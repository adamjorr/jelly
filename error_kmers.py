#!/usr/bin/env python3
import argparse
import io
import screed
import sys
import khmer
import pysam
import vcf
import pybedtools
from pybedtools import BedTool
from khmer import khmer_args
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import time
import datetime
#import dna_jellyfish as jellyfish
from itertools import repeat
from multiprocessing import Pool
import scipy
import scipy.stats
import scipy.signal
import scipy.optimize as op
import os.path

def get_confident_regions(bedfilename):
    '''
    get 1-based region strings
    '''
    bed = BedTool(bedfilename)
    return ['{}:{}-{}'.format(interval.chrom, interval.start+1, interval.stop) for interval in bed]


def load_vcf(vcffilename, conf_regions):
    d = dict()
    for regionstr in conf_regions:
        vcf = pysam.VariantFile(vcffilename, drop_samples = True)
        for record in vcf.fetch(region = regionstr):
            d.setdefault(record.chrom, list()).append(record.pos)
    return d

def get_ref_dict(reffilename):
    reffile = pysam.FastaFile(reffilename)
    return {r : reffile.fetch(region=r) for r in reffile.references}

def get_ref_base(refdict, chrom, pos):
    """pos is 1 based"""
    return refdict[chrom][pos - 1]

def split_into_ranges(read, ksize):
    return [range(start,start+ksize) for start in range(len(read)-ksize+1)]

def count_mers(samfile, ref, vcf, conf_regions, alltable, errortable):
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            alltable.consume(read.query_sequence)
            refchr = read.reference_name
            refpositions = read.get_reference_positions(full_length=True) #these should be 1-based positions but are actually 0-based
            errorpositions = [i for i,pos in enumerate(refpositions) if pos is None or (read.query_sequence[i] != ref[refchr][pos] and pos+1 not in vcf[refchr])]
            if not errorpositions:
                continue
            else:
                mers = errortable.get_kmers(read.query_sequence)
                mranges = mer_ranges(mers, errortable.ksize())
                errorkmers = [k for i,k in enumerate(mers) if any([p in mranges[i] for p in errorpositions])]
                for k in errorkmers:
                    errortable.count(k)
    return alltable, errortable

def get_abundances(samfile, conf_regions, totaltable, errortable, trackingtable):
    totalabund, errorabund = [], []
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            kmers = totaltable.get_kmers(read.query_sequence)
            for mer in kmers:
                if not trackingtable.get(mer):
                    trackingtable.count(mer)
                    errorcount = errortable.get(mer)
                    totalcount = totaltable.get(mer)
                    assert errorcount <= totalcount, "Mer {} has errorcount {} and totalcount {}.".format(mer, errorcount, totalcount)
                    assert totalcount > 0, "Mer {} has totalcount <= 0. ({})".format(mer, totalcount)
                    totalabund.append(totalcount)
                    errorabund.append(errorcount)
    return totalabund, errorabund

def jellyfish_count(samfile, ref, vcf, conf_regions, alltable, errortable, ksize):
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            kmers = jellyfish.string_canonicals(read.query_sequence)
            for mer in kmers:
                alltable.add(mer,1)
            refchr = read.reference_name
            refpositions = read.get_reference_positions(full_length = True) #these should be 1-based but are actually 0-based
            errorpositions = [i for i, pos in enumerate(refpositions) if pos is None or (read.query_sequence[i] != ref[refchr][pos] and pos+1 not in vcf[refchr])]
            if not errorpositions:
                continue
            else:
                mers = jellyfish.string_canonicals(read.query_sequence)
                mranges = mer_ranges(mers, ksize)
                errorkmers = [k for i,k in enumerate(mers) if any([p in mranges[i] for p in errorpositions])]
                for k in errorkmers:
                    errortable.add(k,1)
    return alltable, errortable

def jellyfish_abundances(samfile, conf_regions, totaltable, errortable):
    totalabund, errorabund = [], []
    counted = initialize_jf_hash() #use this hash so we don't double count any kmers
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            kmers = jellyfish.string_canonicals(read.query_sequence)
            for mer in kmers:
                if counted.get(mer) is None:
                    counted.add(mer,1)
                    errorcount = getcount(errortable, mer)
                    totalcount = getcount(totaltable, mer)
                    assert errorcount <= totalcount, "Mer {} has errorcount {} and totalcount {}.".format(mer, errorcount, totalcount)
                    assert totalcount > 0, "Mer {} has totalcount <= 0. ({})".format(mer, totalcount)
                    errorabund.append(errorcount)
                    totalabund.append(totalcount)
    return totalabund, errorabund

def mer_ranges(mers,ksize):
    return [range(k,k+ksize) for k in range(len(mers))]

def initialize_jf_hash():
    return jellyfish.HashCounter(1024,15)

def newinfo(*kwargs):
    return

def getcount(table, mer):
    val = table.get(mer) 
    v = (val if val is not None else 0)
    return v

def getcount_manytables(tables, mer):
    values = map(getcount, tables, repeat(mer))
    return sum(values)

def get_readinfo(samfile, region):
    r = []
    for read in samfile.fetch(region=region, multiple_iterators=True):
        r.append((read.query_sequence, read.reference_name, read.get_reference_positions(full_length=True)))
    return r

def tstamp():
    return '[ ' + datetime.datetime.today().isoformat(' ', 'seconds') + ' ]'

def init_jf_hashes():
    print(tstamp(), "Preparing JF hashes . . .", file=sys.stderr)
    jellyfish.MerDNA.k(30)
    alltable = initialize_jf_hash()
    errortable = initialize_jf_hash()
    return alltable, errortable

def init_hashes():
    print(tstamp(), "Preparing hashes . . .", file=sys.stderr)
    khmer.khmer_args.info = newinfo
    args = khmer.khmer_args.build_counting_args().parse_args()
    alltable = khmer.khmer_args.create_countgraph(args)
    errortable = khmer.khmer_args.create_countgraph(args)
    trackingtable = khmer.khmer_args.create_countgraph(args)
    return alltable, errortable, trackingtable

def load_files(samfilename, fafilename, bedfilename, vcffilename):
    print(tstamp(), "Loading Files . . .", file=sys.stderr)
    samfile = pysam.AlignmentFile(samfilename)
    refdict = get_ref_dict(fafilename)
    conf_regions = get_confident_regions(bedfilename)
    vcf = load_vcf(vcffilename, conf_regions)
    return samfile, refdict, conf_regions, vcf    

def plot_dists(totalabund, errorweight, filename):
    print(tstamp(), "Making distribution plots . . .", file=sys.stderr)
    sns.set()
    plt.xlim(0,100)
    totalfig = plt.figure()
    # totalax = totalfig.add_subplot(211)
    h, bins = np.histogram(totalabund, bins = 'auto')
    sns.distplot(totalabund, bins=bins, color = "g", hist_kws = {'alpha' : 0.6}, kde = False, norm_hist = False)

    # errorax = totalfig.add_subplot(212)
    sns.distplot(totalabund, bins=bins, hist_kws={'weights' : errorweight, 'alpha' : 0.6}, color = "r", kde = False, norm_hist=False)
    totalfig.savefig(filename)

def plot_perror(perror, est_perror, filename):
    #TODO: plot a fit of the error model too
    print(tstamp(), "Making abundance error probability plot . . .", file=sys.stderr)

    x = np.arange(len(perror))
    #estimate = scipy.stats.poisson.sf(x,mu=lambda_est)

    sns.set()
    plt.xlim(0,100)
    probabilityplot = plt.figure()
    plt.plot(x,perror, '.-')
    plt.plot(x,est_perror,'m.-')
    plt.legend(labels = ("empirical","estimate"))
    probabilityplot.savefig(filename)

def calc_perror(totalabund, errorabund, distplot = None, errorplot = None):
    tabund = np.array(totalabund)
    eabund = np.array(errorabund)
    errorweight = np.true_divide(eabund,tabund)
    tcounts = np.bincount(totalabund)
    ecounts = np.bincount(totalabund, weights = errorweight)
    ecounts[0] = 0
    tcounts[0] = 1
    #print("Tcounts is not 0? Index of zero:", np.nonzero(tcounts == 0))
    #print(tcounts[0:10])
    perror = np.true_divide(ecounts, tcounts)

    lambda_ests = scipy.signal.argrelmax(tcounts)[0]
    first_lambda = lambda_ests[0]
    x = np.arange(len(tcounts))
    r = 2
    est_perror = .5 * scipy.stats.nbinom.pmf(x, r, r/(first_lambda+r)) + .5 * scipy.stats.uniform.pdf(x, scale = len(x))
    peak = int(first_lambda / r)
    est_perror = est_perror/max(est_perror) * .975
    est_perror[:peak] = .975
    
    if distplot is not None:
        plot_dists(totalabund, errorweight, distplot)
    if errorplot is not None:
        plot_perror(perror, est_perror, errorplot)

    return perror, tcounts

def count_qual_scores(samfile, ref, conf_regions, vcf):
    print(tstamp(), "Counting Base Quality Scores . . .", file=sys.stderr)
    numerrors = np.zeros(41, dtype = np.uint64)
    numtotal = np.zeros(41, dtype = np.uint64)
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            refchr = read.reference_name
            refpositions = read.get_reference_positions(full_length=True) #these should be 1-based positions but are actually 0-based
            errorpositions = np.array([i for i,pos in enumerate(refpositions) if pos is None or (read.query_sequence[i] != ref[refchr][pos] and pos+1 not in vcf[refchr])], dtype=np.intp)
            quals = np.array(read.query_qualities, dtype=np.intp)
            np.add.at(numtotal,quals,1)
            np.add.at(numerrors,quals[errorpositions],1)
    return numerrors, numtotal

def plot_qual_scores(numerrors, numtotal, plotname, plottitle = None):
    print(tstamp(), "Making Base Quality Score Plot . . .", file=sys.stderr)
    plottitle = (plottitle if plottitle is not None else plotname)
    numtotal[numtotal == 0] = 1 #don't divide by 0
    p = numerrors/numtotal
    p[p == 0] = 1e-4 #1e-4 is the largest quality score, 40
    q = -10.0*np.log10(p)
    q = np.array(np.rint(q), dtype=np.int)
    q = np.clip(q, 0, 40)
    x = np.arange(len(p))
    diff = np.mean(np.absolute(x - q))
    
    sns.set()
    qualplot = plt.figure()
    ax = qualplot.add_subplot(111)
    qualplot.suptitle(plottitle)
    plt.plot(x,x)
    plt.plot(x,q)
    plt.xlabel("Predicted Quality Score")
    plt.ylabel("Actual Quality Score")
    plt.legend(labels = ["Perfect","Estimated"], loc = "upper left")
    plt.text(0.5,0.01, "Mean absolute difference: " + str(diff),
        horizontalalignment = 'center', verticalalignment = 'bottom',
        transform = ax.transAxes)
    qualplot.savefig(plotname)

"""
Calculate the log likelihood of a probability vector given a new value of x and i
"""
def calc_loglike(x, i, A, E, pi, ksize):
    #update A where x = pe0
    if not 0 < x < 1:
        return -np.inf
    pe0 = x
    pe1 = A[i+1,0,1]
    A[i+1] = np.array([[1 - pe1, pe1],[pe0 - pe0*pe1, 1-pe0+pe0*pe1]])

    #update A where x = pe1
    pe0 = A[i+1-ksize,1,0] / A[i+1-ksize,0,0]
    pe1 = x
    A[i+1-ksize] = np.array([[1 - pe1, pe1],[pe0 - pe0*pe1, 1-pe0+pe0*pe1]])

    _, normalizer = normalized_forward(A, E, pi)
    return np.sum(np.log(normalizer))

def correct_sam_test(samfile, conf_regions, outfile, ksize, modelA, modelE, modelxi, modelgamma):
    print(tstamp(), "Correcting Input Reads . . .", file=sys.stderr)
    np.seterr(all='raise')
    outsam = pysam.AlignmentFile(outfile, "wb", template=samfile)
    
    i = 0
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            A = modelA[i,:]
            xi = modelxi[i,:]
            gamma = modelgamma[i,:]
            E = modelE[i,:]
            quals = np.array(read.query_qualities, dtype=np.int)
            p = np.array(10.0**(-quals/10.0), dtype=np.longdouble)

            #update probabilities with model
            #p[ksize:] = xi[:,0,1]
            #p[:ksize] = xi[:ksize,1,0] / xi[:ksize,0,0]
            #p[ksize:] = A[0,1]

            #this seems pretty good
            #p[ksize:-ksize] = (xi[ksize:,1,0] + xi[:-ksize,0,1]) / (gamma[ksize:,1,0] + gamma[:-ksize,0,0]) #this seems pretty close? to working
            
            for j in range(len(p - 2*ksize)):
                p_obs_given_note = np.prod(E[j:j+ksize, 0])
                p_obs_given_e = np.prod(E[j:j+ksize, 1])
                p[j + ksize - 1] = p_obs_given_e * p[j + ksize - 1] / (p_obs_given_e + p_obs_given_note)

            try:
                assert np.all(p >= 0.0) and np.all(p <= 1.0)
            except AssertionError:
                print("p:",p)
                print("p[p < 0.0]",p[p < 0.0])
                print("p[p > 1.0]",p[p > 1.0])
                raise

            q = -10.0*np.log10(p)
            quals = np.array(np.rint(q), dtype=np.int)
            quals = np.clip(quals, 0, 40)
            read.query_qualities = quals
            
            outsam.write(read)
            i = i + 1

def baum_welch(A, E, pi):
    """
    Return new A
    """
    loglike = np.NINF
    likelihood_delta = np.inf
    while (likelihood_delta > .1):
        alpha, normalizer = normalized_forward(A, E, pi) #shape = (t,2,1) and (t,1,1)
        beta = normalized_backward(A, E, normalizer) #shape = (t,2,1)
        newloglike = np.sum(np.log(normalizer))

        #the new likelihood should be better than the old one
        try:
           assert newloglike >= loglike
        except AssertionError:
           print("Previous log likelihood:", loglike)
           print("New log likelihood:", newloglike)
           raise
        likelihood_delta = newloglike - loglike
        loglike = newloglike

        xi_num = alpha[:-1,] * A * np.transpose(beta[1:,],(0,2,1)) * np.transpose(E[1:,],(0,2,1))
        xi_denom = np.sum(xi_num, axis = (1,2), keepdims = True)

        #make sure we don't underflow
        try:
            xi = np.array(xi_num / xi_denom, dtype = np.longdouble) #p(state t = i and state t+1 = j | Obs, parameters); it makes sense from t=0 to length of the state sequence - 1
        except FloatingPointError:
            print("Xi_num:", xi_num)
            print("Xi_denom:",xi_denom)
            raise

        gamma = np.sum(xi, axis = 2, keepdims = True)
        pi = gamma[0,]
        update = np.sum(xi, axis = 0) / np.sum(gamma, axis = 0)
        try:
            assert np.all(update >= 0.0) and np.all(update <= 1.0)
        except AssertionError:
            print("xi:", xi)
            print("gamma:", gamma)
            raise
        A = np.array(update, copy = True)
    return A, xi, gamma

def train_model(samfile, conf_regions, tcounts, perror, kgraph):
    np.seterr(all='raise')
    ksize = kgraph.ksize()
    ecounts = tcounts * perror
    p_e_given_a = np.array(ecounts/tcounts, dtype = np.longdouble)
    p_a_given_e = np.array(ecounts/np.nansum(ecounts), dtype = np.longdouble)
    p_a_given_note = np.array((tcounts-ecounts) / np.nansum(tcounts-ecounts), dtype = np.longdouble)
    #fix zeros
    p_a_given_e[p_a_given_e == 1.0] = 0.9999
    p_a_given_note[p_a_given_note == 0.0] = 0.0001

    models_A = list()
    models_xi = list()
    models_gamma = list()
    models_E = list()
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            kmers = kgraph.get_kmers(read.query_sequence)
            counts = np.array(list(map(kgraph.get, kmers)), dtype = np.int)
            quals = np.array(read.query_qualities, dtype=np.int)
            
            p = np.array(10.0**(-quals/10.0), dtype=np.longdouble)
            A = np.zeros((2,2), dtype=np.longdouble)
            E = np.zeros((len(counts),2,1), dtype=np.longdouble)
            pi = np.array([[p_e_given_a[counts[0]]],[1.0-p_e_given_a[counts[0]]]], dtype=np.longdouble) #probability for 1st state
            for j, count in enumerate(counts): #the emission matrix is of length counts, the transition matrix is of length counts - 1
                pe0 = p[j-1] #A[0] will not make any sense because of this
                pe1 = p[j+ksize-1]
                try:
                    assert 0.0 <= pe0 <= 1.0
                    assert 0.0 <= pe1 <= 1.0
                except AssertionError:
                    print(z)
                    print("pe0:",pe0)
                    print("pe1:",pe1)
                    print("p:",p)
                    raise
                A += np.array([[1.0 - pe1, pe1],[pe0 - pe0*pe1, 1.0 - pe0+pe0*pe1]], dtype=np.longdouble, copy = True) #A is size len(counts), but A[0] is meaningless
                E[j] = np.array([[p_a_given_note[count]],[p_a_given_e[count]]], dtype=np.longdouble, copy = True) #E is size len(counts)
            A = A / len(counts)
            A, xi, gamma = baum_welch(A, E, pi)
            models_A.append(A)
            models_xi.append(xi)
            models_gamma.append(gamma)
            models_E.append(E)
    Astack = np.stack(models_A)
    xistack = np.stack(models_xi)
    gammastack = np.stack(models_gamma)
    Estack = np.stack(models_E)
    return Astack, Estack, xistack, gammastack


def normalized_forward(A, E, pi):
    """
    This is a normalized forward algorithm, normalized by P(O(t)|theta).
    This is sometimes called a filtering recursion.
    """
    normalizer = np.ones((E.shape[0],1,1), dtype = np.longdouble)
    alpha = np.zeros((E.shape[0],2,1), dtype = np.longdouble)
    normalizer[0] = np.sum(pi * E[0])
    alpha[0,] = pi * E[0] / normalizer[0]
    for t in range(1, E.shape[0]):
        unscaled_alpha = np.matmul(np.transpose(A),alpha[t-1]) * E[t]
        normalizer[t] = np.sum(unscaled_alpha)
        alpha[t] = unscaled_alpha / normalizer[t]
    return alpha, normalizer

def normalized_backward(A, E, normalizer):
    """
    This is a normalized backward algorithm normalized by P(O(t)|theta) computed during the normalized forward algorithm.
    It is convenient to use this normalization factor because it cancels out during gamma and xi calculation.
    """
    beta = np.zeros((E.shape[0],2,1), dtype = np.longdouble)
    beta[-1,] = np.array([[1.0],[1.0]])/normalizer[-1]
    for t in reversed(range(0,E.shape[0]-1)):
        beta[t] = np.matmul(A,E[t+1] * beta[t+1]) / normalizer[t]
    return beta

def t_baum_welch(A, E, pi):
    """
    Return new A
    """
    loglike = np.NINF
    likelihood_delta = np.inf
    while (likelihood_delta > .1):
        alpha, normalizer = normalized_forward(A, E, pi) #shape = (t,2,1) and (t,1,1)
        beta = normalized_backward(A, E, normalizer) #shape = (t,2,1)
        newloglike = np.sum(np.log(normalizer))

        #the new likelihood should be better than the old one
        try:
           assert newloglike >= loglike
        except AssertionError:
           print("Previous log likelihood:", loglike)
           print("New log likelihood:", newloglike)
           raise
        likelihood_delta = newloglike - loglike
        loglike = newloglike

        xi_num = alpha[:-1,] * A[1:,] * np.transpose(beta[1:,],(0,2,1)) * np.transpose(E[1:,],(0,2,1))
        xi_denom = np.sum(xi_num, axis = (1,2), keepdims = True)

        #make sure we don't underflow
        try:
            xi = np.array(xi_num / xi_denom, dtype = np.longdouble) #p(state t = i and state t+1 = j | Obs, parameters); it makes sense from t=0 to length of the state sequence - 1
        except FloatingPointError:
            print("Xi_num:", xi_num)
            print("Xi_denom:",xi_denom)
            raise

        gamma = np.sum(xi, axis = 2, keepdims = True)
        pi = gamma[0,]
        update = xi / gamma
        try:
            assert np.all(update >= 0.0) and np.all(update <= 1.0)
        except AssertionError:
            print("xi:", xi)
            print("gamma:", gamma)
            raise
        A[1:,] = np.array(update, copy = True)
    return A

#numpy arrays are n x m, row by column; x[1,2] is 2nd row 3rd col
def t_normalized_forward(A, E, pi):
    """
    This is a normalized forward algorithm, normalized by P(O(t)|theta).
    This is sometimes called a filtering recursion.
    """
    normalizer = np.ones((E.shape[0],1,1), dtype = np.longdouble)
    alpha = np.zeros((E.shape[0],2,1), dtype = np.longdouble)
    normalizer[0] = np.sum(pi * E[0])
    try:
        alpha[0,] = pi * E[0] / normalizer[0]
    except FloatingPointError:
        print("alpha[0,]:",alpha[0,])
        print("E[0]:",E[0])
        print("normalizer[0]",normalizer[0])
    for t in range(1, E.shape[0]):
        unscaled_alpha = np.matmul(np.transpose(A[t]),alpha[t-1]) * E[t]
        normalizer[t] = np.sum(unscaled_alpha)
        alpha[t] = unscaled_alpha / normalizer[t]
    return alpha, normalizer

def t_normalized_backward(A, E, normalizer):
    """
    This is a normalized backward algorithm normalized by P(O(t)|theta) computed during the normalized forward algorithm.
    It is convenient to use this normalization factor because it cancels out during gamma and xi calculation.
    """
    beta = np.zeros((E.shape[0],2,1), dtype = np.longdouble)
    beta[-1,] = np.array([[1.0],[1.0]])/normalizer[-1]
    for t in reversed(range(0,E.shape[0]-1)):
        beta[t] = np.matmul(A[t+1],E[t+1] * beta[t+1]) / normalizer[t]
    return beta

def forward(A, E, pi):
    alpha = np.zeros((E.shape[0],2,1), dtype = np.longdouble)
    alpha[0,] = pi * E[0] #pi is p(err, 1st kmer), which is not in E
    for t in range(1, E.shape[0]): #the first element of the shape is the number of 2 x 1 matrices we have
        alpha[t] = np.matmul(np.transpose(A[t]),alpha[t-1]) * E[t]
    return alpha #alpha is the number of transitions + 1, = to the number of states

def backward(A, E):
    beta = np.zeros((E.shape[0],2,1), dtype = np.longdouble)
    beta[-1,] = np.array([[1],[1]])
    for t in reversed(range(0,E.shape[0]-1)):
        beta[t] = np.matmul(A[t+1],E[t+1] * beta[t+1])
    return beta

def main():
    # args = argparse() #TODO
    # print(get_kmers_covering("ATCGAA",3,4))
    # print(get_confident_regions('/home/adam/variant-standards/CHM-eval/hg19/chr1/chr1_confident.bed.gz')[0:2])

    #foo = "ATCGT"
    #jellyfish.MerDNA.k(4)
    #bar = jellyfish.HashCounter(1024,15)
    #baz = jellyfish.string_mers(foo)
    #for mer in baz:
    #    print(str(mer))
    #exit()

    #example
    # khmer.khmer_args.info = newinfo
    # args = khmer.khmer_args.build_counting_args().parse_args()
    # htable = khmer.khmer_args.create_countgraph(args, ksize=3)
    # htable.count("ATC")
    # htable.count("ATG")
    # htable.count("ATG")
    # print(htable.get("ATC"))
    # print(htable.get("ATG"))

    #arguments
    fileprefix = '/home/ajorr1/variant-standards/CHM-eval/hg19/chr1/'
    samfilename = fileprefix + 'chr1.bam'
    fafilename = fileprefix + 'chr1.renamed.fa'
    bedfilename = fileprefix + 'chr1_first100.bed.gz'
    vcffilename = fileprefix + 'chr1_in_confident.vcf.gz'
    tabundfile = fileprefix + 'tabund.txt.gz'
    eabundfile = fileprefix + 'eabund.txt.gz'
    numerrsfile = fileprefix + 'numerrs.txt.gz'
    numtotalfile = fileprefix + 'numtotal.txt.gz'
    kmergraphfile = fileprefix + 'kmers.khmer'
    outfile = fileprefix + 'test.bam'
    modelfile = fileprefix + 'model.npz'
    np.set_printoptions(edgeitems=100)
    
    
    if all([ os.path.exists(tabundfile), os.path.exists(eabundfile),
    os.path.exists(numerrsfile), os.path.exists(numtotalfile),
    os.path.exists(kmergraphfile) ]) :
        tabund = np.loadtxt(tabundfile, dtype = np.int64)
        eabund = np.loadtxt(eabundfile, dtype = np.int64)
        numerrors = np.loadtxt(numerrsfile, dtype = np.int64)
        numtotal = np.loadtxt(numtotalfile, dtype = np.int64)
        
        alltable, _, _ = init_hashes()
        alltable.load(kmergraphfile)
    else:
        #set up hashes and load files
        alltable, errortable, trackingtable = init_hashes()
        samfile, refdict, conf_regions, vcf = load_files(samfilename, fafilename, bedfilename, vcffilename)
        
        #count
        print(tstamp(), "Counting . . .", file=sys.stderr)
        alltable, errortable = count_mers(samfile, refdict, vcf, conf_regions, alltable, errortable)
        alltable.save(kmergraphfile)

        print(tstamp(), "Calculating Abundances . . .", file=sys.stderr)
        #each kmer has a position in these arrays; abund[kmer idx] = # occurrences
        tabund, eabund = get_abundances(samfile, conf_regions, alltable, errortable, trackingtable)
        np.savetxt(tabundfile, tabund, fmt = '%d')
        np.savetxt(eabundfile, eabund, fmt = '%d')
        
        numerrors, numtotal = count_qual_scores(samfile, refdict, conf_regions, vcf)
        np.savetxt(numerrsfile, numerrors, fmt = '%d')
        np.savetxt(numtotalfile, numtotal, fmt = '%d')
    

    #tabund[i] = count of kmer with index i
    #tcounts[i] = number of kmers with count i
    
    #perror[1] = observed p(error|x) for abundance of 1
    perror, tcounts = calc_perror(tabund, eabund, distplot = 'distributions.png', errorplot = 'probability.png')
    
    samfile, refdict, conf_regions, vcf = load_files(samfilename, fafilename, bedfilename, vcffilename)
    
    if os.path.exists(modelfile):
        print(tstamp(), "Loading model . . .", file = sys.stderr)
        loadedmodel = np.load(modelfile)
        A = np.array(loadedmodel['A'], dtype = np.longdouble)
        xi = np.array(loadedmodel['xi'], dtype = np.longdouble)
        gamma = np.array(loadedmodel['gamma'], dtype = np.longdouble)
        E = np.array(loadedmodel['E'], dtype = np.longdouble)
    else:
        print(tstamp(), "Training model . . .", file = sys.stderr)
        A, E, xi, gamma = train_model(samfile, conf_regions, tcounts, perror, alltable)
        np.savez_compressed(modelfile, A = A, E = E, xi = xi, gamma = gamma)

    correct_sam_test(samfile, conf_regions, outfile, alltable.ksize(), A, E, xi, gamma) #creates outfile
    pysam.index(outfile)
    correctederrs, correctedtot = count_qual_scores(pysam.AlignmentFile(outfile),refdict, conf_regions, vcf)
    
    plot_qual_scores(numerrors, numtotal, "qualscores.png", "Raw Reads")
    plot_qual_scores(correctederrs, correctedtot, "corrected.png", "After Correction")
    print(tstamp(), "Done", file = sys.stderr)

if __name__ == '__main__':
    main()
