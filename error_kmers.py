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
import datetime
#import dna_jellyfish as jellyfish
from itertools import repeat
from multiprocessing import Pool
import scipy
import scipy.stats
import scipy.signal
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

    return perror

def count_qual_scores(samfile, ref, conf_regions, vcf):
    print(tstamp(), "Counting Base Quality Scores . . .", file=sys.stderr)
    numerrors = np.zeros(40, dtype = np.uint64)
    numtotal = np.zeros(40, dtype = np.uint64)
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            refchr = read.reference_name
            refpositions = read.get_reference_positions(full_length=True) #these should be 1-based positions but are actually 0-based
            errorpositions = np.array([i for i,pos in enumerate(refpositions) if pos is None or (read.query_sequence[i] != ref[refchr][pos] and pos+1 not in vcf[refchr])], dtype=np.intp)
            quals = np.array(read.query_qualities, dtype=np.intp)
            np.add.at(numtotal,quals,1)
            np.add.at(numerrors,quals[errorpositions],1)
    return numerrors, numtotal

def plot_qual_scores(numerrors, numtotal, plotname):
    print(tstamp(), "Making Base Quality Score Plot . . .", file=sys.stderr)
    p = numerrors/numtotal
    y = np.arange(len(p))
    x = 10**(-y/10)
    
    sns.set()
    qualplot = plt.figure()
    plt.plot(x,x)
    plt.plot(x,p)
    qualplot.savefig(plotname)

def correct_sam_test(samfile, conf_regions, outfile, tabund, perror):
    largestidx = 0
    kmerdex = dict()
    khmer.khmer_args.info = newinfo
    args = khmer.khmer_args.build_counting_args().parse_args()
    trackingtable = khmer.khmer_args.create_countgraph(args)
    ksize = trackingtable.ksize()
    outsam = pysam.AlignmentFile(outfile, "wb", template=samfile)
    
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            kmers = trackingtable.get_kmers(read.query_sequence)
            quals = np.array(read.query_qualities, dtype=np.int)
            for j, mer in enumerate(kmers):
                if not trackingtable.get(mer):
                    trackingtable.count(mer)
                    i = largestidx
                    kmerdex[mer] = i
                    largestidx += 1
                else:
                    i = kmerdex[mer]
                pe_given_abund = np.float64(perror[tabund[i]])
                floatinfo = np.finfo(np.float64)
                pe_given_abund = np.clip(pe_given_abund,floatinfo.tiny,1)
                p = 10.0**(-quals[j:j+ksize]/10.0) #convert to probability
                newp = np.true_divide(p,pe_given_abund) #in case pe_given_abund == 0
                newp = np.clip(newp, floatinfo.tiny, 1)
                q = -10.0*np.log10(newp)
                quals[j:j+ksize] = np.rint(q)
                
            print(quals)
            read.query_qualities = quals
            outsam.write(read)
    
    
    
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
    outfile = fileprefix + 'test.bam'
    np.set_printoptions(edgeitems=100)
    
    
    if os.path.exists(tabundfile) and os.path.exists(eabundfile) and os.path.exists(numerrsfile) and os.path.exists(numtotalfile):
        tabund = np.loadtxt(tabundfile, dtype = np.int64)
        eabund = np.loadtxt(eabundfile, dtype = np.int64)
        numerrors = np.loadtxt(numerrsfile, dtype = np.int64)
        numtotal = np.loadtxt(numtotalfile, dtype = np.int64)
    else:
        #set up hashes and load files
        alltable, errortable, trackingtable = init_hashes()
        samfile, refdict, conf_regions, vcf = load_files(samfilename, fafilename, bedfilename, vcffilename)
        
        #count
        print(tstamp(), "Counting . . .", file=sys.stderr)
        alltable, errortable = count_mers(samfile, refdict, vcf, conf_regions, alltable, errortable)

        print(tstamp(), "Calculating Abundances . . .", file=sys.stderr)
        #each kmer has a position in these arrays; abund[kmer idx] = # occurrences
        tabund, eabund = get_abundances(samfile, conf_regions, alltable, errortable, trackingtable)
        np.savetxt(tabundfile, tabund, fmt = '%d')
        np.savetxt(eabundfile, eabund, fmt = '%d')
        
        numerrors, numtotal = count_qual_scores(samfile, refdict, conf_regions, vcf)
        np.savetxt(numerrsfile, numerrors, fmt = '%d')
        np.savetxt(numtotalfile, numtotal, fmt = '%d')
        
    #perror[1] = observed p(error|x) for abundance of 1
    perror = calc_perror(tabund, eabund, distplot = 'distributions.png', errorplot = 'probability.png')
    
    samfile, refdict, conf_regions, vcf = load_files(samfilename, fafilename, bedfilename, vcffilename)
    correct_sam_test(samfile, conf_regions, outfile, tabund, perror) #creates outfile
    correctederrs, correctedtot = count_qual_scores(pysam.AlignmentFile(outfile),refdict, conf_regions, vcf)
    
    plot_qual_scores(numerrors, numtotal, "qualscores.png")
    plot_qual_scores(correctederrs, correctedtot, "corrected.png")
    print(tstamp(), "Done", file = sys.stderr)

if __name__ == '__main__':
    main()
