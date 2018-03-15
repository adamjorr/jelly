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
import jellyfish
from itertools import repeat
from multiprocessing import Pool

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
            allmers = jellyfish.string_canonicals(read.query_sequence)
            lmers = []
            for mer in allmers:
                lmers.append(str(mer))
            for mer in lmers:
                m = jellyfish.MerDNA(mer)
                m.canonicalize()
                foo = str(m)
                alltable.add(m,1)
                bar = str(m)
                assert foo == bar, "Mer has mutated! before: {}, after: {}.".format(foo,bar)
            refchr = read.reference_name
            refpositions = read.get_reference_positions(full_length = True)
            refpositions = [p for p in refpositions if p not in vcf[refchr]] #ignore sites in VCF
            errorpositions = [i for i, pos in enumerate(refpositions) if pos is None or read.query_sequence[i] != get_ref_base(ref, refchr, pos)]
            if not errorpositions:
                continue
            else:
                mranges = mer_ranges(lmers, ksize)
                errorkmers = [k for i,k in enumerate(lmers) if any([p in mranges[i] for p in errorpositions])]
                for k in errorkmers:
                    m = jellyfish.MerDNA(k)
                    m.canonicalize()
                    errortable.add(m,1)
    return alltable, errortable

def jellyfish_abundances(samfile, conf_regions, totaltable, errortable):
    totalabund, errorabund = [], []
    counted = initialize_hash() #use this hash so we don't double count any kmers
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            allmers = jellyfish.string_canonicals(read.query_sequence)
            lmers = []
            for mer in allmers:
                lmers.append(str(mer))
            for mer in lmers:
                m = jellyfish.MerDNA(mer)
                m.canonicalize()
                foo = str(m)
                if counted.get(m) is None:
                    counted.add(m,1)
                    errorcount = getcount(errortable, m)
                    totalcount = getcount(totaltable, m)
                    bar = str(m)
                    assert foo == bar, "Mer has mutated! before: {}, after: {}.".format(foo,bar)
                    assert errorcount <= totalcount, "Mer {} has errorcount {} and totalcount {}.".format(m, errorcount, totalcount)
                    assert totalcount > 0, "Mer {} has totalcount <= 0. ({})".format(m, totalcount)
                    errorabund.append(errorcount)
                    totalabund.append(totalcount)
                else:
                    continue
    return totalabund, errorabund

def mer_ranges(mers,ksize):
    return [range(k,k+ksize) for k in range(len(mers))]

def initialize_hash():
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

def main():
    # args = argparse() #TODO
    # print(get_kmers_covering("ATCGAA",3,4))
    # print(get_confident_regions('/home/adam/variant-standards/CHM-eval/hg19/chr1/chr1_confident.bed.gz')[0:2])

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
    bedfilename = fileprefix + 'chr1_confident.bed.gz'
    vcffilename = fileprefix + 'chr1_in_confident.vcf.gz'

    #set up hashes
    print(tstamp(), "Preparing hashes . . .", file=sys.stderr)
    khmer.khmer_args.info = newinfo
    args = khmer.khmer_args.build_counting_args().parse_args()
    alltable = khmer.khmer_args.create_countgraph(args)
    errortable = khmer.khmer_args.create_countgraph(args)
    trackingtable = khmer.khmer_args.create_countgraph(args)

    #do things
    print(tstamp(), "Loading Files . . .", file=sys.stderr)
    samfile = pysam.AlignmentFile(samfilename)
    refdict = get_ref_dict(fafilename)
    conf_regions = get_confident_regions(bedfilename)
    vcf = load_vcf(vcffilename, conf_regions)

    print(tstamp(), "Counting . . .", file=sys.stderr)
    alltable, errortable = count_mers(samfile, refdict, vcf, conf_regions, alltable, errortable)

    print(tstamp(), "Calculating Abundances . . .", file=sys.stderr)
    tabund, eabund = get_abundances(samfile, conf_regions, alltable, errortable, trackingtable)
    #each kmer has a position in these arrays; abund[kmer] = # occurrences
    totalabund = np.array(tabund)
    errorabund = np.array(eabund)
    errorweight = np.true_divide(errorabund,totalabund)
    np.set_printoptions(edgeitems=100)

    totalcounts = np.bincount(totalabund)
    errorcounts = np.bincount(errorabund, weights = errorweight)
    errorcounts = np.pad(errorcounts,(0,len(totalcounts)-len(errorcounts)),'constant')
    divisorcounts = np.array(totalcounts)
    divisorcounts[0] = 1
    perror = np.true_divide(errorcounts,divisorcounts) #element-wise division gets probability any kmer in a bin is an error
    #perror[1] = p(error) for abundance of 1

    print(tstamp(), "Making plots . . .", file=sys.stderr)

    sns.set()
    plt.xlim(0,100)
    totalfig = plt.figure()
    # totalax = totalfig.add_subplot(211)
    sns.distplot(totalabund, color = "g", hist_kws = {'alpha' : 0.75}, kde = False)

    # errorax = totalfig.add_subplot(212)
    sns.distplot(errorabund, hist_kws = {'weights' : errorweight, 'alpha' : 0.75}, color = "r", kde = False)
    totalfig.savefig('distributions.png')

    probabilityplot = plt.figure()
    sns.regplot(np.arange(len(perror+1)),perror,fit_reg = False)
    probabilityplot.savefig('probability.png')

if __name__ == '__main__':
    main()
