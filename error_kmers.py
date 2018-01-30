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

def get_kmers_covering(read, pos, ksize):
    return [read[start:start+ksize] for start in range(0,len(read)-ksize+1)]

def split_into_kmers(read, ksize):
    return {read[start:start+ksize] : range(start,start+ksize) for start in range(0,len(read)-ksize+1)}

def count_mers(samfile, ref, vcf, conf_regions, alltable, errortable):
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            alltable.consume(read.query_sequence)
            refchr = read.reference_name
            refpositions = read.get_reference_positions(full_length=True) #these are 1-based positions
            errorpositions = [i for i, pos in enumerate(refpositions) if pos is None or read.query_sequence[i] != get_ref_base(ref, refchr, pos)]
            if not errorpositions:
                continue
            else:
                kmerdict = split_into_kmers(read.query_sequence, errortable.ksize())
                errorkmers = [k for k,v in kmerdict.items() if any([p in v for p in errorpositions])]
                for k in errorkmers:
                    errortable.add(k)
    return alltable, errortable

def get_abundances(samfile, conf_regions, totaltable, errortable):
    totalabund, errorabund = [], []
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            totalabund.extend(totaltable.get_kmer_counts(read.query_sequence))
            errorabund.extend(errortable.get_kmer_counts(read.query_sequence))
    return totalabund, errorabund

def jellyfish_count(samfile, ref, vcf, conf_regions, alltable, errortable, ksize):
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            allmers = jellyfish.string_canonicals(read.query_sequence)
            lmers = [mer for mer in allmers] #turn iterator into list
            for mer in lmers:
                alltable.add(mer,1)
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
                    errortable.add(mer,1)
    return alltable, errortable

def jellyfish_countregion(reads, ref, vcf, ksize):
    totaltable = initialize_hash()
    errortable = initialize_hash()
    for read in reads:
        allmers = jellyfish.string_canonicals(read.query_sequence)
        lmers = [mer for mer in allmers]
        for mer in lmers:
            totaltable.add(mer,1)
        refchr = read.reference_name
        refpositions = read.get_reference_positions(full_length = True)
        refpositions = [p for p in refpositions if p not in vcf[refchr]]
        errorpositions = [i for i,pos in enumerate(refpositions) if pos is None or read.query_sequence[i] != get_ref_base(ref,refchr,pos)]
        if not errorpositions:
            continue
        else:
            mranges = mer_ranges(lmers, ksize)
            errorkmers = [k for i,k in enumerate(lmers) if any([p in mranges[i] for p in errorpositions])]
            for k in errorkmers:
                errortable.add(mer,1)
    return (alltable, errortable)

def jellyfish_abundances(samfile, conf_regions, totaltable, errortable):
    totalabund, errorabund = [], []
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            allmers = jellyfish.string_canonicals(read.query_sequence)
            for mer in allmers:
                totalabund.append(getcount(totaltable, mer))
                errorabund.append(getcount(errortable, mer))
    return totalabund, errorabund

def jellyfish_threaded_abundances(samfile, conf_regions, totaltables, errortables):
    totalabund, errorabund = [], []
    for regionstr in conf_regions:
        for read in samfile.fetch(region = regionstr):
            allmers = jellyfish.string_canonicals(read.query_sequence)
            for mer in allmers:
                totalabund.append(getcount_manytables(totaltables, mer))
                errorabund.append(getcount_manytables(errortables, mer))
    return totalabund, errorabund

def mer_ranges(mers,ksize):
    return [range(k,k+ksize) for k in range(len(mers))]

def initialize_hash():
    return jellyfish.HashCounter(1024,15)

def newinfo(*kwargs):
    return

def getcount(table, mer):
    val = table.get(mer)
    return val if val is not None else 0

def getcount_manytables(tables, mer):
    values = map(getcount, tables, repeat(mer))
    return sum(values)

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
    bedfilename = fileprefix + 'chr1_first100.bed.gz'
    vcffilename = fileprefix + 'chr1_in_confident.vcf.gz'

    #set up hashes
    print('[',datetime.datetime.today().isoformat(' ', 'seconds'), ']', "Preparing hashes . . .", file=sys.stderr)
    # khmer.khmer_args.info = newinfo
    # args = khmer.khmer_args.build_counting_args().parse_args()
    # alltable = khmer.khmer_args.create_countgraph(args)
    # errortable = khmer.khmer_args.create_countgraph(args)
    jellyfish.MerDNA.k(30)
    alltable = initialize_hash()
    errortable = initialize_hash()


    #do things
    print('[',datetime.datetime.today().isoformat(' ', 'seconds'), ']', "Loading Files . . .", file=sys.stderr)
    samfile = pysam.AlignmentFile(samfilename)
    refdict = get_ref_dict(fafilename)
    conf_regions = get_confident_regions(bedfilename)
    vcf = load_vcf(vcffilename, conf_regions)


    print('[',datetime.datetime.today().isoformat(' ', 'seconds'), ']', "Counting . . .", file=sys.stderr)
    # alltable, errortable = jellyfish_count(samfile, refdict, vcf, conf_regions, alltable, errortable,jellyfish.MerDNA.k())

    #try threaded
    with Pool(processes=16) as pool:
        results = [pool.apply_async(jellyfish_countregion, reads = samfile.fetch(region=r), ref = refdict, vcf = vcf, ksize = jellyfish.MerDNA.k()) for r in conf_regions]
        pool.close()
        pool.join()
        alltables = [r.get() for r in results]
    totaltables, errortables = zip(*alltables) 

    print('[',datetime.datetime.today().isoformat(' ', 'seconds'), ']', "Calculating Abundances . . .", file=sys.stderr)
    totalabund, errorabund = jellyfish_threaded_abundances(samfile, conf_regions, totaltables, errortables)

    print(totalabund[0:10])
    print(errorabund[0:10])

    plt.xlim(0,1000)
    totalplot = sns.distplot(totalabund, kde=True, color = "g")
    totalplot.get_figure().savefig('totalabund.png')

    errorplot = sns.distplot(errorabund, kde=True, color = "r")
    errorplot.get_figure().savefig('errorabund.png')


if __name__ == '__main__':
    main()
