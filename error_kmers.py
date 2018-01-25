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
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import datetime

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

def get_ref_base(reffile, chrom, pos):
    regionstr = '{0}:{1}-{1}'.format(str(chrom),pos)
    return reffile.fetch(region=regionstr)

def get_kmers_covering(read, pos, ksize):
    return [read[start:start+ksize] for start in range(0,len(read)-ksize+1)]

def split_into_kmers(read, ksize):
    return {read[start:start+ksize] : range(start,start+ksize) for start in range(0,len(read)-ksize+1)}

def count_mers(samfile, ref, vcf, conf_regions, alltable, errortable):
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            alltable.consume(read.query_sequence)
            refchr = read.reference_name
            refpositions = read.get_reference_positions(full_length=True)
            errorpositions = [i for i, pos in enumerate(refpositions) if pos is None or read.query_sequence[i] != get_ref_base(reffile, refchr, pos)]
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

def newinfo(*kwargs):
    return

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
    fileprefix = '/home/adam/variant-standards/CHM-eval/hg19/chr1/'
    samfilename = fileprefix + 'chr1.bam'
    fafilename = fileprefix + 'chr1.fa'
    bedfilename = fileprefix + 'chr1_confident.bed.gz'
    vcffilename = fileprefix + 'chr1_in_confident.vcf.gz'

    #set up hashes
    print(sys.stderr, '[',datetime.datetime.today().isoformat(' ', 'seconds'), ']', "Preparing hashes . . .")
    khmer.khmer_args.info = newinfo
    args = khmer.khmer_args.build_counting_args().parse_args()
    alltable = khmer.khmer_args.create_countgraph(args)
    errortable = khmer.khmer_args.create_countgraph(args)

    #do things
    print(sys.stderr, '[',datetime.datetime.today().isoformat(' ', 'seconds'), ']', "Loading Files . . .")
    samfile = pysam.AlignmentFile(samfilename)
    reffile = pysam.FastaFile(fafilename)
    conf_regions = get_confident_regions(bedfilename)
    vcf = load_vcf(vcffilename, conf_regions)

    print(sys.stderr, '[',datetime.datetime.today().isoformat(' ', 'seconds'), ']', "Counting . . .")
    alltable, errortable = count_mers(samfile, reffile, vcf, conf_regions, alltable, errortable)

    print(sys.stderr, '[',datetime.datetime.today().isoformat(' ', 'seconds'), ']', "Calculating Abundances . . .")
    totalabund, errorabund = get_abundances(samfile, conf_regions, alltable, errortable)

    print(totalabund[0:10])

    totalplot = sns.distplot(totalabund, kde=True, color = "m")
    totalplot.get_figure().savefig('totalabund.png')

    errorplot = sns.distplot(errorabund, kde=True, color = "r")
    errorplot.get_figure().savefig('errorabund.png')    


if __name__ == '__main__':
    main()
