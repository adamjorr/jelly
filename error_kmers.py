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
            lmers = []
            for mer in allmers:
                alltable.add(mer,1)
                lmers.append(mer)
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

# def jellyfish_countregion(readinfo, ref, vcf, ksize):
#     print("This worker has " + len(readinfo) + " reads!")
#     kmers, errorkmers = [], []
#     for read in readinfo:
#         query_sequence = read[0]
#         reference_name = read[1]
#         reference_positions = read[2]
#
#         alrlmers = jellyfish.string_canonicals(query_sequence)
#         lmers = []
#         for mer in allmers:
#             m = str(mer)
#             lmers.append(m)
#             kmers.append(m)
#         refpositions = [p for p in reference_positions if p not in vcf[reference_name]]
#         errorpositions = [i for i,pos in enumerate(refpositions) if pos is None or query_sequence[i] != get_ref_base(ref,reference_name,pos)]
#         if not errorpositions:
#             continue
#         else:
#             mranges = mer_ranges(lmers, ksize)
#             errormers = [k for i,k in enumerate(lmers) if any([p in mranges[i] for p in errorpositions])]
#             errorkmers.extend(errormers)
#     return (kmers, errorkmers)

def jellyfish_abundances(samfile, conf_regions, totaltable, errortable):
    totalabund, errorabund = [], []
    counted = initialize_hash() #use this hash so we don't double count any kmers
    for regionstr in conf_regions:
        for read in samfile.fetch(region=regionstr):
            allmers = jellyfish.string_canonicals(read.query_sequence)
            for mer in allmers:
                if counted.get(mer) is None:
                    counted.add(m,1)
                    errorcount = getcount(errortable, m)
                    totalcount = getcount(totaltable, m)
                    assert errorcount <= totalcount, "Mer {} has errorcount {} and totalcount {}.".format(m, errorcount, totalcount)
                    assert totalcount > 0, "Mer {} has totalcount <= 0; processed {} kmers.".format(m,counter)
                    errorabund.append(errorcount)
                    totalabund.append(totalcount)
                else:
                    continue
    return totalabund, errorabund

# def jellyfish_abundances(samfile, conf_regions, totaltable, errortable):
#     totalabund = [getcount(totaltable, mer) for regionstr in conf_regions for read in samfile.fetch(region = regionstr) for mer in jellyfish.string_canonicals(read.query_sequence)]
#     errorabund = [getcount(errortable, mer) for regionstr in conf_regions for read in samfile.fetch(region = regionstr) for mer in jellyfish.string_canonicals(read.query_sequence)]
#     return totalabund, errorabund

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
    bedfilename = fileprefix + 'chr1_first3.bed.gz'
    vcffilename = fileprefix + 'chr1_in_confident.vcf.gz'

    #set up hashes
    print(tstamp(), "Preparing hashes . . .", file=sys.stderr)
    # khmer.khmer_args.info = newinfo
    # args = khmer.khmer_args.build_counting_args().parse_args()
    # alltable = khmer.khmer_args.create_countgraph(args)
    # errortable = khmer.khmer_args.create_countgraph(args)
    jellyfish.MerDNA.k(30)
    alltable = initialize_hash()
    errortable = initialize_hash()


    #do things
    print(tstamp(), "Loading Files . . .", file=sys.stderr)
    samfile = pysam.AlignmentFile(samfilename)
    refdict = get_ref_dict(fafilename)
    conf_regions = get_confident_regions(bedfilename)
    vcf = load_vcf(vcffilename, conf_regions)


    print(tstamp(), "Counting . . .", file=sys.stderr)
    alltable, errortable = jellyfish_count(samfile, refdict, vcf, conf_regions, alltable, errortable,jellyfish.MerDNA.k())

    #do not do threads, something is wrong with parsing the samfile such that the number of kmers present in a region varies
    # with Pool(processes=16) as pool:
    #     results = [pool.apply_async(jellyfish_countregion, [get_readinfo(samfile, r),refdict,vcf,jellyfish.MerDNA.k()]) for r in conf_regions]
    #     progress = 1
    #     for r in results:
    #         print(tstamp(), "Job {} of {}".format(progress, len(results)), file=sys.stderr)
    #         progress += 1
    #         totalmers, errormers = r.get()
    #         for mer in totalmers:
    #             m = jellyfish.MerDNA(mer)
    #             #m.canonicalize() the mer should already be in canonical form
    #             alltable.add(m,1)
    #         for mer in errormers:
    #             m = jellyfish.MerDNA(mer)
    #             # m.canonicalize()
    #             errortable.add(m,1)

    print(tstamp(), "Calculating Abundances . . .", file=sys.stderr)
    totalabund, errorabund = jellyfish_abundances(samfile, conf_regions, alltable, errortable)
    #each kmer has a position in these arrays; abund[kmer] = # occurrences
    totalabund = np.array(totalabund)
    errorabund = np.array(errorabund)
    errorweight = errorabund / totalabund


    totalcounts = np.bincount(totalabund)
    errorcounts = np.bincount(errorabund, weights = errorweight)
    perror = errorcounts / totalcounts #element-wise division gets probability any kmer in a bin is an error
    #perror[1] = p(error) for abundance of 1
    # print(errorcounts)
    # print(totalcounts)
    # print(perror)

    print(tstamp(), "Making plots . . .", file=sys.stderr)

    plt.xlim(0,100)
    totalfig = plt.figure()
    totalax = totalfig.add_subplot(211)
    sns.distplot(totalabund, ax=totalax, bins = 100, kde=False, hist_kws={'range' : (0,100)}, color = "g")

    errorax = totalfig.add_subplot(212)
    sns.distplot(errorabund, bins = 100, kde=False, hist_kws={'range' : (0,100)}, color = "r", ax = errorax)
    totalfig.savefig('distributions.png')

    probabilityplot = plt.figure()
    pax = probabilityplot.add_axes()
    sns.barplot(np.arange(len(perror+1)),perror, ax=pax)
    probabilityplot.savefig('probability.png')

if __name__ == '__main__':
    main()
