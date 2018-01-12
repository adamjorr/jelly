#!/usr/bin/env python3
import argparse
import screed
import sys
import khmer
import pysam
import vcf
from pybedtools import BedTool
import pybedtools

def subset_vcf(vcf, conf_regions):

def subset_samfile(samfile, conf_regions):

def get_errors_at_position(samfile, chr, position):

def get_confident_regions(bedfilename):
	bed = BedTool(bedfilename)
	d = dict()
	for interval in bed:
		dict.setdefault(interval.chr, list())
		d[interval.chr].append(range(interval.start, interval.stop))
	return d

#for now ignore sites in the vcf
#return dict of read : pos
def get_erroneous_reads(samfile, conf_regions, vcf):
	

def get_kmers_covering(read, pos, ksize):
	return [read[start:start+ksize] for start in range(0,len(read)-ksize+1) if pos >= start and pos < start+ksize]

def get_erroneous_kmers(samfile):
	ereads = get_erroneous_reads(samfile)
	return [kmer for read, pos in ereads for kmer in get_kmers_covering(read, pos, ksize)]

def main():
	args = argparse() #TODO
	# print(get_kmers_covering("ATCGAA",3,4))


if __name__ == '__main__':
    main()
