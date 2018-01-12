#!/usr/bin/env python3
import argparse
import screed
import sys
import khmer
import pysam


def subset_vcf(vcf, conf_regions):

def subset_samfile(samfile, conf_regions):

def get_erroneous_reads(graphfile, samfile, vcf):
	htable = khmer.load_countgraph(graphfile)

def get_kmers_covering(read, pos, ksize):
	return [read[start:start+ksize] for start in range(0,len(read)-ksize+1) if pos >= start and pos < start+ksize]

def get_erroneous_kmers(graphfile, samfile):
	ereads = get_erroneous_reads(graphfile,samfile)
	return [kmer for read, pos in ereads for kmer in get_kmers_covering(read, pos, ksize)]

def main():
	print(get_kmers_covering("ATCGAA",3,4))


if __name__ == '__main__':
    main()
