#!/usr/bin/env python3
import argparse
import screed
import sys
import khmer
import pysam






def get_erroneous_reads(graphfile,samfile):
	htable = khmer.load_countgraph(graphfile)

def get_erroneous_kmers(graphfile, samfile):
	ereads = get_erroneous_reads(graphfile,samfile)
	return [k for kmer in foo(ereads) for k in kmer]



def main():
	pass


if __name__ == '__main__':
    main()
