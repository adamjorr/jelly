#!/usr/bin/env python3
import argparse
import screed
import sys
import khmer
import pysam
import vcf
import pybedtools
from pybedtools import BedTool

from khmer import khmer_args

def subset_vcf(vcf, conf_regions):
	pass

def subset_samfile(samfile, conf_regions):
	pass

def get_errors_at_position(samfile, chr, position):
	pass

def _get_confident_regions(bedfilename):
	'''
	get BED-coordinate regions
	'''
	bed = BedTool(bedfilename)
	d = dict()
	for interval in bed:
		dict.setdefault(interval.chr, list())
		d[interval.chr].append([interval.start, interval.stop])
	return d

def get_confident_regions(bedfilename):
	'''
	get 1-based region strings
	'''
	bed = BedTool(bedfilename)
	return ['{}:{}-{}'.format(interval.chrom, interval.start+1, interval.stop) for interval in bed]

def get_ref_base(reffile, chrom, pos):
	regionstr = '{0}:{1}-{1}'.format(chrom,pos)
	return reffile.fetch(region=regionstr)

#for now ignore sites in the vcf
#return dict of read : pos
def get_erroneous_reads(samfile, ref, conf_regions, vcf, ksize):
	#try this and see if it works
	regionstr = ','.join(get_confident_regions(conf_regions))
	# for regionstr in get_confident_regions(conf_regions):
	for col in samfile.pileup(region = regionstr, truncate = True):
		if col.pos in vcf:
			continue
		else:
			refbase = get_ref_base(ref, col.reference_name, col.pos)
			for read in col.pileups:
				if read.is_del or read.is_refskip or read.alignment.query_sequence[read.query_position] != refbase:
					erroneouskmers = get_kmers_covering(read.alignment.query_sequence, read.query_position, ksize)

def get_kmers_covering(read, pos, ksize):
	return [read[start:start+ksize] for start in range(0,len(read)-ksize+1) if pos >= start and pos < start+ksize]

def main():
	# args = argparse() #TODO
	# print(get_kmers_covering("ATCGAA",3,4))
	# print(get_confident_regions('/home/adam/variant-standards/CHM-eval/hg19/chr1/chr1_confident.bed.gz')[0:2])
	htable = khmer.khmer_args.create_countgraph(3, 100, 1)
	htable.count("ATC")
	print(htable.get("ATC"))
	print(htable.get("ATG"))


if __name__ == '__main__':
    main()
