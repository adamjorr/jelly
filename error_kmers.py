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


# def _get_confident_regions(bedfilename):
# 	'''
# 	get BED-coordinate regions
# 	'''
# 	bed = BedTool(bedfilename)
# 	d = dict()
# 	for interval in bed:
# 		dict.setdefault(interval.chr, list())
# 		d[interval.chr].append([interval.start, interval.stop])
# 	return d

def get_confident_regions(bedfilename):
	'''
	get 1-based region strings
	'''
	bed = BedTool(bedfilename)
	return ['{}:{}-{}'.format(interval.chrom, interval.start+1, interval.stop) for interval in bed]


def load_vcf(vcffilename, conf_regions):
	regionstr = ','.join(conf_regions)
	vcf = pysam.VariantFile(vcffilename, drop_samples = True)
	d = dict()
	for record in vcf.fetch(region = regionstr):
		d.setdefault(record.chrom, list()).append(record.pos)
	return d

def get_ref_base(reffile, chrom, pos):
	regionstr = '{0}:{1}-{1}'.format(chrom,pos)
	return reffile.fetch(region=regionstr)

def get_kmers_covering(read, pos, ksize):
	return [read[start:start+ksize] for start in range(0,len(read)-ksize+1) if pos >= start and pos < start+ksize]

def count_from_plp(pileups, refbase, errortable):
	for read in pileups:
		if read.is_del or read.is_refskip or read.alignment.query_sequence[read.query_position] != refbase:
			fastq.write('@{}\n{}\n+\n{}\n'.format(read.query_name, read.query_sequence, read.query_qualities))
			for k in get_kmers_covering(read.alignment.query_sequence, read.query_position, errortable.ksize()):
				errortable.count(k) 
	return errortable, fastq

#for now ignore sites in the vcf
#if there are two errors in a kmer, it will be counted twice.
def count_erroneous_kmers(samfile, ref, conf_regions, vcf, errortable):
	#try this and see if it works
	regionstr = ','.join(conf_regions)
	fastq = io.StringIO()
	# for regionstr in get_confident_regions(conf_regions):
	for col in samfile.pileup(region = regionstr, truncate = True):
		if col.reference_name in vcf and col.pos in vcf[col.reference_name]:
				continue
		else:
			refbase = get_ref_base(ref, col.reference_name, col.pos)
			errortable, fastq = count_from_plp(col.pileups, refbase, errortable, fastq)
	return errortable, fastq

def count_all(samfile, conf_regions, alltable):
	regionstr = ','.join(conf_regions)
	fastq = io.StringIO()
	for read in samfile.fetch(region=regionstr):
		alltable.consume(read.query_sequence)
		fastq.write('@{}\n{}\n+\n{}\n'.format(read.query_name, read.query_sequence, read.query_qualities))
	return alltable, fastq

def newinfo(*kwargs):
	return

def get_abundances(fastq,htable):
	ksize = htable.ksize()
	hashsizes = htable.hashsizes()
	tracking = khmer.Nodegraph(ksize, 1, 1, primes=hashsizes)
	return htable.abundance_distribution_with_reads_parser(khmer.ReadParser(fastq),tracking)

def print_abundances(abundances):
	pass

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

	khmer.khmer_args.info = newinfo
	args = khmer.khmer_args.build_counting_args().parse_args()
	alltable = khmer.khmer_args.create_countgraph(args, ksize=3)
	errortable = khmer.khmer_args.create_countgraph(args, ksize=3)

	samfile = pysam.AlignmentFile(samfilename)
	reffile = pysam.FastaFile(fafilename)
	conf_regions = get_confident_regions(bedfilename)
	vcf = load_vcf(vcffilename, conf_regions)
	alltable, allfq = count_all(samfile, conf_regions, alltable)
	errortable, errorfq = count_erroneous_kmers(samfile, reffile, conf_regions, vcf, errortable)
	allhashsizes = alltable.hashsizes()
	errorhashsizes = errortable.hashsizes()




if __name__ == '__main__':
    main()
