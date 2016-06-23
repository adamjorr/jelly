/*
A stand-alone C++ implementation of:
	https://github.com/adamjorr/khmer/blob/master/sandbox/slice-paired-reads-by-coverage.py
Support import/export of a jellyfish graph
Specify a min and/or max coverage to filter with
Support Interleaved, and 2-file PE reads
Support singleton output for interleaved/PE reads

If no output is specified, output interleaved to stdout
If no input is specified, assume interleaved input from stdin
If no graph is specified, create one and use it
If graph is specified, but file does not exist, create one, output to file, and use it
						if file does exist, use it
Ensure min or max is specified unless g is specified (ie using to solely output a graph).



Load inputs
Build (or load) Jellyfish graph
	graph g;
	if ! args.graph:
		g = build_graph()
	else:
		if -e args.graph:
			g = load_graph()
		else:
			g = build_graph()
			g.print_graph(args.graph)
iterate through reads (in graph if possible, otherwise in file if provided)
Spit to file if read is good



notable jellyfish code:
	sub_commands/count_main.cc#248-253 //code for writing hash to output
	typedef jellyfish::cooperative::hash_counter<jellyfish::mer_dna> mer_hash //jellyfish.hpp



the khmer function for doing this is:

void Hashtable::get_median_count(const std::string &s,
                                 BoundedCounterType &median,
                                 float &average,
                                 float &stddev)
{
    std::vector<BoundedCounterType> counts;
    this->get_kmer_counts(s, counts);

    if (!counts.size()) {
        throw khmer_exception("no k-mer counts for this string; too short?");
    }

    average = 0;
    for (std::vector<BoundedCounterType>::const_iterator i = counts.begin();
            i != counts.end(); ++i) {
        average += *i;
    }
    average /= float(counts.size());

    stddev = 0;
    for (std::vector<BoundedCounterType>::const_iterator i = counts.begin();
            i != counts.end(); ++i) {
        stddev += (float(*i) - average) * (float(*i) - average);
    }
    stddev /= float(counts.size());
    stddev = sqrt(stddev);

    sort(counts.begin(), counts.end());
    median = counts[counts.size() / 2]; // rounds down
}

void Hashtable::get_kmer_counts(const std::string &s,
                                std::vector<BoundedCounterType> &counts) const
{
    KmerIterator kmers(s.c_str(), _ksize);

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();
        BoundedCounterType c = this->get_count(kmer);
        counts.push_back(c);
    }
}

*/

#include <cstring>
#include <iostream>
#include <getopt.h>
#include <typeinfo>
#include <cstdlib>
#include "argparser.h"

int main (int argc, char **argv){

	ArgParser::arg_struc args = ArgParser::parse_args(argc,argv);

	std::cout << args.input_file << "\n";

	exit(EXIT_SUCCESS);
}















