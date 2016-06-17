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
iterate through reads (in graph if possible, otherwise in file if provided)
Spit to file if read is good




the kmer function for this is:

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
#include "argparser.h"

int main (int argc, char **argv){

	// char graph_file[];
	// char input_file[];
	// char input1_file[];
	// char input2_file[];
	// char out_file[];
	// char out1_file[];
	// char out2_file[];
	// char singletons_file[];
	// unsigned long min;
	// unsigned long max;


	// int option_index = 0;
	// static struct option long_options[] = {
	// 	//these have a short option
	// 	{"graph", required_argument, 0, 'g'},
	// 	{"input", required_argument, 0, 'i'},
	// 	{"input1", required_argument, 0, '1'},
	// 	{"input2", required_argument, 0, '2'},
	// 	{"out", required_argument, 0, 'o'},
	// 	{"out1", required_argument, 0, 'p'},
	// 	{"out2", required_argument, 0, 'q'},
	// 	{"singletons", required_argument, 0, 's'},
	// 	//these do not have a short option
	// 	{"min", required_argument, 0, 'a'},
	// 	{"max", required_argument, 0, 'b'},
		
	// 	{0,0,0,0}
	// };

	// const char * shortopts = "g:i:1:2:o:p:q:s:";
	// int c;
	// while (c = getopt_long_only(argc, argv, shortopts, long_options, &option_index) != -1){
	// 	switch(c){
	// 		case 'g':
	// 			graph_file = optarg;
	// 			break;
	// 		case 'i':
	// 			input_file = optarg;
	// 			break;
	// 		case '1':
	// 			input1_file = optarg;
	// 			break;
	// 		case '2':
	// 			input2_file = optarg;
	// 			break;
	// 		case 'o':
	// 			out_file = optarg;
	// 			break;
	// 		case 'p':
	// 			out1_file = optarg;
	// 			break;
	// 		case 'q':
	// 			out2_file = optarg;
	// 			break;
	// 		case 's':
	// 			singletons_file = optarg;
	// 			break;
	// 		case 'a':
	// 			min = strtoul(optarg,0,0);
	// 			break;
	// 		case 'b':
	// 			max = strtoul(optarg,0,0);
	// 			break;
	// 		case '?':
	// 			exit(EXIT_FAILURE);
	// 			break;
	// 		default:
	// 			std::cerr << "Unable to properly parse options. Aborting . . .";
	// 			exit(EXIT_FAILURE);
	// 	}

	// }

	ArgParser p;
	char** args = p.parse_args(argc,argv);

	// for(char* z : args){
	// 	std::cout << z << "\t" << *z << "\n";
	// }

	exit(0);

}















