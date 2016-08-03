#ifndef __JELLY_OUT_H_INCLUDED__
#define __JELLY_OUT_H_INCLUDED__

#include <cstring>
#include <htslib/kseq.h>
#include <zlib.h>
#include <iostream>
#include <jellyfish/jellyfish.hpp>


//https://github.com/gmarcais/Jellyfish/blob/master/sub_commands/count_main.cc#L252

KSEQ_INIT(gzFile,gzread)

class Jellyout{
	std::string out_filename;
	int num_threads;
	std::ostream* fp;
	std::ofstream fout;

	public:
		Jellyout(std::string out_filename, int num_threads);
		~Jellyout();
		void print_read(const kseq_t* seq);
		void print_hash(mer_hash& ary, int length = 4);

};

std::ostream &operator<< (std::ostream &output, const kseq_t* seq);

#endif