#ifndef __JELLY_OUT_H_INCLUDED__
#define __JELLY_OUT_H_INCLUDED__

#include <cstring>
#include <kseq.h>

//https://github.com/gmarcais/Jellyfish/blob/master/sub_commands/count_main.cc#L252

class Jellyout{
	std::string out_filename;
	int num_threads;

	public:
		Jellyout(std::string out_filename, int num_threads);
		void print_read(const kseq_t& seq);
		void print_hash(const mer_hash& ary, int length = 4);

};

ostream &operator<< (ostream &output, const kseq_t& seq);

#endif