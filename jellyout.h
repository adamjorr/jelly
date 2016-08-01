#ifndef __JELLY_OUT_H_INCLUDED__
#define __JELLY_OUT_H_INCLUDED__

#include <cstring>

//https://github.com/gmarcais/Jellyfish/blob/master/sub_commands/count_main.cc#L252

class Jellyout{
	std::string out_filename;
	int num_threads;

	public:
		Jellyout(std::string out_filename, int num_threads);
		void print_reads();
		void print_hash();

};


#endif