#ifndef __JELLY_IN_H_INCLUDED__
#define __JELLY_IN_H_INCLUDED__

#include <htslib/kseq.h>
#include <zlib.h>
#include <iostream>

KSEQ_INIT(gzFile,gzread)

class Jellyin{
	std::string filename;
	gzFile fp;
	kseq_t* seq;

	public:
		Jellyin(std::string filename = "0");
		~Jellyin();
		kseq_t* read_fastq(); //called repeatedly, each time returning a record from a fastq.
		mer_hash& read_hash(); //called once, returning the entire hash.
};

#endif