#ifndef __JELLY_IN_H_INCLUDED__
#define __JELLY_IN_H_INCLUDED__

#include <htslib/kseq.h>
#include <zlib.h>
#include <iostream>
#include <memory>
#include <jellyfish/jellyfish.hpp>

KSEQ_INIT(gzFile, gzread)

//To read a fastq, create an object and call read_fastq(). EOF = -1
//To read a hash, simply use the static method read_hash(file).
class Jellyin{
	std::string filename;
	gzFile fp;
	kseq_t* seq;

	public:
		Jellyin(std::string filename);
		~Jellyin();
		kseq_t* next(); //called repeatedly, each time returning a record from a fastq.
		static std::unique_ptr<binary_reader> read_hash(std::string fname); //called repeatedly, returning an iterator over the hash
};

#endif