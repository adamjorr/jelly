#ifndef __JELLY_IN_H_INCLUDED__
#define __JELLY_IN_H_INCLUDED__

#include <htslib/kseq.h>
#include <zlib.h>
#include <iostream>
#include <memory>
#include <jellyfish/jellyfish.hpp>

KSEQ_INIT(gzFile, gzread);

//Provides a "next" function to iterate over a fastq file
//next() returns -1 when EOF is reached.
class Jellyit{
	std::string filename;
	gzFile fp;
	kseq_t* seq;

	public:
		Jellyit(std::string filename);
		~Jellyit();
		kseq_t* next();
};

//To read a fastq, use static method read_fastq(file). EOF = -1
//To read a hash, use the static method read_hash(file).
class Jellyin{
	public:
		static std::unique_ptr<binary_reader> read_hash(std::string fname); //called once, returning an iterator over the hash
		static std::unique_ptr<Jellyit> read_fastq(std::string fname); //called once, returning an iterator over the fastq
		static std::unique_ptr<mer_hash> get_hash(std::string fname);//return a hash from FASTQ file(s)

};



#endif