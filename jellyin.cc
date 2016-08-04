
#include "jellyin.h"
#include <htslib/kseq.h>
#include <zlib.h>
#include <iostream>
#include <memory>

Jellyin::Jellyin(std::string filename) :
	filename(filename) {
		fp = gzopen(filename.c_str(),"r");
		seq = kseq_init(fp);
	}

Jellyin::~Jellyin(){
	kseq_destroy(seq);
	gzclose(fp);
}

//called repeatedly, each time returning a record from a fastq.
kseq_t* Jellyin::next(){
	kseq_read(seq);
	return seq;
}

//called once, returning an iterator
std::unique_ptr<binary_reader> Jellyin::read_hash(std::string fname){
	std::unique_ptr<std::ifstream> is(new std::ifstream(fname));
	jellyfish::file_header header;
	header.read(*is);
	std::unique_ptr<binary_reader> reader(new binary_reader(*is,&header));
	return reader;
}




