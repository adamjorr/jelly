
#include "jellyin.h"
#include <htslib/kseq.h>
#include <zlib.h>
#include <iostream>
#include <memory>

//called once, returning an iterator
std::unique_ptr<binary_reader> Jellyin::read_hash(std::string fname){
	std::unique_ptr<std::ifstream> is(new std::ifstream(fname));
	if (! is->good()){
		throw std::ios_base::failure("Failure opening " + fname);
	}
	jellyfish::file_header header;
	header.read(*is);
	std::unique_ptr<binary_reader> reader(new binary_reader(*is,&header));
	return reader;
}

//called once, returning an iterator over the fastq
std::unique_ptr<Jellyit> Jellyin::read_fastq(std::string fname){
	std::unique_ptr<Jellyit> iterator(new Jellyit(fname));
	return iterator;
}


Jellyit::Jellyit(std::string filename) :
	filename(filename) {
		fp = gzopen(filename.c_str(),"r");
		if (fp == NULL){
			throw std::ios_base::failure("Failure opening " + filename);
		}
		seq = kseq_init(fp);
	}

Jellyit::~Jellyit(){
	kseq_destroy(seq);
	gzclose(fp);
}

kseq_t* Jellyit::next(){
	kseq_read(seq);
	return seq;
}



