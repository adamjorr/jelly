
#include "jellyin.h"
#include <htslib/kseq.h>
#include <zlib.h>
#include <iostream>
#include <memory>
#include <jellyfish/hash_counter.hpp>
#include <thread>

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

//return a hash from FASTQ file(s)
std::unique_ptr<mer_hash> Jellyin::create_hash(std::string fname){
	jellyfish::mer_dna::k(25); // Set length of mers (k=25)
	const uint64_t hash_size    = 10000000; // Initial size of hash.
	const uint32_t num_reprobes = 126;
	const uint32_t num_threads  = 16; // Number of concurrent threads
	const uint32_t counter_len  = 7;  // Minimum length of counting field
	const bool     canonical    = true; // Use canonical representation
	std::unique_ptr<mer_hash> hash(new mer_hash(hash_size,jellyfish::mer_dna::k()*2,counter_len,num_threads,num_reprobes));
	mer_counter counter(num_threads,hash,)

}

//return a hash object from a hash file
//TODO: Need to make this multithreaded.
//put a lock on the binary reader and make each thread read?
std::unique_ptr<mer_hash> Jellyin::open_hash(std::string fname, uint32_t num_threads){
	jellyfish::mer_dna::k(25); // Set length of mers (k=25)
	const uint64_t hash_size    = 10000000; // Initial size of hash.
	const uint32_t num_reprobes = 126;
	const uint32_t counter_len  = 7;  // Minimum length of counting field
	const bool     canonical    = true; // Use canonical representation

	std::unique_ptr<mer_hash> hash(new mer_hash(hash_size,jellyfish::mer_dna::k()*2,counter_len,num_threads,num_reprobes));

	std::thread threads[num_threads - 1];
	for (int i = 0; i < num_threads - 1; ++i){
		threads[i] = std::thread()
	}


	
	std::unique_ptr<binary_reader> in = read_hash(fname);
	while (in->next()){
		hash->add(in->key(),in->val());
	}
	return hash;
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



