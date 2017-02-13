
#include "jellyin.h"
#include "kmer_counter.h"
#include <htslib/kseq.h>
#include <zlib.h>
#include <iostream>
#include <memory>
#include <jellyfish/hash_counter.hpp>
#include <thread>
#include <mutex>

std::mutex mtx;

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
std::unique_ptr<mer_hash> Jellyin::create_hash(std::vector<const char*> fname){
	jellyfish::mer_dna::k(25); // Set length of mers (k=25)
	const uint64_t hash_size    = 10000000; // Initial size of hash.
	const uint32_t num_reprobes = 126;
	const uint32_t num_threads  = 16; // Number of concurrent threads
	const uint32_t counter_len  = 7;  // Minimum length of counting field
	const bool     canonical    = true; // Use canonical representation
	mer_hash hash(new mer_hash(hash_size,jellyfish::mer_dna::k()*2,counter_len,num_threads,num_reprobes));
	// kmer_qual_counter counter(jellyfish::mer_dna::k,num_threads,3 * num_threads, 4096, hash, new StreamIterator)
	jellyfish::stream_manager<std::vector<const char*>::const_iterator> streams(fname.front(),fname.back());
	sequence_qual_parser parser(num_threads * 3, 100, streams.nb_streams(), &streams);
	Kmer_counter counter(num_threads, &hash, &parser);

}

//multithreaded function to add to a hash
void add_to_hash(mer_hash* hash, binary_reader* reader){
	mtx.lock();
	if (reader->next()){
		auto k = new reader->key(); //make sure we don't keep k and v around
		auto v = new reader->val(); //since we're doing this recursively
		mtx.unlock();
		hash->add(k,v);
		delete k;
		delete v;
		add_to_hash(hash, reader);
	} else {
		mtx.unlock();
		hash->done();
	}
}

//return a hash object from a hash file
std::unique_ptr<mer_hash> Jellyin::open_hash(std::string fname, uint32_t num_threads){
	jellyfish::mer_dna::k(25); // Set length of mers (k=25)
	const uint64_t hash_size    = 10000000; // Initial size of hash.
	const uint32_t num_reprobes = 126;
	const uint32_t counter_len  = 7;  // Minimum length of counting field
	const bool     canonical    = true; // Use canonical representation

	std::unique_ptr<mer_hash> hash(new mer_hash(hash_size,jellyfish::mer_dna::k()*2,counter_len,num_threads,num_reprobes));
	std::unique_ptr<binary_reader> in = read_hash(fname);

	// while (in->next()){
	// 	hash->add(in->key(),in->val());
	// }

	std::vector<std::thread> threads;
	for (int i = 0; i < num_threads; ++i){
		threads.push_back(std::thread(add_to_hash, hash, in));
		}

	for (std::thread t : threads){
		t.join();
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



