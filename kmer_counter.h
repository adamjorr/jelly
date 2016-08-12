#ifndef __KMER_COUNTER_H__
#define __KMER_COUNTER_H__

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <jellyfish/mer_qual_iterator.hpp>

template<typename MerIteratorType, typename ParserType>
class Kmer_counter : public jellyfish::thread_exec {
	int num_threads;
	mer_hash& ary;
	ParserType& parser;

	public:
		Kmer_counter(int num_threads, mer_hash& ary, ParserType& parser);
		virtual void start(int thid);

};

typedef public jellyfish::whole_sequence_parser<jellyfish::stream_manager<std::vector<const char*>::const_iterator>> sequence_qual_parser;
typedef kmer_counter<file_vector::const_iterator, mer_qual_iterator, sequence_qual_parser> kmer_qual_counter;

#endif