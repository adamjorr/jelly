#include "kmer_counter.h"

Kmer_counter::kmer_counter(int num_threads, mer_hash& ary, ParserType& parser) : 
	num_threads(num_threads),
	ary(ary),
	parser(parser)
	{}

void Kmer_counter::start(int thid){
	MerIteratorType it(parser,0);

	for (; it; ++it){
		ary.add(*it,1);
	}
	ary.done();
}
