
#include "hash.h"

mer_counter::mer_counter(int nb_threads, mer_hash_type& mer_hash, char** file_begin, char** file_end, bool canonical):
	mer_hash_(mer_hash),
	streams_(file_begin,file_end),
	parser_(jellyfish::mer_dna::k(),streams_.nb_streams(), 3 * nb_threads, 4096, streams_),
	canonical_(canonical){}

virtual void mer_counter::start(int thid){
	mer_iterator_type mer_iterator(parser_, canonical_);
	for( ; mers; ++mers){
		mer_hash_.add(*mers,1);
	}
	mer_hash_.done();
}

