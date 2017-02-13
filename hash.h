#ifndef __HASH_H_INCLUDED__
#define __HASH_H_INCLUDED__

#include <jellyfish>

typedef jellyfish::cooperative::hash_counter<jellyfish::mer_dna>                  mer_hash_type;
typedef jellyfish::mer_overlap_sequence_parser<jellyfish::stream_manager<char**>> sequence_parser_type;
typedef jellyfish::mer_iterator<sequence_parser_type, jellyfish::mer_dna>         mer_iterator_type;

class mer_counter : public jellyfish::thread_exec {
	mer_hash_type& mer_hash_;
	jellyfish::stream_manager<char**> streams_;
	sequence_parser_type parser_;
	const bool canonical_;

	public:
		mer_counter(int nb_threads, mer_hash_type& mer_hash,
		          char** file_begin, char** file_end,
		          bool canonical)
		virtual void start(int thid) 
};

#endif
