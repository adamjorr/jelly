
#include "jellyin.h"
#include <htslib/kseq.h>
#include <zlib.h>
#include <iostream>

/*
class Jellyin{
	std::string filename;
	gzFile fp;
	kseq_t* seq;

	public:
		Jellyin();
		~Jellyin();
		kseq_t* read_fastq();
		mer_hash& read_hash();
};
*/

//TODO: Figure out what you need to do to get fastq / hash reading to work. Maybe split into 2 classes? Use a "mode" option?


Jellyin::Jellyin(std::string filename) :
	filename(filename) {
		if (filename != "0"){
			fp = gzopen(filename);
			seq = kseq_init(fp);
		}
	}

Jellyin::~Jellyin(){
	if (filename != "0"){
		kseq_destroy(seq);
		gzclose(fp);
	}
}

//called repeatedly, each time returning a record from a fastq.
kseq_t* read_fastq(){
	kseq_read(seq);
	return seq;
}

//called once, returning the entire hash
mer_hash& read_hash(){
	std::ifstream is(args.db_arg);
	jellyfish::file_header header;
	header.read(is);
	jellyfish::mer_dna::k(header.key_len() / 2);
	binary_reader reader(is, &header);
}




