
#include "jellyfish.hpp"
#include "jellyfish/merge_files.hpp"
#include <kseq.h>
#include <zlib.h>
#include <iostream>

Jellyout::Jellyout(std::string out_filename, int num_threads) :
	out_filename(out_filename),
	num_threads(num_threads) {};

Jellyout::print_hash(const mer_hash& ary, int length = 4){
	jellyfish::file_header header;
	std::unique_ptr<jellyfish::dumper_t<mer_array>> dumper;

	header.fill_standard();
	header.set_cmdline(argc, argv);
	dumper.reset(new binary_dumper(length, ary.key_len(), num_threads, out_filename, &header));
	ary.dumper(dumper.get());
	dumper->dump(ary.ary());
	std::vector<const char*> files = dumper->file_names_cstr();
	try {
		jellyfish::merge_files(files, args.output_arg, header, 0, std::numeric_limits<uint64_t>::max());
	} catch(MergeError e) {
		jellyfish::err::die(jellyfish::err::msg() << e.what());
	}
	for(int i = 0; i < dumper->nb_files(); ++i){
		unlink(files[i]);
    }
}

//optimization here? will c++ flush the buffer when it's full automatically?
//that would make << "\n" better than << std::endl
ostream &operator<< (ostream &output, const kseq_t& seq){
	output << seq->name.s << std::endl
	output << seq->seq.s << std::endl
	output << + << std::endl
	output << seq->qual.s << std::endl

}

Jellyout::print_read(const kseq_t& seq){
	
}

