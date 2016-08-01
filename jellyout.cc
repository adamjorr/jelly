
#include "jellyfish.hpp"
#include "jellyfish/merge_files.hpp"

Jellyout::Jellyout(std::string out_filename, int num_threads) :
	out_filename(out_filename),
	num_threads(num_threads) {};

Jellyout::print_hash(mer_hash ary, int length){
	jellyfish::file_header header;
	std::unique_ptr<jellyfish::dumper_t<mer_array>> dumper;

	header.fill_standard();
	header.set_cmdline(argc, argv);
	dumper.reset(new binary_dumper(length, ary.key_len(), num_threads, out_filename, &header));
	ary.dumper(dumper.get());
	dumper->dump(ary.ary());
	std::vector<const char*> files = dumper->file_names_cstr();
	try {
		merge_files(files, args.output_arg, header, 0, std::numeric_limits<uint64_t>::max());
	} catch(MergeError e) {
		jellyfish::err::die(jellyfish::err::msg() << e.what());
	}
	for(int i =0; i < dumper->nb_files(); ++i){
		unlink(files[i]);
    }
}




