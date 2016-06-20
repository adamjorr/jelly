
#include "argparser.h"
#include <getopt.h>
#include <iostream>
#include <string>

ArgParser::arg_struc ArgParser::parse_args(int argc, char **argv){

	std::string graph_file;
	std::string input_file;
	std::string input1_file;
	std::string input2_file;
	std::string out_file;
	std::string out1_file;
	std::string out2_file;
	std::string singletons_file;
	std::string threads;

	std::string min;
	std::string max;
	arg_struc to_return;

	int option_index = 0;
	static struct option long_options[] = {
		//these have a short option
		{"graph", required_argument, 0, 'g'},
		{"input", required_argument, 0, 'i'},
		{"input1", required_argument, 0, '1'},
		{"input2", required_argument, 0, '2'},
		{"out", required_argument, 0, 'o'},
		{"out1", required_argument, 0, 'p'},
		{"out2", required_argument, 0, 'q'},
		{"singletons", required_argument, 0, 's'},
		{"threads", required_argument, 0, 't'},
		//these do not have a short option
		{"min", required_argument, 0, 'a'},
		{"max", required_argument, 0, 'b'},
		
		{0,0,0,0}
	};

	const char * shortopts = "g:i:1:2:o:p:q:s:t:";
	int c = 0;
	while (c != -1){
		c = getopt_long_only(argc, argv, shortopts, long_options, &option_index);
		switch(c){
			case 'g':
				graph_file = optarg;
				break;
			case 'i':
				input_file = optarg;
				break;
			case '1':
				input1_file = optarg;
				break;
			case '2':
				input2_file = optarg;
				break;
			case 'o':
				out_file = optarg;
				break;
			case 'p':
				out1_file = optarg;
				break;
			case 'q':
				out2_file = optarg;
				break;
			case 's':
				singletons_file = optarg;
				break;
			case 'a':
				min = optarg;
				break;
			case 'b':
				max = optarg;
				break;
			case 't':
				threads = optarg;
				break;
			case '?':
				exit(EXIT_FAILURE);
				break;
			case -1:
				break;
			default:
				std::cerr << "Unable to properly parse options. Aborting . . .\n";
				exit(EXIT_FAILURE);
		}

	}

	to_return = {graph_file, input_file, input1_file, input2_file, out_file, out1_file, out2_file, singletons_file, threads, min, max};

	return to_return;

}
