
#include "parser.h"
#include <getopt.h>
#include <iostream>

char** ArgParser::parse_args(int argc, char **argv){

	char *graph_file;
	char *input_file;
	char *input1_file;
	char *input2_file;
	char *out_file;
	char *out1_file;
	char *out2_file;
	char *singletons_file;
	char *min;
	char *max;
	char *to_return[] = {graph_file, input_file, input1_file, input2_file, out_file, out1_file, out2_file, singletons_file, min, max};
	char **ret = to_return;

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
		//these do not have a short option
		{"min", required_argument, 0, 'a'},
		{"max", required_argument, 0, 'b'},
		
		{0,0,0,0}
	};

	const char * shortopts = "g:i:1:2:o:p:q:s:";
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


	// return to_return;
	return ret;

}
