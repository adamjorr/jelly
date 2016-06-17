
#ifndef __ARGPARSER_H_INCLUDED__
#define __ARGPARSER_H_INCLUDED__

#include <string>

class ArgParser{
	public:
		struct arg_struc {
			std::string graph_file;
			std::string input_file;
			std::string input1_file;
			std::string input2_file;
			std::string out_file;
			std::string out1_file;
			std::string out2_file;
			std::string singletons_file;
			std::string min;
			std::string max;
		};
		static arg_struc parse_args(int argc, char **argv);

};

#endif
