jelly : jelly.o argparser.o jellyout.o jellyin.o
	g++ --std=c++11 jelly.o argparser.o jellyout.o jellyin.o -L/usr/local/lib -ljellyfish-2.0 -lhts -lz -lpthread -Wall -Wextra -pedantic -o jelly

argparser.o : argparser.h argparser.cc
	g++ --std=c++11 -c argparser.cc -Wall -Wextra -pedantic

jellyout.o : jellyout.h jellyout.cc
	g++ --std=c++11 -c jellyout.cc -Wall -Wextra -pedantic

jellyin.o : jellyin.h jellyin.cc
	g++ --std=c++11 -c jellyin.cc -Wall -Wextra -pedantic

kmer_counter.o : kmer_counter.h kmer_counter.cc
	g++ --std=c++11 -c jellyin.cc -Wall -Wextra -pedantic
