jelly : jelly.o argparser.o jellyout.o
	g++ --std=c++11 jelly.o argparser.o -L/usr/local/lib -ljellyfish-2.0 -lhts -o jelly

argparser.o : argparser.h argparser.cc
	g++ --std=c++11 -c argparser.cc

jellyout.o : jellyout.h jellyout.cc
	g++ --std=c++11 -c jellyout.cc
