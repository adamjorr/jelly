jelly : jelly.o argparser.o
	g++ --std=c++11 jelly.o argparser.o -L/usr/local/lib -ljellyfish-2.0 -o jelly

argparser.o : argparser.h argparser.cc
	g++ --std=c++11 argparser.cc
