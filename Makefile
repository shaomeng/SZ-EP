SHELL = /bin/sh

FLAGS = -O2 -std=c++17

all: sz-ep

EP.o: EP.cpp 
	c++ ${FLAGS} -c -o $@ $<


