SHELL = /bin/sh

FLAGS = -Wall -O2 -std=c++17

CXX = g++-11

all: sz_ep

EP.o: EP.cpp 
	${CXX} ${FLAGS} -c -o $@ $<

sz_ep: main.cpp  EP.o
	${CXX} ${FLAGS} -fopenmp -o $@ $^

clean:
	rm -f EP.o sz_ep
