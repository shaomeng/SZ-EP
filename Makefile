SHELL = /bin/sh

FLAGS = -Wall -O2 -std=c++17

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
  CXX = g++-11
  LIBSZ = ./SZ-install/lib/
  RPATH = -Wl,-rpath,${LIBSZ}
endif
ifeq ($(UNAME), Linux)
  CXX = g++
  LIBSZ = ./SZ-install/lib/libSZ.so
  RPATH = -Wl,-rpath=${LIBSZ}
endif

all: sz_ep

EP.o: EP.cpp 
	${CXX} ${FLAGS} -c -o $@ $<

sz_ep: main.cpp  EP.o
	${CXX} ${FLAGS} -fopenmp -o $@ $^ -lSZ -L${LIBSZ} ${RPATH}

clean:
	rm -f EP.o sz_ep
