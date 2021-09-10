SHELL = /bin/sh

FLAGS = -Wall -O2 -std=c++17

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
  CXX = g++-11
  LIBSZ_DIR = ./SZ-install/lib/
  RPATH = -Wl,-rpath,${LIBSZ_DIR}
endif
ifeq ($(UNAME), Linux)
  CXX = g++
  LIBSZ_DIR = ./SZ-2.1.12/build/sz/
  RPATH = -Wl,-rpath=${LIBSZ_DIR}
endif

all: sz_ep

EP.o: EP.cpp 
	${CXX} ${FLAGS} -c -o $@ $<

sz_ep: main.cpp  EP.o
	${CXX} ${FLAGS} -fopenmp -o $@ $^ -lSZ -L${LIBSZ_DIR} ${RPATH}

clean:
	rm -f EP.o sz_ep
