CC=g++
CFLAGS=-g -Wall
IFLAGS=-I/usr/X11R6/include
LFLAGS=-L/usr/X11R6/lib -lX11
all: rmn

rmn: rmn.cc points.cc vector.cc graphics.cc polynomial.cc roots.cc rational_map.cc partial_fraction.cc benchmark.cc
	$(CC) $(CFLAGS) $(IFLAGS) -o rmn rmn.cc $(LFLAGS) -lm

clean: rm rmn
