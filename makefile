CC=g++
CFLAGS=-g -Wall
IFLAGS=-I/usr/X11R6/include
LFLAGS=-L/usr/X11R6/lib -lX11
all: rational

rational: rational.cc linear.cc points.cc polynomial.cc rational_map.cc graphics.cc
	$(CC) $(CFLAGS) $(IFLAGS) -o rational rational.cc $(LFLAGS) -lm

clean: rm rational
