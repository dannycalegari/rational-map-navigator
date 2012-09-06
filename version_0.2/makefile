CC=g++
CFLAGS=-g -Wall
IFLAGS=-I/usr/X11R6/include
LFLAGS=-L/usr/X11R6/lib -lX11
all: rmn

rmn: rmn.cc braid.cc points.cc vector.cc graphics.cc polynomial.cc roots.cc rational_map.cc draw.cc insert.cc magnify.cc read_write.cc user_interface.cc
	$(CC) $(CFLAGS) $(IFLAGS) -o rmn rmn.cc $(LFLAGS) -lm

clean: rm rmn
