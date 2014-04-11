# Makefile for VolView plugins
#
# David Worth STFC May 2012

CC=gcc

VOLVIEW=/home/djw/CCPI/VolView-3.4-Linux-x86_64

CFLAGS=-I$(VOLVIEW)/include/KWVolView/Plugins/ -fPIC -Wall -std=c99 -g

LDFLAGS=-shared

vvQuan3DC.so : vvQuan3D.o vvHelper.o
	$(CC) $(LDFLAGS) -Wl,-soname,$@ -o $@ $^

vvQuan3D.o : vvQuan3D.c
	$(CC) $(CFLAGS) -c $<

vvHelper.o : vvHelper.c
	$(CC) $(CFLAGS) -c $<

clean :
	rm -f *.o *.so*
