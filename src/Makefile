VERSION=1.0.0

CC=gcc

DEBUG_FLAGS=-Wall -Wpedantic -Werror

CFLAGS=-g -fPIC $(DEBUG_FLAGS)

libslater.so: slater_condon.o
	$(CC) $(CFLAGS) *.o -fPIC -shared -Wl,-soname,libslater.so.1 -o libslater.so.$(VERSION) -lc
	ln -s libslater.so.$(VERSION) libslater.so
	ln -s libslater.so.$(VERSION) libslater.so.1

%.o: %.c slater_condon.h
	$(CC) $(CFLAGS) slater_condon.c -c

test: test.c libslater.so
	$(CC) $(CFLAGS) test.c -lslater -L. -o test -Wl,-rpath=.

.PHONY: clean

clean: 
	rm -f libslater.so* *.o
