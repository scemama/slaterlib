VERSION=1.0.0

CC=gcc

DEBUG_FLAGS=-Wall -Werror

CFLAGS=-g -fPIC $(DEBUG_FLAGS)

OBJ=core.o simple.o


libslater.so: $(OBJ)
	$(CC) $(CFLAGS) *.o -fPIC -shared -Wl,-soname,libslater.so.1 -o libslater.so.$(VERSION) -lc
	rm -f libslater.so ; ln -s libslater.so.$(VERSION) libslater.so
	rm -f libslater.so.1 ; ln -s libslater.so.$(VERSION) libslater.so.1

%.o: %.c slater_condon.h
	$(CC) $(CFLAGS) $< -c

test: test.c libslater.so
	$(CC) $(CFLAGS) test.c -lslater -L. -o test -Wl,-rpath=.

.PHONY: clean

clean: 
	rm -f libslater.so* *.o
