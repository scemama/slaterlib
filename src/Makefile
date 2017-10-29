VERSION=1.0.0

CC=gcc

DEBUG_FLAGS=-Wall -Werror 

CFLAGS=-g -fPIC $(DEBUG_FLAGS) -O3 -mavx

OBJ_LIB=core.o 
OBJ=simple.o test.o debug.o


libslater.so: $(OBJ_LIB) 
	$(CC) $(OBJ_LIB) -fPIC -shared -Wl,-soname,libslater.so.1 -o libslater.so.$(VERSION) -lc
	rm -f libslater.so ; ln -s libslater.so.$(VERSION) libslater.so
	rm -f libslater.so.1 ; ln -s libslater.so.$(VERSION) libslater.so.1

%.o: %.c slater_condon.h slater_condon_simple.h debug.h tags
	$(CC) $(CFLAGS) $< -c

test: $(OBJ) libslater.so
	$(CC) $(CFLAGS) $(OBJ) -lslater -L. -o test -Wl,-rpath=.

.PHONY: clean tags

clean: 
	rm -f libslater.so* *.o

tags: 
	ctags *.c *.h

