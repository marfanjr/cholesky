all: main

lib-aux.o: lib-aux.c lib-aux.h 
	gcc -c lib-aux.c

cholesky.o: cholesky.c cholesky.h lib-aux.h
	gcc -c cholesky.c

main.o: main.c lib-aux.h cholesky.h
	gcc -c main.c

main: main.o lib-aux.o cholesky.o
	gcc -o invmat main.o lib-aux.o cholesky.o -lm

clean: 
	rm -rf *.o
