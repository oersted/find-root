all: findroot.o function.o norm.o
	gcc -Wall findroot.o function.o norm.o -o findroot -lgsl -lgslcblas -lm -O2
findroot.o:
	gcc -Wall -c findroot.c
function.o:
	gcc -Wall -c function.c
norm.o:
	gcc -Wall -c norm.c
clean:
	rm *.o findroot
