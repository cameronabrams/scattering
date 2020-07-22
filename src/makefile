libskdataspace.so:    skdataspace.o skdataspace_wrap.o
	gcc -O3 -shared -o $@ skdataspace.o skdataspace_wrap.o -lm

skdataspace.o:		skdataspace.c skdataspace.h
	gcc -O3 -c -fpic -o $@ skdataspace.c

skdataspace_wrap.o: 	skdataspace_wrap.c
	gcc -O3 -c -fpic -o $@ skdataspace_wrap.c

skdataspace_wrap.c:	skdataspace.c skdataspace.h skdataspace.i
	swig -tcl8 skdataspace.i

