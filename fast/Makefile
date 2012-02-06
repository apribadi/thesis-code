
COPTS=-O2 -Wall -Wextra
LOPTS=
all: test

test: test.o rbm.o
	g++ -o test test.o rbm.o -larmadillo -lgsl -lgslcblas $(LOPTS)

.cpp.o:
	g++ -c $(COPTS) $<

clean:
	rm -f *.o test
