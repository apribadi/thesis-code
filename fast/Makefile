OPTS=-O1 -larmadillo -lgsl -lgslcblas

all: test

test: test.o rbm.o
	g++ test.o rbm.o -o test -larmadillo -lgsl -lgslcblas

rbm.o: rbm.cpp
	g++ -c rbm.cpp -O1

test.o: test.cpp
	g++ -c test.cpp -O1

clean:
	rm *.o test
