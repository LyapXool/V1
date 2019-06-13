CC=g++
CFLAGS=-I/ADDRES/TO/armadillo/include -std=c++11 -fopenmp -march=native -Wall -Wextra -m64 -DARMA_64BIT_WORD  -O2
LFLAGS=-L/opt/share/armadillopp/armadillo-9.200.4/ -larmadillo -L/opt/ohpc/pub/libs/gnu/openblas/0.2.19/lib/ -lopenblas -llapack -lgfortran
DEPS=Makefile chainrecurrentsets.hpp generalities.hpp instructions.hpp odesystem.hpp RBF.hpp wendland.hpp
OBJ=chainrecurrentsets.o generalities.o instructions.o odesystem.o problem.o odesystem.o RBF.o wendland.o


%.o: %.cpp $(DEPS)
        $(CC) -c $(CFLAGS) -o $@ $<  

solver: $(OBJ)
        $(CC) $(CFLAGS) $(LFLAGS) -o $@ $^ 

.PHONY: clean

clean:
        rm -f *.o
        rm -f *.m
        


