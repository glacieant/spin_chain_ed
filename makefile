# makefile for exact diagonalisation of random heisenberg spin chain

OBJ = randspin.o randomham.o randomgenerate.o randspin_utility.o randstate.o mt19937ar.o

CC = gcc
LFLAGS = -lm -llapack -lblas

randspin : $(OBJ)
	   $(CC) -o randspin $(OBJ) $(LFLAGS) 

randspin.o : randspin.c
	     $(CC) -c randspin.c

randomham.o : randomham.c
	      $(CC) -c randomham.c

randomgenerate.o : randomgenerate.c
	   	   $(CC) -c randomgenerate.c

randspin_utility.o : randspin_utility.c
		     $(CC) -c randspin_utility.c

mt19937ar.o : mt19937ar.c
	      $(CC) -c mt19937ar.c

randstate.o : randstate.c
	      $(CC) -c randstate.c

clean :
	rm randspin $(OBJ)
