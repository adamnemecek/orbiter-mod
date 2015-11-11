SRC=../..
DISCRETA_CPP = /usr/bin/g++
CPPFLAGS = -Wall -I$(SRC)/LIB -I$(SRC)/LIB/DISCRETA
LIB = \
	$(SRC)/LIB/TOP_LEVEL/libtoplevel.a \
	$(SRC)/LIB/DISCRETA/libdiscreta.a \
	$(SRC)/LIB/SNAKES_AND_LADDERS/libsnakesandladders.a  \
	$(SRC)/LIB/ACTION/libaction.a \
	$(SRC)/LIB/INCIDENCE/libincidence.a \
	$(SRC)/LIB/GALOIS/libgalois.a 
LFLAGS = -lm



all: \
	blt.out 




blt.out:  blt.o blt_generator.o ;
	$(DISCRETA_CPP)  blt.o blt_generator.o $(CPPFLAGS) $(LIB) $(LFLAGS) -o $@


clean:
	- rm *.o *.a *.out *~

blt.o: blt.C blt.h; $(DISCRETA_CPP) $(CPPFLAGS) -c blt.C -o $@
blt_generator.o: blt_generator.C blt.h ; $(DISCRETA_CPP) $(CPPFLAGS) -c blt_generator.C -o $@




