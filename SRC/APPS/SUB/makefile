SRC=../../../../ORBITER/SRC
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



all: sub.out

sub.out:  sub.C sub_generator.o ;
	$(DISCRETA_CPP)  sub.C sub_generator.o $(CPPFLAGS) $(LIB) $(LFLAGS) -g -o $@



clean:
	- rm *.o *.a *.out *~



sub_generator.o: sub_generator.C; $(DISCRETA_CPP) $(CPPFLAGS) -c sub_generator.C -o $@
