
CXXFLAGS = -O2 -g -ggdb -pedantic -Wall -I../Include
OBJS=harmonic_test.o harmonics.o math_comb.o
RM=rm
RMFLAGS=-f
LD=$(CXX)

all : harmonic_test

harmonic_test : $(OBJS)
	$(LD) -o $@ $^
clean:
	$(RM) $(RMFLAGS) $(OBJS)
