SHELL=/bin/sh
include ../config.h
objects=util.o math.o transform.o mod_global.o cluster.o  tram.o coarse_grain.o markov.o fileio.o tpt.o pharsecmd.o msm.o scc.o

cpp_obj = scc.o
FC=mpif90 -g
FFLAGS= -fPIC
FOPTFLAGS= -O3 -mtune=native
FPP=cpp -DMPI
FWARNFLAGS=-Wno-unused-function
LAPACKLIBS=/home/haomiao/server_home/tmp/fsatool/lib/lapack-3.7.0/lib64
FLIBS=-llapack -lblas
CPPLIBS = -lstdc++

install: msm.x

msm.x : $(objects)
	$(FC) $(FOPTFLAGS) $(FWARNFLAGS) -o $@ $^ -L$(LAPACKLIBS) $(FLIBS) $(CPPLIBS)


%.o : %.f90
	$(FC) $(FOPTFLAGS) -c $^

scc.o : scc.cpp
	g++ $(FOPTFLAGS) -c $^

clean:
	rm *.o *.mod *.x
