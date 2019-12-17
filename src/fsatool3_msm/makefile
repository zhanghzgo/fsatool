SHELL=/bin/sh
objects=util.o math.o netcdf_helper.o mpi_shared_memory.o transform.o mod_global.o cluster.o  tram.o coarse_grain.o markov.o fileio.o tpt.o pharsecmd.o msm.o scc.o

cpp_obj = scc.o
FC=mpif90 -g
FFLAGS= -fPIC
FOPTFLAGS= -O3 -mtune=native
FPP=cpp -DMPI
FWARNFLAGS=-Wno-unused-function
LAPACKLIBS=../lib/lapack-3.7.0/lib64
NETCDFLIBS=-lnetcdf -lnetcdff
FLIBS=-llapack -lblas
CPPLIBS = -lstdc++

install: msm.x

msm.x : $(objects)
	$(FC) $(FOPTFLAGS) $(FWARNFLAGS) -o $@ $^ ${AMBERHOME}/src/pmemd/src/AmberNetcdf.o  ${AMBERHOME}/lib/libnetcdf.so ${AMBERHOME}/lib/libnetcdff.so -L$(LAPACKLIBS) $(FLIBS) $(CPPLIBS)

%.o : %.F90
	$(FC) $(FOPTFLAGS) -c $^ -I${AMBERHOME}/include -I${AMBERHOME}/src/pmemd/src

scc.o : scc.cpp
	g++ $(FOPTFLAGS) -c $^

clean:
	rm *.o *.mod *.x
