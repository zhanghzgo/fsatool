FCC=mpif90
include ../config.h
objects = moldata.o colvar.o eds.o smd.o flood.o temd.o

fsatool_lib : $(objects)
	ar -crusv fsatool_lib.a $(objects)

%.o : %.f90
	$(FCC) -c $^ $(FOPTFLAGS) $(NETCDFINC) -I$(AMBER_SOURCE)/src/pmemd/src

clean:
	rm *.o *.mod *.a
