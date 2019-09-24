objects = netcdf_helper.o extract_info.o main.o
pmemd_src = ${AMBERHOME}/src/pmemd/src
FC = gfortran

install: extract

extract: $(objects)
	$(FC) -o $@ $^ $(pmemd_src)/AmberNetcdf.o  -I$(pmemd_src) -L${AMBERHOME}/lib -lnetcdf -lnetcdff

%.o : %.F90
	$(FC) -c $^ -I$(pmemd_src) -I${AMBERHOME}/include

clean:
	rm *.o *.mod extract
