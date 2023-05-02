FC=mpiifort
FFLAGS=-O3 -fpp -xHost -qopenmp -parallel -mkl=parallel -mcmodel=medium -shared-intel
SRCMOD=precconst.f90 modules.f90 
SRC=frc.f90 zbrent.f90
OBJMOD=${SRCMOD:.f90=.o}
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

frc: $(OBJMOD) $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

clean:
	@rm -f *.mod *.o frc

