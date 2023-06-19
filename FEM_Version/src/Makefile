FC=mpiifort
FFLAGS=-O3 -fpp -xHost -qopenmp -parallel -mkl=parallel -mcmodel=medium -shared-intel
SRCMOD=precconst.f90 modules.f90
#SRC=frc.f90 indexing.f90 matrix.f90 quadpack_double.f90 integration.f90 gauss.f90 heavisideintegration.f90 rhs.f90 inputoutput.f90 mesh.f90 solver.f90 tseval.f90
SRC=compression.f90 indexing.f90 matrix.f90 quadpack_double.f90 integration.f90 gauss.f90 heavisideintegration.f90 rhs.f90 inputoutput.f90 mesh.f90 polyfit.f90 solver.f90 tseval.f90 
OBJMOD=${SRCMOD:.f90=.o}
OBJ=${SRC:.f90=.o}
VECLIB  = -L$(MKLROOT)/lib/intel64 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -L$(PETSC_DIR)/lib -lpetsc
PETSC = -I$(PETSC_DIR)/include -I${PETSC_DIR}/${PETSC_ARCH}/include

%.o: %.f90
	$(FC) $(FFLAGS) $(VECLIB) $(PETSC) -o $@ -c $<

frc: $(OBJMOD) $(OBJ)
	$(FC) $(FFLAGS) $(VECLIB) $(PETSC) -o $@ $^

clean:
	@rm -f *.mod *.o frc

