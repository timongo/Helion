program main
  implicit none
  PetscErrorCode :: ierr
  PetscInt :: size

  size = 10

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF = MPI_COMM_SELF  

  call PetscSolve(size)

  call PetscFinalize(PETSC_NULL_CHARACTER,ierr)

end program main

subroutine PetscSolve(size)
#include <petsc/finclude/petscsnes.h>
  use petscsnes
  implicit none

  PetscInt :: size

  SNES :: snes
  SNESLineSearch :: linesearch
  KSP :: ksp
  PC :: pc
  Mat :: Jmat
  Vec :: x,rvec
  PetscInt, dimension(0:size-1) :: ix
  PetscScalar, dimension(0:size-1) :: pvals
  PetscInt :: II
  PetscScalar :: one,zero
  PetscInt :: ione
  PetscErrorCode :: ierr

  external :: FormFunction,FormJacobian

  one = 1.d0
  zero = 0.d0
  ione = 1
  
  call VecCreateSeq(PETSC_COMM_WORLD,size,x,ierr)
  call VecDuplicate(x,rvec,ierr)

  call VecSet(one,x,ierr)

  call MatCreate(PETSC_COMM_WORLD,JMat,ierr)
  call MatSetSizes(JMat,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  call MatSetType(JMat,MATSEQAIJ,ierr)

  call SNESCreate(PETSC_COMM_WORLD,snes,ierr)
  call SNESSetFunction(snes,rvec,FormFunction,PETSC_NULL_INTEGER,ierr)
  call SNESSetJacobian(snes,Jmat,Jmat,FormJacobian,PETSC_NULL_INTEGER,ierr)
  ! call SNESSetFromOptions(snes,ierr)
  ! call SNESGetLineSearch(snes, linesearch, ierr)
  ! call SNESLineSearchSetType(linesearch, 'basic', ierr)
  ! call SNESKSPSetUseEW(snes,PETSC_TRUE,ierr)

  call SNESGetKSP(snes,ksp,ierr)
  call KSPSetFromOptions(ksp,ierr)
  call KSPGetPC(ksp,pc,ierr)
  call PCSetFromOptions(pc,ierr)

  ! call FormFunction(snes,x,rvec,PETSC_NULL_INTEGER,ierr)
  
  ! call SNESSolve(snes,PETSC_NULL_VEC,x,ierr)

  call VecGetValues(x,pnws,ix,pvals,ierr)

  call VecGetValues(x,ione,Cind,Carray,ierr)  
  Cval = Carray(0)

  do II=0,pnws-1
     psisol(II+1) = pvals(II)
  end do
  Csol = Cval

  call VecDestroy(x,ierr)
  call VecDestroy(rvec,ierr)
  call MatDestroy(Jmat,ierr)
  call SNESDestroy(snes,ierr)

  Csol = 0._rkind
  psisol = 0._rkind

end subroutine PetscSolve
