subroutine PetscSolve(psistart,psisol,Lambdasol)
#include <petsc/finclude/petscsnes.h>
  use petscsnes
  use globals
  implicit none
  real(rkind), dimension(nws), intent(in) :: psistart
  real(rkind), dimension(nws), intent(out) :: psisol
  real(rkind), intent(out) :: Lambdasol

  SNES :: snes
  SNESLineSearch :: linesearch
  KSP :: ksp
  PC :: pc
  Mat :: Jmat
  Vec :: x,rvec
  PetscInt :: size,pnws
  PetscInt, dimension(0:nws-1) :: ix
  PetscScalar, dimension(0:nws-1) :: pvals
  PetscInt :: II
  PetscScalar :: lambdaval
  PetscScalar :: one,zero
  PetscInt :: ione
  PetscErrorCode :: ierr
  PetscScalar :: normrvec

  PetscScalar, dimension(0:0) :: Larray
  PetscInt, dimension(0:0) :: Lind

  external :: FormFunction,FormJacobian

  pnws = nws
  size = nws+1
  one = 1.d0
  ione = 1
  zero = 0.d0

  ! Required because PETSc is really picky with types sometimes
  ! Argument of VecSetValues and VecGetValues should be arrays even if only one
  ! value is inserted
  Lind(0) = pnws
  Larray(0) = LambdaIni
  
  do II=0,pnws-1
     ix(II) = II
     pvals(II) = psistart(II+1)
  end do

  ! Create PETSc vectors
  call VecCreateSeq(PETSC_COMM_WORLD,size,x,ierr)
  call VecDuplicate(x,rvec,ierr)

  ! Put guess in x
  call VecSetValues(x,pnws,ix,pvals,INSERT_VALUES,ierr)
  call VecSetValues(x,ione,Lind,Larray,INSERT_VALUES,ierr)
  ! Assemble vector
  call VecAssemblyBegin(x,ierr)
  call VecAssemblyEnd(x,ierr)  

  ! Create Jacobian matrix context
  call MatCreate(PETSC_COMM_WORLD,JMat,ierr)
  call MatSetSizes(JMat,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  call MatSetType(JMat,MATSEQAIJ,ierr)

  ! Create SNES context
  call SNESCreate(PETSC_COMM_WORLD,snes,ierr)
  call SNESSetFunction(snes,rvec,FormFunction,PETSC_NULL_INTEGER,ierr)
  ! call SNESSetJacobian(snes,Jmat,Jmat,FormJacobian,PETSC_NULL_INTEGER,ierr)
  call SNESSetFromOptions(snes,ierr)
  call SNESGetLineSearch(snes, linesearch, ierr)
  call SNESLineSearchSetType(linesearch,'basic', ierr)
  ! call SNESKSPSetUseEW(snes,PETSC_TRUE,ierr)

  call SNESGetKSP(snes,ksp,ierr)
  call KSPSetFromOptions(ksp,ierr)
  call KSPGetPC(ksp,pc,ierr)
  call PCSetFromOptions(pc,ierr)

  ! call FormFunction(snes,x,rvec,PETSC_NULL_INTEGER,ierr)
  ! call VecNorm(rvec,NORM_2,normrvec,ierr)

  ! Solve the nonlinear problem with PETSc
  call SNESSolve(snes,PETSC_NULL_VEC,x,ierr)

  ! Put psi solution in pvals
  call VecGetValues(x,pnws,ix,pvals,ierr)
  ! Get solution for the Lagrange multiplier
  call VecGetValues(x,ione,Lind,Larray,ierr)  
  lambdaval = Larray(0)

  do II=0,pnws-1
     psisol(II+1) = pvals(II)
  end do
  Lambdasol = lambdaval

  ! Destroy the PETSc vectors, matrix and SNES context
  call VecDestroy(x,ierr)
  call VecDestroy(rvec,ierr)
  call MatDestroy(Jmat,ierr)
  call SNESDestroy(snes,ierr)

end subroutine PetscSolve

subroutine Preconditioning(psi,lambdaval,rhs,fvals)
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscvec.h>
  use prec_const
  use globals, only : nws, B_BC, PAMat
  use petscksp
  implicit none

  real(rkind), dimension(nws) :: psi,rhs
  real(rkind) :: lambdaval
  
  Vec :: b,x
  KSP :: ksp
  PC :: pc
  PetscInt :: pnws
  PetscErrorCode :: ierr
  PetscReal :: one
  PetscReal :: ptol
  PetscInt :: ione
  PetscInt :: II
  PetscInt, dimension(0:nws-1) :: ix
  PetscScalar, dimension(0:nws-1) :: fvals

  pnws = nws
  ione = 1
  one = 1.d0
  do II=0,pnws-1
     ix(II) = II
     fvals(II) = (B_BC(II+1)+lambdaval*rhs(II+1))
  end do
  
  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
  call KSPSetOperators(ksp,PAMat,PAMat,ierr)

  call KSPSetOptionsPrefix(ksp,'prec_',ierr)
  call KSPGetPC(ksp,pc,ierr)
  call PCSetOptionsPrefix(pc,'prec_',ierr)
  ptol = 1.e-7
  call KSPSetTolerances(ksp,ptol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
  call KSPSetFromOptions(ksp,ierr)
  call PCSetFromOptions(pc,ierr)

  call VecCreateSeq(PETSC_COMM_WORLD,pnws,b,ierr)
  call VecDuplicate(b,x,ierr)

  ! Set b to (B + L*rhs)
  call VecSetValues(b,pnws,ix,fvals,INSERT_VALUES,ierr)
  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)

  call KSPSolve(ksp,b,x,ierr)

  call VecGetValues(x,pnws,ix,fvals,ierr)

  do II=0,pnws
     fvals(II) = psi(II+1)-fvals(II)
  end do

  call KSPDestroy(ksp,ierr)
  call VecDestroy(x,ierr)
  call VecDestroy(b,ierr)

end subroutine Preconditioning

subroutine FormFunction(snes,x,f,user,ierr)
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscvec.h>
  use petscsnes
  use prec_const
  use globals, only : nws, PAMat, psimax, B_BC
  implicit none

  SNES :: snes
  Vec :: x,f
  Vec :: x_nws,b_nws
  PetscInt :: user
  PetscErrorCode :: ierr

  PetscInt :: pnws
  PetscInt :: size
  PetscInt, dimension(0:nws-1) :: ix
  PetscScalar, dimension(0:nws-1) :: psi,fvals
  PetscInt :: ione
  PetscInt :: II
  PetscScalar :: lambdaval

  PetscScalar, dimension(0:0) :: Larray
  PetscInt, dimension(0:0) :: Lind

  real(rkind) :: rmax,psimaxval
  real(rkind), dimension(0:nws-1) :: rhs
  real(rkind), external :: ppfun

  pnws = nws
  size = pnws+1
  ione = 1
  
  do II=0,nws-1
     ix(II) = II
  end do

  call VecGetValues(x,pnws,ix,psi,ierr)
  call PsiMaximum(psi,rmax,psimaxval,.true.)
  call RightHandSide(psi,ppfun,rhs)

  Lind(0) = pnws
  call VecGetValues(x,ione,Lind,Larray,ierr)
  lambdaval = Larray(0)

  ! call VecCreateSeq(PETSC_COMM_WORLD,pnws,x_nws,ierr)
  ! call VecDuplicate(x_nws,b_nws,ierr)

  ! call VecSetValues(x_nws,pnws,ix,psi,INSERT_VALUES,ierr)
  ! call VecAssemblyBegin(x_nws,ierr)
  ! call VecAssemblyEnd(x_nws,ierr)
  
  ! call MatMult(PAMat,x_nws,b_nws,ierr)

  call Preconditioning(psi,lambdaval,rhs,fvals)

  ! call VecGetValues(b_nws,pnws,ix,fvals,ierr)

  ! call VecDestroy(x_nws,ierr)
  ! call VecDestroy(b_nws,ierr)

  ! fvals(:) = fvals(:) - lambdaval*(B_BC(:)+rhs(:))

  call VecSetValues(f,pnws,ix,fvals,INSERT_VALUES,ierr)
  Larray(0) = psimaxval - psimax
  call VecSetValues(f,ione,Lind,Larray,INSERT_VALUES,ierr)

  call VecAssemblyBegin(f,ierr)
  call VecAssemblyEnd(f,ierr)
  
end subroutine FormFunction

subroutine FormJacobian(snes,x,Jmat,Pmat,flag,user,ierr)
#include <petsc/finclude/petscsnes.h>
  use prec_const
  use globals, only : nws, IndArray, nr,nz,IndArrayInv,PsiBC,B_BC
  use petscsnes
  implicit none

  SNES :: snes
  Vec :: x
  Mat :: Jmat,Pmat
  PetscInt :: flag,user
  PetscErrorCode :: ierr
  
  PetscInt :: pnws,size
  PetscInt :: dummynz
  PetscInt :: ione
  PetscInt :: II,JJ
  PetscInt, dimension(0:nws) :: nnz
  PetscScalar :: pval,one

  integer :: ind,jnd
  real(rkind) :: val,valBC
  integer :: iZo,iRo,ko
  integer :: iZp,iRp,kp
  integer :: i,j

  real(rkind) :: res
  real(rkind) :: lambdaval

  PetscInt, dimension(0:nws-1) :: ix
  PetscScalar, dimension(0:nws-1) :: psi

  PetscScalar, dimension(0:0) :: Larray
  PetscInt, dimension(0:0) :: Lind

  real(rkind) :: rmax,psimaxval
  real(rkind), dimension(nws) :: rhs
  real(rkind), external :: ppfun

  integer :: ir
  integer, dimension(4) :: inds
  real(rkind), dimension(4) :: dpsimaxdpsivals

  ! ione = 1
  ! pnws = nws
  ! size = pnws+1
  
  ! do II=0,pnws-1
  !    ix(II) = II
  ! end do

  ! call VecGetValues(x,pnws,ix,psi,ierr)
  ! call PsiMaximum(psi,rmax,psimaxval,.true.)
  ! call RightHandSide(psi,ppfun,rhs)
  ! ! Lind(0) = pnws
  ! ! call VecGetValues(x,ione,Lind,Larray,ierr)
  ! ! lambdaval = Larray(0)

  ! dummynz = 1
  ! nnz = 0
  ! do ind=1,nws
  !    iZo = IndArray(ind,1)
  !    iRo = IndArray(ind,2)
  !    ko = IndArray(ind,3)
  !    do i=-1,1
  !       iZp = iZo+i
  !       do j=-1,1
  !          iRp = iRo+j
  !          if (iRp.ge.1 .and. iRp.le.nr .and. iZp.ge.1 .and. iZp.le.nz) then
  !             do kp=1,4
  !                jnd = IndArrayInv(iZp,iRp,kp)
  !                if (jnd.le.nws) then
  !                   nnz(ind-1) = nnz(ind-1)+1
  !                end if
  !             end do
  !          end if
  !       end do
  !    end do
  !    ! For the last column of the matrix
  !    nnz(ind-1) = nnz(ind-1)+1
  ! end do
  ! nnz(nws) = 4

  ! ! call MatCreate(PETSC_COMM_WORLD,Pmat,ierr)
  ! ! size = nws+1
  ! ! call MatSetSizes(Pmat,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  ! ! call MatSetType(Pmat,MATSEQAIJ,ierr)
  ! call MatSeqAIJSetPreallocation(Pmat,dummynz,nnz,ierr)

  ! one = 1.d0

  ! do ind=1,nws
  !    II = ind-1
  !    valBC = 0._rkind

  !    iZo = IndArray(ind,1)
  !    iRo = IndArray(ind,2)
  !    ko = IndArray(ind,3)

  !    do i=-1,1
  !       iZp = iZo+i
  !       do j=-1,1
  !          iRp = iRo+j
  !          if (iRp.ge.1 .and. iRp.le.nr .and. iZp.ge.1 .and. iZp.le.nz) then
  !             do kp=1,4
  !                jnd = IndArrayInv(iZp,iRp,kp)                 
  !                call Aij_bulk(ind,jnd,val)
  !                if (jnd.le.nws) then
  !                   JJ = jnd-1
  !                   pval = val
  !                   call MatSetValues(Pmat,ione,II,ione,JJ,pval,INSERT_VALUES,ierr)
  !                else
  !                   valBC = valBC - val*PsiBC(iZp,iRp,kp)
  !                end if
  !             end do
  !          end if
  !       end do
  !    end do
  !    B_BC(ind) = valBC
  ! end do

  ! do ind=1,nws
  !    II = ind-1
  !    pval = -(B_BC(ind) + rhs(ind))
  !    call MatSetValues(Pmat,ione,II,ione,pnws,pval,INSERT_VALUES,ierr)
  ! end do

  ! ! Set how psimaxval varies with psi
  ! call dPsimaxdPsi(psi,rmax,inds,dpsimaxdpsivals)

  ! do II=1,4
  !    JJ = inds(II)-1
  !    pval = dpsimaxdpsivals(II)
  !    call MatSetValues(Pmat,ione,pnws,ione,JJ,pval,INSERT_VALUES,ierr)
  ! end do
  
  call MatAssemblyBegin(Pmat,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Pmat,MAT_FINAL_ASSEMBLY,ierr)
  
  if (Jmat .ne. Pmat) then
     call MatAssemblyBegin(Jmat,MAT_FINAL_ASSEMBLY,ierr)
     call MatAssemblyEnd(Jmat,MAT_FINAL_ASSEMBLY,ierr)
  end if

end subroutine FormJacobian

subroutine dPsimaxdPsi(psi,rmax,inds,dpsimaxdpsivals)
  use globals
  implicit none
  real(rkind), dimension(nws), intent(in) :: psi
  real(rkind), intent(in) :: rmax
  integer, intent(out), dimension(4) :: inds
  real(rkind), dimension(4), intent(out) :: dpsimaxdpsivals

  real(rkind) :: psi0,psi1,psip0,psip1
  real(rkind) :: dpsidpsi0,dpsidpsi1,dpsidpsip0,dpsidpsip1
  integer :: ir,iz
  real(rkind) :: a,b,c,delta
  real(rkind) :: tmax,tmax_plus,tmax_minus

    iz = (nz+1)/2
    ir = min(int(rmax/deltar),nr-2)+1

    inds(1) = IndArrayInv(iz,ir,1)
    inds(2) = IndArrayInv(iz,ir+1,1)
    inds(3) = IndArrayInv(iz,ir,3)
    inds(4) = IndArrayInv(iz,ir+1,3)

    psi0 = psi(inds(1))
    psi1 = psi(inds(2))
    psip0 = psi(inds(3))
    psip1 = psi(inds(4))

    a = 3._rkind*psip0 + 3._rkind*psip1 + 6._rkind*psi0 - 6._rkind*psi1
    b = -4._rkind*psip0 - 2._rkind*psip1 - 6._rkind*psi0 + 6._rkind*psi1
    c = psip0

    delta = b**2-4._rkind*a*c

    tmax_plus = (-b+sqrt(delta))/(2._rkind*a)
    tmax_minus = (-b-sqrt(delta))/(2._rkind*a)
    
    if (tmax_plus.ge.0._rkind .and. tmax_plus.le.1._rkind) then
       tmax = tmax_plus
       psimax = p(tmax)
    else
       tmax = tmax_minus
       psimax = p(tmax)
    end if

    dpsidpsi0 = 2._rkind*tmax**3-3._rkind*tmax**2+1._rkind
    dpsidpsi1 = -2._rkind*tmax**3+3._rkind*tmax**2
    dpsidpsip0 = tmax**3-2._rkind*tmax**2+tmax
    dpsidpsip1 = tmax**3-tmax**2

    dpsimaxdpsivals(1) = dpsidpsi0
    dpsimaxdpsivals(2) = dpsidpsi1
    dpsimaxdpsivals(3) = dpsidpsip0
    dpsimaxdpsivals(4) = dpsidpsip1

contains
  
  function p(t)
    real(rkind) :: t,p

    p =    psi0*(2*t**3 - 3*t**2 + 1) + &
         & psip0*(t**3 - 2*t**2 + t) +  &
         & psi1*(-2*t**3 + 3*t**2) +    &
         & psip1*(t**3 - t**2)
  end function p
  
end subroutine dPsimaxdPsi
