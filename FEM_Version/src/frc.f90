program main

  call Initialize
  
  call Run

  call Save

  call Finalize

end program main

subroutine Run
  use globals
  implicit none
  real(rkind), external :: ppfun
  
  call EvolveBrutal(PsiCur,ppfun,AllPsis,ntsmax)

  ! call PetscSolve(PsiCur,1._rkind,PsiFinal,LambdaFinal)

end subroutine Run

subroutine ReadNamelist
  use prec_const
  use globals
  implicit none
  
  integer :: mp

  nz = 100
  nr = 100
  length = 1._rkind
  psiedge = -0.5_rkind
  
  ntsmax = 500
  tol = 1.e-8_rkind
  deltati = 1.e-2_rkind

  gaussorder = 4
  nboundarypoints = 5
  
  relax = 0.1_rkind

  Itotal = 2._rkind

  psimax = 0.05

  npsi = 50
  ntheta = 64

  namelist /frc/ &
       &    nr,nz,length,psiedge,ntsmax,psimax, &
       &    tol,deltati,gaussorder, nboundarypoints,relax,Itotal, &
       &    npsi,ntheta

  mp = 101
  open(mp, file ='nlfrc', delim = 'apostrophe', &
       FORM = 'formatted', action = 'read', status='old')
  read(mp,frc)
  close(mp)

  ! half length
  hlength = 0.5_rkind*length

end subroutine ReadNamelist

subroutine Initialize
#include <petsc/finclude/petsc.h>
  use globals
  use petsc
  implicit none
  integer :: i
  real(rkind) :: t0,t1
  real(rkind), external  :: ppfun
  real(rkind) :: Itot

  PetscErrorCode :: ierr

  call ReadNamelist

  allocate(Z(nz),R(nr),RR(nz,nr),ZZ(nz,nr),logRp1sR(nr), &
       &   PsiBC(nz,nr,4), &
       &   IndArray(nz*nr*4,3), &
       &   IndArrayInv(nz,nr,4), &
       &   RIntegralArray(0:6,nr-1), &
       &   ContributionMat(0:1,0:1,4,0:1,0:1,4,nr-1))

  deltar = 1._rkind/real(nr-1,rkind)
  do i=1,nr
     R(i) = real(i-1,rkind)/real(nr-1,rkind)
     RR(:,i) = real(i-1,rkind)/real(nr-1,rkind)
  end do

  deltaz = 1._rkind/real(nz-1,rkind)*length
  do i=1,nz
     Z(i) = real(i-1,rkind)/real(nz-1,rkind)*length - length/2._rkind
     ZZ(i,:) = real(i-1,rkind)/real(nz-1,rkind)*length - length/2._rkind
  end do

  dRdZ(1) = 1._rkind
  dRdZ(2) = deltaz
  dRdZ(3) = deltar
  dRdZ(4) = deltar*deltaz

  logRp1sR(1) = 0._rkind

  do i=2,nr
     logRp1sR(i) = log((1._rkind+R(i)/deltar)/(R(i)/deltar))
  end do

  ! (nR-2)*(nZ-2)*4 : bulk values
  ! 2*(nZ-2) : 2 unknowns per R=1 edge point, with a total of nZ points (psi and dpsi/dZ are known)
  ! 4*(nR-2) : 2 unknowns per Z=+-L/2 edge point (psi, dRpsi are known)
  ! 2 : 1 unknown (d2psi/dRdZ) per top corner
  ! No unknowns on bottom edge (including corners)
  nws = (nr-2)*(nz-2)*4 + 2*(nz-2)+4*(nr-2)+2
  nkws = nr*nz*4-nws
  ntot = nws+nkws

  allocate(B_BC(nws),PsiCur(nws),PsiFinal(nws),rhs(nws),AllPsis(nws,ntsmax+1))

  call HermiteMatrix

  call Ind_to_iRiZ
  call iRiZ_to_Ind

  call AllCoeffs

  call PsiBoundaryCondition

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF = MPI_COMM_SELF  

  ! call cpu_time(t0)
  ! call SETUPAB
  ! call cpu_time(t1)

  ! print*, (t1-t0)/24

  call ReadGuess

  ! call ReadPprime

  call InterpolatePsi(PsiGuess,nzguess,nrguess,PsiCur)

  call AllRIntegrals
  call Aij_all  
  call SETUPAB

  call TotalCurrent(PsiCur,ppfun,Itot)

  ilast = 1

end subroutine Initialize

function ppfun(psi)
  use prec_const
  use globals
  implicit none

  real(rkind) :: psi,ppfun
  real(rkind), external :: seval

  integer :: i
  real(rkind) :: a1,b1,d1

  ! a1 = 33.943
  ! b1 = 0.25867
  ! d1 = 0.920*41.244

  a1 = 170.64_rkind
  b1 = 0.69919_rkind
  d1 = 3.0156_rkind

  ! ppfun = 0._rkind
  ! do i=1,npprime
     ! ppfun = ppfun + Apprime(i)*psi**(i-1)
  ppfun = d1/cosh(a1*psi-b1)**2
  ! end do

  ! ppfun = seval(npprime,psi,Xpprime,Ypprime,Bpprime,Cpprime,Dpprime)

end function ppfun

subroutine PsiBoundaryCondition
  use globals
  implicit none

  integer :: iZ,iR

  psiBC = 0._rkind

  do iR=1,nR-1
     PsiBC(1,iR,1) = psiedge*R(iR)**2
     PsiBC(1,iR,3) = 2._rkind*psiedge*R(iR)
     PsiBC(nz,iR,1) = psiedge*R(iR)**2
     PsiBC(nz,iR,3) = 2._rkind*psiedge*R(iR)
  end do

  do iZ=2,nZ-1
     PsiBC(iZ,nR,1) = psiedge
  end do

  PsiBC(1,nR,1) = psiedge
  PsiBC(1,nR,3) = 2._rkind*psiedge
  PsiBC(nZ,nR,1) = psiedge
  PsiBC(nZ,nR,3) = 2._rkind*psiedge

end subroutine PsiBoundaryCondition

subroutine HermiteMatrix
  ! https://en.wikipedia.org/wiki/Bicubic_interpolation
  use globals
  implicit none
  integer, dimension(100) :: row,col
  real(rkind), dimension(100) :: dat
  integer :: i,irow,icol,j

  row = (/0,1,2,2,2,2,3,3,3,3,4,5,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9, &
       & 10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10, &
       & 11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11, &
       & 12,12,12,12, &
       & 13,13,13,13, &
       & 14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14, &
       & 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15/)

  col = (/0,4,0,1,4,5,0,1,4,5,8,12,8,9,12,13,8,9,12,13, &
       0,2,8,10,4,6,12,14, &
       0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15, &
       0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15, &
       0,2,8,10,4,6,12,14, &
       0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15, &
       0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)

  dat = (/1, &
       & 1 ,  &
       & -3,3,-2,-1, &
       & 2,-2,1,1, &
       & 1, &
       & 1, &
       & -3,3,-2,-1, &
       & 2,-2,1,1, &
       & -3,3,-2,-1, &
       & -3,3,-2,-1, &
       & 9,-9,-9,9,6,3,-6,-3,6,-6,3,-3,4,2,2,1, &
       & -6,6,6,-6,-3,-3,3,3,-4,4,-2,2,-2,-2,-1,-1, &
       & 2,-2,1,1, &
       & 2,-2,1,1, &
       & -6,6,6,-6,-4,-2,4,2,-3,3,-3,3,-2,-1,-2,-1, &
       & 4,-4,-4,4,2,2,-2,-2,2,-2,2,-2,1,1,1,1/)

  do i=1,100
     irow = row(i)+1
     icol = col(i)+1
     Hmat(irow,icol) = dat(i)
  end do

end subroutine HermiteMatrix

subroutine Finalize
#include <petsc/finclude/petsc.h>
  use globals
  use petsc
  implicit none

  PetscErrorCode :: ierr

  call MatDestroy(PAMat,ierr)

  call PetscFinalize(PETSC_NULL_CHARACTER,ierr)

  deallocate(R,Z,RR,ZZ,logRp1sR,PsiBC,IndArray,IndArrayInv,RIntegralArray,ContributionMat)
  deallocate(B_BC)
  deallocate(PsiGuess,PsiCur,rhs,AllPsis)
  ! deallocate(Xpprime,Ypprime,Bpprime,Cpprime,Dpprime)
  ! deallocate(Apprime)

end subroutine Finalize

subroutine InterpolatePsi(Psizrg,nzg,nrg,Psi)
  ! Transform a 2D psi array into an array adapted to the finite element
  ! Psizrg is typically a 2D guess z,r arry
  ! Psi is in the form of a size nws array
  use globals
  implicit none
  integer, intent(in) :: nzg,nrg
  real(rkind), dimension(nzg,nrg), intent(in) :: Psizrg
  real(rkind), dimension(nws), intent(out) :: Psi
  real(rkind), dimension(nzg) :: zg
  real(rkind), dimension(nrg) :: rg
  real(rkind) :: dzg,drg
  real(rkind) :: z0,z1,r0,r1
  integer :: ind
  integer :: iZ,iR,k
  integer :: izg,irg
  real(rkind) :: t,u
  real(rkind) :: zind,rind
  real(rkind) :: psi00,psi01,psi10,psi11

  z0 = -0.5_rkind*length
  z1 = 0.5_rkind*length

  r0 = 0._rkind
  r1 = 1._rkind

  call linspace(z0,z1,nzg,zg)
  call linspace(r0,r1,nrg,rg)

  dzg = 1._rkind/real(nzg-1,rkind)*length
  drg = 1._rkind/real(nrg-1,rkind)

  do ind=1,nws
     iZ = IndArray(ind,1)
     iR = IndArray(ind,2)
     k = IndArray(ind,3)

     zind = (Z(iZ)+0.5_rkind*length)/dzg
     rind = R(iR)/drg

     izg = min(int(zind),nzg-2)
     irg = min(int(rind),nrg-2)
     t = zind - real(izg,rkind)
     u = rind - real(irg,rkind)

     izg = izg+1
     irg = irg+1

     psi00 = psizrg(izg  ,irg  )
     psi10 = psizrg(izg+1,irg  )
     psi01 = psizrg(izg  ,irg+1)
     psi11 = psizrg(izg+1,irg+1)

     if (k.eq.1) then
        Psi(ind) = (psi00*(1._rkind-t)*(1._rkind-u) + &
             &      psi10*t*(1._rkind-u) +  &
             &      psi01*(1._rkind-t)*u + &
             &      psi11*t*u)
     else if (k.eq.2) then
        Psi(ind) = (psi10 - psi00)*(1._rkind-u)/dzg + (psi11 - psi01)*u/dzg
     else if (k.eq.3) then
        Psi(ind) = (psi01 - psi00)*(1._rkind-t)/drg + (psi11 - psi10)*t/drg
     else if (k.eq.4) then
        Psi(ind) = (psi11 + psi00 - psi10 - psi01)/(drg*dzg)
     end if
  end do

end subroutine InterpolatePsi

subroutine Evolve(PsiIni,fun,Psis,Cs,dts,nts)
  use globals
  implicit none
  real(rkind), dimension(nws), intent(in) :: PsiIni
  integer, intent(in) :: nts
  real(rkind), dimension(nws,nts+1), intent(out) :: Psis
  real(rkind), dimension(nts+1), intent(out) :: Cs
  real(rkind), dimension(nts+1), intent(out) :: dts
  real(rkind), external :: fun
  
  real(rkind) :: C
  real(rkind), dimension(nws) :: Psi1,Psi2,Psi2_star
  real(rkind) :: ratio,n_,n_star
  real(rkind) :: Itot
  real(rkind) :: dtpsi
  real(rkind) :: deltat
  real(rkind) :: Icur
  logical :: advance
  integer :: i,k

  Psis = 0._rkind
  Cs = 0._rkind
  dts = 0._rkind

  i = 1
  dtpsi = tol
  Psis(:,i) = PsiIni
  Cs(i) = 1._rkind
  C = 1._rkind
  call TotalCurrent(PsiIni,fun,Itot)
  deltat = deltati

  advance = .true.
  k = 0

  do while (i.le.nts .and. dtpsi.ge.tol)
     if (advance) then
        write(*,'(A)') '=================='
        write(*,'(A,I3)') 'Step ',i
        write(*,'(A)') '=================='
     end if
     
     call FEMStep(deltat*2._rkind,Psis(:,i),funC,Psi2_star,.true.)
     call FEMStep(deltat,Psis(:,i),funC,Psi1,.false.)
     call FEMStep(deltat,Psi1,funC,Psi2,.true.)

     call DiffNorm(nws,Psis(:,i),Psi2,n_)
     call DiffNorm(nws,Psi2,Psi2_star,n_star)
     ratio = n_star/n_
     dtpsi = n_/deltat

     write(*,'(I3,4E11.3)') k,deltat,Cs(i),ratio,dtpsi/tol
     
     if (ratio.le.0.05_rkind) then
        dts(i) = deltat
        i = i+1
        k=0
        Psis(:,i) = Psi2(:)
        Cs(i) = C
        deltat = deltat*4._rkind
        call TotalCurrent(Psi2,fun,Icur)
        C = Itot/Icur
        advance = .true.
     else if (ratio.gt.0.05_rkind .and. ratio.le.0.1_rkind) then
        dts(i) = deltat
        i = i+1
        k=0
        Psis(:,i) = Psi2(:)
        Cs(i) = C
        deltat = deltat*2._rkind
        call TotalCurrent(Psi2,fun,Icur)
        C = Itot/Icur
        advance = .true.
     else if (ratio.gt.0.1_rkind .and. ratio.le.0.3_rkind) then
        dts(i) = deltat
        i = i+1
        k=0
        Psis(:,i) = Psi2(:)
        Cs(i) = C
        call TotalCurrent(Psi2,fun,Icur)
        C = Itot/Icur
        advance = .true.
     else if (ratio.gt.0.3_rkind .and. ratio.le.0.4_rkind) then
        dts(i) = deltat
        i = i+1
        k=0
        Psis(:,i) = Psi2(:)
        Cs(i) = C
        deltat = deltat*0.5_rkind
        call TotalCurrent(Psi2,fun,Icur)
        C = Itot/Icur
        advance = .true.
     else if (ratio.gt.0.4_rkind .and. ratio.le.0.6_rkind) then
        dts(i) = deltat
        i = i+1
        k=0
        Psis(:,i) = Psi2(:)
        Cs(i) = C
        deltat = deltat*0.25_rkind
        call TotalCurrent(Psi2,fun,Icur)
        C = Itot/Icur
        advance = .true.
     else
        k=k+1
        deltat = deltat*0.0625_rkind
        advance = .false.
     end if
  end do

contains
  
  function funC(psi)
    implicit none
    real(rkind) :: psi,funC

    funC = fun(psi)*C

  end function funC

end subroutine Evolve

subroutine FEMStep(deltat,Psi,fun,PsiSol,compute_rhs)
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscksp.h>
  use globals
  use petscksp
  implicit none

  real(rkind), intent(in) :: deltat
  logical, intent(in) :: compute_rhs
  real(rkind), dimension(nws), intent(in) :: Psi
  real(rkind), dimension(nws), intent(out) :: PsiSol
  real(rkind), external :: fun
  integer :: ind
  real(rkind) :: norm2,norm2_

  Vec :: u,x,b
  PetscInt :: pnws
  PetscErrorCode :: ierr
  Mat :: PM
  PetscReal :: one,negdt
  PetscReal :: ptol
  PetscInt :: ione
  PetscInt :: II
  PetscInt, dimension(0:nws-1) :: ix
  PetscReal, dimension(0:nws-1) :: pvals

  PC :: pc
  KSP :: ksp

  if (compute_rhs) call RightHandSide(Psi,fun,rhs)

  print*, 'Matrix solve'

  pnws = nws
  negdt = -deltat
  do II=0,nws-1
     ix(II) = II
     pvals(II) = psi(II+1) + negdt*(B_BC(II+1)+rhs(II+1))
  end do
  ione = 1
  one = 1.d0

  call MatCreate(PETSC_COMM_WORLD,PM,ierr)
  call MatSetSizes(PM,PETSC_DECIDE,PETSC_DECIDE,pnws,pnws,ierr)
  call MatSetType(PM,MATSEQAIJ,ierr)
  call MatSeqAIJSetPreallocation(PM,ione,PETSC_NULL_INTEGER,ierr)

  do II=0,pnws-1
     call MatSetValues(PM,ione,II,ione,II,one,INSERT_VALUES,ierr)
  end do
  call MatAssemblyBegin(PM,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(PM,MAT_FINAL_ASSEMBLY,ierr)

  call MatAXPY(PM,negdt,PAMat,DIFFERENT_NONZERO_PATTERN,ierr)

  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
  call KSPSetOperators(ksp,PM,PM,ierr)

  call KSPGetPC(ksp,pc,ierr)
  ptol = 1.e-7
  call KSPSetTolerances(ksp,ptol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)

  call KSPSetFromOptions(ksp,ierr)
  call PCSetFromOptions(pc,ierr)

  call VecCreateSeq(PETSC_COMM_WORLD,pnws,b,ierr)
  call VecDuplicate(b,x,ierr)

  ! Set b to psi - deltat*(B + rhs)
  call VecSetValues(b,pnws,ix,pvals,INSERT_VALUES,ierr)
  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)

  call KSPSolve(ksp,b,x,ierr)

  call VecGetValues(x,pnws,ix,pvals,ierr)

  PsiSol(1:nws) = pvals(0:pnws-1)

  call MatDestroy(PM,ierr)
  call VecDestroy(x,ierr)
  call VecDestroy(b,ierr)
  
end subroutine FEMStep

subroutine EvolveBrutal(PsiIni,fun,Psis,nts)
  use globals
  implicit none
  real(rkind), dimension(nws), intent(in) :: PsiIni
  integer, intent(in) :: nts
  real(rkind), dimension(nws,nts+1), intent(out) :: Psis
  real(rkind), dimension(nts+1) :: Cs
  real(rkind), external :: fun
  
  real(rkind) :: C
  real(rkind), dimension(nws) :: Psi1,Psi2
  real(rkind) :: ratio,n_,n_star
  real(rkind) :: Itot
  real(rkind) :: dtpsi
  real(rkind) :: deltat
  real(rkind) :: Icur
  real(rkind) :: rmax,psimaxcur
  logical :: advance
  integer :: i,k
  logical :: equatorial

  ! if nz is odd, magnetic axis ought to be located in the equatorial plane z=0
  equatorial = .false.
  if (mod(nz,2).eq.1) equatorial = .true.

  Psis = 0._rkind
  Cs = 0._rkind

  i = 1
  dtpsi = tol
  Psis(:,i) = PsiIni
  Cs(i) = 1._rkind
  C = 1._rkind
  C = 29.6195935435900_rkind
  ! call TotalCurrent(PsiIni,fun,Itot)
  ! Itot = Itotal
  deltat = deltati

  Psi1(:) = PsiIni(:)

  do while (i.le.nts .and. dtpsi.ge.tol)
     i = i+1
     call FEMStepBrutal(Psi1,funC,Psi2)
     call DiffNorm(nws,Psi1,Psi2,dtpsi)
     call TotalCurrent(Psi2,fun,Icur)
     call PsiMaximum(Psi2,rmax,psimaxcur,equatorial)
     ! C = Itot/Icur
     C = C*psimax/psimaxcur
     write(*,'(2E12.4)'), C,dtpsi
     Psis(:,i) = Psi2(:)
     Psi1(:) = Psi2(:)
  end do

  ilast = i
  LambdaFinal = C

  call TotalCurrent(Psi2,funC,Icur)

  write(*,'(A,2E12.4)') 'Total current at convergence = ', Icur

contains
  
  function funC(psi)
    implicit none
    real(rkind) :: psi,funC

    funC = fun(psi)*C

  end function funC
  
end subroutine EvolveBrutal

subroutine FEMStepBrutal(Psi,fun,PsiSol)
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscksp.h>
  use globals
  use petscksp
  implicit none

  real(rkind), dimension(nws), intent(in) :: Psi
  real(rkind), dimension(nws), intent(out) :: PsiSol
  real(rkind), external :: fun
  integer :: ind
  real(rkind) :: norm2,norm2_

  Vec :: u,x,b
  PetscInt :: pnws
  PetscErrorCode :: ierr
  Mat :: PM
  PetscReal :: one,negdt
  PetscReal :: ptol
  PetscInt :: ione
  PetscInt :: II
  PetscInt, dimension(0:nws-1) :: ix
  PetscReal, dimension(0:nws-1) :: pvals

  PC :: pc
  KSP :: ksp

  call RightHandSide(Psi,fun,rhs)

  print*, 'Matrix solve'

  pnws = nws
  do II=0,nws-1
     ix(II) = II
     pvals(II) = (B_BC(II+1)+rhs(II+1))
  end do
  ione = 1
  one = 1.d0

  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
  call KSPSetOperators(ksp,PAMat,PAMat,ierr)

  call KSPGetPC(ksp,pc,ierr)
  ptol = 1.e-7
  call KSPSetTolerances(ksp,ptol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)

  call KSPSetFromOptions(ksp,ierr)
  call PCSetFromOptions(pc,ierr)

  call VecCreateSeq(PETSC_COMM_WORLD,pnws,b,ierr)
  call VecDuplicate(b,x,ierr)

  ! Set b to (B + rhs)
  call VecSetValues(b,pnws,ix,pvals,INSERT_VALUES,ierr)
  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)

  call KSPSolve(ksp,b,x,ierr)

  call VecGetValues(x,pnws,ix,pvals,ierr)

  PsiSol(1:nws) = Psi(1:nws)*(1._rkind - relax) + relax*pvals(0:pnws-1)

  call VecDestroy(x,ierr)
  call VecDestroy(b,ierr)
  call KSPDestroy(ksp,ierr)
  
end subroutine FEMStepBrutal

subroutine PetscDiffNorm(x,y,norm,ierr)
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscvec.h>
  use petsc
  implicit none

  Vec  :: x,y
  PetscReal :: norm
  PetscReal :: norm_
  PetscErrorCode :: ierr

  PetscScalar :: negone

  negone = -1.d0

  call VecAXPY(y,negone,x,ierr)
  call VecNorm(y,NORM_2,norm,ierr)
  call VecNorm(x,NORM_2,norm_,ierr)

  norm = norm/norm_

end subroutine PetscDiffNorm

subroutine DiffNorm(n,x,y,norm)
  use prec_const
  implicit none
  integer, intent(in) :: n
  real(rkind), dimension(n), intent(in) :: x,y
  real(rkind), intent(out) :: norm

  norm = norm2(y-x)/norm2(x)

end subroutine DiffNorm
