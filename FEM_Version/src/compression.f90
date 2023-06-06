program main

  call Initialize

  call Compress

  call Save

  call Finalize

end program main

subroutine Compress
  use prec_const
  use sizes_indexing
  use globals
  implicit none
  integer, parameter :: max_iter = 40
  real(rkind), parameter :: tol_comp = 1.0e-4
  real(rkind), dimension(npsi) :: S, psi, P, Vprime, P_oldV_old, PV
  real(rkind), dimension(npsi,ntheta+1) :: ZMesh,RMesh,JacobMesh
  real(rkind) :: delta_psiedge
  integer :: iter

  ! Calculate initial psi, pressure, jacobian, Vpime  
  call Run
  call Mesh(inds_r,LambdaFinal,inds_r%PsiFinal,npsi,ntheta,ZMesh,RMesh,JacobMesh,S,psi,P)
  call calculate_Vprime(JacobMesh, Vprime)

  ! Compute initial PV^(5/3)
  P_oldV_old = P * Vprime**(5.0/3.0)
  
  ! delta_psiedge
  delta_psiedge = -0.2
  ! Increase psiedge
  psiedge = psiedge + delta_psiedge
  
  ! Calculate new psi, pressure, jacobian, Vprime
  call Run
  call Mesh(inds_r,LambdaFinal,inds_r%PsiFinal,npsi,ntheta,ZMesh,RMesh,JacobMesh,S,psi,P)
  call calculate_Vprime(JacobMesh, Vprime)

  ! Compute new PV^(5/3)
  PV = P * Vprime**(5.0/3.0)

  ! Iterate until the relative error is below the tolerance
  iter = 0
  do while (norm2( (PV - P_oldV_old) / P_oldV_old) > tol_comp)
     ! Compute the new pressure based on P_oldV_old^(5/3)
     P = P_oldV_old / Vprime**(5.0/3.0)
     ! Calculate new AP_NL and fit it
     call calculate_AP_NL(psi, P, AP_NL)

     ! Recompute based on the updated AP_NL
     call Run 
     call Mesh(inds_r,LambdaFinal,inds_r%PsiFinal,npsi,ntheta,ZMesh,RMesh,JacobMesh,S,psi,P)
     call calculate_Vprime(JacobMesh, Vprime)

     ! Update PV for the recomputed psi
     PV = P * Vprime**(5.0/3.0)

     ! Check the number of iterations
     if (iter >= max_iter) then
        exit
     end if

     iter = iter + 1
  end do

  if (iter >= max_iter) then
     write(*, *) "Maximum iterations reached without convergence for psiedge =", psiedge
  else
     write(*, *) "Converged after ", iter, " iterations for psiedge =", psiedge
  end if
  
end subroutine Compress

subroutine Run
  use globals
  implicit none
  real(rkind), external :: ppfun  

  if (usepetsc) then
     ! coarse grid solve
     call PetscSolve(inds_c,LambdaFinal)
     
     ! Interpolation of coarse result on fine grid
     inds_c%PsiCur = inds_c%PsiFinal
     call InterpolatePsi2(inds_c,inds_r)
     
     ! Refined grid solve
     LambdaIni = LambdaFinal
     call PetscSolve(inds_r,LambdaFinal)
  else
     call Evolve(inds_c,inds_c%PsiCur,ppfun,inds_c%AllPsis,ntsmax)
     inds_c%PsiFinal = inds_c%AllPsis(:,ilast)

     ! Interpolation of coarse result on fine grid
     inds_c%PsiCur = inds_c%PsiFinal
     call InterpolatePsi2(inds_c,inds_r)
     
     ! Refined grid solve
     LambdaIni = LambdaFinal
     call PetscSolve(inds_r,LambdaFinal)
  end if
end subroutine Run

subroutine Initialize
#include <petsc/finclude/petsc.h>
  use globals
  use petsc
  implicit none
  integer :: i
  real(rkind) :: t0,t1
  real(rkind), external  :: ppfun
  real(rkind) :: Itot
  real(rkind) :: rmax

  PetscErrorCode :: ierr

  call ReadNamelist

  call Arrays(inds_c)
  call Arrays(inds_r)

  ! This matrix transforms the set of psi and its derivatives and the four corners of a patch into the array aij such that psi = sum_(i,j) aij Z**i R**j (Z and R going from 0 to 1 on a patch) 
  call HermiteMatrix

  ! Gives the answer to the question 'What is iR,iZ and field type if I know ind (from 1 to nws)?'
  call Ind_to_iRiZ(inds_c)
  call Ind_to_iRiZ(inds_r)

  ! Gives the answer to the question 'What is ind (from 1 to nws) if I know iR,iZ and field type?'
  call iRiZ_to_Ind(inds_c)
  call iRiZ_to_Ind(inds_r)

  ! Uses the Hermite matrix to compute all possible aij for all 16 base cases (4 locations on a patch times 4 field types)
  call AllCoeffs(inds_c)
  call AllCoeffs(inds_r)
  
  ! This sets the nkws values that are known due to the boundary conditions
  call PsiBoundaryCondition(inds_c)
  call PsiBoundaryCondition(inds_r)

  ! Initialize petsc
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF = MPI_COMM_SELF  

  ! call cpu_time(t0)
  ! call SETUPAB
  ! call cpu_time(t1)

  ! print*, (t1-t0)/24

  ! Read guess psi
  call ReadGuess

  ! call ReadPprime

  ! Interpolate the guess to obtain a guess with size nws containing also the field derivatives 
  if (guesstype.eq.1) then
     call InterpolatePsi1(inds_g,inds_c)
  elseif (guesstype.eq.2) then
     call InterpolatePsi2(inds_g,inds_c)
     call DeallocateArrays(inds_g)
  end if

  call PsiMaximum(inds_c,inds_c%PsiCur,rmax,psimaxval,.true.)  

  call TotalCurrent(inds_c,inds_c%PsiCur,ppfun,Itot_target)
  Itot_target = Itot_target*LambdaIni

  ! Computing beforehand all possible types of R**i/(R+a) integrals for matrix preparation
  call AllRIntegrals(inds_c)
  call AllRIntegrals(inds_r)

  ! Computing all possible values for the matrix elements
  call Aij_all(inds_c)
  call Aij_all(inds_r)
  ! Assembling the matrix as well as the boundary condition part of the right hand side
  call SETUPAB(inds_c)
  call SETUPAB(inds_r)

  ilast = 1

end subroutine Initialize

subroutine ReadNamelist
  use prec_const
  use globals
  implicit none
  
  integer :: nzc,nrc
  integer :: nz,nr
  real(rkind) :: length,hlength

  integer :: mp

  ! sizes for refined problem
  nz = 101
  nr = 100
  ! sizes for coarse problem
  nzc = 53
  nrc = 32

  length = 1._rkind
  ! psi at the boundary in R=1
  psiedge = -0.5_rkind
  
  ! max number of iterations
  ntsmax = 500
  tol = 1.e-8_rkind

  ! order for gaussian quadratures
  gaussorder = 4
  ! Number of detected points for the plasma vacuum boundary in each square patch 
  nboundarypoints = 5
  
  ! How much percentage is updated in new solution
  relax = 0.1_rkind

  ! Total current (when total current is used)
  Itotal = 2._rkind
  
  ! Target value of psi on magnetic axis
  psimax = 0.05
  
  ! Number of points for psi,theta mesh (jacobian)
  npsi = 50
  ntheta = 64

  ! written with python or with fortran
  ! guesstype = 1 --> solution from Finite Differences version
  ! guesstype = 2 --> solution from the present Finite Element Version
  guesstype = 1

  ! By default 0, in that case is not used
  LambdaNL = 0._rkind

  ! Polynomial defining pprime
  AP_NL = 0._rkind
  AP_NL(1) = 1._rkind

  ! Use petsc, by default yes
  usepetsc = .true.
  
  ! define namelist

  namelist /frc/ &
       &    nzc,nrc,nz,nr,length,psiedge,ntsmax,psimax, &
       &    tol,gaussorder, nboundarypoints,relax,Itotal, &
       &    npsi,ntheta, &
       &    guesstype, usepetsc, &
       &    LambdaNL, AP_NL

  ! read namelist
  mp = 101
  open(mp, file ='nlfrc', delim = 'apostrophe', &
       FORM = 'formatted', action = 'read', status='old')
  read(mp,frc)
  close(mp)

  ! half length
  hlength = 0.5_rkind*length

  inds_r%nz = nz
  inds_r%nr = nr
  inds_r%length = length
  inds_r%hlength = hlength

  inds_c%nz = nzc
  inds_c%nr = nrc
  inds_c%length = length
  inds_c%hlength = hlength

end subroutine ReadNamelist

subroutine Arrays(inds)
  use sizes_indexing
  use prec_const
  use globals, only: ntsmax
  implicit none
  type(indices), intent(inout) :: inds
  integer :: nz,nr
  integer :: nws
  integer :: i

  nZ = inds%nZ
  nR = inds%nR

  allocate(inds%Z(nZ),inds%R(nR), &
       &   inds%RR(nZ,nR),inds%ZZ(nZ,nR), &
       &   inds%logRp1sR(nR), &
       &   inds%PsiBC(nZ,nR,4), &
       &   inds%IndArray(nZ*nR*4,3), &
       &   inds%IndArrayInv(nZ,nR,4), &
       &   inds%RIntegralArray(0:6,nR-1), &
       &   inds%ContributionMat(0:1,0:1,4,0:1,0:1,4,nR-1))

  ! Constructing Z,R arrays
  inds%deltar = 1._rkind/real(nR-1,rkind)
  do i=1,nR
     inds%R(i) = real(i-1,rkind)/real(nR-1,rkind)
     inds%RR(:,i) = real(i-1,rkind)/real(nR-1,rkind)
  end do

  inds%deltaz = 1._rkind/real(nZ-1,rkind)*inds%length
  do i=1,nZ
     inds%Z(i) = real(i-1,rkind)/real(nZ-1,rkind)*inds%length - inds%hlength
     inds%ZZ(i,:) = real(i-1,rkind)/real(nZ-1,rkind)*inds%length - inds%hlength
  end do
  
  ! dRdZ array is used to transform quantities like dpsi/dZ into dpsi/dZbar, where Zbar goes from 0 to 1 within a patch
  inds%dRdZ(1) = 1._rkind
  inds%dRdZ(2) = inds%deltaz
  inds%dRdZ(3) = inds%deltar
  inds%dRdZ(4) = inds%deltar*inds%deltaz

  ! log(1+R/dR)/(R/dR) is involved in the radial integrations when building the matrix
  inds%logRp1sR(1) = 0._rkind
  do i=2,nR
     inds%logRp1sR(i) = log((1._rkind+inds%R(i)/inds%deltar)/(inds%R(i)/inds%deltar))
  end do

  ! (nR-2)*(nZ-2)*4 : bulk values
  ! 2*(nZ-2) : 2 unknowns per R=1 edge point, with a total of nZ points (psi and dpsi/dZ are known)
  ! 4*(nR-2) : 2 unknowns per Z=+-L/2 edge point (psi, dRpsi are known)
  ! 2 : 1 unknown (d2psi/dRdZ) per top corner
  ! No unknowns on bottom edge (including corners)
  nws = (nR-2)*(nZ-2)*4 + 2*(nZ-2)+4*(nR-2)+2 ! number of unknowns
  inds%nws = nws
  inds%nkws = nR*nZ*4-nws ! number of known values (boundary conditions)
  inds%ntot = inds%nws+inds%nkws ! total size of array needed to reconstruct psi

  allocate(inds%B_BC(nws),      &
       &   inds%PsiFinal(nws),  &
       &   inds%rhs(nws),       &
       &   inds%AllPsis(nws,ntsmax+1))

  if(.not.allocated(inds%PsiCur)) allocate(inds%PsiCur(nws))

end subroutine Arrays

subroutine DeallocateArrays(inds)
  use sizes_indexing
  implicit none
  type(indices), intent(inout) :: inds
  integer :: nz,nr
  integer :: i

  deallocate(inds%Z,inds%R,inds%ZZ,inds%RR, &
       &     inds%logRp1sR,inds%PsiBC,inds%IndArray,inds%IndArrayInv, &
       &     inds%RIntegralArray,inds%ContributionMat, &
       &     inds%B_BC,inds%PsiCur,inds%PsiFinal,inds%rhs,inds%AllPsis)  

end subroutine DeallocateArrays

function ppfun(psiv) result(result_val)
  ! Pprime function in Grad-Shafranov equation
  use prec_const
  use globals, only : AP_NL,psimaxval
  implicit none
  real(rkind), intent(in) :: psiv
  real(rkind) :: x, result_val
  integer :: i

  x = psiv/psimaxval

  result_val = 0._rkind
  do i=1,10
     result_val = result_val + AP_NL(i)*x**(i-1)
  end do

end function ppfun

subroutine PsiBoundaryCondition(inds)
  ! This sets the nkws values that are known due to the boundary conditions
  ! psi = R**2*psiedge
  use globals, only : psiedge
  use prec_const
  use sizes_indexing
  implicit none
  type(indices), intent(inout) :: inds
  integer :: iZ,iR

  inds%psiBC = 0._rkind

  do iR=1,inds%nR-1
     Inds%PsiBC(1,iR,1) = psiedge*inds%R(iR)**2
     Inds%PsiBC(1,iR,3) = 2._rkind*psiedge*inds%R(iR)
     Inds%PsiBC(inds%nZ,iR,1) = psiedge*inds%R(iR)**2
     Inds%PsiBC(inds%nZ,iR,3) = 2._rkind*psiedge*inds%R(iR)
  end do

  do iZ=2,inds%nZ-1
     Inds%PsiBC(iZ,inds%nR,1) = psiedge
  end do

  Inds%PsiBC(1,inds%nR,1) = psiedge
  Inds%PsiBC(1,inds%nR,3) = 2._rkind*psiedge
  Inds%PsiBC(inds%nZ,inds%nR,1) = psiedge
  Inds%PsiBC(inds%nZ,inds%nR,3) = 2._rkind*psiedge

end subroutine PsiBoundaryCondition

subroutine HermiteMatrix
  ! This matrix transforms the set of psi and its derivatives and the four corners of a patch into the array aij such that psi = sum_(i,j) aij Z**i R**j (Z and R going from 0 to 1 on a patch) 
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

subroutine InterpolatePsi1(inds_g,inds)
  ! Transform a 2D psi array into an array adapted to the finite element
  ! Psizrg is typically a 2D guess z,r array
  ! Psi is in the form of a size nws array
  use prec_const
  use sizes_indexing
  use globals, only : PsiGuess1
  implicit none
  type(indices), intent(inout) :: inds_g,inds
  real(rkind), dimension(inds_g%nZ) :: zg
  real(rkind), dimension(inds_g%nR) :: rg
  real(rkind) :: dzg,drg
  real(rkind) :: z0,z1,r0,r1
  real(rkind) :: length
  integer :: nzg,nrg,nws
  integer :: ind
  integer :: iZ,iR,k
  integer :: izg,irg
  real(rkind) :: t,u
  real(rkind) :: zind,rind
  real(rkind) :: psi00,psi01,psi10,psi11

  length = inds_g%length

  z0 = -0.5_rkind*inds_g%length
  z1 = 0.5_rkind*inds_g%length

  r0 = 0._rkind
  r1 = 1._rkind

  nzg = inds_g%nZ
  nrg = inds_g%nR
  nws = inds%nws

  call linspace(z0,z1,nzg,zg)
  call linspace(r0,r1,nrg,rg)

  dzg = 1._rkind/real(nzg-1,rkind)*length
  drg = 1._rkind/real(nrg-1,rkind)

  do ind=1,nws
     iZ = inds%IndArray(ind,1)
     iR = inds%IndArray(ind,2)
     k = inds%IndArray(ind,3)

     zind = (inds%Z(iZ)+0.5_rkind*inds%length)/dzg
     rind = inds%R(iR)/drg

     izg = min(int(zind),nzg-2)
     irg = min(int(rind),nrg-2)
     t = zind - real(izg,rkind)
     u = rind - real(irg,rkind)

     izg = izg+1
     irg = irg+1

     psi00 = PsiGuess1(izg  ,irg  )
     psi10 = PsiGuess1(izg+1,irg  )
     psi01 = PsiGuess1(izg  ,irg+1)
     psi11 = PsiGuess1(izg+1,irg+1)

     if (k.eq.1) then
        ! psi
        inds%PsiCur(ind) = (psi00*(1._rkind-t)*(1._rkind-u) + &
             &      psi10*t*(1._rkind-u) +  &
             &      psi01*(1._rkind-t)*u + &
             &      psi11*t*u)
     else if (k.eq.2) then
        ! dpsi/dz
        inds%PsiCur(ind) = (psi10 - psi00)*(1._rkind-u)/dzg + (psi11 - psi01)*u/dzg
     else if (k.eq.3) then
        ! dpsi/dr
        inds%PsiCur(ind) = (psi01 - psi00)*(1._rkind-t)/drg + (psi11 - psi10)*t/drg
     else if (k.eq.4) then
        ! d2psi/dzdr
        inds%PsiCur(ind) = (psi11 + psi00 - psi10 - psi01)/(drg*dzg)
     end if
  end do

  deallocate(PsiGuess1)

end subroutine InterpolatePsi1

subroutine InterpolatePsi2(inds1,inds2)
  use prec_const
  use sizes_indexing
  implicit none
  type(indices), intent(inout) :: inds1,inds2
  real(rkind), dimension(inds1%ntot) :: PsiAll
  integer :: ind
  integer :: iZ,iR,k
  real(rkind) :: zv,rv
  real(rkind) ::psival

  call FillPsiAll(inds1,inds1%PsiCur,PsiAll)

  do ind=1,inds2%nws
     iZ = inds2%IndArray(ind,1)
     iR = inds2%IndArray(ind,2)
     k = inds2%IndArray(ind,3)
     
     zv = inds2%Z(iZ)
     rv = inds2%R(iR)

     call EvalPsi(inds1,PsiAll,zv,rv,psival,k)

     inds2%PsiCur(ind) = psival

  end do

end subroutine InterpolatePsi2

subroutine Evolve(inds,PsiIni,fun,Psis,nts)
  ! Simple Picard algorithm to find psi
  use prec_const
  use sizes_indexing
  use globals, only : tol, LambdaIni, LambdaFinal, psimax, ilast
  implicit none
  type(indices), intent(inout) :: inds
  real(rkind), dimension(inds%nws), intent(in) :: PsiIni
  integer, intent(in) :: nts
  real(rkind), dimension(inds%nws,nts+1), intent(out) :: Psis
  real(rkind), dimension(nts+1) :: Lambdas
  real(rkind), external :: fun
  
  real(rkind) :: Lambda
  real(rkind), dimension(inds%nws) :: Psi1,Psi2
  real(rkind) :: ratio,n_,n_star
  real(rkind) :: Itot
  real(rkind) :: dtpsi
  real(rkind) :: Icur
  real(rkind) :: rmax,psimaxcur
  logical :: advance
  integer :: i,k
  logical :: equatorial

  ! if nz is odd, magnetic axis ought to be located in the equatorial plane z=0
  equatorial = .false.
  if (mod(inds%nZ,2).eq.1) equatorial = .true.

  Psis = 0._rkind
  Lambdas = 0._rkind

  i = 1
  dtpsi = tol
  Psis(:,i) = PsiIni
  Psi1(:) = PsiIni(:)

  call PsiMaximum(inds,PsiIni,rmax,psimaxcur,equatorial)
  
  Lambdas(i) = Lambdaini
  ! Lambda is a Lagrange multiplier to respect the constraint of maximum psi
  Lambda = LambdaIni
  write(*,'(A,E12.4)'), 'Initial Lambda =', Lambda
  ! call TotalCurrent(PsiIni,fun,Itot)
  ! Itot = Itotal

  write(*,'(A)'), ' Iteration      Lambda       dtpsi'
  ! Do at least 10 iterations
  do while ((i.le.nts .and. dtpsi.ge.tol) .or. (i.le.10))
     ! One step is to solve A psinew = rhs(psiold) 
     ! and update psi2 = psiold + relax*(psinew-psiold)
     call FEMStep(inds,Psi1,funLambda,Psi2)
     ! Compute norm difference between old and new solutions
     call DiffNorm(inds%nws,Psi1,Psi2,dtpsi)
     ! call TotalCurrent(Psi2,fun,Icur)
     ! Compute Psimax
     call PsiMaximum(inds,Psi2,rmax,psimaxcur,equatorial)
     ! Lambda = Itot/Icur
     ! Update the Lagrange multiplier
     Lambda = Lambda*psimax/psimaxcur
     write(*,'(I10,2E12.4)'), i,Lambda,dtpsi
     ! Update psi
     i = i+1
     Psis(:,i) = Psi2(:)
     Psi1(:) = Psi2(:)
  end do

  if (dtpsi.ge.tol) then
     print*, 'Maximum number of iterations reached'
  end if

  ilast = i
  LambdaFinal = Lambda

  ! call TotalCurrent(Psi2,funLambda,Icur)

  ! write(*,'(A,2E12.4)') 'Total current at convergence = ', Icur

contains
  
  function funLambda(psi)
    ! return ppfun times the Lagrange multiplier Lambda
    implicit none
    real(rkind) :: psi,funLambda

    funLambda = fun(psi)*Lambda

  end function funLambda
  
end subroutine Evolve

subroutine FEMStep(inds,Psi,fun,PsiSol)
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscksp.h>
  use prec_const
  use sizes_indexing
  use petscksp
  use globals, only : relax,psimaxval
  implicit none
  type(indices), intent(inout) :: inds

  real(rkind), dimension(inds%nws), intent(in) :: Psi
  real(rkind), dimension(inds%nws), intent(out) :: PsiSol
  real(rkind), external :: fun
  integer :: ind
  real(rkind) :: norm2,norm2_
  real(rkind) :: rmax

  Vec :: u,x,b
  PetscInt :: pnws
  PetscErrorCode :: ierr
  Mat :: PM
  PetscReal :: one,negdt
  PetscReal :: ptol
  PetscInt :: ione
  PetscInt :: II
  PetscInt, dimension(0:inds%nws-1) :: ix
  PetscReal, dimension(0:inds%nws-1) :: pvals

  PC :: pc
  KSP :: ksp

  call PsiMaximum(inds,Psi,rmax,psimaxval,.true.)
  ! Compute psi dependent part of the right hand side (integrals of R*pprime(psi) )
  call RightHandSide(inds,Psi,fun,inds%rhs)

  ! Set right hand side
  pnws = inds%nws
  do II=0,pnws-1
     ix(II) = II
     pvals(II) = (inds%B_BC(II+1)+inds%rhs(II+1))
  end do
  ione = 1
  one = 1.d0

  ! Create KSP context for linear solve
  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
  call KSPSetOperators(ksp,inds%PAMat,inds%PAMat,ierr)

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

  ! Solve the linear problem
  call KSPSolve(ksp,b,x,ierr)

  ! Get the values of the solution
  call VecGetValues(x,pnws,ix,pvals,ierr)

  ! Use them to update the solution with relax
  PsiSol(1:inds%nws) = Psi(1:inds%nws)*(1._rkind - relax) + relax*pvals(0:pnws-1)

  ! Destroy PETSc vectors and ksp context
  call VecDestroy(x,ierr)
  call VecDestroy(b,ierr)
  call KSPDestroy(ksp,ierr)
  
end subroutine FEMStep

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

subroutine calculate_Vprime(JacobMesh, Vprime)
  use globals, only : npsi, ntheta
  implicit none

  real, dimension(:,:), intent(in) :: JacobMesh
  real, dimension(:), intent(out) :: Vprime

  integer :: i

  do i = 1, npsi
     Vprime(i) = sum(JacobMesh(i, 1:ntheta-1)) / real(ntheta-1)
  end do
end subroutine calculate_Vprime

subroutine calculate_AP_NL(psi, P, AP_NL)
  use globals, only : npsi
  implicit none
  
  ! Input parameters
  real(kind=8), dimension(npsi), intent(in) :: psi
  real(kind=8), dimension(npsi), intent(in) :: P
  
  ! Output parameter
  real(kind=8), dimension(10), intent(out) :: AP_NL
  
  ! Local variables
  real(kind=8), dimension(npsi) :: Pprime
  real(kind=8), dimension(0:9) :: poly_coeffs
  integer :: i, j
  
  ! Calculate derivative of P with respect to psi
  do i = 2, npsi-1
    Pprime(i) = (P(i+1) - P(i-1)) / (psi(i+1) - psi(i-1))
  end do
  
  ! Handle edge cases for Pprime at psi[1] and psi[npsi]
  Pprime(1) = (P(2) - P(1)) / (psi(2) - psi(1))
  Pprime(npsi) = (P(npsi) - P(npsi-1)) / (psi(npsi) - psi(npsi-1))
  
  ! Fit Pprime to a 9th degree polynomial using polyfit
  call polyfit(psi, Pprime, 9, AP_NL)
  
  ! Store polynomial coefficients as AP_NL
  do j = 0, 9
    AP_NL(j+1) = poly_coeffs(j)
  end do
  
end subroutine calculate_AP_NL


subroutine Finalize
#include <petsc/finclude/petsc.h>
  use globals
  use petsc
  implicit none

  PetscErrorCode :: ierr

  call MatDestroy(inds_c%PAMat,ierr)
  call MatDestroy(inds_r%PAMat,ierr)

  call PetscFinalize(PETSC_NULL_CHARACTER,ierr)

  call DeallocateArrays(inds_c)
  call DeallocateArrays(inds_r)

  ! deallocate(Xpprime,Ypprime,Bpprime,Cpprime,Dpprime)
  ! deallocate(Apprime)

end subroutine Finalize
