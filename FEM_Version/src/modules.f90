module globals
#include <petsc/finclude/petscksp.h>
  use prec_const
  use petsc

  integer :: nz,nr
  integer :: nws,nkws
  integer :: ntot
  integer :: nzguess,nrguess
  real(rkind) :: lguess
  real(rkind) :: length,hlength
  real(rkind) :: psiedge
  real(rkind) :: current
  integer :: gaussorder
  integer :: nboundarypoints

  integer :: ntsmax
  integer :: ilast
  real(rkind) :: tol
  real(rkind) :: relax
  real(rkind) :: Itotal
  real(rkind) :: psimax

  integer :: npsi
  integer :: ntheta
  real(rkind) :: LagMulC

  real(rkind) :: deltar,deltaz
  real(rkind) :: drguess,dzguess
  real(rkind) :: theta1,theta2,theta3,theta4

  real(rkind), allocatable, dimension(:,:) :: PsiGuess
  real(rkind), allocatable, dimension(:) :: R,Z
  real(rkind), allocatable, dimension(:) :: logRp1sR
  real(rkind), dimension(4) :: dRdZ
  real(rkind), dimension(16,16) :: Hmat
  real(rkind), allocatable, dimension(:,:) :: RR,ZZ
  integer, allocatable, dimension(:,:) :: IndArray
  integer, allocatable, dimension(:,:,:) :: IndArrayInv
  real(rkind), dimension(4,4,0:1,0:1,4) :: CoeffMat
  real(rkind), allocatable, dimension(:,:,:) :: PsiBC
  real(rkind), allocatable, dimension(:,:) :: RIntegralArray
  real(rkind), allocatable, dimension(:,:,:,:,:,:,:) :: ContributionMat
  real(rkind), allocatable, dimension(:) :: B_BC
  real(rkind), allocatable, dimension(:) :: PsiCur
  real(rkind), allocatable, dimension(:) :: PsiFinal
  real(rkind) :: LambdaFinal
  real(rkind), allocatable, dimension(:,:) :: AllPsis
  real(rkind), allocatable, dimension(:) :: rhs

  integer :: npprime
  real(rkind), allocatable, dimension(:) :: Apprime
  ! real(rkind), allocatable, dimension(:) :: Xpprime
  ! real(rkind), allocatable, dimension(:) :: Ypprime
  ! real(rkind), allocatable, dimension(:) :: Bpprime
  ! real(rkind), allocatable, dimension(:) :: Cpprime
  ! real(rkind), allocatable, dimension(:) :: Dpprime

  Mat      :: PAMat
  
end module globals

