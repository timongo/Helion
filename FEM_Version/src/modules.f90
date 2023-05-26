module sizes_indexing
#include <petsc/finclude/petscmat.h>
  use petsc
  use prec_const
  type :: indices
     integer :: nz,nr
     integer :: nws,nkws
     integer :: ntot
     real(rkind) :: length,hlength
     real(rkind) :: deltar,deltaz
     real(rkind), dimension(4) :: dRdZ
     real(rkind), allocatable, dimension(:) :: Z,R
     real(rkind), allocatable, dimension(:,:) :: ZZ,RR
     integer, allocatable, dimension(:,:) :: IndArray
     integer, allocatable, dimension(:,:,:) :: IndArrayInv
     real(rkind), allocatable, dimension(:) :: logRp1sR
     real(rkind), allocatable, dimension(:,:,:) :: PsiBC
     real(rkind), allocatable, dimension(:,:) :: RIntegralArray
     real(rkind), allocatable, dimension(:,:,:,:,:,:,:) :: ContributionMat
     real(rkind), allocatable, dimension(:) :: B_BC
     real(rkind), allocatable, dimension(:,:) :: AllPsis
     real(rkind), allocatable, dimension(:) :: PsiCur
     real(rkind), allocatable, dimension(:) :: PsiFinal
     real(rkind), allocatable, dimension(:) :: rhs
     real(rkind), dimension(4,4,0:1,0:1,4) :: CoeffMat
     Mat      :: PAMat
  end type indices
end module sizes_indexing

module globals
  use prec_const
  use sizes_indexing

  ! container for indices of the coarse problem
  type(indices) :: inds_c
  ! container for indices for refined problem
  type(indices) :: inds_r
  ! container for indices of the guess
  type(indices) :: inds_g
  
  real(rkind) :: psiedge
  real(rkind) :: current
  integer :: gaussorder
  integer :: nboundarypoints

  integer :: ntsmax
  integer :: ilast
  real(rkind) :: tol
  real(rkind) :: relax
  real(rkind) :: Itotal
  real(rkind) :: psimax,psimaxval,psimaxcur
  real(rkind) :: Itot_target

  integer :: npsi
  integer :: ntheta
  real(rkind) :: LagMulC

  real(rkind) :: theta1,theta2,theta3,theta4

  integer :: guesstype

  real(rkind), allocatable, dimension(:,:) :: PsiGuess1
  real(rkind), dimension(16,16) :: Hmat
  real(rkind) :: LambdaIni,LambdaFinal,LambdaNL
  real(rkind), dimension(10) :: AP_NL

  integer :: npprime
  real(rkind), allocatable, dimension(:) :: Apprime
  ! real(rkind), allocatable, dimension(:) :: Xpprime
  ! real(rkind), allocatable, dimension(:) :: Ypprime
  ! real(rkind), allocatable, dimension(:) :: Bpprime
  ! real(rkind), allocatable, dimension(:) :: Cpprime
  ! real(rkind), allocatable, dimension(:) :: Dpprime

  integer,save :: count=0

  logical :: usepetsc
  
end module globals

