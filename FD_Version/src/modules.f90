module globals
  use prec_const

  integer :: nr,nz
  integer :: nrtot,nztot
  integer :: ntot
  integer :: nrguess,nzguess
  integer :: nrguesstot,nzguesstot
  real(rkind) :: lguess
  real(rkind) :: length
  real(rkind) :: psiedge
  real(rkind) :: current

  real(rkind) :: pprime
  real(rkind) :: omegai
  integer :: ntsmax
  real(rkind) :: tol
  real(rkind) :: psimax

  real(rkind) :: deltar,deltaz
  real(rkind) :: drguess,dzguess
  real(rkind), allocatable, dimension(:) :: Z,R
  real(rkind), allocatable, dimension(:,:) :: ZZ,RR
  real(rkind), allocatable, dimension(:,:) :: psi
  real(rkind), allocatable, dimension(:,:) :: psiguess
  real(rkind), allocatable, dimension(:,:) :: psiguess_read

end module globals

