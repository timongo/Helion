module globals
  use prec_const

  integer :: nr,nz
  integer :: nrtot,nztot
  integer :: ntot
  integer :: nrguess,nzguess
  integer :: nrguesstot,nzguesstot
  real(rkind) :: lguess
  real(rkind) :: length,hlength
  real(rkind) :: psiedge
  real(rkind) :: current,CurrentTarget

  real(rkind) :: omegai
  integer :: ntsmax
  real(rkind) :: tol
  real(rkind) :: psimax,psimaxcur

  real(rkind) :: LambdaIni,LambdaSol,LambdaNL

  integer :: guesstype

  integer :: iplot

  logical :: ControlCurrent

  integer :: npsi,ntheta

  real(rkind) :: deltar,deltaz
  real(rkind) :: drguess,dzguess
  real(rkind) :: theta1,theta2,theta3,theta4
  real(rkind), allocatable, dimension(:) :: Z,R
  real(rkind), allocatable, dimension(:,:) :: ZZ,RR
  real(rkind), allocatable, dimension(:,:) :: psi
  real(rkind), allocatable, dimension(:,:) :: psiguess
  real(rkind), allocatable, dimension(:,:) :: psiguess_read
  real(rkind), dimension(10) :: AP
  real(rkind), dimension(10) :: AP_NL
  real(rkind), dimension(10) :: AP_guess

end module globals

