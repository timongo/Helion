MODULE prec_const
  !
  !   Precision for real and complex
  !
  INTEGER, PARAMETER :: IKIND = kind(1)
  INTEGER, PARAMETER :: IKIND8 = 2*IKIND
  INTEGER, PARAMETER :: RKIND = kind(1.d0)
  INTEGER, PARAMETER :: RKIND4 = kind(1.e0)
  INTEGER, PARAMETER :: CKIND = 2*RKIND
  INTEGER, PARAMETER :: CKIND4 = 2*RKIND4
  !
  !   Some useful constants
  !
  REAL(RKIND), PARAMETER :: PI=3.141592653589793238462643383279502884197_rkind
  REAL(RKIND), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_rkind
  REAL(RKIND), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_rkind
  REAL(RKIND), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_rkind
  !
  ! Physical constants
  !
  real(rkind),parameter :: kb = 1.3806504e-23_rkind
  real(rkind),parameter :: mproton = 1.672623e-27_rkind
  real(rkind),parameter :: qproton =  1.60217653e-19_rkind
  real(rkind),parameter :: melectron   = 9.10e-31_rkind
  real(rkind),parameter :: mu0 = 4.0_rkind*PI*1.e-7_rkind
  !
END MODULE prec_const
