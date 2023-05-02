subroutine PsiIntegral(psi,res)
  use globals
  implicit none
  real(rkind), dimension(4,4), intent(inout) :: psi
  real(rkind), intent(out) :: res
  real(rkind) :: x0,x1
  real(rkind), external :: zero,one

  x0 = 0._rkind
  x1 = 1._rkind

  call dblquad(fzr,x0,x1,zero,one,res)

  contains
    function fzr(z,r)
      real(rkind) :: fzr,z,r
      integer :: i,j
      
      fzr = 0._rkind
      do i=0,3
         do j=0,3
            fzr = fzr + psi(i+1,j+1)*z**i*r**j
         end do
      end do
         
    end function fzr

end subroutine PsiIntegral

function zero(z)
  use prec_const
  implicit none
  real(rkind) :: zero,z
  
  zero = 0._rkind
  
end function zero

function one(z)
  use prec_const
  implicit none
  real(rkind) :: one,z
  
  one = 1._rkind
  
end function one

subroutine dblquad(fxy,x0,x1,y0,y1,q)
  use globals
  implicit none
  real(rkind), intent(in) :: x0,x1
  real(rkind), external :: fxy,y0,y1
  real(rkind), intent(out) :: q
  
  call quad(fx,x0,x1,q)

  contains
    
    function fx(x)
      real(rkind) x,fx
      real(rkind) y0_,y1_,res

      y0_ = y0(x)
      y1_ = y1(x)

      call redquad(fxy,x,y0_,y1_,res)

      fx = res

    end function fx  

end subroutine dblquad

subroutine redquad(fxy,x,y0,y1,q)
  implicit none
  real(kind = 8),external :: fxy
  real(kind = 8), intent(in) :: x,y0,y1
  real(kind = 8), intent(out) :: q

  real ( kind = 8 ), parameter :: epsabs = 0.0E+00
  real ( kind = 8 ), parameter :: epsrel = 0.001E+00
  real ( kind = 8 ) abserr
  integer ( kind = 4 ), parameter :: limit = 500
  integer ( kind = 4 ), parameter :: lenw = limit * 4
  integer ( kind = 4 ) :: last
  integer ( kind = 4 ) :: ier
  integer ( kind = 4 ) :: iwork(limit)
  integer ( kind = 4 ) :: neval
  real ( kind = 8 )    :: work(lenw)

  call dqags(fy,y0,y1,epsabs,epsrel,q,abserr,neval,ier, &
       &     limit, lenw, last, iwork, work)

  contains

    function fy(y)
      implicit none
      real(kind = 8) y,fy

      fy = fxy(x,y)

    end function fy

end subroutine redquad

subroutine quad(f,x0,x1,q)
  implicit none
  real(kind = 8), external :: f
  real(kind = 8) :: x0,x1,q

  real ( kind = 8 ), parameter :: epsabs = 0.0E+00
  real ( kind = 8 ), parameter :: epsrel = 0.000001E+00
  real ( kind = 8 ) abserr
  integer ( kind = 4 ), parameter :: limit = 500
  integer ( kind = 4 ), parameter :: lenw = limit * 4
  integer ( kind = 4 ) :: last
  integer ( kind = 4 ) :: ier
  integer ( kind = 4 ) :: iwork(limit)
  integer ( kind = 4 ) :: neval
  real ( kind = 8 )    :: work(lenw)

  call dqags(f,x0,x1,epsabs,epsrel,q,abserr,neval,ier, &
       &     limit, lenw, last, iwork, work)

end subroutine quad
