subroutine HeavisideQuad(fzr,gzr,n,res)
  ! fzr is the function to integrate on [0,1]x[0,1]
  ! over the domain where gzr is positive
  use prec_const
  implicit none
  real(rkind), external :: fzr,gzr
  integer, intent(in) :: n
  real(rkind), intent(out) :: res
  
  real(rkind) :: g00,g01,g10,g11
  integer :: case

  ! Compute values of gzr at four corners
  g00 = gzr(0._rkind,0._rkind)
  g01 = gzr(0._rkind,1._rkind)
  g10 = gzr(1._rkind,0._rkind)
  g11 = gzr(1._rkind,1._rkind)

  ! Use these values to determine which case is relevant
  call WhichCase(g00,g01,g10,g11,case)
  
  ! Compute the integral using quadpack
  call Case_1to12(fzr,gzr,case,n,res)

end subroutine HeavisideQuad

subroutine WhichCase(g00,g01,g10,g11,case)
  use prec_const
  implicit none
  real(rkind), intent(in) :: g00,g01,g10,g11
  integer, intent(out) :: case

  case = 0

  if (g10.ge.0 .and. g00.le.0 .and. g01.le.0 .and. g11.le.0) then
     case=1
  else if (g11.ge.0 .and. g00.le.0 .and. g01.le.0 .and. g10.le.0) then
     case=2
  else if (g01.ge.0 .and. g00.le.0 .and. g10.le.0 .and. g11.le.0) then
     case=3
  else if (g00.ge.0 .and. g01.le.0 .and. g10.le.0 .and. g11.le.0) then
     case=4
  else if (g00.ge.0 .and. g10.ge.0 .and. g01.le.0 .and. g11.le.0) then
     case=5
  else if (g01.ge.0 .and. g11.ge.0 .and. g00.le.0 .and. g10.le.0) then
     case=6
  else if (g10.ge.0 .and. g11.ge.0 .and. g00.le.0 .and. g01.le.0) then
     case=7
  else if (g00.ge.0 .and. g01.ge.0 .and. g10.le.0 .and. g11.le.0) then
     case=8
  else if (g10.le.0 .and. g01.ge.0 .and. g00.ge.0 .and. g11.ge.0) then
     case=9
  else if (g11.le.0 .and. g01.ge.0 .and. g10.ge.0 .and. g00.ge.0) then
     case=10
  else if (g01.le.0 .and. g00.ge.0 .and. g10.ge.0 .and. g11.ge.0) then
     case=11
  else if (g00.le.0 .and. g01.ge.0 .and. g10.ge.0 .and. g11.ge.0) then
     case=12
  end if

end subroutine WhichCase

subroutine Case_1to12(fzr,gzr,case,n,res)
  use prec_const
  implicit none
  real(rkind), external :: fzr,gzr,zero,one
  integer, intent(in) :: case
  integer, intent(in) :: n
  real(rkind), intent(out) :: res
  real(rkind) :: res2
  real(rkind), dimension(n) :: zs,rs
  real(rkind), external :: zbrent

  real(rkind), parameter :: ztol=1.e-8
  real(rkind) :: z0,z1
  real(rkind) :: z0_,z1_
  real(rkind) :: zval
  real(rkind) :: r0,r1
  real(rkind) :: rval
  integer :: i

  if (case.eq.1) then

     z0 = zbrent(gz0,0._rkind,1._rkind,ztol)
     z1 = 1._rkind

     call linspace(z0,z1,n,zs)
     rs = 0._rkind

     do i=2,n
        zval = zs(i)
        rs(i) = zbrent(gz0r,0._rkind,1._rkind,ztol)
     end do

     call dblquad(fzr,z0,z1,zero,rsinterp,res)

  else if (case.eq.2) then

     z0 = zbrent(gz1,0._rkind,1._rkind,ztol)
     z1 = 1._rkind

     call linspace(z0,z1,n,zs)
     rs = 0._rkind
     rs(1) = 1._rkind

     do i=2,n
        zval = zs(i)
        rs(i) = zbrent(gz0r,0._rkind,1._rkind,ztol)
     end do

     call dblquad(fzr,z0,z1,rsinterp,one,res)

  else if (case.eq.3) then

     z0 = 0._rkind
     z1 = zbrent(gz1,0._rkind,1._rkind,ztol)

     call linspace(z0,z1,n,zs)
     rs = 0._rkind
     rs(n) = 1._rkind

     do i=1,n-1
        zval = zs(i)
        rs(i) = zbrent(gz0r,0._rkind,1._rkind,ztol)
     end do

     call dblquad(fzr,z0,z1,rsinterp,one,res)

  else if (case.eq.4) then

     z0 = 0._rkind
     z1 = zbrent(gz0,0._rkind,1._rkind,ztol)

     call linspace(z0,z1,n,zs)
     rs = 0._rkind

     do i=1,n-1
        zval = zs(i)
        rs(i) = zbrent(gz0r,0._rkind,1._rkind,ztol)
     end do

     call dblquad(fzr,z0,z1,zero,rsinterp,res)

  else if (case.eq.5) then

     z0 = 0._rkind
     z1 = 1._rkind

     call linspace(z0,z1,n,zs)
     rs = 0._rkind

     do i=1,n
        zval = zs(i)
        rs(i) = zbrent(gz0r,0._rkind,1._rkind,ztol)
     end do

     call dblquad(fzr,z0,z1,zero,rsinterp,res)

  else if (case.eq.6) then

     z0 = 0._rkind
     z1 = 1._rkind

     call linspace(z0,z1,n,zs)
     rs = 0._rkind

     do i=1,n
        zval = zs(i)
        rs(i) = zbrent(gz0r,0._rkind,1._rkind,ztol)
     end do

     call dblquad(fzr,z0,z1,rsinterp,one,res)

  else if (case.eq.7) then

     r0 = 0._rkind
     r1 = 1._rkind

     zs = 0._rkind
     call linspace(r0,r1,n,rs)

     do i=1,n
        rval = rs(i)
        zs(i) = zbrent(gzr0,0._rkind,1._rkind,ztol)
     end do

     call dblquad(frz,r0,r1,zsinterp,one,res)

  else if (case.eq.8) then

     r0 = 0._rkind
     r1 = 1._rkind

     zs = 0._rkind
     call linspace(r0,r1,n,rs)

     do i=1,n
        rval = rs(i)
        zs(i) = zbrent(gzr0,0._rkind,1._rkind,ztol)
     end do

     call dblquad(frz,r0,r1,zero,zsinterp,res)

  else if (case.eq.9) then

     z0 = 0._rkind
     z0_ = zbrent(gz0,0._rkind,1._rkind,ztol)
     z1 = 1._rkind

     call linspace(z0_,z1,n,zs)
     rs = 0._rkind

     do i=2,n
        zval = zs(i)
        rs(i) = zbrent(gz0r,0._rkind,1._rkind,ztol)
     end do

     call dblquad(fzr,z0,z0_,zero,one,res)
     call dblquad(fzr,z0_,z1,rsinterp,one,res2)

     res = res+res2

  else if (case.eq.10) then

     z0 = 0._rkind
     z0_ = zbrent(gz1,0._rkind,1._rkind,ztol)
     z1 = 1._rkind

     call linspace(z0_,z1,n,zs)
     rs = 0._rkind
     rs(1) = 1._rkind

     do i=2,n
        zval = zs(i)
        rs(i) = zbrent(gz0r,0._rkind,1._rkind,ztol)
     end do

     call dblquad(fzr,z0,z0_,zero,one,res)
     call dblquad(fzr,z0_,z1,zero,rsinterp,res2)

     res = res+res2

  else if (case.eq.11) then

     z0 = 0._rkind
     z0_ = zbrent(gz1,0._rkind,1._rkind,ztol)
     z1 = 1._rkind

     call linspace(z0,z0_,n,zs)
     rs = 0._rkind
     rs(n) = 1._rkind

     do i=1,n-1
        zval = zs(i)
        rs(i) = zbrent(gz0r,0._rkind,1._rkind,ztol)
     end do

     call dblquad(fzr,z0,z0_,zero,rsinterp,res)
     call dblquad(fzr,z0_,z1,zero,one,res2)

     res = res+res2

  else if (case.eq.12) then

     z0 = 0._rkind
     z0_ = zbrent(gz0,0._rkind,1._rkind,ztol)
     z1 = 1._rkind

     call linspace(z0,z0_,n,zs)
     rs = 0._rkind
     rs(n) = 0._rkind

     do i=1,n-1
        zval = zs(i)
        rs(i) = zbrent(gz0r,0._rkind,1._rkind,ztol)
     end do

     call dblquad(fzr,z0,z0_,rsinterp,one,res)
     call dblquad(fzr,z0_,z1,zero,one,res2)

     res = res+res2

  end if

contains    

  function gz0(z)
    real(rkind) :: z,gz0
    gz0 = gzr(z,0._rkind)
  end function gz0

  function gz1(z)
    real(rkind) :: z,gz1
    gz1 = gzr(z,1._rkind)
  end function gz1

  function gz0r(r)
    real(rkind) :: r,gz0r
    gz0r = gzr(zval,r)
  end function gz0r

  function gzr0(z)
    real(rkind) :: z,gzr0
    gzr0 = gzr(z,rval)
  end function gzr0
  
  function frz(r,z)
    ! simply inverts the arguments (for cases 7 and 8)
    real(rkind) :: r,z,frz
    frz = fzr(z,r)
  end function frz

  function rsinterp(z) result(res)
    real(rkind) :: z,res

    if (n.eq.2) then
       ! if only 2 points, linear interpolation is used
       res = rs(1) + (z-z0)/(z1-z0)*(rs(2)-rs(1))
    else
       ! otherwise, a quadratic interpolation is performed
       call interpolation(zs,rs,n,z,res)
    end if
  end function rsinterp

  function zsinterp(r) result(res)
    real(rkind) :: r,res

    if (n.eq.2) then
       ! if only 2 points, linear interpolation is used
       res = zs(1) + (r-r0)/(r1-r0)*(zs(2)-zs(1))
    else
       ! otherwise, a quadratic interpolation is performed
       call interpolation(rs,zs,n,r,res)
    end if
  end function zsinterp

end subroutine Case_1to12

function zbrent(func,x1,x2,tol) 
  use prec_const
  implicit none
  real(rkind) :: zbrent,tol,x1,x2,func
  EXTERNAL func
  integer, parameter :: itmax=100
  real(rkind), parameter :: eps=3.e-8
  ! Using Brent’s method, find the root of a function func known to lie between x1 and x2. The root, returned as zbrent, will be refined until its accuracy is tol.
  ! Parameters: Maximum allowed number of iterations, and machine floating-point precision.
  integer iter
  real(rkind) a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
  a=x1
  b=x2
  fa=func(a)
  fb=func(b)
  if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))  then
     print*, 'root must be bracketed for zbrent' 
     stop
  end if
  c=b
  fc=fb
  do iter=1,ITMAX
     if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.)) then
        c=a !Rename a, b, c and adjust bounding interval d. 
        fc=fa
        d=b-a
        e=d
     endif
     if(abs(fc).lt.abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     endif
     tol1=2.*EPS*abs(b)+0.5*tol ! Convergence check
     xm=.5*(c-b)
     if(abs(xm).le.tol1 .or. fb.eq.0.)then
        zbrent=b
        return 
     endif
     if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
        s=fb/fa ! Attempt inverse quadratic interpolation. 
        if(a.eq.c) then
           p=2.*xm*s
           q=1.-s
        else
           q=fa/fc
           r=fb/fc 
           p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
           q=(q-1.)*(r-1.)*(s-1.)
        endif
        if(p.gt.0.) q=-q ! Check whether in bounds. 
        p=abs(p)
        if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
           e=d ! Accept interpolation.
           d=p/q
        else
           d=xm ! Interpolation failed, use bisection.
           e=d 
        endif
     else ! Bounds decreasing too slowly, use bisection.
        d=xm
        e=d 
     endif
     a=b ! Move last best guess to a.
     fa=fb
     if(abs(d) .gt. tol1) then ! Evaluate new trial root.
        b=b+d 
     else
        b=b+sign(tol1,xm)
     endif
     fb=func(b)
  enddo

  print*, 'zbrent exceeding maximum iterations'

  zbrent=b
  return
end function zbrent

subroutine linspace(x0,x1,n,xs)
  use prec_const
  implicit none
  real(rkind), intent(in) :: x0,x1
  integer, intent(in) :: n
  real(rkind), dimension(n), intent(out) :: xs

  integer :: i
  real(rkind) :: dx

  dx = (x1-x0)/real(n-1,rkind)

  do i=1,n
     xs(i) = x0 + dx*real(i-1,rkind)
  end do

end subroutine linspace

subroutine interpolation(xs,ys,n,x,res)
  ! constant spacing in xs is assumed 
  use prec_const
  implicit none
  integer, intent(in) :: n
  real(rkind), dimension(n) :: xs,ys
  real(rkind), intent(out) :: res
  real(rkind) :: x,dx
  real(rkind) :: a,b,c,num,den
  real(rkind) :: x0,x1,x2,f0,f1,f2
  integer :: i

  dx = xs(2)-xs(1)

  i = min(int((x-xs(1))/dx)+1,n-2)

  x0 = xs(i)
  x1 = xs(i+1)
  x2 = xs(i+2)
  f0 = ys(i)
  f1 = ys(i+1)
  f2 = ys(i+2)

  den = x0**2*(x1-x2) + x1**2*(x2-x0) + x2**2*(x0-x1)
  a = (f0*x1 - f0*x2 - f1*x0 + f1*x2 + f2*x0 - f2*x1)/den
  b = (- f0*x1**2  + f0*x2**2  + f1*x0**2  - f1*x2**2  - f2*x0**2  + f2*x1**2)/den
  num = -f0*x1*x2*(x1*(x0-x2)+x2*(-x0+x1))+f1*x0**2*x2*(x0-x2)-f2*x0**2*x1*(x0-x1)
  den = -x1*(x0-x1)*(x0**2-x2**2)+x2*(x0-x2)*(x0**2-x1**2)
  c = (    2._rkind*dx**2*f0 + 3._rkind*dx*f0*x0 &
       & -4._rkind*dx*f1*x0 + dx*f2*x0 + f0*x0**2 &
       & -2._rkind*f1*x0**2  + f2*x0**2)/(2._rkind*dx**2)

  res = a*x**2+b*x+c
  
end subroutine interpolation

subroutine TestHeaviside(case)
  use prec_const
  implicit none
  integer :: case
  real(rkind) :: res
  integer :: i


  ! print*, 'Case', case
  if (case.eq.1) then
     do i=2,10
        call HeavisideQuad(fxy1,fxy1,i,res)

        ! print*, 'Check case', case
        ! print*, 'number of points for theboundary', i
        ! print*, 'res =', res
        ! print*, 'res_th =', 0.025718570293422206
        print*, 'rel diff =', case,i,abs((res - 0.025718570293422206)/res)
        print*, ' '
     end do

  else if (case.eq.2) then
     do i=2,10
        call HeavisideQuad(fxy2,fxy2,i,res)

        ! print*, 'Check case', case
        ! print*, 'number of points for theboundary', i
        ! print*, 'res =', res
        ! print*, 'res_th =', 0.025718570293422206
        print*, 'rel diff =', case,i,abs((res - 0.025718570293422206)/res)
        print*, ' '
     end do

  else if (case.eq.3) then
     do i=2,10
        call HeavisideQuad(fxy3,fxy3,i,res)

        ! print*, 'Check case', case
        ! print*, 'number of points for theboundary', i
        ! print*, 'res =', res
        ! print*, 'res_th =', 0.025718570293422206
        print*, 'rel diff =', case,i,abs((res - 0.025718570293422206)/res)
        print*, ' '
     end do
  else if (case.eq.4) then
     do i=2,10
        call HeavisideQuad(fxy4,fxy4,i,res)

        ! print*, 'Check case', case
        ! print*, 'number of points for theboundary', i
        ! print*, 'res =', res
        ! print*, 'res_th =', 0.025718570293422206
        print*, 'rel diff =', case,i,abs((res - 0.025718570293422206)/res)
        print*, ' '
     end do
  else if (case.eq.5) then
     do i=2,10
        call HeavisideQuad(fxy5,fxy5,i,res)

        ! print*, 'Check case', case
        ! print*, 'number of points for theboundary', i
        ! print*, 'res =', res
        ! print*, 'res_th =', 0.643397004039
        print*, 'rel diff =', case,i,abs((res - 0.643397004039)/res)
        print*, ' '
     end do
  else if (case.eq.6) then
     do i=2,10
        call HeavisideQuad(fxy6,fxy6,i,res)

        ! print*, 'Check case', case
        ! print*, 'number of points for theboundary', i
        ! print*, 'res =', res
        ! print*, 'res_th =', 0.0267303356816
        print*, 'rel diff =', case,i,abs((res - 0.0267303356816)/res)
        print*, ' '
     end do
  else if (case.eq.7) then
     do i=2,10
        call HeavisideQuad(fxy7,fxy7,i,res)

        ! print*, 'Check case', case
        ! print*, 'number of points for theboundary', i
        ! print*, 'res =', res
        ! print*, 'res_th =', 0.0267303356816
        print*, 'rel diff =',case,i,abs((res - 0.0267303356816)/res)
        print*, ' '
     end do
  else if (case.eq.8) then
     do i=2,10
        call HeavisideQuad(fxy8,fxy8,i,res)

        ! print*, 'Check case', case
        ! print*, 'number of points for theboundary', i
        ! print*, 'res =', res
        ! print*, 'res_th =', 0.643397002353
        print*, 'rel diff =', case,i,abs((res - 0.643397002353)/res)
        print*, ' '
     end do
  else if (case.eq.9) then
     do i=2,10
        call HeavisideQuad(fxy9,fxy9,i,res)

        ! print*, 'Check case', case
        ! print*, 'number of points for theboundary', i
        ! print*, 'res =', res
        ! print*, 'res_th =', 0.7923848250883889
        print*, 'rel diff =', case,i,abs((res - 0.7923848250883889)/res)
        print*, ' '
     end do
  else if (case.eq.10) then
     do i=2,10
        call HeavisideQuad(fxy10,fxy10,i,res)

        ! print*, 'Check case', case
        ! print*, 'number of points for theboundary', i
        ! print*, 'res =', res
        ! print*, 'res_th =', 0.7923848250883889
        print*, 'rel diff =', case,i,abs((res - 0.7923848250883889)/res)
        print*, ' '
     end do
  else if (case.eq.11) then
     do i=2,10
        call HeavisideQuad(fxy11,fxy11,i,res)

        ! print*, 'Check case', case
        ! print*, 'number of points for theboundary', i
        ! print*, 'res =', res
        ! print*, 'res_th =', 0.7923848250883889
        print*, 'rel diff =', case,i,abs((res - 0.7923848250883889)/res)
        print*, ' '
     end do
  else if (case.eq.12) then
     do i=2,10
        call HeavisideQuad(fxy12,fxy12,i,res)

        ! print*, 'Check case', case
        ! print*, 'number of points for theboundary', i
        ! print*, 'res =', res
        ! print*, 'res_th =', 0.7923848250883889
        print*, 'rel diff =', case,i,abs((res - 0.7923848250883889)/res)
        print*, ' '
     end do
  end if
contains
  
  function fxy1(x,y)
    real(rkind) :: x,y,fxy1

    fxy1 = -(1.6_rkind-(x)**2-1.5_rkind*(y-1)**2)

  end function fxy1

  function fxy2(x,y)
    real(rkind) :: x,y,fxy2

    fxy2 = -(1.6_rkind-(x)**2-1.5_rkind*y**2)

  end function fxy2

  function fxy3(x,y)
    real(rkind) :: x,y,fxy3

    fxy3 = -(1.6_rkind-(x-1)**2-1.5_rkind*y**2)

  end function fxy3

  function fxy4(x,y)
    real(rkind) :: x,y,fxy4

    fxy4 = -(1.6_rkind-(x-1)**2-1.5_rkind*(y-1)**2)

  end function fxy4

  function fxy5(x,y)
    real(rkind) :: x,y,fxy5

    fxy5 = (1.2_rkind-(x-0.5_rkind)**2-1.5_rkind*(y)**2)

  end function fxy5

  function fxy6(x,y)
    real(rkind) :: x,y,fxy6

    fxy6 = -(1.2_rkind-(x-0.5_rkind)**2-1.5_rkind*(y)**2)

  end function fxy6

  function fxy7(x,y)
    real(rkind) :: x,y,fxy7

    fxy7 = -(1.2_rkind-1.5_rkind*(x)**2-(y-0.5_rkind)**2)

  end function fxy7

  function fxy8(x,y)
    real(rkind) :: x,y,fxy8

    fxy8 = (1.2_rkind-1.5_rkind*(x)**2-(y-0.5_rkind)**2)

  end function fxy8

  function fxy9(x,y)
    real(rkind) :: x,y,fxy9

    fxy9 = (1.6_rkind - (x)**2-1.5_rkind*(y-1._rkind)**2)

  end function fxy9

  function fxy10(x,y)
    real(rkind) :: x,y,fxy10

    fxy10 = (1.6_rkind - (x)**2-1.5_rkind*(y)**2)

  end function fxy10

  function fxy11(x,y)
    real(rkind) :: x,y,fxy11

    fxy11 = (1.6_rkind - (x-1._rkind)**2-1.5_rkind*(y)**2)

  end function fxy11

  function fxy12(x,y)
    real(rkind) :: x,y,fxy12

    fxy12 = (1.6_rkind - (x-1._rkind)**2-1.5_rkind*(y-1._rkind)**2)

  end function fxy12

end subroutine TestHeaviside

