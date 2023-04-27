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
