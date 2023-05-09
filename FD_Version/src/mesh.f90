subroutine PsiMaximum(psim,rmax,psimaxval,equatorial)
  use globals
  implicit none
  real(rkind), dimension(nztot,nrtot), intent(in) :: Psim
  real(rkind), intent(out) :: psimaxval
  logical, intent(in) :: equatorial
  real(rkind) :: rmax
  integer :: iz,irmax
  
  integer :: mp

  if (equatorial) then

     iz = (nztot+1)/2

     irmax = maxloc(psim(iz,:),1)

     rmax = R(irmax)

     psimaxval = psim(iz,irmax)
     
  end if

end subroutine PsiMaximum

subroutine Mesh(Lambda,Psim,np,nth,Zmesh,RMesh,JacMesh,Smesh,PsiMesh,PressureMesh)
  use globals
  implicit none
  real(rkind), intent(in) :: Lambda
  real(rkind), dimension(nztot,nrtot), intent(in) :: Psim
  integer, intent(in) :: np,nth
  real(rkind), dimension(np), intent(out) :: SMesh,PsiMesh,PressureMesh
  real(rkind), dimension(np,nth+1), intent(out) :: ZMesh,RMesh,JacMesh

  real(rkind), dimension(nth+1) :: thetas
  real(rkind) :: rmax,psimaxval,theta
  real(rkind) :: psival
  real(rkind) :: cost,sint,tant,cotant
  real(rkind) :: tval
  logical :: equatorial
  integer :: it,ip
  real(rkind) :: deltap
  real(rkind) :: psi0,psi1

  real(rkind) :: zval,rval,jval
  real(rkind), parameter :: ztol = 1.e-8
  real(rkind), parameter :: ztol_jac = 1.e-10
  real(rkind), external :: ppfun

  print*, 'Meshing'
  
  equatorial = .false.
  if (mod(nztot,2).eq.1) equatorial = .true.

  if (equatorial) then
     call PsiMaximum(Psim,rmax,psimaxval,equatorial)

     theta1 = atan(1._rkind-rmax,hlength)
     theta2 = pi - theta1
     theta4 = twopi - atan(rmax,hlength)
     theta3 = pi+twopi-theta4

     call linspace(0._rkind,1._rkind,np,SMesh)
     call linspace(0._rkind,twopi,nth+1,thetas)

     PsiMesh(:) = psimaxval*(1._rkind-SMesh(:)**2)
     PsiMesh(np) = 1.e-16_rkind

     RMesh(1,:) = rmax
     ZMesh(1,:) = 0._rkind

     ! Z,R mesh
     do it=1,nth
        theta = thetas(it)
        do ip=2,np
           psival = PsiMesh(ip)
           call ZRPsi(psim,rmax,psival,theta,zval,rval,ztol)
           ZMesh(ip,it) = zval
           RMesh(ip,it) = rval
        end do
     end do

     ! Pressure mesh
     PressureMesh = 0._rkind
     do ip=np-1,1,-1
        psi0 = PsiMesh(ip+1)
        psi1 = PsiMesh(ip)
        call quad(funLambda,psi0,psi1,deltap)
        PressureMesh(ip) = PressureMesh(ip+1)+deltap
     end do

     ! Jacobian
     JacMesh = 0._rkind
     do it=1,nth
        theta = thetas(it)
        do ip=1,np
           if (ip.eq.1) then
              psival = psimaxval*(1._rkind-(1._rkind*SMesh(2))**2)
              call ZRPsi(psim,rmax,psival,theta,zval,rval,ztol)
           else
              zval = ZMesh(ip,it)
              rval = RMesh(ip,it)
           end if
           call Jacobian(psim,rmax,zval,rval,theta,jval)
           JacMesh(ip,it) = jval
        end do
     end do
     
     do ip=1,np
        ZMesh(ip,nth+1) = ZMesh(ip,1)
        RMesh(ip,nth+1) = RMesh(ip,1)
        JacMesh(ip,nth+1) = JacMesh(ip,1)
     end do

  end if

contains
  
  function funLambda(psiv)
    implicit none
    real(rkind) :: psiv,funLambda

    funLambda = ppfun(psiv)*Lambda

  end function funLambda

end subroutine Mesh

subroutine Jacobian(psim,rmax,zv,rv,theta,jv)
  use globals
  implicit none
  real(rkind), dimension(nztot,nrtot), intent(in) :: psim
  real(rkind), intent(in) :: rmax,zv,rv,theta
  real(rkind), intent(out) :: jv
  
  real(rkind) :: dzpsi,drpsi
  real(rkind) :: sint,cost
  
  call EvalPsi(psim,zv,rv,dzpsi,2)
  call EvalPsi(psim,zv,rv,drpsi,3)

  sint = sin(theta)
  cost = cos(theta)

  jv = rv*sqrt(zv**2+(rv-rmax)**2)/abs(drpsi*sint+dzpsi*cost)

end subroutine Jacobian

subroutine ZRPsi(psim,rmax,psiv,thv,zv,rv,ztol)
  use globals
  implicit none
  real(rkind), dimension(nztot,nrtot), intent(in) :: psim
  real(rkind), intent(in) :: rmax
  real(rkind), intent(in) :: psiv,thv
  real(rkind), intent(in) :: ztol
  real(rkind), intent(out) :: zv,rv

  real(rkind) :: sint,cost
  real(rkind) :: tant,cotant
  real(rkind) :: tval

  real(rkind), external :: zbrent

  sint = sin(thv)
  cost = cos(thv)

  if (thv.lt.theta1 .or. thv.ge.theta4) then
     tant = sint/cost
     tval = zbrent(psit1,0._rkind,1._rkind,ztol)
     zv = hlength*tval
     rv = rmax + hlength*tval*tant
  elseif (thv.ge.theta1 .and. thv.lt.theta2) then
     cotant = cost/sint
     tval = zbrent(psit2,0._rkind,1._rkind,ztol)
     zv = (1._rkind-rmax)*tval*cotant
     rv = rmax + (1._rkind-rmax)*tval
  elseif (thv.ge.theta2 .and. thv.lt.theta3) then
     tant = sint/cost
     tval = zbrent(psit3,0._rkind,1._rkind,ztol)
     zv = -hlength*tval
     rv = rmax - hlength*tval*tant
  elseif (thv.ge.theta3 .and. thv.lt.theta4) then
     cotant = cost/sint
     tval = zbrent(psit4,0._rkind,1._rkind,ztol)
     zv = -rmax*tval*cotant
     rv = rmax*(1._rkind - tval)
  end if
  
contains
  function psit1(t)
    real(rkind) :: t,psit1
    real(rkind) :: zv,rv

    zv = hlength*t
    rv = rmax + zv*tant

    call EvalPsi(psim,zv,rv,psit1,1)
    psit1 = psit1 - psiv

  end function psit1

  function psit2(t)
    real(rkind) :: t,psit2
    real(rkind) :: zv,rv

    rv = rmax + (1._rkind-rmax)*t
    zv = (rv-rmax)*cotant

    call EvalPsi(psim,zv,rv,psit2,1)
    psit2 = psit2 - psiv

  end function psit2
  
  function psit3(t)
    real(rkind) :: t,psit3
    real(rkind) :: zv,rv

    zv = -hlength*t
    rv = rmax + zv*tant

    call EvalPsi(psim,zv,rv,psit3,1)
    psit3 = psit3 - psiv

  end function psit3

  function psit4(t)
    real(rkind) :: t,psit4
    real(rkind) :: zv,rv

    rv = rmax*(1._rkind - t)
    zv = -rmax*t*cotant

    call EvalPsi(psim,zv,rv,psit4,1)
    psit4 = psit4 - psiv

  end function psit4
end subroutine ZRPsi

subroutine EvalPsi(psim,zv,rv,psival,kcase)
  ! kcase=1 --> psi
  ! kcase=2 --> dpsi/dz
  ! kcase=3 --> dpsi/dr
  ! kcase=4 --> d2psi/dzdr
  use globals
  implicit none
  real(rkind), dimension(nztot,nrtot), intent(in) :: psim
  real(rkind), intent(in) :: zv,rv
  integer, intent(in) :: kcase
  real(rkind), intent(out) :: psival
  real(rkind) :: zind,rind
  real(rkind) :: t,u

  integer :: iz,ir
  integer :: i,ii

  real(rkind) :: psi00,psi10,psi01,psi11

  zind = (zv+hlength)/deltaz
  rind = rv/deltar

  iz = min(int(zind),nztot-2)
  ir = min(int(rind),nrtot-2)

  t = zind - real(iz,rkind)
  u = rind - real(ir,rkind)

  iz = iz+1
  ir = ir+1

  psi00 = psim(iz,ir)
  psi10 = psim(iz+1,ir)
  psi01 = psim(iz,ir+1)
  psi11 = psim(iz+1,ir+1)

  select case (kcase)
     case(1)
        call EvalPsiPatch(psi00,psi10,psi01,psi11,t,u,psival)
     case(2)
        call EvaldZPsiPatch(psi00,psi10,psi01,psi11,t,u,deltaz,psival)
     case(3)
        call EvaldRPsiPatch(psi00,psi10,psi01,psi11,t,u,deltar,psival)
     case(4)
        call EvaldZdRPsiPatch(psi00,psi10,psi01,psi11,t,u,deltaz,deltar,psival)
  end select

end subroutine EvalPsi

subroutine EvalPsiPatch(psi00,psi10,psi01,psi11,z,r,psival)
  ! Compute the value of psi within a patch at the point z,r (in [0,1]x[0,1]) of the patch
  use prec_const
  implicit none
  real(rkind), intent(in) :: psi00,psi10,psi01,psi11
  real(rkind), intent(in) :: z,r
  real(rkind), intent(out) :: psival

  psival = (psi00*(1._rkind-z)*(1._rkind-r) + &
       &    psi10*z*(1._rkind-r) +  &
       &    psi01*(1._rkind-z)*r + &
       &    psi11*z*r)

end subroutine EvalPsiPatch

subroutine EvaldZPsiPatch(psi00,psi10,psi01,psi11,z,r,dz,psival)
  ! Compute the value of dpsidZ within a patch at the point z,r (in [0,1]x[0,1]) of the patch
  use prec_const
  implicit none
  real(rkind), intent(in) :: psi00,psi10,psi01,psi11
  real(rkind), intent(in) :: z,r
  real(rkind), intent(in) :: dz
  real(rkind), intent(out) :: psival

  psival = (psi10 - psi00)*(1._rkind-r) + (psi11 - psi01)*r
  psival = psival/dz

end subroutine EvaldZPsiPatch

subroutine EvaldRPsiPatch(psi00,psi10,psi01,psi11,z,r,dr,psival)
  ! Compute the value of dpsidR within a patch at the point z,r (in [0,1]x[0,1]) of the patch
  use prec_const
  implicit none
  real(rkind), intent(in) :: psi00,psi10,psi01,psi11
  real(rkind), intent(in) :: z,r
  real(rkind), intent(in) :: dr
  real(rkind), intent(out) :: psival

  psival = (psi01 - psi00)*(1._rkind-z) + (psi11 - psi10)*z
  psival = psival/dr

end subroutine EvaldRPsiPatch

subroutine EvaldZdRPsiPatch(psi00,psi10,psi01,psi11,z,r,dz,dr,psival)
  ! Compute the value of d2psidZdR within a patch at the point z,r (in [0,1]x[0,1]) of the patch
  use prec_const
  implicit none
  real(rkind), intent(in) :: psi00,psi10,psi01,psi11
  real(rkind), intent(in) :: z,r
  real(rkind), intent(in) :: dz,dr
  real(rkind), intent(out) :: psival

  psival = (psi11 + psi00 - psi10 - psi01)/(dr*dz)

end subroutine EvaldZdRPsiPatch

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
