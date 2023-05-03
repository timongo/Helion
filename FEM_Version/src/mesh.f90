subroutine PsiMaximum(psi,rmax,psimaxval,equatorial)
  use globals
  use ieee_arithmetic
  implicit none
  real(rkind), dimension(nws), intent(in) :: Psi
  real(rkind), dimension(nws+nkws) :: PsiAll
  real(rkind), intent(out) :: psimaxval
  logical, intent(in) :: equatorial
  real(rkind) :: r0,r1,rmax
  real(rkind) :: psimaxtol
  real(rkind), external :: zbrent
  
  integer :: mp

  r0 = 1.e-6_rkind
  r1 = 1._rkind
  psimaxtol = 1.e-8_rkind

  count = count+1

  if (count.eq.18) then
     mp = 101
     call openbin(mp,'Debug.bin','unformatted','write','big_endian')
     open(mp+1, file = 'Debug.dat', FORM = 'formatted', action = 'write')

     call matwrtM(mp,'psi',nws,1,psi)

     close(mp)
     close(mp+1)
  end if

  if (equatorial) then
     call FillPsiAll(Psi,PsiAll)

     rmax = zbrent(PsiPrimeEquator,r0,r1,psimaxtol)

     if (.not.ieee_is_finite(rmax)) rmax=0._rkind

     psimaxval = PsiEquator(rmax)
     
     ! print*, count,rmax,psimaxval

  end if
  
contains
 
  function PsiPrimeEquator(r) result(psip)
    real(rkind) :: r,psip
    real(rkind) :: rind
    integer :: ir,iz
    real(rkind) :: t
    integer :: i
    real(rkind), dimension(4) :: psipatch
    integer :: ind1_1,ind1_3,ind2_1,ind2_3
    
    rind = r/deltar

    ir = min(int(rind),nr-2)
    t = rind - real(ir,rkind)
    
    ir = ir+1
    iz = (nz+1)/2

    ind1_1 = IndArrayInv(iz,ir  ,1)
    ind1_3 = IndArrayInv(iz,ir  ,3)
    ind2_1 = IndArrayInv(iz,ir+1,1)
    ind2_3 = IndArrayInv(iz,ir+1,3)
    psipatch = (PsiAll(ind1_1)*CoeffMat(1,:,0,0,1) + &
         &      PsiAll(ind1_3)*CoeffMat(1,:,0,0,3) + &
         &      PsiAll(ind2_1)*CoeffMat(1,:,0,1,1) + &
         &      PsiAll(ind2_3)*CoeffMat(1,:,0,1,3))


    psip = psipatch(2)   + &
         & 2._rkind*psipatch(3)*t + &
         & 3._rkind*psipatch(4)*t**2
  
    psip = psip/deltar

  end function PsiPrimeEquator

  function PsiEquator(r) result(psival)
    real(rkind) :: r,psival
    real(rkind) :: rind
    integer :: ir,iz
    real(rkind) :: t
    integer :: i
    real(rkind), dimension(4) :: psipatch
    integer :: ind1_1,ind1_3,ind2_1,ind2_3
    
    rind = r/deltar

    ir = min(int(rind),nr-2)
    t = rind - real(ir,rkind)
    
    ir = ir+1
    iz = (nz+1)/2

    ind1_1 = IndArrayInv(iz,ir  ,1)
    ind1_3 = IndArrayInv(iz,ir  ,3)
    ind2_1 = IndArrayInv(iz,ir+1,1)
    ind2_3 = IndArrayInv(iz,ir+1,3)
    psipatch = (PsiAll(ind1_1)*CoeffMat(1,:,0,0,1) + &
         &      PsiAll(ind1_3)*CoeffMat(1,:,0,0,3) + &
         &      PsiAll(ind2_1)*CoeffMat(1,:,0,1,1) + &
         &      PsiAll(ind2_3)*CoeffMat(1,:,0,1,3))


    psival = psipatch(1)      + &
         &   psipatch(2)*t    + &
         &   psipatch(3)*t**2 + &
         &   psipatch(4)*t**3
  
  end function PsiEquator

end subroutine PsiMaximum

subroutine Mesh(Psi,np,nth,Zmesh,RMesh,JacMesh,Smesh,PsiMesh,PressureMesh)
  use globals
  implicit none
  real(rkind), dimension(nws), intent(in) :: Psi
  integer, intent(in) :: np,nth
  real(rkind), dimension(np), intent(out) :: SMesh,PsiMesh,PressureMesh
  real(rkind), dimension(np,nth+1), intent(out) :: ZMesh,RMesh,JacMesh

  real(rkind), dimension(nws+nkws) :: PsiAll
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
  
  call FillPsiAll(Psi,PsiAll)

  equatorial = .false.
  if (mod(nz,2).eq.1) equatorial = .true.

  if (equatorial) then
     call PsiMaximum(Psi,rmax,psimaxval,equatorial)

     theta1 = atan(1._rkind-rmax,hlength)
     theta2 = pi - theta1
     theta4 = twopi - atan(rmax,hlength)
     theta3 = pi+twopi-theta4

     call linspace(0._rkind,1._rkind,np,SMesh)
     call linspace(0._rkind,twopi,nth+1,thetas)

     PsiMesh(:) = psimaxval*(1._rkind-SMesh(:)**2)
     PsiMesh(npsi) = 1.e-16_rkind

     RMesh(1,:) = rmax
     ZMesh(1,:) = 0._rkind

     ! Z,R mesh
     do it=1,nth
        theta = thetas(it)
        do ip=2,np
           psival = PsiMesh(ip)
           call ZRPsi(PsiAll,rmax,psival,theta,zval,rval,ztol)
           ZMesh(ip,it) = zval
           RMesh(ip,it) = rval
        end do
     end do

     ! Pressure mesh
     PressureMesh = 0._rkind
     do ip=np-1,1,-1
        psi0 = PsiMesh(ip+1)
        psi1 = PsiMesh(ip)
        call quad(funC,psi0,psi1,deltap)
        PressureMesh(ip) = PressureMesh(ip+1)+deltap
     end do

     ! Jacobian
     JacMesh = 0._rkind
     do it=1,nth
        theta = thetas(it)
        do ip=1,np
           if (ip.eq.1) then
              psival = psimaxval*(1._rkind-(1._rkind*SMesh(2))**2)
              call ZRPsi(PsiAll,rmax,psival,theta,zval,rval,ztol)
           else
              zval = ZMesh(ip,it)
              rval = RMesh(ip,it)
           end if
           call Jacobian(PsiAll,rmax,zval,rval,theta,jval)
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
  
  function funC(psi)
    implicit none
    real(rkind) :: psi,funC

    funC = ppfun(psi)*LambdaFinal

  end function funC

end subroutine Mesh

subroutine Jacobian(PsiAll,rmax,zv,rv,theta,jv)
  use globals
  implicit none
  real(rkind), dimension(nws+nkws), intent(in) :: PsiAll
  real(rkind), intent(in) :: rmax,zv,rv,theta
  real(rkind), intent(out) :: jv
  
  real(rkind) :: dzpsi,drpsi
  real(rkind) :: sint,cost
  
  call EvalPsi(PsiAll,zv,rv,dzpsi,2)
  call EvalPsi(PsiAll,zv,rv,drpsi,3)

  sint = sin(theta)
  cost = cos(theta)

  jv = rv*sqrt(zv**2+(rv-rmax)**2)/abs(drpsi*sint+dzpsi*cost)

end subroutine Jacobian

subroutine ZRPsi(PsiAll,rmax,psiv,thv,zv,rv,ztol)
  use globals
  implicit none
  real(rkind), dimension(nws+nkws), intent(in) :: PsiAll
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

    call EvalPsi(PsiAll,zv,rv,psit1,1)
    psit1 = psit1 - psiv

  end function psit1

  function psit2(t)
    real(rkind) :: t,psit2
    real(rkind) :: zv,rv

    rv = rmax + (1._rkind-rmax)*t
    zv = (rv-rmax)*cotant

    call EvalPsi(PsiAll,zv,rv,psit2,1)
    psit2 = psit2 - psiv

  end function psit2
  
  function psit3(t)
    real(rkind) :: t,psit3
    real(rkind) :: zv,rv

    zv = -hlength*t
    rv = rmax + zv*tant

    call EvalPsi(PsiAll,zv,rv,psit3,1)
    psit3 = psit3 - psiv

  end function psit3

  function psit4(t)
    real(rkind) :: t,psit4
    real(rkind) :: zv,rv

    rv = rmax*(1._rkind - t)
    zv = -rmax*t*cotant

    call EvalPsi(PsiAll,zv,rv,psit4,1)
    psit4 = psit4 - psiv

  end function psit4
end subroutine ZRPsi

subroutine EvalPsi(nz,nr,nws,nkws,IndArrayInv,PsiAll,zv,rv,psival,kcase)
  ! kcase=1 --> psi
  ! kcase=2 --> dpsi/dz
  ! kcase=3 --> dpsi/dr
  ! kcase=4 --> d2psi/dzdr
  use globals, only :: hlength
  use prec_const
  implicit none
  integer, intent(in) :: nz,nr,nws,nkws
  real(rkind), dimension(nz,nr,3), intent(in) :: IndArrayInv
  real(rkind), dimension(nws+nkws), intent(in) :: PsiAll
  real(rkind), intent(in) :: zv,rv
  integer, intent(in) :: kcase
  real(rkind), intent(out) :: psival
  real(rkind) :: zind,rind
  real(rkind) :: t,u

  integer :: iz,ir
  integer :: i,ii

  integer, dimension(4) :: i00,i10,i01,i11
  real(rkind), dimension(16) :: f,alpha
  real(rkind), dimension(4,4) :: co

  zind = (zv+hlength)/deltaz
  rind = rv/deltar

  iz = min(int(zind),nz-2)
  ir = min(int(rind),nr-2)

  t = zind - real(iz,rkind)
  u = rind - real(ir,rkind)

  iz = iz+1
  ir = ir+1

  i00 = IndArrayInv(iz  ,ir  ,:)
  i10 = IndArrayInv(iz+1,ir  ,:)
  i01 = IndArrayInv(iz  ,ir+1,:)
  i11 = IndArrayInv(iz+1,ir+1,:)

  do i=1,4
     
     ii = 1+(i-1)*4

     f(ii  ) = PsiAll(i00(i))*dRdZ(i)
     f(ii+1) = PsiAll(i10(i))*dRdZ(i)
     f(ii+2) = PsiAll(i01(i))*dRdZ(i)
     f(ii+3) = PsiAll(i11(i))*dRdZ(i)
  
  end do

  alpha = matmul(Hmat,f)

  co(:,1) = alpha(1:4)
  co(:,2) = alpha(5:8)
  co(:,3) = alpha(9:12)
  co(:,4) = alpha(13:16)  

  select case (kcase)
     case(1)
        call EvalPsiPatch(co,t,u,psival)
     case(2)
        call EvaldZPsiPatch(co,t,u,psival)
     case(3)
        call EvaldRPsiPatch(co,t,u,psival)
     case(4)
        call EvaldZdRPsiPatch(co,t,u,psival)
  end select

end subroutine EvalPsi

