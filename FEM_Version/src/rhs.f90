subroutine RightHandSide(Psi,fun,rrhs)
  ! In the rhs construction, boundary refers to the plasma vacuum boundary,
  ! not to the domain boundary
  use globals
  implicit none
  real(rkind), dimension(nws), intent(in) :: Psi
  real(rkind), dimension(nws), intent(out) :: rrhs
  real(rkind), external :: fun
  real(rkind) :: res
  integer :: ipatch,npatches               
  integer :: ind,jnd
  integer, dimension(4,2) :: patches
  integer :: iZo,iRo,io,jo,ko
  integer :: iZ,iR,k
  real(rkind), dimension(4,4) :: omega
  real(rkind) :: psi00,psi01,psi10,psi11
  real(rkind),dimension(nws+nkws) :: PsiAll
  logical :: boundary

  rrhs = 0._rkind

  call FillPsiAll(Psi,Psiall)

  do ind=1,nws
     call FindPatches(ind,npatches,patches)
     iZo = IndArray(ind,1)
     iRo = IndArray(ind,2)
     ko = IndArray(ind,3)

     do ipatch=1,npatches
        iZ = patches(ipatch,1)
        iR = patches(ipatch,2)
        io = iZo - iZ
        jo = iRo - iR
        omega = CoeffMat(:,:,io,jo,ko)

        jnd = IndArrayInv(iZ  ,iR  ,1)
        psi00 = PsiAll(jnd)
        jnd = IndArrayInv(iZ+1,iR  ,1)
        psi10 = PsiAll(jnd)
        jnd = IndArrayInv(iZ  ,iR+1,1)
        psi01 = PsiAll(jnd)
        jnd = IndArrayInv(iZ+1,iR+1,1)
        psi11 = PsiAll(jnd)

        ! Special case iR=1
        if (iR.eq.1) then
           if (sign(1._rkind,psi01)*sign(1._rkind,psi11).lt.0._rkind) then
              boundary = .true.
              call RHSPatch(omega,PsiAll,iZ,iR,fun,boundary,res)
              rrhs(ind) = rrhs(ind) + res
           else if (sign(1._rkind,psi01).ge.0._rkind .and. psi11.ge.0._rkind) then
              boundary = .false.
              call RHSPatch(omega,PsiAll,iZ,iR,fun,boundary,res)
              rrhs(ind) = rrhs(ind) + res
           end if
        else
           if (sign(1._rkind,psi00).lt.0._rkind) then
              if (   sign(1._rkind,psi10).gt.0._rkind .or. &
                   & sign(1._rkind,psi01).gt.0._rkind .or. &
                   & sign(1._rkind,psi11).gt.0._rkind) then
                 ! Patch contains a plasma/vacuum boundary
                 boundary = .true.
                 call RHSPatch(omega,PsiAll,iZ,iR,fun,boundary,res)
                 rrhs(ind) = rrhs(ind) + res
              end if
           else if(sign(1._rkind,psi00).ge.0._rkind) then
              if (   sign(1._rkind,psi10).gt.0._rkind .and. &
                   & sign(1._rkind,psi01).gt.0._rkind .and. &
                   & sign(1._rkind,psi11).gt.0._rkind) then
                 ! Patch does not contain a plasma/vacuum boundary
                 boundary = .false.
                 call RHSPatch(omega,PsiAll,iZ,iR,fun,boundary,res)
                 rrhs(ind) = rrhs(ind) + res
              else              
                 ! Patch contains a plasma/vacuum boundary
                 boundary = .true.
                 call RHSPatch(omega,PsiAll,iZ,iR,fun,boundary,res)
                 rrhs(ind) = rrhs(ind) + res
              end if
           end if
        end if
     end do
  end do

end subroutine RightHandSide

subroutine TotalCurrent(Psi,fun,Itot)
  ! In the rhs construction, boundary refers to the plasma vacuum boundary,
  ! not to the domain boundary
  use globals
  implicit none
  real(rkind), dimension(nws), intent(in) :: Psi
  real(rkind), intent(out) :: Itot
  real(rkind), external :: fun
  real(rkind) :: res
  integer :: ind
  integer :: iZ,iR,k
  real(rkind) :: psi00,psi01,psi10,psi11
  real(rkind),dimension(nws+nkws) :: PsiAll
  logical :: boundary

  print*, 'Total current computation'

  Itot = 0._rkind

  call FillPsiAll(Psi,PsiAll)

  do iZ=1,nz-1
     ! Special case iR=1
     iR = 1
     ind = IndArrayInv(iZ  ,iR+1,1)
     psi01 = PsiAll(ind)
     ind = IndArrayInv(iZ+1,iR+1,1)
     psi11 = PsiAll(ind)

     if (sign(1._rkind,psi01)*sign(1._rkind,psi11).lt.0._rkind) then
        boundary = .true.
        call CurrentPatch(PsiAll,iZ,iR,fun,boundary,res)
        Itot = Itot + res
     else if (sign(1._rkind,psi01).ge.0._rkind .and. psi11.ge.0._rkind) then
        boundary = .false.
        call CurrentPatch(PsiAll,iZ,iR,fun,boundary,res)
        Itot = Itot + res
     end if
     
     ! General case iR>0
     do iR=2,nr-1

        ind = IndArrayInv(iZ  ,iR  ,1)
        psi00 = PsiAll(ind)
        ind = IndArrayInv(iZ+1,iR  ,1)
        psi10 = PsiAll(ind)
        ind = IndArrayInv(iZ  ,iR+1,1)
        psi01 = PsiAll(ind)
        ind = IndArrayInv(iZ+1,iR+1,1)
        psi11 = PsiAll(ind)

        if (sign(1._rkind,psi00).lt.0._rkind) then
           if (   sign(1._rkind,psi10).gt.0._rkind .or. &
                & sign(1._rkind,psi01).gt.0._rkind .or. &
                & sign(1._rkind,psi11).gt.0._rkind) then
              ! Patch contains a plasma/vacuum boundary
              boundary = .true.
              call CurrentPatch(PsiAll,iZ,iR,fun,boundary,res)
              Itot = Itot + res
           end if
        else if(sign(1._rkind,psi00).ge.0._rkind) then
           if (   sign(1._rkind,psi10).gt.0._rkind .and. &
                & sign(1._rkind,psi01).gt.0._rkind .and. &
                & sign(1._rkind,psi11).gt.0._rkind) then
              ! Patch does not contain a plasma/vacuum boundary
              boundary = .false.
              call CurrentPatch(PsiAll,iZ,iR,fun,boundary,res)
              Itot = Itot + res
           else              
              ! Patch contains a plasma/vacuum boundary
              boundary = .true.
              call CurrentPatch(PsiAll,iZ,iR,fun,boundary,res)
              Itot = Itot + res
           end if
        end if
     end do
  end do

end subroutine TotalCurrent

subroutine EvalPsiPatch(psi,z,r,psival)
  use prec_const
  implicit none
  real(rkind), intent(in), dimension(4,4) :: psi
  real(rkind), intent(in) :: z,r
  real(rkind), intent(out) :: psival
  integer :: m,n

  psival = 0._rkind

  do m=0,3
     do n=0,3
        psival = psival + psi(m+1,n+1)*z**m*r**n
     end do
  end do

end subroutine EvalPsiPatch

subroutine EvaldZPsiPatch(psi,z,r,psival)
  use prec_const
  use globals, only : deltaz
  implicit none
  real(rkind), intent(in), dimension(4,4) :: psi
  real(rkind), intent(in) :: z,r
  real(rkind), intent(out) :: psival
  integer :: m,n

  psival = 0._rkind

  do m=1,3
     do n=0,3
        psival = psival + psi(m+1,n+1)*real(m,rkind)*z**(m-1)*r**n
     end do
  end do

  psival = psival/deltaz

end subroutine EvaldZPsiPatch

subroutine EvaldRPsiPatch(psi,z,r,psival)
  use prec_const
  use globals, only : deltar
  implicit none
  real(rkind), intent(in), dimension(4,4) :: psi
  real(rkind), intent(in) :: z,r
  real(rkind), intent(out) :: psival
  integer :: m,n

  psival = 0._rkind

  do m=0,3
     do n=1,3
        psival = psival + psi(m+1,n+1)*real(n,rkind)*z**m*r**(n-1)
     end do
  end do

  psival = psival/deltar

end subroutine EvaldRPsiPatch

subroutine EvaldZdRPsiPatch(psi,z,r,psival)
  use prec_const
  use globals, only : deltar,deltaz
  implicit none
  real(rkind), intent(in), dimension(4,4) :: psi
  real(rkind), intent(in) :: z,r
  real(rkind), intent(out) :: psival
  integer :: m,n

  psival = 0._rkind

  do m=1,3
     do n=1,3
        psival = psival + psi(m+1,n+1)*real(n*m,rkind)*z**(m-1)*r**(n-1)
     end do
  end do

  psival = psival/(deltar*deltaz)

end subroutine EvaldZdRPsiPatch

subroutine Jphi(psi,fun,z,r,jval)
  use prec_const
  implicit none
  real(rkind), intent(in), dimension(4,4) :: psi
  real(rkind), external :: fun
  real(rkind), intent(in) :: z,r
  real(rkind), intent(out) :: jval
  real(rkind) :: psival

  call EvalPsiPatch(psi,z,r,psival)

  jval = fun(psi)

end subroutine Jphi

subroutine PsiPatch(PsiAll,iZ,iR,psi)
  use globals
  implicit none
  real(rkind), dimension(nws+nkws), intent(in) :: PsiAll
  integer, intent(in) :: iZ,iR
  real(rkind), dimension(4,4), intent(out) :: psi
  real(rkind), dimension(16) :: f,alpha
  integer :: i,j,k,l
  integer :: ind

  f = 0._rkind
  do i=0,1
     do j=0,1
        do k=1,4
           l = 1+i+2*j+4*(k-1)
           ind = IndArrayInv(iZ+i,iR+j,k)
           f(l) = dRdZ(k)*PsiAll(ind)
        end do
     end do
  end do

  alpha = matmul(Hmat,f)

  psi(:,1) = alpha(1:4)
  psi(:,2) = alpha(5:8)
  psi(:,3) = alpha(9:12)
  psi(:,4) = alpha(13:16)

end subroutine PsiPatch

subroutine CurrentPatch(PsiAll,iZ,iR,fun,boundary,res)
  use globals
  implicit none
  real(rkind), dimension(nws+nkws), intent(in) :: PsiAll
  integer, intent(in) :: iZ,iR
  real(rkind), external :: fun
  logical, intent(in) :: boundary
  real(rkind), intent(out) :: res
  real(rkind), dimension(4,4) :: psi

  call PsiPatch(PsiAll,iZ,iR,psi)

  if (boundary) then
     call HeavisideQuad(CurInt,PsiFun,nboundarypoints,res)
  else
     call Gauss2DQuad(CurInt,gaussorder,res)
  end if

  res = res*deltar**2*deltaz

contains
  
  function PsiFun(z,r) result(psival)
    real(rkind) :: z,r,psival

    call EvalPsiPatch(psi,z,r,psival)

  end function PsiFun

  function CurInt(z,r) result(I)
    real(rkind) :: z,r,I
    integer :: m,n
    real(rkind) :: jval
    
    call Jphi(psi,fun,z,r,jval)

    I = (real(iR-1,rkind)+r)*jval

  end function CurInt

end subroutine CurrentPatch

subroutine RHSPatch(omega,PsiAll,iZ,iR,fun,boundary,res)
  use globals
  implicit none
  real(rkind), dimension(4,4), intent(in) :: omega
  real(rkind), dimension(nws+nkws), intent(in) :: PsiAll
  integer, intent(in) :: iZ,iR
  real(rkind), external :: fun
  logical, intent(in) :: boundary
  real(rkind), intent(out) :: res
  real(rkind), dimension(4,4) :: psi
  
  call PsiPatch(PsiAll,iZ,iR,psi)

  if (boundary) then
     call HeavisideQuad(CurInt,PsiFun,nboundarypoints,res)
  else
     call Gauss2DQuad(CurInt,gaussorder,res)
  end if

  res = res*deltar**2*deltaz

contains
  
  function PsiFun(z,r) result(psival)
    real(rkind) :: z,r,psival

    call EvalPsiPatch(psi,z,r,psival)

  end function PsiFun

  function CurInt(z,r) result(I)
    real(rkind) :: z,r,I
    integer :: m,n
    real(rkind) :: jval
    
    I = 0._rkind
    
    do m=0,3
       do n=0,3
          I = I + omega(m+1,n+1)*z**m*r**n
       end do
    end do

    call Jphi(psi,fun,z,r,jval)

    I = I*(real(iR-1,rkind)+r)*jval

  end function CurInt

end subroutine RHSPatch

subroutine FillPsiAll(Psi,PsiAll)
  use globals
  implicit none
  real(rkind), dimension(nws), intent(in) :: Psi
  real(rkind), dimension(nws+nkws), intent(out) :: PsiAll
  integer :: ind
  integer :: iZ,iR,k
  
  PsiAll = 0._rkind

  PsiAll(:nws) = Psi(:)
  do ind=nws,nws+nkws
     iZ = IndArray(ind,1)
     iR = IndArray(ind,2)
     k = IndArray(ind,3)
     PsiAll(ind) = PsiBC(iZ,iR,k)
  end do

end subroutine FillPsiAll
