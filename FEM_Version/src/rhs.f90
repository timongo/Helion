subroutine RightHandSide(inds,Psi,fun,rrhs)
  ! In the rhs construction, boundary refers to the plasma vacuum boundary,
  ! not to the domain boundary
  ! For each line of the rhs vector, there is one value for omega
  ! The integral int omega R p'(psi) on the whole domain is then comuted
  use prec_const
  use sizes_indexing
  implicit none
  type(indices), intent(inout) :: inds

  real(rkind), dimension(inds%nws), intent(in) :: Psi
  real(rkind), dimension(inds%nws), intent(out) :: rrhs
  real(rkind), external :: fun
  real(rkind) :: res
  integer :: ipatch,npatches               
  integer :: ind,jnd
  integer, dimension(4,2) :: patches
  integer :: iZo,iRo,io,jo,ko
  integer :: iZ,iR,k
  real(rkind), dimension(4,4) :: omega
  real(rkind) :: psi00,psi01,psi10,psi11
  real(rkind),dimension(inds%ntot) :: PsiAll
  logical :: boundary

  rrhs = 0._rkind

  ! Psi has size nws, PsiAll has size nws+nkws and contains BC information
  call FillPsiAll(inds,Psi,Psiall)

  do ind=1,inds%nws
     call FindPatches(inds,ind,npatches,patches)
     iZo = inds%IndArray(ind,1)
     iRo = inds%IndArray(ind,2)
     ko = inds%IndArray(ind,3)

     do ipatch=1,npatches
        ! Sum the contributions of each patch around the node where omega is nonzero
        iZ = patches(ipatch,1)
        iR = patches(ipatch,2)
        io = iZo - iZ
        jo = iRo - iR
        omega = inds%CoeffMat(:,:,io,jo,ko)

        ! For each patch, the values of psi at the four corner must be retrieved
        ! If they are all negative, the patch lives fully in the vacuum region, there is
        ! no current, therefore the contribution vanishes
        ! If the are all positive, the patch lives fully in the plasma region and the 
        ! integral is performed using gauss quadrature
        ! Otherwise, it means the patch is across the plasma/vacuum boundary, and a special
        ! treatment is required to handle the 'heaviside function' character of the integrand
        jnd = inds%IndArrayInv(iZ  ,iR  ,1)
        psi00 = PsiAll(jnd)
        jnd = inds%IndArrayInv(iZ+1,iR  ,1)
        psi10 = PsiAll(jnd)
        jnd = inds%IndArrayInv(iZ  ,iR+1,1)
        psi01 = PsiAll(jnd)
        jnd = inds%IndArrayInv(iZ+1,iR+1,1)
        psi11 = PsiAll(jnd)

        ! Special case iR=1
        ! required because when iR=1, psi00 and psi10 are always zero
        if (iR.eq.1) then
           if (sign(1._rkind,psi01)*sign(1._rkind,psi11).lt.0._rkind) then
              ! Patch contains a plasma/vacuum boundary
              boundary = .true.
              call RHSPatch(inds,omega,PsiAll,iZ,iR,fun,boundary,res)
              rrhs(ind) = rrhs(ind) + res
           else if (sign(1._rkind,psi01).ge.0._rkind .and. psi11.ge.0._rkind) then
              ! Patch is in plasma and does not contain a plasma/vacuum boundary
              boundary = .false.
              call RHSPatch(inds,omega,PsiAll,iZ,iR,fun,boundary,res)
              rrhs(ind) = rrhs(ind) + res
           end if
        else
           if (sign(1._rkind,psi00).lt.0._rkind) then
              if (   sign(1._rkind,psi10).gt.0._rkind .or. &
                   & sign(1._rkind,psi01).gt.0._rkind .or. &
                   & sign(1._rkind,psi11).gt.0._rkind) then
                 ! Patch contains a plasma/vacuum boundary
                 boundary = .true.
                 call RHSPatch(inds,omega,PsiAll,iZ,iR,fun,boundary,res)
                 rrhs(ind) = rrhs(ind) + res
              end if
           else if(sign(1._rkind,psi00).ge.0._rkind) then
              if (   sign(1._rkind,psi10).gt.0._rkind .and. &
                   & sign(1._rkind,psi01).gt.0._rkind .and. &
                   & sign(1._rkind,psi11).gt.0._rkind) then
                 ! Patch does not contain a plasma/vacuum boundary
                 boundary = .false.
                 call RHSPatch(inds,omega,PsiAll,iZ,iR,fun,boundary,res)
                 rrhs(ind) = rrhs(ind) + res
              else              
                 ! Patch contains a plasma/vacuum boundary
                 boundary = .true.
                 call RHSPatch(inds,omega,PsiAll,iZ,iR,fun,boundary,res)
                 rrhs(ind) = rrhs(ind) + res
              end if
           end if
        end if
     end do
  end do

end subroutine RightHandSide

subroutine TotalCurrent(inds,Psi,fun,Itot)
  ! In the rhs construction, boundary refers to the plasma vacuum boundary,
  ! not to the domain boundary
  use prec_const
  use sizes_indexing
  implicit none
  type(indices), intent(inout) :: inds
  real(rkind), dimension(inds%nws), intent(in) :: Psi
  real(rkind), intent(out) :: Itot
  real(rkind), external :: fun
  real(rkind) :: res
  integer :: ind
  integer :: iZ,iR,k
  real(rkind) :: psi00,psi01,psi10,psi11
  real(rkind),dimension(inds%ntot) :: PsiAll
  logical :: boundary

  print*, 'Total current computation'

  Itot = 0._rkind

  ! Psi has size nws, PsiAll has size nws+nkws and contains BC information
  call FillPsiAll(inds,Psi,PsiAll)

  do iZ=1,inds%nZ-1
     ! Special case iR=1
     ! required because when iR=1, psi00 and psi10 are always zero
     iR = 1
     ind = inds%IndArrayInv(iZ  ,iR+1,1)
     psi01 = PsiAll(ind)
     ind = inds%IndArrayInv(iZ+1,iR+1,1)
     psi11 = PsiAll(ind)

     if (sign(1._rkind,psi01)*sign(1._rkind,psi11).lt.0._rkind) then
        ! Patch contains a plasma/vacuum boundary
        boundary = .true.
        call CurrentPatch(inds,PsiAll,iZ,iR,fun,boundary,res)
        Itot = Itot + res
     else if (sign(1._rkind,psi01).ge.0._rkind .and. psi11.ge.0._rkind) then
        ! Patch is in plasma and does not contain a plasma/vacuum boundary
        boundary = .false.
        call CurrentPatch(inds,PsiAll,iZ,iR,fun,boundary,res)
        Itot = Itot + res
     end if
     
     ! General case iR>0
     do iR=2,inds%nR-1

        ! For each patch, the values of psi at the four corner must be retrieved
        ! If they are all negative, the patch lives fully in the vacuum region, there is
        ! no current, therefore the contribution vanishes
        ! If the are all positive, the patch lives fully in the plasma region and the 
        ! integral is performed using gauss quadrature
        ! Otherwise, it means the patch is across the plasma/vacuum boundary, and a special
        ! treatment is required to handle the 'heaviside function' character of the integrand
        ind = inds%IndArrayInv(iZ  ,iR  ,1)
        psi00 = PsiAll(ind)
        ind = inds%IndArrayInv(iZ+1,iR  ,1)
        psi10 = PsiAll(ind)
        ind = inds%IndArrayInv(iZ  ,iR+1,1)
        psi01 = PsiAll(ind)
        ind = inds%IndArrayInv(iZ+1,iR+1,1)
        psi11 = PsiAll(ind)

        if (sign(1._rkind,psi00).lt.0._rkind) then
           if (   sign(1._rkind,psi10).gt.0._rkind .or. &
                & sign(1._rkind,psi01).gt.0._rkind .or. &
                & sign(1._rkind,psi11).gt.0._rkind) then
              ! Patch contains a plasma/vacuum boundary
              boundary = .true.
              call CurrentPatch(inds,PsiAll,iZ,iR,fun,boundary,res)
              Itot = Itot + res
           end if
        else if(sign(1._rkind,psi00).ge.0._rkind) then
           if (   sign(1._rkind,psi10).gt.0._rkind .and. &
                & sign(1._rkind,psi01).gt.0._rkind .and. &
                & sign(1._rkind,psi11).gt.0._rkind) then
              ! Patch does not contain a plasma/vacuum boundary
              boundary = .false.
              call CurrentPatch(inds,PsiAll,iZ,iR,fun,boundary,res)
              Itot = Itot + res
           else              
              ! Patch contains a plasma/vacuum boundary
              boundary = .true.
              call CurrentPatch(inds,PsiAll,iZ,iR,fun,boundary,res)
              Itot = Itot + res
           end if
        end if
     end do
  end do

end subroutine TotalCurrent

subroutine EvalPsiPatch(psi,z,r,psival)
  ! Compute the value of psi within a patch at the point z,r (in [0,1]x[0,1]) of the patch
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

subroutine EvaldZPsiPatch(psi,z,r,dz,psival)
  ! Compute the value of dpsidZ within a patch at the point z,r (in [0,1]x[0,1]) of the patch
  use prec_const
  implicit none
  real(rkind), intent(in), dimension(4,4) :: psi
  real(rkind), intent(in) :: z,r
  real(rkind), intent(in) :: dz
  real(rkind), intent(out) :: psival
  integer :: m,n

  psival = 0._rkind

  do m=1,3
     do n=0,3
        psival = psival + psi(m+1,n+1)*real(m,rkind)*z**(m-1)*r**n
     end do
  end do

  psival = psival/dz

end subroutine EvaldZPsiPatch

subroutine EvaldRPsiPatch(psi,z,r,dr,psival)
  ! Compute the value of dpsidR within a patch at the point z,r (in [0,1]x[0,1]) of the patch
  use prec_const
  implicit none
  real(rkind), intent(in), dimension(4,4) :: psi
  real(rkind), intent(in) :: z,r
  real(rkind), intent(in) :: dr
  real(rkind), intent(out) :: psival
  integer :: m,n

  psival = 0._rkind

  do m=0,3
     do n=1,3
        psival = psival + psi(m+1,n+1)*real(n,rkind)*z**m*r**(n-1)
     end do
  end do

  psival = psival/dr

end subroutine EvaldRPsiPatch

subroutine EvaldZdRPsiPatch(psi,z,r,dz,dr,psival)
  ! Compute the value of d2psidZdR within a patch at the point z,r (in [0,1]x[0,1]) of the patch
  use prec_const
  implicit none
  real(rkind), intent(in), dimension(4,4) :: psi
  real(rkind), intent(in) :: z,r
  real(rkind), intent(in) :: dz,dr
  real(rkind), intent(out) :: psival
  integer :: m,n

  psival = 0._rkind

  do m=1,3
     do n=1,3
        psival = psival + psi(m+1,n+1)*real(n*m,rkind)*z**(m-1)*r**(n-1)
     end do
  end do

  psival = psival/(dr*dz)

end subroutine EvaldZdRPsiPatch

subroutine JphiPatch(psi,fun,z,r,jval)
  ! Computes pprime(psi(z,r)) with z,r in [0,1]x[0,1] within the considered patch
  use prec_const
  implicit none
  real(rkind), intent(in), dimension(4,4) :: psi
  real(rkind), external :: fun
  real(rkind), intent(in) :: z,r
  real(rkind), intent(out) :: jval
  real(rkind) :: psival

  call EvalPsiPatch(psi,z,r,psival)

  jval = fun(psi)

end subroutine JphiPatch

subroutine PsiPatch(inds,PsiAll,iZ,iR,psi)
  ! extracts the psi(i,j) array from the PsiAll array, using the Hermite matrix
  use prec_const
  use sizes_indexing
  use globals, only : Hmat
  implicit none
  type(indices), intent(inout) :: inds
  real(rkind), dimension(inds%ntot), intent(in) :: PsiAll
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
           ind = inds%IndArrayInv(iZ+i,iR+j,k)
           ! The multiplication by dRdZ is necessary to transform derivatives with respect to
           ! Z and R into derivatives with respect to z,r where z,r are 
           ! the local coordinates in the patch, varying in [0,1]x[0,1]
           f(l) = inds%dRdZ(k)*PsiAll(ind)
        end do
     end do
  end do

  alpha = matmul(Hmat,f)

  psi(:,1) = alpha(1:4)
  psi(:,2) = alpha(5:8)
  psi(:,3) = alpha(9:12)
  psi(:,4) = alpha(13:16)

end subroutine PsiPatch

subroutine CurrentPatch(inds,PsiAll,iZ,iR,fun,boundary,res)
  use prec_const
  use sizes_indexing
  use globals, only : nboundarypoints,gaussorder
  implicit none
  type(indices), intent(inout) :: inds
  real(rkind), dimension(inds%ntot), intent(in) :: PsiAll
  integer, intent(in) :: iZ,iR
  real(rkind), external :: fun
  logical, intent(in) :: boundary
  real(rkind), intent(out) :: res
  real(rkind), dimension(4,4) :: psi

  call PsiPatch(inds,PsiAll,iZ,iR,psi)

  if (boundary) then
     ! Special integration package required if across boundary
     call HeavisideQuad(CurInt,PsiFun,nboundarypoints,res)
  else
     ! Standard 2D gauss quadrature otherwise
     call Gauss2DQuad(CurInt,gaussorder,res)
  end if

  res = res*inds%deltar**2*inds%deltaz

contains
  
  function PsiFun(z,r) result(psival)
    ! Psi(z,r) where z,r are in [0,1]x[0,1] within a patch
    real(rkind) :: z,r,psival

    call EvalPsiPatch(psi,z,r,psival)

  end function PsiFun

  function CurInt(z,r) result(I)
    ! R/dr*Jphi(z,r) where z,r are in [0,1]x[0,1] within a patch
    ! R/dr = (iR-1 + r)
    real(rkind) :: z,r,I
    integer :: m,n
    real(rkind) :: jval
    
    call JphiPatch(psi,fun,z,r,jval)

    I = (real(iR-1,rkind)+r)*jval

  end function CurInt

end subroutine CurrentPatch

subroutine RHSPatch(inds,omega,PsiAll,iZ,iR,fun,boundary,res)
  use prec_const
  use sizes_indexing
  use globals, only : nboundarypoints,gaussorder
  implicit none
  type(indices), intent(inout) :: inds
  real(rkind), dimension(4,4), intent(in) :: omega
  real(rkind), dimension(inds%ntot), intent(in) :: PsiAll
  integer, intent(in) :: iZ,iR
  real(rkind), external :: fun
  logical, intent(in) :: boundary
  real(rkind), intent(out) :: res
  real(rkind), dimension(4,4) :: psi
  
  call PsiPatch(inds,PsiAll,iZ,iR,psi)

  if (boundary) then
     ! Special integration package required if across boundary
     call HeavisideQuad(CurInt,PsiFun,nboundarypoints,res)
  else
     ! Standard 2D gauss quadrature otherwise
     call Gauss2DQuad(CurInt,gaussorder,res)
  end if

  res = res*inds%deltar**2*inds%deltaz

contains
  
  function PsiFun(z,r) result(psival)
    ! Psi(z,r) where z,r are in [0,1]x[0,1] within a patch
    real(rkind) :: z,r,psival

    call EvalPsiPatch(psi,z,r,psival)

  end function PsiFun

  function CurInt(z,r) result(I)
    ! omega(z,r)*R/dr*Jphi(z,r) where z,r are in [0,1]x[0,1] within a patch
    ! R/dr = (iR-1 + r)
    real(rkind) :: z,r,I
    integer :: m,n
    real(rkind) :: jval
    
    I = 0._rkind
    
    do m=0,3
       do n=0,3
          I = I + omega(m+1,n+1)*z**m*r**n
       end do
    end do

    call JphiPatch(psi,fun,z,r,jval)

    I = I*(real(iR-1,rkind)+r)*jval

  end function CurInt

end subroutine RHSPatch

subroutine FillPsiAll(inds,Psi,PsiAll)
  ! Psi has size nws, PsiAll has size nws+nkws and contains BC information
  use prec_const
  use sizes_indexing
  implicit none
  type(indices), intent(inout) :: inds
  real(rkind), dimension(inds%nws), intent(in) :: Psi
  real(rkind), dimension(inds%ntot), intent(out) :: PsiAll
  integer :: ind
  integer :: iZ,iR,k
  
  PsiAll = 0._rkind

  PsiAll(:inds%nws) = Psi(:)
  do ind=inds%nws,inds%ntot
     iZ = inds%IndArray(ind,1)
     iR = inds%IndArray(ind,2)
     k = inds%IndArray(ind,3)
     PsiAll(ind) = inds%PsiBC(iZ,iR,k)
  end do

end subroutine FillPsiAll
