subroutine Coeff(i,j,k,co)
  use globals
  implicit none
  integer, intent(in) :: i,j,k
  real(rkind), dimension(4,4), intent(out) :: co
  real(rkind), dimension(16) :: f,alpha
  integer :: l
  
  l = 1+i+2*j+4*(k-1)

  f = 0._rkind
  f(l) = dRdZ(k)

  alpha = matmul(Hmat,f)

  co(:,1) = alpha(1:4)
  co(:,2) = alpha(5:8)
  co(:,3) = alpha(9:12)
  co(:,4) = alpha(13:16)

end subroutine Coeff

subroutine AllCoeffs
  use globals
  implicit none
  integer :: i,j,k
  real(rkind), dimension(4,4) :: co
  
  CoeffMat = 0._rkind

  do i=0,1
     do j=0,1
        do k=1,4
           call Coeff(i,j,k,co)
           CoeffMat(:,:,i,j,k) = co
        end do
     end do
  end do

end subroutine AllCoeffs

subroutine AllRIntegrals
  use globals
  implicit none
  integer :: i,iR
  real(rkind) :: res

  do iR=1,nR-1
     do i=0,6
        call RIntegral(i,iR,res)
        RIntegralArray(i,iR) = res
     end do
  end do
  
end subroutine AllRIntegrals

recursive subroutine RIntegral(i,iR,res)
  use globals
  implicit none
  integer, intent(in) :: i,iR
  real(rkind), intent(out) :: res

  if (iR.eq.1) then
     if (i.eq.0) then
        res = 0._rkind
     else
        res = 1._rkind/(real(i,rkind)+1.e-18_rkind)
     end if
  else
     if (i.eq.0) then
        res = logRp1sR(iR)
     else if (i.gt.0) then
        call RIntegral(i-1,iR,res)
        res = -real(iR-1,rkind)*res + 1._rkind/real(i,rkind)
     else
        res = 0._rkind
     end if
  end if     

end subroutine RIntegral

subroutine SETUPAB
#include <petsc/finclude/petscksp.h>
  use globals
  use petscksp
  implicit none
    
  integer :: ind,jnd
  real(rkind) :: val,valBC
  integer :: iZo,iRo,ko
  integer :: iZp,iRp,kp
  integer :: i,j

  real(rkind), dimension(4,4) :: psi
  real(rkind) :: res

  PetscErrorCode :: ierr
  PetscInt :: pnws
  PetscInt :: dummynz
  PetscInt :: ione
  PetscInt :: II,JJ
  PetscInt, dimension(0:nws-1) :: nnz
  PetscScalar :: pval,one

  ! call AllRIntegrals
  ! call Aij_all

  ione = 1
  dummynz = 1
  nnz = 0
  do ind=1,nws
     iZo = IndArray(ind,1)
     iRo = IndArray(ind,2)
     ko = IndArray(ind,3)
     do i=-1,1
        iZp = iZo+i
        do j=-1,1
           iRp = iRo+j
           if (iRp.ge.1 .and. iRp.le.nr .and. iZp.ge.1 .and. iZp.le.nz) then
              do kp=1,4
                 jnd = IndArrayInv(iZp,iRp,kp)
                 if (jnd.le.nws) then
                    nnz(ind-1) = nnz(ind-1)+1
                 end if
              end do
           end if
        end do
     end do
  end do

  call MatCreate(PETSC_COMM_WORLD,PAMat,ierr)
  pnws = nws
  call MatSetSizes(PAMat,PETSC_DECIDE,PETSC_DECIDE,pnws,pnws,ierr)
  call MatSetType(PAMat,MATSEQAIJ,ierr)
  call MatSetFromOptions(PAMat,ierr)
  call MatSeqAIJSetPreallocation(PAMat,dummynz,nnz,ierr)

  one = 1.d0

  do ind=1,nws
     II = ind-1
     valBC = 0._rkind

     iZo = IndArray(ind,1)
     iRo = IndArray(ind,2)
     ko = IndArray(ind,3)

     do i=-1,1
        iZp = iZo+i
        do j=-1,1
           iRp = iRo+j
           if (iRp.ge.1 .and. iRp.le.nr .and. iZp.ge.1 .and. iZp.le.nz) then
              do kp=1,4
                 jnd = IndArrayInv(iZp,iRp,kp)                 
                 call Aij_bulk(ind,jnd,val)
                 if (jnd.le.nws) then
                    JJ = jnd-1
                    pval = val
                    call MatSetValues(PAMat,ione,II,ione,JJ,pval,INSERT_VALUES,ierr)
                 else
                    valBC = valBC - val*PsiBC(iZp,iRp,kp)
                 end if
              end do
           end if
        end do
     end do
     B_BC(ind) = valBC
  end do

  call MatAssemblyBegin(PAMat,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(PAMat,MAT_FINAL_ASSEMBLY,ierr)

end subroutine SETUPAB

subroutine Aij_bulk(ind,jnd,res)
  use globals
  implicit none
  integer, intent(in) :: ind,jnd
  real(rkind), intent(out) :: res
  integer :: npatches_o,npatches_p,npatches
  integer, dimension(4,2) :: patches_o,patches_p,patches
  integer :: ipatch
  integer :: iZ,iR
  integer :: iZo,iRo,ko
  integer :: iZp,iRp,kp
  integer :: io,jo,ip,jp

  iZo = IndArray(ind,1)
  iRo = IndArray(ind,2)
  ko = IndArray(ind,3)

  iZp = IndArray(jnd,1)
  iRp = IndArray(jnd,2)
  kp = IndArray(jnd,3)

  res = 0._rkind

  if (abs(iZo-iZp).le.1 .and. abs(iRo-iRp).le.1) then
     call FindPatches(ind,npatches_o,patches_o)
     call FindPatches(jnd,npatches_p,patches_p)

     call PatchIntersection(npatches_o,patches_o,npatches_p,patches_p,npatches,patches)

     do ipatch=1,npatches
        iZ = patches(ipatch,1)
        iR = patches(ipatch,2)
        io = iZo-iZ
        ip = iZp-iZ
        jo = iRo-iR
        jp = iRp-iR

        res = res + ContributionMat(io,jo,ko,ip,jp,kp,iR)
     end do
  end if

end subroutine Aij_bulk

subroutine FindPatches(ind,npatches,patches)
  use globals
  implicit none
  integer, intent(in) :: ind
  integer, intent(out) :: npatches
  integer, dimension(4,2), intent(out) :: patches
  integer :: iZ,iR
  integer :: iZi
  integer :: iRj
  integer :: i,j

  iZ = IndArray(ind,1)
  iR = IndArray(ind,2)

  npatches = 0

  patches = 0._rkind

  do i=-1,0
     iZi=iZ+i
     do j=-1,0
        iRj = iR+j
        if (iZi.ge.1 .and. iZi.lt.nz .and. iRj.ge.1 .and. iRj.lt.nr) then
           npatches = npatches+1
           patches(npatches,1) = iZi
           patches(npatches,2) = iRj
        end if
     end do
  end do

end subroutine FindPatches

subroutine PatchIntersection(npatches_o,patches_o,npatches_p,patches_p,npatches,patches)
  use globals
  implicit none
  integer, intent(in) :: npatches_o,npatches_p
  integer, dimension(4,2), intent(in) :: patches_o,patches_p
  integer, intent(out) :: npatches
  integer, dimension(4,2), intent(out) :: patches
  
  integer :: ipatch_o,ipatch_p
  integer :: iZo,iRo
  integer :: iZp,iRp

  npatches = 0
  patches = 0._rkind

  do ipatch_o=1,npatches_o
     iZo = patches_o(ipatch_o,1)
     iRo = patches_o(ipatch_o,2)
     do ipatch_p=1,npatches_p
        iZp = patches_p(ipatch_p,1)
        iRp = patches_p(ipatch_p,2)
        if (iZo.eq.iZp .and. iRo.eq.iRp) then
           npatches = npatches+1
           patches(npatches,:) = patches_o(ipatch_o,:)
        end if
     end do
  end do


end subroutine PatchIntersection

subroutine Aij_all
  use globals
  implicit none

  integer :: iR,io,jo,ko,ip,jp,kp
  real(rkind), dimension(4,4) :: omega,psi
  real(rkind) :: res
  
  ContributionMat = 0._rkind

  print*, 'Matrix Construction'

  do iR=1,nr-1
     do io=0,1
        do jo=0,1
           do ko=1,4
              omega(:,:) = CoeffMat(:,:,io,jo,ko)
              do ip=0,1
                 do jp=0,1
                    do kp=1,4
                       psi(:,:) = CoeffMat(:,:,ip,jp,kp)
                       call ComputeContribution(omega,psi,jo,jp,iR,res)
                       ContributionMat(io,jo,ko,ip,jp,kp,iR) = res
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

end subroutine Aij_all

subroutine ComputeContribution(omega,psi,jo,jp,iR,res)
  use globals
  implicit none
  integer, intent(in) :: jo,jp,iR
  real(rkind), intent(in), dimension(4,4) :: omega,psi
  real(rkind), intent(out) :: res
  real(rkind) :: rval,zval
  integer :: i,j,m,n

  if (iR.eq.1 .and. (jo.eq.0 .or. jp.eq.0)) then
     res = 0._rkind
  else
     res = 0._rkind
     do i=0,3
        do j=0,3
           do m=0,3
              do n=0,3                 
                 zval = omega(i+1,j+1)*psi(m+1,n+1)*real(i*m,rkind)/(real(i+m-1,rkind)+1.e-18)* &
                      & RIntegralArray(j+n,iR)
                 rval = omega(i+1,j+1)*psi(m+1,n+1)*real(j*n,rkind)/real(i+m+1,rkind)* &
                      & RIntegralArray(j+n-2,iR)
                 res = res + zval/deltaz + rval*deltaz/deltar**2
              end do
           end do
        end do
     end do
  end if

end subroutine ComputeContribution
