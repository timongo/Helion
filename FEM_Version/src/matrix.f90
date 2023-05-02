subroutine Coeff(i,j,k,co)
  ! Uses the Hermite matrix to compute aij for a given base case (4 locations on a patch times 4 field types)
  ! i=0,1 -> Z=0,1
  ! j=0,1 -> R=0,1
  ! k = field type
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
  ! Uses the Hermite matrix to compute all possible aij for all 16 base cases (4 locations on a patch times 4 field types)
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
  ! Computing beforehand all possible types of R**i/(R+a) integrals for matrix preparation
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
  ! int_0^1 R**i/(R+iR) integrals for matrix preparation
  ! The analytical result, given recursively, is used
  use globals
  implicit none
  integer, intent(in) :: i,iR
  real(rkind), intent(out) :: res

  if (iR.eq.1) then
     if (i.eq.0) then
        res = 0._rkind
     else
        ! avoid division by 0
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
  ! Assembling the matrix as well as the boundary condition part of the right hand side
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

  ! Counting the number of non zeros in the matrix for fast implementation with PETSc
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

  ! Create PETSc matrix context as PAMat (PETSc A Matrix)
  call MatCreate(PETSC_COMM_WORLD,PAMat,ierr)
  ! PETSc routines like to see PetscInt types, pnws means PETSc nws
  pnws = nws
  call MatSetSizes(PAMat,PETSC_DECIDE,PETSC_DECIDE,pnws,pnws,ierr)
  call MatSetType(PAMat,MATSEQAIJ,ierr)
  call MatSetFromOptions(PAMat,ierr)
  call MatSeqAIJSetPreallocation(PAMat,dummynz,nnz,ierr)
  
  ! Matrix contributions
  one = 1.d0
  do ind=1,nws
     ! Each line of the matrix corresponds to one choice for omega
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
              ! if the supports of omega and psi are disjoints, there is no contribution
              ! otherwise compute A_ind,jnd and put result in val with Aij_bulk(ind,jnd,val)
              do kp=1,4
                 jnd = IndArrayInv(iZp,iRp,kp)     
                 call Aij_bulk(ind,jnd,val)
                 if (jnd.le.nws) then
                    JJ = jnd-1
                    pval = val
                    ! Put entry in PETSc matrix if jnd corresponds to an unknown
                    call MatSetValues(PAMat,ione,II,ione,JJ,pval,INSERT_VALUES,ierr)
                 else
                    ! otherwise Use entry to fill the ind'th line of the right hand side
                    valBC = valBC - val*PsiBC(iZp,iRp,kp)
                 end if
              end do
           end if
        end do
     end do
     ! Set the right hand side (the boundary conditions depending part of the right hand side)
     B_BC(ind) = valBC
  end do

  ! PETSc matrix assembly
  call MatAssemblyBegin(PAMat,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(PAMat,MAT_FINAL_ASSEMBLY,ierr)

end subroutine SETUPAB

subroutine Aij_bulk(ind,jnd,res)
  ! Compute one contribution A_ind,jnd of the matrix. Put the result in res
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

  ! Location and type of ind
  iZo = IndArray(ind,1)
  iRo = IndArray(ind,2)
  ko = IndArray(ind,3)

  ! Location and type of jnd
  iZp = IndArray(jnd,1)
  iRp = IndArray(jnd,2)
  kp = IndArray(jnd,3)

  res = 0._rkind

  ! if the locations are separated by at most one unity, proceed
  if (abs(iZo-iZp).le.1 .and. abs(iRo-iRp).le.1) then
     ! A patch is designated by the values iZ and iR of its bottom left corner
     call FindPatches(ind,npatches_o,patches_o)
     call FindPatches(jnd,npatches_p,patches_p)

     ! Find the patches surrounding ind and jnd that are in common for ind and jnd
     call PatchIntersection(npatches_o,patches_o,npatches_p,patches_p,npatches,patches)

     ! For all such patches, use the contribution matrix to update A_ind,jnd by adding
     ! a contribution to the integral
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
  ! Around a given node, there are four patches, designated by the
  ! Z,R indices of their bottom left corner
  ! This routine returns the intersection of these four patches with the global domain
  ! (that is, if for instance the node is a corner or an edge, some patches are outside the domain)
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
  ! Computing all possible values for the matrix elements
  ! io,jo,ko: location in the patch and field type of omega
  ! ip,jp,kp: location in the patch and field type of psi
  ! io = 0,1: Z=0,1
  ! jo = 0,1: R=0,1
  ! etc
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
  ! Compute one contribution to the matrix
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
                 ! int 1/R domega/dZ dpsi/dZ dZ dR
                 zval = omega(i+1,j+1)*psi(m+1,n+1)*real(i*m,rkind)/(real(i+m-1,rkind)+1.e-18)* &
                      & RIntegralArray(j+n,iR)
                 ! int 1/R domega/dR dpsi/dR dZ dR
                 rval = omega(i+1,j+1)*psi(m+1,n+1)*real(j*n,rkind)/real(i+m+1,rkind)* &
                      & RIntegralArray(j+n-2,iR)
                 ! int grad(omega).grad(psi)/R dZ dR
                 res = res + zval/deltaz + rval*deltaz/deltar**2
              end do
           end do
        end do
     end do
  end if

end subroutine ComputeContribution
