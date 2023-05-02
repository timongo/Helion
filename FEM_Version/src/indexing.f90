subroutine Ind_to_iRiZ
  ! 1 -> psi
  ! 2-> dZpsi
  ! 3 -> dRpsi
  ! 4 -> dRdZpsi
  use globals
  implicit none
  integer :: iZ,iR,k,i,j

  IndArray = 0

  ! Unknowns from i=1 to nws
  i=0
  ! Boundary conditions from nws to nkws
  j=nws

  ! Bottom left corner
  do k=1,4
     j=j+1
     IndArray(j,1) = 1
     IndArray(j,2) = 1
     IndArray(j,3) = k
  end do

  ! Bottom edge, iR=0, iZ=1->nZ
  do iZ=2,nZ-1
     do k=1,4
        j=j+1
        IndArray(j,1) = iZ
        IndArray(j,2) = 1
        IndArray(j,3) = k
     end do
  end do

  ! Bottom right corner, iR=0, iZ=nZ+1
  do k=1,4
     j=j+1
     IndArray(j,1) = nZ
     IndArray(j,2) = 1
     IndArray(j,3) = k
  end do


  do iR=2,nr-1
     ! Left edge, dof = dZpsi and dRdZpsi
     j=j+1
     IndArray(j,1) = 1
     IndArray(j,2) = iR
     IndArray(j,3) = 1
     i=i+1
     IndArray(i,1) = 1
     IndArray(i,2) = iR
     IndArray(i,3) = 2
     j=j+1
     IndArray(j,1) = 1
     IndArray(j,2) = iR
     IndArray(j,3) = 3
     i=i+1
     IndArray(i,1) = 1
     IndArray(i,2) = iR
     IndArray(i,3) = 4

     ! Bulk
     do iZ=2,nz-1
        do k=1,4
           i=i+1
           IndArray(i,1) = iZ
           IndArray(i,2) = iR
           IndArray(i,3) = k
        end do
     end do

     ! Right edge, dof = dZpsi and dRdZpsi
     j=j+1
     IndArray(j,1) = nZ
     IndArray(j,2) = iR
     IndArray(j,3) = 1
     i=i+1
     IndArray(i,1) = nZ
     IndArray(i,2) = iR
     IndArray(i,3) = 2
     j=j+1
     IndArray(j,1) = nZ
     IndArray(j,2) = iR
     IndArray(j,3) = 3
     i=i+1
     IndArray(i,1) = nZ
     IndArray(i,2) = iR
     IndArray(i,3) = 4
  end do

  ! Top left corner
  j=j+1
  IndArray(j,1) = 1
  IndArray(j,2) = nR
  IndArray(j,3) = 1
  j=j+1
  IndArray(j,1) = 1
  IndArray(j,2) = nR
  IndArray(j,3) = 2
  j=j+1
  IndArray(j,1) = 1
  IndArray(j,2) = nR
  IndArray(j,3) = 3
  i=i+1
  IndArray(i,1) = 1
  IndArray(i,2) = nR
  IndArray(i,3) = 4

  ! Top edge, iR=nR+1, iZ=1->nZ. Dof = dRpsi and dRdZpsi
  do iZ=2,nz-1
     j=j+1
     IndArray(j,1) = iZ
     IndArray(j,2) = nR
     IndArray(j,3) = 1
     j=j+1
     IndArray(j,1) = iZ
     IndArray(j,2) = nR
     IndArray(j,3) = 2
     i = i+1
     IndArray(i,1) = iZ
     IndArray(i,2) = nR
     IndArray(i,3) = 3
     i = i+1
     IndArray(i,1) = iZ
     IndArray(i,2) = nR
     IndArray(i,3) = 4
  end do

  ! Top right corner, iR=nR+1, iZ=nZ+1
  j=j+1
  IndArray(j,1) = nZ
  IndArray(j,2) = nR
  IndArray(j,3) = 1
  j=j+1
  IndArray(j,1) = nZ
  IndArray(j,2) = nR
  IndArray(j,3) = 2
  j=j+1
  IndArray(j,1) = nZ
  IndArray(j,2) = nR
  IndArray(j,3) = 3
  i=i+1
  IndArray(i,1) = nZ
  IndArray(i,2) = nR
  IndArray(i,3) = 4

end subroutine Ind_to_iRiZ

subroutine iRiZ_to_Ind
  use globals
  implicit none
  integer :: i,j,iZ,iR,k

  IndArrayInv = 0

  ! First corner, iR=1,iZ=1
  i=0
  j=nws
  iR=1
  iZ=1

  j=j+1
  IndArrayInv(iZ,iR,1) = j
  j=j+1
  IndArrayInv(iZ,iR,2) = j
  j=j+1
  IndArrayInv(iZ,iR,3) = j
  j=j+1
  IndArrayInv(iZ,iR,4) = j

  ! Bottom edge, iR=1, iZ=2->nZ-1
  do iZ=2,nz-1
     j=j+1
     IndArrayInv(iZ,iR,1) = j
     j=j+1
     IndArrayInv(iZ,iR,2) = j
     j=j+1
     IndArrayInv(iZ,iR,3) = j
     j=j+1
     IndArrayInv(iZ,iR,4) = j
  end do
  ! Bottom right corner, iR=1, iZ=nZ
  iR=1
  iZ=nZ
  j=j+1
  IndArrayInv(iZ,iR,1) = j
  j=j+1
  IndArrayInv(iZ,iR,2) = j
  j=j+1
  IndArrayInv(iZ,iR,3) = j
  j=j+1
  IndArrayInv(iZ,iR,4) = j

  do iR=2,nr-1
     ! Left edge, dof = dZpsi and dRdZpsi
     iZ=1
     j=j+1
     IndArrayInv(iZ,iR,1) = j
     i=i+1
     IndArrayInv(iZ,iR,2) = i
     j=j+1
     IndArrayInv(iZ,iR,3) = j
     i=i+1
     IndArrayInv(iZ,iR,4) = i

     ! Bulk
     do iZ=2,nz-1
        do k=1,4
           i=i+1
           IndArrayInv(iZ,iR,k) = i
        end do
     end do

     ! Right edge, dof = dZpsi and dRdZpsi
     iZ=nz
     j=j+1
     IndArrayInv(iZ,iR,1) = j
     i=i+1
     IndArrayInv(iZ,iR,2) = i
     j=j+1
     IndArrayInv(iZ,iR,3) = j
     i=i+1
     IndArrayInv(iZ,iR,4) = i
  end do

  ! Top left corner
  iR=nR
  iZ=1
  j=j+1
  IndArrayInv(iZ,iR,1) = j
  j=j+1
  IndArrayInv(iZ,iR,2) = j
  j=j+1
  IndArrayInv(iZ,iR,3) = j
  i=i+1
  IndArrayInv(iZ,iR,4) = i


  ! Top edge, iR=nR+1, iZ=1->nZ. Dof = dRpsi and dRdZpsi
  do iZ=2,nz-1
     j=j+1
     IndArrayInv(iZ,iR,1) = j
     j=j+1
     IndArrayInv(iZ,iR,2) = j
     i=i+1
     IndArrayInv(iZ,iR,3) = i
     i=i+1
     IndArrayInv(iZ,iR,4) = i
  end do

  ! Top right corner, iR=nR+1, iZ=nZ+1
  iZ=nZ
  j=j+1
  IndArrayInv(iZ,iR,1) = j
  j=j+1
  IndArrayInv(iZ,iR,2) = j
  j=j+1
  IndArrayInv(iZ,iR,3) = j
  i=i+1
  IndArrayInv(iZ,iR,4) = i

end subroutine iRiZ_to_Ind
