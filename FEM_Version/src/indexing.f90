subroutine Ind_to_iRiZ(inds)
  ! Gives the answer to the question 'What is iR,iZ and field type if I know ind (from 1 to nws)?'
  ! 1 -> psi
  ! 2-> dZpsi
  ! 3 -> dRpsi
  ! 4 -> dRdZpsi
  use prec_const
  use sizes_indexing
  implicit none
  type(indices), intent(inout) :: inds
  integer :: iZ,iR,k,i,j

  inds%IndArray = 0

  ! Unknowns from i=1 to nws
  i=0
  ! Boundary conditions from nws to nkws
  j=inds%nws

  ! Bottom left corner
  do k=1,4
     j=j+1
     inds%IndArray(j,1) = 1
     inds%IndArray(j,2) = 1
     inds%IndArray(j,3) = k
  end do

  ! Bottom edge, iR=0, iZ=1->nZ
  do iZ=2,inds%nZ-1
     do k=1,4
        j=j+1
        inds%IndArray(j,1) = iZ
        inds%IndArray(j,2) = 1
        inds%IndArray(j,3) = k
     end do
  end do

  ! Bottom right corner, iR=0, iZ=nZ+1
  do k=1,4
     j=j+1
     inds%IndArray(j,1) = inds%nZ
     inds%IndArray(j,2) = 1
     inds%IndArray(j,3) = k
  end do


  do iR=2,inds%nR-1
     ! Left edge, dof = dZpsi and dRdZpsi
     j=j+1
     inds%IndArray(j,1) = 1
     inds%IndArray(j,2) = iR
     inds%IndArray(j,3) = 1
     i=i+1
     inds%IndArray(i,1) = 1
     inds%IndArray(i,2) = iR
     inds%IndArray(i,3) = 2
     j=j+1
     inds%IndArray(j,1) = 1
     inds%IndArray(j,2) = iR
     inds%IndArray(j,3) = 3
     i=i+1
     inds%IndArray(i,1) = 1
     inds%IndArray(i,2) = iR
     inds%IndArray(i,3) = 4

     ! Bulk
     do iZ=2,inds%nZ-1
        do k=1,4
           i=i+1
           inds%IndArray(i,1) = iZ
           inds%IndArray(i,2) = iR
           inds%IndArray(i,3) = k
        end do
     end do

     ! Right edge, dof = dZpsi and dRdZpsi
     j=j+1
     inds%IndArray(j,1) = inds%nZ
     inds%IndArray(j,2) = iR
     inds%IndArray(j,3) = 1
     i=i+1
     inds%IndArray(i,1) = inds%nZ
     inds%IndArray(i,2) = iR
     inds%IndArray(i,3) = 2
     j=j+1
     inds%IndArray(j,1) = inds%nZ
     inds%IndArray(j,2) = iR
     inds%IndArray(j,3) = 3
     i=i+1
     inds%IndArray(i,1) = inds%nZ
     inds%IndArray(i,2) = iR
     inds%IndArray(i,3) = 4
  end do

  ! Top left corner
  j=j+1
  inds%IndArray(j,1) = 1
  inds%IndArray(j,2) = inds%nR
  inds%IndArray(j,3) = 1
  j=j+1
  inds%IndArray(j,1) = 1
  inds%IndArray(j,2) = inds%nR
  inds%IndArray(j,3) = 2
  j=j+1
  inds%IndArray(j,1) = 1
  inds%IndArray(j,2) = inds%nR
  inds%IndArray(j,3) = 3
  i=i+1
  inds%IndArray(i,1) = 1
  inds%IndArray(i,2) = inds%nR
  inds%IndArray(i,3) = 4

  ! Top edge, iR=nR+1, iZ=1->nZ. Dof = dRpsi and dRdZpsi
  do iZ=2,inds%nZ-1
     j=j+1
     inds%IndArray(j,1) = iZ
     inds%IndArray(j,2) = inds%nR
     inds%IndArray(j,3) = 1
     j=j+1
     inds%IndArray(j,1) = iZ
     inds%IndArray(j,2) = inds%nR
     inds%IndArray(j,3) = 2
     i = i+1
     inds%IndArray(i,1) = iZ
     inds%IndArray(i,2) = inds%nR
     inds%IndArray(i,3) = 3
     i = i+1
     inds%IndArray(i,1) = iZ
     inds%IndArray(i,2) = inds%nR
     inds%IndArray(i,3) = 4
  end do

  ! Top right corner, iR=nR+1, iZ=nZ+1
  j=j+1
  inds%IndArray(j,1) = inds%nZ
  inds%IndArray(j,2) = inds%nR
  inds%IndArray(j,3) = 1
  j=j+1
  inds%IndArray(j,1) = inds%nZ
  inds%IndArray(j,2) = inds%nR
  inds%IndArray(j,3) = 2
  j=j+1
  inds%IndArray(j,1) = inds%nZ
  inds%IndArray(j,2) = inds%nR
  inds%IndArray(j,3) = 3
  i=i+1
  inds%IndArray(i,1) = inds%nZ
  inds%IndArray(i,2) = inds%nR
  inds%IndArray(i,3) = 4

end subroutine Ind_to_iRiZ

subroutine iRiZ_to_Ind(inds)
  ! Gives the answer to the question 'What is ind (from 1 to nws) if I know iR,iZ and field type?'
  use prec_const
  use sizes_indexing
  implicit none
  type(indices), intent(inout) :: inds
  integer :: i,j,iZ,iR,k

  inds%IndArrayInv = 0

  ! First corner, iR=1,iZ=1
  i=0
  j=inds%nws
  iR=1
  iZ=1

  j=j+1
  inds%IndArrayInv(iZ,iR,1) = j
  j=j+1
  inds%IndArrayInv(iZ,iR,2) = j
  j=j+1
  inds%IndArrayInv(iZ,iR,3) = j
  j=j+1
  inds%IndArrayInv(iZ,iR,4) = j

  ! Bottom edge, iR=1, iZ=2->nZ-1
  do iZ=2,inds%nZ-1
     j=j+1
     inds%IndArrayInv(iZ,iR,1) = j
     j=j+1
     inds%IndArrayInv(iZ,iR,2) = j
     j=j+1
     inds%IndArrayInv(iZ,iR,3) = j
     j=j+1
     inds%IndArrayInv(iZ,iR,4) = j
  end do
  ! Bottom right corner, iR=1, iZ=nZ
  iR=1
  iZ=inds%nZ
  j=j+1
  inds%IndArrayInv(iZ,iR,1) = j
  j=j+1
  inds%IndArrayInv(iZ,iR,2) = j
  j=j+1
  inds%IndArrayInv(iZ,iR,3) = j
  j=j+1
  inds%IndArrayInv(iZ,iR,4) = j

  do iR=2,inds%nR-1
     ! Left edge, dof = dZpsi and dRdZpsi
     iZ=1
     j=j+1
     inds%IndArrayInv(iZ,iR,1) = j
     i=i+1
     inds%IndArrayInv(iZ,iR,2) = i
     j=j+1
     inds%IndArrayInv(iZ,iR,3) = j
     i=i+1
     inds%IndArrayInv(iZ,iR,4) = i

     ! Bulk
     do iZ=2,inds%nZ-1
        do k=1,4
           i=i+1
           inds%IndArrayInv(iZ,iR,k) = i
        end do
     end do

     ! Right edge, dof = dZpsi and dRdZpsi
     iZ=inds%nZ
     j=j+1
     inds%IndArrayInv(iZ,iR,1) = j
     i=i+1
     inds%IndArrayInv(iZ,iR,2) = i
     j=j+1
     inds%IndArrayInv(iZ,iR,3) = j
     i=i+1
     inds%IndArrayInv(iZ,iR,4) = i
  end do

  ! Top left corner
  iR=inds%nR
  iZ=1
  j=j+1
  inds%IndArrayInv(iZ,iR,1) = j
  j=j+1
  inds%IndArrayInv(iZ,iR,2) = j
  j=j+1
  inds%IndArrayInv(iZ,iR,3) = j
  i=i+1
  inds%IndArrayInv(iZ,iR,4) = i


  ! Top edge, iR=nR+1, iZ=1->nZ. Dof = dRpsi and dRdZpsi
  do iZ=2,inds%nZ-1
     j=j+1
     inds%IndArrayInv(iZ,iR,1) = j
     j=j+1
     inds%IndArrayInv(iZ,iR,2) = j
     i=i+1
     inds%IndArrayInv(iZ,iR,3) = i
     i=i+1
     inds%IndArrayInv(iZ,iR,4) = i
  end do

  ! Top right corner, iR=nR+1, iZ=nZ+1
  iZ=inds%nZ
  j=j+1
  inds%IndArrayInv(iZ,iR,1) = j
  j=j+1
  inds%IndArrayInv(iZ,iR,2) = j
  j=j+1
  inds%IndArrayInv(iZ,iR,3) = j
  i=i+1
  inds%IndArrayInv(iZ,iR,4) = i

end subroutine iRiZ_to_Ind
