subroutine Save
  use globals
  implicit none
  integer :: mp
  real(rkind) :: psimaxval

  mp = 101

  call openbin(mp,'FRC.bin','unformatted','write','big_endian')
  open(mp+1, file = 'FRC.dat', FORM = 'formatted', action = 'write')

  call matwrtI1(mp,'nz',nz)
  call matwrtI1(mp,'nr',nr)
  call matwrtI1(mp,'nws',nws)
  call matwrtI1(mp,'nkws',nkws)
  call matwrtM1(mp,'length',length)
  call matwrtM1(mp,'psiedge',psiedge)

  call matwrtM(mp,'Z',nz,1,Z)
  call matwrtM(mp,'R',nr,1,R)
  call matwrtM(mp,'RR',nz,nr,RR)
  call matwrtM(mp,'ZZ',nz,nr,ZZ)

  call matwrtM(mp,'Coeffs',1,256,CoeffMat)
  call matwrtM(mp,'CM',1,256*(nr-1),ContributionMat)
  call matwrtM(mp,'B',nws,1,B_BC)

  call matwrtI1(mp,'nzguess',nzguess)
  call matwrtI1(mp,'nrguess',nrguess)
  call matwrtM(mp,'PsiGuess',nzguess,nrguess,PsiGuess)
  call matwrtM(mp,'PsiCur',nws,1,PsiCur)

  call matwrtI1(mp,'ilast',ilast)
  call matwrtM(mp,'AllPsis',nws,ilast,AllPsis(:,1:ilast))
  PsiFinal = AllPsis(:,ilast)

  call matwrtM1(mp,'CSol',LambdaFinal)
  call matwrtM(mp,'PsiSol',nws,1,PsiFinal)

  call SaveMesh(mp,PsiFinal)

  close(mp)
  close(mp+1)

  call PsiMaximum(AllPsis(:,ilast),psimaxval,.true.)

end subroutine Save

subroutine SaveMesh(mp,Psi)
  use globals
  implicit none
  integer :: mp
  real(rkind), dimension(nws), intent(in) :: Psi
  real(rkind), dimension(npsi) :: SMesh,PsiMesh,PressureMesh
  real(rkind), dimension(npsi,ntheta+1) :: ZMesh,RMesh,JacobMesh

  call Mesh(Psi,npsi,ntheta,ZMesh,RMesh,JacobMesh,Smesh,PsiMesh,PressureMesh)
  
  call matwrtM(mp,'s',npsi,1,Smesh)
  call matwrtM(mp,'psi',npsi,1,PsiMesh)
  call matwrtM(mp,'pressure',npsi,1,PressureMesh)
  call matwrtM1(mp,'psimax',Psimesh(1))
  call matwrtM(mp,'ZMesh',npsi,ntheta+1,ZMesh)
  call matwrtM(mp,'RMesh',npsi,ntheta+1,RMesh)
  call matwrtM(mp,'JacobMesh',npsi,ntheta+1,JacobMesh)

end subroutine SaveMesh

subroutine ReadGuess
  use globals
  implicit none
  real(rkind), dimension(3) :: sizes_dpr
  ! real(rkind), allocatable, dimension(:,:) :: psig
  
  integer :: mp

  integer :: i

  mp = 101
  
  open(mp, file='guess.dat', status='old', access='stream', form='unformatted')
  read(mp) sizes_dpr
  close(mp)
  nzguess = int(sizes_dpr(1))
  nrguess = int(sizes_dpr(2))
  lguess = sizes_dpr(3)

  dzguess = 1._rkind/real(nzguess-1,rkind)*length
  drguess = 1._rkind/real(nrguess-1,rkind)

  if (lguess.ne.length) then
     print*, 'length must be equal to that of the guess'
     stop
  end if

  ! allocate(psig(nrguess,nzguess))
  allocate(PsiGuess(nzguess,nrguess))

  open(mp, file='guess.bin', status='old', access='stream', form='unformatted')
  read(mp) PsiGuess
  close(mp)

  ! PsiGuess = transpose(psig)

  ! deallocate(psig)

end subroutine ReadGuess

subroutine ReadPprime
  use globals
  implicit none

  integer :: mp
  integer :: i
  real(rkind), external :: ppfun
  mp = 101

  open(unit=mp,access='sequential',form='formatted', file='PPRIME')
  rewind(mp)
  
  read(mp,'(I5)') npprime

  ! allocate(Xpprime(npprime), &
  !      &   Ypprime(npprime), &
  !      &   Bpprime(npprime), &
  !      &   Cpprime(npprime), &
  !      &   Dpprime(npprime))

  ! do i=1,npprime
  !    read(mp,*) Xpprime(i),Ypprime(i)
  ! end do
  
  allocate(Apprime(npprime))

  do i=1,npprime
     read(mp,*) Apprime(i)
  end do

  CLOSE(UNIT=mp)

  ! call SPLINE(npprime,Xpprime,Ypprime,Bpprime,Cpprime,Dpprime)

end subroutine ReadPprime

subroutine openbin(mp,filename,oformat,oaction,oconvert)
  !
  implicit none
  integer, intent(in)          :: mp
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: oformat
  character(len=*), intent(in) :: oaction
  character(len=*), intent(in) :: oconvert
  !
  open(mp, file=filename,FORM=oformat,action=oaction,convert=oconvert)
  !      
  return
end subroutine openbin

subroutine matwrtI1(mf,name,intvalue)
  use prec_const
  implicit none

  integer, intent(in)                 :: mf
  integer, intent(in)                 :: intvalue
  character(len=*), intent(in)        :: name
  integer                 :: m,n
  character(1)            :: info
  !
  info='*'
  n=1
  m=1
  write(mf) real(intvalue, 4)
  write(mf+1,'(A, 2X, A, 1X, I10, 1X, I10)') name, info, n, m
  !
  return
end subroutine matwrtI1

subroutine matwrtI(mf,name,n,m,intvalue)
  use prec_const
  implicit none

  integer, intent(in)                 :: mf, n, m
  integer, dimension(n*m), intent(in) :: intvalue
  character(len=*), intent(in)        :: name
  character(1)            :: info
  !
  info='*'
  write(mf) real(intvalue, 4)
  write(mf+1,'(A, 2X, A, 1X, I10, 1X, I10)') name, info, n, m
  !
  return
end subroutine matwrtI

subroutine matwrtM1(mf,name,realvalue)
  use prec_const
  implicit none

  integer, intent(in)                     :: mf
  character(len=*), intent(in)            :: name
  real(rkind), intent(in)                 :: realvalue
  integer                 :: m,n
  character(1)            :: info
  !
  info='*'
  n=1
  m=1
  write(mf) real(realvalue,4)
  write(mf+1,'(A, 2X, A, 1X, I10, 1X, I10)') name, info, n, m
  !
  return
end subroutine matwrtM1

subroutine matwrtM(mf,name,n,m,realvalue)
  use prec_const
  implicit none

  integer, intent(in)                      :: mf, n, m
  character(len=*), intent(in)             :: name
  real(rkind), dimension(n*m), intent(in)  :: realvalue
  character(1)                :: info
  !
  info='*'
  write(mf) real(realvalue,4)
  write(mf+1,'(A, 2X, A, 1X, I10, 1X, I10)') name, info, n, m
  !
  return
end subroutine matwrtM

