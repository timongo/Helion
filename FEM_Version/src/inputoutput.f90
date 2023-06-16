subroutine Save
  use globals
  implicit none
  integer :: mp
  ! real(rkind) :: psimaxval
  real(rkind), dimension(inds_c%nz,inds_c%nr) :: psi2d_c
  real(rkind), dimension(inds_r%nz,inds_r%nr) :: psi2d_r

  mp = 101

  call openbin(mp,'FRC.bin','unformatted','write','big_endian')
  open(mp+1, file = 'FRC.dat', FORM = 'formatted', action = 'write')

  call matwrtI1(mp,'nz',inds_r%nZ)
  call matwrtI1(mp,'nr',inds_r%nR)
  call matwrtI1(mp,'nws',inds_r%nws)
  call matwrtI1(mp,'nkws',inds_r%nkws)
  call matwrtI1(mp,'nzc',inds_c%nZ)
  call matwrtI1(mp,'nrc',inds_c%nR)
  call matwrtI1(mp,'nwsc',inds_c%nws)
  call matwrtI1(mp,'nkwsc',inds_c%nkws)
  call matwrtM1(mp,'length',inds_r%length)
  call matwrtM1(mp,'psiedge',psiedge)

  call matwrtM(mp,'Z_c',inds_c%nz,1,inds_c%Z)
  call matwrtM(mp,'R_c',inds_c%nr,1,inds_c%R)
  call matwrtM(mp,'RR_c',inds_c%nz,inds_c%nr,inds_c%RR)
  call matwrtM(mp,'ZZ_c',inds_c%nz,inds_c%nr,inds_c%ZZ)
  call ReconstructPsi(inds_c,inds_c%PsiFinal,psi2d_c)
  call matwrtM(mp,'psi2d_c',inds_c%nz,inds_c%nr,psi2d_c)  

  call matwrtM(mp,'Z',inds_r%nz,1,inds_r%Z)
  call matwrtM(mp,'R',inds_r%nr,1,inds_r%R)
  call matwrtM(mp,'RR',inds_r%nz,inds_r%nr,inds_r%RR)
  call matwrtM(mp,'ZZ',inds_r%nz,inds_r%nr,inds_r%ZZ)
  call ReconstructPsi(inds_r,inds_r%PsiFinal,psi2d_r)
  call matwrtM(mp,'psi2d',inds_r%nz,inds_r%nr,psi2d_r)

  ! call matwrtM(mp,'PsiCur',nws,1,PsiCur)

  ! call matwrtI1(mp,'ilast',ilast)
  ! call matwrtM(mp,'AllPsis',nws,ilast,AllPsis(:,1:ilast))
  ! PsiFinal = AllPsis(:,ilast)

  ! if (usepetsc) then
  call matwrtM1(mp,'LambdaSol',LambdaFinal)
  call matwrtM(mp,'PsiSol_c',inds_c%nws,1,inds_c%PsiFinal)
  call matwrtM(mp,'PsiSol',inds_r%nws,1,inds_r%PsiFinal)
  call SaveMesh(mp,inds_r,inds_r%PsiFinal)
  ! else
  !    call matwrtM1(mp,'LambdaSol',LambdaFinal)
  !    call matwrtM(mp,'PsiSol',inds_c%nws,1,inds_c%PsiFinal)
  !    call SaveMesh(mp,inds_c,inds_c%PsiFinal)
  ! end if
  

  close(mp)
  close(mp+1)

  ! write guess
  call openbin(mp,'solution_guess_write_FEM_coarse.bin','unformatted','write','big_endian')
  write(mp) inds_c%nz
  write(mp) inds_c%nr
  write(mp) inds_c%length
  write(mp) inds_c%nws
  write(mp) inds_c%nkws
  write(mp) LambdaFinal
  write(mp) inds_c%PsiFinal
  close(mp)

  ! write guess
  call openbin(mp,'solution_guess_write_FEM_refined.bin','unformatted','write','big_endian')
  write(mp) inds_r%nz
  write(mp) inds_r%nr
  write(mp) inds_r%length
  write(mp) inds_r%nws
  write(mp) inds_r%nkws
  write(mp) LambdaFinal
  write(mp) inds_r%PsiFinal
  close(mp)

end subroutine Save

subroutine ReconstructPsi(inds,Psi,psi2d)
  use prec_const
  use sizes_indexing
  implicit none
  type(indices), intent(in) :: inds
  real(rkind), dimension(inds%nkws), intent(in) :: Psi
  real(rkind), dimension(inds%nz,inds%nr) :: psi2d
  real(rkind), dimension(inds%ntot) :: PsiAll
  integer :: iz,ir
  real(rkind) :: zv,rv
  real(rkind) :: psiv

  psi2d = 0._rkind

  call FillPsiAll(inds,Psi,PsiAll)
  
  do iz=1,inds%nz
     zv = inds%Z(iz)
     do ir=1,inds%nr
        rv = inds%R(ir)
        call EvalPsi(inds,PsiAll,zv,rv,psiv,1)
        psi2d(iz,ir) = psiv
     end do
  end do

end subroutine ReconstructPsi

subroutine SaveMesh(mp,inds,Psi)
  use prec_const
  use sizes_indexing
  use globals, only : npsi,ntheta,LambdaFinal
  implicit none
  type(indices), intent(inout) :: inds
  integer :: mp
  real(rkind), dimension(inds%nws), intent(in) :: Psi
  real(rkind), dimension(npsi) :: SMesh,PsiMesh,PressureMesh
  real(rkind), dimension(npsi,ntheta+1) :: ZMesh,RMesh,JacobMesh

  call Mesh(inds,LambdaFinal,Psi,npsi,ntheta,ZMesh,RMesh,JacobMesh,Smesh,PsiMesh,PressureMesh)
  
  call matwrtM(mp,'S',npsi,1,Smesh)
  call matwrtM(mp,'PsiMesh',npsi,1,PsiMesh)
  call matwrtM(mp,'Pressure',npsi,1,PressureMesh)
  call matwrtM1(mp,'psimax',Psimesh(1))
  call matwrtM(mp,'ZMesh',npsi,ntheta+1,ZMesh)
  call matwrtM(mp,'RMesh',npsi,ntheta+1,RMesh)
  call matwrtM(mp,'JacobMesh',npsi,ntheta+1,JacobMesh)

end subroutine SaveMesh

subroutine ReadGuess
  use globals
  implicit none
  integer :: nzguess,nrguess
  integer :: nwsguess,nkwsguess
  real(rkind) :: lguess
  ! real(rkind), allocatable, dimension(:,:) :: psig
  
  integer :: mp

  integer :: i

  mp = 101

  if (guesstype.eq.1) then
     call openbin(mp,'solution_guess_read_FD.bin','unformatted','read','big_endian')
     read(mp) nzguess
     read(mp) nrguess
     read(mp) lguess
     read(mp) LambdaIni

     allocate(PsiGuess1(nzguess,nrguess))

     read(mp) PsiGuess1

     close(mp)

     inds_g%nZ = nzguess
     inds_g%nR = nrguess
     inds_g%length = lguess
     inds_g%hlength = 0.5_rkind*lguess

  elseif (guesstype.eq.2) then
     
     call openbin(mp,'solution_guess_read_FEM.bin','unformatted','read','big_endian')
     read(mp) nzguess
     read(mp) nrguess
     read(mp) lguess
     read(mp) nwsguess
     read(mp) nkwsguess
     read(mp) LambdaIni

     allocate(inds_g%PsiCur(nwsguess))

     read(mp) inds_g%PsiCur

     close(mp)

     inds_g%nZ = nzguess
     inds_g%nR = nrguess
     inds_g%length = lguess
     inds_g%hlength = 0.5_rkind*lguess     

     call Arrays(inds_g)
     call Ind_to_iRiZ(inds_g)
     call iRiZ_to_Ind(inds_g)
     call PsiBoundaryCondition(inds_g)

  end if
  
  inds_g%deltaz = 1._rkind/real(nzguess-1,rkind)*lguess
  inds_g%deltar = 1._rkind/real(nrguess-1,rkind)
  
  if (lguess.ne.inds_c%length) then
     print*, "length is not equal to that of the guess. Are you sure you know what you're doing?"
     write(*,'(A,E12.4,A,E12.4)') 'Length = ', inds_c%length, 'while guess length is', lguess
  end if
  
  ! if LambdaNL (NL=namelist) is non zero, replace LambdaIni with LambdaNL
  if (LambdaNL.gt.0_rkind) then
     LambdaIni = LambdaNL
  end if

end subroutine ReadGuess

subroutine ReadGuess_Local
  use prec_const
  use globals
  implicit none
  real(rkind) :: lguess
  
  inds_g%nZ = inds_r%nZ
  inds_g%nR = inds_r%nR
  lguess = inds_r%length
  inds_g%length = lguess
  inds_g%hlength = 0.5_rkind*lguess     
  LambdaIni = LambdaFinal
  inds_g%PsiCur = inds_r%PsiFinal
  inds_g%deltaz = 1._rkind/real(inds_g%nZ-1,rkind)*lguess
  inds_g%deltar = 1._rkind/real(inds_g%nR-1,rkind)

  call Arrays(inds_g)
  call Ind_to_iRiZ(inds_g)
  call iRiZ_to_Ind(inds_g)
  call PsiBoundaryCondition(inds_g)

  ! if LambdaNL (NL=namelist) is non zero, replace LambdaIni with LambdaNL
  if (LambdaNL.gt.0_rkind) then
     LambdaIni = LambdaNL
  end if

end subroutine ReadGuess_Local

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

