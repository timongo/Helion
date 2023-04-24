program main
  use prec_const
  use globals
  implicit none

  call readnamelist

  call initialization

  call ADI_Solve(psiguess,psi)

  call save

  call deallocate_arrays

end program main

function ppfun(psiv)
  use globals
  real(rkind) :: psiv,ppfun

  real(rkind) :: a1,b1,d1

  ! a1 = 3.204_rkind
  ! b1 = 0.93685_rkind
  ! d1 = 1.278_rkind

  ! a1 = 170.64_rkind
  ! b1 = 0.69919_rkind
  ! d1 = 3.0156_rkind

  ! ppfun = pprime
  ! ppfun = d1/cosh(a1*psiv-b1)**2
  ppfun = 1._rkind

end function ppfun

subroutine ADI_Solve(psig,psi2)
  use prec_const
  use globals
  implicit none

  real(rkind), dimension(nrtot,nztot), intent(inout) :: psig
  real(rkind), dimension(nrtot,nztot), intent(out) :: psi2
  ! real(rkind), intent(in) :: PP

  real(rkind), dimension(nrtot,nztot) :: psi1
  real(rkind), dimension(nrtot,nztot) :: psi2_star
  real(rkind), dimension(nrtot,nztot) :: aux
  real(rkind) :: omega,omegamax
  real(rkind) :: start,finish
  integer :: k
  real(rkind) :: I
  real(rkind) :: I0
  real(rkind) :: C
  real(rkind) :: dtpsi

  real(rkind) :: N
  real(rkind) :: Nstar
  real(rkind) :: norm
  real(rkind) :: ratio
  real(rkind) :: dnrm2
  real(rkind) :: psimaxcur

  omegamax = 1.e-1_rkind

  omega = omegai

  norm = dnrm2(ntot,psig,1)
  ! print*, norm

  C = pprime
  !C = 1._rkind
  call TotalCurrent(psig,C,I0)

  k = 0
  dtpsi = tol
  
  print*, 'current = ', I0
  ! print*, 'psitarget = ', psimax

  do while (k.lt.ntsmax .and. dtpsi.ge.tol)
     call ADI_Step_omp(omega,C,psig,psi1)
     call ADI_Step_omp(omega,C,psi1,psi2)
     call ADI_Step_omp(omega*2._rkind,C,psig,psi2_star)

     call dcopy(ntot,psig,1,aux,1)
     call daxpy(ntot,-1._rkind,psi1,1,aux,1)
     N = dnrm2(ntot,aux,1)

     call dcopy(ntot,psi2,1,aux,1)
     call daxpy(ntot,-1._rkind,psi2_star,1,aux,1)
     Nstar = dnrm2(ntot,aux,1)

     ratio = Nstar/N

     dtpsi = N/(norm*omega)

     write(*,'(I6,4E12.4)') k+1,omega,C,ratio,dtpsi/tol

     if (ratio.le.0.05_rkind) then
        k=k+1
        call dcopy(ntot,psi2,1,psig,1)
        omega = min(omega*4._rkind,omegamax)
        call TotalCurrent(psi2,1._rkind,I)
        C = I0/I
        ! psimaxcur = maxval(psi2)
        ! C = C/(psimaxcur/psimax)
     elseif (ratio.gt.0.05_rkind .and. ratio.le.0.1_rkind) then
        k = k+1
        call dcopy(ntot,psi2,1,psig,1)
        omega = min(omega*2._rkind,omegamax)
        call TotalCurrent(psi2,1._rkind,I)
        C = I0/I        
        ! psimaxcur = maxval(psi2)
        ! C = C/(psimaxcur/psimax)
     elseif (ratio.gt.0.1_rkind .and. ratio.le.0.3_rkind) then
        k = k+1
        call dcopy(ntot,psi2,1,psig,1)
        call TotalCurrent(psi2,1._rkind,I)
        C = I0/I        
        ! psimaxcur = maxval(psi2)
        ! C = C/(psimaxcur/psimax)
     elseif (ratio.gt.0.3_rkind .and. ratio.le.0.4_rkind) then
        k = k+1
        call dcopy(ntot,psi2,1,psig,1)
        omega = omega/2._rkind
        call TotalCurrent(psi2,1._rkind,I)
        C = I0/I        
        ! psimaxcur = maxval(psi2)
        ! C = C/(psimaxcur/psimax)
     elseif (ratio.gt.0.4_rkind .and. ratio.le.0.6_rkind) then
        k = k+1
        call dcopy(ntot,psi2,1,psig,1)
        omega = omega/4._rkind
        call TotalCurrent(psi2,1._rkind,I)
        C = I0/I        
        ! psimaxcur = maxval(psi2)
        ! C = C/(psimaxcur/psimax)
     elseif (ratio.gt.0.6_rkind) then
        k = k+1
        call dcopy(ntot,psi2,1,psig,1)
        omega = omega/16._rkind
     end if
  end do

  print*, 'C = ', C
  call TotalCurrent(psi2,C,I)
  print*, 'current = ', I

end subroutine ADI_Solve

subroutine TotalCurrent(psi,PP,I)
  use prec_const
  use globals, only: nrtot,nztot,R,deltar,deltaz
  implicit none
  real(rkind), dimension(nrtot,nztot), intent(in) :: psi
  real(rkind), intent(in) :: PP
  real(rkind), intent(out) :: I
  real(rkind) :: psiv
  integer :: ir,iz
  real(rkind), external :: ppfun

  I = 0._rkind
  
  !$OMP PARALLEL PRIVATE(ir,iz,psiv) &
  !$OMP          SHARED(PP,R)
  !$OMP DO REDUCTION(+:I)
  do ir=1,nrtot
     do iz=1,nztot
        psiv = psi(ir,iz)
        if (psiv.gt.0._rkind) then
           I = I - PP*ppfun(psiv)*R(ir)*deltar*deltaz
        end if
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

end subroutine TotalCurrent

subroutine ADI_Step_omp(omega,PP,psiold,psinew)
  use globals
  implicit none
  real(rkind),intent(in) :: omega
  real(rkind),intent(in) :: PP
  real(rkind), dimension(nrtot,nztot), intent(in) :: psiold
  real(rkind), dimension(nrtot,nztot), intent(out) :: psinew

  real(rkind), dimension(nr) :: DZpsiold
  real(rkind), dimension(nrtot,nztot) :: psistar
  real(rkind), dimension(nz) :: DRpsistar
  real(rkind), dimension(nrtot,nztot) :: cur

  real(rkind), dimension(nr-1) :: ar
  real(rkind), dimension(nr  ) :: br
  real(rkind), dimension(nr-1) :: cr
  real(rkind), dimension(nr  ) :: dr
  real(rkind), dimension(nr  ) :: xr

  real(rkind), dimension(nz-1) :: az
  real(rkind), dimension(nz  ) :: bz
  real(rkind), dimension(nz-1) :: cz
  real(rkind), dimension(nz  ) :: dz
  real(rkind), dimension(nz  ) :: xz

  integer :: ir,iz
  real(rkind) :: dr1,dz1
  real(rkind) :: dr2,dz2
  real(rkind) :: psiv
  real(rkind), external :: ppfun

  dr1 = 1._rkind/deltar
  dz1 = 1._rkind/deltaz
  dr2 = dr1**2
  dz2 = dz1**2

  ! Psistar

  cur = 0._rkind

  do ir=1,nrtot
     do iz=1,nztot
        psiv = psiold(ir,iz)
        if (psiv.gt.0) then
           cur(ir,iz) = -RR(ir,iz)**2*PP*ppfun(psiv)
        end if
     end do
  end do

  ar = 0._rkind
  br = 1._rkind+omega*2._rkind*dr2
  cr = 0._rkind

  do ir=1,nr-1
     ar(ir) = -omega*(dr2 + 0.5_rkind*dr1/R(ir+2))
     cr(ir) = -omega*(dr2 - 0.5_rkind*dr1/R(ir+1))
  end do

  psistar = 0._rkind  
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP          PRIVATE(iz,DZpsiold,dr,xr)
  !$OMP DO
  do iz=2,nz+1
     DZpsiold = (psiold(2:nr+1,iz+1) + psiold(2:nr+1,iz-1) - 2._rkind*psiold(2:nr+1,iz))*dz2
     dr(:) = psiold(2:nr+1,iz) + omega*DZpsiold - omega*cur(2:nr+1,iz)
     dr(nr) = dr(nr) + omega*psiedge*dr1*(dr1 - 0.5/R(nr+1))
     call TDMA_Solver(ar,br,cr,dr,xr,nr)
     psistar(2:nr+1,iz) = xr(:)
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  psistar(1,:) = 0._rkind
  psistar(nr+2,:) = psiedge
  psistar(:,1) = psiedge*R**2
  psistar(:,nz+2) = psiedge*R**2

  ! Psinew

  cur = 0._rkind

  do ir=1,nrtot
     do iz=1,nztot
        psiv = psistar(ir,iz)
        if (psiv.gt.0) then
           cur(ir,iz) = -RR(ir,iz)**2*PP*ppfun(psiv)
        end if
     end do
  end do

  az = -omega*dz2
  bz = 1._rkind+omega*2._rkind*dz2
  cz = -omega*dz2

  psinew = 0._rkind
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP          PRIVATE(ir,DRpsistar,dz,xz)
  !$OMP DO
  do ir=2,nr+1
     DRpsistar = (psistar(ir+1,2:nz+1) + psistar(ir-1,2:nz+1) - 2._rkind*psistar(ir,2:nz+1))*dr2 - &
          &      0.5_rkind*(psistar(ir+1,2:nz+1) - psistar(ir-1,2:nz+1))*dr1/R(ir)
     dz(:) = psistar(ir,2:nz+1) + omega*DRpsistar - omega*cur(ir,2:nz+1)
     dz(1) = dz(1) + omega*psiedge*dz2*R(ir)**2
     dz(nz) = dz(nz) + omega*psiedge*dz2*R(ir)**2
     call TDMA_Solver(az,bz,cz,dz,xz,nz)
     psinew(ir,2:nz+1) = xz(:)
  end do
  !$OMP END DO
  !$OMP END PARALLEL


  psinew(1,:) = 0._rkind
  psinew(nr+2,:) = psiedge
  psinew(:,1) = psiedge*R**2
  psinew(:,nz+2) = psiedge*R**2

end subroutine ADI_Step_omp

subroutine TDMA_Solver(a,b,c,d,x,n)
  use prec_const
  implicit none
  integer, intent(in) :: n
  real(rkind), dimension(n-1), intent(in) :: a
  real(rkind), dimension(n  ), intent(in) :: b
  real(rkind), dimension(n-1), intent(in) :: c
  real(rkind), dimension(n  ), intent(in) :: d
  real(rkind), dimension(n  ), intent(out) :: x
  
  real(rkind), dimension(n-1) :: ac
  real(rkind), dimension(n  ) :: bc
  real(rkind), dimension(n-1) :: cc
  real(rkind), dimension(n  ) :: dc

  real(rkind) :: mc

  integer :: il,it

  x = 0._rkind
  ac(:) = a(:)
  bc(:) = b(:)
  cc(:) = c(:)
  dc(:) = d(:)
  
  do it=2,n
     mc = ac(it-1)/bc(it-1)
     bc(it) = bc(it) - mc*cc(it-1)
     dc(it) = dc(it) - mc*dc(it-1)     
  end do

  x(:) = bc(:)
  x(n) = dc(n)/bc(n)

  do il=n-1,1,-1
     x(il) = (dc(il) - cc(il)*x(il+1))/bc(il)
  end do

end subroutine TDMA_Solver

subroutine readnamelist
  use prec_const
  use globals
  implicit none
  
  integer :: mp

  nr = 100
  nz = 100
  length = 1._rkind
  psiedge = -0.5_rkind
  
  xs = 0.5_rkind
  b  = 0.5_rkind

  pprime = 19.823216814706463_rkind
  ntsmax = 500
  tol = 1.e-8_rkind
  omegai = 1.e-2_rkind

  psimax = 0.05_rkind

  namelist /frc/ &
       nr,nz,length,psiedge,xs,b,pprime,ntsmax,tol,omegai,psimax

  mp = 101
  open(mp, file ='nlfrc', delim = 'apostrophe', &
       FORM = 'formatted', action = 'read', status='old')
  read(mp,frc)
  close(mp)

  nrtot = nr+2
  nztot = nz+2
  ntot = nrtot*nztot
  
end subroutine readnamelist

subroutine initialization
  use globals
  implicit none
  integer :: i

  allocate(R(nrtot),Z(nztot),RR(nrtot,nztot),ZZ(nrtot,nztot), &
       &   psi(nrtot,nztot),psiguess(nrtot,nztot))

  deltar = 1._rkind/real(nrtot-1,rkind)
  do i=1,nrtot
     R(i) = real(i-1,rkind)/real(nrtot-1,rkind)
     RR(i,:) = real(i-1,rkind)/real(nrtot-1,rkind)
  end do

  deltaz = 1._rkind/real(nztot-1,rkind)*length
  do i=1,nztot
     Z(i) = real(i-1,rkind)/real(nztot-1,rkind)*length - length/2._rkind
     ZZ(:,i) = real(i-1,rkind)/real(nztot-1,rkind)*length - length/2._rkind
  end do

  call read_guess

  ! call Steinhauer(xs,b,psiguess,pprime)

end subroutine initialization

subroutine read_guess
  use globals
  implicit none
  real(rkind), dimension(3) :: sizes_dpr
  real(rkind), allocatable, dimension(:,:) :: psig
  
  integer :: mp

  integer :: i

  mp = 101
  
  open(mp, file='guess.dat', status='old', access='stream', form='unformatted')
  read(mp) sizes_dpr
  close(mp)
  nzguesstot = int(sizes_dpr(1))
  nrguesstot = int(sizes_dpr(2))
  lguess = sizes_dpr(3)
  ! nrguesstot = nrguess+2
  ! nzguesstot = nzguess+2

  print*, nrguesstot,nzguesstot

  drguess = 1._rkind/real(nrguesstot-1,rkind)
  dzguess = 1._rkind/real(nzguesstot-1,rkind)*length

  if (lguess.ne.length) then
     print*, 'length must be equal to that of the guess'
     print*, lguess, length
     ! stop
  end if

  allocate(psig(nzguesstot,nrguesstot))

  open(mp, file='guess.bin', status='old', access='stream', form='unformatted')
  read(mp) psig
  close(mp)

  if (nr.eq.nrguess .and. nz.eq.nzguess) then
     psiguess = transpose(psig)
  else
     allocate(psiguess_read(nrguesstot,nzguesstot))
     psiguess_read = transpose(psig)
     call interpolate_guess
     deallocate(psiguess_read)
  end if

  
  deallocate(psig)

end subroutine read_guess

subroutine interpolate_guess
  use globals
  implicit none

  integer :: ir,iz,irguess,izguess
  real(rkind) :: x,y,t,u

  do ir=2,nr+1
     x = R(ir)/drguess
     irguess = floor(x)+1
     t = x - real(irguess-1,rkind)
     do iz=2,nz+1
        y = (Z(iz)+length/2._rkind)/dzguess
        izguess = floor(y)+1
        u = y - real(izguess-1,rkind)

        psiguess(ir,iz) = (1._rkind-t)*(1._rkind-u)*psiguess_read(irguess  ,izguess  ) + &
             &            (1._rkind-t)*      u     *psiguess_read(irguess  ,izguess+1) + &
             &                  t     *(1._rkind-u)*psiguess_read(irguess+1,izguess  ) + &
             &                  t     *      u     *psiguess_read(irguess+1,izguess+1)
     end do
  end do

  psiguess(1,:) = 0._rkind
  psiguess(nr+2,:) = psiedge
  psiguess(:,1) = psiedge*R**2
  psiguess(:,nz+2) = psiedge*R**2

end subroutine interpolate_guess

subroutine Steinhauer(xs_SH,b_SH,psi_SH,PP_SH)
  use globals
  implicit none
  real(rkind), intent(in) :: xs_SH,b_SH
  real(rkind), dimension(nrtot,nztot), intent(out) :: psi_SH
  real(rkind), intent(out) :: PP_SH
  real(rkind) :: N
  real(rkind) :: eps
  real(rkind) :: D0,D1
  real(rkind), dimension(nrtot,nztot) :: D0term,D1term
  real(rkind) :: Bw,B0,B1

  N = 1._rkind/(1._rkind+(0.8_rkind+xs_SH**2)*(b_SH-1._rkind)**2)
  eps = 1._rkind/b_SH
  D0 = (-8._rkind+eps**4)/(4._rkind*(2._rkind+eps**2))
  D1 = -eps**2*(4._rkind+eps**2)/(4._rkind*(2._rkind+eps**2))
  Bw = 1._rkind/(1._rkind-xs_SH**2)
  B0 = Bw*sqrt(1.5_rkind)*xs_SH*(1+(2.56_rkind-3.2_rkind*xs_SH)/(xs_SH**2*(9.65_rkind-15.2_rkind*xs_SH+b_SH**4)))
  B1 = B0*eps**2*(N-1._rkind)/(D0*(4._rkind+eps**2*N)+D1*(12._rkind+2._rkind*eps**2*N))
  
  D0term = D0*(RR**2/xs_SH**2 - 4._rkind*ZZ**2/xs_SH**2)
  D1term = D1*(RR**4/xs_SH**4 - 12._rkind*RR**2*ZZ**2/xs_SH**4 + 8._rkind*ZZ**4/xs_SH**5)
  psi_SH = 0.5_rkind*B0*RR**2*(1._rkind-RR**2/xs_SH**2 - ZZ**2/b_SH**2) + 0.5_rkind*B1*RR**2*(1._rkind+D0term+D1term)

  PP_SH = 4._rkind*B0/xs_SH**2

end subroutine Steinhauer

subroutine deallocate_arrays
  use globals
  implicit none

  deallocate(R,Z,RR,ZZ,psiguess,psi)

end subroutine deallocate_arrays

subroutine save
  use globals
  implicit none
  integer :: mp

  mp = 101

  call openbin(mp,'FRC.bin','unformatted','write','big_endian')
  open(mp+1, file = 'FRC.dat', FORM = 'formatted', action = 'write')

  call matwrtI1(mp,'nr',nr)
  call matwrtI1(mp,'nz',nz)
  call matwrtM1(mp,'length',length)
  call matwrtM1(mp,'psiedge',psiedge)
  call matwrtM1(mp,'xs',xs)
  call matwrtM1(mp,'b',b)

  call matwrtM(mp,'R',nrtot,1,R)
  call matwrtM(mp,'Z',nztot,1,Z)
  call matwrtM(mp,'RR',nrtot,nztot,RR)
  call matwrtM(mp,'ZZ',nrtot,nztot,ZZ)
  call matwrtM(mp,'psiguess',nrtot,nztot,psiguess)
  call matwrtM(mp,'psi',nrtot,nztot,psi)

  close(mp)
  close(mp+1)

end subroutine save

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

! subroutine ADI_Step_mpi(omega,PP,psiold,psinew)
!   use globals
!   use mpi
!   implicit none
!   real(rkind),intent(in) :: omega
!   real(rkind),intent(in) :: PP
!   real(rkind), dimension(nrtot,nztot), intent(in) :: psiold
!   real(rkind), dimension(nrtot,nztot), intent(out) :: psinew

!   real(rkind), dimension(nr) :: DZpsiold
!   real(rkind), dimension(nrtot,nztot) :: psi_l
!   real(rkind), dimension(nz,nr) :: psi_t
!   real(rkind), dimension(nrtot,nztot) :: psistar
!   real(rkind), dimension(nz) :: DRpsistar
!   real(rkind), dimension(nrtot,nztot) :: cur

!   real(rkind), dimension(nr-1) :: ar
!   real(rkind), dimension(nr  ) :: br
!   real(rkind), dimension(nr-1) :: cr
!   real(rkind), dimension(nr  ) :: dr
!   real(rkind), dimension(nr  ) :: xr

!   real(rkind), dimension(nz-1) :: az
!   real(rkind), dimension(nz  ) :: bz
!   real(rkind), dimension(nz-1) :: cz
!   real(rkind), dimension(nz  ) :: dz
!   real(rkind), dimension(nz  ) :: xz

!   integer :: ir,iz
!   real(rkind) :: dr1,dz1
!   real(rkind) :: dr2,dz2

!   integer :: istart,iend
!   integer :: id,Nid,Nidsend
!   integer, dimension(0:numprocs-1) :: Nids,deplts

!   dr1 = 1._rkind/deltar
!   dz1 = 1._rkind/deltaz
!   dr2 = dr1**2
!   dz2 = dz1**2

!   ! Psistar

!   cur = 0._rkind

!   do ir=1,nrtot
!      do iz=1,nztot
!         if (psiold(ir,iz).gt.0) then
!            cur(ir,iz) = -RR(ir,iz)**2*PP
!         end if
!      end do
!   end do

!   ar = 0._rkind
!   br = 1._rkind+omega*2._rkind*dr2
!   cr = 0._rkind

!   do ir=1,nr-1
!      ar(ir) = -omega*(dr2 + 0.5_rkind*dr1/R(ir+2))
!      cr(ir) = -omega*(dr2 - 0.5_rkind*dr1/R(ir+1))
!   end do

!   call endstart(nz,istart,iend)
!   istart = istart+1
!   iend = iend+1
!   Nid = (iend-istart+1)*nr

!   do id=0,numprocs-1
!      Nidsend = Nid
!      call MPI_BCAST(Nidsend,1,MPI_INTEGER,id,MPI_COMM_WORLD,code)
!      Nids(id) = Nidsend
!   end do
!   deplts(0) = 0
!   do id=1,numprocs-1
!      deplts(id) = deplts(id-1) + Nids(id-1)
!   end do

!   psi_l = 0._rkind
!   do iz=istart,iend
!      DZpsiold = (psiold(2:nr+1,iz+1) + psiold(2:nr+1,iz-1) - 2._rkind*psiold(2:nr+1,iz))*dz2
!      dr(:) = psiold(2:nr+1,iz) + omega*DZpsiold - omega*cur(2:nr+1,iz)
!      dr(nr) = dr(nr) + omega*psiedge*dr1*(dr1 - 0.5/R(nr+1))
!      call TDMA_Solver(ar,br,cr,dr,xr,nr)
!      psi_l(2:nr+1,iz) = xr(:)
!   end do  

!   psistar = 0._rkind
!   call MPI_ALLGATHERV(psi_l(2:nr+1,istart:iend),Nid,MPI_DOUBLE_PRECISION, &
!        &              psistar(2:nr+1,2:nz+1),Nids,deplts,MPI_DOUBLE_PRECISION, &
!        &              MPI_COMM_WORLD,code)

!   psistar(1,:) = 0._rkind
!   psistar(nr+2,:) = psiedge
!   psistar(:,1) = psiedge*R**2
!   psistar(:,nz+2) = psiedge*R**2

!   ! Psinew

!   cur = 0._rkind

!   do ir=1,nrtot
!      do iz=1,nztot
!         if (psistar(ir,iz).gt.0) then
!            cur(ir,iz) = -RR(ir,iz)**2*PP
!         end if
!      end do
!   end do

!   az = -omega*dz2
!   bz = 1._rkind+omega*2._rkind*dz2
!   cz = -omega*dz2

!   call endstart(nr,istart,iend)
!   istart = istart+1
!   iend = iend+1
!   Nid = (iend-istart+1)*nz

!   do id=0,numprocs-1
!      Nidsend = Nid
!      call MPI_BCAST(Nidsend,1,MPI_INTEGER,id,MPI_COMM_WORLD,code)
!      Nids(id) = Nidsend
!   end do
!   deplts(0) = 0
!   do id=1,numprocs-1
!      deplts(id) = deplts(id-1) + Nids(id-1)
!   end do

!   psi_l = 0._rkind
!   do ir=istart,iend
!   !do ir=2,nr+1
!      DRpsistar = (psistar(ir+1,2:nz+1) + psistar(ir-1,2:nz+1) - 2._rkind*psistar(ir,2:nz+1))*dr2 - &
!           &      0.5_rkind*(psistar(ir+1,2:nz+1) - psistar(ir-1,2:nz+1))*dr1/R(ir)
!      dz(:) = psistar(ir,2:nz+1) + omega*DRpsistar - omega*cur(ir,2:nz+1)
!      dz(1) = dz(1) + omega*psiedge*dz2*R(ir)**2
!      dz(nz) = dz(nz) + omega*psiedge*dz2*R(ir)**2
!      call TDMA_Solver(az,bz,cz,dz,xz,nz)
!      psi_l(ir,2:nz+1) = xz(:)
!   end do

!   psi_t = 0._rkind
!   psinew = 0._rkind
!   call MPI_ALLGATHERV(transpose(psi_l(istart:iend,2:nz+1)),Nid,MPI_DOUBLE_PRECISION, &
!        &              psi_t,Nids,deplts,MPI_DOUBLE_PRECISION, &
!        &              MPI_COMM_WORLD,code)

!   psinew(2:nr+1,2:nz+1) = transpose(psi_t)
  
!   psinew(1,:) = 0._rkind
!   psinew(nr+2,:) = psiedge
!   psinew(:,1) = psiedge*R**2
!   psinew(:,nz+2) = psiedge*R**2

! end subroutine ADI_Step_mpi


! subroutine endstart(ntasks,istart,iend)
!   use globals
!   implicit none
!   integer, intent(in) :: ntasks
!   integer, intent(out) :: istart,iend

!   integer :: ntasks_per_proc
!   integer :: ntasks_remainder
!   integer :: nt

!   integer :: i,k

!   ntasks_per_proc = ntasks/numprocs
!   ntasks_remainder = ntasks - ntasks_per_proc*numprocs

!   istart = 1 + ntasks_per_proc*myid
!   iend = istart + ntasks_per_proc - 1
  
!   k=0                      
!   do nt=1,ntasks_remainder 
!      if (myid.ge.k) then   
!         iend = iend+1      
!      end if
!      k=k+1                 
!      if (myid.ge.k) then   
!         istart = istart+1  
!      end if
!   end do

! end subroutine endstart
