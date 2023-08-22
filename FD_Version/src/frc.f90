program main
  use prec_const
  use globals
  implicit none

  call readnamelist
  call initialization  
  call run
  !call run_psimax
  call save
  call deallocate_arrays

end program main

subroutine run_psimax
  use globals
  implicit none

  real(rkind) :: psimaxerror
  integer :: iter, max_iter
  real(rkind) :: step_factor
  ! 'psimax', the target value, is read from nlfrc
  
  AP = AP_NL
  step_factor = 0.8 ! multiplicative factor for step size adjustment
  
  iter = 0
  max_iter = 50
  psimaxerror = 2*tol
  do while (iter < max_iter .and. abs(psimaxerror) >= tol)
     call ADI_Solve_Current(psiguess,psi)
     ! Update psiguess
     psiguess = psi

     ! Calculate the error in the maximum value of psi wrt target psimax
     psimaxerror = maxval(psi) - psimax

     ! Update CurrentTarget based on the error and step factor
     CurrentTarget = CurrentTarget - (psimaxerror / psimax) * CurrentTarget * step_factor
     
     iter = iter + 1
  end do

  if (iter >= max_iter) then
     print*, 'Maximum iterations reached without convergence for target psimax = ', psimax
  else
     print*, 'Converged after ', iter, ' iterations for target psimax = ', psimax
  end if
  
end subroutine run_psimax

subroutine run
  use globals
  implicit none
  
  AP = AP_NL
  call ADI_Solve_Current(psiguess,psi)

  ! real(rkind), dimension(10) :: dAP
  ! real(rkind) :: norm
  ! integer :: i,n
  
  ! dAP = AP_NL - AP_guess

  ! call DiffNorm(10,AP_NL,AP_guess,norm)

  ! if (norm.lt.0.05_rkind) then
  !    AP = AP_NL
  !    call ADI_Solve_Current(psiguess,psi)
  ! else
  !    n = int(norm*20._rkind)+1
  !    write(*,'(A,I3,A)') 'Will try to reach the target in ',n,' steps'
  !    do i=1,n
  !       AP = AP_guess + dAP*real(i,rkind)/real(n,rkind)
  !       write(*,'(A,10E12.4)') 'Solving with AP=',AP
  !       call ADI_Solve_Current(psiguess,psi)
  !       psiguess = psi
  !    end do
  ! end if
  
end subroutine run

function ppfun(psiv) result(result_val)
  ! Pprime function in Grad-Shafranov equation
  use prec_const
  use globals, only : AP, psimax, pp_p
  implicit none
  real(rkind), intent(in) :: psiv
  real(rkind) :: x, s, result_val
  integer :: i
  
  result_val = 0._rkind

  if (0 <= psiv .and. psiv <= psimax) then 
     s = sqrt(1._rkind - psiv/psimax) ! roughly a function of r
     ! x = psiv/psimax
     do i=1,10
        result_val = result_val + AP(i)*s**(i-1)
     end do  
  end if

  ! result_val = 2.*pp_p(1)* tanh(pp_p(1)*psiv+pp_p(2)) / cosh(pp_p(1)*psiv+pp_p(2))**2
end function ppfun

subroutine ADI_Solve_Current(psig,psi2)
  use prec_const
  use globals
  implicit none

  real(rkind), dimension(nztot,nrtot), intent(inout) :: psig
  real(rkind), dimension(nztot,nrtot), intent(out) :: psi2
  ! real(rkind), intent(in) :: PP

  real(rkind), dimension(nztot,nrtot) :: psi1
  real(rkind), dimension(nztot,nrtot) :: psi2_star
  real(rkind), dimension(nztot,nrtot) :: aux
  real(rkind) :: omega,omegamax
  real(rkind) :: start,finish
  integer :: k
  real(rkind) :: I
  real(rkind) :: I0
  real(rkind) :: Lambda
  real(rkind) :: dtpsi

  real(rkind) :: N
  real(rkind) :: Nstar
  real(rkind) :: norm
  real(rkind) :: ratio
  real(rkind) :: dnrm2
  
  real(rkind), external :: ppfun

  omegamax = 1.e-1_rkind

  omega = omegai

  norm = dnrm2(ntot,psig,1)
  ! print*, norm

  Lambda = LambdaIni
  !C = 1._rkind
  psimaxcur = maxval(psig)
  I0 = CurrentTarget

  k = 0
  dtpsi = tol
  
  do while (k.lt.ntsmax .and. dtpsi.ge.tol)   
     call ADI_Step_omp(omega,Lambda,psig,psi1)
     call ADI_Step_omp(omega,Lambda,psi1,psi2)
     call ADI_Step_omp(omega*2._rkind,Lambda,psig,psi2_star)

     call dcopy(ntot,psig,1,aux,1)
     call daxpy(ntot,-1._rkind,psi1,1,aux,1)
     N = dnrm2(ntot,aux,1)

     call dcopy(ntot,psi2,1,aux,1)
     call daxpy(ntot,-1._rkind,psi2_star,1,aux,1)
     Nstar = dnrm2(ntot,aux,1)

     ratio = Nstar/N

     dtpsi = N/(norm*omega)
     
     if (mod(k+1,iplot).eq.0 .or. dtpsi.lt.tol) then
        write(*,'(I6,5E12.4)') k+1,omega,Lambda,ratio,dtpsi/tol
     end if

     if (ratio.le.0.05_rkind) then
        k=k+1
        call dcopy(ntot,psi2,1,psig,1)
        psimaxcur = maxval(psig)
        omega = min(omega*4._rkind,omegamax)
        call TotalCurrent(psi2,1._rkind,I)
        Lambda = I0/I
     elseif (ratio.gt.0.05_rkind .and. ratio.le.0.1_rkind) then
        k = k+1
        call dcopy(ntot,psi2,1,psig,1)
        psimaxcur = maxval(psig)
        omega = min(omega*2._rkind,omegamax)
        call TotalCurrent(psi2,1._rkind,I)
        Lambda = I0/I        
     elseif (ratio.gt.0.1_rkind .and. ratio.le.0.3_rkind) then
        k = k+1
        psimaxcur = maxval(psig)
        call dcopy(ntot,psi2,1,psig,1)
        call TotalCurrent(psi2,1._rkind,I)
        Lambda = I0/I        
     elseif (ratio.gt.0.3_rkind .and. ratio.le.0.4_rkind) then
        k = k+1
        call dcopy(ntot,psi2,1,psig,1)
        psimaxcur = maxval(psig)
        omega = omega/2._rkind
        call TotalCurrent(psi2,1._rkind,I)
        Lambda = I0/I        
     elseif (ratio.gt.0.4_rkind .and. ratio.le.0.6_rkind) then
        k = k+1
        call dcopy(ntot,psi2,1,psig,1)
        psimaxcur = maxval(psig)
        omega = omega/4._rkind
        call TotalCurrent(psi2,1._rkind,I)
        Lambda = I0/I        
     elseif (ratio.gt.0.6_rkind) then
        k = k+1
        call dcopy(ntot,psi2,1,psig,1)
        psimaxcur = maxval(psig)
        omega = omega/16._rkind
     end if
  end do
  
  call TotalCurrent(psi2,Lambda,I)
  if (k >= ntsmax) then
     print*, 'Maximum iterations reached without convergence for target I = ', I
  else
     print*, 'Converged after ', k, ' iterations for target I = ', I
  end if

  LambdaSol = Lambda
  ! print*, 'Lambda = ', LambdaSol
  print*, 'psimax = ', maxval(psi2)

end subroutine ADI_Solve_Current

subroutine ADI_Solve(psig,psi2)
  ! Alternate direction implicit algorithm
  ! psig is guess
  ! psi2 is converged solution
  ! omega is the time step
  ! lambda is a Lagrange multiplier
  ! I,I0 are total currents
  ! N, Nstar, norm are norms
  use prec_const
  use globals
  implicit none

  real(rkind), dimension(nztot,nrtot), intent(inout) :: psig
  real(rkind), dimension(nztot,nrtot), intent(out) :: psi2
  ! real(rkind), intent(in) :: PP

  real(rkind), dimension(nztot,nrtot) :: psi1
  real(rkind), dimension(nztot,nrtot) :: psi2_star
  real(rkind), dimension(nztot,nrtot) :: aux
  real(rkind) :: omega,omegamax
  real(rkind) :: start,finish
  integer :: k
  real(rkind) :: I
  real(rkind) :: I0
  real(rkind) :: Lambda, LambdaNew, dLambda
  real(rkind) :: relax
  real(rkind) :: dtpsi

  real(rkind) :: N
  real(rkind) :: Nstar
  real(rkind) :: norm
  real(rkind) :: ratio
  real(rkind) :: dnrm2

  relax = 0.8_rkind

  ! Maximum time step
  omegamax = 1.e-1_rkind

  omega = omegai

  norm = dnrm2(ntot,psig,1)

  Lambda = Lambdaini

  k = 0
  dtpsi = tol
  dLambda = tol
  
  print*, 'psitarget = ', psimax
  psimaxcur = maxval(psig)

  ! Do at least 2 iterations
  do while ((k.lt.ntsmax .and. dtpsi.ge.tol) .or. (k.lt.2))
     ! Carrying out two steps of t = omega
     call ADI_Step_omp(omega,Lambda,psig,psi1)
     call ADI_Step_omp(omega,Lambda,psi1,psi2)
     ! Carrying out a double step of t = 2*omega
     call ADI_Step_omp(omega*2._rkind,Lambda,psig,psi2_star)
     
     ! Computing |psig - psi2| in N
     call dcopy(ntot,psig,1,aux,1)
     call daxpy(ntot,-1._rkind,psi2,1,aux,1)
     N = dnrm2(ntot,aux,1)

     ! Computing |psi2 - psi2_star| in Nstar
     call dcopy(ntot,psi2,1,aux,1)
     call daxpy(ntot,-1._rkind,psi2_star,1,aux,1)
     Nstar = dnrm2(ntot,aux,1)

     ! ratio will we used to determine whether omega increases or decreases
     ratio = Nstar/N

     ! tracking convergence in terms of |dpsi/dt|
     dtpsi = N/(norm*omega)

     if (mod(k+1,iplot).eq.0 .or. dtpsi.lt.tol) then
        write(*,'(I6,5E12.4)') k+1,omega,Lambda,ratio,dtpsi/tol, dLambda/tol
     end if

     if (ratio.le.0.05_rkind) then
        k=k+1
        call dcopy(ntot,psi2,1,psig,1)
        ! updating time step
        omega = min(omega*4._rkind,omegamax)
        psimaxcur = maxval(psi2)
        ! LambdaNew = Lambda/(psimaxcur/psimax)
        ! dLambda = LambdaNew - Lambda
        ! Lambda = Lambda + dLambda*relax
        Lambda = Lambda - ((psimaxcur-psimax) / psimax) * Lambda * relax
     elseif (ratio.gt.0.05_rkind .and. ratio.le.0.1_rkind) then
        k = k+1
        call dcopy(ntot,psi2,1,psig,1)
        ! updating time step
        omega = min(omega*2._rkind,omegamax)
        psimaxcur = maxval(psi2)
        ! LambdaNew = Lambda/(psimaxcur/psimax)
        ! dLambda = LambdaNew - Lambda
        ! Lambda = Lambda + dLambda*relax
        Lambda = Lambda - ((psimaxcur-psimax) / psimax) * Lambda * relax
     elseif (ratio.gt.0.1_rkind .and. ratio.le.0.3_rkind) then
        k = k+1
        call dcopy(ntot,psi2,1,psig,1)
        psimaxcur = maxval(psi2)
        ! LambdaNew = Lambda/(psimaxcur/psimax)
        ! dLambda = LambdaNew - Lambda
        ! Lambda = Lambda + dLambda*relax
        Lambda = Lambda - ((psimaxcur-psimax) / psimax) * Lambda * relax
     elseif (ratio.gt.0.3_rkind .and. ratio.le.0.4_rkind) then
        k = k+1
        call dcopy(ntot,psi2,1,psig,1)
        ! updating time step
        omega = omega/2._rkind
        psimaxcur = maxval(psi2)
        ! LambdaNew = Lambda/(psimaxcur/psimax)
        ! dLambda = LambdaNew - Lambda
        ! Lambda = Lambda + dLambda*relax
        Lambda = Lambda - ((psimaxcur-psimax) / psimax) * Lambda * relax
     elseif (ratio.gt.0.4_rkind .and. ratio.le.0.6_rkind) then
        k = k+1
        call dcopy(ntot,psi2,1,psig,1)
        ! updating time step
        omega = omega/4._rkind
        psimaxcur = maxval(psi2)
        ! LambdaNew = Lambda/(psimaxcur/psimax)
        ! dLambda = LambdaNew - Lambda
        ! Lambda = Lambda + dLambda*relax
        Lambda = Lambda - ((psimaxcur-psimax) / psimax) * Lambda * relax
     elseif (ratio.gt.0.6_rkind) then
        ! updating time step without changing the starting field
        omega = omega/16._rkind
     end if
  end do

  LambdaSol = Lambda

  print*, 'Lambda = ', LambdaSol

end subroutine ADI_Solve

subroutine TotalCurrent(psi,PP,I)
  ! Computing total current with a not very refined algorithm
  ! If a node has psi>0, a contribution R*pprime*dr*dz is added
  use prec_const
  use globals, only: nrtot,nztot,R,deltar,deltaz
  implicit none
  real(rkind), dimension(nztot,nrtot), intent(in) :: psi
  real(rkind), intent(in) :: PP
  real(rkind), intent(out) :: I
  real(rkind) :: psiv
  integer :: ir,iz
  real(rkind), external :: ppfun

  I = 0._rkind
  
  !$OMP PARALLEL PRIVATE(ir,iz,psiv) &
  !$OMP          SHARED(PP,R)
  !$OMP DO REDUCTION(+:I)
  do iz=1,nztot
     do ir=1,nrtot
        psiv = psi(iz,ir)
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
  real(rkind), dimension(nztot,nrtot), intent(in) :: psiold
  real(rkind), dimension(nztot,nrtot), intent(out) :: psinew

  real(rkind), dimension(nr) :: DZpsiold
  real(rkind), dimension(nztot,nrtot) :: psistar
  real(rkind), dimension(nz) :: DRpsistar
  real(rkind), dimension(nztot,nrtot) :: cur

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
  
  ! !!!!!!!!!!!!!
  ! Psistar
  ! !!!!!!!!!!!!!

  ! Current term, rhs of GS equation, R jphi = R**2 pprime
  call CurrentTerm(PP,psiold,cur)

  ! ar is left of diagnoal of matrix
  ar = 0._rkind
  ! br is diagnoal of matrix
  br = 1._rkind+omega*2._rkind*dr2
  ! cr is right of diagnoal of matrix
  cr = 0._rkind

  do ir=1,nr-1
     ar(ir) = -omega*(dr2 + 0.5_rkind*dr1/R(ir+2))
     cr(ir) = -omega*(dr2 - 0.5_rkind*dr1/R(ir+1))
  end do
  
  ! psistar corresponds to the first half step where only radial operator is implicit
  psistar = 0._rkind  
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP          PRIVATE(iz,DZpsiold,dr,xr)
  !$OMP DO
  do iz=2,nz+1
     DZpsiold = (psiold(iz+1,2:nr+1) + psiold(iz-1,2:nr+1) - 2._rkind*psiold(iz,2:nr+1))*dz2
     dr(:) = psiold(iz,2:nr+1) + omega*DZpsiold - omega*cur(iz,2:nr+1)
     dr(nr) = dr(nr) + omega*psiedge*dr1*(dr1 - 0.5/R(nr+1))
     call TDMA_Solver(ar,br,cr,dr,xr,nr)
     psistar(iz,2:nr+1) = xr(:)
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Applying boundary conditions
  psistar(:,1) = 0._rkind
  psistar(:,nr+2) = psiedge
  psistar(1,:) = psiedge*R**2
  psistar(nz+2,:) = psiedge*R**2

  ! !!!!!!!!!!!!!
  ! Psinew
  ! !!!!!!!!!!!!!

  ! Current term
  call CurrentTerm(PP,psistar,cur)

  ! az is left of diagnoal of matrix
  az = -omega*dz2
  ! bz is diagnoal of matrix
  bz = 1._rkind+omega*2._rkind*dz2
  ! cz is right of diagnoal of matrix
  cz = -omega*dz2

  ! psinew corresponds to the second half step where only z-part of the operator is implicit
  psinew = 0._rkind
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP          PRIVATE(ir,DRpsistar,dz,xz)
  !$OMP DO
  do ir=2,nr+1
     DRpsistar = (psistar(2:nz+1,ir+1) + psistar(2:nz+1,ir-1) - 2._rkind*psistar(2:nz+1,ir))*dr2 - &
          &      0.5_rkind*(psistar(2:nz+1,ir+1) - psistar(2:nz+1,ir-1))*dr1/R(ir)
     dz(:) = psistar(2:nz+1,ir) + omega*DRpsistar - omega*cur(2:nz+1,ir)
     dz(1) = dz(1) + omega*psiedge*dz2*R(ir)**2
     dz(nz) = dz(nz) + omega*psiedge*dz2*R(ir)**2
     call TDMA_Solver(az,bz,cz,dz,xz,nz)
     psinew(2:nz+1,ir) = xz(:)
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! boundary conditions
  psinew(:,1) = 0._rkind
  psinew(:,nr+2) = psiedge
  psinew(1,:) = psiedge*R**2
  psinew(nz+2,:) = psiedge*R**2

end subroutine ADI_Step_omp

subroutine CurrentTerm(PP,psio,cur)
  use globals
  implicit none
  real(rkind), intent(in) :: PP
  real(rkind), dimension(nztot,nrtot), intent(in) :: psio
  real(rkind), dimension(nztot,nrtot), intent(out) :: cur
  real(rkind), external :: ppfun
  real(rkind) :: psiv
  integer :: iz,ir
  
  cur = 0._rkind

  do iz=1,nztot
     do ir=1,nrtot
        psiv = psio(iz,ir)
        if (psiv.gt.0) then
           cur(iz,ir) = -RR(iz,ir)**2*PP*ppfun(psiv)
        end if
     end do
  end do

end subroutine CurrentTerm

subroutine TDMA_Solver(a,b,c,d,x,n)
  ! Thomas TriDiagonal Matrix Algorithm
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

  nz = 100
  nr = 100
  ! aspect ratio length/radius
  length = 1._rkind
  
  ! maximum number of iterations
  ntsmax = 500
  ! tolerance for final convergence
  tol = 1.e-8_rkind
  ! initial value of timestep in ADI algorithm
  omegai = 1.e-2_rkind

  ! Psi boundary condition (default corresponds to unit field)
  psiedge = -0.5_rkind

  ! psimax target on axis (not used yet)
  psimax = 0.05_rkind

  ! written with python or with fortran
  guesstype = 1

  ! By default 0, in that case is not used
  LambdaNL = 0._rkind

  ! Parameters Defining Pprime as a function of psi
  pp_p(1) = 130._rkind
  pp_p(2) = 2.3_rkind

  ! Polynomial defining Pprime as a function of s
  AP_NL = 0._rkind
  AP_NL(1) = 0._rkind
  AP_NL(2) = 0.5_rkind
  AP_NL(3) = -1._rkind
  AP_NL(4) = 0.5_rkind

  ! iplot for not printing everytime
  iplot = 50

  ! Current Target
  CurrentTarget = -30._rkind

  ! Number of points for psi,theta mesh (jacobian)
  npsi = 50
  ntheta = 64

  ! define namelist
  namelist /frc/ &
       &  nr,nz,length,psiedge,ntsmax,tol,omegai, &
       &  psimax,guesstype,LambdaNL,pp_p,AP_NL,iplot, &
       &  CurrentTarget, npsi, ntheta

  ! read namelist
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
  ! allocate arrays
  use globals
  implicit none
  integer :: i

  allocate(R(nrtot),Z(nztot),RR(nztot,nrtot),ZZ(nztot,nrtot), &
       &   psi(nztot,nrtot),psiguess(nztot,nrtot))

  ! define R arrays and increment
  deltar = 1._rkind/real(nrtot-1,rkind)
  do i=1,nrtot
     R(i) = real(i-1,rkind)/real(nrtot-1,rkind)
     RR(:,i) = real(i-1,rkind)/real(nrtot-1,rkind)
  end do

  ! define Z arrays and increment
  deltaz = 1._rkind/real(nztot-1,rkind)*length
  do i=1,nztot
     Z(i) = real(i-1,rkind)/real(nztot-1,rkind)*length - length/2._rkind
     ZZ(i,:) = real(i-1,rkind)/real(nztot-1,rkind)*length - length/2._rkind
  end do

  call read_guess

end subroutine initialization

subroutine read_guess
  ! read guess array
  ! if guesstype = 1, it is assumed that it was generated in the routine save of the present program
  ! if guesstype = 2, it is assumed that 
  ! the guess was generated with np.array([...]).tofile('guess.bin') command in python
  use globals
  implicit none
  real(rkind), dimension(3) :: sizes_dpr
  
  integer :: mp

  integer :: i

  if (guesstype.eq.1) then
     
     call openbin(mp,'solution_guess_read_FD.bin','unformatted','read','big_endian')
     read(mp) nzguesstot
     read(mp) nrguesstot
     read(mp) lguess
     read(mp) LambdaIni

     allocate(psiguess_read(nzguesstot,nrguesstot))

     read(mp) psiguess_read
     read(mp) AP_guess

     close(mp)

     dzguess = 1._rkind/real(nzguesstot-1,rkind)*lguess
     drguess = 1._rkind/real(nrguesstot-1,rkind)

     if (lguess.ne.length) then
        print*, "length is not equal to that of the guess. Are you sure you know what you're doing?"
        write(*,'(A,E12.4,A,E12.4)') 'Length = ', length, 'while guess length is', lguess
     end if

  elseif (guesstype.eq.2) then
     mp = 101

     open(mp, file='guess.dat', status='old', access='stream', form='unformatted')
     read(mp) sizes_dpr
     close(mp)
     nzguesstot = int(sizes_dpr(1))
     nrguesstot = int(sizes_dpr(2))
     lguess = sizes_dpr(3)

     print*, nzguesstot,nrguesstot

     dzguess = 1._rkind/real(nzguesstot-1,rkind)*lguess
     drguess = 1._rkind/real(nrguesstot-1,rkind)

     if (lguess.ne.length) then
        print*, "length is not equal to that of the guess. Are you sure you know what you're doing?"
        write(*,'(A,E12.4,A,E12.4)') 'Length = ', length, 'while guess length is', lguess
     end if

     allocate(psiguess_read(nzguesstot,nrguesstot))

     open(mp, file='guess.bin', status='old', access='stream', form='unformatted')
     read(mp) psiguess_read
     close(mp)
  end if

  if (nrtot.eq.nrguesstot .and. nztot.eq.nzguesstot) then
     ! if sizes are the same just reding is enough
     psiguess(:,:) = psiguess_read(:,:)
  else
     ! else interpolation is required
     call interpolate_guess
  end if
  
  deallocate(psiguess_read)

  ! if LambdaNL (NL=namelist) is non zero, replace LambdaIni with LambdaNL
  if (LambdaNL.gt.0_rkind) then
     LambdaIni = LambdaNL
  end if

  hlength = 0.5_rkind*length

end subroutine read_guess

subroutine interpolate_guess
  use globals
  implicit none

  integer :: ir,iz,irguess,izguess
  real(rkind) :: x,y,t,u

  ! simple bilinear interpolation
  do iz=2,nz+1
     x = (Z(iz)+length/2._rkind)/dzguess
     izguess = floor(x)+1
     t = x - real(izguess-1,rkind)
     do ir=2,nr+1
        y = R(ir)/drguess
        irguess = floor(y)+1
        u = y - real(irguess-1,rkind)

        psiguess(iz,ir) = (1._rkind-t)*(1._rkind-u)*psiguess_read(izguess  ,irguess  ) + &
             &                  t     *(1._rkind-u)*psiguess_read(izguess+1,irguess  ) + &
             &            (1._rkind-t)*      u     *psiguess_read(izguess  ,irguess+1) + &
             &                  t     *      u     *psiguess_read(izguess+1,irguess+1)
     end do
  end do

  psiguess(:,1) = 0._rkind
  psiguess(:,nr+2) = psiedge
  psiguess(1,:) = psiedge*R**2
  psiguess(nz+2,:) = psiedge*R**2

end subroutine interpolate_guess

subroutine deallocate_arrays
  ! free memory
  use globals
  implicit none

  deallocate(R,Z,RR,ZZ,psiguess,psi)

end subroutine deallocate_arrays

subroutine save
  ! save results
  use globals
  implicit none
  integer :: mp

  mp = 101

  call openbin(mp,'FRC.bin','unformatted','write','big_endian')
  ! bin file has the data
  open(mp+1, file = 'FRC.dat', FORM = 'formatted', action = 'write')
  ! dat file names and shapes of arrays

  call matwrtI1(mp,'nz',nz) ! I --> Integer
  call matwrtI1(mp,'nr',nr)
  call matwrtM1(mp,'length',length) ! M --> Real
  call matwrtM1(mp,'psiedge',psiedge) 

  call matwrtM(mp,'Z',nztot,1,Z)
  call matwrtM(mp,'R',nrtot,1,R)
  call matwrtM(mp,'ZZ',nztot,nrtot,ZZ)
  call matwrtM(mp,'RR',nztot,nrtot,RR)
  call matwrtM(mp,'psiguess',nztot,nrtot,psiguess)
  call matwrtM(mp,'psi',nztot,nrtot,psi)

  call SaveMesh(mp,psi)

  close(mp)
  close(mp+1)

  ! write guess
  call openbin(mp,'solution_guess_write_FD.bin','unformatted','write','big_endian')
  write(mp) nztot
  write(mp) nrtot
  write(mp) length ! Z/R (R is 1 at the top)
  write(mp) LambdaSol 
  write(mp) psi
  write(mp) AP_NL ! Polynomial of psi/psi_max
  close(mp)

end subroutine save

subroutine SaveMesh(mp,psim)
  use globals
  implicit none
  integer :: mp
  real(rkind), dimension(nztot,nrtot), intent(in) :: psim
  real(rkind), dimension(npsi) :: SMesh,PsiMesh,PressureMesh
  real(rkind), dimension(npsi,ntheta+1) :: ZMesh,RMesh,JacobMesh

  call Mesh(LambdaSol,psim,npsi,ntheta,ZMesh,RMesh,JacobMesh,Smesh,PsiMesh,PressureMesh)
  
  call matwrtM(mp,'S',npsi,1,Smesh)
  call matwrtM(mp,'PsiMesh',npsi,1,PsiMesh)
  call matwrtM(mp,'Pressure',npsi,1,PressureMesh)
  call matwrtM1(mp,'psimax',Psimesh(1))
  call matwrtM(mp,'ZMesh',npsi,ntheta+1,ZMesh)
  call matwrtM(mp,'RMesh',npsi,ntheta+1,RMesh)
  call matwrtM(mp,'JacobMesh',npsi,ntheta+1,JacobMesh)

end subroutine SaveMesh

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

subroutine DiffNorm(n,x,y,norm)
  use prec_const
  implicit none
  integer, intent(in) :: n
  real(rkind), dimension(n), intent(in) :: x,y
  real(rkind), intent(out) :: norm

  norm = norm2(y-x)/norm2(x)

end subroutine DiffNorm

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
