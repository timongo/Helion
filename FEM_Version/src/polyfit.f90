subroutine polyfit(vx, vy, d, coeff)
  use prec_const 
  implicit none

  integer, intent(in) :: d
  real(rkind), dimension(:), intent(in) :: vx, vy
  real(rkind), dimension(d+1), intent(out) :: coeff

  real(rkind), dimension(:,:), allocatable :: X
  real(rkind), dimension(:,:), allocatable :: XT_X
  real(rkind), dimension(:), allocatable :: XT_y
  real(rkind), dimension(d+1, d+1) :: A
  real(rkind), dimension(d+1) :: b

  integer :: i, j, n

  n = size(vx)

  allocate(X(n, d+1))
  allocate(XT_X(d+1, d+1))
  allocate(XT_y(d+1))

  ! Prepare the design matrix X
  do i = 0, d
     X(:, i+1) = vx**i
  end do

  XT_X = matmul(transpose(X), X)
  XT_y = matmul(transpose(X), vy)

  ! Compute the coefficients using the method of least squares
  A = 0.0
  b = 0.0

  do i = 1, d+1
     do j = 1, d+1
        A(i, j) = sum(XT_X(:, i) * XT_X(:, j))
     end do

     b(i) = sum(XT_y * XT_X(:, i))
  end do

  ! Solve the linear system A * coeff = b
  call gauss_elimination(A, b, d+1, coeff)

  ! Deallocate memory
  deallocate(X)
  deallocate(XT_X)
  deallocate(XT_y)

end subroutine polyfit

subroutine gauss_elimination(A, b, n, x)
  use prec_const
  implicit none

  integer, intent(in) :: n
  real(rkind), dimension(n, n) :: A
  real(rkind), dimension(n), intent(in) :: b
  real(rkind), dimension(n), intent(out) :: x

  real(rkind) :: factor
  integer :: i, j, k

  x = b

  do k = 1, n-1
     do i = k+1, n
        factor = A(i, k) / A(k, k)
        x(i) = x(i) - factor * x(k)
        do j = k+1, n
           A(i, j) = A(i, j) - factor * A(k, j)
        end do
     end do
  end do

  do k = n, 1, -1
     x(k) = x(k) / A(k, k)
     do i = k-1, 1, -1
        x(i) = x(i) - A(i, k) * x(k)
     end do
  end do

end subroutine gauss_elimination
