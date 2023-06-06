subroutine polyfit(x, y, deg, coeffs)
  use prec_const
  use globals
  implicit none
  
  ! Input parameters
  real(rkind), dimension(:), intent(in) :: x
  real(rkind), dimension(:), intent(in) :: y
  integer, intent(in) :: deg
  
  ! Output parameter
  real(rkind), dimension(0:deg), intent(out) :: coeffs
  
  ! Local variables
  integer :: n, m, i, j, k
  real(rkind), dimension(:,:), allocatable :: A, Q, R
  real(rkind), dimension(:), allocatable :: b, c

  n = size(x)
  m = deg + 1
  
  ! Allocate memory for matrices
  allocate(A(1:n, 1:m))
  allocate(Q(1:n, 1:m))
  allocate(R(1:m, 1:m))
  allocate(b(1:n))
  allocate(c(1:m))
  
  ! Build the Vandermonde matrix
  do i = 1, n
    do j = 1, m
      A(i, j) = x(i)**(j-1)
    end do
  end do
  
  ! Perform QR factorization
  call qr_factorization(A, Q, R)
  
  ! Calculate the least squares solution
  b = matmul(transpose(Q), y)
  call qr_solve(R, b, c)
  
  ! Store the coefficients
  coeffs = c
  
  ! Deallocate memory
  deallocate(A, Q, R, b, c)
  
end subroutine polyfit

subroutine qr_factorization(A, Q, R)
  implicit none
  
  ! Input parameters
  real(rkind), dimension(:,:), intent(in) :: A
  
  ! Output parameters
  real(rkind), dimension(:,:), intent(out) :: Q
  real(rkind), dimension(:,:), intent(out) :: R
  
  ! Local variables
  integer :: n, m, i, j, k
  real(rkind), dimension(:), allocatable :: u, v
  
  n = size(A, 1)
  m = size(A, 2)
  
  ! Allocate memory for matrices
  allocate(Q(1:n, 1:m))
  allocate(R(1:m, 1:m))
  allocate(u(1:n))
  allocate(v(1:n))
  
  ! Initialize Q as identity matrix
  Q = 0.0
  do i = 1, min(n, m)
    Q(i, i) = 1.0
  end do
  
  ! Perform QR factorization
  R = A
  do k = 1, min(n, m)
    u = R(k:n, k)
    u(1) = u(1) + sign(sqrt(sum(u**2)), u(1))
    u = u / sqrt(sum(u**2))
    R(k:n, k:m) = R(k:n, k:m) - 2.0 * matmul(reshape(u, shape(u, 1), 1), matmul(transpose(reshape(u, shape(u, 1), 1)), R(k:n, k:m)))
    if (k < n) then
      v = R(k, k+1:m)
      v = v - 2.0 * matmul(reshape(u, shape(u, 1), 1), matmul(transpose(reshape(u, shape(u, 1), 1)), v))
      R(k, k+1:m) = v
    end if
    Q(:, k:n) = Q(:, k:n) - 2.0 * matmul(Q(:, k:n), matmul(transpose(reshape(u, shape(u, 1), 1)), R(k:n, k:m)))
  end do
  
  ! Deallocate memory
  deallocate(u, v)
  
end subroutine qr_factorization

subroutine qr_solve(A, b, x)
  implicit none
  
  ! Input parameters
  real(rkind), dimension(:,:), intent(in) :: A
  real(rkind), dimension(:), intent(in) :: b
  
  ! Output parameter
  real(rkind), dimension(:), intent(out) :: x
  
  ! Local variables
  integer :: n, m, i, j, k
  real(rkind), dimension(:), allocatable :: y
  
  n = size(A, 1)
  m = size(A, 2)
  
  ! Allocate memory for vectors
  allocate(y(1:m))
  
  ! Solve the upper triangular system R * y = Q^T * b
  y = matmul(transpose(A(:, 1:m)), b)
  do k = m, 1, -1
    y(k) = y(k) / A(k, k)
    y(1:k-1) = y(1:k-1) - A(1:k-1, k) * y(k)
  end do
  
  ! Assign solution
  x = y
  
  ! Deallocate memory
  deallocate(y)
  
end subroutine qr_solve
