program test

use gmres_module

implicit none

integer, parameter :: n=10 ! dimension of the problem
real(8) :: matrix(n,n), b(n), x(n), invA(n,n)
real(8) :: metric, tol
logical :: print_matrix
integer :: m, max_it, i, j, nout, seed(1)


! make A and b for equation A*x=b
seed(1) = 123
call random_seed(put=seed)
call random_number(matrix)
call random_number(b)

! print matrix to check
print_matrix = .true.
if (print_matrix) then
  print *,'b=',b
  print *,' ------- matrix ------'
  do i = 1, 10
    print *,'column ',i
    do j = 1, 10
      print *, matrix(j, i)
    end do
    print *, ""
  end do
  print *,' --------- end of matrix ------'
endif


! call GMRES()
max_it = 10      ! number of restarts
tol    = 1e-10
metric = 0.01
m      = 20      ! number of GMRES iterations
x      = 0.d0   ! initial x
nout   = 6      ! file number, 6 for stdout
call gmres(n, matvec, x, b, m, max_it, tol, nout)
print *,'gmres x=',x


! benchmark
invA = matrix
call mat_inv_lapack(n,invA)
print *,'benchmark=',matmul(invA,b)/metric


contains

! user defined function to compute A*x
subroutine matvec(n,x,Ax)
  integer :: n
  real(dp) :: x(n)
  real(dp) :: Ax(n)
  Ax = matmul(matrix, x)*metric
end subroutine

end program
