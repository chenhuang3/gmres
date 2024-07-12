!
! Converted from m-file: https://www.netlib.org/templates/matlab/gmres.m
!
! Major modifications:
!
!  1. A*x is replaced by an external user-defined subroutine (matvec)
!  2. Norm and dot product of vectors can take a metric
!
! Chen Huang (Department of Scientific Computing, Florida State University)
!
module gmres_module
  implicit none
  integer, parameter :: dp = kind(1.0d0)

contains

  subroutine gmres(n, matvec, x, b, metric, m, max_it, tol, nout)
    !
    ! INPUTS:
    !  n       --  integer, dimension of x
    !  x       --  real(8) x(n), initial guess, on exit x is the solution
    !  m       --  integer, number of iterations for each GMRES bewteen restarts
    !  max_it  --  integer, number of restarts
    !  matvec  --  external function computes A*x, defined as matvec(n, x, z), where z=A*x.
    !  metric  --  real(8), for dot product of two vectors, sum(v1*v2)*metric, for simple linear algebra, set metric=1.d0
    !  nout    --  integer, file number for print information
    !
    !--------------------------------------------
    ! original description from gmres.m file
    !--------------------------------------------
    !%  -- Iterative template routine --
    !%     Univ. of Tennessee and Oak Ridge National Laboratory
    !%     October 1, 1993
    !%     Details of this algorithm are described in "Templates for the
    !%     Solution of Linear Systems: Building Blocks for Iterative
    !%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
    !%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
    !%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
    !%
    !% [x, error, iter, flag] = gmres( A, x, b, M, restrt, max_it, tol )
    !%
    !% gmres.m solves the linear system Ax=b
    !% using the Generalized Minimal residual ( GMRESm ) method with restarts .
    !%
    !% input   A        REAL nonsymmetric positive definite matrix
    !%         x        REAL initial guess vector
    !%         b        REAL right hand side vector
    !%         M        REAL preconditioner matrix
    !%         restrt   INTEGER number of iterations between restarts
    !%         max_it   INTEGER maximum number of iterations
    !%         tol      REAL error tolerance
    !%
    !% output  x        REAL solution vector
    !%         error    REAL error norm
    !%         iter     INTEGER number of iterations performed
    !%         flag     INTEGER: 0 = solution found to tolerance
    !%                           1 = no convergence given max_it

    implicit none

    ! interface of the external subroutine to compute A*x
    interface
      subroutine matvec(n, x, Ax)
        integer :: n
        real(8) :: x(n), Ax(n)
      end subroutine matvec
    end interface

    integer  :: n, m, max_it, iter, nout
    real(dp) :: tol, metric
    real(dp) :: x(n), b(n), vec(n), error
    integer  :: i, j, k
    real(dp) :: bnrm2, rnorm, temp

    ! note that dim of e1 should be m+1, rather than n in the m-file
    real(dp) :: r(n), w(n), y(m), s(m+1), cs(m), sn(m), e1(m+1)
    real(dp) :: V(n,m+1), H(m+1,m)

    bnrm2 = norm2(b, metric)
    if (bnrm2 == 0.0_dp) bnrm2 = 1.0_dp

    call matvec(n, x, vec)
    r = b - vec
    rnorm = norm2(r, metric)
    error = rnorm / bnrm2
    if (error < tol) return

    e1    = 0.0_dp
    e1(1) = 1.0_dp

    !------------------------------------
    !------------- Restart --------------
    !------------------------------------
    do iter = 1, max_it
      write(nout,'(a,i4)') 'restart iteration: ',iter
      call matvec(n, x, vec)
      r = b - vec
      V(:,1) = r / norm2(r, metric)
      s = norm2(r, metric) * e1

      !-------------------------------------------
      !------------ GMRES iterations -------------
      !-------------------------------------------
      do i = 1, m
        call matvec(n, V(:,i), w)
        do k = 1, i
          H(k,i) = dot_product(w, V(:,k)) * metric
          w = w - H(k,i) * V(:,k)
        end do
        H(i+1,i) = norm2(w, metric)
        V(:,i+1) = w / H(i+1,i)

        do k = 1, i-1
          temp     =  cs(k) * H(k,i) + sn(k) * H(k+1,i)
          H(k+1,i) = -sn(k) * H(k,i) + cs(k) * H(k+1,i)
          H(k,i)   = temp
        end do

        call rotmat(H(i,i), H(i+1,i), cs(i), sn(i))
        temp   = cs(i) * s(i)
        s(i+1) = -sn(i) * s(i)
        s(i)   = temp
        H(i,i) = cs(i) * H(i,i) + sn(i) * H(i+1,i)
        H(i+1,i) = 0.d0

        error = abs(s(i+1)) / bnrm2
        write(nout,'(a,i5,a,es11.3)')'gmres_iter:',i,'  error=',error

        if (error <= tol) then
          y(1:i) = s(1:i)
          call solve(i, H(1:i,1:i), y(1:i))
          x = x + matmul(V(:,1:i), y(1:i))
          exit
        end if
      end do

      if (error <= tol) exit

      ! compute solution x
      y = s(1:m)
      call solve(m, H(1:m,1:m), y(1:m))
      x = x + matmul(V(:,1:m), y(1:m))

      ! compute error
      call matvec(n, x, vec)
      r = b - vec
      error = norm2(r, metric) / bnrm2
      if (error <= tol) exit

    end do

  end subroutine gmres


  ! Based on code: http://www.netlib.org/templates/matlab/rotmat.m
  subroutine rotmat(a,b,c,s)
    implicit none
    real(dp) :: a, b, c, s, temp
    if (b==0.0d0) then
      c = 1.d0
      s = 0.d0
    elseif (abs(b) > abs(a)) then
      temp = a / b
      s = 1.d0 / sqrt(1.d0 + temp**2)
      c = temp * s
    else
      temp = b / a
      c = 1.d0 / sqrt(1.d0 + temp**2)
      s = temp * c
    endif
  end subroutine rotmat


  real(dp) function norm2(x, metric)
    implicit none
    real(dp), dimension(:), intent(in) :: x
    real(dp) :: metric
    norm2 = sqrt(sum(x * x) * metric)
  end function norm2


  ! compute H^{-1} * y
  subroutine solve(n, H, y)
    implicit none
    integer  :: n
    real(dp) :: H(n,n), y(n), H_inv(n,n)

    H_inv = H
    call mat_inv_lapack(n,H_inv)
    y = matmul(H_inv, y)
  end subroutine solve


  ! compute A^{-1}
  subroutine mat_inv_lapack(n,A)
    implicit none
    integer  :: n, info
    integer  :: ipiv(n)
    real(dp) :: A(n,n), work(n)

    call dgetrf(n, n, A, n, ipiv, info)
    if (info /= 0) then
      print *, "Error in LU decomposition: ", info
      stop
    end if

    call dgetri(n, A, n, ipiv, work, n, info)
    if (info /= 0) then
      print *, "Error in matrix inversion: ", info
      stop
    end if
  end subroutine mat_inv_lapack

end module gmres_module
