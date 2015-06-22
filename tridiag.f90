! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution

! Methods for solving tridiagonal systems
module tridiag
  use progvars, only: dp
  
  implicit none

  private
  public :: tridiag_sym_cnst

contains

  ! Solves an n-dimensional tridiagional matrix equation A x = b, where the
  ! upper and lower diagional matrix elements are equal and constant.
  ! This uses a simplified version of the Thomas algorithm for a general
  ! n-dimensional tridiagonal system.
  subroutine tridiag_sym_cnst(diag, sym_cnst, vec, res)
    implicit none
    ! todo: size checking - len(diag) = len(vec) = len(res)
    complex(dp), intent(in) :: diag(:), sym_cnst, vec(:)
    complex(dp), intent(out) :: res(:)

    complex(dp), allocatable :: mat_coeff(:), vec_coeff(:)
    integer(dp) :: n, i

    n = size(diag)
    allocate(mat_coeff(n), vec_coeff(n))

    ! Forward sweep
    mat_coeff(1) = sym_cnst / diag(1)
    vec_coeff(1) = vec(1) / diag(1)

    do i = 2, n - 1
       mat_coeff(i) = sym_cnst / (diag(i) - sym_cnst * mat_coeff(i - 1))
       vec_coeff(i) = (vec(i) - sym_cnst * vec_coeff(i - 1)) / &
            (diag(i) - sym_cnst * mat_coeff(i - 1))
    end do

    vec_coeff(n) = (vec(n) - sym_cnst * vec_coeff(n - 1)) / &
         (diag(n) - sym_cnst * mat_coeff(n - 1))

    ! Backward sweep
    res(n) = vec_coeff(n)

    do i = n - 1, 1, -1
       res(i) = vec_coeff(i) - mat_coeff(i) * res(i + 1)
    end do

    deallocate(mat_coeff, vec_coeff)
    
  end subroutine 
end module tridiag
