! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution

! Methods for solving tridiagonal systems
module tridiag
  use progvars, only: dp
  
  implicit none

  private
  public :: tridiag_sym_cnst

  complex(dp), allocatable :: mat_coeff(:), vec_coeff(:)
  integer(dp) :: n
  
contains

  ! Initialize tridiagonal solver
  ! Allocates work arrays for Thomas algorithm
  subroutine tridiag_init(diag_arr)
    implicit none

    complex(dp) :: diag_arr(:)

    n = size(diag_arr)
    allocate(mat_coeff(n - 1))
    allocate(vec_coeff(n))

  end subroutine tridiag_init

  ! Cleanup tridiagonal solver
  subroutine tridiag_cleanup()
    implicit none

    deallocate(mat_coeff)
    deallocate(vec_coeff)

  end subroutine tridiag_cleanup

  ! General backwards sweep for Thomas algorithm
  subroutine tridiag_backsweep(res, vec_coeff)
    implicit none

    complex(dp), intent(inout) :: res(n)
    complex(dp), intent(in) :: vec_coeff(n)
    integer(dp) :: i
    
    res(n) = vec_coeff(n)

    do i = n - 1, 1, -1
       res(i) = vec_coeff(i) - mat_coeff(i) * res(i + 1)
    end do
    
  end subroutine tridiag_backsweep

  ! Solves an n-dimensional tridiagional matrix equation A x = b, where the
  ! upper and lower diagional matrix elements are equal and constant.
  ! This uses a simplified version of the Thomas algorithm for a general
  ! n-dimensional tridiagonal system.
  subroutine tridiag_sym_cnst(diag, sym_cnst, vec, res)
    implicit none

    complex(dp), intent(in) :: diag(:), sym_cnst, vec(:)
    complex(dp), intent(out) :: res(:)

    integer(dp) :: i

    call tridiag_init(diag)

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
    call tridiag_backsweep(res, vec_coeff)
    call tridiag_cleanup()
    
  end subroutine
end module tridiag
