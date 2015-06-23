! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution

! Methods for solving tridiagonal systems using the Thomas algorithm
module tridiag
  use progvars, only: dp

  implicit none

  private
  public :: tridiag_general
  public :: tridiag_cnst
  public :: tridiag_init
  public :: tridiag_cleanup

contains

  ! Initialize tridiagonal solver work arrays
  subroutine tridiag_init(n, mat_coeff, vec_coeff)
    implicit none

    integer(dp) :: n
    complex(dp), allocatable :: mat_coeff(:), vec_coeff(:)

    allocate(mat_coeff(n - 1))
    allocate(vec_coeff(n))

  end subroutine tridiag_init

  ! Cleanup tridiagonal solver work arrays
  subroutine tridiag_cleanup(mat_coeff, vec_coeff)
    implicit none

    complex(dp), allocatable :: mat_coeff(:), vec_coeff(:)

    deallocate(mat_coeff)
    deallocate(vec_coeff)

  end subroutine tridiag_cleanup

  ! General backwards sweep for Thomas algorithm
  subroutine tridiag_backsweep(res, mat_coeff, vec_coeff)
    implicit none

    complex(dp), intent(inout) :: res(:)
    complex(dp), intent(in) :: mat_coeff(:), vec_coeff(:)
    integer(dp) :: i, n

    n = size(res)

    res(n) = vec_coeff(n)

    do i = n - 1, 1, -1
       res(i) = vec_coeff(i) - mat_coeff(i) * res(i + 1)
    end do

  end subroutine tridiag_backsweep

  ! Solves an n-dimensional tridiagional matrix equation A x = b.
  subroutine tridiag_general(diag, u_diag, l_diag, vec, res, &
       mat_coeff, vec_coeff)
    implicit none

    complex(dp), intent(in) :: diag(:), u_diag(:), l_diag(:), vec(:)
    complex(dp), intent(inout) :: mat_coeff(:), vec_coeff(:)
    complex(dp), intent(out) :: res(:)

    integer(dp) :: i, n

    n = size(res)

    ! Forward sweep
    mat_coeff(1) = u_diag(1) / diag(1)
    vec_coeff(1) = vec(1) / diag(1)

    do i = 2, n - 1
       mat_coeff(i) = u_diag(i) / (diag(i) - l_diag(i) * mat_coeff(i - 1))
       vec_coeff(i) = (vec(i) - l_diag(i) * vec_coeff(i - 1)) / &
            (diag(i) - l_diag(i) * mat_coeff(i - 1))
    end do

    vec_coeff(n) = (vec(n) - l_diag(n) * vec_coeff(n - 1)) / &
         (diag(n) - l_diag(n) * mat_coeff(n - 1))

    ! Backward sweep
    call tridiag_backsweep(res, mat_coeff, vec_coeff)

  end subroutine tridiag_general

  ! Solves an n-dimensional tridiagional matrix equation A x = b, where the
  ! band elements are constant.
  subroutine tridiag_cnst(diag_cnst, u_diag_cnst, l_diag_cnst, vec, res, &
       mat_coeff, vec_coeff)
    implicit none

    complex(dp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst, vec(:)
    complex(dp), intent(inout) :: mat_coeff(:), vec_coeff(:)
    complex(dp), intent(out) :: res(:)

    integer(dp) :: i, n

    n = size(res)

    ! Forward sweep
    mat_coeff(1) = u_diag_cnst / diag_cnst
    vec_coeff(1) = vec(1) / diag_cnst

    do i = 2, n -1
       mat_coeff(i) = u_diag_cnst / (diag_cnst - l_diag_cnst * mat_coeff(i - 1))
       vec_coeff(i) = (vec(i) - l_diag_cnst * vec_coeff(i - 1)) / &
            (diag_cnst - l_diag_cnst * mat_coeff(i - 1))
    end do
    vec_coeff(n) = (vec(n) - l_diag_cnst * vec_coeff(n - 1)) / &
         (diag_cnst - l_diag_cnst * mat_coeff(n - 1))

    ! Backward sweep
    call tridiag_backsweep(res, mat_coeff, vec_coeff)

  end subroutine tridiag_cnst

end module tridiag
