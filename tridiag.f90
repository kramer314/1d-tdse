! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution

! Methods for solving tridiagonal systems using the Thomas algorithm (modified
! Gaussian elimination)
! Note: This is a program-independent module
module tridiag
  use globvars, only: dp

  implicit none

  private
  ! Generic module functions
  public :: tridiag_init
  public :: tridiag_cleanup

  public :: tridiag_general
  public :: tridiag_constant

  ! Private module variables

  ! Temporary storage for Gaussian elimination coefficients
  complex(dp), allocatable :: mat_coeff(:), vec_coeff(:)

contains

  ! Module initialization
  !
  ! n :: dimension of linear system
  subroutine tridiag_init(n)

    integer :: n

    allocate(mat_coeff(n - 1))
    allocate(vec_coeff(n))

  end subroutine tridiag_init


  ! Module cleanup
  subroutine tridiag_cleanup()

    deallocate(mat_coeff)
    deallocate(vec_coeff)

  end subroutine tridiag_cleanup


  ! Thomas algorithm back-substitution to obtain solution from Gaussian
  ! elimination coefficients stored in `vec_coeff(:)` and `mat_coeff(:)`
  !
  ! res(:) :: x (unknown), where A x = b for the linear system
  subroutine tridiag_backsweep(res)

    complex(dp), intent(inout) :: res(:)
    integer :: i, n

    n = size(res)
    res(n) = vec_coeff(n)

    do i = n - 1, 1, -1
       res(i) = vec_coeff(i) - mat_coeff(i) * res(i + 1)
    end do

  end subroutine tridiag_backsweep


  ! Solve an n-dimensional tridiagional matrix equation A x = b.
  !
  ! diag(:) :: main diagonal entries
  ! u_diag(:) :: upper diagonal entries
  ! l_diag(:) :: lower diagonal entries
  ! vec(:) :: b (known), where A x = b for the linear system
  ! res(:) :: x (unknown), where A x = b for the linear system
  subroutine tridiag_general(diag, u_diag, l_diag, vec, res)

    complex(dp), intent(in) :: diag(:), u_diag(:), l_diag(:), vec(:)
    complex(dp), intent(out) :: res(:)

    integer :: i, n

    n = size(res)

    ! Forward sweep to calculate Gaussian elimination coefficients
    mat_coeff(1) = u_diag(1) / diag(1)
    vec_coeff(1) = vec(1) / diag(1)

    do i = 2, n - 1
       mat_coeff(i) = u_diag(i) / (diag(i) - l_diag(i) * mat_coeff(i - 1))
       vec_coeff(i) = (vec(i) - l_diag(i) * vec_coeff(i - 1)) / &
            (diag(i) - l_diag(i) * mat_coeff(i - 1))
    end do

    vec_coeff(n) = (vec(n) - l_diag(n) * vec_coeff(n - 1)) / &
         (diag(n) - l_diag(n) * mat_coeff(n - 1))

    ! Backward sweep to solve the system
    call tridiag_backsweep(res)

  end subroutine tridiag_general


  ! Solves an n-dimensional tridiagional matrix equation A x = b, where the
  ! band elements are constant.
  !
  ! diag_cnst :: main diagonal entry
  ! u_diag_cnst :: upper diagonal entry
  ! l_diag_cnst :: lower diagonal entry
  ! vec(:) :: b (known), where A x = b for the linear system
  ! res(:) :: x (unknown), where A x = b for the linear system
  subroutine tridiag_constant(diag_cnst, u_diag_cnst, l_diag_cnst, vec, res)

    complex(dp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst, vec(:)
    complex(dp), intent(out) :: res(:)

    integer :: i, n

    n = size(res)

    ! Forward sweep to calculate Gaussian elimination coefficients
    mat_coeff(1) = u_diag_cnst / diag_cnst
    vec_coeff(1) = vec(1) / diag_cnst

    do i = 2, n - 1
       mat_coeff(i) = u_diag_cnst / (diag_cnst - l_diag_cnst * mat_coeff(i - 1))
       vec_coeff(i) = (vec(i) - l_diag_cnst * vec_coeff(i - 1)) / &
            (diag_cnst - l_diag_cnst * mat_coeff(i - 1))
    end do

    vec_coeff(n) = (vec(n) - l_diag_cnst * vec_coeff(n - 1)) / &
         (diag_cnst - l_diag_cnst * mat_coeff(n - 1))

    ! Backward sweep to solve the system
    call tridiag_backsweep(res)

  end subroutine tridiag_constant

end module tridiag
