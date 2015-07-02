! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

! One-dimensional wavefunction math
! Note: This is a program-independent module
module wfmath
  use globvars, only: dp

  implicit none

  private

  ! Generic module functions
  public :: wfmath_init
  public :: wfmath_cleanup

  ! Inner products and expectation values
  public :: wfmath_iprod
  public :: wfmath_norm
  public :: wfmath_expec_op
  public :: wfmath_expec_op2
  public :: wfmath_stdev_op

  ! Wrapper / helper functions
  public :: wfmath_expec_x
  public :: wfmath_stdev_x

  ! Other functions
  public :: wfmath_autocorr

  ! Private module variables

  complex(dp), allocatable :: work_arr(:)
  ! Internal copy of spatial grid parameters for modularity
  real(dp), allocatable :: x_arr(:)
  real(dp) :: dx

contains

  ! Module initialization
  !
  ! x_range :: spatial numerical grid
  ! d_x :: grid spacing
  subroutine wfmath_init(x_range, d_x)

    real(dp), intent(in) :: x_range(:), d_x
    integer(dp) :: n_x

    n_x = size(x_range)
    dx = d_x

    allocate(x_arr(n_x))
    x_arr(:) = x_range(:)

    allocate(work_arr(n_x))

  end subroutine wfmath_init

  ! Module cleanup
  subroutine wfmath_cleanup()

    deallocate(x_arr)
    deallocate(work_arr)

  end subroutine wfmath_cleanup

  ! Calculate inner product <psi_1 | func | psi_2>
  !
  ! psi1_arr :: bra array
  ! func_arr :: operator array
  ! psi2_arr :: ket array
  complex(dp) function wfmath_iprod(psi1_arr, func_arr, psi2_arr) result(val)
    implicit none

    complex(dp), intent(in) :: psi1_arr(:), psi2_arr(:), func_arr(:)

    val = sum(conjg(psi1_arr(:)) * func_arr(:) * psi2_arr(:)) * dx

  end function wfmath_iprod

  ! Calculate wavefunction norm ||psi||^2 = <psi | psi>
  !
  ! psi_arr :: ket wavefunction array
  real(dp) function wfmath_norm(psi_arr) result(val)

    complex(dp), intent(in) :: psi_arr(:)

    work_arr(:) = 1.0_dp
    val = real(wfmath_iprod(psi_arr, work_arr, psi_arr))

  end function wfmath_norm

  ! Calculate operator expectation value <op> = <psi | op | psi>
  ! We assume that the operator is Hermitian, hence a real result.
  !
  ! psi_arr :: ket wavefunction array
  ! op_arr :: operator array
  real(dp) function wfmath_expec_op(psi_arr, op_arr) result(val)

    complex(dp), intent(in) :: psi_arr(:), op_arr(:)

    val = real(wfmath_iprod(psi_arr, op_arr, psi_arr))

  end function wfmath_expec_op

  ! Calculate operator squared expectation value <op^2> = <psi | op^2 | psi>
  ! We assume that the operator is Hermitian, hence a real result.
  !
  ! psi_arr :: ket wavefunction array
  ! op_arr :: operator array
  real(dp) function wfmath_expec_op2(psi_arr, op_arr) result(val)

    complex(dp), intent(in) :: psi_arr(:)
    complex(dp), intent(in) :: op_arr(:)

    work_arr(:) = op_arr(:) * op_arr(:)
    val = real(wfmath_iprod(psi_arr, work_arr, psi_arr))

  end function wfmath_expec_op2

  ! Calculate operator standard deviation stdev(op) = sqrt(<op^2> - <op>^2)
  ! We assume that the operator is Hermitian, hence a real result.
  !
  ! psi_arr :: ket wavefunction array
  ! op_arr :: operator array
  real(dp) function wfmath_stdev_op(psi_arr, op_arr) result(val)

    complex(dp), intent(in) :: psi_arr(:)
    complex(dp), intent(in) :: op_arr(:)
    real(dp) :: expec_op, expec_op2

    expec_op = wfmath_expec_op(psi_arr, op_arr)
    expec_op2 = wfmath_expec_op2(psi_arr, op_arr)
    val = sqrt(expec_op2 - expec_op**2)

  end function wfmath_stdev_op

  ! Calculate <x> = <psi | x | psi>
  !
  ! psi_arr :: ket wavefunction array
  real(dp) function wfmath_expec_x(psi_arr) result(val)

    complex(dp), intent(in) :: psi_arr(:)
    val = wfmath_expec_op(psi_arr, cmplx(x_arr, kind=dp))

  end function wfmath_expec_x

  ! Calculate stdev(x) = sqrt(<x^2> - <x>^2)
  !
  ! psi_arr :: ket wavefunction array
  real(dp) function wfmath_stdev_x(psi_arr) result(val)

    complex(dp), intent(in) :: psi_arr(:)
    val = wfmath_stdev_op(psi_arr, cmplx(x_arr, kind=dp))

  end function wfmath_stdev_x

  ! Calculate autocorrelation function C(t) = <psi(x,t) | psi(x,0)>
  !
  ! psi_arr :: ket wavefunction array psi(x,t)
  ! psi0_arr :: initial ket wavefunction array psi(x,0)
  complex(dp) function wfmath_autocorr(psi_arr, psi0_arr) result(val)

    complex(dp), intent(in) :: psi_arr(:), psi0_arr(:)

    work_arr(:) = 1.0_dp
    val = wfmath_iprod(psi_arr, work_arr, psi0_arr)

  end function wfmath_autocorr

end module wfmath
