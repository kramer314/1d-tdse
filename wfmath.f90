! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

! Wavefunction math
module wfmath
  use progvars

  implicit none

  private

  complex(dp), allocatable :: work_arr_x(:)

  public :: wfmath_init
  public :: wfmath_cleanup

  public :: wfmath_iprod_x
  public :: wfmath_norm_x
  public :: wfmath_expec_x
  public :: wfmath_stdev_x
  public :: wfmath_autocorr_x

contains

  ! Setup wfmath work array
  subroutine wfmath_init()
    implicit none

    allocate(work_arr_x(n_x))
  end subroutine wfmath_init

  ! Cleanup wfmath work array
  subroutine wfmath_cleanup()
    implicit none

    deallocate(work_arr_x)
  end subroutine wfmath_cleanup

  ! Inner product <psi_1 | func | psi_2>
  complex(dp) function wfmath_iprod_x(psi1_arr, func_arr, psi2_arr) result(val)
    implicit none

    complex(dp), intent(in) :: psi1_arr(:), psi2_arr(:), func_arr(:)

    val = sum(conjg(psi1_arr(:)) * func_arr(:) * psi2_arr(:)) * dx

  end function wfmath_iprod_x

  ! ||psi||^2 = <psi | psi>
  real(dp) function wfmath_norm_x(psi_arr) result(val)
    implicit none

    complex(dp), intent(in) :: psi_arr(:)

    work_arr_x(:) = 1.0_dp
    val = real(wfmath_iprod_x(psi_arr, work_arr_x, psi_arr))

  end function wfmath_norm_x

  ! <x> = <psi | x | psi>
  real(dp) function wfmath_expec_x(psi_arr) result(val)
    implicit none

    complex(dp), intent(in) :: psi_arr(:)

    work_arr_x(:) = x_range(:)
    val = real(wfmath_iprod_x(psi_arr, work_arr_x, psi_arr))

  end function wfmath_expec_x

  ! <x^2> = <psi | x^2 | psi>
  real(dp) function wfmath_expec_x2(psi_arr) result(val)
    implicit none

    complex(dp), intent(in) :: psi_arr(:)

    work_arr_x(:) = x_range(:) * x_range(:)
    val = real(wfmath_iprod_x(psi_arr, work_arr_x, psi_arr))

  end function wfmath_expec_x2

  ! stdev(x) = sqrt(<x^2> - <x>^2)
  real function wfmath_stdev_x(psi_arr) result(val)
    implicit none

    complex(dp), intent(in) :: psi_arr(:)

    val = sqrt(wfmath_expec_x2(psi_arr) - wfmath_expec_x(psi_arr)**2)

  end function wfmath_stdev_x

  ! C(t) = <psi(x,t) | psi(x,0)>
  complex(dp) function wfmath_autocorr_x(psi_arr) result(val)
    implicit none

    complex(dp), intent(in) :: psi_arr(:)

    work_arr_x(:) = 1.0_dp
    val = wfmath_iprod_x(psi_arr, work_arr_x, psi0_arr)

  end function wfmath_autocorr_x

end module wfmath
