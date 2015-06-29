! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

! Various numerical functions
module numerics
  use progvars

  implicit none

  private

  public :: numerics_linspace
  public :: numerics_norm
  public :: numerics_expec_x
  public :: numerics_expec_x2
  public :: numerics_stdev_x
  public :: numerics_autocorr
  public :: numerics_autocorr2

contains

  ! Populate an array with linearly-spaced values between x_min and x_max
  subroutine numerics_linspace(x_min, x_max, x_arr, dx)
    implicit none

    real(dp), intent(in) :: x_min, x_max
    real(dp), intent(out) :: x_arr(:)
    real(dp), intent(out) :: dx

    integer :: i_x, n_x

    n_x = size(x_arr)
    dx = (x_max - x_min) / n_x

    do i_x = 1, n_x
       x_arr(i_x) = x_min + i_x * dx
    end do

  end subroutine numerics_linspace

  ! Calculate <psi | psi> = ||psi||^2
  real(dp) function numerics_norm(psi_arr) result(val)
    implicit none

    complex(dp), intent(in) :: psi_arr(:)
    integer(dp) :: i_x

    val = 0.0_dp
    do i_x = 1, n_x
       val = val + abs(psi_arr(i_x))**2 * dx
    end do
  end function numerics_norm

  ! Calculate <x>
  real(dp) function numerics_expec_x(psi_arr) result(val)
    implicit none

    complex(dp), intent(in) :: psi_arr(:)
    real(dp) :: x
    integer(dp) :: i_x

    val = 0.0_dp
    do i_x = 1, n_x
       x = x_range(i_x)
       val = val + abs(psi_arr(i_x))**2 * x * dx
    end do

  end function numerics_expec_x

  ! Calculate <x^2>
  real(dp) function numerics_expec_x2(psi_arr) result(val)
    implicit none

    complex(dp), intent(in) :: psi_arr(:)
    real(dp) :: x
    integer(dp) :: i_x

    val = 0.0_dp
    do i_x = 1, n_x
       x = x_range(i_x)
       val = val + abs(psi_arr(i_x))**2 * x**2 * dx
    end do

  end function numerics_expec_x2

  ! Calculate stdev(x)
  real(dp) function numerics_stdev_x(psi_arr) result(val)
    implicit none

    complex(dp), intent(in) :: psi_arr(:)

    ! todo: optimize so loops aren't repeated
    val = sqrt(numerics_expec_x2(psi_arr) - numerics_expec_x(psi_arr)**2)

  end function numerics_stdev_x

  ! Calculate autocorrelation
  complex(dp) function numerics_autocorr(psi_arr) result(val)
    implicit none

    complex(dp), intent(in) :: psi_arr(:)
    real(dp) :: x
    integer(dp) :: i_x

    val = 0.0_dp
    do i_x = 1, n_x
      x = x_range(i_x)
      val = val + conjg(psi_arr(i_x)) * psi0_arr(i_x) * dx
    end do
  end function numerics_autocorr

  ! Calculate norm of autocorrelation
  real(dp) function numerics_autocorr2(psi_arr) result(val)
    implicit none

    complex(dp), intent(in) :: psi_arr(:)
    val = abs(numerics_autocorr(psi_arr))**2
  end function numerics_autocorr2

end module numerics
