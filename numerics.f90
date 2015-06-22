! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

! Various numerical functions
module numerics
  use progvars

  implicit none

  private

  public :: linspace
  public :: norm
  public :: expec_x

contains

  ! Populate an array with linearly-spaced values between x_min and x_max
  subroutine linspace(x_min, x_max, x_arr, dx)
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

  end subroutine linspace

  ! Calculate <psi | psi> = ||psi||^2
  real(dp) function norm(psi_arr) result(val)
    implicit none

    complex(dp), intent(in) :: psi_arr(:)
    integer(dp) :: i_x

    val = 0.0_dp
    do i_x = 1, n_x
       val = val + abs(psi_arr(i_x))**2 * dx
    end do
  end function norm

  ! Calculate <x>
  real(dp) function expec_x(psi_arr) result(val)
    implicit none

    complex(dp), intent(in) :: psi_arr(:)
    real(dp) :: x
    integer(dp) :: i_x
    
    val = 0.0_dp
    do i_x = 1, n_x
       x = x_range(i_x)
       val = val + abs(psi_arr(i_x))**2 * x * dx
    end do
    
  end function expec_x

end module numerics
