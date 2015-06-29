! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

! Various numerical functions
module numerics
  use progvars

  implicit none

  private

  public :: numerics_linspace

contains

  ! Populate an array with linearly-spaced values between x_min and x_max
  subroutine numerics_linspace(x_min, x_max, x_arr, dx)
    implicit none

    real(dp), intent(in) :: x_min, x_max
    real(dp), intent(out) :: x_arr(:)
    real(dp), intent(out) :: dx

    integer(dp) :: i_x, n_x

    n_x = size(x_arr)
    dx = (x_max - x_min) / n_x

    do i_x = 1, n_x
       x_arr(i_x) = x_min + i_x * dx
    end do

  end subroutine numerics_linspace

end module numerics
