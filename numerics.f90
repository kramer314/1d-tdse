! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

! Assorted numerical helper functions
! Note: This is a program-independent module
module numerics
  use globvars, only: dp

  implicit none

  private

  public :: numerics_linspace

contains

  ! Populate an array with linearly-spaced values between x_min and x_max
  !
  ! x_min :: minimum array value
  ! x_max :: maximum array value
  ! x_arr(:) :: array to populate, must already be allocated
  ! dx :: step size of `x_arr(:)`, given its size and `x_min`, `x_max`
  subroutine numerics_linspace(x_min, x_max, x_arr, dx)

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
