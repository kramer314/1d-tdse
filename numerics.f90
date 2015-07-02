! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

! Assorted numerical helper functions
! Note: This is a program-independent module
module numerics
  use globvars, only: dp

  implicit none

  private

  public :: numerics_linspace
  public :: numerics_d1

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

  ! Calculate the 1st derivative of a complex array-valued function with
  ! respect to a real-valued coordinate.
  !
  ! arr :: function array
  ! d_arr :: array to populate
  ! dx :: coordinate spacing
  subroutine numerics_d1(arr, d_arr, dx)

    complex(dp), intent(in) :: arr(:)
    complex(dp), intent(inout) :: d_arr(:)
    real(dp), intent(in) :: dx
    integer(dp) :: i_x, n_x
    real(dp) :: scale

    n_x = size(arr)
    scale = 1.0_dp / (2.0_dp * dx)

    d_arr(1) = (arr(2) - arr(1))
    do i_x = 2, n_x - 1
       d_arr(i_x) = (arr(i_x + 1) - arr(i_x - 1))
    end do
    d_arr(n_x) = (arr(n_x) - arr(n_x - 1))

    d_arr(:) = scale * d_arr(:)

  end subroutine numerics_d1

end module numerics
