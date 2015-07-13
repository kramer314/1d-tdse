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
  public :: numerics_d2

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
  ! respect to a real-valued coordinate with constant spacing:
  !
  !   f'(x) ~= (f(x + dx) - f(x - dx)) / (2 dx)
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

    ! Special edge cases - we assume the grid is fine enough and the function
    ! is smooth enough that this is accurate enough. For the left case:
    !
    !   f(x + h) = arr(2), f(x - h) = arr(1)
    !
    ! For the right case:
    !
    !   f(x + h) = arr(n), f(x - h) = arr(n - 1)
    d_arr(1) = arr(2) - arr(1)
    do i_x = 2, n_x - 1
       d_arr(i_x) = arr(i_x + 1) - arr(i_x - 1)
    end do
    d_arr(n_x) = arr(n_x) - arr(n_x - 1)

    d_arr(:) = scale * d_arr(:)

  end subroutine numerics_d1


  ! Calculate the 2nd derivative of a complex array-valued function with
  ! respect to a real-valued coordinate with constant spacing:
  !
  !   f''(x) ~= (f(x - dx) - 2 f(x) + f(x + dx)) / (dx**2)
  !
  ! arr :: function array
  ! d2_arr :: array to populate
  ! dx :: coordinate spacing
  subroutine numerics_d2(arr, d2_arr, dx)

    complex(dp), intent(in) :: arr(:)
    complex(dp), intent(inout) :: d2_arr(:)
    real(dp), intent(in) :: dx
    integer(dp) :: i_x, n_x
    real(dp) :: scale

    n_x = size(arr)
    scale = 1.0_dp / dx**2

    ! Special edge cases - we assume the grid is fine enough and the function
    ! is smooth enough that this is accurate enough. For the left case:
    !
    !   f(x + h) = arr(2), f(x - h) = arr(1)
    !
    ! For the right case:
    !
    !   f(x + h) = arr(n), f(x - h) = arr(n - 1)
    d2_arr(1) = arr(1) + arr(2)
    do i_x = 2, n_x - 1
       d2_arr(i_x) = arr(i_x - 1) - 2 * arr(i_x) + arr(i_x + 1)
    end do
    d2_arr(n_x) = arr(n_x) + arr(n_x - 1)

    d2_arr(:) = scale * d2_arr(:)

  end subroutine numerics_d2

end module numerics
