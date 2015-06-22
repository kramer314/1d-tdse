! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution

! Global program variables
module progvars
  implicit none

  ! Double precision type / formatting
  integer, parameter :: dp = kind(0.d0)
  character(*), parameter :: dp_format = "(es25.16e3)"

  ! Numerical constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: e = exp(1.0_dp)
  complex(dp), parameter :: j = (0.0_dp, 1.0_dp)

  ! Grid parameters
  real(dp) :: x_min, x_max, dx, t_min, t_max, dt
  integer(dp) :: n_x, n_t

  ! Arrays
  real(dp), allocatable :: x_range(:), t_range(:)
  complex(dp), allocatable :: psi_arr(:), phi_arr(:), pot_arr(:)

end module progvars
