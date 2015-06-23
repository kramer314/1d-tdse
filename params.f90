! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

! Run parameters
module params
  use progvars

  implicit none

  private
  public :: params_init
  public :: params_psi0
  public :: params_pot

  ! Problem specific variables
  real(dp) :: x_0, d_x, k_0

contains

  ! Assign parameters
  subroutine params_init()
    implicit none

    ! Grid parameters
    x_min = -10_dp
    x_max = 10_dp
    n_x = 1e3
    dx = (x_max - x_min) / n_x

    ! Propagation parameters
    t_min = 0.0_dp
    t_max = 10.0_dp
    n_t = 1e3
    dt = (t_max - t_min) / n_t

    ! psi0 parameters
    x_0 = 2.0_dp
    d_x = 1.0_dp
    k_0 = 10.0_dp

    ! output parameters
    output_dir = "./output/"
    x_range_fname = "xrange.dat"
    t_range_fname = "trange.dat"
    psi0_fname = "psi0.dat"
    pot_fname = "pot.dat"
    wfunc_math_fname = "wfunc_math.dat"
    psi_xt_fname = "psi_xt.dat"

    output_grids = .true.
    output_psi0 = .true.
    output_pot = .true.
    output_wfunc_math = .true.
    output_psi_xt = .true.

    ! print_filter
    print_mod_t = 10
    print_mod_x = 10
  end subroutine params_init

  ! Initial wavefunction
  complex(dp) function params_psi0(x) result(val)
    implicit none

    real(dp), intent(in) :: x
    real(dp) :: norm

    norm = (pi * d_x**2) ** (-1.0_dp / 4.0_dp)
    val = norm * exp(-0.5_dp * ((x - x_0) / d_x)**2) * exp(j * k_0 * x)

  end function params_psi0

  ! Potential energy function (time independent)
  complex(dp) function params_pot(x) result(val)
    implicit none
    real(dp), intent(in) :: x

    val = 0.0_dp

    val = 10_dp * x**2

  end function params_pot

end module params
