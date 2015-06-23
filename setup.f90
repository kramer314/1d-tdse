! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

! Initialization module
module setup
  use progvars
  use params, only: params_init, params_psi0, params_pot
  use numerics, only: numerics_linspace
  use output, only: output_init, output_cleanup
  use propagate, only: propagate_init, propagate_cleanup
  use tridiag, only: tridiag_init, tridiag_cleanup

  implicit none

  private
  public :: setup_init
  public :: setup_cleanup

contains
  ! Public pre-execution routine
  subroutine setup_init()
    implicit none

    call params_init
    call output_init
    call allocate_arrays
    call setup_grids
    call setup_psi
    call setup_potential
    call propagate_init
    call tridiag_init(n_x, tridiag_mat_coeff, tridiag_vec_coeff)
  end subroutine setup_init

  ! Public post-execution routine
  subroutine setup_cleanup()
    implicit none

    call output_cleanup
    call deallocate_arrays
    call propagate_cleanup
    call tridiag_cleanup(tridiag_mat_coeff, tridiag_vec_coeff)

  end subroutine setup_cleanup

  ! Allocate arrays
  subroutine allocate_arrays()
    implicit none

    allocate(x_range(n_x))
    allocate(psi_arr(n_x))
    allocate(pot_arr(n_x))
    allocate(phi_arr(n_x))
    allocate(t_range(n_t))

  end subroutine allocate_arrays

  ! Deallocate arrays
  subroutine deallocate_arrays()
    implicit none

    deallocate(x_range)
    deallocate(t_range)
    deallocate(psi_arr)
    deallocate(pot_arr)

  end subroutine deallocate_arrays

  ! Setup numerical grids
  subroutine setup_grids()
    implicit none

    call numerics_linspace(x_min, x_max, x_range, dx)
    call numerics_linspace(t_min, t_max, t_range, dt)

  end subroutine setup_grids

  ! Initialize psi(x, t = 0)
  subroutine setup_psi()
    implicit none

    real(dp) :: x
    integer(dp) :: i_x

    do i_x = 1, n_x
       x = x_range(i_x)
       psi_arr(i_x) = params_psi0(x)
    end do

  end subroutine setup_psi

  ! Initialize potential V(x)
  subroutine setup_potential()
    implicit none

    real(dp) :: x
    integer(dp) :: i_x

    do i_x = 1, n_x
       x = x_range(i_x)
       pot_arr(i_x) = params_pot(x)
    end do

  end subroutine setup_potential

end module setup
