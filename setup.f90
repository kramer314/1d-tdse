! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

! Initialization module
module setup
  use progvars
  use params, only: assign_params, psi0, pot
  use numerics, only: linspace
  use propagate, only: init_propagate

  implicit none

  private
  public :: init, cleanup

contains
  ! Public pre-execution routine
  subroutine init()
    implicit none
    call assign_params
    call allocate_arrays
    call init_grids
    call init_psi
    call init_pot
    call init_propagate
  end subroutine init

  ! Public post-execution routine
  subroutine cleanup()
    implicit none
    deallocate(x_range, psi_arr, phi_arr, pot_arr, t_range)
  end subroutine cleanup
  
  ! Allocate arrays
  subroutine allocate_arrays()
    implicit none
    allocate(x_range(n_x))
    allocate(psi_arr(n_x))
    allocate(pot_arr(n_x))
    allocate(phi_arr(n_x))
    allocate(t_range(n_t))
  end subroutine allocate_arrays

  ! Setup numerical grids
  subroutine init_grids
    implicit none

    call linspace(x_min, x_max, x_range, dx)
    call linspace(t_min, t_max, t_range, dt)
  end subroutine init_grids

  ! Initialize psi(x, t = 0)
  subroutine init_psi
    implicit none

    real(dp) :: x
    integer(dp) :: i_x

    do i_x = 1, n_x
       x = x_range(i_x)
       psi_arr(i_x) = psi0(x)
    end do

  end subroutine init_psi

  ! Initialize potential V(x)
  subroutine init_pot
    implicit none

    real(dp) :: x
    integer(dp) :: i_x

    do i_x = 1, n_x
       x = x_range(i_x)
       pot_arr(i_x) = pot(x)
    end do
  end subroutine init_pot
end module setup
