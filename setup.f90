! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

! Program initialization / cleanup module
module setup
  use progvars
  use params, only: params_init, params_psi0
  use numerics, only: numerics_linspace
  use output, only: output_init, output_cleanup
  use propagate, only: propagate_init, propagate_cleanup
  use tridiag, only: tridiag_init, tridiag_cleanup
  use wfmath, only: wfmath_init, wfmath_cleanup

  implicit none

  private

  ! Generic module functions
  public :: setup_init
  public :: setup_cleanup

contains

  ! Module initialization
  subroutine setup_init()

    ! Initialize parameters
    call params_init()

    ! Allocate main arrays
    allocate(x_range(n_x))
    allocate(psi0_arr(n_x))
    allocate(psi_arr(n_x))
    allocate(flux_arr(n_x))
    allocate(t_range(n_t))

    ! Setup numerical grids
    call numerics_linspace(x_min, x_max, x_range, dx)
    call numerics_linspace(t_min, t_max, t_range, dt)

    ! Setup initial wavefunction array
    call setup_calc_psi0()

    ! Initialize program-independent modules
    call wfmath_init(x_range, dx, m, hbar)
    call tridiag_init(n_x)

    ! Initialize program-dependent modules
    call propagate_init()
    call output_init()

  end subroutine setup_init


  ! Module cleanup
  subroutine setup_cleanup()

    ! Deallocate main arrays
    deallocate(x_range)
    deallocate(t_range)
    deallocate(psi0_arr)
    deallocate(psi_arr)
    deallocate(flux_arr)

    ! Cleanup program-independent modules
    call wfmath_cleanup
    call tridiag_cleanup

    ! Cleanup program-dependent modules
    call output_cleanup
    call propagate_cleanup

  end subroutine setup_cleanup


  ! Calculate psi(x, t=0) and populate relevant arrays
  subroutine setup_calc_psi0()

    real(dp) :: x
    integer :: i_x

    do i_x = 1, n_x
       x = x_range(i_x)
       psi0_arr(i_x) = params_psi0(x)
    end do
    psi_arr(:) = psi0_arr(:)

  end subroutine setup_calc_psi0

end module setup
