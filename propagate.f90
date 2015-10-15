! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

! Wavefunction propagation using a three-point finite difference
! Crank-Nicolson scheme
module propagate
  use progvars
  use tridiag, only: tridiag_constant
  use params, only: params_pot

  implicit none

  private

  ! Generic module functions
  public :: propagate_init
  public :: propagate_cleanup

  public :: propagate_cn_splitop

  ! Private module variables

  complex(dp), allocatable :: exp_pot_arr(:)
  complex(dp), allocatable :: old_pot_arr(:)

  ! auxillary wavefunction array
  complex(dp), allocatable :: phi_arr(:)
  ! tridiagonal matrix elements for auxillary wavefunction calculation
  complex(dp) :: sym_cnst, diag_cnst

contains

  ! Module initialization
  subroutine propagate_init()

    allocate(phi_arr(n_x))
    allocate(exp_pot_arr(n_x))

    allocate(old_pot_arr(n_x))

    ! Ensures the potential is calculated correctly for the first timestep
    old_pot_arr = 1.0_dp

    sym_cnst = - (j * hbar**2 * dt) / (8.0_dp * m * dx**2)
    diag_cnst = (0.5_dp - 2.0_dp * sym_cnst)

  end subroutine propagate_init


  ! Module cleanup
  subroutine propagate_cleanup()

    deallocate(phi_arr)
    deallocate(exp_pot_arr)
    deallocate(old_pot_arr)

  end subroutine propagate_cleanup


  ! Setup time-dependent potential propagation
  !
  ! i_t :: time index
  subroutine propagate_calc_pot(i_t)

    integer, intent(in) :: i_t

    real(dp) :: x, t
    integer :: i_x
    complex(dp) :: pot_xt

    t = t_range(i_t)

    do i_x = 1, n_x
       x = x_range(i_x)

       pot_xt = params_pot(x, t)

       ! Calculating exp(-j V_xt dt) is expensive; we should only update it
       ! as needed.
       if (abs((pot_xt - old_pot_arr(i_x))) .ge. eps_dp) then
          exp_pot_arr(i_x) = exp(-j * pot_xt * dt)
          old_pot_arr(i_x) = pot_xt
       end if

    end do

  end subroutine propagate_calc_pot


  ! Crank-Nicolson propagation for the kinetic energy portion of the
  ! Hamiltonian
  !
  ! psi_arr :: ket wavefunction array
  subroutine propagate_cn_ke(psi_arr)

    complex(dp), intent(inout) :: psi_arr(:)

    ! Solve for auxillary wavefunction, then propagate
    call tridiag_constant(diag_cnst, sym_cnst, sym_cnst, psi_arr, phi_arr)
    psi_arr(:) = phi_arr(:) - psi_arr(:)

  end subroutine propagate_cn_ke


  ! One-dimensional Crank-Nicolson split-operator propagation
  !
  ! psi_arr :: ket wavefunction array
  ! i_t :: time index
  subroutine propagate_cn_splitop(psi_arr, i_t)

    complex(dp), intent(inout) :: psi_arr(:)
    integer, intent(in) :: i_t

    call propagate_calc_pot(i_t)
    psi_arr(:) = exp_pot_arr(:) * psi_arr(:)

    call propagate_cn_ke(psi_arr)

  end subroutine propagate_cn_splitop

end module propagate
