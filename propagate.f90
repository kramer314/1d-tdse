! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

! Wavefunction propagation
module propagate
  use progvars
  use tridiag, only: tridiag_sym_cnst

  implicit none

  complex(dp), allocatable :: exp_pot_arr(:), diag_arr(:)
  complex(dp) :: sym_cnst

  private
  public :: propagate_init
  public :: propagate_cleanup
  public :: propagate_cn_splitop

contains

  ! Initialize propagation variables / arrays
  subroutine propagate_init()
    implicit none

    allocate(diag_arr(n_x))
    allocate(exp_pot_arr(n_x))

    sym_cnst = - (j * dt) / (8.0_dp * dx**2)
    diag_arr(:) = (0.5_dp - 2.0_dp * sym_cnst)
    exp_pot_arr(:) = exp(-j * pot_arr(:) * dt)

  end subroutine propagate_init

  ! Cleanup propagation arrays
  subroutine propagate_cleanup()
    implicit none

    deallocate(diag_arr)
    deallocate(exp_pot_arr)
  end subroutine propagate_cleanup

  ! Crank-Nicolson propagation
  subroutine propagate_cn(psi_arr)
    implicit none

    complex(dp), intent(inout) :: psi_arr(:)

    ! Solve for auxillary wavefunction, then propagate
    call tridiag_sym_cnst(diag_arr, sym_cnst, psi_arr, phi_arr, &
         tridiag_mat_coeff, tridiag_vec_coeff)
    psi_arr(:) = phi_arr(:) - psi_arr(:)

  end subroutine propagate_cn

  ! One-dimensional Crank-Nicolson split-operator propagation
  subroutine propagate_cn_splitop(psi_arr)
    implicit none

    complex(dp), intent(inout) :: psi_arr(:)

    psi_arr(:) = exp_pot_arr(:) * psi_arr(:)
    call propagate_cn(psi_arr)

  end subroutine propagate_cn_splitop

end module propagate
