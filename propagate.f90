! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

! Wavefunction propagation
module propagate
  use progvars
  use tridiag, only: tridiag_sym_cnst

  implicit none

  private
  public :: cn_splitop
  public :: init_propagate
  
contains

  subroutine init_propagate()
    
    implicit none
    allocate(diag_arr(n_x))
    allocate(exp_pot_arr(n_x))

    sym_cnst = - (j * dt) / (8.0_dp * dx**2)
    diag_arr(:) = (0.5_dp - 2.0_dp * sym_cnst)
    exp_pot_arr(:) = exp(-j * pot_arr(:) * dt)
  end subroutine init_propagate
  
  ! Crank-Nicholson propagation
  subroutine cn_propagate(psi_arr)
    implicit none
    complex(dp), intent(inout) :: psi_arr(:)

    ! Solve for auxillary wavefunction
    call tridiag_sym_cnst(diag_arr, sym_cnst, psi_arr, phi_arr)
    ! Propagated wavefunction
    psi_arr(:) = phi_arr(:) - psi_arr(:)

  end subroutine cn_propagate

  subroutine cn_splitop(psi_arr)
    implicit none
    complex(dp), intent(inout) :: psi_arr(:)

    psi_arr(:) = exp_pot_arr(:) * psi_arr(:)
    call cn_propagate(psi_arr)
    
  end subroutine cn_splitop
end module propagate
