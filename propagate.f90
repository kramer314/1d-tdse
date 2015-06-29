! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

! Wavefunction propagation
module propagate
  use progvars
  use tridiag, only: tridiag_cnst
  use params, only: params_pot

  implicit none

  complex(dp), allocatable :: exp_pot_arr(:)
  complex(dp) :: sym_cnst, diag_cnst

  private
  public :: propagate_init
  public :: propagate_cleanup
  public :: propagate_cn_splitop

contains

  ! Initialize propagation variables / arrays
  subroutine propagate_init()
    implicit none

    allocate(exp_pot_arr(n_x))

    sym_cnst = - (j * dt) / (8.0_dp * dx**2)
    diag_cnst = (0.5_dp - 2.0_dp * sym_cnst)
    exp_pot_arr(:) = exp(-j * pot_arr(:) * dt)

  end subroutine propagate_init

  subroutine propagate_calc_pot(i_t)
    implicit none

    integer(dp), intent(in) :: i_t

    real(dp) :: x, t
    integer(dp) :: i_x

    t = t_range(i_t)

    do i_x = 1, n_x
       x = x_range(i_x)
       pot_arr(i_x) = params_pot(x, t)
    end do

    exp_pot_arr(:) = exp(-j * pot_arr(:) * dt)
    
  end subroutine propagate_calc_pot

  ! Cleanup propagation arrays
  subroutine propagate_cleanup()
    implicit none

    deallocate(exp_pot_arr)
  end subroutine propagate_cleanup

  ! Crank-Nicolson propagation
  subroutine propagate_cn(psi_arr)
    implicit none

    complex(dp), intent(inout) :: psi_arr(:)

    ! Solve for auxillary wavefunction, then propagate
    call tridiag_cnst(diag_cnst, sym_cnst, sym_cnst, psi_arr, phi_arr, &
         tridiag_mat_coeff, tridiag_vec_coeff)
    psi_arr(:) = phi_arr(:) - psi_arr(:)

  end subroutine propagate_cn

  ! One-dimensional Crank-Nicolson split-operator propagation
  subroutine propagate_cn_splitop(psi_arr, i_t)
    implicit none

    complex(dp), intent(inout) :: psi_arr(:)
    integer(dp), intent(in) :: i_t

    call propagate_calc_pot(i_t)
    
    psi_arr(:) = exp_pot_arr(:) * psi_arr(:)
    
    call propagate_cn(psi_arr)

  end subroutine propagate_cn_splitop

end module propagate
