! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution

! Program-specific variables
module progvars
  use globvars

  implicit none

  ! Grid parameters
  real(dp) :: x_min, x_max, dx, t_min, t_max, dt
  integer(dp) :: n_x, n_t

  ! Output parameters
  character(120) :: output_dir
  character(120) :: x_range_fname, t_range_fname, pot_fname, psi0_fname, &
       psi_xt_fname, wfunc_math_fname
  logical :: output_grids, output_psi0, output_psixt, output_wfunc_math, &
       output_psi_xt
  integer(dp) :: print_mod_t, print_mod_x

  ! Arrays
  real(dp), allocatable :: x_range(:), t_range(:)
  complex(dp), allocatable :: psi_arr(:), psi0_arr(:)

end module progvars
