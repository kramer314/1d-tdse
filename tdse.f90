! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

program tdse
  use progvars
  use setup, only: setup_init, setup_cleanup
  use propagate, only: propagate_cn_splitop
  use output, only: output_time_dep, output_time_indep

  implicit none

  integer(dp) :: i_t

  call setup_init

  call output_time_indep

  call output_time_dep(psi_arr)
  do i_t = 1, n_t

     call propagate_cn_splitop(psi_arr)

     if (mod(i_t, print_mod_t) .eq. 0) then

        write(*,*) i_t, n_t
        call output_time_dep(psi_arr)

     end if

  end do

  call setup_cleanup

end program tdse
