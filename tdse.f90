! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

program tdse
  use progvars
  use setup, only: setup_init, setup_cleanup
  use propagate, only: propagate_cn_splitop
  use numerics, only: norm, expec_x, stdev_x

  implicit none

  integer(dp) :: i_x, i_t

  call setup_init

  open(unit=1, file="psi0.dat")
  open(unit=2, file="psit.dat")
  open(unit=3, file="norm.dat")
  open(unit=4, file="expec_x.dat")
  open(unit=5, file="stdev_x.dat")
  open(unit=6, file="pot.dat")

  do i_x = 1, n_x
     write(1, dp_format) real(psi_arr(i_x))**2
     write(6, dp_format) real(pot_arr(i_x))
     write(2, dp_format, advance="no") abs(psi_arr(i_x))**2
  end do
  write(2,*)

  do i_t = 1, n_t
     call propagate_cn_splitop(psi_arr)
     write(3, dp_format) norm(psi_arr)
     write(4, dp_format) expec_x(psi_arr)
     write(5, dp_format) stdev_x(psi_arr)
     do i_x = 1, n_x
        write(2, dp_format, advance="no") abs(psi_arr(i_x))**2
     end do
     write(2,*)
  end do

  call setup_cleanup

  close(unit=1)
  close(unit=2)
  close(unit=3)
  close(unit=4)
  close(unit=5)
  close(unit=6)

end program tdse
