! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

module output
  use progvars
  use files, only: files_ensure_dir
  use numerics, only: numerics_expec_x, numerics_norm, numerics_stdev_x

  implicit none

  private
  public :: output_init
  public :: output_cleanup
  public :: output_time_indep
  public :: output_time_dep
contains

  ! Initialize output
  subroutine output_init()
    implicit none

    call files_ensure_dir(output_dir)

    open(unit=1, file=trim(output_dir)//trim(psi_xt_fname))

    if (output_grids) then
       open(unit=2, file=trim(output_dir)//trim(x_range_fname))
       open(unit=3, file=trim(output_dir)//trim(t_range_fname))
    end if

    if (output_psi0) then
       open(unit=4, file=trim(output_dir)//trim(psi0_fname))
    end if

    if (output_pot) then
       open(unit=5, file=trim(output_dir)//trim(pot_fname))
    end if

    if (output_wfunc_math) then
       open(unit=6, file=trim(output_dir)//trim(wfunc_math_fname))
    end if

  end subroutine output_init

  ! Cleanup output
  subroutine output_cleanup()
    implicit none

    close(unit=1)
    close(unit=2)
    close(unit=3)
    close(unit=4)
    close(unit=5)
    close(unit=6)

  end subroutine output_cleanup

  ! Time independent output
  subroutine output_time_indep()
    implicit none

    integer(dp) :: i_x

    do i_x = 1, n_x
       if (output_grids) then
          write(2, dp_format) x_range(i_x)
          write(3, dp_format) t_range(i_x)
       end if
       if (output_psi0) then
          write(4, dp_format) abs(psi_arr(i_x))**2
       end if

       if (output_pot) then
          write(5, dp_format) pot_arr(i_x)
       end if
    end do

  end subroutine output_time_indep

  ! Time-dependent output
  subroutine output_time_dep(psi_arr)
    implicit none

    complex(dp), intent(in) :: psi_arr(:)
    integer(dp) :: i_x
    real(dp) :: norm, e_x, stdev_x

    do i_x = 1, n_x
       write(1, dp_format, advance="no") abs(psi_arr(i_x))**2
    end do
    write(1, *)

    if (output_wfunc_math) then
       norm = numerics_norm(psi_arr)
       e_x = numerics_expec_x(psi_arr)
       stdev_x = numerics_stdev_x(psi_arr)

       write(6, "(3"//dp_format_raw//")") norm, e_x, stdev_x
    end if

  end subroutine output_time_dep

end module output
