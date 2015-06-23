! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

module output
  use progvars
  use files, only: files_ensure_dir
  use numerics, only: numerics_expec_x, numerics_norm, numerics_stdev_x

  implicit none

  ! Output file unit numbers
  integer(dp), parameter :: psi_xt_unit = 99
  integer(dp), parameter :: x_range_unit = 98
  integer(dp), parameter :: t_range_unit = 98
  integer(dp), parameter :: psi0_unit = 96
  integer(dp), parameter :: pot_unit = 95
  integer(dp), parameter :: wfunc_math_unit = 94

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
    if (output_psi_xt) then
       open(unit=psi_xt_unit, file=trim(output_dir)//trim(psi_xt_fname))
    end if

    if (output_grids) then
       open(unit=x_range_unit, file=trim(output_dir)//trim(x_range_fname))
       open(unit=t_range_unit, file=trim(output_dir)//trim(t_range_fname))
    end if

    if (output_psi0) then
       open(unit=psi0_unit, file=trim(output_dir)//trim(psi0_fname))
    end if

    if (output_pot) then
       open(unit=pot_unit, file=trim(output_dir)//trim(pot_fname))
    end if

    if (output_wfunc_math) then
       open(unit=wfunc_math_unit, file=trim(output_dir)//trim(wfunc_math_fname))
    end if

  end subroutine output_init

  ! Cleanup output
  subroutine output_cleanup()
    implicit none

    close(unit=psi_xt_unit)
    close(unit=x_range_unit)
    close(unit=t_range_unit)
    close(unit=psi0_unit)
    close(unit=pot_unit)
    close(unit=wfunc_math_unit)

  end subroutine output_cleanup

  ! Time independent output
  subroutine output_time_indep()
    implicit none

    integer(dp) :: i_x

    do i_x = 1, n_x
       if (output_grids) then
          write(x_range_unit, dp_format) x_range(i_x)
          write(t_range_unit, dp_format) t_range(i_x)
       end if
       if (output_psi0) then
          write(psi0_unit, dp_format) abs(psi_arr(i_x))**2
       end if

       if (output_pot) then
          write(pot_unit, dp_format) pot_arr(i_x)
       end if
    end do

  end subroutine output_time_indep

  ! Time-dependent output
  subroutine output_time_dep(psi_arr)
    implicit none

    complex(dp), intent(in) :: psi_arr(:)
    integer(dp) :: i_x
    real(dp) :: norm, e_x, stdev_x

    if (output_psi_xt) then
       do i_x = 1, n_x
          write(psi_xt_unit, dp_format, advance="no") abs(psi_arr(i_x))**2
       end do
       write(psi_xt_unit, *)
    end if

    if (output_wfunc_math) then
       norm = numerics_norm(psi_arr)
       e_x = numerics_expec_x(psi_arr)
       stdev_x = numerics_stdev_x(psi_arr)

       write(wfunc_math_unit, "(3"//dp_format_raw//")") norm, e_x, stdev_x
    end if

  end subroutine output_time_dep

end module output
