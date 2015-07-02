! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

module output
  use progvars
  use files, only: files_ensure_dir
  use wfmath, only: wfmath_norm, wfmath_expec_x, wfmath_stdev_x, &
       wfmath_autocorr

  implicit none

  private

  ! Generic module functions
  public :: output_init
  public :: output_cleanup

  public :: output_time_indep
  public :: output_time_dep

  ! Private module variables

  ! Output file unit numbers
  integer(dp), parameter :: psi_xt_unit = 99
  integer(dp), parameter :: x_range_unit = 98
  integer(dp), parameter :: t_range_unit = 97
  integer(dp), parameter :: psi0_unit = 96
  integer(dp), parameter :: wfunc_math_unit = 95

contains

  ! Module initialization
  subroutine output_init()

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

    if (output_wfunc_math) then
       open(unit=wfunc_math_unit, file=trim(output_dir)// &
            trim(wfunc_math_fname))
    end if

  end subroutine output_init

  ! Module cleanup
  subroutine output_cleanup()

    close(unit=psi_xt_unit)
    close(unit=x_range_unit)
    close(unit=t_range_unit)
    close(unit=psi0_unit)
    close(unit=wfunc_math_unit)

  end subroutine output_cleanup

  ! Output time-independent quantities
  subroutine output_time_indep()

    integer(dp) :: i_x

    do i_x = 1, n_x

       if (mod(i_x, print_mod_x) .eq. 0) then

          if (output_grids) then
             write(x_range_unit, dp_format) x_range(i_x)
             write(t_range_unit, dp_format) t_range(i_x)
          end if

          if (output_psi0) then
             write(psi0_unit, dp_format) abs(psi_arr(i_x))**2
          end if

       end if

    end do

  end subroutine output_time_indep

  ! Output time-dependent quantities at a given propagation step
  !
  ! psi_arr :: wavefunction array at a particular time step
  subroutine output_time_dep(psi_arr)

    complex(dp), intent(in) :: psi_arr(:)
    integer(dp) :: i_x
    real(dp) :: norm, e_x, stdev_x
    complex(dp) :: autocorr

    if (output_psi_xt) then

       do i_x = 1, n_x

          if (mod(i_x, print_mod_x) .eq. 0) then
             write(psi_xt_unit, dp_format, advance="no") abs(psi_arr(i_x))**2
          end if

       end do
       write(psi_xt_unit, *)

    end if

    if (output_wfunc_math) then

       ! Calculate and output norm, <x>, stdev(x), and autocorrelation
       norm = wfmath_norm(psi_arr)
       e_x = wfmath_expec_x(psi_arr)
       stdev_x = wfmath_stdev_x(psi_arr)
       autocorr = wfmath_autocorr(psi_arr, psi0_arr)

       write(wfunc_math_unit, "(4"//dp_format_raw//")") norm, e_x, stdev_x, &
            abs(autocorr)**2

    end if

  end subroutine output_time_dep

end module output
