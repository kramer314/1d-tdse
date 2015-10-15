! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution.

module output
  use progvars
  use files, only: files_ensure_dir
  use wfmath, only: wfmath_norm, wfmath_expec_x, wfmath_stdev_x, &
       wfmath_autocorr, wfmath_prob_flux

  implicit none

  private

  ! Generic module functions
  public :: output_init
  public :: output_cleanup

  public :: output_time_indep
  public :: output_time_dep

  ! Private module variables

  ! Output file unit numbers
  integer, parameter :: psi_xt_unit = 99
  integer, parameter :: flux_xt_unit = 98
  integer, parameter :: autocorr_unit = 97
  integer, parameter :: wfunc_checks_unit = 96

  integer, parameter :: x_range_unit = 89
  integer, parameter :: t_range_unit = 88
  integer, parameter :: psi0_unit = 87

contains

  ! Module initialization
  subroutine output_init()

    call files_ensure_dir(output_dir)
    if (output_psi_xt) then
       open(unit=psi_xt_unit, file=trim(output_dir)//trim(psi_xt_fname))
    end if

    if (output_flux_xt) then
       open(unit=flux_xt_unit, file=trim(output_dir)//trim(flux_xt_fname))
    end if

    if (output_autocorr) then
       open(unit=autocorr_unit, file=trim(output_dir)//trim(autocorr_fname))
    end if

    if (output_wfunc_checks) then
       open(unit=wfunc_checks_unit, file=trim(output_dir)// &
            trim(wfunc_checks_fname))
    end if

    if (output_grids) then
       open(unit=x_range_unit, file=trim(output_dir)//trim(x_range_fname))
       open(unit=t_range_unit, file=trim(output_dir)//trim(t_range_fname))
    end if

    if (output_psi0) then
       open(unit=psi0_unit, file=trim(output_dir)//trim(psi0_fname))
    end if

  end subroutine output_init


  ! Module cleanup
  subroutine output_cleanup()

    close(unit=psi_xt_unit)
    close(unit=flux_xt_unit)
    close(unit=autocorr_unit)
    close(unit=wfunc_checks_unit)
    close(unit=x_range_unit)
    close(unit=t_range_unit)
    close(unit=psi0_unit)

  end subroutine output_cleanup


  ! Output time-independent quantities
  subroutine output_time_indep()

    integer :: i_x

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
    integer :: i_x
    real(dp) :: norm, expec_x, stdev_x
    complex(dp) :: autocorr

    if (output_psi_xt) then

       do i_x = 1, n_x
          if (mod(i_x, print_mod_x) .eq. 0) then
             write(psi_xt_unit, dp_format, advance="no") abs(psi_arr(i_x))**2
          end if
       end do
       write(psi_xt_unit, *)

    end if

    if (output_flux_xt) then

       call wfmath_prob_flux(psi_arr, flux_arr)

       do i_x = 1, n_x
          if (mod(i_x, print_mod_x) .eq. 0) then
             write(flux_xt_unit, dp_format, advance="no") abs(flux_arr(i_x))**2
          end if
       end do
       write(flux_xt_unit, *)

    end if

    if (output_autocorr) then

       autocorr = wfmath_autocorr(psi_arr, psi0_arr)
       write(autocorr_unit, dp_format) abs(autocorr)**2

    end if

    if (output_wfunc_checks) then

       norm = wfmath_norm(psi_arr)
       expec_x = wfmath_expec_x(psi_arr)
       stdev_x = wfmath_stdev_x(psi_arr)

       write(wfunc_checks_unit, "(3"//dp_format_raw//")") norm, expec_x, stdev_x

    end if

  end subroutine output_time_dep

end module output
