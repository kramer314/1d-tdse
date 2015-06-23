! Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level directory of this distribution.

! File operations
module files
  implicit none

  private
  public :: files_ensure_dir

contains

  ! Ensure that a directory exists
  subroutine files_ensure_dir(dir)
    implicit none

    character(*), intent(in) :: dir
    call system("mkdir -p "//trim(dir))

  end subroutine files_ensure_dir

end module files
