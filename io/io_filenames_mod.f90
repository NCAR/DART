module io_filenames_mod

!> \defgroup io_filenames io_filenames
!> Aim is to store the io filenames
!>  * Restarts
!>  * Diagnostics
!>  * Inflation files
!>
!> Any module can set the filenames here, then state_vector_io_mod
!> can read from this module to get the filenames.
!> Maybe this is a bit lazy, but I just want a way for the 
!> different modules to set the filenames
!>
!> Diagnostic files could have different netcdf variable ids
!> @{

use utilities_mod, only : do_nml_file, nmlfileunit, do_nml_term, check_namelist_read, &
                          find_namelist_in_file

implicit none

private

! These should probably be set and get functions rather than 
! direct access

public :: io_filenames_init, restart_files_in, restart_files_out, extras_in, extras_out

! How do people name there restart files?
! What about domains?
integer, parameter :: max_num_files = 500
integer, parameter :: max_num_domains = 5


! Why not have an in and an out list?
character(len=2048) :: restart_files_in(max_num_files)  = '' ! list of input restart files
character(len=2048) :: restart_files_out(max_num_files) = '' ! list of output restart files

character(len=2048) :: extras_in(6*max_num_domains) ! list of extras in
character(len=2048) :: extras_out(6*max_num_domains) ! list of extras out

namelist / io_filenames_nml / restart_files_in, restart_files_out, extras_in, extras_out

contains

!----------------------------------
subroutine io_filenames_init()

integer :: iunit, io

!call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "io_filenames_nml", iunit)
read(iunit, nml = io_filenames_nml, iostat = io)
call check_namelist_read(iunit, io, "io_filenames_nml")

! Write the namelist values to the log file
if (do_nml_file()) write(nmlfileunit, nml=io_filenames_nml)
if (do_nml_term()) write(     *     , nml=io_filenames_nml)

end subroutine io_filenames_init
!----------------------------------

!> @}
end module io_filenames_mod