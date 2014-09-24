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
use model_mod,     only : construct_file_name

implicit none

private

! These should probably be set and get functions rather than 
! direct access

public :: io_filenames_init, restart_files_in, restart_files_out

! How do people name there restart files?
! What about domains?
integer, parameter :: max_num_files = 5000

! public arrays of filenames. Do we need arrays for restarts AND extras?
character(len=2048), allocatable :: restart_files_in(:,:), restart_files_out(:,:)

! Namelist options
character(len=512) :: restart_in_stub  = 'wrfinput.nc'
character(len=512) :: restart_out_stub = '/Output/wrfinput.nc'

! Should probably get num_domains, num_restarts from elsewhere. In here for now
namelist / io_filenames_nml / restart_in_stub, restart_out_stub

contains

!----------------------------------
!> read namelist and set up filename arrays
subroutine io_filenames_init(ens_size, num_domains)

integer, intent(in) :: ens_size
integer, intent(in) :: num_domains
integer :: iunit, io
integer :: dom, num_files, i

!call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "io_filenames_nml", iunit)
read(iunit, nml = io_filenames_nml, iostat = io)
call check_namelist_read(iunit, io, "io_filenames_nml")

! Write the namelist values to the log file
if (do_nml_file()) write(nmlfileunit, nml=io_filenames_nml)
if (do_nml_term()) write(     *     , nml=io_filenames_nml)

num_files = ens_size + 6

allocate(restart_files_in(num_files, num_domains))
allocate(restart_files_out(num_files, num_domains))

do dom = 1, num_domains
   do i = 1, num_files  ! restarts
      restart_files_in(i, dom)  = construct_file_name(restart_in_stub, dom, i)
      restart_files_out(i, dom) = construct_file_name(restart_out_stub, dom, i)
   enddo
enddo

end subroutine io_filenames_init
!----------------------------------

!----------------------------------
subroutine end_io_filenames()

deallocate(restart_files_in, restart_files_out)

end subroutine end_io_filenames

!----------------------------------

!> @}
end module io_filenames_mod