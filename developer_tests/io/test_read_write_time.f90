! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!>@todo  FIXME - add more tests ... wrong calendars, etc.

program test_read_write_time

use            types_mod, only : r8, i8
use        utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                                 find_namelist_in_file, check_namelist_read,   &
                                 do_nml_file, do_nml_term, nmlfileunit, to_upper, &
                                 initialize_utilities, finalize_utilities

use netcdf_utilities_mod, only : nc_open_file_readwrite, nc_close_file
use     dart_time_io_mod, only : read_model_time, write_model_time
use     time_manager_mod, only : time_type, set_calendar_type, get_calendar_type, &
                                 set_time, print_time, operator(+)

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'test_read_write_time.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

character(len=512) :: msgstring

integer :: iunit, io, ncid, i
integer, parameter :: MAXFILES = 100

type(time_type) :: mytime

! namelist variables
character(len=256) :: input_file(MAXFILES)  = ""
logical            :: verbose = .false.

! namelist items we are going to create/overwrite
namelist /test_read_write_time_nml/ input_file, verbose 


! main code here
 
! initialize the dart libs
call initialize_module()

! Read back the namelist entry
call find_namelist_in_file("input.nml", "test_read_write_time_nml", iunit)
read(iunit, nml = test_read_write_time_nml, iostat = io)
call check_namelist_read(iunit, io, "test_read_write_time_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=test_read_write_time_nml)
if (do_nml_term()) write(     *     , nml=test_read_write_time_nml)

call error_handler(E_MSG, "", "")

! intent is to open a list of netcdf files with various permutations
! of time variable, dimension, size, etc and see if read/write model time
! routines (the default ones) work or error out correctly

do i = 1, MAXFILES
   if (input_file(i) == "") exit

   ! to test:
   ! function read_model_time(filename)
   ! subroutine write_model_time(ncid, dart_time)
   
   mytime = read_model_time(input_file(i))
   
   call print_time(mytime,'read_model_time first')
   
   mytime = mytime + set_time(0, 1)
   
   ncid = nc_open_file_readwrite(input_file(i))
   call write_model_time(ncid, mytime)
   call nc_close_file(ncid)
   
   mytime = read_model_time(input_file(i))
   call print_time(mytime,'read_model_time second')

enddo


call finalize_module()

! end of main code


contains

!----------------------------------------------------------------------

subroutine initialize_module

call initialize_utilities('test_read_write_time')
call register_module(source, revision, revdate)

end subroutine initialize_module

!----------------------------------------------------------------------

subroutine finalize_module

call finalize_utilities('test_read_write_time')

end subroutine finalize_module

!----------------------------------------------------------------------

end program

