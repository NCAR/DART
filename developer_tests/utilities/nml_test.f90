! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program nml_test

use     types_mod, only : r8
use utilities_mod, only : register_module, error_handler, E_ERR, E_MSG,       &
                          open_file, close_file,                              &
                          find_namelist_in_file, check_namelist_read,         &
                          do_nml_file, do_nml_term, nmlfileunit,              &
                          initialize_utilities, finalize_utilities
implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

character(len=128) :: msgstring, msgstring2, tmpstring
integer :: iunit, io
integer  :: test1
logical  :: test2
real(r8) :: test3

! namelist items we are going to create/overwrite
namelist /nml_test_nml/  test1, test2, test3

! main code here
 
!@! option to create a dummy input.nml at runtime.
!@iunit = 11
!@open (iunit, file='input.nml', form='formatted',     &
!@      position='rewind', delim='apostrophe',    &
!@      action='write', status='replace', iostat=io)
!@write(iunit, "(A)") "&utilities_nml"
!@write(iunit, "(A)") "/"
!@close(iunit)

! initialize the dart libs
call initialize_utilities('nml_test')
call initialize_module()

! Create and write a namelist
test1 = 2
test2 = .true.
test3 = 3.14159

iunit = open_file("new_test.nml", "formatted", "write")
write(iunit, nml = nml_test_nml, iostat = io)
call close_file(iunit)

! Read back the namelist entry
call find_namelist_in_file("new_test.nml", "nml_test_nml", iunit)
read(iunit, nml = nml_test_nml, iostat = io)
call check_namelist_read(iunit, io, "nml_test_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=nml_test_nml)
if (do_nml_term()) write(     *     , nml=nml_test_nml)

call finalize_utilities('nml_test')

! end of main code


contains

!----------------------------------------------------------------------

subroutine initialize_module

  call register_module(source, revision, revdate)
  module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
