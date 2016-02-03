! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_state_structure

use        types_mod, only : r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG,       &
                             open_file, close_file, nc_check, get_next_filename, &
                             find_namelist_in_file, check_namelist_read,         &
                             do_nml_file, do_nml_term, nmlfileunit,              &
                             initialize_utilities, finalize_utilities
use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities,       &
                              task_sync, my_task_id

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

integer :: iunit, io
integer :: test1

! namelist items we are going to create/overwrite
namelist /test_state_structure_nml/  test1

! main code here
 
! initialize the dart libs
call initialize_module()
! Read back the namelist entry
call find_namelist_in_file("input.nml", "test_state_structure_nml", iunit)
read(iunit, nml = test_state_structure_nml, iostat = io)
call check_namelist_read(iunit, io, "test_state_structure_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=test_state_structure_nml)
if (do_nml_term()) write(     *     , nml=test_state_structure_nml)

! finalize test_state_structure
call error_handler(E_MSG,'test_state_structure','Finished successfully.',source,revision,revdate)
call finalize_mpi_utilities()

! end of main code


contains

!----------------------------------------------------------------------

subroutine initialize_module

  call initialize_mpi_utilities('test_state_structure')
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
