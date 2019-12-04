! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program error_handler_test

use     types_mod,     only : r8
use utilities_mod,     only : register_module, error_handler, E_ERR, E_MSG, E_ALLMSG, &
                              find_namelist_in_file, check_namelist_read,             &
                              do_nml_file, do_nml_term, nmlfileunit
use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities,       &
                              task_sync, my_task_id

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

character(len=128) :: msgstring
integer :: iunit, io
character(len=8) :: task_id
integer  :: test1
logical  :: test2
real(r8) :: test3

! namelist items we are going to create/overwrite
namelist /error_handler_test_nml/  test1, test2, test3

! main code here
 
!----------------------------------------------------------------------

! initialize the dart libs
call initialize_mpi_utilities('error_handler_test')
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "error_handler_test_nml", iunit)
read(iunit, nml = error_handler_test_nml, iostat = io)
call check_namelist_read(iunit, io, "error_handler_test_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=error_handler_test_nml)
if (do_nml_term()) write(     *     , nml=error_handler_test_nml)

write(task_id, '(i5)' ) my_task_id()

write(msgstring, *) 'test E_MSG from task_id = ', adjustl(trim(task_id))
call error_handler(E_MSG,'E_MSG message', msgstring, source, revision, revdate)

write(msgstring, *) 'test E_ALLMSG from task_id = ', task_id
call error_handler(E_ALLMSG,'E_ALLMSG message', msgstring, source, revision, revdate)

call task_sync()

! write(msgstring, *) 'test E_ERR'
! call error_handler(E_ERR,'E_ERR message', msgstring, source, revision, revdate)

! finalize error_handler_test
call error_handler(E_MSG,'error_handler_test','Finished successfully.',source,revision,revdate)
call finalize_mpi_utilities()

! end of main code

!----------------------------------------------------------------------

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
