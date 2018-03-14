! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program trans_dart_to_coamps

! trans_dart_to_coamps
! --------------------
! This program disassembles the DART state vector then inserts
! the fields into a pre-existing COAMPS restart file.  Since 
! contained in the DART file is a target time for the integration,
! write out a file that can be read by the COAMPS run scripts so
! they can modify their namelist accordingly.

use coamps_translate_mod, only : initialize_translator,     &
                                 open_dart_file,            &
                                 set_dart_current_time,     &
                                 get_dtg,                   &
                                 generate_coamps_filenames, &
                                 open_coamps_files,         &
                                 coamps_write_all_fields,   &
                                 write_pickup_file,         &
                                 finalize_translator

use coamps_hdf5_mod, only : initialize_hdf5_wfile, &
                            open_hdf5_file, &
                            close_hdf5_file, &
                            hdf5_file_write

use utilities_mod, only : E_ERR, error_handler, &
                          initialize_utilities, &
                          finalize_utilities

implicit none

character(len=512) :: string1

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

character(len=*), parameter :: routine = 'trans_dart_to_coamps'

! The translation module uses internal flags for whether it's
! reading or writing - these are just aliases so it's clearer
! what the calls are saying
logical, parameter :: READING_DART = .false.     ! false means read dart
logical, parameter :: WRITING_COAMPS  = .true.   ! true means write coamps

integer :: ierr
character(len=10) :: cdtgp1

call initialize_utilities(routine)

call initialize_translator()
!call initalize_hdf5()

call open_dart_file(READING_DART)
call set_dart_current_time()

call generate_coamps_filenames(WRITING_COAMPS)
call open_coamps_files(WRITING_COAMPS)

cdtgp1 = get_dtg('next')
call initialize_hdf5_wfile(hdf5_file_write, cdtg=cdtgp1, h5_prefix='CoampsUpdate')
call open_hdf5_file(hdf5_file_write, ierr)

if (ierr /= 0) then
   write(string1,*)'open_hdf_file error code ',ierr
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
endif

!>@todo still have to apply clamping

call coamps_write_all_fields()

call close_hdf5_file(hdf5_file_write, ierr)
if (ierr /= 0) then
   write(string1,*)'close_hdf_file error code ',ierr
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
endif

! Script uses a file to get the start/target time
call write_pickup_file()

call finalize_translator()
call finalize_utilities(routine)

end program trans_dart_to_coamps

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
