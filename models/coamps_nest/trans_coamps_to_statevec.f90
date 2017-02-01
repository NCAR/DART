! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program trans_coamps_to_statevec

! trans_coamps_to_statevec
! ------------------------
! This program pulls pieces out of the large COAMPS restart file,
! then assembles them into a state vector that can be used to show
! what DART is seeing for easier reading into programs like MATLAB
! or debugging.  This program is almost identical to the
! trans_coamps_to_dart program but does not output time.

  use coamps_translate_mod, only : initialize_translator,         &
                                   generate_coamps_filenames,     &
                                   open_coamps_files,             &
                                   coamps_read_all_fields,        &
                                   convert_coamps_state_to_dart,  &
                                   set_dart_current_time,         &
                                   open_dart_file, dart_write,    &
                                   finalize_translator

  implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

  ! The translation module uses internal flags for whether it's
  ! reading or writing - help me remember what they are
  logical, parameter :: READING_COAMPS  = .false.
  logical, parameter :: WRITING_DART    = .true.
  logical, parameter :: DONT_WRITE_TIME = .true.     
 
  ! Make sure translation module has everything it needs
  call initialize_translator()
  call generate_coamps_filenames(READING_COAMPS)

  ! Open up the COAMPS file(s)
  call open_coamps_files(READING_COAMPS)
  
  ! Read in and COAMPS convert state vector
  call coamps_read_all_fields()
  call convert_coamps_state_to_dart()

  ! Put this in DART, but don't write the file time
  call open_dart_file(WRITING_DART)
  call dart_write(DONT_WRITE_TIME)

  ! Clean up
  call finalize_translator()

end program trans_coamps_to_statevec

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
