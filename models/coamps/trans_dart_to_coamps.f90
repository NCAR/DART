! trans_dart_to_coamps
! --------------------
! This program disassembles the DART state vector then inserts
! the fields into a pre-existing COAMPS restart file.  Since 
! contained in the DART file is a target time for the integration,
! write out a file that can be read by the COAMPS run scripts so
! they can modify their namelist accordingly.
program trans_dart_to_coamps

  use coamps_translate_mod, only : initialize_translator,          &
                                   open_dart_file, dart_read,      &
                                   fix_negative_values,            &
                                   convert_dart_state_to_coamps,   &
                                   get_dart_current_time,          &
                                   get_dart_target_time,           &
                                   generate_restart_filenames,     &
                                   open_coamps_files,              &
                                   coamps_restart_write_all_fields,&
                                   write_pickup_file,              &
                                   finalize_translator

  implicit none

  ! Modified automatically by Subversion 
  character(len=128) :: &
       source = "$URL$",&
       revision = "$Revision$", &
       revdate = "$Date$"

  ! The translation module uses internal flags for whether it's
  ! reading or writing - these are just aliases so it's clearer
  ! what the calls are saying
  logical, parameter :: READING_DART_RESTART   = .false.
  logical, parameter :: WRITING_COAMPS_RESTART = .true.

  call initialize_translator()

  call open_dart_file(READING_DART_RESTART)

  call dart_read()

  call fix_negative_values()

  call convert_dart_state_to_coamps()

  call get_dart_current_time()
  call get_dart_target_time()

  call generate_restart_filenames(WRITING_COAMPS_RESTART)

  call open_coamps_files(WRITING_COAMPS_RESTART)

  call coamps_restart_write_all_fields()

  ! Script uses a file to get the start/target time
  call write_pickup_file()

  call finalize_translator()
end program trans_dart_to_coamps
