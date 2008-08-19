! trans_coamps_to_dart
! --------------------
! This program pulls pieces out of the large COAMPS restart file,
! then assembles them into a state vector that can be used by DART.
! This includes two pieces of information - the current time and
! the actual state
program trans_coamps_to_dart

  use coamps_translate_mod, only : initialize_translator,         &
                                   generate_restart_filenames,    &
                                   open_coamps_files,             &
                                   coamps_restart_read_all_fields,&
                                   convert_coamps_state_to_dart,  &
                                   set_dart_current_time,         &
                                   open_dart_file, dart_write,    &
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
  logical, parameter :: READING_COAMPS_RESTART = .false.
  logical, parameter :: WRITING_DART_RESTART   = .true.

  call initialize_translator()

  call generate_restart_filenames(READING_COAMPS_RESTART)
  call open_coamps_files(READING_COAMPS_RESTART)
  
  call coamps_restart_read_all_fields()
  call convert_coamps_state_to_dart()

  call set_dart_current_time()

  call open_dart_file(WRITING_DART_RESTART)
  call dart_write()

  call finalize_translator()
end program trans_coamps_to_dart
