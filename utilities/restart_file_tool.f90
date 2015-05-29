! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program restart_file_tool

! Program to overwrite the time on each ensemble in a restart file.

use time_manager_mod,    only : time_type, operator(<), operator(==), &
                                set_time_missing, set_time,           &
                                operator(/=), print_time, print_date, &
                                set_calendar_type, GREGORIAN

use utilities_mod,       only : register_module, error_handler, nmlfileunit, &
                                E_MSG, E_ERR, find_namelist_in_file,         &
                                check_namelist_read, logfileunit,            &
                                do_nml_file, do_nml_term
                                
use assim_model_mod,     only : static_init_assim_model, get_model_size,   &
                                open_restart_read, open_restart_write,     &
                                awrite_state_restart, aread_state_restart, &
                                close_restart

use ensemble_manager_mod, only : init_ensemble_manager, ensemble_type, &
                                 prepare_to_write_to_vars

use mpi_utilities_mod,    only : initialize_mpi_utilities, task_count, &
                                 finalize_mpi_utilities

use types_mod,            only : i8


implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer                 :: iunit, io, member
integer(i8)             :: model_size
type(ensemble_type)     :: ens_handle
character(len = 128)    :: ifile, ofile, msgstring
logical                 :: one_by_one, has_cal
character(len=16)       :: write_format

!----------------------------------------------------------------
! Most of these variables are namelist-controllable.
!
character(len = 128)  :: input_file_name  = "filter_restart",        &
                         output_file_name = "filter_updated_restart"
integer               :: ens_size                     = 1
logical               :: single_restart_file_in       = .true.
logical               :: single_restart_file_out      = .true.
logical               :: write_binary_restart_files   = .true.
logical               :: overwrite_data_time          = .false.
type(time_type)       :: data_time, old_data_time
integer               :: new_data_days = -1, new_data_secs = -1
logical               :: input_is_model_advance_file  = .false.
logical               :: output_is_model_advance_file = .false.
logical               :: overwrite_advance_time       = .false.
type(time_type)       :: advance_time, old_advance_time
integer               :: new_advance_days = -1, new_advance_secs = -1
logical               :: gregorian_cal = .true.
logical               :: print_only = .false.

namelist /restart_file_tool_nml/  &
   input_file_name,              &
   output_file_name,             &
   ens_size,                     &
   single_restart_file_in,       &
   single_restart_file_out,      &
   write_binary_restart_files,   &
   overwrite_data_time,          &
   new_data_days,                &
   new_data_secs,                &
   input_is_model_advance_file,  &
   output_is_model_advance_file, &
   overwrite_advance_time,       &
   new_advance_days,             &
   new_advance_secs,             &
   gregorian_cal,                &
   print_only


!----------------------------------------------------------------
!----------------------------------------------------------------

! This program should only be run with a single process
call initialize_mpi_utilities('restart_file_tool')

if(task_count() > 1) &
   call error_handler(E_ERR,'restart_file_tool','Only use single process', &
                      source,revision,revdate)

call register_module(source,revision,revdate)

! Read the namelist entry and print it
call find_namelist_in_file("input.nml", "restart_file_tool_nml", iunit)
read(iunit, nml = restart_file_tool_nml, iostat = io)
call check_namelist_read(iunit, io, "restart_file_tool_nml")

if (do_nml_file()) write(nmlfileunit, nml=restart_file_tool_nml)
if (do_nml_term()) write(     *     , nml=restart_file_tool_nml)

! if you are not using a gregorian cal, set this to false
! in the namelist.
! NOTE: the namelist entry should probably be: calendar = 'string_name'
! and then the time manager should have a: call set_calendar_by_name('name')
! and then this is much more flexible.  it should probably have a call
! to get the calendar name as well so it could be used in print statements.
if (gregorian_cal) then
   call set_calendar_type(GREGORIAN)
   has_cal = .true.
else
   has_cal = .false.
endif



! ens_size is in the filter namelist, and the single restart file flags
! are in the ensemble manager namelist.  how do i get access to them
! here without replication?  write accessor routines for those parms
! in filter and the ens_mod code?  but filter is a main program, so ?

old_advance_time = set_time_missing()
advance_time     = set_time_missing()
old_data_time    = set_time_missing()
data_time        = set_time_missing()

! Time to reset model/data time to - if set in namelist
if (overwrite_data_time) then
   if (new_data_days >= 0 .and. new_data_secs >= 0) then
      data_time = set_time(new_data_secs, new_data_days)
   else
      call error_handler(E_ERR,'restart_file_tool','must specify data days and times', &
                         source,revision,revdate)
   endif
endif

! Time to reset advance/target time to - if set in namelist
if (overwrite_advance_time) then
   if (.not. output_is_model_advance_file) then
      call error_handler(E_ERR,'restart_file_tool','output_is_model_advance_file must be true to set advance time',&
                         source,revision,revdate)
   endif
   if (new_advance_days >= 0 .and. new_advance_secs >= 0) then
      advance_time = set_time(new_advance_secs, new_advance_days)
   else
      call error_handler(E_ERR,'restart_file_tool','must specify advance days and times',&
                         source,revision,revdate)
   endif
else
! Cannot output an advance time if input is not already one and if the user 
!  has not given us a time to use.
if (.not. input_is_model_advance_file .and. output_is_model_advance_file) then
      call error_handler(E_ERR,'restart_file_tool','overwrite_advance_time must be true if output file has advance time',&
                         source,revision,revdate)
   endif
endif

! Figure out the output format string.  (Input format automatically detected.)
if (write_binary_restart_files) then
   write_format = "unformatted"
else
   write_format = "formatted"
endif

! Initialize the model so we can get the size.
call static_init_assim_model()
model_size = get_model_size()

! If either the restart in or out data is all together in a single file, 
! then there has to be room for the entire dataset in memory at once --
! which is the size of the state vector times the ensemble count.
one_by_one = .not. (single_restart_file_in .or. single_restart_file_out) 

! if both in and out are single files, we can loop here over all ens members
! doing each one piecemeal; no need for space for all members in memory.
if (one_by_one) then

   ! Initialize the ens manager with enough room for a single ensemble member.
   call init_ensemble_manager(ens_handle, num_copies=1, num_vars=model_size)
   call prepare_to_write_to_vars(ens_handle)

   do member=1, ens_size
 
      ! add member number as a suffix: e.g. base.0000
      write(ifile, "(a,a,i4.4)") trim(input_file_name), '.', member
      write(ofile, "(a,a,i4.4)") trim(output_file_name), '.', member

      !------------------- Read restart from file ----------------------
      iunit = open_restart_read(ifile)
      ! Read in the advance time if present
      if (input_is_model_advance_file) then
         call aread_state_restart(ens_handle%time(1), ens_handle%vars(:, 1), iunit, old_advance_time)
      else
         call aread_state_restart(ens_handle%time(1), ens_handle%vars(:, 1), iunit)
      endif
      call close_restart(iunit)
      !------------------- Read restart from file ----------------------
      
      !-----------------  If only printing, print and exit --------------------
      if (print_only) then
         call print_info(model_size, has_cal, ens_handle%time(1), &
                         input_is_model_advance_file, old_advance_time)
         goto 10
      endif
      !-----------------  If only printing, print and exit --------------------

      !-----------------  Update data, advance times if requested -------------
      old_data_time = ens_handle%time(1)
      if (overwrite_data_time) &
         ens_handle%time(1) = data_time
      if (.not. overwrite_advance_time) &
          advance_time = old_advance_time 
      !-----------------  Update data, advance times if requested -------------

      !------------------- Write restart to file -----------------------
      ! Output the restart file if requested; Force to binary for bitwise reproducing
      ! use in filter and perfect_model obs with shell advance options
      iunit = open_restart_write(ofile, write_format)
      if (output_is_model_advance_file) then
         call awrite_state_restart(ens_handle%time(1), ens_handle%vars(:, 1), iunit, advance_time)
      else
         call awrite_state_restart(ens_handle%time(1), ens_handle%vars(:, 1), iunit)
      endif
      call close_restart(iunit)
      !------------------- Write restart to file -----------------------

   enddo

else

   ! Initialize the ens manager with enough room for all ensemble members.
   ! Either read, write, or both will need this.
   call init_ensemble_manager(ens_handle, num_copies=ens_size, num_vars=model_size)
   call prepare_to_write_to_vars(ens_handle)

   ! make the defaults be a single filename, and overwrite them below if
   ! there are individual files.
   ifile = trim(input_file_name)
   ofile = trim(output_file_name)
   
   !------------------- Read restart from file ----------------------
   ! If only one restart file on input, read it all up front.  For individuals
   ! read each one in the loop below.
   if (single_restart_file_in) then
      iunit = open_restart_read(ifile)
      ! Read in the advance time if not present
      do member=1, ens_size
         if (input_is_model_advance_file) then
            call aread_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit, old_advance_time)
         else
            call aread_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit)
         endif
      enddo
      call close_restart(iunit)
    else
      ! loop over each ensemble member, reading in each individually.
      do member=1, ens_size
    
         ! add member number as a suffix: e.g. base.0000
         write(ifile, "(a,a,i4.4)") trim(input_file_name), '.', member
   
         iunit = open_restart_read(ifile)
         ! Read in the advance time if present
         if (input_is_model_advance_file) then
            call aread_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit, old_advance_time)
         else
            call aread_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit)
         endif
         call close_restart(iunit)

      enddo
   endif   
   !------------------- Read restart from file ----------------------
         
   !-----------------  If only printing, print and exit --------------------
   if (print_only) then
      call print_info(model_size, has_cal, ens_handle%time(1), &
                      input_is_model_advance_file, old_advance_time)
      goto 10
   endif
   !-----------------  If only printing, print and exit --------------------
  
   !-----------------  Update data, advance times if requested -----------
   do member=1, ens_size
      old_data_time = ens_handle%time(member)
      if (overwrite_data_time) &
         ens_handle%time(member) = data_time
      if (.not. overwrite_advance_time) &
          advance_time = old_advance_time 
   enddo
   !-----------------  Update data, advance times if requested -----------

   !------------------- Write restart to file -----------------------
   ! If only one restart file on output, write it all in one go.
   ! Otherwise, write each ens in a loop below.
   if (single_restart_file_out) then
      ! Output the restart file if requested; Force to binary for bitwise reproducing
      ! use in filter and perfect_model obs with shell advance options
      iunit = open_restart_write(ofile, write_format)
      do member=1, ens_size
         if (output_is_model_advance_file) then
            call awrite_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit, advance_time)
         else
            call awrite_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit)
         endif
      enddo
      call close_restart(iunit)
    else
      ! loop over each ensemble member, reading in each individually.
      do member=1, ens_size
    
         ! add member number as a suffix: e.g. base.0000
         write(ofile, "(a,a,i4.4)") trim(output_file_name), '.', member
   
         ! Output the restart file if requested; Force to binary for bitwise reproducing
         ! use in filter and perfect_model obs with shell advance options
         iunit = open_restart_write(ofile, write_format)
         if (output_is_model_advance_file) then
            call awrite_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit, advance_time)
         else
            call awrite_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit)
         endif
         call close_restart(iunit)

      enddo
   endif
   !------------------- Write restart to file -----------------------

endif

write(msgstring, '(A, i8)') 'Model size/restart data length =', model_size
call error_handler(E_MSG,'',msgstring)

if (old_advance_time .ne. set_time_missing()) then
   call print_temporal(old_advance_time, has_cal, "input file had an advance_time, which was ")
endif
if ((advance_time .ne. set_time_missing()) .and. output_is_model_advance_file) then
   call print_temporal(advance_time, has_cal, "output file advance_time is now set to ")
endif
if (old_data_time .ne. set_time_missing()) then
   call print_temporal(old_data_time, has_cal, "input file data_time was ")
endif
if ((data_time .ne. set_time_missing()) .or. overwrite_data_time) then
   call print_temporal(data_time, has_cal, "output file data_time is now set to ")
endif

! Early exit for the 'print_only' option
10 continue

call finalize_mpi_utilities()   ! now closes log file, too

contains

subroutine print_info(model_size, has_cal, data_time, is_advance, advance_time)
! called when the 'print_only' namelist item is true; print information
! about the timestamps in the file and the model size.
 integer(i8), intent(in)               :: model_size
 logical, intent(in)                   :: has_cal
 type(time_type), intent(in)           :: data_time
 logical, intent(in)                   :: is_advance
 type(time_type), intent(in), optional :: advance_time


call error_handler(E_MSG,'','')

write(msgstring, '(A, i8)') 'Model size/restart data length =', model_size
call error_handler(E_MSG,'',msgstring)

call error_handler(E_MSG,'','')

call print_temporal(data_time, has_cal, "input file data time is ")
if (is_advance) then
   call print_temporal(advance_time, has_cal, "input file has an advance_time, which is ")
endif

end subroutine print_info

subroutine print_temporal(timeval, has_cal, label)
! the print_date/time routines print only one place.  rather than replicating
! all the code, make a single call which prints both time, and date format if
! there is a calendar, and print both to stdout and to the log file.
 type(time_type), intent(in)  :: timeval
 logical, intent(in)          :: has_cal
 character(len=*), intent(in) :: label

call print_time(timeval, label)
call print_time(timeval, label, logfileunit)

if (has_cal) then
   call print_date(timeval, label)
   call print_date(timeval, label, logfileunit)
endif

end subroutine print_temporal

end program restart_file_tool

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
