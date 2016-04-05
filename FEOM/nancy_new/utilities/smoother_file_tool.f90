! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program smoother_file_tool

! Program to copy the last restart files to Lag 1, Lag 1 to Lag 2, etc
! up to nlags.  has to replicate several namelist settings because 
! it's impossible for this tool to include namelists from other modules.

use types_mod,           only : r8
use time_manager_mod,    only : time_type, operator(<), operator(==),      &
                                set_time_missing, set_time,                &
                                operator(/=), print_time, print_date,      &
                                set_calendar_type, GREGORIAN, NO_CALENDAR, &
                                get_calendar_type

use utilities_mod,       only : register_module, do_output,                &
                                error_handler, nmlfileunit, E_MSG, E_ERR,  &
                                find_namelist_in_file, file_exist,         &
                                check_namelist_read, logfileunit,          &
                                do_nml_file, do_nml_term
                                
use assim_model_mod,     only : static_init_assim_model, get_model_size,   &
                                open_restart_read, open_restart_write,     &
                                awrite_state_restart, aread_state_restart, &
                                close_restart

use ensemble_manager_mod, only : init_ensemble_manager, ensemble_type,     &
                                 put_copy, prepare_to_write_to_vars, prepare_to_read_from_vars

use mpi_utilities_mod,    only : initialize_mpi_utilities, task_count,     &
                                 finalize_mpi_utilities


implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
character(len=128), parameter :: id  = "$Id$"

integer                 :: iunit, model_size, io, member, lag
type(ensemble_type)     :: ens_handle
character(len = 128)    :: ifile, ofile, msgstring
logical                 :: one_by_one, has_cal, seen_lag
character(len=16)       :: write_format

!----------------------------------------------------------------
! Most of these variables are namelist-controllable.
!
character(len = 128)  :: filter_restart_file_name   = "filter_restart", &
                         smoother_restart_file_name = "smoother_restart"
integer               :: ens_size                     = 1
integer               :: num_lags                     = 0
logical               :: single_restart_file_in       = .true.
logical               :: single_restart_file_out      = .true.
logical               :: write_binary_restart_files   = .true.
type(time_type)       :: data_time
logical               :: gregorian_cal = .true.
logical               :: print_only = .false.

namelist /smoother_file_tool_nml/  &
   filter_restart_file_name,     &
   smoother_restart_file_name,   &
   ens_size,                     &
   num_lags,                     &
   single_restart_file_in,       &
   single_restart_file_out,      &
   write_binary_restart_files,   &
   gregorian_cal,                &
   print_only


!----------------------------------------------------------------
!----------------------------------------------------------------

! This program should only be run with a single process
call initialize_mpi_utilities('smoother_file_tool')

if(task_count() > 1) &
   call error_handler(E_ERR,'smoother_file_tool','Only use single process', &
                      source,revision,revdate)

call register_module(id)

! Read the namelist entry and print it
call find_namelist_in_file("input.nml", "smoother_file_tool_nml", iunit)
read(iunit, nml = smoother_file_tool_nml, iostat = io)
call check_namelist_read(iunit, io, "smoother_file_tool_nml")

if (do_nml_file()) write(nmlfileunit, nml=smoother_file_tool_nml)
if (do_nml_term()) write(     *     , nml=smoother_file_tool_nml)

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

data_time        = set_time_missing()


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

! Handle the case of spinning up a smoother run - if the Lag files
! don't exist yet, allow that.  But once you've seen a Lag file, then
! it's an error if Lag N exists but any of Lag 1 to N-1 don't exist.
seen_lag = .false.

! if both in and out are single files, we can loop here over all ens members
! doing each one piecemeal; no need for space for all members in memory.
if (one_by_one) then

   ! Initialize the ens manager with enough room for a single ensemble member.
   call init_ensemble_manager(ens_handle, num_copies=1, num_vars=model_size)
   call prepare_to_write_to_vars(ens_handle)

   LAGLOOP: do lag=num_lags-1, 0, -1

      MEMLOOP: do member=1, ens_size
    
         ! if lag == 0, we're reading a restart file from filter.
         ! otherwise we're reading a previous lag
         ! add member number as a suffix: e.g. base.0000
         if (lag == 0) then
            write(ifile, "(A,A,I4.4)") trim(filter_restart_file_name), '.', member
         else
            write(ifile, '("Lag_",I5.5,"_",A,".",I4.4)') lag, trim(smoother_restart_file_name), member
         endif
   
         ! only check once during the member loop
         if (member == 1 .and. .not. seen_lag) then
            if (.not. file_exist(ifile)) then
                cycle LAGLOOP
            endif
         endif

         ! in all cases, the output is a Lag file.
         write(ofile, '("Lag_",I5.5,"_",A,".",I4.4)') lag+1, trim(smoother_restart_file_name), member
        
         if (.not. file_exist(ifile)) then
            write(msgstring, '(A)') 'Unable to open input file ', trim(ifile)
            call error_handler(E_ERR,'smoother_file_tool: ',msgstring)
         endif

         !------------------- Read restart from file ----------------------
         iunit = open_restart_read(ifile)
         call aread_state_restart(ens_handle%time(1), ens_handle%vars(:, 1), iunit)
         call close_restart(iunit)
         !------------------- Read restart from file ----------------------
         
         !-----------------  If only printing, print and exit --------------------
         if (print_only) then
            call print_info(model_size, has_cal, ens_handle%time(1))
            goto 10
         else
            if (member == 1) &
               call print_smoother_info(ifile, ofile, has_cal, ens_handle%time(1))
         endif
         !-----------------  If only printing, print and exit --------------------
   
   
         !------------------- Write restart to file -----------------------
         ! Output the restart file if requested
         iunit = open_restart_write(ofile, write_format)
         call awrite_state_restart(ens_handle%time(1), ens_handle%vars(:, 1), iunit)
         call close_restart(iunit)
         !------------------- Write restart to file -----------------------
   
      enddo MEMLOOP
   enddo LAGLOOP

else

   ! Initialize the ens manager with enough room for all ensemble members.
   ! Either read, write, or both will need this.
   call init_ensemble_manager(ens_handle, num_copies=ens_size, num_vars=model_size)
   call prepare_to_write_to_vars(ens_handle)

   ! make the defaults be a single filename, and overwrite them below if
   ! there are individual files.
   ifile = trim(filter_restart_file_name)
   ofile = trim(smoother_restart_file_name)
   
   !------------------- Read restart from file ----------------------
   ! If only one restart file on input, read it all up front.  For individuals
   ! read each one in the loop below.
   if (single_restart_file_in) then
      iunit = open_restart_read(ifile)
      ! Read in the advance time if not present
      do member=1, ens_size
         call aread_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit)
      enddo
      call close_restart(iunit)
    else
      ! loop over each ensemble member, reading in each individually.
      do member=1, ens_size
    
         ! add member number as a suffix: e.g. base.0000
         write(ifile, "(a,a,i4.4)") trim(filter_restart_file_name), '.', member
   
         iunit = open_restart_read(ifile)
         call aread_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit)
         call close_restart(iunit)

      enddo
   endif   
   !------------------- Read restart from file ----------------------
         
   !-----------------  If only printing, print and exit --------------------
   if (print_only) then
      call print_info(model_size, has_cal, ens_handle%time(1))
      goto 10
   endif
   !-----------------  If only printing, print and exit --------------------
  
   !------------------- Write restart to file -----------------------
   ! If only one restart file on output, write it all in one go.
   ! Otherwise, write each ens in a loop below.
   if (single_restart_file_out) then
      ! Output the restart file if requested
      iunit = open_restart_write(ofile, write_format)
      do member=1, ens_size
         call awrite_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit)
      enddo
      call close_restart(iunit)
    else
      ! loop over each ensemble member, reading in each individually.
      do member=1, ens_size
    
         ! add member number as a suffix: e.g. base.0000
         write(ofile, "(a,a,i4.4)") trim(smoother_restart_file_name), '.', member
   
         ! Output the restart file if requested
         iunit = open_restart_write(ofile, write_format)
         call awrite_state_restart(ens_handle%time(member), ens_handle%vars(:, member), iunit)
         call close_restart(iunit)

      enddo
   endif
   !------------------- Write restart to file -----------------------

endif

write(msgstring, '(A, I10)') 'Model size/restart data length =', model_size
call error_handler(E_MSG,'',msgstring)

! Early exit for the 'print_only' option
10 continue

call finalize_mpi_utilities()   ! now closes log file, too

contains

subroutine print_info(model_size, has_cal, data_time)
! called when the 'print_only' namelist item is true; print information
! about the timestamps in the file and the model size.
 integer, intent(in)                   :: model_size
 logical, intent(in)                   :: has_cal
 type(time_type), intent(in)           :: data_time

call error_handler(E_MSG,'','')

write(msgstring, '(A, I10)') 'Model size/restart data length =', model_size
call error_handler(E_MSG,'',msgstring)

call error_handler(E_MSG,'','')

call print_temporal(data_time, has_cal, "input file data time is ")

end subroutine print_info


subroutine print_smoother_info(ifile, ofile, has_cal, data_time)
 character(len=*), intent(in)          :: ifile, ofile
 logical, intent(in)                   :: has_cal
 type(time_type), intent(in)           :: data_time

call error_handler(E_MSG,'','')

call print_temporal(data_time, has_cal, "input file data time is ")

write(msgstring, '(A,A)')  'Input  filename: ', trim(ifile)
call error_handler(E_MSG,'',msgstring)

write(msgstring, '(A,A)')  'Output filename: ', trim(ofile)
call error_handler(E_MSG,'',msgstring)

call error_handler(E_MSG,'','')


end subroutine print_smoother_info


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


end program smoother_file_tool

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
