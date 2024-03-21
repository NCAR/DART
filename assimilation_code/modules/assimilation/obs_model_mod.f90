! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module obs_model_mod

use utilities_mod,        only : error_handler, E_ERR, E_MSG, E_WARN, &
                                 get_unit, file_exist, set_output
use assim_model_mod,      only : get_closest_state_time_to,         &
                                 get_model_time_step,  adv_1step

use state_vector_io_mod,  only : read_state, write_state

use obs_sequence_mod,     only : obs_sequence_type, obs_type,  &
                                 get_obs_def, init_obs, destroy_obs, get_num_copies, &
                                 get_num_qc, get_first_obs, get_next_obs_from_key,   &
                                 get_obs_time_range
use obs_def_mod,          only : obs_def_type, get_obs_def_time
use time_manager_mod,     only : time_type, set_time, get_time,                       &
                                 operator(/=), operator(>), operator(-),              &
                                 operator(/), operator(+), operator(<), operator(==), &
                                 operator(<=), operator(>=)
use ensemble_manager_mod, only : get_ensemble_time, ensemble_type, map_task_to_pe

use mpi_utilities_mod,    only : my_task_id, task_sync, block_task, &
                                 sum_across_tasks, shell_execute, my_task_id
use io_filenames_mod,     only : file_info_type

implicit none
private

public :: move_ahead, advance_state, set_obs_model_trace, have_members

character(len=*), parameter :: source = 'obs_model_mod.f90'

logical :: module_initialized  = .false.
integer :: print_timestamps    = 0
integer :: print_trace_details = 0

logical :: debug = .false.   ! set to true to get more status msgs

! Module storage for writing error messages
character(len=512) :: errstring, errstring1, errstring2

contains

!======================================================================

subroutine initialize_module

module_initialized = .true.

end subroutine initialize_module

!-------------------------------------------------------------------------

subroutine move_ahead(ens_handle, ens_size, seq, last_key_used, window_time, &
   key_bounds, num_obs_in_set, curr_ens_time, next_ens_time)

! Based on the current ens time and the time of the next available
! observation, compute whether and how far the ensemble copies of the
! state need to be advanced.  This returns the number of obs in the
! next assimilation window, and both the current and expected next data
! times.  It NO LONGER advances the model - the calling code must call
! advance_state() itself in order to do that if next time is /= curr time.

implicit none

type(ensemble_type),     intent(inout) :: ens_handle
integer,                 intent(in)    :: ens_size
type(obs_sequence_type), intent(in)    :: seq
integer,                 intent(in)    :: last_key_used
type(time_type),         intent(in)    :: window_time
integer,                 intent(out)   :: key_bounds(2), num_obs_in_set
type(time_type),         intent(inout) :: curr_ens_time
type(time_type),         intent(out)   :: next_ens_time

type(time_type)    :: next_time, time2, start_time, end_time, delta_time, ens_time
type(obs_type)     :: observation
type(obs_def_type) :: obs_def
logical            :: is_this_last, is_there_one, out_of_range, leaving_early

! Initialize if needed
if(.not. module_initialized) then
   call initialize_module
   module_initialized = .true.
endif

! Prepare values for 'early' returns
num_obs_in_set  =   0
key_bounds(1:2) = -99
leaving_early   =   .false.
!curr_ens_time   = set_time(0,0)
next_ens_time   = set_time(0,0)


! If none of my copies are regular ensemble members no need to advance. 
! This is true either if I have no copies or if the copy #s I have are
! all larger than the ensemble count.  This does assume that ensemble 
! copies are in the first ens_size global copies of ensemble and that 
! global indices are monotonically increasing in each pes local copies
! (note: Violating the private boundaries on ens_handle by direct access)
! This was the original code:
!   if((ens_handle%my_num_copies < 1) .or. (ens_handle%my_copies(1) > ens_size)) then
! but if the compiler tried to evaluate both halves at the same time,
! (which it's allowed to do), it got an access violation if num copies was 0.

! if this task has no ensembles to advance, return directly from here.
! the calling code will sort out syncing up in the parallel case.
if (.not. have_members(ens_handle, ens_size)) return


! if you get here, this task has at least one ensemble to try to advance;
! it is possible we are at the end of the observations and there in fact
! is no need to advance.  if so, can return.

! ens_handle%my_pe 0 does the output.
! Don't want two pes outputing if task 0 also has a copy
! FIMXE Commment 
if ( map_task_to_pe(ens_handle, 0) >= ens_handle%num_copies .and. &
   ens_handle%my_pe == 0 .and. my_task_id() /= 0) then
  call set_output(.true.)
endif


! Initialize a temporary observation type to use
! after here, must delete observation before returning.
call init_obs(observation, get_num_copies(seq), get_num_qc(seq))

! Get the next observation in the sequence that hasn't already been assimilated
if(last_key_used > 0) then
   call get_next_obs_from_key(seq, last_key_used, observation, is_this_last)
   ! If the end of the observation sequence has been passed, return
   if(is_this_last) leaving_early = .true.
else
   is_there_one = get_first_obs(seq, observation)
   if(.not. is_there_one) leaving_early = .true.
endif

if (leaving_early) then
   ! need to destroy obs here before returning
   call destroy_obs(observation)

   if (ens_handle%my_pe == 0 .and. my_task_id() /= 0) then
    call set_output(.false.)
   endif

   return
endif

! Get the time of this observation
call get_obs_def(observation, obs_def)
next_time = get_obs_def_time(obs_def)

! Get the time of the ensemble, assume consistent across all ensemble copies
call get_ensemble_time(ens_handle, 1, ens_time)

! Compute the model time step and center a window around the closest time
delta_time = get_model_time_step()

! print out current window, if requested
if (print_trace_details > 0) then
   if(delta_time / 2 > ens_time) then
      start_time = set_time(0, 0)
   else
      start_time = ens_time - delta_time / 2 + set_time(1, 0)
   endif
   end_time = ens_time + delta_time / 2
   call timechat(ens_time,    'move_ahead', .false.,        'Current model data time            is: ')
   call timechat(start_time,  'move_ahead', .false.,        'Current assimilation window starts at: ')
   call timechat(next_time,   'move_ahead', .false.,        'Next available observation time    is: ')
   call timechat(end_time,    'move_ahead', .false.,        'Current assimilation window ends   at: ')
  !call timechat(delta_time,  'move_ahead', .false.,        'Width of assimilation window       is: ')
endif

! now recompute for the next window, so the code below can remain unchanged.
! Figure out what time to advance the model to.

! More control over time window use of observations would come in here
time2 = get_closest_state_time_to(ens_time, next_time)

! WATCH OUT FOR USING BOUNDARY TIME OBS TWICE; add one second to bottom time
! ALSO, avoid having a negative time for the start (for low-order models in general)
if(delta_time / 2 > time2) then
   start_time = set_time(0, 0)
else
   start_time = time2 - delta_time / 2 + set_time(1, 0)
endif
end_time = time2 + delta_time / 2

! Output very brief current start and end time at the message level
if (print_trace_details == 0) then
   call timechat(start_time,  'move_ahead', .false.,        'Next assimilation window starts    at: ')
   call timechat(end_time,    'move_ahead', .false.,        'Next assimilation window ends      at: ')
endif

! This block of code gets called either if the next obs is not in the current window,
! or if you're asking for the details of the assimilation window and obs times.
if(next_time < start_time .or. next_time > end_time .or. print_trace_details > 0) then

   if (time2 /= ens_time) then
      !call timechat(next_time,   'move_ahead', .false.,     'Next available observation time    is: ')
      call timechat(time2,       'move_ahead', .false.,     'Next data time should be           at: ', &
         'Not within current window, model will be called to advance state.')
      call timechat(start_time,  'move_ahead', .false.,     'Next assimilation window starts    at: ')
      call timechat(end_time,    'move_ahead', .false.,     'Next assimilation window ends      at: ')
   else 
      if (next_time >= start_time .and. next_time <= end_time) then
         call timechat(next_time,   'move_ahead', .false.,  'Next available observation time    is: ', &
            'Within current assimilation window, model does not need advance.')
      else 
         call timechat(next_time,   'move_ahead', .false.,  'Next available observation time    is: ', &
            'Next obs outside current assimilation window.')
      endif
   endif

   ! if different mpi tasks have different times, the default is only process 0
   ! will print messages.  in this case we're headed towards a fatal error and
   ! just trying to give the most info possible before exiting.  so make all
   ! mpi tasks which get into this block output.  in the worst case you'll get
   ! N sets of these messages which is messy, but probably better than having
   ! the case where process 0 works but some other tasks fail and you get no
   ! helpful info from them.
   if (next_time < start_time .or. next_time > end_time) then
      call set_output(.true.) 
      call error_handler(E_MSG, ' ', ' ')
      call error_handler(E_MSG, 'move_ahead', 'Inconsistent model state/observation times. ')

      if (next_time < start_time) then
         call error_handler(E_MSG, 'move_ahead', &
            'Next observation cannot be earlier than start of new time window')
         call error_handler(E_MSG, ' ', ' ')
         call timechat(next_time,   'move_ahead', .false.,  'Next available observation         at: ')
         call timechat(start_time,  'move_ahead', .false.,  'Next assimilation window starts    at: ')
         call timechat(ens_time,    'move_ahead', .false.,  'Current model data time            is: ')
         errstring1 = 'If this is the start of the obs_seq file, '
         errstring2 = 'can use filter namelist to set first obs or initial data time.'
      else
         call error_handler(E_MSG, 'move_ahead', &
            'Next observation is later than end of new time window')
         call error_handler(E_MSG, ' ', ' ')
         call timechat(next_time,   'move_ahead', .false.,  'Next available observation         at: ')
         call timechat(end_time,    'move_ahead', .false.,  'Next assimilation window ends      at: ')
         call timechat(ens_time,    'move_ahead', .false.,  'Current model data time            is: ')
         errstring1 = 'should not happen; code has miscomputed how far to advance'
         errstring2 = ''
      endif
      call error_handler(E_MSG, ' ', ' ')
      call error_handler(E_ERR, 'move_ahead', &
        'Inconsistent model state/observation times, cannot continue', &
         source, text2=errstring1, text3=errstring2)
   endif
endif

! WANT BETTER WINDOW CONTROL, TIME INTERPOLATION, TOO.
! Get all the observations that are in the observation window - do this with the
! window times, but set the returns to the data times.
call get_obs_time_range(seq, start_time, end_time, key_bounds, num_obs_in_set, &
   out_of_range, observation)

! ok, not really a time detail, but if turned off, the output is pretty much
! the same as the original.
if (print_trace_details > 0) then
   write (errstring, '(A,I8,A)') 'Next assimilation window contains up to ', &
                                  num_obs_in_set, ' observations'
   call error_handler(E_MSG, 'move_ahead', errstring)
endif

curr_ens_time = ens_time
next_ens_time = time2

! Release the storage associated with the observation temp variable
call destroy_obs(observation)

if (ens_handle%my_pe == 0 .and. my_task_id() /= 0) then
  call set_output(.false.)
endif

end subroutine move_ahead


!------------------------------------------------------------------------------

subroutine advance_state(ens_handle, ens_size, target_time, async, adv_ens_command, &
                         tasks_per_model_advance, file_info_output, file_info_input)

! Advances all ensemble copies of the state to the target_time.  Note that 
! there may be more than ens_size copies in the ens_handle storage. Copies other
! than the ensemble copies need not be advanced in time.

implicit none

! ens_size is the number of ensembles that are model state and need to be advanced
type(ensemble_type),  intent(inout) :: ens_handle
type(time_type),      intent(in)    :: target_time
integer,              intent(in)    :: ens_size, async
character(len=*),     intent(in)    :: adv_ens_command
integer,              intent(in)    :: tasks_per_model_advance
type(file_info_type), intent(inout) :: file_info_output
type(file_info_type), intent(inout) :: file_info_input

character(len = 129), dimension(ens_handle%num_copies) :: ic_file_name, ud_file_name 
character(len = 129)                                   :: control_file_name
! WARNING: ARE THESE LENGTHS OKAY?
character(len = 129)                                   :: system_string

type(time_type) :: time_step, ens_time
integer         :: is1, is2, id1, id2, my_num_state_copies, global_ens_index
integer         :: i, control_unit, rc
integer         :: need_advance, any_need_advance
logical         :: read_time_from_file = .true.

! Initialize if needed
if(.not. module_initialized) then
   call initialize_module
   module_initialized = .true.
endif

! Tasks without ensemble members do not know the ensemble time, so assume
! no advance until we can read the time on at least some tasks and confirm
! the advance is needed (which it should be).
need_advance = 0

! Determine model time_step
time_step = get_model_time_step()

! Loop through each model state and advance; Count number of copies that are state ensembles
my_num_state_copies = 0
ENSEMBLE_MEMBERS: do i = 1, ens_handle%my_num_copies
   global_ens_index = ens_handle%my_copies(i)

   ! If global_ens_index is not one of the model state ensembles, don't need to advance
   ! Assumes that state ensemble copies have smallest global index and are stored 
   ! contiguously on local processor
   if(global_ens_index > ens_size) exit ENSEMBLE_MEMBERS

   ! No need to advance if already at target time 
   if(ens_handle%time(i) == target_time) exit ENSEMBLE_MEMBERS

   ! Check for time error
   if(ens_handle%time(i) > target_time) then
      call get_time(ens_handle%time(i), is1, id1)
      call get_time(target_time, is2, id2)
      write(errstring,*)'target time ',is2,id2,' is before model_time ',is1,id1
      call error_handler(E_ERR,'advance_state', errstring, source)
   endif

   ! Ok, this task does need to advance something. 
   need_advance = 1

   ! Increment number of ensemble member copies I have.
   my_num_state_copies = my_num_state_copies + 1

   !------------- Block for subroutine callable adv_1step interface -----------
   if(async == 0) then

      do while(ens_handle%time(i) < target_time)
         call adv_1step(ens_handle%vars(:, i), ens_handle%time(i))
         ens_handle%time(i) = ens_handle%time(i) + time_step
      end do

   !-------------- End single subroutine callable adv_1step interface ---------
   else
   !-------------- Block for calling shell to advance model -------------------

      ! Can only handle up to 10000 ensemble members
      if(global_ens_index >= 10000) then
         write(errstring,*)'Trying to use ', ens_size,' model states -- too many.'
         call error_handler(E_WARN,'advance_state',errstring,source)
         call error_handler(E_ERR,'advance_state','Use less than 10000 member ensemble.', &
            source)
      endif

      !>@todo FIXME write the netCDF restart names into the control file
      ! Create file names for input and output state files for this copy
!     write(ic_file_name(i), '("assim_model_state_ic.", i4.4)') global_ens_index
!     write(ud_file_name(i), '("assim_model_state_ud.", i4.4)') global_ens_index

      ! Open a restart file and write state following a target time
      call write_state(ens_handle, file_info_output)

   endif

   !-------------- End of block for calling shell to advance model -------------

end do ENSEMBLE_MEMBERS


! Following is for async options that use shell to advance model
SHELL_ADVANCE_METHODS: if(async /= 0) then

   ! If no one needs advance, get out now.  This is a global communication routine.
   call sum_across_tasks(need_advance, any_need_advance)
   if (any_need_advance == 0) then
      call error_handler(E_MSG, 'advance_state', 'Model time already at requested time')
      return
   endif

   ! Get a unique name for the control file; use process id
   if(my_task_id() >= 10000) call error_handler(E_ERR, 'advance_state', &
      'Can only have 10000 processes', source)
   write(control_file_name, '("filter_control", i5.5)') my_task_id()

   if (async == 2) then

      if (have_members(ens_handle, ens_size)) then
         ! Everyone writes the ic and ud file names to the control file
         control_unit = get_unit()
         open(unit = control_unit, file = trim(control_file_name))
         ! Write out the file names to a control file
         ! Also write out the global ensemble member number for the script 
         ! if it needs to create unique filenames for its own use.
         do i = 1, my_num_state_copies
            write(control_unit, '(i5)') ens_handle%my_copies(i)
            write(control_unit, '(a)' ) trim(ic_file_name(i)) !!@todo netCDF
            write(control_unit, '(a)' ) trim(ud_file_name(i)) !!@todo netCDF
         end do
         close(control_unit)
   
         ! Arguments to advance model script are unique id and number of copies
         write(system_string, '(i10, 1x, i10)') my_task_id(), my_num_state_copies
        
         if (debug) write(*,*) 'iam ', my_task_id(), ' ready to execute shell command'

         ! Issue a system command with arguments my_task, my_num_copies, and control_file
         rc = shell_execute(trim(adv_ens_command)//' '//trim(system_string)//' '//trim(control_file_name))
         if (debug) write(*,*) 'iam ', my_task_id(), ' shell execute returns ', rc
         if (debug) write(*,*) 'iam ', my_task_id(), ' checking existance of file ', trim(control_file_name)
      
         ! if control file is still here, the advance failed
         if(file_exist(control_file_name)) then
           write(errstring, "(A)")'If advance script finishes ok it removes '//trim(control_file_name)
           write(errstring1,"(A)")'It still exists, so 1+ members listed in that file failed to run'
           write(errstring2,"(A)")'Check the output of the model or script to find the error.'
           call error_handler(E_ERR,'advance_state', errstring, source, &
                              text2=errstring1, text3=errstring2)
         endif
      endif

   else if (async == 4) then

      ! Only task 0 writes all ens names to file.
      if (my_task_id() == 0) then
         ! Write the ic and ud file names to the control file
         control_unit = get_unit()
         open(unit = control_unit, file = trim(control_file_name))
         ! Write out the file names to a control file
         ! Also write out the global ensemble member number for the script 
         ! if it needs to create unique filenames for its own use.
         do i = 1, ens_size
            write(control_unit, '(i5)') i
            write(control_unit, '("assim_model_state_ic.", i4.4)') i
            write(control_unit, '("assim_model_state_ud.", i4.4)') i
         end do
         close(control_unit)
      endif

      ! make sure all tasks have finished writing their initial condition files
      ! (in the previous block above) before letting proc 0 start the model advances.
      call task_sync()

      ! PE0 tells the run script to all the model advance script here, then sleeps.  
      ! the other tasks just sleep.  TODO: this code is too intertwined with the
      ! run script - it assumes that the script has already created 2 fifo/pipe 
      ! files with specific names:  mkfifo filter_to_model.lock model_to_filter.lock
      ! eventually it should probably be yet another separate executable which
      ! keeps the code all together in the mpi module.

      ! all tasks now get here, even those which do not have ensemble members,
      ! so this bookkeeping is not needed, i think.
      !! tell any waiting tasks that we are indeed doing a model advance
      !! here and they need to block.
      !call sum_across_tasks(1, i)

      ! block, not spin, until the model advance is done.
      ! the code does not resume running until the separate program 'wakeup_filter'
      ! calls the companion 'resume_task()'.
      call block_task()


      ! if control file is still here, the advance failed
      if(file_exist(control_file_name)) then
        write(errstring, "(A)")'If advance script finishes ok it removes '//trim(control_file_name)
        write(errstring1,"(A)")'It still exists, so 1+ members listed in that file failed to run'
        write(errstring2,"(A)")'Check the output of the model or script to find the error.'
        call error_handler(E_ERR,'advance_state', errstring, source, &
                           text2=errstring1, text3=errstring2)
      endif

   else
      ! Unsupported option for async error
      write(errstring,*)'input.nml - async is ',async,' must be 0, 2, or 4'
      call error_handler(E_ERR,'advance_state', errstring, source)
   endif

   ! All should be done, read in the states and proceed
   ! NOTE: if async 2 and we add support for tasks_per_model_advance, then we need
   ! another sync here before trying to read in from a task which might not have been
   ! the one which advanced it.
   call read_state(ens_handle, file_info_input, read_time_from_file, ens_time)

end if SHELL_ADVANCE_METHODS


! any process that owns members, check that they advanced ok.
if (have_members(ens_handle, ens_size)) then

   do i = 1, my_num_state_copies

      ! Get times of all ensemble members, one by one, to catch errors
      ! when some advanced ok but not all even if they appeared to.
      call get_ensemble_time(ens_handle, i, ens_time)
      
      ! error out if model state did not advance to when requested.
      if (ens_time /= target_time) then
         ! make sure times print out on whatever task this is (it uses msgs)
         call set_output(.true.)
         call timechat(ens_time,    'advance_state', .true.,   'Model advance time is now:             ')
         call timechat(target_time, 'advance_state', .false.,  'Model advance time NOT what requested: ')
       
         write(errstring2, *) 'in timestamp on ensemble member ', ens_handle%my_copies(i)
         call error_handler(E_ERR,  'advance_state', 'Model advance complete but model time not correct', &
                            source, text2=errstring2)
      endif
   enddo

   if (print_trace_details > 0) then
      call error_handler(E_MSG, 'advance_state', 'Model advance complete, model time updated to requested time')
   endif

endif


end subroutine advance_state

!--------------------------------------------------------------------
! subroutine wait_if_needed(async) ! REMOVED Wed Feb 17 09:17:44 MST 2016 r9776ish
!--------------------------------------------------------------------

!--------------------------------------------------------------------

subroutine timechat(a_time, label, blank, string1, string2, string3)
 type(time_type),  intent(in)           :: a_time
 character(len=*), intent(in)           :: label, string1
 logical,          intent(in)           :: blank
 character(len=*), intent(in), optional :: string2, string3

! prettyprint a time, with an optional preceeding blank line, a label
! before the time, and up to 2 additional text lines afterwards.
 
integer :: a_day, a_sec

call get_time(a_time, a_sec, a_day)
write (errstring, '(A,A,I8,A,I6)') string1, ' day=', a_day, ' sec=', a_sec

if (blank) call error_handler(E_MSG, ' ', ' ')
call error_handler(E_MSG, label, trim(errstring))

if (present(string2)) then
   call error_handler(E_MSG, label, string2)
endif

if (present(string3)) then
   call error_handler(E_MSG, label, string3)
endif

end subroutine timechat

!--------------------------------------------------------------------

subroutine set_obs_model_trace(execution_level, timestamp_level)
 integer, intent(in) :: execution_level
 integer, intent(in) :: timestamp_level

! set module local vars from the calling code to indicate how much
! output we should generate from this code.  execution level is
! intended to make it easier to figure out where in the code a crash
! is happening; timestamp level is intended to help with gross levels
! of overall performance profiling.  eventually, a level of 1 will
! print out only basic info; level 2 will be more detailed.  
! (right now, only > 0 prints anything and it doesn't matter how 
! large the value is.)
 

print_trace_details = execution_level
print_timestamps    = timestamp_level

end subroutine set_obs_model_trace

!--------------------------------------------------------------------

function have_members(ens_handle, ens_size)

! return true if this task contains ensemble copies.
! mean, sd, etc do not count. 
!
! note this depends on how filter (or the calling code) has laid out
! the actual ensemble members vs the other state-length arrays
! carried around in the ensemble handle.

type(ensemble_type), intent(in) :: ens_handle
integer, intent(in)             :: ens_size
logical                         :: have_members

integer :: my_first_copy

! assumes copies are stored in the lowest numbered slots.
if(ens_handle%my_num_copies >= 1) then
   my_first_copy = ens_handle%my_copies(1)
else
   my_first_copy = 0
endif

! if i have at least one copy, and it is a real member, 
! not aux data (mean, inflation, etc), return true.
have_members = ((ens_handle%my_num_copies >= 1) .and. (my_first_copy <= ens_size)) 

end function have_members

!--------------------------------------------------------------------

end module obs_model_mod

