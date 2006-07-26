! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_model_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

use types_mod,            only : r8
use utilities_mod,        only : register_module, error_handler, E_ERR, E_MSG, E_WARN, &
                                 logfileunit, get_unit, file_exist
use assim_model_mod,      only : aget_closest_state_time_to, get_model_time_step, &
                                 open_restart_write, open_restart_read,           &
                                 awrite_state_restart, close_restart, adv_1step,  &
                                 aread_state_restart
use obs_sequence_mod,     only : obs_sequence_type, obs_type, get_obs_from_key,      &
                                 get_obs_def, init_obs, destroy_obs, get_num_copies, &
                                 get_num_qc, get_first_obs, get_next_obs, get_obs_time_range
use obs_def_mod,          only : obs_def_type, get_obs_def_time
use time_manager_mod,     only : time_type, operator(/=), operator(>), set_time, get_time, &
                                 operator(-), operator(/), operator(+), print_time,        &
                                 operator(<), operator(==)
use ensemble_manager_mod, only : get_ensemble_time, ensemble_type
use mpi_utilities_mod,    only : my_task_id

implicit none
private

public :: move_ahead, advance_state

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

logical :: module_initialized = .false.

! Module storage for writing error messages
character(len = 129) :: errstring

contains

!======================================================================

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!-------------------------------------------------------------------------

subroutine move_ahead(ens_handle, ens_size, seq, last_key_used, &
   key_bounds, num_obs_in_set, async, adv_ens_command)

! Advances all ensemble copies of the state to the time appropriate for
! the next unused observation in the observation sequence. Note that there
! may be more than ens_size copies in the ens_handle storage. Copies other
! than the ensemble copies need not be advanced in time.

implicit none

type(ensemble_type),     intent(inout) :: ens_handle
integer,                 intent(in)    :: ens_size
type(obs_sequence_type), intent(in)    :: seq
integer,                 intent(in)    :: last_key_used, async
integer,                 intent(out)   :: key_bounds(2), num_obs_in_set
character(len = *),      intent(in)    :: adv_ens_command

type(time_type)    :: next_time, time2, start_time, end_time, delta_time, ens_time
type(obs_type)     :: observation
type(obs_def_type) :: obs_def
logical            :: is_this_last, is_there_one, out_of_range, i
integer            :: sec, day

! Initialize if needed
if(.not. module_initialized) then
   call initialize_module
   module_initialized = .true.
endif

! Violating the private boundaries on ens_handle by direct access
! If none of my copies are regular ensemble members no need to advance
if(ens_handle%my_num_copies < 1) return
! Next assumes that ensemble copies are in the first ens_size global copies of ensemble
! And that global indices are monotonically increasing in each pes local copies
if(ens_handle%my_copies(1) > ens_size) return

! Initialize a temporary observation type to use
call init_obs(observation, get_num_copies(seq), get_num_qc(seq))

! Get the next observation in the sequence that hasn't already been assimilated
if(last_key_used > 0) then
   call get_obs_from_key(seq, last_key_used, observation)
   call get_next_obs(seq, observation, observation, is_this_last)
   ! If the end of the observation sequence has been passed, return
   if(is_this_last) then
      key_bounds(1:2) = -99
      return
   endif
else
   is_there_one = get_first_obs(seq, observation)
   if(.not. is_there_one) then
      key_bounds(1:2) = -99
      return
   endif
endif

! Get the time of this observation
call get_obs_def(observation, obs_def)
next_time = get_obs_def_time(obs_def)

! Get the time of the ensemble, assume consistent across all ensemble copies
call get_ensemble_time(ens_handle, 1, ens_time)

! Figure out time to which to advance model 
! More control over time window use of observations would come in here
time2 = aget_closest_state_time_to(ens_time, next_time)

! Compute the model time step and center a window around the closest time
delta_time = get_model_time_step()
! WATCH OUT FOR USING BOUNDARY TIME OBS TWICE; add one second to bottom time
! ALSO, avoid having a negative time for the start (for low-order models in general)
if(delta_time / 2 > time2) then
   start_time = set_time(0, 0)
else
   start_time = time2 - delta_time / 2 + set_time(1, 0)
endif
end_time = time2 + delta_time / 2

! Output the start and end time at the message level
call get_time(start_time, sec, day)
write (errstring, *) 'Start time of obs range day=', day, ', sec=', sec
call error_handler(E_MSG, 'move_ahead', errstring, '', '', '')
call get_time(end_time,   sec, day)
write (errstring, *) 'End time of obs range day=  ', day, ', sec=', sec
call error_handler(E_MSG, 'move_ahead', errstring, source, revision, revdate)

! If the next observation is not in the window, then have an error
if(next_time < start_time .or. next_time > end_time) then
   ! Is this test still really needed?
   write(errstring, *) 'next obs time not in model time window: day=', day, ', sec=', sec
   call error_handler(E_ERR, 'move_ahead', errstring, source, revision, revdate)
endif

! WANT BETTER WINDOW CONTROL, TIME INTERPOLATION, TOO.
! Get all the observations that are in the observation window
call get_obs_time_range(seq, start_time, end_time, key_bounds, num_obs_in_set, &
   out_of_range, observation)

! Advance all ensembles to the time for the assimilation
if(time2 /= ens_time) call advance_state(ens_handle, ens_size, time2, async, adv_ens_command)

! Release the storage associated with the observation temp variable
call destroy_obs(observation)

end subroutine move_ahead


!------------------------------------------------------------------------------

subroutine advance_state(ens_handle, ens_size, target_time, async, adv_ens_command) 

! Advances the model extended state until time is equal to the target_time.

implicit none

! ens_size is the number of ensembles that are model state and need to be advanced
type(ensemble_type), intent(inout) :: ens_handle
type(time_type),     intent(in)    :: target_time
integer,             intent(in)    :: ens_size, async
character(len=*),    intent(in)    :: adv_ens_command

character(len = 129), dimension(ens_handle%num_copies) :: ic_file_name, ud_file_name 
character(len = 129)                                   :: control_file_name
! WARNING: ARE THESE LENGTHS OKAY?
character(len = 129)                                   :: system_string

type(time_type) :: time_step
integer         :: is1, is2, id1, id2, my_num_state_copies, global_ens_index
integer         :: i, control_unit, ic_file_unit, ud_file_unit

! Initialize if needed
if(.not. module_initialized) then
   call initialize_module
   module_initialized = .true.
endif

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

   ! Increment number of ensemble member copies I have
   my_num_state_copies = my_num_state_copies + 1

   ! If no need to advance return
   if(ens_handle%time(i) == target_time) return

   ! Check for time error
   if(ens_handle%time(i) > target_time) then
      call get_time(ens_handle%time(i), is1, id1)
      call get_time(target_time, is2, id2)
      write(errstring,*)'target time ',is2,id2,' is before model_time ',is1,id1
      call error_handler(E_ERR,'advance_state', errstring, source, revision, revdate)
   endif

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
      if(global_ens_index > 10000) then
         write(errstring,*)'Trying to use ', ens_size,' model states -- too many.'
         call error_handler(E_WARN,'advance_state',errstring,source,revision,revdate)
         call error_handler(E_ERR,'advance_state','Use less than 10000 member ensemble.', &
            source,revision,revdate)
      endif

      ! Create file names for input and output state files for this copy
      write(ic_file_name(i), '("assim_model_state_ic", i5.5)') global_ens_index
      write(ud_file_name(i), '("assim_model_state_ud", i5.5)') global_ens_index


      ! Open a restart file and write state following a target time
      ! Force writing to binary to have bit-wise reproducing advances
      ic_file_unit = open_restart_write(trim(ic_file_name(i)), "unformatted")
      call awrite_state_restart(ens_handle%time(i), ens_handle%vars(:, i), ic_file_unit, target_time)
      call close_restart(ic_file_unit)

   endif

   !-------------- End of block for calling shell to advance model -------------

end do ENSEMBLE_MEMBERS


! Following block needed for WRF only : PLEASE CONFIRM IF THEY STILL NEED THIS
! Write out the ensemble mean for the calculation of flow-dependent BC of WRF.
!call update_ens_mean(ens_handle)
!ic_file_unit = open_restart_write('assim_model_state_ic_mean')
!call awrite_state_restart(smodel_time, ens_mean, ic_file_unit, target_time)
!call close_restart(ic_file_unit)


! Following is for async options that use shell to advance model
SHELL_ADVANCE_METHODS: if(async /= 0) then
   ! Get a unique name for the control file; use process id
   if(my_task_id() > 10000) call error_handler(E_ERR, 'advance_state', &
      'Can only have 10000 processes', source, revision, revdate)
   write(control_file_name, '("filter_control", i5.5)') my_task_id()


   ! Write the ic and ud file names to the control file
   control_unit = get_unit()
   open(unit = control_unit, file = trim(control_file_name))
   ! Write out the file names to a control file
   do i = 1, my_num_state_copies
      write(control_unit, '(a)' ) trim(ic_file_name(i))
      write(control_unit, '(a)' ) trim(ud_file_name(i))
   end do
   close(control_unit)

   if(async == 2) then

      ! Arguments to advance model script are unique id and number of copies
      write(system_string, '(i10, 1x, i10)') my_task_id(), my_num_state_copies
     
      ! Issue a system command with arguments my_task, my_num_copies, and control_file
      call system(trim(adv_ens_command)//' '//trim(system_string)//' '//trim(control_file_name))
 
      FOREVER: do
         ! When control file is deleted by shell script (adv_ens_command) 
         ! we are free to proceed to read in updated state
         if(file_exist(control_file_name)) then
            ! Adjusting sleep time can modify performance
            ! Use shorter time for shorter expected model advances
            ! Might want to add this to namelist if it's important ???????
            call system ('sleep 0.1')
         else
            exit FOREVER
         endif
      end do FOREVER

   else
   ! Unsupported option for async error
      write(errstring,*)'input.nml - async is ',async,' must be 0, or 2'
      call error_handler(E_ERR,'advance_state', errstring, source, revision, revdate)
   endif

   ! All should be done, read in the states and proceed
   do i = 1, my_num_state_copies
      ud_file_unit = open_restart_read(trim(ud_file_name(i)))
      call aread_state_restart(ens_handle%time(i), ens_handle%vars(:, i), ud_file_unit)
      call close_restart(ud_file_unit)
   end do

end if SHELL_ADVANCE_METHODS

end subroutine advance_state

!--------------------------------------------------------------------

end module obs_model_mod
