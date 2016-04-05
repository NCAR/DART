! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_model_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

use        types_mod, only : r8
use    utilities_mod, only : register_module, error_handler, E_ERR
use     location_mod, only : location_type
use  assim_model_mod, only : aget_closest_state_time_to, &
                             am_get_close_states => get_close_states, &
                             get_model_time_step
use obs_sequence_mod, only : obs_sequence_type, obs_type, get_obs_from_key, &
                             get_obs_def, init_obs, destroy_obs, get_num_copies, &
                             get_num_qc, get_first_obs, get_next_obs, get_obs_time_range
use      obs_def_mod, only : obs_def_type, get_obs_def_location, get_obs_def_time
use time_manager_mod, only : time_type, operator(/=), operator(>), set_time, &
                             operator(-), operator(/), operator(+), print_time, operator(<)

use ensemble_manager_mod, only : get_ensemble_time, Aadvance_state, ensemble_type

implicit none
private

public :: get_close_states, move_ahead

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"


logical, save :: module_initialized = .false.

contains

!======================================================================

subroutine initialize_module

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module



subroutine get_close_states(seq, key, radius, numinds, indices, dist, x)
!------------------------------------------------------------------------
!

type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: key
real(r8), intent(in) :: radius
integer, intent(out) :: numinds, indices(:)
real(r8), intent(out) :: dist(:)
real(r8), intent(in)  :: x(:)

type(obs_type) :: obs
type(obs_def_type) :: obs_def
type(location_type) :: location

if ( .not. module_initialized ) call initialize_module

! Initialize the observation type
!!! Can actually init with the correct size here if wanted
call init_obs(obs, 0, 0)

! Get the location of the obs and pass through to assim_model interface
call get_obs_from_key(seq, key, obs)
call get_obs_def(obs, obs_def)
location = get_obs_def_location(obs_def)

call am_get_close_states(location, radius, numinds, indices, dist, x)

!Need to destroy this obs to free up allocated space
call destroy_obs(obs)

end subroutine get_close_states


!-------------------------------------------------------------------------

subroutine move_ahead(ens_handle, ens_size, model_size, seq, last_key_used, &
   key_bounds, num_obs_in_set, async, adv_ens_command)

! First version of this assumes that observations come in discrete chunks that are exactly
! at times to which the model can advance as per historical filters. Variants should be
! added and routine should be moved to a separate module, probably obs_model_mod.

implicit none

type(ensemble_type),     intent(inout) :: ens_handle
integer,                 intent(in)  :: ens_size, model_size
type(obs_sequence_type), intent(in)  :: seq
integer,                 intent(in)  :: last_key_used, async
integer,                 intent(out) :: key_bounds(2), num_obs_in_set
character(len = 129),    intent(in)  :: adv_ens_command

type(time_type) :: next_time, time2, start_time, end_time, delta_time, ens_time
type(obs_type)  :: observation
type(obs_def_type) :: obs_def
logical :: is_this_last, is_there_one, out_of_range

! Initialize a temporary observation type to use
call init_obs(observation, get_num_copies(seq), get_num_qc(seq))

! Get the next observation in the sequence that hasn't already been assimilated
if(last_key_used > 0) then
   call get_obs_from_key(seq, last_key_used, observation)
   call get_next_obs(seq, observation, observation, is_this_last)
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

!write(*, *) 'next_time is time of next observation'
!call print_time(next_time)

! Get the time of the ensemble, assume consistent across all
call get_ensemble_time(ens_handle, 1, ens_time)

! Figure out time to which to advance model 
time2 = aget_closest_state_time_to(ens_time, next_time)

! Compute the model time step and center a window around the closest time
delta_time = get_model_time_step()
! WATCH OUT FOR USING BOUNDARY TIME OBS TWICE; add one second to bottom time
! ALSO, avoid having a negative time for the start for low-order models
if(delta_time / 2 > time2) then
   start_time = set_time(0, 0)
else
   start_time = time2 - delta_time / 2 + set_time(1, 0)
endif
end_time = time2 + delta_time / 2

call print_time(start_time,str=' start time of obs range')
call print_time(end_time,  str='   end time of obs range')

! If the next observation is not in the window, then have an error
if(next_time < start_time .or. next_time > end_time) then
   call error_handler(E_ERR, 'move_ahead', 'next obs time not in model time window', &
      source, revision, revdate)
endif

! Get all the observations at exactly this time
call get_obs_time_range(seq, start_time, end_time, key_bounds, num_obs_in_set, &
   out_of_range, observation)

! Advance all ensembles (to the time of the first ensemble)
if(time2 /= ens_time) call Aadvance_state(ens_handle, time2, async, adv_ens_command)

! Release the storage associated with the observation temp variable
call destroy_obs(observation)

end subroutine move_ahead

end module obs_model_mod
