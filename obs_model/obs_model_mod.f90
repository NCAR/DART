! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_model_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

use       types_mod,  only : r8
use   utilities_mod,  only : register_module
use    location_mod,  only : location_type, interactive_location
use assim_model_mod,  only : interpolate, aget_closest_state_time_to, &
                             am_get_close_states => get_close_states, &
                             aadvance_state
use    obs_kind_mod,  only : obs_kind_type, interactive_kind, get_obs_kind
use obs_sequence_mod, only : obs_sequence_type, obs_type, get_obs_from_key, &
                             get_obs_def, init_obs, destroy_obs, get_num_copies, &
                             get_num_qc, get_first_obs, get_next_obs, get_obs_time_range
use obs_def_mod,      only : obs_def_type, get_obs_def_location, get_obs_def_kind, &
                             get_obs_def_time
use utilities_mod,    only : error_handler, E_ERR
use time_manager_mod, only : time_type, operator(/=), operator(>), get_time


implicit none
private

public take_obs, interactive_def, get_expected_obs, get_close_states, move_ahead

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





subroutine take_obs(state_vector, location, obs_kind, obs_vals, istatus, rstatus)
!--------------------------------------------------------------------
!

implicit none

real(r8),            intent(in) :: state_vector(:)
type(location_type), intent(in) :: location
type(obs_kind_type), intent(in) :: obs_kind
real(r8),           intent(out) :: obs_vals
integer,  optional, intent(out) :: istatus
real(r8), optional, intent(out) :: rstatus

if ( .not. module_initialized ) call initialize_module

! Initially, have only raw state observations implemented so obs_kind is
! irrelevant. Just do interpolation and return.

if(present(rstatus)) then
   call interpolate(state_vector, location, get_obs_kind(obs_kind), obs_vals, istatus, rstatus)
elseif(present(istatus)) then
   call interpolate(state_vector, location, get_obs_kind(obs_kind), obs_vals, istatus)
else
   call interpolate(state_vector, location, get_obs_kind(obs_kind), obs_vals)
endif

end subroutine take_obs



subroutine interactive_def(location, obs_kind)
!---------------------------------------------------------------------
!
! Used for interactive creation of an obs_def. For now, an observation
! is uniquely defined by a location and a kind. The obs_model level
! knows how kinds and locations need to be combined if there are any
! intricacies (although for models to date there are none, 17 April, 2002).

implicit none

type(location_type), intent(out) :: location
type(obs_kind_type), intent(out) :: obs_kind

if ( .not. module_initialized ) call initialize_module

call interactive_location(location)
call interactive_kind(obs_kind)

end subroutine interactive_def


subroutine get_expected_obs(seq, keys, state, obs_vals, istatus, rstatus)
!---------------------------------------------------------------------------
! Compute forward operator for set of obs in sequence

type(obs_sequence_type), intent(in) :: seq
integer,                 intent(in) :: keys(:)
real(r8),                intent(in) :: state(:)
real(r8),               intent(out) :: obs_vals(:)
integer,  optional,     intent(out) :: istatus
real(r8), optional,     intent(out) :: rstatus(:,:)

integer             :: num_obs, i
type(location_type) :: location
type(obs_kind_type) :: obs_kind
type(obs_type)      :: obs
type(obs_def_type)  :: obs_def
integer             :: obs_kind_ind

if ( .not. module_initialized ) call initialize_module

num_obs = size(keys)

! Initialize the observation type
!!! Can actually init with the correct size here if wanted
call init_obs(obs, 0, 0)

do i = 1, num_obs
   call get_obs_from_key(seq, keys(i), obs)
   call get_obs_def(obs, obs_def)
   location = get_obs_def_location(obs_def)
   obs_kind = get_obs_def_kind(obs_def)
! Check in kind for negative for identity obs
   obs_kind_ind = get_obs_kind(obs_kind)
   if(obs_kind_ind < 0) then
      if ( -obs_kind_ind > size(state) ) call error_handler(E_ERR, &
         'get_expected_obs', &
         'identity obs is outside of state vector ', &
         source, revision, revdate)
      obs_vals(i) = state(-1 * obs_kind_ind)
      if(present(rstatus)) rstatus(i,1) = 0.0_r8
   else
! Otherwise do forward operator
      if(present(rstatus)) then
         call take_obs(state, location, obs_kind, obs_vals(i), istatus, rstatus(i,1))
      elseif(present(istatus)) then
         call take_obs(state, location, obs_kind, obs_vals(i), istatus)
      else
         call take_obs(state, location, obs_kind, obs_vals(i))
      endif
   endif
end do

! Need to destroy this obs to free up allocated space
call destroy_obs(obs)

end subroutine get_expected_obs



subroutine get_close_states(seq, key, radius, numinds, indices, dist)
!------------------------------------------------------------------------
!

type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: key
real(r8), intent(in) :: radius
integer, intent(out) :: numinds, indices(:)
real(r8), intent(out) :: dist(:)

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

call am_get_close_states(location, radius, numinds, indices, dist)

!Need to destroy this obs to free up allocated space
call destroy_obs(obs)

end subroutine get_close_states

!-------------------------------------------------------------------------

subroutine move_ahead(ens, ens_time, ens_size, model_size, seq, last_key_used, &
   key_bounds, num_obs_in_set, async, adv_ens_command)

! First version of this assumes that observations come in discrete chunks that are exactly
! at times to which the model can advance as per historical filters. Variants should be
! added and routine should be moved to a separate module, probably obs_model_mod.

implicit none

integer, intent(in) :: ens_size, model_size
real(r8), intent(inout) :: ens(ens_size, model_size)
type(time_type), intent(inout) :: ens_time(ens_size)
type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: last_key_used, async
integer, intent(out) :: key_bounds(2), num_obs_in_set
character(len = 129), intent(in) :: adv_ens_command

type(time_type) :: next_time, time2
type(obs_type) :: observation
type(obs_def_type) :: obs_def
logical :: is_this_last, is_there_one, out_of_range
integer :: days, secs

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

! Get all the obsevations at exactly this time
call get_obs_time_range(seq, next_time, next_time, key_bounds, num_obs_in_set, out_of_range, observation)

call get_time(next_time, secs, days)
write(*, *) 'time of obs set  is (d, s) = ', days, secs

! If the model time is past the obs set time, should be an error
if(ens_time(1) > next_time) then
   write(*, *) 'model time is already past time of obs_set in move_ahead in filter'
   stop
endif

! Figure out time to which to advance model 
time2 = aget_closest_state_time_to(ens_time(1), next_time)

! Advance all ensembles (to the time of the first ensemble)
if(time2 /= ens_time(1)) call Aadvance_state(ens_time, ens, ens_size, time2, async, adv_ens_command)

! Release the storage associated with the observation temp varialbe
call destroy_obs(observation)

end subroutine move_ahead

end module obs_model_mod
