! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_model_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

use        types_mod, only : r8, t_kelvin, es_alpha, es_beta, es_gamma, gas_constant_v, &
                             gas_constant, L_over_Rv, missing_r8, ps0
use    utilities_mod, only : register_module, error_handler, E_ERR
use     location_mod, only : location_type, interactive_location, get_location
use  assim_model_mod, only : interpolate, aget_closest_state_time_to, &
                             am_get_close_states => get_close_states, &
                             Aadvance_state, get_model_time_step
use     obs_kind_mod, only : obs_kind_type, interactive_kind, get_obs_kind, &
                             KIND_U, KIND_V, KIND_PS, KIND_T, KIND_QV, &
                             KIND_P, KIND_W, KIND_QR, KIND_TD, KIND_VR, &
                             KIND_REF
use obs_sequence_mod, only : obs_sequence_type, obs_type, get_obs_from_key, &
                             get_obs_def, init_obs, destroy_obs, get_num_copies, &
                             get_num_qc, get_first_obs, get_next_obs, get_obs_time_range
use      obs_def_mod, only : obs_def_type, get_obs_def_location, get_obs_def_kind, &
                             get_obs_def_time
!WRF use      obs_def_mod, only : get_obs_def_platform
!WRF use     platform_mod, only : platform_type, get_platform_orientation
use time_manager_mod, only : time_type, operator(/=), operator(>), get_time, set_time, &
                             operator(-), operator(/), operator(+)

implicit none
private

public take_obs, interactive_def, get_expected_obs, get_close_states, move_ahead, take_vr, take_td

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



subroutine take_obs(state_vector, location, obs_kind, obs_vals, istatus)
!--------------------------------------------------------------------
!

implicit none

real(r8),            intent(in) :: state_vector(:)
type(location_type), intent(in) :: location
type(obs_kind_type), intent(in) :: obs_kind
real(r8),           intent(out) :: obs_vals
integer,            intent(out) :: istatus

if ( .not. module_initialized ) call initialize_module

! Initially, have only raw state observations implemented so obs_kind is
! irrelevant. Just do interpolation and return.
   call interpolate(state_vector, location, get_obs_kind(obs_kind), obs_vals, istatus)

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


subroutine get_expected_obs(seq, keys, state, obs_vals, istatus)
!---------------------------------------------------------------------------
! Compute forward operator for set of obs in sequence

type(obs_sequence_type), intent(in)  :: seq
integer,                 intent(in)  :: keys(:)
real(r8),                intent(in)  :: state(:)
real(r8),                intent(out) :: obs_vals(:)
integer,                 intent(out) :: istatus

integer             :: num_obs, i
type(location_type) :: location
type(obs_kind_type) :: obs_kind
type(obs_type)      :: obs
type(obs_def_type)  :: obs_def
integer             :: obs_kind_ind

if ( .not. module_initialized ) call initialize_module

num_obs = size(keys)

! NEED to initialize istatus  to okay value
istatus = 0

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
! Otherwise do forward operator
   elseif(obs_kind_ind == KIND_TD) then ! Dew-point temperature
      call take_td(state, obs_def, obs_vals(i), istatus)
   elseif(obs_kind_ind == KIND_VR) then ! Radial velocity
      call take_vr(state, obs_def, obs_vals(i), istatus)
   else
      call take_obs(state, location, obs_kind, obs_vals(i), istatus)
   endif
end do

! Need to destroy this obs to free up allocated space
call destroy_obs(obs)

end subroutine get_expected_obs


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


subroutine take_vr(state_vector, obs_def, vr, istatus)
!--------------------------------------------------------------------
!

implicit none

real(r8),            intent(in) :: state_vector(:)
type(obs_def_type),  intent(in) :: obs_def
real(r8),           intent(out) :: vr
integer,            intent(out) :: istatus

type(location_type) :: location

real(r8) :: u, v, w, p, qr, alpha, wt, direction(3)

if ( .not. module_initialized ) call initialize_module

location = get_obs_def_location(obs_def)
!WRF direction = get_platform_orientation(get_obs_def_platform(obs_def))

call interpolate(state_vector, location, KIND_U, u, istatus)
call interpolate(state_vector, location, KIND_V, v, istatus)
call interpolate(state_vector, location, KIND_W, w, istatus)
call interpolate(state_vector, location, KIND_P, p, istatus)
call interpolate(state_vector, location, KIND_QR, qr, istatus)

alpha=(ps0/p)**0.4_r8
if(qr <= 0.0_r8) then
   wt=0.0_r8
else
   wt=5.4_r8*alpha*qr**0.125_r8
endif

print*,'Terminal velocity = ',wt

! direction(1) = sin(az)cos(elev)
! direction(2) = cos(az)cos(elev)
! direction(3) = sin(elev)
! az and elev are angles at the observation location.

vr = direction(1)*u + direction(2)*v + direction(3)*(w+wt)

end subroutine take_vr


subroutine take_td(state_vector, obs_def, td, istatus)
!--------------------------------------------------------------------
!

implicit none

real(r8),           intent(in)  :: state_vector(:)
type(obs_def_type), intent(in)  :: obs_def
real(r8),           intent(out) :: td
integer,            intent(out) :: istatus

type(location_type) :: location

real(r8) :: t, qv, t_c, es, qvs, p, INVTD, rd_over_rv, rd_over_rv1

if ( .not. module_initialized ) call initialize_module

location = get_obs_def_location(obs_def)

call interpolate(state_vector, location, KIND_P, p, istatus)
call interpolate(state_vector, location, KIND_QV, qv, istatus)
call interpolate(state_vector, location, KIND_T, t, istatus)

if( p /= missing_r8 .and. qv /= missing_r8 .and. t /= missing_r8 ) then

   rd_over_rv = gas_constant / gas_constant_v
   rd_over_rv1 = 1.0_r8 - rd_over_rv

   t_c = t - t_kelvin

!------------------------------------------------------------------------------
!  Calculate saturation vapour pressure:
!------------------------------------------------------------------------------

   es = es_alpha * exp( es_beta * t_c / ( t_c + es_gamma ) )

!------------------------------------------------------------------------------
!  Calculate saturation specific humidity:
!------------------------------------------------------------------------------

   qvs = rd_over_rv * es / ( p - rd_over_rv1 * es )

   INVTD = 1.0_r8/t  - LOG (qv / qvs) / L_over_Rv

   td = 1.0_r8 / INVTD

else

   td = missing_r8

endif

print*,INVTD, t, qv, qvs, es, p

end subroutine take_td


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

type(time_type) :: next_time, time2, start_time, end_time, delta_time
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

! Figure out time to which to advance model 
time2 = aget_closest_state_time_to(ens_time(1), next_time)

! Compute the model time step and center a window around the closest time
delta_time = get_model_time_step()
! WATCH OUT FOR USING BOUNDARY TIME OBS TWICE; add one second to bottom time
start_time = time2 - delta_time / 2 + set_time(1, 0)
end_time = time2 + delta_time / 2

! Get all the obsevations at exactly this time
call get_obs_time_range(seq, start_time, end_time, key_bounds, num_obs_in_set, &
   out_of_range, observation)

call get_time(start_time, secs, days)
write(*, *) 'start time of obs range is (d, s) = ', days, secs
call get_time(end_time, secs, days)
write(*, *) 'end time of obs range   is (d, s) = ', days, secs

! Advance all ensembles (to the time of the first ensemble)
if(time2 /= ens_time(1)) call Aadvance_state(ens_time, ens, ens_size, time2, async, adv_ens_command)

! Release the storage associated with the observation temp varialbe
call destroy_obs(observation)

end subroutine move_ahead

end module obs_model_mod
