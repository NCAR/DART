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
use assim_model_mod,  only : interpolate, &
                             am_get_close_states => get_close_states
use    obs_kind_mod,  only : obs_kind_type, interactive_kind, get_obs_kind
use obs_sequence_mod, only : obs_sequence_type, obs_type, get_obs_from_key, &
                             get_obs_def, init_obs, destroy_obs
use obs_def_mod,      only : obs_def_type, get_obs_def_location, get_obs_def_kind
use utilities_mod,    only : error_handler, E_ERR


implicit none
private

public take_obs, interactive_def, get_expected_obs, get_close_states

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

stop 'here'
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

end module obs_model_mod
