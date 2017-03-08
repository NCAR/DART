! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module obs_model_mod

use       types_mod, only : r8
use   utilities_mod, only : register_module
use    location_mod, only : location_type, interactive_location
use assim_model_mod, only : interpolate
use    obs_kind_mod, only : obs_kind_type, interactive_kind, get_obs_def_type_of_obs

implicit none
private

public :: take_obs, interactive_def

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


logical, save :: module_initialized = .false.

contains

!======================================================================

subroutine initialize_module

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module





function take_obs(state_vector, location, obs_kind)
!--------------------------------------------------------------------
!

implicit none

real(r8) :: take_obs
real(r8), intent(in) :: state_vector(:)
type(location_type), intent(in) :: location
type(obs_kind_type), intent(in) :: obs_kind

if ( .not. module_initialized ) call initialize_module

! Initially, have only raw state observations implemented so obs_kind is
! irrelevant. Just do interpolation and return.

take_obs = interpolate(state_vector, location, get_obs_def_type_of_obs(obs_kind))

end function take_obs



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


end module obs_model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
