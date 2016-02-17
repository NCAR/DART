! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> This module is used to wrap around the basic portions of existing dynamical models to
!> add capabilities needed by the standard assimilation methods.
module assim_model_mod

use    types_mod, only : r8, digits12
use location_mod, only : location_type
use time_manager_mod, only : time_type,                                            &
                             operator(<), operator(>), operator(+), operator(-),   &
                             operator(/), operator(*), operator(==), operator(/=)
use utilities_mod, only : register_module, error_handler,    &
                          E_ERR, E_WARN, E_MSG, E_DBG, nmlfileunit,                &
                          dump_unit_attributes, find_namelist_in_file,             &
                          check_namelist_read, do_nml_file, do_nml_term, &
                          set_output,                          &
                          ascii_file_format, set_output
use     model_mod, only : get_model_size, static_init_model, get_state_meta_data,  &
                          get_model_time_step, init_conditions,                    &
                          init_time, adv_1step, end_model,                         &
                          nc_write_model_vars,                                     &
                          get_close_maxdist_init, get_close_obs_init,              &
                          model_interpolate,                                       &
                          get_close_obs, pert_model_copies

use ensemble_manager_mod, only : ensemble_type

implicit none
private

public :: static_init_assim_model, &
          get_model_size,  &
          get_closest_state_time_to, &
          get_initial_condition, &
          get_state_meta_data, &
          get_model_time, copy_assim_model, &
          end_assim_model, &
          assim_model_type, &
          init_assim_model, &
          aget_closest_state_time_to,&
          get_model_time_step, &
          adv_1step, &
          aget_initial_condition, &
          get_close_maxdist_init, &
          get_close_obs_init, interpolate, &
          get_close_obs, &
          pert_model_copies

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


! Type to keep model state and time together
type assim_model_type
   !private
   real(r8), pointer :: state_vector(:)
   type(time_type) :: time
   integer                           :: model_size
   integer                           :: copyID
! would like to include character string to indicate which netCDF variable --
! replace "state" in output_diagnostics ...
end type assim_model_type

! Permanent class storage for model_size
integer :: model_size

! Ensure init code is called exactly once
logical :: module_initialized = .false.

!-------------------------------------------------------------

contains

!======================================================================


subroutine init_assim_model(state)
!----------------------------------------------------------------------
!
! Allocates storage for an instance of an assim_model_type. With this
! implementation, need to be VERY careful about assigment and maintaining
! permanent storage locations. Need to revisit the best way to do 
! assim_model_copy below.

type(assim_model_type), intent(inout) :: state

! Get the model_size from the model
model_size = get_model_size()

allocate(state%state_vector(model_size))
state%model_size = model_size

end subroutine init_assim_model


  subroutine static_init_assim_model()
!----------------------------------------------------------------------
! subroutine static_init_assim_model()
!
! Initializes class data for the assim_model. Also calls the static
! initialization for the underlying model. So far, this simply 
! is initializing the position of the state variables as location types.

implicit none

! only execute this code once, even if called multiple times.
if (module_initialized) return

! First thing to do is echo info to logfile ... 
call register_module(source, revision, revdate)
module_initialized = .true.

! Call the underlying model's static initialization
call static_init_model()

end subroutine static_init_assim_model


function get_closest_state_time_to(assim_model, time)
!----------------------------------------------------------------------
!
! Returns the time closest to the given time that the model can reach
! with its state. Initial implementation just assumes fixed timestep.
! Need to describe potentially more general time-stepping capabilities
! from the underlying model in the long run.

implicit none

type(time_type)                    :: get_closest_state_time_to
type(assim_model_type), intent(in) :: assim_model
type(time_type), intent(in) :: time

type(time_type) :: model_time

model_time = assim_model%time

get_closest_state_time_to = aget_closest_state_time_to(model_time, time)

end function get_closest_state_time_to



function aget_closest_state_time_to(model_time, time)
!----------------------------------------------------------------------
!
! Returns the time closest to the given time that the model can reach
! with its state. Initial implementation just assumes fixed timestep.
! Need to describe potentially more general time-stepping capabilities
! from the underlying model in the long run.

implicit none

type(time_type) :: aget_closest_state_time_to
type(time_type), intent(in) :: model_time, time

type(time_type) :: time_step

! Get the model time step capabilities
time_step = get_model_time_step()

if(model_time > time) then
   ! If model_time is past start of obs window, don't advance it
   aget_closest_state_time_to = model_time
   return
endif

aget_closest_state_time_to = model_time

do while((time_step + 2*aget_closest_state_time_to) < 2*time)
   aget_closest_state_time_to = aget_closest_state_time_to + time_step
enddo

end function aget_closest_state_time_to


subroutine get_initial_condition(x)
!----------------------------------------------------------------------
! function get_initial_condition()
!
! Initial conditions. This returns an initial assim_model_type
! which includes both a state vector and a time. Design of exactly where this 
! stuff should come from is still evolving (12 July, 2002) but for now can 
! start at time offset 0 with the initial state.
! Need to carefully coordinate this with the times for observations.

implicit none

type(assim_model_type), intent(inout) :: x

call aget_initial_condition(x%time, x%state_vector)

end subroutine get_initial_condition



subroutine aget_initial_condition(time, x)
!----------------------------------------------------------------------
! function get_initial_condition()
!
! Initial conditions. This returns an initial state vector and a time
! for use in an assim_model_type.  Design of exactly where this 
! stuff should come from is still evolving (12 July, 2002) but for now can 
! start at time offset 0 with the initial state.
! Need to carefully coordinate this with the times for observations.

implicit none

type(time_type), intent(out) :: time
real(r8),        intent(out) :: x(:)

call init_conditions(x)

call init_time(time)

end subroutine aget_initial_condition



function get_model_time(assim_model)
!-----------------------------------------------------------------------
!
! Returns the time component of a assim_model extended state.

implicit none

type(time_type)                    :: get_model_time
type(assim_model_type), intent(in) :: assim_model

get_model_time = assim_model%time

end function get_model_time


subroutine copy_assim_model(model_out, model_in)
!-------------------------------------------------------------------------
!
! Does a copy of assim_model, should be overloaded to =? Still need to be
! very careful about trying to limit copies of the potentially huge state
! vectors for big models.  Interaction with pointer storage?

implicit none

type(assim_model_type), intent(inout) :: model_out
type(assim_model_type), intent(in)    :: model_in

integer :: i

! Need to make sure to copy the actual storage and not just the pointer (verify)
model_out%time       = model_in%time
model_out%model_size = model_in%model_size

do i = 1, model_in%model_size
   model_out%state_vector(i) = model_in%state_vector(i)
end do

end subroutine copy_assim_model

!> Pass through routine to model interpolate
subroutine interpolate(state_handle, ens_size, location, loctype, expected_obs, istatus)
!---------------------------------------------------------------------
!
! Interpolates from the state vector in an assim_model_type to the
! location. Will need to be generalized for more complex state vector
! types. It might be better to be passing an assim_model_type with
! the associated time through here, but that requires changing the
! entire observation side of the class tree. Reconsider this at a 
! later date (JLA, 15 July, 2002). loctype for now is an integer that
! specifies what sort of variable from the model should be interpolated.

implicit none

type(ensemble_type),   intent(in)    :: state_handle
integer,               intent(in)    :: ens_size
type(location_type),   intent(in)    :: location
integer,               intent(in)    :: loctype
real(r8),              intent(out)   :: expected_obs(ens_size)
integer,               intent(out)   :: istatus(ens_size)

istatus = 0

call model_interpolate(state_handle, ens_size, location, loctype, expected_obs, istatus)

end subroutine interpolate

!-------------------------------------------------------------------

subroutine end_assim_model()
!--------------------------------------------------------------------
!
! Closes down assim_model. For now, only thing to do is tell model to end.

implicit none

call end_model()

end subroutine end_assim_model

!
!===================================================================
! End of assim_model_mod
!===================================================================
!
end module assim_model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
