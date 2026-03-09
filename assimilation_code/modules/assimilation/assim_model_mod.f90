! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> This module is used to wrap around the basic portions of existing dynamical models to
!> add capabilities needed by the standard assimilation methods.
!>
!>@todo FIXME explore redundant routines , especially 'diagnostic' ones
!>@todo the documentation for this module is out-of-date.

module assim_model_mod

use    types_mod, only : r8, i8

use time_manager_mod, only : time_type,                                            &
                             operator(<), operator(>), operator(+), operator(-),   &
                             operator(/), operator(*), operator(==), operator(/=)

use ensemble_manager_mod, only : ensemble_type

! these are the required model_mod interfaces:
use     model_mod, only : get_model_size, static_init_model, get_state_meta_data,  &
                          get_model_time_step => shortest_time_between_assimilations, &
                          shortest_time_between_assimilations,                     &
                          init_conditions, init_time, adv_1step, end_model,        &
                          pert_model_copies, get_close_obs,                        &
                          model_get_close_state => get_close_state,                &
                          convert_vertical_obs, convert_vertical_state,            &
                          interpolate => model_interpolate,                        &
                          read_model_time, write_model_time

use ensemble_manager_mod, only : ensemble_type

use location_mod, only  : get_close_type, location_type,    &
                          do_strongly_coupled_localization, &
                          get_close_state_strongly_coupled

implicit none
private

public :: static_init_assim_model, &
          end_assim_model, &
          get_initial_condition, &
          get_closest_state_time_to, &
          get_state_meta_data, &
          get_model_time_step, &
          get_model_size,  &
          adv_1step, &
          interpolate, &
          pert_model_copies, &
          get_close_obs, &
          get_close_state, &
          convert_vertical_obs, &
          convert_vertical_state, &
          read_model_time, &
          write_model_time

! Ensure init code is called exactly once
logical :: module_initialized = .false.

!-------------------------------------------------------------

contains

!======================================================================


  subroutine static_init_assim_model()
!----------------------------------------------------------------------
!
! Initializes class data for the assim_model. Also calls the static
! initialization for the underlying model. So far, this simply 
! is initializing the position of the state variables as location types.

! only execute this code once, even if called multiple times.
if (module_initialized) return

! First thing to do is echo info to logfile ... 
module_initialized = .true.

! give the model a chance to initialize itself once
call static_init_model()

end subroutine static_init_assim_model



subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)
!---------------------------------------------------------------------
!
! Allows assim_model to handle get_close computations by modifying the results
! of the model_mod get_close_state or by doing something completely independent.
! Initially used for localization across Earth system components for 
! strongly-coupled localization in threed_sphere models. 
   
type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(in)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:)
integer(i8),                   intent(in)  :: loc_indx(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ens_handle

! Call the model_mod get_close_state to get horizontal distance
call model_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)

! If strongly coupled localization is being applied, compute it
if(do_strongly_coupled_localization()) &
   call get_close_state_strongly_coupled(gc, base_loc, base_type, locs, &
      loc_qtys, loc_indx, num_close, close_ind, dist, ens_handle)

! Vertical localization MUST be off in threed-sphere when doing strongly-coupled localization
! Can force vertical localization to be off if doing two different components
! Will have to specify the different components in the threed-sphere localization
! module. Intercept there.

end subroutine get_close_state


function get_closest_state_time_to(model_time, given_time)
!----------------------------------------------------------------------
!
! Returns the time closest to the given time that the model can reach
! with its state. Initial implementation just assumes fixed timestep.
! Need to describe potentially more general time-stepping capabilities
! from the underlying model in the long run.

type(time_type) :: get_closest_state_time_to
type(time_type), intent(in) :: model_time, given_time

type(time_type) :: time_step

! Get the model time step capabilities - this is the
! shortest amount of time you can/want to ask the model
! to advance the state
time_step = shortest_time_between_assimilations()

if(model_time > given_time) then
   ! If model_time is past start of obs window, don't advance it
   get_closest_state_time_to = model_time
   return
endif

get_closest_state_time_to = model_time

do while((time_step + 2*get_closest_state_time_to) < 2*given_time)
   get_closest_state_time_to = get_closest_state_time_to + time_step
enddo

end function get_closest_state_time_to


subroutine get_initial_condition(initial_time, x)
!----------------------------------------------------------------------
!

type(time_type), intent(out) :: initial_time
real(r8),        intent(out) :: x(:)

call init_conditions(x)

call init_time(initial_time)

end subroutine get_initial_condition


!-------------------------------------------------------------------

subroutine end_assim_model()
!--------------------------------------------------------------------
!
! Closes down assim_model. For now, only thing to do is tell model to end.

call end_model()

end subroutine end_assim_model

!
!===================================================================
! End of assim_model_mod
!===================================================================
!
end module assim_model_mod

