! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module default_model_mod

!> bypass routines for all required entry points.
!> if a model has no need of a routine, use these instead.
!>
!> see below for brief notes about the purpose of each routine.

use        types_mod,      only : r8, i8, i4, MISSING_R8

use time_manager_mod,      only : time_type, set_time

use     location_mod,      only : location_type, set_location, set_location_missing, &
                                  get_close_type, get_close_obs, get_close_state, &
                                  convert_vertical_obs, convert_vertical_state

use    utilities_mod,      only : error_handler, E_ERR, E_MSG, nmlfileunit, &
                                  do_output, find_namelist_in_file, check_namelist_read,     &
                                  do_nml_file, do_nml_term

use netcdf_utilities_mod,  only : nc_check

use ensemble_manager_mod,  only : ensemble_type

use dart_time_io_mod,      only : read_model_time, write_model_time

implicit none
private

public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          shortest_time_between_assimilations, &
          end_model,              &
          static_init_model,      &
          init_time,              &  ! use this OR fail_
          fail_init_time,         &
          init_conditions,        &  ! use this or fail_
          fail_init_conditions,   &
          nc_write_model_atts,    &
          nc_write_model_vars,    &  ! currently unused
          pert_model_copies,      &
          get_close_obs,          &  ! from the location module
          get_close_state,        &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time,        &  ! from the dart_time_io module
          write_model_time

character(len=*), parameter :: source = 'utilities/default_model_mod.f90'

contains

!==================================================================

!> normal contents:  read the namelist, initialize module variables,
!> call add_domain() to set what data should be read into the state,
!> initialize grid information.

subroutine static_init_model()

end subroutine static_init_model

!------------------------------------------------------------------

!> set hardcoded values in the state vector for a 'cold start' if not 
!> reading from a file.  normally unused for large models because
!> it's not possible to set values in a full state vector from thin air.
!> note that a model_mod should either supply an init_conditions routine,
!> use this one, or set:  init_conditions => fail_init_conditions
!> in the module use list if it isn't supported.

subroutine init_conditions(x)
real(r8), intent(out) :: x(:)

! default
x = 0.0_r8

end subroutine init_conditions

!------------------------------------------------------------------

subroutine fail_init_conditions(x)
real(r8), intent(out) :: x(:)

call error_handler(E_ERR, 'init_conditions', &
           'this model cannot provide initial conditions', source)

! default
x = 0.0_r8

end subroutine fail_init_conditions

!------------------------------------------------------------------

!> for models which can be called as subroutines, this routine should
!> take state vector x() and advance it to the requested time.
!>
!> for any model which is not subroutine-callable or needs to be advanced
!> outside of filter, the model_mod should just use this routine.

subroutine adv_1step(x, time)
real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

call error_handler(E_ERR, 'adv_1step', 'unable to advance model', source)

end  subroutine adv_1step

!------------------------------------------------------------------

!> must compute the size of the state vector, which is the total number 
!> of items in all variables.  this is usually done by reading a file
!> that gives the grid size and then multiplying by the number of variables
!> on that grid.

function get_model_size()
integer :: get_model_size

get_model_size = 1

end function get_model_size

!------------------------------------------------------------------

!> set initial time of the run if not reading the time from a restart file.
!> normally unused for large models because the time comes with the
!> data in a restart file.
!>
!> note that a model_mod should either supply an init_time routine,
!> use this one, or set:  init_time => fail_init_time
!> in the module use list if it isn't supported.

subroutine init_time(time)
type(time_type), intent(out) :: time

time = set_time(0, 0)

end subroutine init_time

!------------------------------------------------------------------

subroutine fail_init_time(time)
type(time_type), intent(out) :: time

call error_handler(E_ERR, 'init_time', &
           'this model cannot provide an initial time', source)

time = set_time(0, 0)

end subroutine fail_init_time

!------------------------------------------------------------------

!> the main interpolation routine.  this routine needs to return an
!> ensemble-sized array of 'expected_obs' values, given a location
!> and a quantity.  previous versions of dart passed in a full state
!> vector (the array x(:)) and the interpolate routine returned only
!> the expected value for that single member.  this version gets
!> an ensemble handle (state_handle), along with a location and quantity.
!> it is expected to return an expected value for all ensemble members.
!> 
!> if converting code from a previous version of dart, look for code
!> where a = x(n) is referenced and change it to a(:) = get_state(n, state_handle)

subroutine model_interpolate(state_handle, ens_size, location, obs_quantity, expected_obs, istatus)

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type),  intent(in) :: location
integer,              intent(in) :: obs_quantity
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

! all fail
expected_obs(:) = MISSING_R8
istatus(:) = 1

end subroutine model_interpolate

!------------------------------------------------------------------

!> set both the current assimilation window size, and the minimum time
!> filter could ask the model to advance, if it's advancing inside a 
!> filter run.
!>
!> this replaces get_model_time_step() in previous versions of dart.
!> the function is the same but the name was confusing.  models often
!> have an internal time step that's very short and is unrelated to
!> what this routine returns.
!>
!> most large models have become too complex for filter to manage,
!> so the model advances are controlled by an external script.
!> in that case, this time only controls the assimilation window size.

function shortest_time_between_assimilations()
type(time_type) :: shortest_time_between_assimilations

! default to 1 day
shortest_time_between_assimilations = set_time(0, 1)

end function shortest_time_between_assimilations

!------------------------------------------------------------------

!> for any given offset into the state vector, return the location
!> and optionally the quantity of that item.  for models using a
!> grid, the state_structure routines can assist in converting a 1D
!> offset into an (i,j,k) index of a grid.

subroutine get_state_meta_data(state_handle, index_in, location, var_type)

type(ensemble_type), intent(in)  :: state_handle !< some large models need this
integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

location = set_location_missing()
if (present(var_type)) var_type = 0    ! default variable type

end subroutine get_state_meta_data

!------------------------------------------------------------------

!> release any allocated storage, close any open file units,
!> clean up in any way needed.  if no cleanup is needed, the
!> model_mod can just use this routine.

subroutine end_model()

end subroutine end_model

!------------------------------------------------------------------

!> if filter is creating an output file (instead of opening an existing
!> file), this routine is called to allow the model_mod to add any additional
!> information such as the grid, or non-state variables.  the state values
!> will be written out by other dart routines so this code SHOULD NOT
!> create or write state variables values.

subroutine nc_write_model_atts(ncid, domain_id) 

integer, intent(in) :: ncid
integer, intent(in) :: domain_id

end subroutine nc_write_model_atts

!------------------------------------------------------------------

!> this routine is currently unused, so model_mods should just use
!> this routine to satisfy the compiler.

subroutine nc_write_model_vars(ncid, domain_id, state_ens_handle, memberindex, timeindex)

integer,             intent(in) :: ncid      
integer,             intent(in) :: domain_id
type(ensemble_type), intent(in) :: state_ens_handle
integer, optional,   intent(in) :: memberindex
integer, optional,   intent(in) :: timeindex

end subroutine nc_write_model_vars

!--------------------------------------------------------------------

!> create an ensemble of states from a single state.  if this routine is
!> not specialized by the model_mod, the default is to add gaussian noise
!> to all parts of the state vector.  unless the model_mod wants to add
!> different amounts of noise to different parts of the state, the model_mod
!> can use this routine.

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,   intent(in) :: ens_size
real(r8),  intent(in) :: pert_amp
logical,  intent(out) :: interf_provided

interf_provided = .false.

end subroutine pert_model_copies

!===================================================================
! End of model_mod
!===================================================================
end module default_model_mod

