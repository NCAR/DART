! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

! This is the interface between pywatershed and DART

use types_mod,             only : r8, i8, i4, MISSING_R8, obstypelength, &
                                  earth_radius

use time_manager_mod,      only : time_type, set_time

use location_mod,          only : location_type, set_location, get_location,  &
                                  get_close_obs, get_close_state,             &
                                  convert_vertical_obs, convert_vertical_state

use utilities_mod,         only : do_nml_file, do_nml_term,    &
                                  nmlfileunit, find_namelist_in_file,           &
                                  check_namelist_read

use location_io_mod,      only :  nc_write_location_atts, nc_write_location

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, nc_begin_define_mode, &
                                 nc_end_define_mode

use         obs_kind_mod,  only : QTY_STATE_VARIABLE

use ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain

use default_model_mod,     only : adv_1step, end_model, pert_model_copies, & 
                                  nc_write_model_vars, MAX_STATE_VARIABLE_FIELDS_CLAMP

use dart_time_io_mod,      only : read_model_time, write_model_time

implicit none
private

! required by DART code - will be called from filter and other
! DART executables.  interfaces to these routines are fixed and
! cannot be changed in any way.
public :: get_model_size,       &
          get_state_meta_data,  &
          model_interpolate,    &
          shortest_time_between_assimilations, &
          static_init_model,    &
          init_conditions,      &
          init_time,            &
          nc_write_model_atts

! required by DART but passed through from another module. 
! To write model specific versions of these routines
! remove the routine from  use statement above and add your code to
! this the file.
public :: pert_model_copies,      &
          nc_write_model_vars,    &
          get_close_obs,          &
          get_close_state,        &
          end_model,              &
          adv_1step,              &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time,        &
          write_model_time

character(len=256), parameter :: source   = "pywatershed/model_mod.f90"

type(location_type), allocatable :: state_loc(:)  ! state locations, compute once and store for speed

type(time_type) :: time_step
integer(i8)     :: model_size

! Model namelist declarations with defaults
integer            :: assimilation_period_days     = 1           ! model instantaneous streamflow at the hour 23 of the day
integer            :: assimilation_period_seconds  = 0           ! NA
integer            :: debug                        = 0           ! debugging flag
real(r8)           :: model_perturbation_amplitude = 0.1         ! perturb parameter for initial ensemble generation
real(r8)           :: max_link_distance            = 10000.0     ! Max distance along the stream: 10 km
character(len=256) :: perturb_distribution         = 'lognormal' ! distribution needed for initial ensemble generation
character(len=256) :: domain_order(3)              = ''          ! Domains: channel, hru, params??
character(len=256) :: domain_shapefiles(3)         = ''          ! template restart files for each domain

character(len=obstypelength) :: channel_variables(MAX_STATE_VARIABLE_FIELDS_CLAMP) = '' ! channel state variables
character(len=obstypelength) :: hru_variables(MAX_STATE_VARIABLE_FIELDS_CLAMP)     = '' ! hru (hydrologic response unit) variables
character(len=obstypelength) :: parameters(MAX_STATE_VARIABLE_FIELDS_CLAMP)        = '' ! parameters

namelist /model_nml/ assimilation_period_days,     &
                     assimilation_period_seconds,  &
                     model_perturbation_amplitude, &
                     perturb_distribution,         &
                     max_link_distance,            &
                     domain_order,                 &
                     domain_shapefiles,            &
                     debug,                        &
                     channel_variables,            &
                     hru_variables,                & 
                     parameters

contains


!---------------------------------------------------------------------
! Initialize the model, read namelist, configure domains and variables

subroutine static_init_model()

integer  :: iunit, io
real(r8) :: x_loc
integer  :: i, dom_id

! Read the namelist 
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Output the namelist values if requested
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

model_size = 0

end subroutine static_init_model


!-------------------------------------------------------------
! The following 2 routines are only used if namelist parameter
! start_from_restart is set to .false.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

x = MISSING_R8

end subroutine init_conditions

subroutine init_time(time)

type(time_type), intent(out) :: time

! for now, just set to 0
time = set_time(0,0)

end subroutine init_time


!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 
! This interface is required for all applications.

function get_model_size()

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size



!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.
! This interface is required for all applications.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations



!------------------------------------------------------------------
! Given a state handle, a location, and a model state variable quantity,
! interpolates the state variable fields to that location and returns
! the values in expected_obs. The istatus variables should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is an integer that specifies the quantity of field (for
! instance temperature, zonal wind component, etc.). In low order
! models that have no notion of types of variables this argument can
! be ignored. For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

subroutine model_interpolate(state_handle, ens_size, location, iqty, expected_obs, istatus)

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type),  intent(in) :: location
integer,              intent(in) :: iqty
real(r8),            intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,             intent(out) :: istatus(ens_size)

! This should be the result of the interpolation of a
! given quantity (iqty) of variable at the given location.
expected_obs(:) = MISSING_R8

! The return code for successful return should be 0. 
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.
istatus(:) = 1

end subroutine model_interpolate



!------------------------------------------------------------------
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument quantity
! (qty) can be returned if the model has more than one type of field
! (for instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

subroutine get_state_meta_data(index_in, location, qty_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: qty_type

! these should be set to the actual location and state quantity
location = state_loc(index_in)
if (present(qty_type)) qty_type = QTY_STATE_VARIABLE 

end subroutine get_state_meta_data


!------------------------------------------------------------------
! Writes model-specific attributes to a netCDF file

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid
integer, intent(in) :: domain_id

! put file into define mode.

integer :: msize

msize = int(model_size, i4)

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )

call nc_add_global_attribute(ncid, "model", "template")

call nc_write_location_atts(ncid, msize)
call nc_end_define_mode(ncid)
call nc_write_location(ncid, state_loc, msize)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

