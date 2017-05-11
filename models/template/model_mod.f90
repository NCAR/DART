! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This is a template showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. The public interfaces listed
! must all be supported with the argument lists as indicated. Many of the interfaces
! are not required for minimal implementation (see the discussion of each
! interface and look for NULL INTERFACE). 

! Modules that are absolutely required for use are listed
use        types_mod, only : r8, i8, MISSING_R8
use time_manager_mod, only : time_type, set_time, set_time_missing
use     location_mod, only : location_type, get_close_type, &
                             get_close_obs, get_close_state, &
                             convert_vertical_obs, convert_vertical_state, &
                             set_location, set_location_missing
use    utilities_mod, only : register_module, error_handler, &
                             E_ERR, E_MSG
                             ! nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             ! find_namelist_in_file, check_namelist_read
use netcdf_utilities_mod, only : nc_add_global_attribute, nc_sync, &
                                 nc_add_global_creation_time, nc_redef, nc_enddef
use state_structure_mod, only : add_domain
use ensemble_manager_mod, only : ensemble_type
use dart_time_io_mod, only  : read_model_time, write_model_time
use default_model_mod, only : pert_model_copies, nc_write_model_vars

implicit none
private

! required by DART code - will be called from filter and other
! DART executables.  interfaces to these routines are fixed and
! cannot be changed in any way.
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          shortest_time_between_assimilations, &
          end_model,              &
          static_init_model,      &
          nc_write_model_atts,    &
          init_time,              &
          init_conditions

! public but in another module
public :: nc_write_model_vars,    &
          pert_model_copies,      &
          get_close_obs,          &
          get_close_state,        &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time, &
          write_model_time


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! EXAMPLE: define model parameters here
integer, parameter               :: model_size = 3
type(time_type)                  :: time_step
type(location_type), allocatable :: state_loc(:)

! EXAMPLE: perhaps a namelist here 
integer  :: time_step_days      = 0
integer  :: time_step_seconds   = 3600
namelist /model_nml/ time_step_days, time_step_seconds

contains

!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.
! Can be a NULL INTERFACE for the simplest models.

subroutine static_init_model()

 real(r8) :: x_loc
 integer  :: i, dom_id
!integer  :: iunit, io

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! This is where you would read a namelist, for example.
!call find_namelist_in_file("input.nml", "model_nml", iunit)
!read(iunit, nml = model_nml, iostat = io)
!call check_namelist_read(iunit, io, "model_nml")


! Record the namelist values used for the run ...
!if (do_nml_file()) write(nmlfileunit, nml=model_nml)
!if (do_nml_term()) write(     *     , nml=model_nml)

! Create storage for locations
allocate(state_loc(model_size))

! Define the locations of the model state variables
! naturally, this can be done VERY differently for more complicated models.
! set_location() is different for 1D vs. 3D models, not surprisingly.
do i = 1, model_size
   x_loc = (i - 1.0_r8) / model_size
   ! must do one of these:
   !state_loc(i) =  set_location(x_loc)
   !state_loc(i) =  set_location(x_loc,y_loc,v_loc,v_type)
end do

! The time_step in terms of a time type must also be initialized.
time_step = set_time(time_step_seconds, time_step_days)

! tell dart the size of the model
dom_id = add_domain(int(model_size,i8))

end subroutine static_init_model




!------------------------------------------------------------------
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

x = MISSING_R8

end subroutine init_conditions



!------------------------------------------------------------------
! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to 
! compute a timestep, for instance for radiation computations.
! This interface is only called if the namelist parameter
! async is set to 0 in perfect_model_obs of filter or if the 
! program integrate_model is to be used to advance the model
! state as a separate executable. If one of these options
! is not going to be used (the model will only be advanced as
! a separate model-specific executable), this can be a 
! NULL INTERFACE.

subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

end subroutine adv_1step



!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 
! This interface is required for all applications.

function get_model_size()

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size



!------------------------------------------------------------------
! Companion interface to init_conditions. Returns a time that is somehow 
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

subroutine init_time(time)

type(time_type), intent(out) :: time

! for now, just set to 0
time = set_time(0,0)

end subroutine init_time




!------------------------------------------------------------------
! Given a state handle, a location, and a model state variable type,
! interpolates the state variable fields to that location and returns
! the values in expected_obs. The istatus variables should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is a model specific integer that specifies the kind of field (for
! instance temperature, zonal wind component, etc.). In low order
! models that have no notion of types of variables this argument can
! be ignored. For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

subroutine model_interpolate(state_handle, ens_size, location, obs_type, expected_obs, istatus)


type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_type
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

! This should be the result of the interpolation of a
! given kind (itype) of variable at the given location.
expected_obs(:) = MISSING_R8

! The return code for successful return should be 0. 
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.
istatus(:) = 1

end subroutine model_interpolate



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
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument kind
! can be returned if the model has more than one type of field (for
! instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

! these should be set to the actual location and state quantity
location = set_location_missing()
if (present(var_type)) var_type = 0  

end subroutine get_state_meta_data



!------------------------------------------------------------------
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

subroutine end_model()

! deallocate(state_loc)

end subroutine end_model


!------------------------------------------------------------------
! write any additional attributes to the output and diagnostic files

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

! put file into define mode.

call nc_redef(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model_revision", revision )
call nc_add_global_attribute(ncid, "model_revdate", revdate )

call nc_add_global_attribute(ncid, "model", "template")

! Flush the buffer and leave netCDF file open
call nc_sync(ncid)

end subroutine nc_write_model_atts

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
