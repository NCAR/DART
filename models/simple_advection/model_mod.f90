! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

use        types_mod, only : r8, PI, i4, i8

use time_manager_mod, only : time_type, set_time, get_time

use     location_mod, only : location_type, set_location, get_location,  &
                             LocationDims, LocationName, LocationLName,  &
                             get_close_maxdist_init, get_close_obs_init, &
                             loc_get_close_obs => get_close_obs, get_close_type

use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             nmlfileunit, find_namelist_in_file,           &
                             check_namelist_read, nc_check, do_output,     &
                             do_nml_file, do_nml_term

use     obs_kind_mod, only : QTY_VELOCITY, QTY_TRACER_CONCENTRATION, &
                             QTY_TRACER_SOURCE, QTY_MEAN_SOURCE, QTY_SOURCE_PHASE

use random_seq_mod,   only : random_seq_type, init_random_seq, random_gaussian

use ensemble_manager_mod,  only : ensemble_type, get_allow_transpose, &
                                  all_vars_to_all_copies, all_copies_to_all_vars, &
                                  init_ensemble_manager, end_ensemble_manager

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain, state_structure_info, &
                                  add_dimension_to_variable, finished_adding_domain

use dart_time_io_mod,      only : read_model_time, write_model_time

implicit none
private

public :: get_model_size, &
          adv_1step, &
          get_state_meta_data, &
          model_interpolate, &
          get_model_time_step, &
          end_model, &
          static_init_model, &
          init_time, &
          init_conditions, &
          nc_write_model_atts, &
          nc_write_model_vars, &
          pert_model_copies, &
          get_close_maxdist_init, &
          get_close_obs_init, &
          get_close_obs, &
          vert_convert, &
          query_vert_localization_coord, &
          read_model_time, &
          write_model_time

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: string1

! Simplest 1D advection model with spatially-constant wind

! Module storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq
logical :: random_seq_init = .false.

! Permanent static initialized information about domain width
real(r8) :: domain_width_meters 

! Base tracer source rate in amount (unitless) / second
real(r8), allocatable :: mean_source(:)
real(r8), allocatable :: source_phase_offset(:)
real(r8), allocatable :: source_random_amp(:)

!---------------------------------------------------------------

! Namelist with default values
!
integer(i8)  :: num_grid_points         = 10
!!!integer  :: grid_spacing_meters     = 100 * 1000
integer  :: grid_spacing_meters     = 100000
integer  :: time_step_days          = 0
integer  :: time_step_seconds       = 3600

! if you have a template file, set the name here.
! if you are not starting from a restart, leave this as ''
! and the code should construct a domain from spec.
character(len=256) :: template_file = ''

! Namelist parameters associated with wind
! Base velocity (expected value over time), in meters/second
real(r8)  :: mean_wind              = 20.0_r8
! Random walk amplitude for wind in meters / second**2
!!!real(r8)  :: wind_random_amp        = 1.0_r8 / 3600.0_r8
real(r8)  :: wind_random_amp        = 0.0002777778
! Rate of damping towards mean wind value in fraction/second
!!!real(r8)  :: wind_damping_rate      = 0.01_r8 / 3600.0_r8
real(r8)  :: wind_damping_rate      = 0.000002777778
! Can use Eulerian (unstable) or lagrangian (stable)
logical   :: lagrangian_for_wind    = .true.

! Namelist parameters associated with tracer
! Tracer destruction rate in fraction per second
!!!real(r8)  :: destruction_rate       = 0.2_r8 / 3600.0_r8
real(r8)  :: destruction_rate       = 0.00005555556

! Namelist parameters associated with tracer source model
! Random walk amplitude for source as fraction of mean source in .../second**2
real(r8)  :: source_random_amp_frac = 0.00001_r8
! Damping toward mean source rate in fraction/second
!!!real(r8)  :: source_damping_rate    = 0.01_r8 / 3600.0_r8
real(r8)  :: source_damping_rate    = 0.000002777778
! Relative amplitude of diurnal cycle of source (dimensionless)
real(r8)  :: source_diurnal_rel_amp = 0.05_r8
real(r8)  :: source_phase_noise     = 0.0_r8

integer, parameter :: NVARS = 5
integer :: my_ens_size = 1

namelist /model_nml/ num_grid_points, grid_spacing_meters, &
                     time_step_days, time_step_seconds, &
                     mean_wind, wind_random_amp, wind_damping_rate, &
                     lagrangian_for_wind, destruction_rate, &
                     source_random_amp_frac, source_damping_rate, &
                     source_diurnal_rel_amp, source_phase_noise, &
                     template_file

!----------------------------------------------------------------

! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
type(time_type)                  :: time_step

!==================================================================

contains

subroutine static_init_model()
!------------------------------------------------------------------
! Initializes class data for this model. For now, simply outputs the
! identity info, sets the location of the state variables, and initializes
! the time type for the time stepping

real(r8) :: x_loc
integer  :: i, iunit, io, j, dom_id, var_id

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Create storage for locations
allocate(state_loc(NVARS*num_grid_points))

! Define the locations of the model state variables
do i = 1, num_grid_points
   x_loc = (i - 1.0_r8) / num_grid_points
   do j = 0, 4
      state_loc(num_grid_points * j + i) = set_location(x_loc)
   enddo
enddo

! The time_step in terms of a time type must also be initialized.
time_step = set_time(time_step_seconds, time_step_days)

! Allocate storage for the base tracer source and initialize it
! Might want to put the source distribution into the namelist somehow?
allocate(mean_source(num_grid_points), &
 source_phase_offset(num_grid_points), &
   source_random_amp(num_grid_points))

! Set the base value for the mean source; This case has a single enhanced source
mean_source    = 0.1_r8
mean_source(1) = 1.0_r8

! A rapidly varying source
!!!mean_source = 0.1_r8
!!!do i = 1, 5, 2
   !!!mean_source(i) = 1.0_r8
!!!enddo

! Set the base value for the diurnal phase offset; initially on diurnal cycle
source_phase_offset = 0.0_r8

! Set the amplitude of the random walk for the source in unit/second^2
source_random_amp = source_random_amp_frac * mean_source

! Compute the width of the domain in meters
domain_width_meters = num_grid_points * grid_spacing_meters

! For any routines that want to use the random number
! generator later on, initialize it.
if(.not. random_seq_init) then
   call init_random_seq(random_seq, 1)
   random_seq_init = .true.
endif

! Tell the DART I/O routines how large the model data is so they
! can read/write it.
if (template_file /= '') then
   dom_id = add_domain(template_file, NVARS, &
                       (/ 'concentration', &
                          'mean_source  ', &
                          'source       ', &
                          'source_phase ', &
                          'wind         ' /))
else 
   !>@todo FIXME : should not need a template file if initializing members from code

   write(string1, *) 'template file is required for now'
   call error_handler(E_ERR,'static_init_model',string1, source, revision, revdate)

   dom_id = add_domain(NVARS, (/ 'concentration', &
                                 'mean_source  ', &
                                 'source       ', &
                                 'source_phase ', &
                                 'wind         ' /))

   do var_id=1, NVARS
      call add_dimension_to_variable(dom_id, var_id, 'time', 1)
      call add_dimension_to_variable(dom_id, var_id, 'member', my_ens_size)
      call add_dimension_to_variable(dom_id, var_id, 'location', int(num_grid_points, i4))
   enddo
   
   call finished_adding_domain(dom_id)
endif


end subroutine static_init_model



subroutine init_conditions(x)
!------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Initial conditions for simplest 1d Advection.

real(r8), intent(out) :: x(:)

! Start by zeroing all
x = 0.0_r8

! First set of variables is the concentration; Everybody starts at 0.0
x(1:num_grid_points) = 0.0_r8

! Set initial source to mean_source
x(num_grid_points + 1 : 2*num_grid_points) = mean_source

! Third set of variables is the u wind; its units are meters/second
x(2*num_grid_points + 1: 3*num_grid_points) = mean_wind

! Fourth set of variables is the time mean source
x(3*num_grid_points + 1 : 4*num_grid_points) = mean_source

! Fifth set of variables is the source phase offset 
x(4*num_grid_points + 1 : 5*num_grid_points) = source_phase_offset

end subroutine init_conditions



subroutine adv_1step(x, time)
!------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
! Does single time step advance for simplest 1D advection model.
! The Time argument is needed for compatibility with more complex models
! that need to know the time to compute their time tendency. 

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

integer :: next, prev, seconds, days, ens_size
type(location_type) :: source_loc
real(r8) :: lctn, source_location, old_u_mean, new_u_mean
real(r8) :: du_dx, dt_seconds, t_phase, phase, random_src
type(ensemble_type) :: temp_handle

integer(i8) :: i
integer  :: istatus(1)
real(r8) :: new_x(size(x),1)

!>@todo  model_interpolate must have a real ensemble_handle
ens_size = 1
call init_ensemble_manager(temp_handle, ens_size, int(NVARS*num_grid_points,i8))
temp_handle%copies(1,:) = x(:)

! State is concentrations (num_grid_points), source(num_grid_points), u(num_grid_points),
! mean_source(num_grid_points), and source_phase_offset(num_grid_points)
! all dimensioned num_grid_points.


! For the concentration do a linear interpolated upstream distance
! Also have an option to do upstream for wind (controlled by namelist; see below).
! Velocity is in meters/second. Need time_step in seconds
dt_seconds = time_step_seconds + time_step_days * 24.0_r8 * 3600.0_r8

do i = 1, num_grid_points

   ! Find the location of the target point
   lctn = get_location(state_loc(i))

   ! Find source point location: Velocity is in meters/second
   ! Figure out meters to move and then convert to fraction of domain
    
   source_location = lctn - x(2*num_grid_points + i) * dt_seconds / &
      domain_width_meters

   if(source_location > 1.0_r8) &
      source_location = source_location - int(source_location)

   if(source_location < 0.0_r8) &
      source_location = source_location - int(source_location) + 1.0_r8

   source_loc = set_location(source_location)
   call model_interpolate(temp_handle, ens_size, source_loc, QTY_TRACER_CONCENTRATION, new_x(i,:), istatus, x)  

   ! Following line does lagangian du

   if(lagrangian_for_wind) &
      call model_interpolate(temp_handle, ens_size, source_loc, QTY_VELOCITY, new_x(2*num_grid_points + i,:), istatus, x)  
enddo


!----- Next block is eulerian du for wind, controlled by namelist
if(.not. lagrangian_for_wind) then
   ! Use centered for du/dx
   do i = 1, num_grid_points
      ! Compute centered difference for du/dx
      next = i + 1
      if(next > num_grid_points) next = 1
      prev = i - 1
      if(prev < 1) prev = num_grid_points
      du_dx = (x(2*num_grid_points + next) - x(2*num_grid_points + prev)) / &
                                       (2.0_r8 * grid_spacing_meters)
      new_x(2*num_grid_points + i,1) = x(2*num_grid_points + i) + &
                                       x(2*num_grid_points + i) * du_dx * dt_seconds
   enddo
endif
!---- End Eulerian block


! Now add in the source contribution and put concentration and wind back into inout x
! Source is in units of .../second
do i = 1, num_grid_points
   x(i) = new_x(i,1) + x(num_grid_points + i) * dt_seconds
   ! Also copy over the new velocity
   x(2*num_grid_points + i) = new_x(2*num_grid_points + i,1)
enddo

! Now do the destruction rate: Units are fraction destroyed per second
do i = 1, num_grid_points
   if(destruction_rate * dt_seconds > 1.0_r8) then
      x(i) = 0.0_r8
   else
      x(i) = x(i) * (1.0_r8 - destruction_rate*dt_seconds)
   endif
enddo

!----- Following block is random walk plus damping to mean for wind
! Random walk for the spatial mean velocity
old_u_mean = sum(x(2*num_grid_points + 1 : 3*num_grid_points)) / num_grid_points

! Add in a random walk to the mean
new_u_mean = random_gaussian(random_seq, old_u_mean, wind_random_amp * dt_seconds)

! Damp the mean back towards base value
! Need an error handler call here or earlier
if(wind_damping_rate*dt_seconds > 1.0_r8) stop
new_u_mean = new_u_mean - wind_damping_rate * dt_seconds * (new_u_mean - mean_wind)

! Substitute the new mean wind
x(2*num_grid_points  + 1 : 3*num_grid_points) = &
   x(2*num_grid_points + 1 : 3*num_grid_points) - old_u_mean + new_u_mean

! Add some noise into each wind element
do i = 1, num_grid_points
   x(2*num_grid_points + i) = random_gaussian(random_seq, x(2*num_grid_points + i), &
      wind_random_amp * dt_seconds)
enddo

!----- End forced damped wind section


!----- Time tendency with diurnal cycle for sources ---
! Figure out what the time in the model is
call get_time(time, seconds, days)
t_phase = (2.0_r8 * PI * seconds / 86400.0_r8) 
do i = 1, num_grid_points
   phase = t_phase + x(4*num_grid_points + i)
   x(num_grid_points + i) = x(num_grid_points + i) &
      + x(3*num_grid_points + i) * source_diurnal_rel_amp * cos(phase)
      !!!+ x(3*num_grid_points + i) * source_diurnal_rel_amp * cos(phase) * dt_seconds
   ! Also add in some random walk
   random_src = random_gaussian(random_seq, 0.0_r8, source_random_amp(i))
   x(num_grid_points + i) = x(num_grid_points + i) + random_src * dt_seconds
   if(x(num_grid_points + i) < 0.0_r8) x(num_grid_points + i) = 0.0_r8
   ! Finally, damp back towards the base value
   ! Need an error handler call here
   if(source_damping_rate*dt_seconds > 1.0_r8) stop
   x(num_grid_points + i) = x(num_grid_points + i) - &
      source_damping_rate*dt_seconds * (x(num_grid_points + i) - x(3*num_grid_points + i))
enddo

!----- End sources time tendency -----

! Time tendency for source model variables mean_source and source_phase_offset 
! is 0 for now. Might want to add in some process noise for filter advance.

! Process noise test for source_phase_offset
do i = 1, num_grid_points
   x(4*num_grid_points + i) = random_gaussian(random_seq, x(4*num_grid_points + i), &
      source_phase_noise*dt_seconds)
enddo

call end_ensemble_manager(temp_handle)

end subroutine adv_1step




function get_model_size()
!------------------------------------------------------------------
! function get_model_size()
!
! Returns size of model

integer(i8) :: get_model_size

get_model_size = NVARS*num_grid_points

end function get_model_size



subroutine init_time(time)
!------------------------------------------------------------------
!
! Gets the initial time for a state from the model. Where should this info
! come from in the most general case?

type(time_type), intent(out) :: time

! For now, just set to 0
time = set_time(0, 0)

end subroutine init_time



!------------------------------------------------------------------
!> this is a non-standard interface to model_interpolate.  it allows
!> an optional last argument, which is present should be a model state
!> and is used instead of the state handle.  this model_mod uses interpolate
!> during the adv_1step() when there's no state handle available; a naked
!> state is passed in to be advanced.    (having optional args is ok with
!> filter - it just never passes in that final arg.)

subroutine model_interpolate(state_handle, ens_size, location, itype, expected_val, istatus, x)
!
! Interpolates from state vector x to the location. 
! This code supports three obs types: concentration, source, and u

type(ensemble_type),  intent(in)  :: state_handle
integer,              intent(in)  :: ens_size
type(location_type),  intent(in)  :: location
integer,              intent(in)  :: itype
real(r8),             intent(out) :: expected_val(ens_size)
integer,              intent(out) :: istatus(ens_size)
real(r8), optional,   intent(in)  :: x(:)  ! old format state vector, not distributed

integer(i8) :: lower_index, upper_index

integer :: offset
real(r8) :: lctn, lctnfrac

! All forward operators supported
istatus(:) = 0

! Convert location to real
lctn = get_location(location)
! Multiply by model size assuming domain is [0, 1] cyclic
lctn = num_grid_points * lctn

lower_index = int(lctn) + 1
upper_index = lower_index + 1
if(lower_index > num_grid_points) lower_index = lower_index - num_grid_points
if(upper_index > num_grid_points) upper_index = upper_index - num_grid_points

lctnfrac = lctn - int(lctn)

! Now figure out which type
if(itype == QTY_TRACER_CONCENTRATION) then
   offset = 0
else if(itype == QTY_TRACER_SOURCE) then
   offset = num_grid_points
else if(itype == QTY_VELOCITY) then
   offset = 2*num_grid_points
else
   write(string1, *) 'itype is not supported in model_interpolate', itype
   call error_handler(E_ERR,'model_interpolate',string1, source, revision, revdate)
endif

! Add the offsets into the lower and upper indices
lower_index = lower_index + offset
upper_index = upper_index + offset

! Do the interpolation and return.  first case is the normal RMA version
! where expected_val is an ens_size array.  the second is a single state.
if (.not. present(x)) then
  expected_val = (1.0_r8 - lctnfrac) *  get_state(lower_index, state_handle) + &
                           lctnfrac  *  get_state(upper_index, state_handle)
else
  expected_val = (1.0_r8 - lctnfrac) *  x(lower_index) + &
                           lctnfrac  *  x(upper_index)
endif

end subroutine model_interpolate



function get_model_time_step()
!------------------------------------------------------------------
! function get_model_time_step()
!
! Returns the the time step of the model. In the long run should be replaced
! by a more general routine that returns details of a general time-stepping
! capability.

type(time_type) :: get_model_time_step

get_model_time_step = time_step

end function get_model_time_step



subroutine get_state_meta_data(state_handle, index_in, location, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?

type(ensemble_type), intent(in)  :: state_handle !< some large models need this
integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

integer :: var_type_index, var_loc_index

! Three variable types
var_type_index = (index_in - 1) / num_grid_points + 1
var_loc_index = index_in - (var_type_index - 1)*num_grid_points

if(present(var_type)) then
   if(var_type_index == 1) then
      var_type = QTY_TRACER_CONCENTRATION
   else if(var_type_index == 2) then
      var_type = QTY_TRACER_SOURCE
   else if(var_type_index == 3) then
      var_type = QTY_VELOCITY
   else if(var_type_index == 4) then
      var_type = QTY_MEAN_SOURCE
   else if(var_type_index == 5) then
      var_type = QTY_SOURCE_PHASE
   endif
endif

location = state_loc(var_loc_index)

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Nothing for now.


end subroutine end_model



subroutine model_get_close_states(o_loc, radius, inum, indices, dist, x)
!------------------------------------------------------------------
! 
! Stub for computation of get close states

type(location_type), intent(in)  :: o_loc
real(r8),            intent(in)  :: radius
integer,             intent(out) :: inum
integer,             intent(in)  :: indices(:)
real(r8),            intent(in)  :: dist(:)
real(r8),            intent(in)  :: x(:)

! Because of F90 limits this stub must be here telling assim_model
! to do exhaustive search (inum = -1 return)
inum = -1

end subroutine model_get_close_states



function nc_write_model_atts( ncFileID, model_mod_writes_state_variables ) result (ierr)
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file

use typeSizes
use netcdf

integer, intent(in)  :: ncFileID      ! netCDF file identifier
logical, intent(out) :: model_mod_writes_state_variables
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: LocationDimID, LocationVarID
integer :: i
type(location_type) :: lctn 

ierr = 0                             ! assume normal termination

! dart code will write the state into the file
! so this routine just needs to write any model-specific
! attributes it wants to record.

model_mod_writes_state_variables = .false.

! are the inq and sync needed?  can we just redef here?
call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID), &
                           'nc_write_model_atts', 'nf90_Inquire')
call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'nf90_sync') ! Ensure netCDF file is current
call nc_check(nf90_Redef(ncFileID), 'nc_write_model_atts', 'nf90_Redef')

! Write Global Attributes 

call put_global_creation_time(ncFileID)

call put_global_char_att(ncFileID, "model_source"           , source )
call put_global_char_att(ncFileID, "model_revision"         , revision )
call put_global_char_att(ncFileID, "model_revdate"          , revdate )
call put_global_char_att(ncFileID, "model"                  , "simple_advection" )
call put_global_real_att(ncFileID, "mean_wind"              , mean_wind )
call put_global_real_att(ncFileID, "wind_random_amp"        , wind_random_amp )
call put_global_real_att(ncFileID, "wind_damping_rate"      , wind_damping_rate )
call put_global_real_att(ncFileID, "destruction_rate"       , destruction_rate )
call put_global_real_att(ncFileID, "source_random_amp_frac" , source_random_amp_frac )
call put_global_real_att(ncFileID, "source_damping_rate"    , source_damping_rate )
call put_global_real_att(ncFileID, "source_diurnal_rel_amp" , source_diurnal_rel_amp )
call put_global_real_att(ncFileID, "source_phase_noise"     , source_phase_noise )

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!>@todo
!> this same location attr and write code is in the oned locations module.
!> this code should call it to fill in the locations attrs and data.

! call nc_write_location_atts(ncFileID)
! call nc_get_location_varids(ncFileID, LocationVarID)
! do i=1, num_grid_points
!    call nc_write_location(ncFileID, LocationVarID, loc, index)
! enddo

call nc_check(nf90_def_dim(ncid=ncFileID, name="location", len=int(num_grid_points,i4), dimid = LocationDimID), &
              'nc_write_model_atts', 'nf90_def_dim location') 
call nc_check(NF90_def_var(ncFileID, name="location", xtype=nf90_double, dimids = LocationDimID, varid=LocationVarID) , &
             'nc_write_model_atts', 'nf90_def_var location')

call nc_check(nf90_put_att(ncFileID, LocationVarID, "short_name", trim(adjustl(LocationLName))), &
              'nc_write_model_atts', 'nf90_put_att short_name')
call nc_check(nf90_put_att(ncFileID, LocationVarID, "long_name", "location on a unit circle"), &
              'nc_write_model_atts', 'nf90_put_att long_name')
call nc_check(nf90_put_att(ncFileID, LocationVarID, "dimension", LocationDims ), &
              'nc_write_model_atts', 'nf90_put_att dimension')
call nc_check(nf90_put_att(ncFileID, LocationVarID, "valid_range", (/ 0.0_r8, 1.0_r8 /)), &
              'nc_write_model_atts', 'nf90_put_att valid_range')

! Leave define mode so we can fill
call nc_check(nf90_enddef(ncfileID), 'nc_write_model_atts', 'nf90_enddef')

! Fill the location variable
do i = 1,num_grid_points
   lctn = state_loc(i)
   call nc_check(nf90_put_var(ncFileID, LocationVarID, get_location(lctn), (/ i /) ), &
                 'nc_write_model_atts', 'nf90_put_var LocationVarID 2')
enddo
!--------------------------------------------------------------------
!--------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'nf90_sync')

end function nc_write_model_atts

!--------------------------------------------------------------------

subroutine put_global_char_att(ncid, name, val)

use netcdf

integer,          intent(in) :: ncid
character(len=*), intent(in) :: name
character(len=*), intent(in) :: val

integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, name, val)
call nc_check(ret, 'put_global_char_att', 'adding the global attribute: '//trim(name))

end subroutine put_global_char_att

!--------------------------------------------------------------------

subroutine put_global_real_att(ncid, name, val)

use netcdf

integer,          intent(in) :: ncid
character(len=*), intent(in) :: name
real(r8),         intent(in) :: val

integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, name, val)
call nc_check(ret, 'put_global_real_att', 'adding the global attribute: '//trim(name))

end subroutine put_global_real_att

!--------------------------------------------------------------------

subroutine put_global_creation_time(ncid)
integer, intent(in) :: ncid

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

character(len=128) :: str1

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call put_global_char_att(ncid, "creation_date",str1)

end subroutine put_global_creation_time


function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH 30 Apr 2007
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
use typeSizes
use netcdf

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

call error_handler(E_ERR,'nc_write_model_vars','never gets used',source, revision, revdate)
ierr = -1  ! never reached

end function nc_write_model_vars



subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)
!------------------------------------------------------------------
! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with uniform
! small independent perturbations.

type(ensemble_type), intent(inout) :: state_ens_handle
integer,   intent(in) :: ens_size
real(r8),  intent(in)  :: pert_amp
logical,  intent(out) :: interf_provided

integer :: i,j
real(r8) :: avg_wind

interf_provided = .true.

! allocating storage space in ensemble manager
if(.not. allocated(state_ens_handle%vars)) &
        allocate(state_ens_handle%vars(state_ens_handle%num_vars, state_ens_handle%my_num_copies))
call all_copies_to_all_vars(state_ens_handle)

do i=1,num_grid_points

   ! Perturb the tracer concentration
   do j=1,state_ens_handle%my_num_copies
      state_ens_handle%vars(i,j) = random_gaussian(random_seq, state_ens_handle%vars(i,j), state_ens_handle%vars(i,j))
   enddo
   where(state_ens_handle%vars(i,:) < 0.0_r8) state_ens_handle%vars(i,:) = 0.0_r8 

   ! Perturb the source
   do j=1,state_ens_handle%my_num_copies
      state_ens_handle%vars(num_grid_points + i,j) = random_gaussian(random_seq, &
                                                     state_ens_handle%vars(num_grid_points + i,j), &
                                                     state_ens_handle%vars(num_grid_points + i,j))
   enddo
   where(state_ens_handle%vars(num_grid_points + i,:) < 0.0_r8) state_ens_handle%vars(num_grid_points + i,:) = 0.0_r8 

   ! Perturb the u field
   do j=1,state_ens_handle%my_num_copies

      ! Find the average value of the wind field for the base
      avg_wind = sum(state_ens_handle%vars(2*num_grid_points + i:3*num_grid_points,j)) / num_grid_points
      ! Get a random draw to get 
      state_ens_handle%vars(2*num_grid_points + i,j) = random_gaussian(random_seq, 0.05_r8, avg_wind)
   enddo
   where(state_ens_handle%vars(2*num_grid_points + i,:) < 0.0_r8) state_ens_handle%vars(2*num_grid_points + i,:) = 0.0_r8 


   ! NOT Perturbing the mean_source field
   ! OLD pert_state(3*num_grid_points + i) = random_gaussian(random_seq, state(3*num_grid_points + i), 0.2_r8)
   ! NEW state_ens_handle%vars(3*num_grid_points + i,j) = &
   !                        random_gaussian(random_seq, state_ens_handle%vars(3*num_grid_points + i,j), 0.2_r8)


   ! NOT Perturbing the source_phase_offset field ONLY if the mean_source is non-zero
   ! OLD if(mean_source(i) /= 0.0_r8) pert_state(4*num_grid_points + i) = &
   !             random_gaussian(random_seq, state(4*num_grid_points + i), 0.4_r8)

enddo 

call all_vars_to_all_copies(state_ens_handle)
! deallocate whole state storage
if(.not. get_allow_transpose(state_ens_handle)) deallocate(state_ens_handle%vars)

end subroutine pert_model_copies

!--------------------------------------------------------------------

!> Pass through to the code in the locations module ... 
!> state_handle not needed in this application

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, &
                         obs_kind, num_close, close_ind, dist, state_handle)

type(ensemble_type),         intent(in)     :: state_handle
type(get_close_type),        intent(in)     :: gc
type(location_type),         intent(inout)  :: base_obs_loc, obs_loc(:)
integer,                     intent(in)     :: base_obs_kind, obs_kind(:)
integer,                     intent(out)    :: num_close, close_ind(:)
real(r8),                    intent(out)    :: dist(:)

call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                          num_close, close_ind, dist)

end subroutine get_close_obs


!--------------------------------------------------------------------

!> Unused in this model.

subroutine vert_convert(state_handle, location, obs_kind, istatus)

type(ensemble_type), intent(in)  :: state_handle
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_kind
integer,             intent(out) :: istatus

istatus = 0

end subroutine vert_convert

!--------------------------------------------------------------------

!> pass the vertical localization coordinate to assim_tools_mod

function query_vert_localization_coord()

integer :: query_vert_localization_coord

!> @TODO should define some parameters including something
!> like HAS_NO_VERT for this use.

query_vert_localization_coord = -1

end function query_vert_localization_coord

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
