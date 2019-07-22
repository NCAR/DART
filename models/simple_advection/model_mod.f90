! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> @brief Simple Advection model
!> 
!>  This model is on a periodic one-dimensional domain. A wind field is
!>  modeled using Burger's Equation with an upstream semi-lagrangian
!>  differencing. This diffusive numerical scheme is stable and forcing
!>  is provided by adding in random gaussian noise to each wind grid
!>  variable independently at each timestep. The domain mean value of the
!>  wind is relaxed to a constant fixed value set by the namelist parameter
!>  mean_wind. The random forcing magnitude is set by namelist parameter
!>  wind_random_amp and the damping of the mean wind is controlled by
!>  parameter wind_damping_rate. An Eulerian option with centered in
!>  space differencing is also provided and can be used by setting namelist
!>  parameter lagrangian_for_wind to .false. The Eulerian differencing is
!>  both numerically unstable and subject to shock formation. However, it
!>  can sometimes be made stable in assimilation mode (see recent work by
!>  Majda and collaborators).
!>
!>  The model state includes a single passive tracer that is advected by
!>  the wind field using semi-lagrangian upstream differencing. The state
!>  also includes a tracer source value at each gridpoint. At each time
!>  step, the source is added into the concentration at each gridpoint.
!>  There is also a constant global destruction of tracer that is controlled
!>  by the namelist parameter destruction_rate. The appropriate percentage
!>  of tracer is destroyed at each gridpoint at each timestep.
!>
!>  The model also includes an associated model for the tracer source rate.
!>  At each gridpoint, there is a value of the time mean source rate and
!>  a value of the phase offset for a diurnal component of the source rate.
!>  The diurnal source rate has an amplitude that is proportional to the
!>  source rate (this proportion is controlled by namelist parameter
!>  source_diurnal_rel_amp). At each grid point, the source is the sum
!>  of the source rate plus the appropriate diurnally varying component.
!>  The phase_offset at the gridpoint controls the diurnal phase. The
!>  namelist parameter source_phase_noise controls the amplitude of
!>  random gaussian noise that is added into the source phase at each
!>  time step. If source_phase_noise is zero then the phase offset is
!>  fixed. Finally, the time mean source rate is constant in time in the
!>  present model version. The time mean source rate controls the
!>  amplitude of the diurnal cycle of the tracer source.
!>
!>

module model_mod

use        types_mod, only : r8, PI, i4, i8

use time_manager_mod, only : time_type, set_time, get_time

use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             nmlfileunit, find_namelist_in_file,           &
                             check_namelist_read, do_output,               &
                             do_nml_file, do_nml_term

use  mpi_utilities_mod, only : sum_across_tasks, my_task_id

use     location_mod,      only : location_type, set_location, get_location, &
                                  get_close_obs, get_close_state, &
                                  convert_vertical_obs, convert_vertical_state

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, &
                                 nc_begin_define_mode, nc_end_define_mode

use location_io_mod,      only :  nc_write_location_atts, nc_get_location_varids, &
                                  nc_write_location

use default_model_mod,     only : end_model, nc_write_model_vars, init_time

use     obs_kind_mod, only : QTY_VELOCITY, QTY_TRACER_CONCENTRATION, &
                             QTY_TRACER_SOURCE, QTY_MEAN_SOURCE, QTY_SOURCE_PHASE, &
                             get_name_for_quantity

use random_seq_mod,   only : random_seq_type, init_random_seq, random_gaussian

use ensemble_manager_mod,  only : ensemble_type, init_ensemble_manager, end_ensemble_manager, &
                                  get_my_num_vars, get_my_vars

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain, state_structure_info, &
                                  add_dimension_to_variable, finished_adding_domain

use dart_time_io_mod,      only : read_model_time, write_model_time

implicit none
private


! these routines must be public and you cannot change the
! arguments because they will be called *from* other DART code.

!> required routines with code in this module
public :: get_model_size, &
          get_state_meta_data,  &
          model_interpolate, &
          shortest_time_between_assimilations, &
          static_init_model, &
          init_conditions,    &
          adv_1step, &
          pert_model_copies, &
          nc_write_model_atts

!> required routines where code is in other modules
public :: nc_write_model_vars, &
          init_time, &
          get_close_obs, &
          get_close_state, &
          end_model, &
          convert_vertical_obs, &
          convert_vertical_state, &
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

! if you change NVARS or the order of any of
! these, you must change the 'add_domain()' call
! in static_init_model() below.
integer, parameter :: NVARS = 5
integer, parameter :: CONC      = 1
integer, parameter :: TSOURCE   = 2
integer, parameter :: WIND      = 3
integer, parameter :: MEAN_SRC  = 4
integer, parameter :: SRC_PHASE = 5
integer :: conc_offset, source_offset, wind_offset
integer :: mean_src_offset, src_phase_offset
integer :: model_size

namelist /model_nml/ num_grid_points, grid_spacing_meters, &
                     time_step_days, time_step_seconds, &
                     mean_wind, wind_random_amp, wind_damping_rate, &
                     lagrangian_for_wind, destruction_rate, &
                     source_random_amp_frac, source_damping_rate, &
                     source_diurnal_rel_amp, source_phase_noise, &
                     template_file


! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
type(time_type)                  :: time_step


contains

!------------------------------------------------------------------
!> initialize anything needed. this routine is called once at startup.

subroutine static_init_model()

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

! NVARS items at each grid location
model_size = NVARS * num_grid_points

! Create storage for locations
allocate(state_loc(model_size))

! Define the locations of the model state variables
do i = 1, num_grid_points
   x_loc = (i - 1.0_r8) / num_grid_points
   do j = 1, NVARS
      state_loc(num_grid_points * (j-1) + i) = set_location(x_loc)
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
   call init_random_seq(random_seq, my_task_id()+1)
   random_seq_init = .true.
endif

! Tell the DART I/O routines how large the model data is so they
! can read/write it.
if (template_file /= '') then
   dom_id = add_domain(template_file, NVARS, &
                       (/ 'concentration', &
                          'source       ', &
                          'wind         ', &
                          'mean_source  ', &
                          'source_phase ' /))
else 
   !>@todo FIXME : should not need a template file if initializing members from code

   write(string1, *) 'template file is required for now'
   call error_handler(E_ERR,'static_init_model',string1, source, revision, revdate)

   dom_id = add_domain(NVARS, (/ 'concentration', &
                                 'source       ', &
                                 'wind         ', &
                                 'mean_source  ', &
                                 'source_phase ' /))

   do var_id=1, NVARS
      call add_dimension_to_variable(dom_id, var_id, 'time', 1)
      call add_dimension_to_variable(dom_id, var_id, 'member', 1)
      call add_dimension_to_variable(dom_id, var_id, 'location', int(num_grid_points, i4))
   enddo
   
   call finished_adding_domain(dom_id)
endif

! set the offsets to the start of each variable
! in the state vector.  the indices for each variable
! quantity are x(offset+1 : offset+num_grid_points).

conc_offset      = (CONC-1)      * num_grid_points 
source_offset    = (TSOURCE-1)   * num_grid_points
wind_offset      = (WIND-1)      * num_grid_points
mean_src_offset  = (MEAN_SRC-1)  * num_grid_points
src_phase_offset = (SRC_PHASE-1) * num_grid_points

end subroutine static_init_model


!------------------------------------------------------------------
!> Initial conditions for simplest 1d Advection.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

! Start by zeroing all
x = 0.0_r8

! Tracer concentration
x(conc_offset     +1 : conc_offset     +num_grid_points) = 0.0_r8

! Initial source 
x(source_offset   +1 : source_offset   +num_grid_points) = mean_source

! U wind velocity; units are meters/second
x(wind_offset     +1 : wind_offset     +num_grid_points) = mean_wind

! Time mean source
x(mean_src_offset +1 : mean_src_offset +num_grid_points) = mean_source

! Source phase offset 
x(src_phase_offset+1 : src_phase_offset+num_grid_points) = source_phase_offset

end subroutine init_conditions


!------------------------------------------------------------------
!> Single time step advance for simplest 1D advection model.
!> time is for models that need to know it to compute their time tendency. 

subroutine adv_1step(x, time)

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

integer :: next, prev, seconds, days
integer, parameter :: ens_size = 1
type(location_type) :: source_loc
real(r8) :: lctn, source_location, old_u_mean, new_u_mean
real(r8) :: du_dx, dt_seconds, t_phase, phase, random_src
type(ensemble_type) :: temp_handle

integer(i8) :: i
integer  :: istatus(ens_size)
real(r8) :: new_x(size(x),ens_size)

!> model_interpolate requires an ensemble handle as one of the
!> parameters to the call.  create a dummy one here.
! the data in this handle is never used. the data passed in as
! the final optional "x" argument is used instead.
call init_ensemble_manager(temp_handle, ens_size, int(model_size,i8))

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
    
   source_location = lctn - x(wind_offset + i) * dt_seconds / &
      domain_width_meters

   if(source_location > 1.0_r8) &
      source_location = source_location - int(source_location)

   if(source_location < 0.0_r8) &
      source_location = source_location - int(source_location) + 1.0_r8

   source_loc = set_location(source_location)
   call model_interpolate(temp_handle, ens_size, source_loc, QTY_TRACER_CONCENTRATION, new_x(conc_offset + i,:), istatus, x)  

   ! Following line does lagangian du

   if(lagrangian_for_wind) &
      call model_interpolate(temp_handle, ens_size, source_loc, QTY_VELOCITY, new_x(wind_offset + i,:), istatus, x)  
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
      du_dx = (x(wind_offset + next) - x(wind_offset + prev)) / &
                                       (2.0_r8 * grid_spacing_meters)
      new_x(wind_offset + i,1) = x(wind_offset + i) + &
                                 x(wind_offset + i) * du_dx * dt_seconds
   enddo
endif
!---- End Eulerian block


! Now add in the source contribution and put concentration and wind back into inout x
! Source is in units of .../second
do i = 1, num_grid_points
   x(conc_offset + i) = new_x(conc_offset + i,1) + x(source_offset + i) * dt_seconds
   ! Also copy over the new velocity
   x(wind_offset + i) = new_x(wind_offset + i,1)
enddo

! Now do the destruction rate: Units are fraction destroyed per second
do i = 1, num_grid_points
   if(destruction_rate * dt_seconds > 1.0_r8) then
      x(conc_offset + i) = 0.0_r8
   else
      x(conc_offset + i) = x(conc_offset + i) * (1.0_r8 - destruction_rate*dt_seconds)
   endif
enddo

!----- Following block is random walk plus damping to mean for wind
! Random walk for the spatial mean velocity
old_u_mean = sum(x(wind_offset + 1 : wind_offset + num_grid_points)) / num_grid_points

! Add in a random walk to the mean
new_u_mean = random_gaussian(random_seq, old_u_mean, wind_random_amp * dt_seconds)

! Damp the mean back towards base value
! Need an error handler call here or earlier
if(wind_damping_rate*dt_seconds > 1.0_r8) stop
new_u_mean = new_u_mean - wind_damping_rate * dt_seconds * (new_u_mean - mean_wind)

! Substitute the new mean wind
x(wind_offset + 1 : wind_offset + num_grid_points) = &
   x(wind_offset + 1 : wind_offset + num_grid_points) - old_u_mean + new_u_mean

! Add some noise into each wind element
do i = 1, num_grid_points
   x(wind_offset + i) = random_gaussian(random_seq, x(wind_offset + i), &
                                        wind_random_amp * dt_seconds)
enddo

!----- End forced damped wind section


!----- Time tendency with diurnal cycle for sources ---
! Figure out what the time in the model is
call get_time(time, seconds, days)
t_phase = (2.0_r8 * PI * seconds / 86400.0_r8) 
do i = 1, num_grid_points
   phase = t_phase + x(src_phase_offset + i)
   x(source_offset + i) = x(source_offset + i) &
      + x(mean_src_offset + i) * source_diurnal_rel_amp * cos(phase)
      !!!+ x(mean_src_offset + i) * source_diurnal_rel_amp * cos(phase) * dt_seconds
   ! Also add in some random walk
   random_src = random_gaussian(random_seq, 0.0_r8, source_random_amp(i))
   x(source_offset + i) = x(source_offset + i) + random_src * dt_seconds
   if(x(source_offset + i) < 0.0_r8) x(source_offset + i) = 0.0_r8
   ! Finally, damp back towards the base value
   ! Need an error handler call here
   if(source_damping_rate*dt_seconds > 1.0_r8) stop
   x(source_offset + i) = x(source_offset + i) - &
      source_damping_rate*dt_seconds * (x(source_offset + i) - x(mean_src_offset + i))
enddo

!----- End sources time tendency -----

! Time tendency for source model variables mean_source and source_phase_offset 
! is 0 for now. Might want to add in some process noise for filter advance.

! Process noise test for source_phase_offset
do i = 1, num_grid_points
   x(src_phase_offset + i) = random_gaussian(random_seq, x(src_phase_offset + i), &
      source_phase_noise*dt_seconds)
enddo

call end_ensemble_manager(temp_handle)

end subroutine adv_1step



!------------------------------------------------------------------
!> Return the number of items in the state vector

function get_model_size()

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size


!------------------------------------------------------------------
!> this is a non-standard interface to model_interpolate.  it allows
!> an optional last argument, which if present should be a model state
!> and is used instead of the state handle.  this model_mod uses interpolate
!> during the adv_1step() when there's no state handle available; a naked
!> state is passed in to be advanced.    (having optional args is ok with
!> filter - it just never passes in that final arg.)
!>
!> Interpolates from ensemble of states to the location.
!> If x is present, Interpolates from state vector x to the location. 
!> This code supports three obs types: concentration, source, and u

subroutine model_interpolate(state_handle, ens_size, location, iqty, expected_val, istatus, x)

type(ensemble_type),  intent(in)  :: state_handle
integer,              intent(in)  :: ens_size
type(location_type),  intent(in)  :: location
integer,              intent(in)  :: iqty
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
if(iqty == QTY_TRACER_CONCENTRATION) then
   offset = conc_offset
else if(iqty == QTY_TRACER_SOURCE) then
   offset = source_offset
else if(iqty == QTY_VELOCITY) then
   offset = wind_offset
else
   ! technically this should just set istatus to a positive value and return
   ! because it shouldn't be a fatal error to try to interpolate a quantity
   ! that is not in the state vector. however for this model we only expect
   ! to get observations of types in the state, so end the program if not.
   write(string1, *) 'quantity ', iqty, ' ('//trim(get_name_for_quantity(iqty))// &
                     ') is not supported in model_interpolate'
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



!------------------------------------------------------------------
!> the smallest time the model should be called to advance

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations



!------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location and optionally the state quantity.

subroutine get_state_meta_data(index_in, location, var_qty)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_qty

integer :: var_qty_index, var_loc_index

! Variable types
var_qty_index = (index_in - 1) / num_grid_points + 1
var_loc_index = index_in - (var_qty_index - 1)*num_grid_points

if(present(var_qty)) then
   if(var_qty_index == CONC) then
      var_qty = QTY_TRACER_CONCENTRATION
   else if(var_qty_index == TSOURCE) then
      var_qty = QTY_TRACER_SOURCE
   else if(var_qty_index == WIND) then
      var_qty = QTY_VELOCITY
   else if(var_qty_index == MEAN_SRC) then
      var_qty = QTY_MEAN_SOURCE
   else if(var_qty_index == SRC_PHASE) then
      var_qty = QTY_SOURCE_PHASE
   endif
endif

location = state_loc(var_loc_index)

end subroutine get_state_meta_data


!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

integer :: msize

! other parts of the dart system will write the state into the file
! so this routine just needs to write any model-specific
! attributes it wants to record.

msize = int(num_grid_points, i4)

! Write Global Attributes

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model_revision", revision )
call nc_add_global_attribute(ncid, "model_revdate", revdate )

call nc_add_global_attribute(ncid, "model", "simple_advection")

call nc_add_global_attribute(ncid, "mean_wind"              , mean_wind )
call nc_add_global_attribute(ncid, "wind_random_amp"        , wind_random_amp )
call nc_add_global_attribute(ncid, "wind_damping_rate"      , wind_damping_rate )
call nc_add_global_attribute(ncid, "destruction_rate"       , destruction_rate )
call nc_add_global_attribute(ncid, "source_random_amp_frac" , source_random_amp_frac )
call nc_add_global_attribute(ncid, "source_damping_rate"    , source_damping_rate )
call nc_add_global_attribute(ncid, "source_diurnal_rel_amp" , source_diurnal_rel_amp )
call nc_add_global_attribute(ncid, "source_phase_noise"     , source_phase_noise )

call nc_write_location_atts(ncid, msize)
call nc_end_define_mode(ncid)
call nc_write_location(ncid, state_loc, msize)

call nc_synchronize_file(ncid)


end subroutine nc_write_model_atts


!------------------------------------------------------------------
!> Perturbs a model state for generating initial ensembles.
!> Returning interf_provided .true. means this code has 
!> added uniform small independent perturbations to a
!> single ensemble member to generate the full ensemble.

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,   intent(in) :: ens_size
real(r8),  intent(in) :: pert_amp
logical,  intent(out) :: interf_provided

integer :: i,j, num_my_grid_points, my_qty
integer(i8), allocatable :: my_grid_points(:)
real(r8) :: avg_wind, localsum

interf_provided = .true.

! if we are running with more than 1 task, then
! we have all the ensemble members for a subset of
! the model state.  which variables we have are determined
! by looking at the global index number into the state vector.

! how many grid points does my task have to work on?
! and what are their indices into the full state vector?
num_my_grid_points = get_my_num_vars(state_ens_handle)
allocate(my_grid_points(num_my_grid_points))
call get_my_vars(state_ens_handle, my_grid_points)

! we also want to compute the average wind field before we start.
! find all the wind values on our task and then sum
! them with wind values on other tasks.

localsum = 0.0_r8
do i=1,num_my_grid_points

   ! Variable quantities
   my_qty = (my_grid_points(i) - 1) / num_grid_points + 1

   if(my_qty == WIND) then
      localsum = localsum + state_ens_handle%copies(1, i)
   endif
enddo

call sum_across_tasks(localsum, avg_wind)
avg_wind = avg_wind / num_grid_points
   
! and now we're ready to perturb the ensemble members
! use the global index into the state vector to see what
! quantities we have.
do i=1,num_my_grid_points

   ! Variable quantities
   my_qty = (my_grid_points(i) - 1) / num_grid_points + 1

   if(my_qty == CONC) then
      ! Perturb the tracer concentration
      do j=1,ens_size
         state_ens_handle%copies(j, i) = random_gaussian(random_seq, &
                                         state_ens_handle%copies(j, i), &
                                         state_ens_handle%copies(j, i))
      enddo
      where(state_ens_handle%copies(1:ens_size, i) < 0.0_r8) &
            state_ens_handle%copies(1:ens_size, i) = 0.0_r8 
   
   else if(my_qty == TSOURCE) then
      ! Perturb the source
      do j=1,ens_size
         state_ens_handle%copies(j, i) = random_gaussian(random_seq, &
                                         state_ens_handle%copies(j, i), &
                                         state_ens_handle%copies(j, i))
      enddo
      where(state_ens_handle%copies(1:ens_size, i) < 0.0_r8) &
            state_ens_handle%copies(1:ens_size, i) = 0.0_r8 

   else if(my_qty == WIND) then
      ! Perturb the u field using the average wind computed above
      do j=1,ens_size
         state_ens_handle%copies(j, i) = random_gaussian(random_seq, 0.05_r8, avg_wind)
      enddo
      where(state_ens_handle%copies(1:ens_size, i) < 0.0_r8) &
            state_ens_handle%copies(1:ens_size, i) = 0.0_r8 

   else if(my_qty == MEAN_SRC) then
      ! NOT Perturbing the mean_source field. 
      ! comment in the following lines to perturb the source mean.
      ! do j=1,ens_size
      !   state_ens_handle%copies(j, i) = random_gaussian(random_seq, &
      !                                   state_ens_handle%copies(j, i), 0.2_r8)
      ! enddo

   else if(my_qty == SRC_PHASE) then

      ! NOT Perturbing the source_phase_offset field 
      !  even if commented in, only perturb if the mean_source is non-zero
      ! if (state_ens_handle%copies(j, i) /= 0.0_r8) then
      !    do j=1,ens_size
      !       state_ens_handle%copies(j, i) = random_gaussian(random_seq, &
      !                                       state_ens_handle%copies(j, i), 0.4_r8)
      !    enddo
      ! endif

   endif

enddo 

deallocate(my_grid_points)

end subroutine pert_model_copies


!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
