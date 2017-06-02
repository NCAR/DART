! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

use        types_mod, only : r8, PI, i4, i8

use time_manager_mod, only : time_type, set_time, get_time

use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             nmlfileunit, find_namelist_in_file,           &
                             check_namelist_read, do_output,               &
                             do_nml_file, do_nml_term

use     location_mod,      only : location_type, set_location, get_location, &
                                  get_close_obs, get_close_state, &
                                  convert_vertical_obs, convert_vertical_state

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_sync, &
                                 nc_add_global_creation_time, nc_redef, nc_enddef

use location_io_mod,      only :  nc_write_location_atts, nc_get_location_varids, &
                                  nc_write_location

use default_model_mod,     only : end_model, nc_write_model_vars, init_time

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

integer, parameter :: NVARS = 5
integer :: my_ens_size = 1

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
      call add_dimension_to_variable(dom_id, var_id, 'member', my_ens_size)
      call add_dimension_to_variable(dom_id, var_id, 'location', int(num_grid_points, i4))
   enddo
   
   call finished_adding_domain(dom_id)
endif


end subroutine static_init_model


!------------------------------------------------------------------
!> Initial conditions for simplest 1d Advection.

subroutine init_conditions(x)

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


!------------------------------------------------------------------
!> Single time step advance for simplest 1D advection model.
!> time is for models that need to know it to compute their time tendency. 

subroutine adv_1step(x, time)

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

!>@todo  model_interpolate requires an ensemble handle as one of the
!> parameters to the call.  create a dummy one here.

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



!------------------------------------------------------------------
!> Return the number of items in the state vector

function get_model_size()

integer(i8) :: get_model_size

get_model_size = NVARS*num_grid_points

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

subroutine model_interpolate(state_handle, ens_size, location, itype, expected_val, istatus, x)

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



!------------------------------------------------------------------
!> the smallest time the model should be called to advance

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations



!------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> assoicated location and optionally the state quantity.

subroutine get_state_meta_data(index_in, location, var_type)

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

call nc_redef(ncid)

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
call nc_enddef(ncid)
call nc_write_location(ncid, state_loc, msize)

call nc_sync(ncid)


end subroutine nc_write_model_atts


!------------------------------------------------------------------
!> Perturbs a model state for generating initial ensembles
!> Returning interf_provided means go ahead and do this with uniform
!> small independent perturbations.

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

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


!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
