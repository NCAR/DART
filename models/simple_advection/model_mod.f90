! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

use        types_mod, only : r8, PI
use time_manager_mod, only : time_type, set_time, get_time
use     location_mod, only : location_type, set_location, get_location,  &
                             LocationDims, LocationName, LocationLName,  &
                             get_close_maxdist_init, get_close_obs_init, &
                             get_close_obs

use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             nmlfileunit, find_namelist_in_file,           &
                             check_namelist_read, nc_check, do_output,     &
                             do_nml_file, do_nml_term

use     obs_kind_mod, only : KIND_VELOCITY, KIND_TRACER_CONCENTRATION, &
                             KIND_TRACER_SOURCE, KIND_MEAN_SOURCE, KIND_SOURCE_PHASE

use random_seq_mod,   only : random_seq_type, init_random_seq, random_gaussian

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
          pert_model_state, &
          get_close_maxdist_init, &
          get_close_obs_init, &
          get_close_obs, &
          ens_mean_for_model

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! Simplest 1D advection model with spatially-constant wind

! Module storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq
logical :: random_seq_init = .false.
logical :: verbose = .false.

! Permanent static initialized information about domain width
real(r8) :: domain_width_meters 


! Base tracer source rate in amount (unitless) / second
real(r8), allocatable :: mean_source(:)
real(r8), allocatable :: source_phase_offset(:)
real(r8), allocatable :: source_random_amp(:)

!---------------------------------------------------------------

! Namelist with default values
!
integer  :: num_grid_points         = 10
!!!integer  :: grid_spacing_meters     = 100 * 1000
integer  :: grid_spacing_meters     = 100000
integer  :: time_step_days          = 0
integer  :: time_step_seconds       = 3600

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

logical   :: output_state_vector    = .false.

namelist /model_nml/ num_grid_points, grid_spacing_meters, &
                     time_step_days, time_step_seconds, &
                     mean_wind, wind_random_amp, wind_damping_rate, &
                     lagrangian_for_wind, destruction_rate, &
                     source_random_amp_frac, source_damping_rate, &
                     source_diurnal_rel_amp, source_phase_noise, &
                     output_state_vector

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
integer  :: i, iunit, io, j

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
allocate(state_loc(5*num_grid_points))

! Define the locations of the model state variables
do i = 1, num_grid_points
   x_loc = (i - 1.0_r8) / num_grid_points
   do j = 0, 4
      state_loc(num_grid_points * j + i) = set_location(x_loc)
   end do
end do

! The time_step in terms of a time type must also be initialized.
time_step = set_time(time_step_seconds, time_step_days)

! Allocate storage for the base tracer source and initialize it
! Might want to put the source distribution into the namelist somehow?
allocate(mean_source(num_grid_points), source_phase_offset(num_grid_points), &
   source_random_amp(num_grid_points))

! Set the base value for the mean source; This case has a single enhanced source
mean_source    = 0.1_r8
mean_source(1) = 1.0_r8

! A rapidly varying source
!!!mean_source = 0.1_r8
!!!do i = 1, 5, 2
   !!!mean_source(i) = 1.0_r8
!!!end do

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

integer :: i, istatus, next, prev, seconds, days
type(location_type) :: this_loc, source_loc
real(r8) :: lctn, source_location, new_x(size(x)), old_u_mean, new_u_mean
real(r8) :: du_dx, dt_seconds, t_phase, phase, random_src

! State is concentrations (num_grid_points), source(num_grid_points), u(num_grid_points),
! mean_source(num_grid_points), and source_phase_offset(num_grid_points)
! all dimensioned num_grid_points.


! For the concentration do a linear interpolated upstream distance
! Also have an option to do upstream for wind (controlled by namelist; see below).
! Velocity is in meters/second. Need time_step in seconds
dt_seconds = time_step_seconds + time_step_days * 24.0_r8 * 3600.0_r8


do i = 1, num_grid_points

   ! Find the location of the target point
   call get_state_meta_data(i, this_loc)
   lctn = get_location(this_loc)

   ! Find source point location: Velocity is in meters/second
   ! Figure out meters to move and then convert to fraction of domain
    
   source_location = lctn - x(2*num_grid_points + i) * dt_seconds / &
      domain_width_meters

   if(source_location > 1.0_r8) &
      source_location = source_location - int(source_location)

   if(source_location < 0.0_r8) &
      source_location = source_location - int(source_location) + 1.0_r8

   source_loc = set_location(source_location)
   call model_interpolate(x, source_loc, KIND_TRACER_CONCENTRATION, &
      new_x(i), istatus)  
   ! Following line does lagangian du
   if(lagrangian_for_wind) call model_interpolate(x, source_loc, KIND_VELOCITY, &
      new_x(2*num_grid_points + i), istatus)  
end do


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
      new_x(2*num_grid_points + i) = x(2*num_grid_points + i) + &
         x(2*num_grid_points + i) * du_dx * dt_seconds
   end do
endif
!---- End Eulerian block


! Now add in the source contribution and put concentration and wind back into inout x
! Source is in units of .../second
do i = 1, num_grid_points
   x(i) = new_x(i) + x(num_grid_points + i) * dt_seconds
   ! Also copy over the new velocity
   x(2*num_grid_points + i) = new_x(2*num_grid_points + i)
end do

! Now do the destruction rate: Units are fraction destroyed per second
do i = 1, num_grid_points
   if(destruction_rate * dt_seconds > 1.0_r8) then
      x(i) = 0.0_r8
   else
      x(i) = x(i) * (1.0_r8 - destruction_rate*dt_seconds)
   endif
end do

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
end do

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
end do

!----- End sources time tendency -----

! Time tendency for source model variables mean_source and source_phase_offset 
! is 0 for now. Might want to add in some process noise for filter advance.

! Process noise test for source_phase_offset
do i = 1, num_grid_points
   x(4*num_grid_points + i) = random_gaussian(random_seq, x(4*num_grid_points + i), &
      source_phase_noise*dt_seconds)
end do


end subroutine adv_1step



function get_model_size()
!------------------------------------------------------------------
! function get_model_size()
!
! Returns size of model

integer :: get_model_size

get_model_size = 5*num_grid_points

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



subroutine model_interpolate(x, location, itype, obs_val, istatus)
!------------------------------------------------------------------
!
! Interpolates from state vector x to the location. 
! Have three types, concentration, source, and u


real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

integer :: lower_index, upper_index, offset
real(r8) :: lctn, lctnfrac

! All forward operators supported
istatus = 0

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
if(itype == KIND_TRACER_CONCENTRATION) then
   offset = 0
else if(itype == KIND_TRACER_SOURCE) then
   offset = num_grid_points
else if(itype == KIND_VELOCITY) then
   offset = 2*num_grid_points
else
   write(*, *) 'itype is not supported in model_interpolate', itype
   stop
endif

! Add the offsets into the lower and upper indices
lower_index = lower_index + offset
upper_index = upper_index + offset

! Do the interpolation and return
obs_val = (1.0_r8 - lctnfrac) * x(lower_index) + lctnfrac * x(upper_index)

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



subroutine get_state_meta_data(index_in, location, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?


integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

integer :: var_type_index, var_loc_index

! Three variable types
var_type_index = (index_in - 1) / num_grid_points + 1
var_loc_index = index_in - (var_type_index - 1)*num_grid_points

if(present(var_type)) then
   if(var_type_index == 1) then
      var_type = KIND_TRACER_CONCENTRATION
   else if(var_type_index == 2) then
      var_type = KIND_TRACER_SOURCE
   else if(var_type_index == 3) then
      var_type = KIND_VELOCITY
   else if(var_type_index == 4) then
      var_type = KIND_MEAN_SOURCE
   else if(var_type_index == 5) then
      var_type = KIND_SOURCE_PHASE
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



function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH 30 Apr 2007
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode 
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

use typeSizes
use netcdf

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

!--------------------------------------------------------------------
! General netCDF variables
!--------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

!--------------------------------------------------------------------
! netCDF variables for Location
!--------------------------------------------------------------------

integer :: LocationDimID, LocationVarID
integer :: StateVarDimID, StateVarVarID
integer :: StateVarID, MemberDimID, TimeDimID
integer :: ConcVarId, SourceVarID, WindVarID, Mean_SourceVarID, Source_PhaseVarID

!--------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer             :: i
type(location_type) :: lctn 
character(len=128)  :: filename

ierr = 0                      ! assume normal termination

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file 
!--------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID), &
              'nc_write_model_atts', 'inquire, '//trim(filename))
call nc_check(nf90_sync(ncFileID), & ! Ensure netCDF file is current
              'nc_write_model_atts', 'sync, '//trim(filename))
call nc_check(nf90_Redef(ncFileID), &
              'nc_write_model_atts', 'redef, '//trim(filename))

!--------------------------------------------------------------------
! Determine ID's from stuff already in the netCDF file
!--------------------------------------------------------------------

! make sure time is unlimited dimid

call nc_check(nf90_inq_dimid(ncFileID,"copy",dimid=MemberDimID), &
              'nc_write_model_atts', 'inq_dimid copy, '//trim(filename))
call nc_check(nf90_inq_dimid(ncFileID,"time",dimid=TimeDimID), &
              'nc_write_model_atts', 'inq_dimid time, '//trim(filename))

!--------------------------------------------------------------------
! Write Global Attributes 
!--------------------------------------------------------------------
call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1), &
              'nc_write_model_atts', 'put_att creation_date, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source", source ), &
              'nc_write_model_atts', 'put_att model_source, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision", revision ), &
              'nc_write_model_atts', 'put_att model_revision, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate", revdate ), &
              'nc_write_model_atts', 'put_att model_revdate, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model", "simple_advection"), &
              'nc_write_model_atts', 'put_att model, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "destruction_rate", destruction_rate ), &
              'nc_write_model_atts', 'put_att destruction_rate, '//trim(filename))

! Define the model size, state variable dimension ... whatever ...
call nc_check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                  len=get_model_size(), dimid = StateVarDimID), &
                'nc_write_model_atts', 'def_dim StateVariable, '//trim(filename))

!----------------------------------------------------------------------------
! Define either the "state vector" variables -OR- the "prognostic" variables.
!----------------------------------------------------------------------------

if ( output_state_vector ) then

   ! Define the Location Variable and add Attributes
   ! Some of the atts come from location_mod (via the USE: stmnt)
   ! CF standards for Locations:
   ! http://www.cgd.ucar.edu/cms/eaton/cf-metadata/index.html

   call nc_check(NF90_def_var(ncFileID, name=trim(adjustl(LocationName)), &
                 xtype=nf90_double, dimids = StateVarDimID, varid=LocationVarID), &
              'nc_write_model_atts', 'check, '//trim(LocationName)//', '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, LocationVarID, "long_name", &
                                     trim(adjustl(LocationLName))), &
               'nc_write_model_atts', 'put_att long_name, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, LocationVarID, "dimension", LocationDims), &
              'nc_write_model_atts', 'put_att dimension, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, LocationVarID, "units", "nondimensional"), &
              'nc_write_model_atts', 'put_att units, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, LocationVarID, "valid_range", (/ 0.0_r8, 1.0_r8 /)), &
              'nc_write_model_atts', 'put_att valid_range, '//trim(filename))

   ! Define the state vector coordinate variable
   call nc_check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
                 dimids=StateVarDimID, varid=StateVarVarID), &
                'nc_write_model_atts', 'def_var StateVariable, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"), &
                 'nc_write_model_atts', 'put_att long_name, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical"), &
                 'nc_write_model_atts', 'put_att units, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, num_grid_points*3 /)), &
                 'nc_write_model_atts', 'put_att valid_range, '//trim(filename)) 
   
   ! Define the actual state vector
   call nc_check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_double, &
              dimids = (/ StateVarDimID, MemberDimID, TimeDimID /), varid=StateVarID), &
              'nc_write_model_atts', 'def_var state, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"), &
                 'nc_write_model_atts', 'put_att long_name, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, StateVarID, "_FillValue", NF90_FILL_DOUBLE), &
                 'nc_write_model_atts', 'put_att fillValue, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, StateVarID, "missing_value", NF90_FILL_DOUBLE), &
                 'nc_write_model_atts', 'put_att missing, '//trim(filename))
   
   !--------------------------------------------------------------------
   ! Leave define mode so we can fill
   !--------------------------------------------------------------------
   call nc_check(nf90_enddef(ncfileID), &
                 'nc_write_model_atts', 'enddef, '//trim(filename))
   
   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,num_grid_points*3) /) ), &
                 'nc_write_model_atts', 'put_var state variable coordinate, '//trim(filename))
   
   ! Fill the location variable
   do i = 1,get_model_size()
      call get_state_meta_data(i,lctn)
      call nc_check(nf90_put_var(ncFileID, LocationVarID, get_location(lctn), (/ i /) ), &
                 'nc_write_model_atts', 'check locationVarId, '//trim(filename))
   enddo

else

   !----------------------------------------------------------------------------
   ! We need to map the state vector into the prognostic variables.
   !----------------------------------------------------------------------------
   ! concentration    (                  1 : 1*num_grid_points)
   ! source           (  num_grid_points+1 : 2*num_grid_points)
   ! wind             (2*num_grid_points+1 : 3*num_grid_points)
   ! mean_source      (3*num_grid_points+1 : 4*num_grid_points)
   ! source_phase     (4*num_grid_points+1 : 5*num_grid_points)
   !----------------------------------------------------------------------------
   ! Create the (empty) Variables and the Attributes
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_dim(ncid=ncFileID, name=trim(adjustl(LocationName)), &
                           len=num_grid_points, dimid = LocationDimID), &
                          'nc_write_model_atts', 'def_dim locations, '//trim(filename))

   ! Define the Location Variable and add Attributes
   ! Some of the atts come from location_mod (via the USE: stmnt)
   ! CF standards for Locations:
   ! http://www.cgd.ucar.edu/cms/eaton/cf-metadata/index.html

   call nc_check(NF90_def_var(ncFileID, name=trim(adjustl(LocationName)), &
                 xtype=nf90_double, dimids = LocationDimID, varid=LocationVarID), &
              'nc_write_model_atts', 'check, '//trim(LocationName)//', '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, LocationVarID, "long_name", &
                                     trim(adjustl(LocationLName))), &
               'nc_write_model_atts', 'put_att long_name, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, LocationVarID, "dimension", LocationDims), &
              'nc_write_model_atts', 'put_att dimension, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, LocationVarID, "units", "nondimensional"), &
              'nc_write_model_atts', 'put_att units, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, LocationVarID, "valid_range", (/ 0.0_r8, 1.0_r8 /)), &
              'nc_write_model_atts', 'put_att valid_range, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, LocationVarID, "axis", "X"), &
                 'nc_write_model_atts', 'locations:axis, '//trim(filename)) 

   ! Define the prognostic variables - concentration

   call nc_check(nf90_def_var(ncid=ncFileID,name="concentration", xtype=nf90_double, &
                 dimids = (/ LocationDimID, MemberDimID, TimeDimID /), varid=ConcVarID), &
                'nc_write_model_atts', 'def_var concentration, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ConcVarID, "long_name", "tracer concentration"), &
                 'nc_write_model_atts', 'conc:long_name, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ConcVarID, "units",     "mass"), &
                 'nc_write_model_atts', 'conc:units, '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID, ConcVarID, "valid_range", (/ 0.0_r8, 1.0_r8 /)), &
!                'nc_write_model_atts', 'conc:valid_range, '//trim(filename)) 
   call nc_check(nf90_put_att(ncFileID, ConcVarID, "_FillValue", NF90_FILL_DOUBLE), &
                 'nc_write_model_atts', 'conc:fillValue, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ConcVarID, "missing_value", NF90_FILL_DOUBLE), &
                 'nc_write_model_atts', 'conc:missing, '//trim(filename))

   ! Define the prognostic variables - source

   call nc_check(nf90_def_var(ncid=ncFileID,name="source", xtype=nf90_double, &
                 dimids = (/ LocationDimID, MemberDimID, TimeDimID /), varid=SourceVarID), &
                'nc_write_model_atts', 'def_var source, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SourceVarID, "long_name", "source"), &
                 'nc_write_model_atts', 'source:long_name, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SourceVarID, "units",     "mass/timestep"), &
                 'nc_write_model_atts', 'source:units, '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID, SourceVarID, "valid_range", (/ 0.0_r8, 1.0_r8 /)), &
!                'nc_write_model_atts', 'source:valid_range, '//trim(filename)) 
   call nc_check(nf90_put_att(ncFileID, SourceVarID, "_FillValue", NF90_FILL_DOUBLE), &
                 'nc_write_model_atts', 'source:fillValue, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SourceVarID, "missing_value", NF90_FILL_DOUBLE), &
                 'nc_write_model_atts', 'source:missing, '//trim(filename))

   ! Define the prognostic variables - wind

   call nc_check(nf90_def_var(ncid=ncFileID,name="wind", xtype=nf90_double, &
                 dimids = (/ LocationDimID, MemberDimID, TimeDimID /), varid=WindVarID), &
                'nc_write_model_atts', 'def_var wind, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, WindVarID, "long_name", "wind"), &
                 'nc_write_model_atts', 'wind:long_name, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, WindVarID, "units",  "gridpoints/timestep"), &
                 'nc_write_model_atts', 'wind:units, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, WindVarID, "_FillValue", NF90_FILL_DOUBLE), &
                 'nc_write_model_atts', 'wind:fillValue, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, WindVarID, "missing_value", NF90_FILL_DOUBLE), &
                 'nc_write_model_atts', 'wind:missing, '//trim(filename))

   ! Define the prognostic variables - mean_source

   call nc_check(nf90_def_var(ncid=ncFileID,name="mean_source", xtype=nf90_double, &
                 dimids = (/ LocationDimID, MemberDimID, TimeDimID /), varid=Mean_SourceVarID), &
                'nc_write_model_atts', 'def_var mean_source, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, Mean_SourceVarID, "long_name", "mean_source"), &
                 'nc_write_model_atts', 'mean_source:long_name, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, Mean_SourceVarID, "units",  "mass/timestep"), &
                 'nc_write_model_atts', 'mean_source:units, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, Mean_SourceVarID, "_FillValue", NF90_FILL_DOUBLE), &
                 'nc_write_model_atts', 'mean_source:fillValue, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, Mean_SourceVarID, "missing_value", NF90_FILL_DOUBLE), &
                 'nc_write_model_atts', 'mean_source:missing, '//trim(filename))

   ! Define the prognostic variables - source_phase

   call nc_check(nf90_def_var(ncid=ncFileID,name="source_phase", xtype=nf90_double, &
                 dimids = (/ LocationDimID, MemberDimID, TimeDimID /), varid=Source_PhaseVarID), &
                'nc_write_model_atts', 'def_var source_phase, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, Source_PhaseVarID, "long_name", "source phase"), &
                 'nc_write_model_atts', 'source_phase:long_name, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, Source_PhaseVarID, "units",  "radians"), &
                 'nc_write_model_atts', 'source_phase:units, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, Source_PhaseVarID, "_FillValue", NF90_FILL_DOUBLE), &
                 'nc_write_model_atts', 'source_phase:fillValue, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, Source_PhaseVarID, "missing_value", NF90_FILL_DOUBLE), &
                 'nc_write_model_atts', 'source_phase:missing, '//trim(filename))

   !--------------------------------------------------------------------
   ! Leave define mode so we can fill
   !--------------------------------------------------------------------
   call nc_check(nf90_enddef(ncfileID), &
                 'nc_write_model_atts', 'enddef, '//trim(filename))
   
   ! Fill the locations coordinate variable
   do i = 1,num_grid_points
      call get_state_meta_data(i,lctn)
      call nc_check(nf90_put_var(ncFileID, LocationVarID, get_location(lctn), (/ i /) ), &
                 'nc_write_model_atts', 'put_var locations, '//trim(filename))
   enddo
endif

!--------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!--------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID), &
              'nc_write_model_atts', 'sync, '//trim(filename))

write (*,*)'Model attributes written, netCDF file synched ...'

end function nc_write_model_atts



function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH 30 Apr 2007
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

use typeSizes
use netcdf

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

!--------------------------------------------------------------------
! General netCDF variables
!--------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID, ConcVarID, SourceVarID, WindVarID, Mean_SourceVarID, Source_PhaseVarID

!--------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------

character(len=128) :: filename
integer :: first, last

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID


ierr = 0                      ! assume normal termination

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!--------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID), &
              'nc_write_model_vars', 'inquire, '//trim(filename))

!------------------------------------------------------------------------
! Branch for state vector OR prognostic variables
!------------------------------------------------------------------------

if ( output_state_vector ) then

   call nc_check(NF90_inq_varid(ncFileID, "state", StateVarID), & 
                 'nc_write_model_vars', 'inq_varid state, '//trim(filename))
   call nc_check(NF90_put_var(ncFileID, StateVarID, statevec(1:num_grid_points*3),  &
                 start=(/ 1, copyindex, timeindex /)),  &
                 'nc_write_model_vars', 'put_var state vector, '//trim(filename))

else

   first = 0*num_grid_points+1
   last  = 1*num_grid_points
   call nc_check(NF90_inq_varid(ncFileID, "concentration", ConcVarID), & 
                 'nc_write_model_vars', 'inq_varid concentration, '//trim(filename))
   call nc_check(NF90_put_var(ncFileID, ConcVarID, statevec(first:last),  &
                 start=(/ 1, copyindex, timeindex /)),  &
                 'nc_write_model_vars', 'put_var concentration, '//trim(filename))

   first = 1*num_grid_points+1
   last  = 2*num_grid_points
   call nc_check(NF90_inq_varid(ncFileID, "source", SourceVarID), & 
                 'nc_write_model_vars', 'inq_varid source, '//trim(filename))
   call nc_check(NF90_put_var(ncFileID, SourceVarID, statevec(first:last),  &
                 start=(/ 1, copyindex, timeindex /)),  &
                 'nc_write_model_vars', 'put_var source, '//trim(filename))

   first = 2*num_grid_points+1
   last  = 3*num_grid_points
   call nc_check(NF90_inq_varid(ncFileID, "wind", WindVarID), & 
                 'nc_write_model_vars', 'inq_varid wind, '//trim(filename))
   call nc_check(NF90_put_var(ncFileID, WindVarID, statevec(first:last),  &
                 start=(/ 1, copyindex, timeindex /)),  &
                 'nc_write_model_vars', 'put_var wind, '//trim(filename))

   first = 3*num_grid_points+1
   last  = 4*num_grid_points
   call nc_check(NF90_inq_varid(ncFileID, "mean_source", Mean_SourceVarID), & 
                 'nc_write_model_vars', 'inq_varid mean_source, '//trim(filename))
   call nc_check(NF90_put_var(ncFileID, Mean_SourceVarID, statevec(first:last),  &
                 start=(/ 1, copyindex, timeindex /)),  &
                 'nc_write_model_vars', 'put_var mean_source, '//trim(filename))

   first = 4*num_grid_points+1
   last  = 5*num_grid_points
   call nc_check(NF90_inq_varid(ncFileID, "source_phase", Source_PhaseVarID), & 
                 'nc_write_model_vars', 'inq_varid source_phase, '//trim(filename))
   call nc_check(NF90_put_var(ncFileID, Source_PhaseVarID, statevec(first:last),  &
                 start=(/ 1, copyindex, timeindex /)),  &
                 'nc_write_model_vars', 'put_var source_phase, '//trim(filename))

endif

! write (*,*)'Finished filling variables ...'
call nc_check(nf90_sync(ncFileID), 'nc_write_model_vars', 'sync, '//trim(filename))
! write (*,*)'netCDF file is synched ...'

end function nc_write_model_vars



subroutine pert_model_state(state, pert_state, interf_provided)
!------------------------------------------------------------------
! subroutine pert_model_state(state, pert_state, interf_provided)
!
! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with uniform
! small independent perturbations.

real(r8), intent(in)  :: state(:)
real(r8), intent(out)  :: pert_state(:)
logical,  intent(out) :: interf_provided

integer :: i
real(r8) :: avg_wind

interf_provided = .true.

write(*, *) 'in pert_model_state'

! Need to make sure perturbed states are not centered on true value
! This model is too forgiving in such circumstances
do i = 1, num_grid_points
   ! Perturb the tracer concentration
   pert_state(i) = random_gaussian(random_seq, state(i), state(i))
   if(pert_state(i) < 0.0_r8) pert_state(i) = 0.0_r8
   ! Perturb the source
   pert_state(num_grid_points + i) = random_gaussian(random_seq, &
      state(num_grid_points + i), state(num_grid_points + i))
   if(pert_state(num_grid_points + i) < 0.0_r8) &
      pert_state(num_grid_points + i) = 0.0_r8

   ! Perturb the u field
   ! Find the average value of the wind field for the base
   avg_wind = sum(state(2*num_grid_points + i:3*num_grid_points)) / num_grid_points
   ! Get a random draw to get 
   pert_state(2*num_grid_points + i) = random_gaussian(random_seq, &
      0.05_r8, avg_wind)
   if(pert_state(2*num_grid_points + i) < 0.0_r8) &
      pert_state(2*num_grid_points + i) = 0.0_r8

   ! Perturb the mean_source field
   ! Need to get this into nameslist, or control stuff at the top
   ! Make sure that this is positive definite here
   pert_state(3*num_grid_points + i) = &
      !!!random_gaussian(random_seq, state(3*num_grid_points + i), 0.2_r8)
   !!!if(pert_state(3*num_grid_points + i) < 0.0_r8) pert_state(3*num_grid_points + i) = 0.0_r8
      state(3*num_grid_points + i)

   ! Perturb the source_phase_offset field ONLY if the mean_source is non-zero
   !!!if(mean_source(i) /= 0.0_r8) pert_state(4*num_grid_points + i) = &
         !!!random_gaussian(random_seq, state(4*num_grid_points + i), 0.4_r8)
   pert_state(4*num_grid_points + i) = state(4*num_grid_points + i)
   

end do 

end subroutine pert_model_state



subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! Not used in low-order models

real(r8), intent(in) :: ens_mean(:)

end subroutine ens_mean_for_model

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
