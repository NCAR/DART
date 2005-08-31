! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module model_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 

! Model mod for the single column climate model for Yuqiong Liu. The model
! is a single column of CCM 3.6 combined with a simple land surface model
! column under it. The land surface column has a number of tiles and
! a number of levels. An arbitrary number of atmospheric tracers, the first
! one being q, are carried.

! Modules that are absolutely required for use are listed
use        types_mod, only : r8
use time_manager_mod, only : time_type, set_time, get_time, get_date, &
                             set_date, operator(+)
use     location_mod, only : location_type, set_location
use    utilities_mod, only : file_exist, open_file, check_nml_error, close_file, &
                             register_module, error_handler, E_ERR, E_MSG, logfileunit

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
          model_get_close_states, &
          nc_write_model_atts, &
          nc_write_model_vars, &
          pert_model_state, &
          write_sccm_state, &
          read_sccm_state

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"


! Basic model parameters controlled by nameslist; have defaults

!---------------------------------------------------------------
! Namelist with default values
!

integer  :: kpt        = 2          ! Number of tiles per land surface box
integer  :: msl        = 6          ! Number of soil levels
integer  :: plev       = 18         ! Number of atmospheric levels
integer  :: pcnst      = 1          ! Number of tracers
integer  :: time_step_days = 0
integer  :: time_step_seconds = 1200

namelist /model_nml/ kpt, msl, plev, pcnst, time_step_days, time_step_seconds
!----------------------------------------------------------------

! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
type(time_type) :: time_step
! Global storage for size of model
integer :: model_size


contains

!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------
!
! Called to do one time initialization of the model. 
! For sccm, this entails reading the namelist to get the size and 
! configuration details of the model.

real(r8) :: x_loc
integer  :: i, iunit, ierr, io

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Begin by reading the namelist input
if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(iunit, nml = model_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'model_nml')
   enddo
 11 continue
   call close_file(iunit)
endif

! Record the namelist values used for the run ...
call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
write(logfileunit, nml=model_nml)
write(     *     , nml=model_nml)

! Compute the model state vector size
model_size = (3 + pcnst) * plev + 1 + 4*kpt + 2*kpt*msl
! Create storage for locations
allocate(state_loc(model_size))

! Define the locations of the model state variables
! For now, just make all locations the same using 3D location
! For small model like this assume no localization will be used
do i = 1, model_size
   state_loc(i) =  set_location(0.0_r8, 0.0_r8, 1000.0_r8, 2)
end do

! The time_step in terms of a time type must also be initialized. 
time_step = set_time(time_step_seconds, time_step_days)

end subroutine static_init_model



subroutine init_conditions(x)
!------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

real(r8), intent(out) :: x(:)


end subroutine init_conditions



subroutine adv_1step(x, time)
!------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
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

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

end subroutine adv_1step



function get_model_size()
!------------------------------------------------------------------
!
! Returns the size of the model as an integer.  Model size is computed
! in static_init_model from namelist parameters.

integer :: get_model_size

get_model_size = model_size 

end function get_model_size



subroutine init_time(time)
!------------------------------------------------------------------
!
! Companion interface to init_conditions. Returns a time that is somehow 
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

type(time_type), intent(out) :: time


end subroutine init_time



subroutine model_interpolate(x, location, itype, obs_val, istatus)
!------------------------------------------------------------------
!
! Given a state vector, a location, and a model state variable type,
! interpolates the state variable field to that location and returns
! the value in obs_val. The istatus variable should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is a model specific integer that specifies the type of field (for
! instance temperature, zonal wind component, etc.). In low order
! models that have no notion of types of variables, this argument can
! be ignored. For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observerd), this can be a NULL INTERFACE.

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

! Default for successful return
istatus = 0

end subroutine model_interpolate



function get_model_time_step()
!------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

get_model_time_step = time_step

end function get_model_time_step



subroutine get_state_meta_data(index_in, location, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location.  For the SCCM, we are assuming for now that no
! localization will be used and all the locations are the same. 
! The type is returned to distinguish the different variable types
! with the following values (these should later be changed to
! use the global kind module: JLA):
!
!   1 = u,  2 = v,  3 = t,   4 = tracer,   5 = ps,  
!   6 = h2osno,   7 = h2ocan,   8 = tv,   9 = tg,   
!   10 = h2osoi,   11 = tsoi
!   

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

integer :: g_index, s_index

location = state_loc(index_in)
if(present(var_type)) then
   ! Set indices for the ground part of the model
   g_index = index_in - (3 + pcnst)*plev + 1
   s_index = g_index - 4*kpt 
   if(index_in <= plev) then
      var_type = 1    ! u column comes first
   else if(index_in <= 2*plev) then
      var_type = 2    ! Then v column
   else if(index_in <= 3*plev) then
      var_type = 3    ! Then t column
   else if(index_in <= (3 + pcnst) *plev) then
      var_type = 4    ! Then pcnst tracer columns
   else if(index_in == 4*plev + 1) then
      var_type = 5    ! ps
   else if(g_index <= kpt) then
      var_type = 6    ! Snow water per tile
   else if(g_index <= 2*kpt) then
      var_type = 7    ! Water in canopy
   else if(g_index <= 3*kpt) then
      var_type = 8    ! Temperature of vegetation
   else if(g_index <= 4*kpt) then
      var_type = 9    ! Surface skin temperature
   else if(s_index <= kpt*msl) then
      var_type = 10   ! Soil moisture, by level and tile
   else if(s_index <= 2*kpt*msl) then
      var_type = 11   ! Temperature of soil by level and tile
   else
      call error_handler(E_ERR,'get_state_meta_data',&   
         'index_in is not in range of model size', source, revision, revdate)
   endif

endif


end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

end subroutine end_model



subroutine model_get_close_states(o_loc, radius, inum, indices, dist, x)
!------------------------------------------------------------------
! 
! Computes a list of model state variable indices that are within 
! distance radius of a given location, o_loc. The units of the radius
! and the metric for computing distances is defined by the location module
! that is in  use. The number of state variables that are within radius
! of o_loc is returned in inum. The indices of each of these is 
! returned in indices and the corresponding distance in dist. The model
! state is available in x because it is sometimes required to determine
! the distance (for instance, the current model surface pressure field
! is required to compute the location of state variables in a sigma
! vertical coordinate model). A model can choose to do no computation
! here and return a value of -1 in inum. If this happens, the filter
! will do a naive search through ALL state variables for close states.
! This can work fine in low-order models, but can be far too expensive
! in large models.

type(location_type), intent(in) :: o_loc
real(r8), intent(in) :: radius
integer, intent(out) :: inum, indices(:)
real(r8), intent(out) :: dist(:)
real(r8), intent(in) :: x(:)

! Simplest interface just sets inum to -1 and returns
inum = -1

end subroutine model_get_close_states



function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH Jan 24 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_96 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
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

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

end function nc_write_model_atts



function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
! Writes the model variables to a netCDF file
! TJH 23 May 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_96 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
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

end function nc_write_model_vars



subroutine pert_model_state(state, pert_state, interf_provided)
!------------------------------------------------------------------
!
! Perturbs a model state for generating initial ensembles.
! The perturbed state is returned in pert_state.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding an O(0.1) magnitude perturbation to each
! model state variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

interf_provided = .false.

end subroutine pert_model_state


subroutine write_sccm_state(state, time, delta_seconds, file_id)
!------------------------------------------------------------------
!

real(r8), intent(in) :: state(:)
integer, intent(in) :: delta_seconds, file_id
type(time_type), intent(in) :: time

real(r8) :: atmos(plev), sfc(kpt), soil(msl * kpt)
integer :: g_offset, s_offset, i, date_stamp
integer :: year, month, day, hour, minute, second, days, seconds

! Write the delta seconds
write(file_id, *) "advance-time ", delta_seconds

! First write day as date integer and seconds
call get_date(time, year, month, day, hour, minute, second)
! Create and write the date stamp for the sccm
year = year - year / 100 * 100
date_stamp = year * 10000 + month * 100 + day
write(file_id, *) 'date ', date_stamp
call get_time(time, seconds, days)
write(file_id, *) 'sec ', seconds

! Write out the model size descriptors
write(file_id, *) 'KPT ', kpt
write(file_id, *) 'MSL ', msl
write(file_id, *) 'PLEV ', plev
write(file_id, *) 'PCNST ', pcnst

! Write out the atmospheric t
atmos(:) = state(2*plev + 1 : 3*plev)
write(file_id, *) 't3(plev) ', atmos
! Write out atmospheric u
atmos(:) = state(1 : plev) 
write(file_id, *) 'u3(plev) ', atmos
! Write out atmospheric v
atmos(:) = state(plev + 1 : 2*plev)
write(file_id, *) 'v3(plev) ', atmos
! Write out the tracer fields
do i = 1, pcnst
   atmos(:) = state((2 + i)*plev + 1: (3 + i) * plev)
   write(file_id, *) 'q3(pcnst,plev) ', atmos
end do
! Write out the ps
write(file_id, *) 'ps ', state((3 + pcnst)* plev + 1)


! Offsets for where the ground and soil state variables start in dart
g_offset = (3 + pcnst) * plev + 2
s_offset = g_offset + 4 * kpt
! Snow moisture content is first
sfc = state(g_offset : g_offset + kpt - 1)
write(file_id, *) 'h2osno(kpt) ', sfc
! Canopy moisture next
sfc = state(g_offset + kpt : g_offset + 2*kpt - 1)
write(file_id, *) 'h20can(kpt) ', sfc
! Now soil moisture
soil = state(s_offset : s_offset + msl*kpt - 1)
write(file_id, *) 'h2osoi(kpt, msl) ', soil
! Now temperature of vegetation
sfc = state(g_offset + 2*kpt : g_offset + 3*kpt - 1)
write(file_id, *) 'tv(kpt) ', sfc
! Now skin temperature
sfc = state(g_offset + 3*kpt : g_offset + 4*kpt - 1)
write(file_id, *) 'tg(kpt) ', sfc
! Finally, soil temperature
soil = state(s_offset + msl*kpt : s_offset + 2*msl*kpt - 1)
write(file_id, *) 'tsoi(kpt,msl) ', soil

end subroutine write_sccm_state



subroutine read_sccm_state(state, time, file_id)
!------------------------------------------------------------------
!

real(r8), intent(out) :: state(:)
type(time_type), intent(out) :: time
integer, intent(in) :: file_id

real(r8) :: atmos(plev), sfc(kpt), soil(msl * kpt)
integer :: g_offset, s_offset, i, date_stamp
integer :: year, month, day, hour, minute, second, days, seconds
integer :: kpt_in, msl_in, plev_in, pcnst_in, delta_seconds
character * 40 :: temp_string, temp_string2

! Read the delta time entry which is not used
read(file_id, *) temp_string, delta_seconds

! First read day as date integer and seconds
read(file_id, *) temp_string, date_stamp
read(file_id, *) temp_string, seconds
year = date_stamp / 10000
month = (date_stamp - year * 10000) / 100
day = (date_stamp - year * 10000 - month * 100)
! Year needs a 19 or 20 in front
if(year > 50) then
   year = year + 1900
else
   year = year + 2000
endif

! First, set a date for just the day part with 0 seconds
time = set_date(year, month, day, 0, 0, 0)
! Then add in the seconds
time = time + set_time(seconds, 0)

! Read the model size descriptors, check for errors
read(file_id, *) temp_string, kpt_in
read(file_id, *) temp_string, msl_in
read(file_id, *) temp_string, plev_in
read(file_id, *) temp_string, pcnst_in
if(kpt_in /= kpt .or. msl_in /= msl .or. plev_in /= plev .or. pcnst_in /= pcnst)then
   write(*, *) 'ERROR ON MODELS SIZE FOR READ'
   stop
endif

! Read the atmospheric t
read(file_id, *) temp_string, atmos
state(2*plev + 1 : 3*plev) = atmos
write(*, *) 't is ', atmos
! Read atmospheric u
read(file_id, *) temp_string, atmos
state(1 : plev) = atmos
write(*, *) 'u is ', atmos
! Read atmospheric v
read(file_id, *) temp_string, atmos
state(plev + 1 : 2*plev) = atmos
write(*, *) 'v is ', atmos
! Read the tracer fields
do i = 1, pcnst
   read(file_id, *) temp_string, temp_string2, atmos
   state((2 + i)*plev + 1: (3 + i) * plev) = atmos
   write(*, *) 'tracer ', i, '  is ', atmos
end do
! Read the ps
read (file_id, *) temp_string, state((3 + pcnst)* plev + 1)
write(*, *) 'ps is ', state((3 + pcnst)*plev + 1)


! Offsets for where the ground and soil state variables start in dart
g_offset = (3 + pcnst) * plev + 2
s_offset = g_offset + 4 * kpt
! Snow moisture content is first
read(file_id, *) temp_string, sfc
state(g_offset : g_offset + kpt - 1) = sfc
! Canopy moisture next
read(file_id, *) temp_string, sfc
state(g_offset + kpt : g_offset + 2*kpt - 1) = sfc
! Now soil moisture
read(file_id, *) temp_string, temp_string2, soil
state(s_offset : s_offset + msl*kpt - 1) = soil
! Now temperature of vegetation
read(file_id, *) temp_string, sfc
state(g_offset + 2*kpt : g_offset + 3*kpt - 1) = sfc
! Now skin temperature
read(file_id, *) temp_string, sfc
state(g_offset + 3*kpt : g_offset + 4*kpt - 1) = sfc
! Finally, soil temperature
read(file_id, *) temp_string, temp_string2, soil
state(s_offset + msl*kpt : s_offset + 2*msl*kpt - 1) = soil

end subroutine read_sccm_state


!===================================================================
! End of model_mod
!===================================================================
end module model_mod
