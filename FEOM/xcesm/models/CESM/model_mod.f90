! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This is the interface between the CESM whole-system climate model and DART.

! Modules that are absolutely required for use are listed
use        types_mod, only : r8, missing_r8
use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                             print_time, print_date,                           &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=), operator(==)
use     location_mod, only : location_type, get_dist, get_close_maxdist_init,  &
                             get_close_obs_init, get_close_state_init,         &
                             set_location, get_location, get_close_type,       &
                             loc_get_close_obs => get_close_obs,               &
                             loc_get_close_state => get_close_state,           &
                             query_location
use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, to_upper, file_to_text, do_output,      &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, file_exist, find_textfile_dims,        &
                             do_nml_file, do_nml_term, nmlfileunit
use cross_comp_localization_mod, only : get_vert_alpha, get_xcomp_dist


use     obs_kind_mod     ! for now, include all


use pop_model_mod, only :                         &
  pop_get_model_size      => get_model_size,      &
  pop_adv_1step           => adv_1step,           &
  pop_get_state_meta_data => get_state_meta_data, &
  pop_model_interpolate   => model_interpolate,   &
  pop_get_model_time_step => get_model_time_step, &
  pop_static_init_model   => static_init_model,   &
  pop_end_model           => end_model,           &
  pop_init_time           => init_time,           &
  pop_init_conditions     => init_conditions,     &
  pop_nc_write_model_atts => nc_write_model_atts, &
  pop_nc_write_model_vars => nc_write_model_vars, &
  pop_pert_model_state    => pert_model_state,    &
  pop_ens_mean_for_model  => ens_mean_for_model,  &
  pop_get_close_obs       => get_close_obs,       &
  pop_get_close_state     => get_close_state,     &
  pop_restart_file_to_sv  => restart_file_to_sv,  &
  pop_sv_to_restart_file  => sv_to_restart_file,  &
  pop_get_gridsize        => get_gridsize,        &
  pop_model_convert_vert_obs    => model_convert_vert_obs,   &
  pop_model_convert_vert_state  => model_convert_vert_state, &
  pop_test_interpolation  => test_interpolation

use cam_model_mod, only :                            &
  cam_get_model_size       => get_model_size,        &
  cam_adv_1step            => adv_1step,             &
  cam_get_state_meta_data  => get_state_meta_data,   &
  cam_model_interpolate    => model_interpolate,     &
  cam_get_model_time_step  => get_model_time_step,   &
  cam_static_init_model    => static_init_model,     &
  cam_end_model            => end_model,             &
  cam_init_time            => init_time,             &
  cam_init_conditions      => init_conditions,       &
  cam_nc_write_model_atts  => nc_write_model_atts,   &
  cam_nc_write_model_vars  => nc_write_model_vars,   &
  cam_pert_model_state     => pert_model_state,      &
  cam_ens_mean_for_model   => ens_mean_for_model,    &
  cam_get_close_obs        => get_close_obs,         &
  cam_get_close_state      => get_close_state,       &
  cam_model_type           => model_type,            &
  cam_prog_var_to_vector   => prog_var_to_vector,    &
  cam_vector_to_prog_var   => vector_to_prog_var,    &
  cam_read_cam_init        => read_cam_init,         &
  cam_init_model_instance  => init_model_instance,   &
  cam_end_model_instance   => end_model_instance,    &
  cam_write_cam_init       => write_cam_init,        &
  cam_write_cam_times      => write_cam_times,       &
  cam_state_vector_to_model => state_vector_to_model, &
  cam_model_to_state_vector => model_to_state_vector, &
  cam_model_convert_vert_obs      => model_convert_vert_obs,  &
  cam_model_convert_vert_state    => model_convert_vert_state


use clm_model_mod, only :                             &
  clm_get_model_size       => get_model_size,         &
  clm_adv_1step            => adv_1step,              &
  clm_get_state_meta_data  => get_state_meta_data,    &
  clm_model_interpolate    => model_interpolate,      &
  clm_get_model_time_step  => get_model_time_step,    &
  clm_static_init_model    => static_init_model,      &
  clm_end_model            => end_model,              &
  clm_init_time            => init_time,              &
  clm_init_conditions      => init_conditions,        &
  clm_nc_write_model_atts  => nc_write_model_atts,    &
  clm_nc_write_model_vars  => nc_write_model_vars,    &
  clm_pert_model_state     => pert_model_state,       &
  clm_ens_mean_for_model   => ens_mean_for_model,     &
  clm_sv_to_restart_file   => sv_to_restart_file,     &
  clm_model_convert_vert_obs     => model_convert_vert_obs,   &
  clm_model_convert_vert_state   => model_convert_vert_state, &
  clm_to_dart_state_vector,  &
   get_gridsize,             &
   get_clm_restart_filename, &
   get_state_time,           &
   get_grid_vertval,         &
   compute_gridcell_value,   &
   DART_get_var,             &
   get_model_time



use typesizes
use netcdf 

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          static_init_model,      &
          end_model,              &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_state,       &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_state_init,   &
          get_close_obs,          &
          get_close_state,        &
          ens_mean_for_model,     &
          model_convert_vert_obs, &
          model_convert_vert_state, &
  restart_file_to_sv, &
  sv_to_restart_file, &
  get_cesm_restart_filename

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer :: cam_model_size, clm_model_size, pop_model_size

character(len=256) :: msgstring
logical, save :: module_initialized = .false.

! FIXME: make the number and names extensible, eventually.
! for the first implementation we only support CAM, POP, and CLM
! explicitly.

! Default to cam and pop
logical :: include_CAM = .true.
logical :: include_POP = .true.
logical :: include_CLM = .false.

integer  :: debug = 0   ! turn up for more and more debug messages


namelist /model_nml/  &
   include_CAM, &
   include_POP, &
   include_CLM, &
   debug

! FIXME: add the input and output filenames here for the
! restart files?

type(time_type) :: model_time, model_timestep
integer :: model_size    ! the state vector length


contains

!------------------------------------------------------------------
!------------------------------------------------------------------

subroutine static_init_model()

integer :: iunit, io, days, ss, dd, model_size

if ( module_initialized ) return ! only need to do this once.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized right away.
module_initialized = .true.

! Read the DART namelist for the CESM pseudo-model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

if (.not. include_CAM .and. .not. include_POP .and. .not. include_CLM) then
   write(msgstring,*)'at least one component must be selected in the &model_nml namelist'
   call error_handler(E_ERR,'static_init_model',msgstring,source,revision,revdate)
endif

if (include_CAM) call cam_static_init_model()
if (include_POP) call pop_static_init_model()
if (include_CLM) call clm_static_init_model()

model_timestep = get_model_time_step()
call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)
write(msgstring,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',msgstring,source,revision,revdate)

model_size = get_model_size()
write(msgstring,*)'model_size = ', model_size
call error_handler(E_MSG,'static_init_model',msgstring,source,revision,revdate)

end subroutine static_init_model

!------------------------------------------------------------

subroutine init_conditions(x)
 real(r8), intent(out) :: x(:)

! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.

character(len=128) :: msgstring2, msgstring3

msgstring2 = "cannot run perfect_model_obs with 'start_from_restart = .false.' "
msgstring3 = "use cesm_to_dart to generate an initial state"
call error_handler(E_ERR,'init_conditions', &
                  'CESM model interface has no built-in default state', &
                  source, revision, revdate, &
                  text2=msgstring2, text3=msgstring3)

! this code never reached - just here to avoid compiler warnings
! about an intent(out) variable not being set to a value.
x = 0.0_r8

end subroutine init_conditions

!------------------------------------------------------------------

subroutine adv_1step(x, time)
 real(r8),        intent(inout) :: x(:)
 type(time_type), intent(in)    :: time

! If the model could be called as a subroutine, does a single
! timestep advance.  CESM cannot be called this way, so fatal error
! if this routine is called.

character(len=128) :: msgstring2, msgstring3

msgstring2 = "DART cannot advance the CESM model; all observations must be within"
msgstring3 = "the current assimilation window.  Check the first/last obs times."
call error_handler(E_ERR,'adv_1step', &
                  'CESM model cannot be called as a subroutine.', &
                  source, revision, revdate, &
                  text2=msgstring2, text3=msgstring3)

end subroutine adv_1step

!------------------------------------------------------------------

function get_model_size()
 integer :: get_model_size

! Returns the size of the model as an integer. Required for all
! applications.

if ( .not. module_initialized ) call static_init_model

! These sizes are module globals - set them so later we can
! get the right offsets for different sections of the state vector.
cam_model_size = 0
pop_model_size = 0
clm_model_size = 0

if (include_CAM) cam_model_size = cam_get_model_size()
if (include_POP) pop_model_size = pop_get_model_size()
if (include_CLM) clm_model_size = clm_get_model_size()

get_model_size = cam_model_size + clm_model_size + pop_model_size

end function get_model_size

!------------------------------------------------------------------

subroutine init_time(time)
 type(time_type), intent(out) :: time

! Companion interface to init_conditions. Returns a time that is
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.

character(len=128) :: msgstring2, msgstring3

msgstring2 = "cannot run perfect_model_obs with 'start_from_restart = .false.' "
msgstring3 = 'use cesm_to_dart to generate an initial state which contains a timestamp'
call error_handler(E_ERR,'init_time', &
                  'CESM model interface has no built-in default time', &
                  source, revision, revdate, &
                  text2=msgstring2, text3=msgstring3)

! this code never reached - just here to avoid compiler warnings
! about an intent(out) variable not being set to a value.
time = set_time(0,0)

end subroutine init_time

!------------------------------------------------------------------

subroutine model_interpolate(x, location, obs_kind, interp_val, istatus)
 real(r8),            intent(in) :: x(:)
 type(location_type), intent(in) :: location
 integer,             intent(in) :: obs_kind
 real(r8),           intent(out) :: interp_val
 integer,            intent(out) :: istatus

! Model interpolate will interpolate any state variable to
! the given location given a state vector. The type of the variable being
! interpolated is obs_type since normally this is used to find the expected
! value of an observation at some location. The interpolated value is 
! returned in interp_val and istatus is 0 for success.

real(r8) :: llon, llat, lvert, loc_array(3)
integer  :: x_start, x_end
character(len=32) :: componentname

if ( .not. module_initialized ) call static_init_model


! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the 
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

interp_val = MISSING_R8     ! the DART bad value flag
istatus = 99                ! unknown error

! Get the individual locations values
loc_array = get_location(location)
llon  = loc_array(1)
llat  = loc_array(2)
lvert = loc_array(3)

if (debug > 1) print *, 'requesting interpolation of ', obs_kind, ' at ', llon, llat, lvert

! based on generic kind, decide which component should compute the forward operator.
! if we add support for more complicated forward operators that span components
! add code at a higher level in an obs_def mod and do the computation there.
call which_model_obs(obs_kind, componentname)
call set_start_end(componentname, x_start, x_end)

if (componentname == 'CAM') then
   call cam_model_interpolate(x(x_start:x_end), location, obs_kind, interp_val, istatus)
if (istatus > 0) print *, 'CAM obs_kind, istatus = ', obs_kind, istatus
else if (componentname == 'POP') then
   call pop_model_interpolate(x(x_start:x_end), location, obs_kind, interp_val, istatus)
if (istatus > 0) print *, 'POP obs_kind, istatus = ', obs_kind, istatus
else if (componentname == 'CLM') then
   call clm_model_interpolate(x(x_start:x_end), location, obs_kind, interp_val, istatus)
if (istatus > 0) print *, 'CLM obs_kind, istatus = ', obs_kind, istatus
else
   return
endif

if (debug > 1) print *, 'interp val, istatus = ', interp_val, istatus

end subroutine model_interpolate

!------------------------------------------------------------------

function get_model_time_step()
 type(time_type) :: get_model_time_step

! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: cam_time, pop_time, clm_time, base_time, zero_time

if ( .not. module_initialized ) call static_init_model

zero_time = set_time(0,0)

cam_time = zero_time
pop_time = zero_time
clm_time = zero_time

if (include_CAM) cam_time = cam_get_model_time_step()
if (include_POP) pop_time = pop_get_model_time_step()
if (include_CLM) clm_time = clm_get_model_time_step()

! make sure they are compatible here and error out if they are not.
if (cam_time > zero_time) then
   base_time = cam_time
else if (pop_time > zero_time) then
   base_time = pop_time
else if (clm_time > zero_time) then
   base_time = clm_time
else
   call error_handler(E_ERR, 'get_model_time_step', 'no models returned a good time step')
endif

if ((cam_time > zero_time .and. cam_time /= base_time) .or. &
    (pop_time > zero_time .and. pop_time /= base_time) .or. &
    (clm_time > zero_time .and. clm_time /= base_time)) then

   call print_time(cam_time, 'CAM time step')
   call print_time(cam_time, 'CAM time step', logfileunit)
   call print_time(pop_time, 'POP time step')
   call print_time(pop_time, 'POP time step', logfileunit)
   call print_time(clm_time, 'CLM time step')
   call print_time(clm_time, 'CLM time step', logfileunit)
   call error_handler(E_ERR, 'get_model_time_step', 'all model time steps must be the same', &
                      source, revision, revdate)
endif

if (base_time == zero_time) then
   call error_handler(E_ERR, 'get_model_time_step', 'all model time steps were zero', &
                      source, revision, revdate)
endif

get_model_time_step = base_time

end function get_model_time_step

!------------------------------------------------------------------

subroutine get_state_meta_data(index_in, location, var_type)
 integer,             intent(in)  :: index_in
 type(location_type), intent(out) :: location
 integer,             intent(out), optional :: var_type

! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument kind
! can be returned if the model has more than one type of field (for
! instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

real(r8) :: lat, lon, vert
integer  :: x_start, x_end
character(len=32) :: componentname

! figure out what offset in the state vector we are at, and then
! call the right sub-model

call which_model_state(index_in, componentname)
call set_start_end(componentname, x_start, x_end)

select case (componentname)
   case ('CAM')
      call cam_get_state_meta_data(index_in - x_start + 1, location, var_type)

   case ('POP')
      call pop_get_state_meta_data(index_in - x_start + 1, location, var_type)

   case ('CLM')
      call clm_get_state_meta_data(index_in - x_start + 1, location, var_type)

   case default
      ! this should not happen
      call error_handler(E_ERR, 'get_state_meta_data', 'error determining right model to use', &
                      source, revision, revdate)
end select


end subroutine get_state_meta_data

!------------------------------------------------------------------

subroutine end_model()

! Shutdown and clean-up.

if (include_CAM) call cam_end_model()
if (include_POP) call pop_end_model()
if (include_CLM) call clm_end_model()

end subroutine end_model

!------------------------------------------------------------------

function nc_write_model_atts(ncFileID, model)
 integer, intent(in)  :: ncFileID            ! netCDF file identifier
 character(len=*), optional, intent(in) :: model
 integer              :: nc_write_model_atts ! function return value

integer :: rc

if ( .not. module_initialized ) call static_init_model

if (present(model)) then
   if (model == 'cam' .and. include_CAM) then
      nc_write_model_atts = cam_nc_write_model_atts(ncFileID)
      return
   else if (model == 'pop' .and. include_POP) then
      nc_write_model_atts = pop_nc_write_model_atts(ncFileID)
      return
   else if (model == 'clm' .and. include_CLM) then
      nc_write_model_atts = clm_nc_write_model_atts(ncFileID)
      return
   endif
endif

nc_write_model_atts = -1 ! If we got here, bad things happened

end function nc_write_model_atts

!------------------------------------------------------------------

function nc_write_model_vars(ncFileID, statevec, copyindex, timeindex, model) 
 integer,                intent(in) :: ncFileID            ! netCDF file identifier
 real(r8), dimension(:), intent(in) :: statevec
 integer,                intent(in) :: copyindex
 integer,                intent(in) :: timeindex
 character(len=*), optional, intent(in) :: model
 integer                            :: nc_write_model_vars ! function return value

integer :: rc, x_start, x_end

if ( .not. module_initialized ) call static_init_model

if (present(model)) then
   if (model == 'cam' .and. include_CAM) then
      call set_start_end('CAM', x_start, x_end)
      nc_write_model_vars = cam_nc_write_model_vars(ncFileID, statevec(x_start:x_end), copyindex, timeindex) 
      return
   else if (model == 'pop' .and. include_POP) then
      call set_start_end('POP', x_start, x_end)
      nc_write_model_vars = pop_nc_write_model_vars(ncFileID, statevec(x_start:x_end), copyindex, timeindex) 
   else if (model == 'clm' .and. include_CLM) then
      call set_start_end('clm', x_start, x_end)
      nc_write_model_vars = clm_nc_write_model_vars(ncFileID, statevec(x_start:x_end), copyindex, timeindex) 
   endif

endif

nc_write_model_vars = -1 ! If we got here, bad things happened

end function nc_write_model_vars

!------------------------------------------------------------------

subroutine pert_model_state(state, pert_state, interf_provided)
 real(r8), intent(in)  :: state(:)
 real(r8), intent(out) :: pert_state(:)
 logical,  intent(out) :: interf_provided

! Perturbs a model state for generating initial ensembles.
! The perturbed state is returned in pert_state.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding a perturbation to each model state 
! variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.

integer :: x_start, x_end
logical :: had_interface, minv, maxv
integer :: tristate(3) 


if ( .not. module_initialized ) call static_init_model

! key to 'tristate' values:
!   0 = this model isn't being used
!   1 = this model has a perturb interface
!   2 = this model doesn't have a perturb interface
tristate(:) = 0

if (include_CAM) then
   call set_start_end('CAM', x_start, x_end)
   call cam_pert_model_state(state(x_start:x_end), pert_state(x_start:x_end), had_interface)
   if (had_interface) then
      tristate = 1
   else
      tristate = 2
   endif
endif

if (include_POP) then
   call set_start_end('POP', x_start, x_end)
   call pop_pert_model_state(state(x_start:x_end), pert_state(x_start:x_end), had_interface)
   if (had_interface) then
      tristate = 1
   else
      tristate = 2
   endif
endif

if (include_CLM) then
   call set_start_end('CLM', x_start, x_end)
   call clm_pert_model_state(state(x_start:x_end), pert_state(x_start:x_end), had_interface)
   if (had_interface) then
      tristate = 1
   else
      tristate = 2
   endif
endif

! FIXME: we cannot handle the case where some models want to
! perturb on their own and some want filter to do it.  it has to
! be every model or no models at this point.

minv = minval(tristate(minloc(tristate, tristate /= 0)))
maxv = maxval(tristate(maxloc(tristate, tristate /= 0)))

if (minv /= maxv) then
   call error_handler(E_MSG, 'pert_model_state', &
      'if any models use a perturb routine, all models must use a perturb routine', &
      source, revision, revdate, text2='overriding local pert routines and using filter code')
   interf_provided = .false.
else
   if (minv == 1) then
      interf_provided = .true.
   else if (minv == 2) then
      interf_provided = .false.
   endif
endif


end subroutine pert_model_state

!------------------------------------------------------------------

subroutine ens_mean_for_model(ens_mean)
 real(r8), intent(in) :: ens_mean(:)

! If needed by the model interface, this is the current mean
! for all state vector items across all ensembles. It is up to this
! code to allocate space and save a copy if it is going to be used
! later on.  For now, we are ignoring it.

integer :: x_start, x_end

if ( .not. module_initialized ) call static_init_model

if (include_CAM) then
   call set_start_end('CAM', x_start, x_end)
   call cam_ens_mean_for_model(ens_mean(x_start:x_end))
endif

if (include_POP) then
   call set_start_end('POP', x_start, x_end)
   call pop_ens_mean_for_model(ens_mean(x_start:x_end))
endif

if (include_CLM) then
   call set_start_end('CLM', x_start, x_end)
   call clm_ens_mean_for_model(ens_mean(x_start:x_end))
endif

end subroutine ens_mean_for_model

!------------------------------------------------------------------

subroutine model_convert_vert_obs(item_count, loc_list, type_list, vertical_localization_coordinate)

integer,             intent(in)    :: item_count
type(location_type), intent(inout) :: loc_list(item_count)
integer,             intent(in)    :: type_list(item_count)
integer,             intent(in)    :: vertical_localization_coordinate

integer :: i, sublist_count
type(location_type), allocatable :: loc_sublist(:)
integer, allocatable :: type_sublist(:)

if ( .not. module_initialized ) call static_init_model

if (include_CAM) then
   ! tricky - here we have to extract the obs types for atm only
   sublist_count = 0
   do i=1, item_count
      if (type_list(i) == RADIOSONDE_TEMPERATURE) sublist_count = sublist_count + 1
   enddo
   if (sublist_count > 0) then
      allocate(loc_sublist(sublist_count), type_sublist(sublist_count))
      sublist_count = 0
      do i=1, item_count
         if (type_list(i) == RADIOSONDE_TEMPERATURE) then
             sublist_count = sublist_count + 1
             loc_sublist(sublist_count)  = loc_list(i)
             type_sublist(sublist_count) = type_list(i)
         endif
      enddo
 
      call cam_model_convert_vert_obs(sublist_count, loc_sublist, type_sublist, vertical_localization_coordinate)

      do i=1, item_count
         if (type_list(i) == RADIOSONDE_TEMPERATURE) then
             sublist_count = sublist_count + 1
             loc_list(i) = loc_sublist(sublist_count) 
         endif
      enddo
      deallocate(loc_sublist, type_sublist)
   endif
endif

if (include_POP) then
   ! tricky - here we have to extract the obs types for ocn only
   sublist_count = 0
   do i=1, item_count
      if (type_list(i) == XBT_TEMPERATURE) sublist_count = sublist_count + 1
   enddo
   if (sublist_count > 0) then
      allocate(loc_sublist(sublist_count), type_sublist(sublist_count))
      sublist_count = 0
      do i=1, item_count
         if (type_list(i) == XBT_TEMPERATURE) then
             sublist_count = sublist_count + 1
             loc_sublist(sublist_count)  = loc_list(i)
             type_sublist(sublist_count) = type_list(i)
         endif
      enddo
 
      call pop_model_convert_vert_obs(sublist_count, loc_sublist, type_sublist, vertical_localization_coordinate)

      do i=1, item_count
         if (type_list(i) == XBT_TEMPERATURE) then
             sublist_count = sublist_count + 1
             loc_list(i) = loc_sublist(sublist_count) 
         endif
      enddo
      deallocate(loc_sublist, type_sublist)
   endif
endif

if (include_CLM) then
   call clm_model_convert_vert_obs(item_count, loc_list, type_list, vertical_localization_coordinate)
endif

end subroutine model_convert_vert_obs

!------------------------------------------------------------------

subroutine model_convert_vert_state(item_count, loc_list, kind_list, indx_list, vertical_localization_coordinate)

integer,             intent(in)    :: item_count
type(location_type), intent(inout) :: loc_list(item_count)
integer,             intent(in)    :: kind_list(item_count)
integer,             intent(in)    :: indx_list(item_count)
integer,             intent(in)    :: vertical_localization_coordinate

integer :: x_start, x_end
integer :: i, sublist_count
type(location_type), allocatable :: loc_sublist(:)
integer, allocatable :: kind_sublist(:), indx_sublist(:)

if ( .not. module_initialized ) call static_init_model

if (include_CAM) then
   ! tricky - here we have to extract only the locs for CAM
   call set_start_end('CAM', x_start, x_end)
   ! FIXME: put this in a subr since we're going to do it 3 times
   sublist_count = 0
   do i=1, item_count
      if (indx_list(i) >= x_start .and. indx_list(i) <= x_end) sublist_count = sublist_count + 1
   enddo
   if (sublist_count > 0) then
      allocate(loc_sublist(sublist_count), kind_sublist(sublist_count), indx_sublist(sublist_count))
      sublist_count = 0
      do i=1, item_count
         if (indx_list(i) >= x_start .and. indx_list(i) <= x_end) then
             sublist_count = sublist_count + 1
             loc_sublist(sublist_count)  = loc_list(i)
             kind_sublist(sublist_count) = kind_list(i)
             indx_sublist(sublist_count) = indx_list(i)
         endif
      enddo
 
      call cam_model_convert_vert_state(sublist_count, loc_sublist, kind_sublist, indx_sublist, vertical_localization_coordinate)
 
      do i=1, item_count
         if (indx_list(i) >= x_start .and. indx_list(i) <= x_end) then
             sublist_count = sublist_count + 1
             loc_list(i) = loc_sublist(sublist_count) 
         endif
      enddo
      deallocate(loc_sublist, kind_sublist, indx_sublist)
   endif
endif

if (include_POP) then
   ! tricky - here we have to extract only the locs for POP
   call set_start_end('POP', x_start, x_end)
   ! FIXME: put this in a subr since we're going to do it 3 times
   sublist_count = 0
   do i=1, item_count
      if (indx_list(i) >= x_start .and. indx_list(i) <= x_end) sublist_count = sublist_count + 1
   enddo
   if (sublist_count > 0) then
      allocate(loc_sublist(sublist_count), kind_sublist(sublist_count), indx_sublist(sublist_count))
      sublist_count = 0
      do i=1, item_count
         if (indx_list(i) >= x_start .and. indx_list(i) <= x_end) then
             sublist_count = sublist_count + 1
             loc_sublist(sublist_count)  = loc_list(i)
             kind_sublist(sublist_count) = kind_list(i)
             indx_sublist(sublist_count) = indx_list(i)
         endif
      enddo
 
      call pop_model_convert_vert_state(sublist_count, loc_sublist, kind_sublist, indx_sublist, vertical_localization_coordinate)

      do i=1, item_count
         if (indx_list(i) >= x_start .and. indx_list(i) <= x_end) then
             sublist_count = sublist_count + 1
             loc_list(i) = loc_sublist(sublist_count) 
         endif
      enddo
      deallocate(loc_sublist, kind_sublist, indx_sublist)
   endif
endif

if (include_CLM) then
   call set_start_end('CLM', x_start, x_end)
   ! FIXME: ditto
   call clm_model_convert_vert_state(item_count, loc_list, kind_list, indx_list, vertical_localization_coordinate)
endif

end subroutine model_convert_vert_state

!------------------------------------------------------------------

subroutine restart_file_to_sv(filenames, state_vector, model_time)
 character(len=*), intent(in)    :: filenames(:)
 real(r8),         intent(inout) :: state_vector(:)
 type(time_type),  intent(out)   :: model_time

! FIXME: we can figure out the parts of the state vector, but
! the filenames are going to be separate so what do we do?

integer :: x_start, x_end, used

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8
used = 1

if (include_CAM) then
   call set_start_end('CAM', x_start, x_end)
   call cam_model_to_state_vector(filenames(used), state_vector(x_start:x_end), model_time)
   used = used + 1
endif

if (include_POP) then
   call set_start_end('POP', x_start, x_end)
   call pop_restart_file_to_sv(filenames(used), state_vector(x_start:x_end), model_time)
   used = used + 1
endif

if (include_CLM) then
   call set_start_end('CLM', x_start, x_end)
   call clm_to_dart_state_vector(state_vector(x_start:x_end), model_time)
   used = used + 1
endif

end subroutine restart_file_to_sv

!------------------------------------------------------------------

subroutine get_cesm_restart_filename(filenames)
 character(len=*), intent(out) :: filenames(:)

! FIXME: this is a hack and should be in the namelist somewhere!
integer :: used

used = 1
if (include_CAM) then
   filenames(used) = 'caminput.nc'
   used = used + 1
endif

if (include_POP) then
   filenames(used) = 'pop.r.nc'
   used = used + 1
endif

if (include_CLM) then
   filenames(used) = 'clminput.nc'
   used = used + 1
endif

end subroutine get_cesm_restart_filename

!------------------------------------------------------------------

subroutine sv_to_restart_file(state_vector, filenames, state_time)
 real(r8),         intent(in) :: state_vector(:)
 character(len=*), intent(in) :: filenames(:)
 type(time_type),  intent(in) :: state_time

integer :: x_start, x_end, used

if ( .not. module_initialized ) call static_init_model

used = 1
if (include_CAM) then
   call set_start_end('CAM', x_start, x_end)
   call cam_state_vector_to_model(state_vector(x_start:x_end), filenames(used), state_time)
   used = used + 1
endif

if (include_POP) then
   call set_start_end('POP', x_start, x_end)
   call pop_sv_to_restart_file(state_vector(x_start:x_end), filenames(used), state_time)
   used = used + 1
endif

if (include_CLM) then
   call set_start_end('CLM', x_start, x_end)
   call clm_sv_to_restart_file(state_vector(x_start:x_end), filenames(used), state_time)
   used = used + 1
endif

end subroutine sv_to_restart_file

!------------------------------------------------------------------
!------------------------------------------------------------------


subroutine get_close_obs(gc, base_obs_loc, base_obs_type, &
                         locs, loc_kind, num_close, close_ind, dist)

 type(get_close_type), intent(in)    :: gc
 type(location_type),  intent(in)    :: base_obs_loc
 integer,              intent(in)    :: base_obs_type
 type(location_type),  intent(inout) :: locs(:)
 integer,              intent(in)    :: loc_kind(:)
 integer,              intent(out)   :: num_close
 integer,              intent(out)   :: close_ind(:)
 real(r8),  optional,  intent(out)   :: dist(:)

! Given a DART location (referred to as "base") and a set of candidate
! locations & kinds (obs, obs_kind), returns the subset close to the
! "base", their indices, and their distances to the "base" ...

! For vertical distance computations, general philosophy is to convert all
! vertical coordinates to a common coordinate. This coordinate type is defined
! in the namelist with the variable "vert_localization_coord".

character(len=32) :: componentname
integer :: i, this, base_obs_kind
real(r8) :: horiz_dist, vert_dist_proxy
integer :: vert1, vert2
integer, save :: n = 1

! Initialize variables to missing status
num_close = 0
close_ind(:) = -1
if (present(dist)) dist = 1.0e9   !something big and positive (far away)

! FIXME: in a real unified model_mod, these would all be called
! and any state vector items from any model are potentially close.
! the vertical conversions are one issue; the other is whether the
! distances should be the same or different in different mediums
! (e.g air vs water vs soil/snow)

! this needs to call all the get close routines for active components
! and merge the lists when done.

! if we go back to generic kinds, then which_model_obs() needs to
! take a kind and this is a real specific type

! we have to allow vertical conversions per component - and then
! what, sum them?  this needs help.

! @TODO:  for now, bypass the model-specific routines and call
! the raw location mod code.   this may mess up land points in pop
! and won't allow top-of-model changes in cam, and no vertical conversion
! will be done.  but it may run with horizontal only.

! NO DISTANCES YET!
call loc_get_close_obs(gc, base_obs_loc, base_obs_type, &
                       locs, loc_kind, num_close, close_ind)

if (.not. present(dist)) return

! get these once before looping below
base_obs_kind = get_obs_kind_var_type(base_obs_type)
vert1 = query_location(base_obs_loc)



do i=1, num_close
   ! compute the horizontal here and then if cross component do something
   ! different in the vertical.
   this = close_ind(i)
   vert2 = query_location(locs(this))
!if (mod(n, 1000) == 1 .and. vert1 /= vert2) print *, 'vert1, 2: ', vert1, vert2
n = n + 1
   if ((base_obs_kind == KIND_AIR_TEMPERATURE   .and. loc_kind(this) == KIND_WATER_TEMPERATURE) .or. &
       (base_obs_kind == KIND_WATER_TEMPERATURE .and. loc_kind(this) == KIND_AIR_TEMPERATURE)) then
      dist(i) = get_xcomp_dist(base_obs_loc, locs(this), base_obs_kind, loc_kind(this))
if (dist(i) == 0) then
print *, '0 distance: ', base_obs_kind, loc_kind(this), vert1, vert2
endif
      ! or
      !horiz_dist = get_dist(base_obs_loc, locs(this), &
      !                      base_obs_type, loc_kind(this), no_vert = .true.)
      !vert_dist_proxy = get_vert_alpha(base_loc_obs, loc(this), base_obs_kind, loc_kind(this))
      !dist(i) = sqrt(horiz_dist**2 + vert_dist_proxy**2)
   else
      vert2 = query_location(locs(this))
      ! FIXME: we should be converting in the vertical for those we know how to convert
      ! ( e.g. height obs in the atmosphere -> pressure)
      if (vert1 /= vert2) then
         dist(i) = get_dist(base_obs_loc, locs(this), base_obs_type, loc_kind(this), no_vert=.true.)
      else
         dist(i) = get_dist(base_obs_loc, locs(this), base_obs_type, loc_kind(this))
      endif
   endif
enddo

return

! NOT REACHED NOT REACHED NOT REACHED

!  FIXME:  @TODO
if (include_CAM) then
      call cam_get_close_obs(gc, base_obs_loc, base_obs_type, &
                             locs, loc_kind, num_close, close_ind, dist)
endif
! save original close_ind list?
if (include_POP) then
      call pop_get_close_obs(gc, base_obs_loc, base_obs_type, &
                             locs, loc_kind, num_close, close_ind, dist)
endif
! merge lists & dists together
if (include_CLM) then
      call loc_get_close_obs(gc, base_obs_loc, base_obs_type, &
                             locs, loc_kind, num_close, close_ind, dist)
endif
! merge lists & dists together

 
end subroutine get_close_obs

!------------------------------------------------------------------

subroutine get_close_state(gc, base_obs_loc, base_obs_type, &
                          locs, loc_kind, num_close, close_ind, dist)

 type(get_close_type), intent(in)    :: gc
 type(location_type),  intent(in)    :: base_obs_loc
 integer,              intent(in)    :: base_obs_type
 type(location_type),  intent(inout) :: locs(:)
 integer,              intent(in)    :: loc_kind(:)
 integer,              intent(out)   :: num_close
 integer,              intent(out)   :: close_ind(:)
 real(r8),  optional,  intent(out)   :: dist(:)

! Given a DART location (referred to as "base") and a set of candidate
! locations & kinds (obs, obs_kind), returns the subset close to the
! "base", their indices, and their distances to the "base" ...

! For vertical distance computations, general philosophy is to convert all
! vertical coordinates to a common coordinate. This coordinate type is defined
! in the namelist with the variable "vert_localization_coord".

integer :: whole_list, cur_start, cur_end
character(len=32) :: componentname

integer,  allocatable :: unified_close_ind(:)
real(r8), allocatable :: unified_dist(:)

! Initialize variables to missing status
whole_list = size(close_ind)
allocate(unified_close_ind(whole_list), unified_dist(whole_list))
num_close = 0
close_ind(:) = -1
unified_close_ind = -1
if (present(dist)) dist = 1.0e9   !something big and positive (far away)
unified_dist = 1.0e9

! FIXME: in a real unified model_mod, these would all be called
! and any state vector items from any model are potentially close.
! the vertical conversions are one issue; the other is whether the
! distances should be the same or different in different mediums
! (e.g air vs water vs soil/snow)

! this needs to call all the get close routines for active components
! and merge the lists when done.

! if we go back to generic kinds, then which_model_obs() needs to
! take a kind and this is a real specific type

! @TODO - need to sort and uniq the list in case models agree on
! the same locations.  and what if they give different distances?
 
! @TODO:  for now, bypass the model-specific routines and call
! the raw location mod code.   this may mess up land points in pop
! and won't allow top-of-model changes in cam, and no vertical conversion
! will be done.  but it may run with horizontal only.

call get_close_obs(gc, base_obs_loc, base_obs_type, &
                   locs, loc_kind, num_close, close_ind, dist)

deallocate(unified_close_ind, unified_dist)
return

! NOT REACHED NOT REACHED NOT REACHED

call loc_get_close_obs(gc, base_obs_loc, base_obs_type, &
                       locs, loc_kind, num_close, close_ind, dist)

deallocate(unified_close_ind, unified_dist)
return


cur_start = 1
cur_end = 1

!  concatinate each list as it is computed
if (include_CAM) then
      call cam_get_close_state(gc, base_obs_loc, base_obs_type, &
                               locs, loc_kind, num_close, close_ind, dist)
      cur_end = cur_start + num_close - 1
      unified_close_ind(cur_start:cur_end) = close_ind(1:num_close)       
      if (present(dist)) &
         unified_dist(cur_start:cur_end) = dist(1:num_close)
      cur_start = cur_start + num_close
endif

if (include_POP) then
      call pop_get_close_state(gc, base_obs_loc, base_obs_type, &
                               locs, loc_kind, num_close, close_ind, dist)
      cur_end = cur_start + num_close - 1
      unified_close_ind(cur_start:cur_end) = close_ind(1:num_close)       
      if (present(dist)) &
         unified_dist(cur_start:cur_end) = dist(1:num_close)
      cur_start = cur_start + num_close
endif

if (include_CLM) then
      call loc_get_close_state(gc, base_obs_loc, base_obs_type, &
                               locs, loc_kind, num_close, close_ind, dist)
      cur_end = cur_start + num_close - 1
      unified_close_ind(cur_start:cur_end) = close_ind(1:num_close)       
      if (present(dist)) &
         unified_dist(cur_start:cur_end) = dist(1:num_close)
      cur_start = cur_start + num_close
endif

! and return combined lists and counts
num_close = cur_start
close_ind = unified_close_ind
if (present(dist)) dist = unified_dist

deallocate(unified_close_ind, unified_dist)
 
end subroutine get_close_state

!------------------------------------------------------------------
!------------------------------------------------------------------
! additional worker routines which figure out which of the other
! components should be called.
!------------------------------------------------------------------

subroutine set_start_end(componentname, x_start, x_end)
 character(len=*), intent(in)  :: componentname
 integer,          intent(out) :: x_start, x_end

! @TODO this assumes a fixed order
select case (componentname)
   case ('CAM')
      x_start = 1
      x_end = x_start + cam_model_size - 1

   case ('POP') 
      x_start = cam_model_size + 1
      x_end = x_start + pop_model_size - 1

   case ('CLM')
      x_start = cam_model_size + pop_model_size + 1
      x_end = x_start + clm_model_size - 1

   case default
      x_start = -1
      x_end = -1

end select

end subroutine set_start_end

!------------------------------------------------------------------

subroutine which_model_state(x_offset, componentname)
 integer,          intent(in)  :: x_offset
 character(len=*), intent(out) :: componentname

integer :: x_start, x_end

call set_start_end('CAM', x_start, x_end)
if (x_offset >= x_start .and. x_offset <= x_end) then
   componentname = 'CAM'
   return
endif
 
call set_start_end('POP', x_start, x_end)
if (x_offset >= x_start .and. x_offset <= x_end) then
   componentname = 'POP'
   return
endif
 
call set_start_end('CLM', x_start, x_end)
if (x_offset >= x_start .and. x_offset <= x_end) then
   componentname = 'CLM'
   return
endif
 
! unknown
componentname = 'NULL'

end subroutine which_model_state

!------------------------------------------------------------------

subroutine which_model_obs(obs_kind, componentname)
 integer,          intent(in)  :: obs_kind
 character(len=*), intent(out) :: componentname

! FIXME: this needs to be beefed up with all possible obs kinds
! and which model would have the best forward operator for it.

select case (obs_kind)
   case (KIND_AIR_TEMPERATURE)
      componentname = 'CAM'
   case (KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT)
      componentname = 'CAM'
 
   case (KIND_WATER_TEMPERATURE)
      componentname = 'POP'
   case (KIND_U_CURRENT_COMPONENT, KIND_V_CURRENT_COMPONENT)
      componentname = 'POP'

   case (KIND_CARBON)
      componentname = 'CLM'

   case default
      ! unknown
      componentname = 'NULL'

end select

if (.not. valid_component(componentname)) componentname = 'NULL'

end subroutine which_model_obs

!------------------------------------------------------------------

function valid_component(componentname)
 character(len=*), intent(in)  :: componentname
 logical :: valid_component

! assume failure
valid_component = .false.

select case (componentname)
   case ('CAM')
      if (include_CAM) valid_component = .TRUE.
   case ('POP')
      if (include_POP) valid_component = .TRUE.
   case ('CLM')
      if (include_CLM) valid_component = .TRUE.
end select

end function valid_component

!------------------------------------------------------------------
! End of model_mod
!------------------------------------------------------------------

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
