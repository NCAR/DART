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
                             operator(/=), operator(<=)
use     location_mod, only : location_type, get_dist, get_close_maxdist_init,  &
                             get_close_obs_init, set_location,                 &
                             get_location, loc_get_close_obs => get_close_obs, &
                             get_close_type
use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper,                    &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, file_exist, find_textfile_dims,        &
                             file_to_text, do_output
use     obs_kind_mod     ! for now, include all


use pop_model_mod
use clm_model_mod
use cam_model_mod


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
          get_close_obs,          &
          ens_mean_for_model, &
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

! FIXME: for now make cam the only default
logical :: include_CAM = .true.
logical :: include_POP = .false.
logical :: include_CLM = .false.

integer  :: debug = 0   ! turn up for more and more debug messages


namelist /model_nml/  &
   include_CAM, &
   include_POP, &
   include_CLM, &
   debug

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
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
if (do_output()) write(logfileunit, nml=model_nml)
if (do_output()) write(     *     , nml=model_nml)


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
msgstring3 = 'use cesm_to_dart to generate an initial state'
call error_handler(E_ERR,'init_conditions', &
                  'WARNING!!  CESM model has no built-in default state', &
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

call error_handler(E_ERR,'adv_1step', &
                  'CESM model cannot be called as a subroutine; async cannot = 0', &
                  source, revision, revdate)

end subroutine adv_1step

!------------------------------------------------------------------

function get_model_size()
 integer :: get_model_size

! Returns the size of the model as an integer. Required for all
! applications.

if ( .not. module_initialized ) call static_init_model

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
                  'WARNING!!  CESM model has no built-in default time', &
                  source, revision, revdate, &
                  text2=msgstring2, text3=msgstring3)

! this code never reached - just here to avoid compiler warnings
! about an intent(out) variable not being set to a value.
time = set_time(0,0)

end subroutine init_time

!------------------------------------------------------------------

subroutine model_interpolate(x, location, obs_type, interp_val, istatus)
 real(r8),            intent(in) :: x(:)
 type(location_type), intent(in) :: location
 integer,             intent(in) :: obs_type
 real(r8),           intent(out) :: interp_val
 integer,            intent(out) :: istatus

! Model interpolate will interpolate any state variable to
! the given location given a state vector. The type of the variable being
! interpolated is obs_type since normally this is used to find the expected
! value of an observation at some location. The interpolated value is 
! returned in interp_val and istatus is 0 for success.

real(r8) :: llon, llat, lvert, loc_array(3)
integer  :: x_start, x_end
character(len=32) :: modelname

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

if (debug > 1) print *, 'requesting interpolation of ', obs_type, ' at ', llon, llat, lvert

call which_model_obs(obs_type, modelname)
call set_start_end(modelname, x_start, x_end)

if (modelname == 'CAM') then
   call cam_model_interpolate(x(x_start:x_end), location, obs_type, interp_val, istatus)
else if (modelname == 'POP') then
   call pop_model_interpolate(x(x_start:x_end), location, obs_type, interp_val, istatus)
else if (modelname == 'CLM') then
   call clm_model_interpolate(x(x_start:x_end), location, obs_type, interp_val, istatus)
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

type(time_type) :: cam_time, clm_time, pop_time

if ( .not. module_initialized ) call static_init_model

if (include_CAM) cam_time = cam_get_model_time_step()
if (include_POP) pop_time = pop_get_model_time_step()
if (include_CLM) clm_time = clm_get_model_time_step()

! FIXME:
! make sure they are compatible here

get_model_time_step = cam_time

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
character(len=32) :: modelname

! figure out what offset in the state vector we are at, and then
! call the right sub-model

call which_model_state(index_in, modelname)
call set_start_end(modelname, x_start, x_end)

if (modelname == 'CAM') then
   call cam_get_state_meta_data(index_in - x_start, location, var_type)

else if (modelname == 'POP') then
   call pop_get_state_meta_data(index_in - x_start, location, var_type)

else if (modelname == 'CLM') then
   call clm_get_state_meta_data(index_in - x_start, location, var_type)

else
   call error_handler(E_ERR, 'get_state_meta_data', 'offset beyond state vector length', &
                      source, revision, revdate)
endif


end subroutine get_state_meta_data

!------------------------------------------------------------------

subroutine end_model()

! Shutdown and clean-up.

if (include_CAM) call cam_end_model()
if (include_POP) call pop_end_model()
if (include_CLM) call clm_end_model()

end subroutine end_model

!------------------------------------------------------------------

function nc_write_model_atts(ncFileID)
 integer, intent(in)  :: ncFileID            ! netCDF file identifier
 integer              :: nc_write_model_atts ! function return value

integer :: rc

if ( .not. module_initialized ) call static_init_model

if (include_CAM) then
   rc = cam_nc_write_model_atts(ncFileID)
   if (rc /= 0) then
      nc_write_model_atts = rc
      return
   endif
endif

if (include_POP) then
   rc = pop_nc_write_model_atts(ncFileID)
   if (rc /= 0) then
      nc_write_model_atts = rc
      return
   endif
endif

if (include_CLM) then
   rc = clm_nc_write_model_atts(ncFileID)
   if (rc /= 0) then
      nc_write_model_atts = rc
      return
   endif
endif

nc_write_model_atts = 0 ! If we got here, things went well.

end function nc_write_model_atts

!------------------------------------------------------------------

function nc_write_model_vars(ncFileID, statevec, copyindex, timeindex) 
 integer,                intent(in) :: ncFileID            ! netCDF file identifier
 real(r8), dimension(:), intent(in) :: statevec
 integer,                intent(in) :: copyindex
 integer,                intent(in) :: timeindex
 integer                            :: nc_write_model_vars ! function return value

integer :: rc, x_start, x_end

if ( .not. module_initialized ) call static_init_model

if (include_CAM) then
   call set_start_end('CAM', x_start, x_end)
   rc = cam_nc_write_model_vars(ncFileID, statevec(x_start:x_end), copyindex, timeindex) 
   if (rc /= 0) then
      nc_write_model_vars = rc
      return
   endif
endif

if (include_POP) then
   call set_start_end('POP', x_start, x_end)
   rc = pop_nc_write_model_vars(ncFileID, statevec(x_start:x_end), copyindex, timeindex) 
   if (rc /= 0) then
      nc_write_model_vars = rc
      return
   endif
endif

if (include_CLM) then
   call set_start_end('CLM', x_start, x_end)
   rc = clm_nc_write_model_vars(ncFileID, statevec(x_start:x_end), copyindex, timeindex) 
   if (rc /= 0) then
      nc_write_model_vars = rc
      return
   endif
endif

nc_write_model_vars = 0 ! If we got here, things went well.

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

if ( .not. module_initialized ) call static_init_model

if (include_CAM) then
   call set_start_end('CAM', x_start, x_end)
   call cam_pert_model_state(state(x_start:x_end), pert_state(x_start:x_end), interf_provided)
endif

if (include_POP) then
   call set_start_end('POP', x_start, x_end)
   call pop_pert_model_state(state(x_start:x_end), pert_state(x_start:x_end), interf_provided)
endif

if (include_CLM) then
   call set_start_end('CLM', x_start, x_end)
   call clm_pert_model_state(state(x_start:x_end), pert_state(x_start:x_end), interf_provided)
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

subroutine restart_file_to_sv(filename, state_vector, model_time)
 character(len=*), intent(in)    :: filename 
 real(r8),         intent(inout) :: state_vector(:)
 type(time_type),  intent(out)   :: model_time

! FIXME: we can figure out the parts of the state vector, but
! the filenames are going to be separate so what do we do?

integer :: x_start, x_end

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

if (include_CAM) then
   call set_start_end('CAM', x_start, x_end)
   call cam_restart_file_to_sv(filename, state_vector(x_start:x_end), model_time)
endif

if (include_POP) then
   call set_start_end('POP', x_start, x_end)
   call pop_restart_file_to_sv(filename, state_vector(x_start:x_end), model_time)
endif

if (include_CLM) then
   call set_start_end('CLM', x_start, x_end)
   call clm_restart_file_to_sv(filename, state_vector(x_start:x_end), model_time)
endif

end subroutine restart_file_to_sv

!------------------------------------------------------------------

subroutine get_cesm_restart_filename(filename)
 character(len=*), intent(out) :: filename

! FIXME:
filename = 'dummy'

end subroutine get_cesm_restart_filename

!------------------------------------------------------------------

subroutine sv_to_restart_file(state_vector, filename, statedate)
 real(r8),         intent(in) :: state_vector(:)
 character(len=*), intent(in) :: filename 
 type(time_type),  intent(in) :: statedate

integer :: x_start, x_end

if ( .not. module_initialized ) call static_init_model

call set_start_end('CAM', x_start, x_end)
call cam_sv_to_restart_file(state_vector(x_start:x_end), filename, statedate)

call set_start_end('POP', x_start, x_end)
call pop_sv_to_restart_file(state_vector(x_start:x_end), filename, statedate)

call set_start_end('CLM', x_start, x_end)
call clm_sv_to_restart_file(state_vector(x_start:x_end), filename, statedate)

end subroutine sv_to_restart_file

!------------------------------------------------------------------
!------------------------------------------------------------------

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, &
                         obs, obs_kind, num_close, close_ind, dist)

 type(get_close_type),              intent(in) :: gc
 type(location_type),               intent(in) :: base_obs_loc
 integer,                           intent(in) :: base_obs_kind
 type(location_type), dimension(:), intent(in) :: obs
 integer,             dimension(:), intent(in) :: obs_kind
 integer,                           intent(out):: num_close
 integer,             dimension(:), intent(out):: close_ind
 real(r8),  optional, dimension(:), intent(out):: dist

! Given a DART location (referred to as "base") and a set of candidate
! locations & kinds (obs, obs_kind), returns the subset close to the
! "base", their indices, and their distances to the "base" ...

! For vertical distance computations, general philosophy is to convert all
! vertical coordinates to a common coordinate. This coordinate type is defined
! in the namelist with the variable "vert_localization_coord".

integer :: t_ind, k
character(len=32) :: modelname

! Initialize variables to missing status
num_close = 0
close_ind(:) = -1
if (present(dist)) dist = 1.0e9   !something big and positive (far away)

call which_model_obs(base_obs_kind, modelname)
if (modelname == 'CAM') then
   call cam_get_close_obs(gc, base_obs_loc, base_obs_kind, &
                          obs, obs_kind, num_close, close_ind, dist)

else if (modelname == 'POP') then
   call pop_get_close_obs(gc, base_obs_loc, base_obs_kind, &
                          obs, obs_kind, num_close, close_ind, dist)

else if (modelname == 'CLM') then
   call clm_get_close_obs(gc, base_obs_loc, base_obs_kind, &
                          obs, obs_kind, num_close, close_ind, dist)

endif


end subroutine get_close_obs

!------------------------------------------------------------------
!------------------------------------------------------------------
! additional worker routines which figure out which of the other
! components should be called.
!------------------------------------------------------------------

subroutine set_start_end(modelname, x_start, x_end)
 character(len=*), intent(in)  :: modelname
 integer,          intent(out) :: x_start, x_end

select case (modelname)
   case ('CAM')
      x_start = 1
      x_end = cam_model_size

   case ('POP') 
      x_start = cam_model_size + 1
      x_end = x_start + pop_model_size

   case ('CLM')
      x_start = cam_model_size + pop_model_size + 1
      x_end = x_start + clm_model_size 

   case default
      x_start = -1
      x_end = -1

end select

end subroutine set_start_end

!------------------------------------------------------------------

subroutine which_model_state(x_offset, modelname)
 integer,          intent(in)  :: x_offset
 character(len=*), intent(out) :: modelname

integer :: x_start, x_end

call set_start_end('CAM', x_start, x_end)
if (x_offset >= x_start .and. x_offset <= x_end) then
   modelname = 'CAM'
   return
endif
 
call set_start_end('POP', x_start, x_end)
if (x_offset >= x_start .and. x_offset <= x_end) then
   modelname = 'POP'
   return
endif
 
call set_start_end('CLM', x_start, x_end)
if (x_offset >= x_start .and. x_offset <= x_end) then
   modelname = 'CLM'
   return
endif
 
! unknown
modelname = 'NULL'

end subroutine which_model_state

!------------------------------------------------------------------

subroutine which_model_obs(obs_type, modelname)
 integer,          intent(in)  :: obs_type
 character(len=*), intent(out) :: modelname


select case (obs_type)
   case (RADIOSONDE_TEMPERATURE)
      modelname = 'CAM'
   case (RADIOSONDE_U_WIND_COMPONENT, RADIOSONDE_V_WIND_COMPONENT)
      modelname = 'CAM'
 
   case (XBT_TEMPERATURE)
      modelname = 'POP'

   case (LEAF_CARBON)
      modelname = 'CLM'

   case default
      ! unknown
      modelname = 'NULL'

end select

end subroutine which_model_obs

!------------------------------------------------------------------
! End of model_mod
!------------------------------------------------------------------

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
