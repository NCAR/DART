! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
!----------------------------------------------------------------
!>
!> This is the interface between the ROMS ocean model and DART.
!> There are 16 required public interfaces whose arguments CANNOT be changed.
!> There are potentially many more public routines that are typically
!> used by the converter programs. As the converter programs get phased out
!> with the impending native netCDF read/write capability, these extra
!> public interfaces may not need to be public.
!>
!> author: PENG XIU 12/2013 @ University of Maine
!>         peng.xiu@maine.edu
!>
!> based on subroutines from others work in the DART package
!> NOTE: For now, this code does not rotate ROMS vectors,
!>       so, rotate them in the observation data if needed
!>
!> subsequently modified by the DART team.
!>
!> \todo
!> @todo really check the land masking and _FillValue processing
!----------------------------------------------------------------

module model_mod

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, digits12, SECPERDAY, DEG2RAD, rad2deg, PI, &
                             MISSING_I, MISSING_R4, MISSING_R8, i4, i8

use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time, &
                             print_time, print_date,set_calendar_type,          &
                             operator(*),  operator(+), operator(-),            &
                             operator(>),  operator(<), operator(/),            &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, get_dist, get_close_maxdist_init,   &
                             get_close_obs_init, set_location,                  &
                             get_location, vert_is_height,write_location,      &
                             set_location_missing,query_location,               &
                             vert_is_level, vert_is_surface,                    &
                             loc_get_close_obs => get_close_obs, get_close_type,&
                             VERTISHEIGHT,VERTISSURFACE,VERTISUNDEF,VERTISLEVEL,&
                             VERTISPRESSURE,VERTISSCALEHEIGHT,horiz_dist_only

use    utilities_mod, only : register_module, error_handler,                    &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,       &
                             nc_check, do_output, to_upper,                     &
                             find_namelist_in_file, check_namelist_read,        &
                             open_file, file_exist, find_textfile_dims,         &
                             file_to_text, do_output,close_file

use     obs_kind_mod, only : KIND_TEMPERATURE, KIND_SALINITY, KIND_DRY_LAND,    &
                             KIND_U_CURRENT_COMPONENT,                          &
                             KIND_V_CURRENT_COMPONENT, KIND_SEA_SURFACE_HEIGHT, &
                             KIND_SEA_SURFACE_PRESSURE,                         &
                             KIND_POTENTIAL_TEMPERATURE,                        &
                             paramname_length,get_raw_obs_kind_index,           &
                             get_raw_obs_kind_name,get_obs_kind_var_type

use mpi_utilities_mod, only: my_task_id

use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use ensemble_manager_mod,  only : ensemble_type, map_pe_to_task, get_copy_owner_index, &
                                  get_var_owner_index

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain, get_model_variable_indices, &
                                  get_num_variables, get_index_start, &
                                  get_num_dims, get_domain_size, get_varid_from_kind, &
                                  get_dart_vector_index, state_structure_info

use typesizes
use netcdf

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
public :: get_model_size,                &
          adv_1step,                     &
          get_state_meta_data,           &
          model_interpolate,             &
          get_model_time_step,           &
          static_init_model,             &
          end_model,                     &
          init_time,                     &
          init_conditions,               &
          nc_write_model_atts,           &
          nc_write_model_vars,           &
          pert_model_copies,             &
          get_close_maxdist_init,        &
          get_close_obs_init,            &
          query_vert_localization_coord, &
          vert_convert,                  &
          construct_file_name_in,        &
          write_model_time,              &
          read_model_time,               &
          get_close_obs

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.
public :: restart_file_to_sv,         &
          sv_to_restart_file,         &
          get_model_restart_filename, &
          get_time_from_namelist,     &
          print_variable_ranges

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2, string3
logical, save :: module_initialized = .false.

! things which can/should be in the model_nml
logical  :: output_state_vector          = .false.
integer  :: assimilation_period_days     = 1
integer  :: assimilation_period_seconds  = 0
character(len=19) :: analysis_time       = '2001-01-01 06:00:00'
real(r8) :: model_perturbation_amplitude = 0.2
logical  :: update_dry_cell_walls        = .false.
integer  :: vert_localization_coord      = VERTISHEIGHT
integer  :: debug = 0   ! turn up for more and more debug messages
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: model_restart_filename = 'roms_restart.nc'
character(len=256) :: grid_definition_filename = 'roms_grid.nc'
real(r8) :: hc=50.0_r8

namelist /model_nml/  &
   analysis_time,               &
   output_state_vector,         &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   hc,                          &
   model_restart_filename,      &
   grid_definition_filename,    &
   vert_localization_coord,     &
   debug,                       &
   variables

! DART state vector contents are specified in the input.nml:&mpas_vars_nml namelist.
integer, parameter :: max_state_variables = 80
integer, parameter :: num_state_table_columns = 2
integer, parameter :: num_bounds_table_columns = 4
character(len=NF90_MAX_NAME) :: variables(max_state_variables * num_state_table_columns ) = ' '
character(len=NF90_MAX_NAME) :: roms_state_bounds(num_bounds_table_columns, max_state_variables ) = ' '
character(len=NF90_MAX_NAME) :: variable_table(max_state_variables, num_state_table_columns )

integer :: nfields   ! This is the number of variables in the DART state vector.

integer :: domain_id ! global variable for state_structure_mod routines

!> Everything needed to describe a variable. Basically all the metadata from
!> a netCDF file is stored here as well as all the information about where
!> the variable is stored in the DART state vector.
!>
!> @todo FIXME ... add a field for what kind of mask to use
!>

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=NF90_MAX_NAME), dimension(NF90_MAX_VAR_DIMS) :: dimname
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer :: xtype         ! netCDF variable type (NF90_double, etc.)
   integer :: numdims       ! number of dims - excluding TIME
   integer :: numvertical   ! number of vertical levels in variable
   integer :: numxi         ! number of horizontal locations (cell centers)
   integer :: numeta        ! number of horizontal locations (edges for velocity components)
   integer :: varsize       ! prod(dimlens(1:numdims))
   integer :: index1        ! location in dart state vector of first occurrence
   integer :: indexN        ! location in dart state vector of last  occurrence
   integer :: dart_kind
   character(len=paramname_length) :: kind_string
   logical  :: clamping     ! does variable need to be range-restricted before
   real(r8) :: range(2)     ! being stuffed back into analysis file.
   logical  :: out_of_range_fail  ! is out of range fatal if range-checking?
   real(r4) :: fill_value_r4
   real(r8) :: fill_value_r8
   logical  :: has_fill_value_r4
   logical  :: has_fill_value_r8
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar


! Grid parameters - the values will be read from a
! standard ROMS namelist and filled in here.

! nx, ny and nz are the size of the rho grids.
integer :: Nx = -1, Ny = -1, Nz = -1

integer :: Nxi_rho
integer :: Nxi_u
integer :: Nxi_v
integer :: Nxi_psi
integer :: Neta_rho
integer :: Neta_u
integer :: Neta_v
integer :: Neta_psi
integer :: Nxi_vert
integer :: Neta_vert
integer :: Ns_rho
integer :: Ns_w

real(r8), allocatable :: ULAT(:,:), ULON(:,:), &
                         TLAT(:,:), TLON(:,:), &
                         VLAT(:,:), VLON(:,:), &
                           PM(:,:),   PN(:,:), &
                         ANGL(:,:),   HT(:,:), ZC(:,:,:)

integer, parameter :: i2 = SELECTED_INT_KIND(2) ! need something to coerce to NF90_SHORT
integer(i2), allocatable :: mask_rho(:,:), &
                            mask_psi(:,:), &
                              mask_u(:,:), &
                              mask_v(:,:)

real(r8)        :: ocean_dynamics_timestep = 900.0_r4
type(time_type) :: model_timestep

integer :: model_size    ! the state vector length
real(r8), allocatable :: ens_mean(:)

!> Reshapes a part of the DART vector back to the original variable shape.
!> @todo FIXME Replaces the DART MISSING value with the original _FillValue value.

INTERFACE vector_to_prog_var
      MODULE PROCEDURE vector_to_1d_prog_var
      MODULE PROCEDURE vector_to_2d_prog_var
      MODULE PROCEDURE vector_to_3d_prog_var
END INTERFACE

!> Packs a ROMS variable into the DART vector.
!> @todo FIXME Replaces the original _FillValue value with the DART MISSING value.

INTERFACE prog_var_to_vector
      MODULE PROCEDURE prog_var_1d_to_vector
      MODULE PROCEDURE prog_var_2d_to_vector
      MODULE PROCEDURE prog_var_3d_to_vector
      MODULE PROCEDURE prog_var_4d_to_vector
END INTERFACE

!> Return the first and last index into the DART array for a specific variable.

INTERFACE get_index_range
      MODULE PROCEDURE get_index_range_int
      MODULE PROCEDURE get_index_range_string
END INTERFACE


contains


!-----------------------------------------------------------------------
! All the REQUIRED interfaces come first - by convention.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!>
!> Returns the size of the DART state vector (i.e. model) as an integer.
!> Required for all applications.
!>

function get_model_size()

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size


!-----------------------------------------------------------------------
!>
!> Does a single timestep advance of the model in a subroutine call.
!> This interface is only called if the namelist parameter
!> async is set to 0 in perfect_model_obs of filter or if the
!> program integrate_model is to be used to advance the model
!> state as a separate executable. If one of these options
!> is not going to be used (the model will only be advanced as
!> a separate model-specific executable), this can be a
!> NULL INTERFACE.
!>
!> NOTE: not supported for ROMS. Will intentionally generate a fatal error.
!>
!> @param x the model state before and after the model advance.
!> @param time the desired time at the end of the model advance.
!>

subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

if ( .not. module_initialized ) call static_init_model

write(string1,*)'Cannot advance ROMS with a subroutine call; async cannot equal 0'
write(string2,*)'Unsupported method for ROMS.'
call error_handler(E_ERR, 'adv_1step:', string1, &
                   source, revision, revdate, text2=string2)

end subroutine adv_1step


!-----------------------------------------------------------------------
!>
!> Given an integer index into the state vector structure, returns the
!> associated location. A second intent(out) optional argument kind
!> can be returned if the model has more than one type of field (for
!> instance temperature and zonal wind component). This interface is
!> required for all filter applications as it is required for computing
!> the distance between observations and state variables.
!>
!> @param state_handle DART ensemble handle
!> @param index_in the index into the DART state vector
!> @param location the location at that index
!> @param var_type the DART KIND at that index
!>

subroutine get_state_meta_data(state_handle, index_in, location, var_type)

type(ensemble_type), intent(in)  :: state_handle
integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

! Local variables

integer  :: nxp, nzp, iloc, vloc, nf, n,nyp,jloc
integer  :: myindx
integer  :: ivar
real(r8) :: depth

if ( .not. module_initialized ) call static_init_model

myindx = -1
nf     = -1

! Determine the right variable
FindIndex : do n = 1,nfields
    if( (progvar(n)%index1 <= index_in) .and. (index_in <= progvar(n)%indexN) ) THEN
      nf = n
      myindx = index_in - progvar(n)%index1 + 1
      exit FindIndex
    endif
enddo FindIndex

if (present(var_type)) then
   var_type = progvar(nf)%dart_kind
endif

if( myindx == -1 ) then
     write(string1,*) 'Problem, cannot find base_offset, index_in is: ', index_in
     call error_handler(E_ERR,'get_state_meta_data:',string1,source,revision,revdate)
endif

! Now that we know the variable, find the cell or edge

if (     progvar(nf)%numxi /= MISSING_I .AND. progvar(nf)%numeta /= MISSING_I ) then
   nyp = progvar(nf)%numxi
   nxp = progvar(nf)%numeta
else
     write(string1,*) 'ERROR, ',trim(progvar(nf)%varname),' is not defined on xi or eta'
     call error_handler(E_ERR,'get_state_meta_data:',string1,source,revision,revdate)
endif

call get_state_indices(progvar(nf)%dart_kind, myindx, iloc, jloc, vloc)

nzp  = progvar(nf)%numvertical
if(nzp==1) then
  depth=0.0
else
  depth= ZC(iloc,jloc,vloc)
endif

ivar = get_progvar_index_from_kind(progvar(nf)%dart_kind)

if (progvar(ivar)%kind_string=='KIND_SEA_SURFACE_HEIGHT') then
      location = set_location(TLON(iloc,jloc),TLAT(iloc,jloc), depth, VERTISSURFACE)
elseif (progvar(ivar)%kind_string=='KIND_U_CURRENT_COMPONENT') then
      location = set_location(ULON(iloc,jloc),ULAT(iloc,jloc), depth, VERTISHEIGHT)
elseif (progvar(ivar)%kind_string=='KIND_V_CURRENT_COMPONENT') then
      location = set_location(VLON(iloc,jloc),VLAT(iloc,jloc), depth, VERTISHEIGHT)
else
      location = set_location(TLON(iloc,jloc),TLAT(iloc,jloc), depth, VERTISHEIGHT)
endif

end subroutine get_state_meta_data


!-----------------------------------------------------------------------
!>
!> Model interpolate will interpolate any DART state variable
!> (i.e. S, T, U, V, Eta) to the given location given a state vector.
!> The type of the variable being interpolated is obs_type since
!> normally this is used to find the expected value of an observation
!> at some location. The interpolated value is returned in interp_vals
!> and istatus is 0 for success. NOTE: This is a workhorse routine and is
!> the basis for all the forward observation operator code.
!>
!> @param state_handle DART ensemble handle
!> @param ens_size DART ensemble size
!> @param location the location of interest
!> @param obs_type the DART KIND of interest
!> @param interp_val the estimated value of the DART state at the location
!>          of interest (the interpolated value).
!> @param istatus interpolation status ... 0 == success, /=0 is a failure
!>
!> @todo FIXME use some unique error code if the location is technically
!> outside the domain, i.e. an extrapolation. At some point it will be
!> useful to know if the interpolation failed because of some illegal
!> state as opposed to simply being outside the domain.

subroutine model_interpolate(state_handle, ens_size, location, obs_type, interp_val, istatus)

 type(ensemble_type), intent(in) :: state_handle
 integer,             intent(in) :: ens_size
 type(location_type), intent(in) :: location
 integer,             intent(in) :: obs_type
 integer,            intent(out) :: istatus(ens_size)
 real(r8),           intent(out) :: interp_val(ens_size) !< array of interpolated values

! Local storage
real(r8)       :: loc_array(3), llon, llat, lheight
integer        :: base_offset, top_offset
integer        :: ivar,obs_kind
integer        :: x_ind, y_ind,lat_bot, lat_top, lon_bot, lon_top
real(r8)       :: x_corners(4), y_corners(4), p(4)

if ( .not. module_initialized ) call static_init_model

! Successful istatus is 0
interp_val = MISSING_R8
istatus = 99

! Get the individual locations values
loc_array = get_location(location)
llon    = loc_array(1)
llat    = loc_array(2)
lheight = loc_array(3)

if( vert_is_height(location) ) then
   ! Nothing to do
elseif ( vert_is_surface(location) ) then
   ! Nothing to do
elseif (vert_is_level(location)) then
   ! convert the level index to an actual depth
   !ind = nint(loc_array(3))
   !if ( (ind < 1) .or. (ind > size(zc)) ) then
   !   lheight = zc(ind)
   !else
   !   istatus = 1
   !   return
   !endif
else   ! if pressure or undefined, we don't know what to do
   istatus = 7
   return
endif

obs_kind = obs_type
ivar = get_progvar_index_from_kind(obs_kind)

! Do horizontal interpolations for the appropriate levels
! Find the basic offset of this field

base_offset = progvar(ivar)%index1
top_offset  = progvar(ivar)%indexN

!print *, 'base offset now for ',trim(progvar(ivar)%varname),base_offset,top_offset

! For Sea Surface Height or Pressure don't need the vertical coordinate
!
if( vert_is_surface(location) ) then
   call lon_lat_interpolate(state_handle, ens_size, llon, llat, &
   obs_type, 1, interp_val, istatus)
   ! call lon_lat_interpolate(x(base_offset:top_offset), llon, llat, &
   ! obs_type, 1, interp_val, istatus)
   return
endif

!for 3D interpolation,
!1. find four corners
!2. do vertical interpolation at four corners
!3. do a spatial interpolation

! write(*,*)'TJH ', trim(progvar(ivar)%varname)//' '//trim(progvar(ivar)%kind_string), &
!                   llon, llat, lheight
! do istatus=1,progvar(ivar)%numdims
!    write(*,*)'TJH ',istatus, trim(progvar(ivar)%dimname(istatus)), &
!                                   progvar(ivar)%dimlens(istatus)
! enddo
istatus = 0

if(progvar(ivar)%kind_string == 'KIND_U_CURRENT_COMPONENT') then
   call get_reg_box_indices(llon, llat, ULON, ULAT, x_ind, y_ind)
   ! Getting corners for accurate interpolation
   call get_quad_corners(ULON, x_ind, y_ind, x_corners)
   call get_quad_corners(ULAT, x_ind, y_ind, y_corners)
elseif (progvar(ivar)%kind_string == 'KIND_V_CURRENT_COMPONENT') then
   call get_reg_box_indices(llon, llat, VLON, VLAT, x_ind, y_ind)
   ! Getting corners for accurate interpolation
   call get_quad_corners(VLON, x_ind, y_ind, x_corners)
   call get_quad_corners(VLAT, x_ind, y_ind, y_corners)
else
   call get_reg_box_indices(llon, llat, TLON, TLAT, x_ind, y_ind)
   ! Getting corners for accurate interpolation
   call get_quad_corners(TLON, x_ind, y_ind, x_corners)
   call get_quad_corners(TLAT, x_ind, y_ind, y_corners)
endif

lon_bot=x_ind
lat_bot=y_ind

! Find the indices to get the values for interpolating
lat_top = lat_bot + 1
if(lat_top > progvar(ivar)%numeta) then
   istatus = 2
   return
endif

lon_top = lon_bot + 1
if(lon_top > progvar(ivar)%numxi) then
   istatus = 2
   return
endif

!write(*,*)'model_mod: obs loc vs. model loc ', llon, llat, TLON(x_ind,y_ind),TLAT(x_ind,y_ind)

! If any of these fail, we can exit (fail) immediately.
! The interp_val and istatus values are initially set to a failed value, so just use those.

call vert_interpolate(state_handle,ens_size,lon_bot,lat_bot,obs_kind,lheight, p(1))
if (p(1) == MISSING_R8) return

call vert_interpolate(state_handle,ens_size,lon_top,lat_bot,obs_kind,lheight, p(2))
if (p(2) == MISSING_R8) return

call vert_interpolate(state_handle,ens_size,lon_top,lat_top,obs_kind,lheight, p(3))
if (p(3) == MISSING_R8) return

call vert_interpolate(state_handle,ens_size,lon_bot,lat_top,obs_kind,lheight, p(4))
if (p(4) == MISSING_R8) return

call quad_bilinear_interp(llon, llat, x_corners, y_corners, p, ens_size, interp_val)
istatus = 0  ! If we get this far, all good.

!write(*,*) 'Ending model interpolate ...'

end subroutine model_interpolate


!-----------------------------------------------------------------------
!>
!> Returns the the time step of the model; the smallest increment in
!> time that the model is capable of advancing the ROMS state.
!>

function get_model_time_step()

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = model_timestep

end function get_model_time_step


!-----------------------------------------------------------------------
!>
!> Called to do one time initialization of the model.
!> In this case, it reads in the grid information, the namelist
!> containing the variables of interest, where to get them, their size,
!> their associated DART KIND, etc.
!>
!> In addition to harvesting the model metadata (grid,
!> desired model advance step, etc.), it also fills a structure
!> containing information about what variables are where in the DART
!> framework.

subroutine static_init_model()

integer :: iunit, io,index1,indexN,ivar
integer :: ss, dd,i,TimeDimID
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname,dimname
character(len=paramname_length)       :: kind_string
integer :: ncid, VarID, numdims, varsize, dimlen

if ( module_initialized ) return

! The Plan:
!
! * read in the grid sizes from grid file
! * allocate space, and read in actual grid values
! * figure out model timestep
! * Compute the model size.
! * set the index numbers where the field types change

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
call error_handler(E_MSG,'static_init_model:','model_nml values are',' ',' ',' ')
if (do_output()) write(logfileunit, nml=model_nml)
if (do_output()) write(     *     , nml=model_nml)

call set_calendar_type( calendar )

! Set the time step ... causes ROMS namelists to be read.
! Ensures model_timestep is multiple of 'ocean_dynamics_timestep'
!> @todo FIXME 'ocean_dynamics_timestep' could/should be gotten from ROMS

model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model:',string1,source,revision,revdate)

! Get the ROMS grid -- sizes and variables.
call get_grid_dimensions()
call get_grid()

call nc_check( nf90_open(trim(model_restart_filename), NF90_NOWRITE, ncid), &
                  'static_init_model', 'open '//trim(model_restart_filename))

call verify_variables( variables, ncid, model_restart_filename, &
                            nfields, variable_table )

domain_id = add_domain(model_restart_filename, nfields, &
                    var_names = variable_table(:,1))

call state_structure_info(domain_id)

TimeDimID = find_time_dimension( ncid, model_restart_filename )

call nc_check( nf90_close(ncid), &
                  'static_init_model', 'close '//trim(model_restart_filename))

model_size = get_domain_size(domain_id)

allocate( ens_mean(model_size) )

end subroutine static_init_model


!-----------------------------------------------------------------------
!>
!> Does any shutdown and clean-up needed for model.
!>

subroutine end_model()

! good style ... perhaps you could deallocate stuff (from static_init_model?).
! deallocate(state_loc)

if (allocated(ULAT)) deallocate(ULAT)
if (allocated(ULON)) deallocate(ULON)
if (allocated(TLAT)) deallocate(TLAT)
if (allocated(TLON)) deallocate(TLON)
if (allocated(VLAT)) deallocate(VLAT)
if (allocated(VLON)) deallocate(VLON)
if (allocated(HT))   deallocate(HT)
if (allocated(PM))   deallocate(PM)
if (allocated(PN))   deallocate(PN)
if (allocated(ANGL)) deallocate(ANGL)
if (allocated(ZC))   deallocate(ZC)

end subroutine end_model


!-----------------------------------------------------------------------
!>
!> Companion interface to init_conditions. Returns a time that is somehow
!> appropriate for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter
!> start_from_restart is set to .false. in the program perfect_model_obs.
!> If this option is not to be used in perfect_model_obs, or if no
!> synthetic data experiments using perfect_model_obs are planned,
!> this can be a NULL INTERFACE.
!>
!> NOTE: Since ROMS cannot start in this manner,
!> DART will intentionally generate a fatal error.
!>
!> @param time the time to associate with the initial state
!>

subroutine init_time(time)

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

time = set_time(0,0)

write(string1,*) 'Cannot initialize ROMS time via subroutine call; start_from_restart cannot be F'
write(string2,*)'Unsupported method for ROMS.'
call error_handler(E_ERR, 'init_time:', string1, &
                   source, revision, revdate, text2=string2)

end subroutine init_time


!-----------------------------------------------------------------------
!>
!> Returns a model state vector, x, that is some sort of appropriate
!> initial condition for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter
!> start_from_restart is set to .false. in the program perfect_model_obs.
!> If this option is not to be used in perfect_model_obs, or if no
!> synthetic data experiments using perfect_model_obs are planned,
!> this can be a NULL INTERFACE.
!>
!> NOTE: This is not supported for ROMS and will generate a FATAL ERROR.
!>       However, this is a required interface - so it must be present.
!>
!> @param x the ROMS initial conditions
!>

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

x = MISSING_R8

write(string1,*)'Cannot initialize ROMS state via subroutine call.'
write(string2,*)'namelist "start_from_restart" cannot be F'
call error_handler(E_ERR, 'init_conditions:', string1, &
                   source, revision, revdate, text2=string2)

end subroutine init_conditions


!-----------------------------------------------------------------------
!>
!> Writes the model-specific attributes to a DART 'diagnostic' netCDF file.
!> This includes coordinate variables and some metadata, but NOT the
!> actual DART state. That may be done multiple times in nc_write_model_vars()
!>
!> @param ncFileID the netCDF handle of the DART diagnostic file opened by
!>                 assim_model_mod:init_diag_output
!> @param model_writes_state have the state structure write out all of the
!>                 state variables
!> @param ierr status ... 0 == all went well, /= 0 failure

function nc_write_model_atts( ncFileID, model_writes_state ) result (ierr)

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

integer, intent(in)  :: ncFileID          ! netCDF file identifier
logical, intent(out) :: model_writes_state ! netCDF file identifier
integer              :: ierr              ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

! variables if we just blast out one long state vector

integer :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)

integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

! variables if we parse the state vector into prognostic variables.

! for the dimensions and coordinate variables
integer :: nxirhoDimID, nxiuDimID, nxivDimID, nxipsiDimID, nxivertDimID
integer :: netarhoDimID, netauDimID, netavDimID, netapsiDimID, netavertDimID
integer :: nsrhoDimID, nswDimID

! for the prognostic variables
integer :: ivar, VarID

! local variables

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1
character(len=NF90_MAX_NAME) :: varname
integer, dimension(NF90_MAX_VAR_DIMS) :: mydimids
integer :: myndims

character(len=256) :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly
model_writes_state = .true.

! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.

write(filename,*) 'ncFileID', ncFileID

! make sure ncFileID refers to an open netCDF file,
! and then put into define mode.

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
                                   'nc_write_model_atts', 'inquire '//trim(filename))
call nc_check(nf90_Redef(ncFileID),'nc_write_model_atts',   'redef '//trim(filename))

! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension.
! Our job is create the 'model size' dimension.

call nc_check(nf90_inq_dimid(ncid=ncFileID, name='copy', dimid=MemberDimID), &
                           'nc_write_model_atts', 'copy dimid '//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid=  TimeDimID), &
                           'nc_write_model_atts', 'time dimid '//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)'Time Dimension ID ',TimeDimID, &
             ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts:', string1, source, revision, revdate)
endif

! Define the model size / state variable dimension / whatever ...

call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable', len=model_size, &
        dimid = StateVarDimID),'nc_write_model_atts', 'state def_dim '//trim(filename))

! Write Global Attributes

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'creation_date' ,str1    ), &
           'nc_write_model_atts', 'creation put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_source'  ,source  ), &
           'nc_write_model_atts', 'source put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revision',revision), &
           'nc_write_model_atts', 'revision put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revdate' ,revdate ), &
           'nc_write_model_atts', 'revdate put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'ROMS' ), &
           'nc_write_model_atts', 'model put '//trim(filename))

! Here is the extensible part. The simplest scenario is to output the state vector,
! parsing the state vector into model-specific parts is complicated, and you need
! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
! complicated part.

if ( output_state_vector ) then

   ! Create a variable for the state vector
   ! Define the actual (3D) state vector, which gets filled as time goes on ...

   call nc_check(nf90_def_var(ncid=ncFileID, name='state', xtype=nf90_real, &
                 dimids=(/StateVarDimID,MemberDimID,unlimitedDimID/),varid=StateVarID),&
                 'nc_write_model_atts','state def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarID,'long_name','model state or fcopy'),&
                 'nc_write_model_atts', 'state long_name '//trim(filename))

   call nc_check(nf90_enddef(ncFileID),'nc_write_model_atts','state enddef '//trim(filename))

else

   ! We need to output the prognostic variables.
   ! Define the new dimensions IDs

   call nc_check(nf90_def_dim(ncid=ncFileID, name='xi_rho',  len = Nxi_rho, &
        dimid = nxirhoDimID),'nc_write_model_atts', 'xi_rho def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='eta_rho', len = Neta_rho,&
        dimid = netarhoDimID),'nc_write_model_atts', 'eta_rho def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='s_rho',   len = Ns_rho,&
        dimid = nsrhoDimID),'nc_write_model_atts', 's_rho def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='s_w',     len = Ns_w,&
        dimid = nswDimID),'nc_write_model_atts','s_w def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='xi_u',    len = Nxi_u,&
        dimid = nxiuDimID),'nc_write_model_atts', 'xi_u def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='xi_v',    len = Nxi_v,&
        dimid = nxivDimID),'nc_write_model_atts', 'xi_v def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='eta_u',   len = Neta_u,&
        dimid = netauDimID),'nc_write_model_atts', 'eta_u def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='eta_v',   len = Neta_v,&
        dimid = netavDimID),'nc_write_model_atts', 'eta_v def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='xi_psi',   len = Nxi_psi,&
        dimid = nxipsiDimID),'nc_write_model_atts', 'xi_psi def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='eta_psi',   len = Neta_psi,&
        dimid = netapsiDimID),'nc_write_model_atts', 'eta_psi def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='xi_vert',   len = Nxi_vert,&
        dimid = nxivertDimID),'nc_write_model_atts', 'xi_vert def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='eta_vert',   len = Neta_vert,&
        dimid = netavertDimID),'nc_write_model_atts', 'eta_vert def_dim '//trim(filename))

   ! Create the (empty) Coordinate Variables and the Attributes

   call nc_check(nf90_def_var(ncFileID,name='lon_rho', xtype=nf90_double, &
                 dimids=(/ nxirhoDimID, netarhoDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lon_rho def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'rho longitudes'), &
                 'nc_write_model_atts', 'lon_rho long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'lon_rho units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='lat_rho', xtype=nf90_double, &
                 dimids=(/ nxirhoDimID, netarhoDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lat_rho def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'rho latitudes'), &
                 'nc_write_model_atts', 'lat_rho long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_north'), &
                 'nc_write_model_atts', 'lat_rho units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='lon_u', xtype=nf90_double, &
                 dimids=(/ nxiuDimID, netauDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lon_u def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'u longitudes'), &
                 'nc_write_model_atts', 'lon_u long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'lon_u units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='lat_u', xtype=nf90_double, &
                 dimids=(/ nxiuDimID, netauDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lat_u def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'u latitudes'), &
                 'nc_write_model_atts', 'lat_u long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_north'), &
                 'nc_write_model_atts', 'lat_u units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='lon_v', xtype=nf90_double, &
                 dimids=(/ nxivDimID, netavDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lon_v def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'v longitudes'), &
                 'nc_write_model_atts', 'lon_v long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'lon_v units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='lat_v', xtype=nf90_double, &
                 dimids=(/ nxivDimID, netavDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lat_v def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'v latitudes'), &
                 'nc_write_model_atts', 'lat_v long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_north'), &
                 'nc_write_model_atts', 'lat_v units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='z_c', xtype=nf90_double, &
                 dimids=(/ nxirhoDimID, netarhoDimID, nsrhoDimID /), varid=VarID),&
                 'nc_write_model_atts', 'z_c def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'z at rho'), &
                 'nc_write_model_atts', 'z_c long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'm'), &
                 'nc_write_model_atts', 'z_c units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='mask_rho', xtype=nf90_short, &
                 dimids=(/ nxirhoDimID, netarhoDimID /), varid=VarID),&
                 'nc_write_model_atts', 'mask_rho def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'mask on RHO-points'), &
                 'nc_write_model_atts', 'mask_rho long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'comment', '0==land, 1==water'), &
                 'nc_write_model_atts', 'mask_rho comment '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='mask_psi', xtype=nf90_short, &
                 dimids=(/ nxipsiDimID, netapsiDimID /), varid=VarID),&
                 'nc_write_model_atts', 'mask_psi def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'mask on psi-points'), &
                 'nc_write_model_atts', 'mask_psi long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'comment', '0==land, 1==water'), &
                 'nc_write_model_atts', 'mask_psi units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='mask_u', xtype=nf90_short, &
                 dimids=(/ nxiuDimID, netauDimID /), varid=VarID),&
                 'nc_write_model_atts', 'mask_u def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'mask on U-points'), &
                 'nc_write_model_atts', 'mask_u long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'comment', '0==land, 1==water'), &
                 'nc_write_model_atts', 'mask_u units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='mask_v', xtype=nf90_short, &
                 dimids=(/ nxivDimID, netavDimID /), varid=VarID),&
                 'nc_write_model_atts', 'mask_v def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'mask on V-points'), &
                 'nc_write_model_atts', 'mask_v long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'comment', '0==land, 1==water'), &
                 'nc_write_model_atts', 'mask_v units '//trim(filename))

   ! Create the (empty) Prognostic Variables and the Attributes

   do ivar=1, nfields

      varname = trim(progvar(ivar)%varname)
      string1 = trim(filename)//' '//trim(varname)

      ! match shape of the variable to the dimension IDs

      call define_var_dims(ncFileID, ivar, MemberDimID, unlimitedDimID, myndims, mydimids)

      ! define the variable and set the attributes
      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(varname), &
              xtype = progvar(ivar)%xtype, dimids = mydimids(1:myndims), varid=VarID), &
              'nc_write_model_atts', trim(string1)//' def_var' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', &
              trim(progvar(ivar)%long_name)), &
              'nc_write_model_atts', trim(string1)//' put_att long_name' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'DART_kind', &
              trim(progvar(ivar)%kind_string)), &
              'nc_write_model_atts', trim(string1)//' put_att dart_kind' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'units', &
              trim(progvar(ivar)%units)), &
              'nc_write_model_atts', trim(string1)//' put_att units' )

      if     (progvar(ivar)%has_fill_value_r4) then
         call nc_check(nf90_put_att(ncFileID, VarID, '_FillValue', &
                                         progvar(ivar)%fill_value_r4), &
              'nc_write_model_atts', trim(string1)//' put_att _FillValue' )
      elseif (progvar(ivar)%has_fill_value_r8) then
         call nc_check(nf90_put_att(ncFileID, VarID, '_FillValue', &
                                         progvar(ivar)%fill_value_r8), &
              'nc_write_model_atts', trim(string1)//' put_att _FillValue' )
      endif

   enddo

   ! Finished with dimension/variable definitions, must end 'define' mode to fill.

   call nc_check(nf90_enddef(ncFileID), 'prognostic enddef '//trim(filename))

   ! Fill the coordinate variables that DART needs and has locally

   call nc_check(NF90_inq_varid(ncFileID, 'lon_rho', VarID), &
                 'nc_write_model_atts', 'lon_rho inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, TLON ), &
                'nc_write_model_atts', 'lon_rho put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'lat_rho', VarID), &
                 'nc_write_model_atts', 'lat_rho inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, TLAT ), &
                'nc_write_model_atts', 'lat_rho put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'lon_u', VarID), &
                 'nc_write_model_atts', 'lon_u inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, ULON ), &
                'nc_write_model_atts', 'lon_u put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'lat_u', VarID), &
                 'nc_write_model_atts', 'lat_u inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, ULAT ), &
                'nc_write_model_atts', 'lat_u put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'lon_v', VarID), &
                 'nc_write_model_atts', 'lon_v inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, VLON ), &
                'nc_write_model_atts', 'lon_v put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'lat_v', VarID), &
                 'nc_write_model_atts', 'lat_v inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, VLAT ), &
                'nc_write_model_atts', 'lat_v put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'z_c', VarID), &
                 'nc_write_model_atts', 'z_c inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, ZC ), &
                'nc_write_model_atts', 'z_c put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'mask_rho', VarID), &
                 'nc_write_model_atts', 'mask_rho inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, mask_rho ), &
                'nc_write_model_atts', 'mask_rho put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'mask_psi', VarID), &
                 'nc_write_model_atts', 'mask_psi inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, mask_psi ), &
                'nc_write_model_atts', 'mask_psi put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'mask_u', VarID), &
                 'nc_write_model_atts', 'mask_u inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, mask_u ), &
                'nc_write_model_atts', 'mask_u put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'mask_v', VarID), &
                 'nc_write_model_atts', 'mask_v inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, mask_v ), &
                'nc_write_model_atts', 'mask_v put_var '//trim(filename))

   endif

! Flush the buffer and leave netCDF file open
call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts


!-----------------------------------------------------------------------
!>
!> With each assimilation cycle, the DART prior and posterior files get
!> inserted into the DART diagnostic files. This routine appends the new
!> states into the unlimited dimension slot.
!>
!> @param ncFileID the netCDF file ID of the DART diagnostic file in question
!> @param state_vec the DART state to insert into the diagnostic file
!> @param copyindex the 'copy' index ... ensemble mean, member 23, etc.
!> @param timeindex the index into the unlimited (time) dimension
!> @param ierr error code. All errors are fatal. 0 == success.

function nc_write_model_vars( ncFileID, state_vec, copyindex, timeindex ) result (ierr)

! TJH 29 Aug 2011 -- all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
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

integer,  intent(in) :: ncFileID
real(r8), intent(in) :: state_vec(:)
integer,  intent(in) :: copyindex
integer,  intent(in) :: timeindex
integer              :: ierr

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME)          :: varname
integer :: i, ivar, VarID, ncNdims, dimlen
integer :: TimeDimID, CopyDimID

real(r8), allocatable, dimension(:)       :: data_1d_array
real(r8), allocatable, dimension(:,:)     :: data_2d_array
real(r8), allocatable, dimension(:,:,:)   :: data_3d_array

character(len=256) :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.

write(filename,*) 'ncFileID', ncFileID

! make sure ncFileID refers to an open netCDF file,

call nc_check(nf90_inq_dimid(ncFileID, 'copy', dimid=CopyDimID), &
            'nc_write_model_vars', 'inq_dimid copy '//trim(filename))

call nc_check(nf90_inq_dimid(ncFileID, 'time', dimid=TimeDimID), &
            'nc_write_model_vars', 'inq_dimid time '//trim(filename))

if ( output_state_vector ) then

   call nc_check(NF90_inq_varid(ncFileID, 'state', VarID), &
                 'nc_write_model_vars', 'state inq_varid '//trim(filename))
   call nc_check(NF90_put_var(ncFileID,VarID,state_vec,start=(/1,copyindex,timeindex/)),&
                 'nc_write_model_vars', 'state put_var '//trim(filename))

else

   ! We need to process the prognostic variables.

   do ivar = 1,nfields

      varname = trim(progvar(ivar)%varname)
      string2 = trim(filename)//' '//trim(varname)

      ! Ensure netCDF variable is conformable with progvar quantity.
      ! The TIME and Copy dimensions are intentionally not queried
      ! by looping over the dimensions stored in the progvar type.

      call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'nc_write_model_vars', 'inq_varid '//trim(string2))

      call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
            'nc_write_model_vars', 'inquire '//trim(string2))

      mystart = 1   ! These are arrays, actually
      mycount = 1
      DimCheck : do i = 1,progvar(ivar)%numdims

         write(string1,'(a,i2,A)') 'inquire dimension ',i,trim(string2)
         call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
               'nc_write_model_vars', trim(string1))

         if ( dimlen /= progvar(ivar)%dimlens(i) ) then
            write(string1,*) trim(string2),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
            write(string2,*)' but it should be.'
            call error_handler(E_ERR, 'nc_write_model_vars:', trim(string1), &
                            source, revision, revdate, text2=trim(string2))
         endif

         mycount(i) = dimlen

      enddo DimCheck

     ! FIXME - wouldn't hurt to make sure each of these match something.
     !         could then eliminate the if ncndims /= xxx checks below.

      where(dimIDs == CopyDimID) mystart = copyindex
      where(dimIDs == CopyDimID) mycount = 1
      where(dimIDs == TimeDimID) mystart = timeindex
      where(dimIDs == TimeDimID) mycount = 1

      if (     progvar(ivar)%numdims == 1 ) then

         if ( ncNdims /= 3 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 3 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars:', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_1d_array( progvar(ivar)%dimlens(1) ) )
         call vector_to_prog_var(state_vec, ivar, data_1d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_1d_array)

      elseif ( progvar(ivar)%numdims == 2 ) then

         if ( ncNdims /= 4 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 4 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars:', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_2d_array( progvar(ivar)%dimlens(1),  &
                                 progvar(ivar)%dimlens(2) ))
         call vector_to_prog_var(state_vec, ivar, data_2d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_2d_array)

      elseif ( progvar(ivar)%numdims == 3) then

         if ( ncNdims /= 5 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 5 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars:', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_3d_array( progvar(ivar)%dimlens(1), &
                                 progvar(ivar)%dimlens(2), &
                                 progvar(ivar)%dimlens(3)))
         call vector_to_prog_var(state_vec, ivar, data_3d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_3d_array)

      else

         write(string1,*)'no support (yet) for 4d fields'
         call error_handler(E_ERR, 'nc_write_model_vars:', string1, &
                            source, revision, revdate)

      endif
   enddo
endif

! Flush the buffer and leave netCDF file open

call nc_check(nf90_sync(ncFileID), 'nc_write_model_vars', 'sync '//trim(filename))

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars

!-----------------------------------------------------------------------
!>
!> Perturbs a model state for generating initial ensembles.
!> The perturbed state is returned in pert_state.
!> A model may choose to provide a NULL INTERFACE by returning
!> .false. for the interf_provided argument. This indicates to
!> the filter that if it needs to generate perturbed states, it
!> may do so by adding a perturbation to each model state
!> variable independently. The interf_provided argument
!> should be returned as .true. if the model wants to do its own
!> perturbing of states.
!>
!> @param state_ens_handle the DART state ensemble handle
!> @param ens_size ensemble size
!> @param pert_amp perturbation amplitude
!> @param interf_provided logical flag that indicates that this routine
!>               is unique for ROMS. TRUE means this routine will
!>               somehow create the perturbed state, FALSE means
!>               the default perturb routine will be used.
!>

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

 type(ensemble_type), intent(inout) :: state_ens_handle
 integer,             intent(in)    :: ens_size
 real(r8),            intent(in)    :: pert_amp
 logical,             intent(out)   :: interf_provided

! Perturbs state copies for generating initial ensembles.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding a perturbation to each model state 
! variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.

! Storage for a random sequence for perturbing a single initial state

if ( .not. module_initialized ) call static_init_model

interf_provided = .false.

!>@ todo FIXME : do we need to take into consideration the land mask when
!>               perturbing the state copies??
write(string1,*)'..  WARNING: pert_model_copies() not finished.'
write(string2,*)'WARNING: Does not respect MISSING values.'
write(string3,*)'WARNING: Fix before using.'
call error_handler(E_ERR,'pert_model_copies:', string1, text2=string2, text3=string3)

end subroutine pert_model_copies


!-----------------------------------------------------------------------
!>
!> Given a DART location (referred to as "base") and a set of candidate
!> locations and kinds (obs, obs_kind); returns the subset close to the
!> "base", their indices, and their distances to the "base" ...
!>
!> @param gc precomputed 'get_close_type' to speed up candidate selection
!> @param base_obs_loc location of the observation in question
!> @param base_obs_kind DART KIND of observation in question
!> @param locs array of comparison locations
!> @param loc_kind matching array of KINDs for the comparison locations
!> @param num_close the number of locs locations that are within the prespecified distance (information contained in 'gc')
!> @param close_ind the indices of the locs locations that are 'close'
!> @param dist the distances of each of the close locations.
!> @param state_handle handle to the dart state
!>

subroutine get_close_obs(gc, base_obs_loc, base_obs_type, &
                         locs, loc_kind, num_close, close_ind, dist, state_handle)
 type(ensemble_type),               intent(in)    :: state_handle
 type(get_close_type),              intent(in)    :: gc
 type(location_type),               intent(inout) :: base_obs_loc
 integer,                           intent(in)    :: base_obs_type
 type(location_type), dimension(:), intent(in)    :: locs
 integer,             dimension(:), intent(in)    :: loc_kind
 integer,                           intent(out)   :: num_close
 integer,             dimension(:), intent(out)   :: close_ind
 real(r8), OPTIONAL,  dimension(:), intent(out)   :: dist !does this need to be optional? It is not in WRF

! Note that both base_obs_loc and locs are intent(inout), meaning that these
! locations are possibly modified here and returned as such to the calling routine.
! The calling routine is always filter_assim and these arrays are local arrays
! within filter_assim. In other words, these modifications will only matter within
! filter_assim, but will not propagate backwards to filter.

integer                :: base_obs_kind ! for sanity
integer                :: ztypeout
integer                :: t_ind, istatus1, istatus2, k
integer                :: base_which, local_obs_which
real(r8), dimension(3) :: base_llv, local_obs_llv   ! lon/lat/vert
type(location_type)    :: local_obs_loc

real(r8) ::  hor_dist
hor_dist = 1.0e9_r8

if ( .not. module_initialized ) call static_init_model

! Initialize variables to missing status

base_obs_kind = base_obs_type

num_close = 0
close_ind = -99
dist      = 1.0e9_r8   !something big and positive (far away) in radians
istatus1  = 0
istatus2  = 0

base_llv = get_location(base_obs_loc)             ! lon/lat/vert
base_which = nint(query_location(base_obs_loc))   ! base_obs_loc%which_vert

ztypeout = vert_localization_coord

if (.not. horiz_dist_only) then
  if (base_llv(3) == MISSING_R8) then
     istatus1 = 1
  else if (base_which /= vert_localization_coord) then
      !>@todo FIXME : JPH. vert_convert needs some work... it's not doing anything
      !>              at the moment.
      call vert_convert(state_handle, base_obs_loc, base_obs_kind, istatus1)
      if(debug > 1) then
         call write_location(0,base_obs_loc,charstring=string1)
         call error_handler(E_MSG, 'get_close_obs: base_obs_loc', string1, &
                            source, revision, revdate)
      endif
  endif
endif

if (istatus1 == 0) then
   call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, locs, loc_kind, &
                          num_close, close_ind)

    do k = 1, num_close

      t_ind = close_ind(k)
      local_obs_loc   = locs(t_ind)
      local_obs_which = nint(query_location(local_obs_loc))

      if (.not. horiz_dist_only) then
          if (local_obs_which /= vert_localization_coord) then
              call vert_convert(state_handle, local_obs_loc, loc_kind(t_ind), istatus2)
          else
              istatus2 = 0
          endif
      endif

      local_obs_llv = get_location(local_obs_loc)

      if (( (.not. horiz_dist_only)           .and. &
            (local_obs_llv(3) == MISSING_R8)) .or.  &
            (istatus2 /= 0)                   ) then
            dist(k) = 1.0e9_r8
      else

       ! if (get_obs_kind_var_type(base_obs_kind) == loc_kind(t_ind)) then
           dist(k) = get_dist(base_obs_loc, local_obs_loc, base_obs_kind, loc_kind(t_ind))
       ! else
       !    dist(k) = 1.0e9_r8
       ! endif

      endif

    enddo

endif


end subroutine get_close_obs


!--------------------------------------------------------------------
!>
!> This is used to pass the vertical localization coordinate 
!> to assim_tools_mod
!>

function query_vert_localization_coord()

integer :: query_vert_localization_coord

query_vert_localization_coord = VERTISHEIGHT

end function query_vert_localization_coord

!-----------------------------------------------------------------------
!>
!> This subroutine converts a given ob/state vertical coordinate to
!> the vertical localization coordinate type requested through the
!> model_mod namelist. This is used in filter_assim(). The vertical
!> conversion is done using the mean state.
!>
!> Notes: (1) obs_kind is only necessary to check whether the ob
!>            is an identity observation.
!>
!>        (2) This subroutine can convert both obs' and state points'
!>            vertical coordinates. Remember that state points get
!>            their DART location information from get_state_meta_data
!>            which is called by filter_assim during the assimilation
!>            process.
!>
!>        (3) state_handle contains relevant DART state information for 
!>            carrying out computations necessary for the vertical coordinate
!>            transformations. As the vertical coordinate is only used
!>            in distance computations, this is actually the "expected"
!>            vertical coordinate, so that computed distance is the
!>            "expected" distance. Thus, under normal circumstances,
!>            state_handle that is supplied to convert_vert should be the
!>            ensemble mean. Nevertheless, the subroutine has the
!>            functionality to operate on any DART state vector that
!>            is supplied to it.
!>
!> @param state_handle handle to the dart state
!> @param location
!> @param obs_kind 
!> @param istatus
!>

subroutine vert_convert(state_handle, location, obs_kind, istatus)

type(ensemble_type), intent(in)  :: state_handle
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_kind
integer,             intent(out) :: istatus

!>@ todo FIXME : JPH. We need to make sure we are converting the observations
!>               vertical values to be the same as the state properly.

istatus = 0


end subroutine vert_convert


!--------------------------------------------------------------------
!>
!> construct restart file name for reading
!>
!> @param stub stub of restart file
!> @param domain for multiple domains in your state
!> @param copy copy number in ensemble handle
!>

function construct_file_name_in(stub, domain, copy)

character(len=256), intent(in) :: stub
integer,            intent(in) :: domain
integer,            intent(in) :: copy
character(len=256) :: construct_file_name_in

! stub is found in input.nml io_filename_nml

!>@todo FIXME : JPH. Need to construct filename properly!!!!

! write(construct_file_name_in, '(A, i4.4, A)') trim(stub), copy, ".nc"


end function construct_file_name_in


!-----------------------------------------------------------------------
!>
!> writes the time of the current state and (optionally) the time
!> to be conveyed to ROMS to dictate the length of the forecast.
!> This file is then used by scripts to modify the ROMS run.
!> The format in the time information is totally at your discretion.
!>
!> @param ncfile_out name of the file (name is set by dart_to_roms_nml:ncfile_out)
!> @param model_time the current time of the model state
!> @param adv_to_time the time in the future of the next assimilation.
!>

subroutine write_model_time(ncfile_out, model_time, adv_to_time)
integer,         intent(in)           :: ncfile_out
type(time_type), intent(in)           :: model_time
type(time_type), intent(in), optional :: adv_to_time

! integer :: iunit
! character(len=19) :: timestring
! type(time_type)   :: deltatime

if ( .not. module_initialized ) call static_init_model

!>@todo FIXME : JPH. Need to write the model time properly to
!>                   the netcdf file

! iunit = open_file(ncfile_out, action='write')
! 
! timestring = time_to_string(model_time)
! write(iunit, '(A)') timestring
! 
! if (present(adv_to_time)) then
!    timestring = time_to_string(adv_to_time)
!    write(iunit, '(A)') timestring
! 
!    deltatime = adv_to_time - model_time
!    timestring = time_to_string(deltatime, interval=.true.)
!    write(iunit, '(A)') timestring
! endif
! 
! call close_file(iunit)

end subroutine write_model_time

!--------------------------------------------------------------------
!>
!> read the time from the input file
!>
!> @param filename name of file that contains the time
!>

function read_model_time(filename)

character(len=256) :: filename
type(time_type) :: read_model_time

!>@todo FIXME : JPH. Need to read the model time properly from
!>                   the netcdf file. This is an example from POP

! integer :: ncid !< netcdf file id
! integer :: iyear, imonth, iday, ihour, iminute, isecond
! 
! if ( .not. module_initialized ) call static_init_model
! 
! if ( .not. file_exist(filename) ) then
!    write(msgstring,*) 'cannot open file ', trim(filename),' for reading.'
!    call error_handler(E_ERR,'read_model_time',msgstring,source,revision,revdate)
! endif
! 
! call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
!                   'read_model_time', 'open '//trim(filename))
! call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iyear'  , iyear), &
!                   'read_model_time', 'get_att iyear')
! call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'imonth' , imonth), &
!                   'read_model_time', 'get_att imonth')
! call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iday'   , iday), &
!                   'read_model_time', 'get_att iday')
! call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'ihour'  , ihour), &
!                   'read_model_time', 'get_att ihour')
! call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iminute', iminute), &
!                   'read_model_time', 'get_att iminute')
! call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'isecond', isecond), &
!                   'read_model_time', 'get_att isecond')
! 
! ! FIXME: we don't allow a real year of 0 - add one for now, but
! ! THIS MUST BE FIXED IN ANOTHER WAY!
! if (iyear == 0) then
!   call error_handler(E_MSG, 'read_model_time', &
!                      'WARNING!!!   year 0 not supported; setting to year 1')
!   iyear = 1
! endif
! 
! read_model_time = set_date(iyear, imonth, iday, ihour, iminute, isecond)


end function read_model_time


!-----------------------------------------------------------------------
! The remaining PUBLIC interfaces come next.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!>
!> Reads the last timestep of the ROMS variables and packs them into a
!> DART state vector. At present, only the last timestep is considered.
!> @todo determining the timestep of the DESIRED time might be a nice extension.
!>
!> @param filename the ROMS file that contains the variables of interest.
!> @param state_vector the DART state vector
!> @param last_file_time the last time in the file.
!>
!> @todo FIXME last_file_time is provided as a convenience, not sure
!>       I am decoding the information correctly given the test file
!>       that I had. It would be nice to get an accurate time and ensure
!>       that the desired time is being read.

subroutine restart_file_to_sv(filename, state_vector, last_file_time)

character(len=*), intent(in)  :: filename
real(r8),         intent(out) :: state_vector(:)
type(time_type),  intent(out) :: last_file_time

! temp space to hold data while we are reading it
integer  :: ndim1, ndim2, ndim3,ndim4
integer  :: indx, iostatus, ivar
real(r8), allocatable, dimension(:)       :: data_1d_array
real(r8), allocatable, dimension(:,:)     :: data_2d_array
real(r8), allocatable, dimension(:,:,:)   :: data_3d_array
real(r8), allocatable, dimension(:,:,:,:) :: data_4d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME) :: varname
integer :: VarID, ncNdims, dimlen
integer :: ncid, nDimensions, nVariables, nAttributes, unlimitedDimID

integer :: TimeDimID, TimeDimLength
character(len=512) :: msgstring

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

! Check that the input file exists ...

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'restart_file_to_sv:',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncid), &
             'restart_file_to_sv','open '//trim(filename))

unlimitedDimID = -1
iostatus = nf90_Inquire(ncid, nDimensions, nVariables, nAttributes, unlimitedDimID)
call nc_check(iostatus, 'restart_file_to_sv', 'inquire '//trim(filename))

! Start counting and filling the state vector one item at a time,
! repacking the Nd arrays into a single 1d list of numbers.
!              (4dimension)            (1dimension)

! The DART prognostic variables are only defined for a single time.
! We already checked the assumption that variables are xy2d or xyz3d ...
! IF the netCDF variable has an 'ocean_time' dimension,
! we need to read the LAST timestep  ...

TimeDimID = find_time_dimension(ncid, filename, last_file_time)

if (TimeDimID /= unlimitedDimID) then
   write(string1,*)'Time dimension is not the unlimited dimension in '//trim(filename)
   write(string2,*)'Time      dimension is dimension ID ',TimeDimID
   write(string3,*)'Unlimited dimension is dimension ID ',unlimitedDimID
   call error_handler(E_ERR, 'restart_file_to_sv:', string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=TimeDimLength), &
         'restart_file_to_sv', 'inquire timedimlength '//trim(filename))

do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   msgstring = trim(filename)//' '//trim(varname)

   ! determine the shape of the netCDF variable

   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            'restart_file_to_sv', 'inq_varid '//trim(msgstring))

   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
            'restart_file_to_sv', 'inquire '//trim(msgstring))

   mystart = 1   ! These are arrays, actually.
   mycount = 1

   ! Only checking the shape of the variable

   DimCheck : do indx = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') indx,trim(msgstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(indx), len=dimlen), &
            'restart_file_to_sv', string1)

      if ( dimlen /= progvar(ivar)%dimlens(indx) ) then
         write(string1,*) trim(msgstring),' dim/dimlen ',&
                          indx, dimlen,' not ',progvar(ivar)%dimlens(indx)
         call error_handler(E_ERR,'restart_file_to_sv:',string1,source,revision,revdate)
      endif

      mycount(indx) = dimlen

   enddo DimCheck

   where(dimIDs == TimeDimID) mystart = TimeDimLength  ! pick the latest time
   where(dimIDs == TimeDimID) mycount = 1              ! only use one time

   if (debug > 2) then
      write(string2,*)'netCDF "start" is ',mystart(1:ncNdims)
      write(string3,*)'netCDF "count" is ',mycount(1:ncNdims)
      call error_handler(E_MSG, 'restart_file_to_sv:', msgstring, &
                        text2=string2, text3=string3)
   endif

   if (ncNdims == 1) then

      ndim1 = mycount(1)
      allocate(data_1d_array(ndim1))

      call nc_check(nf90_get_var(ncid, VarID, data_1d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'restart_file_to_sv', 'get_var '//trim(varname))

      call prog_var_to_vector(data_1d_array, state_vector, ivar)

      deallocate(data_1d_array)

   elseif (ncNdims == 2) then

      ndim1 = mycount(1)
      ndim2 = mycount(2)
      allocate(data_2d_array(ndim1, ndim2))

      call nc_check(nf90_get_var(ncid, VarID, data_2d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'restart_file_to_sv', 'get_var '//trim(varname))
      call prog_var_to_vector(data_2d_array, state_vector, ivar)

      deallocate(data_2d_array)

   elseif (ncNdims == 3) then

      ndim1 = mycount(1)
      ndim2 = mycount(2)
      ndim3 = mycount(3)
      allocate(data_3d_array(ndim1, ndim2, ndim3))

      call nc_check(nf90_get_var(ncid, VarID, data_3d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'restart_file_to_sv', 'get_var '//trim(varname))
      call prog_var_to_vector(data_3d_array, state_vector, ivar)

      deallocate(data_3d_array)

   elseif (ncNdims == 4) then

      ndim1 = mycount(1)
      ndim2 = mycount(2)
      ndim3 = mycount(3)
      ndim4 = mycount(4)
      allocate(data_4d_array(ndim1, ndim2, ndim3, ndim4))

      call nc_check(nf90_get_var(ncid, VarID, data_4d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'restart_file_to_sv', 'get_var '//trim(varname))
      call prog_var_to_vector(data_4d_array, state_vector, ivar)

      deallocate(data_4d_array)

   else
      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'restart_file_to_sv:', string1, &
                        source,revision,revdate)
   endif

enddo

call nc_check(nf90_close(ncid), &
             'restart_file_to_sv','close '//trim(filename))

end subroutine restart_file_to_sv


!-----------------------------------------------------------------------
!>
!> Writes the (posterior) DART state into a ROMS netCDF analysis file.
!> Just to be clear: only the ROMS variables that are part of the DART
!> state are overwritten.
!>
!> @param state_vector the DART posterior state.
!> @param filename the ROMS netCDF analysis file.
!> @param statetime the DART time of the posterior.
!>
!> @todo FIXME ... use 'statetime' to find the proper time slot in the
!>       restart file. At present, we just stuff it into the last slot,
!>       regardless.
!>

subroutine sv_to_restart_file(state_vector, filename, statetime)

real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: filename
type(time_type),  intent(in) :: statetime

! temp space to hold data while we are writing it
integer :: i, ivar
real(r8), allocatable, dimension(:)     :: data_1d_array
real(r8), allocatable, dimension(:,:)   :: data_2d_array
real(r8), allocatable, dimension(:,:,:) :: data_3d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME) :: varname
integer :: VarID, ncNdims, dimlen
integer :: ncFileID, TimeDimID, TimeDimLength
character(len=512) :: msgstring

if ( .not. module_initialized ) call static_init_model

! Check that the output file exists ...

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for writing.'
   call error_handler(E_ERR,'sv_to_restart_file:',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_WRITE, ncFileID), &
             'sv_to_restart_file','open '//trim(filename))

! make sure the time in the file is the same as the time on the data
! we are trying to insert.  we are only updating part of the contents
! of the analysis file, and state vector contents from a different
! time won't be consistent with the rest of the file.

! The DART prognostic variables are only defined for a single time.
! We already checked the assumption that variables are xy2d or xyz3d ...
! IF the netCDF variable has a TIME dimension, it must be the last dimension,
! and we need to read the LAST timestep and effectively squeeze out the
! singleton dimension when we stuff it into the DART state vector.

TimeDimID = find_time_dimension( ncFileID, filename )

call nc_check(nf90_inquire_dimension(ncFileID, TimeDimID, len=TimeDimLength), &
              'sv_to_restart_file', 'inquire timedimlength '//trim(filename))

PROGVARLOOP : do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   msgstring = trim(filename)//' '//trim(varname)

   ! Ensure netCDF variable is conformable with progvar quantity.
   ! The TIME and Copy dimensions are intentionally not queried
   ! by looping over the dimensions stored in the progvar type.

   call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'sv_to_restart_file', 'inq_varid '//trim(msgstring))

   call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
            'sv_to_restart_file', 'inquire '//trim(msgstring))

   mystart = 1   ! These are arrays, actually.
   mycount = 1
   DimCheck : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(msgstring)
      call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
            'sv_to_restart_file', string1)

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(msgstring),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         write(string2,*)' but it should be.'
         call error_handler(E_ERR, 'sv_to_restart_file:', string1, &
                         source, revision, revdate, text2=string2)
      endif

      mycount(i) = dimlen

   enddo DimCheck

   ! FIXME ... this is where statetime would be useful.
   where(dimIDs == TimeDimID) mystart = TimeDimLength
   where(dimIDs == TimeDimID) mycount = 1   ! only the latest timestep

   if (debug > 2) then
      write(string2,*)'netCDF "start" is ',mystart(1:ncNdims)
      write(string3,*)'netCDF "count" is ',mycount(1:ncNdims)
      call error_handler(E_MSG, 'sv_to_restart_file:', msgstring, &
                        text2=string2, text3=string3)

      write(string2,*)'         shape is ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      call error_handler(E_MSG, 'sv_to_restart_file:', msgstring, &
                        text2=string2)
   endif

   if (progvar(ivar)%numdims == 1) then
      allocate(data_1d_array(mycount(1)))

      call vector_to_prog_var(state_vector, ivar, data_1d_array)

      ! did the user specify lower and/or upper bounds for this variable?
      ! if so, follow the instructions to either fail on out-of-range values,
      ! or set out-of-range values to the given min or max vals

      call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array, &
            start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'sv_to_restart_file', 'put_var '//trim(varname))
      deallocate(data_1d_array)

   elseif (progvar(ivar)%numdims == 2) then

      allocate(data_2d_array(mycount(1), mycount(2)))

      call vector_to_prog_var(state_vector, ivar, data_2d_array)

      ! did the user specify lower and/or upper bounds for this variable?
      ! if so, follow the instructions to either fail on out-of-range values,
      ! or set out-of-range values to the given min or max vals
      if ( progvar(ivar)%clamping ) then
         call do_clamping(progvar(ivar)%out_of_range_fail, progvar(ivar)%range, &
                          progvar(ivar)%numdims, varname, array_2d = data_2d_array)
      endif

      call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
            start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'sv_to_restart_file', 'put_var '//trim(varname))
      deallocate(data_2d_array)

   elseif (progvar(ivar)%numdims == 3) then

      allocate(data_3d_array(mycount(1), mycount(2), mycount(3)))

      call vector_to_prog_var(state_vector, ivar, data_3d_array)

      ! did the user specify lower and/or upper bounds for this variable?
      ! if so, follow the instructions to either fail on out-of-range values,
      ! or set out-of-range values to the given min or max vals
      if ( progvar(ivar)%clamping ) then
         call do_clamping(progvar(ivar)%out_of_range_fail, progvar(ivar)%range, &
                          progvar(ivar)%numdims, varname, array_3d = data_3d_array)
      endif

      call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
            start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'sv_to_restart_file', 'put_var '//trim(varname))
      deallocate(data_3d_array)

   else
      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'statevector_to_analysis_file:', string1, &
                        source,revision,revdate)
   endif

enddo PROGVARLOOP

call nc_check(nf90_close(ncFileID), &
             'sv_to_restart_file','close '//trim(filename))

end subroutine sv_to_restart_file


!-----------------------------------------------------------------------
!>
!> return the name of the ROMS analysis filename that was set
!> in the model_nml namelist
!>

subroutine get_model_restart_filename( filename )

character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(model_restart_filename)

end subroutine get_model_restart_filename


!-----------------------------------------------------------------------
!>
!> parse the model_nml "analysis_time" string to extract a DART time.
!> the expected format of the string is YYYY-MM-DD HH:MM:SS
!> The time is expected to come from the Gregorian calendar.

function get_time_from_namelist()

type(time_type) :: get_time_from_namelist

integer :: iyear, imonth, iday, ihour, imin, isec
integer :: io

if ( .not. module_initialized ) call static_init_model

read(analysis_time,'(i4,5(1x,i2))',iostat=io) iyear, imonth, iday, ihour, imin, isec
if (io /= 0) then
   write(string1, *)'unable to parse model_mod_nml:analysis_time string'
   write(string2, *)'string is >'//trim(analysis_time)//'<'
   write(string3, *)'format is >YYYY-MM-DD HH:MM:SS<'
   call error_handler(E_ERR,'get_time_from_namelist:', string1, &
                     source, revision, revdate, text2=string2, text3=string3)
endif

get_time_from_namelist = set_date(iyear, imonth, iday, ihour, imin, isec)

end function get_time_from_namelist


!-----------------------------------------------------------------------
!>
!> given a DART state vector; print out the min and max
!> data values for all the variables in the vector.
!>
!> @param x the vector.
!>

subroutine print_variable_ranges(x)

real(r8), intent(in) :: x(:)

integer :: ivar

if ( .not. module_initialized ) call static_init_model

do ivar = 1, nfields
   call print_minmax(ivar, x)
enddo

end subroutine print_variable_ranges


!-----------------------------------------------------------------------
! The remaining (private) interfaces come last.
! None of the private interfaces need to call static_init_model()
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!>
!> Set the desired minimum model advance time. This is generally NOT the
!> dynamical timestep of the model, but rather the shortest forecast length
!> you are willing to make. This impacts how frequently the observations
!> may be assimilated.
!>

function set_model_time_step()

type(time_type) :: set_model_time_step

! assimilation_period_seconds, assimilation_period_days are from the namelist

!> @todo FIXME make sure set_model_time_step is an integer multiple of
!> the dynamical timestep or whatever strategy ROMS employs.

set_model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end function set_model_time_step


!-----------------------------------------------------------------------
!>
!> Read the grid dimensions from the ROMS grid netcdf file.
!> By reading the dimensions first, we can use them in variable
!> declarations later - which is faster than using allocatable arrays.
!>

subroutine get_grid_dimensions()

integer  :: ncid

call nc_check(nf90_open(trim(grid_definition_filename), nf90_nowrite, ncid), &
              'get_grid_dimensions', 'open '//trim(grid_definition_filename))

Nxi_rho   = get_dimension_length(ncid, 'xi_rho',   grid_definition_filename)
Nxi_u     = get_dimension_length(ncid, 'xi_u',     grid_definition_filename)
Nxi_v     = get_dimension_length(ncid, 'xi_v',     grid_definition_filename)
Nxi_psi   = get_dimension_length(ncid, 'xi_psi',   grid_definition_filename)
Neta_rho  = get_dimension_length(ncid, 'eta_rho',  grid_definition_filename)
Neta_u    = get_dimension_length(ncid, 'eta_u',    grid_definition_filename)
Neta_v    = get_dimension_length(ncid, 'eta_v',    grid_definition_filename)
Neta_psi  = get_dimension_length(ncid, 'eta_psi',  grid_definition_filename)
Nxi_vert  = get_dimension_length(ncid, 'xi_vert',  grid_definition_filename)
Neta_vert = get_dimension_length(ncid, 'eta_vert', grid_definition_filename)
Ns_rho    = get_dimension_length(ncid, 's_rho',    grid_definition_filename)
Ns_w      = get_dimension_length(ncid, 's_w',      grid_definition_filename)

Nx =  Nxi_rho  ! Setting the nominal value of the 'global' variables
Ny = Neta_rho  ! Setting the nominal value of the 'global' variables
Nz =   Ns_rho  ! Setting the nominal value of the 'global' variables

call nc_check(nf90_close(ncid), &
              'get_grid_dimensions','close '//trim(grid_definition_filename))

end subroutine get_grid_dimensions


!-----------------------------------------------------------------------
!>
!> Read the actual grid values from the ROMS netcdf file.
!>
!> @todo FIXME If the grid variables do not exist, they are calculated.
!> It is not clear to me if this code block is correct in all circumstances,
!> as I do not have access to a good test file. Minimally, instead of requiring
!> setting a hardcoded logical, a simple test to see if the required variables
!> exist and then taking appropriate action would be greatly preferred.

subroutine get_grid()

integer  :: k,ncid, VarID, stat

! The following 'automatic' arrays are more efficient than allocatable arrays.
! This is, in part, why the grid dimensions were determined previously.

real(r8) :: dzt0(Nx,Ny)
real(r8) :: s_rho(Ns_rho), Cs_r(Ns_rho), SSH(Nx,Ny)

real(r8), parameter :: all_land = 0.001_r8

if (debug > 1) then
   write(string1,*)'..  NX is ',Nx
   write(string2,*)'NY is ',Ny
   write(string3,*)'NZ is ',Nz
   call error_handler(E_MSG,'get_grid:',string1,text2=string2,text3=string3)
endif

! Open the netcdf file data

call nc_check(nf90_open(trim(grid_definition_filename), nf90_nowrite, ncid), &
      'get_grid', 'open '//trim(grid_definition_filename))

if (.not. allocated(ULAT)) allocate(ULAT(Nxi_u,Neta_u))
if (.not. allocated(ULON)) allocate(ULON(Nxi_u,Neta_u))
if (.not. allocated(VLAT)) allocate(VLAT(Nxi_v,Neta_v))
if (.not. allocated(VLON)) allocate(VLON(Nxi_v,Neta_v))
if (.not. allocated(TLAT)) allocate(TLAT(Nxi_rho,Neta_rho))
if (.not. allocated(TLON)) allocate(TLON(Nxi_rho,Neta_rho))
if (.not. allocated(HT))   allocate(  HT(Nxi_rho,Neta_rho))
if (.not. allocated(ANGL)) allocate(ANGL(Nxi_rho,Neta_rho))

if (.not. allocated(ZC))   allocate(  ZC(Nx,Ny,Nz))

!> todo are PM, PN used for anything?

if (.not. allocated(PM))   allocate(  PM(Nx,Ny))
if (.not. allocated(PN))   allocate(  PN(Nx,Ny))

! If there is no lat and lon information in the file JUST DIE

call nc_check(nf90_inq_varid(ncid, 'lon_rho', VarID), &
   'get_grid', 'inq_varid lon_rho '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, TLON), &
      'get_grid', 'get_var lon_rho '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'lat_rho', VarID), &
      'get_grid', 'inq_varid lat_rho '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, TLAT), &
      'get_grid', 'get_var lat_rho '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'lon_u', VarID), &
      'get_grid', 'inq_varid lon_u '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, ULON), &
      'get_grid', 'get_var lon_u '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'lat_u', VarID), &
      'get_grid', 'inq_varid lat_u '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, ULAT), &
      'get_grid', 'get_var lat_u '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'lon_v', VarID), &
      'get_grid', 'inq_varid lon_v '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, VLON), &
      'get_grid', 'get_var lon_v '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'lat_v', VarID), &
      'get_grid', 'inq_varid lat_v '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, VLAT), &
      'get_grid', 'get_var lat_v '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'h', VarID), &
      'get_grid', 'inq_varid h '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, HT), &
      'get_grid', 'get_var h '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'pm', VarID), &
      'get_grid', 'inq_varid pm '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, PM), &
      'get_grid', 'get_var pm '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'pn', VarID), &
      'get_grid', 'inq_varid pn '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, PN), &
      'get_grid', 'get_var pn '//trim(grid_definition_filename))

stat=nf90_inq_varid(ncid, 'angle', VarID)

if (stat /= nf90_noerr) then
   ANGL(:,:)=0.0
else
   call nc_check(nf90_get_var( ncid, VarID, ANGL), &
         'get_grid', 'get_var angle '//trim(grid_definition_filename))
endif

call nc_check(nf90_inq_varid(ncid, 's_rho', VarID), &
      'get_grid', 'inq_varid s_rho '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, s_rho), &
      'get_grid', 'get_var s_rho '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'Cs_r', VarID), &
      'get_grid', 'inq_varid Cs_r '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, Cs_r), &
      'get_grid', 'get_var Cs_r '//trim(grid_definition_filename))

! Allocate the masks and set them all to water.
! 1 is water, 0 is land (apparently) There are seemingly no intermediate values.
! In the ROMS netCDF file, these are typed 'double' (overkill),
! but we are implementing as short integer.
! (Could even do logicals if netCDF (output) supported them.)

if (.not. allocated(mask_rho)) allocate(mask_rho(Nxi_rho, Neta_rho))
if (.not. allocated(mask_psi)) allocate(mask_psi(Nxi_psi, Neta_psi))
if (.not. allocated(mask_u))   allocate(mask_u(  Nxi_u  , Neta_u  ))
if (.not. allocated(mask_v))   allocate(mask_v(  Nxi_v  , Neta_v  ))

mask_rho = 0
mask_psi = 0
mask_u   = 0
mask_v   = 0

! Read mask on RHO-points

call nc_check(nf90_inq_varid(ncid, 'mask_rho', VarID), &
      'get_grid', 'inq_varid mask_rho '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, mask_rho), &
      'get_grid', 'get_var mask_rho '//trim(grid_definition_filename))
where(mask_rho > all_land) mask_rho = 1

! Read mask on PSI-points

call nc_check(nf90_inq_varid(ncid, 'mask_psi', VarID), &
      'get_grid', 'inq_varid mask_psi '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, mask_psi), &
      'get_grid', 'get_var mask_psi '//trim(grid_definition_filename))
where(mask_psi > all_land) mask_psi = 1

! Read mask on U-points

call nc_check(nf90_inq_varid(ncid, 'mask_u', VarID), &
      'get_grid', 'inq_varid mask_u '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, mask_u), &
      'get_grid', 'get_var mask_u '//trim(grid_definition_filename))
where(mask_u > all_land) mask_u = 1

! Read mask on V-points

call nc_check(nf90_inq_varid(ncid, 'mask_v', VarID), &
      'get_grid', 'inq_varid mask_v '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, mask_v), &
      'get_grid', 'get_var mask_v '//trim(grid_definition_filename))
where(mask_v > all_land) mask_v = 1

call nc_check(nf90_close(ncid), &
             'get_var','close '//trim(grid_definition_filename))

!     s_rho=(/-0.975, -0.925, -0.875, -0.825, -0.775, -0.725, -0.675, &
!             -0.625, -0.575, -0.525, -0.475, -0.425, -0.375, -0.325, &
!             -0.275, -0.225, -0.175, -0.125, -0.075, -0.025/)
!     Cs_r=(/-0.882485522505154, -0.687254423585503, -0.535200908367393, &
!            -0.416772032329648, -0.324527359249072, -0.252671506867986, &
!            -0.196690045050832, -0.153065871294677, -0.11905824454479,  &
!            -0.0925305948496471, -0.0718162907903607, -0.0556139313622304, &
!            -0.042905583895255, -0.0328928312128658, -0.0249466101148999, &
!            -0.0185676897273263, -0.0133553067236276, -0.00898198688799211, &
!            -0.00517297115481568, -0.00168895354075999/)


! Also get zeta values in order to calculate ZC

call nc_check(nf90_open(trim(model_restart_filename), nf90_nowrite, ncid), &
      'get_grid', 'open '//trim(model_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'zeta', VarID), &
      'get_grid', 'inq_varid zeta '//trim(model_restart_filename))
call nc_check(nf90_get_var( ncid, VarID, SSH), &
      'get_grid', 'get_var zeta '//trim(model_restart_filename))

call nc_check(nf90_close(ncid), &
             'get_var','close '//trim(grid_definition_filename))

!this depends on the transform in ROMS, try to modify it if you are using a different one

do k=1,Nz
     dzt0(:,:) = (s_rho(k)-Cs_r(k))*hc + Cs_r(k) * HT(:,:)
     ZC(:,:,k) = -( dzt0(:,:) + SSH(:,:)*(1.0_r8 + dzt0(:,:)/HT(:,:)) )
end do

if (do_output() .and. debug > 0) then
    write(string1,*)'    min/max ULON ',minval(ULON), maxval(ULON)
    write(string2,*)    'min/max ULAT ',minval(ULAT), maxval(ULAT)
    call error_handler(E_MSG,'get_grid',string1, text2=string2)

    write(string1,*)'    min/max VLON ',minval(VLON), maxval(VLON)
    write(string2,*)    'min/max VLAT ',minval(VLAT), maxval(VLAT)
    call error_handler(E_MSG,'get_grid',string1, text2=string2)

    write(string1,*)'    min/max TLON ',minval(TLON), maxval(TLON)
    write(string2,*)    'min/max TLAT ',minval(TLAT), maxval(TLAT)
    call error_handler(E_MSG,'get_grid',string1, text2=string2)
endif

end subroutine get_grid


!-----------------------------------------------------------------------
!>
!> Determines if the variables and DART KINDs specified by the 'variables'
!> namelist variable are available in the ROMS analysis file and that the
!> corresponding DART KIND is also supported.
!>
!>@param state_variables the list of variables and kinds from model_mod_nml
!>@param ncid the netCDF handle to the ROMS analysis file
!>@param filename the name of the ROMS analysis file (for error messages, mostly)
!>@param ngood the number of variable/KIND pairs specified
!>@param table a more convenient form for the variable/KIND pairs. Each row
!>       is a pair, column 1 is the variable, column 2 is the DART KIND

subroutine verify_variables( state_variables, ncid, filename, ngood, table )

character(len=*), intent(in)  :: state_variables(:)
integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: filename
integer,          intent(out) :: ngood
character(len=*), intent(out) :: table(:,:)

integer :: nrows, ncols, i, VarID
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr
logical :: failure

failure = .FALSE. ! perhaps all with go well

nrows = size(table,1)
ncols = size(table,2)

ngood = 0
MyLoop : do i = 1, nrows

   varname    = trim(state_variables(2*i -1))
   dartstr    = trim(state_variables(2*i   ))
   table(i,1) = trim(varname)
   table(i,2) = trim(dartstr)
   if ( table(i,1) == ' ' .and. table(i,2) == ' ' ) exit MyLoop ! Found end of list.

   if ( table(i,1) == ' ' .or. table(i,2) == ' ' ) then
      string1 = 'model_nml:model "variables" not fully specified'
      call error_handler(E_ERR,'verify_variables:',string1,source,revision,revdate)
   endif

   ! Make sure variable exists in model analysis variable list

   write(string1,'(''variable '',a,'' in '',a)') trim(varname), trim(filename)
   write(string2,'(''there is no '',a)') trim(string1)
   call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
                 'verify_variables', trim(string2))

   ! Make sure DART kind is valid

   if( get_raw_obs_kind_index(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_variables:',string1,source,revision,revdate)
   endif

   ngood = ngood + 1
enddo MyLoop


if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'verify_variables:',string1,source,revision,revdate,text2=string2)
endif


end subroutine verify_variables


!-----------------------------------------------------------------------
!>
!> matches variable name in bounds table to assign
!> the bounds if they exist.  otherwise sets the bounds
!> to missing_r8 which means they are unbounded.
!> @todo FIXME the roms_state_bounds module variable that specifies the
!>       bounds is not actually set nor checked nor ... so all variables
!>       are unbounded.
!>
!> @param bounds_table table specifying the numeric bounds for each variable.
!>                     as currently envisioned, a table specifying the variable
!>                     name and its boundaries must be supplied ... ?namelist?
!> @param ivar the handle to the variable in question
!>

subroutine get_variable_bounds(bounds_table, ivar)

character(len=*), intent(in) :: bounds_table(num_bounds_table_columns, max_state_variables)
integer,          intent(in) :: ivar

! local variables
character(len=50) :: bounds_varname, bound
character(len=10) :: clamp_or_fail
real(r8)          :: lower_bound, upper_bound
integer           :: n

write(string1,*)'WARNING routine not tested.'
write(string2,*)'"roms_state_bounds" from namelist is not set, not tested.'
call error_handler(E_MSG,'get_variable_bounds:',string1, text2=string2)

n = 1
do while ( trim(bounds_table(1,n)) /= 'NULL' .and. trim(bounds_table(1,n)) /= '' )

   bounds_varname = trim(bounds_table(1,n))

   if ( bounds_varname == trim(progvar(ivar)%varname) ) then

        bound = trim(bounds_table(2,n))
        if ( bound /= 'NULL' .and. bound /= '' ) then
             read(bound,'(d16.8)') lower_bound
        else
             lower_bound = missing_r8
        endif

        bound = trim(bounds_table(3,n))
        if ( bound /= 'NULL' .and. bound /= '' ) then
             read(bound,'(d16.8)') upper_bound
        else
             upper_bound = missing_r8
        endif

        ! How do we want to handle out of range values?
        ! Set them to predefined limits (clamp) or simply fail (fail).
        clamp_or_fail = trim(bounds_table(4,n))
        if ( clamp_or_fail == 'NULL' .or. clamp_or_fail == '') then
             write(string1, *) 'instructions for CLAMP_or_FAIL on ', &
                                trim(bounds_varname), ' are required'
             call error_handler(E_ERR,'get_variable_bounds:',string1, &
                                source,revision,revdate)
        else if ( clamp_or_fail == 'CLAMP' ) then
             progvar(ivar)%out_of_range_fail = .FALSE.
        else if ( clamp_or_fail == 'FAIL' ) then
             progvar(ivar)%out_of_range_fail = .TRUE.
        else
             write(string1, *) 'last column must be "CLAMP" or "FAIL" for ', &
                  trim(bounds_varname)
             call error_handler(E_ERR,'get_variable_bounds:',string1, &
                  source,revision,revdate, text2='found '//trim(clamp_or_fail))
        endif

        ! Assign the clamping information into the variable
        progvar(ivar)%clamping = .true.
        progvar(ivar)%range    = (/ lower_bound, upper_bound /)

        return
   endif

   n = n + 1

enddo !n

! we got through all the entries in the bounds table and did not
! find any instructions for this variable.  set the values to indicate
! we are not doing any processing when we write the updated state values
! back to the model restart file.

progvar(ivar)%clamping = .false.
progvar(ivar)%range = missing_r8
progvar(ivar)%out_of_range_fail = .false.  ! should be unused so setting shouldn't matter

return

end subroutine get_variable_bounds


!-----------------------------------------------------------------------
!>
!> Pack the values from a 1d array into the DART array
!>
!> @param data_1d_array the array containing the ROMS variable
!> @param x the array containing the DART state
!> @param ivar handle to the DART structure relating what variable is where.

subroutine prog_var_1d_to_vector(data_1d_array, x, ivar)

real(r8), intent(in)    :: data_1d_array(:)
real(r8), intent(inout) :: x(:)
integer,  intent(in)    :: ivar

integer :: idim1,ii

ii = progvar(ivar)%index1

do idim1 = 1, size(data_1d_array, 1)
   x(ii) = data_1d_array(idim1)
   ii = ii + 1
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' read wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'prog_var_1d_to_vector:', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine prog_var_1d_to_vector


!-----------------------------------------------------------------------
!>
!> Pack the values from a 2D array into the 1D DART array
!>
!> @param data_2d_array the array containing the ROMS variable
!> @param x the array containing the DART state
!> @param ivar handle to the DART structure relating what variable is where.
!>

subroutine prog_var_2d_to_vector(data_2d_array, x, ivar)

real(r8), intent(in)    :: data_2d_array(:,:)
real(r8), intent(inout) :: x(:)
integer,  intent(in)    :: ivar

integer :: idim1,idim2,ii

ii = progvar(ivar)%index1

do idim2 = 1,size(data_2d_array, 2)
   do idim1 = 1,size(data_2d_array, 1)
      x(ii) = data_2d_array(idim1,idim2)
      ii = ii + 1
   enddo
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' read wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'prog_var_2d_to_vector:', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine prog_var_2d_to_vector


!-----------------------------------------------------------------------
!>
!> Pack the values from a 3D array into the 1D DART array
!>
!> @param data_3d_array the array containing the ROMS variable
!> @param x the array containing the DART state
!> @param ivar handle to the DART structure relating what variable is where.
!>

subroutine prog_var_3d_to_vector(data_3d_array, x, ivar)

real(r8), intent(in)    :: data_3d_array(:,:,:)
real(r8), intent(inout) :: x(:)
integer,  intent(in)    :: ivar

integer :: idim1,idim2,idim3,ii

ii = progvar(ivar)%index1

do idim3 = 1,size(data_3d_array, 3)
   do idim2 = 1,size(data_3d_array, 2)
      do idim1 = 1,size(data_3d_array, 1)
         x(ii) = data_3d_array(idim1,idim2,idim3)
         ii = ii + 1
      enddo
   enddo
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' read wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'prog_var_3d_to_vector:', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine prog_var_3d_to_vector


!-----------------------------------------------------------------------
!>
!> Pack the values from a 4D array into the 1D DART array
!>
!> @param data_4d_array the array containing the ROMS variable
!> @param x the array containing the DART state
!> @param ivar handle to the DART structure relating what variable is where.
!>

subroutine prog_var_4d_to_vector(data_4d_array, x, ivar)

real(r8), intent(in)    :: data_4d_array(:,:,:,:)
real(r8), intent(inout) :: x(:)
integer,  intent(in)    :: ivar

integer :: idim1,idim2,idim3,idim4,ii

ii = progvar(ivar)%index1

do idim4 = 1,size(data_4d_array, 4)
   do idim3 = 1,size(data_4d_array, 3)
      do idim2 = 1,size(data_4d_array, 2)
         do idim1 = 1,size(data_4d_array, 1)
            x(ii) = data_4d_array(idim1,idim2,idim3,idim4)
            ii = ii + 1
         enddo
      enddo
   enddo
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' read wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'prog_var_4d_to_vector:', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine prog_var_4d_to_vector


!-----------------------------------------------------------------------
!>
!> Extract the values from the 1D DART array into a 1D array
!>
!> @param x the array containing the DART state
!> @param ivar handle to the DART structure relating what variable is where.
!> @param data_1d_array the array containing the ROMS variable
!>

subroutine vector_to_1d_prog_var(x, ivar, data_1d_array)

real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: ivar
real(r8), intent(out) :: data_1d_array(:)

integer :: idim1,ii

ii = progvar(ivar)%index1

do idim1 = 1, size(data_1d_array, 1)
   data_1d_array(idim1) = x(ii)
   ii = ii + 1
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_1d_prog_var:', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_1d_prog_var


!-----------------------------------------------------------------------
!>
!> Extract the values from the 1D DART array into a 2D array
!>
!> @param x the array containing the DART state
!> @param ivar handle to the DART structure relating what variable is where.
!> @param data_2d_array the array containing the ROMS variable
!>

subroutine vector_to_2d_prog_var(x, ivar, data_2d_array)

real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: ivar
real(r8), intent(out) :: data_2d_array(:,:)

integer :: idim1,idim2,ii

ii = progvar(ivar)%index1

do idim2 = 1,size(data_2d_array, 2)
   do idim1 = 1,size(data_2d_array, 1)
      data_2d_array(idim1,idim2) = x(ii)
      ii = ii + 1
   enddo
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_2d_prog_var:', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_2d_prog_var


!-----------------------------------------------------------------------
!>
!> Extract the values from the 1D DART array into a 3D array
!>
!> @param x the array containing the DART state
!> @param ivar handle to the DART structure relating what variable is where.
!> @param data_3d_array the array containing the ROMS variable
!>

subroutine vector_to_3d_prog_var(x, ivar, data_3d_array)

real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: ivar
real(r8), intent(out) :: data_3d_array(:,:,:)

integer :: idim1,idim2,idim3,ii

ii = progvar(ivar)%index1

do idim3 = 1,size(data_3d_array, 3)
   do idim2 = 1,size(data_3d_array, 2)
      do idim1 = 1,size(data_3d_array, 1)
         data_3d_array(idim1,idim2,idim3) = x(ii)
         ii = ii + 1
      enddo
   enddo
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_3d_prog_var:', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_3d_prog_var


!-----------------------------------------------------------------------
!>
!> given a DART KIND string,
!> return the first and last index into the DART array of that KIND.
!>
!> @param string the variable name
!> @param index1 first index into the DART state
!> @param indexN last index into the DART state
!>

subroutine get_index_range_string(string,index1,indexN)

character(len=*),  intent(in)  :: string
integer,           intent(out) :: index1
integer, optional, intent(out) :: indexN

integer :: indx

index1 = 0
if (present(indexN)) indexN = 0

FieldLoop : do indx=1,nfields
   if (progvar(indx)%kind_string /= trim(string)) cycle FieldLoop
   index1 = progvar(indx)%index1
   if (present(indexN)) indexN = progvar(indx)%indexN
   exit FieldLoop
enddo FieldLoop

if (index1 == 0) then
   write(string1,*) 'Problem, cannot find indices for '//trim(string)
   call error_handler(E_ERR,'get_index_range_string:',string1,source,revision,revdate)
endif

end subroutine get_index_range_string


!-----------------------------------------------------------------------
!>
!> given a DART KIND integer,
!> return the first and last index into the DART array of that KIND.
!>
!> @param dartkind the integer describing the DART KIND
!> @param index1 first index into the DART state
!> @param indexN last index into the DART state
!>

subroutine get_index_range_int(dartkind,index1,indexN)

integer,           intent(in)  :: dartkind
integer,           intent(out) :: index1
integer, optional, intent(out) :: indexN

integer :: indx
character(len=paramname_length) :: string

index1 = 0
if (present(indexN)) indexN = 0

FieldLoop : do indx=1,nfields
   if (progvar(indx)%dart_kind /= dartkind) cycle FieldLoop
   index1 = progvar(indx)%index1
   if (present(indexN)) indexN = progvar(indx)%indexN
   exit FieldLoop
enddo FieldLoop

string = get_raw_obs_kind_name(dartkind)

if (index1 == 0) then
   write(string1,*) 'Problem, cannot find indices for kind ',dartkind,trim(string)
   call error_handler(E_ERR,'get_index_range_int:',string1,source,revision,revdate)
endif

end subroutine get_index_range_int


!-----------------------------------------------------------------------
!>
!> Find the 'ocean_time' dimension in a ROMS netCDF file.
!> If there is none - it is a fatal error.
!> If the optional argument is present, the time of the last timestep
!> is also returned.
!>
!> @param TimeDimID the netCDF dimension ID for 'ocean_time'
!> @param ncid the netCDF handle to the ROMS netCDF file.
!> @param filename the name of the ROMS netCDF file
!>                 (used to generate useful error messages).
!> @param last_time the time/date of the last timestep in the file.
!>
!> @todo FIXME Make sure the calculation is correct.
!>       The metadata claims to have a julian calendar, the 64bit real
!>       can support whole numbers that overflow a 32 bit integer. The
!>       only test file I have is for year 223 or something like that.
!>       DART normally uses a Gregorian calendar since our observations
!>       use that same calendar.

function find_time_dimension(ncid, filename, last_time) result(TimeDimId)

integer                                :: TimeDimId
integer,                   intent(in)  :: ncid
character(len=*),          intent(in)  :: filename
type(time_type), optional, intent(out) :: last_time

integer :: ios, VarID, dimlen
character(len=256) :: unitstring
character(len=256) :: calendarstring

integer :: year, month, day, hour, minute, second
real(digits12), allocatable :: ocean_time(:)
real(digits12)  :: dstart
type(time_type) :: time_offset, base_time, relative_time, spinup_end
integer(i8) :: big_integer
integer :: some_seconds, some_days

call nc_check(nf90_inq_dimid(ncid,'ocean_time',dimid=TimeDimId), &
       'find_time_dimension','cannot find "ocean_time" dimension in '//trim(filename))

! If you want the last time from the 'ocean_time' variable,
! you have to decode it a bit from the netCDF file metadata.

if (present(last_time)) then

   call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=dimlen), &
          'find_time_dimension', 'inquire_dimension ocean_time'//trim(filename))

   call nc_check(nf90_inq_varid(ncid, 'ocean_time', VarID), &
          'find_time_dimension', 'inq_varid ocean_time from '//trim(filename))

   allocate(ocean_time(dimlen))

   call nc_check(nf90_get_var( ncid, VarID, ocean_time), &
          'find_time_dimension', 'get_var ocean_time from '//trim(filename))

   ! Make sure the calendar is expected form
   ! ocean-time:calendar = "julian" ;
   ! ocean-time:units    = "seconds since 0001-01-01 00:00:00" ;
   !                        1234567890123

   call nc_check(nf90_get_att(ncid, VarID, 'units', unitstring), &
          'find_time_dimension', 'get_att ocean_time units '//trim(filename))
   call nc_check(nf90_get_att(ncid, VarID, 'calendar', calendarstring), &
          'find_time_dimension', 'get_att ocean_time units '//trim(filename))

   if (trim(calendarstring) /= 'julian') then
      write(string1,*)'expecting ocean_time calendar of "julian"'
      write(string2,*)'got '//trim(calendarstring)
      call error_handler(E_ERR,'find_time_dimension:', string1, &
                source, revision, revdate, text2=string2, text3=string3)
   endif

   call set_calendar_type( trim(calendarstring) )

   if (unitstring(1:13) /= 'seconds since') then
      write(string1,*)'expecting ocean_time units of "seconds since ..."'
      write(string2,*)'got '//trim(unitstring)
      write(string3,*)'Cannot proceed. Stopping.'
      call error_handler(E_ERR,'find_time_dimension:', string1, &
                source, revision, revdate, text2=string2, text3=string3)
   endif

   read(unitstring,'(14x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
   if (ios /= 0) then
      write(string1,*)'Unable to read ocean_time units. Error status was ',ios
      write(string2,*)'expected "seconds since YYYY-MM-DD HH:MM:SS"'
      write(string3,*)'was      "'//trim(unitstring)//'"'
      call error_handler(E_ERR, 'find_time_dimension:', string1, &
             source, revision, revdate, text2=string2, text3=string3)
   endif

   ! big_integer may overflow a 32bit integer, so declare it 64bit
   ! and parse it into an integer number of days and seconds, both
   ! of which can be 32bit. Our set_time, set_date routines need 32bit.

   big_integer  = int(ocean_time(dimlen),i8)
   some_days    = big_integer / (24*60*60)
   big_integer  = big_integer - some_days * (24*60*60)
   some_seconds = int(big_integer,i4)

   time_offset   = set_time(some_seconds, some_days)
   base_time     = set_date(year, month, day, hour, minute, second)
   relative_time = base_time + time_offset

   ! Must get the time stamp assigned to model intitialization
   ! DART is interpreting this as the number of days of the 'origin'

   call nc_check(nf90_inq_varid(ncid, 'dstart', VarID), &
          'find_time_dimension', 'inq_varid dstart from '//trim(filename))

   call nc_check(nf90_get_var( ncid, VarID, dstart), &
          'find_time_dimension', 'get_var dstart from '//trim(filename))

   call nc_check(nf90_get_att(ncid, VarID, 'units', unitstring), &
          'find_time_dimension', 'get_att dstart units '//trim(filename))
   call nc_check(nf90_get_att(ncid, VarID, 'calendar', calendarstring), &
          'find_time_dimension', 'get_att dstart units '//trim(filename))

   if (unitstring(1:10) /= 'days since') then
      write(string1,*)'expecting dstart units of "days since ..."'
      write(string2,*)'got '//trim(unitstring)
      write(string3,*)'Cannot proceed. Stopping.'
      call error_handler(E_ERR,'find_time_dimension:', string1, &
                source, revision, revdate, text2=string2, text3=string3)
   endif

   if (trim(calendarstring) /= 'julian') then
      write(string1,*)'expecting dstart calendar of "julian"'
      write(string2,*)'got '//trim(calendarstring)
      call error_handler(E_ERR,'find_time_dimension:', string1, &
                source, revision, revdate, text2=string2, text3=string3)
   endif

   read(unitstring,'(11x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
   if (ios /= 0) then
      write(string1,*)'Unable to read dtime units. Error status was ',ios
      write(string2,*)'expected "days since YYYY-MM-DD HH:MM:SS"'
      write(string3,*)'was      "'//trim(unitstring)//'"'
      call error_handler(E_ERR, 'find_time_dimension:', string1, &
             source, revision, revdate, text2=string2, text3=string3)
   endif

   some_days    = floor(dstart)
   some_seconds = int(fraction(dstart) * (24*60*60))

   time_offset  = set_time(some_seconds, some_days)  ! seconds, days
   base_time    = set_date(year, month, day, hour, minute, second)
   spinup_end   = base_time + time_offset

   last_time = spinup_end + relative_time

   if (do_output() .and. debug > 3) then
      call print_time(last_time, str='last roms time is ',iunit=logfileunit)
      call print_time(last_time, str='last roms time is ')
      call print_date(last_time, str='last roms date is ',iunit=logfileunit)
      call print_date(last_time, str='last roms date is ')
   endif

   call set_calendar_type( calendar )  ! set calendar back

   deallocate(ocean_time)
endif

end function find_time_dimension


!-----------------------------------------------------------------------
!>
!> for a given directive and range, do the data clamping for the given
!> input array.  only one of the optional array args should be specified - the
!> one which matches the given dimsize.  this still has replicated sections for
!> each possible dimensionality (which so far is only 1 to 3 - add 4-7 only
!> if needed) but at least it is isolated to this subroutine.
!>
!> @param out_of_range_fail switch to prescribe action. if TRUE and the data
!>                          value is outside the expected range, issue a FATAL
!>                          error. if FALSE, set the offending value to the
!>                          appropriate max or min described by the 'range' variable.
!> @param range the expected min and max describing the allowable range.
!> @param dimsize the rank of the variable to be clamped.
!> @param varname the name of the variable to be clamped.
!> @param array_1d if dimsize == 1, this contains the data to be clamped.
!> @param array_2d if dimsize == 1, this contains the data to be clamped.
!> @param array_3d if dimsize == 3, this contains the data to be clamped.
!>

subroutine do_clamping(out_of_range_fail, range, dimsize, varname, &
                       array_1d, array_2d, array_3d)
logical,                      intent(in)    :: out_of_range_fail
real(r8),                     intent(in)    :: range(2)
integer,                      intent(in)    :: dimsize
character(len=NF90_MAX_NAME), intent(in)    :: varname
real(r8), optional,           intent(inout) :: array_1d(:)
real(r8), optional,           intent(inout) :: array_2d(:,:)
real(r8), optional,           intent(inout) :: array_3d(:,:,:)

! these sections should all be identical except for the array_XX specified.
! if anyone can figure out a way to defeat fortran's strong typing for arrays
! so we don't have to replicate each of these sections, i'll buy you a cookie.
! (sorry, you can't suggest using the preprocessor, which is the obvious
! solution.  up to now we have avoided any preprocessed code in the entire
! system.  if we cave at some future point this routine is a prime candidate
! to autogenerate.)

if (dimsize == 1) then
   if (.not. present(array_1d)) then
      call error_handler(E_ERR, 'do_clamping:', 'Internal error.  Should not happen', &
                         source,revision,revdate, text2='array_1d not present for 1d case')
   endif

   ! is lower bound set
   if ( range(1) /= missing_r8 ) then

      if ( out_of_range_fail ) then
         if ( minval(array_1d) < range(1) ) then
            write(string1, *) 'min data val = ', minval(array_1d), &
                              'min data bounds = ', range(1)
            call error_handler(E_ERR, 'statevector_to_analysis_file:', &
                        'Variable '//trim(varname)//' failed lower bounds check.', &
                         source,revision,revdate)
         endif
      else
         where ( array_1d < range(1) ) array_1d = range(1)
      endif

   endif ! min range set

   ! is upper bound set
   if ( range(2) /= missing_r8 ) then

      if ( out_of_range_fail ) then
         if ( maxval(array_1d) > range(2) ) then
            write(string1, *) 'max data val = ', maxval(array_1d), &
                              'max data bounds = ', range(2)
            call error_handler(E_ERR, 'statevector_to_analysis_file:', &
                        'Variable '//trim(varname)//' failed upper bounds check.', &
                         source,revision,revdate, text2=string1)
         endif
      else
         where ( array_1d > range(2) ) array_1d = range(2)
      endif

   endif ! max range set

   write(string1, '(A,A32,2F16.7)') 'BOUND min/max ', trim(varname), &
                      minval(array_1d), maxval(array_1d)
   call error_handler(E_MSG, '', string1, source,revision,revdate)

else if (dimsize == 2) then
   if (.not. present(array_2d)) then
      call error_handler(E_ERR, 'do_clamping:', 'Internal error.  Should not happen', &
                         source,revision,revdate, text2='array_2d not present for 2d case')
   endif

   ! is lower bound set
   if ( range(1) /= missing_r8 ) then

      if ( out_of_range_fail ) then
         if ( minval(array_2d) < range(1) ) then
            write(string1, *) 'min data val = ', minval(array_2d), &
                              'min data bounds = ', range(1)
            call error_handler(E_ERR, 'statevector_to_analysis_file:', &
                        'Variable '//trim(varname)//' failed lower bounds check.', &
                         source,revision,revdate)
         endif
      else
         where ( array_2d < range(1) ) array_2d = range(1)
      endif

   endif ! min range set

   ! is upper bound set
   if ( range(2) /= missing_r8 ) then

      if ( out_of_range_fail ) then
         if ( maxval(array_2d) > range(2) ) then
            write(string1, *) 'max data val = ', maxval(array_2d), &
                              'max data bounds = ', range(2)
            call error_handler(E_ERR, 'statevector_to_analysis_file:', &
                        'Variable '//trim(varname)//' failed upper bounds check.', &
                         source,revision,revdate, text2=string1)
         endif
      else
         where ( array_2d > range(2) ) array_2d = range(2)
      endif

   endif ! max range set

   write(string1, '(A,A32,2F16.7)') 'BOUND min/max ', trim(varname), &
                      minval(array_2d), maxval(array_2d)
   call error_handler(E_MSG, '', string1, source,revision,revdate)

else if (dimsize == 3) then
   if (.not. present(array_3d)) then
      call error_handler(E_ERR, 'do_clamping:', 'Internal error.  Should not happen', &
                         source,revision,revdate, text2='array_3d not present for 3d case')
   endif

   ! is lower bound set
   if ( range(1) /= missing_r8 ) then

      if ( out_of_range_fail ) then
         if ( minval(array_3d) < range(1) ) then
            write(string1, *) 'min data val = ', minval(array_3d), &
                              'min data bounds = ', range(1)
            call error_handler(E_ERR, 'statevector_to_analysis_file:', &
                        'Variable '//trim(varname)//' failed lower bounds check.', &
                         source,revision,revdate)
         endif
      else
         where ( array_3d < range(1) ) array_3d = range(1)
      endif

   endif ! min range set

   ! is upper bound set
   if ( range(2) /= missing_r8 ) then

      if ( out_of_range_fail ) then
         if ( maxval(array_3d) > range(2) ) then
            write(string1, *) 'max data val = ', maxval(array_3d), &
                              'max data bounds = ', range(2)
            call error_handler(E_ERR, 'statevector_to_analysis_file:', &
                        'Variable '//trim(varname)//' failed upper bounds check.', &
                         source,revision,revdate, text2=string1)
         endif
      else
         where ( array_3d > range(2) ) array_3d = range(2)
      endif

   endif ! max range set

   write(string1, '(A,A32,2F16.7)') 'BOUND min/max ', trim(varname), &
                      minval(array_3d), maxval(array_3d)
   call error_handler(E_MSG, '', string1, source,revision,revdate)

else
   write(string1, *) 'dimsize of ', dimsize, ' found where only 1-3 expected'
   call error_handler(E_MSG, 'do_clamping:', 'Internal error, should not happen', &
                      source,revision,revdate, text2=string1)
endif   ! dimsize

end subroutine do_clamping


!-----------------------------------------------------------------------
!>
!> convert DART time type into a character string with the
!> format of YYYYMMDDhh ... or DDhh
!>
!> @param time_to_string the character string containing the time
!> @param t the time
!> @param interval logical flag describing if the time is to be
!>                 interpreted as a calendar date or a time increment.
!>                 If the flag is merely present, the time is to be
!>                 interpreted as an increment and the format is simply
!>                 DDhh. If the flag is not present, the time is a full
!>                 calendar (Gregorian) date and will be renedered with
!>                 the YYYYMMDDhh format.
!>

function time_to_string(t, interval)

character(len=19)              :: time_to_string
type(time_type),   intent(in) :: t
logical, optional, intent(in) :: interval

! local variables

integer :: iyear, imonth, iday, ihour, imin, isec
integer :: ndays, nsecs
logical :: dointerval

if (present(interval)) then
   dointerval = interval
else
   dointerval = .false.
endif

! for interval output, output the number of days, then hours, mins, secs
! for date output, use the calendar routine to get the year/month/day hour:min:sec
if (dointerval) then
   call get_time(t, nsecs, ndays)
   if (ndays > 99) then
      write(string1, *) 'interval number of days is ', ndays
      call error_handler(E_ERR,'time_to_string:', 'interval days cannot be > 99', &
                         source, revision, revdate, text2=string1)
   endif
   ihour = nsecs / 3600
   nsecs = nsecs - (ihour * 3600)
   imin  = nsecs / 60
   nsecs = nsecs - (imin * 60)
   isec  = nsecs
!   write(time_to_string, '(I2.2,3(A1,I2.2))') &
!                        ndays, '_', ihour, ':', imin, ':', isec
   write(time_to_string, '(I2.2,I2.2)') &
                        ndays, ihour
else
   call get_date(t, iyear, imonth, iday, ihour, imin, isec)
!   write(time_to_string, '(I4.4,5(A1,I2.2))') &
!                        iyear, '-', imonth, '-', iday, '_', ihour, ':', imin, ':', isec
   write(time_to_string, '(I4.4,I2.2,I2.2,I2.2)') &
             iyear, imonth, iday, ihour

endif

end function time_to_string


!-----------------------------------------------------------------------
!>
!> given a DART variable index and a state vector; print out the variable
!> name, min, and max data values for that DART variable.
!>

subroutine print_minmax(ivar, x)

integer,  intent(in) :: ivar
real(r8), intent(in) :: x(:)

write(string1, '(A,A32,2ES16.7)') 'data  min/max ', &
           trim(progvar(ivar)%varname), &
           minval(x(progvar(ivar)%index1:progvar(ivar)%indexN)), &
           maxval(x(progvar(ivar)%index1:progvar(ivar)%indexN))

call error_handler(E_MSG, '', string1, source,revision,revdate)

end subroutine print_minmax



!-----------------------------------------------------------------------
!>
!> Determine the DART variable index of a particular DART KIND (integer)
!>
!> @param get_progvar_index_from_kind the DART variable index
!> @param dartkind the DART KIND of interest
!>

function get_progvar_index_from_kind(dartkind)

integer             :: get_progvar_index_from_kind
integer, intent(in) :: dartkind

integer :: i

FieldLoop : do i=1,nfields
   if (progvar(i)%dart_kind /= dartkind) cycle FieldLoop
   get_progvar_index_from_kind = i
   return
enddo FieldLoop

get_progvar_index_from_kind = -1

end function get_progvar_index_from_kind


!-----------------------------------------------------------------------
!>
!> use four corner points to interpolate for one specific layer.
!> This needs to be careful since ROMS uses terrain-following coordinates.
!>
!> @param state_handle handle to the dart state
!> @param lon the longitude of interest
!> @param lat the latitude of interest
!> @param var_type the DART KIND of interest
!> @param height the layer index
!> @param interp_val the estimated value of the DART state at the location
!>          of interest (the interpolated value).
!> @param istatus interpolation status ... 0 == success, /=0 is a failure
!>

subroutine lon_lat_interpolate(state_handle, ens_size, lon, lat, var_type, height, interp_val, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in)  :: ens_size
real(r8),            intent(in)  :: lon
real(r8),            intent(in)  :: lat
integer,             intent(in)  :: var_type
integer,             intent(in)  :: height
real(r8),            intent(out) :: interp_val(ens_size)
integer,             intent(out) :: istatus(ens_size)

!real(r8), intent(in)  :: x(:)

! Local storage
integer  :: lat_bot, lat_top, lon_bot, lon_top
integer  :: x_ind, y_ind,ivar
real(r8) :: p(ens_size,4), x_corners(ens_size,4), y_corners(ens_size,4)

integer(i8) :: state_index
integer     :: var_id

! Succesful return has istatus of 0
istatus    = 99
interp_val = MISSING_R8

ivar = get_progvar_index_from_kind(var_type)

   ! Is this on the U or T grid?
   if(progvar(ivar)%kind_string == 'KIND_U_CURRENT_COMPONENT') then
      call get_reg_box_indices(lon, lat, ULON, ULAT, x_ind, y_ind)
      ! Getting corners for accurate interpolation
      call get_quad_corners(ULON, x_ind, y_ind, x_corners)
      call get_quad_corners(ULAT, x_ind, y_ind, y_corners)
   elseif (progvar(ivar)%kind_string == 'KIND_V_CURRENT_COMPONENT') then
      call get_reg_box_indices(lon, lat, VLON, VLAT, x_ind, y_ind)
      ! Getting corners for accurate interpolation
      call get_quad_corners(VLON, x_ind, y_ind, x_corners)
      call get_quad_corners(VLAT, x_ind, y_ind, y_corners)
   else
      call get_reg_box_indices(lon, lat, TLON, TLAT, x_ind, y_ind)
      ! Getting corners for accurate interpolation
      call get_quad_corners(TLON, x_ind, y_ind, x_corners)
      call get_quad_corners(TLAT, x_ind, y_ind, y_corners)
   endif

lon_bot=x_ind
lat_bot=y_ind

! Find the indices to get the values for interpolating
lat_top = lat_bot + 1
if(lat_top > progvar(ivar)%numeta) then
   istatus = 2
   return
endif

lon_top = lon_bot + 1
if(lon_top > progvar(ivar)%numxi) then
   istatus = 2
   return
endif

! Get the values at the four corners of the box or quad
! Corners go around counterclockwise from lower left
! If any one of these fail, go no further.

!>@ todo FIXME : make sure that we are using the mask
!>               properly for locations that have dry land

istatus = 3
var_id  = get_varid_from_kind(domain_id, var_type)

       state_index = get_dart_vector_index(lon_bot, lat_bot, height, domain_id, var_id)
       p(:,1) = get_state(state_index, state_handle)
if(any(p(:,1) == MISSING_R8)) return

       state_index = get_dart_vector_index(lon_top, lat_bot, height, domain_id, var_id)
       p(:,2) = get_state(state_index, state_handle)
if(any(p(:,2) == MISSING_R8)) return

       state_index = get_dart_vector_index(lon_top, lat_top, height, domain_id, var_id)
       p(:,2) = get_state(state_index, state_handle)
if(any(p(:,3) == MISSING_R8)) return

       state_index = get_dart_vector_index(lon_bot, lat_top, height, domain_id, var_id)
       p(:,4) = get_state(state_index, state_handle)
if(any(p(:,4) == MISSING_R8)) return

! Full bilinear interpolation for quads
! istatus = 0 is good, whatever the interp_val is, it is ...

istatus = 0
call quad_bilinear_interp(lon, lat, x_corners, y_corners, p, ens_size, interp_val)

end subroutine lon_lat_interpolate


!-----------------------------------------------------------------------
!>
!> given a state vector and horizontal location, linearly interpolate
!> to some desired vertical location. If the location is deeper than the
!> deepest level, the deepest value is returned. If the location is shallower
!> than the shallowest level, the shallowest level value is returned.
!>
!> @param state_handle handle to the dart state
!> @param ens_size number of ensemble members
!> @param height the vertical location of interest
!> @param lonid the index of the longitude identifying the column.
!> @param latid the index of the latitude identifying the column.
!> @param var_type the DART KIND of interest
!> @param interp_vals the interpolated value at the desired height.
!>

subroutine vert_interpolate(state_handle, ens_size, lonid, latid, var_type, height, interp_vals)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
real(r8),            intent(in)  :: height
integer,             intent(in)  :: lonid
integer,             intent(in)  :: latid
integer,             intent(in)  :: var_type
real(r8),            intent(out) :: interp_vals(ens_size)

! 'depth' or 'height' is positive; thus ROMS vertical depth. i.e., ZC should be positive.

integer    :: i,iidd
real(r8)   :: tp(ens_size,Nz), zz(Nz)

integer(i8) :: state_index
integer     :: var_id

tp(:,:) = 0.0_r8
iidd    = 0

!>@ todo FIXME ... check ... get_state() can return MISSING_R8 ... is this possible

var_id = get_varid_from_kind(domain_id, var_type)

do i=1,Nz
    state_index = get_dart_vector_index(lonid, latid, i, domain_id, var_id)
    tp(:,i)=get_state(state_index, state_handle)
    zz(i)=ZC(lonid,latid,i)
enddo

doloop: do i=Nz,1,-1
     if (zz(i) >= height) then
        iidd=i
        exit doloop
     endif
enddo doloop

if (iidd==0) then
    interp_vals(:)=tp(:,1)
    return
endif

if (iidd==Nz) then
    interp_vals(:)=tp(:,Nz)
    return
endif

interp_vals(:) = tp(:,iidd) +(tp(:,iidd+1) - tp(:,iidd)) * &
   ((height -zz(iidd))/(zz(iidd+1) - zz(iidd)))

!write(*,*) 'model_mod: height, zz(iidd), val ', height, zz(iidd), interp_vals

end subroutine vert_interpolate


!-----------------------------------------------------------------------
!>
!> Given an integer index into the state vector structure, returns the
!> associated array indices for lat, lon, and depth, as well as the type.
!>
!> @param var_type the DART KIND of interest
!> @param offset relative (to the start of the KIND) index into the DART
!>               vector for the KIND of interest.
!> @param x_index the index of the longitude gridcell
!> @param y_index the index of the latitude gridcell
!> @param z_index the index of the vertical gridcell
!>
!> @todo FIXME Check to make sure that this routine is robust for all
!> grid staggers, etc. Seems too simple given the staggers possible.
!>

subroutine get_state_indices(var_type, offset, x_index, y_index, z_index)

integer, intent(in)  :: var_type
integer, intent(in)  :: offset
integer, intent(out) :: x_index
integer, intent(out) :: y_index
integer, intent(out) :: z_index

integer :: ivar, numxi, numeta

ivar   = get_progvar_index_from_kind(var_type)
numxi  = progvar(ivar)%numxi
numeta = progvar(ivar)%numeta

if (progvar(ivar)%kind_string=='KIND_SEA_SURFACE_HEIGHT') then
  z_index = 1
else
  z_index = ceiling( float(offset) / (numxi * numeta))
endif
y_index = ceiling( float(offset - ((z_index-1)*numxi*numeta)) / numxi )
x_index =  offset - ((z_index-1)*numxi*numeta) - ((y_index-1)*numxi)

if (do_output() .and. debug > 3) then
   write(string1,*) 'checking ',trim(progvar(ivar)%varname), ' offset = ', offset
   write(string2,*) 'lon, lat, depth index = ', x_index, y_index, z_index
   call error_handler(E_MSG,'get_state_indices:', string1, &
              source, revision, revdate, text2=string2)
endif

end subroutine get_state_indices


!-----------------------------------------------------------------------
!>
!> Given a longitude and latitude (lon_in, lat), the longitude and
!> latitude of the 4 corners of a quadrilateral and the values at the
!> four corners, interpolates to (lon_in, lat) which is assumed to
!> be in the quad. This is done by bilinear interpolation, fitting
!> a function of the form a + bx + cy + dxy to the four points and
!> then evaluating this function at (lon, lat). The fit is done by
!> solving the 4x4 system of equations for a, b, c, and d. The system
!> is reduced to a 3x3 by eliminating a from the first three equations
!> and then solving the 3x3 before back substituting. There is concern
!> about the numerical stability of this implementation. Implementation
!> checks showed accuracy to seven decimal places on all tests.
!>
!> @param lon_in longitude of interest
!> @param lat latitude of interest
!> @param x_corners_in the longitudes of the 4 surrounding corners
!> @param y_corners the latitudes of the 4 surrounding corners
!> @param p the values at the 4 surrounding corners
!> @param ens_size number of ensemble members
!> @param interp_vals the desired interpolated value
!>

subroutine quad_bilinear_interp(lon_in, lat, x_corners_in, y_corners, &
                                p, ens_size, interp_vals)

real(r8), intent(in)  :: lon_in
real(r8), intent(in)  :: lat
real(r8), intent(in)  :: x_corners_in(4)
real(r8), intent(in)  :: y_corners(4)
real(r8), intent(in)  :: p(ens_size,4)
integer , intent(in)  :: ens_size
real(r8), intent(out) :: interp_vals(ens_size)

! Given a longitude and latitude (lon_in, lat), the longitude and
! latitude of the 4 corners of a quadrilateral and the values at the
! four corners, interpolates to (lon_in, lat) which is assumed to
! be in the quad. This is done by bilinear interpolation, fitting
! a function of the form a + bx + cy + dxy to the four points and
! then evaluating this function at (lon, lat). The fit is done by
! solving the 4x4 system of equations for a, b, c, and d. The system
! is reduced to a 3x3 by eliminating a from the first three equations
! and then solving the 3x3 before back substituting. There is concern
! about the numerical stability of this implementation. Implementation
! checks showed accuracy to seven decimal places on all tests.

integer :: i
real(r8) :: m(3, 3), v(ens_size,3), r(3), a(ens_size), x_corners(4), lon
! real(r8) :: lon_mean

! Watch out for wraparound on x_corners.
lon = lon_in
x_corners = x_corners_in

!*******
! Problems with extremes in polar cell interpolation can be reduced
! by this block, but it is not clear that it is needed for actual
! ocean grid data
! Find the mean longitude of corners and remove
!!!lon_mean = sum(x_corners) / 4.0_r8
!!!x_corners = x_corners - lon_mean
!!!lon = lon - lon_mean
! Multiply everybody by the cos of the latitude
!!!do i = 1, 4
   !!!x_corners(i) = x_corners(i) * cos(y_corners(i) * deg2rad)
!!!enddo
!!!lon = lon * cos(lat * deg2rad)

!*******


! Fit a surface and interpolate; solve for 3x3 matrix
do i = 1, 3
   ! Eliminate a from the first 3 equations
   m(i, 1) = x_corners(i) - x_corners(i + 1)
   m(i, 2) = y_corners(i) - y_corners(i + 1)
   m(i, 3) = x_corners(i)*y_corners(i) - x_corners(i + 1)*y_corners(i + 1)
   v(:,i) = p(:,i) - p(:,i + 1)
enddo

! Solve the matrix for b, c and d
call mat3x3(m, v, r)

! r contains b, c, and d; solve for a
a(:) = p(:,4) - r(1) * x_corners(4) - r(2) * y_corners(4) - &
   r(3) * x_corners(4)*y_corners(4)


!----------------- Implementation test block
! When interpolating on dipole x3 never exceeded 1e-9 error in this test
!!!do i = 1, 4
   !!!interp_vals = a + r(1)*x_corners(i) + r(2)*y_corners(i)+ r(3)*x_corners(i)*y_corners(i)
   !!!if(abs(interp_vals - p(i)) > 1e-9) then
      !!!write(*, *) 'large interp residual ', interp_vals - p(i)
   !!!endif
!!!enddo
!----------------- Implementation test block


! Now do the interpolation
interp_vals(:) = a(:) + r(1)*lon + r(2)*lat + r(3)*lon*lat

!********
! Avoid exceeding maxima or minima as stopgap for poles problem
! When doing bilinear interpolation in quadrangle, can get interpolated
! values that are outside the range of the corner values
do i = 1,ens_size
   if(interp_vals(i) > maxval(p(i,:))) then
      interp_vals(i) = maxval(p(i,:))
   else if(interp_vals(i) < minval(p(i,:))) then
      interp_vals(i) = minval(p(i,:))
   endif
enddo
!********

end subroutine quad_bilinear_interp


!-----------------------------------------------------------------------
!>
!> Solves rank 3 linear system mr = v for r using Cramer's rule.
!>
!> @ todo This isn't the best choice for speed or numerical stability so
!> might want to replace this at some point.
!>

subroutine mat3x3(m, v, r)
real(r8), intent(in)  :: m(3, 3)
real(r8), intent(in)  :: v(3)
real(r8), intent(out) :: r(3)

real(r8) :: m_sub(3, 3), numer, denom
integer  :: i

! Compute the denominator, det(m)
denom = deter3(m)

! Loop to compute the numerator for each component of r
do i = 1, 3
   m_sub = m
   m_sub(:, i) = v
   numer = deter3(m_sub)
   r(i) = numer / denom
enddo

end subroutine mat3x3


!-----------------------------------------------------------------------
!>
!> Computes determinant of 3x3 matrix m
!>

function deter3(m)

real(r8), intent(in) :: m(3, 3)
real(r8)             :: deter3

deter3 = m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) + &
         m(1,3)*m(2,1)*m(3,2) - m(3,1)*m(2,2)*m(1,3) - &
         m(1,1)*m(2,3)*m(3,2) - m(3,3)*m(2,1)*m(1,2)

end function deter3


!-----------------------------------------------------------------------
!>
!> Given a longitude and latitude in degrees returns the index of the
!> regular lon-lat box that contains the point.
!>
!> @param lon the longitude of interest
!> @param lat the latitude of interest
!> @param LON_al the matrix of longitudes defining the boxes (grids).
!> @param LAT_al the matrix of latitudes defining the boxes (grids).
!> @param x_ind the 'x' (longitude) index of the location of interest.
!> @param y_ind the 'y' (latitude) index of the location of interest.
!>
!> @todo FIXME should this really error out or silently fail?
!> At this point it may silently fail.
!>
!> ROMS has a structured horizontal grid.
!>
!> We group the structured grid cells into regular boxes,
!> e.g., 70 x 40 boxes.
!>
!> (1)calculate distances from a given point to the center cell of each box
!> (2)after we find the box which is closest to the given point,
!>    we can calculate distances to the cells within the box.
!> bjchoi 2014/08/07
!-------------------------------------------------------------

!>@ todo FIXME : check whole routine, not sure that is is doing what
!>               we expect...
subroutine get_reg_box_indices(lon, lat, LON_al, LAT_al, x_ind, y_ind)

real(r8), intent(in)  :: lon
real(r8), intent(in)  :: lat
real(r8), intent(in)  :: LON_al(:,:)
real(r8), intent(in)  :: LAT_al(:,:)
integer,  intent(out) :: x_ind
integer,  intent(out) :: y_ind

real(r8),allocatable  :: boxLON(:,:), boxLAT(:,:)
integer               :: i, j, ii, jj, Nx, Ny, nxbox, nybox
integer               :: num_cell, loc_box, i_start,i_end,j_start,j_end
real(r8)              :: lon_dif, lat_dif, tp1, tp2
real(r8),allocatable  :: rtemp(:)
integer, allocatable  :: itemp(:),jtemp(:)
integer               :: t, loc1
integer               :: min_location(1) ! stupid parser

! Divide the model grid cells into a group of boxes.
! For example, divide the model domain into 47 x 70 subdomains (or boxes).
! Each box contains num_cell * num_cell grid cells.
num_cell = 20 ! please use even number,e.g., 4,6,8,10,12..., for num_cell
Nx = size(LON_al,1)
Ny = size(LON_al,2)
nxbox = ceiling( float(Nx) / num_cell )
nybox = ceiling( float(Ny) / num_cell )

! we define the center cell of each box:
! boxLON and boxLAT
allocate( boxLON(nxbox,nybox) )
allocate( boxLAT(nxbox,nybox) )
allocate( rtemp(Nx*Ny) )
allocate( itemp(Nx*Ny) )
allocate( jtemp(Nx*Ny) )

do ii = 1, nxbox-1
  do jj = 1, nybox-1
     i = (num_cell/2) + (ii-1)*num_cell
     j = (num_cell/2) + (jj-1)*num_cell
     boxLON(ii,jj) = LON_al(i,j)
     boxLAT(ii,jj) = LAT_al(i,j)
  enddo
  jj = nybox
  i = (num_cell/2) + (ii-1)*num_cell
  j = (num_cell/2) + (jj-1)*num_cell
  if (j .GT. Ny) j = Ny
  boxLON(ii,jj) = LON_al(i,j)
  boxLAT(ii,jj) = LAT_al(i,j)
enddo

ii = nxbox
i = (num_cell/2) + (ii-1)*num_cell
if (i .GT. Nx) i = Nx
do jj = 1, nybox-1
   j = (num_cell/2) + (jj-1)*num_cell
   boxLON(ii,jj) = LON_al(i,j)
   boxLAT(ii,jj) = LAT_al(i,j)
enddo
jj = nybox
j = (num_cell/2) + (jj-1)*num_cell
if (j .GT. Ny) j = Ny

! find the box which the obs data belongs to
t=0
do ii=1,nxbox
 do jj=1,nybox
    t=t+1
    lon_dif = (lon - boxLON(ii,jj))*111.41288*cos(lat*DEG2RAD) - &
                 0.09350 * cos(3.0 * lat*DEG2RAD) + &
                 0.00012 * cos(5.0 * lat*DEG2RAD)
    lat_dif = (lat - boxLat(ii,jj))*111.13295
    rtemp(t) = sqrt(lon_dif*lon_dif+lat_dif*lat_dif)
    itemp(t) = ii
    jtemp(t) = jj
 enddo
enddo

! find minium distance
min_location = minloc(rtemp)
loc_box = min_location(1)

! locate the box which contains the obs data (ii, jj)
! set the ranges of i and j around the box.
! the neighboring 4 model grid points are within this range.
ii = itemp(loc_box)
jj = jtemp(loc_box)
i_start = max((ii-1)*num_cell - num_cell/2, 1) ! span 2 cells
i_end   = min(    ii*num_cell + num_cell/2,Nx) ! span 2 cells
j_start = max((jj-1)*num_cell - num_cell/2, 1)
j_end   = min(    jj*num_cell + num_cell/2,Ny)

! now, calculate the distance between the obs and close candiate points
t=0
do i=i_start,i_end
 do j=j_start,j_end
    t=t+1
    lon_dif = (lon - LON_al(i,j))*111.41288*cos(lat*DEG2RAD) - &
              0.09350 * cos(3.0 * lat*DEG2RAD) + &
              0.00012 * cos(5.0 * lat*DEG2RAD)
    lat_dif = (lat - Lat_al(i,j))*111.13295
    rtemp(t) = sqrt(lon_dif*lon_dif+lat_dif*lat_dif)
    itemp(t) = i
    jtemp(t) = j
  enddo
enddo

!find minimum distances
min_location = minloc(rtemp)
loc1 = min_location(1)

!> @todo FIXME dies here with bounds checking turned on and you attempt
!>       to interpolate something outside the domain. This WHOLE ROUTINE
!>       can be replaced with something similar from POP.

!check point location
tp1 = checkpoint(LON_al(itemp(loc1)  ,jtemp(loc1)  ), &
                 LAT_al(itemp(loc1)  ,jtemp(loc1)  ), &
                 LON_al(itemp(loc1)+1,jtemp(loc1)+1), &
                 LAT_al(itemp(loc1)+1,jtemp(loc1)+1),lon,lat)

tp2 = checkpoint(LON_al(itemp(loc1)  ,jtemp(loc1)  ), &
                 LAT_al(itemp(loc1)  ,jtemp(loc1)  ), &
                 LON_al(itemp(loc1)  ,jtemp(loc1)+1), &
                 LAT_al(itemp(loc1)  ,jtemp(loc1)+1),lon,lat)

if(     (tp1  > 0.0_r8) .and. (tp2  > 0.0_r8)) then
      x_ind=itemp(loc1)-1
      y_ind=jtemp(loc1)
elseif( (tp1  > 0.0_r8) .and. (tp2  < 0.0_r8)) then
      x_ind=itemp(loc1)
      y_ind=jtemp(loc1)
elseif( (tp1  < 0.0_r8) .and. (tp2  > 0.0_r8)) then
      x_ind=itemp(loc1)-1
      y_ind=jtemp(loc1)-1
elseif( (tp1  < 0.0_r8) .and. (tp2  < 0.0_r8)) then
      x_ind=itemp(loc1)
      y_ind=jtemp(loc1)-1
elseif( (tp1 == 0.0_r8) .and. (tp2  > 0.0_r8)) then
      x_ind=itemp(loc1)-1
      y_ind=jtemp(loc1)
elseif( (tp1 == 0.0_r8) .and. (tp2  < 0.0_r8)) then
      x_ind=itemp(loc1)
      y_ind=jtemp(loc1)
elseif( (tp1  > 0.0_r8) .and. (tp2 == 0.0_r8)) then
      x_ind=itemp(loc1)
      y_ind=jtemp(loc1)
elseif( (tp1  < 0.0_r8) .and. (tp2 == 0.0_r8)) then
      x_ind=itemp(loc1)-1
      y_ind=jtemp(loc1)-1
else
   write(string1,*)'ERROR in finding matching grid point for'
   write(string2,*)'longitude ',lon
   write(string3,*)'latitude  ',lat
   call error_handler(E_MSG,'get_reg_box_indices:', string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

if (do_output() .and. debug > 0) then

   write(string1,*)'checking'
   write(string2,*)'longitude ',lon,' is index ',x_ind
   write(string3,*)'latitude  ',lat,' is index ',y_ind
   call error_handler(E_MSG,'get_reg_box_indices:', string1, &
              source, revision, revdate, text2=string2, text3=string3)

endif

deallocate(boxLON, boxLAT)
deallocate(rtemp)
deallocate(itemp)
deallocate(jtemp)

end subroutine get_reg_box_indices

!-----------------------------------------------------------------------
!>
!> Given a longitude and latitude in degrees returns the index of the
!> regular lon-lat box that contains the point. TJH ... this was the original
!> version that dies when trying to interpolate something outside the domain.
!> It dies with a run-time error when you have bounds-checking turned on.
!>
!> @param lon the longitude of interest
!> @param lat the latitude of interest
!> @param LON_al the matrix of longitudes defining the boxes (grids).
!> @param LAT_al the matrix of latitudes defining the boxes (grids).
!> @param x_ind the 'x' (longitude) index of the location of interest.
!> @param y_ind the 'y' (latitude) index of the location of interest.
!>
!> @todo FIXME should this really error out or silently fail?
!> At this point it may silently fail.
!>

subroutine get_reg_box_indices_org(lon, lat, LON_al, LAT_al, x_ind, y_ind)

real(r8), intent(in)  :: lon
real(r8), intent(in)  :: lat
real(r8), intent(in)  :: LON_al(:,:)
real(r8), intent(in)  :: LAT_al(:,:)
integer,  intent(out) :: x_ind
integer,  intent(out) :: y_ind

real(r8)              :: lon_dif,lat_dif,tp1,tp2
real(r8), allocatable :: rtemp(:)
integer,  allocatable :: itemp(:),jtemp(:)
integer               :: i,j,t,loc1
integer               :: min_location(1) ! stupid parser

! tricky for rotated domain

! calculate degree distance
! just to find corners

allocate(rtemp(size(LON_al,1)*size(LON_al,2)))
allocate(itemp(size(LON_al,1)*size(LON_al,2)))
allocate(jtemp(size(LON_al,1)*size(LON_al,2)))

t=0
do i=1,size(LON_al,1)
   do j=1,size(LON_al,2)
      t=t+1
      lon_dif = (lon - LON_al(i,j)) * 111.41288_r8 * cos(lat * DEG2RAD) - &
                0.09350_r8 * cos(3.0_r8 * lat * DEG2RAD) + &
                0.00012_r8 * cos(5.0_r8 * lat * DEG2RAD)
      lat_dif = (lat - Lat_al(i,j))*111.13295_r8
      rtemp(t) = sqrt(lon_dif*lon_dif+lat_dif*lat_dif)
      itemp(t) = i
      jtemp(t) = j
   enddo
enddo

!find minimum distances

min_location = minloc(rtemp)
loc1 = min_location(1)

!check point location

tp1 = checkpoint(LON_al(itemp(loc1)  , jtemp(loc1)  ), &
                 LAT_al(itemp(loc1)  , jtemp(loc1)  ), &
                 LON_al(itemp(loc1)+1, jtemp(loc1)+1), &
                 LAT_al(itemp(loc1)+1, jtemp(loc1)+1), lon, lat)

tp2 = checkpoint(LON_al(itemp(loc1), jtemp(loc1)  ), &
                 LAT_al(itemp(loc1), jtemp(loc1)  ), &
                 LON_al(itemp(loc1), jtemp(loc1)+1), &
                 LAT_al(itemp(loc1), jtemp(loc1)+1), lon, lat)

if(    (tp1  > 0.0_r8) .and. (tp2  > 0.0_r8)) then
   x_ind = itemp(loc1)-1
   y_ind = jtemp(loc1)
elseif((tp1  > 0.0_r8) .and. (tp2  < 0.0_r8)) then
   x_ind = itemp(loc1)
   y_ind = jtemp(loc1)
elseif((tp1  < 0.0_r8) .and. (tp2  > 0.0_r8)) then
   x_ind = itemp(loc1)-1
   y_ind = jtemp(loc1)-1
elseif((tp1  < 0.0_r8) .and. (tp2  < 0.0_r8)) then
   x_ind = itemp(loc1)
   y_ind = jtemp(loc1)-1
elseif((tp1 == 0.0_r8) .and. (tp2  > 0.0_r8)) then
   x_ind = itemp(loc1)-1
   y_ind = jtemp(loc1)
elseif((tp1 == 0.0_r8) .and. (tp2  < 0.0_r8)) then
   x_ind = itemp(loc1)
   y_ind = jtemp(loc1)
elseif((tp1  > 0.0_r8) .and. (tp2 == 0.0_r8)) then
   x_ind = itemp(loc1)
   y_ind = jtemp(loc1)
elseif((tp1  < 0.0_r8) .and. (tp2 == 0.0_r8)) then
   x_ind = itemp(loc1)-1
   y_ind = jtemp(loc1)-1
else
   write(string1,*)'ERROR in finding matching grid point for'
   write(string2,*)'longitude ',lon
   write(string3,*)'latitude  ',lat
   call error_handler(E_MSG,'get_reg_box_indices:', string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

deallocate(rtemp)
deallocate(itemp)
deallocate(jtemp)

end subroutine get_reg_box_indices_org


!-----------------------------------------------------------------------
!>
!> check (xp,yp) with line (x1,y1)-->(x2,y2)
!> if checkpoint>0, (xp,yp) is on the left side of the line
!> if checkpoint<0, (xp,yp) is on the right side of the line
!> counterclockwise direction
!>

real(r8) FUNCTION  checkpoint(x1,y1,x2,y2,xp,yp)

real(r8), intent(in) :: x1
real(r8), intent(in) :: y1
real(r8), intent(in) :: x2
real(r8), intent(in) :: y2
real(r8), intent(in) :: xp
real(r8), intent(in) :: yp

real(r8) :: A,B,C

A = -(y2-y1)
B = x2-x1
C = -(A*x1+B*y1)
checkpoint = A*xp+B*yp+C

END FUNCTION checkpoint


!-----------------------------------------------------------------------
!>
!> Grabs the corners for a given quadrilateral from the global array of lower
!> right corners.
!>
!> @param x the matrix of lower right corner locations
!> @param i the 'x' or longitude index
!> @param j the 'y' or latitude index
!> @param corners the 2x2 matrix of the corner values

subroutine get_quad_corners(x, i, j, corners)

real(r8), intent(in)  :: x(:, :)
integer,  intent(in)  :: i
integer,  intent(in)  :: j
real(r8), intent(out) :: corners(4)

integer :: ip1

! Note that corners go counterclockwise around the quad.
! @todo Have to worry about wrapping in longitude but not in latitude

ip1 = i + 1
!if(ip1 > nx) ip1 = 1

corners(1) = x(i,   j  )
corners(2) = x(ip1, j  )
corners(3) = x(ip1, j+1)
corners(4) = x(i,   j+1)

end subroutine get_quad_corners


!-----------------------------------------------------------------------
!>
!> set a netCDF dimension ID array needed to augment the natural shape
!> of the variable with the two additional dimids needed by the DART
!> diagnostic output. DART variables have an additional 'copy' and
!> 'unlimited' (time) dimension.
!>
!> @param ncid the netCDF handle
!> @param ivar the DART variable index
!> @param memberdimid the netCDF dimension ID of the 'copy' coordinate
!> @param unlimiteddimid the netCDF dimension ID of the 'time' coordinate
!> @param ndims the (extended) number of dimensions for this DART variable
!> @param dimids the (extended) dimension IDs
!>

subroutine define_var_dims(ncid, ivar, memberdimid, unlimiteddimid, ndims, dimids)

integer, intent(in)  :: ncid
integer, intent(in)  :: ivar
integer, intent(in)  :: memberdimid
integer, intent(in)  :: unlimiteddimid
integer, intent(out) :: ndims
integer, intent(out) :: dimids(:)

integer :: i,mydimid

ndims  = 0
dimids = 0

do i = 1,progvar(ivar)%numdims

   ! Each of these dimension names
   ! must exist in the DART diagnostic netcdf files.

   call nc_check(nf90_inq_dimid(ncid, trim(progvar(ivar)%dimname(i)), mydimid), &
              'define_var_dims','inq_dimid '//trim(progvar(ivar)%dimname(i)))

   ndims = ndims + 1

   dimids(ndims) = mydimid

enddo

ndims         = ndims + 1
dimids(ndims) = memberdimid
ndims         = ndims + 1
dimids(ndims) = unlimiteddimid

end subroutine define_var_dims


!-----------------------------------------------------------------------
!>
!> gets the length of a netCDF dimension given the dimension name.
!> This bundles the nf90_inq_dimid and nf90_inquire_dimension routines
!> into a slightly easier-to-use function.
!>
!> @param dimlen the length of the netCDF dimension in question
!> @param ncid the netCDF file handle
!> @param dimension_name the character string of the dimension name
!> @param filename the name of the netCDF file (for error message purposes)
!>

function get_dimension_length(ncid, dimension_name, filename) result(dimlen)

integer                      :: dimlen
integer,          intent(in) :: ncid
character(len=*), intent(in) :: dimension_name
character(len=*), intent(in) :: filename

integer :: DimID

write(string1,*)'inq_dimid '//trim(dimension_name)//' '//trim(filename)
write(string2,*)'inquire_dimension '//trim(dimension_name)//' '//trim(filename)

call nc_check(nf90_inq_dimid(ncid, trim(dimension_name), DimID), &
              'get_grid',string1)
call nc_check(nf90_inquire_dimension(ncid, DimID, len=dimlen), &
              'get_grid', string2)

end function get_dimension_length

!===================================================================
! End of model_mod
!===================================================================

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
