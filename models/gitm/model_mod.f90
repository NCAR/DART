! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This is the interface between the gitm model and DART.

! Modules that are absolutely required for use are listed

use        types_mod, only : r4, r8, digits12, SECPERDAY, MISSING_R8, MISSING_I, &
                             rad2deg, deg2rad, PI, obstypelength, i8

use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, get_dist, query_location,          &
                             get_close_type, VERTISHEIGHT, VERTISLEVEL,        &
                             set_location, get_location, VERTISUNDEF,          &
                             loc_get_close_obs => get_close_obs, is_vertical,  &
                             loc_get_close_state => get_close_state,           &
                             vertical_localization_on, set_vertical_localization_coord

use    utilities_mod, only : register_module, error_handler, string_to_logical,&
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             do_output, to_upper, string_to_real,              &
                             find_namelist_in_file, check_namelist_read,       &
                             nmlfileunit, do_nml_file, do_nml_term,            &
                             open_file, file_exist, find_textfile_dims,        &
                             file_to_text, close_file, find_enclosing_indices

use netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file, nc_create_file, &
                                 nc_synchronize_file, nc_add_attribute_to_variable,    &
                                 nc_begin_define_mode, nc_end_define_mode,             &
                                 nc_add_global_attribute,                              &
                                 nc_put_variable, nc_get_dimension_size,               &
                                 nc_define_double_variable, nc_define_dimension,       &
                                 nc_get_dimension_size, nc_get_variable,               &
                                 nc_define_double_scalar

use  ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod,  only : get_state

use   state_structure_mod,  only : add_domain, get_dart_vector_index, get_domain_size, &
                                   get_dim_name, get_kind_index, get_num_dims, &
                                   get_num_variables, get_varid_from_kind, &
                                   get_model_variable_indices, state_structure_info, &
                                   get_index_start

use     obs_kind_mod  ! all for now - fixme  !only : get_index_for_quantity,  &
                                             !get_name_for_quantity,   &
                                             !QTY_GEOPOTENTIAL_HEIGHT

use        quad_utils_mod,  only : quad_interp_handle, init_quad_interp, &
                                   set_quad_coords, finalize_quad_interp, &
                                   quad_lon_lat_locate, quad_lon_lat_evaluate, &
                                   GRID_QUAD_IRREG_SPACED_REGULAR,  &
                                   QUAD_LOCATED_CELL_CENTERS

use     default_model_mod,  only : adv_1step, nc_write_model_vars, &
                                   pert_model_copies, &
                                   read_model_time, write_model_time, &
                                   convert_vertical_obs, convert_vertical_state, &
                                   init_time => fail_init_time, &
                                   init_conditions => fail_init_conditions

use     dart_time_io_mod,   only : read_model_time, write_model_time

use     dart_gitm_mod,      only : get_nSpecies, get_nSpeciesTotal, get_nIons, &
                                   get_nSpeciesAll, get_nLonsPerBlock,         &
                                   get_nLatsPerBlock, get_nAltsPerBlock,       &
                                   decode_gitm_indices

use     ModPlanet,          only : ie_

use netcdf

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
! this list has code in this module.
public :: get_model_size,         &
          get_state_meta_data,    &
          model_interpolate,      &
          shortest_time_between_assimilations, &
          static_init_model,      &
          end_model,              &
          get_close_obs,          &
          get_close_state,        &
          convert_vertical_obs,   &
          convert_vertical_state, &
          nc_write_model_atts     

! these routines also must be public.
! this list are names of routines where the code
! is passed through from other modules
public :: init_time,              &
          init_conditions,        &
          adv_1step,              &
          nc_write_model_vars,    &
          pert_model_copies,      &
          read_model_time,        &
          write_model_time 

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.

public :: get_state_time,          &
          restart_files_to_netcdf, &  ! only used in converter
          netcdf_to_restart_files     ! only used in converter

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2
logical, save :: module_initialized = .false.

character(len=*), parameter :: LON_DIM_NAME = 'LON'
character(len=*), parameter :: LAT_DIM_NAME = 'LAT'
character(len=*), parameter :: ALT_DIM_NAME = 'ALT'

character(len=*), parameter :: LON_VAR_NAME = 'LON'
character(len=*), parameter :: LAT_VAR_NAME = 'LAT'
character(len=*), parameter :: ALT_VAR_NAME = 'ALT'

integer, parameter :: MAX_NAME_LEN = 256

!------------------------------------------------------------------
! things which can/should be in the model_nml
!
!  The DART state vector may consist of things like:
!
!  U    long_name = "X-WIND COMPONENT"      float   U(TIME, ALT, LAT, XE)
!  V    long_name = "Y-WIND COMPONENT"      float   V(TIME, ALT, YE, LON)
!  W    long_name = "Z-WIND COMPONENT"      float   W(TIME, ZE, LAT, LON)
!  TH   long_name = "POTENTIAL TEMPERATURE" float  TH(TIME, ALT, LAT, LON)
!  DBZ  long_name = "RADAR REFLECTIVITY"    float DBZ(TIME, ALT, LAT, LON)
!  WZ   long_name = "VERTICAL VORTICITY"    float  WZ(TIME, ALT, LAT, LON)
!  PI   long_name = "PERT. EXNER"	    float  PI(TIME, ALT, LAT, LON)
!  QV   long_name = "VAPOR MIXING RATIO"    float  QV(TIME, ALT, LAT, LON)
!  QC   long_name = "CLOUD MIXING RATIO"    float  QC(TIME, ALT, LAT, LON)
!  QR   long_name = "RAIN MIXING RATIO"     float  QR(TIME, ALT, LAT, LON)
!  QI   long_name = "ICE MIXING RATIO"      float  QI(TIME, ALT, LAT, LON)
!  QS   long_name = "SNOW MIXING RATIO"     float  QS(TIME, ALT, LAT, LON)
!  QH   long_name = "GRAUPEL MIXING RATIO"  float  QH(TIME, ALT, LAT, LON)
!
!  The variables in the gitm restart file that are used to create the
!  DART state vector are specified in input.nml:model_nml: gitm_state_variables
!
!------------------------------------------------------------------

integer, parameter :: max_state_variables = 80
integer, parameter :: num_state_table_columns = 5
character(len=MAX_NAME_LEN)  :: gitm_state_variables(max_state_variables * num_state_table_columns ) = ' '
logical            :: single_file_in = .false.
integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 60
real(r8)           :: model_perturbation_amplitude = 0.2
integer            :: debug = 0   ! turn up for more and more debug messages
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: template_filename = 'no_file_specified.nc'

namelist /model_nml/            &
   single_file_in,              &
   template_filename,           &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   calendar,                    &
   debug,                       &
   !nLonsPerBlock,               &
   !nLatsPerBlock,               &
   !nAltsPerBlock,               &
   gitm_state_variables

character(len=MAX_NAME_LEN)  :: gitm_block_variables(max_state_variables) = ' '

namelist /gitm_blocks_nml/      &
   gitm_block_variables

integer :: nfields

! Everything needed to describe a GITM variable
! NOTE: these are used to convert the GITM blocks to netCDF files
! They are not used during the assimilation step
type gitmvartype
   private
   character(len=NF90_MAX_NAME) :: varname       ! crazy species name
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=NF90_MAX_NAME) :: gitm_varname  ! NDensityS, IDensityS, ...
   integer :: gitm_dim                           ! dimension defining species
   integer :: gitm_index                         ! 'iSpecies' or u,v,w ...
end type gitmvartype

type(gitmvartype), dimension(max_state_variables) :: gitmvar

type(quad_interp_handle) :: quad_interp

! this id allows us access to all of the state structure
! info and is required for getting state variables.
integer :: domain_id

! nLonsPerBlock, nLatsPerBlock are the number of lons/lats PER block
!                  the number of blocks comes from UAM.in
! nAltsPerBlock  is the number of altitudes, which does not depend on block

integer :: nLonsPerBlock, nLatsPerBlock, nAltsPerBlock

! "... keep in mind that if the model resolution is 5 deg latitude,
!  the model will actually go from -87.5 to 87.5 latitude
! (even though you specify -90 to 90 in the UAM.in file),
! since the latitudes/longitudes are at cell centers,
! while the edges are at the boundaries." -- Aaron Ridley

integer  :: NgridLon=-1, NgridLat=-1, NgridAlt=-1    ! scalar grid counts
integer  :: nBlocksLon=-1, nBlocksLat=-1             ! number of blocks along each dim
real(r8) :: LatStart=MISSING_R8, LatEnd=MISSING_R8, LonStart=MISSING_R8
integer  :: nSpeciesTotal=-1, nSpecies=-1, nIons=-1, nSpeciesAll=-1

! scalar grid positions

real(r8), allocatable :: LON(:)   ! longitude centers
real(r8), allocatable :: LAT(:)   ! latitude  centers
real(r8), allocatable :: ALT(:)   ! vertical level centers

type(time_type)       :: model_time      ! valid time of the model state
type(time_type)       :: model_advance_time  ! smallest time to adv model

! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.

logical :: print_timestamps = .false.
integer :: print_every_Nth  = 10000

integer, parameter :: nGhost = 2   ! number of ghost cells on all edges

!       ERROR codes:
!
!       99:  general error in case something terrible goes wrong...
!       15:  dont know what to do with vertical coord supplied
!       16:  lat/lon illegal
!       17:  altitude illegal
!       20:  asking to interpolate an unknown obs quantity
integer, parameter :: GENERAL_ERROR_CODE = 99
integer, parameter :: INVALID_VERT_COORD_ERROR_CODE = 15
integer, parameter :: INVALID_LATLON_VAL_ERROR_CODE = 16
integer, parameter :: INVALID_ALTITUDE_VAL_ERROR_CODE = 17
integer, parameter :: UNKNOWN_OBS_QTY_ERROR_CODE = 20

!------------------------------------------------------------------
! The gitm restart manager namelist variables
!------------------------------------------------------------------

contains


!==================================================================
! All the REQUIRED interfaces come first - just by convention.
!==================================================================

!==================================================================
!>
!> Returns the size of the DART state vector (i.e. model) as an integer.

function get_model_size()

integer(i8) :: get_model_size

character(len=*), parameter :: routine = 'get_model_size'

if ( .not. module_initialized ) then
      write(string1, *) 'Invalid state - the module has not yet been initialized through static_init_model'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
end if

get_model_size = get_domain_size(domain_id)

end function get_model_size


!==================================================================
!>
!> given an index into the state vector, return its location and
!> if given, the var kind.   despite the name, var_type is a generic
!> kind, like those in obs_kind/obs_kind_mod.f90, starting with QTY_

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

character(len=*), parameter :: routine = 'get_state_meta_data'

! Local variables

integer :: lat_index, lon_index, alt_index
integer :: myvarid, myqty

if ( .not. module_initialized ) then
      write(string1, *) 'Invalid state - the module has not yet been initialized through static_init_model'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
end if

call get_model_variable_indices(index_in, lon_index, lat_index, alt_index, &
   var_id=myvarid, kind_index=myqty)

location = set_location(LON(lon_index), LAT(lat_index), ALT(alt_index), VERTISHEIGHT)

if (present(var_type)) then
   var_type = myqty
endif

end subroutine get_state_meta_data


!==================================================================
!
!>     PURPOSE:
!>
!>     For a given lat, lon, and height, interpolate the correct state value
!>     to that location for the filter from the gitm state vectors
!>
!>     Variables needed to be stored in the MODEL_MODULE data structure
!>
!>       LON   = 1D array storing the local grid center coords (degrees)
!>       LAT   = 1D array storing the local grid center coords (degrees)
!>       ALT   = 1D array storing the local grid center coords (meters)
!>
!> Passed variables

subroutine model_interpolate(state_handle, ens_size, location, obs_qty, interp_vals, status_array)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_qty
real(r8),            intent(out) :: interp_vals(ens_size)
integer,             intent(out) :: status_array(ens_size)

! Local storage

character(len=*), parameter :: routine = 'model_interpolate'

real(r8) :: loc_array(3), llon, llat, lvert, lon_fract, lat_fract
integer  :: four_lons(4), four_lats(4)
integer  :: status1, which_vert, varid
real(r8) :: quad_vals(4, ens_size)

if ( .not. module_initialized ) then
      write(string1, *) 'Invalid state - the module has not yet been initialized through static_init_model'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
end if

! Assume failure.  Set return val to missing, then the code can
! just set status_array to something indicating why it failed, and return.
! If the interpolation is good, interp_vals will be set to the
! good values, and the last line here sets status_array to 0.
! make any error codes set here be in the 10s

interp_vals = MISSING_R8          ! the DART bad value flag
status_array = GENERAL_ERROR_CODE ! unknown error

! Get the individual locations values

loc_array  = get_location(location)
llon       = loc_array(1)
llat       = loc_array(2)
lvert      = loc_array(3)
which_vert = nint(query_location(location))

IF (debug > 85) then
   write(string1,*)  'requesting interpolation at ', llon, llat, lvert
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
end if

! Only height and level for vertical location type is supported at this point
if (.not. is_vertical(location, "HEIGHT") .and. .not. is_vertical(location, "LEVEL")) THEN
     status_array = INVALID_VERT_COORD_ERROR_CODE
     return
endif

if (obs_qty == QTY_GEOMETRIC_HEIGHT .and. is_vertical(location, "LEVEL")) then
   if (nint(lvert) < 1 .or. nint(lvert) > size(ALT,1)) then
      interp_vals = MISSING_R8
      status_array = 1
   else
      interp_vals = ALT(nint(lvert))
      status_array = 0
   endif
   return ! Early Return
endif

! do we know how to interpolate this quantity?
call ok_to_interpolate(obs_qty, varid, status1)

if (status1 /= 0) then
   if(debug > 12) then
      write(string1,*) 'Did not find observation quantity ', obs_qty, ' in the state vector'
      call error_handler(E_WARN,routine,string1,source,revision,revdate)
   endif
   status_array(:) = status1   ! this quantity not in the state vector
   return
endif

! get the indices for the 4 corners of the quad in the horizontal, plus
! the fraction across the quad for the obs location
call quad_lon_lat_locate(quad_interp, llon, llat, & 
                         four_lons, four_lats, lon_fract, lat_fract, status1)
if (status1 /= 0) then
   status_array(:) = INVALID_LATLON_VAL_ERROR_CODE  ! cannot locate enclosing horizontal quad
   return
endif

call get_quad_vals(state_handle, ens_size, varid, four_lons, four_lats, &
                   loc_array, which_vert, quad_vals, status_array)
if (any(status_array /= 0)) return

! do the horizontal interpolation for each ensemble member
call quad_lon_lat_evaluate(quad_interp, lon_fract, lat_fract, ens_size, &
                           quad_vals, interp_vals, status_array)

! All good.
status_array(:) = 0

end subroutine model_interpolate


!==================================================================
!>
!> Returns the the time step of the model; the smallest increment
!> in time that the model is capable of advancing the state in a given
!> implementation. This interface is required for all applications.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

character(len=*), parameter :: routine = 'shortest_time_between_assimilations'

if ( .not. module_initialized ) then
      write(string1, *) 'Invalid state - the module has not yet been initialized through static_init_model'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
end if

shortest_time_between_assimilations = model_advance_time

end function shortest_time_between_assimilations


!==================================================================
!>
!> Called to do one time initialization of the model.
!>
!> All the grid information comes from the initialization of
!> the dart_gitm_mod module.

subroutine static_init_model()

integer :: ss, dd

character(len=*), parameter :: routine = 'static_init_model'

if ( module_initialized ) return ! only need to do this once.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! routines which are used here should *not* call this routine.
! it can mess up the stack unless you say this is a recursive routine.  
! rather than go there just don't call this from anything that is used by it.

module_initialized = .true.

call read_model_namelist()

! Get the GITM variables in a restricted scope setting.

!---------------------------------------------------------------
! Set the time step ... causes gitm namelists to be read.
! Ensures model_advance_time is multiple of 'dynamics_timestep'

call set_calendar_type( calendar )   ! comes from model_mod_nml

model_advance_time = set_model_time_step()

call get_time(model_advance_time,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)

!---------------------------------------------------------------
! get grid dimensions and values

write(string1,'(3A)') "Now reading template file ",trim(template_filename),&
   " for grid information"
call error_handler(E_MSG,routine,string1,source,revision,revdate)

call get_grid_info_from_netcdf(template_filename, NgridLon, NgridLat, NgridAlt)

allocate(LON(NgridLon))
allocate(LAT(NgridLat))
allocate(ALT(NgridAlt))

!---------------------------------------------------------------
! get grid dimensions and values
call get_grid_from_netcdf(template_filename, LON, LAT, ALT)

!---------------------------------------------------------------

! mass points at cell centers
call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, nGridLon, nGridLat, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.false., spans_lon_zero=.false., pole_wrap=.false., &
                      interp_handle=quad_interp)

call set_quad_coords(quad_interp, LON, LAT)

! this calls add_domain() for us.   we should stop cutting and pasting
! this code and make a common routine if it has the same number of columns
! as other models.  replicated code in multiple model_mods is a maintenance
! overhead we should try to avoid.

call set_gitm_variable_info(gitm_state_variables)

if ( debug > 0 ) then
   write(string1,'("grid: NgridLon, NgridLat, NgridAlt =",3(1x,i5))') NgridLon, NgridLat, NgridAlt
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
endif

! needs to set the vertical localization coordinate, too.
call set_vertical_localization_coord(VERTISHEIGHT)   ! not sure which?

end subroutine static_init_model


!==================================================================

subroutine read_model_namelist()

integer :: iunit, io

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

end subroutine

!==================================================================

!> Does any shutdown and clean-up needed for model.

subroutine end_model()

if (allocated(LON)) deallocate(LON, LAT, ALT)
call finalize_quad_interp(quad_interp)

end subroutine end_model


!==================================================================


subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid 
integer, intent(in) :: domain_id

call nc_begin_define_mode(ncid)

call add_nc_definitions(ncid)

!----------------------------------------------------------------------------
! Finished with dimension/variable definitions, must end 'define' mode to fill.
!----------------------------------------------------------------------------

call nc_end_define_mode(ncid)

call add_nc_dimvars(ncid)

end subroutine nc_write_model_atts

!==================================================================

subroutine add_nc_definitions(ncid)

integer, intent(in) :: ncid 

call nc_add_global_attribute(ncid, 'model', 'gitm')

!-------------------------------------------------------------------------------
! Determine shape of most important namelist
!-------------------------------------------------------------------------------
!
!call find_textfile_dims('gitm_vars.nml', nlines, linelen)
!if (nlines > 0) then
!   has_gitm_namelist = .true.
!
!   allocate(textblock(nlines))
!   textblock = ''
!
!   call nc_define_dimension(ncid, 'nlines',  nlines)
!   call nc_define_dimension(ncid, 'linelen', linelen)
!   call nc_define_character_variable(ncid, 'gitm_in', (/ 'nlines ', 'linelen' /))
!   call nc_add_attribute_to_variable(ncid, 'gitm_in', 'long_name', 'contents of gitm_in namelist')
!
!else
!  has_gitm_namelist = .false.
!endif
!
!----------------------------------------------------------------------------
! output only grid info - state vars will be written by other non-model_mod code
!----------------------------------------------------------------------------

call nc_define_dimension(ncid, LON_DIM_NAME, NgridLon)
call nc_define_dimension(ncid, LAT_DIM_NAME, NgridLat)
call nc_define_dimension(ncid, ALT_DIM_NAME, NgridAlt)
call nc_define_dimension(ncid, 'WL',  1)  ! wavelengths - currently only 1?

!----------------------------------------------------------------------------
! Create the (empty) Coordinate Variables and the Attributes
!----------------------------------------------------------------------------

! Grid Longitudes
call nc_define_double_variable(ncid, LON_VAR_NAME, (/ LON_DIM_NAME /) )
call nc_add_attribute_to_variable(ncid, LON_VAR_NAME, 'type',           'x1d')
call nc_add_attribute_to_variable(ncid, LON_VAR_NAME, 'long_name',      'grid longitudes')
call nc_add_attribute_to_variable(ncid, LON_VAR_NAME, 'cartesian_axis', 'X')
call nc_add_attribute_to_variable(ncid, LON_VAR_NAME, 'units',          'degrees_east')
call nc_add_attribute_to_variable(ncid, LON_VAR_NAME, 'valid_range',     (/ 0.0_r8, 360.0_r8 /) )

! Grid Latitudes
call nc_define_double_variable(ncid, LAT_VAR_NAME, (/ LAT_DIM_NAME /) )
call nc_add_attribute_to_variable(ncid, LAT_VAR_NAME, 'type',           'y1d')
call nc_add_attribute_to_variable(ncid, LAT_VAR_NAME, 'long_name',      'grid latitudes')
call nc_add_attribute_to_variable(ncid, LAT_VAR_NAME, 'cartesian_axis', 'Y')
call nc_add_attribute_to_variable(ncid, LAT_VAR_NAME, 'units',          'degrees_north')
call nc_add_attribute_to_variable(ncid, LAT_VAR_NAME, 'valid_range',     (/ -90.0_r8, 90.0_r8 /) )

! Grid Altitudes
call nc_define_double_variable(ncid, ALT_VAR_NAME, (/ ALT_DIM_NAME /) )
call nc_add_attribute_to_variable(ncid, ALT_VAR_NAME, 'type',           'z1d')
call nc_add_attribute_to_variable(ncid, ALT_VAR_NAME, 'long_name',      'grid altitudes')
call nc_add_attribute_to_variable(ncid, ALT_VAR_NAME, 'cartesian_axis', 'Z')
call nc_add_attribute_to_variable(ncid, ALT_VAR_NAME, 'units',          'meters')
call nc_add_attribute_to_variable(ncid, ALT_VAR_NAME, 'positive',       'up')

! Grid wavelengths
call nc_define_double_variable(ncid, 'WL', (/ 'WL' /) )
call nc_add_attribute_to_variable(ncid, 'WL', 'type',           'x1d')
call nc_add_attribute_to_variable(ncid, 'WL', 'long_name',      'grid wavelengths')
call nc_add_attribute_to_variable(ncid, 'WL', 'cartesian_axis', 'X')
call nc_add_attribute_to_variable(ncid, 'WL', 'units',          'wavelength_index')
call nc_add_attribute_to_variable(ncid, 'WL', 'valid_range',     (/ 0.9_r8, 38.1_r8 /) )
end subroutine add_nc_definitions


!==================================================================

subroutine add_nc_dimvars(ncid)

integer, intent(in) :: ncid 

!----------------------------------------------------------------------------
! Fill the coordinate variables
!----------------------------------------------------------------------------

call nc_put_variable(ncid, LON_VAR_NAME, LON)
call nc_put_variable(ncid, LAT_VAR_NAME, LAT)
call nc_put_variable(ncid, ALT_VAR_NAME, ALT)
! what about WL?

!if (has_gitm_namelist) then
!   call file_to_text('gitm_vars.nml', textblock)
!   call nc_put_variable(ncid, 'gitm_in', textblock)
!   deallocate(textblock)
!endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_synchronize_file(ncid)

end subroutine add_nc_dimvars



!==================================================================

!> Given a DART location (referred to as "base") and a set of candidate
!> locations & kinds (obs, obs_kind), returns the subset close to the
!> "base", their indices, and their distances to the "base" ...
!>
!> For vertical distance computations, general philosophy is to convert all
!> vertical coordinates to a common coordinate. This coordinate type is defined
!> in the namelist with the variable "vert_localization_coord".

subroutine get_close_obs(gc, base_obs_loc, base_obs_type, &
                         locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

type(get_close_type),              intent(in)    :: gc
type(location_type),               intent(inout) :: base_obs_loc
integer,                           intent(in)    :: base_obs_type
type(location_type),               intent(inout) :: locs(:)
integer,                           intent(in)    :: loc_qtys(:)
integer,                           intent(in)    :: loc_types(:)
integer,                           intent(out)   :: num_close
integer,                           intent(out)   :: close_ind(:)
real(r8),            optional,     intent(out)   :: dist(:)
type(ensemble_type), optional,     intent(in)    :: ens_handle

integer                :: t_ind, istatus1, istatus2, k, is_in_close_ind, is_in_obs_kind
integer                ::  f107_ind
integer                :: base_which, local_obs_which, vert_type
real(r8), dimension(3) :: base_array, local_obs_array
type(location_type)    :: local_obs_loc

! Initialize variables to missing status

num_close       = 0
close_ind       = -99
dist            = 1.0e9_r8   !something big and positive (far away)
istatus1        = 0
istatus2        = 0
is_in_obs_kind  = 0
is_in_close_ind = 0
f107_ind        = -37 !a bad index, hopefully out of bounds of obs_kind


! if absolute distances aren't needed, or vertical localization isn't on,
! the default version works fine since no conversion will be needed and
! there won't be any damping since there are no vert distances.
if (.not. present(dist) .or. .not. vertical_localization_on()) then
   call loc_get_close_obs(gc, base_obs_loc, base_obs_type, locs, loc_qtys, loc_types, &
                            num_close, close_ind, dist, ens_handle)
   return
endif

! Convert base_obs vertical coordinate to requested vertical coordinate if necessary
! does the base obs need conversion first?
vert_type = query_location(base_obs_loc)

!> @todo FIXME: if we start supporting different vert coords, 
!> this code needs to get smarter.
if (vertical_localization_on() .and. vert_type /= VERTISUNDEF .and. vert_type /= VERTISHEIGHT) then
   call error_handler(E_ERR, 'get_close_state', 'cannot do vertical conversion of base obs', &
                      source, revision, revdate)
endif


! Loop over potentially close subset of obs priors or state variables
! This way, we are decreasing the number of distance computations that will follow.
! This is a horizontal-distance operation and we don't need to have the relevant
! vertical coordinate information yet (for obs_loc).

call loc_get_close_obs(gc, base_obs_loc, base_obs_type, locs, loc_qtys, loc_types, &
                            num_close, close_ind, dist, ens_handle)


!!!! THE following 20-ish+ lines are implementing the search (if f107's dist to obs is to be calculated)
!!!! Alex 03/07/2012
!   do i = 1, size(obs_kind) !have to go over the whole size because these are all the candidates
!      if (obs_kind(i) .eq. get_index_for_quantity('QTY_1D_PARAMETER')) then !so right now any QTY_1D_PARAMETER will match.
!+ right now the only parameter is f107, but if you add more parameters, you might want to change their localizations, as
!+ right now they will be either all at the meas. location or all far (depending on est_f107 setting in pbs_file.sh)
!         is_in_obs_kind = 1 !true
!         f107_ind = i !its index
!      endif
!   enddo
!   if (is_in_obs_kind == 1) then !only check the close_ind if f107 needs to be added
!      do k = 1, num_close !go only as far as the data is reasonable (not -99 = data missing)
!         if (close_ind(k) .eq. f107_ind) then !if is already in close_ind, take note of it
!            is_in_close_ind = 1
!         endif
!      enddo
!   endif
!   if ((is_in_obs_kind == 1) .and. (is_in_close_ind == 0)) then !if it needs to be added (is in obs_kind), but is not added yet
!      num_close = num_close + 1
!      close_ind(num_close) = f107_ind
!      write(*,*) "F107 ADDED, n_c, f107_i ", num_close, f107_ind
!   endif


! i don't think we need any of this but i'm unsure how the F10.7 code
! affects this.

!   do k = 1, num_close
!
!      t_ind = close_ind(k)
!      local_obs_loc   = obs_loc(t_ind)
!      local_obs_which = nint(query_location(local_obs_loc))
!
!      ! Convert vertical coordinate to requested vertical coordinate if necessary.
!      ! Should only be necessary for obs priors, as state location information already
!      ! contains the correct vertical coordinate
!      ! (filter_assim's call to get_state_meta_data).
!      if (vertical_localization_on()) then
! !fixme       if (local_obs_which /= wrf%dom(1)%vert_coord) then
! !fixme           call vert_interpolate(ens_mean, local_obs_loc, obs_kind(t_ind), istatus2)
!            ! Store the "new" location into the original full local array
!            obs_loc(t_ind) = local_obs_loc
! !fixme        endif
!      endif
!
!      ! Compute distance - set distance to a very large value if vert coordinate is
!      ! missing or vert_interpolate returned error (istatus2=1)
!      local_obs_array = get_location(local_obs_loc)
!
!
!      ! TJH FIXME ... the pbs_file script actually modifies the value of dist in
!      ! TJH FIXME ... this file and recompiles. NOT APPROPRIATE.
!      if (((vertical_localization_on())        .and. &
!           (local_obs_array(3) == MISSING_R8)) .or.  &
!           (istatus2 == 1)                   ) then
!         dist(k) = 1.0e9_r8
!      else
!         !if (close_ind(k) .eq. f107_ind) then !check if we came across the parameter
!         !   dist(k) = 0 !changed by pbs_file script
!         !else
!            dist(k) = get_dist(base_obs_loc, local_obs_loc, base_obs_type, obs_kind(t_ind))
!         !endif
!      endif
!   enddo

end subroutine get_close_obs


!==================================================================

! Given a DART location (referred to as "base") and a set of candidate
! locations & kinds (obs, obs_kind), returns the subset close to the
! "base", their indices, and their distances to the "base" ...
!
! For vertical distance computations, general philosophy is to convert all
! vertical coordinates to a common coordinate. This coordinate type is defined
! in the namelist with the variable "vert_localization_coord".

subroutine get_close_state(gc, base_obs_loc, base_obs_type, &
                           locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)

type(get_close_type), intent(in)    :: gc
type(location_type),  intent(inout) :: base_obs_loc
integer,              intent(in)    :: base_obs_type
type(location_type),  intent(inout) :: locs(:)
integer,              intent(in)    :: loc_qtys(:)
integer(i8),          intent(in)    :: loc_indx(:)
integer,              intent(out)   :: num_close
integer,              intent(out)   :: close_ind(:)
real(r8),             optional, intent(out)   :: dist(:)
type(ensemble_type),  optional, intent(in)    :: ens_handle

integer                :: t_ind, istatus1, istatus2, k, is_in_close_ind, is_in_obs_kind
integer                ::  f107_ind
integer                :: base_which, local_obs_which, vert_type
real(r8), dimension(3) :: base_array, local_obs_array
type(location_type)    :: local_obs_loc

! Initialize variables to missing status

num_close       = 0
close_ind       = -99
dist            = 1.0e9_r8   !something big and positive (far away)
istatus1        = 0
istatus2        = 0
is_in_obs_kind  = 0
is_in_close_ind = 0
f107_ind        = -37 !a bad index, hopefully out of bounds of obs_kind


! if absolute distances aren't needed, or vertical localization isn't on,
! the default version works fine since no conversion will be needed and
! there won't be any damping since there are no vert distances.
if (.not. present(dist) .or. .not. vertical_localization_on()) then
   call loc_get_close_state(gc, base_obs_loc, base_obs_type, locs, loc_qtys, loc_indx, &
                            num_close, close_ind, dist, ens_handle)
   return
endif

! Convert base_obs vertical coordinate to requested vertical coordinate if necessary
! does the base obs need conversion first?
vert_type = query_location(base_obs_loc)

!> @todo FIXME: if we start supporting different vert coords, 
!> this code needs to get smarter.
if (vertical_localization_on() .and. vert_type /= VERTISUNDEF .and. vert_type /= VERTISHEIGHT) then
   call error_handler(E_ERR, 'get_close_state', 'cannot do vertical conversion of base obs', &
                      source, revision, revdate)
endif


! Loop over potentially close subset of obs priors or state variables
! This way, we are decreasing the number of distance computations that will follow.
! This is a horizontal-distance operation and we don't need to have the relevant
! vertical coordinate information yet (for obs_loc).

call loc_get_close_state(gc, base_obs_loc, base_obs_type, locs, loc_qtys, loc_indx, &
                            num_close, close_ind, dist, ens_handle)


!!!! THE following 20-ish+ lines are implementing the search (if f107's dist to obs is to be calculated)
!!!! Alex 03/07/2012
!   do i = 1, size(obs_kind) !have to go over the whole size because these are all the candidates
!      if (obs_kind(i) .eq. get_index_for_quantity('QTY_1D_PARAMETER')) then !so right now any QTY_1D_PARAMETER will match.
!+ right now the only parameter is f107, but if you add more parameters, you might want to change their localizations, as
!+ right now they will be either all at the meas. location or all far (depending on est_f107 setting in pbs_file.sh)
!         is_in_obs_kind = 1 !true
!         f107_ind = i !its index
!      endif
!   enddo
!   if (is_in_obs_kind == 1) then !only check the close_ind if f107 needs to be added
!      do k = 1, num_close !go only as far as the data is reasonable (not -99 = data missing)
!         if (close_ind(k) .eq. f107_ind) then !if is already in close_ind, take note of it
!            is_in_close_ind = 1
!         endif
!      enddo
!   endif
!   if ((is_in_obs_kind == 1) .and. (is_in_close_ind == 0)) then !if it needs to be added (is in obs_kind), but is not added yet
!      num_close = num_close + 1
!      close_ind(num_close) = f107_ind
!      write(*,*) "F107 ADDED, n_c, f107_i ", num_close, f107_ind
!   endif


! i don't think we need any of this but i'm unsure how the F10.7 code
! affects this.

!   do k = 1, num_close
!
!      t_ind = close_ind(k)
!      local_obs_loc   = obs_loc(t_ind)
!      local_obs_which = nint(query_location(local_obs_loc))
!
!      ! Convert vertical coordinate to requested vertical coordinate if necessary.
!      ! Should only be necessary for obs priors, as state location information already
!      ! contains the correct vertical coordinate
!      ! (filter_assim's call to get_state_meta_data).
!      if (vertical_localization_on()) then
! !fixme       if (local_obs_which /= wrf%dom(1)%vert_coord) then
! !fixme           call vert_interpolate(ens_mean, local_obs_loc, obs_kind(t_ind), istatus2)
!            ! Store the "new" location into the original full local array
!            obs_loc(t_ind) = local_obs_loc
! !fixme        endif
!      endif
!
!      ! Compute distance - set distance to a very large value if vert coordinate is
!      ! missing or vert_interpolate returned error (istatus2=1)
!      local_obs_array = get_location(local_obs_loc)
!
!
!      ! TJH FIXME ... the pbs_file script actually modifies the value of dist in
!      ! TJH FIXME ... this file and recompiles. NOT APPROPRIATE.
!      if (((vertical_localization_on())        .and. &
!           (local_obs_array(3) == MISSING_R8)) .or.  &
!           (istatus2 == 1)                   ) then
!         dist(k) = 1.0e9_r8
!      else
!         !if (close_ind(k) .eq. f107_ind) then !check if we came across the parameter
!         !   dist(k) = 0 !changed by pbs_file script
!         !else
!            dist(k) = get_dist(base_obs_loc, local_obs_loc, base_obs_type, obs_kind(t_ind))
!         !endif
!      endif
!   enddo

end subroutine get_close_state


!==================================================================
! The remaining PUBLIC interfaces come next
!==================================================================


subroutine get_quad_vals(state_handle, ens_size, varid, four_lons, four_lats, &
                         lon_lat_vert, which_vert, quad_vals, status_array)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: varid
integer,             intent(in)  :: four_lons(4), four_lats(4)
real(r8),            intent(in)  :: lon_lat_vert(3)
integer,             intent(in)  :: which_vert
real(r8),            intent(out) :: quad_vals(4, ens_size)
integer,             intent(out) :: status_array(ens_size)

real(r8) :: vert_val
integer  :: lev1, lev2, stat, integer_level
real(r8) :: vert_fract

character(len=*), parameter :: routine = 'get_quad_vals:'

quad_vals(:,:) = MISSING_R8
status_array(:) = GENERAL_ERROR_CODE

vert_val = lon_lat_vert(3)

select case (which_vert)

   case(VERTISHEIGHT)
      call find_enclosing_indices(NgridAlt, ALT(:), vert_val, lev1, lev2, &
         vert_fract, stat, log_scale = .false.)

      if (stat /= 0) then
         status_array     = INVALID_ALTITUDE_VAL_ERROR_CODE
      end if
      
   case(VERTISLEVEL)
      if (vert_val < 1.0_r8 .or. vert_val > NgridAlt) then
         status_array(:) = INVALID_ALTITUDE_VAL_ERROR_CODE
         return
      else
         integer_level = floor(vert_val)
         if (abs(vert_val - NgridAlt) > 1d-6) then
            lev1 = integer_level
            lev2 = lev1 + 1
            vert_fract = vert_val - lev1
         else
            lev2 = NgridAlt
            lev1 = NgridAlt - 1 
            vert_fract = 1.0_r8
         end if
      end if

      ! because we're given a model level as input, all the ensemble
      ! members have the same outgoing values.
      status_array     = 0

   case default
      status_array(:) = INVALID_VERT_COORD_ERROR_CODE
      write(string1, *) 'unsupported vertical type: ', which_vert
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
end select

! we have all the indices and fractions we could ever want.
! now get the data values at the bottom levels, the top levels, 
! and do vertical interpolation to get the 4 values in the columns.
! the final horizontal interpolation will happen later.
   
if (varid > 0) then

   call get_four_state_values(state_handle, ens_size, four_lons, four_lats, &
                              lev1, lev2, vert_fract, varid, quad_vals, &
                              status_array)

else 
      write(string1, *) 'unsupported variable: ', varid
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
endif

if (any(status_array /= 0)) return

! when you get here, status_array() was set either by passing it to a
! subroutine, or setting it explicitly here.
end subroutine get_quad_vals

!-----------------------------------------------------------------------
!> interpolate in the vertical between 2 arrays of items.
!>
!> vert_fracts: 0 is 100% of the first level and 
!>              1 is 100% of the second level

subroutine vert_interp(nitems, levs1, levs2, vert_fract, out_vals)

integer,  intent(in)  :: nitems
real(r8), intent(in)  :: levs1(nitems)
real(r8), intent(in)  :: levs2(nitems)
real(r8), intent(in)  :: vert_fract
real(r8), intent(out) :: out_vals(nitems)

out_vals(:) = (levs1(:) * (1.0_r8-vert_fract)) + &
              (levs2(:) *         vert_fract)

end subroutine vert_interp


!------------------------------------------------------------------
!> Converts gitm restart files to a netCDF file
!>
!> This routine needs:
!>
!> 1.  A base dirname for the restart files (restart_dirname).
!> they will have the format 'dirname/bNNNN.rst'  where NNNN has
!> leading 0s and is the block number.   Blocks start in the
!> southwest corner of the lat/lon grid and go west first, then
!> north and end in the northeast corner. The other info is in 
!> 'dirname/header.rst'
!> 
!> 2.  The name of the output file to store the netCDF variables 
!> (netcdf_output_file)
!>
!> In the process, the routine will find:
!>
!> 1. The overall grid size, lon/lat/alt when you've read in all
!>    the blocks.  (nGridLon, nGridLat, nGridAlt)
!>
!> 2. The number of blocks in Lon and Lat (nBlocksLon, nBlocksLat)
!>
!> 3. The number of lon/lats in a single grid block  (nLonsPerBlock, 
!>    nLatsPerBlock, nAltsPerBlock)
!>
!> 4. The number of neutral species (and probably a mapping between
!>    the species number and the variable name)  (nSpeciesTotal, nSpecies)
!>
!> 5. The number of ion species (ditto - numbers <-> names) (nIons)
!>
!> We assume that the 'UseTopography' flag is false - that all columns
!> have the same altitude arrays.  This is true on earth but not on
!> other planets.
!>
!> In addition to reading in the state data, it fills Longitude,
!> Latitude, and Altitude arrays with the grid spacing.  This grid
!> is orthogonal and rectangular but can have irregular spacing along
!> any or all of the three dimensions.

subroutine restart_files_to_netcdf(restart_dirname,netcdf_output_file)

character(len=*), intent(in)  :: restart_dirname
character(len=*), intent(in)  :: netcdf_output_file

integer :: ncid

character(len=*), parameter :: routine = 'restart_files_to_netcdf'

if (module_initialized ) then
    write(string1,*)'The gitm mod was already initialized but ',trim(routine),&
      ' uses a separate initialization procedure'
    call error_handler(E_ERR,routine,string1,source,revision,revdate)
end if

call static_init_blocks(restart_dirname)

ncid = nc_create_file(netcdf_output_file)

call add_nc_definitions(ncid)

call get_data(restart_dirname, ncid, define=.true.)

call nc_end_define_mode(ncid)

call add_nc_dimvars(ncid)

call get_data(restart_dirname, ncid, define=.false.)

call print_time(model_time)

call write_model_time(ncid, model_time)

call nc_close_file(ncid)

end subroutine restart_files_to_netcdf


!------------------------------------------------------------------

subroutine static_init_blocks(restart_dirname)

character(len=*), intent(in)  :: restart_dirname

character(len=*), parameter :: routine = 'static_init_blocks'

character(len=NF90_MAX_NAME)    :: varname
integer :: iunit, io, ivar
!logical :: has_gitm_namelist

call read_model_namelist()

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'gitm_blocks_nml', iunit)
read(iunit, nml = gitm_blocks_nml, iostat = io)
call check_namelist_read(iunit, io, 'gitm_blocks_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=gitm_blocks_nml)
if (do_nml_term()) write(     *     , nml=gitm_blocks_nml)

! Get the GITM variables in a restricted scope setting.

nSpecies      = get_nSpecies()
nSpeciesTotal = get_nSpeciesTotal()
nIons         = get_nIons()
nSpeciesAll   = get_nSpeciesAll()
nLonsPerBlock = get_nLonsPerBlock()
nLatsPerBlock = get_nLatsPerBlock()
nAltsPerBlock = get_nAltsPerBlock()

!---------------------------------------------------------------
! Set the time step ... causes gitm namelists to be read.
! Ensures model_advance_time is multiple of 'dynamics_timestep'

call set_calendar_type( calendar )   ! comes from model_mod_nml

!---------------------------------------------------------------
! 1) get grid dimensions
! 2) allocate space for the grids
! 3) read them from the block restart files, could be stretched ...

call get_grid_info_from_blocks(restart_dirname, NgridLon, NgridLat, NgridAlt, nBlocksLon, &
               nBlocksLat, LatStart, LatEnd, LonStart)

if( debug  > 0 ) then
    write(string1,*) 'grid dims are ',NgridLon,NgridLat,NgridAlt
    call error_handler(E_MSG,routine,string1,source,revision,revdate)
endif

allocate( LON( NgridLon ))
allocate( LAT( NgridLat ))
allocate( ALT( NgridAlt ))

call get_grid_from_blocks(restart_dirname, nBlocksLon, nBlocksLat, &
   nLonsPerBlock, nLatsPerBlock, nAltsPerBlock, LON, LAT, ALT )

! this is going to have to loop over all the blocks, both to get
! the data values and to get the full grid spacings.

model_time = get_state_time(restart_dirname)

if (do_output()) &
    call print_time(model_time,'time in restart file '//trim(restart_dirname)//'/header.rst')
if (do_output()) &
    call print_date(model_time,'date in restart file '//trim(restart_dirname)//'/header.rst')

call verify_block_variables( gitm_block_variables, nfields )

do ivar = 1, nfields

   varname                   = trim(gitm_block_variables(ivar))
   gitmvar(ivar)%varname     = varname

   ! This routine also checks to make sure user specified accurate GITM variables
   call decode_gitm_indices( varname,                    &
                             gitmvar(ivar)%gitm_varname, &
                             gitmvar(ivar)%gitm_dim,     &
                             gitmvar(ivar)%gitm_index,   &
                             gitmvar(ivar)%long_name,    &
                             gitmvar(ivar)%units)
   if ( debug > 0 ) then
      call print_gitmvar_info(ivar,routine)
   endif
enddo

if ( debug > 0 ) then
  write(string1,'("grid: NgridLon, NgridLat, NgridAlt =",3(1x,i5))') NgridLon, NgridLat, NgridAlt
  call error_handler(E_MSG,routine,string1,source,revision,revdate)
endif

end subroutine static_init_blocks

!------------------------------------------------------------------

subroutine print_gitmvar_info(ivar,routine)

   integer,          intent(in) :: ivar
   character(len=*), intent(in) :: routine

   call error_handler(E_MSG,routine,'',source,revision,revdate)
   write(string1,*) trim(gitmvar(ivar)%varname),' variable number ',ivar
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) ' long_name    ',trim(gitmvar(ivar)%long_name)
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) ' units        ',trim(gitmvar(ivar)%units)
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) ' gitm_varname ',trim(gitmvar(ivar)%gitm_varname)
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) ' gitm_dim     ',gitmvar(ivar)%gitm_dim
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) ' gitm_index   ',gitmvar(ivar)%gitm_index
   call error_handler(E_MSG,routine,string1,source,revision,revdate)

end subroutine


!------------------------------------------------------------------
! Writes the current time and state variables from a dart state
! vector (1d array) into a gitm netcdf restart file.

subroutine netcdf_to_restart_files(nc_file, output_dirname, input_dirname)

character(len=*), intent(in) :: nc_file
character(len=*), intent(in) :: output_dirname
character(len=*), intent(in) :: input_dirname

integer :: ncid

character(len=*), parameter :: routine = 'netcdf_to_restart_files:'

! sort the required fields into the order they exist in the
! binary restart files and write out the state vector data
! field by field.  when this routine returns all the data has
! been written.

if (module_initialized ) then
    write(string1,*)'The gitm mod was already initialized but ',trim(routine),&
      ' uses a separate initialization procedure'
    call error_handler(E_ERR,routine,string1,source,revision,revdate)
end if

call static_init_blocks(input_dirname)

ncid = nc_open_file_readonly(nc_file, routine)

call put_data(input_dirname, output_dirname, ncid)

end subroutine netcdf_to_restart_files


!------------------------------------------------------------------
! the static_init_model ensures that the gitm namelists are read.
!

function get_state_time( dirname )
type(time_type)              :: get_state_time
character(len=*), intent(in) :: dirname

type(time_type) :: model_offset, base_time

integer  :: iunit, i, ios
integer  :: istep
real(r8) :: tsimulation
integer  :: iyear, imonth, iday, ihour, imin, isec
integer  :: ndays,nsec

character(len=256) :: filename
character(len=100) :: cLine

character(len=*), parameter :: routine = 'get_state_time'

tsimulation = MISSING_R8
istep       = -1
iyear       = -1
imonth      = -1
iday        = -1
ihour       = -1
imin        = -1
isec        = -1

write(filename,'(a,''/header.rst'')') trim(dirname)

iunit = open_file(trim(filename), action='read')

FILEREAD : do i = 1, 100

   read(iunit,'(a)',iostat=ios) cLine

   if (ios < 0) exit FILEREAD  ! end of file

   if (ios /= 0) then
      write(string1,*) 'cannot read ',trim(filename)
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   endif

   select case( cLine(1:6) )
      case('#ISTEP')
         read(iunit,*)istep
      case('#TSIMU')
         read(iunit,*)tsimulation
      case('#TIMES')
         read(iunit,*)iyear
         read(iunit,*)imonth
         read(iunit,*)iday
         read(iunit,*)ihour
         read(iunit,*)imin
         read(iunit,*)isec
      case default
   end select

enddo FILEREAD

call close_file(iunit)

base_time      = set_date(iyear, imonth, iday, ihour, imin, isec)
ndays          = tsimulation/86400
nsec           = tsimulation - ndays*86400
model_offset   = set_time(nsec,ndays)
get_state_time = base_time + model_offset

if (debug > 8) then
   write(string1,*)'iyear       ',iyear
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'imonth      ',imonth
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'iday        ',iday
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ihour       ',ihour
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'imin        ',imin
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'isec        ',isec
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'tsimulation ',tsimulation
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ndays       ',ndays
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'nsec        ',nsec
   call error_handler(E_MSG,routine,string1,source,revision,revdate)

   call print_date(     base_time, 'get_state_time:model base date')
   call print_time(     base_time, 'get_state_time:model base time')
   call print_time(  model_offset, 'get_state_time:model offset')
   call print_date(get_state_time, 'get_state_time:model date')
   call print_time(get_state_time, 'get_state_time:model time')
endif

end function get_state_time


!==================================================================
! The remaining private interfaces come last
!==================================================================


!-----------------------------------------------------------------------
!>
!> Fill the array of requested variables, dart kinds, possible min/max
!> values and whether or not to update the field in the output file.
!> Then calls 'add_domain()' to tell the DART code which variables to
!> read into the state vector after this code returns.
!>
!>@param variable_array  the list of variables and kinds from model_mod_nml
!>@param nfields         the number of variable/Quantity pairs specified

subroutine set_gitm_variable_info(variable_array)

character(len=*), intent(in)  :: variable_array(:)

character(len=*), parameter :: routine = 'set_gitm_variable_info:'

integer :: i, nfields
integer, parameter :: MAX_STRING_LEN = 128

character(len=MAX_STRING_LEN) :: varname    ! column 1, NetCDF variable name
character(len=MAX_STRING_LEN) :: dartstr    ! column 2, DART Quantity
character(len=MAX_STRING_LEN) :: minvalstr  ! column 3, Clamp min val
character(len=MAX_STRING_LEN) :: maxvalstr  ! column 4, Clamp max val
character(len=MAX_STRING_LEN) :: updatestr  ! column 5, Update output or not

integer, parameter :: MAX_STATE_VARIABLES = 100
integer, parameter :: vtablenamelength = 32

character(len=vtablenamelength) :: var_names(MAX_STATE_VARIABLES)
logical  :: update_list(MAX_STATE_VARIABLES)   = .FALSE.
integer  ::   kind_list(MAX_STATE_VARIABLES)   = MISSING_I
real(r8) ::  clamp_vals(MAX_STATE_VARIABLES,2) = MISSING_R8

var_names(:) = ' ' 

nfields = 0
ParseVariables : do i = 1, MAX_STATE_VARIABLES

   varname   = variable_array(num_state_table_columns*i-4)
   dartstr   = variable_array(num_state_table_columns*i-3)
   minvalstr = variable_array(num_state_table_columns*i-2)
   maxvalstr = variable_array(num_state_table_columns*i-1)
   updatestr = variable_array(num_state_table_columns*i  )

   if ( varname == ' ' .and. dartstr == ' ' ) exit ParseVariables ! Found end of list.

   if ( varname == ' ' .or.  dartstr == ' ' ) then
      string1 = 'model_nml:model "state_variables" not fully specified'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   endif

   ! Make sure DART kind is valid

   if( get_index_for_quantity(dartstr) < 0 ) then
      write(string1,'(3A)') 'there is no obs_kind "', trim(dartstr), '" in obs_kind_mod.f90'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   endif

   call to_upper(minvalstr)
   call to_upper(maxvalstr)
   call to_upper(updatestr)

   var_names(   i) = varname
   kind_list(   i) = get_index_for_quantity(dartstr)
   clamp_vals(i,1) = string_to_real(minvalstr)
   clamp_vals(i,2) = string_to_real(maxvalstr)
   update_list( i) = string_to_logical(updatestr, 'UPDATE')

   nfields = nfields + 1

enddo ParseVariables

if (nfields == MAX_STATE_VARIABLES) then
   write(string1,'(2A)') 'WARNING: There is a possibility you need to increase ', &
                         'MAX_STATE_VARIABLES in the global variables in model_mod.f90'

   write(string2,'(A,i4,A)') 'WARNING: you have specified at least ', nfields, &
                             ' perhaps more'

   call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)
endif

! gitm only has a single domain (only a single grid, no nests or multiple grids)

domain_id = add_domain(template_filename, nfields, var_names, kind_list, &
                       clamp_vals, update_list)
!domain_id = add_domain(nfields, var_names, kind_list, &
!                       clamp_vals, update_list)

if (debug > 1) call state_structure_info(domain_id)

end subroutine set_gitm_variable_info


!-----------------------------------------------------------------------
!>

subroutine get_four_state_values(state_handle, ens_size, four_lons, four_lats, &
                                 lev1, lev2, vert_fract, varid, quad_vals, &
                                 my_status)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: four_lons(4), four_lats(4)
integer,             intent(in) :: lev1, lev2
real(r8),            intent(in) :: vert_fract
integer,             intent(in) :: varid
real(r8),           intent(out) :: quad_vals(4, ens_size) !< array of interpolated values
integer,            intent(out) :: my_status(ens_size)

integer  :: icorner
real(r8) :: vals1(ens_size), vals2(ens_size)
real(r8) :: qvals(ens_size)

integer(i8) :: state_indx

character(len=*), parameter :: routine = 'get_four_state_values:'

do icorner=1, 4

   state_indx = get_dart_vector_index(four_lons(icorner), four_lats(icorner), &
         lev1, domain_id, varid)

   if (state_indx < 0) then
      write(string1,*) 'Could not find dart state index from '
      write(string2,*) 'lon, lat, and lev index :', four_lons(icorner), four_lats(icorner), lev2
      call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2)
      return
   endif

   vals1(:) = get_state(state_indx, state_handle)    ! all the ensemble members for level (i)

   state_indx = get_dart_vector_index(four_lons(icorner), four_lats(icorner), &
         lev2, domain_id, varid)

   if (state_indx < 0) then
      write(string1,*) 'Could not find dart state index from '
      write(string2,*) 'lon, lat, and lev index :', four_lons(icorner), four_lats(icorner), lev2
      call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2)
      return
   endif

   vals2(:) = get_state(state_indx, state_handle)    ! all the ensemble members for level (i)

   ! if directly using quad_vals here, it would create a temporary array and give a warning
   call vert_interp(ens_size, vals1, vals2, vert_fract, qvals)
   quad_vals(icorner, :) = qvals
enddo

my_status = 0

end subroutine get_four_state_values


!==================================================================


subroutine get_grid_info_from_netcdf(template_filename, nLon, nLat, nAlt )

character(len=*), intent(in)  :: template_filename
integer,          intent(out) :: nlat, nlon, nalt

character(len=*), parameter :: routine = 'get_grid_info_from_netcdf'

integer :: ncid

ncid = nc_open_file_readonly(template_filename, routine)

nLat = nc_get_dimension_size(ncid, LAT_DIM_NAME, routine)
nLon = nc_get_dimension_size(ncid, LON_DIM_NAME, routine)
nAlt = nc_get_dimension_size(ncid, ALT_DIM_NAME, routine)

call nc_close_file(ncid)

end subroutine


!==================================================================

! Read the lon, lat, and alt arrays from the ncid

subroutine get_grid_from_netcdf(template_filename, LON, LAT, ALT )

character(len=*), intent(in)    :: template_filename
real(r8),         intent(inout) :: LON(:)
real(r8),         intent(inout) :: LAT(:)
real(r8),         intent(inout) :: ALT(:)

character(len=*), parameter :: routine = 'get_grid_from_netcdf'

integer :: ncid

ncid = nc_open_file_readonly(template_filename, routine)

call nc_get_variable(ncid, LAT_VAR_NAME, LAT, routine)
call nc_get_variable(ncid, LON_VAR_NAME, LON, routine)
call nc_get_variable(ncid, ALT_VAR_NAME, ALT, routine)

call nc_close_file(ncid)

end subroutine get_grid_from_netcdf


!==================================================================

!> Read the grid dimensions from the restart netcdf file.
!>
!> The file name comes from module storage ... namelist.

subroutine get_grid_info_from_blocks(gitm_restart_dirname, NgridLon, NgridLat, &
                NgridAlt, nBlocksLon, nBlocksLat, LatStart, LatEnd, LonStart)

character(len=*), intent(in) :: gitm_restart_dirname
integer,  intent(out) :: NgridLon   ! Number of Longitude centers
integer,  intent(out) :: NgridLat   ! Number of Latitude  centers
integer,  intent(out) :: NgridAlt   ! Number of Vertical grid centers
integer,  intent(out) :: nBlocksLon, nBlocksLat
real(r8), intent(out) :: LatStart, LatEnd, LonStart

character(len=*), parameter :: filename = 'UAM.in'

character(len=100) :: cLine  ! iCharLen_ == 100
character(len=256) :: fileloc

integer :: i, iunit, ios

character(len=*), parameter :: routine = 'get_grid_info_from_blocks'

! get the ball rolling ...

nBlocksLon = 0
nBlocksLat = 0
LatStart   = 0.0_r8
LatEnd     = 0.0_r8
LonStart   = 0.0_r8

write(fileloc,'(a,''/'',a)') trim(gitm_restart_dirname),trim(filename)

if (debug > 4) then
   write(string1,*) 'Now opening GITM restart file: ',trim(fileloc)
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
end if


iunit = open_file(trim(fileloc), action='read')

UAMREAD : do i = 1, 1000000

   read(iunit,'(a)',iostat=ios) cLine

   if (ios /= 0) then
      ! If we get to the end of the file or hit a read error without
      ! finding what we need, die.
      write(string1,*) 'cannot find #GRID in ',trim(fileloc)
      call error_handler(E_ERR,'get_grid_info_from_blocks',string1,source,revision,revdate)
   endif

   if (cLine(1:5) .ne. "#GRID") cycle UAMREAD

   nBlocksLon = read_in_int( iunit,'NBlocksLon',trim(fileloc))
   nBlocksLat = read_in_int( iunit,'NBlocksLat',trim(fileloc))
   LatStart   = read_in_real(iunit,'LatStart',  trim(fileloc))
   LatEnd     = read_in_real(iunit,'LatEnd',    trim(fileloc))
   LonStart   = read_in_real(iunit,'LonStart',  trim(fileloc))

   exit UAMREAD

enddo UAMREAD

if (debug > 4) then
   write(string1,*) 'Successfully read GITM restart file:',trim(fileloc)
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   nLonsPerBlock:',nLonsPerBlock
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   nLatsPerBlock:',nLatsPerBlock
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   nBlocksLon:',nBlocksLon
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   nBlocksLat:',nBlocksLat
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   LatStart:',LatStart
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   LatEnd:',LatEnd
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   LonStart:',LonStart
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
end if

call close_file(iunit)

NgridLon = nBlocksLon * nLonsPerBlock
NgridLat = nBlocksLat * nLatsPerBlock
NgridAlt = nAltsPerBlock

write(string1,*)  'NgridLon = ', NgridLon
call error_handler(E_MSG,routine,string1,source,revision,revdate)
write(string1,*)  'NgridLat = ', NgridLat
call error_handler(E_MSG,routine,string1,source,revision,revdate)
write(string1,*)  'NgridAlt = ', NgridAlt
call error_handler(E_MSG,routine,string1,source,revision,revdate)

end subroutine get_grid_info_from_blocks


!==================================================================

! open enough of the restart files to read in the lon, lat, alt arrays

subroutine get_grid_from_blocks(dirname, nBlocksLon, nBlocksLat, &
                  nLonsPerBlock, nLatsPerBlock, nAltsPerBlock,   &
                  LON, LAT, ALT )

character(len=*), intent(in) :: dirname
integer, intent(in) :: nBlocksLon    ! Number of Longitude blocks
integer, intent(in) :: nBlocksLat    ! Number of Latitude  blocks
integer, intent(in) :: nLonsPerBlock ! Number of Longitude centers per block
integer, intent(in) :: nLatsPerBlock ! Number of Latitude  centers per block
integer, intent(in) :: nAltsPerBlock ! Number of Vertical grid centers

real(r8), dimension( : ), intent(inout) :: LON, LAT, ALT

integer :: ios, nb, offset, iunit, nboff
character(len=256) :: filename
real(r8), allocatable :: temp(:)

character(len=*), parameter :: routine = 'get_grid_from_blocks'

! a temp array large enough to hold any of the
! Lon,Lat or Alt array from a block plus ghost cells
allocate(temp(1-nGhost:max(nLonsPerBlock,nLatsPerBlock,nAltsPerBlock)+nGhost))

! go across the south-most block row picking up all longitudes
do nb = 1, nBlocksLon

   iunit = open_block_file(dirname, nb, 'read', filename)

   read(iunit,iostat=ios) temp(1-nGhost:nLonsPerBlock+nGhost)
   if ( ios /= 0 ) then
      print *,'size:',size(temp(1-nGhost:nLonsPerBlock+nGhost))
      print *,'IO error code:',ios
      write(string1,*)'ERROR reading file ', trim(filename)
      write(string2,*)'longitude block ',nb,' of ',nBlocksLon
      call error_handler(E_ERR,'get_grid',string1, &
                 source,revision,revdate,text2=string2)
   endif

   offset = (nLonsPerBlock * (nb - 1))
   LON(offset+1:offset+nLonsPerBlock) = temp(1:nLonsPerBlock)

   call close_file(iunit)
enddo

! go up west-most block row picking up all latitudes
do nb = 1, nBlocksLat

   nboff = ((nb - 1) * nBlocksLon) + 1
   iunit = open_block_file(dirname, nboff, 'read', filename)

   ! get past lon array and read in lats
   read(iunit) temp(1-nGhost:nLonsPerBlock+nGhost)

   read(iunit,iostat=ios) temp(1-nGhost:nLatsPerBlock+nGhost)
   if ( ios /= 0 ) then
      write(string1,*)'ERROR reading file ', trim(filename)
      write(string2,*)'latitude block ',nb,' of ',nBlocksLat
      call error_handler(E_ERR,'get_grid',string1, &
                 source,revision,revdate,text2=string2)
   endif

   offset = (nLatsPerBlock * (nb - 1))
   LAT(offset+1:offset+nLatsPerBlock) = temp(1:nLatsPerBlock)

   call close_file(iunit)
enddo

! this code assumes UseTopography is false - that all columns share
! the same altitude array, so we can read it from the first block.
! if this is not the case, this code has to change.

iunit = open_block_file(dirname, 1, 'read', filename)

! get past lon and lat arrays and read in alt array
read(iunit) temp(1-nGhost:nLonsPerBlock+nGhost)
read(iunit) temp(1-nGhost:nLatsPerBlock+nGhost)
read(iunit) temp(1-nGhost:nAltsPerBlock+nGhost)

ALT(1:nAltsPerBlock) = temp(1:nAltsPerBlock)

call close_file(iunit)

deallocate(temp)

! convert from radians into degrees
LON = LON * rad2deg
LAT = LAT * rad2deg

if (debug > 4) then
   print *, 'All LONs ', LON
   print *, 'All LATs ', LAT
   print *, 'All ALTs ', ALT
endif

if ( debug > 1 ) then ! A little sanity check
   write(string1,*)'LON range ',minval(LON),maxval(LON)
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'LAT range ',minval(LAT),maxval(LAT)
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ALT range ',minval(ALT),maxval(ALT)
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
endif

end subroutine get_grid_from_blocks



!==================================================================

!> open the requested block number restart file and return the
!> file unit

function open_block_file(dirname, blocknum, rw, filename)

integer                       :: open_block_file
character(len=*), intent(in)  :: dirname
integer,          intent(in)  :: blocknum
character(len=*), intent(in)  :: rw   ! 'read' or 'readwrite'
character(len=*), intent(out) :: filename

write(filename, '(A,i4.4,A)') trim(dirname)//'/b', blocknum, '.rst'

if ( rw == 'read' .and. .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'open_block_file',string1,source,revision,revdate)
endif

if (debug > 0) then
   write(string1,*) 'Opening file ', trim(filename), ' for ', trim(rw)
   call error_handler(E_MSG,'open_block_file',string1,source,revision,revdate)
end if

open_block_file = open_file(filename, 'unformatted', rw)

if (debug > 80) then
   write(string1,*) 'Returned file descriptor is ', open_block_file
   call error_handler(E_MSG,'open_block_file',string1,source,revision,revdate)
end if

end function open_block_file


!------------------------------------------------------------------
! open all restart files and read in the requested data item

subroutine get_data(dirname, ncid, define)

character(len=*), intent(in)  :: dirname
integer,          intent(in)  :: ncid
logical,          intent(in)  :: define

integer :: ibLoop, jbLoop
integer :: ib, jb, nb, iunit

character(len=256) :: filename

! get the dirname, construct the filenames inside open_block_file

if (define) then
   ! if define, run one block.
   ! the read_data_from_block call defines the variables in the netCDF file.
   ibLoop = 1
   jbLoop = 1
else
   ! if not define, run all blocks.
   ! the read_data_from_block call adds the (ib,jb) block to a netCDF variable 
   ! in order to make a file containing the data for all the blocks.
   ibLoop = nBlocksLon
   jbLoop = nBlocksLat
end if

do jb = 1, jbLoop
   do ib = 1, ibLoop
      nb = (jb-1) * nBlocksLon + ib

      iunit = open_block_file(dirname, nb, 'read', filename)

      call read_data_from_block(iunit, ib, jb, ncid, define)

      call close_file(iunit)
   enddo
enddo

end subroutine get_data


!==================================================================

! open all restart files and write out the requested data item

subroutine put_data(dirname, dirnameout, ncid)

character(len=*), intent(in) :: dirname, dirnameout
integer,          intent(in) :: ncid

integer :: ib, jb, nb, iunit, ounit
character(len=256) :: readfilename, writefilename

character(len=*), parameter :: routine = 'put_data'

! get the dirname, construct the filenames inside open_block_file

if (debug > 0) then
   write(string1,'(A,I0,A,I0,A)') 'Now putting the data for ',nBlocksLon,' blocks lon by ',nBlocksLat,' blocks lat'
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
end if



do jb = 1, nBlocksLat
   do ib = 1, nBlocksLon

      nb = (jb-1) * nBlocksLon + ib

      iunit = open_block_file(dirname,    nb, 'read',  readfilename)
      ounit = open_block_file(dirnameout, nb, 'write', writefilename)

      if (readfilename == writefilename) then
         write(string1,*) 'Cannot use the same input and output restart files: ',trim(readfilename)
         write(string2,*) 'Please specify different input/output directories in the input.nml file'
         call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2)
      end if

      call write_data(iunit, ounit, ib, jb,  ncid, readfilename, writefilename)

      call close_file(iunit)
      call close_file(ounit)
   enddo
enddo

end subroutine put_data


!==================================================================

! put the requested data into a netcdf variable

subroutine unpack_data(data3d, ivar, block, ncid, define)

real(r8), intent(in)    :: data3d(1-nGhost:nLonsPerBlock+nGhost, &
                                  1-nGhost:nLatsPerBlock+nGhost, &
                                  1-nGhost:nAltsPerBlock+nGhost)

integer,  intent(in)    :: ivar         ! variable index
integer,  intent(in)    :: block(2)
integer,  intent(in)    :: ncid
logical,  intent(in)    :: define

integer :: ib, jb
integer :: starts(3)
character(len=*), parameter :: routine = 'unpack_data'

if (define) then
  
   if (debug > 10) then 
      write(string1,'(A,I0,2A)') 'Defining ivar = ', ivar,':',trim(gitmvar(ivar)%varname)
      call error_handler(E_MSG,routine,string1,source,revision,revdate)
   end if

   call nc_define_double_variable(ncid, gitmvar(ivar)%varname, (/ LON_DIM_NAME, LAT_DIM_NAME, ALT_DIM_NAME /) )
   call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'long_name',      gitmvar(ivar)%long_name)
   call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'units',          gitmvar(ivar)%units)
   !call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'storder',        gitmvar(ivar)%storder)
   call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'gitm_varname',   gitmvar(ivar)%gitm_varname)
   call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'gitm_dim',       gitmvar(ivar)%gitm_dim)
   call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'gitm_index',     gitmvar(ivar)%gitm_index)

else

   ib = block(1)
   jb = block(2)

   ! to compute the start, consider (ib-1)*nLonsPerBlock+1
   starts(1) = (ib-1)*nLonsPerBlock+1
   starts(2) = (jb-1)*nLatsPerBlock+1
   starts(3) = 1

   call nc_put_variable(ncid, gitmvar(ivar)%varname, &
      data3d(1:nLonsPerBlock,1:nLatsPerBlock,1:nAltsPerBlock), &
      context=routine, nc_start=starts, &
      nc_count=(/nLonsPerBlock,nLatsPerBlock,nAltsPerBlock/))
end if

end subroutine unpack_data


!==================================================================

! put the requested data into a netcdf variable

subroutine unpack_data2d(data2d, ivar, block, ncid, define)

real(r8), intent(in)    :: data2d(1-nGhost:nLonsPerBlock+nGhost, &
                                  1-nGhost:nLatsPerBlock+nGhost)

integer,  intent(in)    :: ivar         ! variable index
integer,  intent(in)    :: block(2)
integer,  intent(in)    :: ncid
logical,  intent(in)    :: define

integer :: ib, jb
integer :: starts(2)
character(len=*), parameter :: routine = 'unpack_data2d'

if (define) then
  
   if (debug > 10) then 
      write(string1,'(A,I0,2A)') 'Defining ivar = ', ivar,':',trim(gitmvar(ivar)%varname)
      call error_handler(E_MSG,routine,string1,source,revision,revdate)
   end if

   call nc_define_double_variable(ncid, gitmvar(ivar)%varname, (/ LON_DIM_NAME, LAT_DIM_NAME /) )
   call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'long_name',      gitmvar(ivar)%long_name)
   call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'units',          gitmvar(ivar)%units)
   !call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'storder',        gitmvar(ivar)%storder)
   call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'gitm_varname',   gitmvar(ivar)%gitm_varname)
   call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'gitm_dim',       gitmvar(ivar)%gitm_dim)
   call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'gitm_index',     gitmvar(ivar)%gitm_index)

else
   ib = block(1)
   jb = block(2)

   ! to compute the start, consider (ib-1)*nLonsPerBlock+1
   starts(1) = (ib-1)*nLonsPerBlock+1
   starts(2) = (jb-1)*nLatsPerBlock+1

   call nc_put_variable(ncid, gitmvar(ivar)%varname, &
      data2d(1:nLonsPerBlock,1:nLatsPerBlock), &
      context=routine, nc_start=starts, &
      nc_count=(/nLonsPerBlock,nLatsPerBlock/))
end if

end subroutine unpack_data2d


!==================================================================


!> put the f107 estimate (a scalar, hence 0d) into the state vector.
!> Written specifically
!> for f107 since f107 is the same for all blocks. So what it does
!> is take f107 from the first block (block = 0) and disregard
!> f107 values from all other blocks (hopefully they are the same).
!> written by alex

subroutine unpack_data0d(data0d, ivar, ncid, define)

real(r8), intent(in)    :: data0d
integer,  intent(in)    :: ivar         ! index into state structure
integer,  intent(in)    :: ncid
logical,  intent(in)    :: define


character(len=*), parameter :: routine = 'unpack_data0d'

if (define) then
  
   if (debug > 10) then 
      write(string1,'(A,I0,2A)') 'Defining ivar = ', ivar,':',trim(gitmvar(ivar)%varname)
      call error_handler(E_MSG,routine,string1,source,revision,revdate)
   end if

   call nc_define_double_scalar(ncid,   gitmvar(ivar)%varname)
   call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'long_name',      gitmvar(ivar)%long_name)
   call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'units',          gitmvar(ivar)%units)
   call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'gitm_varname',   gitmvar(ivar)%gitm_varname)
   call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'gitm_dim',       gitmvar(ivar)%gitm_dim)
   call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'gitm_index',     gitmvar(ivar)%gitm_index)

else

   call nc_put_variable(ncid, gitmvar(ivar)%varname, data0d, context=routine)

end if

end subroutine unpack_data0d


!==================================================================

! put the state vector data into a 3d array

subroutine pack_data(ncid, ivar, block, data3d)

integer,  intent(in)    :: ncid
integer,  intent(in)    :: ivar         ! index into state structure
integer,  intent(in)    :: block(2)
real(r8), intent(inout) :: data3d(1-nghost:nLonsPerBlock+nghost,&
                                  1-nghost:nLatsPerblock+nghost,&
                                  1-nGhost:nAltsPerBlock+nghost)
integer :: ib, jb
integer :: starts(3)
integer :: ends(3)
integer :: local_starts(3)
integer :: local_ends(3)
integer :: counts(3)
integer :: maxvals(3)
integer :: i, j

character(len=*), parameter :: routine = 'pack_data'

ib = block(1)
jb = block(2)

! to compute the start, consider (ib-1)*nLonsPerBlock+1
starts(1) = (ib-1)*nLonsPerBlock + 1 - nghost
  ends(1) =     ib*nLonsPerBlock +     nghost
starts(2) = (jb-1)*nLatsPerBlock + 1 - nghost
  ends(2) =     jb*nLatsPerBlock +     nghost
starts(3)=       0               + 1 - nghost
  ends(3)=nAltsPerBlock              + nghost

maxvals = (/NgridLon,NgridLat,NgridAlt/)

do i=1,3
   if (starts(i) < 1) then
      starts(i)       = 1
      local_starts(i) = 1
   else
      local_starts(i) = 1-nghost
   end if
end do

do i=1,3
   if (ends(i) > maxvals(i)) then
      ends(i)         = maxvals(i)
   end if
end do

counts = ends-starts+1

local_ends = local_starts + counts - 1

if (debug > 10) then
   if (ivar == 1) then
      write(string1,'(12(A,I0),A)') 'Now reading netCDF indices (',starts(1),':',ends(1),') and (',starts(2),':',ends(2),') ' // &
         'in the block (',ib,',',jb,') for local (',local_starts(1),':',local_ends(1),') and (',&
         local_starts(2),':',local_ends(2),') for size (',counts(1),',',counts(2),')'
      call error_handler(E_MSG,routine,string1,source,revision,revdate)
   end if
end if

if (debug > 200) then
   if (ivar == 3) then
      print *,'before reading:'
      do i=1-nghost,nLonsPerBlock+nghost
         do j=1-nghost,nLatsPerBlock+nghost
            write(*,'(A,A,I0,A,I0,A,F8.3)') trim(gitmvar(ivar)%varname), &
               ' at level 3 (',i,',',j,'): ',data3D(i,j,3)
         end do
      end do
   end if
end if

call nc_get_variable(ncid, gitmvar(ivar)%varname, data3d(local_starts(1):local_ends(1),&
      local_starts(2):local_ends(2),local_starts(3):local_ends(3)), &
      context="pack_data", nc_start=starts, nc_count=counts)

if (debug > 200) then
   if (ivar == 3) then
      print *,'after reading:'
      do i=1-nghost,nLonsPerBlock+nghost
         do j=1-nghost,nLatsPerBlock+nghost
            write(*,'(A,A,I0,A,I0,A,F8.3)') trim(gitmvar(ivar)%varname), &
               ' at level 3 (',i,',',j,'): ',data3D(i,j,3)
         end do
      end do
   end if
end if

end subroutine pack_data


!==================================================================

!> put the f107 estimate (scalar) from the statevector into a 0d container
!> the only trick this routine does is give all blocks the same f107 (the
!> f107 value from block 1 state vector goes to block 1,2,3,4 restart files)
!> so no matter what, always grab the f107 from block 1 (manipulate
!> the block variable).
!> written by alex

subroutine pack_data0d(ncid, ivar, data0d)

integer,  intent(in)    :: ncid
integer,  intent(in)    :: ivar         ! index into state structure
real(r8), intent(inout) :: data0d

call nc_get_variable(ncid, gitmvar(ivar)%varname, data0d,&
   context="pack_data0d")

end subroutine pack_data0d


!==================================================================

!> open all restart files and read in the requested data items
!>
!> This is a two-pass method: first run through to define the NC variables
!> (define = .true.), then run again to write the data to the NC file
!> (define = .false.)

subroutine read_data_from_block(iunit, ib, jb, ncid, define)

integer,  intent(in) :: iunit
integer,  intent(in) :: ib, jb
integer,  intent(in) :: ncid
logical,  intent(in) :: define

real(r8), allocatable :: temp1d(:), temp2d(:,:), temp3d(:,:,:), temp4d(:,:,:,:)
real(r8), allocatable :: alt1d(:), density_ion_e(:,:,:)
real(r8) :: temp0d !Alex: single parameter has "zero dimensions"
integer :: i, j, inum, maxsize, ivals(NSpeciesTotal)
integer :: block(2) = 0

logical :: no_idensity

character(len=*), parameter :: routine = 'read_data_from_block'

block(1) = ib
block(2) = jb

! a temp array large enough to hold any of the
! Lon,Lat or Alt array from a block plus ghost cells
allocate(temp1d(1-nGhost:max(nLonsPerBlock,nLatsPerBlock,nAltsPerBlock)+nGhost))
! treat alt specially since we want to derive TEC here
allocate( alt1d(1-nGhost:max(nLonsPerBlock,nLatsPerBlock,nAltsPerBlock)+nGhost))

! temp array large enough to hold any 2D field 
allocate(temp2d(1-nGhost:nLonsPerBlock+nGhost, &
                1-nGhost:nLatsPerBlock+nGhost))

! temp array large enough to hold 1 species, temperature, etc
allocate(temp3d(1-nGhost:nLonsPerBlock+nGhost, &
                1-nGhost:nLatsPerBlock+nGhost, &
                1-nGhost:nAltsPerBlock+nGhost))

! save density_ion_e to compute TEC
allocate(density_ion_e(1-nGhost:nLonsPerBlock+nGhost, &
                       1-nGhost:nLatsPerBlock+nGhost, &
                       1-nGhost:nAltsPerBlock+nGhost))

! temp array large enough to hold velocity vect, etc
maxsize = max(3, nSpecies)
allocate(temp4d(1-nGhost:nLonsPerBlock+nGhost, &
                1-nGhost:nLatsPerBlock+nGhost, &
                1-nGhost:nAltsPerBlock+nGhost, maxsize))

! get past lon and lat arrays and read in alt array
read(iunit) temp1d(1-nGhost:nLonsPerBlock+nGhost)
read(iunit) temp1d(1-nGhost:nLatsPerBlock+nGhost)
! save the alt1d for later TEC computation
read(iunit)  alt1d(1-nGhost:nAltsPerBlock+nGhost)

! Read the index from the first species
call get_index_from_gitm_varname('NDensityS', inum, ivals)

if (inum > 0) then
   ! if i equals ival, use the data from the state vect
   ! otherwise read/write what's in the input file
   j = 1
   do i = 1, nSpeciesTotal
      if (debug > 80) then
         write(string1,'(A,I0,A,I0,A,I0,A,I0,A)') 'Now reading species ',i,' of ',nSpeciesTotal, &
            ' for block (',ib,',',jb,')' 
         call error_handler(E_MSG,routine,string1,source,revision,revdate)
      end if
      read(iunit)  temp3d
      if (j <= inum) then
         if (i == gitmvar(ivals(j))%gitm_index) then
            call unpack_data(temp3d, ivals(j), block, ncid, define)
            j = j + 1
         endif
      endif
   enddo
else
   if (debug > 80) then
      write(string1,'(A)') 'Not writing the NDensityS variables to file'
      call error_handler(E_MSG,routine,string1,source,revision,revdate)
   end if
   ! nothing at all from this variable in the state vector.
   ! copy all data over from the input file to output file
   do i = 1, nSpeciesTotal
      read(iunit)  temp3d
   enddo
endif

call get_index_from_gitm_varname('IDensityS', inum, ivals)

! assume we could not find the electron density for VTEC calculations
no_idensity = .true.

if (inum > 0) then
   ! one or more items in the state vector need to replace the
   ! data in the output file.  loop over the index list in order.
   j = 1
   do i = 1, nIons
      if (debug > 80) then
         write(string1,'(A,I0,A,I0,A,I0,A,I0,A)') 'Now reading ion ',i,' of ',nIons, &
            ' for block (',ib,',',jb,')' 
         call error_handler(E_MSG,routine,string1,source,revision,revdate)
      end if
      read(iunit)  temp3d
      if (j <= inum) then
         if (i == gitmvar(ivals(j))%gitm_index) then
            ! ie_, the gitm index for electron density, comes from ModEarth 
            if (gitmvar(ivals(j))%gitm_index == ie_) then
               ! save the electron density for TEC computation
               density_ion_e(:,:,:) = temp3d(:,:,:)
               no_idensity = .false.
            end if
            ! read from input but write from state vector
            call unpack_data(temp3d, ivals(j), block, ncid, define)
            j = j + 1
         endif
      endif
   enddo
else
   ! nothing at all from this variable in the state vector.
   ! read past this variable
   if (debug > 80) then
      write(string1,'(A)') 'Not writing the IDensityS variables to file'
      call error_handler(E_MSG,routine,string1,source,revision,revdate)
   end if
   do i = 1, nIons
      read(iunit)  temp3d
   enddo
endif

read(iunit)  temp3d
call get_index_from_gitm_varname('Temperature', inum, ivals)

if (inum > 0) then
   call unpack_data(temp3d, ivals(1), block, ncid, define)
endif

read(iunit) temp3d
call get_index_from_gitm_varname('ITemperature', inum, ivals)
if (inum > 0) then
   call unpack_data(temp3d, ivals(1), block, ncid, define)
endif

read(iunit) temp3d
call get_index_from_gitm_varname('eTemperature', inum, ivals)
if (inum > 0) then
   call unpack_data(temp3d, ivals(1), block, ncid, define)
endif

read(iunit) temp4d(:,:,:,1:3)
call get_index_from_gitm_varname('Velocity', inum, ivals)
if (inum > 0) then
   ! copy out any requested bits into state vector
   j = 1
   do i = 1, 3
      if (j <= inum) then
         if (i == gitmvar(ivals(j))%gitm_index) then
            temp3d = temp4d(:,:,:,i)
            call unpack_data(temp3d, ivals(j), block, ncid, define)
            j = j + 1
         endif
      endif
   enddo
endif

read(iunit) temp4d(:,:,:,1:3)
call get_index_from_gitm_varname('IVelocity', inum, ivals)
if (inum > 0) then
   ! copy out any requested bits into state vector
   j = 1
   do i = 1, 3
      if (j <= inum) then
         if (i == gitmvar(ivals(j))%gitm_index) then
            ! read from input but write from state vector
            temp3d = temp4d(:,:,:,i)
            call unpack_data(temp3d, ivals(j), block, ncid, define)
            j = j + 1
         endif
      endif
   enddo
endif

!print *, 'reading in temp4d for vvel'
read(iunit) temp4d(:,:,:,1:nSpecies)
call get_index_from_gitm_varname('VerticalVelocity', inum, ivals)
if (inum > 0) then
   ! copy out any requested bits into state vector
   j = 1
   do i = 1, nSpecies
      if (j <= inum) then
         if (i == gitmvar(ivals(j))%gitm_index) then
            temp3d = temp4d(:,:,:,i)
            call unpack_data(temp3d, ivals(j), block, ncid, define)
            j = j + 1
         endif
      endif
   enddo
endif

! add the VTEC as an extended-state variable
! NOTE: This variable will *not* be written out to the GITM blocks to netCDF program
call get_index_from_gitm_varname('TEC', inum, ivals)

if (inum > 0 .and. no_idensity) then
   write(string1,*) 'Cannot compute the VTEC without the electron density'
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
end if

if (inum > 0) then
   if (.not. define) then
      temp2d = 0._r8
      ! comptue the TEC integral
      do i =1,nAltsPerBlock-1 ! approximate the integral over the altitude as a sum of trapezoids
         ! area of a trapezoid: A = (h2-h1) * (f2+f1)/2
         temp2d(:,:) = temp2d(:,:) + ( alt1d(i+1)-alt1d(i) )  * ( density_ion_e(:,:,i+1)+density_ion_e(:,:,i) ) /2.0_r8
      end do  
      ! convert temp2d to TEC units
      temp2d = temp2d/1e16_r8
   end if
   call unpack_data2d(temp2d, ivals(1), block, ncid, define) 
end if

!alex begin
read(iunit)  temp0d
!gitm_index = get_index_start(domain_id, 'VerticalVelocity')
call get_index_from_gitm_varname('f107', inum, ivals)
if (inum > 0) then
  call unpack_data0d(temp0d, ivals(1), ncid, define) !see comments in the body of the subroutine
endif

read(iunit)  temp3d
call get_index_from_gitm_varname('Rho', inum, ivals)
if (inum > 0) then
   call unpack_data(temp3d, ivals(1), block, ncid, define)
endif
!alex end

!print *, 'calling dealloc'
deallocate(temp1d, temp2d, temp3d, temp4d)
deallocate(alt1d, density_ion_e)

end subroutine read_data_from_block


!==================================================================

! open all restart files and write out the requested data item

subroutine write_data(iunit, ounit, ib, jb, ncid, infile, outfile)

integer,          intent(in) :: iunit, ounit
integer,          intent(in) :: ib, jb, ncid
character(len=*), intent(in) :: infile, outfile

real(r8), allocatable :: temp1d(:), temp3d(:,:,:), temp4d(:,:,:,:), data3d(:,:,:)
real(r8) :: data0d, temp0d !Alex !parameter is technically zero-dimensional
integer :: ios
integer :: i, j, inum, maxsize, ivals(NSpeciesTotal)
integer :: block(2)

block(1) = ib
block(2) = jb

! a temp array large enough to hold any of the
! Lon,Lat or Alt array from a block plus ghost cells
allocate(temp1d(1-nGhost:max(nLonsPerBlock,nLatsPerBlock,nAltsPerBlock)+nGhost))

! temp array large enough to hold 1 species, temperature, etc
allocate(temp3d(1-nGhost:nLonsPerBlock+nGhost, 1-nGhost:nLatsPerBlock+nGhost, &
   1-nGhost:nAltsPerBlock+nGhost))
allocate(data3d(1-nGhost:nLonsPerBlock+nGhost, 1-nGhost:nLatsPerBlock+nGhost, &
   1-nGhost:nAltsPerBlock+nGhost))

! temp array large enough to hold velocity vect, etc
maxsize = max(3, nSpecies)
allocate(temp4d(1-nGhost:nLonsPerBlock+nGhost, 1-nGhost:nLatsPerBlock+nGhost, &
   1-nGhost:nAltsPerBlock+nGhost, maxsize))

! copy over lat, lon, alt arrays verbatim
read(iunit,iostat=ios) temp1d(1-nGhost:nLonsPerBlock+nGhost)
if (ios /= 0) then
   write(string1,*)'unable to read lons from ',trim(infile)
   call error_handler(E_ERR,'write_data',string1,source,revision,revdate)
endif
print *,'now writing: lons'
write(ounit,iostat=ios) temp1d(1-nGhost:nLonsPerBlock+nGhost)
if (ios /= 0) then
   write(string1,*)'unable to write lons to ',trim(outfile)
   call error_handler(E_ERR,'write_data',string1,source,revision,revdate)
endif

read(iunit,iostat=ios) temp1d(1-nGhost:nLatsPerBlock+nGhost)
if (ios /= 0) then
   write(string1,*)'unable to read lats from ',trim(infile)
   call error_handler(E_ERR,'write_data',string1,source,revision,revdate)
endif
print *,'now writing: lat'
write(ounit,iostat=ios) temp1d(1-nGhost:nLatsPerBlock+nGhost)
if (ios /= 0) then
   write(string1,*)'unable to write lats to ',trim(outfile)
   call error_handler(E_ERR,'write_data',string1,source,revision,revdate)
endif

read(iunit,iostat=ios) temp1d(1-nGhost:nAltsPerBlock+nGhost)
if (ios /= 0) then
   write(string1,*)'unable to read alts from ',trim(infile)
   call error_handler(E_ERR,'write_data',string1,source,revision,revdate)
endif
print *,'now writing: alt'
write(ounit,iostat=ios) temp1d(1-nGhost:nAltsPerBlock+nGhost)
if (ios /= 0) then
   write(string1,*)'unable to write alts to ',trim(infile)
   call error_handler(E_ERR,'write_data',string1,source,revision,revdate)
endif

call get_index_from_gitm_varname('NDensityS', inum, ivals)
if (inum > 0) then
   ! if i equals ival, use the data from the state vect
   ! otherwise read/write what's in the input file
   j = 1
   do i = 1, nSpeciesTotal
      read(iunit)  temp3d

      if (j <= inum) then
         if (i == gitmvar(ivals(j))%gitm_index) then

            ! FIXME: if the program restart is really resetting the ghost zones
            ! correctly, then we shouldn't need to initialize the array with the
            ! temp3d data (which has pre-assimilation values in it).  but alexey
            ! says this causes problems, which is suspicious and should be looked
            ! at more.  this line might make it run, but the ghost zones were not
            ! updated by the assimilation.
            !alex: the horizontal ghost cells at middle altitudes (just not in the
            !altitude top and bottom ghost cells) should be fine as they get overwritten
            !in GITM. The top and bottom altitude ghost cells are the culprits.
            !Ideally, they should be extrapolated to once DART provides the posterior estimates,
            !but this is not implemented yet. This lack should not affect the
            !assimilation too much if the observations come from middle altitudes
            ! (as is the case with CHAMP and GRACE).
            data3d = temp3d

            call pack_data(ncid, ivals(j), block, data3d)

            ! FIXME: also needs fixing.  if we have made some value negative
            ! where the model doesn't support it, and it can't be 0 either, then
            ! this should be the smallest positive value that the model will accept.
            ! the original data divided by 2 is going to change the distribution of
            ! values and is certainly not right.  leave it here for now to get the
            ! assimilation running, but this needs looking at and changing soon.
            !alex: fixed on 5/20/13. How? Well, the limits of the variables are as
            !follows (taken from a GITM initialization on 12/1/2002 via gitm/matlab/rst2mat.m):
! MINIMA AND MAXIMA
! _
! LonT 0.17453 6.1087
! LatT -1.4835 1.4835
! AltT 100000 630038.9261
! TempT                  163.0163 1223.5239
! ITempT                 163.0154 1967.9977
! eTempT                 184.665 2710.9351
! _
! NDST_1,iO_3P_ 607871671694.2802 624953710511309568
! NDST_2,iO2_        1554285.5124 2977090934472271872
! NDST_3,iN2_      261275675.713  12920995180058857472
! NDST_4,iN_4S_   2725259865.9408 51174404943040.94
! NDST_5,iNO_             91.5983 137340019620842.8
! NDST_6,iN_2D_       490627.9878 656673564758.9565
! NDST_7,iN_2P_       135692.3204 2582963359.5952
! NDST_8,iH_    297189394877.7289 160285753329765.5
! NDST_9,iHe_    39396601335.7323 31530483811658.93
! NDST_10,iCO2_           51.3449 5237123737981628
! NDST_11,iO_1D_          32.2011 26604279344.3065
! _
! IDST_1,iO_4SP_         100 2345587357569.55
! IDST_2,iO2P_             4.0622 121043204145.427
! IDST_3,iN2P_             2.3259e-05 6408254083.7036
! IDST_4,iNP_              1.6073e-05 725487968.9667
! IDST_5,iNOP_            15.9515 182204005544.7968
! IDST_6,iO_2DP_           2.6996e-11 798313237.9133
! IDST_7,iO_2PP_           9.5018e-11 365561613.5574
! IDST_8,iHP_              1 250583438981.8537
! IDST_9,iHeP_             1 13445167727.3174
! IDST_10,ie_      543075134.2391 2346712140512.865
! _
! VT -253.1168 236.8601
! IVT -809.0201 1382.4808
! VVT -82.1731 633.1406
! RhoT 1.6352e-14 7.7801e-07
            !
            !So it makes sense to saturate NDST, IDST and Rho by 1.0e-16 from below,
            ! Temp, ITemp and eTemp by 100.0 from below, and F107 by 60.0 from below

!            where (data3d < 0.0_r8) data3d = temp3d/2 !alex, old - bad because might change distr
            where (data3d < 0.0_r8) data3d = 1.0e-16_r8 !alex, new

            print *,'now writing: ',trim(gitmvar(ivals(j))%varname)
            write(ounit) data3d
            j = j + 1
         else
            write(ounit) temp3d
         endif
      else
         ! this one not in state vector, copy over from input
         write(ounit) temp3d
      endif
   enddo
else
   ! nothing at all from this variable in the state vector.
   ! copy all data over from the input file to output file
   do i = 1, nSpeciesTotal
      read(iunit)  temp3d
      write(ounit) temp3d
   enddo
endif

call get_index_from_gitm_varname('IDensityS', inum, ivals)
if (inum > 0) then
   ! one or more items in the state vector need to replace the
   ! data in the output file.  loop over the index list in order.
   j = 1
   do i = 1, nIons
      read(iunit)  temp3d
      if (j <= inum) then
         if (i == gitmvar(ivals(j))%gitm_index) then
            ! read from input but write from state vector
            data3d = temp3d
            call pack_data(ncid, ivals(j), block, data3d)
            where (data3d < 0.0_r8) data3d = 1.0e-16_r8 !alex
            print *,'now writing: ',trim(gitmvar(ivals(j))%varname)
            write(ounit) data3d
            j = j + 1
         else
            write(ounit) temp3d
         endif
      else
         ! this one not in state vector, copy over from input
         write(ounit) temp3d
      endif
   enddo
else
   ! nothing at all from this variable in the state vector.
   ! copy all data over from the input file to output file
   do i = 1, nIons
      read(iunit)  temp3d
      write(ounit) temp3d
   enddo
endif

read(iunit)  temp3d
data3d = temp3d
call get_index_from_gitm_varname('Temperature', inum, ivals)
if (inum > 0) then
   call pack_data(ncid, ivals(1), block, data3d)
   where (data3d < 0.0_r8) data3d = 100.0_r8 !alex
   print *,'now writing: Temperature'
   write(ounit) data3d
else
   write(ounit) temp3d
endif


read(iunit) temp3d
data3d = temp3d
call get_index_from_gitm_varname('ITemperature', inum, ivals)
if (inum > 0) then
   call pack_data(ncid, ivals(1), block, data3d)
   where (data3d < 0.0_r8) data3d = 100.0_r8 !alex
   print *,'now writing: ITemperature'
   write(ounit) data3d
else
   write(ounit) temp3d
endif

read(iunit) temp3d
data3d = temp3d
call get_index_from_gitm_varname('eTemperature', inum, ivals)
if (inum > 0) then
   call pack_data(ncid, ivals(1), block, data3d)
   print *,'now writing: eTemperature'
   where (data3d < 0.0_r8) data3d = 100.0_r8 !alex
   write(ounit) data3d
else
   write(ounit) temp3d
endif

!print *, 'reading in temp4d for vel'
read(iunit) temp4d(:,:,:,1:3)
call get_index_from_gitm_varname('Velocity', inum, ivals)
if (inum > 0) then
   ! one or more items in the state vector need to replace the
   ! data in the output file.  loop over the index list in order.
   j = 1
   do i = 1, 3
      if (j <= inum) then
         if (i == gitmvar(ivals(j))%gitm_index) then
            print *,'now writing:',trim(gitmvar(ivals(j))%varname)
            ! read from input but write from state vector
            data3d = temp4d(:,:,:,i)
            call pack_data(ncid, ivals(j), block, data3d)
            temp4d(:,:,:,i) = data3d
            j = j + 1
         endif
      endif
   enddo
endif
write(ounit) temp4d(:,:,:,1:3)

!print *, 'reading in temp4d for ivel'
read(iunit) temp4d(:,:,:,1:3)
call get_index_from_gitm_varname('IVelocity', inum, ivals)
if (inum > 0) then
   ! one or more items in the state vector need to replace the
   ! data in the output file.  loop over the index list in order.
   j = 1
   do i = 1, 3
      if (j <= inum) then
         if (i == gitmvar(ivals(j))%gitm_index) then
            print *,'now writing:',trim(gitmvar(ivals(j))%varname)
            ! read from input but write from state vector
            data3d = temp4d(:,:,:,i)
            call pack_data(ncid, ivals(j), block, data3d)
            temp4d(:,:,:,i) = data3d
            j = j + 1
         endif
      endif
   enddo
endif
write(ounit) temp4d(:,:,:,1:3)

!print *, 'reading in temp4d for vvel'
read(iunit) temp4d(:,:,:,1:nSpecies)
call get_index_from_gitm_varname('VerticalVelocity', inum, ivals)
if (inum > 0) then
   ! one or more items in the state vector need to replace the
   ! data in the output file.  loop over the index list in order.
   j = 1
   do i = 1, nSpecies
      if (j <= inum) then
         if (i == gitmvar(ivals(j))%gitm_index) then
            print *,'now writing:',trim(gitmvar(ivals(j))%varname)
            ! read from input but write from state vector
            data3d = temp4d(:,:,:,i)
            call pack_data(ncid, ivals(j), block, data3d)
            temp4d(:,:,:,i) = data3d
            j = j + 1
         endif
      endif
   enddo
endif
write(ounit) temp4d(:,:,:,1:nSpecies)


!alex begin: added f107 and Rho to the restart files:
read(iunit) temp0d
data0d = temp0d
call get_index_from_gitm_varname('f107', inum, ivals)
if (inum > 0) then
   call pack_data0d(ncid, ivals(1), data0d)
   if (data0d < 0.0_r8) data0d = 60.0_r8 !alex
   write(ounit) data0d
else
   write(ounit) temp0d
endif

read(iunit)  temp3d
data3d = temp3d
call get_index_from_gitm_varname('Rho', inum, ivals)
if (inum > 0) then
   call pack_data(ncid, ivals(1), block, data3d)
   where (data3d < 0.0_r8) data3d = 1.0e-16_r8 !alex
   write(ounit) data3d
else
   write(ounit) temp3d
endif
!alex end

deallocate(temp1d, temp3d, temp4d, data3d)

end subroutine write_data


!==================================================================

!> return 0 (ok) if we know how to interpolate this quantity.
!> if it is a field in the state, return the variable id from
!> the state structure.  if not in the state, varid will return -1

subroutine ok_to_interpolate(obs_qty, varid, my_status)

integer, intent(in)  :: obs_qty
integer, intent(out) :: varid
integer, intent(out) :: my_status

! See if the state contains the obs quantity
varid = get_varid_from_kind(domain_id, obs_qty)

! in the state vector
if (varid > 0) then
   my_status = 0
   return
endif

! add any quantities that can be interpolated to this list if they
! are not in the state vector.
select case (obs_qty)
   case (QTY_GEOMETRIC_HEIGHT,  &
         QTY_VERTLEVEL)
      my_status = 0
   case default
      my_status = UNKNOWN_OBS_QTY_ERROR_CODE
end select


end subroutine ok_to_interpolate


!------------------------------------------------------------------
! Determine where any data from a given gitm_varname lies in the
! DART state vector.

subroutine get_index_from_gitm_varname(gitm_varname, inum, ivals)

character(len=*), intent(in) :: gitm_varname
integer, intent(out) :: inum, ivals(:)

integer :: gindex(nfields)
integer :: i, limit

inum = 0
limit = size(ivals)

! GITM handles variables in a way that might seem strange at first.
! It uses the same name but multiple indices. For example, the U, V, 
! and W components of wind are index = 1, 2, 3 for the variable velocity.
! This is why the code below looks the way it does. 
FieldLoop : do i=1,nfields
   if (gitmvar(i)%gitm_varname /= gitm_varname) cycle FieldLoop
   inum = inum + 1
   if (inum > limit) then
      write(string1,*) 'found too many matches, ivals needs to be larger than ', limit
      call error_handler(E_ERR,'get_index_from_gitm_varname',string1,source,revision,revdate)
   endif
   ! i is index into gitmvar array - the order of the fields in the sv
   ! gitm_index is index into the specific variable in the gitm restarts
   ivals(inum) = i
   gindex(inum) = gitmvar(i)%gitm_index
enddo FieldLoop

!if (inum > 0) then
!   print *, 'before sort, inum: ', inum
!   print *, 'before sort, gindex: ', gindex(1:inum)
!   print *, 'before sort, ivals: ', ivals(1:inum)
!endif

! return the vals sorted by gitm_index order if more than 1
if (inum > 1) call sortindexlist(gindex, ivals, inum)

!if (inum > 0) then
!   print *, 'after  sort, inum: ', inum
!   print *, 'after  sort, gindex: ', gindex(1:inum)
!   print *, 'after  sort, ivals: ', ivals(1:inum)
!endif

end subroutine get_index_from_gitm_varname


!------------------------------------------------------------------
! the static_init_model ensures that the gitm namelists are read.

function set_model_time_step()
type(time_type) :: set_model_time_step

set_model_time_step = set_time(assimilation_period_seconds, &
                               assimilation_period_days) ! (seconds, days)

end function set_model_time_step



!------------------------------------------------------------------
subroutine verify_state_variables( variable_array, ngood, table )

character(len=*), dimension(:),   intent(in)  :: variable_array
integer,                          intent(out) :: ngood
character(len=*), dimension(:,:), intent(out) :: table

integer :: nrows, i
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr

character(len=*), parameter :: routine = 'verify_state_variables'

nrows = size(table,1)

ngood = 0
MyLoop : do i = 1, nrows

   varname   = variable_array(num_state_table_columns*i-4)
   dartstr   = variable_array(num_state_table_columns*i-3)
   !minvalstr = variable_array(num_state_table_columns*i-2)
   !maxvalstr = variable_array(num_state_table_columns*i-1)
   !updatestr = variable_array(num_state_table_columns*i  )

   table(i,1) = trim(varname)
   table(i,2) = trim(dartstr)

   if ( table(i,1) == ' ' .and. table(i,2) == ' ' ) exit MyLoop ! Found end of list.

   if ( table(i,1) == ' ' .or. table(i,2) == ' ' ) then
      print *,'i is:',i,'nrows:',nrows
      string1 = 'model_nml:gitm_state_variables not fully specified'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   endif

   ! Because of the way the blocks to netCDF is coded, we are normally sure the 
   ! variable exists in GITM restart variable list

   ! Make sure DART kind is valid
   if( get_index_for_quantity(dartstr) < 0 ) then
      write(string1,'(3A)') 'there is no obs_kind "', trim(dartstr), '" in obs_kind_mod.f90'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   end if

   ! Record the contents of the DART state vector

   if ( debug > 0 ) then
      write(string1,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
      call error_handler(E_MSG,routine,string1,source,revision,revdate)
   endif

   ngood = ngood + 1
enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)
endif

end subroutine verify_state_variables


!------------------------------------------------------------------
subroutine verify_block_variables( variable_array, ngood)

character(len=*), dimension(:),   intent(in)  :: variable_array
integer,                          intent(out) :: ngood

integer :: nrows, i
character(len=NF90_MAX_NAME) :: varname

character(len=*), parameter :: routine = 'verify_state_variables'

nrows = size(variable_array,1)

ngood = 0
MyLoop : do i = 1, nrows

   varname   = variable_array(i)

   if ( varname  == ' ') exit MyLoop ! Found end of list.

   ngood = ngood + 1
enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)
endif

end subroutine verify_block_variables


!------------------------------------------------------------------
function read_in_real(iunit,varname,filename)

integer,          intent(in) :: iunit
character(len=*), intent(in) :: varname,filename
real(r8)                     :: read_in_real

character(len=100) :: cLine
integer :: i, ios

! Read a line and remove anything after a space or TAB
read(iunit,'(a)',iostat=ios) cLine
if (ios /= 0) then
   write(string1,*) 'cannot find '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'get_grid_dims',string1,source,revision,revdate)
endif

i=index(cLine,' ');     if( i > 0 ) cLine(i:len(cLine))=' '
i=index(cLine,char(9)); if( i > 0 ) cLine(i:len(cLine))=' '

! Now that we have a line with nothing else ... parse it
read(cLine,*,iostat=ios)read_in_real

if(ios /= 0) then
   write(string1,*)'unable to read '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'read_in_real',string1,source,revision,revdate)
endif

end function read_in_real


!------------------------------------------------------------------
function read_in_int(iunit,varname,filename)

integer,          intent(in) :: iunit
character(len=*), intent(in) :: varname,filename
integer                      :: read_in_int

character(len=100) :: cLine
integer :: i, ios

! Read a line and remove anything after a space or TAB
read(iunit,'(a)',iostat=ios) cLine
if (ios /= 0) then
   write(string1,*) 'cannot find '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'get_grid_dims',string1,source,revision,revdate)
endif

i=index(cLine,' ');     if( i > 0 ) cLine(i:len(cLine))=' '
i=index(cLine,char(9)); if( i > 0 ) cLine(i:len(cLine))=' '

read(cLine,*,iostat=ios)read_in_int

if(ios /= 0) then
   write(string1,*)'unable to read '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'read_in_int',string1,source,revision,revdate,&
             text2=cLine)
endif

end function read_in_int


!------------------------------------------------------------------
!> sort list x into order based on values in list.
!> should only be called on short ( < hundreds) of values or will be slow
!> @todo FIXME this should be using the sort module routine instead.

subroutine sortindexlist(list, x, inum)

integer, intent(inout) :: list(:)
integer, intent(inout) :: x(:)
integer, intent(in)    :: inum

integer :: tmp
integer :: j, k

!  DO A N^2 SORT - only use for short lists
do j = 1, inum - 1
   do k = j + 1, inum
      ! if list() is in wrong order, exchange both list items and
      ! items in x array.
      if(list(j) .gt. list(k)) then
         tmp = list(k)
         list(k) = list(j)
         list(j) = tmp
         tmp = x(k)
         x(k) = x(j)
         x(j) = tmp
      end if
   end do
end do
end subroutine sortindexlist


!===================================================================
end module model_mod
!===================================================================

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
