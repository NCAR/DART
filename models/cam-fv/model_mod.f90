! DART software - Copyright UCAR. This open source software is provided
! by ucar, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/dares/dart/dart_download
!
! $Id$
!----------------------------------------------------------------
!>
!> this is the interface between the cam-fv atmosphere model and dart.
!> the required public interfaces and arguments cannot be changed.
!>
!----------------------------------------------------------------

module model_mod

use             types_mod,  only : MISSING_R8, MISSING_I, i8, r8, vtablenamelength, &
                                   gravity, DEG2RAD
use      time_manager_mod,  only : set_time, time_type, set_date, &
                                   set_calendar_type, get_date
use          location_mod,  only : location_type, set_vertical, set_location, &
                                   get_location, write_location, is_vertical, &
                                   VERTISUNDEF, VERTISSURFACE, VERTISLEVEL, &
                                   VERTISPRESSURE, VERTISHEIGHT, &
                                   VERTISSCALEHEIGHT, query_location, &
                                   set_vertical_localization_coord, get_dist, &
                                   loc_get_close_obs => get_close_obs, &
                                   loc_get_close_state => get_close_state, &
                                   vertical_localization_on, get_close_type, get_maxdist
use         utilities_mod,  only : find_namelist_in_file, check_namelist_read, &
                                   string_to_logical, string_to_real,& 
                                   nmlfileunit, do_nml_file, do_nml_term, &
                                   register_module, error_handler, &
                                   file_exist, to_upper, E_ERR, E_MSG, array_dump, &
                                   find_enclosing_indices
use          obs_kind_mod,  only : QTY_SURFACE_ELEVATION, QTY_PRESSURE, &
                                   QTY_GEOMETRIC_HEIGHT, QTY_VERTLEVEL, &
                                   QTY_SURFACE_PRESSURE, &
                                   QTY_TEMPERATURE, QTY_SPECIFIC_HUMIDITY, &
                                   QTY_MOLEC_OXYGEN_MIXING_RATIO, &
                                   QTY_ION_O_MIXING_RATIO, QTY_ATOMIC_H_MIXING_RATIO, &
                                   QTY_ATOMIC_OXYGEN_MIXING_RATIO, QTY_NITROGEN, &
                                   get_index_for_quantity, get_num_quantities, &
                                   get_name_for_quantity, get_quantity_for_type_of_obs
use     mpi_utilities_mod,  only : my_task_id
use        random_seq_mod,  only : random_seq_type, init_random_seq, random_gaussian
use  ensemble_manager_mod,  only : ensemble_type, get_my_num_vars, get_my_vars
use distributed_state_mod,  only : get_state
use   state_structure_mod,  only : add_domain, get_dart_vector_index, get_domain_size, &
                                   get_dim_name, get_kind_index, get_num_dims, &
                                   get_num_variables, get_varid_from_kind, &
                                   get_model_variable_indices, state_structure_info
use  netcdf_utilities_mod,  only : nc_get_variable, nc_get_variable_size, &
                                   nc_add_attribute_to_variable, &
                                   nc_define_integer_variable, &
                                   nc_define_real_variable, &
                                   nc_define_real_scalar, &
                                   nc_add_global_creation_time, &
                                   nc_add_global_attribute, &
                                   nc_define_dimension, nc_put_variable, &
                                   nc_synchronize_file, nc_end_define_mode, &
                                   nc_begin_define_mode, nc_open_file_readonly, &
                                   nc_close_file, nc_variable_exists
use        chem_tables_mod, only : init_chem_tables, finalize_chem_tables, &
                                   get_molar_mass, get_volume_mixing_ratio
use        quad_utils_mod,  only : quad_interp_handle, init_quad_interp, &
                                   set_quad_coords, finalize_quad_interp, &
                                   quad_lon_lat_locate, quad_lon_lat_evaluate, &
                                   GRID_QUAD_IRREG_SPACED_REGULAR,  &
                                   QUAD_LOCATED_CELL_CENTERS
use     default_model_mod,  only : adv_1step, nc_write_model_vars, &
                                   init_time => fail_init_time,    &
                                   init_conditions => fail_init_conditions

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the dart code.

! routines in this list have code in this module
public :: static_init_model,                   &
          get_model_size,                      &
          get_state_meta_data,                 &
          model_interpolate,                   & 
          shortest_time_between_assimilations, &
          nc_write_model_atts,                 &
          write_model_time,                    & 
          read_model_time,                     &
          end_model,                           &
          pert_model_copies,                   & 
          convert_vertical_obs,                & 
          convert_vertical_state,              & 
          get_close_obs,                       &
          get_close_state

! code for these routines are in other modules
public :: nc_write_model_vars,           &
          adv_1step,                     &
          init_time,                     &
          init_conditions

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! maximum number of fields you can list to be perturbed
! to generate an ensemble if starting from a single state.
integer, parameter :: MAX_PERT = 100

! model_nml namelist variables and default values
character(len=256) :: cam_template_filename           = 'caminput.nc'
character(len=256) :: cam_phis_filename               = 'cam_phis.nc'
character(len=32)  :: vertical_localization_coord     = 'PRESSURE'
logical            :: use_log_vertical_scale          = .false.
integer            :: assimilation_period_days        = 0
integer            :: assimilation_period_seconds     = 21600
! proposed changes:
integer            :: no_obs_assim_above_level       = -1      ! model levels
integer            :: model_damping_ends_at_level    = -1      ! model levels
! end proposed changes
integer            :: debug_level                     = 0
logical            :: suppress_grid_info_in_output    = .false.
logical            :: custom_routine_to_generate_ensemble = .true.
character(len=32)  :: fields_to_perturb(MAX_PERT)     = ""
real(r8)           :: perturbation_amplitude(MAX_PERT)= 0.0_r8
logical            :: using_chemistry                 = .false.
logical            :: use_variable_mean_mass          = .false.

! in converting to scale height for the vertical: 
!  set this to .true. to compute the log of the pressure.  
!  set this to .false. to additionally normalize by the surface
!    pressure at this location.  this option is backwards-compatible 
!    with previous versions of this code.
logical :: no_normalization_of_scale_heights = .true.

! state_variables defines the contents of the state vector.
! each line of this input should have the form:
!
!    netcdf_variable_name, dart_quantity, clamp_min, clamp_max, update_variable
!
! all items must be strings (even if numerical values).
! for no clamping, use the string 'NA'
! to have the assimilation change the variable use 'UPDATE', else 'NO_UPDATE'

integer, parameter :: MAX_STATE_VARIABLES = 100
integer, parameter :: num_state_table_columns = 5
character(len=vtablenamelength) :: state_variables(MAX_STATE_VARIABLES * &
                                                   num_state_table_columns ) = ' '

namelist /model_nml/  &
   cam_template_filename,               &
   cam_phis_filename,                   &
   vertical_localization_coord,         &
   state_variables,                     &
   assimilation_period_days,            &
   assimilation_period_seconds,         &
   use_log_vertical_scale,              &
   no_obs_assim_above_level,            & 
   model_damping_ends_at_level,         &
   suppress_grid_info_in_output,        &
   custom_routine_to_generate_ensemble, &
   fields_to_perturb,                   &
   perturbation_amplitude,              &
   no_normalization_of_scale_heights,   &
   use_variable_mean_mass,              &
   using_chemistry,                     &
   debug_level

! global variables
character(len=512) :: string1, string2, string3
logical, save      :: module_initialized = .false.

! this id allows us access to all of the state structure
! info and is required for getting state variables.
integer :: domain_id

!> Metadata from the template netCDF file that describes 
!> where the variable data is located and what size it is.

type cam_1d_array
   integer  :: nsize
   real(r8), allocatable :: vals(:)
end type

type cam_grid
   type(cam_1d_array) :: lon
   type(cam_1d_array) :: lat
   type(cam_1d_array) :: slon
   type(cam_1d_array) :: slat
   type(cam_1d_array) :: lev
   type(cam_1d_array) :: ilev
   type(cam_1d_array) :: gw
   type(cam_1d_array) :: hyai
   type(cam_1d_array) :: hybi
   type(cam_1d_array) :: hyam
   type(cam_1d_array) :: hybm
   type(cam_1d_array) :: P0
end type

type(cam_grid) :: grid_data


integer, parameter :: STAGGER_NONE = -1
integer, parameter :: STAGGER_U    =  1
integer, parameter :: STAGGER_V    =  2
integer, parameter :: STAGGER_W    =  3 
integer, parameter :: STAGGER_UV   =  4

type cam_stagger
   integer, allocatable :: qty_stagger(:)
end type

type(cam_stagger) :: grid_stagger

! Surface potential; used for calculation of geometric heights.
real(r8), allocatable :: phis(:, :)

! default to localizing in pressure.  override with namelist
integer :: vertical_localization_type = VERTISPRESSURE

! flag used to know if the vertical unit system has numbers
! that get larger as you move away from the earth's surface
! (e.g. height) or smaller (e.g. pressure)
logical :: higher_is_smaller

! commonly used numbers that we'll set in static_init_model
real(r8) :: ref_model_top_pressure
real(r8) :: ref_surface_pressure
integer  :: ref_nlevels

!>@todo FIXME ask kevin if this threshold value is small enough
! to distinguish cam from waccm configurations?

! an arbitrary value to test the model top against to see
! if we're running cam vs waccm or waccm-x.  it changes the
! standard atmosphere table we use to convert pressure to height, 
! and changes the formatting of numbers in dart_log output.
real(r8), parameter :: high_top_threshold = 0.3_r8  ! pascals

! things related to damping at the model top
logical  :: are_damping = .false.
real(r8) :: ramp_end         ! fixed top of ramp; the start (bottom) varies
logical  :: discarding_high_obs = .false.
real(r8) :: no_assim_above_height    = -1.0_r8 
real(r8) :: no_assim_above_level     = -1.0_r8 
real(r8) :: no_assim_above_pressure  = -1.0_r8 

!> build a pressure/height conversion column based on a
!> standard atmosphere.  this can only be used when we
!> don't have a real ensemble to use, or we don't care
!> about absolute accuracy.

interface single_pressure_value
 module procedure single_pressure_value_int
 module procedure single_pressure_value_real
end interface

! Precompute pressure <-> height map once based on either a low-top or
! high-top table depending on what the model top is.
! Used only to discard obs on heights above the user-defined top threshold.
integer, parameter :: HIGH_TOP_TABLE = 1
integer, parameter ::  LOW_TOP_TABLE = 2
integer :: std_atm_table_len
real(r8), allocatable :: std_atm_hgt_col(:)
real(r8), allocatable :: std_atm_pres_col(:)

! Horizontal interpolation code.  Need a handle for nonstaggered, U and V.
type(quad_interp_handle) :: interp_nonstaggered, &
                            interp_u_staggered, &
                            interp_v_staggered


contains


!-----------------------------------------------------------------------
! All the required interfaces are first.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Called to do one time initialization of the model.
!> In this case, it reads in the grid information, the namelist
!> containing the variables of interest, where to get them, their size,
!> their associated DART Quantity, etc.
!>
!> In addition to harvesting the model metadata (grid,
!> desired model advance step, etc.), it also fills a structure
!> containing information about what variables are where in the DART
!> framework.

subroutine static_init_model()

integer :: iunit, io
integer :: nfields

character(len=*), parameter :: routine = 'static_init_model'

if ( module_initialized ) return

! Record version info
call register_module(source, revision, revdate)

module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

call set_calendar_type('GREGORIAN')

call read_grid_info(cam_template_filename, grid_data)

! initialize global values that are used frequently
call init_globals()

! read the namelist &model_nml :: state_variables
! to set up what will be read into the cam state vector
call set_cam_variable_info(state_variables, nfields)

! convert from string in namelist to integer (e.g. VERTISxxx)
! and tell the dart code which vertical type we want to localize in.
call set_vert_localization(vertical_localization_coord)

! if you have chemistry variables in the model state, set
! this namelist variable so we can initialize the proper tables
if (using_chemistry) call init_chem_tables()

! set top limit where obs impacts are diminished to 0.
! only allowed if doing vertical localization.  error if
! computing horizontal distances only (odd case, intentionally
! choosing not to support this.)
if (model_damping_ends_at_level > 0) then
   if (vertical_localization_on()) then
      call init_damping_ramp_info()
      are_damping = .true.
   else
      string1='cannot support model top damping unless also using vertical localization'
      string2='set "model_damping_ends_at_level = -1" in &model_nml, OR' 
      string3='set "horiz_dist_only = .false." in &location_nml'
      call error_handler(E_ERR, routine, string1, source, revision, revdate, &
                         text2=string2, text3=string3)
   endif
endif

! set top limit where obs are discarded.  -1 to disable.
if (no_obs_assim_above_level > 0) then
   call init_discard_high_obs()
   discarding_high_obs = .true.
endif

! set a flag based on the vertical localization coordinate selected
call init_sign_of_vert_units()

end subroutine static_init_model


!-----------------------------------------------------------------------
!> Returns the size of the DART state vector (i.e. model) as an integer.
!>

function get_model_size()

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = get_domain_size(domain_id)

end function get_model_size



!-----------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location. A second intent(out) optional argument quantity
!> can be returned if the model has more than one type of field (for
!> instance temperature and zonal wind component). This interface is
!> required for all filter applications as it is required for computing
!> the distance between observations and state variables.
!>
!> @param index_in the index into the DART state vector
!> @param location the location at that index
!> @param var_type the DART Quantity at that index
!>

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

! Local variables

integer  :: iloc, vloc, jloc
integer  :: myvarid, myqty, nd

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, iloc, jloc, vloc, var_id=myvarid, kind_index=myqty)

nd = get_num_dims(domain_id, myvarid)

location = get_location_from_index(iloc, jloc, vloc, myqty, nd)

! return state quantity for this index if requested
if (present(var_type)) var_type = myqty

end subroutine get_state_meta_data

!-----------------------------------------------------------------------
!> given the (i,j,k) indices into a field in the state vector,
!> and the quantity, and the dimensionality of the field (2d, 3d),
!> compute the location of that item.  

function get_location_from_index(i, j, k, q, nd)
integer, intent(in) :: i
integer, intent(in) :: j
integer, intent(in) :: k
integer, intent(in) :: q
integer, intent(in) :: nd
type(location_type) :: get_location_from_index

real(r8) :: slon_val
real(r8) :: use_vert_val
integer  :: use_vert_type

! full 3d fields are returned with lon/lat/level.
! 2d fields are either surface fields, or if they
! are column integrated values then they are 'undefined'
! in the vertical.

if (nd == 3) then
   use_vert_type = VERTISLEVEL
   use_vert_val  = real(k,r8)
else
   if (q == QTY_SURFACE_ELEVATION .or. q == QTY_SURFACE_PRESSURE) then
      use_vert_type = VERTISSURFACE
      use_vert_val  = MISSING_R8  
      ! setting the vertical value to missing matches what the previous
      ! version of this code did.  other models choose to set the vertical
      ! value to the actual surface elevation at this location:
      !   use_vert_val  = phis(lon_index, lat_index) / gravity
   else
      ! assume other 2d fields are integrated quantities with no vertical
      ! location. if there are other real surface fields in the state
      ! add their quantitys to the if() test above.
      use_vert_type = VERTISUNDEF
      use_vert_val  = MISSING_R8
   endif
endif

! the horizontal location depends on whether this quantity is on the
! mass point grid or staggered in either lat or lon.  

select case (grid_stagger%qty_stagger(q))
  case (STAGGER_U)
   get_location_from_index = set_location(grid_data%lon%vals(i), &
                                          grid_data%slat%vals(j), &
                                          use_vert_val, use_vert_type)

  case (STAGGER_V)
   ! the first staggered longitude is negative.  dart requires lons
   ! be between 0 and 360.
   slon_val = grid_data%slon%vals(i)
   if (slon_val < 0) slon_val = slon_val + 360.0_r8
   get_location_from_index = set_location(slon_val, &
                                          grid_data%lat%vals(j), &
                                          use_vert_val, use_vert_type)
   
  !>@todo not sure what to do yet. ? +-1/2 ?
  case (STAGGER_W)
   get_location_from_index = set_location(grid_data%lon%vals(i), &
                                          grid_data%lat%vals(j), &
                                          use_vert_val - 0.5_r8, use_vert_type)
  ! no stagger - cell centers
  case default
   get_location_from_index = set_location(grid_data%lon%vals(i), &
                                          grid_data%lat%vals(j), &
                                          use_vert_val, use_vert_type)

end select

end function get_location_from_index

!-----------------------------------------------------------------------
!> this routine should be called to compute a value that comes from an
!> unstaggered grid but needs to correspond to a staggered grid.
!> e.g. you need the surface pressure under a V wind point.

subroutine get_staggered_values_from_qty(ens_handle, ens_size, qty, lon_index, lat_index, &
                                         lev_index, stagger_qty, vals, my_status)
type(ensemble_type), intent(in) :: ens_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: qty
integer,             intent(in) :: lon_index
integer,             intent(in) :: lat_index
integer,             intent(in) :: lev_index
integer,             intent(in) :: stagger_qty
real(r8),            intent(out) :: vals(ens_size)
integer,             intent(out) :: my_status

integer :: next_lat, prev_lon, stagger
real(r8) :: vals1(ens_size), vals2(ens_size)

vals(:) = MISSING_R8
stagger = grid_stagger%qty_stagger(stagger_qty)

!> latitudes:  staggered value N is between N and (N + 1) on the unstaggered grid
!> longitudes: staggered value N is between N and (N - 1) on the unstaggered grid

select case (stagger)
  case (STAGGER_U)
   call quad_index_neighbors(lon_index, lat_index, prev_lon, next_lat)

   call get_values_from_single_level(ens_handle, ens_size, qty, lon_index, lat_index, lev_index, &
                                     vals1, my_status)
   if (my_status /= 0) return
   call get_values_from_single_level(ens_handle, ens_size, qty, lon_index, next_lat,  lev_index, &
                                     vals2, my_status)
   if (my_status /= 0) return

   vals = (vals1 + vals2) * 0.5_r8

  case (STAGGER_V)
   call quad_index_neighbors(lon_index, lat_index, prev_lon, next_lat)

   call get_values_from_single_level(ens_handle, ens_size, qty, lon_index, lat_index, lev_index, &
                                     vals1, my_status)
   if (my_status /= 0) return
   call get_values_from_single_level(ens_handle, ens_size, qty, prev_lon,  lat_index, lev_index, &
                                     vals2, my_status)
   if (my_status /= 0) return

   vals = (vals1 + vals2) * 0.5_r8

  ! no stagger - cell centers, or W stagger
  case default
   call get_values_from_single_level(ens_handle, ens_size, qty, lon_index, lat_index, lev_index, &
                                     vals, my_status)
   if (my_status /= 0) return

end select

! when you reach here, my_status has been to 0 by the last call
! to get_values_from_single_level().  if it was anything else
! it would have already returned.

end subroutine get_staggered_values_from_qty


!-----------------------------------------------------------------------
!> this routine converts the 3 index values and a quantity into a state vector
!> offset and gets the ensemble of state values for that offset.  this only
!> gets a single vertical location - if you need to get values which might 
!> have different vertical locations in different ensemble members
!> see get_values_from_varid() below.

subroutine get_values_from_single_level(ens_handle, ens_size, qty, lon_index, lat_index, lev_index, &
                                        vals, my_status)
type(ensemble_type), intent(in) :: ens_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: qty
integer,             intent(in) :: lon_index
integer,             intent(in) :: lat_index
integer,             intent(in) :: lev_index
real(r8),            intent(out) :: vals(ens_size)
integer,             intent(out) :: my_status

character(len=*), parameter :: routine = 'get_values_from_single_level:'

integer :: varid
integer(i8) :: state_indx

varid = get_varid_from_kind(domain_id, qty)
if (varid < 0) then
   vals(:) = MISSING_R8
   my_status = 12
   return
endif

state_indx = get_dart_vector_index(lon_index, lat_index, lev_index, domain_id, varid)
if (state_indx < 1 .or. state_indx > get_domain_size(domain_id)) then
   write(string1, *) 'state_index out of range: ', state_indx, ' not between ', 1, get_domain_size(domain_id)
   call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2,text3='should not happen')
endif
vals(:) = get_state(state_indx, ens_handle)

my_status = 0

end subroutine get_values_from_single_level


!-----------------------------------------------------------------------
!> this routine takes care of getting the actual state values.  get_state()
!> communicates with other MPI tasks and can be expensive.
!>
!> all ensemble members have the same horizontal location, but different 
!> ensemble members could have different vertical locations and
!> so be between different vertical layers.  this code tries to do the fewest
!> calls to get_state by only calling it for levels that are actually needed
!> and setting all members with those same levels in a single pass.
!> 

subroutine get_values_from_varid(ens_handle, ens_size, lon_index, lat_index, lev_index, varid, &
                                 vals, my_status)
type(ensemble_type), intent(in)  :: ens_handle
integer,  intent(in)  :: ens_size
integer,  intent(in)  :: lon_index
integer,  intent(in)  :: lat_index
integer,  intent(in)  :: lev_index(ens_size)
integer,  intent(in)  :: varid
real(r8), intent(out) :: vals(ens_size)
integer,  intent(out) :: my_status(ens_size)

integer(i8) :: state_indx
integer  :: i, j
real(r8) :: temp_vals(ens_size) 
logical  :: member_done(ens_size)

character(len=*), parameter :: routine = 'get_values_from_varid:'

! as we get the values for each ensemble member, we set the 'done' flag
! and a good return code. 
my_status(:) = 12
member_done(:) = .false.

! start with lev_index(1).  get the vals into a temp var.  
! run through 2-N. any other member that has the same level 
! set the outgoing values.  keep a separate flag for which 
! member(s) have been done.  skip to the next undone member 
! and get the state for that level.  repeat until all levels done.

do i=1, ens_size

   if (member_done(i)) cycle

   state_indx = get_dart_vector_index(lon_index, lat_index, lev_index(i), domain_id, varid)

   if (state_indx < 0) then
      write(string1,*) 'Should not happen: could not find dart state index from '
      write(string2,*) 'lon, lat, and lev index :', lon_index, lat_index, lev_index
      call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2)
      return
   endif

   temp_vals(:) = get_state(state_indx, ens_handle)    ! all the ensemble members for level (i)

   ! start at i, because my ensemble member is clearly at this level.
   ! then continue on to see if any other members are also at this level.
   do j=i, ens_size
      if (member_done(j)) cycle

      if (lev_index(j) == lev_index(i)) then
         vals(j) = temp_vals(j)
         member_done(j) = .true.
         my_status(j) = 0
      endif
         
   enddo
enddo

end subroutine get_values_from_varid

!-----------------------------------------------------------------------
!> this is just for 3d fields

subroutine get_values_from_nonstate_fields(ens_handle, ens_size, lon_index, lat_index, &
                                           lev_index, obs_quantity, vals, my_status)
type(ensemble_type),  intent(in)  :: ens_handle
integer,              intent(in)  :: ens_size
integer,              intent(in)  :: lon_index
integer,              intent(in)  :: lat_index
integer,              intent(in)  :: lev_index(ens_size)
integer,              intent(in)  :: obs_quantity
real(r8),             intent(out) :: vals(ens_size)
integer,              intent(out) :: my_status(ens_size)

integer  :: imember
real(r8) :: vals_array(ref_nlevels,ens_size)

character(len=*), parameter :: routine = 'get_values_from_nonstate_fields:'

vals(:) = MISSING_R8
my_status(:) = 99

select case (obs_quantity) 
   case (QTY_PRESSURE)
      call cam_pressure_levels(ens_handle, ens_size, &
                               lon_index, lat_index, ref_nlevels, &
                               obs_quantity, vals_array, my_status)
      if (any(my_status /= 0)) return

      do imember=1,ens_size
         vals(imember) = vals_array(lev_index(imember), imember)
      enddo

   case (QTY_VERTLEVEL)
      vals(:)      = lev_index(:)
      my_status(:) = 0

   case default
      write(string1,*)'contact dart support. unexpected error for quantity ', obs_quantity
      call error_handler(E_MSG,routine,string1,source,revision,revdate)

end select

end subroutine get_values_from_nonstate_fields

!-----------------------------------------------------------------------
!>
!> Model interpolate will interpolate any DART state variable
!> to the given location.
!>
!> @param state_handle DART ensemble handle
!> @param ens_size DART ensemble size
!> @param location the location of interest
!> @param obs_qty the DART Quantity of interest
!> @param interp_vals the estimated value of the DART state at the location
!>          of interest (the interpolated value).
!> @param istatus interpolation status ... 0 == success, /=0 is a failure
!>
!> istatus = 2    asked to interpolate an unknown/unsupported quantity
!> istatus = 3    cannot locate horizontal quad
!> istatus = 4    cannot locate enclosing vertical levels
!> istatus = 5    cannot retrieve state vector values
!> istatus = 6    cannot get values at quad corners
!> istatus = 7    unused (error code available)
!> istatus = 8    cannot interpolate in the quad to get the values
!> istatus = 9    unused (error code available)
!> istatus = 10   cannot get vertical levels for an obs on pressure levels
!> istatus = 11   cannot get vertical levels for an obs on height levels
!> istatus = 12   cannot get values from obs quantity
!> istatus = 13   can not interpolate values of this quantity
!> istatus = 14   obs above user-defined assimilation top pressure
!> istatus = 15   can not get indices from given state vector index
!> istatus = 16   cannot do vertical interpolation for bottom layer
!> istatus = 17   cannot do vertical interpolation for top layer
!> istatus = 98   unknown error - shouldn't happen
!> istatus = 99   unknown error - shouldn't happen
!>

subroutine model_interpolate(state_handle, ens_size, location, obs_qty, interp_vals, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_qty
real(r8),           intent(out) :: interp_vals(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

character(len=*), parameter :: routine = 'model_interpolate:'

integer  :: varid, which_vert, status1
integer  :: four_lons(4), four_lats(4)
integer  :: status_array(ens_size)
real(r8) :: lon_fract, lat_fract
real(r8) :: lon_lat_vert(3)
real(r8) :: quad_vals(4, ens_size)
type(quad_interp_handle) :: interp_handle   ! should this be a pointer?? 
                                            ! is it replicating the internal arrays on assignment?

if ( .not. module_initialized ) call static_init_model


! Successful istatus is 0
interp_vals(:) = MISSING_R8
istatus(:)     = 99

! do we know how to interpolate this quantity?
call ok_to_interpolate(obs_qty, varid, status1)

if (status1 /= 0) then  
   if(debug_level > 12) then
      write(string1,*)'did not find observation quantity ', obs_qty, ' in the state vector'
      call error_handler(E_MSG,routine,string1,source,revision,revdate)
   endif
   istatus(:) = status1   ! this quantity not in the state vector
   return
endif

! get the grid handle for the right staggered grid
interp_handle = get_interp_handle(obs_qty)

! unpack the location type into lon, lat, vert, vert_type
lon_lat_vert = get_location(location)
which_vert   = nint(query_location(location)) 

! get the indices for the 4 corners of the quad in the horizontal, plus
! the fraction across the quad for the obs location
call quad_lon_lat_locate(interp_handle, lon_lat_vert(1), lon_lat_vert(2), &
                         four_lons, four_lats, lon_fract, lat_fract, status1)
if (status1 /= 0) then
   istatus(:) = 3  ! cannot locate enclosing horizontal quad
   return
endif

! if we are avoiding assimilating obs above a given pressure, test here and return.
if (discarding_high_obs) then
   call obs_too_high(lon_lat_vert(3), which_vert, status1)
   if (status1 /= 0) then
      istatus(:) = status1
      return
   endif
endif

call get_quad_vals(state_handle, ens_size, varid, obs_qty, four_lons, four_lats, &
                   lon_lat_vert, which_vert, quad_vals, status_array)

!>@todo FIXME : Here we are failing if any ensemble member fails. Instead
!>              we should be using track status...
if (any(status_array /= 0)) then
   istatus(:) = maxval(status_array)   ! cannot get the state values at the corners
   return
endif

! do the horizontal interpolation for each ensemble member
call quad_lon_lat_evaluate(interp_handle, lon_fract, lat_fract, ens_size, &
                           quad_vals, interp_vals, status_array)

if (any(status_array /= 0)) then
   istatus(:) = 8   ! cannot evaluate in the quad
   return
endif

if (using_chemistry) &
   interp_vals = interp_vals * get_volume_mixing_ratio(obs_qty)

! all interp values should be set by now.  set istatus
istatus(:) = 0

end subroutine model_interpolate

!-----------------------------------------------------------------------
!> internal only version of model interpolate. 
!> does not check for locations too high - return all actual values.

subroutine interpolate_values(state_handle, ens_size, location, obs_qty, varid, &
                              interp_vals, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_qty
integer,             intent(in) :: varid
real(r8),           intent(out) :: interp_vals(ens_size) 
integer,            intent(out) :: istatus(ens_size)

character(len=*), parameter :: routine = 'interpolate_values:'

integer  :: which_vert, four_lons(4), four_lats(4)
real(r8) :: lon_fract, lat_fract
real(r8) :: lon_lat_vert(3), quad_vals(4, ens_size)
type(quad_interp_handle) :: interp_handle

interp_vals(:) = MISSING_R8
istatus(:)     = 99

interp_handle = get_interp_handle(obs_qty)
lon_lat_vert  = get_location(location)
which_vert    = nint(query_location(location)) 

call quad_lon_lat_locate(interp_handle, lon_lat_vert(1), lon_lat_vert(2), &
                         four_lons, four_lats, lon_fract, lat_fract, istatus(1))
if (istatus(1) /= 0) then
   istatus(:) = 3  ! cannot locate enclosing horizontal quad
   return
endif

call get_quad_vals(state_handle, ens_size, varid, obs_qty, four_lons, four_lats, &
                   lon_lat_vert, which_vert, quad_vals, istatus)
if (any(istatus /= 0)) return

call quad_lon_lat_evaluate(interp_handle, lon_fract, lat_fract, ens_size, &
                           quad_vals, interp_vals, istatus)
if (any(istatus /= 0)) then
   istatus(:) = 8   ! cannot evaluate in the quad
   return
endif

end subroutine interpolate_values

!-----------------------------------------------------------------------
!> return my_status /= 0 if obs is above a user-defined threshold.
!> intended to be quick (low-cost) and not exact. 

subroutine obs_too_high(vert_value, which_vert, my_status)
real(r8), intent(in) :: vert_value
integer,  intent(in) :: which_vert
integer, intent(out) :: my_status

! assume ok to begin with
my_status = 0

if (which_vert == VERTISPRESSURE) then
   ! lower pressures are higher; watch the less than/greater than tests
   if (vert_value < no_assim_above_pressure) my_status = 14
   return
endif

! these are always ok
if (which_vert == VERTISSURFACE .or. which_vert == VERTISUNDEF) return

if (which_vert == VERTISHEIGHT) then
   if (vert_value > no_assim_above_height) my_status = 14
   return
endif

if (which_vert == VERTISLEVEL) then
   ! level 1 is top; watch less than/greater than in tests
   if (vert_value < no_assim_above_level) my_status = 14
   return
endif

! for now we haven't run into observations where the vertical coordinate
! (of the OBS) is in scale height - but if we do it will fall into here.

write(string2, *) 'vertical type: ', which_vert
call error_handler(E_ERR, 'obs_too_high', 'unrecognized vertical type', &
                   source, revision, revdate, text2=string2)

end subroutine obs_too_high

!-----------------------------------------------------------------------
!>

subroutine get_quad_vals(state_handle, ens_size, varid, obs_qty, four_lons, four_lats, &
                         lon_lat_vert, which_vert, quad_vals, my_status)
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: varid
integer,             intent(in) :: obs_qty
integer,             intent(in) :: four_lons(4)
integer,             intent(in) :: four_lats(4)
real(r8),            intent(in) :: lon_lat_vert(3)
integer,             intent(in) :: which_vert
real(r8),           intent(out) :: quad_vals(4, ens_size) !< array of interpolated values
integer,            intent(out) :: my_status(ens_size)

integer  :: icorner, numdims
integer  :: level_one_array(ens_size)
integer  :: four_levs1(4, ens_size), four_levs2(4, ens_size)
real(r8) :: four_vert_fracts(4, ens_size)

character(len=*), parameter :: routine = 'get_quad_vals:'

quad_vals(:,:) = MISSING_R8
my_status(:) = 99

! need to consider the case for 2d vs 3d variables
numdims = get_dims_from_qty(obs_qty, varid)

! now here potentially we have different results for different
! ensemble members.  the things that can vary are dimensioned by ens_size.

if (numdims == 3) then

   ! build 4 columns to find vertical level numbers
   do icorner=1, 4
      call find_vertical_levels(state_handle, ens_size, &
                                four_lons(icorner), four_lats(icorner), lon_lat_vert(3), &
                                which_vert, obs_qty, varid, &
                                four_levs1(icorner, :), four_levs2(icorner, :), & 
                                four_vert_fracts(icorner, :), my_status)
      if (any(my_status /= 0)) return
  
   enddo
   
   ! we have all the indices and fractions we could ever want.
   ! now get the data values at the bottom levels, the top levels, 
   ! and do vertical interpolation to get the 4 values in the columns.
   ! the final horizontal interpolation will happen later.
      
   if (varid > 0) then

      call get_four_state_values(state_handle, ens_size, four_lons, four_lats, &
                                four_levs1, four_levs2, four_vert_fracts, &   
                                varid, quad_vals, my_status)

   else ! get 3d special variables in another ways ( like QTY_PRESSURE )
      call get_four_nonstate_values(state_handle, ens_size, four_lons, four_lats, &
                                   four_levs1, four_levs2, four_vert_fracts, & 
                                   obs_qty, quad_vals, my_status)

   endif

   if (any(my_status /= 0)) return

else if (numdims == 2) then

   if (varid > 0) then
      level_one_array(:) = 1
      do icorner=1, 4
         call get_values_from_varid(state_handle,  ens_size, & 
                                    four_lons(icorner), four_lats(icorner), &
                                    level_one_array, varid, quad_vals(icorner,:),my_status)
         if (any(my_status /= 0)) return

      enddo

   else ! special 2d case
      do icorner=1, 4
         call get_quad_values(ens_size, four_lons(icorner), four_lats(icorner), &
                               obs_qty, obs_qty, quad_vals(icorner,:))
      enddo
      ! apparently this can't fail
      my_status(:) = 0
      
   endif

else
   write(string1, *) trim(get_name_for_quantity(obs_qty)), ' has dimension ', numdims
   call error_handler(E_ERR, routine, 'only supports 2D or 3D fields', &
                      source, revision, revdate, text2=string1)
endif

! when you get here, my_status() was set either by passing it to a
! subroutine, or setting it explicitly here.  if this routine returns
! the default value of 99 something went wrong in this logic.

end subroutine get_quad_vals

!-----------------------------------------------------------------------
!>

subroutine get_four_state_values(state_handle, ens_size, four_lons, four_lats, &
                                 four_levs1, four_levs2, four_vert_fracts, &
                                 varid, quad_vals, my_status)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: four_lons(4), four_lats(4)
integer,             intent(in) :: four_levs1(4, ens_size), four_levs2(4, ens_size)
real(r8),            intent(in) :: four_vert_fracts(4, ens_size)
integer,             intent(in) :: varid
real(r8),           intent(out) :: quad_vals(4, ens_size) !< array of interpolated values
integer,            intent(out) :: my_status(ens_size)

integer  :: icorner
real(r8) :: vals1(ens_size), vals2(ens_size)

character(len=*), parameter :: routine = 'get_four_state_values:'

do icorner=1, 4
   call get_values_from_varid(state_handle,  ens_size, &
                              four_lons(icorner), four_lats(icorner), &
                              four_levs1(icorner, :), varid, vals1, &
                              my_status)

   if (any(my_status /= 0)) then
      my_status(:) = 16   ! cannot retrieve vals1 values
      return
   endif

   call get_values_from_varid(state_handle,  ens_size, &
                              four_lons(icorner), four_lats(icorner), &
                              four_levs2(icorner, :), varid, vals2, my_status)
   if (any(my_status /= 0)) then
      my_status(:) = 17   ! cannot retrieve top values
      return
   endif

   call vert_interp(ens_size, vals1, vals2, four_vert_fracts(icorner, :), & 
                    quad_vals(icorner, :))

enddo


end subroutine get_four_state_values

!-----------------------------------------------------------------------
!>

subroutine get_four_nonstate_values(state_handle, ens_size, four_lons, four_lats, &
                                 four_levs1, four_levs2, four_vert_fracts, &
                                 obs_qty, quad_vals, my_status)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: four_lons(4), four_lats(4)
integer,             intent(in) :: four_levs1(4, ens_size), four_levs2(4, ens_size)
real(r8),            intent(in) :: four_vert_fracts(4, ens_size)
integer,             intent(in) :: obs_qty
real(r8),           intent(out) :: quad_vals(4, ens_size) !< array of interpolated values
integer,            intent(out) :: my_status(ens_size)

integer  :: icorner
real(r8) :: vals1(ens_size), vals2(ens_size)

character(len=*), parameter :: routine = 'get_four_nonstate_values:'

do icorner=1, 4
   call get_values_from_nonstate_fields(state_handle,  ens_size, &
                              four_lons(icorner), four_lats(icorner), &
                              four_levs1(icorner, :), obs_qty, vals1, my_status)
   if (any(my_status /= 0)) then
      my_status(:) = 16   ! cannot retrieve vals1 values
      return
   endif

   call get_values_from_nonstate_fields(state_handle,  ens_size, &
                              four_lons(icorner), four_lats(icorner), &
                              four_levs2(icorner, :), obs_qty, vals2, my_status)
   if (any(my_status /= 0)) then
      my_status(:) = 17   ! cannot retrieve top values
      return
   endif

   call vert_interp(ens_size, vals1, vals2, four_vert_fracts(icorner, :), &
                    quad_vals(icorner, :))

enddo

end subroutine get_four_nonstate_values

!-----------------------------------------------------------------------
!> figure out whether this is a 2d or 3d field based on the quantity.
!> if this field is in the state vector, use the state routines.
!> if it's not, there are cases for known other quantities we can
!> interpolate and return.  add any new non-state fields here.

function get_dims_from_qty(obs_quantity, var_id)
integer, intent(in) :: obs_quantity
integer, intent(in) :: var_id
integer :: get_dims_from_qty

character(len=*), parameter :: routine = 'get_dims_from_qty:'

if (var_id > 0) then
   get_dims_from_qty = get_num_dims(domain_id,var_id)
else
   select case (obs_quantity)
      case (QTY_SURFACE_ELEVATION)
         get_dims_from_qty = 2
      case (QTY_PRESSURE, QTY_GEOMETRIC_HEIGHT)
         get_dims_from_qty = 3
      case default 
         write(string1, *) 'we can not interpolate qty "', get_name_for_quantity(obs_quantity), &
                           '" if the dimension is not known'
         call error_handler(E_ERR,routine, string1,source,revision,revdate)
    end select
endif

end function get_dims_from_qty

!-----------------------------------------------------------------------
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
   case (QTY_SURFACE_ELEVATION, &
         QTY_PRESSURE,          &
         QTY_GEOMETRIC_HEIGHT,  &
         QTY_VERTLEVEL) 
      my_status = 0
   case default
      my_status = 2
end select


end subroutine ok_to_interpolate


!-----------------------------------------------------------------------
!>
!>  This is for 2d special observations quantities not in the state

subroutine get_quad_values(ens_size, lon_index, lat_index, obs_quantity, stagger_qty, vals)
integer,  intent(in) :: ens_size
integer,  intent(in) :: lon_index
integer,  intent(in) :: lat_index
integer,  intent(in) :: obs_quantity
integer,  intent(in) :: stagger_qty
real(r8), intent(out) :: vals(ens_size) 

character(len=*), parameter :: routine = 'get_quad_values'

integer :: stagger, prev_lon, next_lat
real(r8) :: vals1(ens_size), vals2(ens_size)

stagger = grid_stagger%qty_stagger(stagger_qty)

select case (obs_quantity)
   case (QTY_SURFACE_ELEVATION)

     select case (stagger)
       case (STAGGER_U)
          call quad_index_neighbors(lon_index, lat_index, prev_lon, next_lat)
          vals1(:) = phis(lon_index, lat_index) 
          vals2(:) = phis(lon_index, next_lat) 
     
        vals = (vals1 + vals2) * 0.5_r8 
     
       case (STAGGER_V)
          call quad_index_neighbors(lon_index, lat_index, prev_lon, next_lat)
          vals1(:) = phis(lon_index, lat_index) 
          vals2(:) = phis(prev_lon,  lat_index) 
     
        vals = (vals1 + vals2) * 0.5_r8
     
       ! no stagger - cell centers, or W stagger
       case default
  
        vals = phis(lon_index, lat_index)
  
     end select
    
     !>@todo FIXME:
     ! should this be using gravity at the given latitude? 
     vals = vals / gravity

   case default 
      write(string1, *) 'we can not interpolate qty', obs_quantity
      call error_handler(E_ERR,routine,string1,source,revision,revdate)

end select

end subroutine get_quad_values


!-----------------------------------------------------------------------
!> interpolate in the vertical between 2 arrays of items.
!>
!> vert_fracts: 0 is 100% of the first level and 
!>              1 is 100% of the second level

subroutine vert_interp(nitems, levs1, levs2, vert_fracts, out_vals)
integer,  intent(in)  :: nitems
real(r8), intent(in)  :: levs1(nitems)
real(r8), intent(in)  :: levs2(nitems)
real(r8), intent(in)  :: vert_fracts(nitems)
real(r8), intent(out) :: out_vals(nitems)

out_vals(:) = (levs1(:) * (1.0_r8-vert_fracts(:))) + &
              (levs2(:) *         vert_fracts(:))

end subroutine vert_interp

!-----------------------------------------------------------------------
!> given lon/lat indices, add one to lat and subtract one from lon
!> check for wraparound in lon, and north pole at lat.
!> intent is that you give the indices into the staggered grid
!> and the return values are the indices in the original unstaggered
!> grid that you need to compute the midpoints for the staggers.
!>@todo FIXME this needs a picture or ascii art

subroutine quad_index_neighbors(lon_index, lat_index, prev_lon, next_lat)
integer, intent(in)  :: lon_index
integer, intent(in)  :: lat_index
integer, intent(out) :: prev_lon
integer, intent(out) :: next_lat

next_lat = lat_index+1
if (next_lat > grid_data%lat%nsize) next_lat = grid_data%lat%nsize

prev_lon = lon_index-1
if (prev_lon < 1) prev_lon = grid_data%lon%nsize

end subroutine quad_index_neighbors


!-----------------------------------------------------------------------
!> given a lon/lat index number, a quantity and a vertical value and type,
!> return which two levels these are between and the fraction across.
!> 

subroutine find_vertical_levels(ens_handle, ens_size, lon_index, lat_index, vert_val, &
                                which_vert, obs_qty, var_id, levs1, levs2, vert_fracts, my_status)
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: lon_index 
integer,             intent(in)  :: lat_index
real(r8),            intent(in)  :: vert_val
integer,             intent(in)  :: which_vert
integer,             intent(in)  :: obs_qty
integer,             intent(in)  :: var_id
integer,             intent(out) :: levs1(ens_size)
integer,             intent(out) :: levs2(ens_size)
real(r8),            intent(out) :: vert_fracts(ens_size)
integer,             intent(out) :: my_status(ens_size)

character(len=*), parameter :: routine = 'find_vertical_levels:'

integer  :: l1, l2, imember, level_one, status1, k
real(r8) :: fract1
real(r8) :: surf_pressure (  ens_size )
real(r8) :: pressure_array( ref_nlevels, ens_size )
real(r8) :: height_array  ( ref_nlevels, ens_size )

! assume the worst
levs1(:)    = MISSING_I
levs2(:)    = MISSING_I
vert_fracts(:) = MISSING_R8
my_status(:)   = 98

! ref_nlevels is the number of vertical levels (midlayer points)

level_one = 1

select case (which_vert)

   case(VERTISPRESSURE)
      ! construct a pressure column here and find the model levels
      ! that enclose this value
      call get_staggered_values_from_qty(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                                         lon_index, lat_index, level_one, obs_qty, &
                                         surf_pressure, status1)
      if (status1 /= 0) then
         my_status(:) = status1
         return
      endif

      call build_cam_pressure_columns(ens_size, surf_pressure, ref_nlevels, pressure_array)

      do imember=1, ens_size
         call pressure_to_level(ref_nlevels, pressure_array(:, imember), vert_val, & 
                                levs1(imember), levs2(imember), &
                                vert_fracts(imember), my_status(imember))

      enddo

      if (debug_level > 100) then
         do k = 1,ens_size
            print*, 'ISPRESSURE levs1(k), levs2(k), vert_fracts(k), vert_val', &
                     levs1(k), levs2(k), vert_fracts(k), vert_val
          enddo
      endif

   case(VERTISHEIGHT)
      ! construct a height column here and find the model levels
      ! that enclose this value
      call cam_height_levels(ens_handle, ens_size, lon_index, lat_index, ref_nlevels, obs_qty, &
                             height_array, my_status)

      !>@todo FIXME let successful members continue?
      if (any(my_status /= 0)) return

      if (debug_level > 400) then
         do k = 1,ref_nlevels
            print*, 'ISHEIGHT: ', k, height_array(k,1)
         enddo
      endif

      do imember=1, ens_size
         call height_to_level(ref_nlevels, height_array(:, imember), vert_val, & 
                             levs1(imember), levs2(imember), vert_fracts(imember), &
                             my_status(imember))
      enddo

      !>@todo FIXME let successful members continue?
      if (any(my_status /= 0)) return

      if (debug_level > 100) then
         do k = 1,ens_size
            print*, 'ISHEIGHT ens#, levs1(#), levs2(#), vert_fracts(#), top/bot height(#)', &
                     k, levs1(k), levs2(k), vert_fracts(k), height_array(levs2(k),k), height_array(levs1(k), k)
         enddo
      endif
      
   case(VERTISLEVEL)
      ! this routine returns false if the level number is out of range.
      if (.not. check_good_levels(vert_val, ref_nlevels, l1, l2, fract1)) then
         my_status(:) = 8
         return
      endif

      ! because we're given a model level as input, all the ensemble
      ! members have the same outgoing values.
      levs1(:) = l1
      levs2(:) = l2
      vert_fracts(:) = fract1
      my_status(:) = 0

      if (debug_level > 100) then
         do k = 1,ens_size
            print*, 'ISLEVEL levs1(k), levs2(k), vert_fracts(k), vert_val', &
                     levs1(k), levs2(k), vert_fracts(k), vert_val
         enddo
      endif

   ! 2d fields
   case(VERTISUNDEF, VERTISSURFACE)
      if (get_dims_from_qty(obs_qty, var_id) == 2) then
         levs1(:) = ref_nlevels - 1
         levs2(:) = ref_nlevels
         vert_fracts(:) = 1.0_r8
         my_status(:) = 0
      else
         my_status(:) = 4 ! can not get vertical levels
      endif

   case default
      write(string1, *) 'unsupported vertical type: ', which_vert
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
      
end select

! by this time someone has already set my_status(), good or bad.

end subroutine find_vertical_levels

!-----------------------------------------------------------------------
!> Compute the heights at pressure midpoints
!>
!> this version does all ensemble members at once.

subroutine cam_height_levels(ens_handle, ens_size, lon_index, lat_index, nlevels, qty, height_array, my_status) 
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: lon_index
integer,             intent(in)  :: lat_index
integer,             intent(in)  :: nlevels
integer,             intent(in)  :: qty
real(r8),            intent(out) :: height_array(nlevels, ens_size)
integer,             intent(out) :: my_status(ens_size)

integer  :: k, level_one, imember, status1
real(r8) :: surface_elevation(1)
real(r8) :: surface_pressure(ens_size), mbar(nlevels, ens_size)
real(r8) :: tv(nlevels, ens_size)  ! Virtual temperature, top to bottom

! this is for surface obs
level_one = 1

! get the surface pressure from the ens_handle
call get_staggered_values_from_qty(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                                   lon_index, lat_index, level_one, qty, surface_pressure, status1)

! get the surface elevation from the phis, including stagger if needed
call get_quad_values(1, lon_index, lat_index, QTY_SURFACE_ELEVATION, qty, surface_elevation)

call compute_virtual_temperature(ens_handle, ens_size, lon_index, lat_index, nlevels, qty, tv, status1)

if (status1 /= 0) then
   my_status = status1
   return
endif

if (use_variable_mean_mass) then
   call compute_mean_mass(ens_handle, ens_size, lon_index, lat_index, nlevels, qty, mbar, status1)

   if (status1 /= 0) then
      my_status = status1
      return
   endif

   ! compute the height columns for each ensemble member - passing mbar() array in.
   do imember = 1, ens_size
      call build_heights(nlevels, surface_pressure(imember), surface_elevation(1), &
                         tv(:, imember), height_array(:, imember), mbar=mbar(:, imember))
   enddo

else
   ! compute the height columns for each ensemble member - no mbar() argument here.
   ! (you cannot just pass 1.0 in for the mbar() array; it uses a different gas constant
   ! in the variable mean mass case.)
   do imember = 1, ens_size
      call build_heights(nlevels, surface_pressure(imember), surface_elevation(1), &
                         tv(:, imember), height_array(:, imember))
   enddo
endif


if (debug_level > 100) then
 do imember = 1, ens_size
  print *, ''
  print *, 'geopotential, member: ', imember
  do k = 1, nlevels
    print*, 'tv(level)    ', k, tv(k, imember)
  enddo
  do k = 1, nlevels
    print*, 'height(level)', k, height_array(k, imember)
  enddo
 enddo
endif

! convert entire array to geometric height (from potential height)
call gph2gmh(height_array, grid_data%lat%vals(lat_index))

if (debug_level > 100) then
 do imember = 1, ens_size
  print *, ''
  print *, 'geometric, member: ', imember
  do k = 1, nlevels
    print*, 'height(level)', k, height_array(k, imember)
  enddo
 enddo
endif

my_status(:) = 0

end subroutine cam_height_levels

!-----------------------------------------------------------------------
!> Compute the pressures at the layer midpoints for multiple columns

subroutine build_cam_pressure_columns(ens_size, surface_pressure, n_levels, pressure_array)

integer,            intent(in)  :: ens_size
real(r8),           intent(in)  :: surface_pressure(:)   ! in pascals
integer,            intent(in)  :: n_levels
real(r8),           intent(out) :: pressure_array(:,:)

integer :: j, k

! Set midpoint pressures.  This array mirrors the order of the
! cam model levels: 1 is the model top, N is the bottom.

do j=1, ens_size
   do k=1,n_levels
      pressure_array(k, j) = ref_surface_pressure    * grid_data%hyam%vals(k) + &
                                 surface_pressure(j) * grid_data%hybm%vals(k) 
   enddo
enddo

end subroutine build_cam_pressure_columns


!-----------------------------------------------------------------------
!> Compute column of pressures at the layer midpoints for the given 
!> surface pressure.  
!>
!> to get pressure on layer interfaces, the computation would be identical
!> but use hyai, hybi.  (also have n_levels+1)

subroutine single_pressure_column(surface_pressure, n_levels, pressure_array)

real(r8),           intent(in)  :: surface_pressure   ! in pascals
integer,            intent(in)  :: n_levels
real(r8),           intent(out) :: pressure_array(n_levels)

integer :: k

! Set midpoint pressures.  This array mirrors the order of the
! cam model levels: 1 is the model top, N is the bottom.

do k=1, n_levels
   pressure_array(k) = ref_surface_pressure * grid_data%hyam%vals(k) + &
                           surface_pressure * grid_data%hybm%vals(k)
enddo

end subroutine single_pressure_column

!-----------------------------------------------------------------------
!> Compute pressure at one level given the surface pressure
!> cam model levels: 1 is the model top, N is the bottom.
!> in this version of the routine level is integer/whole value

function single_pressure_value_int(surface_pressure, level)

real(r8), intent(in)  :: surface_pressure   ! in pascals
integer,  intent(in)  :: level
real(r8) :: single_pressure_value_int

! cam model levels: 1 is the model top, N is the bottom.

single_pressure_value_int = ref_surface_pressure * grid_data%hyam%vals(level) + &
                                surface_pressure * grid_data%hybm%vals(level)

end function single_pressure_value_int

!-----------------------------------------------------------------------
!> Compute pressure at one level given the surface pressure
!> cam model levels: 1 is the model top, N is the bottom.
!> fraction = 0 is full level 1, fraction = 1 is full level 2
!> level is real/fractional value


function single_pressure_value_real(surface_pressure, level)

real(r8), intent(in)  :: surface_pressure   ! in pascals
real(r8), intent(in)  :: level
real(r8) :: single_pressure_value_real

integer :: k
real(r8) :: fract, pres1, pres2

k = int(level)
fract = level - int(level)

if (k /= ref_nlevels) then
   pres1 = single_pressure_value_int(surface_pressure, k)
   pres2 = single_pressure_value_int(surface_pressure, k+1)
else
   pres1 = single_pressure_value_int(surface_pressure, k-1)
   pres2 = single_pressure_value_int(surface_pressure, k)
   fract = 1.0_r8
endif

single_pressure_value_real = (pres1 * (1.0_r8 - fract)) + &
                              pres2 * (fract)

end function single_pressure_value_real

!-----------------------------------------------------------------------
!> return the level indices and fraction across the level.
!> level 1 is model top, level N is model bottom. 
!> pressure is smallest at the top, so the values are not inverted
!> in the array.
!> fract = 0 means full lev1 value,
!> fract = 1 means full lev2 value. 
!> return non-zero if value outside valid range.

subroutine pressure_to_level(nlevels, pressures, p_val, &
                              lev1, lev2, fract, my_status)

integer,  intent(in)  :: nlevels
real(r8), intent(in)  :: pressures(:)
real(r8), intent(in)  :: p_val
integer,  intent(out) :: lev1
integer,  intent(out) :: lev2
real(r8), intent(out) :: fract  
integer,  intent(out) :: my_status

call find_enclosing_indices(nlevels, pressures, p_val, lev1, lev2, fract, my_status, & 
                            inverted = .false., log_scale = use_log_vertical_scale)

if (my_status /= 0) my_status = 10

end subroutine pressure_to_level

!-----------------------------------------------------------------------
!> return the level indices and fraction across the level.
!> level 1 is model top, level N is model bottom. 
!> height is largest at the top, so the values *are* inverted
!> in the array.
!> fract = 0 means full lev1 value,
!> fract = 1 means full lev2 value. 
!> return non-zero if value outside valid range.

subroutine height_to_level(nlevels, heights, h_val, &
                            lev1, lev2, fract, my_status)

integer,  intent(in)  :: nlevels
real(r8), intent(in)  :: heights(:)
real(r8), intent(in)  :: h_val
integer,  intent(out) :: lev1
integer,  intent(out) :: lev2
real(r8), intent(out) :: fract  
integer,  intent(out) :: my_status

character(len=*), parameter :: routine = 'height_to_level:'

call find_enclosing_indices(nlevels, heights, h_val, lev1, lev2, fract, my_status, &
                            inverted = .true., log_scale = .false.)

if (my_status /= 0) my_status = 11

end subroutine height_to_level

!-----------------------------------------------------------------------
!> in cam level 1 is at the model top, level N is the lowest level
!> our convention in this code is:  between levels a fraction of 0
!> is 100% level 1, and fraction of 1 is 100% level 2.

function check_good_levels(vert_value, valid_range, l1, l2, fract)
real(r8), intent(in)  :: vert_value
integer,  intent(in)  :: valid_range
integer,  intent(out) :: l1
integer,  intent(out) :: l2
real(r8), intent(out) :: fract
logical               :: check_good_levels

integer :: integer_level
real(r8) :: fract_level

! be a pessimist, then you're never disappointed
check_good_levels = .false.
l1 = MISSING_I
l2 = MISSING_I
fract = MISSING_R8

! out of range checks
if (vert_value < 1.0_r8 .or. vert_value > valid_range) return

integer_level = floor(vert_value) 
fract_level = vert_value - integer_level 

! cam levels start at the top so level 1 is
! the highest level and increases on the way down.

!>might want to allow extrapolation - which means
!>allowing out of range values here and handling
!>them correctly in the calling and vert_interp() code.

if (vert_value /= valid_range) then
   l1 = integer_level
   l2 = integer_level + 1
   fract = fract_level
else   
   ! equal to the largest level number
   l1 = integer_level - 1
   l2 = integer_level
   fract = 1.0_r8
endif

check_good_levels = .true.

end function check_good_levels


!-----------------------------------------------------------------------
!> based on the stagger that corresponds to the given quantity,
!> return the handle to the interpolation grid


function get_interp_handle(obs_quantity)
integer, intent(in)      :: obs_quantity
type(quad_interp_handle) :: get_interp_handle

character(len=*), parameter :: routine = 'get_interp_handle:'

select case (grid_stagger%qty_stagger(obs_quantity))
   case ( STAGGER_U ) 
      get_interp_handle = interp_u_staggered
   case ( STAGGER_V ) 
      get_interp_handle = interp_v_staggered
   case ( STAGGER_NONE )
      get_interp_handle = interp_nonstaggered
   case ( STAGGER_W ) 
      write(string1,*) 'w stagger -- not supported yet'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   case ( STAGGER_UV ) 
      write(string1,*) 'uv stagger -- not supported yet'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   case default
      write(string1,*) 'unknown stagger -- this should never happen'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
end select
                      
end function get_interp_handle

!-----------------------------------------------------------------------
!>
!> Set the desired minimum model advance time. This is generally NOT the
!> dynamical timestep of the model, but rather the shortest forecast length
!> you are willing to make. This impacts how frequently the observations
!> may be assimilated.
!>
!>

function shortest_time_between_assimilations()

character(len=*), parameter :: routine = 'shortest_time_between_assimilations:'

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = set_time(assimilation_period_seconds, &
                                               assimilation_period_days)

write(string1,*)'assimilation period is ',assimilation_period_days,   ' days ', &
                                          assimilation_period_seconds,' seconds'
call error_handler(E_MSG,routine,string1,source,revision,revdate)

end function shortest_time_between_assimilations




!-----------------------------------------------------------------------
!>
!> Does any shutdown and clean-up needed for model.
!>

subroutine end_model()

! deallocate arrays from grid and anything else

call free_cam_grid(grid_data)

call free_std_atm_tables()

call finalize_quad_interp(interp_nonstaggered)
call finalize_quad_interp(interp_u_staggered)
call finalize_quad_interp(interp_v_staggered)

if (using_chemistry) call finalize_chem_tables()

end subroutine end_model


!-----------------------------------------------------------------------
!>
!> Writes the model-specific attributes to a DART 'diagnostic' netCDF file.
!> This includes coordinate variables and some metadata, but NOT the
!> actual DART state.
!>
!> @param ncid    the netCDF handle of the DART diagnostic file opened by
!>                assim_model_mod:init_diag_output

subroutine nc_write_model_atts(ncid, dom_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: dom_id    ! not used since there is only one domain

!----------------------------------------------------------------------
! local variables 
!----------------------------------------------------------------------

character(len=*), parameter :: routine = 'nc_write_model_atts'

if ( .not. module_initialized ) call static_init_model

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call nc_begin_define_mode(ncid, routine)

call nc_add_global_creation_time(ncid, routine)

call nc_add_global_attribute(ncid, "model_source", source, routine)
call nc_add_global_attribute(ncid, "model_revision", revision, routine)
call nc_add_global_attribute(ncid, "model_revdate", revdate, routine)

call nc_add_global_attribute(ncid, "model", "CAM", routine)

! this option is for users who want the smallest output
! or diagnostic files - only the state vector data will
! be written.   otherwise, if you want to plot this data
! the rest of this routine writes out enough grid info
! to make the output file look like the input.
if (suppress_grid_info_in_output) then
   call nc_end_define_mode(ncid, routine)
   return
endif

!----------------------------------------------------------------------------
! Output the grid variables.
!----------------------------------------------------------------------------
! Define the new dimensions IDs
!----------------------------------------------------------------------------

call nc_define_dimension(ncid, 'lon',  grid_data%lon%nsize,  routine)
call nc_define_dimension(ncid, 'lat',  grid_data%lat%nsize,  routine)
call nc_define_dimension(ncid, 'slon', grid_data%slon%nsize, routine)
call nc_define_dimension(ncid, 'slat', grid_data%slat%nsize, routine)
call nc_define_dimension(ncid, 'lev',  grid_data%lev%nsize,  routine)
call nc_define_dimension(ncid, 'ilev', grid_data%ilev%nsize, routine)
call nc_define_dimension(ncid, 'gw',   grid_data%gw%nsize,   routine)
call nc_define_dimension(ncid, 'hyam', grid_data%hyam%nsize, routine)
call nc_define_dimension(ncid, 'hybm', grid_data%hybm%nsize, routine)
call nc_define_dimension(ncid, 'hyai', grid_data%hyai%nsize, routine)
call nc_define_dimension(ncid, 'hybi', grid_data%hybi%nsize, routine)

!----------------------------------------------------------------------------
! Create the Coordinate Variables and the Attributes
! The contents will be written in a later block of code.
!----------------------------------------------------------------------------

! U,V Grid Longitudes
call nc_define_real_variable(     ncid, 'lon', (/ 'lon' /),                 routine)
call nc_add_attribute_to_variable(ncid, 'lon', 'long_name', 'longitude',    routine)
call nc_add_attribute_to_variable(ncid, 'lon', 'units',     'degrees_east', routine)


call nc_define_real_variable(     ncid, 'slon', (/ 'slon' /),                       routine)
call nc_add_attribute_to_variable(ncid, 'slon', 'long_name', 'staggered longitude', routine)
call nc_add_attribute_to_variable(ncid, 'slon', 'units',     'degrees_east',        routine)

! U,V Grid Latitudes
call nc_define_real_variable(     ncid, 'lat', (/ 'lat' /),                  routine)
call nc_add_attribute_to_variable(ncid, 'lat', 'long_name', 'latitude',      routine)
call nc_add_attribute_to_variable(ncid, 'lat', 'units',     'degrees_north', routine)


call nc_define_real_variable(     ncid, 'slat', (/ 'slat' /),                      routine)
call nc_add_attribute_to_variable(ncid, 'slat', 'long_name', 'staggered latitude', routine)
call nc_add_attribute_to_variable(ncid, 'slat', 'units',     'degrees_north',      routine)

! Vertical Grid Latitudes
call nc_define_real_variable(     ncid, 'lev', (/ 'lev' /),                                                     routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'long_name',      'hybrid level at midpoints (1000*(A+B))',      routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'units',          'hPa',                                         routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'positive',       'down',                                        routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'standard_name',  'atmosphere_hybrid_sigma_pressure_coordinate', routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'formula_terms',  'a: hyam b: hybm p0: P0 ps: PS',               routine)


call nc_define_real_variable(     ncid, 'ilev', (/ 'ilev' /),                                                    routine)
call nc_add_attribute_to_variable(ncid, 'ilev', 'long_name',      'hybrid level at interfaces (1000*(A+B))',     routine)
call nc_add_attribute_to_variable(ncid, 'ilev', 'units',          'hPa',                                         routine)
call nc_add_attribute_to_variable(ncid, 'ilev', 'positive',       'down',                                        routine)
call nc_add_attribute_to_variable(ncid, 'ilev', 'standard_name',  'atmosphere_hybrid_sigma_pressure_coordinate', routine)
call nc_add_attribute_to_variable(ncid, 'ilev', 'formula_terms',  'a: hyai b: hybi p0: P0 ps: PS',               routine)

! Hybrid Coefficients
call nc_define_real_variable(     ncid, 'hyam', (/ 'lev' /),                                            routine)
call nc_add_attribute_to_variable(ncid, 'hyam', 'long_name', 'hybrid A coefficient at layer midpoints', routine)

call nc_define_real_variable(     ncid, 'hybm', (/ 'lev' /),                                            routine)
call nc_add_attribute_to_variable(ncid, 'hybm', 'long_name', 'hybrid B coefficient at layer midpoints', routine)


call nc_define_real_variable(     ncid, 'hyai', (/ 'ilev' /),                                            routine)
call nc_add_attribute_to_variable(ncid, 'hyai', 'long_name', 'hybrid A coefficient at layer interfaces', routine)


call nc_define_real_variable(     ncid, 'hybi', (/ 'ilev' /),                                            routine)
call nc_add_attribute_to_variable(ncid, 'hybi', 'long_name', 'hybrid B coefficient at layer interfaces', routine)

! Gaussian Weights
call nc_define_real_variable(     ncid, 'gw', (/ 'lat' /),                  routine)
call nc_add_attribute_to_variable(ncid, 'gw', 'long_name', 'gauss weights', routine)

call nc_define_real_scalar(       ncid, 'P0', routine)
call nc_add_attribute_to_variable(ncid, 'P0', 'long_name', 'reference pressure', routine)
call nc_add_attribute_to_variable(ncid, 'P0', 'units',     'Pa',                 routine)

! Finished with dimension/variable definitions, must end 'define' mode to fill.

call nc_end_define_mode(ncid, routine)

!----------------------------------------------------------------------------
! Fill the coordinate variables
!----------------------------------------------------------------------------

call nc_put_variable(ncid, 'lon',  grid_data%lon%vals,  routine)
call nc_put_variable(ncid, 'lat',  grid_data%lat%vals,  routine)
call nc_put_variable(ncid, 'slon', grid_data%slon%vals, routine)
call nc_put_variable(ncid, 'slat', grid_data%slat%vals, routine)
call nc_put_variable(ncid, 'lev',  grid_data%lev%vals,  routine)
call nc_put_variable(ncid, 'ilev', grid_data%ilev%vals, routine)
call nc_put_variable(ncid, 'gw',   grid_data%gw%vals,   routine)
call nc_put_variable(ncid, 'hyam', grid_data%hyam%vals, routine)
call nc_put_variable(ncid, 'hybm', grid_data%hybm%vals, routine)
call nc_put_variable(ncid, 'hyai', grid_data%hyai%vals, routine)
call nc_put_variable(ncid, 'hybi', grid_data%hybi%vals, routine)
call nc_put_variable(ncid, 'P0',   grid_data%P0%vals,   routine)

! flush any pending i/o to disk
call nc_synchronize_file(ncid, routine)

end subroutine nc_write_model_atts

!-----------------------------------------------------------------------
!> writes CAM's model date and time of day into file.  CAM uses
!> integer date values and interger time of day measured in seconds
!>
!> @param ncid         name of the file
!> @param model_time   the current time of the model state
!>

subroutine write_model_time(ncid, model_time)
integer,         intent(in) :: ncid
type(time_type), intent(in) :: model_time

integer :: iyear, imonth, iday, ihour, iminute, isecond
integer :: cam_date(1), cam_tod(1)

character(len=*), parameter :: routine = 'write_model_time'

if ( .not. module_initialized ) call static_init_model

call get_date(model_time, iyear, imonth, iday, ihour, iminute, isecond)

cam_date = iyear*10000 + imonth*100 + iday
cam_tod  = ihour*3600  + iminute*60 + isecond

! if the file doesn't already have a "date" variable make one
if (.not. nc_variable_exists(ncid, "date")) then
   call nc_begin_define_mode(ncid, routine)
   call nc_define_integer_variable(ncid, 'date', (/ 'time' /), routine)
   call nc_end_define_mode(ncid, routine)
   call nc_put_variable(ncid, 'date', cam_date, routine)
endif

! if the file doesn't already have a "datesec" variable make one
if (.not. nc_variable_exists(ncid, "datesec")) then
   call nc_begin_define_mode(ncid, routine)
   call nc_define_integer_variable(ncid, 'datesec', (/ 'time' /), routine)
   call nc_end_define_mode(ncid, routine)
   call nc_put_variable(ncid, 'datesec', cam_tod,  routine)
endif

end subroutine write_model_time

!--------------------------------------------------------------------
!>
!> Read the time from the input file
!>
!> @param filename name of file that contains the time
!>

function read_model_time(filename)

character(len=*), intent(in) :: filename
type(time_type)              :: read_model_time

integer :: ncid
integer :: cam_date, cam_tod
integer :: iyear, imonth, iday, ihour, imin, isec, rem

character(len=*), parameter :: routine = 'read_model_time'

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) trim(filename), ' does not exist.'
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
endif

ncid = nc_open_file_readonly(filename, routine)

! CAM initial files have two variables of length 
! 'time' (the unlimited dimension): date, datesec

call nc_get_variable(ncid, 'date',    cam_date, routine)
call nc_get_variable(ncid, 'datesec', cam_tod,  routine)

! 'date' is YYYYMMDD 
! 'cam_tod' is seconds of current day
iyear  = cam_date / 10000
rem    = cam_date - iyear*10000
imonth = rem / 100
iday   = rem - imonth*100

ihour  = cam_tod / 3600
rem    = cam_tod - ihour*3600
imin   = rem / 60
isec   = rem - imin*60

! some cam files are from before the start of the gregorian calendar.
! since these are 'arbitrary' years, just change the offset.
if (iyear < 1601) then
   write(string1,*)' '
   write(string2,*)'WARNING - ',trim(filename),' changing year from ', &
                   iyear,'to',iyear+1601

   call error_handler(E_MSG, routine, string1, source, revision, &
                      revdate, text2=string2,text3='to make it a valid Gregorian date.')

   write(string1,*)' '
   call error_handler(E_MSG, routine, string1, source, revision)
   iyear = iyear + 1601
endif

read_model_time = set_date(iyear,imonth,iday,ihour,imin,isec)

call nc_close_file(ncid, routine)

end function read_model_time

!--------------------------------------------------------------------
!> if the namelist is set to not use this custom routine, the default
!> dart routine will add 'pert_amp' of noise to every field in the state
!> to generate an ensemble from a single member.  if it is set to true
!> this routine will be called.  the pert_amp will be ignored, and the
!> given list of quantities will be perturbed by the given amplitude
!> (which can be different for each field) to generate an ensemble.

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)
type(ensemble_type), intent(inout) :: state_ens_handle 
integer,             intent(in)    :: ens_size
real(r8),            intent(in)    :: pert_amp   ! ignored in this version
logical,             intent(out)   :: interf_provided

type(random_seq_type) :: seq

integer :: iloc, jloc, vloc, myqty
integer :: max_qtys, j

integer(i8) :: i, state_items
integer(i8), allocatable :: my_vars(:)

logical,  allocatable :: do_these_qtys(:)
real(r8), allocatable :: perturb_by(:)

character(len=*), parameter :: routine = 'pert_model_copies:'

! set by namelist to select using the default routine in filter
! (adds the same noise to all parts of the state vector)
! or the code here that lets you specify which fields get perturbed.
if (custom_routine_to_generate_ensemble) then
   interf_provided = .true.
else
   interf_provided = .false.
   return
endif

! make sure each task is using a different random sequence
call init_random_seq(seq, my_task_id())

max_qtys = get_num_quantities()
allocate(do_these_qtys(0:max_qtys), perturb_by(0:max_qtys))

do_these_qtys(:) = .false.
perturb_by(:)    = 0.0_r8

! this loop is over the number of field names/perturb values
! in the namelist.  it quits when it finds a blank field name.
do i=1, MAX_PERT
   if (fields_to_perturb(i) == '') exit
 
   myqty = get_index_for_quantity(fields_to_perturb(i))
   if (myqty < 0) then
      string1 = 'unrecognized quantity name in "fields_to_perturb" list: ' // &
                trim(fields_to_perturb(i))
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   endif

   do_these_qtys(myqty) = .true.
   perturb_by(myqty)    = perturbation_amplitude(i)
enddo

! get the global index numbers of the part of the state that 
! we have in this task.  here is an example of how to work with
! just the part of the state that is on the current task.
state_items = get_my_num_vars(state_ens_handle)
allocate(my_vars(state_items))
call get_my_vars(state_ens_handle, my_vars)

! this loop is over all the subset of the state items 
! that are on this MPI task.
do i=1, state_items

   ! for each global index number in the state vector find
   ! what quantity it is. (iloc,jloc,vloc are unused here)
   call get_model_variable_indices(my_vars(i), iloc, jloc, vloc, kind_index=myqty)

   ! if myqty is in the namelist, perturb it.  otherwise cycle
   if (.not. do_these_qtys(myqty)) cycle
  
   ! this loop is over the number of ensembles
   do j=1, ens_size
      state_ens_handle%copies(j, i) = random_gaussian(seq, state_ens_handle%copies(j, i), perturb_by(myqty))
   enddo

enddo

deallocate(my_vars)
deallocate(do_these_qtys, perturb_by)

end subroutine pert_model_copies


!-----------------------------------------------------------------------
! The remaining (private) interfaces come last.
! None of the private interfaces need to call static_init_model()
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!>
!> Fill the array of requested variables, dart kinds, possible min/max
!> values and whether or not to update the field in the output file.
!> Then calls 'add_domain()' to tell the DART code which variables to
!> read into the state vector after this code returns.
!>
!>@param variable_array  the list of variables and kinds from model_mod_nml
!>@param nfields         the number of variable/Quantity pairs specified

subroutine set_cam_variable_info( variable_array, nfields )

character(len=*), intent(in)  :: variable_array(:)
integer,          intent(out) :: nfields

character(len=*), parameter :: routine = 'set_cam_variable_info:'

integer :: i
integer, parameter :: MAX_STRING_LEN = 128

character(len=MAX_STRING_LEN) :: varname    ! column 1, NetCDF variable name
character(len=MAX_STRING_LEN) :: dartstr    ! column 2, DART Quantity
character(len=MAX_STRING_LEN) :: minvalstr  ! column 3, Clamp min val
character(len=MAX_STRING_LEN) :: maxvalstr  ! column 4, Clamp max val
character(len=MAX_STRING_LEN) :: updatestr  ! column 5, Update output or not

character(len=vtablenamelength) :: var_names(MAX_STATE_VARIABLES) = ' '
logical  :: update_list(MAX_STATE_VARIABLES)   = .FALSE.
integer  ::   kind_list(MAX_STATE_VARIABLES)   = MISSING_I
real(r8) ::  clamp_vals(MAX_STATE_VARIABLES,2) = MISSING_R8


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

! CAM only has a single domain (only a single grid, no nests or multiple grids)

domain_id = add_domain(cam_template_filename, nfields, var_names, kind_list, &
                       clamp_vals, update_list)

call fill_cam_stagger_info(grid_stagger)

if (debug_level > 100) call state_structure_info(domain_id)

end subroutine set_cam_variable_info


!-----------------------------------------------------------------------
!>
!> Fill the qty_stagger array to tell what type of stagger each variable 
!> has. This will be useful for interpolating observations.
!> This currently doesn't support both slon/slat stagger - but cam-fv 
!> doesn't have any fields like that.
!>

subroutine fill_cam_stagger_info(stagger)
type(cam_stagger), intent(inout) :: stagger

integer :: ivar, jdim, qty_index

allocate(stagger%qty_stagger(0:get_num_quantities()))

stagger%qty_stagger = STAGGER_NONE

do ivar = 1, get_num_variables(domain_id)
   do jdim = 1, get_num_dims(domain_id, ivar)

      if (get_dim_name(domain_id, ivar, jdim) == 'slat') then
         qty_index = get_kind_index(domain_id, ivar) 
         stagger%qty_stagger(qty_index) = STAGGER_U
      endif

      if (get_dim_name(domain_id, ivar, jdim) == 'slon') then
         qty_index = get_kind_index(domain_id, ivar)
         stagger%qty_stagger(qty_index) = STAGGER_V
      endif

      if (get_dim_name(domain_id, ivar, jdim) == 'ilev') then
         qty_index = get_kind_index(domain_id, ivar)
         stagger%qty_stagger(qty_index) = STAGGER_W
      endif

   enddo
enddo

end subroutine fill_cam_stagger_info


!-----------------------------------------------------------------------
!> Read in the grid information from the given CAM restart file.
!> Note that none of the data will be used from this file; just the
!> grid size and locations.  Also read in the elevation information
!> from the "PHIS' file.

subroutine read_grid_info(grid_file, grid)
character(len=*), intent(in)  :: grid_file
type(cam_grid),   intent(out) :: grid

! Get the grid info plus additional non-state arrays
call get_cam_grid(grid_file, grid)

! This non-state variable is used to compute surface elevation.
call read_cam_phis_array(cam_phis_filename)

! Set up the interpolation structures for later 
call setup_interpolation(grid)

end subroutine read_grid_info


!-----------------------------------------------------------------------
!> Read the data from the various cam grid arrays 
!>
!>@todo FIXME not all of these are used.  can we either
!> not read them in, or make them optional?  this does affect
!> what we can write out in the diagnostic file.  if we have
!> to have them in the diag files then we have to read them all
!> even if we never use them.  both ilev and gw currently fall
!> into this category.
!> 

subroutine get_cam_grid(grid_file, grid)
character(len=*), intent(in)  :: grid_file
type(cam_grid), intent(out) :: grid

character(len=*), parameter :: routine = 'get_cam_grid:'

integer :: ncid

! put this in a subroutine that deals with the grid
ncid = nc_open_file_readonly(grid_file, routine)

call fill_cam_1d_array(ncid, 'lon',  grid%lon)
call fill_cam_1d_array(ncid, 'lat',  grid%lat)
call fill_cam_1d_array(ncid, 'lev',  grid%lev)
call fill_cam_1d_array(ncid, 'ilev', grid%ilev) ! for staggered vertical grid
call fill_cam_1d_array(ncid, 'slon', grid%slon)
call fill_cam_1d_array(ncid, 'slat', grid%slat)
call fill_cam_1d_array(ncid, 'gw',   grid%gw)   ! gauss weights
call fill_cam_1d_array(ncid, 'hyai', grid%hyai)
call fill_cam_1d_array(ncid, 'hybi', grid%hybi)
call fill_cam_1d_array(ncid, 'hyam', grid%hyam)
call fill_cam_1d_array(ncid, 'hybm', grid%hybm)

! P0 is a scalar with no dimensionality
call fill_cam_0d_array(ncid, 'P0',   grid%P0) 

call nc_close_file(ncid, routine)

end subroutine get_cam_grid


!-----------------------------------------------------------------------
!>
!> allocate space for a scalar variable and read values into the grid_array
!>   


subroutine fill_cam_1d_array(ncid, varname, grid_array)
integer,            intent(in)    :: ncid
character(len=*),   intent(in)    :: varname
type(cam_1d_array), intent(inout) :: grid_array

character(len=*), parameter :: routine = 'fill_cam_1d_array'

!>@todo do we need to check that this exists?  if all cam input
!> files will have all the arrays we are asking for, then no.

call nc_get_variable_size(ncid, varname, grid_array%nsize)
allocate(grid_array%vals(grid_array%nsize))

call nc_get_variable(ncid, varname, grid_array%vals, routine)

if (debug_level > 80) call array_dump(grid_array%vals, label=varname)

end subroutine fill_cam_1d_array


!-----------------------------------------------------------------------
!>
!> allocate space for a scalar variable and read values into the grid_array
!>   


subroutine fill_cam_0d_array(ncid, varname, grid_array)
integer,            intent(in)    :: ncid
character(len=*),   intent(in)    :: varname
type(cam_1d_array), intent(inout) :: grid_array

character(len=*), parameter :: routine = 'fill_cam_0d_array'

grid_array%nsize = 1
allocate(grid_array%vals(grid_array%nsize))

call nc_get_variable(ncid, varname, grid_array%vals, routine)

if (debug_level > 80) print*, 'variable name ', trim(varname), grid_array%vals

end subroutine fill_cam_0d_array

!-----------------------------------------------------------------------
!>
!> free space in the various grid arrays
!> 

subroutine free_cam_grid(grid)

type(cam_grid), intent(inout) :: grid

call free_cam_1d_array(grid%lon)
call free_cam_1d_array(grid%lat)
call free_cam_1d_array(grid%lev)
call free_cam_1d_array(grid%ilev) 
call free_cam_1d_array(grid%slon)
call free_cam_1d_array(grid%slat)
call free_cam_1d_array(grid%gw)
call free_cam_1d_array(grid%hyai)
call free_cam_1d_array(grid%hybi)
call free_cam_1d_array(grid%hyam)
call free_cam_1d_array(grid%hybm)

call free_cam_1d_array(grid%P0) 

deallocate(phis)

end subroutine free_cam_grid


!-----------------------------------------------------------------------
!>
!> 
!>   

subroutine free_cam_1d_array(grid_array)
type(cam_1d_array), intent(inout) :: grid_array

deallocate(grid_array%vals)
grid_array%nsize = -1

end subroutine free_cam_1d_array

!-----------------------------------------------------------------------
!> convert from string to integer, and set in the dart code the
!> vertical type we are going to want to localize in.
!> 

subroutine set_vert_localization(typename)
character(len=*), intent(in)  :: typename

character(len=*), parameter :: routine = 'set_vert_localization'

character(len=32) :: ucasename
integer :: vcoord

ucasename = typename
call to_upper(ucasename)

select case (ucasename)
  case ("PRESSURE")
     vcoord = VERTISPRESSURE
  case ("HEIGHT")
     vcoord = VERTISHEIGHT
  case ("SCALEHEIGHT", "SCALE_HEIGHT", "SCALE HEIGHT")
     vcoord = VERTISSCALEHEIGHT
  case ("LEVEL", "MODEL_LEVEL", "MODEL LEVEL")
     vcoord = VERTISLEVEL
  case default
     write(string1,*)'unrecognized vertical localization coordinate type: '//trim(typename)
     write(string2,*)'valid values are: PRESSURE, HEIGHT, SCALEHEIGHT, LEVEL'
     call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2)
end select
 
! during assimilation, when get_close() is called to compute the separation distance
! between items, convert all state and obs to use this vertical type if vertical localization 
! is enabled (usually true for cam).

call set_vertical_localization_coord(vcoord)

! save in module global for later use.
vertical_localization_type = vcoord

end subroutine set_vert_localization

!-----------------------------------------------------------------------
!>
!> 
!>   

subroutine setup_interpolation(grid)
type(cam_grid), intent(in) :: grid

!>@todo FIXME the cam fv grid is really evenly spaced in lat and lon,
!>even though they provide full lon() and lat() arrays.  providing the deltas
!>between each pair would be slightly faster inside the interp code.

!print *, 'setting up interpolation: lon/lat sizes = ', grid%lon%nsize, grid%lat%nsize,  &
!                                                       grid%slon%nsize, grid%slat%nsize

! mass points at cell centers
call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, grid%lon%nsize, grid%lat%nsize, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_nonstaggered)
call set_quad_coords(interp_nonstaggered, grid%lon%vals, grid%lat%vals)

! U stagger
call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, grid%lon%nsize, grid%slat%nsize, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_u_staggered)
call set_quad_coords(interp_u_staggered, grid%lon%vals, grid%slat%vals)

! V stagger
call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, grid%slon%nsize, grid%lat%nsize, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_v_staggered)
call set_quad_coords(interp_v_staggered, grid%slon%vals, grid%lat%vals)

end subroutine setup_interpolation

!-----------------------------------------------------------------------
!>
!> 

subroutine read_cam_phis_array(phis_filename)
character(len=*),   intent(in)    :: phis_filename

character(len=*), parameter :: routine = 'read_cam_phis_array'

integer :: ncid, nsize(3)   ! lon, lat, time

ncid = nc_open_file_readonly(phis_filename, routine)

call nc_get_variable_size(ncid, 'PHIS', nsize(:), routine)
allocate( phis(nsize(1), nsize(2)) )

call nc_get_variable(ncid, 'PHIS', phis, routine)

call nc_close_file(ncid, routine)

end subroutine read_cam_phis_array


!-----------------------------------------------------------------------
!> Compute the virtual temperature at the midpoints
!>
!> this version does all ensemble members at once.
!>

subroutine compute_virtual_temperature(ens_handle, ens_size, lon_index, lat_index, nlevels, qty, tv, istatus)

type(ensemble_type), intent(in)   :: ens_handle
integer,             intent(in)   :: ens_size
integer,             intent(in)   :: lon_index
integer,             intent(in)   :: lat_index
integer,             intent(in)   :: nlevels
integer,             intent(in)   :: qty
real(r8),            intent(out)  :: tv(nlevels, ens_size)
integer,             intent(out)  :: istatus

integer :: k
real(r8) :: temperature(ens_size), specific_humidity(ens_size)

!>@todo this should come from a model specific constant module.
!> the forward operators and model_mod should use it.
real(r8), parameter :: rd = 287.05_r8 ! dry air gas constant
real(r8), parameter :: rv = 461.51_r8 ! wet air gas constant
real(r8), parameter :: rr_factor = (rv/rd) - 1.0_r8


! construct a virtual temperature column, one for each ensemble member
do k = 1, nlevels
   ! temperature
   call get_staggered_values_from_qty(ens_handle, ens_size, QTY_TEMPERATURE, &
                                     lon_index, lat_index, k, qty, temperature, istatus)

   if (istatus < 0) return

   ! specific humidity
   call get_staggered_values_from_qty(ens_handle, ens_size, QTY_SPECIFIC_HUMIDITY, &
                                     lon_index, lat_index, k, qty, specific_humidity, istatus)
   if (istatus < 0) return

   !>tv == virtual temperature.
   tv(k,:) = temperature(:)*(1.0_r8 + rr_factor*specific_humidity(:))
   !print*, 'tv(levels)', k,tv(k,1), temperature(1), specific_humidity(1)
enddo


end subroutine compute_virtual_temperature


!-----------------------------------------------------------------------
!> loop through all levels to get the mean mass.
!>

subroutine compute_mean_mass(ens_handle, ens_size, lon_index, lat_index, nlevels, qty, mbar, istatus)
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: lon_index
integer,             intent(in)  :: lat_index
integer,             intent(in)  :: nlevels
integer,             intent(in)  :: qty
real(r8),            intent(out) :: mbar(nlevels, ens_size)
integer,             intent(out) :: istatus

integer :: k, this_qty
real(r8) :: mmr_o1(ens_size, nlevels), &
            mmr_o2(ens_size, nlevels), &
            mmr_h1(ens_size, nlevels), &
            mmr_n2(ens_size, nlevels)
real(r8) :: O_molar_mass, O2_molar_mass, H_molar_mass, N2_molar_mass 

! do this outside the subroutine?  it never changes throughout the
! run of the program
O_molar_mass  = get_molar_mass(QTY_ATOMIC_OXYGEN_MIXING_RATIO)
O2_molar_mass = get_molar_mass(QTY_MOLEC_OXYGEN_MIXING_RATIO)
H_molar_mass  = get_molar_mass(QTY_ATOMIC_H_MIXING_RATIO)
N2_molar_mass = get_molar_mass(QTY_NITROGEN)
   


! High topped models (WACCM-X) need to account for the changing composition 
! of the atmosphere with height.  This requires several variables from the
! initial file, which may not be available from low topped models.
do k = 1, nlevels

   this_qty = QTY_ATOMIC_OXYGEN_MIXING_RATIO
   call get_staggered_values_from_qty(ens_handle, ens_size, this_qty, &
                                      lon_index, lat_index, k, qty, mmr_o1(:, k), istatus)
   if (istatus /= 0) return
   !print *, 'mmr: ', trim(get_name_for_quantity(this_qty)), mmr_o1(1, k)
   
   this_qty = QTY_MOLEC_OXYGEN_MIXING_RATIO
   call get_staggered_values_from_qty(ens_handle, ens_size, this_qty, & 
                                      lon_index, lat_index, k, qty, mmr_o2(:, k), istatus)
   if (istatus /= 0) return
   !print *, 'mmr: ', trim(get_name_for_quantity(this_qty)), mmr_o2(1, k)
   
   this_qty = QTY_ATOMIC_H_MIXING_RATIO
   call get_staggered_values_from_qty(ens_handle, ens_size, this_qty, &
                                      lon_index, lat_index, k, qty, mmr_h1(:, k), istatus)
   if (istatus /= 0) return
   !print *, 'mmr: ', trim(get_name_for_quantity(this_qty)), mmr_h1(1, k)
   
   mmr_n2(:,k) = 1.0_r8 - (mmr_o1(:,k) + mmr_o2(:,k) + mmr_h1(:,k))
   mbar(k,:) = 1.0_r8/( mmr_o1(:,k)/O_molar_mass  &
                      + mmr_o2(:,k)/O2_molar_mass &
                      + mmr_h1(:,k)/H_molar_mass  &
                      + mmr_n2(:,k)/N2_molar_mass)
enddo

end subroutine compute_mean_mass

!-----------------------------------------------------------------------
!> This code is using a finite difference method to evaluate an 
!> integral to solve the hydrostatic equation. 
!>
!> The details are in the reference given below.
!> Don't change this code until you have read the paper and
!> understand what they're doing.  The paper uses a matrix
!> while this code gets away with ignoring 'l' and evaluating
!> the 'k' vector directly. 
!>
!> Equation references are to "Hybrid Coordinates for CCM1"
!>    https://opensky.ucar.edu/islandora/object/technotes%3A149/datastream/PDF/view
!>
!> Here is a comment from the NCL function that does the
!> same thing for them.
!>
!> Purpose:
!>   To compute geopotential height using the CCM2 hybrid coordinate
!>   vertical slice.  Since the vertical integration matrix is a
!>   function of latitude and longitude, it is not explicitly
!>   computed as for sigma coordinates.  The integration algorithm
!>   is derived from Boville's mods in the ibm file hybrid 1mods
!>   (6/17/88).  All vertical slice arrays are oriented top to
!>   bottom as in CCM2.  This field is on full model levels (aka
!>   "midpoints") not half levels.
!>
!> careful - if the calling code passes in the mbar() parameter a different gas
!>           constant is used instead.  an mbar() array of 1.0 is not the same 
!>           as no parameter specified.

subroutine build_heights(nlevels,p_surf,h_surf,virtual_temp,height_midpts,height_interf,mbar)

integer,  intent(in)  :: nlevels                            ! Number of vertical levels
real(r8), intent(in)  :: p_surf                             ! Surface pressure (pascals)
real(r8), intent(in)  :: h_surf                             ! Surface height (m)
real(r8), intent(in)  :: virtual_temp( nlevels)             ! Virtual Temperature
real(r8), intent(out) :: height_midpts(nlevels)             ! Geopotential height at midpoints, top to bottom
real(r8), intent(out), optional :: height_interf(nlevels+1) ! Geopotential height at interfaces, top to bottom
real(r8), intent(in),  optional :: mbar(nlevels)            ! Factor to support for variable gas constant

! Local variables
!>@todo FIXME can we use the types_mod values here?  or have a model constants module?
real(r8), parameter :: const_r = 287.04_r8    ! Different than model_heights (dry air gas constant)
real(r8), parameter :: universal_gas_constant = 8314.0_r8 ! [J/K/kmol]
real(r8), parameter :: g0 = 9.80616_r8        ! Different than model_heights (gph2gmh:G) !

integer  :: k,l

! an array now: real(r8), parameter :: rbyg=r/g0
real(r8) :: pterm(nlevels)   ! vertical scratch space, to improve computational efficiency
real(r8) :: r_g0_tv(nlevels) ! rbyg=r/g0 * tv
real(r8) :: pm_ln(nlevels+1) ! logs of midpoint pressures plus surface interface pressure

! cam uses a uniform gas constant value, but high top
! models like waccm change the gas constant with height.
! allow for the calling code to pass in an array of r.

! if mbar() array is present notice that the units are different
! for the gas constants, so an mbar() array of 1.0s will NOT give
! the same results as if it isn't present.

if (present(mbar)) then
   r_g0_tv(:) = (universal_gas_constant / (mbar(:)*g0)) * virtual_temp(:)
else
   r_g0_tv(:) = (const_r / g0) * virtual_temp(:)
endif

! calculate the log of the pressure column midpoints.
! items 1:nlevels are the midpoints, but NOTICE THAT
! the pressure at nlevels+1 is the pressure of the 
! actual surface interface, not a midpoint!!

call single_pressure_column(p_surf, nlevels, pm_ln)
   
pm_ln(nlevels+1) = p_surf * grid_data%hybi%vals(nlevels+1)   ! surface interface
   
where (pm_ln >  0.0_r8) 
   pm_ln = log(pm_ln) 
else where (pm_ln <= 0.0_r8)
   pm_ln = 0
end where

!debug
!200 format (I3, 6(1X, F24.16))
!201 format (A, 1X, I3, 6(1X, F24.16))
!202 format (A, 6(1X, F24.16))
!203 format (6(1X, F24.16))
!
!print *, 'pm_ln: '
!do i=1, nlevels+1
!  write(*, 200) i, pm_ln(i)
!enddo
!end debug

!        height_midpts(1)=top  ->  height_midpts(nlevels)=bottom
! 
! level
! 1/2    ---------------------------------------------------------------
! 1      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - top
! 3/2    ---------------------------------------------------------------
!
!                 ---------------------------------------------
!        --------/                                             \--------
!                - - - - - - - - - - - - - - - - - - - - - - - - 
! NL     - - - /                                                 \ - - - bottom
!              ---------------------------------------------------
! NL+1/2 -----/|||||||||||||||||||||||||||||||||||||||||||||||||||\-----


! now the finite differences.
! Eq 3.a.109 has 5 piecewise (PW) terms.  The numbers below refer to each of these terms
! in the order they're listed in the paper.

!
! See 2nd PW term here: Eq 3.a.109  where l=K,k<K  h(k,l) = 1/2 * ln [  p(k+1) / p(k) ]
!

do k = 1,nlevels - 1
   height_midpts(k) = h_surf + r_g0_tv(k) * 0.5_r8 * (pm_ln(k+1)-pm_ln(k))
enddo
height_midpts(nlevels) = h_surf + r_g0_tv(nlevels) * (pm_ln(nlevels+1)-pm_ln(nlevels))

!
! See 4th PW term here: Eq 3.a.109  where l=K,k<K  h(k,l) = 1/2*ln[pi*pi/(p(k-1)*p(k))
!

do k = 1,nlevels - 1
    height_midpts(k) = height_midpts(k) + r_g0_tv(nlevels) * &
                       (pm_ln(nlevels+1) - 0.5_r8*(pm_ln(nlevels-1)+pm_ln(nlevels)))
enddo

!
! See 3rd PW term here:  Eqs 1.14 & 3.a.109 where l>K, k<K
!                                h(k,l) = 1/2 * ln [ p(l+1)/p(l-1) ]

! don't recompute the same values multiple times;
! compute once and put into a temporary array.
! (see the double nested loops below with k and l)

! this is really a matrix multiply with a upper triangular
! matrix so it simplifies to a doubly nested loop.

! pterm(1) and (nlevels) are never used, but to prevent
! confusion when debugging set them to 0 so they don't
! look like uninitialized variables.
pterm(1)       = 0.0_r8
do k = 2,nlevels - 1
   pterm(k) = r_g0_tv(k) * 0.5_r8 * (pm_ln(k+1)-pm_ln(k-1))
enddo
pterm(nlevels) = 0.0_r8

do k = 1,nlevels - 2
   do l = k+1, nlevels-1
      height_midpts(k) = height_midpts(k) + pterm(l)
   enddo
enddo

!write(*, 202) 'psurf, hsurf: ', p_surf, h_surf
!do k = 1,nlevels
!   write(*, 201) 'k, height: ', k, height_midpts(k)
!enddo

! not implemented yet.
if (present(height_interf)) then
   height_interf(:) = MISSING_R8
endif

end subroutine build_heights

!-----------------------------------------------------------------------
!>  Convert a 2d array of geopotential altitudes to mean sea level altitudes.
!>  To avoid overflow with very high model tops, convert to km first, compute,
!>  then convert back.  bof.

subroutine gph2gmh(h, lat)
real(r8), intent(inout) :: h(:,:)    ! geopotential altitude in m
real(r8), intent(in)    :: lat       ! latitude in degrees.

real(r8), parameter ::  be = 6356.7516_r8        ! min earth radius, km
real(r8), parameter ::  ae = 6378.1363_r8        ! max earth radius, km
real(r8), parameter ::  G = 0.00980665_r8 ! WMO reference g value, km/s**2, at 45.542N(S)

real(r8) :: g0
real(r8) :: r0
real(r8) :: latr

integer :: i, j

latr = lat * DEG2RAD  ! convert to radians
call compute_surface_gravity(latr, g0)

! compute local earth's radius using ellipse equation

r0 = sqrt( ae**2 * cos(latr)**2 + be**2 * sin(latr)**2)  

! Compute altitude above sea level
do j=1, size(h, 2)
   do i=1, size(h, 1)
      h(i,j) = h(i,j) / 1000.0_r8   ! m to km
      if ( ((g0*r0)/G) - h(i,j) > 0) &
         h(i,j) = (r0 * h(i,j)) / (((g0*r0)/G) - h(i,j))
      h(i,j) = h(i,j) * 1000.0_r8   ! km to m
   enddo
enddo

end subroutine gph2gmh

!-----------------------------------------------------------------------
!> This subroutine computes the Earth's gravity at any latitude.
!> The model assumes the Earth is an oblate spheriod rotating at 
!> the Earth's spin rate.  The model was taken from 
!> "Geophysical Geodesy, Kurt Lambeck, 1988".
!>
!>  input:    xlat, latitude in radians
!>  output:   galt, gravity at the given lat, km/sec**2
!>
!> taken from code from author Bill Schreiner, 5/95
!>
!>

subroutine compute_surface_gravity(xlat, galt)
real(r8), intent(in)  :: xlat
real(r8), intent(out) :: galt

real(r8),parameter :: xmu = 398600.4415_r8         ! km^3/s^2
real(r8),parameter :: ae  = 6378.1363_r8           ! km
real(r8),parameter :: f   = 1.0_r8/298.2564_r8
real(r8),parameter :: xm  = 0.003468_r8            !
real(r8),parameter :: f2  = 5.3481622134089e-03_r8 ! f2 = -f + 5.0* 0.50*xm - 17.0/14.0*f*xm + 15.0/4.0*xm**2
real(r8),parameter :: f4  = 2.3448248012911e-05_r8 ! f4 = -f**2* 0.50 + 5.0* 0.50*f*xm

real(r8) :: g
!real(r8) :: alt = 0.0_r8

! gravity at the equator, km/s2
real(r8), parameter :: ge = xmu/ae**2/(1.0_r8 - f + 1.5_r8*xm - 15.0_r8/14.0_r8*xm*f)


! compute gravity at any latitude, km/s2
g = ge*(1.0_r8 + f2*(sin(xlat))**2 - 1.0_r8/4.0_r8*f4*(sin(2.0_r8*xlat))**2)

! at a fixed altitude of 0.0, g and galt are the same
galt = g

! FIXME: if alt is hardcoded to 0.0, none of this code is needed.
!
! keep it for now in case we want gravity to vary with height.
!
!! compute gravity at any latitude and at any height, km/s2
!galt = g - 2.0_r8*ge*alt/ae*(1.0_r8 + f + xm + (-3.0_r8*f + 5.0_r8* 0.50_r8*xm)*  &
!                          (sin(xlat))**2) + 3.0_r8*ge*alt**2/ae**2
!
!if (g /= galt) print *, 'g, galt: ', g, galt
!
!!! compute gravity at any latitude, km/s2
!!galt = ge*(1.0_r8 + f2*(sin(xlat))**2 - 1.0_r8/4.0_r8*f4*(sin(2.0_r8*xlat))**2)
!
!! convert to meters/s2
!!galt = galt*1000.0_r8

end subroutine compute_surface_gravity

!-----------------------------------------------------------------------
!> This subroutine computes converts vertical state
!>
!>  in:    ens_handle  - mean ensemble handle
!>  in:    num         - number of locations
!>  inout: locs(:)     - locations
!>  in:    loc_qtys(:) - location quantities
!>  in:    loc_indx(:) - location index
!>  in:    which_vert  - vertical location to convert
!>  out:   istatus     - return status 0 is a successful conversion
!>

subroutine convert_vertical_state(ens_handle, num, locs, loc_qtys, loc_indx, &
                                  which_vert, istatus)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer(i8),         intent(in)    :: loc_indx(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus

character(len=*), parameter :: routine = 'convert_vertical_state'

integer :: current_vert_type, ens_size, i

ens_size = 1

do i=1,num
   current_vert_type = nint(query_location(locs(i)))

   if ( current_vert_type == which_vert ) cycle
   if ( current_vert_type == VERTISUNDEF) cycle

   select case (which_vert)
      case (VERTISPRESSURE)
         call state_vertical_to_pressure(    ens_handle, ens_size, locs(i), loc_indx(i), loc_qtys(i) )
      case (VERTISHEIGHT)
         call state_vertical_to_height(      ens_handle, ens_size, locs(i), loc_indx(i), loc_qtys(i) )
      case (VERTISLEVEL)
         call state_vertical_to_level(                   ens_size, locs(i), loc_indx(i), loc_qtys(i) )
      case (VERTISSCALEHEIGHT)
         call state_vertical_to_scaleheight( ens_handle, ens_size, locs(i), loc_indx(i), loc_qtys(i) )
      case default
         write(string1,*)'unable to convert vertical state "', which_vert, '"'
         call error_handler(E_MSG,routine,string1,source,revision,revdate)
   end select
enddo

istatus = 0

end subroutine convert_vertical_state

!--------------------------------------------------------------------

subroutine state_vertical_to_pressure(ens_handle, ens_size, location, location_indx, qty)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx
integer,             intent(in)    :: qty

integer  :: iloc, jloc, vloc, myqty, level_one, status1
integer  :: my_status(ens_size)
real(r8) :: pressure_array(ref_nlevels), surface_pressure(ens_size)


call get_model_variable_indices(location_indx, iloc, jloc, vloc, kind_index=myqty)

if (is_surface_field(myqty)) then
   
   level_one = 1
   call get_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                                     iloc, jloc, level_one, surface_pressure, status1)

   if (status1 /= 0) then
      return
   endif
   call set_vertical(location, surface_pressure(1), VERTISPRESSURE)
else
   call cam_pressure_levels(ens_handle, ens_size, iloc, jloc, ref_nlevels, &
                            qty, pressure_array, my_status)

   call set_vertical(location, pressure_array(vloc), VERTISPRESSURE)
endif

end subroutine state_vertical_to_pressure

!--------------------------------------------------------------------

subroutine state_vertical_to_height(ens_handle, ens_size, location, location_indx, qty)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx
integer,             intent(in)    :: qty

integer  :: iloc, jloc, vloc, my_status(ens_size)
real(r8) :: height_array(ref_nlevels, ens_size)

! build a height column and a pressure column and find the levels
call get_model_variable_indices(location_indx, iloc, jloc, vloc)

call cam_height_levels(ens_handle, ens_size, iloc, jloc, ref_nlevels, &
                       qty, height_array, my_status) 

!>@todo FIXME this can only be used if ensemble size is 1
call set_vertical(location, height_array(vloc,1), VERTISHEIGHT)

end subroutine state_vertical_to_height

!--------------------------------------------------------------------

subroutine state_vertical_to_scaleheight(ens_handle, ens_size, location, location_indx, qty)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx
integer,             intent(in)    :: qty

integer  :: iloc, jloc, vloc, level_one, status1, my_status(ens_size)
real(r8) :: pressure_array(ref_nlevels)
real(r8) :: surface_pressure(1), scaleheight_val

!> this is currently only called with an ensemble size of 1 for
!> vertical conversion.  since it is working only on state variables
!> we don't expect it to ever fail.

level_one = 1
scaleheight_val = MISSING_R8

if (no_normalization_of_scale_heights) then


   if (query_location(location) == VERTISSURFACE) then

      ! get the surface pressure from the ens_handle
      call get_model_variable_indices(location_indx, iloc, jloc, vloc)

      call get_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                                        iloc, jloc, level_one, surface_pressure, status1)
      if (status1 /= 0) goto 200

      scaleheight_val = log(surface_pressure(1))

   else

      ! build a pressure column and and find the levels
      call get_model_variable_indices(location_indx, iloc, jloc, vloc)

      call cam_pressure_levels(ens_handle, ens_size, iloc, jloc, ref_nlevels, &
                               qty, pressure_array, my_status)
      if (any(my_status /= 0)) goto 200
   
      scaleheight_val = log(pressure_array(vloc))

   endif

else

   ! handle surface obs separately here.
   if (query_location(location) == VERTISSURFACE) then

      scaleheight_val = 0.0_r8   ! log(1.0)

   else

      ! build a pressure column and and find the levels
      call get_model_variable_indices(location_indx, iloc, jloc, vloc)

      call cam_pressure_levels(ens_handle, ens_size, iloc, jloc, ref_nlevels, &
                               qty, pressure_array, my_status)
      if (any(my_status /= 0)) goto 200

      ! get the surface pressure from the ens_handle
      call get_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                                        iloc, jloc, level_one, surface_pressure, status1)
      if (status1 /= 0) goto 200
   
      scaleheight_val = scale_height(pressure_array(vloc), surface_pressure(1), no_normalization_of_scale_heights)

   endif

endif
   
200 continue   ! done

call set_vertical(location, scaleheight_val, VERTISSCALEHEIGHT)

end subroutine state_vertical_to_scaleheight

!--------------------------------------------------------------------

subroutine state_vertical_to_level(ens_size, location, location_indx, qty)
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx
integer,             intent(in)    :: qty

integer  :: iloc, jloc, vloc

!>@todo FIXME qty is currently unused.  if we need it, its here.
!>if we really don't need it, we can remove it.  all the other
!>corresponding routines like this use it.  (not clear what to
!>return if field is W or something else with a vertical stagger.)

call get_model_variable_indices(location_indx, iloc, jloc, vloc)

call set_vertical(location, real(vloc,r8), VERTISLEVEL)

end subroutine state_vertical_to_level

!--------------------------------------------------------------------
!> using a standard atmosphere pressure column, convert a height directly to pressure

function generic_height_to_pressure(height, status)
real(r8), intent(in)  :: height
integer,  intent(out) :: status
real(r8) :: generic_height_to_pressure

integer :: lev1, lev2
real(r8) :: fract

generic_height_to_pressure = MISSING_R8

call height_to_level(std_atm_table_len, std_atm_hgt_col, height, & 
                     lev1, lev2, fract, status)
if (status /= 0) return

generic_height_to_pressure = std_atm_pres_col(lev1) * (1.0_r8-fract) + &
                             std_atm_pres_col(lev2) * (fract)

end function generic_height_to_pressure

!--------------------------------------------------------------------
!> using a standard atmosphere pressure column, convert a pressure directly to height

function generic_pressure_to_height(pressure, status)
real(r8), intent(in)  :: pressure
integer,  intent(out) :: status
real(r8) :: generic_pressure_to_height

integer :: lev1, lev2
real(r8) :: fract

generic_pressure_to_height = MISSING_R8

call pressure_to_level(std_atm_table_len, std_atm_pres_col, pressure, &
                       lev1, lev2, fract, status)
if (status /= 0) return

generic_pressure_to_height = std_atm_hgt_col(lev1) * (1.0_r8 - fract) + &
                             std_atm_hgt_col(lev2) * (fract)

end function generic_pressure_to_height

!--------------------------------------------------------------------
!> using the cam eta arrays, convert a pressure directly to model level
!> use P0 as surface, ignore elevation.

function generic_cam_pressure_to_cam_level(pressure, status)
real(r8), intent(in)  :: pressure
integer,  intent(out) :: status
real(r8) :: generic_cam_pressure_to_cam_level

integer :: lev1, lev2
real(r8) :: fract
real(r8) :: pressure_array(ref_nlevels)

generic_cam_pressure_to_cam_level = MISSING_R8

call single_pressure_column(ref_surface_pressure, ref_nlevels, pressure_array)

call pressure_to_level(ref_nlevels, pressure_array, pressure, &
                       lev1, lev2, fract, status)
if (status /= 0) return

generic_cam_pressure_to_cam_level = lev1 + fract

end function generic_cam_pressure_to_cam_level

!-----------------------------------------------------------------------
!> Compute the pressure values at midpoint levels
!>
!> this version does all ensemble members at once.

subroutine cam_pressure_levels(ens_handle, ens_size, lon_index, lat_index, nlevels, qty, &
                               pressure_array, my_status) 
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: lon_index
integer,             intent(in)  :: lat_index
integer,             intent(in)  :: nlevels
integer,             intent(in)  :: qty
real(r8),            intent(out) :: pressure_array(nlevels, ens_size)
integer,             intent(out) :: my_status(ens_size)

integer     :: level_one, status1
real(r8)    :: surface_pressure(ens_size)

! this is for surface obs
level_one = 1

! get the surface pressure from the ens_handle
call get_staggered_values_from_qty(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
        lon_index, lat_index, level_one, qty, surface_pressure, status1)

if (status1 /= 0) then
   my_status(:) = status1
   return
endif

call build_cam_pressure_columns(ens_size, surface_pressure, ref_nlevels, &
                               pressure_array)
my_status(:) = 0

end subroutine cam_pressure_levels

!--------------------------------------------------------------------

subroutine convert_vertical_obs(ens_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, my_status)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer,             intent(in)    :: loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: my_status(:)

character(len=*), parameter :: routine = 'convert_vertical_obs'

integer :: current_vert_type, i

do i=1,num
   current_vert_type = nint(query_location(locs(i)))

   if (( current_vert_type == which_vert ) .or. &
       ( current_vert_type == VERTISUNDEF)) then
      my_status(i) = 0
      cycle
   endif

   select case (which_vert)
      case (VERTISPRESSURE)
         call obs_vertical_to_pressure(   ens_handle, locs(i), my_status(i))
      case (VERTISHEIGHT)
         call obs_vertical_to_height(     ens_handle, locs(i), my_status(i))
      case (VERTISLEVEL)
         call obs_vertical_to_level(      ens_handle, locs(i), my_status(i))
      case (VERTISSCALEHEIGHT)
         call obs_vertical_to_scaleheight(ens_handle, locs(i), my_status(i))
      case default
         write(string1,*)'unable to convert vertical obs "', which_vert, '"'
         call error_handler(E_ERR,routine,string1,source,revision,revdate)
   end select
enddo

end subroutine convert_vertical_obs

!--------------------------------------------------------------------

subroutine obs_vertical_to_pressure(ens_handle, location, my_status)

type(ensemble_type), intent(in)    :: ens_handle
type(location_type), intent(inout) :: location
integer,             intent(out)   :: my_status

integer  :: varid, ens_size, status(1), qty
real(r8) :: pressure_array(ref_nlevels)

character(len=*), parameter :: routine = 'obs_vertical_to_pressure'

ens_size = 1

qty = QTY_PRESSURE
if (query_location(location) == VERTISSURFACE) then
   qty = QTY_SURFACE_PRESSURE
endif

call ok_to_interpolate(qty, varid, my_status)
if (my_status /= 0) return

call interpolate_values(ens_handle, ens_size, location, &
                        qty, varid, pressure_array(:), status(:))


if (status(1) /= 0) then
   my_status = status(1)
   return
endif

call set_vertical(location, pressure_array(1), VERTISPRESSURE)

my_status = 0

end subroutine obs_vertical_to_pressure

!--------------------------------------------------------------------

subroutine obs_vertical_to_height(ens_handle, location, my_status)
type(ensemble_type), intent(in)    :: ens_handle
type(location_type), intent(inout) :: location
integer,             intent(out)   :: my_status

integer  :: varid, ens_size, status(1)
real(r8) :: height_array(1)

character(len=*), parameter :: routine = 'obs_vertical_to_height'

ens_size = 1

call ok_to_interpolate(QTY_GEOMETRIC_HEIGHT, varid, my_status)
if (my_status /= 0) return

call interpolate_values(ens_handle, ens_size, location, &
                        QTY_GEOMETRIC_HEIGHT, varid, height_array(:), status(:))
if (status(1) /= 0) then
   my_status = status(1)
   return
endif

call set_vertical(location, height_array(1), VERTISHEIGHT)

my_status = 0

end subroutine obs_vertical_to_height

!--------------------------------------------------------------------

subroutine obs_vertical_to_level(ens_handle, location, my_status)
type(ensemble_type), intent(in)    :: ens_handle
type(location_type), intent(inout) :: location
integer,             intent(out)   :: my_status

integer  :: varid, ens_size, status(1)
real(r8) :: level_array(1)

ens_size = 1
varid = -1

call interpolate_values(ens_handle, ens_size, location, &
                        QTY_VERTLEVEL, varid, level_array(:), status(:))
if (status(1) /= 0) then
   my_status = status(1)
   return
endif

call set_vertical(location, level_array(1), VERTISLEVEL)

my_status = 0

end subroutine obs_vertical_to_level

!--------------------------------------------------------------------

subroutine obs_vertical_to_scaleheight(ens_handle, location, my_status)
type(ensemble_type), intent(in)    :: ens_handle
type(location_type), intent(inout) :: location
integer,             intent(out)   :: my_status

integer  :: varid1, varid2, ens_size, status(1), ptype
real(r8) :: pressure_array(1), surface_pressure_array(1)
real(r8) :: scaleheight_val

character(len=*), parameter :: routine = 'obs_vertical_to_scaleheight'

ens_size = 1

! there are 4 cases here.

if (no_normalization_of_scale_heights) then

   ! take log of pressure, either surface pressure or regular pressure
   
   if (query_location(location) == VERTISSURFACE) then
      ptype = QTY_SURFACE_PRESSURE
   else
      ptype = QTY_PRESSURE
   endif

   call ok_to_interpolate(ptype, varid1, my_status)
   if (my_status /= 0) return
      
   !>@todo FIXME IFF the obs location is already pressure, we can take it at
   !> face value here and not interpolate it.  however it won't fail if the
   !> pressure here is less than the ensemble mean pressure at this point.
   !> is that ok?
   
   if (ptype == QTY_PRESSURE .and. is_vertical(location, "PRESSURE")) then
      pressure_array(:) = query_location(location, "VLOC")
      my_status = 0
   else
      call interpolate_values(ens_handle, ens_size, location, ptype, varid1, &
                              pressure_array(:), status(:))
      if (status(1) /= 0) then
         my_status = status(1)
         return
      endif
   endif
   
   scaleheight_val = log(pressure_array(1))

else

   ! handle surface obs separately here.
   if (query_location(location) == VERTISSURFACE) then

      scaleheight_val = 0.0_r8  ! -log(1.0)

   else

      call ok_to_interpolate(QTY_PRESSURE, varid1, my_status)
      if (my_status /= 0) return
   
      !>@todo FIXME IFF the obs location is already pressure, we can take it at
      !> face value here and not interpolate it.  however, it can result in negative
      !> scale height values if the pressure is larger than the surface pressure at
      !> that location.  that's what the original cam model_mod did.  is that ok?
   
      if (ptype == QTY_PRESSURE .and. is_vertical(location, "PRESSURE")) then
         pressure_array(:) = query_location(location, "VLOC")
         my_status = 0
      else
         call interpolate_values(ens_handle, ens_size, location, QTY_PRESSURE, varid1, &
                                    pressure_array(:), status(:))
         if (status(1) /= 0) then
            my_status = status(1)
            return
         endif
      endif
                                 
      call ok_to_interpolate(QTY_SURFACE_PRESSURE, varid2, my_status)
      if (my_status /= 0) return
      
      call interpolate_values(ens_handle, ens_size, location, QTY_SURFACE_PRESSURE, varid2, &
                                    surface_pressure_array(:), status(:))
      if (status(1) /= 0) then
         my_status = status(1)
         return
      endif
      
      scaleheight_val = scale_height(pressure_array(1),surface_pressure_array(1), no_normalization_of_scale_heights)

   endif

endif

call set_vertical(location, scaleheight_val, VERTISSCALEHEIGHT)

my_status = 0

end subroutine obs_vertical_to_scaleheight

!--------------------------------------------------------------------

subroutine convert_vert_one_obs(ens_handle, loc, otype, vert_type, status1)
type(ensemble_type), intent(in)    :: ens_handle
type(location_type), intent(inout) :: loc
integer,             intent(in)    :: otype
integer,             intent(in)    :: vert_type
integer,             intent(out)   :: status1

type(location_type) :: bl(1)
integer :: bq(1), bt(1), status(1)

! these need to be arrays.  kinda a pain.
bl(1) = loc
bt(1) = otype
bq(1) = get_quantity_for_type_of_obs(otype)

call convert_vertical_obs(ens_handle, 1, bl, bq, bt, &
                             vert_type, status)
status1 = status(1)
if (status1 /= 0) return

loc = bl(1)

end subroutine convert_vert_one_obs

!--------------------------------------------------------------------

subroutine init_discard_high_obs()

! compute a conversion table between height and pressure based on
! a surface pressure of 1010 mb.  this is a fixed table and does not
! vary with temperature, humidity or surface elevation. 
! use only for quick conversions when absolute accuracy 
! isn't a primary concern.

character(len=*), parameter :: routine = 'init_discard_high_obs'
integer :: my_status

integer :: table_type
character(len=16) :: out_fmt, out_fmt1, pres_fmt
real(r8) :: no_assim_above_scaleh

! pick the better table: 
!  one is more accurate for the lower atmosphere, 
!  and the other has a very high top.
table_type = store_std_atm_tables(ref_model_top_pressure)

! set formatting which is easiest to read in the log.
! the very high top table has very small numbers that need
! exponential notation.
out_fmt  = '(A,F12.5,A)'
out_fmt1 = '(A,I5)'
pres_fmt = out_fmt
if (table_type == HIGH_TOP_TABLE) pres_fmt = '(A,E12.5,A)'

! levels can be fractional but the namelist only allows integer, so simplify the formatting
write(string1, out_fmt1) &
   'Discarding observations higher than model level ', no_obs_assim_above_level
call error_handler(E_MSG, 'init_discard_high_obs', string1, source, revision, revdate)

no_assim_above_pressure = single_pressure_value(ref_surface_pressure, no_obs_assim_above_level)
write(string1, pres_fmt) &
   ' ... which is equivalent to pressure level ', no_assim_above_pressure, ' Pascals' 
call error_handler(E_MSG, 'init_discard_high_obs', string1, source, revision, revdate)

no_assim_above_height = generic_pressure_to_height(no_assim_above_pressure, my_status)
if (my_status /= 0) then
   call error_handler(E_ERR, routine, 'error converting pressure to height', &
                      source, revision, revdate, text2='"no_assim_above_pressure" invalid value')
endif
 
write(string1, out_fmt) &
   ' ... which is equivalent to height         ', no_assim_above_height, ' meters' 
call error_handler(E_MSG, 'init_discard_high_obs', string1, source, revision, revdate)

! special for this - normalize by Ps for printing out
no_assim_above_scaleh = scale_height(no_assim_above_pressure, ref_surface_pressure, .false.)
write(string1, out_fmt) &
   ' ... which is equivalent to scale height   ', no_assim_above_scaleh
call error_handler(E_MSG, 'init_discard_high_obs', string1, source, revision, revdate)

end subroutine init_discard_high_obs

!--------------------------------------------------------------------
! initialize what we can here.  the highest end of the ramp is fixed;
! the start depends on the cutoff distance which can be observation
! type dependent.  at the time the ramping adjustment is applied all
! vertical coordinates will have already been converted to the
! vertical localization type.

subroutine init_damping_ramp_info()

real(r8) :: model_top

character(len=*), parameter :: routine = 'init_damping_ramp_info'

integer :: table_type
character(len=16) :: out_fmt

! pick the better table: 
!  one is more accurate for the lower atmosphere, 
!  and the other has a very high top.
table_type = store_std_atm_tables(ref_model_top_pressure)

! set formatting which is easiest to read in the log.
! the very high top table has very small numbers that need
! exponential notation.
out_fmt = '(A,F12.5,A)'
if (table_type == HIGH_TOP_TABLE .and. &
    vertical_localization_type == VERTISPRESSURE) out_fmt = '(A,E12.5,A)'

! convert to vertical localization units
call convert_vertical_level_generic(real(model_damping_ends_at_level, r8), &
                                         vertical_localization_type, ramp_end, string3, no_norm=.false.)

! check for conversion errors
if (ramp_end == MISSING_R8) then
   write(string1, *) 'error converting ramp_end to vertical localization units'
   call error_handler(E_MSG, routine, 'unexpected error', &
                      source, revision, revdate, text2=string1)
endif

! this value only used for print statement, unused otherwise
call convert_vertical_level_generic(1.0_r8, vertical_localization_type, &
                                    model_top, string3, no_norm=.false.)

! check for conversion errors
if (model_top == MISSING_R8) then
   write(string1, *) 'error converting model_top to vertical localization units'
   call error_handler(E_MSG, routine, 'unexpected error', &
                      source, revision, revdate, text2=string1)
endif

! at this point, ramp_end and model_top are in the localization units

! let the log know what we're doing
write(string1, '(A,I5)') 'Increments will go to 0.0 at model level ', model_damping_ends_at_level
write(string2, out_fmt) 'which is ', ramp_end, ' '//trim(string3)
call error_handler(E_MSG, routine, &
   'Decreasing increments in region damped in the model', &
   string1, source, revision, revdate, text2=string1, text3=string2)

write(string1, out_fmt) 'For reference, model top is ', model_top, ' '//trim(string3)
call error_handler(E_MSG, routine, string1, source, revision, revdate)

end subroutine init_damping_ramp_info

!--------------------------------------------------------------------
!> pressure gets smaller as you go up, everything else gets larger.
!> return true if this value is above the start of the ramp.
!> test_value and ramp_end need to already be in vert localization units

! FIXME: test this new code section carefully.
!
! right now the calling code is expecting extra_dist to be added
! to the original get_dist() value, so any scaling or modifications
! should happen in this routine.
!
! do we need the 2 locations here to compute the horizontal distance?
! or is having the total dist and the vertical separation enough?

function above_ramp_start(test_value, gc, obs_type, ramp_end, total_dist, extra_dist)
real(r8),             intent(in)  :: test_value
type(get_close_type), intent(in)  :: gc
integer,              intent(in)  :: obs_type
real(r8),             intent(in)  :: ramp_end
real(r8),             intent(in)  :: total_dist
real(r8),             intent(out) :: extra_dist
logical :: above_ramp_start

real(r8) :: vert_localize_dist, ramp_start, norm, vert_norm, vert_only_dist
real(r8) :: horiz_dist, ramp_dist, ramp_width
type(location_type) :: this_loc, ramp_start_loc, loc1, loc2
logical, save :: onetime = .true.


! do the easy cases first - either above the ramp end
! or below the ramp start.  leave the middle ground for
! last because we have to then compute a damping factor.

! FIXME: test this!!!
! is it above the ramp end? set damp dist to something
! large enough to turn off all impacts.  is vert_localize_dist enough?
vert_localize_dist = get_maxdist(gc, obs_type)
if (.false. .and. onetime) then
   print *, 'vert_localize_dist = ', vert_localize_dist
   onetime = .false.
endif

if (v_above(test_value, ramp_end)) then
   extra_dist = vert_localize_dist
   above_ramp_start = .true.
   return
endif

! compute ramp start and see if we're lower than that.

! vert norm for this obs type
loc1 = set_location(0.0_r8, 0.0_r8, 0.0_r8, vertical_localization_type)
loc2 = set_location(0.0_r8, 0.0_r8, 1.0_r8, vertical_localization_type)
norm = get_dist(loc1, loc2, obs_type)   ! units: rad/loc units
vert_norm = 1.0_r8 / norm               ! units now: loc units/rad

ramp_start = v_down(ramp_end, vert_norm * vert_localize_dist)

!print *, 'computing ramp start: ramp_end, vert_norm, vert_localize_dist', &
!                    ramp_start, ramp_end, vert_norm, vert_localize_dist

if (.not. v_above(test_value, ramp_start)) then
   extra_dist = 0.0_r8
   above_ramp_start = .false.
   return
endif


! ok, we're somewhere inbetween.  compute horiz and vert distances
! and see what the ramping factor needs to be.

!print *, 'test value within ramp range: ', ramp_start, test_value, ramp_end
above_ramp_start = .true.

! see what the vertical separation is from obs to start of ramp
this_loc       = set_location(0.0_r8, 0.0_r8, test_value, vertical_localization_type)
ramp_start_loc = set_location(0.0_r8, 0.0_r8, ramp_start, vertical_localization_type)

! do we need this?  i think so.   radians
vert_only_dist = get_dist(ramp_start_loc, this_loc, obs_type)

! we need this to compute what?
if (vert_only_dist > total_dist) then
   !print *, 'unexpected, vert larger than total:  ', vert_only_dist, total_dist
   !print *, 'obs_type, vert_norm = ', obs_type, vert_norm
   horiz_dist = 0.0_r8
else
   horiz_dist = sqrt(total_dist**2 - vert_only_dist**2)
endif

ramp_dist  = v_difference(test_value, ramp_start)
ramp_width = v_difference(ramp_end,   ramp_start)
extra_dist = (ramp_dist / ramp_width) * vert_localize_dist

! DEBUG - disable for now
if (.false. .and. above_ramp_start) then
   print *, 'ramp s/v/e: ', ramp_start, test_value, ramp_end
   print *, 'v, h:       ', vert_only_dist, horiz_dist
   print *, 'rampd, tot: ', ramp_dist, ramp_width
   print *, 'ed, return: ', extra_dist, above_ramp_start
endif

end function above_ramp_start

!--------------------------------------------------------------------
! vertical functions - these deal with the fact that pressure,
! scale height, and model levels all get larger as you go from 
! higher in the atmosphere to lower in the atmosphere, but height 
! is the opposite.  these all depend on the global setting of the
! vertical localization type.


!--------------------------------------------------------------------
!> for pressure, level, and one flavor of scale height
!> smaller numbers are further away from the surface.
!> for height and the other flavor of scale height 
!> the opposite is true.  set this once at init time.

subroutine init_sign_of_vert_units()

if (vertical_localization_type == VERTISHEIGHT) then 
   higher_is_smaller = .false.

else if (vertical_localization_type == VERTISSCALEHEIGHT) then 
   ! FIXME: note from nick on scale height:
   !  If no_normalization_of_scale_heights is true, then SH=log(pressure), 
   !  and scale height will decrease with increasing height. 
   !  However, if it is false then SH= -1*log(pressure/surface_pressure) 
   !  and it will increase with increasing height. 

   if (no_normalization_of_scale_heights) then 
      higher_is_smaller = .true.
   else
      higher_is_smaller = .false.
   endif

else
   higher_is_smaller = .true.

endif

end subroutine init_sign_of_vert_units

!--------------------------------------------------------------------
! returns true if a is above b (higher in the atmosphere, 
! further from the surface of the earth).  

pure function v_above(a, b)
real(r8), intent(in) :: a, b
logical :: v_above

if (higher_is_smaller) then
   v_above = (a < b)
else
   v_above = (a > b)
endif

end function v_above

!--------------------------------------------------------------------
! returns new value of moving b distance down in the atmosphere
! starting at a.  for height, this results in a smaller value
! (also one flavor of scale height), but for other vertical types 
! this results in a larger value.

pure function v_down(a, b)
real(r8), intent(in) :: a, b
real(r8) :: v_down

if (higher_is_smaller) then
   v_down = (a + b)
else
   v_down = (a - b)
endif

end function v_down

!--------------------------------------------------------------------
! returns difference of a and b
! (doesn't depend on the vertical_localization_type)

pure function v_difference(a, b)
real(r8), intent(in) :: a, b
real(r8) :: v_difference

v_difference = abs(a - b)

end function v_difference

!--------------------------------------------------------------------
!> this should only be used for converting vertical values which
!> are the same for all ensemble members at all locations. 
!> it uses generic values to do a vertical conversion.

subroutine convert_vertical_level_generic(level_value, want_vert_type, out_value, out_label, no_norm)
real(r8),         intent(in)            :: level_value
integer,          intent(in)            :: want_vert_type
real(r8),         intent(out)           :: out_value
character(len=*), intent(out), optional :: out_label
logical,          intent(in),  optional :: no_norm

character(len=*), parameter :: routine = 'convert_vertical_level_generic'

integer  :: status
real(r8) :: tmp_val
logical  :: no_norm_flag

if (present(no_norm)) then
   no_norm_flag = no_norm
else
   no_norm_flag = no_normalization_of_scale_heights
endif

if (want_vert_type == VERTISLEVEL) then
    out_value = level_value
    if (present(out_label)) out_label = 'levels'
else
   ! convert to the requested units.  start by going to pressure
   tmp_val = single_pressure_value(ref_surface_pressure, level_value)

   select case (want_vert_type)
     case (VERTISPRESSURE)
       out_value = tmp_val
       if (present(out_label)) out_label = 'pascals'
   
     case (VERTISSCALEHEIGHT)
       out_value = scale_height(tmp_val, ref_surface_pressure, no_norm_flag)
       if (present(out_label)) out_label = 'scale heights'
   
     case (VERTISHEIGHT)
       out_value = generic_pressure_to_height(tmp_val, status)
       if (status /= 0) out_value = MISSING_R8
       if (present(out_label)) out_label = 'meters'
   
     case default
       write(string1, *) 'unknown requested vertical type ', want_vert_type
       call error_handler(E_MSG, routine, 'unexpected error', &
                          source, revision, revdate, text2=string1)
   end select
endif

end subroutine convert_vertical_level_generic

!--------------------------------------------------------------------

subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

! The specific type of the base observation, plus the generic kinds list
! for either the state or obs lists are available if a more sophisticated
! distance computation is needed.

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(inout) :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:), loc_types(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ens_handle

character(len=*), parameter :: routine = 'get_close_obs'

integer :: i, status(1), this, vert_type
real(r8) :: vert_value, extra_damping_dist
real(r8), parameter :: LARGE_DIST = 999999.0  ! positive and large

! if absolute distances aren't needed, or vertical localization isn't on,
! the default version works fine since no conversion will be needed and
! there won't be any damping since there are no vert distances.
if (.not. present(dist) .or. .not. vertical_localization_on()) then
   call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                          num_close, close_ind, dist, ens_handle)
   return
endif

if (.not. present(ens_handle)) then
   call error_handler(E_ERR, routine,  &
           'unexpected error: cannot convert distances without an ensemble handle', &
           source, revision, revdate)
endif

! does the base obs need conversion first?
vert_type = query_location(base_loc)

if (vert_type /= vertical_localization_type) then
   call convert_vert_one_obs(ens_handle, base_loc, base_type, &
                             vertical_localization_type, status(1))
   if (status(1) /= 0) then
      num_close = 0
      return
   endif
endif

! FIXME: is here where we need to compute start of ramp for this
! obs type?  should we cache these?  start with doing the computation
! each time, then make an array indexed by obs types with the
! start of the ramp and fill it in on demand.  have to call for
! maxdist(obs_type) and do the math, but just once per type.

! ok, distance is needed and we are localizing in the vertical.
! call default get close to get potentically close locations
! but call without distance so it doesn't do extra work.
call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                       num_close, close_ind)

! compute distances, converting vertical first if need be.
do i=1, num_close
   this = close_ind(i)

   vert_type = query_location(locs(this))
!print *, 'close_o, vval, vtype = ', i, query_location(locs(this), 'VLOC'), vert_type

   if (vert_type /= vertical_localization_type) then
      call convert_vertical_obs(ens_handle, 1, locs(this:this), &
                                loc_qtys(this:this), loc_types(this:this), &
                                vertical_localization_type, status)
      if (status(1) /= 0) then
         dist(i) = LARGE_DIST
         cycle
      endif

   endif

   dist(i) = get_dist(base_loc, locs(this))

   ! do not try to damp impacts when obs has "vert is undefined".  
   ! the impact will go all the way to the model top.  
   ! this is the agreed-on functionality.
   if (.not. are_damping .or. vert_type == VERTISUNDEF) cycle

   vert_value = query_location(locs(this), 'VLOC')
   if (above_ramp_start(vert_value, gc, base_type, ramp_end, dist(i), extra_damping_dist)) then
      dist(i) = dist(i) + extra_damping_dist
   endif
enddo

end subroutine get_close_obs

!----------------------------------------------------------------------------

subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)

! The specific type of the base observation, plus the generic kinds list
! for either the state or obs lists are available if a more sophisticated
! distance computation is needed.

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:)
integer(i8),                   intent(in)  :: loc_indx(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ens_handle

character(len=*), parameter :: routine = 'get_close_state'

integer :: i, status, this, vert_type
real(r8) :: vert_value, extra_damping_dist
real(r8), parameter :: LARGE_DIST = 999999.0  ! positive and large

! if absolute distances aren't needed, or vertical localization isn't on,
! the default version works fine since no conversion will be needed and
! there won't be any damping since there are no vert distances.
if (.not. present(dist) .or. .not. vertical_localization_on()) then
   call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                            num_close, close_ind, dist, ens_handle)
   return
endif

if (.not. present(ens_handle)) then
   call error_handler(E_ERR, routine,  &
           'unexpected error: cannot convert distances without an ensemble handle', &
           source, revision, revdate)
endif

! does the base obs need conversion first?
vert_type = query_location(base_loc)

if (vert_type /= vertical_localization_type) then
   call convert_vert_one_obs(ens_handle, base_loc, base_type, &
                             vertical_localization_type, status)
   if (status /= 0) then
      num_close = 0
      return
   endif
endif

! ok, distance is needed and we are localizing in the vertical.
! call default get close to get potentically close locations
! but call without distance so it doesn't do extra work.
call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                         num_close, close_ind)

! compute distances, converting vertical first if need be.
do i=1, num_close
   this = close_ind(i)

   vert_type = query_location(locs(this))
!print *, 'close_s, vval, vtype = ', i, query_location(locs(this), 'VLOC'), vert_type

   if (vert_type /= vertical_localization_type) then
      call convert_vertical_state(ens_handle, 1, locs(this:this), &
                                 loc_qtys(this:this), loc_indx(this:this), &
                                 vertical_localization_type, status)
      if (status /= 0) then
         dist(i) = LARGE_DIST
         cycle
      endif

   endif

   dist(i) = get_dist(base_loc, locs(this))

   ! do not try to damp impacts when obs has "vert is undefined".  
   ! the impact will go all the way to the model top.  
   ! this is the agreed-on functionality.
   if (.not. are_damping .or. vert_type == VERTISUNDEF) cycle

   vert_value = query_location(locs(this), 'VLOC')
   if (above_ramp_start(vert_value, gc, base_type, ramp_end, dist(i), extra_damping_dist)) then
      dist(i) = dist(i) + extra_damping_dist
   endif
enddo


end subroutine get_close_state

!--------------------------------------------------------------------
!> set values that are used by many routines here and which do not
!> change during the execution of filter.

subroutine init_globals()

ref_surface_pressure = grid_data%P0%vals(1)
ref_model_top_pressure = grid_data%hyai%vals(1) * ref_surface_pressure
ref_nlevels = grid_data%lev%nsize

end subroutine init_globals

!--------------------------------------------------------------------
! Function to calculate scale height given a pressure and optionally
! a surface pressure.  (See the namelist item which controls whether to
! normalize the pressure value aloft with the surface pressure or not.
! We currently only use scale height for computing distances between
! two locations, so the surface pressure terms cancel out - exactly if
! the two locations are co-located horizontally, almost if they are not.
! Normalizing by the surface pressure means in areas of high orography
! the surface differences propagate all the way to the model top.  
! To be backwards-compatible, do this normalization; the current thinking 
! is we shouldn't do it both for scientific reasons and because it 
! doubles the work if it's expensive to find the correct horizontal 
! location, i.e. mpas irregular grids. In this model we always have
! the surface pressure at a location so it's not a performance issue.)
!
! Watch out for unusual cases that could crash the log() function
! We pass in the surface pressure here even if it isn't going to be
! used because in all the cases above we seem to have it (or the standard
! reference pressure) everywhere we are going to compute this value.
! The "skip_norm" parameter controls whether this code uses the 
! surface pressure or not. 

function scale_height(p_above, p_surface, skip_norm)
real(r8), intent(in) :: p_above
real(r8), intent(in) :: p_surface
logical,  intent(in) :: skip_norm
real(r8)             :: scale_height

real(r8), parameter :: tiny = epsilon(1.0_r8)
real(r8) :: diff

if (skip_norm) then
   scale_height = log(p_above)
   return
endif

diff = p_surface - p_above  ! should be positive

if (abs(diff) < tiny) then
   ! surface obs will have (almost) identical values
   scale_height = 0.0_r8   ! -log(1.0_r8)

else if (diff <= tiny .or. p_above <= 0.0_r8) then
   ! weed out bad cases
   scale_height = MISSING_R8

else
   ! normal computation - should be safe now
   scale_height = -log(p_above / p_surface )

endif

end function scale_height

!--------------------------------------------------------------------

function is_surface_field(qty)
integer, intent(in) :: qty
logical :: is_surface_field

is_surface_field = (qty == QTY_SURFACE_PRESSURE .or. qty == QTY_SURFACE_ELEVATION)
   
end function is_surface_field

!-----------------------------------------------------------------------
!> Store a table of pressures and heights. based on a std atmosphere.
!>  not precise - use only when rough numbers are good enough.
!>  return which table was used.
!>
!> table from: http://www.pdas.com/atmos.html
!> and also see:  http://www.pdas.com/upatmos.html
!> for a good explanation of why you can't use the standard
!> equations at high altitudes.   the low tables came from
!> tables.c, and the high one came from bigtables.out.
!> (all found in the atmos.zip file from that web site.)


function store_std_atm_tables(this_model_top)
real(r8), intent(in) :: this_model_top
integer :: store_std_atm_tables

logical, save :: table_initialized = .false.

if (this_model_top < high_top_threshold) then
   if (.not. table_initialized) call load_high_top_table()
   store_std_atm_tables = HIGH_TOP_TABLE
else
   if (.not. table_initialized) call load_low_top_table()
   store_std_atm_tables = LOW_TOP_TABLE
endif

table_initialized = .true.

end function store_std_atm_tables

!-----------------------------------------------------------------------
!> Free arrays associated with generic tables

subroutine free_std_atm_tables()

if (allocated(std_atm_hgt_col))  deallocate(std_atm_hgt_col)
if (allocated(std_atm_pres_col)) deallocate(std_atm_pres_col)

end subroutine free_std_atm_tables

!--------------------------------------------------------------------

subroutine load_low_top_table()
	
std_atm_table_len = 45
allocate(std_atm_hgt_col(std_atm_table_len), std_atm_pres_col(std_atm_table_len))
	
std_atm_hgt_col(1)  = 86.0_r8 ; std_atm_pres_col(1)  = 3.732E-01_r8
std_atm_hgt_col(2)  = 84.0_r8 ; std_atm_pres_col(2)  = 5.308E-01_r8
std_atm_hgt_col(3)  = 82.0_r8 ; std_atm_pres_col(3)  = 7.498E-01_r8
std_atm_hgt_col(4)  = 80.0_r8 ; std_atm_pres_col(4)  = 1.052E+00_r8
std_atm_hgt_col(5)  = 78.0_r8 ; std_atm_pres_col(5)  = 1.467E+00_r8
std_atm_hgt_col(6)  = 76.0_r8 ; std_atm_pres_col(6)  = 2.033E+00_r8
std_atm_hgt_col(7)  = 74.0_r8 ; std_atm_pres_col(7)  = 2.800E+00_r8
std_atm_hgt_col(8)  = 72.0_r8 ; std_atm_pres_col(8)  = 3.835E+00_r8
std_atm_hgt_col(9)  = 70.0_r8 ; std_atm_pres_col(9)  = 5.220E+00_r8
std_atm_hgt_col(10) = 68.0_r8 ; std_atm_pres_col(10) = 7.051E+00_r8
std_atm_hgt_col(11) = 66.0_r8 ; std_atm_pres_col(11) = 9.459E+00_r8
std_atm_hgt_col(12) = 64.0_r8 ; std_atm_pres_col(12) = 1.260E+01_r8
std_atm_hgt_col(13) = 62.0_r8 ; std_atm_pres_col(13) = 1.669E+01_r8
std_atm_hgt_col(14) = 60.0_r8 ; std_atm_pres_col(14) = 2.196E+01_r8
std_atm_hgt_col(15) = 58.0_r8 ; std_atm_pres_col(15) = 2.872E+01_r8
std_atm_hgt_col(16) = 56.0_r8 ; std_atm_pres_col(16) = 3.736E+01_r8
std_atm_hgt_col(17) = 54.0_r8 ; std_atm_pres_col(17) = 4.833E+01_r8
std_atm_hgt_col(18) = 52.0_r8 ; std_atm_pres_col(18) = 6.221E+01_r8
std_atm_hgt_col(19) = 50.0_r8 ; std_atm_pres_col(19) = 7.977E+01_r8
std_atm_hgt_col(20) = 48.0_r8 ; std_atm_pres_col(20) = 1.023E+02_r8
std_atm_hgt_col(21) = 46.0_r8 ; std_atm_pres_col(21) = 1.313E+02_r8
std_atm_hgt_col(22) = 44.0_r8 ; std_atm_pres_col(22) = 1.695E+02_r8
std_atm_hgt_col(23) = 42.0_r8 ; std_atm_pres_col(23) = 2.200E+02_r8
std_atm_hgt_col(24) = 40.0_r8 ; std_atm_pres_col(24) = 2.871E+02_r8
std_atm_hgt_col(25) = 38.0_r8 ; std_atm_pres_col(25) = 3.771E+02_r8
std_atm_hgt_col(26) = 36.0_r8 ; std_atm_pres_col(26) = 4.985E+02_r8
std_atm_hgt_col(27) = 34.0_r8 ; std_atm_pres_col(27) = 6.634E+02_r8
std_atm_hgt_col(28) = 32.0_r8 ; std_atm_pres_col(28) = 8.890E+02_r8
std_atm_hgt_col(29) = 30.0_r8 ; std_atm_pres_col(29) = 1.197E+03_r8
std_atm_hgt_col(30) = 28.0_r8 ; std_atm_pres_col(30) = 1.616E+03_r8
std_atm_hgt_col(31) = 26.0_r8 ; std_atm_pres_col(31) = 2.188E+03_r8
std_atm_hgt_col(32) = 24.0_r8 ; std_atm_pres_col(32) = 2.972E+03_r8
std_atm_hgt_col(33) = 22.0_r8 ; std_atm_pres_col(33) = 4.047E+03_r8
std_atm_hgt_col(34) = 20.0_r8 ; std_atm_pres_col(34) = 5.529E+03_r8
std_atm_hgt_col(35) = 18.0_r8 ; std_atm_pres_col(35) = 7.565E+03_r8
std_atm_hgt_col(36) = 16.0_r8 ; std_atm_pres_col(36) = 1.035E+04_r8
std_atm_hgt_col(37) = 14.0_r8 ; std_atm_pres_col(37) = 1.417E+04_r8
std_atm_hgt_col(38) = 12.0_r8 ; std_atm_pres_col(38) = 1.940E+04_r8
std_atm_hgt_col(39) = 10.0_r8 ; std_atm_pres_col(39) = 2.650E+04_r8
std_atm_hgt_col(40) =  8.0_r8 ; std_atm_pres_col(40) = 3.565E+04_r8
std_atm_hgt_col(41) =  6.0_r8 ; std_atm_pres_col(41) = 4.722E+04_r8
std_atm_hgt_col(42) =  4.0_r8 ; std_atm_pres_col(42) = 6.166E+04_r8
std_atm_hgt_col(43) =  2.0_r8 ; std_atm_pres_col(43) = 7.950E+04_r8
std_atm_hgt_col(44) =  0.0_r8 ; std_atm_pres_col(44) = 1.013E+05_r8
std_atm_hgt_col(45) = -2.0_r8 ; std_atm_pres_col(45) = 1.278E+05_r8

! convert km to m
std_atm_hgt_col(:) = std_atm_hgt_col(:) * 1000.0_r8
	
end subroutine load_low_top_table

!--------------------------------------------------------------------

subroutine load_high_top_table()

std_atm_table_len = 201
allocate(std_atm_hgt_col(std_atm_table_len), std_atm_pres_col(std_atm_table_len))

std_atm_hgt_col(1)   = 1000.0_r8  ;  std_atm_pres_col(1)   = 7.518E-09_r8
std_atm_hgt_col(2)   =  995.0_r8  ;  std_atm_pres_col(2)   = 7.651E-09_r8
std_atm_hgt_col(3)   =  990.0_r8  ;  std_atm_pres_col(3)   = 7.790E-09_r8
std_atm_hgt_col(4)   =  985.0_r8  ;  std_atm_pres_col(4)   = 7.931E-09_r8
std_atm_hgt_col(5)   =  980.0_r8  ;  std_atm_pres_col(5)   = 8.075E-09_r8
std_atm_hgt_col(6)   =  975.0_r8  ;  std_atm_pres_col(6)   = 8.222E-09_r8
std_atm_hgt_col(7)   =  970.0_r8  ;  std_atm_pres_col(7)   = 8.371E-09_r8
std_atm_hgt_col(8)   =  965.0_r8  ;  std_atm_pres_col(8)   = 8.524E-09_r8
std_atm_hgt_col(9)   =  960.0_r8  ;  std_atm_pres_col(9)   = 8.680E-09_r8
std_atm_hgt_col(10)  =  955.0_r8  ;  std_atm_pres_col(10)  = 8.839E-09_r8
std_atm_hgt_col(11)  =  950.0_r8  ;  std_atm_pres_col(11)  = 9.001E-09_r8
std_atm_hgt_col(12)  =  945.0_r8  ;  std_atm_pres_col(12)  = 9.168E-09_r8
std_atm_hgt_col(13)  =  940.0_r8  ;  std_atm_pres_col(13)  = 9.338E-09_r8
std_atm_hgt_col(14)  =  935.0_r8  ;  std_atm_pres_col(14)  = 9.513E-09_r8
std_atm_hgt_col(15)  =  930.0_r8  ;  std_atm_pres_col(15)  = 9.692E-09_r8
std_atm_hgt_col(16)  =  925.0_r8  ;  std_atm_pres_col(16)  = 9.875E-09_r8
std_atm_hgt_col(17)  =  920.0_r8  ;  std_atm_pres_col(17)  = 1.006E-08_r8
std_atm_hgt_col(18)  =  915.0_r8  ;  std_atm_pres_col(18)  = 1.026E-08_r8
std_atm_hgt_col(19)  =  910.0_r8  ;  std_atm_pres_col(19)  = 1.046E-08_r8
std_atm_hgt_col(20)  =  905.0_r8  ;  std_atm_pres_col(20)  = 1.066E-08_r8
std_atm_hgt_col(21)  =  900.0_r8  ;  std_atm_pres_col(21)  = 1.087E-08_r8
std_atm_hgt_col(22)  =  895.0_r8  ;  std_atm_pres_col(22)  = 1.109E-08_r8
std_atm_hgt_col(23)  =  890.0_r8  ;  std_atm_pres_col(23)  = 1.132E-08_r8
std_atm_hgt_col(24)  =  885.0_r8  ;  std_atm_pres_col(24)  = 1.155E-08_r8
std_atm_hgt_col(25)  =  880.0_r8  ;  std_atm_pres_col(25)  = 1.179E-08_r8
std_atm_hgt_col(26)  =  875.0_r8  ;  std_atm_pres_col(26)  = 1.203E-08_r8
std_atm_hgt_col(27)  =  870.0_r8  ;  std_atm_pres_col(27)  = 1.229E-08_r8
std_atm_hgt_col(28)  =  865.0_r8  ;  std_atm_pres_col(28)  = 1.255E-08_r8
std_atm_hgt_col(29)  =  860.0_r8  ;  std_atm_pres_col(29)  = 1.283E-08_r8
std_atm_hgt_col(30)  =  855.0_r8  ;  std_atm_pres_col(30)  = 1.311E-08_r8
std_atm_hgt_col(31)  =  850.0_r8  ;  std_atm_pres_col(31)  = 1.340E-08_r8
std_atm_hgt_col(32)  =  845.0_r8  ;  std_atm_pres_col(32)  = 1.371E-08_r8
std_atm_hgt_col(33)  =  840.0_r8  ;  std_atm_pres_col(33)  = 1.402E-08_r8
std_atm_hgt_col(34)  =  835.0_r8  ;  std_atm_pres_col(34)  = 1.435E-08_r8
std_atm_hgt_col(35)  =  830.0_r8  ;  std_atm_pres_col(35)  = 1.469E-08_r8
std_atm_hgt_col(36)  =  825.0_r8  ;  std_atm_pres_col(36)  = 1.504E-08_r8
std_atm_hgt_col(37)  =  820.0_r8  ;  std_atm_pres_col(37)  = 1.541E-08_r8
std_atm_hgt_col(38)  =  815.0_r8  ;  std_atm_pres_col(38)  = 1.579E-08_r8
std_atm_hgt_col(39)  =  810.0_r8  ;  std_atm_pres_col(39)  = 1.619E-08_r8
std_atm_hgt_col(40)  =  805.0_r8  ;  std_atm_pres_col(40)  = 1.660E-08_r8
std_atm_hgt_col(41)  =  800.0_r8  ;  std_atm_pres_col(41)  = 1.704E-08_r8
std_atm_hgt_col(42)  =  795.0_r8  ;  std_atm_pres_col(42)  = 1.749E-08_r8
std_atm_hgt_col(43)  =  790.0_r8  ;  std_atm_pres_col(43)  = 1.795E-08_r8
std_atm_hgt_col(44)  =  785.0_r8  ;  std_atm_pres_col(44)  = 1.844E-08_r8
std_atm_hgt_col(45)  =  780.0_r8  ;  std_atm_pres_col(45)  = 1.896E-08_r8
std_atm_hgt_col(46)  =  775.0_r8  ;  std_atm_pres_col(46)  = 1.949E-08_r8
std_atm_hgt_col(47)  =  770.0_r8  ;  std_atm_pres_col(47)  = 2.006E-08_r8
std_atm_hgt_col(48)  =  765.0_r8  ;  std_atm_pres_col(48)  = 2.064E-08_r8
std_atm_hgt_col(49)  =  760.0_r8  ;  std_atm_pres_col(49)  = 2.126E-08_r8
std_atm_hgt_col(50)  =  755.0_r8  ;  std_atm_pres_col(50)  = 2.191E-08_r8
std_atm_hgt_col(51)  =  750.0_r8  ;  std_atm_pres_col(51)  = 2.260E-08_r8
std_atm_hgt_col(52)  =  745.0_r8  ;  std_atm_pres_col(52)  = 2.331E-08_r8
std_atm_hgt_col(53)  =  740.0_r8  ;  std_atm_pres_col(53)  = 2.407E-08_r8
std_atm_hgt_col(54)  =  735.0_r8  ;  std_atm_pres_col(54)  = 2.487E-08_r8
std_atm_hgt_col(55)  =  730.0_r8  ;  std_atm_pres_col(55)  = 2.571E-08_r8
std_atm_hgt_col(56)  =  725.0_r8  ;  std_atm_pres_col(56)  = 2.660E-08_r8
std_atm_hgt_col(57)  =  720.0_r8  ;  std_atm_pres_col(57)  = 2.755E-08_r8
std_atm_hgt_col(58)  =  715.0_r8  ;  std_atm_pres_col(58)  = 2.854E-08_r8
std_atm_hgt_col(59)  =  710.0_r8  ;  std_atm_pres_col(59)  = 2.960E-08_r8
std_atm_hgt_col(60)  =  705.0_r8  ;  std_atm_pres_col(60)  = 3.072E-08_r8
std_atm_hgt_col(61)  =  700.0_r8  ;  std_atm_pres_col(61)  = 3.191E-08_r8
std_atm_hgt_col(62)  =  695.0_r8  ;  std_atm_pres_col(62)  = 3.317E-08_r8
std_atm_hgt_col(63)  =  690.0_r8  ;  std_atm_pres_col(63)  = 3.451E-08_r8
std_atm_hgt_col(64)  =  685.0_r8  ;  std_atm_pres_col(64)  = 3.594E-08_r8
std_atm_hgt_col(65)  =  680.0_r8  ;  std_atm_pres_col(65)  = 3.746E-08_r8
std_atm_hgt_col(66)  =  675.0_r8  ;  std_atm_pres_col(66)  = 3.908E-08_r8
std_atm_hgt_col(67)  =  670.0_r8  ;  std_atm_pres_col(67)  = 4.080E-08_r8
std_atm_hgt_col(68)  =  665.0_r8  ;  std_atm_pres_col(68)  = 4.264E-08_r8
std_atm_hgt_col(69)  =  660.0_r8  ;  std_atm_pres_col(69)  = 4.459E-08_r8
std_atm_hgt_col(70)  =  655.0_r8  ;  std_atm_pres_col(70)  = 4.668E-08_r8
std_atm_hgt_col(71)  =  650.0_r8  ;  std_atm_pres_col(71)  = 4.892E-08_r8
std_atm_hgt_col(72)  =  645.0_r8  ;  std_atm_pres_col(72)  = 5.130E-08_r8
std_atm_hgt_col(73)  =  640.0_r8  ;  std_atm_pres_col(73)  = 5.385E-08_r8
std_atm_hgt_col(74)  =  635.0_r8  ;  std_atm_pres_col(74)  = 5.659E-08_r8
std_atm_hgt_col(75)  =  630.0_r8  ;  std_atm_pres_col(75)  = 5.951E-08_r8
std_atm_hgt_col(76)  =  625.0_r8  ;  std_atm_pres_col(76)  = 6.264E-08_r8
std_atm_hgt_col(77)  =  620.0_r8  ;  std_atm_pres_col(77)  = 6.600E-08_r8
std_atm_hgt_col(78)  =  615.0_r8  ;  std_atm_pres_col(78)  = 6.961E-08_r8
std_atm_hgt_col(79)  =  610.0_r8  ;  std_atm_pres_col(79)  = 7.349E-08_r8
std_atm_hgt_col(80)  =  605.0_r8  ;  std_atm_pres_col(80)  = 7.765E-08_r8
std_atm_hgt_col(81)  =  600.0_r8  ;  std_atm_pres_col(81)  = 8.213E-08_r8
std_atm_hgt_col(82)  =  595.0_r8  ;  std_atm_pres_col(82)  = 8.695E-08_r8
std_atm_hgt_col(83)  =  590.0_r8  ;  std_atm_pres_col(83)  = 9.214E-08_r8
std_atm_hgt_col(84)  =  585.0_r8  ;  std_atm_pres_col(84)  = 9.774E-08_r8
std_atm_hgt_col(85)  =  580.0_r8  ;  std_atm_pres_col(85)  = 1.038E-07_r8
std_atm_hgt_col(86)  =  575.0_r8  ;  std_atm_pres_col(86)  = 1.103E-07_r8
std_atm_hgt_col(87)  =  570.0_r8  ;  std_atm_pres_col(87)  = 1.173E-07_r8
std_atm_hgt_col(88)  =  565.0_r8  ;  std_atm_pres_col(88)  = 1.249E-07_r8
std_atm_hgt_col(89)  =  560.0_r8  ;  std_atm_pres_col(89)  = 1.330E-07_r8
std_atm_hgt_col(90)  =  555.0_r8  ;  std_atm_pres_col(90)  = 1.418E-07_r8
std_atm_hgt_col(91)  =  550.0_r8  ;  std_atm_pres_col(91)  = 1.514E-07_r8
std_atm_hgt_col(92)  =  545.0_r8  ;  std_atm_pres_col(92)  = 1.617E-07_r8
std_atm_hgt_col(93)  =  540.0_r8  ;  std_atm_pres_col(93)  = 1.728E-07_r8
std_atm_hgt_col(94)  =  535.0_r8  ;  std_atm_pres_col(94)  = 1.849E-07_r8
std_atm_hgt_col(95)  =  530.0_r8  ;  std_atm_pres_col(95)  = 1.979E-07_r8
std_atm_hgt_col(96)  =  525.0_r8  ;  std_atm_pres_col(96)  = 2.120E-07_r8
std_atm_hgt_col(97)  =  520.0_r8  ;  std_atm_pres_col(97)  = 2.273E-07_r8
std_atm_hgt_col(98)  =  515.0_r8  ;  std_atm_pres_col(98)  = 2.439E-07_r8
std_atm_hgt_col(99)  =  510.0_r8  ;  std_atm_pres_col(99)  = 2.618E-07_r8
std_atm_hgt_col(100) =  505.0_r8  ;  std_atm_pres_col(100) = 2.813E-07_r8
std_atm_hgt_col(101) =  500.0_r8  ;  std_atm_pres_col(101) = 3.024E-07_r8
std_atm_hgt_col(102) =  495.0_r8  ;  std_atm_pres_col(102) = 3.252E-07_r8
std_atm_hgt_col(103) =  490.0_r8  ;  std_atm_pres_col(103) = 3.501E-07_r8
std_atm_hgt_col(104) =  485.0_r8  ;  std_atm_pres_col(104) = 3.770E-07_r8
std_atm_hgt_col(105) =  480.0_r8  ;  std_atm_pres_col(105) = 4.063E-07_r8
std_atm_hgt_col(106) =  475.0_r8  ;  std_atm_pres_col(106) = 4.382E-07_r8
std_atm_hgt_col(107) =  470.0_r8  ;  std_atm_pres_col(107) = 4.728E-07_r8
std_atm_hgt_col(108) =  465.0_r8  ;  std_atm_pres_col(108) = 5.104E-07_r8
std_atm_hgt_col(109) =  460.0_r8  ;  std_atm_pres_col(109) = 5.514E-07_r8
std_atm_hgt_col(110) =  455.0_r8  ;  std_atm_pres_col(110) = 5.960E-07_r8
std_atm_hgt_col(111) =  450.0_r8  ;  std_atm_pres_col(111) = 6.445E-07_r8
std_atm_hgt_col(112) =  445.0_r8  ;  std_atm_pres_col(112) = 6.974E-07_r8
std_atm_hgt_col(113) =  440.0_r8  ;  std_atm_pres_col(113) = 7.550E-07_r8
std_atm_hgt_col(114) =  435.0_r8  ;  std_atm_pres_col(114) = 8.179E-07_r8
std_atm_hgt_col(115) =  430.0_r8  ;  std_atm_pres_col(115) = 8.864E-07_r8
std_atm_hgt_col(116) =  425.0_r8  ;  std_atm_pres_col(116) = 9.612E-07_r8
std_atm_hgt_col(117) =  420.0_r8  ;  std_atm_pres_col(117) = 1.043E-06_r8
std_atm_hgt_col(118) =  415.0_r8  ;  std_atm_pres_col(118) = 1.132E-06_r8
std_atm_hgt_col(119) =  410.0_r8  ;  std_atm_pres_col(119) = 1.229E-06_r8
std_atm_hgt_col(120) =  405.0_r8  ;  std_atm_pres_col(120) = 1.336E-06_r8
std_atm_hgt_col(121) =  400.0_r8  ;  std_atm_pres_col(121) = 1.452E-06_r8
std_atm_hgt_col(122) =  395.0_r8  ;  std_atm_pres_col(122) = 1.579E-06_r8
std_atm_hgt_col(123) =  390.0_r8  ;  std_atm_pres_col(123) = 1.718E-06_r8
std_atm_hgt_col(124) =  385.0_r8  ;  std_atm_pres_col(124) = 1.870E-06_r8
std_atm_hgt_col(125) =  380.0_r8  ;  std_atm_pres_col(125) = 2.037E-06_r8
std_atm_hgt_col(126) =  375.0_r8  ;  std_atm_pres_col(126) = 2.220E-06_r8
std_atm_hgt_col(127) =  370.0_r8  ;  std_atm_pres_col(127) = 2.421E-06_r8
std_atm_hgt_col(128) =  365.0_r8  ;  std_atm_pres_col(128) = 2.641E-06_r8
std_atm_hgt_col(129) =  360.0_r8  ;  std_atm_pres_col(129) = 2.884E-06_r8
std_atm_hgt_col(130) =  355.0_r8  ;  std_atm_pres_col(130) = 3.151E-06_r8
std_atm_hgt_col(131) =  350.0_r8  ;  std_atm_pres_col(131) = 3.445E-06_r8
std_atm_hgt_col(132) =  345.0_r8  ;  std_atm_pres_col(132) = 3.769E-06_r8
std_atm_hgt_col(133) =  340.0_r8  ;  std_atm_pres_col(133) = 4.126E-06_r8
std_atm_hgt_col(134) =  335.0_r8  ;  std_atm_pres_col(134) = 4.521E-06_r8
std_atm_hgt_col(135) =  330.0_r8  ;  std_atm_pres_col(135) = 4.957E-06_r8
std_atm_hgt_col(136) =  325.0_r8  ;  std_atm_pres_col(136) = 5.440E-06_r8
std_atm_hgt_col(137) =  320.0_r8  ;  std_atm_pres_col(137) = 5.975E-06_r8
std_atm_hgt_col(138) =  315.0_r8  ;  std_atm_pres_col(138) = 6.568E-06_r8
std_atm_hgt_col(139) =  310.0_r8  ;  std_atm_pres_col(139) = 7.226E-06_r8
std_atm_hgt_col(140) =  305.0_r8  ;  std_atm_pres_col(140) = 7.957E-06_r8
std_atm_hgt_col(141) =  300.0_r8  ;  std_atm_pres_col(141) = 8.770E-06_r8
std_atm_hgt_col(142) =  295.0_r8  ;  std_atm_pres_col(142) = 9.676E-06_r8
std_atm_hgt_col(143) =  290.0_r8  ;  std_atm_pres_col(143) = 1.069E-05_r8
std_atm_hgt_col(144) =  285.0_r8  ;  std_atm_pres_col(144) = 1.181E-05_r8
std_atm_hgt_col(145) =  280.0_r8  ;  std_atm_pres_col(145) = 1.308E-05_r8
std_atm_hgt_col(146) =  275.0_r8  ;  std_atm_pres_col(146) = 1.449E-05_r8
std_atm_hgt_col(147) =  270.0_r8  ;  std_atm_pres_col(147) = 1.608E-05_r8
std_atm_hgt_col(148) =  265.0_r8  ;  std_atm_pres_col(148) = 1.787E-05_r8
std_atm_hgt_col(149) =  260.0_r8  ;  std_atm_pres_col(149) = 1.989E-05_r8
std_atm_hgt_col(150) =  255.0_r8  ;  std_atm_pres_col(150) = 2.218E-05_r8
std_atm_hgt_col(151) =  250.0_r8  ;  std_atm_pres_col(151) = 2.476E-05_r8
std_atm_hgt_col(152) =  245.0_r8  ;  std_atm_pres_col(152) = 2.770E-05_r8
std_atm_hgt_col(153) =  240.0_r8  ;  std_atm_pres_col(153) = 3.105E-05_r8
std_atm_hgt_col(154) =  235.0_r8  ;  std_atm_pres_col(154) = 3.488E-05_r8
std_atm_hgt_col(155) =  230.0_r8  ;  std_atm_pres_col(155) = 3.927E-05_r8
std_atm_hgt_col(156) =  225.0_r8  ;  std_atm_pres_col(156) = 4.432E-05_r8
std_atm_hgt_col(157) =  220.0_r8  ;  std_atm_pres_col(157) = 5.015E-05_r8
std_atm_hgt_col(158) =  215.0_r8  ;  std_atm_pres_col(158) = 5.690E-05_r8
std_atm_hgt_col(159) =  210.0_r8  ;  std_atm_pres_col(159) = 6.476E-05_r8
std_atm_hgt_col(160) =  205.0_r8  ;  std_atm_pres_col(160) = 7.394E-05_r8
std_atm_hgt_col(161) =  200.0_r8  ;  std_atm_pres_col(161) = 8.474E-05_r8
std_atm_hgt_col(162) =  195.0_r8  ;  std_atm_pres_col(162) = 9.749E-05_r8
std_atm_hgt_col(163) =  190.0_r8  ;  std_atm_pres_col(163) = 1.127E-04_r8
std_atm_hgt_col(164) =  185.0_r8  ;  std_atm_pres_col(164) = 1.308E-04_r8
std_atm_hgt_col(165) =  180.0_r8  ;  std_atm_pres_col(165) = 1.527E-04_r8
std_atm_hgt_col(166) =  175.0_r8  ;  std_atm_pres_col(166) = 1.794E-04_r8
std_atm_hgt_col(167) =  170.0_r8  ;  std_atm_pres_col(167) = 2.121E-04_r8
std_atm_hgt_col(168) =  165.0_r8  ;  std_atm_pres_col(168) = 2.528E-04_r8
std_atm_hgt_col(169) =  160.0_r8  ;  std_atm_pres_col(169) = 3.039E-04_r8
std_atm_hgt_col(170) =  155.0_r8  ;  std_atm_pres_col(170) = 3.693E-04_r8
std_atm_hgt_col(171) =  150.0_r8  ;  std_atm_pres_col(171) = 4.542E-04_r8
std_atm_hgt_col(172) =  145.0_r8  ;  std_atm_pres_col(172) = 5.669E-04_r8
std_atm_hgt_col(173) =  140.0_r8  ;  std_atm_pres_col(173) = 7.203E-04_r8
std_atm_hgt_col(174) =  135.0_r8  ;  std_atm_pres_col(174) = 9.357E-04_r8
std_atm_hgt_col(175) =  130.0_r8  ;  std_atm_pres_col(175) = 1.250E-03_r8
std_atm_hgt_col(176) =  125.0_r8  ;  std_atm_pres_col(176) = 1.736E-03_r8
std_atm_hgt_col(177) =  120.0_r8  ;  std_atm_pres_col(177) = 2.537E-03_r8
std_atm_hgt_col(178) =  115.0_r8  ;  std_atm_pres_col(178) = 4.004E-03_r8
std_atm_hgt_col(179) =  110.0_r8  ;  std_atm_pres_col(179) = 7.149E-03_r8
std_atm_hgt_col(180) =  105.0_r8  ;  std_atm_pres_col(180) = 1.442E-02_r8
std_atm_hgt_col(181) =  100.0_r8  ;  std_atm_pres_col(181) = 3.201E-02_r8
std_atm_hgt_col(182) =   95.0_r8  ;  std_atm_pres_col(182) = 7.577E-02_r8
std_atm_hgt_col(183) =   90.0_r8  ;  std_atm_pres_col(183) = 1.844E-01_r8
std_atm_hgt_col(184) =   85.0_r8  ;  std_atm_pres_col(184) = 4.457E-01_r8
std_atm_hgt_col(185) =   80.0_r8  ;  std_atm_pres_col(185) = 1.052E+00_r8
std_atm_hgt_col(186) =   75.0_r8  ;  std_atm_pres_col(186) = 2.388E+00_r8
std_atm_hgt_col(187) =   70.0_r8  ;  std_atm_pres_col(187) = 5.221E+00_r8
std_atm_hgt_col(188) =   65.0_r8  ;  std_atm_pres_col(188) = 1.093E+01_r8
std_atm_hgt_col(189) =   60.0_r8  ;  std_atm_pres_col(189) = 2.196E+01_r8
std_atm_hgt_col(190) =   55.0_r8  ;  std_atm_pres_col(190) = 4.253E+01_r8
std_atm_hgt_col(191) =   50.0_r8  ;  std_atm_pres_col(191) = 7.978E+01_r8
std_atm_hgt_col(192) =   45.0_r8  ;  std_atm_pres_col(192) = 1.491E+02_r8
std_atm_hgt_col(193) =   40.0_r8  ;  std_atm_pres_col(193) = 2.871E+02_r8
std_atm_hgt_col(194) =   35.0_r8  ;  std_atm_pres_col(194) = 5.746E+02_r8
std_atm_hgt_col(195) =   30.0_r8  ;  std_atm_pres_col(195) = 1.197E+03_r8
std_atm_hgt_col(196) =   25.0_r8  ;  std_atm_pres_col(196) = 2.549E+03_r8
std_atm_hgt_col(197) =   20.0_r8  ;  std_atm_pres_col(197) = 5.529E+03_r8
std_atm_hgt_col(198) =   15.0_r8  ;  std_atm_pres_col(198) = 1.211E+04_r8
std_atm_hgt_col(199) =   10.0_r8  ;  std_atm_pres_col(199) = 2.650E+04_r8
std_atm_hgt_col(200) =    5.0_r8  ;  std_atm_pres_col(200) = 5.405E+04_r8
std_atm_hgt_col(201) =    0.0_r8  ;  std_atm_pres_col(201) = 1.013E+05_r8

! convert km to m
std_atm_hgt_col(:) = std_atm_hgt_col(:) * 1000.0_r8

end subroutine load_high_top_table
	
!===================================================================
! End of model_mod
!===================================================================

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
