! DART software - Copyright UCAR. This open source software is provided
! by ucar, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/dares/dart/dart_download
!
!----------------------------------------------------------------
!>
!> this is the interface between the cam-fv atmosphere model and dart.
!> the required public interfaces and arguments cannot be changed.
!>
!----------------------------------------------------------------

module model_mod

use             types_mod,  only : MISSING_R8, MISSING_I, i8, r8, vtablenamelength, &
                                   gravity, DEG2RAD
use      time_manager_mod,  only : set_calendar_type
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

! examples of additional quantities that cam-chem might need defined from the obs_kind_mod
!                                   ! GASES
!                                   QTY_CO, QTY_SFCO, QTY_SFCO01, QTY_SFCO02, QTY_SFCO03, &
!                                   QTY_O3, QTY_OH, QTY_NO, QTY_NO2, QTY_NO3, QTY_CH2O, &
!                                  ! AEROSOLS
!                                   QTY_AOD, QTY_NUM_A1, QTY_NUM_A2, QTY_NUM_A3, QTY_NUM_A4, & ! AOD and Numbers
!                                   QTY_SFNUM_A1, QTY_SFNUM_A2, QTY_SFNUM_A3, QTY_SFNUM_A4, & ! SF / Numbers
!                                   QTY_POM_A1, QTY_POM_A4, QTY_BC_A1, QTY_BC_A4, &
!                                   QTY_SFPOM_A4, QTY_SFBC_A4, & ! Carbon
!                                   QTY_SO4_A1, QTY_SO4_A2, QTY_SO4_A3, QTY_SFSO4_A1, QTY_SFSO4_A2, & ! Sulfates
!                                   QTY_DST_A1, QTY_DST_A2, QTY_DST_A3, QTY_NCL_A1, QTY_NCL_A2, QTY_NCL_A3, &
!                                   QTY_SOA1_A1, QTY_SOA1_A2, QTY_SOA2_A1, QTY_SOA2_A2, QTY_SOA3_A1, QTY_SOA3_A2, & ! SOA
!                                   QTY_SOA4_A1, QTY_SOA4_A2, QTY_SOA5_A1, QTY_SOA5_A2, & ! SOA

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

use    cam_common_code_mod, only : above_ramp_start, are_damping, build_cam_pressure_columns, build_heights, &
                                   cam_grid, cdebug_level, check_good_levels, cno_normalization_of_scale_heights, &
                                   pert_model_copies, cuse_log_vertical_scale, discarding_high_obs, &
                                   free_cam_grid, free_std_atm_tables, generic_height_to_pressure, &
                                   gph2gmh, grid_data, height_to_level, init_damping_ramp_info, &
                                   init_discard_high_obs, init_globals, init_sign_of_vert_units, &
                                   is_surface_field, obs_too_high, ok_to_interpolate, pressure_to_level, ramp_end, &
                                   read_model_time, ref_model_top_pressure, ref_nlevels, scale_height, &
                                   set_vert_localization, vert_interp, vertical_localization_type, write_model_time

use cam_common_code_mod, only : nc_write_model_atts, grid_data, read_grid_info, &
                                set_cam_variable_info, MAX_STATE_VARIABLES, &
                                num_state_table_columns, MAX_PERT, &
                                shortest_time_between_assimilations, domain_id, &
                                cuse_log_vertical_scale, &
                                cno_normalization_of_scale_heights, &
                                cdebug_level, &
                                ccustom_routine_to_generate_ensemble, &
                                cfields_to_perturb, &
                                cperturbation_amplitude, &
                                cassimilation_period_days, &
                                cassimilation_period_seconds, &
                                csuppress_grid_info_in_output, &
                                common_initialized

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

! public routine for utility program
public :: obs_too_high

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'cam-fv/model_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''


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

!> build a pressure/height conversion column based on a
!> standard atmosphere.  this can only be used when we
!> don't have a real ensemble to use, or we don't care
!> about absolute accuracy.

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

character(len=*), parameter :: routine = 'static_init_model'

if ( module_initialized ) return

! Record version info
call register_module(source, revision, revdate)

module_initialized = .true.
common_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Set values from namelist in cam_common_code_mod
cuse_log_vertical_scale = use_log_vertical_scale
cno_normalization_of_scale_heights = no_normalization_of_scale_heights
cdebug_level = debug_level
ccustom_routine_to_generate_ensemble = custom_routine_to_generate_ensemble
ccustom_routine_to_generate_ensemble = custom_routine_to_generate_ensemble
cperturbation_amplitude = perturbation_amplitude
cassimilation_period_days = assimilation_period_days
cassimilation_period_seconds = assimilation_period_seconds
csuppress_grid_info_in_output = suppress_grid_info_in_output


call set_calendar_type('GREGORIAN')

call read_grid_info(cam_template_filename)
! This non-state variable is used to compute surface elevation.
call read_cam_phis_array(cam_phis_filename)
call setup_interpolation() !grid is global

! initialize global values that are used frequently
call init_globals()

! read the namelist &model_nml :: state_variables
! to set up what will be read into the cam state vector
call set_cam_variable_info(cam_template_filename, state_variables)

call fill_cam_stagger_info(grid_stagger)

if (debug_level > 100) call state_structure_info(domain_id)


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
      call init_damping_ramp_info(model_damping_ends_at_level)
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
   call init_discard_high_obs(no_obs_assim_above_level)
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

character(len=*), parameter :: routine = 'get_location_from_index'
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
else if (nd == 2) then
   ! add any 2d surface fields to this function
   if (is_surface_field(q)) then
      use_vert_type = VERTISSURFACE
      use_vert_val  = MISSING_R8  
      ! setting the vertical value to missing matches what the previous
      ! version of this code did.  other models choose to set the vertical
      ! value to the model surface elevation at this location:
      !   use_vert_val  = phis(lon_index, lat_index) / gravity
   else
      ! any 2d field not listed as a surface field (in is_surface_field() function) 
      ! is assumed to be an integrated quantity with a vert type of VERTISUNDEF.
      use_vert_type = VERTISUNDEF
      use_vert_val  = MISSING_R8
   endif
else
   write(string1, *) 'state vector field not 2D or 3D and no code to handle other dimensionity'
   write(string2, *) 'dimensionality = ', nd, ' quantity type = ', trim(get_name_for_quantity(q))
   call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2)
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
real(r8) :: quad_vals(ens_size,4)
type(quad_interp_handle) :: interp_handle   ! should this be a pointer?? 
                                            ! is it replicating the internal arrays on assignment?

integer :: imem

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

do imem=1,ens_size
   ! do the horizontal interpolation for each ensemble member.
   ! call one at a time to avoid creating temporary arrays
   call quad_lon_lat_evaluate(interp_handle, lon_fract, lat_fract, &
                              quad_vals(imem,:), interp_vals(imem), status_array(imem))
end do

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
real(r8),           intent(out) :: quad_vals(ens_size,4) !< array of interpolated values
integer,            intent(out) :: my_status(ens_size)

integer  :: icorner, numdims
integer  :: level_one_array(ens_size)
integer  :: four_levs1(ens_size, 4), four_levs2(ens_size, 4)
real(r8) :: four_vert_fracts(ens_size, 4)

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
                                four_levs1(:,icorner), four_levs2(:,icorner), & 
                                four_vert_fracts(:,icorner), my_status)
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
                                    level_one_array, varid, quad_vals(:,icorner),my_status)
         if (any(my_status /= 0)) return

      enddo

   else ! special 2d case
      do icorner=1, 4
         call get_quad_values(ens_size, four_lons(icorner), four_lats(icorner), &
                               obs_qty, obs_qty, quad_vals(:,icorner))
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
integer,             intent(in) :: four_levs1(ens_size,4), four_levs2(ens_size,4)
real(r8),            intent(in) :: four_vert_fracts(ens_size,4)
integer,             intent(in) :: varid
real(r8),           intent(out) :: quad_vals(ens_size,4) !< array of interpolated values
integer,            intent(out) :: my_status(ens_size)

integer  :: icorner
real(r8) :: vals1(ens_size), vals2(ens_size)

character(len=*), parameter :: routine = 'get_four_state_values:'

do icorner=1, 4
   call get_values_from_varid(state_handle,  ens_size, &
                              four_lons(icorner), four_lats(icorner), &
                              four_levs1(:,icorner), varid, vals1, &
                              my_status)

   if (any(my_status /= 0)) then
      my_status(:) = 16   ! cannot retrieve vals1 values
      return
   endif

   call get_values_from_varid(state_handle,  ens_size, &
                              four_lons(icorner), four_lats(icorner), &
                              four_levs2(:,icorner), varid, vals2, my_status)
   if (any(my_status /= 0)) then
      my_status(:) = 17   ! cannot retrieve top values
      return
   endif

   call vert_interp(ens_size, vals1, vals2, four_vert_fracts(:,icorner), & 
                    quad_vals(:,icorner))

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
integer,             intent(in) :: four_levs1(ens_size,4), four_levs2(ens_size,4)
real(r8),            intent(in) :: four_vert_fracts(ens_size,4)
integer,             intent(in) :: obs_qty
real(r8),           intent(out) :: quad_vals(ens_size,4) !< array of interpolated values
integer,            intent(out) :: my_status(ens_size)

integer  :: icorner
real(r8) :: vals1(ens_size), vals2(ens_size)

character(len=*), parameter :: routine = 'get_four_nonstate_values:'

do icorner=1, 4
   call get_values_from_nonstate_fields(state_handle,  ens_size, &
                              four_lons(icorner), four_lats(icorner), &
                              four_levs1(:,icorner), obs_qty, vals1, my_status)
   if (any(my_status /= 0)) then
      my_status(:) = 16   ! cannot retrieve vals1 values
      return
   endif

   call get_values_from_nonstate_fields(state_handle,  ens_size, &
                              four_lons(icorner), four_lats(icorner), &
                              four_levs2(:,icorner), obs_qty, vals2, my_status)
   if (any(my_status /= 0)) then
      my_status(:) = 17   ! cannot retrieve top values
      return
   endif

   call vert_interp(ens_size, vals1, vals2, four_vert_fracts(:,icorner), &
                    quad_vals(:,icorner))

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
real(r8) :: pressure(nlevels, ens_size)
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

! Build the pressure columns for the entire ensemble
call build_cam_pressure_columns(ens_size, surface_pressure, nlevels, pressure)

if (use_variable_mean_mass) then
   call compute_mean_mass(ens_handle, ens_size, lon_index, lat_index, nlevels, qty, mbar, status1)

   if (status1 /= 0) then
      my_status = status1
      return
   endif

   ! compute the height columns for each ensemble member - passing mbar() array in.
   do imember = 1, ens_size
      call build_heights(nlevels, surface_pressure(imember), surface_elevation(1), &
                         pressure(:, imember), tv(:, imember), height_array(:, imember), mbar=mbar(:, imember))
   enddo

else
   ! compute the height columns for each ensemble member - no mbar() argument here.
   ! (you cannot just pass 1.0 in for the mbar() array; it uses a different gas constant
   ! in the variable mean mass case.)
   do imember = 1, ens_size
      call build_heights(nlevels, surface_pressure(imember), surface_elevation(1), &
                         pressure(:, imember), tv(:, imember), height_array(:, imember))
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
!> Does any shutdown and clean-up needed for model.
!>

subroutine end_model()

! deallocate arrays from grid and anything else

call free_cam_grid(grid_data)

deallocate(phis)

call free_std_atm_tables()

call finalize_quad_interp(interp_nonstaggered)
call finalize_quad_interp(interp_u_staggered)
call finalize_quad_interp(interp_v_staggered)

if (using_chemistry) call finalize_chem_tables()

end subroutine end_model

!-----------------------------------------------------------------------
! The remaining (private) interfaces come last.
! None of the private interfaces need to call static_init_model()
!-----------------------------------------------------------------------

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
!>
!> 
!>   

subroutine setup_interpolation()

!>@todo FIXME the cam fv grid is really evenly spaced in lat and lon,
!>even though they provide full lon() and lat() arrays.  providing the deltas
!>between each pair would be slightly faster inside the interp code.

!print *, 'setting up interpolation: lon/lat sizes = ', grid%lon%nsize, grid%lat%nsize,  &
!                                                       grid%slon%nsize, grid%slat%nsize

! mass points at cell centers
call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, grid_data%lon%nsize, grid_data%lat%nsize, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_nonstaggered)
call set_quad_coords(interp_nonstaggered, grid_data%lon%vals, grid_data%lat%vals)

! U stagger
call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, grid_data%lon%nsize, grid_data%slat%nsize, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_u_staggered)
call set_quad_coords(interp_u_staggered, grid_data%lon%vals, grid_data%slat%vals)

! V stagger
call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, grid_data%slon%nsize, grid_data%lat%nsize, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_v_staggered)
call set_quad_coords(interp_v_staggered, grid_data%slon%vals, grid_data%lat%vals)

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

! if this location is on the surface, use the surface pressure field
! in the computations below.  otherwise use the 3d pressure field.
if (query_location(location) == VERTISSURFACE) then
   ptype = QTY_SURFACE_PRESSURE
else
   ptype = QTY_PRESSURE
endif


! there are 4 cases here.

if (no_normalization_of_scale_heights) then

   ! take log of pressure, either surface pressure or regular pressure
   
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

!===================================================================
! End of model_mod
!===================================================================

end module model_mod

