! DART software - Copyright UCAR. This open source software is provided
! by ucar, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/dares/dart/dart_download
!
! $Id$
!----------------------------------------------------------------
!>
!> this is the interface between the cam-fv atmosphere model and dart.
!> the required public interfaces and arguments cannot be changed.


! This is a prototype version of CAM-SE with Manhattan using common code
!>
!----------------------------------------------------------------

module model_mod

use             types_mod,  only : MISSING_R8, MISSING_I, i8, r8, vtablenamelength, &
                                   gravity, DEG2RAD, PI, earth_radius
use      time_manager_mod,  only : set_time, time_type, set_date, &
                                   set_calendar_type, get_date
use          location_mod,  only : location_type, set_vertical, set_location, &
                                   get_location, write_location, is_vertical, &
                                   VERTISUNDEF, VERTISSURFACE, VERTISLEVEL, &
                                   VERTISPRESSURE, VERTISHEIGHT, &
                                   VERTISSCALEHEIGHT, query_location, &
                                   set_vertical_localization_coord, get_dist, &
                                   loc_get_close_obs => get_close_obs, &
                                   get_close, &
                                   loc_get_close_state => get_close_state, &
                                   vertical_localization_on, get_close_type, get_maxdist, &
                                   get_close_init

!SENote using 3d space location for now
! NEED to understand how this relates to the standard threed_cartesian which is in the same directory
! IN Addiion, the xyz_location_mod in the threed_cartesian directory in Manhattan is NOT consistent with
! the one from CLASSIC. For now, I've added in the classic module in that directory but eventually this will
! need reconciled. 
use xyz_location_mod,       only : xyz_location_type, xyz_get_close_maxdist_init,          &
                                   xyz_get_close_type, xyz_set_location, xyz_get_location, &
                                   xyz_get_close_obs_init, xyz_get_close_obs_destroy,      &
                                   xyz_find_nearest



use         utilities_mod,  only : find_namelist_in_file, check_namelist_read, &
                                   string_to_logical, string_to_real,& 
                                   nmlfileunit, do_nml_file, do_nml_term, &
                                   register_module, error_handler, &
                                   file_exist, to_upper, E_ERR, E_MSG, E_WARN, array_dump, &
                                   find_enclosing_indices, nc_check
use          obs_kind_mod,  only : QTY_SURFACE_ELEVATION, QTY_PRESSURE, &
                                   QTY_GEOMETRIC_HEIGHT, QTY_VERTLEVEL, &
                                   QTY_SURFACE_PRESSURE, &
                                   QTY_TEMPERATURE, QTY_SPECIFIC_HUMIDITY, &
                                   QTY_MOLEC_OXYGEN_MIXING_RATIO, &
                                   QTY_ION_O_MIXING_RATIO, QTY_ATOMIC_H_MIXING_RATIO, &
                                   QTY_ATOMIC_OXYGEN_MIXING_RATIO, QTY_NITROGEN, &
                                   ! SENote: Added in for tests
                                   QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT, QTY_CLOUD_LIQUID_WATER, QTY_CLOUD_ICE, &
         
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
                                   get_model_variable_indices, state_structure_info, get_short_name, &
                                   get_long_name, get_dim_lengths, get_variable_name
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
                                   nc_close_file, nc_variable_exists, nc_get_global_attribute
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

use    cam_common_code_mod, only : scale_height, HIGH_TOP_TABLE, LOW_TOP_TABLE, std_atm_table_len, &
                                   std_atm_hgt_col, std_atm_pres_col, load_low_top_table, &
                                   load_high_top_table, store_std_atm_tables, free_std_atm_tables, &
                                   high_top_threshold, is_surface_field, init_globals, ref_model_top_pressure, &
                                   ref_surface_pressure, ref_nlevels, cam_1d_array, cam_grid, grid_data, &
                                   are_damping, ramp_end, discarding_high_obs, no_assim_above_height, &
                                   no_assim_above_level, no_assim_above_pressure, vertical_localization_type, &
                                   above_ramp_start, v_above, v_down, v_difference, higher_is_smaller, &
                                   pressure_to_level, cuse_log_vertical_scale, generic_pressure_to_height, &
                                   single_pressure_value, convert_vertical_level_generic, &
                                   cno_normalization_of_scale_heights, init_sign_of_vert_units, &
                                   init_damping_ramp_info, init_discard_high_obs, single_pressure_column, &
                                   build_cam_pressure_columns, height_to_level, check_good_levels, &
                                   generic_cam_pressure_to_cam_level, compute_surface_gravity, gph2gmh, &
                                   build_heights, set_vert_localization, ok_to_interpolate, obs_too_high, &
                                   cdebug_level, get_cam_grid, free_cam_1d_array, free_cam_grid

!SENote the routines to read in the grid geometry use the old netcdf calls
! Tim's new calls are much better and should switch ASAP
! Note also the failure to have use, only here which is annoying
use netcdf
use typeSizes

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
character(len=256), parameter :: source   = "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! maximum number of fields you can list to be perturbed
! to generate an ensemble if starting from a single state.
integer, parameter :: MAX_PERT = 100

! model_nml namelist variables and default values
character(len=256) :: cam_template_filename           = 'caminput.nc'
character(len=256) :: cam_phis_filename               = 'cam_phis.nc'

!SENote added in namelist strings to identify the CS grid mapping files
! NOTE THAT THES lenghts may need to be 256 when we move to new netcdf code
character(len=256) :: homme_map_file                  = 'SEMapping.nc'            ! Corners of each cubed sphere cell.
character(len=256) :: cs_grid_file                    = 'SEMapping_cs_grid.nc'    ! Relationships among corners/nodes.

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
   homme_map_file,                      &
   cs_grid_file,                        &
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

! SENote: the stagger info in the next block will not be required
integer, parameter :: STAGGER_NONE = -1
integer, parameter :: STAGGER_U    =  1
integer, parameter :: STAGGER_V    =  2
integer, parameter :: STAGGER_W    =  3 
integer, parameter :: STAGGER_UV   =  4

type cam_stagger
   integer, allocatable :: qty_stagger(:)
end type

type(cam_stagger) :: grid_stagger
! SENote; end of stagger block that will be deleted

! Surface potential; used for calculation of geometric heights.
! SENote: phis will only have one dimension meaningful dimension. Could it be put in state structure?
! SENote: right now every process has their own complete copy of this
! SENote: Initially appears that SE does not use phis
real(r8), allocatable :: phis(:, :)

!> build a pressure/height conversion column based on a
!> standard atmosphere.  this can only be used when we
!> don't have a real ensemble to use, or we don't care
!> about absolute accuracy.

!SENote This type is not needed for SE
! Horizontal interpolation code.  Need a handle for nonstaggered, U and V.
type(quad_interp_handle) :: interp_nonstaggered, &
                            interp_u_staggered, &
                            interp_v_staggered


! SENote: A veriety of module storage data structures for geometry of grid
! Verify that all of these are still being used
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CS Variables holding relationships among cubed sphere nodes.
! Read in and/or set in static_init_model (and routines it calls) at the beginning of an assimilation.
! Used by model_interpolate.
logical :: l_rectang = .true.        ! Flag to tell whether grid is logically rectangular
                                     ! (Eul, FV) or not (cubed-sphere, ...?)
                                     ! Will be set .false. if cam_template_filename has dimension 'ncol'.
logical :: l_refined = .false.       ! Flag to tell whether grid is a refined mesh or not.

! ne = global metadata from cam_template_filename, giving the number of 'elements' per face edge of the 'cube'.
! np = # of nodes/edge of each element.  Edges are shared with adjacent elements.
! FIXME? put these in a derived type to prevent accidentally using them as local variables.
integer :: ne, np

! Number of columns, or nodes, in the cubed-sphere grid.
integer :: ncol

! The nominal resolution is (30 degrees/ne), assuming np = 4 (3x3 cells per element).
real(r8) :: coarse_grid

! Dimensions of array 'corners', from homme_map_file.
integer :: ncorners, ncenters

! Maximum number of neighbors a node can have (6 in refined, 4 otherwise)
! Get from namelist, to reduce file & array sizes?
! Or derive from l_refined, after ne is read from caminput.nc?
integer, parameter :: max_neighbors = 6

! Array from homme_map_file.
integer,  allocatable :: corners(:,:)    ! The 4 corners (nodes) of each cell, from HommeMapping.nc

! 5 arrays from cs_grid_file
integer,  allocatable :: num_nghbrs(:)   ! Number of neighbors of each node/column in the cubed sphere grid.
integer,  allocatable :: centers(:,:)    ! The names of the cells that use each node as a corner.
real(r8), allocatable :: a(:,:,:)        ! Coefficients of mapping from planar to unit square space for 'x'
real(r8), allocatable :: b(:,:,:)        ! Coefficients of mapping from planar to unit square space for 'y'
real(r8), allocatable :: x_ax_bearings(:,:)  ! The directions from each node to its neighbors,
                                         ! measured from the vector pointing north.  (-PI <= bearing <= PI)

! Locations of cubed sphere nodes, in DART's location_type format.
type(location_type), allocatable :: cs_locs(:)

! Location of cubed sphere nodes, in cartesian coordinates
type(xyz_location_type), allocatable :: cs_locs_xyz(:)
type(xyz_get_close_type)             :: cs_gc_xyz

! Structure containing grid point locations, etc.,
! defined in static_init_mod after reading in CS lons, lats, and levels,
! needed in model_interpolate:interp_cubed_sphere.
type(get_close_type) :: cs_gc

! Array of KINDs of cubed sphere grid points.
! As of 2014-3-28 this is only used by location_mod, which doesn't actually use it.
integer, allocatable :: cs_kinds(:)


! Other useful 1D grid arrays (for cubed sphere)
real(r8), allocatable :: lon_rad(:), lat_rad(:)   ! longitude and latitude in radians, used by bearings()

! This integer is a reminder that some shared calls take 3 dimensions, but SE has only 2: Value is irrelevant
integer :: no_third_dimension = -99

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!SENote: This variable gives extra output for the locations. More global way to do this?:w
! set to .true. to get more details about the state vector and the
! CAM fields and sizes in the init code.
logical :: print_details = .true.




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

!SENote
type(location_type) :: test_loc
integer(i8) :: test_index_in
integer :: test_var_type, i, nc_file_ID
real(r8) :: test_loc_vals(3)

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

! Set values from namelist in cam_common_code_mod
cuse_log_vertical_scale = use_log_vertical_scale
cno_normalization_of_scale_heights = no_normalization_of_scale_heights
cdebug_level = debug_level

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


!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! SENote: This section reads and/or builds the tables needed for model_interpolate for SE
! Read in or create a file containing the relationships among cubed sphere nodes,
! such as neighbors, centers, and bearings, which will be used to identify the cell
! which contains an observation.
! Fields will be stored in global storage.
! Write the cubed sphere grid arrays to a new NetCDF file.

if (file_exist(cs_grid_file)) then
   call nc_read_cs_grid_file()
elseif (file_exist(homme_map_file)) then
   call create_cs_grid_arrays()
   if (my_task_id() == 0) call nc_write_cs_grid_file( cs_grid_file, homme_map_file )
else
   write(string1, *)'No cs_grid_file "',trim(cs_grid_file), &
                 '" nor homme_map_file "',trim(homme_map_file),'"'
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif


!SENote: Following indented block is also used for search for close corners
! Need to make sure that all of it is really needed
   ! Read some attributes from the cubed sphere model_config_file.
   ! ne is the number of elements/cube edge.  Usually 0 for refined grids.
   ! np is the number of nodes/element edge (shared with adjacent element.
   nc_file_ID = nc_open_file_readonly(cam_template_filename, 'Reading ne and np from cam template file')
   call nc_get_global_attribute(nc_file_ID, 'ne', ne, 'Reading ne from cam template file', cam_template_filename)
   call nc_get_global_attribute(nc_file_ID, 'np', np, 'Reading np from cam template file', cam_template_filename)
   call nc_close_file(nc_file_ID, 'Reading ne and np from cam template file', cam_template_filename)

   ! Calculate the nominal resolution of the (coarse) grid,
   ! for use by model_interpolate's call to get_close_obs.
   if (ne == 0) then
      ! Refined mesh; assume the coarsest grid is the default '1-degree'.
      ! Need factor of 1.5 to make sure that there are at least 2 nodes 'close' to any location.
      ! There seems to be a tricky interplay between the lon-lat boxes used in the quick search
      ! for potentially close nodes, and the cubed sphere grid, so that a coarse_grid of only
      ! slightly more than 1.0 degrees can yield 0 close nodes.
      coarse_grid = 1.2_r8 * DEG2RAD
      l_refined = .true.
   else
      ! Standard cubed sphere; there are 3x num_elements/face_edge x 4 nodes
      ! around the equator.  ne = 30 -> 3x4x30 = 360 nodes -> '1-degree'
      ! Yielded a location with only 1 close ob, but need 2.
      ! coarse_grid = (30.01_r8/ne) * DEG2RAD
      coarse_grid = 1.2_r8*(30.0_r8/ne) * DEG2RAD
   endif
   if (print_details) then
      write(string1, *) 'Cubed sphere coarse_grid resolution (rad) used in cs_gc definition = ',&
                      coarse_grid,' because ne = ',ne
      call error_handler(E_MSG, 'static_init_model', string1,source,revision,revdate)
   endif

   ! Fill cs_gc for use by model_mod.  Inputs and outputs are in global storage.
   ! ncol is already defined before this call.
   call fill_gc()

   ! Fill arrays that are useful for bearings and distances.
   allocate(lon_rad(ncol), lat_rad(ncol))
   do i=1,ncol
      lon_rad(i) = grid_data%lon%vals(i)*DEG2RAD
      lat_rad(i) = grid_data%lat%vals(i)*DEG2RAD
   enddo

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! height
! Get dimensions and surface geopotential from a new netcdf file and test for consistency.
! Open file and read PHIS from it.
! Allocate global variables which will be used in vertical interpolations
! Check for pressures on vertically staggered grid, as well as standard grid.

!SENote call read_cam_2Dreal(cam_phis, 'PHIS')

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!SENote; Looking at metadata
!write(*, *) 'end of static init model'
!write(*, *) 'model size ', get_model_size()
!do i = 1, 9500000, 48602
   !test_index_in = i
   !call get_state_meta_data(test_index_in, test_loc, test_var_type)
   !test_loc_vals = get_location(test_loc)
   !write(*, *) 'variable ', i, ' location is ', test_loc_vals(1:3), test_var_type
   !! The test_var_type is a quantity?
!
!
!
   !write(*, *) 'DART quantity ', trim(get_name_for_quantity(test_var_type))
  !! write(*, *) 'long name ', trim(get_long_name(domain_id, test_var_type))
!enddo

end subroutine static_init_model

!-----------------------------------------------------------------------

subroutine fill_gc()

! Subroutine to generate location_types of the cubed sphere grid
! and put them into get_close_type cs_gc, with other derived components.

integer :: c

!SENote: Really don't like all this use of global module storage for communicating among routines
! May want to eliminate some of this.
allocate(cs_locs(ncol), cs_kinds(ncol))

!SENote
write(*, *) 'in fill_gc'
write(*, *) 'ncol ', ncol
write(*, *) 'coarse_grid ', coarse_grid


! CS inputs in degrees.
do c=1,ncol
   cs_locs(c)  = set_location(grid_data%lon%vals(c), grid_data%lat%vals(c), MISSING_R8, VERTISUNDEF)
   cs_kinds(c) = 0
enddo

!SENote: These interfaces to get_close have been changed in Manhattan. Now a single call.
! Initialize cs_gc%maxdist using the maximum grid spacing.
! There will always be at least 2 nodes within 1 coarse_grid in all directions.
!SENotecall get_close_maxdist_init(cs_gc, coarse_grid)
! Use cs_gc%maxdist and node locations to define the rest of cs_gc.
!SENotecall get_close_obs_init(cs_gc, ncol, cs_locs)

call get_close_init(cs_gc, ncol, coarse_grid, cs_locs)
!SENote Test: call get_close_init(cs_gc, ncol, 1000 * coarse_grid, cs_locs)

end subroutine fill_gc

!-----------------------------------------------------------------------

subroutine nc_read_cs_grid_file()

! Read the number of neighbors, corners, centers, a and b coefficients, and x_ax_bearings
! from a netCDF file once for this grid at the beginning of the assimilation.

integer :: nc_file_ID, nc_var_ID, nc_size, n_dims, max_nghbrs, shp(2)
character(len=NF90_MAX_NAME) :: nc_name

! Open the cubed sphere grid relationships file
call nc_check(nf90_open(path=trim(cs_grid_file), mode=nf90_nowrite, ncid=nc_file_ID), &
      'nc_read_cs_grid_file', 'opening '//trim(cs_grid_file))

! learn how many dimensions are defined in this file.
call nc_check(nf90_inquire(nc_file_ID, n_dims), 'nc_read_cs_grid_file', 'inquire n_dims')

! Dimensions written out:
!              name="ncenters",      len = ncenters,      dimid = ncenters_ID), &
!              name="ncorners",      len = ncorners,      dimid = ncorners_ID), &
!              name="max_neighbors", len = max_neighbors, dimid = max_neighbors_ID), &
!              name="ncol",          len = ncol,          dimid = ncol_ID), &
!              name="ncoef_a",       len = 3,             dimid = a_ID), &
!              name="ncoef_b",       len = 2,             dimid = b_ID), &
call nc_check(nf90_inquire_dimension(nc_file_ID, 1, nc_name, ncenters), &
              'nc_read_cs_grid_file', 'inquire for '//trim(nc_name))

call nc_check(nf90_inquire_dimension(nc_file_ID, 2, nc_name, ncorners), &
              'nc_read_cs_grid_file', 'inquire for '//trim(nc_name))

call nc_check(nf90_inquire_dimension(nc_file_ID, 3, nc_name, max_nghbrs), &
              'nc_read_cs_grid_file', 'inquire for '//trim(nc_name))
if (trim(nc_name) /= 'max_neighbors') then
   write(string1, *) trim(cs_grid_file),' max_nghbrs does not match ', trim(cam_template_filename)
   call error_handler(E_ERR,'nc_read_cs_grid_file',string1,source,revision,revdate)
endif
! Check value against the namelist/parameter value.
if (max_nghbrs /= max_neighbors) then
   write(string1, *) trim(cs_grid_file),' max_nghbrs does not match max_neighbors', &
         max_nghbrs,max_neighbors
   call error_handler(E_ERR,'nc_read_cs_grid_file',string1,source,revision,revdate)
endif

call nc_check(nf90_inquire_dimension(nc_file_ID, 4, nc_name, nc_size), &
              'nc_read_cs_grid_file', 'inquire for '//trim(nc_name))

if (nc_size == ncol .and. trim(nc_name) == 'ncol') then
   allocate (corners(ncenters,ncorners),  &
             num_nghbrs           (ncol), &
             centers(max_neighbors,ncol), &
             x_ax_bearings  (ncorners,ncenters), &
             a            (3,ncorners,ncenters), &
             b            (2,ncorners,ncenters)  )
   ! Initialize the grid variables
   num_nghbrs    = MISSING_I
   centers       = MISSING_I
   a             = MISSING_R8
   b             = MISSING_R8
   x_ax_bearings = MISSING_R8

!SENote: I have substituted my_task_id == 0 for the old "output_task0"; confirm
   if (allocated(centers) .and. my_task_id() == 0 .and. print_details) then
      shp = shape(centers)
      write(string1,*) 'Shape of centers = ',shp
      call error_handler(E_MSG,'nc_read_cs_grid_file',string1,source,revision,revdate)
   endif

!SENote: I have substituted my_task_id == 0 for the old "output_task0"; confirm
   if (allocated(corners) .and. my_task_id() == 0 .and. print_details) then
      shp = shape(corners)
      write(string1,*) 'Shape of corners = ',shp
      call error_handler(E_MSG,'nc_read_cs_grid_file',string1,source,revision,revdate)
   endif
else
   write(string1,*) trim(cs_grid_file),' ncol does not match ', trim(cam_template_filename)
   call error_handler(E_ERR,'nc_read_cs_grid_file',string1,source,revision,revdate)
endif

call nc_check(nf90_inq_varid(nc_file_ID, 'num_nghbrs', nc_var_ID), &
                'nc_read_cs_grid_file', 'inq_varid num_nghbrs')
call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, num_nghbrs ), &
                'nc_read_cs_grid_file', 'get_var num_nghbrs')

call nc_check(nf90_inq_varid(nc_file_ID, 'centers', nc_var_ID), &
                'nc_read_cs_grid_file', 'inq_varid centers')
call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, centers ), &
                'nc_read_cs_grid_file', 'get_var centers')

call nc_check(nf90_inq_varid(nc_file_ID, 'corners', nc_var_ID), &
                'nc_read_cs_grid_file', 'inq_varid corners')
call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, corners, &
                           start=(/ 1, 1 /),count=(/ ncenters, ncorners /) ), &
                'nc_read_cs_grid_file', 'get_var corners')

call nc_check(nf90_inq_varid(nc_file_ID, 'a', nc_var_ID), &
                'nc_read_cs_grid_file', 'inq_varid a')
call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, a ), &
                'nc_read_cs_grid_file', 'get_var a')

call nc_check(nf90_inq_varid(nc_file_ID, 'b', nc_var_ID), &
                'nc_read_cs_grid_file', 'inq_varid b')
call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, b ), &
                'nc_read_cs_grid_file', 'get_var b')

call nc_check(nf90_inq_varid(nc_file_ID, 'x_ax_bearings', nc_var_ID), &
                'nc_read_cs_grid_file', 'inq_varid x_ax_bearings')
call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, x_ax_bearings ), &
                'nc_read_cs_grid_file', 'get_var x_ax_bearings')

call nc_check(nf90_close(nc_file_ID), 'nc_read_cs_grid_file', 'closing '//trim(cs_grid_file))

end subroutine nc_read_cs_grid_file



!-----------------------------------------------------------------------

subroutine nc_write_cs_grid_file(cs_grid_file, homme_map_file)

! Write out the number of neighbors, the neighbors, corners, centers, and bearings
! to a netCDF file once for this grid at the beginning of the assimilation.
! Called by static_init_model.

character(len=*), intent(in) :: cs_grid_file
character(len=*), intent(in) :: homme_map_file

integer :: nc_file_ID
integer ::                               &
        ncenters_ID,       centers_var_ID,  &
        ncorners_ID,       corners_var_ID,  &
               a_ID,             a_var_ID,  &
               b_ID,             b_var_ID,  &
            ncol_ID, x_ax_bearings_var_ID,  &
   max_neighbors_ID,    num_nghbrs_var_ID


! Create the file
call nc_check(nf90_create(path=trim(cs_grid_file), cmode=NF90_SHARE, ncid=nc_file_ID), &
              'nc_write_cs_grid_file', 'create '//trim(cs_grid_file))

write(string1,*) trim(cs_grid_file),' is nc_file_ID ',nc_file_ID
call error_handler(E_MSG,'nc_write_cs_grid_file',string1,source,revision,revdate)

! Define the dimensions
call nc_check(nf90_def_dim(ncid=nc_file_ID,                                          &
              name="ncenters",      len = ncenters,      dimid = ncenters_ID), &
              'nc_write_cs_grid_file', 'def_dim ncenters '//trim(cs_grid_file))
call nc_check(nf90_def_dim(ncid=nc_file_ID,                                          &
              name="ncorners",      len = ncorners,      dimid = ncorners_ID), &
              'nc_write_cs_grid_file', 'def_dim ncorners '//trim(cs_grid_file))
call nc_check(nf90_def_dim(ncid=nc_file_ID,                                  &
              name="max_neighbors", len = max_neighbors, dimid = max_neighbors_ID), &
              'nc_write_cs_grid_file', 'def_dim max_neighbors'//trim(cs_grid_file))
call nc_check(nf90_def_dim(ncid=nc_file_ID,                                  &
              name="ncol",          len = ncol,          dimid = ncol_ID), &
              'nc_write_cs_grid_file', 'def_dim ncol '//trim(cs_grid_file))
call nc_check(nf90_def_dim(ncid=nc_file_ID,                                  &
              name="ncoef_a",          len = 3,          dimid = a_ID), &
              'nc_write_cs_grid_file', 'def_dim a '//trim(cs_grid_file))
call nc_check(nf90_def_dim(ncid=nc_file_ID,                                  &
              name="ncoef_b",          len = 2,          dimid = b_ID), &
              'nc_write_cs_grid_file', 'def_dim b '//trim(cs_grid_file))

! Write Global Attributes
call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "title", trim(cs_grid_file)), &
              'nc_write_cs_grid_file',   'put_att title '//trim(cs_grid_file))
call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "model_mod_source", source ), &
              'nc_write_cs_grid_file',   'put_att model_mod_source '//trim(cs_grid_file))
call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "model_mod_revision", revision ), &
              'nc_write_cs_grid_file',   'put_att model_mod_revision '//trim(cs_grid_file))
call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "model_mod_revdate", revdate ), &
              'nc_write_cs_grid_file',   'put_att model_mod_revdate '//trim(cs_grid_file))

call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "elements_per_cube_edge", ne ), &
              'nc_write_cs_grid_file',   'put_att elements_per_cube_edge '//trim(cs_grid_file))
call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "nodes_per_element_edge", np ), &
              'nc_write_cs_grid_file',   'put_att nodes_per_elements_edge '//trim(cs_grid_file))
call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "HommeMapping_file", homme_map_file ), &
              'nc_write_cs_grid_file',   'put_att HommeMapping_file '//trim(cs_grid_file))

! Create variables and attributes.
call nc_check(nf90_def_var(ncid=nc_file_ID, name="num_nghbrs", xtype=nf90_int, &
              dimids=(/ ncol_ID /), varid=num_nghbrs_var_ID),  &
              'nc_write_cs_grid_file', 'def_var num_nghbrs')
call nc_check(nf90_put_att(nc_file_ID, num_nghbrs_var_ID, "long_name", &
              "number of neighbors of each node/column"), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(nc_file_ID, num_nghbrs_var_ID, "units",     "nondimensional"), &
              'nc_write_cs_grid_file', 'units')
call nc_check(nf90_put_att(nc_file_ID, num_nghbrs_var_ID, "valid_range", &
              (/ 1, max_neighbors /)), 'nc_write_cs_grid_file', 'put_att valid_range')

call nc_check(nf90_def_var(ncid=nc_file_ID, name="centers", xtype=nf90_int, &
              dimids=(/ max_neighbors_ID, ncol_ID /), varid=centers_var_ID),  &
              'nc_write_cs_grid_file', 'def_var centers')
call nc_check(nf90_put_att(nc_file_ID, centers_var_ID, "long_name", &
              "cells which use node/column as a corner"), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(nc_file_ID, centers_var_ID, "units",     "nondimensional"), &
              'nc_write_cs_grid_file', 'units')
call nc_check(nf90_put_att(nc_file_ID, centers_var_ID, "valid_range", &
              (/ 1, ncenters /)), 'nc_write_cs_grid_file', 'put_att valid_range')
call nc_check(nf90_put_att(nc_file_ID, centers_var_ID, "missing_value", &
              (/ MISSING_I /)), 'nc_write_cs_grid_file', 'put_att missing_value')

call nc_check(nf90_def_var(ncid=nc_file_ID, name="corners", xtype=nf90_int, &
              dimids=(/ ncenters_ID, ncorners_ID /), varid=corners_var_ID),  &
              'nc_write_cs_grid_file', 'def_var corners')
call nc_check(nf90_put_att(nc_file_ID, corners_var_ID, "long_name", &
              "corners/nodes of each cell "), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(nc_file_ID, corners_var_ID, "units",     "nondimensional"), &
              'nc_write_cs_grid_file', 'units')
call nc_check(nf90_put_att(nc_file_ID, corners_var_ID, "valid_range", &
              (/ 1, ncol /)), 'nc_write_cs_grid_file', 'put_att valid_range')
call nc_check(nf90_put_att(nc_file_ID, corners_var_ID, "missing_value", &
              (/ MISSING_I /)), 'nc_write_cs_grid_file', 'put_att missing_value')

call nc_check(nf90_def_var(ncid=nc_file_ID, name="a", xtype=nf90_double, &
              dimids=(/ a_ID, ncorners_ID, ncenters_ID /), varid=a_var_ID),  &
              'nc_write_cs_grid_file', 'def_var a')
call nc_check(nf90_put_att(nc_file_ID, a_var_ID, "long_name",  &
              "Coefficients of mapping from planar x coord to unit square"), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(nc_file_ID, a_var_ID, "units",     "nondimensional"), &
              'nc_write_cs_grid_file', 'units')
!call nc_check(nf90_put_att(nc_file_ID, a_var_ID, "valid_range", &
!              (/ 1, ncol /)), 'nc_write_cs_grid_file', 'put_att valid_range')
call nc_check(nf90_put_att(nc_file_ID, a_var_ID, "missing_value", &
              (/ MISSING_R8 /)), 'nc_write_cs_grid_file', 'put_att missing_value')


call nc_check(nf90_def_var(ncid=nc_file_ID, name="b", xtype=nf90_double, &
              dimids=(/ b_ID, ncorners_ID, ncenters_ID /), varid=b_var_ID),  &
              'nc_write_cs_grid_file', 'def_var b')
call nc_check(nf90_put_att(nc_file_ID, b_var_ID, "long_name", &
              "Coefficients of mapping from planar y coord to unit square"), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(nc_file_ID, b_var_ID, "units",     "nondimensional"), &
              'nc_write_cs_grid_file', 'units')
!call nc_check(nf90_put_att(nc_file_ID, b_var_ID, "valid_range", &
!              (/ 1, ncol /)), 'nc_write_cs_grid_file', 'put_att valid_range')
call nc_check(nf90_put_att(nc_file_ID, b_var_ID, "missing_value", &
              (/ MISSING_R8 /)), 'nc_write_cs_grid_file', 'put_att missing_value')

call nc_check(nf90_def_var(ncid=nc_file_ID, name="x_ax_bearings", xtype=nf90_double, &
              dimids=(/ ncorners_ID, ncenters_ID /), varid=x_ax_bearings_var_ID),  &
              'nc_write_cs_grid_file', 'def_var x_ax_bearings')
call nc_check(nf90_put_att(nc_file_ID, x_ax_bearings_var_ID, "long_name", &
              "bearing (clockwise from North) from origin node(corner 4) of each mapping to corner 3"), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(nc_file_ID, x_ax_bearings_var_ID, "units",     "radians"), &
              'nc_write_cs_grid_file', 'units')
call nc_check(nf90_put_att(nc_file_ID, x_ax_bearings_var_ID, "valid_range", &
              (/ -PI, PI /)), 'nc_write_cs_grid_file', 'put_att valid_range')
call nc_check(nf90_put_att(nc_file_ID, x_ax_bearings_var_ID, "missing_value", &
              (/ MISSING_R8 /)), 'nc_write_cs_grid_file', 'put_att missing_value')

! Leave define mode so we can fill
call nc_check(nf90_enddef(nc_file_ID), 'nc_write_cs_grid_file', 'enddef '//trim(cs_grid_file))

! sync to disk, but leave open
call nc_check(nf90_sync(nc_file_ID), 'nc_write_cs_grid_file', 'sync '//trim(cs_grid_file))

! Fill the variables
call nc_check(nf90_put_var(nc_file_ID, num_nghbrs_var_ID, num_nghbrs),  &
              'nc_write_cs_grid_file ','put_var num_nghbrs ')
call nc_check(nf90_put_var(nc_file_ID, centers_var_ID, centers),        &
              'nc_write_cs_grid_file ','put_var centers ')
call nc_check(nf90_put_var(nc_file_ID, corners_var_ID, corners),        &
              'nc_write_cs_grid_file ','put_var centers ')
call nc_check(nf90_put_var(nc_file_ID, a_var_ID, a),    &
              'nc_write_cs_grid_file ','put_var a ')
call nc_check(nf90_put_var(nc_file_ID, b_var_ID, b),    &
              'nc_write_cs_grid_file ','put_var b ')
call nc_check(nf90_put_var(nc_file_ID, x_ax_bearings_var_ID, x_ax_bearings),      &
              'nc_write_cs_grid_file ','put_var x_ax_bearings ')

call nc_check(nf90_close(nc_file_ID), 'nc_write_cs_grid_file', 'closing '//trim(cs_grid_file))

end subroutine nc_write_cs_grid_file


!-----------------------------------------------------------------------

subroutine read_cam_2Dint(file_name, cfield, field, num_dim1, num_dim2)

! Read 2d integer field from, e.g., HommeMapping.nc

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

character(len=*),     intent(in)  :: file_name
character(len=*),     intent(in)  :: cfield
integer, allocatable, intent(out) :: field(:,:)
integer,              intent(out) :: num_dim1     !The dimension(s) of cfield
integer,              intent(out) :: num_dim2

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer :: nc_file_ID, nc_var_ID                     !NetCDF variables
integer :: field_dim_IDs(2)                          !Array of dimension IDs for cfield
character(len=NF90_MAX_NAME) :: name_dim1,name_dim2  !Names of dimensions of cfield

field_dim_IDs = MISSING_I                  !Array of dimension IDs for cfield

if (file_exist(file_name)) then
   call nc_check(nf90_open(path=trim(file_name), mode=nf90_nowrite, ncid=nc_file_ID), &
              'read_cam_2Dint', 'opening '//trim(file_name))
!SENote: I have substituted my_task_id == 0 for the old "output_task0"; confirm
   if (print_details .and. my_task_id() == 0) then
      write(string1,*) 'file_name for ',cfield,' is ', trim(file_name)
      call error_handler(E_MSG, 'read_cam_2Dint', string1,source,revision,revdate)
   endif

   ! get field id
   call nc_check(nf90_inq_varid(nc_file_ID, trim(cfield), nc_var_ID), &
              'read_cam_2Dint', 'inq_varid: '//cfield)

   ! get dimension 'id's
   call nc_check(nf90_inquire_variable(nc_file_ID, nc_var_ID, dimids=field_dim_IDs), &
              'read_cam_2Dint', 'inquire_variable: '//cfield)

   ! get dimension sizes
   ! The first spatial dimension is always present.
   call nc_check(nf90_inquire_dimension(nc_file_ID, field_dim_IDs(1), name_dim1, num_dim1 ), &
                 'read_cam_2Dint', 'inquire_dimension: '//name_dim1)
   if (field_dim_IDs(2) /= MISSING_I)  then
      call nc_check(nf90_inquire_dimension(nc_file_ID, field_dim_IDs(2), name_dim2, num_dim2 ), &
                    'read_cam_2Dint', 'inquire_dimension: '//name_dim2)
   else
      num_dim2 = 1
      name_dim2 = 'no2ndDim'
   endif

!SENote: I have substituted my_task_id == 0 for the old "output_task0"; confirm
   if (print_details .and. my_task_id() == 0) then
      write(string1,*) cfield,' dimensions num_dim1, num_dim2 = ',num_dim1, num_dim2
      call error_handler(E_MSG, 'read_cam_2Dint', string1,source,revision,revdate)
   endif
else
   write(string1,'(3A)') 'Required file "',trim(file_name),'" is missing.'
   call error_handler(E_ERR, 'read_cam_2Dint', string1, source, revision, revdate)
endif

! Allocate array, based on size of this variable on the file.
allocate(field(num_dim1,num_dim2))

if (field_dim_IDs(2) /= MISSING_I)  then
   call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, field, start=(/ 1, 1 /), &
                 count=(/ num_dim1, num_dim2 /)), 'read_cam_2Dint', trim(cfield))
else
   call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, field),  &
                  'read_cam_2Dint', trim(cfield))
endif

call nc_check(nf90_close(nc_file_ID), 'read_cam_2Dint', 'closing '//trim(file_name))

end subroutine read_cam_2Dint

!-----------------------------------------------------------------------

subroutine create_cs_grid_arrays()

! Subroutine to create arrays of relationships between cubed sphere nodes (corners)
! and cell centers, including bearings between nodes.
! These will be used to identify the cell containing an observation.
! The relationships read from HommeMapping.nc will be augmented.
! All will be stored in global storage, and written to a new file for
! subsequent use.

! Local variables
integer  :: sh_corn(4), n(4)     ! Shifted corners to put closest at the origin.
integer  :: col, nbr, c, cent            ! Indices for loops.
integer  :: num_n, min_ind(1)
real(r8) :: dist, angle
real(r8) :: bearings(3), x_planar(3), y_planar(3)

! ncol = number of nodes/corners/grid points.  Global storage.
! corners = the names of the corners associated with each cell center
! neighbors = the nodes around each node which partner to make the sides of the cells
!             which may contain an observation.

! Get array of corner nodes/columns which define the cells (identified by 'center').
if (file_exist(homme_map_file)) then
   call read_cam_2Dint(homme_map_file, 'element_corners', corners,ncenters,ncorners)

   if (ncenters /= (ncol -2) ) then
      write(string1, *) trim(homme_map_file),' ncenters inconsistent with ncol-2 ', ncenters, ncol
      call error_handler(E_ERR,'create_cs_grid_arrays',string1,source,revision,revdate)
   endif

      allocate(num_nghbrs           (ncol), &
               centers(max_neighbors,ncol), &
               a          (3,ncorners,ncenters),   &
               b          (2,ncorners,ncenters),   &
               x_ax_bearings(ncorners,ncenters))

   num_nghbrs    = 0
   centers       = MISSING_I
   a             = MISSING_R8
   b             = MISSING_R8
   x_planar      = MISSING_R8
   y_planar      = MISSING_R8
   x_ax_bearings = MISSING_R8
else
   write(string1, *) 'CAM-SE grid file "',trim(homme_map_file),'" can not be found '
   call error_handler(E_ERR,'create_cs_grid_arrays',string1,source,revision,revdate)
endif

! Invert the element_corners array to compile all of the neighbors of each node (corner).
! Loop over HommeMapping cell centers.
Quads: do cent = 1,ncenters
   Corns: do c = 1,4
      ! Get the node numbers that define this cell
      ! and shift (rotate) them to create a separate mapping for each corner/node.
      ! Shift the section of corners 1 place to the 'left'/lower for the first corner,
      ! 2 for the 2nd, etc.  This will put the node closest to the ob in position 4
      ! (of the shifted corners). Then the (x,y) origin will the the closest node,
      ! and the indexing of the a,b,x_ax_bearing arrays will be easy.
      ! Shifting preserves the order of the corners as we go around the cell (clockwise).
      sh_corn = cshift(corners(cent,:), c)

      ! Check a few cells for corner consistency.
      if (print_details .and. sh_corn(4) < 10) then
         write(string1,'(A,5I7,/31X,4I7)') 'c, 4 corners, shifted = ', &
              c,(corners(cent,nbr),nbr=1,4),(sh_corn(nbr),nbr=1,4)
         call error_handler(E_MSG,'create_cs_grid_arrays',string1,source,revision,revdate)
      endif

      ! Increment the number of neighbors (and centers ) this corner(node) has.
      ! sh_corn(4) is used for all cases because the corner we're working on always
      ! ends up in that position, when c is incremented, then the corners are shifted.
      n(c) = num_nghbrs(sh_corn(4)) + 1

      ! Update the number of neighbors of each corner of the cell,
      num_nghbrs(sh_corn(4)) = n(c)

      ! Store the info that this center is associated with 4 node->neighbor pairs.
      centers(n(c),sh_corn(4)) = cent

      ! Define the planar coordinates for this center/cell and this corner/node.
      ! The 4th corner is the origin, and the cell side from the 4th to the 3rd is
      ! the x-axis of this cell's coordinate system for this corner.
      ! This is established in the definition of bearings().
      ! This choice makes mapping coefficients a(0) and b(0) = 0 (see below).
      ! It also helps make the indexing of bearings easy to use and store.

      ! Check a few cells for corner consistency
      if (print_details .and. sh_corn(4) < 10) then
         write(string1,'(A,3F10.6)') 'lon1, lat1 = ', lon_rad(sh_corn(4)), lat_rad(sh_corn(4))
         call error_handler(E_MSG, 'create_cs_grid_arrays', string1, source, revision, revdate)
      endif

      ! Descend through neighbors so that bearings(3) is already defined when needed at loop end.
      do nbr = 3,1,-1
         ! Bearings from the current origin node of the cell to the other 3.
         bearings(nbr) = bearing(lon_rad(sh_corn(4)),   lat_rad(sh_corn(4)),  &
                                 lon_rad(sh_corn(nbr)), lat_rad(sh_corn(nbr)) )

         dist = get_dist(cs_locs(sh_corn(4)), cs_locs(sh_corn(nbr)), 0, 0, .true.)

         if (sh_corn(4) < 10) then
            write(string1,'(A,3F10.6)') 'create_cs_grid:    lon2, lat2, bearing = ', &
                 lon_rad(sh_corn(nbr)), lat_rad(sh_corn(nbr)), bearings(nbr)
            call error_handler(E_MSG, 'create_cs_grid_arrays', string1, source, revision, revdate)
         endif

         ! This difference order looks wrong, but we need to change the sign of angles from the
         ! clockwise direction used by bearings to the counterclockwise direction used by
         ! trig functions.
         angle = bearings(3) - bearings(nbr)

         ! Normalize to -PI < angle <= PI.
         angle = mod(angle,PI) - PI*int(angle/PI)

         ! Set the planar location of this corner/node.
         x_planar(nbr) = dist * cos(angle)
         y_planar(nbr) = dist * sin(angle)

      enddo

      ! Store the baseline for use when interpolating to an ob location.
      x_ax_bearings(c,cent) = bearings(3)

      ! Define another bearings array to allow coord_ind_cs to find the right
      ! cell around the closest node by using a search through bearings,
      ! rather than a call to unit_square_location.
      !   Propagate sort_bearings to writing and reading of HommeMapping_cs_grid.nc file.
      !   real(r8), allocatable, :: sort_bearings(max_neighbors,ncol)
      !   sort_bearings(n(c),sh_corn(4)) = bearings(3)
      ! But see ordering of bearings, commented out below.

      ! Define quantities used to map the planar coordinate system of each cell
      ! to the unit square coordinate system.

      ! I'll use the mapping from planar space (x,y) to unit square space (l,m):
      ! x = a0 + a1*l*m + a2*m + a3*l
      ! y = b0 + b1*l*m + b2*m + b3*l
      ! The 4 corners (x,y) can be mapped to the four corners (l,m) to yield 4 equations.
      ! This can be written as vec_x = mat_A * vec_a^T.
      ! Then AI is the inverse of the mapping from physical space to the unit square space,
      ! and the coefficients of the mapping, aN and bN, can be calculated from:
      !    a = matmul(AI,x_planar(0:3))
      !    b = matmul(AI,y_planar(0:3))
      ! But the mapping from (lon,lat) to (x,y) space put corner "4" of the cell at (x4,y4) = (0,0)
      ! and corner "3" at (x3,y3) = (d3,0)
      ! This ends up making a0 = b0 = b3 = 0, and the equations simplify to the point
      ! that it doesn't make sense to encode this tranformation in a matrix.
      ! Replace matrix and matmul with simpler direct equations:
      a(3,c,cent) = x_planar(3)
      a(2,c,cent) = x_planar(1)
      a(1,c,cent) = x_planar(2) - x_planar(1) - x_planar(3)
      b(2,c,cent) = y_planar(1)
      b(1,c,cent) = y_planar(2) - y_planar(1)

      if (cent < 10) then
         write(string1,'(A,1p4E12.4)') 'create_cs_grid_arrays: a = ',(a(nbr,c,cent),nbr=1,3)
         write(string2,'(A,1p4E12.4)') 'create_cs_grid_arrays: b = ',(b(nbr,c,cent),nbr=1,2)
         call error_handler(E_MSG, 'create_cs_grid_arrays', string1, source, revision, revdate,text2=string2)
      endif

      if (a(3,c,cent)* a(2,c,cent) *a(1,c,cent) == 0.0_r8) then
         write(string1,'(A,2I8,A,1p3E12.4)') 'a(:,',c,cent,') = ',(a(nbr,c,cent),nbr=1,3)
         write(string2,'(A,(6X,2F10.6))') 'create_cs_grid:    lon, lat for nghr=1-3  = ', &
              (lon_rad(sh_corn(nbr)), lat_rad(sh_corn(nbr)), nbr=1,3)
         write(string3,'(A,5I7,/31X,4I7)') 'c, 4 corners, shifted = ', &
               c,(corners(cent,nbr),nbr=1,4),(sh_corn(nbr),nbr=1,4)
         call error_handler(E_MSG, 'create_cs_grid_arrays', string1, source, revision, revdate,text2=string2,text3=string3)
      endif

      if (b(2,c,cent) *b(1,c,cent) == 0.0_r8) then
         write(string1,'(A,2I8,A,1p3E12.4)') 'b(:,',c,cent,') = ',(b(nbr,c,cent),nbr=1,2)
         write(string2,'(A,(6X,2F10.6))') 'create_cs_grid:    lon, lat for nghr=1-3  = ', &
              (lon_rad(sh_corn(nbr)), lat_rad(sh_corn(nbr)), nbr=1,3)
         write(string3,'(A,5I7,/31X,4I7)') 'c, 4 corners, shifted = ', &
               c,(corners(cent,nbr),nbr=1,4),(sh_corn(nbr),nbr=1,4)
         call error_handler(E_MSG, 'create_cs_grid_arrays', string1, source, revision, revdate,text2=string2,text3=string3)
      endif

   enddo Corns
enddo Quads


! Check that all nodes have at least 3 neighbors and no more than 6.
do col = 1,ncol
   if (num_nghbrs(col) < 3 .or. num_nghbrs(col) > max_neighbors) then
      write(string1,'(A,I6,A,6I8)') 'num_nghbrs(',col,') <3 or >6: ', num_nghbrs(col)
      call error_handler(E_ERR,'create_cs_grid_arrays',string1,source,revision,revdate)
   endif
enddo

! There's code in earlier versions of model_mod to
! reorder the neighbors so that they are sequential around each node
! to make the search for the cell containing an ob faster.

return

end subroutine create_cs_grid_arrays


!-----------------------------------------------------------------------

real function bearing(lon1,lat1,lon2,lat2)

! Calculate the direction along the great circle from point 1 on a sphere
! to point 2, relative to north.
! All inputs should have units of radians.
! Output is radians.
! From http://www.movable-type.co.uk/scripts/latlong.html

real(r8), intent(in)    :: lon1,lat1, lon2,lat2

real(r8) :: lon1c,lon2c, cos_lat2, del_lon

real(r8), parameter :: half_PI = PI*0.5_r8

! Make sure the poles are handled consistently:
! If the pole point is the origin point, and the longitude of the pole point is
! defined as 0.0, then the bearing to a nearby point will = the longitude of the point.
! This is consistent/continuous with the bearing from points extremely near
! the pole.
if (half_PI - abs(lat1) < epsilon(lat1)) then
   lon1c = 0.0_r8
else
   lon1c = lon1
endif
if (half_PI - abs(lat2) < epsilon(lat2)) then
   lon2c = 0.0_r8
else
   lon2c = lon2
endif

cos_lat2 = cos(lat2)
del_lon  = lon2c - lon1c

! Normalize del_lon to -pi<=angle<=pi.
del_lon = mod(del_lon,PI) - PI*int(del_lon/PI)
bearing = atan2(cos_lat2*sin(del_lon),  &
                cos(lat1)*sin(lat2) - sin(lat1)*cos_lat2*cos(del_lon) )

end function bearing

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
! SENote: this is totally different for SE.

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

!SENote variable declarations follow
real(r8) :: my_lon, my_lat, my_vert

!SENote: Following implemented for SE

! full 2d fields are returned with lon/lat/level.
! 1d fields are either surface fields, or if they
! are column integrated values then they are 'undefined'
! in the vertical.

! All fields share the same first coordinate into the column list
my_lon = grid_data%lon%vals(i)
my_lat = grid_data%lat%vals(i)
! For SE 3d spatial fields have a 2d storage
if(nd == 2) then
   my_vert = j
   get_location_from_index = set_location(my_lon, my_lat, my_vert, VERTISLEVEL)
elseif(nd == 1) then
   ! setting the vertical value to missing matches what the previous
   ! version of this code did.  other models choose to set the vertical
   ! value to the model surface elevation at this location:
   !   use_vert_val  = phis(lon_index, lat_index) / gravity not available in SE
   my_vert = MISSING_R8
   ! Add any 2d surface fields to this function
   if(is_surface_field(q)) then
      get_location_from_index = set_location(my_lon, my_lat, my_vert, VERTISSURFACE)
   else
      get_location_from_index = set_location(my_lon, my_lat, my_vert, VERTISUNDEF)
   endif
else
   write(string1, *) 'state vector field not 1D or 2D and no code to handle other dimensionity'
   write(string2, *) 'dimensionality = ', nd, ' quantity type = ', trim(get_name_for_quantity(q))
   call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2)
endif
     
end function get_location_from_index


!-----------------------------------------------------------------------
!> this routine converts the column and level index values and a quantity into a state vector
!> offset and gets the ensemble of state values for that offset.  this only
!> gets a single vertical location - if you need to get values which might 
!> have different vertical locations in different ensemble members
!> see get_se_values_from_varid() below.

subroutine get_se_values_from_single_level(ens_handle, ens_size, qty, column_index, lev_index, &
                                        vals, my_status)
type(ensemble_type), intent(in) :: ens_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: qty
integer,             intent(in) :: column_index
integer,             intent(in) :: lev_index
real(r8),            intent(out) :: vals(ens_size)
integer,             intent(out) :: my_status

character(len=*), parameter :: routine = 'get_se_values_from_single_level:'

integer :: varid
integer(i8) :: state_indx

varid = get_varid_from_kind(domain_id, qty)
if (varid < 0) then
   vals(:) = MISSING_R8
   my_status = 12
   return
endif

state_indx = get_dart_vector_index(column_index, lev_index, no_third_dimension, domain_id, varid)

!SENote: Not clear we need error checks like this for things that should never happen.
if (state_indx < 1 .or. state_indx > get_domain_size(domain_id)) then
   write(string1, *) 'state_index out of range: ', state_indx, ' not between ', 1, get_domain_size(domain_id)
   call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2,text3='should not happen')
endif
vals(:) = get_state(state_indx, ens_handle)

my_status = 0

end subroutine get_se_values_from_single_level

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

!SENote modified for se form fv base
subroutine get_se_values_from_varid(ens_handle, ens_size, corner_index, lev_index, varid, &
                                 vals, my_status)
type(ensemble_type), intent(in)  :: ens_handle
integer,  intent(in)  :: ens_size
integer,  intent(in)  :: corner_index
integer,  intent(in)  :: lev_index(ens_size)
integer,  intent(in)  :: varid
real(r8), intent(out) :: vals(ens_size)
integer,  intent(out) :: my_status(ens_size)

integer(i8) :: state_indx
integer  :: i, j
real(r8) :: temp_vals(ens_size) 
logical  :: member_done(ens_size)

character(len=*), parameter :: routine = 'get_se_values_from_varid:'

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
   state_indx = get_dart_vector_index(corner_index, lev_index(i), no_third_dimension, domain_id, varid)

   !SENote: Do we need error checks like this?
   if (state_indx < 0) then
      write(string1,*) 'Should not happen: could not find dart state index from '
      write(string2,*) 'corner and lev index :', corner_index, lev_index
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

end subroutine get_se_values_from_varid


!-----------------------------------------------------------------------
!> this is just for 3d fields

subroutine get_se_values_from_nonstate_fields(ens_handle, ens_size, column_index, &
                                           lev_index, obs_quantity, vals, my_status)
type(ensemble_type),  intent(in)  :: ens_handle
integer,              intent(in)  :: ens_size
integer,              intent(in)  :: column_index
integer,              intent(in)  :: lev_index(ens_size)
integer,              intent(in)  :: obs_quantity
real(r8),             intent(out) :: vals(ens_size)
integer,              intent(out) :: my_status(ens_size)

integer  :: imember
real(r8) :: vals_array(ref_nlevels,ens_size)

character(len=*), parameter :: routine = 'get_se_values_from_nonstate_fields:'

vals(:) = MISSING_R8
my_status(:) = 99

select case (obs_quantity) 
   case (QTY_PRESSURE)
      call cam_se_pressure_levels(ens_handle, ens_size, column_index, ref_nlevels, &
                               vals_array, my_status)
      if (any(my_status /= 0)) return

      do imember=1,ens_size
         vals(imember) = vals_array(lev_index(imember), imember)
      enddo

   case (QTY_VERTLEVEL)
      vals(:)      = lev_index(:)
      my_status(:) = 0

!SENote: Turns out there was no height localization for non-height vertical obs in Manhattan of Classic
!SENote: At present there is no QTY_GEOMETRIC_HEIGHT here as needed to convert to Height
!SENote, what did this look like in get_values_from_nonstate_fields?

   case default
      write(string1,*)'contact dart support. unexpected error for quantity ', obs_quantity
      call error_handler(E_ERR,routine,string1,source,revision,revdate)

end select

end subroutine get_se_values_from_nonstate_fields

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
!SENote: need additional local variables
integer :: closest, cell_corners(4), i
!SENote: these legacy variable names for fractions in the interp need to be more informative
real(r8) :: l, m
type(location_type) :: location_copy


if ( .not. module_initialized ) call static_init_model


! Successful istatus is 0
interp_vals(:) = MISSING_R8
istatus(:)     = 99

! do we know how to interpolate this quantity?
call ok_to_interpolate(obs_qty, varid, domain_id, status1)

if (status1 /= 0) then  
   if(debug_level > 12) then
      write(string1,*)'did not find observation quantity ', obs_qty, ' in the state vector'
      call error_handler(E_MSG,routine,string1,source,revision,revdate)
   endif
   istatus(:) = status1   ! this quantity not in the state vector
   return
endif

! unpack the location type into lon, lat, vert, vert_type
lon_lat_vert = get_location(location)
which_vert   = nint(query_location(location)) 

! if we are avoiding assimilating obs above a given pressure, test here and return.
if (discarding_high_obs) then
   call obs_too_high(lon_lat_vert(3), which_vert, status1)
   if (status1 /= 0) then
      istatus(:) = status1
      return
   endif
endif


!SENote
! Do the interpolation here
! First step, find the columns of the four 'corners' containing the location
! CONFIRM that a qty has replaced a kind 
! SENote2: In the CLASSIC, there is a possibility that the cell_corner was already found and this call can
! be skipped. Understand that and implement as needed.
! Also note that cannot pass location directly because it is intent(inout) in coord_ind_cs.
! Not clear that this will eventually be relevant, so look at changing it ther
location_copy = location
call coord_ind_cs(location_copy, obs_qty, .false., closest , cell_corners, l, m)


!SENote: Now work on the vertical conversions and getting the vertical index for each ensemble member
! Can model this on get_quad_vals which does something similar for the FV
! Just need to send the indices of the corners instead of pair of lon and lat indices
call get_se_quad_vals(state_handle, ens_size, varid, obs_qty, cell_corners, &
                   lon_lat_vert, which_vert, quad_vals, istatus)


! Then interpolate horizontally to the (lon,lat) of the ob.
! The following uses Jeff's recommended 'generalized quadrilateral interpolation', as in
! http://www.particleincell.com/2012/quad-interpolation/.
! Most of the work is done in create_cs_grid_arrays() and coord_ind_cs().

! Interpolate from the cell's corners to the ob location on the unit square.
! This is done by weighting the field at each corner by the rectangular area
! ((l,m) space) diagonally across the ob location from the corner.
! AKA 'linear area weighting'.

! SENote: The status block needs to be worked on to be consistent with CLASSIC
!SENote: It is commented out for now but should be put back in.
! The vals indices are consistent with how mapping of corners was done,
! and how cell_corners was assigned.
!SENoteif (vstatus == 1) then
   !SENoteif (print_details) then
      !SENotewrite(string1,'(A,2F10.6,1pE20.12)') 'istatus = 1, no interpolation'
      !SENotecall error_handler(E_MSG, 'interp_cubed_sphere', string1)
   !SENoteendif
   !SENotereturn
!SENoteelse
   !SENoteif (abs(lon_lat_lev(2)) > max_obs_lat_degree) then
      !SENote! Define istatus to be a combination of vstatus (= 0 or 2 (for higher than highest_obs...))
      !SENote! and whether the ob is poleward of the limits set in the namelist (+ 4).
      !SENoteistatus = 10*vstatus + 4
   !SENoteelse
      !SENoteistatus = vstatus
   !SENoteendif
!SENoteendif

! SENote: could this be done in vector notation?
do i = 1, ens_size
interp_vals(i) = quad_vals(2, i) *           l *          m &
           + quad_vals(1, i) * (1.0_r8 - l)*          m &
           + quad_vals(4, i) * (1.0_r8 - l)*(1.0_r8 - m) &
           + quad_vals(3, i) *           l *(1.0_r8 - m)
end do

if (print_details ) then
   !SENote: Originally for scalar, rather than ensemble size output. Should modify
   write(string1,'(A,2F10.6,1pE20.12)') ' l,m, interpolated vals = ', &
         l,m,interp_vals(1)
   call error_handler(E_MSG, 'interp_cubed_sphere', string1)
endif




if (using_chemistry) &
   interp_vals = interp_vals * get_volume_mixing_ratio(obs_qty)

! all interp values should be set by now.  set istatus
istatus(:) = 0

end subroutine model_interpolate

!-----------------------------------------------------------------------

subroutine coord_ind_cs(obs_loc, obs_kind, closest_only, closest , cell_corners, l, m)

! Find the node closest to a location, and the possibly the corners of the cell which contains 
! the location.

! Variables needed by loc_get_close_obs:
!SENote: Had to change this to intent(inout) because Manhattan loc_get_close_obs requires it
type(location_type),  intent(inout)  :: obs_loc
integer,              intent(in)  :: obs_kind
logical,              intent(in)  :: closest_only
integer,              intent(out) :: closest
integer,              intent(out) :: cell_corners(4)
real(r8),             intent(out) :: l
real(r8),             intent(out) :: m

! Output from loc_get_close_obs
integer  :: num_close

! It would be nice if these could be smaller, but I don't know what number would work.
! It has to be large enough to accommodate all of the grid points that might lie
! within 2xcutoff; resolution and location dependent.
! The size must be specified here; (:) yields an error, and 'allocatable' doesn't help.
integer, allocatable  :: close_ind(:)
real(r8), allocatable :: dist(:)

! Local Variables
! dist_# in radians (Can't be initialized here or they will get the 'save' property,
! and will not be reset during subsequent entries to this subroutine.)
real(r8) :: dist_1, dist_2
real(r8) :: lon_lat_lev(3)
integer  :: k, k1, k2, closest2, origin
logical  :: found_cell

lon_lat_lev = get_location(obs_loc)

! See whether this obs_ is a state variable.
! This could be done by 2 calls to minloc(dist), with the 2nd call using a mask
! to prevent finding the closest, which was found in the first call.
! But would those 2 intrinsic searches through dist be faster than my 1 explicit search?

if (closest_only) then
   ! Use xyz/cartesian coordinates to quickly find the closest node.
   ! If convert_vert only needs the closest node, don't find the l,m weights.
   closest = find_closest_node(lon_lat_lev(2), lon_lat_lev(1))

   ! Can return without deallocating close_ind and dist
   ! because they haven't been allocated yet.
   return
endif

! Allocate space for the potentially close nodes.
allocate(close_ind(ncol), dist(ncol))

! Look for the 2 closest nodes, using slower way of getting all of the close obs
! and searching for the 2 closest.
! --------------
! FIXME: Nancy has a location_xyz:find_closest_???? which will return the N closest points,
! which may be significantly faster than threed_sphere/location_mod.f90:get_close_obs.
! --------------
! FIXME; can the closest node not be a corner of the containing cell in grids generated by SQuadGen?
! --------------
! For a refined grid (from 1 degree to 1/8 degree) loc_get_close_obs is going to return lists
! that are 64x larger in the refined region than in the coarse region

!   obs_'kind' is passed to location.f90:get_close_obs.
!   There it is passed to only get_dist, which only uses it if special_vert_norm is used,
!   and gc%special_maxdist.
!   Model_mod is not using either of those.
!SENote IMPORTANT: This only works with approximate_distance = .false. Somehow this must be overridden from
! Namelist or documented.
call get_close(cs_gc, obs_loc, obs_kind, cs_locs, cs_kinds, &
                       num_close, close_ind, dist)


dist_1 = 10.0_r8
dist_2 = 10.0_r8
closest = MISSING_I
k1 = MISSING_I

! Keep track of k1, k2, and distances in this search.
! Assign closest and closest2 afterwards.
if (num_close <= 0) then
   write(string1,*) 'Unusable num_close, obs_kind : ',num_close, obs_kind
   call write_location(0, obs_loc, charstring=string2)
   write(string3,*) 'dist(1) = ',dist(1)
   call error_handler(E_ERR, 'coord_ind_cs', string1,source,revision,revdate,text2=string2, text3=string3)
endif

do k = 1,num_close
   if (dist(k) < dist_2) then
      ! Replace 2nd with new one.
      k2 = k
      dist_2 = dist(k)
      if (dist_2 <= dist_1) then
         ! Switch 1st and new 2nd.   '<=' To make sure k1 is filled, even for the first k.
         k2 = k1
         k1 = k
         dist_2 = dist_1
         dist_1 = dist(k)
      endif
   endif
enddo
closest  = close_ind(k1)

if (k2 == MISSING_I) then
   write(string1,'(A)') 'Did not find a second closest node to ob:'
   write(string2,'(A,3F10.2,3I6,1p2E12.4)')                    &
        'lon_lat_lev, obs_kind, num_close, closest, dist_1, dist_2 = ', &
         lon_lat_lev, obs_kind, num_close, closest, dist_1, dist_2
   call write_location(0, cs_locs(closest), charstring=string3)
   string3 = 'closest node location = '//string3
   call error_handler(E_ERR, 'coord_ind_cs', string1,source,revision,revdate,text2=string2,text3=string3)
else
   closest2 = close_ind(k2)
endif

! Find the cell which contains the ob.
! First search the cells which have 'closest' as 1 corner.
! If that fails, search the cells around closest2.
! The search consists of passing the ob location to unit_square_location
! and letting it determine whether the ob location maps into the unit square.

! Initial value of success flag.
found_cell = .false.

! FIXME; debug in verify_namelist
! write(string1,*) 'STARTING Cloop num_nghbrs = ',num_nghbrs(closest)
! call error_handler(E_MSG,'coord_ind_cs',string1,source,revision,revdate)

Cloop: do k=1,num_nghbrs(closest)
   ! centers(k,closest) refers to the cell center name associated with neighboring node k
   ! of the closest node.  It is used to retrieve mapping coefficients for the cell being tested.

   call unit_square_location(centers(k,closest), closest, obs_loc,          &
                             lon_lat_lev(1),lon_lat_lev(2), found_cell, origin, l,m)
   if (found_cell) exit Cloop
enddo Cloop

! Try the 2nd closest point, if the first failed.
if ((.not.found_cell) .and. closest2 /= MISSING_I) then

   Second_closest: do k=1,num_nghbrs(closest2)
      call unit_square_location(centers(k,closest2), closest2, obs_loc,         &
                                lon_lat_lev(1),lon_lat_lev(2), found_cell, origin, l,m)
      if (found_cell) then
         ! Put '2nd closest' information into 'closest'.
         dist_1 = dist_2
         closest = closest2

         write(string1,'(A,2F10.7,2I8,1p2E12.4)') &
              'Using 2nd closest node to the ob: l, m, closest2, origin2 = ', &
              l, m, closest, origin
         call error_handler(E_MSG, 'coord_ind_cs', string1,source,revision,revdate)

         exit Second_closest
      endif
   enddo Second_closest
endif

if (found_cell) then
   ! Need to shift corners according to which was chosen as the origin corner
   ! in num_nghbrs loop, above.  The weighted interp calculation assumes, as in
   ! the create_cs_grid_arrays mapping scheme, that the origin node is corner 4.
   cell_corners(1:4) = cshift(corners(centers(k,closest),1:4), origin)

else
   ! Both closest nodes failed; abort
   write(string1, '(A,2I8,A,2F10.4)') &
         'Neither of the 2 closest nodes ',  closest,closest2, &
         ' is a corner of the cell containing ob at ', lon_lat_lev(1),lon_lat_lev(2)
   call error_handler(E_ERR, 'coord_ind_cs', string1,source,revision,revdate)
endif

deallocate(close_ind, dist)

end subroutine coord_ind_cs


!-----------------------------------------------------------------------

subroutine unit_square_location(cell, closest, location, lon_o,lat_o, found_cell, origin, l,m)

! Subroutine based on http://www.particleincell.com/2012/quad-interpolation/.
! The idea is to derive a mapping from any convex quadrilateral(x,y) onto a unit square (l,m).
! Also map the location of the ob onto that square.
! This is a bilinear interpolation;
! x = a0 + a1*l*m + a2*m + a3*l
! y = b0 + b1*l*m + b2*m + b3*l
! so does not take into account the curvature of the quadrilateral on the sphere.
!
! That has been handled by the intermediate mapping from (lon,lat) to a flat planar
! coordinate system (x,y). The locations of the corners/nodes are converted to
! the distances and directions from one node to the other three.  See create_cs_grid_arrays.
! Distances and directions relative to the origin node are preserved, but distances and
! directions between 2 non-origin points are slightly distorted.
! Even these small errors are avoided by defining a planar coordinates system for each corner
! of each cell.
! Then the ob is never near the 'far edges', where distortion could be a problem.
!
! A higher order method exists (Nagata 2005: Simple Local Interpolation of Surfaces
! Using Normal Vectors) to map curved quadrilaterals onto the unit square,
! but the inverse map cannot be done analytically(?), so is not developed here.

integer,             intent(in)    :: cell
integer,             intent(in)    :: closest
type(location_type), intent(in)    :: location
real(r8),            intent(in)    :: lon_o
real(r8),            intent(in)    :: lat_o
logical,             intent(inout) :: found_cell
integer,             intent(out)   :: origin
real(r8),            intent(out)   :: l
real(r8),            intent(out)   :: m

! Observation location in the planar space.
real(r8) :: x_o, y_o

real(r8) :: angle, d, bearing_o  ! Locations in polar coordinate space (bearing,distance).
real(r8) :: aa, bb, cc           ! Coefficients of quadratic equation for m.
real(r8) :: det, m1, m2          ! Determinant and roots.
logical  :: neg_root             ! helpful logical variable to store usefulness of the -root.
real(r8) :: m_neg, l_neg         ! Potential alternate solutions to the m quadratic equation
integer  :: oc(1)

m1    = MISSING_R8  ! first  root returned by solve_quadratic
m2    = MISSING_R8  ! second root returned by solve_quadratic
l     = MISSING_R8  ! unit square abscissa ('x' coord)
m     = MISSING_R8  ! unit square ordinate ('y' coord)
l_neg = MISSING_R8  ! same but for the negative root of the m quadratic equation.
m_neg = MISSING_R8  ! same
neg_root = .false.

! Map the location of the ob into the planar space

! Figure out which corner (1,2,3 or 4) of cell is the closest to the ob,
! by comparing the names of the corners to the name of the node/corner closest
! to the ob, which was passed in.
! Used to get the correct x_ax_bearing and a and b coeffs (from the cs_grid_file).

oc = minloc(corners(cell,:), mask = (corners(cell,:) == closest))
origin = oc(1)

! The bearing of the observation relative to the origin/closest corner.
bearing_o = bearing(lon_rad(closest),lat_rad(closest),lon_o*DEG2RAD,lat_o*DEG2RAD )

! Calculate the difference of the ob bearing from x_axis of this cell.
! The order is opposite of what might be expected because bearings are measured clockwise,
! while angles are measured counterclockwise.
angle = x_ax_bearings(origin,cell) - bearing_o

! Normalize angle to -pi<angle<pi.
angle = mod(angle,PI) - PI*int(angle/PI)

! Calculate the distance from this cell's origin
! and then the (x,y) coordinates of the observation.
d = get_dist(cs_locs(closest), location, no_vert=.true.)
x_o = d * cos(angle)
y_o = d * sin(angle)

! Coefficients of the quadratic equation for m for this cell.
aa =   a(1,origin,cell)*b(2,origin,cell) &
     - a(2,origin,cell)*b(1,origin,cell)

bb =   a(3,origin,cell)*b(2,origin,cell) &
     - a(1,origin,cell)*y_o              &
     + b(1,origin,cell)*x_o

cc = - a(3,origin,cell)*y_o

! Calculate m from the binomial equation, given the quadratic equation coefficients.
call solve_quadratic(aa,bb,cc,m1,m2)
! newFIXME: Simplify this subroutine?
! There can only be one mapping from the (x,y) space to the unit square.
! One of the (potentially) 2 'm's generated here should be eliminated by the l calculation.

if (m1 == MISSING_R8 .and. m2 == MISSING_R8) then
   ! determinant was < 0.
   write(string1,'(A,I6,1X,1p4E12.4)') 'm b^2-4ac <0: cell, angle, d, x_o, y_o',cell, angle, d, x_o, y_o
   call error_handler(E_MSG, 'unit_square_location', string1,source,revision,revdate)
   write(string1,'(A,1X,1p4E12.4)') '   a: a(1)*b(2) - a(2)*b(1) : ',  &
                            a(1,origin,cell),b(2,origin,cell), &
                            a(2,origin,cell),b(1,origin,cell)
   write(string2,'(A,1X,1p4E12.4)') '   b: ', bb
   write(string3,'(A,1X,1p4E12.4)') '   c: a(3)*y_o : ', a(3,origin,cell), y_o
   call error_handler(E_MSG, 'unit_square_location', string1,source,revision,revdate, &
                      text2=string2, text3=string3)
else
   if (aa > 0.0_r8) then
      ! Only m values (roots) between 0 and 1 mean that the ob is in this cell.
      if (m1 >=0 .and. m1 <= 1) then
         m = m1
      elseif (m2 >=0 .and. m2 <= 1) then
         m = m2
      else
         ! Neither root is a map.  Leave m as MISSING_R8
      endif
   elseif (m1 /= MISSING_R8 .and. m2 == MISSING_R8 ) then
      ! Cell is square; solved the linear equation    m*bb + cc = 0
      m = m1
   else
      ! aa < 0; Either both or neither roots yield m>0.
      ! Start with the +root.
!       m = (-bb + sqrt_det)/(2.0_r8*aa)
      m = m1

      if (bb > 0.0_r8) then
         ! Both roots yield m > 0.
         if (m > 1.0_r8) then
            ! The +root didn't yield a usable m.  Try the -root.
            ! m = (-bb - sqrt_det)/(2.0_r8*aa)
            m = m2
         else
            ! It could be that both roots yield a usable m.  Keep track of both
            ! (for testing/debugging only).
            m_neg = m2
         endif

      elseif (bb < 0.0_r8) then
         ! aa < 0 and bb < 0 yields no roots with m>0.
         write(string1,'(A,I6,A)') 'aa < 0 and bb < 0: It appears that cell ',cell,           &
              ' is a highly distorted quadrilateral'
         write(string2,'(A)')                                                                 &
              'and no mapping is possible.  bb = a(3)*b(2) - a(1)*y_o + b(1)*x_o: '
         write(string3,'(1p,(3X,2E12.4))')                                                     &
              a(3,origin,cell),b(2,origin,cell),                                              &
              a(1,origin,cell),y_o,                                                           &
              b(1,origin,cell),x_o
         call error_handler(E_ERR, 'unit_square_location', string1,source,revision,revdate,  &
                            text2=string2, text3=string3)

      elseif (bb == 0.0_r8) then
         ! aa < 0 and bb = 0  should be excluded by the non-negativeness test on det, above.
         write(string1,'(A,1p,2(1x,E12.4))') &
              'aa < 0 and bb = 0 should have been excluded ',aa,bb
         call error_handler(E_ERR, 'unit_square_location', string1,source,revision,revdate)

      endif

   endif

endif

! If m (and maybe m_neg) is out of the possible range, return to calling program
! with found_cell still false.
if (m < 0.0_r8 .or. m > 1.0_r8) then
   if (.not.found_cell) then
      if (m_neg < 0.0_r8 .or. m_neg > 1.0_r8) then
         ! This includes m_neg == MISSING_R8, due to only m being assigned above.
         return
      endif
   ! ? Can these 2 sections ever be entered?
   else
      ! Exit with error if m is outside valid range.
      write(string1, *) 'location of ob in unit square is out of bounds m = [0,1] ',m, &
                        'but status is "found"'
      call error_handler(E_ERR, 'unit_square_location', string1, source, revision, revdate)
   endif
endif

! Use m to calculate the other unit square coordinate value, 'l'.
det = a(3,origin,cell) + a(1,origin,cell) * m
if (det /= 0.0_r8) then
   l = (x_o - a(2,origin,cell)*m) / det
else
   write(string1,'(A,I6,1X,1p4E12.4)') 'l denominator = 0: cell, angle, d, x_o, y_o',cell, angle, d, x_o, y_o
   write(string2,'(A,1X,1p4E12.4)') '  a(3) + a(1)*m = 0 : ', a(3,origin,cell), a(1,origin,cell),m
   call error_handler(E_ERR, 'unit_square_location', string1,source,revision,revdate, text2=string2)
endif

! Repeat for the -root, if it is a possibility.
if (m_neg /= MISSING_R8) then
   det = (a(3,origin,cell) + a(1,origin,cell)*m_neg)
   if (det /= 0.0_r8) then
      l_neg = (x_o -a(2,origin,cell)*m_neg) / det
   else
      write(string1,'(A,I6,1X,1p4E12.4)') 'l_neg denominator = 0: cell, angle, d, x_o, y_o', &
           cell, angle, d, x_o, y_o
      write(string2,'(A,1X,1p4E12.4)') '  a(3) + a(1)*m = 0 : ', a(3,origin,cell), a(1,origin,cell),m
      call error_handler(E_ERR, 'unit_square_location', string1,source,revision,revdate, text2=string2)
   endif

   ! Informational output, if the observation is exactly on the m-axis
!SENote: I have substituted my_task_id == 0 for the old "output_task0"; confirm
   if (l_neg == 0.0_r8 .and. my_task_id() == 0) then
      write(string1,'(A,I6,1X,1p4E12.4)') 'l_neg cell, x_o - a(2)*m = ',cell, x_o ,a(2,origin,cell),m
      call error_handler(E_MSG, 'unit_square_location', string1,source,revision,revdate)
   endif

endif

! Informational output, if the observation is exactly on the m-axis
!SENote: I have substituted my_task_id == 0 for the old "output_task0"; confirm
if (l == 0.0_r8 .and. my_task_id() == 0) then
   write(string1,'(A,I6,1X,1p4E12.4)') 'Ob is on x-axis: l-cell, x_o - a(2)*m = ',cell, x_o ,a(2,origin,cell),m
   call error_handler(E_MSG, 'unit_square_location', string1,source,revision,revdate)
endif

! If l (and maybe l_neg) is out of the possible range, return to calling program
! with found_cell still false.
if (l < 0.0_r8 .or. l > 1.0_r8) then
   if (.not.found_cell) then
      if (l_neg < 0.0_r8 .or. l_neg > 1.0_r8) then
         ! This includes m_neg == MISSING_R8, due to only m being assigned above
         ! Return with found_cell still = failure (0) to test the next cell.
         return
      endif
   ! ? Can these 2 sections ever be entered?
      ! Exit with error if l is outside valid range.
   else
      ! Exit with error if l is outside valid range.
      write(string1, *) 'location of ob in unit square is out of bounds l = [0,1] ',l, &
                        'but status is "found"'
      call error_handler(E_ERR, 'unit_square_location', string1, source, revision, revdate)
   endif
endif

! If we get this far, then this cell contains the ob.

! But which root(s) of the m quadratic equation led to the mapping?
! Put the right values in l and m.
neg_root = m_neg >= 0.0_r8 .and. m_neg <= 1.0_r8 .and. &
           l_neg >= 0.0_r8 .and. l_neg <= 1.0_r8
if (m >= 0.0_r8 .and. m <= 1.0_r8 .and. &
    l >= 0.0_r8 .and. l <= 1.0_r8 ) then
   ! Both roots yield a good mapping.
   if (neg_root) then
      write(string1, *) 'BOTH roots of the m quadratic yield usable mappings.  The +root is being used.'
      call error_handler(E_MSG, 'unit_square_location', string1, source, revision, revdate)
   endif

elseif (neg_root) then
   ! The -root yields a good mapping.  Pass along the -root m and l.
   m = m_neg
   l = l_neg
   write(string1, *) 'The negative root of the m quadratic yielded the only usable mapping.'
   call error_handler(E_MSG, 'unit_square_location', string1, source, revision, revdate)
endif

! Return with found_cell = true; success.
found_cell = .true.

end subroutine unit_square_location

!-----------------------------------------------------------------------


!SENote: THIS IS NOT A NUMERICALLY STABLE QUADRATIC SOLVER: REPLACE
subroutine solve_quadratic(a, b, c, r1, r2)

real(r8), intent(in)  :: a
real(r8), intent(in)  :: b
real(r8), intent(in)  :: c
real(r8), intent(out) :: r1
real(r8), intent(out) :: r2

real(r8) :: scaling, as, bs, cs, disc

r1 = MISSING_R8
r2 = MISSING_R8

! Scale the coefficients to get better round-off tolerance
scaling = max(abs(a), abs(b), abs(c))
as = a / scaling
bs = b / scaling
cs = c / scaling

if (abs(as) < epsilon(as)) then
   ! Solve the linear equation bs*r + cs = 0
   r1 = -cs / bs
else
   ! Get discriminant of scaled equation
   disc = bs * bs - 4.0_r8 * as * cs
   if (disc >= 0.0_r8) then

      ! Calculate the largest root (+ or - determined by sign of bs)
      ! Handling of bs = 0 different from pre-review code
      !    if(bs > 0.0_r8) then
      if(bs >= 0.0_r8) then
         r1 = (-bs - sqrt(disc)) / (2.0_r8 * as)
      else
         r1 = (-bs + sqrt(disc)) / (2.0_r8 * as)
      endif

      ! Compute the second root given the larger (not most positive) one
      if (r1 == 0.0_r8) then
         ! The b AND c must have been 0: solved the equation a*r1^2 = 0 above
         ! and there's no 2nd root.
         r2 = 0.0_r8
      else
         ! 'as' and 'r1' have been tested for 0.
         r2 = cs / (as * r1)
      endif
   endif
endif

end subroutine solve_quadratic

!------------------------------------------------------------
! Subroutines from mpas_atm/model_mod.f90, for using cartesian coordinates to
! find closest node to an ob.

subroutine init_closest_node()

! use ncol, lats and lons of nodes (corners) to initialize a get_close structure
! to be used later in find_closest_node().

integer :: i

allocate(cs_locs_xyz(ncol))

do i=1, ncol
   ! SENote: the lon and lat 1D arrays are now an element in the cam_grid type, declared grid_data
   cs_locs_xyz(i) = xyz_set_location(grid_data%lon%vals(i), grid_data%lat%vals(i), 0.0_r8, earth_radius)
enddo

! the width (2nd arg of ...init) really isn't used anymore, but it's part of the
! interface so we have to pass some number in.
call xyz_get_close_maxdist_init(cs_gc_xyz, 1.0_r8)
call xyz_get_close_obs_init    (cs_gc_xyz, ncol, cs_locs_xyz)

end subroutine init_closest_node

!------------------------------------------------------------

function find_closest_node(lat, lon)

! Determine the index for the closest node to the given point
! 2D calculation only.

real(r8), intent(in)  :: lat
real(r8), intent(in)  :: lon
integer               :: find_closest_node

type(xyz_location_type) :: pointloc
integer                 :: closest_node, rc
! This 'save' is redundant with initializing the variable here in the declaration statement.
logical, save           :: search_initialized = .false.

! do this exactly once.
if (.not. search_initialized) then
   call init_closest_node()
   search_initialized = .true.
endif

pointloc = xyz_set_location(lon, lat, 0.0_r8, earth_radius)

call xyz_find_nearest(cs_gc_xyz, pointloc, cs_locs_xyz, closest_node, rc)

! decide what to do if we don't find anything.
if (rc /= 0 .or. closest_node < 0) then
   !SENote: I have substituted my_task_id == 0 for the old "output_task0"; confirm
   if (my_task_id() == 0) then
      write(string1,*) 'cannot find a nearest node to lon, lat: ', lon, lat
      call error_handler(E_WARN, 'find_closest_node', string1,source,revision,revdate)
      ! newFIXME; should this be E_ERR instead?
   endif
   find_closest_node = -1
   return
endif

! this is the cell index for the closest center
find_closest_node = closest_node

end function find_closest_node

!-----------------------------------------------------------------------
!> internal only version of model interpolate. 
!> does not check for locations too high - return all actual values.

subroutine interpolate_se_values(state_handle, ens_size, location, obs_qty, varid, &
                              interp_vals, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_qty
integer,             intent(in) :: varid
real(r8),           intent(out) :: interp_vals(ens_size) 
integer,            intent(out) :: istatus(ens_size)

character(len=*), parameter :: routine = 'interpolate_se_values:'

integer  :: which_vert
integer :: cell_corners(4), closest, i
type(location_type) :: location_copy
!SENote : These are terrible names
real(r8) :: l, m
real(r8) :: lon_lat_vert(3), quad_vals(4, ens_size)
type(quad_interp_handle) :: interp_handle

interp_vals(:) = MISSING_R8
istatus(:)     = 99

lon_lat_vert  = get_location(location)
which_vert    = nint(query_location(location)) 

location_copy = location
call coord_ind_cs(location_copy, obs_qty, .false., closest , cell_corners, l, m)


!SENote: Now work on the vertical conversions and getting the vertical index for each ensemble member
!SENOte: switching the order on the varid and obs_qty arguments from the call isn't a great idea
call get_se_quad_vals(state_handle, ens_size, varid, obs_qty, cell_corners, &
                   lon_lat_vert, which_vert, quad_vals, istatus)

!SENOte: Again,need to sort out the error returns on all of these
if(istatus(1) /= 0) then
   istatus(:) = 3  ! cannot locate enclosing horizontal quad
   return
endif

if (any(istatus /= 0)) return


! The following uses Jeff's recommended 'generalized quadrilateral interpolation', as in
! http://www.particleincell.com/2012/quad-interpolation/.
! Most of the work is done in create_cs_grid_arrays() and coord_ind_cs().

! Interpolate from the cell's corners to the ob location on the unit square.
! This is done by weighting the field at each corner by the rectangular area
! ((l,m) space) diagonally across the ob location from the corner.
! AKA 'linear area weighting'.

! SENote: The status block needs to be worked on to be consistent with CLASSIC
!SENote: It is commented out for now but should be put back in.
! The vals indices are consistent with how mapping of corners was done,
! and how cell_corners was assigned.
!SENoteif (vstatus == 1) then
   !SENoteif (print_details) then
      !SENotewrite(string1,'(A,2F10.6,1pE20.12)') 'istatus = 1, no interpolation'
      !SENotecall error_handler(E_MSG, 'interp_cubed_sphere', string1)
   !SENoteendif
   !SENotereturn
!SENoteelse
   !SENoteif (abs(lon_lat_lev(2)) > max_obs_lat_degree) then
      !SENote! Define istatus to be a combination of vstatus (= 0 or 2 (for higher than highest_obs...))
      !SENote! and whether the ob is poleward of the limits set in the namelist (+ 4).
      !SENoteistatus = 10*vstatus + 4
   !SENoteelse
      !SENoteistatus = vstatus
   !SENoteendif
!SENoteendif

! SENote: could this be done in vector notation?
do i = 1, ens_size
interp_vals(i) = quad_vals(2, i) *           l *          m &
           + quad_vals(1, i) * (1.0_r8 - l)*          m &
           + quad_vals(4, i) * (1.0_r8 - l)*(1.0_r8 - m) &
           + quad_vals(3, i) *           l *(1.0_r8 - m)
end do

!SENote Probably  not right
if (any(istatus /= 0)) then
   istatus(:) = 8   ! cannot evaluate in the quad
   return
endif

end subroutine interpolate_se_values

!-----------------------------------------------------------------------
!>

subroutine get_se_quad_vals(state_handle, ens_size, varid, obs_qty, corners, &
                         lon_lat_vert, which_vert, quad_vals, my_status)
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: varid
integer,             intent(in) :: obs_qty
integer,             intent(in) :: corners(4)
real(r8),            intent(in) :: lon_lat_vert(3)
integer,             intent(in) :: which_vert
real(r8),           intent(out) :: quad_vals(4, ens_size) !< array of interpolated values
integer,            intent(out) :: my_status(ens_size)

integer  :: icorner, numdims
integer  :: level_one_array(ens_size)
integer  :: four_levs1(4, ens_size), four_levs2(4, ens_size)
real(r8) :: four_vert_fracts(4, ens_size)
!SENote: Need unused integer(4) storage for call to get the four state values
integer :: lats_unused_array(4)

!SENote: four_lons and four_lats just around for legacy during development
integer :: four_lons(4), four_lats(4)

character(len=*), parameter :: routine = 'get_se_quad_vals:'

!SENote Initialize the latitude arrays to 1 since they are not used
lats_unused_array = 1
quad_vals(:,:) = MISSING_R8
my_status(:) = 99

! need to consider the case for 2d vs 3d variables
numdims = get_dims_from_qty(obs_qty, varid)

!SENote: The dimensions are one less than for the FV. Look for ways to share code
! now here potentially we have different results for different
! ensemble members.  the things that can vary are dimensioned by ens_size.
!SENote: It's the first 2 indices (or 1) that are actually in use when fetching data?

!SENote
write(*, *) 'in get_se_quad_vals varid, obs_qty, numdims', varid, obs_qty, numdims

if (numdims == 2) then

   ! build 4 columns to find vertical level numbers
   do icorner=1, 4
      call find_se_vertical_levels(state_handle, ens_size, &
                                corners(icorner), lon_lat_vert(3), &
                                which_vert, obs_qty, varid, &
                                four_levs1(icorner, :), four_levs2(icorner, :), & 
                                four_vert_fracts(icorner, :), my_status)

      if (any(my_status /= 0)) return
      write(*, *) 'finding vertical levels in get_se_quad_vals ', &
         icorner, four_levs1(icorner, :), four_levs2(icorner, :), four_vert_fracts(icorner, :)
  
   enddo
   
   ! we have all the indices and fractions we could ever want.
   ! now get the data values at the bottom levels, the top levels, 
   ! and do vertical interpolation to get the 4 values in the columns.
   ! the final horizontal interpolation will happen later.
      
   if (varid > 0) then

      call get_se_four_state_values(state_handle, ens_size, corners, &
                                four_levs1, four_levs2, four_vert_fracts, &   
                                varid, quad_vals, my_status)

   else ! get 2d special variables in another ways ( like QTY_PRESSURE )
      call get_se_four_nonstate_values(state_handle, ens_size, corners, &
                                   four_levs1, four_levs2, four_vert_fracts, & 
                                   obs_qty, quad_vals, my_status)

   endif

   if (any(my_status /= 0)) return

else if (numdims == 1) then

   if (varid > 0) then
      level_one_array(:) = 1
      do icorner=1, 4
         !SENote: just passing in corner indices
         call get_se_values_from_varid(state_handle,  ens_size, corners(icorner), & 
                                    level_one_array, varid, quad_vals(icorner,:),my_status)

         if (any(my_status /= 0)) return

      enddo

   else ! special 2d case
      do icorner=1, 4
         call get_se_quad_values(ens_size, four_lons(icorner), obs_qty, quad_vals(icorner,:))
      enddo
      ! apparently this can't fail
      my_status(:) = 0
      
   endif

else
   write(string1, *) trim(get_name_for_quantity(obs_qty)), ' has dimension ', numdims
   call error_handler(E_ERR, routine, 'only supports 1D or 2D fields', &
                      source, revision, revdate, text2=string1)
endif

! when you get here, my_status() was set either by passing it to a
! subroutine, or setting it explicitly here.  if this routine returns
! the default value of 99 something went wrong in this logic.

end subroutine get_se_quad_vals

!-----------------------------------------------------------------------
!>

subroutine get_se_four_state_values(state_handle, ens_size, four_corners, &
                                 four_levs1, four_levs2, four_vert_fracts, &
                                 varid, quad_vals, my_status)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: four_corners(4)
integer,             intent(in) :: four_levs1(4, ens_size), four_levs2(4, ens_size)
real(r8),            intent(in) :: four_vert_fracts(4, ens_size)
integer,             intent(in) :: varid
real(r8),           intent(out) :: quad_vals(4, ens_size) !< array of interpolated values
integer,            intent(out) :: my_status(ens_size)

integer  :: icorner
real(r8) :: vals1(ens_size), vals2(ens_size)

character(len=*), parameter :: routine = 'get_se_four_state_values:'

do icorner=1, 4
   call get_se_values_from_varid(state_handle,  ens_size, four_corners(icorner), &
      four_levs1(icorner, :), varid, vals1, my_status)

   if (any(my_status /= 0)) then
      my_status(:) = 16   ! cannot retrieve vals1 values
      return
   endif

   call get_se_values_from_varid(state_handle,  ens_size, four_corners(icorner), &
      four_levs2(icorner, :), varid, vals2, my_status)
   if (any(my_status /= 0)) then
      my_status(:) = 17   ! cannot retrieve top values
      return
   endif

   call vert_interp(ens_size, vals1, vals2, four_vert_fracts(icorner, :), & 
                    quad_vals(icorner, :))

enddo


end subroutine get_se_four_state_values

!-----------------------------------------------------------------------
!>

subroutine get_se_four_nonstate_values(state_handle, ens_size, four_corners, &
                                 four_levs1, four_levs2, four_vert_fracts, &
                                 obs_qty, quad_vals, my_status)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: four_corners(4)
integer,             intent(in) :: four_levs1(4, ens_size), four_levs2(4, ens_size)
real(r8),            intent(in) :: four_vert_fracts(4, ens_size)
integer,             intent(in) :: obs_qty
real(r8),           intent(out) :: quad_vals(4, ens_size) !< array of interpolated values
integer,            intent(out) :: my_status(ens_size)

integer  :: icorner
real(r8) :: vals1(ens_size), vals2(ens_size)

character(len=*), parameter :: routine = 'get_se_four_nonstate_values:'

do icorner=1, 4
   call get_se_values_from_nonstate_fields(state_handle,  ens_size, four_corners(icorner), &
                              four_levs1(icorner, :), obs_qty, vals1, my_status)
   if (any(my_status /= 0)) then
      my_status(:) = 16   ! cannot retrieve vals1 values
      return
   endif

   call get_se_values_from_nonstate_fields(state_handle,  ens_size, four_corners(icorner), &
                              four_levs2(icorner, :), obs_qty, vals2, my_status)
   if (any(my_status /= 0)) then
      my_status(:) = 17   ! cannot retrieve top values
      return
   endif

   call vert_interp(ens_size, vals1, vals2, four_vert_fracts(icorner, :), &
                    quad_vals(icorner, :))

enddo

end subroutine get_se_four_nonstate_values


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
         !SENote: in SE this is a 1 dimensional field, but this is okay
         get_dims_from_qty = 1
      case (QTY_PRESSURE, QTY_GEOMETRIC_HEIGHT)
         !SENote: in SE these are 2d fields
         get_dims_from_qty = 2
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

!SENote: THIS HAS NOT YET BEEN IMPLEMENTED
subroutine get_se_quad_values(ens_size, corner_index, obs_quantity, vals)
integer,  intent(in) :: ens_size
integer,  intent(in) :: corner_index
integer,  intent(in) :: obs_quantity
real(r8), intent(out) :: vals(ens_size) 

character(len=*), parameter :: routine = 'get_se_quad_values'

integer :: prev_lon, next_lat
real(r8) :: vals1(ens_size), vals2(ens_size)


select case (obs_quantity)
   case (QTY_SURFACE_ELEVATION)
      !SENote: Just return phis for this corner_index
      vals = phis(corner_index, 1)

     !>@todo FIXME:
     ! should this be using gravity at the given latitude? 
     vals = vals / gravity

   case default 
      write(string1, *) 'we can not interpolate qty', obs_quantity
      call error_handler(E_ERR,routine,string1,source,revision,revdate)

end select

end subroutine get_se_quad_values

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
!SENote: Changed from original FV
!> given a corner index number, a quantity and a vertical value and type,
!> return which two levels these are between and the fraction across.
!> 

subroutine find_se_vertical_levels(ens_handle, ens_size, corner_index, vert_val, &
                                which_vert, obs_qty, var_id, levs1, levs2, vert_fracts, my_status)
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: corner_index 
real(r8),            intent(in)  :: vert_val
integer,             intent(in)  :: which_vert
integer,             intent(in)  :: obs_qty
integer,             intent(in)  :: var_id
integer,             intent(out) :: levs1(ens_size)
integer,             intent(out) :: levs2(ens_size)
real(r8),            intent(out) :: vert_fracts(ens_size)
integer,             intent(out) :: my_status(ens_size)

character(len=*), parameter :: routine = 'find_se_vertical_levels:'

integer  :: l1, l2, imember, level_one, status1, k
real(r8) :: fract1
real(r8) :: surf_pressure (  ens_size )
real(r8) :: pressure_array( ref_nlevels, ens_size )
real(r8) :: height_array  ( ref_nlevels, ens_size )

!SENote: Temps to replace old arguments
integer :: lon_index, lat_index

!SENote; Debug temps
integer :: i, j

! assume the worst
levs1(:)    = MISSING_I
levs2(:)    = MISSING_I
vert_fracts(:) = MISSING_R8
my_status(:)   = 98


! ref_nlevels is the number of vertical levels (midlayer points)

level_one = 1

select case (which_vert)

   !SENote: Classic has an option for scale_height. Not found here. No observations are reported in scale_height?
   ! Need to see how this works with localization.
   case(VERTISPRESSURE)
      ! construct a pressure column here and find the model levels
      ! that enclose this value
      !SENote: No staggering, so just get the surface pressure at this column
      call get_se_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, corner_index, level_one, &
         surf_pressure, status1)

      !SENote: Need to make a careful study of the error value propagation throughout this call tree
      if (status1 /= 0) then
         my_status(:) = status1
         return
      endif

      call build_cam_pressure_columns(ens_size, surf_pressure, ref_nlevels, pressure_array)
      !SENote: Need to make sure that there are no residual vertical stagger issues here

      do imember = 1, ens_size
         call pressure_to_level(ref_nlevels, pressure_array(:, imember), vert_val, & 
                                levs1(imember), levs2(imember), &
                                vert_fracts(imember), my_status(imember))
      enddo

      if (debug_level > 100) then
         do k = 1,ens_size
            print*, 'ISPRESSURE levs1(k), levs2(k), vert_fracts(k), vert_val', &
                     levs1(k), levs2(k), vert_fracts(k), vert_val, pressure_array(levs1(k) , k), pressure_array(levs2(k), k) 
          enddo
      endif

   case(VERTISHEIGHT)
      ! construct a height column here and find the model levels
      ! that enclose this value
      call cam_se_height_levels(ens_handle, ens_size, corner_index, ref_nlevels, &
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
!SENote, not implemented? It is not clear that anything ever comes through this code block. Need to understand.
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

end subroutine find_se_vertical_levels


!-----------------------------------------------------------------------
!> Compute the heights at pressure midpoints
!>
!> this version does all ensemble members at once.

subroutine cam_se_height_levels(ens_handle, ens_size, column_index, nlevels, height_array, my_status) 
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: column_index
integer,             intent(in)  :: nlevels
real(r8),            intent(out) :: height_array(nlevels, ens_size)
integer,             intent(out) :: my_status(ens_size)

integer  :: k, level_one, imember, status1
real(r8) :: surface_elevation(1)
real(r8) :: surface_pressure(ens_size), mbar(nlevels, ens_size)
real(r8) :: tv(nlevels, ens_size)  ! Virtual temperature, top to bottom

! this is for surface obs
level_one = 1

! Get the surface pressure at this column
call get_se_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, column_index, level_one, &
   surface_pressure, status1)

! get the surface elevation from the phis, including stagger if needed
call get_se_quad_values(1, column_index, QTY_SURFACE_ELEVATION, surface_elevation)

call compute_se_virtual_temperature(ens_handle, ens_size, column_index, nlevels, tv, status1)

if (status1 /= 0) then
   my_status = status1
   return
endif

!SENote: need to change this call CHECK THIS CAREFULLY
if (use_variable_mean_mass) then
   call compute_se_mean_mass(ens_handle, ens_size, column_index, nlevels, mbar, status1)
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
call gph2gmh(height_array, grid_data%lat%vals(column_index))

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

end subroutine cam_se_height_levels

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

! SENote: No phis available in SE for now
!SENote not available, deallocate(phis)

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

!SENote: Added temp variable to get ncol
integer :: ncol_temp(1)


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

! SENote: This seems like a good place to load up the module storage variable ncol
! Try to get it from state structure in a nice way
! The size of the only surface pressure dimension is the number of columns
ncol_temp = get_dim_lengths(domain_id,  get_varid_from_kind(domain_id, QTY_SURFACE_PRESSURE))
ncol = ncol_temp(1)

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

!SENote; We will need to do initialization of interpolation for SE, but not with this
! Set up the interpolation structures for later 
!SENotecall setup_interpolation(grid)

end subroutine read_grid_info

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

!SENote: phis array is really only 1D in the file, but keep a second dimension of 1 for shared code
allocate( phis(nsize(1), 1) )

call nc_get_variable(ncid, 'PHIS', phis, routine)

call nc_close_file(ncid, routine)

end subroutine read_cam_phis_array


!-----------------------------------------------------------------------
!> Compute the virtual temperature at the midpoints
!>
!> this version does all ensemble members at once.
!>

subroutine compute_se_virtual_temperature(ens_handle, ens_size, column_index, nlevels, tv, istatus)

type(ensemble_type), intent(in)   :: ens_handle
integer,             intent(in)   :: ens_size
integer,             intent(in)   :: column_index
integer,             intent(in)   :: nlevels
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
   call get_se_values_from_single_level(ens_handle, ens_size, QTY_TEMPERATURE, column_index, k, &
      temperature, istatus)

   if (istatus < 0) return

   ! specific humidity
   call get_se_values_from_single_level(ens_handle, ens_size, QTY_SPECIFIC_HUMIDITY, column_index, k, &
      specific_humidity, istatus)
   
   if (istatus < 0) return

   !>tv == virtual temperature.
   tv(k,:) = temperature(:)*(1.0_r8 + rr_factor*specific_humidity(:))
   !print*, 'tv(levels)', k,tv(k,1), temperature(1), specific_humidity(1)
enddo

end subroutine compute_se_virtual_temperature


!-----------------------------------------------------------------------
!> loop through all levels to get the mean mass.
!>


!SENote: No test available? Need to work with Nick at some point to test WACCM-X configurations.
subroutine compute_se_mean_mass(ens_handle, ens_size, col_index, nlevels, mbar, istatus)
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: col_index
integer,             intent(in)  :: nlevels
real(r8),            intent(out) :: mbar(nlevels, ens_size)
integer,             intent(out) :: istatus

integer :: k, this_qty
real(r8) :: mmr_o1(ens_size, nlevels), &
            mmr_o2(ens_size, nlevels), &
            mmr_h1(ens_size, nlevels), &
            mmr_n2(ens_size, nlevels)
real(r8) :: O_molar_mass, O2_molar_mass, H_molar_mass, N2_molar_mass 

character(len=*), parameter :: routine = 'compute_se_mean_mass'

!SENote: This subroutine has not been tested yet for SE. Need to work with WACCM folks to test.
call error_handler(E_ERR, routine, 'Subroutine has not been tested', source, revision, revdate)


! do this outside the subroutine?  it never changes throughout the
! run of the program
!SENote: could do an initialization with save storage
O_molar_mass  = get_molar_mass(QTY_ATOMIC_OXYGEN_MIXING_RATIO)
O2_molar_mass = get_molar_mass(QTY_MOLEC_OXYGEN_MIXING_RATIO)
H_molar_mass  = get_molar_mass(QTY_ATOMIC_H_MIXING_RATIO)
N2_molar_mass = get_molar_mass(QTY_NITROGEN)
   


! High topped models (WACCM-X) need to account for the changing composition 
! of the atmosphere with height.  This requires several variables from the
! initial file, which may not be available from low topped models.
do k = 1, nlevels

   call get_se_values_from_single_level(ens_handle, ens_size, QTY_ATOMIC_OXYGEN_MIXING_RATIO, &
      col_index, k, mmr_o1(:, k), istatus)
   if (istatus /= 0) return
   !print *, 'mmr: ', trim(get_name_for_quantity(this_qty)), mmr_o1(1, k)
   
   call get_se_values_from_single_level(ens_handle, ens_size, QTY_MOLEC_OXYGEN_MIXING_RATIO, &
      col_index, k, mmr_o2(:, k), istatus)
   if (istatus /= 0) return
   !print *, 'mmr: ', trim(get_name_for_quantity(this_qty)), mmr_o2(1, k)
   
   call get_se_values_from_single_level(ens_handle, ens_size, QTY_ATOMIC_H_MIXING_RATIO, &
      col_index, k, mmr_h1(:, k), istatus)
   if (istatus /= 0) return
   !print *, 'mmr: ', trim(get_name_for_quantity(this_qty)), mmr_h1(1, k)
   
   mmr_n2(:,k) = 1.0_r8 - (mmr_o1(:,k) + mmr_o2(:,k) + mmr_h1(:,k))
   mbar(k,:) = 1.0_r8/( mmr_o1(:,k)/O_molar_mass  &
                      + mmr_o2(:,k)/O2_molar_mass &
                      + mmr_h1(:,k)/H_molar_mass  &
                      + mmr_n2(:,k)/N2_molar_mass)
enddo

end subroutine compute_se_mean_mass

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
!SENote: This argument is not used here, but is required to support external calls. Should reexamine
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
         call state_vertical_to_pressure(    ens_handle, ens_size, locs(i), loc_indx(i))
      case (VERTISHEIGHT)
         call state_vertical_to_height(      ens_handle, ens_size, locs(i), loc_indx(i))
      case (VERTISLEVEL)
         call state_vertical_to_level(                   ens_size, locs(i), loc_indx(i))
      case (VERTISSCALEHEIGHT)
         call state_vertical_to_scaleheight( ens_handle, ens_size, locs(i), loc_indx(i))
      case default
         write(string1,*)'unable to convert vertical state "', which_vert, '"'
         call error_handler(E_MSG,routine,string1,source,revision,revdate)
   end select
enddo

istatus = 0

end subroutine convert_vertical_state

!--------------------------------------------------------------------

subroutine state_vertical_to_pressure(ens_handle, ens_size, location, location_indx)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx

integer  :: column, level, no_third_dimension, myqty, level_one, status1
integer  :: my_status(ens_size)
real(r8) :: pressure_array(ref_nlevels), surface_pressure(ens_size)


!SENote
!Need to clean up the index naming here. See obs_vertical_to_pressure
!Need to do the other state_vertical conversion routines correctly
call get_model_variable_indices(location_indx, column, level, no_third_dimension, kind_index=myqty)

if (is_surface_field(myqty)) then
   
   level_one = 1
   call get_se_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                                     column, level_one, surface_pressure, status1)

   if (status1 /= 0) then
      return
   endif
   call set_vertical(location, surface_pressure(1), VERTISPRESSURE)
else
   call cam_se_pressure_levels(ens_handle, ens_size, column, ref_nlevels, &
                            pressure_array, my_status)

   call set_vertical(location, pressure_array(level), VERTISPRESSURE)
endif

end subroutine state_vertical_to_pressure

!--------------------------------------------------------------------

subroutine state_vertical_to_height(ens_handle, ens_size, location, location_indx)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx

integer  :: column, level, no_third_dimension, my_status(ens_size)
real(r8) :: height_array(ref_nlevels, ens_size)

! build a height column and a pressure column and find the levels
call get_model_variable_indices(location_indx, column, level, no_third_dimension)

call cam_se_height_levels(ens_handle, ens_size, column, ref_nlevels, &
                       height_array, my_status) 

!>@todo FIXME this can only be used if ensemble size is 1
!SENote: Confirm this
call set_vertical(location, height_array(level, 1), VERTISHEIGHT)

end subroutine state_vertical_to_height

!--------------------------------------------------------------------

subroutine state_vertical_to_scaleheight(ens_handle, ens_size, location, location_indx)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx

integer  :: column, level, no_third_dimension, level_one, status1, my_status(ens_size)
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
      call get_model_variable_indices(location_indx, column, level, no_third_dimension)

      call get_se_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                                        column, level, surface_pressure, status1)
      if (status1 /= 0) goto 200

      scaleheight_val = log(surface_pressure(1))

   else

      ! build a pressure column and and find the levels
      call get_model_variable_indices(location_indx, column, level, no_third_dimension)

      call cam_se_pressure_levels(ens_handle, ens_size, column,  ref_nlevels, &
                               pressure_array, my_status)
      if (any(my_status /= 0)) goto 200
   
      scaleheight_val = log(pressure_array(level))

   endif

else

   ! handle surface obs separately here.
   if (query_location(location) == VERTISSURFACE) then

      scaleheight_val = 0.0_r8   ! log(1.0)

   else

      ! build a pressure column and and find the levels
      call get_model_variable_indices(location_indx, column, level, no_third_dimension)

      call cam_se_pressure_levels(ens_handle, ens_size, column, ref_nlevels, &
                               pressure_array, my_status)
      if (any(my_status /= 0)) goto 200

      ! get the surface pressure from the ens_handle
      call get_se_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                                        column, level, surface_pressure, status1)
      if (status1 /= 0) goto 200
   
      scaleheight_val = scale_height(pressure_array(level), surface_pressure(1), no_normalization_of_scale_heights)

   endif

endif
   
200 continue   ! done

call set_vertical(location, scaleheight_val, VERTISSCALEHEIGHT)

end subroutine state_vertical_to_scaleheight

!--------------------------------------------------------------------

subroutine state_vertical_to_level(ens_size, location, location_indx)
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx

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

subroutine cam_se_pressure_levels(ens_handle, ens_size, col_index, nlevels, &
                               pressure_array, my_status) 
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: col_index
integer,             intent(in)  :: nlevels
real(r8),            intent(out) :: pressure_array(nlevels, ens_size)
integer,             intent(out) :: my_status(ens_size)

integer     :: level_one, status1
real(r8)    :: surface_pressure(ens_size)

! this is for surface obs
level_one = 1

! get the surface pressure from the ens_handle
      call get_se_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, col_index, level_one, &
         surface_pressure, status1)

if (status1 /= 0) then
   my_status(:) = status1
   return
endif

call build_cam_pressure_columns(ens_size, surface_pressure, ref_nlevels, pressure_array)

! No error returns available if we get here, so all good
my_status(:) = 0

end subroutine cam_se_pressure_levels

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

call ok_to_interpolate(qty, varid, domain_id, my_status)
if (my_status /= 0) return

!SENote
write(*, *) 'in obs_vertical_to_pressure'
call interpolate_se_values(ens_handle, ens_size, location, &
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

!SENote Does this work for the FV? 
! Doesn't actually appear to work right for the FV which just blasts through a failed search 
! for the height field. This could be fixed. Confirm with Kevin.
ens_size = 1

call ok_to_interpolate(QTY_GEOMETRIC_HEIGHT, varid, domain_id, my_status)
if (my_status /= 0) return

call interpolate_se_values(ens_handle, ens_size, location, &
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
!SENote
write(*, *) 'In obs_vertical_to_level'

!SENote: This has not been checked, just swapped the call to interpolate_values
call interpolate_se_values(ens_handle, ens_size, location, &
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

   call ok_to_interpolate(ptype, varid1, domain_id, my_status)
   if (my_status /= 0) return
      
   !>@todo FIXME IFF the obs location is already pressure, we can take it at
   !> face value here and not interpolate it.  however it won't fail if the
   !> pressure here is less than the ensemble mean pressure at this point.
   !> is that ok?
   
   if (ptype == QTY_PRESSURE .and. is_vertical(location, "PRESSURE")) then
      pressure_array(:) = query_location(location, "VLOC")
      my_status = 0
   else
      call interpolate_se_values(ens_handle, ens_size, location, ptype, varid1, &
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

      call ok_to_interpolate(QTY_PRESSURE, varid1, domain_id, my_status)
      if (my_status /= 0) return
   
      !>@todo FIXME IFF the obs location is already pressure, we can take it at
      !> face value here and not interpolate it.  however, it can result in negative
      !> scale height values if the pressure is larger than the surface pressure at
      !> that location.  that's what the original cam model_mod did.  is that ok?
   
      if (ptype == QTY_PRESSURE .and. is_vertical(location, "PRESSURE")) then
         pressure_array(:) = query_location(location, "VLOC")
         my_status = 0
      else
         call interpolate_se_values(ens_handle, ens_size, location, QTY_PRESSURE, varid1, &
                                    pressure_array(:), status(:))
         if (status(1) /= 0) then
            my_status = status(1)
            return
         endif
      endif
                                 
      call ok_to_interpolate(QTY_SURFACE_PRESSURE, varid2, domain_id, my_status)
      if (my_status /= 0) return
      
      call interpolate_se_values(ens_handle, ens_size, location, QTY_SURFACE_PRESSURE, varid2, &
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

!SENote
write(*, *) 'in get_close_state before convert_vert_one_obs'
call write_location(0, base_loc)
if (vert_type /= vertical_localization_type) then
   call convert_vert_one_obs(ens_handle, base_loc, base_type, &
                             vertical_localization_type, status)
!SENote
write(*, *) 'back form convert_vert_one_obs status ', status
stop
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

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
