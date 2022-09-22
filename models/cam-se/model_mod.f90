! DART software - Copyright UCAR. This open source software is provided
! by ucar, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/dares/dart/dart_download
!
!----------------------------------------------------------------
!>
!> this is the interface between the cam-se atmosphere model and dart.
!> the required public interfaces and arguments cannot be changed.
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

use         utilities_mod,  only : find_namelist_in_file, check_namelist_read, &
                                   string_to_logical, string_to_real,& 
                                   nmlfileunit, do_nml_file, do_nml_term, &
                                   register_module, error_handler, &
                                   file_exist, to_upper, E_ERR, E_MSG, E_WARN, array_dump
use          obs_kind_mod,  only : QTY_SURFACE_ELEVATION, QTY_PRESSURE, &
                                   QTY_GEOMETRIC_HEIGHT, QTY_VERTLEVEL, &
                                   QTY_SURFACE_PRESSURE, &
                                   QTY_TEMPERATURE, QTY_SPECIFIC_HUMIDITY, &
                                   QTY_MOLEC_OXYGEN_MIXING_RATIO, &
                                   QTY_ION_O_MIXING_RATIO, QTY_ATOMIC_H_MIXING_RATIO, &
                                   QTY_ATOMIC_OXYGEN_MIXING_RATIO, QTY_NITROGEN, &
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
use  netcdf_utilities_mod,  only : nc_get_variable, nc_get_variable_size, nc_create_file, &
                                   nc_add_attribute_to_variable, &
                                   nc_define_integer_variable, nc_define_double_variable, &
                                   nc_define_real_variable, &
                                   nc_define_real_scalar, &
                                   nc_add_global_creation_time, &
                                   nc_add_global_attribute, &
                                   nc_define_dimension, nc_put_variable, &
                                   nc_synchronize_file, nc_end_define_mode, &
                                   nc_begin_define_mode, nc_open_file_readonly, &
                                   nc_close_file, nc_variable_exists, nc_get_global_attribute, &
                                   nc_get_dimension_size
use        chem_tables_mod, only : init_chem_tables, finalize_chem_tables, &
                                   get_molar_mass, get_volume_mixing_ratio
use     default_model_mod,  only : adv_1step, nc_write_model_vars, &
                                   init_time => fail_init_time,    &
                                   init_conditions => fail_init_conditions

use    cam_common_code_mod, only : above_ramp_start, are_damping, build_cam_pressure_columns, build_heights, &
                                   cam_grid, cdebug_level, check_good_levels, cno_normalization_of_scale_heights, &
                                   pert_model_copies, cuse_log_vertical_scale, discarding_high_obs, &
                                   free_cam_grid, free_std_atm_tables, generic_height_to_pressure,  &
                                   gph2gmh, height_to_level, init_damping_ramp_info, &
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

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'cam-se/model_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

! model_nml namelist variables and default values
! Which vertical coordinate: Dry mass if for versions with CESM2 and later
logical            :: dry_mass_vertical_coordinate    = .true.
! If false, uses less precise but vastly cheaper traditional hybrid vertical coordinate for get_close
logical            :: precise_dry_mass_get_close      = .false.

character(len=256) :: cam_template_filename           = 'caminput.nc'
character(len=256) :: cam_phis_filename               = 'cam_phis.nc'

! Identify the CS grid mapping files
character(len=256) :: homme_map_file                  = 'SEMapping.nc'            ! Corners of each cubed sphere cell.
character(len=256) :: cs_grid_file                    = 'SEMapping_cs_grid.nc'    ! Relationships among corners/nodes.

character(len=32)  :: vertical_localization_coord     = 'PRESSURE'
logical            :: use_log_vertical_scale          = .false.
integer            :: assimilation_period_days        = 0
integer            :: assimilation_period_seconds     = 21600
integer            :: no_obs_assim_above_level        = -1      ! model levels
integer            :: model_damping_ends_at_level     = -1      ! model levels
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
   dry_mass_vertical_coordinate,        &
   precise_dry_mass_get_close,          &
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

! Surface potential; used for calculation of geometric heights.
! SENote: right now every process has their own complete copy of this
real(r8), allocatable :: phis(:)

logical :: l_refined = .false.       ! Flag to tell whether grid is a refined mesh or not.

! A veriety of module storage data structures for geometry of grid
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CS Variables holding relationships among cubed sphere nodes.

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

! Used for finding horizontal bounding grid cells
type(get_close_type) :: cs_gc

! Array of KINDs of cubed sphere grid points.
! As of 2014-3-28 this is only used by location_mod, which doesn't actually use it.
integer, allocatable :: cs_kinds(:)

! Other useful 1D grid arrays (for cubed sphere)
real(r8), allocatable :: lon_rad(:), lat_rad(:)   ! longitude and latitude in radians, used by bearings()

! This integer is a reminder that some shared calls take 3 dimensions, but SE has only 2: Value is irrelevant
integer :: no_third_dimension = -99

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!SENote: This variable gives extra output for the locations. More global way to do this?
! set to .true. to get more details about the state vector and the
! CAM fields and sizes in the init code.
logical :: print_details = .true.

contains

!-----------------------------------------------------------------------
! All the public interfaces are first.
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

integer :: iunit, io, i 
integer :: nc_file_ID
integer :: ncol_temp(1)

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
cperturbation_amplitude = perturbation_amplitude
cassimilation_period_days = assimilation_period_days
cassimilation_period_seconds = assimilation_period_seconds
csuppress_grid_info_in_output = suppress_grid_info_in_output

call set_calendar_type('GREGORIAN')

call read_grid_info(cam_template_filename)
! This non-state variable is used to compute surface elevation.
call read_cam_phis_array(cam_phis_filename)

! initialize global values that are used frequently
call init_globals()

! read the namelist &model_nml :: state_variables
! to set up what will be read into the cam state vector
call set_cam_variable_info(cam_template_filename, state_variables)

! The size of the only surface pressure dimension is the number of columns
ncol_temp = get_dim_lengths(domain_id,  get_varid_from_kind(domain_id, QTY_SURFACE_PRESSURE))
ncol = ncol_temp(1)

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


!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! This section reads and/or builds the tables needed for model_interpolate for SE
! Read in or create a file containing the relationships among cubed sphere nodes,
! such as neighbors, centers, and bearings, which will be used to identify the cell
! which contains an observation.
! Fields will be stored in global storage.
! Write the cubed sphere grid arrays to a new NetCDF file.

! Fill arrays that are useful for bearings and distances.
allocate(lon_rad(ncol), lat_rad(ncol))
lon_rad(:) = grid_data%lon%vals(:)*DEG2RAD
lat_rad(:) = grid_data%lat%vals(:)*DEG2RAD

! Following block is also used for search for close corners
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
call fill_gc()

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

integer  :: column, level
integer  :: myvarid, myqty, nd

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, column, level, no_third_dimension, var_id=myvarid, kind_index=myqty)

nd = get_num_dims(domain_id, myvarid)

location = get_location_from_index(column, level, myqty, nd)


! return state quantity for this index if requested
if (present(var_type)) var_type = myqty

end subroutine get_state_meta_data

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
! Many of these error status returns cannot actually happen. Need to verify these.
!> istatus = 2    asked to interpolate an unknown/unsupported quantity
!> istatus = 8    cannot interpolate level, out of range
!> istatus = 10   cannot interpolate in pressure
!> istatus = 11   cannot interpolate in height
!> istatus = 12   cannot get values from obs quantity
!> istatus = 14   obs above user-defined assimilation top pressure
!> istatus = 16   cannot do vertical interpolation for bottom layer
!> istatus = 17   cannot do vertical interpolation for top layer
!>

subroutine model_interpolate(state_handle, ens_size, location, obs_qty, interp_vals, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_qty
real(r8),           intent(out) :: interp_vals(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

character(len=*), parameter :: routine = 'model_interpolate:'

! Should dry mass vertical coordinate be used to build vertical pressure columns? 
! This is expensive but probably necessary to get unbiased forward operators
logical :: precise = .true.

integer  :: varid, which_vert, status1
real(r8) :: lon_lat_vert(3)
real(r8) :: quad_vals(4, ens_size)
integer :: cell_corners(4)
real(r8) :: l_weight, m_weight
type(location_type) :: location_copy

if ( .not. module_initialized ) call static_init_model

! Successful istatus is 0
interp_vals(:) = MISSING_R8
istatus(:)     = 99

! do we know how to interpolate this quantity? Returns status1 = 0 if OK, status1 = 2 if not OK.
call ok_to_interpolate(obs_qty, varid, status1)

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
   ! Returns status1 = 0 if OK, status1 = 14 if too high.
   call obs_too_high(lon_lat_vert(3), which_vert, status1)
   if (status1 /= 0) then
      istatus(:) = status1
      return
   endif
endif


! Do the interpolation here
! First step, find the columns of the four 'corners' containing the location
! SENote2: In the CLASSIC, there is a possibility that the cell_corner was already found and this call can
! be skipped. Understand that and implement as needed.
! Note that cannot pass location directly because it is intent(inout) in coord_ind_cs.
location_copy = location
call coord_ind_cs(location_copy, obs_qty, cell_corners, l_weight, m_weight)

! Now do vertical conversions and get the vertical index for each ensemble member
call get_se_quad_vals(state_handle, ens_size, varid, obs_qty, cell_corners, &
                   lon_lat_vert, which_vert, precise, quad_vals, istatus)

!SENote Do further study of how we want to return istatus for various failures
! For now return istatus 12 for any of the failure modes
if (any(istatus /= 0)) then
   istatus = 12
   return
endif


! Then interpolate horizontally to the (lon,lat) of the ob.
! The following uses JLA's recommended 'generalized quadrilateral interpolation', as in
! http://www.particleincell.com/2012/quad-interpolation/.
! Most of the work is done in create_cs_grid_arrays() and coord_ind_cs().

! Interpolate from the cell's corners to the ob location on the unit square.
! This is done by weighting the field at each corner by the rectangular area
! ((l,m) space) diagonally across the ob location from the corner.
! AKA 'linear area weighting'.

interp_vals(:) = quad_vals(2, :) *       l_weight *          m_weight &
           + quad_vals(1, :) * (1.0_r8 - l_weight)*          m_weight &
           + quad_vals(4, :) * (1.0_r8 - l_weight)*(1.0_r8 - m_weight) &
           + quad_vals(3, :) *           l_weight *(1.0_r8 - m_weight)

if (using_chemistry) &
   interp_vals = interp_vals * get_volume_mixing_ratio(obs_qty)

! all interp values should be set by now.  set istatus
istatus(:) = 0

end subroutine model_interpolate

!-----------------------------------------------------------------------
!>
!> Does any shutdown and clean-up needed for model.
!>

subroutine end_model()

! deallocate arrays from grid and anything else

call free_cam_grid(grid_data)

deallocate(phis)

call free_std_atm_tables()

if (using_chemistry) call finalize_chem_tables()

end subroutine end_model


!--------------------------------------------------------------------
!> Does an conversion to localization vertical coordinate for a set of obs
!> Returns my_status 2 in not able to interp this quantity, 3 if get_se_quad_vals fails.

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
         call obs_vertical_to_pressure(   ens_handle, precise_dry_mass_get_close, locs(i), my_status(i))
      case (VERTISHEIGHT)
         call obs_vertical_to_height(     ens_handle, precise_dry_mass_get_close, locs(i), my_status(i))
      case (VERTISLEVEL)
         call obs_vertical_to_level(      ens_handle, precise_dry_mass_get_close, locs(i), my_status(i))
      case (VERTISSCALEHEIGHT)
         call obs_vertical_to_scaleheight(ens_handle, precise_dry_mass_get_close, locs(i), my_status(i))
      case default
         write(string1,*)'unable to convert vertical obs "', which_vert, '"'
         call error_handler(E_ERR,routine,string1,source,revision,revdate)
   end select
enddo

end subroutine convert_vertical_obs


!-----------------------------------------------------------------------
!> This subroutine converts vertical state
!>
!>  in:    ens_handle  - mean ensemble handle
!>  in:    num         - number of locations
!>  inout: locs(:)     - locations
!>  in:    loc_qtys(:) - location quantities
!>  in:    loc_indx(:) - location index
!>  in:    which_vert  - vertical location to convert
!>  out:   istatus     - return status 0 is a successful conversion
!> At present there is no way for this routine to fail. !HK todo FV also
!>  has no fail in this routine. Is this ok?

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

!SENote: A general note. If many states in the same column are being converted this is a remarkably 
!inefficient way to do it.

do i=1,num
   current_vert_type = nint(query_location(locs(i)))

   if ( current_vert_type == which_vert ) cycle
   if ( current_vert_type == VERTISUNDEF) cycle

   select case (which_vert)
      case (VERTISPRESSURE)
         call state_vertical_to_pressure(    ens_handle, ens_size, precise_dry_mass_get_close, locs(i), loc_indx(i))
      case (VERTISHEIGHT)
         call state_vertical_to_height(      ens_handle, ens_size, precise_dry_mass_get_close, locs(i), loc_indx(i))
      case (VERTISLEVEL)
         call state_vertical_to_level(                   ens_size, locs(i), loc_indx(i))
      case (VERTISSCALEHEIGHT)
         call state_vertical_to_scaleheight( ens_handle, ens_size, precise_dry_mass_get_close, locs(i), loc_indx(i))
      case default
         write(string1,*)'unable to convert vertical state "', which_vert, '"'
         call error_handler(E_MSG,routine,string1,source,revision,revdate)
   end select
enddo

istatus = 0

end subroutine convert_vertical_state

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


!-----------------------------------------------------------------------

subroutine fill_gc()

! Subroutine to generate location_types of the cubed sphere grid
! and put them into get_close_type cs_gc, with other derived components.

integer :: c

!SENote: Really don't like all this use of global module storage for communicating among routines
! May want to eliminate some of this.
allocate(cs_locs(ncol), cs_kinds(ncol))

! CS inputs in degrees.
do c=1,ncol
   cs_locs(c)  = set_location(grid_data%lon%vals(c), grid_data%lat%vals(c), MISSING_R8, VERTISUNDEF)
   cs_kinds(c) = 0
enddo

call get_close_init(cs_gc, ncol, coarse_grid, cs_locs)

end subroutine fill_gc

!-----------------------------------------------------------------------
!> given the column and level in the state vector,
!> and the quantity, and the dimensionality of the field (1d, 2d),
!> compute the location of that item.  

function get_location_from_index(column, level, qty, nd)
integer, intent(in) :: column
integer, intent(in) :: level
integer, intent(in) :: qty
integer, intent(in) :: nd
type(location_type) :: get_location_from_index

character(len=*), parameter :: routine = 'get_location_from_index'
real(r8) :: use_vert_val
real(r8) :: my_lon, my_lat, my_vert

! full 2d fields are returned with column/level.
! 1d fields are either surface fields, or if they
! are column integrated values then they are 'undefined'
! in the vertical.

! All fields share the same first coordinate into the column list
my_lon = grid_data%lon%vals(column)
my_lat = grid_data%lat%vals(column)
! For SE 3d spatial fields have a 2d storage
if(nd == 2) then
   my_vert = level
   get_location_from_index = set_location(my_lon, my_lat, my_vert, VERTISLEVEL)
elseif(nd == 1) then
   ! setting the vertical value to missing matches what the previous
   ! version of this code did.  other models choose to set the vertical
   ! value to the model surface elevation at this location:
   !   use_vert_val  = phis(lon_index, lat_index) / gravity not available in SE
   my_vert = MISSING_R8
   ! Add any 2d surface fields to this function
   if(is_surface_field(qty)) then
      get_location_from_index = set_location(my_lon, my_lat, my_vert, VERTISSURFACE)
   else
      get_location_from_index = set_location(my_lon, my_lat, my_vert, VERTISUNDEF)
   endif
else
   write(string1, *) 'state vector field not 1D or 2D and no code to handle other dimensionity'
   write(string2, *) 'dimensionality = ', nd, ' quantity type = ', trim(get_name_for_quantity(qty))
   call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2)
endif
     
end function get_location_from_index


!-----------------------------------------------------------------------
!> this routine converts the column and level index values and a quantity into a state vector
!> offset and gets the ensemble of state values for that offset.  this only
!> gets a single vertical location - if you need to get values which might 
!> have different vertical locations in different ensemble members
!> see get_se_values_from_varid() below.
!> Returns a 0 for OK, returns 12 for my_status for unable to find.

subroutine get_se_values_from_single_level(ens_handle, ens_size, qty, column, level, &
                                        vals, my_status)
type(ensemble_type), intent(in) :: ens_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: qty
integer,             intent(in) :: column
integer,             intent(in) :: level
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

state_indx = get_dart_vector_index(column, level, no_third_dimension, domain_id, varid)

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

subroutine get_se_values_from_varid(ens_handle, ens_size, column, levels, varid, &
                                 vals, my_status)
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: column
integer,             intent(in)  :: levels(ens_size)
integer,             intent(in)  :: varid
real(r8),            intent(out) :: vals(ens_size)
integer,             intent(out) :: my_status(ens_size)

integer(i8) :: state_indx
integer  :: i, j
real(r8) :: temp_vals(ens_size) 
logical  :: member_done(ens_size)

character(len=*), parameter :: routine = 'get_se_values_from_varid:'

! as we get the values for each ensemble member, we set the 'done' flag
! and a good return code. 
my_status(:) = 12
member_done(:) = .false.

! start with levels(1).  get the vals into a temp var.  
! run through 2-N. any other member that has the same level 
! set the outgoing values.  keep a separate flag for which 
! member(s) have been done.  skip to the next undone member 
! and get the state for that level.  repeat until all levels done.

do i=1, ens_size

   if (member_done(i)) cycle
   state_indx = get_dart_vector_index(column, levels(i), no_third_dimension, domain_id, varid)

   !SENote: Do we need error checks like this? Watch out for the ensemble size with levels being too much
   if (state_indx < 0) then
      write(string1,*) 'Should not happen: could not find dart state index from '
      write(string2,*) 'column and level index :', column, levels
      call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2)
      return
   endif

   temp_vals(:) = get_state(state_indx, ens_handle)    ! all the ensemble members for level (i)

   ! start at i, because my ensemble member is clearly at this level.
   ! then continue on to see if any other members are also at this level.
   do j=i, ens_size
      if (member_done(j)) cycle

      if (levels(j) == levels(i)) then
         vals(j) = temp_vals(j)
         member_done(j) = .true.
         my_status(j) = 0
      endif
         
   enddo
enddo

end subroutine get_se_values_from_varid


!-----------------------------------------------------------------------
!> this is just for 2d fields

subroutine get_se_values_from_nonstate_fields(ens_handle, ens_size, column, &
                                           levels, obs_quantity, precise, vals, my_status)
type(ensemble_type),  intent(in)  :: ens_handle
integer,              intent(in)  :: ens_size
integer,              intent(in)  :: column
integer,              intent(in)  :: levels(ens_size)
integer,              intent(in)  :: obs_quantity
logical,              intent(in)  :: precise
real(r8),             intent(out) :: vals(ens_size)
integer,              intent(out) :: my_status(ens_size)

integer  :: imember
real(r8) :: vals_array(ref_nlevels,ens_size)

character(len=*), parameter :: routine = 'get_se_values_from_nonstate_fields:'

vals(:) = MISSING_R8
! This 99 status value can never be returned. Left for backwards consistency
! with FV get_values_from_nonstate_fields but should be verified and removed
my_status(:) = 99

select case (obs_quantity) 
   case (QTY_PRESSURE)
      call cam_se_pressure_levels(ens_handle, ens_size, column, ref_nlevels, &
                               precise, vals_array, my_status)
      if (any(my_status /= 0)) return

      do imember=1,ens_size
         vals(imember) = vals_array(levels(imember), imember)
      enddo

   case (QTY_VERTLEVEL)
      vals(:)      = levels(:)
      my_status(:) = 0

!SENote: Turns out there was no height localization for non-height vertical obs in Manhattan or Classic
!SENote: At present there is no QTY_GEOMETRIC_HEIGHT here as needed to convert to Height

   case default
      write(string1,*)'contact dart support. unexpected error for quantity ', obs_quantity
      call error_handler(E_ERR,routine,string1,source,revision,revdate)

end select

end subroutine get_se_values_from_nonstate_fields

!-----------------------------------------------------------------------
!> internal only version of model interpolate. 
!> does not check for locations too high - return all actual values.

subroutine interpolate_se_values(state_handle, ens_size, location, obs_qty, varid, &
                              precise, interp_vals, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_qty
integer,             intent(in) :: varid
logical,             intent(in) :: precise
real(r8),           intent(out) :: interp_vals(ens_size) 
integer,            intent(out) :: istatus(ens_size)

character(len=*), parameter :: routine = 'interpolate_se_values:'

integer  :: which_vert
integer :: cell_corners(4), i
type(location_type) :: location_copy
real(r8) :: l_weight, m_weight
real(r8) :: lon_lat_vert(3), quad_vals(4, ens_size)

interp_vals(:) = MISSING_R8
istatus(:)     = 99

lon_lat_vert  = get_location(location)
which_vert    = nint(query_location(location)) 

! Do not want to propagate changes to the location back up the calling tree
location_copy = location
call coord_ind_cs(location_copy, obs_qty, cell_corners, l_weight, m_weight)

!Now work on the vertical conversions and getting the vertical index for each ensemble member
call get_se_quad_vals(state_handle, ens_size, varid, obs_qty, cell_corners, &
                   lon_lat_vert, which_vert, precise, quad_vals, istatus)

!SENote: For now return a 12 for all istatus members if there is any failure from get_se_quad_vals
if (any(istatus /= 0)) then
   istatus = 12
   return
endif

! The following uses Jeff's recommended 'generalized quadrilateral interpolation', as in
! http://www.particleincell.com/2012/quad-interpolation/.
! Most of the work is done in create_cs_grid_arrays() and coord_ind_cs().

! Interpolate from the cell's corners to the ob location on the unit square.
! This is done by weighting the field at each corner by the rectangular area
! ((l,m) space) diagonally across the ob location from the corner.
! AKA 'linear area weighting'.

interp_vals(:) = quad_vals(2, :) *       l_weight *          m_weight &
           + quad_vals(1, :) * (1.0_r8 - l_weight)*          m_weight &
           + quad_vals(4, :) * (1.0_r8 - l_weight)*(1.0_r8 - m_weight) &
           + quad_vals(3, :) *           l_weight *(1.0_r8 - m_weight)

end subroutine interpolate_se_values

!-----------------------------------------------------------------------
!>
!> Finds the values at the quad corners for each ensemble member
!>  Returns all ensemble size my_status as 12 if can't find values.
!> Returns 10 for any ensemble member that cannot be interpolated in presure.
!> Returns 11 for any ensemble member that cannot be interpolated in height.
!> Returns 8 for all ensemble members if level is out of range.
!> Returns for all ensemble members 16 for unable to find lower values, 17 for unable to find upper values.

subroutine get_se_quad_vals(state_handle, ens_size, varid, obs_qty, corners, &
                         lon_lat_vert, which_vert, precise, quad_vals, my_status)
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: varid
integer,             intent(in) :: obs_qty
integer,             intent(in) :: corners(4)
real(r8),            intent(in) :: lon_lat_vert(3)
integer,             intent(in) :: which_vert
logical,             intent(in) :: precise
real(r8),           intent(out) :: quad_vals(4, ens_size) !< array of interpolated values
integer,            intent(out) :: my_status(ens_size)

integer  :: icorner, numdims
integer  :: level_one_array(ens_size)
integer  :: four_levs1(4, ens_size), four_levs2(4, ens_size)
real(r8) :: four_vert_fracts(4, ens_size)

character(len=*), parameter :: routine = 'get_se_quad_vals:'

quad_vals(:,:) = MISSING_R8
my_status(:) = 99

! need to consider the case for 2d vs 2d variables
numdims = get_dims_from_qty(obs_qty, varid)

! Now here potentially we have different results for different
! ensemble members.  the things that can vary are dimensioned by ens_size.

if (numdims == 2) then

   ! build 4 columns to find vertical level numbers
   do icorner=1, 4
      call find_se_vertical_levels(state_handle, ens_size, corners(icorner), lon_lat_vert(3), &
                                which_vert, precise, four_levs1(icorner, :), four_levs2(icorner, :), & 
                                four_vert_fracts(icorner, :), my_status)

      if (any(my_status /= 0)) then
         my_status = 12
         return
      endif
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
                                   obs_qty, precise, quad_vals, my_status)

   endif

   !SENote Technically nothing happens after this point anyway? So is this statement needed.
   if (any(my_status /= 0)) return

else if (numdims == 1) then

   if (varid > 0) then
      level_one_array(:) = 1
      do icorner=1, 4
         call get_se_values_from_varid(state_handle,  ens_size, corners(icorner), & 
                                    level_one_array, varid, quad_vals(icorner,:),my_status)

         if (any(my_status /= 0)) return

      enddo

   else ! special 1d case
      !SENote: Is this ever used at present?
      do icorner=1, 4
         call get_se_quad_values(ens_size, corners(icorner), obs_qty, quad_vals(icorner,:))
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
! Returns my_status 0 for success, 16 if unable to find values at lower level
! and 17 if unable to find values at upper level.

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


! SENote This is inefficient since get_se_values_from_varid is messy. Once we know the lower level for an ensmble
! member, we know the level above and could minimize the work.
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

   !SENote: this is more general. The vertical interpolation is linear in level, but this may be biased for 
   !doing interpolation in height, scale_height, or pressure. Is it worth thinking about this?
   call vert_interp(ens_size, vals1, vals2, four_vert_fracts(icorner, :), & 
                    quad_vals(icorner, :))

enddo


end subroutine get_se_four_state_values

!-----------------------------------------------------------------------
!>

! Returns my_status 0 for success, 16 if unable to find values at lower level
! and 17 if unable to find values at upper level.

subroutine get_se_four_nonstate_values(state_handle, ens_size, four_corners, &
                                 four_levs1, four_levs2, four_vert_fracts, &
                                 obs_qty, precise, quad_vals, my_status)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: four_corners(4)
integer,             intent(in) :: four_levs1(4, ens_size), four_levs2(4, ens_size)
real(r8),            intent(in) :: four_vert_fracts(4, ens_size)
integer,             intent(in) :: obs_qty
logical,             intent(in) :: precise
real(r8),           intent(out) :: quad_vals(4, ens_size) !< array of interpolated values
integer,            intent(out) :: my_status(ens_size)

integer  :: icorner
real(r8) :: vals1(ens_size), vals2(ens_size)

character(len=*), parameter :: routine = 'get_se_four_nonstate_values:'

do icorner=1, 4
   call get_se_values_from_nonstate_fields(state_handle, ens_size, four_corners(icorner), &
                              four_levs1(icorner, :), obs_qty, precise, vals1, my_status)
   if (any(my_status /= 0)) then
      my_status(:) = 16   ! cannot retrieve vals1 values
      return
   endif

   call get_se_values_from_nonstate_fields(state_handle,  ens_size, four_corners(icorner), &
                              four_levs2(icorner, :), obs_qty, precise, vals2, my_status)
   if (any(my_status /= 0)) then
      my_status(:) = 17   ! cannot retrieve top values
      return
   endif

   call vert_interp(ens_size, vals1, vals2, four_vert_fracts(icorner, :), &
                    quad_vals(icorner, :))

enddo

end subroutine get_se_four_nonstate_values


!-----------------------------------------------------------------------
!> figure out whether this is a 1d or 2d field based on the quantity.
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
         ! In SE this is a 1 dimensional field
         get_dims_from_qty = 1
      case (QTY_PRESSURE, QTY_GEOMETRIC_HEIGHT)
         ! In SE these are 2d fields
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
!>  This is for 1d special observations quantities not in the state

! For now this can onlu get surface elevation (phis)

subroutine get_se_quad_values(ens_size, column, obs_quantity, vals)
integer,  intent(in) :: ens_size
integer,  intent(in) :: column
integer,  intent(in) :: obs_quantity
real(r8), intent(out) :: vals(ens_size) 

character(len=*), parameter :: routine = 'get_se_quad_values'

integer :: prev_lon, next_lat
real(r8) :: vals1(ens_size), vals2(ens_size)


select case (obs_quantity)
   case (QTY_SURFACE_ELEVATION)
      ! Just return phis for this column
      vals = phis(column)

     !>@todo FIXME:
     ! should this be using gravity at the given latitude? 
     vals = vals / gravity

   case default 
      write(string1, *) 'we can not interpolate qty', obs_quantity
      call error_handler(E_ERR,routine,string1,source,revision,revdate)

end select

end subroutine get_se_quad_values

!-----------------------------------------------------------------------
!> given a column index number, a quantity and a vertical value and type,
!> return which two levels these are between and the fraction across.
!> 
! my_status is 0 for success, 12 for all ensemble members if values cannot be found, 
! 10 for any ensemble member that cannot be interpolated in pressure.
! 11 for any ensemble member that cannot be interpolated in height.
! 8 for all ensemble members if cannot be interpolated in level.

subroutine find_se_vertical_levels(ens_handle, ens_size, column, vert_val, &
                                which_vert, precise, levs1, levs2, vert_fracts, my_status)
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: column
real(r8),            intent(in)  :: vert_val
integer,             intent(in)  :: which_vert
logical,             intent(in)  :: precise
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

! assume the worst
levs1(:)    = MISSING_I
levs2(:)    = MISSING_I
vert_fracts(:) = MISSING_R8
my_status(:)   = 98


! ref_nlevels is the number of vertical levels (midlayer points)

level_one = 1

select case (which_vert)

   case(VERTISPRESSURE)
      ! construct a pressure column here and find the model levels that enclose this value
      call get_se_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, column, level_one, &
         surf_pressure, status1)

      ! Returns all my_status members as 12 if unable to find the value
      if (status1 /= 0) then
         my_status(:) = status1
         return
      endif

      if(dry_mass_vertical_coordinate .and. precise) then
         call build_dry_mass_pressure_columns(ens_handle, ens_size, ref_nlevels, column, surf_pressure, &
            pressure_array, status1)
         ! All output variables set to missing already, just return with all my_status as 12
         if(status1 /= 0) then
            my_status = 12
            return
         endif
      else
         call build_cam_pressure_columns(ens_size, surf_pressure, ref_nlevels, pressure_array)
      endif

      do imember = 1, ens_size
         ! Returns my_status 10 if unable to interpolate in this column
         call pressure_to_level(ref_nlevels, pressure_array(:, imember), vert_val, & 
                                levs1(imember), levs2(imember), vert_fracts(imember), my_status(imember))
      enddo

      !SENote: Can we somehow get all of these disruptive debug statements out of the code? Preprocess?
      if (debug_level > 100) then
         do k = 1,ens_size
            print*, 'ISPRESSURE levs1(k), levs2(k), vert_fracts(k), vert_val', &
                     levs1(k), levs2(k), vert_fracts(k), vert_val, pressure_array(levs1(k) , k), pressure_array(levs2(k), k) 
          enddo
      endif

   case(VERTISHEIGHT)
      ! construct a height column here and find the model levels that enclose this value
      ! All my_status are returned as 12 if a failure
      call cam_se_height_levels(ens_handle, ens_size, column, ref_nlevels, &
                             precise, height_array, my_status)

      !>@todo FIXME let successful members continue?
      if (any(my_status /= 0)) return

      if (debug_level > 400) then
         do k = 1,ref_nlevels
            print*, 'ISHEIGHT: ', k, height_array(k,1)
         enddo
      endif

      do imember=1, ens_size
         ! Returns 11 if unable to interpolate in this column
         call height_to_level(ref_nlevels, height_array(:, imember), vert_val, & 
                             levs1(imember), levs2(imember), vert_fracts(imember), &
                             my_status(imember))
      enddo

      if (debug_level > 100) then
         do k = 1,ens_size
            print*, 'ISHEIGHT ens#, levs1(#), levs2(#), vert_fracts(#), top/bot height(#)', &
                     k, levs1(k), levs2(k), vert_fracts(k), height_array(levs2(k),k), height_array(levs1(k), k)
         enddo
      endif

      !>@todo FIXME let successful members continue?
      
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

!SENote: This subroutine is only called from one place and only for 2d fields so this next block can't be reached
! This allows removal of the obs_qty and varid arguments from the call
   ! 1d fields for SE
   !case(VERTISUNDEF, VERTISSURFACE)
      !if (get_dims_from_qty(obs_qty, var_id) == 2) then
         !levs1(:) = ref_nlevels - 1
         !levs2(:) = ref_nlevels
         !vert_fracts(:) = 1.0_r8
         !my_status(:) = 0
      !else
         !my_status(:) = 4 ! can not get vertical levels
      !endif

   case default
      write(string1, *) 'unsupported vertical type: ', which_vert
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
      
end select

! by this time someone has already set my_status(), good or bad.

end subroutine find_se_vertical_levels


!-----------------------------------------------------------------------
!> Compute pressure column for the dry mass vertical coordinate option
!>
!> this version does all ensemble members at once.


subroutine build_dry_mass_pressure_columns(ens_handle, ens_size, nlevels, column, surf_pressure, pressure, status)
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: nlevels
integer,             intent(in)  :: column
real(r8),            intent(in)  :: surf_pressure(ens_size)
real(r8),            intent(out) :: pressure(nlevels, ens_size)
integer,             intent(out) :: status

real(r8) :: specific_humidity(ens_size), cldliq(ens_size), cldice(ens_size), sum_specific_water_ratios(ens_size)
real(r8) :: a_width(nlevels), b_width(nlevels)
real(r8) :: sum_dry_mix_ratio(nlevels, ens_size), dry_mass_top, mass_diff_term, denom, numer
real(r8) :: dry_mass_sfc(ens_size), half_pressure(nlevels + 1)
integer  :: k, n, istatus

! Building pressure columns for dry mass vertical coordinate; Add in references to Lauritzen and pointer
! to the document on the algorithm

! Begin by getting column values for specific mixing ratio for water vapor (specific humidity), cloud liquid 
! and cloud ice. For more accuracy could also include rain, snow, and any other tracers.
! Note that other tracers in the dry mass cam/se have a dry mixing ratio, not specific mixing ratio (although
! Peter plans to change this for consistency).

!SENote: The following should all be 0, but A's at levels 2, 3, 4, and 5 are not.
! This is confirmed in the caminput.nc file and is a PROBLEM.
! Some tests on the A and B coefficients
!do k = 1, nlevels
   !write(*, *) k, (grid_data%hyai%vals(k) + grid_data%hyai%vals(k+1))/2.0_r8 - grid_data%hyam%vals(k)
   !write(*, *) k, (grid_data%hybi%vals(k) + grid_data%hybi%vals(k+1))/2.0_r8 - grid_data%hybm%vals(k)
!enddo
!stop

! For now, will fail all ensemble members and levels if any level/member fails

! Need the water tracer specific mixing ratios for ever level in the column to compute their mass sum
do k = 1, nlevels

   ! Specific Humidity
   call get_se_values_from_single_level(ens_handle, ens_size, QTY_SPECIFIC_HUMIDITY, column, k, &
      specific_humidity, status)
   if (status /= 0) then
      pressure = MISSING_R8
      return
   endif
   
   ! Cloud liquid
   call get_se_values_from_single_level(ens_handle, ens_size, QTY_CLOUD_LIQUID_WATER, column, k, &
      cldliq, status)
   if (status /= 0) then
      pressure = MISSING_R8
      return
   endif

   ! Cloud ice
   call get_se_values_from_single_level(ens_handle, ens_size, QTY_CLOUD_ICE, column, k, &
      cldice, status)
   if (status /= 0) then
      pressure = MISSING_R8
      return
   endif

   ! Compute the sum of the dry mixing ratio of dry air plus all the water tracers (ref. to notes)
   sum_specific_water_ratios = specific_humidity(:) + cldliq(:) + cldice(:) 
   sum_dry_mix_ratio(k, :) = 1.0_r8 + sum_specific_water_ratios(:) / (1.0_r8 - sum_specific_water_ratios(:))

   ! Compute the A 'width' and B 'width' of each level
   a_width(k) = grid_data%hyai%vals(k + 1) - grid_data%hyai%vals(k)
   b_width(k) = grid_data%hybi%vals(k + 1) - grid_data%hybi%vals(k)

enddo

! Compute the dry mass at the bottom of the column for each enseble member
! Do we need to worry about latitudinal variation in g?
! Nothing but dry air above the model top
dry_mass_top = ref_model_top_pressure / gravity
do n = 1, ens_size
   mass_diff_term = (surf_pressure(n) - ref_model_top_pressure) / gravity
   numer = mass_diff_term - dry_mass_top * sum(a_width(:) * sum_dry_mix_ratio(:, n)) / grid_data%hyai%vals(1)
   denom = sum(b_width(:) * sum_dry_mix_ratio(:, n))
   dry_mass_sfc(n) = numer / denom

   ! Now compute the pressure columns
   half_pressure(1) = ref_model_top_pressure
   do k = 1, nlevels
      half_pressure(k + 1) = half_pressure(k) + &
         !gravity * (a_width(k)*dry_mass_top + b_width(k)*dry_mass_sfc(n)) * sum_dry_mix_ratio(k, n)
         ! SENote: NEXT LINE IS BELIEVED TO BE CORRECT BUT NEEDS TO BE VETTED WITH CGD
         gravity * (a_width(k)*dry_mass_top / grid_data%hyai%vals(1) + b_width(k)*dry_mass_sfc(n)) * sum_dry_mix_ratio(k, n)
      pressure(k, n) = (half_pressure(k) + half_pressure(k + 1)) / 2
   end do
end do

status = 0

end subroutine build_dry_mass_pressure_columns


!-----------------------------------------------------------------------
!> Compute the heights at pressure midpoints
!>
!> this version does all ensemble members at once.
!> Returns my_status 12 for all members if unable to compute levels.

subroutine cam_se_height_levels(ens_handle, ens_size, column, nlevels, precise, height_array, my_status) 
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: column
integer,             intent(in)  :: nlevels
logical,             intent(in)  :: precise
real(r8),            intent(out) :: height_array(nlevels, ens_size)
integer,             intent(out) :: my_status(ens_size)

integer  :: k, level_one, imember, status1
real(r8) :: surface_elevation(1)
real(r8) :: surface_pressure(ens_size), mbar(nlevels, ens_size)
real(r8) :: pressure(nlevels, ens_size)
real(r8) :: tv(nlevels, ens_size)  ! Virtual temperature, top to bottom

! this is for surface obs
level_one = 1

! Get the surface pressure at this column; Returns status1 12 if value cannot be found
call get_se_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, column, level_one, &
   surface_pressure, status1)
if(status1 /= 0) then 
   my_status = status1
   return
endif

! get the surface elevation from the phis
call get_se_quad_values(1, column, QTY_SURFACE_ELEVATION, surface_elevation)

! Returns status1 12 if unsuccessful in getting needed quantities
call compute_se_virtual_temperature(ens_handle, ens_size, column, nlevels, tv, status1)

if (status1 /= 0) then
   my_status = status1
   return
endif

! Build the pressure columns for the entire ensemble
if(dry_mass_vertical_coordinate .and. precise) then
   call build_dry_mass_pressure_columns(ens_handle, ens_size, nlevels, column, surface_pressure, &
      pressure, status1)
   if(status1 /= 0) then
      my_status = 12
      height_array = MISSING_R8
      return
   endif
else
   call build_cam_pressure_columns(ens_size, surface_pressure, nlevels, pressure)
endif

if (use_variable_mean_mass) then
   call compute_se_mean_mass(ens_handle, ens_size, column, nlevels, mbar, status1)
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
call gph2gmh(height_array, grid_data%lat%vals(column))

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
! The remaining (private) interfaces come last.
! None of the private interfaces need to call static_init_model()
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!>
!>

subroutine read_cam_phis_array(phis_filename)
character(len=*),   intent(in)    :: phis_filename

character(len=*), parameter :: routine = 'read_cam_phis_array'

integer :: ncid, nsize(3)   ! lon, lat, time !HK todo nope? space filling curve nsize=1

ncid = nc_open_file_readonly(phis_filename, routine)

call nc_get_variable_size(ncid, 'PHIS', nsize(:), routine)

allocate( phis(nsize(1)))

call nc_get_variable(ncid, 'PHIS', phis, routine)

call nc_close_file(ncid, routine)

end subroutine read_cam_phis_array


!-----------------------------------------------------------------------
!> Compute the virtual temperature at the midpoints
!>
!> this version does all ensemble members at once.
!>

subroutine compute_se_virtual_temperature(ens_handle, ens_size, column, nlevels, tv, istatus)

type(ensemble_type), intent(in)   :: ens_handle
integer,             intent(in)   :: ens_size
integer,             intent(in)   :: column
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
   call get_se_values_from_single_level(ens_handle, ens_size, QTY_TEMPERATURE, column, k, &
      temperature, istatus)

   if (istatus < 0) return

   ! specific humidity
   call get_se_values_from_single_level(ens_handle, ens_size, QTY_SPECIFIC_HUMIDITY, column, k, &
      specific_humidity, istatus)
   
   if (istatus < 0) return

   !>tv == virtual temperature.
   tv(k,:) = temperature(:)*(1.0_r8 + rr_factor*specific_humidity(:))
enddo

end subroutine compute_se_virtual_temperature


!-----------------------------------------------------------------------
!> loop through all levels to get the mean mass.
!>


subroutine compute_se_mean_mass(ens_handle, ens_size, column, nlevels, mbar, istatus)
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: column
integer,             intent(in)  :: nlevels
real(r8),            intent(out) :: mbar(nlevels, ens_size)
integer,             intent(out) :: istatus

integer :: k, this_qty
real(r8) :: mmr_o1(ens_size, nlevels), &
            mmr_o2(ens_size, nlevels), &
            mmr_h1(ens_size, nlevels), &
            mmr_n2(ens_size, nlevels)
real(r8) :: O_molar_mass, O2_molar_mass, H_molar_mass, N2_molar_mass 
integer  :: my_status, varid

character(len=*), parameter :: routine = 'compute_se_mean_mass'

!SENote: This subroutine has not been tested yet for SE. Need to work with WACCM folks to test.
call error_handler(E_ERR, routine, 'Subroutine has not been tested', source, revision, revdate)

! Default is successful return
istatus = 0

! do this outside the subroutine?  it never changes throughout the
! run of the program
!SENote: could do an initialization with save storage
! See if the quantities can be interpolated. 
call ok_to_interpolate(QTY_ATOMIC_OXYGEN_MIXING_RATIO, varid, my_status)
if(my_status /= 0) call error_handler(E_ERR, routine, 'Cannot get QTY_ATOMIC_OXYGEN_MIXING_RATIO', source, revision, revdate)
O_molar_mass  = get_molar_mass(QTY_ATOMIC_OXYGEN_MIXING_RATIO)

call ok_to_interpolate(QTY_MOLEC_OXYGEN_MIXING_RATIO, varid, my_status)
if(my_status /= 0) call error_handler(E_ERR, routine, 'Cannot get QTY_MOLEC_OXYGEN_MIXING_RATIO', source, revision, revdate)
O2_molar_mass = get_molar_mass(QTY_MOLEC_OXYGEN_MIXING_RATIO)

call ok_to_interpolate(QTY_ATOMIC_H_MIXING_RATIO, varid,my_status)
if(my_status /= 0) call error_handler(E_ERR, routine, 'Cannot get QTY_ATOMIC_H_MIXING_RATIO', source, revision, revdate)
H_molar_mass  = get_molar_mass(QTY_ATOMIC_H_MIXING_RATIO)

call ok_to_interpolate(QTY_NITROGEN, varid, my_status)
if(my_status /= 0) call error_handler(E_ERR, routine, 'Cannot get QTY_NITROGEN', source, revision, revdate)
N2_molar_mass = get_molar_mass(QTY_NITROGEN)
   


! High topped models (WACCM-X) need to account for the changing composition 
! of the atmosphere with height.  This requires several variables from the
! initial file, which may not be available from low topped models.
do k = 1, nlevels

   call get_se_values_from_single_level(ens_handle, ens_size, QTY_ATOMIC_OXYGEN_MIXING_RATIO, &
      column, k, mmr_o1(:, k), istatus)
   if (istatus /= 0) return
   !print *, 'mmr: ', trim(get_name_for_quantity(this_qty)), mmr_o1(1, k)
   
   call get_se_values_from_single_level(ens_handle, ens_size, QTY_MOLEC_OXYGEN_MIXING_RATIO, &
      column, k, mmr_o2(:, k), istatus)
   if (istatus /= 0) return
   !print *, 'mmr: ', trim(get_name_for_quantity(this_qty)), mmr_o2(1, k)
   
   call get_se_values_from_single_level(ens_handle, ens_size, QTY_ATOMIC_H_MIXING_RATIO, &
      column, k, mmr_h1(:, k), istatus)
   if (istatus /= 0) return
   !print *, 'mmr: ', trim(get_name_for_quantity(this_qty)), mmr_h1(1, k)
   
   mmr_n2(:,k) = 1.0_r8 - (mmr_o1(:,k) + mmr_o2(:,k) + mmr_h1(:,k))
   mbar(k,:) = 1.0_r8/( mmr_o1(:,k)/O_molar_mass  &
                      + mmr_o2(:,k)/O2_molar_mass &
                      + mmr_h1(:,k)/H_molar_mass  &
                      + mmr_n2(:,k)/N2_molar_mass)
enddo

end subroutine compute_se_mean_mass

!--------------------------------------------------------------------

subroutine state_vertical_to_pressure(ens_handle, ens_size, precise, location, location_indx)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
logical,             intent(in)    :: precise
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx

integer  :: column, level, myqty, level_one, status1
integer  :: my_status(ens_size)
real(r8) :: pressure_array(ref_nlevels), surface_pressure(ens_size)


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
                            precise, pressure_array, my_status)

   call set_vertical(location, pressure_array(level), VERTISPRESSURE)
endif

end subroutine state_vertical_to_pressure

!--------------------------------------------------------------------

subroutine state_vertical_to_height(ens_handle, ens_size, precise, location, location_indx)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
logical,             intent(in)    :: precise
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx

integer  :: column, level, my_status(ens_size)
real(r8) :: height_array(ref_nlevels, ens_size)

! build a height column and a pressure column and find the levels
call get_model_variable_indices(location_indx, column, level, no_third_dimension)

call cam_se_height_levels(ens_handle, ens_size, column, ref_nlevels, &
                       precise, height_array, my_status) 

!>@todo FIXME this can only be used if ensemble size is 1
call set_vertical(location, height_array(level, 1), VERTISHEIGHT)

end subroutine state_vertical_to_height

!--------------------------------------------------------------------

subroutine state_vertical_to_scaleheight(ens_handle, ens_size, precise, location, location_indx)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
logical,             intent(in)    :: precise
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
                               precise, pressure_array, my_status)
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
                               precise, pressure_array, my_status)
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

integer  :: column, level

!>@todo FIXME qty is currently unused.  if we need it, its here.
!>if we really don't need it, we can remove it.  all the other
!>corresponding routines like this use it. 

call get_model_variable_indices(location_indx, column, level, no_third_dimension)

call set_vertical(location, real(level, r8), VERTISLEVEL)

end subroutine state_vertical_to_level


!-----------------------------------------------------------------------
!> Compute the pressure values at midpoint levels
!>
!> this version does all ensemble members at once.

subroutine cam_se_pressure_levels(ens_handle, ens_size, column, nlevels, &
                               precise, pressure_array, my_status) 
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: column
integer,             intent(in)  :: nlevels
logical,             intent(in)  :: precise
real(r8),            intent(out) :: pressure_array(nlevels, ens_size)
integer,             intent(out) :: my_status(ens_size)

integer     :: level_one, status1
real(r8)    :: surface_pressure(ens_size)

! this is for surface obs
level_one = 1

! get the surface pressure from the ens_handle
call get_se_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, column, level_one, &
                                     surface_pressure, status1)

if (status1 /= 0) then
   my_status(:) = status1
   return
endif

if(dry_mass_vertical_coordinate .and. precise) then
   call build_dry_mass_pressure_columns(ens_handle, ens_size, ref_nlevels, column, surface_pressure, &
      pressure_array, status1)
   if(status1 /= 0) then
      my_status = 12
      pressure_array = MISSING_R8
      return
   endif
else
   call build_cam_pressure_columns(ens_size, surface_pressure, ref_nlevels, pressure_array)
endif

! No error returns available if we get here, so all good
my_status(:) = 0

end subroutine cam_se_pressure_levels

!--------------------------------------------------------------------

subroutine obs_vertical_to_pressure(ens_handle, precise, location, my_status)

type(ensemble_type), intent(in)    :: ens_handle
logical,             intent(in)    :: precise
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

call interpolate_se_values(ens_handle, ens_size, location, &
                        qty, varid, precise, pressure_array(:), status(:))

if (status(1) /= 0) then
   my_status = status(1)
   return
endif

call set_vertical(location, pressure_array(1), VERTISPRESSURE)

my_status = 0

end subroutine obs_vertical_to_pressure

!--------------------------------------------------------------------

subroutine obs_vertical_to_height(ens_handle, precise, location, my_status)
type(ensemble_type), intent(in)    :: ens_handle
logical,             intent(in)    :: precise
type(location_type), intent(inout) :: location
integer,             intent(out)   :: my_status

integer  :: varid, ens_size, status(1)
real(r8) :: height_array(1)

character(len=*), parameter :: routine = 'obs_vertical_to_height'

!SENote Does this work for the FV? 
! Doesn't actually appear to work right for the FV which just blasts through a failed search 
! for the height field. This could be fixed.
ens_size = 1

call ok_to_interpolate(QTY_GEOMETRIC_HEIGHT, varid, my_status)
if (my_status /= 0) return

call interpolate_se_values(ens_handle, ens_size, location, &
                        QTY_GEOMETRIC_HEIGHT, varid, precise, height_array(:), status(:))
if (status(1) /= 0) then
   my_status = status(1)
   return
endif

call set_vertical(location, height_array(1), VERTISHEIGHT)

my_status = 0

end subroutine obs_vertical_to_height

!--------------------------------------------------------------------

subroutine obs_vertical_to_level(ens_handle, precise, location, my_status)
type(ensemble_type), intent(in)    :: ens_handle
logical,             intent(in)    :: precise
type(location_type), intent(inout) :: location
integer,             intent(out)   :: my_status

integer  :: varid, ens_size, status(1)
real(r8) :: level_array(1)

ens_size = 1
varid = -1
!SENote; This does not work yet so can't do vertical localization in level
! of observations that are not on level.
!SENote: could be implemented but is there any demand?
call error_handler(E_ERR, 'obs_vertical_to_level',  &
           'Localization in level for obs not on levels is not implemented', &
           source, revision, revdate)

!SENote: This has not been checked, just swapped the call to interpolate_values
call interpolate_se_values(ens_handle, ens_size, location, &
                        QTY_VERTLEVEL, varid, precise, level_array(:), status(:))
if (status(1) /= 0) then
   my_status = status(1)
   return
endif

call set_vertical(location, level_array(1), VERTISLEVEL)

my_status = 0

end subroutine obs_vertical_to_level

!--------------------------------------------------------------------

subroutine obs_vertical_to_scaleheight(ens_handle, precise, location, my_status)
type(ensemble_type), intent(in)    :: ens_handle
logical,             intent(in)    :: precise
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
      call interpolate_se_values(ens_handle, ens_size, location, ptype, varid1, &
                              precise, pressure_array(:), status(:))
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
         call interpolate_se_values(ens_handle, ens_size, location, QTY_PRESSURE, varid1, &
                                    precise, pressure_array(:), status(:))
         if (status(1) /= 0) then
            my_status = status(1)
            return
         endif
      endif
                                 
      call ok_to_interpolate(QTY_SURFACE_PRESSURE, varid2, my_status)
      if (my_status /= 0) return
      
      call interpolate_se_values(ens_handle, ens_size, location, QTY_SURFACE_PRESSURE, varid2, &
                                    precise, surface_pressure_array(:), status(:))
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

type(location_type) :: base_loc(1)
integer :: base_qty(1), base_type(1), status(1)

! SENote: Only reason this is needed is to do the conversion from a scalar to a 1-element array. Annoying.

! these need to be arrays.  kinda a pain.
base_loc(1) = loc
base_type(1) = otype
base_qty(1) = get_quantity_for_type_of_obs(otype)

call convert_vertical_obs(ens_handle, 1, base_loc, base_qty, base_type, vert_type, status)

status1 = status(1)

if (status1 /= 0) return

loc = base_loc(1)

end subroutine convert_vert_one_obs


!-----------------------------------------------------------------------------------------------------
! Routines for computing horizontal grid box location with cubed sphere spectral element grids follows
!-----------------------------------------------------------------------------------------------------


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
integer :: nc_file_ID
real(r8) :: dist, angle
real(r8) :: bearings(3), x_planar(3), y_planar(3)

! ncol = number of nodes/corners/grid points.  Global storage.
! corners = the names of the corners associated with each cell center
! neighbors = the nodes around each node which partner to make the sides of the cells
!             which may contain an observation.

! Get array of corner nodes/columns which define the cells (identified by 'center').
if (file_exist(homme_map_file)) then
   nc_file_ID = nc_open_file_readonly(homme_map_file)
   ncorners = nc_get_dimension_size(nc_file_ID, 'ncorners')
   ncenters = nc_get_dimension_size(nc_file_ID, 'ncenters')

   if (ncenters /= (ncol -2) ) then
      write(string1, *) trim(homme_map_file),' ncenters inconsistent with ncol-2 ', ncenters, ncol
      call error_handler(E_ERR,'create_cs_grid_arrays',string1,source,revision,revdate)
   endif

   ! Allocate array for the homme mapping file contents
   allocate(corners(ncenters, ncorners))

   ! Read it in
   call nc_get_variable(nc_file_ID, 'element_corners', corners)
   call nc_close_file(nc_file_ID)

   allocate(num_nghbrs           (ncol), centers(max_neighbors,ncol),  a(3,ncorners,ncenters),   &
               b(2,ncorners,ncenters),  x_ax_bearings(ncorners,ncenters))

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

subroutine nc_read_cs_grid_file()

! Read the number of neighbors, corners, centers, a and b coefficients, and x_ax_bearings
! from a netCDF file once for this grid at the beginning of the assimilation.

integer :: nc_file_ID,  max_nghbrs, local_ncol

! Open the file for reading
nc_file_ID = nc_open_file_readonly(trim(cs_grid_file), 'reading the cs_grid_file')

! Get the number of centers and number of corners and check the number of columns
ncenters = nc_get_dimension_size(nc_file_ID, 'ncenters')
ncorners = nc_get_dimension_size(nc_file_ID, 'ncorners')
local_ncol = nc_get_dimension_size(nc_file_ID, 'ncol')
max_nghbrs = nc_get_dimension_size(nc_file_ID, 'max_neighbors')

! Check value against the namelist/parameter value.
if (max_nghbrs /= max_neighbors) then
   write(string1, *) trim(cs_grid_file),' max_nghbrs does not match max_neighbors', &
         max_nghbrs,max_neighbors
   call error_handler(E_ERR,'nc_read_cs_grid_file',string1,source,revision,revdate)
endif

! Check value against the namelist/parameter value.
if (local_ncol /= ncol) then
   write(string1, *) trim(cs_grid_file),' ncol in cs_grid_file does not match the one in caminput.nc', &
         local_ncol, ncol
   call error_handler(E_ERR,'nc_read_cs_grid_file',string1,source,revision,revdate)
endif

! Allocate space for all the cs geometry variables
allocate (corners(ncenters, ncorners),  num_nghbrs(ncol), centers(max_neighbors, ncol), &
             x_ax_bearings(ncorners, ncenters), a(3, ncorners, ncenters), b(2, ncorners, ncenters))

! Read in the values for these fields
call nc_get_variable(nc_file_ID, 'corners', corners)
call nc_get_variable(nc_file_ID, 'num_nghbrs', num_nghbrs)
call nc_get_variable(nc_file_ID, 'centers', centers)
call nc_get_variable(nc_file_ID, 'x_ax_bearings', x_ax_bearings)
call nc_get_variable(nc_file_ID, 'a', a)
call nc_get_variable(nc_file_ID, 'b', b)
call nc_close_file(nc_file_ID, 'closing cs_grid_file')

end subroutine nc_read_cs_grid_file

!-----------------------------------------------------------------------

subroutine nc_write_cs_grid_file(cs_grid_file, homme_map_file)

! Write out the number of neighbors, the neighbors, corners, centers, and bearings
! to a netCDF file once for this grid at the beginning of the assimilation.

character(len=*), intent(in) :: cs_grid_file
character(len=*), intent(in) :: homme_map_file

integer :: nc_file_ID

! Create the file
nc_file_ID = nc_create_file(trim(cs_grid_file), 'creating cs_grid_file')

! Define the dimensions
call nc_define_dimension(nc_file_ID, 'ncenters', ncenters)
call nc_define_dimension(nc_file_ID, 'ncorners', ncorners)
call nc_define_dimension(nc_file_ID, 'max_neighbors', max_neighbors)
call nc_define_dimension(nc_file_ID, 'ncol', ncol)
call nc_define_dimension(nc_file_ID, 'ncoef_a', 3)
call nc_define_dimension(nc_file_ID, 'ncoef_b', 2)

! Write Global Attributes
call nc_add_global_attribute(nc_file_ID, 'title', trim(cs_grid_file))
call nc_add_global_attribute(nc_file_ID, 'model_mod_source', source)
call nc_add_global_attribute(nc_file_ID, 'model_mod_revision', revision)
call nc_add_global_attribute(nc_file_ID, 'model_mod_revdate', revdate)
call nc_add_global_attribute(nc_file_ID, 'elements_per_cube_edge', ne)
call nc_add_global_attribute(nc_file_ID, 'nodes_per_element_edge', np)
call nc_add_global_attribute(nc_file_ID, 'HommeMapping_file', homme_map_file)

! Create variables and attributes.
call nc_define_integer_variable(nc_file_ID, 'num_nghbrs', 'ncol')
call nc_add_attribute_to_variable(nc_file_ID, 'num_nghbrs', 'long_name', 'number of neighbors of each node/column')
call nc_add_attribute_to_variable(nc_file_ID, 'num_nghbrs', 'units', 'nondimensional')
call nc_add_attribute_to_variable(nc_file_ID, 'num_nghbrs', 'valid_range', (/1, max_neighbors/))

call nc_define_integer_variable(nc_file_ID, 'centers', (/'max_neighbors', 'ncol         '/))
call nc_add_attribute_to_variable(nc_file_ID, 'centers', 'long_name', 'cells which use node/column as a corner')
call nc_add_attribute_to_variable(nc_file_ID, 'centers', 'units', 'nondimensional')
call nc_add_attribute_to_variable(nc_file_ID, 'centers', 'valid_range', (/1, ncenters/))
call nc_add_attribute_to_variable(nc_file_ID, 'centers', 'missing_value', MISSING_I)

call nc_define_integer_variable(nc_file_ID, 'corners', (/'ncenters', 'ncorners'/))
call nc_add_attribute_to_variable(nc_file_ID, 'corners', 'long_name', 'corners/nodes of each cell')
call nc_add_attribute_to_variable(nc_file_ID, 'corners', 'units', 'nondimensional')
call nc_add_attribute_to_variable(nc_file_ID, 'corners', 'valid_range', (/1, ncol/))
call nc_add_attribute_to_variable(nc_file_ID, 'corners', 'missing_value', MISSING_I)

call nc_define_double_variable(nc_file_ID, 'a', (/'ncoef_a ', 'ncorners', 'ncenters'/))
call nc_add_attribute_to_variable(nc_file_ID, 'a', 'long_name', &
   'Coefficients of mapping from planar x coord to unit square')
call nc_add_attribute_to_variable(nc_file_ID, 'a', 'units', 'nondimensional')
call nc_add_attribute_to_variable(nc_file_ID, 'a', 'missing_value', MISSING_R8)

call nc_define_double_variable(nc_file_ID, 'b', (/'ncoef_b ', 'ncorners', 'ncenters'/))
call nc_add_attribute_to_variable(nc_file_ID, 'b', 'long_name', &
   'Coefficients of mapping from planar y coord to unit square')
call nc_add_attribute_to_variable(nc_file_ID, 'b', 'units', 'nondimensional')
call nc_add_attribute_to_variable(nc_file_ID, 'b', 'missing_value', MISSING_R8)

call nc_define_double_variable(nc_file_ID, 'x_ax_bearings', (/'ncorners', 'ncenters'/))
call nc_add_attribute_to_variable(nc_file_ID, 'x_ax_bearings', 'long_name', &
   'bearing (clockwise from North) from origin node(corner 4) of each mapping to corner 3')
call nc_add_attribute_to_variable(nc_file_ID, 'x_ax_bearings', 'units', 'radians')
call nc_add_attribute_to_variable(nc_file_ID, 'x_ax_bearings', 'valid_range', (/-PI, PI/))
call nc_add_attribute_to_variable(nc_file_ID, 'x_ax_bearings', 'missing_value', MISSING_R8)

! Fill 'em up
call nc_end_define_mode(nc_file_ID)
call nc_put_variable(nc_file_ID, 'num_nghbrs', num_nghbrs)
call nc_put_variable(nc_file_ID, 'centers', centers)
call nc_put_variable(nc_file_ID, 'corners', corners)
call nc_put_variable(nc_file_ID, 'a', a)
call nc_put_variable(nc_file_ID, 'b', b)
call nc_put_variable(nc_file_ID, 'x_ax_bearings', x_ax_bearings)
call nc_close_file(nc_file_ID)

end subroutine nc_write_cs_grid_file

!-----------------------------------------------------------------------

subroutine coord_ind_cs(obs_loc, obs_kind, cell_corners, l_weight, m_weight)

! Find the corners of the cell which contains the location.

! Variables needed by loc_get_close_obs:
type(location_type),  intent(inout)  :: obs_loc
integer,              intent(in)  :: obs_kind
integer,              intent(out) :: cell_corners(4)
real(r8),             intent(out) :: l_weight
real(r8),             intent(out) :: m_weight

! Output from loc_get_close_obs
integer  :: num_close

! SENote: Need to reduce the memory usage for these for standard configurations. Classic notes follow:
! It would be nice if these could be smaller, but I don't know what number would work.
! It has to be large enough to accommodate all of the grid points that might lie
! within 2xcutoff; resolution and location dependent.
! The size must be specified here; (:) yields an error, and 'allocatable' doesn't help.
integer, allocatable  :: close_ind(:)
real(r8), allocatable :: dist(:)

! dist_# in radians (Can't be initialized here or they will get the 'save' property,
! and will not be reset during subsequent entries to this subroutine.)
real(r8) :: dist_1, dist_2
real(r8) :: lon_lat_lev(3)
integer  :: k, k1, k2, closest, closest2, origin
logical  :: found_cell

lon_lat_lev = get_location(obs_loc)

! See whether this obs_ is a state variable.
! This could be done by 2 calls to minloc(dist), with the 2nd call using a mask
! to prevent finding the closest, which was found in the first call.
! But would those 2 intrinsic searches through dist be faster than my 1 explicit search?

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
! Namelist or documented. Note that there is an error trap for failure that recommends changing to 
! approximate_distance false below.
call get_close(cs_gc, obs_loc, 1, cs_locs, cs_kinds, &
                       num_close, close_ind, dist)

dist_1 = 10.0_r8
dist_2 = 10.0_r8
closest = MISSING_I
k1 = MISSING_I

! Keep track of k1, k2, and distances in this search.
! Assign closest and closest2 afterwards.
if (num_close <= 0) then
   write(string1,*) "Can't find enclosing quadrilatersl. Unusable num_close, obs_kind : ",num_close, obs_kind
   call write_location(0, obs_loc, charstring=string2)
   write(string3,*) 'Setting namelist approximate_distance = .false. might help: dist(1) = ',dist(1)
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
   string3 = 'Setting namelist approximate_distance = .false. might help: closest node location = '//string3
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
                             lon_lat_lev(1), lon_lat_lev(2), found_cell, origin, l_weight, m_weight)
   if (found_cell) exit Cloop
enddo Cloop

! Try the 2nd closest point, if the first failed.
if ((.not.found_cell) .and. closest2 /= MISSING_I) then

   Second_closest: do k=1,num_nghbrs(closest2)
      call unit_square_location(centers(k,closest2), closest2, obs_loc,         &
                                lon_lat_lev(1), lon_lat_lev(2), found_cell, origin, l_weight,m_weight)
      if (found_cell) then
         ! Put '2nd closest' information into 'closest'.
         dist_1 = dist_2
         closest = closest2

         write(string1,'(A,2F10.7,2I8,1p2E12.4)') &
              'Using 2nd closest node to the ob: l, m, closest2, origin2 = ', &
              l_weight, m_weight, closest, origin
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
   string2 = 'Setting namelist approximate_distance = .false. might help.'
   call error_handler(E_ERR, 'coord_ind_cs', string1,source,revision,revdate, text2 = string2)
endif

deallocate(close_ind, dist)

end subroutine coord_ind_cs

!-----------------------------------------------------------------------------------------------------

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
   if (l_neg == 0.0_r8 .and. my_task_id() == 0) then
      write(string1,'(A,I6,1X,1p4E12.4)') 'l_neg cell, x_o - a(2)*m = ',cell, x_o ,a(2,origin,cell),m
      call error_handler(E_MSG, 'unit_square_location', string1,source,revision,revdate)
   endif

endif

! Informational output, if the observation is exactly on the m-axis
!SENote: Why does this message get printed a billion times in CLASSIC?
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

subroutine solve_quadratic(a, b, c, r1, r2)
!SENote: This is similar to the version in adaptive_inflation. Should put in utilities.

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

!-----------------------------------------------------------------------


!-----------------------------------------------------------------------------------------------------
! End of routines for computing horizontal grid box location with cubed sphere spectral element grids 
!-----------------------------------------------------------------------------------------------------


!===================================================================
! End of model_mod
!===================================================================

end module model_mod
