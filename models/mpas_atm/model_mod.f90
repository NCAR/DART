! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module model_mod

! MPAS Atmosphere model interface to the DART data assimilation system.

! This revision of the model_mod supports both a global MPAS grid and
! a regional grid.  For the regional grid only observations which 
! are completely inside the interior will be assimilated, meaning obs
! which need interpolation information from the boundary cells 
! (in any of the 7 boundary layers) will be rejected.  However, during the
! assimilation phase all locations in the local grid will be impacted, 
! even locations in the boundary layers if there are obs close to the
! boundaries.  A post-processing step will smooth the GFS external
! values with the values updated by the assimilation in the boundary layers.

! Note that to reject obs during interpolation requires the model_interpolate()
! routine to check and return an error, but during the vertical conversion and
! get_close routines all state point operations must succeed, even those
! in the boundary layers.  Pay close attention to which internal routines are used
! by each to make sure the intended actions are what happens.

use        types_mod, only : r4, r8, i8, digits12, SECPERDAY, MISSING_R8,      &
                             rad2deg, deg2rad, PI, MISSING_I, obstypelength

use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, get_dist, query_location,          &
                             get_close_type, set_location, get_location,       &
                             write_location, vertical_localization_on,         &
                             VERTISUNDEF, VERTISSURFACE, VERTISLEVEL,          &
                             VERTISPRESSURE, VERTISHEIGHT, VERTISSCALEHEIGHT,  &
                             loc_get_close_obs => get_close_obs,               &
                             loc_get_close_state => get_close_state,           &
                             is_vertical, set_vertical_localization_coord
                             
use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, nc_check, &
                                 nc_begin_define_mode, nc_end_define_mode,                       &
                                 nc_open_file_readonly, nc_close_file,                           &
                                 nc_add_attribute_to_variable, nc_define_dimension,              &
                                 nc_define_unlimited_dimension, nc_define_character_variable,    &
                                 nc_define_real_variable, nc_get_variable, nc_put_variable,      &
                                 nc_get_dimension_size, nc_variable_exists, nc_dimension_exists, &
                                 nc_define_integer_variable

use location_io_mod,      only :  nc_write_location_atts, nc_write_location

use default_model_mod,     only : nc_write_model_vars, adv_1step,          &
                                  init_time => fail_init_time,             &
                                  init_conditions => fail_init_conditions

use xyz_location_mod, only : xyz_location_type, xyz_set_location, xyz_get_location,         &
                             xyz_get_close_type, xyz_get_close_init, xyz_get_close_destroy, &
                             xyz_find_nearest, xyz_use_great_circle_dist

use    utilities_mod, only : error_handler, get_unit,         &
                             E_ERR, E_WARN, E_MSG, E_ALLMSG, logfileunit,      &
                             to_upper, nmlfileunit,                 &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, file_exist, find_textfile_dims,        &
                             file_to_text, close_file, do_nml_file,            &
                             do_nml_term, scalar

use     obs_kind_mod, only : get_index_for_quantity,       &
                             get_name_for_quantity,        &
                             get_quantity_for_type_of_obs, &
                             QTY_SURFACE_ELEVATION,        &
                             QTY_SURFACE_PRESSURE,         &
                             QTY_10M_U_WIND_COMPONENT,     &
                             QTY_10M_V_WIND_COMPONENT,     &
                             QTY_2M_TEMPERATURE,           &
                             QTY_2M_SPECIFIC_HUMIDITY,     &
                             QTY_VERTICAL_VELOCITY,        &
                             QTY_POTENTIAL_TEMPERATURE,    &
                             QTY_EDGE_NORMAL_SPEED,        &
                             QTY_TEMPERATURE,              &
                             QTY_U_WIND_COMPONENT,         &
                             QTY_V_WIND_COMPONENT,         &
                             QTY_PRESSURE,                 &
                             QTY_DENSITY,                  &
                             QTY_VAPOR_MIXING_RATIO,       &
                             QTY_CLOUDWATER_MIXING_RATIO,  &
                             QTY_RAINWATER_MIXING_RATIO,   &
                             QTY_ICE_MIXING_RATIO,         &
                             QTY_SNOW_MIXING_RATIO,        &
                             QTY_GRAUPEL_MIXING_RATIO,     &
                             QTY_SPECIFIC_HUMIDITY,        &
                             QTY_GEOPOTENTIAL_HEIGHT,      &
                             QTY_PRECIPITABLE_WATER,       &
                             QTY_SKIN_TEMPERATURE,         &  ! for rttov
                             QTY_SURFACE_TYPE,             &  ! for rttov
                             QTY_CLOUD_FRACTION               ! for rttov

use mpi_utilities_mod, only: my_task_id, broadcast_minmax

use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use ensemble_manager_mod, only : ensemble_type, get_my_num_vars, get_my_vars

use distributed_state_mod, only : get_state, get_state_array

use netcdf

use state_structure_mod, only :  add_domain, get_model_variable_indices, &
                                 state_structure_info, get_index_start, get_index_end, &
                                 get_num_variables, get_domain_size, get_varid_from_varname, &
                                 get_variable_name, get_num_dims, get_dim_lengths, &
                                 get_dart_vector_index, get_num_varids_from_kind, &
                                 get_dim_name, get_varid_from_kind

use get_reconstruct_mod, only : get_reconstruct_init, get_reconstruct
use get_geometry_mod,    only : get_geometry


implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.

! routines in this list have code in this module
public :: get_model_size,                 &
          get_state_meta_data,            &
          model_interpolate,              &
          shortest_time_between_assimilations, &
          static_init_model,              &
          end_model,                      &
          nc_write_model_atts,            &
          pert_model_copies,              &
          get_close_obs,                  &
          get_close_state,                &
          convert_vertical_obs,           &
          convert_vertical_state,         &
          read_model_time,                &
          write_model_time

! code for these routines are in other modules
public :: init_time,                      &
          init_conditions,                &
          adv_1step,                      &
          nc_write_model_vars

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.

public :: get_init_template_filename,   &
          get_analysis_time,            &
          get_grid_dims,                &
          get_xland,                    &
          get_surftype,                 &
          get_cell_center_coords,       &
          get_bdy_mask,                 &
          find_closest_cell_center,     &
          find_triangle,                &
          read_2d_from_nc_file,         &
          find_height_bounds,           &
          cell_ok_to_interpolate,       &
          is_global_grid,               &
          uv_cell_to_edges

! try adjusting what static_init_model does, before it is called.
! the main update_bc program, for example, can call these before
! calling static_init_model() and that will modify what it does.
public :: set_lbc_variables, &
          force_u_into_state

! set_lbc_variables sets the lbc_variables string array
! force_u_into_state sets a logical add_u_to_state_list that forces u to be in state

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'models/mpas_atm/model_mod.f90'

! error codes: 
integer, parameter :: GENERAL_RTTOV_ERROR = 1 ! general error for rttov - at least one of the input variables goes wrong
integer, parameter :: CRITICAL_ERROR = 99 ! general error in case something terrible goes wrong
integer, parameter :: VERTICAL_TOO_HIGH = 81 ! Vertical location too high
integer, parameter :: VERTICAL_TOO_LOW = 80 ! Vertical location too low
integer, parameter :: QTY_NOT_IN_STATE_VECTOR = 88 ! qty is not in the state vector
integer, parameter :: TK_COMPUTATION_ERROR = 89 ! tk cannot be computed
integer, parameter :: CELL_CENTER_NOT_FOUND = 11 ! Could not find the closest cell center that contains this lat/lon
integer, parameter :: SURFACE_OBS_TOO_FAR = 12 ! Surface obs too far away from model elevation
integer, parameter :: INTERPOLATION_MISSING_VALUE = 13 ! Missing value in interpolation
integer, parameter :: TRIANGLE_CELL_CENTER_NOT_FOUND = 14 ! Could not find the other two cell centers of the triangle that contains this lat/lon
integer, parameter :: TRIANGLE_CELL_CENTER_IN_BOUNDARY = 15 ! Cell centers of the triangle fall in the lateral boundary zone
integer, parameter :: VERTICAL_VELOCITY_NOT_AVAIL = 16 ! Do not know how to do vertical velocity for now
integer, parameter :: PRESSURE_COMPUTATION_ERROR = 17 ! Unable to compute pressure values
integer, parameter :: ILLEGAL_ALTITUDE = 18
integer, parameter :: RBF_U_COMPUTATION_ERROR = 19 ! could not compute u using RBF code
integer, parameter :: INTERNAL_ERROR = 101 ! reached end of subroutine without finding an applicable case.
integer, parameter :: REJECT_OBS_USER_PRESSURE_LEVEL = 201 ! Reject observation from user specified pressure level
integer, parameter :: PRESSURE_NOT_MONOTONIC = 988 ! Pressure is not monotonically descreased with level


! module global storage; maintains values between calls, accessible by
! any subroutine
character(len=512) :: string1, string2, string3
logical, save :: module_initialized = .false.

! length of an mpas (also wrf) time string:  YYYY-MM-DD_hh:mm:ss
integer, parameter :: TIMELEN = 19

! Real (physical) constants as defined exactly in MPAS.
! redefined here for consistency with the model.
real(r8), parameter :: rgas = 287.0_r8  ! Constant: Gas constant for dry air [J kg-1 K-1]
real(r8), parameter :: rv = 461.6_r8    ! Constant: Gas constant for water vapor [J kg-1 K-1]
real(r8), parameter :: cp = 7.0_r8*rgas/2.0_r8 ! = 1004.5 Constant: Specific heat of dry air at constant pressure [J kg-1 K-1] 
real(r8), parameter :: cv = cp - rgas          ! = 717.5  Constant: Specific heat of dry air at constant volume [J kg-1 K-1]
real(r8), parameter :: p0 = 100000.0_r8
real(r8), parameter :: rcv = rgas/(cp-rgas)
real(r8), parameter :: rvord = rv/rgas    

! earth radius; needed to convert lat/lon to x,y,z cartesian coords.
! for the highest accuracy this should match what the model uses.
real(r8), parameter :: radius = 6371229.0_r8 ! meters

! roundoff error tolerance
! set in static_init_model() to 1e-5 or 1e-12 depending on compiled precision
real(r8) :: roundoff

! Storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq

! Structure for computing distances to cell centers, and assorted arrays
! needed for the get_close code.
type(xyz_get_close_type)             :: cc_gc
type(xyz_location_type), allocatable :: cell_locs(:)
logical :: search_initialized = .false.

! compile-time control over whether grid information is written to the
! diagnostic files or not.  if it is, the files are self-contained (e.g. for
! ease of plotting), but are also much larger than they would be otherwise.
! change this private variable to control whether the grid information is
! written or not.
logical :: add_static_data_to_diags = .false.

! variables which are in the module namelist
character(len=256) :: init_template_filename = 'mpas_init.nc'
character(len=256) :: bdy_template_filename = ''
integer            :: vert_localization_coord = VERTISHEIGHT
integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 21600
real(r8)           :: model_perturbation_amplitude = 0.0001   ! tiny amounts
logical            :: log_p_vert_interp = .true.     ! if true, interpolate vertical pressure in log space
character(len=32)  :: calendar = 'Gregorian'
real(r8)           :: highest_obs_pressure_mb   = 100.0_r8    ! do not assimilate obs higher than this level.
real(r8)           :: sfc_elev_max_diff = -1.0_r8    ! do not assimilate if |model - station| height is larger than this [m].

! this is not in the namelist or supported generally.
! (setting this to true avoids the surface elevation max diff 
! test for elevation and surface pressure.)
logical :: always_assim_surf_altimeters = .false.

integer :: anl_domid = -1 ! For state_structure_mod access
integer :: lbc_domid = -1

! if .false. use U/V reconstructed winds tri interp at centers for wind forward ops
! if .true.  use edge normal winds (u) with RBF functs for wind forward ops
logical :: use_u_for_wind = .false.

! if using rbf, options 1,2,3 control how many points go into the rbf.
! larger numbers use more points from further away
integer :: use_rbf_option = 2

! if .false. edge normal winds (u) should be in state vector and are written out directly
! if .true.  edge normal winds (u) are updated based on U/V reconstructed winds
logical :: update_u_from_reconstruct = .true.

! only if update_u_from_reconstruct is true,
! if .false. use the cell center u,v reconstructed values to update edge winds
! if .true., read in the original u,v winds, compute the increments after the
! assimilation, and use only the increments to update the edge winds
logical :: use_increments_for_u_update = .true.

! if set by: call force_u_into_state()  *BEFORE* calling static_init_model(),
! the code will add U (edge normal winds) to the mpas state vector even if it
! isn't in the namelist.
logical :: add_u_to_state_list = .false.

! if > 0, amount of distance in fractional model levels that a vertical
! point can be above or below the top or bottom of the grid and still be
! evaluated without error. 
real(r8) :: outside_grid_level_tolerance = -1.0_r8

! if the calling code updates an existing file it simply writes the state variable
! data.  if it needs to create a file from scratch it calls nc_write_model_atts()
! to optionally add grid info or any other non-state variables or attributes.
! setting this to .true. adds the mpas grid info to the file.  .false. does
! not and results in smaller diag/output files.
logical :: write_grid_to_diag_files = .false. 

! in converting to scale height for the vertical, set this to .false. to
! use simply the log of the pressure.  to normalize by the surface pressure
! (backwards compatible with previous code), set this to .true.
logical :: no_normalization_of_scale_heights = .true. 

! for regional MPAS
real(r8) :: dxmax  ! max distance between two adjacent cell centers in the mesh (in meters)

! when updating boundary files for regional mpas, note whether the boundary
! file has the reconstructed winds (lbc_ur, lbc_vr) or not. (this is set by
! looking at the bdy template file and seeing if those variables are there.)  
! if not, the other two options are ignored.  
! if they are in the lbc file, then the other logicals control whether to use 
! them instead of updating the U edge winds directly, and whether to use the 
! reconstructed increments or not. 
! the latter two options could be added to the namelist if someone wanted to
! explore the options for how the edge winds are updated in the boundary file.
! for now they're not - they're hardcoded true.
! note that these are for the boundary file update only - there are separate
! options for how to update the U winds in the main assimilation state.

logical :: lbc_file_has_reconstructed_winds     = .false.  


namelist /model_nml/             &
   init_template_filename,       &
   vert_localization_coord,      &
   assimilation_period_days,     &
   assimilation_period_seconds,  &
   model_perturbation_amplitude, &
   log_p_vert_interp,            &
   calendar,                     &
   use_u_for_wind,               &
   use_rbf_option,               &
   update_u_from_reconstruct,    &
   use_increments_for_u_update,  &
   highest_obs_pressure_mb,      &
   outside_grid_level_tolerance, &
   sfc_elev_max_diff,            &
   write_grid_to_diag_files,     &
   no_normalization_of_scale_heights

! DART state vector contents are specified in the input.nml:&mpas_vars_nml namelist.
integer, parameter :: MAX_STATE_VARIABLES = 80
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 2
integer, parameter :: NUM_BOUNDS_TABLE_COLUMNS = 4 !HK @todo get rid of clamp or fail
character(len=NF90_MAX_NAME) :: mpas_state_variables(MAX_STATE_VARIABLES * NUM_STATE_TABLE_COLUMNS ) = ' '
character(len=NF90_MAX_NAME) :: mpas_state_bounds(NUM_BOUNDS_TABLE_COLUMNS, MAX_STATE_VARIABLES ) = ' '
character(len=NF90_MAX_NAME) :: variable_table(MAX_STATE_VARIABLES, NUM_STATE_TABLE_COLUMNS )

! this is special and not in the namelist.  boundary files have a fixed
! set of variables with fixed names.
character(len=NF90_MAX_NAME) :: lbc_variables(MAX_STATE_VARIABLES) = ''

namelist /mpas_vars_nml/ mpas_state_variables, mpas_state_bounds

integer :: nfields

! Grid parameters - the values will be read from an mpas analysis file.

integer :: nCells        = -1  ! Total number of cells making up the grid
integer :: nVertices     = -1  ! Unique points in grid that are corners of cells
integer :: nEdges        = -1  ! Straight lines between vertices making up cells
integer :: maxEdges      = -1  ! Largest number of edges a cell can have
integer :: nVertLevels   = -1  ! Vertical levels; count of vert cell centers
integer :: nVertLevelsP1 = -1  ! Vert levels plus 1; count of vert cell faces
integer :: vertexDegree  = -1  ! Max number of cells/edges that touch any vertex
integer :: nSoilLevels   = -1  ! Number of soil layers

! scalar grid positions

! FIXME: we read in a lot of metadata about the grids.  if space becomes an
! issue we could consider reading in only the x,y,z arrays for all the items
! plus the radius, and then compute the lat/lon for locations needed by
! get_state_meta_data() on demand.  most of the computations we need to do
! are actually easier in xyz coords (no issue with poles).

! FIXME: it may be desirable to read in xCell(:), yCell(:), zCell(:)
! to keep from having to compute them on demand, especially since we
! have converted the radian lat/lon of the cell centers into degrees.
! we have to convert back, then take a few sin and cos to get xyz.
! time/space/accuracy tradeoff here.

real(r8), allocatable :: xVertex(:), yVertex(:), zVertex(:)
real(r8), allocatable :: xEdge(:), yEdge(:), zEdge(:)
real(r8), allocatable :: lonEdge(:) ! edge longitudes (degrees, original radians in file)
real(r8), allocatable :: latEdge(:) ! edge longitudes (degrees, original radians in file)
real(r8), allocatable :: lonCell(:) ! cell center longitudes (degrees, original radians in file)
real(r8), allocatable :: latCell(:) ! cell center latitudes  (degrees, original radians in file)
real(r8), allocatable :: dcEdge(:)  ! distance between two adjacent cell centers (in meters)
real(r8), allocatable :: xland(:)   ! land-ocean mask (1=land including sea-ice ; 2=ocean)  
real(r8), allocatable :: seaice(:)  ! sea-ice flag (0=no seaice; =1 otherwise) - for rttov
real(r8), allocatable :: skintemp(:)! ground or water surface temperature      - for rttov
real(r8), allocatable :: zGridFace(:,:)   ! geometric height at cell faces   (nVertLevelsP1,nCells)
real(r8), allocatable :: zGridCenter(:,:) ! geometric height at cell centers (nVertLevels,  nCells)
real(r8), allocatable :: zGridEdge(:,:)   ! geometric height at edge centers (nVertLevels,  nEdges)
!real(r8), allocatable :: zEdgeFace(:,:)   ! geometric height at edges faces  (nVertLevelsP1,nEdges)
!real(r8), allocatable :: zEdgeCenter(:,:) ! geometric height at edges faces  (nVertLevels  ,nEdges)

integer,  allocatable :: cellsOnVertex(:,:) ! list of cell centers defining a triangle
integer,  allocatable :: verticesOnCell(:,:)

integer,  allocatable :: edgesOnCell(:,:) ! list of edges that bound each cell
integer,  allocatable :: cellsOnEdge(:,:) ! list of cells that bound each edge
integer,  allocatable :: nedgesOnCell(:) ! list of edges that bound each cell
real(r8), allocatable :: edgeNormalVectors(:,:)

! Boundary information might be needed ... regional configuration?
! Read if available.

integer,  allocatable :: bdyMaskCell(:)
integer,  allocatable :: bdyMaskEdge(:)
integer,  allocatable :: maxLevelCell(:)

integer,  allocatable :: surftype(:)   !  ! surface type (land=0, water=1, seaice = 2) - for rttov

integer         :: model_size          ! the state vector length
type(time_type) :: model_timestep      ! smallest time to adv model

! useful flags in making decisions when searching for points, etc
logical :: global_grid = .true.        ! true = the grid covers the sphere with no holes
logical :: all_levels_exist_everywhere = .true. ! true = cells defined at all levels
logical :: has_edge_u = .false.        ! true = has original normal u on edges
logical :: has_uvreconstruct = .false. ! true = has reconstructed at centers

! Do we have any state vector items located on the cell edges?
! If not, avoid reading in or using the edge arrays to save space.
! FIXME: this should be set after looking at the fields listed in the
! namelist which are to be read into the state vector - if any of them
! are located on the edges then this flag should be changed to .true.
! however, the way the code is structured these arrays are allocated
! before the details of field list is examined.  since right now the
! only possible field array that is on the edges is the 'u' edge normal
! winds, search specifically for that in the state field list and set
! this based on that.  if any other data might be on edges, this routine
! will need to be updated: is_edgedata_in_state_vector()
logical :: data_on_edges = .false.

! currently unused; for a regional model it is going to be necessary to know
! if the grid is continuous around longitudes (wraps in east-west) or not,
! and if it covers either of the poles.
character(len= 64) :: ew_boundary_type, ns_boundary_type


INTERFACE get_analysis_time
      MODULE PROCEDURE get_analysis_time_ncid
      MODULE PROCEDURE get_analysis_time_fname
END INTERFACE

!------------------------------------------------

! The regular grid used for triangle interpolation divides the sphere into
! a set of regularly spaced lon-lat boxes. The number of boxes in
! longitude and latitude are set by num_reg_x and num_reg_y. Making the
! number of regular boxes smaller decreases the computation required for
! doing each interpolation but increases the static storage requirements
! and the initialization computation (which seems to be pretty small).
integer, parameter :: num_reg_x = 90, num_reg_y = 90

! The max_reg_list_num controls the size of temporary storage used for
! initializing the regular grid. Two arrays
! of size num_reg_x*num_reg_y*max_reg_list_num are needed. The initialization
! fails and returns an error if max_reg_list_num is too small. A value of
! ??? is sufficient for ???
integer, parameter :: max_reg_list_num = 100

! The triangle interpolation keeps a list of how many and which triangles
! overlap each regular lon-lat box. The number is stored in
! array triangle_num. The allocatable array
! triangle_list lists the uniquen index
! of each overlapping triangle. The entry in
! triangle_start for a given regular lon-lat box indicates
! where the list of triangles begins in the triangle_list.

integer :: triangle_start(num_reg_x, num_reg_y)
integer :: triangle_num  (num_reg_x, num_reg_y) = 0
integer, allocatable :: triangle_list(:)


contains

!------------------------------------------------------------------

subroutine static_init_model()
!>@todo FIXME - can we add an optional arg here for update_bc use? !HK No.

! Called to do one time initialization of the model.

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname,dimname
character(len=obstypelength)       :: kind_string
integer :: ncid, VarID, numdims, varsize, dimlen
integer :: iunit, io, ivar, i
integer :: ss, dd, z1, m1
integer :: nDimensions, nVariables, nAttributes, unlimitedDimID, TimeDimID
integer :: lbc_nfields
logical :: both
real(r8) :: variable_bounds(MAX_STATE_VARIABLES, 2)
integer :: variable_qtys(MAX_STATE_VARIABLES)
character(len=*), parameter :: routine = 'static_init_model'

if ( module_initialized ) return ! only need to do this once.

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Read the MPAS variable list to populate DART state vector
! Intentionally do not try to dump them to the nml unit because
! they include large character arrays which output pages of white space.
! The routine that reads and parses this namelist will output what
! values it found into the log.
call find_namelist_in_file('input.nml', 'mpas_vars_nml', iunit)
read(iunit, nml = mpas_vars_nml, iostat = io)
call check_namelist_read(iunit, io, 'mpas_vars_nml')

!---------------------------------------------------------------
! Set the time step ... causes mpas namelists to be read - HK does it?
! Ensures model_timestep is multiple of 'dynamics_timestep'

call set_calendar_type( calendar )   

model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation window is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,routine,string1,source)

ncid = nc_open_file_readonly(init_template_filename, routine)

call verify_state_variables(ncid, init_template_filename, &
                            nfields, variable_qtys, variable_bounds)
call read_grid()

call nc_close_file(ncid, routine)

anl_domid = add_domain(init_template_filename, nfields,           &
                       var_names  = variable_table (1:nfields,1), &
                       kind_list  = variable_qtys(1:nfields),     &
                       clamp_vals = variable_bounds(1:nfields,:) )

model_size = get_domain_size(anl_domid)

lbc_nfields = 0

! if we have a lateral boundary file, add it to the domain
! so we have access to the corresponding lbc_xxx fields.
!>@todo FIXME: if we want to do increments, we could also add a
! third domain which is the original forecast fields before
! the assimilation (so we can compute increments)  !HK why would you need to add a third domain?
if (.not. global_grid .and. lbc_variables(1) /= '') then
   ! regional: count number of lbc fields to read in
   COUNTUP: do i=1, MAX_STATE_VARIABLES
      if (lbc_variables(i) /= '') then
         lbc_nfields = lbc_nfields + 1
      else
         exit COUNTUP
      endif
   enddo COUNTUP
   lbc_domid = add_domain(bdy_template_filename, lbc_nfields,    &
                          var_names  = lbc_variables)
                        ! FIXME clamp_vals = variable_bounds(1:nfields,:) )
   model_size = model_size + get_domain_size(lbc_domid)
else
   lbc_domid = -1
endif

! do some sanity checking here:

! if you have at least one of these wind components in the state vector,
! you have to have them both.  the subroutine will error out if only one
! is found and not both.
if (has_uvreconstruct) then
   call winds_present(z1,m1,both)
endif

! if you set the namelist to use the reconstructed cell center winds,
! they have to be in the state vector.
if (update_u_from_reconstruct .and. .not. has_uvreconstruct) then
   write(string1,*) 'update_u_from_reconstruct cannot be True'
   write(string2,*) 'because state vector does not contain U/V reconstructed winds'
   call error_handler(E_ERR,'static_init_model',string1,source, &
                      text2=string2)
endif

! if you set the namelist to update the edge normal winds based on the
! updated cell center winds, and you also have the edge normal winds in
! the state vector, warn that the edge normal winds will be ignored
! when going back to the mpas_analysis.nc file.  not an error, but the
! updates to the edge normal vectors won't affect the results.
if (update_u_from_reconstruct .and. has_edge_u) then
   write(string1,*) 'edge normal winds (u) in MPAS file will be updated based on U/V reconstructed winds'
   write(string2,*) 'and not from the updated edge normal values in the state vector'
   write(string3,*) 'because update_u_from_reconstruct is True'
   call error_handler(E_MSG,'static_init_model',string1,source, &
                      text2=string2, text3=string3)
endif

! there are two choices when updating the edge normal winds based on the
! winds at the cell centers.  one is a direct interpolation of the values;
! the other is to read in the original cell centers, compute the increments
! changed by the interpolation, and then add or substract only the increments
! from the original edge normal wind values.
if (update_u_from_reconstruct) then
   if (use_increments_for_u_update) then
      write(string1,*) 'edge normal winds (u) in MPAS file will be updated based on the difference'
      write(string2,*) 'between the original U/V reconstructed winds and the updated values'
      write(string3,*) 'because use_increment_for_u_update is True'
   else
      write(string1,*) 'edge normal winds (u) in MPAS file will be updated by averaging the'
      write(string2,*) 'updated U/V reconstructed winds at the corresponding cell centers'
      write(string3,*) 'because use_increment_for_u_update is False'
   endif
   call error_handler(E_MSG,'static_init_model',string1,source, &
                      text2=string2, text3=string3)
endif

! basically we cannot do much without having at least these
! three fields in the state vector.  refuse to go further
! if these are not present:
if ((get_num_varids_from_kind(anl_domid, QTY_POTENTIAL_TEMPERATURE) < 0) .or. &
    (get_num_varids_from_kind(anl_domid, QTY_DENSITY) < 0) .or. &
    (get_num_varids_from_kind(anl_domid, QTY_VAPOR_MIXING_RATIO) < 0)) then
   write(string1, *) 'State vector is missing one or more of the following fields:'
   write(string2, *) 'Potential Temperature (theta), Density (rho), Vapor Mixing Ratio (qv).'
   write(string3, *) 'Cannot convert between height/pressure nor compute sensible temps.'
   call error_handler(E_ERR,'static_init_model',string1,source, &
                      text2=string2, text3=string3)
endif

! tell the location module how we want to localize in the vertical
call set_vertical_localization_coord(vert_localization_coord)

! set an appropriate value for roundoff tests based
! on this code being compiled single or double precision.
! set to 1e-5 (for single) or 1e-12 (for double precision).
if (r8 == digits12) then
   roundoff = 1.0e-12_r8
else
   roundoff = 1.0e-5_r8
endif

end subroutine static_init_model

!------------------------------------------------------------------
subroutine read_grid()

integer :: cel1, cel2
integer :: iloc, kloc
integer :: ncid
character(len=*), parameter :: routine = 'read_grid'

ncid = nc_open_file_readonly(init_template_filename, routine)

nCells        = nc_get_dimension_size(ncid, 'nCells',         routine)
nVertices     = nc_get_dimension_size(ncid, 'nVertices',      routine)
nEdges        = nc_get_dimension_size(ncid, 'nEdges',         routine)
maxEdges      = nc_get_dimension_size(ncid, 'maxEdges',       routine)
nVertLevels   = nc_get_dimension_size(ncid, 'nVertLevels',    routine)
nVertLevelsP1 = nc_get_dimension_size(ncid, 'nVertLevelsP1',  routine)
vertexDegree  = nc_get_dimension_size(ncid, 'vertexDegree',   routine)
nSoilLevels   = nc_get_dimension_size(ncid, 'nSoilLevels',    routine)

allocate(latCell(nCells), lonCell(nCells))
allocate(zGridFace(nVertLevelsP1, nCells))
allocate(zGridCenter(nVertLevels, nCells))

allocate(cellsOnVertex(vertexDegree, nVertices))
allocate(dcEdge(nEdges))
allocate(nEdgesOnCell(nCells))
allocate(xland(nCells))
allocate(seaice(nCells))              ! for rttov
allocate(skintemp(nCells))            ! for rttov
allocate(surftype(nCells))            ! for rttov
allocate(edgesOnCell(maxEdges, nCells))
allocate(cellsOnEdge(2, nEdges))
allocate(verticesOnCell(maxEdges, nCells))
allocate(edgeNormalVectors(3, nEdges))
allocate(xVertex(nVertices), yVertex(nVertices), zVertex(nVertices))

! see if U is in the state vector list.  if not, don't read in or
! use any of the Edge arrays to save space.
data_on_edges = is_edgedata_in_state_vector(variable_table, lbc_variables)

if(data_on_edges) then
   allocate(zGridEdge(nVertLevels, nEdges))
   allocate(xEdge(nEdges), yEdge(nEdges), zEdge(nEdges))
   allocate(latEdge(nEdges), lonEdge(nEdges))
endif

! is this a global or regional grid?
! if regional, allocate and read in the boundary info here.
call set_global_grid(ncid)

! fill in the grid values
call get_grid(ncid)

! vertical faces are in the input file.  compute vertical center locations here.
do kloc=1, nCells
   do iloc=1, nVertLevels
      zGridCenter(iloc,kloc) = (zGridFace(iloc,kloc) + zGridFace(iloc+1,kloc))*0.5_r8
   enddo
enddo

if(data_on_edges) then
   ! FIXME: This code is supposed to check whether an edge has 2 neighbours or 1 neighbour and then
   !        compute the height accordingly.  HOWEVER, the array cellsOnEdge does not change with
   !        depth, but it should as an edge may have 2 neighbour cells at the top but not at depth.
   do kloc=1, nEdges
      do iloc=1, nVertLevels
         cel1 = cellsOnEdge(1,kloc)
         cel2 = cellsOnEdge(2,kloc)
         if (cel1>0 .and. cel2>0) then
            zGridEdge(iloc,kloc) = (zGridCenter(iloc,cel1) + zGridCenter(iloc,cel2))*0.5_r8
         else if (cel1>0) then
            zGridEdge(iloc,kloc) = zGridCenter(iloc,cel1)
         else if (cel2>0) then
            zGridEdge(iloc,kloc) = zGridCenter(iloc,cel2)
         else  !this is bad...
            write(string1,*)'Edge ',kloc,' at vertlevel ',iloc,' has no neighbouring cells!'
            call error_handler(E_ERR,'static_init_model', string1, source)
         endif
      enddo
   enddo
endif

call nc_close_file(ncid, routine)


end subroutine read_grid


!------------------------------------------------------------------

subroutine get_state_meta_data(index_in, location, qty)

! Given an index into the state vector, return its location and
! optionally the qty

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: qty

integer  :: nzp, iloc, vloc, nf, ndim
integer  :: index, i, j, k, dom_id, var_id, qty_local, vert, cell, level
real(r8) :: lon, lat, height
type(location_type) :: new_location

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, i, j, k, dom_id=dom_id, var_id=var_id, kind_index=qty_local)

! Need to find if the variable is 2d (1d netcdf) or 3d (2d netcdf) 
! Need to find if variable is on cells or edges
! Need to find if variable is on nVertLevels or nVertLevelsP1

if (get_num_dims(dom_id, var_id) == 1) then  ! 2d (1d netcdf) variable: (cells)
   cell = i
   level = 1
   vert = VERTISSURFACE
   if (on_edge(dom_id, var_id)) then
      lon = lonEdge(cell)
      lat = latEdge(cell)
      height = zGridEdge(level, cell)
   else
      lon = lonCell(cell)
      lat = latCell(cell)
      height = zGridFace(level, cell) !HK 2D on face in original code
   endif

else ! 2d (1d netcdf) variable: (levels, cells)
   cell = j
   level = i
   vert = VERTISHEIGHT
   if (on_edge(dom_id, var_id)) then
      lon = lonEdge(cell)
      lat = latEdge(cell)
      height = zGridEdge(level, cell)
   else
      lon = lonCell(cell)
      lat = latCell(cell)
      if (get_dim_name(dom_id, var_id, 1) == 'nVertLevels') then
         height = zGridCenter(level, cell)
      else
         height = zGridFace(level, cell) 
      endif
   endif
endif

location = set_location(lon, lat, height, vert)

if (present(qty)) then
   qty = qty_local
endif

end subroutine get_state_meta_data

!------------------------------------------------------------------
function on_edge(dom, var)

integer, intent(in) :: dom, var
logical :: on_edge

if (get_num_dims(dom, var) == 1) then 
   on_edge = ( get_dim_name(dom, var, 1) == 'nEdges' ) 
else
   on_edge = ( get_dim_name(dom, var, 2) == 'nEdges' ) 
endif

end function on_edge


!------------------------------------------------------------------
!> given a state ensemble_handle, a location, and a QTY_xxx, return the
!> interpolated value at that location, and an error code.  0 is success,
!> anything positive is an error.  (negative reserved for system use)
subroutine model_interpolate(state_handle, ens_size, location, qty, expected_obs, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: qty
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)


type(location_type) :: location_tmp(ens_size)
integer  :: ivar, obs_kind
integer  :: tvars(3)
integer  :: cellid
logical  :: goodkind, surface_obs
real(r8) :: lpres(ens_size), values(3, ens_size)
real(r8) :: llv(3)    ! lon/lat/vert
integer  :: e, verttype
real(r8) :: single_expected_obs
integer  :: single_istatus

if ( .not. module_initialized ) call static_init_model

expected_obs = MISSING_R8
istatus      = INTERNAL_ERROR

llv = get_location(location)
verttype = nint(query_location(location))
surface_obs = (verttype == VERTISSURFACE) 

cellid = cell_ok_to_interpolate(location)
if (cellid < 1) then
   istatus = CELL_CENTER_NOT_FOUND
   return
endif

! HK @todo elevation check (but not for rttov)

! Explicit fail for qty_vertical_velocity
if (qty == QTY_VERTICAL_VELOCITY) then
   istatus = VERTICAL_VELOCITY_NOT_AVAIL
   return
endif

if (.not. qty_ok_to_interpolate(qty)) then
   istatus = QTY_NOT_IN_STATE_VECTOR
   return
endif

! HK @todo reject obs above a user specified pressure level

select case (qty)

   case (QTY_U_WIND_COMPONENT)
      if (use_u_for_wind .and. has_edge_u) then
         call compute_u_with_rbf(state_handle, ens_size, location, .TRUE., expected_obs, istatus)
      else
         tvars(1) = get_varid_from_kind(anl_domid, qty)
         call compute_scalar_with_barycentric(state_handle, ens_size, location, 1, tvars(1), expected_obs, istatus)
      endif
      
   case (QTY_V_WIND_COMPONENT)
      if (use_u_for_wind .and. has_edge_u) then
         call compute_u_with_rbf(state_handle, ens_size, location, .FALSE., expected_obs, istatus)
      else
         tvars(1) = get_varid_from_kind(anl_domid, qty)
         call compute_scalar_with_barycentric(state_handle, ens_size, location, 1, tvars(1), expected_obs, istatus)
      endif

   case (QTY_TEMPERATURE) 
      ! need to get potential temp, pressure, qv here, but can
      ! use same weights, so push all three types into the subroutine.
      tvars(1) = get_varid_from_kind(anl_domid, QTY_POTENTIAL_TEMPERATURE)
      tvars(2) = get_varid_from_kind(anl_domid, QTY_DENSITY)
      tvars(3) = get_varid_from_kind(anl_domid, QTY_VAPOR_MIXING_RATIO)

      call compute_scalar_with_barycentric(state_handle, ens_size, location, 3, tvars, values, istatus)
      if (fails(istatus)) then
         istatus = TK_COMPUTATION_ERROR
         return ! early return to not pass nonsense values to theta_to_tk
      endif
      ! convert pot_temp, density, vapor mixing ratio into sensible temperature
      expected_obs(:) = theta_to_tk(ens_size, values(1, :), values(2, :), values(3, :), istatus) 

   case (QTY_PRESSURE) 
      call compute_pressure_at_loc(state_handle, ens_size, location, expected_obs, istatus)
      if (fails(istatus)) istatus = PRESSURE_COMPUTATION_ERROR

   case (QTY_GEOPOTENTIAL_HEIGHT)
      location_tmp = location
      call convert_vert(state_handle, ens_size, location_tmp, QTY_GEOPOTENTIAL_HEIGHT, VERTISHEIGHT, istatus)

      do e = 1, ens_size
      if(istatus(e) == 0) expected_obs(e) = query_location(location_tmp(e), 'VLOC')
      enddo

   ! HK @todo things that cannot be negative
   case (QTY_VAPOR_MIXING_RATIO, &
         QTY_CLOUDWATER_MIXING_RATIO, &
         QTY_ICE_MIXING_RATIO, &
         QTY_GRAUPEL_MIXING_RATIO, &
         QTY_2M_SPECIFIC_HUMIDITY, &
         QTY_RAINWATER_MIXING_RATIO, &
         QTY_SNOW_MIXING_RATIO, &
         QTY_CLOUD_FRACTION  )

   case(QTY_SPECIFIC_HUMIDITY)
      tvars(1) = get_varid_from_kind(anl_domid, QTY_VAPOR_MIXING_RATIO)
      call compute_scalar_with_barycentric(state_handle, ens_size, location, 1, tvars(1), expected_obs, istatus)

      ! compute vapor pressure, then: sh = vp / (1.0 + vp)
      do e = 1, ens_size
         if (istatus(e) == 0) then
            if (expected_obs(e) >= 0.0_r8) then
               expected_obs(e) = expected_obs(e) / (1.0_r8 + expected_obs(e))
            else
               expected_obs(e) = 1.0e-12_r8  
            endif
         endif 
      enddo

   ! qtys not included in the dart state vector, but static across the ensemble
   case (QTY_SURFACE_ELEVATION)
      call compute_surface_data_with_barycentric(zGridFace(1,:), location, single_expected_obs, single_istatus)
      expected_obs(:) = single_expected_obs
      istatus(:) = single_istatus

   case (QTY_SKIN_TEMPERATURE)
      call compute_surface_data_with_barycentric(skintemp(:), location, single_expected_obs, single_istatus)
      expected_obs(:) = single_expected_obs
      istatus(:) = single_istatus

   case (QTY_SURFACE_TYPE)
      call get_surftype(nCells,surftype)
      call compute_surface_data_with_barycentric(surftype(:)*1.0_r8, location, single_expected_obs, single_istatus)
      expected_obs(:) = single_expected_obs
      istatus(:) = single_istatus
   
   case default ! regular interpolation
      tvars(1) = get_varid_from_kind(anl_domid, qty)
      call compute_scalar_with_barycentric(state_handle, ens_size, location, 1, tvars(1), expected_obs, istatus)

end select



end subroutine model_interpolate

!------------------------------------------------------------------
function fails(e) result(fail)

integer, intent(in) :: e(:)
logical :: fail

fail = (any(e /= 0))

end function fails

!------------------------------------------------------------------
function qty_ok_to_interpolate(qty) result(qty_ok)

integer, intent(in) :: qty
logical :: qty_ok

integer :: varid

qty_ok = .false.
varid = get_varid_from_kind(anl_domid, qty)

if (varid > 0) then ! in the state vector
   qty_ok = .true.
   return
endif

select case (qty)
case (QTY_TEMPERATURE, QTY_2M_TEMPERATURE)
   qty_ok = .true.
case (QTY_SURFACE_ELEVATION, QTY_GEOPOTENTIAL_HEIGHT)
   qty_ok = .true.
case (QTY_PRESSURE)   ! surface pressure should be in the state !HK @todo check "should"
   qty_ok = .true.
case (QTY_SKIN_TEMPERATURE, QTY_SURFACE_TYPE, QTY_CLOUD_FRACTION)   
   qty_ok = .true.
case (QTY_SPECIFIC_HUMIDITY, QTY_2M_SPECIFIC_HUMIDITY)
   qty_ok = .true.
case (QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT)
   if (get_varid_from_kind(anl_domid, QTY_EDGE_NORMAL_SPEED) > 0 .and. use_u_for_wind) qty_ok = .true.
case default
   qty_ok = .false.
end select

end function qty_ok_to_interpolate
!------------------------------------------------------------------

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid
integer, intent(in) :: domain_id

character(len=*), parameter :: routine = 'nc_write_model_atts'

real(r8), allocatable :: data1d(:)
integer :: ncid2


!-------------------------------------------------------------------------------
! put file into define mode.
!-------------------------------------------------------------------------------

call nc_begin_define_mode(ncid)

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source)
call nc_add_global_attribute(ncid, "model", "MPAS_ATM")


!----------------------------------------------------------------------------
! if not adding grid info, return here
!----------------------------------------------------------------------------

if (.not. write_grid_to_diag_files) then
   call nc_end_define_mode(ncid)
   call nc_synchronize_file(ncid)
   return
endif

!----------------------------------------------------------------------------
!  Everything below here is static grid info
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
! Dimensions 
!----------------------------------------------------------------------------

call nc_define_dimension(ncid, 'nCells',        nCells,        routine)
call nc_define_dimension(ncid, 'nEdges',        nEdges,        routine)
call nc_define_dimension(ncid, 'nVertLevels',   nVertLevels,   routine)
call nc_define_dimension(ncid, 'nVertLevelsP1', nVertLevelsP1, routine)
call nc_define_dimension(ncid, 'nSoilLevels',   nSoilLevels,   routine)

call nc_define_dimension(ncid, 'maxEdges',      maxEdges,     routine)
call nc_define_dimension(ncid, 'nVertices',     nVertices,    routine)
call nc_define_dimension(ncid, 'VertexDegree',  VertexDegree, routine)

!----------------------------------------------------------------------------
! Coordinate Variables and the Attributes
!----------------------------------------------------------------------------
! Cell Longitudes
call      nc_define_real_variable(ncid, 'lonCell', 'nCells', routine)
call nc_add_attribute_to_variable(ncid, 'lonCell', 'long_name',   'cell center longitudes', routine)
call nc_add_attribute_to_variable(ncid, 'lonCell', 'units',       'degrees_east',           routine)
call nc_add_attribute_to_variable(ncid, 'lonCell', 'valid_range', (/ 0.0_r8, 360.0_r8 /),   routine)

! Cell Latitudes
call      nc_define_real_variable(ncid, 'latCell', 'nCells', routine)
call nc_add_attribute_to_variable(ncid, 'latCell', 'long_name',   'cell center latitudes',  routine)
call nc_add_attribute_to_variable(ncid, 'latCell', 'units',       'degrees_north',          routine)
call nc_add_attribute_to_variable(ncid, 'latCell', 'valid_range', (/ -90.0_r8, 90.0_r8 /),  routine)

! Grid vertical information
call      nc_define_real_variable(ncid, 'zgrid', (/ 'nVertLevelsP1', 'nCells       ' /), routine)
call nc_add_attribute_to_variable(ncid, 'zgrid', 'long_name',        'grid zgrid',     routine)
call nc_add_attribute_to_variable(ncid, 'zgrid', 'units',            'meters',         routine)
call nc_add_attribute_to_variable(ncid, 'zgrid', 'positive',         'up',             routine)
call nc_add_attribute_to_variable(ncid, 'zgrid', 'cartesian_axis',   'Z',              routine)

! Grid relationship information
call   nc_define_integer_variable(ncid, 'nEdgesOnCell', 'nCells', routine)
call nc_add_attribute_to_variable(ncid, 'nEdgesOnCell', 'long_name', 'grid nEdgesOnCell', routine)

call   nc_define_integer_variable(ncid, 'cellsOnVertex', (/ 'VertexDegree', 'nVertices   ' /), routine)
call nc_add_attribute_to_variable(ncid, 'cellsOnVertex', 'long_name', 'grid cellsOnVertex', routine)

call   nc_define_integer_variable(ncid, 'verticesOnCell', (/ 'maxEdges', 'nCells  ' /), routine)
call nc_add_attribute_to_variable(ncid, 'verticesOnCell', 'long_name', 'grid verticesOnCell', routine)

! Cartesian coordinates for the same cells
call      nc_define_real_variable(ncid, 'xCell', 'nCells', routine)
call nc_add_attribute_to_variable(ncid, 'xCell', 'long_name', 'cell center x cartesian coordinates',   routine)

call      nc_define_real_variable(ncid, 'yCell', 'nCells', routine)
call nc_add_attribute_to_variable(ncid, 'yCell', 'long_name', 'cell center y cartesian coordinates',   routine)

call      nc_define_real_variable(ncid, 'zCell', 'nCells', routine)
call nc_add_attribute_to_variable(ncid, 'zCell', 'long_name', 'cell center z cartesian coordinates',   routine)

! Vertex Longitudes and latitudes
call      nc_define_real_variable(ncid, 'lonVertex', 'nVertices', routine)
call nc_add_attribute_to_variable(ncid, 'lonVertex', 'long_name', 'vertex longitudes', routine)

call      nc_define_real_variable(ncid, 'latVertex', 'nVertices', routine)
call nc_add_attribute_to_variable(ncid, 'latVertex', 'long_name', 'vertex latitudes', routine)

! Edge Longitudes and latitudes
if(data_on_edges) then
   call      nc_define_real_variable(ncid, 'lonEdge', 'nEdges', routine)
   call nc_add_attribute_to_variable(ncid, 'lonEdge', 'long_name', 'edge longitudes', routine)

   call      nc_define_real_variable(ncid, 'latEdge', 'nEdges', routine)
   call nc_add_attribute_to_variable(ncid, 'latEdge', 'long_name', 'edge latitudes', routine)
endif

call nc_define_real_variable(ncid, 'areaCell', 'nCells', routine)

!----------------------------------------------------------------------------
! Finished with dimension/variable definitions, must end 'define' mode to fill.
!----------------------------------------------------------------------------

call nc_end_define_mode(ncid)

!----------------------------------------------------------------------------
! Fill the coordinate variables
!----------------------------------------------------------------------------

call nc_put_variable(ncid, 'lonCell', lonCell, routine)
call nc_put_variable(ncid, 'latCell', latCell, routine)

call nc_put_variable(ncid, 'zgrid', zGridFace, routine)

if(data_on_edges) then
   call nc_put_variable(ncid, 'lonEdge', lonEdge, routine)
   call nc_put_variable(ncid, 'latEdge', latEdge, routine)
endif

call nc_put_variable(ncid, 'nEdgesOnCell', nEdgesOnCell, routine)
call nc_put_variable(ncid, 'verticesOnCell', verticesOnCell, routine)
call nc_put_variable(ncid, 'cellsOnVertex', cellsOnVertex, routine)

!----------------------------------------------------------------------------
! DART has not read these in, so we have to read them from the input file
! and copy them to the DART output file.
!----------------------------------------------------------------------------

ncid2 = nc_open_file_readonly(init_template_filename, routine)

allocate(data1d(nCells))

call nc_get_variable(ncid2, 'xCell', data1d, routine)
call nc_put_variable(ncid,  'xCell', data1d, routine)

call nc_get_variable(ncid2, 'yCell', data1d, routine)
call nc_put_variable(ncid,  'yCell', data1d, routine)

call nc_get_variable(ncid2, 'zCell', data1d, routine)
call nc_put_variable(ncid,  'zCell', data1d, routine)

call nc_get_variable(ncid2, 'areaCell', data1d, routine)
call nc_put_variable(ncid,  'areaCell', data1d, routine)

deallocate(data1d)

allocate(data1d(nVertices))

call nc_get_variable(ncid2, 'lonVertex', data1d, routine)
call nc_put_variable(ncid,  'lonVertex', data1d, routine)

call nc_get_variable(ncid2, 'latVertex', data1d, routine)
call nc_put_variable(ncid,  'latVertex', data1d, routine)

deallocate(data1d)

call nc_close_file(ncid2)

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!------------------------------------------------------------------

function get_model_size()

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size

!------------------------------------------------------------------
function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = model_timestep

end function shortest_time_between_assimilations

!------------------------------------------------------------------

subroutine end_model()

! Does any shutdown and clean-up needed for model.

if (allocated(latCell))        deallocate(latCell)
if (allocated(lonCell))        deallocate(lonCell)
if (allocated(zGridFace))      deallocate(zGridFace)
if (allocated(zGridCenter))    deallocate(zGridCenter)
if (allocated(dcEdge))         deallocate(dcEdge)
if (allocated(cellsOnVertex))  deallocate(cellsOnVertex)
if (allocated(xland))          deallocate(xland)
if (allocated(seaice))         deallocate(seaice)
if (allocated(skintemp))       deallocate(skintemp)
if (allocated(surftype))       deallocate(surftype)
if (allocated(nEdgesOnCell))   deallocate(nEdgesOnCell)
if (allocated(edgesOnCell))    deallocate(edgesOnCell)
if (allocated(cellsOnEdge))    deallocate(cellsOnEdge)
if (allocated(verticesOnCell)) deallocate(verticesOnCell)
if (allocated(edgeNormalVectors)) deallocate(edgeNormalVectors)
if (allocated(xVertex))        deallocate(xVertex)
if (allocated(yVertex))        deallocate(yVertex)
if (allocated(zVertex))        deallocate(zVertex)
if (allocated(zGridEdge))      deallocate(zGridEdge)
if (allocated(xEdge))          deallocate(xEdge)
if (allocated(yEdge))          deallocate(yEdge)
if (allocated(zEdge))          deallocate(zEdge)
if (allocated(latEdge))        deallocate(latEdge)
if (allocated(lonEdge))        deallocate(lonEdge)

end subroutine end_model

!------------------------------------------------------------------

subroutine pert_model_copies(ens_handle, ens_size, pert_amp, interf_provided)

 type(ensemble_type), intent(inout) :: ens_handle
 integer,                intent(in) :: ens_size
 real(r8),               intent(in) :: pert_amp
 logical,               intent(out) :: interf_provided

logical, allocatable  :: within_range(:)
real(r8), allocatable :: min_var(:), max_var(:)
integer  :: start_ind, end_ind
real(r8) :: pert_val, range
integer  :: copy
integer  :: num_variables
integer  :: i, j
integer(i8), allocatable :: var_list(:)


! Perturbs a model state copies for generating initial ensembles.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, 
! it may do so by adding a perturbation to each model state 
! variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.

interf_provided = .true.

num_variables = get_num_variables(anl_domid) !HK what about the other domain?

! Get min and max of each variable in each domain
allocate(var_list(get_my_num_vars(ens_handle)))
call get_my_vars(ens_handle, var_list)

allocate(min_var(num_variables), max_var(num_variables))
allocate(within_range(ens_handle%my_num_vars))

do i = 1, get_num_variables(anl_domid)

   start_ind = get_index_start(anl_domid, i)
   end_ind = get_index_end(anl_domid, i)

   within_range = (var_list >= start_ind .and. var_list <= end_ind)
   min_var(i) = minval(ens_handle%copies(1,:), MASK=within_range)
   max_var(i) = maxval(ens_handle%copies(1,:), MASK=within_range)

enddo

! get global min/max for each variable
call broadcast_minmax(min_var, max_var, num_variables)
deallocate(within_range)

call init_random_seq(random_seq, my_task_id()+1)

do i = 1, num_variables

   start_ind = get_index_start(anl_domid, i)
   end_ind = get_index_end(anl_domid, i)

   ! make the perturbation amplitude a fraction of the
   ! entire variable range.
   range = max_var(i) - min_var(i)
   pert_val = model_perturbation_amplitude * range   ! this is a namelist item

   do j=1, ens_handle%my_num_vars
      if (ens_handle%my_vars(j) >= start_ind .and. ens_handle%my_vars(j) <= end_ind) then
         do copy = 1, ens_size
            ens_handle%copies(copy, j) = random_gaussian(random_seq, ens_handle%copies(copy, j), pert_val)
         enddo

         ! keep variable from exceeding the original range
         ens_handle%copies(1:ens_size,j) = max(min_var(i), ens_handle%copies(1:ens_size,j))
         ens_handle%copies(1:ens_size,j) = min(max_var(i), ens_handle%copies(1:ens_size,j))

      endif
   enddo

enddo

deallocate(var_list, min_var, max_var)

end subroutine pert_model_copies

!------------------------------------------------------------------
! Given a DART location (referred to as "base") and a set of candidate
! locations & qtys/types (locs, loc_qtys/types), returns the subset close 
! to the "base", their indices, and their distances to the "base" 

subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, state_handle)

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:), loc_types(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: state_handle


integer                :: ztypeout
integer                :: t_ind, istatus(1), k, istat_arr(1)
integer                :: base_which, local_obs_which, base_qty
real(r8), dimension(3) :: base_llv, local_obs_llv   ! lon/lat/vert
type(location_type)    :: local_obs_loc, location_arr(1)

num_close = 0

! Convert base_obs vertical coordinate to requested vertical coordinate if necessary
base_llv = get_location(base_loc)
base_which = nint(query_location(base_loc))

ztypeout = vert_localization_coord

if (vertical_localization_on()) then
  if (base_llv(3) == MISSING_R8) then !HK @todo what about VERTISUNDEF?
     return
  else if (base_which /= vert_localization_coord .and. base_which /= VERTISUNDEF) then
      base_qty = get_quantity_for_type_of_obs(base_type)
      location_arr(1) = base_loc
      call convert_vert(state_handle, 1, location_arr, base_qty, vert_localization_coord, istat_arr)
      if (istat_arr(1) /= 0) return
      base_loc = location_arr(1)
   endif
endif

! Loop over potentially close subset of obs priors or state variables
call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, num_close, close_ind)

do k = 1, num_close

   t_ind = close_ind(k)
   local_obs_loc   = locs(t_ind)
   local_obs_which = nint(query_location(local_obs_loc))

   ! Convert local_obs vertical coordinate to requested vertical coordinate if necessary.
   if ((base_which /= VERTISUNDEF) .and.  (vertical_localization_on())) then
         if (local_obs_which /= vert_localization_coord .and. local_obs_which /= VERTISUNDEF) then
            location_arr(1) = local_obs_loc
            call convert_vert(state_handle, 1, location_arr, loc_qtys(t_ind), vert_localization_coord, istatus)
            locs(t_ind) = location_arr(1)
         else
            istatus = 0
         endif
   endif

   if (present(dist)) then
      if (istatus(1) == 0) then
         dist(k) = get_dist(base_loc, locs(t_ind), base_type, loc_qtys(t_ind))
      else
         dist(k) = 1.0e9_r8
      endif
   endif

enddo

end subroutine get_close_obs

!------------------------------------------------------------------
! Given a DART location (referred to as "base") and a set of candidate
! locations & qtys/indices (locs, loc_qtys/loc_indx), returns the subset close 
! to the "base", their indices, and their distances to the "base" 

subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, state_handle)

!>@todo FIXME this is working on state vector items.  if a vertical
!>conversion is needed, it doesn't need to interpolate.  it can compute
!>the location using the logic that get_state_meta_data() uses.

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:)
integer(i8),                   intent(in)  :: loc_indx(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: state_handle


integer                :: ztypeout
integer                :: t_ind, istatus1, istatus2, k, istat_arr(1)
integer                :: base_which, local_obs_which, base_qty
real(r8), dimension(3) :: base_llv, local_obs_llv   ! lon/lat/vert
type(location_type)    :: local_obs_loc, location_arr(1)
! timing
real(digits12) :: t_base, t_base2, interval


real(r8) ::  hor_dist
hor_dist = 1.0e9_r8

! Initialize variables to missing status

num_close = 0
istatus1  = 0
istatus2  = 0

! Convert base_obs vertical coordinate to requested vertical coordinate if necessary

base_llv = get_location(base_loc)
base_which = nint(query_location(base_loc))

ztypeout = vert_localization_coord

if (vertical_localization_on()) then
  if (base_llv(3) == MISSING_R8) then
     istatus1 = 1
  else if (base_which /= vert_localization_coord .and. base_which /= VERTISUNDEF) then
      base_qty = get_quantity_for_type_of_obs(base_type)
      location_arr(1) = base_loc
      call convert_vert(state_handle, 1, location_arr, base_qty, vert_localization_coord, istat_arr)
      istatus1 = istat_arr(1)
      base_loc = location_arr(1)
   endif
endif

if (istatus1 == 0) then

   ! Loop over potentially close subset of obs priors or state variables
   ! This way, we are decreasing the number of distance computations that will follow.
   ! This is a horizontal-distance operation and we don't need to have the relevant vertical
   ! coordinate information yet (for locs).
   call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                            num_close, close_ind)

   do k = 1, num_close

      t_ind = close_ind(k)
      local_obs_loc   = locs(t_ind)
      local_obs_which = nint(query_location(local_obs_loc))

      ! Convert local_obs vertical coordinate to requested vertical coordinate if necessary.
      ! This should only be necessary for obs priors, as state location information already
      ! contains the correct vertical coordinate (filter_assim's call to get_state_meta_data).
      if ((base_which /= VERTISUNDEF) .and.  (vertical_localization_on())) then
          if (local_obs_which /= vert_localization_coord .and. local_obs_which /= VERTISUNDEF) then
              location_arr(1) = local_obs_loc
              call convert_vert_state(state_handle, 1, location_arr, loc_qtys(t_ind), &
                                              loc_indx(t_ind), vert_localization_coord, istat_arr)
              istatus2 = istat_arr(1)
              locs(t_ind) = location_arr(1)
          else
              istatus2 = 0
          endif
      endif

      if (present(dist)) then
         ! Compute distance - set distance to a very large value if vert coordinate is missing
         ! or vert_interpolate returned error (istatus2=1)
         local_obs_llv = get_location(local_obs_loc)
         if ( (vertical_localization_on() .and. &
              (local_obs_llv(3) == MISSING_R8)) .or. (istatus2 /= 0) ) then
               dist(k) = 1.0e9_r8
         else
               dist(k) = get_dist(base_loc, locs(t_ind), base_type, loc_qtys(t_ind))
         endif
      endif

   enddo
endif

end subroutine get_close_state


!==================================================================
! The (model-specific) additional public interfaces come next
!  (these are not required by dart but are used by other programs)
!==================================================================

subroutine get_init_template_filename( filename )  !HK why n

! return the name of the template filename that was set
! in the model_nml namelist (ex. init.nc)

character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(init_template_filename)

end subroutine get_init_template_filename


!-------------------------------------------------------------------
! modify what static_init_model does.  this *must* be called before !HK Nope this is not good.
! calling static_init_model().
! the boundary file variables are fixed by the model and so we
! don't allow the user to set them via namelist

subroutine set_lbc_variables(template_filename)

character(len=*), intent(in) :: template_filename

integer :: ncid

bdy_template_filename = template_filename

! this initial list always exists.  hardcode them for now,
! and query to see if the reconstructed winds are there or not.

lbc_variables(1) = 'lbc_qc'
lbc_variables(2) = 'lbc_qr'
lbc_variables(3) = 'lbc_qv'
lbc_variables(4) = 'lbc_rho'
lbc_variables(5) = 'lbc_theta'
lbc_variables(6) = 'lbc_u'
lbc_variables(7) = 'lbc_w'

ncid = nc_open_file_readonly(template_filename, 'set_lbc_variables')
if (nc_variable_exists(ncid, 'lbc_ur')) then
   lbc_variables(8) = 'lbc_ur'  
   lbc_variables(9) = 'lbc_vr' 
   lbc_file_has_reconstructed_winds = .true.
endif
call nc_close_file(ncid)

end subroutine set_lbc_variables


!-------------------------------------------------------------------
! modify what static_init_model does.  this *must* be called before
! calling static_init_model().
! set a logical add_u_to_state_list that forces u to be in state

subroutine force_u_into_state()

add_u_to_state_list = .true.

end subroutine force_u_into_state


!------------------------------------------------------------------

function get_analysis_time_ncid( ncid, filename )

! The analysis netcdf files have the start time of the experiment.
! The time array contains the time trajectory since then.
! This routine returns the start time of the experiment.

integer,          intent(in) :: ncid
character(len=*), intent(in) :: filename
type(time_type) :: get_analysis_time_ncid

! local variables
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, idims
integer           :: VarID, numdims

character(len=64) :: timestring

if ( .not. module_initialized ) call static_init_model

call nc_check( nf90_inq_varid(ncid, 'xtime', VarID), &
              'get_analysis_time', 'inquire xtime '//trim(filename))

call nc_check( nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
              'get_analysis_time', 'inquire TIME '//trim(filename))

if (numdims /= 2) then
   write(string1,*) 'xtime variable has unknown shape in ', trim(filename)
   call error_handler(E_ERR,'get_analysis_time',string1,source)
endif

call nc_check( nf90_inquire_dimension(ncid, dimIDs(1), len=idims(1)), &
                 'get_analysis_time', 'inquire time dimension length '//trim(filename))
call nc_check( nf90_inquire_dimension(ncid, dimIDs(2), len=idims(2)), &
                 'get_analysis_time', 'inquire time dimension length '//trim(filename))

if (idims(2) /= 1) then
   write(string1,*) 'multiple timesteps (',idims(2),') in file ', trim(filename)
   write(string2,*) 'We are using the LAST one, presumably, the LATEST timestep.'
   call error_handler(E_MSG,'get_analysis_time',string1,source,text2=string2)
endif

! Get the highest ranking time ... the last one, basically.

call nc_check( nf90_get_var(ncid, VarID, timestring, start = (/ 1, idims(2) /)), &
              'get_analysis_time', 'get_var xtime '//trim(filename))

get_analysis_time_ncid = string_to_time(timestring)

end function get_analysis_time_ncid


!------------------------------------------------------------------

function get_analysis_time_fname(filename)

! The analysis netcdf files have the start time of the experiment.
! The time array contains the time trajectory since then.
! This routine returns the start time of the experiment.

type(time_type) :: get_analysis_time_fname

character(len=*), intent(in) :: filename

integer :: i

if ( .not. module_initialized ) call static_init_model

! why do we care?  we aren't opening this file, just taking the time
! string from the filename itself.
if ( .not. file_exist(filename) ) then
   write(string1,*) 'file ', trim(filename),' does not exist.' 
   call error_handler(E_ERR,'get_analysis_time',string1,source)
endif

! find the first digit and use that as the start of the string conversion
i = scan(filename, "0123456789")
if (i <= 0) then
   write(string1,*) 'cannot find time string in name ', trim(filename)
   call error_handler(E_ERR,'get_analysis_time',string1,source)
endif

get_analysis_time_fname = string_to_time(filename(i:i+TIMELEN-1))

end function get_analysis_time_fname


!-----------------------------------------------------------------------
subroutine write_model_time(ncid, dart_time)

integer,             intent(in) :: ncid !< netcdf file handle
type(time_type),     intent(in) :: dart_time

integer :: year, month, day, hour, minute, second
character(len=64) :: timestring
character(len=*), parameter :: routine = 'write_model_time'

timestring = ''
call get_date(dart_time, year, month, day, hour, minute, second)
call set_wrf_date(timestring, year, month, day, hour, minute, second)

! Define xtime variable if it does not already exist
if (.not. nc_variable_exists(ncid, 'xtime')) then

   call nc_begin_define_mode(ncid)

   ! check to see if there are Time and date_string_length dimensions
   if (.not. nc_dimension_exists(ncid, 'Time')) &
      call nc_define_unlimited_dimension(ncid, 'Time', routine)

   if (.not. nc_dimension_exists(ncid, 'StrLen')) &
      call nc_define_dimension(ncid, 'StrLen', len(timestring), routine)

   ! make xtime(Time, StrLen)
   call nc_define_character_variable(ncid, 'xtime', (/ 'StrLen', 'Time  ' /), routine)
   call nc_add_attribute_to_variable(ncid, 'xtime', 'units', "YYYY-MM-DD_hh:mm:ss", routine)
   call nc_add_attribute_to_variable(ncid, 'xtime', 'long_name', 'Model valid time', routine)

call nc_end_define_mode(ncid)

endif

call nc_put_variable(ncid, 'xtime', timestring, routine)

end subroutine write_model_time

!------------------------------------------------------------------

subroutine get_grid_dims(Cells, Vertices, Edges, VertLevels, VertexDeg, SoilLevels)

! public routine for returning the counts of various things in the grid
!

integer, intent(out) :: Cells         ! Total number of cells making up the grid
integer, intent(out) :: Vertices      ! Unique points in grid which are corners of cells
integer, intent(out) :: Edges         ! Straight lines between vertices making up cells
integer, intent(out) :: VertLevels    ! Vertical levels; count of vert cell centers
integer, intent(out) :: VertexDeg     ! Max number of edges that touch any vertex
integer, intent(out) :: SoilLevels    ! Number of soil layers

if ( .not. module_initialized ) call static_init_model

Cells      = nCells
Vertices   = nVertices
Edges      = nEdges
VertLevels = nVertLevels
VertexDeg  = vertexDegree
SoilLevels = nSoilLevels

end subroutine get_grid_dims

!------------------------------------------------------------------

subroutine get_xland(Cells,LandOrNot)

! public routine for returning land mask
! Later, we may want to add more variables such as sfc_albedo.

integer,  intent(in)  :: Cells
real(r8), allocatable, intent(out) :: LandOrNot(:)

if ( .not. module_initialized ) call static_init_model()

allocate(LandOrNot(Cells))

LandOrNot  = xland

end subroutine get_xland

!------------------------------------------------------------------

subroutine get_surftype(Cells,surface_type)

! public routine for returning surface type (for rttov)
! As defined in atmos_profile_type in obs_def_rttov_mod.f90
! surface type (land=0, water=1, seaice = 2)

integer,  intent(in)  :: Cells
integer,  allocatable, intent(out) :: surface_type(:)

if ( .not. module_initialized ) call static_init_model()

allocate(surface_type(Cells))

! xland(:)   ! land-ocean mask (1=land including sea-ice ; 2=ocean)  
! seaice(:)  ! sea-ice flag (0=no seaice; =1 seaice) - for rttov

surface_type = 0                          ! land
where (seaice == 1.0_r8) surface_type  = 2     ! seaice
where ( xland == 2.0_r8) surface_type  = 1     ! ocean

end subroutine get_surftype

!------------------------------------------------------------------

subroutine get_cell_center_coords(Cells,Lats,Lons)

! public routine for returning cell center coordinates

integer,  intent(in)  :: Cells
real(r8), allocatable, intent(out) :: Lats(:), Lons(:)

if ( .not. module_initialized ) call static_init_model()

allocate(Lats(Cells), Lons(Cells))

Lats = latCell
Lons = lonCell

end subroutine get_cell_center_coords

!------------------------------------------------------------------

subroutine get_bdy_mask(Cells,Mask)

! public routine for returning mask for boundary cells

integer, intent(in)  :: Cells
integer, allocatable, intent(out) :: Mask(:)

if ( .not. module_initialized ) call static_init_model()

allocate(Mask(Cells))

Mask = bdyMaskCell

end subroutine get_bdy_mask


!==================================================================
! The (model-specific) private interfaces come last
!==================================================================


!------------------------------------------------------------------
!> convert time type into a character string with the
!> format of YYYY-MM-DD_hh:mm:ss

function time_to_string(t, interval)

 type(time_type), intent(in) :: t
 logical, intent(in), optional :: interval
character(len=TIMELEN) :: time_to_string

integer :: iyear, imonth, iday, ihour, imin, isec
integer :: ndays, nsecs
logical :: dointerval

   dointerval = .false.
if (present(interval)) dointerval = interval

! for interval output, output the number of days, then hours, mins, secs
! for date output, use the calendar routine to get the year/month/day hour:min:sec
if (dointerval) then
   call get_time(t, nsecs, ndays)
   if (ndays > 99) then
      write(string1, *) 'interval number of days is ', ndays
      call error_handler(E_ERR,'time_to_string', 'interval days cannot be > 99', &
                         source, text2=string1)
   endif
   ihour = nsecs / 3600
   nsecs = nsecs - (ihour * 3600)
   imin  = nsecs / 60
   nsecs = nsecs - (imin * 60)
   isec  = nsecs
   write(time_to_string, '(I2.2,3(A1,I2.2))') &
                        ndays, '_', ihour, ':', imin, ':', isec
else
   call get_date(t, iyear, imonth, iday, ihour, imin, isec)
   write(time_to_string, '(I4.4,5(A1,I2.2))') &
                        iyear, '-', imonth, '-', iday, '_', ihour, ':', imin, ':', isec
endif

end function time_to_string


!------------------------------------------------------------------

function string_to_time(s)

! parse a string to extract time.  the expected format of
! the string is YYYY-MM-DD_hh:mm:ss  (although the exact
! non-numeric separator chars are skipped and not validated.)

type(time_type) :: string_to_time
character(len=*), intent(in) :: s

integer :: iyear, imonth, iday, ihour, imin, isec

read( s ,'(i4,5(1x,i2))') iyear, imonth, iday, ihour, imin, isec
string_to_time = set_date(iyear, imonth, iday, ihour, imin, isec)

end function string_to_time


!------------------------------------------------------------------

function set_model_time_step()

! the static_init_model ensures that the model namelists are read.

type(time_type) :: set_model_time_step

if ( .not. module_initialized ) call static_init_model

! these are from the namelist
set_model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end function set_model_time_step

!------------------------------------------------------------------
!> Read the grid values in from the MPAS netcdf file.
!>

subroutine get_grid(ncid)
integer, intent(in) :: ncid

character(len=*), parameter :: routine = 'get_grid'


call nc_get_variable(ncid, 'lonCell',       lonCell,       routine)
call nc_get_variable(ncid, 'latCell',       latCell,       routine)

! MPAS locations are in radians - at this point DART needs degrees.
! watch out for tiny rounding errors and clamp to exactly +/- 90

latCell = latCell * rad2deg
lonCell = lonCell * rad2deg
where (latCell >  90.0_r8) latCell = 90.0_r8
where (latCell < -90.0_r8) latCell = -90.0_r8

call nc_get_variable(ncid, 'dcEdge',        dcEdge,        routine)
call nc_get_variable(ncid, 'zgrid',         zGridFace,     routine)
call nc_get_variable(ncid, 'cellsOnVertex', cellsOnVertex, routine)
call nc_get_variable(ncid, 'xland',         xland,         routine)
call nc_get_variable(ncid, 'seaice',        seaice,        routine)
call nc_get_variable(ncid, 'skintemp',      skintemp,      routine)

dxmax = maxval(dcEdge)  ! max grid resolution in meters

call nc_get_variable(ncid, 'edgeNormalVectors', edgeNormalVectors, routine)
call nc_get_variable(ncid, 'nEdgesOnCell',      nEdgesOnCell,      routine)
call nc_get_variable(ncid, 'edgesOnCell',       edgesOnCell,       routine)
call nc_get_variable(ncid, 'cellsOnEdge',       cellsOnEdge,       routine)
call nc_get_variable(ncid, 'verticesOnCell',    verticesOnCell,    routine)

if(data_on_edges) then
   call nc_get_variable(ncid, 'lonEdge', lonEdge, routine)
   call nc_get_variable(ncid, 'latEdge', latEdge, routine)

   latEdge = latEdge * rad2deg
   lonEdge = lonEdge * rad2deg
   where (latEdge >  90.0_r8) latEdge =  90.0_r8
   where (latEdge < -90.0_r8) latEdge = -90.0_r8

   call nc_get_variable(ncid, 'xEdge', xEdge, routine)
   call nc_get_variable(ncid, 'yEdge', yEdge, routine)
   call nc_get_variable(ncid, 'zEdge', zEdge, routine)
endif

call nc_get_variable(ncid, 'xVertex', xVertex, routine)
call nc_get_variable(ncid, 'yVertex', yVertex, routine)
call nc_get_variable(ncid, 'zVertex', zVertex, routine)

if (nc_variable_exists(ncid, 'maxLevelCell')) then
   allocate(maxLevelCell(nCells))
   call nc_get_variable(ncid, 'maxLevelCell', maxlevelCell, routine)
   all_levels_exist_everywhere = .false.
endif

end subroutine get_grid

!------------------------------------------------------------------

subroutine read_2d_from_nc_file(ncid, filename, varname, data)
 integer,          intent(in)  :: ncid
 character(len=*), intent(in)  :: filename
 character(len=*), intent(in)  :: varname
 real(r8),         intent(out) :: data(:,:)

!
! Read the values for all dimensions but the time dimension.
! Only read the last time (if more than 1 present)
!

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME) :: dimname
integer :: VarID, numdims, dimlen, i

call nc_check(nf90_inq_varid(ncid, varname, VarID), &
              'read_from_nc_file', &
              'inq_varid '//trim(varname)//' '//trim(filename))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
              'read_from_nc_file', &
              'inquire '//trim(varname)//' '//trim(filename))

do i=1, numdims
   write(string1,*)'inquire length for dimension ',i
   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen, name=dimname), &
                 'read_2d_from_nc_file', &
                  trim(string1)//' '//trim(filename))
   if (trim(dimname) == 'Time') then
      mystart(i)       = dimlen
      mycount(numdims) = 1
   else
      mystart(i)       = 1
      mycount(i)       = dimlen
   endif
enddo

call nc_check( nf90_get_var(ncid, VarID, data, &
               start=mystart(1:numdims), count=mycount(1:numdims)), &
              'read_2d_from_nc_file', &
              'get_var u '//trim(filename))

end subroutine read_2d_from_nc_file


!------------------------------------------------------------------

subroutine put_u(ncid, filename, u, vchar)
 character(len=*) :: vchar
 integer,  intent(in) :: ncid
 character(len=*), intent(in) :: filename
 real(r8), intent(in) :: u(:,:)       ! u(nVertLevels, nEdges)

! Put the newly updated 'u' field back into the netcdf file.

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount, numu
integer :: VarID, numdims, nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: ntimes, i

if ( .not. module_initialized ) call static_init_model

string2 = trim(filename)//' '//trim(vchar)

call nc_check(nf90_Inquire(ncid,nDimensions,nVariables,nAttributes,unlimitedDimID), &
              'put_u', 'inquire '//trim(filename))

call nc_check(nf90_inquire_dimension(ncid, unlimitedDimID, len=ntimes), &
              'put_u', 'inquire time dimension length '//trim(filename))

call nc_check(nf90_inq_varid(ncid, trim(vchar), VarID), &
              'put_u', 'inq_varid '//trim(string2))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
              'put_u', 'inquire '//trim(string2))

do i=1, numdims
   write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=numu(i)), &
                 'put_u', string1)
enddo

! for all but the time dimension, read all the values.
! for time read only the last one (if more than 1 present)
mystart = 1
mystart(numdims) = ntimes
mycount = numu
mycount(numdims) = 1

call nc_check(nf90_put_var(ncid, VarID, u, start=mystart, count=mycount), &
              'put_u', 'put_var '//trim(string2))

end subroutine put_u


!------------------------------------------------------------------
subroutine verify_state_variables(ncid, filename, ngood, qty_list, variable_bounds)

integer,                          intent(in)  :: ncid
character(len=*),                 intent(in)  :: filename
integer,                          intent(out) :: ngood
integer,                          intent(out) :: qty_list(:)
real(r8),                         intent(out) :: variable_bounds(:,:)

integer :: i, j, VarID
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: varname, dimname
character(len=NF90_MAX_NAME) :: dartstr
real(r8) :: lower_bound, upper_bound
character(len=10) :: bound
integer :: dimlen, numdims
logical :: u_already_in_list
logical :: failure

if ( .not. module_initialized ) call static_init_model

failure = .FALSE. ! this is for unsupported dimensions
u_already_in_list = .FALSE.

!HK @todo varible bounds
variable_bounds(:,:) = MISSING_R8   

ngood = 0
MyLoop : do i = 1, MAX_STATE_VARIABLES

   varname    = trim(mpas_state_variables(2*i -1))
   dartstr    = trim(mpas_state_variables(2*i   ))
   variable_table(i,1) = trim(varname)
   variable_table(i,2) = trim(dartstr) !HK @todo to_upper

   if (varname == 'u') has_edge_u = .true.
   if (varname == 'uReconstructZonal' .or. &
       varname == 'uReconstructMeridional') has_uvreconstruct = .true.

   if ( variable_table(i,1) == ' ' .and. variable_table(i,2) == ' ' ) exit MyLoop ! Found end of list.

   if ( variable_table(i,1) == ' ' .or. variable_table(i,2) == ' ' ) then
      string1 = 'mpas_vars_nml:model state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source)
   endif

   ! Make sure variable exists in model analysis variable list

   write(string1,'(''variable '',a,'' in '',a)') trim(varname), trim(filename)
   write(string2,'(''there is no '',a)') trim(string1)
   call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
                 'verify_state_variables', trim(string2))

   ! Make sure variable is defined by (Time,nCells) or (Time,nCells,vertical)
   ! unable to support Edges or Vertices at this time.

   call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
                 'verify_state_variables', 'inquire '//trim(string1))

   DimensionLoop : do j = 1,numdims

      write(string2,'(''inquire dimension'',i2,'' of '',a)') j,trim(string1)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(j), len=dimlen, name=dimname), &
                                          'verify_state_variables', trim(string2))
      select case ( trim(dimname) )
         case ('Time')
            ! supported - do nothing
         case ('nCells')
            ! supported - do nothing
         case ('nEdges')
            ! supported - do nothing
         case ('nVertLevels')
            ! supported - do nothing
         case ('nVertLevelsP1')
            ! supported - do nothing
         case ('nSoilLevels')
            ! supported - do nothing
         case default
            write(string2,'(''unsupported dimension '',a,'' in '',a)') trim(dimname),trim(string1)
            call error_handler(E_MSG,'verify_state_variables',string2,source)
            failure = .TRUE.
      end select

   enddo DimensionLoop

   if (failure) then
       string2 = 'unsupported dimension(s) are fatal'
       call error_handler(E_ERR,'verify_state_variables',string2,source)
   endif

   qty_list(i) = get_index_for_quantity(dartstr)
   if( qty_list(i) < 0 ) then
      write(string1,'(''there is no qty <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source)
   endif

   ! Keep track of whether U (edge normal winds) is part of the user-specified field list
   if ((variable_table(i, 1) == 'u') .or. (variable_table(i, 1) == 'U')) u_already_in_list = .true.

   ! check variable bounds
   do j = 1, MAX_STATE_VARIABLES
       if (trim(mpas_state_bounds(1, j)) == ' ') exit
         if (trim(mpas_state_bounds(1, j)) == trim(varname)) then

            bound = trim(mpas_state_bounds(2, j))
            if ( bound /= 'NULL' .and. bound /= '' ) then
               read(bound,'(d16.8)') lower_bound
             else
               lower_bound = MISSING_R8
             endif

             bound = trim(mpas_state_bounds(3, j))
             if ( bound /= 'NULL' .and. bound /= '' ) then
                  read(bound,'(d16.8)') upper_bound
             else
                  upper_bound = MISSING_R8
             endif
            variable_bounds(i,1) = lower_bound
            variable_bounds(i,2) = upper_bound
         endif
   enddo

   ngood = ngood + 1
enddo MyLoop

if (ngood == MAX_STATE_VARIABLES) then
   string1 = 'WARNING: There is a possibility you need to increase ''MAX_STATE_VARIABLES'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'verify_state_variables',string1,source,text2=string2)
endif

! if this flag is true and the user hasn't said U should be in the state,
! add it to the list.
if (add_u_to_state_list .and. .not. u_already_in_list) then
   ngood = ngood + 1
   variable_table(ngood,1) = "u"
   variable_table(ngood,2) = "QTY_EDGE_NORMAL_SPEED"
endif

! state bounds


end subroutine verify_state_variables


!------------------------------------------------------------------

!> @todo FIXME if you can call this *after* add_domain() has been
!> called, then we could use state structure calls for this.
!>
!> but right now, it's called first.  so pass in the namelist 
!> and boundary lists and base the decision on those.
!>
!> this routine can only be called after the namelist is read.  
!> also, it's called BY static_init_model() so it can't call
!> back into it.
!>

function is_edgedata_in_state_vector(state_list, bdy_list)
  character(len=*), intent(in) :: state_list(MAX_STATE_VARIABLES, NUM_STATE_TABLE_COLUMNS)
  character(len=*), intent(in) :: bdy_list(MAX_STATE_VARIABLES)
  logical :: is_edgedata_in_state_vector

integer :: i

StateLoop : do i = 1, MAX_STATE_VARIABLES

   if (state_list(i, 1) == ' ') exit StateLoop ! Found end of list.

   if (state_list(i, 1) == 'u') then
      is_edgedata_in_state_vector = .true.
      return
   endif

enddo StateLoop

! if U is not in the state, does it matter if U is
! in the boundary file?  yes, return true if so.
BdyLoop : do i = 1, MAX_STATE_VARIABLES

   if (bdy_list(i) == ' ') exit BdyLoop ! Found end of list.

   if (bdy_list(i) == 'lbc_u') then
      is_edgedata_in_state_vector = .true.
      return
   endif

enddo BdyLoop

is_edgedata_in_state_vector = .false.

end function is_edgedata_in_state_vector


!------------------------------------------------------------------

subroutine winds_present(zonal,meridional,both)

integer, intent(out) :: zonal, meridional
logical, intent(out) :: both

! if neither of uReconstructZonal or uReconstructMeridional are in the
!   state vector, set both to .false. and we're done.
! if both are there, return the ivar indices for each
! if only one is there, it's an error.

zonal      = get_varid_from_varname(anl_domid, 'uReconstructZonal')
meridional = get_varid_from_varname(anl_domid, 'uReconstructMeridional')

if (zonal > 0 .and. meridional > 0) then
  both = .true.
  return
else if (zonal < 0 .and. meridional < 0) then
  both = .false.
  return
endif

! only one present - error.
write(string1,*) 'both components for U/V reconstructed winds must be in state vector'
write(string2,*) 'cannot have only one of uReconstructMeridional and uReconstructZonal'
call error_handler(E_ERR,'winds_present',string1,source,text2=string2)

end subroutine winds_present

!------------------------------------------------------------------

subroutine compute_pressure_at_loc(state_handle, ens_size, location, ploc, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
real(r8),            intent(out) :: ploc(ens_size)
integer,             intent(out) :: istatus(ens_size)

! convert the vertical coordinate at the given location
! to a pressure.  if the vertical type is "undefined" then
! return 2001 mb.  this routine is currently only used to
! test for and reject obs which are above the upper limit
! (and for a 'VERTISUNDEF' obs it needs to return a large
! value with no error code).  it's also used if the 
! interpolation type is QTY_PRESSURE.  we don't explicitly
! prevent someone from creating a synthetic obs that has
! a VERTISUNDEF type, but it doesn't make sense in that case.

real(r8), dimension(3) :: llv, llv_new      ! lon/lat/vert
real(r8), dimension(3, ens_size) :: values
real(r8), dimension(ens_size)    :: tk
type(location_type), dimension(ens_size) :: new_location
integer :: ivars(3)
integer :: e ! loop index

! default is failure
ploc = MISSING_R8
istatus = 99

! base location
llv = get_location(location)

if(is_vertical(location, "PRESSURE")) then
   ploc(:) = llv(3)
   istatus(:) = 0

else if(is_vertical(location, "HEIGHT") .or. is_vertical(location, "LEVEL")) then
   new_location = location

   ! the quantity doesn't matter for this call but we have to pass in something.
   call convert_vert(state_handle, ens_size, new_location, QTY_TEMPERATURE, VERTISPRESSURE, istatus)

   do e = 1, ens_size
      if(istatus(e) == 0) then
         llv_new = get_location(new_location(e))
         ploc(e) = llv_new(3)
      endif
   enddo

else if(is_vertical(location, "SURFACE")) then


else if(is_vertical(location, "UNDEFINED")) then    ! not error, but no exact vert loc either
   ploc(:) = 200100.0_r8    ! this is an unrealistic pressure value to indicate no known pressure.
   istatus(:) = 0           ! see comment at top of this routine for why this is ok.

else
   call error_handler(E_ERR, 'compute_pressure:', 'internal error: unknown type of vertical', &
        source)
endif

end subroutine compute_pressure_at_loc

!------------------------------------------------------------------
!> convert the vertical type for one or more observation locations.

subroutine convert_vertical_obs(state_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, status)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:), loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: status(:)

integer :: i

do i=1, num
   call convert_vert(state_handle, 1, locs(i:i), loc_qtys(i), which_vert, status(i:i))
enddo

end subroutine convert_vertical_obs

!--------------------------------------------------------------------
!> convert the vertical type for one or more state locations.

subroutine convert_vertical_state(state_handle, num, locs, loc_qtys, loc_indx, &
                                  which_vert, istatus)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer(i8),         intent(in)    :: loc_indx(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus

integer :: i, status(1)

!> we are using code from get_state_meta_data to get
!> the indices for cell id and vert level.  this calls a
!> modified convert_vert() routine that doesn't
!> have to search for the cell centers. 

istatus = 0

do i=1, num
   call convert_vert_state(state_handle, 1, locs(i:i), loc_qtys(i), &
                                   loc_indx(i), which_vert, status)

   ! save the first error we see - but continue to convert the rest
   if (istatus == 0 .and. status(1) /= 0) istatus = status(1)
enddo


end subroutine convert_vertical_state


!------------------------------------------------------------------
!> code to convert an observation location's vertical coordinate type. !HK also state locations

subroutine convert_vert(state_handle, ens_size, location, obs_kind, ztypeout, istatus)

! This subroutine converts a given obs vertical coordinate to
! the vertical localization coordinate type requested through the
! model_mod namelist.
!
! Notes: (1) obs_kind is only necessary to check whether the ob
!            is an identity ob.
!        (2) This subroutine can convert both obs' and state points'
!            vertical coordinates. Remember that state points get
!            their DART location information from get_state_meta_data
!            which is called by filter_assim during the assimilation
!            process.
!        (3) state_handle is the relevant DART state vector for carrying out
!            computations necessary for the vertical coordinate
!            transformations. As the vertical coordinate is only used
!            in distance computations, this is actually the "expected"
!            vertical coordinate, so that computed distance is the
!            "expected" distance. Thus, under normal circumstances,
!            the state that is supplied to convert_vert should be the
!            ensemble mean. Nevertheless, the subroutine has the
!            functionality to operate on any DART state vector that
!            is supplied to it.

type(ensemble_type),    intent(in)    :: state_handle
integer,                intent(in)    :: ens_size
type(location_type),    intent(inout) :: location(ens_size)  ! because the verticals may differ
integer,                intent(in)    :: obs_kind
integer,                intent(in)    :: ztypeout
integer,                intent(out)   :: istatus(ens_size)

! zin and zout are the vert values coming in and going out.
! ztype{in,out} are the vert types as defined by the 3d sphere
! locations mod (location/threed_sphere/location_mod.f90)
real(r8), dimension(3, ens_size) :: llv_loc
real(r8), dimension(3,ens_size)  :: zk_mid, values, fract, fdata
integer,  dimension(3,ens_size)  :: k_low, k_up
real(r8), dimension(ens_size)    :: zin, zout
real(r8), dimension(ens_size)    :: tk, fullp, surfp
logical :: at_surf, do_norm
type(location_type), dimension(ens_size) :: surfloc

real(r8) :: weights(3)
integer  :: ztypein, i, e, n
integer  :: c(3), ivars(3) 

! assume failure.
istatus = 1

! initialization
k_low = 0.0_r8
k_up = 0.0_r8
weights = 0.0_r8

! first off, check if ob is identity ob.  if so get_state_meta_data() will
! return the location.
if (obs_kind < 0) then
   call get_state_meta_data(-int(obs_kind,i8),location(1)) ! will be the same across the ensemble
   location(:) = location(1)
endif

! if the existing coord is already in the requested vertical units
! or if the vert is 'undef' which means no specifically defined
! vertical coordinate, return now.
ztypein  = nint(query_location(location(1), 'which_vert'))
if ((ztypein == ztypeout) .or. (ztypein == VERTISUNDEF)) then
   istatus = 0
   return
endif

! we do need to convert the vertical.  start by
! extracting the location lon/lat/vert values.
do e = 1, ens_size
   llv_loc(:, e) = get_location(location(e))
enddo

! the routines below will use zin as the incoming vertical value
! and zout as the new outgoing one.  start out assuming failure
! (zout = missing) and wait to be pleasantly surprised when it works.
zin(:)     = llv_loc(3, :)
zout(:)    = missing_r8

! if the vertical is missing to start with, return it the same way
! with the requested type as out.
do e = 1, ens_size
   if (zin(e) == missing_r8) then
      location(e) = set_location(llv_loc(1, e),llv_loc(2, e),missing_r8,ztypeout)
   endif
enddo
! if the entire ensemble has missing vertical values we can return now.
! otherwise we need to continue to convert the members with good vertical values.
if (all(zin == missing_r8)) then
   istatus(:) = 0 
   return
endif

! Convert the incoming vertical type (ztypein) into the vertical
! localization coordinate given in the namelist (ztypeout).
! Various incoming vertical types (ztypein) are taken care of
! inside find_vert_level. So we only check ztypeout here.

! convert into:
select case (ztypeout)

   ! ------------------------------------------------------------
   ! outgoing vertical coordinate should be 'model level number'
   ! ------------------------------------------------------------
   case (VERTISLEVEL)

   ! Identify the three cell ids (c) in the triangle enclosing the obs and
   ! the vertical indices for the triangle at two adjacent levels (k_low and k_up)
   ! and the fraction (fract) for vertical interpolation.

   call find_triangle (location(1), n, c, weights, istatus(1))
   if(istatus(1) /= 0) then
      istatus(:) = istatus(1)
      location(:) = set_location(llv_loc(1, 1),llv_loc(2, 1),missing_r8,ztypeout)
      return
   endif
   call find_vert_indices (state_handle, ens_size, location(1), n, c, k_low, k_up, fract, istatus)
   if( all(istatus /= 0) ) then
      location(:) = set_location(llv_loc(1, 1),llv_loc(2, 1),missing_r8,ztypeout)
      return
   endif

   zk_mid = k_low + fract
   do e = 1, ens_size
      if(istatus(e) == 0) then
         zout(e) = sum(weights(:) * zk_mid(:, e))
      endif
   enddo

   ! ------------------------------------------------------------
   ! outgoing vertical coordinate should be 'pressure' in Pa
   ! ------------------------------------------------------------
   case (VERTISPRESSURE)

   ! Need to get base offsets for the potential temperature, density, and water
   ! vapor mixing fields in the state vector
   !ivars(1) = get_progvar_index_from_kind(QTY_POTENTIAL_TEMPERATURE)
   !ivars(2) = get_progvar_index_from_kind(QTY_DENSITY)
   !ivars(3) = get_progvar_index_from_kind(QTY_VAPOR_MIXING_RATIO)

   if (any(ivars(1:3) < 0)) then
      write(string1,*) 'Internal error, cannot find one or more of: theta, rho, qv'
      call error_handler(E_ERR, 'convert_vert',string1,source)
   endif

   ! Get theta, rho, qv at the interpolated location
   call compute_scalar_with_barycentric (state_handle, ens_size, location(1), 3, ivars, values, istatus)
   if( all(istatus /= 0) ) then
      location(:) = set_location(llv_loc(1, 1),llv_loc(2, 1),missing_r8,ztypeout)
      return
   endif

   ! Convert theta, rho, qv into pressure
   call compute_full_pressure(ens_size, values(1,:), values(2,:), values(3,:), zout, tk, istatus)

   ! ------------------------------------------------------------
   ! outgoing vertical coordinate should be 'height' in meters
   ! ------------------------------------------------------------
   case (VERTISHEIGHT)

   call find_triangle (location(1), n, c, weights, istatus(1))
   if(istatus(1) /= 0) then
      istatus(:) = istatus(1)
      location(:) = set_location(llv_loc(1, 1),llv_loc(2, 1),missing_r8,ztypeout)
      return
   endif
   call find_vert_indices (state_handle, ens_size, location(1), n, c, k_low, k_up, fract, istatus)
   if( all(istatus /= 0) ) then
      location(:) = set_location(llv_loc(1, 1),llv_loc(2, 1),missing_r8,ztypeout)
      return
   endif

   fdata = 0.0_r8
   do i = 1, n
      where (istatus == 0)
         fdata(i, :) = zGridFace(k_low(i, :),c(i))*(1.0_r8 - fract(i, :)) + &
                       zGridFace(k_up (i, :),c(i))*fract(i, :)
      end where
   enddo

   ! now have vertically interpolated values at cell centers.
   ! use horizontal weights to compute value at interp point.
   do e = 1, ens_size
      if(istatus(e) == 0) then
         zout(e) = sum(weights(:) * fdata(:,e))
      else
         zout(e) = missing_r8
      endif
   enddo

   ! ------------------------------------------------------------
   ! outgoing vertical coordinate should be 'scale height' (a ratio)
   ! ------------------------------------------------------------
   case (VERTISSCALEHEIGHT)

     ! Scale Height is defined as:  log(pressure) 
     ! if namelist item:  no_normalization_of_scale_heights = .true. 
     ! otherwise it is defined as: -log(pressure / surface_pressure)

     ! set logicals here so we can do the minimum amount of work.
     ! finding gridcells and computing pressure is expensive in this model.
     ! logic table is:
     !  surf T, norm T:  return 0.0 by definition
     !  surf T, norm F:  need surfp only
     !  surf F, norm F:  need fullp only
     !  surf F, norm T:  need both surfp and fullp

     at_surf = (ztypein == VERTISSURFACE)
     do_norm = .not. no_normalization_of_scale_heights

     ! if normalizing pressure and we're on the surface, by definition scale height 
     ! is log(1.0) so skip the rest of these computations.
     if (at_surf .and. do_norm) then
        zout = 0.0_r8  
        istatus(:) = 0
        return
     endif

     ! Base offsets for the potential temperature, density, and water
     ! vapor mixing fields in the state vector.
       !ivars(1) = get_progvar_index_from_kind(QTY_POTENTIAL_TEMPERATURE)
       !ivars(2) = get_progvar_index_from_kind(QTY_DENSITY)
       !ivars(3) = get_progvar_index_from_kind(QTY_VAPOR_MIXING_RATIO)

     if (at_surf .or. do_norm) then  ! we will need surface pressure

       ! Get theta, rho, qv at the surface corresponding to the interpolated location
       surfloc(1) = set_location(llv_loc(1, 1), llv_loc(2, 1), 1.0_r8, VERTISLEVEL)
       call compute_scalar_with_barycentric (state_handle, ens_size, surfloc(1), 3, ivars, values, istatus)
       if( all(istatus /= 0) ) then
           location(:) = set_location(llv_loc(1, 1),llv_loc(2, 1),missing_r8,ztypeout)
           return
       endif

       ! Convert surface theta, rho, qv into pressure
       call compute_full_pressure(ens_size, values(1,:), values(2,:), values(3,:), surfp, tk, istatus)

     endif

     if (.not. at_surf) then   ! we will need full pressure

       ! Get theta, rho, qv at the interpolated location
       call compute_scalar_with_barycentric (state_handle, ens_size, location(1), 3, ivars, values, istatus)

       ! Convert theta, rho, qv into pressure
       call compute_full_pressure(ens_size, values(1,:), values(2,:), values(3,:), fullp, tk, istatus)

     endif

     ! we have what we need now.  figure out the case and set zout to the right values.
     ! we've already taken care of the (at_surf .and. do_norm) case, so that simplifies
     ! the tests here.
     if (at_surf) then
        where (surfp /= MISSING_R8)
           zout = log(surfp)
       else where
         zout = MISSING_R8
       end where

     else if (.not. do_norm) then
        where (fullp /= MISSING_R8)
           zout = log(fullp)
        else where
           zout = MISSING_R8
        end where

     else  ! not at surface, and normalizing by surface pressure
        where (surfp /= MISSING_R8 .and. surfp > 0.0_r8 .and. fullp /= MISSING_R8)   
           zout = -log(fullp / surfp)
        else where
           zout = MISSING_R8
        end where
     endif

   ! -------------------------------------------------------
   ! outgoing vertical coordinate is unrecognized
   ! -------------------------------------------------------
   case default
      write(string1,*) 'Requested vertical coordinate not recognized: ', ztypeout
      call error_handler(E_ERR,'convert_vert', string1, &
                         source)

end select   ! outgoing vert type

! Returned location
do e = 1, ens_size
   if(istatus(e) == 0) then
      location(e) = set_location(llv_loc(1, e),llv_loc(2, e),zout(e),ztypeout)
   else
      location(e) = set_location(llv_loc(1, e),llv_loc(2, e),missing_r8,ztypeout)
   endif
enddo


end subroutine convert_vert

!------------------------------------------------------------------
!> code to convert an state location's vertical coordinate type.

subroutine convert_vert_state(state_handle, ens_size, location, quantity, state_indx, ztypeout, istatus)

! This subroutine converts a given state vertical coordinate to
! the vertical localization coordinate type requested through the
! model_mod namelist.
!
!        (3) state_handle is the relevant DART state vector for carrying out
!            computations necessary for the vertical coordinate
!            transformations. As the vertical coordinate is only used
!            in distance computations, this is actually the "expected"
!            vertical coordinate, so that computed distance is the
!            "expected" distance. Thus, under normal circumstances,
!            state_handle that is supplied to convert_vert should be the
!            ensemble mean. Nevertheless, the subroutine has the
!            functionality to operate on any DART state vector that
!            is supplied to it.

type(ensemble_type),    intent(in)    :: state_handle
integer,                intent(in)    :: ens_size
type(location_type),    intent(inout) :: location(ens_size)  ! because the verticals may differ
integer,                intent(in)    :: quantity
integer(i8),            intent(in)    :: state_indx
integer,                intent(in)    :: ztypeout
integer,                intent(out)   :: istatus(ens_size)

! zin and zout are the vert values coming in and going out.
! ztype{in,out} are the vert types as defined by the 3d sphere
! locations mod (location/threed_sphere/location_mod.f90)
real(r8), dimension(3, ens_size) :: llv_loc
real(r8), dimension(3,ens_size)  :: zk_mid, values, fract, fdata
integer,  dimension(3,ens_size)  :: k_low, k_up
real(r8), dimension(ens_size)    :: zin, zout
real(r8), dimension(ens_size)    :: tk, fullp, surfp
logical :: at_surf, do_norm, on_bound
type(location_type), dimension(ens_size) :: surfloc

real(r8) :: weights(3)
integer  :: ztypein, i, e, n, ndim
integer  :: c(3), ivars(3), vert_level
integer  :: cellid              !SYHA

! assume failure.
istatus = 1

! initialization
k_low = 0.0_r8
k_up = 0.0_r8
weights = 0.0_r8


! if the existing coord is already in the requested vertical units
! or if the vert is 'undef' which means no specifically defined
! vertical coordinate, return now.
ztypein  = nint(query_location(location(1), 'which_vert'))
if ((ztypein == ztypeout) .or. (ztypein == VERTISUNDEF)) then
   istatus = 0
   return
endif

!> assume that all locations have the same incoming lat/lon and level.
!> depending on the output vert type each member might have a different
!> vertical value.

! unpack the incoming location(s)
do e = 1, ens_size
   llv_loc(:, e) = get_location(location(e))
enddo

call find_mpas_indices(state_indx, cellid, vert_level, ndim)

! the routines below will use zin as the incoming vertical value
! and zout as the new outgoing one.  start out assuming failure
! (zout = missing) and wait to be pleasantly surprised when it works.
zin(:)     = vert_level
ztypein    = VERTISLEVEL
zout(:)    = missing_r8

! if the vertical is missing to start with, return it the same way
! with the requested type as out.
do e = 1, ens_size
   if (zin(e) == missing_r8) then 
      location(e) = set_location(llv_loc(1, e),llv_loc(2, e),missing_r8,ztypeout)
   endif
enddo
! if the entire ensemble has missing vertical values we can return now.
! otherwise we need to continue to convert the members with good vertical values.
! boundary cells will be updated by the assimilation.
! if all the vertical localization coord values are missing, 
! we don't call this routine again, and return.
if (all(zin == missing_r8)) then ! .or. on_bound) then
   istatus(:) = 0
   return
endif

! Convert the incoming vertical type (ztypein) into the vertical
! localization coordinate given in the namelist (ztypeout).
! Because this is only for state vector locations, we have already
! computed the vertical level, so all these conversions are from
! model level to something.

! convert into:
select case (ztypeout)

   ! ------------------------------------------------------------
   ! outgoing vertical coordinate should be 'model level number'
   ! ------------------------------------------------------------
   case (VERTISLEVEL)

   ! we have the vert_level and cellid - no need to call find_triangle or find_vert_indices

   zout(:) = vert_level

   ! ------------------------------------------------------------
   ! outgoing vertical coordinate should be 'pressure' in Pa
   ! ------------------------------------------------------------
   case (VERTISPRESSURE)

   !>@todo FIXME - this is the original code from the
   !> observation version which does horizontal interpolation.
   !> in this code we know we are on a state vector location
   !> so no interp is needed.  we should be able to make this
   !> more computationally efficient.

   ! Need to get base offsets for the potential temperature, density, and water
   ! vapor mixing fields in the state vector
   !ivars(1) = get_progvar_index_from_kind(QTY_POTENTIAL_TEMPERATURE)
   !ivars(2) = get_progvar_index_from_kind(QTY_DENSITY)
   !ivars(3) = get_progvar_index_from_kind(QTY_VAPOR_MIXING_RATIO)


   ! Get theta, rho, qv at the interpolated location - pass in cellid we have already located 
   ! to save the search time.
   call compute_scalar_with_barycentric (state_handle, ens_size, location(1), 3, ivars, values, istatus, cellid)
   if( all(istatus /= 0) ) then
      location(:) = set_location(llv_loc(1, 1),llv_loc(2, 1),missing_r8,ztypeout)
      return
   endif

   ! Convert theta, rho, qv into pressure
   call compute_full_pressure(ens_size, values(1,:), values(2,:), values(3,:), zout, tk, istatus)

   ! ------------------------------------------------------------
   ! outgoing vertical coordinate should be 'height' in meters
   ! ------------------------------------------------------------
   case (VERTISHEIGHT)

   ! surface obs should use the lower face of the first level.  the rest
   ! of the quantities should use the level centers.
   if ( ndim == 1 )  then
      zout(:) = zGridFace(1, cellid)
   else
      zout(:) = zGridCenter(vert_level, cellid)
      if ( quantity == QTY_VERTICAL_VELOCITY ) zout(:) = zGridFace(vert_level, cellid)
      if ( quantity == QTY_EDGE_NORMAL_SPEED ) zout(:) = zGridEdge(vert_level, cellid)
   endif

   ! ------------------------------------------------------------
   ! outgoing vertical coordinate should be 'scale height' (a ratio)
   ! ------------------------------------------------------------
   case (VERTISSCALEHEIGHT)

   !>@todo FIXME
   !> whatever we do for pressure, something similar here

     ! Scale Height is defined as:  log(pressure) 
     ! if namelist item:  no_normalization_of_scale_heights = .true. 
     ! otherwise it is defined as: -log(pressure / surface_pressure)

     ! set logicals here so we can do the minimum amount of work.
     ! finding gridcells and computing pressure is expensive in this model.
     ! logic table is:
     !  surf T, norm T:  return 0.0 by definition
     !  surf T, norm F:  need surfp only
     !  surf F, norm F:  need fullp only
     !  surf F, norm T:  need both surfp and fullp

     at_surf = (ztypein == VERTISSURFACE)
     do_norm = .not. no_normalization_of_scale_heights

     ! if normalizing pressure and we're on the surface, by definition scale height 
     ! is log(1.0) so skip the rest of these computations.
     if (at_surf .and. do_norm) then
        zout = 0.0_r8  
        istatus(:) = 0
        return
     endif

     ! Base offsets for the potential temperature, density, and water
     ! vapor mixing fields in the state vector.
       !ivars(1) = get_progvar_index_from_kind(QTY_POTENTIAL_TEMPERATURE)
       !ivars(2) = get_progvar_index_from_kind(QTY_DENSITY)
       !ivars(3) = get_progvar_index_from_kind(QTY_VAPOR_MIXING_RATIO)

     if (at_surf .or. do_norm) then  ! we will need surface pressure

       ! Get theta, rho, qv at the surface corresponding to the interpolated location
       surfloc(1) = set_location(llv_loc(1, 1), llv_loc(2, 1), 1.0_r8, VERTISLEVEL)
       call compute_scalar_with_barycentric (state_handle, ens_size, surfloc(1), 3, ivars, values, istatus)
       if( all(istatus /= 0) ) then
           location(:) = set_location(llv_loc(1, 1),llv_loc(2, 1),missing_r8,ztypeout)
           return
       endif

       ! Convert surface theta, rho, qv into pressure
       call compute_full_pressure(ens_size, values(1,:), values(2,:), values(3,:), surfp, tk, istatus)

     endif

     if (.not. at_surf) then   ! we will need full pressure

        ! Get theta, rho, qv at the interpolated location
        call compute_scalar_with_barycentric (state_handle, ens_size, location(1), 3, ivars, values, istatus)

        ! Convert theta, rho, qv into pressure
        call compute_full_pressure(ens_size, values(1,:), values(2,:), values(3,:), fullp, tk, istatus)

     endif

     ! we have what we need now.  figure out the case and set zout to the right values.
     ! we've already taken care of the (at_surf .and. do_norm) case, so that simplifies
     ! the tests here.
     if (at_surf) then
        where (surfp /= MISSING_R8)
           zout = log(surfp)
        else where
           zout = MISSING_R8
        end where

     else if (.not. do_norm) then
        where (fullp /= MISSING_R8)
           zout = log(fullp)
        else where
           zout = MISSING_R8
        end where

     else  ! not at surface, and normalizing by surface pressure
        where (surfp /= MISSING_R8 .and. surfp > 0.0_r8 .and. fullp /= MISSING_R8)   
           zout = -log(fullp / surfp)
        else where
           zout = MISSING_R8
        end where
     endif

   ! -------------------------------------------------------
   ! outgoing vertical coordinate is unrecognized
   ! -------------------------------------------------------
   case default
      write(string1,*) 'Requested vertical coordinate not recognized: ', ztypeout
      call error_handler(E_ERR,'convert_vert_state', string1, &
                         source)

end select   ! outgoing vert type

! Returned location
do e = 1, ens_size
   if(istatus(e) == 0) then
      location(e) = set_location(llv_loc(1, e),llv_loc(2, e),zout(e),ztypeout)
   else
      location(e) = set_location(llv_loc(1, e),llv_loc(2, e),missing_r8,ztypeout)
   endif
enddo


end subroutine convert_vert_state

!-------------------------------------------------------------------


!------------------------------------------------------------------
subroutine find_mpas_indices(index_in, cellid, vert_level, ndim)

integer(i8), intent(in)  :: index_in
integer,     intent(out) :: cellid
integer,     intent(out) :: vert_level
integer,     intent(out) :: ndim

integer  :: i, j, k ! Indices into variable (note k is not used in MPAS)
integer  :: nzp, iloc, vloc, nnf

call get_model_variable_indices(index_in, i, j, k, var_id=nnf)
ndim = get_num_dims(anl_domid, nnf)

if ( ndim == 2) then   ! 3d (2d netcdf ) variable(vcol, iloc)
   vert_level = i
   cellid = j
else  ! 2d (1d netcdf) variable: (cells)
   cellid = i
   vert_level = 1
endif

end subroutine find_mpas_indices

!==================================================================
! The following (private) interfaces are used for triangle interpolation
!==================================================================


!------------------------------------------------------------------
!> Finds position of a given height in an array of height grid points and returns
!> the index of the lower and upper bounds and the fractional offset.  ier
!> returns 0 unless there is an error. Could be replaced with a more efficient
!> search if there are many vertical levels.

subroutine find_height_bounds(height, nbounds, bounds, lower, upper, fract, ier)

real(r8), intent(in)  :: height
integer,  intent(in)  :: nbounds
real(r8), intent(in)  :: bounds(nbounds)
integer,  intent(out) :: lower(:), upper(:)
real(r8), intent(out) :: fract(:)
integer,  intent(out) :: ier(:)

! Assume that the spacing on altitudes is arbitrary and do the simple thing
! which is a linear search. Probably not worth any fancier searching unless
! models get to be huge.

integer :: i
real(r8) :: hlevel

! Initialization
ier = 0
fract = -1.0_r8
lower = -1
upper = -1

! FIXME:
! difficult to allow extrapolation because this routine only returns
! an upper and lower level, plus a fraction.  it would have to return
! a fraction < 0 or > 1 and then the calling code would have to test for
! and support extrapolation at some later point in the code.  so ignore
! extrapolation for now.
!
! do try to convert heights to grid levels and if within the limits
! return the top or bottom level as the value.  convert lower locations
! based on bounds(1)/(2) and upper on bounds(nbounds-1)/bounds(nbounds)

! Bounds checks:

! too low?
if(height < bounds(1)) then
   if(outside_grid_level_tolerance <= 0.0_r8) then
      ier = 998
   else
      ! let points outside the grid but "close enough" use the
      ! bottom level as the height.  hlevel should come back < 1
      hlevel = extrapolate_level(height, bounds(1), bounds(2))
      if (hlevel + outside_grid_level_tolerance >= 1.0_r8) then
         fract = 0.0_r8
         lower = 1
         upper = 2
      endif
   endif
endif

! too high?
if(height > bounds(nbounds)) then
   if(outside_grid_level_tolerance <= 0.0_r8) then
      ier = 9998
   else
      ! let points outside the grid but "close enough" use the
      ! top level as the height.  hlevel should come back > nbounds
      hlevel = extrapolate_level(height, bounds(nbounds-1), bounds(nbounds))
      if (hlevel - outside_grid_level_tolerance <= real(nbounds, r8)) then
         fract = 1.0_r8
         lower = nbounds-1
         upper = nbounds
      endif
   endif
endif


if (all(ier /= 0)) return

! Search two adjacent layers that enclose the given point
do i = 2, nbounds
   if(height <= bounds(i)) then
      lower = i - 1
      upper = i
      fract = (height - bounds(lower)) / (bounds(upper) - bounds(lower))
      return
   endif
end do

! should never get here because heights above and below the grid
! are tested for at the start of this routine.  if you get here
! there is a coding error.
ier = 3

end subroutine find_height_bounds


!------------------------------------------------------------------
!> given a location and 3 cell ids, return three sets of:
!> the two level numbers that enclose the given vertical value
!> plus the fraction between them for each of the 3 cell centers.
!> If the requested location uses a vertical coordinate that only
!> has one level, both level numbers are identical ... 1 and the
!> fractions are identical ... 0.0

subroutine find_vert_level(state_handle, ens_size, loc, nc, ids, oncenters, lower, upper, fract, ier)

! note that this code handles data at cell centers, at edges, but not
! data on faces.  so far we don't have any on faces.

! loc is an intrisic funtion
type(ensemble_type), intent(in)  :: state_handle
type(location_type), intent(in)  :: loc
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: nc, ids(:)
logical,             intent(in)  :: oncenters
integer,             intent(out) :: lower(nc, ens_size), upper(nc, ens_size)
real(r8),            intent(out) :: fract(nc, ens_size)
integer,             intent(out) :: ier(ens_size)

real(r8) :: lat, lon, vert, llv(3)
real(r8) :: vert_array(ens_size)
integer  :: track_ier(ens_size)
integer     :: verttype, i
integer     :: e

! Initialization
ier = 0 
lower(1:nc, :) = -1
upper(1:nc, :) = -1
fract(1:nc, :) = -1.0_r8

! the plan is to take in: whether this var is on cell centers or edges,
! and the location so we can extract the vert value and which vert flag.
! compute and return the lower and upper index level numbers that
! enclose this vert, along with the fract between them.  ier is set
! in case of error (e.g. outside the grid, on dry land, etc).

! kinds we have to handle:
!is_vertical:    VERTISUNDEF
!is_vertical:  VERTISSURFACE
!is_vertical:    VERTISLEVEL
!is_vertical: VERTISPRESSURE
!is_vertical:   VERTISHEIGHT

! unpack the location into local vars
! I think you can do this with the first ensemble member? 
! Because they are the same horizontally?
llv = get_location(loc)
lon  = llv(1)
lat  = llv(2)
vert = llv(3)
verttype = nint(query_location(loc))

! VERTISSURFACE describes variables on A surface (not necessarily THE surface)
! VERTISUNDEF describes no defined vertical location (e.g. vertically integrated vals)
if(verttype == VERTISSURFACE .or. verttype == VERTISUNDEF) then  ! same across the ensemble
   lower(1:nc, :) = 1
   upper(1:nc, :) = 1
   fract(1:nc, :) = 0.0_r8 
   ier = 0
   return
endif

! model level numbers (supports fractional levels)
if(verttype == VERTISLEVEL) then
   ! FIXME: if this is W, the top is nVertLevels+1
   if (vert > nVertLevels) then 
      ier(:) = 81
      return
   endif

   ! at the top we have to make the fraction 1.0; all other
   ! integral levels can be 0.0 from the lower level.
   if (vert == nVertLevels) then
      lower(1:nc, :) = nint(vert) - 1   ! round down
      upper(1:nc, :) = nint(vert)
      fract(1:nc, :) = 1.0_r8
      ier = 0
      return
   endif

   lower(1:nc, :) = aint(vert)   ! round down
   upper(1:nc, :) = lower+1
   fract(1:nc, :) = vert - lower
   ier = 0
   return

endif

if(verttype == VERTISPRESSURE ) then

   track_ier = 0
   vert_array = vert

   do i=1, nc
      call find_pressure_bounds(state_handle, ens_size, vert_array, ids(i), nVertLevels, &
            lower(i, :), upper(i, :), fract(i, :), ier)

      ! we are inside a loop over each corner. consolidate error codes
      ! so that we return an error for that ensemble member if any
      ! of the corners fails the pressure bounds test.
      where (ier /= 0 .and. track_ier == 0) track_ier = ier

   enddo

   ier = track_ier
   return
endif

if(verttype == VERTISHEIGHT) then
   ! For height, can do simple vertical search for interpolation for now
   ! Get the lower and upper bounds and fraction for each column
   track_ier = 0

   do i=1, nc
      if (oncenters) then
         call find_height_bounds(vert, nVertLevels, zGridCenter(:, ids(i)), &
                                 lower(i, :), upper(i, :), fract(i, :), ier)
      else

         if (.not. data_on_edges) then
            call error_handler(E_ERR, 'find_vert_level', &
                              'Internal error: oncenters false but data_on_edges false', &
                              source)
         endif
         call find_height_bounds(vert, nVertLevels, zGridEdge(:, ids(i)), &
                                 lower(i, :), upper(i, :), fract(i, :), ier)
      endif

      ! we are inside a loop over each corner. consolidate error codes
      ! so that we return an error for that ensemble member if any
      ! of the corners fails the pressure bounds test.
      where (ier /= 0 .and. track_ier == 0) track_ier = ier

   enddo

   ier = track_ier
   return

endif

end subroutine find_vert_level

!------------------------------------------------------------------

subroutine find_pressure_bounds(state_handle, ens_size, p, cellid, nbounds, &
   lower, upper, fract, ier)

! Finds vertical interpolation indices and fraction for a quantity with
! pressure vertical coordinate. Loops through the height levels and
! computes the corresponding pressure at the horizontal point.  nbounds is
! the number of vertical levels in the potential temperature, density,
! and water vapor grids.

type(ensemble_type), intent(in) :: state_handle
integer,     intent(in)  :: ens_size
real(r8),    intent(in)  :: p(ens_size)
integer,     intent(in)  :: cellid
integer,     intent(in)  :: nbounds ! number of vertical levels !HK @todo why not call it num_levels?
integer,     intent(out) :: lower(ens_size), upper(ens_size) 
real(r8),    intent(out) :: fract(ens_size) 
integer,     intent(out) :: ier(ens_size)

integer  :: i, ier2
real(r8) :: pr
real(r8) :: pressure(nbounds, ens_size)
integer  ::  e
integer  :: temp_ier(ens_size)
logical  :: found_level(ens_size)


! Initialize to bad values so unset returns will be caught.
fract = -1.0_r8
lower = -1
upper = -1

! this must start out success (0) and then as ensemble members fail
! we will record the first encountered error code.
ier = 0

! Find the lowest pressure
call get_interp_pressure(state_handle, ens_size, cellid, 1, pressure(1, :), temp_ier)

where(ier(:) == 0) ier(:) = temp_ier(:)

! Get the highest pressure level
call get_interp_pressure(state_handle, ens_size, cellid, nbounds, pressure(nbounds, :), temp_ier)

where(ier(:) == 0) ier(:) = temp_ier(:)

! Check for out of the column range
where(p(:) < pressure(nbounds, :)) ier(:) = 81   ! too high
!where(p(:) > pressure(      1, :)) ier(:) = 80  ! too low -> just take the lowest level

if(all(ier /= 0)) return

found_level(:) = .false.

! Check if the obs is located below the lowest model level in each member.
do e = 1, ens_size
   if(p(e) >= pressure(1, e)) then
      found_level(e) = .true.
            fract(e) = 0.0_r8
            lower(e) = 1
            upper(e) = 1
   endif
enddo !e = 1, ens_size
if(all(found_level)) return

! Loop through the rest of the column from the bottom up
do i = 2, nbounds
   ! we've already done this call for level == nbounds
   if (i /= nbounds) then
      call get_interp_pressure(state_handle, ens_size, cellid, i, pressure(i, :), temp_ier)
      where (ier(:) == 0) ier(:) = temp_ier(:)
   endif
   
   ! Check if pressure is not monotonically descreased with level.
   if(any(pressure(i, :) > pressure(i-1, :))) then
     where(pressure(i, :) > pressure(i-1, :)) ier(:) = PRESSURE_NOT_MONOTONIC
   endif

   ! each ensemble member could have a vertical between different levels,
   ! and more likely a different fract across a level.
   do e = 1, ens_size

      ! if we've already run into an error, or we've already found the
      ! level for this ensemble member, skip the rest of this loop.
      if (ier(e) /= 0)    cycle
      if (found_level(e)) cycle

      ! Is pressure between levels i-1 and i?  
      ! if so, set the lower and upper level numbers and fraction across.
      ! fraction is 0 at level (i-1) and 1 at level(i).
      if(p(e) >= pressure(i, e)) then
         found_level(e) = .true.
         lower(e) = i - 1
         upper(e) = i
         if (pressure(i, e) == pressure(i-1, e)) then
            fract(e) = 0.0_r8
         else if (log_p_vert_interp) then
            fract(e) = (log(p(e))          - log(pressure(i-1,e))) / &
                       (log(pressure(i,e)) - log(pressure(i-1,e)))
         else
            fract(e) = (p(e) - pressure(i-1,e)) / (pressure(i,e) - pressure(i-1,e))
         endif
      endif
   enddo
   if(all(found_level)) return
enddo

end subroutine find_pressure_bounds

!------------------------------------------------------------------

subroutine get_interp_pressure(state_handle, ens_size, cellid, lev, pressure, ier)

! Finds the value of pressure at a given point at model level lev

type(ensemble_type), intent(in) :: state_handle
integer,      intent(in)  :: ens_size
integer,      intent(in)  :: cellid
integer,      intent(in)  :: lev
real(r8),     intent(out) :: pressure(ens_size)
integer,      intent(out) :: ier(ens_size)

integer(i8) :: pt_idx, density_idx, qv_idx
integer :: dummyk
real(r8) :: pt(ens_size), density(ens_size), qv(ens_size), tk(ens_size)
integer :: e

! Get the values of potential temperature, density, and vapor

pt_idx      = get_dart_vector_index(lev, cellid, dummyk, anl_domid, get_varid_from_kind(anl_domid, QTY_POTENTIAL_TEMPERATURE))
density_idx = get_dart_vector_index(lev, cellid, dummyk, anl_domid, get_varid_from_kind(anl_domid, QTY_DENSITY))
qv_idx      = get_dart_vector_index(lev, cellid, dummyk, anl_domid, get_varid_from_kind(anl_domid, QTY_VAPOR_MIXING_RATIO))

pt      =  get_state(pt_idx, state_handle)
density =  get_state(density_idx, state_handle)
qv      =  get_state(qv_idx, state_handle)

! Initialization
ier = 0

! Error if any of the values are missing; probably will be all or nothing
!HK @todo why would the state be missing?
do e = 1, ens_size
   if(pt(e) == MISSING_R8 .or. density(e) == MISSING_R8 .or. qv(e) == MISSING_R8) then
      ier(e) = 2
      pressure(e) = MISSING_R8
   endif
enddo

! Convert theta, rho, qv into pressure
call compute_full_pressure(ens_size, pt(:), density(:), qv(:), pressure(:), tk(:), ier) !HK where is tk used?

end subroutine get_interp_pressure

!------------------------------------------------------------

subroutine get_3d_weights(p, v1, v2, v3, lat, lon, weights)

! Given a point p (x,y,z) inside a triangle, and the (x,y,z)
! coordinates of the triangle corner points (v1, v2, v3),
! find the weights for a barycentric interpolation.  this
! computation only needs two of the three coordinates, so figure
! out which quadrant of the sphere the triangle is in and pick
! the 2 axes which are the least planar:
!  (x,y) near the poles,
!  (y,z) near 0 and 180 longitudes near the equator,
!  (x,z) near 90 and 270 longitude near the equator.
! (lat/lon are the coords of p. we could compute them here
! but since in all cases we already have them, pass them
! down for efficiency)

real(r8), intent(in)  :: p(3)
real(r8), intent(in)  :: v1(3), v2(3), v3(3)
real(r8), intent(in)  :: lat, lon
real(r8), intent(out) :: weights(3)

real(r8) :: cxs(3), cys(3)

! above or below 45 in latitude, where -90 < lat < 90:
if (lat >= 45.0_r8 .or. lat <= -45.0_r8) then
   cxs(1) = v1(1)
   cxs(2) = v2(1)
   cxs(3) = v3(1)
   cys(1) = v1(2)
   cys(2) = v2(2)
   cys(3) = v3(2)
   call get_barycentric_weights(p(1), p(2), cxs, cys, weights)
   return
endif

! nearest 0 or 180 in longitude, where 0 < lon < 360:
if ( lon <= 45.0_r8 .or. lon >= 315.0_r8 .or. &
    (lon >= 135.0_r8 .and. lon <= 225.0_r8)) then
   cxs(1) = v1(2)
   cxs(2) = v2(2)
   cxs(3) = v3(2)
   cys(1) = v1(3)
   cys(2) = v2(3)
   cys(3) = v3(3)
   call get_barycentric_weights(p(2), p(3), cxs, cys, weights)
   return
endif

! last option, nearest 90 or 270 in lon:
cxs(1) = v1(1)
cxs(2) = v2(1)
cxs(3) = v3(1)
cys(1) = v1(3)
cys(2) = v2(3)
cys(3) = v3(3)
call get_barycentric_weights(p(1), p(3), cxs, cys, weights)

end subroutine get_3d_weights


!------------------------------------------------------------

subroutine get_barycentric_weights(x, y, cxs, cys, weights)

! Computes the barycentric weights for a 2d interpolation point
! (x,y) in a 2d triangle with the given (cxs,cys) corners.

real(r8), intent(in)  :: x, y, cxs(3), cys(3)
real(r8), intent(out) :: weights(3)

real(r8) :: denom

! Get denominator
denom = (cys(2) - cys(3)) * (cxs(1) - cxs(3)) + &
   (cxs(3) - cxs(2)) * (cys(1) - cys(3))

weights(1) = ((cys(2) - cys(3)) * (x - cxs(3)) + &
   (cxs(3) - cxs(2)) * (y - cys(3))) / denom

weights(2) = ((cys(3) - cys(1)) * (x - cxs(3)) + &
   (cxs(1) - cxs(3)) * (y - cys(3))) / denom

weights(3) = 1.0_r8 - weights(1) - weights(2)

if (any(abs(weights) < roundoff)) then
   where (abs(weights) < roundoff) weights = 0.0_r8
   where (abs(1.0_r8 - abs(weights)) < roundoff) weights = 1.0_r8
endif

end subroutine get_barycentric_weights


!------------------------------------------------------------
!> this routine computes 1 or more values at a single location,
!> for each ensemble member.  "n" is the number of different
!> quantities to compute, ival is an array of 'n' progval() indices
!> to indicate which quantities to compute.
!> dval(n, ens_size) are the output values, and ier(ens_size) are 
!> the success/error returns for each ensemble member.

subroutine compute_scalar_with_barycentric(state_handle, ens_size, loc, n, ival, dval, ier, this_cellid)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: loc !(ens_size)
integer,             intent(in)  :: n
integer,             intent(in)  :: ival(n)
real(r8),            intent(out) :: dval(n, ens_size)
integer,             intent(out) :: ier(ens_size)
integer, optional,   intent(in)  :: this_cellid

real(r8), dimension(3, ens_size) :: fract, fdata
real(r8), dimension(ens_size) :: lowval, uppval
real(r8)    :: weights(3)
integer     :: c(3), nvert, k, i, nc
integer(i8) :: low_state_indx(ens_size), upp_state_indx(ens_size)
integer     :: jdummy, kdummy
integer     :: lower(3, ens_size), upper(3,ens_size)
integer     :: e, e2, thislower, thisupper
logical     :: did_member(ens_size)

dval = MISSING_R8
ier = QTY_NOT_IN_STATE_VECTOR

! make sure we have all good field indices first
if (any(ival < 0)) return

call find_triangle (loc, nc, c, weights, ier(1), this_cellid) !HK c is cell centers
if(ier(1) /= 0) then
   ier(:) = ier(1)
   return
endif

! If the field is on a single level, lower and upper are both 1
call find_vert_indices (state_handle, ens_size, loc, nc, c, lower, upper, fract, ier)

if(fails(ier)) return

! for each field to compute at this location: !HK @todo why loop in here?  why not outside?
do k = 1, n

   if( get_num_dims(anl_domid, ival(k)) == 1 ) then  ! 1-D field
      do i = 1, nc
         ! go around triangle and interpolate in the vertical
         ! c(3) are the cell ids
         low_state_indx(1) =  get_dart_vector_index(c(i), jdummy, kdummy, anl_domid, ival(k))
         fdata(i,:) = get_state(low_state_indx(1), state_handle)
      enddo

   else ! 2-D field

      do i = 1, nc

         members: do e = 1, ens_size

            low_state_indx(e) =  get_dart_vector_index(lower(i,e), c(i), kdummy, anl_domid, ival(k))
            upp_state_indx(e) =  get_dart_vector_index(upper(i,e), c(i), kdummy, anl_domid, ival(k))
            
         enddo members

         call get_state_array(lowval(:), low_state_indx, state_handle)
         call get_state_array(uppval(:), upp_state_indx, state_handle)

         fdata(i, :) = lowval(:)*(1.0_r8 - fract(i, :)) + uppval(:)*fract(i,:)

      enddo  ! corners
 
   endif
   ! now have vertically interpolated values at cell centers.
   ! use weights to compute value at interp point.
   do e = 1, ens_size
    dval(k, e) = sum(weights(1:nc) * fdata(1:nc, e))
   enddo
enddo


end subroutine compute_scalar_with_barycentric

!------------------------------------------------------------
!> Interpolates metadata variables to an arbitrary location.
!> This routine can only be used with variables in module storage.

subroutine compute_surface_data_with_barycentric(var1d, loc, dval, ier, this_cellid)

type(location_type), intent(in)  :: loc
real(r8),            intent(in)  :: var1d(:) ! whole static variable. HK why pass this in?
real(r8),            intent(out) :: dval
integer,             intent(out) :: ier
integer, optional,   intent(in)  :: this_cellid

real(r8)    :: weights(3), fdata(3)
integer     :: c(3), i, nc

! assume failure
dval = MISSING_R8

call find_triangle (loc, nc, c, weights, ier, this_cellid)
if(ier /= 0) return

do i = 1, nc
   fdata(i) = var1d(c(i))    ! selected cell ID number
enddo

! use weights to compute value at interp point.
dval = sum(weights(1:nc) * fdata(1:nc))

end subroutine compute_surface_data_with_barycentric

!------------------------------------------------------------

subroutine find_triangle(loc, nc, c, weights, ier, this_cellid)

type(location_type), intent(in)  :: loc 
integer,             intent(out) :: nc
integer,             intent(out) :: c(:) ! single value - cell id
real(r8),            intent(out) :: weights(:)
integer,             intent(out) :: ier
integer, optional,   intent(in)  :: this_cellid

! compute the values at the correct vertical level for each
! of the 3 cell centers defining a triangle that encloses the
! the interpolation point, then interpolate once in the horizontal
! using barycentric weights to get the value at the interpolation point.

integer, parameter :: listsize = 30
integer  :: nedges, i, neighborcells(maxEdges), edgeid
real(r8) :: xdata(listsize), ydata(listsize), zdata(listsize)
real(r8) :: t1(3), t2(3), t3(3), r(3)
integer  :: cellid, verts(listsize), closest_vert
real(r8) :: lat, lon, vert, llv(3)
integer  :: verttype, vindex, v, vp1
logical  :: inside, foundit

! initialization
      c = MISSING_I
weights = 0.0_r8
    ier = 0
     nc = 1

! unpack the location into local vars
llv = get_location(loc)
lon  = llv(1)
lat  = llv(2)
vert = llv(3)
verttype = nint(query_location(loc))

! if we already know the closest cell center, pass it in instead
! of searching for it again.
if (present(this_cellid)) then
   cellid = this_cellid
else
   cellid = find_closest_cell_center(lat, lon)
endif



if (cellid < 1) then
   ier = 11
   return
endif

c(1) = cellid

! closest vertex to given point.
closest_vert = closest_vertex_ll(cellid, lat, lon)

! collect the neighboring cell ids and vertex numbers
! this 2-step process avoids us having to read in the
! cellsOnCells() array which i think we only need here.
! if it comes up in more places, we can give up the space
! and read it in and then this is a direct lookup.
! also note which index is the closest vert and later on
! we can start the triangle search there.
vindex = 1
nedges = nEdgesOnCell(cellid)
do i=1, nedges
   edgeid = edgesOnCell(i, cellid)
   if (.not. global_grid .and. &
      (cellsOnEdge(1, edgeid) <= 0 .or. cellsOnEdge(2, edgeid) <= 0)) then
      ier = TRIANGLE_CELL_CENTER_NOT_FOUND 
      return
   endif
   if (cellsOnEdge(1, edgeid) /= cellid) then
      neighborcells(i) = cellsOnEdge(1, edgeid)
   else
      neighborcells(i) = cellsOnEdge(2, edgeid)
   endif
   verts(i) = verticesOnCell(i, cellid)
   if (verts(i) == closest_vert) vindex = i
   call latlon_to_xyz(latCell(neighborcells(i)), lonCell(neighborcells(i)), &
      xdata(i), ydata(i), zdata(i))
enddo


! get the cartesian coordinates in the cell plane for the closest center
call latlon_to_xyz(latCell(cellid), lonCell(cellid), t1(1), t1(2), t1(3))

! and the observation point
call latlon_to_xyz(lat, lon, r(1), r(2), r(3))

if (all(abs(t1-r) < roundoff)) then   ! Located at a grid point (counting roundoff errors)

   ! need vert index for the vertical level
   nc = 1

   ! This is a grid point - horiz interpolation NOT needed
   weights(1) = 1.0_r8

else                       ! an arbitrary point

   ! find the cell-center-tri that encloses the obs point
   ! figure out which way vertices go around cell?
   foundit = .false.
   findtri: do i=vindex, vindex+nedges
      v = mod(i-1, nedges) + 1
      vp1 = mod(i, nedges) + 1
      t2(1) = xdata(v)
      t2(2) = ydata(v)
      t2(3) = zdata(v)
      t3(1) = xdata(vp1)
      t3(2) = ydata(vp1)
      t3(3) = zdata(vp1)
      call inside_triangle(t1, t2, t3, r, lat, lon, inside, weights)
      if (inside) then
         ! weights are the barycentric weights for the point r
         ! in the triangle formed by t1, t2, t3.
         ! v and vp1 are vert indices which are same indices
         ! for cell centers
         nc = 3
         c(2) = neighborcells(v)
         c(3) = neighborcells(vp1)
         foundit = .true.
         exit findtri
      endif
   enddo findtri
   if (.not. foundit) then
      ier = TRIANGLE_CELL_CENTER_NOT_FOUND
      return
   endif

endif     ! horizontal index search is done now.

end subroutine find_triangle

!------------------------------------------------------------
!> This routine adds some error checking because find_vert_level
!> has some early return statements that prevent having the debugging
!> information in it.

subroutine find_vert_indices (state_handle, ens_size, loc, nc, c, lower, upper, fract, ier)

type(ensemble_type), intent(in)  :: state_handle
type(location_type), intent(in)  :: loc
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: nc
integer,             intent(in)  :: c(nc)
integer,             intent(out) :: lower(nc, ens_size), upper(nc, ens_size)
real(r8),            intent(out) :: fract(nc, ens_size) ! ens_size
integer,             intent(out) :: ier(ens_size)

! initialization
lower = MISSING_I
upper = MISSING_I
fract = 0.0_r8

! need vert index for the vertical level
call find_vert_level(state_handle, ens_size, loc, nc, c, .true., lower, upper, fract, ier)

end subroutine find_vert_indices

!------------------------------------------------------------

subroutine compute_u_with_rbf(state_handle, ens_size, loc, zonal, uval, ier)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: loc
logical,             intent(in)  :: zonal
real(r8),            intent(out) :: uval(ens_size)
integer,             intent(out) :: ier(ens_size) 

! max edges we currently use is ~50, but overestimate for now
integer, parameter :: listsize = 200
logical, parameter :: on_a_sphere = .true.
integer     :: nedges, edgelist(listsize), i, j, nvert
real(r8)    :: xdata(listsize), ydata(listsize), zdata(listsize)
real(r8)    :: edgenormals(3, listsize)
real(r8)    :: veldata(listsize, ens_size)
real(r8)    :: xreconstruct, yreconstruct, zreconstruct
real(r8)    :: ureconstructx, ureconstructy, ureconstructz
real(r8)    :: ureconstructzonal, ureconstructmeridional
real(r8)    :: datatangentplane(3,2)
real(r8)    :: coeffs_reconstruct(3,listsize)
integer(i8) :: upindx(ens_size), lowindx(ens_size)
integer     :: cellid, vertexid
real(r8)    :: lat, lon, vert, llv(3), fract(listsize, ens_size), lowval(ens_size), uppval(ens_size)
integer     :: verttype, lower(listsize, ens_size), upper(listsize, ens_size), ncells, celllist(listsize)
integer     :: var_id, dummy

integer :: e ! loop index

ier  = 0
uval = MISSING_R8

var_id = get_varid_from_varname(anl_domid, 'u')
if (var_id < 0 .or. .not. data_on_edges) then
   ! cannot compute u if it is not in the state vector, or if we
   ! have not read in the edge data (which shouldn't happen if
   ! u is in the state vector.
   ier = RBF_U_COMPUTATION_ERROR
   return
endif

! unpack the location into local vars
llv = get_location(loc) ! I believe this is the same across the ensemble
lon = llv(1)
lat = llv(2)
vert = llv(3)
verttype = nint(query_location(loc))

call find_surrounding_edges(lat, lon, nedges, edgelist, cellid, vertexid)
if (nedges <= 0) then ! we are on a boundary, no interpolation
    ier = 18
    return
 endif

if (verttype == VERTISPRESSURE) then
   ! get all cells which share any edges with the 3 cells which share
   ! the closest vertex.
   call make_cell_list(vertexid, 3, ncells, celllist)

   call find_vert_level(state_handle, ens_size, loc, ncells, celllist, .true., &
                         lower, upper, fract, ier)

   if (fails(ier)) return

   ! now have pressure at all cell centers - need to interp to get pressure
   ! at edge centers.
   call move_pressure_to_edges(ncells, celllist, lower, upper, fract, &
                                nedges, edgelist, ier)
   if (fails(ier)) return

else
   ! need vert index for the vertical level
   call find_vert_level(state_handle, ens_size, loc, nedges, edgelist, .false., &
                         lower, upper, fract, ier)
    if (fails(ier)) return
endif

! the rbf code needs (their names == our names):
! nData == nedges
! xyz data == xyzEdge
! normalDirectionData == edgeNormalVectors
! velocitydata = U field

do i = 1, nedges
   !>@todo FIXME:
   ! ryan has size(xEdge, 1) but this is a 1d array, so a dimension shouldn't be needed
   ! was possibly required for ifort v17 on cheyenne
   if (edgelist(i) > size(xEdge,1)) then
       write(string1, *) 'edgelist has index larger than edge count', &
                           i, edgelist(i), size(xEdge,1)
      call error_handler(E_ERR, 'compute_u_with_rbf', 'internal error', &
                          source, text2=string1)
   endif

   xdata(i) = xEdge(edgelist(i))
   ydata(i) = yEdge(edgelist(i))
   zdata(i) = zEdge(edgelist(i))

   do j=1, 3
      edgenormals(j, i) = edgeNormalVectors(j, edgelist(i))
   enddo

  ! Lower/upper could be different levels in pressure
   do e = 1, ens_size
      lowindx = get_dart_vector_index(lower(i, e), edgelist(i), dummy,anl_domid, var_id) !HK @todo check x,y,z
      upindx = get_dart_vector_index(upper(i, e), edgelist(i), dummy,anl_domid, var_id) !HK @todo check x,y,z
   enddo

   call get_state_array(lowval, lowindx, state_handle)
   call get_state_array(uppval, upindx, state_handle)
   veldata(i, :) = lowval*(1.0_r8 - fract(i, :)) + uppval*fract(i, :)


enddo


! this intersection routine assumes the points are on the surface
! of a sphere.  they do NOT try to make a plane between the three
! points.  the difference is small, granted.
call latlon_to_xyz(lat, lon, xreconstruct,yreconstruct,zreconstruct)

! call a simple subroutine to define vectors in the tangent plane
call get_geometry(nedges, xdata, ydata, zdata, &
               xreconstruct, yreconstruct, zreconstruct, edgenormals, &
               on_a_sphere, datatangentplane)

! calculate coeffs_reconstruct
call get_reconstruct_init(nedges, xdata, ydata, zdata, &
              xreconstruct, yreconstruct, zreconstruct, edgenormals, &
              datatangentplane, coeffs_reconstruct)

! do the reconstruction
do e = 1, ens_size
   call get_reconstruct(nedges, lat*deg2rad, lon*deg2rad, &
                coeffs_reconstruct, on_a_sphere, veldata(:, e), &
                ureconstructx, ureconstructy, ureconstructz, &
                ureconstructzonal, ureconstructmeridional)

   if (zonal) then
       uval(e) = ureconstructzonal
   else
       uval(e) = ureconstructmeridional
   endif

enddo

end subroutine compute_u_with_rbf

! !------------------------------------------------------------

subroutine find_surrounding_edges(lat, lon, nedges, edge_list, cellid, vertexid)
 real(r8), intent(in)  :: lat, lon
 integer,  intent(out) :: nedges, edge_list(:)
 integer,  intent(out) :: cellid, vertexid

! ! given an arbitrary lat/lon location, find the edges of the
! ! cells that share the nearest vertex.  return the ids for the
! ! closest cell center and closest vertex to save having to redo
! ! the search.

! integer :: vertex_list(30), nedges1, edge_list1(300), c, i
! integer :: ncells, cell_list(150), nedgec, edgeid, nverts

! ! find the cell id that has a center point closest
! ! to the given point.
! cellid = find_closest_cell_center(lat, lon)
! if (cellid < 1) then
!    nedges = 0
!    edge_list(:) = -1
!    return
! endif

! ! inside this cell, find the vertex id that the point
! ! is closest to.  this is a return from this subroutine.
! vertexid = closest_vertex_ll(cellid, lat, lon)
! if (vertexid <= 0) then
!    ! call error handler?  unexpected
!    nedges = 0
!    edge_list(:) = -1
!    return
! endif

! nedges = 0
! edge_list = 0

! select case (use_rbf_option)
!   case (0)
!       ! use edges on the closest cell only.
!       nedges = nEdgesOnCell(cellid)
!       do i=1, nedges
!          edge_list(i) = edgesOnCell(i, cellid)
!       enddo
!   case (1)
!       ! fill in the number of unique edges and fills the
!       ! edge list with the edge ids.  the code that detects
!       ! boundary edges for the ocean or regional atmosphere
!       ! is incorporated here.  nedges can come back 0 in that case.
!       vertex_list(1) = vertexid
!       call make_edge_list_from_verts(1, vertex_list, nedges, edge_list)

!    case (2)
!       ! for all vertices of the enclosing cell, add the
!       ! edges which share this vertex.
!       nverts = nEdgesOnCell(cellid)
!       do i=1, nverts
!          vertex_list(i) = verticesOnCell(i, cellid)
!       enddo

!       ! fill in the number of unique edges and fills the
!       ! edge list with the edge ids.  the code that detects
!       ! boundary edges for the ocean or regional atmosphere
!       ! is incorporated here.  nedges can come back 0 in that case.
!       call make_edge_list_from_verts(nverts, vertex_list, nedges, edge_list)

!    case (3)
!       call make_cell_list(vertexid, 1, ncells, cell_list)

!       ! for all cells:
!       do c=1, ncells
!          ! for all vertices of the enclosing cell, add the
!          ! edges which share this vertex.
!          nverts = nEdgesOnCell(cell_list(c))
!          do i=1, nverts
!             vertex_list(i) = verticesOnCell(i, cell_list(c))
!          enddo

!          ! fill in the number of unique edges and fills the
!          ! edge list with the edge ids.  the code that detects
!          ! boundary edges for the ocean or regional atmosphere
!          ! is incorporated here.  nedges can come back 0 in that case.
!          call make_edge_list_from_verts(nverts, vertex_list, nedges1, edge_list1)
!          call merge_edge_lists(nedges, edge_list, nedges1, edge_list1)
!       enddo

!    case (4)
!       ! this one gives the same edge list as case (2).  see what's faster.
!       nedgec = nEdgesOnCell(cellid)
!       do i=1, nedgec
!          edgeid = edgesOnCell(i, cellid)
!          if (cellsOnEdge(1, edgeid) /= cellid) then
!             cell_list(i) = cellsOnEdge(1, edgeid)
!          else
!             cell_list(i) = cellsOnEdge(2, edgeid)
!          endif
!       enddo

!       call make_edge_list_from_cells(nedgec, cell_list, nedges, edge_list)

!    case default
!       call error_handler(E_ERR, 'find_surrounding_edges', 'bad use_rbf_option value', &
!                          source)
! end select

! ! Check if any of edges are located in the boundary zone.
! ! (We will skip the obs if any edges are located there.)
! if (on_boundary_edgelist(edge_list)) then
!    nedges = -1
!    edge_list(:) = -1
! endif

end subroutine find_surrounding_edges

!------------------------------------------------------------

subroutine init_closest_center()

! use nCells, latCell, lonCell to initialize a GC structure
! to be used later in find_closest_cell_center().

! set up a GC in the locations mod

integer :: i

allocate(cell_locs(nCells))

do i=1, nCells
   cell_locs(i) = xyz_set_location(lonCell(i), latCell(i), 0.0_r8, radius)
enddo

! get 2nd arg from max dcEdge.
call xyz_get_close_init(cc_gc, dxmax, nCells, cell_locs)

if (.not. global_grid) &
   call xyz_use_great_circle_dist(radius)

end subroutine init_closest_center

!------------------------------------------------------------
!> Determine the cell index for the closest center to the given point
!> 2D calculation only.   If the closest cell is on the boundary for
!> the regional case, the location is considered outside the region.
!> for the global case this can't happen.

function find_closest_cell_center(lat, lon)

real(r8), intent(in)  :: lat, lon
integer               :: find_closest_cell_center

type(xyz_location_type) :: pointloc
integer :: closest_cell, rc

! do this exactly once.
if (.not. search_initialized) then
   call init_closest_center()
   search_initialized = .true.
endif

pointloc = xyz_set_location(lon, lat, 0.0_r8, radius)

call xyz_find_nearest(cc_gc, pointloc, cell_locs, closest_cell, rc)

!> updated xyz_find_nearest to return -1 if outside the volume.
!> this is a code change so allow -1 returns here.  make sure
!> calling code is ready to handle a -1 cellid.  should only
!> happen in the regional case.

if (rc /= 0 .or. closest_cell < 0) then
   find_closest_cell_center = -1
   return
endif


! do allow boundary cells to be returned.
find_closest_cell_center = closest_cell

end function find_closest_cell_center

!------------------------------------------------------------

subroutine finalize_closest_center()

! get rid of storage associated with GC for cell centers if
! they were used.

if (search_initialized) call xyz_get_close_destroy(cc_gc)

end subroutine finalize_closest_center

!------------------------------------------------------------
!> we determine whether this is a local or global grid by
!> whether the 'bdyMaskXXX' arrays have values > 0 in them.
!> unfortunately both global and regional grids have these
!> arrays so we have to read them in and check the values
!> before knowing whether this is a global or regional case.
!> if regional, keep these arrays in memory. 
!> if global, discard them.

subroutine set_global_grid(ncid)
integer, intent(in) :: ncid

character(len=*), parameter :: routine = 'set_global_grid'

global_grid = .true.

if (nc_variable_exists(ncid, 'bdyMaskCell')) then
   allocate(bdyMaskCell(nCells))
   call nc_get_variable(ncid, 'bdyMaskCell', bdyMaskCell, routine)
   if(maxval(bdyMaskCell) > 0) then
      global_grid = .false.
   else
     deallocate(bdyMaskCell)
   endif
endif

if (nc_variable_exists(ncid, 'bdyMaskEdge')) then
   allocate(bdyMaskEdge(nEdges))
   call nc_get_variable(ncid, 'bdyMaskEdge', bdyMaskEdge, routine)
   if(maxval(bdyMaskEdge) > 0) then
      global_grid = .false.
   else
      deallocate(bdyMaskEdge)
   endif
endif

if (global_grid) then
   string1 = 'MPAS running in global mode'
else
   string1 = 'MPAS running in regional (limited-area) mode'
endif
call error_handler(E_MSG,'set_global_grid',string1,source)

end subroutine set_global_grid

!------------------------------------------------------------
!> accessor function for global_grid module variable

function is_global_grid()
logical :: is_global_grid

is_global_grid = global_grid 

end function is_global_grid

!------------------------------------------------------------
!> find the closest cell center.  if global grid, return it.
!> if regional grid, continue to be sure that all three corners
!> of the enclosing triangle are also inside the main part of
!> the regional grid (none are in the boundary layers)
!> return -1 if not ok.

function cell_ok_to_interpolate(location)

type(location_type), intent(in)  :: location
integer :: cell_ok_to_interpolate

integer :: cellid, nc, c(3), istatus, i
real(r8) :: llv(3), lat, lon
real(r8) :: weights(3)

llv = get_location(location)
lat = llv(2)
lon = llv(1)

! if outside a regional grid or the grid is
! global, return here.
cellid = find_closest_cell_center(lat, lon)
if (cellid < 1 .or. global_grid) then
   cell_ok_to_interpolate = cellid
   return
endif

! quick check to see if the current cell is
! a regional boundary cell
if (on_boundary_cell(cellid)) then
   cell_ok_to_interpolate = -1
   return
endif

! for regional, continue on and find the other 2 cell centers
! that create a triangle enclosing this location.  verify the
! other vertices are also fully inside the grid and not in the
! boundary layers.

call find_triangle(location, nc, c, weights, istatus, cellid)
if (istatus /= 0) then
   cell_ok_to_interpolate = -1
   return
endif

do i=1, nc
   if (on_boundary_cell(c(i))) then
      cell_ok_to_interpolate = -1
      return
   endif
enddo

cell_ok_to_interpolate = cellid

end function cell_ok_to_interpolate

!------------------------------------------------------------
!> Determine if this cell is on the boundary, and return true if so.   
!> if the global flag is set, skip all code and return false immediately.
!> Unlike the previous version of on_boundary, we do not return true
!> if any surrounding edges belong to the boundary zone. We will take care
!> of those edges either in uv_cell_to_edges or in find_surrounding_edges
!> individually.

function on_boundary_cell(cellid)

integer,  intent(in)  :: cellid
logical               :: on_boundary_cell

on_boundary_cell = .false.

if (global_grid .or. .not. allocated(bdyMaskCell)) return

if (bdyMaskCell(cellid) > 0) on_boundary_cell = .true.

end function on_boundary_cell

!------------------------------------------------------------
!> Determine if this edge is on the boundary

function on_boundary_edge(edgeid)

integer,  intent(in)  :: edgeid
logical               :: on_boundary_edge

on_boundary_edge = .false.

if (global_grid .or. .not. allocated(bdyMaskEdge)) return

if (bdyMaskEdge(edgeid) > 0) on_boundary_edge = .true.

end function on_boundary_edge

!------------------------------------------------------------
!> Determine if this edge is the outermost edge (bdyMaskEdge = 7)

function on_outermost_edge(edgeid)

integer,  intent(in)  :: edgeid
logical               :: on_outermost_edge

on_outermost_edge = .false.

if (global_grid .or. .not. allocated(bdyMaskEdge)) return

if (bdyMaskEdge(edgeid) > 6) on_outermost_edge = .true.

end function on_outermost_edge

!------------------------------------------------------------
!> Determine if this cell is the outermost cell (bdyMaskCell = 7)

function on_outermost_cell(cellid)

integer,  intent(in)  :: cellid
logical               :: on_outermost_cell

on_outermost_cell = .false.

if (global_grid .or. .not. allocated(bdyMaskCell)) return

if (bdyMaskCell(cellid) > 6) on_outermost_cell = .true.

end function on_outermost_cell

!------------------------------------------------------------
!> Determine if any of these edges are on the boundary

function on_boundary_edgelist(edgeids)

integer,  intent(in)  :: edgeids(:)
logical               :: on_boundary_edgelist

on_boundary_edgelist = .false.

if (global_grid .or. .not. allocated(bdyMaskEdge)) return

if (any(bdyMaskEdge(edgeids) > 0)) on_boundary_edgelist = .true.

end function on_boundary_edgelist

!------------------------------------------------------------

!> based on boundary cell mask value, return the weight
!> given to the analysis values.  interior values are 1.0
!> (full analysis); exterior values are 0.0 (full boundary)

function get_analysis_weight(cellid, is_it_cell)
integer, intent(in) :: cellid
real(r8)            :: get_analysis_weight

logical, optional, intent(in) :: is_it_cell

! Adopted from MPAS/src/core_init_atmosphere/mpas_init_atm_cases.F
integer, parameter  :: nBdyLayers = 7   ! The number of relaxation layers plus the number of specified layers
integer, parameter  :: nSpecLayers = 2  ! The number of specified layers

! local variables
integer :: imask
logical :: if_cell

if_cell = .true.    ! default is cells, not edges.
imask   = bdyMaskCell(cellid)

if (present(is_it_cell)) if_cell = is_it_cell

if (global_grid) then
    get_analysis_weight = 1.0_r8
    return
endif

if(.not. if_cell) then
   imask = bdyMaskEdge(cellid)
endif

if (imask > (nBdyLayers - nSpecLayers)) then  ! Specified Zone (6,7)
    get_analysis_weight = 0.0_r8
else        ! Relexation Zone (1-5)
    get_analysis_weight = 1.0_r8 - real(imask,kind=r8)/real(nBdyLayers - nSpecLayers,kind=r8)
end if

end function get_analysis_weight

!------------------------------------------------------------


function closest_vertex_ll(cellid, lat, lon)

! Return the vertex id of the closest one to the given point
! this version uses lat/lon.  see closest_vertex_xyz for the
! cartesian version.

integer,  intent(in)  :: cellid
real(r8), intent(in)  :: lat, lon
integer               :: closest_vertex_ll

real(r8) :: px, py, pz

! use the same radius as MPAS for computing this
call latlon_to_xyz(lat, lon, px, py, pz)

closest_vertex_ll = closest_vertex_xyz(cellid, px, py, pz)

end function closest_vertex_ll

!------------------------------------------------------------

function closest_vertex_xyz(cellid, px, py, pz)

! Return the vertex id of the closest one to the given point
! see closest_vertex_ll for the lat/lon version (which calls this)

integer,  intent(in)  :: cellid
real(r8), intent(in)  :: px, py, pz
integer               :: closest_vertex_xyz

integer :: nverts, i, vertexid
real(r8) :: distsq, closest_dist, dx, dy, dz

! nedges and nverts is same in a closed figure
nverts = nEdgesOnCell(cellid)

closest_dist = 1.0e38_r8   ! something really big; these are meters not radians
closest_vertex_xyz = -1

do i=1, nverts
   vertexid = verticesOnCell(i, cellid)
   dx = xVertex(vertexid) - px
   dy = yVertex(vertexid) - py
   dz = zVertex(vertexid) - pz
   distsq = (dx * dx) + (dy * dy) + (dz * dz)
   if (distsq < closest_dist) then
      closest_dist = distsq
      closest_vertex_xyz = vertexid
   endif
enddo

end function closest_vertex_xyz

!------------------------------------------------------------

subroutine make_edge_list_from_verts(nverts, vertex_list, nedges, edge_list)

! FIXME: will need a vertical level number input arg here to detect
! the boundary edges/vertices correctly.

! given a vertexid, look up the N cells which share it as
! a vertex, and then for those cells look up the edge ids.
! return a list with all edges listed exactly once (shared edges
! are detected and not replicated in the output list).
! the edge_list output should be at least 10x the Ncells to
! guarentee it will be large enough if all cells are disjoint.

integer, intent(in)  :: nverts, vertex_list(:)
integer, intent(out) :: nedges, edge_list(:)

integer :: edgecount, c, e, listlen, l, nextedge, v
integer :: ncells, cellid_list(3*nverts), v_base
logical :: found

! use the cellsOnVertex() array to find the three cells
! which share each input vertex.  the "edge list from cells"
! subroutine has an input 'degree' which controls
! how far from each vertex we go when collecting edges.
! if we need it, this subroutine could be changed to have
! the same input option.  right now it is hardcoded to have
! degree = 2
ncells = 0
do v=1, nverts
   ncells = ncells + 3
   v_base = (v-1) * 3
   cellid_list(v_base+1) = cellsOnVertex(1, vertex_list(v))
   cellid_list(v_base+2) = cellsOnVertex(2, vertex_list(v))
   cellid_list(v_base+3) = cellsOnVertex(3, vertex_list(v))
enddo

! use nEdgesOnCell(nCells) and edgesOnCell(nCells, 10) to
! find the edge numbers.  add them to the list, skipping if
! the edge is already there.  increment nedges each time a
! new edge is added.  check arrays for enough length before
! starting to work.

listlen = 0
do c=1, ncells
   edgecount = nEdgesOnCell(cellid_list(c))
   do e=1, edgecount
      nextedge = edgesOnCell(e, cellid_list(c))
      found = .false.
      addloop: do l=1, listlen
         if (edge_list(l) == nextedge) then
            found = .true.
            exit addloop
         endif
      enddo addloop
      if ( .not. found) then
         listlen = listlen + 1
         edge_list(listlen) = nextedge
      endif
   enddo
enddo

nedges = listlen

end subroutine make_edge_list_from_verts

!------------------------------------------------------------

subroutine make_edge_list_from_cells(ncells, cellids, nedges, edge_list)

! given a list of cellids, return a unique list of edges
! the edge_list output should be at least ncells * maxedges long

integer, intent(in)  :: ncells, cellids(:)
integer, intent(out) :: nedges, edge_list(:)

integer :: edgecount, c, e, listlen, l, nextedge
logical :: found

! use nEdgesOnCell(nCells) and edgesOnCell(nCells, 10) to
! find the edge numbers.  add them to the list, skipping if
! the edge is already there.  increment nedges each time a
! new edge is added.  check arrays for enough length before
! starting to work.

listlen = 0
do c=1, ncells
   edgecount = nEdgesOnCell(cellids(c))
   do e=1, edgecount
      nextedge = edgesOnCell(e, cellids(c))
      found = .false.
      addloop: do l=1, listlen
         if (edge_list(l) == nextedge) then
            found = .true.
            exit addloop
         endif
      enddo addloop
      if ( .not. found) then
         listlen = listlen + 1
         edge_list(listlen) = nextedge
      endif
   enddo
enddo

nedges = listlen

end subroutine make_edge_list_from_cells

!------------------------------------------------------------

subroutine merge_edge_lists(nedges1, edge_list1, nedges2, edge_list2)

! FIXME: will need a vertical level number input arg here to detect
! the boundary edges/vertices correctly.

! given 2 edge lists, add any edges in list 2 that are not already
! in list 2. edge_list1 needs to be at least nedges1+nedges2 in size
! to guarentee it will be large enough if edge sets are disjoint.

integer, intent(inout) :: nedges1, edge_list1(:)
integer, intent(in)    :: nedges2, edge_list2(:)

integer :: e, l, listlen
logical :: found


listlen = nedges1
do e=1, nedges2
   found = .false.
   addloop: do l=1, listlen
      if (edge_list1(l) == edge_list2(e)) then
         found = .true.
         exit addloop
      endif
   enddo addloop
   if ( .not. found) then
      listlen = listlen + 1
      edge_list1(listlen) = edge_list2(e)
   endif
enddo

nedges1 = listlen

end subroutine merge_edge_lists

!------------------------------------------------------------

subroutine make_cell_list(vertexid, degree, ncells, cell_list)

! FIXME: will need a vertical level number input arg here to detect
! the boundary edges/vertices correctly.

! given a vertexid and a degree, look up the N cells which
! share it as a vertex.   for degree 1, return the 3 cells
! which share this vertex.  for degree 2, return the second-nearest
! cells as well.  only handles degree 1 and 2.
! returns a unique list of cells; duplicates are removed.
! the cell_list output should be at least 10x the Ncells times
! the degree to guarentee it will be large enough.

integer, intent(in)  :: vertexid, degree
integer, intent(out) :: ncells, cell_list(:)

integer :: vertcount, c, c2, v, listlen, l, nextvert, nextcell
logical :: found

if (degree < 1 .or. degree > 3) then
   call error_handler(E_ERR, 'make_cell_list', 'internal error, must have 1 <= degree <= 3', &
                      source)
endif

! use the cellsOnVertex() array to find the three cells
! which share this vertex.  for degree 1, we are done.
! for degree 2, search the vertices of these cells as well.

cell_list(1) = cellsOnVertex(1, vertexid)
cell_list(2) = cellsOnVertex(2, vertexid)
cell_list(3) = cellsOnVertex(3, vertexid)

ncells = 3
if (degree == 1) return

listlen = ncells
do c=1, ncells
   vertcount = nEdgesOnCell(cell_list(c))
   do v=1, vertcount
      nextvert = verticesOnCell(v, cell_list(c))
      do c2=1, 3
         nextcell = cellsOnVertex(c2, nextvert)
         found = .false.
         addloop: do l=1, listlen
            if (cell_list(l) == nextcell) then
               found = .true.
               exit addloop
            endif
         enddo addloop
         if ( .not. found) then
            listlen = listlen + 1
            cell_list(listlen) = nextcell
         endif
      enddo
   enddo
enddo

ncells = listlen
if (degree == 2) return


! degree 3:
do c=1, ncells
   vertcount = nEdgesOnCell(cell_list(c))
   do v=1, vertcount
      nextvert = verticesOnCell(v, cell_list(c))
      do c2=1, 3
         nextcell = cellsOnVertex(c2, nextvert)
         found = .false.
         addloop2: do l=1, listlen
            if (cell_list(l) == nextcell) then
               found = .true.
               exit addloop2
            endif
         enddo addloop2
         if ( .not. found) then
            listlen = listlen + 1
            cell_list(listlen) = nextcell
         endif
      enddo
   enddo
enddo

ncells = listlen

end subroutine make_cell_list

!------------------------------------------------------------

subroutine move_pressure_to_edges(ncells, celllist, lower, upper, fract, nedges, edgelist, ier)

! given a list of n cell ids, and a list of edge ids, figure out which
! cells are on either side of each edge and compute new lower, upper and
! fract values.  the lower, upper and fract lists match the cells on the
! way in, they match the edges on the way out.

integer,  intent(in)    :: ncells, celllist(:)
integer,  intent(inout) :: lower(:, :), upper(:, :) ! ens_size
real(r8), intent(inout) :: fract(:, :) ! ens_size
integer,  intent(in)    :: nedges, edgelist(:)
integer,  intent(out)   :: ier(:)

!>@todo this is set in different places to different values
integer, parameter :: listsize = 30 
!real(r8) :: o_lower(listsize), o_fract(listsize)
real(r8), allocatable :: o_lower(:,:), o_fract(:,:)
real(r8) :: x1, x2, x
integer  :: i, c1, c2, c
integer :: length

length = size(upper, 2)

allocate(o_lower(listsize, length), o_fract(listsize, length))

! save the originals; we are going to overwrite these arrays
o_lower = lower
o_fract = fract

do i=1, nedges
   c1 = cellsOnEdge(1, edgelist(i))
   c2 = cellsOnEdge(2, edgelist(i))
   x1 = -1.0_r8
   x2 = -1.0_r8
   do c=1, ncells
      if (celllist(c) == c1) then
         x1 = o_lower(c, 1) + o_fract(c, 1)
      endif
      if (celllist(c) == c2) then
         x2 = o_lower(c, 1) + o_fract(c, 1)
      endif
   enddo
   if (x1 < 0.0_r8 .or. x2 < 0.0_r8) then
      write(string1, *) 'neither cell found in celllist', c1, c2, ' edge is ',  edgelist(i)
      call error_handler(E_ERR, 'move_pressure_to_edges', 'cell id not found in list', &
                         source, text2=string1)
   endif
   ! FIXME: should this be a log interpolation since we know this is going
   ! to be used for pressure only?  it's in the horizontal so the values
   ! shouldn't be too different, but still for consistency...
   !if (log_p_vert_interp) then
   !   x = exp((log(x1) + log(x2)) / 2.0_r8)
   ! else
       x = (x1 + x2) / 2.0_r8
   ! endif
   lower(i, :) = aint(x)
   upper(i, :) = lower(i, 1) + 1
   fract(i, :) = x - lower(i, 1)
enddo

ier = 0

end subroutine move_pressure_to_edges

!------------------------------------------------------------

function extrapolate_level(h, lb, ub)

! Given a test height and two level heights, extrapolate the
! actual level number which would correspond to that height.

real(r8), intent(in) :: h
real(r8), intent(in) :: lb
real(r8), intent(in) :: ub
real(r8)             :: extrapolate_level

real(r8) :: thickness, fract

thickness = ub - lb

if (h < lb) then
   ! extrapolate down - not in data values but in level
   ! relative to the height difference in levels 1-2
   extrapolate_level = lb - (((lb - h) + lb) / thickness)
   
else if (h > ub) then
   ! extrapolate up - not in data values but in level
   ! relative to the height difference in levels nb-1 to nb
   extrapolate_level = ub + ((ub - (h - ub)) / thickness)

else
   write(string1, *) h, ' not outside range of ', lb, ' and ', ub
   call error_handler(E_ERR, 'extrapolate_level', &
                      'extrapolate code called for height not out of range', &
                      source, text2=string1)
endif

end function extrapolate_level

!------------------------------------------------------------

subroutine latlon_to_xyz(lat, lon, x, y, z)

! Given a lat, lon in degrees, return the cartesian x,y,z coordinate
! on the surface of a specified radius relative to the origin
! at the center of the earth.  (this radius matches the one
! used at MPAS grid generation time and must agree in order
! to be consistent with the cartisian coordinate arrays in
! the MPAS data files.)

real(r8), intent(in)  :: lat, lon
real(r8), intent(out) :: x, y, z

real(r8) :: rlat, rlon

rlat = lat * deg2rad
rlon = lon * deg2rad

x = radius * cos(rlon) * cos(rlat)
y = radius * sin(rlon) * cos(rlat)
z = radius * sin(rlat)

end subroutine latlon_to_xyz

!------------------------------------------------------------

subroutine inside_triangle(t1, t2, t3, r, lat, lon, inside, weights)

! given 3 corners of a triangle and an xyz point, compute whether
! the point is inside the triangle.  this assumes r is coplanar
! with the triangle - the caller must have done the lat/lon to
! xyz conversion with a constant radius and then this will be
! true (enough).  sets inside to true/false, and returns the
! weights if true.  weights are set to 0 if false.

real(r8), intent(in)  :: t1(3), t2(3), t3(3)
real(r8), intent(in)  :: r(3), lat, lon
logical,  intent(out) :: inside
real(r8), intent(out) :: weights(3)

! check for degenerate cases first - is the test point located
! directly on one of the vertices?  (this case may be common
! if we're computing on grid point locations.)
if (all(abs(r - t1) < roundoff)) then
   inside = .true.
   weights = (/ 1.0_r8, 0.0_r8, 0.0_r8 /)
   return
else if (all(abs(r - t2) < roundoff)) then
   inside = .true.
   weights = (/ 0.0_r8, 1.0_r8, 0.0_r8 /)
   return
else if (all(abs(r - t3) < roundoff)) then
   inside = .true.
   weights = (/ 0.0_r8, 0.0_r8, 1.0_r8 /)
   return
endif

! not a vertex. compute the weights.  if any are
! negative, the point is outside.  since these are
! real valued computations define a lower bound for
! numerical roundoff error and be sure that the
! weights are not just *slightly* negative.
call get_3d_weights(r, t1, t2, t3, lat, lon, weights)

if (any(weights < -roundoff)) then
   inside = .false.
   weights = 0.0_r8
   return
endif

! truncate barely negative values to 0
inside = .true.
where (weights < 0.0_r8) weights = 0.0_r8
return

end subroutine inside_triangle

!------------------------------------------------------------

subroutine uv_cell_to_edges(zonal_wind, meridional_wind, uedge, full_u)

! Project u, v wind (increments) at cell centers onto the edges.
! FIXME:
!        we can hard-code R3 here since it comes from the (3d) x/y/z cartesian coordinate.
!        We read cellsOnEdge and edgeNormalVectors in get_grid.
!        Here "uedge" is an edge normal wind which can be either a full field or an increment
!        depending on the optional input (full_u).
!        if full_u = .true., then the edge wind is replaced by averaging zonal and meridional 
!        winds at cell centers.
!        if full_u does not exist, uedge is the analysis increment from the updated cell-center winds.
!        We do not update/compute uedge in the outermost edge in the regional MPAS.
! This routine followed the updating part in tend_toEdges in 
! MPAS/src/core_atmosphere/physics/mpas_atmphys_todynamics.F.

real(r8), intent(in) :: zonal_wind(:,:)             ! u wind updated from filter
real(r8), intent(in) :: meridional_wind(:,:)        ! v wind updated from filter
real(r8), intent(inout):: uedge(:,:)                ! normal velocity (increment) on the edges
logical,  intent(in), optional :: full_u            ! compute a full field, not an increment

! Local variables
integer, parameter :: R3 = 3
real(r8) :: east(R3,nCells), north(R3,nCells)
real(r8) :: lonCell_rad(nCells), latCell_rad(nCells)
integer  :: iCell, iEdge, cell1, cell2

if ( .not. module_initialized ) call static_init_model

! Initialization for increments
if ( .not. present(full_u)) then
   uedge(:,:) = 0.0_r8
endif

! Back to radians (locally)
lonCell_rad = lonCell*deg2rad
latCell_rad = latCell*deg2rad

! Compute unit vectors in east and north directions for each cell:
do iCell = 1, nCells

    east(1,iCell) = -sin(lonCell_rad(iCell))
    east(2,iCell) =  cos(lonCell_rad(iCell))
    east(3,iCell) =  0.0_r8
    call r3_normalize(east(1,iCell), east(2,iCell), east(3,iCell))

    north(1,iCell) = -cos(lonCell_rad(iCell))*sin(latCell_rad(iCell))
    north(2,iCell) = -sin(lonCell_rad(iCell))*sin(latCell_rad(iCell))
    north(3,iCell) =  cos(latCell_rad(iCell))
    call r3_normalize(north(1,iCell), north(2,iCell), north(3,iCell))

enddo

! Project reconstructed winds from the cell centers to the edges
do iEdge = 1, nEdges
   if(.not.on_outermost_edge(iEdge)) then ! we do not update the outermost edge
      cell1 = cellsOnEdge(1,iEdge)
      cell2 = cellsOnEdge(2,iEdge)
      if((.not.on_outermost_cell(cell1)) .and. (.not.on_outermost_cell(cell2))) then 
      uedge(:,iEdge) = zonal_wind(:,cell1)      * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell1)   &
                                                      +  edgeNormalVectors(2,iEdge) * east(2,cell1)   &
                                                      +  edgeNormalVectors(3,iEdge) * east(3,cell1))  &
                     + meridional_wind(:,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell1)   &
                                                      +  edgeNormalVectors(2,iEdge) * north(2,cell1)   &
                                                      +  edgeNormalVectors(3,iEdge) * north(3,cell1))  &
                     + zonal_wind(:,cell2)      * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell2)   &
                                                      +  edgeNormalVectors(2,iEdge) * east(2,cell2)   &
                                                      +  edgeNormalVectors(3,iEdge) * east(3,cell2))  &
                     + meridional_wind(:,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell2)   &
                                                      +  edgeNormalVectors(2,iEdge) * north(2,cell2)   &
                                                      +  edgeNormalVectors(3,iEdge) * north(3,cell2))
      endif ! if((.not.on_outermost_cell(cell1)) .and. (.not.on_outermost_cell(cell2)))
   endif    ! if(.not.on_outermost_edge(iEdge)) then
enddo       ! do iEdge = 1, nEdges

end subroutine uv_cell_to_edges


!------------------------------------------------------------------

!==================================================================
! The following (private) routines were borrowed from the MPAS code
!==================================================================

!------------------------------------------------------------------

subroutine r3_normalize(ax, ay, az)

!normalizes the vector (ax, ay, az)

real(r8), intent(inout) :: ax, ay, az
real(r8) :: mi

 mi = 1.0_r8 / sqrt(ax**2 + ay**2 + az**2)
 ax = ax * mi
 ay = ay * mi
 az = az * mi

end subroutine r3_normalize


!------------------------------------------------------------------

function theta_to_tk (ens_size, theta, rho, qv, istatus) 

! Compute sensible temperature [K] from potential temperature [K].
! code matches computation done in MPAS model

integer,                       intent(in)  :: ens_size
real(r8), dimension(ens_size), intent(in)  :: theta    ! potential temperature [K]
real(r8), dimension(ens_size), intent(in)  :: rho      ! dry air density [kg/m3]
real(r8), dimension(ens_size), intent(in)  :: qv       ! water vapor mixing ratio [kg/kg]
integer,  dimension(ens_size), intent(inout) :: istatus
real(r8), dimension(ens_size) :: theta_to_tk          ! sensible temperature [K]

! Local variables
real(r8), dimension(ens_size) :: theta_m    ! potential temperature modified by qv
real(r8), dimension(ens_size) :: exner      ! exner function
real(r8), dimension(ens_size) :: qv_nonzero ! qv >= 0
integer :: e

qv_nonzero = max(qv,0.0_r8)
theta_to_tk = missing_r8

where (istatus == 0)

   theta_m = (1.0_r8 + rvord * qv_nonzero)*theta
   
   where (theta_m > 0.0_r8 .and. rho > 0.0_r8)  ! Check if all the input are positive

      exner = ( (rgas/p0) * (rho*theta_m) )**rcv
   
      ! Temperature [K]
      theta_to_tk = theta * exner

   elsewhere

      istatus = TK_COMPUTATION_ERROR

   endwhere

endwhere


end function theta_to_tk


!------------------------------------------------------------------

subroutine compute_full_pressure(ens_size, theta, rho, qv, pressure, tk, istatus)

! Compute full pressure from the equation of state.
! since it has to compute sensible temp along the way,
! make temp one of the return values rather than having
! to call theta_to_tk() separately.
! code matches computation done in MPAS model

integer,  intent(in)  :: ens_size
real(r8), dimension(ens_size), intent(in)  :: theta    ! potential temperature [K]
real(r8), dimension(ens_size), intent(in)  :: rho      ! dry air density [kg/m3]
real(r8), dimension(ens_size), intent(in)  :: qv       ! water vapor mixing ratio [kg/kg]
real(r8), dimension(ens_size), intent(out) :: pressure ! full pressure [Pa]
real(r8), dimension(ens_size), intent(out) :: tk       ! return sensible temperature to caller
integer,  dimension(ens_size), intent(inout):: istatus

! Local variables
real(r8), dimension(ens_size) :: qv_nonzero
integer   :: e

pressure = missing_r8
qv_nonzero = max(qv,0.0_r8)

tk = theta_to_tk(ens_size, theta, rho, qv_nonzero, istatus)

where (istatus == 0)       ! We only take non-missing tk here
   pressure = rho * rgas * tk * (1.0_r8 + rvord * qv_nonzero)
end where

end subroutine compute_full_pressure

!--------------------------------------------------------------------
!> pass the vertical localization coordinate to assim_tools_mod
function query_vert_localization_coord()

integer :: query_vert_localization_coord

query_vert_localization_coord = vert_localization_coord

end function query_vert_localization_coord

!--------------------------------------------------------------------
!> read the time from the input file
function read_model_time(filename)

character(len=256), intent(in) :: filename

type(time_type) :: read_model_time
integer         :: ncid  ! netcdf file id
integer         :: ret ! return code for netcdf

ret = nf90_open(filename, NF90_NOWRITE, ncid)
call nc_check(ret, 'opening', filename)

read_model_time = get_analysis_time_ncid(ncid, filename)

ret = nf90_close(ncid)
call nc_check(ret, 'closing', filename)


end function read_model_time

!----------------------------------------------------------------------
! Returns integers taken from tstring
! It is assumed that the tstring char array is as YYYY-MM-DD_hh:mm:ss

subroutine set_wrf_date (tstring, year, month, day, hour, minute, second)

integer,           intent(in) :: year, month, day, hour, minute, second
character(len=TIMELEN), intent(out)  :: tstring

character(len=4)  :: ch_year
character(len=2)  :: ch_month, ch_day, ch_hour, ch_minute, ch_second

write(ch_year,'(i4)') year
write(ch_month,'(i2)') month
if (ch_month(1:1) == " ") ch_month(1:1) = "0"
write(ch_day,'(i2)') day
if (ch_day(1:1) == " ") ch_day(1:1) = "0"
write(ch_hour,'(i2)') hour
if (ch_hour(1:1) == " ") ch_hour(1:1) = "0"
write(ch_minute,'(i2)') minute
if (ch_minute(1:1) == " ") ch_minute(1:1) = "0"
write(ch_second,'(i2)') second
if (ch_second(1:1) == " ") ch_second(1:1) = "0"
tstring(1:4)   = ch_year
tstring(5:5)   = "-"
tstring(6:7)   = ch_month
tstring(8:8)   = "-"
tstring(9:10)  = ch_day
tstring(11:11) = "_"
tstring(12:13) = ch_hour
tstring(14:14) = ":"
tstring(15:16) = ch_minute
tstring(17:17) = ":"
tstring(18:19) = ch_second

end subroutine set_wrf_date

end module model_mod
