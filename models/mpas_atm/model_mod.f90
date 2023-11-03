! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module model_mod

! MPAS Atmosphere model interface to the DART data assimilation system.
! Code in this module is compiled with the DART executables.  It isolates
! all information about the MPAS grids, model variables, and other details.
! There are a set of 16 subroutine interfaces that are required by DART;
! these cannot be changed.  Additional public routines in this file can
! be used by converters and utilities and those interfaces can be anything
! that is useful to other pieces of code.

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

use    utilities_mod, only : register_module, error_handler, get_unit,         &
                             E_ERR, E_WARN, E_MSG, E_ALLMSG, logfileunit,      &
                             do_output, to_upper, nmlfileunit,                 &
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

use distributed_state_mod

! netcdf modules
use typesizes
use netcdf

! RBF (radial basis function) modules, donated by LANL. currently deprecated
! in this version.  they did the job but not as well as other techniques and
! at a much greater execution-time code.  they were used to interpolate
! values at arbitrary locations, not just at cell centers.  with too small
! a set of basis points the values were discontinuous at cell boundaries;
! with too many the values were too smoothed out.  we went back to
! barycentric interpolation in triangles formed by the three cell centers
! that enclosed the given point.
use get_geometry_mod
use get_reconstruct_mod

use state_structure_mod, only :  add_domain, get_model_variable_indices, &
                                 state_structure_info, get_index_start, get_index_end, &
                                 get_num_variables, get_domain_size, get_varid_from_varname, &
                                 get_variable_name, get_num_dims, get_dim_lengths, &
                                 get_dart_vector_index, get_num_varids_from_kind


implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.

! routines in this list have code in this module
public :: get_model_size,                 &
          get_num_vars,                   &
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
          analysis_file_to_statevector, &
          statevector_to_analysis_file, &
          statevector_to_boundary_file, &
          get_analysis_time,            &
          get_grid_dims,                &
          get_xland,                    &
          get_surftype,                 &
          get_cell_center_coords,       &
          get_bdy_mask,                 &
          print_variable_ranges,        &
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
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

! module global storage; maintains values between calls, accessible by
! any subroutine
character(len=512) :: string1, string2, string3
character(len=256) :: locstring
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
integer            :: xyzdebug = 0
integer            :: debug = 0   ! turn up for more and more debug messages

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
! evaluated without error.  if extrapolate is true, extrapolate from the
! first or last model level.  if extrapolate is false, use level 1 or N.
real(r8) :: outside_grid_level_tolerance = -1.0_r8
logical  :: extrapolate = .false.

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
   debug,                        &
   xyzdebug,                     &
   use_u_for_wind,               &
   use_rbf_option,               &
   update_u_from_reconstruct,    &
   use_increments_for_u_update,  &
   highest_obs_pressure_mb,      &
   outside_grid_level_tolerance, &
   extrapolate,                  &
   sfc_elev_max_diff,            &
   write_grid_to_diag_files,     &
   no_normalization_of_scale_heights

! DART state vector contents are specified in the input.nml:&mpas_vars_nml namelist.
integer, parameter :: max_state_variables = 80
integer, parameter :: num_state_table_columns = 2
integer, parameter :: num_bounds_table_columns = 4
character(len=NF90_MAX_NAME) :: mpas_state_variables(max_state_variables * num_state_table_columns ) = ' '
character(len=NF90_MAX_NAME) :: mpas_state_bounds(num_bounds_table_columns, max_state_variables ) = ' '
character(len=NF90_MAX_NAME) :: variable_table(max_state_variables, num_state_table_columns )

! this is special and not in the namelist.  boundary files have a fixed
! set of variables with fixed names.
character(len=NF90_MAX_NAME) :: lbc_variables(max_state_variables) = ''

namelist /mpas_vars_nml/ mpas_state_variables, mpas_state_bounds

integer :: nfields

!>@todo FIXME - REMOVE AS MUCH OF THIS AS POSSIBLE.
!> some of this information is in the state structure now.
!> the duplicate progvar stuff should be removed and the 
!> state routines used instead.  this duplicates work and 
!> makes us keep up code in 2 different places.

! original code:
! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME), dimension(NF90_MAX_VAR_DIMS) :: dimname
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer :: xtype         ! netCDF variable type (NF90_double, etc.)
   integer :: numdims       ! number of dims - excluding TIME
   integer :: numvertical   ! number of vertical levels in variable
   integer :: numcells      ! number of horizontal locations (cell centers)
   integer :: numedges      ! number of horizontal locations (edges for velocity components)
   logical :: ZonHalf       ! vertical coordinate has dimension nVertLevels
   integer :: varsize       ! prod(dimlens(1:numdims))
   integer :: index1        ! location in dart state vector of first occurrence
   integer :: indexN        ! location in dart state vector of last  occurrence
   integer :: dart_kind
   character(len=obstypelength) :: kind_string
   logical  :: clamping     ! does variable need to be range-restricted before
   real(r8) :: range(2)     ! being stuffed back into MPAS analysis file.
   logical  :: out_of_range_fail  ! is out of range fatal if range-checking?
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

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

! common names that call specific subroutines based on the arg types
INTERFACE vector_to_prog_var
      MODULE PROCEDURE vector_to_1d_prog_var
      MODULE PROCEDURE vector_to_2d_prog_var
      MODULE PROCEDURE vector_to_3d_prog_var
END INTERFACE

INTERFACE prog_var_to_vector
      MODULE PROCEDURE prog_var_1d_to_vector
      MODULE PROCEDURE prog_var_2d_to_vector
      MODULE PROCEDURE prog_var_3d_to_vector
END INTERFACE

INTERFACE get_analysis_time
      MODULE PROCEDURE get_analysis_time_ncid
      MODULE PROCEDURE get_analysis_time_fname
END INTERFACE

INTERFACE get_index_range
      MODULE PROCEDURE get_index_range_int
      MODULE PROCEDURE get_index_range_string
END INTERFACE


interface write_model_time
   module procedure write_model_time_file
   module procedure write_model_time_restart
end interface


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

!==================================================================
! All the public REQUIRED interfaces come first - just by convention.
!==================================================================


!------------------------------------------------------------------

subroutine static_init_model()
!>@todo FIXME - can we add an optional arg here for update_bc use?

! Called to do one time initialization of the model.
!
! All the grid information comes from the initialization of
! the dart_model_mod module.

! Local variables - all the important ones have module scope


integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname,dimname
character(len=obstypelength)       :: kind_string
integer :: ncid, VarID, numdims, varsize, dimlen
integer :: iunit, io, ivar, i, index1, indexN, iloc, kloc
integer :: ss, dd, z1, m1
integer :: nDimensions, nVariables, nAttributes, unlimitedDimID, TimeDimID
integer :: cel1, cel2, lbc_nfields
logical :: both
real(r8) :: variable_bounds(max_state_variables, 2)
integer :: variable_qtys(max_state_variables)
character(len=*), parameter :: routine = 'static_init_model'

if ( module_initialized ) return ! only need to do this once.

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

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
! Set the time step ... causes mpas namelists to be read.
! Ensures model_timestep is multiple of 'dynamics_timestep'

call set_calendar_type( calendar )   ! comes from model_mod_nml

model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation window is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,routine,string1,source,revision,revdate)

!---------------------------------------------------------------
! 1) get grid dimensions
! 2) allocate space for the grids
! 3) read them from the analysis file

ncid = nc_open_file_readonly(init_template_filename, routine)

! move this up - some of the routines below depend on having
! the variable table parsed already, and U added if needed.
call verify_state_variables( mpas_state_variables, ncid, init_template_filename, &
                             nfields, variable_table)

! get sizes
call read_grid_dims(ncid)

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
            call error_handler(E_ERR,'static_init_model', string1, source, revision, revdate)
         endif
      enddo
   enddo
endif

!---------------------------------------------------------------
! Compile the list of model variables to use in the creation
! of the DART state vector.
!
! THIS CODE SHOULD BE REMOVED - it is done by the add_domain code.
!


TimeDimID = FindTimeDimension( ncid )

if (TimeDimID < 0 ) then
   write(string1,*)'unable to find a dimension named Time.'
   call error_handler(E_MSG,'static_init_model', string1, source, revision, revdate)
endif

call nc_check(nf90_Inquire(ncid,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                    'static_init_model', 'inquire '//trim(init_template_filename))

if ( (TimeDimID > 0) .and. (unlimitedDimID > 0) .and. (TimeDimID /= unlimitedDimID)) then
   write(string1,*)'IF Time is not the unlimited dimension, I am lost.'
   call error_handler(E_MSG,'static_init_model', string1, source, revision, revdate)
endif

index1  = 1;
indexN  = 0;
do ivar = 1, nfields

   varname                   = trim(variable_table(ivar,1))
   kind_string               = trim(variable_table(ivar,2))
   progvar(ivar)%varname     = varname
   progvar(ivar)%kind_string = kind_string
   progvar(ivar)%dart_kind   = get_index_for_quantity( progvar(ivar)%kind_string )
   progvar(ivar)%numdims     = 0
   progvar(ivar)%numvertical = 1
   progvar(ivar)%dimlens     = MISSING_I
   progvar(ivar)%numcells    = MISSING_I
   progvar(ivar)%numedges    = MISSING_I

   string2 = trim(init_template_filename)//' '//trim(varname)

   call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), &
            'static_init_model', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid, VarID, xtype=progvar(ivar)%xtype, &
           dimids=dimIDs, ndims=numdims), 'static_init_model', 'inquire '//trim(string2))

   ! Since we are not concerned with the TIME dimension, we need to skip it.
   ! When the variables are read, only a single timestep is ingested into
   ! the DART state vector.

   varsize = 1
   dimlen  = 1
   DimensionLoop : do i = 1,numdims

      if (dimIDs(i) == TimeDimID) cycle DimensionLoop

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen, name=dimname), &
                                          'static_init_model', string1)

      progvar(ivar)%numdims    = progvar(ivar)%numdims + 1
      progvar(ivar)%dimlens(i) = dimlen
      progvar(ivar)%dimname(i) = trim(dimname)
      varsize = varsize * dimlen

      select case ( dimname(1:6) )
         case ('nCells')
            progvar(ivar)%numcells = dimlen
         case ('nEdges')
            progvar(ivar)%numedges = dimlen
         case ('nVertL')  ! nVertLevels, nVertLevelsP1, nVertLevelsP2
            progvar(ivar)%numvertical = dimlen
         case ('nSoilL')  ! nSoilLevels
            progvar(ivar)%numvertical = dimlen
      end select

   enddo DimensionLoop

   ! this call sets: clamping, bounds, and out_of_range_fail in the progvar entry
   call get_variable_bounds(mpas_state_bounds, ivar)

   if (progvar(ivar)%numvertical == nVertLevels) then
      progvar(ivar)%ZonHalf = .TRUE.
   else
      progvar(ivar)%ZonHalf = .FALSE.
   endif

   if (varname == 'u') has_edge_u = .true.
   if (varname == 'uReconstructZonal' .or. &
       varname == 'uReconstructMeridional') has_uvreconstruct = .true.

   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1
   index1                    = index1 + varsize      ! sets up for next variable

   if ( debug > 11 .and. do_output()) call dump_progvar(ivar)

enddo

call nc_close_file(ncid, routine)


! FIXME: moved below to after add_domain() calls
!model_size = progvar(nfields)%indexN

if ( debug > 0 .and. do_output()) then
  write(logfileunit,*)
  write(     *     ,*)
  write(logfileunit,'(" static_init_model: nCells, nEdges, nVertices, nVertLevels =",4(1x,i9))') &
                                          nCells, nEdges, nVertices, nVertLevels
  write(     *     ,'(" static_init_model: nCells, nEdges, nVertices, nVertLevels =",4(1x,i9))') &
                                          nCells, nEdges, nVertices, nVertLevels
!  write(logfileunit, *)'static_init_model: model_size = ', model_size
!  write(     *     , *)'static_init_model: model_size = ', model_size
  if ( global_grid ) then
     write(logfileunit, *)'static_init_model: grid is a global grid '
     write(     *     , *)'static_init_model: grid is a global grid '
  else
     write(logfileunit, *)'static_init_model: grid is NOT a global grid. Lateral boundaries exist '
     write(     *     , *)'static_init_model: grid is NOT a global grid. Lateral boundaries exist '
  endif
  if ( all_levels_exist_everywhere ) then
     write(logfileunit, *)'static_init_model: all cells have same number of vertical levels '
     write(     *     , *)'static_init_model: all cells have same number of vertical levels '
  else
     write(logfileunit, *)'static_init_model: cells have varying number of vertical levels '
     write(     *     , *)'static_init_model: cells have varying number of vertical levels '
  endif
endif


! set the domain(s) in the state structure here

variable_bounds(1:nfields, 1) = progvar(1:nfields)%range(1)
variable_bounds(1:nfields, 2) = progvar(1:nfields)%range(2)
variable_qtys(1:nfields) = progvar(1:nfields)%dart_kind

anl_domid = add_domain(init_template_filename, nfields,           &
                       var_names  = variable_table (1:nfields,1), &
                       kind_list  = variable_qtys(1:nfields),     &
                       clamp_vals = variable_bounds(1:nfields,:) )

model_size = get_domain_size(anl_domid)
if ( debug > 4 .and. do_output()) print*,'model_size(anl_domid)=',model_size  ! HA

lbc_nfields = 0

! if we have a lateral boundary file, add it to the domain
! so we have access to the corresponding lbc_xxx fields.
!>@todo FIXME: if we want to do increments, we could also add a
! third domain which is the original forecast fields before
! the assimilation (so we can compute increments)
if (.not. global_grid .and. lbc_variables(1) /= '') then
   ! regional: count number of lbc fields to read in
   COUNTUP: do i=1, max_state_variables
      if (lbc_variables(i) /= '') then
         lbc_nfields = lbc_nfields + 1
      else
         exit COUNTUP
      endif
   enddo COUNTUP
   if( debug > 4 .and. do_output()) print*, 'Regional: number of lbc fields to read in = ', lbc_nfields
   lbc_domid = add_domain(bdy_template_filename, lbc_nfields,    &
                          var_names  = lbc_variables)
                        ! FIXME clamp_vals = variable_bounds(1:nfields,:) )
   if( debug > 4 .and. do_output()) print*, 'model_size, lbc_domid =',model_size,lbc_domid

   model_size = model_size + get_domain_size(lbc_domid)
else
   lbc_domid = -1
endif

if ( debug > 4 .and. do_output()) then
     call state_structure_info(anl_domid)
     print*, 'model_size =',model_size
     if(lbc_domid >= 0) print*, 'get_domain_size(lbc_domid):',get_domain_size(lbc_domid)
endif
!if ( debug > 4 .and. do_output() .and. lbc_domid >= 0) call state_structure_info(lbc_domid)

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
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate, &
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
   call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate, &
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
   call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate, &
                      text2=string2, text3=string3)
endif

! basically we cannot do much without having at least these
! three fields in the state vector.  refuse to go further
! if these are not present:
!print *, get_num_varids_from_kind(anl_domid, QTY_POTENTIAL_TEMPERATURE)
!print *, get_num_varids_from_kind(anl_domid, QTY_DENSITY)
!print *, get_num_varids_from_kind(anl_domid, QTY_VAPOR_MIXING_RATIO)
if ((get_num_varids_from_kind(anl_domid, QTY_POTENTIAL_TEMPERATURE) < 0) .or. &
    (get_num_varids_from_kind(anl_domid, QTY_DENSITY) < 0) .or. &
    (get_num_varids_from_kind(anl_domid, QTY_VAPOR_MIXING_RATIO) < 0)) then
   write(string1, *) 'State vector is missing one or more of the following fields:'
   write(string2, *) 'Potential Temperature (theta), Density (rho), Vapor Mixing Ratio (qv).'
   write(string3, *) 'Cannot convert between height/pressure nor compute sensible temps.'
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate, &
                      text2=string2, text3=string3)
endif


if (extrapolate) then
   call error_handler(E_MSG,'static_init_model',&
                      'extrapolate not supported yet; will use level 1 or N values', &
                      source, revision, revdate)
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

subroutine get_state_meta_data(index_in, location, var_type)

! given an index into the state vector, return its location and
! if given, the var kind.   despite the name, var_type is a generic
! kind, like those in obs_kind/obs_kind_mod.f90, starting with QTY_

! passed variables
integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

! Local variables

integer  :: nzp, iloc, vloc, nf, ndim
real(r8) :: height
type(location_type) :: new_location

if ( .not. module_initialized ) call static_init_model

! get the local indicies and type from dart index. kloc is a dummy variable for this subroutine

call find_mpas_indices(index_in, iloc, vloc, ndim, nf)

nzp  = progvar(nf)%numvertical

! the zGrid array contains the location of the cell top and bottom faces, so it has one
! more value than the number of cells in each column.  for locations of cell centers
! you have to take the midpoint of the top and bottom face of the cell.
if (progvar(nf)%numedges /= MISSING_I) then
   if (.not. data_on_edges) then
      call error_handler(E_ERR, 'get_state_meta_data', &
                        'Internal error: numedges present but data_on_edges false', &
                        source, revision, revdate, text2='variable '//trim(progvar(nf)%varname))
   endif
   if ( progvar(nf)%ZonHalf ) then
      height = zGridEdge(vloc,iloc)
   else
      call error_handler(E_ERR, 'get_state_meta_data', 'no support for edges at face heights', &
                         source, revision, revdate)
   endif

   if (nzp <= 1) then
      location = set_location(lonEdge(iloc),latEdge(iloc), height, VERTISSURFACE)
   else
      location = set_location(lonEdge(iloc),latEdge(iloc), height, VERTISHEIGHT)
   endif
else
   if ( progvar(nf)%ZonHalf ) then
      height = zGridCenter(vloc,iloc)
   else if (nzp <= 1) then
      height = zGridFace(1,iloc)
   else
      height = zGridFace(vloc,iloc)
   endif

   if (nzp <= 1) then
      location = set_location(lonCell(iloc),latCell(iloc), height, VERTISSURFACE)
   else
      location = set_location(lonCell(iloc),latCell(iloc), height, VERTISHEIGHT)
   endif
endif

! we are no longer passing in an ensemble handle to this routine, so we
! cannot do vertical conversion here.  assim_tools will call vertical conversion
! on the obs and on the state.

if (debug > 20 .and. do_output()) then

    write(*,'("INDEX_IN / IVAR : ",(i10,2x),(i5,2x))') index_in, nf
    write(*,'("                 ILOC, VLOC: ",2(i5,2x))') iloc, vloc
    write(*,'("                 LON/LAT/HGT: ",3(f12.3,2x))') lonCell(iloc), latCell(iloc), height

endif

if (present(var_type)) then
   var_type = progvar(nf)%dart_kind
endif

end subroutine get_state_meta_data

!------------------------------------------------------------------
!> given an index into the state vector, return with the cellid
!> and the vertical level if this is a 2d variable.  also return
!> the dimensionality, and optionally the progvar index.

subroutine find_mpas_indices(index_in, cellid, vert_level, ndim, nf)
integer(i8), intent(in)  :: index_in
integer,     intent(out) :: cellid
integer,     intent(out) :: vert_level
integer,     intent(out), optional :: ndim
integer,     intent(out), optional :: nf

integer  :: i, j, k ! Indices into variable (note k is not used in MPAS)
integer  :: nzp, iloc, vloc, nnf

! get the local indicies and type from dart index. 'k' is a dummy variable for this subroutine

call get_model_variable_indices(index_in, i, j, k, var_id=nnf)

if (progvar(nnf)%numdims == 2) then     ! variable(vcol, iloc)
   vert_level = i
   cellid = j
elseif (progvar(nnf)%numdims == 1) then ! variable(iloc)
   cellid = i
   vert_level = 1
else
   call error_handler(E_ERR, 'find_mpas_indices ', 'expecting 1D or 2D variable')
endif

if (present(ndim)) ndim = progvar(nnf)%numdims
if (present(nf)) nf = nnf

end subroutine find_mpas_indices

!------------------------------------------------------------------
!> given a domain and varid, return 3 dimensions - setting to 1 if
!> not present.

subroutine find_mpas_dims(domid, ivar, ndims, dims)
integer, intent(in)  :: domid
integer, intent(in)  :: ivar
integer, intent(out) :: ndims
integer, intent(out) :: dims(3)

! start with everything 1, then set values which are different
dims(:) = 1

! get the domain info
ndims = get_num_dims(domid, ivar)
dims(1:ndims) = get_dim_lengths(domid, ivar)

end subroutine find_mpas_dims

!------------------------------------------------------------------
!> given a state vector, a location, and a QTY_xxx, return the
!> interpolated value at that location, and an error code.  0 is success,
!> anything positive is an error.  (negative reserved for system use)
!>
!>       ERROR codes:
!>
!>       ISTATUS =  1:  general error for rttov - at least one of the input variables goes wrong
!>       ISTATUS = 99:  general error in case something terrible goes wrong...
!>       ISTATUS = 81:  Vertical location too high
!>       ISTATUS = 80:  Vertical location too low
!>       ISTATUS = 88:  this kind is not in the state vector
!>       ISTATUS = 89:  tk cannot be computed.
!>       ISTATUS = 11:  Could not find the closest cell center that contains this lat/lon
!>       ISTATUS = 12:  Surface obs too far away from model elevation
!>       ISTATUS = 13:  Missing value in interpolation.
!>       ISTATUS = 14:  Could not find the other two cell centers of the triangle that contains this lat/lon
!>       ISTATUS = 15:  Cell centers of the triangle fall in the lateral boundary zone
!>       ISTATUS = 16:  Don't know how to do vertical velocity for now
!>       ISTATUS = 17:  Unable to compute pressure values
!>       ISTATUS = 18:  altitude illegal
!>       ISTATUS = 19:  could not compute u using RBF code
!>       ISTATUS = 101: Internal error; reached end of subroutine without
!>                      finding an applicable case.
!>       ISTATUS = 201: Reject observation from user specified pressure level
!>       ISTATUS = 988: pressure is not monotonically descreased with level.
!>
!> Debugging options:
!> 0: No prints, whatsoever.
!> 1: Print only when each obs is rejected.
!> 2: Print the basic info on each obs.
!> 3: Print the base location of the localization.
!> 5: Print info on the localized obs and the model grids
!> 9: Print info on wind updates
!> 10: Print info on each obs in detail.
!> 11: Print info on localized states
!> 12: Print info on each ensemble member in detail.
!> 13: Print detailed info on triangle and level search


subroutine model_interpolate(state_handle, ens_size, location, obs_type, expected_obs, istatus)

! passed variables

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_type
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

! local storage

type(location_type) :: location_tmp(ens_size)
integer  :: ivar, obs_kind
integer  :: tvars(3)
integer  :: cellid
logical  :: goodkind, surface_obs
real(r8) :: lpres(ens_size), values(3, ens_size)
real(r8) :: llv(3)    ! lon/lat/vert
integer  :: e, verttype

if ( .not. module_initialized ) call static_init_model

expected_obs = MISSING_R8
istatus      = 0          ! must be positive (and integer)

! rename for sanity - we can't change the argument names
! to this subroutine, but this really is a kind.
obs_kind = obs_type

! write the location information into a string for error messages and debugging
call write_location(0,location,charstring=locstring)

llv = get_location(location)
verttype = nint(query_location(location))
surface_obs = (verttype == VERTISSURFACE) 

! this routine returns the cellid for a global mpas grid, same as
! find_closest_cell_center().
! for a regional grid it only returns a good cellid if the closest cell center
! AND the other 2 triangle points surrounding this location are completely inside 
! the grid and none of the vertices are in the boundary region.
cellid = cell_ok_to_interpolate(location)
if (cellid < 1) then
   !print *, 'model_interpolate: lon/lat is outside the domain: ', llv(1), llv(2)
   istatus = 11
   goto 100
endif

if (debug > 1 .and. do_output()) &
   print *, 'model_interpolate for obs_kind:',&
             trim(get_name_for_quantity(obs_kind)),' at ',trim(locstring),' cellid:',cellid

! FIXME see issue #96 - remove all but surface elevation
! pass surface variables for rttov - it should be ok as these are diagnostic (except for surface elevation).
if((obs_kind == QTY_SURFACE_PRESSURE) .or. (obs_kind == QTY_SURFACE_ELEVATION)    .or. &
   (obs_kind == QTY_2M_TEMPERATURE)   .or. (obs_kind == QTY_2M_SPECIFIC_HUMIDITY) .or. &
   (obs_kind == QTY_SKIN_TEMPERATURE) .or. (obs_kind == QTY_SURFACE_TYPE)         .or. &
   (obs_kind == QTY_CLOUD_FRACTION)   .or. &
   (obs_kind == QTY_10M_U_WIND_COMPONENT) .or. (obs_kind == QTY_10M_V_WIND_COMPONENT)) then
    istatus = 0
   if ( debug > 1 .and. do_output()) print *, &
   'model_interpolate: pass sfc_elev_max_diff check for ',trim(get_name_for_quantity(obs_kind))
else
! Reject obs if the station height is far way from the model terrain.
   if(is_vertical(location, "SURFACE").and. sfc_elev_max_diff >= 0) then
   if(abs(llv(3) - zGridFace(1,cellid)) > sfc_elev_max_diff) then
      istatus = 12
      if (debug > 0 .and. do_output()) then
         print*, 'model_interpolate: Skip due to abs(dz) > sfc_elev_max_diff:', &
         llv(3)-zGridFace(1,cellid), trim(get_name_for_quantity(obs_kind)),' at ',trim(locstring)
         goto 100
      endif
   endif
   endif
endif

! see if observation quantity is in the state vector.  this sets an
! error code and returns without a fatal error if answer is no.
! exceptions:  the state vector has potential temp, but we can
! compute sensible temperature from pot_temp, rho, and qv.
! also there are options for the winds because mpas has both
! winds on the cell edges (normal only) and reconstructed winds
! at the cell centers (U,V).  there are namelist options to control
! which to use if both are in the state vector.  
! As another note, mpas defines the model variable qv as water vapor 
! mixing ratio while it defines q2 as 2-meter specific humidity,
! not 2-m water vapor mixing ratio, so q2 should be specified as 
! QTY_2M_SPECIFIC_HUMIDITY in mpas_state_variables in &mpas_vars_nml.

! is this field in the state?
ivar = get_progvar_index_from_kind(obs_kind)
if (ivar > 0) then
   goodkind = .true.     ! yes

else
   goodkind = .false.    ! but check for exceptions

   ! exceptions if the kind isn't directly
   ! a field in the state vector:
   select case (obs_kind)
      case (QTY_TEMPERATURE, QTY_2M_TEMPERATURE)
         goodkind = .true.
      case (QTY_SURFACE_ELEVATION, QTY_GEOPOTENTIAL_HEIGHT)
         goodkind = .true.
      case (QTY_PRESSURE)   ! surface pressure should be in the state
         goodkind = .true.
      case (QTY_SKIN_TEMPERATURE, QTY_SURFACE_TYPE, QTY_CLOUD_FRACTION)   
         goodkind = .true.
      case (QTY_SPECIFIC_HUMIDITY, QTY_2M_SPECIFIC_HUMIDITY)
         goodkind = .true.
      case (QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT)
         ! if the reconstructed winds at the cell centers aren't there,
         ! we can use the edge normal winds, if the user allows it.
         if (get_progvar_index_from_kind(QTY_EDGE_NORMAL_SPEED) > 0 &
             .and. use_u_for_wind) goodkind = .true.
   end select
endif

if (debug > 10 .and. do_output()) &
   print *, 'model_interpolate: ivar, goodkind? ', ivar, goodkind,' ',&
             trim(get_name_for_quantity(obs_kind)),' at ',trim(locstring)

! this kind is not in the state vector and it isn't one of the exceptions
! that we know how to handle.
if (.not. goodkind) then
   istatus(:) = 88
   if (debug > 4 .and. do_output()) print *, 'model_interpolate: kind rejected', obs_kind
   goto 100
endif

! Reject obs above a user specified pressure level.
! this is expensive - only do it if users want to reject observations
! at the top of the model.  negative values mean ignore this test.

if (.not.surface_obs) then

if (highest_obs_pressure_mb > 0.0) then
   if (.not.(is_vertical(location, "SURFACE")) .and. (.not.(is_vertical(location, "UNDEFINED")))) then

   call compute_pressure_at_loc(state_handle, ens_size, location, lpres, istatus)
   where (lpres < highest_obs_pressure_mb * 100.0_r8)
      ! Exclude from assimilation the obs above a user specified level
      istatus(:) = 201
   end where

   if (debug > 10 .and. do_output()) then
      do e = 1, ens_size
         if (istatus(e) == 201) print *, 'ens ', e, ' rejected, pressure < upper limit', lpres(e), highest_obs_pressure_mb
      enddo
   endif
   if (all(istatus /= 0)) goto 100   ! if everyone has failed, we can quit

   endif
endif

if (debug > 9 .and. do_output()) then
  print *, 'high pressure test is passed, ready to interpolate kind ', obs_kind
endif

endif !(.not.surface_obs) then

! Not prepared to do W interpolation at this time
if(obs_kind == QTY_VERTICAL_VELOCITY) then
   if (debug > 0 .and. do_output()) print *, 'model_interpolate: code does not handle vertical velocity yet'
   istatus(:) = 16
   goto 100
endif

! winds
if ((obs_kind == QTY_U_WIND_COMPONENT .or. &
     obs_kind == QTY_V_WIND_COMPONENT) .and. has_edge_u .and. use_u_for_wind) then
   if (obs_kind == QTY_U_WIND_COMPONENT) then
      ! return U
      call compute_u_with_rbf(state_handle, ens_size, location, .TRUE., expected_obs, istatus)
   else
      ! return V
      call compute_u_with_rbf(state_handle, ens_size, location, .FALSE., expected_obs, istatus)
   endif
   if (debug > 9 .and. do_output()) print *, 'model_interpolate: u_with_rbf ', obs_kind, istatus(1), expected_obs(1)
   if (all(istatus /= 0)) goto 100   ! if any member has failed, we can exit

else if (obs_kind == QTY_TEMPERATURE) then
   ! need to get potential temp, pressure, qv here, but can
   ! use same weights, so push all three types into the subroutine.
   tvars(1) = get_progvar_index_from_kind(QTY_POTENTIAL_TEMPERATURE)
   tvars(2) = get_progvar_index_from_kind(QTY_DENSITY)
   tvars(3) = get_progvar_index_from_kind(QTY_VAPOR_MIXING_RATIO)

   call compute_scalar_with_barycentric(state_handle, ens_size, location, 3, tvars, values, istatus)
   if (all(istatus /= 0)) goto 100

   ! convert pot_temp, density, vapor mixing ratio into sensible temperature
   expected_obs(:) = theta_to_tk(ens_size, values(1, :), values(2, :), values(3, :), istatus) 

   if (debug > 9 .and. do_output()) &
      print *, 'model_interpolate: TEMPERATURE ', istatus(1), expected_obs(1), trim(locstring)
   if ( all(istatus /= 0 ) ) goto 100

else if (obs_kind == QTY_PRESSURE) then
   call  compute_pressure_at_loc(state_handle, ens_size, location, expected_obs, istatus)

   if (debug > 9 .and. do_output()) &
      print *, 'model_interpolate: PRESSURE ', istatus(1), expected_obs(1), trim(locstring)
   if ( all(istatus /= 0 ) ) goto 100

else if (obs_kind == QTY_GEOPOTENTIAL_HEIGHT) then
   location_tmp = location
   call convert_vert_distrib(state_handle, ens_size, location_tmp, QTY_GEOPOTENTIAL_HEIGHT, VERTISHEIGHT, istatus)

   do e = 1, ens_size
     if(istatus(e) == 0) expected_obs(e) = query_location(location_tmp(e), 'VLOC')
   enddo

   if ( all(istatus /= 0 ) ) goto 100

else if (obs_kind == QTY_VAPOR_MIXING_RATIO      .or. obs_kind == QTY_2M_SPECIFIC_HUMIDITY   .or. &
         obs_kind == QTY_CLOUDWATER_MIXING_RATIO .or. obs_kind == QTY_RAINWATER_MIXING_RATIO .or. &
         obs_kind == QTY_ICE_MIXING_RATIO        .or. obs_kind == QTY_SNOW_MIXING_RATIO      .or. &
         obs_kind == QTY_GRAUPEL_MIXING_RATIO    .or. obs_kind == QTY_CLOUD_FRACTION       ) then
   tvars(1) = get_progvar_index_from_kind(obs_kind)
   call compute_scalar_with_barycentric(state_handle, ens_size, location, 1, tvars(1), values(1,:), istatus)
   expected_obs = values(1, :)

   ! Don't accept negative hydrometeors
   do e = 1, ens_size
      if(istatus(e) == 0) then
         expected_obs(e) = max(values(1, e),0.0_r8)
         if(obs_kind == QTY_CLOUD_FRACTION) expected_obs(e) = min(values(1, e),1.0_r8)
      endif
   enddo

   if (debug > 9 .and. do_output()) then
      print *, 'model_interpolate: obs_kind,name,istatus,expected_obs,location =', &
      obs_kind,trim(get_name_for_quantity(obs_kind)),istatus(1),expected_obs(1),trim(locstring)
   endif
   if ( all(istatus /= 0 ) ) goto 100

else if (obs_kind == QTY_SPECIFIC_HUMIDITY) then
   tvars(1) = get_progvar_index_from_kind(QTY_VAPOR_MIXING_RATIO)
   call compute_scalar_with_barycentric(state_handle, ens_size, location, 1, tvars(1), values(1,:), istatus)
   expected_obs = values(1, :)

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
   if ( all(istatus /= 0 ) ) goto 100

   if (debug > 9 .and. do_output()) &
      print *, 'model_interpolate: ', trim(get_name_for_quantity(obs_kind)), istatus(1), &
                                      expected_obs(1), trim(locstring)

! Only for the variables NOT included in the dart state vector.
! Anything included in the state_handle should not go into the case right down here.
! ex. QTY_SURFACE_PRESSURE will be taken care of in the else statement (one further down)
! as a generic interpolation case thru compute_scalar_with_barycentric
! because that obs_kind is included in mpas_state_variables in &mpas_vars_nml (e.g. state vector).
else if (obs_kind == QTY_SURFACE_ELEVATION .or. &
         obs_kind == QTY_SKIN_TEMPERATURE  .or. obs_kind == QTY_SURFACE_TYPE ) then

   if (obs_kind == QTY_SURFACE_ELEVATION) then
       call compute_surface_data_with_barycentric(zGridFace(1,:), location, expected_obs(1), istatus(1))
   else if (obs_kind == QTY_SKIN_TEMPERATURE) then
       call compute_surface_data_with_barycentric(skintemp(:), location, expected_obs(1), istatus(1))
   else if (obs_kind == QTY_SURFACE_TYPE) then
       call get_surftype(nCells,surftype)
       call compute_surface_data_with_barycentric(surftype(:)*1.0_r8, location, expected_obs(1), istatus(1))
   endif

   expected_obs(2:ens_size) = expected_obs(1)
   istatus(2:ens_size) = istatus(1)

   if (debug > 9 .and. do_output()) &
      print *, 'model_interpolate: ', trim(get_name_for_quantity(obs_kind)), ' ', istatus(1), &
                                      expected_obs(1), trim(locstring)
   if ( all(istatus /= 0 ) ) goto 100

else
   ! all other kinds come here.
   ! direct interpolation: kind is in the state vector and no clamping or other conversions needed

   tvars(1) = ivar
   call compute_scalar_with_barycentric(state_handle, ens_size, location, 1, tvars(1), values(1,:), istatus)
   expected_obs = values(1, :)

   if (debug > 9 .and. do_output()) &
       print*, 'generic interpolation in compute_scalar_with_barycentric: ', obs_kind, &
               trim(get_name_for_quantity(obs_kind))
   if (debug > 11 .and. do_output()) then
       do e = 1, ens_size
          print*, 'member ',e, ' istatus(e)=',istatus(e), ', expected_obs(e)=', expected_obs(e)
       enddo
   endif
   if ( all(istatus /= 0 ) ) goto 100

endif


100 continue

! this is for debugging - when we're confident the code is returning
! consistent values and rc codes, both these tests can be removed for speed. 
! also optionally check for generated NaNs for now.

do e = 1, ens_size

   if ((istatus(e) < 0) .or. &
       (istatus(e) /= 0 .and. expected_obs(e) /= MISSING_R8) .or. &
       (istatus(e) == 0 .and. expected_obs(e) == MISSING_R8)) then

      write(string2,*) 'member ',e,' obs_kind', obs_kind,' value = ', expected_obs(e), &
                       ' istatus = ', istatus(e), ' cellid:',cellid
      write(string3,*) trim(locstring)

      if (istatus(e) < 0) then
         write(string1,*) 'interp routine returned a negative status which is an illegal value'
      else if (istatus(e) /= 0 .and. expected_obs(e) /= MISSING_R8) then
         write(string1,*) 'interp routine returned a bad status but not a MISSING_R8 value'
         expected_obs(e) = MISSING_R8  
      else
         write(string1,*) 'interp routine returned a good status but set value to MISSING_R8'
      endif

      call error_handler(E_ALLMSG,'model_interpolate', string1, source,revision,revdate, &
                         text2=string2, text3=string3)
   endif

   ! the only portable, reliable test for NaNs we know - if some number is neither
   ! less than nor equal to/greater than 0, it must be a NaN.  all numerical comparisons
   ! fail if one or more of the operands are NaN.

   if (.not. expected_obs(e) < 0 .and. .not. expected_obs(e) >= 0) then
      write(string1,*) 'Skip member ', e, ' expected obs may be NaN: ', expected_obs(e), &
      ' obs_kind', obs_kind,' istatus = ', istatus(e)
      write(string2,*) 'cellid:',cellid,' location ', trim(locstring)
      expected_obs(e) = MISSING_R8
      istatus(e) = 99
      call error_handler(E_ERR,'model_interpolate', string1, source,revision,revdate,text2=string2)
   endif

   if (debug > 11 .and. do_output()) then
      write(string2,*) 'Completed for member ',e,' obs_kind', obs_kind,' expected_obs = ', expected_obs(e)
      write(string3,*) 'istatus = ', istatus(e), ' at ', trim(locstring)
      call error_handler(E_ALLMSG, 'model_interpolate', string2, source, revision, revdate, text2=string3)
   endif
enddo

end subroutine model_interpolate


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
call nc_add_global_attribute(ncid, "model_revision", revision)
call nc_add_global_attribute(ncid, "model_revdate", revdate)

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

! Returns the size of the model as an integer.
! Required for all applications.

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

!print *, 'model_size = ', model_size
get_model_size = model_size

end function get_model_size

!------------------------------------------------------------------
!> Returns the number of variables as an integer.

function get_num_vars()

integer :: get_num_vars

if ( .not. module_initialized ) call static_init_model

get_num_vars = nfields   

end function get_num_vars

!------------------------------------------------------------------
!> Returns the the time step of the model; the smallest increment
!> in time that the model is capable of advancing the state in a given
!> implementation. This interface is required for all applications.

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

!>@todo If MPAS ever supports more than a single domain then
!>look at the wrf model_mod code for how to change this.  you
!>have to separate out the total number of variables across
!>all domains for the min/max part, and then loop over only
!>the number of variables in each domain in the second part.

num_variables = get_num_variables(anl_domid)

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
      call convert_vert_distrib(state_handle, 1, location_arr, base_qty, vert_localization_coord, istat_arr)
      istatus1 = istat_arr(1)
      base_loc = location_arr(1)
      if(debug > 4 .and. do_output()) then
         call write_location(0,base_loc,charstring=string1)
         call error_handler(E_ALLMSG, 'get_close_obs: base_loc',string1,source, revision, revdate)
     endif
   endif
endif

if (istatus1 == 0) then

   ! Loop over potentially close subset of obs priors or state variables
   ! This way, we are decreasing the number of distance computations that will follow.
   ! This is a horizontal-distance operation and we don't need to have the relevant vertical
   ! coordinate information yet (for locs).
   call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
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
              call convert_vert_distrib(state_handle, 1, location_arr, loc_qtys(t_ind), vert_localization_coord, istat_arr)
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

if ((debug > 2) .and. do_output()) then
   call write_location(0,base_loc,charstring=string2)
   print *, 'get_close_obs: nclose, base_loc ', num_close, trim(string2)
endif

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
      call convert_vert_distrib(state_handle, 1, location_arr, base_qty, vert_localization_coord, istat_arr)
      istatus1 = istat_arr(1)
      base_loc = location_arr(1)
      if(debug > 9 .and. do_output()) then
         print*, 'get_close_state: istatus1 from convert_vert_distrib for base_loc, base_which = ',&
                  istatus1,base_which
         call write_location(0,base_loc,charstring=string1)
         call error_handler(E_ALLMSG, 'get_close_state: base_loc',string1,source, revision, revdate)
      endif
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
              call convert_vert_distrib_state(state_handle, 1, location_arr, loc_qtys(t_ind), &
                                              loc_indx(t_ind), vert_localization_coord, istat_arr)
              istatus2 = istat_arr(1)
              locs(t_ind) = location_arr(1)
              if(istatus2 /= 0 .and. debug > 9 .and. do_output()) then
                 call write_location(0,local_obs_loc,charstring=string1)
                 call error_handler(E_ALLMSG, 'get_close_state: local_obs_loc',string1,source, revision, revdate)
              endif
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

if ((debug > 2) .and. do_output()) then
   call write_location(0,base_loc,charstring=string2)
   print *, 'get_close_state: nclose, base_loc ', num_close, trim(string2)
endif

end subroutine get_close_state


!==================================================================
! The (model-specific) additional public interfaces come next
!  (these are not required by dart but are used by other programs)
!==================================================================

subroutine get_init_template_filename( filename )

! return the name of the template filename that was set
! in the model_nml namelist (ex. init.nc)

character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(init_template_filename)

end subroutine get_init_template_filename


!-------------------------------------------------------------------
! modify what static_init_model does.  this *must* be called before
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


!-------------------------------------------------------------------

!>@todo FIXME: 
!>  i believe no one is calling this routine anymore.  we should remove
!>  it from the public list and see what breaks.  if nothing, remove it.
!>
!> these need to be replaced by calls to the state structure.
!> at the moment they are only used by the postprocessing program that
!> moves the wind increments from cell centers to edges - it's not needed
!> by anything in filter.  it's a holdover from code before the state structure
!> was a general facility.

subroutine analysis_file_to_statevector(filename, state_vector, model_time)

! Reads the current time and state variables from a mpas analysis
! file and packs them into a dart state vector.

character(len=*), intent(in)    :: filename
real(r8),         intent(inout) :: state_vector(:)
type(time_type),  intent(out)   :: model_time

! temp space to hold data while we are reading it
integer  :: i, ivar
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array
real(r8), allocatable, dimension(:,:,:)     :: data_3d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME) :: varname
integer :: VarID, ncNdims, dimlen
integer :: ncid, TimeDimID, TimeDimLength
character(len=256) :: myerrorstring
character(len=*), parameter :: routine = 'analysis_file_to_statevector'

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

! Check that the input file exists ...

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'analysis_file_to_statevector',string1,source,revision,revdate)
endif

ncid = nc_open_file_readonly(filename, routine)

model_time = get_analysis_time(ncid, filename)

! Start counting and filling the state vector one item at a time,
! repacking the Nd arrays into a single 1d list of numbers.

! The DART prognostic variables are only defined for a single time.
! We already checked the assumption that variables are xy2d or xyz3d ...
! IF the netCDF variable has a TIME dimension, it must be the last dimension,
! and we need to read the LAST timestep and effectively squeeze out the
! singleton dimension when we stuff it into the DART state vector.

TimeDimID = FindTimeDimension( ncid )

if ( TimeDimID > 0 ) then
   call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=TimeDimLength), &
            'analysis_file_to_statevector', 'inquire timedimlength '//trim(filename))
else
   TimeDimLength = 0
endif

do ivar=1, nfields

   varname = progvar(ivar)%varname
   myerrorstring = trim(filename)//' '//trim(varname)

   ! determine the shape of the netCDF variable

   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            'analysis_file_to_statevector', 'inq_varid '//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
            'analysis_file_to_statevector', 'inquire '//trim(myerrorstring))

   mystart = 1   ! These are arrays, actually.
   mycount = 1

   ! Only checking the shape of the variable - sans TIME
   DimCheck : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'analysis_file_to_statevector', string1)

      mycount(i) = dimlen

   enddo DimCheck

   where(dimIDs == TimeDimID) mystart = TimeDimLength  ! pick the latest time
   where(dimIDs == TimeDimID) mycount = 1              ! only use one time

   if ((debug > 0) .and. do_output()) then
      write(*,*)'analysis_file_to_statevector '//trim(varname)//' start = ',mystart(1:ncNdims)
      write(*,*)'analysis_file_to_statevector '//trim(varname)//' count = ',mycount(1:ncNdims)
   endif

   if (ncNdims == 1) then

      ! If the single dimension is TIME, we only need a scalar.
      ! Pretty sure this cannot happen ...
      allocate(data_1d_array(mycount(1)))
      call nc_check(nf90_get_var(ncid, VarID, data_1d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'analysis_file_to_statevector', 'get_var '//trim(varname))

      call prog_var_to_vector(data_1d_array, state_vector, ivar)
      deallocate(data_1d_array)

   elseif (ncNdims == 2) then

      allocate(data_2d_array(mycount(1), mycount(2)))
      call nc_check(nf90_get_var(ncid, VarID, data_2d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'analysis_file_to_statevector', 'get_var '//trim(varname))

      call prog_var_to_vector(data_2d_array, state_vector, ivar)
      deallocate(data_2d_array)

   elseif (ncNdims == 3) then

      allocate(data_3d_array(mycount(1), mycount(2), mycount(3)))
      call nc_check(nf90_get_var(ncid, VarID, data_3d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'analysis_file_to_statevector', 'get_var '//trim(varname))

      call prog_var_to_vector(data_3d_array, state_vector, ivar)
      deallocate(data_3d_array)

   else
      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'analysis_file_to_statevector', string1, &
                        source,revision,revdate)
   endif

enddo

call nc_close_file(ncid, routine)

end subroutine analysis_file_to_statevector


!-------------------------------------------------------------------

!> @todo FIXME  this routine is a holdover from the days before
!> we had a state structure module.  it should be replaced by
!> something like write_state and other calls to the state structure
!> module.
!>
!> but there is one twist here.  this code is only being called
!> by update_mpas_states, which does a simple ncks copy from the
!> output of filter (netcdf with only the state vector fields in it)
!> to a full mpas restart file with grid info, other fields, etc.
!> 
!> however it also averages the cell center winds and projects them
!> onto the cell edges and updates 'u' - which is usually not in the
!> model state, so it isn't in the input file, only the output file.
!>
!> the normal i/o routines wouldn't write it.  they would write 
!> everything else, so read_state()/write_state() plus one additional
!> 'fix_u' routine would seem to suffice.  TBD.
!>

subroutine statevector_to_analysis_file(state_vector, ncid, filename)

! Writes the current time and state variables from a dart state
! vector (1d array) into a mpas netcdf analysis file.

real(r8),         intent(in) :: state_vector(:)
integer,          intent(in) :: ncid
character(len=*), intent(in) :: filename

! temp space to hold data while we are writing it
integer :: i, ivar
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array
real(r8), allocatable, dimension(:,:,:)     :: data_3d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME) :: varname
integer :: VarID, ncNdims, dimlen
integer :: TimeDimID, TimeDimLength
logical :: done_winds
type(time_type) :: model_time

if ( .not. module_initialized ) call static_init_model

TimeDimID = FindTimeDimension( ncid )

if ( TimeDimID > 0 ) then
   call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=TimeDimLength), &
            'statevector_to_analysis_file', 'inquire timedimlength '//trim(filename))
else
   TimeDimLength = 0
endif

done_winds = .false.
PROGVARLOOP : do ivar=1, nfields

   varname = progvar(ivar)%varname
   string2 = trim(filename)//' '//trim(varname)

   if (( varname == 'uReconstructZonal' .or. &
         varname == 'uReconstructMeridional' ) .and. update_u_from_reconstruct ) then
      if (done_winds) cycle PROGVARLOOP

      ! this routine updates the edge winds from both the zonal and meridional
      ! fields, so only call it once.
      call update_wind_components(ncid, filename, state_vector, use_increments_for_u_update)
      done_winds = .true.
      cycle PROGVARLOOP
   endif
   if ( varname == 'u' .and. update_u_from_reconstruct ) then
      write(string1, *) 'skipping update of edge normal winds (u) because'
      write(string2, *) 'update_u_from_reconstruct is True'
      call error_handler(E_MSG,'statevector_to_analysis_file',string1,&
                         source,revision,revdate, text2=string2)
      cycle PROGVARLOOP
   endif

   ! Ensure netCDF variable is conformable with progvar quantity.
   ! The TIME and Copy dimensions are intentionally not queried
   ! by looping over the dimensions stored in the progvar type.

   call nc_check(nf90_inq_varid(ncid, varname, VarID), &
            'statevector_to_analysis_file', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
            'statevector_to_analysis_file', 'inquire '//trim(string2))

   mystart = 1   ! These are arrays, actually.
   mycount = 1
   DimCheck : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'statevector_to_analysis_file', string1)

      mycount(i) = dimlen

   enddo DimCheck


   where(dimIDs == TimeDimID) mystart = TimeDimLength
   where(dimIDs == TimeDimID) mycount = 1   ! only the latest one

   if ((debug > 0) .and. do_output()) then
      write(*,*)'statevector_to_analysis_file '//trim(varname)//' start is ',mystart(1:ncNdims)
      write(*,*)'statevector_to_analysis_file '//trim(varname)//' count is ',mycount(1:ncNdims)
   endif

!> @todo FIXME the clamping can be done on the 1d array by getting
!> the start/end index from the state vector.  you can also call nf90_put_var() 
!> with a 1d conformable array without allocating, copying, writing, and freeing
!> extra space.   or call write_state() which does all this for us.

   if (progvar(ivar)%numdims == 1) then
      allocate(data_1d_array(mycount(1)))
      call vector_to_prog_var(state_vector, ivar, data_1d_array)

      ! did the user specify lower and/or upper bounds for this variable?
      ! if so, follow the instructions to either fail on out-of-range values,
      ! or set out-of-range values to the given min or max vals
      if ( progvar(ivar)%clamping ) then
         call do_clamping(progvar(ivar)%out_of_range_fail, progvar(ivar)%range, &
                          progvar(ivar)%numdims, varname, array_1d = data_1d_array)
      endif

      call nc_check(nf90_put_var(ncid, VarID, data_1d_array, &
            start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'statevector_to_analysis_file', 'put_var '//trim(varname))
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

      call nc_check(nf90_put_var(ncid, VarID, data_2d_array, &
            start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'statevector_to_analysis_file', 'put_var '//trim(varname))
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

      call nc_check(nf90_put_var(ncid, VarID, data_3d_array, &
            start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'statevector_to_analysis_file', 'put_var '//trim(varname))
      deallocate(data_3d_array)

   else
      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'statevector_to_analysis_file', string1, &
                        source,revision,revdate)
   endif

enddo PROGVARLOOP


end subroutine statevector_to_analysis_file

!------------------------------------------------------------------

!> regional mpas only:
!> update the boundary fields based on the analysis file values
!> in the boundary region, which were updated by blending the analysis
!> and the original (e.g., prior) boundary values.
!> There are options to update edge winds directly (by blending 'u'),
!> or to update reconstructed winds at cell centers first then project them onto edges.
!> When the edge wind is updated, it can be replaced by the blended value,
!> or modified with the increments (from cell-center winds).
!> state_vector - read from both ncid_a and ncid_b.
!> ncid_b = lbc to be overwritten with blended fields in the boundary zone.
!> ncid_a = analysis - either init or restart.

subroutine statevector_to_boundary_file(state_vector, ncid_b, ncid_a, &
    lbc_update_from_reconstructed_winds, lbc_update_winds_from_increments, idebug)

real(r8), intent(inout) :: state_vector(:)
integer,  intent(in)    :: ncid_b, ncid_a
logical,  intent(in)    :: lbc_update_from_reconstructed_winds
logical,  intent(in)    :: lbc_update_winds_from_increments
integer,  intent(in)    :: idebug

integer :: i, a_ivar, b_ivar, ivar, ivar_u, avar_u
integer(i8) :: a_index, b_index, l, sb_index, eb_index
integer :: cellid, vert_level, ndims, nvars, dims(3), col
integer :: adims, dima(3), iup   ! HA
integer :: edgeid
real(r8) :: weight
real(r8), allocatable :: lbc_u(:,:), lbc_ucell(:,:), lbc_vcell(:,:)
real(r8), allocatable :: delta_u(:,:), old_lbc_ucell(:,:), old_lbc_vcell(:,:)
real(r8), allocatable :: inc_lbc_ucell(:,:), inc_lbc_vcell(:,:)

character(len=NF90_MAX_NAME) :: avarname, bvarname
character(len=*), parameter :: routine = 'statevector_to_boundary_file'

nvars = get_num_variables(lbc_domid)

write(string1, *) 'lbc_update_from_reconstructed_winds = ', lbc_update_from_reconstructed_winds
write(string2, *) 'lbc_update_winds_from_increments = ',lbc_update_winds_from_increments
call error_handler(E_MSG,'statevector_to_boundary_file',string1,&
                         source,revision,revdate, text2=string2)

! save a copy of the reconstructed cell winds in separate arrays
! if we are doing an incremental update of the edge normal winds.
if (lbc_update_winds_from_increments) then
   if (.not. lbc_file_has_reconstructed_winds) then
      write(string1, *) 'Cannot update edge winds from increments because the boundary file does not contain the reconstructed winds (lbc_ur, lbc_vr)'
      write(string2, *) 'lbc_update_winds_from_increments should be .false.'
      call error_handler(E_MSG,'statevector_to_boundary_file',string1,&
                         source,revision,revdate, text2=string2)
   endif

   allocate(old_lbc_ucell(nVertLevels, nCells))
   allocate(old_lbc_vcell(nVertLevels, nCells))
   allocate(      delta_u(nVertLevels, nEdges))

   ivar = get_varid_from_varname(lbc_domid, 'lbc_ur')
   call bdy_vector_to_prog_var(state_vector, ivar, old_lbc_ucell)

   ivar = get_varid_from_varname(lbc_domid, 'lbc_vr')
   call bdy_vector_to_prog_var(state_vector, ivar, old_lbc_vcell)
endif

! for each cell in the grid, find the analysis in the
! boundary region and blend them with prior lbc values.
CELLS: do cellid = 1, nCells

   ! Soyoung: We blend the analysis in the boundary zone only.
   if (.not. on_boundary_cell(cellid)) cycle CELLS

   ! 1.0 is interior, 0.0 is exterior boundary
   weight = get_analysis_weight(cellid)

   ! do all variables associated with this cellid.
   
   VARLOOP: do b_ivar = 1, nvars

      bvarname = get_variable_name(lbc_domid, b_ivar)
      if (bvarname(1:4) /= 'lbc_') then
         write(string1, *) 'skipping update of boundary variable ', trim(bvarname)
         write(string2, *) 'because the name does not start with "lbc"'
         call error_handler(E_MSG,'statevector_to_boundary_file',string1,&
                            source,revision,revdate, text2=string2)
         cycle VARLOOP
      endif

      ! skip edge normal 'U' winds here - they will
      ! be handled in a separate code section below.
      if (bvarname == 'lbc_u') cycle VARLOOP

      ! get corresponding field in analysis domain
      avarname = trim(bvarname(5:))

      ! reconstructed cell-center winds have different names in the lbc file.
      if (bvarname == 'lbc_ur') avarname = 'uReconstructZonal'
      if (bvarname == 'lbc_vr') avarname = 'uReconstructMeridional'

      a_ivar = get_varid_from_varname(anl_domid, avarname)

      !if(do_output().and.idebug > 4) print*, 'statevector_to_boundary_file: ', &
      !                trim(bvarname), b_ivar, trim(avarname), a_ivar

      call find_mpas_dims(lbc_domid, b_ivar, ndims, dims)

      ! HA: double-check if dimensions are the same between lbc_domid and anl_domid.
      call find_mpas_dims(anl_domid, a_ivar, adims, dima)
      if(dims(1) /= dima(1) .or. dims(2) /= dima(2)) then
         write(string1, *) 'Dimension mismatches:',dims,' vs.',dima
         call error_handler(E_ERR,'statevector_to_boundary_file',string1,&
                            source,revision,revdate)
         exit
      endif

      ! loop over vert_levels.
      THISCOL: do col=1, dims(1)
         a_index = get_dart_vector_index(col, cellid, 1, anl_domid, a_ivar)
         b_index = get_dart_vector_index(col, cellid, 1, lbc_domid, b_ivar)
   
         ! compute (1-w)*x_lbc + w*x_anl
         state_vector(a_index) = (1.0_r8 - weight) * state_vector(b_index) + &
                                           weight  * state_vector(a_index)

      enddo THISCOL

   enddo VARLOOP

enddo CELLS

!> here is where we fix up the U edge normal winds
!> two options - do them directly, or compute increments from the
!> reconstructed winds and update from them.

if (.not. lbc_update_from_reconstructed_winds) then

   ! this is the prior u, not updated yet
   a_ivar = get_varid_from_varname(anl_domid, 'u')       ! analysis edge winds
   b_ivar = get_varid_from_varname(lbc_domid, 'lbc_u')   ! prior edge winds in the lbc file

   call find_mpas_dims(lbc_domid, b_ivar, ndims, dims)

   ! for each edge in the grid, find the ones which are in the
   ! boundary region and blend their values.
   EDGES: do edgeid = 1, nEdges
   
      if (.not. on_boundary_edge(edgeid)) cycle EDGES
   
      ! 1.0 is interior, 0.0 is exterior boundary
      weight = get_analysis_weight(edgeid,.false.)
   
      ! loop over vert_levels.
      THATCOL: do col=1, dims(1)
            a_index = get_dart_vector_index(col, edgeid, 1, anl_domid, a_ivar)
            b_index = get_dart_vector_index(col, edgeid, 1, lbc_domid, b_ivar)
      
            ! compute (1-w)*x_lbc + w*x_anl
            state_vector(a_index) = (1.0_r8 - weight) * state_vector(b_index) + &
                                              weight  * state_vector(a_index)
      enddo THATCOL
   
   enddo EDGES

else  ! do the increment process

   ! We only blended cell-center fields (e.g. looping over all the variables in the CELLS loop above, but not over 'u' in nEdges).
   ! Now we compute diffs (or increments) between the blended ur (vr) and the prior lbc_ur (vr).
   
   allocate(        lbc_u(nVertLevels, nEdges))
   allocate(    lbc_ucell(nVertLevels, nCells))
   allocate(    lbc_vcell(nVertLevels, nCells))
   
   ! these analyses have been blended in the boundary zone already.
   ivar = get_varid_from_varname(anl_domid, 'uReconstructZonal')
   call vector_to_prog_var(state_vector, ivar, lbc_ucell)
   
   ivar = get_varid_from_varname(anl_domid, 'uReconstructMeridional')
   call vector_to_prog_var(state_vector, ivar, lbc_vcell)
   
   ! this is the analysis u, not blended in the boundary zone yet.
   avar_u = get_varid_from_varname(anl_domid, 'u')
   call vector_to_prog_var(state_vector, avar_u, lbc_u)

   if (idebug > 9 .and. do_output()) print *, 'MIN/MAX lbc_u before update:',MINVAL(lbc_u),MAXVAL(lbc_u)
   
   if (lbc_update_winds_from_increments) then
   
      ! project analysis increments at cell centers onto the edges.
   
      allocate(inc_lbc_ucell(nVertLevels, nCells))
      allocate(inc_lbc_vcell(nVertLevels, nCells))
   
      inc_lbc_ucell = lbc_ucell - old_lbc_ucell 
      inc_lbc_vcell = lbc_vcell - old_lbc_vcell 
   
      call uv_cell_to_edges(inc_lbc_ucell, inc_lbc_vcell, delta_u)
   
      ! Soyoung: Add the blended u increments back to lbc_u in the boundary zone.
      !          We should not change the analysis u in the interior domain, but
      !          We should also check bdyMaskCell for the two adjacent cells as
      !          bdyMaskEdge is assigned with the lower mask value between the two
      !          cells.
      ! Ex) An edge between cell1 (w/ bdyMaskCell = 0) and cell2 (w/ bdyMaskCell = 1)
      ! has bdyMaskEdge = 0. In this case, even if bdyMaskEdge of the edge is zero,
      ! cell2 has been updated in the CELLS loop above, so the edge has to be updated.

      iup = 0
      IEDGE: do edgeid = 1, nEdges

      if (.not. on_boundary_edge(edgeid) .and. &
          .not. on_boundary_cell(cellsOnEdge(1,edgeid)) .and. &
          .not. on_boundary_cell(cellsOnEdge(2,edgeid)) ) cycle IEDGE
   
           lbc_u(:,edgeid) = lbc_u(:,edgeid) + delta_u(:,edgeid)

           if (idebug > 9 .and. .not.on_boundary_edge(edgeid)) &
               print*, iup, edgeid, delta_u(1,edgeid)

      enddo IEDGE

      if (idebug > 9 .and. do_output()) print*, 'MIN/MAX delta_u:            ',MINVAL(delta_u),MAXVAL(delta_u)

      deallocate(old_lbc_ucell, old_lbc_vcell, delta_u)

   else
   
      ! just replace, no increments
      ! project the diffs (or increments) onto the edges by calling uv_cell_to_edges.
      call uv_cell_to_edges(lbc_ucell, lbc_vcell, lbc_u, .true.)
      
   endif
   if (idebug > 9 .and. do_output()) print *, 'MIN/MAX lbc_u  after update:',MINVAL(lbc_u),MAXVAL(lbc_u)
   
   ! put lbc_u array data back into the state_vector
   
   sb_index = get_index_start(anl_domid, avar_u)
   eb_index = get_index_end  (anl_domid, avar_u)
   state_vector(sb_index:eb_index) = reshape(lbc_u, (/eb_index-sb_index+1/) )
      
   deallocate(lbc_u, lbc_ucell, lbc_vcell)

endif   ! U updates
   

! for each boundary variable, write it to the output file.
VARLOOP2: do b_ivar = 1, nvars

   bvarname = get_variable_name(lbc_domid, b_ivar)
   if (bvarname(1:4) /= 'lbc_') then
      write(string1, *) 'skipping update of boundary variable ', trim(bvarname)
      write(string2, *) 'because the name does not start with "lbc"'
      call error_handler(E_MSG,'statevector_to_boundary_file',string1,&
                         source,revision,revdate, text2=string2)
      cycle VARLOOP2
   endif

   ! by default strip off 'lbc_' from boundary file name
   ! and update that field name for the analysis file,
   ! unless it doesn't follow the pattern
   avarname = trim(bvarname(5:))

   if (bvarname == 'lbc_ur') avarname = 'uReconstructZonal'
   if (bvarname == 'lbc_vr') avarname = 'uReconstructMeridional'

   ! Soyoung - we blended the analysis vector in the boundary zone.
   a_ivar = get_varid_from_varname(anl_domid, avarname)

   ! nsc - it's possible we could remove this line and the
   !       reshape()s from put_variable lines below.
   call find_mpas_dims(lbc_domid, b_ivar, ndims, dims)

   sb_index = get_index_start(anl_domid, a_ivar)
   eb_index = get_index_end  (anl_domid, a_ivar)

   if (idebug > 9 .and. do_output()) print *, 'updating ', trim(bvarname), ' min, max:',&
       minval(state_vector(sb_index:eb_index)), maxval(state_vector(sb_index:eb_index))
       !, ' with length ', eb_index - sb_index + 1

   ! Soyoung - Now, the lbc file has the analysis in the interior and the blended values
   !           in the boundary zone. In other words, the lbc file is updated not only in
   !           the boundary, but over the entire domain (although it is used only for the
   !           boundary zone).
   call nc_put_variable(ncid_b, bvarname, reshape(state_vector(sb_index:eb_index), dims), routine)
!   call nc_put_variable(ncid_a, avarname, reshape(state_vector(sb_index:eb_index), dims), routine)

enddo VARLOOP2

end subroutine statevector_to_boundary_file

!------------------------------------------------------------------

subroutine do_clamping(out_of_range_fail, range, dimsize, varname, array_1d, array_2d, array_3d)
 logical,          intent(in)    :: out_of_range_fail
 real(r8),         intent(in)    :: range(2)
 integer,          intent(in)    :: dimsize
 character(len=*), intent(in)    :: varname
 real(r8),optional,intent(inout) :: array_1d(:), array_2d(:,:), array_3d(:,:,:)

! FIXME: This should be replaced by the clamp_variable in direct_netcdf_mpi.
! clamp_variable currently works on one dimensional arrays for variables. This
! only becomes an issue in statevector_to_analysis_file which requires the dimension
! of the variable to output to a netcdf file, however, this information is
! already stored in the state_structure.

! for a given directive and range, do the data clamping for the given
! input array.  only one of the optional array args should be specified - the
! one which matches the given dimsize.  this still has replicated sections for
! each possible dimensionality (which so far is only 1 to 3 - add 4-7 only
! if needed) but at least it is isolated to this subroutine.

! these sections should all be identical except for the array_XX specified.
! if anyone can figure out a way to defeat fortran's strong typing for arrays
! so we don't have to replicate each of these sections, i'll buy you a cookie.
! (sorry, you can't suggest using the preprocessor, which is the obvious
! solution.  up to now we have avoided any preprocessed code in the entire
! system.  if we cave at some future point this routine is a prime candidate
! to autogenerate.)

if (dimsize == 1) then
   if (.not. present(array_1d)) then
      call error_handler(E_ERR, 'do_clamping', 'Internal error.  Should not happen', &
                         source,revision,revdate, text2='array_1d not present for 1d case')
   endif

   ! is lower bound set
   if ( range(1) /= missing_r8 ) then

      if ( out_of_range_fail ) then
         if ( minval(array_1d) < range(1) ) then
            write(string1, *) 'min data val = ', minval(array_1d), &
                              'min data bounds = ', range(1)
            call error_handler(E_ERR, 'statevector_to_analysis_file', &
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
            call error_handler(E_ERR, 'statevector_to_analysis_file', &
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
      call error_handler(E_ERR, 'do_clamping', 'Internal error.  Should not happen', &
                         source,revision,revdate, text2='array_2d not present for 2d case')
   endif

   ! is lower bound set
   if ( range(1) /= missing_r8 ) then

      if ( out_of_range_fail ) then
         if ( minval(array_2d) < range(1) ) then
            write(string1, *) 'min data val = ', minval(array_2d), &
                              'min data bounds = ', range(1)
            call error_handler(E_ERR, 'statevector_to_analysis_file', &
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
            call error_handler(E_ERR, 'statevector_to_analysis_file', &
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
      call error_handler(E_ERR, 'do_clamping', 'Internal error.  Should not happen', &
                         source,revision,revdate, text2='array_3d not present for 3d case')
   endif

   ! is lower bound set
   if ( range(1) /= missing_r8 ) then

      if ( out_of_range_fail ) then
         if ( minval(array_3d) < range(1) ) then
            write(string1, *) 'min data val = ', minval(array_3d), &
                              'min data bounds = ', range(1)
            call error_handler(E_ERR, 'statevector_to_analysis_file', &
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
            call error_handler(E_ERR, 'statevector_to_analysis_file', &
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
   call error_handler(E_MSG, 'do_clamping', 'Internal error, should not happen', &
                      source,revision,revdate, text2=string1)
endif   ! dimsize

end subroutine do_clamping

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
   call error_handler(E_ERR,'get_analysis_time',string1,source,revision,revdate)
endif

call nc_check( nf90_inquire_dimension(ncid, dimIDs(1), len=idims(1)), &
                 'get_analysis_time', 'inquire time dimension length '//trim(filename))
call nc_check( nf90_inquire_dimension(ncid, dimIDs(2), len=idims(2)), &
                 'get_analysis_time', 'inquire time dimension length '//trim(filename))

if (idims(2) /= 1) then
   write(string1,*) 'multiple timesteps (',idims(2),') in file ', trim(filename)
   write(string2,*) 'We are using the LAST one, presumably, the LATEST timestep.'
   call error_handler(E_MSG,'get_analysis_time',string1,source,revision,revdate,text2=string2)
endif

! Get the highest ranking time ... the last one, basically.

call nc_check( nf90_get_var(ncid, VarID, timestring, start = (/ 1, idims(2) /)), &
              'get_analysis_time', 'get_var xtime '//trim(filename))

get_analysis_time_ncid = string_to_time(timestring)

if ((debug > 0) .and. do_output()) then
   call print_date(get_analysis_time_ncid, 'get_analysis_time:model date')
   call print_time(get_analysis_time_ncid, 'get_analysis_time:model time')
endif

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
   call error_handler(E_ERR,'get_analysis_time',string1,source,revision,revdate)
endif

! find the first digit and use that as the start of the string conversion
i = scan(filename, "0123456789")
if (i <= 0) then
   write(string1,*) 'cannot find time string in name ', trim(filename)
   call error_handler(E_ERR,'get_analysis_time',string1,source,revision,revdate)
endif

get_analysis_time_fname = string_to_time(filename(i:i+TIMELEN-1))

end function get_analysis_time_fname


!------------------------------------------------------------------

subroutine write_model_time_file(time_filename, model_time, adv_to_time)
 character(len=*), intent(in)           :: time_filename
 type(time_type),  intent(in)           :: model_time
 type(time_type),  intent(in), optional :: adv_to_time

integer :: iunit
character(len=TIMELEN) :: timestring
type(time_type)   :: deltatime

iunit = open_file(time_filename, action='write')

timestring = time_to_string(model_time)
write(iunit, '(A)') timestring

if (present(adv_to_time)) then
   timestring = time_to_string(adv_to_time)
   write(iunit, '(A)') timestring

   deltatime = adv_to_time - model_time
   timestring = time_to_string(deltatime, interval=.true.)
   write(iunit, '(A)') timestring
endif

call close_file(iunit)

end subroutine write_model_time_file

!-----------------------------------------------------------------------
subroutine write_model_time_restart(ncid, dart_time)

integer,             intent(in) :: ncid !< netcdf file handle
type(time_type),     intent(in) :: dart_time

integer :: year, month, day, hour, minute, second
character(len=64) :: timestring
character(len=*), parameter :: routine = 'write_model_time_restart'

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

end subroutine write_model_time_restart

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
                         source, revision, revdate, text2=string1)
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
!> Read the grid dimensions from the MPAS netcdf file.
!>

subroutine read_grid_dims(ncid)
integer, intent(in) :: ncid

character(len=*), parameter :: routine = 'read_grid_dims'

nCells        = nc_get_dimension_size(ncid, 'nCells',         routine)
nVertices     = nc_get_dimension_size(ncid, 'nVertices',      routine)
nEdges        = nc_get_dimension_size(ncid, 'nEdges',         routine)
maxEdges      = nc_get_dimension_size(ncid, 'maxEdges',       routine)
nVertLevels   = nc_get_dimension_size(ncid, 'nVertLevels',    routine)
nVertLevelsP1 = nc_get_dimension_size(ncid, 'nVertLevelsP1',  routine)
vertexDegree  = nc_get_dimension_size(ncid, 'vertexDegree',   routine)
nSoilLevels   = nc_get_dimension_size(ncid, 'nSoilLevels',    routine)

if (debug > 4 .and. do_output()) then
   write(*,*)
   write(*,*)'read_grid_dims: nCells        is ', nCells
   write(*,*)'read_grid_dims: nVertices     is ', nVertices
   write(*,*)'read_grid_dims: nEdges        is ', nEdges
   write(*,*)'read_grid_dims: maxEdges      is ', maxEdges
   write(*,*)'read_grid_dims: nVertLevels   is ', nVertLevels
   write(*,*)'read_grid_dims: nVertLevelsP1 is ', nVertLevelsP1
   write(*,*)'read_grid_dims: vertexDegree  is ', vertexDegree
   write(*,*)'read_grid_dims: nSoilLevels   is ', nSoilLevels
endif

end subroutine read_grid_dims


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


! A little sanity check

if ( debug > 9 .and. do_output() ) then

   write(*,*)
   write(*,*)'latCell           range ',minval(latCell),           maxval(latCell)
   write(*,*)'lonCell           range ',minval(lonCell),           maxval(lonCell)
   write(*,*)'zgrid             range ',minval(zGridFace),         maxval(zGridFace)
   write(*,*)'cellsOnVertex     range ',minval(cellsOnVertex),     maxval(cellsOnVertex)
   write(*,*)'edgeNormalVectors range ',minval(edgeNormalVectors), maxval(edgeNormalVectors)
   write(*,*)'nEdgesOnCell      range ',minval(nEdgesOnCell),      maxval(nEdgesOnCell)
   write(*,*)'EdgesOnCell       range ',minval(EdgesOnCell),       maxval(EdgesOnCell)
   write(*,*)'cellsOnEdge       range ',minval(cellsOnEdge),       maxval(cellsOnEdge)
   write(*,*)'xland             range ',minval(xland),             maxval(xland)
   write(*,*)'seaice            range ',minval(seaice),            maxval(seaice)       ! for rttov
   write(*,*)'skintemp          range ',minval(skintemp),          maxval(skintemp)     ! for rttov
   if(data_on_edges) then
      write(*,*)'latEdge        range ',minval(latEdge),           maxval(latEdge)
      write(*,*)'lonEdge        range ',minval(lonEdge),           maxval(lonEdge)
      write(*,*)'xEdge          range ',minval(xEdge),             maxval(xEdge)
      write(*,*)'yEdge          range ',minval(yEdge),             maxval(yEdge)
      write(*,*)'zEdge          range ',minval(zEdge),             maxval(zEdge)
   endif
   write(*,*)'xVertex           range ',minval(xVertex),           maxval(xVertex)
   write(*,*)'yVertex           range ',minval(yVertex),           maxval(yVertex)
   write(*,*)'zVertex           range ',minval(zVertex),           maxval(zVertex)
   write(*,*)'verticesOnCell    range ',minval(verticesOnCell),    maxval(verticesOnCell)
   if (allocated(bdyMaskCell)) &
   write(*,*)'bdyMaskCell       range ',minval(bdyMaskCell),       maxval(bdyMaskCell)
   if (allocated(bdyMaskEdge)) &
   write(*,*)'bdyMaskEdge       range ',minval(bdyMaskEdge),       maxval(bdyMaskEdge)
   if (allocated(maxLevelCell)) &
   write(*,*)'maxLevelCell      range ',minval(maxLevelCell),      maxval(maxLevelCell)

endif

end subroutine get_grid


!------------------------------------------------------------------

subroutine update_wind_components(ncid, filename, state_vector, use_increments_for_u_update)

 integer,  intent(in)  :: ncid
 character(len=*), intent(in) :: filename
 real(r8), intent(in)  :: state_vector(:)
 logical,  intent(in)  :: use_increments_for_u_update

! the winds pose a special problem because the model uses the edge-normal component
! of the winds at the center of the cell edges at half levels in the vertical ('u').
! the output files from the model can include interpolated Meridional and Zonal
! winds as prognostic fields, which are easier for us to use when computing
! forward operator values.  but in the end we need to update 'u' in the output
! model file.

! this routine is only called when 'u' is not being directly updated by the
! assimilation, and the updated cell center values need to be converted back
! to update 'u'.   there are several choices for how to do this and most are
! controlled by namelist settings.

! If 'use_increments_for_u_update' is .true.:
!  Read in the previous reconstructed winds from the original mpas netcdf file
!  and compute what increments (changes in values) were added by the assimilation.
!  Read in the original edge normal wind 'u' field from that same mpas netcdf
!  file and add the interpolated increments to compute the updated 'u' values.
!  (note that we can't use the DART Prior_Diag.nc file to get the previous
!  values if we're using Prior inflation, because the diagnostic values are
!  written out after inflation is applied.)

! If 'use_increments_for_u_update' is .false.:
!  use the Zonal/Meridional cell center values directly. The edge normal winds
!  are each directly between 2 cell centers, so average the components normal
!  to the edge direction.  don't read in the previous values at the cell centers
!  or the edge normal winds.

! there are several changes here from previous versions:
!  1. it requires both zonal and meridional fields to be there.  it doesn't
!  make sense to assimilate with only one component of the winds.
!  2. i removed the 'return if this has been called already' flag.
!  this is a generic utility routine.  if someone wrote a main program
!  that wanted to cycle over multiple files in a loop, this would have
!  only updated the first file and silently returned for all the rest.
!  3. i moved the netcdf code for 3 identical operations into a subroutine.
!  the code is easier to read and it's less likely to make a mistake if
!  the code needs to be changed (one place vs three).

! space to hold existing data from the analysis file
real(r8), allocatable :: u(:,:)              ! u(nVertLevels, nEdges)
real(r8), allocatable :: ucell(:,:)          ! uReconstructZonal(nVertLevels, nCells)
real(r8), allocatable :: vcell(:,:)          ! uReconstructMeridional(nVertLevels, nCells)
real(r8), allocatable :: data_2d_array(:,:)  ! temporary
logical :: both
integer :: zonal, meridional

if ( .not. module_initialized ) call static_init_model

! get the ivar values for the zonal and meridional wind fields
call winds_present(zonal,meridional,both)
if (.not. both) call error_handler(E_ERR, 'update_wind_components', &
   'internal error: wind fields not found', source, revision, revdate)

allocate(    u(nVertLevels, nEdges))
allocate(ucell(nVertLevels, nCells))
allocate(vcell(nVertLevels, nCells))

! if doing increments, read in 'u' (edge normal winds), plus the uReconstructZonal
! and uReconstructMeridional fields from the mpas analysis netcdf file.

if (use_increments_for_u_update) then
   call read_2d_from_nc_file(ncid, filename, 'u', u)
   call read_2d_from_nc_file(ncid, filename, 'uReconstructZonal', ucell)
   call read_2d_from_nc_file(ncid, filename, 'uReconstructMeridional', vcell)

   if ((debug > 8) .and. do_output()) then
      write(*,*)
      write(*,*)'update_winds: org u          range ',minval(u),     maxval(u)
      write(*,*)'update_winds: org zonal      range ',minval(ucell), maxval(ucell)
      write(*,*)'update_winds: org meridional range ',minval(vcell), maxval(vcell)
   endif

   ! compute the increments compared to the updated values in the state vector

   allocate(data_2d_array(nVertLevels, nCells))

   ! write updated reconstructed winds to restart file for later use.
   call vector_to_prog_var(state_vector, zonal, data_2d_array)
   call put_u(ncid, filename, data_2d_array, 'uReconstructZonal')
   ucell = data_2d_array - ucell

   call vector_to_prog_var(state_vector, meridional, data_2d_array)
   call put_u(ncid, filename, data_2d_array, 'uReconstructMeridional')
   vcell = data_2d_array - vcell

   deallocate(data_2d_array)

   ! this is by nedges, not ncells as above
   allocate(data_2d_array(nVertLevels, nEdges))
   call uv_cell_to_edges(ucell, vcell, data_2d_array)
   u(:,:) = u(:,:) + data_2d_array(:,:)
   deallocate(data_2d_array)

   if ((debug > 8) .and. do_output()) then
      write(*,*)
      write(*,*)'update_winds: u increment    range ',minval(ucell), maxval(ucell)
      write(*,*)'update_winds: v increment    range ',minval(vcell), maxval(vcell)
   endif

else

   ! The state vector has updated zonal and meridional wind components.
   ! put them directly into the arrays.  these are the full values, not
   ! just increments.
   call vector_to_prog_var(state_vector, zonal, ucell)
   call vector_to_prog_var(state_vector, meridional, vcell)

   call uv_cell_to_edges(ucell, vcell, u, .true.)

   if ((debug > 8) .and. do_output()) then
      write(*,*)
      write(*,*)'update_winds: u values    range ',minval(ucell), maxval(ucell)
      write(*,*)'update_winds: v values    range ',minval(vcell), maxval(vcell)
   endif

   call put_u(ncid, filename, ucell, 'uReconstructZonal')
   call put_u(ncid, filename, vcell, 'uReconstructMeridional')
endif

if ((debug > 8) .and. do_output()) then
   write(*,*)
   write(*,*)'update_winds: u after update:',minval(u), maxval(u)
endif

! Write back to the mpas analysis file.

call put_u(ncid, filename, u, 'u')

deallocate(ucell, vcell, u)

end subroutine update_wind_components


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


! A little sanity check

if ((debug > 8) .and. do_output()) then

   write(*,*)
   write(*,*) trim(vchar),'       range ',minval(u),     maxval(u)

endif

end subroutine put_u


!------------------------------------------------------------------

subroutine vector_to_1d_prog_var(x, ivar, data_1d_array)

! convert the values from a 1d array, starting at an offset,
! into a 1d array.

real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:),   intent(out) :: data_1d_array

integer :: start_offset, end_offset

start_offset = progvar(ivar)%index1
end_offset   = start_offset + size(data_1d_array) - 1

data_1d_array = x(start_offset:end_offset)

end subroutine vector_to_1d_prog_var


!------------------------------------------------------------------

subroutine vector_to_2d_prog_var(x, ivar, data_2d_array)

! convert the values from a 1d array, starting at an offset,
! into a 2d array.

real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:,:), intent(out) :: data_2d_array

integer :: start_offset, end_offset

start_offset = get_index_start(anl_domid, ivar)
end_offset   = get_index_end(anl_domid, ivar)

data_2d_array = reshape(x(start_offset:end_offset), shape(data_2d_array))

end subroutine vector_to_2d_prog_var


!------------------------------------------------------------------

subroutine bdy_vector_to_prog_var(x, ivar, data_2d_array)

! convert the values from a 1d array, starting at an offset,
! into a 2d array.

real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:,:), intent(out) :: data_2d_array

integer :: start_offset, end_offset

start_offset = get_index_start(lbc_domid, ivar)
end_offset   = get_index_end(lbc_domid, ivar)

data_2d_array = reshape(x(start_offset:end_offset), shape(data_2d_array))

end subroutine bdy_vector_to_prog_var


!------------------------------------------------------------------

subroutine vector_to_3d_prog_var(x, ivar, data_3d_array)

! convert the values from a 1d array, starting at an offset,
! into a 3d array.

real(r8), dimension(:),     intent(in)  :: x
integer,                    intent(in)  :: ivar
real(r8), dimension(:,:,:), intent(out) :: data_3d_array

integer :: start_offset, end_offset

start_offset = progvar(ivar)%index1
end_offset   = start_offset + size(data_3d_array) - 1

data_3d_array = reshape(x(start_offset:end_offset), shape(data_3d_array))

end subroutine vector_to_3d_prog_var


!------------------------------------------------------------------

subroutine prog_var_1d_to_vector(data_1d_array, x, ivar)

! convert the values from a 1d array into a 1d array
! starting at an offset.

real(r8), dimension(:),   intent(in)    :: data_1d_array
real(r8), dimension(:),   intent(inout) :: x
integer,                  intent(in)    :: ivar

integer :: start_offset, end_offset

start_offset = progvar(ivar)%index1
end_offset   = start_offset + size(data_1d_array) - 1

x(start_offset:end_offset) = data_1d_array

end subroutine prog_var_1d_to_vector


!------------------------------------------------------------------

subroutine prog_var_2d_to_vector(data_2d_array, x, ivar)

! convert the values from a 2d array into a 1d array
! starting at an offset.

real(r8), dimension(:,:), intent(in)    :: data_2d_array
real(r8), dimension(:),   intent(inout) :: x
integer,                  intent(in)    :: ivar

integer :: start_offset, end_offset

start_offset = progvar(ivar)%index1
end_offset   = start_offset + size(data_2d_array) - 1

x(start_offset:end_offset) = reshape(data_2d_array, (/ size(data_2d_array) /) )

end subroutine prog_var_2d_to_vector


!------------------------------------------------------------------

subroutine prog_var_3d_to_vector(data_3d_array, x, ivar)

! convert the values from a 2d array into a 1d array
! starting at an offset.

real(r8), dimension(:,:,:), intent(in)    :: data_3d_array
real(r8), dimension(:),     intent(inout) :: x
integer,                    intent(in)    :: ivar

integer :: start_offset, end_offset

start_offset = progvar(ivar)%index1
end_offset   = start_offset + size(data_3d_array) - 1

x(start_offset:end_offset) = reshape(data_3d_array, (/ size(data_3d_array) /))

end subroutine prog_var_3d_to_vector


!------------------------------------------------------------------

!>@todo fill in the inputs we need for the add_domain() routine

subroutine verify_state_variables( state_variables, ncid, filename, ngood, table )

character(len=*), dimension(:),   intent(in)  :: state_variables
integer,                          intent(in)  :: ncid
character(len=*),                 intent(in)  :: filename
integer,                          intent(out) :: ngood
character(len=*), dimension(:,:), intent(out) :: table

integer :: nrows, ncols, i, j, VarID
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: varname, dimname
character(len=NF90_MAX_NAME) :: dartstr
integer :: dimlen, numdims
logical :: u_already_in_list
logical :: failure

if ( .not. module_initialized ) call static_init_model

failure = .FALSE. ! perhaps all with go well
u_already_in_list = .FALSE.

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
      string1 = 'mpas_vars_nml:model state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
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
            call error_handler(E_MSG,'verify_state_variables',string2,source,revision,revdate)
            failure = .TRUE.
      end select

   enddo DimensionLoop

   if (failure) then
       string2 = 'unsupported dimension(s) are fatal'
       call error_handler(E_ERR,'verify_state_variables',string2,source,revision,revdate)
   endif

   ! Make sure DART kind is valid

   if( get_index_for_quantity(dartstr) < 0 ) then
      write(string1,'(''there is no qty <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector

   if ((debug > 0) .and. do_output()) then
      write(logfileunit,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
      write(     *     ,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
   endif

   ! Keep track of whether U (edge normal winds) is part of the user-specified field list
   if ((table(i, 1) == 'u') .or. (table(i, 1) == 'U')) u_already_in_list = .true.

   ngood = ngood + 1
enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'verify_state_variables',string1,source,revision,revdate,text2=string2)
endif

! if this flag is true and the user hasn't said U should be in the state,
! add it to the list.
if (add_u_to_state_list .and. .not. u_already_in_list) then
   ngood = ngood + 1
   table(ngood,1) = "u"
   table(ngood,2) = "QTY_EDGE_NORMAL_SPEED"
endif

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
  character(len=*), intent(in) :: state_list(max_state_variables, num_state_table_columns)
  character(len=*), intent(in) :: bdy_list(max_state_variables)
  logical :: is_edgedata_in_state_vector

integer :: i

StateLoop : do i = 1, max_state_variables

   if (state_list(i, 1) == ' ') exit StateLoop ! Found end of list.

   if (state_list(i, 1) == 'u') then
      is_edgedata_in_state_vector = .true.
      return
   endif

enddo StateLoop

! if U is not in the state, does it matter if U is
! in the boundary file?  yes, return true if so.
BdyLoop : do i = 1, max_state_variables

   if (bdy_list(i) == ' ') exit BdyLoop ! Found end of list.

   if (bdy_list(i) == 'lbc_u') then
      is_edgedata_in_state_vector = .true.
      return
   endif

enddo BdyLoop

is_edgedata_in_state_vector = .false.

end function is_edgedata_in_state_vector


!------------------------------------------------------------------

subroutine dump_progvar(ivar)

! dump the contents of the metadata for an individual entry.
! expected to be called in a loop or called for entries of interest.

integer,  intent(in)           :: ivar

!%! type progvartype
!%!    private
!%!    character(len=NF90_MAX_NAME) :: varname
!%!    integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
!%!    integer :: xtype         ! netCDF variable type (NF90_double, etc.)
!%!    integer :: numdims       ! number of dims - excluding TIME
!%!    integer :: numvertical   ! number of vertical levels in variable
!%!    integer :: numcells      ! number of horizontal locations (typically cell centers)
!%!    integer :: numedges
!%!    logical :: ZonHalf       ! vertical coordinate has dimension nVertLevels
!%!    integer :: varsize       ! prod(dimlens(1:numdims))
!%!    integer :: index1        ! location in dart state vector of first occurrence
!%!    integer :: indexN        ! location in dart state vector of last  occurrence
!%!    integer :: dart_kind
!%!    character(len=obstypelength) :: kind_string
!%!    logical  :: clamping     ! does variable need to be range-restricted before
!%!    real(r8) :: range(2)     ! being stuffed back into MPAS analysis file.
!%! end type progvartype

integer :: i

! take care of parallel runs where we only want a single copy of
! the output.
if (.not. do_output()) return

write(logfileunit,*)
write(     *     ,*)
write(logfileunit,*) 'variable number ',ivar,' is ',trim(progvar(ivar)%varname)
write(     *     ,*) 'variable number ',ivar,' is ',trim(progvar(ivar)%varname)
write(logfileunit,*) '  xtype       ',progvar(ivar)%xtype
write(     *     ,*) '  xtype       ',progvar(ivar)%xtype
write(logfileunit,*) '  dimlens     ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
write(     *     ,*) '  dimlens     ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
write(logfileunit,*) '  numdims     ',progvar(ivar)%numdims
write(     *     ,*) '  numdims     ',progvar(ivar)%numdims
write(logfileunit,*) '  numvertical ',progvar(ivar)%numvertical
write(     *     ,*) '  numvertical ',progvar(ivar)%numvertical
write(logfileunit,*) '  numcells    ',progvar(ivar)%numcells
write(     *     ,*) '  numcells    ',progvar(ivar)%numcells
write(logfileunit,*) '  numedges    ',progvar(ivar)%numedges
write(     *     ,*) '  numedges    ',progvar(ivar)%numedges
write(logfileunit,*) '  ZonHalf     ',progvar(ivar)%ZonHalf
write(     *     ,*) '  ZonHalf     ',progvar(ivar)%ZonHalf
write(logfileunit,*) '  varsize     ',progvar(ivar)%varsize
write(     *     ,*) '  varsize     ',progvar(ivar)%varsize
write(logfileunit,*) '  index1      ',progvar(ivar)%index1
write(     *     ,*) '  index1      ',progvar(ivar)%index1
write(logfileunit,*) '  indexN      ',progvar(ivar)%indexN
write(     *     ,*) '  indexN      ',progvar(ivar)%indexN
write(logfileunit,*) '  dart_kind   ',progvar(ivar)%dart_kind
write(     *     ,*) '  dart_kind   ',progvar(ivar)%dart_kind
write(logfileunit,*) '  kind_string ',progvar(ivar)%kind_string
write(     *     ,*) '  kind_string ',progvar(ivar)%kind_string
write(logfileunit,*) '  clamping    ',progvar(ivar)%clamping
write(     *     ,*) '  clamping    ',progvar(ivar)%clamping
write(logfileunit,*) '  clamp range  ',progvar(ivar)%range
write(     *     ,*) '  clamp range  ',progvar(ivar)%range
write(logfileunit,*) '  clamp fail   ',progvar(ivar)%out_of_range_fail
write(     *     ,*) '  clamp fail   ',progvar(ivar)%out_of_range_fail
do i = 1,progvar(ivar)%numdims
   write(logfileunit,*) '  dimension/length/name ',i,progvar(ivar)%dimlens(i),trim(progvar(ivar)%dimname(i))
   write(     *     ,*) '  dimension/length/name ',i,progvar(ivar)%dimlens(i),trim(progvar(ivar)%dimname(i))
enddo

end subroutine dump_progvar

!------------------------------------------------------------------

subroutine print_variable_ranges(x,title)

! given a state vector, print out the min and max
! data values for the variables in the vector.

real(r8), intent(in) :: x(:)
character(len=*), optional, intent(in)  :: title

integer :: ivar

if (present(title)) write(*,*) trim(title)

do ivar = 1, nfields
   call print_minmax(ivar, x)
enddo

end subroutine print_variable_ranges

!------------------------------------------------------------------

subroutine print_minmax(ivar, x)

! given an index and a state vector, print out the min and max
! data values for the items corresponding to that progvar index.

integer,  intent(in) :: ivar
real(r8), intent(in) :: x(:)

write(string1, '(A,A32,2F16.7)') 'data  min/max ', trim(progvar(ivar)%varname), &
           minval(x(progvar(ivar)%index1:progvar(ivar)%indexN)), &
           maxval(x(progvar(ivar)%index1:progvar(ivar)%indexN))

call error_handler(E_MSG, '', string1, source,revision,revdate)

end subroutine print_minmax


!------------------------------------------------------------------

function FindTimeDimension(ncid) result(timedimid)

! Find the Time Dimension ID in a netCDF file.
! If there is none - (spelled the obvious way) - the routine
! returns a negative number. You don't HAVE to have a TIME dimension.

integer                      :: timedimid
integer,          intent(in) :: ncid

integer :: nc_rc

TimeDimID = -1 ! same as the netCDF library routines.
nc_rc = nf90_inq_dimid(ncid,'Time',dimid=TimeDimID)

end function FindTimeDimension


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
call error_handler(E_ERR,'winds_present',string1,source,revision,revdate,text2=string2)

end subroutine winds_present


!------------------------------------------------------------
subroutine get_variable_bounds(bounds_table, ivar)

! matches MPAS variable name in bounds table to assign
! the bounds if they exist.  otherwise sets the bounds
! to missing_r8
!
! SYHA (May-30-2013)
! Adopted from wrf/model_mod.f90 after adding mpas_state_bounds in mpas_vars_nml.

! bounds_table is the global mpas_state_bounds
character(len=*), intent(in)  :: bounds_table(num_bounds_table_columns, max_state_variables)
integer,          intent(in)  :: ivar

! local variables
character(len=50)             :: bounds_varname, bound
character(len=10)             :: clamp_or_fail
real(r8)                      :: lower_bound, upper_bound
integer                       :: n

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
             call error_handler(E_ERR,'get_variable_bounds',string1, &
                                source,revision,revdate)
        else if ( clamp_or_fail == 'CLAMP' ) then
             progvar(ivar)%out_of_range_fail = .FALSE.
        else if ( clamp_or_fail == 'FAIL' ) then
             progvar(ivar)%out_of_range_fail = .TRUE.
        else
             write(string1, *) 'last column must be "CLAMP" or "FAIL" for ', &
                  trim(bounds_varname)
             call error_handler(E_ERR,'get_variable_bounds',string1, &
                  source,revision,revdate, text2='found '//trim(clamp_or_fail))
        endif

        ! Assign the clamping information into the variable
        progvar(ivar)%clamping = .true.
        progvar(ivar)%range    = (/ lower_bound, upper_bound /)

        if ((debug > 0) .and. do_output()) then
           write(*,*) 'In get_variable_bounds assigned ', trim(progvar(ivar)%varname)
           write(*,*) ' clamping range  ',progvar(ivar)%range
        endif

        ! we found the progvar entry and set the values.  return here.
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

!------------------------------------------------------------

subroutine get_index_range_string(string,index1,indexN)

! Determine where a particular DART kind (string) exists in the
! DART state vector.

character(len=*),  intent(in)  :: string
integer,           intent(out) :: index1
integer, optional, intent(out) :: indexN

integer :: i

index1 = 0
if (present(indexN)) indexN = 0

FieldLoop : do i=1,nfields
   if (progvar(i)%kind_string /= trim(string)) cycle FieldLoop
   index1 = progvar(i)%index1
   if (present(indexN)) indexN = progvar(i)%indexN
   exit FieldLoop
enddo FieldLoop

if (index1 == 0) then
   write(string1,*) 'Problem, cannot find indices for '//trim(string)
   call error_handler(E_ERR,'get_index_range_string',string1,source,revision,revdate)
endif
end subroutine get_index_range_string


!------------------------------------------------------------------

subroutine get_index_range_int(dartkind,index1,indexN)

! Determine where a particular DART kind (integer) exists in the
! DART state vector.

integer    ,           intent(in)  :: dartkind
integer(i8),           intent(out) :: index1
integer(i8), optional, intent(out) :: indexN

integer :: i
character(len=obstypelength) :: string

index1 = 0
if (present(indexN)) indexN = 0

FieldLoop : do i=1,nfields
   if (progvar(i)%dart_kind /= dartkind) cycle FieldLoop
   index1 = progvar(i)%index1
   if (present(indexN)) indexN = progvar(i)%indexN
   exit FieldLoop
enddo FieldLoop

string = get_name_for_quantity(dartkind)

if (index1 == 0) then
   write(string1,*) 'Problem, cannot find indices for kind ',dartkind,trim(string)
   call error_handler(E_ERR,'get_index_range_int',string1,source,revision,revdate)
endif

end subroutine get_index_range_int


!------------------------------------------------------------------

function get_progvar_index_from_kind(dartkind)

! Determine what index a particular DART kind (integer) is in the
! progvar array.
integer :: get_progvar_index_from_kind
integer, intent(in) :: dartkind

integer :: i

FieldLoop : do i=1,nfields
   if (progvar(i)%dart_kind /= dartkind) cycle FieldLoop
   get_progvar_index_from_kind = i
   return
enddo FieldLoop

get_progvar_index_from_kind = -1

end function get_progvar_index_from_kind


!------------------------------------------------------------------

function get_index_from_varname(varname)

! Determine what index corresponds to the given varname
! if name not in state vector, return -1 -- not an error.

integer :: get_index_from_varname
character(len=*), intent(in) :: varname

integer :: i

FieldLoop : do i=1,nfields
   if (trim(progvar(i)%varname) == trim(varname)) then
      get_index_from_varname = i
      return
   endif
enddo FieldLoop

get_index_from_varname = -1
return

end function get_index_from_varname

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
   call convert_vert_distrib(state_handle, ens_size, new_location, QTY_TEMPERATURE, VERTISPRESSURE, istatus)

   do e = 1, ens_size
      if(istatus(e) == 0) then
         llv_new = get_location(new_location(e))
         ploc(e) = llv_new(3)
      endif
   enddo

else if(is_vertical(location, "SURFACE")) then

   ivars(1) = get_progvar_index_from_kind(QTY_SURFACE_PRESSURE)
   if ( ivars(1) >= 0 ) then

     call compute_scalar_with_barycentric(state_handle, ens_size, location, 1, ivars(1), values(1,:), istatus)
     where (istatus == 0) ploc(:) = values(1,:)

   else
     istatus = 88    ! required quantity not in state vector
     return

!%!     !>original code:
!%!     !>@todo FIXME: do we really want to do this if the vert is surface and
!%!     !> the surface pressure field is not in the state?  this is going to return
!%!     !> the pressure at the midpoint of the first level, is it not?
!%!
!%!     new_location(1) = set_location(llv(1), llv(2), 1.0_r8, VERTISLEVEL)
!%!
!%!     ! Need to get base offsets for the potential temperature, density, and water
!%!     ! vapor mixing fields in the state vector
!%!     ivars(1) = get_progvar_index_from_kind(QTY_POTENTIAL_TEMPERATURE)
!%!     ivars(2) = get_progvar_index_from_kind(QTY_DENSITY)
!%!     ivars(3) = get_progvar_index_from_kind(QTY_VAPOR_MIXING_RATIO)
!%!
!%!     call compute_scalar_with_barycentric (state_handle, ens_size, new_location(1), 3, ivars, values, istatus)
!%!     if ( all(istatus /= 0) ) return
!%!
!%!     ! Convert surface theta, rho, qv into pressure
!%!     call compute_full_pressure(ens_size, values(1, :), values(2, :), values(3, :), ploc(:), tk(:), istatus(:) )

   endif

else if(is_vertical(location, "UNDEFINED")) then    ! not error, but no exact vert loc either
   ploc(:) = 200100.0_r8    ! this is an unrealistic pressure value to indicate no known pressure.
   istatus(:) = 0           ! see comment at top of this routine for why this is ok.

else
   call error_handler(E_ERR, 'compute_pressure:', 'internal error: unknown type of vertical', &
        source, revision, revdate)
endif

if(debug > 9 .and. do_output()) then
   print *, 'compute_pressure_at_loc: base location ',llv(1:3)
   print *, 'compute_pressure_at_loc: istatus, pressure ', istatus(1), ploc(1)
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
   call convert_vert_distrib(state_handle, 1, locs(i:i), loc_qtys(i), which_vert, status(i:i))
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
!> modified convert_vert_distrib() routine that doesn't
!> have to search for the cell centers. 

istatus = 0

do i=1, num
   call convert_vert_distrib_state(state_handle, 1, locs(i:i), loc_qtys(i), &
                                   loc_indx(i), which_vert, status)

   ! save the first error we see - but continue to convert the rest
   if (istatus == 0 .and. status(1) /= 0) istatus = status(1)
enddo


end subroutine convert_vertical_state


!------------------------------------------------------------------
!> code to convert an observation location's vertical coordinate type.

subroutine convert_vert_distrib(state_handle, ens_size, location, obs_kind, ztypeout, istatus)

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
!            the state that is supplied to convert_vert_distrib should be the
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
else
   if (debug > 9 .and. do_output()) then
      write(string1,'(A,3X,2I3)') 'ztypein, ztypeout:',ztypein,ztypeout
      call error_handler(E_MSG, 'convert_vert_distrib',string1,source, revision, revdate)
   endif
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

   if (debug > 9 .and. do_output()) then
      write(string2,'("Zk for member 1:", 3F8.2," => ",F8.2)') zk_mid(:,1),zout(1)
      call error_handler(E_MSG, 'convert_vert_distrib',string2,source, revision, revdate)
   endif

   ! ------------------------------------------------------------
   ! outgoing vertical coordinate should be 'pressure' in Pa
   ! ------------------------------------------------------------
   case (VERTISPRESSURE)

   ! Need to get base offsets for the potential temperature, density, and water
   ! vapor mixing fields in the state vector
   ivars(1) = get_progvar_index_from_kind(QTY_POTENTIAL_TEMPERATURE)
   ivars(2) = get_progvar_index_from_kind(QTY_DENSITY)
   ivars(3) = get_progvar_index_from_kind(QTY_VAPOR_MIXING_RATIO)

   if (any(ivars(1:3) < 0)) then
      write(string1,*) 'Internal error, cannot find one or more of: theta, rho, qv'
      call error_handler(E_ERR, 'convert_vert_distrib',string1,source, revision, revdate)
   endif

   ! Get theta, rho, qv at the interpolated location
   call compute_scalar_with_barycentric (state_handle, ens_size, location(1), 3, ivars, values, istatus)
   if( all(istatus /= 0) ) then
      location(:) = set_location(llv_loc(1, 1),llv_loc(2, 1),missing_r8,ztypeout)
      return
   endif

   ! Convert theta, rho, qv into pressure
   call compute_full_pressure(ens_size, values(1,:), values(2,:), values(3,:), zout, tk, istatus)
   if( all(istatus /= 0) ) then
      if (debug > 1 .and. do_output()) then
       write(string2,'("Failed in compute_full_pressure in VERTISPRESSURE")')
       call error_handler(E_ALLMSG, 'convert_vert_distrib',string2,source, revision, revdate)
      endif
   endif

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

   if (debug > 9 .and. do_output()) then
      write(string2,'("zout_in_height [m] for member 1:",F10.2)') zout(1)
      call error_handler(E_MSG, 'convert_vert_distrib',string2,source, revision, revdate)
   endif


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
        goto 101
     endif

     ! Base offsets for the potential temperature, density, and water
     ! vapor mixing fields in the state vector.
       ivars(1) = get_progvar_index_from_kind(QTY_POTENTIAL_TEMPERATURE)
       ivars(2) = get_progvar_index_from_kind(QTY_DENSITY)
       ivars(3) = get_progvar_index_from_kind(QTY_VAPOR_MIXING_RATIO)

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
       if( all(istatus /= 0) ) then
          if (debug > 1 .and. do_output()) then
              write(string2,'("Failed in compute_full_pressure in surfp")')
              call error_handler(E_MSG, 'convert_vert_distrib',string2,source, revision, revdate)
          endif
       endif

     endif

     if (.not. at_surf) then   ! we will need full pressure

       ! Get theta, rho, qv at the interpolated location
       call compute_scalar_with_barycentric (state_handle, ens_size, location(1), 3, ivars, values, istatus)

       ! Convert theta, rho, qv into pressure
       call compute_full_pressure(ens_size, values(1,:), values(2,:), values(3,:), fullp, tk, istatus)
       if( all(istatus /= 0) ) then
          if (debug > 1 .and. do_output()) then
              write(string2,'("Failed in compute_full_pressure in fullp")')
              call error_handler(E_MSG, 'convert_vert_distrib',string2,source, revision, revdate)
          endif
       endif

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

101 continue

     if (debug > 9 .and. do_output()) then
       write(string2,'("zout_in_scaleheight for member 1:",F10.2)') zout(1)
       call error_handler(E_MSG, 'convert_vert_distrib',string2,source, revision, revdate)
     endif

   ! -------------------------------------------------------
   ! outgoing vertical coordinate is unrecognized
   ! -------------------------------------------------------
   case default
      write(string1,*) 'Requested vertical coordinate not recognized: ', ztypeout
      call error_handler(E_ERR,'convert_vert_distrib', string1, &
                         source, revision, revdate)

end select   ! outgoing vert type

! Returned location
do e = 1, ens_size
   if(istatus(e) == 0) then
      location(e) = set_location(llv_loc(1, e),llv_loc(2, e),zout(e),ztypeout)
   else
      location(e) = set_location(llv_loc(1, e),llv_loc(2, e),missing_r8,ztypeout)
   endif
enddo


end subroutine convert_vert_distrib

!------------------------------------------------------------------
!> code to convert an state location's vertical coordinate type.

subroutine convert_vert_distrib_state(state_handle, ens_size, location, quantity, state_indx, ztypeout, istatus)

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
!            state_handle that is supplied to convert_vert_distrib should be the
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
else
   if (debug > 9 .and. do_output()) then
      write(string1,'(A,3X,2I3)') 'ztypein, ztypeout:',ztypein,ztypeout
      call error_handler(E_MSG, 'convert_vert_distrib_state',string1,source, revision, revdate)
   endif
endif

!> assume that all locations have the same incoming lat/lon and level.
!> depending on the output vert type each member might have a different
!> vertical value.

! unpack the incoming location(s)
do e = 1, ens_size
   llv_loc(:, e) = get_location(location(e))
enddo

!> state_indx is i8, intent(in).  
!> cellid and vert_level are integer, intent(out)
call find_mpas_indices(state_indx, cellid, vert_level, ndim)
if (debug > 9 .and. do_output()) print*,'convert_vert_distrib_state: ndim=', ndim

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

   if (debug > 9 .and. do_output()) then
      write(string2,'("zout_in_level for member 1:",F10.2)') zout(1)
      call error_handler(E_MSG, 'convert_vert_distrib_state',string2,source, revision, revdate)
   endif

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
   ivars(1) = get_progvar_index_from_kind(QTY_POTENTIAL_TEMPERATURE)
   ivars(2) = get_progvar_index_from_kind(QTY_DENSITY)
   ivars(3) = get_progvar_index_from_kind(QTY_VAPOR_MIXING_RATIO)

   if (any(ivars(1:3) < 0)) then
      write(string1,*) 'Internal error, cannot find one or more of: theta, rho, qv'
      call error_handler(E_ERR, 'convert_vert_distrib_state',string1,source, revision, revdate)
   endif

   ! Get theta, rho, qv at the interpolated location - pass in cellid we have already located 
   ! to save the search time.
   call compute_scalar_with_barycentric (state_handle, ens_size, location(1), 3, ivars, values, istatus, cellid)
   if( all(istatus /= 0) ) then
      location(:) = set_location(llv_loc(1, 1),llv_loc(2, 1),missing_r8,ztypeout)
      return
   endif

   ! Convert theta, rho, qv into pressure
   call compute_full_pressure(ens_size, values(1,:), values(2,:), values(3,:), zout, tk, istatus)
   if( all(istatus /= 0) ) then
       if (debug > 1 .and. do_output()) then
           write(string2,'("Failed in compute_full_pressure")')
           call error_handler(E_MSG, 'convert_vert_distrib_state',string2,source, revision, revdate)
       endif
   else
       if (debug > 9 .and. do_output()) then
           write(string2,'("zout[Pa] for member 1,theta,rho,qv,ier:",3F10.2,F12.8,I5)') zout(1),values(1:3,1),istatus(1)
           call error_handler(E_MSG, 'convert_vert_distrib_state',string2,source, revision, revdate)
       endif
   endif      ! istatus

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

   if (debug > 9 .and. do_output()) then
      write(string2,'("zout[m] for member 1:",F10.2)') zout(1)
      call error_handler(E_MSG, 'convert_vert_distrib_state',string2,source, revision, revdate)
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
        goto 101
     endif

     ! Base offsets for the potential temperature, density, and water
     ! vapor mixing fields in the state vector.
       ivars(1) = get_progvar_index_from_kind(QTY_POTENTIAL_TEMPERATURE)
       ivars(2) = get_progvar_index_from_kind(QTY_DENSITY)
       ivars(3) = get_progvar_index_from_kind(QTY_VAPOR_MIXING_RATIO)

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
       if( all(istatus /= 0) ) then
           if (debug > 1 .and. do_output()) then
           write(string2,'("Failed in compute_full_pressure")')
           call error_handler(E_MSG, 'convert_vert_distrib_state',string2,source, revision, revdate)
           endif
       else
           if (debug > 9 .and. do_output()) then
           write(string2,'("zout_psfc,theta,rho,qv:",3F10.2,F12.8)') surfp(1), values(1:3,1)
           call error_handler(E_MSG, 'convert_vert_distrib_state',string2,source, revision, revdate)
           endif
       endif      ! istatus

     endif

     if (.not. at_surf) then   ! we will need full pressure

        ! Get theta, rho, qv at the interpolated location
        call compute_scalar_with_barycentric (state_handle, ens_size, location(1), 3, ivars, values, istatus)

        ! Convert theta, rho, qv into pressure
        call compute_full_pressure(ens_size, values(1,:), values(2,:), values(3,:), fullp, tk, istatus)
        if( all(istatus /= 0) ) then
           if (debug > 1 .and. do_output()) then
           write(string2,'("Failed in compute_full_pressure")')
           call error_handler(E_MSG, 'convert_vert_distrib_state',string2,source, revision, revdate)
           endif
        else
           if (debug > 9 .and. do_output()) then
           write(string2,'("zout [Pa] for member 1,theta,rho,qv:",3F10.2,F12.8)') fullp(1), values(1:3,1)
           call error_handler(E_MSG, 'convert_vert_distrib_state',string2,source, revision, revdate)
           endif
        endif      ! istatus

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

101 continue

     if (debug > 9 .and. do_output()) then
        write(string2,'("zout_in_scaleheight for member 1:",F10.2)') zout(1)
        call error_handler(E_MSG, 'convert_vert_distrib_state',string2,source, revision, revdate)
     endif

   ! -------------------------------------------------------
   ! outgoing vertical coordinate is unrecognized
   ! -------------------------------------------------------
   case default
      write(string1,*) 'Requested vertical coordinate not recognized: ', ztypeout
      call error_handler(E_ERR,'convert_vert_distrib_state', string1, &
                         source, revision, revdate)

end select   ! outgoing vert type

! Returned location
do e = 1, ens_size
   if(istatus(e) == 0) then
      location(e) = set_location(llv_loc(1, e),llv_loc(2, e),zout(e),ztypeout)
   else
      location(e) = set_location(llv_loc(1, e),llv_loc(2, e),missing_r8,ztypeout)
   endif
enddo


end subroutine convert_vert_distrib_state

!-------------------------------------------------------------------


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
      if ((debug > 1) .and. do_output()) print*,'find_height_bounds: height < bounds(1) ', &
      height, bounds(1)
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
integer,             intent(out) :: lower(:, :), upper(:, :) ! ens_size
real(r8),            intent(out) :: fract(:, :)
integer,             intent(out) :: ier(:)

real(r8) :: lat, lon, vert, llv(3)
real(r8) :: vert_array(ens_size)
integer  :: track_ier(ens_size)
integer(i8) :: pt_base_offset, density_base_offset, qv_base_offset
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

! these first 3 types need no cell/edge location information.
if ((debug > 9) .and. do_output()) then
   write(string2,'("vert, which_vert:",3F20.12,I5)') lon,lat,vert,verttype
   call error_handler(E_MSG, 'find_vert_level',string2,source, revision, revdate)
endif

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

! ok, now we need to know where we are in the grid for heights or pressures
! as the vertical coordinate.

! Vertical interpolation for pressure coordinates
if(verttype == VERTISPRESSURE ) then

   track_ier = 0
   vert_array = vert

   ! Need to get base offsets for the potential temperature, density, and water
   ! vapor mixing fields in the state vector
   call get_index_range(QTY_POTENTIAL_TEMPERATURE, pt_base_offset)
   call get_index_range(QTY_DENSITY, density_base_offset)
   call get_index_range(QTY_VAPOR_MIXING_RATIO, qv_base_offset)

   do i=1, nc
      call find_pressure_bounds(state_handle, ens_size, vert_array, ids(i), nVertLevels, &
            pt_base_offset, density_base_offset, qv_base_offset,  &
            lower(i, :), upper(i, :), fract(i, :), ier)

      !if(debug > 9) print '(A,5I5,F10.4)', &
      !   '    after find_pressure_bounds: ier, i, cellid, lower, upper, fract = ', &
      !      ier, i, ids(i), lower(i, :), upper(i, :), fract(i, :)
      !if(debug > 5 .and. ier(e) /= 0) print '(A,4I5,F10.4,I3,F10.4)', &
      !'fail in find_pressure_bounds: ier, nc, i, id, vert, lower, fract: ', &
      !   ier, nc, i, ids(i), vert_array, lower(i, :), fract(i, :)

      ! we are inside a loop over each corner. consolidate error codes
      ! so that we return an error for that ensemble member if any
      ! of the corners fails the pressure bounds test.
      where (ier /= 0 .and. track_ier == 0) track_ier = ier

   enddo

   ier = track_ier
   return
endif

!HK I believe height is the same across the ensemble, so you can use
! vert and not vert_array:

! grid is in height, so this needs to know which cell to index into
! for the column of heights and call the bounds routine.
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
                              source, revision, revdate)
         endif
         call find_height_bounds(vert, nVertLevels, zGridEdge(:, ids(i)), &
                                 lower(i, :), upper(i, :), fract(i, :), ier)
      endif

      ! we are inside a loop over each corner. consolidate error codes
      ! so that we return an error for that ensemble member if any
      ! of the corners fails the pressure bounds test.
      where (ier /= 0 .and. track_ier == 0) track_ier = ier

   enddo

   if ((debug > 9) .and. do_output()) print*,'find_vert_level is done'

   ier = track_ier
   return

endif

end subroutine find_vert_level

!------------------------------------------------------------------

subroutine find_pressure_bounds(state_handle, ens_size, p, cellid, nbounds, &
   pt_base_offset, density_base_offset, qv_base_offset, &
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
integer,     intent(in)  :: nbounds ! number of vertical levels?
integer(i8), intent(in)  :: pt_base_offset, density_base_offset, qv_base_offset
integer,     intent(out) :: lower(:), upper(:) ! ens_size
real(r8),    intent(out) :: fract(:) ! ens_size
integer,     intent(out) :: ier(:) ! ens_size

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
call get_interp_pressure(state_handle, ens_size, pt_base_offset, density_base_offset, &
   qv_base_offset, cellid, 1, nbounds, pressure(1, :), temp_ier)
if(debug > 10 .and. do_output()) &
   print *, 'find_pressure_bounds k=1: p, temp_ier, ier', 1, pressure(1,:), temp_ier, ier

where(ier(:) == 0) ier(:) = temp_ier(:)

! Get the highest pressure level
call get_interp_pressure(state_handle, ens_size, pt_base_offset, density_base_offset, &
   qv_base_offset, cellid, nbounds, nbounds, pressure(nbounds, :), temp_ier)
if(debug > 10 .and. do_output()) &
   print *, 'find_pressure_bounds k=N: p, temp_ier, ier', nbounds, pressure(nbounds,:), temp_ier, ier

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
      call get_interp_pressure(state_handle, ens_size, pt_base_offset, density_base_offset, &
         qv_base_offset, cellid, i, nbounds, pressure(i, :), temp_ier)
      if(debug > 11 .and. do_output()) &
         print *, 'find_pressure_bounds k=?: p, temp_ier, ier ', i, pressure(i,:), temp_ier, ier

      where (ier(:) == 0) ier(:) = temp_ier(:)
   endif
   
   ! Check if pressure is not monotonically descreased with level.
   if(any(pressure(i, :) > pressure(i-1, :))) then
      if (debug > 0 .and. do_output()) then
         write(*, *) 'pressure at the level is larger than the level below at cellid',  cellid
         do e=1, ens_size
            write(*, '(A,3I4,F9.2,A,F9.2)') &
               'ens#, level, level-1, p(level), p(level-1) ', e,i,i-1,pressure(i,e),' ?>? ',pressure(i-1,e)
         enddo
      endif
     where(pressure(i, :) > pressure(i-1, :)) ier(:) = 988
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

         if ((debug > 9) .and. do_output()) print '(A,3F10.2,2I4,F10.2)', &
            "find_pressure_bounds: p_in, pr(i-1), pr(i), lower, upper, fract = ", &
            p(e), pressure(i-1,e), pressure(i,e), lower(e), upper(e), fract(e)

      endif
   enddo
   if(all(found_level)) return
enddo

end subroutine find_pressure_bounds

!------------------------------------------------------------------

subroutine get_interp_pressure(state_handle, ens_size, pt_offset, density_offset, qv_offset, &
   cellid, lev, nlevs, pressure, ier)

! Finds the value of pressure at a given point at model level lev

type(ensemble_type), intent(in) :: state_handle
integer,      intent(in)  :: ens_size
integer(i8),  intent(in)  :: pt_offset, density_offset, qv_offset
integer,      intent(in)  :: cellid
integer,      intent(in)  :: lev, nlevs
real(r8),     intent(out) :: pressure(:)
integer,      intent(out) :: ier(:)

integer  :: offset
real(r8) :: pt(ens_size), density(ens_size), qv(ens_size), tk(ens_size)
integer :: e

if(debug > 9 .and. do_output()) print*, 'get_interp_pressure for k = ',lev

! Get the values of potential temperature, density, and vapor
offset = (cellid - 1) * nlevs + lev - 1
pt      =  get_state(pt_offset      + offset, state_handle)
density =  get_state(density_offset + offset, state_handle)
qv      =  get_state(qv_offset      + offset, state_handle)

! Initialization
ier = 0

! Error if any of the values are missing; probably will be all or nothing
do e = 1, ens_size
   if(pt(e) == MISSING_R8 .or. density(e) == MISSING_R8 .or. qv(e) == MISSING_R8) then
      ier(e) = 2
      pressure(e) = MISSING_R8
   endif
enddo

! Convert theta, rho, qv into pressure
call compute_full_pressure(ens_size, pt(:), density(:), qv(:), pressure(:), tk(:), ier) !HK where is tk used?
if( all(ier/= 0) ) then
    if (debug > 1 .and. do_output()) then
        write(string2,'("Failed in compute_full_pressure")')
        call error_handler(E_ALLMSG, 'get_interp_pressure',string2,source, revision, revdate)
    endif
else
    if((debug > 11) .and. do_output()) then
    write(*,*) 'get_interp_pressure: cellid, lev', cellid, lev
    write(*,*) 'get_interp_pressure: theta,rho,qv,p,tk', pt(1), density(1), qv(1), pressure(1), tk(1)
    endif
endif   !( all(istatus /= 0) )

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
integer(i8) :: index1, low_offset, upp_offset, low_state_indx, upp_state_indx
integer     :: lower(3, ens_size), upper(3,ens_size)
integer     :: e, e2, thislower, thisupper
logical     :: did_member(ens_size)

! assume failure
dval = MISSING_R8
ier = 88   ! field not in state vector

! make sure we have all good field indices first
if (any(ival < 0)) return

call find_triangle (loc, nc, c, weights, ier(1), this_cellid)
if(ier(1) /= 0) then
   ier(:) = ier(1)
   return
endif

! If the field is on a single level, lower and upper are both 1
call find_vert_indices (state_handle, ens_size, loc, nc, c, lower, upper, fract, ier)
if(all(ier /= 0)) return

! for each field to compute at this location:
do k=1, n
   ! get the starting index in the state vector
   index1 = progvar(ival(k))%index1
   nvert  = progvar(ival(k))%numvertical

   ! for each corner: could be 1 if location is directly at a vertex, but
   ! normally is 3 for the enclosing triangle made up of cell centers.
   do i = 1, nc
      ! go around triangle and interpolate in the vertical
      ! c(3) are the cell ids

      low_offset = (c(i)-1) * nvert !FIXME low_offset and upp_offset are the same?
      upp_offset = (c(i)-1) * nvert

      if( nvert == 1 ) then         ! fields on a surface (1-D, one level, ...)

          ! Because 'lower(i,1)-1' evaluates to zero for 1-D, no need to add it
          fdata(i,:) = get_state(index1 + low_offset, state_handle)

      else                          ! 2-D fields

         did_member(:) = .false.

         members: do e = 1, ens_size

            if (did_member(e)) cycle members

            !> minimize the number of times we call get_state() by
            !> doing all the ensemble members which are between the same
            !> two vertical levels.  this is true most of the time. 
            !> in some cases it could be 2 or 3 different pairs of levels because 
            !> of differences in vertical conversion that depends on per-member fields.
            low_state_indx = index1 + low_offset + lower(i,e)-1
            lowval(:) = get_state(low_state_indx, state_handle)
            upp_state_indx = index1 + upp_offset + upper(i,e)-1
            uppval(:) = get_state(upp_state_indx, state_handle)

            thislower = lower(i, e)
            thisupper = upper(i, e)

            ! for all remaining ensemble members, use these values if the lower and
            ! upper level numbers are the same.  fract() will vary with member.
            do e2=e, ens_size
               if (thislower == lower(i, e2) .and. thisupper == upper(i, e2)) then
                  fdata(i, e2) = lowval(e2)*(1.0_r8 - fract(i, e2)) + uppval(e2)*fract(i, e2)
                  did_member(e2) = .true.
               endif
            enddo
         enddo members

      endif                         ! 2-D fields

   enddo  ! corners

   ! now have vertically interpolated values at cell centers.
   ! use weights to compute value at interp point.
   do e = 1, ens_size
      if (ier(e) /= 0) cycle
      dval(k, e) = sum(weights(1:nc) * fdata(1:nc, e))
   enddo

enddo

end subroutine compute_scalar_with_barycentric

!------------------------------------------------------------
!> Interpolates metadata variables to an arbitrary location.
!> This routine can only be used with variables in module storage.

subroutine compute_surface_data_with_barycentric(var1d, loc, dval, ier, this_cellid)

type(location_type), intent(in)  :: loc
real(r8),            intent(in)  :: var1d(:)
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

if(debug > 9 .and. do_output()) &
   print '(A,7f12.5)','compute_surface_data_with_barycentric: corner vals, weights, result: ',&
                      fdata(:),weights(:),dval

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

if ((xyzdebug > 1) .and. do_output()) &
   print *, 'ft: closest cell center for lon/lat: ', lon, lat, cellid

if (cellid < 1) then
   if(xyzdebug > 0) print *, 'ft: closest cell center for lon/lat: ', lon, lat, cellid
   ier = 11
   return
endif

c(1) = cellid

! closest vertex to given point.
closest_vert = closest_vertex_ll(cellid, lat, lon)
if ((xyzdebug > 5) .and. do_output()) &
   print *, 'ft: closest vertex for lon/lat: ', lon, lat, closest_vert

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
!print *, 'ft: i: ', i
   edgeid = edgesOnCell(i, cellid)
!print *, 'ft: edgeid: ', edgeid
   if (.not. global_grid .and. &
      (cellsOnEdge(1, edgeid) <= 0 .or. cellsOnEdge(2, edgeid) <= 0)) then
      ier = 14
      return
   endif
   if (cellsOnEdge(1, edgeid) /= cellid) then
      neighborcells(i) = cellsOnEdge(1, edgeid)
   else
      neighborcells(i) = cellsOnEdge(2, edgeid)
   endif
   verts(i) = verticesOnCell(i, cellid)
!print *, 'ft: verts: ', verts(i), closest_vert
   if (verts(i) == closest_vert) vindex = i
   call latlon_to_xyz(latCell(neighborcells(i)), lonCell(neighborcells(i)), &
      xdata(i), ydata(i), zdata(i))
enddo


! get the cartesian coordinates in the cell plane for the closest center
call latlon_to_xyz(latCell(cellid), lonCell(cellid), t1(1), t1(2), t1(3))

! and the observation point
call latlon_to_xyz(lat, lon, r(1), r(2), r(3))
!print *, 'ft: lat/lon: ', lat, lon

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
      ier = 14     ! 11?
      return
   endif

endif     ! horizontal index search is done now.
if (ier /= 0) return

if (debug > 12 .and. do_output()) then
   write(string3,*) 'ier = ',ier, ' triangle = ',c(1:nc), ' weights = ',weights(1:nc)
   call error_handler(E_MSG, 'find_triangle', string3, source, revision, revdate)
endif

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
integer,             intent(in)  :: c(:)
integer,             intent(out) :: lower(:, :), upper(:, :) ! ens_size
real(r8),            intent(out) :: fract(:, :) ! ens_size
integer,             intent(out) :: ier(:) ! ens_size

integer :: e

! initialization
lower = MISSING_I
upper = MISSING_I
fract = 0.0_r8

! need vert index for the vertical level
call find_vert_level(state_handle, ens_size, loc, nc, c, .true., lower, upper, fract, ier)

if (debug > 9 .and. do_output()) then
   write(string3,*) 'ier = ',ier(1), ' triangle = ',c(1:nc), ' vert_index = ',lower(1:nc, 1)+fract(1:nc, 1), ' nc = ', nc
   call error_handler(E_MSG, 'find_vert_indices', string3, source, revision, revdate)
endif

if (debug > 12 .and. do_output()) then
  do e = 1, ens_size
   if(ier(e) /= 0) then   
      print *, 'find_vert_indices: e = ', e, ' nc = ', nc, ' ier = ', ier(e)
      print *, 'find_vert_indices: c = ', c
      print *, 'find_vert_indices: lower = ', lower(1:nc, e)
      print *, 'find_vert_indices: upper = ', upper(1:nc, e)
      print *, 'find_vert_indices: fract = ', fract(1:nc, e)
   endif
  enddo
endif

end subroutine find_vert_indices

!------------------------------------------------------------

subroutine compute_u_with_rbf(state_handle, ens_size, loc, zonal, uval, ier)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: loc
logical,             intent(in)  :: zonal
real(r8),            intent(out) :: uval(:) ! ens_size
integer,             intent(out) :: ier(:) ! ens_size

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
integer(i8) :: upindx, lowindx
integer     :: index1, progindex, cellid, vertexid
real(r8)    :: lat, lon, vert, llv(3), fract(listsize, ens_size), lowval(ens_size), uppval(ens_size)
integer     :: verttype, lower(listsize, ens_size), upper(listsize, ens_size), ncells, celllist(listsize)

integer :: e ! loop index

! Initialization
ier  = 0
uval = MISSING_R8

! FIXME: it would be great to make this cache the last value and 
! if the location is the same as before and it's asking for V now 
! instead of U, skip the expensive computation.  however, given
! how we currently distribute observations the V wind obs will
! almost certainly be given to a different task.  if that changes
! in some future version of dart, revisit this code as well.

progindex = get_index_from_varname('u')
if (progindex < 0 .or. .not. data_on_edges) then
   ! cannot compute u if it isn't in the state vector, or if we
   ! haven't read in the edge data (which shouldn't happen if
   ! u is in the state vector.
   ier = 18
   return
endif
index1 = progvar(progindex)%index1
nvert = progvar(progindex)%numvertical

! unpack the location into local vars
llv = get_location(loc) ! I believe this is the same across the ensemble
lon = llv(1)
lat = llv(2)
vert = llv(3)
verttype = nint(query_location(loc))

call find_surrounding_edges(lat, lon, nedges, edgelist, cellid, vertexid)
if (nedges <= 0) then
   ! we are on a boundary, no interpolation
   ier = 18
   return
endif

if (verttype == VERTISPRESSURE) then
   ! get all cells which share any edges with the 3 cells which share
   ! the closest vertex.
   call make_cell_list(vertexid, 3, ncells, celllist)

   call find_vert_level(state_handle, ens_size, loc, ncells, celllist, .true., &
                        lower, upper, fract, ier)

   if (all(ier /= 0)) return

   ! now have pressure at all cell centers - need to interp to get pressure
   ! at edge centers.
   call move_pressure_to_edges(ncells, celllist, lower, upper, fract, &
                               nedges, edgelist, ier)
   if (all(ier /= 0)) return

else
   ! need vert index for the vertical level
   call find_vert_level(state_handle, ens_size, loc, nedges, edgelist, .false., &
                        lower, upper, fract, ier)
   if (all(ier /= 0)) return
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
                         source, revision, revdate, text2=string1)
   endif
   xdata(i) = xEdge(edgelist(i))
   ydata(i) = yEdge(edgelist(i))
   zdata(i) = zEdge(edgelist(i))

   do j=1, 3
      edgenormals(j, i) = edgeNormalVectors(j, edgelist(i))
   enddo

   !>@todo lower (upper) could be different levels in pressure
   !lowval = x(index1 + (edgelist(i)-1) * nvert + lower(i)-1)
   lowindx = int(index1,i8) + int((edgelist(i)-1) * nvert,i8) + int(lower(i,1)-1,i8)
   lowval =  get_state(lowindx, state_handle)

   !uppval = x(index1 + (edgelist(i)-1) * nvert + upper(i)-1)
   upindx = int(index1,i8) + int((edgelist(i)-1) * nvert,i8) + int(upper(i,1)-1,i8)
   uppval =  get_state(upindx, state_handle)

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
      uval = ureconstructzonal
   else
      uval = ureconstructmeridional
   endif

enddo

end subroutine compute_u_with_rbf

!------------------------------------------------------------

subroutine find_surrounding_edges(lat, lon, nedges, edge_list, cellid, vertexid)
real(r8), intent(in)  :: lat, lon
integer,  intent(out) :: nedges, edge_list(:)
integer,  intent(out) :: cellid, vertexid

! given an arbitrary lat/lon location, find the edges of the
! cells that share the nearest vertex.  return the ids for the
! closest cell center and closest vertex to save having to redo
! the search.

integer :: vertex_list(30), nedges1, edge_list1(300), c, i
integer :: ncells, cell_list(150), nedgec, edgeid, nverts

! find the cell id that has a center point closest
! to the given point.
cellid = find_closest_cell_center(lat, lon)
if (cellid < 1) then
   nedges = 0
   edge_list(:) = -1
   return
endif

! inside this cell, find the vertex id that the point
! is closest to.  this is a return from this subroutine.
vertexid = closest_vertex_ll(cellid, lat, lon)
if (vertexid <= 0) then
   ! call error handler?  unexpected
   nedges = 0
   edge_list(:) = -1
   return
endif

nedges = 0
edge_list = 0

select case (use_rbf_option)
  case (0)
      ! use edges on the closest cell only.
      nedges = nEdgesOnCell(cellid)
      do i=1, nedges
         edge_list(i) = edgesOnCell(i, cellid)
      enddo
  case (1)
      ! fill in the number of unique edges and fills the
      ! edge list with the edge ids.  the code that detects
      ! boundary edges for the ocean or regional atmosphere
      ! is incorporated here.  nedges can come back 0 in that case.
      vertex_list(1) = vertexid
      call make_edge_list_from_verts(1, vertex_list, nedges, edge_list)

   case (2)
      ! for all vertices of the enclosing cell, add the
      ! edges which share this vertex.
      nverts = nEdgesOnCell(cellid)
      do i=1, nverts
         vertex_list(i) = verticesOnCell(i, cellid)
      enddo

      ! fill in the number of unique edges and fills the
      ! edge list with the edge ids.  the code that detects
      ! boundary edges for the ocean or regional atmosphere
      ! is incorporated here.  nedges can come back 0 in that case.
      call make_edge_list_from_verts(nverts, vertex_list, nedges, edge_list)

   case (3)
      call make_cell_list(vertexid, 1, ncells, cell_list)

      ! for all cells:
      do c=1, ncells
         ! for all vertices of the enclosing cell, add the
         ! edges which share this vertex.
         nverts = nEdgesOnCell(cell_list(c))
         do i=1, nverts
            vertex_list(i) = verticesOnCell(i, cell_list(c))
         enddo

         ! fill in the number of unique edges and fills the
         ! edge list with the edge ids.  the code that detects
         ! boundary edges for the ocean or regional atmosphere
         ! is incorporated here.  nedges can come back 0 in that case.
         call make_edge_list_from_verts(nverts, vertex_list, nedges1, edge_list1)
         call merge_edge_lists(nedges, edge_list, nedges1, edge_list1)
      enddo

   case (4)
      ! this one gives the same edge list as case (2).  see what's faster.
      nedgec = nEdgesOnCell(cellid)
      do i=1, nedgec
         edgeid = edgesOnCell(i, cellid)
         if (cellsOnEdge(1, edgeid) /= cellid) then
            cell_list(i) = cellsOnEdge(1, edgeid)
         else
            cell_list(i) = cellsOnEdge(2, edgeid)
         endif
      enddo

      call make_edge_list_from_cells(nedgec, cell_list, nedges, edge_list)

   case default
      call error_handler(E_ERR, 'find_surrounding_edges', 'bad use_rbf_option value', &
                         source, revision, revdate)
end select

! Check if any of edges are located in the boundary zone.
! (We will skip the obs if any edges are located there.)
if (on_boundary_edgelist(edge_list)) then
   nedges = -1
   edge_list(:) = -1
   if(debug > 0 .and. do_output()) call error_handler(E_MSG, 'find_surrounding_edges', 'edges in the boundary', &
                      source, revision, revdate)
   return
endif

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
   if (debug > 0 .or. global_grid) then
       print *, 'cannot find nearest cell to lon, lat: ', lon, lat
   endif
   find_closest_cell_center = -1
   return
endif

if (debug > 9 .and. do_output()) print *, 'find_closest_cell_center: lat/lon closest to cellid, with lat/lonCell: ', &
                          lat, lon, closest_cell, latCell(closest_cell), lonCell(closest_cell)

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
call error_handler(E_MSG,'set_global_grid',string1,source,revision,revdate)

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
if ((closest_vertex_ll < 0)  .and. &
    (debug > 0) .and. do_output()) &
   print *, 'cannot find nearest vertex to lon, lat: ', lon, lat

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
                      source, revision, revdate)
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
                         source, revision, revdate, text2=string1)
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
                      source, revision, revdate, text2=string1)
endif

if (debug > 1 .and. do_output()) print *, 'extrapolate: ', h, lb, ub, extrapolate_level
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

subroutine xyz_to_latlon(x, y, z, lat, lon)

! Given a cartesian x, y, z coordinate relative to the origin
! at the center of the earth, using a fixed radius specified
! by MPAS (in the grid generation step), return the corresponding
! lat, lon location in degrees.

real(r8), intent(in)  :: x, y, z
real(r8), intent(out) :: lat, lon

real(r8) :: rlat, rlon

! right now this is only needed for debugging messages.
! the arc versions of routines are expensive.

rlat = PI/2.0_r8 - acos(z/radius)
rlon = atan2(y,x)
if (rlon < 0) rlon = rlon + PI*2

lat = rlat * rad2deg
lon = rlon * rad2deg

end subroutine xyz_to_latlon

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

function vector_magnitude(a)

! Given a cartesian vector, compute the magnitude

real(r8), intent(in)  :: a(3)
real(r8) :: vector_magnitude

vector_magnitude = sqrt(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))

end function vector_magnitude

!------------------------------------------------------------

subroutine vector_cross_product(a, b, r)

! Given 2 cartesian vectors, compute the cross product of a x b

real(r8), intent(in)  :: a(3), b(3)
real(r8), intent(out) :: r(3)

r(1) = a(2)*b(3) - a(3)*b(2)
r(2) = a(3)*b(1) - a(1)*b(3)
r(3) = a(1)*b(2) - a(2)*b(1)

end subroutine vector_cross_product

!------------------------------------------------------------

function vector_dot_product(a, b)

! Given 2 cartesian vectors, compute the dot product of a . b

real(r8), intent(in)  :: a(3), b(3)
real(r8) :: vector_dot_product

vector_dot_product = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

end function vector_dot_product

!------------------------------------------------------------

subroutine vector_projection(a, b, r)

! Given 2 cartesian vectors, project a onto b

real(r8), intent(in)  :: a(3), b(3)
real(r8), intent(out) :: r(3)

real(r8) :: ab_over_bb

ab_over_bb = vector_dot_product(a, b) / vector_dot_product(b, b)
r = (ab_over_bb) * b

end subroutine vector_projection

!------------------------------------------------------------

subroutine determinant3(a, r)

! Given a 3x3 matrix, compute the determinant

real(r8), intent(in)  :: a(3,3)
real(r8), intent(out) :: r

r = a(1,1)*(a(2,2)*a(3,3) - (a(3,2)*a(2,3))) + &
    a(2,1)*(a(3,2)*a(1,3) - (a(3,3)*a(1,2))) + &
    a(3,1)*(a(1,2)*a(2,3) - (a(2,2)*a(1,3)))

end subroutine determinant3

!------------------------------------------------------------

subroutine invert3(a, r)

! Given a 3x3 matrix, compute the inverse

real(r8), intent(in)  :: a(3,3)
real(r8), intent(out) :: r(3,3)

real(r8) :: det, b(3,3)

call determinant3(a, det)
if (det == 0.0_r8) then
   print *, 'matrix cannot be inverted'
   r = 0.0_r8
   return
endif

b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
b(2,1) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
b(3,1) = a(2,1)*a(3,2) - a(3,1)*a(2,2)

b(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
b(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
b(3,2) = a(3,1)*a(1,2) - a(1,1)*a(3,2)

b(1,3) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
b(2,3) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)

r = b / det

end subroutine invert3

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


if ( debug > 0 .and. do_output()) then
   if (any(istatus /= 0)) then
     print *, 'theta_to_tk - nonzero istatus coming in'
     do e = 1, ens_size
      if (istatus(e) /= 0) then
          write(string2, *) 'member ', e, ' incoming istatus = ', istatus(e)
          call error_handler(E_ALLMSG, 'theta_to_tk',string2,source, revision, revdate)
      endif  !(istatus(e) /= 0) then
     enddo   ! ens_size
   endif
endif
where (istatus == 0)

   theta_m = (1.0_r8 + rvord * qv_nonzero)*theta
   
   where (theta_m > 0.0_r8 .and. rho > 0.0_r8)  ! Check if all the input are positive

      exner = ( (rgas/p0) * (rho*theta_m) )**rcv
   
      ! Temperature [K]
      theta_to_tk = theta * exner

   elsewhere

      istatus = 89

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

if ( debug > 1 ) then
   if( any(istatus /= 0) ) then   
      do e = 1, ens_size
       if (istatus(e) /= 0) then
           write(string2,'("Failed in member,istatus,P[Pa],tk,theta,rho,qv:",2I4,4F12.2,F15.6)') &
           e, istatus(e), pressure(e), tk(e), theta(e), rho(e), qv(e)
           call error_handler(E_ALLMSG, 'compute_full_pressure',string2,source, revision, revdate)
       endif  !(istatus(e) /= 0) then
      enddo   ! ens_size
  endif       ! debug
else
  if ((debug > 12) .and. do_output()) print *, 'compute_full_pressure: theta,r,q,p,tk,istatus =', &
                                                theta, rho, qv, pressure, tk, istatus
endif

end subroutine compute_full_pressure

!--------------------------------------------------------------------
!> pass the vertical localization coordinate to assim_tools_mod
function query_vert_localization_coord()

integer :: query_vert_localization_coord

query_vert_localization_coord = vert_localization_coord

end function query_vert_localization_coord

!--------------------------------------------------------------------
!> read the time from the input file
!> stolen get_analysis_time_fname
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

!===================================================================
! End of model_mod
!===================================================================
end module model_mod
