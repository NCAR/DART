! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! MPAS ocean model interface to the DART data assimilation system.
! Code in this module is compiled with the DART executables.  It isolates
! all information about the MPAS grids, model variables, and other details.
! There are a set of 16 subroutine interfaces that are required by DART;
! these cannot be changed.  Additional public routines in this file can
! be used by converters and utilities and those interfaces can be anything
! that is useful to other pieces of code.

! Units on everything are MKS:
!
! u   (velocity):  meter / second
! h   (depth)   :  meter
! rho (density) :  kilograms / meter^3
! temperature   :  *potential temperature* degrees C
! salinity      :  PSU
!
! Note:  the 'temperature' variable is *potential* temperature.


! Routines in other modules that are used here.

use        types_mod, only : r4, r8, digits12, SECPERDAY, MISSING_R8,          &
                             rad2deg, deg2rad, PI, MISSING_I, obstypelength

use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, get_dist, query_location,          &
                             get_close_maxdist_init, get_close_type,           &
                             set_location, get_location, horiz_dist_only,      &
                             write_location,                                   &
                             vert_is_undef,        VERTISUNDEF,                &
                             vert_is_surface,      VERTISSURFACE,              &
                             vert_is_level,        VERTISLEVEL,                &
                             vert_is_pressure,     VERTISPRESSURE,             &
                             vert_is_height,       VERTISHEIGHT,               &
                             vert_is_scale_height, VERTISSCALEHEIGHT,          &
                             get_close_obs_init, get_close_obs_destroy,        &
                             loc_get_close_obs => get_close_obs

use xyz_location_mod, only : xyz_location_type, xyz_get_close_maxdist_init,    &
                             xyz_get_close_type, xyz_set_location, xyz_get_location, &
                             xyz_get_close_obs_init, xyz_get_close_obs_destroy, &
                             xyz_find_nearest

use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper, nmlfileunit,       &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, file_exist, find_textfile_dims,        &
                             file_to_text, close_file, do_nml_file, do_nml_term

use     obs_kind_mod, only : get_index_for_quantity,  &
                             get_name_for_quantity,   &
                             QTY_VERTICAL_VELOCITY,  &
                             QTY_POTENTIAL_TEMPERATURE, &
                             QTY_TEMPERATURE,        &
                             QTY_SALINITY,              &
                             QTY_DRY_LAND,              &
                             QTY_EDGE_NORMAL_SPEED,     &
                             QTY_U_CURRENT_COMPONENT,   &
                             QTY_V_CURRENT_COMPONENT,   &
                             QTY_SEA_SURFACE_HEIGHT,    &
                             QTY_SEA_SURFACE_PRESSURE,  &
                             QTY_TRACER_CONCENTRATION

use mpi_utilities_mod, only: my_task_id

use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

! netcdf modules
use typesizes
use netcdf

! RBF (radial basis function) modules, donated by LANL. currently deprecated
! in this version.  they did the job but not as well as other techniques and
! at a much greater execution-time code.  they were used to interpolate
! values at arbitrary locations, not just at cell centers.  with too small
! a set of basis points, the values were discontinuous at cell boundaries;
! with too many the values were too smoothed out.  we went back to
! barycentric interpolation in triangles formed by the three cell centers
! that enclosed the given point.
use get_geometry_mod
use get_reconstruct_mod


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
          ens_mean_for_model

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.

public :: get_model_analysis_filename,  &
          get_grid_definition_filename, &
          analysis_file_to_statevector, &
          statevector_to_analysis_file, &
          get_analysis_time,            &
          write_model_time,             &
          get_grid_dims,                &
          print_variable_ranges

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! module global storage; maintains values between calls, accessible by
! any subroutine
character(len=256) :: string1, string2, string3
logical, save :: module_initialized = .false.

! Real (physical) constants as defined exactly in MPAS.
! redefined here for consistency with the model.
real(r8), parameter :: rgas = 287.0_r8
real(r8), parameter :: cp = 1003.0_r8
real(r8), parameter :: cv = 716.0_r8
real(r8), parameter :: p0 = 100000.0_r8
real(r8), parameter :: rcv = rgas/(cp-rgas)

! earth radius; needed to convert lat/lon to x,y,z cartesian coords.
! FIXME: the example ocean files has a global attr with 6371220.0
! the actual mpas code may have hardwired values (it did in the atmosphere)
! need to check what is really going on.
real(r8), parameter :: radius = 6371229.0 ! meters

! roundoff error
real(r8), parameter :: roundoff = 1.0e-12_r8

! Storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq

! Structure for computing distances to cell centers, and assorted arrays
! needed for the get_close code.
type(xyz_get_close_type)             :: cc_gc
type(xyz_location_type), allocatable :: cell_locs(:)

! not part of the namelist because in the ocean it's not clear
! that we need to worry about log pressure in the vertical.
! if this makes a difference, set this to .true.
logical :: log_p_vert_interp = .false.  ! if true, interpolate vertical pressure in log space

! variables which are in the module namelist
integer            :: vert_localization_coord = VERTISHEIGHT
integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 60
real(r8)           :: model_perturbation_amplitude = 0.0001   ! tiny amounts
logical            :: output_state_vector = .false.  ! output prognostic variables (if .false.)
integer            :: debug = 0   ! turn up for more and more debug messages
integer            :: xyzdebug = 0
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: model_analysis_filename = 'mpas_analysis.nc'
character(len=256) :: grid_definition_filename = 'mpas_analysis.nc'

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

namelist /model_nml/             &
   model_analysis_filename,      &
   grid_definition_filename,     &
   vert_localization_coord,      &
   assimilation_period_days,     &
   assimilation_period_seconds,  &
   model_perturbation_amplitude, &
   calendar,                     &
   debug,                        &
   xyzdebug,                     &
   use_u_for_wind,               &
   use_rbf_option,               &
   update_u_from_reconstruct,    &
   use_increments_for_u_update

! DART state vector contents are specified in the input.nml:&mpas_vars_nml namelist.
integer, parameter :: max_state_variables = 80
integer, parameter :: num_state_table_columns = 2
integer, parameter :: num_bounds_table_columns = 4
character(len=NF90_MAX_NAME) :: mpas_state_variables(max_state_variables * num_state_table_columns ) = ' '
character(len=NF90_MAX_NAME) :: mpas_state_bounds(num_bounds_table_columns, max_state_variables ) = ' '
character(len=NF90_MAX_NAME) :: variable_table(max_state_variables, num_state_table_columns )

namelist /mpas_vars_nml/ mpas_state_variables, mpas_state_bounds

! FIXME: this shouldn't be a global.  the progvar array
! should be allocated at run time and nfields should be part
! of a larger derived type that includes nfields and an array
! of progvartypes.
integer :: nfields

! Everything needed to describe a variable

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
real(r8), allocatable :: zGridFace(:,:)   ! geometric depth at cell faces   (nVertLevelsP1,nCells)
real(r8), allocatable :: zGridCenter(:,:) ! geometric depth at cell centers (nVertLevels,  nCells)
real(r8), allocatable :: zGridEdge(:,:)   ! geometric depth at edge centers (nVertLevels,  nEdges)

real(r8), allocatable :: zMid(:,:)    ! depths at midpoints - may be able to be computed instead
                                      !  of requiring it in the input file (save file space). FIXME
real(r8), allocatable :: hZLevel(:)   ! layer thicknesses - maybe - FIXME
!real(r8), allocatable :: zEdgeCenter(:,:) ! geometric height at edges faces  (nVertLevels  ,nEdges)

integer,  allocatable :: cellsOnVertex(:,:) ! list of cell centers defining a triangle
integer,  allocatable :: verticesOnCell(:,:)

integer,  allocatable :: edgesOnCell(:,:) ! list of edges that bound each cell
integer,  allocatable :: cellsOnEdge(:,:) ! list of cells that bound each edge
integer,  allocatable :: nedgesOnCell(:) ! list of edges that bound each cell
real(r8), allocatable :: edgeNormalVectors(:,:)

! Boundary information might be needed ... regional configuration?
! Read if available.

integer,  allocatable :: boundaryCell(:,:) ! logical, cells that are on boundaries
integer,  allocatable :: maxLevelEdgeTop(:) !
integer,  allocatable :: boundaryEdge(:,:) ! logical, edges that are boundaries
integer,  allocatable :: boundaryVertex(:,:) ! logical, vertices that are on boundaries
integer,  allocatable :: maxLevelCell(:) ! list of maximum (deepest) level for each cell

real(r8), allocatable :: ens_mean(:)   ! needed to convert vertical distances consistently

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
integer :: cel1, cel2
logical :: both

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

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)

!---------------------------------------------------------------
! 1) get grid dimensions
! 2) allocate space for the grids
! 3) read them from the analysis file

! read_grid_dims() fills in the following module global variables:
!  nCells, nVertices, nEdges, maxEdges, nVertLevels, nVertLevelsP1, vertexDegree
call read_grid_dims()

allocate(latCell(nCells), lonCell(nCells))
allocate(zGridFace(nVertLevelsP1, nCells))
allocate(zGridCenter(nVertLevels, nCells))

allocate(cellsOnVertex(vertexDegree, nVertices))
allocate(nEdgesOnCell(nCells))
allocate(edgesOnCell(maxEdges, nCells))
allocate(cellsOnEdge(2, nEdges))
allocate(verticesOnCell(maxEdges, nCells))
allocate(edgeNormalVectors(3, nEdges))
allocate(xVertex(nVertices), yVertex(nVertices), zVertex(nVertices))

! see if U is in the state vector list.  if not, don't read in or
! use any of the Edge arrays to save space.
data_on_edges = is_edgedata_in_state_vector(mpas_state_variables)

if(data_on_edges) then
   allocate(zGridEdge(nVertLevels, nEdges))
   allocate(xEdge(nEdges), yEdge(nEdges), zEdge(nEdges))
   allocate(latEdge(nEdges), lonEdge(nEdges))
endif

! this reads in latCell, lonCell, hZLevel, zMid, cellsOnVertex
call get_grid()

! determine which edges are boundaries
! (this requires 'maxLevelEdgeTop' to exist in the restart file)
boundaryEdge(:,1:nEdges) = 1
do iloc=1, nEdges
   if(maxLevelEdgeTop(iloc) > 0) then
      boundaryEdge(1:maxLevelEdgeTop(iloc),iloc) = 0
   endif
end do

! determine which cells are on boundaries
boundaryCell(:,1:nCells) = 0
do kloc=1, nEdges
   do iloc=1, nVertLevels
      if(boundaryEdge(iloc,kloc) .eq. 1) then
         if(cellsOnEdge(1,kloc) > 0) then
            boundaryCell(iloc,cellsOnEdge(1,kloc)) = 1
         endif
         if(cellsOnEdge(2,kloc) > 0) then
            boundaryCell(iloc,cellsOnEdge(2,kloc)) = 1
         endif
     endif
   end do
end do

if(data_on_edges) then
   ! FIXME: This code is supposed to check whether an edge has 2 neighbours or 1 neighbour and then
   !        compute the depth accordingly.  HOWEVER, the array cellsOnEdge does not change with
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
! of the DART state vector. Required to determine model_size.
!
! Verify all variables are in the model analysis file
!
! Compute the offsets into the state vector for the start of each
! different variable type. Requires reading shapes from the model
! analysis file. As long as TIME is the LAST dimension, we're OK.
!
! Record the extent of the data type in the state vector.


call nc_check( nf90_open(trim(model_analysis_filename), NF90_NOWRITE, ncid), &
                  'static_init_model', 'open '//trim(model_analysis_filename))

call verify_state_variables( mpas_state_variables, ncid, model_analysis_filename, &
                             nfields, variable_table )

TimeDimID = FindTimeDimension( ncid )

if (TimeDimID < 0 ) then
   write(string1,*)'unable to find a dimension named Time.'
   call error_handler(E_MSG,'static_init_model', string1, source, revision, revdate)
endif

call nc_check(nf90_Inquire(ncid,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                    'static_init_model', 'inquire '//trim(model_analysis_filename))

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

   string2 = trim(model_analysis_filename)//' '//trim(varname)

   call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), &
            'static_init_model', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid, VarID, xtype=progvar(ivar)%xtype, &
           dimids=dimIDs, ndims=numdims), 'static_init_model', 'inquire '//trim(string2))

   ! If the long_name and/or units attributes are set, get them.
   ! They are not REQUIRED to exist but are nice to use if they are present.

   if( nf90_inquire_attribute(    ncid, VarID, 'long_name') == NF90_NOERR ) then
      call nc_check( nf90_get_att(ncid, VarID, 'long_name' , progvar(ivar)%long_name), &
                  'static_init_model', 'get_att long_name '//trim(string2))
   else
      progvar(ivar)%long_name = varname
   endif

   if( nf90_inquire_attribute(    ncid, VarID, 'units') == NF90_NOERR )  then
      call nc_check( nf90_get_att(ncid, VarID, 'units' , progvar(ivar)%units), &
                  'static_init_model', 'get_att units '//trim(string2))
   else
      progvar(ivar)%units = 'unknown'
   endif

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

   if ( debug > 4 ) call dump_progvar(ivar)

enddo

call nc_check( nf90_close(ncid), &
                  'static_init_model', 'close '//trim(model_analysis_filename))

model_size = progvar(nfields)%indexN

if ( debug > 0 .and. do_output()) then
  write(logfileunit,*)
  write(     *     ,*)
  write(logfileunit,'(" static_init_model: nCells, nEdges, nVertices, nVertLevels =",4(1x,i6))') &
                                          nCells, nEdges, nVertices, nVertLevels
  write(     *     ,'(" static_init_model: nCells, nEdges, nVertices, nVertLevels =",4(1x,i6))') &
                                          nCells, nEdges, nVertices, nVertLevels
  write(logfileunit, *)'static_init_model: model_size = ', model_size
  write(     *     , *)'static_init_model: model_size = ', model_size
  if ( global_grid ) then
     write(logfileunit, *)'static_init_model: grid is a global grid '
     write(     *     , *)'static_init_model: grid is a global grid '
  else
     write(logfileunit, *)'static_init_model: grid has boundaries '
     write(     *     , *)'static_init_model: grid has boundaries '
  endif
  if ( all_levels_exist_everywhere ) then
     write(logfileunit, *)'static_init_model: all cells have same number of vertical levels '
     write(     *     , *)'static_init_model: all cells have same number of vertical levels '
  else
     write(logfileunit, *)'static_init_model: cells have varying number of vertical levels '
     write(     *     , *)'static_init_model: cells have varying number of vertical levels '
  endif
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
!
! TJH if ((get_progvar_index_from_kind(QTY_POTENTIAL_TEMPERATURE) < 0) .or. &
! TJH     (get_progvar_index_from_kind(QTY_DENSITY) < 0) .or. &
! TJH     (get_progvar_index_from_kind(QTY_VAPOR_MIXING_RATIO) < 0)) then
! TJH    write(string1, *) 'State vector is missing one or more of the following fields:'
! TJH    write(string2, *) 'Potential Temperature (theta), Density (rho), Vapor Mixing Ratio (qv).'
! TJH    write(string3, *) 'Cannot convert between height/pressure nor compute sensible temps.'
! TJH    call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate, &
! TJH                       text2=string2, text3=string3)
! TJH endif

string1 = 'WARNING: fix block of required variables - detritus from atmosphere'
call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)

allocate( ens_mean(model_size) )

end subroutine static_init_model


!------------------------------------------------------------------

subroutine get_state_meta_data(index_in, location, var_type)

! given an index into the state vector, return its location and
! if given, the var kind.   despite the name, var_type is a generic
! kind, like those in obs_kind/obs_kind_mod.f90, starting with QTY_

! passed variables

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

! Local variables

integer  :: nxp, nzp, iloc, vloc, nf, n
integer  :: myindx
integer  :: istatus
real(r8) :: depth
type(location_type) :: new_location

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

if( myindx == -1 ) then
     write(string1,*) 'Problem, cannot find base_offset, index_in is: ', index_in
     call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
endif

! Now that we know the variable, find the cell or edge

if (     progvar(nf)%numcells /= MISSING_I) then
   nxp = progvar(nf)%numcells
elseif ( progvar(nf)%numedges /= MISSING_I) then
   nxp = progvar(nf)%numedges
else
     write(string1,*) 'ERROR, ',trim(progvar(nf)%varname),' is not defined on edges or cells'
     call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
endif

nzp  = progvar(nf)%numvertical
iloc = 1 + (myindx-1) / nzp    ! cell index
vloc = myindx - (iloc-1)*nzp   ! vertical level index

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
      depth = zGridEdge(vloc,iloc)
   else
      call error_handler(E_ERR, 'get_state_meta_data', 'no support for edges at face depths', &
                         source, revision, revdate)
   endif
else
   if ( progvar(nf)%ZonHalf ) then
      depth = zGridCenter(vloc,iloc)
   else if (nzp <= 1) then
      depth = zGridFace(1,iloc)
   else
      depth = zGridFace(vloc,iloc)
   endif
endif

if (progvar(nf)%numedges /= MISSING_I) then
   if (.not. data_on_edges) then
      call error_handler(E_ERR, 'get_state_meta_data', &
                        'Internal error: numedges present but data_on_edges false', &
                        source, revision, revdate, text2='variable '//trim(progvar(nf)%varname))
   endif
   if (nzp <= 1) then
      location = set_location(lonEdge(iloc),latEdge(iloc), depth, VERTISSURFACE)
   else
      location = set_location(lonEdge(iloc),latEdge(iloc), depth, VERTISHEIGHT)
   endif
else ! must be on cell centers
   if (nzp <= 1) then
      location = set_location(lonCell(iloc),latCell(iloc), depth, VERTISSURFACE)
   else
      location = set_location(lonCell(iloc),latCell(iloc), depth, VERTISHEIGHT)
   endif
endif

! Let us return the vert location with the requested vertical localization coordinate
! hoping that the code can run faster when same points are close to many obs
! since vert_convert does not need to be called repeatedly in get_close_obs any more.
! FIXME: we should test this to see if it's a win (in computation time).  for obs
! which are not dense relative to the grid, this might be slower than doing the
! conversions on demand in the localization code (in get_close_obs()).

if ( .not. horiz_dist_only .and. vert_localization_coord /= VERTISHEIGHT ) then
     new_location = location
     call vert_convert(ens_mean, new_location, progvar(nf)%dart_kind, vert_localization_coord, istatus)
     if(istatus == 0) location = new_location
endif

if (debug > 12) then

    write(*,'("INDEX_IN / myindx / IVAR / NX, NZ: ",2(i10,2x),3(i5,2x))') index_in, myindx, nf, nxp, nzp
    write(*,'("                       ILOC, KLOC: ",2(i5,2x))') iloc, vloc
    write(*,'("                      LON/LAT/HGT: ",3(f12.3,2x))') lonCell(iloc), latCell(iloc), depth

endif

if (present(var_type)) then
   var_type = progvar(nf)%dart_kind
endif

end subroutine get_state_meta_data


!------------------------------------------------------------------

subroutine model_interpolate(x, location, obs_type, interp_val, istatus)

! given a state vector, a location, and a QTY_xxx, return the
! interpolated value at that location, and an error code.  0 is success,
! anything positive is an error.  (negative reserved for system use)
!
!       ERROR codes:
!
!       ISTATUS = 99:  general error in case something terrible goes wrong...
!       ISTATUS = 88:  this kind is not in the state vector
!       ISTATUS = 11:  Could not find a triangle that contains this lat/lon
!       ISTATUS = 12:  depth vertical coordinate out of model range.
!       ISTATUS = 13:  Missing value in interpolation.
!       ISTATUS = 16:  Don't know how to do vertical velocity for now
!       ISTATUS = 17:  Unable to compute pressure values
!       ISTATUS = 18:  altitude illegal
!       ISTATUS = 19:  could not compute u using RBF code
!       ISTATUS = 101: Internal error; reached end of subroutine without
!                      finding an applicable case.
!

! passed variables

real(r8),            intent(in)  :: x(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_type
real(r8),            intent(out) :: interp_val
integer,             intent(out) :: istatus

! local storage

integer  :: ivar, obs_kind
integer  :: tvars(3)
logical  :: goodkind
real(r8) :: values(3), lpres

if ( .not. module_initialized ) call static_init_model

interp_val = MISSING_R8
istatus    = 99           ! must be positive (and integer)

! rename for sanity - we can't change the argument names
! to this subroutine, but this really is a kind.
obs_kind = obs_type

if (debug > 0) then
   call write_location(0,location,charstring=string1)
   print *, my_task_id(), 'stt x(1), kind, loc ', x(1), obs_kind, trim(string1)
endif

! if we can interpolate any other kinds that are not directly in the
! state vector, then add those cases here.
ivar = get_progvar_index_from_kind(obs_kind)
if (ivar > 0) then
   goodkind = .true.
else
   goodkind = .false.
   ! exceptions if the kind isn't directly
   ! a field in the state vector:
   select case (obs_kind)
      !case (QTY_TEMPERATURE)  ! potential temperature is in state vector not in-situ temp
      !   goodkind = .true.
      !case (QTY_PRESSURE)
      !   goodkind = .true.
   end select
endif

! this kind is not in the state vector and it isn't one of the exceptions
! that we know how to handle.
if (.not. goodkind) then
   if (debug > 4) print *, 'kind rejected', obs_kind
   istatus = 88
   goto 100
endif

! if you add exceptions above, then you need code like this:
! if (obs_kind == QTY_xxx) then
!   add code here for how to interpolate it.
!
!   compute_scalar_with_barycentric() can take up to 3 kinds at one location
!   and return up to 3 values which can be combined or computed on here.
!
!   to use the RBF functions, here is an example call.  the third arg
!   controls whether to return the meridional (.false.) or zonal (.true.) component.
!   call compute_u_with_rbf(x, location, .TRUE., interp_val, istatus)
!
!   on error, goto 100 to get the final error checking code
!
! else

   ! direct interpolation, kind is in the state vector
   tvars(1) = ivar
   call compute_scalar_with_barycentric(x, location, 1, tvars, values, istatus)
   interp_val = values(1)
   if (debug > 4) print *, 'called generic compute_w_bary, kind, val, istatus: ', obs_kind, interp_val, istatus
   if (istatus /= 0) goto 100

!endif

100 continue

! this is for debugging - when we're confident the code is
! returning consistent values and rc codes, both these tests can
! be removed for speed.  FIXME.
if (istatus /= 0 .and. interp_val /= MISSING_R8) then
   write(string1,*) 'interp routine returned a bad status but not a MISSING_R8 value'
   write(string2,*) 'value = ', interp_val, ' istatus = ', istatus
   call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate, &
                      text2=string2)
endif
if (istatus == 0 .and. interp_val == MISSING_R8) then
   write(string1,*) 'interp routine returned a good status but set value to MISSING_R8'
   call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)
endif

if (debug > 0) then
   call write_location(0,location,charstring=string1)
   print *, my_task_id(), 'end x(1), kind, loc, val, rc ', x(1), obs_kind, trim(string1), interp_val, istatus
endif

end subroutine model_interpolate


!------------------------------------------------------------------

function nc_write_model_atts( ncFileID ) result (ierr)

! TJH -- Writes the model-specific attributes to a netCDF file.
!     This includes coordinate variables and some metadata, but NOT
!     the model state vector.
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

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

!----------------------------------------------------------------------
! variables if we just blast out one long state vector
!----------------------------------------------------------------------

integer :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)

integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

!----------------------------------------------------------------------
! variables if we parse the state vector into prognostic variables.
!----------------------------------------------------------------------

! for the dimensions and coordinate variables
integer :: nCellsDimID
integer :: nEdgesDimID, maxEdgesDimID
integer :: nVerticesDimID
integer :: VertexDegreeDimID
integer :: nVertLevelsDimID
integer :: nVertLevelsP1DimID


! for the prognostic variables
integer :: ivar, VarID, mpasFileID

!----------------------------------------------------------------------
! local variables
!----------------------------------------------------------------------

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

character(len=128) :: filename

real(r8), allocatable, dimension(:)   :: data1d

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file,
! and then put into define mode.
!-------------------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
                                   'nc_write_model_atts', 'inquire '//trim(filename))
call nc_check(nf90_Redef(ncFileID),'nc_write_model_atts',   'redef '//trim(filename))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension.
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name='copy', dimid=MemberDimID), &
                           'nc_write_model_atts', 'copy dimid '//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid=  TimeDimID), &
                           'nc_write_model_atts', 'time dimid '//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)'Time Dimension ID ',TimeDimID, &
             ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', string1, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------

call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable', len=model_size, &
        dimid = StateVarDimID),'nc_write_model_atts', 'state def_dim '//trim(filename))

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------

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
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'MPAS_OCN' ), &
           'nc_write_model_atts', 'model put '//trim(filename))

!-------------------------------------------------------------------------------
! Here is the extensible part. The simplest scenario is to output the state vector,
! parsing the state vector into model-specific parts is complicated, and you need
! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
! complicated part.
!-------------------------------------------------------------------------------

if ( output_state_vector ) then

   !----------------------------------------------------------------------------
   ! Create a variable for the state vector
   !----------------------------------------------------------------------------

   ! Define the actual (3D) state vector, which gets filled as time goes on ...
   call nc_check(nf90_def_var(ncid=ncFileID, name='state', xtype=nf90_real, &
                 dimids=(/StateVarDimID,MemberDimID,unlimitedDimID/),varid=StateVarID),&
                 'nc_write_model_atts','state def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarID,'long_name','model state or fcopy'),&
                 'nc_write_model_atts', 'state long_name '//trim(filename))

   ! Leave define mode.
   call nc_check(nf90_enddef(ncFileID),'nc_write_model_atts','state enddef '//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to output the prognostic variables.
   !----------------------------------------------------------------------------
   ! Define the new dimensions IDs
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nCells', &
          len = nCells, dimid = nCellsDimID),'nc_write_model_atts', 'nCells def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nEdges', &
          len = nEdges, dimid = nEdgesDimID),'nc_write_model_atts', 'nEdges def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='maxEdges', &
          len = maxEdges, dimid = maxEdgesDimID),'nc_write_model_atts', 'maxEdges def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nVertices', &
          len = nVertices, dimid = nVerticesDimID),'nc_write_model_atts', &
               'nVertices def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='VertexDegree', &
          len = VertexDegree, dimid = VertexDegreeDimID),'nc_write_model_atts', &
               'VertexDegree def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nVertLevels', &
          len = nVertLevels, dimid = NVertLevelsDimID),'nc_write_model_atts', &
                                      'nVertLevels def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nVertLevelsP1', &
          len = nVertLevelsP1, dimid = NVertLevelsP1DimID),'nc_write_model_atts', &
                                      'nVertLevelsP1 def_dim '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Coordinate Variables and the Attributes
   !----------------------------------------------------------------------------

   ! Cell Longitudes
   call nc_check(nf90_def_var(ncFileID,name='lonCell', xtype=nf90_double, &
                 dimids=nCellsDimID, varid=VarID),&
                 'nc_write_model_atts', 'lonCell def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'cell center longitudes'), &
                 'nc_write_model_atts', 'lonCell long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'lonCell units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'lonCell valid_range '//trim(filename))

   ! Cell Latitudes
   call nc_check(nf90_def_var(ncFileID,name='latCell', xtype=nf90_double, &
                 dimids=nCellsDimID, varid=VarID),&
                 'nc_write_model_atts', 'latCell def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'cell center latitudes'), &
                 'nc_write_model_atts', 'latCell long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_north'),  &
                 'nc_write_model_atts', 'latCell units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'latCell valid_range '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='xCell', xtype=nf90_double, &
                 dimids=nCellsDimID, varid=VarID),&
                 'nc_write_model_atts', 'xCell def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'cell center x cartesian coordinates'), &
                 'nc_write_model_atts', 'xCell long_name '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='yCell', xtype=nf90_double, &
                 dimids=nCellsDimID, varid=VarID),&
                 'nc_write_model_atts', 'yCell def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'cell center y cartesian coordinates'), &
                 'nc_write_model_atts', 'yCell long_name '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='zCell', xtype=nf90_double, &
                 dimids=nCellsDimID, varid=VarID),&
                 'nc_write_model_atts', 'zCell def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'cell center z cartesian coordinates'), &
                 'nc_write_model_atts', 'zCell long_name '//trim(filename))

   ! Grid vertical information
!  call nc_check(nf90_def_var(ncFileID,name='zgrid',xtype=nf90_double, &
!                dimids=(/ nVertLevelsP1DimID, nCellsDimID /) ,varid=VarID), &
!                'nc_write_model_atts', 'zgrid def_var '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'grid zgrid'), &
!                'nc_write_model_atts', 'zgrid long_name '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'meters'),  &
!                'nc_write_model_atts', 'zgrid units '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID,  VarID, 'positive', 'up'),  &
!                'nc_write_model_atts', 'zgrid units '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Z'),   &
!                'nc_write_model_atts', 'zgrid cartesian_axis '//trim(filename))

   ! Vertex Longitudes
   call nc_check(nf90_def_var(ncFileID,name='lonVertex', xtype=nf90_double, &
                 dimids=nVerticesDimID, varid=VarID),&
                 'nc_write_model_atts', 'lonVertex def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'vertex longitudes'), &
                 'nc_write_model_atts', 'lonVertex long_name '//trim(filename))

   ! Vertex Latitudes
   call nc_check(nf90_def_var(ncFileID,name='latVertex', xtype=nf90_double, &
                 dimids=nVerticesDimID, varid=VarID),&
                 'nc_write_model_atts', 'latVertex def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'vertex latitudes'), &
                 'nc_write_model_atts', 'latVertex long_name '//trim(filename))

   if(data_on_edges) then
      ! Edge Longitudes
      call nc_check(nf90_def_var(ncFileID,name='lonEdge', xtype=nf90_double, &
                    dimids=nEdgesDimID, varid=VarID),&
                    'nc_write_model_atts', 'lonEdge def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'edge longitudes'), &
                    'nc_write_model_atts', 'lonEdge long_name '//trim(filename))

      ! Edge Latitudes
      call nc_check(nf90_def_var(ncFileID,name='latEdge', xtype=nf90_double, &
                    dimids=nEdgesDimID, varid=VarID),&
                    'nc_write_model_atts', 'latEdge def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'edge latitudes'), &
                    'nc_write_model_atts', 'latEdge long_name '//trim(filename))
   endif

   ! Grid relationship information
   call nc_check(nf90_def_var(ncFileID,name='nEdgesOnCell',xtype=nf90_int, &
                 dimids=nCellsDimID ,varid=VarID), &
                 'nc_write_model_atts', 'nEdgesOnCell def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'grid nEdgesOnCell'), &
                 'nc_write_model_atts', 'nEdgesOnCell long_name '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='cellsOnVertex',xtype=nf90_int, &
                 dimids=(/ VertexDegreeDimID, nVerticesDimID /) ,varid=VarID), &
                 'nc_write_model_atts', 'cellsOnVertex def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'grid cellsOnVertex'), &
                 'nc_write_model_atts', 'cellsOnVertex long_name '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='verticesOnCell',xtype=nf90_int, &
                 dimids=(/ maxEdgesDimID, nCellsDimID /) ,varid=VarID), &
                 'nc_write_model_atts', 'verticesOnCell def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'grid verticesOnCell'), &
                 'nc_write_model_atts', 'verticesOnCell long_name '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='areaCell', xtype=nf90_double, &
                 dimids=nCellsDimID, varid=VarID),&
                 'nc_write_model_atts', 'areaCell def_var '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------

   do ivar=1, nfields

      varname = trim(progvar(ivar)%varname)
      string1 = trim(filename)//' '//trim(varname)

      ! match shape of the variable to the dimension IDs

      call define_var_dims(ncFileID, ivar, MemberDimID, unlimitedDimID, myndims, mydimids)

      ! define the variable and set the attributes

      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(varname), xtype=progvar(ivar)%xtype, &
                    dimids = mydimids(1:myndims), varid=VarID),&
                    'nc_write_model_atts', trim(string1)//' def_var' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', trim(progvar(ivar)%long_name)), &
           'nc_write_model_atts', trim(string1)//' put_att long_name' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'DART_kind', trim(progvar(ivar)%kind_string)), &
           'nc_write_model_atts', trim(string1)//' put_att dart_kind' )
      call nc_check(nf90_put_att(ncFileID, VarID, 'units', trim(progvar(ivar)%units)), &
           'nc_write_model_atts', trim(string1)//' put_att units' )

   enddo

   !----------------------------------------------------------------------------
   ! Finished with dimension/variable definitions, must end 'define' mode to fill.
   !----------------------------------------------------------------------------

   call nc_check(nf90_enddef(ncFileID), 'prognostic enddef '//trim(filename))

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables that DART needs and has locally
   !----------------------------------------------------------------------------

   call nc_check(NF90_inq_varid(ncFileID, 'lonCell', VarID), &
                 'nc_write_model_atts', 'lonCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, lonCell ), &
                'nc_write_model_atts', 'lonCell put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'latCell', VarID), &
                 'nc_write_model_atts', 'latCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, latCell ), &
                'nc_write_model_atts', 'latCell put_var '//trim(filename))

   if(data_on_edges) then
      call nc_check(NF90_inq_varid(ncFileID, 'lonEdge', VarID), &
                    'nc_write_model_atts', 'lonEdge inq_varid '//trim(filename))
      call nc_check(nf90_put_var(ncFileID, VarID, lonEdge ), &
                   'nc_write_model_atts', 'lonEdge put_var '//trim(filename))

      call nc_check(NF90_inq_varid(ncFileID, 'latEdge', VarID), &
                    'nc_write_model_atts', 'latEdge inq_varid '//trim(filename))
      call nc_check(nf90_put_var(ncFileID, VarID, latEdge ), &
                   'nc_write_model_atts', 'latEdge put_var '//trim(filename))
   endif

   call nc_check(NF90_inq_varid(ncFileID, 'hZLevel', VarID), &
                 'nc_write_model_atts', 'hZLevel inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, hZLevel ), &
                'nc_write_model_atts', 'hZLevel put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'zMid', VarID), &
                 'nc_write_model_atts', 'zMid inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, zMid ), &
                'nc_write_model_atts', 'zMid put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'nEdgesOnCell', VarID), &
                 'nc_write_model_atts', 'nEdgesOnCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, nEdgesOnCell ), &
                'nc_write_model_atts', 'nEdgesOnCell put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'verticesOnCell', VarID), &
                 'nc_write_model_atts', 'verticesOnCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, verticesOnCell ), &
                'nc_write_model_atts', 'verticesOnCell put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'cellsOnVertex', VarID), &
                 'nc_write_model_atts', 'cellsOnVertex inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, cellsOnVertex ), &
                'nc_write_model_atts', 'cellsOnVertex put_var '//trim(filename))

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables needed for plotting only.
   ! DART has not read these in, so we have to read them from the input file
   ! and parrot them to the DART output file.
   !----------------------------------------------------------------------------

   call nc_check(nf90_open(trim(grid_definition_filename), NF90_NOWRITE, mpasFileID), &
                 'nc_write_model_atts','open '//trim(grid_definition_filename))

   allocate(data1d(nCells))
   call nc_check(nf90_inq_varid(mpasFileID, 'xCell', VarID), &
                 'nc_write_model_atts',     'xCell inq_varid ')
   call nc_check(nf90_get_var(mpasFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'xCell get_var ')
   call nc_check(nf90_inq_varid(ncFileID,   'xCell', VarID), &
                 'nc_write_model_atts',     'xCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'xCell put_var '//trim(filename))

   call nc_check(nf90_inq_varid(mpasFileID, 'yCell', VarID), &
                 'nc_write_model_atts',     'yCell inq_varid ')
   call nc_check(nf90_get_var(mpasFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'yCell get_var ')
   call nc_check(nf90_inq_varid(ncFileID,   'yCell', VarID), &
                 'nc_write_model_atts',     'yCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'yCell put_var '//trim(filename))

   call nc_check(nf90_inq_varid(mpasFileID, 'zCell', VarID), &
                 'nc_write_model_atts',     'zCell inq_varid ')
   call nc_check(nf90_get_var(mpasFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'zCell get_var ')
   call nc_check(nf90_inq_varid(ncFileID,   'zCell', VarID), &
                 'nc_write_model_atts',     'zCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'zCell put_var '//trim(filename))

   call nc_check(nf90_inq_varid(mpasFileID, 'areaCell', VarID), &
                 'nc_write_model_atts',     'areaCell inq_varid ')
   call nc_check(nf90_get_var(mpasFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'areaCell get_var ')
   call nc_check(nf90_inq_varid(ncFileID,   'areaCell', VarID), &
                 'nc_write_model_atts',     'areaCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'areaCell put_var '//trim(filename))
   deallocate(data1d)

   allocate(data1d(nVertices))
   call nc_check(nf90_inq_varid(mpasFileID, 'lonVertex', VarID), &
                 'nc_write_model_atts',     'lonVertex inq_varid ')
   call nc_check(nf90_get_var(mpasFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'lonVertex get_var ')
   call nc_check(nf90_inq_varid(ncFileID,   'lonVertex', VarID), &
                 'nc_write_model_atts',     'lonVertex inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'lonVertex put_var '//trim(filename))

   call nc_check(nf90_inq_varid(mpasFileID, 'latVertex', VarID), &
                 'nc_write_model_atts',     'latVertex inq_varid ')
   call nc_check(nf90_get_var(mpasFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'latVertex get_var ')
   call nc_check(nf90_inq_varid(ncFileID,   'latVertex', VarID), &
                 'nc_write_model_atts',     'latVertex inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'latVertex put_var '//trim(filename))
   deallocate(data1d)

   call nc_check(nf90_close(mpasFileID),'nc_write_model_atts','close '//trim(grid_definition_filename))
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts


!------------------------------------------------------------------

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

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: state_vec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME)          :: varname
integer :: i, ivar, VarID, ncNdims, dimlen
integer :: TimeDimID, CopyDimID

real(r8), allocatable, dimension(:)       :: data_1d_array
real(r8), allocatable, dimension(:,:)     :: data_2d_array
real(r8), allocatable, dimension(:,:,:)   :: data_3d_array

character(len=128) :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file,
!-------------------------------------------------------------------------------

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

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------

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
            call error_handler(E_ERR, 'nc_write_model_vars', trim(string1), &
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

      if ((debug > 9) .and. do_output()) then
         write(*,*)'nc_write_model_vars '//trim(varname)//' start is ',mystart(1:ncNdims)
         write(*,*)'nc_write_model_vars '//trim(varname)//' count is ',mycount(1:ncNdims)
      endif

      if (     progvar(ivar)%numdims == 1 ) then

         if ( ncNdims /= 3 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 3 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
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
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
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
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
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

         ! FIXME put an error message here
         write(string1,*)'no support (yet) for 4d fields'
         call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate)

      endif

   enddo


endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), 'nc_write_model_vars', 'sync '//trim(filename))

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars


!------------------------------------------------------------------

function get_model_size()

! Returns the size of the model as an integer.
! Required for all applications.

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size


!------------------------------------------------------------------

function get_model_time_step()

! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = model_timestep

end function get_model_time_step


!------------------------------------------------------------------

subroutine ens_mean_for_model(filter_ens_mean)

! If needed by the model interface, this is the current mean
! for all state vector items across all ensembles.

real(r8), intent(in) :: filter_ens_mean(:)

if ( .not. module_initialized ) call static_init_model

ens_mean = filter_ens_mean

if ((debug > 3) .and. do_output()) then
   if (do_output()) print *, 'resetting ensemble mean: '
   call print_variable_ranges(ens_mean)
endif

end subroutine ens_mean_for_model


!------------------------------------------------------------------

subroutine end_model()

! Does any shutdown and clean-up needed for model.

if (allocated(latCell))        deallocate(latCell)
if (allocated(lonCell))        deallocate(lonCell)
if (allocated(zGridFace))      deallocate(zGridFace)
if (allocated(zGridCenter))    deallocate(zGridCenter)
if (allocated(cellsOnVertex))  deallocate(cellsOnVertex)
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

call finalize_closest_center()

end subroutine end_model


!------------------------------------------------------------------

subroutine pert_model_state(state, pert_state, interf_provided)

! Perturbs a model state for generating initial ensembles.
! The perturbed state is returned in pert_state.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding a perturbation to each model state
! variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

real(r8)              :: pert_ampl
real(r8)              :: minv, maxv, temp
type(random_seq_type) :: random_seq
integer               :: i, j, s, e
integer, save         :: counter = 1


! generally you do not want to perturb a single state
! to begin an experiment - unless you make minor perturbations
! and then run the model free for long enough that differences
! develop which contain actual structure.
!
! the subsequent code is a pert routine which
! can be used to add minor perturbations which can be spun up.
!
! if all values in a field are identical (i.e. 0.0) this
! routine will not change those values since it won't
! make a new value outside the original min/max of that
! variable in the state vector.  to handle this case you can
! remove the min/max limit lines below.


! start of pert code

if ( .not. module_initialized ) call static_init_model

interf_provided = .true.

! the first time through get the task id (0:N-1)
! and set a unique seed per task.  this won't
! be consistent between different numbers of mpi
! tasks, but at least it will reproduce with
! multiple runs with the same task count.
! best i can do since this routine doesn't have
! the ensemble member number as an argument
! (which i think it needs for consistent seeds).
!
! this only executes the first time since counter
! gets incremented after the first use and the value
! is saved between calls.
if (counter == 1) counter = counter + (my_task_id() * 1000)

call init_random_seq(random_seq, counter)
counter = counter + 1

do i=1, nfields
   ! starting and ending indices in the linear state vect
   ! for each different state kind.
   s = progvar(i)%index1
   e = progvar(i)%indexN
   ! original min/max data values of each type
   minv = minval(state(s:e))
   maxv = maxval(state(s:e))
   do j=s, e
      ! once you change pert_state, state is changed as well
      ! since they are the same storage as called from filter.
      ! you have to save it if you want to use it again.
      temp = state(j)  ! original value
      ! perturb each value individually
      ! make the perturbation amplitude N% of this value
      pert_ampl = model_perturbation_amplitude * temp
      pert_state(j) = random_gaussian(random_seq, state(j), pert_ampl)
      ! keep it from exceeding the original range
      pert_state(j) = max(minv, pert_state(j))
      pert_state(j) = min(maxv, pert_state(j))
   enddo
enddo

end subroutine pert_model_state


!------------------------------------------------------------------

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, &
                         obs_loc, obs_kind, num_close, close_ind, dist)
!
! Given a DART location (referred to as "base") and a set of candidate
! locations & kinds (obs, obs_kind), returns the subset close to the
! "base", their indices, and their distances to the "base" ...

! For vertical distance computations, general philosophy is to convert all
! vertical coordinates to a common coordinate. This coordinate type can be
! defined in the namelist with the variable "vert_localization_coord".
! But we first try a single coordinate type as the model level here.
! FIXME: We need to add more options later.

! Vertical conversion is carried out by the subroutine vert_convert.

! Note that both base_obs_loc and obs_loc are intent(inout), meaning that these
! locations are possibly modified here and returned as such to the calling routine.
! The calling routine is always filter_assim and these arrays are local arrays
! within filter_assim. In other words, these modifications will only matter within
! filter_assim, but will not propagate backwards to filter.

type(get_close_type),              intent(in)    :: gc
type(location_type),               intent(inout) :: base_obs_loc
integer,                           intent(in)    :: base_obs_kind
type(location_type), dimension(:), intent(inout) :: obs_loc
integer,             dimension(:), intent(in)    :: obs_kind
integer,                           intent(out)   :: num_close
integer,             dimension(:), intent(out)   :: close_ind
real(r8),            dimension(:), intent(out)   :: dist

integer                :: ztypeout
integer                :: t_ind, istatus1, istatus2, k
integer                :: base_which, local_obs_which
real(r8), dimension(3) :: base_llv, local_obs_llv   ! lon/lat/vert
type(location_type)    :: local_obs_loc

real(r8) ::  hor_dist
hor_dist = 1.0e9_r8

! Initialize variables to missing status

num_close = 0
close_ind = -99
dist      = 1.0e9_r8   !something big and positive (far away) in radians
istatus1  = 0
istatus2  = 0

! Convert base_obs vertical coordinate to requested vertical coordinate if necessary

base_llv = get_location(base_obs_loc)
base_which = nint(query_location(base_obs_loc))

ztypeout = vert_localization_coord

if (.not. horiz_dist_only) then
  if (base_llv(3) == MISSING_R8) then
     istatus1 = 1
  else if (base_which /= vert_localization_coord) then
      call vert_convert(ens_mean, base_obs_loc, base_obs_kind, ztypeout, istatus1)
      if(debug > 5) then
      call write_location(0,base_obs_loc,charstring=string1)
      call error_handler(E_MSG, 'get_close_obs: base_obs_loc',string1,source, revision, revdate)
  endif
endif
endif

if (istatus1 == 0) then

   ! Loop over potentially close subset of obs priors or state variables
   ! This way, we are decreasing the number of distance computations that will follow.
   ! This is a horizontal-distance operation and we don't need to have the relevant vertical
   ! coordinate information yet (for obs_loc).
   call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                          num_close, close_ind)

   do k = 1, num_close

      t_ind = close_ind(k)
      local_obs_loc   = obs_loc(t_ind)
      local_obs_which = nint(query_location(local_obs_loc))

      ! Convert local_obs vertical coordinate to requested vertical coordinate if necessary.
      ! This should only be necessary for obs priors, as state location information already
      ! contains the correct vertical coordinate (filter_assim's call to get_state_meta_data).
      if (.not. horiz_dist_only) then
          if (local_obs_which /= vert_localization_coord) then
              call vert_convert(ens_mean, local_obs_loc, obs_kind(t_ind), ztypeout, istatus2)
          else
              istatus2 = 0
          endif
      endif

      ! Compute distance - set distance to a very large value if vert coordinate is missing
      ! or vert_interpolate returned error (istatus2=1)
      local_obs_llv = get_location(local_obs_loc)
      if (( (.not. horiz_dist_only)           .and. &
            (local_obs_llv(3) == MISSING_R8)) .or.  &
            (istatus2 /= 0)                   ) then
            dist(k) = 1.0e9_r8
      else
            dist(k) = get_dist(base_obs_loc, local_obs_loc, base_obs_kind, obs_kind(t_ind))
!       if ((debug > 4) .and. (k < 100) .and. do_output()) then
!           print *, 'calling get_dist'
!           call write_location(0,base_obs_loc,charstring=string2)
!           call error_handler(E_MSG, 'get_close_obs: base_obs_loc',string2,source, revision, revdate)
!           call write_location(0,local_obs_loc,charstring=string2)
!           call error_handler(E_MSG, 'get_close_obs: local_obs_loc',string2,source, revision, revdate)
!           hor_dist = get_dist(base_obs_loc, local_obs_loc, base_obs_kind, obs_kind(t_ind), no_vert=.true.)
!           print *, 'hor/3d_dist for k =', k, ' is ', hor_dist,dist(k)
!       endif
      endif

   enddo
endif

if ((debug > 2) .and. do_output()) then
   call write_location(0,base_obs_loc,charstring=string2)
   print *, 'get_close_obs: nclose, base_obs_loc ', num_close, trim(string2)
endif

end subroutine get_close_obs


!------------------------------------------------------------------

subroutine init_time(time)

! Companion interface to init_conditions. Returns a time that is somehow
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

! this shuts up the compiler warnings about unused variables
time = set_time(0, 0)

write(string1,*) 'Cannot initialize MPAS time via subroutine call; start_from_restart cannot be F'
call error_handler(E_ERR,'init_time',string1,source,revision,revdate)

end subroutine init_time


!------------------------------------------------------------------

subroutine init_conditions(x)

! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

! this shuts up the compiler warnings about unused variables
x = 0.0_r8

write(string1,*) 'Cannot initialize MPAS state via subroutine call; start_from_restart cannot be F'
call error_handler(E_ERR,'init_conditions',string1,source,revision,revdate)

end subroutine init_conditions


!------------------------------------------------------------------

subroutine adv_1step(x, time)

! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to
! compute a timestep, for instance for radiation computations.
! This interface is only called IF the namelist parameter
! async is set to 0 in perfect_model_obs or filter -OR- if the
! program integrate_model is to be used to advance the model
! state as a separate executable. If none of these options
! are used (the model will only be advanced as a separate
! model-specific executable), this can be a NULL INTERFACE.

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

if ( .not. module_initialized ) call static_init_model

if (do_output()) then
   call print_time(time,'NULL interface adv_1step (no advance) DART time is')
   call print_time(time,'NULL interface adv_1step (no advance) DART time is',logfileunit)
endif

write(string1,*) 'Cannot advance MPAS with a subroutine call; async cannot equal 0'
call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate)

end subroutine adv_1step



!==================================================================
! The (model-specific) additional public interfaces come next
!  (these are not required by dart but are used by other programs)
!==================================================================


subroutine get_model_analysis_filename( filename )

! return the name of the analysis filename that was set
! in the model_nml namelist

character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(model_analysis_filename)

end subroutine get_model_analysis_filename


!-------------------------------------------------------------------

subroutine get_grid_definition_filename( filename )

! return the name of the grid_definition filename that was set
! in the model_nml namelist

character(len=*), intent(out) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(grid_definition_filename)

end subroutine get_grid_definition_filename


!-------------------------------------------------------------------

subroutine analysis_file_to_statevector(filename, state_vector, model_time)

! Reads the current time and state variables from a mpas analysis
! file and packs them into a dart state vector.

character(len=*), intent(in)    :: filename
real(r8),         intent(inout) :: state_vector(:)
type(time_type),  intent(out)   :: model_time

! temp space to hold data while we are reading it
integer  :: ndim1, ndim2, ndim3
integer  :: i, ivar
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array
real(r8), allocatable, dimension(:,:,:)     :: data_3d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME) :: varname
integer :: VarID, ncNdims, dimlen
integer :: ncid, TimeDimID, TimeDimLength
character(len=256) :: myerrorstring

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

! Check that the input file exists ...

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'analysis_file_to_statevector',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncid), &
             'analysis_file_to_statevector','open '//trim(filename))

model_time = get_analysis_time(ncid, filename)

! let the calling program print out the time information it wants.
!if (do_output()) &
!    call print_time(model_time,'time in restart file '//trim(filename))
!if (do_output()) &
!    call print_date(model_time,'date in restart file '//trim(filename))

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

   varname = trim(progvar(ivar)%varname)
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

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(myerrorstring),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         call error_handler(E_ERR,'analysis_file_to_statevector',string1,source,revision,revdate)
      endif

      mycount(i) = dimlen

   enddo DimCheck

   where(dimIDs == TimeDimID) mystart = TimeDimLength  ! pick the latest time
   where(dimIDs == TimeDimID) mycount = 1              ! only use one time

   if ((debug > 10) .and. do_output()) then
      write(*,*)'analysis_file_to_statevector '//trim(varname)//' start = ',mystart(1:ncNdims)
      write(*,*)'analysis_file_to_statevector '//trim(varname)//' count = ',mycount(1:ncNdims)
   endif

   if (ncNdims == 1) then

      ! If the single dimension is TIME, we only need a scalar.
      ! Pretty sure this cannot happen ...
      ndim1 = mycount(1)
      allocate(data_1d_array(ndim1))
      call nc_check(nf90_get_var(ncid, VarID, data_1d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'analysis_file_to_statevector', 'get_var '//trim(varname))

      call prog_var_to_vector(data_1d_array, state_vector, ivar)
      deallocate(data_1d_array)

   elseif (ncNdims == 2) then

      ndim1 = mycount(1)
      ndim2 = mycount(2)
      allocate(data_2d_array(ndim1, ndim2))
      call nc_check(nf90_get_var(ncid, VarID, data_2d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'analysis_file_to_statevector', 'get_var '//trim(varname))

      call prog_var_to_vector(data_2d_array, state_vector, ivar)
      deallocate(data_2d_array)

   elseif (ncNdims == 3) then

      ndim1 = mycount(1)
      ndim2 = mycount(2)
      ndim3 = mycount(3)
      allocate(data_3d_array(ndim1, ndim2, ndim3))
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

call nc_check(nf90_close(ncid), &
             'analysis_file_to_statevector','close '//trim(filename))

end subroutine analysis_file_to_statevector


!-------------------------------------------------------------------

subroutine statevector_to_analysis_file(state_vector, filename, statetime)

! Writes the current time and state variables from a dart state
! vector (1d array) into a mpas netcdf analysis file.

real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: filename
type(time_type),  intent(in) :: statetime

! temp space to hold data while we are writing it
integer :: i, ivar
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array
real(r8), allocatable, dimension(:,:,:)     :: data_3d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME) :: varname
integer :: VarID, ncNdims, dimlen
integer :: ncFileID, TimeDimID, TimeDimLength
logical :: done_winds
type(time_type) :: model_time

if ( .not. module_initialized ) call static_init_model

! Check that the output file exists ...

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for writing.'
   call error_handler(E_ERR,'statevector_to_analysis_file',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_WRITE, ncFileID), &
             'statevector_to_analysis_file','open '//trim(filename))

! make sure the time in the file is the same as the time on the data
! we are trying to insert.  we are only updating part of the contents
! of the mpas analysis file, and state vector contents from a different
! time won't be consistent with the rest of the file.

model_time = get_analysis_time(ncFileID, filename)

if ( model_time /= statetime ) then
   call print_time( statetime,'DART current time',logfileunit)
   call print_time(model_time,'mpas current time',logfileunit)
   call print_time( statetime,'DART current time')
   call print_time(model_time,'mpas current time')
   write(string1,*)trim(filename),' current time must equal model time'
   call error_handler(E_ERR,'statevector_to_analysis_file',string1,source,revision,revdate)
endif

! let the calling program print out the time information it wants.
!if (do_output()) &
!    call print_time(statetime,'time of DART file '//trim(filename))
!if (do_output()) &
!    call print_date(statetime,'date of DART file '//trim(filename))

! The DART prognostic variables are only defined for a single time.
! We already checked the assumption that variables are xy2d or xyz3d ...
! IF the netCDF variable has a TIME dimension, it must be the last dimension,
! and we need to read the LAST timestep and effectively squeeze out the
! singleton dimension when we stuff it into the DART state vector.

TimeDimID = FindTimeDimension( ncFileID )

if ( TimeDimID > 0 ) then
   call nc_check(nf90_inquire_dimension(ncFileID, TimeDimID, len=TimeDimLength), &
            'statevector_to_analysis_file', 'inquire timedimlength '//trim(filename))
else
   TimeDimLength = 0
endif

done_winds = .false.
PROGVARLOOP : do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   string2 = trim(filename)//' '//trim(varname)

   if (( varname == 'uReconstructZonal' .or. &
         varname == 'uReconstructMeridional' ) .and. update_u_from_reconstruct ) then
      if (done_winds) cycle PROGVARLOOP

      ! this routine updates the edge winds from both the zonal and meridional
      ! fields, so only call it once.
      call update_wind_components(ncFileID, state_vector, use_increments_for_u_update)
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

   call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'statevector_to_analysis_file', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
            'statevector_to_analysis_file', 'inquire '//trim(string2))

   mystart = 1   ! These are arrays, actually.
   mycount = 1
   DimCheck : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
            'statevector_to_analysis_file', string1)

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string2),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         write(string2,*)' but it should be.'
         call error_handler(E_ERR, 'statevector_to_analysis_file', string1, &
                         source, revision, revdate, text2=string2)
      endif

      mycount(i) = dimlen

   enddo DimCheck


   where(dimIDs == TimeDimID) mystart = TimeDimLength
   where(dimIDs == TimeDimID) mycount = 1   ! only the latest one

   if ((debug > 9) .and. do_output()) then
      write(*,*)'statevector_to_analysis_file '//trim(varname)//' start is ',mystart(1:ncNdims)
      write(*,*)'statevector_to_analysis_file '//trim(varname)//' count is ',mycount(1:ncNdims)
   endif


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

      call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array, &
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

      call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
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

      call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
            start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'statevector_to_analysis_file', 'put_var '//trim(varname))
      deallocate(data_3d_array)

   else
      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'statevector_to_analysis_file', string1, &
                        source,revision,revdate)
   endif

enddo PROGVARLOOP

call nc_check(nf90_close(ncFileID), &
             'statevector_to_analysis_file','close '//trim(filename))

end subroutine statevector_to_analysis_file


!------------------------------------------------------------------

subroutine do_clamping(out_of_range_fail, range, dimsize, varname, array_1d, array_2d, array_3d)
 logical,          intent(in)    :: out_of_range_fail
 real(r8),         intent(in)    :: range(2)
 integer,          intent(in)    :: dimsize
 character(len=*), intent(in)    :: varname
 real(r8),optional,intent(inout) :: array_1d(:), array_2d(:,:), array_3d(:,:,:)

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

if ((debug > 6) .and. do_output()) then
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

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'get_analysis_time',string1,source,revision,revdate)
endif

! find the first number and use that as the start of the string conversion
i = scan(filename, "0123456789")
if (i <= 0) then
   write(string1,*) 'cannot find time string in name ', trim(filename)
   call error_handler(E_ERR,'get_analysis_time',string1,source,revision,revdate)
endif

get_analysis_time_fname = string_to_time(filename(i:i+19))

end function get_analysis_time_fname


!------------------------------------------------------------------

subroutine write_model_time(time_filename, model_time, adv_to_time)
 character(len=*), intent(in)           :: time_filename
 type(time_type),  intent(in)           :: model_time
 type(time_type),  intent(in), optional :: adv_to_time

integer :: iunit
character(len=19) :: timestring
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

end subroutine write_model_time


!------------------------------------------------------------------

subroutine get_grid_dims(Cells, Vertices, Edges, VertLevels, VertexDeg)

! public routine for returning the counts of various things in the grid
!

integer, intent(out) :: Cells         ! Total number of cells making up the grid
integer, intent(out) :: Vertices      ! Unique points in grid which are corners of cells
integer, intent(out) :: Edges         ! Straight lines between vertices making up cells
integer, intent(out) :: VertLevels    ! Vertical levels; count of vert cell centers
integer, intent(out) :: VertexDeg     ! Max number of edges that touch any vertex

if ( .not. module_initialized ) call static_init_model

Cells      = nCells
Vertices   = nVertices
Edges      = nEdges
VertLevels = nVertLevels
VertexDeg  = vertexDegree

end subroutine get_grid_dims


!==================================================================
! The (model-specific) private interfaces come last
!==================================================================


!------------------------------------------------------------------

function time_to_string(t, interval)

! convert time type into a character string with the
! format of YYYY-MM-DD_hh:mm:ss

! passed variables
 character(len=19) :: time_to_string
 type(time_type), intent(in) :: t
 logical, intent(in), optional :: interval

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
!
! TJH Seems like the MPAS OCN developers start their runs from year 1.
! This just doesn't work for DA - but I need to exercise the code,
! so I am modifying it. FIXME ... should be a hard error.

type(time_type) :: string_to_time
character(len=*), intent(in) :: s

integer :: iyear, imonth, iday, ihour, imin, isec

read( s ,'(i4,5(1x,i2))') iyear, imonth, iday, ihour, imin, isec

if (iyear < 1601) then
   write(string1,*)'WARNING: Converting YEAR ',iyear,' to ',iyear+1601
   write(string2,*)'original time (string) is <',trim(s),'>'
   call error_handler(E_MSG, 'string_to_time', string1, &
               source, revision, revdate, text2=string2)
   iyear = iyear + 1601
endif

string_to_time = set_date(iyear, imonth, iday, ihour, imin, isec)

end function string_to_time


!------------------------------------------------------------------

function set_model_time_step()

! the static_init_model ensures that the model namelists are read.

type(time_type) :: set_model_time_step

if ( .not. module_initialized ) call static_init_model

! these are from the namelist
!FIXME: sanity check these for valid ranges?
set_model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end function set_model_time_step


!------------------------------------------------------------------

subroutine read_grid_dims()

! Read the grid dimensions from the MPAS netcdf file.
!
! The file name comes from module storage ... namelist.

integer :: grid_id, dimid

if ( .not. module_initialized ) call static_init_model

! get the ball rolling ...

call nc_check( nf90_open(trim(grid_definition_filename), NF90_NOWRITE, grid_id), &
              'read_grid_dims', 'open '//trim(grid_definition_filename))

! nCells : get dimid for 'nCells' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'nCells', dimid), &
              'read_grid_dims','inq_dimid nCells '//trim(grid_definition_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nCells), &
            'read_grid_dims','inquire_dimension nCells '//trim(grid_definition_filename))

! nVertices : get dimid for 'nVertices' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'nVertices', dimid), &
              'read_grid_dims','inq_dimid nVertices '//trim(grid_definition_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nVertices), &
            'read_grid_dims','inquire_dimension nVertices '//trim(grid_definition_filename))

! nEdges : get dimid for 'nEdges' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'nEdges', dimid), &
              'read_grid_dims','inq_dimid nEdges '//trim(grid_definition_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nEdges), &
            'read_grid_dims','inquire_dimension nEdges '//trim(grid_definition_filename))

! maxEdges : get dimid for 'maxEdges' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'maxEdges', dimid), &
              'read_grid_dims','inq_dimid maxEdges '//trim(grid_definition_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=maxEdges), &
            'read_grid_dims','inquire_dimension maxEdges '//trim(grid_definition_filename))

! nVertLevels : get dimid for 'nVertLevels' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'nVertLevels', dimid), &
              'read_grid_dims','inq_dimid nVertLevels '//trim(grid_definition_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nVertLevels), &
            'read_grid_dims','inquire_dimension nVertLevels '//trim(grid_definition_filename))

! nVertLevelsP1 : get dimid for 'nVertLevelsP1' and then get value
! This exists in the atmosphere, but not in the ocean - at this point.
! In the hopes we can consolidate this function to be useful in both instances,
! I am extending this to work in both scenarios.

if( nf90_inq_dimid(grid_id, 'nVertLevelsP1', dimid) == NF90_NOERR ) then
   call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nVertLevelsP1), &
            'read_grid_dims','inquire_dimension nVertLevelsP1 '//trim(grid_definition_filename))
else
   nVertLevelsP1 = nVertLevels + 1
endif

! vertexDegree : get dimid for 'vertexDegree' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'vertexDegree', dimid), &
              'read_grid_dims','inq_dimid vertexDegree '//trim(grid_definition_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=vertexDegree), &
            'read_grid_dims','inquire_dimension vertexDegree '//trim(grid_definition_filename))

! tidy up

call nc_check(nf90_close(grid_id), &
         'read_grid_dims','close '//trim(grid_definition_filename) )

if (debug > 4) then
   write(*,*)
   write(*,*)'read_grid_dims: nCells        is ', nCells
   write(*,*)'read_grid_dims: nVertices     is ', nVertices
   write(*,*)'read_grid_dims: nEdges        is ', nEdges
   write(*,*)'read_grid_dims: maxEdges      is ', maxEdges
   write(*,*)'read_grid_dims: nVertLevels   is ', nVertLevels
   write(*,*)'read_grid_dims: nVertLevelsP1 is ', nVertLevelsP1
   write(*,*)'read_grid_dims: vertexDegree  is ', vertexDegree
endif

end subroutine read_grid_dims


!------------------------------------------------------------------

subroutine get_grid()

! Read the actual grid values in from the MPAS netcdf file.
!
! The file name comes from module storage ... namelist.
! This reads in the following arrays:
!   latCell, lonCell, zGridFace, cellsOnVertex (all in module global storage)


integer  :: ncid, VarID, numdims, dimlen, i

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: dimname

! Read the netcdf file data

call nc_check(nf90_open(trim(grid_definition_filename), nf90_nowrite, ncid), 'get_grid', 'open '//trim(grid_definition_filename))

! Read the variables

call nc_check(nf90_inq_varid(ncid, 'latCell', VarID), &
      'get_grid', 'inq_varid latCell '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, latCell), &
      'get_grid', 'get_var latCell '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'lonCell', VarID), &
      'get_grid', 'inq_varid lonCell '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, lonCell), &
      'get_grid', 'get_var lonCell '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'hZLevel', VarID), &
      'get_grid', 'inq_varid hZLevel '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, hZLevel), &
      'get_grid', 'get_var hZLevel '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'zMid', VarID), &
      'get_grid', 'inq_varid zMid '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, zMid), &
      'get_grid', 'get_var zMid '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'cellsOnVertex', VarID), &
      'get_grid', 'inq_varid cellsOnVertex '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, cellsOnVertex), &
      'get_grid', 'get_var cellsOnVertex '//trim(grid_definition_filename))

! MPAS analysis files are in radians - at this point DART needs degrees.

latCell = latCell * rad2deg
lonCell = lonCell * rad2deg

! Read the variables

! FIXME ...   inq_varid edgeNormalVectors ../data/mpas_out.nc: NetCDF: Variable not found
edgeNormalVectors(:,:) = MISSING_R8
! call nc_check(nf90_inq_varid(ncid, 'edgeNormalVectors', VarID), &
!       'get_grid', 'inq_varid edgeNormalVectors '//trim(grid_definition_filename))
! call nc_check(nf90_get_var( ncid, VarID, edgeNormalVectors), &
!       'get_grid', 'get_var edgeNormalVectors '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'nEdgesOnCell', VarID), &
      'get_grid', 'inq_varid nEdgesOnCell '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, nEdgesOnCell), &
      'get_grid', 'get_var nEdgesOnCell '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'edgesOnCell', VarID), &
      'get_grid', 'inq_varid edgesOnCell '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, edgesOnCell), &
      'get_grid', 'get_var edgesOnCell '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'cellsOnEdge', VarID), &
      'get_grid', 'inq_varid cellsOnEdge '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, cellsOnEdge), &
      'get_grid', 'get_var cellsOnEdge '//trim(grid_definition_filename))

if(data_on_edges) then
   call nc_check(nf90_inq_varid(ncid, 'latEdge', VarID), &
         'get_grid', 'inq_varid latEdge '//trim(grid_definition_filename))
   call nc_check(nf90_get_var( ncid, VarID, latEdge), &
         'get_grid', 'get_var latEdge '//trim(grid_definition_filename))

   call nc_check(nf90_inq_varid(ncid, 'lonEdge', VarID), &
         'get_grid', 'inq_varid lonEdge '//trim(grid_definition_filename))
   call nc_check(nf90_get_var( ncid, VarID, lonEdge), &
         'get_grid', 'get_var lonEdge '//trim(grid_definition_filename))

   latEdge = latEdge * rad2deg
   lonEdge = lonEdge * rad2deg

   call nc_check(nf90_inq_varid(ncid, 'xEdge', VarID), &
         'get_grid', 'inq_varid xEdge '//trim(grid_definition_filename))
   call nc_check(nf90_get_var( ncid, VarID, xEdge), &
         'get_grid', 'get_var xEdge '//trim(grid_definition_filename))

   call nc_check(nf90_inq_varid(ncid, 'yEdge', VarID), &
         'get_grid', 'inq_varid yEdge '//trim(grid_definition_filename))
   call nc_check(nf90_get_var( ncid, VarID, yEdge), &
         'get_grid', 'get_var yEdge '//trim(grid_definition_filename))

   call nc_check(nf90_inq_varid(ncid, 'zEdge', VarID), &
         'get_grid', 'inq_varid zEdge '//trim(grid_definition_filename))
   call nc_check(nf90_get_var( ncid, VarID, zEdge), &
         'get_grid', 'get_var zEdge '//trim(grid_definition_filename))
endif

call nc_check(nf90_inq_varid(ncid, 'xVertex', VarID), &
      'get_grid', 'inq_varid xVertex '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, xVertex), &
      'get_grid', 'get_var xVertex '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'yVertex', VarID), &
      'get_grid', 'inq_varid yVertex '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, yVertex), &
      'get_grid', 'get_var yVertex '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'zVertex', VarID), &
      'get_grid', 'inq_varid zVertex '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, zVertex), &
      'get_grid', 'get_var zVertex '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'verticesOnCell', VarID), &
      'get_grid', 'inq_varid verticesOnCell '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, verticesOnCell), &
      'get_grid', 'get_var verticesOnCell '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'maxLevelEdgeTop', VarID), &
      'get_grid', 'inq_varid maxLevelEdgeTop '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, maxLevelEdgeTop), &
      'get_grid', 'get_var maxLevelEdgeTop '//trim(grid_definition_filename))

! Get the boundary information if available.
! Assuming the existence of this variable is sufficient to determine if
! the grid is defined everywhere or not.

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
              'get_grid', 'inquire boundaryEdge'//trim(model_analysis_filename))

!do i=1, numdims
!   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen, name=dimname),'get_grid')
!   write(*,*)'boundaryEdge inquire length for dimension (',i,') is ',dimlen,' '//trim(dimname)
!enddo

if ( nf90_inq_varid(ncid, 'boundaryEdge', VarID) == NF90_NOERR ) then
   allocate(boundaryEdge(nVertLevels,nEdges))
   call nc_check(nf90_get_var( ncid, VarID, boundaryEdge), &
      'get_grid', 'get_var boundaryEdge '//trim(grid_definition_filename))
   boundaryEdge(:,:) = MISSING_R8
   global_grid = .false.
endif

if ( nf90_inq_varid(ncid, 'boundaryVertex', VarID) == NF90_NOERR ) then
   allocate(boundaryVertex(nVertLevels,nVertices))
   call nc_check(nf90_get_var( ncid, VarID, boundaryVertex), &
      'get_grid', 'get_var boundaryVertex '//trim(grid_definition_filename))
   boundaryVertex(:,:) = MISSING_R8
   global_grid = .false.
endif

if ( nf90_inq_varid(ncid, 'maxLevelCell', VarID) == NF90_NOERR ) then
   allocate(maxLevelCell(nCells))
   call nc_check(nf90_get_var( ncid, VarID, maxLevelCell), &
      'get_grid', 'get_var maxLevelCell '//trim(grid_definition_filename))
   all_levels_exist_everywhere = .false.
endif

call nc_check(nf90_close(ncid), 'get_grid','close '//trim(grid_definition_filename) )

! FIXME : Use layer thicknesses (perhaps hZLevel) to determine all layer information

zGridFace(:,:)   = MISSING_R8
zGridCenter(:,:) = -1.0_r8 * zMid(:,:)    ! all depths are negative in the input file; obs are +
zGridEdge(:,:)   = MISSING_R8

! A little sanity check
if ((debug > 9) .and. do_output()) then

   write(*,*)
   write(*,*)'latCell           range ',minval(latCell),           maxval(latCell)
   write(*,*)'lonCell           range ',minval(lonCell),           maxval(lonCell)
   write(*,*)'hZLevel           range ',minval(hZLevel),           maxval(hZLevel)
   write(*,*)'cellsOnVertex     range ',minval(cellsOnVertex),     maxval(cellsOnVertex)
   write(*,*)'edgeNormalVectors range ',minval(edgeNormalVectors), maxval(edgeNormalVectors)
   write(*,*)'nEdgesOnCell      range ',minval(nEdgesOnCell),      maxval(nEdgesOnCell)
   write(*,*)'EdgesOnCell       range ',minval(EdgesOnCell),       maxval(EdgesOnCell)
   write(*,*)'cellsOnEdge       range ',minval(cellsOnEdge),       maxval(cellsOnEdge)
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
   write(*,*)'maxLevelEdgeTop   range ',minval(maxLevelEdgeTop),   maxval(maxLevelEdgeTop)
   if (allocated(boundaryEdge)) &
   write(*,*)'boundaryEdge      range ',minval(boundaryEdge),      maxval(boundaryEdge)
   if (allocated(boundaryVertex)) &
   write(*,*)'boundaryVertex    range ',minval(boundaryVertex),    maxval(boundaryVertex)
   if (allocated(maxLevelCell)) &
   write(*,*)'maxLevelCell      range ',minval(maxLevelCell),      maxval(maxLevelCell)

endif

end subroutine get_grid


!------------------------------------------------------------------

subroutine update_wind_components(ncid, state_vector, use_increments_for_u_update)

 integer,  intent(in)  :: ncid                  ! netCDF handle for model_analysis_filename
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
   call read_2d_from_nc_file(ncid, 'u', u)
   call read_2d_from_nc_file(ncid, 'uReconstructZonal', ucell)
   call read_2d_from_nc_file(ncid, 'uReconstructMeridional', vcell)

   if ((debug > 8) .and. do_output()) then
      write(*,*)
      write(*,*)'update_winds: org u          range ',minval(u),     maxval(u)
      write(*,*)'update_winds: org zonal      range ',minval(ucell), maxval(ucell)
      write(*,*)'update_winds: org meridional range ',minval(vcell), maxval(vcell)
   endif

   ! compute the increments compared to the updated values in the state vector

   allocate(data_2d_array(nVertLevels, nCells))

   call vector_to_prog_var(state_vector, zonal, data_2d_array)
   ucell = data_2d_array - ucell
   call vector_to_prog_var(state_vector, meridional, data_2d_array)
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

   call uv_cell_to_edges(ucell, vcell, u)

   if ((debug > 8) .and. do_output()) then
      write(*,*)
      write(*,*)'update_winds: u values    range ',minval(ucell), maxval(ucell)
      write(*,*)'update_winds: v values    range ',minval(vcell), maxval(vcell)
   endif

endif

if ((debug > 8) .and. do_output()) then
   write(*,*)
   write(*,*)'update_winds: u after update:',minval(u), maxval(u)
endif

! Write back to the mpas analysis file.

call put_u(ncid, u)

deallocate(ucell, vcell, u)

end subroutine update_wind_components


!------------------------------------------------------------------

subroutine read_2d_from_nc_file(ncid, varname, data)
 integer,          intent(in)  :: ncid
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
              'inq_varid '//trim(varname)//' '//trim(model_analysis_filename))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
              'read_from_nc_file', &
              'inquire '//trim(varname)//' '//trim(model_analysis_filename))

do i=1, numdims
   write(string1,*)'inquire length for dimension ',i
   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen, name=dimname), &
                 'read_2d_from_nc_file', &
                  trim(string1)//' '//trim(model_analysis_filename))
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
              'get_var u '//trim(model_analysis_filename))

end subroutine read_2d_from_nc_file


!------------------------------------------------------------------

subroutine put_u(ncid, u)
 integer,  intent(in) :: ncid
 real(r8), intent(in) :: u(:,:)       ! u(nVertLevels, nEdges)

! Put the newly updated 'u' field back into the netcdf file.

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount, numu
integer :: VarID, numdims, nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: ntimes, i

if ( .not. module_initialized ) call static_init_model


call nc_check(nf90_Inquire(ncid,nDimensions,nVariables,nAttributes,unlimitedDimID), &
              'put_u', 'inquire '//trim(model_analysis_filename))

call nc_check(nf90_inquire_dimension(ncid, unlimitedDimID, len=ntimes), &
              'put_u', 'inquire time dimension length '//trim(model_analysis_filename))

call nc_check(nf90_inq_varid(ncid, 'u', VarID), &
              'put_u', 'inq_varid u '//trim(model_analysis_filename))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
              'put_u', 'inquire u '//trim(model_analysis_filename))

do i=1, numdims
   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=numu(i)), &
                 'put_u', 'inquire U dimension length '//trim(model_analysis_filename))
enddo

! for all but the time dimension, read all the values.
! for time read only the last one (if more than 1 present)
mystart = 1
mystart(numdims) = ntimes
mycount = numu
mycount(numdims) = 1

call nc_check(nf90_put_var(ncid, VarID, u, start=mystart, count=mycount), &
              'put_u', 'get_var u '//trim(model_analysis_filename))


! A little sanity check

if ((debug > 9) .and. do_output()) then

   write(*,*)
   write(*,*)'u       range ',minval(u),     maxval(u)

endif

end subroutine put_u


!------------------------------------------------------------------

subroutine vector_to_1d_prog_var(x, ivar, data_1d_array)

! convert the values from a 1d array, starting at an offset,
! into a 1d array.

real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:),   intent(out) :: data_1d_array

integer :: idim1,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar(ivar)%index1

do idim1 = 1, size(data_1d_array, 1)
   data_1d_array(idim1) = x(ii)
   ii = ii + 1
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_1d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_1d_prog_var


!------------------------------------------------------------------

subroutine vector_to_2d_prog_var(x, ivar, data_2d_array)

! convert the values from a 1d array, starting at an offset,
! into a 2d array.

real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:,:), intent(out) :: data_2d_array

integer :: idim1,idim2,ii

if ( .not. module_initialized ) call static_init_model

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
   call error_handler(E_ERR,'vector_to_2d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_2d_prog_var


!------------------------------------------------------------------

subroutine vector_to_3d_prog_var(x, ivar, data_3d_array)

! convert the values from a 1d array, starting at an offset,
! into a 3d array.

real(r8), dimension(:),     intent(in)  :: x
integer,                    intent(in)  :: ivar
real(r8), dimension(:,:,:), intent(out) :: data_3d_array

integer :: idim1,idim2,idim3,ii

if ( .not. module_initialized ) call static_init_model

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
   call error_handler(E_ERR,'vector_to_3d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_3d_prog_var


!------------------------------------------------------------------

subroutine prog_var_1d_to_vector(data_1d_array, x, ivar)

! convert the values from a 1d array into a 1d array
! starting at an offset.

real(r8), dimension(:),   intent(in)    :: data_1d_array
real(r8), dimension(:),   intent(inout) :: x
integer,                  intent(in)    :: ivar

integer :: idim1,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar(ivar)%index1

do idim1 = 1, size(data_1d_array, 1)
   x(ii) = data_1d_array(idim1)
   ii = ii + 1
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' read wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'prog_var_1d_to_vector', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine prog_var_1d_to_vector


!------------------------------------------------------------------

subroutine prog_var_2d_to_vector(data_2d_array, x, ivar)

! convert the values from a 2d array into a 1d array
! starting at an offset.

real(r8), dimension(:,:), intent(in)    :: data_2d_array
real(r8), dimension(:),   intent(inout) :: x
integer,                  intent(in)    :: ivar

integer :: idim1,idim2,ii

if ( .not. module_initialized ) call static_init_model

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
   call error_handler(E_ERR,'prog_var_2d_to_vector', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine prog_var_2d_to_vector


!------------------------------------------------------------------

subroutine prog_var_3d_to_vector(data_3d_array, x, ivar)

! convert the values from a 2d array into a 1d array
! starting at an offset.

real(r8), dimension(:,:,:), intent(in)    :: data_3d_array
real(r8), dimension(:),     intent(inout) :: x
integer,                    intent(in)    :: ivar

integer :: idim1,idim2,idim3,ii

if ( .not. module_initialized ) call static_init_model

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
   call error_handler(E_ERR,'prog_var_3d_to_vector', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine prog_var_3d_to_vector


!------------------------------------------------------------------

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
logical :: failure

if ( .not. module_initialized ) call static_init_model

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
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector

   if ((debug > 0) .and. do_output()) then
      write(logfileunit,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
      write(     *     ,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
   endif

   ngood = ngood + 1
enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'verify_state_variables',string1,source,revision,revdate,text2=string2)
endif

! TJH FIXME need to add check so they cannot have both normal winds and reconstructed winds in
! DART state vector.   nsc - not sure that should be illegal.  which one is
! updated in the restart file is controlled by namelist and both could be in
! state vector for testing.

end subroutine verify_state_variables


!------------------------------------------------------------------

function is_edgedata_in_state_vector(state_variables)
  character(len=*), dimension(:),   intent(in)  :: state_variables
  logical :: is_edgedata_in_state_vector

! before we actually read and parse the input table,
! do a quick check to see if 'u' is one of the state variable names.
! add any other fields here if their data values are based on edges
! and not cell centers.

integer :: nrows, ncols, i
character(len=NF90_MAX_NAME) :: mpasname

if ( .not. module_initialized ) call static_init_model

nrows = max_state_variables
ncols = num_state_table_columns

MyLoop : do i = 1, nrows

   mpasname = trim(state_variables(ncols*i - 1))

   if (mpasname == ' ') exit MyLoop ! Found end of list.

   if (mpasname == 'u') then
      is_edgedata_in_state_vector = .true.
      return
   endif

enddo MyLoop

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
!%!    character(len=NF90_MAX_NAME) :: long_name
!%!    character(len=NF90_MAX_NAME) :: units
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
write(logfileunit,*) '  long_name   ',trim(progvar(ivar)%long_name)
write(     *     ,*) '  long_name   ',trim(progvar(ivar)%long_name)
write(logfileunit,*) '  units       ',trim(progvar(ivar)%units)
write(     *     ,*) '  units       ',trim(progvar(ivar)%units)
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
write(logfileunit,*) '  clmp range  ',progvar(ivar)%range
write(     *     ,*) '  clmp range  ',progvar(ivar)%range
write(logfileunit,*) '  clmp fail   ',progvar(ivar)%out_of_range_fail
write(     *     ,*) '  clmp fail   ',progvar(ivar)%out_of_range_fail
do i = 1,progvar(ivar)%numdims
   write(logfileunit,*) '  dimension/length/name ',i,progvar(ivar)%dimlens(i),trim(progvar(ivar)%dimname(i))
   write(     *     ,*) '  dimension/length/name ',i,progvar(ivar)%dimlens(i),trim(progvar(ivar)%dimname(i))
enddo

end subroutine dump_progvar

!------------------------------------------------------------------

subroutine print_variable_ranges(x)

! given a state vector, print out the min and max
! data values for the variables in the vector.

real(r8), intent(in) :: x(:)

integer :: ivar

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

zonal      = get_index_from_varname('uReconstructZonal')
meridional = get_index_from_varname('uReconstructMeridional')

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

subroutine define_var_dims(ncid,ivar, memberdimid, unlimiteddimid, ndims, dimids)

! set the dimids array needed to augment the natural shape of the variable
! with the two additional dimids needed by the DART diagnostic output.
integer,               intent(in)  :: ncid
integer,               intent(in)  :: ivar
integer,               intent(in)  :: memberdimid, unlimiteddimid
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: dimids

integer :: i,mydimid

ndims  = 0
dimids = 0

do i = 1,progvar(ivar)%numdims

   ! Each of these dimension names (originally from the MPAS analysis file)
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

integer,           intent(in)  :: dartkind
integer,           intent(out) :: index1
integer, optional, intent(out) :: indexN

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

function compute_pressure_at_loc(x, location)

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
real(r8) :: compute_pressure_at_loc

! convert the vertical coordinate at the given location
! to a pressure.  if the vertical is undefined
! return a fixed 1000 mb.

real(r8) :: loc_array(3), tk, values(3)
type(location_type) :: new_location
integer :: istatus, ivars(3)

string1 = 'routine not written for ocean.'
string2 = 'if needed, borrow from mpas_atm/model_mod.f90'
call error_handler(E_ERR, 'compute_pressure', string1, &
        source, revision, revdate, text2=string2)

compute_pressure_at_loc = MISSING_R8

end function compute_pressure_at_loc

!------------------------------------------------------------------

subroutine vert_convert(x, location, obs_kind, ztypeout, istatus)

! This subroutine converts a given ob/state vertical coordinate to
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
!        (3) x is the relevant DART state vector for carrying out
!            computations necessary for the vertical coordinate
!            transformations. As the vertical coordinate is only used
!            in distance computations, this is actually the "expected"
!            vertical coordinate, so that computed distance is the
!            "expected" distance. Thus, under normal circumstances,
!            x that is supplied to vert_convert should be the
!            ensemble mean. Nevertheless, the subroutine has the
!            functionality to operate on any DART state vector that
!            is supplied to it.

real(r8), dimension(:), intent(in)    :: x
type(location_type),    intent(inout) :: location
integer,                intent(in)    :: obs_kind
integer,                intent(in)    :: ztypeout
integer,                intent(out)   :: istatus

! zin and zout are the vert values coming in and going out.
! ztype{in,out} are the vert types as defined by the 3d sphere
! locations mod (location/threed_sphere/location_mod.f90)
real(r8) :: llv_loc(3)
real(r8) :: zin, zout, tk, fullp, surfp
real(r8) :: weights(3), zk_mid(3), values(3), fract(3), fdata(3)
integer  :: ztypein, i
integer  :: k_low(3), k_up(3), c(3), n
integer  :: ivars(3)
type(location_type) :: surfloc

! assume failure.
istatus = 1

! initialization
k_low = 0.0_r8
k_up = 0.0_r8
weights = 0.0_r8

! first off, check if ob is identity ob.  if so get_state_meta_data() will
! have returned location information already in the requested vertical type.
if (obs_kind < 0) then
   call get_state_meta_data(obs_kind,location)
   istatus = 0
   return
endif

! if the existing coord is already in the requested vertical units
! or if the vert is 'undef' which means no specifically defined
! vertical coordinate, return now.
ztypein  = nint(query_location(location, 'which_vert'))
if ((ztypein == ztypeout) .or. (ztypein == VERTISUNDEF)) then
   istatus = 0
   return
else
   if ((debug > 9) .and. do_output()) then
      write(string1,'(A,3X,2I3)') 'ztypein, ztypeout:',ztypein,ztypeout
      call error_handler(E_MSG, 'vert_convert',string1,source, revision, revdate)
   endif
endif

! we do need to convert the vertical.  start by
! extracting the location lon/lat/vert values.
llv_loc = get_location(location)

! the routines below will use zin as the incoming vertical value
! and zout as the new outgoing one.  start out assuming failure
! (zout = missing) and wait to be pleasantly surprised when it works.
zin     = llv_loc(3)
zout    = missing_r8

! if the vertical is missing to start with, return it the same way
! with the requested type as out.
if (zin == missing_r8) then
   location = set_location(llv_loc(1),llv_loc(2),missing_r8,ztypeout)
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

   call find_triangle_vert_indices (x, location, n, c, k_low, k_up, fract, weights, istatus)
   if(istatus /= 0) return

   zk_mid = k_low + fract
   zout = sum(weights * zk_mid)

   if ((debug > 9) .and. do_output()) then
      write(string2,'("Zk:",3F8.2," => ",F8.2)') zk_mid,zout
      call error_handler(E_MSG, 'vert_convert',string2,source, revision, revdate)
   endif

   ! ------------------------------------------------------------
   ! outgoing vertical coordinate should be 'pressure' in Pa
   ! ------------------------------------------------------------
   case (VERTISPRESSURE)

   ! Need to get base offsets for the potential temperature, density, and water
   ! vapor mixing fields in the state vector
! TJH   ivars(1) = get_progvar_index_from_kind(QTY_POTENTIAL_TEMPERATURE)
! TJH   ivars(2) = get_progvar_index_from_kind(QTY_DENSITY)
! TJH   ivars(3) = get_progvar_index_from_kind(QTY_VAPOR_MIXING_RATIO)

string1 = 'fix VERTISPRESSURE get base offsets - detritus from atmosphere'
call error_handler(E_ERR,'vert_convert',string1,source,revision,revdate)

   if (any(ivars(1:3) < 0)) then
      write(string1,*) 'Internal error, cannot find one or more of: theta, rho, qv'
      call error_handler(E_ERR, 'vert_convert',string1,source, revision, revdate)
   endif

   ! Get theta, rho, qv at the interpolated location
   call compute_scalar_with_barycentric (x, location, 3, ivars, values, istatus)
   if (istatus /= 0) return

   ! Convert theta, rho, qv into pressure
   call compute_full_pressure(values(1), values(2), values(3), zout, tk)
   if ((debug > 9) .and. do_output()) then
      write(string2,'("zout_in_pressure, theta, rho, qv:",3F10.2,F15.10)') zout, values
      call error_handler(E_MSG, 'vert_convert',string2,source, revision, revdate)
   endif

   ! ------------------------------------------------------------
   ! outgoing vertical coordinate should be 'depth' in meters
   ! ------------------------------------------------------------
   case (VERTISHEIGHT)

   call find_triangle_vert_indices (x, location, n, c, k_low, k_up, fract, weights, istatus)
   if (istatus /= 0) return

   fdata = 0.0_r8
   do i = 1, n
      fdata(i) = zGridFace(k_low(i),c(i))*(1.0_r8 - fract(i)) + &
                 zGridFace(k_up (i),c(i))*fract(i)
   enddo

   ! now have vertically interpolated values at cell centers.
   ! use horizontal weights to compute value at interp point.
   zout = sum(weights * fdata)

   ! ------------------------------------------------------------
   ! outgoing vertical coordinate should be 'scale height' (a ratio)
   ! ------------------------------------------------------------
   case (VERTISSCALEHEIGHT)

   ! Scale Height is defined here as: -log(pressure / surface_pressure)

   ! Need to get base offsets for the potential temperature, density, and water
   ! vapor mixing fields in the state vector
! TJH   ivars(1) = get_progvar_index_from_kind(QTY_POTENTIAL_TEMPERATURE)
! TJH   ivars(2) = get_progvar_index_from_kind(QTY_DENSITY)
! TJH   ivars(3) = get_progvar_index_from_kind(QTY_VAPOR_MIXING_RATIO)

string1 = 'fix vertisscaleheight get base offsets - detritus from atmosphere'
call error_handler(E_ERR,'vert_convert',string1,source,revision,revdate)

   ! Get theta, rho, qv at the interpolated location
   call compute_scalar_with_barycentric (x, location, 3, ivars, values, istatus)
   if (istatus /= 0) return

   ! Convert theta, rho, qv into pressure
   call compute_full_pressure(values(1), values(2), values(3), fullp, tk)
   if ((debug > 9) .and. do_output()) then
      write(string2,'("zout_full_pressure, theta, rho, qv:",3F10.2,F15.10)') fullp, values
      call error_handler(E_MSG, 'vert_convert',string2,source, revision, revdate)
   endif

   ! Get theta, rho, qv at the surface corresponding to the interpolated location
   surfloc = set_location(llv_loc(1), llv_loc(2), 1.0_r8, VERTISLEVEL)
   call compute_scalar_with_barycentric (x, surfloc, 3, ivars, values, istatus)
   if (istatus /= 0) return

   ! Convert surface theta, rho, qv into pressure
   call compute_full_pressure(values(1), values(2), values(3), surfp, tk)
   if ((debug > 9) .and. do_output()) then
      write(string2,'("zout_surf_pressure, theta, rho, qv:",3F10.2,F15.10)') surfp, values
      call error_handler(E_MSG, 'vert_convert',string2,source, revision, revdate)
   endif

   ! and finally, convert into scale height
   if (surfp /= 0.0_r8) then
      zout = -log(fullp / surfp)
   else
      zout = MISSING_R8
   endif

   if ((debug > 9) .and. do_output()) then
      write(string2,'("zout_in_pressure:",F10.2)') zout
      call error_handler(E_MSG, 'vert_convert',string2,source, revision, revdate)
   endif

   ! -------------------------------------------------------
   ! outgoing vertical coordinate is unrecognized
   ! -------------------------------------------------------
   case default
      write(string1,*) 'Requested vertical coordinate not recognized: ', ztypeout
      call error_handler(E_ERR,'vert_convert', string1, &
                         source, revision, revdate)

end select   ! outgoing vert type

! Returned location
location = set_location(llv_loc(1),llv_loc(2),zout,ztypeout)

! Set successful return code only if zout has good value
if(zout /= missing_r8) istatus = 0


end subroutine vert_convert

!==================================================================
! The following (private) interfaces are used for triangle interpolation
!==================================================================


!------------------------------------------------------------------

subroutine vert_interp(x, base_offset, cellid, nlevs, lower, fract, val, ier)

! Interpolates in vertical in column indexed by tri_index for a field
! with base_offset.  Vertical index is varying fastest here. Returns ier=0
! unless missing value is encounterd.

real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: base_offset
integer,  intent(in)  :: cellid
integer,  intent(in)  :: nlevs
integer,  intent(in)  :: lower
real(r8), intent(in)  :: fract
real(r8), intent(out) :: val
integer,  intent(out) :: ier

integer  :: offset
real(r8) :: lx, ux

! Default return is good
ier = 0

! Get the value at the lower and upper points
offset = base_offset + (cellid - 1) * nlevs + lower - 1
lx = x(offset)
ux = x(offset + 1)

! Check for missing value
if(lx == MISSING_R8 .or. ux == MISSING_R8) then
   ier = 2
   return
endif

! Interpolate
val = (1.0_r8 - fract)*lx + fract*ux

end subroutine vert_interp


!------------------------------------------------------------------

subroutine find_depth_bounds(depth, nbounds, bounds, lower, upper, fract, ier)

! Finds position of a given height in an array of height grid points and returns
! the index of the lower and upper bounds and the fractional offset.  ier
! returns 0 unless there is an error. Could be replaced with a more efficient
! search if there are many vertical levels.

real(r8), intent(in)  :: depth
integer,  intent(in)  :: nbounds
real(r8), intent(in)  :: bounds(nbounds)
integer,  intent(out) :: lower, upper
real(r8), intent(out) :: fract
integer,  intent(out) :: ier

! Assume that the spacing on altitudes is arbitrary and do the simple thing
! which is a linear search. Probably not worth any fancier searching unless
! models get to be huge.

integer :: i

! Initialization
ier = 0
fract = -1.0_r8
lower = -1
upper = -1

if (debug > 7) then
   print *, 'ready to check height bounds'
   print *, 'array ranges from 1 to ', nbounds
   print *, 'depth to check is ' , depth
   print *, 'array(1) = ', bounds(1)
   print *, 'array(nbounds) = ', bounds(nbounds)
endif

! Bounds check
if(depth < bounds(1)) ier = 998
if(depth > bounds(nbounds)) ier = 9998
if (ier /= 0) then
   if (debug > 7) print *, 'could not find vertical location'
   return
endif

! Search two adjacent layers that enclose the given point
do i = 2, nbounds
   if(depth <= bounds(i)) then
      lower = i - 1
      upper = i
      fract = (depth - bounds(lower)) / (bounds(upper) - bounds(lower))
      if (debug > 7) print *, 'found it.  lower, upper, fract = ', lower, upper, fract
      return
   endif
end do

! should never get here because depths above and below the grid
! are tested for at the start of this routine.  if you get here
! there is a coding error.
ier = 3
if (debug > 7) print *, 'internal code inconsistency: could not find vertical location'

end subroutine find_depth_bounds

!------------------------------------------------------------------

subroutine find_vert_level(x, loc, nc, ids, oncenters, lower, upper, fract, ier)

! given a location and 3 cell ids, return three sets of:
!  the two level numbers that enclose the given vertical
!  value plus the fraction between them for each of the 3 cell centers.

! FIXME:  this handles data at cell centers, at edges, but not
! data on faces.

real(r8),            intent(in)  :: x(:)
type(location_type), intent(in)  :: loc
integer,             intent(in)  :: nc, ids(:)
logical,             intent(in)  :: oncenters
integer,             intent(out) :: lower(:), upper(:)
real(r8),            intent(out) :: fract(:)
integer,             intent(out) :: ier

real(r8) :: lat, lon, vert, llv(3)
integer  :: verttype, i
integer  :: pt_base_offset, density_base_offset, qv_base_offset

! the plan is to take in: whether this var is on cell centers or edges,
! and the location so we can extract the vert value and which vert flag.
! compute and return the lower and upper index level numbers that
! enclose this vert, along with the fract between them.  ier is set
! in case of error (e.g. outside the grid, on dry land, etc).

! kinds we have to handle:
!vert_is_undef,    VERTISUNDEF
!vert_is_surface,  VERTISSURFACE
!vert_is_level,    VERTISLEVEL
!vert_is_pressure, VERTISPRESSURE
!vert_is_height,   VERTISHEIGHT

! unpack the location into local vars
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

! no defined vertical location (e.g. vertically integrated vals)
if (vert_is_undef(loc)) then
   ier = 12
   return
endif

! vertical is defined to be on the surface (level 1 here)
if(vert_is_surface(loc)) then
   lower(1:nc) = 1
   upper(1:nc) = 2
   fract(1:nc) = 0.0_r8
   ier = 0
   return
endif

! model level numbers (supports fractional levels)
if(vert_is_level(loc)) then
   ! FIXME: if this is W, the top is nVertLevels+1
   if (vert > nVertLevels) then
      ier = 12
      return
   endif

   ! at the top we have to make the fraction 1.0; all other
   ! integral levels can be 0.0 from the lower level.
   if (vert == nVertLevels) then
      lower(1:nc) = nint(vert) - 1   ! round down
      upper(1:nc) = nint(vert)
      fract(1:nc) = 1.0_r8
      ier = 0
      return
   endif

   lower(1:nc) = aint(vert)   ! round down
   upper(1:nc) = lower+1
   fract(1:nc) = vert - lower
   ier = 0
   return

endif

! ok, now we need to know where we are in the grid for depths or pressures
! as the vertical coordinate.


! Vertical interpolation for pressure coordinates
if(vert_is_pressure(loc) ) then
   ! Need to get base offsets for the potential temperature, density, and water
   ! vapor mixing fields in the state vector
! TJH   call get_index_range(QTY_POTENTIAL_TEMPERATURE, pt_base_offset)
! TJH   call get_index_range(QTY_DENSITY,          density_base_offset)
! TJH   call get_index_range(QTY_VAPOR_MIXING_RATIO,    qv_base_offset)

string1 = 'fix VERTISPRESSURE vertical interpolation for pressure - detritus from atmosphere'
call error_handler(E_ERR,'find_vert_level',string1,source,revision,revdate)

   do i=1, nc
      call find_pressure_bounds(x, vert, ids(i), nVertLevels, &
            pt_base_offset, density_base_offset, qv_base_offset,  &
            lower(i), upper(i), fract(i), ier)
      if ((debug > 9) .and. do_output()) &
         print '(A,5(1x,I5),1x,F10.4)', &
         '    after find_pressure_bounds: ier, i, cellid, lower, upper, fract = ', &
              ier, i, ids(i), lower(i), upper(i), fract(i)
      if ((debug > 5) .and. (ier /= 0) .and. do_output()) &
         print '(A,4(1x,I5),1x,F10.4,1x,I3,1x,F10.4)', &
        'fail in find_pressure_bounds: ier, nc, i, id, vert, lower, fract: ', &
         ier, nc, i, ids(i), vert, lower(i), fract(i)

      if (ier /= 0) return
   enddo

   return
endif

! grid is in depth, so this needs to know which cell to index into
! for the column of depths and call the bounds routine.
if(vert_is_height(loc)) then
   ! For depth, can do simple vertical search for interpolation for now
   ! Get the lower and upper bounds and fraction for each column
   do i=1, nc
      if (oncenters) then
         call find_depth_bounds(vert, nVertLevels, zGridCenter(:, ids(i)), &
                                 lower(i), upper(i), fract(i), ier)
      else

         if (.not. data_on_edges) then
            call error_handler(E_ERR, 'find_vert_level', &
                              'Internal error: oncenters false but data_on_edges false', &
                              source, revision, revdate)
         endif
         call find_depth_bounds(vert, nVertLevels, zGridEdge(:, ids(i)), &
                                 lower(i), upper(i), fract(i), ier)
      endif
      if (ier /= 0) return
   enddo

   return
endif

! If we get here, the vertical type is not understood.  Should not
! happen.
ier = 3

end subroutine find_vert_level

!------------------------------------------------------------------

subroutine find_pressure_bounds(x, p, cellid, nbounds, &
   pt_base_offset, density_base_offset, qv_base_offset, &
   lower, upper, fract, ier)

! Finds vertical interpolation indices and fraction for a quantity with
! pressure vertical coordinate. Loops through the depth levels and
! computes the corresponding pressure at the horizontal point.  nbounds is
! the number of vertical levels in the potential temperature, density,
! and water vapor grids.

real(r8),  intent(in)  :: x(:)
real(r8),  intent(in)  :: p
integer,   intent(in)  :: cellid
integer,   intent(in)  :: nbounds
integer,   intent(in)  :: pt_base_offset, density_base_offset, qv_base_offset
integer,   intent(out) :: lower, upper
real(r8),  intent(out) :: fract
integer,   intent(out) :: ier

integer  :: i, j, ier2
real(r8) :: pressure(nbounds), pr

! Initialize to bad values so unset returns will be caught.
fract = -1.0_r8
lower = -1
upper = -1

! Find the lowest pressure
call get_interp_pressure(x, pt_base_offset, density_base_offset, &
   qv_base_offset, cellid, 1, nbounds, pressure(1), ier)
!print *, 'find p bounds2, pr(1) = ', pressure(1), ier
if(ier /= 0) return

! Get the highest pressure level
call get_interp_pressure(x, pt_base_offset, density_base_offset, &
   qv_base_offset, cellid, nbounds, nbounds, pressure(nbounds), ier)
!print *, 'find p bounds2, pr(n) = ', pressure(nbounds), ier
if(ier /= 0) return

! Check for out of the column range
ier = 0
if(p > pressure(1)) ier = 88
if(p < pressure(nbounds)) ier = 888
if (ier /= 0) return

! Loop through the rest of the column from the bottom up
do i = 2, nbounds
   call get_interp_pressure(x, pt_base_offset, density_base_offset, &
      qv_base_offset, cellid, i, nbounds, pressure(i), ier)
!print *, 'find p bounds i, pr(i) = ', i, pressure(i), ier
   if (ier /= 0) return

   ! Check if pressure at lower level is higher than at upper level.
   if(pressure(i) > pressure(i-1)) then
      if ((debug > 0) .and. do_output()) then
      write(string1, *) 'lower pressure larger than upper pressure at cellid',  cellid
      write(string2, *) 'level nums, pressures: ', i-1,i,pressure(i-1),pressure(i)
      write(*,*) 'find_pressure_bounds: ', trim(string1), trim(string2)
      endif
      if ((debug > 5) .and. do_output())  then
      do j = 1, nbounds
         call get_interp_pressure(x, pt_base_offset, density_base_offset, &
            qv_base_offset, cellid, j, nbounds, pr, ier2, .true.)
      enddo
      endif

      ier = 988
      return

      !call error_handler(E_ERR, "find_pressure_bounds", string1, &
      !   source, revision, revdate, text2=string2)
   endif

   ! Is pressure between i-1 and i level?
   if(p > pressure(i)) then
      lower = i - 1
      upper = i
      if (pressure(i) == pressure(i-1)) then
         fract = 0.0_r8
      else if (log_p_vert_interp) then
         fract = exp((log(p) - log(pressure(i-1))) / &
                    (log(pressure(i)) - log(pressure(i-1))))
      else
         fract = (p - pressure(i-1)) / (pressure(i) - pressure(i-1))
      endif

      if ((debug > 9) .and. do_output()) print '(A,3F26.18,2I4,F22.18)', &
         "find_pressure_bounds: p_in, pr(i-1), pr(i), lower, upper, fract = ", &
         p, pressure(i-1), pressure(i), lower, upper, fract

      return
   endif

end do

! should never get here because pressures above and below the column
! were tested for at the start of this routine.  if you get here
! there is a coding error.
ier = 3

end subroutine find_pressure_bounds

!------------------------------------------------------------------

subroutine get_interp_pressure(x, pt_offset, density_offset, qv_offset, &
   cellid, lev, nlevs, pressure, ier, debug)

! Finds the value of pressure at a given point at model level lev

real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: pt_offset, density_offset, qv_offset
integer,  intent(in)  :: cellid
integer,  intent(in)  :: lev, nlevs
real(r8), intent(out) :: pressure
integer,  intent(out) :: ier
logical,  intent(in), optional :: debug

integer  :: offset
real(r8) :: pt, density, qv, tk


! Get the values of potential temperature, density, and vapor
offset = (cellid - 1) * nlevs + lev - 1
pt = x(pt_offset + offset)
density = x(density_offset + offset)
qv = x(qv_offset + offset)

! Error if any of the values are missing; probably will be all or nothing
if(pt == MISSING_R8 .or. density == MISSING_R8 .or. qv == MISSING_R8) then
   ier = 2
   return
endif

! Convert theta, rho, qv into pressure
call compute_full_pressure(pt, density, qv, pressure, tk)

if (present(debug)) then
   if (debug) then
      write(*,*) 'get_interp_pressure: cellid, lev', cellid, lev
      write(*,*) 'get_interp_pressure: pt,rho,qv,p,tk', pt, density, qv, pressure, tk
   endif
endif

! Default is no error
ier = 0

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

! FIXME: i want to remove this code.  does it affect the answers?
if (any(abs(weights) < roundoff)) then
   !print *, 'get_barycentric_weights due to roundoff errors: ', weights
   where (abs(weights) < roundoff) weights = 0.0_r8
   where (abs(1.0_r8 - abs(weights)) < roundoff) weights = 1.0_r8
endif
if(abs(sum(weights)-1.0_r8) > roundoff) &
   print *, 'fail in get_barycentric_weights: sum(weights) = ',sum(weights)
!end FIXME section

end subroutine get_barycentric_weights


!------------------------------------------------------------

subroutine compute_scalar_with_barycentric(x, loc, n, ival, dval, ier)
real(r8),            intent(in)  :: x(:)
type(location_type), intent(in)  :: loc
integer,             intent(in)  :: n
integer,             intent(in)  :: ival(:)
real(r8),            intent(out) :: dval(:)
integer,             intent(out) :: ier

real(r8) :: fract(3), lowval(3), uppval(3), fdata(3), weights(3)
integer  :: lower(3), upper(3), c(3), nvert, index1, k, i, nc

dval = MISSING_R8

call find_triangle_vert_indices (x, loc, nc, c, lower, upper, fract, weights, ier)
if(ier /= 0) return

dval = 0.0_r8
do k=1, n
   ! get the starting index in the state vector
   index1 = progvar(ival(k))%index1
   nvert = progvar(ival(k))%numvertical

   ! go around triangle and interpolate in the vertical
   ! t1, t2, t3 are the xyz of the cell centers
   ! c(3) are the cell ids
   do i = 1, nc
      lowval(i) = x(index1 + (c(i)-1) * nvert + lower(i)-1)
      uppval(i) = x(index1 + (c(i)-1) * nvert + upper(i)-1)
      fdata(i) = lowval(i)*(1.0_r8 - fract(i)) + uppval(i)*fract(i)
!      if((debug > 9) .and. do_output()) &
!      print '(A,I2,A,I2,5f12.5)','compute_scalar_with_barycentric: nv=',k,' ic =',i, &
!                                  lowval(i),uppval(i),fdata(i),fract(i),weights(i)

   enddo

   ! now have vertically interpolated values at cell centers.
   ! use weights to compute value at interp point.
   dval(k) = sum(weights(1:nc) * fdata(1:nc))
!print *, 'k, dval(k) = ', k, dval(k)
enddo

ier = 0

end subroutine compute_scalar_with_barycentric

!------------------------------------------------------------

subroutine find_triangle_vert_indices (x, loc, nc, c, lower, upper, fract, weights, ier)
real(r8),            intent(in)  :: x(:)
type(location_type), intent(in)  :: loc
integer,             intent(out) :: nc
integer,             intent(out) :: c(:)
integer,             intent(out) :: lower(:), upper(:)
real(r8),            intent(out) :: fract(:)
real(r8),            intent(out) :: weights(:)
integer,             intent(out) :: ier

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
      c = MISSING_R8
  lower = MISSING_R8
  upper = MISSING_R8
  fract = 0.0_r8
weights = 0.0_r8
    ier = 0
     nc = 1

! unpack the location into local vars
llv = get_location(loc)
lon  = llv(1)
lat  = llv(2)
vert = llv(3)
verttype = nint(query_location(loc))

cellid = find_closest_cell_center(lat, lon)
if ((xyzdebug > 5) .and. do_output()) &
   print *, 'closest cell center for lon/lat: ', lon, lat, cellid
if (cellid < 1) then
   if(xyzdebug > 0) print *, 'closest cell center for lon/lat: ', lon, lat, cellid
   ier = 11
   return
endif

c(1) = cellid

if (on_boundary(cellid)) then
   ier = 12
   return
endif

if (.not. inside_cell(cellid, lat, lon)) then
   ier = 13
   return
endif

! closest vertex to given point.
closest_vert = closest_vertex_ll(cellid, lat, lon)
if ((xyzdebug > 5) .and. do_output()) &
   print *, 'closest vertex for lon/lat: ', lon, lat, closest_vert

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
! FIXME: i want to remove this code.  does it affect the answers?
      if(any(weights == 1.0_r8)) then
         nc = 1
      else
!end FIXME section
         nc = 3
         c(2) = neighborcells(v)
         c(3) = neighborcells(vp1)
! FIXME: i want to remove this code.  does it affect the answers?
      endif
!end FIXME section
      foundit = .true.
      exit findtri
   endif
enddo findtri
if (.not. foundit) then
   ier = 14     ! 11
   return
endif

endif     ! horizontal index search is done now.

! need vert index for the vertical level
call find_vert_level(x, loc, nc, c, .true., lower, upper, fract, ier)
if(ier /= 0) return

if ((debug > 9) .and. do_output()) then
   write(string3,*) 'ier = ',ier, ' triangle = ',c(1:nc), ' vert_index = ',lower(1:nc)+fract(1:nc)
   call error_handler(E_MSG, 'find_triangle_vert_indices', string3, source, revision, revdate)
endif

if ((debug > 8) .and. do_output()) then
   print *, 'nc = ', nc
   print *, 'c = ', c
   print *, 'lower = ', lower(1:nc)
   print *, 'upper = ', upper(1:nc)
   print *, 'fract = ', fract(1:nc)
   print *, 'weights = ', weights
   print *, 'ier = ', ier
endif

end subroutine find_triangle_vert_indices

!------------------------------------------------------------

subroutine compute_u_with_rbf(x, loc, zonal, uval, ier)
real(r8),            intent(in)  :: x(:)
type(location_type), intent(in) :: loc
logical,             intent(in)  :: zonal
real(r8),            intent(out) :: uval
integer,             intent(out) :: ier

! max edges we currently use is ~50, but overestimate for now
integer, parameter :: listsize = 200
logical, parameter :: on_a_sphere = .true.
integer  :: nedges, edgelist(listsize), i, j, nvert
real(r8) :: xdata(listsize), ydata(listsize), zdata(listsize)
real(r8) :: edgenormals(3, listsize)
real(r8) :: veldata(listsize)
real(r8) :: xreconstruct, yreconstruct, zreconstruct
real(r8) :: ureconstructx, ureconstructy, ureconstructz
real(r8) :: ureconstructzonal, ureconstructmeridional
real(r8) :: datatangentplane(3,2)
real(r8) :: coeffs_reconstruct(3,listsize)
integer  :: index1, progindex, cellid, vertexid
real(r8) :: lat, lon, vert, llv(3), fract(listsize), lowval, uppval
integer  :: verttype, lower(listsize), upper(listsize), ncells, celllist(listsize)


! FIXME: make this cache the last value and if the location is
! the same as before and it's asking for V now instead of U,
! skip the expensive computation.

progindex = get_index_from_varname('u')
if (progindex < 0 .or. .not. data_on_edges) then
   ! cannot compute u if it isn't in the state vector, or if we
   ! haven't read in the edge data (which shouldn't happen if
   ! u is in the state vector.
   uval = MISSING_R8
   ier = 18
   return
endif
index1 = progvar(progindex)%index1
nvert = progvar(progindex)%numvertical

! unpack the location into local vars
llv = get_location(loc)
lon = llv(1)
lat = llv(2)
vert = llv(3)
verttype = nint(query_location(loc))

call find_surrounding_edges(lat, lon, nedges, edgelist, cellid, vertexid)
if (nedges <= 0) then
   ! we are on a boundary, no interpolation
   uval = MISSING_R8
   ier = 18
   return
endif

if (verttype == VERTISPRESSURE) then
   ! get all cells which share any edges with the 3 cells which share
   ! the closest vertex.
   call make_cell_list(vertexid, 3, ncells, celllist)

   call find_vert_level(x, loc, ncells, celllist, .true., &
                        lower, upper, fract, ier)
   if (ier /= 0) return

   ! now have pressure at all cell centers - need to interp to get pressure
   ! at edge centers.
   call move_pressure_to_edges(ncells, celllist, lower, upper, fract, &
                               nedges, edgelist, ier)
   if (ier /= 0) return

else
   ! need vert index for the vertical level
   call find_vert_level(x, loc, nedges, edgelist, .false., &
                        lower, upper, fract, ier)
   if (ier /= 0) return
endif

! the rbf code needs (their names == our names):
! nData == nedges
! xyz data == xyzEdge
! normalDirectionData == edgeNormalVectors
! velocitydata = U field

do i = 1, nedges
   if (edgelist(i) > size(xEdge)) then
      write(string1, *) 'edgelist has index larger than edge count', &
                          i, edgelist(i), size(xEdge)
      call error_handler(E_ERR, 'compute_u_with_rbf', 'internal error', &
                         source, revision, revdate, text2=string1)
   endif
   xdata(i) = xEdge(edgelist(i))
   ydata(i) = yEdge(edgelist(i))
   zdata(i) = zEdge(edgelist(i))

   do j=1, 3
      edgenormals(j, i) = edgeNormalVectors(j, edgelist(i))
   enddo

   lowval = x(index1 + (edgelist(i)-1) * nvert + lower(i)-1)
   uppval = x(index1 + (edgelist(i)-1) * nvert + upper(i)-1)
   veldata(i) = lowval*(1.0_r8 - fract(i)) + uppval*fract(i)
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
call get_reconstruct(nedges, lat*deg2rad, lon*deg2rad, &
              coeffs_reconstruct, on_a_sphere, veldata, &
              ureconstructx, ureconstructy, ureconstructz, &
              ureconstructzonal, ureconstructmeridional)


! FIXME: it would be nice to return both and not have to call this
! code twice.  crap.
if (zonal) then
   uval = ureconstructzonal
else
   uval = ureconstructmeridional
endif

ier = 0

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

if (.not. inside_cell(cellid, lat, lon)) then
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

! the width really isn't used anymore, but it's part of the
! interface so we have to pass some number in.
call xyz_get_close_maxdist_init(cc_gc, 1.0_r8)
call xyz_get_close_obs_init(cc_gc, nCells, cell_locs)

end subroutine init_closest_center

!------------------------------------------------------------

function find_closest_cell_center(lat, lon)

! Determine the cell index for the closest center to the given point
! 2D calculation only.

real(r8), intent(in)  :: lat, lon
integer               :: find_closest_cell_center

type(xyz_location_type) :: pointloc
integer :: closest_cell, rc
logical, save :: search_initialized = .false.

! do this exactly once.
if (.not. search_initialized) then
   call init_closest_center()
   search_initialized = .true.
endif

pointloc = xyz_set_location(lon, lat, 0.0_r8, radius)

call xyz_find_nearest(cc_gc, pointloc, cell_locs, closest_cell, rc)

! decide what to do if we don't find anything.
if (rc /= 0 .or. closest_cell < 0) then
   if ((debug > 8) .and. do_output()) &
       print *, 'cannot find nearest cell to lon, lat: ', lon, lat
   find_closest_cell_center = -1
   return
endif

! this is the cell index for the closest center
find_closest_cell_center = closest_cell

end function find_closest_cell_center

!------------------------------------------------------------

subroutine finalize_closest_center()

! get rid of storage associated with GC for cell centers.

call xyz_get_close_obs_destroy(cc_gc)

end subroutine finalize_closest_center

!------------------------------------------------------------

function find_closest_edge(cellid, lat, lon)

! given a cellid and a lat/lon, find which edge is closest
! to the location.  2D calculation only.

integer,  intent(in)  :: cellid
real(r8), intent(in)  :: lat, lon
integer               :: find_closest_edge

type(location_type) :: pointloc, edgeloc
integer :: i, closest_edge, edgeid
real(r8) :: closest_dist, dist

pointloc = set_location(lon, lat, 0.0_r8, VERTISSURFACE)

closest_edge = -1
closest_dist = 9.0e9_r8  ! something large in radians
nedges = nEdgesOnCell(cellid)
do i=1, nedges
   edgeid = edgesOnCell(i, cellid)
   edgeloc = set_location(lonEdge(edgeid), latEdge(edgeid), 0.0_r8, VERTISSURFACE)
   dist = get_dist(pointloc, edgeloc, no_vert = .true.)
   if (dist < closest_dist) then
      closest_dist = dist
      closest_edge = edgeid
   endif
enddo

! decide what to do if we don't find anything.
if (closest_edge < 0) then
   if ((debug > 8) .and. do_output()) &
      print *, 'cannot find nearest edge to lon, lat: ', lon, lat
   find_closest_edge = -1
   return
endif

! this is the edge index for the closest edge center
find_closest_edge = edgeid

end function find_closest_edge

!------------------------------------------------------------

function on_boundary(cellid)

! use the surface (level 1) to determine if any edges (or vertices?)
! are on the boundary, and return true if so.   if the global flag
! is set, skip all code and return false immediately.

integer,  intent(in)  :: cellid
logical               :: on_boundary

integer :: vertical

vertical = 1

if (global_grid) then
   on_boundary = .false.
   return
endif

write(*,*) 'boundaryCell ', cellid, boundaryCell(vertical,cellid)

on_boundary = boundaryCell(vertical,cellid) .eq. 1

end function on_boundary

!------------------------------------------------------------

function inside_cell(cellid, lat, lon)

! this function no longer really determines if we are inside
! the cell or not.  what it does do is determine if the nearest
! cell is on the grid boundary in any way and says no if it is
! a boundary.  if we have a flag saying this a global grid, we
! can avoid doing any work and immediately return true.  for a
! global atmosphere this is always so; for a regional atmosphere
! and for the ocean (which does not have cells on land) this is
! necessary test.

integer,  intent(in)  :: cellid
real(r8), intent(in)  :: lat, lon
logical               :: inside_cell

! do this completely with topology of the grid.  if any of
! the cell edges are marked as boundary edges, return no.
! otherwise return yes.

integer :: nedges, i, edgeid, vert

! if we're on a global grid, skip all this code
if (global_grid) then
   inside_cell = .true.
   return
endif

nedges = nEdgesOnCell(cellid)

! go around the edges and check the boundary array.
! if any are true, return false.  even if we are inside
! this cell, we aren't going to be able to interpolate it
! so shorten the code path.

! FIXME: at some point we can be more selective and try to
! interpolate iff the edges of the three cells which are
! going to contribute edges to the RBF exist, even if some
! of the other cell edges are on the boundary.  so this
! decision means we won't be interpolating some obs that in
! theory we have enough information to interpolate.  but it
! is conservative for now - we certainly won't try to interpolate
! outside the existing grid.

! FIXME: can this loop over edges be replaced by a single test of the 'boundaryCell' array?
do i=1, nedges
   edgeid = edgesOnCell(i, cellid)

   ! FIXME: this is an int array.  is it 0=false,1=true?
   ! BOTHER - we need the vert for this and we don't have it
   ! and in fact can't compute it if the interpolation point
   ! has pressure or depth as its vertical coordinate.
   vert = 1

   if (boundaryEdge(vert, edgeid) > 0) then
      inside_cell = .false.
      return
   endif

enddo

inside_cell = .true.

end function inside_cell

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
    (debug > 8) .and. do_output()) &
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


! FIXME: the ocean files have:
!  integer boundaryEdge(nVertLevels, nEdges)
!  integer boundaryVertex(nVertLevels, nVertices)
! as a first pass, if ANY of the edges or vertices are on
! the boundary, punt and return 0 as the edge count.  later
! once this is working, decide if a single boundary vertex or
! edge is ok if it's the exterior of the edges we are including
! and if it has good data values.

listlen = 0
do c=1, ncells
   edgecount = nEdgesOnCell(cellid_list(c))
   do e=1, edgecount
      nextedge = edgesOnCell(e, cellid_list(c))
      ! FIXME:
      ! if (boundaryEdge(nextedge, vert)) then
      !    nedges = 0
      !    edge_list(:) = -1
      !    return
      ! endif
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

! FIXME: will need a vertical level number input arg here to detect
! the boundary edges/vertices correctly.

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


! FIXME: the ocean files have:
!  integer boundaryEdge(nVertLevels, nEdges)
!  integer boundaryVertex(nVertLevels, nVertices)
! as a first pass, if ANY of the edges or vertices are on
! the boundary, punt and return 0 as the edge count.  later
! once this is working, decide if a single boundary vertex or
! edge is ok if it's the exterior of the edges we are including
! and if it has good data values.

listlen = 0
do c=1, ncells
   edgecount = nEdgesOnCell(cellids(c))
   do e=1, edgecount
      nextedge = edgesOnCell(e, cellids(c))
      ! FIXME:
      ! if (boundaryEdge(nextedge, vert)) then
      !    nedges = 0
      !    edge_list(:) = -1
      !    return
      ! endif
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


! FIXME: the ocean files have:
!  integer boundaryEdge(nVertLevels, nEdges)
!  integer boundaryVertex(nVertLevels, nVertices)
! as a first pass, if ANY of the edges or vertices are on
! the boundary, punt and return 0 as the edge count.  later
! once this is working, decide if a single boundary vertex or
! edge is ok if it's the exterior of the edges we are including
! and if it has good data values.

listlen = ncells
do c=1, ncells
   vertcount = nEdgesOnCell(cell_list(c))
   do v=1, vertcount
      nextvert = verticesOnCell(v, cell_list(c))
      ! FIXME:
      ! if (boundaryVertex(nextvert, vert)) then
      !    ncells = 0
      !    cell_list(:) = -1
      !    return
      ! endif
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
      ! FIXME:
      ! if (boundaryVertex(nextvert, vert)) then
      !    ncells = 0
      !    cell_list(:) = -1
      !    return
      ! endif
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
integer,  intent(inout) :: lower(:), upper(:)
real(r8), intent(inout) :: fract(:)
integer,  intent(in)    :: nedges, edgelist(:)
integer,  intent(out)   :: ier

integer, parameter :: listsize = 30
real(r8) :: o_lower(listsize), o_fract(listsize)
real(r8) :: x1, x2, x
integer  :: i, c1, c2, c

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
         x1 = o_lower(c) + o_fract(c)
      endif
      if (celllist(c) == c2) then
         x2 = o_lower(c) + o_fract(c)
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
   lower(i) = aint(x)
   upper(i) = lower(i) + 1
   fract(i) = x - lower(i)
enddo

ier = 0

end subroutine move_pressure_to_edges

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

subroutine latlon_to_xyz_on_plane(lat, lon, cellid, x, y, z)

! FIXME: currently unused.  unless needed, could be removed.

! Given a lat, lon in degrees, and the id of a cell in the
! MPAS grid, return the cartesian x,y,z coordinate of that
! location ON THE PLANE defined by the vertices of that cell.
! This will be different from the x,y,z of the surface of the
! sphere.  Uses the parametric form description from
! http://en.wikipedia.org/wiki/Line-plane_intersection

real(r8), intent(in)  :: lat, lon
integer,  intent(in)  :: cellid
real(r8), intent(out) :: x, y, z

integer  :: nverts, i, vertexid
real(r8) :: s(3)         ! location of point on surface
real(r8) :: p(3,3)       ! first 3 vertices of cell, xyz
real(r8) :: intp(3)      ! intersection point with plane

call latlon_to_xyz(lat, lon, s(1), s(2), s(3))

! get the first 3 vertices to define plane
! intersect with sx,sy,sz to get answer

! nedges and nverts is same
nverts = nEdgesOnCell(cellid)
if (nverts < 3) then
   print *, 'nverts is < 3', nverts
   stop
endif

! use first 3 verts to define plane
do i=1, 3
   vertexid = verticesOnCell(i, cellid)

   p(1,i) = xVertex(vertexid)
   p(2,i) = yVertex(vertexid)
   p(3,i) = zVertex(vertexid)
enddo

x = intp(1)
y = intp(2)
z = intp(3)

end subroutine latlon_to_xyz_on_plane

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

!==================================================================
! The following (private) routines were borrowed from the MPAS code
!==================================================================

!------------------------------------------------------------------

subroutine uv_cell_to_edges(zonal_wind, meridional_wind, du)

! Project u, v wind increments at cell centers onto the edges.
! FIXME:
!        we can hard-code R3 here since it comes from the (3d) x/y/z cartesian coordinate.
!        We define nEdgesOnCell in get_grid_dims, and read edgesOnCell in get_grid.
!        We read edgeNormalVectors in get_grid to use this subroutine.
!        Here "U" is the prognostic variable in MPAS, and we update it with the wind
!        increments at cell centers.

real(r8), intent(in) :: zonal_wind(:,:)             ! u wind updated from filter
real(r8), intent(in) :: meridional_wind(:,:)        ! v wind updated from filter
real(r8), intent(out):: du(:,:)                     ! normal velocity increment on the edges

! Local variables
integer, parameter :: R3 = 3
real(r8) :: east(R3,nCells), north(R3,nCells)
real(r8) :: lonCell_rad(nCells), latCell_rad(nCells)
integer  :: iCell, iEdge, jEdge, k

if ( .not. module_initialized ) call static_init_model

! Initialization
du(:,:) = 0.0_r8

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

! Project analysis increments from the cell centers to the edges
do iCell = 1, nCells
   do jEdge = 1, nEdgesOnCell(iCell)
      iEdge = edgesOnCell(jEdge, iCell)
      do k = 1, nVertLevels
         du(k,iEdge) = du(k,iEdge) + 0.5_r8 * zonal_wind(k,iCell)   &
                     * (edgeNormalVectors(1,iEdge) * east(1,iCell)  &
                     +  edgeNormalVectors(2,iEdge) * east(2,iCell)  &
                     +  edgeNormalVectors(3,iEdge) * east(3,iCell)) &
                     + 0.5_r8 * meridional_wind(k,iCell)            &
                     * (edgeNormalVectors(1,iEdge) * north(1,iCell) &
                     +  edgeNormalVectors(2,iEdge) * north(2,iCell) &
                     +  edgeNormalVectors(3,iEdge) * north(3,iCell))
      enddo
   enddo
enddo

end subroutine uv_cell_to_edges


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

function theta_to_tk (theta, rho, qv)

! Compute sensible temperature [K] from potential temperature [K].
! code matches computation done in MPAS model

real(r8), intent(in)  :: theta    ! potential temperature [K]
real(r8), intent(in)  :: rho      ! dry density
real(r8), intent(in)  :: qv       ! water vapor mixing ratio [kg/kg]
real(r8)  :: theta_to_tk          ! sensible temperature [K]

! Local variables
real(r8) :: theta_m               ! potential temperature modified by qv
real(r8) :: exner                 ! exner function
real(r8) :: qv_nonzero            ! qv >= 0

qv_nonzero = max(qv,0.0_r8)
theta_m = (1.0_r8 + 1.61_r8 * qv_nonzero)*theta

!theta_m = (1.0_r8 + 1.61_r8 * (max(qv, 0.0_r8)))*theta
exner = ( (rgas/p0) * (rho*theta_m) )**rcv

! Temperature [K]
theta_to_tk = theta * exner

end function theta_to_tk


!------------------------------------------------------------------

subroutine compute_full_pressure(theta, rho, qv, pressure, tk)

! Compute full pressure from the equation of state.
! since it has to compute sensible temp along the way,
! make temp one of the return values rather than having
! to call theta_to_tk() separately.
! code matches computation done in MPAS model

real(r8), intent(in)  :: theta    ! potential temperature [K]
real(r8), intent(in)  :: rho      ! dry density
real(r8), intent(in)  :: qv       ! water vapor mixing ratio [kg/kg]
real(r8), intent(out) :: pressure ! full pressure [Pa]
real(r8), intent(out) :: tk       ! return sensible temperature to caller

! Local variables
real(r8) :: qv_nonzero            ! qv >= 0

qv_nonzero = max(qv,0.0_r8)
tk = theta_to_tk(theta, rho, qv_nonzero)

!tk = theta_to_tk(theta, rho, max(qv,0.0_r8))
pressure = rho * rgas * tk * (1.0_r8 + 1.61_r8 * qv)
!if ((debug > 9) .and. do_output()) print *, 't,r,q,p,tk =', theta, rho, qv, pressure, tk

end subroutine compute_full_pressure


!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
