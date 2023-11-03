! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module model_mod

! This is the interface between the Community Land Model (CLM) and DART.
! The following is the hierarchy as I see it:
! top ... a gridcell has one or more more landunits
!     ... ... landunits have one or more columns
!     ... ... ... columns may have one or more layers (snow, for instance)
!     ... ... ... columns have one or more PFTs
!     ... ... ... ... some PFTs have layers (radiance bands, for instance)
!
! landunit types
! 1 soil (natural vegation/bare ground)
! 2 glacier
! 3 lake
! 4 not used (shallow lake)
! 5 wetland
! 6 urban
! 7 ice (new glacier model)
! 8 crop (if using crop model)

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, i8, MISSING_R8, MISSING_I, MISSING_R4,    &
                             obstypelength

use time_manager_mod, only : time_type, set_time, get_time, set_date, get_date,&
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+),  operator(-),          &
                             operator(>),  operator(<),  operator(/),          &
                             operator(/=), operator(<=), operator(==)

use     location_mod, only : location_type, set_location, get_location,        &
                             write_location, is_vertical, VERTISLEVEL,         &
                             VERTISHEIGHT, LocationDims, get_close_type,       &
                             convert_vertical_obs, convert_vertical_state,     &
                             loc_get_close_obs   => get_close_obs,             &
                             loc_get_close_state => get_close_state

use    utilities_mod, only : error_handler, E_ALLMSG, E_ERR, E_WARN, E_MSG,    &
                             logfileunit, get_unit, do_output, to_upper,       &
                             find_namelist_in_file, check_namelist_read,       &
                             file_exist, find_textfile_dims, file_to_text,     &
                             open_file, close_file, do_nml_file, do_nml_term,  &
                             nmlfileunit

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_check, nc_add_global_creation_time,        &
                                 nc_begin_define_mode, nc_end_define_mode,     &
                                 nc_open_file_readonly, nc_open_file_readwrite,&
                                 nc_close_file, nc_add_attribute_to_variable,  &
                                 nc_get_dimension_size, nc_get_variable_size,  &
                                 nc_get_variable, nc_put_variable, &
                                 nc_get_variable_dimension_names, &
                                 nc_define_dimension, nc_variable_exists, &
                                 nc_define_integer_variable, &
                                 nc_define_real_variable, &
                                 nc_define_double_variable

use     obs_kind_mod, only : QTY_SOIL_TEMPERATURE,       &
                             QTY_SOIL_MOISTURE,          &
                             QTY_SOIL_LIQUID_WATER,      &
                             QTY_SOIL_ICE,               &
                             QTY_SNOWCOVER_FRAC,         &
                             QTY_WATER_TABLE_DEPTH,      &
                             QTY_LEAF_CARBON,   QTY_LIVE_STEM_CARBON,   QTY_DEAD_STEM_CARBON, &
                             QTY_LEAF_NITROGEN, QTY_LIVE_STEM_NITROGEN, QTY_DEAD_STEM_NITROGEN, &
                             QTY_NET_CARBON_PRODUCTION,  &
                             QTY_LEAF_AREA_INDEX,        &
                             QTY_WATER_TABLE_DEPTH,      &
                             QTY_GEOPOTENTIAL_HEIGHT,    &
                             QTY_VEGETATION_TEMPERATURE, &
                             QTY_PAR_DIRECT,             &
                             QTY_PAR_DIFFUSE,            &
                             QTY_ABSORBED_PAR,           &
                             QTY_FRACTION_ABSORBED_PAR,  &
                             QTY_SOLAR_INDUCED_FLUORESCENCE, &
                             QTY_SOIL_CARBON,            &
                             QTY_LATENT_HEAT_FLUX,       &
                             QTY_LANDMASK,               &
                             get_index_for_quantity,     &
                             get_name_for_quantity,      &
                             get_quantity_for_type_of_obs

 use ensemble_manager_mod, only : ensemble_type, &
                                  map_pe_to_task, &
                                  get_var_owner_index, &
                                  all_copies_to_all_vars, &
                                  all_vars_to_all_copies

use distributed_state_mod, only : get_state

use   state_structure_mod, only : add_domain, state_structure_info,   &
                                  get_index_start, get_index_end,     &
                                  get_num_domains, get_num_variables, &
                                  get_num_dims, get_dim_name,         &
                                  get_dim_length, get_variable_name,  &
                                  do_io_update, get_variable_size,    &
                                  get_model_variable_indices,         &
                                  get_domain_size, get_varid_from_kind

use obs_def_utilities_mod, only : track_status

use     mpi_utilities_mod, only : my_task_id

use     default_model_mod, only : adv_1step, init_time, init_conditions, &
                                  nc_write_model_vars

use typesizes

use netcdf

implicit none
private

! These required routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
! Routines in this list have code in this module.
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          shortest_time_between_assimilations, &
          static_init_model,      &
          nc_write_model_atts,    &
          get_close_obs,          &
          get_close_state,        &
          pert_model_copies,      &
          write_model_time,       &
          read_model_time,        &
          end_model

! These required routines are in other modules.
public::  init_time,              &
          init_conditions,        &
          nc_write_model_vars,    &
          convert_vertical_obs,   &
          convert_vertical_state

! Routines for support purposes, these interfaces can be changed as appropriate.
public :: get_gridsize,                 &
          get_clm_restart_filename,     &
          compute_gridcell_value,       &
          gridcell_components

character(len=*), parameter :: source   = 'clm/model_mod.f90'
character(len=*), parameter :: revdate  = ''

character(len=512) :: string1, string2, string3

logical, save :: module_initialized = .false.

! 'Handles' for the different domains.
integer :: dom_restart = -1
integer :: dom_history = -1
integer :: dom_vector  = -1

!------------------------------------------------------------------
!
!  The DART state vector may consist of things like:
!
!  SNOWDP       aka  "snow depth"
!  frac_sno     aka  "snow cover fraction"
!  leafc        aka  "leaf carbon"
!  T_SOISNO     aka  "temperature of soil & snow"
!  H2OSOI_LIQ   aka  "liquid water in soil & snow"
!  H2OSOI_ICE   aka  "water equivalent of ice in soil & snow"
!  DZSNO        aka  "snow layer thickness"
!
!  The variables in the clm restart file that are used to create the
!  DART state vector are specified in the input.nml:model_nml namelist.
!
!------------------------------------------------------------------

! Codes for COLUMN variables
integer :: icol_vegetated_or_bare_soil
integer :: icol_crop
! character(:) :: icol_crop_noncompete    NO IDEA 
integer :: icol_landice
! character(:) :: icol_landice_multiple_elevation_classes    NO IDEA 
integer :: icol_deep_lake
integer :: icol_wetland
integer :: icol_urban_roof
integer :: icol_urban_sunwall
integer :: icol_urban_shadewall
integer :: icol_urban_impervious_road
integer :: icol_urban_pervious_road

! Codes for LANDUNIT variables
integer :: ilun_vegetated_or_bare_soil
integer :: ilun_crop
integer :: ilun_UNUSED
integer :: ilun_landice
integer :: ilun_deep_lake
integer :: ilun_wetland
integer :: ilun_urban_tbd
integer :: ilun_urban_hd
integer :: ilun_urban_md

! Codes for restricting the range of a variable
integer, parameter :: BOUNDED_NONE  = 0 ! ... unlimited range
integer, parameter :: BOUNDED_BELOW = 1 ! ... minimum, but no maximum
integer, parameter :: BOUNDED_ABOVE = 2 ! ... maximum, but no minimum
integer, parameter :: BOUNDED_BOTH  = 3 ! ... minimum and maximum

integer :: nfields
integer, parameter :: max_state_variables = 40
integer, parameter :: num_state_table_columns = 6
character(len=obstypelength) :: variable_table(max_state_variables, num_state_table_columns)

! Codes for interpreting the columns of the variable_table
integer, parameter :: VT_VARNAMEINDX  = 1 ! ... variable name
integer, parameter :: VT_KINDINDX     = 2 ! ... DART kind
integer, parameter :: VT_MINVALINDX   = 3 ! ... minimum value if any
integer, parameter :: VT_MAXVALINDX   = 4 ! ... maximum value if any
integer, parameter :: VT_ORIGININDX   = 5 ! ... file of origin
integer, parameter :: VT_STATEINDX    = 6 ! ... update (state) or not

integer :: domain_count = 0

! things which can/should be in the model_nml

integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 60
integer            :: debug = 0   ! turn up for more and more debug messages
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: clm_restart_filename = 'clm_restart.nc'
character(len=256) :: clm_history_filename = 'clm_history.nc'
character(len=256) :: clm_vector_history_filename = 'clm_vector_history.nc'

character(len=obstypelength) :: clm_variables(max_state_variables*num_state_table_columns) = ' '

namelist /model_nml/            &
   clm_restart_filename,        &
   clm_history_filename,        &
   clm_vector_history_filename, &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   calendar,                    &
   debug,                       &
   clm_variables

!----------------------------------------------------------------------
! how many and which columns are in each gridcell
!----------------------------------------------------------------------

type gridcellcolumns !  given a gridcell, which columns contribute
   private
   integer  :: ncols
   integer, pointer, dimension(:) :: columnids
   integer  :: npfts
   integer, pointer, dimension(:) :: pftids
end type gridcellcolumns
type(gridcellcolumns), allocatable, dimension(:,:), target :: gridCellInfo

!------------------------------------------------------------------------------
! Things that come from the CLM history file.
!
! The LON,LAT arrays store the longitude and latitude of grid cell centers.
! For the FV cores, there actually _are_ cells CENTERED at the poles.

integer :: nlon     = -1
integer :: nlat     = -1
integer :: nlevgrnd = -1 ! Number of 'ground' levels
integer :: nlevsoi  = -1 ! Number of 'soil' levels
integer :: nlevdcmp = -1 ! Number of 'decomposition' levels

real(r8), allocatable :: LEVGRND(:)  ! in meters
real(r8), allocatable :: LEVSOI(:)   ! in meters
real(r8), allocatable :: LEVDCMP(:)  ! in meters
real(r8), allocatable ::     LON(:)
real(r8), allocatable ::     LAT(:)
real(r8), allocatable ::  AREA1D(:),   LANDFRAC1D(:)   ! unstructured grid
real(r8), allocatable ::  AREA2D(:,:), LANDFRAC2D(:,:) ! 2D grid

logical :: unstructured = .false.

!------------------------------------------------------------------------------
! Things that come from the CLM restart file.
!
! These are the 'sparse' type arrays pertaining to the land gridcells.
!
! ZSNO  contains the height of the middle of each snow layer;
! ZISNO contains the height of the top of each snow layer;
! DZSNO tells the snow thickness of each layer.
! snow heights are stored as negative values

! Unlike the soil layers, snow layer thickness, as well as snow layer depth,
! may change as time goes on (due to snow metamorphism, overburden and the like).
! So there's no uniform levsno as levgrnd coordinate variable.

! The packing order is: Z-LON-LAT, Z moving fastest.
! all levels at a location, then
! scan along longitudes, then
! move to next latitude.

integer :: ngridcell = -1 ! Number of gridcells containing land
integer :: nlandunit = -1 ! Number of land units
integer :: ncolumn   = -1 ! Number of columns
integer :: npft      = -1 ! Number of plant functional types
integer :: nlevlak   = -1 ! Number of
integer :: nlevsno   = -1 ! Number of snow levels
integer :: nlevsno1  = -1 ! Number of snow level ... interfaces?
integer :: nlevtot   = -1 ! Number of total levels
integer :: nnumrad   = -1 ! Number of
integer :: nlevcan   = -1 ! Number of canopy layers (*XY*)

integer,  allocatable, dimension(:)  :: grid1d_ixy, grid1d_jxy ! 2D lon/lat index of corresponding gridcell
integer,  allocatable, dimension(:)  :: land1d_ixy, land1d_jxy ! 2D lon/lat index of corresponding gridcell
integer,  allocatable, dimension(:)  :: cols1d_ixy, cols1d_jxy ! 2D lon/lat index of corresponding gridcell
integer,  allocatable, dimension(:)  :: pfts1d_ixy, pfts1d_jxy ! 2D lon/lat index of corresponding gridcell
real(r8), allocatable, dimension(:)  :: land1d_wtxy    ! landunit weight relative to corresponding gridcell
real(r8), allocatable, dimension(:)  :: cols1d_wtxy    ! column   weight relative to corresponding gridcell
real(r8), allocatable, dimension(:)  :: pfts1d_wtxy    ! pft      weight relative to corresponding gridcell
integer,  allocatable, dimension(:)  :: land1d_ityplun ! landunit type
integer,  allocatable, dimension(:)  :: cols1d_ityplun ! column landunit type
integer,  allocatable, dimension(:)  :: pfts1d_ityplun ! pft landunit type
real(r8), allocatable, dimension(:)  :: cols1d_lon, cols1d_lat ! column longitude, column latitude
real(r8), allocatable, dimension(:)  :: pfts1d_lon, pfts1d_lat ! pft longitude, pft latitude
real(r8), allocatable, dimension(:)  :: levtot
real(r8), allocatable, dimension(:,:):: zsno   ! (column,levsno) ... snow layer midpoint
real(r8), allocatable, dimension(:,:):: zisno  ! (column,LEVSNO) yes, levsno ... snow layer tops

!------------------------------------------------------------------------------
! These are the metadata arrays that are the same size as the state vector.

integer,  allocatable, dimension(:) :: lonixy       ! longitude index of parent gridcell
integer,  allocatable, dimension(:) :: latjxy       ! latitude  index of parent gridcell
real(r8), allocatable, dimension(:) :: levels       ! depth
real(r8), allocatable, dimension(:) :: landarea     ! land area ... 'support' ... 'weight'

!------------------------------------------------------------------------------
! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.

logical :: print_timestamps = .false.
integer :: print_every_Nth  = 10000

!------------------------------------------------------------------
! module storage
!------------------------------------------------------------------

integer(i8)     :: model_size      ! the state vector length
type(time_type) :: model_time      ! valid time of the model state
type(time_type) :: model_timestep  ! smallest time to adv model


INTERFACE get_state_time
      MODULE PROCEDURE get_state_time_ncid
      MODULE PROCEDURE read_model_time
END INTERFACE


contains

!==================================================================
! All the REQUIRED interfaces come first - just by convention.
!==================================================================


!------------------------------------------------------------------
!> Returns the size of the model as an integer.
!> Required for all applications.

function get_model_size()

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size


!------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated array indices for lat, lon, and height, as well as the type.

subroutine get_state_meta_data(indx, location, var_type)

integer(i8),         intent(in)  :: indx
type(location_type), intent(out) :: location
integer, OPTIONAL,   intent(out) :: var_type

! Local variables

integer :: var_id, dom_id
integer :: ip, jp, kp   ! required, but unused in this context

! Module variables

! LON
! LAT
! lonixy
! latjxy
! levels

if ( .not. module_initialized ) call static_init_model

location = set_location( LON(lonixy(indx)), LAT(latjxy(indx)), levels(indx), VERTISHEIGHT)

if (present(var_type)) then

   var_type = MISSING_I

   ! from the dart index get the local variables indices
   call get_model_variable_indices(indx, ip, jp, kp, &
            var_id=var_id, dom_id=dom_id, kind_index=var_type)

   if( var_type == MISSING_I ) then
      write(string1,*) 'Cannot find DART QTY  for indx ', indx
      write(string2,*) 'variable "'//trim(get_variable_name(dom_id, var_id))//'"'
      call error_handler(E_ERR,'get_state_meta_data',string1,source,text2=string2)
   endif

endif

end subroutine get_state_meta_data


!------------------------------------------------------------------
!> Returns the the time step of the model; the smallest increment
!> in time that the model is capable of advancing the state in a given
!> implementation. This interface is required for all applications.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations= model_timestep

end function shortest_time_between_assimilations


!------------------------------------------------------------------
!> Called to do one time initialization of the model.

subroutine static_init_model()

! Local variables - all the important ones have module scope

character(len=*), parameter :: routine = 'static_init_model'

integer                      :: dimlens(NF90_MAX_VAR_DIMS)
character(len=obstypelength) :: dimnames(NF90_MAX_VAR_DIMS)

integer :: ncid
integer :: iunit, io, idom, ivar
integer :: j, rank
integer(i8) :: indx
integer :: ss, dd

integer :: nvars
character(len=obstypelength) :: var_names(max_state_variables)
real(r8) :: var_ranges(max_state_variables,2)
logical  :: var_update(max_state_variables)
integer  :: var_qtys(  max_state_variables)

character(len=obstypelength) :: varname

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

!---------------------------------------------------------------
! Set the time step ... causes clm namelists to be read.
! Ensures model_timestep is multiple of 'dynamics_timestep'

call set_calendar_type( calendar )   ! comes from model_mod_nml

model_time     = get_state_time(clm_restart_filename)
model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,routine,string1,source)

!---------------------------------------------------------------
! The CLM history file (h0?) has the 'superset' of information.
! The CLM restart files are intentionally lean and, in so doing,
! do not have the full lat/lon arrays nor any depth information.
!
ncid = 0; ! signal that netcdf file is closed
call get_history_dims(ncid, clm_history_filename, 'open', nlon, nlat, &
                      nlevgrnd, nlevsoi, nlevdcmp)

if (unstructured) then
   allocate( AREA1D(nlon),      LANDFRAC1D(nlon) )
else
   allocate( AREA2D(nlon,nlat), LANDFRAC2D(nlon,nlat) )
endif

allocate(LON(nlon), LAT(nlat),  LEVGRND(nlevgrnd), LEVSOI(nlevsoi), LEVDCMP(nlevdcmp))

call get_clm_codes(ncid, clm_history_filename, 'open')
call get_full_grid(ncid, clm_history_filename, 'close')

ncid = 0; ! signal that netcdf file is closed

!---------------------------------------------------------------
! The CLM grid in a restart file is fundamentally a sparse matrix
! representation that lacks the native grid dimensions.
! The full lat/lon arrays do not exist in the restart files.
! only the grid cells that contain land are preserved.

call get_sparse_dims(ncid, clm_restart_filename, 'open')

allocate(grid1d_ixy(ngridcell), grid1d_jxy(ngridcell))
allocate(land1d_ixy(nlandunit),     land1d_jxy(nlandunit),   land1d_wtxy(nlandunit))
allocate(cols1d_ixy(ncolumn),       cols1d_jxy(ncolumn),     cols1d_wtxy(ncolumn))
allocate(pfts1d_ixy(npft),          pfts1d_jxy(npft),        pfts1d_wtxy(npft))
allocate(land1d_ityplun(nlandunit), cols1d_ityplun(ncolumn), pfts1d_ityplun(npft))
allocate(cols1d_lon(ncolumn),       cols1d_lat(ncolumn))
allocate(pfts1d_lon(npft),          pfts1d_lat(npft))
allocate(levtot(nlevtot))

if (nlevsno > 0) then
   allocate( zsno(nlevsno,ncolumn))
   allocate(zisno(nlevsno1,ncolumn))
endif

call get_sparse_geog(ncid, clm_restart_filename, 'close')

!---------------------------------------------------------------
! Generate list of columns in each gridcell

allocate(gridCellInfo(nlon,nlat))
call SetLocatorArrays()

!---------------------------------------------------------------
! Compile the list of clm variables to use in the creation
! of the DART state vector.

call parse_variable_table( clm_variables, nfields, variable_table )

! Must group the variables according to the file they come from.
!>@todo FIXME ... io_filenames_nml:rpointer_file order must somehow match 
!> - or be insensitive to - the add_domain() calls below (if nvars == 0) ...

model_size = 0_i8

call cluster_variables(variable_table, 'RESTART', nvars, var_names, &
                       var_qtys, var_ranges, var_update)
if (nvars > 0) then
   dom_restart = add_domain(clm_restart_filename, nvars, var_names(1:nvars),  &
                  kind_list   = var_qtys(  1:nvars),   &
                  clamp_vals  = var_ranges(1:nvars,:), &
                  update_list = var_update(1:nvars) )
   model_size = model_size + get_domain_size(dom_restart)
   if (debug > 1)  call state_structure_info(dom_restart)
endif

call cluster_variables(variable_table, 'HISTORY', nvars, var_names, &
                       var_qtys, var_ranges, var_update)
if (nvars > 0) then
   dom_history = add_domain(clm_history_filename, nvars, var_names(1:nvars), &
                  kind_list   = var_qtys(  1:nvars),   &
                  clamp_vals  = var_ranges(1:nvars,:), &
                  update_list = var_update(1:nvars) )
   model_size = model_size + get_domain_size(dom_history)
   if (debug > 1)  call state_structure_info(dom_history)
endif

call cluster_variables(variable_table, 'VECTOR', nvars, var_names, &
                       var_qtys, var_ranges, var_update)
if (nvars > 0) then
   dom_vector = add_domain(clm_vector_history_filename, nvars, var_names(1:nvars), &
                  kind_list   = var_qtys(  1:nvars),   &
                  clamp_vals  = var_ranges(1:nvars,:), &
                  update_list = var_update(1:nvars) )
   model_size = model_size + get_domain_size(dom_vector)
   if (debug > 1)  call state_structure_info(dom_vector)
endif

if ((debug > 0) .and. do_output()) then
  write(string1,'(" grid: nlon, nlat =",2(1x,i6))') nlon, nlat
  write(string2, *)'model_size = ', model_size
  call say('')
  call say(string1)
  call say(string2)
endif

! Create the metadata arrays that are the same shape as the state vector.
! The metadata arrays will provide the ability to determine what grid cell is the parent
! of the state vector index in question ... as well as the actual surface area.
! By using the state structure, we are guaranteed to stride through the state vector 
! the same way the state vector is filled. This makes get_state_meta_data() quite simple.
!>@todo remove these huge allocations

allocate(lonixy(model_size), latjxy(model_size), levels(model_size), landarea(model_size))

! Initialize all levels to surface. If there is a level, we will explicitly specify it.
levels(:) = 0.0_r8

DOMAINS : do idom = 1,get_num_domains()
VARIABLES : do ivar=1,get_num_variables(idom)

   varname = get_variable_name(idom,ivar)
   indx    = get_index_start(  idom,ivar)
   rank    = get_num_dims(     idom,ivar)
   do j = 1,rank
      dimnames(j) = get_dim_name(  idom,ivar,j)
      dimlens( j) = get_dim_length(idom,ivar,j)
   enddo

   ! 'time' dimensions do not count toward 'rank'

   if (rank == 1) then
      call fill_rank1_metadata(varname, dimnames(1), dimlens(1), indx)
   elseif (rank == 2) then
      call fill_rank2_metadata(varname, dimnames(1:2), dimlens(1:2), indx)
   elseif (rank == 3) then
      call fill_rank3_metadata(varname, dimnames(1:3), dimlens(1:3), indx)
   else
      write(string1,*)'variables of rank ',rank,' are unsupported.'
      write(string2,*)trim(varname),' is dimensioned ', dimlens(1:rank)
      call error_handler(E_ERR,routine,string1,source,text2=string2)
   endif

enddo VARIABLES
enddo DOMAINS

end subroutine static_init_model


!------------------------------------------------------------------
!> Does any shutdown and clean-up needed for model. Can be a NULL
!> INTERFACE if the model has no need to clean up storage, etc.

subroutine end_model()

if (unstructured) then
   deallocate(AREA1D, LANDFRAC1D)
else
   deallocate(AREA2D, LANDFRAC2D)
endif

deallocate(LAT, LON, LEVGRND, LEVSOI, LEVDCMP)
deallocate(grid1d_ixy, grid1d_jxy)
deallocate(land1d_ixy, land1d_jxy, land1d_wtxy, land1d_ityplun)
deallocate(cols1d_ixy, cols1d_jxy, cols1d_wtxy, cols1d_ityplun)
deallocate(pfts1d_ixy, pfts1d_jxy, pfts1d_wtxy, pfts1d_ityplun)
deallocate(cols1d_lon, cols1d_lat, pfts1d_lon, pfts1d_lat)
deallocate(lonixy, latjxy, landarea)

end subroutine end_model


!------------------------------------------------------------------
!> Writes the model-specific attributes to a netCDF file.
!> This includes coordinate variables and some metadata, but NOT
!> the model state vector.

subroutine nc_write_model_atts( ncid, domain_id ) 

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

character(len=*), parameter :: routine = 'nc_write_model_atts'

if ( .not. module_initialized ) call static_init_model

! Put file into define mode.
! Write Global Attributes

call nc_begin_define_mode(ncid)
call nc_add_global_creation_time(ncid)
call nc_add_global_attribute(ncid, "model_source", source)
call nc_add_global_attribute(ncid, "model", "CLM")

!----------------------------------------------------------------------------
! Define the new dimensions IDs
!----------------------------------------------------------------------------

call nc_define_dimension(ncid,'lon', nlon, routine)
call nc_define_dimension(ncid,'lat', nlat, routine)

call nc_define_dimension(ncid,'gridcell', ngridcell, routine)
call nc_define_dimension(ncid,'landunit', nlandunit, routine)
call nc_define_dimension(ncid,'column',   ncolumn,   routine)
call nc_define_dimension(ncid,'pft',      npft,      routine)
call nc_define_dimension(ncid,'levgrnd',  nlevgrnd,  routine)
call nc_define_dimension(ncid,'levsoi',   nlevsoi,   routine)
call nc_define_dimension(ncid,'levdcmp',  nlevdcmp,  routine)
call nc_define_dimension(ncid,'levlak',   nlevlak,   routine)
call nc_define_dimension(ncid,'levsno',   nlevsno,   routine)
call nc_define_dimension(ncid,'levsno1',  nlevsno1,  routine)
call nc_define_dimension(ncid,'levtot',   nlevtot,   routine)
call nc_define_dimension(ncid,'numrad',   nnumrad,   routine)

if (unstructured) call nc_define_dimension(ncid,'lndgrid', ngridcell, routine)

if (nlevcan > 0) call nc_define_dimension(ncid,'levcan',nlevcan,routine)

!----------------------------------------------------------------------------
! Create the (empty) Coordinate Variables and the Attributes
!----------------------------------------------------------------------------

! Grid Longitudes

call nc_define_real_variable(     ncid,'lon','lon',routine)
call nc_add_attribute_to_variable(ncid,'lon','long_name','coordinate longitude',routine)
call nc_add_attribute_to_variable(ncid,'lon','cartesian_axis', 'X',routine)
call nc_add_attribute_to_variable(ncid,'lon','units', 'degrees_east',routine)
call nc_add_attribute_to_variable(ncid,'lon','valid_range', (/ 0.0_r8, 360.0_r8 /),routine)

! Grid Latitudes

call nc_define_real_variable(     ncid,'lat','lat',routine)
call nc_add_attribute_to_variable(ncid,'lat','long_name','coordinate latitude',routine)
call nc_add_attribute_to_variable(ncid,'lat','cartesian_axis', 'Y',routine)
call nc_add_attribute_to_variable(ncid,'lat','units', 'degrees_north',routine)
call nc_add_attribute_to_variable(ncid,'lat','valid_range', (/ -90.0_r8, 90.0_r8 /),routine)

! subsurface levels

call nc_define_real_variable(     ncid,'levgrnd','levgrnd',routine)
call nc_add_attribute_to_variable(ncid,'levgrnd','long_name','coordinate ground levels',routine)
call nc_add_attribute_to_variable(ncid,'levgrnd','cartesian_axis', 'Z',routine)
call nc_add_attribute_to_variable(ncid,'levgrnd','units', 'm',routine)

call nc_define_real_variable(     ncid,'levsoi','levsoi',routine) 
call nc_add_attribute_to_variable(ncid,'levsoi','long_name', &
          'coordinate soil levels (equivalent to top nlevsoi levels of levgrnd)',routine)
call nc_add_attribute_to_variable(ncid,'levsoi','cartesian_axis', 'Z',routine)
call nc_add_attribute_to_variable(ncid,'levsoi','units', 'm',routine)

call nc_define_real_variable(     ncid,'levdcmp','levdcmp',routine) 
call nc_add_attribute_to_variable(ncid,'levdcmp','long_name', &
                  'coordinate levels for soil decomposition variables',routine)
call nc_add_attribute_to_variable(ncid,'levdcmp','cartesian_axis', 'Z',routine)
call nc_add_attribute_to_variable(ncid,'levdcmp','units', 'm',routine)

call nc_define_real_variable(     ncid,'levlak','levlak',routine) 
call nc_add_attribute_to_variable(ncid,'levlak','long_name','coordinate lake levels',routine)
call nc_add_attribute_to_variable(ncid,'levlak','cartesian_axis', 'Z',routine)
call nc_add_attribute_to_variable(ncid,'levlak','units', 'm',routine)

! grid cell areas

if (unstructured) then
   call nc_define_real_variable(ncid,'area','lon',routine)
else
   call nc_define_real_variable(ncid,'area',(/'lon', 'lat'/),routine)
endif
call nc_add_attribute_to_variable(ncid,'area','long_name','grid cell areas',routine)
call nc_add_attribute_to_variable(ncid,'area','units','km^2',routine)

! grid cell land fractions

if (unstructured) then
   call nc_define_real_variable(ncid,'landfrac','lon',routine) 
else
   call nc_define_real_variable(ncid,'landfrac',(/ 'lon', 'lat' /),routine) 
endif
call nc_add_attribute_to_variable(ncid,'landfrac','long_name','land fraction',routine)
call nc_add_attribute_to_variable(ncid,'landfrac','units','km^2',routine)

! ---------
! landunits

! longitude grid index for each landunit
call nc_define_integer_variable(  ncid,'land1d_ixy','landunit',routine) 
call nc_add_attribute_to_variable(ncid,'land1d_ixy','long_name', &
                '2d longitude index of corresponding landunit',routine)

! latitude grid index for each land
call nc_define_integer_variable(  ncid,'land1d_jxy','landunit',routine) 
call nc_add_attribute_to_variable(ncid,'land1d_jxy','long_name', &
                '2d latitude index of corresponding landunit',routine)

! land weight relative to corresponding gridcell
call nc_define_double_variable(   ncid,'land1d_wtxy','landunit',routine) 
call nc_add_attribute_to_variable(ncid,'land1d_wtxy','long_name', &
           'landunit weight relative to corresponding gridcell',routine)

! land type of each land
call nc_define_integer_variable(  ncid,'land1d_ityplun','landunit',routine) 
call nc_add_attribute_to_variable(ncid,'land1d_ityplun','long_name', 'landunit type',routine)

! ---------
! columns

! longitude grid index for each column
call nc_define_integer_variable(  ncid,'cols1d_ixy','column',routine) 
call nc_add_attribute_to_variable(ncid,'cols1d_ixy','long_name', &
                '2d longitude index of corresponding column',routine)

! latitude grid index for each column
call nc_define_integer_variable(  ncid,'cols1d_jxy','column',routine) 
call nc_add_attribute_to_variable(ncid,'cols1d_jxy','long_name', &
                '2d latitude index of corresponding column',routine)

! column weight relative to corresponding gridcell
call nc_define_double_variable(   ncid,'cols1d_wtxy','column',routine) 
call nc_add_attribute_to_variable(ncid,'cols1d_wtxy','long_name', &
           'column weight relative to corresponding gridcell',routine)

! column longitude
call nc_define_double_variable(   ncid,'cols1d_lon','column',routine) 
call nc_add_attribute_to_variable(ncid,'cols1d_lon','long_name', 'column longitude',routine)

! column latitude
call nc_define_double_variable(   ncid,'cols1d_lat','column',routine) 
call nc_add_attribute_to_variable(ncid,'cols1d_lat','long_name', 'column latitude',routine)

! column land type
call nc_define_integer_variable(  ncid,'cols1d_ityplun','column',routine) 
call nc_add_attribute_to_variable(ncid,'cols1d_ityplun','long_name', &
              'column landunit type (vegetated,urban,lake,wetland or glacier)',routine)

! ---------
! patches

! longitude grid index for each pft
call nc_define_integer_variable(  ncid,'pfts1d_ixy','pft',routine) 
call nc_add_attribute_to_variable(ncid,'pfts1d_ixy','long_name', &
                '2d longitude index of corresponding pft',routine)

! latitude grid index for each pft
call nc_define_integer_variable(  ncid,'pfts1d_jxy','pft',routine) 
call nc_add_attribute_to_variable(ncid,'pfts1d_jxy','long_name', &
                '2d latitude index of corresponding pft',routine)

! pft weight relative to corresponding gridcell
call nc_define_double_variable(   ncid,'pfts1d_wtxy','pft',routine) 
call nc_add_attribute_to_variable(ncid,'pfts1d_wtxy','long_name', &
           'pft weight relative to corresponding gridcell',routine)

! pft longitude
call nc_define_double_variable(   ncid,'pfts1d_lon','pft',routine) 
call nc_add_attribute_to_variable(ncid,'pfts1d_lon','long_name', 'pft longitude',routine)

! pft latitude
call nc_define_double_variable(   ncid,'pfts1d_lat','pft',routine) 
call nc_add_attribute_to_variable(ncid,'pfts1d_lat','long_name', 'pft latitude',routine)

! pft land type
call nc_define_integer_variable(  ncid,'pfts1d_ityplun','pft',routine) 
call nc_add_attribute_to_variable(ncid,'pfts1d_ityplun','long_name', &
                                          'pft landunit type',routine)

!----------------------------------------------------------------------------
! Finished with dimension/variable definitions, must end 'define' mode to fill.
!----------------------------------------------------------------------------

call nc_end_define_mode(ncid)

!----------------------------------------------------------------------------
! Fill the coordinate variables
!----------------------------------------------------------------------------

call nc_put_variable(ncid,'lon',LON,routine)
call nc_put_variable(ncid,'lat',LAT,routine)
call nc_put_variable(ncid,'levgrnd',LEVGRND,routine)
call nc_put_variable(ncid,'levsoi',LEVSOI,routine)
call nc_put_variable(ncid,'levdcmp',LEVDCMP,routine)

if (unstructured) then
   call nc_put_variable(ncid,'area',    AREA1D,    routine)
   call nc_put_variable(ncid,'landfrac',LANDFRAC1D,routine)
else
   call nc_put_variable(ncid,'area',    AREA2D,    routine)
   call nc_put_variable(ncid,'landfrac',LANDFRAC2D,routine)
endif

call nc_put_variable(ncid,'cols1d_ixy',     cols1d_ixy,     routine)
call nc_put_variable(ncid,'cols1d_jxy',     cols1d_jxy,     routine)
call nc_put_variable(ncid,'cols1d_wtxy',    cols1d_wtxy,    routine)
call nc_put_variable(ncid,'cols1d_lon',     cols1d_lon,     routine)
call nc_put_variable(ncid,'cols1d_lat',     cols1d_lat,     routine)
call nc_put_variable(ncid,'cols1d_ityplun', cols1d_ityplun, routine)

call nc_put_variable(ncid,'pfts1d_ixy',     pfts1d_ixy,     routine)
call nc_put_variable(ncid,'pfts1d_jxy',     pfts1d_jxy,     routine)
call nc_put_variable(ncid,'pfts1d_wtxy',    pfts1d_wtxy,    routine)
call nc_put_variable(ncid,'pfts1d_lon',     pfts1d_lon,     routine)
call nc_put_variable(ncid,'pfts1d_lat',     pfts1d_lat,     routine)
call nc_put_variable(ncid,'pfts1d_ityplun', pfts1d_ityplun, routine)

! Flush the buffer and leave netCDF file open

call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts


!-------------------------------------------------------------------------------
!> Actually, nothing special about get_close_obs for CLM. Including it here
!> because we may need a unique one in the future. Right now it is simply
!> a pass-through to the default routine.

subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc
type(location_type),           intent(inout) :: base_loc
type(location_type),           intent(inout) :: locs(:)
integer,                       intent(in)    :: base_type
integer,                       intent(in)    :: loc_qtys(:)
integer,                       intent(in)    :: loc_types(:)
integer,                       intent(out)   :: num_close
integer,                       intent(out)   :: close_ind(:)
real(r8),            optional, intent(out)   :: dist(:)
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_obs'

! default version 

call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                          num_close, close_ind, dist, ens_handle)

end subroutine get_close_obs


!----------------------------------------------------------------------------
!> Routine to enforce a-priori assumptions about correlations. Manipulate the
!> distance between two objects to put them outside the localization radius if
!> we expect or know them to be uncorrelated.

subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc           !< handle to a get_close structure
type(location_type),           intent(inout) :: base_loc     !< location of interest
integer,                       intent(in)    :: base_type    !< observation TYPE
type(location_type),           intent(inout) :: locs(:)      !< locations on my task
integer,                       intent(in)    :: loc_qtys(:)  !< QTYs for locations on my task
integer(i8),                   intent(in)    :: loc_indx(:)  !< indices into DART state on my task
integer,                       intent(out)   :: num_close    !< how many are close
integer,                       intent(out)   :: close_ind(:) !< indices into the locs array
real(r8),            optional, intent(out)   :: dist(:)      !< distances (in radians)
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_state'

integer     :: iloc, obsqty
integer(i8) :: full_index

integer :: iunit

! Despite the fact dist is an optional argument, the only routine to
! call get_close_state() is in assim_tools_mod.f90 and in both places
! it has the dist argument. Consequently, we do not have to modify
! the num_close, close_ind, or dist variables - we can simply set the
! distance to be really far away.
! (although it feels better to modify num_close, etc.)

call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                            num_close, close_ind, dist, ens_handle)

! Convert the observation type to a DART Quantity.
obsqty = get_quantity_for_type_of_obs(base_type)

! Check to see if the Quantity is supposed to modify whatever
! is pointed to by the current state index.
do iloc=1, num_close

   full_index = loc_indx(close_ind(iloc))

   if (.not. related(obsqty, full_index) ) then
      if (present(dist)) dist(iloc) = huge(1.0_r8) ! really far away
   endif

enddo

! The close_ind indices are only the subset of indices that are on this particular
! task. Use these indices to reference the absolute index into the state vector to
! call get_state_meta_data() or access one of the arrays that indicates
! if this is a column or pft or ... and then we can determine if it is a lake
! or an urban area or ...

if (debug > 99) then
   iunit = my_task_id() + 100
   call write_location(iunit,base_loc)
   write(string1,'("PE ",I4)') my_task_id()
   do iloc = 1,num_close
      write(iunit,*)trim(string1),' loc_get_close_state: ',iloc, &
                dist(iloc), &
                close_ind(iloc), &
                loc_indx(close_ind(iloc))
   enddo
   write(iunit,*)''
endif

end subroutine get_close_state


!------------------------------------------------------------------
!> Perturbs a single model state for generating initial ensembles.
!> This (required interface) is unsupported in CLM and any attempt
!> to use it will cause DART to terminate.
!>
!> So far, we have generated intial ensembles by taking a single
!> instance and replicating it N times - and pairing each of the
!> instances with a unique atmospheric forcing file and integrating
!> for some period of time till the ensemble demonstrates sufficient
!> spread. This is an area of active research.
!>
!> The naive approach does not work -- it generates negative
!> snow cover fractions, for example.  Must check for out-of-range
!> values specific to each type.
!> The WRF model mod has something that might be useful.

subroutine pert_model_copies(ens_handle, ens_size, pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: ens_size
real(r8),            intent(in)    :: pert_amp
logical,             intent(out)   :: interf_provided

if ( .not. module_initialized ) call static_init_model

call error_handler(E_WARN,'pert_model_copies', &
                  'CLM cannot be started from a single vector', &
                  source, &
                  text2='see comments in clm/model_mod.f90::pert_model_copies()')

! Should provide a minimal pert routine that only perturbs some variables on columns
! or pfts that we like ... 'vegetated or bare ground, crop' ... but not lake, glacier, etc.

interf_provided = .false.

end subroutine pert_model_copies


!------------------------------------------------------------------
!> reads the valid time of the CLM model state from the restart file
!> variables timemgr_rst_curr_ymd, timemgr_rst_curr_tod

function read_model_time(filename)
type(time_type) :: read_model_time
character(len=*), intent(in) :: filename

integer :: ncid

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'read_model_time',string1,source)
endif

ncid = nc_open_file_readonly(filename,'read_model_time')

read_model_time = get_state_time_ncid(ncid)

call nc_close_file(ncid, 'read_model_time', filename)

end function read_model_time


!-----------------------------------------------------------------------
!> this routine writes the model time when creating files from scratch
!>
!>	int timemgr_rst_start_ymd ;
!>		timemgr_rst_start_ymd:long_name = "start date" ;
!>		timemgr_rst_start_ymd:units = "YYYYMMDD" ;
!>		timemgr_rst_start_ymd:interpinic_flag = 3 ;
!>		timemgr_rst_start_ymd:interpinic_flag_meanings = "1=nearest neighbor, 2=copy directly, 3=skip" ;
!>		timemgr_rst_start_ymd:varnames_on_old_files = "timemgr_rst_start_ymd" ;
!>		timemgr_rst_start_ymd:_FillValue = -999999999 ;
!>		timemgr_rst_start_ymd:missing_value = -9999 ;
!>	int timemgr_rst_start_tod ;
!>		timemgr_rst_start_tod:long_name = "start time of day" ;
!>		timemgr_rst_start_tod:units = "sec" ;
!>		timemgr_rst_start_tod:interpinic_flag = 3 ;
!>		timemgr_rst_start_tod:interpinic_flag_meanings = "1=nearest neighbor, 2=copy directly, 3=skip" ;
!>		timemgr_rst_start_tod:varnames_on_old_files = "timemgr_rst_start_tod" ;
!>		timemgr_rst_start_tod:_FillValue = -999999999 ;
!>		timemgr_rst_start_tod:missing_value = -9999 ;

subroutine write_model_time(ncid, dart_time)

integer,             intent(in) :: ncid
type(time_type),     intent(in) :: dart_time

character(len=*), parameter :: routine = 'write_model_time'
integer :: iyear, imonth, iday, ihour, imin, isec
integer :: rst_curr_ymd, rst_curr_tod
integer :: ymdVarID, todVarID
integer :: io, io1, io2
logical :: defining

call get_date(dart_time, iyear, imonth, iday, ihour, imin, isec)

rst_curr_ymd = iyear*10000 + imonth*100 + iday
rst_curr_tod = ihour*3600  + imin*60    + isec

io1 = nf90_inq_varid(ncid, 'timemgr_rst_curr_ymd', ymdVarID)
io2 = nf90_inq_varid(ncid, 'timemgr_rst_curr_tod', todVarID)

defining = .false.

if (io1 /= NF90_NOERR .or. io2 /= NF90_NOERR) defining = .true.

if ( defining ) call nc_check(nf90_redef(ncid), routine, 'redef')

if ( io1 /= NF90_NOERR ) then
   io = nf90_def_var(ncid,'timemgr_rst_curr_ymd',NF90_INT,ymdVarID)
   call nc_check(io, routine, 'defining timemgr_rst_curr_ymd')
   call nc_add_attribute_to_variable(ncid,'timemgr_rst_curr_ymd','long_name','start date',routine)
   call nc_add_attribute_to_variable(ncid,'timemgr_rst_curr_ymd','units','YYYYMMDD',routine)
endif

if ( io2 /= NF90_NOERR ) then
   io = nf90_def_var(ncid,'timemgr_rst_curr_tod',NF90_INT,todVarID)
   call nc_check(io, routine, 'defining timemgr_rst_curr_tod')
   call nc_add_attribute_to_variable(ncid,'timemgr_rst_curr_tod','long_name','start time of day',routine)
   call nc_add_attribute_to_variable(ncid,'timemgr_rst_curr_tod','units','sec',routine)
endif

if ( defining ) call nc_check(nf90_enddef(ncid), routine, 'enddef')

io = nf90_put_var(ncid, ymdVarID, rst_curr_ymd)
call nc_check(io, routine, 'put_var timemgr_rst_curr_ymd')

io = nf90_put_var(ncid, todVarID, rst_curr_tod)
call nc_check(io, routine, 'put_var timemgr_rst_curr_tod')

!>@todo maybe we also want to write the time in a DART-like way.

end subroutine write_model_time

!==================================================================
! The remaining PUBLIC interfaces come next
!==================================================================


!------------------------------------------------------------------
!>

subroutine get_gridsize(num_lon, num_lat, num_lev)

integer, intent(out) :: num_lon, num_lat, num_lev

if ( .not. module_initialized ) call static_init_model

 num_lon = nlon
 num_lat = nlat
 num_lev = nlevtot

end subroutine get_gridsize


!==================================================================
! The remaining required interfaces
!==================================================================


!-----------------------------------------------------------------------
!> Interpolate any QUANTITY of the model state to any arbitrary location.
!> A status of 0 indicates a successful interpolation.

subroutine model_interpolate(state_handle, ens_size, location, obs_kind, &
                             expected_obs, istatus)

! Reconstructing the vertical profile of the gridcell is complicated.
! Each land unit/column can have a different number of vertical levels.
! Do I just try to recontruct the whole profile and then interpolate?
! Impossible to know which elements are 'above' and 'below' without
! finding all the elements in the first place. The vertical information
! is in the levels() array for each state vector component.

! Interpolates from state vector x to the location. It's not particularly
! happy dumping all of this straight into the model. Eventually some
! concept of a grid underlying models but above locations is going to
! be more general. May want to wait on external infrastructure projects
! for this?

type(ensemble_type),    intent(in)  :: state_handle
integer,                intent(in)  :: ens_size
type(location_type),    intent(in)  :: location
integer,                intent(in)  :: obs_kind
real(r8),               intent(out) :: expected_obs(ens_size)
integer,                intent(out) :: istatus(ens_size)

character(len=*), parameter  :: routine = 'model_interpolate'

! Local storage
character(len=obstypelength) :: qty_string
logical  :: return_now
real(r8) :: loc_array(LocationDims)

real(r8) :: llon, llat, lheight
integer  ::    istatus_liq(ens_size),    istatus_ice(ens_size)
real(r8) :: interp_val_liq(ens_size), interp_val_ice(ens_size)
integer  :: i

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the expected_obs will be set to the
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

expected_obs   = MISSING_R8 ! the DART bad value flag
interp_val_liq = MISSING_R8 ! the DART bad value flag
interp_val_ice = MISSING_R8 ! the DART bad value flag
istatus     = 0             ! will be updated with failed values
istatus_liq = 0
istatus_ice = 0

! Get the individual locations values

loc_array = get_location(location)
llon      = loc_array(1)
llat      = loc_array(2)
lheight   = loc_array(3)

call write_location(0,location,charstring=string2)
if (debug > 6) then
   write (string3,'(3(2x,f18.12))') llon, llat, lheight
   call error_handler(E_ALLMSG,routine,'requesting interpolation at',source, &
                      text2=string2,text3=string3)
endif

! Some applications just need to know the number of vertical levels.
! This is done by trying to 'interpolate' height on a large number of levels.
! When the interpolation fails, you've gone one level too far. 

if ((obs_kind == QTY_GEOPOTENTIAL_HEIGHT) .and. is_vertical(location, "LEVEL")) then
   if (nint(lheight) > nlevgrnd) then
      expected_obs = MISSING_R8
      istatus = 1
   else
      expected_obs = LEVGRND(nint(lheight))
      istatus = 0
   endif
   return ! Early Return
endif

! Standard method of interpolation.
! get_grid_vertval()       for quantities that have a vertical profile.
! compute_gridcell_value() for quantities computed for the entire gridcell.

select case( obs_kind )

   case ( QTY_SOIL_MOISTURE )

      ! only clm (history) variable 'H2OSOI' is actually soil moisture.
      ! If this is part of the state - get it, if not, construct it from
      ! QTY_SOIL_LIQUID_WATER and QTY_SOIL_ICE ... beware units.
      ! In the history file, the units are mm3/mm3
      ! In the restart file, the units are kg/m2 ... thanks ...

      call get_grid_vertval(state_handle, ens_size, location, QTY_SOIL_MOISTURE, &
                            interp_val_liq, istatus)
      if (any(istatus == 0)) then
         where(istatus == 0) expected_obs = interp_val_liq
         return
      endif

      istatus = 0

      ! TJH FIXME : make sure this is consistent with the COSMOS operator
      ! This is terrible ... the COSMOS operator wants m3/m3 ... CLM is kg/m2
      call get_grid_vertval(state_handle, ens_size, location, QTY_SOIL_LIQUID_WATER, &
                            interp_val_liq, istatus_liq)
      call track_status(ens_size, istatus_liq, interp_val_liq, istatus, return_now)
      if (return_now) return
    
      call get_grid_vertval(state_handle, ens_size, location, QTY_SOIL_ICE, &
                            interp_val_ice, istatus_ice)
      call track_status(ens_size, istatus_ice, interp_val_ice, istatus, return_now)
      if (return_now) return

      where (istatus == 0) expected_obs = interp_val_liq + interp_val_ice

   case ( QTY_SOIL_TEMPERATURE, QTY_SOIL_LIQUID_WATER, QTY_SOIL_ICE, QTY_SOIL_CARBON )

      call get_grid_vertval(state_handle, ens_size, location, obs_kind, &
                            expected_obs, istatus)

   case ( QTY_SNOWCOVER_FRAC, QTY_LEAF_AREA_INDEX, &
          QTY_LEAF_CARBON,   QTY_LIVE_STEM_CARBON,   QTY_DEAD_STEM_CARBON, &
          QTY_LEAF_NITROGEN, QTY_LIVE_STEM_NITROGEN, QTY_DEAD_STEM_NITROGEN, &
          QTY_WATER_TABLE_DEPTH, QTY_VEGETATION_TEMPERATURE, &
          QTY_FRACTION_ABSORBED_PAR, QTY_NET_CARBON_PRODUCTION, &
          QTY_PAR_DIRECT, QTY_PAR_DIFFUSE, QTY_ABSORBED_PAR, &
          QTY_SOLAR_INDUCED_FLUORESCENCE, QTY_LATENT_HEAT_FLUX)

      call compute_gridcell_value(state_handle, ens_size, location, obs_kind, &
                                  expected_obs, istatus)

   case default

      qty_string = get_name_for_quantity(obs_kind)

      write(string1,*)'not written for (integer) kind ',obs_kind
      write(string2,*)'AKA '//trim(qty_string)
      call error_handler(E_ERR,routine,string1,source,text2=string2)
      expected_obs = MISSING_R8
      istatus = 5

end select

if ((debug > 6) .and. do_output()) then
   do i = 1,ens_size
      write(     *     ,*)'ensemble member ',i,'; expected_obs value =',expected_obs(i)
      write(logfileunit,*)'ensemble member ',i,'; expected_obs value =',expected_obs(i)
   enddo
endif

! istatus is set by the calls to get_grid_vertval() or compute_gridcell_value() above.
! leave it with the value it has - don't override it here.

end subroutine model_interpolate


!------------------------------------------------------------------
!> Each gridcell may contain values for several land units, each land unit may contain
!> several columns, each column may contain several pft's.
!> aggregates across multiple pft's. So, each gridcell value
!> is an area-weighted value of an unknown number of column-based quantities.

subroutine compute_gridcell_value(state_handle, ens_size, location, qty_index, &
                                  interp_val, istatus)

! Passed variables

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location     ! location somewhere in a grid cell
integer,             intent(in)  :: qty_index    ! QTY in DART state needed for interpolation 
real(r8),            intent(out) :: interp_val(ens_size)   ! area-weighted result
integer,             intent(out) :: istatus(ens_size)      ! error code (0 == good)

character(len=*), parameter :: routine = 'compute_gridcell_value:'

! Local storage
integer  :: domain, varid, counter(ens_size)
integer(i8) :: index1, indexi, indexN 
integer  :: gridloni,gridlatj
real(r8) :: loc_lat, loc_lon
real(r8) :: state(ens_size)
real(r8) :: total(ens_size)
real(r8) :: total_area(ens_size)
real(r8), dimension(1) :: loninds,latinds
real(r8), dimension(LocationDims) :: loc
integer :: imem
character(len=obstypelength) :: varname

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

interp_val = MISSING_R8  ! the DART bad value flag
istatus    = 0

loc        = get_location(location)  ! loc is in DEGREES
loc_lon    = loc(1)
loc_lat    = loc(2)

! determine the portion of interest of the state vector
call FindVarOfInterest(qty_index, routine, domain, varid, varname)

if (varid < 1) then
   istatus = 11
   return
endif

! BOMBPROOFING - check for a vertical dimension for this variable
!>@todo check for unsupported dimensions 

! determine the grid cell for the location
latinds  = minloc(abs(LAT - loc_lat))   ! these return 'arrays' ...
loninds  = minloc(abs(LON - loc_lon))   ! these return 'arrays' ...
gridlatj = latinds(1)
gridloni = loninds(1)

if (debug > 0) then
   write(string1,*)'Working on "',trim(varname),'"'
   write(string2,*)'targetlon, lon, lon index is ',loc_lon,LON(gridloni),gridloni
   write(string3,*)'targetlat, lat, lat index is ',loc_lat,LAT(gridlatj),gridlatj
   call error_handler(E_ALLMSG,routine,string1,source,text2=string2,text3=string3)
endif

!>@todo simplify ... check for dimensionality, etc.
! If the variable is dimensioned lon,lat ... this is easy ...
!if (get_dim_name(domain,varid,1) == 'lon' .and. 
!    get_dim_name(domain,varid,2) == 'lat')
!
!endif

! Since there is no vertical component, the problem is greatly simplified.
! Area-weight an average of all pieces in the grid cell.
!>@todo  FIXME ... this is the loop that can exploit the knowledge of what 
! columnids or pftids are needed for any particular gridcell.
! gridCellInfo%pftids, gridCellInfo%columnids
! gridCellInfo(gridlatj,gridloni)%gridcell

index1 = get_index_start(domain, varid)
indexN = get_index_end(  domain, varid)

counter    = 0
total      = 0.0_r8      ! temp storage for state vector
total_area = 0.0_r8      ! temp storage for area
ELEMENTS : do indexi = index1, indexN

   if (   lonixy(indexi) /=  gridloni ) cycle ELEMENTS
   if (   latjxy(indexi) /=  gridlatj ) cycle ELEMENTS
   if ( landarea(indexi) ==   0.0_r8  ) cycle ELEMENTS

   state = get_state(indexi, state_handle)

   MEMBERS: do imem = 1, ens_size

      if(state(imem) == MISSING_R8) cycle MEMBERS

      counter(imem)    = counter(imem)    + 1
      total(imem)      = total(imem)      + state(imem)*landarea(indexi)
      total_area(imem) = total_area(imem) +             landarea(indexi)

!#    if ((debug > 99) .and. do_output()) then
!#       write(*,*)
!#       write(*,*)'gridcell location match',counter(1),'at statevector index',indexi
!#       write(*,*)'area is              (',landarea(indexi),')'
!#       write(*,*)'LON index is         (',lonixy(indexi),')'
!#       write(*,*)'LAT index is         (',latjxy(indexi),')'
!#       write(*,*)'closest LON is       (',LON(gridloni),')'
!#       write(*,*)'closest LAT is       (',LAT(gridlatj),')'
!#       write(*,*)'closest lev is       (',levels(indexi),')'
!#       write(*,*)'counter     value is (',counter(imem),')'
!#       write(*,*)'total_area  value is (',total_area(imem),')'
!#       write(*,*)'statevector value is (',state(imem),')'
!#       write(*,*)'wgtd partial  sum is (',total(imem),')'
!#    endif

   enddo MEMBERS
enddo ELEMENTS

do imem = 1,ens_size
   if (total_area(imem) > 0.0_r8 .and. istatus(imem) == 0) then
      interp_val(imem) = total(imem) / total_area(imem)
   else
      interp_val(imem) = MISSING_R8
      istatus(imem)    = 32
   endif
enddo

!# if( any(istatus == 32) ) then
!#    if (debug > 4) then
!#       write(string1, *)'Variable '//trim(varname)//' had no viable data'
!#       write(string2, *)'at gridcell ilon/jlat = (',gridloni,',',gridlatj,')'
!#       write(string3, *)'obs lon/lat = (',loc_lon,',',loc_lat,')'
!#       call error_handler(E_ALLMSG,routine,string1,source,text2=string2,text3=string3)
!#    endif
!# endif

!>@todo FIXME Need to print debugging info for any task, not just task 0
! Print more information for the really curious
!# if ((debug > 99)) then
!#    write(string1,*)'counter, total_area', counter(1), total_area(1)
!#    write(string2,*)'interp_val, istatus', interp_val(1), istatus(1)
!#    call error_handler(E_ALLMSG,routine,string1,source,text2=string2)
!# endif

end subroutine compute_gridcell_value


!------------------------------------------------------------------
!> Calculate the expected vertical value for the gridcell.
!> Each gridcell value is an area-weighted value of an unknown number
!> of column-based quantities.

subroutine get_grid_vertval(state_handle, ens_size, location, qty_index, interp_val, istatus)

type(ensemble_type), intent(in)  :: state_handle ! state vector
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location     ! location somewhere in a grid cell
integer,             intent(in)  :: qty_index    ! QTY in DART state needed for interpolation
real(r8),            intent(out) :: interp_val(ens_size)   ! area-weighted result
integer,             intent(out) :: istatus(ens_size)      ! error code (0 == good)

character(len=*), parameter :: routine = 'get_grid_vertval:'

! Local storage
integer  :: domain, ivar, index1, indexN, counter1, counter2
integer(i8) :: indexi
integer  :: gridloni,gridlatj
real(r8), dimension(LocationDims) :: loc
real(r8) :: loc_lat, loc_lon, loc_lev
real(r8) :: value_below(ens_size), value_above(ens_size), total_area(ens_size)
real(r8) :: depthbelow, depthabove
real(r8) :: topwght, botwght
real(r8), dimension(1) :: loninds,latinds

real(r8), allocatable :: above(:,:), area_above(:,:)
real(r8), allocatable :: below(:,:), area_below(:,:)
integer,  allocatable :: counter_above(:), counter_below(:)
integer :: counter
integer :: levelabove, levelbelow
integer :: imem
real(r8) :: state(ens_size)
character(len=obstypelength) :: varname
logical :: matched

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 20s

interp_val = MISSING_R8  ! the DART bad value flag
istatus    = 0           ! presume that it will work

loc        = get_location(location)  ! loc is in DEGREES
loc_lon    = loc(1)
loc_lat    = loc(2)
loc_lev    = loc(3)

!>@todo might need to relax this if we get observations IN the snow layers.
!> what happens with canopy?
if ( loc_lev < 0.0_r8 ) then
   write(string1,*)'Cannot support above-ground vertical interpolation.'
   write(string2,*)'requested a value at a depth of ',loc_lev
   write(string3,*)'CLM has negative depths to indicate above-ground values.'
   call error_handler(E_ERR,routine,string1,source,text2=string2,text3=string3)
endif

! determine the portion of interest of the state vector
call FindVarOfInterest(qty_index, routine, domain, ivar, varname)
if (ivar < 1) then
   istatus = 20
   return
endif

index1 = get_index_start(domain, varname)
indexN = get_index_end(  domain, varname)

!>@todo check applicability
!if ( variable does not have levels ) then
!   write(string1, *)'Variable '//trim(varname)//' should not use this routine.'
!   write(string2, *)'use compute_gridcell_value() instead.'
!   call error_handler(E_ERR,routine,string1,source,text2=string2)
!endif

! determine the grid cell for the location
latinds  = minloc(abs(LAT - loc_lat))   ! these return 'arrays' ...
loninds  = minloc(abs(LON - loc_lon))   ! these return 'arrays' ...
gridlatj = latinds(1)
gridloni = loninds(1)

if (debug > 4) then
   write(string1,*)'Working on "',trim(varname),'" at depth ',loc_lev
   write(string2,*)'targetlon, lon, lon index is ',loc_lon,LON(gridloni),gridloni
   write(string3,*)'targetlat, lat, lat index is ',loc_lat,LAT(gridlatj),gridlatj
   call error_handler(E_ALLMSG,routine,string1,source,text2=string2,text3=string3)
endif

! Determine the level 'above' and 'below' the desired vertical
! The above-ground 'depths' are calculated from ZISNO and are negative.
! The 'depths' are all positive numbers, increasingly positive is deeper.
! The variables currently supported use the subsurface definitions in
! the module variable LEVNGRND -- number of layers is nlevgrnd
!@todo FIXME other variables may use coordinate variables other than LEVGRND
!(despite the fact that levsoi and levdcmp use the depths defined in LEVGRND)

if (loc_lev  <= LEVGRND(1)) then  ! the top level is so close to the surface
   depthabove = LEVGRND(1)        ! just use the top level
   depthbelow = LEVGRND(1)
   levelabove = 1
   levelbelow = 1
elseif (loc_lev >= maxval(LEVGRND)) then  ! at depth, however ... do we
   depthabove    = maxval(LEVGRND)        ! fail or just use the deepest
   depthbelow    = maxval(LEVGRND)        ! I am using the deepest.
   levelabove    = nlevgrnd
   levelbelow    = nlevgrnd
else
   LAYERS : do indexi = 2,nlevgrnd
      if (loc_lev   < LEVGRND(indexi)) then
         levelabove = indexi-1
         levelbelow = indexi
         depthabove = LEVGRND(levelabove)
         depthbelow = LEVGRND(levelbelow)
         exit LAYERS
      endif
   enddo LAYERS
endif

if (debug > 4) then
   write(string1,*)'..  depth       above ',depthabove, 'level above', levelabove
   write(string2,*)'depth of interest ',loc_lev
   write(string3,*)'depth       below ',depthbelow, 'level below', levelbelow
   call error_handler(E_ALLMSG,routine,string1,source,text2=string2,text3=string3)
endif

! Determine how many elements can contribute to the gridcell value.
! There are multiple column-based contributors, each column has a
! separate area-based weight. There are multiple levels.
! I believe I have to keep track of all of them to sort out how to
! calculate the gridcell value at a particular depth.
! FIXME: can we leverage the gridCellInfo structure?

! These are just to get the max size
counter1 = 0
counter2 = 0

GRIDCELL : do indexi = index1, indexN

   if ( lonixy(indexi) /=  gridloni )  cycle GRIDCELL
   if ( latjxy(indexi) /=  gridlatj )  cycle GRIDCELL

   if (levels(indexi) == depthabove) counter1 = counter1 + 1
   if (levels(indexi) == depthbelow) counter2 = counter2 + 1

enddo GRIDCELL

counter = max(counter1, counter2)

if (debug > 1) then
   write(string1, *)'"'//trim(varname)//'" had ',counter,' viable data'
   write(string2, *)'at gridcell lon/lat = (',gridloni,',',gridlatj,')'
   call write_location(0,location,charstring=string3)
   call error_handler(E_ALLMSG,routine,string1,source,text2=string2,text3=string3)
endif

if ( counter == 0 ) then
   istatus = 21
   return
endif

allocate(counter_above(ens_size), above(ens_size,counter), area_above(ens_size,counter), &
         counter_below(ens_size), below(ens_size,counter), area_below(ens_size,counter))

counter_above = 0
counter_below = 0
above         = 0.0_r8
below         = 0.0_r8
area_below    = 0.0_r8
area_above    = 0.0_r8
value_above   = MISSING_R8
value_below   = MISSING_R8

ELEMENTS : do indexi = index1, indexN

   if ( lonixy(indexi) /=  gridloni )  cycle ELEMENTS
   if ( latjxy(indexi) /=  gridlatj )  cycle ELEMENTS
   if ( levels(indexi) /=  depthabove .and. &
        levels(indexi) /=  depthbelow) cycle ELEMENTS

   state = get_state(indexi, state_handle)

   matched = .false.

   MEMBERS : do imem = 1,ens_size

      if (state(imem) == MISSING_R8) cycle MEMBERS

      ! sometimes it could be right ON a level, so it could 
      ! match both depths

      if (levels(indexi) == depthabove) then
         counter_above(imem)             = counter_above(imem) + 1
         above(     imem, counter_above(imem)) = state(imem)
         area_above(imem, counter_above(imem)) = landarea(indexi)
         matched = .true.
      endif

      if(levels(indexi) == depthbelow) then
         counter_below(imem)             = counter_below(imem) + 1
         below(     imem, counter_below(imem)) = state(imem)
         area_below(imem, counter_below(imem)) = landarea(indexi)
         matched = .true.
      endif

      if (debug > 99 .and. do_output() .and. matched) then
         write(*,*)routine,'Variable "'//trim(varname)//'" matched ',matched
         write(*,*)routine,'gridcell location match at statevector index',indexi
         write(*,*)routine,'area is ',landarea(indexi), ' depth is  ',levels(indexi)
         write(*,*)routine,'LON index is     ',lonixy(indexi),' LON is ',LON(gridloni)
         write(*,*)routine,'LAT index is     ',latjxy(indexi),' LAT is ',LAT(gridlatj)
         write(*,*)routine,'state value is (',state(imem),')'
      endif

   enddo MEMBERS

enddo ELEMENTS

! could arise if the above or below was 'missing' ... but the mate was not.

if ( any(counter_above /= counter_below) ) then
   write(string1, *)'Variable "'//trim(varname)//'" has peculiar interpolation problems.'
   write(string2, *)'uneven number of values "above" and "below"'
   write(string3, *)'counter_above == ',counter_above,' /= ',counter_below,' == counter_below'
   call error_handler(E_ERR,routine,string1,source,text2=string2,text3=string3)
   deallocate(counter_above, counter_below, above, below, area_above, area_below)
   istatus = 22
   return
endif

do imem = 1, ens_size

   ! Determine the area using the level above the depth of interest.
   total_area(imem) = sum(area_above(imem, :))

   if ( total_area(imem) /= 0.0_r8 .and. istatus(imem) == 0) then
      ! normalize the area-based weights
      area_above(imem, :) = area_above(imem, :) / total_area(imem)
      value_above(imem) = sum(above(imem, :) * area_above(imem, :))
   else
      if (debug > 1) then
         write(string1, *)'Variable '//trim(varname)//' had no viable data above'
         write(string2, *)'at gridcell lon/lat/level = (',gridloni,',',gridlatj,',',levelabove,')'
         call write_location(0,location,charstring=string3)
         call error_handler(E_ALLMSG,routine,string1,source,text2=string2,text3=string3)
      endif
   endif

   ! Determine the value for the level below the depth of interest.
   total_area(imem) = sum(area_below(imem, :))

   if ( total_area(imem) /= 0.0_r8 .and. istatus(imem) == 0 ) then
      ! normalize the area-based weights
      area_below(imem, :) = area_below(imem, :) / total_area(imem)
      value_below(imem) = sum(below(imem, :) * area_below(imem, :))
   else
      if (debug > 1) then
         write(string1, *)'Variable '//trim(varname)//' had no viable data below'
         write(string2, *)'at gridcell lon/lat/lev = (',gridloni,',',gridlatj,',',levelbelow,')'
         call write_location(0,location,charstring=string3)
         call error_handler(E_ALLMSG,routine,string1,source,text2=string2,text3=string3)
      endif
   endif

enddo

if (depthbelow == depthabove) then
   topwght = 1.0_r8
   botwght = 0.0_r8
else
   topwght = (depthbelow - loc_lev) / (depthbelow - depthabove)
   botwght = (loc_lev - depthabove) / (depthbelow - depthabove)
endif

!>@todo The istatus should be poor if the state was missing to begin with
do imem = 1, ens_size
   if (istatus(imem)     == 0          .and. &
       value_above(imem) /= MISSING_R8 .and. &
       value_below(imem) /= MISSING_R8) then
      interp_val(imem) = value_above(imem)*topwght + value_below(imem)*botwght
   else
      istatus(imem) = 23
      interp_val(imem) = MISSING_R8
   endif
enddo

deallocate(counter_above, counter_below, above, below, area_above, area_below)

end subroutine get_grid_vertval


!------------------------------------------------------------------
!> Read the dimensions from the history netcdf file.

subroutine get_history_dims(ncid, fname, cstat, lon, lat, levgrnd, levsoi, levdcmp, &
                            lonatm, latatm, lonrof, latrof)

integer,           intent(inout) :: ncid
character(len=*),  intent(in)    :: fname
character(len=*),  intent(in)    :: cstat ! how do you want to leave the netcdf file
integer,           intent(out)   :: lon
integer,           intent(out)   :: lat
integer,           intent(out)   :: levgrnd
integer,           intent(out)   :: levsoi
integer,           intent(out)   :: levdcmp
integer, OPTIONAL, intent(out)   :: lonatm, latatm
integer, OPTIONAL, intent(out)   :: lonrof, latrof

character(len=*), parameter :: routine = 'get_history_dims'

integer :: dimid

! get the ball rolling ...

if (ncid == 0) ncid = nc_open_file_readonly(fname, routine)

! The SingleColumn (and unstructured grid) configurations
! do not have a 'lon' and 'lat' dimension. There is only 'lndgrid'

if ( nf90_inq_dimid(ncid, 'lndgrid', dimid) == NF90_NOERR ) unstructured = .true.

if (unstructured) then ! use the lndgrid dimension for both lon and lat

   lon = nc_get_dimension_size(ncid,'lndgrid',routine)
   lat = lon

else

   lon = nc_get_dimension_size(ncid,'lon',routine)
   lat = nc_get_dimension_size(ncid,'lat',routine)

endif

levgrnd = nc_get_dimension_size(ncid,'levgrnd',routine)
levsoi  = nc_get_dimension_size(ncid,'levsoi', routine)
levdcmp = nc_get_dimension_size(ncid,'levdcmp',routine)

if (present(lonatm)) lonatm = nc_get_dimension_size(ncid,'lonatm',routine)
if (present(latatm)) latatm = nc_get_dimension_size(ncid,'latatm',routine)
if (present(lonrof)) lonrof = nc_get_dimension_size(ncid,'lonrof',routine)
if (present(latrof)) latrof = nc_get_dimension_size(ncid,'latrof',routine)

if ((debug > 8) .and. do_output()) then
   write(logfileunit,*)
   write(logfileunit,*)'get_history_dims output follows:'
   write(logfileunit,*)'nlon = ',nlon
   write(logfileunit,*)'nlat = ',nlat
   write(     *     ,*)
   write(     *     ,*)'get_history_dims output follows:'
   write(     *     ,*)'nlon = ',nlon
   write(     *     ,*)'nlat = ',nlat

   if (present(lonatm)) write(logfileunit,*)'lonatm   = ',lonatm
   if (present(latatm)) write(logfileunit,*)'latatm   = ',latatm
   if (present(lonrof)) write(logfileunit,*)'lonrof   = ',lonrof
   if (present(latrof)) write(logfileunit,*)'latrof   = ',latrof
   if (present(lonatm)) write(     *     ,*)'lonatm   = ',lonatm
   if (present(latatm)) write(     *     ,*)'latatm   = ',latatm
   if (present(lonrof)) write(     *     ,*)'lonrof   = ',lonrof
   if (present(latrof)) write(     *     ,*)'latrof   = ',latrof
endif

if (cstat == 'close') then
   call nc_close_file(ncid, routine, fname)
   ncid = 0
endif

end subroutine get_history_dims


!------------------------------------------------------------------
!> Read the grid dimensions from the CLM history netcdf file.
!> LON,LAT,AREA,LANDFRAC,LEVGRN all have module scope

subroutine get_full_grid(ncid, fname, cstat)

integer,                  intent(inout) :: ncid
character(len=*),         intent(in)    :: fname
character(len=*),         intent(in)    :: cstat

character(len=*), parameter :: routine = 'get_full_grid'

! Make sure the variables are the right size ...
! at some point in the future ...

if (ncid == 0) ncid = nc_open_file_readonly(fname, routine)

! The lat/lon matrices in the history file have been masked by
! the land values such that the wet cells are 'missing' values.
! This makes it less than useful for us. Thankfully, the 1D
! lat/lon arrays have no such mask applied. We use these.

call nc_get_variable(ncid, 'lon'    , LON,     routine)
call nc_get_variable(ncid, 'lat'    , LAT,     routine)
call nc_get_variable(ncid, 'levgrnd', LEVGRND, routine)

if ( nc_variable_exists(ncid,'levsoi') ) then
   call nc_get_variable(ncid,'levsoi', LEVSOI, routine)
else
   LEVSOI(1:nlevsoi) = LEVGRND(1:nlevsoi)
endif
    
if ( nc_variable_exists(ncid,'levdcmp') ) then
   call nc_get_variable(ncid,'levdcmp', LEVDCMP, routine)
else
   LEVDCMP(1:nlevdcmp) = LEVGRND(1:nlevdcmp)
endif


if (unstructured) then
   call nc_get_variable(ncid,'area'    ,AREA1D,     routine)
   call nc_get_variable(ncid,'landfrac',LANDFRAC1D, routine)
   where(AREA1D     == MISSING_R8) AREA1D     = 0.0_r8
   where(LANDFRAC1D == MISSING_R8) LANDFRAC1D = 0.0_r8
else
   call nc_get_variable(ncid,'area'    ,AREA2D,     routine)
   call nc_get_variable(ncid,'landfrac',LANDFRAC2D, routine)
   where(AREA2D     == MISSING_R8) AREA2D     = 0.0_r8
   where(LANDFRAC2D == MISSING_R8) LANDFRAC2D = 0.0_r8
endif

! just to make sure we are [0,360] and [-90,90]

where (LON <   0.0_r8) LON = LON + 360.0_r8
where (LON > 360.0_r8) LON = LON - 360.0_r8

if (any(LON < 0.0_r8)) then
   write(string1,*)'longitudes in history file variable "lon" still negative.'
   call error_handler(E_ERR,routine,string1,source)
endif

where (LAT < -90.0_r8) LAT = -90.0_r8
where (LAT >  90.0_r8) LAT =  90.0_r8

if (cstat == 'close') then
   call nc_close_file(ncid, routine, fname)
   ncid = 0
endif

! A little sanity check

if ((debug > 7) .and. do_output()) then

   write(logfileunit,*)
   write(logfileunit,*)'history_file grid information as interpreted ...'
   write(logfileunit,*)'lon      range ',minval(LON     ),maxval(LON     )
   write(logfileunit,*)'lat      range ',minval(LAT     ),maxval(LAT     )
   write(logfileunit,*)'levgrnd  range ',minval(LEVGRND ),maxval(LEVGRND )
   write(logfileunit,*)'levgrnd  is ',LEVGRND
   write(     *     ,*)
   write(     *     ,*)'history_file grid information as interpreted ...'
   write(     *     ,*)'lon      range ',minval(LON     ),maxval(LON     )
   write(     *     ,*)'lat      range ',minval(LAT     ),maxval(LAT     )
   write(     *     ,*)'levgrnd  range ',minval(LEVGRND ),maxval(LEVGRND )
   write(     *     ,*)'levgrnd  is ',LEVGRND

   if (unstructured) then
      write(logfileunit,*)'area     range ',minval(AREA1D    ),maxval(AREA1D    )
      write(logfileunit,*)'landfrac range ',minval(LANDFRAC1D),maxval(LANDFRAC1D)
      write(     *     ,*)'area     range ',minval(AREA1D    ),maxval(AREA1D    )
      write(     *     ,*)'landfrac range ',minval(LANDFRAC1D),maxval(LANDFRAC1D)
   else
      write(logfileunit,*)'area     range ',minval(AREA2D    ),maxval(AREA2D    )
      write(logfileunit,*)'landfrac range ',minval(LANDFRAC2D),maxval(LANDFRAC2D)
      write(     *     ,*)'area     range ',minval(AREA2D    ),maxval(AREA2D    )
      write(     *     ,*)'landfrac range ',minval(LANDFRAC2D),maxval(LANDFRAC2D)
   endif

endif

end subroutine get_full_grid


!------------------------------------------------------------------
! From "main/landunit_varcon.F90" svn rev 76577 
!
! integer, parameter, public :: istsoil    = 1  !soil         landunit type (natural vegetation)
! integer, parameter, public :: istcrop    = 2  !crop         landunit type
! integer, parameter, public :: istice     = 3  !land ice     landunit type (glacier)
! integer, parameter, public :: istice_mec = 4  !land ice (multiple elevation classes) landunit type
! integer, parameter, public :: istdlak    = 5  !deep lake    landunit type (now used for all lakes)
! integer, parameter, public :: istwet     = 6  !wetland      landunit type (swamp, marsh, etc.)
! 
! integer, parameter, public :: isturb_tbd = 7  !urban tbd    landunit type
! integer, parameter, public :: isturb_hd  = 8  !urban hd     landunit type
! integer, parameter, public :: isturb_md  = 9  !urban md     landunit type
!
! Were these from CLM 4.0?
! landunit types
! 1 soil (natural vegation/bare ground)
! 2 glacier
! 3 lake
! 4 not used (shallow lake)
! 5 wetland
! 6 urban
! 7 ice (new glacier model)
! 8 crop (if using crop model)

subroutine get_clm_codes(ncid, fname, cstat)

integer,          intent(inout) :: ncid
character(len=*), intent(in)    :: fname
character(len=*), intent(in)    :: cstat

character(len=*), parameter :: routine = 'get_clm_codes'
integer :: ret

if (ncid == 0) ncid = nc_open_file_readonly(fname, routine)

! If the read is successful, use it - if not, use the default from CLM 4.5

! Codes for column-level variables

ret = nf90_get_att(ncid, NF90_GLOBAL, 'icol_vegetated_or_bare_soil', icol_vegetated_or_bare_soil)
if ( ret /= NF90_NOERR ) icol_vegetated_or_bare_soil = 1

ret = nf90_get_att(ncid, NF90_GLOBAL, 'icol_crop', icol_crop)
if ( ret /= NF90_NOERR ) icol_crop = 2

! ret = nf90_get_att(ncid, NF90_GLOBAL, 'icol_crop_noncompete', icol_crop_noncompete)
! if ( ret /= NF90_NOERR ) icol_crop_noncompete = "no idea what this is"

ret = nf90_get_att(ncid, NF90_GLOBAL, 'icol_landice', icol_landice)
if ( ret /= NF90_NOERR ) icol_landice = 3

! ret = nf90_get_att(ncid, NF90_GLOBAL, 'icol_landice_multiple_elevation_classes', icol_landice_multiple_elevation_classes)
! if ( ret /= NF90_NOERR ) icol_landice_multiple_elevation_classes = "no idea what this is"

ret = nf90_get_att(ncid, NF90_GLOBAL, 'icol_deep_lake', icol_deep_lake)
if ( ret /= NF90_NOERR ) icol_deep_lake = 5

ret = nf90_get_att(ncid, NF90_GLOBAL, 'icol_wetland', icol_wetland)
if ( ret /= NF90_NOERR ) icol_wetland = 6

ret = nf90_get_att(ncid, NF90_GLOBAL, 'icol_urban_roof', icol_urban_roof)
if ( ret /= NF90_NOERR ) icol_urban_roof = 71

ret = nf90_get_att(ncid, NF90_GLOBAL, 'icol_urban_sunwall', icol_urban_sunwall)
if ( ret /= NF90_NOERR ) icol_urban_sunwall = 72

ret = nf90_get_att(ncid, NF90_GLOBAL, 'icol_urban_shadewall', icol_urban_shadewall)
if ( ret /= NF90_NOERR ) icol_urban_shadewall = 73

ret = nf90_get_att(ncid, NF90_GLOBAL, 'icol_urban_impervious_road', icol_urban_impervious_road)
if ( ret /= NF90_NOERR ) icol_urban_impervious_road = 74

ret = nf90_get_att(ncid, NF90_GLOBAL, 'icol_urban_pervious_road', icol_urban_pervious_road)
if ( ret /= NF90_NOERR ) icol_urban_pervious_road = 75

! Codes for landunit-level variables

ret = nf90_get_att(ncid, NF90_GLOBAL, 'ilun_vegetated_or_bare_soil', ilun_vegetated_or_bare_soil)
if ( ret /= NF90_NOERR ) ilun_vegetated_or_bare_soil = 1

ret = nf90_get_att(ncid, NF90_GLOBAL, 'ilun_crop', ilun_crop)
if ( ret /= NF90_NOERR ) ilun_crop = 2

ret = nf90_get_att(ncid, NF90_GLOBAL, 'ilun_UNUSED', ilun_UNUSED)
if ( ret /= NF90_NOERR ) ilun_UNUSED = 3

ret = nf90_get_att(ncid, NF90_GLOBAL, 'ilun_landice_multiple_elevation_classes', ilun_landice)
if ( ret /= NF90_NOERR ) ilun_landice = 4

ret = nf90_get_att(ncid, NF90_GLOBAL, 'ilun_deep_lake', ilun_deep_lake)
if ( ret /= NF90_NOERR ) ilun_deep_lake = 5

ret = nf90_get_att(ncid, NF90_GLOBAL, 'ilun_wetland', ilun_wetland)
if ( ret /= NF90_NOERR ) ilun_wetland = 6

ret = nf90_get_att(ncid, NF90_GLOBAL, 'ilun_urban_tbd', ilun_urban_tbd)
if ( ret /= NF90_NOERR ) ilun_urban_tbd = 7

ret = nf90_get_att(ncid, NF90_GLOBAL, 'ilun_urban_hd', ilun_urban_hd)
if ( ret /= NF90_NOERR ) ilun_urban_hd = 8

ret = nf90_get_att(ncid, NF90_GLOBAL, 'ilun_urban_md', ilun_urban_md)
if ( ret /= NF90_NOERR ) ilun_urban_md = 9

! Summary

if ((debug > 1) .and. do_output()) then

   call say('Reporting column codes: ')
   write(string1,*) 'icol_vegetated_or_bare_soil = ', icol_vegetated_or_bare_soil
   call say(string1)
   write(string1,*) 'icol_crop = ', icol_crop
   call say(string1)
   write(string1,*) 'icol_landice = ', icol_landice
   call say(string1)
   write(string1,*) 'icol_deep_lake = ', icol_deep_lake
   call say(string1)
   write(string1,*) 'icol_wetland = ', icol_wetland
   call say(string1)
   write(string1,*) 'icol_urban_roof = ', icol_urban_roof
   call say(string1)
   write(string1,*) 'icol_urban_sunwall = ', icol_urban_sunwall
   call say(string1)
   write(string1,*) 'icol_urban_shadewall = ', icol_urban_shadewall
   call say(string1)
   write(string1,*) 'icol_urban_impervious_road = ', icol_urban_impervious_road
   call say(string1)
   write(string1,*) 'icol_urban_pervious_road = ', icol_urban_pervious_road
   call say(string1)

   call say('Reporting landunit codes: ')
   write(string1,*) 'ilun_vegetated_or_bare_soil = ',ilun_vegetated_or_bare_soil
   call say(string1)
   write(string1,*) 'ilun_crop = ',ilun_crop
   call say(string1)
   write(string1,*) 'ilun_UNUSED = ',ilun_UNUSED
   call say(string1)
   write(string1,*) 'ilun_landice_multiple_elevation_classes = ',ilun_landice
   call say(string1)
   write(string1,*) 'ilun_deep_lake = ',ilun_deep_lake
   call say(string1)
   write(string1,*) 'ilun_wetland = ',ilun_wetland
   call say(string1)
   write(string1,*) 'ilun_urban_tbd = ',ilun_urban_tbd
   call say(string1)
   write(string1,*) 'ilun_urban_hd = ',ilun_urban_hd
   call say(string1)
   write(string1,*) 'ilun_urban_md = ',ilun_urban_md
   call say(string1)

endif

end subroutine get_clm_codes


!------------------------------------------------------------------
!> Read the dimensions from the CLM restart netcdf file.


subroutine get_sparse_dims(ncid, fname, cstat)

integer,          intent(inout) :: ncid
character(len=*), intent(in)    :: fname
character(len=*), intent(in)    :: cstat

character(len=*), parameter :: routine = 'get_sparse_dims'

integer :: mylevgrnd

if (ncid == 0) ncid = nc_open_file_readonly(fname, routine)

ngridcell = nc_get_dimension_size(ncid, 'gridcell', routine)
nlandunit = nc_get_dimension_size(ncid, 'landunit', routine)
ncolumn   = nc_get_dimension_size(ncid, 'column',   routine)
npft      = nc_get_dimension_size(ncid, 'pft',      routine)
mylevgrnd = nc_get_dimension_size(ncid, 'levgrnd',  routine)

if (mylevgrnd /= nlevgrnd) then
   write(string1,*)'Number of ground levels in restart file is',mylevgrnd
   write(string2,*)'Number of ground levels in history file is', nlevgrnd
   call error_handler(E_ERR,routine,string1,source,text2=string2)
endif

nlevlak  = nc_get_dimension_size(ncid, 'levlak',  routine)
nlevtot  = nc_get_dimension_size(ncid, 'levtot',  routine)
nnumrad  = nc_get_dimension_size(ncid, 'numrad',  routine)
nlevcan  = nc_get_dimension_size(ncid, 'levcan',  routine)
nlevsno  = nc_get_dimension_size(ncid, 'levsno',  routine)
nlevsno1 = nc_get_dimension_size(ncid, 'levsno1', routine)

! CLM4 does not have a multi-level canopy.
! CLM4.5 has a multi-level canopy.

if (cstat == 'close') then
   call nc_close_file(ncid, routine, fname)
   ncid = 0
endif

! Echo what we know if desired.
if ((debug > 1) .and. do_output()) then
   write(logfileunit,*)
   write(logfileunit,*)'get_sparse_dims output follows:'
   write(logfileunit,*)'ngridcell = ',ngridcell
   write(logfileunit,*)'nlandunit = ',nlandunit
   write(logfileunit,*)'ncolumn   = ',ncolumn
   write(logfileunit,*)'npft      = ',npft
   write(logfileunit,*)'nlevgrnd  = ',nlevgrnd
   write(logfileunit,*)'nlevlak   = ',nlevlak
   write(logfileunit,*)'nlevsno   = ',nlevsno
   write(logfileunit,*)'nlevsno1  = ',nlevsno1
   write(logfileunit,*)'nlevtot   = ',nlevtot
   write(logfileunit,*)'nnumrad   = ',nnumrad
   write(logfileunit,*)'nlevcan   = ',nlevcan
   write(     *     ,*)
   write(     *     ,*)'get_sparse_dims output follows:'
   write(     *     ,*)'ngridcell = ',ngridcell
   write(     *     ,*)'nlandunit = ',nlandunit
   write(     *     ,*)'ncolumn   = ',ncolumn
   write(     *     ,*)'npft      = ',npft
   write(     *     ,*)'nlevgrnd  = ',nlevgrnd
   write(     *     ,*)'nlevlak   = ',nlevlak
   write(     *     ,*)'nlevsno   = ',nlevsno
   write(     *     ,*)'nlevsno1  = ',nlevsno1
   write(     *     ,*)'nlevtot   = ',nlevtot
   write(     *     ,*)'nnumrad   = ',nnumrad
   write(     *     ,*)'nlevcan   = ',nlevcan
endif

end subroutine get_sparse_dims


!------------------------------------------------------------------



subroutine get_sparse_geog(ncid, fname, cstat)
!------------------------------------------------------------------
!
! Read the geography information from from the restart netcdf file.

integer,          intent(inout) :: ncid
character(len=*), intent(in)    :: fname
character(len=*), intent(in)    :: cstat

character(len=*), parameter :: routine = 'get_sparse_geog'

real(r8), allocatable :: temp2d(:,:)

if (ncid == 0) ncid = nc_open_file_readonly(fname, routine)

! Make sure the variables are the right size ...
! by comparing agains the size of the variable ...

if ( ngridcell < 0 ) then
   write(string1,*)'Unable to read the number of gridcells.'
   call error_handler(E_ERR,routine,string1,source)
endif

if ( nlandunit < 0 ) then
   write(string1,*)'Unable to read the number of land units.'
   call error_handler(E_ERR,routine,string1,source)
endif

if ( ncolumn < 0 ) then
   write(string1,*)'Unable to read the number of columns.'
   call error_handler(E_ERR,routine,string1,source)
endif

if ( npft < 0 ) then
   write(string1,*)'Unable to read the number of pfts.'
   call error_handler(E_ERR,routine,string1,source)
endif

! Read the netcdf file data

call nc_get_variable(ncid, 'grid1d_ixy',     grid1d_ixy,     routine)
call nc_get_variable(ncid, 'grid1d_jxy',     grid1d_jxy,     routine)
call nc_get_variable(ncid, 'land1d_ixy',     land1d_ixy,     routine)
call nc_get_variable(ncid, 'land1d_jxy',     land1d_jxy,     routine)
call nc_get_variable(ncid, 'land1d_wtxy',    land1d_wtxy,    routine)
call nc_get_variable(ncid, 'land1d_ityplun', land1d_ityplun, routine)
call nc_get_variable(ncid, 'cols1d_ixy',     cols1d_ixy,     routine)
call nc_get_variable(ncid, 'cols1d_jxy',     cols1d_jxy,     routine)
call nc_get_variable(ncid, 'cols1d_wtxy',    cols1d_wtxy,    routine)
call nc_get_variable(ncid, 'cols1d_lon',     cols1d_lon,    routine)
call nc_get_variable(ncid, 'cols1d_lat',     cols1d_lat,    routine)
call nc_get_variable(ncid, 'cols1d_ityplun', cols1d_ityplun, routine)
call nc_get_variable(ncid, 'pfts1d_ixy',     pfts1d_ixy,     routine)
call nc_get_variable(ncid, 'pfts1d_jxy',     pfts1d_jxy,     routine)
call nc_get_variable(ncid, 'pfts1d_wtxy',    pfts1d_wtxy,    routine)
call nc_get_variable(ncid, 'pfts1d_lon',     pfts1d_lon,    routine)
call nc_get_variable(ncid, 'pfts1d_lat',     pfts1d_lat,    routine)
call nc_get_variable(ncid, 'pfts1d_ityplun', pfts1d_ityplun, routine)

! zsno is NOT optional ... so it IS a fatal error if it is not present (for now, anyway).
! as read into fortran ... zsno(1,:) is the level closest to the sun.
! as read into fortran ... zsno(N,:) is the level closest to the ground.
!
! ZSNO and ZISNO are dimensioned identically. 
! Here is example when SNLSNO = -10:
!            DZSNO     ZSNO                  ZISNO    
! ( 1,11345) 0,        0,                    0,
! ( 2,11345) 0,        0,                    0,
! ( 3,11345) 0.02,     -22.9724115403977,    -22.9824115403977
! ( 4,11345) 0.05,     -22.9374115403977,    -22.9624115403977
! ( 5,11345) 0.11,     -22.8574115403977,    -22.9124115403977
! ( 6,11345) 0.23,     -22.6874115403977,    -22.8024115403977
! ( 7,11345) 0.47,     -22.3374115403977,    -22.5724115403977
! ( 8,11345) 0.95,     -21.6274115403977,    -22.1024115403977
! ( 9,11345) 1.91,     -20.1974115403977,    -21.1524115403977
! (10,11345) 3.83,     -17.3274115403977,    -19.2424115403977
! (11,11345) 7.67,     -11.5774115403977,    -15.4124115403977
! (12,11345) 7.74+,    -3.87120577019883,    -7.74241154039766

if (nlevsno > 0 ) then
   call nc_get_variable(ncid, 'ZSNO',  zsno,  routine)
   allocate(temp2d(nlevsno,ncolumn))
   call nc_get_variable(ncid, 'ZISNO', temp2D, routine)
   zisno(1:nlevsno,:) = temp2D;
   zisno( nlevsno1,:) = 0.0_r8;
   deallocate(temp2D)
else
   write(string1,*) 'levsno must be in restart file'
   call error_handler(E_ERR,routine,string1,source)
endif

if (cstat == 'close') then
   call nc_close_file(ncid, routine, fname)
   ncid = 0
endif

! A little sanity check

if ((debug > 7) .and. do_output()) then

   write(logfileunit,*)
   write(logfileunit,*)'Raw lat/lon information as read ...'
   write(logfileunit,*)'grid1d_ixy     range ',minval(grid1d_ixy),    maxval(grid1d_ixy)
   write(logfileunit,*)'grid1d_jxy     range ',minval(grid1d_jxy),    maxval(grid1d_jxy)

   write(logfileunit,*)'land1d_ixy     range ',minval(land1d_ixy),    maxval(land1d_ixy)
   write(logfileunit,*)'land1d_jxy     range ',minval(land1d_jxy),    maxval(land1d_jxy)
   write(logfileunit,*)'land1d_wtxy    range ',minval(land1d_wtxy),   maxval(land1d_wtxy)
   write(logfileunit,*)'land1d_ityplun range ',minval(land1d_ityplun),maxval(land1d_ityplun)

   write(logfileunit,*)'cols1d_ixy     range ',minval(cols1d_ixy),    maxval(cols1d_ixy)
   write(logfileunit,*)'cols1d_jxy     range ',minval(cols1d_jxy),    maxval(cols1d_jxy)
   write(logfileunit,*)'cols1d_wtxy    range ',minval(cols1d_wtxy),   maxval(cols1d_wtxy)
   write(logfileunit,*)'cols1d_lon     range ',minval(cols1d_lon),    maxval(cols1d_lon)
   write(logfileunit,*)'cols1d_lat     range ',minval(cols1d_lat),    maxval(cols1d_lat)
   write(logfileunit,*)'cols1d_ityplun range ',minval(cols1d_ityplun),maxval(cols1d_ityplun)

   write(logfileunit,*)'pfts1d_ixy     range ',minval(pfts1d_ixy),    maxval(pfts1d_ixy)
   write(logfileunit,*)'pfts1d_jxy     range ',minval(pfts1d_jxy),    maxval(pfts1d_jxy)
   write(logfileunit,*)'pfts1d_wtxy    range ',minval(pfts1d_wtxy),   maxval(pfts1d_wtxy)
   write(logfileunit,*)'pfts1d_lon     range ',minval(pfts1d_lon),    maxval(pfts1d_lon)
   write(logfileunit,*)'pfts1d_lat     range ',minval(pfts1d_lat),    maxval(pfts1d_lat)
   write(logfileunit,*)'pfts1d_ityplun range ',minval(pfts1d_ityplun),maxval(pfts1d_ityplun)

   if (nlevsno > 0) write(logfileunit,*)'zsno           range ',minval(zsno), maxval(zsno)
   if (nlevsno > 0) write(logfileunit,*)'zisno          range ',minval(zisno),maxval(zisno)

   write(     *     ,*)
   write(     *     ,*)'Raw lat/lon information as read ...'
   write(     *     ,*)'grid1d_ixy     range ',minval(grid1d_ixy),    maxval(grid1d_ixy)
   write(     *     ,*)'grid1d_jxy     range ',minval(grid1d_jxy),    maxval(grid1d_jxy)

   write(     *     ,*)'land1d_ixy     range ',minval(land1d_ixy),    maxval(land1d_ixy)
   write(     *     ,*)'land1d_jxy     range ',minval(land1d_jxy),    maxval(land1d_jxy)
   write(     *     ,*)'land1d_wtxy    range ',minval(land1d_wtxy),   maxval(land1d_wtxy)
   write(     *     ,*)'land1d_ityplun range ',minval(land1d_ityplun),maxval(land1d_ityplun)

   write(     *     ,*)'cols1d_ixy     range ',minval(cols1d_ixy),    maxval(cols1d_ixy)
   write(     *     ,*)'cols1d_jxy     range ',minval(cols1d_jxy),    maxval(cols1d_jxy)
   write(     *     ,*)'cols1d_wtxy    range ',minval(cols1d_wtxy),   maxval(cols1d_wtxy)
   write(     *     ,*)'cols1d_lon     range ',minval(cols1d_lon),    maxval(cols1d_lon)
   write(     *     ,*)'cols1d_lat     range ',minval(cols1d_lat),    maxval(cols1d_lat)
   write(     *     ,*)'cols1d_ityplun range ',minval(cols1d_ityplun),maxval(cols1d_ityplun)

   write(     *     ,*)'pfts1d_ixy     range ',minval(pfts1d_ixy),    maxval(pfts1d_ixy)
   write(     *     ,*)'pfts1d_jxy     range ',minval(pfts1d_jxy),    maxval(pfts1d_jxy)
   write(     *     ,*)'pfts1d_wtxy    range ',minval(pfts1d_wtxy),   maxval(pfts1d_wtxy)
   write(     *     ,*)'pfts1d_lon     range ',minval(pfts1d_lon),    maxval(pfts1d_lon)
   write(     *     ,*)'pfts1d_lat     range ',minval(pfts1d_lat),    maxval(pfts1d_lat)
   write(     *     ,*)'pfts1d_ityplun range ',minval(pfts1d_ityplun),maxval(pfts1d_ityplun)

   if (nlevsno > 0) write(     *     ,*)'zsno           range ',minval(zsno), maxval(zsno)
   if (nlevsno > 0) write(     *     ,*)'zisno          range ',minval(zisno),maxval(zisno)

endif

end subroutine get_sparse_geog


!------------------------------------------------------------------
!> The restart netcdf files have the time of the state.

function get_state_time_ncid( ncid )

type(time_type) :: get_state_time_ncid
integer, intent(in) :: ncid

character(len=*), parameter :: routine = 'get_state_time_ncid'
integer :: io, VarID
integer :: rst_curr_ymd, rst_curr_tod, leftover
integer :: year, month, day, hour, minute, second

if ( .not. module_initialized ) call static_init_model

io = nf90_inq_varid(ncid, 'timemgr_rst_curr_ymd', VarID)
call nc_check(io, routine, 'inq_varid timemgr_rst_curr_ymd', ncid=ncid)

io = nf90_get_var(ncid, VarID, rst_curr_ymd)
call nc_check(io, routine, 'get_var rst_curr_ymd', ncid=ncid)

io = nf90_inq_varid(ncid, 'timemgr_rst_curr_tod', VarID)
call nc_check(io, routine, 'inq_varid timemgr_rst_curr_tod', ncid=ncid)

io = nf90_get_var(ncid, VarID, rst_curr_tod)
call nc_check(io, routine, 'get_var rst_curr_tod', ncid=ncid)

year     = rst_curr_ymd/10000
leftover = rst_curr_ymd - year*10000
month    = leftover/100
day      = leftover - month*100

hour     = rst_curr_tod/3600
leftover = rst_curr_tod - hour*3600
minute   = leftover/60
second   = leftover - minute*60

get_state_time_ncid = set_date(year, month, day, hour, minute, second)

end function get_state_time_ncid


!------------------------------------------------------------------
!> This defines the window used for assimilation.
!> all observations +/- half this timestep are assimilated.

function set_model_time_step()

type(time_type) :: set_model_time_step

set_model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end function set_model_time_step


!------------------------------------------------------------------


subroutine get_clm_restart_filename( filename )

character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(clm_restart_filename)

end subroutine get_clm_restart_filename


!------------------------------------------------------------------
!>  This routine checks the user input against the variables available in the
!>  input netcdf file to see if it is possible to construct the DART state vector
!>  specified by the input.nml:model_nml:clm_variables  variable.
!>  Each variable must have 6 entries.
!>  1: variable name
!>  2: DART KIND
!>  3: minimum value - as a character string - if none, use 'NA'
!>  4: maximum value - as a character string - if none, use 'NA'
!>  5: what file contains the variable - '.r. => restart', '.h0. => h0 history file'
!>  6: does the variable get updated in the restart file or not ...
!>     only variables from restart files may be updated.
!>     'UPDATE'       => update the variable in the restart file
!>     'NO_COPY_BACK' => do not copy the variable back to the restart file
!>     all these variables will be updated INTERNALLY IN DART
!>     only variables marked '.r', 'UPDATE' will be modified for CLM.
!>
!>  The calling code should check to see if the variable exists.

subroutine parse_variable_table( state_variables, ngood, table )

character(len=*), dimension(:),   intent(in)  :: state_variables
integer,                          intent(out) :: ngood
character(len=*), dimension(:,:), intent(out) :: table

character(len=*), parameter :: routine = 'parse_variable_table'

integer :: nrows, ncols, i
character(len=NF90_MAX_NAME) :: varname       ! column 1
character(len=NF90_MAX_NAME) :: dartstr       ! column 2
character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
character(len=NF90_MAX_NAME) :: origin_file   ! column 5
character(len=NF90_MAX_NAME) :: state_or_aux  ! column 6

nrows = size(table,1)
ncols = size(table,2)

! This loop just repackages the 1D array of values into a 2D array.
! We can do some miniminal checking along the way.
! Determining which file to check is going to be more complicated.

ngood = 0
MyLoop : do i = 1, nrows

   varname      = trim(state_variables(ncols*i - 5))
   dartstr      = trim(state_variables(ncols*i - 4))
   minvalstring = trim(state_variables(ncols*i - 3))
   maxvalstring = trim(state_variables(ncols*i - 2))
   origin_file  = trim(state_variables(ncols*i - 1))
   state_or_aux = trim(state_variables(ncols*i    ))

   call to_upper(origin_file)
   call to_upper(state_or_aux)

   table(i,VT_VARNAMEINDX) = trim(varname)
   table(i,VT_KINDINDX)    = trim(dartstr)
   table(i,VT_MINVALINDX)  = trim(minvalstring)
   table(i,VT_MAXVALINDX)  = trim(maxvalstring)
   table(i,VT_ORIGININDX)  = trim(origin_file)
   table(i,VT_STATEINDX)   = trim(state_or_aux)

   ! If the first element is empty, we have found the end of the list.
   if ( table(i,1) == ' ' ) exit MyLoop

   ! Any other condition is an error.
   if ( any(table(i,:) == ' ') ) then
      string1 = 'input.nml &model_nml:clm_variables not fully specified'
      string2 = 'must be 6 entries per variable. Last known variable name is'
      string3 = '['//trim(table(i,1))//'] ... (without the [], naturally)'
      call error_handler(E_ERR,routine,string1,source,text2=string2,text3=string3)
   endif

   ! Make sure DART kind is valid

   if( get_index_for_quantity(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,routine,string1,source)
   endif

   ! Record the contents of the DART state vector

   if ((debug > 8) .and. do_output()) then
      write(logfileunit,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2)),' ', &
                                               trim(table(i,3)), ' ', trim(table(i,4)),' ', &
                                               trim(table(i,5)), ' ', trim(table(i,6))
      write(     *     ,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2)),' ', &
                                               trim(table(i,3)), ' ', trim(table(i,4)),' ', &
                                               trim(table(i,5)), ' ', trim(table(i,6))
   endif

   ngood = ngood + 1
enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_ALLMSG,routine,string1,source,text2=string2)
endif

! Check to see if zsno is part of the requested variables

end subroutine parse_variable_table


!------------------------------------------------------------------
!> icol            ... which CLM 'column' are we in
!> dimname         ... is it dimensioned 'levgrnd' or 'levsno' or 'levtot' ...
!> enlevels        ... the expected number of levels ... varshape
!> levtot          ... the output array of vertical coordinates
!>
!> The total number of levels is defined to be the soil levels (fixed)
!> plus the number of snow levels, which can vary by column.
!> The history file contains the depths of the soil levels (ncdf var 'levgrnd');
!> these are in levtot(1:nlevgrnd).
!> The restart file contains the depths of the snow levels (ncdf var 'ZSNO').
!> The tricky bit is that they are in reverse order ... and differ from one model to another.
!> If you simply grab the netcdf variable zsno,
!> the level closest to the soil is actually the highest index.
!>
!> From Matlab (which indexes like Fortran)
!>> size(zsno) ans = 13693           5
!>> zsno(1,:)  ans = -1.4202   -1.3852   -1.3052   -1.1352   -0.5101
!>                      |          |         |         |        |...... closest to soil surface
!>                      |          |         |         |............... one level 'up'
!>                      |          |         |......................... one level 'up'
!>                      |          |................................... one level 'up'
!>                      |.............................................. closest to sun
!>
!> If there is no snow ... the corresponding zsno is ZERO ...
!>> zsno(508,:) ans = 0   -0.5736   -0.5389   -0.4591   -0.2021
!>
!> The following Matlab code may be used to explore a variable's storage order.
!> (a better version is in the clm/matlab/CheckStorageOrder.m function)
!>
!> h2o = nc_varget(fname,'H2OSOI_LIQ');
!> h2o(h2o > 1.0E30) = NaN;
!> lat = nc_varget(fname,'cols1d_lat');
!> lon = nc_varget(fname,'cols1d_lon');
!> figure(1); plot3(lon,lat,h2o(:,1),'x'); hold on; worldmap; view(0,90)
!> figure(2); plot3(lon,lat,h2o(:,2),'x'); hold on; worldmap; view(0,90)
!> figure(3); plot3(lon,lat,h2o(:,3),'x'); hold on; worldmap; view(0,90)
!> figure(4); plot3(lon,lat,h2o(:,4),'x'); hold on; worldmap; view(0,90)
!> figure(5); plot3(lon,lat,h2o(:,5),'x'); hold on; worldmap; view(0,90)
!> figure(6); plot3(lon,lat,h2o(:,6),'x'); hold on; worldmap; view(0,90)

subroutine fill_levels(varname,icol,dimname,enlevels,levtot)

character(len=*), intent(in)    :: varname
integer,          intent(in)    :: icol
character(len=*), intent(in)    :: dimname
integer,          intent(in)    :: enlevels
real(r8),         intent(out)   :: levtot(:)

integer :: j

if     (dimname == 'levsno') then

   if (nlevsno /= enlevels) then
      write(string1,*) trim(varname),' dimension ', trim(dimname), &
                       ' has declared length ',enlevels
      write(string2,*) 'not the known number of snow levels ',nlevsno
      call error_handler(E_ERR,'fill_levels',string1,source,text2=string2)
   endif
   levtot(1:nlevsno) = zsno(1:nlevsno,icol)

elseif (dimname == 'levsno1') then

   if (nlevsno1 /= enlevels) then
      write(string1,*) trim(varname),' dimension ', trim(dimname), &
                       ' has declared length ',enlevels
      write(string2,*) 'not the known number of snow interfaces ',nlevsno1
      call error_handler(E_ERR,'fill_levels',string1,source,text2=string2)
   endif
   levtot(1:nlevsno1) = zisno(1:nlevsno1,icol)

elseif (dimname == 'levgrnd') then

   if (nlevgrnd /= enlevels) then
      write(string1,*) trim(varname),' dimension ', trim(dimname), &
                       ' has declared length ',enlevels
      write(string2,*) 'not the known number of ground levels ',nlevgrnd
      call error_handler(E_ERR,'fill_levels',string1,source,text2=string2)
   endif
   levtot(1:nlevgrnd) = LEVGRND

elseif (dimname == 'levsoi') then

   if (nlevsoi /= enlevels) then
      write(string1,*) trim(varname),' dimension ', trim(dimname), &
                       ' has declared length ',enlevels
      write(string2,*) 'not the known number of soil levels ',nlevsoi
      call error_handler(E_ERR,'fill_levels',string1,source,text2=string2)
   endif
   levtot(1:nlevsoi) = LEVSOI

elseif (dimname == 'levdcmp') then

   if (nlevdcmp /= enlevels) then
      write(string1,*) trim(varname),' dimension ', trim(dimname), &
                       ' has declared length ',enlevels
      write(string2,*) 'not the known number of decomposition levels ',nlevdcmp
      call error_handler(E_ERR,'fill_levels',string1,source,text2=string2)
   endif
   levtot(1:nlevdcmp) = LEVDCMP

elseif (dimname == 'levtot') then

   ! This block assumes anything dimensioned 'levtot' has the first nlevsno levels
   ! followed by nlevgrnd levels. Dunno what to do with lake stuff ...

   if (nlevtot /= enlevels) then
      write(string1,*) trim(varname),' dimension ', trim(dimname), &
                       ' has declared length ',enlevels
      write(string2,*) 'not the known number of total levels ',nlevtot
      call error_handler(E_ERR,'fill_levels',string1,source,text2=string2)
   endif

   if (nlevtot /= nlevgrnd + nlevsno) then
      write(string1,*) trim(varname),' nlevtot ', nlevtot, &
                       ' is not equal to nlevgrnd + nlevsno'
      write(string2,*) 'nlevgrnd is ',nlevgrnd,' nlevsno is ',nlevsno, &
                       ' total of ',nlevgrnd+nlevsno
      call error_handler(E_ERR,'fill_levels',string1,source,text2=string2)
   endif

   levtot(1:nlevsno) = zsno(1:nlevsno,icol)
   levtot(nlevsno+1:nlevsno+nlevgrnd) = LEVGRND

elseif (dimname == 'numrad') then

   if (nnumrad /= enlevels) then
      write(string1,*) trim(varname),' dimension ', trim(dimname), &
                       ' has declared length ',enlevels
      write(string2,*) 'not the known number of radiation levels ',nnumrad
      call error_handler(E_ERR,'fill_levels',string1,source,text2=string2)
   endif
   levtot(1:nnumrad) = (/ (j,j=1,nnumrad) /)

else
   write(string1,*) 'Unable to determine vertical coordinate for column ',icol
   write(string2,*) 'Variable in question is: "',trim(varname),'"'
   write(string3,*) 'unknown dimension name: "',trim(dimname),'"'
   call error_handler(E_ERR,'fill_levels',string1,source,text2=string2,text3=string3)
endif


end subroutine fill_levels


!------------------------------------------------------------------


subroutine gridcell_components(varstring)

! In order to exercise some of the routines, it is necessary to know
! which gridcells have multiple land units
! which gridcells have multiple columns
! which gridcells have multiple PFTs
!
! This routine simply tells me which gridcells are 'interesting'.
! Each level counts separately. 1 column with 20 levels ... yields a count of 20
!
! It is very similar to SetLocatorArrays(), but only does the landunit,column,pft
! that the variable uses. It is also public - currently only used by
! model_mod_check.

character(len=*), intent(in) :: varstring    ! T_SOISNO, H2OSOI_LIQ

! Local storage

integer :: idom, ivar, indexi, i, j
integer, allocatable, dimension(:,:) :: countmat

if ( .not. module_initialized ) call static_init_model

allocate(countmat(nlon,nlat))
countmat = 0

do idom = 1, get_num_domains()
   VARTYPES : do ivar = 1,get_num_variables(idom)
       ! Skip to the right variable
       if ( get_variable_name(idom,ivar) /= varstring) cycle VARTYPES

       ! Create a count of all the multiples in a gridcell
       do indexi = get_index_start(idom, ivar), &
                   get_index_end(  idom, ivar)
          i = lonixy(indexi)
          j = latjxy(indexi)
          countmat(i,j) = countmat(i,j) + 1
       enddo
   enddo VARTYPES
enddo 

if (do_output() .and. (debug > 0)) then
   write(*,*)'exploring '//trim(varstring)

   do j = 1,nlat
   do i = 1,nlon
      if ( countmat(i,j) > 1) &
         write(*,'(''gridcell'',2(1x,i8),'' has '',i6,'' lon/lat'',2(1x,f12.7))') &
                             i,j,countmat(i,j),LON(i),LAT(j)
   enddo
   enddo
endif

deallocate(countmat)

end subroutine gridcell_components


! notes on using scale_factor and add_offset 
! If _FillValue is defined then it should be scalar and of the same type as the variable.
! If the variable is packed using scale_factor and add_offset attributes (see below),
! the _FillValue attribute should have the data type of the packed data.
!
! missing_value
! When scale_factor and add_offset are used for packing, the value(s) of the missing_value
! attribute should be specified in the domain of the data in the file (the packed data),
! so that missing values can be detected before the scale_factor and add_offset are applied.
!
! scale_factor
! If present for a variable, the data are to be multiplied by this factor after the data
! are read by the application that accesses the data.  If valid values are specified using
! the valid_min, valid_max, valid_range, or _FillValue attributes, those values should be
! specified in the domain of the data in the file (the packed data), so that they can be
! interpreted before the scale_factor and add_offset are applied.
!
! add_offset
! If present for a variable, this number is to be added to the data after it is read by
! the application that accesses the data. If both scale_factor and add_offset attributes
! are present, the data are first scaled before the offset is added.



function get_model_time()
type(time_type) :: get_model_time

if ( .not. module_initialized ) call static_init_model

get_model_time = model_time

end function get_model_time


!-------------------------------------------------------------------------------
!> Find the first occurrence of the desired DART QTY of interest. 
!> The first domain, the first variable. This variable will be used for all
!> forward operator calculations.

subroutine FindVarOfInterest(qty_index, caller, dom_id, var_id, varname)

integer,          intent(in)  :: qty_index
character(len=*), intent(in)  :: caller
integer,          intent(out) :: dom_id
integer,          intent(out) :: var_id
character(len=*), intent(out) :: varname

integer :: idom
character(len=obstypelength) :: kind_string

logical, save :: warned = .false.

dom_id = -1
var_id = -1
varname = "missing_varname"

DOMAINS : do idom = 1,get_num_domains()

   var_id = get_varid_from_kind(idom, qty_index) 

   if (var_id > 0) then
      dom_id = idom
      varname = get_variable_name(dom_id,var_id)
      return
   endif

enddo DOMAINS

if (.not. warned .and. debug > 0) then
   kind_string = get_name_for_quantity( qty_index )
   write(string1,*) trim(caller)//' cannot find "'//trim(kind_string)//'" in list of DART state variables.'
   write(string2,*) trim(caller)//' looking for DART KIND (index) ', qty_index
   call error_handler(E_WARN,'FindVarOfInterest',string1,source,text2=string2)
   warned = .true.
endif

end subroutine FindVarOfInterest


!-------------------------------------------------------------------------------
!> This function will create the relational table that will indicate how many
!> and which columns pertain to the gridcells. A companion function will
!> return the column indices that are needed to recreate the gridcell value.
!>
!> This fills the gridCellInfo(:,:) structure.
!> given a gridcell, the gridCellInfo(:,:) structure will indicate how many and
!> which columns are part of the gridcell.

subroutine SetLocatorArrays()

integer :: ilon, ilat, ij, iunit
integer :: icol, currenticol(nlon,nlat)
integer :: ipft, currentipft(nlon,nlat)

gridCellInfo(:,:)%ncols = 0
gridCellInfo(:,:)%npfts = 0

! Count up how many columns are in each gridcell

do ij = 1,ncolumn
   ilon = cols1d_ixy(ij)
   ilat = cols1d_jxy(ij)
   gridCellInfo(ilon,ilat)%ncols = gridCellInfo(ilon,ilat)%ncols + 1
enddo

! Count up how many pfts are in each gridcell

do ij = 1,npft
   ilon = pfts1d_ixy(ij)
   ilat = pfts1d_jxy(ij)
   gridCellInfo(ilon,ilat)%npfts = gridCellInfo(ilon,ilat)%npfts + 1
enddo

! Create storage for the list of column,pft indices

do ilat = 1,nlat
do ilon = 1,nlon
   if ( gridCellInfo(ilon,ilat)%ncols > 0 ) then
      allocate( gridCellInfo(ilon,ilat)%columnids( gridCellInfo(ilon,ilat)%ncols ))
   endif
   if ( gridCellInfo(ilon,ilat)%npfts > 0 ) then
      allocate( gridCellInfo(ilon,ilat)%pftids( gridCellInfo(ilon,ilat)%npfts ))
   endif
enddo
enddo

! Fill the column pointer arrays

currenticol(:,:) = 0
do ij = 1,ncolumn

   ilon = cols1d_ixy(ij)
   ilat = cols1d_jxy(ij)

   currenticol(ilon,ilat) = currenticol(ilon,ilat) + 1
   icol = currenticol(ilon,ilat)

   if ( icol <= gridCellInfo(ilon,ilat)%ncols ) then
      gridCellInfo(ilon,ilat)%columnids(icol) = ij
   else
      write(string1,'(''gridcell('',i4,'','',i4,'') has at most '',i4,'' columns.'')') &
         ilon, ilat, gridCellInfo(ilon,ilat)%ncols
      write(string2,'(''Found '',i8,'' at dart index '',i12)') icol, ij
      call error_handler(E_ERR,'SetLocatorArrays',string1,source,text2=string2)
   endif
enddo

! Fill the pft pointer arrays

currentipft(:,:) = 0
do ij = 1,npft

   ilon = pfts1d_ixy(ij)
   ilat = pfts1d_jxy(ij)

   currentipft(ilon,ilat) = currentipft(ilon,ilat) + 1
   ipft = currentipft(ilon,ilat)

   if ( ipft <= gridCellInfo(ilon,ilat)%npfts ) then
      gridCellInfo(ilon,ilat)%pftids(ipft) = ij
   else
      write(string1,'(''gridcell('',i4,'','',i4,'') has at most '',i4,'' pfts.'')') &
         ilon, ilat, gridCellInfo(ilon,ilat)%npfts
      write(string2,'(''Found '',i8,'' at dart index '',i12)') ipft, ij
      call error_handler(E_ERR,'SetLocatorArrays',string1,source,text2=string2)
   endif
enddo

! Check block

if ((debug > 0) .and. do_output()) then

   iunit = open_file('gridcell_column_table.txt',form='formatted',action='write')

   do ilon = 1,nlon
   do ilat = 1,nlat
      if (gridCellInfo(ilon,ilat)%ncols > 0) then
         write(iunit,'(''gridcell'',i8,1x,i8,'' has '', i6, '' columns:'')') &
                   ilon,ilat,gridCellInfo(ilon,ilat)%ncols
         write(iunit,*)gridCellInfo(ilon,ilat)%columnids
      endif
   enddo
   enddo

   call close_file(iunit)

   iunit = open_file('gridcell_pft_table.txt',form='formatted',action='write')

   do ilon = 1,nlon
   do ilat = 1,nlat
      if (gridCellInfo(ilon,ilat)%npfts > 0) then
         write(iunit,'(''gridcell'',i8,1x,i8,'' has '', i6, '' pfts : '')') &
                   ilon,ilat,gridCellInfo(ilon,ilat)%npfts
         write(iunit,*)gridCellInfo(ilon,ilat)%pftids
      endif
   enddo
   enddo

   call close_file(iunit)

endif

end subroutine SetLocatorArrays


!-------------------------------------------------------------------------------
!> FindDesiredTimeIndx() returns the index into the time array that matches
!> the model_time from the CLM restart file.

function FindDesiredTimeIndx(ncid, ntimes, varname)
integer,          intent(in) :: ncid
integer,          intent(in) :: ntimes
character(len=*), intent(in) :: varname
integer                      :: FindDesiredTimeIndx

integer :: VarID
real(r8),  dimension(ntimes) :: mytimes
character(len=NF90_MAX_NAME) :: attvalue

type(time_type) :: thistime
integer :: ios, itime, basedays, baseseconds
integer :: iyear, imonth, iday, ihour, imin, isec

FindDesiredTimeIndx = MISSING_I   ! initialize to failure setting

call nc_check(nf90_inq_varid(ncid, 'time', VarID), &
        'FindDesiredTimeIndx:', 'inq_varid time '//varname)
call nc_check(nf90_get_var(  ncid, VarID, mytimes), &
        'FindDesiredTimeIndx:', 'get_var   time '//varname)
call nc_check(nf90_get_att(ncid, VarID, 'units', attvalue), &
        'FindDesiredTimeIndx:', 'time get_att units '//varname)

! time:units = "days since 2004-01-01 00:00:00" ;
!               1234567890

if (attvalue(1:10) /= 'days since') then
   write(string1,*)'expecting time units of [days since ... ]'
   write(string2,*)'read time units of ['//trim(attvalue)//']'
   call error_handler(E_ERR,'FindDesiredTimeIndx:',string1,source,text2=string2)
endif

read(attvalue,'(11x,i4,5(1x,i2))',iostat=ios)iyear,imonth,iday,ihour,imin,isec
if (ios /= 0) then
   write(string1,*)'Unable to read time units. Error status was ',ios
   write(string2,*)'expected "days since YYYY-MM-DD HH:MM:SS"'
   write(string3,*)'was      "'//trim(attvalue)//'"'
   call error_handler(E_ERR,'FindDesiredTimeIndx:',string1, &
          source,text2=string2,text3=string3)
endif

thistime = set_date(iyear, imonth, iday, ihour, imin, isec)
call get_time(thistime, baseseconds, basedays)

! convert each time to a DART time and compare to desired

TIMELOOP : do itime = 1,ntimes

   iday     = int(mytimes(itime))
   isec     = (mytimes(itime) - iday)*86400
   thistime = set_time(baseseconds+isec, basedays+iday)

   if (thistime == model_time) then
      FindDesiredTimeIndx = itime
      exit TIMELOOP
   endif

enddo TIMELOOP

! FIXME ... do we actually need a perfect match ... or do we just use the last
! one. History files may not have a perfect match to the restart file.
if ( FindDesiredTimeIndx == MISSING_I ) then
   call print_time(model_time,str='model time is ',iunit=logfileunit)
   call print_time(model_time,str='model time is ')
   call print_date(model_time,str='model date is ',iunit=logfileunit)
   call print_date(model_time,str='model date is ')
   write(string1,*)'No time matching model_time found for '//trim(varname)
   call error_handler(E_ERR,'FindDesiredTimeIndx:',string1,source)
endif

if (debug > 0) then
   write(string1,*)trim(varname)//' matching time index is ',FindDesiredTimeIndx
   call error_handler(E_ALLMSG,'FindDesiredTimeIndx:',string1,source)
endif

end function FindDesiredTimeIndx


!-----------------------------------------------------------------------
!> Collect all the variables that are from the specific file.
!> Some variables come from the restart file, some from the history file, ...
!> Each file specifies a new 'domain'.

subroutine cluster_variables(table, origin, nvars, var_names, var_qtys, var_ranges, var_update)

character(len=*), intent(in)  :: table(:,:)
character(len=*), intent(in)  :: origin
integer,          intent(out) :: nvars
character(len=*), intent(out) :: var_names(:)
integer,          intent(out) :: var_qtys(:)
real(r8),         intent(out) :: var_ranges(:,:)
logical,          intent(out) :: var_update(:)

character(len=*), parameter :: routine = 'cluster_variables'

integer :: ivar, ios
real(r8) :: minvalue, maxvalue

nvars      = 0
var_names  = 'no_variable_specified'
var_qtys   = MISSING_I
var_ranges = MISSING_R8
var_update = .false.

domain_count = domain_count + 1

VARLOOP : do ivar = 1,size(table,1)

   ! Skip if the variable is for a different domain
   if (table(ivar,VT_ORIGININDX) /= origin) cycle VARLOOP

   nvars = nvars + 1
   var_names( nvars)  = trim(table(ivar,VT_VARNAMEINDX))
   var_qtys(  nvars)  = get_index_for_quantity(table(ivar,VT_KINDINDX))

   ! If the character string can be interpreted as an r8, great.
   ! If not, there is no value to be used.

   var_ranges(nvars,1) = MISSING_R8
   var_ranges(nvars,2) = MISSING_R8

   read(table(ivar,VT_MINVALINDX),*,iostat=ios) minvalue
   if (ios == 0) var_ranges(nvars,1) = minvalue

   read(table(ivar,VT_MAXVALINDX),*,iostat=ios) maxvalue
   if (ios == 0) var_ranges(nvars,2) = maxvalue

   ! You may want to inflate the history file variables used for forward operators.
   ! You can always write to other output files and DART will create new ones
   ! without actually changing the values in the input history files.
   ! [output_state_file_list /= input_state_file_list]

   if ( table(ivar,VT_STATEINDX) == 'UPDATE' ) var_update(nvars) = .true.

   if (debug > 0 .and. do_output()) then
      call say(' ')
      write(string1,*)trim(origin),' defines domain ',domain_count, &
                      ' var "',trim(var_names(nvars)),'"'
      write(string2,*)'variable has declared range ', &
                      var_ranges(nvars,1),var_ranges(nvars,2)
      write(string3,*)'and quantity index ',var_qtys(nvars), &
                      '. Is var on update list:',var_update(nvars)
      call error_handler(E_MSG,routine,string1,source,text2=string2,text3=string3)
   endif

enddo VARLOOP

end subroutine cluster_variables


!------------------------------------------------------------------

subroutine say(what)
character(len=*), intent(in) :: what

write(logfileunit, *) trim(what)
write(  *        , *) trim(what)

end subroutine say


!------------------------------------------------------------------
!> Declare all urban areas and lakes 'unrelated'
!> to whatever observations we are considering. It is actually
!> easier to specify that we only want to update vegetated or crop columns/landunits.
!> At some point, we may want to update land ice or glaciers or ... 
!> but that will probably entail an expanded set of namelist options.

function related(obs_quantity, state_index)

integer,     intent(in) :: obs_quantity
integer(i8), intent(in) :: state_index
logical                 :: related

! land1d_ityplun:long_name = "landunit type (vegetated,urban,lake,wetland or glacier)" ;
! cols1d_ityplun:long_name = "column landunit type (vegetated,urban,lake,wetland or glacier)" ;
! pfts1d_ityplun:long_name = "pft landunit type (vegetated,urban,lake,wetland or glacier)" ;

integer :: jdim
integer :: var_id, dom_id, var_type
character(len=32) :: var_kind_string
character(len=NF90_MAX_NAME) :: dimension_name
integer :: indices(3)

related = .false.

! from the dart index get a whole host of information
call get_model_variable_indices(state_index, indices(1), indices(2), indices(3), &
            var_id=var_id, dom_id=dom_id, kind_index=var_type, kind_string=var_kind_string)

! Determine if state_index is a variable from a column (or whatever is of interest).
! Determine what dimension is of interest, need to know to index into
! cols1d_ityplun(ncolumn) array (for example).

RELATEDLOOP: do jdim = 1, get_num_dims(dom_id, var_id)

   dimension_name = get_dim_name(dom_id, var_id, jdim)
   select case ( trim(dimension_name) )
          case ("gridcell","lon","lat")
             related = .true.
          case ("lndgrid")
             related = .true.
          case ("landunit")

             if ( land1d_ityplun(indices(jdim)) == ilun_vegetated_or_bare_soil ) related = .true.
             if ( land1d_ityplun(indices(jdim)) == ilun_crop                   ) related = .true.

          case ("column")
                       
             if ( cols1d_ityplun(indices(jdim)) == icol_vegetated_or_bare_soil ) related = .true.
             if ( cols1d_ityplun(indices(jdim)) == icol_crop                   ) related = .true.

          case ("pft")
             related = .true.
             !write(*,*)'related: TJH pfts1d_ityplun ', &
             !          pfts1d_ityplun(indices(jdim)),trim(dimension_name)
          case default
   end select

   ! Since variables can use only one of these dimensions,
   ! there is no need to check the other dimensions. 
   if (related) exit RELATEDLOOP

enddo RELATEDLOOP

! write(*,*)'related: ',obs_quantity, state_index, var_type, trim(var_kind_string), related

end function related


!-------------------------------------------------------------------------------
!> horrible way to fill the metadata arrays:
!> lonixy   contains the index to the longitude array
!> latjxy   contains the index to the latitude array
!> landarea contains the weight to use when calculating area-weighted quantities

subroutine fill_rank1_metadata(varname, dimension_name, dimension_length, indx)

character(len=*), intent(in)    :: varname
character(len=*), intent(in)    :: dimension_name
integer,          intent(in)    :: dimension_length
integer(i8),      intent(inout) :: indx

integer :: i, xi, xj

if ((debug > 0) .and. do_output()) then
   write(*,*)
   write(*,*)'variable "'//trim(varname)//'"'
   write(*,*)'dimension 1 (i) ',trim(dimension_name),dimension_length
endif

SELECT CASE ( trim(dimension_name) )
   CASE ("gridcell","lndgrid")
      do i = 1, dimension_length
         xi             = grid1d_ixy(i)
         xj             = grid1d_jxy(i) ! nnnnn_jxy(:) always 1 if unstructured
         if (unstructured) then
            lonixy(  indx) = xi
            latjxy(  indx) = xi
            landarea(indx) = AREA1D(xi) * LANDFRAC1D(xi)
         else
            lonixy(  indx) = xi
            latjxy(  indx) = xj
            landarea(indx) = AREA2D(xi,xj) * LANDFRAC2D(xi,xj)
         endif
         indx = indx + 1
      enddo

   CASE ("landunit")
      do i = 1, dimension_length
         xi             = land1d_ixy(i)
         xj             = land1d_jxy(i) ! nnnnn_jxy(:) always 1 if unstructured
         if (unstructured) then
            lonixy(  indx) = xi
            latjxy(  indx) = xi
            landarea(indx) = AREA1D(xi) * LANDFRAC1D(xi) * land1d_wtxy(i)
         else
            lonixy(  indx) = xi
            latjxy(  indx) = xj
            landarea(indx) = AREA2D(xi,xj) * LANDFRAC2D(xi,xj) * land1d_wtxy(i)
         endif
         indx = indx + 1
      enddo

   CASE ("column")
      do i = 1, dimension_length
         xi             = cols1d_ixy(i)
         xj             = cols1d_jxy(i) ! nnnnn_jxy(:) always 1 if unstructured
         if (unstructured) then
            lonixy(  indx) = xi
            latjxy(  indx) = xi
            landarea(indx) = AREA1D(xi) * LANDFRAC1D(xi) * cols1d_wtxy(i)
         else
            lonixy(  indx) = xi
            latjxy(  indx) = xj
            landarea(indx) = AREA2D(xi,xj) * LANDFRAC2D(xi,xj) * cols1d_wtxy(i)
         endif
         indx = indx + 1
      enddo

   CASE ("pft")
      do i = 1, dimension_length
         xi             = pfts1d_ixy(i)
         xj             = pfts1d_jxy(i) ! nnnnn_jxy(:) always 1 if unstructured
         if (unstructured) then
            lonixy(  indx) = xi
            latjxy(  indx) = xi
            landarea(indx) = AREA1D(xi) * LANDFRAC1D(xi) * pfts1d_wtxy(i)
         else
            lonixy(  indx) = xi
            latjxy(  indx) = xj
            landarea(indx) = AREA2D(xi,xj) * LANDFRAC2D(xi,xj) * pfts1d_wtxy(i)
         endif
         indx = indx + 1
      enddo

   CASE DEFAULT
      write(string1,*)'unknown dimension name "'//trim(dimension_name)//'"'
      write(string2,*)' while trying to create metadata for "'//trim(varname)//'"'
      call error_handler(E_ERR,'fill_rank1_metadata',string1,source,text2=string2)

END SELECT

end subroutine fill_rank1_metadata


!-------------------------------------------------------------------------------
!> Horrible way to fill the metadata arrays: The arrays have to be filled in the
!> same order as the data is stored in the DART state vector, which happens
!> somewhere else in DART.
!> lonixy   contains the index to the longitude array
!> latjxy   contains the index to the latitude array
!> landarea contains the weight to use when calculating area-weighted quantities

subroutine fill_rank2_metadata(varname, dimension_name, dimension_length, indx)

character(len=*), intent(in)    :: varname
character(len=*), intent(in)    :: dimension_name(2)
integer,          intent(in)    :: dimension_length(2)
integer(i8),      intent(inout) :: indx

! RANK 2 VARIABLES (as summarized by 'ncdump') - 'time' dimension omitted
!         RESTART FILES          HISTORY_FILES          VECTOR_FILES
!      X(  column, levtot)       X(lat, lon)            X(lat,lon)
!      X(  column, levgrnd)                             X(levgrnd, column)
!      X(  column, levsno)                              X(levsoi,  column)
!      X(  column, levlak)                              X(levdcmp, column)
!      X(  column, numrad )
!      X(  column, levsno1)
!      X(     pft, levgrnd)
!      X(     pft, levcan)
!      X(     pft, numrad)
!      X(landunit, numrad)

integer :: i, j, xi, xj
integer :: horiz_dimension, other_dimension

! sometimes the horizontal dimension is 1, sometimes it is 2
! In this context, gridcell,lat,colum,pft,landunit are 'horizontal',
! but 'lon' is not ... the 'horizontal dimension' is used to index the
! first dimension of the LANDFRAC2D variable ... which is 'lat' 

call parse_dimensions(varname,dimension_name,horiz_dimension,other_dimension)

! If the variable is dimensioned (*vertical*,column) the order of the looping must
! be reversed to match the storage order in the DART state vector.

if ((debug > 0) .and. do_output()) then
   write(*,*)
   write(*,*)'variable ',trim(varname)
   write(*,*)'dimension 1 (i) ',trim(dimension_name(1)),dimension_length(1)
   write(*,*)'dimension 2 (j) ',trim(dimension_name(2)),dimension_length(2)
   write(*,*)'horizontal dimension : ',trim(dimension_name(horiz_dimension))
   write(*,*)'other      dimension : ',trim(dimension_name(other_dimension))
endif

SELECT CASE ( trim(dimension_name(horiz_dimension)) )
   CASE ("gridcell")
      if ((debug > 8) .and. do_output()) write(*,*)'length grid1d_ixy ',size(grid1d_ixy)
      do j = 1, dimension_length(horiz_dimension)
         xi = grid1d_ixy(j)
         xj = grid1d_jxy(j) ! nnnnn_jxy(:) always 1 if unstructured
         do i = 1, dimension_length(other_dimension)
            if (unstructured) then
               lonixy(  indx) = xi
               latjxy(  indx) = xi
               landarea(indx) = AREA1D(xi) * LANDFRAC1D(xi)
            else
               lonixy(  indx) = xi
               latjxy(  indx) = xj
               landarea(indx) = AREA2D(xi,xj) * LANDFRAC2D(xi,xj)
            endif
            indx = indx + 1
         enddo
      enddo

   CASE ("landunit")
      if ((debug > 8) .and. do_output()) write(*,*)'length land1d_ixy ',size(land1d_ixy)
      do j = 1, dimension_length(horiz_dimension)
         xi = land1d_ixy(j)
         xj = land1d_jxy(j) ! nnnnn_jxy(:) always 1 if unstructured
         do i = 1, dimension_length(other_dimension)
            if (unstructured) then
               lonixy(  indx) = xi
               latjxy(  indx) = xi
               landarea(indx) = AREA1D(xi) * LANDFRAC1D(xi) * land1d_wtxy(j)
            else
               lonixy(  indx) = xi
               latjxy(  indx) = xj
               landarea(indx) = AREA2D(xi,xj) * LANDFRAC2D(xi,xj) * land1d_wtxy(j)
            endif
            indx = indx + 1
         enddo
      enddo

   CASE ("column") ! Column is the only coordinate that has vertical levels.
                   ! Sometimes the storage order is (column,level), sometimes
                   ! it is (level,column). Remember netCDF reporting is opposite
                   ! to Fortran.

      if (horiz_dimension == 2) then
         call loop_other_fastest(varname, indx, dimension_name, dimension_length, &
                                      horiz_dimension, other_dimension)
      else
         call loop_column_fastest(varname, indx, dimension_name, dimension_length, &
                                      horiz_dimension, other_dimension)
      endif


   CASE ("pft")
      if ((debug > 8) .and. do_output()) write(*,*)'length pfts1d_ixy ',size(pfts1d_ixy)
      do j = 1, dimension_length(horiz_dimension)
         xi = pfts1d_ixy(j)
         xj = pfts1d_jxy(j) ! nnnnn_jxy(:) always 1 if unstructured
         do i = 1, dimension_length(other_dimension)
            if (unstructured) then
               lonixy(  indx) = xi
               latjxy(  indx) = xi
               landarea(indx) = AREA1D(xi) * LANDFRAC1D(xi) * pfts1d_wtxy(j)
            else
               lonixy(  indx) = xi
               latjxy(  indx) = xj
               landarea(indx) = AREA2D(xi,xj) * LANDFRAC2D(xi,xj) * pfts1d_wtxy(j)
            endif
            indx = indx + 1
         enddo
      enddo

   CASE ("lat")
      do j = 1, dimension_length(horiz_dimension)
         do i = 1, dimension_length(other_dimension)
            if (unstructured) then
               lonixy(  indx) = i
               latjxy(  indx) = i
               landarea(indx) = AREA1D(i) * LANDFRAC1D(i)
               call error_handler(E_ERR,'fill_rank2_metadata','section untested',source)
            else
               lonixy(  indx) = i
               latjxy(  indx) = j
               landarea(indx) = AREA2D(i,j) * LANDFRAC2D(i,j)
            endif
            indx = indx + 1
         enddo
      enddo

   CASE DEFAULT


      write(string1,*)'unsupported dimension "'//trim(dimension_name(horiz_dimension))//'"'
      write(string2,*)' while trying to create metadata for "'//trim(varname)//'"'
      call error_handler(E_ERR,'fill_rank2_metadata',string1,source,text2=string2)

END SELECT

end subroutine fill_rank2_metadata


!-------------------------------------------------------------------------------
!> horrible way to fill the metadata arrays:
!> lonixy   contains the index to the longitude array
!> latjxy   contains the index to the latitude array
!> landarea contains the weight to use when calculating area-weighted quantities

subroutine fill_rank3_metadata(varname, dimension_name, dimension_length, indx)

character(len=*), intent(in)    :: varname
character(len=*), intent(in)    :: dimension_name(3)
integer,          intent(in)    :: dimension_length(3)
integer(i8),      intent(inout) :: indx

! restart file   variables  never have 3 dimensions
! vector_history variables    may have 3 dimensions [time, lat, lon]
! history file   variables always have 3 dimensions [time, lat, lon]
!     float      H2OSOI(time, levgrnd, lat, lon) ;
!     float    TSOI_ICE(time, levgrnd, lat, lon) ;
!     float PCT_GLC_MEC(time, glc_nec, lat, lon) ;
!     float  PCT_LANDUNIT(time, ltype, lat, lon) ;
!     float  PCT_NAT_PFT(time, natpft, lat, lon) ;
!     float        TLAKE(time, levlak, lat, lon) ;
!     float       VEGWP(time, nvegwcs, lat, lon) ;

integer :: i, j, k

if ((debug > 8) .and. do_output()) then
   write(*,*)
   write(*,*)'variable ',trim(varname)
   write(*,*)'dimension 1 (i) ',dimension_name(1),dimension_length(1)
   write(*,*)'dimension 2 (j) ',dimension_name(2),dimension_length(2)
   write(*,*)'dimension 3 (k) ',dimension_name(3),dimension_length(3)
endif

! Remember the order is reversed from ncdump to fortran
! The messages are in ncdump-order, checking is fortran-order
if (trim(dimension_name(1)) .ne. 'lon' .or.  trim(dimension_name(2)) .ne. 'lat') then
   write(string1,*)'3D variables must be [~,lat,lon] (as reported by ncdump)'
   write(string2,*)'"'//trim(varname)//'" is ', &
                        trim(dimension_name(3)), ' ', &
                        trim(dimension_name(2)), ' ', &
                        trim(dimension_name(1))
   call error_handler(E_ERR,'fill_rank3_metadata',string1,source,text2=string2)
endif

!>@todo  extend fill_levels to support all the vertical levels

SELECT CASE ( trim(dimension_name(3)) )
   CASE ("levgrnd")
      levtot(1:dimension_length(3)) = LEVGRND;
   CASE ("levsoi")
      levtot(1:dimension_length(3)) = LEVSOI;
   CASE ("levdcmp")
      levtot(1:dimension_length(3)) = LEVDCMP;
   CASE DEFAULT
      write(string1,*)'unsupported vertical dimension name "'//trim(dimension_name(3))//'"'
      write(string2,*)' while trying to create metadata for "'//trim(varname)//'"'
      call error_handler(E_ERR,'fill_rank3_metadata',string1,source,text2=string2)
END SELECT

do k = 1, dimension_length(3)
   do j = 1, dimension_length(2)
      do i = 1, dimension_length(1)
         levels(  indx) = levtot(k)
         lonixy(  indx) = i
         latjxy(  indx) = j
         landarea(indx) = AREA2D(i,j) * LANDFRAC2D(i,j)
         indx = indx + 1
      enddo
   enddo
enddo

end subroutine fill_rank3_metadata


!-------------------------------------------------------------------------------
!>
! RANK 2 VARIABLES (as summarized by 'ncdump') - 'time' dimension omitted
!         RESTART FILES          HISTORY_FILES          VECTOR_FILES
!      X(  column, levtot)       X(lat, lon)            X(lat,lon)
!      X(  column, levgrnd)                             X(levgrnd, column)
!      X(  column, levsno)                              X(levsoi,  column)
!      X(  column, levlak)                              X(levdcmp, column)
!      X(  column, numrad )
!      X(  column, levsno1)
!      X(     pft, levgrnd)
!      X(     pft, levcan)
!      X(     pft, numrad)
!      X(landunit, numrad)

subroutine parse_dimensions(varname,dimension_names,h_dimension,o_dimension)

character(len=*), intent(in)  :: varname
character(len=*), intent(in)  :: dimension_names(:)
integer,          intent(out) :: h_dimension, o_dimension

integer :: i

h_dimension = -1
o_dimension = -1

string2 = ''

DIMLOOP : do i = 1,size(dimension_names)

   write(string2,*) trim(string2),'"'//trim(dimension_names(i))//'",'

   SELECT CASE ( trim(dimension_names(i)) )
     CASE ("lat","landunit","column","pft")
        h_dimension = i
     CASE ("lon")
        ! the fill_rank2_metadata() routine counts on finding 'lat' to
        ! correctly index into the shape of LANDFRAC2D varible. ugly.
        continue
     CASE DEFAULT
        continue
   END SELECT

   if (h_dimension > 0) exit DIMLOOP
enddo DIMLOOP

! In the current context, it is a fatal error to call this function and not
! find a horizontal dimension or the variable is a 3D variable.

if (h_dimension < 1) then
   write(string1,*)'Unable to determine horizontal dimension for "'//trim(varname)//'"'
   write(string2,*)'dimensions are '//trim(string2)
   call error_handler(E_ERR,'parse_dimensions',string1,source,text2=string2)
elseif (h_dimension == 1 ) then
   o_dimension = 2
elseif (h_dimension == 2 ) then
   o_dimension = 1
else
   write(string1,*)'Only variables with two dimensiona are supported.'
   write(string2,*)'Unable to determine horizontal dimension for "'//trim(varname)//'"'
   write(string3,*)'dimensions are '//trim(string2)
   call error_handler(E_ERR,'parse_dimensions',string1,source,text2=string2,text3=string3)
endif

end subroutine parse_dimensions


!-------------------------------------------------------------------------------
!>

subroutine loop_other_fastest(varname, indx, dimension_names, dimension_lengths, &
                                      horiz_dimension, other_dimension)
character(len=*), intent(in) :: varname
integer(i8),          intent(inout) :: indx
character(len=*), intent(in) :: dimension_names(:)
integer,          intent(in) :: dimension_lengths(:)
integer,          intent(in) :: horiz_dimension
integer,          intent(in) :: other_dimension

integer :: i,j,xi,xj

COLUMN : do j = 1, dimension_lengths(horiz_dimension)

   call fill_levels(varname, j, dimension_names(  other_dimension), & 
                                dimension_lengths(other_dimension),levtot)

   xi = cols1d_ixy(j)
   xj = cols1d_jxy(j) ! nnnnn_jxy(:) always 1 if unstructured
   OTHER :  do i = 1, dimension_lengths(other_dimension)
      levels(  indx) = levtot(i)
      if (unstructured) then
         lonixy(  indx) = xi
         latjxy(  indx) = xi
         landarea(indx) = AREA1D(xi) * LANDFRAC1D(xi) * cols1d_wtxy(j)
      else
         lonixy(  indx) = xi
         latjxy(  indx) = xj
         landarea(indx) = AREA2D(xi,xj) * LANDFRAC2D(xi,xj) * cols1d_wtxy(j)
      endif
      indx = indx + 1
   enddo OTHER
enddo COLUMN

end subroutine loop_other_fastest


!-------------------------------------------------------------------------------
!>

subroutine loop_column_fastest(varname, indx, dimension_names, dimension_lengths, &
                                      horiz_dimension, other_dimension)
character(len=*), intent(in) :: varname
integer(i8),          intent(inout) :: indx
character(len=*), intent(in) :: dimension_names(:)
integer,          intent(in) :: dimension_lengths(:)
integer,          intent(in) :: horiz_dimension
integer,          intent(in) :: other_dimension

integer :: i,j,xi,xj

OTHER : do j = 1, dimension_lengths(other_dimension)
   COLUMN :  do i = 1, dimension_lengths(horiz_dimension)

      xi = cols1d_ixy(i)
      xj = cols1d_jxy(i) ! nnnnn_jxy(:) always 1 if unstructured

      call fill_levels(varname, i, dimension_names(  other_dimension), & 
                                   dimension_lengths(other_dimension),levtot)

      levels(  indx) = levtot(j)
      if (unstructured) then
         lonixy(  indx) = xi
         latjxy(  indx) = xi
         landarea(indx) = AREA1D(xi) * LANDFRAC1D(xi) * cols1d_wtxy(i)
      else
         lonixy(  indx) = xi
         latjxy(  indx) = xj
         landarea(indx) = AREA2D(xi,xj) * LANDFRAC2D(xi,xj) * cols1d_wtxy(i)
      endif
      indx = indx + 1
   enddo COLUMN
enddo OTHER

end subroutine loop_column_fastest


!===================================================================
! End of model_mod
!===================================================================

end module model_mod

