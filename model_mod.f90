! DART software - Copyright 2004 - 2015 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! FESOM interface to the DART data assimilation system.
! code in this module is compiled with the DART executables.  It isolates
! all information about the FESOM grids, model variables, and other details.
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

use            types_mod, only : r4, r8, i8, digits12, SECPERDAY, MISSING_R8,       &
                                 rad2deg, deg2rad, PI, MISSING_I, obstypelength

use     time_manager_mod, only : time_type, set_time, set_date, get_date, &
                                 get_time, print_time, print_date,        &
                                 set_calendar_type, increment_time,       &
                                 operator(*),  operator(+), operator(-),  &
                                 operator(>),  operator(<), operator(/),  &
                                 operator(/=), operator(<=)

use         location_mod, only : location_type, get_dist, query_location,    &
                                 set_location, get_location, write_location, &
                                 get_close_type, VERTISHEIGHT,               &
                                 loc_get_close_obs => get_close_obs,         &
                                 loc_get_close_state => get_close_state,     &
                                 is_vertical, set_vertical_localization_coord, &
                                 vertical_localization_on

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, nc_check,        &
                                 nc_begin_define_mode, nc_end_define_mode,     &
                                 nc_open_file_readonly, nc_close_file,         &
                                 nc_define_dimension

use      location_io_mod, only : nc_write_location_atts, nc_write_location

use        utilities_mod, only : register_module, error_handler,                   &
                                 E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                                 do_output, to_upper, nmlfileunit,                 &
                                 find_namelist_in_file, check_namelist_read,       &
                                 open_file, file_exist, find_textfile_dims,        &
                                 file_to_text, close_file, do_nml_file,            &
                                 do_nml_term, scalar

use         obs_kind_mod, only : get_index_for_quantity,     &
                                 get_name_for_quantity,      &
                                 QTY_VERTICAL_VELOCITY,      &
                                 QTY_POTENTIAL_TEMPERATURE,  &
                                 QTY_TEMPERATURE,            &
                                 QTY_SALINITY,               &
                                 QTY_DRY_LAND,               &
                                 QTY_EDGE_NORMAL_SPEED,      &
                                 QTY_U_CURRENT_COMPONENT,    &
                                 QTY_V_CURRENT_COMPONENT,    &
                                 QTY_SEA_SURFACE_HEIGHT,     &
                                 QTY_SEA_SURFACE_PRESSURE,   &
                                 QTY_TRACER_CONCENTRATION

use mpi_utilities_mod, only: my_task_id, broadcast_minmax, task_count

use        random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use         fesom_modules, only: read_node, read_aux3, read_depth, read_namelist, &
                                 nCells => myDim_nod2D, & ! number of surface locations
                                 nVertices => myDim_nod3D, & ! wet points in grid
                                 nVertLevels => max_num_layers, & ! number of vertical levels
                                 layerdepth, & ! depth at each level (m)
                                 coord_nod2D, &
                                 coord_nod3D, &
                                 num_layers_below_nod2d, &
                                 nod2d_corresp_to_nod3D, &
                                 nod3d_below_nod2d

use        random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use  ensemble_manager_mod, only: ensemble_type, get_my_num_vars, get_my_vars

use default_model_mod,    only : nc_write_model_vars, adv_1step,          &
                                 init_time => fail_init_time,             &
                                 init_conditions => fail_init_conditions

use distributed_state_mod, only : get_state

! netcdf modules
use typesizes
use netcdf

use state_structure_mod, only : add_domain, &
                                get_domain_size, &
                                get_model_variable_indices, &
                                state_structure_info, &
                                get_index_start, &
                                get_index_end, &
                                get_num_variables, &
                                get_variable_name, &
                                get_varid_from_kind

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
public :: get_model_size,                      & !>!DoNe
          get_num_vars,                        & !>!DoNe
          get_state_meta_data,                 & !>!DoNe
          model_interpolate,                   & !>!ToDo
          shortest_time_between_assimilations, & !>!DoNe
          static_init_model,                   & !>!ToDo add_domain...
          end_model,                           & !>!Check if other variables should be deallocated
          nc_write_model_atts,                 & !>!DoNe
          pert_model_copies,                   & !>!DoNe
          get_close_obs,                       & !>!ToDo vert_convert
          get_close_state,                     & !>!ToDo
          convert_vertical_obs,                & !>!ToDo
          convert_vertical_state,              & !>!ToDo
          read_model_time,                     & !>!ToDo how to read model time
          write_model_time                       !>!ToDo

! code for these mandatory routines are in other modules
public :: init_time,           &
          init_conditions,     &
          adv_1step,           &
          nc_write_model_vars

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.

public :: get_grid_dims,                       &
          read_2d_from_nc_file,                &
          print_variable_ranges,               &
          get_model_time_step

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

! module global storage; maintains values between calls, accessible by
! any subroutine
character(len=512) :: string1, string2, string3
logical, save :: module_initialized = .false.

! Real (physical) constants as defined exactly in FESOM.
! redefined here for consistency with the model.
real(r8), parameter :: rgas = 287.0_r8
real(r8), parameter :: cp = 1003.0_r8
real(r8), parameter :: cv = 716.0_r8
real(r8), parameter :: p0 = 100000.0_r8
real(r8), parameter :: rcv = rgas/(cp-rgas)

! earth radius; needed to convert lat/lon to x,y,z cartesian coords.
! FIXME: the example ocean files has a global attr with 6371220.0
! the actual FESOM code may have hardwired values (it did in the atmosphere)
! need to check what is really going on.
real(r8), parameter :: radius = 6371229.0 ! meters

! roundoff error
real(r8), parameter :: roundoff = 1.0e-12_r8

! Structure for computing distances to cell centers, and assorted arrays
! needed for the get_close code.

! Storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq

integer :: domid ! For state_structure_mod access

type(location_type), allocatable :: cell_locations(:)
integer,             allocatable :: cell_kinds(:)
integer,             allocatable :: close_cell_inds(:)
real(r8),            allocatable :: depths(:)

type(get_close_type)  :: cc_gc

integer, parameter :: MAX_STATE_VARIABLES = 3
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 5
integer, parameter :: VARNAME_INDEX = 1
integer, parameter ::    KIND_INDEX = 2
integer, parameter ::  MINVAL_INDEX = 3
integer, parameter ::  MAXVAL_INDEX = 4
integer, parameter :: REPLACE_INDEX = 5

character(len=obstypelength) :: var_names(max_state_variables)
real(r8) :: var_ranges(max_state_variables,2)
logical  :: var_update(max_state_variables)
integer  :: var_qtys(  max_state_variables)

! variables which are in the module namelist
integer            :: vert_localization_coord = VERTISHEIGHT
integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 21600     ! 86400 | 43200 | 21600
real(r8)           :: model_perturbation_amplitude = 0.0001   ! tiny amounts
logical            :: output_state_vector = .false.  ! output prognostic variables (if .false.)
logical            :: diagnostic_metadata = .false.
integer            :: debug = 0   ! turn up for more and more debug messages
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: model_analysis_filename = 'expno.year.oce.nc'
character(len=NF90_MAX_NAME) :: variables(NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) = ' '

namelist /model_nml/             &
   model_analysis_filename,      &
   output_state_vector,          &
   vert_localization_coord,      &
   diagnostic_metadata,          &
   assimilation_period_days,     &
   assimilation_period_seconds,  &
   model_perturbation_amplitude, &
   calendar,                     &
   variables,                    &
   debug

! Everything needed to describe a variable

integer :: nfields

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
   integer :: varsize       ! prod(dimlens(1:numdims))
   integer :: index1        ! location in dart state vector of first occurrence
   integer :: indexN        ! location in dart state vector of last  occurrence
   integer :: dart_kind
   character(len=obstypelength) :: kind_string
   logical  :: clamping     ! does variable need to be range-restricted before
   real(r8) :: range(2)     ! being stuffed back into FESOM analysis file.
   logical  :: out_of_range_fail  ! is out of range fatal if range-checking?
   logical  :: replace      ! does variable need to be stuffed back info FESOM
end type progvartype

type(progvartype), dimension(MAX_STATE_VARIABLES) :: progvar

! Grid parameters - the values will be read from an FESOM analysis file.

real(r8), allocatable :: ens_mean(:)   ! needed to convert vertical distances consistently

!$! integer :: FESOM:num_layers_below_2d(:) ! list of maximum (deepest) level index for each cell

integer         :: model_size          ! the state vector length
type(time_type) :: model_timestep      ! smallest time to adv model

! useful flags in making decisions when searching for points, etc
!$! logical :: global_grid = .true.        ! true = the grid covers the sphere with no holes
!$! logical :: all_levels_exist_everywhere = .true. ! true = cells defined at all levels

! currently unused; for a regional model it is going to be necessary to know
! if the grid is continuous around longitudes (wraps in east-west) or not,
! and if it covers either of the poles.
!#! character(len= 64) :: ew_boundary_type, ns_boundary_type

! common names that call specific subroutines based on the arg types

!TJH interface write_model_time
!TJH    module procedure write_model_time_file
!TJH    module procedure write_model_time_restart
!TJH end interface


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

logical :: state_table_needed = .false.
logical :: close_structure_allocated = .false.

contains

!==================================================================
! All the public REQUIRED interfaces come first - just by convention.
!==================================================================


!------------------------------------------------------------------
!> Called to do one time initialization of the model.

subroutine static_init_model()

! Local variables - all the important ones have module scope

integer :: ncid
integer :: iunit, io
integer :: ss, dd

if ( module_initialized ) return ! only need to do this once.

! Print module information to log file and stdout.

call register_module(source, revision, revdate)

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.

module_initialized = .true.

! Read the DART namelist for this model and
! record the namelist values used for the run

call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

!---------------------------------------------------------------
! 1) get grid dimensions
! 2) allocate space for the grids
! 3) read them from the analysis file

! Get the FESOM run-time configurations from a hardcoded filename
! must be called 'namelist.config' in the current directory.

call read_namelist()

!---------------------------------------------------------------
!>@todo Ensure model_timestep is multiple of "dynamics timestep"
!> FESOM must be able to be stopped at the requested time

call set_calendar_type( calendar )

model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)

!---------------------------------------------------------------
! Compile the list of model variables to use in the creation
! of the DART state vector. Required to determine model_size.
!
! parse and verify all variables are in the model analysis file
!
! Compute the offsets into the state vector for the start of each
! different variable type. Requires reading shapes from the model
! analysis file. As long as TIME is the LAST dimension, we're OK.
!
! Record the extent of the data type in the state vector.

ncid = nc_open_file_readonly(model_analysis_filename,'static_init_model')

call parse_variable_input(variables, ncid, model_analysis_filename, nfields)

call nc_close_file(ncid, 'static_init_model', model_analysis_filename)

!>@ TODO try to move read_grid to some other routine that gets called very
!   early by all the tasks.

call read_grid() ! sets nCells, nVertices, nVertLevels

! dump_tables() is VERY expensive and should only be called
! once to generate sanity checks.
if (state_table_needed .and. debug > 99) call dump_tables()

domid = add_domain( model_analysis_filename, nfields, &
            var_names  = var_names( 1:nfields),       &
            kind_list  = var_qtys(  1:nfields),       &
            clamp_vals = var_ranges(1:nfields,:),     &
            update_list= var_update(1:nfields)        )

if ( debug > 4 .and. do_output()) call state_structure_info(domid)

model_size = get_domain_size(domid)

allocate( ens_mean(model_size) )

! tell the location module how we want to localize in the vertical
call set_vertical_localization_coord(vert_localization_coord)

end subroutine static_init_model


!------------------------------------------------------------------
!>

subroutine get_state_meta_data(index_in, location, var_type)

! given an index into the state vector, return its location and
! if given, the var kind.   despite the name, var_type is a generic
! kind, like those in obs_kind/obs_kind_mod.f90, starting with KIND_

! passed variables

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

! Local variables

integer  :: nf, n, myindx
real(r8) :: lon, lat, depth
integer  :: dim2_j, dim3_k
integer  :: var_id, dom_id, variable_type

if ( .not. module_initialized ) call static_init_model

myindx = -1
nf     = -1

! from the dart index get the local variables indices

call get_model_variable_indices(index_in, myindx, dim2_j, dim3_k, &
            var_id=var_id, dom_id=dom_id, kind_index=variable_type)

if( myindx == -1 ) then
     write(string1,*) 'Problem, cannot find base_offst, index_in is: ', index_in
     call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
endif

! so by now, 'myindx' can range from 1-to-nVertices or 1-to-nCells,
! depending on the variable. Since the first nCells parts of coord_nod3D are
! identical to coord_nod2D, we can just index everything with coord_nod3D

lon   = coord_nod3D(1,myindx)
lat   = coord_nod3D(2,myindx)
depth = coord_nod3D(3,myindx)
location = set_location(lon,lat,depth,VERTISHEIGHT)

if (present(var_type)) then
   var_type = variable_type

   if( var_type == MISSING_I ) then
      write(string1,*) 'Cannot find DART QTY for indx ', index_in
      write(string2,*) 'variable "'//trim(get_variable_name(dom_id, var_id))//'"'
      call error_handler(E_ERR, 'get_state_meta_data', string1, &
                 source, revision, revdate, text2=string2)
   endif
endif

end subroutine get_state_meta_data

!------------------------------------------------------------------

subroutine model_interpolate(state_handle, ens_size, location, obs_type, expected_obs, istatus)

! given a state vector, a location, and a QTY_xxx, return the
! value at from the closest grid location and an error code.  0 is success,
! anything positive is an error.  (negative reserved for system use)
! TODO: possibly interpolate ... 
!
!       ERROR codes:
!
!       ISTATUS = 99:  general error in case something terrible goes wrong...
!       ISTATUS = 88:  this kind is not in the state vector
!       ISTATUS = 11:  Could not find a triangle that contains this lat/lon
!       ISTATUS = 12:  Height vertical coordinate out of model range.
!       ISTATUS = 13:  Missing value in interpolation.
!       ISTATUS = 16:  Don't know how to do vertical velocity for now
!       ISTATUS = 17:  Unable to compute pressure values
!       ISTATUS = 18:  altitude illegal
!       ISTATUS = 19:  could not compute u using RBF code
!       ISTATUS = 101: Internal error; reached end of subroutine without
!                      finding an applicable case.
!

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
logical  :: goodkind
real(r8) :: lpres(ens_size), values(3, ens_size)
real(r8) :: llv(3)    ! lon/lat/vert
integer(i8) :: index1  ! the DART index of the start of the variable of interest
integer  :: verttype

! local storage

integer(i8)  :: surface_index
integer(i8)  :: ilayer, layer_below, layer_above
real(r8) :: lon, lat, vert
real(r8) :: depth_below, depth_above, layer_thick
integer(i8)  :: closest_index_above, closest_index_below
integer(i8)  :: index_above, index_below

type(location_type) :: location_above, location_below

if ( .not. module_initialized ) call static_init_model

expected_obs = MISSING_R8
istatus      = 0          ! must be positive (and integer)

! rename for sanity - we can't change the argument names
! to this subroutine, but this really is a kind.
obs_kind = obs_type

! Make sure the DART state has the type (T,S,U,etc.) that we are asking for.
! If we cannot, simply return and 'fail' with an 88

ivar = get_varid_from_kind(domid, obs_kind)
if (ivar < 1) then
   istatus = 88
   return
endif

! Determine the offset into the DART state vector for this variable
index1 = get_index_start(domid,ivar)

! Decode the location into bits for error messages ...
llv  = get_location(location)
lon  = llv(1)    ! degrees East [0,360)
lat  = llv(2)    ! degrees North [-90,90]
vert = llv(3)    ! depth in meters ... even 2D fields have a value of 0.0

surface_index = find_closest_surface_location(location, obs_kind)

if (surface_index < 1) then ! nothing close
   istatus = 11
   return
endif

! If it is a surface variable, we're done.

if (progvar(ivar)%varsize == nCells) then
   index_above   = index1 + surface_index - 1
   expected_obs  = get_state(index_above, state_handle)
   istatus     = 0
   return
endif

! Which vertical level is closest
! Find the first nominal layer deeper than than the observation.
! 2 <= layer_below <= nVertLevels
! use the nod3D_below_nod2D(nVertLevels,nCells)  array to figure out
! 1 <= surface_index <= nCells

layer_below = 0
LAYER: do ilayer = 1,nVertLevels
   if (depths(ilayer) > vert ) then
        layer_below = ilayer
        exit LAYER
   endif
enddo LAYER

if     (layer_below == 0) then ! below the deepest level, RETURN
   istatus = 19
   return
elseif (layer_below == 1) then ! too shallow, RETURN
   istatus = 18
   return
endif

! somewhere in the water column

layer_above = layer_below - 1

! If there is no water, the return value is a negative number

closest_index_above = nod3d_below_nod2d(layer_above,surface_index)
closest_index_below = nod3d_below_nod2d(layer_below,surface_index)

if ((closest_index_below < 1) .or.  (closest_index_above < 1)) then
   istatus = 17
   return
endif

! observation must be 'wet' as far as the model resolution is concerned

index_above = index1 + closest_index_above - 1
index_below = index1 + closest_index_below - 1

depth_below = depths(layer_below) - vert
depth_above = vert - depths(layer_above)
layer_thick = depths(layer_below) - depths(layer_above)

expected_obs  = ( &
                depth_below*get_state(index_above,state_handle) &
              + depth_above*get_state(index_below,state_handle)) &
              / layer_thick
istatus     = 0

! DEBUG block to confirm that the interpolation is using the state
! at the correct location.

if (do_output() .and. debug > 2) then
   call get_state_meta_data(index_above, location_above)
   call get_state_meta_data(index_below, location_below)

   call write_location(0,location_above,charstring=string1)
   call write_location(0,location      ,charstring=string2)
   call write_location(0,location_below,charstring=string3)

   write(logfileunit,*)
   write(     *     ,*)
   call error_handler(E_MSG,'model_interpolate', '... '//string1, &
              text2=string2, text3=string3)
endif

end subroutine model_interpolate


!------------------------------------------------------------------
!>

!> FIXME -aLi- subroutine model_interpolate(x, location, obs_type, interp_val, istatus)
!> 
!> ! given a state vector, a location, and a KIND_xxx, return the
!> ! interpolated value at that location, and an error code.  0 is success,
!> ! anything positive is an error.  (negative reserved for system use)
!> !
!> ! This version simply returns the value at the closest node.
!> !
!> !       ERROR codes:
!> !
!> !       ISTATUS = 99:  general error in case something terrible goes wrong...
!> !       ISTATUS = 88:  this kind is not in the state vector
!> !       ISTATUS = 11:  Could not find a triangle that contains this lat/lon
!> !       ISTATUS = 12:  depth vertical coordinate out of model range.
!> !       ISTATUS = 13:  Missing value in interpolation.
!> !       ISTATUS = 16:  Don't know how to do vertical velocity for now
!> !       ISTATUS = 17:  Unable to compute pressure values
!> !       ISTATUS = 18:  observation too shallow
!> !       ISTATUS = 19:  observation too deep
!> !       ISTATUS = 101: Internal error; reached end of subroutine without
!> !                      finding an applicable case.
!> !
!> 
!> ! passed variables
!> 
!> real(r8),            intent(in)  :: x(:)
!> type(location_type), intent(in)  :: location
!> integer,             intent(in)  :: obs_type
!> real(r8),            intent(out) :: interp_val
!> integer,             intent(out) :: istatus
!> 
!> ! local storage
!> 
!> integer  :: ivar, obs_kind, surface_index
!> integer  :: ilayer, layer_below, layer_above
!> real(r8) :: llv(3), lon, lat, vert
!> real(r8) :: depth_below, depth_above, layer_thick
!> integer  :: closest_index_above, closest_index_below
!> integer  :: index_above, index_below
!> 
!> type(location_type) :: location_above, location_below
!> 
!> if ( .not. module_initialized ) call static_init_model
!> 
!> interp_val = MISSING_R8
!> istatus    = 99           ! must be positive (and integer)
!> 
!> ! rename for sanity - we can't change the argument names
!> ! to this subroutine, but this really is a kind.
!> obs_kind = obs_type
!> 
!> ! Make sure the DART state has the type (T,S,U,etc.) that we are asking for.
!> ! If we cannot, simply return and 'fail' with an 88
!> 
!> ivar = get_progvar_index_from_kind(obs_kind)
!> if (ivar < 1) then
!>    istatus = 88
!>    return
!> endif
!> 
!> ! Decode the location into bits for error messages ...
!> llv  = get_location(location)
!> lon  = llv(1)    ! degrees East [0,360)
!> lat  = llv(2)    ! degrees North [-90,90]
!> vert = llv(3)    ! depth in meters ... even 2D fields have a value of 0.0
!> 
!> surface_index = find_closest_surface_location(location, obs_kind)
!> 
!> if (surface_index < 1) then ! nothing close
!>    istatus = 11
!>    return
!> endif
!> 
!> ! If it is a surface variable, we're done.
!> 
!> if (progvar(ivar)%varsize == nCells) then
!>    index_above = progvar(ivar)%index1 + surface_index - 1
!>    interp_val  = x( index_above )
!>    istatus     = 0
!>    return
!> endif
!> 
!> ! Which vertical level is closest
!> ! Find the first nominal layer deeper than than the observation.
!> ! 2 <= layer_below <= nVertLevels
!> ! use the nod3D_below_nod2D(nVertLevels,nCells)  array to figure out
!> ! 1 <= surface_index <= nCells
!> 
!> layer_below = 0
!> LAYER: do ilayer = 1,nVertLevels
!>    if (depths(ilayer) > vert ) then
!>         layer_below = ilayer
!>         exit LAYER
!>    endif
!> enddo LAYER
!> 
!> if     (layer_below == 0) then ! below the deepest level
!>    istatus = 19
!> elseif (layer_below == 1) then ! too shallow
!>    istatus = 18
!> else                           ! somewhere in the water column
!> 
!>    layer_above = layer_below - 1
!> 
!>    ! If there is no water, the return value is a negative number
!> 
!>    closest_index_above = nod3d_below_nod2d(layer_above,surface_index)
!>    closest_index_below = nod3d_below_nod2d(layer_below,surface_index)
!> 
!>    if ((closest_index_below < 1) .or.  (closest_index_above < 1)) then
!>       istatus = 17
!>       return
!>    endif
!> 
!>    ! observation must be 'wet' as far as the model resolution is concerned
!> 
!>    index_above = progvar(ivar)%index1 + closest_index_above - 1
!>    index_below = progvar(ivar)%index1 + closest_index_below - 1
!> 
!>    depth_below = depths(layer_below) - vert
!>    depth_above = vert - depths(layer_above)
!>    layer_thick = depths(layer_below) - depths(layer_above)
!> 
!>    interp_val  = (depth_below*x(index_above) + depth_above*x(index_below)) &
!>                  / layer_thick
!>    istatus     = 0
!> 
!>    ! DEBUG block to confirm that the interpolation is using the state
!>    ! at the correct location.
!> 
!>    if (do_output() .and. debug > 2) then
!>       call get_state_meta_data(index_above, location_above)
!>       call get_state_meta_data(index_below, location_below)
!> 
!>       call write_location(0,location_above,charstring=string1)
!>       call write_location(0,location      ,charstring=string2)
!>       call write_location(0,location_below,charstring=string3)
!> 
!>       write(logfileunit,*)
!>       write(     *     ,*)
!>       call error_handler(E_MSG,'model_interpolate', '... '//string1, &
!>                  text2=string2, text3=string3)
!>    endif
!> 
!> endif
!> 
!> end subroutine model_interpolate


!-----------------------------------------------------------------------
!>

subroutine nc_write_model_atts( ncFileID, domain_id )

! TJH -- Writes the model-specific attributes to a netCDF file.
!     This includes coordinate variables and some metadata, but NOT
!     the model state vector.
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...

integer, intent(in) :: domain_id
integer, intent(in)  :: ncFileID

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

!----------------------------------------------------------------------
! variables if we just blast out one long state vector
!----------------------------------------------------------------------

integer :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)

!----------------------------------------------------------------------
! variables if we parse the state vector into prognostic variables.
!----------------------------------------------------------------------

! for the dimensions and coordinate variables
integer :: nodes_3DimID
integer :: nodes_2DimID
integer :: nVertLevelsDimID

! for the prognostic variables
integer :: ivar, VarID

!----------------------------------------------------------------------
! local variables
!----------------------------------------------------------------------

! we are going to need these to record the creation date in the netCDF file.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

character(len=NF90_MAX_NAME) :: str1
character(len=NF90_MAX_NAME) :: varname
integer, dimension(NF90_MAX_VAR_DIMS) :: mydimids
integer :: io, myndims

character(len=256) :: filename

if ( .not. module_initialized ) call static_init_model

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file,
! and then put into define mode.
!-------------------------------------------------------------------------------

io = nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID)
call nc_check(io, 'nc_write_model_atts', 'inquire')

call nc_begin_define_mode(ncFileID)

call nc_add_global_creation_time(ncFileID)
call nc_add_global_attribute(ncFileID, "model", "FESOM")

!----------------------------------------------------------------------------
! Define useful geometry variables.
!----------------------------------------------------------------------------

if (diagnostic_metadata) then

   call nc_define_dimension(ncFileID, 'nodes_3d', nVertices, 'nc_write_model_atts')
   call nc_define_dimension(ncFileID, 'nodes_2d', nCells,    'nc_write_model_atts')

   io = nf90_inq_dimid(ncFileID, 'nodes_3d', nodes_3DimID)
   call nc_check(io,'nc_write_model_atts', 'inquire nodes_3d')

   io = nf90_inq_dimid(ncFileID, 'nodes_2d', nodes_2DimID)
   call nc_check(io,'nc_write_model_atts', 'inquire nodes_2d')

   io = nf90_def_var(ncFileID, 'longitudes', NF90_DOUBLE, (/ nodes_3DimID /), VarID)
   call nc_check(io,'nc_write_model_atts', 'longitudes def_var')
   call nc_check(nf90_put_att(ncFileID,VarID,'units','degrees East'),&
                 'nc_write_model_atts', 'longitude units')

   io = nf90_def_var(ncid=ncFileID, name='latitudes', &
                    xtype=NF90_DOUBLE, dimids = (/ nodes_3DimID /), varid=VarID)
   call nc_check(io,'nc_write_model_atts', 'latitudes def_var')
   call nc_check(nf90_put_att(ncFileID,VarID,'units','degrees North'),&
                 'nc_write_model_atts', 'latitudes units')

   io = nf90_def_var(ncid=ncFileID, name='depths', &
                    xtype=NF90_REAL, dimids = (/ nodes_3DimID /), varid=VarID)
   call nc_check(io,'nc_write_model_atts', 'depths def_var')
   call nc_check(nf90_put_att(ncFileID,VarID,'units','meters'),&
                 'nc_write_model_atts', 'depths units')

!TJH    io = nf90_def_var(ncid=ncFileID, name='node_table', &
!TJH              xtype=NF90_INT, dimids = (/ nVertLevelsDimID, nodes_2DimID /), varid=VarID)
!TJH    call nc_check(io,'nc_write_model_atts', 'node_table def_var')
!TJH    call nc_check(nf90_put_att(ncFileID,VarID,'long_name','nod3D_below_nod2D'),&
!TJH                  'nc_write_model_atts', 'node_table long_name')
!TJH    call nc_check(nf90_put_att(ncFileID,VarID,'valid_range',(/ 1, nVertices /)),&
!TJH                  'nc_write_model_atts', 'node_table long_name')
!TJH    call nc_check(nf90_put_att(ncFileID,VarID,'description', &
!TJH                 'table defining packing order. Given a level and a &
!TJH            &horizontal cell, return the vertex index (between 1 and nodes_3d). &
!TJH            &A value outside the valid_range means there is no wet location.'),&
!TJH                  'nc_write_model_atts', 'node_table description')

endif

!----------------------------------------------------------------------------
! Finished with dimension/variable definitions, must end 'define' mode to fill.
!----------------------------------------------------------------------------

call nc_end_define_mode(ncFileID, 'nc_write_model_atts', 'enddef')

!----------------------------------------------------------------------------
! Fill the coordinate variables that DART needs and has locally
!----------------------------------------------------------------------------

if (diagnostic_metadata) then

   call nc_check(NF90_inq_varid(ncFileID, 'longitudes', VarID), &
                 'nc_write_model_atts', 'longitudes inq_varid')
   call nc_check(nf90_put_var(ncFileID, VarID, coord_nod3D(1,:) ), &
                'nc_write_model_atts', 'longitudes put_var')

   call nc_check(NF90_inq_varid(ncFileID, 'latitudes', VarID), &
                 'nc_write_model_atts', 'latitudes inq_varid')
   call nc_check(nf90_put_var(ncFileID, VarID, coord_nod3D(2,:) ), &
                'nc_write_model_atts', 'latitudes put_var')

   call nc_check(NF90_inq_varid(ncFileID, 'depths', VarID), &
                 'nc_write_model_atts', 'depths inq_varid')
   call nc_check(nf90_put_var(ncFileID, VarID, coord_nod3D(3,:) ), &
                'nc_write_model_atts', 'depths put_var')

!TJH    call nc_check(NF90_inq_varid(ncFileID, 'node_table', VarID), &
!TJH                  'nc_write_model_atts', 'node_table inq_varid')
!TJH    call nc_check(nf90_put_var(ncFileID, VarID, nod3D_below_nod2D ), &
!TJH                 'nc_write_model_atts', 'node_table put_var')

endif


!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')


end subroutine nc_write_model_atts


!------------------------------------------------------------------
!>

function get_model_size()

! Returns the size of the model as an integer.
! Required for all applications.

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

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
!>

function get_model_time_step()

! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = model_timestep

end function get_model_time_step


!------------------------------------------------------------------
!>
function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = model_timestep

end function shortest_time_between_assimilations

!------------------------------------------------------------------


subroutine end_model()

! Does any shutdown and clean-up needed for model.

if (allocated(cell_locations))  deallocate(cell_locations)
if (allocated(cell_kinds))      deallocate(cell_kinds)
if (allocated(close_cell_inds)) deallocate(close_cell_inds)
if (allocated(depths))          deallocate(depths)

!> if (close_structure_allocated) call finalize_closest_center()

end subroutine end_model


!------------------------------------------------------------------
!>

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

num_variables = get_num_variables(domid)

! Get min and max of each variable in each domain
allocate(var_list(get_my_num_vars(ens_handle)))
call get_my_vars(ens_handle, var_list)

allocate(min_var(num_variables), max_var(num_variables))
allocate(within_range(ens_handle%my_num_vars))

do i = 1, get_num_variables(domid)

   start_ind = get_index_start(domid, i)
   end_ind = get_index_end(domid, i)

   within_range = (var_list >= start_ind .and. var_list <= end_ind)
   min_var(i) = minval(ens_handle%copies(1,:), MASK=within_range)
   max_var(i) = maxval(ens_handle%copies(1,:), MASK=within_range)

enddo

! get global min/max for each variable
call broadcast_minmax(min_var, max_var, num_variables)
deallocate(within_range)

call init_random_seq(random_seq, my_task_id()+1)

do i = 1, num_variables

   start_ind = get_index_start(domid, i)
   end_ind = get_index_end(domid, i)

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
!>

subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, state_handle)
!
! Given a DART location (referred to as "base") and a set of candidate
! locations & kinds (obs, obs_kind), returns the subset close to the
! "base", their indices, and their distances to the "base" ...

! Note that both base_obs_loc and obs_loc are intent(inout), meaning that these
! locations are possibly modified here and returned as such to the calling routine.
! The calling routine is always filter_assim and these arrays are local arrays
! within filter_assim. In other words, these modifications will only matter within
! filter_assim, but will not propagate backwards to filter.

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:), loc_types(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: state_handle


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
if (present(dist)) dist = 1.0e9_r8   !something big and positive (far away) in radians
istatus1  = 0
istatus2  = 0

! If you want to impose some sort of special localization, you can key
! off things like the obs_kind and make things 'infinitely' far away.
! Otherwise, this does nothing. Take a look at the POP model_mod.f90 for an example.

! Convert base_obs vertical coordinate to requested vertical coordinate if necessary

base_llv = get_location(base_loc)
base_which = nint(query_location(base_loc))

ztypeout = vert_localization_coord

if (vertical_localization_on()) then
  if (base_llv(3) == MISSING_R8) then
     istatus1 = 1
  else if (base_which /= vert_localization_coord) then
!>      call vert_convert(state_handle, base_loc, base_type, istatus1)
      if(debug > 5) then
         call write_location(0,base_loc,charstring=string1)
         call error_handler(E_MSG, 'get_close_obs: base_loc',string1,source, revision, revdate)
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
      if (vertical_localization_on()) then
          if (local_obs_which /= vert_localization_coord) then
!>              call vert_convert(state_handle, local_obs_loc, loc_qtys(t_ind), istatus2)
              locs(t_ind) = local_obs_loc
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
               dist(k) = get_dist(base_loc, local_obs_loc, base_type, loc_qtys(t_ind))
         ! if ((debug > 4) .and. (k < 100) .and. do_output()) then
         !     print *, 'calling get_dist'
         !     call write_location(0,base_loc,charstring=string2)
         !     call error_handler(E_MSG, 'get_close_obs: base_loc',string2,source, revision, revdate)
         !     call write_location(0,local_obs_loc,charstring=string2)
         !     call error_handler(E_MSG, 'get_close_obs: local_obs_loc',string2,source, revision, revdate)
         !     hor_dist = get_dist(base_loc, local_obs_loc, base_type, loc_qtys(t_ind), no_vert=.true.)
         !     print *, 'hor/3d_dist for k =', k, ' is ', hor_dist,dist(k)
         ! endif
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

integer :: istatus1, istatus2
real(r8) :: hor_dist

hor_dist = 1.0e9_r8

call error_handler(E_ERR,'get_close_state','routine not written')

! Initialize variables to missing status

num_close = 0
close_ind = -99
if (present(dist)) dist = 1.0e9_r8   !something big and positive (far away) in radians
istatus1  = 0
istatus2  = 0


end subroutine get_close_state

!------------------------------------------------------------------
!>

subroutine convert_vertical_obs(state_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, status)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:), loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: status(:)

call error_handler(E_ERR,'convert_vertical_obs','routine not written, not needed?')

end subroutine convert_vertical_obs

!--------------------------------------------------------------------

subroutine convert_vertical_state(state_handle, num, locs, loc_qtys, loc_indx, &
                                  which_vert, istatus)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer(i8),         intent(in)    :: loc_indx(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus

call error_handler(E_ERR,'convert_vertical_state','routine not written, not needed?')

end subroutine convert_vertical_state


!==================================================================
! The (model-specific) additional public interfaces come next
!  (these are not required by dart but are used by other programs)
!==================================================================


!------------------------------------------------------------------
!>

function get_analysis_time( ncid, filename )

! The analysis netcdf files have the start time of the experiment.
! The time array contains the time trajectory since then.
! This routine returns the start time of the experiment.

integer,          intent(in) :: ncid
character(len=*), intent(in) :: filename
type(time_type)              :: get_analysis_time

! local variables
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, idims
integer           :: VarID, numdims
integer           :: model_year
type(time_type)   :: model_advance_forecast, model_start_time

real(r8), allocatable :: timearray(:)

if ( .not. module_initialized ) call static_init_model

call nc_check( nf90_inq_varid(ncid, 'time', VarID), &
              'get_analysis_time', 'inquire time '//trim(filename))

call nc_check( nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
              'get_analysis_time', 'inquire TIME '//trim(filename))

if (numdims /= 1) then
   write(string1,*) 'time variable has unknown shape in ', trim(filename)
   call error_handler(E_ERR,'get_analysis_time',string1,source,revision,revdate)
endif

call nc_check( nf90_inquire_dimension(ncid, dimIDs(1), len=idims(1)), &
                 'get_analysis_time', 'inquire time dimension length '//trim(filename))

if (idims(1) /= 1) then
   write(string1,*) 'multiple timesteps (',idims(1),') in file ', trim(filename)
   write(string2,*) 'We are using the LAST one, presumably, the LATEST timestep.'
   call error_handler(E_MSG,'get_analysis_time',string1,source,revision,revdate,text2=string2)
endif

allocate(timearray(idims(1)))
! Get the highest ranking time ... the last one, basically.

call nc_check( nf90_get_var(ncid, VarID, timearray), &
              'get_analysis_time', 'get_var time '//trim(filename))

model_year = year_from_filename(filename)

model_advance_forecast = set_time(int(sum(timearray)-timearray(idims(1))), 0)

model_start_time = set_date(model_year,1,1)

get_analysis_time = model_start_time + model_advance_forecast
if (do_output() .and. debug > 0) then
   call print_date(get_analysis_time, 'get_analysis_time:model date')
   call print_time(get_analysis_time, 'get_analysis_time:model time')
endif

deallocate(timearray)

end function get_analysis_time


!------------------------------------------------------------------
!>

function year_from_filename(filename)

! The analysis netcdf files have the start time of the experiment.
! The time array contains the time trajectory since then.
! This routine returns the start time of the experiment.

integer :: year_from_filename

character(len=*), intent(in) :: filename

integer :: i

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'year_from_filename',string1,source,revision,revdate)
endif

! find the first "period" and use that as the start of the string conversion
i = scan(filename, ".")
if (i <= 0) then
   write(string1,*) 'cannot find time string in name ', trim(filename)
   call error_handler(E_ERR,'year_from_filename',string1,source,revision,revdate)
endif

year_from_filename = string_to_time(filename(i+1:i+4))

end function year_from_filename


!--------------------------------------------------------------------
!> get the model time from the ****.clock   file

function read_model_time(filename)

!  There is a filename with the following format:
!  previous-seconds     day-of-year      year
!  current-seconds      day-of-year      year

character(len=256), intent(in) :: filename
type(time_type) :: read_model_time

real(r8) :: fseconds
integer  :: seconds, dayofyear, year
integer  :: io, iunit
type(time_type) :: january1, forecast

! ALI ... make a namelist variable specifying the name of the file we 
! SHOULD be reading ... and open it here instead of argument.

iunit = open_file(filename, action = 'read')

! skip the first record
read(iunit,*,iostat=io) fseconds, dayofyear, year
if (io /= 0) call error_handler(E_ERR,'read_model_time', &
             'bomb reading "'//trim(filename)//'"')

! read the time of interest
read(iunit,*,iostat=io) fseconds, dayofyear, year
if (io /= 0) call error_handler(E_ERR,'read_model_time','BOMB')

january1 = set_date(year,1,1)

seconds = nint(fseconds)

read_model_time = january1 + set_time(seconds,dayofyear) - set_time(0,1)

call print_time(read_model_time,'read_model_time:')
call print_date(read_model_time,'read_model_time:')

call close_file(iunit)

end function read_model_time

!--------------------------------------------------------------------
!> not written

subroutine write_model_time(ncid, dart_time)

integer,         intent(in) :: ncid
type(time_type), intent(in) :: dart_time

call error_handler(E_ERR,'write_model_time','routine not written')

end subroutine write_model_time

!-----------------------------------------------------------------------

!TJH subroutine write_model_time_file(time_filename, model_time, adv_to_time)
!TJH  character(len=*), intent(in)           :: time_filename
!TJH  type(time_type),  intent(in)           :: model_time
!TJH  type(time_type),  intent(in), optional :: adv_to_time
!TJH 
!TJH integer :: iunit
!TJH character(len=19) :: timestring
!TJH type(time_type)   :: deltatime
!TJH 
!TJH iunit = open_file(time_filename, action='write')
!TJH 
!TJH timestring = time_to_string(model_time)
!TJH write(iunit, '(A)') timestring
!TJH 
!TJH if (present(adv_to_time)) then
!TJH    timestring = time_to_string(adv_to_time)
!TJH    write(iunit, '(A)') timestring
!TJH 
!TJH    deltatime = adv_to_time - model_time
!TJH    timestring = time_to_string(deltatime, interval=.true.)
!TJH    write(iunit, '(A)') timestring
!TJH endif
!TJH 
!TJH call close_file(iunit)
!TJH 
!TJH end subroutine write_model_time_file
!TJH 
!TJH 
!TJH !-----------------------------------------------------------------------
!TJH 
!TJH subroutine write_model_time_restart(ncid, dart_time)
!TJH 
!TJH integer,             intent(in) :: ncid !< netcdf file handle
!TJH type(time_type),     intent(in) :: dart_time
!TJH 
!TJH call error_handler(E_MSG, 'write_model_time', 'no routine for fesom write model time')
!TJH 
!TJH end subroutine write_model_time_restart

!> FIXME: -aLi-----------------------------------------------------
!> subroutine write_model_time(time_filename, model_time, adv_to_time)
!>  character(len=*), intent(in)           :: time_filename
!>  type(time_type),  intent(in)           :: model_time
!>  type(time_type),  intent(in), optional :: adv_to_time
!> 
!> integer :: iunit
!> type(time_type) :: deltatime
!> type(time_type) :: new_model_time, assim_start, assim_end
!> integer :: seconds, days
!> 
!> iunit = open_file(time_filename, form='formatted', action='write')
!> 
!> deltatime = set_time(assimilation_period_seconds, assimilation_period_days)
!> 
!> new_model_time = model_time + deltatime
!> assim_start    = new_model_time - deltatime/2 + set_time(1,0)
!> assim_end      = new_model_time + deltatime/2
!> 
!> ! By writing the days,seconds to strings with a free-format write,
!> ! you can avoid any decision about formatting precision.
!> ! The '(A)' syntax avoids printing the quotes delimiting the string
!> 
!> call get_time(new_model_time, seconds, days)
!> write(string1,*) days
!> write(string2,*) seconds
!> write(iunit,'(A,A)') ' init_time_days     = ',trim(string1)
!> write(iunit,'(A,A)') ' init_time_seconds  = ',trim(string2)
!> 
!> call get_time(assim_start, seconds, days)
!> write(string1,*) days
!> write(string2,*) seconds
!> write(iunit,'(A,A)') ' first_obs_days     = ',trim(string1)
!> write(iunit,'(A,A)') ' first_obs_seconds  = ',trim(string2)
!> 
!> call get_time(assim_end, seconds, days)
!> write(string1,*) days
!> write(string2,*) seconds
!> write(iunit,'(A,A)') ' last_obs_days      = ',trim(string1)
!> write(iunit,'(A,A)') ' last_obs_seconds   = ',trim(string2)
!> 
!> string2 = time_to_string(new_model_time)
!> 
!> call error_handler(E_MSG,'write_model_time:','next model time should be '//trim(string2))
!> 
!> write(iunit, '(A)') trim(string2)
!> 
!> call print_time(new_model_time,'FESOM    stop at :',  iunit)
!> call print_time(    model_time,'FESOM current at :',  iunit)
!> call print_date(new_model_time,'stop    date :',  iunit)
!> call print_date(    model_time,'current date :',  iunit)
!> 
!> if (present(adv_to_time)) then
!>    string1 = time_to_string(adv_to_time)
!>    write(iunit, '(A)') trim(string1)
!> 
!>    deltatime = adv_to_time - model_time
!>    string1 = time_to_string(deltatime, interval=.true.)
!>    write(iunit, '(A)') trim(string1)
!> endif
!> 
!> call close_file(iunit)
!> 
!> end subroutine write_model_time


!------------------------------------------------------------------
!>

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
Edges      = missing_I
VertLevels = nVertLevels
VertexDeg  = 60

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
   write(time_to_string, '(I2.2,A1,I2.2)') &
                        ndays, ' ', nsecs
   ihour = nsecs / 3600
   nsecs = nsecs - (ihour * 3600)
   imin  = nsecs / 60
   nsecs = nsecs - (imin * 60)
   isec  = nsecs
   write(time_to_string, '(I2.2,3(A1,I2.2))') &
                        ndays, ' ', ihour, ' ', imin, ' ', isec
else
   call get_date(t, iyear, imonth, iday, ihour, imin, isec)
   write(time_to_string, '(I4.4,5(A1,I2.2))') &
                        iyear, ' ', imonth, ' ', iday, ' ', ihour, ' ', imin, ' ', isec

endif

end function time_to_string


!------------------------------------------------------------------
!>

function string_to_time(s)

character(len=*), intent(in) :: s
integer                      :: string_to_time

integer :: iyear

read(s,'(i4)') iyear

if (iyear < 1601) then
   write(string1,*)'WARNING: Converting YEAR ',iyear,' to ',iyear+1601
   write(string2,*)'original time (string) is <',trim(s),'>'
   call error_handler(E_MSG, 'string_to_time', string1, &
               source, revision, revdate, text2=string2)
   iyear = iyear + 1601
endif

string_to_time = iyear

end function string_to_time


!------------------------------------------------------------------
!>

function set_model_time_step()

! the static_init_model ensures that the model namelists are read.

type(time_type) :: set_model_time_step

if ( .not. module_initialized ) call static_init_model

! these are from the namelist
!FIXME: sanity check these for valid ranges?
set_model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end function set_model_time_step


!------------------------------------------------------------------
!> Read the grid from the FESOM metadata files.
!>@ TODO If this can be moved to some routine that gets called by
!   every task - but not when called by the converter routines,
!   it would be faster.

subroutine read_grid()

! This routine does not need to be called by the model_to_dart and
! dart_to_model routines, but DOES need to be called very early by
! filter and perfect_model_obs.

logical, save :: grid_read = .false.
integer :: i
real(r8) :: maxdist_km

if (grid_read) return ! only read the grid once.

! referenced by module 'use'
! nCells      = myDim_nod2D
! nVertices   = myDim_nod3D
! nVertLevels = max_num_layers

! the grid coordinates (coord_nod3D) come from the FESOM module
! lon   = coord_nod3D(1,:)
! lat   = coord_nod3D(2,:)
! depth = coord_nod3D(3,:)

call read_node()  ! sets myDim_nod2D, myDim_nod3D
call read_aux3()  ! sets nVertLevels and node locating arrays
call read_depth()

! convert to [0, 360.0) longitude base
where(coord_nod3D(1,:) < 0.0_r8) &
      coord_nod3D(1,:) = coord_nod3D(1,:) + 360.0_r8

! convert depths from negative to positive.
! the 1D depths array has positive depths.
! the observations in the WOD have positive depths.
coord_nod3D(3,:) = -1.0_r8 * coord_nod3D(3,:)

!>@ TODO there are a lot of variables allocated in the fesom_modules.f90
!        that are not needed. For memory-efficiency, we should deallocate
!        all those (large) variables. coord_nod2D, index_nod2D, myList_nod2D, ...

allocate( cell_locations(nCells), cell_kinds(nCells), close_cell_inds(nCells) )
allocate(depths(nVertLevels))

depths = real(layerdepth,r8)
cell_kinds = 0 ! not used

! The first nCells locations in coord_nod3D are the same as those in coord_nod2D
! and define the first layer.

do i=1,nCells
   cell_locations(i) = set_location(coord_nod3D(1,i), &
                                    coord_nod3D(2,i), &
                                    coord_nod3D(3,i), VERTISHEIGHT)
enddo

! Initialize a 'get_close' structure that sets up the lookup table
! to speed up the identification of the grid location closest to
! any arbitrary location.  Note: there are 40,000km in 2PI radians.

maxdist_km = 2.5_r8 ! more than the largest separation between vertices

!for NANCY ... how do we initialize a get_close structure for _part_ of the state

close_structure_allocated = .true.

if (debug > 5) &
   call error_handler(E_MSG,'read_grid','get close lookup table initialized')

if (debug > 1) then
   write(string1,*)'nCells      is ', nCells
   write(string2,*)'nVertices   is ', nVertices
   write(string3,*)'nVertLevels is ', nVertLevels
   call error_handler(E_MSG,'read_grid','... '//string1,text2=string2,text3=string3)

   write(string1,*)'latitude  range ',minval(coord_nod3D(2,:)),maxval(coord_nod3D(2,:))
   write(string2,*)'longitude range ',minval(coord_nod3D(1,:)),maxval(coord_nod3D(1,:))
   write(string3,*)'depths    range ',minval(depths),          maxval(depths)
   call error_handler(E_MSG,'read_grid','... '//string1,text2=string2,text3=string3)
endif

grid_read = .true.

end subroutine read_grid


!------------------------------------------------------------------
!>

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
              'read_2d_from_nc_file', &
              'inq_varid '//trim(varname)//' '//trim(model_analysis_filename))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
              'read_2d_from_nc_file', &
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
              'get_var '//trim(varname)//' '//trim(model_analysis_filename))

end subroutine read_2d_from_nc_file


!------------------------------------------------------------------
!>

subroutine parse_variable_input( state_variables, ncid, filename, ngood )

character(len=*), intent(in)  :: state_variables(:,:)
integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: filename
integer,          intent(out) :: ngood

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: dimname
integer :: nrows, i, j, VarID, dimlen, numdims
logical :: failure

character(len=NF90_MAX_NAME) :: varname       ! column 1
character(len=NF90_MAX_NAME) :: dartstr       ! column 2
character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
character(len=NF90_MAX_NAME) :: state_or_aux  ! column 5

real(r8) :: minvalue, maxvalue
integer  :: ios

if ( .not. module_initialized ) call static_init_model

var_names  = 'no_variable_specified'
var_qtys   = MISSING_I
var_ranges = MISSING_R8
var_update = .false.

failure = .FALSE. ! perhaps all with go well

ngood = 0
MyLoop : do i = 1, MAX_STATE_VARIABLES

   if ( variables(1,i) == ' ' .and. variables(2,i) == ' ' ) exit MyLoop ! Found end of list.

   if ( any(state_variables(:,i) == ' ') ) then
      string1 = '...  model_nml:"variables" not fully specified'
      write(string2,*)'failing on line ',i
      call error_handler(E_ERR, 'parse_variable_input', string1, &
                 source, revision, revdate, text2=string2)
   endif

   varname      = trim(state_variables(VARNAME_INDEX,i))
   dartstr      = trim(state_variables(   KIND_INDEX,i))
   minvalstring = trim(state_variables( MINVAL_INDEX,i))
   maxvalstring = trim(state_variables( MAXVAL_INDEX,i))
   state_or_aux = trim(state_variables(REPLACE_INDEX,i))
   call to_upper(state_or_aux)

   var_names(i) = trim(varname)
   var_qtys(i)  = get_index_for_quantity(dartstr) 

   read(minvalstring,*,iostat=ios) minvalue
   if (ios == 0) var_ranges(i,1) = minvalue

   read(maxvalstring,*,iostat=ios) maxvalue
   if (ios == 0) var_ranges(i,1) = maxvalue

   if (state_or_aux == 'UPDATE' ) var_update(i) = .true.

   ! Make sure DART kind is valid

   if( var_qtys(i) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') &
            trim(dartstr)
      call error_handler(E_ERR,'parse_variable_input',string1,source,revision,revdate)
   endif

   ! Make sure variable exists in model analysis variable list

   write(string1,'(''variable '',a,'' in '',a)') trim(varname), trim(filename)
   write(string2,'(''there is no '',a)') trim(string1)
   call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
                 'parse_variable_input', trim(string2))

   ! Make sure variable is defined by (Time,nCells) or (Time,nCells,vertical)
   ! unable to support Edges or Vertices at this time.

   call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
                 'parse_variable_input', 'inquire '//trim(string1))

   DimensionLoop : do j = 1,numdims

      write(string2,'(''inquire dimension'',i2,'' of '',a)') j,trim(string1)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(j), len=dimlen, name=dimname), &
                                          'parse_variable_input', trim(string2))
      select case ( trim(dimname) )
         case ('T')
            ! supported - do nothing
         case ('nodes_2d')
            ! supported - do nothing
         case ('nodes_3d')
            ! supported - do nothing
         case default
            write(string2,'(''unsupported dimension '',a,'' in '',a)') trim(dimname),trim(string1)
            call error_handler(E_MSG,'parse_variable_input',string2,source,revision,revdate)
            failure = .TRUE.
      end select

   enddo DimensionLoop

   if (failure) then
       string2 = 'unsupported dimension(s) are fatal'
       call error_handler(E_ERR,'parse_variable_input',string2,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector

   if (debug > 0) then
      write(string1,*)'variable ',i,' is ',trim(varname), ' ', trim(dartstr)
      call error_handler(E_MSG,'parse_variable_input',string1)
   endif

   ngood = ngood + 1
enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''MAX_STATE_VARIABLES'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'parse_variable_input',string1,text2=string2)
endif

end subroutine parse_variable_input


!------------------------------------------------------------------
!>

subroutine dump_progvar(ivar)

! dump the contents of the metadata for an individual entry.
! expected to be called in a loop or called for entries of interest.

integer,  intent(in)           :: ivar
integer                        :: i

! take care of parallel runs where we only want a single copy of
! the output.
if (.not. do_output()) return

write(logfileunit,*)
write(     *     ,*)
write(logfileunit,*) 'variable number ',ivar,' is ',trim(progvar(ivar)%varname)
write(     *     ,*) 'variable number ',ivar,' is ',trim(progvar(ivar)%varname)
write(logfileunit,*) '  replace     ',progvar(ivar)%replace
write(     *     ,*) '  replace     ',progvar(ivar)%replace
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
write(logfileunit,*)
write(     *     ,*)

end subroutine dump_progvar


!------------------------------------------------------------------
!>

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
!>

subroutine print_minmax(ivar, x)

! given an index and a state vector, print out the min and max
! data values for the items corresponding to that progvar index.

integer,  intent(in) :: ivar
real(r8), intent(in) :: x(:)

write(string1, '(A,A32,2F16.7)') 'data  min/max ', trim(progvar(ivar)%varname), &
           minval(x(progvar(ivar)%index1:progvar(ivar)%indexN)), &
           maxval(x(progvar(ivar)%index1:progvar(ivar)%indexN))

call error_handler(E_MSG,'print_minmax', string1, source,revision,revdate)

end subroutine print_minmax


!------------------------------------------------------------------
!>

function FindTimeDimension(ncid) result(timedimid)

! Find the Time Dimension ID in a netCDF file.
! If there is none - (spelled the obvious way) - the routine
! returns a negative number. You don't HAVE to have a TIME dimension.

integer                      :: timedimid
integer,          intent(in) :: ncid

integer :: nc_rc

TimeDimID = -1 ! same as the netCDF library routines.
nc_rc = nf90_inq_dimid(ncid,'T',dimid=TimeDimID)

end function FindTimeDimension


!------------------------------------------------------------
!>

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

   ! Each of these dimension names (originally from the FESOM analysis file)
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
!>

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
!>

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


!==================================================================
! The following (private) interfaces are used for triangle interpolation
!==================================================================


!------------------------------------------------------------------
!>

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
!>

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


!------------------------------------------------------------
!> Determine the cell index closest to the given point
!> 2D calculation only.

function find_closest_surface_location(location, obs_kind)

type(location_type), intent(in) :: location
integer,             intent(in) :: obs_kind
integer                         :: find_closest_surface_location

integer :: num_close, surface_index, iclose, indx, i
real(r8) :: closest
real(r8) :: distance

find_closest_surface_location = -1

! Generate the list of indices into the DART vector of the close candidates.

!TJH call loc_get_close_state(cc_gc, location, obs_kind, cell_locations, cell_kinds, &
!TJH                       num_close, close_ind=close_cell_inds)

!TJH just do an exhaustive search for late on a Friday afternoon
num_close = nCells
close_cell_inds = (/ (i, i=1,num_close) /)

! Sometimes the location is outside the model domain.
! In this case, we cannot interpolate.

if (num_close == 0) return

! Loop over close candidates. They come in without regard to what DART KIND they are,
! nor are they sorted by distance. We are only interested in the close locations
! of the right DART KIND.

closest = 1000000.0_r8 ! not very close
surface_index = 0

CLOSE : do iclose = 1, num_close

   indx     = close_cell_inds(iclose)
   distance = get_dist(location, cell_locations(indx))

   if (distance < closest) then
      closest       = distance
      surface_index = indx
   endif

enddo CLOSE

if (do_output() .and. debug > 2) &
   write(*,*)'HORIZONTALLY closest is state index ',surface_index, ' is at distance ',closest

find_closest_surface_location = surface_index

end function find_closest_surface_location

!------------------------------------------------------------

subroutine dump_tables

! TJH DEBUG create the rosetta stone for locations
! separate file for temp, salinity, coord_nod3D, all dart state_vector.
! should be able to cat temp & salinity to recreate all dart_state-vector
! temp == salt == coord_nod3D, that sort of thing.

integer     :: iunit, var_type
integer(i8) :: i
real(r8)    :: llv(3), lon, lat, vert
type(location_type) :: location

if (do_output() .and. state_table_needed ) then
   state_table_needed = .false. ! only do this once, it is expensive as hell

   write(*,*)'writing dart_state_locations.txt'

   iunit = open_file('dart_state_locations.txt',action='write')
   do i = 1,model_size
      call get_state_meta_data(i,location,var_type)
      llv  = get_location(location)
      lon  = llv(1) * deg2rad
      lat  = llv(2) * deg2rad
      vert = llv(3)
      call write_location(0,location,charstring=string1)
      write(iunit,100) i, var_type, trim(string1), lon, lat, vert
   enddo
   call close_file(iunit)

   write(*,*)'writing salinity_locations.txt'

   iunit = open_file('salinity_locations.txt',action='write')
   do i = progvar(1)%index1,progvar(1)%indexN
      call get_state_meta_data(i,location,var_type)
      llv  = get_location(location)
      lon  = llv(1) * deg2rad
      lat  = llv(2) * deg2rad
      vert = llv(3)
      call write_location(0,location,charstring=string1)
      write(iunit,100) i, var_type, trim(string1), lon, lat, vert
   enddo
   call close_file(iunit)

   write(*,*)'writing temperature_locations.txt'

   iunit = open_file('temperature_locations.txt',action='write')
   do i = progvar(2)%index1,progvar(2)%indexN
      call get_state_meta_data(i,location,var_type)
      llv  = get_location(location)
      lon  = llv(1) * deg2rad
      lat  = llv(2) * deg2rad
      vert = llv(3)
   var_type = 50  ! to match to salinity
      call write_location(0,location,charstring=string1)
      write(iunit,100) i-progvar(2)%index1+1, var_type, trim(string1), lon, lat, vert
   enddo
   call close_file(iunit)

   write(*,*)'writing fesom_state_locations.txt'

   iunit = open_file('fesom_state_locations.txt',action='write')
   do i = 1,nVertices
      var_type = 50  ! to match to salinity
      lon  = coord_nod3D(1,i) * deg2rad
      lat  = coord_nod3D(2,i) * deg2rad
      vert = coord_nod3D(3,i)
      location = set_location(coord_nod3D(1,i), &
                              coord_nod3D(2,i), &
                              coord_nod3D(3,i), VERTISHEIGHT)
      call write_location(0,location,charstring=string1)
      write(iunit,100) i, var_type, trim(string1), lon, lat, vert
   enddo
   call close_file(iunit)

   write(*,*)'writing fesom_surface_locations.txt'

   iunit = open_file('fesom_surface_locations.txt',action='write')
   do i = 1,nCells
      var_type = 50  ! to match to salinity
      lon  = coord_nod2D(1,i) * deg2rad
      lat  = coord_nod2D(2,i) * deg2rad
      vert = coord_nod3D(3,i)  ! deliberately using 3D
      location = set_location(coord_nod2D(1,i), &
                              coord_nod2D(2,i), &
                              coord_nod3D(3,i), VERTISHEIGHT)
      call write_location(0,location,charstring=string1)
      write(iunit,100) i, var_type, trim(string1), lon, lat, vert
   enddo
   call close_file(iunit)

 100  format('index ',i7,1x,i3,1x,A,2(f24.12,1x),f13.7)

endif

end subroutine dump_tables

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
