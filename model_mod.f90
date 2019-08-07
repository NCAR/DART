! DART software - Copyright 2004 - 2015 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: model_mod.f90 10268 2016-05-09 17:06:15Z thoar $

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

use        types_mod, only : r4, r8, digits12, SECPERDAY, MISSING_R8,       &
                             rad2deg, deg2rad, PI, MISSING_I, obstypelength

use time_manager_mod, only : time_type, set_time, set_date, get_date, &
                             get_time, print_time, print_date,        &
                             set_calendar_type, increment_time,       &
                             operator(*),  operator(+), operator(-),  &
                             operator(>),  operator(<), operator(/),  &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, get_dist, query_location,    &
                             get_close_maxdist_init, get_close_type,     &
                             set_location, get_location, write_location, &
                             vert_is_height, VERTISHEIGHT,               &
                             get_close_obs_init, get_close_obs_destroy,  &
                             loc_get_close_obs => get_close_obs

use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper, nmlfileunit,       &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, file_exist, find_textfile_dims,        &
                             file_to_text, close_file, do_nml_file, do_nml_term

use     obs_kind_mod, only : paramname_length,           &
                             get_raw_obs_kind_index,     &
                             get_raw_obs_kind_name,      &
                             KIND_VERTICAL_VELOCITY,     &
                             KIND_POTENTIAL_TEMPERATURE, &
                             KIND_TEMPERATURE,           &
                             KIND_SALINITY,              &
                             KIND_DRY_LAND,              &
                             KIND_EDGE_NORMAL_SPEED,     &
                             KIND_U_CURRENT_COMPONENT,   &
                             KIND_V_CURRENT_COMPONENT,   &
                             KIND_SEA_SURFACE_HEIGHT,    &
                             KIND_SEA_SURFACE_PRESSURE,  &
                             KIND_TRACER_CONCENTRATION

use mpi_utilities_mod, only: my_task_id

use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use      fesom_modules, only: read_node, read_aux3, read_depth, read_namelist, &
                             nCells => myDim_nod2D, & ! number of surface locations
                             nVertices => myDim_nod3D, & ! wet points in grid
                             nVertLevels => max_num_layers, & ! number of vertical levels
                             layerdepth, & ! depth at each level (m)
                             coord_nod2D, &
                             coord_nod3D, &
                             num_layers_below_nod2d, &
                             nod2d_corresp_to_nod3D, &
                             nod3d_below_nod2d

! netcdf modules
use typesizes
use netcdf

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
          analysis_file_to_statevector, &
          statevector_to_analysis_file, &
          write_model_time,             &
          get_grid_dims,                &
          print_variable_ranges

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/FEOM/models/FeoM/model_mod.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 10268 $"
character(len=128), parameter :: revdate  = "$Date: 2016-05-09 19:06:15 +0200 (Mon, 09 May 2016) $"

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

type(location_type), allocatable :: cell_locations(:)
integer,             allocatable :: cell_kinds(:)
integer,             allocatable :: close_cell_inds(:)
real(r8),            allocatable :: depths(:)

type(get_close_type)  :: cc_gc

integer, parameter :: MAX_STATE_VARIABLES = 80
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 5
integer, parameter :: VARNAME_INDEX = 1
integer, parameter ::    KIND_INDEX = 2
integer, parameter ::  MINVAL_INDEX = 3
integer, parameter ::  MAXVAL_INDEX = 4
integer, parameter :: REPLACE_INDEX = 5

! variables which are in the module namelist
integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 21600     ! 86400 | 43200 | 21600
real(r8)           :: model_perturbation_amplitude = 0.0001   ! tiny amounts
logical            :: output_state_vector = .false.  ! output prognostic variables (if .false.)
logical            :: diagnostic_metadata = .false.
integer            :: debug = 0   ! turn up for more and more debug messages
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: model_analysis_filename = 'expno.year.oce.nc'
character(len=NF90_MAX_NAME) :: variables(MAX_STATE_VARIABLES * NUM_STATE_TABLE_COLUMNS ) = ' '

namelist /model_nml/             &
   model_analysis_filename,      &
   output_state_vector,          &
   diagnostic_metadata,          &
   assimilation_period_days,     &
   assimilation_period_seconds,  &
   model_perturbation_amplitude, &
   calendar,                     &
   variables,                    &
   debug

character(len=NF90_MAX_NAME) :: variable_table(MAX_STATE_VARIABLES, NUM_STATE_TABLE_COLUMNS )

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
   character(len=paramname_length) :: kind_string
   logical  :: clamping     ! does variable need to be range-restricted before
   real(r8) :: range(2)     ! being stuffed back into FESOM analysis file.
   logical  :: out_of_range_fail  ! is out of range fatal if range-checking?
   logical  :: replace      ! does variable need to be stuffed back info FESOM
end type progvartype

type(progvartype), dimension(MAX_STATE_VARIABLES) :: progvar

! Grid parameters - the values will be read from an FESOM analysis file.

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

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname,dimname
character(len=paramname_length)       :: kind_string
integer :: ncid, VarID, numdims, varsize, dimlen
integer :: iunit, io, ivar, i, index1, indexN
integer :: ss, dd
integer :: nDimensions, nVariables, nAttributes, unlimitedDimID, TimeDimID
real(r8) :: lower_bound, upper_bound

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

call nc_check( nf90_open(trim(model_analysis_filename), NF90_NOWRITE, ncid), &
                  'static_init_model', 'open '//trim(model_analysis_filename))

call parse_variable_input( variables, ncid, model_analysis_filename, &
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

   varname                   = trim(variable_table(ivar,VARNAME_INDEX))
   kind_string               = trim(variable_table(ivar,   KIND_INDEX))
   progvar(ivar)%varname     = varname
   progvar(ivar)%kind_string = kind_string
   progvar(ivar)%dart_kind   = get_raw_obs_kind_index( progvar(ivar)%kind_string )
   progvar(ivar)%numdims     = 0
   progvar(ivar)%numvertical = 1
   progvar(ivar)%dimlens     = MISSING_I
   progvar(ivar)%numcells    = MISSING_I
   progvar(ivar)%replace     = .true.
   progvar(ivar)%clamping    = .false.
   progvar(ivar)%out_of_range_fail = .false.  ! FIXME ... not used

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

      select case ( dimname(1:8) )
         case ('nodes_2d')
            progvar(ivar)%numcells = nCells
         case ('nodes_3d')
            progvar(ivar)%numcells    = nVertices
            progvar(ivar)%numvertical = nVertLevels
      end select

   enddo DimensionLoop

   progvar(ivar)%varsize = varsize
   progvar(ivar)%index1  = index1
   progvar(ivar)%indexN  = index1 + varsize - 1
   index1                = index1 + varsize      ! sets up for next variable

   ! resolve issues with bounded variables.

   read(variable_table(ivar,MINVAL_INDEX),*,iostat=io)lower_bound
   if (io /= 0) lower_bound = MISSING_R8

   read(variable_table(ivar,MAXVAL_INDEX),*,iostat=io)upper_bound
   if (io /= 0) upper_bound = MISSING_R8

   progvar(ivar)%range = (/ lower_bound, upper_bound /)

   if (lower_bound /= MISSING_R8 .or. upper_bound /= MISSING_R8) then
      progvar(ivar)%clamping = .true.
   endif

   ! resolve issues with variables that may not need to be reinserted into FESOM

   string1 = adjustl(variable_table(ivar,REPLACE_INDEX))
   call to_upper(string1)

   if (string1(1:2) /= 'UP') then
      progvar(ivar)%replace = .false.
   endif

   ! print summary if desired.

   if ( debug > 1 ) call dump_progvar(ivar)

enddo

call nc_check( nf90_close(ncid), &
                  'static_init_model', 'close '//trim(model_analysis_filename))

model_size = progvar(nfields)%indexN

!>@ TODO try to move read_grid to some other routine that gets called very
!   early by all the tasks.

call read_grid() ! sets nCells, nVertices, nVertLevels

! dump_tables() is VERY expensive and should only be called
! once to generate sanity checks.
if (state_table_needed .and. debug > 99) call dump_tables()

end subroutine static_init_model


!------------------------------------------------------------------
!>

subroutine get_state_meta_data(index_in, location, var_type)

! given an index into the state vector, return its location and
! if given, the var kind.   despite the name, var_type is a generic
! kind, like those in obs_kind/obs_kind_mod.f90, starting with KIND_

! passed variables

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

! Local variables

integer  :: nf, n, myindx
real(r8) :: lon, lat, depth

if ( .not. module_initialized ) call static_init_model

myindx = -1
nf     = -1

! Determine the right variable
FindIndex : do n = 1,nfields
    if( (progvar(n)%index1 <= index_in) .and. (index_in <= progvar(n)%indexN) ) then
      nf = n
      myindx = index_in - progvar(n)%index1 + 1
      exit FindIndex
    endif
enddo FindIndex

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
   var_type = progvar(nf)%dart_kind
endif

end subroutine get_state_meta_data


!------------------------------------------------------------------
!>

subroutine model_interpolate(x, location, obs_type, interp_val, istatus)

! given a state vector, a location, and a KIND_xxx, return the
! interpolated value at that location, and an error code.  0 is success,
! anything positive is an error.  (negative reserved for system use)
!
! This version simply returns the value at the closest node.
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
!       ISTATUS = 18:  observation too shallow
!       ISTATUS = 19:  observation too deep
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

integer  :: ivar, obs_kind, surface_index
integer  :: ilayer, layer_below, layer_above
real(r8) :: llv(3), lon, lat, vert
real(r8) :: depth_below, depth_above, layer_thick
integer  :: closest_index_above, closest_index_below
integer  :: index_above, index_below

type(location_type) :: location_above, location_below

if ( .not. module_initialized ) call static_init_model

interp_val = MISSING_R8
istatus    = 99           ! must be positive (and integer)

! rename for sanity - we can't change the argument names
! to this subroutine, but this really is a kind.
obs_kind = obs_type

! Make sure the DART state has the type (T,S,U,etc.) that we are asking for.
! If we cannot, simply return and 'fail' with an 88

ivar = get_progvar_index_from_kind(obs_kind)
if (ivar < 1) then
   istatus = 88
   return
endif

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
   index_above = progvar(ivar)%index1 + surface_index - 1
   interp_val  = x( index_above )
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

if     (layer_below == 0) then ! below the deepest level
   istatus = 19
elseif (layer_below == 1) then ! too shallow
   istatus = 18
else                           ! somewhere in the water column

   layer_above = layer_below - 1

   ! If there is no water, the return value is a negative number

   closest_index_above = nod3d_below_nod2d(layer_above,surface_index)
   closest_index_below = nod3d_below_nod2d(layer_below,surface_index)

   if ((closest_index_below < 1) .or.  (closest_index_above < 1)) then
      istatus = 17
      return
   endif

   ! observation must be 'wet' as far as the model resolution is concerned

   index_above = progvar(ivar)%index1 + closest_index_above - 1
   index_below = progvar(ivar)%index1 + closest_index_below - 1

   depth_below = depths(layer_below) - vert
   depth_above = vert - depths(layer_above)
   layer_thick = depths(layer_below) - depths(layer_above)

   interp_val  = (depth_below*x(index_above) + depth_above*x(index_below)) &
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

endif

end subroutine model_interpolate


!-----------------------------------------------------------------------
!>

function nc_write_model_atts( ncFileID ) result (ierr)

! TJH -- Writes the model-specific attributes to a netCDF file.
!     This includes coordinate variables and some metadata, but NOT
!     the model state vector.
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
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'FESOM' ), &
           'nc_write_model_atts', 'model put '//trim(filename))

!-------------------------------------------------------------------------------
! Here is the extensible part. The simplest scenario is to output the state vector,
! parsing the state vector into model-specific parts is complicated, and you need
! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
! complicated part.
!-------------------------------------------------------------------------------

if ( output_state_vector ) then

   !----------------------------------------------------------------------------
   ! Define the DART model size and
   ! create a variable for the state vector.
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable', len=model_size, &
        dimid = StateVarDimID),'nc_write_model_atts', 'state def_dim '//trim(filename))

   ! Define the actual (3D) state vector, which gets filled as time goes on ...
   call nc_check(nf90_def_var(ncid=ncFileID, name='state', xtype=nf90_real, &
                 dimids=(/StateVarDimID,MemberDimID,unlimitedDimID/),varid=VarID),&
                 'nc_write_model_atts','state def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,VarID,'long_name','model state or fcopy'),&
                 'nc_write_model_atts', 'state long_name '//trim(filename))

   ! Leave define mode.
   call nc_check(nf90_enddef(ncFileID),'nc_write_model_atts','state enddef '//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to output the prognostic variables.
   !----------------------------------------------------------------------------
   ! Define the new dimensions IDs
   !----------------------------------------------------------------------------

   io = nf90_def_dim(ncid=ncFileID, name='nodes_2d',len=nCells, dimid=nodes_2DimID)
   call nc_check(io, 'nc_write_model_atts', 'nodes_2D def_dim '//trim(filename))

   io = nf90_def_dim(ncid=ncFileID, name='nodes_3d', len=nVertices, dimid=nodes_3DimID)
   call nc_check(io,'nc_write_model_atts', 'nodes_3D def_dim '//trim(filename))

   io = nf90_def_dim(ncid=ncFileID, name='levels', len=nVertLevels, dimid=nVertLevelsDimID)
   call nc_check(io,'nc_write_model_atts', 'levels def_dim '//trim(filename))

   !----------------------------------------------------------------------------
   ! Define useful geometry variables.
   !----------------------------------------------------------------------------

   if (diagnostic_metadata) then

      io = nf90_def_var(ncid=ncFileID, name='longitudes', &
                       xtype=NF90_DOUBLE, dimids = (/ nodes_3DimID /), varid=VarID)
      call nc_check(io,'nc_write_model_atts', 'longitudes def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,VarID,'units','degrees East'),&
                    'nc_write_model_atts', 'longitude units  '//trim(filename))

      io = nf90_def_var(ncid=ncFileID, name='latitudes', &
                       xtype=NF90_DOUBLE, dimids = (/ nodes_3DimID /), varid=VarID)
      call nc_check(io,'nc_write_model_atts', 'latitudes def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,VarID,'units','degrees North'),&
                    'nc_write_model_atts', 'latitudes units  '//trim(filename))

      io = nf90_def_var(ncid=ncFileID, name='depths', &
                       xtype=NF90_REAL, dimids = (/ nodes_3DimID /), varid=VarID)
      call nc_check(io,'nc_write_model_atts', 'depths def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,VarID,'units','meters'),&
                    'nc_write_model_atts', 'depths units  '//trim(filename))

      io = nf90_def_var(ncid=ncFileID, name='node_table', &
                xtype=NF90_INT, dimids = (/ nVertLevelsDimID, nodes_2DimID /), varid=VarID)
      call nc_check(io,'nc_write_model_atts', 'node_table def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,VarID,'long_name','nod3D_below_nod2D'),&
                    'nc_write_model_atts', 'node_table long_name  '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,VarID,'valid_range',(/ 1, nVertices /)),&
                    'nc_write_model_atts', 'node_table long_name  '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,VarID,'description', &
                   'table defining packing order. Given a level and a &
              &horizontal cell, return the vertex index (between 1 and nodes_3d). &
              &A value outside the valid_range means there is no wet location.'),&
                    'nc_write_model_atts', 'node_table description  '//trim(filename))

   endif

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

   if (diagnostic_metadata) then

      call nc_check(NF90_inq_varid(ncFileID, 'longitudes', VarID), &
                    'nc_write_model_atts', 'longitudes inq_varid '//trim(filename))
      call nc_check(nf90_put_var(ncFileID, VarID, coord_nod3D(1,:) ), &
                   'nc_write_model_atts', 'longitudes put_var '//trim(filename))

      call nc_check(NF90_inq_varid(ncFileID, 'latitudes', VarID), &
                    'nc_write_model_atts', 'latitudes inq_varid '//trim(filename))
      call nc_check(nf90_put_var(ncFileID, VarID, coord_nod3D(2,:) ), &
                   'nc_write_model_atts', 'latitudes put_var '//trim(filename))

      call nc_check(NF90_inq_varid(ncFileID, 'depths', VarID), &
                    'nc_write_model_atts', 'depths inq_varid '//trim(filename))
      call nc_check(nf90_put_var(ncFileID, VarID, coord_nod3D(3,:) ), &
                   'nc_write_model_atts', 'depths put_var '//trim(filename))

      call nc_check(NF90_inq_varid(ncFileID, 'node_table', VarID), &
                    'nc_write_model_atts', 'node_table inq_varid '//trim(filename))
      call nc_check(nf90_put_var(ncFileID, VarID, nod3D_below_nod2D ), &
                   'nc_write_model_atts', 'node_table put_var '//trim(filename))

   endif

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts


!------------------------------------------------------------------
!>

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

      if (do_output() .and. debug > 9) then
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
!>

function get_model_size()

! Returns the size of the model as an integer.
! Required for all applications.

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size


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

subroutine ens_mean_for_model(filter_ens_mean)

! If needed by the model interface, this is the current mean
! for all state vector items across all ensembles.

real(r8), intent(in) :: filter_ens_mean(:)

if ( .not. module_initialized ) call static_init_model

! Not needed by FESOM ... just a stub.

end subroutine ens_mean_for_model


!------------------------------------------------------------------
!>

subroutine end_model()

! Does any shutdown and clean-up needed for model.

if (allocated(cell_locations))  deallocate(cell_locations)
if (allocated(cell_kinds))      deallocate(cell_kinds)
if (allocated(close_cell_inds)) deallocate(close_cell_inds)
if (allocated(depths))          deallocate(depths)

if (close_structure_allocated) call finalize_closest_center()

end subroutine end_model


!------------------------------------------------------------------
!>

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
!>

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, &
                         obs_loc, obs_kind, num_close, close_ind, dist)
!
! Given a DART location (referred to as "base") and a set of candidate
! locations & kinds (obs, obs_kind), returns the subset close to the
! "base", their indices, and their distances to the "base" ...

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
real(r8), optional,  dimension(:), intent(out)   :: dist

! If you want to impose some sort of special localization, you can key
! off things like the obs_kind and make things 'infinitely' far away.
! Otherwise, this does nothing. Take a look at the POP model_mod.f90 for an example.

call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                          num_close, close_ind, dist)

end subroutine get_close_obs


!------------------------------------------------------------------
!>

subroutine init_time(time)

! Companion interface to init_conditions. Returns a time that is somehow
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

! this shuts up the compiler warnings about unused variables
time = set_time(0, 0)

write(string1,*) 'Cannot initialize FESOM time via subroutine call.'
write(string2,*) 'input.nml:start_from_restart cannot be FALSE'
call error_handler(E_ERR, 'init_time', string1, &
           source, revision, revdate, text2=string2)

end subroutine init_time


!------------------------------------------------------------------
!>

subroutine init_conditions(x)

! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'Cannot initialize FESOM time via subroutine call.'
write(string2,*) 'input.nml:start_from_restart cannot be FALSE'
call error_handler(E_ERR, 'init_conditions', string1, &
           source, revision, revdate, text2=string2)

! this shuts up the compiler warnings about unused variables
x = 0.0_r8

end subroutine init_conditions


!------------------------------------------------------------------
!>

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

write(string1,*) 'Cannot advance FESOM with a subroutine call; async cannot equal 0'
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
!>

subroutine analysis_file_to_statevector(filename, state_vector, model_time)

! Reads the current time and state variables from a FESOM analysis
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
if (do_output()) &
    call print_time(model_time,'time in restart file '//trim(filename))
if (do_output()) &
    call print_date(model_time,'date in restart file '//trim(filename))

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
   string2 = trim(filename)//' '//trim(varname)

   ! determine the shape of the netCDF variable

   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            'analysis_file_to_statevector', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
            'analysis_file_to_statevector', 'inquire '//trim(string2))

   mystart = 1   ! These are arrays, actually.
   mycount = 1

   ! Only checking the shape of the variable - sans TIME
   DimCheck : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'analysis_file_to_statevector', string1)

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string2),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         call error_handler(E_ERR,'analysis_file_to_statevector',string1,source,revision,revdate)
      endif

      mycount(i) = dimlen

   enddo DimCheck

   where(dimIDs == TimeDimID) mystart = TimeDimLength  ! pick the latest time
   where(dimIDs == TimeDimID) mycount = 1              ! only use one time

   if (debug > 0) then
      write(string1,*)'..  '//trim(varname)//' start = ',mystart(1:ncNdims)
      write(string3,*)        trim(varname)//' count = ',mycount(1:ncNdims)
      call error_handler(E_MSG,'analysis_file_to_statevector',string1,text2=string3)
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
!>

subroutine statevector_to_analysis_file(state_vector, filename, statetime)

! Writes the posterior state from a dart state
! vector (1d array) into a FESOM file.

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
! of the FESOM analysis file, and state vector contents from a different
! time won't be consistent with the rest of the file.

model_time = get_analysis_time(ncFileID, filename)

if ( model_time /= statetime ) then
   call print_time( statetime,'DART current time',logfileunit)
   call print_time(model_time,'FESOM current time',logfileunit)
   call print_time( statetime,'DART current time')
   call print_time(model_time,'FESOM current time')
   write(string1,*)trim(filename),' current time must be equal to model time'
   call error_handler(E_ERR,'statevector_to_analysis_file',string1,source,revision,revdate)
endif

if (do_output()) call print_time(statetime,'time in DART file '//trim(filename))
if (do_output()) call print_date(statetime,'date in DART file '//trim(filename))

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

   if ( .not. progvar(ivar)%replace ) cycle PROGVARLOOP

   varname = trim(progvar(ivar)%varname)
   string2 = trim(filename)//' '//trim(varname)

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

   if (debug > 2) then
      write(string1,*)'..  '//trim(varname)//' start is ',mystart(1:ncNdims)
      write(string3,*)        trim(varname)//' count is ',mycount(1:ncNdims)
      call error_handler(E_MSG,'statevector_to_analysis_file',string1,text2=string3)
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
!>

subroutine do_clamping(out_of_range_fail, range, dimsize, varname, &
                       array_1d, array_2d, array_3d)
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
   if ( range(1) /= MISSING_R8 ) then

      if ( out_of_range_fail ) then
         if ( minval(array_1d) < range(1) ) then
            write(string1, *) 'min data val = ', minval(array_1d), &
                              'min data bounds = ', range(1)
            call error_handler(E_ERR, 'do_clamping', &
                        'Variable '//trim(varname)//' failed lower bounds check.', &
                         source,revision,revdate)
         endif
      else
         where ( array_1d < range(1) ) array_1d = range(1)
      endif

   endif ! min range set

   ! is upper bound set
   if ( range(2) /= MISSING_R8 ) then

      if ( out_of_range_fail ) then
         if ( maxval(array_1d) > range(2) ) then
            write(string1, *) 'max data val = ', maxval(array_1d), &
                              'max data bounds = ', range(2)
            call error_handler(E_ERR, 'do_clamping', &
                        'Variable '//trim(varname)//' failed upper bounds check.', &
                         source,revision,revdate, text2=string1)
         endif
      else
         where ( array_1d > range(2) ) array_1d = range(2)
      endif

   endif ! max range set

   write(string1, '(A,A32,2F16.7)') 'BOUND min/max ', trim(varname), &
                      minval(array_1d), maxval(array_1d)
   call error_handler(E_MSG,'do_clamping',string1,source,revision,revdate)

else if (dimsize == 2) then
   if (.not. present(array_2d)) then
      call error_handler(E_ERR, 'do_clamping', 'Internal error.  Should not happen', &
                         source,revision,revdate, text2='array_2d not present for 2d case')
   endif

   ! is lower bound set
   if ( range(1) /= MISSING_R8 ) then

      if ( out_of_range_fail ) then
         if ( minval(array_2d) < range(1) ) then
            write(string1, *) 'min data val = ', minval(array_2d), &
                              'min data bounds = ', range(1)
            call error_handler(E_ERR, 'do_clamping', &
                        'Variable '//trim(varname)//' failed lower bounds check.', &
                         source,revision,revdate)
         endif
      else
         where ( array_2d < range(1) ) array_2d = range(1)
      endif

   endif ! min range set

   ! is upper bound set
   if ( range(2) /= MISSING_R8 ) then

      if ( out_of_range_fail ) then
         if ( maxval(array_2d) > range(2) ) then
            write(string1, *) 'max data val = ', maxval(array_2d), &
                              'max data bounds = ', range(2)
            call error_handler(E_ERR, 'do_clamping', &
                        'Variable '//trim(varname)//' failed upper bounds check.', &
                         source,revision,revdate, text2=string1)
         endif
      else
         where ( array_2d > range(2) ) array_2d = range(2)
      endif

   endif ! max range set

   write(string1, '(A,A32,2F16.7)') 'BOUND min/max ', trim(varname), &
                      minval(array_2d), maxval(array_2d)
   call error_handler(E_MSG,'do_clamping',string1,source,revision,revdate)

else if (dimsize == 3) then
   if (.not. present(array_3d)) then
      call error_handler(E_ERR, 'do_clamping', 'Internal error.  Should not happen', &
                         source,revision,revdate, text2='array_3d not present for 3d case')
   endif

   ! is lower bound set
   if ( range(1) /= MISSING_R8 ) then

      if ( out_of_range_fail ) then
         if ( minval(array_3d) < range(1) ) then
            write(string1, *) 'min data val = ', minval(array_3d), &
                              'min data bounds = ', range(1)
            call error_handler(E_ERR, 'do_clamping', &
                        'Variable '//trim(varname)//' failed lower bounds check.', &
                         source,revision,revdate)
         endif
      else
         where ( array_3d < range(1) ) array_3d = range(1)
      endif

   endif ! min range set

   ! is upper bound set
   if ( range(2) /= MISSING_R8 ) then

      if ( out_of_range_fail ) then
         if ( maxval(array_3d) > range(2) ) then
            write(string1, *) 'max data val = ', maxval(array_3d), &
                              'max data bounds = ', range(2)
            call error_handler(E_ERR, 'do_clamping', &
                        'Variable '//trim(varname)//' failed upper bounds check.', &
                         source,revision,revdate, text2=string1)
         endif
      else
         where ( array_3d > range(2) ) array_3d = range(2)
      endif

   endif ! max range set

   write(string1, '(A,A32,2F16.7)') 'BOUND min/max ', trim(varname), &
                      minval(array_3d), maxval(array_3d)
   call error_handler(E_MSG,'do_clamping',string1,source,revision,revdate)

else
   write(string1, *) 'dimsize of ', dimsize, ' found where only 1-3 expected'
   call error_handler(E_MSG,'do_clamping','Internal error, should not happen', &
                      source,revision,revdate, text2=string1)
endif   ! dimsize

end subroutine do_clamping


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


!------------------------------------------------------------------
!>

subroutine write_model_time(time_filename, model_time, adv_to_time)
 character(len=*), intent(in)           :: time_filename
 type(time_type),  intent(in)           :: model_time
 type(time_type),  intent(in), optional :: adv_to_time

integer :: iunit
type(time_type) :: deltatime
type(time_type) :: new_model_time, assim_start, assim_end
integer :: seconds, days

iunit = open_file(time_filename, form='formatted', action='write')

deltatime = set_time(assimilation_period_seconds, assimilation_period_days)

new_model_time = model_time + deltatime
assim_start    = new_model_time - deltatime/2 + set_time(1,0)
assim_end      = new_model_time + deltatime/2

! By writing the days,seconds to strings with a free-format write,
! you can avoid any decision about formatting precision.
! The '(A)' syntax avoids printing the quotes delimiting the string

call get_time(new_model_time, seconds, days)
write(string1,*) days
write(string2,*) seconds
write(iunit,'(A,A)') ' init_time_days     = ',trim(string1)
write(iunit,'(A,A)') ' init_time_seconds  = ',trim(string2)

call get_time(assim_start, seconds, days)
write(string1,*) days
write(string2,*) seconds
write(iunit,'(A,A)') ' first_obs_days     = ',trim(string1)
write(iunit,'(A,A)') ' first_obs_seconds  = ',trim(string2)

call get_time(assim_end, seconds, days)
write(string1,*) days
write(string2,*) seconds
write(iunit,'(A,A)') ' last_obs_days      = ',trim(string1)
write(iunit,'(A,A)') ' last_obs_seconds   = ',trim(string2)

string2 = time_to_string(new_model_time)

call error_handler(E_MSG,'write_model_time:','next model time should be '//trim(string2))

write(iunit, '(A)') trim(string2)

call print_time(new_model_time,'FESOM    stop at :',  iunit)
call print_time(    model_time,'FESOM current at :',  iunit)
call print_date(new_model_time,'stop    date :',  iunit)
call print_date(    model_time,'current date :',  iunit)

if (present(adv_to_time)) then
   string1 = time_to_string(adv_to_time)
   write(iunit, '(A)') trim(string1)

   deltatime = adv_to_time - model_time
   string1 = time_to_string(deltatime, interval=.true.)
   write(iunit, '(A)') trim(string1)
endif

call close_file(iunit)

end subroutine write_model_time


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

call get_close_maxdist_init(cc_gc, maxdist = maxdist_km*PI/20000.0_r8)
call get_close_obs_init(cc_gc, nCells, cell_locations)

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
!>

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
!>

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
!>

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
!>

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
!>

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
!>

subroutine parse_variable_input( state_variables, ncid, filename, ngood, table )

character(len=*), intent(in)  :: state_variables(:)
integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: filename
integer,          intent(out) :: ngood
character(len=*), intent(out) :: table(:,:)

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: dimname
integer :: nrows, ncols, i, j, VarID, dimlen, numdims
logical :: failure

character(len=NF90_MAX_NAME) :: varname       ! column 1
character(len=NF90_MAX_NAME) :: dartstr       ! column 2
character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
character(len=NF90_MAX_NAME) :: state_or_aux  ! column 5

if ( .not. module_initialized ) call static_init_model

failure = .FALSE. ! perhaps all with go well

nrows = size(table,1)
ncols = size(table,2)

ngood = 0
MyLoop : do i = 1, nrows

   varname      = trim(state_variables(NUM_STATE_TABLE_COLUMNS*i-4))
   dartstr      = trim(state_variables(NUM_STATE_TABLE_COLUMNS*i-3))
   minvalstring = trim(state_variables(NUM_STATE_TABLE_COLUMNS*i-2))
   maxvalstring = trim(state_variables(NUM_STATE_TABLE_COLUMNS*i-1))
   state_or_aux = trim(state_variables(NUM_STATE_TABLE_COLUMNS*i  ))

   table(i,VARNAME_INDEX) = trim(varname)
   table(i,   KIND_INDEX) = trim(dartstr)
   table(i, MINVAL_INDEX) = trim(minvalstring)
   table(i, MAXVAL_INDEX) = trim(maxvalstring)
   table(i,REPLACE_INDEX) = trim(state_or_aux)

   if ( table(i,1) == ' ' .and. table(i,2) == ' ' ) exit MyLoop ! Found end of list.

   if ( any(table(i,:) == ' ') ) then
      string1 = '...  model_nml:"variables" not fully specified'
      write(string2,*)'failing on line ',i
      call error_handler(E_ERR, 'parse_variable_input', string1, &
                 source, revision, revdate, text2=string2)
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

   ! Make sure DART kind is valid

   if( get_raw_obs_kind_index(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'parse_variable_input',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector

   if (debug > 0) then
      write(string1,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
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
!%!    integer :: varsize       ! prod(dimlens(1:numdims))
!%!    integer :: index1        ! location in dart state vector of first occurrence
!%!    integer :: indexN        ! location in dart state vector of last  occurrence
!%!    integer :: dart_kind
!%!    character(len=paramname_length) :: kind_string
!%!    logical  :: clamping     ! does variable need to be range-restricted before
!%!    real(r8) :: range(2)     ! being stuffed back into FESOM analysis file.
!%!    logical  :: replace      ! does it need to be stuffed back into FESOM file.
!%! end type progvartype

integer :: i

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


!------------------------------------------------------------------
!>

!#! !>@todo FIXME : threed_cartesian/location assumes everything is in height
!#!
!#! subroutine vert_convert(x, location, obs_kind, ztypeout, istatus)
!#!
!#! ! This subroutine converts a given ob/state vertical coordinate to
!#! ! the vertical localization coordinate type requested through the
!#! ! model_mod namelist.
!#! !
!#! ! Notes: (1) obs_kind is only necessary to check whether the ob
!#! !            is an identity ob.
!#! !        (2) This subroutine can convert both obs' and state points'
!#! !            vertical coordinates. Remember that state points get
!#! !            their DART location information from get_state_meta_data
!#! !            which is called by filter_assim during the assimilation
!#! !            process.
!#! !        (3) x is the relevant DART state vector for carrying out
!#! !            computations necessary for the vertical coordinate
!#! !            transformations. As the vertical coordinate is only used
!#! !            in distance computations, this is actually the "expected"
!#! !            vertical coordinate, so that computed distance is the
!#! !            "expected" distance. Thus, under normal circumstances,
!#! !            x that is supplied to vert_convert should be the
!#! !            ensemble mean. Nevertheless, the subroutine has the
!#! !            functionality to operate on any DART state vector that
!#! !            is supplied to it.
!#!
!#! real(r8), dimension(:), intent(in)    :: x
!#! type(location_type),    intent(inout) :: location
!#! integer,                intent(in)    :: obs_kind
!#! integer,                intent(in)    :: ztypeout
!#! integer,                intent(out)   :: istatus
!#!
!#!  ! zin and zout are the vert values coming in and going out.
!#!  ! ztype{in,out} are the vert types as defined by the 3d sphere
!#!  ! locations mod (location/threed_sphere/location_mod.f90)
!#!  real(r8) :: llv_loc(3)
!#!  real(r8) :: zin, zout, tk, fullp, surfp
!#!  real(r8) :: weights(3), zk_mid(3), values(3), fract(3), fdata(3)
!#!  integer  :: ztypein, i
!#!  integer  :: k_low(3), k_up(3), c(3), n
!#!  integer  :: ivars(3)
!#!  type(location_type) :: surfloc
!#! !$!
!#! !$! ! assume failure.
!#! !$! istatus = 1
!#! !$!
!#! !$! ! initialization
!#! !$! k_low   = 0.0_r8
!#! !$! k_up    = 0.0_r8
!#! !$! weights = 0.0_r8
!#! !$!
!#! !$! ! first off, check if ob is identity ob.  if so get_state_meta_data() will
!#! !$! ! have returned location information already in the requested vertical type.
!#! !$! if (obs_kind < 0) then
!#! !$!    call get_state_meta_data(obs_kind,location)
!#! !$!    istatus = 0
!#! !$!    return
!#! !$! endif
!#! !$!
!#! !$! !#! ! if the existing coord is already in the requested vertical units
!#! !$! !#! ! or if the vert is 'undef' which means no specifically defined
!#! !$! !#! ! vertical coordinate, return now.
!#! !$!
!#! !>@todo FIXME : threed_cartesian/location assumes everything is in height
!#! !$!
!#! !$! !#! ztypein  = nint(query_location(location, 'which_vert'))
!#! !$! !#! if ((ztypein == ztypeout) .or. (ztypein == VERTISUNDEF)) then
!#! !$! !#!    istatus = 0
!#! !$! !#!    return
!#! !$! !#! else
!#! !$! !#!    if (do_output() .and. debug > 9) then
!#! !$! !#!       write(string1,'(A,3X,2I3)') 'ztypein, ztypeout:',ztypein,ztypeout
!#! !$! !#!       call error_handler(E_MSG, 'vert_convert',string1,source, revision, revdate)
!#! !$! !#!    endif
!#! !$! !#! endif
!#! !$!
!#! !$! ! we do need to convert the vertical.  start by
!#! !$! ! extracting the location lon/lat/vert values.
!#! !$! llv_loc = get_location(location)
!#! !$!
!#! !$! ! the routines below will use zin as the incoming vertical value
!#! !$! ! and zout as the new outgoing one.  start out assuming failure
!#! !$! ! (zout = missing) and wait to be pleasantly surprised when it works.
!#! !$! zin     = llv_loc(3)
!#! !$! !@>todo FIXME : this probably does not matter since everything is always in
!#! !$! !> hight, some models have to do vertival conversion for incomming observations
!#! !$! zout    = MISSING_R8
!#! !$!
!#! !$! ! if the vertical is missing to start with, return it the same way
!#! !$! ! with the requested type as out.
!#! !$! if (zin == MISSING_R8) then
!#! !$!    location = set_location(llv_loc(1),llv_loc(2),MISSING_R8,ztypeout)
!#! !$!    return
!#! !$! endif
!#! !$!
!#! !$! !>@todo FIXME : threed_cartesian/location assumes everything is in height
!#! !$!
!#! !$! !#! ! Convert the incoming vertical type (ztypein) into the vertical
!#! !$! !#! ! localization coordinate given in the namelist (ztypeout).
!#! !$! !#! ! Various incoming vertical types (ztypein) are taken care of
!#! !$! !#! ! inside find_vert_level. So we only check ztypeout here.
!#! !$! !#!
!#! !$! !#! ! convert into:
!#! !$! !#! select case (ztypeout)
!#! !$! !#!
!#! !$! !#!    ! ------------------------------------------------------------
!#! !$! !#!    ! outgoing vertical coordinate should be 'model level number'
!#! !$! !#!    ! ------------------------------------------------------------
!#! !$! !#!    case (VERTISLEVEL)
!#! !$! !#!
!#! !$! !#!    ! Identify the three cell ids (c) in the triangle enclosing the obs and
!#! !$! !#!    ! the vertical indices for the triangle at two adjacent levels (k_low and k_up)
!#! !$! !#!    ! and the fraction (fract) for vertical interpolation.
!#! !$! !#!
!#! !$! !#!   call find_triangle_vert_indices (x, location, n, c, k_low, k_up, fract, weights, istatus)
!#! !$! !#!   if(istatus /= 0) return
!#! !$! !#!
!#! !$! !#!   zk_mid = k_low + fract
!#! !$! !#!   zout = sum(weights * zk_mid)
!#! !$! !#!
!#! !$! !#!   if (do_output() .and. debug > 9) then
!#! !$! !#!      write(string2,'("Zk:",3F8.2," => ",F8.2)') zk_mid,zout
!#! !$! !#!      call error_handler(E_MSG, 'vert_convert',string2,source, revision, revdate)
!#! !$! !#!   endif
!#! !$! !#!
!#! !$! !#!    ! ------------------------------------------------------------
!#! !$! !#!    ! outgoing vertical coordinate should be 'pressure' in Pa
!#! !$! !#!    ! ------------------------------------------------------------
!#! !$! !#!    case (VERTISPRESSURE)
!#! !$! !#!
!#! !$! !#!    ! Need to get base offsets for the potential temperature, density, and water
!#! !$! !#!    ! vapor mixing fields in the state vector
!#! !$! !#! ! TJH   ivars(1) = get_progvar_index_from_kind(KIND_POTENTIAL_TEMPERATURE)
!#! !$! !#! ! TJH   ivars(2) = get_progvar_index_from_kind(KIND_DENSITY)
!#! !$! !#! ! TJH   ivars(3) = get_progvar_index_from_kind(KIND_VAPOR_MIXING_RATIO)
!#! !$! !#!
!#! !$! !#! string1 = 'fix VERTISPRESSURE get base offsets - detritus from atmosphere'
!#! !$! !#! call error_handler(E_ERR,'vert_convert',string1,source,revision,revdate)
!#! !$! !#!
!#! !$! !#!    if (any(ivars(1:3) < 0)) then
!#! !$! !#!       write(string1,*) 'Internal error, cannot find one or more of: theta, rho, qv'
!#! !$! !#!       call error_handler(E_ERR, 'vert_convert',string1,source, revision, revdate)
!#! !$! !#!    endif
!#! !$! !#!
!#! !$! !#!    ! Get theta, rho, qv at the interpolated location
!#! !$! !#!    call compute_scalar_with_barycentric (x, location, 3, ivars, values, istatus)
!#! !$! !#!    if (istatus /= 0) return
!#! !$! !#!
!#! !$! !#!    ! Convert theta, rho, qv into pressure
!#! !$! !#!    call compute_full_pressure(values(1), values(2), values(3), zout, tk)
!#! !$! !#!    if (do_output() .and. debug > 9) then
!#! !$! !#!       write(string2,'("zout_in_pressure, theta, rho, qv:",3F10.2,F15.10)') zout, values
!#! !$! !#!       call error_handler(E_MSG, 'vert_convert',string2,source, revision, revdate)
!#! !$! !#!    endif
!#! !$! !#!
!#! !$! !#!    ! ------------------------------------------------------------
!#! !$! !#!    ! outgoing vertical coordinate should be 'depth' in meters
!#! !$! !#!    ! ------------------------------------------------------------
!#! !$! !#!    case (VERTISHEIGHT)
!#!
!#! !>@todo FIXME : might still want something similar to this code
!#!
!#!      call find_triangle_vert_indices (x, location, n, c, k_low, k_up, fract, weights, istatus)
!#!      if (istatus /= 0) return
!#!
!#!      ! now have vertically interpolated values at cell centers.
!#!      ! use horizontal weights to compute value at interp point.
!#!      zout = sum(weights * fdata)
!#!
!#! !$! !#!
!#! !$! !#!    ! ------------------------------------------------------------
!#! !$! !#!    ! outgoing vertical coordinate should be 'scale height' (a ratio)
!#! !$! !#!    ! ------------------------------------------------------------
!#! !$! !#!    case (VERTISSCALEHEIGHT)
!#! !$! !#!
!#! !$! !#!    ! Scale Height is defined here as: -log(pressure / surface_pressure)
!#! !$! !#!
!#! !$! !#!    ! Need to get base offsets for the potential temperature, density, and water
!#! !$! !#!    ! vapor mixing fields in the state vector
!#! !$! !#! ! TJH   ivars(1) = get_progvar_index_from_kind(KIND_POTENTIAL_TEMPERATURE)
!#! !$! !#! ! TJH   ivars(2) = get_progvar_index_from_kind(KIND_DENSITY)
!#! !$! !#! ! TJH   ivars(3) = get_progvar_index_from_kind(KIND_VAPOR_MIXING_RATIO)
!#! !$! !#!
!#! !$! !#! string1 = 'fix vertisscaleheight get base offsets - detritus from atmosphere'
!#! !$! !#! call error_handler(E_ERR,'vert_convert',string1,source,revision,revdate)
!#! !$! !#!
!#! !$! !#!    ! Get theta, rho, qv at the interpolated location
!#! !$! !#!    call compute_scalar_with_barycentric (x, location, 3, ivars, values, istatus)
!#! !$! !#!    if (istatus /= 0) return
!#! !$! !#!
!#! !$! !#!    ! Convert theta, rho, qv into pressure
!#! !$! !#!    call compute_full_pressure(values(1), values(2), values(3), fullp, tk)
!#! !$! !#!    if (do_output() .and. debug > 9) then
!#! !$! !#!       write(string2,'("zout_full_pressure, theta, rho, qv:",3F10.2,F15.10)') fullp, values
!#! !$! !#!       call error_handler(E_MSG, 'vert_convert',string2,source, revision, revdate)
!#! !$! !#!    endif
!#! !$! !#!
!#! !$! !#!    ! Get theta, rho, qv at the surface corresponding to the interpolated location
!#! !$! !#!    surfloc = set_location(llv_loc(1), llv_loc(2), 1.0_r8, VERTISLEVEL)
!#! !$! !#!    call compute_scalar_with_barycentric (x, surfloc, 3, ivars, values, istatus)
!#! !$! !#!    if (istatus /= 0) return
!#! !$! !#!
!#! !$! !#!    ! Convert surface theta, rho, qv into pressure
!#! !$! !#!    call compute_full_pressure(values(1), values(2), values(3), surfp, tk)
!#! !$! !#!    if (do_output() .and. debug > 9) then
!#! !$! !#!       write(string2,'("zout_surf_pressure, theta, rho, qv:",3F10.2,F15.10)') surfp, values
!#! !$! !#!       call error_handler(E_MSG, 'vert_convert',string2,source, revision, revdate)
!#! !$! !#!    endif
!#! !$! !#!
!#! !$! !#!    ! and finally, convert into scale height
!#! !$! !#!    if (surfp /= 0.0_r8) then
!#! !$! !#!       zout = -log(fullp / surfp)
!#! !$! !#!    else
!#! !$! !#!       zout = MISSING_R8
!#! !$! !#!    endif
!#! !$! !#!
!#! !$! !#!    if (do_output() .and. debug > 9) then
!#! !$! !#!       write(string2,'("zout_in_pressure:",F10.2)') zout
!#! !$! !#!       call error_handler(E_MSG, 'vert_convert',string2,source, revision, revdate)
!#! !$! !#!    endif
!#! !$! !#!
!#! !$! !#!    ! -------------------------------------------------------
!#! !$! !#!    ! outgoing vertical coordinate is unrecognized
!#! !$! !#!    ! -------------------------------------------------------
!#! !$! !#!    case default
!#! !$! !#!       write(string1,*) 'Requested vertical coordinate not recognized: ', ztypeout
!#! !$! !#!       call error_handler(E_ERR,'vert_convert', string1, &
!#! !$! !#!                          source, revision, revdate)
!#! !$! !#!
!#! !$! !#! end select   ! outgoing vert type
!#! !$!
!#! !$! ! Returned location
!#! !$! location = set_location(llv_loc(1),llv_loc(2),zout,ztypeout)
!#! !$!
!#! !$! ! Set successful return code only if zout has good value
!#! !$! if(zout /= MISSING_R8) istatus = 0
!#! !$!
!#! !$!
!#! end subroutine vert_convert

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

integer :: num_close, surface_index, iclose, indx
real(r8) :: closest
real(r8) :: distance

find_closest_surface_location = -1

! Generate the list of indices into the DART vector of the close candidates.

call loc_get_close_obs(cc_gc, location, obs_kind, cell_locations, cell_kinds, &
                       num_close, close_ind=close_cell_inds)

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
!>

subroutine finalize_closest_center()

! get rid of storage associated with get close lookup table

!>@ TODO make this initialization a global T/F variable and only call if ...

call get_close_obs_destroy(cc_gc)

end subroutine finalize_closest_center


!$! !------------------------------------------------------------
!$!
!$! function on_boundary(cellid)
!$!
!$! ! use the surface (level 1) to determine if any edges (or vertices?)
!$! ! are on the boundary, and return true if so.   if the global flag
!$! ! is set, skip all code and return false immediately.
!$!
!$! integer,  intent(in)  :: cellid
!$! logical               :: on_boundary
!$!
!$! integer :: vertical
!$!
!$! vertical = 1
!$!
!$! if (global_grid) then
!$!    on_boundary = .false.
!$!    return
!$! endif
!$!
!$! write(*,*) 'boundaryCell ', cellid, boundaryCell(vertical,cellid)
!$!
!$! on_boundary = boundaryCell(vertical,cellid) .eq. 1
!$!
!$! end function on_boundary
!$!
!$! !------------------------------------------------------------
!$!
!$! function inside_cell(cellid, lat, lon)
!$!
!$! ! this function no longer really determines if we are inside
!$! ! the cell or not.  what it does do is determine if the nearest
!$! ! cell is on the grid boundary in any way and says no if it is
!$! ! a boundary.  if we have a flag saying this a global grid, we
!$! ! can avoid doing any work and immediately return true.  for a
!$! ! global atmosphere this is always so; for a regional atmosphere
!$! ! and for the ocean (which does not have cells on land) this is
!$! ! necessary test.
!$!
!$! integer,  intent(in)  :: cellid
!$! real(r8), intent(in)  :: lat, lon
!$! logical               :: inside_cell
!$!
!$! ! do this completely with topology of the grid.  if any of
!$! ! the cell edges are marked as boundary edges, return no.
!$! ! otherwise return yes.
!$!
!$! integer :: nedges, i, edgeid, vert
!$!
!$! ! if we're on a global grid, skip all this code
!$! if (global_grid) then
!$!    inside_cell = .true.
!$!    return
!$! endif
!$!
!$! nedges = nEdgesOnCell(cellid)
!$!
!$! ! go around the edges and check the boundary array.
!$! ! if any are true, return false.  even if we are inside
!$! ! this cell, we aren't going to be able to interpolate it
!$! ! so shorten the code path.
!$!
!$! ! FIXME: at some point we can be more selective and try to
!$! ! interpolate iff the edges of the three cells which are
!$! ! going to contribute edges to the RBF exist, even if some
!$! ! of the other cell edges are on the boundary.  so this
!$! ! decision means we won't be interpolating some obs that in
!$! ! theory we have enough information to interpolate.  but it
!$! ! is conservative for now - we certainly won't try to interpolate
!$! ! outside the existing grid.
!$!
!$! ! FIXME: can this loop over edges be replaced by a single test of the 'boundaryCell' array?
!$! do i=1, nedges
!$!    edgeid = edgesOnCell(i, cellid)
!$!
!$!    ! FIXME: this is an int array.  is it 0=false,1=true?
!$!    ! BOTHER - we need the vert for this and we don't have it
!$!    ! and in fact can't compute it if the interpolation point
!$!    ! has pressure or depth as its vertical coordinate.
!$!    vert = 1
!$!
!$!    if (boundaryEdge(vert, edgeid) > 0) then
!$!       inside_cell = .false.
!$!       return
!$!    endif
!$!
!$! enddo
!$!
!$! inside_cell = .true.
!$!
!$! end function inside_cell

!$-----------------------------------------------------------

!$! function closest_vertex_ll(cellid, lat, lon)
!$!
!$! ! Return the vertex id of the closest one to the given point
!$! ! this version uses lat/lon.  see closest_vertex_xyz for the
!$! ! cartesian version.
!$!
!$! integer,  intent(in)  :: cellid
!$! real(r8), intent(in)  :: lat, lon
!$! integer               :: closest_vertex_ll
!$!
!$! real(r8) :: px, py, pz
!$!

!>@todo FIXME : should not need to convert between latlon and xyz
!$! ! use the same radius as FESOM for computing this
!$! call latlon_to_xyz(lat, lon, px, py, pz)
!$!
!$! closest_vertex_ll = closest_vertex_xyz(cellid, px, py, pz)
!$! if ((closest_vertex_ll < 0)  .and. &
!$!     (debug > 8) .and. do_output()) &
!$!    print *, 'cannot find nearest vertex to lon, lat: ', lon, lat, &
!$!             'cellid', cellid,'px: ', px,'py: ', py,'pz: ', pz
!$!
!$! end function closest_vertex_ll
!$!
!$! !------------------------------------------------------------
!$!
!$! function closest_vertex_xyz(cellid, px, py, pz)
!$!
!$! ! Return the vertex id of the closest one to the given point
!$! ! see closest_vertex_ll for the lat/lon version (which calls this)
!$!
!$! integer,  intent(in)  :: cellid
!$! real(r8), intent(in)  :: px, py, pz
!$! integer               :: closest_vertex_xyz
!$!
!$! integer :: nverts, i, vertexid
!$! real(r8) :: distsq, closest_dist, dx, dy, dz
!$!
!$! ! nedges and nverts is same in a closed figure
!$! nverts = nEdgesOnCell(cellid)
!$!
!$! closest_dist = 1.0e38_r8   ! something really big; these are meters not radians
!$! closest_vertex_xyz = -1
!$!
!$! do i=1, nverts
!$!    vertexid = verticesOnCell(i, cellid)
!$!    distsq = (dx * dx) + (dy * dy) + (dz * dz)
!$!    if (distsq < closest_dist) then
!$!       closest_dist = distsq
!$!       closest_vertex_xyz = vertexid
!$!    endif
!$! enddo
!$!
!$! end function closest_vertex_xyz

!------------------------------------------------------------
!----DON'T CHANGE THE REST-----------------------------------
!------------------------------------------------------------
!------------------------------------------------------------

!#! subroutine latlon_to_xyz(lat, lon, x, y, z)
!#!
!#! ! Given a lat, lon in degrees, return the cartesian x,y,z coordinate
!#! ! on the surface of a specified radius relative to the origin
!#! ! at the center of the earth.  (this radius matches the one
!#! ! used at FESOM grid generation time and must agree in order
!#! ! to be consistent with the cartisian coordinate arrays in
!#! ! the FESOM data files.)
!#!
!#! real(r8), intent(in)  :: lat, lon
!#! real(r8), intent(out) :: x, y, z
!#!
!#! real(r8) :: rlat, rlon
!#!
!#! rlat = lat * deg2rad
!#! rlon = lon * deg2rad
!#!
!#! x = radius * cos(rlon) * cos(rlat)
!#! y = radius * sin(rlon) * cos(rlat)
!#! z = radius * sin(rlat)
!#!
!#! end subroutine latlon_to_xyz
!#!
!#! !------------------------------------------------------------
!#!
!#! subroutine xyz_to_latlon(x, y, z, lat, lon)
!#!
!#! ! Given a cartesian x, y, z coordinate relative to the origin
!#! ! at the center of the earth, using a fixed radius specified
!#! ! by FESOM (in the grid generation step), return the corresponding
!#! ! lat, lon location in degrees.
!#!
!#! real(r8), intent(in)  :: x, y, z
!#! real(r8), intent(out) :: lat, lon
!#!
!#! real(r8) :: rlat, rlon
!#!
!#! ! right now this is only needed for debugging messages.
!#! ! the arc versions of routines are expensive.
!#!
!#! rlat = PI/2.0_r8 - acos(z/radius)
!#! rlon = atan2(y,x)
!#! if (rlon < 0) rlon = rlon + PI*2
!#!
!#! lat = rlat * rad2deg
!#! lon = rlon * rad2deg
!#!
!#! end subroutine xyz_to_latlon

!------------------------------------------------------------

!#! subroutine inside_triangle(t1, t2, t3, r, lat, lon, inside, weights)
!#!
!#! ! given 3 corners of a triangle and an xyz point, compute whether
!#! ! the point is inside the triangle.  this assumes r is coplanar
!#! ! with the triangle - the caller must have done the lat/lon to
!#! ! xyz conversion with a constant radius and then this will be
!#! ! true (enough).  sets inside to true/false, and returns the
!#! ! weights if true.  weights are set to 0 if false.
!#!
!#! real(r8), intent(in)  :: t1(3), t2(3), t3(3)
!#! real(r8), intent(in)  :: r(3), lat, lon
!#! logical,  intent(out) :: inside
!#! real(r8), intent(out) :: weights(3)
!#!
!#! ! check for degenerate cases first - is the test point located
!#! ! directly on one of the vertices?  (this case may be common
!#! ! if we're computing on grid point locations.)
!#! if (all(abs(r - t1) < roundoff)) then
!#!    inside = .true.
!#!    weights = (/ 1.0_r8, 0.0_r8, 0.0_r8 /)
!#!    return
!#! else if (all(abs(r - t2) < roundoff)) then
!#!    inside = .true.
!#!    weights = (/ 0.0_r8, 1.0_r8, 0.0_r8 /)
!#!    return
!#! else if (all(abs(r - t3) < roundoff)) then
!#!    inside = .true.
!#!    weights = (/ 0.0_r8, 0.0_r8, 1.0_r8 /)
!#!    return
!#! endif
!#!
!#! ! not a vertex. compute the weights.  if any are
!#! ! negative, the point is outside.  since these are
!#! ! real valued computations define a lower bound for
!#! ! numerical roundoff error and be sure that the
!#! ! weights are not just *slightly* negative.
!#! call get_3d_weights(r, t1, t2, t3, lat, lon, weights)
!#!
!#! if (any(weights < -roundoff)) then
!#!    inside = .false.
!#!    weights = 0.0_r8
!#!    return
!#! endif
!#!
!#! ! truncate barely negative values to 0
!#! inside = .true.
!#! where (weights < 0.0_r8) weights = 0.0_r8
!#! return
!#!
!#! end subroutine inside_triangle

!>@todo FIXME : function no used
!#! !------------------------------------------------------------
!#!
!#! function vector_magnitude(a)
!#!
!#! ! Given a cartesian vector, compute the magnitude
!#
!#! real(r8), intent(in)  :: a(3)
!#! real(r8) :: vector_magnitude
!#!
!#! vector_magnitude = sqrt(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))
!#!
!#! end function vector_magnitude
!#!
!#! !------------------------------------------------------------
!#!
!#! subroutine vector_cross_product(a, b, r)
!#!
!#! ! Given 2 cartesian vectors, compute the cross product of a x b
!#!
!#! real(r8), intent(in)  :: a(3), b(3)
!#! real(r8), intent(out) :: r(3)
!#!
!#! r(1) = a(2)*b(3) - a(3)*b(2)
!#! r(2) = a(3)*b(1) - a(1)*b(3)
!#! r(3) = a(1)*b(2) - a(2)*b(1)
!#!
!#! end subroutine vector_cross_product
!#!
!#! !------------------------------------------------------------
!#!
!#! function vector_dot_product(a, b)
!#!
!#! ! Given 2 cartesian vectors, compute the dot product of a . b
!#!
!#! real(r8), intent(in)  :: a(3), b(3)
!#! real(r8) :: vector_dot_product
!#!
!#! vector_dot_product = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
!#!
!#! end function vector_dot_product
!#!
!#! !------------------------------------------------------------
!#!
!#! subroutine vector_projection(a, b, r)
!#!
!#! ! Given 2 cartesian vectors, project a onto b
!#!
!#! real(r8), intent(in)  :: a(3), b(3)
!#! real(r8), intent(out) :: r(3)
!#!
!#! real(r8) :: ab_over_bb
!#!
!#! ab_over_bb = vector_dot_product(a, b) / vector_dot_product(b, b)
!#! r = (ab_over_bb) * b
!#!
!#! end subroutine vector_projection
!#!
!#! !------------------------------------------------------------
!#!
!#! subroutine determinant3(a, r)
!#!
!#! ! Given a 3x3 matrix, compute the determinant
!#!
!#! real(r8), intent(in)  :: a(3,3)
!#! real(r8), intent(out) :: r
!#!
!#! r = a(1,1)*(a(2,2)*a(3,3) - (a(3,2)*a(2,3))) + &
!#!     a(2,1)*(a(3,2)*a(1,3) - (a(3,3)*a(1,2))) + &
!#!     a(3,1)*(a(1,2)*a(2,3) - (a(2,2)*a(1,3)))
!#!
!#! end subroutine determinant3
!#!
!#! !------------------------------------------------------------
!#!
!#! subroutine invert3(a, r)
!#!
!#! ! Given a 3x3 matrix, compute the inverse
!#!
!#! real(r8), intent(in)  :: a(3,3)
!#! real(r8), intent(out) :: r(3,3)
!#!
!#! real(r8) :: det, b(3,3)
!#!
!#! call determinant3(a, det)
!#! if (det == 0.0_r8) then
!#!    print *, 'matrix cannot be inverted'
!#!    r = 0.0_r8
!#!    return
!#! endif
!#!
!#! b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
!#! b(2,1) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
!#! b(3,1) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
!#!
!#! b(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
!#! b(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
!#! b(3,2) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
!#!
!#! b(1,3) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
!#! b(2,3) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
!#! b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
!#!
!#! r = b / det
!#!
!#! end subroutine invert3
!#!
!#! !------------------------------------------------------------
!#!
!#! !==================================================================
!#! ! The following (private) routines were borrowed from the FESOM code
!#! !==================================================================
!#!
!#! !------------------------------------------------------------------
!#!
!#! subroutine r3_normalize(ax, ay, az)
!#!
!#! !normalizes the vector (ax, ay, az)
!#!
!#! real(r8), intent(inout) :: ax, ay, az
!#! real(r8) :: mi
!#!
!#!  mi = 1.0_r8 / sqrt(ax**2 + ay**2 + az**2)
!#!  ax = ax * mi
!#!  ay = ay * mi
!#!  az = az * mi
!#!
!#! end subroutine r3_normalize


!------------------------------------------------------------------

!#! function theta_to_tk (theta, rho, qv)
!#!
!#! ! Compute sensible temperature [K] from potential temperature [K].
!#! ! code matches computation done in FESOM model
!#!
!#! real(r8), intent(in)  :: theta    ! potential temperature [K]
!#! real(r8), intent(in)  :: rho      ! dry density
!#! real(r8), intent(in)  :: qv       ! water vapor mixing ratio [kg/kg]
!#! real(r8)  :: theta_to_tk          ! sensible temperature [K]
!#!
!#! ! Local variables
!#! real(r8) :: theta_m               ! potential temperature modified by qv
!#! real(r8) :: exner                 ! exner function
!#! real(r8) :: qv_nonzero            ! qv >= 0
!#!
!#! qv_nonzero = max(qv,0.0_r8)
!#! theta_m = (1.0_r8 + 1.61_r8 * qv_nonzero)*theta
!#!
!#! !theta_m = (1.0_r8 + 1.61_r8 * (max(qv, 0.0_r8)))*theta
!#! exner = ( (rgas/p0) * (rho*theta_m) )**rcv
!#!
!#! ! Temperature [K]
!#! theta_to_tk = theta * exner
!#!
!#! end function theta_to_tk


!#! !------------------------------------------------------------------
!#!
!#! subroutine compute_full_pressure(theta, rho, qv, pressure, tk)
!#!
!#! ! Compute full pressure from the equation of state.
!#! ! since it has to compute sensible temp along the way,
!#! ! make temp one of the return values rather than having
!#! ! to call theta_to_tk() separately.
!#! ! code matches computation done in FESOM model
!#!
!#! real(r8), intent(in)  :: theta    ! potential temperature [K]
!#! real(r8), intent(in)  :: rho      ! dry density
!#! real(r8), intent(in)  :: qv       ! water vapor mixing ratio [kg/kg]
!#! real(r8), intent(out) :: pressure ! full pressure [Pa]
!#! real(r8), intent(out) :: tk       ! return sensible temperature to caller
!#!
!#! ! Local variables
!#! real(r8) :: qv_nonzero            ! qv >= 0
!#!
!#! qv_nonzero = max(qv,0.0_r8)
!#! tk = theta_to_tk(theta, rho, qv_nonzero)
!#!
!#! !tk = theta_to_tk(theta, rho, max(qv,0.0_r8))
!#! pressure = rho * rgas * tk * (1.0_r8 + 1.61_r8 * qv)
!#! !if (do_output() .and. debug > 9) print *, 't,r,q,p,tk =', theta, rho, qv, pressure, tk
!#!
!#! end subroutine compute_full_pressure


!------------------------------------------------------------------
!>

subroutine dump_tables

! TJH DEBUG create the rosetta stone for locations
! separate file for temp, salinity, coord_nod3D, all dart state_vector.
! should be able to cat temp & salinity to recreate all dart_state-vector
! temp == salt == coord_nod3D, that sort of thing.

integer :: iunit, i, var_type
real(r8) :: llv(3), lon, lat, vert
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
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/FEOM/models/FeoM/model_mod.f90 $
! $Id: model_mod.f90 10268 2016-05-09 17:06:15Z thoar $
! $Revision: 10268 $
! $Date: 2016-05-09 19:06:15 +0200 (Mon, 09 May 2016) $
