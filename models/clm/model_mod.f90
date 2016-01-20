! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

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
use        types_mod, only : r4, r8, SECPERDAY, MISSING_R8,                    &
                             MISSING_I, MISSING_R4, rad2deg, deg2rad, PI,      &
                             obstypelength, i8
use time_manager_mod, only : time_type, set_time, get_time, set_date, get_date,&
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+),  operator(-),          &
                             operator(>),  operator(<),  operator(/),          &
                             operator(/=), operator(<=), operator(==)

use     location_mod, only : location_type, get_dist, query_location,          &
                             get_close_maxdist_init, get_close_type,           &
                             set_location, get_location, horiz_dist_only,      &
                             vert_is_undef,    VERTISUNDEF,                    &
                             vert_is_surface,  VERTISSURFACE,                  &
                             vert_is_level,    VERTISLEVEL,                    &
                             vert_is_pressure, VERTISPRESSURE,                 &
                             vert_is_height,   VERTISHEIGHT,                   &
                             get_close_obs_init, LocationDims,                 &
                             loc_get_close_obs => get_close_obs

use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper,                    &
                             find_namelist_in_file, check_namelist_read,       &
                             file_exist, find_textfile_dims, file_to_text,     &
                             open_file, close_file, do_nml_file, do_nml_term,  &
                             nmlfileunit

use     obs_kind_mod, only : KIND_SOIL_TEMPERATURE,       &
                             KIND_SOIL_MOISTURE,          &
                             KIND_LIQUID_WATER,           &
                             KIND_ICE,                    &
                             KIND_SNOWCOVER_FRAC,         &
                             KIND_SNOW_THICKNESS,         &
                             KIND_LEAF_CARBON,            &
                             KIND_LEAF_AREA_INDEX,        &
                             KIND_WATER_TABLE_DEPTH,      &
                             KIND_GEOPOTENTIAL_HEIGHT,    &
                             KIND_VEGETATION_TEMPERATURE, &
                             KIND_FPAR,                   &
                             KIND_FPAR_SUNLIT_DIRECT,     &
                             KIND_FPAR_SUNLIT_DIFFUSE,    &
                             KIND_FPAR_SHADED_DIRECT,     &
                             KIND_FPAR_SHADED_DIFFUSE,    &
                             KIND_FPAR_SHADED_DIRECT,     &
                             KIND_FPAR_SHADED_DIFFUSE,    &
                             paramname_length,            &
                             get_raw_obs_kind_index,      &
                             get_raw_obs_kind_name

 use ensemble_manager_mod, only : ensemble_type, &
                                  map_pe_to_task, &
                                  get_var_owner_index

use distributed_state_mod, only : get_state
use state_structure_mod,   only : add_domain

use mpi_utilities_mod, only: my_task_id

use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

!> @todo FIXME write a write_model_time for CLM
use dart_time_io_mod,      only : write_model_time

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
          pert_model_copies,      &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs, &
          write_model_time, &
          read_model_time, &
          query_vert_localization_coord, &
          vert_convert

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.

public :: get_gridsize,                 &
          clm_to_dart_state_vector,     &
          sv_to_restart_file,           &
          get_clm_restart_filename,     &
          get_grid_vertval,             &
          compute_gridcell_value,       &
          gridcell_components,          &
          DART_get_var,                 &
          get_model_time,               &
          construct_file_name_in

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=256) :: string1, string2, string3
logical, save :: module_initialized = .false.

! Storage for a random sequence for perturbing a single initial state

type(random_seq_type) :: random_seq

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

integer, parameter :: LAKE = 3

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

! things which can/should be in the model_nml

integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 60
real(r8)           :: model_perturbation_amplitude = 0.2
logical            :: output_state_vector = .true.
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
   output_state_vector,         &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   calendar,                    &
   debug,                       &
   clm_variables

! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=obstypelength), dimension(NF90_MAX_VAR_DIMS) :: dimnames
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer  :: numdims
   integer  :: maxlevels
   integer  :: xtype
   integer(i8)  :: varsize     ! prod(dimlens(1:numdims))
   integer(i8)  :: index1      ! location in dart state vector of first occurrence
   integer(i8)  :: indexN      ! location in dart state vector of last  occurrence
   integer  :: dart_kind
   integer  :: rangeRestricted
   real(r8) :: minvalue
   real(r8) :: maxvalue
   integer  :: spvalINT, missingINT
   real(r4) :: spvalR4, missingR4
   real(r8) :: spvalR8, missingR8
   logical  :: has_fill_value      ! intended for future use
   logical  :: has_missing_value   ! intended for future use
   character(len=paramname_length) :: kind_string
   character(len=512) :: origin    ! the file it came from
   logical  :: update
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

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
integer :: nlonatm  = -1
integer :: nlonrof  = -1
integer :: nlat     = -1
integer :: nlatatm  = -1
integer :: nlatrof  = -1
integer :: nlevgrnd = -1 ! Number of soil levels

real(r8), allocatable :: LEVGRND(:)
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
integer :: nrtmlon   = -1 ! Number of river transport model longitudes
integer :: nrtmlat   = -1 ! Number of river transport model latitudes
integer :: nlevcan   = -1 ! Number of canopy layers (*XY*)

integer,  allocatable, dimension(:)  :: grid1d_ixy, grid1d_jxy ! 2D lon/lat index of corresponding gridcell
integer,  allocatable, dimension(:)  :: land1d_ixy, land1d_jxy ! 2D lon/lat index of corresponding gridcell
integer,  allocatable, dimension(:)  :: cols1d_ixy, cols1d_jxy ! 2D lon/lat index of corresponding gridcell
integer,  allocatable, dimension(:)  :: pfts1d_ixy, pfts1d_jxy ! 2D lon/lat index of corresponding gridcell
real(r8), allocatable, dimension(:)  :: land1d_wtxy    ! landunit weight relative to corresponding gridcell
real(r8), allocatable, dimension(:)  :: cols1d_wtxy    ! column   weight relative to corresponding gridcell
real(r8), allocatable, dimension(:)  :: pfts1d_wtxy    ! pft      weight relative to corresponding gridcell
real(r8), allocatable, dimension(:)  :: levtot
real(r8), allocatable, dimension(:,:):: zsno           ! (column,levsno) ... snow layer midpoint
integer,  allocatable, dimension(:)  :: cols1d_ityplun ! columntype ... lake, forest, city ...

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

integer         :: model_size      ! the state vector length
type(time_type) :: model_time      ! valid time of the model state
type(time_type) :: model_timestep  ! smallest time to adv model


INTERFACE vector_to_prog_var
      MODULE PROCEDURE vector_to_1d_prog_var
      MODULE PROCEDURE vector_to_2d_prog_var
END INTERFACE

INTERFACE DART_get_var
      MODULE PROCEDURE get_var_1d
      MODULE PROCEDURE get_var_2d
      MODULE PROCEDURE get_var_3d
END INTERFACE

INTERFACE get_state_time
      MODULE PROCEDURE get_state_time_ncid
      MODULE PROCEDURE read_model_time
END INTERFACE


contains

!==================================================================
! All the REQUIRED interfaces come first - just by convention.
!==================================================================


function get_model_size()
!------------------------------------------------------------------
! Returns the size of the model as an integer.
! Required for all applications.

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size



subroutine adv_1step(x, time)
!------------------------------------------------------------------
! CLM is never advanced by DART.
! If we get here, something is wrong and we should stop right away.

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

if ( .not. module_initialized ) call static_init_model

if (do_output()) then
   call print_time(time,'NULL interface adv_1step (no advance) DART time is')
   call print_time(time,'NULL interface adv_1step (no advance) DART time is',logfileunit)
endif

write(string1,*)'DART should not be trying to advance CLM'
call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate)

! just so suppress compiler warnings. code unreachable
x(:) = MISSING_R8

end subroutine adv_1step



subroutine get_state_meta_data(state_handle, indx, location, var_type)
!------------------------------------------------------------------
! Given an integer index into the state vector structure, returns the
! associated array indices for lat, lon, and height, as well as the type.

! Passed variables

type(ensemble_type), intent(in)    :: state_handle
integer(i8), intent(in)            :: indx
type(location_type), intent(out)   :: location
integer, OPTIONAL, intent(out)     :: var_type

! Local variables

integer  :: n

! Module variables

! LON
! LAT
! lonixy
! latjxy
! levels
! progvar

if ( .not. module_initialized ) call static_init_model

location = set_location( LON(lonixy(indx)), LAT(latjxy(indx)), levels(indx), VERTISHEIGHT)

if (present(var_type)) then

   var_type = MISSING_I

   FINDTYPE : do n = 1,nfields
      if((indx >= progvar(n)%index1) .and. &
         (indx <= progvar(n)%indexN) ) then
         var_type = progvar(n)%dart_kind
         exit FINDTYPE
      endif
   enddo FINDTYPE

   if( var_type == MISSING_I ) then
      write(string1,*) 'Problem, cannot find base_offset, indx is: ', indx
      call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
   endif

endif

return
end subroutine get_state_meta_data



function get_model_time_step()
!------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = model_timestep

end function get_model_time_step



subroutine static_init_model()
!------------------------------------------------------------------
!
! Called to do one time initialization of the model.
!
! All the grid information comes from the initialization of
! the dart_clm_mod module.

! Local variables - all the important ones have module scope

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=obstypelength)          :: dimname
integer :: ncid, TimeDimID, VarID, dimlen, varsize
integer :: iunit, io, ivar
integer :: i, j, xi, xj, index1, indexN, indx
integer :: ss, dd

integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: spvalR8

character(len=paramname_length) :: var_names(max_state_variables)
real(r8) :: var_ranges(max_state_variables,2)
integer :: nvars, domid

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

!---------------------------------------------------------------
! Set the time step ... causes clm namelists to be read.
! Ensures model_timestep is multiple of 'dynamics_timestep'

call set_calendar_type( calendar )   ! comes from model_mod_nml

model_time     = get_state_time(clm_restart_filename)
model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',string1)

!---------------------------------------------------------------
! The CLM history file (h0?) has the 'superset' of information.
! The CLM restart files are intentionally lean and, in so doing,
! do not have the full lat/lon arrays nor any depth information.
!
ncid = 0; ! signal that netcdf file is closed
call get_history_dims(ncid, clm_history_filename, 'open', nlon, nlat, nlevgrnd )

if (unstructured) then
   allocate( AREA1D(nlon),      LANDFRAC1D(nlon) )
else
   allocate( AREA2D(nlon,nlat), LANDFRAC2D(nlon,nlat) )
endif

allocate(LON(nlon), LAT(nlat),  LEVGRND(nlevgrnd))

call get_full_grid(ncid, clm_history_filename, 'close')

ncid = 0; ! signal that netcdf file is closed

!---------------------------------------------------------------
! The CLM grid in a restart file is fundamentally a sparse matrix
! representation that lacks the native grid dimensions.
! The full lat/lon arrays do not exist in the restart files.
! only the grid cells that contain land are preserved.

call get_sparse_dims(ncid, clm_restart_filename, 'open')

allocate(grid1d_ixy(ngridcell), grid1d_jxy(ngridcell))
allocate(land1d_ixy(nlandunit), land1d_jxy(nlandunit), land1d_wtxy(nlandunit))
allocate(cols1d_ixy(ncolumn),   cols1d_jxy(ncolumn))
allocate(cols1d_wtxy(ncolumn),  cols1d_ityplun(ncolumn))
allocate(pfts1d_ixy(npft),      pfts1d_jxy(npft)     , pfts1d_wtxy(npft))
allocate(levtot(nlevtot))
if (nlevsno > 0) allocate(zsno(nlevsno,ncolumn))

call get_sparse_geog(ncid, clm_restart_filename, 'close')

!---------------------------------------------------------------
! Generate list of columns in each gridcell

allocate(gridCellInfo(nlon,nlat))
call SetLocatorArrays()

!---------------------------------------------------------------
! Compile the list of clm variables to use in the creation
! of the DART state vector.

call parse_variable_table( clm_variables, nfields, variable_table )

! Compute the offsets into the state vector for the start of each
! variable type. Requires reading shapes from the clm restart file.
! Record the extent of the data type in the state vector.

index1  = 1;
indexN  = 0;
do ivar = 1, nfields

   ! convey the information in the variable_table to each progvar. 

   call SetVariableAttributes(ivar)

   ! Open the file for each variable and get dimensions, etc.

   call nc_check(nf90_open(trim(progvar(ivar)%origin), NF90_NOWRITE, ncid), &
              'static_init_model','open '//trim(progvar(ivar)%origin))

   ! File is not required to have a time dimension
   io = nf90_inq_dimid(ncid, 'time', TimeDimID)
   if (io /= NF90_NOERR) TimeDimID = MISSING_I

   string2 = trim(progvar(ivar)%origin)//' '//trim(progvar(ivar)%varname)

   call nc_check(nf90_inq_varid(ncid, trim(progvar(ivar)%varname), VarID), &
            'static_init_model', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, &
                 ndims=progvar(ivar)%numdims, xtype=progvar(ivar)%xtype), &
            'static_init_model', 'inquire '//trim(string2))

   ! If the long_name and/or units attributes are set, get them.
   ! They are not REQUIRED to exist but are nice to use if they are present.

   if( nf90_inquire_attribute(    ncid, VarID, 'long_name') == NF90_NOERR ) then
      call nc_check( nf90_get_att(ncid, VarID, 'long_name' , progvar(ivar)%long_name), &
                  'static_init_model', 'get_att long_name '//trim(string2))
   else
      progvar(ivar)%long_name = progvar(ivar)%varname
   endif

   if( nf90_inquire_attribute(    ncid, VarID, 'units') == NF90_NOERR )  then
      call nc_check( nf90_get_att(ncid, VarID, 'units' , progvar(ivar)%units), &
                  'static_init_model', 'get_att units '//trim(string2))
   else
      progvar(ivar)%units = 'unknown'
   endif

   ! Saving any FillValue, missing_value attributes so I can use it when I read and write ...
   ! CESM1_1_1 ... no attributes in the restart file for rank1 or greater
   ! variables.

   if (progvar(ivar)%xtype == NF90_INT) then
       if (nf90_get_att(ncid, VarID, '_FillValue'    , spvalINT) == NF90_NOERR) then
          progvar(ivar)%spvalINT       = spvalINT
          progvar(ivar)%has_fill_value = .true.
       endif
       if (nf90_get_att(ncid, VarID, 'missing_value' , spvalINT) == NF90_NOERR) then
          progvar(ivar)%missingINT        = spvalINT
          progvar(ivar)%has_missing_value = .true.
       endif

   elseif (progvar(ivar)%xtype == NF90_FLOAT) then
       if (nf90_get_att(ncid, VarID, '_FillValue'    , spvalR4) == NF90_NOERR) then
          progvar(ivar)%spvalR4        = spvalR4
          progvar(ivar)%has_fill_value = .true.
       endif
       if (nf90_get_att(ncid, VarID, 'missing_value' , spvalR4) == NF90_NOERR) then
          progvar(ivar)%missingR4         = spvalR4
          progvar(ivar)%has_missing_value = .true.
       endif

   elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
       if (nf90_get_att(ncid, VarID, '_FillValue'    , spvalR8) == NF90_NOERR) then
          progvar(ivar)%spvalR8        = spvalR8
          progvar(ivar)%has_fill_value = .true.
       endif
       if (nf90_get_att(ncid, VarID, 'missing_value' , spvalR8) == NF90_NOERR) then
          progvar(ivar)%missingR8         = spvalR8
          progvar(ivar)%has_missing_value = .true.
       endif
   endif

   varsize = 1
   dimlen  = 1
   DimensionLoop : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), name=dimname, len=dimlen), &
                                          'static_init_model', string1)

      ! Only reserve space for a single time slice 
      if (dimIDs(i) == TimeDimID) dimlen = 1

      progvar(ivar)%dimlens( i) = dimlen
      progvar(ivar)%dimnames(i) = dimname
      varsize = varsize * dimlen

   enddo DimensionLoop

   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1
   index1                    = index1 + varsize      ! sets up for next variable

   if ((debug > 0) .and. do_output()) then
      write(logfileunit,*)
      write(logfileunit,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(logfileunit,*) '  filename    ',trim(progvar(ivar)%origin)
      write(logfileunit,*) '  update      ',progvar(ivar)%update
      write(logfileunit,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(logfileunit,*) '  units       ',trim(progvar(ivar)%units)
      write(logfileunit,*) '  xtype       ',progvar(ivar)%xtype
      write(logfileunit,*) '  dimnames    ',(/ (trim(progvar(ivar)%dimnames(i))//' ', i=1,progvar(ivar)%numdims ) /)
      write(logfileunit,*) '  dimlens     ',progvar(ivar)%dimlens( 1:progvar(ivar)%numdims)
      write(logfileunit,*) '  numdims     ',progvar(ivar)%numdims
      write(logfileunit,*) '  varsize     ',progvar(ivar)%varsize
      write(logfileunit,*) '  index1      ',progvar(ivar)%index1
      write(logfileunit,*) '  indexN      ',progvar(ivar)%indexN
      write(logfileunit,*) '  dart_kind   ',progvar(ivar)%dart_kind
      write(logfileunit,*) '  kind_string ',progvar(ivar)%kind_string
      write(logfileunit,*) '  spvalINT    ',progvar(ivar)%spvalINT
      write(logfileunit,*) '  spvalR4     ',progvar(ivar)%spvalR4
      write(logfileunit,*) '  spvalR8     ',progvar(ivar)%spvalR8
      write(logfileunit,*) '  missingINT  ',progvar(ivar)%missingINT
      write(logfileunit,*) '  missingR4   ',progvar(ivar)%missingR4
      write(logfileunit,*) '  missingR8   ',progvar(ivar)%missingR8
      write(logfileunit,*) '  has_fill_value    ',progvar(ivar)%has_fill_value
      write(logfileunit,*) '  has_missing_value ',progvar(ivar)%has_missing_value
      write(logfileunit,*)'   rangeRestricted   ',progvar(ivar)%rangeRestricted
      write(logfileunit,*)'   minvalue          ',progvar(ivar)%minvalue
      write(logfileunit,*)'   maxvalue          ',progvar(ivar)%maxvalue

      write(     *     ,*)
      write(     *     ,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(     *     ,*) '  filename    ',trim(progvar(ivar)%origin)
      write(     *     ,*) '  update      ',progvar(ivar)%update
      write(     *     ,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(     *     ,*) '  units       ',trim(progvar(ivar)%units)
      write(     *     ,*) '  xtype       ',progvar(ivar)%xtype
      write(     *     ,*) '  dimnames    ',(/ (trim(progvar(ivar)%dimnames(i))//' ', i=1,progvar(ivar)%numdims ) /)
      write(     *     ,*) '  dimlens     ',progvar(ivar)%dimlens( 1:progvar(ivar)%numdims)
      write(     *     ,*) '  numdims     ',progvar(ivar)%numdims
      write(     *     ,*) '  varsize     ',progvar(ivar)%varsize
      write(     *     ,*) '  index1      ',progvar(ivar)%index1
      write(     *     ,*) '  indexN      ',progvar(ivar)%indexN
      write(     *     ,*) '  dart_kind   ',progvar(ivar)%dart_kind
      write(     *     ,*) '  kind_string ',progvar(ivar)%kind_string
      write(     *     ,*) '  spvalINT    ',progvar(ivar)%spvalINT
      write(     *     ,*) '  spvalR4     ',progvar(ivar)%spvalR4
      write(     *     ,*) '  spvalR8     ',progvar(ivar)%spvalR8
      write(     *     ,*) '  missingINT  ',progvar(ivar)%missingINT
      write(     *     ,*) '  missingR4   ',progvar(ivar)%missingR4
      write(     *     ,*) '  missingR8   ',progvar(ivar)%missingR8
      write(     *     ,*) '  has_fill_value    ',progvar(ivar)%has_fill_value
      write(     *     ,*) '  has_missing_value ',progvar(ivar)%has_missing_value
      write(     *     ,*)'   rangeRestricted   ',progvar(ivar)%rangeRestricted
      write(     *     ,*)'   minvalue          ',progvar(ivar)%minvalue
      write(     *     ,*)'   maxvalue          ',progvar(ivar)%maxvalue
   endif

   call nc_check(nf90_close(ncid),'static_init_model','close '//trim(string2))
   ncid = 0

enddo

model_size = progvar(nfields)%indexN

if ((debug > 0) .and. do_output()) then
  write(logfileunit, *)
  write(logfileunit,'("grid: nlon, nlat, nz =",2(1x,i6))') nlon, nlat
  write(logfileunit, *)'model_size = ', model_size
  write(     *     , *)
  write(     *     ,'("grid: nlon, nlat, nz =",2(1x,i6))') nlon, nlat
  write(     *     , *)'model_size = ', model_size
endif

!---------------------------------------------------------------
! Create the metadata arrays that are the same shape as the state vector.
! The metadata arrays will provide the ability to determine what grid cell is the parent
! of the state vector index in question ... as well as the actual surface area.
! This MUST stride through the state vector the same way the state vector is filled.

allocate(lonixy(model_size), latjxy(model_size), levels(model_size), landarea(model_size))

! Initialize all levels to surface. If there is a level, we will explicitly specify it.
levels(:) = 0.0_r8

do ivar=1, nfields

   ! All variables are at the surface until proven otherwise.
   progvar(ivar)%maxlevels = 1

   indx = progvar(ivar)%index1

   ! 1D variables are usually from the restart file
   ! FIXME this same logic is used for 2D variables from the 'vector' file which
   ! has a singleton dimension of 'time'
   ! What if I check on the rank of the variable instead of the number of dimensions ...
   ! If I require 'time' to be the unlimited dimension, it will always be 'last',
   ! so the first N dimensions and the first N ranks are identical ...

   if (progvar(ivar)%numdims == 1) then

      if ((debug > 8) .and. do_output()) then
         write(*,*)
         write(*,*)'variable ',trim(progvar(ivar)%varname)
         write(*,*)'dimension 1 (i) ',progvar(ivar)%dimnames(1),progvar(ivar)%dimlens(1)
      endif

      SELECT CASE ( trim(progvar(ivar)%dimnames(1)) )
         CASE ("gridcell","lndgrid")
            do i = 1, progvar(ivar)%dimlens(1)
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
            do i = 1, progvar(ivar)%dimlens(1)
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
            do i = 1, progvar(ivar)%dimlens(1)
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
            do i = 1, progvar(ivar)%dimlens(1)
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
            write(string1,*)'(1d) unknown Dimension name '//trim(progvar(ivar)%dimnames(1))// &
             ' while trying to create metadata for '//trim(progvar(ivar)%varname)
            call error_handler(E_ERR,'static_init_model', string1, source, revision, revdate)

      END SELECT

   elseif (progvar(ivar)%numdims == 2) then

      ! restart file variables may have 2 dimensions, one of which is layers.
      ! vector_history files may have 2 dimensions, one of which is time
      ! history file variables always have 3 dimension [time,lat,lon]

      ! In the ncdump output, dimension 2 is the leftmost dimension.
      ! Only dimension 2 matters for the weights.

      if ((debug > 8) .and. do_output()) then
         write(*,*)
         write(*,*)'variable ',trim(progvar(ivar)%varname)
         write(*,*)'dimension 1 (i) ',progvar(ivar)%dimnames(1),progvar(ivar)%dimlens(1)
         write(*,*)'dimension 2 (j) ',progvar(ivar)%dimnames(2),progvar(ivar)%dimlens(2)
      endif

      SELECT CASE ( trim(progvar(ivar)%dimnames(2)) )
         CASE ("gridcell")
            if ((debug > 8) .and. do_output()) write(*,*)'length grid1d_ixy ',size(grid1d_ixy)
            do j = 1, progvar(ivar)%dimlens(2)
               xi = grid1d_ixy(j)
               xj = grid1d_jxy(j) ! nnnnn_jxy(:) always 1 if unstructured
               do i = 1, progvar(ivar)%dimlens(1)
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
            do j = 1, progvar(ivar)%dimlens(2)
               xi = land1d_ixy(j)
               xj = land1d_jxy(j) ! nnnnn_jxy(:) always 1 if unstructured
               do i = 1, progvar(ivar)%dimlens(1)
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
                         ! The vertical levels are fully defined by the levgrnd and
                         ! levsno variables. levgrnd is static, levsno varies by column.

            progvar(ivar)%maxlevels = progvar(ivar)%dimlens(2)

            if ((debug > 8) .and. do_output()) write(*,*)'length cols1d_ixy ',size(cols1d_ixy)
            if ((debug > 8) .and. do_output()) write(*,*)'size zsno ',size(zsno,1), size(zsno,2)
            if ((debug > 8) .and. do_output()) write(*,*)'nlevsno ',nlevsno

            LANDCOLUMN : do j = 1, progvar(ivar)%dimlens(2)

               call fill_levels(progvar(ivar)%dimnames(1),j,progvar(ivar)%dimlens(1),levtot)

               xi = cols1d_ixy(j)
               xj = cols1d_jxy(j) ! nnnnn_jxy(:) always 1 if unstructured
               VERTICAL :  do i = 1, progvar(ivar)%dimlens(1)
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
               enddo VERTICAL
            enddo LANDCOLUMN

         CASE ("pft")
            if ((debug > 8) .and. do_output()) write(*,*)'length pfts1d_ixy ',size(pfts1d_ixy)
            do j = 1, progvar(ivar)%dimlens(2)
               xi = pfts1d_ixy(j)
               xj = pfts1d_jxy(j) ! nnnnn_jxy(:) always 1 if unstructured
               do i = 1, progvar(ivar)%dimlens(1)
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

         CASE ("time")

            ! The vector history files can have things 'pft,time' or 'column,time'
            ! The single-column history files can have things 'lndgrid,time', 'pft,time' or 'column,time'

            SELECT CASE ( trim(progvar(ivar)%dimnames(1)) )
               CASE ("pft")
                  do i = 1, progvar(ivar)%dimlens(1)
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

               CASE ("column")
                  do i = 1, progvar(ivar)%dimlens(1)
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

               CASE ("lndgrid")
                  do i = 1, progvar(ivar)%dimlens(1)
                     xi             = cols1d_ixy(i)
                     xj             = cols1d_jxy(i) ! nnnnn_jxy(:) always 1 if unstructured
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

               CASE DEFAULT
                  write(string1,*)'(2d) unknown first Dimension name '//trim(progvar(ivar)%dimnames(1))// &
                   ' while trying to create metadata for '//trim(progvar(ivar)%varname)
                  call error_handler(E_ERR,'static_init_model', string1, source, revision, revdate)
            END SELECT

         CASE DEFAULT
            write(string1,*)'(2d) unknown Dimension name '//trim(progvar(ivar)%dimnames(2))// &
             ' while trying to create metadata for '//trim(progvar(ivar)%varname)
            call error_handler(E_ERR,'static_init_model', string1, source, revision, revdate)

      END SELECT


   elseif (progvar(ivar)%numdims == 3) then

      ! restart file variables never have 3 dimensions
      ! vector_history variables   may have 3 dimensions [time, lat, lon]
      ! history file variables always have 3 dimensions [time, lat, lon]
      !     exception is float H2OSOI(time, levgrnd, lat, lon) ... but we
      !     have access to restart file h2osoi_[liq,ice]

      if ((debug > 8) .and. do_output()) then
         write(*,*)
         write(*,*)'variable ',trim(progvar(ivar)%varname)
         write(*,*)'dimension 1 (i) ',progvar(ivar)%dimnames(1),progvar(ivar)%dimlens(1)
         write(*,*)'dimension 2 (j) ',progvar(ivar)%dimnames(2),progvar(ivar)%dimlens(2)
         write(*,*)'dimension 3 (k) ',progvar(ivar)%dimnames(3),progvar(ivar)%dimlens(3)
      endif

      ! Remember the order is reversed from ncdump to fortran
      ! The messages are in ncdump-order, checking is fortran-order
      if ((trim(progvar(ivar)%dimnames(1)) .ne. 'lon') .or. &
          (trim(progvar(ivar)%dimnames(2)) .ne. 'lat' ) .or. &
          (trim(progvar(ivar)%dimnames(3)) .ne. 'time' )) then
         write(string1,*)'3D variables must be [time,lat,lon] (as reported by ncdump)'
         write(string2,*)trim(progvar(ivar)%varname),' is ', &
           trim(progvar(ivar)%dimnames(3)), ' ', &
           trim(progvar(ivar)%dimnames(2)), ' ', &
           trim(progvar(ivar)%dimnames(1))
         call error_handler(E_ERR,'static_init_model', string1, &
                 source, revision, revdate, text2=string2)
      endif

      ! The get_var_3d() routine ensures there is only 1 timestep.
      ! So there is no need to loop over dimlens(3) 

      do j = 1, progvar(ivar)%dimlens(2)     ! time-invariant
         do i = 1, progvar(ivar)%dimlens(1)  ! time-invariant
            lonixy(  indx) = i
            latjxy(  indx) = j
            landarea(indx) = AREA2D(i,j) * LANDFRAC2D(i,j)
            indx = indx + 1
         enddo
      enddo

   else

      ! Cannot support 4D variables yet.
      ! These only occurr in 2D history files for variables with levels (and time)
      ! i.e. H2OSOI(time, levgrnd, lat, lon)
      !
      ! You can get the same information from the restart file, usually.

      write(string1,*)'variables of rank ',progvar(ivar)%numdims,' are unsupported.'
      write(string2,*)trim(progvar(ivar)%varname),' is dimensioned ',&
                           progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      call error_handler(E_ERR,'static_init_model', string1, source, revision, &
                         revdate, text2=string2)

   endif

   indx = indx - 1
   if (indx /= progvar(ivar)%indexN ) then
      write(string1,*)'variable ',trim(progvar(ivar)%varname), &
       ' is supposed to end at index ',progvar(ivar)%indexN
      write(string2,*)'it ends at index ',indx
      call error_handler(E_ERR,'static_init_model', string1, source, revision, &
                         revdate, text2=string2)
   endif

enddo

! Must group the variables according to the file they come from.
!> @TODO FIXME ... add_domain() must do the right thing if nvars == 0
!> @TODO FIXME ... io_filenames_nml:rpointer_file order must somehow match 
!> - or be insensitive to - the add_domain() calls below (if nvars == 0) ...

call cluster_variables(clm_restart_filename, nvars, var_names,            var_ranges)
domid =     add_domain(clm_restart_filename, nvars, var_names, clamp_vals=var_ranges)

call cluster_variables(clm_history_filename, nvars, var_names,            var_ranges)
domid =     add_domain(clm_history_filename, nvars, var_names, clamp_vals=var_ranges)

call cluster_variables(clm_vector_history_filename, nvars, var_names,            var_ranges)
domid =     add_domain(clm_vector_history_filename, nvars, var_names, clamp_vals=var_ranges)

end subroutine static_init_model



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

! if ( .not. module_initialized ) call static_init_model

if (unstructured) then
   deallocate(AREA1D, LANDFRAC1D)
else
   deallocate(AREA2D, LANDFRAC2D)
endif

deallocate(LAT, LON, LEVGRND)
deallocate(grid1d_ixy, grid1d_jxy)
deallocate(land1d_ixy, land1d_jxy, land1d_wtxy)
deallocate(cols1d_ixy, cols1d_jxy, cols1d_wtxy, cols1d_ityplun)
deallocate(pfts1d_ixy, pfts1d_jxy, pfts1d_wtxy)

deallocate(lonixy, latjxy, landarea)

end subroutine end_model



subroutine init_time(time)
!------------------------------------------------------------------
!
! Companion interface to init_conditions. Returns a time that is somehow
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no
! synthetic data experiments using perfect_model_obs are planned,
! this can be a NULL INTERFACE.

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

! for now, just set to 0
time = set_time(0,0)

end subroutine init_time



subroutine init_conditions(x)
!------------------------------------------------------------------
!
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no
! synthetic data experiments using perfect_model_obs are planned,
! this can be a NULL INTERFACE.

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

x = 0.0_r8

end subroutine init_conditions



function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
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

integer :: StateVarVarID   ! netCDF pointer to state variable coordinate array
integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

!----------------------------------------------------------------------
! variables if we parse the state vector into prognostic variables.
!----------------------------------------------------------------------

! for the dimensions and coordinate variables
integer ::     nlonDimID
integer ::     nlatDimID
integer ::  lndgridDimID
integer :: gridcellDimID
integer :: landunitDimID
integer ::   columnDimID
integer ::      pftDimID
integer ::  levgrndDimID
integer ::   levlakDimID
integer ::   levsnoDimID
integer ::  levsno1DimID
integer ::   levtotDimID
integer ::   numradDimID
integer ::   levcanDimID

! for the prognostic variables
integer :: ivar, VarID

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
integer :: i, myndims

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
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'clm' ), &
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

  ! Define the state vector coordinate variable and some attributes.
   call nc_check(nf90_def_var(ncid=ncFileID,name='StateVariable', xtype=nf90_int, &
                 dimids=StateVarDimID, varid=StateVarVarID), 'nc_write_model_atts', &
                 'statevariable def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarVarID,'long_name','State Variable ID'),&
                 'nc_write_model_atts','statevariable long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, 'units','indexical'), &
                 'nc_write_model_atts', 'statevariable units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarVarID,'valid_range',(/ 1,model_size /)),&
                 'nc_write_model_atts', 'statevariable valid_range '//trim(filename))

   ! Define the actual (3D) state vector, which gets filled as time goes on ...
   call nc_check(nf90_def_var(ncid=ncFileID, name='state', xtype=nf90_real, &
                 dimids=(/StateVarDimID,MemberDimID,unlimitedDimID/),varid=StateVarID),&
                 'nc_write_model_atts','state def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarID,'long_name','model state or fcopy'),&
                 'nc_write_model_atts', 'state long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarID,'_FillValue',MISSING_R8),&
                 'nc_write_model_atts', 'state FillValue '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarID,'missing_value',MISSING_R8),&
                 'nc_write_model_atts', 'state missing_value '//trim(filename))

   ! Leave define mode so we can fill the coordinate variable.
   call nc_check(nf90_enddef(ncfileID),'nc_write_model_atts','state enddef '//trim(filename))

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ), &
                 'nc_write_model_atts', 'state put_var '//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to output the prognostic variables.
   !----------------------------------------------------------------------------
   ! Define the new dimensions IDs
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_dim(ncid=ncFileID, name='lon', len = nlon, &
             dimid =     nlonDimID),'nc_write_model_atts', 'lon def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='lat', len = nlat, &
             dimid =     nlatDimID),'nc_write_model_atts', 'lat def_dim '//trim(filename))
   if (unstructured) then
      call nc_check(nf90_def_dim(ncid=ncFileID, name='lndgrid', len = ngridcell, &
             dimid = lndgridDimID),'nc_write_model_atts', 'lndgrid def_dim '//trim(filename))
   endif

   call nc_check(nf90_def_dim(ncid=ncFileID, name='gridcell', len = ngridcell, &
             dimid = gridcellDimID),'nc_write_model_atts', 'gridcell def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='landunit', len = nlandunit, &
             dimid = landunitDimID),'nc_write_model_atts', 'landunit def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='column', len = ncolumn, &
             dimid =   columnDimID),'nc_write_model_atts', 'column def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='pft', len = npft, &
             dimid =      pftDimID),'nc_write_model_atts', 'pft def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='levgrnd',  len = nlevgrnd, &
             dimid =  levgrndDimID),'nc_write_model_atts', 'levgrnd def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='levlak', len = nlevlak, &
             dimid =   levlakDimID),'nc_write_model_atts', 'levlak def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='levsno', len = nlevsno, &
             dimid =   levsnoDimID),'nc_write_model_atts', 'levsno def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='levsno1', len = nlevsno1, &
             dimid =  levsno1DimID),'nc_write_model_atts', 'levsno1 def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='levtot', len = nlevtot, &
             dimid =   levtotDimID),'nc_write_model_atts', 'levtot def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='numrad', len = nnumrad, &
             dimid =   numradDimID),'nc_write_model_atts', 'numrad def_dim '//trim(filename))
   if (nlevcan > 0) &
   call nc_check(nf90_def_dim(ncid=ncFileID, name='levcan', len = nlevcan, &
             dimid =   levcanDimID),'nc_write_model_atts', 'levcan def_dim'//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Coordinate Variables and the Attributes
   !----------------------------------------------------------------------------

   ! Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='lon', xtype=nf90_real, &
                 dimids=(/ nlonDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lon def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'coordinate longitude'), &
                 'nc_write_model_atts', 'lon long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'X'),  &
                 'nc_write_model_atts', 'lon cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'lon units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'lon valid_range '//trim(filename))

   ! Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='lat', xtype=nf90_real, &
                 dimids=(/ nlatDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lat def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'coordinate latitude'), &
                 'nc_write_model_atts', 'lat long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Y'),   &
                 'nc_write_model_atts', 'lat cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_north'),  &
                 'nc_write_model_atts', 'lat units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'lat valid_range '//trim(filename))

   ! subsurface levels
   call nc_check(nf90_def_var(ncFileID,name='levgrnd', xtype=nf90_real, &
                 dimids=(/ levgrndDimID /), varid=VarID),&
                 'nc_write_model_atts', 'levgrnd def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'coordinate soil levels'), &
                 'nc_write_model_atts', 'levgrnd long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Z'),   &
                 'nc_write_model_atts', 'levgrnd cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'm'),  &
                 'nc_write_model_atts', 'levgrnd units '//trim(filename))

   ! grid cell areas
   if (unstructured) then
      call nc_check(nf90_def_var(ncFileID,name='area', xtype=nf90_real, &
                 dimids=(/ nlonDimID /), varid=VarID),&
                 'nc_write_model_atts', 'area def_var '//trim(filename))
   else
      call nc_check(nf90_def_var(ncFileID,name='area', xtype=nf90_real, &
                 dimids=(/ nlonDimID,nlatDimID /), varid=VarID),&
                 'nc_write_model_atts', 'area def_var '//trim(filename))
   endif
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'grid cell areas'), &
                 'nc_write_model_atts', 'area long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'km^2'),  &
                 'nc_write_model_atts', 'area units '//trim(filename))

   ! grid cell land fractions
   if (unstructured) then
      call nc_check(nf90_def_var(ncFileID,name='landfrac', xtype=nf90_real, &
                 dimids=(/ nlonDimID /), varid=VarID),&
                 'nc_write_model_atts', 'landfrac def_var '//trim(filename))
   else
      call nc_check(nf90_def_var(ncFileID,name='landfrac', xtype=nf90_real, &
                 dimids=(/ nlonDimID,nlatDimID /), varid=VarID),&
                 'nc_write_model_atts', 'landfrac def_var '//trim(filename))
   endif
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'land fraction'), &
                 'nc_write_model_atts', 'landfrac long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'km^2'),  &
                 'nc_write_model_atts', 'landfrac units '//trim(filename))

   ! longitude grid index for each column
   call nc_check(nf90_def_var(ncFileID,name='cols1d_ixy', xtype=nf90_int, &
                 dimids=(/ columnDimID /), varid=VarID),&
                 'nc_write_model_atts', 'cols1d_ixy def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', &
                 '2d longitude index of corresponding column'), &
                 'nc_write_model_atts', 'cols1d_ixy long_name '//trim(filename))

   ! latitude grid index for each column
   call nc_check(nf90_def_var(ncFileID,name='cols1d_jxy', xtype=nf90_int, &
                 dimids=(/ columnDimID /), varid=VarID),&
                 'nc_write_model_atts', 'cols1d_jxy def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', &
                 '2d latitude index of corresponding column'), &
                 'nc_write_model_atts', 'cols1d_jxy long_name '//trim(filename))

   ! column weight relative to corresponding gridcell
   call nc_check(nf90_def_var(ncFileID,name='cols1d_wtxy', xtype=nf90_double, &
                 dimids=(/ columnDimID /), varid=VarID),&
                 'nc_write_model_atts', 'cols1d_wtxy def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', &
                 'column weight relative to corresponding gridcell'), &
                 'nc_write_model_atts', 'cols1d_wtxy long_name '//trim(filename))

   ! land type to corresponding gridcell
   call nc_check(nf90_def_var(ncFileID,name='cols1d_ityplun', xtype=nf90_int, &
                 dimids=(/ columnDimID /), varid=VarID),&
                 'nc_write_model_atts', 'cols1d_ityplun def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', &
                 'column landunit type (vegetated,urban,lake,wetland or glacier)'), &
                 'nc_write_model_atts', 'cols1d_ityplun long_name '//trim(filename))

   ! longitude grid index for each pft
   call nc_check(nf90_def_var(ncFileID,name='pfts1d_ixy', xtype=nf90_int, &
                 dimids=(/ pftDimID /), varid=VarID),&
                 'nc_write_model_atts', 'pfts1d_ixy def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', &
                 '2d longitude index of corresponding column'), &
                 'nc_write_model_atts', 'pfts1d_ixy long_name '//trim(filename))

   ! latitude grid index for each pft
   call nc_check(nf90_def_var(ncFileID,name='pfts1d_jxy', xtype=nf90_int, &
                 dimids=(/ pftDimID /), varid=VarID),&
                 'nc_write_model_atts', 'pfts1d_jxy def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', &
                 '2d latitude index of corresponding column'), &
                 'nc_write_model_atts', 'pfts1d_jxy long_name '//trim(filename))

   ! pft weight relative to corresponding gridcell
   call nc_check(nf90_def_var(ncFileID,name='pfts1d_wtxy', xtype=nf90_double, &
                 dimids=(/ pftDimID /), varid=VarID),&
                 'nc_write_model_atts', 'pfts1d_wtxy def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', &
                 'pft weight relative to corresponding gridcell'), &
                 'nc_write_model_atts', 'pfts1d_wtxy long_name '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------

   do ivar=1, nfields

      varname = trim(progvar(ivar)%varname)
      string1 = trim(filename)//' '//trim(varname)

      ! match shape of the variable to the dimension IDs

      call define_var_dims(ivar, ncFileID, MemberDimID, unlimitedDimID, myndims, mydimids)

      ! define the variable and set the attributes

      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(varname), &
                    xtype=progvar(ivar)%xtype, &
                    dimids = mydimids(1:myndims), varid=VarID),&
                    'nc_write_model_atts', trim(string1)//' def_var' )

      call nc_check(nf90_put_att(ncFileID, VarID, &
              'long_name', trim(progvar(ivar)%long_name)), &
              'nc_write_model_atts', trim(string1)//' put_att long_name' )
      call nc_check(nf90_put_att(ncFileID, VarID, &
              'DART_kind', trim(progvar(ivar)%kind_string)), &
              'nc_write_model_atts', trim(string1)//' put_att dart_kind' )
      call nc_check(nf90_put_att(ncFileID, VarID, &
              'units', trim(progvar(ivar)%units)), &
              'nc_write_model_atts', trim(string1)//' put_att units' )

      ! Preserve the original missing_value/_FillValue code.

      if (  progvar(ivar)%xtype == NF90_INT ) then
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 'missing_value', progvar(ivar)%spvalINT), &
                 'nc_write_model_atts', trim(string1)//' put_att missing_value' )
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 '_FillValue',  progvar(ivar)%spvalINT), &
                 'nc_write_model_atts', trim(string1)//' put_att _FillValue' )

      elseif (  progvar(ivar)%xtype == NF90_FLOAT ) then
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 'missing_value', progvar(ivar)%spvalR4), &
                 'nc_write_model_atts', trim(string1)//' put_att missing_value' )
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 '_FillValue',  progvar(ivar)%spvalR4), &
                 'nc_write_model_atts', trim(string1)//' put_att _FillValue' )

      elseif (  progvar(ivar)%xtype == NF90_DOUBLE ) then
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 'missing_value', progvar(ivar)%spvalR8), &
                 'nc_write_model_atts', trim(string1)//' put_att missing_value' )
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 '_FillValue',  progvar(ivar)%spvalR8), &
                 'nc_write_model_atts', trim(string1)//' put_att _FillValue' )
      endif

   enddo

   !----------------------------------------------------------------------------
   ! Finished with dimension/variable definitions, must end 'define' mode to fill.
   !----------------------------------------------------------------------------

   call nc_check(nf90_enddef(ncfileID), 'prognostic enddef '//trim(filename))

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables
   !----------------------------------------------------------------------------

   call nc_check(nf90_inq_varid(ncFileID, 'lon', VarID), &
                'nc_write_model_atts', 'put_var lon '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, lon ), &
                'nc_write_model_atts', 'lon put_var '//trim(filename))

   call nc_check(nf90_inq_varid(ncFileID, 'lat', VarID), &
                'nc_write_model_atts', 'put_var lat '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, lat ), &
                'nc_write_model_atts', 'lat put_var '//trim(filename))

   call nc_check(nf90_inq_varid(ncFileID, 'levgrnd', VarID), &
                'nc_write_model_atts', 'put_var levgrnd '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, LEVGRND ), &
                'nc_write_model_atts', 'levgrnd put_var '//trim(filename))

   ! AREA can be 1D or 2D
   call nc_check(nf90_inq_varid(ncFileID, 'area', VarID), &
                'nc_write_model_atts', 'put_var area '//trim(filename))
   if (unstructured) then
      call nc_check(nf90_put_var(ncFileID, VarID, AREA1D ), &
                'nc_write_model_atts', 'area put_var '//trim(filename))
   else
      call nc_check(nf90_put_var(ncFileID, VarID, AREA2D ), &
                'nc_write_model_atts', 'area put_var '//trim(filename))
   endif


   ! LANDFRAC can be 1D or 2D
   call nc_check(nf90_inq_varid(ncFileID, 'landfrac', VarID), &
                'nc_write_model_atts', 'put_var landfrac '//trim(filename))
   if (unstructured) then
      call nc_check(nf90_put_var(ncFileID, VarID, LANDFRAC1D ), &
                'nc_write_model_atts', 'landfrac put_var '//trim(filename))
   else
      call nc_check(nf90_put_var(ncFileID, VarID, LANDFRAC2D ), &
                'nc_write_model_atts', 'landfrac put_var '//trim(filename))
   endif

   call nc_check(nf90_inq_varid(ncFileID, 'cols1d_ixy', VarID), &
                'nc_write_model_atts', 'put_var cols1d_ixy '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, cols1d_ixy ), &
                'nc_write_model_atts', 'cols1d_ixy put_var '//trim(filename))

   call nc_check(nf90_inq_varid(ncFileID, 'cols1d_jxy', VarID), &
                'nc_write_model_atts', 'put_var cols1d_jxy '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, cols1d_jxy ), &
                'nc_write_model_atts', 'cols1d_jxy put_var '//trim(filename))

   call nc_check(nf90_inq_varid(ncFileID, 'cols1d_wtxy', VarID), &
                'nc_write_model_atts', 'put_var cols1d_wtxy '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, cols1d_wtxy ), &
                'nc_write_model_atts', 'cols1d_wtxy put_var '//trim(filename))

   call nc_check(nf90_inq_varid(ncFileID, 'cols1d_ityplun', VarID), &
                'nc_write_model_atts', 'put_var cols1d_ityplun '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, cols1d_ityplun ), &
                'nc_write_model_atts', 'cols1d_ityplun put_var '//trim(filename))

   call nc_check(nf90_inq_varid(ncFileID, 'pfts1d_ixy', VarID), &
                'nc_write_model_atts', 'put_var pfts1d_ixy '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, pfts1d_ixy ), &
                'nc_write_model_atts', 'pfts1d_ixy put_var '//trim(filename))

   call nc_check(nf90_inq_varid(ncFileID, 'pfts1d_jxy', VarID), &
                'nc_write_model_atts', 'put_var pfts1d_jxy '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, pfts1d_jxy ), &
                'nc_write_model_atts', 'pfts1d_jxy put_var '//trim(filename))

   call nc_check(nf90_inq_varid(ncFileID, 'pfts1d_wtxy', VarID), &
                'nc_write_model_atts', 'put_var pfts1d_wtxy '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, pfts1d_wtxy ), &
                'nc_write_model_atts', 'pfts1d_wtxy put_var '//trim(filename))

endif

!-------------------------------------------------------------------------------
! Fill the variables we can
!-------------------------------------------------------------------------------

! none

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts



function nc_write_model_vars( ncFileID, state_vec, copyindex, timeindex ) result (ierr)
!------------------------------------------------------------------
! Writes the model variables to a netCDF file.
!
! All errors are fatal, so the
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

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, ncstart, nccount
character(len=NF90_MAX_NAME)          :: varname
integer :: i, ivar, VarID, ncNdims, dimlen
integer :: TimeDimID, CopyDimID

real(r8), allocatable, dimension(:)       :: data_1d_array
real(r8), allocatable, dimension(:,:)     :: data_2d_array

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

   do ivar = 1,nfields  ! Very similar to loop in sv_to_restart_file

      varname = trim(progvar(ivar)%varname)
      string2 = trim(filename)//' '//trim(varname)

      ! Ensure netCDF variable is conformable with progvar quantity.
      ! The TIME and Copy dimensions are intentionally not queried
      ! by looping over the dimensions stored in the progvar type.

      call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'nc_write_model_vars', 'inq_varid '//trim(string2))

      call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
            'nc_write_model_vars', 'inquire '//trim(string2))

      ncstart = 1   ! These are arrays, actually
      nccount = 1
      DimCheck : do i = 1,progvar(ivar)%numdims

         write(string1,'(a,i2,A)') 'inquire dimension ',i,trim(string2)
         call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
               'nc_write_model_vars', trim(string1))

         if (progvar(ivar)%dimnames(i) == 'time') cycle DimCheck

         if ( dimlen /= progvar(ivar)%dimlens(i) ) then
            write(string1,*)trim(string2),' dim/dimlen ',i,dimlen, &
                            ' not ',progvar(ivar)%dimlens(i)
            write(string2,*)' but it should be.'
            call error_handler(E_ERR, 'nc_write_model_vars', trim(string1), &
                            source, revision, revdate, text2=trim(string2))
         endif

         nccount(i) = dimlen

      enddo DimCheck

      where(dimIDs == CopyDimID) ncstart = copyindex
      where(dimIDs == CopyDimID) nccount = 1
      where(dimIDs == TimeDimID) ncstart = timeindex
      where(dimIDs == TimeDimID) nccount = 1

      if ((debug > 0) .and. do_output()) then
         write(*,*)'nc_write_model_vars '//trim(varname)//' start is ',ncstart(1:ncNdims)
         write(*,*)'nc_write_model_vars '//trim(varname)//' count is ',nccount(1:ncNdims)
      endif

      ! Since 'time' is a singleton dimension, we can use the same logic
      ! as if it the variable had one less dimension.

      if (     (progvar(ivar)%numdims == 1) .or. &
              ((progvar(ivar)%numdims == 2) .and. &
               (progvar(ivar)%dimnames(2) == 'time')) )then

         if ( ncNdims /= 3 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 3 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_1d_array( progvar(ivar)%dimlens(1) ) )
         call vector_to_prog_var(state_vec, ivar, data_1d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array, &
             start = ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_1d_array)

      elseif ( (progvar(ivar)%numdims == 2) .or. &
              ((progvar(ivar)%numdims == 3) .and. &
               (progvar(ivar)%dimnames(3) == 'time')) )then

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
             start = ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_2d_array)

      else

         write(string1,*)'do not know how to handle CLM variables with more than 2 dimensions'
         write(string2,*)trim(progvar(ivar)%varname),'has shape', &
                              progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
         call error_handler(E_ERR,'nc_write_model_vars',string1,source,revision,revdate)

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

subroutine pert_model_copies(ens_handle, ens_size, pert_amp, interf_provided)

!
! Perturbs a single model state for generating initial ensembles.
! This (required interface) is unsupported in CLM and any attempt
! to use it will cause DART to terminate.
!
! So far, we have generated intial ensembles by taking a single
! instance and replicating it N times - and pairing each of the
! instances with a unique atmospheric forcing file and integrating
! for some period of time till the ensemble demonstrates sufficient
! spread. This is an area of active research.
!
! The naieve approach does not work -- it generates negative
! snow cover fractions, for example.  Must check for out-of-range
! values specific to each type.
! The WRF model mod has something that might be useful.

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: ens_size
real(r8),            intent(in)    :: pert_amp
logical,             intent(out)   :: interf_provided


if ( .not. module_initialized ) call static_init_model

call error_handler(E_ERR,'pert_model_copies', &
                  'CLM cannot be started from a single vector', &
                  source, revision, revdate, &
                  text2='see comments in clm/model_mod.f90::pert_model_copies()')

interf_provided = .true.

end subroutine pert_model_copies

!--------------------------------------------------------------------
!> pass the vertical localization coordinate to assim_tools_mod
function query_vert_localization_coord()

integer :: query_vert_localization_coord

query_vert_localization_coord = 1 ! any old value

end function query_vert_localization_coord

!--------------------------------------------------------------------
!> This is used in the filter_assim. The vertical conversion is done using the 
!> mean state.
subroutine vert_convert(state_handle, location, obs_kind, istatus)

type(ensemble_type), intent(in)  :: state_handle
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_kind
integer,             intent(out) :: istatus

istatus = 0

end subroutine vert_convert


!> @TODO dist should be optional ... should state_handle (intent in) be before the intent(out) 
!> @TODO state_handle is only needed for vertical conversion,  should not be here at all ...

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, &
                         obs_loc, obs_kind, num_close, close_ind, dist, state_handle)
!
! Given a DART location (referred to as "base") and a set of candidate
! locations & kinds (obs, obs_kind), returns the subset close to the
! "base", their indices, and their distances to the "base" ...

! For vertical distance computations, general philosophy is to convert all
! vertical coordinates to a common coordinate. This coordinate type can be
! defined in the namelist with the variable "vert_localization_coord".
! But we first try a single coordinate type as the model level here.
! FIXME: We need to add more options later.

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
type(ensemble_type),               intent(in)    :: state_handle

! Initialize variables to missing status

num_close = 0
close_ind = -99
dist      = 1.0e9_r8   !something big and positive (far away) in radians

! coordinate information yet (for obs_loc).
call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                          num_close, close_ind, dist)

end subroutine get_close_obs


!--------------------------------------------------------------------
! New RMA rountines
!--------------------------------------------------------------------
!> construct restart file name for reading
!> model time for CESM format?
function construct_file_name_in(stub, domain, copy)

character(len=512), intent(in) :: stub
integer,            intent(in) :: domain
integer,            intent(in) :: copy
character(len=1024)            :: construct_file_name_in

write(construct_file_name_in, '(A, i4.4)') TRIM(stub), copy

end function construct_file_name_in



function read_model_time(filename)
!------------------------------------------------------------------
! the static_init_model ensures that the clm namelists are read.
!
type(time_type) :: read_model_time
character(len=*), intent(in) :: filename

integer         :: ncid

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'read_model_time',string1,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'read_model_time', 'open '//trim(filename))

read_model_time = get_state_time_ncid(ncid)

call nc_check(nf90_close(ncid),'read_model_time', 'close '//trim(filename))

end function read_model_time


!==================================================================
! The remaining PUBLIC interfaces come next
!==================================================================

subroutine get_gridsize(num_lon, num_lat, num_lev)
integer, intent(out) :: num_lon, num_lat, num_lev
!------------------------------------------------------------------
! public utility routine.

if ( .not. module_initialized ) call static_init_model

 num_lon = nlon
 num_lat = nlat
 num_lev = nlevtot

end subroutine get_gridsize



subroutine clm_to_dart_state_vector(state_vector, restart_time)
!------------------------------------------------------------------
! Reads the current time and state variables from a clm restart
! file and packs them into a dart state vector. This better happen
! in the same fashion as the metadata arrays are built.

real(r8),         intent(out) :: state_vector(:)
type(time_type),  intent(out) :: restart_time

! temp space to hold data while we are reading it
integer  :: i, j, ni, nj, ivar, indx, numsnowlevels
integer,  allocatable, dimension(:)         :: snlsno
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array
real(r8), allocatable, dimension(:,:,:)     :: data_3d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: varname
integer :: io, TimeDimID, VarID, ncNdims, dimlen
integer :: ncid
character(len=256) :: myerrorstring

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

! Must check anything with a dimension of 'levtot' or 'levsno' and manually
! set the values to DART missing. If only it were that easy ...
!
! The treatment of snow-related variables is complicated.
! The SNLSNO variable defines the number of snow layers with valid values.
! HOWEVER, if the snow depth is < 0.01 m, the snow is not represented by a layer,
! so the SNLSNO(i) is zero even though there is a trace of snow.
! Even a trace amount of snow results in some sort of snow cover fraction.
!
! Lakes are treated differently.
! The SNLSNO(i) is always zero, even though there is snow.
! The snow over lakes is wholly contained in the bulk formulation variables
! as opposed to the snow layer variables.

! Some things must be gotten explicitly from the restart file - ONCE.
! number of snow layers
! time of restart file

allocate(snlsno(ncolumn))
call nc_check(nf90_open(trim(clm_restart_filename), NF90_NOWRITE, ncid), &
              'clm_to_dart_state_vector', 'open SNLSNO'//clm_restart_filename)
call nc_check(nf90_inq_varid(ncid,'SNLSNO', VarID), &
              'clm_to_dart_state_vector', 'inq_varid SNLSNO'//clm_restart_filename)
call nc_check(nf90_get_var(ncid, VarID, snlsno), &
              'clm_to_dart_state_vector', 'get_var SNLSNO'//clm_restart_filename)

restart_time = get_state_time(ncid)

if (do_output()) call print_time(restart_time,'time in restart file '//clm_restart_filename)
if (do_output()) call print_date(restart_time,'date in restart file '//clm_restart_filename)

call nc_check(nf90_close(ncid),'clm_to_dart_state_vector','close '//clm_restart_filename)

! Start counting and filling the state vector one item at a time,
! repacking the Nd arrays into a single 1d list of numbers.

do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   myerrorstring = trim(progvar(ivar)%origin)//' '//trim(progvar(ivar)%varname)

   call nc_check(nf90_open(trim(progvar(ivar)%origin), NF90_NOWRITE, ncid), &
              'clm_to_dart_state_vector','open '//trim(myerrorstring))

   ! File is not required to have a time dimension
   io = nf90_inq_dimid(ncid, 'time', TimeDimID)
   if (io /= NF90_NOERR) TimeDimID = MISSING_I

   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            'clm_to_dart_state_vector', 'inq_varid '//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
            'clm_to_dart_state_vector', 'inquire '//trim(myerrorstring))

   ! Check the rank of the variable

   if ( ncNdims /= progvar(ivar)%numdims ) then
      write(string1, *) 'netCDF rank of '//trim(varname)//' does not match derived type knowledge'
      write(string2, *) 'netCDF rank is ',ncNdims,' expected ',progvar(ivar)%numdims
      call error_handler(E_ERR,'clm_to_dart_state_vector', string1, &
                        source,revision,revdate,text2=string2)
   endif

   ! Check the shape of the variable

   do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'clm_to_dart_state_vector', string1)

      ! Time dimension will be 1 in progvar, but not necessarily
      ! in origin file. We only want a single matching time.
      ! static_init_model() only reserves space for a single time.
      
      if ( dimIDs(i) == TimeDimID ) dimlen = 1
          
      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(myerrorstring),' dim/dimlen ',i,dimlen, &
                              ' not ',progvar(ivar)%dimlens(i)
         call error_handler(E_ERR,'clm_to_dart_state_vector',string1,source,revision,revdate)
      endif

   enddo

   ! Pack the variable into the DART state vector
   ! Could/should fill metadata arrays at the same time ...
   ! As of 24 Aug 2011, CLM was not consistent about using a single fill_value
   ! or missing value, and the restart files didn't use the right attributes anyway ...
   ! (bugzilla report 1401)

   indx = progvar(ivar)%index1

   if (ncNdims == 1) then

      ni = progvar(ivar)%dimlens(1)
      allocate(data_1d_array(ni))
      call DART_get_var(ncid, varname, data_1d_array)

      do i = 1, ni
         state_vector(indx) = data_1d_array(i)
         indx = indx + 1
      enddo
      deallocate(data_1d_array)

   elseif (ncNdims == 2) then

      ni = progvar(ivar)%dimlens(1)
      nj = progvar(ivar)%dimlens(2)
      allocate(data_2d_array(ni, nj))
      call DART_get_var(ncid, varname, data_2d_array)

      ! README: The values in unused snow layers must be assumed to be
      ! indeterminate. If the layer is not in use, fill with a missing value.
      ! (levsno,column) and (levtot,column) variables may be treated identically.
      ! abs(snlsno(j)) defines the number of valid levels in each column -
      ! even over lakes. Lakes use a 'bulk' formula, so all the pertinent
      ! values are in the 1D variables, SNOWDP and frac_sno.

      ! FIXME: Question, what happens to unused levels below ground? Are those
      ! values 'special'?

      if     ( (trim(progvar(ivar)%dimnames(1)) == 'levsno')   .and. &
               (trim(progvar(ivar)%dimnames(2)) == 'column') ) then

         do j = 1, nj  ! loop over columns
            numsnowlevels = abs(snlsno(j))
            do i = 1, nlevsno - numsnowlevels  ! loop over layers
               data_2d_array(i,j) = MISSING_R8
            enddo
         enddo

      elseif ( (trim(progvar(ivar)%dimnames(1)) == 'levtot') .and. &
               (trim(progvar(ivar)%dimnames(2)) == 'column') ) then

         do j = 1, nj  ! loop over columns
            numsnowlevels = abs(snlsno(j))
            do i = 1, nlevsno - numsnowlevels  ! loop over layers
               data_2d_array(i,j) = MISSING_R8
            enddo
         enddo

      endif

      ! Block of checks that will hopefully be corrected in the
      ! core CLM code. There are some indeterminate values being
      ! used instead of the missing_value code - and even then,
      ! the missing_value code is not reliably implemented.

      if (progvar(ivar)%varname == 'T_SOISNO') then
         where(data_2d_array < 1.0_r8) data_2d_array = MISSING_R8
         do j = 1,nj  ! T_SOISNO has missing data in lake columns
           if (cols1d_ityplun(j) == LAKE) then
           !  write(*,*)'Found a lake column resetting the following:'
           !  write(*,*)data_2d_array(:,j)
              data_2d_array(:,j) = MISSING_R8
           endif
         enddo
      endif
      if ((progvar(ivar)%varname == 'H2OSOI_LIQ')  .or. &
          (progvar(ivar)%varname == 'H2OSOI_ICE')) then
         where(data_2d_array < 0.0_r8) data_2d_array = MISSING_R8
         do j = 1,nj  ! missing data in lake columns
           if (cols1d_ityplun(j) == LAKE) then
              data_2d_array(:,j) = MISSING_R8
           endif
         enddo
      endif

      ! Finally, pack it into the DART state vector.

      do j = 1, nj
      do i = 1, ni
         state_vector(indx) = data_2d_array(i, j)
         indx = indx + 1
      enddo
      enddo
      deallocate(data_2d_array)

   elseif (ncNdims == 3) then

      ! restart file variables never have 3 dimensions
      ! vector_history variables  may have 3 dimensions [time, lat, lon]
      ! history file variables always have 3 dimensions [time, lat, lon]
      !     exception is float H2OSOI(time, levgrnd, lat, lon) ... but we
      !     have access to restart file h2osoi_[liq,ice]

      if     ( (trim(progvar(ivar)%dimnames(1)) == 'lon')   .and. &
               (trim(progvar(ivar)%dimnames(2)) == 'lat')   .and. &
               (trim(progvar(ivar)%dimnames(3)) == 'time') ) then

         ni = progvar(ivar)%dimlens(1)
         nj = progvar(ivar)%dimlens(2)
       ! nk = progvar(ivar)%dimlens(3) not needed ... time is always a singleton

         allocate(data_3d_array(ni, nj, 1))
         call DART_get_var(ncid, varname, data_3d_array)

         ! In the CLM history files, the _missing_value_ flag seems to be
         ! applied correctly for PBOT, TBOT ... so there is no need for the
         ! extra processing that is present in the previous loops.

         do j = 1, nj
         do i = 1, ni
            state_vector(indx) = data_3d_array(i, j, 1)
            indx = indx + 1
         enddo
         enddo
         deallocate(data_3d_array)
      else

         write(string1, *) '3D variable unexpected shape -- only support nlon, nlat, time(=1)'
         write(string2, *) 'variable [',trim(progvar(ivar)%varname),']'
         write(string3, *) 'file [',trim(progvar(ivar)%origin),']'
         call error_handler(E_ERR,'clm_to_dart_state_vector', string1, &
                           source, revision, revdate, text2=string2, text3=string3)

      endif

   else

      write(string1, *) 'no support for data array of dimension ', ncNdims
      write(string2, *) 'variable [',trim(progvar(ivar)%varname),']'
      write(string3, *) 'file [',trim(progvar(ivar)%origin),']'
      call error_handler(E_ERR,'clm_to_dart_state_vector', string1, &
                        source, revision, revdate, text2=string2, text3=string3)
   endif

   indx = indx - 1
   if ( indx /= progvar(ivar)%indexN ) then
      write(string1, *)'Variable '//trim(varname)//' filled wrong.'
      write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',indx
      call error_handler(E_ERR,'clm_to_dart_state_vector', string1, &
                        source,revision,revdate,text2=string2)
   endif

   call nc_check(nf90_close(ncid),'clm_to_dart_state_vector','close '//progvar(ivar)%origin)
   ncid = 0

enddo

deallocate(snlsno)

end subroutine clm_to_dart_state_vector



subroutine sv_to_restart_file(state_vector, filename, dart_time)
!------------------------------------------------------------------
! Writes the current time and state variables from a dart state
! vector (1d array) into a clm netcdf restart file.
!
real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: filename
type(time_type),  intent(in) :: dart_time

! temp space to hold data while we are writing it
integer :: i, ni, nj, ivar
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname
integer         :: VarID, ncNdims, dimlen
integer         :: ncFileID
type(time_type) :: file_time

if ( .not. module_initialized ) call static_init_model

! Check that the output file exists ...

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for writing.'
   call error_handler(E_ERR,'sv_to_restart_file',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_WRITE, ncFileID), &
             'sv_to_restart_file','open '//trim(filename))

! make sure the time in the file is the same as the time on the data
! we are trying to insert.  we are only updating part of the contents
! of the clm restart file, and state vector contents from a different
! time won't be consistent with the rest of the file.

file_time = get_state_time(ncFileID)

if ( file_time /= dart_time ) then
   call print_time(dart_time,'DART current time',logfileunit)
   call print_time(file_time,'clm  current time',logfileunit)
   call print_time(dart_time,'DART current time')
   call print_time(file_time,'clm  current time')
   write(string1,*)trim(filename),' current time /= model time. FATAL error.'
   call error_handler(E_ERR,'sv_to_restart_file',string1,source,revision,revdate)
endif

if (do_output()) call print_time(file_time,'time of restart file '//trim(filename))
if (do_output()) call print_date(file_time,'date of restart file '//trim(filename))

! The DART prognostic variables are only defined for a single time.
! We already checked the assumption that variables are xy2d or xyz3d ...
! IF the netCDF variable has a TIME dimension, it must be the last dimension,
! and we need to read the LAST timestep and effectively squeeze out the
! singleton dimension when we stuff it into the DART state vector.

! The snow water equivalent (H2OSNO) cannot be zero since H2OSNO is used to calculate the
! bulk snow density, which in turn is a parameter in the equation of Snow Cover Fraction.
! In order to avoid the negative values of H2OSNO produced by DART, I added some "if" conditions
! to set the value of H2OSNO back to the value before assimilation if negative value is found.

UPDATE : do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   string2 = trim(filename)//' '//trim(varname)

   if ( .not. progvar(ivar)%update ) then
      write(string1,*)'intentionally not updating '//trim(string2) 
      write(string3,*)'as per namelist control in model_nml:clm_variables'
      call error_handler(E_MSG, 'sv_to_restart_file', string1, text2=string3)
      cycle UPDATE
   endif

   ! Ensure netCDF variable is conformable with progvar quantity.
   ! The TIME and Copy dimensions are intentionally not queried
   ! by looping over the dimensions stored in the progvar type.

   call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'sv_to_restart_file', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
            'sv_to_restart_file', 'inquire '//trim(string2))

   DimCheck : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
            'sv_to_restart_file', string1)

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string2),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         write(string2,*)' but it should be.'
         call error_handler(E_ERR, 'sv_to_restart_file', string1, &
                         source, revision, revdate, text2=string2)
      endif

   enddo DimCheck

   ! When called with a 4th argument, vector_to_prog_var() replaces the DART
   ! missing code with the value in the corresponding variable in the netCDF file.
   ! Any clamping to physically meaningful values occurrs in vector_to_prog_var.

   if (progvar(ivar)%numdims == 1) then

      ni = progvar(ivar)%dimlens(1)
      allocate(data_1d_array(ni))
      call vector_to_prog_var(state_vector, ivar, data_1d_array, ncFileID)

      call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array), &
            'sv_to_restart_file', 'put_var '//trim(varname))
      deallocate(data_1d_array)

   elseif (progvar(ivar)%numdims == 2) then

      ni = progvar(ivar)%dimlens(1)
      nj = progvar(ivar)%dimlens(2)
      allocate(data_2d_array(ni, nj))
      call vector_to_prog_var(state_vector, ivar, data_2d_array, ncFileID)

      call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array), &
            'sv_to_restart_file', 'put_var '//trim(varname))
      deallocate(data_2d_array)

   else
      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'sv_to_restart_file', string1, &
                        source,revision,revdate)
   endif

   ! TJH FIXME ... this works perfectly if it were not for a bug in netCDF.
   ! When they fix the bug, this will be a useful thing to restore.
   ! Make note that the variable has been updated by DART
!  call nc_check(nf90_Redef(ncFileID),'sv_to_restart_file', 'redef '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID, VarID,'DART','variable modified by DART'),&
!                'sv_to_restart_file', 'modified '//trim(varname))
!  call nc_check(nf90_enddef(ncfileID),'sv_to_restart_file','state enddef '//trim(filename))

enddo UPDATE

call nc_check(nf90_close(ncFileID),'sv_to_restart_file','close '//trim(filename))
ncFileID = 0

end subroutine sv_to_restart_file



!==================================================================
! The remaining interfaces come last
!==================================================================


subroutine model_interpolate(state_handle, ens_size, location, obs_kind, expected_obs, istatus)


! PURPOSE:
!
! For a given lat, lon, and height, interpolate the correct state value
! to that location for the filter from the clm state vectors
!
! Reconstructing the vertical profile of the gridcell is complicated.
! Each land unit/column can have a different number of vertical levels.
! Do I just try to recontruct the whole profile and then interpolate?
! Impossible to know which elements are 'above' and 'below' without
! finding all the elements in the first place. The vertical information
! is in the levels() array for each state vector component.

! Passed variables

! Interpolates from state vector x to the location. It's not particularly
! happy dumping all of this straight into the model. Eventually some
! concept of a grid underlying models but above locations is going to
! be more general. May want to wait on external infrastructure projects
! for this?

! Argument itype is not used here because there is only one type of variable.
! Type is needed to allow swap consistency with more complex models.

type(ensemble_type),    intent(in)  :: state_handle
integer,                intent(in)  :: ens_size
type(location_type),    intent(in)  :: location
integer,                intent(in)  :: obs_kind
real(r8),               intent(out) :: expected_obs(ens_size)
integer,                intent(out) :: istatus(ens_size)


! Local storage

real(r8), dimension(LocationDims) :: loc_array
real(r8) :: llon, llat, lheight
integer  :: imem
integer  :: istatus_2(ens_size)
real(r8) :: interp_val_2(ens_size)
character(len=paramname_length) :: kind_string

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

expected_obs = MISSING_R8     ! the DART bad value flag
interp_val_2 = MISSING_R8     ! the DART bad value flag
istatus = 99                ! unknown error

! Get the individual locations values

loc_array = get_location(location)
llon      = loc_array(1)
llat      = loc_array(2)
lheight   = loc_array(3)

if ((debug > 6) .and. do_output()) print *, 'requesting interpolation at ', llon, llat, lheight

! Standard method of interpolation.
! get_grid_vertval()       for quantities that have a vertical profile.
! compute_gridcell_value() for quantities computed for the entire gridcell.

select case( obs_kind )

   case ( KIND_SOIL_MOISTURE )

      ! TJH FIXME - actually ROLAND FIXME
      ! This is terrible ... the COSMOS operator wants m3/m3 ... CLM is kg/m2
      call get_grid_vertval(state_handle, ens_size, location, KIND_LIQUID_WATER, expected_obs,  istatus)
      call get_grid_vertval(state_handle, ens_size, location, KIND_ICE,          interp_val_2, istatus_2)

      do imem = 1,ens_size
         if ((istatus(imem) == 0) .and. (istatus_2(imem) == 0)) then
            expected_obs(imem) = expected_obs(imem) + interp_val_2(imem)
         else
            expected_obs(imem) = MISSING_R8
            istatus(imem)      = 6
         endif
      enddo

   case ( KIND_SOIL_TEMPERATURE, KIND_LIQUID_WATER, KIND_ICE )

      call get_grid_vertval(state_handle, ens_size, location, obs_kind, expected_obs, istatus)

   case ( KIND_SNOWCOVER_FRAC, KIND_LEAF_AREA_INDEX, KIND_LEAF_CARBON, KIND_WATER_TABLE_DEPTH, &
          KIND_VEGETATION_TEMPERATURE)

      call compute_gridcell_value(state_handle, ens_size, location, obs_kind, expected_obs, istatus)

   case default

      kind_string = get_raw_obs_kind_name(obs_kind)

      write(string1,*)'not written for (integer) kind ',obs_kind
      write(string2,*)'AKA '//trim(kind_string)
      call error_handler(E_ERR, 'model_interpolate', string1, &
             source, revision, revdate, text2=string2)
      istatus = 5

end select

if ((debug > 6) .and. do_output()) write(*,*)'expected_obs ',expected_obs

end subroutine model_interpolate


!------------------------------------------------------------------


subroutine compute_gridcell_value(state_handle, ens_size, location, kind_index, interp_val, istatus)
!
! Each gridcell may contain values for several land units, each land unit may contain
! several columns, each column may contain several pft's. BUT this routine never
! aggregates across multiple pft's. So, each gridcell value
! is an area-weighted value of an unknown number of column-based quantities.

! Passed variables

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location     ! location somewhere in a grid cell
integer,             intent(in)  :: kind_index   ! KIND in DART state needed for interpolation 
real(r8),            intent(out) :: interp_val(ens_size)   ! area-weighted result
integer,             intent(out) :: istatus(ens_size)      ! error code (0 == good)

! Local storage

integer  :: ivar, index1, indexN, counter(ens_size)
integer(i8) :: indexi
integer  :: gridloni,gridlatj
real(r8) :: loc_lat, loc_lon
real(r8) :: state(ens_size)
real(r8) :: total(ens_size)
real(r8) :: total_area(ens_size)
real(r8), dimension(1) :: loninds,latinds
real(r8), dimension(LocationDims) :: loc
integer :: imem
character(len=paramname_length) :: varstring

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

interp_val = MISSING_R8  ! the DART bad value flag
istatus    = 99          ! unknown error

loc        = get_location(location)  ! loc is in DEGREES
loc_lon    = loc(1)
loc_lat    = loc(2)

! determine the portion of interest of the state vector
ivar   = findKindIndex(kind_index, 'compute_gridcell_value')
index1 = progvar(ivar)%index1 ! in the DART state vector, start looking here
indexN = progvar(ivar)%indexN ! in the DART state vector, stop  looking here

varstring = progvar(ivar)%varname ! used in error messages

! BOMBPROOFING - check for a vertical dimension for this variable
if (progvar(ivar)%maxlevels > 1) then
   write(string1, *)'Variable '//trim(varstring)//' cannot use this routine.'
   write(string2, *)'use get_grid_vertval() instead.'
   call error_handler(E_ERR,'compute_gridcell_value', string1, &
                  source, revision, revdate, text2=string2)
endif

! determine the grid cell for the location
latinds  = minloc(abs(LAT - loc_lat))   ! these return 'arrays' ...
loninds  = minloc(abs(LON - loc_lon))   ! these return 'arrays' ...
gridlatj = latinds(1)
gridloni = loninds(1)

if ((debug > 5) .and. do_output()) then
   write(*,*)'compute_gridcell_value:targetlon, lon, lon index is ',&
                  loc_lon,LON(gridloni),gridloni
   write(*,*)'compute_gridcell_value:targetlat, lat, lat index is ',&
                  loc_lat,LAT(gridlatj),gridlatj
endif

! If there is no vertical component, the problem is greatly simplified.
! Simply area-weight an average of all pieces in the grid cell.
! FIXME ... this is the loop that can exploit the knowledge of what 
! columnids or pftids are needed for any particular gridcell.
! gridCellInfo%pftids, gridCellInfo%columnids

counter    = 0
total      = 0.0_r8      ! temp storage for state vector
total_area = 0.0_r8      ! temp storage for area
ELEMENTS : do indexi = index1, indexN

   if (   lonixy(indexi) /=  gridloni ) cycle ELEMENTS
   if (   latjxy(indexi) /=  gridlatj ) cycle ELEMENTS
   if ( landarea(indexi) ==   0.0_r8  ) cycle ELEMENTS

   state = get_state(indexi, state_handle)

   do imem = 1, ens_size
      if(state(imem) /= MISSING_R8) then
         counter(imem)    = counter(imem)    + 1
         total(imem)      = total(imem)      + state(imem)*landarea(indexi)
         total_area(imem) = total_area(imem) +       landarea(indexi)
      endif
   enddo

   if ((debug > 5) .and. do_output()) then
      write(*,*)
      write(*,*)'gridcell location match',counter(1),'at statevector index',indexi
  !    write(*,*)'statevector value is (',state(indexi),')'
      write(*,*)'area is              (',landarea(indexi),')'
      write(*,*)'LON index is         (',lonixy(indexi),')'
      write(*,*)'LAT index is         (',latjxy(indexi),')'
      write(*,*)'closest LON is       (',LON(gridloni),')'
      write(*,*)'closest LAT is       (',LAT(gridlatj),')'
      write(*,*)'closest lev is       (',levels(indexi),')'
   endif

enddo ELEMENTS


if (all(total_area /= 0.0_r8)) then ! All good.
   interp_val = total/total_area
   istatus    = 0
else
   if ((debug > 4) .and. do_output()) then
      write(string1, *)'Variable '//trim(varstring)//' had no viable data'
      write(string2, *)'at gridcell ilon/jlat = (',gridloni,',',gridlatj,')'
      write(string3, *)'obs lon/lat = (',loc_lon,',',loc_lat,')'
      call error_handler(E_MSG,'compute_gridcell_value', string1, &
                     text2=string2,text3=string3)
   endif
endif

!> @todo FIXME Need to print debugging info for any task, not just task 0
! Print more information for the really curious
if ((debug > 5)) then
   ! write(string1,*)'counter, total, total_area', counter, total, total_area
   write(string1,*)'counter, total_area', counter(1), total_area(1)
   write(string2,*)'interp_val, istatus', interp_val(1), istatus(1)
   call error_handler(E_MSG,'compute_gridcell_value', string1, text2=string2)
endif

end subroutine compute_gridcell_value


!------------------------------------------------------------------


subroutine get_grid_vertval(state_handle, ens_size, location, kind_index, interp_val, istatus)

! Calculate the expected vertical value for the gridcell.
! Each gridcell value is an area-weighted value of an unknown number of
! column-based quantities.

! Passed variables

type(ensemble_type), intent(in)  :: state_handle ! state vector
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location     ! location somewhere in a grid cell
integer,             intent(in)  :: kind_index   ! KIND in DART state needed for interpolation
real(r8),            intent(out) :: interp_val(ens_size)   ! area-weighted result
integer,             intent(out) :: istatus(ens_size)      ! error code (0 == good)

! Local storage

integer  :: ivar, index1, indexN, counter1, counter2
integer(i8) :: indexi
integer  :: gridloni,gridlatj
real(r8), dimension(LocationDims) :: loc
real(r8) :: loc_lat, loc_lon, loc_lev
real(r8) :: value_below(ens_size), value_above(ens_size), total_area(ens_size)
real(r8) :: depthbelow, depthabove
real(r8) :: topwght, botwght
real(r8), dimension(1) :: loninds,latinds

real(r8), allocatable, dimension(:, :) :: above, below
real(r8), allocatable, dimension(:, :) :: area_above
real(r8), allocatable, dimension(:, :) :: area_below
integer :: counter, counter_above, counter_below
integer :: imem
real(r8) :: state(ens_size)
character(len=paramname_length) :: varstring

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

interp_val = MISSING_R8  ! the DART bad value flag
istatus    = 99          ! unknown error

loc        = get_location(location)  ! loc is in DEGREES
loc_lon    = loc(1)
loc_lat    = loc(2)
loc_lev    = loc(3)

if ( loc_lev < 0.0_r8 ) then
   write(string1,*)'Cannot support above-ground vertical interpolation.'
   write(string2,*)'requested a value at a depth of ',loc_lev
   write(string3,*)'CLM has negative depths to indicate above-ground values.'
   call error_handler(E_ERR,'get_grid_vertval', string1, &
      source, revision, revdate, text2=string2, text3=string3)
endif

! determine the portion of interest of the state vector
ivar   = findKindIndex(kind_index, 'get_grid_vertval')
index1 = progvar(ivar)%index1 ! in the DART state vector, start looking here
indexN = progvar(ivar)%indexN ! in the DART state vector, stop  looking here

varstring = progvar(ivar)%varname ! used in a lot of error messages

! BOMBPROOFING - check for a vertical dimension for this variable
if (progvar(ivar)%maxlevels < 2) then
   write(string1, *)'Variable '//trim(varstring)//' should not use this routine.'
   write(string2, *)'use compute_gridcell_value() instead.'
   call error_handler(E_ERR,'get_grid_vertval', string1, &
                  source, revision, revdate, text2=string2)
endif

! determine the grid cell for the location
latinds  = minloc(abs(LAT - loc_lat))   ! these return 'arrays' ...
loninds  = minloc(abs(LON - loc_lon))   ! these return 'arrays' ...
gridlatj = latinds(1)
gridloni = loninds(1)

if ((debug > 4) .and. do_output()) then
   write(*,*)'get_grid_vertval:targetlon, lon, lon index, level is ', &
              loc_lon,LON(gridloni),gridloni,loc_lev
   write(*,*)'get_grid_vertval:targetlat, lat, lat index, level is ', &
              loc_lat,LAT(gridlatj),gridlatj,loc_lev
endif

! Determine the level 'above' and 'below' the desired vertical
! The above-ground 'depths' are calculated from ZISNO and are negative.
! The 'depths' are all positive numbers, increasingly positive is deeper.
! The variables currently supported use the subsurface definitions in
! the module variable LEVNGRND.

if (loc_lev  <= LEVGRND(1)) then  ! the top level is so close to the surface
   depthabove = LEVGRND(1)        ! just use the top level
   depthbelow = LEVGRND(1)
elseif (loc_lev >= maxval(LEVGRND)) then  ! at depth, however ... do we
   depthabove    = maxval(LEVGRND)        ! fail or just use the deepest
   depthbelow    = maxval(LEVGRND)        ! I am using the deepest.
else

   LAYERS : do indexi = 2,size(LEVGRND)
      if (loc_lev < LEVGRND(indexi)) then
         depthabove = LEVGRND(indexi-1)
         depthbelow = LEVGRND(indexi  )
         exit LAYERS
      endif
   enddo LAYERS

endif

if ((debug > 4) .and. do_output()) then
   write(*,*)'get_grid_vertval:depthbelow ',depthbelow,'>= loc_lev', &
                   loc_lev,'>= depthabove',depthabove
endif

! Determine how many elements can contribute to the gridcell value.
! There are multiple column-based contributors, each column has a
! separate area-based weight. There are multiple levels.
! I believe I have to keep track of all of them to sort out how to
! calculate the gridcell value at a particular depth.

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

if ( counter == 0 ) then
   if ((debug > 0) .and. do_output()) then
      write(string1, *)'statevector variable '//trim(varstring)//' had no viable data'
      write(string2, *)'at gridcell lon/lat = (',gridloni,',',gridlatj,')'
      write(string3, *)'obs lon/lat/lev (',loc_lon,',',loc_lat,',',loc_lev,')'
      call error_handler(E_MSG,'get_grid_vertval', string1, &
                     text2=string2,text3=string3)
   endif
   return
endif


allocate(above(ens_size, counter),&
         below(ens_size, counter),&
         area_below(ens_size, counter), &
         area_above(ens_size, counter))

above(:, :)  = 0.0_r8
below(:, :)  = 0.0_r8
area_below(:, :) = 0.0_r8
area_above(:, :) = 0.0_r8

counter_above = 0
counter_below = 0
ELEMENTS : do indexi = index1, indexN

   if ( lonixy(indexi) /=  gridloni )  cycle ELEMENTS
   if ( latjxy(indexi) /=  gridlatj )  cycle ELEMENTS

!  write(*,*)'level ',indexi,' is ',levels(indexi),' location depth is ',loc_lev
   state = get_state(indexi, state_handle)

   if (levels(indexi) == depthabove) then

      counter_above           = counter_above + 1
      above(:, counter_above) = state
      where(state /= MISSING_R8) area_above(:, counter_above) = landarea(indexi)

   elseif(levels(indexi) == depthbelow) then

      counter_below            = counter_below + 1
      below(:, counter_below)    =  state
      where(state /= MISSING_R8) area_below(:, counter_below) = landarea(indexi)

   else
      cycle ELEMENTS
   endif

   if ((debug > 4) .and. do_output()) then
   write(*,*)
   write(*,*)'gridcell location match at statevector index',indexi
   !write(*,*)'statevector value is (',state,')'
   write(*,*)'area is          (',landarea(indexi),')'
   write(*,*)'LON index is     (',lonixy(indexi),')'
   write(*,*)'LAT index is     (',latjxy(indexi),')'
   write(*,*)'gridcell LON is  (',LON(gridloni),')'
   write(*,*)'gridcell LAT is  (',LAT(gridlatj),')'
   write(*,*)'depth        is  (',levels(indexi),')'
   endif

enddo ELEMENTS

! could arise if the above or below was 'missing' ... but the mate was not.

if ( counter_above /= counter_below ) then
   write(string1, *)'Variable '//trim(varstring)//' has peculiar interpolation problems.'
   write(string2, *)'uneven number of values "above" and "below"'
   write(string3, *)'counter_above == ',counter_above,' /= ',counter_below,' == counter_below'
   call error_handler(E_MSG,'get_grid_vertval', string1, &
                  text2=string2,text3=string3)
   return
endif

do imem = 1, ens_size

   ! Determine the value for the level above the depth of interest.
   total_area(imem) = sum(area_above(imem, :))

   if ( total_area(imem) /= 0.0_r8 ) then
      ! normalize the area-based weights
      area_above(imem, :) = area_above(imem, :) / total_area(imem)
      value_above(imem) = sum(above(imem, :) * area_above(imem, :))
   else
      write(string1, *)'Variable '//trim(varstring)//' had no viable data above'
      write(string2, *)'at gridcell lon/lat/lev = (',gridloni,',',gridlatj,',',depthabove,')'
      write(string3, *)'obs lon/lat/lev (',loc_lon,',',loc_lat,',',loc_lev,')'
      call error_handler(E_ERR,'get_grid_vertval', string1, &
                  source, revision, revdate, text2=string2,text3=string3)
   endif

   ! Determine the value for the level below the depth of interest.
   total_area(imem) = sum(area_below(imem, :))

   if ( total_area(imem) /= 0.0_r8 ) then
      ! normalize the area-based weights
      area_below(imem, :) = area_below(imem, :) / total_area(imem)
      value_below(imem) = sum(below(imem, :) * area_below(imem, :))
   else
      write(string1, *)'Variable '//trim(varstring)//' had no viable data below'
      write(string2, *)'at gridcell lon/lat/lev = (',gridloni,',',gridlatj,',',depthbelow,')'
      write(string3, *)'obs lon/lat/lev (',loc_lon,',',loc_lat,',',loc_lev,')'
      call error_handler(E_ERR,'get_grid_vertval', string1, &
                  source, revision, revdate, text2=string2,text3=string3)
   endif

enddo

if (depthbelow == depthabove) then
   topwght = 1.0_r8
   botwght = 0.0_r8
else
   topwght = (depthbelow - loc_lev) / (depthbelow - depthabove)
   botwght = (loc_lev - depthabove) / (depthbelow - depthabove)
endif

interp_val = value_above*topwght + value_below*botwght
istatus      = 0

deallocate(above, below, area_above, area_below)

end subroutine get_grid_vertval


!------------------------------------------------------------------


subroutine vector_to_1d_prog_var(x, ivar, data_1d_array, ncid)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset, into a 1d array.
!
! If the optional argument (ncid) is specified, some additional
! processing takes place. The variable in the netcdf is read.
! This must be the same shape as the intended output array.
! Anywhere the DART MISSING code is encountered in the input array,
! the corresponding (i.e. original) value from the netCDF file is
! used.

real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:),   intent(out) :: data_1d_array
integer, OPTIONAL,        intent(in)  :: ncid

integer :: i,ii, VarID
real(r8), allocatable, dimension(:) :: org_array

! unpack the right part of the DART state vector into a 1D array.

ii = progvar(ivar)%index1

do i = 1, progvar(ivar)%dimlens(1)
   data_1d_array(i) = x(ii)
   ii = ii + 1
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_1d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

! Apply the min/max values, if applicable
! This should only be true when converting to a variable that will
! be reinserted into the CLM restart file. This is indicated
! by the presence of the ncid variable.

if (present(ncid)) then

   if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_1d_array /= MISSING_R8) .and. &
             (data_1d_array > progvar(ivar)%maxvalue)) &
              data_1d_array = progvar(ivar)%maxvalue
   endif

   if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_1d_array /= MISSING_R8) .and. &
             (data_1d_array < progvar(ivar)%minvalue)) &
              data_1d_array = progvar(ivar)%minvalue
   endif

   ! Replace the DART fill value with the original value and apply any clamping.
   ! Get the 'original' variable from the netcdf file.

   allocate(org_array(size(data_1d_array)))

   call nc_check(nf90_inq_varid(ncid, progvar(ivar)%varname, VarID), &
            'vector_to_1d_prog_var', 'inq_varid '//trim(progvar(ivar)%varname))

   call nc_check(nf90_get_var(ncid, VarID, org_array), &
            'vector_to_1d_prog_var', 'get_var '//trim(progvar(ivar)%varname))

   ! restoring the indeterminate original values

   where(data_1d_array == MISSING_R8) data_1d_array = org_array

   ! clamping the assimilated values to physically meaningful ranges.

   if (trim(progvar(ivar)%varname) == 'SNOWDP') &
      where((data_1d_array < 0.0_r8)) data_1d_array = org_array

   if (trim(progvar(ivar)%varname) == 'H2OSNO') &
      where((data_1d_array <= 0.0_r8)) data_1d_array = org_array

   deallocate(org_array)

else

   if     (progvar(ivar)%xtype == NF90_INT) then
      where(data_1d_array == MISSING_I) data_1d_array = progvar(ivar)%spvalINT
   elseif (progvar(ivar)%xtype == NF90_FLOAT) then
      where(data_1d_array == MISSING_R4) data_1d_array = progvar(ivar)%spvalR4
   elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
      where(data_1d_array == MISSING_R8) data_1d_array = progvar(ivar)%spvalR8
   endif

endif

end subroutine vector_to_1d_prog_var


!------------------------------------------------------------------


subroutine vector_to_2d_prog_var(x, ivar, data_2d_array, ncid)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset,
! into a 2d array.
!
real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:,:), intent(out) :: data_2d_array
integer, OPTIONAL,        intent(in)  :: ncid

integer :: i,j,ii, VarID
real(r8), allocatable, dimension(:,:) :: org_array

! unpack the right part of the DART state vector into a 1D array.

ii = progvar(ivar)%index1

do j = 1,progvar(ivar)%dimlens(2)
do i = 1,progvar(ivar)%dimlens(1)
   data_2d_array(i,j) = x(ii)
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

! Apply the min/max values, if applicable
! This should only be true when converting to a variable that will
! be reinserted into the CLM restart file. This is indicated
! by the presence of the ncid variable.

if (present(ncid)) then

   if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_2d_array /= MISSING_R8) .and. &
             (data_2d_array > progvar(ivar)%maxvalue)) &
              data_2d_array = progvar(ivar)%maxvalue
   endif

   if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_2d_array /= MISSING_R8) .and. &
             (data_2d_array < progvar(ivar)%minvalue)) &
              data_2d_array = progvar(ivar)%minvalue
   endif

   ! Replace the DART fill value with the original value and apply any clamping.
   ! Get the 'original' variable from the netcdf file if need be.

   allocate(org_array(size(data_2d_array,1),size(data_2d_array,2)))

   call nc_check(nf90_inq_varid(ncid, progvar(ivar)%varname, VarID), &
            'vector_to_2d_prog_var', 'inq_varid '//trim(progvar(ivar)%varname))

   call nc_check(nf90_get_var(ncid, VarID, org_array), &
            'vector_to_2d_prog_var', 'get_var '//trim(progvar(ivar)%varname))

   ! restoring the indeterminate original values

   where(data_2d_array == MISSING_R8 ) data_2d_array = org_array

   ! clamping the assimilated values to physically meaningful ranges.

   if     (trim(progvar(ivar)%varname) == 'DZSNO') then
      where((data_2d_array < 0.0_r8)) data_2d_array = org_array
   elseif (trim(progvar(ivar)%varname) == 'ZSNO') then
      where((data_2d_array < 0.0_r8)) data_2d_array = org_array
   elseif (trim(progvar(ivar)%varname) == 'ZISNO') then
      where((data_2d_array < 0.0_r8)) data_2d_array = org_array
   elseif (trim(progvar(ivar)%varname) == 'H2OSOI_LIQ') then
      where((data_2d_array < 0.0_r8)) data_2d_array = org_array
   elseif (trim(progvar(ivar)%varname) == 'H2OSOI_ICE') then
      where((data_2d_array < 0.0_r8)) data_2d_array = org_array
   elseif (trim(progvar(ivar)%varname) == 'T_SOISNO') then
      where((data_2d_array < 0.0_r8)) data_2d_array = org_array
   elseif (trim(progvar(ivar)%varname) == 'T_LAKE') then
      where((data_2d_array < 0.0_r8)) data_2d_array = org_array
   elseif (trim(progvar(ivar)%varname) == 'leafc') then
      where((data_2d_array < 0.0_r8)) data_2d_array = 0.0_r8
   endif

   deallocate(org_array)

else

   if     (progvar(ivar)%xtype == NF90_INT) then
      where(data_2d_array == MISSING_I) data_2d_array = progvar(ivar)%spvalINT
   elseif (progvar(ivar)%xtype == NF90_FLOAT) then
      where(data_2d_array == MISSING_R4) data_2d_array = progvar(ivar)%spvalR4
   elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
      where(data_2d_array == MISSING_R8) data_2d_array = progvar(ivar)%spvalR8
   endif

endif

end subroutine vector_to_2d_prog_var


!------------------------------------------------------------------


subroutine get_history_dims(ncid,fname,cstat,lon,lat,levgrnd,lonatm,latatm,lonrof,latrof)
!------------------------------------------------------------------
!
! Read the dimensions from the history netcdf file.
!
! The file name comes from module storage ... namelist.

integer,           intent(inout) :: ncid
character(len=*),  intent(in)    :: fname
character(len=*),  intent(in)    :: cstat ! how do you want to leave the netcdf file
integer,           intent(out)   :: lon
integer,           intent(out)   :: lat
integer,           intent(out)   :: levgrnd
integer, OPTIONAL, intent(out)   :: lonatm, latatm
integer, OPTIONAL, intent(out)   :: lonrof, latrof

integer :: dimid

! get the ball rolling ...

if (ncid == 0) then ! we need to open it
   call nc_check(nf90_open(trim(fname), nf90_nowrite, ncid), &
       'get_history_dims','open '//trim(fname))
endif

! The new SingleColumn (and unstructured grid) configurations
! do not have a 'lon' and 'lat' dimension. There is only 'lndgrid'

if ( nf90_inq_dimid(ncid, 'lndgrid', dimid) == NF90_NOERR ) unstructured = .true.

if (unstructured) then ! use the lndgrid dimension for both lon and lat

      call nc_check(nf90_inq_dimid(ncid, 'lndgrid', dimid), &
                  'get_history_dims','inq_dimid lndgrid '//trim(fname))
      call nc_check(nf90_inquire_dimension(ncid, dimid, len=lon), &
                  'get_history_dims','inquire_dimension lndgrid '//trim(fname))
      lat = lon

else

    call nc_check(nf90_inq_dimid(ncid, 'lon', dimid), &
                'get_history_dims','inq_dimid lon '//trim(fname))
    call nc_check(nf90_inquire_dimension(ncid, dimid, len=lon), &
                'get_history_dims','inquire_dimension lon '//trim(fname))

    call nc_check(nf90_inq_dimid(ncid, 'lat', dimid), &
                'get_history_dims','inq_dimid lat '//trim(fname))
    call nc_check(nf90_inquire_dimension(ncid, dimid, len=lat), &
                'get_history_dims','inquire_dimension lat '//trim(fname))

endif

call nc_check(nf90_inq_dimid(ncid, 'levgrnd', dimid), &
            'get_history_dims','inq_dimid levgrnd '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=levgrnd), &
            'get_history_dims','inquire_dimension levgrnd '//trim(fname))

if (present(lonatm)) then
   call nc_check(nf90_inq_dimid(ncid, 'lonatm', dimid), &
               'get_history_dims','inq_dimid lonatm '//trim(fname))
   call nc_check(nf90_inquire_dimension(ncid, dimid, len=lonatm), &
               'get_history_dims','inquire_dimension lonatm '//trim(fname))
endif

if (present(latatm)) then
   call nc_check(nf90_inq_dimid(ncid, 'latatm', dimid), &
               'get_history_dims','inq_dimid latatm '//trim(fname))
   call nc_check(nf90_inquire_dimension(ncid, dimid, len=latatm), &
               'get_history_dims','inquire_dimension latatm '//trim(fname))
endif

if (present(lonrof)) then
   call nc_check(nf90_inq_dimid(ncid, 'lonrof', dimid), &
               'get_history_dims','inq_dimid lonrof '//trim(fname))
   call nc_check(nf90_inquire_dimension(ncid, dimid, len=lonrof), &
               'get_history_dims','inquire_dimension lonrof '//trim(fname))
endif

if (present(latrof)) then
   call nc_check(nf90_inq_dimid(ncid, 'latrof', dimid), &
               'get_history_dims','inq_dimid latrof '//trim(fname))
   call nc_check(nf90_inquire_dimension(ncid, dimid, len=latrof), &
               'get_history_dims','inquire_dimension latrof '//trim(fname))
endif

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
   call nc_check(nf90_close(ncid),'get_history_dims','close '//trim(fname) )
   ncid = 0
endif

end subroutine get_history_dims


!------------------------------------------------------------------


subroutine get_full_grid(ncid, fname, cstat)
!------------------------------------------------------------------
!
! Read the grid dimensions from the CLM history netcdf file.
! LON,LAT,AREA,LANDFRAC,LEVGRN all have module scope

integer,                  intent(inout) :: ncid
character(len=*),         intent(in)    :: fname
character(len=*),         intent(in)    :: cstat

! Make sure the variables are the right size ...
! at some point in the future ...

if (ncid == 0) then ! we need to open it
   call nc_check(nf90_open(trim(fname), nf90_nowrite, ncid), &
       'get_full_grid','open '//trim(fname))
endif

! The lat/lon matrices in the history file have been masked by
! the land values such that the wet cells are 'missing' values.
! This makes it less than useful for us. Thankfully, the 1D
! lat/lon arrays have no such mask applied. We use these.

call DART_get_var(ncid,'lon'     ,LON)
call DART_get_var(ncid,'lat'     ,LAT)
call DART_get_var(ncid,'levgrnd' ,LEVGRND)
if (unstructured) then
   call DART_get_var(ncid,'area'    ,AREA1D)
   call DART_get_var(ncid,'landfrac',LANDFRAC1D)
   where(AREA1D     == MISSING_R8) AREA1D     = 0.0_r8
   where(LANDFRAC1D == MISSING_R8) LANDFRAC1D = 0.0_r8
else
   call DART_get_var(ncid,'area'    ,AREA2D)
   call DART_get_var(ncid,'landfrac',LANDFRAC2D)
   where(AREA2D     == MISSING_R8) AREA2D     = 0.0_r8
   where(LANDFRAC2D == MISSING_R8) LANDFRAC2D = 0.0_r8
endif


! just to make sure we are [0,360] and [-90,90]

where (LON <   0.0_r8) LON = LON + 360.0_r8
where (LON > 360.0_r8) LON = LON - 360.0_r8

if (any(LON < 0.0_r8)) then
   write(string1,*)'longitudes in history file variable "lon" still negative.'
   call error_handler(E_ERR,'get_full_grid',string1,source,revision,revdate)
endif

where (LAT < -90.0_r8) LAT = -90.0_r8
where (LAT >  90.0_r8) LAT =  90.0_r8

if (cstat == 'close') then
   call nc_check(nf90_close(ncid),'get_full_grid','close '//trim(fname) )
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

return
end subroutine get_full_grid


!------------------------------------------------------------------


subroutine get_sparse_dims(ncid, fname, cstat)
!------------------------------------------------------------------
!
! Read the dimensions from the CLM restart netcdf file.

integer,          intent(inout) :: ncid
character(len=*), intent(in)    :: fname
character(len=*), intent(in)    :: cstat

integer :: dimid, istatus, mylevgrnd

if (ncid == 0) then
   call nc_check(nf90_open(trim(fname), nf90_nowrite, ncid), &
               'get_sparse_dims','open '//trim(fname))
endif

! get dimid for 'gridcell' and then get value ...

call nc_check(nf90_inq_dimid(ncid, 'gridcell', dimid), &
            'get_sparse_dims','inq_dimid gridcell '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=ngridcell), &
            'get_sparse_dims','inquire_dimension gridcell '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'landunit', dimid), &
            'get_sparse_dims','inq_dimid landunit '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=nlandunit), &
            'get_sparse_dims','inquire_dimension landunit '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'column', dimid), &
            'get_sparse_dims','inq_dimid column '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=ncolumn), &
            'get_sparse_dims','inquire_dimension column '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'pft', dimid), &
            'get_sparse_dims','inq_dimid pft '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=npft), &
            'get_sparse_dims','inquire_dimension pft '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'levgrnd', dimid), &
            'get_sparse_dims','inq_dimid levgrnd '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=mylevgrnd), &
            'get_sparse_dims','inquire_dimension levgrnd '//trim(fname))

if (mylevgrnd /= nlevgrnd) then
   write(string1,*)'Number of ground levels in restart file is',mylevgrnd
   write(string2,*)'Number of ground levels in history file is', nlevgrnd
   call error_handler(E_ERR,'get_sparse_dims',string1,source,revision,revdate,text2=string2)
endif

call nc_check(nf90_inq_dimid(ncid, 'levlak', dimid), &
            'get_sparse_dims','inq_dimid levlak '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=nlevlak), &
            'get_sparse_dims','inquire_dimension levlak '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'levtot', dimid), &
            'get_sparse_dims','inq_dimid levtot '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=nlevtot), &
            'get_sparse_dims','inquire_dimension levtot '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'numrad', dimid), &
            'get_sparse_dims','inq_dimid numrad '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=nnumrad), &
            'get_sparse_dims','inquire_dimension numrad '//trim(fname))

! CLM4 does not have a multi-level canopy.
! CLM4.5 has a multi-level canopy.
istatus = nf90_inq_dimid(ncid, 'levcan', dimid)
if (istatus == nf90_noerr) then
   call nc_check(nf90_inquire_dimension(ncid, dimid, len=nlevcan), &
               'get_sparse_dims','inquire_dimension levcan '//trim(fname))
endif

! levsno is presently required, but I can envision a domain/experiment that
! will not have snow levels. How this relates to variables dimensioned 'levtot'
! is unclear. For that reason, levsno is presently required.

istatus = nf90_inq_dimid(ncid, 'levsno', dimid)
if (istatus == nf90_noerr) then
   call nc_check(nf90_inquire_dimension(ncid, dimid, len=nlevsno), &
               'get_sparse_dims','inquire_dimension levsno '//trim(fname))
endif

! levsno1, rtmlon, rtmlat are optional ... so it is not a fatal error if they are not present.

istatus = nf90_inq_dimid(ncid, 'levsno1', dimid)
if (istatus == nf90_noerr) then
   call nc_check(nf90_inquire_dimension(ncid, dimid, len=nlevsno1), &
               'get_sparse_dims','inquire_dimension levsno1 '//trim(fname))
endif

istatus = nf90_inq_dimid(ncid, 'rtmlon', dimid)
if (istatus == nf90_noerr) then
   call nc_check(nf90_inquire_dimension(ncid, dimid, len=nrtmlon), &
               'get_sparse_dims','inquire_dimension rtmlon '//trim(fname))
endif

istatus = nf90_inq_dimid(ncid, 'rtmlat', dimid)
if (istatus == nf90_noerr) then
   call nc_check(nf90_inquire_dimension(ncid, dimid, len=nrtmlat), &
               'get_sparse_dims','inquire_dimension rtmlat '//trim(fname))
endif

if (cstat == 'close') then
   call nc_check(nf90_close(ncid),'get_sparse_dims','close '//trim(fname) )
   ncid = 0
endif

! Echo what we know if desired.
if ((debug > 7) .and. do_output()) then
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
   write(logfileunit,*)'nrtmlon   = ',nrtmlon
   write(logfileunit,*)'nrtmlat   = ',nrtmlat
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
   write(     *     ,*)'nrtmlon   = ',nrtmlon
   write(     *     ,*)'nrtmlat   = ',nrtmlat
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

integer  :: VarID

if (ncid == 0) then
   call nc_check(nf90_open(trim(fname), nf90_nowrite, ncid), &
               'get_sparse_geog','open '//trim(fname))
endif

! Make sure the variables are the right size ...
! by comparing agains the size of the variable ...

if ( ngridcell < 0 ) then
   write(string1,*)'Unable to read the number of gridcells.'
   call error_handler(E_ERR,'get_sparse_geog',string1,source,revision,revdate)
endif

if ( nlandunit < 0 ) then
   write(string1,*)'Unable to read the number of land units.'
   call error_handler(E_ERR,'get_sparse_geog',string1,source,revision,revdate)
endif

if ( ncolumn < 0 ) then
   write(string1,*)'Unable to read the number of columns.'
   call error_handler(E_ERR,'get_sparse_geog',string1,source,revision,revdate)
endif

if ( npft < 0 ) then
   write(string1,*)'Unable to read the number of pfts.'
   call error_handler(E_ERR,'get_sparse_geog',string1,source,revision,revdate)
endif

! Read the netcdf file data

call nc_check(nf90_inq_varid(ncid, 'grid1d_ixy', VarID),     'get_sparse_geog', &
                         'inq_varid grid1d_ixy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   grid1d_ixy),     'get_sparse_geog', &
                                   'get_var grid1d_ixy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'grid1d_jxy', VarID),     'get_sparse_geog', &
                         'inq_varid grid1d_jxy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   grid1d_jxy),     'get_sparse_geog', &
                                   'get_var grid1d_jxy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'land1d_ixy', VarID),     'get_sparse_geog', &
                         'inq_varid land1d_ixy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   land1d_ixy),     'get_sparse_geog', &
                                   'get_var land1d_ixy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'land1d_jxy', VarID),     'get_sparse_geog', &
                         'inq_varid land1d_jxy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   land1d_jxy),     'get_sparse_geog', &
                                   'get_var land1d_jxy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'land1d_wtxy', VarID),    'get_sparse_geog', &
                         'inq_varid land1d_wtxy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   land1d_wtxy),    'get_sparse_geog', &
                                   'get_var land1d_wtxy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'cols1d_ixy', VarID),     'get_sparse_geog', &
                         'inq_varid cols1d_ixy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   cols1d_ixy),     'get_sparse_geog', &
                                   'get_var cols1d_ixy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'cols1d_jxy', VarID),     'get_sparse_geog', &
                         'inq_varid cols1d_jxy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   cols1d_jxy),     'get_sparse_geog', &
                                   'get_var cols1d_jxy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'cols1d_wtxy', VarID),    'get_sparse_geog', &
                         'inq_varid cols1d_wtxy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   cols1d_wtxy),    'get_sparse_geog', &
                                   'get_var cols1d_wtxy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'cols1d_ityplun', VarID), 'get_sparse_geog', &
                         'inq_varid cols1d_ityplun '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   cols1d_ityplun), 'get_sparse_geog', &
                                   'get_var cols1d_ityplun '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'pfts1d_ixy', VarID),     'get_sparse_geog', &
                         'inq_varid pfts1d_ixy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   pfts1d_ixy),     'get_sparse_geog', &
                                   'get_var pfts1d_ixy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'pfts1d_jxy', VarID),     'get_sparse_geog', &
                         'inq_varid pfts1d_jxy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   pfts1d_jxy),     'get_sparse_geog', &
                                   'get_var pfts1d_jxy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'pfts1d_wtxy', VarID),    'get_sparse_geog', &
                         'inq_varid pfts1d_wtxy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   pfts1d_wtxy),    'get_sparse_geog', &
                                   'get_var pfts1d_wtxy '//trim(clm_restart_filename))

! zsno is NOT optional ... so it IS a fatal error if it is not present (for now, anyway).
! as read into fortran ... zsno(:,1) is the level closest to the sun.
! as read into fortran ... zsno(:,5) is the level closest to the ground.

if (nlevsno > 0 ) then
   call nc_check(nf90_inq_varid(ncid,   'ZSNO', VarID), &
        &    'get_sparse_geog', 'inq_varid ZSNO '//trim(fname))
   call nc_check(nf90_get_var(  ncid, VarID,   zsno), &
        &    'get_sparse_geog',   'get_var ZSNO '//trim(fname))
else
   write(string1,*) 'levsno must be in restart file'
   call error_handler(E_ERR,'get_sparse_geog',string1,source,revision,revdate)
endif

if (cstat == 'close') then
   call nc_check(nf90_close(ncid),'get_sparse_geog','close '//trim(fname) )
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

   write(logfileunit,*)'cols1d_ixy     range ',minval(cols1d_ixy),    maxval(cols1d_ixy)
   write(logfileunit,*)'cols1d_jxy     range ',minval(cols1d_jxy),    maxval(cols1d_jxy)
   write(logfileunit,*)'cols1d_wtxy    range ',minval(cols1d_wtxy),   maxval(cols1d_wtxy)
   write(logfileunit,*)'cols1d_ityplun range ',minval(cols1d_ityplun),maxval(cols1d_ityplun)

   write(logfileunit,*)'pfts1d_ixy     range ',minval(pfts1d_ixy),    maxval(pfts1d_ixy)
   write(logfileunit,*)'pfts1d_jxy     range ',minval(pfts1d_jxy),    maxval(pfts1d_jxy)
   write(logfileunit,*)'pfts1d_wtxy    range ',minval(pfts1d_wtxy),   maxval(pfts1d_wtxy)
   if (nlevsno > 0) write(logfileunit,*)'zsno           range ',minval(zsno),maxval(zsno)

   write(     *     ,*)
   write(     *     ,*)'Raw lat/lon information as read ...'
   write(     *     ,*)'grid1d_ixy     range ',minval(grid1d_ixy),    maxval(grid1d_ixy)
   write(     *     ,*)'grid1d_jxy     range ',minval(grid1d_jxy),    maxval(grid1d_jxy)

   write(     *     ,*)'land1d_ixy     range ',minval(land1d_ixy),    maxval(land1d_ixy)
   write(     *     ,*)'land1d_jxy     range ',minval(land1d_jxy),    maxval(land1d_jxy)
   write(     *     ,*)'land1d_wtxy    range ',minval(land1d_wtxy),   maxval(land1d_wtxy)

   write(     *     ,*)'cols1d_ixy     range ',minval(cols1d_ixy),    maxval(cols1d_ixy)
   write(     *     ,*)'cols1d_jxy     range ',minval(cols1d_jxy),    maxval(cols1d_jxy)
   write(     *     ,*)'cols1d_wtxy    range ',minval(cols1d_wtxy),   maxval(cols1d_wtxy)
   write(     *     ,*)'cols1d_ityplun range ',minval(cols1d_ityplun),maxval(cols1d_ityplun)

   write(     *     ,*)'pfts1d_ixy     range ',minval(pfts1d_ixy),    maxval(pfts1d_ixy)
   write(     *     ,*)'pfts1d_jxy     range ',minval(pfts1d_jxy),    maxval(pfts1d_jxy)
   write(     *     ,*)'pfts1d_wtxy    range ',minval(pfts1d_wtxy),   maxval(pfts1d_wtxy)
   if (nlevsno > 0) write(     *     ,*)'zsno           range ',minval(zsno),maxval(zsno)

endif

return
end subroutine get_sparse_geog


!------------------------------------------------------------------


function get_state_time_ncid( ncid )
!------------------------------------------------------------------
! The restart netcdf files have the time of the state.

type(time_type) :: get_state_time_ncid
integer, intent(in) :: ncid

integer :: VarID
integer :: rst_curr_ymd, rst_curr_tod, leftover
integer :: year, month, day, hour, minute, second

if ( .not. module_initialized ) call static_init_model

call nc_check(nf90_inq_varid(ncid, 'timemgr_rst_curr_ymd', VarID), 'get_state_time_ncid', &
                      &  'inq_varid timemgr_rst_curr_ymd '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   rst_curr_ymd), 'get_state_time_ncid', &
                      &            'get_var rst_curr_ymd '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'timemgr_rst_curr_tod', VarID), 'get_state_time_ncid', &
                      &  'inq_varid timemgr_rst_curr_tod '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   rst_curr_tod), 'get_state_time_ncid', &
                      &            'get_var rst_curr_tod '//trim(clm_restart_filename))

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




function set_model_time_step()
!------------------------------------------------------------------
! This defines the window used for assimilation.
! all observations +/- half this timestep are assimilated.

type(time_type) :: set_model_time_step

set_model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end function set_model_time_step


!------------------------------------------------------------------


subroutine get_clm_restart_filename( filename )

character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(clm_restart_filename)

end subroutine get_clm_restart_filename



subroutine get_clm_history_filename( filename )

character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(clm_history_filename)

end subroutine get_clm_history_filename


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
      call error_handler(E_ERR, 'parse_variable_table', string1, &
         source, revision, revdate, text2=string2, text3=string3)
   endif

   ! Make sure DART kind is valid

   if( get_raw_obs_kind_index(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'parse_variable_table',string1,source,revision,revdate)
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
   call error_handler(E_MSG,'parse_variable_table',string1,text2=string2)
endif

! Check to see if zsno is part of the requested variables

end subroutine parse_variable_table


!------------------------------------------------------------------


function FindTimeDimension(ncid) result(timedimid)

! Find the Time Dimension ID in a netCDF file.
! If there is none - (spelled the obvious way) - the routine
! returns a negative number. You don't HAVE to have a TIME dimension.

integer                      :: timedimid
integer,          intent(in) :: ncid

integer :: nc_rc

TimeDimID = -1 ! same as the netCDF library routines.
nc_rc = nf90_inq_dimid(ncid,'TIME',dimid=TimeDimID)
if ( nc_rc /= NF90_NOERR ) then ! did not find it - try another spelling
   nc_rc = nf90_inq_dimid(ncid,'Time',dimid=TimeDimID)
   if ( nc_rc /= NF90_NOERR ) then ! did not find it - try another spelling
      nc_rc = nf90_inq_dimid(ncid,'time',dimid=TimeDimID)
   endif
endif

end function FindTimeDimension


!------------------------------------------------------------------


!>  define_var_dims() takes the N-dimensional variable and appends the DART
!>  dimensions of 'copy' and 'time'. If the variable initially had a 'time'
!>  dimension, it is ignored because (by construction) it is a singleton
!>  dimension.

subroutine define_var_dims(ivar, ncid, memberdimid, unlimiteddimid, ndims, dimids)

integer,               intent(in)  :: ivar, ncid, memberdimid, unlimiteddimid
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: dimids

character(len=NF90_MAX_NAME),dimension(NF90_MAX_VAR_DIMS) :: dimnames

integer :: i, mydimid

ndims = 0

DIMLOOP : do i = 1,progvar(ivar)%numdims

   if (progvar(ivar)%dimnames(i) == 'time') cycle DIMLOOP

   call nc_check(nf90_inq_dimid(ncid=ncid, name=progvar(ivar)%dimnames(i), dimid=mydimid), &
                           'define_var_dims','inq_dimid '//trim(progvar(ivar)%dimnames(i)))

   ndims         = ndims + 1
   dimids(ndims) = mydimid
   dimnames(ndims) = progvar(ivar)%dimnames(i)

enddo DIMLOOP

! The last two dimensions are always 'copy' and 'time'
ndims           = ndims + 1
dimids(ndims)   = memberdimid
dimnames(ndims) = 'copy'
ndims           = ndims + 1
dimids(ndims)   = unlimitedDimid
dimnames(ndims) = 'time'

if ((debug > 8) .and. do_output()) then

   write(logfileunit,*)
   write(logfileunit,*)'define_var_dims knowledge'

   write(logfileunit,*)trim(progvar(ivar)%varname),' has original     dimnames: ', &
                   (/( trim(progvar(ivar)%dimnames(i))//' ',i=1,progvar(ivar)%numdims) /)
   write(logfileunit,*)trim(progvar(ivar)%varname),' repackaging into dimnames: ', &
                       (/ (trim(dimnames(i))//' ',i=1,ndims) /)

   write(logfileunit,*)'thus dimids ',dimids(1:ndims)
   write(     *     ,*)
   write(     *     ,*)'define_var_dims knowledge'
   write(     *     ,*)trim(progvar(ivar)%varname),' has original     dimnames: ', &
                   (/( trim(progvar(ivar)%dimnames(i))//' ',i=1,progvar(ivar)%numdims) /)
   write(     *     ,*)trim(progvar(ivar)%varname),' repackaging into dimnames: ', &
                       (/ (trim(dimnames(i))//' ',i=1,ndims) /)
   write(     *     ,*)'thus dimids ',dimids(1:ndims)

endif

return
end subroutine define_var_dims


!------------------------------------------------------------------


  subroutine fill_levels(dimname,icol,enlevels,levtot)
! subroutine fill_levels(dimname,icol,enlevels,levtot)
!
! dimname         ... is it dimensioned 'levgrnd' or 'levsno' or 'levtot' ...
! icol            ... which CLM 'column' are we in
! enlevels        ... the expected number of levels ... varshape
! levtot          ... the output array of vertical coordinates
!
! The total number of levels is defined to be the soil levels (fixed)
! plus the number of snow levels, which can vary by column.
! The history file contains the depths of the soil levels (ncdf var 'levgrnd');
! these are in levtot(1:nlevgrnd).
! The restart file contains the depths of the snow levels (ncdf var 'ZSNO').
! The tricky bit is that they are in reverse order ... and differ from one model to another.
! If you simply grab the netcdf variable zsno,
! the level closest to the soil is actually the highest index.
!
! From Matlab (which indexes like Fortran)
!> size(zsno) ans = 13693           5
!> zsno(1,:)  ans = -1.4202   -1.3852   -1.3052   -1.1352   -0.5101
!                      |          |         |         |        |...... closest to soil surface
!                      |          |         |         |............... one level 'up'
!                      |          |         |......................... one level 'up'
!                      |          |................................... one level 'up'
!                      |.............................................. closest to sun
!
! If there is no snow ... the corresponding zsno is ZERO ...
!> zsno(508,:) ans = 0   -0.5736   -0.5389   -0.4591   -0.2021
!
! The following Matlab code may be used to explore a variable's storage order.
! (a better version is in the clm/matlab/CheckStorageOrder.m function)
!
! h2o = nc_varget(fname,'H2OSOI_LIQ');
! h2o(h2o > 1.0E30) = NaN;
! lat = nc_varget(fname,'cols1d_lat');
! lon = nc_varget(fname,'cols1d_lon');
! figure(1); plot3(lon,lat,h2o(:,1),'x'); hold on; worldmap; view(0,90)
! figure(2); plot3(lon,lat,h2o(:,2),'x'); hold on; worldmap; view(0,90)
! figure(3); plot3(lon,lat,h2o(:,3),'x'); hold on; worldmap; view(0,90)
! figure(4); plot3(lon,lat,h2o(:,4),'x'); hold on; worldmap; view(0,90)
! figure(5); plot3(lon,lat,h2o(:,5),'x'); hold on; worldmap; view(0,90)
! figure(6); plot3(lon,lat,h2o(:,6),'x'); hold on; worldmap; view(0,90)

character(len=*),          intent(in) :: dimname
integer,                   intent(in) :: icol
integer,                   intent(in) :: enlevels
real(r8), dimension(:), intent(inout) :: levtot

if     (dimname == 'levsno') then

   if (nlevsno /= enlevels) then
      write(string1,*) 'dimension ', trim(dimname),' has declared length ',enlevels
      write(string2,*) 'not the known number of snow levels ',nlevsno
      call error_handler(E_ERR,'fill_levels', string1, &
                             source, revision, revdate, text2=string2)
   endif
   levtot(1:nlevsno) = zsno(1:nlevsno,icol)

elseif (dimname == 'levgrnd') then

   if (nlevgrnd /= enlevels) then
      write(string1,*) 'dimension ', trim(dimname),' has declared length ',enlevels
      write(string2,*) 'not the known number of soil levels ',nlevgrnd
      call error_handler(E_ERR,'fill_levels', string1, &
                             source, revision, revdate, text2=string2)
   endif
   levtot(1:nlevgrnd) = LEVGRND

elseif (dimname == 'levtot') then

   ! This block assumes anything dimensioned 'levtot' has the first nlevsno levels
   ! followed by nlevgrnd levels. Dunno what to do with lake stuff ...

   if (nlevtot /= enlevels) then
      write(string1,*) 'dimension ', trim(dimname),' has declared length ',enlevels
      write(string2,*) 'not the known number of total levels ',nlevtot
      call error_handler(E_ERR,'fill_levels', string1, &
                             source, revision, revdate, text2=string2)
   endif

   if (nlevtot /= nlevgrnd + nlevsno) then
      write(string1,*) 'nlevtot ', nlevtot,' is not equal to nlevgrnd + nlevsno'
      write(string2,*) 'nlevgrnd is ',nlevgrnd,' nlevsno is ',nlevsno,' total of ',nlevgrnd+nlevsno
      call error_handler(E_ERR,'fill_levels', string1, &
                             source, revision, revdate, text2=string2)
   endif

   levtot(1:nlevsno) = zsno(1:nlevsno,icol)
   levtot(nlevsno+1:nlevsno+nlevgrnd) = LEVGRND

else
   write(string1,*) 'Unable to determine vertical coordinate for column ',icol
   write(string2,*) 'unknown dimension name: ',trim(dimname)
   call error_handler(E_ERR,'fill_levels', string1, &
                             source, revision, revdate, text2=string2)
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

integer :: ivar, indexi, i, j
integer, allocatable, dimension(:,:) :: countmat

if ( .not. module_initialized ) call static_init_model

allocate(countmat(nlon,nlat))
countmat = 0

VARTYPES : do ivar = 1,nfields

    ! Skip to the right variable
    if ( trim(progvar(ivar)%varname) /= varstring) cycle VARTYPES

    ! Create a count of all the multiples in a gridcell
    do indexi = progvar(ivar)%index1, progvar(ivar)%indexN
       i = lonixy(indexi)
       j = latjxy(indexi)
       countmat(i,j) = countmat(i,j) + 1
    enddo

enddo VARTYPES

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



subroutine get_var_1d(ncid, varname, var1d)
! This function will return a R8 array with the netCDF attributes applied.
! scale_factor, offset will be applied,
! missing_value, _FillValue will be replaced by the DART missing value ...

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

integer,                intent(in)  :: ncid
character(len=*),       intent(in)  :: varname
real(r8), dimension(:), intent(out) :: var1d

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens, ncstart, nccount
integer  :: VarID, numdims, xtype, io1, io2
integer  :: TimeDimID, time_dimlen, timeindex
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8

integer,  allocatable, dimension(:) :: intarray
real(r4), allocatable, dimension(:) :: r4array

if ( .not. module_initialized ) call static_init_model

! a little whitespace makes this a lot more readable
if (do_output() .and. (debug > 1)) then
   write(*,*)
   write(logfileunit,*)
endif

io1 = nf90_inq_dimid(ncid, 'time', TimeDimID)
if (io1 /= NF90_NOERR) TimeDimID = MISSING_I

call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), 'get_var_1d', 'inq_varid '//varname)
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_1d', 'inquire_variable '//varname)
call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlens(1)), &
              'get_var_1d', 'inquire_dimension '//varname)

if ((numdims /= 1) .or. (size(var1d) /= dimlens(1)) ) then
   write(string1,*) trim(varname)//' is not the expected shape/length of ', size(var1d)
   call error_handler(E_ERR,'get_var_1d',string1,source,revision,revdate)
endif

ncstart = 1
nccount = dimlens(1)

if (dimIDs(1) == TimeDimID) then
   call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), &
         'get_var_1d', 'inquire_dimension time '//trim(varname))
   timeindex  = FindDesiredTimeIndx(ncid, time_dimlen, varname)
   ncstart(1) = timeindex
   nccount(1) = 1
endif

if (do_output() .and. (debug > 1)) then
   write(*,*)'get_var_1d: variable ['//trim(varname)//']'
   write(*,*)'get_var_1d: start ',ncstart(1:numdims)
   write(*,*)'get_var_1d: count ',nccount(1:numdims)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1)))
   call nc_check(nf90_get_var(ncid, VarID, values=intarray, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_1d', 'get_var '//varname)
   var1d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d + add_offset
   endif

   deallocate(intarray)

elseif (xtype == NF90_FLOAT) then

   allocate(r4array(dimlens(1)))
   call nc_check(nf90_get_var(ncid, VarID, values=r4array, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_1d', 'get_var '//varname)
   var1d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d + add_offset
   endif

   deallocate(r4array)

elseif (xtype == NF90_DOUBLE) then

   call nc_check(nf90_get_var(ncid, VarID, values=var1d, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_1d', 'get_var '//varname)

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var1d == spvalR8) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var1d == spvalR8) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d + add_offset
   endif

else
   write(string1,*) trim(varname)//' has unsupported (by DART) xtype of', xtype
   call error_handler(E_ERR,'get_var_1d',string1,source,revision,revdate)
endif

end subroutine get_var_1d



subroutine get_var_2d(ncid, varname, var2d)
! This function will return a R8 array with the netCDF attributes applied.
! scale_factor, offset will be applied,
! missing_value, _FillValue will be replaced by the DART missing value ...

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

integer,                  intent(in)  :: ncid
character(len=*),         intent(in)  :: varname
real(r8), dimension(:,:), intent(out) :: var2d

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens, ncstart, nccount
integer  :: VarID, numdims, xtype, io1, io2, i
integer  :: TimeDimID, time_dimlen, timeindex
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8

integer,  allocatable, dimension(:,:) :: intarray
real(r4), allocatable, dimension(:,:) :: r4array

if ( .not. module_initialized ) call static_init_model

! a little whitespace makes this a lot more readable
if (do_output() .and. (debug > 1)) then
   write(*,*)
   write(logfileunit,*)
endif

io1 = nf90_inq_dimid(ncid, 'time', TimeDimID)
if (io1 /= NF90_NOERR) TimeDimID = MISSING_I

call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), 'get_var_2d', 'inq_varid')
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_2d', 'inquire_variable')

if ( (numdims /= 2)  ) then
   write(string1,*) trim(varname)//' is not a 2D variable as expected.'
   call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate)
endif

ncstart(:) = 1
nccount(:) = 1

DimCheck : do i = 1,numdims

   write(string1,'(a,i2,1x,A)') 'inquire dimension ',i,trim(varname)

   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
              'get_var_2d', string1)

   if ( dimIDs(i) == TimeDimID ) then

      call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen), &
            'get_var_2d', 'inquire_dimension time '//trim(varname))

      timeindex  = FindDesiredTimeIndx(ncid, time_dimlen, varname)
      ncstart(i) = timeindex
      dimlens(i) = 1

   elseif ( size(var2d,i) /= dimlens(i) ) then
      write(string1,*) trim(varname)//' has shape ', dimlens(1:numdims)
      write(string2,*) 'which is not the expected shape of ', size(var2d,1),size(var2d,2)
      call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate,text2=string2)
   endif

   nccount(i) = dimlens(i)

enddo DimCheck

if (do_output() .and. (debug > 1)) then
   write(*,*)'get_var_2d: variable ['//trim(varname)//']'
   write(*,*)'get_var_2d: start ',ncstart(1:numdims)
   write(*,*)'get_var_2d: count ',nccount(1:numdims)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1),dimlens(2)))
   call nc_check(nf90_get_var(ncid, VarID, values=intarray, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_2d', 'get_var '//varname)
   var2d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d + add_offset
   endif

   deallocate(intarray)

elseif (xtype == NF90_FLOAT) then

   allocate(r4array(dimlens(1),dimlens(2)))
   call nc_check(nf90_get_var(ncid, VarID, values=r4array, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_2d', 'get_var '//varname)
   var2d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d + add_offset
   endif

   deallocate(r4array)

elseif (xtype == NF90_DOUBLE) then

   call nc_check(nf90_get_var(ncid, VarID, values=var2d, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_2d', 'get_var '//varname)

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var2d == spvalR8) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var2d == spvalR8) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d + add_offset
   endif

else
   write(string1,*) trim(varname)//' has unsupported (by DART) xtype of', xtype
   call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate)
endif

end subroutine get_var_2d



subroutine get_var_3d(ncid, varname, var3d)
! This function will return a R8 array with the netCDF attributes applied.
! scale_factor, offset will be applied,
! missing_value, _FillValue will be replaced by the DART missing value ...

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

integer,                    intent(in)  :: ncid
character(len=*),           intent(in)  :: varname
real(r8), dimension(:,:,:), intent(out) :: var3d

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens, ncstart, nccount
integer  :: i, TimeDimID, time_dimlen, timeindex
integer  :: VarID, numdims, xtype, io1, io2
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8

integer,  allocatable, dimension(:,:,:) :: intarray
real(r4), allocatable, dimension(:,:,:) :: r4array

if ( .not. module_initialized ) call static_init_model

! a little whitespace makes this a lot more readable
if (do_output() .and. (debug > 1)) then
   write(*,*)
   write(logfileunit,*)
endif

! 3D fields must have a time dimension.
! Need to know the Time Dimension ID and length

call nc_check(nf90_inq_dimid(ncid, 'time', TimeDimID), &
         'get_var_3d', 'inq_dimid time '//varname)
call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), &
         'get_var_3d', 'inquire_dimension time '//varname)

call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), 'get_var_3d', 'inq_varid')
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_3d', 'inquire_variable '//varname)

if ( (numdims /= 3)  ) then
   write(string1,*) trim(varname)//' is not a 3D variable as expected.'
   call error_handler(E_ERR,'get_var_3d',string1,source,revision,revdate)
endif

! only expecting [nlon,nlat,time]

ncstart(:) = 1
nccount(:) = 1
DimCheck : do i = 1,numdims

   write(string1,'(a,i2,1x,A)') 'inquire dimension ',i,trim(varname)

   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
                 'get_var_3d', trim(string1))

   if ( dimIDs(i) == TimeDimID ) then

       timeindex  = FindDesiredTimeIndx(ncid, time_dimlen, varname)
       ncstart(i) = timeindex
       dimlens(i) = 1

   elseif ( size(var3d,i) /= dimlens(i) ) then
      write(string1,*) trim(varname)//' has shape ', dimlens(1:numdims)
      write(string2,*) 'which is not the expected shape of ', size(var3d,i)
      call error_handler(E_ERR,'get_var_3d',string1,source,revision,revdate,text2=string2)
   endif

   nccount(i) = dimlens(i)

enddo DimCheck

if (do_output() .and. (debug > 1)) then
   write(*,*)'get_var_3d: variable ['//trim(varname)//']'
   write(*,*)'get_var_3d: start ',ncstart(1:numdims)
   write(*,*)'get_var_3d: count ',nccount(1:numdims)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1),dimlens(2),dimlens(3)))
   call nc_check(nf90_get_var(ncid, VarID, values=intarray, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_3d', 'get_var '//varname)
   var3d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var3d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var3d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d + add_offset
   endif

   deallocate(intarray)

elseif (xtype == NF90_FLOAT) then

   allocate(r4array(dimlens(1),dimlens(2),dimlens(3)))
   call nc_check(nf90_get_var(ncid, VarID, values=r4array, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_3d', 'get_var '//varname)
   var3d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var3d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var3d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d + add_offset
   endif

   deallocate(r4array)

elseif (xtype == NF90_DOUBLE) then

   call nc_check(nf90_get_var(ncid, VarID, values=var3d, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_3d', 'get_var '//varname)

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var3d == spvalR8) var3d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var3d == spvalR8) var3d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d + add_offset
   endif

else
   write(string1,*) trim(varname)//' has unsupported (by DART) xtype of', xtype
   call error_handler(E_ERR,'get_var_3d',string1,source,revision,revdate)
endif

end subroutine get_var_3d



function get_model_time()
type(time_type) :: get_model_time

if ( .not. module_initialized ) call static_init_model

get_model_time = model_time

end function get_model_time



function findKindIndex(kind_index, caller)
integer,          intent(in) :: kind_index
character(len=*), intent(in) :: caller
integer                      :: findKindIndex

integer :: i
character(len=paramname_length) :: kind_string
findKindIndex = -1

! Skip to the right variable
VARTYPES : do i = 1,nfields
    if (progvar(i)%dart_kind == kind_index) then
       findKindIndex = i
       exit VARTYPES
    endif
enddo VARTYPES

if (findKindIndex < 1) then
   kind_string = get_raw_obs_kind_name( kind_index )
   write(string1,*) trim(caller)//' cannot find "'//trim(kind_string)//'" in list of DART state variables.'
   write(string2,*) trim(caller)//' looking for DART KIND (index) ', kind_index
   call error_handler(E_ERR,'findKindIndex',string1,source,revision,revdate, text2=string2)
endif

end function findKindIndex



subroutine SetLocatorArrays()
! This function will create the relational table that will indicate how many
! and which columns pertain to the gridcells. A companion function will
! return the column indices that are needed to recreate the gridcell value.
!
! This fills the gridCellInfo(:,:) structure.
! given a gridcell, the gridCellInfo(:,:) structure will indicate how many and
! which columns are part of the gridcell.

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

do ilon = 1,nlon
do ilat = 1,nlat
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
      call error_handler(E_ERR, 'SetLocatorArrays', string1, &
                   source, revision, revdate, text2=string2)
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
      call error_handler(E_ERR, 'SetLocatorArrays', string1, &
                   source, revision, revdate, text2=string2)
   endif
enddo

! Check block

if ((debug > 99) .and. do_output()) then

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


!> SetVariableAttributes() converts the information in the variable_table 
!> to the progvar structure for each variable.
!> If the numerical limit does not apply, it is set to MISSING_R8, even if
!> it is the maximum that does not apply.

subroutine SetVariableAttributes(ivar)

integer, intent(in) :: ivar

integer  :: ios
real(r8) :: minvalue, maxvalue

progvar(ivar)%varname     = trim(variable_table(ivar,VT_VARNAMEINDX))
progvar(ivar)%kind_string = trim(variable_table(ivar,VT_KINDINDX))
progvar(ivar)%dart_kind   = get_raw_obs_kind_index( progvar(ivar)%kind_string )
progvar(ivar)%maxlevels   = 0
progvar(ivar)%dimlens     = 0
progvar(ivar)%dimnames    = ' '
progvar(ivar)%spvalINT    = -9999        ! from CESM/CLM clm_varcon.F90
progvar(ivar)%spvalR4     = 1.e36_r4     ! from CESM/CLM clm_varcon.F90
progvar(ivar)%spvalR8     = 1.e36_r8     ! from CESM/CLM clm_varcon.F90
progvar(ivar)%missingINT  = MISSING_I
progvar(ivar)%missingR4   = MISSING_R4
progvar(ivar)%missingR8   = MISSING_R8
progvar(ivar)%rangeRestricted   = BOUNDED_NONE
progvar(ivar)%minvalue          = MISSING_R8
progvar(ivar)%maxvalue          = MISSING_R8
progvar(ivar)%has_fill_value    = .true.
progvar(ivar)%has_missing_value = .true.
progvar(ivar)%update            = .false.

if (variable_table(ivar,VT_ORIGININDX) == 'VECTOR') then
   progvar(ivar)%origin = trim(clm_vector_history_filename)
elseif (variable_table(ivar,VT_ORIGININDX) == 'HISTORY') then
   progvar(ivar)%origin = trim(clm_history_filename)
else
   variable_table(ivar,VT_ORIGININDX) = 'RESTART'
   progvar(ivar)%origin = trim(clm_restart_filename)
endif

if ((variable_table(ivar,VT_STATEINDX)  == 'UPDATE') .and. &
    (variable_table(ivar,VT_ORIGININDX) == 'RESTART')) progvar(ivar)%update = .true.

! set the default values

minvalue = MISSING_R8
maxvalue = MISSING_R8
progvar(ivar)%minvalue = MISSING_R8
progvar(ivar)%maxvalue = MISSING_R8

! If the character string can be interpreted as an r8, great.
! If not, there is no value to be used.

read(variable_table(ivar,VT_MINVALINDX),*,iostat=ios) minvalue
if (ios == 0) progvar(ivar)%minvalue = minvalue

read(variable_table(ivar,VT_MAXVALINDX),*,iostat=ios) maxvalue
if (ios == 0) progvar(ivar)%maxvalue = maxvalue

! rangeRestricted == BOUNDED_NONE  == 0 ... unlimited range
! rangeRestricted == BOUNDED_BELOW == 1 ... minimum, but no maximum
! rangeRestricted == BOUNDED_ABOVE == 2 ... maximum, but no minimum
! rangeRestricted == BOUNDED_BOTH  == 3 ... minimum and maximum

if (   (progvar(ivar)%minvalue /= MISSING_R8) .and. &
       (progvar(ivar)%maxvalue /= MISSING_R8) ) then
   progvar(ivar)%rangeRestricted = BOUNDED_BOTH

elseif (progvar(ivar)%maxvalue /= MISSING_R8) then
   progvar(ivar)%rangeRestricted = BOUNDED_ABOVE

elseif (progvar(ivar)%minvalue /= MISSING_R8) then
   progvar(ivar)%rangeRestricted = BOUNDED_BELOW

else
   progvar(ivar)%rangeRestricted = BOUNDED_NONE

endif

! Check to make sure min is less than max if both are specified.

if ( progvar(ivar)%rangeRestricted == BOUNDED_BOTH ) then
   if (maxvalue < minvalue) then
      write(string1,*)'&model_nml state_variable input error for ',trim(progvar(ivar)%varname)
      write(string2,*)'minimum value (',minvalue,') must be less than '
      write(string3,*)'maximum value (',maxvalue,')'
      call error_handler(E_ERR,'SetVariableAttributes',string1, &
         source,revision,revdate,text2=trim(string2),text3=trim(string3))
   endif
endif

end subroutine SetVariableAttributes




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
   call error_handler(E_ERR, 'FindDesiredTimeIndx:', string1, &
          source, revision, revdate, text2=string2)
endif

read(attvalue,'(11x,i4,5(1x,i2))',iostat=ios)iyear,imonth,iday,ihour,imin,isec
if (ios /= 0) then
   write(string1,*)'Unable to read time units. Error status was ',ios
   write(string2,*)'expected "days since YYYY-MM-DD HH:MM:SS"'
   write(string3,*)'was      "'//trim(attvalue)//'"'
   call error_handler(E_ERR, 'FindDesiredTimeIndx:', string1, &
          source, revision, revdate, text2=string2, text3=string3)
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
   call error_handler(E_ERR, 'FindDesiredTimeIndx:', string1, &
          source, revision, revdate )
endif

if ((debug > 0) .and. do_output()) then
   write(string1,*)trim(varname)//' matching time index is ',FindDesiredTimeIndx
   call error_handler(E_MSG, 'FindDesiredTimeIndx:', string1)
endif

end function FindDesiredTimeIndx



subroutine cluster_variables(filename, nvars, var_names, var_ranges)
character(len=*), intent(in)  :: filename
integer,          intent(out) :: nvars
character(len=*), intent(out) :: var_names(:)
real(r8),         intent(out) :: var_ranges(:,:)

integer :: ivar

nvars      = 0
var_names  = 'no_variable_specified'
var_ranges = MISSING_R8

do ivar = 1,nfields

   if (trim(progvar(ivar)%origin) == trim(filename)) then
      nvars = nvars + 1
      var_names(nvars)    = progvar(ivar)%varname
      var_ranges(nvars,1) = progvar(ivar)%minvalue
      var_ranges(nvars,2) = progvar(ivar)%maxvalue
   endif

enddo

!>@todo FIXME ... this summary message is not finished
if (do_output() .and. debug > 99) then
   do ivar = 1,nvars
      write(string1,*)trim(filename),' has ',trim(var_names(ivar)),var_ranges(ivar,1),var_ranges(ivar,2)
      call error_handler(E_MSG,'cluster_variables',string1)
   enddo
endif

end subroutine cluster_variables

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
