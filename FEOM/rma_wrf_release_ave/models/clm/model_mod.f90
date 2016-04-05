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
                             obstypelength
use time_manager_mod, only : time_type, set_time, set_date, get_time,          &
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, get_dist, query_location,          &
                             get_close_maxdist_init, get_close_type,           &
                             set_location, get_location, horiz_dist_only,      &
                             vert_is_undef,    VERTISUNDEF,                    &
                             vert_is_surface,  VERTISSURFACE,                  &
                             vert_is_level,    VERTISLEVEL,                    &
                             vert_is_pressure, VERTISPRESSURE,                 &
                             vert_is_height,   VERTISHEIGHT,                   &
                             get_close_obs_init, get_close_obs

use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper,                    &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, file_exist, find_textfile_dims,        &
                             file_to_text

use     obs_kind_mod, only : KIND_SOIL_TEMPERATURE,   &
                             KIND_SOIL_MOISTURE,      &
                             KIND_LIQUID_WATER,       &
                             KIND_ICE,                &
                             KIND_SNOWCOVER_FRAC,     &
                             KIND_SNOW_THICKNESS,     &
                             KIND_LEAF_CARBON,        &
                             KIND_WATER_TABLE_DEPTH,  &
                             KIND_GEOPOTENTIAL_HEIGHT,&
                             paramname_length,        &
                             get_raw_obs_kind_index

use mpi_utilities_mod, only: my_task_id

use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

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

public :: get_gridsize,                 &
          restart_file_to_sv,           &
          sv_to_restart_file,           &
          get_clm_restart_filename,     &
          get_state_time,               &
          get_grid_vertval,             &
          compute_gridcell_value,       &
          find_gridcell_Npft,           &
          DART_get_var,                 &
          get_model_time

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

integer :: nfields
integer, parameter :: max_state_variables = 40
integer, parameter :: num_state_table_columns = 2
character(len=obstypelength) :: variable_table(max_state_variables, num_state_table_columns)

! things which can/should be in the model_nml

integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 60
real(r8)           :: model_perturbation_amplitude = 0.2
logical            :: output_state_vector = .true.
integer            :: debug = 0   ! turn up for more and more debug messages
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: clm_restart_filename = 'clm_restart.nc'
character(len=256) :: clm_history_filename = 'clm_history.nc'
character(len=obstypelength) :: clm_state_variables(max_state_variables*num_state_table_columns) = ' '

namelist /model_nml/            &
   clm_restart_filename,        &
   clm_history_filename,        &
   output_state_vector,         &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   calendar,                    &
   debug,                       &
   clm_state_variables

! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=obstypelength), dimension(NF90_MAX_VAR_DIMS) :: dimnames
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer :: numdims
   integer :: maxlevels
   integer :: xtype
   integer :: varsize     ! prod(dimlens(1:numdims))
   integer :: index1      ! location in dart state vector of first occurrence
   integer :: indexN      ! location in dart state vector of last  occurrence
   integer :: dart_kind
   integer :: spvalINT, missingINT
   real(r4) :: spvalR4, missingR4
   real(r8) :: spvalR8, missingR8
   logical  :: has_fill_value
   logical  :: has_missing_value
   character(len=paramname_length) :: kind_string
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

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

integer :: Ngridcell = -1 ! Number of gridcells containing land
integer :: Nlandunit = -1 ! Number of land units
integer :: Ncolumn   = -1 ! Number of columns
integer :: Npft      = -1 ! Number of plant functional types
integer :: Nlevlak   = -1 ! Number of
integer :: Nlevsno   = -1 ! Number of snow levels
integer :: Nlevsno1  = -1 ! Number of snow level ... interfaces?
integer :: Nlevtot   = -1 ! Number of total levels
integer :: Nnumrad   = -1 ! Number of
integer :: Nrtmlon   = -1 ! Number of river transport model longitudes
integer :: Nrtmlat   = -1 ! Number of river transport model latitudes

integer,  allocatable, dimension(:)  :: grid1d_ixy, grid1d_jxy ! 2D lon/lat index of corresponding gridcell
integer,  allocatable, dimension(:)  :: land1d_ixy, land1d_jxy ! 2D lon/lat index of corresponding gridcell
integer,  allocatable, dimension(:)  :: cols1d_ixy, cols1d_jxy ! 2D lon/lat index of corresponding gridcell
integer,  allocatable, dimension(:)  :: pfts1d_ixy, pfts1d_jxy ! 2D lon/lat index of corresponding gridcell
real(r8), allocatable, dimension(:)  :: land1d_wtxy    ! landunit weight relative to corresponding gridcell
real(r8), allocatable, dimension(:)  :: cols1d_wtxy    ! column   weight relative to corresponding gridcell
real(r8), allocatable, dimension(:)  :: pfts1d_wtxy    ! pft      weight relative to corresponding gridcell
real(r8), allocatable, dimension(:)  :: levtot
real(r8), allocatable, dimension(:,:):: zsno           ! (column,levsno) ... snow layer midpoint

!------------------------------------------------------------------------------
! These are the metadata arrays that are the same size as the state vector.

real(r8), allocatable, dimension(:) :: ens_mean     ! may be needed for forward ops
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
END INTERFACE

INTERFACE get_state_time
      MODULE PROCEDURE get_state_time_ncid
      MODULE PROCEDURE get_state_time_fname
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

end subroutine adv_1step



subroutine get_state_meta_data(indx, location, var_type)
!------------------------------------------------------------------
! Given an integer index into the state vector structure, returns the
! associated array indices for lat, lon, and height, as well as the type.

! Passed variables

integer, intent(in)            :: indx
type(location_type)            :: location
integer, optional, intent(out) :: var_type

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
character(len=NF90_MAX_NAME)          :: varname
character(len=obstypelength)          :: dimname
character(len=paramname_length)       :: kind_string
integer :: ncid, VarID, dimlen, varsize
integer :: iunit, io, ivar
integer :: i, j, xi, xj, index1, indexN, indx
integer :: ss, dd

integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: spvalR8

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
if (do_output()) call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
if (do_output()) write(logfileunit, nml=model_nml)
if (do_output()) write(     *     , nml=model_nml)

!---------------------------------------------------------------
! Set the time step ... causes clm namelists to be read.
! Ensures model_timestep is multiple of 'dynamics_timestep'

call set_calendar_type( calendar )   ! comes from model_mod_nml

model_time     = get_state_time(clm_restart_filename)
model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)

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

allocate(grid1d_ixy(Ngridcell), grid1d_jxy(Ngridcell))
allocate(land1d_ixy(Nlandunit), land1d_jxy(Nlandunit), land1d_wtxy(Nlandunit))
allocate(cols1d_ixy(Ncolumn),   cols1d_jxy(Ncolumn)  , cols1d_wtxy(Ncolumn) )
allocate(pfts1d_ixy(Npft),      pfts1d_jxy(Npft)     , pfts1d_wtxy(Npft))
allocate(levtot(Nlevtot))
if (Nlevsno > 0) allocate(zsno(Nlevsno,Ncolumn))

call get_sparse_geog(ncid, clm_restart_filename, 'open')

!---------------------------------------------------------------
! Compile the list of clm variables to use in the creation
! of the DART state vector. Required to determine model_size.
!
! Verify all variables are in the clm restart file
!
! Compute the offsets into the state vector for the start of each
! variable type. Requires reading shapes from the clm restart file.
! Record the extent of the data type in the state vector.

call verify_state_variables( clm_state_variables, ncid, clm_restart_filename, &
                             nfields, variable_table )

index1  = 1;
indexN  = 0;
do ivar = 1, nfields

   varname                   = trim(variable_table(ivar,1))
   kind_string               = trim(variable_table(ivar,2))
   progvar(ivar)%varname     = varname
   progvar(ivar)%kind_string = kind_string
   progvar(ivar)%dart_kind   = get_raw_obs_kind_index( progvar(ivar)%kind_string )
   progvar(ivar)%dimlens     = 0
   progvar(ivar)%dimnames    = ' '
   progvar(ivar)%spvalINT    = MISSING_I
   progvar(ivar)%spvalR4     = MISSING_R4
   progvar(ivar)%spvalR8     = MISSING_R8
   progvar(ivar)%missingINT  = MISSING_I
   progvar(ivar)%missingR4   = MISSING_R4
   progvar(ivar)%missingR8   = MISSING_R8
   progvar(ivar)%has_fill_value    = .false.
   progvar(ivar)%has_missing_value = .false.
   progvar(ivar)%maxlevels   = 0

   string2 = trim(clm_restart_filename)//' '//trim(varname)

   call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), &
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
      progvar(ivar)%long_name = varname
   endif

   if( nf90_inquire_attribute(    ncid, VarID, 'units') == NF90_NOERR )  then
      call nc_check( nf90_get_att(ncid, VarID, 'units' , progvar(ivar)%units), &
                  'static_init_model', 'get_att units '//trim(string2))
   else
      progvar(ivar)%units = 'unknown'
   endif

   ! Saving any FillValue so I can use it when I read and write ...

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
      progvar(ivar)%dimlens( i) = dimlen
      progvar(ivar)%dimnames(i) = dimname
      varsize = varsize * dimlen

   enddo DimensionLoop

   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1
   index1                    = index1 + varsize      ! sets up for next variable

   if ((debug > 8) .and. do_output()) then
      write(logfileunit,*)
      write(logfileunit,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(logfileunit,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(logfileunit,*) '  units       ',trim(progvar(ivar)%units)
      write(logfileunit,*) '  xtype       ',progvar(ivar)%xtype
      write(logfileunit,*) '  dimnames    ',progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
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

      write(     *     ,*)
      write(     *     ,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(     *     ,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(     *     ,*) '  units       ',trim(progvar(ivar)%units)
      write(     *     ,*) '  xtype       ',progvar(ivar)%xtype
      write(     *     ,*) '  dimnames    ',progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
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
   endif

enddo

call nc_check(nf90_close(ncid),'static_init_model','close '//trim(clm_restart_filename))
ncid = 0

model_size = progvar(nfields)%indexN

if ((debug > 8) .and. do_output()) then
  write(logfileunit, *)
  write(logfileunit,'("grid: nlon, nlat, nz =",2(1x,i6))') nlon, nlat
  write(logfileunit, *)'model_size = ', model_size
  write(     *     , *)
  write(     *     ,'("grid: nlon, nlat, nz =",2(1x,i6))') nlon, nlat
  write(     *     , *)'model_size = ', model_size
endif

allocate(ens_mean(model_size))

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

   if (progvar(ivar)%numdims == 1) then

      if ((debug > 8) .and. do_output()) then
         write(*,*)
         write(*,*)'variable ',trim(progvar(ivar)%varname)
         write(*,*)'dimension 1 (i) ',progvar(ivar)%dimnames(1),progvar(ivar)%dimlens(1)
      endif

      SELECT CASE ( trim(progvar(ivar)%dimnames(1)) )
         CASE ("gridcell")
            do i = 1, progvar(ivar)%dimlens(1)
               xi             = grid1d_ixy(i)
               xj             = grid1d_jxy(i) ! always unity if unstructured
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
               xj             = land1d_jxy(i) ! always unity if unstructured
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
               xj             = cols1d_jxy(i) ! always unity if unstructured
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
               xj             = pfts1d_jxy(i) ! always unity if unstructured
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
             ' while trying to create metadata arrays.'
            call error_handler(E_ERR,'static_init_model', string1, source, revision, revdate)

      END SELECT

   elseif (progvar(ivar)%numdims == 2) then

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
               xj = grid1d_jxy(j) ! always unity if unstructured
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
               xj = land1d_jxy(j) ! always unity if unstructured
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
               xj = cols1d_jxy(j) ! always unity if unstructured
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
               xj = pfts1d_jxy(j) ! always unity if unstructured
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

         CASE DEFAULT
            write(string1,*)'(2d) unknown Dimension name '//trim(progvar(ivar)%dimnames(2))// &
             ' while trying to create metadata arrays.'
            call error_handler(E_ERR,'static_init_model', string1, source, revision, revdate)

      END SELECT

   else

      write(string1,*)'variables of rank ',progvar(ivar)%numdims,' are unsupported.'
      write(string2,*)trim(progvar(ivar)%varname),' is dimensioned ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      call error_handler(E_ERR,'static_init_model', string1, source, revision, revdate, text2=string2)

   endif

   indx = indx - 1
   if (indx /= progvar(ivar)%indexN ) then
      write(string1,*)'variable ',trim(progvar(ivar)%varname),' is supposed to end at index ',progvar(ivar)%indexN
      write(string2,*)'it ends at index ',indx
      call error_handler(E_ERR,'static_init_model', string1, source, revision, revdate, text2=string2)
   endif

enddo

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
deallocate(cols1d_ixy, cols1d_jxy, cols1d_wtxy)
deallocate(pfts1d_ixy, pfts1d_jxy, pfts1d_wtxy)

deallocate(ens_mean)
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

   call nc_check(nf90_def_dim(ncid=ncFileID, name='lon', &
          len = nlon, dimid = nlonDimID),'nc_write_model_atts', 'lon def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='lat', &
          len = nlat, dimid = NlatDimID),'nc_write_model_atts', 'lat def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='gridcell', len = ngridcell, &
             dimid = gridcellDimID),'nc_write_model_atts', 'gridcell def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='landunit', len = nlandunit, &
             dimid = landunitDimID),'nc_write_model_atts', 'landunit def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='column',   len = ncolumn,   &
             dimid =   columnDimID),'nc_write_model_atts',   'column def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='pft',      len = npft,      &
             dimid =      pftDimID),'nc_write_model_atts',      'pft def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='levgrnd',  len = nlevgrnd,  &
             dimid =  levgrndDimID),'nc_write_model_atts',  'levgrnd def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='levlak',   len = nlevlak,   &
             dimid =   levlakDimID),'nc_write_model_atts',   'levlak def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='levsno',   len = nlevsno,   &
             dimid =   levsnoDimID),'nc_write_model_atts',   'levsno def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='levsno1',  len = nlevsno1,  &
             dimid =  levsno1DimID),'nc_write_model_atts',  'levsno1 def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='levtot',   len = nlevtot,   &
             dimid =   levtotDimID),'nc_write_model_atts',   'levtot def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='numrad',   len = nnumrad,   &
             dimid =   numradDimID),'nc_write_model_atts',   'numrad def_dim '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Coordinate Variables and the Attributes
   !----------------------------------------------------------------------------

   ! Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='lon', xtype=nf90_real, &
                 dimids=(/ NlonDimID /), varid=VarID),&
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
                 dimids=(/ NlonDimID /), varid=VarID),&
                 'nc_write_model_atts', 'area def_var '//trim(filename))
   else
      call nc_check(nf90_def_var(ncFileID,name='area', xtype=nf90_real, &
                 dimids=(/ NlonDimID,nlatDimID /), varid=VarID),&
                 'nc_write_model_atts', 'area def_var '//trim(filename))
   endif
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'grid cell areas'), &
                 'nc_write_model_atts', 'area long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'km^2'),  &
                 'nc_write_model_atts', 'area units '//trim(filename))

   ! grid cell land fractions
   if (unstructured) then
      call nc_check(nf90_def_var(ncFileID,name='landfrac', xtype=nf90_real, &
                 dimids=(/ NlonDimID /), varid=VarID),&
                 'nc_write_model_atts', 'landfrac def_var '//trim(filename))
   else
      call nc_check(nf90_def_var(ncFileID,name='landfrac', xtype=nf90_real, &
                 dimids=(/ NlonDimID,nlatDimID /), varid=VarID),&
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
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', '2d longitude index of corresponding column'), &
                 'nc_write_model_atts', 'cols1d_ixy long_name '//trim(filename))

   ! latitude grid index for each column
   call nc_check(nf90_def_var(ncFileID,name='cols1d_jxy', xtype=nf90_int, &
                 dimids=(/ columnDimID /), varid=VarID),&
                 'nc_write_model_atts', 'cols1d_jxy def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', '2d latitude index of corresponding column'), &
                 'nc_write_model_atts', 'cols1d_jxy long_name '//trim(filename))

   ! column weight relative to corresponding gridcell
   call nc_check(nf90_def_var(ncFileID,name='cols1d_wtxy', xtype=nf90_double, &
                 dimids=(/ columnDimID /), varid=VarID),&
                 'nc_write_model_atts', 'cols1d_wtxy def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'column weight relative to corresponding gridcell'), &
                 'nc_write_model_atts', 'cols1d_wtxy long_name '//trim(filename))

   ! longitude grid index for each pft
   call nc_check(nf90_def_var(ncFileID,name='pfts1d_ixy', xtype=nf90_int, &
                 dimids=(/ pftDimID /), varid=VarID),&
                 'nc_write_model_atts', 'pfts1d_ixy def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', '2d longitude index of corresponding column'), &
                 'nc_write_model_atts', 'pfts1d_ixy long_name '//trim(filename))

   ! latitude grid index for each pft
   call nc_check(nf90_def_var(ncFileID,name='pfts1d_jxy', xtype=nf90_int, &
                 dimids=(/ pftDimID /), varid=VarID),&
                 'nc_write_model_atts', 'pfts1d_jxy def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', '2d latitude index of corresponding column'), &
                 'nc_write_model_atts', 'pfts1d_jxy long_name '//trim(filename))

   ! pft weight relative to corresponding gridcell
   call nc_check(nf90_def_var(ncFileID,name='pfts1d_wtxy', xtype=nf90_double, &
                 dimids=(/ pftDimID /), varid=VarID),&
                 'nc_write_model_atts', 'pfts1d_wtxy def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'pft weight relative to corresponding gridcell'), &
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

      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(varname), xtype=progvar(ivar)%xtype, &
                    dimids = mydimids(1:myndims), varid=VarID),&
                    'nc_write_model_atts', trim(string1)//' def_var' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', trim(progvar(ivar)%long_name)), &
           'nc_write_model_atts', trim(string1)//' put_att long_name' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'DART_kind', trim(progvar(ivar)%kind_string)), &
           'nc_write_model_atts', trim(string1)//' put_att dart_kind' )
      call nc_check(nf90_put_att(ncFileID, VarID, 'units', trim(progvar(ivar)%units)), &
           'nc_write_model_atts', trim(string1)//' put_att units' )

      ! Preserve the original missing_value/_FillValue code.

      if (  progvar(ivar)%xtype == NF90_INT ) then
         call nc_check(nf90_put_att(ncFileID, VarID, 'missing_value', progvar(ivar)%spvalINT), &
              'nc_write_model_atts', trim(string1)//' put_att missing_value' )
         call nc_check(nf90_put_att(ncFileID, VarID, '_FillValue',  progvar(ivar)%spvalINT), &
              'nc_write_model_atts', trim(string1)//' put_att _FillValue' )

      elseif (  progvar(ivar)%xtype == NF90_FLOAT ) then
         call nc_check(nf90_put_att(ncFileID, VarID, 'missing_value', progvar(ivar)%spvalR4), &
              'nc_write_model_atts', trim(string1)//' put_att missing_value' )
         call nc_check(nf90_put_att(ncFileID, VarID, '_FillValue',  progvar(ivar)%spvalR4), &
              'nc_write_model_atts', trim(string1)//' put_att _FillValue' )

      elseif (  progvar(ivar)%xtype == NF90_DOUBLE ) then
         call nc_check(nf90_put_att(ncFileID, VarID, 'missing_value', progvar(ivar)%spvalR8), &
              'nc_write_model_atts', trim(string1)//' put_att missing_value' )
         call nc_check(nf90_put_att(ncFileID, VarID, '_FillValue',  progvar(ivar)%spvalR8), &
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

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
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

      where(dimIDs == CopyDimID) mystart = copyindex
      where(dimIDs == CopyDimID) mycount = 1
      where(dimIDs == TimeDimID) mystart = timeindex
      where(dimIDs == TimeDimID) mycount = 1

      if ((debug > 7) .and. do_output()) then
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

      else

         write(string1,*)'do not know how to handle CLM variables with more than 2 dimensions'
         write(string2,*)trim(progvar(ivar)%varname),'has shape', progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
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



subroutine pert_model_state(state, pert_state, interf_provided)
!------------------------------------------------------------------
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

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

integer :: i
logical, save :: random_seq_init = .false.

if ( .not. module_initialized ) call static_init_model

call error_handler(E_ERR,'pert_model_state', &
                  'CLM cannot be started from a single vector', &
                  source, revision, revdate, &
                  text2='see comments in clm/model_mod.f90::pert_model_state()')

interf_provided = .true.

! Initialize my random number sequence
if(.not. random_seq_init) then
   call init_random_seq(random_seq, my_task_id())
   random_seq_init = .true.
endif

! This does keep the compiler error messages down, but this
! section of code must never be reached.
do i=1,size(state)
   pert_state(i) = random_gaussian(random_seq, state(i), &
                                   model_perturbation_amplitude)
enddo

end subroutine pert_model_state



subroutine ens_mean_for_model(filter_ens_mean)
!------------------------------------------------------------------
! If needed by the model interface, this is the current mean
! for all state vector items across all ensembles.

real(r8), intent(in) :: filter_ens_mean(:)

if ( .not. module_initialized ) call static_init_model

ens_mean = filter_ens_mean

end subroutine ens_mean_for_model



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



subroutine restart_file_to_sv(filename, state_vector, restart_time)
!------------------------------------------------------------------
! Reads the current time and state variables from a clm restart
! file and packs them into a dart state vector. This better happen
! in the same fashion as the metadata arrays are built.

character(len=*), intent(in)    :: filename
real(r8),         intent(inout) :: state_vector(:)
type(time_type),  intent(out)   :: restart_time

! temp space to hold data while we are reading it
integer  :: i, j, ni, nj, ivar, indx, numsnowlevels
integer,  allocatable, dimension(:)         :: snlsno, cols1d_ityplun
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: varname
integer :: VarID, ncNdims, dimlen
integer :: ncid
character(len=256) :: myerrorstring
integer :: fixme = 7000   ! intentionally declaring here so value
                          ! persists across calls.

integer, parameter :: LAKE = 3

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

! Check that the input file exists ...

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'restart_file_to_sv',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncid), &
              'restart_file_to_sv','open '//trim(filename))

restart_time = get_state_time(ncid)

if (do_output()) call print_time(restart_time,'time in restart file '//trim(filename))
if (do_output()) call print_date(restart_time,'date in restart file '//trim(filename))

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

! read number of snow layers

allocate(snlsno(Ncolumn))
call nc_check(nf90_inq_varid(ncid,'SNLSNO', VarID), 'restart_file_to_sv', 'inq_varid SNLSNO')
call nc_check(nf90_get_var(ncid, VarID, snlsno), 'restart_file_to_sv', 'get_var SNLSNO')

! Read in the column landunit type. Knowing what columns are lakes is useful.

allocate(cols1d_ityplun(Ncolumn))
call nc_check(nf90_inq_varid(ncid,'cols1d_ityplun', VarID), 'restart_file_to_sv', 'inq_varid cols1d_ityplun')
call nc_check(nf90_get_var(ncid, VarID, cols1d_ityplun), 'restart_file_to_sv', 'get_var cols1d_ityplun')

! Start counting and filling the state vector one item at a time,
! repacking the Nd arrays into a single 1d list of numbers.

do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   myerrorstring = trim(filename)//' '//trim(varname)

   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            'restart_file_to_sv', 'inq_varid '//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
            'restart_file_to_sv', 'inquire '//trim(myerrorstring))

   ! Check the rank of the variable

   if ( ncNdims /= progvar(ivar)%numdims ) then
      write(string1, *) 'netCDF rank of '//trim(varname)//' does not match derived type knowledge'
      write(string2, *) 'netCDF rank is ',ncNdims,' expected ',progvar(ivar)%numdims
      call error_handler(E_ERR,'restart_file_to_sv', string1, &
                        source,revision,revdate,text2=string2)
   endif

   ! Check the shape of the variable

   do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'restart_file_to_sv', string1)

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(myerrorstring),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         call error_handler(E_ERR,'restart_file_to_sv',string1,source,revision,revdate)
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

      ! TJH FIXME DEBUG TEST
      ! Trying to force a gridpoint that will be affected to have
      ! MISSING values so we can test the block that aborts the update
      ! if some or any of the ensemble members have a missing value.

 !    write(*,*)trim(progvar(ivar)%varname),'replacing ',data_1d_array(fixme),'with',MISSING_R8, 'element',fixme
 !    data_1d_array(fixme) = MISSING_R8
 !    fixme = fixme + 1

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
            do i = 1, Nlevsno - numsnowlevels  ! loop over layers
               data_2d_array(i,j) = MISSING_R8
            enddo
         enddo

      elseif ( (trim(progvar(ivar)%dimnames(1)) == 'levtot') .and. &
               (trim(progvar(ivar)%dimnames(2)) == 'column') ) then

         do j = 1, nj  ! loop over columns
            numsnowlevels = abs(snlsno(j))
            do i = 1, Nlevsno - numsnowlevels  ! loop over layers
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

   else
      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'restart_file_to_sv', string1, &
                        source,revision,revdate)
   endif

   indx = indx - 1
   if ( indx /= progvar(ivar)%indexN ) then
      write(string1, *)'Variable '//trim(varname)//' filled wrong.'
      write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',indx
      call error_handler(E_ERR,'restart_file_to_sv', string1, &
                        source,revision,revdate,text2=string2)
   endif

enddo

call nc_check(nf90_close(ncid),'restart_file_to_sv','close '//trim(filename))
ncid = 0

deallocate(snlsno,cols1d_ityplun)

end subroutine restart_file_to_sv



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

if (do_output()) &
    call print_time(file_time,'time of restart file '//trim(filename))
if (do_output()) &
    call print_date(file_time,'date of restart file '//trim(filename))

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

   if (trim(varname) == 'frac_sno') then
      ! frac_sno (snow cover fraction) is a diagnosed field.
      ! Simply update SNOWDP (snow depth) and CLM recalculates frac_sno.
      ! Actually updating frac_sno causes restart problems for CLM.
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



subroutine model_interpolate(x, location, obs_kind, interp_val, istatus)

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

real(r8),            intent(in)  :: x(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_kind
real(r8),            intent(out) :: interp_val
integer,             intent(out) :: istatus

! Local storage

real(r8), dimension(3) :: loc_array
real(r8) :: llon, llat, lheight
real(r8) :: interp_val_2
integer  :: istatus_2

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

interp_val   = MISSING_R8     ! the DART bad value flag
interp_val_2 = MISSING_R8     ! the DART bad value flag
istatus = 99                ! unknown error

! Get the individual locations values

loc_array = get_location(location)
llon      = loc_array(1)
llat      = loc_array(2)
lheight   = loc_array(3)

if ((debug > 6) .and. do_output()) print *, 'requesting interpolation at ', llon, llat, lheight

! FIXME may be better to check the %maxlevels and kick the interpolation to the
! appropriate routine based on that ... or check the dimnames for the
! vertical coordinate  ...

if (obs_kind == KIND_SOIL_TEMPERATURE) then
   call get_grid_vertval(x, location, 'T_SOISNO',  interp_val, istatus )

elseif (obs_kind == KIND_SOIL_MOISTURE) then
   ! TJH FIXME - actually ROLAND FIXME
   ! This is terrible ... the COSMOS operator wants m3/m3 ... CLM is kg/m2
   call get_grid_vertval(x, location, 'H2OSOI_LIQ',interp_val  , istatus   )
   call get_grid_vertval(x, location, 'H2OSOI_ICE',interp_val_2, istatus_2 )
   if ((istatus == 0) .and. (istatus_2 == 0)) then
      interp_val = interp_val + interp_val_2
   else
      interp_val = MISSING_R8
      istatus = 6
   endif

elseif (obs_kind == KIND_LIQUID_WATER ) then
   call get_grid_vertval(x, location, 'H2OSOI_LIQ',interp_val, istatus )
elseif (obs_kind == KIND_ICE ) then
   call get_grid_vertval(x, location, 'H2OSOI_ICE',interp_val, istatus )
elseif (obs_kind == KIND_SNOWCOVER_FRAC ) then
   call compute_gridcell_value(x, location, 'frac_sno', interp_val, istatus)
elseif (obs_kind == KIND_LEAF_CARBON ) then
   call compute_gridcell_value(x, location, 'leafc',    interp_val, istatus)
elseif (obs_kind == KIND_WATER_TABLE_DEPTH ) then
   call compute_gridcell_value(x, location, 'ZWT',    interp_val, istatus)
elseif (obs_kind == KIND_SNOW_THICKNESS ) then
   write(string1,*)'model_interpolate for DZSNO not written yet.'
   call error_handler(E_ERR,'compute_gridcell_value',string1,source,revision,revdate)
   istatus = 5
elseif ((obs_kind == KIND_GEOPOTENTIAL_HEIGHT) .and. vert_is_level(location)) then
   if (nint(lheight) > nlevgrnd) then
      interp_val = MISSING_R8
      istatus = 1
   else
      interp_val = LEVGRND(nint(lheight))
      istatus = 0
   endif
else
   write(string1,*)'model_interpolate not written for (integer) kind ',obs_kind
   call error_handler(E_ERR,'compute_gridcell_value',string1,source,revision,revdate)
   istatus = 5
endif

if ((debug > 6) .and. do_output()) write(*,*)'interp_val ',interp_val

end subroutine model_interpolate


!------------------------------------------------------------------


subroutine compute_gridcell_value(x, location, varstring, interp_val, istatus)
!
! Each gridcell may contain values for several land units, each land unit may contain
! several columns, each column may contain several pft's. BUT this routine never
! aggregates across multiple pft's. So, each gridcell value
! is an area-weighted value of an unknown number of column-based quantities.

! Passed variables

real(r8),            intent(in)  :: x(:)         ! state vector
type(location_type), intent(in)  :: location     ! location somewhere in a grid cell
character(len=*),    intent(in)  :: varstring    ! frac_sno, leafc
real(r8),            intent(out) :: interp_val   ! area-weighted result
integer,             intent(out) :: istatus      ! error code (0 == good)

! Local storage

integer  :: ivar, index1, indexN, indexi, counter
integer  :: gridloni,gridlatj
real(r8) :: loc_lat, loc_lon
real(r8) :: total, total_area
real(r8), dimension(1) :: loninds,latinds
real(r8), dimension(3) :: loc

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
ivar   = findVarindex(varstring, 'compute_gridcell_value')
index1 = progvar(ivar)%index1 ! in the DART state vector, start looking here
indexN = progvar(ivar)%indexN ! in the DART state vector, stop  looking here

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

if ((debug > 6) .and. do_output()) then
   write(*,*)'compute_gridcell_value:targetlon, lon, lon index is ',loc_lon,LON(gridloni),gridloni
   write(*,*)'compute_gridcell_value:targetlat, lat, lat index is ',loc_lat,LAT(gridlatj),gridlatj
endif

! If there is no vertical component, the problem is greatly simplified.
! Simply area-weight an average of all pieces in the grid cell.

counter    = 0
total      = 0.0_r8      ! temp storage for state vector
total_area = 0.0_r8      ! temp storage for area
ELEMENTS : do indexi = index1, indexN

   if (   lonixy(indexi) /=  gridloni ) cycle ELEMENTS
   if (   latjxy(indexi) /=  gridlatj ) cycle ELEMENTS
   if (        x(indexi) == MISSING_R8) cycle ELEMENTS
   if ( landarea(indexi) ==   0.0_r8  ) cycle ELEMENTS

   counter    = counter    + 1
   total      = total      + x(indexi)*landarea(indexi)
   total_area = total_area +           landarea(indexi)

   if ((debug > 6) .and. do_output()) then
      write(*,*)
      write(*,*)'gridcell location match',counter,'at statevector index',indexi
      write(*,*)'statevector value is (',x(indexi),')'
      write(*,*)'area is              (',landarea(indexi),')'
      write(*,*)'LON index is         (',lonixy(indexi),')'
      write(*,*)'LAT index is         (',latjxy(indexi),')'
      write(*,*)'closest LON is       (',LON(gridloni),')'
      write(*,*)'closest LAT is       (',LAT(gridlatj),')'
      write(*,*)'closest lev is       (',levels(indexi),')'
   endif

enddo ELEMENTS

if (total_area /= 0.0_r8) then ! All good.
   interp_val = total/total_area
   istatus    = 0
else
   if ((debug > 0) .and. do_output()) then
      write(string1, *)'Variable '//trim(varstring)//' had no viable data'
      write(string2, *)'at gridcell ilon/jlat = (',gridloni,',',gridlatj,')'
      write(string3, *)'obs lon/lat = (',loc_lon,',',loc_lat,')'
      call error_handler(E_MSG,'compute_gridcell_value', string1, &
                     source, revision, revdate, text2=string2,text3=string3)
   endif
endif

! Print more information for the really curious
if ((debug > 6) .and. do_output()) then
   write(string1,*)'counter, total, total_area', counter, total, total_area
   write(string2,*)'interp_val, istatus', interp_val, istatus
   call error_handler(E_MSG,'compute_gridcell_value', string1, &
                     source, revision, revdate, text2=string2)
endif

end subroutine compute_gridcell_value


!------------------------------------------------------------------


subroutine get_grid_vertval(x, location, varstring, interp_val, istatus)
!
! Calculate the expected vertical value for the gridcell.
! Each gridcell value is an area-weighted value of an unknown number of
! column-based quantities.

! Passed variables

real(r8),            intent(in)  :: x(:)         ! state vector
type(location_type), intent(in)  :: location     ! location somewhere in a grid cell
character(len=*),    intent(in)  :: varstring    ! T_SOISNO, H2OSOI_LIQ, H2OSOI_ICE
real(r8),            intent(out) :: interp_val   ! area-weighted result
integer,             intent(out) :: istatus      ! error code (0 == good)

! Local storage

integer  :: ivar, index1, indexN, indexi, counter1, counter2
integer  :: gridloni,gridlatj
real(r8), dimension(3) :: loc
real(r8) :: loc_lat, loc_lon, loc_lev
real(r8) :: value_below, value_above, total_area
real(r8) :: depthbelow, depthabove
real(r8) :: topwght, botwght
real(r8), dimension(1) :: loninds,latinds

real(r8), allocatable, dimension(:)   :: above, below
real(r8), allocatable, dimension(:,:) :: myarea

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
ivar   = findVarindex(varstring, 'get_grid_vertval')
index1 = progvar(ivar)%index1 ! in the DART state vector, start looking here
indexN = progvar(ivar)%indexN ! in the DART state vector, stop  looking here

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

if ((debug > 6) .and. do_output()) then
   write(*,*)'get_grid_vertval:targetlon, lon, lon index, level is ',loc_lon,LON(gridloni),gridloni,loc_lev
   write(*,*)'get_grid_vertval:targetlat, lat, lat index, level is ',loc_lat,LAT(gridlatj),gridlatj,loc_lev
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

if ((debug > 6) .and. do_output()) then
   write(*,*)'get_grid_vertval:depthbelow ',depthbelow,'>= loc_lev',loc_lev,'>= depthabove',depthabove
endif

! Determine how many elements can contribute to the gridcell value.
! There are multiple column-based contributors, each column has a
! separate area-based weight. There are multiple levels.
! I believe I have to keep track of all of them to sort out how to
! calculate the gridcell value at a particular depth.

counter1 = 0
counter2 = 0
GRIDCELL : do indexi = index1, indexN
   if ( lonixy(indexi) /=  gridloni )  cycle GRIDCELL
   if ( latjxy(indexi) /=  gridlatj )  cycle GRIDCELL
   if (      x(indexi) == MISSING_R8)  cycle GRIDCELL

   if (levels(indexi) == depthabove) counter1 = counter1 + 1
   if (levels(indexi) == depthbelow) counter2 = counter2 + 1

enddo GRIDCELL

if ( (counter1+counter2) == 0 ) then
   if ((debug > 0) .and. do_output()) then
      write(string1, *)'statevector variable '//trim(varstring)//' had no viable data'
      write(string2, *)'at gridcell lon/lat = (',gridloni,',',gridlatj,')'
      write(string3, *)'obs lon/lat/lev (',loc_lon,',',loc_lat,',',loc_lev,')'
      call error_handler(E_MSG,'get_grid_vertval', string1, &
                     source, revision, revdate, text2=string2,text3=string3)
   endif
   return
endif

allocate(above(counter1),below(counter2),myarea(max(counter1,counter2),2))
above  = MISSING_R8
below  = MISSING_R8
myarea = 0.0_r8

counter1 = 0
counter2 = 0
ELEMENTS : do indexi = index1, indexN

   if ( lonixy(indexi) /=  gridloni )  cycle ELEMENTS
   if ( latjxy(indexi) /=  gridlatj )  cycle ELEMENTS
   if (      x(indexi) == MISSING_R8)  cycle ELEMENTS

!  write(*,*)'level ',indexi,' is ',levels(indexi),' location depth is ',loc_lev

   if (levels(indexi)     == depthabove) then
      counter1            = counter1 + 1
      above( counter1)    =        x(indexi)
      myarea(counter1,1)  = landarea(indexi)
   endif
   if (levels(indexi)     == depthbelow) then
      counter2            = counter2 + 1
      below( counter2)    =        x(indexi)
      myarea(counter2,2)  = landarea(indexi)
   endif

   if ((levels(indexi) /= depthabove) .and. &
       (levels(indexi) /= depthbelow)) then
      cycle ELEMENTS
   endif

   if ((debug > 6) .and. do_output()) then
   write(*,*)
   write(*,*)'gridcell location match at statevector index',indexi
   write(*,*)'statevector value is (',x(indexi),')'
   write(*,*)'area is          (',landarea(indexi),')'
   write(*,*)'LON index is     (',lonixy(indexi),')'
   write(*,*)'LAT index is     (',latjxy(indexi),')'
   write(*,*)'gridcell LON is  (',LON(gridloni),')'
   write(*,*)'gridcell LAT is  (',LAT(gridlatj),')'
   write(*,*)'depth        is  (',levels(indexi),')'
   endif

enddo ELEMENTS

! could arise if the above or below was 'missing' ... but the mate was not.

if ( counter1 /= counter2 ) then
   write(string1, *)'Variable '//trim(varstring)//' has peculiar interpolation problems.'
   write(string2, *)'uneven number of values "above" and "below"'
   write(string3, *)'counter1 == ',counter1,' /= ',counter2,' == counter2'
   call error_handler(E_MSG,'get_grid_vertval', string1, &
                  source, revision, revdate, text2=string2,text3=string3)
   return
endif

! Determine the value for the level above the depth of interest.

total_area = sum(myarea(1:counter1,1))

if ( total_area /= 0.0_r8 ) then
   ! normalize the area-based weights
   myarea(1:counter1,1) = myarea(1:counter1,1) / total_area
   value_above = sum(above(1:counter1) * myarea(1:counter1,1))
else
   write(string1, *)'Variable '//trim(varstring)//' had no viable data above'
   write(string2, *)'at gridcell lon/lat/lev = (',gridloni,',',gridlatj,',',depthabove,')'
   write(string3, *)'obs lon/lat/lev (',loc_lon,',',loc_lat,',',loc_lev,')'
   call error_handler(E_ERR,'get_grid_vertval', string1, &
                  source, revision, revdate, text2=string2,text3=string3)
endif

! Determine the value for the level below the depth of interest.

total_area = sum(myarea(1:counter2,2))

if ( total_area /= 0.0_r8 ) then
   ! normalize the area-based weights
   myarea(1:counter2,2) = myarea(1:counter2,2) / total_area
   value_below = sum(below(1:counter2) * myarea(1:counter2,2))
else
   write(string1, *)'Variable '//trim(varstring)//' had no viable data below'
   write(string2, *)'at gridcell lon/lat/lev = (',gridloni,',',gridlatj,',',depthbelow,')'
   write(string3, *)'obs lon/lat/lev (',loc_lon,',',loc_lat,',',loc_lev,')'
   call error_handler(E_ERR,'get_grid_vertval', string1, &
                  source, revision, revdate, text2=string2,text3=string3)
endif

if (depthbelow == depthabove) then
   topwght = 1.0_r8
   botwght = 0.0_r8
else
   topwght = (depthbelow - loc_lev) / (depthbelow - depthabove)
   botwght = (loc_lev - depthabove) / (depthbelow - depthabove)
endif

interp_val = value_above*topwght + value_below*botwght
istatus    = 0

deallocate(above, below, myarea)

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
integer, optional,        intent(in)  :: ncid

integer :: i,ii, VarID
real(r8), allocatable, dimension(:) :: org_array

if ( .not. module_initialized ) call static_init_model

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

! Replace the DART fill value with the original value and apply any clamping.
! Get the 'original' variable from the netcdf file.

if (present(ncid)) then

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
integer, optional,        intent(in)  :: ncid

integer :: i,j,ii, VarID
real(r8), allocatable, dimension(:,:) :: org_array

if ( .not. module_initialized ) call static_init_model

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

! Replace the DART fill value with the original value and apply any clamping.
! Get the 'original' variable from the netcdf file if need be.

if (present(ncid)) then

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
integer, optional, intent(out)   :: lonatm, latatm
integer, optional, intent(out)   :: lonrof, latrof

integer :: dimid

if ( .not. module_initialized ) call static_init_model

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

if ( .not. module_initialized ) call static_init_model

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

if ( .not. module_initialized ) call static_init_model

if (ncid == 0) then
   call nc_check(nf90_open(trim(fname), nf90_nowrite, ncid), &
               'get_sparse_dims','open '//trim(fname))
endif

! get dimid for 'gridcell' and then get value ...

call nc_check(nf90_inq_dimid(ncid, 'gridcell', dimid), &
            'get_sparse_dims','inq_dimid gridcell '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Ngridcell), &
            'get_sparse_dims','inquire_dimension gridcell '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'landunit', dimid), &
            'get_sparse_dims','inq_dimid landunit '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nlandunit), &
            'get_sparse_dims','inquire_dimension landunit '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'column', dimid), &
            'get_sparse_dims','inq_dimid column '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Ncolumn), &
            'get_sparse_dims','inquire_dimension column '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'pft', dimid), &
            'get_sparse_dims','inq_dimid pft '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Npft), &
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
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nlevlak), &
            'get_sparse_dims','inquire_dimension levlak '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'levtot', dimid), &
            'get_sparse_dims','inq_dimid levtot '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nlevtot), &
            'get_sparse_dims','inquire_dimension levtot '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'numrad', dimid), &
            'get_sparse_dims','inq_dimid numrad '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nnumrad), &
            'get_sparse_dims','inquire_dimension numrad '//trim(fname))

! levsno is presently required, but I can envision a domain/experiment that
! will not have snow levels. How this relates to variables dimensioned 'levtot'
! is unclear. For that reason, levsno is presently required.

istatus = nf90_inq_dimid(ncid, 'levsno', dimid)
if (istatus == nf90_noerr) then
   call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nlevsno), &
               'get_sparse_dims','inquire_dimension levsno '//trim(fname))
endif

! levsno1, rtmlon, rtmlat are optional ... so it is not a fatal error if they are not present.

istatus = nf90_inq_dimid(ncid, 'levsno1', dimid)
if (istatus == nf90_noerr) then
   call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nlevsno1), &
               'get_sparse_dims','inquire_dimension levsno1 '//trim(fname))
endif

istatus = nf90_inq_dimid(ncid, 'rtmlon', dimid)
if (istatus == nf90_noerr) then
   call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nrtmlon), &
               'get_sparse_dims','inquire_dimension rtmlon '//trim(fname))
endif

istatus = nf90_inq_dimid(ncid, 'rtmlat', dimid)
if (istatus == nf90_noerr) then
   call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nrtmlat), &
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
   write(logfileunit,*)'Ngridcell = ',Ngridcell
   write(logfileunit,*)'Nlandunit = ',Nlandunit
   write(logfileunit,*)'Ncolumn   = ',Ncolumn
   write(logfileunit,*)'Npft      = ',Npft
   write(logfileunit,*)'Nlevgrnd  = ',Nlevgrnd
   write(logfileunit,*)'Nlevlak   = ',Nlevlak
   write(logfileunit,*)'Nlevsno   = ',Nlevsno
   write(logfileunit,*)'Nlevsno1  = ',Nlevsno1
   write(logfileunit,*)'Nlevtot   = ',Nlevtot
   write(logfileunit,*)'Nnumrad   = ',Nnumrad
   write(logfileunit,*)'Nrtmlon   = ',Nrtmlon
   write(logfileunit,*)'Nrtmlat   = ',Nrtmlat
   write(     *     ,*)
   write(     *     ,*)'get_sparse_dims output follows:'
   write(     *     ,*)'Ngridcell = ',Ngridcell
   write(     *     ,*)'Nlandunit = ',Nlandunit
   write(     *     ,*)'Ncolumn   = ',Ncolumn
   write(     *     ,*)'Npft      = ',Npft
   write(     *     ,*)'Nlevgrnd  = ',Nlevgrnd
   write(     *     ,*)'Nlevlak   = ',Nlevlak
   write(     *     ,*)'Nlevsno   = ',Nlevsno
   write(     *     ,*)'Nlevsno1  = ',Nlevsno1
   write(     *     ,*)'Nlevtot   = ',Nlevtot
   write(     *     ,*)'Nnumrad   = ',Nnumrad
   write(     *     ,*)'Nrtmlon   = ',Nrtmlon
   write(     *     ,*)'Nrtmlat   = ',Nrtmlat
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

if ( .not. module_initialized ) call static_init_model

if (ncid == 0) then
   call nc_check(nf90_open(trim(fname), nf90_nowrite, ncid), &
               'get_sparse_geog','open '//trim(fname))
endif

! Make sure the variables are the right size ...
! by comparing agains the size of the variable ...

if ( Ngridcell < 0 ) then
   write(string1,*)'Unable to read the number of gridcells.'
   call error_handler(E_ERR,'get_sparse_geog',string1,source,revision,revdate)
endif

if ( Nlandunit < 0 ) then
   write(string1,*)'Unable to read the number of land units.'
   call error_handler(E_ERR,'get_sparse_geog',string1,source,revision,revdate)
endif

if ( Ncolumn < 0 ) then
   write(string1,*)'Unable to read the number of columns.'
   call error_handler(E_ERR,'get_sparse_geog',string1,source,revision,revdate)
endif

if ( Npft < 0 ) then
   write(string1,*)'Unable to read the number of pfts.'
   call error_handler(E_ERR,'get_sparse_geog',string1,source,revision,revdate)
endif

! Read the netcdf file data

call nc_check(nf90_inq_varid(ncid, 'grid1d_ixy', VarID), 'get_sparse_geog', 'inq_varid grid1d_ixy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   grid1d_ixy), 'get_sparse_geog',   'get_var grid1d_ixy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'grid1d_jxy', VarID), 'get_sparse_geog', 'inq_varid grid1d_jxy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   grid1d_jxy), 'get_sparse_geog',   'get_var grid1d_jxy '//trim(clm_restart_filename))


call nc_check(nf90_inq_varid(ncid, 'land1d_ixy', VarID), 'get_sparse_geog', 'inq_varid land1d_ixy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   land1d_ixy), 'get_sparse_geog',   'get_var land1d_ixy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'land1d_jxy', VarID), 'get_sparse_geog', 'inq_varid land1d_jxy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   land1d_jxy), 'get_sparse_geog',   'get_var land1d_jxy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'land1d_wtxy', VarID), 'get_sparse_geog', 'inq_varid land1d_wtxy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   land1d_wtxy), 'get_sparse_geog',   'get_var land1d_wtxy '//trim(clm_restart_filename))


call nc_check(nf90_inq_varid(ncid, 'cols1d_ixy', VarID), 'get_sparse_geog', 'inq_varid cols1d_ixy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   cols1d_ixy), 'get_sparse_geog',   'get_var cols1d_ixy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'cols1d_jxy', VarID), 'get_sparse_geog', 'inq_varid cols1d_jxy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   cols1d_jxy), 'get_sparse_geog',   'get_var cols1d_jxy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'cols1d_wtxy', VarID), 'get_sparse_geog', 'inq_varid cols1d_wtxy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   cols1d_wtxy), 'get_sparse_geog',   'get_var cols1d_wtxy '//trim(clm_restart_filename))


call nc_check(nf90_inq_varid(ncid, 'pfts1d_ixy', VarID), 'get_sparse_geog', 'inq_varid pfts1d_ixy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   pfts1d_ixy), 'get_sparse_geog',   'get_var pfts1d_ixy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'pfts1d_jxy', VarID), 'get_sparse_geog', 'inq_varid pfts1d_jxy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   pfts1d_jxy), 'get_sparse_geog',   'get_var pfts1d_jxy '//trim(clm_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'pfts1d_wtxy', VarID), 'get_sparse_geog', 'inq_varid pfts1d_wtxy '//trim(clm_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   pfts1d_wtxy), 'get_sparse_geog',   'get_var pfts1d_wtxy '//trim(clm_restart_filename))

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
   write(logfileunit,*)'grid1d_ixy  range ',minval(grid1d_ixy),maxval(grid1d_ixy)
   write(logfileunit,*)'grid1d_jxy  range ',minval(grid1d_jxy),maxval(grid1d_jxy)

   write(logfileunit,*)'land1d_ixy  range ',minval(land1d_ixy),maxval(land1d_ixy)
   write(logfileunit,*)'land1d_jxy  range ',minval(land1d_jxy),maxval(land1d_jxy)
   write(logfileunit,*)'land1d_wtxy range ',minval(land1d_wtxy),maxval(land1d_wtxy)

   write(logfileunit,*)'cols1d_ixy  range ',minval(cols1d_ixy),maxval(cols1d_ixy)
   write(logfileunit,*)'cols1d_jxy  range ',minval(cols1d_jxy),maxval(cols1d_jxy)
   write(logfileunit,*)'cols1d_wtxy range ',minval(cols1d_wtxy),maxval(cols1d_wtxy)

   write(logfileunit,*)'pfts1d_ixy  range ',minval(pfts1d_ixy),maxval(pfts1d_ixy)
   write(logfileunit,*)'pfts1d_jxy  range ',minval(pfts1d_jxy),maxval(pfts1d_jxy)
   write(logfileunit,*)'pfts1d_wtxy range ',minval(pfts1d_wtxy),maxval(pfts1d_wtxy)
   if (nlevsno > 0) write(logfileunit,*)'zsno        range ',minval(zsno),maxval(zsno)

   write(     *     ,*)
   write(     *     ,*)'Raw lat/lon information as read ...'
   write(     *     ,*)'grid1d_ixy  range ',minval(grid1d_ixy),maxval(grid1d_ixy)
   write(     *     ,*)'grid1d_jxy  range ',minval(grid1d_jxy),maxval(grid1d_jxy)

   write(     *     ,*)'land1d_ixy  range ',minval(land1d_ixy),maxval(land1d_ixy)
   write(     *     ,*)'land1d_jxy  range ',minval(land1d_jxy),maxval(land1d_jxy)
   write(     *     ,*)'land1d_wtxy range ',minval(land1d_wtxy),maxval(land1d_wtxy)

   write(     *     ,*)'cols1d_ixy  range ',minval(cols1d_ixy),maxval(cols1d_ixy)
   write(     *     ,*)'cols1d_jxy  range ',minval(cols1d_jxy),maxval(cols1d_jxy)
   write(     *     ,*)'cols1d_wtxy range ',minval(cols1d_wtxy),maxval(cols1d_wtxy)

   write(     *     ,*)'pfts1d_ixy  range ',minval(pfts1d_ixy),maxval(pfts1d_ixy)
   write(     *     ,*)'pfts1d_jxy  range ',minval(pfts1d_jxy),maxval(pfts1d_jxy)
   write(     *     ,*)'pfts1d_wtxy range ',minval(pfts1d_wtxy),maxval(pfts1d_wtxy)
   if (nlevsno > 0) write(     *     ,*)'zsno        range ',minval(zsno),maxval(zsno)

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


function get_state_time_fname(filename)
!------------------------------------------------------------------
! the static_init_model ensures that the clm namelists are read.
!
type(time_type) :: get_state_time_fname
character(len=*), intent(in) :: filename

integer         :: ncid

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'get_state_time_fname',string1,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'get_state_time_fname', 'open '//trim(filename))

get_state_time_fname = get_state_time_ncid(ncid)

call nc_check(nf90_close(ncid),'get_state_time_fname', 'close '//trim(filename))

end function get_state_time_fname


!------------------------------------------------------------------


function set_model_time_step()
!------------------------------------------------------------------
! This defines the window used for assimilation.
! all observations +/- half this timestep are assimilated.

type(time_type) :: set_model_time_step

if ( .not. module_initialized ) call static_init_model

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


subroutine verify_state_variables( state_variables, ncid, filename, ngood, table )

character(len=*), dimension(:),   intent(in)  :: state_variables
integer,                          intent(in)  :: ncid
character(len=*),                 intent(in)  :: filename
integer,                          intent(out) :: ngood
character(len=*), dimension(:,:), intent(out) :: table

integer :: nrows, ncols, i, VarID
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr

if ( .not. module_initialized ) call static_init_model

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
      string1 = 'model_nml:clm_state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Make sure variable exists in netCDF file

   write(string1,'(''there is no variable '',a,'' in '',a)') trim(varname), trim(filename)
   call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
                 'verify_state_variables', trim(string1))

   ! Make sure DART kind is valid

   if( get_raw_obs_kind_index(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector

   if ((debug > 3) .and. do_output()) then
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

! Check to see if zsno is part of the requested variables

end subroutine verify_state_variables


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


subroutine define_var_dims(ivar, ncid, memberdimid, unlimiteddimid, ndims, dimids)

integer,               intent(in)  :: ivar, ncid, memberdimid, unlimiteddimid
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: dimids

integer :: i, mydimid

ndims = 0

do i = 1,progvar(ivar)%numdims

   call nc_check(nf90_inq_dimid(ncid=ncid, name=progvar(ivar)%dimnames(i), dimid=mydimid), &
                           'define_var_dims','inq_dimid '//trim(progvar(ivar)%dimnames(i)))

   ndims         = ndims + 1
   dimids(ndims) = mydimid

enddo

ndims = ndims + 1
dimids(ndims) = memberdimid
ndims = ndims + 1
dimids(ndims) = unlimitedDimid

if ((debug > 8) .and. do_output()) then

   write(logfileunit,*)
   write(logfileunit,*)'define_var_dims knowledge'
   write(logfileunit,*)trim(progvar(ivar)%varname),' has dimnames ', &
                       progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
   write(logfileunit,*)' thus dimids ',dimids(1:ndims)
   write(     *     ,*)
   write(     *     ,*)'define_var_dims knowledge'
   write(     *     ,*)trim(progvar(ivar)%varname),' has dimnames ', &
                       progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
   write(     *     ,*)' thus dimids ',dimids(1:ndims)

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


subroutine find_gridcell_Npft(varstring)

! In order to exercise some of the routines, it is necessary to know
! which gridcells have multiple land units
! which gridcells have multiple columns
! which gridcells have multiple PFTs
!
! This routine simply tells me which gridcells are 'interesting'

character(len=*), intent(in) :: varstring    ! T_SOISNO, H2OSOI_LIQ

! Local storage

integer                       :: ivar, indexi, i, j
integer, dimension(nlon,nlat) :: countmat

if ( .not. module_initialized ) call static_init_model

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

write(*,*)'exploring '//trim(varstring)

do j = 1,nlat
do i = 1,nlon

   if ( countmat(i,j) > 1) then
      write(*,'(''gridcell'',2(1x,i8),'' has '',i6,'' lon/lat'',2(1x,f12.7))') &
                          i,j,countmat(i,j),LON(i),LAT(j)
   endif

enddo
enddo

end subroutine find_gridcell_Npft



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

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs,dimlens
integer  :: VarID, numdims, xtype, io1, io2
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8

integer,  allocatable, dimension(:) :: intarray
real(r4), allocatable, dimension(:) :: r4array

if ( .not. module_initialized ) call static_init_model

call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), 'get_var_1d', 'inq_varid')
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, xtype=xtype), 'get_var_1d', 'inquire_variable')
call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlens(1)), 'get_var_1d', 'inquire_dimension')

if ((numdims /= 1) .or. (size(var1d) /= dimlens(1)) ) then
   write(string1,*) trim(varname)//' is not the expected shape/length of ', size(var1d)
   call error_handler(E_ERR,'get_var_1d',string1,source,revision,revdate)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1)))
   call nc_check(nf90_get_var(ncid, VarID, intarray), 'get_var_1d', 'get_var')
   var1d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1,source,revision,revdate)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1,source,revision,revdate)
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
   call nc_check(nf90_get_var(ncid, VarID, r4array), 'get_var_1d', 'get_var')
   var1d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1,source,revision,revdate)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1,source,revision,revdate)
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

   call nc_check(nf90_get_var(ncid, VarID, var1d), 'get_var_1d', 'get_var')

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var1d == spvalR8) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1,source,revision,revdate)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var1d == spvalR8) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1,source,revision,revdate)
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

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens
integer  :: VarID, numdims, xtype, io1, io2
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8

integer,  allocatable, dimension(:,:) :: intarray
real(r4), allocatable, dimension(:,:) :: r4array

if ( .not. module_initialized ) call static_init_model

call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), 'get_var_2d', 'inq_varid')
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, xtype=xtype), 'get_var_2d', 'inquire_variable')

if ( (numdims /= 2)  ) then
   write(string1,*) trim(varname)//' is not a 2D variable as expected.'
   call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate)
endif

call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlens(1)), 'get_var_2d', 'inquire_dimension 1')
call nc_check(nf90_inquire_dimension(ncid, dimIDs(2), len=dimlens(2)), 'get_var_2d', 'inquire_dimension 2')

if ( (size(var2d,1) /= dimlens(1))  .or. (size(var2d,2) /= dimlens(2)) ) then
   write(string1,*) trim(varname)//' has shape ', dimlens(1), dimlens(2)
   write(string2,*) 'which is not the expected shape of ', size(var2d,1),size(var2d,2)
   call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate,text2=string2)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1),dimlens(2)))
   call nc_check(nf90_get_var(ncid, VarID, intarray), 'get_var_2d', 'get_var')
   var2d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1,source,revision,revdate)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1,source,revision,revdate)
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
   call nc_check(nf90_get_var(ncid, VarID, r4array), 'get_var_2d', 'get_var')
   var2d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1,source,revision,revdate)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1,source,revision,revdate)
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

   call nc_check(nf90_get_var(ncid, VarID, var2d), 'get_var_2d', 'get_var')

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var2d == spvalR8) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1,source,revision,revdate)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var2d == spvalR8) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1,source,revision,revdate)
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



function get_model_time()
type(time_type) :: get_model_time

if ( .not. module_initialized ) call static_init_model

get_model_time = model_time

end function get_model_time



function findVarindex(varstring, caller)
character(len=*), intent(in) :: varstring
character(len=*), intent(in) :: caller
integer                      :: findVarindex

integer :: i

findVarindex = -1

! Skip to the right variable
VARTYPES : do i = 1,nfields
    findVarindex = i
    if ( trim(progvar(i)%varname) == varstring) exit VARTYPES
enddo VARTYPES

if (findVarindex < 1) then
   write(string1,*) trim(caller)//' cannot find '//trim(varstring)
   call error_handler(E_ERR,'findVarindex',string1,source,revision,revdate)
endif

end function findVarindex

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
