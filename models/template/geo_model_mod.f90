! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> this version of the model_mod is probably more appropriate to start with
!> if you are porting a geophysical model with locations specified in lat/lon
!> and some kind of vertical.  it assumes a gridded model.

!> THIS FILE HAS NOT BEEN CONVERTED to read netcdf directly by
!> defining a domain - a list of netcdf variables per input file.
!> THIS WILL BE UPDATED SOON.

module model_mod

! This is the interface between the model model and DART.

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, digits12, SECPERDAY, MISSING_R8,          &
                             rad2deg, deg2rad, PI
use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
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
                             loc_get_close_obs => get_close_obs,               &
                             get_close_state => get_close_state


use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper,                    &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, file_exist, find_textfile_dims,        &
                             file_to_text, close_file

use     obs_kind_mod, only : obstypelength,        &
                             get_index_for_quantity,  &
                             get_name_for_quantity

use mpi_utilities_mod, only: my_task_id

use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use     dart_model_mod, only: get_model_nLons, get_model_nLats, get_model_nAlts, &
                             get_nSpecies, get_nSpeciesTotal, get_nIons,     &
                             get_nSpeciesAll, decode_model_indices

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
          get_close_state_init,   &
          get_close_state

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.

public :: get_gridsize,                &
          restart_file_to_statevector, &
          statevector_to_restart_file, &
          get_model_restart_dirname,    &
          get_base_time,               &
          get_state_time

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=256) :: string1, string2
logical, save :: module_initialized = .false.

! Storage for a random sequence for perturbing a single initial state

type(random_seq_type) :: random_seq

! things which can/should be in the model_nml

integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 60
real(r8)           :: model_perturbation_amplitude = 0.2
logical            :: output_state_vector = .true.
integer            :: debug = 0   ! turn up for more and more debug messages
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: model_restart_dirname = 'model_restartdir'

namelist /model_nml/  &
   model_restart_dirname,        &
   output_state_vector,         &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   calendar,                    &
   debug

!------------------------------------------------------------------
!
!  The DART state vector may consist of things like:  
!
!  U    long_name = "X-WIND COMPONENT"      float   U(TIME, ALT, LAT, XE) 
!  V    long_name = "Y-WIND COMPONENT"      float   V(TIME, ALT, YE, LON)
!  W    long_name = "Z-WIND COMPONENT"      float   W(TIME, ZE, LAT, LON)
!  TH   long_name = "POTENTIAL TEMPERATURE" float  TH(TIME, ALT, LAT, LON)
!  DBZ  long_name = "RADAR REFLECTIVITY"    float DBZ(TIME, ALT, LAT, LON)
!  WZ   long_name = "VERTICAL VORTICITY"    float  WZ(TIME, ALT, LAT, LON)
!  PI   long_name = "PERT. EXNER"	    float  PI(TIME, ALT, LAT, LON)
!  QV   long_name = "VAPOR MIXING RATIO"    float  QV(TIME, ALT, LAT, LON)
!  QC   long_name = "CLOUD MIXING RATIO"    float  QC(TIME, ALT, LAT, LON)
!  QR   long_name = "RAIN MIXING RATIO"     float  QR(TIME, ALT, LAT, LON)
!  QI   long_name = "ICE MIXING RATIO"      float  QI(TIME, ALT, LAT, LON)
!  QS   long_name = "SNOW MIXING RATIO"     float  QS(TIME, ALT, LAT, LON)
!  QH   long_name = "GRAUPEL MIXING RATIO"  float  QH(TIME, ALT, LAT, LON)
!
!  The variables in the model restart file that are used to create the 
!  DART state vector are specified in the input.nml:model_vars_nml namelist.
!
!------------------------------------------------------------------

integer, parameter :: max_state_variables = 80
integer, parameter :: num_state_table_columns = 2
character(len=NF90_MAX_NAME) :: model_state_variables(max_state_variables * num_state_table_columns ) = ' '
character(len=NF90_MAX_NAME) :: variable_table(max_state_variables, num_state_table_columns )

namelist /model_vars_nml/ model_state_variables

integer :: nfields

! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname       ! crazy species name
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=NF90_MAX_NAME) :: storder
   character(len=NF90_MAX_NAME) :: model_varname  ! NDensityS, IDensityS, ...
   integer :: model_dim                           ! dimension defining species
   integer :: model_index                         ! 'iSpecies' or u,v,w ...
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens ! nlons, nlats, nalts [, nspecies]
   integer :: posdef
   integer :: numdims
   integer :: varsize     ! prod(dimlens(1:numdims))
   integer :: index1      ! location in dart state vector of first occurrence
   integer :: indexN      ! location in dart state vector of last  occurrence
   integer :: dart_kind
   character(len=obstypelength) :: kind_string
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

! These are statically defined in ModSize.f90 ...
! nAlts  is the one and only number of altitudes ... no block-dependence
! nLons, nLats are the number of lons/lats PER block
! the number of blocks comes from UAM.in 

integer :: nLons, nLats, nAlts

! "... keep in mind that if the model resolution is 5 deg latitude,
!  the model will actually go from -87.5 to 87.5 latitude
! (even though you specify -90 to 90 in the UAM.in file), 
! since the latitudes/longitudes are at cell centers, 
! while the edges are at the boundaries." -- Aaron Ridley

integer  :: NgridLon=-1, NgridLat=-1, NgridAlt=-1    ! scalar grid counts
integer  :: nBlocksLon=-1, nBlocksLat=-1             ! number of blocks along each dim
real(r8) :: LatStart=MISSING_R8, LatEnd=MISSING_R8, LonStart=MISSING_R8
integer  :: nSpeciesTotal=-1, nSpecies=-1, nIons=-1, nSpeciesAll=-1

! scalar grid positions

real(r8), allocatable :: LON(:)   ! longitude centers
real(r8), allocatable :: LAT(:)   ! latitude  centers
real(r8), allocatable :: ALT(:)   ! vertical level centers

integer               :: model_size      ! the state vector length
type(time_type)       :: model_time      ! valid time of the model state
type(time_type)       :: model_timestep  ! smallest time to adv model
real(r8), allocatable :: ens_mean(:)     ! may be needed for forward ops

! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.

logical :: print_timestamps = .false.
integer :: print_every_Nth  = 10000

integer, parameter :: nGhost = 2   ! number of ghost cells on all edges

!------------------------------------------------------------------
! The model restart manager namelist variables
!------------------------------------------------------------------

character(len= 64) :: ew_boundary_type, ns_boundary_type

INTERFACE vector_to_prog_var
      MODULE PROCEDURE vector_to_1d_prog_var
      MODULE PROCEDURE vector_to_2d_prog_var
      MODULE PROCEDURE vector_to_3d_prog_var
      MODULE PROCEDURE vector_to_4d_prog_var
END INTERFACE

INTERFACE get_base_time
      MODULE PROCEDURE get_base_time_ncid
      MODULE PROCEDURE get_base_time_fname
END INTERFACE

INTERFACE get_index_range
      MODULE PROCEDURE get_index_range_string
      MODULE PROCEDURE get_index_range_int
END INTERFACE

contains

!==================================================================
! All the REQUIRED interfaces come first - just by convention.
!==================================================================


function get_model_size()
!------------------------------------------------------------------
! Done - TJH.
! Returns the size of the model as an integer. 
! Required for all applications.

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size



subroutine adv_1step(x, time)
!------------------------------------------------------------------
! Done - TJH.
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

! FIXME: put an error handler call here - we cannot advance the model
! this way and it would be an error if filter called it.

end subroutine adv_1step



subroutine get_state_meta_data(index_in, location, var_type)
!------------------------------------------------------------------
! Done - JLA.
! given an index into the state vector, return its location and
! if given, the var kind.   despite the name, var_type is a generic
! kind, like those in obs_kind/obs_kind_mod.f90, starting with QTY_

integer, intent(in)            :: index_in
type(location_type)            :: location
integer, optional, intent(out) :: var_type
  
! Local variables

integer :: lat_index, lon_index, alt_index
integer :: n, nf, myindx, remainder, remainder2

if ( .not. module_initialized ) call static_init_model

! Find out which of the 3D fields index_in is part of 
nf     = -1
  
FindIndex : do n = 1,nfields
   if( (progvar(n)%index1 <= index_in) .and. (index_in <= progvar(n)%indexN) ) then
      nf = n
      myindx = index_in - progvar(n)%index1 + 1
      exit FindIndex
   endif
enddo FindIndex 
  
if( myindx == -1 ) then
   write(string1,*) 'Problem, cannot find base_offset, index_in is: ', index_in
   call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
endif
  
alt_index = 1 + (myindx - 1) / (NgridLon * NgridLat)
remainder = myindx - (alt_index-1) * NgridLon * NgridLat
lat_index = 1 + (remainder - 1) / NgridLon
remainder2 = remainder - (lat_index - 1) * NgridLon
lon_index = remainder2

location = set_location(LON(lon_index), LAT(lat_index), ALT(alt_index), VERTISHEIGHT)
  
if (present(var_type)) then
   var_type = progvar(nf)%dart_kind
endif

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
! the dart_model_mod module.

! Local variables - all the important ones have module scope

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname
character(len=obstypelength)       :: kind_string
integer :: ncid, VarID, numdims, dimlen, varsize
integer :: iunit, io, ivar, i, index1, indexN
integer :: ss, dd
integer :: nDimensions, nVariables, nAttributes, unlimitedDimID, TimeDimID
logical :: shapeok

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
call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
if (do_output()) write(logfileunit, nml=model_nml)
if (do_output()) write(     *     , nml=model_nml)

! Get the model variables in a restricted scope setting.

nLons         = get_model_nLons()
nLats         = get_model_nLats()
nAlts         = get_model_nAlts()
nSpecies      = get_nSpecies()
nSpeciesTotal = get_nSpeciesTotal()
nIons         = get_nIons()
nSpeciesAll   = get_nSpeciesAll()

if ((debug > 0) .and.  do_output() ) then
   write(*,*)
   write(*,*)'nLons         is ',nLons
   write(*,*)'nLats         is ',nLats
   write(*,*)'nAlts         is ',nAlts
   write(*,*)'nSpecies      is ',nSpecies
   write(*,*)'nSpeciesTotal is ',nSpeciesTotal
   write(*,*)'nIons         is ',nIons
   write(*,*)'nSpeciesAll   is ',nSpeciesAll
endif

! Read the model variable list to populate DART state vector
! Once parsed, the values will be recorded for posterity
call find_namelist_in_file('input.nml', 'model_vars_nml', iunit)
read(iunit, nml = model_vars_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_vars_nml')

!---------------------------------------------------------------
! Set the time step ... causes model namelists to be read.
! Ensures model_timestep is multiple of 'dynamics_timestep'

call set_calendar_type( calendar )   ! comes from model_mod_nml

model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)

!---------------------------------------------------------------
! 1) get grid dimensions
! 2) allocate space for the grids 
! 3) read them from the block restart files, could be stretched ...

call get_grid_info(NgridLon, NgridLat, NgridAlt, nBlocksLon, nBlocksLat, &
                   LatStart, LatEnd, LonStart)

if( debug  > 0 ) then
   write(*,*)'grid dims are ',NgridLon,NgridLat,NgridAlt
endif

allocate( LON( NgridLon ))
allocate( LAT( NgridLat ))
allocate( ALT( NgridAlt ))

call get_grid(model_restart_dirname, nBlocksLon, nBlocksLat, &
              nLons, nLats, nAlts, LON, LAT, ALT )
              
!---------------------------------------------------------------
! Compile the list of model variables to use in the creation
! of the DART state vector. Required to determine model_size.
!
! Verify all variables are in the model restart file
!
! Compute the offsets into the state vector for the start of each
! different variable type. Requires reading shapes from the model
! restart file. As long as TIME is the LAST dimension, we're OK.
!
! Record the extent of the data type in the state vector.

call verify_state_variables( model_state_variables, ncid, model_restart_dirname, &
                             nfields, variable_table )

index1  = 1;
indexN  = 0;

do ivar = 1, nfields 

   varname                   = trim(variable_table(ivar,1))
   kind_string               = trim(variable_table(ivar,2))
   progvar(ivar)%varname     = varname
   progvar(ivar)%kind_string = kind_string
   progvar(ivar)%dart_kind   = get_index_for_quantity( progvar(ivar)%kind_string ) 
   progvar(ivar)%dimlens     = 0

   ! I would really like decode_model_indices to set the following (on a per-variable basis)
   ! progvar(ivar)%storder
   ! progvar(ivar)%numdims
   ! progvar(ivar)%dimlens

   call decode_model_indices( varname, progvar(ivar)%model_varname, progvar(ivar)%model_dim, &
                             progvar(ivar)%model_index, progvar(ivar)%long_name, &
                             progvar(ivar)%units)

   varsize = NgridLon * NgridLat * NgridAlt

   progvar(ivar)%storder     = 'xyz3d'
   progvar(ivar)%numdims     = 3
   progvar(ivar)%dimlens(1:progvar(ivar)%numdims) = (/ NgridLon, NgridLat, NgridAlt /)
   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1 
   index1                    = index1 + varsize      ! sets up for next variable

   if ( debug > 0 ) then
      write(logfileunit,*)
      write(logfileunit,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(logfileunit,*) ' storage      ',trim(progvar(ivar)%storder)
      write(logfileunit,*) ' long_name    ',trim(progvar(ivar)%long_name)
      write(logfileunit,*) ' units        ',trim(progvar(ivar)%units)
      write(logfileunit,*) ' numdims      ',progvar(ivar)%numdims
      write(logfileunit,*) ' dimlens      ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      write(logfileunit,*) ' varsize      ',progvar(ivar)%varsize
      write(logfileunit,*) ' index1       ',progvar(ivar)%index1
      write(logfileunit,*) ' indexN       ',progvar(ivar)%indexN
      write(logfileunit,*) ' dart_kind    ',progvar(ivar)%dart_kind
      write(logfileunit,*) ' kind_string  ',trim(progvar(ivar)%kind_string)
      write(logfileunit,*) ' model_varname ',trim(progvar(ivar)%model_varname)
      write(logfileunit,*) ' model_dim     ',progvar(ivar)%model_dim
      write(logfileunit,*) ' model_index   ',progvar(ivar)%model_index

      write(     *     ,*)
      write(     *     ,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(     *     ,*) ' storage      ',trim(progvar(ivar)%storder)
      write(     *     ,*) ' long_name    ',trim(progvar(ivar)%long_name)
      write(     *     ,*) ' units        ',trim(progvar(ivar)%units)
      write(     *     ,*) ' numdims      ',progvar(ivar)%numdims
      write(     *     ,*) ' dimlens      ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      write(     *     ,*) ' varsize      ',progvar(ivar)%varsize
      write(     *     ,*) ' index1       ',progvar(ivar)%index1
      write(     *     ,*) ' indexN       ',progvar(ivar)%indexN
      write(     *     ,*) ' dart_kind    ',progvar(ivar)%dart_kind
      write(     *     ,*) ' kind_string  ',trim(progvar(ivar)%kind_string)
      write(     *     ,*) ' model_varname ',trim(progvar(ivar)%model_varname)
      write(     *     ,*) ' model_dim     ',progvar(ivar)%model_dim
      write(     *     ,*) ' model_index   ',progvar(ivar)%model_index
   endif

enddo

model_size = progvar(nfields)%indexN

if ( debug > 0 ) then
  write(logfileunit,'("grid: NgridLon, NgridLat, NgridAlt =",3(1x,i5))') NgridLon, NgridLat, NgridAlt
  write(     *     ,'("grid: NgridLon, NgridLat, NgridAlt =",3(1x,i5))') NgridLon, NgridLat, NgridAlt
  write(logfileunit, *)'model_size = ', model_size
  write(     *     , *)'model_size = ', model_size
endif

allocate( ens_mean(model_size) )

end subroutine static_init_model



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

if (allocated(LON)) deallocate(LON, LAT, ALT)

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

! FIXME: put an error handler call here - we cannot initialize the model
! this way and it would be an error if filter called it.

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

! FIXME: put an error handler call here - we cannot initialize the model
! this way and it would be an error if filter called it.

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
integer :: NLONDimID
integer :: NLATDimID
integer :: NALTDimID

! for the prognostic variables
integer :: ivar, VarID

!----------------------------------------------------------------------
! variables for the namelist output
!----------------------------------------------------------------------

character(len=129), allocatable, dimension(:) :: textblock
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen
logical :: has_model_namelist

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

call nc_check(nf90_inq_dimid(ncid=ncFileID, name='NMLlinelen', dimid=LineLenDimID), &
                           'nc_write_model_atts','inq_dimid NMLlinelen')
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
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'model' ), &
           'nc_write_model_atts', 'model put '//trim(filename))

!-------------------------------------------------------------------------------
! Determine shape of most important namelist
!-------------------------------------------------------------------------------

call find_textfile_dims('model_vars.nml', nlines, linelen)
if (nlines > 0) then
  has_model_namelist = .true.
else
  has_model_namelist = .false.
endif

if (debug > 0)    print *, 'model namelist: nlines, linelen = ', nlines, linelen
  
if (has_model_namelist) then 
   allocate(textblock(nlines))
   textblock = ''

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nlines', &
                 len = nlines, dimid = nlinesDimID), &
                 'nc_write_model_atts', 'def_dim nlines ')

   call nc_check(nf90_def_var(ncFileID,name='model_in', xtype=nf90_char,    &
                 dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID), &
                 'nc_write_model_atts', 'def_var model_in')
   call nc_check(nf90_put_att(ncFileID, nmlVarID, 'long_name',       &
                 'contents of model_in namelist'), 'nc_write_model_atts', 'put_att model_in')

endif

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
   
   call nc_check(nf90_def_dim(ncid=ncFileID, name='LON', &
          len = NgridLon, dimid = NLONDimID),'nc_write_model_atts', 'LON def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='LAT', &
          len = NgridLat, dimid = NLATDimID),'nc_write_model_atts', 'LAT def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='ALT', &
          len = NgridAlt, dimid = NALTDimID),'nc_write_model_atts', 'ALT def_dim '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Coordinate Variables and the Attributes
   !----------------------------------------------------------------------------

   ! Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='LON', xtype=nf90_real, &
                 dimids=NLONDimID, varid=VarID),&
                 'nc_write_model_atts', 'LON def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VarID, 'type', 'x1d'),  &
                 'nc_write_model_atts', 'LON type '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'grid longitudes'), &
                 'nc_write_model_atts', 'LON long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'X'),  &
                 'nc_write_model_atts', 'LON cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'LON units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'LON valid_range '//trim(filename))

   ! Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='LAT', xtype=nf90_real, &
                 dimids=NLATDimID, varid=VarID),&
                 'nc_write_model_atts', 'LAT def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VarID, 'type', 'y1d'),  &
                 'nc_write_model_atts', 'LAT type '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'grid latitudes'), &
                 'nc_write_model_atts', 'LAT long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Y'),   &
                 'nc_write_model_atts', 'LAT cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_north'),  &
                 'nc_write_model_atts', 'LAT units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'LAT valid_range '//trim(filename))

   ! Grid Altitudes
   call nc_check(nf90_def_var(ncFileID,name='ALT',xtype=nf90_real, &
                 dimids=NALTDimID,varid=VarID), &
                 'nc_write_model_atts', 'ALT def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VarID, 'type', 'z1d'),  &
                 'nc_write_model_atts', 'ALT type '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'grid altitudes'), &
                 'nc_write_model_atts', 'ALT long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'meters'),  &
                 'nc_write_model_atts', 'ALT units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'positive', 'up'),  &
                 'nc_write_model_atts', 'ALT units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Z'),   &
                 'nc_write_model_atts', 'ALT cartesian_axis '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------

   do ivar=1, nfields

      varname = trim(progvar(ivar)%varname)
      string1 = trim(filename)//' '//trim(varname)

      ! match shape of the variable to the dimension IDs

      call define_var_dims(progvar(ivar), myndims, mydimids, MemberDimID, unlimitedDimID, &
                      NLONDimID, NLATDimID, NALTDimID, NgridLon, NgridLat, NgridAlt)

      ! define the variable and set the attributes

      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(varname), xtype=nf90_real, &
                    dimids = mydimids(1:myndims), varid=VarID),&
                    'nc_write_model_atts', trim(string1)//' def_var' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'type', trim(progvar(ivar)%storder)), &
           'nc_write_model_atts', trim(string1)//' put_att storage type' )

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

   call nc_check(nf90_enddef(ncfileID), 'prognostic enddef '//trim(filename))

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables
   !----------------------------------------------------------------------------

   call nc_check(NF90_inq_varid(ncFileID, 'LON', VarID), &
                 'nc_write_model_atts', 'LON inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, LON ), &
                'nc_write_model_atts', 'LON put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'LAT', VarID), &
                 'nc_write_model_atts', 'LAT inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, LAT ), &
                'nc_write_model_atts', 'LAT put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'ALT', VarID), &
                 'nc_write_model_atts', 'ALT inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, ALT ), &
                'nc_write_model_atts', 'ALT put_var '//trim(filename))

endif

!-------------------------------------------------------------------------------
! Fill the variables we can
!-------------------------------------------------------------------------------

if (has_model_namelist) then
   call file_to_text('model_vars.nml', textblock)
   call nc_check(nf90_put_var(ncFileID, nmlVarID, textblock ), &
                 'nc_write_model_atts', 'put_var nmlVarID')
   deallocate(textblock)
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts



function nc_write_model_vars( ncFileID, state_vec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
!
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
real(r8), allocatable, dimension(:,:,:,:) :: data_4d_array

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

      call nc_check(nf90_inquire_variable(ncFileID,VarId,dimids=dimIDs,ndims=ncNdims), &
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

      if ( debug > 1 ) then
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
         call vector_to_prog_var(state_vec, progvar(ivar), data_1d_array)
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
         call vector_to_prog_var(state_vec, progvar(ivar), data_2d_array)
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
         call vector_to_prog_var(state_vec, progvar(ivar), data_3d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_3d_array)

      elseif ( progvar(ivar)%numdims == 4) then

         if ( ncNdims /= 6 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 6 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_4d_array( progvar(ivar)%dimlens(1), &
                                 progvar(ivar)%dimlens(2), &
                                 progvar(ivar)%dimlens(3), &
                                 progvar(ivar)%dimlens(4)))
         call vector_to_prog_var(state_vec, progvar(ivar), data_4d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_4d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_4d_array)

      else

         ! FIXME put an error message here

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

integer :: i
logical, save :: random_seq_init = .false.

if ( .not. module_initialized ) call static_init_model

interf_provided = .true.

! Initialize my random number sequence
if(.not. random_seq_init) then
   call init_random_seq(random_seq, my_task_id())
   random_seq_init = .true.
endif

! add some uncertainty to each ...
do i=1,size(state)
   pert_state(i) = random_gaussian(random_seq, state(i), &
                                   model_perturbation_amplitude)
enddo

end subroutine pert_model_state



subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, &
                         obs_loc, obs_kind, num_close, close_ind, dist)
!------------------------------------------------------------------

! Given a DART location (referred to as "base") and a set of candidate
! locations & kinds (obs, obs_kind), returns the subset close to the
! "base", their indices, and their distances to the "base" ...

! For vertical distance computations, general philosophy is to convert all
! vertical coordinates to a common coordinate. This coordinate type is defined
! in the namelist with the variable "vert_localization_coord".

type(get_close_type),              intent(in)    :: gc
type(location_type),               intent(inout) :: base_obs_loc
integer,                           intent(in)    :: base_obs_kind
type(location_type), dimension(:), intent(inout) :: obs_loc
integer,             dimension(:), intent(in)    :: obs_kind
integer,                           intent(out)   :: num_close
integer,             dimension(:), intent(out)   :: close_ind
real(r8),            dimension(:), intent(out)   :: dist

integer                :: t_ind, istatus1, istatus2, k
integer                :: base_which, local_obs_which
real(r8), dimension(3) :: base_array, local_obs_array
type(location_type)    :: local_obs_loc

! Initialize variables to missing status

num_close = 0
close_ind = -99
dist      = 1.0e9   !something big and positive (far away)
istatus1  = 0
istatus2  = 0

! Convert base_obs vertical coordinate to requested vertical coordinate if necessary

base_array = get_location(base_obs_loc)
base_which = nint(query_location(base_obs_loc))

! fixme ... 
if (.not. horiz_dist_only) then
!  if (base_which /= wrf%dom(1)%vert_coord) then
!     call vert_interpolate(ens_mean, base_obs_loc, base_obs_kind, istatus1)
!  elseif (base_array(3) == MISSING_R8) then
!     istatus1 = 1
!  endif
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
 !fixme       if (local_obs_which /= wrf%dom(1)%vert_coord) then
 !fixme           call vert_interpolate(ens_mean, local_obs_loc, obs_kind(t_ind), istatus2)
            ! Store the "new" location into the original full local array
            obs_loc(t_ind) = local_obs_loc
 !fixme        endif
      endif

      ! Compute distance - set distance to a very large value if vert coordinate is missing
      ! or vert_interpolate returned error (istatus2=1)
      local_obs_array = get_location(local_obs_loc)
      if (( (.not. horiz_dist_only)             .and. &
            (local_obs_array(3) == MISSING_R8)) .or.  &
            (istatus2 == 1)                   ) then
            dist(k) = 1.0e9
      else
            dist(k) = get_dist(base_obs_loc, local_obs_loc, base_obs_kind, obs_kind(t_ind))
      endif

   enddo
endif

end subroutine get_close_obs


!==================================================================
! The remaining PUBLIC interfaces come next
!==================================================================



subroutine get_gridsize(num_LON, num_xe, num_LAT, num_ye, num_ALT, num_ze )
 integer, intent(out) :: num_LON, num_LAT, num_ALT
 integer, intent(out) :: num_xe, num_ye, num_ze
!------------------------------------------------------------------
! public utility routine.

if ( .not. module_initialized ) call static_init_model

 num_LON = NgridLon
 num_LAT = NgridLat
 num_ALT = NgridAlt

end subroutine get_gridsize



! FIXME:
!  this routine needs:  
!  1.  a base dirname for the restart files.
!  they will have the format 'dirname/bNNNN.rst'  where NNNN has
!  leading 0s and is the block number.   blocks start in the
!  southwest corner of the lat/lon grid and go west first, then
!  north and end in the northeast corner.   (assuming var 'dirname')
!  the other info is in 'dirname/header.rst'
!
!  2. the overall grid size, lon/lat/alt when you've read in all
!  the blocks.  (nGridLon, nGridLat, nGridAlt, will compute totalVarSize)
!
!  3. the number of blocks in Lon and Lat (nBlocksLon, nBlocksLat)
!
!  4. the number of lon/lats in a single grid block  (nLons, nLats, nAlts)
!
!  5. the number of neutral species (and probably a mapping between
!  the species number and the variable name)  (nSpeciesTotal, nSpecies)
!
!  6. the number of ion species (ditto - numbers <-> names) (nIons)
!
! we assume that the 'UseTopography' flag is false - that all columns
! have the same altitude arrays.  this is true on earth but not on
! other planets.
! 
!  in addition to reading in the state data, it fills Longitude,
!  Latitude, and Altitude arrays with the grid spacing.  this grid
!  is orthogonal and rectangular but can have irregular spacing along
!  any or all of the three dimensions.
!

subroutine restart_file_to_statevector(dirname, state_vector, model_time)
!------------------------------------------------------------------
! Reads the current time and state variables from a model restart
! file and packs them into a dart state vector.

character(len=*), intent(in)  :: dirname 
real(r8),         intent(out) :: state_vector(:)
type(time_type),  intent(out) :: model_time

integer :: ivar, i
character(len=NF90_MAX_NAME) :: varname

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

! this is going to have to loop over all the blocks, both to get
! the data values and to get the full grid spacings.

model_time = get_state_time(dirname)

if (do_output()) &
    call print_time(model_time,'time in restart file '//trim(dirname)//'/header.rst')
if (do_output()) &
    call print_date(model_time,'date in restart file '//trim(dirname)//'/header.rst')

! sort the required fields into the order they exist in the
! binary restart files and fill in the state vector as you
! read each field.  when this routine returns all the data has
! been read.

call get_data(dirname, state_vector)


end subroutine restart_file_to_statevector



subroutine statevector_to_restart_file(state_vector, dirname, statedate)
!------------------------------------------------------------------
! Writes the current time and state variables from a dart state
! vector (1d array) into a model netcdf restart file.
!
real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: dirname 
type(time_type),  intent(in) :: statedate


integer :: ivar 
character(len=NF90_MAX_NAME) :: varname
character(len=128) :: dirnameout

if ( .not. module_initialized ) call static_init_model

print *, 'in statevector_to_restart_file, debug, nfields = ', debug, nfields

! sort the required fields into the order they exist in the
! binary restart files and write out the state vector data
! field by field.  when this routine returns all the data has
! been written.

dirnameout = trim(dirname) // '.out'
call put_data(dirname, dirnameout, state_vector)

! FIXME:
! write out model_time to a text file?
   call print_time(model_time)

 
if (do_output()) &
    call print_time(model_time,'time in restart file '//trim(dirname)//'/header.rst')
if (do_output()) &
    call print_date(model_time,'date in restart file '//trim(dirname)//'/header.rst')

end subroutine statevector_to_restart_file



!==================================================================
! The remaining interfaces come last
!==================================================================



!############################################################################
!
!     ##################################################################
!     ######                                                      ######
!     ######            SUBROUTINE MODEL_INTERPOLATE              ######
!     ######                                                      ######
!     ##################################################################
!
!
!     PURPOSE:
!
!     For a given lat, lon, and height, interpolate the correct state value
!     to that location for the filter from the model state vectors
!
!############################################################################
!     Modified from NCOMMAS model_mod.f90 by Tim Hoar and Jeff Anderson
!     Revision Date: August 2011
!
!     Original routine information
!
!     Author:  Tim Hoar, Ted Mansell, Lou Wicker
!
!     Creation Date:  August 2010
!     
!     Variables needed to be stored in the MODEL_MODULE data structure
!
!       LON   = 1D array storing the local grid center coords (degrees)
!       LAT   = 1D array storing the local grid center coords (degrees)
!       ALT   = 1D array storing the local grid center coords (meters)
!
!       ERROR codes:
!
!       ISTATUS = 99:  general error in case something terrible goes wrong...
!       ISTATUS = 15:  dont know what to do with vertical coord supplied
!       ISTATUS = 16:  longitude illegal
!       ISTATUS = 17:  latitude illegal
!       ISTATUS = 18:  altitude illegal
!       
!
!############################################################################


subroutine model_interpolate(x, location, obs_type, interp_val, istatus)

! Passed variables

  real(r8),            intent(in)  :: x(:)
  type(location_type), intent(in)  :: location
  integer,             intent(in)  :: obs_type
  real(r8),            intent(out) :: interp_val
  integer,             intent(out) :: istatus

! Local storage

  real(r8)         :: loc_array(3), llon, llat, lheight, lon_fract, lat_fract, alt_fract
  integer          :: base_offset, end_offset, blon(2), blat(2), balt(2), i, j, k, ier
  real(r8)         :: cube(2, 2, 2), square(2, 2), line(2)

  IF ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the 
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

  interp_val = MISSING_R8     ! the DART bad value flag
  istatus = 99                ! unknown error

! Get the individual locations values

  loc_array = get_location(location)
  llon      = loc_array(1)
  llat      = loc_array(2)
  lheight   = loc_array(3)

  IF (debug > 2) print *, 'requesting interpolation at ', llon, llat, lheight

! Only height for vertical location type is supported at this point
  IF(.not. vert_is_height(location) ) THEN 
     istatus = 15
     return
  ENDIF

! Find the start and end offsets for this field in the state vector x(:)

  call get_index_range(obs_type, base_offset, end_offset)

  IF (debug > 2) print *, 'base offset now ', base_offset

! Need to find bounding center indices and fractional offsets for lat, lon, alt
call find_lon_bounds(llon, blon(1), blon(2), lon_fract, ier)
if(ier /= 0) then
   istatus = 16
   return
endif

call find_lat_or_alt_bounds(llat, NgridLat, LAT, blat(1), blat(2), lat_fract, ier)
if(ier /= 0) then
   istatus = 17
   return
endif

call find_lat_or_alt_bounds(lheight, NgridAlt, ALT, balt(1), balt(2), alt_fract, ier)
if(ier /= 0) then
   istatus = 18
   return
endif

! Get the grid values for the first 
do i = 1, 2
   do j = 1, 2
      do k = 1, 2
         cube(i, j, k) = get_grid_value(base_offset, blon(i), blat(j), balt(k), x)
      end do
   end do
end do

! Interpolate to the given altitude
do i = 1, 2
   do j = 1, 2
      square(i, j) = (1 - alt_fract) * cube(i, j, 1) + alt_fract * cube(i, j, 2)
   end do
end do

! Interpolate to the given latitude
do i = 1, 2
   line(i) = (1 - lat_fract) * square(i, 1) + lat_fract * square(i, 2)
end do

! Interpolate to the given longitude
interp_val = (1 - lon_fract) * line(1) + lon_fract * line(2)

! All good.
  istatus = 0

return
end subroutine model_interpolate


!------------------------------------------------------------------

function get_grid_value(base_offset, ilon, ilat, ialt, x)

real(r8)             :: get_grid_value
integer, intent(in)  :: base_offset, ilon, ilat, ialt
real(r8), intent(in) :: x(:)

! Returns the value for the given lon,lat,alt point in the field that 
! starts at offset base_offset

integer :: offset

offset = (ilon - 1) + (ilat - 1) * NgridLon + (ialt - 1) * (NgridLon * NgridLat)
get_grid_value = x(base_offset + offset)

end function get_grid_value

!------------------------------------------------------------------


subroutine find_lon_bounds(llon, lower, upper, fract, ier)

! Finds position of a given longitude in an array of longitude grid points and returns
! the index of the lower and upper bounds and the fractional offset. Assumes that the
! first longitude in the list is the smallest and that the largest is less than
! 360 degrees.  ier returns 0 unless there is an error.

real(r8), intent(in)  :: llon
integer,  intent(out) :: lower, upper
real(r8), intent(out) :: fract
integer,  intent(out) :: ier

! For now, assume that the spacing on longitudes is arbitrary.
! Do a silly linear search. Probably not worth any fancier searching unless
! models get to be huge.

integer  :: i
real(r8) :: width

if(llon < 0.0_r8 .or. llon > 360.0_r8) then
   ier = 2
   return
endif

! Look for case where longitude is less than smallest in list
if(llon <= LON(1)) then
   lower = NgridLon
   upper = 1
   width = 360.0_r8 - LON(NgridLon) + LON(1)
   fract = (llon + 360.0_r8 - LON(NgridLon)) / width
   ier = 0
   return
endif

! Look for case where longitude is greater than largest in list
if(llon >= LON(NgridLon)) then
  lower = NgridLon
  upper = 1
  width = 360.0 - LON(NgridLon) + LON(1)
  fract = (llon - LON(NgridLon)) / width
  ier = 0
  return
endif

! Otherwise in the interior
do i = 2, NgridLon
   if(llon < LON(i)) then
      lower = i - 1
      upper = i
      fract = (llon - LON(i-1)) / (LON(i) - LON(i - 1))
      ier = 0
      return
   endif
end do

! Shouldn't ever fall off end of loop
ier = 2

end subroutine find_lon_bounds

!------------------------------------------------------------------

subroutine find_lat_or_alt_bounds(llat, nbounds, bounds, lower, upper, fract, ier)

! Finds position of a given latitude in an array of latitude grid points and returns
! the index of the lower and upper bounds and the fractional offset. Used for both
! latitude and altitude which have similar linear arrays. ier returns 0 unless there
! is an error.

real(r8), intent(in)  :: llat
integer,  intent(in)  :: nbounds
real(r8), intent(in)  :: bounds(nbounds)
integer,  intent(out) :: lower, upper
real(r8), intent(out) :: fract
integer,  intent(out) :: ier

! For now, assume that the spacing on latitudes or altitudes is arbitrary
! Do a silly linear search. Probably not worth any fancier searching unless
! models get to be huge.

integer :: i

if(llat < bounds(1) .or. llat > bounds(nbounds)) then
   ier = 2
   return
endif

do i = 2, nbounds
   if(llat <= bounds(i)) then
      lower = i - 1
      upper = i
      fract = (llat - bounds(i-1)) / (bounds(i) - bounds(i - 1))
      ier = 0
      return
   endif
end do

! Shouldn't ever fall off end of loop
ier = 2

end subroutine find_lat_or_alt_bounds

!------------------------------------------------------------------


subroutine vector_to_1d_prog_var(x, progvar, data_1d_array)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset,
! into a 1d array.
!
real(r8), dimension(:),   intent(in)  :: x
type(progvartype),        intent(in)  :: progvar
real(r8), dimension(:),   intent(out) :: data_1d_array

integer :: i,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar%index1

do i = 1, progvar%dimlens(1)
   data_1d_array(i) = x(ii)
   ii = ii + 1
enddo

ii = ii - 1
if ( ii /= progvar%indexN ) then
   write(string1, *)'Variable '//trim(progvar%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_1d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_1d_prog_var


!------------------------------------------------------------------


subroutine vector_to_2d_prog_var(x, progvar, data_2d_array)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset,
! into a 2d array.
!
real(r8), dimension(:),   intent(in)  :: x
type(progvartype),        intent(in)  :: progvar
real(r8), dimension(:,:), intent(out) :: data_2d_array

integer :: i,j,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar%index1

do j = 1,progvar%dimlens(2)
do i = 1,progvar%dimlens(1)
   data_2d_array(i,j) = x(ii)
   ii = ii + 1
enddo
enddo

ii = ii - 1
if ( ii /= progvar%indexN ) then
   write(string1, *)'Variable '//trim(progvar%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_2d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_2d_prog_var


!------------------------------------------------------------------


subroutine vector_to_3d_prog_var(x, progvar, data_3d_array)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset,
! into a 3d array.
!
real(r8), dimension(:),     intent(in)  :: x
type(progvartype),          intent(in)  :: progvar
real(r8), dimension(:,:,:), intent(out) :: data_3d_array

integer :: i,j,k,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar%index1

do k = 1,progvar%dimlens(3)
do j = 1,progvar%dimlens(2)
do i = 1,progvar%dimlens(1)
   data_3d_array(i,j,k) = x(ii)
   ii = ii + 1
enddo
enddo
enddo

ii = ii - 1
if ( ii /= progvar%indexN ) then
   write(string1, *)'Variable '//trim(progvar%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_3d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_3d_prog_var


!------------------------------------------------------------------


subroutine vector_to_4d_prog_var(x, progvar, data_4d_array)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset,
! into a 4d array.
!
real(r8), dimension(:),       intent(in)  :: x
type(progvartype),            intent(in)  :: progvar
real(r8), dimension(:,:,:,:), intent(out) :: data_4d_array

integer :: i,j,k,m,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar%index1

do m = 1,progvar%dimlens(4)
do k = 1,progvar%dimlens(3)
do j = 1,progvar%dimlens(2)
do i = 1,progvar%dimlens(1)
   data_4d_array(i,j,k,m) = x(ii)
   ii = ii + 1
enddo
enddo
enddo
enddo

ii = ii - 1
if ( ii /= progvar%indexN ) then
   write(string1, *)'Variable '//trim(progvar%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_4d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_4d_prog_var


!------------------------------------------------------------------


subroutine get_grid_info(NgridLon, NgridLat, NgridAlt, &
                nBlocksLon, nBlocksLat, LatStart, LatEnd, LonStart)
!------------------------------------------------------------------
!
! Read the grid dimensions from the restart netcdf file.
!
! The file name comes from module storage ... namelist.

integer,  intent(out) :: NgridLon   ! Number of Longitude centers
integer,  intent(out) :: NgridLat   ! Number of Latitude  centers
integer,  intent(out) :: NgridAlt   ! Number of Vertical grid centers
integer,  intent(out) :: nBlocksLon, nBlocksLat
real(r8), intent(out) :: LatStart, LatEnd, LonStart

integer :: grid_id, dimid
character(len=obstypelength) :: filename = 'UAM.in'

character(len=100) :: cLine  ! iCharLen_ == 100
character(len=128) :: fileloc

integer :: i, iunit, ios

if ( .not. module_initialized ) call static_init_model

! get the ball rolling ...

nBlocksLon = 0
nBlocksLat = 0
LatStart   = 0.0_r8
LatEnd     = 0.0_r8
LonStart   = 0.0_r8

write(fileloc,'(a,''/'',a)') trim(model_restart_dirname),trim(filename)

iunit = open_file(fileloc, action='read')

UAMREAD : do i = 1, 1000000

   read(iunit,'(a)',iostat=ios) cLine

   if (ios /= 0) then
      write(string1,*) 'cannot find #GRID in ',trim(filename)
      call error_handler(E_ERR,'get_grid_info',string1,source,revision,revdate)
   endif

   if (cLine(1:5) .ne. "#GRID") cycle UAMREAD

      nBlocksLon = read_in_int( iunit,'NBlocksLon',filename)
      nBlocksLat = read_in_int( iunit,'NBlocksLat',filename)
      LatStart   = read_in_real(iunit,'LatStart',  filename)
      LatEnd     = read_in_real(iunit,'LatEnd',    filename)
      LonStart   = read_in_real(iunit,'LonStart',  filename)

   exit UAMREAD

enddo UAMREAD

call close_file(iunit)

NgridLon = nBlocksLon * nLons
NgridLat = nBlocksLat * nLats
NgridAlt = nAlts

end subroutine get_grid_info


!------------------------------------------------------------------


subroutine get_grid(dirname, nBlocksLon, nBlocksLat, &
                  nLons, nLats, nAlts, LON, LAT, ALT )
!------------------------------------------------------------------
! open enough of the restart files to read in the lon, lat, alt arrays
!
character(len=*), intent(in) :: dirname
integer, intent(in) :: nBlocksLon ! Number of Longitude blocks
integer, intent(in) :: nBlocksLat ! Number of Latitude  blocks
integer, intent(in) :: NLons      ! Number of Longitude centers per block
integer, intent(in) :: NLats      ! Number of Latitude  centers per block
integer, intent(in) :: NAlts      ! Number of Vertical grid centers

real(r8), dimension( : ), intent(inout) :: LON, LAT, ALT

integer :: nb, offset, iunit, nboff
character(len=128) :: filename
real(r8), allocatable :: temp(:)

! a temp array large enough to hold any of the 
! Lon,Lat or Alt array from a block plus ghost cells
allocate(temp(1-nGhost:max(nLons,nLats,nAlts)+nGhost))

! get the dirname, construct the filenames inside 

! go across the south-most block row picking up all longitudes
do nb = 1, nBlocksLon

   iunit = open_block_file(dirname, nb, 'read')

   read(iunit) temp(1-nGhost:nLons+nGhost)

   offset = (nLons * (nb - 1)) 
   LON(offset+1:offset+nLons) = temp(1:nLons)

   call close_file(iunit)
enddo

! go up west-most block row picking up all latitudes
do nb = 1, nBlocksLat

   nboff = ((nb - 1) * nBlocksLon) + 1
   iunit = open_block_file(dirname, nboff, 'read')

   ! get past lon array and read in lats
   read(iunit) temp(1-nGhost:nLons+nGhost)

   read(iunit) temp(1-nGhost:nLats+nGhost)

   offset = (nLats * (nb - 1)) 
   LAT(offset+1:offset+nLats) = temp(1:nLats)

   call close_file(iunit)
enddo

! this code assumes UseTopography is false - that all columns share
! the same altitude array, so we can read it from the first block.  
! if this is not the case, this code has to change.

iunit = open_block_file(dirname, 1, 'read')

! get past lon and lat arrays and read in alt array
read(iunit) temp(1-nGhost:nLons+nGhost)
read(iunit) temp(1-nGhost:nLats+nGhost)
read(iunit) temp(1-nGhost:nAlts+nGhost)

ALT(1:nAlts) = temp(1:nAlts)

call close_file(iunit)

deallocate(temp)

! convert from radians into degrees
LON = LON * rad2deg
LAT = LAT * rad2deg

if (debug > 4) then
   print *, 'All LONs ', LON
   print *, 'All LATs ', LAT
   print *, 'All ALTs ', ALT
endif

if ( debug > 1 ) then ! A little sanity check
   write(*,*)'LON range ',minval(LON),maxval(LON)
   write(*,*)'LAT range ',minval(LAT),maxval(LAT)
   write(*,*)'ALT range ',minval(ALT),maxval(ALT)
endif

end subroutine get_grid


!------------------------------------------------------------------


function open_block_file(dirname, blocknum, rw)
!------------------------------------------------------------------
! open the requested block number restart file and return the
! open file unit

integer                      :: open_block_file
character(len=*), intent(in) :: dirname
integer,          intent(in) :: blocknum
character(len=*), intent(in) :: rw   ! 'read' or 'readwrite'

character(len=128) :: filename

write(filename, '(A,i4.4,A)') trim(dirname)//'/b', blocknum, '.rst'

if ( rw == 'read' .and. .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'open_block_file',string1,source,revision,revdate)
endif

print *, 'opening file ', trim(filename), ' for ', trim(rw)

open_block_file = open_file(filename, 'unformatted', rw)
 
!print *, 'returned file descriptor is ', open_block_file

end function open_block_file



subroutine get_data(dirname, statevector)
!------------------------------------------------------------------
! open all restart files and read in the requested data item
!
 character(len=*), intent(in)  :: dirname
 real(r8),         intent(out) :: statevector(:)

integer :: ib, jb, nb, iunit, blockoffset, i

! get the dirname, construct the filenames inside open_block_file

do jb = 1, nBlocksLat
 do ib = 1, nBlocksLon
   
   nb = (jb-1) * nBlocksLon + ib

   blockoffset = nLats * ngridLon * (jb-1) + &
                 nLons * (ib-1)

!print *, 'ib,jb = ', ib, jb
!print *, 'blockoffset, nb = ', blockoffset, nb

   iunit = open_block_file(dirname, nb, 'read')

   call read_data(iunit, blockoffset, statevector)

   call close_file(iunit)
 enddo
enddo

if ( debug > 4 ) then ! A little sanity check
   write (*,*) 'variable data after read: '
   do i = 1, nfields
      write(*,*) trim(progvar(i)%varname), ' range ', &
                 minval(statevector(progvar(i)%index1:progvar(i)%indexN)), &
                 maxval(statevector(progvar(i)%index1:progvar(i)%indexN))
   enddo
endif

end subroutine get_data

subroutine put_data(dirname, dirnameout, statevector)
!------------------------------------------------------------------
! open all restart files and write out the requested data item
!
 character(len=*), intent(in) :: dirname, dirnameout
 real(r8),         intent(in) :: statevector(:)

integer :: ib, jb, nb, iunit, ounit
integer :: i, j, k, blockoffset

! get the dirname, construct the filenames inside open_block_file

if ( debug > 4 ) then ! A little sanity check
   write (*,*) 'variable data to be written: '
   do i = 1, nfields
      write(*,*) trim(progvar(i)%varname), 'range ', &
                 minval(statevector(progvar(i)%index1:progvar(i)%indexN)), &
                 maxval(statevector(progvar(i)%index1:progvar(i)%indexN))
   enddo
endif

do jb = 1, nBlocksLat
 do ib = 1, nBlocksLon
   
   nb = (jb-1) * nBlocksLon + ib

   blockoffset = nLats * ngridLon * (jb-1) + &
                 nLons * (ib-1)

!print *, 'ib,jb = ', ib, jb
!print *, 'blockoffset, nb = ', blockoffset, nb

   iunit = open_block_file(dirname,    nb, 'read' )
   ounit = open_block_file(dirnameout, nb, 'write')

   call write_data(iunit, ounit, blockoffset, statevector)

   call close_file(iunit)
   call close_file(ounit)
 enddo
enddo

end subroutine put_data


subroutine unpack_data(data3d, ivar, blockoffset, statevector)
!------------------------------------------------------------------
! put the requested data into the state vector
!
real(r8), intent(in)    :: data3d(:,:,:)
integer,  intent(in)    :: ivar         ! index into progvar struct
integer,  intent(in)    :: blockoffset
real(r8), intent(inout) :: statevector(:)

integer :: i, j, k, offset, base

!print *, 'ivar = ', ivar
base = progvar(ivar)%index1 - 1   ! FIXME: -1?
!print *, 'blockoffset, base = ', blockoffset, base

do k=1,nAlts
 do j=1,nLats
  do i=1,nLons

      offset = ((k-1) * ngridLat * ngridLon) +  &
               ((j-1) * ngridLon) +             &
               i
    if (base+blockoffset+offset < 1 .or. &
        base+blockoffset+offset > model_size) then
      print *, 'i,j,k, index: ', i, j, k, base+blockoffset+offset
    else
      statevector(base + blockoffset + offset) = data3d(i, j, k)
      !print *, 'i,j,k,varoffset = ', i,j,k,blockoffset + offset
    endif

  enddo
 enddo
enddo

end subroutine unpack_data


subroutine pack_data(statevector, ivar, blockoffset, data3d)
!------------------------------------------------------------------
! put the state vector data into a 3d array
!
real(r8), intent(in)    :: statevector(:)
integer,  intent(in)    :: ivar         ! index into progvar struct
integer,  intent(in)    :: blockoffset
real(r8), intent(inout) :: data3d(:,:,:)

integer :: i, j, k, offset, base

base = progvar(ivar)%index1 - 1  ! FIXME: -1?

do k=1,nAlts
 do j=1,nLats
  do i=1,nLons

      offset = ((k-1) * ngridLat * ngridLon) +  &
               ((j-1) * ngridLon) +             &
               i
      data3d(i,j,k) = statevector(base + blockoffset + offset) 
      !print *, 'i,j,k,varoffset = ', i,j,k,blockoffset + offset

  enddo
 enddo
enddo

end subroutine pack_data


subroutine read_data(iunit, blockoffset, statevector)
!------------------------------------------------------------------
! open all restart files and read in the requested data items
!
integer,  intent(in)    :: iunit
integer,  intent(in)    :: blockoffset
real(r8), intent(inout) :: statevector(:)

real(r8), allocatable :: temp1d(:), temp3d(:,:,:), temp4d(:,:,:,:)
integer :: i, j, count, maxsize, ivals(NSpeciesTotal)
logical :: done

! a temp array large enough to hold any of the 
! Lon,Lat or Alt array from a block plus ghost cells
allocate(temp1d(1-nGhost:max(nLons,nLats,nAlts)+nGhost))

! temp array large enough to hold 1 species, temperature, etc
allocate(temp3d(1-nGhost:nLons+nGhost, 1-nGhost:nLats+nGhost, 1-nGhost:nAlts+nGhost))

! temp array large enough to hold velocity vect, etc
maxsize = max(3, nSpecies)
allocate(temp4d(1-nGhost:nLons+nGhost, 1-nGhost:nLats+nGhost, 1-nGhost:nAlts+nGhost, maxsize))

! get past lon and lat arrays and read in alt array
print *, 'iunit in read_data = ', iunit
read(iunit) temp1d(1-nGhost:nLons+nGhost)
read(iunit) temp1d(1-nGhost:nLats+nGhost)
read(iunit) temp1d(1-nGhost:nAlts+nGhost)


call get_index_from_model_varname('NDensityS', count, ivals)
if (count > 0) then
   ! if i equals ival, use the data from the state vect
   ! otherwise read/write what's in the input file
   j = 1
   do i = 1, nSpeciesTotal
      read(iunit)  temp3d
      if (j <= count) then
         if (i == progvar(ivals(j))%model_index) then
            call unpack_data(temp3d, ivals(j), blockoffset, statevector)
            j = j + 1
         endif
      endif
   enddo
else
   ! nothing at all from this variable in the state vector.
   ! copy all data over from the input file to output file
   do i = 1, nSpeciesTotal
      read(iunit)  temp3d
   enddo
endif

call get_index_from_model_varname('IDensityS', count, ivals)
if (count > 0) then
   ! one or more items in the state vector need to replace the 
   ! data in the output file.  loop over the index list in order.
   j = 1
   do i = 1, nIons
      read(iunit)  temp3d
      if (j <= count) then
         if (i == progvar(ivals(j))%model_index) then
            ! read from input but write from state vector
            call unpack_data(temp3d, ivals(j), blockoffset, statevector)
            j = j + 1
         endif
      endif
   enddo
else
   ! nothing at all from this variable in the state vector.
   ! read past this variable
   do i = 1, nIons
      read(iunit)  temp3d
   enddo
endif

read(iunit)  temp3d
call get_index_from_model_varname('Temperature', count, ivals)
if (count > 0) then
   call unpack_data(temp3d, ivals(1), blockoffset, statevector)
endif

read(iunit) temp3d
call get_index_from_model_varname('ITemperature', count, ivals)
if (count > 0) then
   call unpack_data(temp3d, ivals(1), blockoffset, statevector)
endif

read(iunit) temp3d
call get_index_from_model_varname('eTemperature', count, ivals)
if (count > 0) then
   call unpack_data(temp3d, ivals(1), blockoffset, statevector)
endif

!print *, 'reading in temp4d for vel'
read(iunit) temp4d(:,:,:,1:3)
call get_index_from_model_varname('Velocity', count, ivals)
if (count > 0) then
   ! copy out any requested bits into state vector
   j = 1
   do i = 1, 3
      if (j <= count) then
         if (i == progvar(ivals(j))%model_index) then
            temp3d = temp4d(:,:,:,i)
            call unpack_data(temp3d, ivals(j), blockoffset, statevector)
            j = j + 1
         endif
      endif
   enddo
endif

!print *, 'reading in temp4d for ivel'
read(iunit) temp4d(:,:,:,1:3)
call get_index_from_model_varname('IVelocity', count, ivals)
if (count > 0) then
   ! copy out any requested bits into state vector
   j = 1
   do i = 1, 3
      if (j <= count) then
         if (i == progvar(ivals(j))%model_index) then
            ! read from input but write from state vector
            temp3d = temp4d(:,:,:,i)
            call unpack_data(temp3d, ivals(j), blockoffset, statevector)
            j = j + 1
         endif
      endif
   enddo
endif

!print *, 'reading in temp4d for vvel'
read(iunit) temp4d(:,:,:,1:nSpecies)
call get_index_from_model_varname('VerticalVelocity', count, ivals)
if (count > 0) then
   ! copy out any requested bits into state vector
   j = 1
   do i = 1, nSpecies
      if (j <= count) then
         if (i == progvar(ivals(j))%model_index) then
            temp3d = temp4d(:,:,:,i)
            call unpack_data(temp3d, ivals(j), blockoffset, statevector)
            j = j + 1
         endif
      endif
   enddo
endif

!print *, 'calling dealloc'
deallocate(temp1d, temp3d, temp4d)

end subroutine read_data


subroutine write_data(iunit, ounit, blockoffset, statevector)
!------------------------------------------------------------------
! open all restart files and write out the requested data item
!
integer,  intent(in) :: iunit, ounit
integer,  intent(in) :: blockoffset
real(r8), intent(in) :: statevector(:)

real(r8), allocatable :: temp1d(:), temp3d(:,:,:), temp4d(:,:,:,:), data3d(:,:,:)
integer :: i, j, count, maxsize, ivals(NSpeciesTotal)
logical :: done

! a temp array large enough to hold any of the 
! Lon,Lat or Alt array from a block plus ghost cells
allocate(temp1d(1-nGhost:max(nLons,nLats,nAlts)+nGhost))

! temp array large enough to hold 1 species, temperature, etc
allocate(temp3d(1-nGhost:nLons+nGhost, 1-nGhost:nLats+nGhost, 1-nGhost:nAlts+nGhost))
allocate(data3d(1-nGhost:nLons+nGhost, 1-nGhost:nLats+nGhost, 1-nGhost:nAlts+nGhost))

! temp array large enough to hold velocity vect, etc
maxsize = max(3, nSpecies)
allocate(temp4d(1-nGhost:nLons+nGhost, 1-nGhost:nLats+nGhost, 1-nGhost:nAlts+nGhost, maxsize))

! copy over lat, lon, alt arrays verbatim
 read(iunit) temp1d(1-nGhost:nLons+nGhost)
write(ounit) temp1d(1-nGhost:nLons+nGhost)
 read(iunit) temp1d(1-nGhost:nLats+nGhost)
write(ounit) temp1d(1-nGhost:nLats+nGhost)
 read(iunit) temp1d(1-nGhost:nAlts+nGhost)
write(ounit) temp1d(1-nGhost:nAlts+nGhost)


call get_index_from_model_varname('NDensityS', count, ivals)
if (count > 0) then
   ! if i equals ival, use the data from the state vect
   ! otherwise read/write what's in the input file
   j = 1
   do i = 1, nSpeciesTotal
      read(iunit)  temp3d
      if (j <= count) then
         if (i == progvar(ivals(j))%model_index) then
            call pack_data(statevector, ivals(j), blockoffset, data3d)
            write(ounit) data3d
            j = j + 1
         else
            write(ounit) temp3d
         endif
      else
         ! this one not in state vector, copy over from input
         write(ounit) temp3d
      endif
   enddo
else
   ! nothing at all from this variable in the state vector.
   ! copy all data over from the input file to output file
   do i = 1, nSpeciesTotal
      read(iunit)  temp3d
      write(ounit) temp3d
   enddo
endif

call get_index_from_model_varname('IDensityS', count, ivals)
if (count > 0) then
   ! one or more items in the state vector need to replace the 
   ! data in the output file.  loop over the index list in order.
   j = 1
   do i = 1, nIons
      read(iunit)  temp3d
      if (j <= count) then
         if (i == progvar(ivals(j))%model_index) then
            ! read from input but write from state vector
            call pack_data(statevector, ivals(j), blockoffset, data3d)
            write(ounit) data3d
            j = j + 1
         else
            write(ounit) temp3d
         endif
      else
         ! this one not in state vector, copy over from input
         write(ounit) temp3d
      endif
   enddo
else
   ! nothing at all from this variable in the state vector.
   ! copy all data over from the input file to output file
   do i = 1, nIons
      read(iunit)  temp3d
      write(ounit) temp3d
   enddo
endif

read(iunit)  temp3d
call get_index_from_model_varname('Temperature', count, ivals)
if (count > 0) then
   call pack_data(statevector, ivals(1), blockoffset, data3d)
   write(ounit) data3d
else
   write(ounit) temp3d
endif


read(iunit) temp3d
call get_index_from_model_varname('ITemperature', count, ivals)
if (count > 0) then
   call pack_data(statevector, ivals(1), blockoffset, data3d)
   write(ounit) data3d
else
   write(ounit) temp3d
endif

read(iunit) temp3d
call get_index_from_model_varname('eTemperature', count, ivals)
if (count > 0) then
   call pack_data(statevector, ivals(1), blockoffset, data3d)
   write(ounit) data3d
else
   write(ounit) temp3d
endif

!print *, 'reading in temp4d for vel'
read(iunit) temp4d(:,:,:,1:3)
call get_index_from_model_varname('Velocity', count, ivals)
if (count > 0) then
   ! one or more items in the state vector need to replace the 
   ! data in the output file.  loop over the index list in order.
   j = 1
   do i = 1, 3
      if (j <= count) then
         if (i == progvar(ivals(j))%model_index) then
            ! read from input but write from state vector
            call pack_data(statevector, ivals(j), blockoffset, data3d)
            temp4d(:,:,:,i) = data3d
            j = j + 1
         endif
      endif
   enddo
endif
write(ounit) temp4d(:,:,:,1:3)

!print *, 'reading in temp4d for ivel'
read(iunit) temp4d(:,:,:,1:3)
call get_index_from_model_varname('IVelocity', count, ivals)
if (count > 0) then
   ! one or more items in the state vector need to replace the 
   ! data in the output file.  loop over the index list in order.
   j = 1
   do i = 1, 3
      if (j <= count) then
         if (i == progvar(ivals(j))%model_index) then
            ! read from input but write from state vector
            call pack_data(statevector, ivals(j), blockoffset, data3d)
            temp4d(:,:,:,i) = data3d
            j = j + 1
         endif
      endif
   enddo
endif
write(ounit) temp4d(:,:,:,1:3)

!print *, 'reading in temp4d for vvel'
read(iunit) temp4d(:,:,:,1:nSpecies)
call get_index_from_model_varname('VerticalVelocity', count, ivals)
if (count > 0) then
   ! one or more items in the state vector need to replace the 
   ! data in the output file.  loop over the index list in order.
   j = 1
   do i = 1, nSpecies
      if (j <= count) then
         if (i == progvar(ivals(j))%model_index) then
            ! read from input but write from state vector
            call pack_data(statevector, ivals(j), blockoffset, data3d)
            temp4d(:,:,:,i) = data3d
            j = j + 1
         endif
      endif
   enddo
endif
write(ounit) temp4d(:,:,:,1:nSpecies)

deallocate(temp1d, temp3d, temp4d, data3d)

end subroutine write_data


function get_base_time_ncid( ncid )
!------------------------------------------------------------------
! The restart netcdf files have the start time of the experiment.
! The time array contains the time trajectory since then.
! This routine returns the start time of the experiment.

type(time_type) :: get_base_time_ncid

integer, intent(in) :: ncid

integer :: year, month, day, hour, minute, second

if ( .not. module_initialized ) call static_init_model

call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'YEAR'  , year), &
                  'get_base_time', 'get_att year')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'MONTH' , month), &
                  'get_base_time', 'get_att month')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'DAY'   , day), &
                  'get_base_time', 'get_att day')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'HOUR'  , hour), &
                  'get_base_time', 'get_att hour')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'MINUTE', minute), &
                  'get_base_time', 'get_att minute')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'SECOND', second), &
                  'get_base_time', 'get_att second')

get_base_time_ncid = set_date(year, month, day, hour, minute, second)

end function get_base_time_ncid


!------------------------------------------------------------------


function get_base_time_fname(filename)
!------------------------------------------------------------------
! The restart netcdf files have the start time of the experiment.
! The time array contains the time trajectory since then.
! This routine returns the start time of the experiment.

type(time_type) :: get_base_time_fname

character(len=*), intent(in) :: filename

integer :: ncid, year, month, day, hour, minute, second

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'get_base_time',string1,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'get_base_time', 'open '//trim(filename))
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'YEAR'  , year), &
                  'get_base_time', 'get_att year')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'MONTH' , month), &
                  'get_base_time', 'get_att month')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'DAY'   , day), &
                  'get_base_time', 'get_att day')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'HOUR'  , hour), &
                  'get_base_time', 'get_att hour')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'MINUTE', minute), &
                  'get_base_time', 'get_att minute')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'SECOND', second), &
                  'get_base_time', 'get_att second')
call nc_check(nf90_close(ncid), 'get_base_time', 'close '//trim(filename))

get_base_time_fname = set_date(year, month, day, hour, minute, second)

end function get_base_time_fname


!------------------------------------------------------------------


function get_state_time( dirname )
!------------------------------------------------------------------
! the static_init_model ensures that the model namelists are read.
!
type(time_type)              :: get_state_time
character(len=*), intent(in) :: dirname

type(time_type) :: model_offset, base_time

integer  :: iunit, i, ios
integer  :: istep
real(r8) :: tsimulation
integer  :: iyear, imonth, iday, ihour, imin, isec
integer  :: ndays,nsec

character(len=256) :: filename
character(len=100) :: cLine

if ( .not. module_initialized ) call static_init_model

tsimulation = MISSING_R8
istep       = -1
iyear       = -1
imonth      = -1
iday        = -1
ihour       = -1
imin        = -1
isec        = -1

write(filename,'(a,''/header.rst'')') trim(dirname)

iunit = open_file(trim(filename), action='read')

FILEREAD : do i = 1, 100

   read(iunit,'(a)',iostat=ios) cLine

   if (ios < 0) exit FILEREAD  ! end of file

   if (ios /= 0) then
      write(string1,*) 'cannot read ',trim(filename)
      call error_handler(E_ERR,'get_grid_info',string1,source,revision,revdate)
   endif

   select case( cLine(1:6) ) 
      case('#ISTEP')
         read(iunit,*)istep
      case('#TSIMU')
         read(iunit,*)tsimulation
      case('#TIMES')
         read(iunit,*)iyear
         read(iunit,*)imonth
         read(iunit,*)iday
         read(iunit,*)ihour
         read(iunit,*)imin
         read(iunit,*)isec
      case default
   end select

enddo FILEREAD

call close_file(iunit)

base_time      = set_date(iyear, imonth, iday, ihour, imin, isec)
ndays          = tsimulation/86400
nsec           = tsimulation - ndays*86400
model_offset   = set_time(nsec,ndays)
get_state_time = base_time + model_offset

if (debug > 8) then
   write(*,*)'get_state_time : iyear       ',iyear
   write(*,*)'get_state_time : imonth      ',imonth
   write(*,*)'get_state_time : iday        ',iday
   write(*,*)'get_state_time : ihour       ',ihour
   write(*,*)'get_state_time : imin        ',imin
   write(*,*)'get_state_time : isec        ',isec
   write(*,*)'get_state_time : tsimulation ',tsimulation
   write(*,*)'get_state_time : ndays       ',ndays
   write(*,*)'get_state_time : nsec        ',nsec

   call print_date(     base_time, 'get_state_time:model base date')
   call print_time(     base_time, 'get_state_time:model base time')
   call print_time(  model_offset, 'get_state_time:model offset')
   call print_date(get_state_time, 'get_state_time:model date')
   call print_time(get_state_time, 'get_state_time:model time')
endif

end function get_state_time


!------------------------------------------------------------------


subroutine get_index_range_string(string,index1,indexN)
!------------------------------------------------------------------
! Determine where a particular DART kind (string) exists in the 
! DART state vector.

character(len=*), intent(in)  :: string
integer,          intent(out) :: index1,indexN

integer :: i

index1 = 0
indexN = 0

write(*,*)'FIXME ... actually, just test get_index_range_string ... '

FieldLoop : do i=1,nfields
   if (progvar(i)%kind_string /= trim(string)) cycle FieldLoop
   index1 = progvar(i)%index1
   indexN = progvar(i)%indexN
   exit FieldLoop
enddo FieldLoop

if ((index1 == 0) .or. (indexN == 0)) then
   write(string1,*) 'Problem, cannot find indices for '//trim(string)
   call error_handler(E_ERR,'get_index_range_string',string1,source,revision,revdate)
endif
end subroutine get_index_range_string


subroutine get_index_range_int(dartkind,index1,indexN)
!------------------------------------------------------------------
! Determine where a particular DART kind (integer) exists in the 
! DART state vector.

integer, intent(in) :: dartkind
integer, intent(out) :: index1,indexN

integer :: i
character(len=obstypelength) :: string

index1 = 0
indexN = 0

write(*,*)'FIXME ... actually, just test get_index_range_int ... '

FieldLoop : do i=1,nfields
   if (progvar(i)%dart_kind /= dartkind) cycle FieldLoop
   index1 = progvar(i)%index1
   indexN = progvar(i)%indexN
   exit FieldLoop
enddo FieldLoop

string = get_name_for_quantity(dartkind)

if ((index1 == 0) .or. (indexN == 0)) then
   write(string1,*) 'Problem, cannot find indices for kind ',dartkind,trim(string)
   call error_handler(E_ERR,'get_index_range_int',string1,source,revision,revdate)
endif

end subroutine get_index_range_int


subroutine get_index_from_model_varname(model_varname, count, ivals)
!------------------------------------------------------------------
! Determine where any data from a given model_varname lies in the
! DART state vector.

character(len=*), intent(in) :: model_varname
integer, intent(out) :: count, ivals(:)

integer :: gindex(nfields)
integer :: i, limit

count = 0
limit = size(ivals)

FieldLoop : do i=1,nfields
   if (progvar(i)%model_varname /= model_varname) cycle FieldLoop
   count = count + 1
   if (count > limit) then
      write(string1,*) 'found too many matches, ivals needs to be larger than ', limit
      call error_handler(E_ERR,'get_index_from_model_varname',string1,source,revision,revdate)
   endif
   ! i is index into progvar array - the order of the fields in the sv
   ! model_index is index into the specific variable in the model restarts
   ivals(count) = i
   gindex(count) = progvar(i)%model_index
enddo FieldLoop

!if (count > 0) then
!   print *, 'before sort, count: ', count
!   print *, 'before sort, gindex: ', gindex(1:count)
!   print *, 'before sort, ivals: ', ivals(1:count)
!endif

! return the vals sorted by model_index order if more than 1
if (count > 1) call sortindexlist(gindex, ivals, count)

!if (count > 0) then
!   print *, 'after  sort, count: ', count
!   print *, 'after  sort, gindex: ', gindex(1:count)
!   print *, 'after  sort, ivals: ', ivals(1:count)
!endif


end subroutine get_index_from_model_varname


!------------------------------------------------------------------


function set_model_time_step()
!------------------------------------------------------------------
! the static_init_model ensures that the model namelists are read.
!
type(time_type) :: set_model_time_step

if ( .not. module_initialized ) call static_init_model

! FIXME - determine when we can stop the model

   set_model_time_step = set_time(0, 1) ! (seconds, days)

end function set_model_time_step


!------------------------------------------------------------------


subroutine get_model_restart_dirname( dirname )

character(len=*), intent(OUT) :: dirname

if ( .not. module_initialized ) call static_init_model

dirname = trim(model_restart_dirname)

end subroutine get_model_restart_dirname


!------------------------------------------------------------------


subroutine verify_state_variables( state_variables, ncid, filename, ngood, table )

character(len=*), dimension(:),   intent(in)  :: state_variables
integer,                          intent(in)  :: ncid
character(len=*),                 intent(in)  :: filename
integer,                          intent(out) :: ngood
character(len=*), dimension(:,:), intent(out) :: table

integer :: nrows, ncols, i, varid
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
      string1 = 'model_vars_nml:model state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Make sure variable exists in model restart variable list

!%!   write(string1,'(''there is no variable '',a,'' in '',a)') trim(varname), trim(filename)
!%!   call nc_check(NF90_inq_varid(ncid, trim(varname), varid), &
!%!                 'verify_state_variables', trim(string1))

   ! Make sure DART kind is valid

   if( get_index_for_quantity(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector 

   if ( debug > 0 ) then
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



!###########################################################################
!
!     ##################################################################
!     ######                                                      ######
!     ######             INTEGER FUNCTION FIND_INDEX              ######
!     ######                                                      ######
!     ##################################################################
!
!     PURPOSE:
!
!     This function returns the array index (here, the value returned by
!     find_index is designated as i) such that x is between xa(i) and xa(i+1).
!     If x is less than xa(1), then i=-1 is returned.  If x is greater than
!     xa(n), then i=-1 is returned.  It is assumed that the values of
!    xa increase monotonically with increasing i.
!
!############################################################################
!
!     Author:  David Dowell (based on "locate" algorithm in Numerical Recipes)
!
!     Creation Date:  17 November 2004
!
!############################################################################

function find_index(x, xa, n)

    integer              :: find_index
    integer,  intent(in) :: n             ! array size
    real(r8), intent(in) :: xa(n)         ! array of locations
    real(r8), intent(in) :: x             ! location of interest

    integer :: lower, upper, mid          ! lower and upper limits, and midpoint
    integer :: order

    lower = 0
    upper = n + 1

    IF( xa(1) .lt. xa(n) ) THEN
      order = 1
    ELSE
      order = -1
    ENDIF

    IF ( x .gt. maxval(xa) ) THEN
      mid = -1
    ELSEIF ( x .lt. minval(xa) ) THEN
      mid = -1
    ELSE

10    IF ((upper-lower).gt.1) THEN
        mid=(lower+upper)/2
        IF( order .eq. 1 ) THEN
          IF (x .ge. xa(mid)) THEN
            lower = mid
          ELSE
            upper = mid
          ENDIF
          go to 10
        ELSE
          IF (x .lt. xa(mid)) THEN
            lower = mid
          ELSE
            upper = mid
          ENDIF
          go to 10
        ENDIF
      ENDIF

    ENDIF

    find_index = lower

return
end function find_index




subroutine define_var_dims(myprogvar, ndims, dimids, memberdimid, unlimiteddimid, &
                       nLONdimid, nLATdimid, nALTdimid, NgridLon, NgridLat, NgridAlt)

type(progvartype),     intent(in)  :: myprogvar
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: dimids
integer,               intent(in)  :: memberdimid, unlimiteddimid
integer,               intent(in)  :: nLONdimid, nLATdimid, nALTdimid
integer,               intent(in)  :: NgridLon, NgridLat, NgridAlt

select case( myprogvar%storder ) 
case('xyz3d')

      ndims = 5

      dimids(1) = nLONdimid
      dimids(2) = nLATdimid
      dimids(3) = nALTdimid
      dimids(4) = memberdimid
      dimids(5) = unlimitedDimid

case('xy2d')

      ndims = 4

      dimids(1) = nLONdimid
      dimids(2) = nLATdimid
      dimids(3) = memberdimid
      dimids(4) = unlimitedDimid

case('x1d')

      ndims = 3

      dimids(1) = nLONdimid
      dimids(2) = memberdimid
      dimids(3) = unlimitedDimid

case('y1d')

      ndims = 3

      dimids(1) = nLATdimid
      dimids(2) = memberdimid
      dimids(3) = unlimitedDimid

case('z1d')

      ndims = 3

      dimids(1) = nALTdimid
      dimids(2) = memberdimid
      dimids(3) = unlimitedDimid

case default

      write(string1,*)'unknown storage order '//trim(myprogvar%storder)//& 
                              ' for variable '//trim(myprogvar%varname)
      call error_handler(E_ERR,'define_var_dims',string1,source,revision,revdate)

end select

return
end subroutine define_var_dims




function read_in_real(iunit,varname,filename)
integer,          intent(in) :: iunit
character(len=*), intent(in) :: varname,filename
real(r8)                     :: read_in_real

character(len=100) :: cLine
integer :: i, ios

! Read a line and remove anything after a space or TAB
read(iunit,'(a)',iostat=ios) cLine
if (ios /= 0) then 
   write(string1,*) 'cannot find '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'get_grid_dims',string1,source,revision,revdate)
endif

i=index(cLine,' ');     if( i > 0 ) cLine(i:len(cLine))=' '
i=index(cLine,char(9)); if( i > 0 ) cLine(i:len(cLine))=' '

! Now that we have a line with nothing else ... parse it
read(cLine,*,iostat=ios)read_in_real

if(ios /= 0) then
   write(string1,*)'unable to read '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'read_in_real',string1,source,revision,revdate)
endif

end function read_in_real



function read_in_int(iunit,varname,filename)
integer,          intent(in) :: iunit
character(len=*), intent(in) :: varname,filename
integer                      :: read_in_int

character(len=100) :: cLine
integer :: i, ios

! Read a line and remove anything after a space or TAB
read(iunit,'(a)',iostat=ios) cLine
if (ios /= 0) then 
   write(string1,*) 'cannot find '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'get_grid_dims',string1,source,revision,revdate)
endif

i=index(cLine,' ');     if( i > 0 ) cLine(i:len(cLine))=' '
i=index(cLine,char(9)); if( i > 0 ) cLine(i:len(cLine))=' '

read(cLine,*,iostat=ios)read_in_int

if(ios /= 0) then
   write(string1,*)'unable to read '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'read_in_int',string1,source,revision,revdate,&
             text2=cLine)
endif

end function read_in_int


subroutine sortindexlist(list, x, count)
! sort list x into order based on values in list.
! should only be called on short ( < hundreds) of values or will be slow

integer, intent(inout) :: list(:)
integer, intent(inout) :: x(:)
integer, intent(in)    :: count

integer :: tmp
integer :: j, k

!  DO A N^2 SORT - only use for short lists
do j = 1, count - 1
   do k = j + 1, count
      ! if list() is in wrong order, exchange both list items and 
      ! items in x array.
      if(list(j) .gt. list(k)) then
         tmp = list(k)
         list(k) = list(j)
         list(j) = tmp
         tmp = x(k)
         x(k) = x(j)
         x(j) = tmp
      end if
   end do
end do
end subroutine sortindexlist

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
