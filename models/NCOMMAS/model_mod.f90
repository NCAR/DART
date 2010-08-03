! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! This is the interface between the ncommas model and DART.

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, SECPERDAY, MISSING_R8, rad2deg, PI
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
                             get_close_obs_init, loc_get_close_obs => get_close_obs

use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper,                    &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, file_exist, find_textfile_dims,        &
                             file_to_text

use     obs_kind_mod, only : KIND_U_WIND_COMPONENT,        &   ! index 1
                             KIND_V_WIND_COMPONENT,        &   ! index 2
                             KIND_VERTICAL_VELOCITY,       &   ! index 3
                             KIND_POTENTIAL_TEMPERATURE,   &   ! index 4
                             KIND_RADAR_REFLECTIVITY,      &   ! index 5
                             KIND_VERTICAL_VORTICITY,      &   ! index 6
                             KIND_EXNER_FUNCTION,          &   ! index 7
                             KIND_VAPOR_MIXING_RATIO,      &   ! index 8
                             KIND_CLOUDWATER_MIXING_RATIO, &   ! index 9 
                             KIND_RAINWATER_MIXING_RATIO,  &   ! index 10
                             KIND_ICE_MIXING_RATIO,        &   ! index 11
                             KIND_SNOW_MIXING_RATIO,       &   ! index 12
                             KIND_GRAUPEL_MIXING_RATIO,    &   ! index 13
                             get_raw_obs_kind_name

use mpi_utilities_mod, only: my_task_id

use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use  dart_ncommas_mod, only: set_model_time_step, grid_type, get_grid, &
                             get_grid_dims, get_base_time, get_state_time, &
                             get_ncommas_restart_filename, write_ncommas_namelist

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
public :: get_gridsize, restart_file_to_sv, sv_to_restart_file

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = '$URL$', &
   revision = '$Revision$', &
   revdate  = '$Date$'

character(len=256) :: string1, string2
logical, save :: module_initialized = .false.

character(len=256) :: ncommas_filename

! Storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq

! things which can/should be in the model_nml
integer  :: assimilation_period_days = 1
integer  :: assimilation_period_seconds = 0
real(r8) :: model_perturbation_amplitude = 0.2
logical  :: output_state_vector = .true.
integer  :: debug = 0   ! turn up for more and more debug messages
character(len=32):: calendar

namelist /model_nml/  &
   output_state_vector,         &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   calendar,                    &
   debug

!------------------------------------------------------------------
!
! The DART state vector will consist of:  
!
! scalar PSFC long_name = "SURFACE PRESSURE"
! scalar TSFC long_name = "SURFACE TEMPERATURE AT GROUND"
! scalar QSFC long_name = "SURFACE MIXING RATIO AT GROUND"
!  U    long_name = "X-WIND COMPONENT"      float   U(TIME, ZC, YC, XE) 
!  V    long_name = "Y-WIND COMPONENT"      float   V(TIME, ZC, YE, XC)
!  W    long_name = "Z-WIND COMPONENT"      float   W(TIME, ZE, YC, XC)
!  TH   long_name = "POTENTIAL TEMPERATURE" float  TH(TIME, ZC, YC, XC)
!  DBZ  long_name = "RADAR REFLECTIVITY"    float DBZ(TIME, ZC, YC, XC)
!  WZ   long_name = "VERTICAL VORTICITY"    float  WZ(TIME, ZC, YC, XC)
!  PI   long_name = "PERT. EXNER"	    float  PI(TIME, ZC, YC, XC)
!  QV   long_name = "VAPOR MIXING RATIO"    float  QV(TIME, ZC, YC, XC)
!  QC   long_name = "CLOUD MIXING RATIO"    float  QC(TIME, ZC, YC, XC)
!  QR   long_name = "RAIN MIXING RATIO"     float  QR(TIME, ZC, YC, XC)
!  QI   long_name = "ICE MIXING RATIO"      float  QI(TIME, ZC, YC, XC)
!  QS   long_name = "SNOW MIXING RATIO"     float  QS(TIME, ZC, YC, XC)
!  QH   long_name = "GRAUPEL MIXING RATIO"  float  QH(TIME, ZC, YC, XC)
!
! FIXME: make this completely namelist driven,
!        both contents and order of vars.
!        Example: WRF input.nml sets kind_string, etc.
!------------------------------------------------------------------

! FIXME: this ought to be set by the length of the namelist.
integer, parameter :: n3dfields = 13
integer, parameter :: n2dfields = 0
integer, parameter :: nfields   = n3dfields + n2dfields

! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer :: numdims
   integer :: varsize     ! prod(dimlens(1:numdims))
   integer :: index1      ! location in dart state vector of first occurrence
   integer :: indexN      ! location in dart state vector of last  occurrence
   integer :: dart_kind
   character(len=32) :: kind_string
end type progvartype

type(progvartype), dimension(nfields) :: progvar

character(len=128) :: progvarnames(nfields) = &
                         (/ 'U   ', 'V   ', 'W   ', 'TH  ', 'DBZ ', &
                            'WZ  ', 'PI  ', 'QV  ', 'QC  ', 'QR  ', &
                            'QI  ', 'QS  ', 'QH  ' /)

integer :: progvarkinds(nfields) = (/ &
   KIND_U_WIND_COMPONENT, &
   KIND_V_WIND_COMPONENT, &
   KIND_VERTICAL_VELOCITY, &
   KIND_POTENTIAL_TEMPERATURE, &
   KIND_RADAR_REFLECTIVITY, &
   KIND_VERTICAL_VORTICITY, &
   KIND_EXNER_FUNCTION, &
   KIND_VAPOR_MIXING_RATIO, &
   KIND_CLOUDWATER_MIXING_RATIO, &
   KIND_RAINWATER_MIXING_RATIO, &
   KIND_ICE_MIXING_RATIO, &
   KIND_SNOW_MIXING_RATIO, &
   KIND_GRAUPEL_MIXING_RATIO  /)

integer :: start_index(nfields)

! Grid parameters - the values will be read from a
! standard ncommas namelist and filled in here.

! Each spatial dimension has a staggered counterpart.
integer :: nxc=-1, nyc=-1, nzc=-1    ! scalar grid positions
integer :: nxe=-1, nye=-1, nze=-1    ! staggered grid positions

! locations of cell centers (C) and edges (E) for each axis.
real(r8), allocatable :: ZC(:), ZE(:)

! These arrays store the longitude and latitude of the lower left corner
real(r8), allocatable :: ULAT(:,:), ULON(:,:)  ! XE,YC
real(r8), allocatable :: VLAT(:,:), VLON(:,:)  ! XC,YE
real(r8), allocatable :: WLAT(:,:), WLON(:,:)  ! XC,YC

integer               :: model_size      ! the state vector length
type(time_type)       :: model_time      ! valid time of the model state
type(time_type)       :: model_timestep  ! smallest time to adv model
real(r8), allocatable :: ens_mean(:)     ! may be needed for forward ops

INTERFACE vector_to_prog_var
      MODULE PROCEDURE vector_to_2d_prog_var
      MODULE PROCEDURE vector_to_3d_prog_var
END INTERFACE

!------------------------------------------------
! These bits are left over from the POP dipole grid.
!------------------------------------------------

! The regular grid used for dipole interpolation divides the sphere into
! a set of regularly spaced lon-lat boxes. The number of boxes in
! longitude and latitude are set by num_reg_x and num_reg_y. Making the
! number of regular boxes smaller decreases the computation required for
! doing each interpolation but increases the static storage requirements
! and the initialization computation (which seems to be pretty small).
integer, parameter :: num_reg_x = 90, num_reg_y = 90

! The max_reg_list_num controls the size of temporary storage used for
! initializing the regular grid. Four arrays
! of size num_reg_x*num_reg_y*max_reg_list_num are needed. The initialization
! fails and returns an error if max_reg_list_num is too small. A value of
! 30 is sufficient for the x3 ncommas grid with 180 regular lon and lat boxes 
! and a value of 80 is sufficient for for the x1 grid.
integer, parameter :: max_reg_list_num = 80

! The dipole interpolation keeps a list of how many and which dipole quads
! overlap each regular lon-lat box. The number for the u and t grids are 
! stored in u_dipole_num and t_dipole_num. The allocatable arrays
! u_dipole_lon(lat)_list and t_dipole_lon(lat)_list list the longitude 
! and latitude indices for the overlapping dipole quads. The entry in
! u_dipole_start and t_dipole_start for a given regular lon-lat box indicates
! where the list of dipole quads begins in the u_dipole_lon(lat)_list and
! t_dipole_lon(lat)_list arrays.

integer :: u_dipole_start(num_reg_x, num_reg_y)
integer :: u_dipole_num  (num_reg_x, num_reg_y) = 0
integer :: t_dipole_start(num_reg_x, num_reg_y)
integer :: t_dipole_num  (num_reg_x, num_reg_y) = 0
integer, allocatable :: u_dipole_lon_list(:), t_dipole_lon_list(:)
integer, allocatable :: u_dipole_lat_list(:), t_dipole_lat_list(:)

! Need to check for pole quads: for now we are not interpolating in them
integer :: pole_x, t_pole_y, u_pole_y


! Have a global variable saying whether this is dipole or regular lon-lat grid
! This should be initialized static_init_model. Code to do this is below.
logical :: dipole_grid

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

end subroutine adv_1step



subroutine get_state_meta_data(index_in, location, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument kind
! can be returned if the model has more than one type of field (for
! instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

real(r8):: lat, lon, height
integer :: lon_index, lat_index, height_index, local_var

call get_state_indices(index_in, lat_index, lon_index, height_index, local_var)

if     (is_on_ugrid(local_var)) then
   lon = ULON(lon_index, lat_index)
   lat = ULAT(lon_index, lat_index)
elseif (is_on_vgrid(local_var)) then
   lon = VLON(lon_index, lat_index)
   lat = VLAT(lon_index, lat_index)
else 
   lon = WLON(lon_index, lat_index)
   lat = WLAT(lon_index, lat_index)
endif

if (debug > 5) print *, 'lon, lat, height = ', lon, lat, height

location = set_location(lon, lat, height, VERTISHEIGHT)

if (present(var_type)) then
   var_type = local_var
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
! Harvest a ton of information from the NCOMMAS restart file
! about grid sizes, grid contents, variable sizes, etc..

! Local variables - all the important ones have module scope
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname
integer :: ncid, VarID, numdims, dimlen, varsize
integer :: iunit, io, ivar, i, index1, indexN
integer :: ss, dd

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

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

!---------------------------------------------------------------
! Set the time step ... causes ncommas namelists to be read.
! Ensures model_timestep is multiple of 'dynamics_timestep'

call set_calendar_type( calendar )   ! comes from model_mod_nml

model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)

call get_ncommas_restart_filename( ncommas_filename )

!---------------------------------------------------------------
! 1) get grid dimensions
! 2) allocate space for the grids 
! 3) read them, convert them from X-Y-Z to lat-lon-z

call get_grid_dims(nxc, nxe, nyc, nye, nzc, nze )

allocate(ULAT(nxe,nyc), ULON(nxe,nyc))
allocate(VLAT(nxc,nye), VLON(nxc,nye))
allocate(WLAT(nxc,nyc), WLON(nxc,nyc))
allocate(  ZC(  nzc  ),   ZE(  nze  ))

call get_grid(nxc, nxe, nyc, nye, nzc, nze, &
              ULAT, ULON, VLAT, VLON, WLAT, WLON, ZC, ZE)

!---------------------------------------------------------------
! compute the offsets into the state vector for the start of each
! different variable type. Requires reading shapes from the NCOMMAS
! restart file.
!
! FIXME - this should go in dart_ncommas_mod.f90 
! as well as the progvartype declaration, should be query routines.
!
! Record where in the state vector the data type changes
! from one type to another, by computing the starting
! index for each block of data.

call nc_check( nf90_open(trim(ncommas_filename), NF90_NOWRITE, ncid), &
                  'static_init_model', 'open '//trim(ncommas_filename))

! Find the Time (Unlimited) dimension - so we can skip it.
call nc_check(nf90_Inquire(ncid,nDimensions,nVariables,nAttributes,unlimitedDimID),&
                    'static_init_model', 'inquire '//trim(ncommas_filename))

index1  = 1;
indexN  = 0;
do ivar = 1, nfields 

   varname = adjustl(progvarnames(ivar))
   string2 = trim(ncommas_filename)//' '//trim(varname)

   progvar(ivar)%varname = trim(varname)
   progvar(ivar)%dimlens = 0

   call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), &
            'static_init_model', 'inq_varid '//trim(string2))

   call nc_check( nf90_get_att(ncid, VarId, 'long_name' , progvar(ivar)%long_name), &
            'static_init_model', 'get_att long_name '//trim(string2))

   call nc_check( nf90_get_att(ncid, VarId, 'units' , progvar(ivar)%units), &
            'static_init_model', 'get_att units '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid, VarId, dimids=dimIDs, ndims=numdims), &
            'static_init_model', 'inquire '//trim(string2))

   progvar(ivar)%numdims = numdims

   varsize = 1
   DimensionLoop : do i = 1,numdims

      if (dimIDs(i) == unlimitedDimID) then
         dimlen = 1
      else   
         write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
         call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
                                          'static_init_model', string1)
      endif
      progvar(ivar)%dimlens(i) = dimlen
      varsize = varsize * dimlen

   enddo DimensionLoop

   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1 
   index1                    = index1 + varsize      ! sets up for next variable
   progvar(ivar)%dart_kind   = progvarkinds(ivar)
   progvar(ivar)%kind_string = get_raw_obs_kind_name(progvar(ivar)%dart_kind)

   if (do_output()) then
      write(logfileunit,*)
      write(logfileunit,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(logfileunit,*) '  long_name ',trim(progvar(ivar)%long_name)
      write(logfileunit,*) '  units     ',trim(progvar(ivar)%units)
      write(logfileunit,*) '  numdims   ',progvar(ivar)%numdims
      write(logfileunit,*) '  dimlens   ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      write(logfileunit,*) '  varsize   ',progvar(ivar)%varsize
      write(logfileunit,*) '  index1    ',progvar(ivar)%index1
      write(logfileunit,*) '  indexN    ',progvar(ivar)%indexN
      write(logfileunit,*) '  dart_kind ',progvar(ivar)%dart_kind
      write(logfileunit,*) '  kind_string ',progvar(ivar)%kind_string

      write(     *     ,*)
      write(     *     ,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(     *     ,*) '  long_name ',trim(progvar(ivar)%long_name)
      write(     *     ,*) '  units     ',trim(progvar(ivar)%units)
      write(     *     ,*) '  numdims   ',progvar(ivar)%numdims
      write(     *     ,*) '  dimlens   ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      write(     *     ,*) '  varsize   ',progvar(ivar)%varsize
      write(     *     ,*) '  index1    ',progvar(ivar)%index1
      write(     *     ,*) '  indexN    ',progvar(ivar)%indexN
      write(     *     ,*) '  dart_kind ',progvar(ivar)%dart_kind
      write(     *     ,*) '  kind_string ',progvar(ivar)%kind_string
   endif

enddo

model_size = progvar(nfields)%indexN

if (do_output()) then
  write(logfileunit, *)'grid: nx[ce], ny[ce], nz[ce] = ', nxc, nxe, nyc, nye, nzc, nze
  write(     *     , *)'grid: nx[ce], ny[ce], nz[ce] = ', nxc, nxe, nyc, nye, nzc, nze
  write(logfileunit, *)'model_size = ', model_size
  write(     *     , *)'model_size = ', model_size
endif

allocate( ens_mean(model_size) )

! Initialize the interpolation routines
call init_interp()

end subroutine static_init_model



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

! if ( .not. module_initialized ) call static_init_model

deallocate(ULAT, ULON, VLAT, VLON, WLAT, WLON)
deallocate(ZC, ZE)

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
integer :: NxcDimID, NycDimID, NzcDimID
integer :: NxeDimID, NyeDimID, NzeDimID

integer :: ulonVarID, ulatVarID
integer :: vlonVarID, vlatVarID
integer :: wlonVarID, wlatVarID
integer :: ZEVarID, ZCVarID

! for the prognostic variables
integer :: SVarID, TVarID, UVarID, VVarID, PSURFVarID 

!----------------------------------------------------------------------
! variables for the namelist output
!----------------------------------------------------------------------

character(len=129), allocatable, dimension(:) :: textblock
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen
logical :: has_ncommas_namelist

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

integer :: i
character(len=128)  :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename, '(a, i3)') 'ncFileID', ncFileID

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
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'ncommas' ), &
           'nc_write_model_atts', 'model put '//trim(filename))

!-------------------------------------------------------------------------------
! Determine shape of most important namelist
!-------------------------------------------------------------------------------

call find_textfile_dims('ncommas_in', nlines, linelen)
if (nlines > 0) then
  has_ncommas_namelist = .true.
else
  has_ncommas_namelist = .false.
endif

if (debug > 0)    print *, 'ncommas namelist: nlines, linelen = ', nlines, linelen
  
if (has_ncommas_namelist) then 
   allocate(textblock(nlines))
   textblock = ''

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nlines', &
                 len = nlines, dimid = nlinesDimID), &
                 'nc_write_model_atts', 'def_dim nlines ')

   call nc_check(nf90_def_var(ncFileID,name='ncommas_in', xtype=nf90_char,    &
                 dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID), &
                 'nc_write_model_atts', 'def_var ncommas_in')
   call nc_check(nf90_put_att(ncFileID, nmlVarID, 'long_name',       &
                 'contents of ncommas_in namelist'), 'nc_write_model_atts', 'put_att ncommas_in')

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
   
   call nc_check(nf90_def_dim(ncid=ncFileID, name='XC', &
          len = nxc, dimid = NxcDimID),'nc_write_model_atts', 'xc def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='XE', &
          len = nxe, dimid = NxeDimID),'nc_write_model_atts', 'xe def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='YC', &
          len = nyc, dimid = NycDimID),'nc_write_model_atts', 'yc def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='YE', &
          len = nye, dimid = NyeDimID),'nc_write_model_atts', 'ye def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='ZC', &
          len = nzc, dimid = NzcDimID),'nc_write_model_atts', 'zc def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='ZE', &
          len = nze, dimid = NzeDimID),'nc_write_model_atts', 'ze def_dim '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Coordinate Variables and the Attributes
   !----------------------------------------------------------------------------


   ! U Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='ULON', xtype=nf90_real, &
                 dimids=(/ NxeDimID, NycDimID /), varid=ulonVarID),&
                 'nc_write_model_atts', 'ULON def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'long_name', 'longitudes of U grid'), &
                 'nc_write_model_atts', 'ULON long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'cartesian_axis', 'X'),  &
                 'nc_write_model_atts', 'ULON cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'ULON units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'ULON valid_range '//trim(filename))

   ! U Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='ULAT', xtype=nf90_real, &
                 dimids=(/ NxeDimID, NycDimID /), varid=ulatVarID),&
                 'nc_write_model_atts', 'ULAT def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'long_name', 'latitudes of U grid'), &
                 'nc_write_model_atts', 'ULAT long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'cartesian_axis', 'Y'),   &
                 'nc_write_model_atts', 'ULAT cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'units', 'degrees_north'),  &
                 'nc_write_model_atts', 'ULAT units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'ULAT valid_range '//trim(filename))

   ! V Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='VLON', xtype=nf90_real, &
                 dimids=(/ NxcDimID, NyeDimID /), varid=vlonVarID),&
                 'nc_write_model_atts', 'vlon def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, vlonVarID, 'long_name', 'longitudes of V grid'), &
                 'nc_write_model_atts', 'vlon long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, vlonVarID, 'cartesian_axis', 'X'),   &
                 'nc_write_model_atts', 'vlon cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, vlonVarID, 'units', 'degrees_east'),  &
                 'nc_write_model_atts', 'vlon units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, vlonVarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'vlon valid_range '//trim(filename))

   ! V Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='VLAT', xtype=nf90_real, &
                 dimids= (/ NxcDimID, NyeDimID /), varid=vlatVarID), &
                 'nc_write_model_atts', 'vlat def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, vlatVarID, 'long_name', 'latitudes of V grid'), &
                 'nc_write_model_atts', 'vlat long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, vlatVarID, 'cartesian_axis', 'Y'),   &
                 'nc_write_model_atts', 'vlat cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, vlatVarID, 'units', 'degrees_north'),  &
                 'nc_write_model_atts', 'vlat units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, vlatVarID, 'valid_range', (/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'vlat valid_range '//trim(filename))

   ! W Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='WLON', xtype=nf90_real, &
                 dimids=(/ NxcDimID, NycDimID /), varid=wlonVarID),&
                 'nc_write_model_atts', 'wlon def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, wlonVarID, 'long_name', 'longitudes of all others... grid'), &
                 'nc_write_model_atts', 'wlon long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, wlonVarID, 'cartesian_axis', 'X'),   &
                 'nc_write_model_atts', 'wlon cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, wlonVarID, 'units', 'degrees_east'),  &
                 'nc_write_model_atts', 'wlon units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, wlonVarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'wlon valid_range '//trim(filename))

   ! V Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='WLAT', xtype=nf90_real, &
                 dimids= (/ NxcDimID, NycDimID /), varid=wlatVarID), &
                 'nc_write_model_atts', 'wlat def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, wlatVarID, 'long_name', 'latitudes of all others ... grid'), &
                 'nc_write_model_atts', 'wlat long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, wlatVarID, 'cartesian_axis', 'Y'),   &
                 'nc_write_model_atts', 'wlat cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, wlatVarID, 'units', 'degrees_north'),  &
                 'nc_write_model_atts', 'wlat units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, wlatVarID, 'valid_range', (/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'wlat valid_range '//trim(filename))

   ! heights
   call nc_check(nf90_def_var(ncFileID,name='ZE', xtype=nf90_real, &
                 dimids=NzeDimID, varid= ZEVarID), &
                 'nc_write_model_atts', 'ZE def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZEVarID, 'long_name', 'height at grid edges'), &
                 'nc_write_model_atts', 'ZE long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZEVarID, 'cartesian_axis', 'Z'),   &
                 'nc_write_model_atts', 'ZE cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZEVarID, 'units', 'meters'),  &
                 'nc_write_model_atts', 'ZE units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZEVarID, 'positive', 'down'),  &
                 'nc_write_model_atts', 'ZE units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZEVarID, 'comment', &
                 'more positive is closer to the center of the earth'),  &
                 'nc_write_model_atts', 'ZE comment '//trim(filename))

   ! heights
   call nc_check(nf90_def_var(ncFileID,name='ZC',xtype=nf90_real, &
                 dimids=NzcDimID,varid=ZCVarID), &
                 'nc_write_model_atts', 'ZC def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, 'long_name', 'height at grid centroids'), &
                 'nc_write_model_atts', 'ZC long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, 'cartesian_axis', 'Z'),   &
                 'nc_write_model_atts', 'ZC cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, 'units', 'meters'),  &
                 'nc_write_model_atts', 'ZC units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, 'positive', 'down'),  &
                 'nc_write_model_atts', 'ZC units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, 'comment', &
                 'more positive is closer to the center of the earth'),  &
                 'nc_write_model_atts', 'ZC comment '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------

!    call nc_check(nf90_def_var(ncid=ncFileID, name='SALT', xtype=nf90_real, &
!          dimids = (/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=SVarID),&
!          'nc_write_model_atts', 'S def_var '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, SVarID, 'long_name', 'salinity'), &
!          'nc_write_model_atts', 'S long_name '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, SVarID, 'units', 'kg/kg'), &
!          'nc_write_model_atts', 'S units '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, SVarID, 'missing_value', NF90_FILL_REAL), &
!          'nc_write_model_atts', 'S missing '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, SVarID, '_FillValue', NF90_FILL_REAL), &
!          'nc_write_model_atts', 'S fill '//trim(filename))
! 
! 
!    call nc_check(nf90_def_var(ncid=ncFileID, name='TEMP', xtype=nf90_real, &
!          dimids=(/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=TVarID),&
!          'nc_write_model_atts', 'T def_var '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, TVarID, 'long_name', 'Potential Temperature'), &
!          'nc_write_model_atts', 'T long_name '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, TVarID, 'units', 'deg C'), &
!          'nc_write_model_atts', 'T units '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, TVarID, 'units_long_name', 'degrees celsius'), &
!          'nc_write_model_atts', 'T units_long_name '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, TVarID, 'missing_value', NF90_FILL_REAL), &
!          'nc_write_model_atts', 'T missing '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, TVarID, '_FillValue', NF90_FILL_REAL), &
!          'nc_write_model_atts', 'T fill '//trim(filename))
! 
! 
!    call nc_check(nf90_def_var(ncid=ncFileID, name='UVEL', xtype=nf90_real, &
!          dimids=(/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=UVarID),&
!          'nc_write_model_atts', 'U def_var '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, UVarID, 'long_name', 'U velocity'), &
!          'nc_write_model_atts', 'U long_name '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, UVarID, 'units', 'cm/s'), &
!          'nc_write_model_atts', 'U units '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, UVarID, 'units_long_name', 'centimeters per second'), &
!          'nc_write_model_atts', 'U units_long_name '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, UVarID, 'missing_value', NF90_FILL_REAL), &
!          'nc_write_model_atts', 'U missing '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, UVarID, '_FillValue', NF90_FILL_REAL), &
!          'nc_write_model_atts', 'U fill '//trim(filename))
! 
! 
!    call nc_check(nf90_def_var(ncid=ncFileID, name='VVEL', xtype=nf90_real, &
!          dimids=(/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=VVarID),&
!          'nc_write_model_atts', 'V def_var '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, VVarID, 'long_name', 'V Velocity'), &
!          'nc_write_model_atts', 'V long_name '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, VVarID, 'units', 'cm/s'), &
!          'nc_write_model_atts', 'V units '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, VVarID, 'units_long_name', 'centimeters per second'), &
!          'nc_write_model_atts', 'V units_long_name '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, VVarID, 'missing_value', NF90_FILL_REAL), &
!          'nc_write_model_atts', 'V missing '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, VVarID, '_FillValue', NF90_FILL_REAL), &
!          'nc_write_model_atts', 'V fill '//trim(filename))
! 
! 
!    call nc_check(nf90_def_var(ncid=ncFileID, name='PSURF', xtype=nf90_real, &
!          dimids=(/NlonDimID,NlatDimID,MemberDimID,unlimitedDimID/),varid=PSURFVarID), &
!          'nc_write_model_atts', 'PSURF def_var '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, PSURFVarID, 'long_name', 'surface pressure'), &
!          'nc_write_model_atts', 'PSURF long_name '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, PSURFVarID, 'units', 'dyne/cm2'), &
!          'nc_write_model_atts', 'PSURF units '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, PSURFVarID, 'missing_value', NF90_FILL_REAL), &
!          'nc_write_model_atts', 'PSURF missing '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, PSURFVarID, '_FillValue', NF90_FILL_REAL), &
!          'nc_write_model_atts', 'PSURF fill '//trim(filename))

   ! Finished with dimension/variable definitions, must end 'define' mode to fill.

   call nc_check(nf90_enddef(ncfileID), 'prognostic enddef '//trim(filename))

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables
   !----------------------------------------------------------------------------

   call nc_check(nf90_put_var(ncFileID, ulonVarID, ULON ), &
                'nc_write_model_atts', 'ULON put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, ulatVarID, ULAT ), &
                'nc_write_model_atts', 'ULAT put_var '//trim(filename))

   call nc_check(nf90_put_var(ncFileID, vlonVarID, VLON ), &
                'nc_write_model_atts', 'VLON put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, vlatVarID, VLAT ), &
                'nc_write_model_atts', 'VLAT put_var '//trim(filename))

   call nc_check(nf90_put_var(ncFileID, wlonVarID, WLON ), &
                'nc_write_model_atts', 'WLON put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, wlatVarID, WLAT ), &
                'nc_write_model_atts', 'WLAT put_var '//trim(filename))

   call nc_check(nf90_put_var(ncFileID, ZCVarID, ZC ), &
                'nc_write_model_atts', 'ZC put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, ZEVarID, ZE ), &
                'nc_write_model_atts', 'ZE put_var '//trim(filename))

endif

!-------------------------------------------------------------------------------
! Fill the variables we can
!-------------------------------------------------------------------------------

if (has_ncommas_namelist) then
   call file_to_text('ncommas_in', textblock)
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



function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
! TJH 24 Oct 2006 -- Writes the model variables to a netCDF file.
!
! TJH 29 Jul 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_96 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
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
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname 
integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: i, ivar, VarID, numdims, dimlen

real(r8), allocatable, dimension(:)       :: data_1d
real(r8), allocatable, dimension(:,:)     :: data_2d
real(r8), allocatable, dimension(:,:,:)   :: data_3d
real(r8), allocatable, dimension(:,:,:,:) :: data_4d

character(len=128) :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename, '(a, i3)') 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
!-------------------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
              'nc_write_model_vars', 'inquire '//trim(filename))

if ( output_state_vector ) then

   call nc_check(NF90_inq_varid(ncFileID, 'state', VarID), &
                 'nc_write_model_vars', 'state inq_varid '//trim(filename))
   call nc_check(NF90_put_var(ncFileID,VarID,statevec,start=(/1,copyindex,timeindex/)),&
                 'nc_write_model_vars', 'state put_var '//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   ! Replace missing values (0.0) with netcdf missing value.
   !----------------------------------------------------------------------------

   do ivar = 1,nfields

      varname = trim(progvar(ivar)%varname)
      string2 = trim(filename)//' '//trim(varname)

      ! ensure netCDF variable is conformable 
      ! the TIME (unlimited) dimension will be skipped

      call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'nc_write_model_vars', 'inq_varid '//trim(string2))

      call nc_check(nf90_inquire_variable(ncFileID,VarId,dimids=dimIDs,ndims=numdims), &
            'nc_write_model_vars', 'inquire '//trim(string2))

      ConformableDimensions : do i = 1,numdims
         if ( dimIDs(i) == unlimitedDimID ) cycle ConformableDimensions

         write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
         call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
               'nc_write_model_vars', string1)

         if ( dimlen /= progvar(ivar)%dimlens(i) ) then
            write(string1,*) trim(string2),'dim/dimlen',i,dimlen,'not',progvar(ivar)%dimlens(i)
            call error_handler(E_ERR,'nc_write_model_vars',string1,source,revision,revdate)
         endif
      enddo ConformableDimensions

!     if ( progvar(ivar)%numdims == 1) then
!     elseif ( progvar(ivar)%numdims == 1) then
!     else
!     endif

!     call vector_to_prog_var(statevec, S_index, data_3d)
!     where (data_3d == 0.0_r8) data_3d = NF90_FILL_REAL
!     call nc_check(NF90_inq_varid(ncFileID, 'SALT', VarID), &
!                  'nc_write_model_vars', 'S inq_varid '//trim(filename))
!     call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
!                  'nc_write_model_vars', 'S put_var '//trim(filename))

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

integer :: i, var_type
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
!  elseif (base_array(3) == missing_r8) then
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
            (local_obs_array(3) == missing_r8)) .or.  &
            (istatus2 == 1)                   ) then
            dist(k) = 1.0e9
      else
            dist(k) = get_dist(base_obs_loc, local_obs_loc, base_obs_kind, obs_kind(t_ind))
      endif

   enddo
endif

end subroutine get_close_obs



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



subroutine get_gridsize(num_x, num_y, num_z)
 integer, intent(out) :: num_x, num_y, num_z
!------------------------------------------------------------------
! public utility routine.

if ( .not. module_initialized ) call static_init_model

 num_x = nxc
 num_y = nyc
 num_z = nzc

end subroutine get_gridsize



subroutine restart_file_to_sv(filename, state_vector, model_time)
!------------------------------------------------------------------
! Reads the current time and state variables from a ncommas restart
! file and packs them into a dart state vector.

character(len=*), intent(in)    :: filename 
real(r8),         intent(inout) :: state_vector(:)
type(time_type),  intent(out)   :: model_time

! temp space to hold data while we are reading it
integer  :: mystart(1), mycount(1)
integer  :: i, j, k, l, ni, nj, nk, nl, ivar, indx
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array
real(r8), allocatable, dimension(:,:,:)     :: data_3d_array
real(r8), allocatable, dimension(:,:,:,:)   :: data_4d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: varname 
integer :: VarID, numdims, dimlen
integer :: ncid, year, month, day, hour, minute, second, nc_rc
character(len=256) :: myerrorstring 

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

! Check that the input file exists ... 

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'restart_file_to_sv',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncid), &
             'restart_file_to_sv','open '//trim(filename))

model_time = get_state_time(ncid, filename)

if (do_output()) &
    call print_time(model_time,'time for restart file '//trim(filename))
if (do_output()) &
    call print_date(model_time,'date for restart file '//trim(filename))

! Start counting and filling the state vector one item at a time,
! repacking the Nd arrays into a single 1d list of numbers.

indx = 1

! If these arrays have an extra dimension (like TIME), and it is only 1,
! we might have to make the query code smarter, and might have to set
! a start and count on the get_var() calls.  If it is more than 1, then
! we have a problem.
do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   myerrorstring = trim(filename)//' '//trim(varname)

   ! determine the shape of the netCDF variable 

   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            'restart_file_to_sv', 'inq_varid '//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid,VarId,dimids=dimIDs,ndims=numdims), &
            'restart_file_to_sv', 'inquire '//trim(myerrorstring))

   do i = 1,numdims
      write(string1,'(''inquire dimension'',i2,A)') i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'restart_file_to_sv', string1)

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(myerrorstring),'dim/dimlen',i,dimlen,'not',progvar(ivar)%dimlens(i)
         call error_handler(E_ERR,'restart_file_to_sv',string1,source,revision,revdate)
      endif
   enddo

   if (numdims == 1) then
      ni = progvar(ivar)%dimlens(1)
      allocate(data_1d_array(ni))
      call nc_check(nf90_get_var(ncid, VarID, data_1d_array), &
            'restart_file_to_sv', 'get_var '//trim(varname))
      do i = 1, ni   ! size(data_1d_array,1)
         state_vector(indx) = data_1d_array(i)
         indx = indx + 1
      enddo
      deallocate(data_1d_array)
   elseif (numdims == 2) then
      ni = progvar(ivar)%dimlens(1)
      nj = progvar(ivar)%dimlens(2)
      allocate(data_2d_array(ni, nj))
      call nc_check(nf90_get_var(ncid, VarID, data_2d_array), &
            'restart_file_to_sv', 'get_var '//trim(varname))
      do j = 1, nj   ! size(data_2d_array,2)
      do i = 1, ni   ! size(data_2d_array,1)
         state_vector(indx) = data_2d_array(i, j)
         indx = indx + 1
      enddo
      enddo
      deallocate(data_2d_array)
   elseif (numdims == 3) then
      ni = progvar(ivar)%dimlens(1)
      nj = progvar(ivar)%dimlens(2)
      nk = progvar(ivar)%dimlens(3)
      allocate(data_3d_array(ni, nj, nk))
      call nc_check(nf90_get_var(ncid, VarID, data_3d_array), &
            'restart_file_to_sv', 'get_var '//trim(varname))
   
      do k = 1, nk   ! size(data_3d_array,3)
      do j = 1, nj   ! size(data_3d_array,2)
      do i = 1, ni   ! size(data_3d_array,1)
         state_vector(indx) = data_3d_array(i, j, k)
         indx = indx + 1
      enddo
      enddo
      enddo

      deallocate(data_3d_array)
   elseif (numdims == 4) then
      ni = progvar(ivar)%dimlens(1)
      nj = progvar(ivar)%dimlens(2)
      nk = progvar(ivar)%dimlens(3)
      nl = progvar(ivar)%dimlens(4)
      allocate(data_4d_array(ni, nj, nk, nl))
      call nc_check(nf90_get_var(ncid, VarID, data_4d_array), &
            'restart_file_to_sv', 'get_var '//trim(varname))
      do l = 1, nl   ! size(data_4d_array,3)
      do k = 1, nk   ! size(data_4d_array,3)
      do j = 1, nj   ! size(data_4d_array,2)
      do i = 1, ni   ! size(data_4d_array,1)
         state_vector(indx) = data_4d_array(i, j, k, l)
         indx = indx + 1
      enddo
      enddo
      enddo
      enddo
      deallocate(data_4d_array)
   else
      write(string1, *) 'no support for data array of dimension ', numdims
      call error_handler(E_ERR,'restart_file_to_sv', string1, &
                        source,revision,revdate)
   endif

enddo

call nc_check(nf90_close(ncid), &
             'restart_file_to_sv','close '//trim(filename))

end subroutine restart_file_to_sv



subroutine sv_to_restart_file(state_vector, filename, statedate)
!------------------------------------------------------------------
! Writes the current time and state variables from a dart state
! vector (1d fortran array) into a ncommas netcdf restart file.
!
real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: filename 
type(time_type),  intent(in) :: statedate

integer :: year, month, day, hour, minute, second, nowseconds
type(time_type) :: ncommas_time, ncommas_time0

! temp space to hold data while we are writing it
real(r8) :: data_2d_array(nxc,nyc), data_3d_array(nxc,nyc,nzc)

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname 
character(len=256)                    :: myerrorstring 

integer :: i, ivar, ncid, VarID, numdims, dimlen

!----------------------------------------------------------------------
! Get the show underway
!----------------------------------------------------------------------

if ( .not. module_initialized ) call static_init_model

! Check that the input file exists. 
! make sure the time tag in the restart file matches 
! the current time of the DART state ...

if ( .not. file_exist(filename)) then
   write(string1,*)trim(filename),' does not exist. FATAL error.'
   call error_handler(E_ERR,'sv_to_restart_file',string1,source,revision,revdate) 
endif

call nc_check( nf90_open(trim(filename), NF90_WRITE, ncid), &
                  'sv_to_restart_file', 'open '//trim(filename))
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'YEAR'  , year), &
                  'sv_to_restart_file', 'get_att YEAR')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'MONTH' , month), &
                  'sv_to_restart_file', 'get_att MONTH')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'DAY'   , day), &
                  'sv_to_restart_file', 'get_att DAY')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'HOUR'  , hour), &
                  'sv_to_restart_file', 'get_att HOUR')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'MINUTE', minute), &
                  'sv_to_restart_file', 'get_att MINUTE')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'SECOND', second), &
                  'sv_to_restart_file', 'get_att SECOND')

! initial time
ncommas_time0 = set_date(year, month, day, hour, minute, second)

! have to open TIME variable (not attribute) to get number of seconds
! since time 0 for current time.  put that into nowseconds
nowseconds = 300   ! FIXME - get this from netcdf
ncommas_time = ncommas_time0 + set_time(nowseconds)

if ( ncommas_time /= statedate ) then
   call print_time(statedate,'DART current time',logfileunit) 
   call print_time( ncommas_time,'ncommas  current time',logfileunit) 
   call print_time(statedate,'DART current time') 
   call print_time( ncommas_time,'ncommas  current time') 
   write(string1,*)trim(filename),' current time /= model time. FATAL error.'
   call error_handler(E_ERR,'sv_to_restart_file',string1,source,revision,revdate) 
endif

if (do_output()) &
    call print_time(ncommas_time,'time of restart file '//trim(filename))
if (do_output()) &
    call print_date(ncommas_time,'date of restart file '//trim(filename))

! FIXME: this needs to change.   read the namelist to see what variables
! are in the state vector and write them out in a loop.
! fill S, T, U, V in that order
do ivar=1, n3dfields

   varname = trim(progvarnames(ivar))//'_CUR'
   myerrorstring = trim(filename)//' '//trim(varname)

   ! Is the netCDF variable the right shape?
   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            'sv_to_restart_file', 'inq_varid '//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid,VarId,dimids=dimIDs,ndims=numdims), &
            'sv_to_restart_file', 'inquire '//trim(myerrorstring))

   if (numdims /= 3) then
      write(string1,*) trim(myerrorstring),' does not have exactly 3 dimensions'
      call error_handler(E_ERR,'sv_to_restart_file',string1,source,revision,revdate)
   endif

   do i = 1,numdims
      write(string1,'(''inquire dimension'',i2,A)') i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'sv_to_restart_file', string1)

      if (dimlen /= size(data_3d_array,i)) then
         write(string1,*) trim(myerrorstring),'dim/dimlen',i,dimlen,'not',size(data_3d_array,i)
         call error_handler(E_ERR,'sv_to_restart_file',string1,source,revision,revdate)
      endif
   enddo

   call vector_to_prog_var(state_vector, ivar, data_3d_array)

   ! Actually stuff it into the netcdf file
   call nc_check(nf90_put_var(ncid, VarID, data_3d_array), &
            'sv_to_restart_file', 'put_var '//trim(myerrorstring))

enddo

! and finally, PSURF (and any other 2d fields)
do ivar=(n3dfields+1), (n3dfields+n2dfields)

   varname = trim(progvarnames(ivar))//'_CUR'
   myerrorstring = trim(varname)//' '//trim(filename)

   ! Is the netCDF variable the right shape?

   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            'sv_to_restart_file', 'inq_varid '//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid,VarId,dimids=dimIDs,ndims=numdims), &
            'sv_to_restart_file', 'inquire '//trim(myerrorstring))

   if (numdims /= 2) then
      write(string1,*) trim(myerrorstring),' does not have exactly 2 dimensions'
      call error_handler(E_ERR,'sv_to_restart_file',string1,source,revision,revdate)
   endif

   do i = 1,numdims
      write(string1,'(''inquire dimension'',i2,A)') i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'sv_to_restart_file', string1)

      if (dimlen /= size(data_2d_array,i)) then
         write(string1,*) trim(myerrorstring),'dim/dimlen',i,dimlen,'not',size(data_2d_array,i)
         call error_handler(E_ERR,'sv_to_restart_file',string1,source,revision,revdate)
      endif
   enddo

   call vector_to_prog_var(state_vector, ivar, data_2d_array)

   call nc_check(nf90_put_var(ncid, VarID, data_2d_array), &
            'sv_to_restart_file', 'put_var '//trim(myerrorstring))

enddo

call nc_check(nf90_close(ncid), 'sv_to_restart_file', 'close '//trim(filename))

end subroutine sv_to_restart_file



!==================================================================
! The remaining interfaces come last
!==================================================================



subroutine init_interp()

! Initializes data structures needed for ncommas interpolation for
! either dipole or irregular grid.
! This should be called at static_init_model time to avoid 
! having all this temporary storage in the middle of a run.

integer :: i

! Determine whether this is a irregular lon-lat grid or a dipole.
! Do this by seeing if the lons have the same values at both
! the first and last latitude row; this is not the case for dipole. 

dipole_grid = .false.
do i = 1, nxc
   if(ulon(i, 1) /= ulon(i, nyc)) then
      dipole_grid = .true.
!     call init_dipole_interp()
      return
   endif
enddo

end subroutine init_interp


!------------------------------------------------------------


subroutine get_reg_box_indices(lon, lat, x_ind, y_ind)

! Given a longitude and latitude in degrees returns the index of the regular
! lon-lat box that contains the point.

real(r8), intent(in)  :: lon, lat
integer,  intent(out) :: x_ind, y_ind

call get_reg_lon_box(lon, x_ind)
call get_reg_lat_box(lat, y_ind)

end subroutine get_reg_box_indices

!------------------------------------------------------------

subroutine get_reg_lon_box(lon, x_ind)

! Determine which regular longitude box a longitude is in.

real(r8), intent(in)  :: lon
integer,  intent(out) :: x_ind

x_ind = int(num_reg_x * lon / 360.0_r8) + 1

! Watch out for exactly at top; assume all lats and lons in legal range
if(lon == 360.0_r8) x_ind = num_reg_x

end subroutine get_reg_lon_box

!------------------------------------------------------------

subroutine get_reg_lat_box(lat, y_ind)

! Determine which regular latitude box a latitude is in.

real(r8), intent(in)  :: lat
integer,  intent(out) :: y_ind

y_ind = int(num_reg_y * (lat + 90.0_r8) / 180.0_r8) + 1

! Watch out for exactly at top; assume all lats and lons in legal range
if(lat == 90.0_r8)  y_ind = num_reg_y

end subroutine get_reg_lat_box

!------------------------------------------------------------

subroutine reg_box_overlap(x_corners, y_corners, is_pole, reg_lon_ind, reg_lat_ind)

real(r8), intent(in)  :: x_corners(4), y_corners(4)
logical,  intent(in)  :: is_pole
integer,  intent(out) :: reg_lon_ind(2), reg_lat_ind(2)

! Find a set of regular lat lon boxes that covers all of the area covered by 
! a dipole grid qaud whose corners are given by the dimension four x_corners 
! and y_corners arrays.  The two dimensional arrays reg_lon_ind and reg_lat_ind
! return the first and last indices of the regular boxes in latitude and
! longitude respectively. These indices may wraparound for reg_lon_ind.  
! A special computation is needed for a dipole quad that has the true north 
! pole in its interior. The logical is_pole is set to true if this is the case.
! This can only happen for the t grid.  If the longitude boxes overlap 0
! degrees, the indices returned are adjusted by adding the total number of
! boxes to the second index (e.g. the indices might be 88 and 93 for a case
! with 90 longitude boxes).

real(r8) :: lat_min, lat_max, lon_min, lon_max
integer  :: i

!  A quad containing the pole is fundamentally different
if(is_pole) then
   ! Need all longitude boxes
   reg_lon_ind(1) = 1
   reg_lon_ind(2) = num_reg_x
   ! Need to cover from lowest latitude to top box
   lat_min = minval(y_corners)
   reg_lat_ind(1) = int(num_reg_y * (lat_min + 90.0_r8) / 180.0_r8) + 1
   call get_reg_lat_box(lat_min, reg_lat_ind(1))
   reg_lat_ind(2) = num_reg_y
else
   ! All other quads do not contain pole (pole could be on edge but no problem)
   ! This is specific to the dipole ncommas grids that do not go to the south pole
   ! Finding the range of latitudes is cake
   lat_min = minval(y_corners)
   lat_max = maxval(y_corners)

   ! Figure out the indices of the regular boxes for min and max lats
   call get_reg_lat_box(lat_min, reg_lat_ind(1))
   call get_reg_lat_box(lat_max, reg_lat_ind(2))

   ! Lons are much trickier. Need to make sure to wraparound the
   ! right way. There is no guarantee on direction of lons in the
   ! high latitude dipole rows.
   ! All longitudes for non-pole rows have to be within 180 degrees
   ! of one another. 
   lon_min = minval(x_corners)
   lon_max = maxval(x_corners)
   if((lon_max - lon_min) > 180.0_r8) then
      ! If the max longitude value is more than 180 
      ! degrees larger than the min, then there must be wraparound.
      ! Then, find the smallest value > 180 and the largest < 180 to get range.
      lon_min = 360.0_r8
      lon_max = 0.0_r8
      do i=1, 4
         if(x_corners(i) > 180.0_r8 .and. x_corners(i) < lon_min) lon_min = x_corners(i)
         if(x_corners(i) < 180.0_r8 .and. x_corners(i) > lon_max) lon_max = x_corners(i)
      enddo
   endif

   ! Get the indices for the extreme longitudes
   call get_reg_lon_box(lon_min, reg_lon_ind(1))
   call get_reg_lon_box(lon_max, reg_lon_ind(2))

   ! Watch for wraparound again; make sure that second index is greater than first
   if(reg_lon_ind(2) < reg_lon_ind(1)) reg_lon_ind(2) = reg_lon_ind(2) + num_reg_x
endif

end subroutine reg_box_overlap

!------------------------------------------------------------

subroutine get_quad_corners(x, i, j, corners)

real(r8), intent(in)  :: x(:, :)
integer,  intent(in)  :: i, j
real(r8), intent(out) :: corners(4)

! Grabs the corners for a given quadrilateral from the global array of lower
! right corners. Note that corners go counterclockwise around the quad.

integer :: ip1

! Have to worry about wrapping in longitude but not in latitude
ip1 = i + 1
if(ip1 > nxc) ip1 = 1

corners(1) = x(i,   j  ) 
corners(2) = x(ip1, j  )
corners(3) = x(ip1, j+1)
corners(4) = x(i,   j+1)

end subroutine get_quad_corners

!------------------------------------------------------------

subroutine update_reg_list(reg_list_num, reg_list_lon, reg_list_lat, &
   reg_lon_ind, reg_lat_ind, dipole_lon_index, dipole_lat_index)

integer, intent(inout) :: reg_list_num(:, :), reg_list_lon(:, :, :), reg_list_lat(:, :, :)
integer, intent(inout) :: reg_lon_ind(2), reg_lat_ind(2)
integer, intent(in)    :: dipole_lon_index, dipole_lat_index
 
! Updates the data structure listing dipole quads that are in a given regular box
integer :: ind_x, index_x, ind_y

! Loop through indices for each possible regular cell
! Have to watch for wraparound in longitude
if(reg_lon_ind(2) < reg_lon_ind(1)) reg_lon_ind(2) = reg_lon_ind(2) + num_reg_x

do ind_x = reg_lon_ind(1), reg_lon_ind(2)
   ! Inside loop, need to go back to wraparound indices to find right box
   index_x = ind_x
   if(index_x > num_reg_x) index_x = index_x - num_reg_x
   
   do ind_y = reg_lat_ind(1), reg_lat_ind(2)
      ! Make sure the list storage isn't full
      if(reg_list_num(index_x, ind_y) >= max_reg_list_num) then
         write(string1,*) 'max_reg_list_num (',max_reg_list_num,') is too small ... increase'
         call error_handler(E_ERR, 'update_reg_list', string1, source, revision, revdate)
      endif

      ! Increment the count
      reg_list_num(index_x, ind_y) = reg_list_num(index_x, ind_y) + 1
      ! Store this quad in the list for this regular box
      reg_list_lon(index_x, ind_y, reg_list_num(index_x, ind_y)) = dipole_lon_index
      reg_list_lat(index_x, ind_y, reg_list_num(index_x, ind_y)) = dipole_lat_index
   enddo
enddo

end subroutine update_reg_list

!------------------------------------------------------------

subroutine model_interpolate(x, location, obs_type, interp_val, istatus)

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_type
real(r8),           intent(out) :: interp_val
integer,            intent(out) :: istatus

! Model interpolate will interpolate any state variable (S, T, U, V, PSURF) to
! the given location given a state vector. The type of the variable being
! interpolated is obs_type since normally this is used to find the expected
! value of an observation at some location. The interpolated value is 
! returned in interp_val and istatus is 0 for success.

! Local storage
real(r8)       :: loc_array(3), llon, llat, lheight
integer        :: base_offset, offset, ind
integer        :: hgt_bot, hgt_top
real(r8)       :: hgt_fract
real(r8)       :: top_val, bot_val
integer        :: hstatus

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the 
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

interp_val = MISSING_R8     ! the DART bad value flag
istatus = 99                ! unknown error

! Get the individual locations values
loc_array = get_location(location)
llon    = loc_array(1)
llat    = loc_array(2)
lheight = loc_array(3)

if (debug > 1) print *, 'requesting interpolation at ', llon, llat, lheight

if( vert_is_height(location) ) then
   ! Nothing to do 
elseif ( vert_is_surface(location) ) then
   ! Nothing to do 
elseif (vert_is_level(location)) then
   ! convert the level index to an actual height 
   ind = nint(loc_array(3))
   if ( (ind < 1) .or. (ind > size(zc)) ) then 
      istatus = 11
      return
   else
      lheight = zc(ind)
   endif
else   ! if pressure or undefined, we don't know what to do
   istatus = 17
   return
endif

! Do horizontal interpolations for the appropriate levels
! Find the basic offset of this field

if(obs_type == KIND_U_WIND_COMPONENT) then
   base_offset = start_index(1)
else if(obs_type == KIND_V_WIND_COMPONENT) then
   base_offset = start_index(2)
else if(obs_type == KIND_VERTICAL_VELOCITY) then
   base_offset = start_index(3)
else if(obs_type == KIND_POTENTIAL_TEMPERATURE) then
   base_offset = start_index(4)
else if(obs_type == KIND_RADAR_REFLECTIVITY) then
   base_offset = start_index(5)
else if(obs_type == KIND_VERTICAL_VORTICITY) then
   base_offset = start_index(6)
else if(obs_type == KIND_EXNER_FUNCTION) then
   base_offset = start_index(7)
else if(obs_type == KIND_VAPOR_MIXING_RATIO) then
   base_offset = start_index(8)
else if(obs_type == KIND_CLOUDWATER_MIXING_RATIO) then
   base_offset = start_index(9)
else if(obs_type == KIND_RAINWATER_MIXING_RATIO) then
   base_offset = start_index(10)
else if(obs_type == KIND_ICE_MIXING_RATIO) then
   base_offset = start_index(11)
else if(obs_type == KIND_SNOW_MIXING_RATIO) then
   base_offset = start_index(12)
else if(obs_type == KIND_GRAUPEL_MIXING_RATIO) then
   base_offset = start_index(13)
else
   ! Not a legal type for interpolation, return istatus error
   istatus = 15
   return
endif


if (debug > 1) print *, 'base offset now ', base_offset

! surface variables are simpler
if( vert_is_surface(location) ) then
   call lon_lat_interpolate(x(base_offset:), llon, llat, obs_type, 1, interp_val, istatus)
   return
endif

! Get the bounding vertical levels and the fraction between bottom and top
call height_bounds(lheight, nzc, zc, hgt_bot, hgt_top, hgt_fract, hstatus)
if(hstatus /= 0) then
   istatus = 12
   return
endif

! Find the base location for the bottom height and interpolate horizontally 
!  on this level.  Do bottom first in case it is below the ocean floor; can
!  avoid the second horizontal interpolation.
offset = base_offset + (hgt_bot - 1) * nxc * nyc
if (debug > 1) print *, 'relative bot height offset = ', offset - base_offset
if (debug > 1) print *, 'absolute bot height offset = ', offset
call lon_lat_interpolate(x(offset:), llon, llat, obs_type, hgt_bot, bot_val, istatus)
! Failed istatus from interpolate means give up
if(istatus /= 0) return

! Find the base location for the top height and interpolate horizontally 
!  on this level.
offset = base_offset + (hgt_top - 1) * nxc * nyc
if (debug > 1) print *, 'relative top height offset = ', offset - base_offset
if (debug > 1) print *, 'absolute top height offset = ', offset
call lon_lat_interpolate(x(offset:), llon, llat, obs_type, hgt_top, top_val, istatus)
! Failed istatus from interpolate means give up
if(istatus /= 0) return


! Then weight them by the vertical fraction and return
interp_val = bot_val + hgt_fract * (top_val - bot_val)
if (debug > 1) print *, 'model_interp: interp val = ',interp_val

! All good.
istatus = 0

end subroutine model_interpolate

!------------------------------------------------------------

subroutine lon_lat_interpolate(x, lon, lat, var_type, height, interp_val, istatus)

! Subroutine to interpolate to a lon lat location given the state vector 
! for that level, x. This works just on one horizontal slice.
! NOTE: Using array sections to pass in the x array may be inefficient on some
! compiler/platform setups. Might want to pass in the entire array with a base
! offset value instead of the section if this is an issue.
! This routine works for either the dipole or a regular lat-lon grid.
! Successful interpolation returns istatus=0.

real(r8),            intent(in) :: x(:)
real(r8),            intent(in) :: lon, lat
integer,             intent(in) :: var_type, height
real(r8),           intent(out) :: interp_val
integer,            intent(out) :: istatus

! Local storage
integer  :: lat_bot, lat_top, lon_bot, lon_top, num_inds, start_ind
integer  :: x_ind, y_ind, i
real(r8) :: p(4), x_corners(4), y_corners(4), xbot, xtop
real(r8) :: lon_fract, lat_fract
logical  :: masked

if ( .not. module_initialized ) call static_init_model

! Succesful return has istatus of 0
istatus = 0

!FIXME

end subroutine lon_lat_interpolate

!------------------------------------------------------------


function get_val(lon_index, lat_index, nlon, x, var_type, height)

! Returns the value from a single level array given the lat and lon indices
! 'masked' returns true if this is NOT a valid grid location (e.g. land, or
! below the ocean floor in shallower areas).

integer,     intent(in) :: lon_index, lat_index, nlon, var_type, height
real(r8),    intent(in) :: x(:)
real(r8)                :: get_val

if ( .not. module_initialized ) call static_init_model

! Layout has lons varying most rapidly
get_val = x((lat_index - 1) * nlon + lon_index)

end function get_val

!------------------------------------------------------------

subroutine get_irreg_box(lon, lat, lon_array, lat_array, &
   found_x, found_y, lon_fract, lat_fract, istatus)

! Given a longitude and latitude of a point (lon and lat) and the
! longitudes and latitudes of the lower left corner of the regular grid
! boxes, gets the indices of the grid box that contains the point and
! the fractions along each directrion for interpolation.

real(r8),            intent(in) :: lon, lat
real(r8),            intent(in) :: lon_array(nxc, nyc), lat_array(nxc, nyc)
real(r8),           intent(out) :: lon_fract, lat_fract
integer,            intent(out) :: found_x, found_y, istatus

! Local storage
integer  :: lon_status, lat_status, lon_top, lat_top

! Succesful return has istatus of 0
istatus = 0

! Get latitude box boundaries
call lat_bounds(lat, nyc, lat_array, found_y, lat_top, lat_fract, lat_status)

! Check for error on the latitude interpolation
if(lat_status /= 0) then
   istatus = 1
   return
endif

! Find out what longitude box and fraction
call lon_bounds(lon, nxc, lon_array, found_x, lon_top, lon_fract)

end subroutine get_irreg_box

!------------------------------------------------------------

subroutine lon_bounds(lon, nlons, lon_array, bot, top, fract)

! Given a longitude lon, the array of longitudes for grid boundaries, and the
! number of longitudes in the grid, returns the indices of the longitude
! below and above the location longitude and the fraction of the distance
! between. It is assumed that the longitude wraps around for a global grid. 
! Since longitude grids are going to be regularly spaced, this could be made more efficient.
! Algorithm fails for a silly grid that has only two longitudes separated by 180 degrees.

real(r8),          intent(in) :: lon
integer,           intent(in) :: nlons
real(r8),          intent(in) :: lon_array(:, :)
integer,          intent(out) :: bot, top
real(r8),         intent(out) :: fract

! Local storage
integer  :: i
real(r8) :: dist_bot, dist_top

if ( .not. module_initialized ) call static_init_model

! This is inefficient, someone could clean it up since longitudes are regularly spaced
! But note that they don't have to start at 0
do i = 2, nlons
   dist_bot = lon_dist(lon, lon_array(i - 1, 1))
   dist_top = lon_dist(lon, lon_array(i, 1))
   if(dist_bot <= 0 .and. dist_top > 0) then
      bot = i - 1
      top = i
      fract = abs(dist_bot) / (abs(dist_bot) + dist_top)
      return
   endif
enddo

! Falling off the end means it's in between; wraparound
bot = nlons
top = 1
dist_bot = lon_dist(lon, lon_array(bot, 1))
dist_top = lon_dist(lon, lon_array(top, 1)) 
fract = abs(dist_bot) / (abs(dist_bot) + dist_top)

end subroutine lon_bounds

!-------------------------------------------------------------

subroutine lat_bounds(lat, nlats, lat_array, bot, top, fract, istatus)

! Given a latitude lat, the array of latitudes for grid boundaries, and the
! number of latitudes in the grid, returns the indices of the latitude
! below and above the location latitude and the fraction of the distance
! between. istatus is returned as 0 unless the location latitude is 
! south of the southernmost grid point (1 returned) or north of the 
! northernmost (2 returned). If one really had lots of polar obs would 
! want to worry about interpolating around poles.

real(r8),          intent(in) :: lat
integer,           intent(in) :: nlats
real(r8),          intent(in) :: lat_array(:, :)
integer,          intent(out) :: bot, top
real(r8),         intent(out) :: fract
integer,          intent(out) :: istatus

! Local storage
integer :: i

if ( .not. module_initialized ) call static_init_model

! Success should return 0, failure a positive number.
istatus = 0

! Check for too far south or north
if(lat < lat_array(1, 1)) then
   istatus = 1
   return
else if(lat > lat_array(1, nlats)) then
   istatus = 2
   return
endif

! In the middle, search through
do i = 2, nlats
   if(lat <= lat_array(1, i)) then
      bot = i - 1
      top = i
      fract = (lat - lat_array(1, bot)) / (lat_array(1, top) - lat_array(1, bot))
      return
   endif
enddo

! Shouldn't get here. Might want to fail really hard through error handler
istatus = 40

end subroutine lat_bounds

!------------------------------------------------------------

function lon_dist(lon1, lon2)

! Returns the smallest signed distance between lon1 and lon2 on the sphere
! If lon1 is less than 180 degrees east of lon2 the distance is negative
! If lon1 is less than 180 degrees west of lon2 the distance is positive

real(r8), intent(in) :: lon1, lon2
real(r8)             :: lon_dist

if ( .not. module_initialized ) call static_init_model

lon_dist = lon2 - lon1
if(lon_dist >= -180.0_r8 .and. lon_dist <= 180.0_r8) then
   return
else if(lon_dist < -180.0_r8) then
   lon_dist = lon_dist + 360.0_r8
else
   lon_dist = lon_dist - 360.0_r8
endif

end function lon_dist

!------------------------------------------------------------

subroutine get_dipole_quad(lon, lat, qlons, qlats, num_inds, start_ind, &
   x_inds, y_inds, found_x, found_y, istatus)

real(r8), intent(in)  :: lon, lat, qlons(:, :), qlats(:, :)
integer,  intent(in)  :: num_inds, start_ind, x_inds(:), y_inds(:)
integer,  intent(out) :: found_x, found_y, istatus

! Given the lon and lat of a point, and a list of the
! indices of the quads that might contain a point at (lon, lat), determines
! which quad contains the point.  istatus is returned as 0 if all went 
! well and 1 if the point was not found to be in any of the quads.

integer :: i, my_index
real(r8) :: x_corners(4), y_corners(4)

istatus = 0

! Loop through all the quads and see if the point is inside
do i = 1, num_inds
   my_index = start_ind + i - 1
   call get_quad_corners(qlons, x_inds(my_index), y_inds(my_index), x_corners)
   call get_quad_corners(qlats, x_inds(my_index), y_inds(my_index), y_corners)

   ! Ssearch in this individual quad
   if(in_quad(lon, lat, x_corners, y_corners)) then
      found_x = x_inds(my_index)
      found_y = y_inds(my_index)
      return
   endif
enddo

! Falling off the end means search failed, return istatus 1
istatus = 1

end subroutine get_dipole_quad

!------------------------------------------------------------

function in_quad(lon, lat, x_corners, y_corners)

! Return in_quad true if the point (lon, lat) is in the quad with 
! the given corners.
! Do this by line tracing in latitude for now. For non-pole point, want a vertical
! line from the lon, lat point to intersect a side of the quad both above
! and below the point.

real(r8), intent(in)  :: lon, lat, x_corners(4), y_corners(4)
logical               :: in_quad

real(r8) :: x(2), y(2)
logical  :: cant_be_in_box, in_box
integer  :: intercepts_above(4), intercepts_below(4), i
integer  :: num_above, num_below

! Default answer is point is not in quad
in_quad = .false.

! Loop through the sides and compute intercept (if any) with a vertical line
! from the point. This line has equation x=lon.
do i = 1, 4
   ! Load up the sides endpoints
   if(i <= 3) then
      x(1:2) = x_corners(i:i+1)
      y(1:2) = y_corners(i:i+1)
   else
      x(1) = x_corners(4)
      x(2) = x_corners(1)
      y(1) = y_corners(4)
      y(2) = y_corners(1)
   endif

   ! Check to see how a vertical line from the point is related to this side
   call line_intercept(x, y, lon, lat, cant_be_in_box, in_box, intercepts_above(i), &
      intercepts_below(i))

   ! If cant_be_in_box is true, can return right away
   if(cant_be_in_box) then
      in_quad = .false.
      return
   ! Return true if on a side
   else if(in_box) then
      in_quad = .true.
      return
   endif

enddo

! See if the line intercepted a side of the quad both above and below
num_above = sum(intercepts_above)
num_below = sum(intercepts_below)

if(num_above > 0 .and. num_below > 0) then
   in_quad = .true.
endif

end function in_quad

!------------------------------------------------------------

subroutine line_intercept(side_x_in, side_y, x_point_in, y_point, &
   cant_be_in_box, in_box, intercept_above, intercept_below)

real(r8), intent(in)  :: side_x_in(2), side_y(2), x_point_in, y_point
logical,  intent(out) :: cant_be_in_box, in_box
integer,  intent(out) :: intercept_above, intercept_below

! Find the intercept of a vertical line from point (x_point, y_point) and
! a line segment with endpoints side_x and side_y.
! For a given side have endpoints (side_x1, side_y1) and (side_x2, side_y2)
! so equation of segment is y = side_y1 + m(x-side_x1) for y 
! between side_y1 and side_y2.
! Intersection of vertical line and line containing side 
! occurs at y = side_y1 + m(x_point - side_x1); need this
! y to be between side_y1 and side_y2.
! If the vertical line is colinear with the side but the point is not on the side, return
! cant_be_in_box as true. If the point is on the side, return in_box true.
! If the intersection of the vertical line and the side occurs at a point above
! the given point, return 1 for intercept_above. If the intersection occurs 
! below, return 1 for intercept_below. If the vertical line does not intersect
! the segment, return false and 0 for all intent out arguments.

! WARNING: CERTAINLY A PROBLEM FOR THE POLE BOX!!! POLE BOX COULD
! HAVE SIDES THAT ARE LONGER THAN 180. For now pole boxes are excluded.

! This can probably be made much cleaner and more efficient.

real(r8) :: slope, y_intercept, side_x(2), x_point

! May have to adjust the longitude intent in values, so copy
side_x = side_x_in
x_point = x_point_in

! See if the side wraps around in longitude
if(maxval(side_x) - minval(side_x) > 180.0_r8) then
   if(side_x(1) < 180.0_r8)  side_x(1) =  side_x(1) + 360.0_r8
   if(side_x(2) < 180.0_r8)  side_x(2) =  side_x(2) + 360.0_r8
   if(x_point < 180.0_r8) x_point = x_point + 360.0_r8
endif

! Initialize the default returns 
cant_be_in_box   = .false.
in_box           = .false.
intercept_above = 0
intercept_below = 0

! First easy check, if x_point is not between endpoints of segment doesn't intersect
if(x_point < minval(side_x) .or. x_point > maxval(side_x)) return

! Otherwise line must intersect the segment

! First subblock, slope is undefined
if(side_x(2) == side_x(1)) then
   ! The line is colinear with the side
   ! If y_point is between endpoints then point is on this side
   if(y_point <= maxval(side_y) .and. y_point >= minval(side_y)) then
      in_box = .true.
      return
   ! If not on side but colinear with side, point cant be in quad
   else
      cant_be_in_box = .true.
      return
   endif

else

   ! Second possibility; slope is defined
   ! FIXME: watch out for numerical instability.
   ! near-zero x's and large side_y's may cause overflow
   slope = (side_y(2) - side_y(1)) / (side_x(2) - side_x(1))

   ! Intercept of vertical line through is at x_point and...
   y_intercept = side_y(1) + slope * (x_point - side_x(1))

   ! Intersects the segment, is it above, below, or at the point
   if(y_intercept == y_point) then
      in_box = .true.
      return
   else if(y_intercept > y_point) then
      intercept_above = 1
      return
   else
      intercept_below = 1
      return
   endif
endif

end subroutine line_intercept

!------------------------------------------------------------

subroutine quad_bilinear_interp(lon_in, lat, x_corners_in, y_corners, &
   p, interp_val)

real(r8), intent(in)  :: lon_in, lat, x_corners_in(4), y_corners(4), p(4)
real(r8), intent(out) :: interp_val

! Given a longitude and latitude (lon_in, lat), the longitude and
! latitude of the 4 corners of a quadrilateral and the values at the
! four corners, interpolates to (lon_in, lat) which is assumed to
! be in the quad. This is done by bilinear interpolation, fitting
! a function of the form a + bx + cy + dxy to the four points and 
! then evaluating this function at (lon, lat). The fit is done by
! solving the 4x4 system of equations for a, b, c, and d. The system
! is reduced to a 3x3 by eliminating a from the first three equations
! and then solving the 3x3 before back substituting. There is concern
! about the numerical stability of this implementation. Implementation
! checks showed accuracy to seven decimal places on all tests.

integer :: i
real(r8) :: m(3, 3), v(3), r(3), a, x_corners(4), lon, lon_mean

! Watch out for wraparound on x_corners.
lon = lon_in
x_corners = x_corners_in

! See if the side wraps around in longitude. If the corners longitudes
! wrap around 360, then the corners and the point to interpolate to
! must be adjusted to be in the range from 180 to 540 degrees.
if(maxval(x_corners) - minval(x_corners) > 180.0_r8) then
   if(lon < 180.0_r8) lon = lon + 360.0_r8
   do i = 1, 4
      if(x_corners(i) < 180.0_r8) x_corners(i) = x_corners(i) + 360.0_r8
   enddo
endif


!*******
! Problems with extremes in polar cell interpolation can be reduced
! by this block, but it is not clear that it is needed for actual
! ocean grid data
! Find the mean longitude of corners and remove
!!!lon_mean = sum(x_corners) / 4.0_r8
!!!x_corners = x_corners - lon_mean
!!!lon = lon - lon_mean
! Multiply everybody by the cos of the latitude
!!!do i = 1, 4
   !!!x_corners(i) = x_corners(i) * cos(y_corners(i) * deg2rad)
!!!enddo
!!!lon = lon * cos(lat * deg2rad)

!*******


! Fit a surface and interpolate; solve for 3x3 matrix
do i = 1, 3
   ! Eliminate a from the first 3 equations
   m(i, 1) = x_corners(i) - x_corners(i + 1)
   m(i, 2) = y_corners(i) - y_corners(i + 1)
   m(i, 3) = x_corners(i)*y_corners(i) - x_corners(i + 1)*y_corners(i + 1)
   v(i) = p(i) - p(i + 1)
enddo

! Solve the matrix for b, c and d
call mat3x3(m, v, r)

! r contains b, c, and d; solve for a
a = p(4) - r(1) * x_corners(4) - r(2) * y_corners(4) - &
   r(3) * x_corners(4)*y_corners(4)


!----------------- Implementation test block
! When interpolating on dipole x3 never exceeded 1e-9 error in this test
!!!do i = 1, 4
   !!!interp_val = a + r(1)*x_corners(i) + r(2)*y_corners(i)+ r(3)*x_corners(i)*y_corners(i)
   !!!if(abs(interp_val - p(i)) > 1e-9) then
      !!!write(*, *) 'large interp residual ', interp_val - p(i)
   !!!endif
!!!enddo
!----------------- Implementation test block


! Now do the interpolation
interp_val = a + r(1)*lon + r(2)*lat + r(3)*lon*lat

!********
! Avoid exceeding maxima or minima as stopgap for poles problem
! When doing bilinear interpolation in quadrangle, can get interpolated
! values that are outside the range of the corner values
if(interp_val > maxval(p)) then 
   interp_val = maxval(p)
else if(interp_val < minval(p)) then
   interp_val = minval(p)
endif
!********

end subroutine quad_bilinear_interp

!------------------------------------------------------------

subroutine mat3x3(m, v, r)

real(r8), intent(in)  :: m(3, 3), v(3)
real(r8), intent(out) :: r(3)

! Solves rank 3 linear system mr = v for r
! using Cramer's rule. This isn't the best choice
! for speed or numerical stability so might want to replace
! this at some point.

real(r8) :: m_sub(3, 3), numer, denom
integer  :: i

! Compute the denominator, det(m)
denom = deter3(m)

! Loop to compute the numerator for each component of r
do i = 1, 3
   m_sub = m
   m_sub(:, i) = v   
   numer = deter3(m_sub)
   r(i) = numer / denom
enddo

end subroutine mat3x3

!------------------------------------------------------------
function deter3(m)

real(r8) :: deter3
real(r8), intent(in) :: m(3, 3)

! Computes determinant of 3x3 matrix m

deter3 = m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) + &
         m(1,3)*m(2,1)*m(3,2) - m(3,1)*m(2,2)*m(1,3) - &
         m(1,1)*m(2,3)*m(3,2) - m(3,3)*m(2,1)*m(1,2)

end function deter3

!------------------------------------------------------------

subroutine height_bounds(lheight, nheights, hgt_array, bot, top, fract, istatus)

real(r8), intent( in) :: lheight
integer,  intent( in) :: nheights
real(r8), intent( in) :: hgt_array(nheights)
integer,  intent(out) :: bot, top
real(r8), intent(out) :: fract
integer,  intent(out) :: istatus

! Local variables
integer   :: i

if ( .not. module_initialized ) call static_init_model

! Succesful istatus is 0
! Make any failure here return istatus in the 20s
istatus = 0

! The zc array contains the heights of the center of the vertical grid boxes
! In this case (unlike how we handle the MIT heights), positive is really down.
! FIXME: in the MIT model, we're given box widths and we compute the centers,
! and we computed them with larger negative numbers being deeper.  Here,
! larger positive numbers are deeper.

! It is assumed that the top box is shallow and any observations shallower
! than the height of this box's center are just given the value of the
! top box.
if(lheight <= hgt_array(1)) then
   top = 1
   bot = 2
   ! NOTE: the fract definition is the relative distance from bottom to top
   fract = 1.0_r8 
if (debug > 1) print *, 'above first level in height'
if (debug > 1) print *, 'hgt_array, top, bot, fract=', hgt_array(1), top, bot, fract
   return
endif

! Search through the boxes
do i = 2, nheights
   ! If the location is shallower than this entry, it must be in this box
   if(lheight < hgt_array(i)) then
      top = i -1
      bot = i
      fract = (hgt_array(bot) - lheight) / (hgt_array(bot) - hgt_array(top))
if (debug > 1) print *, 'i, hgt_array, top, bot, fract=', i, hgt_array(i), top, bot, fract
      return
   endif
enddo

! Falling off the end means the location is lower than the deepest height
istatus = 20

end subroutine height_bounds



subroutine get_state_indices(index_in, lat_index, lon_index, height_index, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated array indices for lat, lon, and height, as well as the type.

integer, intent(in)  :: index_in
integer, intent(out) :: lat_index, lon_index, height_index
integer, intent(out) :: var_type

integer :: startind, offset, var_index

if ( .not. module_initialized ) call static_init_model

if (debug > 5) print *, 'asking for meta data about index ', index_in

call get_state_kind(index_in, var_index, var_type, startind, offset)

! FIXME:  this should be using progvar(var_index)%numdims, %dimlens(:), %index1, etc.

! wrong.
height_index = 1

! old code
!lat_index = (offset - ((height_index-1)*nxc*nyc)) / nxc + 1
!lon_index =  offset - ((height_index-1)*nxc*nyc) - ((lat_index-1)*nxc) + 1
lat_index = 1
lon_index = 1

if (debug > 5) print *, 'lon, lat, height index = ', lon_index, lat_index, height_index

end subroutine get_state_indices



subroutine get_state_kind(index_in, var_index, var_type, startind, offset)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the kind,
! and both the starting offset for this kind, as well as the offset into
! the block of this kind.

integer, intent(in)  :: index_in
integer, intent(out) :: var_index, var_type, startind, offset

integer :: i

if ( .not. module_initialized ) call static_init_model

if (debug > 5) print *, 'asking for meta data about index ', index_in

do i = 1, nfields
   if ((index_in >= progvar(i)%index1) .and. (index_in <= progvar(i)%indexN)) then
      var_index = i
      var_type = progvar(i)%dart_kind
      startind = progvar(i)%index1
      offset = index_in - startind

      if (debug > 5) print *, 'var type = ', var_type
      if (debug > 5) print *, 'startind = ', startind
      if (debug > 5) print *, 'offset = ', offset

      return
   endif
enddo

! shouldn't get here.

end subroutine get_state_kind



subroutine vector_to_2d_prog_var(x, varindex, data_2d_array)
!------------------------------------------------------------------
! convert the values from a 1d fortran array, starting at an offset,
! into a 2d fortran array.  the 2 dims are taken from the array size.
!
real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: varindex
real(r8), dimension(:,:), intent(out) :: data_2d_array

integer :: i,j,ii
integer :: dim1,dim2
character(len=128) :: varname

if ( .not. module_initialized ) call static_init_model

dim1 = size(data_2d_array,1)
dim2 = size(data_2d_array,2)

varname = progvarnames(varindex)

if (dim1 /= nxc) then
   write(string1,*)trim(varname),' 2d array dim 1 ',dim1,' /= ',nxc
   call error_handler(E_ERR,'vector_to_2d_prog_var',string1,source,revision,revdate) 
endif
if (dim2 /= nyc) then
   write(string1,*)trim(varname),' 2d array dim 2 ',dim2,' /= ',nyc
   call error_handler(E_ERR,'vector_to_2d_prog_var',string1,source,revision,revdate) 
endif

ii = start_index(varindex)

do j = 1,nyc   ! latitudes
do i = 1,nxc   ! longitudes
   data_2d_array(i,j) = x(ii)
   ii = ii + 1
enddo
enddo

end subroutine vector_to_2d_prog_var



subroutine vector_to_3d_prog_var(x, varindex, data_3d_array)
!------------------------------------------------------------------
! convert the values from a 1d fortran array, starting at an offset,
! into a 3d fortran array.  the 3 dims are taken from the array size.
!
real(r8), dimension(:),     intent(in)  :: x
integer,                    intent(in)  :: varindex
real(r8), dimension(:,:,:), intent(out) :: data_3d_array

integer :: i,j,k,ii
integer :: dim1,dim2,dim3
character(len=128) :: varname

if ( .not. module_initialized ) call static_init_model

dim1 = size(data_3d_array,1)
dim2 = size(data_3d_array,2)
dim3 = size(data_3d_array,3)

varname = progvarnames(varindex)

if (dim1 /= nxc) then
   write(string1,*)trim(varname),' 3d array dim 1 ',dim1,' /= ',nxc
   call error_handler(E_ERR,'vector_to_3d_prog_var',string1,source,revision,revdate) 
endif
if (dim2 /= nyc) then
   write(string1,*)trim(varname),' 3d array dim 2 ',dim2,' /= ',nyc
   call error_handler(E_ERR,'vector_to_3d_prog_var',string1,source,revision,revdate) 
endif
if (dim3 /= nzc) then
   write(string1,*)trim(varname),' 3d array dim 3 ',dim3,' /= ',nzc
   call error_handler(E_ERR,'vector_to_3d_prog_var',string1,source,revision,revdate) 
endif

ii = start_index(varindex)

do k = 1,nzc   ! vertical
do j = 1,nyc   ! latitudes
do i = 1,nxc   ! longitudes
   data_3d_array(i,j,k) = x(ii)
   ii = ii + 1
enddo
enddo
enddo

end subroutine vector_to_3d_prog_var



  function is_on_ugrid(obs_type)
!------------------------------------------------------------------
!  returns true if U 
integer, intent(in) :: obs_type
logical             :: is_on_ugrid

if ( .not. module_initialized ) call static_init_model

is_on_ugrid = .FALSE.

if (obs_type == KIND_U_WIND_COMPONENT) is_on_ugrid = .TRUE.

end function is_on_ugrid



  function is_on_vgrid(obs_type)
!------------------------------------------------------------------
!  returns true if V
integer, intent(in) :: obs_type
logical             :: is_on_vgrid

if ( .not. module_initialized ) call static_init_model

is_on_vgrid = .FALSE.

if (obs_type == KIND_V_WIND_COMPONENT) is_on_vgrid = .TRUE.

end function is_on_vgrid



  function is_on_wgrid(obs_type)
!------------------------------------------------------------------
!  returns true if W
integer, intent(in) :: obs_type
logical             :: is_on_wgrid

if ( .not. module_initialized ) call static_init_model

is_on_wgrid = .FALSE.

if (obs_type == KIND_VERTICAL_VELOCITY) is_on_wgrid = .TRUE.

end function is_on_wgrid



!===================================================================
! End of model_mod
!===================================================================
end module model_mod


