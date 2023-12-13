! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This is the interface between the ncommas model and DART.

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, digits12, SECPERDAY, MISSING_R8,          &
                             rad2deg, deg2rad, PI, obstypelength

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
                             do_output, to_upper,                              &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, file_exist, find_textfile_dims,        &
                             file_to_text

use     obs_kind_mod, only : QTY_U_WIND_COMPONENT,   &
                             QTY_V_WIND_COMPONENT,   &
                             QTY_VERTICAL_VELOCITY,  &
                             get_index_for_quantity

use mpi_utilities_mod, only: my_task_id

use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use    netcdf_utilities_mod, only : nc_check

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
          get_ncommas_restart_filename, &
          get_base_time,                &
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
logical            :: output_state_vector = .false.
integer            :: debug = 0   ! turn up for more and more debug messages
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: ncommas_restart_filename = 'ncommas_restart.nc'

namelist /model_nml/  &
   ncommas_restart_filename,    &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   calendar,                    &
   debug

!------------------------------------------------------------------
!
!  The DART state vector may consist of things like:  
!
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
!  The variables in the ncommas restart file that are used to create the 
!  DART state vector are specified in the input.nml:ncommas_vars_nml namelist.
!
!------------------------------------------------------------------

integer, parameter :: max_state_variables = 80
integer, parameter :: num_state_table_columns = 2
character(len=NF90_MAX_NAME) :: ncommas_state_variables(max_state_variables * num_state_table_columns ) = ' '
character(len=NF90_MAX_NAME) :: variable_table(max_state_variables, num_state_table_columns )

namelist /ncommas_vars_nml/ ncommas_state_variables

integer :: nfields

! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=NF90_MAX_NAME) :: storder
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer :: originalnumdims
   integer :: posdef
   integer :: numdims
   integer :: varsize     ! prod(dimlens(1:numdims))
   integer :: index1      ! location in dart state vector of first occurrence
   integer :: indexN      ! location in dart state vector of last  occurrence
   integer :: dart_kind
   character(len=obstypelength) :: kind_string
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

! little stuff here and there..

real(digits12), parameter :: rearth=1000.0_digits12 * 6367.0_digits12 ! radius of earth (m)

! Grid parameters - the values will be read from an ncommas restart file. 
! Each spatial dimension has a staggered counterpart.

integer :: nxc=-1, nyc=-1, nzc=-1    ! scalar grid positions
integer :: nxe=-1, nye=-1, nze=-1    ! staggered grid positions

! Global storage for the grid's reference lat/lon

real(r8) :: ref_lat, ref_lon
real(r8) :: xg_pos, yg_pos, hgt_offset

! locations of cell centers (C) and edges (E) for each axis.

real(r8), allocatable :: XC(:), XE(:)
real(r8), allocatable :: YC(:), YE(:)
real(r8), allocatable :: ZC(:), ZE(:)

! These arrays store the longitude and latitude of the lower left corner

real(r8), allocatable :: ULAT(:,:), ULON(:,:)  ! XE,YC
real(r8), allocatable :: VLAT(:,:), VLON(:,:)  ! XC,YE
real(r8), allocatable :: WLAT(:,:), WLON(:,:)  ! XC,YC

integer               :: model_size      ! the state vector length
type(time_type)       :: model_time      ! valid time of the model state
type(time_type)       :: model_timestep  ! smallest time to adv model
real(r8), allocatable :: ens_mean(:)     ! may be needed for forward ops

! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.

integer :: print_every_Nth  = 10000

!------------------------------------------------------------------
! The ncommas restart manager namelist variables
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



!############################################################################
!
!     ##################################################################
!     ######                                                      ######
!     ######            SUBROUTINE GET_STATE_META_DATA            ######
!     ######                                                      ######
!     ##################################################################
!
!     PURPOSE:
!
!     Given an integer index into the state vector structure, returns the
!     associated array indices for lat, lon, and height, as well as the type.
!
!############################################################################
!
!     Author:  Tim Hoar, Ted Mansell, Lou Wicker
!
!     Creation Date:  August 2010
!     
!     Variables needed to be stored in the MODEL_MODULE data structure
!
!       PROGVAR => Module's data structure
!
!############################################################################

subroutine get_state_meta_data(index_in, location, var_type)

! Passed variables

  integer, intent(in)            :: index_in
  type(location_type)            :: location
  integer, optional, intent(out) :: var_type
  
! Local variables

  integer  :: nxp, nyp, nzp, iloc, jloc, kloc, nf, n
  integer  :: myindx, lat_index, lon_index, index2
  real(r8) :: height
  real(r8) :: x1,y1

  if ( .not. module_initialized ) call static_init_model
  
  myindx = -1
  nf     = -1
  
  FindIndex : DO n = 1,nfields
    IF( (progvar(n)%index1 <= index_in) .and. (index_in <= progvar(n)%indexN) ) THEN
      nf = n
      myindx = index_in - progvar(n)%index1 + 1
      EXIT FindIndex
    ENDIF
  ENDDO FindIndex 
  
  IF( myindx == -1 ) THEN
     write(string1,*) 'Problem, cannot find base_offset, index_in is: ', index_in
     call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
  ENDIF
  
  nxp = progvar(nf)%dimlens(1)
  nyp = progvar(nf)%dimlens(2)
  nzp = progvar(nf)%dimlens(3)

  index2 = myindx
  kloc   = 1 + (myindx-1) / (nxp*nyp)
  myindx = myindx - (kloc-1)*nyp*nxp
  jloc   = 1 + (myindx-1) / nxp
  myindx = myindx - (jloc-1)*nxp
  iloc   = myindx
  
  lat_index = jloc
  lon_index = iloc

  IF ( progvar(nf)%dart_kind == QTY_VERTICAL_VELOCITY ) THEN
        height = ze(kloc)
  ELSE
        height = zc(kloc)
  ENDIF
  
  IF (debug > 5) THEN

    IF ( progvar(nf)%dart_kind == QTY_U_WIND_COMPONENT ) THEN
       x1 = xe(lon_index)
    ELSE
       x1 = xc(lon_index)
    ENDIF
    
    IF ( progvar(nf)%dart_kind == QTY_V_WIND_COMPONENT ) THEN
       y1 = ye(lat_index)
    ELSE
       y1 = yc(lat_index)
    ENDIF

    write(*,FMT='("INDEX_IN / INDEX / NVAR / NXP, NYP, NZP: ",2(i10,2x),4(i5,2x))') index_in, index2, nf, nxp, nyp, nzp
    write(*,FMT='("                       ILOC, JLOC, KLOC: ",3(i5,2x))') lon_index, lat_index, kloc
    write(*,FMT='("                                  X/Y/Z: ",3(f10.1,2x))') x1,y1, height

  ENDIF
  
! Here we assume:
! That anything not a velocity is zone centered.
  
  IF(progvar(nf)%dart_kind == QTY_U_WIND_COMPONENT) THEN
    location = set_location(ulon(iloc,jloc), ulat(iloc,jloc), height, VERTISHEIGHT)
  ELSEIF(progvar(nf)%dart_kind == QTY_V_WIND_COMPONENT) THEN
    location = set_location(vlon(iloc,jloc), vlat(iloc,jloc), height, VERTISHEIGHT)
  ELSEIF(progvar(nf)%dart_kind == QTY_VERTICAL_VELOCITY) THEN
    height   = ze(kloc)
    location = set_location(wlon(iloc,jloc), wlat(iloc,jloc), height, VERTISHEIGHT)
  ELSE
    location = set_location(wlon(iloc,jloc), wlat(iloc,jloc), height, VERTISHEIGHT)
  ENDIF
  
  IF (present(var_type)) THEN
     var_type = progvar(nf)%dart_kind
  ENDIF
  
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
! the dart_ncommas_mod module.

! Local variables - all the important ones have module scope

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname
character(len=obstypelength)          :: kind_string
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

! Read the NCOMMAS variable list to populate DART state vector
! Once parsed, the values will be recorded for posterity
call find_namelist_in_file('ncommas_vars.nml', 'ncommas_vars_nml', iunit)
read(iunit, nml = ncommas_vars_nml, iostat = io)
call check_namelist_read(iunit, io, 'ncommas_vars_nml')

!---------------------------------------------------------------
! Set the time step ... causes ncommas namelists to be read.
! Ensures model_timestep is multiple of 'dynamics_timestep'

call set_calendar_type( calendar )   ! comes from model_mod_nml

model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)

!---------------------------------------------------------------
! 1) get grid dimensions
! 2) allocate space for the grids 
! 3) read them, convert them from X-Y-Z to lat-lon-z

call get_grid_dims(nxc, nxe, nyc, nye, nzc, nze )

allocate(ULAT(nxe,nyc), ULON(nxe,nyc))
allocate(VLAT(nxc,nye), VLON(nxc,nye))
allocate(WLAT(nxc,nyc), WLON(nxc,nyc))
allocate(  XC(  nxc  ),   XE(  nxe  ))
allocate(  YC(  nyc  ),   YE(  nye  ))
allocate(  ZC(  nzc  ),   ZE(  nze  ))

CALL GET_GRID(nxc, nxe, nyc, nye, nzc, nze,       &
              XC, XE, YC, YE, ZC, ZE,             &
              ULAT, ULON, VLAT, VLON, WLAT, WLON, &
              xg_pos, yg_pos, ref_lat, ref_lon, hgt_offset)
              
!---------------------------------------------------------------
! Compile the list of ncommas variables to use in the creation
! of the DART state vector. Required to determine model_size.
!
! Verify all variables are in the ncommas restart file
!
! Compute the offsets into the state vector for the start of each
! different variable type. Requires reading shapes from the NCOMMAS
! restart file. As long as TIME is the LAST dimension, we're OK.
!
! Record the extent of the data type in the state vector.

call nc_check( nf90_open(trim(ncommas_restart_filename), NF90_NOWRITE, ncid), &
                  'static_init_model', 'open '//trim(ncommas_restart_filename))

call verify_state_variables( ncommas_state_variables, ncid, ncommas_restart_filename, &
                             nfields, variable_table )

TimeDimID = FindTimeDimension( ncid )

if (TimeDimID < 0 ) then
   write(string1,*)'unable to find a dimension named TIME Time or time.'
   call error_handler(E_MSG,'static_init_model', string1, source, revision, revdate)
endif

call nc_check(nf90_Inquire(ncid,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                    'static_init_model', 'inquire '//trim(ncommas_restart_filename))

if ( (TimeDimID > 0) .and. (unlimitedDimID > 0) .and. (TimeDimID /= unlimitedDimID)) then
   write(string1,*)'IF TIME is not the unlimited dimension, I am lost.'
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
   progvar(ivar)%dimlens     = 0

   string2 = trim(ncommas_restart_filename)//' '//trim(varname)

   call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), &
            'static_init_model', 'inq_varid '//trim(string2))

   call nc_check( nf90_get_att(ncid, VarId, 'type' , progvar(ivar)%storder), &
            'static_init_model', 'get_att type '//trim(string2))

   call nc_check( nf90_get_att(ncid, VarId, 'long_name' , progvar(ivar)%long_name), &
            'static_init_model', 'get_att long_name '//trim(string2))

   call nc_check( nf90_get_att(ncid, VarId, 'posdef' , progvar(ivar)%posdef), &
            'static_init_model', 'get_att posdef '//trim(string2))

   call nc_check( nf90_get_att(ncid, VarId, 'units' , progvar(ivar)%units), &
            'static_init_model', 'get_att units '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid, VarId, dimids=dimIDs, ndims=numdims), &
            'static_init_model', 'inquire '//trim(string2))
   progvar(ivar)%originalnumdims = numdims

   ! TIME is the unlimited dimension, and in Fortran, this is always the
   ! LAST dimension in this loop. Since we are not concerned with it, we
   ! need to skip it. 
   varsize = 1
   dimlen  = 1
   progvar(ivar)%numdims = 0
   DimensionLoop : do i = 1,numdims

      if (dimIDs(i) == TimeDimID) cycle DimensionLoop
      
      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
                                          'static_init_model', string1)

      progvar(ivar)%numdims = progvar(ivar)%numdims + 1 
      progvar(ivar)%dimlens(i) = dimlen
      varsize = varsize * dimlen

   enddo DimensionLoop

   ! Make sure the storage order is as we expect.
   shapeok = .false.
   if (progvar(ivar)%numdims == 1) then
      if (trim(progvar(ivar)%storder) ==   'x1d') shapeok = .true.
      if (trim(progvar(ivar)%storder) ==   'y1d') shapeok = .true.
      if (trim(progvar(ivar)%storder) ==   'z1d') shapeok = .true.
   elseif (progvar(ivar)%numdims == 2) then
      if (trim(progvar(ivar)%storder) ==  'xy2d') shapeok = .true.
   elseif (progvar(ivar)%numdims == 3) then
      if (trim(progvar(ivar)%storder) == 'xyz3d') shapeok = .true.
   endif
   if ( .not. shapeok ) then
      write(string1,*)'unable to handle storage order of '//trim(progvar(ivar)%storder)//' numdims = ',progvar(ivar)%numdims
      call error_handler(E_ERR,'static_init_model', string1, source, revision, revdate, &
                                              text2=string2)
   endif

   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1 
   index1                    = index1 + varsize      ! sets up for next variable

   if ( debug > 0 ) then
      write(logfileunit,*)
      write(logfileunit,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(logfileunit,*) '  type        ',trim(progvar(ivar)%storder)
      write(logfileunit,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(logfileunit,*) '  units       ',trim(progvar(ivar)%units)
      write(logfileunit,*) '  orgnalndims ',progvar(ivar)%originalnumdims
      write(logfileunit,*) '  numdims     ',progvar(ivar)%numdims
      write(logfileunit,*) '  dimlens     ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      write(logfileunit,*) '  varsize     ',progvar(ivar)%varsize
      write(logfileunit,*) '  index1      ',progvar(ivar)%index1
      write(logfileunit,*) '  indexN      ',progvar(ivar)%indexN
      write(logfileunit,*) '  dart_kind   ',progvar(ivar)%dart_kind
      write(logfileunit,*) '  kind_string ',progvar(ivar)%kind_string

      write(     *     ,*)
      write(     *     ,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(     *     ,*) '  type        ',trim(progvar(ivar)%storder)
      write(     *     ,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(     *     ,*) '  units       ',trim(progvar(ivar)%units)
      write(     *     ,*) '  orgnalndims ',progvar(ivar)%originalnumdims
      write(     *     ,*) '  numdims     ',progvar(ivar)%numdims
      write(     *     ,*) '  dimlens     ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      write(     *     ,*) '  varsize     ',progvar(ivar)%varsize
      write(     *     ,*) '  index1      ',progvar(ivar)%index1
      write(     *     ,*) '  indexN      ',progvar(ivar)%indexN
      write(     *     ,*) '  dart_kind   ',progvar(ivar)%dart_kind
      write(     *     ,*) '  kind_string ',progvar(ivar)%kind_string
   endif

enddo

call nc_check( nf90_close(ncid), &
                  'static_init_model', 'close '//trim(ncommas_restart_filename))

model_size = progvar(nfields)%indexN

if ( debug > 0 ) then
  write(logfileunit,'("grid: nx[ce], ny[ce], nz[ce] =",6(1x,i5))') nxc, nxe, nyc, nye, nzc, nze
  write(     *     ,'("grid: nx[ce], ny[ce], nz[ce] =",6(1x,i5))') nxc, nxe, nyc, nye, nzc, nze
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
integer :: XEVarID, XCVarID
integer :: YEVarID, YCVarID
integer :: ZEVarID, ZCVarID

! for the prognostic variables
integer :: ivar, VarID

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
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'ncommas' ), &
           'nc_write_model_atts', 'model put '//trim(filename))

!-------------------------------------------------------------------------------
! Determine shape of most important namelist
!-------------------------------------------------------------------------------

call find_textfile_dims('ncommas_vars.nml', nlines, linelen)
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

   ! X-Grid Centers
   call nc_check(nf90_def_var(ncFileID,name='XC',xtype=nf90_real, &
                       dimids=NxcDimID,varid=XCVarID), &
                     'nc_write_model_atts', 'XC def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, XCVarID, 'type', 'x1d'),  &
                 'nc_write_model_atts','XC type '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, XCVarID, 'long_name', 'SCALAR GRID POSITION IN X'), &
                 'nc_write_model_atts','XC long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, XCVarID, 'units', 'meters'),  &
                 'nc_write_model_atts','XC units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, XCVarID, 'cartesian_axis', 'X'),   &
                 'nc_write_model_atts','XC cartesian_axis '//trim(filename))

   ! X-Grid Edges
   call nc_check(nf90_def_var(ncFileID,name='XE', xtype=nf90_real, &
                     dimids=NxeDimID, varid= XEVarID), &
                     'nc_write_model_atts', 'XE def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, XEVarID, 'type', 'x1d'),  &
                 'nc_write_model_atts','XE type '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, XEVarID, 'long_name', 'STAGGERED GRID POSITION IN X'), &
                 'nc_write_model_atts','XE long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, XEVarID, 'units', 'meters'),  &
                 'nc_write_model_atts','XE units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, XEVarID, 'cartesian_axis', 'X'),   &
                 'nc_write_model_atts','XE cartesian_axis '//trim(filename))

   ! Y-Grid Centers
   call nc_check(nf90_def_var(ncFileID,name='YC', xtype=nf90_real, &
                     dimids=NycDimID, varid= YCVarID), &
                     'nc_write_model_atts', 'YC def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YCVarID, 'type', 'y1d'),  &
                 'nc_write_model_atts','YC type '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YCVarID, 'long_name', 'SCALAR GRID POSITION IN Y'), &
                 'nc_write_model_atts','YC long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YCVarID, 'units', 'meters'),  &
                 'nc_write_model_atts','YC units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YCVarID, 'cartesian_axis', 'Y'),   &
                 'nc_write_model_atts','YC cartesian_axis '//trim(filename))

   ! Y-Grid Edges
   call nc_check(nf90_def_var(ncFileID,name='YE',xtype=nf90_real, &
                       dimids=NyeDimID,varid=YEVarID), &
                     'nc_write_model_atts', 'YE def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YEVarID, 'type', 'y1d'),  &
                 'nc_write_model_atts','YE type '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YEVarID, 'long_name', 'STAGGERED GRID POSITION IN Y'), &
                 'nc_write_model_atts','YE long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YEVarID, 'units', 'meters'),  &
                 'nc_write_model_atts','YE units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YEVarID, 'cartesian_axis', 'Y'),   &
                 'nc_write_model_atts','YE cartesian_axis '//trim(filename))

   ! Z-Grid Centers
   call nc_check(nf90_def_var(ncFileID,name='ZC',xtype=nf90_real, &
                       dimids=NzcDimID,varid=ZCVarID), &
                     'nc_write_model_atts', 'ZC def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, 'type', 'z1d'),  &
                 'nc_write_model_atts','ZC type '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, 'long_name', 'SCALAR GRID POSITION IN Z'), &
                 'nc_write_model_atts','ZC long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, 'units', 'meters'),  &
                 'nc_write_model_atts','ZC units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, 'positive', 'up'),  &
                 'nc_write_model_atts','ZC units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, 'cartesian_axis', 'Z'),   &
                 'nc_write_model_atts','ZC cartesian_axis '//trim(filename))

   ! Z-Grid Edges
   call nc_check(nf90_def_var(ncFileID,name='ZE', xtype=nf90_real, &
                     dimids=NzeDimID, varid= ZEVarID), &
                     'nc_write_model_atts', 'ZE def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZEVarID, 'type', 'z1d'),  &
                 'nc_write_model_atts','ZE type '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZEVarID, 'long_name', 'STAGGERED GRID POSITION IN Z'), &
                 'nc_write_model_atts','ZE long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZEVarID, 'units', 'meters'),  &
                 'nc_write_model_atts','ZE units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZEVarID, 'positive', 'up'),  &
                 'nc_write_model_atts','ZE units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZEVarID, 'cartesian_axis', 'Z'),   &
                 'nc_write_model_atts','ZE cartesian_axis '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------

   do ivar=1, nfields

      varname = trim(progvar(ivar)%varname)
      string1 = trim(filename)//' '//trim(varname)

      ! match shape of the variable to the dimension IDs

      call define_var_dims(progvar(ivar), myndims, mydimids, MemberDimID, unlimitedDimID, &
                      NxcDimID, NycDimID, NzcDimID, NxeDimID, NyeDimID, NzeDimID, & 
                      nxc     , nyc     , nzc     , nxe     , nye     , nze      )

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

   call nc_check(nf90_put_var(ncFileID, XCVarID, XC ), &
                'nc_write_model_atts', 'XC put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, XEVarID, XE ), &
                'nc_write_model_atts', 'XE put_var '//trim(filename))

   call nc_check(nf90_put_var(ncFileID, YCVarID, YC ), &
                'nc_write_model_atts', 'YC put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, YEVarID, YE ), &
                'nc_write_model_atts', 'YE put_var '//trim(filename))

   call nc_check(nf90_put_var(ncFileID, ZCVarID, ZC ), &
                'nc_write_model_atts', 'ZC put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, ZEVarID, ZE ), &
                'nc_write_model_atts', 'ZE put_var '//trim(filename))

endif

!-------------------------------------------------------------------------------
! Fill the variables we can
!-------------------------------------------------------------------------------

if (has_ncommas_namelist) then
   call file_to_text('ncommas_vars.nml', textblock)
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

   do ivar = 1,nfields  ! Very similar to loop in sv_to_restart_file

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



subroutine get_gridsize(num_xc, num_xe, num_yc, num_ye, num_zc, num_ze )
 integer, intent(out) :: num_xc, num_yc, num_zc
 integer, intent(out) :: num_xe, num_ye, num_ze
!------------------------------------------------------------------
! public utility routine.

if ( .not. module_initialized ) call static_init_model

 num_xc = nxc
 num_yc = nyc
 num_zc = nzc

 num_xe = nxe
 num_ye = nye
 num_ze = nze

end subroutine get_gridsize



subroutine restart_file_to_sv(filename, state_vector, model_time)
!------------------------------------------------------------------
! Reads the current time and state variables from a ncommas restart
! file and packs them into a dart state vector.

character(len=*), intent(in)    :: filename 
real(r8),         intent(inout) :: state_vector(:)
type(time_type),  intent(out)   :: model_time

! temp space to hold data while we are reading it
integer  :: i, j, k, l, ni, nj, nk, nl, ivar, indx
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array
real(r8), allocatable, dimension(:,:,:)     :: data_3d_array
real(r8), allocatable, dimension(:,:,:,:)   :: data_4d_array

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
   call error_handler(E_ERR,'restart_file_to_sv',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncid), &
             'restart_file_to_sv','open '//trim(filename))

model_time = get_state_time(ncid, filename)

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
            'restart_file_to_sv', 'inquire timedimlength '//trim(filename))
else
   TimeDimLength = 0
endif

do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   myerrorstring = trim(filename)//' '//trim(varname)

   ! determine the shape of the netCDF variable 

   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            'restart_file_to_sv', 'inq_varid '//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid,VarId,dimids=dimIDs,ndims=ncNdims), &
            'restart_file_to_sv', 'inquire '//trim(myerrorstring))

   mystart = 1   ! These are arrays, actually.
   mycount = 1
   ! Only checking the shape of the variable - sans TIME
   DimCheck : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'restart_file_to_sv', string1)

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(myerrorstring),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         call error_handler(E_ERR,'restart_file_to_sv',string1,source,revision,revdate)
      endif

      mycount(i) = dimlen

   enddo DimCheck

   where(dimIDs == TimeDimID) mystart = TimeDimLength
   where(dimIDs == TimeDimID) mycount = 1   ! only the latest one

   if ( debug > 1 ) then
      write(*,*)'restart_file_to_sv '//trim(varname)//' start = ',mystart(1:ncNdims)
      write(*,*)'restart_file_to_sv '//trim(varname)//' count = ',mycount(1:ncNdims)
   endif

   indx = progvar(ivar)%index1

   if (ncNdims == 1) then

      ! If the single dimension is TIME, we only need a scalar.
      ! Pretty sure this cant happen given the test for x1d,y1d,z1d. 
      ni = mycount(1)
      allocate(data_1d_array(ni))
      call nc_check(nf90_get_var(ncid, VarID, data_1d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'restart_file_to_sv', 'get_var '//trim(varname))
      do i = 1, ni
         state_vector(indx) = data_1d_array(i)
         indx = indx + 1
      enddo
      deallocate(data_1d_array)

   elseif (ncNdims == 2) then

      ni = mycount(1)
      nj = mycount(2)
      allocate(data_2d_array(ni, nj))
      call nc_check(nf90_get_var(ncid, VarID, data_2d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'restart_file_to_sv', 'get_var '//trim(varname))
      do j = 1, nj
      do i = 1, ni
         state_vector(indx) = data_2d_array(i, j)
         indx = indx + 1
      enddo
      enddo
      deallocate(data_2d_array)

   elseif (ncNdims == 3) then

      ni = mycount(1)
      nj = mycount(2)
      nk = mycount(3)
      allocate(data_3d_array(ni, nj, nk))
      call nc_check(nf90_get_var(ncid, VarID, data_3d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'restart_file_to_sv', 'get_var '//trim(varname))
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
         state_vector(indx) = data_3d_array(i, j, k)
         indx = indx + 1
      enddo
      enddo
      enddo
      deallocate(data_3d_array)

   elseif (ncNdims == 4) then

      ni = mycount(1)
      nj = mycount(2)
      nk = mycount(3)
      nl = mycount(4)
      allocate(data_4d_array(ni, nj, nk, nl))
      call nc_check(nf90_get_var(ncid, VarID, data_4d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'restart_file_to_sv', 'get_var '//trim(varname))
      do l = 1, nl
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
         state_vector(indx) = data_4d_array(i, j, k, l)
         indx = indx + 1
      enddo
      enddo
      enddo
      enddo
      deallocate(data_4d_array)

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

call nc_check(nf90_close(ncid), &
             'restart_file_to_sv','close '//trim(filename))

end subroutine restart_file_to_sv



subroutine sv_to_restart_file(state_vector, filename, statedate)
!------------------------------------------------------------------
! Writes the current time and state variables from a dart state
! vector (1d array) into a ncommas netcdf restart file.
!
real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: filename 
type(time_type),  intent(in) :: statedate

! temp space to hold data while we are writing it
integer :: i, ni, nj, nk, nl, ivar
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array
real(r8), allocatable, dimension(:,:,:)     :: data_3d_array
real(r8), allocatable, dimension(:,:,:,:)   :: data_4d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME) :: varname 
integer :: VarID, ncNdims, dimlen
integer :: ncFileID, TimeDimID, TimeDimLength

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
! of the ncommas restart file, and state vector contents from a different
! time won't be consistent with the rest of the file.

model_time = get_state_time(ncFileID, filename)

if ( model_time /= statedate ) then
   call print_time( statedate,'DART current time',logfileunit) 
   call print_time(model_time,'ncommas  current time',logfileunit) 
   call print_time( statedate,'DART current time') 
   call print_time(model_time,'ncommas  current time') 
   write(string1,*)trim(filename),' current time /= model time. FATAL error.'
   call error_handler(E_ERR,'sv_to_restart_file',string1,source,revision,revdate) 
endif

if (do_output()) &
    call print_time(statedate,'time of restart file '//trim(filename))
if (do_output()) &
    call print_date(statedate,'date of restart file '//trim(filename))

! The DART prognostic variables are only defined for a single time.
! We already checked the assumption that variables are xy2d or xyz3d ...
! IF the netCDF variable has a TIME dimension, it must be the last dimension,
! and we need to read the LAST timestep and effectively squeeze out the
! singleton dimension when we stuff it into the DART state vector. 

TimeDimID = FindTimeDimension( ncFileID )

if ( TimeDimID > 0 ) then
   call nc_check(nf90_inquire_dimension(ncFileID, TimeDimID, len=TimeDimLength), &
            'sv_to_restart_file', 'inquire timedimlength '//trim(filename))
else
   TimeDimLength = 0
endif

do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   string2 = trim(filename)//' '//trim(varname)

   ! Ensure netCDF variable is conformable with progvar quantity.
   ! The TIME and Copy dimensions are intentionally not queried
   ! by looping over the dimensions stored in the progvar type.

   call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'sv_to_restart_file', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncFileID,VarId,dimids=dimIDs,ndims=ncNdims), &
            'sv_to_restart_file', 'inquire '//trim(string2))

   mystart = 1   ! These are arrays, actually.
   mycount = 1
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

      mycount(i) = dimlen

   enddo DimCheck

   where(dimIDs == TimeDimID) mystart = TimeDimLength
   where(dimIDs == TimeDimID) mycount = 1   ! only the latest one

   if ( debug > 1 ) then
      write(*,*)'sv_to_restart_file '//trim(varname)//' start is ',mystart(1:ncNdims)
      write(*,*)'sv_to_restart_file '//trim(varname)//' count is ',mycount(1:ncNdims)
   endif

   if (progvar(ivar)%numdims == 1) then
      ni = mycount(1)
      allocate(data_1d_array(ni))
      call vector_to_prog_var(state_vector, progvar(ivar), data_1d_array)
      call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array, &
            start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'sv_to_restart_file', 'put_var '//trim(varname))
      deallocate(data_1d_array)

   elseif (progvar(ivar)%numdims == 2) then

      ni = mycount(1)
      nj = mycount(2)
      allocate(data_2d_array(ni, nj))
      call vector_to_prog_var(state_vector, progvar(ivar), data_2d_array)
      
      if ( progvar(ivar)%posdef == 1 ) then
        where ( data_2d_array < 0 ) data_2d_array = 0
      endif

      call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'sv_to_restart_file', 'put_var '//trim(varname))
      deallocate(data_2d_array)

   elseif (progvar(ivar)%numdims == 3) then

      ni = mycount(1)
      nj = mycount(2)
      nk = mycount(3)
      allocate(data_3d_array(ni, nj, nk))
      call vector_to_prog_var(state_vector, progvar(ivar), data_3d_array)

      if ( progvar(ivar)%posdef == 1 ) then
        where ( data_3d_array < 0 ) data_3d_array = 0
      endif

      call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'sv_to_restart_file', 'put_var '//trim(varname))
      deallocate(data_3d_array)

   elseif (progvar(ivar)%numdims == 4) then

      ni = mycount(1)
      nj = mycount(2)
      nk = mycount(3)
      nl = mycount(4)
      allocate(data_4d_array(ni, nj, nk, nl))
      call vector_to_prog_var(state_vector, progvar(ivar), data_4d_array)

      call nc_check(nf90_put_var(ncFileID, VarID, data_4d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'sv_to_restart_file', 'put_var '//trim(varname))
      deallocate(data_4d_array)

   else
      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'sv_to_restart_file', string1, &
                        source,revision,revdate)
   endif

enddo

call nc_check(nf90_close(ncFileID), &
             'sv_to_restart_file','close '//trim(filename))

end subroutine sv_to_restart_file



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
!     to that location for the filter from the NCOMMAS state vectors
!
!############################################################################
!
!     Author:  Tim Hoar, Ted Mansell, Lou Wicker
!
!     Creation Date:  August 2010
!     
!     Variables needed to be stored in the MODEL_MODULE data structure
!
!       XG_POS   = time dependent X (meters) offset of model from its lat/lon point
!       YG_POS   = time dependent Y (meters) offset of model from its lat/lon point
!       XE, XC   = 1D arrays storing the local grid edge and center coords (meters)
!       YE, YC   = 1D arrays storing the local grid edge and center coords (meters)
!       ZE, ZC   = 1D arrays storing the local grid edge and center coords (meters)
!       REF_LAT  = grid reference latitude
!       REF_LON  = grid reference longitude
!
!       ERROR codes:
!
!       ISTATUS = 99:  general error in case something terrible goes wrong...
!       ISTATUS = 15:  dont know what to do with vertical coord supplied
!       ISTATUS = 16:  Obs_type is not found
!       ISTATUS = 11:  index from from xi is outside domain WRT XE
!       ISTATUS = 12:  index from from xi is outside domain WRT XC
!       ISTATUS = 21:  index from from yi is outside domain WRT YE
!       ISTATUS = 22:  index from from yi is outside domain WRT YC
!       ISTATUS = 31:  index from from zi is outside domain WRT ZE
!       ISTATUS = 32:  index from from zi is outside domain WRT ZC
!       ISTATUS = 33:  index from location(3) is not bounded by 1 --> nz-1
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

  real(r8)         :: loc_array(3), llon, llat
  real(r8)         :: lheight
  integer          :: base_offset
  integer          :: nf, i, n0, n1
  integer          :: iloc, jloc, kloc, nxp, nyp, nzp
  real(r8)         :: xi, yi, xf, yf, zf, q1, q2, vt, vb

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

  IF( vert_is_height(location) ) THEN        ! Nothing to do 
  
  ELSEIF ( vert_is_surface(location) ) THEN  ! Nothing to do
 
  ELSEIF (vert_is_level(location)) THEN      ! convert the level index to an actual height
                                             ! This code is wrong (I think) if variable == "W"
     kloc = nint(loc_array(3))
     IF( (kloc < 1) .or. (kloc > size(zc)) ) THEN 
        istatus = 33
        return
     ELSE
        lheight = zc(kloc)
     ENDIF
  ELSE   ! if we don't know what to do
     istatus = 15
     return
  ENDIF

! Find what field this is....

  nf = -1
  DO i = 1,nfields
     IF (progvar(i)%dart_kind == obs_type) THEN
       nf = i
       exit
     ENDIF
  ENDDO

! Not a legal type for interpolation, return istatus error

  IF (nf < 0) then
     istatus = 16
     return
  ENDIF

  base_offset = progvar(nf)%index1

  IF (debug > 2) print *, 'base offset now ', base_offset

! Not implemented...

  IF( vert_is_surface(location) ) THEN
     istatus = 15
     return
  ENDIF
  
! From input lat/lon position, find x/y distance from grid lat/lon
! Assumed lat/lon input is in degrees (.true. flag is used to convert inside routine)

  CALL ll_to_xy(xi, yi, 0, ref_lat, ref_lon, llat, llon, .true.)
  
  IF( debug > 2 ) THEN
    print *, 'Ref_LON / OBS_LON / XMIN / X-LOC from lat/lon / XMAX:  ', ref_lon, llon, xe(1), xi, xe(nxe)
    print *, 'Ref_LAT / OBS_LAT / YMIN / Y-LOC from lat/lon / YMAX:  ', ref_lat, llat, ye(1), yi, ye(nye)
  ENDIF
  
! Using nf, get the true dimensions of the variable

   nxp = progvar(nf)%dimlens(1)
   nyp = progvar(nf)%dimlens(2)
   nzp = progvar(nf)%dimlens(3)   
   
! Get X (LON) index

  IF(obs_type == 1) THEN
    iloc = find_index(xi, xe, nxe)
    IF( iloc == -1) THEN  ! Error occurred, observation is outside of the domain
     istatus = 11
     return
    ENDIF
    xf = (xe(iloc+1) - xi) / (xe(iloc+1) - xe(iloc))
  ELSE
    iloc = find_index(xi, xc, nxc)
    IF( iloc == -1) THEN  ! Error occurred, observation is outside of the domain
     istatus = 12
     return
    ENDIF
    xf = (xc(iloc+1) - xi) / (xc(iloc+1) - xc(iloc))
  ENDIF
  
! Get Y (LAT) index

   IF(obs_type == 2) THEN
    jloc = find_index(yi, ye, nye)
    IF( jloc == -1) THEN  ! Error occurred, observation is outside of the vertical domain
     istatus = 21
     return
    ENDIF
    yf = (ye(jloc+1) - yi) / (ye(jloc+1) - ye(jloc))
  ELSE
    jloc = find_index(yi, yc, nyc) 
    IF( jloc == -1) THEN  ! Error occurred, observation is outside of the vertical domain
     istatus = 22
     return
    ENDIF
    yf = (yc(jloc+1) - yi) / (yc(jloc+1) - yc(jloc))
  ENDIF
  
! Get Z (HGT) index

  IF(obs_type == 3) THEN
    kloc = find_index(lheight, ze, nze)
    IF( kloc == -1) THEN  ! Error occurred, observation is outside of the vertical domain
     istatus = 31
     return
    ENDIF
    zf = (ze(kloc+1) - lheight) / (ze(kloc+1) - ze(kloc))
  ELSE
    kloc = find_index(lheight, zc, nzc)
    IF( kloc == -1) THEN  ! Error occurred, observation is outside of the vertical domain
     istatus = 32
     return
    ENDIF
    zf = (zc(kloc+1) - lheight) / (zc(kloc+1) - zc(kloc))
  ENDIF

  IF( debug > 2 ) THEN
    print *, "ILOC:  ", iloc
    print *, "XI = ", xi
    print *, "XF = ", xf
    print *, "JLOC:  ", jloc
    print *, "YI = ", yi 
    print *, "YF = ", yf
    print *, "ZLOC = ", kloc
    print *, "ZI = ", lheight 
    print *, "ZF = ", zf
  ENDIF
  
! Find all the corners of the trilinear cube and then form the 1D->2D->3D interpolation
! Since base_offset is always the first point (e.g., (1,1,1)), then subtract "1" to get that loc
  
  n0 = base_offset + (jloc-1)*nxp + (kloc-1)*nxp*nyp + iloc - 1  ! Location of (i,  j,k) ?                         
  n1 = base_offset + (jloc-1)*nxp + (kloc-1)*nxp*nyp + iloc      ! Location of (i+1,j,k) ?                    
  q1 = (1.0_r8-xf)*x(n1) + xf*x(n0)                              ! Value along "j,k" grid line
  
  n0 = base_offset + (jloc)*nxp + (kloc-1)*nxp*nyp + iloc - 1    ! Location of (i,  j+1,k) ?                         
  n1 = base_offset + (jloc)*nxp + (kloc-1)*nxp*nyp + iloc        ! Location of (i+1,j+1,k) ?                    
  q2 = (1.0_r8-xf)*x(n1) + xf*x(n0)                              ! Value along "j+1,k" grid line
  
  vb = (1.0_r8-yf)*q2 + yf*q1                                    ! Binlinear value on the bottom plane
  
  n0 = base_offset + (jloc-1)*nxp + (kloc)*nxp*nyp + iloc - 1    ! Location of (i,  j,k+1) ?                         
  n1 = base_offset + (jloc-1)*nxp + (kloc)*nxp*nyp + iloc        ! Location of (i+1,j,k+1) ?                    
  q1 = (1.0_r8-xf)*x(n1) + xf*x(n0)                              ! Value along "j,k+1" grid line
  
  n0 = base_offset + (jloc)*nxp + (kloc)*nxp*nyp + iloc - 1      ! Location of (i,  j+1,k+1) ?                         
  n1 = base_offset + (jloc)*nxp + (kloc)*nxp*nyp + iloc          ! Location of (i+1,j+1,k+1) ?                    
  q2 = (1.0_r8-xf)*x(n1) + xf*x(n0)                              ! Value along "j+1,k+1" grid line
  
  vt = (1.0_r8-yf)*q2 + yf*q1                                    ! Binlinear value on the bottom plane
  
  interp_val = (1.0_r8-zf)*vt + zf*vb                            ! Trilinear value
  
! All good.
  istatus = 0

return
end subroutine model_interpolate


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


subroutine get_grid_dims(NXC, NXE, NYC, NYE, NZC, NZE)
!------------------------------------------------------------------
!
! Read the grid dimensions from the restart netcdf file.
!
! The file name comes from module storage ... namelist.

integer, intent(out) :: NXC   ! Number of Longitude centers
integer, intent(out) :: NXE   ! Number of Longitude edges
integer, intent(out) :: NYC   ! Number of Latitude  centers
integer, intent(out) :: NYE   ! Number of Latitude  edges
integer, intent(out) :: NZC   ! Number of Vertical grid centers
integer, intent(out) :: NZE   ! Number of Vertical grid edges

integer :: grid_id, dimid

if ( .not. module_initialized ) call static_init_model

! get the ball rolling ...

call nc_check(nf90_open(trim(ncommas_restart_filename), nf90_nowrite, grid_id), &
            'get_grid_dims','open '//trim(ncommas_restart_filename))

! Longitudes : get dimid for 'XC' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'XC', dimid), &
            'get_grid_dims','inq_dimid XC '//trim(ncommas_restart_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=NXC), &
            'get_grid_dims','inquire_dimension XC '//trim(ncommas_restart_filename))

! Longitudes : get dimid for 'XE and then get value

call nc_check(nf90_inq_dimid(grid_id, 'XE', dimid), &
            'get_grid_dims','inq_dimid XE '//trim(ncommas_restart_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=NXE), &
            'get_grid_dims','inquire_dimension XE '//trim(ncommas_restart_filename))

! Latitudes : get dimid for 'YC' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'YC', dimid), &
            'get_grid_dims','inq_dimid YC '//trim(ncommas_restart_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=NYC), &
            'get_grid_dims','inquire_dimension YC '//trim(ncommas_restart_filename))

! Latitudes : get dimid for 'YE' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'YE', dimid), &
            'get_grid_dims','inq_dimid YE '//trim(ncommas_restart_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=NYE), &
            'get_grid_dims','inquire_dimension YE '//trim(ncommas_restart_filename))

! Vertical Levels : get dimid for 'ZC' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'ZC', dimid), &
            'get_grid_dims','inq_dimid ZC '//trim(ncommas_restart_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=NZC), &
            'get_grid_dims','inquire_dimension ZC '//trim(ncommas_restart_filename))

! Vertical Levels : get dimid for 'ZE' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'ZE', dimid), &
            'get_grid_dims','inq_dimid ZE '//trim(ncommas_restart_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=NZE), &
            'get_grid_dims','inquire_dimension ZE '//trim(ncommas_restart_filename))

! tidy up

call nc_check(nf90_close(grid_id), &
         'get_grid_dims','close '//trim(ncommas_restart_filename) )

end subroutine get_grid_dims


!------------------------------------------------------------------


subroutine get_grid(nxc, nxe, nyc, nye, nzc, nze,       &
                    xc, xe, yc, ye, zc, ze,             &
                    ulat, ulon, vlat, vlon, wlat, wlon, &
                    xg_pos, yg_pos, ref_lat, ref_lon, hgt_offset)
!------------------------------------------------------------------
!
! Read the grid dimensions from the restart netcdf file.
!
! The file name comes from module storage ... namelist.

integer, intent(in) :: NXC   ! Number of Longitude centers
integer, intent(in) :: NXE   ! Number of Longitude edges
integer, intent(in) :: NYC   ! Number of Latitude  centers
integer, intent(in) :: NYE   ! Number of Latitude  edges
integer, intent(in) :: NZC   ! Number of Vertical grid centers
integer, intent(in) :: NZE   ! Number of Vertical grid edges

real(r8), dimension(:,:), intent(out) :: ULAT, ULON, VLAT, VLON, WLAT, WLON
real(r8), dimension( : ), intent(out) :: XC, XE, YC, YE, ZC, ZE

real(r8),                 intent(out) :: ref_lat, ref_lon
real(r8),                 intent(out) :: xg_pos, yg_pos, hgt_offset

integer  :: ncid, VarID
integer  :: i,j
real(r8) :: x,y,lat,lon

! Read the netcdf file data

call nc_check(nf90_open(trim(ncommas_restart_filename), nf90_nowrite, ncid), 'get_grid', 'open '//trim(ncommas_restart_filename))

! fixme - in a perfect world - 
! Get the variable ID
! Check to make sure it is the right shape
! Read it

call nc_check(nf90_inq_varid(ncid, 'XC', VarID), 'get_grid', 'inq_varid XC '//trim(ncommas_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   XC), 'get_grid',   'get_var XC '//trim(ncommas_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'XE', VarID), 'get_grid', 'inq_varid XE '//trim(ncommas_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   XE), 'get_grid',   'get_var XE '//trim(ncommas_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'YC', VarID), 'get_grid', 'inq_varid YC '//trim(ncommas_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   YC), 'get_grid',   'get_var YC '//trim(ncommas_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'YE', VarID), 'get_grid', 'inq_varid YE '//trim(ncommas_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   YE), 'get_grid',   'get_var YE '//trim(ncommas_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'ZC', VarID), 'get_grid', 'inq_varid ZC '//trim(ncommas_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   ZC), 'get_grid',   'get_var ZC '//trim(ncommas_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'ZE', VarID), 'get_grid', 'inq_varid ZE '//trim(ncommas_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   ZE), 'get_grid',   'get_var ZE '//trim(ncommas_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'LAT',   VarID), 'get_grid', 'inq_varid LAT '//trim(ncommas_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID, ref_lat), 'get_grid',   'get_var LAT '//trim(ncommas_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'LON',   VarID), 'get_grid', 'inq_varid LON '//trim(ncommas_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID, ref_lon), 'get_grid',   'get_var LON '//trim(ncommas_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'XG_POS', VarID), 'get_grid', 'inq_varid XG_POS '//trim(ncommas_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   xg_pos), 'get_grid',   'get_var XG_POS '//trim(ncommas_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'YG_POS', VarID), 'get_grid', 'inq_varid YG_POS '//trim(ncommas_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   yg_pos), 'get_grid',   'get_var YG_POS '//trim(ncommas_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'HGT',      VarID), 'get_grid', 'inq_varid HGT '//trim(ncommas_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID, hgt_offset), 'get_grid',   'get_var HGT '//trim(ncommas_restart_filename))

call nc_check(nf90_close(ncid), 'get_grid','close '//trim(ncommas_restart_filename) )

! A little sanity check

if ( debug > 1 ) then

   write(*,*)'xc range ',minval(xc),maxval(xc)
   write(*,*)'xe range ',minval(xe),maxval(xe)
   write(*,*)'yc range ',minval(yc),maxval(yc)
   write(*,*)'ye range ',minval(ye),maxval(ye)
   write(*,*)'zc range ',minval(zc),maxval(zc)
   write(*,*)'ze range ',minval(ze),maxval(ze)
   write(*,*)'ref_lat/LAT    ',ref_lat
   write(*,*)'ref_lon/LON    ',ref_lon
   write(*,*)'xg_pos/XG_POS  ',xg_pos
   write(*,*)'yg_pos/YG_POS  ',yg_pos
   write(*,*)'hgt_offset/HGT ',hgt_offset

endif

! Here convert model vertical height to height above sea level

zc(:) = zc(:) + hgt_offset
ze(:) = ze(:) + hgt_offset

! Do the same with the horizontal grid offsets so that all the coordinates are now relative to lat/lon reference point

xc(:) = xc(:) + xg_pos
xe(:) = xe(:) + xg_pos
yc(:) = yc(:) + yg_pos
ye(:) = ye(:) + yg_pos
  
! Create lat lons

DO j = 1,nyc
  DO i = 1,nxc
    call xy_to_ll(wlat(i,j), wlon(i,j), 0, xc(i), yc(j), ref_lat, ref_lon, .true.)
  ENDDO  
ENDDO  
      
DO j = 1,nyc
  DO i = 1,nxe
    call xy_to_ll(ulat(i,j), ulon(i,j), 0, xe(i), yc(j), ref_lat, ref_lon, .true.)
  ENDDO  
ENDDO  

DO j = 1,nye
  DO i = 1,nxc
    call xy_to_ll(vlat(i,j), vlon(i,j), 0, xc(i), ye(j), ref_lat, ref_lon, .true.)
  ENDDO  
ENDDO  

! check:

IF( debug > 2 ) THEN
  DO j = 1,nyc,4
    DO i = 1,nxc,4
      call ll_to_xy(x, y, 0, ref_lat, ref_lon, wlat(i,j), wlon(i,j),.true.)
      call xy_to_ll(lat, lon, 0, xc(i), yc(j), ref_lat, ref_lon)
      write(*,*) 'i,j,x,y,x1,y1: ',i,j,xc(i), yc(j),x,y,wlat(i,j),wlon(i,j),lat,lon
    ENDDO  
  ENDDO 
ENDIF

! IF WE DO THIS, THEN CAN ASSUME THAT LON > 0 IN LL_TO_XY
where (ULON <   0.0_r8) ULON = ULON + 360.0_r8
where (ULON > 360.0_r8) ULON = ULON - 360.0_r8
where (VLON <   0.0_r8) VLON = VLON + 360.0_r8
where (VLON > 360.0_r8) VLON = VLON - 360.0_r8
where (WLON <   0.0_r8) WLON = WLON + 360.0_r8
where (WLON > 360.0_r8) WLON = WLON - 360.0_r8

where (ULAT < -90.0_r8) ULAT = -90.0_r8
where (ULAT >  90.0_r8) ULAT =  90.0_r8
where (VLAT < -90.0_r8) VLAT = -90.0_r8
where (VLAT >  90.0_r8) VLAT =  90.0_r8
where (WLAT < -90.0_r8) WLAT = -90.0_r8
where (WLAT >  90.0_r8) WLAT =  90.0_r8

! Print a little summary.
if ( debug > 1 ) then
   write(*,*)'ulon longitude range ',minval(ulon),maxval(ulon)
   write(*,*)'ulat latitude  range ',minval(ulat),maxval(ulat)

   write(*,*)'vlon longitude range ',minval(vlon),maxval(vlon)
   write(*,*)'vlat latitude  range ',minval(vlat),maxval(vlat)

   write(*,*)'wlon longitude range ',minval(wlon),maxval(wlon)
   write(*,*)'wlat latitude  range ',minval(wlat),maxval(wlat)
endif

return
end subroutine get_grid


!------------------------------------------------------------------


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


function get_state_time_ncid( ncid, filename )
!------------------------------------------------------------------
! the static_init_model ensures that the ncommas namelists are read.
!
type(time_type)              :: get_state_time_ncid
integer,          intent(in) :: ncid
character(len=*), intent(in) :: filename

integer         :: VarID, numdims, dimlen
type(time_type) :: model_offset, base_time

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
integer, allocatable, dimension(:)    :: mytimes

if ( .not. module_initialized ) call static_init_model

base_time = get_base_time(ncid)

call nc_check( nf90_inq_varid(ncid, 'TIME', VarID), &
                  'get_state_time', 'inq_varid TIME '//trim(filename))

call nc_check(nf90_inquire_variable(ncid, VarId, dimids=dimIDs, ndims=numdims), &
                  'get_state_time', 'inquire TIME '//trim(filename))

if ( numdims > 1 ) then
   write(string1,*) 'TIME is not expected to have ',numdims,' dimensions.'
   call error_handler(E_ERR,'get_state_time',string1,source,revision,revdate)
endif

call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlen), &
            'get_state_time', 'inquire time dimension length '//trim(filename))

allocate(mytimes(dimlen))

call nc_check( nf90_get_var(ncid, VarID, mytimes ), &
                  'get_state_time', 'get_var TIME '//trim(filename))

model_offset = set_time(maxval(mytimes))

get_state_time_ncid = base_time + model_offset

deallocate(mytimes)

end function get_state_time_ncid


!------------------------------------------------------------------


function get_state_time_fname(filename)
!------------------------------------------------------------------
! the static_init_model ensures that the ncommas namelists are read.
!
type(time_type) :: get_state_time_fname
character(len=*), intent(in) :: filename

integer         :: ncid, VarID, numdims, dimlen
type(time_type) :: model_offset, base_time

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
integer, allocatable, dimension(:)    :: mytimes

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'get_state_time',string1,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'get_base_time', 'open '//trim(filename))

base_time = get_base_time(ncid)

call nc_check( nf90_inq_varid(ncid, 'TIME', VarID), &
                  'get_state_time', 'inq_varid TIME '//trim(filename))

call nc_check(nf90_inquire_variable(ncid, VarId, dimids=dimIDs, ndims=numdims), &
                  'get_state_time', 'inquire TIME '//trim(filename))

if ( numdims > 1 ) then
   write(string1,*) 'TIME is not expected to have ',numdims,' dimensions.'
   call error_handler(E_ERR,'get_state_time',string1,source,revision,revdate)
endif

call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlen), &
            'get_state_time', 'inquire time dimension length '//trim(filename))

allocate(mytimes(dimlen))

call nc_check( nf90_get_var(ncid, VarID, mytimes ), &
                  'get_state_time', 'get_var TIME '//trim(filename))
call nc_check(nf90_close(ncid), 'get_state_time', 'close '//trim(filename))

write(*,*)' temporal offset is (in seconds) is ',maxval(mytimes)
model_offset = set_time(maxval(mytimes))

get_state_time_fname = base_time + model_offset

deallocate(mytimes)

end function get_state_time_fname


!------------------------------------------------------------------


function set_model_time_step()
!------------------------------------------------------------------
! the static_init_model ensures that the ncommas namelists are read.
!
type(time_type) :: set_model_time_step

if ( .not. module_initialized ) call static_init_model

! FIXME - determine when we can stop the model

   set_model_time_step = set_time(0, 1) ! (seconds, days)

end function set_model_time_step


!------------------------------------------------------------------


subroutine get_ncommas_restart_filename( filename )

character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(ncommas_restart_filename)

end subroutine get_ncommas_restart_filename


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
      string1 = 'ncommas_vars_nml:ncommas_state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Make sure variable exists in netCDF file

   write(string1,'(''there is no variable '',a,'' in '',a)') trim(varname), trim(filename)
   call nc_check(NF90_inq_varid(ncid, trim(varname), varid), &
                 'verify_state_variables', trim(string1))

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


!###########################################################################
!
!     ##################################################################
!     ######                                                      ######
!     ######             REAL FUNCTION LINTERP_1D                 ######
!     ######                                                      ######
!     ##################################################################
!
!     PURPOSE:
!
!############################################################################
!
!     Author:  Lou Wicker for DART interp
!
!############################################################################

function linterp1D(x, x0, x1, y0, y1)

  real(r8)             :: linterp1D
  real(r8), intent(in) :: x, x0, x1, y0, y1
  
  real(r8) :: f
  
  f = (x1 - x) / (x1 - x0)
  
  IF( f < 0.0_r8 ) THEN
    write(*,*) "LINTERP ERROR, X outside X0->X1:  x/x0/x1", x, x0, x1
    linterp1d = -999.0_r8
    return
  ENDIF
  
  linterp1d = y0*f + (1.0_r8-f)*y1
  
return
end function linterp1D


!############################################################################
!
!     ##################################################################
!     ######                                                      ######
!     ######                   SUBROUTINE LL_TO_XY                ######
!     ######                                                      ######
!     ##################################################################
!
!
!     PURPOSE:
!
!     This subroutine computes the projected (x, y) coordinates of the
!     point (lat2, lon2) relative to (lat1, lon1).
!
!############################################################################
!
!     Author:  David Dowell
!
!     Creation Date:  25 February 2005
!     Modified:  August 2010:  Number of mods to add optional arg to convert 
!                from and to degrees for input/output, and some a check
!                to see if the inputs are in degrees (Lou Wicker)
!
!                Note, the automatic checks will FAIL with small latitudes
!                or longitudes - we cannot differntiate between deg and rad
!                
!############################################################################

  subroutine ll_to_xy(x, y, map_proj, lat1, lon1, lat2, lon2, degrees)

! Passed variables

    real(r8),     intent(out)     :: x, y           ! distance (m)  RETURNED FIELDS
    real(r8),     intent(in)      :: lat1, lon1     ! coordinates of first point (in deg or radians)
    real(r8),     intent(in)      :: lat2, lon2     ! coordinates of second point (in deg or radians)
    integer ,     intent(in)      :: map_proj       ! map projection:
                                                    !   0 = flat earth
                                                    !   1 = oblique azimuthal (not supported)
                                                    !   2 = Lambert conformal (not supported)

    logical, intent(in), optional :: degrees        ! If lat/lon inputs are in degrees, convert to radians

! Local variables

    real(r8)      :: llat1, llon1     ! local coordinates of first point (in radians)
    real(r8)      :: llat2, llon2     ! local coordinates of second point (in radians)
    real(r8)      :: llon1p, llon2p   ! local positive longitudes (in radians)

! In case they are already in radians..

    llat1 = lat1
    llat2 = lat2 
    llon1 = lon1 
    llon2 = lon2 

! If not, check things

    IF( present(degrees) ) THEN
      IF( degrees ) THEN
       llat1 = lat1 * deg2rad
       llat2 = lat2 * deg2rad
       llon1 = lon1 * deg2rad
       llon2 = lon2 * deg2rad
      ENDIF
     ELSE
! Try and be smart, in case user is stupid (like me!) - Note, this check will not work with lats ~ equator,
! nor lons near in western Europe, or mid-Pacific

      IF( abs(llat1) >  2.0_r8*PI ) llat1 = llat1 * deg2rad
      IF( abs(llat2) >  2.0_r8*PI ) llat2 = llat2 * deg2rad
      IF( abs(llon1) >  2.0_r8*PI ) llon1 = llon1 * deg2rad
      IF( abs(llon2) >  2.0_r8*PI ) llon2 = llon2 * deg2rad
    ENDIF

    if (llon1 < 0.0_r8) then
      llon1p = llon1+2.0_r8*PI
    else
      llon1p = llon1
    endif
    if (llon2 < 0.0_r8) then
      llon2p = llon2+2.0_r8*PI
    else
      llon2p = llon2
    endif

    if (map_proj == 0) then
      x = rearth * cos(0.5_r8*(llat1+llat2)) * (llon2p-llon1p)
      y = rearth * (llat2-llat1)
    else
       write(string1,*) 'Requested map projection unavailable: ', map_proj
       call error_handler(E_ERR,'ll_to_xy',string1,source,revision,revdate)
    endif

  return
  end subroutine ll_to_xy


!############################################################################
!
!     ##################################################################
!     ######                                                                                                                              ######
!     ######                  SUBROUTINE XY_TO_LL                 ######
!     ######                                                                                                                              ######
!     ##################################################################
!
!
!     PURPOSE:
!
!     This subroutine computes the projected (lat, lon) coordinates of the
!     point (x, y) relative to (lat0, lon0).  Various map projections
!     are possible.
!
!############################################################################
!
!     Author:  David Dowell
!
!     Creation Date:  25 February 2005
!     Modified:  August 2010:  Number of mods to add optional arg to convert 
!                from and to degrees for input/output, and some a check
!                to see if the inputs are in degrees - outputs are then
!                automatically converted back.  (Lou Wicker)
!
!                Note, the automatic checks will FAIL with small latitudes
!                or longitudes - we cannot differntiate between deg and rad
!############################################################################

  subroutine xy_to_ll(lat, lon, map_proj, x, y, lat1, lon1, degrees)
    
! Passed variables
  
    real(r8),     intent(in)     :: x, y         ! distance (in meters) of point
    real(r8),     intent(out)    :: lat, lon     ! coordinates of first point (in deg or radians)
    real(r8),     intent(in)     :: lat1, lon1   ! coordinates of second point (in deg or radians)
    integer ,     intent(in)     :: map_proj     ! map projection:
                                                 !   0 = flat earth
                                                 !   1 = oblique azimuthal (not supported)
                                                 !   2 = Lambert conformal (not supported)
                                
    logical, intent(in), optional :: degrees     ! If TRUE and lat/lon inputs are in degrees, return degrees

! Local variables

    real(r8)      :: llat1, llon1     ! local coordinates of reference lat/lon (in radians)

! In case they are already in radians..

    llat1 = lat1
    llon1 = lon1 

! If not, check things

    IF( present(degrees) ) THEN
      IF( degrees ) THEN
       llat1 = lat1 * deg2rad
       llon1 = lon1 * deg2rad
      ENDIF
    ENDIF

! Try and be smart, in case user is stupid (like me!) - Note, this check will not work with lats ~ equator,
! nor lons near in western Europe, or mid-Pacific

   IF ( .not. present(degrees) ) THEN
    IF( abs(llat1) .gt. 2.0_r8*PI ) llat1 = llat1 * deg2rad
    IF( abs(llon1) .gt. 2.0_r8*PI ) llon1 = llon1 * deg2rad
   ENDIF
   
    IF (map_proj.eq.0) THEN
      lat = llat1 + y / rearth
      lon = llon1 + x / ( rearth * cos(0.5_r8*(llat1+lat)) )
    ELSE
      write(string1,*) 'Requested map projection unavailable: ', map_proj
      call error_handler(E_ERR,'xy_to_ll',string1,source,revision,revdate)
    ENDIF

    IF( present(degrees) ) THEN      ! USER said convert them back to degrees
      IF( degrees ) THEN
       lat = lat * rad2deg ! / deg2rad
       lon = lon * rad2deg ! / deg2rad
      ENDIF
    ENDIF

   IF ( .not. present(degrees) ) THEN
    IF( abs(lat1) > 2.0_r8*PI ) lat = lat * rad2deg ! / deg2rad   ! User inputs were in degrees (we think)
    IF( abs(lon1) > 2.0_r8*PI ) lon = lon * rad2deg ! / deg2rad
   ENDIF
  return
  end subroutine xy_to_ll




subroutine define_var_dims(myprogvar, ndims, dimids, memberdimid, unlimiteddimid, &
                       nxcdimid, nycdimid, nzcdimid, nxedimid, nyedimid, nzedimid, & 
                       nxc     , nyc     , nzc     , nxe     , nye     , nze      )

type(progvartype),     intent(in)  :: myprogvar
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: dimids
integer,               intent(in)  :: memberdimid, unlimiteddimid
integer,               intent(in)  :: nxcdimid, nycdimid, nzcdimid, nxedimid, nyedimid, nzedimid
integer,               intent(in)  :: nxc     , nyc     , nzc     , nxe     , nye     , nze

select case( myprogvar%storder ) 
case('xyz3d')

      ndims = 5

      dimids(1) = nxcdimid
      dimids(2) = nycdimid
      dimids(3) = nzcdimid
      dimids(4) = memberdimid
      dimids(5) = unlimitedDimid

      if (myprogvar%dimlens(1) == nxe) dimids(1) = nxedimid
      if (myprogvar%dimlens(2) == nye) dimids(2) = nyedimid
      if (myprogvar%dimlens(3) == nze) dimids(3) = nzedimid

case('xy2d')

      ndims = 4

      dimids(1) = nxcdimid
      dimids(2) = nycdimid
      dimids(3) = memberdimid
      dimids(4) = unlimitedDimid

      if (myprogvar%dimlens(1) == nxe) dimids(1) = nxedimid
      if (myprogvar%dimlens(2) == nye) dimids(2) = nyedimid

case('x1d','y1d','z1d')

      ndims = 3

      dimids(1) = nxcdimid
      dimids(2) = memberdimid
      dimids(3) = unlimitedDimid

      if (myprogvar%dimlens(1) == nxe) dimids(1) = nxedimid

case default

      write(string1,*)'unknown storage order '//trim(myprogvar%storder)//& 
                              ' for variable '//trim(myprogvar%varname)
      call error_handler(E_ERR,'define_var_dims',string1,source,revision,revdate)

end select

return
end subroutine define_var_dims


!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
