! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This is the interface between the MITgcm ocean model and DART.

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, SECPERDAY
use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time, &
                             set_calendar_type, GREGORIAN, print_time, print_date, &
                             operator(*),  operator(+), operator(-), &
                             operator(>),  operator(<), operator(/), &
                             operator(/=), operator(<=)
use     location_mod, only : location_type,      get_close_maxdist_init, &
                             get_close_obs_init, get_close_obs, set_location, &
                             VERTISHEIGHT, get_location, vert_is_height, &
                             vert_is_level, vert_is_surface
use    utilities_mod, only : register_module, error_handler, E_ERR, E_WARN, E_MSG, &
                             logfileunit, get_unit, nc_check, do_output, to_upper, &
                             find_namelist_in_file, check_namelist_read, &
                             open_file, file_exist, find_textfile_dims, file_to_text
use     obs_kind_mod, only : QTY_TEMPERATURE, QTY_SALINITY, QTY_U_CURRENT_COMPONENT, &
                             QTY_V_CURRENT_COMPONENT, QTY_SEA_SURFACE_HEIGHT
use mpi_utilities_mod, only: my_task_id
use random_seq_mod,   only : random_seq_type, init_random_seq, random_gaussian


implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          end_model,              &
          static_init_model,      &
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
public :: prog_var_to_vector, vector_to_prog_var, &
          MIT_meta_type, read_meta, write_meta, &
          read_snapshot, write_snapshot, get_gridsize, &
          write_data_namelistfile, set_model_end_time, &
          snapshot_files_to_sv, sv_to_snapshot_files, &
          timestep_to_DARTtime, DARTtime_to_MITtime, &
          DARTtime_to_timestepindex 

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=129) :: msgstring
logical, save :: module_initialized = .false.

! Storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq


!! FIXME: This is horrid ... 'reclen' is compiler-dependent.
!! IBM XLF  -- item_size_direct_access == 4,8
!! IFORT    -- needs compile switch "-assume byterecl"
integer, parameter :: item_size_direct_access = 4

!------------------------------------------------------------------
!
! MITgcm namelist section:  we want to share the 'data' namelist file
! with the model, so we must declare all possible namelist entries to
! avoid getting an error from a valid namelist file.  Most of these
! values are unused in this model_mod code; only a few are needed and
! those are indicated in comments below.
!
!------------------------------------------------------------------
! The time manager namelist variables
! some types/etc come from   <mitsource>/pkg/cal/cal.h
! some useful insight from   cal_set.F, cal_convdate.F
!
! startDate_1 (integer) yyyymmdd   "start date of the integration"
! startDate_2 (integer) hhmmss
!------------------------------------------------------------------

character(len=9) :: TheCalendar = 'gregorian'
! integration start date follows: yyyymmddhhmmss
integer          :: startDate_1 = 19530101
integer          :: startDate_2 =          60000
logical          :: calendarDumps = .false.

NAMELIST /CAL_NML/ TheCalendar, startDate_1, startDate_2, calendarDumps

! FIXME: these namelists should probably be in a separate file, and only
! the actual values needed should be made public, so this isn't so messy.

! must match the value in EEPARAMS.h
integer, parameter :: MAX_LEN_FNAM = 512

! these come from SIZE.h and can vary.  
! should they be namelist controlled?
integer, parameter :: max_nx = 1024
integer, parameter :: max_ny = 1024
integer, parameter :: max_nz = 512
integer, parameter :: max_nr = 512

!-- FIXME: i have not been able to find all of these in the
!-- original source code.  the ones we're using have been
!-- checked that they are the right type.  the others might
!-- cause problems since i don't have a namelist file i can
!-- examine that has the full set of values specified.
!-- 
!-- must match lists declared in ini_parms.f

!--   Time stepping parameters variable declarations
real(r8) :: pickupSuff, &
      deltaT, deltaTClock, deltaTmom, &
      deltaTtracer, dTtracerLev(max_nr), deltaTfreesurf, &
      abEps, alph_AB, beta_AB, &
      tauCD, rCD, &
      baseTime, startTime, endTime, chkPtFreq, &
      dumpFreq, dumpInitAndLast, adjDumpFreq, taveFreq, tave_lastIter, &
      diagFreq, monitorFreq, adjMonitorFreq, pChkPtFreq, cAdjFreq, &
      outputTypesInclusive, &
      tauThetaClimRelax, tauSaltClimRelax, latBandClimRelax, &
      tauThetaClimRelax3Dim, tauSaltClimRelax3Dim, tauTr1ClimRelax, &
      periodicExternalForcing, externForcingPeriod, externForcingCycle
integer :: nIter0, nTimeSteps, nEndIter, momForcingOutAB, tracForcingOutAB
logical :: forcing_In_AB, &
      momDissip_In_AB, doAB_onGtGs, &
      startFromPickupAB2

!--   Gridding parameters variable declarations 
logical :: usingCartesianGrid, usingCylindricalGrid, &
           usingSphericalPolarGrid, usingCurvilinearGrid, &
           deepAtmosphere
real(r8) :: dxSpacing, dySpacing, delX(max_nx), delY(max_ny), &
            ygOrigin, xgOrigin, rSphere, &
            Ro_SeaLevel, delZ(max_nz), delP, delR(max_nr), delRc(max_nr+1), &
            rkFac, groundAtK1
character(len=MAX_LEN_FNAM) :: delXFile, delYFile, &
                      delRFile, delRcFile, &
                      horizGridFile

!--   Input files variable declarations 
character(len=MAX_LEN_FNAM) :: &
      bathyFile, topoFile, shelfIceFile, &
      hydrogThetaFile, hydrogSaltFile, diffKrFile, &
      zonalWindFile, meridWindFile, &
      thetaClimFile, saltClimFile, &
      surfQfile, surfQnetFile, surfQswFile, EmPmRfile, saltFluxFile, &
      lambdaThetaFile, lambdaSaltFile, &
      uVelInitFile, vVelInitFile, pSurfInitFile, &
      dQdTFile, ploadFile,tCylIn,tCylOut, &
      eddyTauxFile, eddyTauyFile, &
      mdsioLocalDir, &
      the_run_name

!--   Time stepping parameters namelist
NAMELIST /PARM03/ &
      nIter0, nTimeSteps, nEndIter, pickupSuff, &
      deltaT, deltaTClock, deltaTmom, &
      deltaTtracer, dTtracerLev, deltaTfreesurf, &
      forcing_In_AB, momForcingOutAB, tracForcingOutAB, &
      momDissip_In_AB, doAB_onGtGs, &
      abEps, alph_AB, beta_AB, startFromPickupAB2, &
      tauCD, rCD, &
      baseTime, startTime, endTime, chkPtFreq, &
      dumpFreq, dumpInitAndLast, adjDumpFreq, taveFreq, tave_lastIter, &
      diagFreq, monitorFreq, adjMonitorFreq, pChkPtFreq, cAdjFreq, &
      outputTypesInclusive, &
      tauThetaClimRelax, tauSaltClimRelax, latBandClimRelax, &
      tauThetaClimRelax3Dim, tauSaltClimRelax3Dim, tauTr1ClimRelax, &
      periodicExternalForcing, externForcingPeriod, externForcingCycle, &
      calendarDumps

!--   Gridding parameters namelist
NAMELIST /PARM04/ &
      usingCartesianGrid, usingCylindricalGrid, &
      dxSpacing, dySpacing, delX, delY, delXFile, delYFile, &
      usingSphericalPolarGrid, ygOrigin, xgOrigin, rSphere, &
      usingCurvilinearGrid, horizGridFile, deepAtmosphere, &
      Ro_SeaLevel, delZ, delP, delR, delRc, delRFile, delRcFile, &
      rkFac, groundAtK1

!--   Input files namelist
NAMELIST /PARM05/ &
      bathyFile, topoFile, shelfIceFile, &
      hydrogThetaFile, hydrogSaltFile, diffKrFile, &
      zonalWindFile, meridWindFile, &
      thetaClimFile, saltClimFile, &
      surfQfile, surfQnetFile, surfQswFile, EmPmRfile, saltFluxFile, &
      lambdaThetaFile, lambdaSaltFile, &
      uVelInitFile, vVelInitFile, pSurfInitFile, &
      dQdTFile, ploadFile,tCylIn,tCylOut, &
      eddyTauxFile, eddyTauyFile, &
      mdsioLocalDir, &
      the_run_name

!------------------------------------------------------------------
!
! The DART state vector (control vector) will consist of:  S, T, U, V, Eta
! (Salinity, Temperature, U velocity, V velocity, Sea Surface Height).
! S, T are 3D arrays, located at cell centers.  U is staggered in X
! and V is staggered in Y (meaning the points are located on the cell
! faces) but the grids are offset by half a cell, so there are actually
! the same number of points in each grid. 
! Eta is a 2D field (X,Y only).  The Z direction is downward.
!
!------------------------------------------------------------------

integer, parameter :: n3dfields = 4
integer, parameter :: n2dfields = 1
integer, parameter :: nfields   = n3dfields + n2dfields

integer, parameter :: S_index   = 1
integer, parameter :: T_index   = 2
integer, parameter :: U_index   = 3
integer, parameter :: V_index   = 4
integer, parameter :: Eta_index = 5

! (the absoft compiler likes them to all be the same length during declaration)
! we trim the blanks off before use anyway, so ...
character(len=128) :: progvarnames(nfields) = (/'S  ','T  ','U  ','V  ','Eta'/)

integer :: start_index(nfields)

! Grid parameters - the values will be read from a
! standard MITgcm namelist and filled in here.

integer :: Nx=-1, Ny=-1, Nz=-1    ! grid counts for each field

! locations of cell centers (C) and edges (G) for each axis.
real(r8), allocatable :: XC(:), XG(:), YC(:), YG(:), ZC(:), ZG(:)

! location information - these grids can either be regularly
! spaced or the spacing along each axis can vary.

!real(r8) :: lat_origin, lon_origin
!logical  :: regular_lat, regular_lon, regular_depth
!real(r8) :: delta_lat, delta_lon, delta_depth
!real(r8), allocatable :: lat_grid(:), lon_grid(:), depth_grid(:)

real(r8)        :: ocean_dynamics_timestep = 900.0_r4
integer         :: timestepcount = 0
type(time_type) :: model_time, model_timestep

integer :: model_size    ! the state vector length

! Skeleton of a model_nml that would be in input.nml
! This is where dart-related model parms could be set.
logical  :: output_state_vector = .false.
integer  :: assimilation_period_days = 7
integer  :: assimilation_period_seconds = 0
real(r8) :: model_perturbation_amplitude = 0.2

namelist /model_nml/ assimilation_period_days,    &
                     assimilation_period_seconds, &
                     model_perturbation_amplitude

! /pkg/mdsio/mdsio_write_meta.F writes the .meta files 
type MIT_meta_type
!  private
   integer :: nDims
   integer :: dimList(3)
   character(len=32) :: dataprec
   integer :: reclen
   integer :: nrecords
   integer :: timeStepNumber    ! optional
end type MIT_meta_type

INTERFACE read_snapshot
      MODULE PROCEDURE read_2d_snapshot
      MODULE PROCEDURE read_3d_snapshot
END INTERFACE

INTERFACE write_snapshot
      MODULE PROCEDURE write_2d_snapshot
      MODULE PROCEDURE write_3d_snapshot
END INTERFACE

INTERFACE vector_to_prog_var
      MODULE PROCEDURE vector_to_2d_prog_var
      MODULE PROCEDURE vector_to_3d_prog_var
END INTERFACE


contains

!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------
!
! Called to do one time initialization of the model. In this case,
! it reads in the grid information and then the model data.

integer :: i, iunit, io
integer :: ss, dd

! The Plan:
!
!   read the standard MITgcm namelist file 'data.cal' for the calendar info
!
!   read the standard MITgcm namelist file 'data' for the
!   time stepping info and the grid info.
!
!   open the grid data files to get the actual grid coordinates
!
!   Compute the model size.
!
!   set the index numbers where the field types change
!
!   set the grid location info
!

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.


! Read the DART namelist for this model
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run
call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
if (do_output()) write(logfileunit, nml=model_nml)
if (do_output()) write(     *     , nml=model_nml)

! MIT calendar information
call find_namelist_in_file("data.cal", "CAL_NML", iunit)
read(iunit, nml = CAL_NML, iostat = io)
call check_namelist_read(iunit, io, "CAL_NML")

if (index(TheCalendar,'g') > 0 ) then
   call set_calendar_type(GREGORIAN)
elseif (index(TheCalendar,'G') > 0 )then
   call set_calendar_type(GREGORIAN)
else
   write(msgstring,*)"namelist data.cal indicates a ",trim(TheCalendar)," calendar."
   call error_handler(E_MSG,"static_init_model", msgstring, source, revision, revdate)
   write(msgstring,*)"only have support for Gregorian"
   call error_handler(E_ERR,"static_init_model", msgstring, source, revision, revdate)
endif
if (do_output()) write(*,*)'model_mod:namelist cal_NML',startDate_1,startDate_2
if (do_output()) write(*,nml=CAL_NML)

! Time stepping parameters are in PARM03
call find_namelist_in_file("data", "PARM03", iunit)
read(iunit, nml = PARM03, iostat = io)
call check_namelist_read(iunit, io, "PARM03")

if ((deltaTmom   == deltaTtracer) .and. &
    (deltaTmom   == deltaTClock ) .and. &
    (deltaTClock == deltaTtracer)) then
   ocean_dynamics_timestep = deltaTmom                    ! need a time_type version
else
   write(msgstring,*)"namelist PARM03 has deltaTmom /= deltaTtracer /= deltaTClock"
   call error_handler(E_MSG,"static_init_model", msgstring, source, revision, revdate)
   write(msgstring,*)"values were ",deltaTmom, deltaTtracer, deltaTClock
   call error_handler(E_MSG,"static_init_model", msgstring, source, revision, revdate)
   write(msgstring,*)"At present, DART only supports equal values."
   call error_handler(E_ERR,"static_init_model", msgstring, source, revision, revdate)
endif

! Define the assimilation period as the model_timestep
! Ensure model_timestep is multiple of ocean_dynamics_timestep

model_time     = timestep_to_DARTtime(timestepcount)
model_timestep = set_model_time_step(assimilation_period_seconds, &
                                     assimilation_period_days,    &
                                     ocean_dynamics_timestep)

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(msgstring,*)"assimilation period is ",dd," days ",ss," seconds"
call error_handler(E_MSG,'static_init_model',msgstring,source,revision,revdate)
if (do_output()) write(logfileunit,*)msgstring

! Grid-related variables are in PARM04
delX(:) = 0.0_r4
delY(:) = 0.0_r4
delZ(:) = 0.0_r4
call find_namelist_in_file("data", "PARM04", iunit)
read(iunit, nml = PARM04, iostat = io)
call check_namelist_read(iunit, io, "PARM04")

! Input datasets are in PARM05
call find_namelist_in_file("data", "PARM05", iunit)
read(iunit, nml = PARM05, iostat = io)
call check_namelist_read(iunit, io, "PARM05")

! The only way I know to compute the number of
! levels/lats/lons is to set the default value of delZ to 0.0
! before reading the namelist.  now loop until you get back
! to zero and that is the end of the list.
! Not a very satisfying/robust solution ...

Nx = -1
do i=1, size(delX)
 if (delX(i) == 0.0_r4) then
    Nx = i-1
    exit
 endif
enddo
if (Nx == -1) then
   write(msgstring,*)'could not figure out number of longitudes from delX in namelist'
   call error_handler(E_ERR,"static_init_model", msgstring, source, revision, revdate)
endif

Ny = -1
do i=1, size(delY)
 if (delY(i) == 0.0_r4) then
    Ny = i-1
    exit
 endif
enddo
if (Ny == -1) then
   write(msgstring,*)'could not figure out number of latitudes from delY in namelist'
   call error_handler(E_ERR,"static_init_model", msgstring, source, revision, revdate)
endif

Nz = -1
do i=1, size(delZ)
 if (delZ(i) == 0.0_r4) then
    Nz = i-1
    exit
 endif
enddo
if (Nz == -1) then
   write(msgstring,*)'could not figure out number of depth levels from delZ in namelist'
   call error_handler(E_ERR,"static_init_model", msgstring, source, revision, revdate)
endif

! We know enough to allocate grid variables. 

allocate(XC(Nx), YC(Ny), ZC(Nz))
allocate(XG(Nx), YG(Ny), ZG(Nz))

! XG (the grid edges) and XC (the grid centroids) must be computed.

XG(1) = xgOrigin
XC(1) = xgOrigin + 0.5_r8 * delX(1)
do i=2, Nx
 XG(i) = XG(i-1) + delX(i-1)
 XC(i) = XC(i-1) + 0.5_r8 * delX(i-1) + 0.5_r8 * delX(i) 
enddo

! YG (the grid edges) and YC (the grid centroids) must be computed.

YG(1) = ygOrigin
YC(1) = ygOrigin + 0.5_r8 * delY(1)
do i=2, Ny
 YG(i) = YG(i-1) + delY(i-1)
 YC(i) = YC(i-1) + 0.5_r8 * delY(i-1) + 0.5_r8 * delY(i) 
enddo

! the namelist contains a list of thicknesses of each depth level (delZ)
! ZG (the grid edges) and ZC (the grid centroids) must be computed.

ZG(1) = 0.0_r8
ZC(1) = -0.5_r8 * delZ(1)
do i=2, Nz
 ZG(i) = ZG(i-1) - delZ(i-1)
 ZC(i) = ZC(i-1) - 0.5_r8 * delZ(i-1) - 0.5_r8 * delZ(i) 
enddo

! record where in the state vector the data type changes
! from one type to another, by computing the starting
! index for each block of data.
start_index(S_index)   = 1
start_index(T_index)   = start_index(S_index) + (Nx * Ny * Nz)
start_index(U_index)   = start_index(T_index) + (Nx * Ny * Nz)
start_index(V_index)   = start_index(U_index) + (Nx * Ny * Nz)
start_index(Eta_index) = start_index(V_index) + (Nx * Ny * Nz)

! in spite of the staggering, all grids are the same size
! and offset by half a grid cell.  4 are 3D and 1 is 2D.
!  e.g. S,T,U,V = 256 x 225 x 70
!  e.g. Eta = 256 x 225

if (do_output()) write(logfileunit, *) 'Using grid size : '
if (do_output()) write(logfileunit, *) '  Nx, Ny, Nz = ', Nx, Ny, Nz
if (do_output()) write(     *     , *) 'Using grid size : '
if (do_output()) write(     *     , *) '  Nx, Ny, Nz = ', Nx, Ny, Nz
!print *, ' 3d field size: ', n3dfields * (Nx * Ny * Nz)
!print *, ' 2d field size: ', n2dfields * (Nx * Ny)
model_size = (n3dfields * (Nx * Ny * Nz)) + (n2dfields * (Nx * Ny))
if (do_output()) write(*,*) 'model_size = ', model_size

end subroutine static_init_model




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



subroutine adv_1step(x, time)
!------------------------------------------------------------------
!
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



function get_model_size()
!------------------------------------------------------------------
!
! Returns the size of the model as an integer. Required for all
! applications.

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size



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



subroutine model_interpolate(x, location, obs_type, interp_val, istatus)
!=======================================================================
!

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_type
real(r8),           intent(out) :: interp_val
integer,            intent(out) :: istatus

! Model interpolate will interpolate any state variable (S, T, U, V, Eta) to
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

! Successful istatus is 0
interp_val = 0.0_r8
istatus = 0

! Get the individual locations values
loc_array = get_location(location)
llon    = loc_array(1)
llat    = loc_array(2)
lheight = loc_array(3)

if( vert_is_height(location) ) then
   ! Nothing to do 
elseif ( vert_is_surface(location) ) then
   ! Nothing to do 
elseif (vert_is_level(location)) then
   ! convert the level index to an actual depth 
   ind = nint(loc_array(3))
   if ( (ind < 1) .or. (ind > size(zc)) ) then 
      lheight = zc(ind)
   else
      istatus = 1
      return
   endif
else   ! if pressure or undefined, we don't know what to do
   istatus = 7
   return
endif

! Do horizontal interpolations for the appropriate levels
! Find the basic offset of this field
if(obs_type == QTY_SALINITY) then
   base_offset = start_index(1)
else if(obs_type == QTY_TEMPERATURE) then
   base_offset = start_index(2)
else if(obs_type == QTY_U_CURRENT_COMPONENT) then
   base_offset = start_index(3)
else if(obs_type == QTY_V_CURRENT_COMPONENT) then
   base_offset = start_index(4)
else if(obs_type == QTY_SEA_SURFACE_HEIGHT) then
   base_offset = start_index(5)
else
   ! Not a legal type for interpolation, return istatus error
   istatus = 5
   return
endif

!print *, 'base offset now ', base_offset

! For Sea Surface Height don't need the vertical coordinate
if( vert_is_surface(location) ) then
   call lat_lon_interpolate(x(base_offset:), llon, llat, obs_type, interp_val, istatus)
   return
endif

! Get the bounding vertical levels and the fraction between bottom and top
call height_bounds(lheight, nz, zc, hgt_bot, hgt_top, hgt_fract, hstatus)
if(hstatus /= 0) then
   istatus = 2
   return
endif

! Find the base location for the top height and interpolate horizontally on this level
offset = base_offset + (hgt_top - 1) * nx * ny
!print *, 'relative top height offset = ', offset - base_offset
!print *, 'absolute top height offset = ', offset
call lat_lon_interpolate(x(offset:), llon, llat, obs_type, top_val, istatus)
! Failed istatus from interpolate means give up
if(istatus /= 0) return

! Find the base location for the bottom height and interpolate horizontally on this level
offset = base_offset + (hgt_bot - 1) * nx * ny
!print *, 'relative bot height offset = ', offset - base_offset
!print *, 'absolute bot height offset = ', offset
call lat_lon_interpolate(x(offset:), llon, llat, obs_type, bot_val, istatus)
! Failed istatus from interpolate means give up
if(istatus /= 0) return

! Then weight them by the fraction and return
interp_val = bot_val + hgt_fract * (top_val - bot_val)
!print *, 'model_interp: interp val = ',interp_val

end subroutine model_interpolate



subroutine height_bounds(lheight, nheights, hgt_array, bot, top, fract, istatus)
!=======================================================================
!
real(r8),             intent(in) :: lheight
integer,              intent(in) :: nheights
real(r8),             intent(in) :: hgt_array(nheights)
integer,             intent(out) :: bot, top
real(r8),            intent(out) :: fract
integer,             intent(out) :: istatus

! Local variables
integer   :: i

if ( .not. module_initialized ) call static_init_model

! Succesful istatus is 0
istatus = 0

! The zc array contains the depths of the center of the vertical grid boxes

! It is assumed that the top box is shallow and any observations shallower
! than the depth of this boxes center are just given the value of the
! top box.
if(lheight > hgt_array(1)) then
   top = 1
   bot = 2
   ! NOTE: the fract definition is the relative distance from bottom to top
   ! ??? Make sure this is consistent with the interpolation
   fract = 1.0_r8 
endif

! Search through the boxes
do i = 2, nheights
   ! If the location is shallower than this entry, it must be in this box
   if(lheight >= hgt_array(i)) then
      top = i -1
      bot = i
      fract = (lheight - hgt_array(bot)) / (hgt_array(top) - hgt_array(bot))
      return
   endif
end do

! Falling off the end means the location is lower than the deepest height
! Fail with istatus 2 in this case
istatus = 2

end subroutine height_bounds


subroutine lat_lon_interpolate(x, llon, llat, var_type, interp_val, istatus)
!=======================================================================
!

! Subroutine to interpolate to a lat lon location given that portion of state
! vector. 
! NOTE: Using array sections to pass in the x array may be inefficient on some
! compiler/platform setups. Might want to pass in the entire array with a base
! offset value instead of the section if this is an issue.

real(r8),            intent(in) :: x(:)
real(r8),            intent(in) :: llon, llat
integer,             intent(in) :: var_type
integer,            intent(out) :: istatus
real(r8),           intent(out) :: interp_val

! Local storage
real(r8) :: lat_array(ny), lon_array(nx)
integer  :: lat_bot, lat_top, lon_bot, lon_top
real(r8) :: lat_fract, lon_fract
real(r8) :: pa, pb, pc, pd, xbot, xtop
integer  :: lat_status, lon_status
logical  :: masked

if ( .not. module_initialized ) call static_init_model

! Succesful return has istatus of 0
istatus = 0

! Failed return values for istatus are ?????

! Find out what latitude box and fraction
! The latitude grid being used depends on the variable type
! V is on the YG latitude grid
if(var_type == QTY_V_CURRENT_COMPONENT) then
   lat_array = yg
   call lat_bounds(llat, ny, lat_array, lat_bot, lat_top, lat_fract, lat_status)
else 
   ! Eta, U, T and S are on the YC latitude grid
   lat_array = yc
   call lat_bounds(llat, ny, lat_array, lat_bot, lat_top, lat_fract, lat_status)
endif

! Check for error on the latitude interpolation
if(lat_status /= 0) then 
   istatus = 1
   return
endif

! Find out what longitude box and fraction
if(var_type == QTY_U_CURRENT_COMPONENT) then
   ! U velocity is on the XG grid
   lon_array = xg
   call lon_bounds(llon, nx, lon_array, lon_bot, lon_top, lon_fract, lon_status)
else
   ! Eta, V, T, and S are on the XC grid
   lon_array = xc
   call lon_bounds(llon, nx, lon_array, lon_bot, lon_top, lon_fract, lon_status)
endif

! Check for error on the longitude interpolation
if(lat_status /= 0) then 
   istatus = 2
   return
endif

! Vector is laid out with lat outermost loop, lon innermost loop
! Find the bounding points for the lat lon box
! NOTE: For now, it is assumed that a real(r8) value of exactly 0.0 indicates
! that a particular gridded quantity is masked and not available. This is not
! the most robust way to do this, but may be sufficient since exact 0's are
! expected to happen rarely. Jeff Anderson believes that the only implication
! will be that an observation whos forward operator requires interpolating
! from a point that has exactly 0.0 (but is not masked) will not be 
! assimilated.
!print *, 'lon_bot, lon_top = ', lon_bot, lon_top
!print *, 'lat_bot, lat_top = ', lat_bot, lat_top

pa = get_val(lon_bot, lat_bot, nx, x, masked)
!print *, 'pa = ', pa
if(masked) then
   istatus = 3
   return
endif
pb = get_val(lon_top, lat_bot, nx, x, masked)
!print *, 'pb = ', pb
if(masked) then
   istatus = 3
   return
endif
pc = get_val(lon_bot, lat_top, nx, x, masked)
!print *, 'pc = ', pc
if(masked) then
   istatus = 3
   return
endif
pd = get_val(lon_top, lat_top, nx, x, masked)
!print *, 'pd = ', pd
if(masked) then
   istatus = 3
   return
endif

!print *, 'pa,b,c,d = ', pa, pb, pc, pd

! Finish bi-linear interpolation 
! First interpolate in longitude
!print *, 'bot lon_fract, delta = ', lon_fract, (pb-pa)
xbot = pa + lon_fract * (pb - pa)
!print *, 'xbot = ', xbot
!print *, 'top lon_fract, delta = ', lon_fract, (pd-pc)
xtop = pc + lon_fract * (pd - pc)
!print *, 'xtop = ', xtop
! Now interpolate in latitude
!print *, 'lat_fract, delta = ', lat_fract, (xtop - xbot)
interp_val = xbot + lat_fract * (xtop - xbot)
!print *, 'lat_lon_interp: interp_val = ', interp_val

end subroutine lat_lon_interpolate




subroutine lat_bounds(llat, nlats, lat_array, bot, top, fract, istatus)

!=======================================================================
!

! Given a latitude llat, the array of latitudes for grid boundaries, and the
! number of latitudes in the grid, returns the indices of the latitude
! below and above the location latitude and the fraction of the distance
! between. istatus is returned as 0 unless the location latitude is 
! south of the southernmost grid point (1 returned) or north of the 
! northernmost (2 returned). If one really had lots of polar obs would 
! want to worry about interpolating around poles.

real(r8),          intent(in) :: llat
integer,           intent(in) :: nlats
real(r8),          intent(in) :: lat_array(nlats)
integer,          intent(out) :: bot, top
real(r8),         intent(out) :: fract
integer,          intent(out) :: istatus

! Local storage
integer    :: i

if ( .not. module_initialized ) call static_init_model

! Default is success
istatus = 0

! Check for too far south or north
if(llat < lat_array(1)) then
   istatus = 1
   return
else if(llat > lat_array(nlats)) then
   istatus = 2
   return
endif

! In the middle, search through
do i = 2, nlats
   if(llat <= lat_array(i)) then
      bot = i - 1
      top = i
      fract = (llat - lat_array(bot)) / (lat_array(top) - lat_array(bot))
      return
   endif
end do

end subroutine lat_bounds



subroutine lon_bounds(llon, nlons, lon_array, bot, top, fract, istatus)

!=======================================================================
!

! Given a longitude llon, the array of longitudes for grid boundaries, and the
! number of longitudes in the grid, returns the indices of the longitude
! below and above the location longitude and the fraction of the distance
! between. istatus is returned as 0 unless the location longitude is 
! not between any of the longitude box boundaries. This should be modified
! for global wrap-around grids.
! Algorithm fails for a silly grid that
! has only two longitudes separated by 180 degrees.

real(r8),          intent(in) :: llon
integer,           intent(in) :: nlons
real(r8),          intent(in) :: lon_array(nlons)
integer,          intent(out) :: bot, top
real(r8),         intent(out) :: fract
integer,          intent(out) :: istatus

! Local storage
integer  :: i
real(r8) :: dist_bot, dist_top

if ( .not. module_initialized ) call static_init_model

! Default is success
istatus = 0

!print *, 'computing bounds for = ', llon
! This is inefficient, someone could clean it up
! Plus, it doesn't work for a global model that wraps around
do i = 2, nlons
   dist_bot = lon_dist(llon, lon_array(i - 1))
   dist_top = lon_dist(llon, lon_array(i))
!print *, 'dist top, bot = ', dist_top, dist_bot
   if(dist_bot >= 0 .and. dist_top < 0) then
      bot = i - 1
      top = i
!print *, 'bot, top = ', bot, top
!print *, 'numerator = ',  dist_bot
!print *, 'denomenator = ', dist_bot + abs(dist_top)
      fract = dist_bot / (dist_bot + abs(dist_top))
      ! orig: fract = abs(dist_bot) / (abs(dist_bot) + dist_top)
!print *, 'fract = ', fract
      return
   endif
end do

! Falling off the end means its in between. Add the wraparound check.
! For now, return istatus 1
istatus = 1

end subroutine lon_bounds



function lon_dist(lon1, lon2)
!=======================================================================
!

! Returns the smallest signed distance between lon1 and lon2 on the sphere
! If lon1 is less than 180 degrees east of lon2 the distance is negative
! If lon1 is less than 180 degrees west of lon2 the distance is positive

real(r8), intent(in) :: lon1, lon2
real(r8)             :: lon_dist

if ( .not. module_initialized ) call static_init_model

lon_dist = lon1 - lon2
if(lon_dist >= -180.0_r8 .and. lon_dist <= 180.0_r8) then 
   return
else if(lon_dist < -180.0_r8) then
   lon_dist = lon_dist + 360.0_r8
else
   lon_dist = lon_dist - 360.0_r8
endif

end function lon_dist


function get_val(lon_index, lat_index, nlon, x, masked)
!=======================================================================
!

! Returns the value from a single level array given the lat and lon indices
integer,     intent(in) :: lon_index, lat_index, nlon
real(r8),    intent(in) :: x(:)
logical,    intent(out) :: masked
real(r8)                :: get_val

if ( .not. module_initialized ) call static_init_model

! Layout has lons varying most rapidly
!print *, 'lat_index, lon_index, nlon', lat_index, lon_index, nlon
!print *, 'computing index val ', (lat_index - 1) * nlon + lon_index
get_val = x((lat_index - 1) * nlon + lon_index)
!print *, 'get_val = ', get_val

! Masked returns false if the value is masked
! A grid variable is assumed to be masked if its value is exactly 0.
! See discussion in lat_lon_interpolate.
if(get_val == 0.0_r8) then
   masked = .true.
else
   masked = .false.
endif

!print *, 'masked is ', masked
end function get_val



function set_model_time_step(ss, dd, dt)
!------------------------------------------------------------------
!
! Sets the model 'timestep' AKA 'assimilation period'.
! Must make sure the assimilation period is a multiple of the 
! model's dynamical timestepping requirement.

integer,  intent(in) :: ss    ! assimilation_period_seconds
integer,  intent(in) :: dd    ! assimilation_period_days
real(r8), intent(in) :: dt    ! ocean_dynamics_timestep

type(time_type) :: set_model_time_step

integer :: assim_period, ndt

if ( .not. module_initialized ) call static_init_model

assim_period = ss + dd*SECPERDAY   ! in seconds 
ndt = max(nint(assim_period / dt),1)
assim_period = nint(ndt * dt)

set_model_time_step = set_time(assim_period, 0) ! works seconds > 86400

end function set_model_time_step



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



subroutine set_model_end_time(adv_to_offset)
!------------------------------------------------------------------
!
! sets PARM03:endTime to reflect the time when the model will stop. 
! endTime is in module storage 

type(time_type), intent(in) :: adv_to_offset

integer :: secs, days

if ( .not. module_initialized ) call static_init_model

call get_time(adv_to_offset, secs, days)

endTime = (secs + days*SECPERDAY)

end subroutine set_model_end_time



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

real(R8) :: lat, lon, depth
integer :: var_num, offset, lon_index, lat_index, depth_index

if ( .not. module_initialized ) call static_init_model

!print *, 'asking for meta data about index ', index_in

if (index_in < start_index(S_index+1)) then
   if (present(var_type)) var_type = QTY_SALINITY  
   var_num = S_index
else if (index_in < start_index(T_index+1)) then
   if (present(var_type)) var_type = QTY_TEMPERATURE  
   var_num = T_index
else if (index_in < start_index(U_index+1)) then
   if (present(var_type)) var_type = QTY_U_CURRENT_COMPONENT
   var_num = U_index
else if (index_in < start_index(V_index+1)) then
   if (present(var_type)) var_type = QTY_V_CURRENT_COMPONENT
   var_num = V_index
else 
   if (present(var_type)) var_type = QTY_SEA_SURFACE_HEIGHT
   var_num = Eta_index
endif

!print *, 'var num = ', var_num

! local offset into this var array
offset = index_in - start_index(var_num)

!print *, 'offset = ', offset

if (var_num == Eta_index) then
  depth = 0.0
  depth_index = 1
else
  depth_index = (offset / (Nx * Ny)) + 1
  depth = ZC(depth_index)
endif

lat_index = (offset - ((depth_index-1)*Nx*Ny)) / Nx + 1
lon_index = offset - ((depth_index-1)*Nx*Ny) - ((lat_index-1)*Nx) + 1

!print *, 'lon, lat, depth index = ', lon_index, lat_index, depth_index
lon = XC(lon_index)
lat = YC(lat_index)

!print *, 'lon, lat, depth = ', lon, lat, depth

location = set_location(lon, lat, depth, VERTISHEIGHT)

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

! good style ... perhaps you could deallocate stuff (from static_init_model?).
! deallocate(state_loc)

if ( .not. module_initialized ) call static_init_model

end subroutine end_model



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

use typeSizes
use netcdf

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
integer :: XGDimID, XCDimID, YGDimID, YCDimID, ZGDimID, ZCDimID
integer :: XGVarID, XCVarID, YGVarID, YCVarID, ZGVarID, ZCVarID

! for the prognostic variables
integer :: SVarID, TVarID, UVarID, VVarID, EtaVarID 

!----------------------------------------------------------------------
! variables for the namelist output
!----------------------------------------------------------------------

character(len=129), allocatable, dimension(:) :: textblock
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen

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

write(filename,*) 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! and then put into define mode.
!-------------------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
                                   "nc_write_model_atts", "inquire "//trim(filename))
call nc_check(nf90_Redef(ncFileID),"nc_write_model_atts",   "redef "//trim(filename))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension. 
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name="NMLlinelen", dimid=LineLenDimID), &
                           "nc_write_model_atts","inq_dimid NMLlinelen")
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID), &
                           "nc_write_model_atts", "copy dimid "//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID), &
                           "nc_write_model_atts", "time dimid "//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(msgstring,*)"Time Dimension ID ",TimeDimID, &
             " should equal Unlimited Dimension ID",unlimitedDimID
   call error_handler(E_ERR,"nc_write_model_atts", msgstring, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------
call nc_check(nf90_def_dim(ncid=ncFileID, name="StateVariable", len=model_size, &
        dimid = StateVarDimID),"nc_write_model_atts", "state def_dim "//trim(filename))

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date" ,str1    ), &
           "nc_write_model_atts", "creation put "//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source"  ,source  ), &
           "nc_write_model_atts", "source put "//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision), &
           "nc_write_model_atts", "revision put "//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate" ,revdate ), &
           "nc_write_model_atts", "revdate put "//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model",  "MITgcm_ocean" ), &
           "nc_write_model_atts", "model put "//trim(filename))

!-------------------------------------------------------------------------------
! Determine shape of most important namelist
!-------------------------------------------------------------------------------

call find_textfile_dims("data", nlines, linelen)

allocate(textblock(nlines))
textblock = ''

call nc_check(nf90_def_dim(ncid=ncFileID, name="nlines", &
              len = nlines, dimid = nlinesDimID), &
              'nc_write_model_atts', 'def_dim nlines ')

call nc_check(nf90_def_var(ncFileID,name="datanml", xtype=nf90_char,    &
              dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID), &
              'nc_write_model_atts', 'def_var datanml')
call nc_check(nf90_put_att(ncFileID, nmlVarID, "long_name",       &
              "contents of data namelist"), 'nc_write_model_atts', 'put_att datanml')

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
   call nc_check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
                 dimids=StateVarDimID, varid=StateVarVarID), "nc_write_model_atts", &
                 "statevariable def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarVarID,"long_name","State Variable ID"),&
                 "nc_write_model_atts","statevariable long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "units","indexical"), &
                 "nc_write_model_atts", "statevariable units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarVarID,"valid_range",(/ 1,model_size /)),&
                 "nc_write_model_atts", "statevariable valid_range "//trim(filename))

   ! Define the actual (3D) state vector, which gets filled as time goes on ... 
   call nc_check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_real, &
                 dimids=(/StateVarDimID,MemberDimID,unlimitedDimID/),varid=StateVarID),&
                 "nc_write_model_atts","state def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarID,"long_name","model state or fcopy"),&
                 "nc_write_model_atts", "state long_name "//trim(filename))

   ! Leave define mode so we can fill the coordinate variable.
   call nc_check(nf90_enddef(ncfileID),"nc_write_model_atts","state enddef "//trim(filename))

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ), &
                 "nc_write_model_atts", "state put_var "//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to output the prognostic variables.
   !----------------------------------------------------------------------------
   ! Define the new dimensions IDs
   !----------------------------------------------------------------------------
   
   call nc_check(nf90_def_dim(ncid=ncFileID, name="XG", &
          len = Nx, dimid = XGDimID),"nc_write_model_atts", "XG def_dim "//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name="XC", &
          len = Nx, dimid = XCDimID),"nc_write_model_atts", "XC def_dim "//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name="YG", &
          len = Ny, dimid = YGDimID),"nc_write_model_atts", "YG def_dim "//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name="YC", &
          len = Ny, dimid = YCDimID),"nc_write_model_atts", "YC def_dim "//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name="ZG", &
          len = Nz, dimid = ZGDimID),"nc_write_model_atts", "ZG def_dim "//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name="ZC", &
          len = Nz, dimid = ZCDimID),"nc_write_model_atts", "ZC def_dim "//trim(filename))
   
   !----------------------------------------------------------------------------
   ! Create the (empty) Coordinate Variables and the Attributes
   !----------------------------------------------------------------------------

   ! U Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name="XG",xtype=nf90_real,dimids=XGDimID,varid=XGVarID),&
                 "nc_write_model_atts", "XG def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  XGVarID, "long_name", "longitude grid edges"), &
                 "nc_write_model_atts", "XG long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  XGVarID, "cartesian_axis", "X"),  &
                 "nc_write_model_atts", "XG cartesian_axis "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  XGVarID, "units", "degrees_east"), &
                 "nc_write_model_atts", "XG units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  XGVarID, "valid_range", (/ 0.0_r8, 360.0_r8 /)), &
                 "nc_write_model_atts", "XG valid_range "//trim(filename))

   ! S,T,V,Eta Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name="XC",xtype=nf90_real,dimids=XCDimID,varid=XCVarID),&
                 "nc_write_model_atts", "XC def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, XCVarID, "long_name", "longitude grid centroids"), &
                 "nc_write_model_atts", "XC long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, XCVarID, "cartesian_axis", "X"),   &
                 "nc_write_model_atts", "XC cartesian_axis "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, XCVarID, "units", "degrees_east"),  &
                 "nc_write_model_atts", "XC units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, XCVarID, "valid_range", (/ 0.0_r8, 360.0_r8 /)), &
                 "nc_write_model_atts", "XC valid_range "//trim(filename))

   ! V Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name="YG",xtype=nf90_real,dimids=YGDimID,varid=YGVarID),&
                 "nc_write_model_atts", "YG def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YGVarID, "long_name", "latitude grid edges"), &
                 "nc_write_model_atts", "YG long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YGVarID, "cartesian_axis", "Y"),   &
                 "nc_write_model_atts", "YG cartesian_axis "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YGVarID, "units", "degrees_north"),  &
                 "nc_write_model_atts", "YG units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,YGVarID,"valid_range",(/-90.0_r8,90.0_r8 /)), &
                 "nc_write_model_atts", "YG valid_range "//trim(filename))

   ! S,T,U,Eta Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name="YC",xtype=nf90_real,dimids=YCDimID,varid=YCVarID), &
                 "nc_write_model_atts", "YC def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YCVarID, "long_name", "latitude grid centroids"), &
                 "nc_write_model_atts", "YC long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YCVarID, "cartesian_axis", "Y"),   &
                 "nc_write_model_atts", "YC cartesian_axis "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YCVarID, "units", "degrees_north"),  &
                 "nc_write_model_atts", "YC units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YCVarID, "valid_range", (/ -90.0_r8, 90.0_r8 /)), &
                 "nc_write_model_atts", "YC valid_range "//trim(filename))

   ! Depths
   call nc_check(nf90_def_var(ncFileID,name="ZG",xtype=nf90_real,dimids=ZGDimID,varid=ZGVarID), &
                 "nc_write_model_atts", "ZG def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZGVarID, "long_name", "depth at grid edges"), &
                 "nc_write_model_atts", "ZG long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZGVarID, "cartesian_axis", "Z"),   &
                 "nc_write_model_atts", "ZG cartesian_axis "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZGVarID, "units", "meters"),  &
                 "nc_write_model_atts", "ZG units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZGVarID, "positive", "up"),  &
                 "nc_write_model_atts", "ZG units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZGVarID, "comment", &
                  "more negative is closer to the center of the earth"),  &
                 "nc_write_model_atts", "ZG comment "//trim(filename))

   ! Depths
   call nc_check(nf90_def_var(ncFileID,name="ZC",xtype=nf90_real,dimids=ZCDimID,varid=ZCVarID), &
                 "nc_write_model_atts", "ZC def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, "long_name", "depth at grid centroids"), &
                 "nc_write_model_atts", "ZC long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, "cartesian_axis", "Z"),   &
                 "nc_write_model_atts", "ZC cartesian_axis "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, "units", "meters"),  &
                 "nc_write_model_atts", "ZC units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, "positive", "down"),  &
                 "nc_write_model_atts", "ZC units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, "comment", &
                  "more negative is closer to the center of the earth"),  &
                 "nc_write_model_atts", "ZC comment "//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncid=ncFileID, name="S", xtype=nf90_real, &
         dimids = (/XCDimID,YCDimID,ZCDimID,MemberDimID,unlimitedDimID/),varid=SVarID),&
         "nc_write_model_atts", "S def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SVarID, "long_name", "salinity"), &
         "nc_write_model_atts", "S long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SVarID, "_FillValue", NF90_FILL_REAL), &
         "nc_write_model_atts", "S fill "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SVarID, "units", "psu"), &
         "nc_write_model_atts", "S units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SVarID, "units_long_name", "practical salinity units"), &
         "nc_write_model_atts", "S units_long_name "//trim(filename))

   call nc_check(nf90_def_var(ncid=ncFileID, name="T", xtype=nf90_real, &
         dimids=(/XCDimID,YCDimID,ZCDimID,MemberDimID,unlimitedDimID/),varid=TVarID),&
         "nc_write_model_atts", "T def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, "long_name", "Temperature"), &
         "nc_write_model_atts", "T long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, "_FillValue", NF90_FILL_REAL), &
         "nc_write_model_atts", "T fill "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, "units", "C"), &
         "nc_write_model_atts", "T units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, "units_long_name", "degrees celsius"), &
         "nc_write_model_atts", "T units_long_name "//trim(filename))

   call nc_check(nf90_def_var(ncid=ncFileID, name="U", xtype=nf90_real, &
         dimids=(/XGDimID,YCDimID,ZCDimID,MemberDimID,unlimitedDimID/),varid=UVarID),&
         "nc_write_model_atts", "U def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, "long_name", "Zonal Velocity"), &
         "nc_write_model_atts", "U long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, "_FillValue", NF90_FILL_REAL), &
         "nc_write_model_atts", "U fill "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, "units", "m/s"), &
         "nc_write_model_atts", "U units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, "units_long_name", "meters per second"), &
         "nc_write_model_atts", "U units_long_name "//trim(filename))

   call nc_check(nf90_def_var(ncid=ncFileID, name="V", xtype=nf90_real, &
         dimids=(/XCDimID,YGDimID,ZCDimID,MemberDimID,unlimitedDimID/),varid=VVarID),&
         "nc_write_model_atts", "V def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, "long_name", "Meridional Velocity"), &
         "nc_write_model_atts", "V long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, "_FillValue", NF90_FILL_REAL), &
         "nc_write_model_atts", "V fill "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, "units", "m/s"), &
         "nc_write_model_atts", "V units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, "units_long_name", "meters per second"), &
         "nc_write_model_atts", "V units_long_name "//trim(filename))

   call nc_check(nf90_def_var(ncid=ncFileID, name="Eta", xtype=nf90_real, &
         dimids=(/XCDimID,YCDimID,MemberDimID,unlimitedDimID/),varid=EtaVarID), &
         "nc_write_model_atts", "Eta def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, EtaVarID, "long_name", "sea surface height"), &
         "nc_write_model_atts", "Eta long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, EtaVarID, "_FillValue", NF90_FILL_REAL), &
         "nc_write_model_atts", "Eta fill "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, EtaVarID, "units", "meters"), &
         "nc_write_model_atts", "Eta units "//trim(filename))

   ! Finished with dimension/variable definitions, must end 'define' mode to fill.

   call nc_check(nf90_enddef(ncfileID), "prognostic enddef "//trim(filename))

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables
   !----------------------------------------------------------------------------

   call nc_check(nf90_put_var(ncFileID, XGVarID, XG ), &
                "nc_write_model_atts", "XG put_var "//trim(filename))
   call nc_check(nf90_put_var(ncFileID, XCVarID, XC ), &
                "nc_write_model_atts", "XC put_var "//trim(filename))
   call nc_check(nf90_put_var(ncFileID, YGVarID, YG ), &
                "nc_write_model_atts", "YG put_var "//trim(filename))
   call nc_check(nf90_put_var(ncFileID, YCVarID, YC ), &
                "nc_write_model_atts", "YC put_var "//trim(filename))
   call nc_check(nf90_put_var(ncFileID, ZGVarID, ZG ), &
                "nc_write_model_atts", "ZG put_var "//trim(filename))
   call nc_check(nf90_put_var(ncFileID, ZCVarID, ZC ), &
                "nc_write_model_atts", "ZC put_var "//trim(filename))

endif

!-------------------------------------------------------------------------------
! Fill the variables we can
!-------------------------------------------------------------------------------

call file_to_text("data", textblock)
call nc_check(nf90_put_var(ncFileID, nmlVarID, textblock ), &
              'nc_write_model_atts', 'put_var nmlVarID')
deallocate(textblock)

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID), "nc_write_model_atts", "atts sync")

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

use typeSizes
use netcdf

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: VarID

real(r4), dimension(Nx,Ny,Nz) :: data_3d
real(r4), dimension(Nx,Ny)    :: data_2d
character(len=128)  :: filename

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

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
              "nc_write_model_vars", "inquire "//trim(filename))

if ( output_state_vector ) then

   call nc_check(NF90_inq_varid(ncFileID, "state", VarID), &
                 "nc_write_model_vars", "state inq_varid "//trim(filename))
   call nc_check(NF90_put_var(ncFileID,VarID,statevec,start=(/1,copyindex,timeindex/)),&
                 "nc_write_model_vars", "state put_var "//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   ! Replace missing values (0.0) with netcdf missing value.
   ! Staggered grid causes some logistical problems.
   ! Hopefully, the conversion between r8 and r4 still preserves 'hard' zeros.
   !----------------------------------------------------------------------------

   call vector_to_prog_var(statevec,S_index,data_3d)
   where (data_3d == 0.0_r4) data_3d = NF90_FILL_REAL
   call nc_check(NF90_inq_varid(ncFileID, "S", VarID), &
                "nc_write_model_vars", "S inq_varid "//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                "nc_write_model_vars", "S put_var "//trim(filename))

   call vector_to_prog_var(statevec,T_index,data_3d)
   where (data_3d == 0.0_r4) data_3d = NF90_FILL_REAL
   call nc_check(NF90_inq_varid(ncFileID, "T", VarID), &
                "nc_write_model_vars", "T inq_varid "//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                "nc_write_model_vars", "T put_var "//trim(filename))

   call vector_to_prog_var(statevec,U_index,data_3d)
   where (data_3d == 0.0_r4) data_3d = NF90_FILL_REAL
   call nc_check(NF90_inq_varid(ncFileID, "U", VarID), &
                "nc_write_model_vars", "U inq_varid "//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                "nc_write_model_vars", "U put_var "//trim(filename))

   call vector_to_prog_var(statevec,V_index,data_3d)
   where (data_3d == 0.0_r4) data_3d = NF90_FILL_REAL
   call nc_check(NF90_inq_varid(ncFileID, "V", VarID), &
                "nc_write_model_vars", "V inq_varid "//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                "nc_write_model_vars", "V put_var "//trim(filename))

   call vector_to_prog_var(statevec,Eta_index,data_2d)
   where (data_2d == 0.0_r4) data_2d = NF90_FILL_REAL
   call nc_check(NF90_inq_varid(ncFileID, "Eta", VarID), &
                "nc_write_model_vars", "Eta inq_varid "//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,data_2d,start=(/1,1,copyindex,timeindex/)),&
                "nc_write_model_vars", "Eta put_var "//trim(filename))

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), "nc_write_model_vars", "sync "//trim(filename))

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

! only perturb the non-zero values.  0 is a flag for missing
! ocean cells (e.g. land or under the sea floor)
do i=1,size(state)
   if (state(i) /= 0.0_r8) &
      pert_state(i) = random_gaussian(random_seq, state(i), &
                                      model_perturbation_amplitude)
enddo


end subroutine pert_model_state




subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! If needed by the model interface, this is the current mean
! for all state vector items across all ensembles. It is up to this
! code to allocate space and save a copy if it is going to be used
! later on.  For now, we are ignoring it.

real(r8), intent(in) :: ens_mean(:)

if ( .not. module_initialized ) call static_init_model

end subroutine ens_mean_for_model




function read_meta(fbase, vartype)
!------------------------------------------------------------------
!
! Reads the meta files associated with each 'snapshot'
! and fills the appropriate parts of the output structure.
!
! I believe  pkg/mdsio/mdsio_write_meta.F writes the .meta files 
!
! The files look something like: 
!
! nDims = [   2 ];
! dimList = [
!   256,    1,  256,
!   225,    1,  225
! ];
! dataprec = [ 'float32' ];
! nrecords = [     1 ];
! timeStepNumber = [          0 ];
!
! USAGE:
! metadata = read_meta('U.0000000024')
! ... or ...
! metadata = read_meta('0000000024','U')

character(len=*),           intent(in) :: fbase
character(len=*), OPTIONAL, intent(in) :: vartype
type(MIT_meta_type)                    :: read_meta

character(len=128) :: filename, charstring
integer :: iunit, io
integer :: i, j, indx, nlines, dim1, dimN
logical :: fexist

if ( .not. module_initialized ) call static_init_model

if (present(vartype)) then
   filename = vartype//'.'//trim(fbase)//'.meta'
else
   filename = trim(fbase)//'.meta'
endif

! Initialize to (mostly) bogus values

read_meta%nDims = 0
read_meta%dimList = (/ Nx, Ny, Nz /)
read_meta%dataprec = 'float32'
read_meta%reclen = 0
read_meta%nrecords = 0
read_meta%timeStepNumber = timestepcount

iunit = get_unit()

! See if the file exists ... typically not the case for a cold start.
! If the file does not exist, we must use the namelist values as parsed.

inquire(file=filename, exist=fexist, iostat=io)
if ((io /= 0) .or. (.not. fexist)) then 
   write(msgstring,*) trim(filename), ' does not exist, using namelist values'
   call error_handler(E_MSG,'model_mod:read_meta',msgstring,source,revision,revdate)
   return
endif

! Get next available unit number and open the file

open(unit=iunit, file=filename, action='read', form='formatted', iostat = io)
if (io /= 0) then
   write(msgstring,*) 'cannot open ', trim(filename), ', using namelist values'
   call error_handler(E_MSG,'model_mod:read_meta',msgstring,source,revision,revdate)
   return
endif

! Read every line looking for the nDims entry
! Count the lines just to make future loops more sensible.
! nDims = [   2 ];

nlines = 0
ReadnDims: do i = 1,1000 
   read(iunit,'(a)', iostat = io)charstring
   if (io /= 0) exit ReadnDims
   nlines = nlines + 1

   indx = index(charstring,'nDims = [')
   if (indx > 0) then
      read(charstring(indx+9:),*,iostat=io)read_meta%nDims
      if (io /= 0 )then
         write(msgstring,*)'unable to parse line ',nlines,' from ', trim(filename)
         call error_handler(E_ERR,'model_mod:read_meta:nDims',msgstring,source,revision,revdate)
      endif
   endif
enddo ReadnDims

if (read_meta%nDims < 1) then
   write(msgstring,*) 'unable to determine nDims from ', trim(filename)
   call error_handler(E_ERR,'model_mod:read_meta',msgstring,source,revision,revdate)
endif

! Read every line looking for the dimList entry
! dimList = [
!   256,    1,  256,
!   225,    1,  225
! ];

rewind(iunit)
ReaddimList: do i = 1,nlines 
   
   read(iunit,'(a)', iostat = io)charstring
   if (io /= 0) then
      write(msgstring,*) 'unable to read line ',i,' of ', trim(filename)
      call error_handler(E_ERR,'model_mod:read_meta:dimList',msgstring,source,revision,revdate)
   endif

   indx = index(charstring,'dimList = [')

   if (indx > 0) then
      do j = 1,read_meta%nDims
         read(iunit,*,iostat=io)read_meta%dimList(j),dim1,dimN
         if (io /= 0) then
            write(msgstring,*)'unable to parse dimList(',j, ') from ', trim(filename)
            call error_handler(E_ERR,'model_mod:read_meta',msgstring,source,revision,revdate)
         endif
      enddo
      exit ReaddimList
   endif
enddo ReaddimList

if (all(read_meta%dimList < 1)) then
   write(msgstring,*) 'unable to determine dimList from ', trim(filename)
   call error_handler(E_ERR,'model_mod:read_meta',msgstring,source,revision,revdate)
endif


! Read every line looking for the dataprec entry
! dataprec = [ 'float32' ];

rewind(iunit)
Readdataprec: do i = 1,nlines 

   read(iunit,'(a)', iostat = io)charstring
   if (io /= 0) then
      write(msgstring,*) 'unable to read line ',i,' of ', trim(filename)
      call error_handler(E_ERR,'model_mod:read_meta:dataprec',msgstring,source,revision,revdate)
   endif

   indx = index(charstring,'dataprec = [')

   if (indx > 0) then
      read(charstring(indx+12:),*,iostat=io)read_meta%dataprec
      if (io /= 0) then
         write(msgstring,*)'unable to parse dataprec from ', trim(filename)
         call error_handler(E_ERR,'model_mod:read_meta',msgstring,source,revision,revdate)
      endif
      exit Readdataprec
   endif
enddo Readdataprec

if (index(read_meta%dataprec,'null') > 0) then
   write(msgstring,*) 'unable to determine dataprec from ', trim(filename)
   call error_handler(E_ERR,'model_mod:read_meta',msgstring,source,revision,revdate)
endif


! Read every line looking for the nrecords entry
! nrecords = [     1 ];

rewind(iunit) 
Readnrecords: do i = 1,nlines 
   read(iunit,'(a)', iostat = io)charstring
   if (io /= 0) then
      call error_handler(E_ERR,'model_mod:read_meta','message',source,revision,revdate)
   endif

   indx = index(charstring,'nrecords = [')

   if (indx > 0) then
      read(charstring(indx+12:),*,iostat=io)read_meta%nrecords
      if (io /= 0) then
         write(msgstring,*)'unable to parse nrecords from ', trim(filename)
         call error_handler(E_ERR,'model_mod:read_meta',msgstring,source,revision,revdate)
      endif
      exit Readnrecords
   endif
enddo Readnrecords

if (read_meta%nrecords < 1) then
   write(msgstring,*) 'unable to determine nrecords from ', trim(filename)
   call error_handler(E_ERR,'model_mod:read_meta',msgstring,source,revision,revdate)
endif


! Read every line looking for the timeStepNumber entry
! timeStepNumber = [          0 ];

rewind(iunit)
ReadtimeStepNumber: do i = 1,nlines 
   read(iunit,'(a)', iostat = io)charstring
   if (io /= 0) then
      call error_handler(E_ERR,'model_mod:read_meta','message',source,revision,revdate)
   endif

   indx = index(charstring,'timeStepNumber = [')
   if (indx > 0) then
      read(charstring(indx+18:),*,iostat=io)read_meta%timeStepNumber
      if (io /= 0) then
         write(msgstring,*)'unable to parse timeStepNumber from ', trim(filename)
         call error_handler(E_ERR,'model_mod:read_meta',msgstring,source,revision,revdate)
      endif
      exit ReadtimeStepNumber

   endif
enddo ReadtimeStepNumber

if (read_meta%timeStepNumber < 0) then
   write(msgstring,*) 'unable to determine timeStepNumber from ', trim(filename)
   call error_handler(E_MSG,'model_mod:read_meta',msgstring,source,revision,revdate)
endif

close(iunit)

end function read_meta


subroutine write_meta(metadata, filebase)
!------------------------------------------------------------------
!
type(MIT_meta_type), intent(in) :: metadata
character(len=*), intent(in) :: filebase

integer :: iunit, io, i
character(len=128) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(filebase)//'.meta'

iunit = get_unit()
open(unit=iunit, file=filename, action='write', form='formatted', iostat = io)
if (io /= 0) then
   write(msgstring,*) 'unable to open file ', trim(filename), ' for writing'
   call error_handler(E_ERR,'model_mod:write_meta',msgstring,source,revision,revdate)
endif

write(iunit, "(A,I5,A)") "nDims = [ ", metadata%nDims, " ];"
write(iunit, "(A)")     "dimList = [ "
do i=1, metadata%nDims-1
  write(iunit, "(3(I5,A))") metadata%dimList(i), ',', &
                           1, ',', &
                           metadata%dimList(i), ','
enddo
write(iunit, "(3(I5,A))") metadata%dimList(i), ',', &
                         1, ',', &
                         metadata%dimList(i), ' '

write(iunit, "(A)") "];"

write(iunit, "(3A)") "dataprec = [ '", trim(metadata%dataprec), "' ];"

write(iunit, "(A,I5,A)") "nrecords = [ ", metadata%nrecords, " ];"

write(iunit, "(A,I8,A)") "timeStepNumber = [ ", metadata%timeStepNumber, " ];"

close(iunit)

end subroutine write_meta



subroutine read_2d_snapshot(fbase, x, timestep, vartype)
!------------------------------------------------------------------
!
! This routine reads the Fortran direct access binary files ... eg
! U.nnnnnnnnnn.data    by getting the dimension information from
! U.nnnnnnnnnn.meta
!
! As it stands now ... the .data files appear to be big-endian.

character(len=*),           intent(in)  :: fbase
real(r4), dimension(:,:),   intent(out) :: x
integer,                    intent(out) :: timestep
character(len=*), optional, intent(in)  :: vartype

character(len=128) :: metafilename, datafilename
type(MIT_meta_type):: metadata
integer :: iunit, io
integer :: reclen

if ( .not. module_initialized ) call static_init_model

if (present(vartype)) then
   metafilename = vartype//'.'//trim(fbase)//'.meta'
   datafilename = vartype//'.'//trim(fbase)//'.data'
else
   metafilename = trim(fbase)//'.meta'
   datafilename = trim(fbase)//'.data'
endif

! If the companion ".meta" file exists, it is used. 
! Otherwise, the namelist variables are used.

metadata = read_meta(fbase,vartype)
timestep = metadata%timeStepNumber

! check to make sure storage modes match
! somehow have to pair the string with the fortran kind
! and protect against the redefinition of 'r4' ...
! This code is not robust as it stands ...

if     ( index(metadata%dataprec,'float32') > 0 ) then
   ! r4 is probably OK
else if( index(metadata%dataprec,'real*4') > 0 ) then
   ! r4 is probably OK
else
   write(msgstring,*) 'storage mode mismatch for ', trim(datafilename)
   call error_handler(E_ERR,'model_mod:read_2d_snapshot',msgstring,source,revision,revdate)
endif

! Check dimensions

if (size(x, 1) /= metadata%dimList(1) ) then
   write(msgstring,*)trim(metafilename),' dim 1 does not match delX grid size from namelist'
   call error_handler(E_MSG,"model_mod:read_2d_snapshot",msgstring,source,revision,revdate)
   write(msgstring,*)'expected ',size(x,1),' got ',metadata%dimList(1)
   call error_handler(E_ERR,"model_mod:read_2d_snapshot",msgstring,source,revision,revdate)
endif

if (size(x, 2) /= metadata%dimList(2)) then
   write(msgstring,*)trim(metafilename),' dim 2 does not match delY grid size from namelist'
   call error_handler(E_MSG,"model_mod:read_2d_snapshot",msgstring,source,revision,revdate)
   write(msgstring,*)'expected ',size(x,2),' got ',metadata%dimList(2)
   call error_handler(E_ERR,"model_mod:read_2d_snapshot",msgstring,source,revision,revdate)
endif

reclen = product(shape(x)) * item_size_direct_access

if (do_output()) write(logfileunit,*)'item_size is ',item_size_direct_access, ' reclen is ',reclen
if (do_output()) write(     *     ,*)'item_size is ',item_size_direct_access, ' reclen is ',reclen

! Get next available unit number, read file.

iunit = get_unit()
open(unit=iunit, file=datafilename, action='read', access='direct', recl=reclen, iostat=io)
if (io /= 0) then
   write(msgstring,*) 'cannot open (',io,') file ', trim(datafilename),' for reading.'
   call error_handler(E_ERR,'model_mod:read_2d_snapshot',msgstring,source,revision,revdate)
endif

read(iunit, rec=1, iostat = io) x
if (io /= 0) then
   write(msgstring,*) 'unable to read (',io,') snapshot file ', trim(datafilename)
   call error_handler(E_ERR,'model_mod:read_2d_snapshot',msgstring,source,revision,revdate)
endif

close(iunit)


end subroutine read_2d_snapshot



subroutine read_3d_snapshot(fbase, x, timestep, vartype)
!------------------------------------------------------------------
!
! This routine reads the Fortran direct access binary files ... eg
! U.nnnnnnnnnn.data    by getting the dimension information from
! U.nnnnnnnnnn.meta
!
! As it stands now ... the .data files appear to be big-endian.

character(len=*),           intent(in)  :: fbase
real(r4), dimension(:,:,:), intent(out) :: x
integer,                    intent(out) :: timestep
character(len=*), optional, intent(in)  :: vartype

character(len=128) :: metafilename, datafilename
type(MIT_meta_type):: metadata
integer :: iunit, io
integer :: reclen

if ( .not. module_initialized ) call static_init_model

if (present(vartype)) then
   metafilename = vartype//'.'//trim(fbase)//'.meta'
   datafilename = vartype//'.'//trim(fbase)//'.data'
else
   metafilename = trim(fbase)//'.meta'
   datafilename = trim(fbase)//'.data'
endif

! If the companion ".meta" file exists, it is used. 
! Otherwise, the namelist variables are used.

metadata = read_meta(fbase,vartype)
timestep = metadata%timeStepNumber

! check to make sure storage modes match
! somehow have to pair the string with the fortran kind
! and protect against the redefinition of 'r4' ...
! This code is not robust as it stands ...

if     ( index(metadata%dataprec,'float32') > 0 ) then
   ! r4 is probably OK
else if( index(metadata%dataprec,'real*4') > 0 ) then
   ! r4 is probably OK
else
   write(msgstring,*) 'storage mode mismatch for ', trim(datafilename)
   call error_handler(E_ERR,'model_mod:read_3d_snapshot',msgstring,source,revision,revdate)
endif

! Check dimensions

if (size(x, 1) /= metadata%dimList(1) ) then
   write(msgstring,*)trim(metafilename),' dim 1 does not match delX grid size from namelist'
   call error_handler(E_MSG,"model_mod:read_3d_snapshot",msgstring,source,revision,revdate)
   write(msgstring,*)'expected ',size(x,1),' got ',metadata%dimList(1)
   call error_handler(E_ERR,"model_mod:read_3d_snapshot",msgstring,source,revision,revdate)
endif

if (size(x, 2) /= metadata%dimList(2)) then
   write(msgstring,*)trim(metafilename),' dim 2 does not match delY grid size from namelist'
   call error_handler(E_MSG,"model_mod:read_3d_snapshot",msgstring,source,revision,revdate)
   write(msgstring,*)'expected ',size(x,2),' got ',metadata%dimList(2)
   call error_handler(E_ERR,"model_mod:read_3d_snapshot",msgstring,source,revision,revdate)
endif

if (size(x, 3) /= metadata%dimList(3)) then
   write(msgstring,*)trim(metafilename),' dim 3 does not match delZ grid size from namelist'
   call error_handler(E_MSG,"model_mod:read_3d_snapshot",msgstring,source,revision,revdate)
   write(msgstring,*)'expected ',size(x,3),' got ',metadata%dimList(3)
   call error_handler(E_ERR,"model_mod:read_3d_snapshot",msgstring,source,revision,revdate)
endif

reclen = product(shape(x)) * item_size_direct_access

! Get next available unit number, read file.

iunit = get_unit()
open(unit=iunit, file=datafilename, action='read', access='direct', recl=reclen, iostat=io)
if (io /= 0) then
   write(msgstring,*) 'cannot open file ', trim(datafilename),' for reading.'
   call error_handler(E_ERR,'model_mod:read_3d_snapshot',msgstring,source,revision,revdate)
endif

read(iunit, rec=1, iostat = io) x
if (io /= 0) then
   write(msgstring,*) 'unable to read snapshot file ', trim(datafilename)
   call error_handler(E_ERR,'model_mod:read_3d_snapshot',msgstring,source,revision,revdate)
endif

close(iunit)

end subroutine read_3d_snapshot



subroutine write_2d_snapshot(x, fbase, timestepindex)
!------------------------------------------------------------------
!
! This routine writes the Fortran direct access binary files ... eg
! U.nnnnnnnnnn.data    and
! U.nnnnnnnnnn.meta
!
! As it stands now ... the .data files appear to be big-endian.

real(r4), dimension(:,:), intent(in) :: x
character(len=*),         intent(in) :: fbase
integer, optional,        intent(in) :: timestepindex

character(len=128) :: metafilename, datafilename
type(MIT_meta_type) :: metadata
integer :: iunit, io
integer :: reclen

if ( .not. module_initialized ) call static_init_model

metafilename = trim(fbase)//'.meta'
datafilename = trim(fbase)//'.data'

metadata%nDims = 2
metadata%dimList(:) = (/ Nx, Ny, 0 /)
metadata%dataprec = "float32"
metadata%reclen = Nx * Ny 
metadata%nrecords = 1
if (present(timestepindex)) then
   metadata%timeStepNumber = timestepindex
else
   metadata%timeStepNumber = 0
endif

call write_meta(metadata, fbase)

reclen = metadata%reclen * item_size_direct_access

iunit = get_unit()
open(unit=iunit, file=datafilename, action='write', access='direct', recl=reclen, iostat=io)
if (io /= 0) then
   write(msgstring,*) 'cannot open file ', trim(datafilename),' for reading.'
   call error_handler(E_ERR,'model_mod:write_2d_snapshot',msgstring,source,revision,revdate)
endif

write(iunit, rec=1, iostat = io) x
if (io /= 0) then
   write(msgstring,*) 'unable to read snapshot file ', trim(datafilename)
   call error_handler(E_ERR,'write_2d_snapshot',msgstring,source,revision,revdate)
endif

close(iunit)

end subroutine write_2d_snapshot


subroutine write_3d_snapshot(x, fbase, timestepindex)
!------------------------------------------------------------------
!
! This routine writes the Fortran direct access binary files ... eg
! U.nnnnnnnnnn.data  and the associated
! U.nnnnnnnnnn.meta 
!
! As it stands now ... the .data files appear to be big-endian.

real(r4), dimension(:,:,:), intent(in) :: x
character(len=*),           intent(in) :: fbase
integer, optional,          intent(in) :: timestepindex

character(len=128) :: metafilename, datafilename
type(MIT_meta_type) :: metadata
integer :: iunit, io
integer :: reclen

if ( .not. module_initialized ) call static_init_model

metafilename = trim(fbase)//'.meta'
datafilename = trim(fbase)//'.data'

metadata%nDims = 3
metadata%dimList(:) = (/ Nx, Ny, Nz /)
metadata%dataprec = "float32"     ! FIXME depends on defn of 'r4' 
metadata%reclen = Nx * Ny * Nz
metadata%nrecords = 1
if (present(timestepindex)) then
   metadata%timeStepNumber = timestepindex
else
   metadata%timeStepNumber = 0
endif

call write_meta(metadata, fbase)

reclen = metadata%reclen * item_size_direct_access

! Get next available unit number, write file.

iunit = get_unit()
open(unit=iunit, file=datafilename, action='write', access='direct', recl=reclen, iostat=io)
if (io /= 0) then
   write(msgstring,*) 'cannot open file ', trim(datafilename),' for writing.'
   call error_handler(E_ERR,'model_mod:write_3d_snapshot',msgstring,source,revision,revdate)
endif

write(iunit, rec=1, iostat = io) x
if (io /= 0) then
   write(msgstring,*) 'unable to write snapshot file ', trim(datafilename)
   call error_handler(E_ERR,'write_3d_snapshot',msgstring,source,revision,revdate)
endif

close(iunit)

end subroutine write_3d_snapshot


subroutine snapshot_files_to_sv(timestepcount, state_vector)
!------------------------------------------------------------------
!
integer,  intent(in)    :: timestepcount 
real(r8), intent(inout) :: state_vector(:)

! temp space to hold data while we are reading it
real(r4) :: data_2d_array(Nx,Ny), data_3d_array(Nx,Ny,Nz)
integer :: i, j, k, l, indx, timestepcount_out

! These must be a fixed number and in a fixed order.
character(len=128)  :: prefixstring

if ( .not. module_initialized ) call static_init_model

! start counting up and filling the state vector
! one item at a time, linearizing the 3d arrays into a single
! 1d list of numbers.
indx = 1

! fill S, T, U, V in that order
do l=1, n3dfields

   ! the filenames are constructed here and assumed to be:
   !  Variable.Timestep.[data,.meta]
   ! e.g. S.0000000672.data
   !      S.0000000672.meta
   write(prefixstring, '(A,''.'',I10.10)') trim(progvarnames(l)),timestepcount

   call read_snapshot(prefixstring, data_3d_array, timestepcount_out)
   do k = 1, Nz
   do j = 1, Ny
   do i = 1, Nx
      state_vector(indx) = data_3d_array(i, j, k)
      indx = indx + 1
   enddo
   enddo
   enddo

enddo

! and finally, Eta (and any other 2d fields)
do l=(n3dfields+1), (n3dfields+n2dfields)

   write(prefixstring, '(A,''.'',I10.10)') trim(progvarnames(l)), timestepcount

   call read_snapshot(prefixstring, data_2d_array, timestepcount_out)
   do j = 1, Ny
   do i = 1, Nx
      state_vector(indx) = data_2d_array(i, j)
      indx = indx + 1
   enddo
   enddo

enddo

end subroutine snapshot_files_to_sv



subroutine sv_to_snapshot_files(state_vector, date1, date2)
!------------------------------------------------------------------
!
real(r8), intent(in) :: state_vector(:)
integer,  intent(in) :: date1, date2 

! temp space to hold data while we are writing it
real(r4) :: data_2d_array(Nx,Ny), data_3d_array(Nx,Ny,Nz)
integer :: i, j, k, l, indx

! These must be a fixed number and in a fixed order.
character(len=128)  :: prefixstring

if ( .not. module_initialized ) call static_init_model

! start counting up and filling the state vector
! one item at a time, linearizing the 3d arrays into a single
! 1d list of numbers.
indx = 1

! fill S, T, U, V in that order
do l=1, n3dfields

   ! the filenames are going to be constructed here and assumed to be:
   !  Variable.Basename.Timestep.data  and .meta
   ! e.g. S.Prior.0000000672.data
   !      S.Prior.0000000672.meta
   write(prefixstring, '(A,''.'',I8.8,''.'',I6.6)') trim(progvarnames(l)),date1,date2

   do k = 1, Nz
   do j = 1, Ny
   do i = 1, Nx
      data_3d_array(i, j, k) = state_vector(indx)
      indx = indx + 1
   enddo
   enddo
   enddo
   call write_snapshot(data_3d_array, prefixstring, timestepcount)

enddo

! and finally, Eta (and any other 2d fields)
do l=(n3dfields+1), (n3dfields+n2dfields)

   write(prefixstring, '(A,''.'',I8.8,''.'',I6.6)') trim(progvarnames(l)),date1,date2

   do j = 1, Ny
   do i = 1, Nx
      data_2d_array(i, j) = state_vector(indx)
      indx = indx + 1
   enddo
   enddo
   call write_snapshot(data_2d_array, prefixstring, timestepcount)

enddo

end subroutine sv_to_snapshot_files



subroutine prog_var_to_vector(s,t,u,v,eta,x)
!------------------------------------------------------------------
! deprecated in favor of snapshot_files_to_sv

real(r4), dimension(:,:,:), intent(in)  :: s,t,u,v
real(r4), dimension(:,:),   intent(in)  :: eta
real(r8), dimension(:),     intent(out) :: x

integer :: i,j,k,ii

if ( .not. module_initialized ) call static_init_model

! check shapes

if (size(s,1) /= Nx) then
   write(msgstring,*) 'dim 1 of S /= Nx ',size(s,1),Nx
   call error_handler(E_ERR,'model_mod:prog_var_to_vector', &
                      msgstring,source,revision,revdate) 
endif

if (size(s,2) /= Ny) then
   write(msgstring,*) 'dim 2 of S /= Ny ',size(s,2),Ny
   call error_handler(E_ERR,'model_mod:prog_var_to_vector', &
                      msgstring,source,revision,revdate) 
endif

if (size(s,3) /= Nz) then
   write(msgstring,*) 'dim 3 of S /= Nz ',size(s,3),Nz
   call error_handler(E_ERR,'model_mod:prog_var_to_vector', &
                      msgstring,source,revision,revdate) 
endif

if (size(eta,1) /= Nx) then
   write(msgstring,*) 'dim 1 of Eta /= Nx ',size(eta,1),Nx
   call error_handler(E_ERR,'model_mod:prog_var_to_vector', &
                      msgstring,source,revision,revdate) 
endif

if (size(eta,2) /= Ny) then
   write(msgstring,*) 'dim 2 of Eta /= Ny ',size(eta,2),Ny
   call error_handler(E_ERR,'model_mod:prog_var_to_vector', &
                      msgstring,source,revision,revdate) 
endif

! Should check sizes of T,U,V against that of S

ii = 0

! Salinity
do k = 1,Nz   ! vertical
do j = 1,Ny   ! latitudes
do i = 1,Nx   ! longitudes
   ii = ii + 1
   x(ii) = s(i,j,k)
enddo
enddo
enddo

! Temperature
do k = 1,Nz   ! vertical
do j = 1,Ny   ! latitudes
do i = 1,Nx   ! longitudes
   ii = ii + 1
   x(ii) = t(i,j,k)
enddo
enddo
enddo

! E-W 
do k = 1,Nz   ! vertical
do j = 1,Ny   ! latitudes
do i = 1,Nx   ! longitudes
   ii = ii + 1
   x(ii) = u(i,j,k)
enddo
enddo
enddo

! N-S
do k = 1,Nz   ! vertical
do j = 1,Ny   ! latitudes
do i = 1,Nx   ! longitudes
   ii = ii + 1
   x(ii) = v(i,j,k)
enddo
enddo
enddo

! Sea Surface Height
do j = 1,Ny   ! latitudes
do i = 1,Nx   ! longitudes
   ii = ii + 1
   x(ii) = eta(i,j)
enddo
enddo

if (ii /= get_model_size()) then
   write(msgstring,*)'data size ',ii,' /= ',get_model_size(),' model size'
   call error_handler(E_ERR,'model_mod:prog_var_to_vector', &
                      msgstring,source,revision,revdate) 
endif

end subroutine prog_var_to_vector



subroutine vector_to_2d_prog_var(x, varindex, data_2d_array)
!------------------------------------------------------------------
!
real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: varindex
real(r4), dimension(:,:), intent(out) :: data_2d_array

integer :: i,j,ii
integer :: dim1,dim2
character(len=128) :: varname

if ( .not. module_initialized ) call static_init_model

dim1 = size(data_2d_array,1)
dim2 = size(data_2d_array,2)

varname = progvarnames(varindex)

if (dim1 /= Nx) then
   write(msgstring,*)trim(varname),' 2d array dim 1 ',dim1,' /= ',Nx
   call error_handler(E_ERR,'model_mod:vector_to_2d_prog_var',msgstring,source,revision,revdate) 
endif
if (dim2 /= Ny) then
   write(msgstring,*)trim(varname),' 2d array dim 2 ',dim2,' /= ',Ny
   call error_handler(E_ERR,'model_mod:vector_to_2d_prog_var',msgstring,source,revision,revdate) 
endif

ii = start_index(varindex)

do j = 1,Ny   ! latitudes
do i = 1,Nx   ! longitudes
   data_2d_array(i,j) = x(ii)
   ii = ii + 1
enddo
enddo

end subroutine vector_to_2d_prog_var



subroutine vector_to_3d_prog_var(x, varindex, data_3d_array)
!------------------------------------------------------------------
!
real(r8), dimension(:),     intent(in)  :: x
integer,                    intent(in)  :: varindex
real(r4), dimension(:,:,:), intent(out) :: data_3d_array

integer :: i,j,k,ii
integer :: dim1,dim2,dim3
character(len=128) :: varname

if ( .not. module_initialized ) call static_init_model

dim1 = size(data_3d_array,1)
dim2 = size(data_3d_array,2)
dim3 = size(data_3d_array,3)

varname = progvarnames(varindex)

if (dim1 /= Nx) then
   write(msgstring,*)trim(varname),' 3d array dim 1 ',dim1,' /= ',Nx
   call error_handler(E_ERR,'model_mod:vector_to_3d_prog_var',msgstring,source,revision,revdate) 
endif
if (dim2 /= Ny) then
   write(msgstring,*)trim(varname),' 3d array dim 2 ',dim2,' /= ',Ny
   call error_handler(E_ERR,'model_mod:vector_to_3d_prog_var',msgstring,source,revision,revdate) 
endif
if (dim3 /= Nz) then
   write(msgstring,*)trim(varname),' 3d array dim 3 ',dim3,' /= ',Nz
   call error_handler(E_ERR,'model_mod:vector_to_3d_prog_var',msgstring,source,revision,revdate) 
endif

ii = start_index(varindex)

do k = 1,Nz   ! vertical
do j = 1,Ny   ! latitudes
do i = 1,Nx   ! longitudes
   data_3d_array(i,j,k) = x(ii)
   ii = ii + 1
enddo
enddo
enddo

end subroutine vector_to_3d_prog_var


function timestep_to_DARTtime(TimeStepIndex)
!
! The MITtime is composed of an offset to a fixed time base.
! The base time is derived from the namelist in 'date.cal',
! the model timestep (deltaT) is from the namelist 'PARM03',
! and the timestepindex is the middle portion of the filename
! of the MIT files   [S,T,U,V,Eta].nnnnnnnnnn.dat 
!
! (namelist) startDate_1  yyyymmdd (year/month/day)
! (namelist) startDate_2    hhmmss (hours/minutes/seconds)
! (namelist) deltaTmom   aka 'timestep' ... r4 ... implies roundoff nuances
!
integer, intent(in) :: TimeStepIndex
type(time_type)     :: timestep_to_DARTtime

integer         :: yy,mn,dd,hh,mm,ss
integer         :: modeloffset, maxtimestep
type(time_type) :: offset

if ( .not. module_initialized ) call static_init_model

offset = set_time(0,0)

! Calculate how many seconds/days are represented by the 
! timestepindex and timestepdeltaT .... the offset.
! the timestepindex can be a 10 digit integer ... potential overflow
! when multiplied by a large deltaT

maxtimestep = HUGE(modeloffset)/ocean_dynamics_timestep

if (TimeStepIndex >= maxtimestep) then
   write(msgstring,*)' timestepindex (',TimeStepIndex, &
                     ') * timestep (',ocean_dynamics_timestep,') overflows.'
   call error_handler(E_ERR,'model_mod:timestep_to_DARTtime',msgstring,source,revision,revdate) 
endif

modeloffset = nint(TimeStepIndex * ocean_dynamics_timestep)
dd          = modeloffset /     SECPERDAY    ! use integer arithmetic 
ss          = modeloffset - (dd*SECPERDAY)
offset      = set_time(ss,dd)

! Calculate the DART time_type for the MIT base time.

yy =     startDate_1/10000
mn = mod(startDate_1/100,100)
dd = mod(startDate_1,100)

hh =     startDate_2/10000 
mm = mod(startDate_2/100,100)
ss = mod(startDate_2,100)

timestep_to_DARTtime =  set_date(yy,mn,dd,hh,mm,ss) + offset

end function timestep_to_DARTtime



subroutine DARTtime_to_MITtime(darttime,date1,date2)
!
! The MITtime is composed of:
! (namelist) startDate_1  yyyymmdd (year/month/day)
! (namelist) startDate_2    hhmmss (hours/minutes/seconds)
!
type(time_type), intent(in)  :: darttime
integer,         intent(out) :: date1, date2

integer         :: yy,mn,dd,hh,mm,ss

if ( .not. module_initialized ) call static_init_model

call get_date(darttime,yy,mn,dd,hh,mm,ss)

date1 = yy*10000 + mn*100 + dd
date2 = hh*10000 + mm*100 + ss

!if (do_output()) call print_time(darttime,'DART2MIT dart model time')
!if (do_output()) call print_date(darttime,'DART2MIT dart model date')
!if (do_output()) write(*,*)'DART2MIT',date1,date2

end subroutine DARTtime_to_MITtime



function DARTtime_to_timestepindex(mytime)
!
! The MITtime is composed of an offset to a fixed time base.
! The base time is derived from the namelist in 'date.cal',
! the model timestep (deltaT) is from the namelist 'PARM03',
! and the timestepindex is the middle portion of the filename
! of the MIT files   [S,T,U,V,Eta].nnnnnnnnnn.dat 
!
! (namelist) startDate_1  yyyymmdd (year/month/day)
! (namelist) startDate_2    hhmmss (hours/minutes/seconds)
! (namelist) deltaTmom   aka 'timestep' ... r4 ... implies roundoff nuances
!
type(time_type), intent(in)  :: mytime
integer                      :: DARTtime_to_timestepindex

integer :: dd,ss
type(time_type) :: timeorigin, offset

if ( .not. module_initialized ) call static_init_model

timeorigin = timestep_to_DARTtime(0)
offset     = mytime - timeorigin
call get_time(offset,ss,dd)

if (dd >= (HUGE(dd)/SECPERDAY)) then   ! overflow situation
   call print_time(mytime,'DART time is',logfileunit)
   write(msgstring,*)'Trying to convert DART time to MIT timestep overflows'
   call error_handler(E_ERR,'model_mod:DARTtime_to_timestepindex',msgstring,source,revision,revdate) 
endif

DARTtime_to_timestepindex = nint((dd*SECPERDAY+ss) / ocean_dynamics_timestep)

end function DARTtime_to_timestepindex



subroutine get_gridsize(num_x, num_y, num_z)
!------------------------------------------------------------------
!
 integer, intent(out) :: num_x, num_y, num_z

 num_x = Nx
 num_y = Ny
 num_z = Nz

end subroutine get_gridsize



subroutine write_data_namelistfile
!------------------------------------------------------------------
! Essentially, we want to set the PARM03:endTime value to tell the
! model when to stop. To do that, we have to write an entirely new
! 'data' file, which unfortunately contains multiple namelists.
!
! The strategy here is to determine where the PARM03 namelist starts
! and stops. The stopping part is generally tricky, since 
! the terminator is not well-defined.
!
! So - once we know where the namelist starts and stops, we can
! hunt for the values we need to change and change them while 
! preserving everything else.
! 

integer :: iunit, ounit
integer :: linenum1, linenumE, linenumN
integer :: io, iline
real(r8) :: MyEndTime

character(len=169) :: nml_string, uc_string

if ( .not. file_exist('data') ) then
   call error_handler(E_ERR,'write_data_namelistfile', &
      'namelist file "data" does not exist',source,revision,revdate) 
endif

iunit = open_file('data',      action = 'read')
ounit = open_file('data.DART', action = 'write')
rewind(ounit)

! Find which line number contains &PARM03 and how many lines total.
! Since static_init_model() has already read this once, we know
! that the data file exists and that it contains a PARM03 namelist.
linenumN = 0
linenum1 = 0

FINDSTART : do

   read(iunit, '(A)', iostat = io) nml_string

   if (io < 0 ) then ! end of file
   !  write(*,*)'FINDSTART ... end-of-file at ',linenumN
      exit FINDSTART
   elseif (io /= 0 ) then ! read error
      write(*,msgstring)'manual namelist read failed at line ',linenum1
      call error_handler(E_ERR,'write_data_namelistfile', &
         msgstring,source,revision,revdate) 
   endif

   linenumN = linenumN + 1

   if('&PARM03' == trim(adjustl(nml_string))) then
      linenum1 = linenumN
   endif

enddo FINDSTART

! write(*,*)'Namelist PARM03 starts at line ',linenum1
! write(*,*)'File has ',linenumN,' lines'

if (linenum1 < 1) then
   write(*,msgstring)'unable to find string PARM03'
   call error_handler(E_ERR,'write_data_namelistfile', &
      msgstring,source,revision,revdate) 
endif

! We must preserve the value of the paramters right now,
! so we can write THEM (and not the values we are about to read!)

MyEndTime  = endTime

! Hopefully, we can read the namelist and stay positioned
! Since static_init_model() has already read this once, 
! it is highly unlikely to fail here ...

rewind(iunit) 
read(iunit, nml = PARM03, iostat = io)
if (io /= 0 ) then
   call error_handler(E_ERR,'write_data_namelistfile', &
   'namelist READ failed somehow',source,revision,revdate) 
endif

endTime  = MyEndTime
dumpFreq = MyEndTime
taveFreq = MyEndTime

! Find how many more lines till the end-of-file 
linenumE = 0

FINDEND : do

   read(iunit, '(A)', iostat = io) nml_string

   if (io < 0 ) then ! end of file
   !  write(*,*)'FINDEND ... end-of-file at ',linenumE
      exit FINDEND
   elseif (io /= 0 ) then ! read error
      write(*,msgstring)'manual namelist read failed at line ',linenumE
      call error_handler(E_ERR,'write_data_namelistfile', &
         msgstring,source,revision,revdate) 
   endif

   linenumE = linenumE + 1

enddo FINDEND

! write(*,*)'There are ',linenumE,' lines after the namelist ends.'
rewind(iunit) 

! Read the original namelistfile, write the new namelistfile

do iline = 1,linenum1
    read(iunit, '(A)', iostat = io) nml_string
   write(ounit, '(A)', iostat = io) trim(nml_string)
enddo 
do iline = 1,(linenumN-linenum1-linenumE+1)

    read(iunit, '(A)', iostat = io) nml_string
   uc_string = nml_string
   call to_upper(uc_string)

   if     (index(uc_string,'ADJDUMPFREQ') > 0) then
      continue

   elseif (index(uc_string,'STARTTIME') > 0) then
      write(nml_string,'('' startTime = '',f12.6,'','')')0.0_r8

   elseif (index(uc_string,'DUMPFREQ') > 0) then
      write(nml_string,'('' dumpFreq = '',f12.6,'','')')dumpFreq

   elseif (index(uc_string,'ENDTIME') > 0) then
      write(nml_string,'('' endTime  = '',f12.6,'','')')endTime

   elseif (index(uc_string,'TAVEFREQ') > 0) then
      write(nml_string,'('' taveFreq = '',f12.8,'','')')taveFreq

   endif

   write(ounit, '(A)', iostat = io) trim(nml_string)
enddo 
do iline = 1,(linenumE-1)
    read(iunit, '(A)', iostat = io) nml_string
   write(ounit, '(A)', iostat = io) trim(nml_string)
enddo 

close(iunit)
close(ounit)

end subroutine write_data_namelistfile



!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
