! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module model_mod

! This is the interface between the MITgcm ocean model and DART.

! Modules that are absolutely required for use are listed
use        types_mod,      only : r4, r8, i8, SECPERDAY, vtablenamelength, &
                                  MISSING_I, MISSING_R4, MISSING_R8

use time_manager_mod,      only : time_type, set_time, set_date, get_date, get_time, &
                                  set_calendar_type, GREGORIAN, print_time, print_date, &
                                  operator(*),  operator(+), operator(-), &
                                  operator(>),  operator(<), operator(/), &
                                  operator(/=), operator(<=)

use     location_mod,      only : location_type, get_close_init, &
                                  get_close_state, get_close_obs, set_location, &
                                  VERTISHEIGHT, get_location, is_vertical, &
                                  convert_vertical_obs, convert_vertical_state
! EL use only nc_check was here, deleted for now for testing
use    utilities_mod,      only : error_handler, E_ERR, E_WARN, E_MSG, &
                                  logfileunit, get_unit, do_output, to_upper, &
                                  find_namelist_in_file, check_namelist_read, &
                                  open_file, file_exist, find_textfile_dims, file_to_text, &
                                  string_to_real, string_to_logical

use  netcdf_utilities_mod, only : nc_check

use     obs_kind_mod,      only : QTY_TEMPERATURE, QTY_SALINITY, QTY_U_CURRENT_COMPONENT,  &
                                  QTY_V_CURRENT_COMPONENT, QTY_SEA_SURFACE_HEIGHT,         &
                                  QTY_NITRATE_CONCENTRATION, QTY_SURFACE_CHLOROPHYLL,      &
                                  QTY_PHOSPHATE_CONCENTRATION, QTY_DISSOLVED_OXYGEN,       &
                                  QTY_PHYTOPLANKTON_BIOMASS, QTY_DISSOLVED_INORGANIC_IRON, &
                                  QTY_DISSOLVED_INORGANIC_CARBON, QTY_DISSOLVED_ORGANIC_P, &
                                  QTY_DISSOLVED_ORGANIC_NITROGEN, QTY_ALKALINITY,          &
                                  get_index_for_quantity

use mpi_utilities_mod,     only : my_task_id

use random_seq_mod,        only : random_seq_type, init_random_seq, random_gaussian

use default_model_mod,     only : nc_write_model_vars, adv_1step, &
                                  init_conditions => fail_init_conditions

use dart_time_io_mod,      only : write_model_time

use ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain, get_model_variable_indices, &
                                  get_varid_from_kind, &
                                  get_index_start, get_index_end, &
                                  get_dart_vector_index, get_num_variables, &
                                  get_domain_size, &
                                  get_io_clamping_minval, get_kind_index
                                  
use netcdf_utilities_mod,  only : nc_open_file_readonly, nc_get_variable, & 
                                  nc_get_dimension_size, nc_close_file

use netcdf

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          shortest_time_between_assimilations, &
          end_model,              &
          static_init_model,      &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_copies,      &
          get_close_init,         &
          get_close_obs,          &
          get_close_state,        &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time,        &
          write_model_time

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.
public :: MIT_meta_type, read_meta, write_meta, &
          get_gridsize, &
          write_data_namelistfile, set_model_end_time, &
          timestep_to_DARTtime, DARTtime_to_MITtime, &
          DARTtime_to_timestepindex, &
          lon_bounds,lat_bounds, lon_dist, max_nx, max_ny, max_nz, max_nr, MAX_LEN_FNAM 

character(len=*), parameter :: source = 'MITgcm_ocean/model_mod.f90'

character(len=512) :: string1, string2, string3
logical, save      :: module_initialized = .false.


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

character(len=9) :: TheCalendar         = 'gregorian'
! integration start date follows: yyyymmddhhmmss
integer          :: startDate_1         = 19530101
integer          :: startDate_2         = 60000
logical          :: calendarDumps       = .false.
logical          :: pickupStrictlyMatch = .false.

NAMELIST /CAL_NML/ TheCalendar, startDate_1, startDate_2, calendarDumps

! FIXME: these namelists should probably be in a separate file, and only
! the actual values needed should be made public, so this isn't so messy.

! must match the value in EEPARAMS.h
integer, parameter :: MAX_LEN_FNAM = 512

! these come from SIZE.h and can vary.  
! should they be namelist controlled?
! MEG increased these for the high-res runs
integer, parameter :: max_nx = 2048
integer, parameter :: max_ny = 2048
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
      startFromPickupAB2, &
      checkIniTemp, checkIniSalt

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
      calendarDumps, &
            pickupStrictlyMatch

!--   Gridding parameters namelist
NAMELIST /PARM04/ &
      usingCartesianGrid, usingCylindricalGrid, &
      dxSpacing, dySpacing, delX, delY, delXFile, delYFile, &
      usingSphericalPolarGrid, ygOrigin, xgOrigin, rSphere, &
      usingCurvilinearGrid, horizGridFile, deepAtmosphere, &
      Ro_SeaLevel, delZ, delP, delR, delRc, delRFile, delRcFile, &
      rkFac, groundAtK1

!--   Input files namelist
! MEG: These parameters are not used anywhere in the mod_mod!
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
      the_run_name, &
      checkIniTemp, checkIniSalt

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

integer :: FVAL=-999.0 !SIVA: The FVAL is the fill value used for input netcdf files.

! Grid parameters - the values will be read from a
! standard MITgcm namelist and filled in here.

integer :: Nx=-1, Ny=-1, Nz=-1    ! grid counts for each field
integer :: comp2d = -1, comp3d=-1, comp3dU = -1, comp3dV = -1 ! size of commpressed variables

! locations of cell centers (C) and edges (G) for each axis.
real(r8), allocatable :: XC(:), XG(:), YC(:), YG(:), ZC(:), ZG(:)
real(r4), allocatable :: XC_sq(:), YC_sq(:), XG_sq(:), YG_sq(:)
real(r8), allocatable :: ZC_sq(:)

integer, allocatable :: Xc_Ti(:), Yc_Ti(:), Zc_Ti(:)
integer, allocatable :: Xc_Ui(:), Yc_Ui(:), Zc_Ui(:)
integer, allocatable :: Xc_Vi(:), Yc_Vi(:), Zc_Vi(:) 

real(r8)        :: ocean_dynamics_timestep = 900.0_r4
integer         :: timestepcount = 0
type(time_type) :: model_time, model_timestep

integer(i8) :: model_size    ! the state vector length

! Model namelist declarations with defaults
! Codes for interpreting the columns of the variable_table

integer, parameter :: VT_VARNAMEINDX  = 1 ! ... variable name
integer, parameter :: VT_KINDINDX     = 2 ! ... DART QUANTITY
integer, parameter :: VT_MINVALINDX   = 3 ! ... minimum value if any
integer, parameter :: VT_MAXVALINDX   = 4 ! ... maximum value if any
integer, parameter :: VT_STATEINDX    = 5 ! ... update (state) or not
integer, parameter :: MAX_STATE_VARIABLES = 20
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 5
character(len=vtablenamelength) :: mitgcm_variables(NUM_STATE_TABLE_COLUMNS, MAX_STATE_VARIABLES ) = ' '

character(len=256) :: model_shape_file = ' '
integer  :: assimilation_period_days = 7
integer  :: assimilation_period_seconds = 0
real(r8) :: model_perturbation_amplitude = 0.2

namelist /model_nml/ assimilation_period_days,     &
                     assimilation_period_seconds,  &
                     model_perturbation_amplitude, &
                     model_shape_file,             &
                     mitgcm_variables

logical :: go_to_dart    = .false.
logical :: do_bgc        = .false.
logical :: log_transform = .false.
logical :: compress      = .false.

namelist /trans_mitdart_nml/ go_to_dart, do_bgc, log_transform, compress

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

integer :: domain_id
integer :: nvars

contains

!==================================================================



!------------------------------------------------------------------
!> Called to do one-time initialization of the model.

subroutine static_init_model()

character(len=vtablenamelength) :: var_names(MAX_STATE_VARIABLES) = ' '
integer  :: quantity_list(MAX_STATE_VARIABLES)   = MISSING_I
real(r8) ::    clamp_vals(MAX_STATE_VARIABLES,2) = MISSING_R8
logical  ::   update_list(MAX_STATE_VARIABLES)   = .FALSE.

integer :: i, iunit, io
integer :: ss, dd
integer :: ncid ! for reading compressed coordinates

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

! mit2dart, dart2mit interface
call find_namelist_in_file("input.nml", "trans_mitdart_nml", iunit)
read(iunit, nml = trans_mitdart_nml, iostat = io)
call check_namelist_read(iunit, io, "trans_mitdart_nml")

! MIT calendar information
call find_namelist_in_file("data.cal", "CAL_NML", iunit)
read(iunit, nml = CAL_NML, iostat = io)
call check_namelist_read(iunit, io, "CAL_NML")

if (index(TheCalendar,'g') > 0 ) then
   call set_calendar_type(GREGORIAN)
elseif (index(TheCalendar,'G') > 0 )then
   call set_calendar_type(GREGORIAN)
else
   write(string1,*)"namelist data.cal indicates a ",trim(TheCalendar)," calendar."
   call error_handler(E_MSG,"static_init_model", string1, source)
   write(string1,*)"only have support for Gregorian"
   call error_handler(E_ERR,"static_init_model", string1, source)
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
   write(string1,*)"namelist PARM03 has deltaTmom /= deltaTtracer /= deltaTClock"
   call error_handler(E_MSG,"static_init_model", string1, source)
   write(string1,*)"values were ",deltaTmom, deltaTtracer, deltaTClock
   call error_handler(E_MSG,"static_init_model", string1, source)
   write(string1,*)"At present, DART only supports equal values."
   call error_handler(E_ERR,"static_init_model", string1, source)
endif

! Define the assimilation period as the model_timestep
! Ensure model_timestep is multiple of ocean_dynamics_timestep

model_time     = timestep_to_DARTtime(timestepcount)
model_timestep = set_model_time_step(assimilation_period_seconds, &
                                     assimilation_period_days,    &
                                     ocean_dynamics_timestep)

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)"assimilation period is ",dd," days ",ss," seconds"
call error_handler(E_MSG,'static_init_model',string1,source)
if (do_output()) write(logfileunit,*)string1

! Grid-related variables are in PARM04
delX(:) = 0.0_r4
delY(:) = 0.0_r4
delZ(:) = 0.0_r4
delR(:) = 0.0_r4
call find_namelist_in_file("data", "PARM04", iunit)
read(iunit, nml = PARM04, iostat = io)
call check_namelist_read(iunit, io, "PARM04")

! we use either delR or delZ in mitgcm
if (delR(1) /= 0.0_r4) then
   delZ = delR
endif

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
   write(string1,*)'could not figure out number of longitudes from delX in namelist'
   call error_handler(E_ERR,"static_init_model", string1, source)
endif

Ny = -1
do i=1, size(delY)
 if (delY(i) == 0.0_r4) then
    Ny = i-1
    exit
 endif
enddo
if (Ny == -1) then
   write(string1,*)'could not figure out number of latitudes from delY in namelist'
   call error_handler(E_ERR,"static_init_model", string1, source)
endif

Nz = -1
do i=1, size(delZ)
 if (delZ(i) == 0.0_r4) then
    Nz = i-1
    exit
 endif
enddo
if (Nz == -1) then
   write(string1,*)'could not figure out number of depth levels from delZ in namelist'
   call error_handler(E_ERR,"static_init_model", string1, source)
endif

! We know enough to allocate grid variables. 

if (.not. allocated(XC)) allocate(XC(Nx))
if (.not. allocated(YC)) allocate(YC(Ny))
if (.not. allocated(ZC)) allocate(ZC(Nz))
if (.not. allocated(XG)) allocate(XG(Nx))
if (.not. allocated(YG)) allocate(YG(Ny))
if (.not. allocated(ZG)) allocate(ZG(Nz))

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

! in spite of the staggering, all grids are the same size
! and offset by half a grid cell.  4 are 3D and 1 is 2D.
!  e.g. S,T,U,V = 256 x 225 x 70
!  e.g. Eta = 256 x 225

if (do_output()) write(logfileunit, *) 'Using grid size : '
if (do_output()) write(logfileunit, *) '  Nx, Ny, Nz = ', Nx, Ny, Nz
if (do_output()) write(     *     , *) 'Using grid size : '
if (do_output()) write(     *     , *) '  Nx, Ny, Nz = ', Nx, Ny, Nz

call parse_variable_input(mitgcm_variables, model_shape_file, nvars, &
                      var_names, quantity_list, clamp_vals, update_list)

domain_id = add_domain(model_shape_file, nvars, &
                    var_names, quantity_list, clamp_vals, update_list )

if (compress) then ! read in compressed coordinates

   ncid    = nc_open_file_readonly(model_shape_file)
   comp2d  = nc_get_dimension_size(ncid, 'comp2d' , 'static_init_model', model_shape_file)
   comp3d  = nc_get_dimension_size(ncid, 'comp3d' , 'static_init_model', model_shape_file)
   comp3dU = nc_get_dimension_size(ncid, 'comp3dU', 'static_init_model', model_shape_file)
   comp3dV = nc_get_dimension_size(ncid, 'comp3dV', 'static_init_model', model_shape_file)

   allocate(XC_sq(comp3d))
   allocate(YC_sq(comp3d))
   allocate(ZC_sq(comp3d))  ! ZC is r8

   allocate(XG_sq(comp3d))
   allocate(YG_sq(comp3d))

   allocate(Xc_Ti(comp3d))
   allocate(Yc_Ti(comp3d))
   allocate(Zc_Ti(comp3d))

   allocate(Xc_Ui(comp3dU))
   allocate(Yc_Ui(comp3dU))
   allocate(Zc_Ui(comp3dU))

   allocate(Xc_Vi(comp3dV))
   allocate(Yc_Vi(comp3dV))
   allocate(Zc_Vi(comp3dV))

   call nc_get_variable(ncid, 'XCcomp', XC_sq)
   call nc_get_variable(ncid, 'YCcomp', YC_sq)
   call nc_get_variable(ncid, 'ZCcomp', ZC_sq)

   call nc_get_variable(ncid, 'XGcomp', XG_sq)
   call nc_get_variable(ncid, 'YGcomp', YG_sq)

   call nc_get_variable(ncid, 'Xcomp_ind', Xc_Ti)
   call nc_get_variable(ncid, 'Ycomp_ind', Yc_Ti)
   call nc_get_variable(ncid, 'Zcomp_ind', Zc_Ti)

   call nc_get_variable(ncid, 'Xcomp_indU', Xc_Ui)
   call nc_get_variable(ncid, 'Ycomp_indU', Yc_Ui)
   call nc_get_variable(ncid, 'Zcomp_indU', Zc_Ui)

   call nc_get_variable(ncid, 'Xcomp_indV', Xc_Vi)
   call nc_get_variable(ncid, 'Ycomp_indV', Yc_Vi)
   call nc_get_variable(ncid, 'Zcomp_indV', Zc_Vi)

   call nc_close_file(ncid)

endif

model_size = get_domain_size(domain_id)

if (do_output()) write(*,*) 'model_size = ', model_size

end subroutine static_init_model

function get_model_size()
!------------------------------------------------------------------
!
! Returns the size of the model as an integer. Required for all
! applications.

integer(i8) :: get_model_size

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
!@>todo read time from netCDF input file if possible

time = set_time(0,0)

end subroutine init_time

!-----------------------------------------------------------------------
!

subroutine model_interpolate(state_handle, ens_size, location, quantity, interp_val, istatus)

type(ensemble_type),   intent(in)  :: state_handle
integer,               intent(in)  :: ens_size
type(location_type),   intent(in)  :: location 
integer,               intent(in)  :: quantity
real(r8),              intent(out) :: interp_val(ens_size)
integer,               intent(out) :: istatus(ens_size)

! model_interpolate will interpolate any variable in the DART vector to the given location.
! The first variable matching the quantity of interest will be used for the interpolation.
!
! istatus =  0 ... success
! istatus =  1 ... unknown model level
! istatus =  2 ... vertical coordinate unsupported
! istatus =  3 ... quantity not in the DART vector
! istatus =  4 ... vertical  interpolation failed
! istatus = 11 ... latitude  interpolation failed
! istatus = 12 ... longitude interpolation failed
! istatus = 13 ... corner a retrieval failed
! istatus = 14 ... corner b retrieval failed
! istatus = 15 ... corner c retrieval failed
! istatus = 16 ... corner d retrieval failed

! Local storage
real(r8)    :: loc_array(3), llon, llat, lheight
integer(i8) :: base_offset, offset
integer     :: ind
integer     :: hgt_bot, hgt_top
real(r8)    :: hgt_fract
real(r8)    :: top_val(ens_size), bot_val(ens_size)
integer     :: hstatus
integer     :: i, varid

if ( .not. module_initialized ) call static_init_model

! Successful istatus is 0
interp_val = MISSING_R8
istatus    = 0

! Get the individual locations values
loc_array = get_location(location)
llon      = loc_array(1)
llat      = loc_array(2)
lheight   = loc_array(3)
   
if( is_vertical(location,"HEIGHT") ) then
   ! Nothing to do 
elseif ( is_vertical(location,"SURFACE") ) then
   ! Nothing to do 
elseif (is_vertical(location,"LEVEL")) then
   ! convert the level index to an actual depth 
   ind = nint(loc_array(3))
   if ( (ind < 1) .or. (ind > size(zc)) ) then 
      lheight = zc(ind)
   else
      istatus = 1
      return
   endif
else   ! if pressure or undefined, we don't know what to do
   istatus = 2
   return
endif

! determine which variable is the desired QUANTITY
varid = get_varid_from_kind(domain_id, quantity)

if (varid < 1) then
   istatus = 3
   return
endif

! Do horizontal interpolations for the appropriate levels

! For Sea Surface Height don't need the vertical coordinate
if( is_vertical(location,"SURFACE") ) then
   call lat_lon_interpolate(state_handle, ens_size, llon, llat, 1, varid, quantity, interp_val, istatus)
   return
endif
   
! Get the bounding vertical levels and the fraction between bottom and top
call height_bounds(lheight, nz, zc, hgt_bot, hgt_top, hgt_fract, hstatus)
if(hstatus /= 0) then
   istatus = 4
   return
endif

call lat_lon_interpolate(state_handle, ens_size, llon, llat, hgt_top, varid, quantity, top_val, istatus)
! Failed istatus from interpolate means give up
do i =1,ens_size
   if(istatus(i) /= 0) return
enddo
   
call lat_lon_interpolate(state_handle, ens_size, llon, llat, hgt_bot, varid, quantity, bot_val, istatus)
! Failed istatus from interpolate means give up
do i =1,ens_size
   if(istatus(i) /= 0) return
enddo
! Then weight them by the fraction and return
interp_val = bot_val + hgt_fract * (top_val - bot_val)

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


subroutine lat_lon_interpolate(state_handle, ens_size, llon, llat, level, var_id, qty, interp_val, istatus)
!=======================================================================
!

! Subroutine to interpolate to a lat lon location for a given level

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size

real(r8),            intent(in) :: llon, llat
integer,             intent(in) :: level
integer,             intent(in) :: var_id, qty
integer,            intent(out) :: istatus(ens_size)
real(r8),           intent(out) :: interp_val(ens_size)

! Local storage
real(r8) :: lat_array(ny), lon_array(nx)
integer  :: lat_bot, lat_top, lon_bot, lon_top
real(r8) :: lat_fract, lon_fract
real(r8),dimension(ens_size) :: pa, pb, pc, pd, xbot, xtop
integer  :: lat_status, lon_status
logical  :: masked

if ( .not. module_initialized ) call static_init_model

! Succesful return has istatus of 0
istatus = 0

! Find out what latitude box and fraction
! The latitude grid being used depends on the variable type
! V is on the YG latitude grid

lat_array = yc
if(qty == QTY_V_CURRENT_COMPONENT) lat_array = yg

call lat_bounds(llat, ny, lat_array, lat_bot, lat_top, lat_fract, lat_status)

! Check for error on the latitude interpolation
if(lat_status /= 0) then 
   istatus = 11
   return
endif

! Find out what longitude box and fraction
lon_array = xc
if(qty == QTY_U_CURRENT_COMPONENT) lon_array = xg

call lon_bounds(llon, nx, lon_array, lon_bot, lon_top, lon_fract, lon_status)

! Check for error on the longitude interpolation
if(lon_status /= 0) then 
   istatus = 12
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

pa = get_val(lon_bot, lat_bot, level, var_id, state_handle, ens_size, masked)
if(masked) then
   istatus = 13
   return
endif
pb = get_val(lon_top, lat_bot, level, var_id, state_handle, ens_size, masked)
if(masked) then
   istatus = 14
   return
endif
pc = get_val(lon_bot, lat_top, level, var_id, state_handle, ens_size, masked)
if(masked) then
   istatus = 15
   return
endif
pd = get_val(lon_top, lat_top, level, var_id, state_handle, ens_size, masked)
if(masked) then
   istatus = 16
   return
endif

xbot = pa + lon_fract * (pb - pa)
xtop = pc + lon_fract * (pd - pc)
interp_val = xbot + lat_fract * (xtop - xbot)

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

! This is inefficient, someone could clean it up
! Plus, it doesn't work for a global model that wraps around
do i = 2, nlons
   dist_bot = lon_dist(llon, lon_array(i - 1))
   dist_top = lon_dist(llon, lon_array(i))
   if(dist_bot >= 0 .and. dist_top < 0) then
      bot = i - 1
      top = i
      fract = dist_bot / (dist_bot + abs(dist_top))
      ! orig: fract = abs(dist_bot) / (abs(dist_bot) + dist_top)
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


function get_compressed_dart_vector_index(iloc, jloc, kloc, dom_id, var_id)
!=======================================================================
!

! returns the dart vector index for the compressed state

integer, intent(in) :: iloc, jloc, kloc
integer, intent(in) :: dom_id, var_id
integer(i8)         :: get_compressed_dart_vector_index

integer     :: i    ! loop counter
integer     :: qty
integer(i8) :: offset

offset = get_index_start(dom_id, var_id)

qty = get_kind_index(dom_id, var_id)

get_compressed_dart_vector_index = -1

! MEG: Using the already established compressed indices
!
! 2D compressed variables
if (qty == QTY_SEA_SURFACE_HEIGHT .or. qty == QTY_SURFACE_CHLOROPHYLL ) then
   do i = 1, comp2d
      if (Xc_Ti(i) == iloc .and. Yc_Ti(i) == jloc .and. Zc_Ti(i) == 1) then
         get_compressed_dart_vector_index = offset + i - 1
      endif
   enddo
   return
endif

! 3D compressed variables
if (qty == QTY_U_CURRENT_COMPONENT) then 
   do i = 1, comp3dU
      if (Xc_Ui(i) == iloc .and. Yc_Ui(i) == jloc .and. Zc_Ui(i) == kloc) then 
         get_compressed_dart_vector_index = offset + i - 1
      endif  
   enddo
elseif (qty == QTY_V_CURRENT_COMPONENT) then
   do i = 1, comp3dV
      if (Xc_Vi(i) == iloc .and. Yc_Vi(i) == jloc .and. Zc_Vi(i) == kloc) then 
         get_compressed_dart_vector_index = offset + i - 1
      endif  
   enddo
else
   do i = 1, comp3d
      if (Xc_Ti(i) == iloc .and. Yc_Ti(i) == jloc .and. Zc_Ti(i) == kloc) then 
         get_compressed_dart_vector_index = offset + i - 1
      endif  
   enddo
endif


end function get_compressed_dart_vector_index


function get_val(lon_index, lat_index, level, var_id, state_handle,ens_size, masked)
!=======================================================================
!

! Returns the value from a single level array given the lat and lon indices
integer,             intent(in)  :: lon_index, lat_index, level
integer,             intent(in)  :: var_id ! state variable
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
logical,             intent(out) :: masked
real(r8)                         :: get_val(ens_size)

integer(i8) :: state_index
integer :: i

if ( .not. module_initialized ) call static_init_model

masked = .false.

if (compress) then

   state_index = get_compressed_dart_vector_index(lon_index, lat_index, level, domain_id, var_id)

   if (state_index .ne. -1) then
      get_val = get_state(state_index,state_handle)
   else
      masked = .true.
   endif

else

   state_index = get_dart_vector_index(lon_index, lat_index, level, domain_id, var_id)
   get_val = get_state(state_index,state_handle)

   do i=1,ens_size  ! HK this is checking the whole ensemble, can you have different masks for each ensemble member?
       if(get_val(i) == FVAL) masked = .true.
   enddo

endif

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



function shortest_time_between_assimilations()
!------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = model_timestep

end function shortest_time_between_assimilations



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



!------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location. A second intent(out) optional argument kind
!> can be returned if the model has more than one type of field (for
!> instance temperature and zonal wind component). This interface is
!> required for all filter applications as it is required for computing
!> the distance between observations and state variables.

subroutine get_state_meta_data(index_in, location, qty)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: qty

real(r8) :: lat, lon, depth
integer  :: iloc, jloc, kloc

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, iloc, jloc, kloc, kind_index = qty)

if (compress) then ! all variables ae 1D
   lon   = XC_sq(iloc)
   lat   = YC_sq(iloc)
   depth = ZC_sq(iloc)
   ! Acounting for variables those on staggered grids
   if (qty == QTY_U_CURRENT_COMPONENT) lon   = XG_sq(iloc)
   if (qty == QTY_V_CURRENT_COMPONENT) lat   = YG_sq(iloc)
else

   lon   = XC(iloc)
   lat   = YC(jloc)
   depth = ZC(kloc)

   ! Acounting for variables those on staggered grids
   if (qty == QTY_U_CURRENT_COMPONENT) lon   = XG(iloc)
   if (qty == QTY_V_CURRENT_COMPONENT) lat   = YG(jloc)

endif

! MEG: check chl's depth here
if (qty == QTY_SEA_SURFACE_HEIGHT .or. &
    qty == QTY_SURFACE_CHLOROPHYLL) depth = 0.0_r8

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



!------------------------------------------------------------------
!> Writes the model-specific attributes to a netCDF file.
!> This includes coordinate variables and some metadata,
!> but NOT the model state.

subroutine nc_write_model_atts(ncFileID, domain_id)

! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!

integer, intent(in)  :: ncFileID
integer, intent(in)  :: domain_id

integer :: ierr
integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

!integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)

!----------------------------------------------------------------------
! variables if we parse the state vector into prognostic variables.
!----------------------------------------------------------------------

! for the dimensions and coordinate variables
integer :: XGDimID, XCDimID, YGDimID, YCDimID, ZGDimID, ZCDimID
integer :: XGVarID, XCVarID, YGVarID, YCVarID, ZGVarID, ZCVarID

!----------------------------------------------------------------------
! local variables 
!----------------------------------------------------------------------

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
                                   "nc_write_model_atts", "inquire "//trim(filename))
call nc_check(nf90_Redef(ncFileID),"nc_write_model_atts",   "redef "//trim(filename))


!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source"  ,source  ), &
           "nc_write_model_atts", "source put "//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model",  "MITgcm_ocean" ), &
           "nc_write_model_atts", "model put "//trim(filename))

!-------------------------------------------------------------------------------
! Determine shape of most important namelist
!-------------------------------------------------------------------------------

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

   !>@todo FIXME ... why are we not defining ZGVarID

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
   call nc_check(nf90_put_var(ncFileID, ZCVarID, ZC ), &
                "nc_write_model_atts", "ZC put_var "//trim(filename))


!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID), "nc_write_model_atts", "atts sync")

ierr = 0 ! If we got here, things went well.

end subroutine nc_write_model_atts



!------------------------------------------------------------------
! Create an ensemble of states from a single state.

! Note if you perturb a compressed state, this will not be bitwise
! with perturbing a non-compressed state.
subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: ens_size
real(r8),            intent(in)    :: pert_amp
logical,             intent(out)   :: interf_provided

integer  :: ivar
integer  :: copy
integer  :: i
real(r8) :: pertval, clamp_min_val 
integer  :: iloc, jloc, kloc ! not used, but required for get_model_variable_indices

type(random_seq_type) :: random_seq

if ( .not. module_initialized ) call static_init_model

interf_provided = .true.

call init_random_seq(random_seq, my_task_id())
 
INDICES : do i = 1, state_ens_handle%my_num_vars

   call get_model_variable_indices(state_ens_handle%my_vars(i), iloc, jloc, kloc, ivar)
   clamp_min_val = get_io_clamping_minval(domain_id, ivar)

   MEMBERS : do copy = 1, ens_size

      ! Only perturb the actual ocean cells;
      ! Leave the land and ocean floor values alone.
      if( state_ens_handle%copies(copy, i) == FVAL ) cycle MEMBERS

      pertval = random_gaussian(random_seq, state_ens_handle%copies(copy, i), &
                                      model_perturbation_amplitude)

      ! Clamping: Samples obtained from truncated Gaussian dist. 
      if (.not. log_transform) pertval = max(clamp_min_val, pertval)

      state_ens_handle%copies(copy, i) = pertval

   enddo MEMBERS
enddo INDICES

end subroutine pert_model_copies


!--------------------------------------------------------------------
!> read the time from the input file

function read_model_time(filename)

character(len=*), intent(in) :: filename
type(time_type)              :: read_model_time

if ( .not. module_initialized ) call static_init_model

! The filename is not actually used for this model.
! The time comes from the MITgcm namelists, which are read in static_init_model

read_model_time = model_time

if (do_output()) then
   call print_time(read_model_time, str='MITgcm_ocean time is ',iunit=logfileunit)
   call print_time(read_model_time, str='MITgcm_ocean time is ')
   call print_date(read_model_time, str='MITgcm_ocean date is ',iunit=logfileunit)
   call print_date(read_model_time, str='MITgcm_ocean date is ')
endif

end function read_model_time


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
   write(string1,*) trim(filename), ' does not exist, using namelist values'
   call error_handler(E_MSG,'model_mod:read_meta',string1,source)
   return
endif

! Get next available unit number and open the file

open(unit=iunit, file=filename, action='read', form='formatted', iostat = io)
if (io /= 0) then
   write(string1,*) 'cannot open ', trim(filename), ', using namelist values'
   call error_handler(E_MSG,'model_mod:read_meta',string1,source)
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
         write(string1,*)'unable to parse line ',nlines,' from ', trim(filename)
         call error_handler(E_ERR,'model_mod:read_meta:nDims',string1,source)
      endif
   endif
enddo ReadnDims

if (read_meta%nDims < 1) then
   write(string1,*) 'unable to determine nDims from ', trim(filename)
   call error_handler(E_ERR,'model_mod:read_meta',string1,source)
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
      write(string1,*) 'unable to read line ',i,' of ', trim(filename)
      call error_handler(E_ERR,'model_mod:read_meta:dimList',string1,source)
   endif

   indx = index(charstring,'dimList = [')

   if (indx > 0) then
      do j = 1,read_meta%nDims
         read(iunit,*,iostat=io)read_meta%dimList(j),dim1,dimN
         if (io /= 0) then
            write(string1,*)'unable to parse dimList(',j, ') from ', trim(filename)
            call error_handler(E_ERR,'model_mod:read_meta',string1,source)
         endif
      enddo
      exit ReaddimList
   endif
enddo ReaddimList

if (all(read_meta%dimList < 1)) then
   write(string1,*) 'unable to determine dimList from ', trim(filename)
   call error_handler(E_ERR,'model_mod:read_meta',string1,source)
endif


! Read every line looking for the dataprec entry
! dataprec = [ 'float32' ];

rewind(iunit)
Readdataprec: do i = 1,nlines 

   read(iunit,'(a)', iostat = io)charstring
   if (io /= 0) then
      write(string1,*) 'unable to read line ',i,' of ', trim(filename)
      call error_handler(E_ERR,'model_mod:read_meta:dataprec',string1,source)
   endif

   indx = index(charstring,'dataprec = [')

   if (indx > 0) then
      read(charstring(indx+12:),*,iostat=io)read_meta%dataprec
      if (io /= 0) then
         write(string1,*)'unable to parse dataprec from ', trim(filename)
         call error_handler(E_ERR,'model_mod:read_meta',string1,source)
      endif
      exit Readdataprec
   endif
enddo Readdataprec

if (index(read_meta%dataprec,'null') > 0) then
   write(string1,*) 'unable to determine dataprec from ', trim(filename)
   call error_handler(E_ERR,'model_mod:read_meta',string1,source)
endif


! Read every line looking for the nrecords entry
! nrecords = [     1 ];

rewind(iunit) 
Readnrecords: do i = 1,nlines 
   read(iunit,'(a)', iostat = io)charstring
   if (io /= 0) then
      call error_handler(E_ERR,'model_mod:read_meta','message',source)
   endif

   indx = index(charstring,'nrecords = [')

   if (indx > 0) then
      read(charstring(indx+12:),*,iostat=io)read_meta%nrecords
      if (io /= 0) then
         write(string1,*)'unable to parse nrecords from ', trim(filename)
         call error_handler(E_ERR,'model_mod:read_meta',string1,source)
      endif
      exit Readnrecords
   endif
enddo Readnrecords

if (read_meta%nrecords < 1) then
   write(string1,*) 'unable to determine nrecords from ', trim(filename)
   call error_handler(E_ERR,'model_mod:read_meta',string1,source)
endif


! Read every line looking for the timeStepNumber entry
! timeStepNumber = [          0 ];

rewind(iunit)
ReadtimeStepNumber: do i = 1,nlines 
   read(iunit,'(a)', iostat = io)charstring
   if (io /= 0) then
      call error_handler(E_ERR,'model_mod:read_meta','message',source)
   endif

   indx = index(charstring,'timeStepNumber = [')
   if (indx > 0) then
      read(charstring(indx+18:),*,iostat=io)read_meta%timeStepNumber
      if (io /= 0) then
         write(string1,*)'unable to parse timeStepNumber from ', trim(filename)
         call error_handler(E_ERR,'model_mod:read_meta',string1,source)
      endif
      exit ReadtimeStepNumber

   endif
enddo ReadtimeStepNumber

if (read_meta%timeStepNumber < 0) then
   write(string1,*) 'unable to determine timeStepNumber from ', trim(filename)
   call error_handler(E_MSG,'model_mod:read_meta',string1,source)
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
   write(string1,*) 'unable to open file ', trim(filename), ' for writing'
   call error_handler(E_ERR,'model_mod:write_meta',string1,source)
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


!-------------------------------------------------------------------------------
!> The MITtime is composed of an offset to a fixed time base.
!> The base time is derived from the namelist in 'date.cal',
!> the model timestep (deltaT) is from the namelist 'PARM03',
!> and the timestepindex is the middle portion of the filename
!> of the MIT files   [S,T,U,V,Eta].nnnnnnnnnn.dat 
!
!> (namelist) startDate_1  yyyymmdd (year/month/day)
!> (namelist) startDate_2    hhmmss (hours/minutes/seconds)
!> (namelist) deltaTmom   aka 'timestep' ... r4 ... implies roundoff nuances

function timestep_to_DARTtime(TimeStepIndex)

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
   write(string1,*)' timestepindex (',TimeStepIndex, &
                     ') * timestep (',ocean_dynamics_timestep,') overflows.'
   call error_handler(E_ERR,'model_mod:timestep_to_DARTtime',string1,source) 
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
   write(string1,*)'Trying to convert DART time to MIT timestep overflows'
   call error_handler(E_ERR,'model_mod:DARTtime_to_timestepindex',string1,source) 
endif

DARTtime_to_timestepindex = nint((dd*SECPERDAY+ss) / ocean_dynamics_timestep)

end function DARTtime_to_timestepindex


!------------------------------------------------------------------
!

subroutine get_gridsize(num_x, num_y, num_z)

integer, intent(out) :: num_x, num_y, num_z

num_x = Nx
num_y = Ny
num_z = Nz

end subroutine get_gridsize



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

subroutine write_data_namelistfile

integer :: iunit, ounit
integer :: linenum1, linenumE, linenumN
integer :: io, iline
real(r8) :: MyEndTime

character(len=169) :: nml_string, uc_string

if ( .not. file_exist('data') ) then
   call error_handler(E_ERR,'write_data_namelistfile', &
      'namelist file "data" does not exist',source) 
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
      write(*,string1)'manual namelist read failed at line ',linenum1
      call error_handler(E_ERR,'write_data_namelistfile',string1,source) 
   endif

   linenumN = linenumN + 1

   if('&PARM03' == trim(adjustl(nml_string))) then
      linenum1 = linenumN
   endif

enddo FINDSTART

! write(*,*)'Namelist PARM03 starts at line ',linenum1
! write(*,*)'File has ',linenumN,' lines'

if (linenum1 < 1) then
   write(*,string1)'unable to find string PARM03'
   call error_handler(E_ERR,'write_data_namelistfile',string1,source) 
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
   call error_handler(E_ERR,'write_data_namelistfile','namelist READ failed somehow',source) 
endif

endTime  = MyEndTime
dumpFreq = 10800
taveFreq = 0

! Find how many more lines till the end-of-file 
linenumE = 0

FINDEND : do

   read(iunit, '(A)', iostat = io) nml_string

   if (io < 0 ) then ! end of file
   !  write(*,*)'FINDEND ... end-of-file at ',linenumE
      exit FINDEND
   elseif (io /= 0 ) then ! read error
      write(*,string1)'manual namelist read failed at line ',linenumE
      call error_handler(E_ERR,'write_data_namelistfile',string1,source) 
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


!-----------------------------------------------------------------------
!>
!> Fill the array of requested variables, dart kinds, possible min/max
!> values and whether or not to update the field in the output file.
!>
!>@param state_variables the list of variables and kinds from model_mod_nml
!>@param ngood the number of variable/KIND pairs specified

subroutine parse_variable_input(state_variables, filename, ngood, &
                     var_names, quantity_list, clamp_vals, update_list)

character(len=*), intent(in)  :: state_variables(:,:)
character(len=*), intent(in)  :: filename
integer,          intent(out) :: ngood
character(len=*), intent(out) :: var_names(:)
integer,          intent(out) :: quantity_list(:)
real(r8),         intent(out) :: clamp_vals(:,:)
logical,          intent(out) :: update_list(:)

integer :: i
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr
character(len=NF90_MAX_NAME) :: minvalstring
character(len=NF90_MAX_NAME) :: maxvalstring
character(len=NF90_MAX_NAME) :: updateable

ngood = 0
MyLoop : do i = 1, MAX_STATE_VARIABLES

   varname      = trim(state_variables(VT_VARNAMEINDX,i))
   dartstr      = trim(state_variables(VT_KINDINDX   ,i))
   minvalstring = trim(state_variables(VT_MINVALINDX ,i))
   maxvalstring = trim(state_variables(VT_MAXVALINDX ,i))
   updateable   = trim(state_variables(VT_STATEINDX  ,i))

   if ( varname == ' ' .and. dartstr == ' ' ) exit MyLoop ! Found end of list.

   if ( varname == ' ' .or. dartstr == ' ' ) then
      string1 = 'model_nml:model "variables" not fully specified'
      string2 = 'reading from "'//trim(filename)//'"'
      call error_handler(E_ERR,'parse_variable_input:',string1,source,text2=string2)
   endif

   ! Make sure DART quantity is valid

   if( get_index_for_quantity(dartstr) < 0 ) then
      write(string1,'(''there is no quantity <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'parse_variable_input:',string1,source)
   endif

   call to_upper(minvalstring)
   call to_upper(maxvalstring)
   call to_upper(updateable)

   var_names(i)     = varname
   quantity_list(i) = get_index_for_quantity(dartstr)
   clamp_vals(i, 1) = string_to_real(minvalstring)
   clamp_vals(i, 2) = string_to_real(maxvalstring)
   update_list(i)   = string_to_logical(updateable, 'UPDATE')

   ! Adjust clamping in case of log-transform
   if (quantity_list(i) == QTY_NITRATE_CONCENTRATION     ) call adjust_clamp(clamp_vals(i, 1))
   if (quantity_list(i) == QTY_PHOSPHATE_CONCENTRATION   ) call adjust_clamp(clamp_vals(i, 1))
   if (quantity_list(i) == QTY_DISSOLVED_OXYGEN          ) call adjust_clamp(clamp_vals(i, 1))
   if (quantity_list(i) == QTY_PHYTOPLANKTON_BIOMASS     ) call adjust_clamp(clamp_vals(i, 1))
   if (quantity_list(i) == QTY_ALKALINITY                ) call adjust_clamp(clamp_vals(i, 1))
   if (quantity_list(i) == QTY_DISSOLVED_INORGANIC_CARBON) call adjust_clamp(clamp_vals(i, 1))
   if (quantity_list(i) == QTY_DISSOLVED_ORGANIC_P       ) call adjust_clamp(clamp_vals(i, 1))
   if (quantity_list(i) == QTY_DISSOLVED_ORGANIC_NITROGEN) call adjust_clamp(clamp_vals(i, 1))
   if (quantity_list(i) == QTY_DISSOLVED_INORGANIC_IRON  ) call adjust_clamp(clamp_vals(i, 1))
   if (quantity_list(i) == QTY_SURFACE_CHLOROPHYLL       ) call adjust_clamp(clamp_vals(i, 1))

   ngood = ngood + 1

enddo MyLoop

if (ngood == MAX_STATE_VARIABLES) then
   string1 = 'WARNING: There is a possibility you need to increase ''MAX_STATE_VARIABLES'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'parse_variable_input:',string1,source,text2=string2)
endif

end subroutine parse_variable_input


!-----------------------------------------------------------------------
!>
!> We may choose to do a log transformation for the bgc tracers
!> Make sure the user has the right lower_bound 

subroutine adjust_clamp(lower_bound)

real(r8), intent(inout) :: lower_bound

if (log_transform) then 
   lower_bound = MISSING_R8
else
   lower_bound = 0.0_r8
endif

end subroutine adjust_clamp


!===================================================================
! End of MITgcm_ocean model_mod
!===================================================================
end module model_mod

