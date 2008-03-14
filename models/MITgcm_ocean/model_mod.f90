! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id: model_mod.f90 2786 2007-04-03 22:44:36Z nancy $
! $Revision$
! $Date: 2007-04-03 16:44:36 -0600 (Tue, 03 Apr 2007) $

! This is the interface between the MITgcm ocean model and DART.

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8
use time_manager_mod, only : time_type, set_time
use     location_mod, only : location_type,      get_close_maxdist_init, &
                             get_close_obs_init, get_close_obs, set_location, &
                             VERTISHEIGHT
use    utilities_mod, only : register_module, error_handler, E_ERR, E_WARN, E_MSG, &
                             logfileunit, get_unit, nc_check, &
                             find_namelist_in_file, check_namelist_read

use     obs_kind_mod, only : KIND_TEMPERATURE

implicit none
private

! FIXME: these belong in obs_kind 
integer, parameter :: KIND_U_CURRENT_COMPONENT = 90
integer, parameter :: KIND_V_CURRENT_COMPONENT = 91
integer, parameter :: KIND_SALINITY = 92
integer, parameter :: KIND_SEA_SURFACE_HEIGHT = 93

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
          ens_mean_for_model,     &
          prog_var_to_vector,     &
          MIT_meta_type, read_meta, read_snapshot, drop_snapshot



! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date: 2007-04-03 16:44:36 -0600 (Tue, 03 Apr 2007) $"

character(len=128) :: msgstring

!------------------------------------------------------------------
!
! MITgcm namelist section:  we want to share the 'data' namelist file
! with the model, so we must declare all possible namelist entries to
! avoid getting an error from a valid namelist file.  Most of these
! values are unused in this model_mod code; only a few are needed and
! those are indicated in comments below.
!
! FIXME: these namelists should probably be in a separate file, and only
! the actual values needed should be made public, so this isn't so messy.
!
!------------------------------------------------------------------

! must match the value in EEPARAMS.h
integer, parameter :: MAX_LEN_FNAM = 512

! these come from SIZE.h and can vary.  
! should they be namelist controlled?
integer, parameter :: max_nx = 1024
integer, parameter :: max_ny = 1024
integer, parameter :: max_nz = 512
integer, parameter :: max_nr = 512

! must match lists declared in ini_parms.f

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
      usingSphericalPolarGrid, phiMin, thetaMin, rSphere, &
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

!-- FIXME: i have not been able to find all of these in the
!-- original source code.  the ones we're using have been
!-- checked that they are the right type.  the others might
!-- cause problems since i don't have a namelist file i can
!-- examine that has the full set of values specified.

!--   Time stepping parameters variable declarations
real(r8) :: nIter0, nTimeSteps, nEndIter, pickupSuff, &
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
      periodicExternalForcing, externForcingPeriod, externForcingCycle, &
      calendarDumps
integer :: momForcingOutAB, tracForcingOutAB
logical :: forcing_In_AB, &
      momDissip_In_AB, doAB_onGtGs, &
      startFromPickupAB2

!--   Gridding parameters variable declarations 
logical :: usingCartesianGrid, usingCylindricalGrid, &
           usingSphericalPolarGrid, usingCurvilinearGrid, &
           deepAtmosphere
real(r8) :: dxSpacing, dySpacing, delX(max_nx), delY(max_ny), &
            phiMin, thetaMin, rSphere, &
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


!------------------------------------------------------------------
!
! The DART state vector (control vector) will consist of:  S, T, U, V, SSH
! (Salinity, Temperature, U velocity, V velocity, Sea Surface Height).
! S, T are 3D arrays, located at cell centers.  U is staggered in X
! and V is staggered in Y (meaning the points are located on the cell
! faces) but the grids are offset by half a cell, so there are actually
! the same number of points in each grid. 
! SSH is a 2D field (X,Y only).  The Z direction is downward.
!
!------------------------------------------------------------------

! Grid parameters that we define - the values will be read from a
! standard MITgcm namelist and filled in here.

integer, parameter :: nfields   = 5
integer, parameter :: n3dfields = 4
integer, parameter :: n2dfields = 1

integer, parameter :: S_index   = 1
integer, parameter :: T_index   = 2
integer, parameter :: U_index   = 3
integer, parameter :: V_index   = 4
integer, parameter :: SSH_index = 5

! FIXME: put the input/restart data filenames into an array
! so we can loop over all 3d arrays and then 2d arrays as we
! fill the state vector.


! grid counts for each field
integer :: Nx, Ny, Nz
integer :: start_index(nfields)

! locations of the cell centers (C) and cell faces (grid?) (G)
! for each axis.
real(r8), allocatable :: XC(:), XG(:), YC(:), YG(:), ZC(:), ZG(:)

! temp space to hold data while we are reading it
real(r4), allocatable :: data_2d_array(:,:), data_3d_array(:,:,:)

! location information - these grids can either be regularly
! spaced or the spacing along each axis can vary.

!real(r8) :: lat_origin, lon_origin
!logical :: regular_lat, regular_lon, regular_depth
!real(r8) :: delta_lat, delta_lon, delta_depth
!real(r8), allocatable :: lat_grid(:), lon_grid(:), depth_grid(:)


! What is the natural model timestep?
real(r8) :: model_timestep
type(time_type) :: time_step
integer :: timestepcount
integer :: time_step_days      = 0
integer :: time_step_seconds   = 900


! The real state vector and its length
real(r8), allocatable :: state_vector(:)
integer :: model_size

! Skeleton of a model_nml that would be in input.nml
! This is where dart-related model parms could be set.
logical  :: output_state_vector = .true.
namelist /model_nml/ output_state_vector

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

INTERFACE drop_snapshot
      MODULE PROCEDURE drop_3d_snapshot
      MODULE PROCEDURE drop_2d_snapshot
END INTERFACE



contains

!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------
!
! Called to do one time initialization of the model. In this case,
! it reads in the grid information and then the model data.

integer :: iunit, io
integer :: i, j, k, indx

! The Plan:
!
!   read the standard MITgcm namelist file 'data' for the
!   file info, the time step size, and maybe some grid info
!   (e.g. projection type)
!
!   open the individual files, one at a time, and read in
!   the meta files for the array sizes.  add them up as you go
!   to compute the total model size
!
!   open the grid data files to get the actual grid coordinates
!
!   open the S,T,U,V,SSH files to read in the data values
!   into the state vector
!
!   set the index numbers where the field types change
!
!   set the grid location info
!

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the DART namelist for this model
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run
call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
write(logfileunit, nml=model_nml)
write(     *     , nml=model_nml)

! Read in the MITgcm namelists from the 'data' file
call find_namelist_in_file("data", "PARM03", iunit)
read(iunit, nml = PARM03, iostat = io)
call check_namelist_read(iunit, io, "PARM03")

delZ(:) = 0.0_r4

! depths are going to be in this namelist
call find_namelist_in_file("data", "PARM04", iunit)
read(iunit, nml = PARM04, iostat = io)
call check_namelist_read(iunit, io, "PARM04")

call find_namelist_in_file("data", "PARM05", iunit)
read(iunit, nml = PARM05, iostat = io)
call check_namelist_read(iunit, io, "PARM05")

! we pass in a filename (without the .meta and the .data)
! and an unallocated array of the right dimensionality.
! it returns the allocated, filled array, along with the timestep count.
! when we are done, drop_snapshot() frees the space.

! first, X (longitude)
call read_snapshot('XC', data_2d_array, timestepcount)
Nx = size(data_2d_array, 1)
Ny = size(data_2d_array, 2)
do i=1, Nx
  XC(i) = data_2d_array(i, 1)
enddo
call drop_snapshot(data_2d_array)

call read_snapshot('XG', data_2d_array, timestepcount)
if (size(data_2d_array, 1) /= Nx) then
  ! call error handler, files are inconsistent sizes
endif
if (size(data_2d_array, 2) /= Ny) then
  ! call error handler, files are inconsistent sizes
endif
do i=1, Nx
  XG(i) = data_2d_array(i, 1)
enddo
call drop_snapshot(data_2d_array)

! then Y (latitude)
call read_snapshot('YC', data_2d_array, timestepcount)
Nx = size(data_2d_array, 1)
Ny = size(data_2d_array, 2)
do i=1, Ny
  YC(i) = data_2d_array(i, 1)
enddo
call drop_snapshot(data_2d_array)

call read_snapshot('YG', data_2d_array, timestepcount)
if (size(data_2d_array, 1) /= Nx) then
  ! call error handler, files are inconsistent sizes
endif
if (size(data_2d_array, 2) /= Ny) then
  ! call error handler, files are inconsistent sizes
endif
do i=1, Ny
  YG(i) = data_2d_array(i, 1)
enddo
call drop_snapshot(data_2d_array)

! compute ZC and ZG here based on delZ from the namelist


! record where in the state vector the data type changes
! from one type to another, by computing the starting
! index for each block of data.
start_index(S_index)   = 1
start_index(T_index)   = start_index(S_index) + (Nx * Ny * Nz)
start_index(U_index)   = start_index(T_index) + (Nx * Ny * Nz)
start_index(V_index)   = start_index(U_index) + (Nx * Ny * Nz)
start_index(SSH_index) = start_index(V_index) + (Nx * Ny * Nz)


! in spite of the staggering, all grids are the same size
! and offset by half a grid cell.  4 are 3D and 1 is 2D.
!  e.g. S,T,U,V = 256 x 225 x 70
!  e.g. SSH = 256 x 225

model_size = (n3dfields * (Nx * Ny * Nz)) + (n2dfields * (Nx * Ny))

! Create storage for model data
allocate(state_vector(model_size))

! The time_step in terms of a time type must also be initialized.
time_step = set_time(time_step_seconds, time_step_days)

end subroutine static_init_model




subroutine init_conditions(x)
!------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

real(r8), intent(out) :: x(:)


end subroutine init_conditions



subroutine adv_1step(x, time)
!------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to 
! compute a timestep, for instance for radiation computations.
! This interface is only called if the namelist parameter
! async is set to 0 in perfect_model_obs of filter or if the 
! program integrate_model is to be used to advance the model
! state as a separate executable. If one of these options
! is not going to be used (the model will only be advanced as
! a separate model-specific executable), this can be a 
! NULL INTERFACE.

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

end subroutine adv_1step



function get_model_size()
!------------------------------------------------------------------
!
! Returns the size of the model as an integer. Required for all
! applications.

integer :: get_model_size

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

! for now, just set to 0
time = set_time(0,0)

end subroutine init_time



subroutine model_interpolate(x, location, itype, obs_val, istatus)
!------------------------------------------------------------------
!
! Given a state vector, a location, and a model state variable type,
! interpolates the state variable field to that location and returns
! the value in obs_val. The istatus variable should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is a model specific integer that specifies the type of field (for
! instance temperature, zonal wind component, etc.). In low order
! models that have no notion of types of variables, this argument can
! be ignored. For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

! FIXME: we need to put interp code here.

! Default for successful return
istatus = 0

end subroutine model_interpolate



function get_model_time_step()
!------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

get_model_time_step = time_step

end function get_model_time_step



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

if (index_in < start_index(S_index+1)) then
   var_type = KIND_SALINITY  
   !location = set_location(lon, lat, depth, VERTISHEIGHT)
else if (index_in < start_index(T_index+1)) then
   var_type = KIND_TEMPERATURE  
   !location =
else if (index_in < start_index(U_index+1)) then
   var_type = KIND_U_CURRENT_COMPONENT
   !location =
else if (index_in < start_index(V_index+1)) then
   var_type = KIND_V_CURRENT_COMPONENT
   !location =
else 
   var_type = KIND_SEA_SURFACE_HEIGHT
   !location =
endif

! something bad happened

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

! good style ... perhaps you could deallocate stuff (from static_init_model?).
! deallocate(state_loc)

end subroutine end_model



function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! TJH 24 Oct 2006 -- Writes the model-specific attributes to a netCDF file.
!     This includes coordinate variables and some metadata, but NOT
!     the model state vector. We do have to allocate SPACE for the model
!     state vector, but that variable gets filled as the model advances.
!
! As it stands, this routine will work for ANY model, with no modification.
!
! The simplest possible netCDF file would contain a 3D field
! containing the state of 'all' the ensemble members. This requires
! three coordinate variables -- one for each of the dimensions 
! [model_size, ensemble_member, time]. A little metadata is useful, 
! so we can also create some 'global' attributes. 
! This is what is implemented here.
!
! Once the simplest case is working, this routine (and nc_write_model_vars)
! can be extended to create a more logical partitioning of the state vector,
! fundamentally creating a netCDF file with variables that are easily 
! plotted. The bgrid model_mod is perhaps a good one to view, keeping
! in mind it is complicated by the fact it has two coordinate systems. 
! There are stubs in this template, but they are only stubs.
!
! TJH 29 Jul 2003 -- for the moment, all errors are fatal, so the
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
integer :: XGDimID, XCDimID, YGDimID, YCDimID, ZGDimID
integer :: XGVarID, XCVarID, YGVarID, YCVarID, ZGVarID

! for the prognostic variables
integer :: SVarID, TVarID, UVarID, VVarID, SSHVarID 

!----------------------------------------------------------------------
! local variables 
!----------------------------------------------------------------------

character(len=129)    :: errstring

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer :: i
character(len=128)  :: filename

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
                                   "nc_write_model_atts", "inquire "//trim(filename))
call nc_check(nf90_Redef(ncFileID),"nc_write_model_atts",   "redef "//trim(filename))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension. 
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID), &
                           "nc_write_model_atts", "copy dimid "//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID), &
                           "nc_write_model_atts", "time dimid "//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(errstring,*)"Time Dimension ID ",TimeDimID, &
             " should equal Unlimited Dimension ID",unlimitedDimID
   call error_handler(E_ERR,"nc_write_model_atts", errstring, source, revision, revdate)
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

   ! S,T,V,SSH Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name="XC",xtype=nf90_real,dimids=XCDimID,varid=XCVarID),&
                 "nc_write_model_atts", "XC def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, XCVarID, "long_name", "longitude grid edges"), &
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

   ! S,T,U,SSH Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name="YC",xtype=nf90_real,dimids=YCDimID,varid=YCVarID), &
                 "nc_write_model_atts", "YC def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YCVarID, "long_name", "latitude grid edges"), &
                 "nc_write_model_atts", "YC long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YCVarID, "cartesian_axis", "Y"),   &
                 "nc_write_model_atts", "YC cartesian_axis "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YCVarID, "units", "degrees_north"),  &
                 "nc_write_model_atts", "YC units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, YCVarID, "valid_range", (/ -90.0_r8, 90.0_r8 /)), &
                 "nc_write_model_atts", "YC valid_range "//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncid=ncFileID, name="S", xtype=nf90_real, &
         dimids = (/XCDimID,YCDimID,ZGDimID,MemberDimID,unlimitedDimID/),varid=SVarID),&
         "nc_write_model_atts", "S def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SVarID, "long_name", "salinity"), &
         "nc_write_model_atts", "S long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SVarID, "units", "psu"), &
         "nc_write_model_atts", "S units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SVarID, "units_long_name", "practical salinity units"), &
         "nc_write_model_atts", "S units_long_name "//trim(filename))

   call nc_check(nf90_def_var(ncid=ncFileID, name="T", xtype=nf90_real, &
         dimids=(/XCDimID,YCDimID,ZGDimID,MemberDimID,unlimitedDimID/),varid=TVarID),&
         "nc_write_model_atts", "T def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, "long_name", "Temperature"), &
         "nc_write_model_atts", "T long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, "units", "C"), &
         "nc_write_model_atts", "T units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, "units_long_name", "degrees celsius"), &
         "nc_write_model_atts", "T units_long_name "//trim(filename))

   call nc_check(nf90_def_var(ncid=ncFileID, name="U", xtype=nf90_real, &
         dimids=(/XGDimID,YCDimID,ZGDimID,MemberDimID,unlimitedDimID/),varid=UVarID),&
         "nc_write_model_atts", "U def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, "long_name", "Zonal Velocity"), &
         "nc_write_model_atts", "U long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, "units", "m/s"), &
         "nc_write_model_atts", "U units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, "units_long_name", "meters per second"), &
         "nc_write_model_atts", "U units_long_name "//trim(filename))

   call nc_check(nf90_def_var(ncid=ncFileID, name="V", xtype=nf90_real, &
         dimids=(/XCDimID,YGDimID,ZGDimID,MemberDimID,unlimitedDimID/),varid=VVarID),&
         "nc_write_model_atts", "V def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, "long_name", "Meridional Velocity"), &
         "nc_write_model_atts", "V long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, "units", "m/s"), &
         "nc_write_model_atts", "V units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, "units_long_name", "meters per second"), &
         "nc_write_model_atts", "V units_long_name "//trim(filename))

   call nc_check(nf90_def_var(ncid=ncFileID, name="SSH", xtype=nf90_real, &
         dimids=(/XCDimID,YCDimID,MemberDimID,unlimitedDimID/),varid=SSHVarID), &
         "nc_write_model_atts", "SSH def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SSHVarID, "long_name", "sea surface height"), &
         "nc_write_model_atts", "SSH long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SSHVarID, "units", "meters"), &
         "nc_write_model_atts", "SSH units "//trim(filename))

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

endif

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

real(r4), allocatable, dimension(:,:,:) :: s,t,u,v
real(r4), allocatable, dimension(:,:)   :: ssh
character(len=128)  :: filename

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
              "nc_write_model_vars", "inquire "//trim(filename))

if ( output_state_vector ) then

   call nc_check(NF90_inq_varid(ncFileID, "state", VarID), &
                 "nc_write_model_vars", "state inq_varid "//trim(filename))
   call nc_check(NF90_put_var(ncFileID,VarID,statevec,start=(/1,copyindex,timeindex/)),&
                 "nc_write_model_vars", "state put_var "//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------

   call vector_to_prog_var(statevec,s,t,u,v,ssh) ! arrays allocated internally

   call nc_check(NF90_inq_varid(ncFileID, "S", VarID), &
                "nc_write_model_vars", "S inq_varid "//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,s,start=(/1,1,1,copyindex,timeindex/)),&
                "nc_write_model_vars", "S put_var "//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, "T", VarID), &
                "nc_write_model_vars", "T inq_varid "//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,t,start=(/1,1,1,copyindex,timeindex/)),&
                "nc_write_model_vars", "T put_var "//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, "U", VarID), &
                "nc_write_model_vars", "U inq_varid "//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,u,start=(/1,1,1,copyindex,timeindex/)),&
                "nc_write_model_vars", "U put_var "//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, "V", VarID), &
                "nc_write_model_vars", "V inq_varid "//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,v,start=(/1,1,1,copyindex,timeindex/)),&
                "nc_write_model_vars", "V put_var "//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, "SSH", VarID), &
                "nc_write_model_vars", "SSH inq_varid "//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,ssh,start=(/1,1,copyindex,timeindex/)),&
                "nc_write_model_vars", "SSH put_var "//trim(filename))

   deallocate(s,t,u,v,ssh)

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
! may do so by adding an O(0.1) magnitude perturbation to each
! model state variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

interf_provided = .false.

end subroutine pert_model_state




subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! Not used in low-order models

real(r8), intent(in) :: ens_mean(:)

end subroutine ens_mean_for_model




  function read_meta(fbase, vartype)
!------------------------------------------------------------------
! function read_meta(fbase,vartype)
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

if (present(vartype)) then
   filename = vartype//'.'//trim(fbase)//'.meta'
else
   filename = trim(fbase)//'.meta'
endif

! Initialize to bogus values

read_meta%nDims = 0
read_meta%dimList = (/ 0, 0, 0 /)
read_meta%dataprec = 'null'
read_meta%reclen = 0
read_meta%nrecords = 0
read_meta%timeStepNumber = 0

! Get next available unit number and open the file

iunit = get_unit()
open(unit=iunit, file=filename, action='read', form='formatted', iostat = io)
if (io /= 0) then
   write(msgstring,*) 'unable to open file ', trim(filename), ' for reading'
   call error_handler(E_ERR,'model_mod:read_meta',msgstring,source,revision,revdate)
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

if (read_meta%nrecords < 1) then
   write(msgstring,*) 'unable to determine timeStepNumber from ', trim(filename)
   call error_handler(E_WARN,'model_mod:read_meta',msgstring,source,revision,revdate)
endif

end function read_meta



  subroutine read_2d_snapshot(fbase, x, timestep, vartype)
!------------------------------------------------------------------
! subroutine read_2d_snapshot(fbase, x, timestep, vartype)
!
! This routine reads the Fortran direct access binary files ... eg
! U.nnnnnnnnnn.data    by getting the dimension information from
! U.nnnnnnnnnn.meta
!
! As it stands now ... the .data files appear to be big-endian.

character(len=*),                        intent(in)  :: fbase
real(r4), allocatable, dimension(:,:),   intent(out) :: x
integer,                                 intent(out) :: timestep
character(len=*), optional,              intent(in)  :: vartype

character(len=128) :: metafilename, datafilename
type(MIT_meta_type):: metadata
integer :: iunit, io, indx
integer :: reclen

if (present(vartype)) then
   metafilename = vartype//'.'//trim(fbase)//'.meta'
   datafilename = vartype//'.'//trim(fbase)//'.data'
else
   metafilename = trim(fbase)//'.meta'
   datafilename = trim(fbase)//'.data'
endif

metadata = read_meta(fbase,vartype)

write(*,*)'nDims   ',metadata%nDims
write(*,*)'dimList ',metadata%dimList

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

! good to go 

allocate(x(metadata%dimList(1), metadata%dimList(2)))

write(*,*)'shape   ',shape(x)

reclen = product(shape(x))

! Get next available unit number, read file.

iunit = get_unit()
open(unit=iunit, file=datafilename, action='read', access='direct', recl=reclen, iostat=io)
if (io /= 0) then
   write(msgstring,*) 'cannot open file ', trim(datafilename),' for reading.'
   call error_handler(E_ERR,'model_mod:read_2d_snapshot',msgstring,source,revision,revdate)
endif

read(iunit, rec=1, iostat = io) x
if (io /= 0) then
   write(msgstring,*) 'unable to read snapshot file ', trim(datafilename)
   call error_handler(E_ERR,'read_2d_snapshot',msgstring,source,revision,revdate)
endif

end subroutine read_2d_snapshot



  subroutine read_3d_snapshot(fbase, x, timestep, vartype)
!------------------------------------------------------------------
! subroutine read_3d_snapshot(fbase, x, timestep, vartype)
!
! This routine reads the Fortran direct access binary files ... eg
! U.nnnnnnnnnn.data    by getting the dimension information from
! U.nnnnnnnnnn.meta
!
! As it stands now ... the .data files appear to be big-endian.

character(len=*),                        intent(in)  :: fbase
real(r4), allocatable, dimension(:,:,:), intent(out) :: x
integer,                                 intent(out) :: timestep
character(len=*), optional,              intent(in)  :: vartype

character(len=128) :: metafilename, datafilename
type(MIT_meta_type):: metadata
integer :: iunit, io, indx
integer :: reclen

if (present(vartype)) then
   metafilename = vartype//'.'//trim(fbase)//'.meta'
   datafilename = vartype//'.'//trim(fbase)//'.data'
else
   metafilename = trim(fbase)//'.meta'
   datafilename = trim(fbase)//'.data'
endif

metadata = read_meta(fbase,vartype)

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

! good to go 

allocate(x(metadata%dimList(1), metadata%dimList(2), metadata%dimList(3)))

reclen = product(shape(x))

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
   call error_handler(E_ERR,'read_3d_snapshot',msgstring,source,revision,revdate)
endif

end subroutine read_3d_snapshot



subroutine drop_3d_snapshot(x)
real(r4), allocatable, intent(inout) :: x(:,:,:)
deallocate(x)
end subroutine drop_3d_snapshot

subroutine drop_2d_snapshot(x)
real(r4), allocatable, intent(inout) :: x(:,:)
deallocate(x)
end subroutine drop_2d_snapshot




subroutine prog_var_to_vector(s,t,u,v,ssh,svector)
real(r4), dimension(:,:,:), intent(in)  :: s,t,u,v
real(r4), dimension(:,:),   intent(in)  :: ssh
real(r8), dimension(:),     intent(out) :: svector

integer :: i,j,k,ii

! check shapes

if (size(s,1) /= Nx) then
   write(msgstring,*),'dim 1 of S /= Nx ',size(s,1),Nx
   call error_handler(E_ERR,'model_mod:prog_var_to_vector', &
                      msgstring,source,revision,revdate) 
endif

if (size(s,2) /= Ny) then
   write(msgstring,*),'dim 2 of S /= Nx ',size(s,2),Nx
   call error_handler(E_ERR,'model_mod:prog_var_to_vector', &
                      msgstring,source,revision,revdate) 
endif

if (size(s,3) /= Nz) then
   write(msgstring,*),'dim 3 of S /= Nx ',size(s,3),Nx
   call error_handler(E_ERR,'model_mod:prog_var_to_vector', &
                      msgstring,source,revision,revdate) 
endif

if (size(ssh,1) /= Nx) then
   write(msgstring,*),'dim 1 of SSH /= Nx ',size(ssh,1),Nx
   call error_handler(E_ERR,'model_mod:prog_var_to_vector', &
                      msgstring,source,revision,revdate) 
endif

if (size(ssh,2) /= Nx) then
   write(msgstring,*),'dim 2 of SSH /= Nx ',size(ssh,2),Nx
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
   svector(ii) = s(i,j,k)
enddo
enddo
enddo

! Temperature
do k = 1,Nz   ! vertical
do j = 1,Ny   ! latitudes
do i = 1,Nx   ! longitudes
   ii = ii + 1
   svector(ii) = t(i,j,k)
enddo
enddo
enddo

! E-W 
do k = 1,Nz   ! vertical
do j = 1,Ny   ! latitudes
do i = 1,Nx   ! longitudes
   ii = ii + 1
   svector(ii) = u(i,j,k)
enddo
enddo
enddo

! N-S
do k = 1,Nz   ! vertical
do j = 1,Ny   ! latitudes
do i = 1,Nx   ! longitudes
   ii = ii + 1
   svector(ii) = v(i,j,k)
enddo
enddo
enddo

! Sea Surface Height
do j = 1,Ny   ! latitudes
do i = 1,Nx   ! longitudes
   ii = ii + 1
   svector(ii) = ssh(i,j)
enddo
enddo

if (ii /= get_model_size()) then
   write(msgstring,*)'data size ',ii,' /= ',get_model_size(),' model size'
   call error_handler(E_ERR,'model_mod:prog_var_to_vector', &
                      msgstring,source,revision,revdate) 
endif

end subroutine prog_var_to_vector



!===================================================================
! End of model_mod
!===================================================================
end module model_mod
