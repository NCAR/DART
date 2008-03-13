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
use        types_mod, only : r8
use time_manager_mod, only : time_type, set_time
use     location_mod, only : location_type,      get_close_maxdist_init, &
                             get_close_obs_init, get_close_obs, set_location, &
                             VERTISHEIGHT
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             logfileunit, &
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
          ens_mean_for_model


! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date: 2007-04-03 16:44:36 -0600 (Tue, 03 Apr 2007) $"


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

! these come from SIZE.h and can vary.  should they be namelist
! controlled?
integer, parameter :: Nx = 1024
integer, parameter :: Ny = 1024
integer, parameter :: Nr = 512

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
      deltaTtracer, dTtracerLev(Nr), deltaTfreesurf, &
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
real(r8) :: dxSpacing, dySpacing, delX(Nx), delY(Ny), &
            phiMin, thetaMin, rSphere, &
            Ro_SeaLevel, delZ, delP, delR(Nr), delRc(Nr+1), &
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

integer, parameter :: S_index   = 1
integer, parameter :: T_index   = 2
integer, parameter :: U_index   = 3
integer, parameter :: V_index   = 4
integer, parameter :: SSH_index = 5

! grid counts for each field
integer :: num_x(nfields), num_y(nfields), num_z(nfields)
integer :: start_index(nfields)

! location information - these grids can either be regularly
! spaced or the spacing along each axis can vary.

real(r8) :: lat_origin, lon_origin
logical :: regular_lat, regular_lon, regular_depth
real(r8) :: delta_lat, delta_lon, delta_depth
real(r8), allocatable :: lat_grid(:), lon_grid(:), depth_grid(:)
real(r8), allocatable :: data_2d_array(:,:), data_3d_array(:,:,:)


! What is the natural model timestep?
real(r8) :: model_timestep
type(time_type) :: time_step
integer :: time_step_days      = 0
integer :: time_step_seconds   = 900


! The real state vector and its length
real(r8), allocatable :: state_vector(:)
integer :: model_size

! Skeleton of a model_nml that would be in input.nml
! This is where dart-related model parms could be set.
logical  :: output_state_vector = .true.
namelist /model_nml/ output_state_vector


contains

!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------
!
! Called to do one time initialization of the model. In this case,
! it reads in the grid information and then the model data.

integer :: iunit, io

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

! depths are going to be in this namelist
call find_namelist_in_file("data", "PARM04", iunit)
read(iunit, nml = PARM04, iostat = io)
call check_namelist_read(iunit, io, "PARM04")

call find_namelist_in_file("data", "PARM05", iunit)
read(iunit, nml = PARM05, iostat = io)
call check_namelist_read(iunit, io, "PARM05")

! for reading/writing the restart files, tim is writing
! some utility routines which will look like this:
!
! call read_snapshot('XC', array(:,:), data_time)
! call read_snapshot('XG', array(:,:), data_time)
! call read_snapshot('YC', array(:,:), data_time)
! call read_snapshot('YG', array(:,:), data_time)
! 
! we pass in a filename (without the .meta and the .data)
! and an unallocated array of the right dimensionality.
! it returns the allocated, filled array, along with the timestamp.

! when we are done, drop_snapshot() frees the space for symmetry


! for now, here's a hard-coded example to get this compiling:
! in spite of the staggering, all grids are the same size
! and offset by half a grid cell.  4 are 3D and 1 is 2D.
!  S,T,U,V = 256 x 225 x 70
!  SSH = 256 x 225

model_size = 4 * (256*225*70) + (256*225)

! Create storage for locations
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

integer :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)

integer :: StateVarVarID   ! netCDF pointer to state variable coordinate array
integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

character(len=129)    :: errstring

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer :: i

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! and then put into define mode.
!-------------------------------------------------------------------------------

ierr = -1 ! assume things go poorly

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, &
                                  nAttributes, unlimitedDimID), "inquire")
call check(nf90_Redef(ncFileID),"redef")

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension. 
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID),"copy dimid")
call check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID),"time dimid")

if ( TimeDimID /= unlimitedDimId ) then
   write(errstring,*)"Time Dimension ID ",TimeDimID, &
                     " should equal Unlimited Dimension ID",unlimitedDimID
   call error_handler(E_ERR,"nc_write_model_atts", errstring, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------
call check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                        len=model_size, dimid = StateVarDimID),"state def_dim")

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date" ,str1    ),"creation put")
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source"  ,source  ),"source put")
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision),"revision put")
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate" ,revdate ),"revdate put")
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","template"       ),"model put")

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
   call check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
              dimids=StateVarDimID, varid=StateVarVarID), "statevariable def_var")
   call check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"), &
                                    "statevariable long_name")
   call check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical"), &
                                    "statevariable units")
   call check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, model_size /)), &
                                    "statevariable valid_range")

   ! Define the actual (3D) state vector, which gets filled as time goes on ... 
   call check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_real, &
              dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), &
              varid=StateVarID), "state def_var")
   call check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"), &
                                           "state long_name")

   ! Leave define mode so we can fill the coordinate variable.
   call check(nf90_enddef(ncfileID),"state enddef")

   ! Fill the state variable coordinate variable
   call check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ), &
                                    "state put_var")

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------

   ! This block is a stub for something more complicated.
   ! Usually, the control for the execution of this block is a namelist variable.
   ! Take a peek at the bgrid model_mod.f90 for a (rather complicated) example.

   call check(nf90_enddef(ncfileID), "prognostic enddef")

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call check(nf90_sync(ncFileID),"atts sync")

ierr = 0 ! If we got here, things went well.

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus, string1)
  !
  ! string1 was added to provide some sense of WHERE things were bombing.
  ! It helps to determine which particular 'nf90_put_att' was generating
  ! the error, for example.

    integer, intent ( in) :: istatus
    character(len=*), intent(in), optional :: string1

    character(len=20)  :: myname = 'nc_write_model_atts '
    character(len=129) :: mystring
    integer            :: indexN

    if( istatus /= nf90_noerr) then

       if (present(string1) ) then
          if ((len_trim(string1)+len(myname)) <= len(mystring) ) then
             mystring = myname // trim(adjustl(string1))
          else
             indexN = len(mystring) - len(myname)
             mystring = myname // string1(1:indexN)
          endif
       else
          mystring = myname
       endif

       call error_handler(E_ERR, mystring, trim(nf90_strerror(istatus)), &
                          source, revision, revdate)
    endif

  end subroutine check

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

integer :: StateVarID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
!-------------------------------------------------------------------------------

ierr = -1 ! assume things go poorly

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, &
                                  nAttributes, unlimitedDimID), "inquire")

if ( output_state_vector ) then

   call check(NF90_inq_varid(ncFileID, "state", StateVarID), "state inq_varid" )
   call check(NF90_put_var(ncFileID, StateVarID, statevec,  &
                start=(/ 1, copyindex, timeindex /)), "state put_var")                   

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------

   ! This block is a stub for something more complicated.
   ! Usually, the control for the execution of this block is a namelist variable.
   ! Take a peek at the bgrid model_mod.f90 for a (rather complicated) example.
   !
   ! Generally, it is necessary to take the statevec and decompose it into 
   ! the separate prognostic variables. In this (commented out) example,
   ! global_Var is a user-defined type that has components like:
   ! global_Var%ps, global_Var%t, ... etc. Each of those can then be passed
   ! directly to the netcdf put_var routine. This may cause a huge storage
   ! hit, so large models may want to avoid the duplication if possible.

   ! call vector_to_prog_var(statevec, get_model_size(), global_Var)

   ! the 'start' array is crucial. In the following example, 'ps' is a 2D
   ! array, and the netCDF variable "ps" is a 4D array [lat,lon,copy,time]

   ! call check(NF90_inq_varid(ncFileID, "ps", psVarID), "ps inq_varid")
   ! call check(nf90_put_var( ncFileID, psVarID, global_Var%ps, &
   !                          start=(/ 1, 1, copyindex, timeindex /) ), "ps put_var")

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call check(nf90_sync(ncFileID), "sync")

ierr = 0 ! If we got here, things went well.

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus, string1)
    integer, intent ( in) :: istatus
    character(len=*), intent(in), optional :: string1

    character(len=20)  :: myname = 'nc_write_model_vars '
    character(len=129) :: mystring
    integer            :: indexN

    if( istatus /= nf90_noerr) then

       if (present(string1) ) then
          if ((len_trim(string1)+len(myname)) <= len(mystring) ) then
             mystring = myname // trim(adjustl(string1))
          else
             indexN = len(mystring) - len(myname)
             mystring = myname // string1(1:indexN)
          endif
       else
          mystring = myname
       endif

       call error_handler(E_ERR, mystring, trim(nf90_strerror(istatus)), &
                          source, revision, revdate)
    endif

  end subroutine check

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



!===================================================================
! End of model_mod
!===================================================================
end module model_mod
