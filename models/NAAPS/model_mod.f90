! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This is a template showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. The public interfaces listed
! must all be supported with the argument lists as indicated. Many of the interfaces
! are not required for minimal implementation (see the discussion of each
! interface and look for NULL INTERFACE). 

! Modules that are absolutely required for use are listed

use        types_mod, only : r8, r4, MISSING_R8, metadatalength, obstypelength
use time_manager_mod, only : time_type, set_calendar_type, operator(/=), &
                             set_time, print_time, set_date, print_date
use     location_mod, only : location_type,      get_close_maxdist_init, &
                             get_close_obs_init, get_close_obs, set_location, &
                             set_location_missing, VERTISUNDEF, vert_is_height, &
                             vert_is_level, vert_is_surface, vert_is_undef, &
                             get_location
use    utilities_mod, only : register_module, error_handler, nc_check, &
                             E_ERR, E_MSG, get_unit, open_file, close_file, &
                             find_namelist_in_file, check_namelist_read, &
                             do_nml_file, do_nml_term, do_output, &
                             nmlfileunit, logfileunit
use     obs_kind_mod, only : get_index_for_quantity

use netcdf
use typesizes

implicit none
private

public :: get_model_size,               &
          adv_1step,                    &
          get_state_meta_data,          &
          model_interpolate,            &
          get_model_time_step,          &
          end_model,                    &
          static_init_model,            &
          init_time,                    &
          init_conditions,              &
          nc_write_model_atts,          &
          nc_write_model_vars,          &
          pert_model_state,             &
          get_close_maxdist_init,       &
          get_close_obs_init,           &
          get_close_obs,                &
          ens_mean_for_model

public :: analysis_file_to_statevector, &
          statevector_to_analysis_file, &
          get_naaps_restart_path,       &
          get_naaps_metadata,           &
          get_naaps_dtg,                &
          get_naaps_ensemble_member,    &
          dtg_to_time

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! EXAMPLE: define model parameters here
integer                          :: model_size   !-Length of state vector 
type(time_type)                  :: time_step
integer                          :: io, iunit, istat
integer                          :: nx, ny, nz, nspecies, nvars 
REAL(r4),ALLOCATABLE             :: xlat(:), xlon(:)
REAL(r8)                         :: dlat, dlon
logical, save                    :: module_initialized  = .false.
character(len=256)               :: string1, string2

! model namelist variables and declaration 
character(len=10)  :: dtg                 = '1999123100'
integer            :: nens                = 20
character(len=256) :: naaps_restart_path  = './'
integer            :: time_step_days      = 0
integer            :: time_step_seconds   = 21600
logical            :: output_state_vector = .false.
logical            :: debug               = .false.
integer            :: member              = 1

namelist /model_nml/ nens, dtg, naaps_restart_path, &
                     time_step_seconds, time_step_days, &
                     debug, member

! Everything needed to describe a variable
integer, parameter :: max_state_variables = 80

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME), dimension(NF90_MAX_VAR_DIMS) :: dimname
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer :: xtype         ! netCDF variable type (NF90_double, etc.) 
   integer :: numdims       ! number of dims - excluding TIME
   integer :: numvertical   ! number of vertical levels in variable
   integer :: varsize       ! prod(dimlens(1:numdims))
   integer :: index1        ! location in dart state vector of first occurrence
   integer :: indexN        ! location in dart state vector of last  occurrence
   integer :: dart_kind
   character(len=obstypelength) :: kind_string
   logical  :: clamping     ! does variable need to be range-restricted before 
   real(r8) :: range(2)     ! being stuffed back into MPAS analysis file.
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

real(r8), allocatable, dimension(:) :: module_ensemble_mean

INTERFACE vector_to_prog_var
      MODULE PROCEDURE vector_to_1d_prog_var
      MODULE PROCEDURE vector_to_2d_prog_var
      MODULE PROCEDURE vector_to_3d_prog_var
END INTERFACE

INTERFACE dtg_to_time
      MODULE PROCEDURE dtg_integer_to_time
      MODULE PROCEDURE dtg_string_to_time
END INTERFACE


contains

!==================================================================

subroutine static_init_model()
!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.
! Can be a NULL INTERFACE for the simplest models.

       INTEGER           :: i
       TYPE(time_type)   :: model_time
       CHARACTER(len=metadatalength),dimension(10) :: species_names = (/ &
                     'sulfate_aod ', &
                     'dust_aod    ', &
                     'smoke_aod   ', &
                     'seasalt_aod ', &
                     'total_aod   ', &
                     'so2_conc    ', &
                     'sulfate_conc', &
                     'dust_conc   ', &
                     'smoke_conc  ', &
                     'seasalt_conc'/)
       CHARACTER(len=metadatalength),dimension(10) :: species_kinds = (/ &
                     'QTY_INTEGRATED_SULFATE', &
                     'QTY_INTEGRATED_DUST   ', &
                     'QTY_INTEGRATED_SMOKE  ', &
                     'QTY_INTEGRATED_SEASALT', &
                     'QTY_INTEGRATED_AOD    ', & 
                     'QTY_SO2               ', &
                     'QTY_SULFATE           ', &
                     'QTY_DUST              ', &
                     'QTY_SMOKE             ', &
                     'QTY_SEASALT           '/)
       !integer  :: iunit, io

       if ( module_initialized ) return

       module_initialized = .true.

       ! Print module information to log file and stdout.
       call register_module(source, revision, revdate)

       !_Update the basics - dtg, path, et cetera
       call find_namelist_in_file("input.nml", "model_nml", iunit)
       read(iunit, nml = model_nml, iostat = io)
       call check_namelist_read(iunit, io, "model_nml")
       call set_calendar_type('Gregorian')

       ! Record the namelist values used for the run ...
       if (do_nml_file()) write(nmlfileunit, nml=model_nml)
       if (do_nml_term()) write(     *     , nml=model_nml)

       !_Provided the path, open 

       !_Open sample concentration file, get metadata
       call get_naaps_metadata( naaps_restart_path, dtg, model_time )
       
       !_Definine 2D aod fields.  3D to follow
       DO i=1, nspecies
           progvar(i)%numdims = 2 !_lat/lon fields first
           progvar(i)%dimlens(1) = nx
           progvar(i)%dimlens(2) = ny
           progvar(i)%dimname(1) = 'LONGITUDE'
           progvar(i)%dimname(2) = 'LATITUDE'
           progvar(i)%numvertical = 1
           progvar(i)%varsize = nx * ny
           progvar(i)%index1 = (i-1)*(nx*ny)+1 !
           progvar(i)%indexN = i*(nx*ny) 
           progvar(i)%varname = species_names(i)
           progvar(i)%kind_string = species_kinds(i)
           progvar(i)%dart_kind = get_index_for_quantity(progvar(i)%kind_string)

           if (debug) WRITE(*,*)
           if (debug) WRITE(*,*) progvar(i)%numdims 
           if (debug) WRITE(*,*) progvar(i)%dimlens(1:progvar(i)%numdims)
           if (debug) WRITE(*,*) progvar(i)%numvertical
           if (debug) WRITE(*,*) progvar(i)%varsize 
           if (debug) WRITE(*,*) progvar(i)%index1 
           if (debug) WRITE(*,*) progvar(i)%indexN
           if (debug) WRITE(*,*) trim(progvar(i)%varname) 
           if (debug) WRITE(*,*) trim(progvar(i)%kind_string)
       END DO

       DO i=nspecies+1, nspecies*2
           progvar(i)%numdims = 3 !_lat/lon fields first
           progvar(i)%dimlens(1) = nx
           progvar(i)%dimlens(2) = ny
           progvar(i)%dimlens(3) = nz
           progvar(i)%dimname(1) = 'LONGITUDE'
           progvar(i)%dimname(2) = 'LATITUDE'
           progvar(i)%dimname(3) = 'SIGMA'
           progvar(i)%numvertical = nz 
           progvar(i)%varsize = nx * ny * nz
           progvar(i)%index1 = progvar(i-1)%indexN + 1 
           progvar(i)%indexN = progvar(i)%index1 + nx*ny*nz -1 
           progvar(i)%varname = species_names(i)
           progvar(i)%kind_string = species_kinds(i)
           progvar(i)%dart_kind = get_index_for_quantity(progvar(i)%kind_string)

           if (debug) WRITE(*,*)
           if (debug) WRITE(*,*) progvar(i)%numdims 
           if (debug) WRITE(*,*) progvar(i)%dimlens(1:progvar(i)%numdims)
           if (debug) WRITE(*,*) progvar(i)%numvertical
           if (debug) WRITE(*,*) progvar(i)%varsize 
           if (debug) WRITE(*,*) progvar(i)%index1 
           if (debug) WRITE(*,*) progvar(i)%indexN
           if (debug) WRITE(*,*) trim(progvar(i)%varname) 
           if (debug) WRITE(*,*) trim(progvar(i)%kind_string)
           
           model_size = progvar(i)%indexN
       END DO
       nvars = nspecies * 2

       ! The time_step in terms of a time type must also be initialized.
       time_step = set_time(time_step_seconds, time_step_days)

       ! Allocate space for an ensemble mean
       allocate(module_ensemble_mean(model_size))

end subroutine static_init_model




subroutine init_conditions(x)
!------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If it is not possible to use some default initial conditions, this
! can be a NULL INTERFACE that should issue a horrible warning and DIE.

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'Cannot initialize NAAPS with a default state.'
write(string2,*) 'perfect_model_obs_nml:start_from_restart cannot be .FALSE.'
call error_handler(E_ERR,'init_conditions',string1,source,revision,revdate,&
                                     text2=string2)

x = MISSING_R8 ! tell compiler to be quiet about unused variables

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
! NULL INTERFACE, just die if it ever gets here.

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

if ( .not. module_initialized ) call static_init_model

if (do_output()) then
   call print_time(time,'NULL interface adv_1step (no advance) DART time is')
   call print_time(time,'NULL interface adv_1step (no advance) DART time is',logfileunit)
endif

write(string1,*) 'Cannot advance NAAPS with a subroutine call; async cannot equal 0'
call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate)

x = MISSING_R8 ! tell compiler to be quiet about unused variables

end subroutine adv_1step



function get_model_size()
!------------------------------------------------------------------
! Returns the size of the model as an integer. Required. 

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size



subroutine init_time(time)
!------------------------------------------------------------------
! Companion interface to init_conditions. Returns a time that is somehow 
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If it is not possible to use some default initial conditions, this
! can be a NULL INTERFACE that should issue a horrible warning and DIE.

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'Cannot initialize NAAPS with a default state.'
write(string2,*) 'perfect_model_obs_nml:start_from_restart cannot be .FALSE.'
call error_handler(E_ERR,'init_time',string1,source,revision,revdate,&
                                     text2=string2)

! it might be reasonable to set it to the dtg in the model_nml
time = dtg_to_time(dtg)

end subroutine init_time



subroutine model_interpolate(x, location, itype, interp_val, istatus)
!------------------------------------------------------------------
!
! Given a state vector, a location, and a model state variable type,
! interpolates the state variable field to that location and returns
! the value in interp_val. The istatus variable should be returned as
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
real(r8),           intent(out) :: interp_val
integer,            intent(out) :: istatus

  real(r8)         :: loc_array(3), llon, llat
  real(r8)         :: lheight
  integer          :: base_offset, rel_offset
  integer          :: i
  integer          :: iloc, jloc, kloc

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

  IF ( llon .LT. 180.0_r8 ) llon = llon + 360.0_r8 

! TJH 27.Oct.2011 This model_interpolate routine only works for '2D' observations.
!                 There is no vertical interpolation capability.
!                 The 'interpolation' is trivial - the assumption is that 
!                 the observations have been preoprocessd to be coincident with
!                 the NAAPS state vector.

  IF( vert_is_undef(location) ) THEN  ! Nothing to do 
  ELSE   ! if we don't know what to do, error out
      istatus = 15
      return
  ENDIF

!  ELSEIF ( vert_is_surface(location) ) THEN  ! Nothing to do
!  ELSEIF (vert_is_level(location)) THEN      ! convert the level index to an actual height
!     kloc = nint(loc_array(3))
!     IF( (kloc < 1) .or. (kloc > size(zc)) ) THEN
!        istatus = 33
!        return
!     ELSE
!        lheight = zc(kloc)
!     ENDIF

  kloc = 1
  jloc = NINT((llat - xlat(1)) / dlat) + 1
  iloc = NINT((llon - xlon(1)) / dlon) + 1
  MY_SPECIES: DO i = 1, nvars
    IF ( progvar(i)%dart_kind == itype ) THEN
      base_offset = progvar(i)%index1
      rel_offset = nx*ny*(kloc - 1) + nx * (jloc - 1) + iloc - 1 
      interp_val = x(base_offset + rel_offset)
      istatus = 0 
      EXIT MY_SPECIES
    ENDIF
  END DO MY_SPECIES


  IF ( 1 /= 1 ) THEN
     print *, 'requesting interpolation at '
     print *, llon, iloc, xlon(iloc) 
     print *, llat, jloc, xlat(jloc)
     print *, base_offset, rel_offset, rel_offset + base_offset
  ENDIF

end subroutine model_interpolate



function get_model_time_step()
!------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model
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

integer,             intent(in)            :: index_in
type(location_type), intent(out)           :: location
integer,             intent(out), optional :: var_type

integer :: nxp, nyp, iloc, jloc, kloc, nf, n
integer :: myindx

if (.not. module_initialized ) call static_init_model

myindx = -1
nf     = -1

FindIndex : DO n = 1, nvars
  IF( (progvar(n)%index1 <= index_in) .and. (index_in <= progvar(n)%indexN) ) THEN
    nf = n
    myindx = index_in - progvar(n)%index1 + 1
    EXIT FindIndex
  ENDIF
ENDDO FindIndex

IF( myindx == -1 ) THEN
   write(string1,*) 'Problem, cannot find base_offset, index_in is: ', & 
                    index_in
   call error_handler(E_ERR,'get_state_meta_data',string1,source,revision, &
                     revdate)
ENDIF

nxp = progvar(nf)%dimlens(1)
nyp = progvar(nf)%dimlens(2)

kloc   = 1 + (myindx-1) / (nxp*nyp)
myindx = myindx - (kloc-1)*nyp*nxp
jloc   = 1 + (myindx-1) / nxp
myindx = myindx - (jloc-1)*nxp
iloc   = myindx

location = set_location(real(xlon(iloc),r8), real(xlat(jloc),r8), 0.0_r8, VERTISUNDEF)

IF (present(var_type)) THEN
   var_type = progvar(nf)%dart_kind
ENDIF

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

if (allocated(module_ensemble_mean)) deallocate(module_ensemble_mean)
if (allocated(xlon))                 deallocate(xlon)
if (allocated(xlat))                 deallocate(xlat)

end subroutine end_model



function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! TJH 24 Oct 2011 -- Writes the model-specific attributes to a netCDF file.
!     This includes coordinate variables and some metadata, but NOT
!     the model state vector. We do have to allocate SPACE for the model
!     state vector, but that variable gets filled as the model advances.
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
integer :: lonDimID        ! netCDF pointer to longitude dimension 
integer :: latDimID        ! netCDF pointer to latitude  dimension 
integer :: levDimID        ! netCDF pointer to vertical  dimension 
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)
integer :: VarID

character(len=128) :: filename

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer :: i, ivar, idim, ndims
integer, dimension(NF90_MAX_VAR_DIMS) :: dims

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

call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                                    "nc_write_model_atts", "inquire")
call nc_check(nf90_redef(ncFileID), "nc_write_model_atts", "redef")

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension. 
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name="copy",       dimid=MemberDimID), &
                            "nc_write_model_atts", "inq_dimid copy")
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="time",       dimid=TimeDimID), &
                            "nc_write_model_atts", "inq_dimid time")

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)"Time Dimension ID ",TimeDimID, &
                     " should equal Unlimited Dimension ID",unlimitedDimID
   call error_handler(E_ERR,"nc_write_model_atts", string1, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------

call nc_check(nf90_def_dim(ncid=ncFileID, name="StateVariable",  &
                           len=model_size, dimid=StateVarDimID), &
                           "nc_write_model_atts", "def_dim state")

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date" ,str1), &
                          "nc_write_model_atts", "put_att creation_date")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source"  ,source), &
                          "nc_write_model_atts", "put_att model_source")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision), &
                          "nc_write_model_atts", "put_att model_revision")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate" ,revdate), &
                          "nc_write_model_atts", "put_att model_revdate")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","template"), &
                          "nc_write_model_atts", "put_att model")

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

   ! Define the actual (3D) state vector, which gets filled as time goes on ... 
   call nc_check(nf90_def_var(ncid=ncFileID, name="state", xtype=NF90_REAL, &
                 dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), &
                 varid=VarID), "nc_write_model_atts", "def_var state")
   call nc_check(nf90_put_att(ncFileID, VarID, "long_name", "model state or fcopy"), &
                             "nc_write_model_atts", "put_att state long_name")

   ! Define the state vector coordinate variable and some attributes.
   call nc_check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=NF90_INT, &
                              dimids=StateVarDimID, varid=VarID), &
                             "nc_write_model_atts", "def_var StateVariable")
   call nc_check(nf90_put_att(ncFileID, VarID,"long_name","State Variable ID"), &
                             "nc_write_model_atts", "put_att StateVariable long_name")
   call nc_check(nf90_put_att(ncFileID, VarID, "units",     "indexical"), &
                             "nc_write_model_atts", "put_att StateVariable units")
   call nc_check(nf90_put_att(ncFileID, VarID, "valid_range", (/ 1, model_size /)), &
                             "nc_write_model_atts", "put_att StateVariable valid_range")

   ! Leave define mode so we can fill the coordinate variable.
   call nc_check(nf90_enddef(ncfileID),"nc_write_model_atts", "state_vector enddef")

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, VarID, (/ (i,i=1,model_size) /)), &
                                    "nc_write_model_atts", "put_var state")

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_dim(ncid=ncFileID, name="lon", len=nx, dimid = lonDimID), &
                     "nc_write_model_atts", "lon def_dim "//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name="lat", len=ny, dimid = latDimID), &
                     "nc_write_model_atts", "lat def_dim "//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name="lev", len=nz, dimid = levDimID), &
                     "nc_write_model_atts", "lev def_dim "//trim(filename))

   ! Standard Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name="LON", xtype=nf90_real, &
                    dimids=(/ lonDimID /), varid=VarID), &
                    "nc_write_model_atts", "LON def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, "long_name", "longitudes of grid"), &
                    "nc_write_model_atts", "LON long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, "cartesian_axis", "X"),  &
                    "nc_write_model_atts", "LON cartesian_axis "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, "units", "degrees_east"), &
                    "nc_write_model_atts", "LON units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, "valid_range", (/ -180.0_r8, 360.0_r8 /)), &
                    "nc_write_model_atts", "LON valid_range "//trim(filename))

   ! Standard Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name="LAT", xtype=nf90_real, &
                    dimids=(/ latDimID /), varid=VarID), &
                    "nc_write_model_atts", "LAT def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, "long_name", "latitudes of grid"), &
                    "nc_write_model_atts", "LAT long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, "cartesian_axis", "Y"),  &
                    "nc_write_model_atts", "LAT cartesian_axis "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, "units", "degrees_north"), &
                    "nc_write_model_atts", "LAT units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, "valid_range", (/ -90.0_r8, 90.0_r8 /)), &
                    "nc_write_model_atts", "LAT valid_range "//trim(filename))

   ! Standard Z Levels
   call nc_check(nf90_def_var(ncFileID,name="LEV", xtype=nf90_real, &
                    dimids=(/ levDimID /), varid=VarID), &
                    "nc_write_model_atts", "LEV def_var "//trim(filename))
!  call nc_check(nf90_put_att(ncFileID,  VarID, "long_name", "standard hybrid model levels"), &
!                   "nc_write_model_atts", "LEV long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, "cartesian_axis", "Z"),  &
                    "nc_write_model_atts", "LEV cartesian_axis "//trim(filename))
!  call nc_check(nf90_put_att(ncFileID,  VarID, "units", "model level"), &
!                   "nc_write_model_atts", "LEV units "//trim(filename))
!  call nc_check(nf90_put_att(ncFileID,  VarID, "valid_range", (/ 1._r8,float(nz)+1._r8 /)), &
!                   "nc_write_model_atts", "LEV valid_range "//trim(filename))

   ! DEBUG block to check shape of netCDF variables
   if ( 1 == 1 ) then
         write(*,*)"lon   dimid is ",lonDimID
         write(*,*)"lat   dimid is ",latDimID
         write(*,*)"lev   dimid is ",levDimID
         write(*,*)"unlim dimid is ",unlimitedDimID
         write(*,*)"copy  dimid is ",MemberDimID
   endif

   ! Define variables that are (lon,lat[,lev],copy,time) for each 'prognostic' variable.

   do ivar = 1, nvars

      string1 = trim(filename)//" "//trim(progvar(ivar)%varname)

      dims(1)=lonDimID
      dims(2)=latDimID

      idim=3
      if (progvar(ivar)%numvertical > 1) then
         dims(idim) = levDimID
         idim       = idim+1
      end if

      ! Create a dimension for the ensemble member/copy
      dims(idim) = memberDimID
      idim       = idim+1

      ! And now for the TIME dimension
      dims(idim) = unlimitedDimID
      ndims      = idim

      ! check shape of netCDF variables

      if ( 1 == 1 ) write(*,*)trim(progvar(ivar)%varname)," has netCDF dimIDs ",dims(1:ndims)

      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(progvar(ivar)%varname), &
                       xtype=nf90_real, dimids = dims(1:ndims), varid=VarID), &
                       "nc_write_model_atts", trim(string1)//" def_var" )

      call nc_check(nf90_put_att(ncFileID, VarID, "DART_kind", progvar(ivar)%dart_kind), &
                       "nc_write_model_atts", trim(string1)//" put_att dart_kind" )

      call nc_check(nf90_put_att(ncFileID, VarID, "DART_kind_str", trim(progvar(ivar)%kind_string)), &
                       "nc_write_model_atts", trim(string1)//" put_att dart_kind_str" )
   enddo

   call nc_check(nf90_enddef(ncfileID), "nc_write_model_atts", "prognostic enddef")

   ! Fill the coordinate variables

   call nc_check(nf90_inq_varid(ncFileID, "LON", VarID), &
                             "nc_write_model_atts", "inq_varid lon" )
   call nc_check(nf90_put_var(ncFileID, VarID, xlon),    &
                             "nc_write_model_atts", "put_var lon")

   call nc_check(nf90_inq_varid(ncFileID, "LAT", VarID), &
                             "nc_write_model_atts", "inq_varid lat" )
   call nc_check(nf90_put_var(ncFileID, VarID, xlat),    &
                             "nc_write_model_atts", "put_var lat")

!  call nc_check(nf90_inq_varid(ncFileID, "LEV", VarID), &
!                            "nc_write_model_atts", "inq_varid lev" )
!  call nc_check(nf90_put_var(ncFileID, VarID, dunnowhatgoeshere),    &
!                            "nc_write_model_atts", "put_var lev")

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID),"nc_write_model_atts", "sync")

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts



function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
! TJH 24 Oct 2011 -- Writes the model variables to a netCDF file.

use typeSizes
use netcdf

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: VarID

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME)          :: varname
integer :: ivar, i, ndims, dimlen, ncNdims, vardims(3)
integer :: TimeDimID, CopyDimID

real(r8), allocatable, dimension(:,:)   :: data_2d_array
real(r8), allocatable, dimension(:,:,:) :: data_3d_array

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

call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                          "nc_write_model_vars", "inquire")

if ( output_state_vector ) then

   call nc_check(nf90_inq_varid(ncFileID, "state", VarID), &
                             "nc_write_model_vars", "inq_varid state" )
   call nc_check(nf90_put_var(ncFileID, VarID, statevec,  &
                              start=(/ 1, copyindex, timeindex /)), &
                             "nc_write_model_vars", "put_var state")                   

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   ! The desire is to create a mystart,mycount array that can be used to
   ! simply blast the entire variable into the right slot in the netCDF file.
   !----------------------------------------------------------------------------

   do ivar=1,nvars

      varname = trim(progvar(ivar)%varname)
      string1 = trim(filename)//" "//trim(varname)

      mystart(:) = 1
      mycount(:) = 1
      dimIDs     = 0

      ! Ensure netCDF variable is conformable with progvar quantity.
      ! The TIME and Copy dimensions are intentionally not queried
      ! by looping over the dimensions stored in the progvar type.

      call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
                    "nc_write_model_vars", "inq_varid "//trim(string1))

      call nc_check(nf90_inquire_variable(ncFileID,VarId,dimids=dimIDs,ndims=ncNdims), &
                    "nc_write_model_vars", "inquire "//trim(string1))

      vardims(1) = progvar(ivar)%dimlens(1) ! num lons
      vardims(2) = progvar(ivar)%dimlens(2) ! num lats
      ndims      = 2                        ! only lon/lat

      if (progvar(ivar)%numvertical > 1) then
         vardims(3) = progvar(ivar)%dimlens(3) ! num vert
         ndims      = 3                        ! lon/lat/vert
      endif

      DimCheck : do i = 1,ndims

         write(string1,'(a,i2,a)') "inquire dimension ",i,trim(varname)
         call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
                       "nc_write_model_vars", trim(string1))

         if ( dimlen /= vardims(i) ) then
            write(string1,*) trim(varname)," dim/dimlen ",i,dimlen," not ",vardims(i)
            write(string2,*)" but it should be."
            call error_handler(E_ERR, "nc_write_model_vars", trim(string1), &
                            source, revision, revdate, text2=trim(string2))
         endif

         mycount(i) = dimlen

      enddo DimCheck

      where(dimIDs == CopyDimID) mystart = copyindex
      where(dimIDs == CopyDimID) mycount = 1
      where(dimIDs == TimeDimID) mystart = timeindex
      where(dimIDs == TimeDimID) mycount = 1

      if (ndims==2) then
         allocate(data_2d_array(vardims(1),vardims(2)))
         call vector_to_prog_var(statevec, ivar, data_2d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
                         start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                         "nc_write_model_vars", "put_var "//trim(string2))
         deallocate(data_2d_array)

      elseif (ndims==3) then
         allocate(data_3d_array(vardims(1),vardims(2),vardims(3)))
         call vector_to_prog_var(statevec, ivar, data_3d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
                       start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                       "nc_write_model_vars", "put_var "//trim(string2))
         deallocate(data_3d_array)

      else
         write(string1, *) "no support for data array of dimension ", ncNdims
         call error_handler(E_ERR,"nc_write_model_vars", string1, &
                       source,revision,revdate)
      endif
   enddo

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), "nc_write_model_vars", "sync")

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
! perturbing of states.  The returned pert_state should in any
! case be valid, since it will be read by filter even if 
! interf_provided is .false.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

if ( .not. module_initialized ) call static_init_model

pert_state      = state
interf_provided = .false.

end subroutine pert_model_state




subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! Store a copy of the ensemble mean in the model_mod module for use
! by routines that might need it. Most commonly, using a single
! definition to find the 'height' of an observation when the model
! uses a hybrid vertical coordinate system.

real(r8), intent(in) :: ens_mean(:)

if ( .not. module_initialized ) call static_init_model

module_ensemble_mean = ens_mean

end subroutine ens_mean_for_model




subroutine analysis_file_to_statevector(path, state_vector, ens_num, model_time)
!------------------------------------------------------------------
! Smooshes ensemble files into a single state vector

character(len=*), intent(in)    :: path 
real(r8),         intent(inout) :: state_vector(:)
integer,          intent(in)    :: ens_num 
type(time_type),  intent(out)   :: model_time

       integer                         :: x, y, s, z, lun, rel_offset, base_offset
       INTEGER                         :: i !_state vector index
       CHARACTER(len=5)                :: member_dir
       CHARACTER(len=256)              :: file_aod, file_conc
       REAL(r8)                        :: f_aod(nspecies+2)
       REAL(r4)                        :: f_conc(nx,ny,nz,nspecies)
       integer             :: icdtg,fhr
       type(time_type)                 :: file_time
       if ( .not. module_initialized ) call static_init_model

       !_Columns of AOD files 
       !_statevector will be in the order of (nens,nx,ny,nz,nspecies)+(nens,nx,ny,nspecies)
       ! Concentrations first, aod second?  Or just work with aod?
       !_CURRENTLY AOD ONLY
       
       !_Is NAVDAS our foward operator?  (and backward operator?)
       ! Take out the 2d<->3d portion of NAVDAS, use that to operate
       !_Currently just going to convert AOD fields.
       ! Loop over ensemble members, check file existence
        
       state_vector = MISSING_R8
       WRITE(member_dir,'(A1,I0.2,A2)') 'E', ens_num, '00'
       file_aod = trim(path) // '/NAAPSAOD/' // trim(member_dir) &
                // '/' // trim(dtg) // '_aod'
       file_conc = trim(path) // '/NAAPS/' // trim(member_dir) &
                // '/' // trim(dtg) // '_conc'
      
       if (debug) write(*,*)'analysis_file_to_statevector:Reading AOD  file',trim(file_aod)
       if (debug) write(*,*)'analysis_file_to_statevector:Reading CONC file',trim(file_conc)
 
       !_Check for file existence
       lun = OPEN_FILE(file_aod, form='formatted', action='read')
       !write(*,*) 'COORDS', nx, ny, nspecies
       !_Loop over x, y, s
       DO x = 1, nx
       DO y = 1, ny          
           READ(lun,*) f_aod 
           DO s = 1, nspecies
               i = (nx*ny*(s-1)) + (nx*(y-1)) + x !_(s,lat,lon) 
               state_vector(i) = f_aod(2 + s)
               !write(*,*) x, y, s, i 
               !write(*,*) state_vector(i)
           ENDDO
       ENDDO
       ENDDO
       call CLOSE_FILE(lun)

       !_Read concentration data in (look to readn.f in NAVDAS)
       !_CURRENTLY AOD ONLY

       !_Read in CONC data
       lun = OPEN_FILE(file_conc, form='unformatted', action='read')
       read(lun) icdtg, fhr
       read(lun)            ! nx, ny, nz, ns 
       read(lun)            ! lats 
       read(lun)            ! lons 
       read(lun)            ! siga, sigb 
       read(lun)            ! height
       read(lun)            ! binrad 
       read(lun)            ! binradw 
       read(lun)            ! rho
       read(lun)            ! refract_r
       read(lun)            ! refract_i
       read(lun)            ! budget_vars
       read(lun) f_conc     ! conc
       call CLOSE_FILE(lun)

       model_time = dtg_to_time(dtg)
       file_time = dtg_to_time(icdtg)
       if ( file_time /= model_time) then
           write(string1, *)'DART valid time does not equal analysis file time '  
           write(string2, *)'DART: ', dtg, ' FILE: ', icdtg
           call error_handler(E_ERR,'analysis_file_to_statevector', string1, &
                    source, revision, revdate, text2=string2)
       endif
      
       DO s = nspecies+1, nvars
          base_offset = progvar(s)%index1
          rel_offset = 0 
          DO z = 1, nz
          DO y = 1, ny
          DO x = 1, nx
               state_vector(base_offset + rel_offset) = f_conc(x,y,z,s-nspecies)
               rel_offset = rel_offset+1
!               write(*,*) 'CONC:', x, y, z, s, base_offset + rel_offset 
               !write(*,*) state_vector(i)
          ENDDO
          ENDDO
          ENDDO
       ENDDO

end subroutine analysis_file_to_statevector



subroutine statevector_to_analysis_file( statevector, naaps_restart_path, ens_num )
       CHARACTER(len=*), INTENT(in)    :: naaps_restart_path 
       INTEGER         , INTENT(in)    :: ens_num 
       REAL(r8),         INTENT(inout) :: statevector(:)
       CHARACTER(len=256)              :: file_concda, file_conc
       CHARACTER(len=5)                :: member_dir
       REAL(r4)                        :: height(nx,ny), binrad(nspecies), binradw(nspecies), &
                                          rho(nspecies), refract_i(nspecies), refract_r(nspecies), &
                                          budget_vars(15), mix(nx,ny), lift(nx,ny,nspecies),       &
                                          sinkd(nx,ny,nspecies), sinkw(nx,ny,nspecies),            &
                                          temperature(nx,ny,nz), sfc_pressure(nx,ny), siga(nz+1),  & 
                                          sigb(nz+1), conc(nx,ny,nz,nspecies)
       INTEGER                         :: n, i, j, k, l, s, lun, icdtg, fhr
       real(r8)                        :: tmp(nx,ny,nz) 
       REAL(r4)                        :: lats(ny), lons(nx)
       type(time_type)                 :: file_time, model_time

       if ( .not. module_initialized ) call static_init_model

       !_read in bootstrap _conc

       WRITE(member_dir,'(A1,I0.2,A2)') 'E', ens_num, '00'
       file_conc = trim(naaps_restart_path) // '/NAAPS/' // trim(member_dir) &
                // '/' // trim(dtg) // '_conc'
       file_concda = trim(naaps_restart_path) // '/NAAPS/' // trim(member_dir) &
                // '/' // trim(dtg) // '_dart'!_temporarily named for debugging
       print *, 'OUTPUT: ', trim(file_concda)

       if (debug) write(*,*)'statevector_to_analysis_file:Reading CONC file',trim(file_conc)
       if (debug) write(*,*)'statevector_to_analysis_file:Writing CONC file',trim(file_concda)

       !_Read in CONC data
       lun = OPEN_FILE(file_conc, form='unformatted', action='read')
       read(lun) icdtg, fhr   ! icdtg, fhr
       read(lun)              ! nx, ny, nz, ns 
       read(lun) lats         ! lats 
       read(lun) lons         ! lons 
       read(lun) siga, sigb   ! siga, sigb 
       read(lun) height       ! height
       read(lun) binrad       ! binrad 
       read(lun) binradw      ! binradw 
       read(lun) rho          ! rho
       read(lun) refract_r    ! refract_r
       read(lun) refract_i    ! refract_i
       read(lun) budget_vars  ! budget_vars
       read(lun)              ! conc
       read(lun) mix          ! mix
       read(lun) lift         ! lift
       read(lun) sinkd        ! sinkd
       read(lun) sinkw        ! sinkw
       read(lun) temperature  ! temperature
       read(lun) sfc_pressure ! sfc_pressure

       !_either rewind or create separate output path - write out
       file_time = dtg_to_time(icdtg)
       model_time = dtg_to_time(dtg)
       if ( file_time /= model_time) then
           write(string1, *)'DART valid time does not equal analysis file time '  
           write(string2, *)'DART: ', dtg, ' FILE: ', icdtg
           call error_handler(E_ERR,'statevector_to_analysis_file', string1, &
                    source, revision, revdate, text2=string2)
       endif
       call CLOSE_FILE(lun)

       !_rehape concentration portion of sv (i,j,k,s)
       DO n = nspecies + 1, nvars
           CALL vector_to_3d_prog_var( statevector, n, tmp )
           conc(:,:,:,n-nspecies) = REAL(tmp,r4)
       ENDDO

       !_shove into conc, write file
       lun = OPEN_FILE(file_concda, form='unformatted', action='write')
       WRITE(lun) icdtg, fhr 
       WRITE(lun) nx, ny, nz, nspecies 
       WRITE(lun) (lats(j),j=1,ny)
       WRITE(lun) (lons(i),i=1,nx)
       WRITE(lun) (siga(k),k=1,nz+1), (sigb(k),k=1,nz+1)
       WRITE(lun) ((height(i,j),i=1,nx),j=1,ny)
       WRITE(lun) binrad
       WRITE(lun) binradw
       WRITE(lun) rho     
       WRITE(lun) refract_r
       WRITE(lun) refract_i 
       WRITE(lun) budget_vars
       WRITE(lun) ((((conc(i,j,k,l),i=1,nx),j=1,ny),k=1,nz),l=1,nspecies)      
       WRITE(lun) ((mix(i,j),i=1,nx),j=1,ny)
       WRITE(lun) (((lift(i,j,s), i=1, nx), j=1, ny), s=1, nspecies)
       WRITE(lun) (((sinkd(i,j,s), i=1, nx), j=1, ny), s=1, nspecies)
       WRITE(lun) (((sinkw(i,j,s), i=1, nx), j=1, ny), s=1, nspecies)
       WRITE(lun) (((temperature(i,j,k),i=1,nx),j=1,ny),k=1,nz)        
       WRITE(lun) ((sfc_pressure(i,j), i=1, nx), j=1, ny)
       CALL CLOSE_FILE(lun)
 
END SUBROUTINE statevector_to_analysis_file



subroutine get_naaps_restart_path(path)
   character(len=*), intent(out) :: path
   if ( .not. module_initialized ) call static_init_model
   path = trim(naaps_restart_path)
end subroutine get_naaps_restart_path



subroutine get_naaps_dtg( mydtg )
   character(len=*), intent(out) :: mydtg 
   if ( .not. module_initialized ) call static_init_model
   mydtg = dtg
end subroutine get_naaps_dtg



subroutine get_naaps_ensemble_member( myne )
   integer, intent(out) :: myne
   if ( .not. module_initialized ) call static_init_model
   myne = member
end subroutine get_naaps_ensemble_member



subroutine get_naaps_metadata(path, dtg, model_time )

       CHARACTER(len=*),  INTENT(in)  :: path
       CHARACTER(len=*),  INTENT(in)  :: dtg
       TYPE(time_type),   INTENT(out) :: model_time
       CHARACTER(len=256)             :: filename 
       INTEGER                        :: dum, lun, icdtg
       REAL(r4), ALLOCATABLE          :: asig(:), bsig(:)
       
       if ( .not. module_initialized ) call static_init_model

       write(filename,'(a,''/NAAPS/E'',i2.2,''00/'',a,''_conc'')') &
                    trim(path),       member,     trim(dtg)
       !filename = trim(path) // '/NAAPS/E2000/' // trim(dtg) // '_conc'
       !filename = '/aerosol_opstmp/native/naaps/conc/201110/2011101000_conc'

       if (debug) write(*,*) 'get_naaps_metadata: opening ',trim(filename)

       lun = get_unit() 
       open(lun,file=filename,status='old',form='unformatted')!,iostat=istat)
       read(lun) icdtg, dum
       read(lun) nx, ny, nz, nspecies
       allocate(xlat(ny), xlon(nx), asig(nz+1), bsig(nz+1))
       !read(lun) !itype, ipack
       read(lun) xlat 
       read(lun) xlon
       where ( xlon < 0.0_r4 ) xlon = xlon + 360.0_r4
       
       !read(lun) asig, bsig !(nz+1, nz+1)
       close(lun)

       !_Assuming evenly spaced lat/lon grid. 
       dlat = ABS(xlat(2) - xlat(1)) 
       dlon = ABS(xlon(2) - xlon(1))
       if (debug) write(*,*) dlat, dlon, 'DLATLON'

       model_time = dtg_to_time(dtg)
       call print_date(model_time,str='get_naaps_metadata: model time')

       deallocate(asig,bsig)
       ! DO NOT deallocate xlat, xlon ... they are to remain 

end subroutine get_naaps_metadata



subroutine vector_to_1d_prog_var(x, ivar, data_1d_array)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset,
! into a 1d array.
!
real(r8), dimension(:), intent(in)  :: x
integer,                intent(in)  :: ivar
real(r8), dimension(:), intent(out) :: data_1d_array

integer :: idim1,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar(ivar)%index1

do idim1 = 1,progvar(ivar)%dimlens(1)
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



subroutine vector_to_2d_prog_var(x, ivar, data_2d_array)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset,
! into a 2d array.
!
real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:,:), intent(out) :: data_2d_array

integer :: idim1,idim2,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar(ivar)%index1

do idim2 = 1,progvar(ivar)%dimlens(2)
do idim1 = 1,progvar(ivar)%dimlens(1)
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



subroutine vector_to_3d_prog_var(x, ivar, data_3d_array)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset,
! into a 3d array.
!
real(r8), dimension(:),     intent(in)  :: x
integer,                    intent(in)  :: ivar
real(r8), dimension(:,:,:), intent(out) :: data_3d_array

integer :: idim1,idim2,idim3,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar(ivar)%index1

do idim3 = 1,progvar(ivar)%dimlens(3)
do idim2 = 1,progvar(ivar)%dimlens(2)
do idim1 = 1,progvar(ivar)%dimlens(1)
   !print *, idim3, idim2, idim1, ii, ivar, 'walter'
   !print *, progvar(ivar)%dimlens(3), 
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


function dtg_integer_to_time(mydtg)
integer, intent(in) :: mydtg
type(time_type)     :: dtg_integer_to_time

integer :: dummy,yyyy,mm,dd,hh,mn,ss

yyyy  =        mydtg/1000000  ! integer arithmetic truncates
dummy = mydtg - yyyy*1000000
mm    =        dummy/10000
dummy = dummy -   mm*10000
dd    =        dummy/100
hh    = dummy -   dd*100
mn    = 0
ss    = 0

dtg_integer_to_time = set_date(yyyy,mm,dd,hh,mn,ss)

if (debug) then
   call print_time(dtg_integer_to_time, str='dtg_integer_to_time')
   call print_date(dtg_integer_to_time, str='dtg_integer_to_time')
endif

end function dtg_integer_to_time



function dtg_string_to_time(mydtg)
character(len=*), intent(in) :: mydtg
type(time_type)              :: dtg_string_to_time

character(len=20) :: lj_dtg
integer           :: yyyy,mm,dd,hh,mn,ss

lj_dtg = adjustl(mydtg)
read(lj_dtg,'(i4,i2,i2,i2)') yyyy, mm, dd, hh
mn = 0
ss = 0
dtg_string_to_time = set_date(yyyy,mm,dd,hh,mn,ss)

if (debug) then
   call print_time(dtg_string_to_time, str='dtg_string_to_time')
   call print_date(dtg_string_to_time, str='dtg_string_to_time')
endif

end function dtg_string_to_time

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
