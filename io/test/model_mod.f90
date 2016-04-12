! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: model_mod.f90 9198 2015-12-09 16:23:52Z nancy $

module model_mod

! This is a template showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. The public interfaces listed
! must all be supported with the argument lists as indicated. Many of the interfaces
! are not required for minimal implementation (see the discussion of each
! interface and look for NULL INTERFACE). 

! Modules that are absolutely required for use are listed
use        types_mod, only : r8, i8, MISSING_R8
use time_manager_mod, only : time_type, set_time, set_time_missing, set_calendar_type, &
                             GREGORIAN
use     location_mod, only : location_type, get_close_type, get_close_maxdist_init, &
                             get_close_obs_init, loc_get_close_obs => get_close_obs, &
                             set_location, set_location_missing
use    utilities_mod, only : register_module, error_handler, nc_check, &
                             E_ERR, E_MSG, to_upper, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, check_namelist_read

use ensemble_manager_mod,  only : ensemble_type

use state_structure_mod,  only : add_domain, get_domain_size, state_structure_info

use obs_kind_mod,  only : get_raw_obs_kind_index

use dart_time_io_mod, only : write_model_time

use netcdf

implicit none
private

! required by DART code - will be called from filter and other
! DART executables.  interfaces to these routines are fixed and
! cannot be changed in any way.
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
          pert_model_copies,      &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          query_vert_localization_coord, &
          vert_convert, &
          construct_file_name_in, &
          read_model_time, &
          write_model_time

! not required by DART but for larger models can be useful for
! utility programs that are tightly tied to the other parts of
! the model_mod code.
public :: model_file_to_dart_vector, &
          dart_vector_to_model_file


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_cf_conventions/models/template/model_mod.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 9198 $"
character(len=128), parameter :: revdate  = "$Date: 2015-12-09 09:23:52 -0700 (Wed, 09 Dec 2015) $"

! global variables
character(len=512) :: string1, string2, string3

integer(i8) :: model_size

! needed for get_model_time_step
type(time_type) :: time_step

! Codes for interpreting the columns of the variable_table
integer, parameter :: VT_VARNAMEINDX  = 1 ! ... variable name
integer, parameter :: VT_KINDINDX     = 2 ! ... DART kind
integer, parameter :: VT_MINVALINDX   = 3 ! ... minimum value if any
integer, parameter :: VT_MAXVALINDX   = 4 ! ... maximum value if any
integer, parameter :: VT_STATEINDX    = 5 ! ... update (state) or not

! dimensions for the namelist state variable table 

integer, parameter :: max_state_variables = 10
integer, parameter :: num_state_table_columns = 5

integer, allocatable :: state_kinds_list(:)

! this defines what variables are in the state, and in what order.                  !------------- 
character(len=13) :: state_variables(num_state_table_columns*max_state_variables) = '             '
logical :: output_state_vector = .true.


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

integer  :: nfields, domid
integer  :: statekindlist(max_state_variables)
real(r8) :: stateclampvals(max_state_variables,2)
logical  :: stateupdatelist(max_state_variables)                                         
character(len=129), dimension(max_state_variables, num_state_table_columns):: var_table = ' '
integer  :: numvar = 4

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! set calendar type
call set_calendar_type(GREGORIAN)

                      !----5----10--     -------------    -------------    -------------    -------------
state_variables(1:numvar*num_state_table_columns)  = [character(len=13) ::  &
                   (/ 'A            ',  'KIND_PRESSURE', 'NA           ', 'NA           ', 'UPDATE       ', &
                      'B            ',  'KIND_PRESSURE', '0.0          ', 'NA           ', 'UPDATE       ', &
                      'C            ',  'KIND_PRESSURE', 'NA           ', '3.0          ', 'UPDATE       ', &
                      'D            ',  'KIND_PRESSURE', '0.0          ', '200.0        ', 'NO_COPY_BACK ' /)]

call parse_variable_table( state_variables, nfields, var_table, statekindlist, stateclampvals, stateupdatelist )

domid = add_domain('cf_test.nc', nfields, var_names = var_table(1:nfields, VT_VARNAMEINDX), &
                   kind_list   = statekindlist(1:nfields), &
                   clamp_vals  = stateclampvals(1:nfields,:), &
                   update_list = stateupdatelist(1:nfields))

call state_structure_info(domid)

model_size = get_domain_size(domid)
! ! The time_step in terms of a time type must also be initialized.
time_step = set_time(0, 0)

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

x = MISSING_R8

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

integer(i8) :: get_model_size

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




subroutine model_interpolate(state_handle, ens_size, location, obs_type, expected_obs, istatus)
!------------------------------------------------------------------
!
! Given a state handle, a location, and a model state variable type,
! interpolates the state variable fields to that location and returns
! the values in expected_obs. The istatus variables should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is a model specific integer that specifies the kind of field (for
! instance temperature, zonal wind component, etc.). In low order
! models that have no notion of types of variables this argument can
! be ignored. For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.



type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_type
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

! This should be the result of the interpolation of a
! given kind (itype) of variable at the given location.
expected_obs(:) = MISSING_R8

! The return code for successful return should be 0. 
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.
istatus(:) = 1

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



subroutine get_state_meta_data(state_handle, index_in, location, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument kind
! can be returned if the model has more than one type of field (for
! instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

type(ensemble_type), intent(in)  :: state_handle
integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

! these should be set to the actual location and obs kind
location = set_location_missing()
if (present(var_type)) var_type = 0  

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

! good style ... perhaps you could deallocate stuff (from static_init_model?).
! deallocate(state_loc)

end subroutine end_model



function nc_write_model_atts( ncFileID, model_mod_writes_state_variables ) result (ierr)
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
logical, intent(out) :: model_mod_writes_state_variables
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)

character(len=129)    :: errstring

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! and then put into define mode.
!-------------------------------------------------------------------------------

ierr = -1 ! assume things go poorly

!-------------------------------------------------------------------------------
! Have dart write out the prognostic variables
!-------------------------------------------------------------------------------
model_mod_writes_state_variables = .false.

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! and then put into define mode.
!-------------------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
                                   'nc_write_model_atts', 'inquire ')
call nc_check(nf90_Redef(ncFileID),'nc_write_model_atts',   'redef ')

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension. 
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID), &
                            "nc_write_model_atts", "inq_dimid copy")
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid= TimeDimID), &
                            "nc_write_model_atts", "inq_dimid time")

if ( TimeDimID /= unlimitedDimId ) then
   write(errstring,*)"Time Dimension ID ",TimeDimID, &
                     " should equal Unlimited Dimension ID",unlimitedDimID
   call error_handler(E_ERR,"nc_write_model_atts", errstring, source, revision, revdate)
endif

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

! Finished with dimension/variable definitions, must end 'define' mode to fill.

call nc_check(NF90_enddef(ncfileID), 'prognostic enddef ')


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

integer :: StateVarID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
!-------------------------------------------------------------------------------

ierr = -1 ! assume things go poorly

call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                          "nc_write_model_vars", "inquire")

if ( output_state_vector ) then

   call nc_check(nf90_inq_varid(ncFileID, "state", StateVarID), &
                               "nc_write_model_vars", "inq_varid state" )
   call nc_check(nf90_put_var(ncFileID, StateVarID, statevec,  &
                              start=(/ 1, copyindex, timeindex /)), &
                             "nc_write_model_vars", "put_var state")                   

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

   ! call nc_check(nf90_inq_varid(ncFileID, "ps", psVarID), &
   !                             "nc_write_model_vars",  "inq_varid ps")
   ! call nc_check(nf90_put_var( ncFileID, psVarID, global_Var%ps, &
   !                             start=(/ 1, 1, copyindex, timeindex /)), &
   !                            "nc_write_model_vars", "put_var ps")

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

pert_state      = state
interf_provided = .false.

end subroutine pert_model_state


!------------------------------------------------------------------

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

 type(ensemble_type), intent(inout) :: state_ens_handle
 integer,             intent(in)    :: ens_size
 real(r8),            intent(in)    :: pert_amp
 logical,             intent(out)   :: interf_provided

! Perturbs a model state copies for generating initial ensembles.
! The perturbed state is returned in pert_state.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding a perturbation to each model state
! variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.

interf_provided = .false.
return

end subroutine pert_model_copies

!--------------------------------------------------------------------

!> Pass through to the code in the locations module

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, &
                         obs_kind, num_close, close_ind, dist, state_handle)

type(ensemble_type),         intent(in)     :: state_handle
type(get_close_type),        intent(in)     :: gc
type(location_type),         intent(inout)  :: base_obs_loc, obs_loc(:)
integer,                     intent(in)     :: base_obs_kind, obs_kind(:)
integer,                     intent(out)    :: num_close, close_ind(:)
real(r8),                    intent(out)    :: dist(:)


call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                          num_close, close_ind, dist)

end subroutine get_close_obs



!--------------------------------------------------------------------
!> pass the vertical localization coordinate to assim_tools_mod
function query_vert_localization_coord()

integer :: query_vert_localization_coord

query_vert_localization_coord = 0 ! any value

end function query_vert_localization_coord

!--------------------------------------------------------------------
!> This is used in the filter_assim. The vertical conversion is done using the
!> mean state.
!> Calling this is a waste of time
subroutine vert_convert(state_handle, location, obs_kind, istatus)

type(ensemble_type), intent(in)  :: state_handle
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_kind
integer,             intent(out) :: istatus

istatus = 0

end subroutine vert_convert

!--------------------------------------------------------------------
!> construct restart file name for reading
function construct_file_name_in(stub, domain, copy)

character(len=512), intent(in) :: stub
integer,            intent(in) :: domain
integer,            intent(in) :: copy
character(len=1024)            :: construct_file_name_in

write(construct_file_name_in, '(2A, i4.4, A)') trim(stub), ".", copy, ".nc"

end function construct_file_name_in

!--------------------------------------------------------------------
!> read the time from the input file
!> Stolen from pop model_mod.f90 restart_to_sv
function read_model_time(filename)

character(len=256) :: filename
type(time_type) :: read_model_time

read_model_time = set_time_missing()


end function read_model_time



!=================================================================
! PUBLIC interfaces that aren't required by the DART code but are
! generally useful for other related utility programs.
! (less necessary for small models; generally used for larger models
! with predefined file formats and control structures.)
!==================================================================


subroutine model_file_to_dart_vector(filename, state_vector, model_time)
!------------------------------------------------------------------
! Reads the current time and state variables from a model data
! file and packs them into a dart state vector.

character(len=*), intent(in)    :: filename
real(r8),         intent(inout) :: state_vector(:)
type(time_type),  intent(out)   :: model_time

! code goes here

end subroutine model_file_to_dart_vector


subroutine dart_vector_to_model_file(state_vector, filename, statedate)
!------------------------------------------------------------------
! Writes the current time and state variables from a dart state
! vector (1d array) into a ncommas netcdf restart file.
!
real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: filename
type(time_type),  intent(in) :: statedate

! code goes here

end subroutine dart_vector_to_model_file

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

subroutine parse_variable_table( state_variables, ngood, table, kindslist, clampvals, updatelist )

character(len=*), dimension(:),   intent(in)  :: state_variables
integer,                          intent(out) :: ngood
character(len=*), dimension(:,:), intent(out) :: table
integer,                          intent(out) :: kindslist(:)
real(r8),                         intent(out) :: clampvals(:,:)
logical,                          intent(out) :: updatelist(:)

integer :: nrows, ncols, i, ios
character(len=NF90_MAX_NAME) :: varname       ! column 1
character(len=NF90_MAX_NAME) :: dartstr       ! column 2
character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
character(len=NF90_MAX_NAME) :: state_or_aux  ! column 6

real(r8) :: minvalue, maxvalue

! initalize clamp values to be missing
clampvals(:,:) = MISSING_R8
updatelist(:)  = .true.

nrows = size(table,1)
ncols = size(table,2)

! This loop just repackages the 1D array of values into a 2D array.
! We can do some miniminal checking along the way.
! Determining which file to check is going to be more complicated.

ngood = 0
MyLoop : do i = 1, nrows

   varname      = trim(state_variables(ncols*i - 4))
   dartstr      = trim(state_variables(ncols*i - 3))
   minvalstring = trim(state_variables(ncols*i - 2))
   maxvalstring = trim(state_variables(ncols*i - 1))
   state_or_aux = trim(state_variables(ncols*i    ))

   call to_upper(state_or_aux)

   table(i,VT_VARNAMEINDX) = trim(varname)
   table(i,VT_KINDINDX)    = trim(dartstr)
   table(i,VT_MINVALINDX)  = trim(minvalstring)
   table(i,VT_MAXVALINDX)  = trim(maxvalstring)
   table(i,VT_STATEINDX)   = trim(state_or_aux)

   ! set values for clamping
   read(table(i,VT_MINVALINDX),*,iostat=ios) minvalue
   if (ios == 0) clampvals(i,1) = minvalue

   read(table(i,VT_MAXVALINDX),*,iostat=ios) maxvalue
   if (ios == 0) clampvals(i,2) = maxvalue

   if ((table(i,VT_STATEINDX)  /= 'UPDATE')) updatelist(i) = .false.

   ! If the first element is empty, we have found the end of the list.
   if ( trim(table(i,1)) == ' ' ) exit MyLoop

   ! Any other condition is an error.
   if ( any(table(i,:) == ' ') ) then
      string1 = 'input.nml & state_variables not fully specified'
      string2 = 'must be 5 entries per variable. Last known variable name is'
      string3 = '['//trim(table(i,1))//'] ... (without the [], naturally)'
      call error_handler(E_ERR, 'parse_variable_table', string1, &
         source, revision, revdate, text2=string2, text3=string3)
   endif

   ! Make sure DART kind is valid
   kindslist(i) = get_raw_obs_kind_index(dartstr)
   if( kindslist(i) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'parse_variable_table',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector

   if (do_output()) then
      write(     *     ,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2)),' ', &
                                               trim(table(i,3)), ' ', trim(table(i,4)),' ', &
                                               trim(table(i,5))
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


!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_cf_conventions/models/template/model_mod.f90 $
! $Id: model_mod.f90 9198 2015-12-09 16:23:52Z nancy $
! $Revision: 9198 $
! $Date: 2015-12-09 09:23:52 -0700 (Wed, 09 Dec 2015) $
