! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
!----------------------------------------------------------------
!>
!> This is the interface between the ROMS ocean model and DART.
!> The required public interfaces arguments CANNOT be changed.
!>
!> Written in collaboration between Hernan Arango, Andy Moore, Chris Edwards
!> and the DART team. Uses/requires precomputed forward operator
!> values output directly from ROMS (i.e. the output of s4dvar.in:MODname).
!> The obs values are computed at the exact time of the obs as the model
!> is advancing (FGAT - First Guess at Appropriate Time). As a result,
!> there is currently NO model_interpolate() routine, since all
!> the observations already have expected values. What is required is
!> the ability to convert -on-demand- a DART observation sequence file
!> from the ROMS observation format. This is done with the convert_roms_obs
!> program in observations/ROMS.
!>
!> If required for other obs, the model_interpolate code will
!> need to be written and tested.
!>
!> The ROMS model uses a mask array to indicate dry land, so inside
!> DART and this model mod we can ignore dry land completely.  there
!> is no need to read the mask or account for it in get_state_meta_data
!> or get_close.
!>
!----------------------------------------------------------------

module model_mod

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, digits12, SECPERDAY, DEG2RAD, rad2deg, PI, &
                             MISSING_I, MISSING_R4, MISSING_R8, i4, i8

use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time, &
                             print_time, print_date,                            &
                             set_calendar_type, get_calendar_type,              &
                             operator(*),  operator(+), operator(-),            &
                             operator(>),  operator(<), operator(/),            &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, get_dist, get_close_maxdist_init,   &
                             get_close_obs_init, set_location,                  &
                             get_location, vert_is_height,write_location,       &
                             set_location_missing,query_location,               &
                             vert_is_level, vert_is_surface,                    &
                             loc_get_close_obs => get_close_obs, get_close_type,&
                             VERTISHEIGHT,VERTISSURFACE,VERTISUNDEF,VERTISLEVEL,&
                             VERTISPRESSURE,VERTISSCALEHEIGHT,horiz_dist_only

use    utilities_mod, only : register_module, error_handler, do_nml_term,       &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,       &
                             nc_check, do_output, to_upper, do_nml_file,        &
                             find_namelist_in_file, check_namelist_read,        &
                             open_file, file_exist, find_textfile_dims,         &
                             file_to_text, do_output, close_file,               &
                             string_to_real, string_to_logical

use     obs_kind_mod, only : QTY_TEMPERATURE,           &
                             QTY_SALINITY,              &
                             QTY_U_CURRENT_COMPONENT,   &
                             QTY_V_CURRENT_COMPONENT,   &
                             QTY_SEA_SURFACE_HEIGHT,    &
                             QTY_SEA_SURFACE_PRESSURE,  &
                             QTY_POTENTIAL_TEMPERATURE, &
                             get_index_for_quantity,     &
                             get_name_for_quantity

use     mpi_utilities_mod, only : my_task_id

use        random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

use  ensemble_manager_mod, only : ensemble_type, map_pe_to_task, get_copy_owner_index, &
                                  get_var_owner_index

use distributed_state_mod, only : get_state

use   state_structure_mod, only : add_domain, get_model_variable_indices, &
                                  get_num_variables, get_index_start, &
                                  get_num_dims, get_domain_size, get_varid_from_kind, &
                                  get_dart_vector_index, state_structure_info, &
                                  get_index_start, get_index_end, get_variable_name, &
                                  get_kind_index, get_kind_string, get_dim_length, &
                                  get_dim_name, get_missing_value, get_units, &
                                  get_long_name, get_xtype, get_has_missing_value, &
                                  get_dim_lengths

use typesizes
use netcdf

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
public :: get_model_size,                &
          adv_1step,                     &
          get_state_meta_data,           &
          model_interpolate,             &
          get_model_time_step,           &
          static_init_model,             &
          end_model,                     &
          init_time,                     &
          init_conditions,               &
          nc_write_model_atts,           &
          nc_write_model_vars,           &
          pert_model_copies,             &
          get_close_maxdist_init,        &
          get_close_obs_init,            &
          query_vert_localization_coord, &
          vert_convert,                  &
          write_model_time,              &
          read_model_time,               &
          get_close_obs

public :: get_time_information,          &
          get_location_from_ijk

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2, string3
logical, save :: module_initialized = .false.

! things which can/should be in the model_nml
logical  :: output_state_vector          = .false.
integer  :: assimilation_period_days     = 1
integer  :: assimilation_period_seconds  = 0
integer  :: vert_localization_coord      = VERTISHEIGHT
integer  :: debug = 0   ! turn up for more and more debug messages
character(len=256) :: roms_filename = 'roms_input.nc'

namelist /model_nml/  &
   output_state_vector,         &
   assimilation_period_days,    &
   assimilation_period_seconds, &
   roms_filename,               &
   vert_localization_coord,     &
   debug,                       &
   variables

! DART contents are specified in the input.nml:&model_nml namelist.
!>@todo  NF90_MAX_NAME is 256 ... this makes the namelist output unreadable
integer, parameter :: MAX_STATE_VARIABLES = 80
integer, parameter :: num_state_table_columns = 5
character(len=NF90_MAX_NAME) :: variables(MAX_STATE_VARIABLES * num_state_table_columns ) = ' '
character(len=NF90_MAX_NAME) :: var_names(MAX_STATE_VARIABLES) = ' '
logical  ::                   update_list(MAX_STATE_VARIABLES) = .FALSE.
integer  ::                     kind_list(MAX_STATE_VARIABLES) = MISSING_I
real(r8) ::                    clamp_vals(MAX_STATE_VARIABLES,2) = MISSING_R8

integer :: nfields   ! This is the number of variables in the DART state vector.

integer :: domain_id ! global variable for state_structure_mod routines

!> Everything needed to describe a variable. Basically all the metadata from
!> a netCDF file is stored here as well as all the information about where
!> the variable is stored in the DART state vector.
!>

! Grid parameters - the values will be read from a
! standard ROMS namelist and filled in here.

! nx, ny and nz are the size of the rho grids.
integer :: Nx = -1, Ny = -1, Nz = -1

integer :: Nxi_rho
integer :: Nxi_u
integer :: Nxi_v
integer :: Neta_rho
integer :: Neta_u
integer :: Neta_v
integer :: Ns_rho
integer :: Ns_w

!>@todo FIXME ... nancy suggested creating pointers for each of these so
!    we could simply use the myvarid as the index in the pointer ...

real(r8), allocatable, target :: ULAT(:,:), ULON(:,:), UDEP(:,:,:), &
                                 TLAT(:,:), TLON(:,:), TDEP(:,:,:), &
                                 VLAT(:,:), VLON(:,:), VDEP(:,:,:), &
                                 WDEP(:,:,:) !>@todo FIXME : JPH may not need this array

type(time_type) :: model_timestep

integer :: model_size    ! the state vector length

!> Reshapes a part of the DART vector back to the original variable shape.
!>@todo FIXME Replaces the DART MISSING value with the original _FillValue value.


contains


!-----------------------------------------------------------------------
! All the REQUIRED interfaces come first - by convention.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!>
!> Returns the size of the DART state vector (i.e. model) as an integer.
!> Required for all applications.
!>

function get_model_size()

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size


!-----------------------------------------------------------------------
!>
!> Does a single timestep advance of the model in a subroutine call.
!> This interface is only called if the namelist parameter
!> async is set to 0 in perfect_model_obs of filter or if the
!> program integrate_model is to be used to advance the model
!> state as a separate executable. If one of these options
!> is not going to be used (the model will only be advanced as
!> a separate model-specific executable), this can be a
!> NULL INTERFACE.
!>
!> NOTE: not supported for ROMS. Will intentionally generate a fatal error.
!>
!> @param x the model state before and after the model advance.
!> @param time the desired time at the end of the model advance.
!>

subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

if ( .not. module_initialized ) call static_init_model

write(string1,*)'Cannot advance ROMS with a subroutine call; async cannot equal 0'
write(string2,*)'Unsupported method for ROMS.'
call error_handler(E_ERR, 'adv_1step:', string1, &
                   source, revision, revdate, text2=string2)

end subroutine adv_1step


!-----------------------------------------------------------------------
!>
!> Given an integer index into the state vector structure, returns the
!> associated location. A second intent(out) optional argument kind
!> can be returned if the model has more than one type of field (for
!> instance temperature and zonal wind component). This interface is
!> required for all filter applications as it is required for computing
!> the distance between observations and state variables.
!>
!> @param state_handle DART ensemble handle
!> @param index_in the index into the DART state vector
!> @param location the location at that index
!> @param var_type the DART KIND at that index
!>

subroutine get_state_meta_data(state_handle, index_in, location, var_type)

type(ensemble_type), intent(in)  :: state_handle
integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

! Local variables

integer  :: iloc, vloc, jloc
integer  :: myindx, myvarid

if ( .not. module_initialized ) call static_init_model

myindx  = -1
myvarid = -1

call get_model_variable_indices(index_in, iloc, jloc, vloc, var_id=myvarid)

! FIXME - get this always and compare below
if (present(var_type)) then
   var_type = get_kind_index(domain_id, myvarid)
endif

if     (get_kind_index(domain_id,myvarid)==QTY_U_CURRENT_COMPONENT) then
      location = set_location(ULON(iloc,jloc),ULAT(iloc,jloc), UDEP(iloc,jloc,vloc), VERTISHEIGHT)
elseif (get_kind_index(domain_id,myvarid)==QTY_V_CURRENT_COMPONENT) then
      location = set_location(VLON(iloc,jloc),VLAT(iloc,jloc), VDEP(iloc,jloc,vloc), VERTISHEIGHT)
elseif (get_kind_index(domain_id,myvarid)==QTY_SEA_SURFACE_HEIGHT) then
      location = set_location(TLON(iloc,jloc),TLAT(iloc,jloc), 0.0_r8, VERTISSURFACE)
else  ! Everything else is assumed to be on the rho points
      location = set_location(TLON(iloc,jloc),TLAT(iloc,jloc), TDEP(iloc,jloc,vloc), VERTISHEIGHT)
endif

end subroutine get_state_meta_data


!-----------------------------------------------------------------------
!>
!> Model interpolate will interpolate any DART state variable
!> (i.e. S, T, U, V, Eta) to the given location given a state vector.
!> The type of the variable being interpolated is obs_type since
!> normally this is used to find the expected value of an observation
!> at some location. The interpolated value is returned in interp_vals
!> and istatus is 0 for success. NOTE: This is a workhorse routine and is
!> the basis for all the forward observation operator code.
!>
!> @param state_handle DART ensemble handle
!> @param ens_size DART ensemble size
!> @param location the location of interest
!> @param obs_type the DART KIND of interest
!> @param interp_val the estimated value of the DART state at the location
!>          of interest (the interpolated value).
!> @param istatus interpolation status ... 0 == success, /=0 is a failure
!>

subroutine model_interpolate(state_handle, ens_size, location, obs_type, interp_val, istatus)

 type(ensemble_type), intent(in) :: state_handle
 integer,             intent(in) :: ens_size
 type(location_type), intent(in) :: location
 integer,             intent(in) :: obs_type
 integer,            intent(out) :: istatus(ens_size)
 real(r8),           intent(out) :: interp_val(ens_size) !< array of interpolated values

if ( .not. module_initialized ) call static_init_model

! Successful istatus is 0
interp_val = MISSING_R8
istatus = 99

write(string1,*)'model_interpolate should not be called.'
write(string2,*)'we are getting forward observations directly from ROMS'
call error_handler(E_MSG,'model_interpolate:',string1,source,revision,revdate, text2=string2)

end subroutine model_interpolate


!-----------------------------------------------------------------------
!>
!> Returns the the time step of the model; the smallest increment in
!> time that the model is capable of advancing the ROMS state.
!>

function get_model_time_step()

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = model_timestep

end function get_model_time_step


!-----------------------------------------------------------------------
!>
!> Called to do one time initialization of the model.
!> In this case, it reads in the grid information, the namelist
!> containing the variables of interest, where to get them, their size,
!> their associated DART KIND, etc.
!>
!> In addition to harvesting the model metadata (grid,
!> desired model advance step, etc.), it also fills a structure
!> containing information about what variables are where in the DART
!> framework.

subroutine static_init_model()

integer :: iunit, io
integer :: ss, dd
integer :: ncid

character(len=32) :: calendar

type(time_type) :: model_time

if ( module_initialized ) return

! The Plan:
!
! * read in the grid sizes from grid file
! * allocate space, and read in actual grid values
! * figure out model timestep
! * Compute the model size.
! * set the index numbers where the field types change

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(logfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model:',string1,source,revision,revdate)

call nc_check( nf90_open(trim(roms_filename), NF90_NOWRITE, ncid), &
                  'static_init_model', 'open '//trim(roms_filename))

call get_time_information(roms_filename, ncid, 'ocean_time', 'ocean_time', &
                          calendar=calendar, last_time=model_time)

call set_calendar_type( trim(calendar) )

! Get the ROMS grid -- sizes and variables.
call get_grid_dimensions()
call get_grid()

! parse_variable_input() fills var_names, kind_list, clamp_vals, update_list
call parse_variable_input(variables, nfields)

domain_id = add_domain(roms_filename, nfields, &
                    var_names, kind_list, clamp_vals, update_list )

if (debug > 2) call state_structure_info(domain_id)

call nc_check( nf90_close(ncid), &
                  'static_init_model', 'close '//trim(roms_filename))

model_size = get_domain_size(domain_id)

call write_roms_time_information(model_time)

end subroutine static_init_model


!-----------------------------------------------------------------------
!>
!> Does any shutdown and clean-up needed for model.
!>

subroutine end_model()

! good style ... perhaps you could deallocate stuff (from static_init_model?).
! deallocate(state_loc)

if (allocated(ULAT)) deallocate(ULAT)
if (allocated(ULON)) deallocate(ULON)
if (allocated(UDEP)) deallocate(UDEP)

if (allocated(VLAT)) deallocate(VLAT)
if (allocated(VLON)) deallocate(VLON)
if (allocated(VDEP)) deallocate(VDEP)

if (allocated(TLAT)) deallocate(TLAT)
if (allocated(TLON)) deallocate(TLON)
if (allocated(TDEP)) deallocate(TDEP)

if (allocated(WDEP)) deallocate(WDEP)

end subroutine end_model


!-----------------------------------------------------------------------
!>
!> Companion interface to init_conditions. Returns a time that is somehow
!> appropriate for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter
!> start_from_restart is set to .false. in the program perfect_model_obs.
!> If this option is not to be used in perfect_model_obs, or if no
!> synthetic data experiments using perfect_model_obs are planned,
!> this can be a NULL INTERFACE.
!>
!> NOTE: Since ROMS cannot start in this manner,
!> DART will intentionally generate a fatal error.
!>
!> @param time the time to associate with the initial state
!>

subroutine init_time(time)

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

time = set_time(0,0)

write(string1,*) 'Cannot initialize ROMS time via subroutine call; start_from_restart cannot be F'
write(string2,*)'Unsupported method for ROMS.'
call error_handler(E_ERR, 'init_time:', string1, &
                   source, revision, revdate, text2=string2)

end subroutine init_time


!-----------------------------------------------------------------------
!>
!> Returns a model state vector, x, that is some sort of appropriate
!> initial condition for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter
!> start_from_restart is set to .false. in the program perfect_model_obs.
!> If this option is not to be used in perfect_model_obs, or if no
!> synthetic data experiments using perfect_model_obs are planned,
!> this can be a NULL INTERFACE.
!>
!> NOTE: This is not supported for ROMS and will generate a FATAL ERROR.
!>       However, this is a required interface - so it must be present.
!>
!> @param x the ROMS initial conditions
!>

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

x = MISSING_R8

write(string1,*)'Cannot initialize ROMS state via subroutine call.'
write(string2,*)'namelist "start_from_restart" cannot be F'
call error_handler(E_ERR, 'init_conditions:', string1, &
                   source, revision, revdate, text2=string2)

end subroutine init_conditions


!-----------------------------------------------------------------------
!>
!> Writes the model-specific attributes to a DART 'diagnostic' netCDF file.
!> This includes coordinate variables and some metadata, but NOT the
!> actual DART state.
!>
!> @param ncFileID the netCDF handle of the DART diagnostic file opened by
!>                 assim_model_mod:init_diag_output
!> @param model_writes_state have the state structure write out all of the
!>                 state variables
!> @param ierr status ... 0 == all went well, /= 0 failure

function nc_write_model_atts( ncFileID, model_writes_state ) result (ierr)

! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

integer, intent(in)  :: ncFileID           ! netCDF file identifier
logical, intent(out) :: model_writes_state ! if true, dart lib writes state info
integer              :: ierr               ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

! variables if we just blast out one long state vector

integer :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)

integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

! variables if we parse the state vector into prognostic variables.

! for the dimensions and coordinate variables
integer :: nxirhoDimID, nxiuDimID, nxivDimID
integer :: netarhoDimID, netauDimID, netavDimID
integer :: nsrhoDimID

! for the prognostic variables
integer :: VarID

! local variables

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

character(len=256) :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly
model_writes_state = .false.

! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.

write(filename,*) 'ncFileID', ncFileID

! make sure ncFileID refers to an open netCDF file,
! and then put into define mode.

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
                                   'nc_write_model_atts', 'inquire '//trim(filename))
call nc_check(nf90_Redef(ncFileID),'nc_write_model_atts',   'redef '//trim(filename))

! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension.
! Our job is create the 'model size' dimension.

call nc_check(nf90_inq_dimid(ncid=ncFileID, name='copy', dimid=MemberDimID), &
                           'nc_write_model_atts', 'copy dimid '//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid=  TimeDimID), &
                           'nc_write_model_atts', 'time dimid '//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)'Time Dimension ID ',TimeDimID, &
             ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts:', string1, source, revision, revdate)
endif

! Define the model size / state variable dimension / whatever ...

call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable', len=model_size, &
        dimid = StateVarDimID),'nc_write_model_atts', 'state def_dim '//trim(filename))

! Write Global Attributes

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
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'ROMS' ), &
           'nc_write_model_atts', 'model put '//trim(filename))

! Here is the extensible part. The simplest scenario is to output the state vector,
! parsing the state vector into model-specific parts is complicated, and you need
! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
! complicated part.

if ( output_state_vector ) then

   ! Create a variable for the state vector
   ! Define the actual (3D) state vector, which gets filled as time goes on ...

   call nc_check(nf90_def_var(ncid=ncFileID, name='state', xtype=nf90_real, &
                 dimids=(/StateVarDimID,MemberDimID,unlimitedDimID/),VarID=StateVarID),&
                 'nc_write_model_atts','state def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarID,'long_name','model state or fcopy'),&
                 'nc_write_model_atts', 'state long_name '//trim(filename))

   call nc_check(nf90_enddef(ncFileID),'nc_write_model_atts','state enddef '//trim(filename))

else

   ! We need to output the prognostic variables.
   ! Define the new dimensions IDs

   call nc_check(nf90_def_dim(ncid=ncFileID, name='xi_rho',  len = Nxi_rho, &
        dimid = nxirhoDimID),'nc_write_model_atts', 'xi_rho def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='eta_rho', len = Neta_rho,&
        dimid = netarhoDimID),'nc_write_model_atts', 'eta_rho def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='s_rho',   len = Ns_rho,&
        dimid = nsrhoDimID),'nc_write_model_atts', 's_rho def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='xi_u',    len = Nxi_u,&
        dimid = nxiuDimID),'nc_write_model_atts', 'xi_u def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='xi_v',    len = Nxi_v,&
        dimid = nxivDimID),'nc_write_model_atts', 'xi_v def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='eta_u',   len = Neta_u,&
        dimid = netauDimID),'nc_write_model_atts', 'eta_u def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='eta_v',   len = Neta_v,&
        dimid = netavDimID),'nc_write_model_atts', 'eta_v def_dim '//trim(filename))

   ! Create the (empty) Coordinate Variables and the Attributes

   call nc_check(nf90_def_var(ncFileID,name='lon_rho', xtype=nf90_double, &
                 dimids=(/ nxirhoDimID, netarhoDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lon_rho def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'rho longitudes'), &
                 'nc_write_model_atts', 'lon_rho long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'lon_rho units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='lat_rho', xtype=nf90_double, &
                 dimids=(/ nxirhoDimID, netarhoDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lat_rho def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'rho latitudes'), &
                 'nc_write_model_atts', 'lat_rho long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_north'), &
                 'nc_write_model_atts', 'lat_rho units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='lon_u', xtype=nf90_double, &
                 dimids=(/ nxiuDimID, netauDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lon_u def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'u longitudes'), &
                 'nc_write_model_atts', 'lon_u long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'lon_u units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='lat_u', xtype=nf90_double, &
                 dimids=(/ nxiuDimID, netauDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lat_u def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'u latitudes'), &
                 'nc_write_model_atts', 'lat_u long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_north'), &
                 'nc_write_model_atts', 'lat_u units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='lon_v', xtype=nf90_double, &
                 dimids=(/ nxivDimID, netavDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lon_v def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'v longitudes'), &
                 'nc_write_model_atts', 'lon_v long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'lon_v units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='lat_v', xtype=nf90_double, &
                 dimids=(/ nxivDimID, netavDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lat_v def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'v latitudes'), &
                 'nc_write_model_atts', 'lat_v long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_north'), &
                 'nc_write_model_atts', 'lat_v units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='z_rho', xtype=nf90_double, &
                 dimids=(/ nxirhoDimID, netarhoDimID /), varid=VarID),&
                 'nc_write_model_atts', 'z_rho def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'z at rho'), &
                 'nc_write_model_atts', 'z_rho long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'm'), &
                 'nc_write_model_atts', 'z_rho units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='z_u', xtype=nf90_double, &
                 dimids=(/ nxiuDimID, netauDimID /), varid=VarID),&
                 'nc_write_model_atts', 'z_u def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'z at rho'), &
                 'nc_write_model_atts', 'z_u long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'm'), &
                 'nc_write_model_atts', 'z_u units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='z_v', xtype=nf90_double, &
                 dimids=(/ nxivDimID, netavDimID /), varid=VarID),&
                 'nc_write_model_atts', 'z_v def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'z at rho'), &
                 'nc_write_model_atts', 'z_v long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'm'), &
                 'nc_write_model_atts', 'z_v units '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='z_w', xtype=nf90_double, &
                 dimids=(/ nxivDimID, netavDimID /), varid=VarID),&
                 'nc_write_model_atts', 'z_w def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'z at rho'), &
                 'nc_write_model_atts', 'z_w long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'm'), &
                 'nc_write_model_atts', 'z_w units '//trim(filename))

   ! Finished with dimension/variable definitions, must end 'define' mode to fill.

   call nc_check(nf90_enddef(ncFileID), 'prognostic enddef '//trim(filename))

   ! Fill the coordinate variables that DART needs and has locally

   ! the RHO grid

   call nc_check(NF90_inq_varid(ncFileID, 'lon_rho', VarID), &
                 'nc_write_model_atts', 'lon_rho inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, TLON ), &
                'nc_write_model_atts', 'lon_rho put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'lat_rho', VarID), &
                 'nc_write_model_atts', 'lat_rho inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, TLAT ), &
                'nc_write_model_atts', 'lat_rho put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'z_rho', VarID), &
                 'nc_write_model_atts', 'z_rho inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, TDEP ), &
                'nc_write_model_atts', 'z_rho put_var '//trim(filename))

   ! the U grid

   call nc_check(NF90_inq_varid(ncFileID, 'lon_u', VarID), &
                 'nc_write_model_atts', 'lon_u inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, ULON ), &
                'nc_write_model_atts', 'lon_u put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'lat_u', VarID), &
                 'nc_write_model_atts', 'lat_u inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, ULAT ), &
                'nc_write_model_atts', 'lat_u put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'z_u', VarID), &
                 'nc_write_model_atts', 'z_u inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, UDEP ), &
                'nc_write_model_atts', 'z_u put_var '//trim(filename))

   ! the V grid

   call nc_check(NF90_inq_varid(ncFileID, 'lon_v', VarID), &
                 'nc_write_model_atts', 'lon_v inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, VLON ), &
                'nc_write_model_atts', 'lon_v put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'lat_v', VarID), &
                 'nc_write_model_atts', 'lat_v inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, VLAT ), &
                'nc_write_model_atts', 'lat_v put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'z_v', VarID), &
                 'nc_write_model_atts', 'z_v inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, VDEP ), &
                'nc_write_model_atts', 'z_v put_var '//trim(filename))

   ! the W grid

   call nc_check(NF90_inq_varid(ncFileID, 'z_w', VarID), &
                 'nc_write_model_atts', 'z_w inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, WDEP ), &
                'nc_write_model_atts', 'z_w put_var '//trim(filename))

   endif

! Flush the buffer and leave netCDF file open
call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts


!-----------------------------------------------------------------------
!>
!> With each assimilation cycle, the DART prior and posterior files get
!> inserted into the DART diagnostic files. This routine appends the new
!> states into the unlimited dimension slot.
!>
!> @param ncFileID the netCDF file ID of the DART diagnostic file in question
!> @param state_vec the DART state to insert into the diagnostic file
!> @param copyindex the 'copy' index ... ensemble mean, member 23, etc.
!> @param timeindex the index into the unlimited (time) dimension
!> @param ierr error code. All errors are fatal. 0 == success.

function nc_write_model_vars( ncFileID, state_vec, copyindex, timeindex ) result (ierr)

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

integer,  intent(in) :: ncFileID
real(r8), intent(in) :: state_vec(:)
integer,  intent(in) :: copyindex
integer,  intent(in) :: timeindex
integer              :: ierr

integer :: VarID
integer :: TimeDimID, CopyDimID

character(len=256) :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.

write(filename,*) 'ncFileID', ncFileID

! make sure ncFileID refers to an open netCDF file,

call nc_check(nf90_inq_dimid(ncFileID, 'copy', dimid=CopyDimID), &
            'nc_write_model_vars', 'inq_dimid copy '//trim(filename))

call nc_check(nf90_inq_dimid(ncFileID, 'time', dimid=TimeDimID), &
            'nc_write_model_vars', 'inq_dimid time '//trim(filename))

if ( output_state_vector ) then

   call nc_check(NF90_inq_varid(ncFileID, 'state', VarID), &
                 'nc_write_model_vars', 'state inq_varid '//trim(filename))
   call nc_check(NF90_put_var(ncFileID,VarID,state_vec,start=(/1,copyindex,timeindex/)),&
                 'nc_write_model_vars', 'state put_var '//trim(filename))

endif

! Flush the buffer and leave netCDF file open

call nc_check(nf90_sync(ncFileID), 'nc_write_model_vars', 'sync '//trim(filename))

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars

!-----------------------------------------------------------------------
!>
!> Perturbs a model state for generating initial ensembles.
!> The perturbed state is returned in pert_state.
!> A model may choose to provide a NULL INTERFACE by returning
!> .false. for the interf_provided argument. This indicates to
!> the filter that if it needs to generate perturbed states, it
!> may do so by adding a perturbation to each model state
!> variable independently. The interf_provided argument
!> should be returned as .true. if the model wants to do its own
!> perturbing of states.
!>
!> @param state_ens_handle the DART state ensemble handle
!> @param ens_size ensemble size
!> @param pert_amp perturbation amplitude
!> @param interf_provided logical flag that indicates that this routine
!>               is unique for ROMS. TRUE means this routine will
!>               somehow create the perturbed state, FALSE means
!>               the default perturb routine will be used.
!>

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

 type(ensemble_type), intent(inout) :: state_ens_handle
 integer,             intent(in)    :: ens_size
 real(r8),            intent(in)    :: pert_amp
 logical,             intent(out)   :: interf_provided

interf_provided = .false.

end subroutine pert_model_copies


!-----------------------------------------------------------------------
!>
!> Given a DART location (referred to as "base") and a set of candidate
!> locations and kinds (obs, obs_kind); returns the subset close to the
!> "base", their indices, and their distances to the "base" ...
!>
!> @param gc precomputed 'get_close_type' to speed up candidate selection
!> @param base_obs_loc location of the observation in question
!> @param base_obs_kind DART KIND of observation in question
!> @param locs array of comparison locations
!> @param loc_kind matching array of KINDs for the comparison locations
!> @param num_close the number of locs locations that are within the prespecified distance (information contained in 'gc')
!> @param close_ind the indices of the locs locations that are 'close'
!> @param dist the distances of each of the close locations.  optional.
!> @param state_handle handle to the dart state
!>

subroutine get_close_obs(gc, base_obs_loc, base_obs_type, &
                         locs, loc_kind, num_close, close_ind, dist, state_handle)
 type(ensemble_type),               intent(in)    :: state_handle
 type(get_close_type),              intent(in)    :: gc
 type(location_type),               intent(inout) :: base_obs_loc
 integer,                           intent(in)    :: base_obs_type
 type(location_type), dimension(:), intent(in)    :: locs
 integer,             dimension(:), intent(in)    :: loc_kind
 integer,                           intent(out)   :: num_close
 integer,             dimension(:), intent(out)   :: close_ind
 real(r8), optional,  dimension(:), intent(out)   :: dist

! Note that both base_obs_loc and locs are intent(inout), meaning that these
! locations are possibly modified here and returned as such to the calling routine.
! The calling routine is always filter_assim and these arrays are local arrays
! within filter_assim. In other words, these modifications will only matter within
! filter_assim, but will not propagate backwards to filter.

!>@todo FIXME implement masking ...

! use the default system routine

call loc_get_close_obs(gc, base_obs_loc, base_obs_type, locs, loc_kind, &
                       num_close, close_ind, dist)

end subroutine get_close_obs


!--------------------------------------------------------------------
!>
!> This is used to pass the vertical localization coordinate
!> to assim_tools_mod
!>

function query_vert_localization_coord()

integer :: query_vert_localization_coord

query_vert_localization_coord = VERTISHEIGHT

end function query_vert_localization_coord

!-----------------------------------------------------------------------
!>
!> This subroutine converts a given ob/state vertical coordinate to
!> the vertical localization coordinate type requested through the
!> model_mod namelist. This is used in filter_assim(). The vertical
!> conversion is done using the mean state.
!>
!> Notes: (1) obs_kind is only necessary to check whether the ob
!>            is an identity observation.
!>
!>        (2) This subroutine can convert both obs' and state points'
!>            vertical coordinates. Remember that state points get
!>            their DART location information from get_state_meta_data
!>            which is called by filter_assim during the assimilation
!>            process.
!>
!>        (3) state_handle contains relevant DART state information for
!>            carrying out computations necessary for the vertical coordinate
!>            transformations. As the vertical coordinate is only used
!>            in distance computations, this is actually the "expected"
!>            vertical coordinate, so that computed distance is the
!>            "expected" distance. Thus, under normal circumstances,
!>            state_handle that is supplied to convert_vert should be the
!>            ensemble mean. Nevertheless, the subroutine has the
!>            functionality to operate on any DART state vector that
!>            is supplied to it.
!>
!> @param state_handle handle to the dart state
!> @param location
!> @param obs_kind
!> @param istatus
!>

subroutine vert_convert(state_handle, location, obs_kind, istatus)

type(ensemble_type), intent(in)  :: state_handle
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_kind
integer,             intent(out) :: istatus

! no need for vertical conversion since the forward observations
! are coming from ROMS and vertical coordinate is always in HEIGHT

istatus = 0

end subroutine vert_convert


!-----------------------------------------------------------------------
!>
!> writes the time of the current state and (optionally) the time
!> to be conveyed to ROMS to dictate the length of the forecast.
!> This file is then used by scripts to modify the ROMS run.
!> The format in the time information is totally at your discretion.
!>
!> @param ncfile_out name of the file
!> @param model_time the current time of the model state
!> @param adv_to_time the time in the future of the next assimilation.
!>

subroutine write_model_time(ncid, model_time, adv_to_time)
integer,         intent(in)           :: ncid
type(time_type), intent(in)           :: model_time
type(time_type), intent(in), optional :: adv_to_time

integer :: io, varid, seconds, days
type(time_type) :: origin_time, deltatime
real(digits12)  :: run_duration

if ( .not. module_initialized ) call static_init_model

if (present(adv_to_time)) then
   string3 = time_to_string(adv_to_time)
   write(string1,*)'ROMS/DART not configured to advance ROMS.'
   write(string2,*)'called with optional advance_to_time of'
   call error_handler(E_ERR, 'write_model_time', string1, &
              source, revision, revdate, text2=string2,text3=string3)
endif

! If the ocean_time variable exists, we are updating a ROMS file,
! if not ... must be updating a DART diagnostic file.

io = nf90_inq_varid(ncid,'ocean_time',varid)
if (io == NF90_NOERR) then
   call get_time_information('unknown', ncid, 'ocean_time', 'ocean_time', &
                myvarid=varid, origin_time=origin_time)
   deltatime = model_time - origin_time
   call get_time(deltatime, seconds, days)
   run_duration = real(days,digits12)*86400.0_digits12 + real(seconds,digits12)
   call nc_check(nf90_put_var(ncid, varid, run_duration), 'write_model_time', 'put_var')
   return
endif

io = nf90_inq_varid(ncid,'time',varid)
if (io == NF90_NOERR) then
   call get_time_information('unknown', ncid, 'time', 'time', &
                myvarid=varid, origin_time=origin_time)
   deltatime = model_time - origin_time
   call get_time(deltatime, seconds, days)
   run_duration = real(days,digits12)*86400.0_digits12 + real(seconds,digits12)
   call nc_check(nf90_put_var(ncid, varid, run_duration), 'write_model_time', 'put_var')
   return
endif

end subroutine write_model_time

!--------------------------------------------------------------------
!>
!> read the time from the input file
!>
!> @param filename name of file that contains the time
!>

function read_model_time(filename)

character(len=*), intent(in) :: filename
type(time_type)              :: read_model_time

integer :: ncid

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'read_model_time',string1,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'read_model_time', 'open '//trim(filename))

call get_time_information(filename, ncid, 'ocean_time', 'ocean_time', last_time=read_model_time)

call nc_check( nf90_close(ncid), 'read_model_time', 'close '//trim(filename))

end function read_model_time



!-----------------------------------------------------------------------
! The remaining (private) interfaces come last.
! None of the private interfaces need to call static_init_model()
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!>
!> Set the desired minimum model advance time. This is generally NOT the
!> dynamical timestep of the model, but rather the shortest forecast length
!> you are willing to make. This impacts how frequently the observations
!> may be assimilated.
!>

function set_model_time_step()

type(time_type) :: set_model_time_step

! assimilation_period_seconds, assimilation_period_days are from the namelist

!>@todo FIXME make sure set_model_time_step is an integer multiple of
!> the dynamical timestep or whatever strategy ROMS employs.
! TJH NHIST*DT ... from the inputfile ... can remove assim_* from DART input namelist

!>@todo FIXME : JPH we should really be getting this from the history file??
set_model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end function set_model_time_step


!-----------------------------------------------------------------------
!>
!> Read the grid dimensions from the ROMS grid netcdf file.
!> By reading the dimensions first, we can use them in variable
!> declarations later - which is faster than using allocatable arrays.
!>

subroutine get_grid_dimensions()

integer :: ncid

! Read the (static) grid dimensions from the ROMS grid file.

call nc_check(nf90_open(trim(roms_filename), nf90_nowrite, ncid), &
              'get_grid_dimensions', 'open '//trim(roms_filename))

Nxi_rho   = get_dimension_length(ncid, 'xi_rho',   roms_filename)
Nxi_u     = get_dimension_length(ncid, 'xi_u',     roms_filename)
Nxi_v     = get_dimension_length(ncid, 'xi_v',     roms_filename)
Neta_rho  = get_dimension_length(ncid, 'eta_rho',  roms_filename)
Neta_u    = get_dimension_length(ncid, 'eta_u',    roms_filename)
Neta_v    = get_dimension_length(ncid, 'eta_v',    roms_filename)

call nc_check(nf90_close(ncid), &
              'get_grid_dimensions','close '//trim(roms_filename))

! Read the vertical dimensions from the dedicated file.

call nc_check(nf90_open(trim(roms_filename), nf90_nowrite, ncid), &
               'get_grid_dimensions', 'open '//trim(roms_filename))

Ns_rho    = get_dimension_length(ncid, 's_rho',    roms_filename)
Ns_w      = get_dimension_length(ncid, 's_w'  ,    roms_filename)

call nc_check(nf90_close(ncid), &
              'get_grid_dimensions','close '//trim(roms_filename))

Nx =  Nxi_rho  ! Setting the nominal value of the 'global' variables
Ny = Neta_rho  ! Setting the nominal value of the 'global' variables
Nz =   Ns_rho  ! Setting the nominal value of the 'global' variables

end subroutine get_grid_dimensions


!-----------------------------------------------------------------------
!>
!> Read the actual grid values from the ROMS netcdf file.
!>
!>@todo FIXME:  the original implementation opened 3 different files
!> to get the grid info - the namelist was:
!>    roms_ini_filename            = '../data/wc13_ini.nc'
!>    grid_definition_filename     = '../data/wc13_grd.nc'
!>    depths_definition_filename   = '../data/wc13_depths.nc'
!>
!> these have been consolidated by hernan for the santa cruz version
!> into a single file.  check with the other rutgers folks to see if
!> they still need to open 3 different files.  if so, we might need
!> to restore the 3 namelist items and we can use the same file for
!> all 3 types of grid info in the first case, and 3 different files
!> for the second case.
!>

subroutine get_grid()

integer  :: ncid, VarID

real(r8), parameter :: all_land = 0.001_r8

if (.not. allocated(ULAT)) allocate(ULAT(Nxi_u, Neta_u))
if (.not. allocated(ULON)) allocate(ULON(Nxi_u, Neta_u))
if (.not. allocated(UDEP)) allocate(UDEP(Nxi_u, Neta_u, Nz))

if (.not. allocated(VLAT)) allocate(VLAT(Nxi_v, Neta_v))
if (.not. allocated(VLON)) allocate(VLON(Nxi_v, Neta_v))
if (.not. allocated(VDEP)) allocate(VDEP(Nxi_v, Neta_v, Nz))

if (.not. allocated(TLAT)) allocate(TLAT(Nxi_rho, Neta_rho))
if (.not. allocated(TLON)) allocate(TLON(Nxi_rho, Neta_rho))
if (.not. allocated(TDEP)) allocate(TDEP(Nxi_rho, Neta_rho, Nz))
if (.not. allocated(WDEP)) allocate(WDEP(Nxi_rho, Neta_rho, Ns_w))

! Read the vertical information from the (separate) roms_filename

call nc_check(nf90_open(trim(roms_filename), nf90_nowrite, ncid), &
      'get_grid', 'open '//trim(roms_filename))

call nc_check(nf90_inq_varid(ncid, 'z_u', VarID), &
      'get_grid', 'inq_varid z_u '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, UDEP), &
      'get_grid', 'get_var z_u '//trim(roms_filename))

call nc_check(nf90_inq_varid(ncid, 'z_w', VarID), &
      'get_grid', 'inq_varid z_w '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, WDEP), &
      'get_grid', 'get_var z_w '//trim(roms_filename))

call nc_check(nf90_inq_varid(ncid, 'z_v', VarID), &
      'get_grid', 'inq_varid z_v '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, VDEP), &
      'get_grid', 'get_var z_v '//trim(roms_filename))

call nc_check(nf90_inq_varid(ncid, 'z_rho', VarID), &
      'get_grid', 'inq_varid z_rho '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, TDEP), &
      'get_grid', 'get_var z_rho '//trim(roms_filename))

call nc_check(nf90_close(ncid), &
             'get_var','close '//trim(roms_filename))

! Read the rest of the grid information from the traditional grid file

call nc_check(nf90_open(trim(roms_filename), nf90_nowrite, ncid), &
      'get_grid', 'open '//trim(roms_filename))

call nc_check(nf90_inq_varid(ncid, 'lon_rho', VarID), &
   'get_grid', 'inq_varid lon_rho '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, TLON), &
      'get_grid', 'get_var lon_rho '//trim(roms_filename))

where (TLON < 0.0_r8) TLON = TLON + 360.0_r8

call nc_check(nf90_inq_varid(ncid, 'lat_rho', VarID), &
      'get_grid', 'inq_varid lat_rho '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, TLAT), &
      'get_grid', 'get_var lat_rho '//trim(roms_filename))

call nc_check(nf90_inq_varid(ncid, 'lon_u', VarID), &
      'get_grid', 'inq_varid lon_u '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, ULON), &
      'get_grid', 'get_var lon_u '//trim(roms_filename))

where (ULON < 0.0_r8) ULON = ULON + 360.0_r8

call nc_check(nf90_inq_varid(ncid, 'lat_u', VarID), &
      'get_grid', 'inq_varid lat_u '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, ULAT), &
      'get_grid', 'get_var lat_u '//trim(roms_filename))

call nc_check(nf90_inq_varid(ncid, 'lon_v', VarID), &
      'get_grid', 'inq_varid lon_v '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, VLON), &
      'get_grid', 'get_var lon_v '//trim(roms_filename))

where (VLON < 0.0_r8) VLON = VLON + 360.0_r8

call nc_check(nf90_inq_varid(ncid, 'lat_v', VarID), &
      'get_grid', 'inq_varid lat_v '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, VLAT), &
      'get_grid', 'get_var lat_v '//trim(roms_filename))

! Be aware that all the depths are negative values.
! The surface of the ocean is 0.0, the deepest is a big negative value.

if (do_output() .and. debug > 0) then
    write(string1,*)'    min/max ULON ',minval(ULON), maxval(ULON)
    write(string2,*)    'min/max ULAT ',minval(ULAT), maxval(ULAT)
    write(string3,*)    'min/max UDEP ',minval(UDEP), maxval(UDEP)
    call error_handler(E_MSG,'get_grid',string1, text2=string2, text3=string3)

    write(string1,*)'    min/max VLON ',minval(VLON), maxval(VLON)
    write(string2,*)    'min/max VLAT ',minval(VLAT), maxval(VLAT)
    write(string3,*)    'min/max VDEP ',minval(VDEP), maxval(VDEP)
    call error_handler(E_MSG,'get_grid',string1, text2=string2, text3=string3)

    write(string1,*)'    min/max TLON ',minval(TLON), maxval(TLON)
    write(string2,*)    'min/max TLAT ',minval(TLAT), maxval(TLAT)
    write(string3,*)    'min/max TDEP ',minval(TDEP), maxval(TDEP)
    call error_handler(E_MSG,'get_grid',string1, text2=string2, text3=string3)
endif

end subroutine get_grid


!-----------------------------------------------------------------------
!>
!> Fill the array of requested variables, dart kinds, possible min/max
!> values and whether or not to update the field in the output file.
!>
!>@param state_variables the list of variables and kinds from model_mod_nml
!>@param ngood the number of variable/KIND pairs specified

subroutine parse_variable_input( state_variables, ngood )

character(len=*), intent(in)  :: state_variables(:)
integer,          intent(out) :: ngood

integer :: i
character(len=NF90_MAX_NAME) :: varname       ! column 1
character(len=NF90_MAX_NAME) :: dartstr       ! column 2
character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
character(len=NF90_MAX_NAME) :: state_or_aux  ! column 5   change to updateable

ngood = 0
MyLoop : do i = 1, MAX_STATE_VARIABLES

   varname      = trim(state_variables(num_state_table_columns*i-4))
   dartstr      = trim(state_variables(num_state_table_columns*i-3))
   minvalstring = trim(state_variables(num_state_table_columns*i-2))
   maxvalstring = trim(state_variables(num_state_table_columns*i-1))
   state_or_aux = trim(state_variables(num_state_table_columns*i  ))

   if ( varname == ' ' .and. dartstr == ' ' ) exit MyLoop ! Found end of list.

   if ( varname == ' ' .or. dartstr == ' ' ) then
      string1 = 'model_nml:model "variables" not fully specified'
      call error_handler(E_ERR,'parse_variable_input:',string1,source,revision,revdate)
   endif

   ! Make sure DART kind is valid

   if( get_index_for_quantity(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'parse_variable_input:',string1,source,revision,revdate)
   endif

   call to_upper(minvalstring)
   call to_upper(maxvalstring)
   call to_upper(state_or_aux)

   var_names(   i) = varname
   kind_list(   i) = get_index_for_quantity(dartstr)
   clamp_vals(i,1) = string_to_real(minvalstring)
   clamp_vals(i,2) = string_to_real(maxvalstring)
   update_list( i) = string_to_logical(state_or_aux, 'UPDATE')

   ngood = ngood + 1

enddo MyLoop

if (ngood == MAX_STATE_VARIABLES) then
   string1 = 'WARNING: There is a possibility you need to increase ''MAX_STATE_VARIABLES'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'parse_variable_input:',string1,source,revision,revdate,text2=string2)
endif

end subroutine parse_variable_input


!-----------------------------------------------------------------------
!
!> Find the named variable (often 'ocean_time') in a ROMS netCDF file.
!> If it is not found, it is a fatal error.
!>
!> @param filename the name of the ROMS netCDF file
!>                 (used to generate useful error messages).
!> @param ncid the netCDF handle to the ROMS netCDF file.
!> @param variable name which contains the time
!> @param calendar the character string indicating the calendar in use
!> @param last_time_index the value of the last time dimension
!> @param last_time the time/date of the last time
!> @param origin_time the base time other times are relative to
!> @param all_times an array of all times in the variable
!>
!>@todo FIXME Make sure the calculation is correct.
!>  A 64bit real can support whole numbers that overflow a 32 bit integer.

subroutine get_time_information(filename, ncid, var_name, dim_name, myvarid, &
                    calendar, last_time_index, last_time, origin_time, all_times)

character(len=*),            intent(in)  :: filename
integer,                     intent(in)  :: ncid
character(len=*),            intent(in)  :: var_name
character(len=*),            intent(in)  :: dim_name
integer,           optional, intent(out) :: myvarid
character(len=32), optional, intent(out) :: calendar
integer,           optional, intent(out) :: last_time_index
type(time_type),   optional, intent(out) :: last_time
type(time_type),   optional, intent(out) :: origin_time
type(time_type),   optional, intent(out) :: all_times(:)

integer :: ios, DimID, VarID, dimlen, i
character(len=64) :: unitstring
character(len=32) :: calendarstring

integer :: year, month, day, hour, minute, second, rc
real(digits12), allocatable :: these_times(:)
type(time_type) :: time_offset, base_time

logical :: offset_in_seconds  ! if .false., assuming offset in days

integer :: original_calendar_type

!>@todo FIXME get the variable length from the varid and remove the need for the dimension name

call nc_check(nf90_inq_dimid(ncid,dim_name,dimid=DimID), &
       'get_time_information','cannot find "'//trim(dim_name)//'" dimension in '//trim(filename))

call nc_check(nf90_inquire_dimension(ncid, DimID, len=dimlen), &
       'get_time_information', 'inquire_dimension '//trim(dim_name)//' from '//trim(filename))
if (present(last_time_index)) last_time_index = dimlen

call nc_check(nf90_inq_varid(ncid, var_name, VarID), &
       'get_time_information', 'inq_varid '//trim(var_name)//' from '//trim(filename))
if (present(myvarid)) myvarid = VarID

! assume gregorian calendar unless there's a calendar attribute saying elsewise
rc = nf90_get_att(ncid, VarID, 'calendar', calendarstring)
if (rc /= nf90_noerr) calendarstring = 'gregorian'
if (present(calendar)) calendar = trim(calendarstring)

if (trim(calendarstring) /= 'gregorian') then
   write(string1,*)'expecting '//trim(var_name)//' calendar of "gregorian"'
   write(string2,*)'got '//trim(calendarstring)
   call error_handler(E_MSG,'get_time_information:', string1, &
             source, revision, revdate, text2=string2, text3=string3)
endif

if (present(last_time) .or. present(origin_time) .or. present(all_times)) then

   ! May need to put the calendar back to some original value
   original_calendar_type = get_calendar_type()

   ! We need to set the calendar to interpret the time values
   ! do we need to preserve the original calendar setting if there is one?
   call set_calendar_type( trim(calendarstring) )

   ! Make sure the calendar is expected form
   ! var_name:units    = "seconds since 1999-01-01 00:00:00" ;
   !                      1234567890123
   !   OR
   ! var_name:units    = "days since 1999-01-01 00:00:00" ;
   !                      1234567890

   call nc_check(nf90_get_att(ncid, VarID, 'units', unitstring), &
          'get_time_information', 'get_att '//trim(var_name)//' units '//trim(filename))

   ! decode the start time of the time variable - expecting time to be coded
   ! as an offset to some base

   if (unitstring(1:13) == 'seconds since') then
      read(unitstring,'(14x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
      if (ios /= 0) then
         write(string1,*)'Unable to read time variable units. Error status was ',ios
         write(string2,*)'expected "seconds since YYYY-MM-DD HH:MM:SS"'
         write(string3,*)'was      "'//trim(unitstring)//'"'
         call error_handler(E_ERR, 'get_time_information:', string1, &
                source, revision, revdate, text2=string2, text3=string3)
      endif
      offset_in_seconds = .true.

   else if (unitstring(1:10) == 'days since') then
      read(unitstring,'(11x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
      if (ios /= 0) then
         write(string1,*)'Unable to read time variable units. Error status was ',ios
         write(string2,*)'expected "days since YYYY-MM-DD HH:MM:SS"'
         write(string3,*)'was      "'//trim(unitstring)//'"'
         call error_handler(E_ERR, 'get_time_information:', string1, &
                source, revision, revdate, text2=string2, text3=string3)
      endif
      offset_in_seconds = .false.

   else
      write(string1,*)'expecting time attribute units of "seconds since ..." -OR-'
      write(string2,*)'                              "days since ..."'
      write(string3,*)'got "'//trim(unitstring)//'"'
      call error_handler(E_ERR,'get_time_information:', string1, &
                source, revision, revdate, text2=string2, text3=string3)
   endif

   base_time = set_date(year, month, day, hour, minute, second)

   if (present(origin_time)) origin_time = base_time

   if (present(last_time) .or. present(all_times)) then

      ! big_integer may overflow a 32bit integer, so declare it 64bit
      ! and parse it into an integer number of days and seconds, both
      ! of which can be 32bit. Our set_time, set_date routines need 32bit integers.

      allocate(these_times(dimlen))

      call nc_check(nf90_get_var( ncid, VarID, these_times), &
             'get_time_information', 'get_var '//trim(var_name)//' from '//trim(filename))

      if (present(last_time)) then
         time_offset = convert_to_time_offset(these_times(dimlen), offset_in_seconds)
         last_time = base_time + time_offset
      endif

      if (present(all_times)) then
         do i=1, dimlen
            time_offset = convert_to_time_offset(these_times(i), offset_in_seconds)
            all_times(i) = base_time + time_offset
         enddo
      endif

      if (do_output() .and. debug > 0 .and. present(last_time)) then
         call print_time(last_time, str='last roms time is ',iunit=logfileunit)
         call print_time(last_time, str='last roms time is ')
         call print_date(last_time, str='last roms date is ',iunit=logfileunit)
         call print_date(last_time, str='last roms date is ')
      endif

      deallocate(these_times)

   endif

   call set_calendar_type(original_calendar_type)

endif

end subroutine get_time_information

!-----------------------------------------------------------------------

!> convert a fractional day to a dart time type

function convert_to_time_offset(offset, offset_in_seconds)

real(digits12), intent(in) :: offset
logical,        intent(in) :: offset_in_seconds
type(time_type) :: convert_to_time_offset

integer(i8) :: big_integer
integer :: some_seconds, some_days

if (offset_in_seconds) then
   big_integer  = int(offset,i8)
   some_days    = big_integer / (24*60*60)
   big_integer  = big_integer - (some_days * (24*60*60))
   some_seconds = int(big_integer,i4)
else
   ! offset in fractional days
   some_days    = int(offset)
   some_seconds = (offset - some_days) * (24*60*60)
endif

convert_to_time_offset = set_time(some_seconds, some_days)

end function convert_to_time_offset

!-----------------------------------------------------------------------
!>
!> convert DART time type into a character string with the
!> format of YYYYMMDDhh ... or DDhh
!>
!> @param time_to_string the character string containing the time
!> @param t the time
!> @param interval logical flag describing if the time is to be
!>                 interpreted as a calendar date or a time increment.
!>                 If the flag is merely present, the time is to be
!>                 interpreted as an increment and the format is simply
!>                 DDhh. If the flag is not present, the time is a full
!>                 calendar (Gregorian) date and will be renedered with
!>                 the YYYYMMDDhh format.
!>

function time_to_string(t, interval)

character(len=19)              :: time_to_string
type(time_type),   intent(in) :: t
logical, optional, intent(in) :: interval

! local variables

integer :: iyear, imonth, iday, ihour, imin, isec
integer :: ndays, nsecs
logical :: dointerval

if (present(interval)) then
   dointerval = interval
else
   dointerval = .false.
endif

! for interval output, output the number of days, then hours, mins, secs
! for date output, use the calendar routine to get the year/month/day hour:min:sec
if (dointerval) then
   call get_time(t, nsecs, ndays)
   if (ndays > 99) then
      write(string1, *) 'interval number of days is ', ndays
      call error_handler(E_ERR,'time_to_string:', 'interval days cannot be > 99', &
                         source, revision, revdate, text2=string1)
   endif
   ihour = nsecs / 3600
   nsecs = nsecs - (ihour * 3600)
   imin  = nsecs / 60
   nsecs = nsecs - (imin * 60)
   isec  = nsecs
!   write(time_to_string, '(I2.2,3(A1,I2.2))') &
!                        ndays, '_', ihour, ':', imin, ':', isec
   write(time_to_string, '(I2.2,I2.2)') &
                        ndays, ihour
else
   call get_date(t, iyear, imonth, iday, ihour, imin, isec)
   write(time_to_string, '(I4.4,5(A1,I2.2))') &
          iyear, '-', imonth, '-', iday, ' ', ihour, ':', imin, ':', isec
endif

end function time_to_string


!-----------------------------------------------------------------------
!>
!> gets the length of a netCDF dimension given the dimension name.
!> This bundles the nf90_inq_dimid and nf90_inquire_dimension routines
!> into a slightly easier-to-use function.
!>
!> @param dimlen the length of the netCDF dimension in question
!> @param ncid the netCDF file handle
!> @param dimension_name the character string of the dimension name
!> @param filename the name of the netCDF file (for error message purposes)
!>

function get_dimension_length(ncid, dimension_name, filename) result(dimlen)

integer                      :: dimlen
integer,          intent(in) :: ncid
character(len=*), intent(in) :: dimension_name
character(len=*), intent(in) :: filename

integer :: DimID

write(string1,*)'inq_dimid '//trim(dimension_name)//' '//trim(filename)
write(string2,*)'inquire_dimension '//trim(dimension_name)//' '//trim(filename)

call nc_check(nf90_inq_dimid(ncid, trim(dimension_name), DimID), &
              'get_dimension_length',string1)
call nc_check(nf90_inquire_dimension(ncid, DimID, len=dimlen), &
              'get_dimension_length', string2)

end function get_dimension_length

!-----------------------------------------------------------------------
!>
!> writes the time of the current state and (optionally) the time
!> to be conveyed to ROMS to dictate the length of the forecast.
!> This file is then used by scripts to modify the ROMS run.
!> The format in the time information is totally at your discretion.
!>
!> @param model_time the current time of the model state
!>

subroutine write_roms_time_information( model_time )
type(time_type), intent(in) :: model_time

integer :: iunit, day, second, ios, ncid
real(digits12) :: dstart
type(time_type) :: base_time, forecast_time

base_time = model_time

call nc_check( nf90_open(trim(roms_filename), NF90_NOWRITE, ncid), &
                  'write_roms_time_information', 'open '//trim(roms_filename))

call get_time_information(roms_filename, ncid, 'ocean_time', 'ocean_time', origin_time=base_time)

call nc_check( nf90_close(ncid), &
                 'write_roms_time_information', 'close '//trim(roms_filename))

forecast_time = model_time - base_time

call get_time(forecast_time,second,day)

dstart = real(day,digits12) + second/86400.0_digits12

iunit = open_file('new_dstart.txt', action='write')

write(iunit,'(A,1x,f25.5)', iostat=ios) 'DSTART ', dstart
if (ios /= 0) then
   write(string1,*)'Unable to write new DSTART. Error status was ',ios
   write(string2,*)'dstart = ', dstart
   call error_handler(E_ERR, 'write_roms_time_information:', string1, &
          source, revision, revdate, text2=string2)
endif

! The next records are totally optional - not checking write status
write(iunit,'(A,1x,f25.5)') 'ROMS_offset ', dstart*86400.0_digits12
call print_date(model_time, str='ROMS_date ',iunit=iunit)
call print_time(model_time, str='DART_time ',iunit=iunit)

string1 = time_to_string(model_time)
write(iunit,'(A,1x,A)') 'YYYYMMDD',trim(string1)

call close_file(iunit)

end subroutine write_roms_time_information


!-----------------------------------------------------------------------
!>
!> Returns the lat,lon,depth given a fractional i,j,k and a specified kind
!>
!> @param filoc fractional x index
!> @param fjloc fractional y index
!> @param fkloc fractional vert index
!> @param dart_kind
!> @param locatiation location at fractional i,j,k
!>
!>  Each grid cell is oriented in a counter clockwise direction
!>  for interpolating locations.  First we interpolate in latitude
!>  and longitude, then interpolate in height.  The hgt of each grid
!>  cell can very on each interpolation, so we have to be careful how
!>  we interpolate in the horizontal.  Using the 4 different heights
!>  and lat_frac, lon_frac, hgt_frac we can do a simple trilinear
!>  interpolation to find the location given fractional indicies.
!>
!>              (i ,j+1) ----- (i+1,j+1)   hgt(1) hgt(2) hgt(3) hgt(4)
!>               hgt(4)         hgt(3)       |      |      |      |
!>                 |               |         |      *      |      |
!>  lon_frac _____ |       X       |         |      |      |      *
!>                 |               |         |             *      |
!>                 |               |         *             |
!>              (i ,j)   ----- (i+1,j)       |
!>               hgt(1)   |     hgt(2)             * - hgt_frac location
!>                        |
!>                     lat_frac
!>
!>    ISTATUS : 10 - bad incoming dart_kind
!>    ISTATUS : 11 - fkloc out of range
!>    ISTATUS : 12 - filoc or fjloc out of range for u grid
!>    ISTATUS : 13 - filoc or fjloc out of range for v grid
!>    ISTATUS : 14 - filoc or fjloc out of range for rho grid
!>    ISTATUS : 99 - initalized istatus, this should not happen

function get_location_from_ijk(filoc, fjloc, fkloc, dart_kind, location) result(istatus)
real(r8),            intent(in)  :: filoc
real(r8),            intent(in)  :: fjloc
real(r8),            intent(in)  :: fkloc
integer,             intent(in)  :: dart_kind
type(location_type), intent(out) :: location
integer :: istatus ! 0 good, else bad

integer  :: var_id, iloc, jloc, vloc
integer  :: vert_type, my_kind
real(r8) :: lon_fract, lat_fract, hgt_fract
real(r8) :: lon_val, lat_val, hgt_val, tmp_hgt(2), hgt(4)
real(r8), pointer :: mylon(:,:), mylat(:,:), mydep(:,:,:)
logical, save :: first_time = .true.


write(string1,*)'Routine not finished.'
call error_handler(E_ERR, 'get_location_from_ijk:', string1, &
                      source, revision, revdate)

! start out assuming bad istatus
istatus  = 99

var_id = get_varid_from_kind(domain_id, dart_kind)
if (var_id < 0) then
  istatus = 10 ! can not find variable id for dart_kind
  return
endif

! check that we have a valid vertical location.
! allow obs above the top rho point but below the
! surface to be considered at the top rho point.
! the commented out code is the test for excluding
! obs above the top rho point in addition to obs
! below the bottom rho point.
!if (fkloc < 0.5 .or. fkloc > Ns_rho - 0.5) then
if (fkloc < 0.5 .or. fkloc > Ns_rho) then
  istatus = 11
  location = set_location_missing()
  return
endif

iloc = FLOOR(filoc)
jloc = FLOOR(fjloc)
vloc = FLOOR(fkloc)

lon_fract = filoc - iloc
lat_fract = fjloc - jloc
hgt_fract = fkloc - vloc

my_kind = get_kind_index(domain_id,var_id)
if (my_kind==QTY_U_CURRENT_COMPONENT) then
   write(string1,*)'Not interpolating ', get_name_for_quantity(my_kind), ' at the moment.'
   write(string2,*)'Need to check that we are using the right grid for location interpolation'
   call error_handler(E_ERR, 'get_location_from_ijk:', string1, &
                      source, revision, revdate, text2=string2)
   if (filoc < 1 .or. filoc > Nxi_u-1 .or. &
       fjloc < 1 .or. fjloc > Neta_u-1 ) then
     istatus = 12
     location = set_location_missing()
     return
   endif
   mylon => ULON
   mylat => ULAT
   mydep => UDEP
elseif (my_kind==QTY_V_CURRENT_COMPONENT) then
   write(string1,*)'Not interpolating ', get_name_for_quantity(my_kind), ' at the moment.'
   write(string2,*)'Need to check that we are using the right grid for location interpolation'
   call error_handler(E_ERR, 'get_location_from_ijk:', string1, &
                      source, revision, revdate, text2=string2)
   if (filoc < 1 .or. filoc > Nxi_v-1 .or. &
       fjloc < 1 .or. fjloc > Neta_v-1 ) then
     istatus = 13
     location = set_location_missing()
     return
   endif
   mylon => VLON
   mylat => VLAT
   mydep => VDEP
else  ! Everything else is assumed to be on the rho points
   if (filoc < 1 .or. filoc > Nxi_rho-1 .or. &
       fjloc < 1 .or. fjloc > Neta_rho-1 ) then

     write(*,*)
     write(*,*)'filoc, Nxi_rho-1       = ',filoc, Nxi_rho-1
     write(*,*)'fjloc, Neta_rho-1      = ',fjloc, Neta_rho-1
     write(*,*)'fkloc, vloc, hgt_fract = ',fkloc,vloc,hgt_fract

     write(logfileunit,*)
     write(logfileunit,*)'filoc, Nxi_rho-1     = ',filoc, Nxi_rho-1
     write(logfileunit,*)'fjloc, Neta_rho-1    = ',fjloc, Neta_rho-1
     write(logfileunit,*)'fkloc,vloc,hgt_fract = ',fkloc,vloc,hgt_fract

     istatus = 14
     location = set_location_missing()
     return
   endif
   mylon => TLON
   mylat => TLAT
   mydep => TDEP
endif

lon_val = (1.0-lon_fract)*(mylon(iloc,jloc)) + (lon_fract)*(mylon(iloc+1,jloc))
lat_val = (1.0-lat_fract)*(mylat(iloc,jloc)) + (lat_fract)*(mylat(iloc  ,jloc+1))

if( get_kind_index(domain_id,var_id) == QTY_SEA_SURFACE_HEIGHT ) then
   hgt_val    = 0.0_r8
   vert_type  = VERTISSURFACE
else
   if (fkloc > Ns_rho - 0.5) then
      ! leave something in the log to remind us we're doing this and
      ! make sure it's ok with everyone.
      if (first_time) then
         call error_handler(E_MSG, 'ROMS model_mod', &
           'NOTE: Locations above the top Rho point are moved to that point')
         first_time = .false.
      endif
      ! special case for obs at or closer to surface than top rho point
      ! make them 100% at the top point
      vloc = vloc - 1
      hgt_fract = 1.0_r8
   endif

   ! fractional heights in the vertical for each corner of the horizontal box
   hgt(1) = (1.0-hgt_fract)*(mydep(iloc  ,jloc,  vloc)) + hgt_fract*(mydep(iloc  ,jloc  ,vloc+1))
   hgt(2) = (1.0-hgt_fract)*(mydep(iloc+1,jloc  ,vloc)) + hgt_fract*(mydep(iloc+1,jloc  ,vloc+1))
   hgt(3) = (1.0-hgt_fract)*(mydep(iloc+1,jloc+1,vloc)) + hgt_fract*(mydep(iloc+1,jloc+1,vloc+1))
   hgt(4) = (1.0-hgt_fract)*(mydep(iloc  ,jloc+1,vloc)) + hgt_fract*(mydep(iloc  ,jloc+1,vloc+1))

   tmp_hgt(1) = (1.0-lon_fract)*hgt(1) + lon_fract*hgt(2)
   tmp_hgt(2) = (1.0-lon_fract)*hgt(4) + lon_fract*hgt(3)

   hgt_val   = (1.0-lat_fract)*tmp_hgt(1)+lat_fract*tmp_hgt(2)
   vert_type = VERTISHEIGHT
endif

istatus  = 0
location = set_location(lon_val, lat_val, hgt_val, vert_type)

if (debug > 5) then
   print*,' i,j,k', filoc, fjloc, fkloc

   print*,' lon(i  ,j  ), lat(i  ,j  ) : (', mylon(iloc  ,jloc  ),',',mylat(iloc  ,jloc  ),')'
   print*,' lon(i+1,j  ), lat(i+1,j  ) : (', mylon(iloc+1,jloc  ),',',mylat(iloc+1,jloc  ),')'
   print*,' lon(i+1,j+1), lat(i+1,j+1) : (', mylon(iloc+1,jloc+1),',',mylat(iloc+1,jloc+1),')'
   print*,' lon(i  ,j+1), lat(i  ,j+1) : (', mylon(iloc  ,jloc+1),',',mylon(iloc  ,jloc+1),')'
   print*,' '
   print*,' tmp_hgt(1)  = ', tmp_hgt(1)
   print*,' tmp_hgt(2)  = ', tmp_hgt(2)
   print*, ' '
   print*,' lon_frac    = ', lon_fract
   print*,' lat_frac    = ', lat_fract
   print*,' hgt_frac    = ', hgt_fract
   print*,' '
   print*,' lon_val     = ', lon_val
   print*,' lat_val     = ', lat_val
   print*,' hgt_val     = ', hgt_val
   print*,' '
   print*,' WDEP(i  , j  , k  )', WDEP(iloc,   jloc  , vloc)
   print*,' WDEP(i  , j+1, k  )', WDEP(iloc,   jloc+1, vloc)
   print*,' WDEP(i+1, j+1, k  )', WDEP(iloc+1, jloc+1, vloc)
   print*,' WDEP(i+1, j  , k  )', WDEP(iloc+1, jloc  , vloc)
   print*,' '
   print*,' WDEP(i  , j  , k+1)', WDEP(iloc,   jloc  , vloc+1)
   print*,' WDEP(i  , j+1, k+1)', WDEP(iloc,   jloc+1, vloc+1)
   print*,' WDEP(i+1, j+1, k+1)', WDEP(iloc+1, jloc+1, vloc+1)
   print*,' WDEP(i+1, j  , k+1)', WDEP(iloc+1, jloc  , vloc+1)
endif

end function get_location_from_ijk


!===================================================================
! End of model_mod
!===================================================================

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
