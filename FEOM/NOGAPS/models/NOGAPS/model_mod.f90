! This code may (or may not) be part of the NOGAPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

module model_mod

! This is a template showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. The public interfaces listed
! must all be supported with the argument lists as indicated. Many of the interfaces
! are not required for minimal implementation (see the discussion of each
! interface and look for NULL INTERFACE). 

! Modules that are absolutely required for use are listed
use        types_mod, only : r8, MISSING_R8, PI, digits12, RAD2DEG, DEG2RAD
use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                             print_time, print_date, set_calendar_type, GREGORIAN, &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=)
use     location_mod, only : location_type, get_dist, get_close_maxdist_init, &
                             get_close_obs_init, set_location, query_location,&
                             VERTISHEIGHT, VERTISSURFACE, VERTISLEVEL,        &  
                             VERTISPRESSURE, VERTISUNDEF, &  
                             vert_is_height, vert_is_level, vert_is_surface,  &
                             get_location, loc_get_close_obs => get_close_obs, &
                             get_close_type, vert_is_pressure
use    utilities_mod, only : register_module, error_handler, nc_check, &
                             E_ERR, E_MSG, open_file, file_exist, logfileunit, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, check_namelist_read, &
                             find_textfile_dims, file_to_text
use     obs_kind_mod, only : KIND_TEMPERATURE,      KIND_U_WIND_COMPONENT,     &
                             KIND_V_WIND_COMPONENT, KIND_SPECIFIC_HUMIDITY,    &
                             KIND_PRESSURE

use typesizes
use netcdf 

implicit none
private

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

public :: test_interpolate

! These interfaces are useful, but not required.

public :: get_gridsize, &
          pres_sounding

! This is from the original assembla server we used during collaboration.
! character(len=128), parameter :: &
!    source   = "$orgURL: https://svn2.assembla.com/svn/ngdart/model_mod.f90 $", &
!    revision = "$orgRevision: 113 $", &
!    revdate  = "$orgDate: 2010-06-11 14:49:52 -0600 (Fri, 11 Jun 2010) $"

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


character(len=256) :: msgstring
logical, save :: module_initialized = .false.

integer, parameter :: NOT_INITIALIZED = -1

! define model parameters here
integer         :: model_size = NOT_INITIALIZED
type(time_type) :: time_step

! define model namelist here 
logical  :: output_state_vector = .true.
integer  :: time_step_days      = 0
integer  :: time_step_seconds   = 900
character(len=128) :: geometry_text_file = "noggeom.txt"
integer  :: debug = 0   ! turn up for more and more debug messages

namelist /model_nml/ output_state_vector, time_step_days, time_step_seconds, &
                     geometry_text_file, debug

! Global storage for prognostic variables
integer, parameter :: n3dfields = 4
integer, parameter :: n2dfields = 1
integer, parameter :: nfields   = n3dfields + n2dfields

! (the absoft compiler likes them to all be the same length during declaration)
! we trim the blanks off before use anyway, so ...
character(len=128) :: progvarnames(nfields) = (/'UVEL ','VVEL ','TEMP ',&
                                                'SHUM ','PSURF'/)

integer, parameter :: U_index     = 1
integer, parameter :: V_index     = 2
integer, parameter :: T_index     = 3
integer, parameter :: Q_index     = 4
integer, parameter :: PSURF_index = 5
integer            :: start_index(nfields)

integer :: numlons = NOT_INITIALIZED
integer :: numlats = NOT_INITIALIZED
integer :: numlevs = NOT_INITIALIZED

real(r8), allocatable :: asig(:), bsig(:), sinl(:)
real(r8), allocatable :: lons(:)
real(r8), allocatable :: lats(:)
real(r8), allocatable :: levs(:)

real(r8), allocatable :: ensemble_mean(:)

! locations of cell centers (C) and edges (G) for each axis.
! These arrays store the longitude and latitude arrays. 
real(r8), allocatable :: zc(:), zg(:)

! Surface Geopotential (m^2/s^2)
! To convert to height, divide phis by the gravitational constant for
! whatever planet you are on.
real(r8), allocatable :: phis(:,:)

INTERFACE vector_to_prog_var
      MODULE PROCEDURE vector_to_2d_prog_var
      MODULE PROCEDURE vector_to_3d_prog_var
END INTERFACE

! More unit conversion fun!
real(r8), parameter :: HPA_TO_PA = 100.0_r8
contains

!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------

integer  :: i, j, k, iunit, io, kx, jx
real(r8) :: dx

if ( module_initialized ) return ! only need to do this once.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! If you're going to use obs from this planet, you're 
! going to use this calendar ...
call set_calendar_type(GREGORIAN)

! NOGAPS has a straight ascii file with the geometry information.
! The model_nml knows the name.

iunit = open_file( geometry_text_file, form="formatted", action="read" )

! This block comes right from the "geometry_text_file" itself. 
! NOGAPS not only writes out the file, it provides code on how
! to read the file AFTER the data ... concise.

! Read dimensions
read(iunit,'(a80)')
read(iunit,'(3i10)') numlons,numlats,numlevs
read(iunit,'(a80)')

! Create storage for locations
allocate(asig(numlevs+1), bsig(numlevs+1), levs(numlevs))
allocate(lons(numlons), lats(numlats), sinl(numlats))
allocate(phis(numlons,numlats))

! Finish reading location information
read(iunit,'(i5,2f16.12)') (kx,asig(k),bsig(k),k=1,numlevs+1)
read(iunit,'(a80)')
read(iunit,'(i5,5x,f16.12)') (jx,sinl(j),j=numlats,1,-1 )

! Now that we know the dimensions, 
! that the longitudes always start at 0.0, 
! that the latitudes start at the south pole, and
! that the level == 1 is at the top of the atmosphere ...
! we can create coordinate variables - in radians.
dx = (2.0_r8 * PI)/numlons
do i = 1,numlons
   lons(i) = (i - 1.0_r8) * dx 
enddo

do i = 1,numlats
   lats(i) = asin(sinl(i))
enddo

do i = 1,numlevs
   levs(i) = real(i,r8)
enddo

!---------------------------------------------------------------
! compute the offsets into the state vector for the start of each
! different variable type.

! record where in the state vector the data type changes
! from one type to another, by computing the starting
! index for each block of data.
start_index(U_index)     = 1
start_index(V_index)     = start_index(U_index) + (numlons * numlats * numlevs)
start_index(T_index)     = start_index(V_index) + (numlons * numlats * numlevs)
start_index(Q_index)     = start_index(T_index) + (numlons * numlats * numlevs)
start_index(PSURF_index) = start_index(Q_index) + (numlons * numlats * numlevs)

if (do_output()) write(logfileunit, *) 'Using grid : numlons, numlats, numlevs = ', &
                                                     numlons, numlats, numlevs
if (do_output()) write(     *     , *) 'Using grid : numlons, numlats, numlevs = ', &
                                                     numlons, numlats, numlevs

model_size = (n3dfields * (numlons * numlats * numlevs)) + &
             (n2dfields * (numlons * numlats))

if (do_output()) write(*,*) 'model_size = ', model_size

allocate(ensemble_mean(model_size))

! The time_step in terms of a time type must also be initialized.
time_step = set_time(time_step_seconds, time_step_days)

! Read in the surface geopotential from some (namelist-specified?) file ... 
! FIXME: netcdf open ... read ... close ...

phis = 1.125_r8    ! just to give me something to write

! Create a netCDF file of what we've got so far.
! if (debug > 0) call write_grid_netcdf()
! if (debug > 0) call write_grid_interptest()

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

x = MISSING_R8

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

! FIXME: for now, just set to 0 ... is this a good idea?
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

use nogaps_interp_mod, only : compute_neighbors, &
                              NUM_NEIGHBORS,     &
                              get_val_at_pressure  

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

integer, parameter :: STATUS_WTF     = 99
integer, parameter :: STATUS_SUCCESS = 0

real(r8) :: loc_array(3)
real(r8) :: obs_lon, obs_lat, obs_lev

integer  :: starting_offset

integer  :: p_surf_position
real(r8) :: p_surf(NUM_NEIGHBORS)
real(r8) :: pressure_profile(numlevs)

integer  :: neighbor_i(NUM_NEIGHBORS)
integer  :: neighbor_j(NUM_NEIGHBORS)
real(r8) :: neighbor_weights(NUM_NEIGHBORS)
real(r8) :: neighbor_vals(NUM_NEIGHBORS)
integer  :: neighbor_indices(numlevs)
integer  :: cur_neighbor

integer :: i

logical :: can_interpolate

! Default for successful return
istatus = STATUS_SUCCESS
obs_val = MISSING_R8

can_interpolate = (vert_is_surface(location) .and. itype == KIND_PRESSURE)
can_interpolate = can_interpolate .or. vert_is_pressure(location)
if (.not. can_interpolate) then
   call error_handler(E_ERR, 'model_interpolate', 'Can only interpolate' // &
                      'obs on pressure levels.', source, revision, revdate)
end if

! Break down the location into components
loc_array = get_location(location)
obs_lon   = loc_array(1) * DEG2RAD
obs_lat   = loc_array(2) * DEG2RAD
obs_lev   = loc_array(3)

if (debug > 1) print *, 'requesting interpolation at ', obs_lon, obs_lat, obs_lev

if (obs_lat > maxval(lats) .or.  obs_lat < minval(lats)) then
    istatus = 128
    return
end if

call compute_neighbors(obs_lat, obs_lon, lats, lons, numlats, numlons, &
                       neighbor_i, neighbor_j, neighbor_weights)

! Find out where the target variable is in the state vector
starting_offset = get_start_index_by_kind(itype)

do cur_neighbor = 1, NUM_NEIGHBORS
    p_surf_position = (start_index(PSURF_index)-1) +                    &
                      index_2d_to_1d(neighbor_i(cur_neighbor),          &
                                     neighbor_j(cur_neighbor), numlons)
    p_surf(cur_neighbor) = x(p_surf_position)

    ! Surface pressure is stored in hPa - observations are in Pa
    if (obs_lev > (p_surf(cur_neighbor)*HPA_TO_PA)) then
        ! Below ground
        istatus = 130
        return
    end if
end do

if (vert_is_surface(location)) then
    ! Assume that the only thing at the surface is pressure (for now)
    neighbor_vals(:) = p_surf(:)
else
    do cur_neighbor = 1, NUM_NEIGHBORS
        call get_column_indices(neighbor_i(cur_neighbor),  &
                                neighbor_j(cur_neighbor),  &
                                numlons, numlats, numlevs, &
                                neighbor_indices) 
        neighbor_indices(:) = neighbor_indices(:) + (starting_offset-1)

        ! surface pressure (hPa) -> column pressure profile (Pa)
        call pres_sounding(numlevs, asig, bsig, p_surf(cur_neighbor), &
                           pressure_profile)
        call get_val_at_pressure(x(neighbor_indices), pressure_profile, &
                                 obs_lev, neighbor_vals(cur_neighbor),  &
                                 istatus)
    end do
end if

obs_val = dot_product(neighbor_weights, neighbor_vals)
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
integer, optional,   intent(out) :: var_type

real(r8) :: lat, lon, level
integer  :: lon_index, lat_index, level_index, local_var

call get_state_indices(index_in, lat_index, lon_index, level_index, local_var)

if (debug > 5) print *, 'g_s_metadata: lonindx, latindx, lvlindx = ', lon_index, lat_index, level_index

lon = lons(lon_index)*RAD2DEG
lat = lats(lat_index)*RAD2DEG

if (local_var == KIND_PRESSURE) then
   ! This is the surface pressure field; vertical coordinate is surface
   level = 0.0_r8
   location = set_location(lon, lat, level, VERTISSURFACE)
else
   ! This is a 3-dimensional field; vertical coordinate is model level
   level = levs(level_index)
   location = set_location(lon, lat, level, VERTISLEVEL)
endif

if (debug > 5) print *, 'g_s_metadata: lon, lat, level = ', lon, lat, level

if (present(var_type)) then
   var_type = local_var
endif

end subroutine get_state_meta_data


subroutine get_state_indices(index_in, lat_index, lon_index, level_index, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated array indices for lat, lon, and level, as well as the type.

integer, intent(in)  :: index_in
integer, intent(out) :: lat_index, lon_index, level_index
integer, intent(out) :: var_type

integer :: startind, offset

if ( .not. module_initialized ) call static_init_model

if (debug > 5) print *, 'g_s_indices:asking for meta data about index ', index_in

call get_state_kind(index_in, var_type, startind, offset)

if (debug > 5) print *, 'g_s_indices: var_type startind = ', var_type, startind 

if (startind == start_index(PSURF_index)) then
  level_index = 1
else
  level_index = (offset / (numlons * numlats)) + 1
endif

lat_index = (offset - ((level_index-1)*numlons*numlats)) / numlons + 1
lon_index =  offset - ((level_index-1)*numlons*numlats) - ((lat_index-1)*numlons) + 1

if (debug > 5) print *, 'g_s_indices:lon, lat, height index = ', lon_index, lat_index, level_index

end subroutine get_state_indices


subroutine get_state_kind(index_in, var_type, startind, offset)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the kind,
! and both the starting offset for this kind, as well as the offset into
! the block of this kind.

integer, intent(in)  :: index_in
integer, intent(out) :: var_type, startind, offset


if ( .not. module_initialized ) call static_init_model

if (debug > 5) print *, 'g_s_kind:asking for meta data about index ', index_in

if (index_in < start_index(U_index+1)) then
   var_type = KIND_U_WIND_COMPONENT
   startind = start_index(U_index)

else if (index_in < start_index(V_index+1)) then
   var_type = KIND_V_WIND_COMPONENT
   startind = start_index(V_index)

else if (index_in < start_index(T_index+1)) then
   var_type = KIND_TEMPERATURE
   startind = start_index(T_index)

else if (index_in < start_index(Q_index+1)) then
   var_type = KIND_SPECIFIC_HUMIDITY
   startind = start_index(Q_index)

else 
   var_type = KIND_PRESSURE
   startind = start_index(PSURF_index)
endif

! local offset into this var array
offset = index_in - startind

if (debug > 5) print *, 'var type = ', var_type
if (debug > 5) print *, 'startind = ', startind
if (debug > 5) print *, 'offset = ', offset

end subroutine get_state_kind



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

! if ( .not. module_initialized ) call static_init_model

deallocate(asig, bsig, sinl, lons, lats, levs, phis, ensemble_mean)

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
integer :: NlonDimID, NlatDimID, NzDimID, NzeDimID
integer :: lonVarID, latVarID, levVarID, levedgeVarID
integer :: sinlVarID, asigVarID, bsigVarID
integer :: PHISVarID

! for the prognostic variables
integer :: UVarID, VVarID, TVarID, QVarID, PSURFVarID

!----------------------------------------------------------------------
! variables for the namelist output
!----------------------------------------------------------------------

character(len=129), allocatable, dimension(:) :: textblock
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen
logical :: has_nogaps_namelist

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
   write(msgstring,*)'Time Dimension ID ',TimeDimID, &
             ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', msgstring, source, revision, revdate)
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
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'NOGAPS' ), &
           'nc_write_model_atts', 'model put '//trim(filename))

!-------------------------------------------------------------------------------
! Determine shape of most important namelist
!-------------------------------------------------------------------------------

call find_textfile_dims('nogaps_in', nlines, linelen)
if (nlines > 0) then
  has_nogaps_namelist = .true.
else
  has_nogaps_namelist = .false.
endif

if (debug > 1) print *, 'NOGAPS namelist: nlines, linelen = ', nlines, linelen
  
if (has_nogaps_namelist) then 
   allocate(textblock(nlines))
   textblock = ''
   
   call nc_check(nf90_def_dim(ncid=ncFileID, name='nlines', &
                 len = nlines, dimid = nlinesDimID), &
                 'nc_write_model_atts', 'def_dim nlines ')
   
   call nc_check(nf90_def_var(ncFileID,name='nogaps_in', xtype=nf90_char,    &
                 dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID), &
                 'nc_write_model_atts', 'def_var nogaps_in')
   call nc_check(nf90_put_att(ncFileID, nmlVarID, 'long_name',       &
                 'contents of nogaps_in namelist'), 'nc_write_model_atts', 'put_att nogaps_in')

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
   
   call nc_check(nf90_def_dim(ncid=ncFileID, name='lon', &
          len = numlons, dimid = NlonDimID),'nc_write_model_atts', 'lon def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='lat', &
          len = numlats, dimid = NlatDimID),'nc_write_model_atts', 'lat def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='lev', &
          len = numlevs, dimid =   NzDimID),'nc_write_model_atts', 'lev def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='leveledges', &
          len = numlevs+1,dimid = NzeDimID),'nc_write_model_atts', 'edges def_dim '//trim(filename))
   
   !----------------------------------------------------------------------------
   ! Create the (empty) Coordinate Variables and the Attributes
   !----------------------------------------------------------------------------

if (has_nogaps_namelist) then 
   call nc_check(nf90_def_var(ncFileID,name='nogapsnml', xtype=nf90_char,    &
                 dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID), &
                 'nc_write_model_atts', 'def_var nogapsnml')
   call nc_check(nf90_put_att(ncFileID, nmlVarID, 'long_name',       &
                 'namelist.input contents'), 'nc_write_model_atts', 'put_att nogapsnml')
endif

   ! Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='lon', xtype=nf90_real, &
                 dimids= NlonDimID, varid=lonVarID),&
                 'nc_write_model_atts', 'lon def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  lonVarID, 'long_name', 'longitudes of grid'), &
                 'nc_write_model_atts', 'lon long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  lonVarID, 'cartesian_axis', 'X'),  &
                 'nc_write_model_atts', 'lon cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  lonVarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'lon units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  lonVarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'lon valid_range '//trim(filename))

   ! Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='lat', xtype=nf90_real, &
                 dimids=NlatDimID, varid=latVarID),&
                 'nc_write_model_atts', 'lat def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  latVarID, 'long_name', 'latitudes of grid'), &
                 'nc_write_model_atts', 'lat long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  latVarID, 'cartesian_axis', 'Y'),   &
                 'nc_write_model_atts', 'lat cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  latVarID, 'units', 'degrees_north'),  &
                 'nc_write_model_atts', 'lat units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  latVarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'lat valid_range '//trim(filename))

   ! Vertical - centers
   call nc_check(nf90_def_var(ncFileID,name='lev', xtype=nf90_real, &
                 dimids=NzDimID, varid= levVarID), &
                 'nc_write_model_atts', 'lev def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, levVarID, 'long_name', 'centers of vertical levels'), &
                 'nc_write_model_atts', 'lev long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, levVarID, 'cartesian_axis', 'Z'),   &
                 'nc_write_model_atts', 'lev cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, levVarID, 'units', 'unknown'),  &
                 'nc_write_model_atts', 'lev units '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID, levVarID, 'positive', 'down'),  &
!                'nc_write_model_atts', 'lev units '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID, levVarID, 'comment', &
!                 'more positive is closer to the center of the earth'),  &
!                'nc_write_model_atts', 'lev comment '//trim(filename))

   ! asig
   call nc_check(nf90_def_var(ncFileID,name='asig', xtype=nf90_real, &
                 dimids=NzeDimID, varid=asigVarID),&
                 'nc_write_model_atts', 'asig def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  asigVarID, 'long_name', 'first coefficient of hybrid vertical coordinates'), &
                 'nc_write_model_atts', 'asig long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  asigVarID, 'formula', 'p(k) = asig(k+1) + bsig(k+1)*ps'),   &
                 'nc_write_model_atts', 'asig formula '//trim(filename))

   ! bsig
   call nc_check(nf90_def_var(ncFileID,name='bsig', xtype=nf90_real, &
                 dimids=NzeDimID, varid=bsigVarID),&
                 'nc_write_model_atts', 'bsig def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  bsigVarID, 'long_name', 'second coefficient of hybrid vertical coordinates'), &
                 'nc_write_model_atts', 'bsig long_name '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID,  bsigVarID, 'cartesian_axis', 'Y'),   &
!                'nc_write_model_atts', 'bsig cartesian_axis '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID,  bsigVarID, 'units', 'degrees_north'),  &
!                'nc_write_model_atts', 'bsig units '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID,  bsigVarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
!                'nc_write_model_atts', 'bsig valid_range '//trim(filename))

   ! sinl
   call nc_check(nf90_def_var(ncFileID,name='sinl', xtype=nf90_real, &
                 dimids=NlatDimID, varid=sinlVarID),&
                 'nc_write_model_atts', 'sinl def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  sinlVarID, 'long_name', 'sine of gaussian latitudes'), &
                 'nc_write_model_atts', 'sinl long_name '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID,  sinlVarID, 'cartesian_axis', 'Y'),   &
!                'nc_write_model_atts', 'sinl cartesian_axis '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID,  sinlVarID, 'units', 'degrees_north'),  &
!                'nc_write_model_atts', 'sinl units '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID,  sinlVarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
!                'nc_write_model_atts', 'sinl vald_range '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncid=ncFileID, name='U', xtype=nf90_real, &
         dimids=(/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=UVarID),&
         'nc_write_model_atts', 'U def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, 'long_name', 'U velocity'), &
         'nc_write_model_atts', 'U long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, 'units', 'm/s'), &
         'nc_write_model_atts', 'U units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, 'units_long_name', 'meters per second'), &
         'nc_write_model_atts', 'U units_long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, 'missing_value', NF90_FILL_REAL), &
         'nc_write_model_atts', 'U missing '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, '_FillValue', NF90_FILL_REAL), &
         'nc_write_model_atts', 'U fill '//trim(filename))

   call nc_check(nf90_def_var(ncid=ncFileID, name='V', xtype=nf90_real, &
         dimids=(/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=VVarID),&
         'nc_write_model_atts', 'V def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, 'long_name', 'V Velocity'), &
         'nc_write_model_atts', 'V long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, 'units', 'm/s'), &
         'nc_write_model_atts', 'V units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, 'units_long_name', 'meters per second'), &
         'nc_write_model_atts', 'V units_long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, 'missing_value', NF90_FILL_REAL), &
         'nc_write_model_atts', 'V missing '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, '_FillValue', NF90_FILL_REAL), &
         'nc_write_model_atts', 'V fill '//trim(filename))

   call nc_check(nf90_def_var(ncid=ncFileID, name='T', xtype=nf90_real, &
         dimids=(/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=TVarID),&
         'nc_write_model_atts', 'T def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, 'long_name', 'Temperature'), &
         'nc_write_model_atts', 'T long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, 'units', 'deg K'), &
         'nc_write_model_atts', 'T units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, 'units_long_name', 'degrees Kelvin'), &
         'nc_write_model_atts', 'T units_long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, 'missing_value', NF90_FILL_REAL), &
         'nc_write_model_atts', 'T missing '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, '_FillValue', NF90_FILL_REAL), &
         'nc_write_model_atts', 'T fill '//trim(filename))

   call nc_check(nf90_def_var(ncid=ncFileID, name='Q', xtype=nf90_real, &
         dimids = (/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=QVarID),&
         'nc_write_model_atts', 'Q def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, QVarID, 'long_name', 'specific humidity'), &
         'nc_write_model_atts', 'Q long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, QVarID, 'units', 'g/g'), &
         'nc_write_model_atts', 'Q units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, QVarID, 'missing_value', NF90_FILL_REAL), &
         'nc_write_model_atts', 'Q missing '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, QVarID, '_FillValue', NF90_FILL_REAL), &
         'nc_write_model_atts', 'Q fill '//trim(filename))

   call nc_check(nf90_def_var(ncid=ncFileID, name='PSURF', xtype=nf90_real, &
         dimids=(/NlonDimID,NlatDimID,MemberDimID,unlimitedDimID/),varid=PSURFVarID), &
         'nc_write_model_atts', 'PSURF def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, PSURFVarID, 'long_name', 'surface pressure'), &
         'nc_write_model_atts', 'PSURF long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, PSURFVarID, 'units', 'hPa'), &
         'nc_write_model_atts', 'PSURF units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, PSURFVarID, 'missing_value', NF90_FILL_REAL), &
         'nc_write_model_atts', 'PSURF missing '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, PSURFVarID, '_FillValue', NF90_FILL_REAL), &
         'nc_write_model_atts', 'PSURF fill '//trim(filename))

   ! Finished with dimension/variable definitions, must end 'define' mode to fill.

   call nc_check(nf90_enddef(ncfileID), 'prognostic enddef '//trim(filename))

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables & module storage static arrays
   !----------------------------------------------------------------------------

   call nc_check(nf90_put_var(ncFileID, lonVarID, lons*RAD2DEG ), &
                'nc_write_model_atts', 'lon put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, latVarID, lats*RAD2DEG ), &
                'nc_write_model_atts', 'lat put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, levVarID, levs ), &
                'nc_write_model_atts', 'lev put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, sinlVarID, sinl ), &
                'nc_write_model_atts', 'sinl put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, asigVarID, asig ), &
                'nc_write_model_atts', 'asig put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, bsigVarID, bsig ), &
                'nc_write_model_atts', 'bsig put_var '//trim(filename))

endif

!-------------------------------------------------------------------------------
! Fill the variables we can
!-------------------------------------------------------------------------------

if (has_nogaps_namelist) then
   call file_to_text('nogaps_in', textblock)
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

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: VarID

real(r8), dimension(numlons,numlats,numlevs) :: data_3d
real(r8), dimension(numlons,numlats)    :: data_2d
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
   ! Staggered grid causes some logistical problems.
   !----------------------------------------------------------------------------

   call vector_to_prog_var(statevec, U_index, data_3d)
   where (data_3d == 0.0_r8) data_3d = NF90_FILL_REAL
   call nc_check(NF90_inq_varid(ncFileID, 'U', VarID), &
                'nc_write_model_vars', 'U inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                'nc_write_model_vars', 'U put_var '//trim(filename))

   call vector_to_prog_var(statevec, V_index, data_3d)
   where (data_3d == 0.0_r8) data_3d = NF90_FILL_REAL
   call nc_check(NF90_inq_varid(ncFileID, 'V', VarID), &
                'nc_write_model_vars', 'V inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                'nc_write_model_vars', 'V put_var '//trim(filename))

   call vector_to_prog_var(statevec, T_index, data_3d)
   where (data_3d == 0.0_r8) data_3d = NF90_FILL_REAL
   call nc_check(NF90_inq_varid(ncFileID, 'T', VarID), &
                'nc_write_model_vars', 'T inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                'nc_write_model_vars', 'T put_var '//trim(filename))

   call vector_to_prog_var(statevec, Q_index, data_3d)
   where (data_3d == 0.0_r8) data_3d = NF90_FILL_REAL
   call nc_check(NF90_inq_varid(ncFileID, 'Q', VarID), &
                'nc_write_model_vars', 'Q inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                'nc_write_model_vars', 'Q put_var '//trim(filename))

   call vector_to_prog_var(statevec, PSURF_index, data_2d)
   where (data_2d == 0.0_r8) data_2d = NF90_FILL_REAL
   call nc_check(NF90_inq_varid(ncFileID, 'PSURF', VarID), &
                'nc_write_model_vars', 'PSURF inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,data_2d,start=(/1,1,copyindex,timeindex/)),&
                'nc_write_model_vars', 'PSURF put_var '//trim(filename))

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
! may do so by adding an O(0.1) magnitude perturbation to each
! model state variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

pert_state      = MISSING_R8
interf_provided = .false.

end subroutine pert_model_state




subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! If needed by the model interface, this is the current mean
! for all state vector items across all ensembles. It is up to this
! code to allocate space and save a copy if it is going to be used
! later on.  For now, we are ignoring it.

real(r8), intent(in) :: ens_mean(:)

if ( .not. module_initialized ) call static_init_model

! Copy into module-level storage
ensemble_mean(:) = ens_mean(:)

end subroutine ens_mean_for_model



subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, &
                         obs_kind, num_close, close_ind, dist)
implicit none

type(get_close_type), intent(in)  :: gc
type(location_type),  intent(in)  :: base_obs_loc, obs_loc(:)
integer,              intent(in)  :: base_obs_kind, obs_kind(:)
integer,              intent(out) :: num_close, close_ind(:)
real(r8),             intent(out) :: dist(:)

! remove some (unused) variables?
integer                :: k, t_ind
integer                :: base_which, local_base_which, obs_which, local_obs_which
real(r8), dimension(3) :: base_array, local_base_array, obs_array, local_obs_array
real(r8)               :: increment, threshold, thresh_wght
type(location_type)    :: local_base_obs_loc, local_obs_loc

real(r8), parameter :: HIGHEST_STATE_PRESSURE = 15000.0_r8 !Pa
! If base_obs vert type is not pressure; convert it to pressure
base_which = nint(query_location(base_obs_loc))
if (base_which == VERTISPRESSURE) then
    local_base_obs_loc = base_obs_loc
elseif (base_which == VERTISLEVEL) then
    local_base_obs_loc = pressure_from_level(base_obs_loc)
elseif (base_which == VERTISSURFACE) then
    local_base_obs_loc = pressure_from_surface(base_obs_loc)
end if

local_base_which   = VERTISPRESSURE
!! DEBUG: comment this in if you want to bypass the top level damping code below.
!call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
!                       num_close, close_ind, dist)
!return

! Get all the potentially close obs but no dist (optional argument dist(:) is not present)
call loc_get_close_obs(gc, local_base_obs_loc, base_obs_kind, obs_loc, &
                       obs_kind, num_close, close_ind)

! this should be a namelist parameter
threshold = HIGHEST_STATE_PRESSURE
if (threshold > 0.0_r8) thresh_wght = 1.0_r8/(threshold * threshold)

do k = 1, num_close
   t_ind = close_ind(k)
   obs_array = get_location(obs_loc(t_ind))
   obs_which = nint(query_location(obs_loc(t_ind)))

    if (obs_which == VERTISPRESSURE ) then
        local_obs_loc = obs_loc(t_ind)
    elseif (obs_which == VERTISLEVEL) then
        local_obs_loc = pressure_from_level(obs_loc(t_ind))
    elseif (obs_which == VERTISSURFACE) then
        local_obs_loc = pressure_from_surface(obs_loc(t_ind))
    else
        call error_handler(E_ERR, 'get_close_obs', 'Unsupported ob ' // &
                           'vertical coordinate type', source, revision, &
                           revdate)
    end if

!  nsc fri, 13mar09
!  allow a namelist specified kind string to restrict the impact of those
!  obs kinds to only other obs and state vars of the same kind.
   if (local_base_which == VERTISUNDEF) then
      ! The last argument, no_vert = .true., makes get_dist calculate horzontal distance only.
      dist(k) = get_dist(local_base_obs_loc, local_obs_loc, base_obs_kind, obs_kind(t_ind),.true.)
      ! Then no damping can be done since vertical distance is undefined.
      ! ? Is this routine called *both* to get model points close to a real obs,
      !   AND ob close to a model point?  I want damping in the latter case,
      !   even if ob has which_vert = VERTISUNDEF.
      !   I think that testing on local_base_which will do that.
   else
      dist(k) = get_dist(local_base_obs_loc, local_obs_loc, base_obs_kind, obs_kind(t_ind))

      ! Damp the influence of obs (below the namelist variable highest_obs_pressure_mb) 
      ! on variables above highest_state_pressure_mb.  
      ! This section could also change the distance based on the KIND_s of the base_obs and obs.
   
      ! dist = 0 for some for synthetic obs.
      ! Additive increase, based on height above threshold, works better than multiplicative
   
      ! See model_mod circa 1/1/2007 for other damping algorithms.
   
      increment = threshold - query_location(local_obs_loc, 'VLOC')
      ! This if-test handles the case where no damping is performed, i.e. 
      ! highest_state_pressure_mb = 0 and threshold = 0.
      if (increment > 0) then
         dist(k) = dist(k) + increment * increment * thresh_wght
      end if
   endif

end do

end subroutine get_close_obs

!--------------------------------------------------------------------
!-HELPER ROUTINES----------------------------------------------------
!--------------------------------------------------------------------

function pressure_from_level(lev_loc) result (p_loc)
    type(location_type), intent(in) :: lev_loc
    type(location_type)             :: p_loc

    real(kind=r8)       :: lev_loc_array(3)
    real(kind=r8)       :: p_surf
    real(kind=r8)       :: p_profile(numlevs)
    type(location_type) :: surface_loc

    integer :: interp_status
    
    lev_loc_array = get_location(lev_loc)
    surface_loc   = set_location(lev_loc_array(1), lev_loc_array(2), &
                                 0.0_r8, VERTISSURFACE)
    call model_interpolate(ensemble_mean, surface_loc, KIND_PRESSURE, &
                           p_surf, interp_status)
    call pres_sounding(numlevs, asig, bsig, p_surf, p_profile)

    p_loc = set_location(lev_loc_array(1), lev_loc_array(2), &
                         p_profile(nint(lev_loc_array(3))),&
                         VERTISPRESSURE)
end function pressure_from_level

function pressure_from_surface(surf_loc) result (p_loc)
    type(location_type), intent(in) :: surf_loc
    type(location_type)             :: p_loc

    real(kind=r8)       :: surf_loc_array(3)
    real(kind=r8)       :: p_surf
    real(kind=r8)       :: p_profile(numlevs)
    type(location_type) :: surface_loc

    integer :: interp_status
    
    surf_loc_array = get_location(surf_loc)
    surface_loc   = set_location(surf_loc_array(1), surf_loc_array(2), &
                                 0.0_r8, VERTISSURFACE)
    call model_interpolate(ensemble_mean, surface_loc, KIND_PRESSURE, &
                           p_surf, interp_status)
    p_loc = set_location(surf_loc_array(1), surf_loc_array(2), &
                         p_surf*HPA_TO_PA, VERTISPRESSURE)
end function pressure_from_surface

subroutine restart_file_to_sv(filename, state_vector, model_time)
!------------------------------------------------------------------
! Reads the current time and state variables from a nogaps restart
! file and packs them into a dart state vector.

character(len=*), intent(in)    :: filename 
real(r8),         intent(inout) :: state_vector(:)
type(time_type),  intent(out)   :: model_time

! temp space to hold data while we are reading it
real(r8) :: data_2d_array(numlons,numlats), data_3d_array(numlons,numlats,numlevs)
integer  :: i, j, k, ivar, indx

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: varname 
integer :: VarID, numdims, dimlen
integer :: ncid, iyear, imonth, iday, ihour, iminute, isecond, nc_rc
character(len=256) :: myerrorstring 

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

! Check that the input file exists ... 
! Read the time data. 
! Note from Nancy Norton as pertains time:
! "The time recorded in the nogaps2 restart files is the current time,
! which corresponds to the time of the XXXX_CUR variables.
!
! current time is determined from iyear, imonth, iday, and *seconds_this_day*
!
! The ihour, iminute, and isecond variables are used for internal
! model counting purposes, but because isecond is rounded to the nearest
! integer, it is possible that using ihour,iminute,isecond information
! on the restart file to determine the exact curtime would give you a 
! slightly wrong answer."
!
! DART only knows about integer number of seconds, so using the rounded one
! is what we would have to do anyway ... and we already have a set_date routine
! that takes ihour, iminute, isecond information.

if ( .not. file_exist(filename) ) then
   write(msgstring,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'restart_file_to_sv',msgstring,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'restart_file_to_sv', 'open '//trim(filename))
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iyear'  , iyear), &
                  'restart_file_to_sv', 'get_att iyear')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'imonth' , imonth), &
                  'restart_file_to_sv', 'get_att imonth')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iday'   , iday), &
                  'restart_file_to_sv', 'get_att iday')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'ihour'  , ihour), &
                  'restart_file_to_sv', 'get_att ihour')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iminute', iminute), &
                  'restart_file_to_sv', 'get_att iminute')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'isecond', isecond), &
                  'restart_file_to_sv', 'get_att isecond')

! FIXME: we don't allow a real year of 0 - add one for now, but
! THIS MUST BE FIXED IN ANOTHER WAY!
if (iyear == 0) then
  call error_handler(E_MSG, 'restart_file_to_sv', &
                     'WARNING!!!   year 0 not supported; setting to year 1')
  iyear = 1
endif

model_time = set_date(iyear, imonth, iday, ihour, iminute, isecond)

if (do_output()) &
    call print_time(model_time,'time for restart file '//trim(filename))
if (do_output()) &
    call print_date(model_time,'date for restart file '//trim(filename))

! Start counting and filling the state vector one item at a time,
! repacking the 3d arrays into a single 1d list of numbers.
! These must be a fixed number and in a fixed order.

indx = 1

! fill SALT, TEMP, UVEL, VVEL in that order
! The nogaps restart files have two time steps for each variable,
! the variables are named SALT_CUR and SALT_OLD ... for example.
! We are only interested in the CURrent time step.

do ivar=1, n3dfields

   varname = trim(progvarnames(ivar))//'_CUR'
   myerrorstring = trim(filename)//' '//trim(varname)

   ! Is the netCDF variable the right shape?

   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            'restart_file_to_sv', 'inq_varid '//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid,VarId,dimids=dimIDs,ndims=numdims), &
            'restart_file_to_sv', 'inquire '//trim(myerrorstring))

   if (numdims /= 3) then
      write(msgstring,*) trim(myerrorstring),' does not have exactly 3 dimensions'
      call error_handler(E_ERR,'restart_file_to_sv',msgstring,source,revision,revdate)
   endif

   do i = 1,numdims
      write(msgstring,'(''inquire dimension'',i2,A)') i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'restart_file_to_sv', msgstring)

      if (dimlen /= size(data_3d_array,i)) then
         write(msgstring,*) trim(myerrorstring),'dim/dimlen',i,dimlen,'not',size(data_3d_array,i)
         call error_handler(E_ERR,'restart_file_to_sv',msgstring,source,revision,revdate)
      endif
   enddo   

   ! Actually get the variable and stuff it into the array

   call nc_check(nf90_get_var(ncid, VarID, data_3d_array), 'restart_file_to_sv', &
                'get_var '//trim(varname))

   do k = 1, numlevs   ! size(data_3d_array,3)
   do j = 1, numlats   ! size(data_3d_array,2)
   do i = 1, numlons   ! size(data_3d_array,1)
      state_vector(indx) = data_3d_array(i, j, k)
      indx = indx + 1
   enddo
   enddo
   enddo

enddo

! and finally, PSURF (and any other 2d fields)
do ivar=(n3dfields+1), (n3dfields+n2dfields)

   varname = trim(progvarnames(ivar))//'_CUR'
   myerrorstring = trim(varname)//' '//trim(filename)

   ! Is the netCDF variable the right shape?

   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            'restart_file_to_sv', 'inq_varid '//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid,VarId,dimids=dimIDs,ndims=numdims), &
            'restart_file_to_sv', 'inquire '//trim(myerrorstring))

   if (numdims /= 2) then
      write(msgstring,*) trim(myerrorstring),' does not have exactly 2 dimensions'
      call error_handler(E_ERR,'restart_file_to_sv',msgstring,source,revision,revdate)
   endif

   do i = 1,numdims
      write(msgstring,'(''inquire dimension'',i2,A)') i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'restart_file_to_sv', msgstring)

      if (dimlen /= size(data_2d_array,i)) then
         write(msgstring,*) trim(myerrorstring),'dim/dimlen',i,dimlen,'not',size(data_2d_array,i)
         call error_handler(E_ERR,'restart_file_to_sv',msgstring,source,revision,revdate)
      endif
   enddo   

   ! Actually get the variable and stuff it into the array

   call nc_check(nf90_get_var(ncid, VarID, data_2d_array), 'restart_file_to_sv', &
                'get_var '//trim(varname))

   do j = 1, numlats   ! size(data_3d_array,2)
   do i = 1, numlons   ! size(data_3d_array,1)
      state_vector(indx) = data_2d_array(i, j)
      indx = indx + 1
   enddo
   enddo

enddo

end subroutine restart_file_to_sv



subroutine sv_to_restart_file(state_vector, statedate, ncid, fname)
!------------------------------------------------------------------
! Writes the current time and state variables from a dart state
! vector (1d fortran array) into a NOGAPS netcdf restart file.
!
real(r8),         intent(in) :: state_vector(:)
type(time_type),  intent(in) :: statedate
integer,          intent(in) :: ncid
character(len=*), intent(in) :: fname 

integer :: ierr

!----------------------------------------------------------------------
! Get the show underway
!----------------------------------------------------------------------

if ( .not. module_initialized ) call static_init_model

!----------------------------------------------------------------------
! reuse the nc_write_model_xxxx routines.
!----------------------------------------------------------------------

ierr = nc_write_model_atts(ncid)
if (ierr /= 0 ) then
   write(msgstring,*)'nc_write_model_atts failed '
   call error_handler(E_ERR,'sv_to_restart_file',msgstring,source,revision,revdate) 
endif

ierr = nc_write_model_vars(ncid, state_vector, 1, 1)
if (ierr /= 0 ) then
   write(msgstring,*)'nc_write_model_vars failed '
   call error_handler(E_ERR,'sv_to_restart_file',msgstring,source,revision,revdate) 
endif

call nc_check(nf90_close(ncid), 'sv_to_restart_file', 'close '//trim(fname))

end subroutine sv_to_restart_file



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

if (dim1 /= numlons) then
   write(msgstring,*)trim(varname),' 2d array dim 1 ',dim1,' /= ',numlons
   call error_handler(E_ERR,'vector_to_2d_prog_var',msgstring,source,revision,revdate) 
endif
if (dim2 /= numlats) then
   write(msgstring,*)trim(varname),' 2d array dim 2 ',dim2,' /= ',numlats
   call error_handler(E_ERR,'vector_to_2d_prog_var',msgstring,source,revision,revdate) 
endif

ii = start_index(varindex)

do j = 1,numlats   ! latitudes
do i = 1,numlons   ! longitudes
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

if (dim1 /= numlons) then
   write(msgstring,*)trim(varname),' 3d array dim 1 ',dim1,' /= ',numlons
   call error_handler(E_ERR,'vector_to_3d_prog_var',msgstring,source,revision,revdate) 
endif
if (dim2 /= numlats) then
   write(msgstring,*)trim(varname),' 3d array dim 2 ',dim2,' /= ',numlats
   call error_handler(E_ERR,'vector_to_3d_prog_var',msgstring,source,revision,revdate) 
endif
if (dim3 /= numlevs) then
   write(msgstring,*)trim(varname),' 3d array dim 3 ',dim3,' /= ',numlevs
   call error_handler(E_ERR,'vector_to_3d_prog_var',msgstring,source,revision,revdate) 
endif

ii = start_index(varindex)

do k = 1,numlevs   ! vertical
do j = 1,numlats   ! latitudes
do i = 1,numlons   ! longitudes
   data_3d_array(i,j,k) = x(ii)
   ii = ii + 1
enddo
enddo
enddo

end subroutine vector_to_3d_prog_var



subroutine get_gridsize(num_x, num_y, num_z)
 integer, intent(out) :: num_x, num_y, num_z
!------------------------------------------------------------------
! public utility routine.

if ( .not. module_initialized ) call static_init_model

 num_x = numlons
 num_y = numlats
 num_z = numlevs

end subroutine get_gridsize



  subroutine write_grid_netcdf()
!------------------------------------------------------------------
!
! Write the grid to a netcdf file for checking.

integer :: ncid, NlonDimID, NlatDimID, NzDimID, NzeDimID
integer :: latVarID, lonVarID, levVarID, sinlVarID, asigVarID, bsigVarID
integer :: phisVarID

integer :: dimids(2);

if ( .not. module_initialized ) call static_init_model

call nc_check(nf90_create('dart_grid.nc', NF90_CLOBBER, ncid),'write_grid_netcdf:create')

! define dimensions

call nc_check(nf90_def_dim(ncid, 'lon',      numlons  , NlonDimID),'write_grid_netcdf:defdim lon')
call nc_check(nf90_def_dim(ncid, 'lat',      numlats  , NlatDimID),'write_grid_netcdf:defdim lat')
call nc_check(nf90_def_dim(ncid, 'lev',      numlevs  ,   NzDimID),'write_grid_netcdf:defdim lev')
call nc_check(nf90_def_dim(ncid, 'nlevedges', numlevs+1, NzeDimID),'write_grid_netcdf:defdim edges')

dimids(1) = NlonDimID 
dimids(2) = NlatDimID 

! define variables

! we should add attributes to say what units the grids are in (degrees).
call nc_check(nf90_def_var(ncid, 'phis', nf90_double,    dimids, phisVarID),'write_grid_netcdf:defvar phis')
call nc_check(nf90_def_var(ncid,  'lon', nf90_double, NlonDimID,  lonVarID),'write_grid_netcdf:defvar lon')
call nc_check(nf90_def_var(ncid,  'lat', nf90_double, NlatDimID,  latVarID),'write_grid_netcdf:defvar lat')
call nc_check(nf90_def_var(ncid, 'sinl', nf90_double, NlatDimID, sinlVarID),'write_grid_netcdf:defvar sinl')
call nc_check(nf90_def_var(ncid,  'lev', nf90_double,   NzDimID,  levVarID),'write_grid_netcdf:defvar lev')
call nc_check(nf90_def_var(ncid, 'asig', nf90_double,  NzeDimID, asigVarID),'write_grid_netcdf:defvar asig')
call nc_check(nf90_def_var(ncid, 'bsig', nf90_double,  NzeDimID, bsigVarID),'write_grid_netcdf:defvar bsig')

call nc_check(nf90_put_att(ncid,phisVarID,'long_name','geopotential height'), 'write_grid_netcdf')
call nc_check(nf90_put_att(ncid,phisVarID,    'units','m^2/s^2'            ), 'write_grid_netcdf')

call nc_check(nf90_put_att(ncid,lonVarID,'long_name','grid longitudes'), 'write_grid_netcdf')
call nc_check(nf90_put_att(ncid,lonVarID,    'units','radians'        ), 'write_grid_netcdf')

call nc_check(nf90_put_att(ncid,latVarID,'long_name','grid latitudes'), 'write_grid_netcdf')
call nc_check(nf90_put_att(ncid,latVarID,    'units','radians'       ), 'write_grid_netcdf')

call nc_check(nf90_put_att(ncid,levVarID,'long_name','grid levels'), 'write_grid_netcdf')
call nc_check(nf90_put_att(ncid,levVarID,    'units','varies'     ), 'write_grid_netcdf')

call nc_check(nf90_enddef(ncid),'write_grid_netcdf:enddef')

! fill variables

call nc_check(nf90_put_var(ncid, phisVarID, phis),'write_grid_netcdf:put phis')
call nc_check(nf90_put_var(ncid,  lonVarID, lons),'write_grid_netcdf:put lons')
call nc_check(nf90_put_var(ncid,  latVarID, lats),'write_grid_netcdf:put lats')
call nc_check(nf90_put_var(ncid,  levVarID, levs),'write_grid_netcdf:put levs')
call nc_check(nf90_put_var(ncid, sinlVarID, sinl),'write_grid_netcdf:put sinl')
call nc_check(nf90_put_var(ncid, asigVarID, asig),'write_grid_netcdf:put asig')
call nc_check(nf90_put_var(ncid, bsigVarID, bsig),'write_grid_netcdf:put bsig')

call nc_check(nf90_close(ncid),'write_grid_netcdf:close')

end subroutine write_grid_netcdf



! Helper functions for interpolation

subroutine get_column_indices(target_i, target_j, i_size, j_size, &
                              column_size, column_indices)
    integer, intent(in)  :: target_i
    integer, intent(in)  :: target_j
    integer, intent(in)  :: i_size
    integer, intent(in)  :: j_size
    integer, intent(in)  :: column_size
    integer, intent(out) :: column_indices(column_size)

    integer :: cur_level
    integer :: index_1d
    integer :: size_2d

    index_1d = index_2d_to_1d(target_i, target_j, i_size)
    size_2d  = i_size*j_size
    do cur_level = 1, column_size
        column_indices(cur_level) = index_1d + size_2d*(cur_level-1)
    end do 

end subroutine get_column_indices



function index_2d_to_1d(i, j, i_size) result (index_1d)
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(in) :: i_size
    integer             :: index_1d

    index_1d = (i_size*(j-1) + i)
end function index_2d_to_1d



subroutine pres_sounding(lm, asig, bsig, ps, pv)
!
!  This is lifted from prexp.f in NOGAPS ... it appears that
!  the "half" levels are traditional hybrid-sigma coordinates
!  defined by p(k) = A(k+1) + B(k+1)*ps
!  The "full" levels are determined by the formulas below...
!
! Input:
!    lm = number of "full" vertical levels
!    asig(lm+1),bsig(lm+1) = hybrid coordinate coefficients
!    ps = surface pressure (in hPa)
!
! Ouput:
!    pv(lm) = pressure on "full" levels (in Pa)
!
!
               
      integer,  intent(in)  :: lm
      real(r8), intent(in)  :: asig(lm+1)
      real(r8), intent(in)  :: bsig(lm+1)
      real(r8), intent(in)  :: ps
      real(r8), intent(out) :: pv(lm)

! Local variables:
      integer  :: k
      real(r8) :: capa
      real(r8) :: capap1
      real(r8) :: opok
      real(r8) :: p0
      real(r8) :: pl2(lm)
      real(r8) :: pk2(lm)
      real(r8) :: pk(lm)
      real(r8) :: ptop
      real(r8) :: ptopk
!
!
      capa   = 1.0_r8/3.5_r8 
      capap1 = 1.0_r8 + capa
      p0     = 1000.0_r8
      opok   = 1.0_r8/p0**capa
      ptop   = asig(1)
      ptopk  = ptop*opok*ptop**capa
!  pressure on the top "full" level
!
      k = 1
      pl2(k) = asig(k+1) + bsig(k+1)*ps
      pk2(k) = opok*pl2(k)**capa
      pk(k)  = ( pl2(k)*pk2(k)-ptopk )/( capap1*( pl2(k)-ptop ) )
      pv(k)  = p0*pk(k)*pk(k)*pk(k)*sqrt( pk(k) )                               
!
!  pressure on the rest of the levels
!
      do k = 2, lm
         pl2(k) = asig(k+1) + bsig(k+1)*ps
         pk2(k) = opok*pl2(k)**capa
         pk(k)  = ( pl2(k)*pk2(k)-pl2(k-1)*pk2(k-1))   &
                  / ( capap1*( pl2(k)-pl2(k-1) ) )
         pv(k)  = p0*pk(k)*pk(k)*pk(k)*sqrt( pk(k) )
      end do

      pv(:) = pv(:) * HPA_TO_PA

      return
end subroutine pres_sounding



function get_start_index_by_kind(var_kind) result (offset)
    integer, intent(in) :: var_kind
    integer             :: offset

    integer :: var_index 
    select case (var_kind)
        case(KIND_U_WIND_COMPONENT)
            var_index = U_index
        case(KIND_V_WIND_COMPONENT)
            var_index = V_index
        case(KIND_TEMPERATURE)
            var_index = T_index
        case(KIND_SPECIFIC_HUMIDITY)
            var_index = Q_index
        case(KIND_PRESSURE)
            var_index = PSURF_index
        case default
            call error_handler(E_ERR, 'get_start_index_by_kind', &
                               'Unknown variable kind!', source, &
                               revision, revdate)
    end select

    offset = start_index(var_index)

end function get_start_index_by_kind


!===================================================================
! Test routines
!===================================================================

subroutine test_interpolate(x, test_pressure, start_lon)
    real(kind=r8), intent(in) :: x(:)
    real(kind=r8), intent(in) :: test_pressure
    real(kind=r8), intent(in) :: start_lon


    integer,       parameter :: TEST_POINTS = 100
    real(kind=r8), parameter :: START_LAT   = -80.0_r8
    real(kind=r8), parameter :: DELTA_LON   = 360.0_r8 / TEST_POINTS
    real(kind=r8), parameter :: DELTA_LAT   = 160.0_r8 / TEST_POINTS

    real(kind=r8)       :: test_lats(TEST_POINTS)
    real(kind=r8)       :: test_lons(TEST_POINTS)
    type(location_type) :: test_loc

    real(kind=r8)       :: interp_vals(TEST_POINTS, TEST_POINTS)

    integer :: i, j
    integer :: interp_status

    do i = 1, TEST_POINTS
        test_lats(i) = START_LAT + (i-1)*DELTA_LAT
        test_lons(i) = start_lon + (i-1)*DELTA_LON
        if (test_lons(i) > 360.0_r8) test_lons(i) = test_lons(i) - 360.0_r8
    end do

    do j = 1, TEST_POINTS
        do i = 1, TEST_POINTS
            test_loc = set_location(test_lons(i), test_lats(j), &
                                    test_pressure, VERTISSURFACE)
            call model_interpolate(x, test_loc, KIND_PRESSURE, &
                                   interp_vals(i,j), interp_status)
            write (42,*) interp_vals(i,j)
            test_loc = set_location(test_lons(i), test_lats(j), &
                                    test_pressure, VERTISPRESSURE)
            call model_interpolate(x, test_loc, KIND_TEMPERATURE, &
                                   interp_vals(i,j), interp_status)
            write (43,*) interp_vals(i,j)
            test_loc = set_location(test_lons(i), test_lats(j), &
                                    test_pressure, VERTISPRESSURE)
            call model_interpolate(x, test_loc, KIND_U_WIND_COMPONENT, &
                                   interp_vals(i,j), interp_status)
            write (44,*) interp_vals(i,j)
            test_loc = set_location(test_lons(i), test_lats(j), &
                                    test_pressure, VERTISPRESSURE)
            call model_interpolate(x, test_loc, KIND_V_WIND_COMPONENT, &
                                   interp_vals(i,j), interp_status)
            write (45,*) interp_vals(i,j)
            test_loc = set_location(test_lons(i), test_lats(j), &
                                    test_pressure, VERTISPRESSURE)
            call model_interpolate(x, test_loc, KIND_SPECIFIC_HUMIDITY, &
                                   interp_vals(i,j), interp_status)
            write (46,*) interp_vals(i,j)
        end do
    end do

end subroutine test_interpolate


!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

