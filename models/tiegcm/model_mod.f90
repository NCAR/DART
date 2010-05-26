! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!------------------------------------------------------------------
!
! Interface for HAO-TIEGCM 
!
!------------------------------------------------------------------

! DART Modules
use        types_mod, only : r8, digits12
use time_manager_mod, only : time_type, set_calendar_type, set_time_missing, &
                             set_time, get_time, print_time, &
                             set_date, get_date, print_date, & 
                             operator(*),  operator(+), operator(-),    &
                             operator(>),  operator(<), operator(/),    &
                             operator(/=), operator(<=)
use     location_mod, only : location_type, get_close_maxdist_init,                 &
                             get_close_obs_init, get_close_obs,                     &
                             set_location, get_location, query_location,            &
                             get_dist, vert_is_height
use    utilities_mod, only : file_exist, open_file, close_file,                     &       
                             error_handler, E_ERR, E_MSG, E_WARN, nmlfileunit,      & 
                             do_output, find_namelist_in_file, check_namelist_read, &
                             do_nml_file, do_nml_term, nc_check,                    &
                             register_module
use     obs_kind_mod, only : KIND_U_WIND_COMPONENT, &! just for definition
                             KIND_V_WIND_COMPONENT, &! just for definition
                             KIND_TEMPERATURE,      &! just for definition
                             KIND_ELECTRON_DENSITY, &! for Ne obs 
                             KIND_DENSITY            ! for neutral density observation
                                                     ! (for now [O]; minor species ignored)
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
  
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
!TIEGCM specific routines
public :: model_type,             &
          init_model_instance,    &
          end_model_instance,     &
          prog_var_to_vector,     &
          vector_to_prog_var,     &
          read_TIEGCM_restart,    &
          update_TIEGCM_restart,  &
          read_TIEGCM_definition, &
          read_TIEGCM_secondary

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!------------------------------------------------------------------
! define model parameters


integer                               :: nilev, nlev, nlon, nlat
real(r8),dimension(:),    allocatable :: lons, lats, levs, ilevs
real(r8),dimension(:,:,:),allocatable :: ZGtiegcm !! auxiliary variable(geopotential height[cm])
real                                  :: TIEGCM_missing_value !! global attribute

integer                               :: time_step_seconds = 120 ! tiegcm time step
integer                               :: time_step_days = 0      ! tiegcm time step
character (len=19)                    :: restart_file = 'tiegcm_restart_p.nc'
character (len=11)                    :: secondary_file = 'tiegcm_s.nc'

integer, parameter                    :: TYPE_local_NE    = 0  
integer, parameter                    :: TYPE_local_TN    = 1  
integer, parameter                    :: TYPE_local_TN_NM = 2  
integer, parameter                    :: TYPE_local_UN    = 3  
integer, parameter                    :: TYPE_local_UN_NM = 4  
integer, parameter                    :: TYPE_local_VN    = 5  
integer, parameter                    :: TYPE_local_VN_NM = 6  
integer, parameter                    :: TYPE_local_O1    = 7
integer, parameter                    :: TYPE_local_O1_NM = 8

type(time_type)                       :: time_step
integer                               :: model_size
                                        
type model_type
  real(r8), pointer                   :: vars_3d(:,:,:,:)
  type(time_type)                     :: valid_time
end type model_type

integer                               :: state_num_3d = 9  
                                          ! -- interface levels --
                                          ! NE 
                                          ! -- midpoint levels --
                                          ! O1 O1_NM (O2 O2_NM NO NO_NM excluded for now)    
                                          ! -- midpoint levels; top slot missing --
                                          ! TN TN_NM UN UN_NM VN VN_NM
logical                               :: output_state_vector = .false.
                                          ! .true.  results in a "state-vector" netCDF file
                                          ! .false. results in a "prognostic-var" netCDF file
namelist /model_nml/ output_state_vector, state_num_3d

logical                               :: first_pert_call = .true.
type(random_seq_type)                 :: random_seq



!------------------------------------------------------------------

character(len = 129) :: msgstring
logical, save :: module_initialized = .false.

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

 integer  :: i
 integer  :: iunit, io

if (module_initialized) return ! only need to do this once
 
! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

!! Read the namelist entry for model_mod from input.nml
!call find_namelist_in_file("input.nml", "model_nml", iunit)
!read(iunit, nml = model_nml, iostat = io)
!call check_namelist_read(iunit, io, "model_nml")

!if (do_nml_file()) write(nmlfileunit, nml=model_nml)
!if (do_nml_term()) write(     *     , nml=model_nml)

! Reading in TIEGCM grid definition etc from TIEGCM restart file
call read_TIEGCM_definition(restart_file)

! Reading in Geopotential Height from TIEGCM secondary output file
call read_TIEGCM_secondary(secondary_file)

! Compute overall model size 
model_size = nlon * nlat * nlev * state_num_3d

if (do_output()) write(*,*)'nlon = ',nlon
if (do_output()) write(*,*)'nlat = ',nlat
if (do_output()) write(*,*)'nlev = ',nlev
if (do_output()) write(*,*)'n3D  = ',state_num_3d

! Might as well use the Gregorian Calendar
call set_calendar_type('Gregorian')

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

if ( .not. module_initialized ) call static_init_model

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

integer  :: i, vstatus, which_vert
integer  :: lat_below, lat_above, lon_below, lon_above
real(r8) :: lon_fract, temp_lon, lat_fract
real(r8) :: lon, lat, height, lon_lat_lev(3)
real(r8) :: bot_lon, top_lon, delta_lon, bot_lat, top_lat, delta_lat
real(r8) :: val(2,2), a(2)

if ( .not. module_initialized ) call static_init_model

! Default for successful return
istatus = 0
vstatus = 0

! Get the position, determine if it is model level or pressure in vertical
lon_lat_lev = get_location(location)
lon = lon_lat_lev(1) ! degree
lat = lon_lat_lev(2) ! degree
if(vert_is_height(location)) then
   height = lon_lat_lev(3)
else
   which_vert = nint(query_location(location))
   write(msgstring,*) 'vertical coordinate type:',which_vert,' cannot be handled'
   call error_handler(E_ERR,'model_interpolate',msgstring,source,revision,revdate)
endif

! Get lon and lat grid specs
bot_lon   = lons(1)                ! 0
top_lon   = lons(nlon)             ! 355 
delta_lon = abs((lons(1)-lons(2))) ! 5
bot_lat   = lats(1)                ! -87.5
top_lat   = lats(nlat)             ! 87.5
delta_lat = abs((lats(1)-lats(2))) ! 5


! Compute bracketing lon indices: DART [0, 360] TIEGCM [0, 355]
if(lon >= bot_lon .and. lon <= top_lon) then !  0 <= lon <= 355 
   lon_below = int((lon - bot_lon) / delta_lon) + 1
   lon_above = lon_below + 1
   lon_fract = (lon - lons(lon_below)) / delta_lon
elseif (lon < bot_lon) then                  !  lon < 0
   temp_lon = lon + 360.0
   lon_below = int((temp_lon - bot_lon) / delta_lon) + 1
   lon_above = lon_below + 1
   lon_fract = (temp_lon - lons(lon_below)) / delta_lon  

elseif (lon >= (top_lon+delta_lon)) then     !  360 <= lon  
   temp_lon = lon - 360.0
   lon_below = int((temp_lon - bot_lon) / delta_lon) + 1
   lon_above = lon_below + 1
   lon_fract = (temp_lon - lons(lon_below)) / delta_lon  
else                                         !  355 < lon < 360 at wraparound point
   lon_below = nlon
   lon_above = 1
   lon_fract = (lon - top_lon) / delta_lon
endif

! compute neighboring lat rows: TIEGCM [-87.5, 87.5] DART [-90, 90]
! NEED TO BE VERY CAREFUL ABOUT POLES; WHAT'S BEING DONE IS NOT GREAT!
if(lat >= bot_lat .and. lat <= top_lat) then ! -87.5 <= lat <= 87.5
   lat_below = int((lat - bot_lat) / delta_lat) + 1
   lat_above = lat_above + 1
   lat_fract = (lat - lats(lat_below) ) / delta_lat
else if(lat < bot_lat) then ! South of bottom lat 
   lat_below = 1
   lat_above = 1
   lat_fract = 1.
else                        ! North of top lat
   lat_below = nlat
   lat_above = nlat
   lat_fract = 1.
endif

! Now, need to find the values for the four corners
                     call get_val(val(1, 1), x, lon_below, lat_below, height, itype, vstatus)
   if (vstatus /= 1) call get_val(val(1, 2), x, lon_below, lat_above, height, itype, vstatus)
   if (vstatus /= 1) call get_val(val(2, 1), x, lon_above, lat_below, height, itype, vstatus)
   if (vstatus /= 1) call get_val(val(2, 2), x, lon_above, lat_above, height, itype, vstatus)


! istatus   meaning                  return expected obs?   assimilate?
! 0         obs and model are fine;  yes                    yes
! 1         fatal problem;           no                     no
! 2         exclude valid obs        yes                    no
!
istatus = vstatus
if(istatus /= 1) then
   do i = 1, 2
      a(i) = lon_fract * val(2, i) + (1.0 - lon_fract) * val(1, i)
   end do
   obs_val = lat_fract * a(2) + (1.0 - lat_fract) * a(1)
else
   obs_val = 0.
endif

!print*, 'model_interpolate', lon, lat, height,obs_val


end subroutine model_interpolate



subroutine get_val(val, x, lon_index, lat_index, height, obs_kind, istatus)
!------------------------------------------------------------------
!
real(r8), intent(out) :: val
real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: lon_index, lat_index
real(r8), intent(in)  :: height
integer,  intent(in)  :: obs_kind
integer,  intent(out) :: istatus

integer               :: var_type
integer               :: k, lev_top, lev_bottom
real(r8)              :: zgrid, delta_z
real(r8)              :: val_top, val_bottom, frac_lev

! No errors to start with
istatus = 0

! To find a layer height: what's the unit of height [m] 
h_loop:do k = 1, nlev
  zgrid = ZGtiegcm(lon_index,lat_index,k)/100.0_r8 ! [m] = ZGtiegcm/100 [cm]
  if (height <= zgrid) then  
    lev_top = k
    lev_bottom = lev_top -1
    delta_z = zgrid - ZGtiegcm(lon_index,lat_index,lev_bottom)
    frac_lev = (zgrid - height)/delta_z
    exit h_loop
  endif
enddo h_loop

if (obs_kind == KIND_DENSITY) then

  var_type   = TYPE_local_O1
  val_top    = x(get_index(lat_index, lon_index, lev_top, var_type))
  val_bottom = x(get_index(lat_index, lon_index, lev_bottom, var_type))

elseif (obs_kind == KIND_ELECTRON_DENSITY) then

  var_type   = TYPE_local_NE
  val_top    = x(get_index(lat_index, lon_index, lev_top, var_type))
  val_bottom = x(get_index(lat_index, lon_index, lev_bottom, var_type))

else
 
  istatus = 1
  val = 0.
  return

end if

val = frac_lev * val_bottom  +  (1.0 - frac_lev) * val_top

end subroutine get_val



function get_index(lat_index, lon_index, lev_index, var_type)
!------------------------------------------------------------------
!
integer,  intent(in) :: lat_index, lon_index, lev_index, var_type
integer              :: get_index

get_index = 1 + var_type + (lev_index -1)*state_num_3d           &
                         + (lat_index -1)*state_num_3d*nlev      &
                         + (lon_index -1)*state_num_3d*nlev*nlat

end function get_index



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

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

integer  :: indx, num_per_col, col_num, col_elem
integer  :: lon_index, lat_index, lev_index
real(r8) :: lon, lat, lev
integer  :: local_var_type, var_type_temp

if ( .not. module_initialized ) call static_init_model

! Easier to compute with a 0 to size -1 index
indx = index_in -1

! Compute number of items per column
num_per_col = nlev * state_num_3d

! What column is this index in
col_num  = indx / num_per_col
col_elem = indx - col_num * num_per_col 

! what lon and lat index for this column
lon_index = col_num /nlat
lat_index = col_num - lon_index * nlat

! Now figure out which beast in column this is
lev_index = col_elem / state_num_3d

! Get actual lon lat values from static_init_model arrays
lon = lons(lon_index + 1)
lat = lats(lat_index + 1)

! Find which var_type this element is
var_type_temp = mod(col_elem, state_num_3d)

if      (var_type_temp == 0) then  !NE
  local_var_type = KIND_ELECTRON_DENSITY
else if (var_type_temp == 1) then  !TN
  local_var_type = KIND_TEMPERATURE
else if (var_type_temp == 2) then  !TN_NM
  local_var_type = KIND_TEMPERATURE
else if (var_type_temp == 3) then  !UN
  local_var_type = KIND_U_WIND_COMPONENT
else if (var_type_temp == 4) then  !UN_NM
  local_var_type = KIND_U_WIND_COMPONENT
else if (var_type_temp == 5) then  !VN
  local_var_type = KIND_V_WIND_COMPONENT
else if (var_type_temp == 6) then  !VN_NM
  local_var_type = KIND_V_WIND_COMPONENT
else if (var_type_temp == 7) then  !O1
  local_var_type = KIND_DENSITY
else if (var_type_temp == 8) then  !O1_NM
  local_var_type = KIND_DENSITY
else
   write(msgstring,*)"unknown var_type for index ",index_in
   call error_handler(E_ERR,"get_state_meta_data", msgstring, source, revision, revdate)
endif

if (var_type_temp == 0) then                !NE defined at interface levels
  lev = ilevs(lev_index + 1)
else
  lev = levs(lev_index + 1)        !TN UN VN O1 defined at midpoints
endif

location = set_location(lon,lat,lev,2) ! 2 == pressure  3 == height

! If the type is wanted, return it
if(present(var_type)) var_type = local_var_type

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

if ( .not. module_initialized ) call static_init_model

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


integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

integer :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)

integer :: StateVarVarID   ! netCDF pointer to state variable coordinate array
integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

integer :: TNVarID, TN_NMVarID, UNVarID, UN_NMVarID, VNVarID, VN_NMVarID
integer :: O1VarID, O1_NMVarID 
integer :: NEVarID
integer :: lonDimID, latDimID, levDimID, ilevDimID
integer :: lonVarID, latVarID, levVarID, ilevVarID

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer :: i

if ( .not. module_initialized ) call static_init_model

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! and then put into define mode.
!-------------------------------------------------------------------------------

ierr = -1 ! assume things go poorly

call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, &
              nAttributes, unlimitedDimID), 'nc_write_model_atts','inquire')
call nc_check(nf90_Redef(ncFileID),'nc_write_model_atts','redef')

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension. 
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID),&
       'nc_write_model_atts', 'copy dimid')
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID),&
       'nc_write_model_atts', 'time dimid')

if ( TimeDimID /= unlimitedDimId ) then
   write(msgstring,*)"Time Dimension ID ",TimeDimID, &
                     " should equal Unlimited Dimension ID",unlimitedDimID
   call error_handler(E_ERR,"nc_write_model_atts", msgstring, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------
call nc_check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                        len=model_size, dimid = StateVarDimID),&
       'nc_write_model_atts', 'state def_dim')

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date" ,str1    ),&
       'nc_write_model_atts', 'creation put')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source"  ,source  ),&
       'nc_write_model_atts', 'source put')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision),&
       'nc_write_model_atts', 'revision put')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate" ,revdate ),&
       'nc_write_model_atts', 'revdate put')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","TIEGCM"         ),&
       'nc_write_model_atts', 'model put')

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
              dimids=StateVarDimID, varid=StateVarVarID), &
              'nc_write_model_atts', 'statevariable def_var')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"), &
              'nc_write_model_atts', 'statevariable long_name')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical"), &
              'nc_write_model_atts', 'statevariable units')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, model_size /)), &
              'nc_write_model_atts', 'statevariable valid_range')

   ! Define the actual (3D) state vector, which gets filled as time goes on ... 
   call nc_check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_real, &
              dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), &
              varid=StateVarID), 'nc_write_model_atts', 'state def_var')
   call nc_check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"), &
              'nc_write_model_atts', 'state long_name')

   ! Leave define mode so we can fill the coordinate variable.
   call nc_check(nf90_enddef(ncfileID), 'nc_write_model_atts', 'state enddef')

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ), &
              'nc_write_model_atts', 'state put_var')

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------

   ! This block is a stub for something more complicated.
   ! Usually, the control for the execution of this block is a namelist variable.
   ! Take a peek at the bgrid model_mod.f90 for a (rather complicated) example.

   !----------------------------------------------------------------------------
   ! Define the dimensions IDs
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_dim(ncid=ncFileID, name="lon", &
             & len = nlon,    dimid =   lonDimID), 'nc_write_model_atts')
   call nc_check(nf90_def_dim(ncid=ncFileID, name="lat", &
             & len = nlat,    dimid =   latDimID), 'nc_write_model_atts')
   call nc_check(nf90_def_dim(ncid=ncFileID, name="lev", &
             & len = nlev, dimid =   levDimID),  'nc_write_model_atts')
   call nc_check(nf90_def_dim(ncid=ncFileID, name="ilev", &
             & len = nilev, dimid =  ilevDimID), 'nc_write_model_atts')
   
   !----------------------------------------------------------------------------
   ! Create the (empty) Variables and the Attributes
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncFileID, name="lon", &
             & xtype=nf90_double, dimids=lonDimID, varid=lonVarID),&
             'nc_write_model_atts') 
   call nc_check(nf90_put_att(ncFileID, lonVarID, &
             & "long_name", "geographic longitude (-west, +east)"),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, lonVarID, "units", "degrees_east"),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, lonVarID, "valid_range", &
             & (/ -180.0_r8, 180.0_r8 /)),'nc_write_model_atts')

   call nc_check(nf90_def_var(ncFileID, name="lat", &
             & xtype=nf90_double, dimids=latDimID, varid=latVarID),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, latVarID, &
             & "long_name", "geographic latitude (-south +north)"),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, latVarID, "units", "degrees_north"),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, latVarID, "valid_range", &
             & (/ -90.0_r8, 90.0_r8 /)),'nc_write_model_atts')

   call nc_check(nf90_def_var(ncFileID, name="lev", &
             & xtype=nf90_double, dimids=levDimID, varid=levVarID),&
             'nc_write_model_atts') 
   call nc_check(nf90_put_att(ncFileID, levVarID, "long_name", "midpoint levels"),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, levVarID, "units", "Pa"),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, levVarID, "positive", "up"),&
             'nc_write_model_atts') 

   call nc_check(nf90_def_var(ncFileID, name="ilev", &
             & xtype=nf90_double, dimids=ilevDimID, varid=ilevVarID),&
             'nc_write_model_atts') 
   call nc_check(nf90_put_att(ncFileID, ilevVarID, "long_name", "interface levels"),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, ilevVarID, "units", "Pa"),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, ilevVarID, "positive", "up"),&
             'nc_write_model_atts') 

   !----------------------------------------------------------------------------
   ! Create attributes for the state vector
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncid=ncFileID, name="NE", xtype=nf90_real, &
       dimids = (/ lonDimID, latDimID, ilevDimID, MemberDimID, unlimitedDimID /), &
       varid  = NEVarID), 'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, NEVarID, "long_name", "electron density"), &
          'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, NEVarID, "units", "cm-3"), &
          'nc_write_model_atts')

   call nc_check(nf90_def_var(ncid=ncFileID, name="TN", xtype=nf90_real, &
       dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
       varid  = TNVarID), 'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, TNVarID, "long_name", "neutral temperature"), &
          'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, TNVarID, "units", "K"), &
          'nc_write_model_atts')

   call nc_check(nf90_def_var(ncid=ncFileID, name="TN_NM", xtype=nf90_real, &
       dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
       varid  = TN_NMVarID), 'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, TN_NMVarID, "long_name", &
       "neutral temperature (time N-1)"), 'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, TN_NMVarID, "units", "K"), &
          'nc_write_model_atts')
 
   call nc_check(nf90_def_var(ncid=ncFileID, name="UN", xtype=nf90_real, &
       dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
       varid  = UNVarID), 'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, UNVarID, "long_name", & 
        "neutral zonal wind (+east)"), 'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, UNVarID, "units", "cm/s"), & 
          'nc_write_model_atts')

   call nc_check(nf90_def_var(ncid=ncFileID, name="UN_NM", xtype=nf90_real, &
       dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
       varid  = UN_NMVarID), 'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, UN_NMVarID, "long_name", &
       "neutral zonal wind (+east) (time N-1)"), 'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, UN_NMVarID, "units", "cm/s"), & 
          'nc_write_model_atts')  

   call nc_check(nf90_def_var(ncid=ncFileID, name="VN", xtype=nf90_real, &
       dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
       varid  = VNVarID), 'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, VNVarID, "long_name", &
       "neutral meridional wind (+north)"), 'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, VNVarID, "units", "cm/s"), & 
          'nc_write_model_atts')

   call nc_check(nf90_def_var(ncid=ncFileID, name="VN_NM", xtype=nf90_real, &
       dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
       varid  = VN_NMVarID),'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, VN_NMVarID, "long_name", &
       "neutral meridional wind (+north) (time N-1)"), 'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, VN_NMVarID, "units", "cm/s"), & 
          'nc_write_model_atts')  

   call nc_check(nf90_def_var(ncid=ncFileID, name="O1", xtype=nf90_real, &
       dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
       varid  = O1VarID), 'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, O1VarID, "long_name", "atomic oxygen"), &
          'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, O1VarID, "units", "mmr"), & 
          'nc_write_model_atts')

   call nc_check(nf90_def_var(ncid=ncFileID, name='O1_NM', xtype=nf90_real, &
       dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
       varid  = O1_NMVarID), 'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, O1_NMVarID, "long_name", "atomic oxygen (time N-1)"), &
          'nc_write_model_atts')
   call nc_check(nf90_put_att(ncFileID, O1_NMVarID, "units", "mmr"), &
          'nc_write_model_atts')
  

   call nc_check(nf90_enddef(ncfileID), 'nc_write_model_atts', 'prognostic enddef')

   !-------------------------------------------------------------------------------
   ! Fill the variables
   !-------------------------------------------------------------------------------

   call nc_check(nf90_put_var(ncFileID, lonVarID, lons), & 
          'nc_write_model_atts', 'put_var lons')
   call nc_check(nf90_put_var(ncFileID, latVarID, lats), & 
          'nc_write_model_atts', 'put_var lats')
   call nc_check(nf90_put_var(ncFileID, levVarID, levs), & 
          'nc_write_model_atts', 'put_var levs')
   call nc_check(nf90_put_var(ncFileID,ilevVarID,ilevs), & 
          'nc_write_model_atts', 'put_var ilevs')

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'sync')
if (do_output()) write (*,*) 'nc_write_model_atts: netCDF file ', ncFileID, ' is synched '

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

integer         :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer         :: StateVarID
integer         :: TNVarID, TN_NMVarID, UNVarID, UN_NMVarID, VNVarID, VN_NMVarID  
integer         :: O1VarID, O1_NMVarID
integer         :: NEVarID

type(model_type):: var 

if ( .not. module_initialized ) call static_init_model
  
!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
!-------------------------------------------------------------------------------

ierr = -1 ! assume things go poorly

call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, &
                  nAttributes, unlimitedDimID), 'nc_write_model_vars', 'inquire')

if ( output_state_vector ) then

   call nc_check(NF90_inq_varid(ncFileID, 'state', StateVarID), &
          'nc_write_model_vars', 'state inq_varid' )
   call nc_check(NF90_put_var(ncFileID, StateVarID, statevec,  &
                start=(/ 1, copyindex, timeindex /)), &
          'nc_write_model_vars', 'state put_var')                   

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

   call init_model_instance(var)
   call vector_to_prog_var(statevec, var)

   ! the 'start' array is crucial. In the following example, 'ps' is a 2D
   ! array, and the netCDF variable "ps" is a 4D array [lat,lon,copy,time]
   !
   ! call check(NF90_inq_varid(ncFileID, "ps", psVarID), "ps inq_varid")
   ! call check(nf90_put_var( ncFileID, psVarID, global_Var%ps, &
   !                          start=(/ 1, 1, copyindex, timeindex /) ), "ps put_var")

   call nc_check(NF90_inq_varid(ncFileID, 'NE',    NEVarID), &
          'nc_write_model_vars', 'NE inq_varid') 
   call nc_check(nf90_put_var( ncFileID, NEVarID, var%vars_3d(:,:,:,1), & 
                 start=(/ 1, 1, 1, copyindex, timeindex /) ), &
          'nc_write_model_vars', 'NE put_var') 
 
   call nc_check(NF90_inq_varid(ncFileID, 'TN',    TNVarID), &
          'nc_write_model_vars', 'TN inq_varid')
   call nc_check(nf90_put_var( ncFileID, TNVarID, var%vars_3d(:,:,:,2), & 
                 start=(/ 1, 1, 1, copyindex, timeindex /) ), &
          'nc_write_model_vars', 'TN put_var')     

   call nc_check(NF90_inq_varid(ncFileID, 'TN_NM', TN_NMVarID), &
          'nc_write_model_vars', 'TN_NM inq_varid')
   call nc_check(nf90_put_var( ncFileID, TN_NMVarID, var%vars_3d(:,:,:,3),& 
                start=(/ 1, 1, 1, copyindex, timeindex /) ), &
          'nc_write_model_vars', 'TN_NM put_var')     

   call nc_check(NF90_inq_varid(ncFileID, 'UN',    UNVarID), &
          'nc_write_model_vars',  'UN inq_varid')
   call nc_check(nf90_put_var( ncFileID, UNVarID, var%vars_3d(:,:,:,4),& 
                start=(/ 1, 1, 1, copyindex, timeindex /) ), & 
          'nc_write_model_vars',  'UN put_var')     

   call nc_check(NF90_inq_varid(ncFileID, 'UN_NM', UN_NMVarID), &
          'nc_write_model_vars', 'UN_NM inq_varid') 
   call nc_check(nf90_put_var( ncFileID, UN_NMVarID, var%vars_3d(:,:,:,5),& 
                start=(/ 1, 1, 1, copyindex, timeindex /) ), &
          'nc_write_model_vars', 'UN_NM put_var')     

   call nc_check(NF90_inq_varid(ncFileID, 'VN',    VNVarID), &
          'nc_write_model_vars', 'VN inq_varid')
   call nc_check(nf90_put_var( ncFileID, VNVarID, var%vars_3d(:,:,:,6),& 
                start=(/ 1, 1, 1, copyindex, timeindex /) ), &
          'nc_write_model_vars', 'VN put_var')     

   call nc_check(NF90_inq_varid(ncFileID, 'VN_NM', VN_NMVarID), &
          'nc_write_model_vars', 'VN_NM inq_varid') 
   call nc_check(nf90_put_var( ncFileID, VN_NMVarID, var%vars_3d(:,:,:,7),& 
                start=(/ 1, 1, 1, copyindex, timeindex /) ), &
          'nc_write_model_vars', 'VN_NM put_var')     

   call nc_check(NF90_inq_varid(ncFileID, 'O1',    O1VarID), &
          'nc_write_model_vars', 'O1 inq_varid') 
   call nc_check(nf90_put_var( ncFileID, O1VarID, var%vars_3d(:,:,:,8),& 
                start=(/ 1, 1, 1, copyindex, timeindex /) ), &
          'nc_write_model_vars', 'O1 put_var')    

   call nc_check(NF90_inq_varid(ncFileID, 'O1_NM', O1_NMVarID), &
          'nc_write_model_vars', 'O1_NM inq_varid') 
   call nc_check(nf90_put_var( ncFileID, O1_NMVarID, var%vars_3d(:,:,:,9),& 
                start=(/ 1, 1, 1, copyindex, timeindex /) ), &
          'nc_write_model_vars', 'O1_NM put_var')

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

if (do_output()) write (*,*) 'nc_write_model_vars: Finished filling variables '
call nc_check(nf90_sync(ncFileID), 'nc_write_model_vars', 'sync')
if (do_output()) write (*,*) 'nc_write_model_vars: netCDF file is synched '

ierr = 0 ! If we got here, things went well.

call end_model_instance(Var)   ! should avoid any memory leaking

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

real(r8), intent(in)    :: state(:)
real(r8), intent(out)   :: pert_state(:)
logical,  intent(out)   :: interf_provided

integer                 :: i, variable_type
type(location_type)     :: temp_loc

if ( .not. module_initialized ) call static_init_model

! An interface is provided
interf_provided = .true.

! If first call initialize random sequence
if(first_pert_call) then
   call init_random_seq(random_seq)
   first_pert_call = .false.
endif

do i = 1, get_model_size()
   call get_state_meta_data(i, temp_loc, variable_type)
   if(variable_type == KIND_U_WIND_COMPONENT .or. &
      variable_type == KIND_V_WIND_COMPONENT .or. &
      variable_type == KIND_TEMPERATURE      .or. & 
      variable_type == KIND_DENSITY          .or. & 
      variable_type == KIND_ELECTRON_DENSITY ) then
      pert_state(i) = &
       & random_gaussian(random_seq, state(i), state(i)*0.01)
   else
      pert_state(i) = state(i)
   endif
end do

end subroutine pert_model_state



subroutine prog_var_to_vector(var, x)
!=======================================================================
!
! Copies fields to straight vector
!

type(model_type), intent(in) :: var
real(r8), intent(out)        :: x(:)

integer :: i, j, k, nf, indx

if ( .not. module_initialized ) call static_init_model

indx = 0
loop_longitude: do i = 1, nlon       
   loop_latitude: do j = 1, nlat     
        loop_level: do k = 1, nlev 
            loop_var: do nf = 1, state_num_3d         
                        indx = indx + 1
                        x(indx) = var%vars_3d(i,j,k,nf)
                      enddo loop_var
                    enddo loop_level
                  enddo loop_latitude
                enddo loop_longitude

if(indx /= model_size) then
   write(msgstring, *) 'indx ',indx,' model_size ',model_size,' must be equal '
   call error_handler(E_ERR, 'prog_var_to_vector', msgstring, source, revision, revdate)
endif

end subroutine prog_var_to_vector



subroutine vector_to_prog_var(x, var) 
!==================================================================
! 
! Copies fields from straight vector 
!

real(r8), intent(in)          :: x(:)
type(model_type), intent(out) :: var
integer                       :: i, j, k, nf, indx

if ( .not. module_initialized ) call static_init_model

indx = 0
loop_longitude: do i = 1, nlon       
   loop_latitude: do j = 1, nlat     
        loop_level: do k = 1, nlev 
            loop_var: do nf = 1, state_num_3d         
                        indx = indx + 1
                        var%vars_3d(i,j,k,nf) = x(indx)
                      enddo loop_var
                    enddo loop_level
                  enddo loop_latitude
                enddo loop_longitude

if(indx /= model_size) then
   write(msgstring, *) 'indx ',indx,' model_size ',model_size,' must be equal '
   call error_handler(E_ERR, 'vector_to_prog_var', msgstring, source, revision, revdate)
endif

end subroutine vector_to_prog_var



subroutine init_model_instance(var,valid_time)
!==================================================================
!
! Initializes an instance of TIEGCM model state variables
!

type(model_type),          intent(out) :: var
type(time_type), optional, intent(in)  :: valid_time

if ( .not. module_initialized ) call static_init_model
  
allocate(var%vars_3d(nlon, nlat, nlev, state_num_3d))
  
if (present(valid_time)) then
   var%valid_time = valid_time
else
   var%valid_time = set_time_missing()
endif

end subroutine init_model_instance



subroutine end_model_instance(var)
!==================================================================
!
! Ends an instance of TIEGCM model state variables
!
                                                                                            
type(model_type), intent(inout) :: var

if ( .not. module_initialized ) call static_init_model
                                                                                      
deallocate(var%vars_3d)
                                                                                  
end subroutine end_model_instance



subroutine update_TIEGCM_restart(file_name, var)
!=======================================================================
!
! Updates TIEGCM restart file fields 
!
  
   character (len = *), intent(in)     :: file_name
   type(model_type),    intent(in)     :: var

   integer                             :: ncerr
   integer                             :: restart_id
   integer                             :: dim_id, dim_len
   integer                             :: dim_time_id, dim_time_len
   integer                             :: var_id
   integer, parameter                  :: nmtime = 3
   integer, dimension(nmtime)          :: mtime  ! day, hour, minute 
   integer                             :: utsec, doy !year

   real(r8), dimension(nlon,nlat,nlev) :: TN, TN_NM, UN, UN_NM, VN, VN_NM 
   real(r8), dimension(nlon,nlat,nlev) :: O1, O1_NM
   real(r8), dimension(nlon,nlat,nilev):: NE
   
   integer                             :: nlevm1
   type(time_type)                     :: jan1, tbase
   integer                             :: year, month, day, hour, mins, sec

   if ( .not. module_initialized ) call static_init_model

   nlevm1 = nlev -1

   if( .not. file_exist(file_name)) then
      write(msgstring,*) trim(adjustl(file_name)),' not available.'
      call error_handler(E_ERR,'update_TIEGCM_restart',msgstring,source,revision,revdate)
   endif


   if (do_output()) print *, 'update_TIEGCM_restart: opening restart'
   ncerr = nf90_open( file_name, NF90_WRITE, restart_id )! open with read/write access
   call nc_check(ncerr, 'update_TIEGCM_restart','open')  ! will die if error
   if (do_output()) print *, 'update_TIEGCM_restart: opened with '//trim(nf90_strerror(ncerr))


   !... check for matching dimensions
   call nc_check( nf90_inq_dimid(restart_id, 'lon', dim_id), &
          'update_TIEGCM_restart', 'inq_dimid lon')
   call nc_check( nf90_inquire_dimension(restart_id, dim_id, len=dim_len), &
          'update_TIEGCM_restart', 'inquire_dimension lon')   
   if (dim_len .ne. nlon) then
     write(msgstring, *) trim(file_name), ' dim_lon = ',dim_len, ' DART expects ',nlon
     call error_handler(E_ERR,'update_TIEGCM_restart',msgstring,source,revision,revdate)
   endif


   call nc_check( nf90_inq_dimid(restart_id, 'lat', dim_id), &
          'update_TIEGCM_restart', 'inq_dimid lat')
   call nc_check( nf90_inquire_dimension(restart_id, dim_id, len=dim_len), &
          'update_TIEGCM_restart', 'inquire_dimension lat')   
   if (dim_len .ne. nlat) then
     write(msgstring, *) trim(file_name), ' dim_lat = ',dim_len, ' DART expects ',nlat
     call error_handler(E_ERR,'update_TIEGCM_restart',msgstring,source,revision,revdate)
   endif


   call nc_check( nf90_inq_dimid(restart_id, 'lev', dim_id), &
          'update_TIEGCM_restart', 'inq_dimid lev')
   call nc_check( nf90_inquire_dimension(restart_id, dim_id, len=dim_len), &
          'update_TIEGCM_restart', 'inquire_dimension lev')   
   if (dim_len .ne. nlev) then
     write(msgstring, *) trim(file_name), ' dim_lev = ',dim_len, ' DART expects ',nlev
     call error_handler(E_ERR,'update_TIEGCM_restart',msgstring,source,revision,revdate)
   endif


   call nc_check( nf90_inq_dimid(restart_id, 'ilev', dim_id), &
          'update_TIEGCM_restart', 'inq_dimid ilev')
   call nc_check( nf90_inquire_dimension(restart_id, dim_id, len=dim_len), &
          'update_TIEGCM_restart', 'inquire_dimension ilev')   
   if (dim_len .ne. nilev) then
     write(msgstring, *) trim(file_name), ' dim_ilev = ',dim_len, ' DART expects ',nilev
     call error_handler(E_ERR,'update_TIEGCM_restart',msgstring,source,revision,revdate)
   endif

    
   call nc_check( nf90_inquire(restart_id, unlimitedDimId = dim_time_id), &
          'update_TIEGCM_restart', 'inquire id of unlimited dimension time')
   call nc_check( nf90_inquire_dimension(restart_id, dim_time_id, len=dim_time_len ), &
          'update_TIEGCM_restart', 'inquire_dimension time')
   

!... put variables into TIEGCM array
!... (order: NE, TN, TN_NM, UN, UN_NM, VN, VN_NM, O1, O1_NM)


   NE                  = var%vars_3d(:,:,:,1)
   TN(:,:,1:nlevm1)    = var%vars_3d(:,:,1:nlevm1,2) 
   TN(:,:,  nlev)      = TIEGCM_missing_value        !fill top slot with missing value
   TN_NM(:,:,1:nlevm1) = var%vars_3d(:,:,1:nlevm1,3) 
   TN_NM(:,:,   nlev ) = TIEGCM_missing_value        !fill top slot with missing value

   UN(:,:,1:nlevm1)    = var%vars_3d(:,:,1:nlevm1,4) 
   UN(:,:,   nlev)     = TIEGCM_missing_value        !fill top slot with missing value
   UN_NM(:,:,1:nlevm1) = var%vars_3d(:,:,1:nlevm1,5) 
   UN_NM(:,:,nlev)     = TIEGCM_missing_value        !fill top slot with missing value  

   VN(:,:,1:nlevm1)    = var%vars_3d(:,:,1:nlevm1,6)  
   VN(:,:,  nlev)      = TIEGCM_missing_value        !fill top slot with missing value
   VN_NM(:,:,1:nlevm1) = var%vars_3d(:,:,1:nlevm1,7) 
   VN_NM(:,:,nlev)     = TIEGCM_missing_value        !fill top slot with missing value
   O1                  = var%vars_3d(:,:,:,8)
   O1_NM               = var%vars_3d(:,:,:,9)
  

   call nc_check( nf90_inq_varid(restart_id, 'NE', var_id), &
          'update_TIEGCM_restart', 'inq_varid NE')
   call nc_check( nf90_put_var(restart_id, var_id, values=NE, &
                  start = (/1,1,1,dim_time_len/), count = (/nlon,nlat,nilev,1/)), &
          'update_TIEGCM_restart', 'put_var NE')


   call nc_check( nf90_inq_varid(restart_id, 'TN', var_id), &
          'update_TIEGCM_restart', 'inq_varid TN')
   call nc_check( nf90_put_var(restart_id, var_id, values=TN, &
                  start = (/1,1,1,dim_time_len/), count = (/nlon,nlat,nlev,1/)), &
          'update_TIEGCM_restart', 'put_var TN')


   call nc_check( nf90_inq_varid(restart_id, 'TN_NM', var_id), &
          'update_TIEGCM_restart', 'inq_varid TN_NM')
   call nc_check( nf90_put_var(restart_id, var_id, values=TN_NM, &
                  start = (/1,1,1,dim_time_len/), count = (/nlon,nlat,nlev,1/)), &
          'update_TIEGCM_restart', 'put_var TN_NM')


   call nc_check( nf90_inq_varid(restart_id, 'UN', var_id), &
          'update_TIEGCM_restart', 'inq_varid UN')
   call nc_check( nf90_put_var(restart_id, var_id, values=UN, &
                  start = (/1,1,1,dim_time_len/), count = (/nlon,nlat,nlev,1/)), &
          'update_TIEGCM_restart', 'put_var UN')


   call nc_check( nf90_inq_varid(restart_id, 'UN_NM', var_id), &
          'update_TIEGCM_restart', 'inq_varid UN_NM')
   call nc_check( nf90_put_var(restart_id, var_id, values=UN_NM, &
                  start = (/1,1,1,dim_time_len/), count = (/nlon,nlat,nlev,1/)), &
          'update_TIEGCM_restart', 'put_var UN_NM')

 
   call nc_check( nf90_inq_varid(restart_id, 'VN', var_id), &
          'update_TIEGCM_restart', 'inq_varid VN')
   call nc_check( nf90_put_var(restart_id, var_id, values=VN, &
                  start = (/1,1,1,dim_time_len/), count = (/nlon,nlat,nlev,1/)), &
          'update_TIEGCM_restart', 'put_var VN')


   call nc_check( nf90_inq_varid(restart_id, 'VN_NM', var_id), &
          'update_TIEGCM_restart', 'inq_varid VN_NM')
   call nc_check( nf90_put_var(restart_id, var_id, values=VN_NM, &
                  start = (/1,1,1,dim_time_len/), count = (/nlon,nlat,nlev,1/)), &
          'update_TIEGCM_restart', 'put_var VN_NM')

  
   call nc_check( nf90_inq_varid(restart_id, 'O1', var_id), &
          'update_TIEGCM_restart', 'inq_varid O1')
   call nc_check( nf90_put_var(restart_id, var_id, values=O1, &
                  start = (/1,1,1,dim_time_len/), count = (/nlon,nlat,nlev,1/)), &
          'update_TIEGCM_restart', 'put_var O1')


   call nc_check( nf90_inq_varid(restart_id, 'O1_NM', var_id), &
          'update_TIEGCM_restart', 'inq_varid O1_NM')
   call nc_check( nf90_put_var(restart_id, var_id, values=O1_NM, &
                  start = (/1,1,1,dim_time_len/), count = (/nlon,nlat,nlev,1/)), &
          'update_TIEGCM_restart', 'put_var O1_NM')

!... mtime and year
   call get_date(var%valid_time, year, month, day, hour, mins, sec )
   jan1  = set_date(year,1,1)
   tbase = var%valid_time - jan1    ! total time since the start of the year.

   call get_time(tbase, utsec, doy)

   mtime(1) = doy + 1  ! Have to add January 1 back in
   mtime(2) = hour
   mtime(3) = mins

   if (do_output()) print *, 'update_TIEGCM_restart: mtime (doy/hour/minute):', mtime  
   call nc_check( nf90_inq_varid(restart_id, 'mtime', var_id), &
          'update_TIEGCM_restart','inq_varid mtime')
   call nc_check( nf90_put_var( restart_id, var_id, values=mtime, & 
                  start = (/1,dim_time_len/), count = (/nmtime,1/)), &
          'update_TIEGCM_restart','get_var mtime')

   if (do_output()) print *, 'update_TIEGCM_restart: year:', year  
   call nc_check( nf90_inq_varid(restart_id, 'year', var_id), &
          'update_TIEGCM_restart','inq_varid year')
   call nc_check( nf90_put_var( restart_id, var_id, values=year, &
                  start = (/dim_time_len/)) , &
          'update_TIEGCM_restart','get_var year')
    
   call nc_check( nf90_sync(restart_id), 'update_TIEGCM_restart', 'sync')

   call nc_check( nf90_close(restart_id), 'update_TIEGCM_restart', 'close')


end subroutine update_TIEGCM_restart



subroutine read_TIEGCM_restart(file_name, var, model_time)
!=======================================================================
!
! Read TIEGCM restart file fields 
!

   character (len = *),           intent(in) :: file_name
   type(model_type),              intent(out):: var
   type(time_type),               intent(out):: model_time

   integer                                   :: ncerr
   integer                                   :: restart_id 
   integer                                   :: dim_id, dim_len
   integer                                   :: dim_time_id, dim_time_len
   integer                                   :: var_Vtmp_id, var_mtime_id, var_year_id
   real(r8), dimension(nlon,nlat,nlev)       :: TN, TN_NM, UN, UN_NM, VN, VN_NM  
   real(r8), dimension(nlon,nlat,nlev)       :: O1, O1_NM
   real(r8), dimension(nlon,nlat,nilev)      :: NE
   integer,  dimension(:), allocatable       :: yeartmp
   integer,  dimension(:,:), allocatable     :: mtimetmp
   integer,  parameter                       :: nmtime = 3
   integer,  dimension(nmtime)               :: mtime  ! day, hour, minute 
   integer                                   :: year, utsec, doy
   integer                                   :: nlevm1

   if ( .not. module_initialized ) call static_init_model

   nlevm1 = nlev - 1

   if( .not. file_exist(file_name)) then
      write(msgstring,*)trim(file_name)//' does not exist.'
      call error_handler(E_ERR,'read_TIEGCM_restart',msgstring,source,revision,revdate)
   endif

   if (do_output()) print *, 'read_TIEGCM_restart:reading restart:', file_name

   call nc_check( nf90_open( file_name, NF90_NOWRITE, restart_id ), &
                                 'read_TIEGCM_restart', 'open')
   
!... check for matching dimensions
   call nc_check( nf90_inq_dimid(restart_id, 'lon', dim_id), &
          'read_TIEGCM_restart', 'inq_dimid lon')
   call nc_check( nf90_inquire_dimension(restart_id, dim_id, len=dim_len), &
          'read_TIEGCM_restart', 'inquire_dimension lon')   
   if (dim_len .ne. nlon) then
     write(msgstring, *) trim(file_name), ' dim_lon = ',dim_len, ' DART expects ',nlon
     call error_handler(E_ERR,'read_TIEGCM_restart',msgstring,source,revision,revdate)
   endif

   call nc_check( nf90_inq_dimid(restart_id, 'lat', dim_id), &
          'read_TIEGCM_restart', 'inq_dimid lat')
   call nc_check( nf90_inquire_dimension(restart_id, dim_id, len=dim_len), &
          'read_TIEGCM_restart', 'inquire_dimension lat')   
   if (dim_len .ne. nlat) then
     write(msgstring, *) trim(file_name), ' dim_lat = ',dim_len, ' DART expects ',nlat
     call error_handler(E_ERR,'read_TIEGCM_restart',msgstring,source,revision,revdate)
   endif

   call nc_check( nf90_inq_dimid(restart_id, 'lev', dim_id), &
          'read_TIEGCM_restart', 'inq_dimid lev')
   call nc_check( nf90_inquire_dimension(restart_id, dim_id, len=dim_len), &
          'read_TIEGCM_restart', 'inquire_dimension lev')   
   if (dim_len .ne. nlev) then
     write(msgstring, *) trim(file_name), ' dim_lev = ',dim_len, ' DART expects ',nlev
     call error_handler(E_ERR,'read_TIEGCM_restart',msgstring,source,revision,revdate)
   endif

   call nc_check( nf90_inq_dimid(restart_id, 'ilev', dim_id), &
          'read_TIEGCM_restart', 'inq_dimid ilev')
   call nc_check( nf90_inquire_dimension(restart_id, dim_id, len=dim_len), &
          'read_TIEGCM_restart', 'inquire_dimension ilev')   
   if (dim_len .ne. nilev) then
     write(msgstring, *) trim(file_name), ' dim_ilev = ',dim_len, ' DART expects ',nilev
     call error_handler(E_ERR,'read_TIEGCM_restart',msgstring,source,revision,revdate)
   endif
    
   call nc_check(nf90_inquire(restart_id, unlimitedDimId = dim_time_id), &
          'read_TIEGCM_restart', 'inquire id of unlimited dimension time')
   call nc_check(nf90_inquire_dimension(restart_id, dim_time_id, len=dim_time_len ), &
          'read_TIEGCM_restart', 'inquire_dimension time')

!... NE
   call nc_check(nf90_inq_varid(restart_id, 'NE', var_Vtmp_id),       &
          'read_TIEGCM_restart', 'inq_varid NE')
   call nc_check(nf90_get_var(restart_id, var_Vtmp_id, values=NE,     &
                            start = (/ 1, 1, 1, dim_time_len /),      &
                            count = (/ nlon, nlat, nilev, 1 /)),      &
                            'read_TIEGCM_restart', 'get_var NE') 

!... TN
   call nc_check(nf90_inq_varid(restart_id, 'TN', var_Vtmp_id),       &
          'read_TIEGCM_restart', 'inq_varid TN')                            
   call nc_check(nf90_get_var(restart_id, var_Vtmp_id, values=TN,     &
                            start = (/ 1, 1, 1, dim_time_len /),      &
                            count = (/ nlon, nlat, nlev, 1 /)),       &
                            'read_TIEGCM_restart', 'get_var TN') 

!... TN_NM
   call nc_check(nf90_inq_varid(restart_id, 'TN_NM', var_Vtmp_id),    &
          'read_TIEGCM_restart', 'inq_varid TN_NM')
   call nc_check(nf90_get_var(restart_id, var_Vtmp_id, values=TN_NM,  &
                            start = (/ 1, 1, 1, dim_time_len /),      &
                            count = (/ nlon, nlat, nlev, 1 /)),       &
                            'read_TIEGCM_restart', 'get_var TN_NM') 

!... UN
   call nc_check(nf90_inq_varid(restart_id, 'UN', var_Vtmp_id),       &
          'read_TIEGCM_restart', 'inq_varid UN')
   call nc_check(nf90_get_var(restart_id, var_Vtmp_id, values=UN,     &
                            start = (/ 1, 1, 1, dim_time_len /),      &
                            count = (/ nlon, nlat, nlev, 1 /)),       &
          'read_TIEGCM_restart', 'get_var UN') 
 
!... UN_NM
   call nc_check(nf90_inq_varid(restart_id, 'UN_NM', var_Vtmp_id),    &
          'read_TIEGCM_restart', 'inq_varid UN_NM')
   call nc_check(nf90_get_var(restart_id, var_Vtmp_id, values=UN_NM,  &
                            start = (/ 1, 1, 1, dim_time_len /),      &
                            count = (/ nlon, nlat, nlev, 1 /)),       &
          'read_TIEGCM_restart', 'get_var UN_NM') 

!... VN
   call nc_check(nf90_inq_varid(restart_id, 'VN', var_Vtmp_id),       &
          'read_TIEGCM_restart', 'inq_varid VN')
   call nc_check(nf90_get_var(restart_id, var_Vtmp_id, values=VN,     &
                            start = (/ 1, 1, 1, dim_time_len /),      &
                            count = (/ nlon, nlat, nlev, 1 /)),       &                            
          'read_TIEGCM_restart', 'get_var VN') 

!... VN_NM
   call nc_check(nf90_inq_varid(restart_id, 'VN_NM', var_Vtmp_id),    &
          'read_TIEGCM_restart', 'inq_varid VN_NM')
   call nc_check(nf90_get_var(restart_id, var_Vtmp_id, values=VN_NM,  &
                            start = (/ 1, 1, 1, dim_time_len /),      &
                            count = (/ nlon, nlat, nlev, 1 /)),       &     
          'read_TIEGCM_restart', 'get_var VN_NM') 

!... O1
   call nc_check(nf90_inq_varid(restart_id, 'O1', var_Vtmp_id),       &
          'read_TIEGCM_restart', 'inq_varid O1')
   call nc_check(nf90_get_var(restart_id, var_Vtmp_id, values=O1,     &
                            start = (/ 1, 1, 1, dim_time_len /),      &
                            count = (/ nlon, nlat, nlev, 1 /)),       &     
          'read_TIEGCM_restart', 'get_var O1') 

!... O1_NM
   call nc_check(nf90_inq_varid(restart_id, 'O1_NM', var_Vtmp_id),    &
          'read_TIEGCM_restart', 'inq_varid O1_NM')
   call nc_check(nf90_get_var(restart_id, var_Vtmp_id, values=O1_NM,  &
                            start = (/ 1, 1, 1, dim_time_len /),      &
                            count = (/ nlon, nlat, nlev, 1 /)),       &     
          'read_TIEGCM_restart', 'get_var O1_NM') 

!... get mtime
   call nc_check(nf90_inq_dimid(restart_id, 'mtimedim', dim_id), &
          'read_TIEGCM_restart', 'inq_dimid mtimedim')
   call nc_check(nf90_inquire_dimension(restart_id, dim_id, len=dim_len), &
          'read_TIEGCM_definition', 'inquire_dimension mtimedim')  
   if (dim_len .ne. nmtime) then
     write(msgstring, *) trim(file_name), ' mtimedim = ',dim_len, ' DART expects ', nmtime
     call error_handler(E_ERR,'read_TIEGCM_restart',msgstring,source,revision,revdate)
   endif

   allocate(mtimetmp(dim_len, dim_time_len))
   call nc_check(nf90_inq_varid(restart_id, 'mtime', var_mtime_id), &
          'read_TIEGCM_restart', 'inq_varid mtime')
   call nc_check(nf90_get_var(restart_id, var_mtime_id, values=mtimetmp), &
          'read_TIEGCM_restart', 'get_var mtime')    
   mtime = mtimetmp(:,dim_time_len)
   deallocate(mtimetmp)

!... get year
   allocate(yeartmp(dim_time_len))
   call nc_check(nf90_inq_varid(restart_id, 'year', var_year_id), &
          'read_TIEGCM_restart', 'inq_varid year')
   call nc_check(nf90_get_var(restart_id, var_year_id, values=yeartmp), &
          'read_TIEGCM_restart', 'get_var year')    
   year = yeartmp(dim_time_len)
   deallocate(yeartmp)

!... close the file
   call nc_check(nf90_close( restart_id),'read_TIEGCM_restart','close')


! ... Now we want to convert the year/doy/hour/minute to a dart time.
! ... We start by finding the dart time of the year and adding the rest to it.

   if (do_output()) print *, 'read_TIEGCM_restart: mtime (doy/hour/minute) and year:', mtime, year

   doy   =  mtime(1)
   utsec = (mtime(2)*60 + mtime(3))*60

   model_time = set_time(utsec, doy-1) + set_date(year, 1, 1)  ! Jan 1 of whatever year.

   var%valid_time = model_time
 
  if (do_output()) call print_time(model_time, str=" read_TIEGCM_restart: model_time ")   
  if (do_output()) call print_date(model_time, str=" read_TIEGCM_restart: model_date ")   

!... fill DART state vector with TIEGCM variables
!    (order: NE  TN TN_NM UN UN_NM VN VN_NM  O1 O1_NM)
   var%vars_3d(:,:,:,1)        = NE(:,:,:)

   var%vars_3d(:,:,1:nlevm1,2) = TN(:,:,1:nlevm1)    
   var%vars_3d(:,:,    nlev,2) = TN(:,:,  nlevm1)    !fill top slot with values at nlev-1  

   var%vars_3d(:,:,1:nlevm1,3) = TN_NM(:,:,1:nlevm1) 
   var%vars_3d(:,:,    nlev,3) = TN_NM(:,:,  nlevm1) !fill top slot with values at nlev-1

   var%vars_3d(:,:,1:nlevm1,4) = UN(:,:,1:nlevm1)    
   var%vars_3d(:,:,    nlev,4) = UN(:,:,  nlevm1)    !fill top slot with values at nlev-1

   var%vars_3d(:,:,1:nlevm1,5) = UN_NM(:,:,1:nlevm1) 
   var%vars_3d(:,:,    nlev,5) = UN_NM(:,:,  nlevm1) !fill top slot with values at nlev-1

   var%vars_3d(:,:,1:nlevm1,6) = VN(:,:,1:nlevm1)    
   var%vars_3d(:,:,    nlev,6) = VN(:,:,  nlevm1)    !fill top slot with values at nlev-1

   var%vars_3d(:,:,1:nlevm1,7) = VN_NM(:,:,1:nlevm1) 
   var%vars_3d(:,:,    nlev,7) = VN_NM(:,:,  nlevm1) !fill top slot with values at nlev-1

   var%vars_3d(:,:,:,8)        = O1
   var%vars_3d(:,:,:,9)        = O1_NM

end subroutine read_TIEGCM_restart



subroutine read_TIEGCM_definition(file_name)
!=======================================================================
!
! Read TIEGCM grid definition and Geopotential from a tiegcm restart file
!
  
   character (len = *), intent(in)           :: file_name

   integer                                   :: ncerr
   integer                                   :: restart_id
   integer                                   :: var_lon_id, var_lat_id, var_lev_id, var_ilev_id
   integer                                   :: dim_lon_id, dim_lat_id, dim_lev_id, dim_ilev_id
   integer                                   :: dim_time_id, dim_time_len, missing_value_len
   integer                                   :: var_p0_id
   integer                                   :: l
   real(r8)                                  :: p0, lon_tmp


   if( .not. file_exist(file_name)) then
     write(msgstring,*) trim(adjustl(file_name)),' not available.'
     call error_handler(E_ERR,'read_TIEGCM_definition',msgstring,source,revision,revdate)
   endif


   if (do_output()) print *, 'read_TIEGCM_definition:reading restart:', file_name
   ncerr = nf90_open(file_name, NF90_NOWRITE, restart_id)
   call nc_check(ncerr, 'read_TIEGCM_definition', 'open')   
   if (do_output()) print *, 'read_TIEGCM_definition:opened with '//trim(nf90_strerror(ncerr))


   call nc_check(nf90_inq_dimid(restart_id, 'lon', dim_lon_id), &
          'read_TIEGCM_definition', 'inq_dimid lon')
   call nc_check(nf90_inquire_dimension(restart_id, dim_lon_id, len=nlon), &
          'read_TIEGCM_definition', 'inquire_dimension lon')
   allocate(lons(nlon))
   call nc_check(nf90_inq_varid(restart_id, 'lon', var_lon_id), &
          'read_TIEGCM_definition', 'inq_varid lon')
   call nc_check(nf90_get_var(restart_id, var_lon_id, values=lons), &
          'read_TIEGCM_definition', 'get_var lon') 

   do l = 1, nlon
      lon_tmp = lons(l)
      if (lon_tmp < 0) lons(l) = lons(l) + 360 ! DART [0, 360] TIEGCM [-180, 180]     
   enddo

   call nc_check(nf90_inq_dimid(restart_id, 'lat', dim_lat_id), &
          'read_TIEGCM_definition', 'inq_dimid lat')
   call nc_check(nf90_inquire_dimension(restart_id, dim_lat_id, len=nlat), &
          'read_TIEGCM_definition', 'inquire_dimension lat')    
   allocate(lats(nlat))
   call nc_check(nf90_inq_varid(restart_id, 'lat', var_lat_id), &
          'read_TIEGCM_definition', 'inq_varid lat')
   call nc_check(nf90_get_var(restart_id, var_lat_id, values=lats), &
          'read_TIEGCM_definition', 'get_var lat') 


   call nc_check(nf90_inq_varid(restart_id, 'p0', var_p0_id), &
          'read_TIEGCM_definition', 'inq_varid p0')
   call nc_check(nf90_get_var(restart_id, var_p0_id, values=p0), &
          'read_TIEGCM_definition', 'get_var p0') 


   call nc_check(nf90_inq_dimid(restart_id, 'lev', dim_lev_id), &
          'read_TIEGCM_definition', 'inq_dimid lev')
   call nc_check(nf90_inquire_dimension(restart_id, dim_lev_id, len=nlev), &
          'read_TIEGCM_definition', 'inquire_dimension lev')    
   allocate(levs(nlev))
   call nc_check(nf90_inq_varid(restart_id, 'lev', var_lev_id), &
          'read_TIEGCM_definition', 'inq_varid lev')
   call nc_check(nf90_get_var(restart_id, var_lev_id, values=levs), &
          'read_TIEGCM_definition', 'get_var lev') 

   levs = p0 * exp(-levs) * 100.0_r8 ![Pa] = 100* [millibars] = 100* [hPa]

   call nc_check(nf90_inq_dimid(restart_id, 'ilev', dim_ilev_id), &
          'read_TIEGCM_definition', 'inq_dimid ilev')
   call nc_check(nf90_inquire_dimension(restart_id, dim_ilev_id, len=nilev), &
          'read_TIEGCM_definition', 'inquire_dimension ilev')      
   allocate(ilevs(nilev)) 
   call nc_check(nf90_inq_varid(restart_id, 'ilev', var_ilev_id), &
          'read_TIEGCM_definition', 'inq_varid ilev')
   call nc_check(nf90_get_var(restart_id, var_ilev_id, values=ilevs), &
          'read_TIEGCM_definition', 'get_var ilev') 

   ilevs = p0 * exp(-ilevs) * 100.0_r8 ! [Pa] = 100* [millibars] = 100* [hPa]

   if (nlev .ne. nilev) then
     write(msgstring, *) ' nlev = ',nlev,' nilev = ',nilev, 'are different; DART assumes them to be the same'
     call error_handler(E_ERR,'read_TIEGCM_definition',msgstring,source,revision,revdate)
   endif
   
   call nc_check(nf90_inquire_attribute(restart_id, nf90_global, 'missing_value', len = missing_value_len), &
          'read_TIEGCM_definition', 'inquire global attribute named missing_value')
   if (missing_value_len .ne. 1) then
     write(msgstring, *) ' global attribute missing_value length is ', missing_value_len, ' DART expects 1'
     call error_handler(E_ERR,'read_TIEGCM_definition',msgstring,source,revision,revdate)
   endif
   
   call nc_check(nf90_get_att(restart_id, nf90_global, 'missing_value', TIEGCM_missing_value), &
          'read_TIEGCM_definition', 'get_att global attribute named missing_value')
   
   call nc_check(nf90_close(restart_id),'read_TIEGCM_definition', 'close')   


end subroutine read_TIEGCM_definition



subroutine read_TIEGCM_secondary(file_name)
!=======================================================================
!
! Read TIEGCM Geopotential from a tiegcm secondary output file
!
  
   character (len = *), intent(in)           :: file_name
   integer                                   :: restart_id
   integer                                   :: dim_time_id, dim_time_len
   integer                                   :: var_ZGtmp_id
   real(r8), dimension(:,:,:,:), allocatable :: ZGtmp


   if( .not. file_exist(file_name)) then
     write(msgstring,*) trim(adjustl(file_name)),' not available.'
     call error_handler(E_ERR,'read_TIEGCM_secondary',msgstring,source,revision,revdate)
   endif

   if (do_output()) print *, 'read_TIEGCM_secondary:reading restart:', file_name

   call nc_check(nf90_open(file_name, NF90_NOWRITE, restart_id), &
                           'read_TIEGCM_secondary', 'open')   

   call nc_check(nf90_inquire(restart_id, unlimitedDimId = dim_time_id), &
          'read_TIEGCM_secondary', 'inquire id of unlimited dimension time')

   call nc_check(nf90_inquire_dimension(restart_id, dim_time_id, len=dim_time_len ), &
          'read_TIEGCM_secondary', 'inquire_dimension time')

   allocate(ZGtmp(nlon,nlat,nilev,dim_time_len))

   call nc_check(nf90_inq_varid(restart_id, 'ZG', var_ZGtmp_id), &
          'read_TIEGCM_secondary', 'inq_varid ZG')

   call nc_check(nf90_get_var(restart_id, var_ZGtmp_id, values=ZGtmp), &
          'read_TIEGCM_secondary', 'get_var ZG') 

   allocate(ZGtiegcm(nlon,nlat,nilev))

   ZGtiegcm = ZGtmp(:,:,:,dim_time_len) ! [cm]

   deallocate(ZGtmp)
 
   call nc_check(nf90_close(restart_id),'read_TIEGCM_secondary', 'close')   

end subroutine read_TIEGCM_secondary



subroutine ens_mean_for_model(ens_mean)
!===================================================================
!
!

real(r8), intent(in) :: ens_mean(:)

if ( .not. module_initialized ) call static_init_model
 
end subroutine ens_mean_for_model


!===================================================================
! End of model_mod
!===================================================================
end module model_mod
