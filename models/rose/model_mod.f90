! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

!-----------------------------------------------------------------------
!
!     Interface for ROSE model
!
!-----------------------------------------------------------------------

! DART Modules 
use        types_mod, only : r8, digits12, pi
use time_manager_mod, only : time_type, set_time, print_time, get_time, &
                             operator(<), operator(>), operator(+), &
                             operator(-), operator(/), operator(*), &
                             operator(==), operator(/=), set_time_missing
use     location_mod, only : location_type, get_close_maxdist_init, &
                             get_close_obs_init, get_close_obs, &
                             set_location, get_location, query_location, &
                             get_dist, vert_is_height
use    utilities_mod, only : file_exist, open_file, close_file, &       
                             error_handler, E_ERR, E_MSG, E_WARN, nmlfileunit, &
                             do_output, find_namelist_in_file, check_namelist_read, &
                             do_nml_file, do_nml_term
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use     obs_kind_mod, only : QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT,&
                             QTY_TEMPERATURE
use    netcdf_utilities_mod, only : nc_check
use typesizes
use netcdf

! ROSE Modules
use params, only : nx, ny, nz, ntime

!-----------------------------------------------------------------------

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

! ROSE-specific routines 
public :: model_type, &
          prog_var_to_vector, &
          vector_to_prog_var, &
          read_ROSE_restart, &
          update_ROSE_restart, &
          update_ROSE_namelist, &
          read_ROSE_tref,&
          init_model_instance, &
          end_model_instance

!-----------------------------------------------------------------------

type model_type
  real(r8), pointer :: vars_3d(:,:,:,:)
  type(time_type)   :: valid_time
end type model_type

integer, parameter :: TYPE_local_U0 = 0, &
                      TYPE_local_V0 = 1, &
                      TYPE_local_T0 = 2, &
                      TYPE_local_U  = 3, &
                      TYPE_local_V  = 4, &
                      TYPE_local_T  = 5

!----------------------------------------------------------------------

! ROSE vertical grid variables
real(r8)                   :: dz = 2100.0_r8, zbot = 16800.0_r8
real, dimension (nz)       :: tref !! auxiliary variable needed for H                           

! Global storage for describing ROSE model class
integer :: model_size 
type(time_type) :: Time_step_ROSE

! Arrays to store lats, lons
real(r8) :: lons(nx), lats(ny), levs(nz)

! Random sequence and init for pert_model_state
logical :: first_pert_call = .true.
type(random_seq_type)   :: random_seq

!----------------------------------------------------------------------
! Nameslist variables with default values follow
! Namelist variables for defining state vector, and default values
integer :: state_num_3d = 6             ! # of 3d fields to read from file
namelist /model_nml/ state_num_3d

logical :: output_prog_diag = .false.
character(len=128)   :: input_dir = '../input_current/'
character(len=50)   :: out_dir   = '/ptmp/tmatsuo/rose/'
character(len=30)   :: restart_file = 'dyn_restart_001-1999.nc'
real(kind=r8)       :: amp_tune = 1.
real(kind=r8)       :: pha_tune = 0.
real(kind=digits12) :: target_time = 0.125 ! [hr] 
integer             :: ens_element = 1

namelist /ROSE_NML/ target_time, &
                    input_dir, out_dir, restart_file,&
                    output_prog_diag, &
                    amp_tune, pha_tune, &
                    ntime, ens_element

integer, parameter :: TYPE_U0 = 0, TYPE_V0 = 1, TYPE_T0 = 2, & 
                      TYPE_U = 3, TYPE_V = 4, TYPE_T = 5, & 
                      TYPE_Q_H = 6, TYPE_Q_OH = 7, TYPE_Q_O = 8

!----------------------------------------------------------------------
! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

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

integer  :: iunit, io, i, j, k
integer  :: Time_step_seconds, Time_step_days
integer  :: seconds_of_day = 86400
real(r8) :: d_lat, d_lon
real(r8) :: z_m
real(r8) :: dz = 2100.0_r8, zbot = 16800.0_r8

! Read the namelist ROSE_NML from the file rose.nml
call find_namelist_in_file("rose.nml", "ROSE_NML", iunit)
read(iunit, nml = ROSE_NML, iostat = io)
call check_namelist_read(iunit, io, "ROSE_NML")

if (do_nml_file()) write(nmlfileunit, nml=ROSE_NML)
if (do_nml_term()) write(     *     , nml=ROSE_NML)

! Read the namelist entry for model_mod from file input.nml
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Compute overall model size and put in global storage
model_size = nx * ny * nz * state_num_3d

! Set the model minimum time step from the namelist seconds and days input
! ROSE control variable "ntime" specify the number of steps per hour

Time_step_seconds = nint(3600.0_r8 / real(ntime))
Time_step_days = 0
write(*, *) 'time step secs days ' , Time_step_seconds, Time_step_days
if (Time_step_seconds > seconds_of_day) then
   Time_step_days = floor(real(Time_step_seconds) / real(seconds_of_day))
   Time_step_seconds = mod(Time_step_seconds, Time_step_days * seconds_of_day) 
endif
Time_step_ROSE = set_time(Time_step_seconds, Time_step_days)

if (do_output()) call print_time(Time_step_ROSE,'ROSE time step')

! lon: long_name = "geographic longitude", units = "degrees" ;
d_lon = 360.0_r8/real(nx)
do i = 1, nx
   lons(i) = (i-1)*d_lon
enddo

! lat: long_name = "geographic latitude",  units = "degrees" ;
d_lat = 180.0_r8/real(ny)
do j = 1, ny
   lats(j) = -90 + d_lat/2. + (j-1)*d_lat
enddo

! Pressure levels ! from msetfix.f
! 1013 [mb] exp(-z/7km) with a fixed scale height (7km) 
! (2.1 km increment from 16.8 km)
do k = 1, nz
z_m = (k-1)*dz + zbot
levs(k) = 1013._r8*exp(-z_m/7.e3)   ![hPa]
enddo

! Read in rose_restart file for tref (global mean temperature) 
! for computation of total temperature needed for H
call read_ROSE_tref(restart_file)

module_initialized = .true.

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
! compute a timestep, for instance for radiation compuations.
! This interface is only called if the namelist parameter
! async is set to 0 in perfect_model_obs of filter or if the 
! program integrate_model is to be used to advance the model
! state as a separate executable. If one of these options
! is not going to be used (the model will only be advanced as
! a separate model-specific executable), this can be a 
! NULL INTERFACE.

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

if ( .not. module_initialized ) call static_init_model

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

! For now, just set to 0
time = set_time(0, 0)
if ( .not. module_initialized ) call static_init_model

end subroutine init_time



subroutine model_interpolate(x, location, itype, obs_val, istatus)
!------------------------------------------------------------------
!
! Given a state vector, a location, and a model state variable type,
! interpolates the state variable field`to that location and returns
! the value in obs_val. The istatus variable should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is a model specific integer that specifies the type of field (for
! instance temperature, zonal wind component, etc.). In low order
! models that have no notion of types of variables, this argument can
! be ingored. For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observerd), this can be a NULL INTERFACE.

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

integer  :: i, vstatus, which_vert
integer  :: lat_below, lat_above, lon_below, lon_above
real(r8) :: lon_fract, temp_lon, lat_fract
real(r8) :: lon, lat, height, lon_lat_lev(3)
real(r8) :: bot_lon, top_lon, delta_lon, bot_lat, top_lat
real(r8) :: val(2,2), a(2)
if ( .not. module_initialized ) call static_init_model

! Default for successful return
istatus = 0
vstatus = 0

! Get the position, determine if it is model level or pressure in vertical
lon_lat_lev = get_location(location)
lon = lon_lat_lev(1); lat = lon_lat_lev(2);
if(vert_is_height(location)) then
   height = lon_lat_lev(3)
else
   which_vert = nint(query_location(location))
   write(msgstring,*) 'vertical coordinate type:',which_vert,' cannot be handled'
   call error_handler(E_ERR,'model_interpolate',msgstring,source,revision,revdate)
endif

! Get lon and lat grid specs, nx, ny, nz from ROSE param mod
   bot_lon = lons(1)
   top_lon = lons(nx)
   delta_lon = lons(2) - lons(1)
   bot_lat = lats(1)
   top_lat = lats(ny)

! Compute bracketing lon indices
if(lon >= bot_lon .and. lon <= top_lon) then
   lon_below = int((lon - bot_lon) / delta_lon) + 1
   lon_above = lon_below + 1
   lon_fract = (lon - ((lon_below - 1) * delta_lon + bot_lon)) / delta_lon
else
! At wraparound point
   lon_below = nx
   lon_above = 1
   if(lon < bot_lon) then
      temp_lon = lon + 360.0
   else
      temp_lon = lon
   endif
   lon_fract = (temp_lon - top_lon) / delta_lon
endif

! Next, compute neighboring lat rows
! NEED TO BE VERY CAREFUL ABOUT POLES; WHAT'S BEING DONE MAY BE WRONG
! Inefficient search used for latitudes in Gaussian grid. Might want to speed up.
if(lat >= bot_lat .and. lat <= top_lat) then

   do i = 2, ny
      if(lat <= lats(i)) then
         lat_above = i
         lat_below = i - 1
         lat_fract = (lat - lats(lat_below)) / (lats(lat_above) - lats(lat_below))
         goto 20
      end if
   end do
else if(lat <= bot_lat) then
! South of bottom lat NEED TO DO BETTER: NOT REALLY REGULAR
   lat_below = 1
   lat_above = 1
   lat_fract = 1.
else
! North of top lat NEED TO DO BETTER: NOT REALLY REGULAR
   lat_below = ny
   lat_above = ny
   lat_fract = 1.
endif

20 continue

! Now, need to find the values for the four corners
                     call get_val(val(1, 1), x, lon_below, lat_below, height, itype, vstatus)
   if (vstatus /= 1) call get_val(val(1, 2), x, lon_below, lat_above, height, itype, vstatus)
   if (vstatus /= 1) call get_val(val(2, 1), x, lon_above, lat_below, height, itype, vstatus)
   if (vstatus /= 1) call get_val(val(2, 2), x, lon_above, lat_above, height, itype, vstatus)


! kdr Guam;s
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
real(r8), intent(in)  :: x(:), height
integer, intent(in)   :: lon_index, lat_index, obs_kind
integer, intent(out)  :: istatus

integer               :: var_type
integer               :: k, lev_top, lev_bottom
real(r8)              :: zgrid
real(r8)              :: val_top, val_bottom, frac_lev

! No errors to start with
istatus = 0

! Pressure levels ! from msetfix.f
! 1013 [mb] exp(-z/7km) with a fixed scale height (7km)
! (2.5 km increment from 17.5 km)
!levs(k) = 1013._r8*exp(-z_m/7.e3)   ![hPa]

height_loop:do k = 1, nz
   zgrid = (k-1)*dz + zbot
   if (height <= zgrid) then
      lev_top = k
      frac_lev = (zgrid - height)/dz
      exit height_loop
   endif
enddo height_loop

lev_bottom = lev_top -1

if (obs_kind == QTY_TEMPERATURE) then

  var_type   = TYPE_local_T
! val_top = x(get_index(lat_index, lon_index, lev_top, var_type))
! val_bottom = x(get_index(lat_index, lon_index, lev_bottom, var_type))

  val_top = x(get_index(lat_index,lon_index,lev_top,var_type)) &
          + tref(lev_top)
  val_bottom = x(get_index(lat_index,lon_index,lev_bottom,var_type)) &
          + tref(lev_bottom)

elseif (obs_kind == QTY_U_WIND_COMPONENT) then

  var_type = TYPE_local_U
  val_top = x(get_index(lat_index, lon_index, lev_top, var_type))
  val_bottom =  x(get_index(lat_index, lon_index, lev_bottom, var_type))

elseif (obs_kind == QTY_V_WIND_COMPONENT) then

  var_type = TYPE_local_V
  val_top = x(get_index(lat_index, lon_index, lev_top, var_type))
  val_bottom =  x(get_index(lat_index, lon_index, lev_bottom, var_type))

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

get_index = 1 + var_type + (lev_index - 1)*state_num_3d      &
                         + (lat_index - 1)*state_num_3d*nz   &
                         + (lon_index - 1)*state_num_3d*nz*ny

end function get_index



function get_model_time_step()
!------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.
!
! Limited by ROSE's fixed time step (also by 2-step restart file)

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step =  Time_step_ROSE

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

integer  :: indx, num_per_col, col_num, col_elem, lon_index,&
            lat_index, lev_index
real(r8) :: lon, lat, lev
integer  :: local_var_type, var_type_temp

if ( .not. module_initialized ) call static_init_model

! Easier to compute with a 0 to size - 1 index
indx = index_in - 1

! Compute number of items per column
num_per_col = nz * state_num_3d 

! What column is this index in
col_num  = indx / num_per_col 
col_elem = indx - col_num * num_per_col

! What lon and lat index for this column
lon_index = col_num / ny      ! ny is number of ROSE latitude grid
lat_index = col_num - lon_index * ny

! Get actual lon lat values from static_init arrays ???
lon = lons(lon_index + 1)
lat = lats(lat_index + 1)

! Now figure out which beast in column this is
lev_index = col_elem / state_num_3d 
lev = levs(lev_index + 1)

! Find which var_type this element is
var_type_temp = mod(col_elem, state_num_3d )
if(var_type_temp == 0) then
  local_var_type = QTY_U_WIND_COMPONENT
else if(var_type_temp == 1) then
  local_var_type = QTY_V_WIND_COMPONENT
else if(var_type_temp == 2) then
  local_var_type = QTY_TEMPERATURE
else if(var_type_temp == 3) then
  local_var_type = QTY_U_WIND_COMPONENT
else if(var_type_temp == 4) then
  local_var_type = QTY_V_WIND_COMPONENT
else if(var_type_temp == 5) then
  local_var_type = QTY_TEMPERATURE
! else if(var_type_temp == 6) then
!  local_var_type = TYPE_Q_H
!else if(var_type_temp == 7) then 
!  local_var_type = TYPE_Q_OH
!else
!  local_var_type = TYPE_Q_O
endif

!write(*, '(1x,3(f6.2,1x),i3)') lon, lat, lev, local_var_type

location = set_location(lon, lat, lev, 2)  ! 2 == pressure (hPa)

! If the type is wanted, return it
if(present(var_type)) var_type = local_var_type

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storgae, etc.

if ( .not. module_initialized ) call static_init_model

end subroutine end_model






function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH Jan 24 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
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

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID


integer :: MemberDimID, TimeDimID
integer :: lonDimID, latDimID, levDimID
integer :: lonVarID, latVarID, levVarID
integer :: u1VarID, v1VarID, t1VarID, uVarID, vVarID, tVarID
! integer :: qnHVarID, qnOHVarID, qnOVarID

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

!------------------------------------------------------------------

if ( .not. module_initialized ) call static_init_model

ierr = 0  ! assume normal termination 

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!--------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID),'nc_write_model_atts','inquire')
call nc_check(nf90_Redef(ncFileID),'nc_write_model_atts','redef')

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies
!-------------------------------------------------------------------------------
                                                                                                           
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID),'nc_write_model_atts')
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID),'nc_write_model_atts')
                                                                                                           
if ( TimeDimID /= unlimitedDimId ) then
   write(msgstring,*)'Time Dimension ID ',TimeDimID,' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', msgstring, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------
                                                                                                           
call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source",source),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate",revdate),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","ROSE"),'nc_write_model_atts')

!-------------------------------------------------------------------------------
! Define the dimensions IDs
!-------------------------------------------------------------------------------

call nc_check(nf90_def_dim(ncid=ncFileID, name="lon",   len = nx,   dimid =   lonDimID),'nc_write_model_atts')
call nc_check(nf90_def_dim(ncid=ncFileID, name="lat",   len = ny,   dimid =   latDimID),'nc_write_model_atts')
call nc_check(nf90_def_dim(ncid=ncFileID, name="lev",   len = nz,   dimid =   levDimID),'nc_write_model_atts')

!-------------------------------------------------------------------------------
! Create the (empty) Variables and the Attributes
!-------------------------------------------------------------------------------

call nc_check(nf90_def_var(ncFileID, name="lon", xtype=nf90_double, &
                                                dimids=lonDimID, varid=lonVarID),'nc_write_model_atts' ) 
call nc_check(nf90_put_att(ncFileID, lonVarID, "long_name", "longitude"),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, lonVarID, "cartesian_axis", "X"),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, lonVarID, "units", "degrees_east"),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, lonVarID, "valid_range", (/ 0.0_r8, 360.0_r8 /)),'nc_write_model_atts')

call nc_check(nf90_def_var(ncFileID, name="lat", xtype=nf90_double, &
                                                dimids=latDimID, varid=latVarID),'nc_write_model_atts' ) 
call nc_check(nf90_put_att(ncFileID, latVarID, "long_name", "latitude"),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, latVarID, "cartesian_axis", "Y"),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, latVarID, "units", "degrees_north"),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, latVarID, "valid_range", (/ -90.0_r8, 90.0_r8 /)),'nc_write_model_atts')

call nc_check(nf90_def_var(ncFileID, name="lev", xtype=nf90_double, &
                                                dimids=levDimID, varid=levVarID),'nc_write_model_atts') 
call nc_check(nf90_put_att(ncFileID, levVarID, "long_name", "level"),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, levVarID, "cartesian_axis", "Z"),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, levVarID, "units", "hPa"),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, levVarID, "positive", "down"),'nc_write_model_atts')

!----------------------------------------------------------------------------
! Create attributes for the state vector
!----------------------------------------------------------------------------
! TJH NOTE: the 'other' 3D models have netcdf variables that
! are allocated  
!      dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
! which is consistent with the Fortran convention of having the 
! fastest-moving index on the left.

call nc_check(nf90_def_var(ncid=ncFileID, name="u1", xtype=nf90_real, &
       dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
       varid  = u1VarID),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, u1VarID, "long_name", "zonal wind component"),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, u1VarID, "units", "m/s"),'nc_write_model_atts')
                                                                                                           
call nc_check(nf90_def_var(ncid=ncFileID, name="v1", xtype=nf90_real, &
       dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
       varid  = v1VarID),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, v1VarID, "long_name", "meridional wind component"),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, v1VarID, "units", "m/s"),'nc_write_model_atts')

call nc_check(nf90_def_var(ncid=ncFileID, name="t1", xtype=nf90_real, &
       dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
       varid  = t1VarID),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, t1VarID, "long_name", "temperature"),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, t1VarID, "units", "degrees Kelvin"),'nc_write_model_atts')


call nc_check(nf90_def_var(ncid=ncFileID, name="u", xtype=nf90_real, &
       dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
       varid  = uVarID),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, uVarID, "long_name", "zonal wind component"),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, uVarID, "units", "m/s"),'nc_write_model_atts')
                                                                                                           
call nc_check(nf90_def_var(ncid=ncFileID, name="v", xtype=nf90_real, &
       dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
       varid  = vVarID),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, vVarID, "long_name", "meridional wind component"),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, vVarID, "units", "m/s"),'nc_write_model_atts')

call nc_check(nf90_def_var(ncid=ncFileID, name="t", xtype=nf90_real, &
       dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
       varid  = tVarID),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, tVarID, "long_name", "temperature"),'nc_write_model_atts')
call nc_check(nf90_put_att(ncFileID, tVarID, "units", "degrees Kelvin"),'nc_write_model_atts')

! call nc_check(nf90_def_var(ncid=ncFileID, name="qnH", xtype=nf90_real, &
!        dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
!        varid  = qnHVarID),'nc_write_model_atts')
! call nc_check(nf90_put_att(ncFileID, qnHVarID, "long_name", "mixing ratio H"),'nc_write_model_atts')
! call nc_check(nf90_put_att(ncFileID, qnHVarID, "units", "?"),'nc_write_model_atts')

! call nc_check(nf90_def_var(ncid=ncFileID, name="qnOH", xtype=nf90_real, &
!        dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
!        varid  = qnOHVarID),'nc_write_model_atts')
! call nc_check(nf90_put_att(ncFileID, qnOHVarID, "long_name", "mixing ratio OH"),'nc_write_model_atts')
! call nc_check(nf90_put_att(ncFileID, qnOHVarID, "units", "?"),'nc_write_model_atts')

! call nc_check(nf90_def_var(ncid=ncFileID, name="qnO", xtype=nf90_real, &
!        dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
!        varid  = qnOVarID),'nc_write_model_atts')
! call nc_check(nf90_put_att(ncFileID, qnOVarID, "long_name", "mixing ratio O"),'nc_write_model_atts')
! call nc_check(nf90_put_att(ncFileID, qnOVarID, "units", "?"),'nc_write_model_atts')

call nc_check(nf90_enddef(ncfileID),'nc_write_model_atts')

!-------------------------------------------------------------------------------
! Fill the variables
!-------------------------------------------------------------------------------

call nc_check(nf90_put_var(ncFileID, lonVarID, lons),'nc_write_model_atts')
call nc_check(nf90_put_var(ncFileID, latVarID, lats),'nc_write_model_atts')
call nc_check(nf90_put_var(ncFileID, levVarID, levs),'nc_write_model_atts')

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID),'nc_write_model_atts')
write (*,*)'nc_write_model_atts: netCDF file ',ncFileID,' is synched ...'

end function nc_write_model_atts


function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH 23 May 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
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

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: u1VarID, v1VarID, t1VarID, uVarID, vVarID, tVarID
! integer :: qnHVarID, qnOHVarID, qnOVarID

type(model_type)                   :: Var

if ( .not. module_initialized ) call static_init_model

ierr = 0  ! assume normal termination

call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID),'nc_write_model_vars','inquire')


call init_model_instance(Var)
call vector_to_prog_var(statevec, Var)

call nc_check(NF90_inq_varid(ncFileID,  "u1",  u1VarID),'nc_write_model_vars','inq_varid u1')
call nc_check(nf90_put_var( ncFileID,  u1VarId, var%vars_3d(:,:,:, 1), &
           start=(/ 1, 1, 1, copyindex, timeindex /) ),'nc_write_model_vars','put_var u1')
                                                                                                                                     
call nc_check(NF90_inq_varid(ncFileID,  "v1",  v1VarID),'nc_write_model_vars','inq_varid v1')
call nc_check(nf90_put_var( ncFileID,  v1VarId, var%vars_3d(:,:,:, 2), &
           start=(/ 1, 1, 1, copyindex, timeindex /) ),'nc_write_model_vars','put_var v1')
                                                                                                                                     
call nc_check(NF90_inq_varid(ncFileID,  "t1",  t1VarID),'nc_write_model_vars','inq_varid t1')
call nc_check(nf90_put_var( ncFileID,  t1VarId, var%vars_3d(:,:,:, 3), &
           start=(/ 1, 1, 1, copyindex, timeindex /) ),'nc_write_model_vars','put_var t1')

call nc_check(NF90_inq_varid(ncFileID,  "u",  uVarID),'nc_write_model_vars','inq_ varid u')
call nc_check(nf90_put_var( ncFileID,  uVarId, var%vars_3d(:,:,:, 4), &
           start=(/ 1, 1, 1, copyindex, timeindex /) ),'nc_write_model_vars','put_var u')

call nc_check(NF90_inq_varid(ncFileID,  "v",  vVarID),'nc_write_model_vars','inq_varid v')
call nc_check(nf90_put_var( ncFileID,  vVarId, var%vars_3d(:,:,:, 5), &
           start=(/ 1, 1, 1, copyindex, timeindex /) ),'nc_write_model_vars','put_var v')

call nc_check(NF90_inq_varid(ncFileID,  "t",  tVarID),'nc_write_model_vars','inq_varid t')
call nc_check(nf90_put_var( ncFileID,  tVarId, var%vars_3d(:,:,:, 6), &
           start=(/ 1, 1, 1, copyindex, timeindex /) ),'nc_write_model_vars','put_var t')

! call nc_check(NF90_inq_varid(ncFileID,  "qnH",  qnHVarID),'nc_write_model_vars','inq_varid qnH')
! call nc_check(nf90_put_var( ncFileID,  qnHVarId, var%vars_3d(:,:,:, 7), &
!            start=(/ 1, 1, 1, copyindex, timeindex /) ),'nc_write_model_vars','put_var qnH')

! call nc_check(NF90_inq_varid(ncFileID,  "qnOH",  qnOHVarID),'nc_write_model_vars','inq_varid qnOH')
! call nc_check(nf90_put_var( ncFileID,  qnOHVarId, var%vars_3d(:,:,:, 8), &
!            start=(/ 1, 1, 1, copyindex, timeindex /) ),'nc_write_model_vars','put_var qnOH')

! call nc_check(NF90_inq_varid(ncFileID,  "qnH",  qnOVarID),'nc_write_model_vars','inq_varid qnH')
! call nc_check(nf90_put_var( ncFileID,  qnOVarId, var%vars_3d(:,:,:, 9), &
!            start=(/ 1, 1, 1, copyindex, timeindex /) ),'nc_write_model_vars','put_var qnH')

if (do_output()) write (*,*)'Finished filling variables ...'
call nc_check(nf90_sync(ncFileID),'nc_write_model_vars','sync')
if (do_output()) write (*,*)'netCDF file is synched ...'

call end_model_instance(Var)   ! should avoid any memory leaking

end function nc_write_model_vars


subroutine pert_model_state(state, pert_state, interf_provided)
!------------------------------------------------------------------
!
! Perturbs a model state for generating initial ensembles.
! The perturbed state is returned in pert_state.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the inter_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding an O(0.1) magnitude perturbation to each
! model state variable independently. The interf_provided argument
! should be returned as .true. if the model want to do its own
! perturbing of states.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

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
   if(variable_type == QTY_U_WIND_COMPONENT .or. &
      variable_type == QTY_V_WIND_COMPONENT .or. &
      variable_type == QTY_TEMPERATURE ) then
      pert_state(i) = &
!       & random_gaussian(random_seq, state(i), state(i)*0.05)
        & state(i)
   else
      pert_state(i) = state(i)
   endif
end do

end subroutine pert_model_state


subroutine prog_var_to_vector(var, x)
!=======================================================================
! based on prog_var_to_vector in model_mod for CAM.
!

type(model_type), intent(in) :: var
real(r8), intent(out) :: x(:)

integer :: i, j, k, nf, indx

if ( .not. module_initialized ) call static_init_model

! Do order as ps, t, u, v, q, tracers to be consistent with b-grid

! Start copying fields to straight vector
! TJH NOTE: the 'other' 3D models have variables that
! are allocated  [lon,lat,lev ...] which is consistent with the Fortran
! convention of having the fastest-moving index on the left.
indx = 0
do i = 1, nx !longitude
   do j = 1, ny !latitude
      ! u,v,t,q, and tracers at successively lower levels
      do k = 1, nz ! height
         do nf= 1, state_num_3d
            indx = indx + 1
            x(indx) = var%vars_3d(k, i, j, nf)
         end do
      end do
   end do
end do

! Temporary check
if(indx /= model_size) then
   write(msgstring, *) 'indx ',indx,' model_size ',model_size,' must be equal '
   call error_handler(E_ERR, 'prog_var_to_vector', msgstring, source, revision, revdate)
endif

end subroutine prog_var_to_vector


subroutine vector_to_prog_var(x, var) 
!=======================================================================
! subroutine vector_to_prog_var(x, var) 
!

real(r8), intent(in) :: x(:)
type(model_type), intent(out) :: var

integer :: i, j, k, nf, indx

if ( .not. module_initialized ) call static_init_model

! Start copying fields from straight vector
! TJH NOTE: the 'other' 3D models have variables that
! are allocated  [lon,lat,lev ...] which is consistent with the Fortran
! convention of having the fastest-moving index on the left.
indx = 0
do i = 1, nx ! longitude
   do j = 1, ny ! latitude
      ! u,v,t,q  and tracers at successive levels
      do k = 1, nz ! height
         do nf = 1, state_num_3d
            indx = indx + 1
            var%vars_3d(k, i, j, nf) = x(indx)
         end do 
      end do
   end do
end do

! Temporary check
if(indx /= model_size) then
   write(msgstring, *) 'indx ',indx,' model_size ',model_size,' must be equal '
   call error_handler(E_ERR, 'vector_to_prog_var', msgstring, source, revision, revdate)
endif

end subroutine vector_to_prog_var



subroutine init_model_instance(var, valid_time)
!=======================================================================
! subroutine init_model_instance(var)
!
! Initializes an instance of a ROSE model state variable

type(model_type),          intent(out) :: var
type(time_type), optional, intent( in) :: valid_time

if ( .not. module_initialized ) call static_init_model

! Initialize the storage space and return
! TJH NOTE: the 'other' 3D models have variables that
! are allocated  [lon,lat,lev ...] which is consistent with the Fortran
! convention of having the fastest-moving index on the left.

allocate(var%vars_3d(nz, nx, ny, state_num_3d))

if (present(valid_time)) then
   var%valid_time = valid_time
else
   var%valid_time = set_time_missing()
endif

end subroutine init_model_instance



subroutine end_model_instance(var)
!=======================================================================
! subroutine end_model_instance(var)
!
! Ends an instance of a rose model state variable
                                                                                                         
type(model_type), intent(inout) :: var
                                                                                                         
if ( .not. module_initialized ) call static_init_model

deallocate(var%vars_3d)
                                                                                                         
end subroutine end_model_instance



subroutine update_ROSE_restart(file_name, var)
!=======================================================================
! update ROSE restart file fields

  character (len = *), intent(in) :: file_name
  type(model_type),    intent(in) :: var

! ROSE restart file prognostic variables
  real, dimension (nz,nx,ny) :: un1, vn1, tn1, un0, vn0, tn0

! real, dimension (nz,nx,ny,nbcon) :: qn1
! real, dimension (nz,ny) :: q_o2, q_n2

  integer :: restart_id
  integer :: idd
  integer :: ncerr,  nlons, nlats, nlevs
  integer :: var_id 
  integer, dimension(3) :: mtime
  integer :: seconds, days

!====================================================================

if ( .not. module_initialized ) call static_init_model

if( .not. file_exist(file_name)) then
  write(msgstring,*) trim(adjustl(file_name)),' not available.'
  call error_handler(E_ERR,'update_ROSE_restart',msgstring,source,revision,revdate)
endif

if (do_output()) print *, 'update_ROSE_restart: reading restart'
ncerr = nf90_open( file_name, NF90_WRITE, restart_id )
call nc_check(ncerr, 'update_ROSE_restart','open')  ! will die if error
if (do_output()) print *, 'update_ROSE_restart: opened with '//trim(nf90_strerror(ncerr))

!... check for matching dimensions

  call nc_check( nf90_inq_dimid(restart_id, 'lon', idd), &
         'update_rose_restart', 'inq_dimid lon')
  call nc_check( nf90_inquire_dimension(restart_id, idd, len=nlons), &
         'update_rose_restart', 'inquire_dimension lon')

  call nc_check( nf90_inq_dimid(restart_id, 'lat', idd), &
         'update_rose_restart', 'inq_dimid lat')
  call nc_check( nf90_inquire_dimension(restart_id, idd, len=nlats), &
         'update_rose_restart', 'inquire_dimension lat')

  call nc_check( nf90_inq_dimid(restart_id, 'lev', idd), &
         'update_rose_restart', 'inq_dimid lev')
  call nc_check( nf90_inquire_dimension(restart_id, idd, len=nlevs), &
         'update_rose_restart', 'inquire_dimension lev')

  if ( nlons /= nx ) then
   write(msgstring,*) trim(file_name),' nlons = ',nlons,' DART expects ',nx
   call error_handler(E_ERR,'update_rose_restart',msgstring,source,revision,revdate)
  endif

  if ( nlats /= ny ) then
   write(msgstring,*) trim(file_name),' nlats = ',nlats,' DART expects ',ny
   call error_handler(E_ERR,'update_rose_restart',msgstring,source,revision,revdate)
  endif

  if ( nlevs /= nz ) then
   write(msgstring,*) trim(file_name),' nlevs = ',nlevs,' DART expects ',nz
   call error_handler(E_ERR,'update_rose_restart',msgstring,source,revision,revdate)
  endif

! unpack the variables into something familiar and then stuff them
! into the existing netCDF variables. 

  call get_time(var%valid_time, seconds, days)

  un1 = var%vars_3d(:,:,:, 1)
  vn1 = var%vars_3d(:,:,:, 2)
  tn1 = var%vars_3d(:,:,:, 3)
  un0 = var%vars_3d(:,:,:, 4)
  vn0 = var%vars_3d(:,:,:, 5)
  tn0 = var%vars_3d(:,:,:, 6)
  !qn1(:,:,:,7)  = var%vars_3d(:,:,:,7) ! H
  !qn1(:,:,:,8)  = var%vars_3d(:,:,:,8) ! OH
  !qn1(:,:,:,18) = var%vars_3d(:,:,:,9) ! O

  call nc_check(  nf90_inq_varid(restart_id, 'un1', var_id), &
         'update_rose_restart', 'inq_varid un1')
  call nc_check(  nf90_put_var(restart_id, var_id, values=un1), &
         'update_rose_restart', 'put_var un1')

  call nc_check(  nf90_inq_varid(restart_id, 'vn1', var_id), &
         'update_rose_restart', 'inq_varid vn1')
  call nc_check(  nf90_put_var(restart_id, var_id, values=vn1), &
         'update_rose_restart', 'put_var vn1')

  call nc_check(  nf90_inq_varid(restart_id, 'tn1', var_id), &
         'update_rose_restart', 'inq_varid tn1')
  call nc_check(  nf90_put_var(restart_id, var_id, values=tn1), &
         'update_rose_restart', 'put_var tn1')

  call nc_check(  nf90_inq_varid(restart_id, 'un0', var_id), &
         'update_rose_restart', 'inq_varid un0')
  call nc_check(  nf90_put_var(restart_id, var_id, values=un0), &
         'update_rose_restart', 'put_var un0')

  call nc_check(  nf90_inq_varid(restart_id, 'vn0', var_id), &
         'update_rose_restart', 'inq_varid vn0')
  call nc_check(  nf90_put_var(restart_id, var_id, values=vn0), &
         'update_rose_restart', 'put_var vn0')

  call nc_check(  nf90_inq_varid(restart_id, 'tn0', var_id), &
         'update_rose_restart', 'inq_varid tn0')
  call nc_check(  nf90_put_var(restart_id, var_id, values=tn0), &
         'update_rose_restart', 'put_var tn0')

  ! FIXME - since we're not handling the year correctly in rose,
  ! we need to capture the existing one and reuse it.

  call nc_check(  nf90_inq_varid(restart_id, 'mtime', var_id), &
            'update_rose_restart','inq_varid mtime')
  call nc_check(nf90_get_var( restart_id, var_id, mtime) , &
                     'update_rose_restart','get_var mtime')

!  mtime(1) = cal_year  FIXME ... this is that I'm talking about
   mtime(2) = days
   mtime(3) = seconds

  call nc_check(nf90_put_var( restart_id, var_id, mtime) , &
                     'update_rose_restart','get_var mtime')

if (do_output()) print *, 'update_ROSE_restart: mtime (year/doy/seconds):', mtime

  call nc_check( nf90_sync( restart_id), 'update_rose_restart', 'sync')
  call nc_check( nf90_close(restart_id), 'update_rose_restart', 'close')

end subroutine update_ROSE_restart


subroutine read_ROSE_restart(file_name, var, model_time)
!=======================================================================
! read ROSE restart file fields that have been updated
!
  character (len = *), intent(in) :: file_name
  type(model_type),   intent(out) :: var
  type(time_type),    intent(out) :: model_time

  integer :: ncerr, idd
  integer :: var_id, restart_id, nlons, nlats, nlevs
  integer, dimension(3) :: mtime

! ROSE restart file prognostic variables
  integer :: doy                   ! day of year
  integer :: utsec
  integer :: cal_year              ! calendar year
  real, dimension (nz,nx,ny) :: un1, vn1, tn1, un0, vn0, tn0
! real, dimension (nz)       :: tref                                 

! from chem.mod.f ... constituent mixing ratios
! real, dimension (nz,nx,ny,nbcon) :: qn1
! real, dimension (nz,ny)          :: q_o2, q_n2
!
!====================================================================

if ( .not. module_initialized ) call static_init_model

if( .not. file_exist(file_name)) then
   write(msgstring,*)trim(file_name)//' does not exist.'
   call error_handler(E_ERR,'read_ROSE_restart',msgstring,source,revision,revdate)
endif

if (do_output()) print *, 'read_ROSE_restart:reading restart:', file_name
ncerr = nf90_open( file_name, NF90_NOWRITE, restart_id )
call nc_check(ncerr, 'read_ROSE_restart', 'open')
if (do_output()) print *, 'read_ROSE_restart:opened with '//trim(nf90_strerror(ncerr))

!... check for matching dimensions

!  call dim_check( 'lon', nx)
!  call dim_check( 'lat', ny)
!  call dim_check( 'lev', nz)

!... check for matching dimensions

   call nc_check( nf90_inq_dimid(restart_id, 'lon', idd), &
          'read_rose_restart', 'inq_dimid lon')
   call nc_check( nf90_inquire_dimension(restart_id, idd, len=nlons), &
          'read_rose_restart', 'inquire_dimension lon')

   call nc_check( nf90_inq_dimid(restart_id, 'lat', idd), &
          'read_rose_restart', 'inq_dimid lat')
   call nc_check( nf90_inquire_dimension(restart_id, idd, len=nlats), &
          'read_rose_restart', 'inquire_dimension lat')

   call nc_check( nf90_inq_dimid(restart_id, 'lev', idd), &
          'read_rose_restart', 'inq_dimid lev')
   call nc_check( nf90_inquire_dimension(restart_id, idd, len=nlevs), &
          'read_rose_restart', 'inquire_dimension lev')

   if ( nlons /= nx ) then
    write(msgstring,*) trim(file_name),' nlons = ',nlons,' DART expects ',nx
    call error_handler(E_ERR,'read_rose_restart',msgstring,source,revision,revdate)
   endif

   if ( nlats /= ny ) then
    write(msgstring,*) trim(file_name),' nlats = ',nlats,' DART expects ',ny
    call error_handler(E_ERR,'read_rose_restart',msgstring,source,revision,revdate)
   endif

   if ( nlevs /= nz ) then
    write(msgstring,*) trim(file_name),' nlevs = ',nlevs,' DART expects ',nz
    call error_handler(E_ERR,'read_rose_restart',msgstring,source,revision,revdate)
   endif

!... get dynamical variables (un1, vn1, tn1, un0, vn0, tn0)

!  call var_read( 'un1', un1)
!  call var_read( 'vn1', vn1)
!  call var_read( 'tn1', tn1)
!  call var_read( 'un0', un0)
!  call var_read( 'vn0', vn0)
!  call var_read( 'tn0', tn0)

   call nc_check(  nf90_inq_varid(restart_id, 'un1', var_id), &
          'read_rose_restart', 'inq_varid un1')
   call nc_check(  nf90_get_var(restart_id, var_id, values=un1), &
          'read_rose_restart', 'get_var un1')

   call nc_check(  nf90_inq_varid(restart_id, 'vn1', var_id), &
          'read_rose_restart', 'inq_varid vn1')
   call nc_check(  nf90_get_var(restart_id, var_id, values=vn1), &
          'read_rose_restart', 'get_var vn1')

   call nc_check(  nf90_inq_varid(restart_id, 'tn1', var_id), &
          'read_rose_restart', 'inq_varid tn1')
   call nc_check(  nf90_get_var(restart_id, var_id, values=tn1), &
          'read_rose_restart', 'get_var tn1')

   call nc_check(  nf90_inq_varid(restart_id, 'un0', var_id), &
          'read_rose_restart', 'inq_varid un0')
   call nc_check(  nf90_get_var(restart_id, var_id, values=un0), &
          'read_rose_restart', 'get_var un0')

   call nc_check(  nf90_inq_varid(restart_id, 'vn0', var_id), &
          'read_rose_restart', 'inq_varid vn0')
   call nc_check(  nf90_get_var(restart_id, var_id, values=vn0), &
          'read_rose_restart', 'get_var vn0')

   call nc_check(  nf90_inq_varid(restart_id, 'tn0', var_id), &
          'read_rose_restart', 'inq_varid tn0')
   call nc_check(  nf90_get_var(restart_id, var_id, values=tn0), &
          'read_rose_restart', 'get_var tn0')

   call nc_check(nf90_inq_varid( restart_id, 'mtime', var_id), &
             'read_rose_restart','inq_varid mtime')
   call nc_check(nf90_get_var( restart_id, var_id, mtime) , &
                      'read_rose_restart','get_var mtime')

   call nc_check(nf90_inq_varid( restart_id, 'tref', var_id), &
             'read_rose_restart','inq_varid tref')
   call nc_check(nf90_get_var(restart_id, var_id, tref) , &
                      'read_rose_restart','get_var tref')

   ! That's it ...

   call nc_check(nf90_close( restart_id),'read_rose_restart','close')

   if (do_output()) print *, 'read_ROSE_restart: mtime (year/doy/seconds):', mtime
   cal_year = mtime(1)
   doy      = mtime(2)
   utsec    = mtime(3)

   ! repack the variables into our structure

   var%vars_3d(:,:,:,1) = un1 
   var%vars_3d(:,:,:,2) = vn1 
   var%vars_3d(:,:,:,3) = tn1 
   var%vars_3d(:,:,:,4) = un0           ! at model_time
   var%vars_3d(:,:,:,5) = vn0           ! at model_time
   var%vars_3d(:,:,:,6) = tn0           ! at model_time
   !var%vars_3d(:,:,:,7) = qn1(:,:,:,7)  ! H
   !var%vars_3d(:,:,:,8) = qn1(:,:,:,8)  ! OH
   !var%vars_3d(:,:,:,9) = qn1(:,:,:,18) ! O

   ! FIXME ... this ignores years - no calendar
   model_time = set_time(utsec, doy)

   var%valid_time = model_time

   if (do_output()) call print_time(model_time, str=" read_ROSE_restart: model_time ")

end subroutine read_ROSE_restart


subroutine update_ROSE_namelist(file_name, time1, timeN, ens_member, &
       atune, ptune )
!=======================================================================
! Update the ROSE namelist - especially the new target_time.
! The target_time is actually an offset - so we need to calculate that.
!=======================================================================

character(len=*),   intent(in) :: file_name
type(time_type),    intent(in) :: time1, timeN
integer,            intent(in) :: ens_member
real(r8), optional, intent(in) :: atune, ptune

type(time_type) :: forecast_length
integer :: second, day
integer :: iunit

forecast_length = timeN - time1
call get_time(forecast_length, second, day)
target_time = real(   day,digits12)*24.0_digits12 + &
              real(second,digits12)/3600.0_digits12

if (do_output()) then
   PRINT*,'update_ROSE_namelist: forecast length [days seconds] ', &
                 day, second, ' = hours ', target_time
endif

! Update the other namelist parameters 
ens_element     = ens_member
if (present(atune)) amp_tune = atune 
if (present(ptune)) pha_tune = ptune 

iunit = open_file(file_name, action='write')
write(iunit, nml=ROSE_NML)
call close_file(iunit)

end subroutine update_ROSE_namelist



subroutine read_ROSE_tref(file_name)
!=======================================================================
! read the reference temperature from a ROSE restart file
!====================================================================

   character (len = *), intent(in) :: file_name

   integer :: dim_id, dim_len
   integer :: var_id, restart_id

   if( .not. file_exist(file_name)) then
     write(msgstring,*) trim(adjustl(file_name)),' not available.'
     call error_handler(E_ERR,'read_ROSE_tref',msgstring,source,revision,revdate)
   endif

   call nc_check(nf90_open(file_name, NF90_NOWRITE, restart_id), 'read_ROSE_tref', 'open')

   call nc_check(nf90_inq_dimid(restart_id, 'lev', dim_id), 'read_ROSE_tref', 'inq_dimid lev')
   call nc_check(nf90_inquire_dimension(restart_id, dim_id, len=dim_len),'read_ROSE_tref','inquire_dimension lev')

   if (dim_len .ne. nz) then
     write(msgstring,'(a,2i5)') 'restart dimension mismatch', dim_len, nz
     call error_handler(E_ERR,'read_ROSE_tref',msgstring,source,revision,revdate)
   endif

   call nc_check(nf90_inq_varid(restart_id, 'tref', var_id),'read_ROSE_tref', 'inq_varid tref')
   call nc_check(nf90_get_var(restart_id, var_id, tref),'read_ROSE_tref', 'get_var tref')

   call nc_check(nf90_close(restart_id),'read_ROSE_tref', 'close')

end subroutine read_ROSE_tref



subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! Not used ????

real(r8), intent(in) :: ens_mean(:)

if ( .not. module_initialized ) call static_init_model

end subroutine ens_mean_for_model


!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
