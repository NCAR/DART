! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> This is the interface between the OpenGGCM space weather model and DART.

module model_mod

! Modules that are absolutely required for use are listed
use        types_mod,    only : r4, r8, i4, i8, SECPERDAY, MISSING_R8, rad2deg, PI, &
                                earth_radius

use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=), GREGORIAN

use     location_mod, only : location_type, get_dist, get_close_maxdist_init,  &
                             get_close_obs_init, set_location,                 &
                             VERTISUNDEF, VERTISHEIGHT, get_location,          &
                             vert_is_height, vert_is_level, vert_is_surface,   &
                             vert_is_undef, get_close_type,                    &
                             loc_get_close_obs => get_close_obs, write_location

use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper,                    &
                             find_namelist_in_file, check_namelist_read,       &
                             file_exist, find_textfile_dims, file_to_text,     &
                             do_nml_file, do_nml_term, nmlfileunit, open_file, &
                             close_file

use     obs_kind_mod, only : KIND_ELECTRON_DENSITY, KIND_ELECTRIC_POTENTIAL,   &
                             get_raw_obs_kind_index, get_raw_obs_kind_name,    &
                             paramname_length 

use     mpi_utilities_mod, only : my_task_id, task_count

use        random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

use ensemble_manager_mod,  only : ensemble_type, map_pe_to_task, get_copy_owner_index, &
                                  get_var_owner_index

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain, get_model_variable_indices,        &
                                  get_num_variables, get_index_start,            &
                                  get_num_dims, get_domain_size, get_kind_index, &
                                  get_varid_from_kind, get_dart_vector_index,    &
                                  get_dim_name, get_variable_name,               &
                                  state_structure_info

use dart_time_io_mod,      only : write_model_time

use cotr_mod,              only : transform, cotr_set, cotr, xyzdeg, degxyz

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
          get_close_obs,                 &
          query_vert_localization_coord, &
          vert_convert,                  &
          construct_file_name_in,        &
          read_model_time,               &
          write_model_time



! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! Return netcdf data array given a variable name
interface get_data
   module procedure get_data_1d
   module procedure get_data_2d
   module procedure get_data_3d
end interface get_data

! Transform type openggcm, this come directly from the 
! openggmc source code that had been modified and put
! as a subroutine in cotr_mod.nml
type(transform) :: openggcm_transform

! Message strings
character(len=512) :: msgstring
character(len=512) :: msgstring2
character(len=512) :: msgstring3

integer, parameter :: VERT_LEVEL_1 = 1

! Grid types
integer, parameter :: GEOGRAPHIC_GRID = 1
integer, parameter :: MAGNETIC_GRID   = 2

! Grid convertion direction
integer, parameter :: SM_TO_GEO = 1
integer, parameter :: GEO_TO_SM = 2

! Logical to keep track of if we have initialized static_init_model
logical, save :: module_initialized = .false.

! Storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq

! DART state vector contents are specified in the input.nml:&model_nml namelist.
integer, parameter :: max_state_variables = 10 
integer, parameter :: num_state_table_columns = 4
! NOTE: may need to increase character length if netcdf variables are
! larger than paramname_length = 32.
character(len=paramname_length) :: variable_table( max_state_variables, num_state_table_columns )
integer :: state_kinds_list( max_state_variables )
logical ::  update_var_list( max_state_variables )
integer ::   grid_info_list( max_state_variables )
integer ::   dim_order_list( max_state_variables, 3 )

! identifiers for variable_table
integer, parameter :: VAR_NAME_INDEX   = 1
integer, parameter :: VAR_KIND_INDEX   = 2
integer, parameter :: VAR_UPDATE_INDEX = 3

! identifiers for LAT, LON and HEIGHT
integer, parameter :: VAR_LAT_INDEX   = 1
integer, parameter :: VAR_LON_INDEX   = 2
integer, parameter :: VAR_HGT_INDEX   = 3

! things which can/should be in the model_nml
character(len=NF90_MAX_NAME) :: openggcm_template
integer  :: assimilation_period_days     = 1
integer  :: assimilation_period_seconds  = 0
real(r8) :: model_perturbation_amplitude = 0.2
character(len=paramname_length) :: model_state_variables(max_state_variables * num_state_table_columns ) = ' '
integer  :: debug = 0   ! turn up for more and more debug messages

namelist /model_nml/  &
   openggcm_template,           &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   model_state_variables,       &
   debug

! Generic grid type to hold CTIM and Magnetic grid information
type grid_type
   ! the size of the grid
   integer :: nlon, nlat, nheight

   ! grid information
   real(r8), allocatable :: longitude(:), latitude(:)
   logical :: uses_colatitude
  
   ! decide if we can get away with 1D heights or if they
   ! have to be per-column.  would prefer 1D for simplicity
   ! if possible.
   !real(r8), allocatable :: heights(:,:,:)
   real(r8), allocatable :: heights(:)

   ! optional conversion grid - 2d arrays: (lon, lat)
   real(r8), allocatable :: conv_2d_lon(:,:), conv_2d_lat(:,:)
end type grid_type


!------------------------------------------------------------------
! Global Variables 
!------------------------------------------------------------------

! Number of fields in the state vector
integer :: nfields

! Geometric (CTIM) Grid and Magnetic Grid
type(grid_type), target :: geo_grid, mag_grid

! Global Time Variables
type(time_type) :: model_time, model_timestep

! The state vector length
integer(i8) :: model_size

! Domain id to be used by routines in state_structure_mod
integer :: domain_id

contains

!------------------------------------------------------------------
!------------------------------------------------------------------

!> Called to do one time initialization of the model. In this case,
!> it reads in the grid information.

subroutine static_init_model()

integer :: iunit, io, ncid
integer :: ss, dd


if ( module_initialized ) return ! only need to do this once.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! set calendar type
call set_calendar_type(GREGORIAN)

! Set the time step ... causes openggcm namelists to be read.
! Ensures model_timestep is multiple of 'ocean_dynamics_timestep'

model_timestep = set_time(assimilation_period_days,assimilation_period_seconds)

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(msgstring,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',msgstring,source,revision,revdate)

! get data dimensions, then allocate space, then open the files
! and actually fill in the arrays.

ncid = get_grid_template_fileid(openggcm_template)

call get_grid_sizes(ncid, geo_grid, 'cg_lon', 'cg_lat','cg_height')
call get_grid_sizes(ncid, mag_grid, 'ig_lon', 'ig_lat')

! allocate space for geographic and magnetic grids

call allocate_grid_space(geo_grid, conv=.false.)
call allocate_grid_space(mag_grid, conv=.true.)

! read in geographic and magnetic grids, and for the mag grid read
! in the 2d conversion arrays to go to geographic coords

call read_horiz_grid(ncid, geo_grid, 'cg_lon',  'cg_lat',  is_conv=.false., is_co_latitude=.false.)
call read_horiz_grid(ncid, mag_grid, 'ig_lon',  'ig_lat',  is_conv=.false., is_co_latitude=.true.) 
call read_horiz_grid(ncid, mag_grid, 'geo_lon', 'geo_lat', is_conv=.true.,  is_co_latitude=.true.)

call read_vert_levels(ncid,geo_grid,'cg_height')
call read_vert_levels(ncid,mag_grid,'ig_height')

! verify that the model_state_variables namelist was filled in correctly.  
! returns variable_table which has variable names, kinds and update strings, 
! and grid information.

call verify_state_variables(model_state_variables, nfields, variable_table, &
                            state_kinds_list, update_var_list, grid_info_list)

! fill up the state structure with information from the model template file

domain_id = add_domain(openggcm_template, nfields, &
                       var_names   = variable_table  (1:nfields , VAR_NAME_INDEX), &
                       kind_list   = state_kinds_list(1:nfields), &
                       update_list = update_var_list (1:nfields))

if (debug > 0) call state_structure_info(domain_id)

! order the dimensions according to lat, lon and height

call make_dim_order_table(nfields)

! only one domain in the model at the moment.

model_size = get_domain_size(domain_id)
write(msgstring,*)'model_size = ', model_size
call error_handler(E_MSG,'model_interpolate',msgstring,source,revision,revdate)

! set the transform geo -> magnetic grid
call initialize_openggcm_transform(openggcm_template)

end subroutine static_init_model

!------------------------------------------------------------------

!> Returns the size of the model state vector as an I8 integer.
!> The size is computed in static_init_model() an stored in a
!> module global variable.

function get_model_size()

integer(i8) :: get_model_size  !< state vector length

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size

!------------------------------------------------------------------

!> Model interpolate will interpolate any state variable
!> the given location given a state vector. The 'generic kind' of the variable being
!> interpolated is obs_kind since normally this is used to find the expected
!> value of an observation at some location. The interpolated value is 
!> returned in interp_val and istatus is 0 for success.

subroutine model_interpolate(state_handle, ens_size, location, obs_kind, expected_obs, istatus)

type(ensemble_type), intent(in) :: state_handle !< ensemble handle for data to interpolate in
integer,             intent(in) :: ens_size !< number of ensembles, sets size of expected_obs and istatus arrays
type(location_type), intent(in) :: location !< dart location to interpolate to
integer,             intent(in) :: obs_kind !< dart kind to interpolate
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size) !< array of returned statuses


! Local storage
real(r8)    :: loc_array(3), llon, llat, lheight
integer     :: ind
integer     :: hgt_bot, hgt_top
real(r8)    :: hgt_fract
integer     :: hstatus
type(grid_type), pointer :: mygrid

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the 
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

expected_obs(:) = MISSING_R8     ! the DART bad value flag
istatus(:)      = 99             ! unknown error

! Get the individual locations values
loc_array = get_location(location)
llon    = loc_array(1)
llat    = loc_array(2)
lheight = loc_array(3)

if( vert_is_undef(location) ) then
   ! this is what we expect and it is ok
elseif ( vert_is_height(location) ) then
   ! this is what we expect and it is ok
   ! once we write the code to search in the vertical
elseif (vert_is_level(location)) then
   write(msgstring,*)'requesting interp of an obs on level, not supported yet'
   call error_handler(E_ERR,'model_interpolate',msgstring,source,revision,revdate)

   !> @todo FIXME something like this
   ! convert the heights index to an actual height 
   ind = nint(loc_array(3))
   if ( (ind < 1) .or. (ind > size(geo_grid%heights,1)) ) then 
      istatus = 11
      return
   else
      !>@todo FIXME : assuming everything is flat at the moment
      lheight = geo_grid%heights(1)
   endif
else   ! if pressure or surface we don't know what to do
   write(msgstring,*)'requesting interp of an obs on pressure or surface, not supported yet'
   call error_handler(E_ERR,'model_interpolate',msgstring,source,revision,revdate)

   istatus = 17
   return
endif

! Determine which grid the incoming observation is on. If the observation
! is on the magnetic grid transform it to the geographic grid.

SELECT CASE (get_grid_type(obs_kind))
   CASE (MAGNETIC_GRID)
      call transform_mag_geo(llon, llat, lheight, GEO_TO_SM)
      mygrid => mag_grid

   CASE (GEOGRAPHIC_GRID)
      mygrid => geo_grid

   CASE DEFAULT
      call error_handler(E_ERR, 'model_interpolate', 'unknown grid type, should not happen', &
            source, revision, revdate)
END SELECT


if( vert_is_undef(location) ) then
   call lon_lat_interpolate(state_handle, ens_size, mygrid, obs_kind, llon, llat, VERT_LEVEL_1, &
                            expected_obs, istatus)

   return
endif

call height_bounds(lheight, mygrid%nheight, mygrid%heights, hgt_bot, hgt_top, hgt_fract, hstatus)

if(hstatus /= 0) then
   istatus = 12
   return
endif

! do a 2d interpolation for the value at the bottom level, then again for
! the top level, then do a linear interpolation in the vertical to get the
! final value.  this sets both interp_val and istatus.
call do_interp(state_handle, ens_size, mygrid, hgt_bot, hgt_top, hgt_fract, &
               llon, llat, obs_kind, expected_obs, istatus)

end subroutine model_interpolate

!------------------------------------------------------------------

!> Subroutine to interpolate to a lon lat location given the state handle.
!> Successful interpolation returns istatus=0.

subroutine lon_lat_interpolate(state_handle, ens_size, grid_handle, var_kind, &
                               lon, lat, height_index, expected_obs, istatus)

type(ensemble_type), intent(in)  :: state_handle !< state ensemble handle
integer,             intent(in)  :: ens_size !< ensemble size
type(grid_type),     intent(in)  :: grid_handle !< geo or mag grid
integer,             intent(in)  :: var_kind !< dart variable kind
real(r8),            intent(in)  :: lon !< longitude to interpolate
real(r8),            intent(in)  :: lat !< latitude to interpolate
integer,             intent(in)  :: height_index !< height index to interpolate
real(r8),            intent(out) :: expected_obs(ens_size) !< returned interpolations
integer,             intent(out) :: istatus(ens_size) !< returned statuses


! Local storage, 
integer  :: lat_bot, lat_top, lon_bot, lon_top
real(r8) :: p(4,ens_size), xbot(ens_size), xtop(ens_size)
real(r8) :: lon_fract, lat_fract

! Succesful return has istatus of 0
istatus = 0

! find the lower and upper indices which enclose the given value
! in this model, the data at lon 0 is replicated at lon 360, so no special
! wrap case is needed.
call lon_bounds(lon, grid_handle, lon_bot, lon_top, lon_fract)

if (grid_handle%uses_colatitude) then
   call colat_bounds(lat, grid_handle, lat_bot, lat_top, lat_fract, istatus(1))
else
   call lat_bounds(lat, grid_handle, lat_bot, lat_top, lat_fract, istatus(1))
endif

if (istatus(1) /= 0) then
   istatus(:) = 18 
   return
endif

! Get the values at the four corners of the box or quad
! Corners go around counterclockwise from lower left
p(1, :) = get_val(lon_bot, lat_bot, height_index, var_kind, state_handle, ens_size)
p(2, :) = get_val(lon_top, lat_bot, height_index, var_kind, state_handle, ens_size)
p(3, :) = get_val(lon_top, lat_top, height_index, var_kind, state_handle, ens_size)
p(4, :) = get_val(lon_bot, lat_top, height_index, var_kind, state_handle, ens_size)

! Rectangular bilinear interpolation
xbot = p(1, :) + lon_fract * (p(2, :) - p(1, :))
xtop = p(4, :) + lon_fract * (p(3, :) - p(4, :))

! Now interpolate in latitude
expected_obs = xbot + lat_fract * (xtop - xbot)

end subroutine lon_lat_interpolate

!------------------------------------------------------------

 !> Returns the value for a single model level given the lat and lon indices

function get_val(lon_index, lat_index, height_index, var_kind, state_handle, ens_size)

integer,             intent(in)  :: lon_index !< longitude index
integer,             intent(in)  :: lat_index !< lattude index
integer,             intent(in)  :: height_index !< height index
integer,             intent(in)  :: var_kind !< dart variable kind
type(ensemble_type), intent(in)  :: state_handle !< ensemble handle for state vector
integer,             intent(in)  :: ens_size !< size of the ensemble


! Local variables
real(r8)    :: get_val(ens_size)
integer(i8) :: state_index
integer     :: var_id

integer, dimension(3) :: dim_index

if (var_kind < 0 ) then
   write(msgstring,*) 'dart kind < 0 which should not happen'
   write(msgstring2,*) 'the dart kind provided is : ', var_kind
   call error_handler(E_ERR, 'get_val', msgstring, &
                      source, revision, revdate, text2=msgstring2)
endif   

var_id = get_varid_from_kind(domain_id, var_kind)

dim_index(dim_order_list(var_id, VAR_LON_INDEX)) = lon_index
dim_index(dim_order_list(var_id, VAR_LAT_INDEX)) = lat_index
dim_index(dim_order_list(var_id, VAR_HGT_INDEX)) = height_index

state_index = get_dart_vector_index(dim_index(1), dim_index(2), dim_index(3), &
                                    domain_id, var_id)

get_val = get_state(state_index, state_handle)

end function get_val

!------------------------------------------------------------

!> Given a longitude lon and a grid handle which contains both the 1D array 
!> of longitudes and the grid longitude size, returns the indices of the grid
!> below and above the location longitude and the fraction of the distance
!> between.  This code assumes that the first and last rows are replicated
!> and identical (e.g. 0 and 360 both have entries in the array)

subroutine lon_bounds(lon, grid_handle, bot, top, fract)

real(r8),        intent(in)  :: lon !< input longitude
type(grid_type), intent(in)  :: grid_handle  !< geo or mag grid
integer,         intent(out) :: bot !< index of bottom layer
integer,         intent(out) :: top !< index of top layer
real(r8),        intent(out) :: fract !< fraction between layers


! Local storage
integer  :: i

do i = 2, grid_handle%nlon
   if (lon <= grid_handle%longitude(i)) then
      bot = i-1
      top = i
      fract = (lon - grid_handle%longitude(bot)) / &
              (grid_handle%longitude(top) - grid_handle%longitude(bot))
      return
   endif
enddo

write(msgstring, *) 'looking for lon ', lon
call error_handler(E_ERR, 'lon_bounds', 'reached end of loop without finding lon', &
                   source, revision, revdate, text2=msgstring)

end subroutine lon_bounds

!-------------------------------------------------------------

!> Given a latitude lat and the grid_handle which contains both the 
!> 1D array of latitudes and the grid latitude count, returns the
!> indices of the grid below and above the location latitude and 
!> the fraction of the distance between. istatus is returned as 0 
!> unless the location latitude is south of the southernmost grid 
!> point (1 returned) or north of the northernmost (2 returned),
!> which may not be possible anymore and possibly could be removed.

subroutine lat_bounds(lat, grid_handle, bot, top, fract, istatus)

real(r8),        intent(in)  :: lat !< input latitude
type(grid_type), intent(in)  :: grid_handle  !< geo or mag grid
integer,         intent(out) :: bot !< index of bottom layer
integer,         intent(out) :: top !< index of top layer
real(r8),        intent(out) :: fract !< fraction between layers
integer,         intent(out) :: istatus !< return status


! Local storage
integer :: i

! Success should return 0, failure a positive number.
istatus = 0

! Check for too far south or north
if(lat < grid_handle%latitude(1)) then
   istatus = 1
   return
else if(lat > grid_handle%latitude(grid_handle%nlat)) then
   istatus = 2
   return
endif

! In the middle, search through
do i = 2, grid_handle%nlat
   if(lat <= grid_handle%latitude(i)) then
      bot = i - 1
      top = i
      fract = (lat - grid_handle%latitude(bot)) / &
              (grid_handle%latitude(top) - grid_handle%latitude(bot))
      return
   endif
enddo

write(msgstring, *) 'looking for lat ', lat
call error_handler(E_ERR, 'lat_bounds', 'reached end of loop without finding lat', &
                   source, revision, revdate, text2=msgstring)

end subroutine lat_bounds

!-------------------------------------------------------------

!> Given a latitude lat, the grid handle which contains the 1d array of 
!> colatitudes for the grid and the grid count, return the indices of
!> the grid below and above the location colatitude and the fraction 
!> of the distance between. colatitudes start at 0 and go to 180, but
!> to be consistent with our locations mod we have already transformed
!> them into 90 to -90.  this routine has to be different because the
!> order of the points is north pole to south, while latitudes are ordered
!> south pole to north.  we have to search in a different order and the
!> test itself is reversed from the lat_bounds() routine.
!> istatus is returned as 0 unless the location latitude is 
!> south of the southernmost grid point (1 returned) or north of the 
!> northernmost (2 returned). given our locations module i believe this
!> test is no longer needed since the grid includes the poles.

subroutine colat_bounds(lat, grid_handle, bot, top, fract, istatus)

real(r8),        intent(in)  :: lat          !< input latitude
type(grid_type), intent(in)  :: grid_handle  !< geo or mag grid
integer,         intent(out) :: bot          !< index of bottom layer
integer,         intent(out) :: top          !< index of top layer
real(r8),        intent(out) :: fract        !< fraction between layers
integer,         intent(out) :: istatus      !< return status


! Local storage
integer :: i

! Success should return 0, failure a positive number.
istatus = 0

! Check for too far south or north
if(lat > grid_handle%latitude(1)) then
   istatus = 1
   return
else if(lat < grid_handle%latitude(grid_handle%nlat)) then
   istatus = 2
   return
endif

! In the middle, search through
do i = 2, grid_handle%nlat
   if(lat >= grid_handle%latitude(i)) then
      bot = i - 1
      top = i
      fract = (lat - grid_handle%latitude(bot)) / &
              (grid_handle%latitude(top) - grid_handle%latitude(bot))
      return
   endif
enddo

write(msgstring, *) 'looking for colat ', lat
call error_handler(E_ERR, 'colat_bounds', 'reached end of loop without finding colat', &
                   source, revision, revdate, text2=msgstring)

end subroutine colat_bounds

!------------------------------------------------------------

!> find the index top and bottom index for a variable given an lheight and an
!> array of heights.

subroutine height_bounds(lheight, nheights, hgt_array, bot, top, fract, istatus)

real(r8),   intent(in)  :: lheight             !< height location
integer,    intent(in)  :: nheights            !< number of total heights
real(r8),   intent(in)  :: hgt_array(nheights) !< array of heights
integer,    intent(out) :: bot                 !< bottom bounding height
integer,    intent(out) :: top                 !< top bounding height
real(r8),   intent(out) :: fract               !< fraction inbetween
integer,    intent(out) :: istatus             !< return status


! Local variables
integer   :: i

! Succesful istatus is 0
! Make any failure here return istatus in the 20s
istatus = 0

if(lheight <= hgt_array(1)) then
   bot = 1
   top = 2
   ! NOTE: the fract definition is the relative distance from bottom to top
   fract = 1.0_r8 
   return
endif

! Search through the boxes
do i = 2, nheights
   ! If the location is lower than this entry, it must be in this box
   if(lheight < hgt_array(i)) then
      top = i
      bot = i -1
      fract = (lheight - hgt_array(bot)) / (hgt_array(top) - hgt_array(bot))
      return
   endif
enddo

! Falling off the end means the location is higher than the model top.
bot   = -1
top   = -1
fract = -1.0_r8

istatus = 20

end subroutine height_bounds

!------------------------------------------------------------------

!> Returns the the time step of the model; the smallest increment
!> in time that the model is capable of advancing the state in a given
!> implementation. In fact this sets the assimilation window size.

function get_model_time_step()

type(time_type) :: get_model_time_step !< returned timestep

if ( .not. module_initialized ) call static_init_model

get_model_time_step = model_timestep

end function get_model_time_step

!------------------------------------------------------------------

!> Given an integer index into the state vector structure, returns the
!> associated location and (optionally) the generic kind.  For this model
!> we must always return geographic coordinates.  For fields on the
!> magnetic grid this requires a coordinate transformation.

subroutine get_state_meta_data(state_handle, index_in, location, var_type)

type(ensemble_type),           intent(in)  :: state_handle  !< state ensemble handle
integer(i8),                   intent(in)  :: index_in !< dart state index of interest
type(location_type),           intent(out) :: location !< locartion of interest
integer,             optional, intent(out) :: var_type !< optional dart kind return


! Local variables
real(r8) :: lat, lon, height
integer  :: lon_index, lat_index, height_index, local_var, var_id
integer  :: state_loc(3)

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, state_loc(1), state_loc(2), state_loc(3), var_id=var_id)
local_var = get_kind_index(domain_id, var_id)

lon_index    = state_loc(dim_order_list(var_id, VAR_LON_INDEX))
lat_index    = state_loc(dim_order_list(var_id, VAR_LAT_INDEX))
height_index = state_loc(dim_order_list(var_id, VAR_HGT_INDEX))

! we are getting a mapping array between magnetic -> geogrid

if ( get_grid_type(local_var) == MAGNETIC_GRID ) then
   lon = mag_grid%conv_2d_lon(lon_index, lat_index)
   lat = mag_grid%conv_2d_lat(lon_index, lat_index)
else
   lon = geo_grid%longitude(lon_index)
   lat = geo_grid%latitude(lat_index)
endif

! Here we are ASSUMING that electric potential is a 2D
! variable.  If this is NOT the case FIX HERE!
if (local_var == KIND_ELECTRIC_POTENTIAL) then
   height   = 0.0_r8
   location = set_location(lon, lat, height, VERTISUNDEF)
else
   height   = geo_grid%heights(1)
   location = set_location(lon, lat, height, VERTISHEIGHT)
endif

if (present(var_type)) then
   var_type = local_var
endif

end subroutine get_state_meta_data

!------------------------------------------------------------------

!> Shutdown and clean-up.

subroutine end_model()

call deallocate_grid_space(geo_grid)
call deallocate_grid_space(mag_grid)

end subroutine end_model

!------------------------------------------------------------------

!> open and return the netcdf file id of template file

function get_grid_template_fileid(filename)

character(len=*), intent(in) :: filename
integer :: get_grid_template_fileid


call nc_check( NF90_open(filename, NF90_NOWRITE, get_grid_template_fileid), &
                  'get_grid_template_fileid', 'open '//trim(filename))

end function get_grid_template_fileid

!------------------------------------------------------------------

!> get grid sizes given netcdf id, lon, lat and height name. results
!> are store in the provided grid handle

subroutine get_grid_sizes(ncFileID, grid_handle, lon_name, lat_name, height_name)

integer,                    intent(in)    :: ncFileID !< netcdf file id
type(grid_type),            intent(inout) :: grid_handle !< geo or mag grid
character(len=*),           intent(in)    :: lon_name !< longitude name
character(len=*),           intent(in)    :: lat_name !< latitude name
character(len=*), optional, intent(in)    :: height_name !< height name


grid_handle%nlon = get_dim(ncFileID,lon_name, 'get_grid_sizes')
grid_handle%nlat = get_dim(ncFileID,lat_name, 'get_grid_sizes')

if (present(height_name)) then
   grid_handle%nheight = get_dim(ncFileID, height_name, 'get_grid_sizes')
else
   grid_handle%nheight = 1
endif

end subroutine get_grid_sizes

!------------------------------------------------------------------

!> is_conv:  if true, read the data into the conversion grid.
!> otherwise read into the normal lat/lon arrays.
!> 
!> is_co_latitude:  if true, subtract 90 from the lat values
!> co_latitudes start at 0 at the north pole and go to 180 at the south.
!> "normal" latitudes for us are -90 at the south pole up to 90 at the north.

subroutine read_horiz_grid(ncFileID, grid_handle, lon_name, lat_name, is_conv, is_co_latitude)

integer,          intent(in)    :: ncFileID !< netcdf file id for grid
type(grid_type),  intent(inout) :: grid_handle !< grid handle to be read into
character(len=*), intent(in)    :: lon_name !< longitude variable name
character(len=*), intent(in)    :: lat_name !< latitude variable name
logical,          intent(in)    :: is_conv !< fill conversion grid
logical,          intent(in)    :: is_co_latitude !< is grid in co-latitude

 
if (is_conv) then

   call get_data(ncFileID, lon_name, grid_handle%conv_2d_lon, 'read_conv_horiz_grid')
   call get_data(ncFileID, lat_name, grid_handle%conv_2d_lat, 'read_conv_horiz_grid')

   if (is_co_latitude) then
      grid_handle%conv_2d_lat(:,:) = 90.0_r8 - grid_handle%conv_2d_lat(:,:)
      grid_handle%uses_colatitude = .true.
   else
      grid_handle%uses_colatitude = .false.
   endif

   where(grid_handle%conv_2d_lon < 0) grid_handle%conv_2d_lon = grid_handle%conv_2d_lon + 360.0_r8

else

   call get_data(ncFileID, lon_name, grid_handle%longitude, 'read_horiz_grid')
   call get_data(ncFileID, lat_name, grid_handle%latitude,  'read_horiz_grid')

   if (is_co_latitude) then
      grid_handle%latitude(:) = 90.0_r8 - grid_handle%latitude(:)
      grid_handle%uses_colatitude = .true.
   else
      grid_handle%uses_colatitude = .false.
   endif

   where(grid_handle%longitude < 0) grid_handle%longitude = grid_handle%longitude + 360.0_r8

endif

end subroutine read_horiz_grid

!------------------------------------------------------------------

!> Read the height array from a netcdf file.

subroutine read_vert_levels(ncFileID, grid_handle, height_name)

integer,          intent(in)    :: ncFileID !< netcdf id with vertical information
type(grid_type),  intent(inout) :: grid_handle !< geo or mag grid
character(len=*), intent(in)    :: height_name !< name of height variable


call get_data(ncFileID, height_name, grid_handle%heights, 'read_vert_levels')

end subroutine read_vert_levels

!------------------------------------------------------------------

!> Writes the model-specific attributes to a netCDF file.
!> This includes coordinate variables for the geometric and
!> magnetic grids.

function nc_write_model_atts( ncFileID, model_mod_writes_state_variables ) result (ierr)

integer, intent(in)  :: ncFileID !> netCDF file identifier
logical, intent(out) :: model_mod_writes_state_variables !< false if you want DART to write the state variables
integer              :: ierr !> return error status


!  Typical sequence for adding new dimensions,variables,attributes:
!  NF90_OPEN             ! open existing netCDF dataset
!     NF90_redef         ! put into define mode 
!     NF90_def_dim       ! define additional dimensions (if any)
!     NF90_def_var       ! define variables: from name, type, and dims
!     NF90_put_att       ! assign attribute values
!  NF90_ENDDEF           ! end definitions: leave define mode
!     NF90_put_var       ! provide values for variable
!  NF90_CLOSE            ! close: save updated netCDF dataset

!----------------------------------------------------------------------
! local variables 
!----------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

! for the dimensions and coordinate variables
integer :: NlonDimID, NlatDimID, NhgtDimID
integer :: geoLonVarID, geoLatVarID, geoHeightVarID
integer :: magLonVarID, magLatVarID, magHeightVarID
integer ::  coLonVarID,  coLatVarID


! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

character(len=128)  :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

!-------------------------------------------------------------------------------
! Have dart write out the prognostic variables
!-------------------------------------------------------------------------------
model_mod_writes_state_variables = .false. 

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
! Write Global Attributes 
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call add_string_att(ncFileID, NF90_GLOBAL, 'creation_date', str1,      filename)
call add_string_att(ncFileID, NF90_GLOBAL, 'model_source'  ,source,    filename)
call add_string_att(ncFileID, NF90_GLOBAL, 'model_revision',revision,  filename)
call add_string_att(ncFileID, NF90_GLOBAL, 'model_revdate' ,revdate,   filename)
call add_string_att(ncFileID, NF90_GLOBAL, 'model',        'openggcm', filename)

!----------------------------------------------------------------------------
! We need to output grid information
!----------------------------------------------------------------------------
! Define the new dimensions IDs
!----------------------------------------------------------------------------

NlonDimID = set_dim(ncFileID, 'geo_lon', geo_grid%nlon,    filename)
NlatDimID = set_dim(ncFileID, 'geo_lat', geo_grid%nlat,    filename)
NhgtDimID = set_dim(ncFileID, 'geo_hgt', geo_grid%nheight, filename)

!----------------------------------------------------------------------------
! Create the (empty) Coordinate Variables and the Attributes
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
! Write out Geographic Grid attributes
!----------------------------------------------------------------------------

call nc_check(NF90_def_var(ncFileID,name='geo_lon', xtype=NF90_real, &
              dimids=NlonDimID, varid=geoLonVarID),&
              'nc_write_model_atts', 'geo_lon def_var '//trim(filename))
call nc_check(NF90_put_att(ncFileID,  geoLonVarID, 'long_name', 'longitudes of geo grid'), &
              'nc_write_model_atts', 'geo_lon long_name '//trim(filename))

! Grid Latitudes
call nc_check(NF90_def_var(ncFileID,name='geo_lat', xtype=NF90_real, &
              dimids=NlatDimID, varid=geoLatVarID),&
              'nc_write_model_atts', 'geo_lat def_var '//trim(filename))
call nc_check(NF90_put_att(ncFileID,  geoLatVarID, 'long_name', 'latitudes of geo grid'), &
              'nc_write_model_atts', 'geo_lat long_name '//trim(filename))

! Heights
call nc_check(NF90_def_var(ncFileID,name='geo_heights', xtype=NF90_real, &
              dimids=NhgtDimID, varid= geoHeightVarID), &
              'nc_write_model_atts', 'geo_heights def_var '//trim(filename))
call nc_check(NF90_put_att(ncFileID, geoHeightVarID, 'long_name', 'height of geo grid'), &
              'nc_write_model_atts', 'geo_heights long_name '//trim(filename))

!----------------------------------------------------------------------------
! Write out Magnetic Grid attributes
!----------------------------------------------------------------------------

NlonDimID = set_dim(ncFileID, 'mag_lon', mag_grid%nlon,    filename)
NlatDimID = set_dim(ncFileID, 'mag_lat', mag_grid%nlat,    filename)
NhgtDimID = set_dim(ncFileID, 'mag_hgt', mag_grid%nheight, filename)

! Grid Longitudes
call nc_check(NF90_def_var(ncFileID,name='mag_lon', xtype=NF90_real, &
              dimids=NlonDimID, varid=magLonVarID),&
              'nc_write_model_atts', 'mag_lon def_var '//trim(filename))
call nc_check(NF90_put_att(ncFileID,  magLonVarID, 'long_name', 'longitudes mag grid'), &
              'nc_write_model_atts', 'mag_lon long_name '//trim(filename))

! Grid Latitudes
call nc_check(NF90_def_var(ncFileID,name='mag_lat', xtype=NF90_real, &
              dimids=NlatDimID, varid=magLatVarID),&
              'nc_write_model_atts', 'mag_lat def_var '//trim(filename))
call nc_check(NF90_put_att(ncFileID,  magLatVarID, 'long_name', 'latitudes of mag grid'), &
              'nc_write_model_atts', 'mag_lat long_name '//trim(filename))

! Heights
call nc_check(NF90_def_var(ncFileID,name='mag_heights', xtype=NF90_real, &
              dimids=NhgtDimID, varid= magHeightVarID), &
              'nc_write_model_atts', 'mag_heights def_var '//trim(filename))
call nc_check(NF90_put_att(ncFileID, magHeightVarID, 'long_name', 'height for mag grid'), &
              'nc_write_model_atts', 'mag_heights long_name '//trim(filename))

!----------------------------------------------------------------------------
! Write out Co-Latitude Grid attributes
!----------------------------------------------------------------------------

NlonDimID = set_dim(ncFileID, 'co_lon', mag_grid%nlon,    filename)
NlatDimID = set_dim(ncFileID, 'co_lat', mag_grid%nlat,    filename)

! Grid Longitudes
call nc_check(NF90_def_var(ncFileID,name='co_lon', xtype=NF90_real, &
              dimids=(/ NlonDimID, NlatDimID /), varid=coLonVarID),&
              'nc_write_model_atts', 'co_lon def_var '//trim(filename))
call nc_check(NF90_put_att(ncFileID,  coLonVarID, 'long_name', 'longitudes co grid'), &
              'nc_write_model_atts', 'co_lon long_name '//trim(filename))

! Grid Latitudes
call nc_check(NF90_def_var(ncFileID,name='co_lat', xtype=NF90_real, &
              dimids=(/ NlonDimID, NlatDimID /), varid=coLatVarID),&
              'nc_write_model_atts', 'co_lat def_var '//trim(filename))
call nc_check(NF90_put_att(ncFileID,  coLatVarID, 'long_name', 'latitudes of co grid'), &
              'nc_write_model_atts', 'co_lat long_name '//trim(filename))

! Finished with dimension/variable definitions, must end 'define' mode to fill.

call nc_check(NF90_enddef(ncfileID), 'prognostic enddef '//trim(filename))

!----------------------------------------------------------------------------
! Fill Geographic Grid values
!----------------------------------------------------------------------------

call nc_check(NF90_put_var(ncFileID, geoLonVarID, geo_grid%longitude ), &
             'nc_write_model_atts', 'geo_lon put_var '//trim(filename))
call nc_check(NF90_put_var(ncFileID, geoLatVarID, geo_grid%latitude ), &
             'nc_write_model_atts', 'geo_lat put_var '//trim(filename))
call nc_check(NF90_put_var(ncFileID, geoHeightVarID, geo_grid%heights ), &
             'nc_write_model_atts', 'geo_heights put_var '//trim(filename))

!----------------------------------------------------------------------------
! Fill Magnetic Grid values
!----------------------------------------------------------------------------

call nc_check(NF90_put_var(ncFileID, magLonVarID, mag_grid%longitude ), &
             'nc_write_model_atts', 'mag_lon put_var '//trim(filename))
call nc_check(NF90_put_var(ncFileID, magLatVarID, mag_grid%latitude ), &
             'nc_write_model_atts', 'mag_lat put_var '//trim(filename))
call nc_check(NF90_put_var(ncFileID, magHeightVarID, mag_grid%heights ), &
             'nc_write_model_atts', 'mag_heights put_var '//trim(filename))

!----------------------------------------------------------------------------
! Fill Co-Latitude Grid values
!----------------------------------------------------------------------------

call nc_check(NF90_put_var(ncFileID, coLonVarID, mag_grid%conv_2d_lon(:,:)), &
             'nc_write_model_atts', 'co_lon put_var '//trim(filename))
call nc_check(NF90_put_var(ncFileID, coLatVarID, mag_grid%conv_2d_lat(:,:) ), &
             'nc_write_model_atts', 'co_lat put_var '//trim(filename))

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(NF90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts

!------------------------------------------------------------------

!> Do not need to to use this since state_space_diag_mod takes care 
!> of writing out the state variables. This is because :
!>
!>   model_mod_writes_state_variables = .false. in nc_write_model_atts

function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
integer,                intent(in) :: ncFileID !< netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec !< dart state vector
integer,                intent(in) :: copyindex !< copy number index
integer,                intent(in) :: timeindex !< time index from model
integer                            :: ierr !< return value of function


if ( .not. module_initialized ) call static_init_model

ierr = 0

end function nc_write_model_vars

!------------------------------------------------------------------

!> Perturbs state copies for generating initial ensembles.

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: state_ens_handle !< state ensemble handle
integer,             intent(in)    :: ens_size !< ensemble size
real(r8),            intent(in)    :: pert_amp !< perterbation amplitude
logical,             intent(out)   :: interf_provided !< have you provided an interface?


! Local Variables
integer     :: i, j
integer(i8) :: dart_index

! Storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq

if ( .not. module_initialized ) call static_init_model

interf_provided = .true.

! Initialize random number sequence
call init_random_seq(random_seq, my_task_id())

! perturb the state using random gaussian noise
do i=1,state_ens_handle%my_num_vars
   dart_index = state_ens_handle%my_vars(i)
   do j=1, ens_size
      ! Since potential and electron density have such radically different values
      ! we weight the standard deviation with the actual state value so that noise 
      ! created is closer to the actual values in the state. NOTE: If the value is
      ! state value is zero the standard deviation will be zero and your values will
      ! not be perturbed.
      state_ens_handle%copies(j,i) = random_gaussian(random_seq, &
         mean               = state_ens_handle%copies(j,i),      &
         standard_deviation = state_ens_handle%copies(j,i)*model_perturbation_amplitude)
   enddo
enddo


end subroutine pert_model_copies

!------------------------------------------------------------------

!> Given a DART location (referred to as "base") and a set of candidate
!> locations & kinds (obs, obs_kind), returns the subset close to the
!> "base", their indices, and their distances to the "base" ...
!>
!> For vertical distance computations, general philosophy is to convert all
!> vertical coordinates to a common coordinate. This coordinate type is defined
!> in the namelist with the variable "vert_localization_coord".

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, &
                         obs, obs_kind, num_close, close_ind, dist, state_handle)

type(ensemble_type),               intent(in) :: state_handle !< state ensemble handle
type(get_close_type),              intent(in) :: gc !< get_close_type handle
type(location_type),               intent(in) :: base_obs_loc !< base observation location
integer,                           intent(in) :: base_obs_kind !< base dart observation kind
type(location_type), dimension(:), intent(in) :: obs !< observation sequence
integer,             dimension(:), intent(in) :: obs_kind !< dart kind
integer,                           intent(out):: num_close !< number of close found
integer,             dimension(:), intent(out):: close_ind !< list of close indicies
real(r8), optional,  dimension(:), intent(out):: dist !< list of distances of close observations


call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs, obs_kind, &
                       num_close, close_ind, dist)

end subroutine get_close_obs

!------------------------------------------------------------------

!> do a 2d horizontal interpolation for the value at the bottom level, 
!> then again for the top level, then do a linear interpolation in the 
!> vertical to get the final value.

subroutine do_interp(state_handle, ens_size, grid_handle, hgt_bot, hgt_top, hgt_fract, &
                     llon, llat, obs_kind, expected_obs, istatus)

type(ensemble_type), intent(in)  :: state_handle !< state ensemble handle
integer,             intent(in)  :: ens_size !< ensemble size
type(grid_type),     intent(in)  :: grid_handle !< geo or mag grid
integer,             intent(in)  :: hgt_bot !< index to bottom bound
integer,             intent(in)  :: hgt_top !< index to top bound
real(r8),            intent(in)  :: hgt_fract !< fraction inbetween top and bottom
real(r8),            intent(in)  :: llon !< longitude to interpolate
real(r8),            intent(in)  :: llat !< latitude to interpolate
integer,             intent(in)  :: obs_kind !< dart kind
real(r8),            intent(out) :: expected_obs(ens_size) !< interpolated value
integer,             intent(out) :: istatus(ens_size) !< status of interpolation


! Local Variables
real(r8)    :: bot_val(ens_size), top_val(ens_size)
integer     :: temp_status(ens_size)
logical     :: return_now

istatus(:) = 0 ! need to start with istatus = 0 for track status to work properly

call lon_lat_interpolate(state_handle, ens_size, grid_handle, obs_kind, &
                         llon, llat, hgt_bot, bot_val, temp_status)
call track_status(ens_size, temp_status, bot_val, istatus, return_now)
if (return_now) return

call lon_lat_interpolate(state_handle, ens_size, grid_handle, obs_kind, &
                         llon, llat, hgt_top, top_val, temp_status)
call track_status(ens_size, temp_status, top_val, istatus, return_now)
if (return_now) return

! Then weight them by the vertical fraction and return
where (istatus == 0) 
   expected_obs = bot_val + hgt_fract * (top_val - bot_val)
elsewhere
   expected_obs = MISSING_R8
endwhere

end subroutine do_interp

!--------------------------------------------------------------------

!> construct restart file name for reading
!>
!> stub is found in input.nml io_filename_nml
!> restart files typically are of the form openggcm3D_0001.nc

function construct_file_name_in(stub, domain, copy)

character(len=512), intent(in) :: stub !< stubname of file
integer,            intent(in) :: domain !< domain of file
integer,            intent(in) :: copy !< copy number (i.e. ensemble number)
character(len=1024)            :: construct_file_name_in !< constructed filename


write(construct_file_name_in, '(A,''_'',I4.4,''.nc'')') trim(stub), copy

end function construct_file_name_in

!--------------------------------------------------------------------

!> read the current model time for this data from template file

function read_model_time(filename)

character(len=*), intent(in) :: filename !< file to get time
type(time_type) :: read_model_time !< returned time from file


! netcdf variables
integer :: ncFileID, VarID

! time variables
integer :: seconds
type(time_type) :: base_time

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(msgstring,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'read_model_time',msgstring,source,revision,revdate)
endif

call nc_check( NF90_open(filename, NF90_NOWRITE, ncFileID), &
                  'read_model_time', 'open '//filename )

call nc_check(NF90_inq_varid(ncFileID, 'time', VarID), &
              'read_model_time', 'time inq_varid')

call nc_check(NF90_get_var(ncFileID, VarID, seconds), &
              'read_model_time', 'time get_var')

!>@todo FIXME : Should be grabbing the base model time from the 
!>              time variable attributes. This is hardcoded for now.
base_time = set_date(1966,1,1,0,0)

read_model_time = base_time + set_time(seconds)

end function read_model_time

!--------------------------------------------------------------------

!> pass the vertical localization coordinate to assim_tools_mod

function query_vert_localization_coord()

integer :: query_vert_localization_coord !< return which height we want to 
                                         !< localize in

query_vert_localization_coord = VERTISHEIGHT

end function query_vert_localization_coord

!--------------------------------------------------------------------

!> This is used in the filter_assim. The vertical conversion is done using the 
!> mean state.  Currently this model is not doing any vertical conversions.

subroutine vert_convert(state_handle, location, obs_kind, istatus)

type(ensemble_type), intent(in)  :: state_handle !< state ensemble handle
type(location_type), intent(in)  :: location !< location to convert
integer,             intent(in)  :: obs_kind !< dart kind
integer,             intent(out) :: istatus !< status of conversion

istatus = 0

end subroutine vert_convert

!------------------------------------------------------------

!> Always an error to call this routine.
!> At present, this is only called if the namelist parameter 
!> start_from_restart is set to .false. in the program perfect_model_obs.
!> Returns a model state vector, x, that is some sort of appropriate
!> initial condition for starting up a long integration of the model.
!> This is not possible for a large geophysical model.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:) !< initalize state from scratch


msgstring2 = "cannot run perfect_model_obs with 'start_from_restart = .false.' "
msgstring3 = 'use openggcm_to_dart to generate an initial state'
call error_handler(E_ERR,'init_conditions', &
                  'ERROR!!  openggcm model has no built-in default state', &
                  source, revision, revdate, &
                  text2=msgstring2, text3=msgstring3)

! this code never reached - just here to avoid compiler warnings
! about an intent(out) variable not being set to a value.
x = 0.0_r8

end subroutine init_conditions

!------------------------------------------------------------------

!> If the model could be called as a subroutine, does a single
!> timestep advance.  openggcm cannot be called this way, so fatal error
!> if this routine is called.

subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:) !< state vector
type(time_type), intent(in)    :: time !< time stamp for state_vector

call error_handler(E_ERR,'adv_1step', &
                  'openggcm model cannot be called as a subroutine; async cannot = 0', &
                  source, revision, revdate)

end subroutine adv_1step

!------------------------------------------------------------------

!> Companion interface to init_conditions. Returns a time that is
!> appropriate for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter 
!> start_from_restart is set to .false. in the program perfect_model_obs.

subroutine init_time(time)

type(time_type), intent(out) :: time !< time restart file


msgstring2 = "cannot run perfect_model_obs with 'start_from_restart = .false.' "
msgstring3 = 'use openggcm_to_dart to generate an initial state which contains a timestamp'
call error_handler(E_ERR,'init_time', &
                  'ERROR!!  openggcm model has no built-in default time', &
                  source, revision, revdate, &
                  text2=msgstring2, text3=msgstring3)

! this code never reached - just here to avoid compiler warnings
! about an intent(out) variable not being set to a value.
time = set_time(0,0)

end subroutine init_time

!------------------------------------------------------------------

!> Verify that the namelist was filled in correctly, and check
!> that there are valid entries for the dart_kind. 
!> Returns a table with columns:  
!>
!>    netcdf_variable_name ; dart_kind_string ; update_string ; grid_id
!>

subroutine verify_state_variables( state_variables, ngood, table, kind_list, update_var, grid_id )

character(len=*), intent(inout) :: state_variables(:) !< list of state variables and attributes from nml
integer,          intent(out)   :: ngood !< number of good namelist values or nfields
character(len=*), intent(out)   :: table(:,:) !< 2d table with information from namelist
integer,          intent(out)   :: kind_list(:) !< dart kind
logical,          intent(out)   :: update_var(:) !< list of logical update information
integer,          intent(out)   :: grid_id(:) !< list of dart kind numbers


! Local Variables
integer :: nrows, i
character(len=NF90_MAX_NAME) :: varname, dartstr, update, gridname

nrows = size(table,1)

ngood = 0

if ( state_variables(1) == ' ' ) then ! no model_state_variables namelist provided
   msgstring = 'model_nml:model_state_variables not specified using default variables'
   call error_handler(E_ERR,'verify_state_variables',msgstring,source,revision,revdate)
endif

MyLoop : do i = 1, nrows

   varname  = trim(state_variables(4*i -3))
   dartstr  = trim(state_variables(4*i -2))
   update   = trim(state_variables(4*i -1))
   gridname = trim(state_variables(4*i   ))
   
   call to_upper(update)

   table(i,1) = trim(varname)
   table(i,2) = trim(dartstr)
   table(i,3) = trim(update)
   table(i,4) = trim(gridname)

   if ( table(i,1) == ' ' .and. &
        table(i,2) == ' ' .and. &
        table(i,3) == ' ' .and. &
        table(i,4) == ' ') exit MyLoop ! Found end of list.

   if ( table(i,1) == ' ' .or. &
        table(i,2) == ' ' .or. &
        table(i,3) == ' ' .or. &
        table(i,4) == ' ') then
      msgstring = 'model_nml:model_state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',msgstring,source,revision,revdate)
   endif

   ! Make sure DART kind is valid

   kind_list(i) = get_raw_obs_kind_index(dartstr)
   if( kind_list(i)  < 0 ) then
      write(msgstring,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',msgstring,source,revision,revdate)
   endif
   
   ! Make sure the update variable has a valid name

   SELECT CASE (update)
      CASE ('UPDATE')
         update_var(i) = .true.
      CASE ('NO_COPY_BACK')
         update_var(i) = .false.
      CASE DEFAULT
         write(msgstring, '(A)')  'only UPDATE or NO_COPY_BACK supported in model_state_variable namelist'
         write(msgstring2,'(6A)') 'you provided : ', trim(varname), ', ', trim(dartstr), ', ', trim(update), ', ', trim(gridname)
         call error_handler(E_ERR,'verify_state_variables',msgstring,source,revision,revdate, text2=msgstring2)
   END SELECT

   ! Make sure the update variable has a valid name

   SELECT CASE (gridname)
      CASE ('GEOGRAPHIC_GRID')
         grid_id(i) = GEOGRAPHIC_GRID
      CASE ('MAGNETIC_GRID')
         grid_id(i) = MAGNETIC_GRID
      CASE DEFAULT
         write(msgstring, '(A)')  'only GEOGRAPHIC_GRID or MAGNETIC_GRID supported in model_state_variable namelist'
         write(msgstring2,'(8A)') 'you provided : ',&
                               trim(varname), ', ', trim(dartstr), ', ', trim(update), ', ', trim(gridname)
         call error_handler(E_ERR,'verify_state_variables',msgstring,source,revision,revdate, text2=msgstring2)
   END SELECT

   ! Record the contents of the DART state vector

   if (do_output()) then
      write(msgstring,'(A,I2,8A)') 'variable ',i,' is ',trim(varname), ', ', &
                                  trim(dartstr), ', ', trim(update), ', ', trim(gridname)
      call error_handler(E_MSG,'verify_state_variables',msgstring,source,revision,revdate)
   endif

   ngood = ngood + 1
enddo MyLoop

end subroutine verify_state_variables

!----------------------------------------------------------------------

!> determine if grid is geographic or magnetic given a dart kind
!> this is based on what the user tells us in the namelist.

function get_grid_type(dart_kind)

integer, intent(in) :: dart_kind !< dart kind
integer :: get_grid_type !< returned grid


! Local Variables
integer :: i

do i = 1,nfields
   if ( dart_kind == state_kinds_list(i) ) then
      get_grid_type = grid_info_list(i)
      return
   endif
enddo

write(msgstring,*)' Can not find grid type for : ', get_raw_obs_kind_name(dart_kind)
call error_handler(E_ERR,'get_grid_type',msgstring,source,revision,revdate)

end function get_grid_type

!----------------------------------------------------------------------

!> track_status can be used to keep track of the status of
!> each ensemble member during multiple calls to model_interpolate
!> for a given obs_def.
!> It assumes that you are starting with istatus(:) = 0
!> If debugging, return_now is only set to true if all istatuses are non-zero
!> If not debugging, return_now is set to true if any istatues are non-zero
!>  and any remaining zero istatues are set to 1.

!> @todo FIXME: this should be in the utilities mod!!!

subroutine track_status(ens_size, val_istatus, val_data, istatus, return_now)

integer,  intent(in)    :: ens_size
integer,  intent(in)    :: val_istatus(ens_size)
real(r8), intent(inout) :: val_data(ens_size) !> expected_obs for obs_def
integer,  intent(inout) :: istatus(ens_size) !> istatus for obs_def
logical,  intent(out)   :: return_now


where (istatus == 0) istatus = val_istatus
where (istatus /= 0) val_data = MISSING_R8

return_now = .false.
if (debug > 0) then
   if( all(istatus /= 0))then
      return_now = .true.
      val_data(:) = missing_r8
   endif
else
   if( any(istatus /= 0)) then
      return_now = .true.
      val_data(:) = missing_r8
      where (istatus == 0) istatus = 1
   endif
endif

end subroutine track_status

!----------------------------------------------------------------------

!> Allocate space for grid variables. 
!> cannot be called until the grid sizes are set.

subroutine allocate_grid_space(grid_handle, conv)

type(grid_type), intent(inout) :: grid_handle !< geo or mag grid handle
logical,         intent(in)    :: conv !< if true, the grid has conversion arrays


allocate(grid_handle%longitude(grid_handle%nlon))
allocate(grid_handle%latitude(grid_handle%nlat))
allocate(grid_handle%heights(grid_handle%nheight))

if (conv) then
   allocate(grid_handle%conv_2d_lon(grid_handle%nlon,grid_handle%nlat))
   allocate(grid_handle%conv_2d_lat(grid_handle%nlon,grid_handle%nlat))
endif

end subroutine allocate_grid_space

!----------------------------------------------------------------------

!> Deallocate space for grid variables. 

subroutine deallocate_grid_space(grid_handle)

type(grid_type), intent(inout) :: grid_handle !< geo or mag grid handle


if (allocated(grid_handle%longitude))  deallocate(grid_handle%longitude)
if (allocated(grid_handle%latitude))   deallocate(grid_handle%latitude)
if (allocated(grid_handle%heights))    deallocate(grid_handle%heights)

if (allocated(grid_handle%conv_2d_lon)) deallocate(grid_handle%conv_2d_lon)
if (allocated(grid_handle%conv_2d_lat)) deallocate(grid_handle%conv_2d_lat)

end subroutine deallocate_grid_space

!------------------------------------------------------------------

!> @todo FIXME: the following routines should be in a netcdf utils 
!> module somewhere

!------------------------------------------------------------------

!< returns the length of the dimension from a given a dimension name
!< and netcdf file id

function get_dim(ncFileID, dim_name, context)

integer,          intent(in)    :: ncFileID !< netcdf file id
character(len=*), intent(in)    :: dim_name !< dimension name of interest
character(len=*), intent(in)    :: context !< routine calling from
integer :: get_dim !< returns the length of the dimension


! netcdf variables
integer :: DimID, rc

rc = NF90_inq_dimid(ncid=ncFileID, name=dim_name, dimid=DimID)
call nc_check(rc, trim(context)//' inquiring for dimension '//trim(dim_name))

rc = NF90_inquire_dimension(ncFileID, DimID, len=get_dim)
call nc_check(rc, trim(context)//' getting length of dimension '//trim(dim_name))

end function get_dim

!------------------------------------------------------------------

!< sets a netcdf file dimension given a dimension name and length.
!< returns the value of the netcdf dimension id.

function set_dim(ncFileID, dim_name, dim_val, context)

integer,          intent(in)    :: ncFileID !< netcdf id
character(len=*), intent(in)    :: dim_name !< dimension name
integer,          intent(in)    :: dim_val !< dimension length
character(len=*), intent(in)    :: context !< routine calling from
integer :: set_dim !< netcdf id to the defined dimension


! netcdf variables
integer :: rc

rc = NF90_def_dim(ncid=ncFileID, name=dim_name, len=dim_val, dimid=set_dim)
call nc_check(rc, trim(context)//' setting dimension '//trim(dim_name))

end function set_dim

!------------------------------------------------------------------

!> read 1d variable data from netcdf file

subroutine get_data_1d(ncFileID, var_name, data_array, context)

integer,          intent(in)    :: ncFileID !< netcdf id
real(r8),         intent(out)   :: data_array(:) !< id array of values
character(len=*), intent(in)    :: var_name !< variable of interest
character(len=*), intent(in)    :: context !< routine called from
 

! netcdf variables
integer :: VarID, rc

rc = NF90_inq_varid(ncFileID, var_name, VarID)
call nc_check(rc, trim(context)//' inquiring for 1d array '//trim(var_name))

rc = NF90_get_var(ncFileID, VarID, data_array)
call nc_check(rc, trim(context)//' getting data for 1d array '//trim(var_name))

end subroutine get_data_1d

!------------------------------------------------------------------

!> read 2d variable data from netcdf file

subroutine get_data_2d(ncFileID, var_name, data_array, context)

integer,          intent(in)    :: ncFileID !< netcdf id
real(r8),         intent(out)   :: data_array(:,:) !< id array of values
character(len=*), intent(in)    :: var_name !< variable of interest
character(len=*), intent(in)    :: context !< routine called from
 

! netcdf variables
integer :: VarID, rc

rc = NF90_inq_varid(ncFileID, var_name, VarID)
call nc_check(rc, trim(context)//' inquiring for 2d array '//trim(var_name))

rc = NF90_get_var(ncFileID, VarID, data_array)
call nc_check(rc, trim(context)//' getting data for 2d array '//trim(var_name))

end subroutine get_data_2d

!----------------------------------------------------------------------

!> read 3d variable data from netcdf file

subroutine get_data_3d(ncFileID, var_name, data_array, context)

integer,          intent(in)    :: ncFileID !< netcdf id
real(r8),         intent(out)   :: data_array(:,:,:) !< id array of values
character(len=*), intent(in)    :: var_name !< variable of interest
character(len=*), intent(in)    :: context !< routine called from
 

! netcdf variables
integer :: VarID, rc

rc = NF90_inq_varid(ncFileID, var_name, VarID)
call nc_check(rc, trim(context)//' inquiring for 3d array '//trim(var_name))

rc = NF90_get_var(ncFileID, VarID, data_array)
call nc_check(rc, trim(context)//' getting data for 3d array '//trim(var_name))

end subroutine get_data_3d

!----------------------------------------------------------------------

!> add a string attribute to a netcdf variable

subroutine add_string_att(ncFileID, varid, attname, attval, context)

integer,          intent(in)    :: ncFileID !< netcdf file id
integer,          intent(in)    :: varid !< netcdf variable id
character(len=*), intent(in)    :: attname !< netcdf attribute name
character(len=*), intent(in)    :: attval !< netcdf attrivute value
character(len=*), intent(in)    :: context !< routine called by


!netcdf variables
integer :: rc

rc = NF90_put_att(ncFileID, varid, attname, attval)
call nc_check(rc, trim(context)//' putting attribute '//trim(attname))

end subroutine add_string_att

!----------------------------------------------------------------------

!> initialize a grid transformation type.  the transform type
!> and conversion routines are in the cort_mod.  the grid transformations
!> change with time, so the initialization routine must know the
!> current model time.

subroutine initialize_openggcm_transform(filename)

character(len=*), intent(inout) :: filename !< file with openggcm time information


! Local Variables
type(time_type) :: dart_time
integer :: yr,mo,dy,hr,mn,se

dart_time = read_model_time(filename)

call get_date(dart_time, yr, mo, dy, hr, mn, se)

call cotr_set(yr,mo,dy,hr,mn,real(se,r4),openggcm_transform)

end subroutine initialize_openggcm_transform

!----------------------------------------------------------------------

!> transform grid from geo->magnetic or vice-versa. 
!> currently we only need it to translate from magnetic to geo.
!> direction options GEO_TO_SM or SM_TO_GEO

subroutine transform_mag_geo(llon, llat, lheight, direction)

real(r8), intent(inout) :: llon !< dart longitude [0,360]
real(r8), intent(inout) :: llat !< dart latitude [-90,90]
real(r8), intent(inout) :: lheight !< height
integer,  intent(in)    :: direction


! Local Variables
real(r4) :: xin, xout, yin, yout, zin, zout
character(len=4) :: dirstrin, dirstrout

if (direction == SM_TO_GEO) then
   dirstrin  = 'sm '
   dirstrout = 'geo'
else if (direction == GEO_TO_SM) then
   dirstrin  = 'geo'
   dirstrout = 'sm '
else
  call error_handler(E_ERR, 'transform_mag_geo', &
          'unexpected direction input, should not happen', &
           source, revision, revdate)
endif

! cort_mod is expecting height from center of earth
if (lheight == MISSING_R8) then
   lheight = 0.0_r8 + earth_radius
else
   lheight = lheight + earth_radius
endif

! degxyz requires longitude to be between -180 and 180
if (llon > 180.0_r8) llon = llon - 360.0_r8
! degxyz requires latitude to be between 0 and 180
llat = 90.0_r8 - llat

! transform spherical coordinates to cartesian
call degxyz(lheight, llon, llat, xin, yin, zin)

! transform from geographic to magnetic grid or back
call cotr(openggcm_transform, dirstrin, dirstrout, &
          xin, yin, zin, xout, yout, zout)

! transform cartesian coordinates to spherical 
call xyzdeg(xout, yout, zout, lheight, llon, llat)

! transform back to dart longitude coordinates [0,360]
if (llon < 0.0_r8) llon = llon + 360.0_r8
! transform back to dart latitude coordinates [-90,90]
llat = 90.0_r8 - llat

end subroutine transform_mag_geo

!----------------------------------------------------------------------

!> recording the storage order of the dimensions for each variable.
!>
!> lon_index is   dim_order_list(VAR_ID, VAR_LON_INDEX)
!> lat_index is   dim_order_list(VAR_ID, VAR_LAT_INDEX)
!> hgt_index is   dim_order_list(VAR_ID, VAR_HGT_INDEX)
!>
!> variables without a height dimension, VAR_HGT_INDEX is set to one.

subroutine make_dim_order_table(ngood)
integer, intent(in) :: ngood !< number of good fields


! Local Variables
integer :: ivar, jdim
character(len=NF90_MAX_NAME) :: dimname

! initialize list
dim_order_list(:,:) = 1

do ivar = 1,ngood
   do jdim = 1,get_num_dims(domain_id, ivar)
      dimname = get_dim_name(domain_id, ivar, jdim)
      SELECT CASE (trim(dimname))
         CASE ('cg_lon','ig_lon')
            dim_order_list(ivar, VAR_LON_INDEX) = jdim
         CASE ('cg_lat','ig_lat')
            dim_order_list(ivar, VAR_LAT_INDEX) = jdim
         CASE ('cg_height','ig_height')
            dim_order_list(ivar, VAR_HGT_INDEX) = jdim
         CASE DEFAULT
            write(msgstring,*) 'cannot find dimension ', trim(dimname),&
                               ' for variable', get_variable_name(domain_id, ivar)
            call error_handler(E_ERR,'make_dim_order_table',msgstring,source,revision,revdate)
      END SELECT
   enddo
enddo

end subroutine make_dim_order_table

!----------------------------------------------------------------------
!------------------------------------------------------------------
! End of model_mod
!------------------------------------------------------------------

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
