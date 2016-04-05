! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!-----------------------------------------------------------------------
!>
!> This is the interface between JULES and DART.
!> There are 16 required public interfaces whose arguments CANNOT be changed.
!> There are potentially many more public routines that are typically
!> used by the converter programs. As the converter programs get phased out
!> with the impending native netCDF read/write capability, these extra
!> public interfaces may not need to be public. 
!>

module model_mod

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, SECPERDAY, MISSING_R8,                    &
                             MISSING_I, MISSING_R4, rad2deg, deg2rad, PI,      &
                             obstypelength
use time_manager_mod, only : time_type, set_time, get_time, set_date, get_date,&
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+),  operator(-),          &
                             operator(>),  operator(<),  operator(/),          &
                             operator(/=), operator(<=), operator(==)

use     location_mod, only : location_type, get_dist, query_location,          &
                             set_location, get_location, horiz_dist_only,      &
                             vert_is_undef,    VERTISUNDEF,                    &
                             vert_is_surface,  VERTISSURFACE,                  &
                             vert_is_level,    VERTISLEVEL,                    &
                             vert_is_pressure, VERTISPRESSURE,                 &
                             vert_is_height,   VERTISHEIGHT,                   &
                             get_close_maxdist_init, get_close_type,           &
                             get_close_obs_init, get_close_obs_destroy,        &
                             loc_get_close_obs => get_close_obs,               &
                             LocationDims, write_location,                     &
                             is_location_in_region

use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper,                    &
                             find_namelist_in_file, check_namelist_read,       &
                             file_exist, find_textfile_dims, file_to_text,     &
                             open_file, close_file

use     obs_kind_mod, only : KIND_SOIL_TEMPERATURE,   &
                             KIND_SOIL_MOISTURE,      &
                             KIND_LIQUID_WATER,       &
                             KIND_ICE,                &
                             KIND_SNOWCOVER_FRAC,     &
                             KIND_SNOW_THICKNESS,     &
                             KIND_LEAF_CARBON,        &
                             KIND_WATER_TABLE_DEPTH,  &
                             KIND_GEOPOTENTIAL_HEIGHT,&
                             paramname_length,        &
                             get_raw_obs_kind_index,  &
                             get_raw_obs_kind_name

use mpi_utilities_mod, only: my_task_id

use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use typesizes
use netcdf

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          static_init_model,      &
          end_model,              &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_state,       &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          ens_mean_for_model

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.

public :: jules_to_dart_state_vector,   &
          dart_to_jules_restart,        &
          get_jules_restart_filename,   &
          get_state_time,               &
          DART_get_var,                 &
          get_model_time

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=256) :: string1, string2, string3
logical, save :: module_initialized = .false.

! Storage for a random sequence for perturbing a single initial state

type(random_seq_type) :: random_seq

! -----------------------------------------------------------------
!
!  The DART state vector may consist of things like:
!
!  The variables in the JULES restart file that are used to create the
!  DART state vector are specified in the input.nml:model_nml namelist.
!
! -----------------------------------------------------------------

! A gridcell may consist of up to 9 tiles.
! 5 of the tiles are PFTs,
! 4 of the tiles are non-vegetated classifications.
! The tiles occupy a fraction of the gridcell.
! The fractions for each tile in each grid cell are 
! contained in a file specified in 'ancillaries.nml'
!
! tile types

integer, parameter ::  BROADLEAF_TREE   = 1
integer, parameter ::  NEEDLE_LEAF_TREE = 2
integer, parameter ::  C3_GRASS         = 3
integer, parameter ::  C4_GRASS         = 4
integer, parameter ::  SHRUBS           = 5
integer, parameter ::  URBAN            = 6
integer, parameter ::  INLAND_WATER     = 7
integer, parameter ::  BARE_SOIL        = 8
integer, parameter ::  LAND_ICE         = 9

! Codes for restricting the range of a variable
integer, parameter :: BOUNDED_NONE  = 0 ! ... unlimited range
integer, parameter :: BOUNDED_BELOW = 1 ! ... minimum, but no maximum
integer, parameter :: BOUNDED_ABOVE = 2 ! ... maximum, but no minimum
integer, parameter :: BOUNDED_BOTH  = 3 ! ... minimum and maximum

integer :: nfields
integer, parameter :: max_state_variables = 20
integer, parameter :: num_state_table_columns = 6
character(len=obstypelength) :: variable_table(max_state_variables, num_state_table_columns)

! Codes for interpreting the columns of the variables
integer, parameter :: VT_VARNAMEINDX  = 1 ! ... variable name
integer, parameter :: VT_KINDINDX     = 2 ! ... DART kind
integer, parameter :: VT_MINVALINDX   = 3 ! ... minimum value, if any
integer, parameter :: VT_MAXVALINDX   = 4 ! ... maximum value, if any
integer, parameter :: VT_ORIGININDX   = 5 ! ... file of origin
integer, parameter :: VT_STATEINDX    = 6 ! ... update (state) or not

! things which can/should be in the model_nml

integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 60
real(r8)           :: model_perturbation_amplitude = 0.2
integer            :: debug = 0   ! turn up for more and more debug messages
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: jules_restart_filename = 'jules_restart.nc'
character(len=256) :: jules_output_filename = 'jules_history.nc'


character(len=obstypelength) :: variables(max_state_variables*num_state_table_columns) = ' '

namelist /model_nml/            &
   jules_restart_filename,      &
   jules_output_filename,       &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   calendar,                    &
   debug,                       &
   variables

! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=NF90_MAX_NAME) :: coordinates
   character(len=obstypelength), dimension(NF90_MAX_VAR_DIMS) :: dimnames
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer  :: numdims
   integer  :: rank        ! excludes the (singleton) time dimension if present
   integer  :: maxlevels
   integer  :: xtype
   integer  :: varsize     ! prod(dimlens(1:numdims))
   integer  :: index1      ! location in dart state vector of first occurrence
   integer  :: indexN      ! location in dart state vector of last  occurrence
   integer  :: dart_kind
   integer  :: rangeRestricted
   real(r8) :: minvalue
   real(r8) :: maxvalue
   integer  :: spvalINT, missingINT
   real(r4) :: spvalR4, missingR4
   real(r8) :: spvalR8, missingR8
   logical  :: has_fill_value      ! intended for future use
   logical  :: has_missing_value   ! intended for future use
   character(len=paramname_length) :: kind_string
   character(len=512) :: origin    ! the file it came from
   logical  :: update
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

! ---------------------------------------------------------------------
! The LatLon_to_dart structure is used when you know the 
! physical latitude and longitude gridcell indices and want to know 
! what indices of the DART state vector need to be queried.

type LatLon_to_dart
    private
    integer :: n
    integer, pointer, dimension(:) :: dartIndex
end type LatLon_to_dart
type(LatLon_to_dart), allocatable, dimension(:,:), target :: ij_to_dart

integer, allocatable, dimension(:,:) :: ij_to_land

! ---------------------------------------------------------------------
! Convert the land index to the 'parent' lat / lon gridcell
! This is like a sparse matrix representation.
! At present, since the dimension of y is forced to be 1 such that
! the size of x is equal to nland, ... the land_to_LatLon structure
! can be used with either 'ix' or 'iland' as the index 

type land_to_LatLon
    private
    integer :: latindex
    integer :: lonindex
end type land_to_LatLon
type(land_to_LatLon), allocatable, dimension(:), target :: land_to_ij

! -----------------------------------------------------------------------------
! ancillaries.nml &jules_frac      file='frac_lai_canht.nc' -or- 'static_input.nc'
! ancillaries.nml &jules_frac      frac_name='frac'
!
! model_grid.nml &jules_latlon     file='grid_lat_lon.nc' -or- 'static_input.nc'
!                                  lat_name='latitude'
!                                  lon_name='longitude'
!
! model_grid.nml &jules_land_frac  file='land_fraction.nc' -or- 'static_input.nc'
!                                  land_frac_name='land_mask'
!
! timesteps.nml  &jules_time       main_run_start='YYYY-MM-DD HH:MM:SS'
! timesteps.nml  &jules_time       main_run_end='YYYY-MM-DD HH:MM:SS'
! timesteps.nml  &jules_time       timestep_len=3600
!
! jules_surface_types.nml &jules_surface_types might be interesting ...
!
character(len=256) :: tile_fraction_file   = 'none'
character(len=256) :: land_fraction_file   = 'none'
character(len=256) :: physical_latlon_file = 'none'
character(len=NF90_MAX_NAME) :: tilefraction_variable = 'none'
character(len=NF90_MAX_NAME) :: landfraction_variable = 'none'
character(len=NF90_MAX_NAME) :: physical_lon_variable = 'none'
character(len=NF90_MAX_NAME) :: physical_lat_variable = 'none'

integer :: Nx_model = -1   ! output file, AKA 'x'
integer :: Ny_model = -1   ! output file, AKA 'y'
integer :: Nsoil    = -1   ! output, restart, jules_soil.nml
integer :: Ntile    = -1   ! output, restart
integer :: Nland    = -1   ! Number of gridcells containing land
integer :: Nscpool  = -1   ! Number of soil carbon pools
integer :: Npft     = -1   ! Number of plant types

real(r8), allocatable :: LONGITUDE(:,:)    ! output file, grid cell centers
real(r8), allocatable ::  LATITUDE(:,:)    ! output file, grid cell centers
real(r8), allocatable :: SOILLEVEL(:)      ! jules_soil.nml, soil interfaces
real(r4), allocatable :: land_tile_fractions(:,:,:)

type(location_type) :: ll_boundary   ! lower  left grid BOUNDARY
type(location_type) :: ur_boundary  ! upper right grid BOUNDARY

! The following variables are 'physically' shaped (i.e. you can plot them easily)

integer :: Nlon = -1   ! shape of physical (input) longitude/latitude matrix
integer :: Nlat = -1   ! shape of physical (input) longitude/latitude matrix

real(r8), allocatable :: physical_tile_fractions(:,:,:) ! land types (Nlon,Nlat,[Ntile,Ntype])
real(r8), allocatable :: land_mask(:,:)      ! gridcell fraction of land (Nlon,Nlat)
real(r8), allocatable :: physical_longitudes(:,:)
real(r8), allocatable :: physical_latitudes( :,:)

type(location_type), allocatable :: statespace_locations(:)
integer,             allocatable :: statespace_kinds(:)
type(get_close_type) :: gc_state

! -----------------------------------------------------------------
! module storage
! -----------------------------------------------------------------

integer         :: model_size      ! the state vector length
type(time_type) :: model_time      ! valid time of the model state
type(time_type) :: model_timestep  ! smallest time to adv model


INTERFACE vector_to_prog_var
      MODULE PROCEDURE vector_to_1d_prog_var
      MODULE PROCEDURE vector_to_2d_prog_var
      MODULE PROCEDURE vector_to_3d_prog_var
END INTERFACE

INTERFACE DART_get_var
      MODULE PROCEDURE get_var_1d
      MODULE PROCEDURE get_var_2d
      MODULE PROCEDURE get_var_3d
      MODULE PROCEDURE get_var_4d
END INTERFACE

INTERFACE get_state_time
      MODULE PROCEDURE get_state_time_ncid
      MODULE PROCEDURE get_state_time_fname
END INTERFACE


contains

!==================================================================
! All the REQUIRED interfaces come first - just by convention.
!==================================================================

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
!> Responsible for advancing JULES as a subroutine call.
!> since JULES is not subroutine callable, this is a stub.
!> any 'async' value other than 0 will result in an error.
!>
!> @param x the model state before and after the model advance.
!> @param time the desired time at the end of the model advance.
!>

subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

call error_handler(E_MSG, 'adv_1step:', 'FIXME RAFAEL routine not tested yet (needs long obs_seq file) ', source, revision, revdate)

if ( .not. module_initialized ) call static_init_model

if (do_output()) then
   call print_time(time,'NULL interface adv_1step (no advance) DART time is')
   call print_time(time,'NULL interface adv_1step (no advance) DART time is',logfileunit)
endif

write(string1,*)'DART should not be trying to advance JULES as a subroutine.'
write(string2,*)'DART can only run JULES as a shell script advance.'
call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate,text2=string2)

! just so suppress compiler warnings. code unreachable
x(:) = MISSING_R8

return
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
!> @param indx the index into the DART state vector
!> @param location the location at that index
!> @param var_type the DART KIND at that index
!>

subroutine get_state_meta_data(indx, location, var_type)

integer, intent(in)            :: indx
type(location_type)            :: location
integer, OPTIONAL, intent(out) :: var_type

! Local variables

integer  ::    lon_index   ! in the statespace framework, not physical space
integer  ::    lat_index   ! in the statespace framework, not physical space
integer  ::   soil_index
integer  ::   tile_index
integer  :: scpool_index
integer  :: vert_coord
integer  :: varindex
real(r8) :: mylon, mylat, mylev

if ( .not. module_initialized ) call static_init_model

varindex  = -1 ! if varindex is negative, get_state_indices will calculate it
call get_state_indices(indx, lon_index, lat_index, soil_index, &
                             tile_index, scpool_index, varindex)

call get_state_lonlatlev(varindex, lon_index, lat_index, soil_index, &
                             mylon, mylat, mylev, vert_coord)

location = set_location( mylon, mylat, mylev, vert_coord)

if (present(var_type)) var_type = progvar(varindex)%dart_kind

if (do_output() .and. debug > 5) then
   write(*,301) trim(progvar(varindex)%varname)
   write(*,300) trim(progvar(varindex)%varname), indx, &
      lon_index, lat_index, soil_index, tile_index, scpool_index, &
      mylon, mylat, mylev, vert_coord
endif

 301 format(A,' DART_index  x_ind  y_ind levind   tile  cpool   longitude    latitude       level coord          value')
 300 format(A,1x,i10,5(1x,i6),3(1x,f11.6),1x,i3)

return
end subroutine get_state_meta_data


!-----------------------------------------------------------------------
!>
!> The basis for all 'forward observation operators'.
!> For a given lat, lon, and height, interpolate the correct state value
!> to that location for the filter from the JULES state vectors
!>
!> The type of the variable being interpolated is obs_type since 
!> normally this is used to find the expected value of an observation 
!> at some location. The interpolated value is returned in interp_val 
!> and istatus is 0 for success. NOTE: This is a workhorse routine and is
!> the basis for all the forward observation operator code.
!>
!> Reconstructing the vertical profile of the gridcell is complicated.
!> Each land unit/column can have a different number of vertical levels.
!> Do I just try to recontruct the whole profile and then interpolate?
!> Impossible to know which elements are 'above' and 'below' without
!> finding all the elements in the first place. The vertical information
!> is in the levels() array for each state vector component.
!>
!> @param x the DART state vector
!> @param location the location of interest
!> @param obs_kind the DART KIND of interest
!> @param interp_val the estimated value of the DART state at the location 
!>          of interest (the interpolated value).
!> @param istatus interpolation status ... 0 == success, /=0 is a failure
!>
!> @todo FIXME use some unique error code if the location is technically
!> outside the domain, i.e. an extrapolation. At some point it will be
!> useful to know if the interpolation failed because of some illegal
!> state as opposed to simply being outside the domain.

subroutine model_interpolate(x, location, obs_kind, interp_val, istatus)

real(r8),            intent(in)  :: x(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_kind
real(r8),            intent(out) :: interp_val
integer,             intent(out) :: istatus

! Local storage

real(r8), dimension(LocationDims) :: loc_array
real(r8) :: llon, llat, lheight

character(len=obstypelength) :: kind_name

integer  :: varindex, testindex, itile
integer  :: i, x_index, y_index, indx, vert_coord
integer  :: lon_index, lat_index, soil_index, tile_index, scpool_index
integer  :: closest_index(1)
integer  :: num_close, ncontrib
integer  :: close_ind(Nlon*Nlat)
real(r8) :: distances(Nlon*Nlat)
real(r8) :: mindistance
real(r8) :: mylon, mylat, mylev
real(r8), allocatable :: profile(:)

integer  :: top, bot
real(r8) :: slope, xtrcp

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

interp_val = MISSING_R8  ! the DART bad value flag
istatus    = 99          ! unknown error

! If identity observation (obs_kind < 0), then no need to interpolate
! identity observation -> -(obs_kind)=DART state vector index
! obtain state value directly from index
! obs_kind == 0 == RAW_STATE_VARIABLE

if ( obs_kind < 0 ) then

   if ((-1*obs_kind > 1) .and. (-1*obs_kind <= model_size)) then
      interp_val = x(-1*obs_kind)
      istatus = 0
   else
      istatus = 98
   endif

   return

endif

! Get the individual locations values - mostly for messages

loc_array = get_location(location)
llon      = loc_array(1)
llat      = loc_array(2)
lheight   = loc_array(3)
kind_name = get_raw_obs_kind_name(obs_kind)

! Check to make sure we can interpolate this kind before we go through any more work.
varindex = 0
KindCheck : do i = 1,nfields
   if (progvar(i)%kind_string == kind_name) then
      varindex = i
      exit KindCheck
   endif
enddo KindCheck

if (varindex == 0) then
   istatus = 97
   return
endif

! Check to ensure the observation is within the physical domain.
! Prevents against extrapolation.
!> @ todo FIXME return some specific error code that will indicate it 
!>        is outside the domain ... See JIRA ticket DARTSUP-294

if( .not. is_location_in_region(location, ll_boundary, ur_boundary)) then

   if (do_output() .and. debug > 2) then
      write(string1,*)'Requesting interpolation for ',trim(kind_name), ' at'
      call write_location(istatus,location,charstring=string2)
      write(string3,*)'which is outside the domain.'
      call error_handler(E_MSG, 'model_interpolate:', string1, &
                                text2=string2, text3=string3)
   endif
   
   istatus = 96
   return
endif

! Generate the list of physical gridcell indices that are 'close' to the observation.

call loc_get_close_obs(gc_state, location, 1, statespace_locations, statespace_kinds, &
                       num_close, close_ind, distances)

if (num_close == 0) then
   call write_location(num_close, location, charstring=string1)
   call error_handler(E_MSG, 'model_interpolate:', string1, text2='nothing close')
   istatus = 95
   return
endif

mindistance   =           minval(distances(1:num_close))
closest_index = close_ind(minloc(distances(1:num_close)))

! We now have the closest index [1, Nlon*Nlat] (i.e. statespace format) and need
! to decompose it into the lat,lon indices of the JULES physical grid.

y_index = 1 + (closest_index(1) - 1)/Nlon
x_index = closest_index(1) - (y_index - 1)*Nlon

if (do_output() .and. debug > 3) then

   write(string1,*)'Requesting interpolation for ',trim(kind_name), ' at'
   call write_location(istatus,location,charstring=string2)
   write(string3,*)'There are ',num_close,' close "sparse" grid locations.'
   call error_handler(E_MSG, 'model_interpolate:', string1, text2=string2, text3=string3)

   write(*,*)'distances       are ',distances(1:num_close)
   write(*,*)'indices         are ',close_ind(1:num_close)
   write(*,*)'minimum distance is ',mindistance
   write(*,*)'statespace index is ',closest_index
   write(*,*)'physical lon  index ',x_index
   write(*,*)'physical lat  index ',y_index

endif

! Closest gridcell is not land.

if( land_mask(x_index, y_index) < 1.0_r8 ) then

   if (do_output() .and. debug > 2) then
      write(string1,*)'Requesting interpolation for ',trim(kind_name), ' at'
      call write_location(istatus,location,charstring=string2)
      write(string3,*)'which is apparently not a land gridcell.'
      call error_handler(E_MSG, 'model_interpolate:', string1, &
                                text2=string2, text3=string3)
   endif

   istatus = 94
   return
endif

! structure to return the DART elements that are at that location.
! This block is just a check.

if (do_output() .and. debug > 2) then

   write(*,*) ij_to_dart(x_index,y_index)%n, &
              'DART elements from parent gridcell ', x_index, y_index
   write(*,301) trim(progvar(varindex)%varname)
   
   do i = 1, ij_to_dart(x_index,y_index)%n
      indx = ij_to_dart(x_index,y_index)%dartIndex(i)
   
      testindex = -1   
      call get_state_indices(indx, lon_index, lat_index, soil_index, &
                                tile_index, scpool_index, testindex)
   
      call get_state_lonlatlev(testindex, lon_index, lat_index, soil_index, &
                                mylon, mylat, mylev, vert_coord)

      if ( testindex == varindex ) then ! this is an item of interest.
         write(*,300) trim(progvar(varindex)%varname), indx, &
            lon_index, lat_index, soil_index, tile_index, scpool_index, &
            mylon, mylat, mylev, vert_coord, x(indx)
      endif
   enddo
endif

 301 format(A,' DART_index  x_ind  y_ind levind   tile  cpool   longitude    latitude       level coord          value')
 300 format(A,1x,i10,5(1x,i6),3(1x,f11.6),1x,i3,1x,f16.7)

! some variables are (land,tile), some are (land,soil), (land,scpool)
! so some results are simply a scalar (land,tile)
! some results are vectors (land,soil), (land,scpool) and may need more 
!
! Since get_state_indices works with any shape variable,
! we can treat the 'x y tile time' case the same as 'land tile'.
! Similarly for other cases.
!
! At this point:
! x_index  is for the ij_to_dart structure
! y_index  is for the ij_to_dart structure
! varindex is the DART variable we need to interpolate

select case (trim(progvar(varindex)%coordinates))

case ('land tile','x y tile time') ! There is no vertical dimension ... 

   if (debug > 2) then
      do itile = 1,Ntile
         write(string1,*)'physical x_i,y_i,tile_fraction ',x_index,y_index,itile, &
                physical_tile_fractions(x_index,y_index,itile)
         call error_handler(E_MSG,'model_interpolate:',string1)
      enddo
   endif

   interp_val = 0.0_r8

   ! Create the weighted average of the gridcell value based on
   ! the tile fractions of that particular gridcell 

   do i = 1, ij_to_dart(x_index,y_index)%n
      indx = ij_to_dart(x_index,y_index)%dartIndex(i)
   
      testindex = -1   
      call get_state_indices(indx, lon_index, lat_index, soil_index, &
                                tile_index, scpool_index, testindex)

      if (testindex == varindex ) then

         interp_val = interp_val + &
                      physical_tile_fractions(x_index,y_index,tile_index) * x(indx)

         if (debug > 2) then
            write(string1,*)'tile_fraction * x += weighted_sum ', & 
            physical_tile_fractions(x_index,y_index,tile_index), x(indx), interp_val
            call error_handler(E_MSG,'model_interpolate:',string1)
         endif

      endif
   enddo

   istatus = 0

case ('land soil','x y soil time')
   ! must vertically interpolate to proper depth
   allocate( profile(Nsoil) )
   profile = 0.0_r8

   ncontrib = 0
WANTED : do i = 1, ij_to_dart(x_index,y_index)%n
      indx = ij_to_dart(x_index,y_index)%dartIndex(i)

      if ( (indx < progvar(varindex)%index1) .or. &
           (indx > progvar(varindex)%indexN) ) cycle WANTED

      call get_state_indices(indx, lon_index, lat_index, soil_index, &
                             tile_index, scpool_index, varindex)

      profile(soil_index) = x(indx)

      ncontrib = ncontrib + 1
   
   enddo WANTED

   if (ncontrib /= Nsoil) then
      write(string1,*)trim(progvar(varindex)%varname),'has coordinates', &
                      trim(progvar(varindex)%coordinates)
      write(string2,*)'but had ',ncontrib,' components instead of ',Nsoil  
      write(string3,*)'Not right - stopping.'
      call error_handler(E_ERR, 'model_interpolate', string1, &
                 source, revision, revdate, text2=string2, text3=string3)
   endif

   ! So now we have the soil profile at the location of interest.
   ! Do the vertical interpolation and we're done.

   if ( lheight <= SOILLEVEL(1) ) then
      interp_val = profile(1)
      istatus    = 0
      return
   endif

   !> @todo FIXME ... may want to put a depth limit on this ...
   if ( lheight >= SOILLEVEL(Nsoil) ) then
      interp_val = profile(Nsoil)
      istatus    = 0
      return
   endif

   ! OK - so we're somewhere in the middle ... do a simple
   ! linear interpolation
   Depth : do i = 2,Nsoil
      if ( lheight < SOILLEVEL(i) ) then
         top   = i-1
         bot   = i
         slope = (profile(bot) - profile(top)) / (SOILLEVEL(bot) - SOILLEVEL(top))
         xtrcp =  profile(bot) - slope*SOILLEVEL(bot)
         interp_val = slope * lheight + xtrcp
         istatus = 0
         return
      endif
   enddo Depth

   deallocate(profile)

case ('land scpool','x y scpool time')

   !> @todo FIXME don't know what to do here ... what is scpool
   !> if scpool is a singleton dimension - I suspect we can just use it
   !> the problem comes up when the scpool dimension is > 1

   write(string1,*)'variable ',trim(progvar(varindex)%varname)
   write(string2,*)'has dimensions <'//trim( progvar(varindex)%coordinates)//'>'
   write(string3,*)'which contains "scpool" - Unable to proceed - stopping.'
   call error_handler(E_ERR, 'model_interpolate', string1, &
              source, revision, revdate, text2=string2, text3=string3)

   ! if the dimension of scpool is just '1', I think the simple case('land')
   ! has the same logic that we want ...

case ('land','x y')

   do i = 1, ij_to_dart(x_index,y_index)%n
      indx = ij_to_dart(x_index,y_index)%dartIndex(i)
   
      testindex = -1   
      call get_state_indices(indx, lon_index, lat_index, soil_index, &
                                tile_index, scpool_index, testindex)

      if (testindex == varindex ) then

         if (debug > 2) then
            write(string1,*)'only land dimension', & 
            x_index, y_index, indx, x(indx), interp_val
            call error_handler(E_MSG,'model_interpolate:',string1)
         endif

         if (interp_val /= MISSING_R8) then
            write(string1,*)'more than 1 value for ',trim(progvar(varindex)%varname)
            write(string2,*)'at the desired location. Unexpected. Stopping.'
            call error_handler(E_ERR, 'model_interpolate', string1, &
                       source, revision, revdate, text2=string2)
         else
            interp_val = x(indx)
         endif

      endif
   enddo

   istatus = 0

case default

   write(string1,*)'only supports specific interpolations.'
   write(string2,*)trim(progvar(varindex)%varname),' has coordinates "', &
                   trim(progvar(varindex)%coordinates),'"'
   write(string3,*)'Unable to proceed - stopping.'
   call error_handler(E_ERR, 'model_interpolate', string1, &
              source, revision, revdate, text2=string2, text3=string3)
end select

return
end subroutine model_interpolate


!-----------------------------------------------------------------------
!>
!> Returns the the time step of the model; the smallest increment in 
!> time that the model is capable of advancing the JULES state.
!>

function get_model_time_step()

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = model_timestep

end function get_model_time_step


!-----------------------------------------------------------------------
!>
!> Called to do one-time initialization of the model.
!> In this case, it reads in the grid information, the namelist
!> containing the variables of interest, where to get them, their size, 
!> their associated DART KIND, etc.
!>
!> In addition to harvesting the model metadata (grid,
!> desired model advance step, etc.), it also fills a structure
!> containing information about what variables are where in the DART
!> framework.
!>

subroutine static_init_model()

! Local variables - all the important ones have module scope

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: dimname, coordinates
integer :: ncid, TimeDimID, VarID, dimlen, varsize
integer :: iunit, io, ivar
integer :: i, index1, indexN
integer :: ss, dd

integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: spvalR8

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
if (do_output()) call error_handler(E_MSG,'static_init_model:','model_nml values are')
if (do_output()) write(logfileunit, nml=model_nml)
if (do_output() .and. debug > 3) write(     *     , nml=model_nml)

! Set the time step ... causes JULES namelists to be read.
! Ensures model_timestep is multiple of 'dynamics_timestep'

call set_calendar_type( calendar )   ! comes from model_mod_nml

model_time     = get_state_time(jules_output_filename)
model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model:',string1)

! The 'Input' grids are always physically-based, the 'Model' grid
! may have all the non-land gridcells squeezed out of it.
! The JULES output file has the model grid metadata.
! The JULES restart files are intentionally lean and, in so doing,
! do not have metadata.
! The JULES grid in an output file can take on two forms.
! If there are no masked (non-land) gridcells, it can be a regular 2D grid.
! If a mask has been applied, the number of longitudes takes on the number
! of useful land gridcells and the number of latitudes is 1.
! DART only supports the 'masked' version, i.e. the 1D version.

call get_jules_output_dimensions()
call get_jules_restart_dimensions()

! need to know the filenames for the physical grid information

call get_physical_filenames()   ! from jules namelist files
call set_physical_grid(physical_latlon_file)
call set_grid_boundary()
call set_land_mask(land_fraction_file)
call set_tile_fractions(tile_fraction_file)

call get_soil_levels()

call Get_Model_Grid(jules_output_filename)

! Compile the list of JULES variables to use in the creation
! of the DART state vector. This just checks to see that the
! DART KIND is valid and that the variable was specified correctly.
! Whether or not the JULES variables exist is checked in the
! rather long loop that follows.

nfields = parse_variable_table()

! Compute the offsets into the state vector for the start of each
! variable type. Requires reading shapes from the JULES restart file.
! Record the extent of the data type in the state vector.

index1  = 1;
indexN  = 0;
do ivar = 1, nfields

   ! convey the information in the variable_table to each progvar.

   call SetVariableAttributes(ivar)

   ! Open the file for each variable and get dimensions, etc.
   ! Since the JULES restart files have no metadata (!?)
   ! we must get all the per-variable metadata from the
   ! (hopefully identical) variables in the output file.

   call nc_check(nf90_open(trim(progvar(ivar)%origin), NF90_NOWRITE, ncid), &
              'static_init_model','open '//trim(progvar(ivar)%origin))

   ! File is not required to have a time dimension
   io = nf90_inq_dimid(ncid, 'time', TimeDimID)
   if (io /= NF90_NOERR) TimeDimID = MISSING_I

   string2 = trim(progvar(ivar)%origin)//' '//trim(progvar(ivar)%varname)

   call nc_check(nf90_inq_varid(ncid, trim(progvar(ivar)%varname), VarID), &
            'static_init_model', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, &
                 ndims=progvar(ivar)%numdims, xtype=progvar(ivar)%xtype), &
            'static_init_model', 'inquire '//trim(string2))

   ! by default, the rank and the number of dimensions are identical.
   ! Since DART treats 'time' as a singleton dimension - the rank is one less
   ! if one of the dimensions is 'time'.

   progvar(ivar)%rank = progvar(ivar)%numdims

   ! If the long_name and/or units attributes are set, get them.
   ! They are not REQUIRED to exist but are nice to use if they are present.
   !> @todo FIXME ... we could check the same variable in the output file
   ! to get the metadata

   if( nf90_inquire_attribute(    ncid, VarID, 'long_name') == NF90_NOERR ) then
      call nc_check( nf90_get_att(ncid, VarID, 'long_name' , progvar(ivar)%long_name), &
                  'static_init_model', 'get_att long_name '//trim(string2))
   else
      progvar(ivar)%long_name = progvar(ivar)%varname
   endif

   if( nf90_inquire_attribute(    ncid, VarID, 'units') == NF90_NOERR )  then
      call nc_check( nf90_get_att(ncid, VarID, 'units' , progvar(ivar)%units), &
                  'static_init_model', 'get_att units '//trim(string2))
   else
      progvar(ivar)%units = 'unknown'
   endif

   ! Save any FillValue, missing_value attributes for reading, writing ...

   if (progvar(ivar)%xtype == NF90_INT) then
       if (nf90_get_att(ncid, VarID, '_FillValue'    , spvalINT) == NF90_NOERR) then
          progvar(ivar)%spvalINT       = spvalINT
          progvar(ivar)%has_fill_value = .true.
       endif
       if (nf90_get_att(ncid, VarID, 'missing_value' , spvalINT) == NF90_NOERR) then
          progvar(ivar)%missingINT        = spvalINT
          progvar(ivar)%has_missing_value = .true.
       endif

   elseif (progvar(ivar)%xtype == NF90_FLOAT) then
       if (nf90_get_att(ncid, VarID, '_FillValue'    , spvalR4) == NF90_NOERR) then
          progvar(ivar)%spvalR4        = spvalR4
          progvar(ivar)%has_fill_value = .true.
       endif
       if (nf90_get_att(ncid, VarID, 'missing_value' , spvalR4) == NF90_NOERR) then
          progvar(ivar)%missingR4         = spvalR4
          progvar(ivar)%has_missing_value = .true.
       endif

   elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
       if (nf90_get_att(ncid, VarID, '_FillValue'    , spvalR8) == NF90_NOERR) then
          progvar(ivar)%spvalR8        = spvalR8
          progvar(ivar)%has_fill_value = .true.
       endif
       if (nf90_get_att(ncid, VarID, 'missing_value' , spvalR8) == NF90_NOERR) then
          progvar(ivar)%missingR8         = spvalR8
          progvar(ivar)%has_missing_value = .true.
       endif
   endif

   varsize = 1
   dimlen  = 1
   coordinates = ''
   DimensionLoop : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), name=dimname, len=dimlen), &
                                          'static_init_model', string1)

      ! The 'land' dimension of this variable must match the number
      ! of elements in the LONGITUDE, LATITUDE arrays.
      ! The 'x' and 'y' dimensions in the output file already match
      ! so there is no need to check those.
      if ( (dimname == 'land') .and. (dimlen /= Nx_model*Ny_model) )then
         write(string1,*)'dimension mismatch between restart and output'
         write(string2,*)trim(progvar(ivar)%varname),' land dimension is ',dimlen
         write(string3,*)'number of active land cells from output is ',Nx_model*Ny_model
         call error_handler(E_ERR, 'static_init_model', string1, &
                    source, revision, revdate, text2=string2, text3=string3)
      endif

      ! The 'soil' dimension of this variable must match the number
      ! of soil layers from the namelist
      if ( (dimname == 'soil') .and. (dimlen /= Nsoil) )then
         write(string1,*)'dimension mismatch between file and namelist'
         write(string2,*)trim(progvar(ivar)%varname),' soil dimension is ',dimlen
         write(string3,*)'number of soil layers from namelist is ',Nsoil
         call error_handler(E_ERR, 'static_init_model', string1, &
                    source, revision, revdate, text2=string2, text3=string3)
      endif

      ! Only reserve space for a single time slice
      if (dimIDs(i) == TimeDimID) then
         dimlen = 1
         progvar(ivar)%rank = progvar(ivar)%numdims - 1
      endif

      progvar(ivar)%dimlens( i) = dimlen
      progvar(ivar)%dimnames(i) = dimname
      varsize = varsize * dimlen

      ! For some purposes, it is nice to have a little coordinate array of all the 
      ! dimnames pasted together with a minimal amount of whitespace separation.

      write(coordinates,*) trim(coordinates)//' '//trim(dimname)

   enddo DimensionLoop

   progvar(ivar)%coordinates = adjustl(coordinates)
   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1
   index1                    = index1 + varsize      ! sets up for next variable

   if (do_output() .and. debug > 0) call dump_structure(ivar)

   call nc_check(nf90_close(ncid),'static_init_model','close '//trim(string2))
   ncid = 0

enddo

model_size = progvar(nfields)%indexN

if (do_output() .and. debug > 99) then
  write(logfileunit, *)
  write(logfileunit,'("grid: Nx_model, Ny_model, Nsoil =",3(1x,i6))') &
                             Nx_model, Ny_model, Nsoil
  write(logfileunit, *)'model_size = ', model_size
  write(     *     , *)
  write(     *     ,'("grid: Nx_model, Ny_model, Nsoil =",3(1x,i6))') &
                             Nx_model, Ny_model, Nsoil
  write(     *     , *)'model_size = ', model_size
endif

! Generate list of dart indices for each gridcell

allocate(ij_to_dart(Nlon,Nlat))
call determine_parent_gridcells()

! Generate speedup table of 'close' items to figure out what gridcell
! is close to some arbitrary (observation) location.

call set_sparse_locations_kinds()
call init_interp()

return
end subroutine static_init_model


!-----------------------------------------------------------------------
!>
!> Does any shutdown and clean-up needed for model.
!>

subroutine end_model()

   if (allocated(LATITUDE)               ) deallocate(LATITUDE) 
   if (allocated(LONGITUDE)              ) deallocate(LONGITUDE) 
   if (allocated(SOILLEVEL)              ) deallocate(SOILLEVEL)
   if (allocated(land_tile_fractions)    ) deallocate(land_tile_fractions)

   if (allocated(land_mask)              ) deallocate(land_mask)
   if (allocated(physical_longitudes)    ) deallocate(physical_longitudes) 
   if (allocated(physical_latitudes)     ) deallocate(physical_latitudes) 
   if (allocated(physical_tile_fractions)) deallocate(physical_tile_fractions) 

   if (allocated(statespace_locations)   ) deallocate(statespace_locations)
   if (allocated(statespace_kinds)       ) deallocate(statespace_kinds)

   if (allocated(land_to_ij)             ) deallocate(land_to_ij)
   if (allocated(ij_to_land)             ) deallocate(ij_to_land)

return
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
!> NOTE: Since JULES cannot start in this manner,
!> DART will intentionally generate a fatal error.
!>
!> @param time the time to associate with the initial state
!>

subroutine init_time(time)

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'No good way to specify an arbitrary initial time for JULES.'
call error_handler(E_ERR,'init_time',string1,source,revision,revdate)

! Just set to 0 to silence the compiler warnings.
time = set_time(0,0)

return
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
!> NOTE: This is not supported for JULES and will generate a FATAL ERROR. 
!>       However, this is a required interface - so it must be present.
!>
!> @param x the JULES initial conditions  

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'No good way to specify an arbitrary initial conditions for JULES.'
call error_handler(E_ERR,'init_conditions',string1,source,revision,revdate)

! Just set to 0 to silence the compiler warnings.
x = 0.0_r8

return
end subroutine init_conditions


!-----------------------------------------------------------------------
!>
!> Writes the model-specific attributes to a DART 'diagnostic' netCDF file.
!> This includes coordinate variables and some metadata, but NOT the
!> actual DART state. That may be done multiple times in nc_write_model_vars()
!>
!> @param ncFileID the netCDF handle of the DART diagnostic file opened by
!>                 assim_model_mod:init_diag_output
!> @param ierr status ... 0 == all went well, /= 0 failure
!>

function nc_write_model_atts( ncFileID ) result (ierr)

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

! variables if we just blast out one long state vector

integer :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)

! variables if we parse the state vector into prognostic variables.

! for the dimensions and coordinate variables
integer ::   NlonDimID   ! Physical grid longitude dimension
integer ::   NlatDimID   ! Physical grid latitude  dimension
integer ::   NxDimID     ! model grid dimension
integer ::   NyDimID     ! model grid dimension
integer ::  nsoilDimID
integer ::  nlandDimID   ! number of land grid cells 
integer ::  ntileDimID   ! number of tiles
integer :: NscpoolDimID

! for the prognostic variables
integer :: ivar, VarID

! local variables

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1
character(len=NF90_MAX_NAME) :: varname
integer, dimension(NF90_MAX_VAR_DIMS) :: mydimids
integer :: myndims

character(len=128) :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

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
                          'nc_write_model_atts', 'inq_dimid copy '//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid=  TimeDimID), &
                          'nc_write_model_atts', 'inq_dimid time '//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)'Time Dimension ID ',TimeDimID, &
             ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', string1, source, revision, revdate)
endif

! Define the model size / state variable dimension / whatever ...
call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable', len=model_size, &
        dimid = StateVarDimID),'nc_write_model_atts', &
       'def_dim StateVariable '//trim(filename))

! Write Global Attributes

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'creation_date' ,str1    ), &
           'nc_write_model_atts', 'put_att creation '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_source'  ,source  ), &
           'nc_write_model_atts', 'put_att source '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revision',revision), &
           'nc_write_model_atts', 'put_att revision '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revdate' ,revdate ), &
           'nc_write_model_atts', 'put_att revdate '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'JULES' ), &
           'nc_write_model_atts', 'put_att model '//trim(filename))

! We need to output the prognostic variables.
! Define the new dimensions IDs

call nc_check(nf90_def_dim(ncid=ncFileID, name='x' , len = Nx_model, &
          dimid=NxDimID),   'nc_write_model_atts', 'def_dim x '//trim(filename))

call nc_check(nf90_def_dim(ncid=ncFileID, name='y' , len = Ny_model, &
          dimid=NyDimID),   'nc_write_model_atts', 'def_dim y '//trim(filename))

call nc_check(nf90_def_dim(ncid=ncFileID, name='soil', len = Nsoil, &
          dimid=nsoilDimID),  'nc_write_model_atts', 'def_dim soil '//trim(filename))

call nc_check(nf90_def_dim(ncid=ncFileID, name='land', len = Nland, &
          dimid=nlandDimID),  'nc_write_model_atts', 'def_dim land '//trim(filename))

call nc_check(nf90_def_dim(ncid=ncFileID, name='tile', len = Ntile, &
          dimid=ntileDimID),  'nc_write_model_atts', 'def_dim tile '//trim(filename))

call nc_check(nf90_def_dim(ncid=ncFileID, name='scpool', len = Nscpool, &
          dimid=NscpoolDimID),'nc_write_model_atts', 'def_dim scpool '//trim(filename))

call nc_check(nf90_def_dim(ncid=ncFileID, name='longitude', len = Nlon, &
          dimid=NlonDimID),'nc_write_model_atts', 'def_dim longitude '//trim(filename))

call nc_check(nf90_def_dim(ncid=ncFileID, name='latitude', len = Nlat, &
          dimid=NlatDimID),'nc_write_model_atts', 'def_dim latitude '//trim(filename))

! Create the (empty) Coordinate Variables and the Attributes

! Model Grid Longitudes
call nc_check(nf90_def_var(ncFileID,name='x', xtype=nf90_real, &
              dimids=(/ NxDimID, NyDimID /), varid=VarID),&
              'nc_write_model_atts', 'def_var x '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'land-only gridcell longitudes'), &
              'nc_write_model_atts', 'put_att longitude long_name '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'X'),  &
              'nc_write_model_atts', 'put_att longitude cartesian_axis '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
              'nc_write_model_atts', 'put_att longitude units '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ 0.0_r4, 360.0_r4 /)), &
              'nc_write_model_atts', 'put_att longitude valid_range '//trim(filename))

! Model Grid Latitudes
call nc_check(nf90_def_var(ncFileID,name='y', xtype=nf90_real, &
              dimids=(/ NxDimID, NyDimID /), varid=VarID),&
              'nc_write_model_atts', 'def_var y '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'land-only gridcell latitudes'), &
              'nc_write_model_atts', 'put_att latitude long_name '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Y'),   &
              'nc_write_model_atts', 'put_att latitude cartesian_axis '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_north'),  &
              'nc_write_model_atts', 'put_att latitude units '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID,'valid_range',(/ -90.0_r4, 90.0_r4 /)), &
              'nc_write_model_atts', 'put_att latitude valid_range '//trim(filename))

! subsurface levels
call nc_check(nf90_def_var(ncFileID,name='soil', xtype=nf90_real, &
              dimids=(/ nsoilDimID /), varid=VarID),&
              'nc_write_model_atts', 'def_var soil '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'coordinate soil levels'), &
              'nc_write_model_atts', 'put_att soil long_name '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Z'),   &
              'nc_write_model_atts', 'put_att soil cartesian_axis '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'm'),  &
              'nc_write_model_atts', 'put_att soil units '//trim(filename))

! Model Grid tile fractions
call nc_check(nf90_def_var(ncFileID,name='land_tile_fractions', xtype=nf90_real, &
              dimids=(/ NxDimID, NyDimID, nTileDimID /), varid=VarID),&
              'nc_write_model_atts', 'def_var land_tile_fractions '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'land-only gridcell tile fractions'), &
              'nc_write_model_atts', 'put_att land_tile_fractions long_name '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'units', ' '),  &
              'nc_write_model_atts', 'put_att land_tile_fractions units '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID,'valid_range',(/ 0.0_r4, 1.0_r4 /)), &
              'nc_write_model_atts', 'put_att land_tile_fractions valid_range '//trim(filename))

! Physical Grid Longitudes
call nc_check(nf90_def_var(ncFileID,name='longitude', xtype=nf90_real, &
              dimids=(/ NlonDimID, NlatDimID /), varid=VarID),&
              'nc_write_model_atts', 'def_var longitude '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'physical-space longitudes'), &
              'nc_write_model_atts', 'put_att longitude long_name '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'X'),  &
              'nc_write_model_atts', 'put_att longitude cartesian_axis '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
              'nc_write_model_atts', 'put_att longitude units '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ 0.0_r4, 360.0_r4 /)), &
              'nc_write_model_atts', 'put_att longitude valid_range '//trim(filename))

! Physical Grid Latitudes
call nc_check(nf90_def_var(ncFileID,name='latitude', xtype=nf90_real, &
              dimids=(/ NlonDimID, NlatDimID /), varid=VarID),&
              'nc_write_model_atts', 'def_var latitude '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'physical-space latitudes'), &
              'nc_write_model_atts', 'put_att latitude long_name '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Y'),   &
              'nc_write_model_atts', 'put_att latitude cartesian_axis '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_north'),  &
              'nc_write_model_atts', 'put_att latitude units '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID,'valid_range',(/ -90.0_r4, 90.0_r4 /)), &
              'nc_write_model_atts', 'put_att latitude valid_range '//trim(filename))

! Physical Grid tile fractions
call nc_check(nf90_def_var(ncFileID,name='physical_tile_fractions', xtype=nf90_real, &
              dimids=(/ NlonDimID, NlatDimID, nTileDimID /), varid=VarID),&
              'nc_write_model_atts', 'def_var physical_tile_fractions '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'physical-space tile fractions'), &
              'nc_write_model_atts', 'put_att physical_tile_fractions long_name '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'units', ' '),  &
              'nc_write_model_atts', 'put_att physical_tile_fractions units '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID,'valid_range',(/ 0.0_r4, 1.0_r4 /)), &
              'nc_write_model_atts', 'put_att physical_tile_fractions valid_range '//trim(filename))

!> @todo FIXME ... do we need a model grid tile fractions ... and if so, what
! dimensions should it use ... nland -or- nx and ny ... what variables use it?

! Physical Grid land fractions
call nc_check(nf90_def_var(ncFileID,name='land_mask', xtype=nf90_real, &
              dimids=(/ NlonDimID, NlatDimID /), varid=VarID),&
              'nc_write_model_atts', 'def_var land_mask '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'physical-space land fractions'), &
              'nc_write_model_atts', 'put_att land_mask long_name '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'units', ' '),  &
              'nc_write_model_atts', 'put_att land_mask units '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID,'valid_range',(/ 0.0_r4, 1.0_r4 /)), &
              'nc_write_model_atts', 'put_att land_mask valid_range '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'comment', '1==land, 0==not land'),  &
              'nc_write_model_atts', 'put_att land_mask units '//trim(filename))

! Create the (empty) Prognostic Variables and the Attributes

do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   string1 = trim(filename)//' '//trim(varname)

   ! match shape of the variable to the dimension IDs

   call define_var_dims(ivar, ncFileID, MemberDimID, unlimitedDimID, myndims, mydimids)

   ! define the variable and set the attributes

   call nc_check(nf90_def_var(ncid=ncFileID, name=trim(varname), &
                 xtype=progvar(ivar)%xtype, &
                 dimids = mydimids(1:myndims), varid=VarID),&
                 'nc_write_model_atts', 'def_var '//trim(string1))

   call nc_check(nf90_put_att(ncFileID, VarID, &
           'long_name', trim(progvar(ivar)%long_name)), &
           'nc_write_model_atts', 'put_att long_name '//trim(string1))
   call nc_check(nf90_put_att(ncFileID, VarID, &
           'DART_kind', trim(progvar(ivar)%kind_string)), &
           'nc_write_model_atts', 'put_att DART_kind '//trim(string1))
   call nc_check(nf90_put_att(ncFileID, VarID, &
           'units', trim(progvar(ivar)%units)), &
           'nc_write_model_atts', 'put_att units '//trim(string1))

   ! Preserve the original missing_value/_FillValue code.

   if (  progvar(ivar)%xtype == NF90_INT ) then
      call nc_check(nf90_put_att(ncFileID, VarID, &
              'missing_value', progvar(ivar)%spvalINT), &
              'nc_write_model_atts', 'put_att missing_value '//trim(string1))
      call nc_check(nf90_put_att(ncFileID, VarID, &
              '_FillValue',  progvar(ivar)%spvalINT), &
              'nc_write_model_atts', 'put_att _FillValue '//trim(string1))

   elseif (  progvar(ivar)%xtype == NF90_FLOAT ) then
      call nc_check(nf90_put_att(ncFileID, VarID, &
              'missing_value', progvar(ivar)%spvalR4), &
              'nc_write_model_atts', 'put_att missing_value '//trim(string1))
      call nc_check(nf90_put_att(ncFileID, VarID, &
              '_FillValue',  progvar(ivar)%spvalR4), &
              'nc_write_model_atts', 'put_att _FillValue '//trim(string1))

   elseif (  progvar(ivar)%xtype == NF90_DOUBLE ) then
      call nc_check(nf90_put_att(ncFileID, VarID, &
              'missing_value', progvar(ivar)%spvalR8), &
              'nc_write_model_atts', 'put_att missing_value '//trim(string1))
      call nc_check(nf90_put_att(ncFileID, VarID, &
              '_FillValue',  progvar(ivar)%spvalR8), &
              'nc_write_model_atts', 'put_att _FillValue '//trim(string1))
   endif

enddo

! Finished with dimension/variable definitions, must end 'define' mode to fill.

call nc_check(nf90_enddef(ncfileID), 'prognostic enddef '//trim(filename))

! Fill the coordinate variables

call nc_check(nf90_inq_varid(ncFileID, 'x', VarID), &
             'nc_write_model_atts', 'inq_varid longitude '//trim(filename))
call nc_check(nf90_put_var(ncFileID, VarID, LONGITUDE ), &
             'nc_write_model_atts', 'put_var longitude '//trim(filename))

call nc_check(nf90_inq_varid(ncFileID, 'y', VarID), &
             'nc_write_model_atts', 'inq_varid latitude '//trim(filename))
call nc_check(nf90_put_var(ncFileID, VarID, LATITUDE ), &
             'nc_write_model_atts', 'put_var latitude '//trim(filename))

call nc_check(nf90_inq_varid(ncFileID, 'soil', VarID), &
             'nc_write_model_atts', 'inq_varid soil '//trim(filename))
call nc_check(nf90_put_var(ncFileID, VarID, SOILLEVEL ), &
             'nc_write_model_atts', 'put_var soil '//trim(filename))

call nc_check(nf90_inq_varid(ncFileID, 'longitude', VarID), &
             'nc_write_model_atts', 'inq_varid longitude '//trim(filename))
call nc_check(nf90_put_var(ncFileID, VarID, physical_longitudes ), &
             'nc_write_model_atts', 'put_var longitude '//trim(filename))

call nc_check(nf90_inq_varid(ncFileID, 'latitude', VarID), &
             'nc_write_model_atts', 'inq_varid latitude '//trim(filename))
call nc_check(nf90_put_var(ncFileID, VarID, physical_latitudes ), &
             'nc_write_model_atts', 'put_var latitude '//trim(filename))

call nc_check(nf90_inq_varid(ncFileID, 'physical_tile_fractions', VarID), &
             'nc_write_model_atts', 'inq_varid physical_tile_fractions '//trim(filename))
call nc_check(nf90_put_var(ncFileID, VarID, physical_tile_fractions ), &
             'nc_write_model_atts', 'put_var physical_tile_fractions '//trim(filename))

call nc_check(nf90_inq_varid(ncFileID, 'land_mask', VarID), &
             'nc_write_model_atts', 'inq_varid land_mask '//trim(filename))
call nc_check(nf90_put_var(ncFileID, VarID, land_mask ), &
             'nc_write_model_atts', 'put_var land_mask '//trim(filename))

call nc_check(nf90_inq_varid(ncFileID, 'land_tile_fractions', VarID), &
             'nc_write_model_atts', 'inq_varid land_tile_fractions '//trim(filename))
call nc_check(nf90_put_var(ncFileID, VarID, land_tile_fractions ), &
             'nc_write_model_atts', 'put_var land_tile_fractions '//trim(filename))

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
!>

function nc_write_model_vars( ncFileID, state_vec, copyindex, timeindex ) result (ierr)

integer,  intent(in) :: ncFileID      ! netCDF file identifier
real(r8), intent(in) :: state_vec(:)
integer,  intent(in) :: copyindex
integer,  intent(in) :: timeindex
integer              :: ierr          ! return value of function

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, ncstart, nccount
character(len=NF90_MAX_NAME)          :: varname
integer :: i, ivar, VarID, ncNdims, dimlen
integer :: TimeDimID, CopyDimID

real(r8), allocatable, dimension(:)       :: data_1d_array
real(r8), allocatable, dimension(:,:)     :: data_2d_array
real(r8), allocatable, dimension(:,:,:)   :: data_3d_array

character(len=128) :: filename

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

! We need to process the prognostic variables.

do ivar = 1,nfields  ! Very similar to loop in dart_to_jules_restart

   varname = trim(progvar(ivar)%varname)
   string2 = trim(filename)//' '//trim(varname)

   ! Ensure netCDF variable is conformable with progvar quantity.
   ! The TIME and Copy dimensions are intentionally not queried
   ! by looping over the dimensions stored in the progvar type.

   call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
         'nc_write_model_vars', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
         'nc_write_model_vars', 'inquire '//trim(string2))

   ncstart = 1   ! These are arrays, actually
   nccount = 1
   DimCheck : do i = 1,progvar(ivar)%numdims

      write(string1,'(a,i2,A)') 'inquire dimension ',i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
            'nc_write_model_vars', trim(string1))

      if (progvar(ivar)%dimnames(i) == 'time') cycle DimCheck

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*)trim(string2),' dim/dimlen ',i,dimlen, &
                         ' not ',progvar(ivar)%dimlens(i)
         write(string2,*)' but it should be.'
         call error_handler(E_ERR, 'nc_write_model_vars', trim(string1), &
                         source, revision, revdate, text2=trim(string2))
      endif

      nccount(i) = dimlen

   enddo DimCheck

   where(dimIDs == CopyDimID) ncstart = copyindex
   where(dimIDs == CopyDimID) nccount = 1
   where(dimIDs == TimeDimID) ncstart = timeindex
   where(dimIDs == TimeDimID) nccount = 1

   if (do_output() .and. debug > 9) then
      write(*,*)'nc_write_model_vars '//trim(varname)//' start is ',ncstart(1:ncNdims)
      write(*,*)'nc_write_model_vars '//trim(varname)//' count is ',nccount(1:ncNdims)
   endif

   ! Since 'time' is a singleton dimension, we can use the same logic
   ! as if it the variable had one less dimension.

   if (     (progvar(ivar)%numdims == 1) .or. &
           ((progvar(ivar)%numdims == 2) .and. &
            (progvar(ivar)%dimnames(2) == 'time')) )then

      if ( ncNdims /= 3 ) then
         write(string1,*)trim(varname),' no room for copy,time dimensions.'
         write(string2,*)'netcdf file should have 3 dimensions, has ',ncNdims
         call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                         source, revision, revdate, text2=string2)
      endif

      allocate(data_1d_array( progvar(ivar)%dimlens(1) ) )
      call vector_to_prog_var(state_vec, ivar, data_1d_array)
      call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array, &
          start = ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
                'nc_write_model_vars', 'put_var '//trim(string2))
      deallocate(data_1d_array)

   elseif ( (progvar(ivar)%numdims == 2) .or. &
           ((progvar(ivar)%numdims == 3) .and. &
            (progvar(ivar)%dimnames(3) == 'time')) )then

      if ( ncNdims /= 4 ) then
         write(string1,*)trim(varname),' no room for copy,time dimensions.'
         write(string2,*)'netcdf file should have 4 dimensions, has ',ncNdims
         call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                         source, revision, revdate, text2=string2)
      endif

      allocate(data_2d_array( progvar(ivar)%dimlens(1),  &
                              progvar(ivar)%dimlens(2) ))
      call vector_to_prog_var(state_vec, ivar, data_2d_array)
      call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
          start = ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
                'nc_write_model_vars', 'put_var '//trim(string2))
      deallocate(data_2d_array)

   elseif ( (progvar(ivar)%numdims == 3) .or. &
           ((progvar(ivar)%numdims == 4) .and. &
            (progvar(ivar)%dimnames(4) == 'time')) )then

      if ( ncNdims /= 5 ) then
         write(string1,*)trim(varname),' no room for copy,time dimensions.'
         write(string2,*)'netcdf file should have 5 dimensions, has ',ncNdims
         call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                         source, revision, revdate, text2=string2)
      endif

      allocate(data_3d_array( progvar(ivar)%dimlens(1),  &
                              progvar(ivar)%dimlens(2),  &
                              progvar(ivar)%dimlens(3) ))
      call vector_to_prog_var(state_vec, ivar, data_3d_array)
      call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
          start = ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
                'nc_write_model_vars', 'put_var '//trim(string2))
      deallocate(data_3d_array)

   else

      write(string1,*)'Do not know how to handle JULES variables with > 3 dimensions'
      write(string2,*)trim(progvar(ivar)%varname),'has shape', &
                           progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      call error_handler(E_ERR,'nc_write_model_vars',string1,source,revision,revdate)

   endif

enddo

! Flush the buffer and leave netCDF file open

call nc_check(nf90_sync(ncFileID), 'nc_write_model_vars', 'sync '//trim(filename))

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars


!-----------------------------------------------------------------------
!>
!> Perturbs a single model state for generating initial ensembles.
!> This (required interface) is unsupported in JULES and any attempt
!> to use it will cause DART to terminate. Initial ensemble members
!> are generated externally for JULES applications.
!>
!> A model may choose to provide a NULL INTERFACE by returning
!> .false. for the interf_provided argument. This indicates to
!> the filter that if it needs to generate perturbed states, it
!> may do so by adding a perturbation to each model state
!> variable independently. The interf_provided argument
!> should be returned as .true. if the model wants to do its own
!> perturbing of states.
!> 
!> @param state the base DART state vector to perturb
!> @param pert_state the (new) perturbed DART state vector
!> @param interf_provided logical flag that indicates that this routine
!>               is unique for JULES. TRUE means this routine will
!>               somehow create the perturbed state, FALSE means
!>               the default perturb routine will be used.
!> 

subroutine pert_model_state(state, pert_state, interf_provided)

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

integer :: i
logical, save :: random_seq_init = .false.

if ( .not. module_initialized ) call static_init_model

call error_handler(E_ERR,'pert_model_state', &
                  'JULES cannot be started from a single vector', &
                  source, revision, revdate, &
                  text2='see comments in jules/model_mod.f90::pert_model_state()')

interf_provided = .true.

! Initialize my random number sequence
if(.not. random_seq_init) then
   call init_random_seq(random_seq, my_task_id())
   random_seq_init = .true.
endif

! This does keep the compiler error messages down, but this
! section of code must never be reached.
do i=1,size(state)
   pert_state(i) = random_gaussian(random_seq, state(i), &
                                   model_perturbation_amplitude)
enddo

return
end subroutine pert_model_state


!-----------------------------------------------------------------------
!>
!> Given a DART location (referred to as "base") and a set of candidate
!> locations & kinds (obs, obs_kind), returns the subset close to the
!> "base", their indices, and their distances to the "base" ...
!>
!> @param gc precomputed 'get_close_type' to speed up candidate selection 
!> @param base_obs_loc location of the observation in question
!> @param base_obs_kind DART KIND of observation in question
!> @param locs array of comparison locations
!> @param loc_kind matching array of KINDs for the comparison locations
!> @param num_close the number of locs locations that are within the prespecified distance (information contained in 'gc')
!> @param close_ind the indices of the locs locations that are 'close'
!> @param dist the distances of each of the close locations.
!>

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, &
                         locs, loc_kind, num_close, close_ind, dist)

type(get_close_type), intent(in) :: gc
type(location_type),  intent(in) :: base_obs_loc
integer,              intent(in) :: base_obs_kind
type(location_type),  intent(in) :: locs(:)
integer,              intent(in) :: loc_kind(:)
integer,              intent(out):: num_close
integer,              intent(out):: close_ind(:)
real(r8),  OPTIONAL,  intent(out):: dist(:)

integer :: k

if ( .not. module_initialized ) call static_init_model

! Initialize variables to missing status

num_close = 0
close_ind = -99
if (present(dist)) dist = 1.0e9   ! something big and positive (far away)

! Get all the potentially close obs but no dist (optional argument dist(:)
! is not present) This way, we are decreasing the number of distance
! computations that will follow.  This is a horizontal-distance operation and
! we don't need to have the relevant vertical coordinate information yet
! (for obs).
!> @todo FIXME - confirm that the close_ind() array does not benefit from having 
! all the dry_land locations pruned out.

call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, locs, loc_kind, &
                       num_close, close_ind)

! Loop over potentially close subset of obs priors or state variables
if (present(dist)) then
do k = 1, num_close

   dist(k) = get_dist(base_obs_loc,       locs(close_ind(k)), &
                      base_obs_kind, loc_kind(close_ind(k)))

enddo
   !> @todo FIXME
endif

return
end subroutine get_close_obs


!-----------------------------------------------------------------------
!>
!> If needed by the model interface, this is the current mean
!> for all state vector items across all ensembles. The ensemble mean
!> may or may not be needed by JULES.
!>
!> @param filter_ens_mean the ensemble mean DART state vector
!>

subroutine ens_mean_for_model(filter_ens_mean)

real(r8), intent(in) :: filter_ens_mean(:)

if ( .not. module_initialized ) call static_init_model

! ens_mean = filter_ens_mean

return
end subroutine ens_mean_for_model


!==================================================================
! The remaining PUBLIC interfaces come next
!==================================================================


!-----------------------------------------------------------------------
!>
!> Reads the current time and state variables from a JULES restart
!> file and packs them into a dart state vector. This better happen
!> in the same fashion as the metadata arrays are built.
!>

subroutine jules_to_dart_state_vector(state_vector, restart_time)

real(r8),         intent(out) :: state_vector(:)
type(time_type),  intent(out) :: restart_time

! temp space to hold data while we are reading it
integer  :: i, j, k, ni, nj, nk, ivar, indx
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array
real(r8), allocatable, dimension(:,:,:)     :: data_3d_array
real(r8), allocatable, dimension(:,:,:,:)   :: data_4d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: varname
integer :: io, TimeDimID, VarID, ncNdims, dimlen
integer :: ncid
character(len=256) :: myerrorstring

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

call nc_check(nf90_open(trim(jules_output_filename), NF90_NOWRITE, ncid), &
              'jules_to_dart_state_vector', 'open for time '//jules_output_filename)

restart_time = get_state_time(ncid)

if (do_output()) then
   call print_time(restart_time,'time of model state '//jules_output_filename)
   call print_date(restart_time,'date of model state '//jules_output_filename)
endif

call nc_check(nf90_close(ncid), 'jules_to_dart_state_vector', &
                                'close '//jules_output_filename)

! Start counting and filling the state vector one item at a time,
! repacking the Nd arrays into a single 1d list of numbers.
! The DART state vector may be composed of variables from either the restart 
! or output file. Sometimes it is useful to have variables from the output 
! file to use for forward ops.

do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   myerrorstring = trim(progvar(ivar)%origin)//' '//trim(progvar(ivar)%varname)

   call nc_check(nf90_open(trim(progvar(ivar)%origin), NF90_NOWRITE, ncid), &
              'jules_to_dart_state_vector','open '//trim(myerrorstring))

   ! File is not required to have a time dimension
   io = nf90_inq_dimid(ncid, 'time', TimeDimID)
   if (io /= NF90_NOERR) TimeDimID = MISSING_I

   call nc_check(nf90_inq_varid(ncid, varname, VarID), &
            'jules_to_dart_state_vector', 'inq_varid '//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=ncNdims), &
            'jules_to_dart_state_vector', 'inquire '//trim(myerrorstring))

   ! Check the rank of the variable

   if ( ncNdims /= progvar(ivar)%numdims ) then
      write(string1, *) 'netCDF rank of '//trim(varname)// &
                        ' does not match derived type knowledge'
      write(string2, *) 'netCDF rank is ',ncNdims,' expected ',progvar(ivar)%numdims
      call error_handler(E_ERR,'jules_to_dart_state_vector', string1, &
                        source,revision,revdate,text2=string2)
   endif

   ! Check the shape of the variable

   do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'jules_to_dart_state_vector', string1)

      ! Time dimension will be unity in progvar, but not necessarily
      ! in origin file. We only want a single matching time.
      ! static_init_model() only reserves space for a single time.

      if ( dimIDs(i) == TimeDimID ) dimlen = 1

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(myerrorstring),' dim/dimlen ',i,dimlen, &
                              ' not ',progvar(ivar)%dimlens(i)
         call error_handler(E_ERR, 'jules_to_dart_state_vector', string1, &
                            source, revision, revdate)
      endif

   enddo

   ! Pack the variable into the DART state vector

   indx = progvar(ivar)%index1

   if (ncNdims == 1) then

      ni = progvar(ivar)%dimlens(1)
      allocate(data_1d_array(ni))
      call DART_get_var(ncid, varname, data_1d_array)

      do i = 1, ni
         state_vector(indx) = data_1d_array(i)
         indx = indx + 1
      enddo
      deallocate(data_1d_array)

   elseif (ncNdims == 2) then

      ! restart file variables are at most 2D, and never have time
      !   float canopy(tile, land) ;
      !   float t_soil(soil, land) ;
      !   float cs(scpool, land) ;
      !   float gs(land) ;

      ni = progvar(ivar)%dimlens(1)
      nj = progvar(ivar)%dimlens(2)
      allocate(data_2d_array(ni, nj))
      call DART_get_var(ncid, varname, data_2d_array)

      do j = 1, nj
      do i = 1, ni
         state_vector(indx) = data_2d_array(i, j)
         indx = indx + 1
      enddo
      enddo
      deallocate(data_2d_array)

   elseif (ncNdims == 3) then

      ! restart file variables never have 3 dimensions
      ! output  file variables may have 3 dimensions
      !   float precip(time, y, x) ;

      if     ( (trim(progvar(ivar)%dimnames(1)) == 'x')   .and. &
               (trim(progvar(ivar)%dimnames(2)) == 'y')   .and. &
               (trim(progvar(ivar)%dimnames(3)) == 'time') ) then

         ni = progvar(ivar)%dimlens(1)
         nj = progvar(ivar)%dimlens(2)
       ! nk = progvar(ivar)%dimlens(3) not needed ... time is always a singleton

         allocate(data_3d_array(ni, nj, 1))
         call DART_get_var(ncid, varname, data_3d_array)

         do j = 1, nj
         do i = 1, ni
            state_vector(indx) = data_3d_array(i, j, 1)
            indx = indx + 1
         enddo
         enddo
         deallocate(data_3d_array)
      else

         write(string1, *) '3D variable unexpected shape -- only support x, y, time(=1)'
         write(string2, *) 'variable [',trim(progvar(ivar)%varname),']'
         write(string3, *) 'file [',trim(progvar(ivar)%origin),']'
         call error_handler(E_ERR,'jules_to_dart_state_vector', string1, &
                           source, revision, revdate, text2=string2, text3=string3)

      endif

   elseif (ncNdims == 4) then

      ! output  file variables may have 4 dimensions
      !   float soil_wet(time, soil, y, x)
      !   float    esoil(time, tile, y, x)

      if     ( (trim(progvar(ivar)%dimnames(1)) == 'x')   .and. &
               (trim(progvar(ivar)%dimnames(2)) == 'y')   .and. &
               (trim(progvar(ivar)%dimnames(4)) == 'time') ) then

         ni = progvar(ivar)%dimlens(1)
         nj = progvar(ivar)%dimlens(2)
         nk = progvar(ivar)%dimlens(3)
       ! nl = progvar(ivar)%dimlens(3) not needed ... time is always a singleton

         allocate(data_4d_array(ni, nj, nk, 1))
         call DART_get_var(ncid, varname, data_4d_array)

         do k = 1, nk
         do j = 1, nj
         do i = 1, ni
            state_vector(indx) = data_4d_array(i, j, k, 1)
            indx = indx + 1
         enddo
         enddo
         enddo
         deallocate(data_4d_array)
      else

         write(string1, *) '4D variable unexpected shape -- only support x, y, xxx, time(=1)'
         write(string2, *) 'variable [',trim(progvar(ivar)%varname),']'
         write(string3, *) 'file [',trim(progvar(ivar)%origin),']'
         call error_handler(E_ERR,'jules_to_dart_state_vector', string1, &
                           source, revision, revdate, text2=string2, text3=string3)

      endif

   else

      write(string1, *) 'no support for data array of dimension ', ncNdims
      write(string2, *) 'variable [',trim(progvar(ivar)%varname),']'
      write(string3, *) 'file [',trim(progvar(ivar)%origin),']'
      call error_handler(E_ERR,'jules_to_dart_state_vector', string1, &
                        source, revision, revdate, text2=string2, text3=string3)
   endif

   indx = indx - 1
   if ( indx /= progvar(ivar)%indexN ) then
      write(string1, *)'Variable '//trim(varname)//' filled wrong.'
      write(string2, *)'Should have ended at ',progvar(ivar)%indexN, &
                       ' actually ended at ',indx
      call error_handler(E_ERR,'jules_to_dart_state_vector', string1, &
                        source,revision,revdate,text2=string2)
   endif

   call nc_check(nf90_close(ncid),'jules_to_dart_state_vector', &
                                  'close '//progvar(ivar)%origin)
   ncid = 0

enddo

return
end subroutine jules_to_dart_state_vector


!-----------------------------------------------------------------------
!>
!> Writes the current time and state variables from a dart state
!> vector (1d array) into a JULES netcdf restart file.
!>

subroutine dart_to_jules_restart( state_vector )

real(r8),         intent(in) :: state_vector(:)

! temp space to hold data while we are writing it
integer :: i, ni, nj, ivar
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname
integer         :: VarID, ncNdims, dimlen
integer         :: ncFileID

call error_handler(E_MSG, 'dart_to_jules_restart:', 'FIXME RAFAEL routine not tested', source, revision, revdate)

if ( .not. module_initialized ) call static_init_model

! Check that the output file exists ...

if ( .not. file_exist(jules_restart_filename) ) then
   write(string1,*) 'cannot open file ', trim(jules_restart_filename),' for writing.'
   call error_handler(E_ERR,'dart_to_jules_restart',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(jules_restart_filename), NF90_WRITE, ncFileID), &
             'dart_to_jules_restart','open '//trim(jules_restart_filename))

!> @todo FIXME ... if/when the restart file gets metadata indicating the time,
! it would be a good idea to check the model_time against the time in
! the restart file.

! The DART prognostic variables are only defined for a single time.
! We already checked the assumption that variables are xy2d or xyz3d ...
! IF the netCDF variable has a TIME dimension, it must be the last dimension,
! and we need to read the LAST timestep and effectively squeeze out the
! singleton dimension when we stuff it into the DART state vector.

UPDATE : do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   string2 = trim(progvar(ivar)%origin)//' '//trim(varname)

   if ( .not. progvar(ivar)%update ) then
      write(string1,*)'intentionally not updating '//trim(string2)
      write(string3,*)'as per namelist control in model_nml:variables'
      call error_handler(E_MSG, 'dart_to_jules_restart:', string1, text2=string3)
      cycle UPDATE
   endif

   ! Ensure netCDF variable is conformable with progvar quantity.
   ! The TIME and Copy dimensions are intentionally not queried
   ! by looping over the dimensions stored in the progvar type.

   call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'dart_to_jules_restart', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
            'dart_to_jules_restart', 'inquire '//trim(string2))

   DimCheck : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
            'dart_to_jules_restart', string1)

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*)trim(string2),' dim/dimlen ',i,dimlen, &
                         ' not ',progvar(ivar)%dimlens(i)
         write(string2,*)' but it should be.'
         call error_handler(E_ERR, 'dart_to_jules_restart', string1, &
                         source, revision, revdate, text2=string2)
      endif

   enddo DimCheck

   ! When called with a 4th argument, vector_to_prog_var()
   ! clamps to physically meaningful values as specified in the model_mod namelist.

   if (progvar(ivar)%numdims == 1) then

      ni = progvar(ivar)%dimlens(1)
      allocate(data_1d_array(ni))
      call vector_to_prog_var(state_vector, ivar, data_1d_array, ncFileID)

      call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array), &
            'dart_to_jules_restart', 'put_var '//trim(varname))
      deallocate(data_1d_array)

   elseif (progvar(ivar)%numdims == 2) then

      ni = progvar(ivar)%dimlens(1)
      nj = progvar(ivar)%dimlens(2)
      allocate(data_2d_array(ni, nj))
      call vector_to_prog_var(state_vector, ivar, data_2d_array, ncFileID)

      call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array), &
            'dart_to_jules_restart', 'put_var '//trim(varname))
      deallocate(data_2d_array)

   else
      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'dart_to_jules_restart', string1, &
                        source,revision,revdate)
   endif

   !> @todo FIXME ... this works perfectly if it were not for a bug in netCDF.
   ! When they fix the bug, this will be a useful thing to restore.
   ! Make note that the variable has been updated by DART
!  call nc_check(nf90_Redef(ncFileID),'dart_to_jules_restart', &
!                                     'redef '//trim(jules_restart_filename))
!  call nc_check(nf90_put_att(ncFileID, VarID,'DART','variable modified by DART'),&
!                'dart_to_jules_restart', 'modified '//trim(varname))
!  call nc_check(nf90_enddef(ncfileID),'dart_to_jules_restart',&
!                                     'state enddef '//trim(jules_restart_filename))

enddo UPDATE

call nc_check(nf90_close(ncFileID),'dart_to_jules_restart', &
                                   'close '//trim(jules_restart_filename))
ncFileID = 0

return
end subroutine dart_to_jules_restart


!-----------------------------------------------------------------------
!>
!> provides access to the filename normally in module storage
!> the filename originally comes from the dart namelist.
!>

subroutine get_jules_restart_filename( filename )

character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(jules_restart_filename)

return
end subroutine get_jules_restart_filename


!-----------------------------------------------------------------------
!>
!> Since the JULES restart files do not have the time as part of the
!> metadata, the times must be read from the companion output file.
!> The last time in the output file is the time used as the valid time
!> of the model state.
!>
!>   float time(time) ;
!>           time:standard_name = "time" ;
!>           time:long_name = "Time of data" ;
!>           time:units = "seconds since 2014-01-01 03:00:00" ;
!>           time:bounds = "time_bounds" ;
!>           time:calendar = "standard" ;
!>

function get_state_time_ncid( ncid )

type(time_type) :: get_state_time_ncid
integer, intent(in) :: ncid

integer :: VarID, numdims, xtype, ios, dimlen
integer :: year, month, day, hour, minute, second

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: attvalue

real(r8), allocatable :: time_array(:)

type(time_type) :: base_time
type(time_type) :: forecast_length

if ( .not. module_initialized ) call static_init_model

call nc_check(nf90_Inquire(ncid, nDimensions, nVariables, nAttributes, unlimitedDimID),&
              'get_state_time_ncid:', 'inquire '//trim(jules_output_filename))

call nc_check(nf90_inq_varid(ncid, 'time', VarID), 'get_state_time_ncid:', &
              'inq_varid time '//trim(jules_output_filename))
call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_state_time_ncid:', &
              'inquire_variable time '//trim(jules_output_filename))

if (numdims /= 1) then
   write(string1,*)'"time" is supposed to be 1D.'
   write(string2,*)'"time" has ',numdims,' dimensions.'
   call error_handler(E_ERR, 'get_state_time_ncid:', string1, &
          source, revision, revdate, text2=string2)
endif

if (dimIDs(1) /= unlimitedDimID) then
   write(string1,*)'"time" is supposed to be the unlimited dimension.'
   call error_handler(E_ERR, 'get_state_time_ncid:', string1, &
          source, revision, revdate)
endif

call nc_check(nf90_inquire_dimension(ncid, unlimitedDimID, len=dimlen), &
        'get_state_time_ncid', 'inquire_dimension time '//trim(jules_output_filename))

! Make sure 'time' is the unlimited dimension and just grab the whole thing,
! use the LAST one ...
! Get the time units attribute so we can add the offset to the base, etc.

call nc_check(nf90_get_att(ncid, VarID, 'units', attvalue), &
        'get_state_time_ncid:', 'time get_att units '//varname)

! Make sure the calendar is 'standard', which implies a gregorian calendar

! time:units = "seconds since 2014-01-01 03:00:00" ;
!               1234567890123

if (attvalue(1:13) /= 'seconds since') then
   write(string1,*)'expecting time units of [seconds since ... ]'
   write(string2,*)'read time units of ['//trim(attvalue)//']'
   call error_handler(E_ERR, 'get_state_time_ncid:', string1, &
          source, revision, revdate, text2=string2)
endif

read(attvalue,'(14x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
if (ios /= 0) then
   write(string1,*)'Unable to read time units. Error status was ',ios
   write(string2,*)'expected "seconds since YYYY-MM-DD HH:MM:SS"'
   write(string3,*)'was      "'//trim(attvalue)//'"'
   call error_handler(E_ERR, 'get_state_time_ncid:', string1, &
          source, revision, revdate, text2=string2, text3=string3)
endif

allocate(time_array(dimlen))

call nc_check(nf90_get_var(ncid, VarID, time_array), 'get_state_time_ncid', &
                      &            'get_var time '//trim(jules_output_filename))

forecast_length = set_time(int(time_array(dimlen)),0)
base_time = set_date(year, month, day, hour, minute, second)

deallocate(time_array)

get_state_time_ncid = base_time + forecast_length

end function get_state_time_ncid


!-----------------------------------------------------------------------
!>
!> sometimes it is useful to use the netCDF file name
!>

function get_state_time_fname(filename)

type(time_type) :: get_state_time_fname
character(len=*), intent(in) :: filename

integer :: ncid

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'get_state_time_fname',string1,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'get_state_time_fname', 'open '//trim(filename))

get_state_time_fname = get_state_time_ncid(ncid)

call nc_check(nf90_close(ncid),'get_state_time_fname', 'close '//trim(filename))

end function get_state_time_fname


!-----------------------------------------------------------------------
!>
!> This function will return a R8 array with the netCDF attributes applied.
!> scale_factor, offset will be applied,
!> missing_value, _FillValue will be replaced by the DART missing value ...
!>
!> If _FillValue is defined then it should be scalar and of the same type as the variable.
!> If the variable is packed using scale_factor and add_offset attributes (see below),
!> the _FillValue attribute should have the data type of the packed data.
!>
!> missing_value
!> When scale_factor and add_offset are used for packing, the value(s) of the missing_value
!> attribute should be specified in the domain of the data in the file (the packed data),
!> so that missing values can be detected before the scale_factor and add_offset are applied.
!>
!> scale_factor
!> If present for a variable, the data are to be multiplied by this factor after the data
!> are read by the application that accesses the data.  If valid values are specified using
!> the valid_min, valid_max, valid_range, or _FillValue attributes, those values should be
!> specified in the domain of the data in the file (the packed data), so that they can be
!> interpreted before the scale_factor and add_offset are applied.
!>
!> add_offset
!> If present for a variable, this number is to be added to the data after it is read by
!> the application that accesses the data. If both scale_factor and add_offset attributes
!> are present, the data are first scaled before the offset is added.
!>

subroutine get_var_1d(ncid, varname, var1d, context)

integer,                    intent(in)  :: ncid
character(len=*),           intent(in)  :: varname
real(r8), dimension(:),     intent(out) :: var1d
character(len=*), OPTIONAL, intent(in)  :: context

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens, ncstart, nccount
integer  :: VarID, numdims, xtype, io1, io2
integer  :: TimeDimID, time_dimlen, timeindex
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8
character(len=512) :: msgstring

integer,  allocatable, dimension(:) :: intarray
real(r4), allocatable, dimension(:) :: r4array

if ( .not. module_initialized ) call static_init_model

! a little whitespace makes this a lot more readable
if (do_output() .and. debug > 8) then
   write(*,*)
   write(logfileunit,*)
endif

if (present(context)) then
   msgstring = context
else
   msgstring = varname
endif

io1 = nf90_inq_dimid(ncid, 'time', TimeDimID)
if (io1 /= NF90_NOERR) TimeDimID = MISSING_I

call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), &
              'get_var_1d', 'inq_varid '//trim(msgstring))
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_1d', 'inquire_variable '//trim(msgstring))
call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlens(1)), &
              'get_var_1d', 'inquire_dimension '//trim(msgstring))

if ((numdims /= 1) .or. (size(var1d) /= dimlens(1)) ) then
   write(string1,*) trim(varname)//' is not the expected shape/length of ', size(var1d)
   call error_handler(E_ERR,'get_var_1d',string1,source,revision,revdate)
endif

ncstart = 1
nccount = dimlens(1)

if (dimIDs(1) == TimeDimID) then
   call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), &
         'get_var_1d', 'inquire_dimension time '//trim(msgstring))
   timeindex  = FindDesiredTimeIndx(ncid, time_dimlen, varname)
   ncstart(1) = timeindex
   nccount(1) = 1
endif

if (do_output() .and. debug > 8) then
   write(*,*)'get_var_1d: variable ['//trim(varname)//']'
   write(*,*)'get_var_1d: start ',ncstart(1:numdims)
   write(*,*)'get_var_1d: count ',nccount(1:numdims)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1)))
   call nc_check(nf90_get_var(ncid, VarID, values=intarray, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_1d', 'get_var '//trim(msgstring))
   var1d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d:',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d:',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d + add_offset
   endif

   deallocate(intarray)

elseif (xtype == NF90_FLOAT) then

   allocate(r4array(dimlens(1)))
   call nc_check(nf90_get_var(ncid, VarID, values=r4array, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_1d', 'get_var '//trim(msgstring))
   var1d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d:',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d:',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d + add_offset
   endif

   deallocate(r4array)

elseif (xtype == NF90_DOUBLE) then

   call nc_check(nf90_get_var(ncid, VarID, values=var1d, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_1d', 'get_var '//trim(msgstring))

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var1d == spvalR8) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d:',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var1d == spvalR8) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d:',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d + add_offset
   endif

else
   write(string1,*) trim(msgstring)//' has unsupported (by DART) xtype of', xtype
   call error_handler(E_ERR,'get_var_1d',string1,source,revision,revdate)
endif

return
end subroutine get_var_1d


!-----------------------------------------------------------------------
!>
!> This function will return a R8 array with the netCDF attributes applied.
!> scale_factor, offset will be applied,
!> missing_value, _FillValue will be replaced by the DART missing value ...
!>
!> If _FillValue is defined then it should be scalar and of the same type as the variable.
!> If the variable is packed using scale_factor and add_offset attributes (see below),
!> the _FillValue attribute should have the data type of the packed data.
!>
!> missing_value
!> When scale_factor and add_offset are used for packing, the value(s) of the missing_value
!> attribute should be specified in the domain of the data in the file (the packed data),
!> so that missing values can be detected before the scale_factor and add_offset are applied.
!>
!> scale_factor
!> If present for a variable, the data are to be multiplied by this factor after the data
!> are read by the application that accesses the data.  If valid values are specified using
!> the valid_min, valid_max, valid_range, or _FillValue attributes, those values should be
!> specified in the domain of the data in the file (the packed data), so that they can be
!> interpreted before the scale_factor and add_offset are applied.
!>
!> add_offset
!> If present for a variable, this number is to be added to the data after it is read by
!> the application that accesses the data. If both scale_factor and add_offset attributes
!> are present, the data are first scaled before the offset is added.
!>

subroutine get_var_2d(ncid, varname, var2d, context)

integer,                    intent(in)  :: ncid
character(len=*),           intent(in)  :: varname
real(r8), dimension(:,:),   intent(out) :: var2d
character(len=*), OPTIONAL, intent(in)  :: context

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens, ncstart, nccount
integer  :: VarID, numdims, xtype, io1, io2, i
integer  :: TimeDimID, time_dimlen, timeindex
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8
character(len=512) :: msgstring

integer,  allocatable, dimension(:,:) :: intarray
real(r4), allocatable, dimension(:,:) :: r4array

if ( .not. module_initialized ) call static_init_model

! a little whitespace makes this a lot more readable
if (do_output() .and. debug > 8) then
   write(*,*)
   write(logfileunit,*)
endif

if (present(context)) then
   msgstring = context
else
   msgstring = varname
endif

io1 = nf90_inq_dimid(ncid, 'time', TimeDimID)
if (io1 /= NF90_NOERR) TimeDimID = MISSING_I

call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), 'get_var_2d', 'inq_varid')
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_2d', 'inquire_variable')

if ( (numdims /= 2)  ) then
   write(string1,*) trim(varname)//' is not a 2D variable as expected.'
   call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate)
endif

ncstart(:) = 1
nccount(:) = 1

DimCheck : do i = 1,numdims

   write(string1,'(a,i2,1x,A)') 'inquire dimension ',i,trim(msgstring)

   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
              'get_var_2d', string1)

   if ( dimIDs(i) == TimeDimID ) then
      call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen), &
            'get_var_2d', 'inquire_dimension time '//trim(msgstring))
      timeindex  = FindDesiredTimeIndx(ncid, time_dimlen, varname)
      ncstart(i) = timeindex
      dimlens(i) = 1

   elseif ( size(var2d,i) /= dimlens(i) ) then
      write(string1,*) trim(varname)//' has shape ', dimlens(1:numdims)
      write(string2,*) 'which is not the expected shape of ', size(var2d,1),size(var2d,2)
      call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate,text2=string2)
   endif

   nccount(i) = dimlens(i)

enddo DimCheck

if (do_output() .and. debug > 8) then
   write(*,*)'get_var_2d: variable ['//trim(varname)//']'
   write(*,*)'get_var_2d: start ',ncstart(1:numdims)
   write(*,*)'get_var_2d: count ',nccount(1:numdims)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1),dimlens(2)))
   call nc_check(nf90_get_var(ncid, VarID, values=intarray, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_2d', 'get_var '//trim(msgstring))
   var2d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d:',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d:',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d + add_offset
   endif

   deallocate(intarray)

elseif (xtype == NF90_FLOAT) then

   allocate(r4array(dimlens(1),dimlens(2)))
   call nc_check(nf90_get_var(ncid, VarID, values=r4array, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_2d', 'get_var '//trim(msgstring))
   var2d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d:',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d:',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d + add_offset
   endif

   deallocate(r4array)

elseif (xtype == NF90_DOUBLE) then

   call nc_check(nf90_get_var(ncid, VarID, values=var2d, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_2d', 'get_var '//trim(msgstring))

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var2d == spvalR8) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d:',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var2d == spvalR8) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d:',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d + add_offset
   endif

else
   write(string1,*) trim(varname)//' has unsupported (by DART) xtype of', xtype
   call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate)
endif

return
end subroutine get_var_2d


!-----------------------------------------------------------------------
!>
!> This function will return a R8 array with the netCDF attributes applied.
!> scale_factor, offset will be applied,
!> missing_value, _FillValue will be replaced by the DART missing value ...
!>
!> If _FillValue is defined then it should be scalar and of the same type as the variable.
!> If the variable is packed using scale_factor and add_offset attributes (see below),
!> the _FillValue attribute should have the data type of the packed data.
!>
!> missing_value
!> When scale_factor and add_offset are used for packing, the value(s) of the missing_value
!> attribute should be specified in the domain of the data in the file (the packed data),
!> so that missing values can be detected before the scale_factor and add_offset are applied.
!>
!> scale_factor
!> If present for a variable, the data are to be multiplied by this factor after the data
!> are read by the application that accesses the data.  If valid values are specified using
!> the valid_min, valid_max, valid_range, or _FillValue attributes, those values should be
!> specified in the domain of the data in the file (the packed data), so that they can be
!> interpreted before the scale_factor and add_offset are applied.
!>
!> add_offset
!> If present for a variable, this number is to be added to the data after it is read by
!> the application that accesses the data. If both scale_factor and add_offset attributes
!> are present, the data are first scaled before the offset is added.
!>

subroutine get_var_3d(ncid, varname, var3d, context)

integer,                    intent(in)  :: ncid
character(len=*),           intent(in)  :: varname
real(r8), dimension(:,:,:), intent(out) :: var3d
character(len=*), OPTIONAL, intent(in)  :: context

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens, ncstart, nccount
integer  :: i, TimeDimID, time_dimlen, timeindex
integer  :: VarID, numdims, xtype, io1, io2
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8
character(len=512) :: msgstring

integer,  allocatable, dimension(:,:,:) :: intarray
real(r4), allocatable, dimension(:,:,:) :: r4array

if ( .not. module_initialized ) call static_init_model

! a little whitespace makes this a lot more readable
if (do_output() .and. debug > 8) then
   write(*,*)
   write(logfileunit,*)
endif

if (present(context)) then
   msgstring = context
else
   msgstring = varname
endif

! 3D fields _may_ have a time dimension.
! If so, we need to know the Time Dimension ID and length

if ( nf90_inq_dimid(ncid, 'time', TimeDimID) == NF90_NOERR ) then
   call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), &
         'get_var_3d', 'inquire_dimension time '//trim(msgstring))
else
   TimeDimID = -999    ! an impossible value 
   time_dimlen = 0
endif

call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), 'get_var_3d', 'inq_varid')
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_3d', 'inquire_variable '//trim(msgstring))

if ( (numdims /= 3)  ) then
   write(string1,*) trim(msgstring)//' is not a 3D variable as expected.'
   call error_handler(E_ERR,'get_var_3d',string1,source,revision,revdate)
endif

ncstart(:) = 1
nccount(:) = 1
DimCheck : do i = 1,numdims

   write(string1,'(a,i2,1x,A)') 'inquire dimension ',i,trim(msgstring)

   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
                 'get_var_3d', trim(string1))

   if ( dimIDs(i) == TimeDimID ) then
       timeindex  = FindDesiredTimeIndx(ncid, time_dimlen, varname)
       ncstart(i) = timeindex
       dimlens(i) = 1

   elseif ( size(var3d,i) /= dimlens(i) ) then
      write(string1,*) trim(msgstring)//' dimension ',i,' is ', dimlens(i)
      write(string2,*) 'which is not the expected size of ', size(var3d,i)
      call error_handler(E_ERR, 'get_var_3d', string1, &
                source, revision, revdate, text2=string2)
   endif

   nccount(i) = dimlens(i)

enddo DimCheck

if (do_output() .and. debug > 8) then
   write(*,*)'get_var_3d: context ['//trim(msgstring)//']'
   write(*,*)'get_var_3d: variable ['//trim(varname)//']'
   write(*,*)'get_var_3d: start ',ncstart(1:numdims)
   write(*,*)'get_var_3d: count ',nccount(1:numdims)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1),dimlens(2),dimlens(3)))
   call nc_check(nf90_get_var(ncid, VarID, values=intarray, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_3d', 'get_var '//trim(msgstring))
   var3d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var3d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(msgstring)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d:',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var3d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(msgstring)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d:',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d + add_offset
   endif

   deallocate(intarray)

elseif (xtype == NF90_FLOAT) then

   allocate(r4array(dimlens(1),dimlens(2),dimlens(3)))
   call nc_check(nf90_get_var(ncid, VarID, values=r4array, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_3d', 'get_var '//trim(msgstring))
   var3d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var3d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(msgstring)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d:',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var3d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(msgstring)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d:',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d + add_offset
   endif

   deallocate(r4array)

elseif (xtype == NF90_DOUBLE) then

   call nc_check(nf90_get_var(ncid, VarID, values=var3d, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_3d', 'get_var '//trim(msgstring))

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var3d == spvalR8) var3d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(msgstring)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d:',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var3d == spvalR8) var3d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(msgstring)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d:',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d + add_offset
   endif

else
   write(string1,*) trim(varname)//' has unsupported (by DART) xtype of', xtype
   call error_handler(E_ERR,'get_var_3d',string1,source,revision,revdate)
endif

return
end subroutine get_var_3d


!-----------------------------------------------------------------------
!>
!> This function will return a R8 array with the netCDF attributes applied.
!> scale_factor, offset will be applied,
!> missing_value, _FillValue will be replaced by the DART missing value ...
!>
!> If _FillValue is defined then it should be scalar and of the same type as the variable.
!> If the variable is packed using scale_factor and add_offset attributes (see below),
!> the _FillValue attribute should have the data type of the packed data.
!>
!> missing_value
!> When scale_factor and add_offset are used for packing, the value(s) of the missing_value
!> attribute should be specified in the domain of the data in the file (the packed data),
!> so that missing values can be detected before the scale_factor and add_offset are applied.
!>
!> scale_factor
!> If present for a variable, the data are to be multiplied by this factor after the data
!> are read by the application that accesses the data.  If valid values are specified using
!> the valid_min, valid_max, valid_range, or _FillValue attributes, those values should be
!> specified in the domain of the data in the file (the packed data), so that they can be
!> interpreted before the scale_factor and add_offset are applied.
!>
!> add_offset
!> If present for a variable, this number is to be added to the data after it is read by
!> the application that accesses the data. If both scale_factor and add_offset attributes
!> are present, the data are first scaled before the offset is added.
!>

subroutine get_var_4d(ncid, varname, var4d, context)

integer,                      intent(in)  :: ncid
character(len=*),             intent(in)  :: varname
real(r8), dimension(:,:,:,:), intent(out) :: var4d
character(len=*), OPTIONAL,   intent(in)  :: context

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens, ncstart, nccount
integer  :: i, TimeDimID, time_dimlen, timeindex
integer  :: VarID, numdims, xtype, io1, io2
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8
character(len=512) :: msgstring

integer,  allocatable, dimension(:,:,:,:) :: intarray
real(r4), allocatable, dimension(:,:,:,:) :: r4array

if ( .not. module_initialized ) call static_init_model

! a little whitespace makes this a lot more readable
if (do_output() .and. debug > 8) then
   write(*,*)
   write(logfileunit,*)
endif

if (present(context)) then
   msgstring = context
else
   msgstring = varname
endif

! 4D fields must have a time dimension.
! Need to know the Time Dimension ID and length

call nc_check(nf90_inq_dimid(ncid, 'time', TimeDimID), &
         'get_var_4d', 'inq_dimid time '//trim(msgstring))
call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), &
         'get_var_4d', 'inquire_dimension time '//trim(msgstring))

call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), 'get_var_4d', 'inq_varid')
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_4d', 'inquire_variable '//trim(msgstring))

if ( (numdims /= 4)  ) then
   write(string1,*) trim(varname)//' is not a 4D variable as expected.'
   call error_handler(E_ERR,'get_var_4d',string1,source,revision,revdate)
endif

! only expecting [x,y,[soil,tile],time]

ncstart(:) = 1
nccount(:) = 1
DimCheck : do i = 1,numdims

   write(string1,'(a,i2,1x,A)') 'inquire dimension ',i,trim(varname)

   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
                 'get_var_4d', trim(string1))

   if ( dimIDs(i) == TimeDimID ) then
       timeindex  = FindDesiredTimeIndx(ncid, time_dimlen, varname)
       ncstart(i) = timeindex
       dimlens(i) = 1

   elseif ( size(var4d,i) /= dimlens(i) ) then
      write(string1,*) trim(varname)//' has shape ', dimlens(1:numdims)
      write(string2,*) 'which is not the expected shape of ', size(var4d,i)
      call error_handler(E_ERR,'get_var_4d',string1,source,revision,revdate,text2=string2)
   endif

   nccount(i) = dimlens(i)

enddo DimCheck

if (do_output() .and. debug > 8) then
   write(*,*)'get_var_4d: variable ['//trim(varname)//']'
   write(*,*)'get_var_4d: start ',ncstart(1:numdims)
   write(*,*)'get_var_4d: count ',nccount(1:numdims)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1),dimlens(2),dimlens(3),1))
   call nc_check(nf90_get_var(ncid, VarID, values=intarray, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_4d', 'get_var '//trim(msgstring))
   var4d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var4d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_4d:',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var4d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_4d:',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var4d /= MISSING_R8) var4d = var4d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var4d /= MISSING_R8) var4d = var4d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var4d /= MISSING_R8) var4d = var4d + add_offset
   endif

   deallocate(intarray)

elseif (xtype == NF90_FLOAT) then

   allocate(r4array(dimlens(1),dimlens(2),dimlens(3),1))
   call nc_check(nf90_get_var(ncid, VarID, values=r4array, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_4d', 'get_var '//trim(msgstring))
   var4d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var4d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_4d:',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var4d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_4d:',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var4d /= MISSING_R8) var4d = var4d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var4d /= MISSING_R8) var4d = var4d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var4d /= MISSING_R8) var4d = var4d + add_offset
   endif

   deallocate(r4array)

elseif (xtype == NF90_DOUBLE) then

   call nc_check(nf90_get_var(ncid, VarID, values=var4d, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_4d', 'get_var '//trim(msgstring))

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var4d == spvalR8) var4d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_4d:',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var4d == spvalR8) var4d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() .and. debug > 8 ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_4d:',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var4d /= MISSING_R8) var4d = var4d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var4d /= MISSING_R8) var4d = var4d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var4d /= MISSING_R8) var4d = var4d + add_offset
   endif

else
   write(string1,*) trim(varname)//' has unsupported (by DART) xtype of', xtype
   call error_handler(E_ERR,'get_var_4d',string1,source,revision,revdate)
endif

return
end subroutine get_var_4d


!-----------------------------------------------------------------------
!>  This routine
!>

function get_model_time()
type(time_type) :: get_model_time

if ( .not. module_initialized ) call static_init_model

get_model_time = model_time

end function get_model_time


!==================================================================
! The remaining (private) interfaces come last
!==================================================================


!-----------------------------------------------------------------------
!>
!> convert the values from a 1d array, starting at an offset, into a 1d array.
!>
!> If the optional argument (ncid) is specified, some additional
!> processing takes place.
!>

subroutine vector_to_1d_prog_var(x, ivar, data_1d_array, ncid)

real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:),   intent(out) :: data_1d_array
integer, OPTIONAL,        intent(in)  :: ncid

integer :: i,ii

! unpack the right part of the DART state vector into a 1D array.

ii = progvar(ivar)%index1

do i = 1, progvar(ivar)%dimlens(1)
   data_1d_array(i) = x(ii)
   ii = ii + 1
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_1d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

! Apply the min/max values, if applicable
! This should only be true when converting to a variable that will
! be reinserted into the JULES restart file. This is indicated
! by the presence of the ncid variable.

if (present(ncid)) then

   if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_1d_array /= MISSING_R8) .and. &
             (data_1d_array > progvar(ivar)%maxvalue)) &
              data_1d_array = progvar(ivar)%maxvalue
   endif

   if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_1d_array /= MISSING_R8) .and. &
             (data_1d_array < progvar(ivar)%minvalue)) &
              data_1d_array = progvar(ivar)%minvalue
   endif

   ! Replace the DART missing value flag with the one JULES uses.
   !> @todo FIXME ... I am not sure if there would ever be any missing values
   ! in a JULES restart file.

   if     (progvar(ivar)%xtype == NF90_INT) then
      where(data_1d_array == MISSING_I) data_1d_array = progvar(ivar)%spvalINT
   elseif (progvar(ivar)%xtype == NF90_FLOAT) then
      where(data_1d_array == MISSING_R4) data_1d_array = progvar(ivar)%spvalR4
   elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
      where(data_1d_array == MISSING_R8) data_1d_array = progvar(ivar)%spvalR8
   endif

endif

return
end subroutine vector_to_1d_prog_var


!-----------------------------------------------------------------------
!>
!> convert the values from a 1d array, starting at an offset,
!> into a 2d array.
!>

subroutine vector_to_2d_prog_var(x, ivar, data_2d_array, ncid)

real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:,:), intent(out) :: data_2d_array
integer, OPTIONAL,        intent(in)  :: ncid

integer :: i,j,ii

! unpack the right part of the DART state vector into a 1D array.

ii = progvar(ivar)%index1

do j = 1,progvar(ivar)%dimlens(2)
do i = 1,progvar(ivar)%dimlens(1)
   data_2d_array(i,j) = x(ii)
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

! Apply the min/max values, if applicable
! This should only be true when converting to a variable that will
! be reinserted into the JULES restart file. This is indicated
! by the presence of the ncid variable.

if (present(ncid)) then

   if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_2d_array /= MISSING_R8) .and. &
             (data_2d_array > progvar(ivar)%maxvalue)) &
              data_2d_array = progvar(ivar)%maxvalue
   endif

   if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_2d_array /= MISSING_R8) .and. &
             (data_2d_array < progvar(ivar)%minvalue)) &
              data_2d_array = progvar(ivar)%minvalue
   endif

   ! replace the missing values with the original missing values.
   !> @todo FIXME ... I am not sure if there would ever be any missing values
   ! in a JULES restart file.

   if     (progvar(ivar)%xtype == NF90_INT) then
      where(data_2d_array == MISSING_I) data_2d_array = progvar(ivar)%spvalINT
   elseif (progvar(ivar)%xtype == NF90_FLOAT) then
      where(data_2d_array == MISSING_R4) data_2d_array = progvar(ivar)%spvalR4
   elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
      where(data_2d_array == MISSING_R8) data_2d_array = progvar(ivar)%spvalR8
   endif

endif

return
end subroutine vector_to_2d_prog_var


!-----------------------------------------------------------------------
!>
!> convert the values from a 1d array, starting at an offset,
!> into a 3d array.
!>

subroutine vector_to_3d_prog_var(x, ivar, data_3d_array, ncid)

real(r8), dimension(:),     intent(in)  :: x
integer,                    intent(in)  :: ivar
real(r8), dimension(:,:,:), intent(out) :: data_3d_array
integer, OPTIONAL,          intent(in)  :: ncid

integer :: i,j,k,ii

! unpack the right part of the DART state vector into a 1D array.

ii = progvar(ivar)%index1

do k = 1,progvar(ivar)%dimlens(3)
do j = 1,progvar(ivar)%dimlens(2)
do i = 1,progvar(ivar)%dimlens(1)
   data_3d_array(i,j,k) = x(ii)
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

! Apply the min/max values, if applicable
! This should only be true when converting to a variable that will
! be reinserted into the JULES restart file. This is indicated
! by the presence of the ncid variable.

if (present(ncid)) then

   if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_3d_array /= MISSING_R8) .and. &
             (data_3d_array > progvar(ivar)%maxvalue)) &
              data_3d_array = progvar(ivar)%maxvalue
   endif

   if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_3d_array /= MISSING_R8) .and. &
             (data_3d_array < progvar(ivar)%minvalue)) &
              data_3d_array = progvar(ivar)%minvalue
   endif

   ! replace the missing values with the original missing values.
   !> @todo FIXME ... I am not sure if there would ever be any missing values
   ! in a JULES restart file.

   if     (progvar(ivar)%xtype == NF90_INT) then
      where(data_3d_array == MISSING_I) data_3d_array = progvar(ivar)%spvalINT
   elseif (progvar(ivar)%xtype == NF90_FLOAT) then
      where(data_3d_array == MISSING_R4) data_3d_array = progvar(ivar)%spvalR4
   elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
      where(data_3d_array == MISSING_R8) data_3d_array = progvar(ivar)%spvalR8
   endif

endif

return
end subroutine vector_to_3d_prog_var


!-----------------------------------------------------------------------
!>
!> Read the dimensions from the history netcdf file.
!> The file name comes from module storage ... namelist.
!>

subroutine get_jules_output_dimensions()

integer :: dimid
integer :: ncid

! sets module variables Nx_model, Ny_model, Nsoil, Ntile

call nc_check(nf90_open(trim(jules_output_filename), nf90_nowrite, ncid), &
            'get_jules_output_dimensions','open '//trim(jules_output_filename))

call nc_check(nf90_inq_dimid(ncid, 'x', dimid), &
            'get_jules_output_dimensions','inq_dimid x '//trim(jules_output_filename))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nx_model), &
            'get_jules_output_dimensions', &
            'inquire_dimension x '//trim(jules_output_filename))

call nc_check(nf90_inq_dimid(ncid, 'y', dimid), &
            'get_jules_output_dimensions', &
            'inq_dimid y '//trim(jules_output_filename))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Ny_model), &
            'get_jules_output_dimensions', &
            'inquire_dimension y '//trim(jules_output_filename))

call nc_check(nf90_inq_dimid(ncid, 'soil', dimid), &
            'get_jules_output_dimensions', &
            'inq_dimid soil '//trim(jules_output_filename))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nsoil), &
            'get_jules_output_dimensions', &
            'inquire_dimension soil '//trim(jules_output_filename))

call nc_check(nf90_inq_dimid(ncid, 'tile', dimid), &
            'get_jules_output_dimensions', &
            'inq_dimid tile '//trim(jules_output_filename))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Ntile), &
            'get_jules_output_dimensions', &
            'inquire_dimension tile '//trim(jules_output_filename))

! At this point, we are forcing people to string out the JULES state into a 1D array

if (Ny_model /= 1) then
   write(string1,*)'Only supporting case where "y" dimension is 1.'
   write(string2,*)trim(jules_output_filename)//' has a "y" dimension of ',Ny_model
   write(string3,*)'Unable to continue. Stopping.'
   call error_handler(E_ERR,'get_jules_output_dimensions',string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

if (do_output() .and. debug > 9) then
   write(logfileunit,*)
   write(logfileunit,*)'get_jules_output_dimensions output follows:'
   write(logfileunit,*)'Nx_model = ',Nx_model
   write(logfileunit,*)'Ny_model = ',Ny_model
   write(logfileunit,*)'Nsoil    = ',Nsoil
   write(logfileunit,*)'Ntile    = ',Ntile
   write(     *     ,*)
   write(     *     ,*)'get_jules_output_dimensions output follows:'
   write(     *     ,*)'Nx_model = ',Nx_model
   write(     *     ,*)'Ny_model = ',Ny_model
   write(     *     ,*)'Nsoil    = ',Nsoil
   write(     *     ,*)'Ntile    = ',Ntile
endif

call nc_check(nf90_close(ncid),'get_jules_output_dimensions', &
            'close '//trim(jules_output_filename) )

return
end subroutine get_jules_output_dimensions


!-----------------------------------------------------------------------
!>
!> Read the dimensions from the history netcdf file.
!> The file name comes from module storage ... namelist.
!>

subroutine get_jules_restart_dimensions()

integer :: dimid
integer :: ncid
integer :: myNsoil, myNtile

! sets module variables Nland, Nscpool, Nsoil, Ntile

call nc_check(nf90_open(trim(jules_restart_filename), nf90_nowrite, ncid), &
            'get_jules_restart_dimensions', &
            'open '//trim(jules_restart_filename))

call nc_check(nf90_inq_dimid(ncid, 'land', dimid), &
            'get_jules_restart_dimensions', &
            'inq_dimid land '//trim(jules_restart_filename))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nland), &
            'get_jules_restart_dimensions', &
            'inquire_dimension land '//trim(jules_restart_filename))

call nc_check(nf90_inq_dimid(ncid, 'scpool', dimid), &
            'get_jules_restart_dimensions', &
            'inq_dimid scpool '//trim(jules_restart_filename))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nscpool), &
            'get_jules_restart_dimensions', &
            'inquire_dimension scpool '//trim(jules_restart_filename))

call nc_check(nf90_inq_dimid(ncid, 'soil', dimid), &
            'get_jules_restart_dimensions', &
            'inq_dimid soil '//trim(jules_restart_filename))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=myNsoil), &
            'get_jules_restart_dimensions', &
            'inquire_dimension soil '//trim(jules_restart_filename))

call nc_check(nf90_inq_dimid(ncid, 'tile', dimid), &
            'get_jules_restart_dimensions', &
            'inq_dimid tile '//trim(jules_restart_filename))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=myNtile), &
            'get_jules_restart_dimensions', &
            'inquire_dimension tile '//trim(jules_restart_filename))

if (myNsoil /= Nsoil) then
   write(string1,*)'Confusion about number of soil levels.'
   write(string2,*)'number of soil levels in output  file is ',Nsoil
   write(string3,*)'number of soil levels in restart file is ',myNsoil
   call error_handler(E_ERR,'get_jules_restart_dimensions',string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

if (myNtile /= Ntile) then
   write(string1,*)'Confusion about number of tiles.'
   write(string2,*)'number of tiles in output  file is ',Ntile
   write(string3,*)'number of tiles in restart file is ',myNtile
   call error_handler(E_ERR,'get_jules_restart_dimensions',string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

if (Nx_model /= Nland) then
   write(string1,*)'Only supporting case where "y" dimension is 1, i.e. nland == nx'
   write(string2,*)trim(jules_restart_filename)//' has a "land" dimension of ',Nx_model
   write(string3,*)trim(jules_output_filename)//' has an "x" dimension of ',Nland
   call error_handler(E_ERR,'get_jules_restart_dimensions',string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif


if (do_output() .and. debug > 9) then
   write(logfileunit,*)
   write(logfileunit,*)'get_jules_restart_dimensions output follows:'
   write(logfileunit,*)'Nland   = ',Nland
   write(logfileunit,*)'Nscpool = ',Nscpool
   write(logfileunit,*)'Nsoil   = ',myNsoil
   write(logfileunit,*)'Ntile   = ',myNtile
   write(     *     ,*)
   write(     *     ,*)'get_jules_restart_dimensions output follows:'
   write(     *     ,*)'Nland   = ',Nland
   write(     *     ,*)'Nscpool = ',Nscpool
   write(     *     ,*)'Nsoil   = ',myNsoil
   write(     *     ,*)'Ntile   = ',Ntile
endif

call nc_check(nf90_close(ncid),'get_jules_restart_dimensions', &
            'close '//trim(jules_restart_filename) )

return
end subroutine get_jules_restart_dimensions


!-----------------------------------------------------------------------
!>
!> Read the model grid dimensions from the JULES output netcdf file.
!> LONGITUDE, LATITUDE ... all have module scope
!>

subroutine Get_Model_Grid(fname)

character(len=*), intent(in) :: fname

integer  :: ncid
integer  :: iland, ix, iy
real(r8) :: total_mismatch_lon
real(r8) :: total_mismatch_lat

call nc_check(nf90_open(trim(fname), nf90_nowrite, ncid), &
        'Get_Model_Grid','open '//trim(fname))

! The lat/lon matrices in the history file have been 'packed' so
! that only land gridcells are represented.
! This means we need arrays to convert between physical lat/lon
! grids and the 'packed' arrays.

allocate(LONGITUDE(Nx_model,Ny_model))
allocate( LATITUDE(Nx_model,Ny_model))
allocate(land_to_ij(Nland))
allocate(ij_to_land(Nlon,Nlat))

call DART_get_var(ncid, 'longitude', LONGITUDE)
call DART_get_var(ncid, 'latitude',  LATITUDE)

! just to make sure we are [0,360] and [-90,90]
where (LONGITUDE <   0.0_r8) LONGITUDE = LONGITUDE + 360.0_r8
where (LONGITUDE > 360.0_r8) LONGITUDE = LONGITUDE - 360.0_r8

if (any(LONGITUDE < 0.0_r8)) then
   write(string1,*)'longitudes in history file variable "lon" still negative.'
   call error_handler(E_ERR,'Get_Model_Grid',string1,source,revision,revdate)
endif

where (LATITUDE < -90.0_r8) LATITUDE = -90.0_r8
where (LATITUDE >  90.0_r8) LATITUDE =  90.0_r8

call nc_check(nf90_close(ncid),'Get_Model_Grid','close '//trim(fname) )

! A little sanity check

if (do_output() .and. debug > 0) then

   write(logfileunit,*)
   write(logfileunit,*)'history_file grid information as interpreted ...'
   write(logfileunit,*)'longitude  range ', minval(LONGITUDE), maxval(LONGITUDE)
   write(logfileunit,*)'latitude   range ', minval( LATITUDE), maxval( LATITUDE)
   write(logfileunit,*)'soillevel  is ', SOILLEVEL
   write(     *     ,*)
   write(     *     ,*)'history_file grid information as interpreted ...'
   write(     *     ,*)'longitude  range ', minval(LONGITUDE), maxval(LONGITUDE)
   write(     *     ,*)'latitude   range ', minval( LATITUDE), maxval(LATITUDE)
   write(     *     ,*)'SOILLEVEL  is ', SOILLEVEL

endif

! Confirm that we strip out the land cells in the same fashion as JULES
! This used to be in a routine called confirm_unwrapping()
! Create the lookup arrays relating land cell to lat/lon gridcell, etc.

total_mismatch_lon = 0.0_r8
total_mismatch_lat = 0.0_r8
iland = 0
do iy = 1,Nlat
do ix = 1,Nlon
   if (land_mask(ix,iy) > 0.0_r8) then

      iland = iland + 1 

      land_tile_fractions(iland,1,:) = physical_tile_fractions(ix,iy,:)

      total_mismatch_lon = total_mismatch_lon + &
                           LONGITUDE(iland,1) - physical_longitudes(ix,iy)
      total_mismatch_lat = total_mismatch_lat + &
                           LATITUDE( iland,1) - physical_latitudes( ix,iy)

      if (do_output() .and. debug > 99) write(*,100) iland, ix, iy, &
          LONGITUDE(iland,1) - physical_longitudes(ix,iy), &
           LATITUDE(iland,1) - physical_latitudes( ix,iy)

      ! construct the lookup tables here since we know the
      ! unrolling will be correct

      land_to_ij(iland)%lonindex = ix
      land_to_ij(iland)%latindex = iy
      ij_to_land(ix,iy)          = iland

   endif
enddo
enddo

 100 format('iland, ix, iy, londiff, latdiff',3(1x,i5),2(1x,f14.9))

if (iland /= Nland) then
   write(string1,*)'Unexpected number of land cells.'
   write(string2,*)'expected ',Nland,' got ',iland
   call error_handler(E_ERR, 'Get_Model_Grid', string1, &
          source, revision, revdate, text2=string2)
endif

if ((total_mismatch_lon > 0.001) .or. (total_mismatch_lat > 0.001)) then
   write(string1,*)'Check mapping between model gridcell layout and physical layout.'
   write(string2,*)'total mismatch of land cell longitudes is ',total_mismatch_lon
   write(string3,*)'total mismatch of land cell latitudes  is ',total_mismatch_lat
   call error_handler(E_ERR, 'Get_Model_Grid', string1, &
          source, revision, revdate, text2=string2, text3=string3)
endif

return
end subroutine Get_Model_Grid


!-----------------------------------------------------------------------
!>
!> This defines the window used for assimilation.
!> all observations +/- half this timestep are assimilated.
!>

function set_model_time_step()

type(time_type) :: set_model_time_step

call error_handler(E_MSG, 'set_model_time_step:', 'FIXME SHAMS routine is not tested')

!> @todo FIXME ... should check to see that time step is attainable given the JULES namelist values.

set_model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end function set_model_time_step


!-----------------------------------------------------------------------
!>
!>  This routine checks the user input against the variables available in the
!>  input netcdf file to see if it is possible to construct the DART state vector
!>  specified by the input.nml:model_nml:variables  variable.
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
!>     only variables marked '.r', 'UPDATE' will be modified for JULES.
!>
!>  The calling code should check to see if the variable exists.
!>

function parse_variable_table() result(ngood)

integer :: ngood
! character variables(:)        is module scope
! character variable_table(:,:) is module scope

integer :: i
character(len=NF90_MAX_NAME) :: varname       ! column 1
character(len=NF90_MAX_NAME) :: dartstr       ! column 2
character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
character(len=NF90_MAX_NAME) :: origin_file   ! column 5
character(len=NF90_MAX_NAME) :: state_or_aux  ! column 6

! This loop just repackages the 1D array of values into a 2D array.
! We can do some miniminal checking along the way.
! Determining which file to check is going to be more complicated.

ngood = 0
MyLoop : do i = 1, max_state_variables

   varname      = trim(variables(num_state_table_columns*i - 5))
   dartstr      = trim(variables(num_state_table_columns*i - 4))
   minvalstring = trim(variables(num_state_table_columns*i - 3))
   maxvalstring = trim(variables(num_state_table_columns*i - 2))
   origin_file  = trim(variables(num_state_table_columns*i - 1))
   state_or_aux = trim(variables(num_state_table_columns*i    ))

   call to_upper(origin_file)
   call to_upper(state_or_aux)

   variable_table(i,VT_VARNAMEINDX) = trim(varname)
   variable_table(i,VT_KINDINDX)    = trim(dartstr)
   variable_table(i,VT_MINVALINDX)  = trim(minvalstring)
   variable_table(i,VT_MAXVALINDX)  = trim(maxvalstring)
   variable_table(i,VT_ORIGININDX)  = trim(origin_file)
   variable_table(i,VT_STATEINDX)   = trim(state_or_aux)

   ! If the first element is empty, we have found the end of the list.
   if ( variable_table(i,1) == ' ' ) exit MyLoop

   ! Any other condition is an error.
   if ( any(variable_table(i,:) == ' ') ) then
      string1 = 'input.nml &model_nml:clm_variables not fully specified'
      string2 = 'must be 6 entries per variable. Last known variable name is'
      string3 = '['//trim(variable_table(i,1))//'] ... (without the [], naturally)'
      call error_handler(E_ERR, 'parse_variable_table', string1, &
         source, revision, revdate, text2=string2, text3=string3)
   endif

   ! Make sure DART kind is valid

   if( get_raw_obs_kind_index(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'parse_variable_table',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector

   if (do_output() .and. debug > 8) then
      write(logfileunit,*)'variable ',i,' is ',trim(variable_table(i,1)), ' ', &
                                               trim(variable_table(i,2)), ' ', &
                                               trim(variable_table(i,3)), ' ', &
                                               trim(variable_table(i,4)), ' ', &
                                               trim(variable_table(i,5)), ' ', &
                                               trim(variable_table(i,6))
      write(     *     ,*)'variable ',i,' is ',trim(variable_table(i,1)), ' ', &
                                               trim(variable_table(i,2)), ' ', &
                                               trim(variable_table(i,3)), ' ', &
                                               trim(variable_table(i,4)), ' ', &
                                               trim(variable_table(i,5)), ' ', &
                                               trim(variable_table(i,6))
   endif

   ngood = ngood + 1
enddo MyLoop

if (ngood == max_state_variables) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'parse_variable_table:',string1,text2=string2)
endif

end function parse_variable_table


!-----------------------------------------------------------------------
!>
!> takes the N-dimensional variable and appends the DART
!> dimensions of 'copy' and 'time'. If the variable initially had a 'time'
!> dimension, it is ignored because (by construction) it is a singleton
!> dimension.
!>

subroutine define_var_dims(ivar, ncid, memberdimid, unlimiteddimid, ndims, dimids)

integer,               intent(in)  :: ivar, ncid, memberdimid, unlimiteddimid
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: dimids

character(len=NF90_MAX_NAME),dimension(NF90_MAX_VAR_DIMS) :: dimnames

integer :: i, mydimid

ndims = 0

DIMLOOP : do i = 1,progvar(ivar)%numdims

   if (progvar(ivar)%dimnames(i) == 'time') cycle DIMLOOP

   call nc_check(nf90_inq_dimid(ncid=ncid, name=progvar(ivar)%dimnames(i), &
                 dimid=mydimid), 'define_var_dims', &
                 'inq_dimid '//trim(progvar(ivar)%dimnames(i)))

   ndims         = ndims + 1
   dimids(ndims) = mydimid
   dimnames(ndims) = progvar(ivar)%dimnames(i)

enddo DIMLOOP

! The last two dimensions are always 'copy' and 'time'
ndims           = ndims + 1
dimids(ndims)   = memberdimid
dimnames(ndims) = 'copy'
ndims           = ndims + 1
dimids(ndims)   = unlimitedDimid
dimnames(ndims) = 'time'

if (do_output() .and. debug > 3) then

   write(logfileunit,*)
   write(logfileunit,*)'define_var_dims knowledge'
   write(logfileunit,*)trim(progvar(ivar)%varname),' repackaging into dimnames: ', &
                       (/ (trim(dimnames(i))//' ',i=1,ndims) /)
   write(logfileunit,*)'thus dimids ',dimids(1:ndims)

   write(     *     ,*)
   write(     *     ,*)'define_var_dims knowledge'
   write(     *     ,*)trim(progvar(ivar)%varname),' repackaging into dimnames: ', &
                       (/ (trim(dimnames(i))//' ',i=1,ndims) /)
   write(     *     ,*)'thus dimids ',dimids(1:ndims)

endif

return
end subroutine define_var_dims


!-----------------------------------------------------------------------
!>
!> finds the index into the progvar structure for the named variable
!>

function findVarIndex(varstring, caller)
character(len=*), intent(in) :: varstring
character(len=*), intent(in) :: caller
integer                      :: findVarIndex

integer :: i

findVarIndex = -1

! Skip to the right variable
VARTYPES : do i = 1,nfields
    findVarIndex = i
    if ( trim(progvar(i)%varname) == varstring) exit VARTYPES
enddo VARTYPES

if (findVarIndex < 1) then
   write(string1,*) trim(caller)//' cannot find "'//trim(varstring)// &
                                  '" in list of DART state variables.'
   call error_handler(E_ERR,'findVarIndex',string1,source,revision,revdate)
endif

end function findVarIndex


!-----------------------------------------------------------------------
!>
!> This fills the ij_to_dart structure.
!> After ij_to_dart is filled, a pair of (physical) lon/lat gridcell indices
!> can be used to address how many - and what - dart state vector components
!> are from that (parent) gridcell.
!>

subroutine determine_parent_gridcells()

integer :: soil_index, tile_index, scpool_index, varindex  ! all dummies
integer :: ij, iunit
integer :: ilon, ilat   ! physical grid indices
integer ::   ix,   iy   ! model grid indices

ij_to_dart(:,:)%n = 0   ! everybody gets nobody

! Since get_state_indices() can return the 'land' cell index ('ix')
! and we have a way to convert the land index back to its parent gridcell
! (using land_to_ij) we can simply traverse the entire model structure once
! and count up how many DART contributions are in each parent gridcell.
! Once we know how many DART contributions are in each gridcell, we can allocate
! space to hold an array that holds the DART indices that came from that gridcell.

do ij = 1,model_size
   varindex  = -1 ! if varindex is negative, get_state_indices will calculate it
   call get_state_indices(ij, ix, iy, soil_index, tile_index, scpool_index, varindex)

   ! remember, iy will always equal 1 ... ix is iland
   ilon = land_to_ij(ix)%lonindex
   ilat = land_to_ij(ix)%latindex

   ij_to_dart(ilon, ilat)%n = ij_to_dart(ilon, ilat)%n + 1
enddo

! Create storage for the list of indices

do ilat = 1,Nlat
do ilon = 1,Nlon
   if ( ij_to_dart(ilon,ilat)%n > 0 ) then
      allocate( ij_to_dart(ilon,ilat)%dartIndex(ij_to_dart(ilon,ilat)%n) )
   endif
enddo
enddo

! Now that we know how many DART elements are in each parent gridcell,
! go back and fill ij_to_dart with the DART indices.

ij_to_dart(:,:)%n = 0   ! everybody gets nobody - again

do ij = 1,model_size
   varindex  = -1 ! if varindex is negative, get_state_indices will calculate it
   call get_state_indices(ij, ix, iy, soil_index, tile_index, scpool_index, varindex)

   ilon = land_to_ij(ix)%lonindex
   ilat = land_to_ij(ix)%latindex

   ij_to_dart(ilon, ilat)%n        = ij_to_dart(ilon, ilat)%n + 1
   ij_to_dart(ilon, ilat)%dartIndex( ij_to_dart(ilon, ilat)%n ) = ij
enddo

! Verification
if (do_output() .and. debug > 1) then

   iunit = open_file('gridcell_lookup_table.txt',form='formatted',action='write')

   do ilat = 1,Nlat
   do ilon = 1,Nlon
      write(iunit,'(''ij_to_dart'',i8,1x,i8,'' has '', i6, '' constituents:'')') &
                   ilon,ilat,ij_to_dart(ilon,ilat)%n
      if (ij_to_dart(ilon,ilat)%n > 0) &
         write(iunit,*) ij_to_dart(ilon,ilat)%dartIndex
   enddo
   enddo

   call close_file(iunit)

endif

return
end subroutine determine_parent_gridcells


!-----------------------------------------------------------------------
!>
!> converts the information in the variable_table
!> to the progvar structure for each variable.
!> If the numerical limit does not apply, it is set to MISSING_R8, even if
!> it is the maximum that does not apply.
!>

subroutine SetVariableAttributes(ivar)

integer, intent(in) :: ivar

integer  :: ios
real(r8) :: minvalue, maxvalue

progvar(ivar)%varname     = trim(variable_table(ivar,VT_VARNAMEINDX))
progvar(ivar)%long_name   = ''
progvar(ivar)%units       = ''
progvar(ivar)%coordinates = ''
progvar(ivar)%dimnames    = ''
progvar(ivar)%dimlens     = 0
progvar(ivar)%numdims     = 0
progvar(ivar)%rank        = 0
progvar(ivar)%maxlevels   = 0
progvar(ivar)%xtype       = 0
progvar(ivar)%varsize     = 0
progvar(ivar)%index1      = 0
progvar(ivar)%indexN      = 0
progvar(ivar)%kind_string = trim(variable_table(ivar,VT_KINDINDX))
progvar(ivar)%dart_kind   = get_raw_obs_kind_index( progvar(ivar)%kind_string )
progvar(ivar)%rangeRestricted   = BOUNDED_NONE
progvar(ivar)%minvalue          = MISSING_R8
progvar(ivar)%maxvalue          = MISSING_R8
progvar(ivar)%spvalINT    = -9999
progvar(ivar)%spvalR4     = -1.e20_r4    ! value in the JULES output file
progvar(ivar)%spvalR8     = -1.e20_r8
progvar(ivar)%missingINT  = MISSING_I
progvar(ivar)%missingR4   = MISSING_R4
progvar(ivar)%missingR8   = MISSING_R8
progvar(ivar)%has_fill_value    = .true.
progvar(ivar)%has_missing_value = .true.
progvar(ivar)%update            = .false.

if (variable_table(ivar,VT_ORIGININDX) == 'RESTART') then
   progvar(ivar)%origin = trim(jules_restart_filename)
else
   variable_table(ivar,VT_ORIGININDX) = 'OUTPUT'
   progvar(ivar)%origin = trim(jules_output_filename)
endif

if ((variable_table(ivar,VT_STATEINDX)  == 'UPDATE') .and. &
    (variable_table(ivar,VT_ORIGININDX) == 'RESTART')) progvar(ivar)%update = .true.

! set the default values

minvalue = MISSING_R8
maxvalue = MISSING_R8
progvar(ivar)%minvalue = MISSING_R8
progvar(ivar)%maxvalue = MISSING_R8

! If the character string can be interpreted as an r8, great.
! If not, there is no value to be used.

read(variable_table(ivar,VT_MINVALINDX),*,iostat=ios) minvalue
if (ios == 0) progvar(ivar)%minvalue = minvalue

read(variable_table(ivar,VT_MAXVALINDX),*,iostat=ios) maxvalue
if (ios == 0) progvar(ivar)%maxvalue = maxvalue

! rangeRestricted == BOUNDED_NONE  == 0 ... unlimited range
! rangeRestricted == BOUNDED_BELOW == 1 ... minimum, but no maximum
! rangeRestricted == BOUNDED_ABOVE == 2 ... maximum, but no minimum
! rangeRestricted == BOUNDED_BOTH  == 3 ... minimum and maximum

if (   (progvar(ivar)%minvalue /= MISSING_R8) .and. &
       (progvar(ivar)%maxvalue /= MISSING_R8) ) then
   progvar(ivar)%rangeRestricted = BOUNDED_BOTH

elseif (progvar(ivar)%maxvalue /= MISSING_R8) then
   progvar(ivar)%rangeRestricted = BOUNDED_ABOVE

elseif (progvar(ivar)%minvalue /= MISSING_R8) then
   progvar(ivar)%rangeRestricted = BOUNDED_BELOW

else
   progvar(ivar)%rangeRestricted = BOUNDED_NONE

endif

! Check to make sure min is less than max if both are specified.

if ( progvar(ivar)%rangeRestricted == BOUNDED_BOTH ) then
   if (maxvalue < minvalue) then
      write(string1,*)'&model_nml state_variable input error for ', &
                      trim(progvar(ivar)%varname)
      write(string2,*)'minimum value (',minvalue,') must be less than '
      write(string3,*)'maximum value (',maxvalue,')'
      call error_handler(E_ERR,'SetVariableAttributes',string1, &
         source,revision,revdate,text2=trim(string2),text3=trim(string3))
   endif
endif

return
end subroutine SetVariableAttributes


!-----------------------------------------------------------------------
!>
!> returns the index into the time array that matches
!> the model_time from the JULES restart file.
!>

function FindDesiredTimeIndx(ncid, ntimes, varname)

integer,          intent(in) :: ncid
integer,          intent(in) :: ntimes
character(len=*), intent(in) :: varname
integer                      :: FindDesiredTimeIndx

integer :: VarID
real(r8),  dimension(ntimes) :: mytimes
character(len=NF90_MAX_NAME) :: attvalue

type(time_type) :: thistime, basetime
integer :: ios, itime
integer :: iyear, imonth, iday, ihour, imin, isec

FindDesiredTimeIndx = MISSING_I   ! initialize to failure setting

call nc_check(nf90_inq_varid(ncid, 'time', VarID), &
        'FindDesiredTimeIndx:', 'inq_varid time '//varname)
call nc_check(nf90_get_var(  ncid, VarID, mytimes), &
        'FindDesiredTimeIndx:', 'get_var   time '//varname)
call nc_check(nf90_get_att(ncid, VarID, 'units', attvalue), &
        'FindDesiredTimeIndx:', 'time get_att units '//varname)

! time:units = "seconds since 2004-01-01 00:00:00" ;
!               12345678901234

if (attvalue(1:13) /= 'seconds since') then
   write(string1,*)'expecting time units of [seconds since ... ]'
   write(string2,*)'read time units of ['//trim(attvalue)//']'
   call error_handler(E_ERR, 'FindDesiredTimeIndx:', string1, &
          source, revision, revdate, text2=string2)
endif

read(attvalue,'(14x,i4,5(1x,i2))',iostat=ios)iyear,imonth,iday,ihour,imin,isec
if (ios /= 0) then
   write(string1,*)'Unable to read time units. Error status was ',ios
   write(string2,*)'expected "seconds since YYYY-MM-DD HH:MM:SS"'
   write(string3,*)'was      "'//trim(attvalue)//'"'
   call error_handler(E_ERR, 'FindDesiredTimeIndx:', string1, &
          source, revision, revdate, text2=string2, text3=string3)
endif

basetime = set_date(iyear, imonth, iday, ihour, imin, isec)

! convert each time to a DART time and compare to desired

TIMELOOP : do itime = 1,ntimes

   isec     = int(mytimes(itime))
   thistime = set_time(isec, 0) + basetime

   if (thistime == model_time) then
      FindDesiredTimeIndx = itime
      exit TIMELOOP
   endif

enddo TIMELOOP

!> @todo FIXME ... do we actually need a perfect match ... or do we just use the last one
if ( FindDesiredTimeIndx == MISSING_I ) then
   call print_time(model_time,str='model time is ',iunit=logfileunit)
   call print_time(model_time,str='model time is ')
   call print_date(model_time,str='model date is ',iunit=logfileunit)
   call print_date(model_time,str='model date is ')
   write(string1,*)'No time matching model_time found for '//trim(varname)
   call error_handler(E_ERR, 'FindDesiredTimeIndx:', string1, &
          source, revision, revdate )
endif

if (do_output() .and. debug > 0) then
   write(string1,*)trim(varname)//' matching time index is ',FindDesiredTimeIndx
   call error_handler(E_MSG, 'FindDesiredTimeIndx:', string1)
endif

end function FindDesiredTimeIndx


!-----------------------------------------------------------------------
!>
!> reads the namelist specifying the number of soil layers and the
!> associated layer thicknesses. DART wants the soil layer depth,
!> so the thicknesses are converted to depths.
!> The soil level values are only available from the namelist, apparently.
!>

subroutine get_soil_levels()

! integer,  intent(out) :: Nsoil        ! were it not for module scope
! real(r8), intent(out) :: SOILLEVEL(:) ! were it not for module scope

integer :: iunit, io, i
integer, PARAMETER :: MAX_SOIL_LEVELS = 200

real(r8) :: confrac
real(r8) :: dzsoil_io(MAX_SOIL_LEVELS)
real(r8) :: soil_interfaces(Nsoil + 1)
logical :: l_bedrock
logical :: l_dpsids_dsdz
logical :: l_soil_sat_down
logical :: l_vg_soil
integer :: sm_levels
integer :: soilhc_method
real(r8) :: zsmc
real(r8) :: zst

namelist /jules_soil/ confrac, dzsoil_io, l_bedrock, l_dpsids_dsdz, &
    l_soil_sat_down, l_vg_soil, sm_levels, soilhc_method, zsmc, zst

call find_namelist_in_file('jules_soil.nml', 'jules_soil', iunit)
read(iunit, nml = jules_soil, iostat = io)
call check_namelist_read(iunit, io, 'jules_soil')

if (Nsoil /= sm_levels) then
   write(string1,*)'Number of soil layers unknown.'
   write(string2,*)'Number of soil layers from namelist is ',sm_levels
   write(string3,*)'Number of soil layers from netCDF   is ',Nsoil
   call error_handler(E_ERR, 'get_soil_levels:', string1, &
          source, revision, revdate, text2=string2, text3=string3 )
endif

allocate(SOILLEVEL(Nsoil))

! We want the midpoints of the soil layers.
! The soil temperature in layer 1 is the average of that layer.
! If the first thickness is 10cm, it is more accurate to
! say the temperature 'depth' is 5cm, for example.
! Could use some higher order function of temperature with depth, but ...

soil_interfaces(1) = 0.0_r8
do i = 2,Nsoil+1
   soil_interfaces(i) = soil_interfaces(i-1) + dzsoil_io(i-1)
enddo

do i = 1,Nsoil
   SOILLEVEL(i) = (soil_interfaces(i) + soil_interfaces(i+1)) / 2.0_r8
enddo

return
end subroutine get_soil_levels


!-----------------------------------------------------------------------
!>
!> reads the latitude and longitude matrices for the physical grid
!>
!>         float latitude(y, x) ;
!>                latitude:units = " " ;
!>                latitude:_FillValue = -9999.f ;
!>        float longitude(y, x) ;
!>                longitude:units = " " ;
!>                longitude:_FillValue = -9999.f ;
!>
!> in this context, 'x' and 'y' are 'Nlon' and 'Nlat', respectively.
!> the model sizes are called Nx_model, Ny_model by DART.
!>

subroutine set_physical_grid(filename)

character(len=*), intent(in) :: filename

integer :: ncid, i, j
integer :: VarID, numdims
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens

if ( .not. file_exist(filename) ) then
   write(string1,*)'cannot open file <', trim(filename),'> to read physical locations.'
   write(string2,*)'filename came from model_grid.nml &jules_latlon "file"'
   call error_handler(E_ERR, 'set_physical_grid', string1, &
       source, revision, revdate, text2=string2)
endif

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncid), &
        'set_physical_grid', 'open '//trim(filename))

! Create a convenient string for messages
string3 = trim(filename)//' '//trim(physical_lon_variable)

! The physical grid dimensions are based on the shape of the longitude variable.

call nc_check(nf90_inq_varid(ncid, trim(physical_lon_variable), VarID), &
         'set_physical_grid', 'inq_varid lon_name '//trim(string3))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
         'set_physical_grid', 'inquire '//trim(string3))

if (numdims /= 2) then
   write(string1,*)trim(string3)//' expected to have 2 dimensions, has ',numdims
   write(string2,*)'unable to proceed.'
   call error_handler(E_ERR, 'set_physical_grid', string1, &
          source, revision, revdate, text2=string2)
endif

do i = 1,numdims
   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
      'set_physical_grid', 'inquire dimension '//trim(string3))
enddo

Nlon = dimlens(1) ! set global variable for number of longitudes of physical grid
Nlat = dimlens(2) ! set global variable for number of latitudes  of physical grid

allocate(physical_longitudes(Nlon,Nlat))
allocate(physical_latitudes( Nlon,Nlat))

! Actually get the longitude variable

call DART_get_var(ncid, trim(physical_lon_variable), physical_longitudes, string3)

where (physical_longitudes < 0.0_r8) &
       physical_longitudes = physical_longitudes + 360.0_r8

! Get the latitude variable

string3 = trim(filename)//' '//trim(physical_lat_variable)

call nc_check(nf90_inq_varid(ncid, trim(physical_lat_variable), VarID), &
         'set_physical_grid', 'inq_varid lat_name '//trim(string3))

call DART_get_var(ncid, trim(physical_lat_variable), physical_latitudes, string3)

! Clean up

call nc_check(nf90_close(ncid),'set_physical_grid','close '//trim(filename))

if (do_output() .and. debug > 99) then
   write(logfileunit,*)
   write(logfileunit,*)'Summary of lat/lon in '//trim(filename)
   write(logfileunit,*)'lat shape      ', shape(physical_latitudes)
   write(logfileunit,*)'lat minval     ',minval(physical_latitudes)
   write(logfileunit,*)'lat maxval     ',maxval(physical_latitudes)
   write(logfileunit,*)'lon shape      ', shape(physical_longitudes)
   write(logfileunit,*)'lon minval     ',minval(physical_longitudes)
   write(logfileunit,*)'lon maxval     ',maxval(physical_longitudes)
   write(     *     ,*)
   write(     *     ,*)'Summary of lat/lon in '//trim(filename)
   write(     *     ,*)'lat shape      ', shape(physical_latitudes)
   write(     *     ,*)'lat minval     ',minval(physical_latitudes)
   write(     *     ,*)'lat maxval     ',maxval(physical_latitudes)
   write(     *     ,*)'lon shape      ', shape(physical_longitudes)
   write(     *     ,*)'lon minval     ',minval(physical_longitudes)
   write(     *     ,*)'lon maxval     ',maxval(physical_longitudes)
endif

if (do_output() .and. debug > 99) then
   write(*,*)'checking physical_longitude(i,j),physical_latitude(i,j)'
   do j = 1,Nlat
   do i = 1,Nlon
      write(*,'(''('',i1,1x,i1,'') ='',2(1x,f14.6) )'), i, j, &
                physical_longitudes(i,j), physical_latitudes(i,j)
   enddo
   enddo
endif

return
end subroutine set_physical_grid


!-----------------------------------------------------------------------
!>
!> reads the land fractions for each physical grid cell
!>

subroutine set_land_mask(filename)

character(len=*), intent(in) :: filename

integer :: ncid, i
integer :: VarID, numdims
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens

! the file that contains the land fraction must also contain the
! latitude and longitude matrices of the physical grid, not the grid
! that just contains the 'land' gridcells. As such, the 'full' shape of
! the land fraction, latitude, and longitude matrices 

allocate(land_mask(Nlon,Nlat))

if ( trim(filename) == 'none' ) then ! implicitly, all land.
   land_mask = 1.0_r8
   return
endif

if ( .not. file_exist(filename) ) then
   write(string1,*)'cannot open file <', trim(filename),'> to read land fractions.'
   write(string2,*)'filename came from model_grid.nml &jules_land_frac "file"'
   call error_handler(E_ERR, 'set_land_mask', string1, &
       source, revision, revdate, text2=string2)
endif

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncid), &
        'set_land_mask', 'open '//trim(filename))

! Create a convenient string for messages
string3 = trim(filename)//' '//trim(landfraction_variable)

! The land fractions come from the input grid, which is always 'physical' 
! (x,y are 'full grid'). The only way to relate which physical-space tile fraction
! to the model-space tile fractions is to skip the non-land physical-space gridcells.
! The packing order is longitude-then-latitude. 

call nc_check(nf90_inq_varid(ncid, trim(landfraction_variable), VarID), &
         'set_land_mask', 'inq_varid '//trim(string3))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
         'set_land_mask', 'inquire '//trim(string3))

if (numdims /= 2) then
   write(string1,*)trim(string3)//' expected to have 2 dimensions, has ',numdims
   write(string2,*)'unable to proceed.'
   call error_handler(E_ERR, 'set_land_mask', string1, &
          source, revision, revdate, text2=string2)
endif

do i = 1,numdims
   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
      'set_land_mask', 'inquire dimension '//trim(string3))
enddo

! Check to see if land mask is the same shape as the physical lat & lon arrays

if (dimlens(1) /= Nlon) then
   write(string1,*)trim(string3)//' expected to have ',Nlon,' longitudes.'
   write(string2,*)'Found ',dimlens(1),'. Unable to proceed.'
   call error_handler(E_ERR, 'set_land_mask', string1, &
          source, revision, revdate, text2=string2)
endif

if (dimlens(2) /= Nlat) then
   write(string1,*)trim(string3)//' expected to have ',Nlat,' latitudes.'
   write(string2,*)'Found ',dimlens(2),'. Unable to proceed.'
   call error_handler(E_ERR, 'set_land_mask', string1, &
          source, revision, revdate, text2=string2)
endif

! DART_get_var() is responsible for making sure variable is expected shape.

call DART_get_var(ncid, trim(landfraction_variable), land_mask, string3)

call nc_check(nf90_close(ncid),'set_land_mask','close '//trim(filename))

if (do_output() .and. debug > 99) then
   write(logfileunit,*)
   write(logfileunit,*)'set_land_mask:Summary of '//trim(string3)
   write(logfileunit,*)'set_land_mask:shape      ',dimlens(1:numdims)
   
   write(     *     ,*)
   write(     *     ,*)'set_land_mask:Summary of '//trim(string3)
   write(     *     ,*)'set_land_mask:shape      ',dimlens(1:numdims)

   !> @todo FIXME ... report on number of non-land & land gridcells
endif

return
end subroutine set_land_mask


!-----------------------------------------------------------------------
!>
!> reads the namelist specifying the name of the file containing the 
!> tile fractions for each gridcell, and then reads the tile fractions.
!>

subroutine set_tile_fractions(filename)

character(len=*), intent(in) :: filename

integer :: ncid, i
integer :: VarID, numdims, dimid
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens

integer :: iland, ix, iy

! Create a convenient string for messages
string3 = trim(filename)//' '//trim(tilefraction_variable)

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncid), &
             'set_tile_fractions','open '//trim(filename))

! Depending on some namelist settings, the shape of the tile fractions may be
! 'physical' (x,y are 'full grid') or 'sparse' (x=nland,y=1). If the restart/output 
! files are 'sparse', then the only way to relate which physical-space tile fraction
! should be used is to read the physical-space land mask.
! The ordering of the land cells is longitude-then-latitude. 

call nc_check(nf90_inq_varid(ncid, trim(tilefraction_variable), VarID), &
         'set_tile_fractions', 'inq_varid '//trim(string3))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
         'set_tile_fractions', 'inquire '//trim(string3))

do i = 1,numdims
   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
      'set_tile_fractions', 'inquire dimension '//trim(string3))
enddo

if ( numdims /= 3) then
   write(string1,*)'<', trim(string3),'> has shape',dimlens(1:numdims)
   write(string2,*)'Expected 3 dimensions. Unable to proceed.'
   call error_handler(E_ERR, 'set_tile_fractions', string1, &
       source, revision, revdate, text2=string2)
endif

if (dimlens(3) /= Ntile) then
   write(string1,*)trim(string3)//' has wrong number of tiles (dimension 3).'
   write(string2,*)'expected ',Ntile,' got ',dimlens(3)
   call error_handler(E_ERR, 'set_tile_fractions', string1, &
          source, revision, revdate, text2=string2)
endif

allocate(physical_tile_fractions(Nlon, Nlat, Ntile))
allocate(land_tile_fractions(Nx_model, Ny_model, Ntile))

call DART_get_var(ncid, trim(tilefraction_variable), physical_tile_fractions, string3)

call nc_check(nf90_inq_dimid(ncid, 'pft', dimid), &
            'set_tile_fractions', 'inq_dimid pft '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Npft), &
            'set_tile_fractions', 'inquire_dimension pft '//trim(filename))

call nc_check(nf90_close(ncid),'set_tile_fractions','close '//trim(string3))

! confirm_unwrapping() has confirmed that if we store the land cells in this order,
! it is consistent as when JULES squeezes out the non-land cells.

iland = 0
do iy = 1,Nlat
do ix = 1,Nlon
   if (land_mask(ix,iy) > 0.0_r8) then
      iland = iland + 1 
      land_tile_fractions(iland,Ny_model,:) = physical_tile_fractions(ix,iy,:)
   endif
enddo
enddo

if (iland /= Nland) then
   write(string1,*)trim(string3)//' does not have expected number of land cells.'
   write(string2,*)'expected ',Nland,' got ',iland
   call error_handler(E_ERR, 'set_tile_fractions', string1, &
          source, revision, revdate, text2=string2)
endif


if (do_output() .and. debug > 99) then
   write(logfileunit,*)
   write(logfileunit,*)'set_tile_fractions:Summary of '//trim(string3)
   write(logfileunit,*)'set_tile_fractions:shape      ',shape(physical_tile_fractions)
   write(logfileunit,*)'set_tile_fractions:Npft = ',Npft
   
   write(     *     ,*)
   write(     *     ,*)'set_tile_fractions:Summary of '//trim(string3)
   write(     *     ,*)'set_tile_fractions:shape      ',shape(physical_tile_fractions)
   write(     *     ,*)'set_tile_fractions:Npft = ',Npft

   !> @todo FIXME ... report min/max?
endif

return
end subroutine set_tile_fractions


!-----------------------------------------------------------------------
!>
!> Given an integer index into the state vector structure, returns the
!> associated array indices for lat, lon, and depth, as well as the type.
!>

subroutine get_state_indices(index_in, x_index, y_index, soil_index, &
                             tile_index, scpool_index, varindex)

integer, intent(in)    :: index_in
integer, intent(out)   ::      x_index
integer, intent(out)   ::      y_index
integer, intent(out)   ::   soil_index
integer, intent(out)   ::   tile_index
integer, intent(out)   :: scpool_index
integer, intent(inout) :: varindex

integer :: n, offset, ndim1, ndim2
integer :: land_index, indx1, indx2, indx3

! make sure this is valid index
if ((index_in < progvar(1)%index1) .or. &
    (index_in > progvar(nfields)%indexN) ) then
   write(string1,*) 'desired index ',index_in
   write(string2,*) 'is not within the bounds of the DART address space of'
   write(string3,*) progvar(1)%index1, progvar(nfields)%indexN
   call error_handler(E_ERR,'get_state_indices',string1, &
         source,revision,revdate,text2=string2,text3=string3)
endif

! If need be, find the DART variable that spans the index of interest.
if (varindex < 0) then
   FindIndex : do n = 1,nfields
      if( (progvar(n)%index1 <= index_in) .and. (index_in <= progvar(n)%indexN) ) then
         varindex = n
         exit FindIndex
      endif
   enddo FindIndex
elseif ((index_in < progvar(varindex)%index1) .or. &
        (index_in > progvar(varindex)%indexN) ) then
   write(string1,*) 'desired index ',index_in
   write(string2,*) 'is not within the bounds of the DART address space for ', &
                     trim(progvar(varindex)%varname)
   write(string3,*) 'must be between',progvar(varindex)%index1, 'and ', &
                                      progvar(varindex)%indexN
   call error_handler(E_ERR,'get_state_indices',string1, &
         source,revision,revdate,text2=string2,text3=string3)
endif

     x_index = -1
     y_index = -1
  soil_index = -1
  tile_index = -1
scpool_index = -1
  land_index = -1

offset = index_in - progvar(varindex)%index1 + 1

! So now we know the variable and its shape. Remember that DART only
! stores a single timestep, so any time dimension is a singleton.
! Furthermore, netCDF requires the time/unlimited dimension be the LAST
! (in Fortran) dimension, so we can just focus on the first N dimensions.
! Relying on integer arithmetic.

if     ( progvar(varindex)%rank == 1) then
   ! This can only be dimensioned 'land' ... possible FIXME ... check

   if (Ny_model /= 1) then ! the X-Y plane needs to be further decomposed

      ! If Ny_model /= 1 then Nx_model * Ny_model is known to equal Nland
      ! So if I have the index into 'Nland', I should be able
      ! to unwrap the 'offset' into Nx_model,Ny_model with the 2D algorithm

      ndim1     = Nx_model
      y_index = 1 + (offset - 1)/ndim1
      x_index = offset - (y_index - 1)*ndim1

   else

      x_index = offset
      y_index = 1

   endif

elseif ( progvar(varindex)%rank == 2) then

   ! Nland always equals Nx_model*Ny_model. Always.
   ! The lon,lat arrays are ALWAYS (Nx_model,Ny_model)
   ! Sometimes Ny_model == 1 (in which case Nx_model == Nland)
   ! So when (for example) Nland = 679, the lat,lon arrays are (679,1)
   ! But when Nland = 800 (all cells land), the lat,lon arrays could be (20,40)
   ! or (40,20) ... So I can correctly get the land_index, but I need to
   ! decompose it into x_index, y_index 
   !
   ! (time does not count toward 'rank') variables like 
   ! float latent_heat([time,] y,    x) 
   ! float      t_soil(     soil, land)
   ! float      canopy(     tile, land)
   ! float          cs(   scpool, land)
   ! Remember that the order in Fortran is reversed from the ncdump order.

   ndim1 = progvar(varindex)%dimlens(1)

   indx2 = 1 + (offset - 1)/ndim1
   indx1 = offset - (indx2 - 1)*ndim1

   ! Depending on the variable, the
   ! 'indx2' is actually either 'soil_index' or 'tile_index' or ...
   ! and 'indx1' is the index into the X-Y plane

   if     (trim(progvar(varindex)%dimnames(2)) ==   'tile') then
      land_index   = indx1
      tile_index   = indx2
   elseif (trim(progvar(varindex)%dimnames(2)) ==   'soil') then
      land_index   = indx1
      soil_index   = indx2
   elseif (trim(progvar(varindex)%dimnames(2)) == 'scpool') then
      land_index   = indx1
      scpool_index = indx2
   else
      x_index    = indx1
      y_index    = indx2
   endif

   if (Ny_model /= 1) then ! the X-Y plane needs to be further decomposed

      ! If Ny_model /= 1 then Nx_model * Ny_model is known to equal Nland
      ! So if I have the index into 'Nland', I should be able
      ! to unwrap the 'land_index' into Nx_model,Ny_model with the same algorithm

      ndim1     = Nx_model
      offset    = land_index

      y_index = 1 + (offset - 1)/ndim1
      x_index = offset - (y_index - 1)*ndim1

   else

      x_index = indx1
      y_index = 1

   endif

elseif ( progvar(varindex)%rank == 3) then

   ! (time does not count toward 'rank') variables like 
   ! float soil_wet(time, soil, y, x)
   ! float    esoil(time, tile, y, x)
   ! Remember that the order in Fortran is reversed from the ncdump order.
   !
   ! Since none of the rank 3 variables have 'land' as the fastest dimension
   ! this case is easier than the rank 2 case.

   ndim1 = progvar(varindex)%dimlens(1)  ! num_fastest_dimension (Fortran leftmost)
   ndim2 = progvar(varindex)%dimlens(2)  ! num_next_fastest

   indx3 = 1 + (offset - 1)/(ndim1*ndim2)
   indx2 = 1 + (offset - (indx3-1)*ndim1*ndim2 -1)/ndim1
   indx1 = offset - (indx3-1)*ndim1*ndim2 - (indx2-1)*ndim1

   ! Depending on the variable, the
   ! 'indx3' is actually either 'soil_index' or 'tile_index' or ...
   ! and 'indx1' is the index into the X-Y plane

   if     (trim(progvar(varindex)%dimnames(3)) ==   'tile') then
      tile_index   = indx3
   elseif (trim(progvar(varindex)%dimnames(3)) ==   'soil') then
      soil_index   = indx3
   elseif (trim(progvar(varindex)%dimnames(3)) == 'scpool') then
      scpool_index = indx3
   else
      write(string1,*) 'No support for coordinates ', &
                       trim(progvar(varindex)%coordinates)
      write(string2,*) 'variable is ',trim(progvar(varindex)%varname)
      call error_handler(E_ERR,'get_state_indices:',string1, &
            source,revision,revdate,text2=string2)
   endif

   x_index = indx1
   y_index = indx2

else
   write(string1,*) 'Does not support variables with rank ',progvar(varindex)%rank
   write(string2,*) 'variable     is ',trim(progvar(varindex)%varname)
   write(string3,*) 'coordinates are ',trim(progvar(varindex)%coordinates)
   call error_handler(E_ERR,'get_state_indices:',string1, &
         source, revision, revdate, text2=string2, text3=string3)
endif

if (do_output() .and. debug > 9) then
   print *, 'get_state_indices: asking for meta data about index ', index_in
   print *, 'get_state_indices: variable      ', trim(progvar(varindex)%varname)
   print *, 'get_state_indices: coordinates   ', trim(progvar(varindex)%coordinates)
   print *, 'get_state_indices: rank          ', progvar(varindex)%rank
   print *, 'get_state_indices: dimensions    ', progvar(varindex)%dimlens(1:progvar(varindex)%numdims)
   print *, 'get_state_indices: Nland, Nx_model, Ny_model ', Nland, Nx_model, Ny_model
   print *, 'get_state_indices: x      index = ',    x_index, 'of ',Nx_model
   print *, 'get_state_indices: y      index = ',    y_index, 'of ',Ny_model
   print *, 'get_state_indices: soil   index = ',   soil_index
   print *, 'get_state_indices: tile   index = ',   tile_index
   print *, 'get_state_indices: scpool index = ', scpool_index
   print *
endif

return
end subroutine get_state_indices


!-----------------------------------------------------------------------
!>
!> Given an integer index into the state vector structure,
!> return the longitude, latitude, and level
!>

subroutine get_state_lonlatlev(varindex, x_index, y_index, soil_index, &
                               mylon, mylat, mylev, mycoordsystem)

integer,  intent(in)  :: varindex
integer,  intent(in)  :: x_index
integer,  intent(in)  :: y_index
integer,  intent(in)  :: soil_index
real(r8), intent(out) :: mylon
real(r8), intent(out) :: mylat
real(r8), intent(out) :: mylev
integer,  intent(out) :: mycoordsystem

if     ((trim(progvar(varindex)%dimnames(1)) == 'x') .and. &
        (trim(progvar(varindex)%dimnames(2)) == 'y') .and. &
        (trim(progvar(varindex)%dimnames(3)) == 'soil')) then
   mylon =  LONGITUDE(x_index, y_index)
   mylat =   LATITUDE(x_index, y_index)
   mylev =  SOILLEVEL(soil_index)
   mycoordsystem = VERTISHEIGHT

elseif ((trim(progvar(varindex)%dimnames(1)) == 'land') .and. &
        (trim(progvar(varindex)%dimnames(2)) == 'soil')) then
   mylon = LONGITUDE(x_index, y_index)
   mylat =  LATITUDE(x_index, y_index)
   mylev =  SOILLEVEL(soil_index)
   mycoordsystem = VERTISHEIGHT

elseif ((trim(progvar(varindex)%dimnames(1)) == 'x') .and. &
        (trim(progvar(varindex)%dimnames(2)) == 'y')) then
   mylon = LONGITUDE(x_index, y_index)
   mylat =  LATITUDE(x_index, y_index)
   mylev = 0.0_r8
   mycoordsystem = VERTISSURFACE

elseif ((trim(progvar(varindex)%dimnames(1)) == 'land') .and. &
        (trim(progvar(varindex)%dimnames(2)) == 'tile')) then
   mylon = LONGITUDE(x_index, y_index)
   mylat =  LATITUDE(x_index, y_index)
   mylev = 0.0_r8
   mycoordsystem = VERTISSURFACE

elseif ((trim(progvar(varindex)%dimnames(1)) == 'land') .and. &
        (trim(progvar(varindex)%dimnames(2)) == 'scpool')) then
   mylon = LONGITUDE(x_index, y_index)
   mylat =  LATITUDE(x_index, y_index)
   mylev = 0.0_r8
   mycoordsystem = VERTISSURFACE

elseif (trim(progvar(varindex)%coordinates) == 'land') then
   mylon = LONGITUDE(x_index, y_index)
   mylat =  LATITUDE(x_index, y_index)
   mylev = 0.0_r8
   mycoordsystem = VERTISSURFACE

else

   !> @todo FIXME ... should check for supported variable shapes much earlier on.

   write(string1,*) 'unsupported variable shape of ['// &
                    &trim(progvar(varindex)%coordinates)//']'
   write(string2,*) 'for variable ',trim(progvar(varindex)%varname)
   call error_handler(E_ERR,'get_state_lonlatlev',string1, &
         source,revision,revdate,text2=string2)

endif

return
end subroutine get_state_lonlatlev


!-----------------------------------------------------------------------
!>
!> set_sparse_locations_kinds() creates an array the size of the physical
!> grid with the location of each of the corresponding gridcell centers.
!>

subroutine set_sparse_locations_kinds()

integer  :: ix, iy, indx
real(r8) :: mylon, mylat

allocate(statespace_locations(Nlon*Nlat))
allocate(statespace_kinds(Nlon*Nlat))

statespace_kinds = 1   ! This is a dummy because we have to satisfy get_close_obs() 
                 ! get_close_obs() requires a kind, but this implementation does not.

indx = 0
do iy = 1,Nlat
do ix = 1,Nlon

    indx  = indx + 1
    mylon = physical_longitudes(ix,iy)
    mylat = physical_latitudes( ix,iy)

    statespace_locations(indx) = set_location(mylon, mylat, 0.0_r8, VERTISUNDEF)
enddo
enddo

return
end subroutine set_sparse_locations_kinds


!-----------------------------------------------------------------------
!>
!> init_interp() Initializes data structures needed for interpolation.
!>
!> General philosophy is to precompute a 'get_close' structure with
!> a list of what state vector elements are within a certain distance
!> of each other.
!>

subroutine init_interp()

! This should be called at static_init_model time to avoid
! having all this temporary storage in the middle of a run.
!
! DART has a 'get_close' type that divides the domain into a set of boxes
! and tracks which locations are in which box. The get_close_obs() routine
! then finds the box and the subsequent distances. The trick is to set up
! the set of boxes parsimoniously.

integer :: i, j

real(r8) :: max_lon, max_lat, maxdist
real(r8) :: lon_dist, lat_dist

! Need to determine the likely maximum size of any of the gridcells.
! the idea is to find a distance that allows us to find the surrounding
! gridcell locations without finding 'too many' neighboring gridcells.
!
! Since this may be a point, regional, or global model, I am going 
! to use the lat & lon values as proxies for the physical distance.
! The maximum grid cell size from one grid should suffice.

max_lon = 0.0_r8
max_lat = 0.0_r8

EW_lat : do j = 1, Nlat
EW_lon : do i = 1, Nlon-1

   lat_dist = abs(physical_latitudes( i,j) - physical_latitudes( i+1,j))
   lon_dist = abs(physical_longitudes(i,j) - physical_longitudes(i+1,j))

   !> @todo FIXME the prime meridian is a problem ... 
   ! hopefully, there are other gridcells 
   if (lon_dist > 180.0_r8) cycle EW_lon

   if(lon_dist > max_lon) max_lon = lon_dist
   if(lat_dist > max_lat) max_lat = lat_dist

enddo EW_lon
enddo EW_lat

NS_lat : do j = 1, Nlat-1
NS_lon : do i = 1, Nlon

   lat_dist = abs(physical_latitudes( i,j) - physical_latitudes( i,j+1))
   lon_dist = abs(physical_longitudes(i,j) - physical_longitudes(i,j+1))

   !> @todo FIXME the prime meridian is a problem ... 
   ! hopefully, there are other gridcells 
   if (lon_dist > 180.0_r8) cycle NS_lon

   if(lon_dist > max_lon) max_lon = lon_dist
   if(lat_dist > max_lat) max_lat = lat_dist
enddo NS_lon
enddo NS_lat

!> @todo FIXME ... what to do if max_lon = 0.0_r8

! convert to radians
max_lon = max_lon / rad2deg
max_lat = max_lat / rad2deg
maxdist = sqrt(max_lon*max_lon + max_lat*max_lat)

if (do_output() .and. debug > 1) then
   write(*,*)
   write(*,*)'init_interp: summary'
   write(*,*)'init_interp: Maximum distance ',maxdist, ' (radians) based on '
   write(*,*)'init_interp: max_dlon = ', max_lon, 'radians', max_lon*rad2deg, 'degrees'
   write(*,*)'init_interp: max_dlat = ', max_lat, 'radians', max_lat*rad2deg, 'degrees'
   write(*,*)
endif

! maxdist unscaled resulted in NNN candidates for a location in the middle of the domain

call get_close_maxdist_init(gc_state, maxdist)
call get_close_obs_init(gc_state, Nlon*Nlat, statespace_locations)

return
end subroutine init_interp


!-----------------------------------------------------------------------
!>
!> Simply prints a summary of everything in the progvar structure.
!>
!> General philosophy is to precompute a 'get_close' structure with
!> a list of what state vector elements are within a certain distance
!> of each other.
!>

subroutine dump_structure(ivar)

integer, intent(in) :: ivar

integer :: i

write(logfileunit,*)
write(logfileunit,*) trim(progvar(ivar)%varname),' variable number ',ivar
write(logfileunit,*) '  filename    :',trim(progvar(ivar)%origin)
write(logfileunit,*) '  update      :',progvar(ivar)%update
write(logfileunit,*) '  long_name   :',trim(progvar(ivar)%long_name)
write(logfileunit,*) '  units       :',trim(progvar(ivar)%units)
write(logfileunit,*) '  xtype       :',progvar(ivar)%xtype
write(logfileunit,*) '  rank        :',progvar(ivar)%rank
write(logfileunit,*) '  numdims     :',progvar(ivar)%numdims
write(logfileunit,*) '  coordinates :',trim(progvar(ivar)%coordinates)

do i = 1,progvar(ivar)%numdims
   write(logfileunit,'(''   dimension ('',i1,'') length '',i10,'' name '',A)') &
              i,progvar(ivar)%dimlens(i),trim(progvar(ivar)%dimnames(i))
enddo

write(logfileunit,*) '  varsize     :',progvar(ivar)%varsize
write(logfileunit,*) '  index1      :',progvar(ivar)%index1
write(logfileunit,*) '  indexN      :',progvar(ivar)%indexN
write(logfileunit,*) '  dart_kind   :',progvar(ivar)%dart_kind
write(logfileunit,*) '  kind_string :',progvar(ivar)%kind_string
write(logfileunit,*) '  spvalINT    :',progvar(ivar)%spvalINT
write(logfileunit,*) '  spvalR4     :',progvar(ivar)%spvalR4
write(logfileunit,*) '  spvalR8     :',progvar(ivar)%spvalR8
write(logfileunit,*) '  missingINT  :',progvar(ivar)%missingINT
write(logfileunit,*) '  missingR4   :',progvar(ivar)%missingR4
write(logfileunit,*) '  missingR8   :',progvar(ivar)%missingR8
write(logfileunit,*) '  has_fill_value    :',progvar(ivar)%has_fill_value
write(logfileunit,*) '  has_missing_value :',progvar(ivar)%has_missing_value
write(logfileunit,*)'   rangeRestricted   :',progvar(ivar)%rangeRestricted
write(logfileunit,*)'   minvalue          :',progvar(ivar)%minvalue
write(logfileunit,*)'   maxvalue          :',progvar(ivar)%maxvalue

write(     *     ,*)
write(     *     ,*) trim(progvar(ivar)%varname),' variable number ',ivar
write(     *     ,*) '  filename    :',trim(progvar(ivar)%origin)
write(     *     ,*) '  update      :',progvar(ivar)%update
write(     *     ,*) '  long_name   :',trim(progvar(ivar)%long_name)
write(     *     ,*) '  units       :',trim(progvar(ivar)%units)
write(     *     ,*) '  xtype       :',progvar(ivar)%xtype
write(     *     ,*) '  rank        :',progvar(ivar)%rank
write(     *     ,*) '  numdims     :',progvar(ivar)%numdims
write(     *     ,*) '  coordinates :',trim(progvar(ivar)%coordinates)

do i = 1,progvar(ivar)%numdims
   write(  *,'(''   dimension ('',i1,'') length '',i10,'' name '',A)') &
              i,progvar(ivar)%dimlens(i),trim(progvar(ivar)%dimnames(i))
enddo

write(     *     ,*) '  varsize     :',progvar(ivar)%varsize
write(     *     ,*) '  index1      :',progvar(ivar)%index1
write(     *     ,*) '  indexN      :',progvar(ivar)%indexN
write(     *     ,*) '  dart_kind   :',progvar(ivar)%dart_kind
write(     *     ,*) '  kind_string :',progvar(ivar)%kind_string
write(     *     ,*) '  spvalINT    :',progvar(ivar)%spvalINT
write(     *     ,*) '  spvalR4     :',progvar(ivar)%spvalR4
write(     *     ,*) '  spvalR8     :',progvar(ivar)%spvalR8
write(     *     ,*) '  missingINT  :',progvar(ivar)%missingINT
write(     *     ,*) '  missingR4   :',progvar(ivar)%missingR4
write(     *     ,*) '  missingR8   :',progvar(ivar)%missingR8
write(     *     ,*) '  has_fill_value    :',progvar(ivar)%has_fill_value
write(     *     ,*) '  has_missing_value :',progvar(ivar)%has_missing_value
write(     *     ,*)'   rangeRestricted   :',progvar(ivar)%rangeRestricted
write(     *     ,*)'   minvalue          :',progvar(ivar)%minvalue
write(     *     ,*)'   maxvalue          :',progvar(ivar)%maxvalue

return
end subroutine dump_structure


!-----------------------------------------------------------------------
!>
!> Get the names of that files that contain the physically-shaped variables. 
!> the land fraction mask,
!> the tile fraction mask, that sort of stuff.
!>

subroutine get_physical_filenames()

! module variables set by this routine:
! physical_latlon_file 
! land_fraction_file 
! landfraction_variable 

integer :: iunit, io

character(len=256) :: file = 'none'
character(len=256) :: lat_name = 'none'
character(len=256) :: lon_name = 'none'
character(len=256) :: land_frac_name = 'none'
character(len=256) :: frac_name = 'none'

namelist /jules_latlon/ file, lat_name, lon_name

namelist /jules_land_frac/ file, land_frac_name

namelist /jules_frac/ file, frac_name

! Read the land fraction filename and variable name first

call find_namelist_in_file('model_grid.nml', 'jules_land_frac', iunit)
read(iunit, nml = jules_land_frac, iostat = io)
call check_namelist_read(iunit, io, 'jules_land_frac')

land_fraction_file = file
landfraction_variable = land_frac_name

if ( (trim(file) /= 'none') .and. (.not. file_exist(file)) ) then
   write(string1,*)'cannot open file <', trim(file),'> to read land fractions.'
   write(string2,*)'filename came from model_grid.nml &jules_land_frac "file"'
   call error_handler(E_ERR, 'get_physical_filenames', string1, &
       source, revision, revdate, text2=string2)
endif

! Read the namelist specifying the physical grid

file = 'none'

call find_namelist_in_file('model_grid.nml', 'jules_latlon', iunit)
read(iunit, nml = jules_latlon, iostat = io)
call check_namelist_read(iunit, io, 'jules_latlon')

if ( trim(file) == 'none' ) then ! single column case, maybe
   write(string1,*)'No latlon file to read input grid latitudes and longitudes.'
   write(string2,*)'filename should come from model_grid.nml &jules_latlon "file"'
   call error_handler(E_ERR, 'get_physical_filenames', string1, &
       source, revision, revdate, text2=string2)
endif

if ( .not. file_exist(file) ) then
   write(string1,*)'cannot open file <', trim(file),'> to read lats and lons.'
   write(string2,*)'filename came from model_grid.nml &jules_latlon "file"'
   call error_handler(E_ERR, 'get_physical_filenames', string1, &
       source, revision, revdate, text2=string2)
endif

physical_latlon_file  = file
physical_lon_variable = lon_name
physical_lat_variable = lat_name

! Read the namelist specifying the tile fractions

file = 'none'

call find_namelist_in_file('ancillaries.nml', 'jules_frac', iunit)
read(iunit, nml = jules_frac, iostat = io)
call check_namelist_read(iunit, io, 'jules_frac')

if ( trim(file) == 'none' ) then
   write(string1,*)'tile fraction file not specified.'
   write(string2,*)'filename should come from ancillaries.nml &jules_frac "file"'
   write(string3,*)'unable to proceed.'
   call error_handler(E_ERR, 'get_physical_filenames', string1, &
       source, revision, revdate, text2=string2, text3=string3)
endif

if ( .not. file_exist(file) ) then
   write(string1,*)'cannot open file <', trim(file),'> to read tile fractions.'
   write(string2,*)'filename came from ancillaries.nml &jules_frac "file"'
   call error_handler(E_ERR, 'get_physical_filenames', string1, &
       source, revision, revdate, text2=string2)
endif

tile_fraction_file  = file
tilefraction_variable = frac_name

! Print a little summary if so desired.

if (do_output() .and. debug > 99) then
   write(     *     ,*)
   write(     *     ,*)'tile_fraction_file     is <'//trim(tile_fraction_file)//'>'
   write(     *     ,*)'tile fraction variable is <'//trim(tilefraction_variable)//'>'
   write(     *     ,*)'land_fraction_file     is <'//trim(land_fraction_file)//'>'
   write(     *     ,*)'land fraction variable is <'//trim(landfraction_variable)//'>'
   write(     *     ,*)'physical_latlon_file   is <'//trim(physical_latlon_file)//'>'
   write(     *     ,*)'longitude variable     is <'//trim(physical_lon_variable)//'>'
   write(     *     ,*)'latitude variable      is <'//trim(physical_lat_variable)//'>'

   write(logfileunit,*)
   write(logfileunit,*)'tile_fraction_file     is <'//trim(tile_fraction_file)//'>'
   write(logfileunit,*)'tile fraction variable is <'//trim(tilefraction_variable)//'>'
   write(logfileunit,*)'land_fraction_file     is <'//trim(land_fraction_file)//'>'
   write(logfileunit,*)'land fraction variable is <'//trim(landfraction_variable)//'>'
   write(logfileunit,*)'physical_latlon_file   is <'//trim(physical_latlon_file)//'>'
   write(logfileunit,*)'longitude variable     is <'//trim(physical_lon_variable)//'>'
   write(logfileunit,*)'latitude variable      is <'//trim(physical_lat_variable)//'>'
endif

return
end subroutine get_physical_filenames


!-----------------------------------------------------------------------
!>
!> Set the outermost edge of the domain. The gridcell locations as given
!> represent the gridcell centers. The edges must be computed.
!> Care is taken to consider spanning the prime meridian.
!>

subroutine set_grid_boundary()

! the physical grid is stored (Nlon,Nlat) ... with the following layout 
! physical_longitudes(  : ,1) specify the westernmost longitudes 
! physical_longitudes(  1 ,1) specifies the NorthWest corner
! physical_longitudes(Nlon,1) specifies the SouthWest corner
!
! physical_longitudes(   :,Nlat) specify the easternmost longitudes 
! physical_longitudes(   1,Nlat) specifies the NorthEast corner
! physical_longitudes(Nlon,Nlat) specifies the SouthEast corner

! type(location_type) :: ll_boundary  ! lower  left grid BOUNDARY
! type(location_type) :: ur_boundary  ! upper right grid BOUNDARY

integer  :: x_i, y_i
real(r8) :: dx, dy, lon_boundary, lat_boundary

! LOWER left ... aka ... southwest

x_i = Nlon
y_i = 1

!> @todo FIXME what if dx wraps ...

dx = physical_longitudes(x_i  , y_i+1) - physical_longitudes(x_i, y_i)
dy = physical_latitudes( x_i-1, y_i  ) - physical_latitudes( x_i, y_i)

!> @todo FIXME what if lon_boundary wraps ...

lon_boundary = physical_longitudes(x_i,y_i) - dx/2.0_r8
lat_boundary = physical_latitudes( x_i,y_i) - dy/2.0_r8

ll_boundary = set_location(lon_boundary, lat_boundary, 0.0_r8, VERTISUNDEF)

! UPPER right ... aka northeast

x_i = 1
y_i = Nlat

dx = abs(physical_longitudes(x_i, y_i) - physical_longitudes(x_i, y_i-1))
dy = abs(physical_latitudes( x_i, y_i) - physical_latitudes( x_i+1, y_i))

lon_boundary = physical_longitudes(x_i,y_i) + dx/2.0_r8
lat_boundary = physical_latitudes( x_i,y_i) + dy/2.0_r8

ur_boundary = set_location(lon_boundary, lat_boundary, 0.0_r8, VERTISUNDEF)

if (debug > 0) then
   call write_location(x_i, ll_boundary, charstring=string3)
   write(string1,*)'..  lower left  boundary ',trim(string3)
   call write_location(x_i, ur_boundary, charstring=string3)
   write(string2,*)'upper right boundary ',trim(string3)
   call error_handler(E_MSG, 'set_grid_boundary:', string1, text2=string2)
endif

!> @todo FIXME check to see what happens if upper right is > 360.0 and location is < 10.0

end subroutine set_grid_boundary

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
