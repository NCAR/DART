! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module location_mod

! Implements location interfaces for a two dimensional spherical shell. 
! The internal representation of the location is currently implemented
! as radians from 0 to 2 PI for longitude and -PI/2 to PI/2 for latitude to
! minimize computational cost for distances. However, the external 
! representation is longitude in degrees from 0 to 360 and latitude 
! from -90 to 90 for consistency with most applications in the field.

use      types_mod, only : r8, DEG2RAD, RAD2DEG, PI, MISSING_R8, i8
use  utilities_mod, only : error_handler, E_ERR, &
                           ascii_file_format, is_longitude_between
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform
use ensemble_manager_mod, only : ensemble_type
use default_location_mod, only : has_vertical_choice, vertical_localization_on, &
                                 get_vertical_localization_coord, &
                                 set_vertical_localization_coord, &
                                 is_vertical, set_vertical, &
                                 convert_vertical_obs, convert_vertical_state

implicit none
private

public :: location_type, get_location, set_location, &
          set_location_missing, is_location_in_region, get_maxdist, &
          write_location, read_location, interactive_location, query_location, &
          LocationDims, LocationName, LocationLName, LocationStorageOrder, LocationUnits, &
          get_close_type, get_close_init, get_close_obs, get_close_state, get_close_destroy, &
          operator(==), operator(/=), get_dist, has_vertical_choice, vertical_localization_on, &
          set_vertical, is_vertical, get_vertical_localization_coord, &
          set_vertical_localization_coord, convert_vertical_obs, convert_vertical_state

character(len=*), parameter :: source = 'twod_sphere/location_mod.f90'

type location_type
   private
   real(r8) :: lon, lat
end type location_type

! Needed as stub but not used in this low-order model
type get_close_type
   private
   integer  :: num
   real(r8) :: maxdist
end type get_close_type

type(random_seq_type) :: ran_seq
logical :: ran_seq_init = .false.
logical, save :: module_initialized = .false.

integer,          parameter :: LocationDims = 2
character(len=*), parameter :: LocationName = "loc2Dsphere"
character(len=*), parameter :: LocationLName = &
                                 "twod sphere locations: lon, lat"
character(len=*), parameter :: LocationStorageOrder = "X Y"
character(len=*), parameter :: LocationUnits = "degrees degrees"

character(len = 512) :: errstring

interface operator(==); module procedure loc_eq; end interface
interface operator(/=); module procedure loc_ne; end interface

interface set_location
   module procedure set_location_single
   module procedure set_location_array
end interface set_location

contains

!----------------------------------------------------------------------------

subroutine initialize_module
 
if (module_initialized) return

module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------------

function get_dist(loc1, loc2, type1, kind2)

! Returns distance in radians (independent of diameter of sphere)

type(location_type), intent(in) :: loc1, loc2
integer, optional,   intent(in) :: type1, kind2
real(r8)                        :: get_dist

real(r8) :: lon_dif

if ( .not. module_initialized ) call initialize_module


! Compute great circle path shortest route between two points
lon_dif = abs(loc1%lon - loc2%lon)
if(lon_dif > PI) lon_dif = 2.0_r8 * PI - lon_dif

if(cos(loc1%lat) == 0) then
   get_dist = abs(loc2%lat - loc1%lat)
else
   get_dist = acos(sin(loc2%lat) * sin(loc1%lat) + &
      cos(loc2%lat) * cos(loc1%lat) * cos(lon_dif))
endif

end function get_dist

!---------------------------------------------------------------------------

function loc_eq(loc1,loc2)
 
! Interface operator used to compare two locations.
! Returns true only if all components are 'the same' to within machine
! precision.

type(location_type), intent(in) :: loc1, loc2
logical                         :: loc_eq

if ( .not. module_initialized ) call initialize_module

loc_eq = .false.

if ( abs(loc1%lon  - loc2%lon ) > epsilon(loc1%lon ) ) return
if ( abs(loc1%lat  - loc2%lat ) > epsilon(loc1%lat ) ) return

loc_eq = .true.

end function loc_eq

!---------------------------------------------------------------------------

function loc_ne(loc1,loc2)
 
! Interface operator used to compare two locations.
! Returns true if locations are not identical to machine precision.

type(location_type), intent(in) :: loc1, loc2
logical                         :: loc_ne

if ( .not. module_initialized ) call initialize_module

loc_ne = (.not. loc_eq(loc1,loc2))

end function loc_ne

!---------------------------------------------------------------------------

function get_location(loc)
 
! Given a location type, return the longitude and latitude

type(location_type), intent(in) :: loc
real(r8), dimension(2) :: get_location

if ( .not. module_initialized ) call initialize_module

get_location(1) = loc%lon * RAD2DEG
get_location(2) = loc%lat * RAD2DEG

end function get_location

!----------------------------------------------------------------------------

function set_location_single(lon, lat)
 
! Given a longitude and latitude, put this value into the location.

real(r8), intent(in) :: lon, lat
type (location_type) :: set_location_single

if ( .not. module_initialized ) call initialize_module

if(lon < 0.0_r8 .or. lon > 360.0_r8) then
   write(errstring,*)'longitude (',lon,') is not within range [0,360]'
   call error_handler(E_ERR, 'set_location', errstring, source)
endif

if(lat < -90.0_r8 .or. lat > 90.0_r8) then
   write(errstring,*)'latitude (',lat,') is not within range [-90,90]'
   call error_handler(E_ERR, 'set_location', errstring, source)
endif

set_location_single%lon = lon * DEG2RAD
set_location_single%lat = lat * DEG2RAD

end function set_location_single

!----------------------------------------------------------------------------

function set_location_array(list)
 
! location semi-independent interface routine
! given 2 float numbers, call the underlying set_location routine

real(r8), intent(in) :: list(:)
type (location_type) :: set_location_array

if ( .not. module_initialized ) call initialize_module

if (size(list) < 2) then
   write(errstring,*)'requires 2 input values'
   call error_handler(E_ERR, 'set_location', errstring, source)
endif

set_location_array = set_location_single(list(1), list(2))

end function set_location_array

!----------------------------------------------------------------------------

function set_location_missing()

! Initialize a location type to indicate the contents are unset.

type (location_type) :: set_location_missing

if ( .not. module_initialized ) call initialize_module

set_location_missing%lon = MISSING_R8
set_location_missing%lat = MISSING_R8

end function set_location_missing

!---------------------------------------------------------------------------

function query_location(loc, attr)
 
! Returns the value of the attribute

type(location_type),        intent(in) :: loc
character(len=*), optional, intent(in) :: attr
real(r8)                               :: query_location

if ( .not. module_initialized ) call initialize_module

! see the long comment in this routine in the threed_sphere
! module for warnings about compiler bugs before you change
! this code.


! the default value
query_location = loc%lon
if (.not. present(attr)) return

select case(attr)
 case ('lon','LON')
   query_location = loc%lon
 case ('lat','LAT')
   query_location = loc%lat
 case default
   call error_handler(E_ERR, 'query_location; twod_sphere', &
         'Only lon or lat are legal attributes to request from location', source)
end select

end function query_location

!----------------------------------------------------------------------------

subroutine write_location(locfile, loc, fform, charstring)
 
! Writes a location to a file.
! additional functionality: if optional argument charstring is specified,
! it must be long enough to hold the string, and the location information is
! written into it instead of to a file.  fform must be ascii (which is the
! default if not specified) to use this option.

integer, intent(in)                        :: locfile
type(location_type), intent(in)            :: loc
character(len = *),  intent(in),  optional :: fform
character(len = *),  intent(out), optional :: charstring

integer             :: charlength
logical             :: writebuf

! 10 format(1x,2(f22.14,1x))  ! old
10 format(1X,2(F22.16,1X)) 

if ( .not. module_initialized ) call initialize_module

! writing to a file (normal use) or to a character buffer?
writebuf = present(charstring)

! output file; test for ascii or binary, write what's asked, and return
if (.not. writebuf) then
   if (ascii_file_format(fform)) then
      write(locfile, '(''loc2s'')' )
      write(locfile, 10) loc%lon, loc%lat
   else
      write(locfile) loc%lon, loc%lat
   endif
   return
endif

! you only get here if you're writing to a buffer and not
! to a file, and you can't have binary format set.
if (.not. ascii_file_format(fform)) then
   call error_handler(E_ERR, 'write_location', &
      'Cannot use string buffer with binary format', source)
endif

! format the location to be more human-friendly; meaning
! degrees instead of radians.

! this must be the sum of the formats below.
charlength = 39

if (len(charstring) < charlength) then
   write(errstring, *) 'charstring buffer must be at least ', charlength, ' chars long'
   call error_handler(E_ERR, 'write_location', errstring, source)
endif

write(charstring, '(A,F12.8,1X,F12.8)') 'Lon/Lat(deg): ',  loc%lon*RAD2DEG, loc%lat*RAD2DEG

end subroutine write_location

!----------------------------------------------------------------------------

function read_location(locfile, fform)
 
! Reads a location from a file that was written by write_location. 
! See write_location for additional discussion.

integer, intent(in)                      :: locfile
character(len = *), intent(in), optional :: fform
type(location_type)                      :: read_location

character(len=5) :: header

if ( .not. module_initialized ) call initialize_module

if (ascii_file_format(fform)) then
   read(locfile, '(a5)' ) header
   if(header /= 'loc2s') then
      write(errstring,*)'Expected location header "loc2s" in input file, got ', header 
      call error_handler(E_ERR, 'read_location', errstring, source)
   endif
   ! Now read the location data value
   read(locfile, *) read_location%lon, read_location%lat
else
   read(locfile) read_location%lon, read_location%lat
endif

end function read_location

!--------------------------------------------------------------------------

subroutine interactive_location(location, set_to_default)
 
! Allows for interactive input of a location. Also gives option of selecting
! a uniformly distributed random location.

type(location_type), intent(out) :: location
logical, intent(in), optional    :: set_to_default

real(r8) :: lon, lat

if ( .not. module_initialized ) call initialize_module

! If set_to_default is true, then just zero out and return
if(present(set_to_default)) then
   if(set_to_default) then
      location%lon = 0.0
      location%lat = 0.0
      return
   endif
endif

write(*, *) 'Input longitude for this obs: value 0 to 360.0 or a negative number for '
write(*, *) 'Uniformly distributed random location'
read(*, *) lon

do while(lon > 360.0_r8)
   write(*, *) 'Input value greater than 360.0 is illegal, please try again'
   read(*, *) lon
end do

if(lon < 0.0_r8) then

   ! Need to make sure random sequence is initialized

   if(.not. ran_seq_init) then
      call init_random_seq(ran_seq)
      ran_seq_init = .TRUE.
   endif

   ! Longitude is random from 0 to 2 PI
   location%lon = random_uniform(ran_seq) * 2.0_r8 * PI

   ! Latitude must be area weighted
   location%lat = asin(random_uniform(ran_seq) * 2.0_r8 - 1.0_r8)

   write(*, *) 'random location is ', location%lon, location%lat

else
   write(*, *) 'Input latitude for this obs: value -90.0 to 90.0'
   read(*, *) lat
   do while(lat < -90.0_r8 .or. lat > 90.0_r8)
      write(*, *) 'Input value < -90.0 or > 90.0 is illegal, please try again'
      read(*, *) lat
   end do
   location = set_location(lon, lat)
end if

end subroutine interactive_location

!----------------------------------------------------------------------------
! Initializes get_close accelerator - unused in this location module

subroutine get_close_init(gc, num, maxdist, locs, maxdist_list)

type(get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
real(r8),             intent(in)    :: maxdist
type(location_type),  intent(in)    :: locs(:)
real(r8), intent(in), optional      :: maxdist_list(:)

! Set the maximum localization distance
gc%maxdist = maxdist

! Save the value of num_locs
gc%num = num

if (present(maxdist_list)) then
   write(errstring,*)'twod_sphere locations does not support different cutoff distances by type'
   call error_handler(E_ERR, 'get_close_init', errstring, source)
endif

end subroutine get_close_init

!----------------------------------------------------------------------------

subroutine get_close_destroy(gc)

type(get_close_type), intent(inout) :: gc

end subroutine get_close_destroy

!----------------------------------------------------------------------------

subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ensemble_handle)

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(in)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:), loc_types(:) 
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ensemble_handle

call get_close(gc, base_loc, base_type, locs, loc_qtys, &
               num_close, close_ind, dist, ensemble_handle)

end subroutine get_close_obs

!----------------------------------------------------------------------------

subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ensemble_handle)

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(in)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:)
integer(i8),                   intent(in)  :: loc_indx(:) 
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ensemble_handle

call get_close(gc, base_loc, base_type, locs, loc_qtys, &
               num_close, close_ind, dist, ensemble_handle)

end subroutine get_close_state

!----------------------------------------------------------------------------

subroutine get_close(gc, base_loc, base_type, locs, loc_qtys, &
                     num_close, close_ind, dist, ensemble_handle)

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(in)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ensemble_handle

integer :: i
real(r8) :: this_dist

! the list of locations in the locs() argument must be the same
! as the list of locations passed into get_close_init(), so
! gc%num and size(locs) better be the same.   if the list changes
! you have to destroy the old gc and init a new one.
if (size(locs) /= gc%num) then
   write(errstring,*)'locs() array must match one passed to get_close_init()'
   call error_handler(E_ERR, 'get_close', errstring, source)
endif

! Return list of obs that are within maxdist and their distances
num_close = 0
do i = 1, gc%num
   this_dist = get_dist(base_loc, locs(i), base_type, loc_qtys(i))
   if(this_dist <= gc%maxdist) then
      ! Add this ob to the list
      num_close = num_close + 1
      close_ind(num_close) = i
      if (present(dist)) dist(num_close) = this_dist 
   endif
end do

end subroutine get_close

!---------------------------------------------------------------------------

function get_maxdist(gc, obs_type)
type(get_close_type), intent(in) :: gc
integer, optional,    intent(in) :: obs_type
real(r8) :: get_maxdist

get_maxdist = gc%maxdist

end function get_maxdist

!----------------------------------------------------------------------------

function is_location_in_region(loc, minl, maxl)
 
! Returns true if the given location is between the other two.
! will handle wrap in longitude, but not in latitude over poles.

logical                          :: is_location_in_region
type(location_type), intent(in)  :: loc, minl, maxl

if ( .not. module_initialized ) call initialize_module

! assume failure and return as soon as we are confirmed right.
! set to success only at the bottom after all tests have passed.
is_location_in_region = .false.

! latitude: we do not allow wrap of rectangular regions over the poles.
if ((loc%lat < minl%lat) .or. (loc%lat > maxl%lat)) return

! use common routine in utilities module to do all the wrapping
if (.not. is_longitude_between(loc%lon, minl%lon, maxl%lon, doradians=.TRUE.)) return
 
is_location_in_region = .true.

end function is_location_in_region


!----------------------------------------------------------------------------
! end of location/twod_sphere/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod

