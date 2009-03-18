! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module location_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! Implements location interfaces for a two dimensional spherical shell. 
! The internal representation of the location is currently implemented
! as radians from 0 to 2 PI for longitude and -PI/2 to PI/2 for latitude to
! minimize computational cost for distances. However, the external 
! representation is longitude in degrees from 0 to 360 and latitude 
! from -90 to 90 for consistency with most applications in the field.

use      types_mod, only : r8, DEG2RAD, RAD2DEG, PI, MISSING_R8
use  utilities_mod, only : register_module, error_handler, E_ERR
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none
private

public :: location_type, get_dist, get_location, set_location, &
          set_location2, set_location_missing, is_location_in_region, &
          write_location, read_location, interactive_location, &
          get_close_obs, alloc_get_close_obs, &
          operator(==), operator(/=)

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type location_type
   private
   real(r8) :: lon, lat
end type location_type

type(random_seq_type) :: ran_seq
logical :: ran_seq_init = .false.
logical, save :: module_initialized = .false.

interface operator(==); module procedure loc_eq; end interface
interface operator(/=); module procedure loc_ne; end interface

contains


  subroutine initialize_module
!----------------------------------------------------------------------------
! subroutine initialize_module

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module



function get_dist(loc1, loc2)
!----------------------------------------------------------------------------

implicit none

type(location_type), intent(in) :: loc1, loc2
real(r8) :: get_dist

real(r8) :: lon_dif

if ( .not. module_initialized ) call initialize_module

! Returns distance in radians (independent of diameter of sphere)

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



function get_location(loc)
!---------------------------------------------------------------------------
!
! Given a location type, return the longitude and latitude

implicit none

type(location_type), intent(in) :: loc
real(r8), dimension(2) :: get_location

if ( .not. module_initialized ) call initialize_module

get_location(1) = loc%lon * RAD2DEG
get_location(2) = loc%lat * RAD2DEG

end function get_location



function loc_eq(loc1,loc2)
!---------------------------------------------------------------------------
!
! interface operator used to compare two locations.
! Returns true only if all components are 'the same' to within machine
! precision.

implicit none

type(location_type), intent(in) :: loc1, loc2
logical :: loc_eq

if ( .not. module_initialized ) call initialize_module

loc_eq = .false.

if ( abs(loc1%lon  - loc2%lon ) > epsilon(loc1%lon ) ) return
if ( abs(loc1%lat  - loc2%lat ) > epsilon(loc1%lat ) ) return

loc_eq = .true.

end function loc_eq



function loc_ne(loc1,loc2)
!---------------------------------------------------------------------------
!
! interface operator used to compare two locations.
! Returns true if locations are not identical to machine precision.

implicit none

type(location_type), intent(in) :: loc1, loc2
logical :: loc_ne

if ( .not. module_initialized ) call initialize_module

loc_ne = (.not. loc_eq(loc1,loc2))

end function loc_ne



function get_location_lon(loc)
!---------------------------------------------------------------------------
!
! Given a location type, return the longitude

implicit none

type(location_type), intent(in) :: loc
real(r8) :: get_location_lon

if ( .not. module_initialized ) call initialize_module

get_location_lon = loc%lon * RAD2DEG    

end function get_location_lon



function get_location_lat(loc)
!---------------------------------------------------------------------------
!
! Given a location type, return the latitude

implicit none

type(location_type), intent(in) :: loc
real(r8) :: get_location_lat

if ( .not. module_initialized ) call initialize_module

get_location_lat = loc%lat * RAD2DEG      

end function get_location_lat



function set_location(lon, lat)
!----------------------------------------------------------------------------
!
! Given a longitude and latitude
! puts this value into the location.

implicit none

type (location_type) :: set_location
real(r8), intent(in) :: lon, lat

if ( .not. module_initialized ) call initialize_module

if(lon < 0.0_r8 .or. lon > 360.0_r8) then
   write(errstring,*)'longitude (',lon,') is not within range [0,360]'
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
endif

if(lat < -90.0_r8 .or. lat > 90.0_r8) then
   write(errstring,*)'latitude (',lat,') is not within range [-90,90]'
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
endif

set_location%lon = lon * DEG2RAD
set_location%lat = lat * DEG2RAD

end function set_location


function set_location2(list)
!----------------------------------------------------------------------------
!
! location semi-independent interface routine
! given 2 float numbers, call the underlying set_location routine

implicit none

type (location_type) :: set_location2
real(r8), intent(in) :: list(:)

character(len=129) :: errstring

if ( .not. module_initialized ) call initialize_module

if (size(list) /= 2) then
   write(errstring,*)'requires 2 input values'
   call error_handler(E_ERR, 'set_location2', errstring, source, revision, revdate)
endif

set_location2 = set_location(list(1), list(2))

end function set_location2




function set_location_missing()
!----------------------------------------------------------------------------
!
! Given a longitude and latitude
! puts this value into the location.

implicit none

type (location_type) :: set_location_missing

if ( .not. module_initialized ) call initialize_module

set_location_missing%lon = MISSING_R8
set_location_missing%lat = MISSING_R8

end function set_location_missing



subroutine write_location(ifile, loc, fform)
!----------------------------------------------------------------------------
!
! Writes a 2D location to the file. Implemented as a subroutine but  could
! rewrite as a function with error control info returned. For initial implementation,
! ifile is just an integer file unit number. Probably want to replace this with ifile
! as a file_type allowing more flexibility for IO at later point. file_type and 
! associated operations would have to be supported. The mpp_io intefaces are a good
! place to head with this, perhaps, when we need to extend to supporting parallel
! platforms. 

implicit none

integer,                    intent(in) :: ifile
type(location_type),        intent(in) :: loc
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! For now, output a character tag followed by the r8 value.

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) loc%lon, loc%lat
   CASE DEFAULT
      write(ifile, '(''loc2s'')' ) 
      write(ifile, *) loc%lon, loc%lat
END SELECT

end subroutine write_location



function read_location(ifile, fform)
!----------------------------------------------------------------------------
!
! Reads a 2D location from ifile that was written by write_location. 
! See write_location for additional discussion.

implicit none

integer, intent(in) :: ifile
type(location_type) :: read_location
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_location%lon, read_location%lat
   CASE DEFAULT
      read(ifile, '(a5)' ) header
      if(header /= 'loc2s') then
         write(errstring,*)'Expected location header "loc2s" in input file, got ', header 
         call error_handler(E_ERR, 'read_location', errstring, source, revision, revdate)
      endif
! Now read the location data value
      read(ifile, *) read_location%lon, read_location%lat
END SELECT

end function read_location



subroutine interactive_location(location, set_to_default)
!--------------------------------------------------------------------------
!
! Allows for interactive input of a location. Also gives option of selecting
! a uniformly distributed random location (what the heck).

implicit none

type(location_type), intent(out) :: location
logical, intent(in), optional :: set_to_default

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

   ! Latitude must be area weightedA
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

subroutine alloc_get_close_obs(num, obs, cutoff, obs_box)

implicit none

integer, intent(in) :: num
type(location_type), intent(in) :: obs(num)
real(r8), intent(in) :: cutoff
integer, intent(out) :: obs_box(num)

! This does pre-computing for close obs; no function needed in one dimension

return

end subroutine alloc_get_close_obs


!----------------------------------------------------------------------------

subroutine get_close_obs(base_ob, num, obs, cutoff, obs_box, num_close, close_ind, dist)

! Default version with no smarts; no need to be smart in 1D

implicit none

integer, intent(in) :: base_ob, num
type(location_type), intent(in) :: obs(num)
real(r8), intent(in) :: cutoff
integer, intent(in) :: obs_box(num)
integer, intent(out) :: num_close, close_ind(num)
real(r8), intent(out) :: dist(num)

integer :: i
real(r8) :: this_dist

! Return list of obs that are within cutoff and their distances
num_close = 0
do i = 1, num
   this_dist = get_dist(obs(base_ob), obs(i))
   if(this_dist <= cutoff) then
      ! Add this ob to the list
      num_close = num_close + 1
      close_ind(num_close) = i
      dist(num_close) = this_dist 
   endif
end do

end subroutine get_close_obs



function is_location_in_region(loc, minl, maxl)
!----------------------------------------------------------------------------
!
! Returns true if the given location is between the other two.

implicit none

logical                          :: is_location_in_region
type(location_type), intent(in)  :: loc, minl, maxl


character(len=129) :: errstring

if ( .not. module_initialized ) call initialize_module


! assume failure and return as soon as we are confirmed right.
! set to success only at the bottom after all tests have passed.
is_location_in_region = .false.

if ((loc%lon < minl%lon) .or. (loc%lon > maxl%lon)) return
if ((loc%lat < minl%lat) .or. (loc%lat > maxl%lat)) return
 
is_location_in_region = .true.

end function is_location_in_region


!----------------------------------------------------------------------------
! end of location/twod_sphere/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod
