module location_mod
!
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
! Implements location interfaces for a two dimensional spherical shell. 
! The internal representation of the location is currently implemented
! as radians from 0 to 2 PI for longitude and -PI/2 to PI/2 for latitude to
! minimize computational cost for distances. However, the external 
! representation is longitude in degrees from 0 to 360 and latitude 
! from -90 to 90 for consistency with most applications in the field.

use      types_mod
use  utilities_mod, only : output_err, E_ERR
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform


private

public location_type, get_dist, get_location, set_location, &
       write_location, read_location, interactive_location

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type location_type
   private
   real(r8) :: lon, lat
end type location_type

type(random_seq_type) :: ran_seq
logical :: ran_seq_init = .false.


! CVS Generated file description for error handling, do not edit
character(len = 129), parameter :: &
   e_src = "$Source$", &
   e_rev = "$Revision$", &
   e_dat = "$Date$", &
   e_aut = "$Author$"


contains


function get_dist(loc1, loc2)
!----------------------------------------------------------------------------

implicit none

type(location_type), intent(in) :: loc1, loc2
real(r8) :: get_dist

real(r8) :: lon_dif

! Returns distance in radians (independent of diameter of sphere)

! Compute great circle path shortest route between two points
lon_dif = abs(loc1%lon - loc2%lon)
if(lon_dif > pi) lon_dif = 2.0_r8 * pi - lon_dif

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

get_location(1) = loc%lon * 360.0_r8 / (2.0_r8 * pi)
get_location(2) = loc%lat * 180.0_r8 / pi

end function get_location



function get_location_lon(loc)
!---------------------------------------------------------------------------
!
! Given a location type, return the longitude

implicit none

type(location_type), intent(in) :: loc
real(r8) :: get_location_lon

get_location_lon = loc%lon * 360.0_r8 / (2.0_r8 * pi)

end function get_location_lon



function get_location_lat(loc)
!---------------------------------------------------------------------------
!
! Given a location type, return the latitude

implicit none

type(location_type), intent(in) :: loc
real(r8) :: get_location_lat

get_location_lat = loc%lat * 180.0_r8 / pi

end function get_location_lat



function set_location(lon, lat)
!----------------------------------------------------------------------------
!
! Given a longitude and latitude
! puts this value into the location.

implicit none

type (location_type) :: set_location
real(r8), intent(in) :: lon, lat

if(lon < 0.0_r8 .or. lon > 360.0_r8) call output_err(E_ERR, e_src, e_rev, e_dat, e_aut, &
   'set_location', 'Longitude is out of 0->360 range')

if(lat < -90.0_r8 .or. lat > 90.0_r8) call output_err(E_ERR, e_src, e_rev, e_dat, e_aut, &
   'set_location', 'Latitude is out of -90->90 range')

set_location%lon = lon * 2.0_r8 * pi / 360.0_r8
set_location%lat = lat * pi / 180.0_r8

end function set_location



subroutine write_location(file, loc)
!----------------------------------------------------------------------------
!
! Writes a oned location to the file. Implemented as a subroutine but  could
! rewrite as a function with error control info returned. For initial implementation,
! file is just an integer file unit number. Probably want to replace this with file
! as a file_type allowing more flexibility for IO at later point. file_type and 
! associated operations would have to be supported. The mpp_io intefaces are a good
! place to head with this, perhaps, when we need to extend to supporting parallel
! platforms. 

implicit none

integer, intent(in) :: file
type(location_type), intent(in) :: loc

! For now, output a character tag followed by the r8 value. Is this written at 
! machine precision ???

write(file, 11) 
11 format('loc2s')
write(file, *) loc%lon, loc%lat

! I need to learn how to do formatted IO with the types package
!21 format(r8)

end subroutine write_location



function read_location(file)
!----------------------------------------------------------------------------
!
! Reads a oned location from file that was written by write_location. See write_location
! for additional discussion.

implicit none

integer, intent(in) :: file
type(location_type) :: read_location

character*5 :: header

! Will want to add additional error checks on the read
read(file, 11) header
11 format(a5)
if(header /= 'loc2s') call output_err(E_ERR, e_src, e_rev, e_dat, e_aut, &
   'read_location', 'Expected location header "loc1d" in input file')

! Now read the location data value
read(file, *) read_location%lon, read_location%lat

end function read_location



subroutine interactive_location(location)
!--------------------------------------------------------------------------
!
! Allows for interactive input of a location. Also gives option of selecting
! a uniformly distributed random location (what the heck).

implicit none

type(location_type), intent(out) :: location

real(r8) :: lon, lat

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

!
!----------------------------------------------------------------------------
! end of location/twod_sphere/location_mod.f90
!----------------------------------------------------------------------------
!
end module location_mod
