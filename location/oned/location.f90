module location_mod
!
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
! Implements location interfaces for a one dimensional periodic domain. Initial 
! implementation has domain 'longitude' running from 0 to 1. May want to investigate
! allowing an arbitrary real domain size at some point.

use types_mod

private

public location_type, get_dist, get_location, set_location, write_location, read_location

type location_type
   private
   real(r8) :: x
end type location_type

contains

!----------------------------------------------------------------------------

function get_dist(loc1, loc2)

implicit none

type(location_type), intent(in) :: loc1, loc2
real(r8) :: get_dist

! Reentrant domain, if distance is greater than half wraparound the other way.
get_dist = abs(loc1%x - loc2%x)
if(get_dist > 0.5_r8) get_dist = 1.0_r8 - get_dist

end function get_dist

!----------------------------------------------------------------------------

function set_location(x)

! Given a location type and a double precision value between 0 and 1
! puts this value into the location.

implicit none

type (location_type) :: set_location
real(r8), intent(in) :: x

if(x < 0.0_r8 .or. x > 1.0_r8) then
! Need to go through error handler at some point
   write(*, *) 'Error in set_loc: value of x is out of 0->1 range'
   stop
endif
set_location%x = x

end function set_location

!---------------------------------------------------------------------------

function get_location(loc)

! Given a location type, return the x coordinate

implicit none

type(location_type), intent(in) :: loc
real(r8) :: get_location

get_location = loc%x

end function get_location

!----------------------------------------------------------------------------

subroutine write_location(file, loc)

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
11 format('loc1d')
write(file, *) loc%x

! I need to learn how to do formatted IO with the types package
!21 format(r8)

end subroutine write_location

!----------------------------------------------------------------------------

function read_location(file)

! Reads a oned location from file that was written by write_location. See write_location
! for additional discussion.

implicit none

integer, intent(in) :: file
type(location_type) :: read_location

character*5 :: header

! Will want to add additional error checks on the read
read(file, 11) header
11 format(a5)
if(header /= 'loc1d') then
   write(*, *) 'Error: Expected location header "loc1d" in input file'
   stop
endif

! Now read the location data value
read(file, *) read_location%x

end function read_location

!----------------------------------------------------------------------------

end module location_mod
