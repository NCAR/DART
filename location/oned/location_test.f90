program location_test

!  $Source$
!  $Revision$
!  $Date$

! Simple test program to exercise oned location module.

use location_mod
use types_mod

implicit none

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type(location_type) :: loc1, loc2
integer :: unit, i
real(r8) :: loc2_val

! Open an output file
unit = 21
open(unit, file = 'location_test_file')

! Set the first location
loc1 = set_location(1.4_r8)
loc2_val = get_location(loc1)
if(loc2_val /= 0.4_r8) then
   write(*, *) 'Error or rounding error'
   stop
endif


! Write this location to the file
call write_location(unit, loc1)

! Loop to set up four other locations
do i = 1, 4
   call write_location(unit, set_location(i * 1.0_r8 / 4))
end do

close(unit)

! Now read them back in and compute the distances from loc1
open(unit, file = 'location_test_file')
loc1 = read_location(unit)

do i = 1, 4
   loc2 = read_location(unit)
   write(*, *) 'distance ', i, ' is ', get_dist(loc1, loc2)
end do

end program location_test

