program location_test

! Simple test program to exercise twod_sphere location module.

use location_mod
use types_mod

implicit none

type(location_type) :: loc1, loc2
integer :: unit, i
real(r8) :: loc2_val, lon, lat

! Test distribution of random locations
call interactive_location(loc2)
do i = 1, 10
   loc1 = loc2
   call interactive_location(loc2)
   write(*, *) 'location is ', get_location(loc2)
   write(*, *) 'distance to previous is ', get_dist(loc1, loc2)
   write(*, *) 'distance to previous is ', get_dist(loc2, loc1)
end do

! Open an output file
unit = 21
open(unit, file = 'location_test_file')

! Set the first location
call interactive_location(loc1)
call interactive_location(loc2)

! Write this location to the file
call write_location(unit, loc1)
call write_location(unit, loc2)

close(unit)

! Now read them back in and compute the distances from loc1
open(unit, file = 'location_test_file')
loc1 = read_location(unit)
loc2 = read_location(unit)

write(*, *) 'location 1 is ', get_location(loc1)
write(*, *) 'location 2 is ', get_location(loc2)

write(*, *) 'distance is ', get_dist(loc1, loc2)

end program location_test

