! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program location_test

! Simple test program to exercise annulus location module.

use     types_mod, only : r8
use utilities_mod, only : open_file, close_file, initialize_utilities, &
                          finalize_utilities
use  location_mod

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type(location_type) :: loc0(6), loc1, loc2, loc3, loc4, locA(7), locB(5)
integer             :: iunit, iunit1, iunit2, i
character(len=102)  :: testbuf

call initialize_utilities('location_test')

print *, ''
print *, 'reading in 6 locations to use in distance tests:'
do i = 1, 6
   call interactive_location(loc0(i))
enddo

print *, ''
print *, 'testing distances between pairs of points:'
do i = 1, 5
   loc1 = loc0(i)
   loc2 = loc0(i+1)
   print *, 'prev location is ', get_location(loc1)
   print *, 'this location is ', get_location(loc2)
   print *, 'distance prev to this is ', get_dist(loc1, loc2)
   print *, 'distance this to prev is ', get_dist(loc2, loc1)
end do

print *, ''
print *, 'reading in 7 locations for file read/write tests:'

! read in 7 different locs for the tests below
do i = 1, 7
   call interactive_location(locA(i))
enddo

loc1 = locA(1)
loc2 = locA(2)

print *, ''
print *, 'location 1 is ', get_location(loc1)
print *, 'location 2 is ', get_location(loc2)
print *, 'distance   is ', get_dist(loc1, loc2)

print *, ''
print *, 'testing ascii file read/write:'

! Open an output file
iunit = open_file('location_test_file.txt', form='formatted', action='write')

! Write these locations to the file
call write_location(iunit, loc1)
call write_location(iunit, loc2)

do i = 1, 5
   call write_location(iunit, locA(i))
enddo

call close_file(iunit)

loc3 = set_location_missing()
loc4 = set_location_missing()

! Now read them back in and compute the distances between them
iunit1 = open_file('location_test_file.txt', form='formatted', action='read')

loc3 = read_location(iunit1)
loc4 = read_location(iunit1)

print *, ''
print *, 'read-back location 1 is ', get_location(loc3)
print *, 'read-back location 2 is ', get_location(loc4)
print *, 'distance   is ', get_dist(loc3, loc4)

print *, ''
print *, 'testing binary file read/write:'

! Open an output file
iunit = open_file('location_test_file.bin', form='unformatted', action='write')

! Write this location to the file
call write_location(iunit, loc1, 'unformatted')
call write_location(iunit, loc2, 'unformatted')

do i = 1, 5
   call write_location(iunit, locA(i), 'unformatted')
enddo

close(iunit)

loc3 = set_location_missing()
loc4 = set_location_missing()

! Now read them back in and compute the distances between them
iunit2 = open_file('location_test_file.bin', form='unformatted', action='read')

loc3 = read_location(iunit2, 'unformatted')
loc4 = read_location(iunit2, 'unformatted')

print *, ''
print *, 'read-back location 1 is ', get_location(loc3)
print *, 'read-back location 2 is ', get_location(loc4)
print *, 'distance   is ', get_dist(loc3, loc4)


! testing write to a character buffer
print *, ''
write(*,*) 'testing write to a char buffer'

print *, ''
print *, 'raw locations: '
do i = 1, 5
   print *, get_location(locA(i))
enddo
print *, ''
print *, 'locations formatted to a char buffer:'
do i = 1, 5
   call write_location(0, locA(i), 'formatted', testbuf)
   print *, 'string length: ', len_trim(testbuf)
enddo
do i = 1, 5
   call write_location(0, locA(i), 'formatted', testbuf)
   print *, trim(testbuf)
enddo

print *, ''
print *, 'ascii readback: '
do i = 1, 5
   locB(i) = set_location_missing()
   locB(i) = read_location(iunit1, 'formatted')
   print *, get_location(locB(i))
enddo
do i = 1, 5
   call write_location(0, locB(i), 'formatted', testbuf)
   print *, trim(testbuf)
enddo

! Now read them back in and compute the distances between them
open(iunit, file = 'location_test_file.bin', form='unformatted', action='read')

print *, ''
print *, 'binary readback: '
do i = 1, 5
   locB(i) = set_location_missing()
   locB(i) = read_location(iunit2, 'unformatted')
   print *, get_location(locB(i))
enddo
do i = 1, 5
   call write_location(0, locB(i), 'formatted', testbuf)
   print *, trim(testbuf)
enddo

call close_file(iunit1)
call close_file(iunit2)

call finalize_utilities()

end program location_test

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
