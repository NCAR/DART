! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program location_test2

! Simple test program to exercise oned location module.

use location_mod
use types_mod,     only : r8
use utilities_mod, only : get_unit, error_handler, E_ERR, &
                          initialize_utilities, finalize_utilities

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type(location_type) :: loc1, loc2
integer  :: iunit, i
real(r8) :: loc2_val

call initialize_utilities('location_test2')

! Open an output file
iunit = get_unit()
open(iunit, file = 'location_test_file')

! Set the first location
loc1     = set_location(0.4_r8)
loc2_val = get_location(loc1)

if(loc2_val /= 0.4_r8) call error_handler(E_ERR,'main program unit',&
       'Error or rounding error',  source, revision, revdate)

! Write this location to the file
call write_location(iunit, loc1)

! Loop to set up four other locations
do i = 1, 4
   call write_location(iunit, set_location(i * 1.0_r8 / 4))
end do

close(iunit)

! Now read them back in and compute the distances from loc1
open(iunit, file = 'location_test_file')
loc1 = read_location(iunit)

do i = 1, 4
   loc2 = read_location(iunit)
   write(*, *) 'distance ', i, ' is ', get_dist(loc1, loc2)
end do

close(iunit)

call finalize_utilities()

end program location_test2

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
