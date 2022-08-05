! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program location_test2

! Simple test program to exercise vertical related routines in location module.

use types_mod,     only : r8
use utilities_mod, only : get_unit, error_handler, E_ERR, &
                          initialize_utilities, finalize_utilities

use location_mod

implicit none


type(location_type) :: loc1, loc2
integer  :: iunit, i, nitems
real(r8) :: loc2_val
real(r8) :: vals(4) = (/ 1.0 1.0 1.0 1.0 /)
logical  :: tf
character(len=128) :: true_false, has, vertlocal

call initialize_utilities('location_test2')

nitems = LocationDims

!! Open an output file
!iunit = get_unit()
!open(iunit, file = 'location_test_file')

! Set a test location

loc1     = set_location(vals(1:nitems))

tf = has_vertical_choice()
call set_tf_string(tf, "has", "does not have", has)
print *, 'Location type '//trim(LocationName)//' '//trim(has)//' choices for vertical coordinate'


tf = vertical_localization_on()
call set_tf_string(tf, "on", "off", vertlocal)
print *, 'Vertical localization in '//trim(LocationName)//' is '//trim(vertlocal)


print *, 'Vertical localization coordinate is ', get_vertical_localization_coord()

call set_vertical(loc1, 1.0, 1)
print*, 'set_vertical() called successfully'

!call convert_vertical_obs()
!call convert_vertical_state()

!close(iunit)

call finalize_utilities()

contains

subrouine set_tf_string(tf_flag, true_choice, false_choice, returnstring)
 logical, intent(in) :: tf_flag
 character(len=*), intent(in) :: true_choice, false_choice
 character(len=*), intent(out) :: returnstring

if (tf_flag) then
   returnstring = trim(true_choice)
else
   returnstring = trim(false_choice)
endif

end subrouine set_tf_string

end program location_test2

