! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program timetest

use time_manager_mod, only : set_calendar_type, GREGORIAN, time_type, &
                             set_date, leap_year

type(time_type) :: mytime
integer :: i

! start of code

call set_calendar_type(GREGORIAN)

print *, ' '
print * ,' testing gregorian calendar 1980 to 1910'
do i = 1890, 1910

 mytime = set_date(i, 1, 1, 0, 0, 0)
 if (leap_year(mytime)) then
   print *, i, 'is a leap year'
 else
   print *, i, 'is not a leap year'
 endif

enddo

print *, ' '
print * ,' testing gregorian calendar 1990 to 2010'
do i = 1990, 2010

 mytime = set_date(i, 1, 1, 0, 0, 0)
 if (leap_year(mytime)) then
   print *, i, 'is a leap year'
 else
   print *, i, 'is not a leap year'
 endif

enddo

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
