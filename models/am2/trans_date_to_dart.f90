! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program trans_date_to_dart

!----------------------------------------------------------------------
! purpose: generate a Gregorian/DART date & time from standard date and time
!
! method: Read ASCII input(/output) file containing yyyy/mm/dd hh:mm:ss . 
!         Reform time and date into form needed by DART.
!         Write out AM2 time and date to i/o file for use in input.nlm . 
!
! author: Kevin Raeder 8/18/03
!
!----------------------------------------------------------------------

use utilities_mod,    only : get_unit, initialize_utilities, finalize_utilities
use time_manager_mod, only : time_type, write_time, &
                             get_time, set_time, get_date, set_date, &
                             set_calendar_type, GREGORIAN, get_calendar_type

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer               :: calendar_type = GREGORIAN
integer               :: file_unit, seconds, &
                         year, month, day, hour, minute, second, &
                         cam_date, cam_tod
type(time_type)       :: dart_time
character (len = 128) :: file_name = 'date_greg'

call initialize_utilities('Trans_date_to_dart')

call set_calendar_type(calendar_type)
! debug
seconds = get_calendar_type()
PRINT*,'calendar type = ',seconds

file_unit = get_unit()
PRINT*,'file_unit = ',file_unit

! read in date and time 
open(unit = file_unit, file = file_name, status='old',form='formatted')
read(file_unit, '(I4,5(1X,I2))') year, month, day, hour, minute, second
PRINT*,'read in date = ',year, month, day, hour, minute, second

! create and write DART date (Gregorian)
dart_time = set_date(year, month, day, hour, minute, second)
call write_time (file_unit,dart_time)

! create and write AM2 date
cam_date = (year)*10000 + month*100 + day
cam_tod  = hour*3600 + minute*60 + second
write (file_unit,'(2I8)') cam_date, cam_tod

close(unit = file_unit)

call finalize_utilities()

end program trans_date_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
