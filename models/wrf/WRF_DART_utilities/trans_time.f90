program trans_time

!----------------------------------------------------------------------
! purpose: interface between CAM and WRF time and date
!
! method: Read DART 'state vector' file (proprietary format).
!         Reform time and date into form needed by WRF.
!         Write out WRF time and date to wrfinput file.
!
! Note: This needs to be done BEFORE update_wrf_bc
!
! author: Kevin Raeder 8/1/03
! mods for WRF: Josh Hacker 4/16/04
!
!----------------------------------------------------------------------


use time_manager_mod, only : time_type, read_time, &
                             operator(-), get_date, &
                             print_time, set_date, operator (+), &
                             set_calendar_type, GREGORIAN
use assim_model_mod,  only : static_init_assim_model, binary_restart_files

use model_mod,        only : get_wrf_date, set_wrf_date, output_wrf_time


use utilities_mod,    only : get_unit

implicit none

integer               :: calendar_type = GREGORIAN
integer               :: file_unit, &
                         year, month, day, hour, minute, second

type(time_type)       :: dart_time(2), wrf_time
character (len = 128) :: file_name = 'dart_wrf_vector'
character (len = 16)  :: file_form

call set_calendar_type(calendar_type)

! Static init assim model calls static_init_model
call static_init_assim_model()

! get form of file output from assim_model_mod
if ( binary_restart_files ) then
   file_form = 'unformatted'
else
   file_form = 'formatted'
endif
file_unit = get_unit()

open(unit = file_unit, file = file_name, form=file_form)

! this should be the init time
call get_wrf_date(year, month, day, hour, minute, second)
wrf_time = set_date(year, month, day, hour, minute, second)
dart_time(1) = read_time(file_unit, file_form)
call print_time(dart_time(1))
dart_time(2) = read_time(file_unit, file_form)

! the current time is in 2, and the next time is in 1, so increment to 2
wrf_time = dart_time(2) + wrf_time

! get new date
call get_date(wrf_time, year, month, day, hour, minute, second)

! put new date into the wrf construct
call set_wrf_date(year, month, day, hour, minute, second)

! tag the wrfinput file with the new date.
call output_wrf_time()

close(file_unit)

end program trans_time

