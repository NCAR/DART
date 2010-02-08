! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program trans_time

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between CAM and DART time and date
!
! method: Read DART 'state vector' file (proprietary format).
!         Reform time and date into form needed by CAM.
!         Write out CAM time and date to file for use by run-pc.csh
!
! author: Kevin Raeder 8/1/03
!
!----------------------------------------------------------------------

use time_manager_mod, only : time_type, read_time, write_time, &
                             get_time, set_time, operator(-), get_date, &
                             set_calendar_type, GREGORIAN, NOLEAP
use  assim_model_mod, only : static_init_assim_model, init_assim_model, open_restart_read, close_restart, &
                             get_model_size, get_model_time, read_state_restart, assim_model_type
use    utilities_mod, only : get_unit, initialize_utilities, finalize_utilities
use        types_mod, only : r8

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

integer :: ntimes = 2, n, nhtfrq, calendar_type = GREGORIAN
integer :: file_unit(2), year, month, day, hour, minute, second, &
           cam_date, cam_tod
type(time_type)       :: dart_time(2), forecast_length
character (len = 128) :: file_name = 'temp_ic', file_out = 'times'
type(assim_model_type)  :: x

call initialize_utilities('Trans_time')

call set_calendar_type(calendar_type)

! Static init assim model calls static_init_model
call static_init_assim_model()
call init_assim_model(x)

file_unit(1) = open_restart_read(file_name)
call read_state_restart(x, file_unit(1), dart_time(1))
dart_time(2) = get_model_time(x)

file_unit(2) = get_unit()
open(unit = file_unit(2), file = file_out)
! end time is first, then beginning time
!  -namelist "&camexp START_YMD=$times[3] START_TOD=$times[4] \
!                     STOP_YMD=$times[1] STOP_TOD=$times[2] NHTFRQ=$times[5] /" \

do n=1,ntimes
   call get_date(dart_time(n), year, month, day, hour, minute, second)
   PRINT*,'date = ',year, month, day, hour, minute, second
   if (calendar_type.eq.GREGORIAN) then
      cam_date = year*10000 + month*100 + day
      cam_tod  = hour*3600 + minute*60 + second
   elseif (calendar_type.eq.NOLEAP) then
      cam_date = (1899 + year)*10000 + month*100 + day
      cam_tod  = hour*3600 + minute*60 + second
   endif
   write (file_unit(2),'(2I8)') cam_date, cam_tod
enddo

call close_restart(file_unit(1))

! calculate number of hours in forecast, and pass to history tape write frequency

forecast_length = dart_time(1) - dart_time(2)

call get_time(forecast_length, second, day)
PRINT*,'forecast length = ', day, second
hour = second/3600
minute = mod(second,3600)
if (minute.ne.0) &
   print*,' not integer number of hours; nhtfrq error in trans_time'

! convert to hours, and negative to signal units are hours

! nhtfrq = -1*((((year-1)*365 + (month-1))*30 + (day-1))*24 + hour)
nhtfrq = -1*(day*24 + hour)
write (file_unit(2),'(I8)') nhtfrq

close(file_unit(2))

call finalize_utilities()

end program trans_time
