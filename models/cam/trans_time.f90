program trans_time

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
use assim_model_mod, only : static_init_assim_model, binary_restart_files


                         
use utilities_mod, only : get_unit

implicit none

integer               :: ntimes = 2, n, nhtfrq, &
                         calendar_type = GREGORIAN
integer               :: file_unit(2), seconds, days, &
                         year, month, day, hour, minute, second, &
                         cam_date, cam_tod
type(time_type)       :: dart_time(2), forecast_length
character (len = 128) :: file_name = 'assim_model_state_ic1', file_out = 'times'
character (len = 16)  :: file_form

call set_calendar_type(calendar_type)

! Static init assim model calls static_init_model
call static_init_assim_model()

! get form of file output from assim_model_mod
if (binary_restart_files == .true.) then
   file_form = 'unformatted'
else
   file_form = 'formatted'
endif
file_unit(1) = get_unit()

! write out a test file
!open(unit = file_unit(1), file = file_name, form=file_form)
!forecast_length = set_time(0, 2)
!call write_time (file_unit(1),forecast_length)
!forecast_length = set_time(43200, 1)
!call write_time (file_unit(1),forecast_length)
!close(unit = file_unit(1))

open(unit = file_unit(1), file = file_name, form=file_form)
file_unit(2) = get_unit()
open(unit = file_unit(2), file = file_out)
! end time is first, then beginning time
!  -namelist "&camexp START_YMD=$times[3] START_TOD=$times[4] \
!                     STOP_YMD=$times[1] STOP_TOD=$times[2] NHTFRQ=$times[5] /" \

do n=1,ntimes
   dart_time(n) = read_time(file_unit(1), file_form)
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

close(file_unit(1))

! calculate number of hours in forecast, and pass to history tape write frequency

forecast_length = dart_time(1) - dart_time(2)

! kdr
!   can't use this for GREG
! call get_date(forecast_length, year, month, day, hour, minute, second)

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

end program trans_time
