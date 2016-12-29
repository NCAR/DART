! This code may (or may not) be part of the NOGAPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

program trans_time

!----------------------------------------------------------------------
! purpose: interface between NOGAPS and DART time and date
!
! method: Read times out of the start of a DART 'state vector' file (proprietary format).
!         Reform time and date into form needed by NOGAPS scripts.
!         Write out time and date to file.
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             open_file, close_file, error_handler, E_ERR
use time_manager_mod, only : time_type, read_time, write_time, &
                             get_time, set_time, operator(-), get_date, &
                             set_calendar_type, GREGORIAN

implicit none

! This is from the original assembla server we used during collaboration.
! character(len=128), parameter :: &
!    source   = "$orgURL: https://svn2.assembla.com/svn/ngdart/trans_time.f90 $", &
!    revision = "$orgRevision: 108 $", &
!    revdate  = "$orgDate: 2010-06-09 14:49:52 -0600 (Wed, 09 Jun 2010) $"

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer :: ntimes = 2, n, nhtfrq
integer :: in_unit, out_unit, year, month, day, hour, minute, second, &
           cam_date, cam_tod, ios_out

type(time_type) :: dart_now_time, dart_advance_time, forecast_length

character(len=128) :: file_name = 'temp_ic'
character(len=128) :: file_out = 'time_info' 
character(len=128) :: formatting = ''
character(len=128) :: read_format, msgstring

logical :: is_advance_file = .true.

!----------------------------------------------------------------------
! read in just the first 1 or 2 times from a dart restart file.
! in theory this should work for any dart restart file because it
! doesn't read in the state vector.  however, the output format
! of the time(s) is going to be specific for any particular model.
! look below for the comment which says it's starting in on the
! model specific output code.

call initialize_utilities('trans_time')

! choose your calendar type from the time manager options.   
! (obvious candidate for a namelist item)
call set_calendar_type(GREGORIAN)

! this is supposed to autodetect whether the restart file is in
! binary or ascii, by trying to read in the initial time and
! trying the other option if it fails.  this has been known to
! have problems on some versions of the absoft compiler.
! could always add a namelist override if necessary.
! whether to read a single time or 2 times should be a namelist
! choice - which is not implemented ...
in_unit = open_restart_for_time(file_name, formatting)
dart_now_time = read_time(in_unit, formatting, ios_out)
if (is_advance_file) then
   dart_advance_time = dart_now_time
   dart_now_time = read_time(in_unit, formatting, ios_out)
else
   dart_advance_time = set_time(0, 0)
endif

! ok, up to here everything is the same for any model.  now, you have
! the time(s) and you need to output them in a format that makes
! the scripting happy.   in this case, it's ascii to a file in the
! very specific format YYYYMMDDHH

! out filename another candidate for a namelist
out_unit = open_file(file_out, 'formatted', 'write')

call get_date(dart_now_time, year, month, day, hour, minute, second)
write (out_unit,'(I4.4,3(I2.2))') year, month, day, hour
! debug
write(*,'(''trans_time:now date = '',6(1x,i4))') year, month, day, hour, minute, second

call get_date(dart_advance_time, year, month, day, hour, minute, second)
write (out_unit,'(I4.4,3(I2.2))') year, month, day, hour
! debug
write(*,'(''trans_time:adv date = '',6(1x,i4))') year, month, day, hour, minute, second

! calculate number of hours in forecast

forecast_length = dart_advance_time - dart_now_time

call get_time(forecast_length, second, day)

hour = day*24 + second/3600
minute = mod(second,3600)
if (minute .ne. 0) then
   write(msgstring,*) 'forecast length ',hour,':',minute,'not integer number of hours.'
   call error_handler(E_ERR,'trans_time',msgstring,source,revision,revdate)
endif

! output hours
write(out_unit, '(I3.3)') hour
! debug
write(*,'(''trans_time:forecast length (hours) = '',i8)') hour

call close_file(in_unit)
call close_file(out_unit)
call finalize_utilities()

contains


function open_restart_for_time(file_name, formatting)
!----------------------------------------------------------------------
!
! Opens a restart file just to read in time, cannot read in data
! because it's avoiding calling the static model init code, which
! calls the model_mod init code.  this code is a condensed version
! of the assim_model/open_restart_read routine.  return which format
! worked as part of output.   function return value is unit number.

integer :: open_restart_for_time
character(len = *), intent(in) :: file_name
character(len = *), intent(out) :: formatting

integer :: ios, ios_out
type(time_type) :: temp_time

! first open formatted and try to read.  if that doesn't work,
! try again unformatted.
read_format = 'formatted'
open_restart_for_time = open_file(file_name, read_format, 'read')

temp_time = read_time(open_restart_for_time, read_format, ios_out)
if(ios_out == 0) then 
   ! It appears to be formatted, proceed
   rewind open_restart_for_time
   formatting = read_format
   return
endif

! Next, try to see if an unformatted read works instead
call close_file(open_restart_for_time)

read_format = 'unformatted'
open_restart_for_time = open_file(file_name, read_format, 'read')

temp_time = read_time(open_restart_for_time, read_format, ios_out)
if(ios_out == 0) then 
   ! It appears to be unformatted, proceed
   rewind open_restart_for_time
   formatting = read_format
   return
endif


end function open_restart_for_time


end program trans_time

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

