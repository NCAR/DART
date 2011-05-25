! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program trans_time

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: extract two times from the start of a DART 'state vector'
!          file (which is a proprietary format).
!          makes year, month, day, hour, minute, second available.
!          typical use is to write the information into a text file
!          which is then picked up by the advance_model.csh script.
!
! nancy collins 3 mar 2010 
!  adapted from the CAM version by Kevin Raeder 8/1/03
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             open_file, close_file, error_handler, E_ERR
use time_manager_mod, only : time_type, read_time, write_time, &
                             get_time, set_time, operator(-), get_date, &
                             set_calendar_type, GREGORIAN

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

integer :: in_unit, out_unit, year, month, day, hour, minute, second
integer :: ios_out
type(time_type)       :: dart_now_time
type(time_type)       :: dart_advance_time
character (len = 128) :: file_name = 'temp_ic'
character (len = 128) :: file_out  = 'times'
character (len = 128) :: read_format = ''
logical               :: is_advance_file = .true.

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
! choice.
in_unit = open_restart_for_time(file_name, read_format)
dart_now_time = read_time(in_unit, read_format, ios_out)
if (is_advance_file) then
   dart_advance_time = dart_now_time
   dart_now_time = read_time(in_unit, read_format, ios_out)
else
   dart_advance_time = set_time(0, 0)
endif

! ok, up to here everything is the same for any model.  now, you have
! the time(s) and you need to output them in a format that makes
! the scripting happy.   in this case, it's ascii to a file in the
! format:  yyyy mm dd hh nn ss

! out filename another candidate for a namelist
out_unit = open_file(file_out, 'formatted', 'write')


call get_date(dart_now_time, year, month, day, hour, minute, second)
write (out_unit,'(I4.4, 5(1X, I2.2))') year, month, day, hour, minute, second
! debug
write (*, '')
write (*,'(A, I4.4, 5(1X, I2.2))') 'Now time: ', year, month, day, hour, minute, second


call get_date(dart_advance_time, year, month, day, hour, minute, second)
write (out_unit,'(I4.4, 5(1X, I2.2))') year, month, day, hour, minute, second
! debug
write (*, '')
write (*,'(A, I4.4, 5(1X, I2.2))') 'Adv time: ', year, month, day, hour, minute, second


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

character(len = *), intent(in)  :: file_name
character(len = *), intent(out) :: formatting
integer                         :: open_restart_for_time

integer         :: ios_out
type(time_type) :: temp_time

! first open formatted and try to read.  if that doesn't work,
! try again unformatted.
formatting = 'formatted'
open_restart_for_time = open_file(file_name, read_format, 'read')

temp_time = read_time(open_restart_for_time, read_format, ios_out)
if(ios_out == 0) then 
   ! It appears to be formatted, proceed
   rewind open_restart_for_time
   return
endif

! Next, try to see if an unformatted read works instead
call close_file(open_restart_for_time)

formatting = 'unformatted'
open_restart_for_time = open_file(file_name, read_format, 'read')

temp_time = read_time(open_restart_for_time, read_format, ios_out)
if(ios_out == 0) then 
   ! It appears to be unformatted, proceed
   rewind open_restart_for_time
   return
endif

! error - can't get a time from the file
formatting = ''
call error_handler(E_ERR, 'open_restart_for_time', 'Does not appear to be a DART restart file', &
                   source, revision, revdate, &
                   text2 = 'Cannot read a DART formatted time from the start of the file')

end function open_restart_for_time


end program trans_time
