! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program trans_time

! <next few lines under version control, do not edit>
! $URL$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between ROSE and DART time and date
!
! method: Read DART 'state vector' file (proprietary format).
!         Reform time and date into form needed by ROSE.
!         Write out ROSE traget time and date to file for use by run-pc.csh
!
!         Based on "trans_time" for CAM.
!
!----------------------------------------------------------------------

use time_manager_mod, only : time_type, get_time, operator(-)
use  assim_model_mod, only : static_init_assim_model, init_assim_model, open_restart_read, close_restart, &
                             get_model_time, read_state_restart, assim_model_type
use    utilities_mod, only : get_unit
use        types_mod, only : r8

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

integer :: file_unit, day, second
real(r8) :: target_hours
type(time_type)       :: dart_time(2), forecast_length
character (len = 128) :: file_name = 'temp_ic', file_out = 'times'
type(assim_model_type)  :: x

! Static init assim model calls static_init_model
call static_init_assim_model()
call init_assim_model(x)

file_unit = open_restart_read(file_name)
call read_state_restart(x, file_unit, dart_time(1)) ! target_time
dart_time(2) = get_model_time(x) ! current model_time
call close_restart(file_unit)

file_unit = get_unit()
open(unit = file_unit, file = file_out)
forecast_length = dart_time(1) - dart_time(2)
call get_time(forecast_length, second, day)
target_hours = real(day)*24._r8 + real(second)/3600._r8     
PRINT*,'trans_time: forecast length = ', day, second, &
       'trans_time: hours =', target_hours
write (file_unit,'(f20.15)') target_hours
close(file_unit)

end program trans_time
