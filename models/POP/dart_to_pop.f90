! Data Assimilation Research Testbed -- DART
! Copyright 2004-2009, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
program dart_to_pop

!----------------------------------------------------------------------
! purpose: interface between DART and the POP model
!
! method: Read DART state vector (in file 'assim_model_state_ic') and 
!         overwrite values in a POP 'restart' file.
!         Must do something with the advance-to time.
!
! author: Tim Hoar 25 Jun 09
!
!----------------------------------------------------------------------

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r4, r8, SECPERDAY
use    utilities_mod, only : E_ERR, E_WARN, E_MSG, error_handler, open_file, &
                             initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit
use        model_mod, only : static_init_model, sv_to_restart_file, &
                             get_model_size, get_model_time_step, &
                             set_model_end_time
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, get_time, print_time, print_date, &
                             operator(-)

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 128) :: dart_to_pop_input_file   = 'assim_model_state_ic'
character (len = 128) :: dart_to_pop_restart_file = 'my_pop_restart_file'

namelist /dart_to_pop_nml/ dart_to_pop_input_file, dart_to_pop_restart_file

!----------------------------------------------------------------------

integer               :: iunit, io, x_size
integer               :: secs, days
type(time_type)       :: model_time, adv_to_time
type(time_type)       :: model_timestep, offset
real(r8), allocatable :: statevector(:)

!----------------------------------------------------------------------

call initialize_utilities('dart_to_pop')

! Call model_mod:static_init_model(), which reads the namelists
! to set calendar type, grid sizes, etc.

call static_init_model()

x_size = get_model_size()
allocate(statevector(x_size))

! Read the namelist to get the input and output filenames. 

call find_namelist_in_file("input.nml", "dart_to_pop_nml", iunit)
read(iunit, nml = dart_to_pop_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_pop_nml")

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_to_pop_input_file)
call aread_state_restart(model_time, statevector, iunit, adv_to_time)
call close_restart(iunit)

!----------------------------------------------------------------------
! update the current POP state vector
!----------------------------------------------------------------------

call sv_to_restart_file(statevector, dart_to_pop_restart_file, &
                        model_time, adv_to_time)

!iunit = open_file('data.cal.DART',form='formatted',action='rewind')
!write(iunit, nml=CAL_NML)
!close(iunit)

!----------------------------------------------------------------------
! convert the adv_to_time to the appropriate number of POP
! time_manager_nml: stop_option, stop_count increments
!----------------------------------------------------------------------

model_timestep = get_model_time_step()
offset         = adv_to_time - model_time

!call set_model_end_time(offset)
!call write_data_namelistfile()

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call get_time(offset, secs, days)

call print_date( model_time,'dart_to_pop:dart model date')
call print_date(adv_to_time,'dart_to_pop:advance_to date')
call print_time( model_time,'dart_to_pop:dart model time')
call print_time(adv_to_time,'dart_to_pop:advance_to time')
call print_time(     offset,'dart_to_pop:a distance of')
write(    *      ,'(''dart_to_pop:PARM03   endTime '',i,'' seconds'')') &
                   (secs + days*SECPERDAY)

call print_date( model_time,'dart_to_pop:dart model date',logfileunit)
call print_date(adv_to_time,'dart_to_pop:advance_to date',logfileunit)
call print_time( model_time,'dart_to_pop:dart model time',logfileunit)
call print_time(adv_to_time,'dart_to_pop:advance_to time',logfileunit)
call print_time(     offset,'dart_to_pop:  a distance of',logfileunit)
write(logfileunit,'(''dart_to_pop:PARM03   endTime '',i,'' seconds'')') &
                   (secs + days*SECPERDAY)

call finalize_utilities()

end program dart_to_pop
