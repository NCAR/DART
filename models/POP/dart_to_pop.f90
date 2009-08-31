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
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, get_time, print_time, print_date, &
                             operator(-), set_time
use        model_mod, only : static_init_model, sv_to_restart_file, &
                             get_model_size, get_model_time_step, &
                             set_model_end_time
use     dart_pop_mod, only : write_pop_namelist

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 128) :: dart_to_pop_input_file   = 'filter_ics'
character (len = 128) :: dart_to_pop_restart_file = 'my_pop_restart_file'
logical               :: advance_time_present     = .TRUE.

namelist /dart_to_pop_nml/ dart_to_pop_input_file, dart_to_pop_restart_file, &
                           advance_time_present

!----------------------------------------------------------------------

integer               :: iunit, io, x_size
integer               :: secs, days
type(time_type)       :: model_time, adv_to_time
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

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! update the current POP state vector
! Convey the amount of time to integrate the model ...
! time_manager_nml: stop_option, stop_count increments
!----------------------------------------------------------------------

call sv_to_restart_file(statevector, dart_to_pop_restart_file, model_time)

if ( advance_time_present ) then
   call write_pop_namelist(model_time, adv_to_time)
endif

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_pop:dart model date')
call print_time( model_time,'dart_to_pop:dart model time')
call print_date( model_time,'dart_to_pop:dart model date',logfileunit)
call print_time( model_time,'dart_to_pop:dart model time',logfileunit)

if ( advance_time_present ) then
call print_time(adv_to_time,'dart_to_pop:advance_to time')
call print_date(adv_to_time,'dart_to_pop:advance_to date')
call print_time(adv_to_time,'dart_to_pop:advance_to time',logfileunit)
call print_date(adv_to_time,'dart_to_pop:advance_to date',logfileunit)
endif

call finalize_utilities()

end program dart_to_pop
