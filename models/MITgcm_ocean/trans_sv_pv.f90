! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
 
program trans_sv_pv

!----------------------------------------------------------------------
! purpose: interface between DART and the MITgcm_ocean model
!
! method: Read DART state vector (in file 'assim_model_state_ic') and 
!         write out MITgcm_ocean "snapshot" files.
!
! author: Tim Hoar 5Apr08
!
!----------------------------------------------------------------------

use        types_mod, only : r4, r8, SECPERDAY
use    utilities_mod, only : E_ERR, E_WARN, E_MSG, error_handler, open_file, &
                             initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit
use        model_mod, only : sv_to_snapshot_files, static_init_model, &
                             get_model_size, DARTtime_to_MITtime, &
                             get_model_time_step, write_data_namelistfile, &
                             set_model_end_time
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, get_time, print_time, print_date, &
                             operator(-)

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character (len = 128) :: file_in  = 'assim_model_state_ic'

!------------------------------------------------------------------
! The time manager namelist variables
! some types/etc come from   <mitsource>/pkg/cal/cal.h
! some useful insight from   cal_set.F, cal_convdate.F
!
! startDate_1 (integer) yyyymmdd   "start date of the integration"
! startDate_2 (integer) hhmmss
!------------------------------------------------------------------

character(len=9) :: TheCalendar = 'gregorian'
integer          :: startDate_1 = 19530101
integer          :: startDate_2 =          60000
logical          :: calendarDumps = .false.

NAMELIST /CAL_NML/ TheCalendar, startDate_1, startDate_2, calendarDumps

!----------------------------------------------------------------------

integer               :: iunit, io, x_size
integer               :: secs, days
type(time_type)       :: model_time, adv_to_time
type(time_type)       :: model_timestep, offset
real(r8), allocatable :: statevector(:)

!----------------------------------------------------------------------

call initialize_utilities('trans_sv_pv')
call static_init_model()

! MIT calendar information. The namelist is already read in 
! static_init_model(), so no further bulletproofing is needed here.
call find_namelist_in_file("data.cal", "CAL_NML", iunit)
read(iunit, nml = CAL_NML, iostat = io)
call check_namelist_read(iunit, io, "CAL_NML")

x_size = get_model_size()
allocate(statevector(x_size))

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(file_in)
call aread_state_restart(model_time, statevector, iunit, adv_to_time)
call close_restart(iunit)

!----------------------------------------------------------------------
! update the CAL_NML variables so we can rewrite that namelist.
! Ultimately, we want to keep data:&PARM03:startTime = 0.,
!----------------------------------------------------------------------

call DARTtime_to_MITtime(  model_time, startDate_1, startDate_2)
call sv_to_snapshot_files(statevector, startDate_1, startDate_2)

iunit = open_file('data.cal.DART',form='formatted',action='rewind')
write(iunit, nml=CAL_NML)
close(iunit)

!----------------------------------------------------------------------
! convert the adv_to_time to the appropriate number of seconds.
!----------------------------------------------------------------------

model_timestep = get_model_time_step()
offset         = adv_to_time - model_time

call set_model_end_time(offset)
call write_data_namelistfile()

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call get_time(offset, secs, days)

call print_date( model_time,'trans_sv_pv:dart model date')
call print_date(adv_to_time,'trans_sv_pv:advance_to date')
call print_time( model_time,'trans_sv_pv:dart model time')
call print_time(adv_to_time,'trans_sv_pv:advance_to time')
call print_time(     offset,'trans_sv_pv:a distance of')
write(    *      ,'(''trans_sv_pv:PARM03   endTime '',i8,'' seconds'')') &
                   (secs + days*SECPERDAY)

call print_date( model_time,'trans_sv_pv:dart model date',logfileunit)
call print_date(adv_to_time,'trans_sv_pv:advance_to date',logfileunit)
call print_time( model_time,'trans_sv_pv:dart model time',logfileunit)
call print_time(adv_to_time,'trans_sv_pv:advance_to time',logfileunit)
call print_time(     offset,'trans_sv_pv:  a distance of',logfileunit)
write(logfileunit,'(''trans_sv_pv:PARM03   endTime '',i8,'' seconds'')') &
                   (secs + days*SECPERDAY)

call finalize_utilities()

end program trans_sv_pv

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
