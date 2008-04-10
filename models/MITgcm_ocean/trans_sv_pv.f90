! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
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

! <next few lines under version control, do not edit>
! $URL: http://subversion.ucar.edu/DAReS/DART/trunk/models/MITgcm_ocean/trans_sv_pv.f90 $
! $Id: trans_sv_pv.f90 3258 2008-03-14 15:58:36Z thoar $
! $Revision: 3258 $
! $Date: 2008-03-14 09:58:36 -0600 (Fri, 14 Mar 2008) $

use        types_mod, only : r4, r8
use    utilities_mod, only : E_ERR, E_WARN, E_MSG, error_handler, open_file, &
                             initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read
use        model_mod, only : sv_to_snapshot_files, static_init_model, &
                             get_model_size, DARTtime_to_MITtime, &
                             get_model_time_step
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, get_time, print_time, print_date, &
                             operator(-)

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL: http://subversion.ucar.edu/DAReS/DART/trunk/models/MITgcm_ocean/trans_sv_pv.f90 $", &
   revision = "$Revision: 3258 $", &
   revdate  = "$Date: 2008-03-14 09:58:36 -0600 (Fri, 14 Mar 2008) $"

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

! call print_date(model_time,'dart model date')
! call print_time(model_time,'dart model time')
! call print_date(adv_to_time,'advance_to date')
! call print_time(adv_to_time,'advance_to time')

!----------------------------------------------------------------------
! update the CAL_NML variables so we can rewrite that namelist.
! Ultimately, we want to keep data:&PARM03:startTime = 0.,
!----------------------------------------------------------------------

call DARTtime_to_MITtime(  model_time, startDate_1, startDate_2)
call sv_to_snapshot_files(statevector, startDate_1, startDate_2)

iunit = open_file('data.cal.new',form='formatted',action='rewind')
write(iunit, nml=CAL_NML)
close(iunit)

!----------------------------------------------------------------------
! convert the adv_to_time to the appropriate number of seconds.
!----------------------------------------------------------------------

model_timestep = get_model_time_step()
offset         = adv_to_time - model_time

call get_time(offset, secs, days)

write(*, nml=CAL_NML)
write(*, *)'PARM03 startTime ',0.0
write(*, *)'PARM03   endTime ',(secs + days*86400)

call finalize_utilities()

end program trans_sv_pv
