! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program dart_to_model

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between TIEGCM and DART
!
! method: Read DART state vector ("proprietary" format)
!         Reform state vector back into TIEGCM fields.
!         Replace those fields on the TIEGCM restart file with the new values,
!
!         Replace 'mtime' variable in the TIEGCM restart file
!         with the 'valid time' of the DART state vector.
!         Write out updated namelist variables (e.g., model_time, adv_to_time) 
!         in a file called 'namelist_update'
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : get_unit, initialize_utilities, E_ERR, &
                             error_handler, timestamp, do_output
use        model_mod, only : model_type, get_model_size, init_model_instance, &
                             vector_to_prog_var, update_TIEGCM_restart, &
                             static_init_model
use  assim_model_mod, only : assim_model_type, aread_state_restart, &
                             open_restart_read, close_restart
use time_manager_mod, only : time_type, get_time, get_date, set_calendar_type, &
                             print_time, print_date, set_date, set_time, &      
                             operator(*),  operator(+), operator(-), &
                             operator(>),  operator(<), operator(/), &
                             operator(/=), operator(<=)


implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type(model_type)       :: var
type(time_type)        :: model_time, adv_to_time, jan1, tbase, target_time
real(r8), allocatable  :: x_state(:)
integer                :: file_unit, x_size, ens_member, io
character (len = 128)  :: file_name = 'tiegcm_restart_p.nc', file_in = 'temp_ic'
character (len = 128)  :: file_namelist_out = 'namelist_update'
integer                :: model_doy, model_hour, model_minute    ! current tiegcm mtime (doy,hour,mitute) 
integer                :: adv_to_doy, adv_to_hour, adv_to_minute ! advance tiegcm mtime
integer                :: target_doy, target_hour, target_minute ! forecast time interval in mtime format 
integer                :: utsec, year, month, day, sec, model_year

!----------------------------------------------------------------------
!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_model', output_flag=.true.)

call static_init_model()        ! reads input.nml, etc., sets the table 
x_size = get_model_size()       ! now that we know how big state vector is ...
allocate(x_state(x_size))       ! allocate space for the (empty) state vector

! Open the DART model state ... 
! Read in the time to which TIEGCM must advance.  
! Read in the valid time for the model state
! Read in state vector from DART

file_unit = open_restart_read(file_in)

call aread_state_restart(model_time, x_state, file_unit, adv_to_time)
call close_restart(file_unit)

if (do_output()) &
    call print_time(model_time,'time for restart file '//trim(file_in))
if (do_output()) &
    call print_date(model_time,'date for restart file '//trim(file_in))

if (do_output()) &
    call print_time(adv_to_time,'time for restart file '//trim(file_in))
if (do_output()) &
    call print_date(adv_to_time,'date for restart file '//trim(file_in))


! Parse the vector into TIEGCM fields (prognostic variables)
call init_model_instance(var, model_time)
call vector_to_prog_var(x_state, var)
deallocate(x_state)

!----------------------------------------------------------------------
! This program writes out parameters to a file called 'namelist_update' 
! for TIEGCM namelist update used in advance_model.csh 
!----------------------------------------------------------------------

file_unit = get_unit()
open(unit = file_unit, file = file_namelist_out)

! write fields to the binary TIEGCM restart file
call update_TIEGCM_restart(file_name, var)


! Get updated TIEGCM namelist variables 
!
call get_date(model_time, model_year, month, day, model_hour, model_minute, sec)
jan1  = set_date(model_year,1,1)
tbase = model_time - jan1    ! total time since the start of the year
call get_time(tbase, utsec, model_doy)
model_doy = model_doy + 1        ! add jan1 back in

call get_date(adv_to_time, year, month, day, adv_to_hour, adv_to_minute, sec)
jan1  = set_date(year,1,1)
tbase = adv_to_time - jan1    ! total time since the start of the year
call get_time(tbase, utsec, adv_to_doy)
adv_to_doy = adv_to_doy + 1      ! add jan1 back in

! Calculate number of hours to advance tiegcm
target_time = adv_to_time - model_time
call get_time(target_time, sec, day)
target_doy    = day
target_hour   = sec/3600
target_minute = (sec - target_hour*3600)/60

!START_YEAR
write(file_unit, *, iostat = io )  model_year
if (io /= 0 )then
   call error_handler(E_ERR,'dart_to_model:','cannot write model_year to STDOUT', &
         source,revision,revdate)
endif
!START_DAY
write(file_unit, *, iostat = io )  model_doy
if (io /= 0 )then
   call error_handler(E_ERR,'dart_to_model:','cannot write model_day to STDOUT', &
         source,revision,revdate)
endif
!SOURCE_START, START, SECSTART
write(file_unit, *, iostat = io )  model_doy,',',model_hour,',',model_minute
if (io /= 0 )then
   call error_handler(E_ERR,'dart_to_model:','cannot write mtime (day/hour/minute) to STDOUT', &
         source,revision,revdate)
endif
!STOP, SECSTOP
write(file_unit, *, iostat = io )  adv_to_doy,',',adv_to_hour,',',adv_to_minute
if (io /= 0 )then
   call error_handler(E_ERR,'dart_to_model:','cannot write adv_to_time mtime (day/hour/minute) to STDOUT', &
         source,revision,revdate)
endif
!HIST, SAVE, SECHIST, SECSAVE
write(file_unit, *, iostat = io )  target_doy,',',target_hour,',',target_minute
if (io /= 0 )then
   call error_handler(E_ERR,'dart_to_model:','cannot write target_time mtime (day/hour/minute) to STDOUT', &
         source,revision,revdate)
endif
!F107
if (size(var%vars_1d) > 0) then
  write(file_unit, *, iostat = io )  var%vars_1d(1)
  if (io /= 0 )then
    call error_handler(E_ERR,'dart_to_model:','cannot write f107 (var%vars_1d) to STDOUT', &
         source,revision,revdate)
  endif
else
  write(file_unit, *, iostat = io ) 'NA'
  if (io /= 0 )then
    call error_handler(E_ERR,'dart_to_model:','cannot write f107 (var%vars_1d) to STDOUT', &
         source,revision,revdate)
  endif
endif

close(file_unit)

!----------------------------------------------------------------------
! When called with 'end', timestamp will also call finalize_utilities()
!----------------------------------------------------------------------
call timestamp(string1=source, pos='end')

end program dart_to_model
