! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!-----------------------------------------------------------------------
!> One of the interface executables between DART and ROMS.
!>
!> Reads a file containing a DART state vector and overwrite the 
!> corresponding values in a ROMS model analysis file.
!> If the DART state vector has an 'advance_to_time' present, a
!> file called model_in.DART is created with a time_manager_nml namelist 
!> appropriate to advance model to the requested time.
!>
!> The dart_to_roms_nml namelist setting for advance_time_present 
!> determines whether or not the input file has an 'advance_to_time'.
!> Typically, only temporary files like 'assim_model_state_ic' have
!> an 'advance_to_time'.
!>
!> author: PENG XIU 12/2013 @ University of Maine
!>         peng.xiu@maine.edu
!>
!> subsequently modified by TJH 1/2015

program dart_to_roms

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file, &
                             error_handler, E_MSG, E_ERR
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, print_time, print_date, operator(-), &
                             get_time, get_date
use        model_mod, only : static_init_model, sv_to_restart_file, &
                             get_model_size, get_model_restart_filename, &
                             write_model_time, print_variable_ranges

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=256) :: dart_to_roms_input_file = 'dart.restart'
logical            :: print_data_ranges       = .true.
logical            :: advance_time_present    = .false.
character(len=256) :: time_filename           = 'roms_time'

namelist /dart_to_roms_nml/ dart_to_roms_input_file, &
                            advance_time_present,    &
                            print_data_ranges,       &
                            time_filename

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: io, iunit, x_size
type(time_type)       :: model_time, adv_to_time
real(r8), allocatable :: statevector(:)
character(len=256)    :: roms_restart_filename
character(len=512)    :: string1, string2

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_roms')

!----------------------------------------------------------------------
! Read the namelist to get the input filename
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "dart_to_roms_nml", iunit)
read(iunit, nml = dart_to_roms_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_roms_nml")

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the model namelists
! to set grid sizes, etc.  The 'ocean_time' times in the ROMS file
! cannot (at present) be correctly decoded to give dates for which we
! have real observations. Either the test files I have are from AD 223,
! or I am doing it wrong. Since I cannot decode them sensibly now,
! I am setting the time from the namelist value.
!----------------------------------------------------------------------

call static_init_model()

x_size = get_model_size()
allocate(statevector(x_size))

call get_model_restart_filename( roms_restart_filename )

write(string1,*) '..  converting  DART file <'//trim(dart_to_roms_input_file)//'>'
write(string2,*) 'to ROMS analysis file <'//trim(roms_restart_filename)//'>' 
call error_handler(E_MSG,'dart_to_roms:',string1,text2=string2)

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_to_roms_input_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

if (print_data_ranges) call print_variable_ranges(statevector)

!----------------------------------------------------------------------
! update the current model state vector
! Convey the amount of time to integrate the model ...
! time_manager_nml: stop_option, stop_count increments
!----------------------------------------------------------------------

call sv_to_restart_file(statevector, roms_restart_filename, model_time)

! write time into in text format (YYYY-MM-DD_hh:mm:ss) into a file.
! if advance time is there, write the current time then advance time.
! otherwise just write current time.
if ( advance_time_present ) then
   call write_model_time(time_filename, model_time, adv_to_time)
else
   call write_model_time(time_filename, model_time)
endif

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_roms:ROMS date')
call print_time( model_time,'dart_to_roms:ROMS time')
call print_date( model_time,'dart_to_roms:ROMS date',logfileunit)
call print_time( model_time,'dart_to_roms:ROMS time',logfileunit)

if ( advance_time_present ) then
call print_time(adv_to_time,'dart_to_roms:advance_to time')
call print_date(adv_to_time,'dart_to_roms:advance_to date')
call print_time(adv_to_time,'dart_to_roms:advance_to time',logfileunit)
call print_date(adv_to_time,'dart_to_roms:advance_to date',logfileunit)
endif

call finalize_utilities()

end program dart_to_roms

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

