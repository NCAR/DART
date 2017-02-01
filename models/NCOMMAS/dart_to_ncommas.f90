! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_ncommas

!----------------------------------------------------------------------
! purpose: interface between DART and the ncommas model
!
! method: Read DART state vector and overwrite values in a ncommas restart file.
!         If the DART state vector has an 'advance_to_time' present, a
!         file called ncommas_in.DART is created with a time_manager_nml namelist 
!         appropriate to advance ncommas to the requested time.
!
!         The dart_to_ncommas_nml namelist setting for advance_time_present 
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!
! author: Tim Hoar 25 Jun 09, revised 12 July 2010
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file, &
                             error_handler, E_MSG
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, print_time, print_date, operator(-), get_time
use        model_mod, only : static_init_model, sv_to_restart_file, &
                             get_model_size, get_base_time, get_ncommas_restart_filename

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 128) :: dart_to_ncommas_input_file = 'dart_restart'
logical               :: advance_time_present       = .false.

namelist /dart_to_ncommas_nml/ dart_to_ncommas_input_file, &
                           advance_time_present

!----------------------------------------------------------------------

character(len=256)    :: ncommas_restart_filename
integer               :: iunit, io, x_size, diff1, diff2
type(time_type)       :: model_time, adv_to_time, base_time
real(r8), allocatable :: statevector(:)
logical               :: verbose              = .FALSE.

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_ncommas', output_flag=verbose)

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the ncommas namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

x_size = get_model_size()
allocate(statevector(x_size))

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "dart_to_ncommas_nml", iunit)
read(iunit, nml = dart_to_ncommas_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_ncommas_nml")

call get_ncommas_restart_filename( ncommas_restart_filename )

write(*,*)
write(*,'(''dart_to_ncommas:converting DART file '',A, &
      &'' to ncommas restart file '',A)') &
     trim(dart_to_ncommas_input_file), trim(ncommas_restart_filename)

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_to_ncommas_input_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! update the current ncommas state vector
! Convey the amount of time to integrate the model ...
! time_manager_nml: stop_option, stop_count increments
!----------------------------------------------------------------------

call sv_to_restart_file(statevector, ncommas_restart_filename, model_time)

if ( advance_time_present ) then
   base_time = get_base_time(ncommas_restart_filename)
   call get_time((model_time  - base_time), diff1)
   call get_time((adv_to_time - base_time), diff2)
   iunit = open_file('times', action='write')
   write(iunit, '(I8, I8)') diff1, diff2
   call close_file(iunit)
endif

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_ncommas:ncommas  model date')
call print_time( model_time,'dart_to_ncommas:DART model time')
call print_date( model_time,'dart_to_ncommas:ncommas  model date',logfileunit)
call print_time( model_time,'dart_to_ncommas:DART model time',logfileunit)

if ( advance_time_present ) then
call print_time(adv_to_time,'dart_to_ncommas:advance_to time')
call print_date(adv_to_time,'dart_to_ncommas:advance_to date')
call print_time(adv_to_time,'dart_to_ncommas:advance_to time',logfileunit)
call print_date(adv_to_time,'dart_to_ncommas:advance_to date',logfileunit)
endif

call error_handler(E_MSG,'dart_to_ncommas','Finished successfully.',source,revision,revdate)
call finalize_utilities()

end program dart_to_ncommas

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
