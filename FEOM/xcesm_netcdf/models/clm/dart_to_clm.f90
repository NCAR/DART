! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_clm

!----------------------------------------------------------------------
! purpose: interface between DART and the CLM model
!
! method: Read DART state vector and overwrite values in a CLM restart file.
!         If the DART state vector has an 'advance_to_time' present, 
!         it is read ... but nothing happens with it at this time.
!         DART is NEVER expected to advance CLM.
!
!         The dart_to_clm_nml namelist setting for advance_time_present 
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!
! author: Tim Hoar 12 July 2011
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file, &
                             error_handler, E_MSG
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, print_time, print_date, operator(-), get_time
use        model_mod, only : static_init_model, sv_to_restart_file, &
                             get_model_size, get_clm_restart_filename

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 128) :: dart_to_clm_input_file = 'dart_restart'
logical               :: advance_time_present   = .false.

namelist /dart_to_clm_nml/ dart_to_clm_input_file, &
                           advance_time_present

!----------------------------------------------------------------------

character(len=256)    :: clm_restart_filename
integer               :: iunit, io, x_size
type(time_type)       :: model_time, adv_to_time
real(r8), allocatable :: statevector(:)

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_clm')

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the clm namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

x_size = get_model_size()
allocate(statevector(x_size))

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "dart_to_clm_nml", iunit)
read(iunit, nml = dart_to_clm_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_clm_nml")

call get_clm_restart_filename( clm_restart_filename )

write(*,*)
write(*,'(''dart_to_clm:converting DART file '',A, &
      &'' to clm restart file '',A)') &
     trim(dart_to_clm_input_file), trim(clm_restart_filename)

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_to_clm_input_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! update the current clm state vector
!----------------------------------------------------------------------

call sv_to_restart_file(statevector, clm_restart_filename, model_time)

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_clm:clm  model date')
call print_time( model_time,'dart_to_clm:DART model time')
call print_date( model_time,'dart_to_clm:clm  model date',logfileunit)
call print_time( model_time,'dart_to_clm:DART model time',logfileunit)

if ( advance_time_present ) then
   call error_handler(E_MSG,'dart_to_clm','warning: DART not configured to advance CLM', &
            source, revision, revdate)
   call print_time(adv_to_time,'dart_to_clm:advance_to time')
   call print_date(adv_to_time,'dart_to_clm:advance_to date')
   call print_time(adv_to_time,'dart_to_clm:advance_to time',logfileunit)
   call print_date(adv_to_time,'dart_to_clm:advance_to date',logfileunit)
endif

call finalize_utilities('dart_to_clm')

end program dart_to_clm

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
