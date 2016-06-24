! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: dart_to_cice.f90 8565 2015-09-11 17:16:08Z hkershaw $

program dart_to_cice

!----------------------------------------------------------------------
! purpose: interface between DART and the CICE model
!
! method: Read DART state vector and overwrite values in a CICE restart file.
!         If the DART state vector has an 'advance_to_time' present, a
!         file called cice_in.DART is created with a time_manager_nml namelist 
!         appropriate to advance CICE to the requested time.
!
!         The dart_to_cice_nml namelist setting for advance_time_present 
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!
! author: C Bitz June 2016
! borrowed heavily from Tim Hoar 25 Jun 09, revised 12 July 2010
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit
use  state_vector_io_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, print_time, print_date, operator(-)
use        model_mod, only : static_init_model, sv_to_restart_file, &
                             get_model_size, get_cice_restart_filename
use     dart_cice_mod, only : write_cice_namelist

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_cice/models/CICE/dart_to_cice.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 8565 $"
character(len=128), parameter :: revdate  = "$Date: 2015-09-11 10:16:08 -0700 (Fri, 11 Sep 2015) $"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 128) :: dart_to_cice_input_file   = 'dart_restart'
logical               :: advance_time_present     = .false.

namelist /dart_to_cice_nml/ dart_to_cice_input_file, &
                           advance_time_present

!----------------------------------------------------------------------

integer               :: iunit, io, x_size
type(time_type)       :: model_time, adv_to_time
real(r8), allocatable :: statevector(:)
character (len = 128) :: cice_restart_filename = 'no_cice_restart_file'

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_cice')

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the CICE namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

x_size = get_model_size()
allocate(statevector(x_size))

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "dart_to_cice_nml", iunit)
read(iunit, nml = dart_to_cice_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_cice_nml")

call get_cice_restart_filename( cice_restart_filename )

write(*,*)
write(*,'(''dart_to_cice:converting DART file '',A, &
      &'' to CICE restart file '',A)') &
     trim(dart_to_cice_input_file), trim(cice_restart_filename)

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_to_cice_input_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! update the current CICE state vector
! Convey the amount of time to integrate the model ...
! time_manager_nml: stop_option, stop_count increments
!----------------------------------------------------------------------

call sv_to_restart_file(statevector, cice_restart_filename, model_time)

if ( advance_time_present ) then
   call write_cice_namelist(model_time, adv_to_time)
endif

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_cice:CICE model date')
call print_time( model_time,'dart_to_cice:DART model time')
call print_date( model_time,'dart_to_cice:CICE model date',logfileunit)
call print_time( model_time,'dart_to_cice:DART model time',logfileunit)

if ( advance_time_present ) then
call print_time(adv_to_time,'dart_to_cice:advance_to time')
call print_date(adv_to_time,'dart_to_cice:advance_to date')
call print_time(adv_to_time,'dart_to_cice:advance_to time',logfileunit)
call print_date(adv_to_time,'dart_to_cice:advance_to date',logfileunit)
endif

call finalize_utilities('dart_to_cice')

end program dart_to_cice

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_cice/models/CICE/dart_to_cice.f90 $
! $Id: dart_to_cice.f90 8565 2015-09-11 17:16:08Z hkershaw $
! $Revision: 8565 $
! $Date: 2015-09-11 10:16:08 -0700 (Fri, 11 Sep 2015) $
