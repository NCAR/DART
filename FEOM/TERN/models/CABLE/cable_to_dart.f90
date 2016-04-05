! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program cable_to_dart

!----------------------------------------------------------------------
! purpose: interface between CABLE and DART
!
! method: Read CABLE "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The CABLE filename is read from the cable_in namelist
!         <edit cable_to_dart_output_file in input.nml:cable_to_dart_nml>
!         cable_to_dart
!
! author: Tim Hoar 20 February 2014
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             E_ERR, error_handler, nc_check, file_exist, logfileunit
use        model_mod, only : get_model_size, cable_state_to_dart_vector, &
                             get_cable_restart_filename, static_init_model
use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart
use time_manager_mod, only : time_type, print_time, print_date, set_date, &
                             get_time, operator(-)

use typesizes
use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=128)    :: cable_to_dart_output_file  = 'dart_ics'

namelist /cable_to_dart_nml/ cable_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)
character(len=256)    :: cable_restart_filename

!======================================================================

call initialize_utilities(progname='cable_to_dart')

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the cable namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "cable_to_dart_nml", iunit)
read(iunit, nml = cable_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "cable_to_dart_nml") ! closes, too.

call get_cable_restart_filename( cable_restart_filename )

write(*,*)
write(*,'(''cable_to_dart:converting cable restart file '',A, &
      &'' to DART file '',A)') &
       trim(cable_restart_filename), trim(cable_to_dart_output_file)

write(logfileunit,*)
write(logfileunit,'(''cable_to_dart:converting cable restart file '',A, &
      &'' to DART file '',A)') &
       trim(cable_restart_filename), trim(cable_to_dart_output_file)

!----------------------------------------------------------------------
! get to work
!----------------------------------------------------------------------

x_size = get_model_size()
allocate(statevector(x_size))

call cable_state_to_dart_vector(cable_restart_filename, statevector, model_time) 

iunit = open_restart_write(cable_to_dart_output_file)

call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

call print_date(model_time, str='cable_to_dart:CABLE model date')
call print_time(model_time, str='cable_to_dart:DART  model time')
call print_date(model_time, str='cable_to_dart:CABLE model date',iunit=logfileunit)
call print_time(model_time, str='cable_to_dart:DART  model time',iunit=logfileunit)

call finalize_utilities('cable_to_dart')

end program cable_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
