! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program ncommas_to_dart

!----------------------------------------------------------------------
! purpose: interface between ncommas and DART
!
! method: Read ncommas "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The ncommas filename is read from the ncommas_in namelist
!         <edit ncommas_to_dart_output_file in input.nml:ncommas_to_dart_nml>
!         ncommas_to_dart
!
! author: Tim Hoar 6/24/09
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, error_handler, E_MSG, &
                             find_namelist_in_file, check_namelist_read, &
                             finalize_utilities
use        model_mod, only : get_model_size, restart_file_to_sv, &
                             get_ncommas_restart_filename
use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart
use time_manager_mod, only : time_type, print_time, print_date

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=128) :: ncommas_to_dart_output_file  = 'dart_ics'

namelist /ncommas_to_dart_nml/ ncommas_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

logical               :: verbose = .TRUE.
integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)
character(len=256)    :: ncommas_restart_filename

!======================================================================

call initialize_utilities(progname='ncommas_to_dart', output_flag=verbose)

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "ncommas_to_dart_nml", iunit)
read(iunit, nml = ncommas_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "ncommas_to_dart_nml") ! closes, too.

write(*,*)
write(*,'(''ncommas_to_dart:converting ncommas restart file '',A, &
      &'' to DART file '',A)') &
       trim(ncommas_restart_filename), trim(ncommas_to_dart_output_file)

!----------------------------------------------------------------------
! get to work
!----------------------------------------------------------------------

x_size = get_model_size()
allocate(statevector(x_size))

call get_ncommas_restart_filename( ncommas_restart_filename )

call restart_file_to_sv(ncommas_restart_filename, statevector, model_time) 

iunit = open_restart_write(ncommas_to_dart_output_file)

call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

call print_date(model_time, str='ncommas_to_dart:ncommas  model date')
call print_time(model_time, str='ncommas_to_dart:DART model time')

call error_handler(E_MSG,'ncommas_to_dart','Finished successfully.',source,revision,revdate)
call finalize_utilities()

end program ncommas_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
