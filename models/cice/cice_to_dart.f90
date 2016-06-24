! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: cice_to_dart.f90 8565 2015-09-11 17:16:08Z hkershaw $

program cice_to_dart

!----------------------------------------------------------------------
! purpose: interface between CICE and DART
!
! method: Read CICE "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The CICE filename is read from the cice_in namelist
!         <edit cice_to_dart_output_file in input.nml:cice_to_dart_nml>
!         cice_to_dart
!
! author: C Bitz June 2016
! borrowed heavily from Tim Hoar 25 Jun 09, revised 12 July 2010
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read
use        model_mod, only : restart_file_to_sv, static_init_model, &
                             get_model_size, get_cice_restart_filename
use state_vector_io_mod, only : awrite_state_restart, open_restart_write, close_restart
use time_manager_mod, only : time_type, print_time, print_date

use netcdf
implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_cice/models/cice/cice_to_dart.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 8565 $"
character(len=128), parameter :: revdate  = "$Date: 2015-09-11 10:16:08 -0700 (Fri, 11 Sep 2015) $"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character (len = 128) :: cice_to_dart_output_file  = 'dart_ics'

namelist /cice_to_dart_nml/ cice_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)
character (len = 128) :: cice_restart_filename = 'no_cice_restart_filename' 

!----------------------------------------------------------------------

call initialize_utilities(progname='cice_to_dart')

!----------------------------------------------------------------------
! Call model_mod:static_init_model(), which reads the namelists
! to set calendar type, starting date, deltaT, etc.
!----------------------------------------------------------------------

call static_init_model()

!----------------------------------------------------------------------
! Read the namelist to get the input and output filenames.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "cice_to_dart_nml", iunit)
read(iunit, nml = cice_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "cice_to_dart_nml") ! closes, too.

call get_cice_restart_filename( cice_restart_filename )

write(*,*)
write(*,'(''cice_to_dart:converting CICE restart file '',A, &
      &'' to DART file '',A)') &
       trim(cice_restart_filename), trim(cice_to_dart_output_file)

!----------------------------------------------------------------------
! Now that we know the names, get to work.
!----------------------------------------------------------------------

x_size = get_model_size()
allocate(statevector(x_size))
call restart_file_to_sv(cice_restart_filename, statevector, model_time) 

iunit = open_restart_write(cice_to_dart_output_file)

call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

call print_date(model_time, str='cice_to_dart:CICE  model date')
call print_time(model_time, str='cice_to_dart:DART model time')
call finalize_utilities('cice_to_dart')

end program cice_to_dart

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_cice/models/cice/cice_to_dart.f90 $
! $Id: cice_to_dart.f90 8565 2015-09-11 17:16:08Z hkershaw $
! $Revision: 8565 $
! $Date: 2015-09-11 10:16:08 -0700 (Fri, 11 Sep 2015) $
