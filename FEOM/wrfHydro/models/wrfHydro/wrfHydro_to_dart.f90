! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program wrfHydro_to_dart

!----------------------------------------------------------------------
! purpose: interface between wrfHydro and DART
!
! method: Read the 2 wrfHydro "restart" files (LSM = noah | noahMP, HYDRO) 
!         of model state.
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
!
! USAGE:  The filename is read from the noah_in namelist
!         wrfHydro_to_dart
!
! author: Tim Hoar 11 July 2012
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit
use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart
use time_manager_mod, only : time_type, print_time, print_date
use        model_mod, only : static_init_model, get_model_size, model_to_dart_vector, &
                             get_lsm_restart_filename, get_hydro_restart_filename, &
                             get_assimOnly_restart_filename, &
                             get_debug_level

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------
character(len=128) :: wrfHydro_to_dart_output_file  = 'dart_ics'
namelist /wrfHydro_to_dart_nml/ wrfHydro_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------
integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)
character(len=256)    :: lsm_restart_filename, hydro_restart_filename, assimOnly_restart_filename

!======================================================================
call initialize_utilities(progname='wrfHydro_to_dart')

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------
call static_init_model()

! Read the namelist to get the input filename.

call find_namelist_in_file("input.nml", "wrfHydro_to_dart_nml", iunit)
read(iunit, nml = wrfHydro_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "wrfHydro_to_dart_nml")

! the output filename comes from the initialization of model_mod

call get_lsm_restart_filename( lsm_restart_filename )
call get_hydro_restart_filename( hydro_restart_filename )
call get_assimOnly_restart_filename( assimOnly_restart_filename )

!----------------------------------------------------------------------
! Read the restart and return the DART state vector.
!----------------------------------------------------------------------
write(*,*)
write(*,'(''wrfHydro_to_dart:converting restart files <'',A, &
      &'' & '', A, ''> to DART file <'',A,''>'')') &
       trim(lsm_restart_filename), trim(hydro_restart_filename), trim(wrfHydro_to_dart_output_file)
write(logfileunit,*)
write(logfileunit,'(''wrfHydro_to_dart:converting noah restart file <'',A, &
      &'' & '', A, ''> to DART file <'',A,''>'')') &
       trim(lsm_restart_filename), trim(hydro_restart_filename), trim(wrfHydro_to_dart_output_file)

x_size = get_model_size()
allocate(statevector(x_size))

call model_to_dart_vector(lsm_restart_filename, &
                          hydro_restart_filename, &
                          assimOnly_restart_filename, &
                          statevector, model_time)

! Output the DART state vector (apriori)
iunit = open_restart_write(wrfHydro_to_dart_output_file)
call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

!----------------------------------------------------------------------
! call finalize_utilities()
!----------------------------------------------------------------------

if (get_debug_level() > 0) then
   call print_date(model_time, str='wrfHydro_to_dart:DART model date',iunit=logfileunit)
   call print_date(model_time, str='wrfHydro_to_dart:DART model date')
   call print_time(model_time, str='wrfHydro_to_dart:DART model time')
   call print_time(model_time, str='wrfHydro_to_dart:DART model time',iunit=logfileunit)
endif

call finalize_utilities()

end program wrfHydro_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
