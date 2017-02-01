! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program model_to_dart

!----------------------------------------------------------------------
! purpose: interface between model and DART
!
! method: Read MPAS "history" files of model state.
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The model filename is read from the model_in namelist
!         <edit model_to_dart_output_file in input.nml:model_to_dart_nml>
!         model_to_dart
!
! author: Tim Hoar 12 Sep 2011
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit
use        model_mod, only : get_model_size, analysis_file_to_statevector, &
                             get_model_analysis_filename, static_init_model, &
                             print_variable_ranges
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

character(len=128) :: model_to_dart_output_file  = 'dart_ics'
logical            :: print_data_ranges          = .true.

namelist /model_to_dart_nml/    &
     model_to_dart_output_file, &
     print_data_ranges

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)
character(len=256)    :: model_analysis_filename

!======================================================================

call initialize_utilities(progname='model_to_dart')

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "model_to_dart_nml", iunit)
read(iunit, nml = model_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "model_to_dart_nml") ! closes, too.

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the model namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------
call static_init_model()

call get_model_analysis_filename(model_analysis_filename)

x_size = get_model_size()
allocate(statevector(x_size))

write(*,*)
write(*,*) 'model_to_dart: converting model analysis file ', &
           "'"//trim(model_analysis_filename)//"'" 
write(*,*) ' to DART file ', "'"//trim(model_to_dart_output_file)//"'"

!----------------------------------------------------------------------
! Read the valid time and the state from the MPAS netcdf file
!----------------------------------------------------------------------
call analysis_file_to_statevector(model_analysis_filename, statevector, model_time) 

!----------------------------------------------------------------------
! if requested, print out the data ranges variable by variable
!----------------------------------------------------------------------
if (print_data_ranges) then
    call print_variable_ranges(statevector)
endif

!----------------------------------------------------------------------
! Write the valid time and the state to the dart restart file
!----------------------------------------------------------------------
iunit = open_restart_write(model_to_dart_output_file)

call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

!----------------------------------------------------------------------
! finish up
!----------------------------------------------------------------------

call print_date(model_time, 'model_to_dart:model date')
call print_time(model_time, 'model_to_dart:model time')
call print_date(model_time, 'model_to_dart:model date',logfileunit)
call print_time(model_time, 'model_to_dart:model time',logfileunit)

call finalize_utilities()

end program model_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
