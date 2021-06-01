! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program clm_to_dart

!-------------------------------------------------------------------------------
! purpose: interface between clm and DART
!
! method: Read CLM restart files of model state
!         Reform variables into a netCDF file with DART missing values
!         in all the right places. CLM restart files have certain
!         variables that - despite having a missing_value or _FillValue
!         attribute - may or may not actually USE the values.
! 
! USAGE:  The clm filename is read from the clm_in namelist
!         <edit clm_to_dart_output_file in input.nml:clm_to_dart_nml>
!         clm_to_dart
!
! author: Tim Hoar 12 July 2011
!         and again 26 May 2021
!-------------------------------------------------------------------------------

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read

use        model_mod, only : get_model_size, mark_missing_r8_values, &
                             read_model_time, static_init_model

use time_manager_mod, only : time_type, print_time, print_date

implicit none

character(len=*), parameter :: source = 'clm_to_dart.f90'

!-------------------------------------------------------------------------------
! namelist parameters ... no default values.
!-------------------------------------------------------------------------------
! There can be at most three input files - restart, history, vector history

character(len=256) :: clm_to_dart_input_files(3) = 'null'
character(len=256) :: clm_to_dart_output_files(3) = 'null'

namelist /clm_to_dart_nml/ clm_to_dart_input_files, clm_to_dart_output_files

!-------------------------------------------------------------------------------
! global storage
!-------------------------------------------------------------------------------

integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)

!===============================================================================

call initialize_utilities(progname='clm_to_dart')

call static_init_model()

! Read the namelist to get the filenames.

call find_namelist_in_file("input.nml", "clm_to_dart_nml", iunit)
read(iunit, nml = clm_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "clm_to_dart_nml") ! closes, too.

! Reads the valid time, the state, and the target time.

model_time = read_model_time(clm_to_dart_output_file)

! Each variable specifies its own file of origin.
call mark_missing_r8_values(clm_to_dart_output_file, model_time) 

call write_clm_state_to_netcdf()

call print_date(model_time, str='clm_to_dart:clm  model date')
call print_time(model_time, str='clm_to_dart:DART model time')

call finalize_utilities('clm_to_dart')

end program clm_to_dart

