! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program dart_to_cosmo

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between DART and the cosmo model
!
! method: Read DART state vector and overwrite values in a cosmo restart file.
!         If the DART state vector has an 'advance_to_time' present, a
!         file called cosmo_in.DART is created with a time_manager_nml namelist 
!         appropriate to advance cosmo to the requested time.
!
!         The dart_to_cosmo_nml namelist setting for advance_time_present 
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, print_time, print_date, operator(-), get_time, set_time
use        model_mod, only : static_init_model, write_grib_file, get_model_size, write_state_times

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 128) :: dart_output_file        = 'dart.ic'
character (len = 128) :: new_cosmo_analysis_file = 'out.grb'
logical               :: advance_time_present    = .true.

namelist /dart_to_cosmo_nml/ dart_output_file,        &
                             new_cosmo_analysis_file, &
                             advance_time_present

!----------------------------------------------------------------------

integer               :: iunit, io, model_size
type(time_type)       :: model_time, adv_to_time
real(r8), allocatable :: state_vector(:)
logical               :: verbose = .false.

!----------------------------------------------------------------------
! Read the namelist to get the output filename. 
!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_cosmo', output_flag=verbose)

call find_namelist_in_file("input.nml", "dart_to_cosmo_nml", iunit)
read(iunit, nml = dart_to_cosmo_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_cosmo_nml")

write(*,*)
write(*,'(''dart_to_cosmo:converting DART file "'',A, &
      &''" to cosmo file "'',A,''"'')') &
     trim(dart_output_file), trim(new_cosmo_analysis_file)

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the cosmo namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

model_size=get_model_size()

allocate(state_vector(1:model_size))

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_output_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, state_vector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, state_vector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! update the current cosmo state vector
! Convey the amount of time to integrate the model ...
! time_manager_nml: stop_option, stop_count increments
!----------------------------------------------------------------------

call write_grib_file(state_vector, new_cosmo_analysis_file)

if ( advance_time_present ) then
   iunit = open_file('times', action='write')
   call write_state_times(iunit, model_time, adv_to_time)
   call close_file(iunit)
endif

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_cosmo:cosmo model date')
call print_time( model_time,'dart_to_cosmo:DART model time')
call print_date( model_time,'dart_to_cosmo:cosmo model date',logfileunit)
call print_time( model_time,'dart_to_cosmo:DART model time',logfileunit)

if ( advance_time_present ) then
call print_time(adv_to_time,'dart_to_cosmo:advance_to time')
call print_date(adv_to_time,'dart_to_cosmo:advance_to date')
call print_time(adv_to_time,'dart_to_cosmo:advance_to time',logfileunit)
call print_date(adv_to_time,'dart_to_cosmo:advance_to date',logfileunit)
endif

call finalize_utilities()

end program dart_to_cosmo
