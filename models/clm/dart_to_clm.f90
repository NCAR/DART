! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_clm

!----------------------------------------------------------------------
! purpose: interface between DART and the CLM model
!
! method: Read DART state vector and overwrite values in a CLM restart file.
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file, &
                             error_handler, E_MSG
use time_manager_mod, only : time_type, print_time, print_date, operator(-), get_time
use        model_mod, only : static_init_model, fill_missing_r8_with_orig, &
                             get_model_size, get_clm_restart_filename, &
                             read_model_time

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=256) :: dart_to_clm_input_file = 'model_restart.nc'
logical            :: advance_time_present   = .false.

namelist /dart_to_clm_nml/ dart_to_clm_input_file, &
                           advance_time_present

!----------------------------------------------------------------------

character(len=256) :: clm_restart_filename
integer            :: iunit, io
type(time_type)    :: model_time, adv_to_time

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_clm')

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the clm namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

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

model_time = read_model_time(clm_restart_filename)

!----------------------------------------------------------------------
! update the current clm state vector
!----------------------------------------------------------------------

call fill_missing_r8_with_orig(clm_restart_filename, dart_to_clm_input_file, model_time)

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
