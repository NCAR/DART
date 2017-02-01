! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_naaps

!----------------------------------------------------------------------
! purpose: interface between DART and the naaps model
!
! method: Read DART state vector and overwrite values in a naaps restart file.
!         If the DART state vector has an 'advance_to_time' present, a
!         file called naaps_in.DART is created with a date-time-group
!         appropriate to advance naaps to the requested time.
!
!         The dart_to_naaps_nml namelist setting for advance_time_present 
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities,   &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, print_time, print_date, get_time, get_date
use        model_mod, only : static_init_model, statevector_to_analysis_file, &
                             get_model_size, get_naaps_restart_path,          &
                             get_naaps_ensemble_member 

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=128) :: dart_to_naaps_input_file = 'dart_restart'
logical            :: advance_time_present     = .true.
logical            :: verbose                  = .false.
namelist /dart_to_naaps_nml/ dart_to_naaps_input_file, &
                             advance_time_present,     &
                             verbose

!----------------------------------------------------------------------
! Local variables
!----------------------------------------------------------------------

integer               :: iunit, io, x_size, member
type(time_type)       :: model_time, adv_to_time
real(r8), allocatable :: statevector(:)
character(len=256)    :: naaps_restart_path
integer               :: iyear,imonth,iday,ihour,imin,isec,dtg

!======================================================================

call initialize_utilities(progname='dart_to_naaps', output_flag=verbose)

!----------------------------------------------------------------------
! Read the namelist to get the output filename, etc.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "dart_to_naaps_nml", iunit)
read(iunit, nml = dart_to_naaps_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_naaps_nml")

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the naaps namelists
! to set grid sizes, restart path, etc.
!----------------------------------------------------------------------

call static_init_model()
call get_naaps_restart_path( naaps_restart_path )
call get_naaps_ensemble_member( member )

write(*,*)
write(*,'(''dart_to_naaps:converting DART file '',a, &
      &'' to naaps restart dir '',a,'' member '',i4)') &
     trim(dart_to_naaps_input_file), trim(naaps_restart_path), member

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

x_size = get_model_size()
allocate(statevector(x_size))

iunit = open_restart_read(dart_to_naaps_input_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! update the current naaps state vector
! write out a tiny file with the advance_to_dtg information ... mebbe
!----------------------------------------------------------------------

call statevector_to_analysis_file(statevector, naaps_restart_path, member)

if ( advance_time_present ) then
   call get_date(adv_to_time,iyear,imonth,iday,ihour,imin,isec)
   dtg = iyear*1000000 + imonth*10000 + iday*100 + ihour
   iunit = open_file('adv_to_dtg', action='write')
   write(iunit,*) dtg
   call close_file(iunit)
endif

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_naaps:NAAPS model date')
call print_time( model_time,'dart_to_naaps:DART  model time')
call print_date( model_time,'dart_to_naaps:NAAPS model date',logfileunit)
call print_time( model_time,'dart_to_naaps:DART  model time',logfileunit)

if ( advance_time_present ) then
call print_time(adv_to_time,'dart_to_naaps:advance_to time')
call print_date(adv_to_time,'dart_to_naaps:advance_to date')
call print_time(adv_to_time,'dart_to_naaps:advance_to time',logfileunit)
call print_date(adv_to_time,'dart_to_naaps:advance_to date',logfileunit)
endif

call finalize_utilities()
end program dart_to_naaps

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
