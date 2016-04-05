! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_gitm

!----------------------------------------------------------------------
! purpose: interface between DART and the GITM model
!
! method: Read DART state vector and overwrite values in a gitm restart file.
!         If the DART state vector has an 'advance_to_time' present, a
!         file called DART_GITM_time_control.txt is created with
!         information appropriate to advance gitm to the requested time.
!
!         The dart_to_gitm_nml namelist setting for advance_time_present
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!----------------------------------------------------------------------

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file

use        model_mod, only : static_init_model, get_model_size, &
                             get_gitm_restart_dirname, statevector_to_restart_file

use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart

use time_manager_mod, only : time_type, print_time, print_date, operator(-), &
                             get_time, get_date

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character (len = 128) :: dart_to_gitm_input_file = 'dart_restart'
logical               :: advance_time_present    = .false.

namelist /dart_to_gitm_nml/ dart_to_gitm_input_file, &
                            advance_time_present

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: iunit, io, x_size
type(time_type)       :: model_time, adv_to_time
real(r8), allocatable :: statevector(:)
character(len=256)    :: gitm_restart_dirname  = 'none'

!======================================================================

call initialize_utilities(progname='dart_to_gitm')

!----------------------------------------------------------------------
! Read the namelist.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "dart_to_gitm_nml", iunit)
read(iunit, nml = dart_to_gitm_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_gitm_nml")

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the gitm namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()
call get_gitm_restart_dirname(gitm_restart_dirname)

write(*,*)
write(*,*) 'dart_to_gitm: converting DART file ', "'"//trim(dart_to_gitm_input_file)//"'"
write(*,*) 'to gitm restart files in directory ', "'"//trim(gitm_restart_dirname)//"'"

x_size = get_model_size()
allocate(statevector(x_size))

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_to_gitm_input_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! update the current gitm state vector
! Convey the amount of time to integrate the model ...
! time_manager_nml: stop_option, stop_count increments
!----------------------------------------------------------------------

call statevector_to_restart_file(statevector, gitm_restart_dirname, model_time)

if ( advance_time_present ) then
   call write_gitm_time_control(model_time, adv_to_time)
endif

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_gitm:gitm model date')
call print_time( model_time,'dart_to_gitm:DART model time')
call print_date( model_time,'dart_to_gitm:gitm model date',logfileunit)
call print_time( model_time,'dart_to_gitm:DART model time',logfileunit)

if ( advance_time_present ) then
call print_time(adv_to_time,'dart_to_gitm:advance_to time')
call print_date(adv_to_time,'dart_to_gitm:advance_to date')
call print_time(adv_to_time,'dart_to_gitm:advance_to time',logfileunit)
call print_date(adv_to_time,'dart_to_gitm:advance_to date',logfileunit)
endif

! end - close the log, etc
call finalize_utilities()

!======================================================================
contains
!======================================================================

subroutine write_gitm_time_control(model_time, adv_to_time)
! The idea is to write a text file with the following structure:
!
!#TIMESTART
!2003            year
!06              month
!21              day
!00              hour
!00              minute
!00              second
!
!#TIMEEND
!2003            year
!07              month
!21              day
!00              hour
!00              minute
!00              second
!

type(time_type), intent(in) :: model_time, adv_to_time
integer :: iyear,imonth,iday,ihour,imin,isec

iunit = open_file('DART_GITM_time_control.txt', action='write')
write(iunit,*)

! the end time comes first.

call get_date(adv_to_time,iyear,imonth,iday,ihour,imin,isec)
write(iunit,'(''#TIMEEND'')')
write(iunit,'(i4.4,10x,''year''  )')iyear
write(iunit,'(i2.2,12x,''month'' )')imonth
write(iunit,'(i2.2,12x,''day''   )')iday
write(iunit,'(i2.2,12x,''hour''  )')ihour
write(iunit,'(i2.2,12x,''minute'')')imin
write(iunit,'(i2.2,12x,''second'')')isec
write(iunit,*)

call get_date(model_time,iyear,imonth,iday,ihour,imin,isec)
write(iunit,'(''#TIMESTART'')')
write(iunit,'(i4.4,10x,''year''  )')iyear
write(iunit,'(i2.2,12x,''month'' )')imonth
write(iunit,'(i2.2,12x,''day''   )')iday
write(iunit,'(i2.2,12x,''hour''  )')ihour
write(iunit,'(i2.2,12x,''minute'')')imin
write(iunit,'(i2.2,12x,''second'')')isec
write(iunit,*)

call close_file(iunit)
end subroutine write_gitm_time_control



end program dart_to_gitm

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
