! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_model

!----------------------------------------------------------------------
! purpose: interface between DART and the model model
!
! method: Read DART state vector and overwrite values in a model restart file.
!         If the DART state vector has an 'advance_to_time' present, a
!         file called model_in.DART is created with a time_manager_nml namelist 
!         appropriate to advance model to the requested time.
!
!         The dart_to_model_nml namelist setting for advance_time_present 
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!
! author: Tim Hoar 25 Jun 09, revised 12 July 2010
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, print_time, print_date, operator(-), &
                             get_time, get_date
use        model_mod, only : static_init_model, dart_vector_to_model_file, &
                             get_model_size

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 128) :: dart_to_model_input_file = 'dart_restart'
logical               :: advance_time_present     = .false.
character(len=256)    :: model_restart_filename   = 'model_restartfile'

namelist /dart_to_model_nml/ dart_to_model_input_file, &
                            advance_time_present,    &
                            model_restart_filename

!----------------------------------------------------------------------

integer               :: iunit, io, x_size, diff1, diff2
type(time_type)       :: model_time, adv_to_time, base_time
real(r8), allocatable :: statevector(:)
logical               :: verbose              = .FALSE.

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_model', output_flag=verbose)

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the model namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

x_size = get_model_size()
allocate(statevector(x_size))

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "dart_to_model_nml", iunit)
read(iunit, nml = dart_to_model_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_model_nml")

write(*,*)
write(*,*) 'dart_to_model: converting DART file ', "'"//trim(dart_to_model_input_file)//"'"
write(*,*) 'to model restart files named        ', "'"//trim(model_restart_filename)//"'" 

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_to_model_input_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

print *, 'read state vector'
!----------------------------------------------------------------------
! update the current model state vector
! Convey the amount of time to integrate the model ...
! time_manager_nml: stop_option, stop_count increments
!----------------------------------------------------------------------

print *, 'calling sv to restart file'
call dart_vector_to_model_file(statevector, model_restart_filename, model_time)

if ( advance_time_present ) then
   call write_model_time_control(model_time, adv_to_time)
endif

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_model:model model date')
call print_time( model_time,'dart_to_model:DART model time')
call print_date( model_time,'dart_to_model:model model date',logfileunit)
call print_time( model_time,'dart_to_model:DART model time',logfileunit)

if ( advance_time_present ) then
call print_time(adv_to_time,'dart_to_model:advance_to time')
call print_date(adv_to_time,'dart_to_model:advance_to date')
call print_time(adv_to_time,'dart_to_model:advance_to time',logfileunit)
call print_date(adv_to_time,'dart_to_model:advance_to date',logfileunit)
endif

call finalize_utilities()

!======================================================================
contains
!======================================================================

subroutine write_model_time_control(model_time, adv_to_time)
! Write a text file that the model can use to figure out how
! far to run until.  Could be as simple as a text file containing: 
!    YYYYMMDD hh:mm:ss
! or for something a bit more complicated, here's another example:
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

iunit = open_file('DART_model_time_control.txt', action='write')
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

call get_date(adv_to_time,iyear,imonth,iday,ihour,imin,isec)
write(iunit,'(''#TIMEEND'')') 
write(iunit,'(i4.4,10x,''year''  )')iyear
write(iunit,'(i2.2,12x,''month'' )')imonth
write(iunit,'(i2.2,12x,''day''   )')iday
write(iunit,'(i2.2,12x,''hour''  )')ihour
write(iunit,'(i2.2,12x,''minute'')')imin
write(iunit,'(i2.2,12x,''second'')')isec
write(iunit,*)

call close_file(iunit)
end subroutine write_model_time_control


end program dart_to_model

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
