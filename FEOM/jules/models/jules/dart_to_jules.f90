! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!-----------------------------------------------------------------------
!> purpose: interface between DART and the jules model
!> 
!> method: Read DART state vector and overwrite values in a jules restart file.
!>         If the DART state vector has an 'advance_to_time' present, 
!>         it is read ... but nothing happens with it at this time.
!>         DART is NEVER expected to advance jules.
!>
!>         The dart_to_jules_nml namelist setting for advance_time_present 
!>         determines whether or not the input file has an 'advance_to_time'.
!>         Typically, only temporary files like 'assim_model_state_ic' have
!>         an 'advance_to_time'.
!> 
!> author: Tim Hoar 10 August 2015

program dart_to_jules

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file, &
                             error_handler, E_MSG, E_ERR, get_unit

use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart

use time_manager_mod, only : time_type, print_time, print_date, get_time, get_date

use        model_mod, only : static_init_model, dart_to_jules_restart, &
                             get_model_size, get_jules_restart_filename

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=256) :: dart_to_jules_input_file = 'dart_restart'
logical            :: advance_time_present     = .false.

namelist /dart_to_jules_nml/ dart_to_jules_input_file, &
                           advance_time_present

!----------------------------------------------------------------------

character(len=256)    :: jules_restart_filename
integer               :: iunit, io, x_size
type(time_type)       :: model_time, adv_to_time
real(r8), allocatable :: statevector(:)

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_jules')

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the jules namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

x_size = get_model_size()
allocate(statevector(x_size))

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "dart_to_jules_nml", iunit)
read(iunit, nml = dart_to_jules_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_jules_nml")

call get_jules_restart_filename( jules_restart_filename )

write(*,*)
write(*,'(''dart_to_jules:converting DART file '',A, &
      &'' to jules restart file '',A)') &
     trim(dart_to_jules_input_file), trim(jules_restart_filename)

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_to_jules_input_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! update the current jules state vector
!----------------------------------------------------------------------

call dart_to_jules_restart( statevector )

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_jules:JULES model date')
call print_time( model_time,'dart_to_jules:DART  model time')
call print_date( model_time,'dart_to_jules:JULES model date',logfileunit)
call print_time( model_time,'dart_to_jules:DART  model time',logfileunit)

if ( advance_time_present ) then
   call print_time(adv_to_time,'dart_to_jules:advance_to time')
   call print_date(adv_to_time,'dart_to_jules:advance_to date')
   call print_time(adv_to_time,'dart_to_jules:advance_to time',logfileunit)
   call print_date(adv_to_time,'dart_to_jules:advance_to date',logfileunit)

   call write_jules_control_file('DART_time_control.txt',model_time, adv_to_time)
endif

call finalize_utilities('dart_to_jules')

!======================================================================
contains
!======================================================================


subroutine write_jules_control_file(filename, model_time, adv_to_time)

character(len=*), intent(in) :: filename
type(time_type),  intent(in) :: model_time
type(time_type),  intent(in) :: adv_to_time

integer :: file_unit, io
integer :: year, month, day, hour, minute, second
character(len=64) :: timestring

! Write updated JULES namelist variables to a text file.
! It is up to advance_model.csh to update the JULES namelist.
! Write the information to a text file so we can grep the desired strings and
! then update the right parts of the timesteps.nml - without having to write the
! whole timesteps.nml. 

file_unit = get_unit()
open(unit = file_unit, file = trim(filename))

call get_date(model_time, year, month, day, hour, minute, second)
write(timestring,100) year, month, day, hour, minute, second

write(*,*)' model time timestring is ',trim(timestring)

write(file_unit, *, iostat=io) 'main_run_start = "',trim(timestring),'"'
if (io /= 0) call error_handler(E_ERR,'dart_to_model:', &
   'cannot write main_run_start to '//trim(filename),source,revision,revdate)

call get_date(adv_to_time, year, month, day, hour, minute, second)
write(timestring,100) year, month, day, hour, minute, second

write(file_unit, *, iostat=io) 'main_run_end   = "',trim(timestring),'"'
if (io /= 0) call error_handler(E_ERR,'dart_to_model:', &
   'cannot write main_run_end to '//trim(filename),source,revision,revdate)

close(file_unit)

100 format (i4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',i2.2)

end subroutine write_jules_control_file


end program dart_to_jules

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
