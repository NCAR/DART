! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_cable

!----------------------------------------------------------------------
! purpose: interface between DART and the CABLE model
!
! method: Read DART state vector and overwrite values in a CABLE restart file.
!         If the DART state vector has an 'advance_to_time' present, 
!         it is read ... but nothing happens with it at this time.
!
!         The dart_to_cable_nml namelist setting for advance_time_present 
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!
! author: Tim Hoar 20 February 2014
!----------------------------------------------------------------------

use        types_mod, only : r8, digits12

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file, nc_check, &
                             file_exist, error_handler, E_MSG, E_ERR

use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart

use time_manager_mod, only : time_type, print_time, print_date, get_time, set_date, &
                             set_time, operator(-), operator(+), operator(==)

use        model_mod, only : static_init_model, dart_vector_to_model_file, &
                             get_model_size, get_cable_restart_filename, &
                             get_time_origin

use typesizes
use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2, string3

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=256) :: dart_to_cable_input_file = 'dart_restart'
character(len=256) :: cable_time_file = 'CABLE_time_file.nc'
logical            :: advance_time_present   = .false.

namelist /dart_to_cable_nml/ dart_to_cable_input_file, &
                           cable_time_file, advance_time_present

!----------------------------------------------------------------------

character(len=256)    :: cable_restart_filename
integer               :: iunit, io, x_size, newunit
type(time_type)       :: model_time, adv_to_time
real(r8), allocatable :: statevector(:)
integer               :: kstart, kend

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_cable')

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the cable namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

x_size = get_model_size()
allocate(statevector(x_size))

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "dart_to_cable_nml", iunit)
read(iunit, nml = dart_to_cable_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_cable_nml")

call get_cable_restart_filename( cable_restart_filename )

write(*,*)
write(*,'(''dart_to_cable:converting DART file '',A, &
      &'' to cable restart file '',A)') &
     trim(dart_to_cable_input_file), trim(cable_restart_filename)

write(logfileunit,*)
write(logfileunit,'(''dart_to_cable:converting DART file '',A, &
      &'' to cable restart file '',A)') &
     trim(dart_to_cable_input_file), trim(cable_restart_filename)

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_to_cable_input_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)

   ! Determine the indices of the forcing file that are needed to run.
   ! Yingping: "There is a limitation to the values of kstart or kend. 
   ! Kstart has to be the first time step of a day and kend has to be the 
   ! last time step of a day. If a model time step is hourly, therefore 
   ! possible values are 1, 25, 37.. for kstart, and 24, 48, 72 .. for kend.

   call get_cable_time_array(model_time, adv_to_time, kstart, kend)

   newunit = open_file('timestep.txt',action='rewind')
   write(newunit,*)kstart, kend
   call print_date( model_time,'dart_to_cable:DART   model date',newunit)
   call print_date(adv_to_time,'dart_to_cable:DART desired date',newunit)
   call close_file(newunit)

else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! update the current cable state vector
!----------------------------------------------------------------------

call dart_vector_to_model_file(statevector, cable_restart_filename, model_time)

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_cable:CABLE model date')
call print_time( model_time,'dart_to_cable:DART  model time')
call print_date( model_time,'dart_to_cable:CABLE model date',logfileunit)
call print_time( model_time,'dart_to_cable:DART  model time',logfileunit)

if ( advance_time_present ) then
   call print_time(adv_to_time,'dart_to_cable:advance_to time')
   call print_date(adv_to_time,'dart_to_cable:advance_to date')
   call print_time(adv_to_time,'dart_to_cable:advance_to time',logfileunit)
   call print_date(adv_to_time,'dart_to_cable:advance_to date',logfileunit)
endif

call finalize_utilities('dart_to_cable')

deallocate(statevector)

!===============================================================================
contains
!===============================================================================

subroutine get_cable_time_array(time1,time2,indx1,indx2)

! The time control for CABLE is a simple pair of integers that
! index into the time array of the gspwfile 'time' variable.
! We have to read that time array to figure out how to start/stop
! the model advance.

type(time_type), intent(in)  :: time1, time2
integer,         intent(out) :: indx1, indx2

integer :: ncid, TimeDimID, VarID, ntimes, i
integer :: year, month, day, hour, minute, second
character(len=256) :: unitstring

real(digits12), dimension(:), allocatable :: time_array, seconds
integer, dimension(:), allocatable :: days

type(time_type) :: time_origin, curr_time, test_time

if ( .not. file_exist(cable_time_file) ) then
   write(string1,*) 'cannot open file ', trim(cable_time_file),' for reading.'
   call error_handler(E_ERR,'get_cable_time_array',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(cable_time_file), NF90_NOWRITE, ncid), &
              'get_cable_time_array','open '//trim(cable_time_file))

! Get the length of the time dimension
call nc_check(nf90_inq_dimid(ncid, 'time', TimeDimID), &
            'get_cable_time_array', 'inq_dimid '//trim(string3))
call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=ntimes), &
            'get_cable_time_array', 'inq_dimid '//trim(string3))

allocate(time_array(ntimes), days(ntimes), seconds(ntimes))

! read the time variable
call nc_check(nf90_inq_varid(ncid, 'time', VarID), &
            'get_cable_time_array', 'inq_varid '//trim(string3))
call nc_check(nf90_get_var(ncid, VarID, time_array), &
            'get_cable_time_array', 'get_var'//trim(string3))

! read the time units/base
if( nf90_inquire_attribute(    ncid, VarID, 'units') == NF90_NOERR )  then
   call nc_check( nf90_get_att(ncid, VarID, 'units' , unitstring), &
           'get_cable_time_array', 'get_att units '//trim(string2))
else
   write(string1,*) 'cannot get time base from <', trim(cable_time_file),'>.'
   call error_handler(E_ERR,'get_cable_time_array',string1,source,revision,revdate)
endif

read(unitstring,'(14x,i4,5(1x,i2))'),year,month,day,hour,minute,second
time_origin = set_date(year, month, day, hour, minute, second)

! make sure that the time array is just 'seconds' and not fractional days'
if (index(unitstring,'seconds since') > 0) then
   write(*,*)'Methinks we have seconds since'
   days    = 0
   seconds = time_array
elseif (index(unitstring,'days since') >0) then
   write(*,*)'Methinks we have fractional days'
   days    = nint(time_array)
   seconds = (time_array - real(days,digits12))*86400.0_digits12
else
   write(string1,*) 'unknown time base from <', trim(cable_time_file),'>.'
   call error_handler(E_ERR, 'get_cable_time_array', string1, &
                 source, revision, revdate, text2=trim(unitstring))
endif

indx1 = -1
indx2 = -1

! indx1 (AKA 'kstart') is the timestep AFTER the time that matches.
FindKstart : do i = 1,ntimes-1

   curr_time = set_time(nint(seconds(i)), days(i))
   test_time = curr_time + time_origin
   if (test_time == time1) then
      indx1 = i + 1
      exit FindKstart
   endif

enddo FindKstart

if (indx1 < 1) then
   call print_date(time1,'Cannot find starting date of',logfileunit)
   call print_date(time1,'Cannot find starting date of')
   write(string1,*) 'cannot find start date in', trim(cable_time_file),'>.'
   call error_handler(E_ERR,'get_cable_time_array',string1,source,revision,revdate)
endif

FindKend : do i = indx1,ntimes

   curr_time = set_time(nint(seconds(i)), days(i))
   test_time = curr_time + time_origin
   if (test_time == time2) then
      indx2 = i
      exit FindKend
   endif

enddo FindKend

if (indx2 < 1) then
   call print_date(time1,'Cannot find ending date of',logfileunit)
   call print_date(time1,'Cannot find ending date of')
   write(string1,*) 'cannot find end date in', trim(cable_time_file),'>.'
   call error_handler(E_ERR,'get_cable_time_array',string1,source,revision,revdate)
endif

deallocate(time_array, days, seconds)

end subroutine get_cable_time_array


end program dart_to_cable

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
