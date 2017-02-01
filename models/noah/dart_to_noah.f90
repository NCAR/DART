! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_noah

!----------------------------------------------------------------------
! purpose: interface between DART and the NOAH model
!
! method: Read DART state vector and overwrite values in a noah restart file.
!         If the DART state vector has an 'advance_to_time' present,
!         it is read ... but nothing happens with it at this time.
!         DART is NEVER expected to advance noah.
!
!         The dart_to_noah_nml namelist setting for advance_time_present
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!
! author: Tim Hoar 12 July 2011
!----------------------------------------------------------------------

use        types_mod, only : r8, obstypelength
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file, &
                             error_handler, E_ERR
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, print_time, print_date, get_date, &
                             set_time, operator(+), operator(-), operator(<=)
use        model_mod, only : static_init_model, dart_vector_to_model_file, &
                             get_model_size, get_noah_restart_filename, &
                             get_noah_timestepping, get_debug_level

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character (len = 128) :: dart_to_noah_input_file = 'dart_restart'
character (len=obstypelength), dimension(40) :: skip_variables = ' '
logical               :: advance_time_present   = .true.

namelist /dart_to_noah_nml/ dart_to_noah_input_file, &
                            skip_variables, &
                            advance_time_present

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

character(len=20)     :: noah_restart_filename
integer               :: ifile, nfiles, iunit, io, x_size
type(time_type)       :: model_time, adv_to_time, mytime
type(time_type)       :: forcingtimestep, nexttimestep
real(r8), allocatable :: statevector(:)
integer               :: kday, khour, noah_timestep, output_timestep
integer               :: forcing_timestep, restart_frequency_seconds

integer :: year,month,day,hour,minute,second
character(len=32)     :: datestring
character(len=128)    :: string1,string2,string3

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_noah')

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the NOAH namelist
! to set location and state vector
!----------------------------------------------------------------------

call static_init_model()

! Read the namelist to get the input filename.

call find_namelist_in_file("input.nml", "dart_to_noah_nml", iunit)
read(iunit, nml = dart_to_noah_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_noah_nml")

! the output filename comes from the initialization of model_mod

call get_noah_restart_filename( noah_restart_filename )

write(*,*)
write(*,'(''dart_to_noah:converting DART file <'',A, &
      &''> to NOAH restart file <'',A,''>'')') &
     trim(dart_to_noah_input_file), trim(noah_restart_filename)
write(logfileunit,*)
write(logfileunit,'(''dart_to_noah:converting DART file <'',A, &
      &''> to NOAH restart file <'',A,''>'')') &
     trim(dart_to_noah_input_file), trim(noah_restart_filename)

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

x_size = get_model_size()
allocate(statevector(x_size))

iunit = open_restart_read(dart_to_noah_input_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! write the updated state to the NOAH restart file.
!----------------------------------------------------------------------

call dart_vector_to_model_file(statevector, noah_restart_filename, model_time, &
                               skip_variables)

!----------------------------------------------------------------------
! Convey adv_to_time to noah by updating kday or khour in the namelist.
! Predict the names of the LDASIN files needed -
!    write these to a file that will be queried by advance_model.csh
!    so they get staged appropriately.
!----------------------------------------------------------------------

if ( advance_time_present ) then

   ! AS long as the output_timestep, forcing_timestep, and noah_timestep are all identical,
   ! the time in the restart file is the time of the next forcing file needed.

   call get_noah_timestepping(kday, khour, noah_timestep, output_timestep, forcing_timestep, &
                              restart_frequency_seconds)

   if ( (noah_timestep == forcing_timestep) .and. &
        (noah_timestep == output_timestep )) then

   else
      write(string2,*) 'noah_timestep ', noah_timestep, ' /= forcing_timestep of ',forcing_timestep
      write(string3,*) 'noah_timestep ', noah_timestep, ' /= output_timestep of ',output_timestep
      write(string1,*) 'temporarily unsupported configuration'
      call error_handler(E_ERR,'dart_to_noah',string1,source,revision,revdate,&
                         text2=string2,text3=string3)
   endif

   ! figure out how many timesteps between adv_to_time and model_time
   mytime          = model_time
   nexttimestep    = set_time(   noah_timestep, 0)
   forcingtimestep = set_time(forcing_timestep, 0)
   nfiles = 0
   TIMELOOP: do while (mytime <= adv_to_time)
      nfiles = nfiles + 1
      mytime = mytime + nexttimestep
   enddo TIMELOOP

   if (get_debug_level() > 0) &
   write(*,*)'needed ',nfiles,' LDASIN files to get from model_time to adv_to_time.'

   iunit = open_file('noah_advance_information.txt',form='formatted',action='write')
   call print_date(  model_time,'dart_to_noah:noah  model      date',iunit)
   call print_date( adv_to_time,'dart_to_noah:noah  advance_to_date',iunit)
   write(iunit,'(''khour  = '',i6)') nfiles-1
   write(iunit,'(''nfiles = '',i6)') nfiles

   mytime = model_time

   do ifile = 1,nfiles
      ! The model time is one noah_timestep (which must be equal to RESTART_FREQUENCY_HOURS)
      ! behind the file names needed.

      mytime = mytime + nexttimestep

      call get_date(mytime,year,month,day,hour,minute,second)

      ! TJH FIXME - what happens to seconds ...
      if ((minute == 0) .and. (second == 0)) then
         write(datestring,'(i4.4,3(i2.2))')year,month,day,hour
      else
         write(datestring,'(i4.4,4(i2.2))')year,month,day,hour,minute
      endif

      write(iunit,'(a)')trim(datestring)//'.LDASIN_DOMAIN1'
   enddo

   call close_file(iunit)

endif

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

if (get_debug_level() > 0) then
   call print_date( model_time,'dart_to_noah:noah model date')
   call print_time( model_time,'dart_to_noah:DART model time')
   call print_date( model_time,'dart_to_noah:noah model date',logfileunit)
   call print_time( model_time,'dart_to_noah:DART model time',logfileunit)

   if ( advance_time_present ) then
      call print_time(adv_to_time,'dart_to_noah:advance_to time')
      call print_date(adv_to_time,'dart_to_noah:advance_to date')
      call print_time(adv_to_time,'dart_to_noah:advance_to time',logfileunit)
      call print_date(adv_to_time,'dart_to_noah:advance_to date',logfileunit)
   endif
endif

call finalize_utilities()

end program dart_to_noah

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
