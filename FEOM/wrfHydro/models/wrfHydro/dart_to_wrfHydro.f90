! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_wrfHydro

!----------------------------------------------------------------------
! purpose: interface between DART and wrfHydro model
!
! method: Read DART state vector and overwrite values in the two wrfHydro restart files.
!         If the DART state vector has an 'advance_to_time' present,
!         it is read ... but nothing happens with it at this time.
!         DART is not expected to advance wrfHydro.
!
!         The dart_to_wrfHydro_nml namelist setting for advance_time_present
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
                             set_time, operator(+), operator(-), operator(<=), & 
                             operator(<)
use        model_mod, only : static_init_model, dart_vector_to_model_files, &
                             get_model_size, get_lsm_restart_filename, &
                             get_lsm_restart_filename, &
                             get_hydro_restart_filename, & 
                             get_assimOnly_restart_filename, & 
                             get_model_timestepping, get_debug_level

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character (len = 128) :: dart_to_wrfHydro_input_file = 'dart_restart'
character (len=obstypelength), dimension(40) :: skip_variables = ' '
logical               :: advance_time_present   = .true.

namelist /dart_to_wrfHydro_nml/ dart_to_wrfHydro_input_file, &
                            skip_variables, &
                            advance_time_present

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

character(len=20)     :: lsm_restart_filename, hydro_restart_filename, assimOnly_restart_filename
integer               :: ifile, nfiles, iunit, io, x_size
type(time_type)       :: model_time, adv_to_time, mytime
type(time_type)       :: forcingtimestep, nexttimestep
real(r8), allocatable :: statevector(:)
integer               :: kday, khour, noah_timestep, output_timestep
integer               :: forcing_timestep, restart_frequency_seconds

integer :: year,month,day,hour,minute,second
character(len=obstypelength) :: datestring
character(len=128)           :: string1,string2,string3

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_wrfHydro')

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the namelists
! to set location and state vector
!----------------------------------------------------------------------

call static_init_model()

! Read the namelist to get the input filename.

call find_namelist_in_file("input.nml", "dart_to_wrfHydro_nml", iunit)
read(iunit, nml = dart_to_wrfHydro_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_wrfHydro_nml")

! the output filename comes from the initialization of model_mod

call get_lsm_restart_filename( lsm_restart_filename )
call get_hydro_restart_filename( hydro_restart_filename )
call get_assimOnly_restart_filename( assimOnly_restart_filename )

write(*,*)
write(*,'(''dart_to_wrfHydro:converting DART file <'',A, &
      &''> to wrfHydro restart files <'',A,'' & '', A, ''>'')') &
      trim(dart_to_wrfHydro_input_file),trim(lsm_restart_filename), trim(hydro_restart_filename)

write(logfileunit,*)
write(logfileunit,'(''dart_to_wrfHydro:converting DART file <'',A, &
      &''> wrfHydro restart file <'',A,'' & '', A, ''>'')') &
      trim(dart_to_wrfHydro_input_file),trim(lsm_restart_filename), trim(hydro_restart_filename)


!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time. 
! JLM fixme - from DART? we need to be explicit about in/out of time 
! and which component it pertains to. model_time seems to be dart_time.
! I need to check this when I have time. 
!----------------------------------------------------------------------
x_size = get_model_size()
allocate(statevector(x_size))
iunit = open_restart_read(dart_to_wrfHydro_input_file)
if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time) 
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! write the updated state to the wrfHydro restart file.
!----------------------------------------------------------------------
call dart_vector_to_model_files(statevector, &
                               lsm_restart_filename, &
                               hydro_restart_filename, &
                               assimOnly_restart_filename, &
                               model_time, skip_variables)

!----------------------------------------------------------------------
! Convey adv_to_time to noah by updating kday or khour in the namelist.
! Predict the names of the LDASIN files needed -
!    write these to a file that will be queried by advance_model.csh
!    so they get staged appropriately.
!----------------------------------------------------------------------
if ( advance_time_present ) then

   ! AS long as the output_timestep, forcing_timestep, and noah_timestep are all identical,
   ! the time in the restart file is the time of the next forcing file needed.

   call get_model_timestepping(kday, khour, noah_timestep, output_timestep, forcing_timestep, &
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
   !! jlm fixme - mytime?  huh?
   nexttimestep    = set_time(   noah_timestep, 0)
   forcingtimestep = set_time(forcing_timestep, 0)
   mytime          = model_time
   nfiles = 0

   TIMELOOP: do while (mytime < adv_to_time)  
         nfiles = nfiles + 1
         mytime = mytime + nexttimestep
   enddo TIMELOOP

   if (get_debug_level() > 0) &
   write(*,*)'needed ',nfiles,' LDASIN files to get from model_time to adv_to_time.'

   iunit = open_file('wrfHydro_advance_information.txt',form='formatted',action='write')
   call print_date(  model_time,'dart_to_wrfHydro:wrfHydro  model      date',iunit)
   call print_date( adv_to_time,'dart_to_wrfHydro:wrfHydro  advance_to_date',iunit)
   write(iunit,'(''khour  = '',i6)') nfiles
   write(iunit,'(''nfiles = '',i6)') nfiles

   mytime = model_time  !! why???? my=hrldas, model=dart???? this needs to be explicit

   do ifile = 1,nfiles
      mytime = mytime + nexttimestep 

      call get_date(mytime,year,month,day,hour,minute,second)

      ! TJH FIXME - what happens to seconds ...
      ! we only seem to run at time resolutions of hours for wrfHydro.... jlm
      if ((minute == 0) .and. (second == 0)) then
         write(datestring,'(i4.4,3(i2.2))')year,month,day,hour
      else
         write(datestring,'(i4.4,4(i2.2))')year,month,day,hour,minute
      endif

      ! we dont know the final character on this file name. I'm going to leave it off 
      ! and let advance_model match against files in the FORCING/ directory.
      write(iunit,'(a)')trim(datestring)//'.LDASIN_DOMAIN'
   enddo

   call close_file(iunit)

end if

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

if (get_debug_level() > 0) then
   call print_date( model_time,'dart_to_wrfHydro:wrfHydro model date')
   call print_time( model_time,'dart_to_wrfHydro:DART model time')
   call print_date( model_time,'dart_to_wrfHydro:wrfHydro model date',logfileunit)
   call print_time( model_time,'dart_to_wrfHydro:DART model time',logfileunit)

   if ( advance_time_present ) then
      call print_time(adv_to_time,'dart_to_wrfHydro:advance_to time')
      call print_date(adv_to_time,'dart_to_wrfHydro:advance_to date')
      call print_time(adv_to_time,'dart_to_wrfHydro:advance_to time',logfileunit)
      call print_date(adv_to_time,'dart_to_wrfHydro:advance_to date',logfileunit)
   endif
endif

call finalize_utilities()

end program dart_to_wrfHydro

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
