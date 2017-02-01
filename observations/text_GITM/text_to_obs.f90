! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program text_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   text_to_obs - a program that only needs minor customization to read
!      in a text-based dataset - either white-space separated values or
!      fixed-width column data.
!
!     created 29 Mar 2010   nancy collins NCAR/IMAGe
!     modified 15 Aug 2012 Alexey Morozov (Univ. of Michigan)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, PI, DEG2RAD

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              open_file, close_file

use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              operator(>=), increment_time, get_time, &
                              operator(-), GREGORIAN, operator(+), print_date

use      location_mod, only : VERTISHEIGHT

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq

use      obs_kind_mod, only : SAT_TEMPERATURE, SAT_RHO

implicit none

character(len=64), parameter :: text_input_file = 'text.txt'
character(len=64), parameter :: obs_out_file    = 'obs_seq.out'

logical, parameter :: debug = .false.  ! set to .true. to print info

character (len=129) :: input_line

integer :: oday, osec, rcio, iunit, otype
integer :: year, month, day, hour, minute, second
integer :: num_copies, num_qc, max_obs

logical  :: first_obs

real(r8) :: temp, terr, qc
real(r8) :: lat, lon, vert

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

! start of executable code

call initialize_utilities('text_to_obs')

! time setup
call set_calendar_type(GREGORIAN)

! each observation in this series will have a single observation value
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
max_obs    = 100000
num_copies = 1
num_qc     = 1

! Set the DART data quality control.   0 is good data.
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8

! call the initialization code, and initialize two empty observation types
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
first_obs = .true.

! create a new, empty obs_seq file.  you must give a max limit
! on number of obs.  increase the size if too small.
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

! the first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call set_qc_meta_data(obs_seq, 1, 'Data QC')

! if you want to append to existing files (e.g. you have a lot of
! small text files you want to combine), you can do it this way,
! or you can use the obs_sequence_tool to merge a list of files
! once they are in DART obs_seq format.

!  ! existing file found, append to it
!  inquire(file=obs_out_file, exist=file_exist)
!  if ( file_exist ) then
!     call read_obs_seq(obs_out_file, 0, 0, max_obs, obs_seq)
!  endif

! open input text file
iunit = open_file(text_input_file, 'formatted', 'read')
if (debug) print *, 'opened input file ' // trim(text_input_file)

obsloop: do    ! no end limit - have the loop break when input ends

   ! read in a line from the text file.   What you need to create an obs:
   !  location: lat, lon, and height in pressure or meters
   !  time: when the observation was taken
   !  type: from the DART list of obs types
   !  error: very important - the instrument error plus representativeness error
   !        (see html file for more info)

   ! assume here a line is a type, location, time, value, obs error

   ! read in entire text line into a buffer
   read(iunit, "(A)", iostat=rcio) input_line
   if (rcio /= 0) then
      if (debug) print *, 'got bad read code from input file, rcio = ', rcio
      exit obsloop
   endif

   ! pull off the first value as an integer, to decode the type
   read(input_line, *, iostat=rcio) otype
   if (rcio /= 0) then
      if (debug) print *, 'got bad read code trying to get obs type, rcio = ', rcio
      exit obsloop
   endif

   if (debug) print *, 'next observation type = ', otype

   ! for this example, assume there is an obs type, where otype=1 is
   ! a temperature, and if otype=37, there's Rho (neutral density)
   ! The input text file has these as their own hardcoded convention.

   if (otype == 1) then !Temperature
      read(input_line, *, iostat=rcio) otype, lat, lon, vert, &
                                 year, month, day, hour, minute, second, &
                                 temp, terr
      if (rcio /= 0) then
         if (debug) print *, 'got bad read code getting rest of temp obs, rcio = ', rcio
         exit obsloop
      endif

   elseif (otype == 37) then !Rho
      read(input_line, *, iostat=rcio) otype, lat, lon, vert, &
                                 year, month, day, hour, minute, second, &
                                 temp, terr
      if (rcio /= 0) then
         if (debug) print *, 'got bad read code getting rest of temp obs, rcio = ', rcio
         exit obsloop
      endif

   else !WHAT?

      ! no method defined to convert this type
      write(*,*) "NOOOOO, DON'T DO THIS: no method defined to convert this type = ", otype

   endif

   if (debug) print *, 'next observation located at lat, lon = ', lat, lon

   ! check the lat/lon values to see if they are ok
   if ( lat >  90.0_r8 .or. lat <  -90.0_r8 ) cycle obsloop
   if ( lon <   0.0_r8 .or. lon >  360.0_r8 ) cycle obsloop

   ! put date into a dart time format
   time_obs = set_date(year, month, day, hour, minute, second)

   if (debug) call print_date(time_obs, 'next obs time is')

   ! extract time of observation into gregorian day, sec.
   call get_time(time_obs, osec, oday)

   ! this example assumes there is an obs type, where otype=1 is
   ! a temperature measured in height, and if otype=38, there's a TEC obs
   if (otype == 1) then !Temperature

      ! height is in meters (gph)
      ! make an obs derived type, and then add it to the sequence
      call create_3d_obs(lat, lon, vert, VERTISHEIGHT, temp, &
                         SAT_TEMPERATURE, terr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
      if (debug) print *, 'added Temperature obs to output seq'

   elseif (otype == 37) then !SAT_RHO=37, but probably only in GITM
      ! height is in meters (gph)
      call create_3d_obs(lat, lon, vert, VERTISHEIGHT, temp, &
                         SAT_RHO, terr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
      if (debug) print *, 'added Rho obs to output seq'


   else !WHAT???

      ! no method defined to convert this type
      write(*,*) 'no method defined to convert this type = ', otype

   endif

enddo obsloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

end program text_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
