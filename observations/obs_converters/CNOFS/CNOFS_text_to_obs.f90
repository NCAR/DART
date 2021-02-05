! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program CNOFS_text_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   CNOFS_text_to_obs - a program that only needs minor customization to read
!      in a text-based dataset - either white-space separated values or
!      fixed-width column data.
!
!     created 29 Mar 2010   nancy collins NCAR/IMAGe
!     modified 15 Aug 2012 Alexey Morozov (Univ. of Michigan)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, PI, DEG2RAD

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              open_file, close_file, find_namelist_in_file, check_namelist_read

use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              operator(>=), increment_time, get_time, &
                              operator(-), GREGORIAN, operator(+), print_date

use      location_mod, only : VERTISHEIGHT

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq

use      obs_kind_mod, only : SAT_DENSITY_ION_E, SAT_TEMPERATURE_ION


implicit none

!alex: things which can/should be in the CNOFS_text_to_obs_nml

character(len=64)  :: text_input_file = 'Density_3deg_02_335_2p2.ascii'
character(len=64)  :: obs_out_file    = 'obs_seq.out'
logical            :: debug = .true.  ! set to .true. to print info

namelist /CNOFS_text_to_obs_nml/  &
     text_input_file, &
     obs_out_file,    &
     debug

character (len=200) :: input_line !alex: GPS TEC data is 153 wide, but give myself some room

integer :: oday, osec, rcio, iunit
integer :: year, month, day, hour, minute, second
integer :: num_copies, num_qc, max_obs

logical  :: file_exist, first_obs

real(r8) :: temp, terr, qc, dens, denserr  !alexey
real(r8) :: lat, lon, vert

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

integer  :: ignore_i        ! integers we read but don't use
real(r8) :: ignore_r, mi,ma ! calculate the min and max of the dataset on the fly

! start of executable code

call initialize_utilities('CNOFS_text_to_obs')

! Read the DART namelist for this model !alex
call find_namelist_in_file('input.nml', 'CNOFS_text_to_obs_nml', iunit)
read(iunit, nml = CNOFS_text_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, 'CNOFS_text_to_obs_nml')

! time setup
call set_calendar_type(GREGORIAN)

! open input text file
iunit = open_file(text_input_file, 'formatted', 'read')
if (debug) print *, 'opened input file ' // trim(text_input_file)

! each observation in this series will have a single observation value
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
max_obs    = 10000000
num_copies = 1
num_qc     = 1

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
call   set_qc_meta_data(obs_seq, 1, 'Data QC')

! if you want to append to existing files (e.g. you have a lot of
! small text files you want to combine), you can do it this way,
! or you can use the obs_sequence_tool to merge a list of files
! once they are in DART obs_seq format.

!  ! existing file found, append to it
!  inquire(file=obs_out_file, exist=file_exist)
!  if ( file_exist ) then
!     call read_obs_seq(obs_out_file, 0, 0, max_obs, obs_seq)
!  endif

! Set the DART data quality control.   0 is good data.
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8

read(iunit, "(A)", iostat=rcio) input_line ! skip the first line as it is a header
read(iunit, "(A)", iostat=rcio) input_line ! skip the second line as it is #START

mi =  1.0e32 !set initial minimum as large positive
ma = -1.0e32 !set initial maximum as large negative

obsloop: do    ! no end limit - have the loop break when input ends

   ! read in a line from the text file.   What you need to create an obs:
   !  location: lat, lon, and height in pressure or meters
   !  time: when the observation was taken
   !  type: from the DART list of obs types
   !  error: very important - the instrument error plus representativeness error
   !        (see html file for more info)

   ! read the whole line
   read(iunit, "(A)", iostat=rcio) input_line

   if (rcio /= 0) then
      if (debug) print *, 'got bad read code from input file, rcio = ', rcio
      exit obsloop
   endif

   !here is the header
   !year mm dd hh mm ss msec long lat alt vpar vpar_var idens idens_var itemp itemp_var
   !here is a line from Angeline's cnofs_jan_2012_data_assim.ivm
   !2012 1 001 1 38 10 291 331.65 13.0664 407.954 -44.3753 1156.34 638705 375394 825.747 6078.96

   read(input_line, *, iostat=rcio) year, month, day, hour, minute, second, &
        ignore_i, lon, lat, vert, ignore_r, ignore_r, dens, denserr, temp, terr

   vert    = vert*1000.0_r8     ! convert km to m
   dens    = dens * 10.0**6     ! convert ion density cm^-3 to m^-3
   denserr = denserr * 10.0**12 ! convert ion density variance cm^-6 to m^-6
   denserr = sqrt(denserr)      ! convert var to std dev (m^-3)
                                ! temp is ion temp in K
   terr    = sqrt(terr)         ! convert var to std dev (K)

   if (rcio /= 0) then
      if (debug) print *, 'got bad read code getting rest of temp obs, rcio = ', rcio
      exit obsloop
   endif

   if (debug) print *, 'next observation located at lat, lon = ', lat, lon

   ! check the lat/lon values to see if they are ok
   if ( lat >  90.0_r8 .or. lat <  -90.0_r8 ) cycle obsloop
   if ( lon <   0.0_r8 .or. lon >  360.0_r8 ) cycle obsloop

   mi = min(mi, temp)
   ma = max(ma, temp)

   ! put date into a dart time format
   time_obs = set_date(year, month, day, hour, minute, second)

   if (debug) call print_date(time_obs, 'next obs time is')

   ! extract time of observation into gregorian day, sec.
   call get_time(time_obs, osec, oday)

   call create_3d_obs(lat, lon, vert, VERTISHEIGHT, dens, &
        SAT_DENSITY_ION_E, denserr, oday, osec, qc, obs)
   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   if (debug) print *, 'added an obs to output seq, now the total count is', get_num_obs(obs_seq)

   call create_3d_obs(lat, lon, vert, VERTISHEIGHT, temp, &
        SAT_TEMPERATURE_ION, terr, oday, osec, qc, obs)
   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   if (debug) print *, 'added an obs to output seq, now the total count is', get_num_obs(obs_seq)

enddo obsloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   print *, 'writing obs_seq, obs_count, min, max = ', get_num_obs(obs_seq), mi, ma
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

end program CNOFS_text_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
