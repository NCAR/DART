! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program CHAMP_density_text_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   CHAMP_density_text_to_obs - a program that only needs minor customization to read
!      in a text-based dataset - either white-space separated values or
!      fixed-width column data.
!
!     created 29 Mar 2010   nancy collins NCAR/IMAGe
!
!+ modified 15 Aug 2012 Alexey Morozov (Univ. of Michigan), alexeymor at google mail
!
!+ It is designed to read CHAMP ascii files
!+ For example of input files, see Density_3deg_02_335.ascii in work folder, which is taken from
!  http://sisko.colorado.edu/sutton/data/ver2.2/champ/density/2002/ascii/
!+ This program reads the name of the text file, obs_seq file, and "debug" from input.nml
!+ APPENDS new observations to existing obs_seq.out - see line 130ish
!+ (but not if you change the obs_out_file in input.nml)
!+ For added convenience, see convert.sh in work folder, which runs this program repeatedly
!+ to convert+append many CHAMP files
!+ Implemented the suggestion about times starting from weird points (like 335th day in 2002)
!+ - see lines 190ish
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, PI, DEG2RAD

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              open_file, close_file, find_namelist_in_file, check_namelist_read

use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              operator(>=), increment_time, get_time, set_time, &
                              operator(-), GREGORIAN, operator(+), print_date

use      location_mod, only : VERTISHEIGHT

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq

use      obs_kind_mod, only : SAT_RHO !this is density

implicit none

! things which can/should be in the text_to_obs_nml

character(len=64)  :: text_input_file = 'Density_3deg_02_335_2p2.ascii'
character(len=64)  :: obs_out_file    = 'obs_seq.out'
logical            :: debug = .true.

namelist /CHAMP_density_text_to_obs_nml/  &
     text_input_file, &
     obs_out_file,    &
     debug


character (len=200) :: input_line !162 is the width of CHAMPdens2.2 files, but to be safe do 200

integer :: oday, osec, rcio, iunit
integer :: year, day, second
integer :: num_copies, num_qc, max_obs

logical  :: file_exist, first_obs

real(r8) :: temp, terr, qc
real(r8) :: lat, lon, vert

real(r8) :: second_r !CHAMP seconds are reals instead of ints

!variables to be discarded (only needed so that the read line works)
integer  :: ignore_i
real     :: ignore_r

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

! start of executable code

call initialize_utilities('CHAMP_density_text_to_obs')

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'CHAMP_density_text_to_obs_nml', iunit)
read(iunit, nml = CHAMP_density_text_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, 'CHAMP_density_text_to_obs_nml')

! time setup
call set_calendar_type(GREGORIAN)


! open input text file

iunit = open_file(text_input_file, 'formatted', 'read')
if (debug) print *, 'opened input file ' // trim(text_input_file)


! each observation in this series will have a single observation value
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
max_obs    = 100000
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
call set_qc_meta_data(obs_seq, 1, 'Data QC')

! if you want to append to existing files (e.g. you have a lot of
! small text files you want to combine), you can do it this way,
! or you can use the obs_sequence_tool to merge a list of files
! once they are in DART obs_seq format.

! existing file found, append to it
inquire(file=obs_out_file, exist=file_exist)
if ( file_exist ) then
  call read_obs_seq(obs_out_file, 0, 0, max_obs, obs_seq)
endif

! Set the DART data quality control.   0 is good data.
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8

! first two lines are just text (description), so just skip them
read(iunit, "(A)", iostat=rcio) input_line
read(iunit, "(A)", iostat=rcio) input_line

obsloop: do    ! no end limit - have the loop break when input ends

   ! read in a line from the text file.   What you need to create an obs:
   !  location: lat, lon, and height in pressure or meters
   !  time: when the observation was taken
   !  type: from the DART list of obs types
   !  error: very important - the instrument error plus representativeness error
   !        (see html file for more info)

   ! read the whole line into a buffer
   read(iunit, "(A)", iostat=rcio) input_line
   if (rcio /= 0) then
      if (debug) print *, 'got bad read code from input file, rcio = ', rcio
      exit obsloop
   endif

   ! assume here is a line from sisko.colorado.edu/sutton/data/ver2.2/champ/density/2002/ascii/,
   !data format is:
   !+ 1)year(2I), 2)day(3I), 3)second(8.3F), 4)round(lat), 5)lat(d,-90 90), 6)lon(d,-180 180), 7)alt(km),
   !+ 8)LT, 9)Mlat, 10)Mlon, 11)MLT, 12)Rho(Density!), 13)MSISRho400, 14)MSISRho410, 15)MSISRhoSat
   !+ 16)Rho uncertainty (I'm guessing std deviation from units: kg/m^3) 17)points averaged over
   !+ 18)points needing interpolation 19)coeff of drag averaged over bin

   read(input_line, *, iostat=rcio) &
        year, day, second_r, ignore_i, lat, lon, vert, &
        ignore_r, ignore_r, ignore_r, ignore_r, temp, ignore_r, ignore_r, ignore_r, &
        terr, ignore_i, &
        ignore_i, ignore_r

   vert=vert*1000 !DART needs alt in m, whereas in champ files it's in km

   if (rcio /= 0) then
      if (debug) print *, 'got bad read code getting rest of temp obs, rcio = ', rcio
      exit obsloop
   endif

   if (debug) print *, 'this observation located at lat, lon = ', lat, lon

   ! if lon comes in between -180 and 180, use these lines instead:
   if ( lat >  90.0_r8 .or. lat <  -90.0_r8 ) cycle obsloop
   if ( lon > 180.0_r8 .or. lon < -180.0_r8 ) cycle obsloop
   if ( lon < 0.0_r8 )  lon = lon + 360.0_r8 ! changes into 0-360

   ! put date into a dart time format

   year = 2000 + year !because year in file is (2I) - 2 digits
   second = nint(second_r)

   !! some times are supplied as number of seconds since some reference
   !! date.  This is an example of how to support that.
   !! put the reference date into DART format
   comp_day0 = set_date(year, 1, 1, 0, 0, 0)
   time_obs  = comp_day0 + set_time(second, day-1)

   ! extract time of observation into gregorian day, sec.
   call get_time(time_obs, osec, oday)

   if (debug) call print_date(time_obs, 'this obs time is')

   ! height is in kilometers (yardstick)
   ! make an obs derived type, and then add it to the sequence

   call create_3d_obs(lat, lon, vert, VERTISHEIGHT, temp, &
        SAT_RHO, terr, oday, osec, qc, obs)
   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

   if (debug) print *, 'added RHO obs to output seq'

end do obsloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   !if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

end program CHAMP_density_text_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
