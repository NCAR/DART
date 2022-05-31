! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program gnd_gps_vtec_text_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   gnd_gps_vtec_text_to_obs - a program that only needs minor customization to read
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

use      obs_kind_mod, only : GND_GPS_VTEC

implicit none

! things which can/should be in the gnd_gps_vtec_text_to_obs_nml

character(len=64)  :: text_input_file = 'Density_3deg_02_335_2p2.ascii'
character(len=64)  :: obs_out_file    = 'obs_seq.out'
logical            :: debug = .true.  ! set to .true. to print info

namelist /gnd_gps_vtec_text_to_obs_nml/  &
     text_input_file, &
     obs_out_file,    &
     debug

character (len=200) :: input_line ! GPS TEC data is 153 wide, but give myself some room

integer :: oday, osec, rcio, iunit
integer :: year, month, day, hour, minute, second
integer :: num_copies, num_qc, max_obs

logical  :: file_exist, first_obs

real(r8) :: temp, terr, qc
real(r8) :: lat, lon, vert

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

integer  :: ignore_i ! integers we read but don't use
real(r8) :: mi,ma    ! calculate the min and max of the dataset on the fly

! start of executable code

call initialize_utilities('gnd_gps_vtec_text_to_obs')

write(*,*)
write(*,*) "WARNING, BE CAREFUL: the TEC files from different dates might have different number of columns!"
write(*,*) "For example, on some dates in 2002 there are 13 columns, whereas in 2006 there are 10 "
write(*,*) "If you are getting some weird data/results, this could be the issue "
write(*,*) "To fix it, change ../gnd_gps_vtec_text_to_obs.f90:line~175 and recompile with ./quickbuild.sh"
write(*,*)

! Set the DART data quality control.   0 is good data.
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8
mi = 1.0e32 ! set initial minimum as large positive
ma =-1.0e32 ! set initial maximum as large negative

! each observation in this series will have a single observation value
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
max_obs    = 10000000
num_copies = 1
num_qc     = 1

! time setup
call set_calendar_type(GREGORIAN)

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'gnd_gps_vtec_text_to_obs_nml', iunit)
read(iunit, nml = gnd_gps_vtec_text_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, 'gnd_gps_vtec_text_to_obs_nml')

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

read(iunit, "(A)", iostat=rcio) input_line ! skip the first line as it is a header

obsloop: do    ! FIXME no end limit - have the loop break when input ends

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

  !BE CAREFUL: the files from different years have different number of columns!
  !for example, on some dates in 2002 it is 13, in 2006 it is 10

   ! assume here is a line from http://madrigal.haystack.mit.edu/cgi-bin/madrigal/displayCatHead.py?expName=experiments/
   !2002/gps/30nov02&fileName=gps021130g.002&expTitle=World-wide%20Vertical%20Total%20Electron%20Content&displayLevel=0
   !here is the header
   !YEAR MONTH DAY HOUR MIN SEC UT1_UNIX UT2_UNIX RECNO GDLAT GLON TEC DTEC
   !here is a sample line (white space removed by hand)
   !2002 12 1 0 2 30 1038700800 1038701100 2 -82.00 -178.00 8.40000e+00 5.00000e+00

   read(input_line, *, iostat=rcio) year, month, day, hour, minute, second, &
        ignore_i, ignore_i, ignore_i, lat, lon, temp, terr

   !temp is the TEC value
   !terr is the STANDARD DEVIATION (measured in TECU) uncertainty associated with that measurement

   ! replacing "harsh" exit with a "soft" cycle because GPS data has a lot of "missing" values
   ! which we simply skip (i.e. throw away)
   if (rcio /= 0) then
      if (debug) print *, 'got bad read code getting rest of temp obs, rcio = ', rcio
      cycle obsloop
   endif

   if (debug) print *, 'next observation located at lat, lon = ', lat, lon

   ! check the lat/lon values to see if they are ok
   if ( lat >= 90.0_r8 .or. lat <= -90.0_r8 ) cycle obsloop
   if ( lon > 180.0_r8 .or. lon < -180.0_r8 ) cycle obsloop
   if ( lon < 0.0_r8 )  lon = lon + 360.0_r8 ! changes into 0-360

   ! after all checks were done on this data point, update the min and max
   mi = min(mi, temp)
   ma = max(ma, temp)

   ! put date into a dart time format
   time_obs = set_date(year, month, day, hour, minute, second)

   if (debug) call print_date(time_obs, 'next obs time is')

   ! extract time of observation into gregorian day, sec.
   call get_time(time_obs, osec, oday)

   !height is tricky since ground TEC measurement is an integral over height, so it has no height. I'll set it to
   !350 000m, but in model_mod I'll make TEC measurement localized only in horizontal (no vertical).

   vert=350000.0
   call create_3d_obs(lat, lon, vert, VERTISHEIGHT, temp, &
                      GND_GPS_VTEC, terr, oday, osec, qc, obs)
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

end program gnd_gps_vtec_text_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
