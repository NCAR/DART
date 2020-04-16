! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program MPD_text_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   MPD/MPD_text_to_obs - a program that only needs minor customization to read
!      in a text-based dataset - either white-space separated values or
!      fixed-width column data.
!
!     created 20 Mar 2020  Michael Ying
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, PI, DEG2RAD

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              open_file, close_file, error_handler, &
                              E_ERR, E_WARN, E_MSG

use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              operator(>=), get_time, &
                              operator(-), GREGORIAN, operator(+), print_date

use      location_mod, only : VERTISHEIGHT

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq

use      obs_kind_mod, only : get_num_quantities, MPD_ABSOLUTE_HUMIDITY

implicit none

character(len=*), parameter :: source = 'MPD_text_to_obs.f90'
character(len=*), parameter :: text_input_file = 'text.txt'
character(len=*), parameter :: obs_out_file    = 'obs_seq.out'

logical, parameter :: debug = .false.  ! set to .true. to print info

character(len=129) :: input_line
character(len=512) :: string1

integer :: oday, osec, rcio, iunit, ilinecount
integer :: year, month, day, hour, minute, second
integer :: num_copies, num_qc, max_obs

logical  :: first_obs

real(r8) :: abs_humid, terr, qc
real(r8) :: lat, lon, vert

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

! start of executable code

call initialize_utilities(source)

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

! This is sort of a do-nothing call to initialize the obs_kind_mod module.
! The only real purpose of this call is to make obs_kind_mod print its
! initialization report at the beginning of the output from MPD_text_to_obs.
! Without this call, the initialization report actually comes AFTER all
! the observations have been processed and right before the output file
! is closed.  We do not actually need to know what is in string1.
write(string1,*)' num_observations understood is ',get_num_quantities()

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

ilinecount = 1

obsloop: do    ! no end limit - have the loop break when input ends

   ! read in a line from the text file.   What you need to create an obs:
   !  location: lat, lon, and height in pressure or meters
   !  time: when the observation was taken
   !  type: from the DART list of obs types
   !  error: very important - the instrument error plus representativeness error
   !        (see html file for more info)

   ! read in entire text line into a buffer
   read(iunit, "(A)", iostat=rcio) input_line
   if (rcio == 0) then
      continue ! line read normally
   elseif (rcio < 0) then ! Normal end-of-file
      exit obsloop
   else
      write(string1,*)'got bad read code from input file, rcio = ', rcio, &
                      ', line number =',ilinecount
      call error_handler(E_ERR,source,string1,text2='Hard Stop. No output.')
   endif

   ! the input text file has its own hardcoded convention for abs humidity
   ! 'lat'  is degrees Latitude [-90,90]
   ! 'lon'  is degrees East Longitude [0,360]
   ! 'vert' is height in meters
   ! 'terr' is the observation error STANDARD DEVIATION, this is converted into
   !        units of variance inside DART. This should include both the instrument
   !        error and an estimate of the representativeness error.  

   read(input_line, *, iostat=rcio) lat, lon, vert, &
                              year, month, day, hour, minute, second, &
                              abs_humid, terr
   if (rcio /= 0) then
      write(string1,*)'Unable to parse line ',ilinecount,', rcio = ', rcio
      call error_handler(E_ERR,source,string1)
   endif

   if (debug) print *, 'observation ',ilinecount,' located at lat, lon = ', lat, lon

   ! skip any observations outside these bounds
   if ( lat >  90.0_r8 .or. lat <  -90.0_r8 ) cycle obsloop
   if ( lon <   0.0_r8 .or. lon >  360.0_r8 ) cycle obsloop

   ! put date into a dart time format
   time_obs = set_date(year, month, day, hour, minute, second)

   if (debug) call print_date(time_obs, 'next obs time is')

   ! extract time of observation into gregorian day, sec.
   call get_time(time_obs, osec, oday)

   ! make an obs derived type, and then add it to the sequence
   call create_3d_obs(lat, lon, vert, VERTISHEIGHT, abs_humid, &
                      MPD_ABSOLUTE_HUMIDITY, terr, oday, osec, qc, obs)
   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

   if (debug) print *, 'added abs humidity obs to output seq'

   ilinecount = ilinecount + 1

enddo obsloop

call close_file(iunit)

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   write(string1,*)'writing "'//trim(obs_out_file)//'", obs_count = ', get_num_obs(obs_seq)
   call error_handler(E_MSG,source,string1)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

end program MPD_text_to_obs
