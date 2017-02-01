! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program tc_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   tc_to_obs - read in ascii data about tropical cyclones and create
!      DART observations which correspond to that information.
!
!     created 21 Jul 2014   nancy collins NCAR/IMAGe
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8
use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              open_file, close_file, find_namelist_in_file, &
                              check_namelist_read, nmlfileunit, &
                              do_nml_file, do_nml_term, &
                              get_next_filename, error_handler, E_ERR, E_MSG
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              operator(>=), increment_time, get_time, &
                              operator(-), GREGORIAN, operator(+), print_date
use      location_mod, only : VERTISUNDEF
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use      obs_kind_mod, only : VORTEX_LAT, VORTEX_LON, VORTEX_PMIN, &
                              VORTEX_WMAX
use obs_utilities_mod, only : add_obs_to_seq, create_3d_obs

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


! IMPORTANT!  These are the observational errors for position (in degrees),
! minimum sea level pressure (in hectopascals/mb), and max wind (in m/s).
! They are hardcoded here; to change their values either edit them here
! or add code to make them dependent on the magnitude of the obs values.

real(r8), parameter :: tc_pos_error  = 0.30_r8,  &
                       tc_pmin_error = 3.00_r8,  &
                       tc_wmax_error = 3.00_r8


logical  :: file_exist, first_obs
real(r8) :: lat, lon, qc
real(r8) :: mslp, mwnd, wsval
integer  :: fhr, istat, io
integer  :: iyear, imonth, iday, ihour
integer  :: oday, osec, rcio, iunit
integer  :: num_copies, num_qc, max_obs
           
character(len=256) :: input_line, msgstring
character(len=10)  :: yyyymmddhh
character(len=4)   :: ftype
character(len=2)   :: basin, tcnum, stypei
character(len=1)   :: ew, ns

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

! things in the namelist
character(len=256) :: input_atcf_file         = 'input.txt'
character(len=128) :: fileformat              = 'b-deck'
character(len=256) :: obs_out_file            = 'obs_seq.out'
logical            :: append_to_existing_file = .false.
logical            :: debug                   = .false.  ! set to .true. to print more info


namelist /tc_to_obs_nml/ &
   input_atcf_file, fileformat, obs_out_file, append_to_existing_file, debug


! start of executable code

call initialize_utilities('tc_to_obs')

! Read the namelist entry
call find_namelist_in_file("input.nml", "tc_to_obs_nml", iunit)
read(iunit, nml = tc_to_obs_nml, iostat = io)
call check_namelist_read(iunit, io, "tc_to_obs_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=tc_to_obs_nml)
if (do_nml_term()) write(     *     , nml=tc_to_obs_nml)

! time setup
call set_calendar_type(GREGORIAN)

! open input text file
iunit = open_file(input_atcf_file, 'formatted', 'read')
if (debug) print *, 'opened input ATCF file ' // trim(input_atcf_file)

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

! if you want to append to existing files (e.g. you have a lot of
! small text files you want to combine), you can do it here,
! or you can use the obs_sequence_tool to merge a list of files 
! once they are in DART obs_seq format.

! existing file found, append to it if user requests this via the namelist
inquire(file=obs_out_file, exist=file_exist)

if ( file_exist .and. append_to_existing_file) then
   if (debug) print *, 'found existing obs_seq file, appending to ' // trim(obs_out_file)
   call read_obs_seq(obs_out_file, 0, 0, max_obs, obs_seq)
else
   ! create a new, empty obs_seq file.  you must give a max limit
   ! on number of obs.  increase the size if too small.
   if (debug) print *, 'creating new obs_seq file for output ' // trim(obs_out_file)
   call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)
   
   ! the first one needs to contain the string 'observation' and the
   ! second needs the string 'QC'.
   call set_copy_meta_data(obs_seq, 1, 'observation')
   call set_qc_meta_data(obs_seq, 1, 'Data QC')
endif

! Set the DART data quality control.   0 is good data. 
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8

OBSLOOP: do    ! no end limit - have the loop break when input ends

   if (debug) print *, ''

   ! read in a line from the text file.   What you need to create an obs:
   !  location: lat, lon, and for these types of obs the vertical is 'UNDEF'
   !  time: when the observation was taken
   !  type: from the DART list of obs types, typically VORTEX_xxx
   !  error: very important - the instrument error plus model representativeness error
   !        (see html doc file for more info)

   ! read the entire text line into a buffer here.  this makes it easier to
   ! parse lines which may have variable numbers of fields.  the later read()
   ! then reads out of this character variable.
   read(iunit, "(A)", iostat=rcio) input_line
   if (rcio /= 0) then 
      if (debug) print *, 'read code from input file, rcio = ', rcio
      if (debug) print *, '(normal if we have reached the end of file)'
      exit OBSLOOP
   endif

   if (debug) print *, 'next input line: ' // trim(input_line)  

   select case(fileformat)
    case ("b-deck")
       ! read the various data fields out of the most recently read input line
       read(input_line, 21, iostat=istat ) basin, tcnum, yyyymmddhh, ftype, fhr, &
                                           lat, ns, lon, ew, mwnd, mslp, stypei, wsval

21     format(a2,1x,a2,1x,a10,1x,a4,1x,i1,1x, &
              f5.2,a1,1x,f6.2,a1,1x,f5.1,1x,f6.1,1x,a2)

       read(yyyymmddhh(1:4),  fmt='(i4)') iyear
       read(yyyymmddhh(5:6),  fmt='(i2)') imonth
       read(yyyymmddhh(7:8),  fmt='(i2)') iday
       read(yyyymmddhh(9:10), fmt='(i2)') ihour

       if ( ns == 'S' )  lat = -1.0_r8 * lat
       if ( ew == 'W' )  lon = -1.0_r8 * lon

       mslp = mslp * 100.0_r8   ! pressures in dart are in pascals

       time_obs = set_date(iyear, imonth, iday, ihour, 0, 0)
       call get_time(time_obs, osec, oday)
   
    ! case ("a-deck")
    ! case ("other")
    !  add different read format lines for other data layouts here

    case default
      write(msgstring, *) 'unrecognized fileformat string: "'//trim(fileformat)//'"'
      call error_handler(E_ERR,'tc_to_obs', msgstring, source, revision, revdate)
   end select


   ! check the lat/lon values to see if they are ok
   if ( lat >  90.0_r8 .or. lat <  -90.0_r8 ) cycle OBSLOOP
   if ( lon > 360.0_r8 .or. lon < -180.0_r8 ) cycle OBSLOOP
   if ( lon < 0.0_r8 )  lon = lon + 360.0_r8 ! dart needs lon 0-360

   if (debug) call print_date(time_obs, ' time of observations: ')

   ! make an obs derived type, and then add it to the sequence
   if (debug) print *, ' adding vortex  latitude observation, value,err (in degrees) = ', lat, tc_pos_error
   call create_3d_obs(lat, lon, 0.0_r8, VERTISUNDEF, lat, &
                      VORTEX_LAT, tc_pos_error, oday, osec, qc, obs)
   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
   ! make an obs derived type, and then add it to the sequence
   if (debug) print *, ' adding vortex longitude observation, value,err (in degrees) = ', lon, tc_pos_error
   call create_3d_obs(lat, lon, 0.0_r8, VERTISUNDEF, lon, &
                      VORTEX_LON, tc_pos_error, oday, osec, qc, obs)
   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

   ! make an obs derived type, and then add it to the sequence
   if (debug) print *, ' adding          min slp observation, value,err (in pascals) = ', mslp, tc_pmin_error*100.0_r8
   call create_3d_obs(lat, lon, 0.0_r8, VERTISUNDEF, mslp, &
                      VORTEX_PMIN, tc_pmin_error*100.0_r8, oday, osec, qc, obs)
   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

   if (mwnd > 0.0_r8) then
      ! make an obs derived type, and then add it to the sequence
      if (debug) print *, ' adding         max wind observation, value,err (in m/sec)   = ', mwnd, tc_wmax_error
      call create_3d_obs(lat, lon, 0.0_r8, VERTISUNDEF, mwnd, &
                         VORTEX_WMAX, tc_wmax_error, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   endif

end do OBSLOOP

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (debug) print *, 'writing output obs_seq file, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()


end program tc_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
