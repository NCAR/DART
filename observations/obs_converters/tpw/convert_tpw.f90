! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_tpw

! convert MODIS observations of Total Precipitable Water into
! DART observation sequence files.  since the original observations
! are very dense, this program has options to bin both in time and
! in space.  THIS CODE ASSUMES ALL ELEVATIONS ARE AT 0.0 meters.
! this is true for all ocean obs.  if you have land obs, you will
! have to get the hght of the surface elevation at the lat/lon and
! use that instead.
!
! assumes data is in text files which match the following read line:
!      read(iunit, '(f11.6, f13.5, f10.4, 4x, i4, 4i3, f7.3)') &
!                lat, lon, tpw, iyear, imonth, iday, ihour, imin, seconds
!
! constructs a input filename based on:
!  ObsBase/InfilePrefix + YYYYMMDD + InfileSuffix
! 
! constructs a output filename based on:
!  ./OutfilePrefix + YYYYMMDD + OutfileSuffix
!
! any of the prefix or suffixes can be '' (blank)
! 
! FIXME:
! cannot loop past month boundaries.  also does not handle time windows
! that span day boundaries, in the sense that it reads one day at a time
! and outputs one obs_seq file for each day.  if a time bin crosses the
! day boundary it should read in the next day, keep binning the obs for
! the rest of the current time bin, and write out an obs_seq file that
! has a full bin of all available observations.  the current code starts
! at 0Z and ends at 0Z and so you should construct your time bins carefully
! around the day boundaries.


use         types_mod, only : r8, metadatalength, missing_r8
use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time,&
                              increment_time, get_time, set_date, operator(-),  &
                              print_date, decrement_time, operator(>), print_time
use     utilities_mod, only : initialize_utilities, find_namelist_in_file,      &
                              check_namelist_read, nmlfileunit, do_nml_file,    &
                              get_next_filename, error_handler, E_ERR, E_MSG,   &
                              find_textfile_dims, do_nml_term,                  &
                              is_longitude_between, finalize_utilities,         &
                              open_file, close_file, register_module
use      location_mod, only : VERTISSURFACE, set_location
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,        &
                              static_init_obs_sequence, init_obs, destroy_obs,  &
                              write_obs_seq, init_obs_sequence, get_num_obs,    &
                              insert_obs_in_seq, destroy_obs_sequence,          &
                              set_copy_meta_data, set_qc_meta_data, set_qc,     &
                              set_obs_values, set_obs_def, insert_obs_in_seq
use       obs_def_mod, only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                              set_obs_def_error_variance, set_obs_def_location, &
                              set_obs_def_key
use      obs_kind_mod, only :  AQUA_TOTAL_PRECIPITABLE_WATER,  &
                              TERRA_TOTAL_PRECIPITABLE_WATER,  &
                               AMSR_TOTAL_PRECIPITABLE_WATER,  &
                              MODIS_TOTAL_PRECIPITABLE_WATER,  &
                              get_index_for_type_of_obs
use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
 "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

!--------------------------------------


type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs 
type(time_type)         :: obs_time, prev_time
type(time_type)         :: time_diff

integer  :: ibin, nr, io, obstype
integer  :: iyear, imonth, iday, ihour, imin
integer  :: nrec, idd, iunit, num_obs, num_avg
integer  :: nbin, kk, nn, oday, osec, sec_diff

integer, allocatable :: this_bin(:)

integer, parameter :: ncol       = 9
integer, parameter :: num_copies = 1
integer, parameter :: num_qc     = 1

real(r8) :: oerr, qc
real(r8) :: bin_beg, bin_end
real(r8) :: lon, lat, tpw, seconds
real(r8) :: lonk, dlat, hght
real(r8) :: nloni, nlonk

logical  :: first_obs, lon_between

character(len = 8)   :: obsdate
character(len = 256) :: msgstring, infile, outfile

type tpw_type
  real(r8) :: lat
  real(r8) :: lon
  real(r8) :: tpw_value
  type(time_type) :: time
  real(r8) :: fract_hours
  logical  :: been_used
end type

type(tpw_type), allocatable :: tpw_all(:)
type(tpw_type)              :: tpw_base, tpw_next, tpw_avg

! items in namelist, along with default values
integer  :: start_year  = 2008
integer  :: start_month = 1
integer  :: start_day   = 1
integer  :: total_days  = 31
integer  :: max_obs     = 150000
real(r8) :: time_bin_start      =  0.00_r8  ! fractional hours
real(r8) :: time_bin_interval   =  0.50_r8  ! fractional hours
real(r8) :: time_bin_half_width =  0.25_r8  ! fractional hours
real(r8) :: time_bin_end        = 24.00_r8  ! fractional hours
real(r8) :: delta_lat_box = 1.0_r8
real(r8) :: delta_lon_box = 1.0_r8
real(r8) :: min_lon = missing_r8
real(r8) :: max_lon = missing_r8
real(r8) :: min_lat = missing_r8
real(r8) :: max_lat = missing_r8
! the date, in 'YYYYMMDD' format, will be inserted between
! the input and output file prefix and suffix.  ObsBase is
! only prepended to the input file
character(len=128) :: ObsBase       = '../data'
character(len=64)  :: InfilePrefix  = 'datafile.'
character(len=64)  :: InfileSuffix  = '.txt'
character(len=64)  :: OutfilePrefix = 'obs_seq.'
character(len=64)  :: OutfileSuffix = ''
character(len=32)  :: observation_name = 'MODIS_TOTAL_PRECIPITABLE_WATER'

namelist /convert_tpw_nml/ start_year, start_month, start_day, &
      total_days, max_obs, delta_lat_box, delta_lon_box,       &
      time_bin_start, time_bin_interval, time_bin_half_width,  &
      min_lon, max_lon, min_lat, max_lat, time_bin_end,        &
      ObsBase, InfilePrefix, InfileSuffix, OutfilePrefix,      &
      OutfileSuffix, observation_name


! ----------------------------------------------------------------------
! start of executable program code
! ----------------------------------------------------------------------

call initialize_utilities('convert_tpw')
call register_module(source,revision,revdate)

! Initialize the obs_sequence module 

call static_init_obs_sequence()

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file('input.nml', 'convert_tpw_nml', iunit)
read(iunit, nml = convert_tpw_nml, iostat = io)
call check_namelist_read(iunit, io, 'convert_tpw_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=convert_tpw_nml)
if (do_nml_term()) write(    *      , nml=convert_tpw_nml)

! some simple error checking - if you specify either of the longitude
! limits, since it's a cyclic value you have to specify both.  the included
! area will start at min_lon and go east until max_lon.   for latitude
! you can let one default if you only want to set a min or a max.
if (((min_lon /= missing_r8) .and. (max_lon == missing_r8)) .or. &
    ((min_lon == missing_r8) .and. (max_lon /= missing_r8))) then
   write(msgstring, *) 'if you set the min or max limit on longitude you must specify both'
   call error_handler(E_ERR,'convert_tpw', msgstring, source, revision, revdate)
endif

! convert a string into an observation type number
obstype = get_index_for_type_of_obs(observation_name)
if (obstype < 0) then
   write(msgstring, *) 'unrecognized observation type ', trim(observation_name)
   call error_handler(E_ERR,'convert_tpw', msgstring, source, revision, revdate, &
                      text2='check the obs_defs in your &preprocess list')
endif

call set_calendar_type(GREGORIAN)

! allocate arrays which depend on input values in the namelist
allocate(tpw_all(max_obs), this_bin(max_obs))

!--------------------------------------

DAYLOOP: do idd = start_day, start_day + total_days - 1

   ! set up for next obs_seq for next day

   call init_obs(obs,      num_copies, num_qc)
   call init_obs(prev_obs, num_copies, num_qc)
   first_obs = .true.
   num_obs = 0

   call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

   call set_copy_meta_data(obs_seq, 1, 'MODIS observation')
   call set_qc_meta_data(obs_seq, 1, 'MODIS QC')

   ! construct the input and output filenames
   write(obsdate, '(i4.4, i2.2, i2.2)') start_year, start_month, idd
   write(infile, '(A)') trim(ObsBase)//'/'//trim(InfilePrefix)//obsdate//trim(InfileSuffix)
   write(*, '(A)') 'reading  file: '//trim(infile)
   write(outfile, '(A)') trim(OutfilePrefix)//obsdate//trim(OutfileSuffix)
   write(*, '(A)') 'creating file: '//trim(outfile)

   iunit = open_file(infile,'formatted',action='read')

   nrec = 0

   ! read in all observation records and fill the tpw type array
   READLOOP: do 

      read(iunit, '(f11.6, f13.5, f10.4, 4x, i4, 4i3, f7.3)', iostat=io) &
         lat, lon, tpw, iyear, imonth, iday, ihour, imin, seconds
      if (io /= 0) exit READLOOP

      nrec = nrec + 1

      if(nrec > max_obs) then
         write(msgstring, *) 'number of input observations larger than "max_obs" limit'
         call error_handler(E_ERR,'convert_tpw', msgstring, source, revision, revdate, &
                            text2='increase parameter "max_obs" in namelist and rerun')
      endif

      tpw_all(nrec)%lat = lat
      tpw_all(nrec)%lon = lon
      tpw_all(nrec)%tpw_value = tpw
      tpw_all(nrec)%time = set_date(iyear, imonth, iday, ihour, imin, int(seconds))
      tpw_all(nrec)%fract_hours = ihour + real(imin)/60.0_r8 + seconds/3600.0_r8
      tpw_all(nrec)%been_used = .false.
     
   enddo READLOOP

   call close_file(iunit)

   write(msgstring, '(A,I12)') 'number of observations in input file is ', nrec
   call error_handler(E_MSG,'', msgstring, source, revision, revdate)

! ---------------------------------------------------
!   select observations for specific time interval
! ---------------------------------------------------

   ibin = 0

   TIMEBINLOOP:  do

      bin_beg = time_bin_start + ibin*time_bin_interval - time_bin_half_width
      bin_end = time_bin_start + ibin*time_bin_interval + time_bin_half_width
      if (bin_beg >= time_bin_end) exit TIMEBINLOOP

      ibin = ibin + 1

      ! make a list of which tpw index numbers are in this time bin

      nbin = 0
      NRECLOOP: do nr = 1, nrec

         if(tpw_all(nr)%fract_hours <= bin_beg .or. tpw_all(nr)%fract_hours > bin_end) cycle NRECLOOP

         nbin = nbin + 1
         this_bin(nbin) = nr

      enddo NRECLOOP

      write(msgstring, '(A,F6.3,A,F6.3,A,I12)') 'num obs in time bin ', &
            bin_beg, ' to ', bin_end, ' hrs is ', nbin
      call error_handler(E_MSG,'', msgstring, source, revision, revdate)

      ! if there are no observations in this bin, go straight to the next bin
      if (nbin == 0) cycle TIMEBINLOOP

      ! main loop over all observations in this time bin
      OBSLOOP2: do kk = 1, nbin 
  
         tpw_base = tpw_all(this_bin(kk))

         ! if we have already averaged this observation in with another
         ! then skip it.
         if (tpw_base%been_used) cycle OBSLOOP2

         lonk = tpw_base%lon

         ! this will always return true, but if the longitude box crosses the
         ! prime meridian it will add 360 to all values on the east of the line. 
         ! when we compute an average the points will have consistent values. 
         ! (e.g. we will never average 355 and 5 and get 180 instead of 0.)
         lon_between = is_longitude_between(lonk, lonk-delta_lon_box, &
                                            lonk+delta_lon_box, newlon=nlonk)

         tpw_base%lon = nlonk

         ! use this observation as the base for averaging any other
         ! close obs nearby
         num_avg = 1
         tpw_avg = tpw_base

         ! average time by tracking the difference, in seconds, between the
         ! times of all the obs being averaged together.
         sec_diff = 0

         ! loop over the remaining observations
         ! do the average of the nearby observations
         OTHEROBS: do nn = (kk+1), nbin

            tpw_next = tpw_all(this_bin(nn))

            ! if we have already averaged this observation in with another
            ! then skip it.
            if (tpw_next%been_used) cycle OTHEROBS

            ! for latitudes you can do a simple subtraction.  since longitudes
            ! are cyclic the computation is more complex.  nloni will be returned
            ! as a value between 0 and 720 to guarentee that values in the delta_lon
            ! region are contiguous (no -180/180 or 359/0 discontinuous vals).
            dlat = abs(tpw_base%lat - tpw_next%lat)
            if (dlat > delta_lat_box) cycle OTHEROBS

            lon_between = is_longitude_between(tpw_next%lon, lonk-delta_lon_box, &
                                               lonk+delta_lon_box, newlon=nloni)
            if (.not. lon_between) cycle OTHEROBS

            ! obs is close enough to average

            num_avg = num_avg + 1
            tpw_avg%lat = tpw_avg%lat + tpw_next%lat
            tpw_avg%lon = tpw_avg%lon + nloni  ! possibly + 360 if near prime meridian
            tpw_avg%tpw_value = tpw_avg%tpw_value + tpw_next%tpw_value

            ! do the average by tracking the difference in seconds between the
            ! original time and all the times being averaged in.  this is to
            ! avoid adding a lot of large numbers of seconds together and running
            ! into roundoff error.  at the end, divide the total difference in secs
            ! by the count and add (or subtract) that from the initial time.

            ! when subtracting two time types, the time mgr always returns 
            ! positive vals, so we have to test to see if we need to add or
            ! subtract from the running total difference.
            time_diff = tpw_next%time - tpw_base%time   
            call get_time(time_diff, osec)

            if (tpw_next%time > tpw_base%time) then
               sec_diff = sec_diff + osec
            else
               sec_diff = sec_diff - osec
            endif

            ! since we have averaged this obs in with the others, remove
            ! it from the original list so it isn't used again.
            tpw_all(this_bin(nn))%been_used = .true.

         enddo OTHEROBS

         ! compute the average values, including location and time

         if (num_avg > 1) then
            lat = tpw_avg%lat / real(num_avg)
            lon = tpw_avg%lon / real(num_avg)
            tpw = tpw_avg%tpw_value / real(num_avg)

            ! sec_diff is the total running difference, in seconds, of all
            ! the obs from the base time.  divide by the obs count to get 
            ! the average diff from the base time.
            sec_diff = nint(real(sec_diff) / real(num_avg))
            if (sec_diff > 0) then
               obs_time = increment_time(tpw_base%time, sec_diff)
            else if (sec_diff < 0) then
               obs_time = decrement_time(tpw_base%time, sec_diff)
            else
               obs_time = tpw_base%time
            endif
         else
            lat = tpw_base%lat
            lon = tpw_base%lon
            tpw = tpw_base%tpw_value
            obs_time = tpw_base%time
         endif
 
         ! at this point longitudes are between 0 and 720; subtract
         ! 360 if > 360 to get back into range of 0 to 360.
         if (lon > 360.0_r8) lon = lon - 360.0_r8


         ! optionally select a range of lat/lon.
         ! now that we have the average lat and lon, test them against the limits.
         ! we do this test last so the location of the average lat/lon is guarenteed
         ! to be within the bounds; possibly individual obs being averaged with this
         ! one may be slightly out of the range depending on the box size.
         if (min_lat /= missing_r8) then
            if (lat .lt. min_lat) cycle OBSLOOP2
         endif
         if (max_lat /= missing_r8) then
            if (lat .gt. max_lat) cycle OBSLOOP2
         endif
         if ((min_lon /= missing_r8) .and. (max_lon /= missing_r8)) then
            if (.not. is_longitude_between(lon, min_lon, max_lon)) cycle OBSLOOP2
         endif


         qc   = 0.0_r8  ! quality control - all good.
         hght = 0.0_r8  ! This only works over the ocean; over land need elevation
         oerr = 0.5     ! observation error in cm
         call get_time(obs_time,  osec, oday)  

         call create_3d_obs(lat, lon, hght, VERTISSURFACE, tpw, &
                           obstype, oerr, oday, osec, qc, obs)
         call add_obs_to_seq(obs_seq, obs, obs_time, prev_obs, prev_time, first_obs)
         num_obs = num_obs + 1

      enddo OBSLOOP2

   enddo TIMEBINLOOP


   if (num_obs > 0) then
      write(msgstring, '(A,I12)') 'number of observations in output file is ', num_obs
      call error_handler(E_MSG,'', msgstring, source, revision, revdate)
      call write_obs_seq(obs_seq, outfile)
   else
      call error_handler(E_MSG,'', 'not creating output file because 0 observations found', &
                         source, revision, revdate)
   endif

   call destroy_obs(obs)

enddo DAYLOOP

deallocate(tpw_all, this_bin)

call finalize_utilities()

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
