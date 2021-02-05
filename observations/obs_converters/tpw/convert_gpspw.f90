! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_gpspw

! convert ground-based GPS observations of Total Precipitable Water 
! into DART observation sequence files, that is archived in the
! SuomiNet website (http://www.suominet.ucar.edu/data/index.html).
! THIS CODE ASSUMES ALL ELEVATIONS ARE AT 0.0 meters.
! this is true for all ocean obs.  if you have land obs, you will
! have to get the hght of the surface elevation at the lat/lon and
! use that instead.
!
! assumes two input text files:
! the first one is the data which match the following read line:
!      read(iunit, '(A4, 1x, i4, 2i2, 1x, 2i2, 5x, f7.1)' &
!                stnid, iyear, imonth, iday, ihour, imin
! 
! the second one is a list of station locations.
!      read(junit, '(2f11.6,f8.2,'(A4)') &
!           lon, lat, hgt, stnid
! 
! constructs a input filename based on:
!  ObsBase/InfilePrefix + YYYYMMDD + InfileSuffix
! 
! constructs a output filename based on:
!  ./OutfilePrefix + YYYYMMDD + OutfileSuffix
!
! any of the prefix or suffixes can be '' (blank)
! 

use         types_mod, only : r8, metadatalength, missing_r8
use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time,&
                              increment_time, get_time, set_date, operator(-),  &
                              print_date, decrement_time, operator(>),          &
                              print_time, julian_day
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
use      obs_kind_mod, only :   GPS_TOTAL_PRECIPITABLE_WATER,  &
                              get_index_for_type_of_obs
use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq


use           netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
 "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!--------------------------------------


type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs 
type(time_type)         :: obs_time, prev_time
type(time_type)         :: base_time, time_diff

integer  :: idate, junit, ihr, irec
integer  :: ibin, nr, io, obstype
integer  :: iyear, imonth, iday, ihour, imin
integer  :: nrec, idd, iunit, num_obs, num_avg
integer  :: nbin, kk, nn, oday, osec, sec_diff

integer, parameter :: ncol       = 9
integer, parameter :: num_copies = 1
integer, parameter :: num_qc     = 1

real(r8) :: oerr, qc
real(r8) :: lon, lat, hght, tpw

logical  :: first_obs, lon_between

character(len = 4)   :: stnid, locid
character(len = 19)  :: obsdate, indate
character(len = 256) :: msgstring, infile, outfile

! items in namelist, along with default values
integer  :: start_year  = 2008
integer  :: start_month = 1
integer  :: start_day   = 1
integer  :: total_days  = 31
integer  :: max_obs     = 150000
real(r8) :: min_lon = missing_r8
real(r8) :: max_lon = missing_r8
real(r8) :: min_lat = missing_r8
real(r8) :: max_lat = missing_r8
! the date, in 'YYYYMMDD' format, will be inserted between
! the input and output file prefix and suffix.  ObsBase is
! only prepended to the input file
character(len=128) :: ObsBase       = '../data'
character(len=128) :: StationFile   = '../data/Suominet_North_America_stations_LatLon.table.txt'
character(len=64)  :: InfilePrefix  = 'SUOh_'
character(len=64)  :: InfileSuffix  = '.PWV'
character(len=64)  :: OutfilePrefix = 'obs_seq.'
character(len=64)  :: OutfileSuffix = ''
character(len=32)  :: observation_name = 'GPS_TOTAL_PRECIPITABLE_WATER'

namelist /convert_gpspw_nml/ start_year, start_month, start_day, &
      total_days, max_obs,                              &
      min_lon, max_lon, min_lat, max_lat,                      &
      ObsBase, InfilePrefix, InfileSuffix, OutfilePrefix,      &
      OutfileSuffix, observation_name, StationFile


! ----------------------------------------------------------------------
! start of executable program code
! ----------------------------------------------------------------------

call initialize_utilities('convert_gpspw')
call register_module(source,revision,revdate)

! Initialize the obs_sequence module 

call static_init_obs_sequence()

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file('input.nml', 'convert_gpspw_nml', iunit)
read(iunit, nml = convert_gpspw_nml, iostat = io)
call check_namelist_read(iunit, io, 'convert_gpspw_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=convert_gpspw_nml)
if (do_nml_term()) write(    *      , nml=convert_gpspw_nml)

! some simple error checking - if you specify either of the longitude
! limits, since it's a cyclic value you have to specify both.  the included
! area will start at min_lon and go east until max_lon.   for latitude
! you can let one default if you only want to set a min or a max.
if (((min_lon /= missing_r8) .and. (max_lon == missing_r8)) .or. &
    ((min_lon == missing_r8) .and. (max_lon /= missing_r8))) then
   write(msgstring, *) 'if you set the min or max limit on longitude you must specify both'
   call error_handler(E_ERR,'convert_gpspw', msgstring, source, revision, revdate)
endif

! convert a string into an observation type number
obstype = get_index_for_type_of_obs(observation_name)
if (obstype < 0) then
   write(msgstring, *) 'unrecognized observation type ', trim(observation_name)
   call error_handler(E_ERR,'convert_gpspw', msgstring, source, revision, revdate, &
                      text2='check the obs_defs in your &preprocess list')
endif

call set_calendar_type(GREGORIAN)

!--------------------------------------
! Read the station list first
junit = open_file(StationFile,'formatted',action='read')

! FIXME: Assume all the input files have the same year (for now).
idate = julian_day(start_year, start_month, start_day)
write(*,*) 'idate = ',idate

DAYLOOP: do idd = 0, total_days-1

   ! set up for next obs_seq for next day

   call init_obs(obs,      num_copies, num_qc)
   call init_obs(prev_obs, num_copies, num_qc)
   first_obs = .true.
   num_obs = 0

   call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

   call set_copy_meta_data(obs_seq, 1, 'GPS PW observation')
   call set_qc_meta_data(obs_seq, 1, 'DATA QC')

   write(obsdate, '(i4.4, 2i2.2)') start_year, start_month, start_day + idd
   write(outfile, '(A)') trim(OutfilePrefix)//trim(obsdate)//trim(OutfileSuffix)
   write(*, '(A)') 'creating file: '//trim(outfile)

   nrec = 0

   ! read in all observation records 
   HOURLOOP: do ihr = 0, 23

   ! construct the input and output filenames
   write(indate, '(i4.4, A1, i3.3, A1, i2.2, A3)') start_year, '.', idate + idd, '.', ihr, '.00'
   write(infile, '(A)') trim(ObsBase)//'/'//trim(InfilePrefix)//trim(indate)//trim(InfileSuffix)
   write(*, '(A)') 'reading  file: '//trim(infile)

   iunit = open_file(infile,'formatted',action='read')
   read(iunit,*) 
   read(iunit,*) 

   READLOOP: do 

      read(iunit, '(A4, 1x, i4, 2i2, 1x, 2i2, 5x, f7.1)', iostat=io) &
                    stnid, iyear, imonth, iday, ihour, imin, tpw
      if (io /= 0) exit READLOOP

      rewind(junit)
      STNLOOP: do 
         read(junit, '(f11.6,f10.6,f8.2,1x,A4)', iostat=io) lon, lat, hght, locid
         if(locid.eq.stnid) exit STNLOOP
         if(io /= 0) exit STNLOOP
      enddo STNLOOP
      if(locid.ne.stnid) then
         write(msgstring, *) 'Could not find the station: ',stnid,' vs. ',locid
         call error_handler(E_ERR,'convert_gpspw', msgstring, source, revision, revdate)
      endif

      nrec = nrec + 1

      if(nrec > max_obs) then
         write(msgstring, *) 'number of input observations larger than "max_obs" limit'
         call error_handler(E_ERR,'convert_gpspw', msgstring, source, revision, revdate, &
                            text2='increase parameter "max_obs" in namelist and rerun')
      endif

      ! adjust longitude ranges to [0,360]
      if (lon < 360.0_r8) lon = lon + 360.0_r8

         ! optionally select a range of lat/lon.
         ! now that we have the average lat and lon, test them against the limits.
         ! we do this test last so the location of the average lat/lon is guarenteed
         ! to be within the bounds; possibly individual obs being averaged with this
         ! one may be slightly out of the range depending on the box size.
         if (min_lat /= missing_r8) then
            if (lat .lt. min_lat) cycle READLOOP
         endif
         if (max_lat /= missing_r8) then
            if (lat .gt. max_lat) cycle READLOOP
         endif
         if ((min_lon /= missing_r8) .and. (max_lon /= missing_r8)) then
            if (.not. is_longitude_between(lon, min_lon, max_lon)) cycle READLOOP
         endif

         qc   = 0.0_r8  ! quality control - all good.
        !hght = 0.0_r8  ! This only works over the ocean; over land need elevation
         oerr = 1.5     ! observation error in cm
         obs_time = set_date(iyear, imonth, iday, ihour, imin, 0)
         call get_time(obs_time,  osec, oday)  

         call create_3d_obs(lat, lon, hght, VERTISSURFACE, tpw, &
                           obstype, oerr, oday, osec, qc, obs)
         call add_obs_to_seq(obs_seq, obs, obs_time, prev_obs, prev_time, first_obs)
         num_obs = num_obs + 1

   enddo READLOOP

   write(msgstring, '(A,I12)') 'number of observations in input file is ', nrec
   call error_handler(E_MSG,'', msgstring, source, revision, revdate)

   call close_file(iunit)

   enddo HOURLOOP


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

call finalize_utilities()

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
