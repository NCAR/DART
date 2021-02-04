! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

!> cice_to_obs - read binary files from the national snow and ice data
!>      center and convert percentages into sea ice concentration 
!>      observation sequence files.
!> 
!>      nsc 28 june 2016 based on conversion programs from Yongfei Zhang, 
!>      now at University of Washington working with Cecelia Bitz.

program cice_to_obs

! ;============data description from CC matlab script=================================
! ; Data are stored as flat two-byte integers representing sea ice concentration
! ; values. The sea ice concentration data values are packed into integer
! ; format by multiplying the original sea ice concentration values by
! ; 10. These values range from 0 to 1,000, with land registered as -800 and
! ; -100 being the northern hemisphere hole poleward of 85 degrees for the
! ; SMMR data and poleward of 87 degrees for the SSM/I data, which the
! ; satellite cannot cover due to orbit inclination.
! ;
! ; (nsc -- what are these? ...)
! ; little-endian byte order no mention of header
! ; One-byte (scaled, unsigned integer) flat binary arrays preceded by a 300-byte header;
! ; *uint8
! ;=====================================================================================

!>@todo right now no code uses netcdf or nc_check, so it's commented out.

use         types_mod, only : r8, i2, i4, missing_r8
use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              open_file, close_file, error_handler, E_ERR, &
                              do_nml_file, do_nml_term, nmlfileunit, &
                              find_namelist_in_file, check_namelist_read, &
                              get_unit
!use  netcdf_utilities_mod, only : nc_check
use  time_manager_mod, only : time_type, set_calendar_type, set_date, set_time, &
                              get_date, get_time, GREGORIAN, &
                              operator(>=), operator(-), operator(+)
use      location_mod, only : VERTISLEVEL, VERTISUNDEF
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, destroy_obs_sequence, & 
                              set_copy_meta_data, set_qc_meta_data
use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq
use      obs_kind_mod, only : SAT_SEAICE_AGREG_CONCENTR, &
                              SAT_SEAICE_AGREG_VOLUME, &
                              SAT_SEAICE_AGREG_SNOWVOLUME, &
                              SAT_SEAICE_AGREG_THICKNESS, &
                              SAT_SEAICE_AGREG_SNOWDEPTH, &
                              SAT_U_SEAICE_COMPONENT, &
                              SAT_V_SEAICE_COMPONENT, &
                              SAT_SEAICE_CONCENTR, &
                              SAT_SEAICE_VOLUME, &
                              SAT_SEAICE_SNOWVOLUME


!>@todo FIXME - originally i thought we needed to read the netcdf file
!> which contained the model land mask.  but creating observations from
!> the binary files is unrelated to the model grid so this code is
!> currently commented out and probably should be removed. BUT: it could 
!> be useful to create an obs_seq.in (synthetic obs) for  
!> a perfect_model_obs experiment, so i've left it for now.   
!> (all related code marked with LANDMASK tag for easier cleanup.)
!> nsc 30jun2016

! netcdf access for land mask LANDMASK
!use netcdf

implicit none

character(len=*), parameter :: routine ='cice_to_obs'

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

! this is based on the model, not on the satellite grid.
! ignore this for now. LANDMASK
!character(len=256) :: land_mask_file  = 'cice_hist.nc'

character(len=256) :: next_file
character(len=512) :: msgstring, msgstring1

integer :: oday, osec, iunit, io, rc
integer :: year, month, day, hour, minute, second
integer :: num_copies, num_qc, max_obs, ilon, ilat
integer :: start_index
! integer :: ncid, varid
           
logical  :: file_exist, first_obs

real(r8), allocatable :: lat(:,:), lon(:,:), percent(:,:)
real(r8) :: perr, qc
! LANDMASK
!integer, allocatable :: tmask(:,:)
integer(i2), allocatable :: rawdata_i2(:)
integer(i4), allocatable :: rawdata_i4(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time, curr_time, end_time, one_day

! these are determined by the format of the binary data files
! and probably shouldn't be changed.
character(len=256) :: cice_lat_file   = 'psn25lats_v3.dat'
character(len=256) :: cice_lon_file   = 'psn25lons_v3.dat'
integer            :: num_latitudes   = 448
integer            :: num_longitudes  = 304
real(r8)           :: grid_scale_factor = 100000.0
real(r8)           :: data_scale_factor = 10.0
integer            :: land_missing_value = -800
integer            :: pole_missing_value = -100

!> the observation error will be set by multiplying
!> the obs value by this factor.  we can also change
!> this to be a fixed value, or ??   but for now the
!> obs error will be 10% of the observed value.
!> (or whatever you specify in the namelist)
real(r8)           :: error_factor    = 0.10

!> science question - like radar obs with no reflectivity,
!> there are many obs with 0 ice.  this is information and
!> these obs could be assimilated to remove ice where it
!> isn't observed.  but it greatly increases the number of
!> obs created.  it should be tested to see if these obs
!> are worth using.
logical            :: ignore_zero_obs = .false.

!> these are alternatives for specifying the input file(s).
!> if pattern is .true. the YYYY, MM, and DD will be
!> replaced by the actual year, month, day and the program
!> will loop from the start year/month/day to the end date. 
!> if pattern is false, it will ignore the start/end values 
!> below and simply open the single cice_data_file.
logical            :: use_data_filename_pattern = .true.
character(len=256) :: data_filename_pattern     = 'bt_YYYYMMDD_n07_v02_n.bin'
character(len=256) :: cice_data_file            = 'bt_19800101_n07_v02_n.bin'

!> these are alternatives for specifying the output file(s).
!> if pattern is .true. the YYYY, MM, and DD will be
!> replaced by the actual year, month, day and the program
!> will create a different output file for each day.
!> if pattern is false, all the obs will be added to a single
!> output obs_seq file with the name given below.
logical            :: use_obsseq_filename_pattern = .true.
character(len=256) :: obsseq_filename_pattern     = 'obs_seq.YYYYMMDD'
character(len=256) :: obsseq_out_file             = 'obs_seq.out'

!> if using an input filename pattern, this program can loop 
!> over multiple days.  this section controls how many days to do.  
!> to do a single day make the start/end values the same or 
!> make the input pattern false.
integer            :: start_year      = 1980
integer            :: start_month     = 1
integer            :: start_day       = 1

integer            :: end_year        = 1980
integer            :: end_month       = 1
integer            :: end_day         = 1

!> probably not as useful here as might be in other converters.
!> if looping outside this program (from a script), you can
!> set this to true to append more obs to an existing obs_seq file.
!> since this program loops already, you normally want to leave
!> this as .false.
logical            :: append_to_existing_file     = .false.

!> set to true to print out more info.
logical            :: debug           = .false.


namelist /cice_to_obs_nml/ &
   cice_lat_file,     &
   cice_lon_file,     &
   num_latitudes,     &
   num_longitudes,    &
   start_year,        &
   start_month,       &
   start_day,         &
   end_year,          &
   end_month,         &
   end_day,           &
   grid_scale_factor, &
   data_scale_factor, &
   error_factor,      &
   land_missing_value,          &
   pole_missing_value,          &
   use_data_filename_pattern,   &
   data_filename_pattern,       &
   cice_data_file,              &
   use_obsseq_filename_pattern, &
   obsseq_filename_pattern,     &
   obsseq_out_file,             &
   append_to_existing_file,     &
   ignore_zero_obs,             &
   debug


! start of executable code

call initialize_utilities(routine)

! Read the namelist entry
call find_namelist_in_file("input.nml", "cice_to_obs_nml", iunit)
read(iunit, nml = cice_to_obs_nml, iostat = io)
call check_namelist_read(iunit, io, "cice_to_obs_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=cice_to_obs_nml)
if (do_nml_term()) write(     *     , nml=cice_to_obs_nml)


! time setup
call set_calendar_type(GREGORIAN)

! files are daily. use time manager routines to loop over days
! between start and end.  initialize the year/month/day so we can
! use it to construct a filename if requested.
curr_time  = set_date(start_year, start_month, start_day, 0, 0, 0)
end_time   = set_date(end_year,   end_month,   end_day,   0, 0, 0)
one_day    = set_time(0, 1)  ! one day

call get_date(curr_time, year, month, day, hour, minute, second)

! LANDMASK
!         tmask(num_latitudes,num_longitudes), &

! make space for the data arrays
allocate(lat(num_latitudes,num_longitudes), &
         lon(num_latitudes,num_longitudes), &
         percent(num_latitudes,num_longitudes), &
         rawdata_i2(num_latitudes*num_longitudes), &
         rawdata_i4(num_latitudes*num_longitudes))


! read in lats/lons first.  applies to all data files.

! data is binary, distributed as scaled 4-byte integers in stream format
! (no record counts).
iunit = open_file(cice_lat_file, 'unformatted', 'read', 'stream')

if (debug) print *, 'opened input latitude file ' // trim(cice_lat_file)
read(iunit) rawdata_i4
lat(:,:) = reshape(real(rawdata_i4(:),r8) / grid_scale_factor, (/ num_latitudes, num_longitudes /) )
call close_file(iunit)

!> see comment in cice_lat_file section
!iunit = open_file(cice_lon_file, 'unformatted', 'read', 'stream')
iunit = get_unit()
open (iunit, file=trim(cice_lon_file), form='unformatted', action='read', &
               position='rewind', access='stream', status='old', iostat=rc)
if (rc /= 0) then
   write(msgstring,*)'Cannot open file "'//trim(cice_lon_file)//'" for reading'
   write(msgstring1,*)'Error code was ', rc
   call error_handler(E_ERR, 'open_file: ', msgstring, source, revision, revdate, &
                      text2=msgstring1)

endif
if (debug) print *, 'opened input longitude file ' // trim(cice_lon_file)
read(iunit) rawdata_i4
lon(:,:) = reshape(real(rawdata_i4(:),r8) / grid_scale_factor, (/ num_latitudes, num_longitudes /) )
call close_file(iunit)

! lon given as -180 and 180 - convert all values at once
where (lon < 0.0_r8)  lon = lon + 360.0_r8  ! changes into 0-360

! ! LANDMASK
! ! read in land mask.  we are going to skip obs over land.
! rc = nf90_open(land_mask_file, nf90_nowrite, ncid)
! call nc_check(rc, routine, 'opening land mask file '//trim(land_mask_file))
! rc = nf90_inq_varid(ncid, "tmask", varid)
! call nc_check(rc, routine, 'inquire var "tmask"')
! rc = nf90_get_var(ncid, varid, tmask)
! call nc_check(rc, routine, 'getting var "tmask"')
! rc = nf90_close(ncid)
! call nc_check(rc, routine, 'closing land mask file '//trim(land_mask_file))
!

! each observation in this series will have a single observation value 
! the max possible number of obs needs to be specified but it will only 
! write out the actual number created.
max_obs    = num_latitudes * num_longitudes
num_copies = 1
num_qc     = 1

! call the initialization code, and initialize two empty observation types
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

! create a new, empty obs_seq.  you must give a max limit
! on number of obs.  increase the size if too small.
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)
first_obs = .true.

! the first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call set_qc_meta_data(obs_seq, 1, 'Data QC')

! if you want to append to existing files (e.g. you have a lot of
! small text files you want to combine), you can do it this way,
! or you can use the obs_sequence_tool to merge a list of files 
! once they are in DART obs_seq format.

! existing file found, append to it
inquire(file=obsseq_out_file, exist=file_exist)
if (file_exist) then
   if (append_to_existing_file .and. .not. use_obsseq_filename_pattern) then
      call read_obs_seq(obsseq_out_file, 0, 0, max_obs, obs_seq)
   else
      !>@todo FIXME tell user we are ignoring existing file and why
   endif
endif

! Set the DART data quality control.   0 is good data. 
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8

obsloop: do    ! no end limit - have the loop break when end time exceeded
   
   if (use_data_filename_pattern) then
      call fix_filename(data_filename_pattern, next_file, year, month, day)
   else
      next_file = cice_data_file
      !>@todo make sure start_time = end_time if you aren't
      !> generating the filenames by pattern.
      curr_time = end_time
   endif

   if (use_obsseq_filename_pattern) then
      call fix_filename(obsseq_filename_pattern, obsseq_out_file, year, month, day)
  
      !>@todo FIXME open the new output file?

      ! create a new, empty obs_seq.  you must give a max limit
      ! on number of obs.  increase the size if too small.
      call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)
      first_obs = .true.

      ! the first one needs to contain the string 'observation' and the
      ! second needs the string 'QC'.
      call set_copy_meta_data(obs_seq, 1, 'observation')
      call set_qc_meta_data(obs_seq, 1, 'Data QC')
      
   else
      ! the obsseq_out_file already has the right value
      !obsseq_out_file = obsseq_data_file
      !>@todo make sure start_time = end_time if you aren't
      !> generating the filenames by pattern.
      curr_time = end_time
   endif

   ! read in concentration data
   !> see comment in cice_lat_file section
   iunit = open_file(next_file, 'unformatted', 'read', 'stream')
   
   if (debug) print *, 'opened data file ' // trim(next_file)
   read(iunit) rawdata_i2
   percent(:,:) = reshape(real(rawdata_i2, r8) / data_scale_factor, (/ num_latitudes, num_longitudes /) )
   call close_file(iunit)

   ! all observations in this file have the same time?  midnight?
   time_obs = curr_time

   ! extract time of observation into gregorian day, sec.
   call get_time(time_obs, osec, oday)

   ! set missing values where no data or bad data, and skip them 
   ! in the following loop.  doc says:
   !    land_missing_value = -800
   !    pole_missing_value = -100
   ! count these to see if they actually exist

   ! make this an error handler call?
   if (debug) print *, 'total obs = ', size(percent), &
              ' out of range = ',  count (percent < 0.0_r8 .or. percent > 100.0_r8) , &
              ' zeros = ', count (percent == 0.0_r8)

   where (percent < 0.0_r8 .or. percent > 100.0_r8) percent = missing_r8

   do ilon=1, num_longitudes
      do ilat=1, num_latitudes
 
         ! here is where we decide to ignore obs or not
         !if (lat(ilat,ilon) < 40 ) cycle

         if (percent(ilat, ilon) == missing_r8) cycle

         if (percent(ilat, ilon) /= 0.0_r8) then
            perr = percent(ilat, ilon) / 100.0_r8 * error_factor   ! percentage from namelist
         else
            if (ignore_zero_obs) cycle
            perr = 0.05 !error_factor ! if obs value is 0, give a non-zero error
         endif

         ! LANDMASK
         !if (tmask(ilat, ilon) /= 1) cycle   ! land point

         ! make an obs derived type, and then add it to the sequence
         call create_3d_obs(lat(ilat, ilon), lon(ilat, ilon), 0.0_r8, VERTISUNDEF, &
                            percent(ilat, ilon)/100.0_r8, SAT_SEAICE_AGREG_CONCENTR, &
                            perr, oday, osec, qc, obs)
         call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

      enddo
   enddo


   ! here is where we decide to loop again for another day or exit if
   ! we have converted the last day's file.
  
   !>@todo FIXME: here if we're looping over the output, write the
   !> obs_seq file and open a new one

   ! if we added any obs to the sequence, write it out to a file now.
   if ( get_num_obs(obs_seq) > 0 .and. use_obsseq_filename_pattern ) then
      if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
      call write_obs_seq(obs_seq, obsseq_out_file)
   endif

   ! destroy the obs_seq so we don't have a mem leak
   call destroy_obs_sequence(obs_seq)

   if (curr_time >= end_time) exit obsloop
   curr_time = curr_time + one_day

   call get_date(curr_time, year, month, day, hour, minute, second)

end do obsloop

! if we haven't already written the obs_seq file inside the loop above, write it now
if ( .not. use_obsseq_filename_pattern .and. get_num_obs(obs_seq) > 0 ) then
   if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obsseq_out_file)
endif

! clean up
call destroy_obs_sequence(obs_seq)

! LANDMASK
!          tmask, 

deallocate(lat, lon, percent, rawdata_i2, rawdata_i4)

! end of main program
call finalize_utilities()

contains

!-----------------------------------------------------------------------
!> take a character string which must contain YYYY, MM, and DD
!> and substitute the year, month and day.  it is an error to
!> have a pattern without all 3 items.

subroutine fix_filename(inpattern, outstring, year, month, day)

character(len=*), intent(in)  :: inpattern
character(len=*), intent(out) :: outstring
integer,          intent(in)  :: year
integer,          intent(in)  :: month
integer,          intent(in)  :: day

character(len=32)  :: buf

outstring = inpattern

! substitute year for YYYY, month for MM, day for DD
!>@todo add error handling

start_index = index(outstring, 'YYYY')
if (start_index < 0) then
   write(msgstring,*)'String "YYYY" was not found in pattern "'//trim(inpattern)//'"'
   write(msgstring1, *)'Trying to replace YYYY with year ', year
   call error_handler(E_ERR, 'cice_to_obs: ', msgstring, source, revision, revdate, &
                      text2=msgstring1)
endif
write(buf, '(I4.4)') year
outstring(start_index:start_index+3) = buf(1:4)

start_index = index(outstring, 'MM')
if (start_index < 0) then
   write(msgstring,*)'String "MM" was not found in pattern "'//trim(inpattern)//'"'
   write(msgstring1, *)'Trying to replace MM with month ', month
   call error_handler(E_ERR, 'cice_to_obs: ', msgstring, source, revision, revdate, &
                      text2=msgstring1)
endif
write(buf, '(I2.2)') month
outstring(start_index:start_index+1) = buf(1:2)

start_index = index(outstring, 'DD')
if (start_index < 0) then
   write(msgstring,*)'String "DD" was not found in pattern "'//trim(inpattern)//'"'
   write(msgstring1, *)'Trying to replace DD with day ', day
   call error_handler(E_ERR, 'cice_to_obs: ', msgstring, source, revision, revdate, &
                      text2=msgstring1)
endif
write(buf, '(I2.2)') day
outstring(start_index:start_index+1) = buf(1:2)

end subroutine fix_filename

end program cice_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
