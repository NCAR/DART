! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> convert snow and ice data center files from binary to obs_seq format.
!> 
!> ;============data description from CC matlab script=================================
!> ; Data are stored as flat two-byte integers representing sea ice concentration
!> ; values. The sea ice concentration data values are packed into integer
!> ; format by multiplying the original sea ice concentration values by
!> ; 10. These values range from 0 to 1,000, with land registered as -800 and
!> ; -100 being the northern hemisphere hole poleward of 85 degrees for the
!> ; SMMR data and poleward of 87 degrees for the SSM/I data, which the
!> ; satellite cannot cover due to orbit inclination.
!> ;
!> ; (nsc -- what are these? ...)
!> ; little-endian byte order no mention of hearder
!> ; One-byte (scaled, unsigned integer) flat binary arrays preceded by a 300-byte header;
!> ; *uint8
!> ;=====================================================================================
!>

program cice_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   cice_to_obs - read binary files from the national snow and ice data
!      center and convert percentages into sea ice concentration observations.
!
!      nsc 28 june 2016 based on conversion programs from yongfei zhang, 
!                       now at UW working with cecelia bitz.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, i2, missing_r8
use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              open_file, close_file, error_handler, E_ERR, &
                              do_nml_file, do_nml_term, nmlfileunit, &
                              find_namelist_in_file, check_namelist_read
use  time_manager_mod, only : time_type, set_calendar_type, set_date, set_time, &
                              operator(>), increment_time, get_date, get_time, &
                              operator(-), GREGORIAN, operator(+), print_date
use      location_mod, only : VERTISLEVEL, VERTISUNDEF
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq
use      obs_kind_mod, only : &
  SAT_SEAICE_AGREG_CONCENTR, &
  SAT_SEAICE_AGREG_VOLUME, &
  SAT_SEAICE_AGREG_SNOWVOLUME, &
  SAT_SEAICE_AGREG_THICKNESS, &
  SAT_SEAICE_AGREG_SNOWDEPTH, &
  SAT_U_SEAICE_COMPONENT, &
  SAT_V_SEAICE_COMPONENT, &
  SAT_SEAICE_CONCENTR, &
  SAT_SEAICE_VOLUME, &
  SAT_SEAICE_SNOWVOLUME


implicit none

!> namelist items
!> @todo : give them reasonable defaults later.
character(len=256) :: cice_lat_file   = 'psn25lats_v3.dat'
character(len=256) :: cice_lon_file   = 'psn25lons_v3.dat'
integer            :: num_latitudes   = 448
integer            :: num_longitudes  = 304
integer            :: start_year      = 1980
integer            :: start_month     = 1
integer            :: start_day       = 1
integer            :: end_year        = 1980
integer            :: end_month       = 1
integer            :: end_day         = 1
real(r8)           :: data_scale_factor = 10.0
real(r8)           :: error_factor    = 0.10
integer            :: land_missing_value = -800
integer            :: pole_missing_value = -100
logical            :: use_data_filename_pattern = .true.
character(len=256) :: data_filename_pattern     = 'bt_YYYYMMDD_n07_v02_n.bin'
character(len=256) :: cice_data_file            = 'bt_19800101_n07_v02_n.bin'
logical            :: use_obsseq_filename_pattern = .true.
character(len=256) :: obsseq_filename_pattern     = 'obs_seq.YYYYMMDD'
character(len=256) :: obsseq_out_file             = 'obs_seq.out'
logical            :: append_to_existing_file     = .false.
logical            :: debug           = .false.


character(len=256) :: input_line, input_filename, next_file, out_file
character(len=32)  :: buf

integer :: oday, osec, rcio, iunit, otype, io
integer :: year, month, day, hour, minute, second
integer :: num_copies, num_qc, max_obs, ilon, ilat, i, j
integer :: start_index
           
logical  :: file_exist, first_obs

real(r8), allocatable :: lat(:,:), lon(:,:), percent(:,:)
real(r8) :: perr, qc
integer(i2), allocatable :: rawdata(:,:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time, curr_time, end_time, one_day

namelist /cice_to_obs_nml/ &
   cice_lat_file,     &
   cice_lon_file,     &
   obsseq_out_file,   &
   num_latitudes,     &
   num_longitudes,    &
   start_year,        &
   start_month,       &
   start_day,         &
   end_year,          &
   end_month,         &
   end_day,           &
   data_scale_factor, &
   use_data_filename_pattern, &
   data_filename_pattern,     &
   cice_data_file,            &
   use_obsseq_filename_pattern, &
   obsseq_filename_pattern,     &
   obsseq_out_file,             &
   append_to_existing_file,   &
   debug


! start of executable code

call initialize_utilities('cice_to_obs')

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

! make space for the data arrays
allocate(lat(num_latitudes,num_longitudes), lon(num_latitudes,num_longitudes), &
         percent(num_latitudes,num_longitudes), rawdata(num_latitudes,num_longitudes))

! read in lats/lons first.  applies to all data files.

! data is binary, distributed as packed short integers
iunit = open_file(cice_lat_file, 'unformatted', 'read')
if (debug) print *, 'opened input file ' // trim(cice_lat_file)
read(iunit) rawdata
lat(:,:) = real(rawdata(:,:),r8) / 10.0_r8
call close_file(iunit)

iunit = open_file(cice_lon_file, 'unformatted', 'read')
if (debug) print *, 'opened input file ' // trim(cice_lon_file)
read(iunit) rawdata
! convert from packed short ints to real values
lon(:,:) = real(rawdata(:,:),r8) / 10.0_r8
call close_file(iunit)

! lon given as -180 and 180 - convert all values at once
where (lon < 0.0_r8)  lon = lon + 360.0_r8  ! changes into 0-360

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
first_obs = .true.

! create a new, empty obs_seq.  you must give a max limit
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
inquire(file=obsseq_out_file, exist=file_exist)
if ( file_exist) then
   if (append_to_existing_file .and. .not. use_obsseq_filename_pattern) then
      call read_obs_seq(obsseq_out_file, 0, 0, max_obs, obs_seq)
   else
      !> @todo FIXME tell user we are ignoring existing file and why
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
      !> @todo make sure start_time = end_time if you aren't
      !> generating the filenames by pattern.
      curr_time = end_time
   endif

   if (use_obsseq_filename_pattern) then
      call fix_filename(obsseq_filename_pattern, obsseq_out_file, year, month, day)
  
      !> @todo FIXME open the new output file?

! create a new, empty obs_seq.  you must give a max limit
! on number of obs.  increase the size if too small.
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

! the first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call set_qc_meta_data(obs_seq, 1, 'Data QC')

      
   else
      ! the obsseq_out_file already has the right value
      !obsseq_out_file = obsseq_data_file
      !> @todo make sure start_time = end_time if you aren't
      !> generating the filenames by pattern.
      curr_time = end_time
   endif

   ! read in concentration data
   iunit =  open_file(next_file, 'unformatted', 'read')
   read(iunit) rawdata
   percent(:,:) = real(rawdata, r8)
   call close_file(iunit)

   ! all observations in this file have the same time?  midnight?
   time_obs = curr_time

   ! extract time of observation into gregorian day, sec.
   call get_time(time_obs, osec, oday)

   ! set missing values where no data or bad data, and skip them 
   ! in the following loop.  doc says:
   !    land_missing_value = -800
   !    pole_missing_value = -100

   where (percent < 0.0_r8 .or. percent > 100.0_r8) percent = missing_r8

   do ilon=1, num_longitudes
      do ilat=1, num_latitudes
 
         if (percent(ilat, ilon) == missing_r8) cycle

         perr = percent(ilat, ilon) * error_factor   ! percentage from namelist

         ! make an obs derived type, and then add it to the sequence
         call create_3d_obs(lat(ilat, ilon), lon(ilat, ilon), 0.0_r8, VERTISUNDEF, &
                            percent(ilat, ilon), SAT_SEAICE_AGREG_CONCENTR, &
                            perr, oday, osec, qc, obs)
         call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

      enddo
   enddo


   ! here is where we decide to loop again for another day or exit if
   ! we have converted the last day's file.
  
   !> @todo FIXME: here if we're looping over the output, write the
   !> obs_seq file and open a new one

 ! if we added any obs to the sequence, write it out to a file now.
 if ( get_num_obs(obs_seq) > 0 ) then
    if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
    call write_obs_seq(obs_seq, obsseq_out_file)
 endif

 !> @todo FIXME destroy the obs_seq so we don't have a mem leak

   curr_time = curr_time + one_day
   if (curr_time > end_time) exit obsloop

   call get_date(curr_time, year, month, day, hour, minute, second)

end do obsloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obsseq_out_file)
endif

!> @todo FIXME destroy the obs_seq so we don't have a mem leak

! end of main program
call finalize_utilities()

contains

!> take a character string which must contain YYYY, MM, and DD
!> and substitute the year, month and day.  it is an error to
!> have a pattern without all 3 items.

subroutine fix_filename(inpattern, outstring, year, month, day)

character(len=*), intent(in)  :: inpattern
character(len=*), intent(out) :: outstring
integer,          intent(in)  :: year
integer,          intent(in)  :: month
integer,          intent(in)  :: day

outstring = inpattern

!> substitute year for YYYY, month for MM, day for DD
!> @todo add error handling

start_index = index(outstring, 'YYYY')
if (start_index < 0) stop
write(buf, '(A)') year
outstring(start_index:start_index+3) = buf(1:4)

start_index = index(outstring, 'MM')
if (start_index < 0) stop
write(buf, '(A)') month
outstring(start_index:start_index+1) = buf(1:2)

start_index = index(outstring, 'DD')
if (start_index < 0) stop
write(buf, '(A)') day
outstring(start_index:start_index+1) = buf(1:2)

end subroutine fix_filename

end program cice_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
