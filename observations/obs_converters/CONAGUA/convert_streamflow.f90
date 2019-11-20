! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_streamflow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! James McCreight jamesmcc at ucar dot edu
! 9/24/2018
!
! Bring in the CONAGUA daily or hourly streamflow (gasto) observations
! 1. Read the meta data file specified in
!       input.nml:convert_streamflow_nml:meta_data_files
! 2. Process the files listed in
!       intput.nml:convert_streamflow_nwml:data_file_list
!    All files must either be daily dd*.csv or hourly hd*.csv
!    If a given gage is not in the meta data file, fatal error.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8

use      location_mod, only : VERTISSURFACE

use     utilities_mod, only : initialize_utilities, finalize_utilities,      &
                              open_file, close_file, find_namelist_in_file,  &
                              check_namelist_read, nmlfileunit, do_nml_file, &
                              do_nml_term

use  time_manager_mod, only : time_type, set_calendar_type, set_date, set_time, &
                              operator(>=), increment_time, get_time,           &
                              operator(-), GREGORIAN, operator(+), print_date

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,     &
                              static_init_obs_sequence, init_obs,            &
                              write_obs_seq, init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data

use      obs_kind_mod, only : STREAM_FLOW

use obs_def_streamflow_mod, only : set_streamflow_metadata, &
                                   missing_gage_ID, missing_link_ID

use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq

implicit none

character(len=64), parameter :: obs_out_file    = 'obs_seq.out'

integer :: oday, osec, rcio, iunit, io
integer :: num_copies, num_qc, max_obs, ix, iy
           
logical  :: file_exist, first_obs

real(r8) :: obs_err, qc
real(r8) :: lat, lon, vert
real(r8) :: dlon, dlat
real(r8), allocatable :: coverage(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

real(r8) :: missing_value = -20.0_r8
logical  :: debug = .false.  ! set to .true. to print info


! -------------------------------------------------------
! Specific CONAGUA additions

! Namelist
character(len=256) :: meta_data_file, data_file_list
real(r8)           :: obs_fraction_for_error = 0.01
real(r8)           :: obs_min_err_std = 0.5

namelist /convert_streamflow_nml/  &
     meta_data_file, data_file_list, missing_value, debug, &
     obs_fraction_for_error, obs_min_err_std

integer :: meta_unit, list_unit, stn_unit
integer :: key

integer, parameter :: max_stations = 5000
integer, parameter :: max_records = 10000
integer :: i_rec, i_stn
character(len=512) :: line_in, stn_name, station_station
integer :: count_stations, count_meta_stations, ii, dd, stn_name_len, wh_meta_stn

logical :: hourly_resolution

integer :: year, month, day
real :: obs_val

real(r8) :: missing_float = -88888.0_r8
integer :: missing_int = -88888

! Meta data derived type/structure.
!FID,LON,LAT,STATION,Name,HOST
!1,-92.75,17.38333,Oxolotan,30111,BANDAS
type :: meta_data_record
   integer :: fid
   real :: lon
   real :: lat
   character(len=64) :: station
   character(len=64) :: name
   character(len=64) :: host
end type meta_data_record
type(meta_data_record), dimension(max_stations) :: meta_data
integer, parameter :: n_meta_fields = 6
character(len=32) :: meta_header(n_meta_fields)

! data_file_list
character(len=256) :: data_files(max_stations)
character(len=256) :: stn_file_name, stn_resolution

! Daily data derived type/structure.
type :: station_daily_data_record
   integer :: year
   integer :: month
   real :: day_of_month_streamflow(31)
end type station_daily_data_record
! Will process one in and one out, so the dimension is only len 1.
type(station_daily_data_record) :: stn_daily_record
type(station_daily_data_record) :: stn_daily_record_null
integer, parameter :: n_stn_daily_fields = 33
character(len=32) :: stn_daily_header(n_stn_daily_fields)

! Hourly data derived type/structure.
type :: station_hourly_data_record
   integer :: day
   integer :: month
   integer :: year
   integer :: hour
   real :: streamflow
end type station_hourly_data_record
! Will process one in and one out, so the dimension is only len 1.
type(station_hourly_data_record) :: stn_hourly_record
type(station_hourly_data_record) :: stn_hourly_record_null
integer, parameter :: n_stn_hourly_fields = 3
character(len=32) :: stn_hourly_header(n_stn_hourly_fields)
! End declarations


! -------------------------------------------------------
! Start work
! Initialize DART modules.
call initialize_utilities('convert_streamflow')

call find_namelist_in_file('input.nml', 'convert_streamflow_nml', iunit)
read(iunit, nml = convert_streamflow_nml, iostat = io)
call check_namelist_read(iunit, io, 'convert_streamflow_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=convert_streamflow_nml)
if (do_nml_term()) write(     *     , nml=convert_streamflow_nml)

! time setup
call set_calendar_type(GREGORIAN)

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

! Several obs_seq quantities which are static for our observations:

! surface obs.  normally we set vert to the actual surface elevation, 
! we do not have it here, so set to 0 for now.
vert = 0.0_r8   

! Set the DART data quality control.   0 is good data. 
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8


! -------------------------------------------------------
! Read meta data file 
meta_data_file = trim(meta_data_file)
meta_unit = open_file(meta_data_file, 'formatted', 'read')
if (debug) print *, 'opened meta data file ' // trim(meta_data_file)

read(meta_unit, *) meta_header

count_meta_stations = 0
meta_loop: do i_rec = 1, max_stations
   read(meta_unit, *, iostat=rcio) meta_data(i_rec)
   if (rcio /= 0) then 
      if (debug) print *, 'meta_loop: got bad read code getting line, rcio = ', rcio
      exit meta_loop
   endif
   count_meta_stations = count_meta_stations + 1
end do meta_loop
call close_file(meta_unit)

print *, ''
print *, '----------------------------------'
print *, 'reading meta data file: ', trim(meta_data_file)
print *, ''
print *, 'The following header fields were found: '
do ii=1, n_meta_fields
    write(*, fmt="(1x,a)", advance="no") trim(meta_header(ii)) // ', '
end do
print *,''
print *,''
print *, "count_meta_stations: ", count_meta_stations
do ii = 1, count_meta_stations
   print*, trim(meta_data(ii)%name), ': ', trim(meta_data(ii)%station)
end do
print *,''   


! -------------------------------------------------------
! Pre-Process the data_file_list and check it's consistent with
! 1) meta_data
! 2) itself (all files at the same time resolution)
list_unit = open_file(data_file_list, 'formatted', 'read')
count_stations = 0
list_loop: do i_rec = 1, max_stations
   read(list_unit, '(a)', iostat=rcio) data_files(i_rec)
   if (rcio /= 0) then 
      if (debug) print *, 'list_loop: got bad read code getting line, rcio = ', rcio
      exit list_loop
   endif
   count_stations = count_stations + 1
end do list_loop
call close_file(list_unit)

! Validate that all the stations in the list are in the meta data.
valid_loop: do ii = 1, count_stations
   stn_file_name = data_files(ii)
   !print *, "stn_file_name: ", stn_file_name 
   stn_name_len = len(trim(stn_file_name))
   !print *, "stn_name_len: ", stn_name_len
   ! As far as we can tell, CONAGUA stations are always 5 digits wide.
   stn_name = stn_file_name((stn_name_len-8):(stn_name_len-4))
   stn_resolution = stn_file_name((stn_name_len-10):(stn_name_len-9))
   !print *,''
   !print *,'stn_name: ', trim(stn_name), '  ;  resolution: ', trim(stn_resolution)

   if (.not. any(meta_data%name .eq. stn_name)) then
      print *, 'Metadata not found for the requested station file: ', stn_file_name
      stop 1
   end if
  
   if (ii .eq. 1) then
      hourly_resolution = (stn_resolution .eq. 'hd')
   else
      if (hourly_resolution .and. stn_resolution .ne. 'hd') then
         print *, "Time resolution does not appear consistent across the data_file_list.txt"
         stop 2
      end if
      if (.not. hourly_resolution .and. stn_resolution .ne. 'dd') then
         print *, "Time resolution does not appear consistent across the data_file_list.txt"
         stop 3
      end if
   end if

end do valid_loop


! -------------------------------------------------------
! Daily observation processing.

stn_daily_record_null%year = missing_int
stn_daily_record_null%month = missing_int
stn_daily_record_null%day_of_month_streamflow = missing_float

station_loop: do i_stn = 1, count_stations

   stn_file_name = data_files(i_stn)
   stn_name_len = len(trim(stn_file_name))
   stn_name = stn_file_name((stn_name_len-8):(stn_name_len-4))
   ! My pack and wheres are rusty
   wh_meta_stn = 0
   wh_loop: do ii = 1, count_meta_stations
      if (trim(meta_data(ii)%name) .eq. trim(stn_name)) then
         wh_meta_stn = ii
         exit wh_loop
      end if
   end do wh_loop

   station_station = meta_data(wh_meta_stn)%station
   lon = meta_data(wh_meta_stn)%lon
   if (lon < 0.0_r8) lon = lon + 360.0_r8
   lat = meta_data(wh_meta_stn)%lat
   !print *, wh_meta_stn
   
   print *, '' 
   print *, '!-----------------------------------'
   print *, 'Processing station ', i_stn, ': ', trim(station_station)
   
   stn_unit = open_file(trim(stn_file_name), 'formatted', 'read')   
   print *, 'opened stn record file: ' // trim(stn_file_name)
   read(stn_unit, *, iostat=rcio) stn_daily_header
   !print *, stn_daily_header
   
   record_loop: do i_rec = 1, max_records
      stn_daily_record = stn_daily_record_null
      read(stn_unit, '(a)', iostat=rcio) line_in
      if (rcio /= 0) then 
         if (debug) print *, 'recordloop: got bad read code getting line_in, rcio = ', rcio, &
                             '; record number = ', i_rec
         exit record_loop
      endif

      ! This is a band-aid for records/lines which end with a missing value.
      ! I'd prefer to keep it to a single comma so any "badly mangled" records might cause an error.
      write(line_in, '(a)') trim(line_in) // ','
      
      read(line_in, *, iostat=rcio) stn_daily_record
      if (rcio /= 0) then
         if (debug) print *, 'recordloop: got bad read code parsing line_in, rcio = ', rcio, &
                             '; record number = ', i_rec
         exit record_loop
      endif

      year = stn_daily_record%year
      month = stn_daily_record%month
      
      day_loop: do dd = 1, 31

         day = dd

         obs_val = stn_daily_record%day_of_month_streamflow(dd)
         ! Set the observation error variance here.
         obs_err = max(obs_val*obs_fraction_for_error, obs_min_err_std)

         if (obs_val .eq. missing_float .or. obs_val <= 0.000001) then
            cycle day_loop ! skip this recrod.
         end if
         
         print *, year, month, day, obs_val
         
         time_obs = set_date(year, month, day, 0, 0, 0)
         ! extract time of observation into gregorian day, sec.
         call get_time(time_obs, osec, oday)

         ! print *,'lat: ', lat
         ! print *,'lon: ', lon
         ! !print
         ! print *,'dble(obs_val): ', dble(obs_val)
         ! print *,'obs_err: ', obs_err
         ! print *,'oday: ', oday
         ! print *,'osec: ', osec
         ! print *,'qc: ', qc

         ! This routine can help tie gages to model spatial elements which can then
         ! be used for more advanced localization.
         ! Currently it is not used for gridded, though a layer in the Fulldom.nc file
         ! could easily facilitate this. The bigger issue with gridded is lack of
         ! connectivity approach currently (only calculated inside the model).
         ! Just supplying dummy values now... 
         call set_streamflow_metadata(key, missing_gage_ID, missing_link_ID)

         call create_3d_obs(                 &
              lat, lon, vert, VERTISSURFACE, &
              dble(obs_val),                 &
              STREAM_FLOW,                   &
              obs_err, oday, osec, qc, obs,  &
              key                            &
         )

         call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs) 

         ! TODO(JLM): This is just demonstrating that the write_obs_seq below fails 
         ! for just one record remove when that is fixed
         ! exit station_loop 
         
      end do day_loop
      
   end do record_loop

   call close_file(stn_unit)

end do station_loop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

end program convert_streamflow

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
