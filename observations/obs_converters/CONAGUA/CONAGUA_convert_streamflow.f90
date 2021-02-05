! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program CONAGUA_convert_streamflow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! James McCreight jamesmcc at ucar dot edu
! 9/24/2018
! Modified by TJH 12/2019 - use error_handler
!
! Bring in the CONAGUA daily or hourly streamflow (gasto) observations
! 1. Read the meta data file specified in
!       input.nml:CONAGUA_convert_streamflow_nml:meta_data_files
! 2. Process the files listed in
!       intput.nml:CONAGUA_convert_streamflow_nwml:data_file_list
!    All files must either be daily dd*.csv or hourly hd*.csv
!    If a given gage is not in the meta data file, fatal error.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8

use      location_mod, only : VERTISSURFACE

use     utilities_mod, only : initialize_utilities, finalize_utilities,      &
                              open_file, close_file, find_namelist_in_file,  &
                              check_namelist_read, nmlfileunit, do_nml_file, &
                              do_nml_term, error_handler, E_ERR, E_MSG

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

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'CONAGUA_convert_streamflow.f90'
character(len=*), parameter :: revision = '$Revision$'
character(len=*), parameter :: revdate  = '$Date$'

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------

character(len=256) :: meta_data_file         = 'unset'
character(len=256) :: data_file_list         = 'unset'
character(len=256) :: obs_out_file           = 'obs_seq.out'
real(r8)           :: obs_fraction_for_error = 0.01_r8
real(r8)           :: obs_min_err_std        = 0.5_r8
logical            :: debug                  = .false.  ! .true. print info

namelist /CONAGUA_convert_streamflow_nml/  &
     meta_data_file, data_file_list, obs_out_file, &
     obs_fraction_for_error, obs_min_err_std, debug

!------------------------------------------------------------------------

integer :: oday, osec, rcio, iunit, io
integer :: num_copies, num_qc, max_obs
           
logical  :: first_obs

real(r8) :: obs_err, qc
real(r8) :: lat, lon, vert

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time
integer :: meta_unit, list_unit, stn_unit
integer :: key

integer, parameter :: MAX_STATIONS = 5000
integer, parameter :: MAX_RECORDS = 10000
integer :: i_rec, i_stn
character(len=512) :: line_in, stn_name, station_station
integer :: count_stations, count_meta_stations, ii, stn_name_len, wh_meta_stn

logical :: hourly_resolution

integer :: year, month, day
real(r8) :: obs_val
real(r8), parameter :: MISSING_FLOAT = -88888.0_r8
integer,  parameter :: MISSING_INT = -88888

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

type(meta_data_record), dimension(MAX_STATIONS) :: meta_data
integer, parameter :: N_META_FIELDS = 6
character(len=32) :: meta_header(N_META_FIELDS)

! data_file_list
character(len=256) :: data_files(MAX_STATIONS)
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
integer, parameter :: N_STN_DAILY_FIELDS = 33
character(len=32) :: stn_daily_header(N_STN_DAILY_FIELDS)

character(len=512) :: string1, string2, string3

! -------------------------------------------------------
! Start work
! Initialize DART modules.
call initialize_utilities('CONAGUA_convert_streamflow')

call find_namelist_in_file('input.nml', 'CONAGUA_convert_streamflow_nml', iunit)
read(iunit, nml = CONAGUA_convert_streamflow_nml, iostat = io)
call check_namelist_read(iunit, io, 'CONAGUA_convert_streamflow_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=CONAGUA_convert_streamflow_nml)
if (do_nml_term()) write(     *     , nml=CONAGUA_convert_streamflow_nml)

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

if (debug)  then
   write(string1,*) 'opened meta data file "'//trim(meta_data_file)//'"'
   call error_handler(E_MSG, source, string1)
endif

read(meta_unit, *) meta_header

count_meta_stations = 0
meta_loop: do i_rec = 1, MAX_STATIONS
   read(meta_unit, *, iostat=rcio) meta_data(i_rec)
   if (rcio /= 0) then 
      if (debug) then
         write(string1,*)'meta_loop: got bad read code getting line, rcio = ', rcio
         call error_handler(E_MSG, source, string1)
      endif
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
do ii=1, N_META_FIELDS
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

if (debug)  then
   write(string1,*) 'opened data_file_list "'//trim(data_file_list)//'"'
   call error_handler(E_MSG, source, string1)
endif

count_stations = 0
list_loop: do i_rec = 1, MAX_STATIONS
   read(list_unit, '(a)', iostat=rcio) data_files(i_rec)
   if (rcio /= 0) then 
      if (debug) then
         write(string1,*)'list_loop: got bad read code getting line, rcio = ', rcio
         call error_handler(E_MSG, source, string1)
      endif
      exit list_loop
   endif
   count_stations = count_stations + 1
end do list_loop
call close_file(list_unit)

! Validate that all the stations in the list are in the meta data.

valid_loop: do ii = 1, count_stations
   stn_file_name = data_files(ii)
   stn_name_len  = len_trim(stn_file_name)

   ! As far as we can tell, CONAGUA stations are always 5 digits wide.
   stn_resolution = stn_file_name((stn_name_len-10):(stn_name_len-9))
   stn_name       = stn_file_name((stn_name_len-8):(stn_name_len-4))

   if (debug) then
      write(string1,*)'stn_file_name substring is ', stn_file_name((stn_name_len-10):(stn_name_len-4))
      write(string2,*)'station resolution is <',trim(stn_resolution),'>'
      write(string3,*)'station name       is <',trim(stn_name),'>'
      call error_handler(E_MSG, source, string1, text2=string2, text3=string3)
   endif

   if (.not. any(meta_data%name == stn_name)) then
      write(string1,*)'Metadata not found for the requested station file: "'//trim(stn_file_name)//'"'
      call error_handler(E_ERR, source, string1)
   end if
  
   if (ii == 1) then
      hourly_resolution = (stn_resolution == 'hd')
   else
      if (hourly_resolution .and. stn_resolution .ne. 'hd') then
         write(string1,*)"hourly resolution does not appear consistent across the data_file_list.txt"
         call error_handler(E_ERR, source, string1)
      end if
      if (.not. hourly_resolution .and. stn_resolution .ne. 'dd') then
         write(string1,*) "daily resolution does not appear consistent across the data_file_list.txt"
         call error_handler(E_ERR, source, string1)
      end if
   end if

end do valid_loop


! -------------------------------------------------------
! Daily observation processing.

stn_daily_record_null%year  = MISSING_INT
stn_daily_record_null%month = MISSING_INT
stn_daily_record_null%day_of_month_streamflow = MISSING_FLOAT

station_loop: do i_stn = 1, count_stations

   stn_file_name = data_files(i_stn)
   stn_name_len = len(trim(stn_file_name))
   stn_name = stn_file_name((stn_name_len-8):(stn_name_len-4))

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
   
   write(string1,*) '' 
   write(string2,*) '!-----------------------------------'
   write(string3,*) 'Processing station ', i_stn, ': "'//trim(station_station)//'"'
   call error_handler(E_MSG,source,string1,text2=string2,text3=string3)
   
   stn_unit = open_file(trim(stn_file_name), 'formatted', 'read')   

   write(string1,*) 'opened stn record file: "'//trim(stn_file_name)//'"'
   call error_handler(E_MSG,source,string1)

   read(stn_unit, *, iostat=rcio) stn_daily_header
   
   record_loop: do i_rec = 1, MAX_RECORDS
      stn_daily_record = stn_daily_record_null
      read(stn_unit, '(a)', iostat=rcio) line_in
      if (rcio /= 0) then 
         if (debug) then
            write(string1,*) 'recordloop: unable to read record for "'//trim(stn_name)//'"'
            write(string2,*) ' rcio = ', rcio, '; record number = ', i_rec
            call error_handler(E_MSG,source,string1,text2=string2)
         endif
         exit record_loop
      endif

      ! This is a band-aid for records/lines which end with a missing value.
      ! I'd prefer to keep it to a single comma so any "badly mangled" records might cause an error.
      write(line_in, '(a)') trim(line_in) // ','
      
      read(line_in, *, iostat=rcio) stn_daily_record
      if (rcio /= 0) then
         if (debug) then
            write(string1,*) 'recordloop: unable to parse stn_daily_record for "'//trim(stn_name)//'"'
            write(string2,*) 'rcio = ', rcio, '; record number = ', i_rec
            call error_handler(E_MSG,source,string1,text2=string2)
         endif
         exit record_loop
      endif

      year  = stn_daily_record%year
      month = stn_daily_record%month
      
      day_loop: do day = 1, 31

         obs_val = stn_daily_record%day_of_month_streamflow(day)

         if (obs_val == MISSING_FLOAT .or. obs_val <= 0.000001_r8) then
            cycle day_loop ! skip this recrod.
         end if
         
         ! Set the observation error variance here.
         obs_err = max(obs_val*obs_fraction_for_error, obs_min_err_std)
         
         time_obs = set_date(year, month, day, 0, 0, 0)
         ! extract time of observation into gregorian day, sec.
         call get_time(time_obs, osec, oday)

         ! This routine can help tie gages to model spatial elements which can then
         ! be used for more advanced localization.
         ! Currently it is not used for gridded, though a layer in the Fulldom.nc file
         ! could easily facilitate this. The bigger issue with gridded is lack of
         ! connectivity approach currently (only calculated inside the model).
         ! Just supplying dummy values now... 

         call set_streamflow_metadata(key, missing_gage_ID, missing_link_ID)

         call create_3d_obs(                 &
              lat, lon, vert, VERTISSURFACE, &
              obs_val,                       &
              STREAM_FLOW,                   &
              obs_err, oday, osec, qc, obs,  &
              key                            )

         call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs) 

      end do day_loop
      
   end do record_loop

   call close_file(stn_unit)

end do station_loop

write(string1,*)'writing "'//trim(obs_out_file)//'"'
write(string2,*)'obs_count = ', get_num_obs(obs_seq)

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   call error_handler(E_MSG,source,string1,text2=string2)
   call write_obs_seq(obs_seq, obs_out_file)
else
   write(string3,*)'no observations in output. Something is probably wrong.'
   call error_handler(E_ERR,source,string1,text2=string2,text3=string3)
endif

! end of main program
call finalize_utilities()

end program CONAGUA_convert_streamflow

