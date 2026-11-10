! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! The following is an obs converter for Copernicus Sea Surface 
! Temperature L3S product. The data is based on a blended global 
! high resolution ODYSSEA SST Multi-sensor L3 composite. 

! The converter reads in the adjusted SST (bias corrected) 
! and its associated error sd from the incoming netcdf file. 
! A lower bound on the SD is applied to avoid over-fitting 
! a relatively raw data product. Observations from multiple
! files (at different times) are collected and added to the
! observation sequence file. 

! The observed temperature is converted from K to C. 


program cmems_sst_to_obs

use types_mod,            only : i4, i8, r8, t_kelvin, MISSING_R8, digits12
use time_manager_mod,     only : time_type, set_calendar_type, GREGORIAN, get_time,   &   
                                 set_date, set_time, print_date, operator(+), operator(-)
use utilities_mod,        only : initialize_utilities, find_namelist_in_file, E_MSG,  &   
                                 nmlfileunit, error_handler, do_nml_term, E_ERR,      &   
                                 finalize_utilities, do_nml_file, get_next_filename,  &
                                 find_textfile_dims, file_exist, check_namelist_read 
use location_mod,         only : VERTISSURFACE
use obs_sequence_mod,     only : obs_type, obs_sequence_type, init_obs, get_num_obs,  &
                                 static_init_obs_sequence, init_obs_sequence,         &   
                                 set_copy_meta_data, set_qc_meta_data, write_obs_seq, &
                                 destroy_obs_sequence, destroy_obs
use obs_utilities_mod,    only : create_3d_obs, add_obs_to_seq
use obs_kind_mod,         only : SATELLITE_BLENDED_SST
use netcdf_utilities_mod, only : nc_check, nc_open_file_readonly, nc_close_file,      &   
                                 nc_get_variable, nc_get_attribute_from_variable,     &   
                                 nc_get_dimension_size   
use netcdf

implicit none

character(len=*), parameter :: source = 'cmems_sst_to_obs'

! Lower bound for obs_error_sd
real(r8), parameter :: OBS_ERROR_SD_MIN   = 0.05_r8 ! Add this to the namelist? 
integer, parameter  :: OBS_QC_LOW_QUALITY = 3

! File variables
character(len=512) :: string1
character(len=256) :: next_infile
integer            :: io, iunit, filenum
integer            :: num_new_obs, nfiles
logical            :: first_obs = .true.

! Obs sequence
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time, prev_time
integer                 :: num_copies = 1, &     ! number of copies in sequence
                           num_qc     = 1        ! number of QC entries
real(r8)                :: obs_qc     = 0.0_r8

! Data arrays
real(r8), allocatable   :: lon(:), lat(:)     ! Location 
real(r8), allocatable   :: sst(:,:), err(:,:) ! Obs and associated error
real(r8), allocatable   :: iqc(:,:)           ! Incoming data QC 

!------------------------------------------------------------------------
!  Declare namelist parameters
character(len=256) :: file_list         = 'sst_file_list.txt'
character(len=256) :: file_out          = 'obs_seq.sst'
integer            :: avg_obs_per_file  = 500000
logical            :: debug             = .true.


namelist /cmems_sst_to_obs_nml/ file_list,        &
                                file_out,         &
                                avg_obs_per_file, &
                                debug

! Start Converter
call initialize_utilities()

! Read the namelist options
call find_namelist_in_file('input.nml', 'cmems_sst_to_obs_nml', iunit)
read(iunit, nml = cmems_sst_to_obs_nml, iostat = io)
call check_namelist_read(iunit, io, 'cmems_sst_to_obs_nml')

if (do_nml_file()) write(nmlfileunit, nml=cmems_sst_to_obs_nml)
if (do_nml_term()) write(     *     , nml=cmems_sst_to_obs_nml)

! Set the calendar kind
call set_calendar_type(GREGORIAN)

! Get number of observations
call find_textfile_dims(file_list, nfiles)
num_new_obs = avg_obs_per_file * nfiles

! Initialize obs seq file
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

if (file_exist(file_out)) then
   string1 = 'Output file: '//trim(adjustl(file_out))//' exists. Replacing it ...'
   call error_handler(E_MSG, source, string1)
else
   string1 = 'Creating "'//trim(adjustl(file_out))//'" file.'
   call error_handler(E_MSG, source, string1)
endif

call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)
call set_copy_meta_data(obs_seq, num_copies, 'observation')
call set_qc_meta_data(obs_seq, num_qc, 'QC')

! Loop over the obs files
do filenum = 1, nfiles
  next_infile = get_next_filename(file_list, filenum)
  
  if (debug) write(*, '(/, A, i0, X, A)') 'Input file: #', filenum, trim(next_infile)

  ! Collect data from netcdf file
  call collect_data(next_infile, obs, prev_obs, prev_time)
  call cleanup()
enddo 

call finish_obs()
call finalize_utilities()

contains 

!---------------------------------------------------------------
! Read SST and associated metadata and add them to the obs
! sequence file 
subroutine collect_data(datafile, obs, prev_obs, prev_time)

character(len=*), parameter :: routine = 'collect_data'

character(len=*), intent(in)    :: datafile
type(obs_type),   intent(inout) :: obs, prev_obs
type(time_type),  intent(inout) :: prev_time

integer           :: ncid, nlon, nlat, ntime 
integer           :: ilon, ilat
real(r8)          :: missing_sst

integer           :: year, month, day, hour, minute, second, ios
real(digits12)    :: time_since_init
character(len=64) :: datestr
type(time_type)   :: base_date
integer(i8)       :: big_integer
integer           :: some_seconds, some_days

integer           :: osec, oday

ncid  = nc_open_file_readonly(datafile)
ntime = nc_get_dimension_size(ncid, 'time', routine)

if (ntime /= 1) then
   string1 = '"time" dimension > 1 not supported.'
   call error_handler(E_ERR, routine, string1, source)
endif

! Dimensions
nlat = nc_get_dimension_size(ncid, 'latitude', routine)
nlon = nc_get_dimension_size(ncid, 'longitude', routine)

! Memory
allocate(lat(nlat), lon(nlon)) 
allocate(sst(nlon, nlat), iqc(nlon, nlat), err(nlon, nlat))

! Location data
call nc_get_variable(ncid, 'latitude', lat, routine)
call nc_get_variable(ncid, 'longitude', lon, routine)
where(lon < 0.0_r8) lon = lon + 360.0_r8

! SST
call nc_get_variable(ncid, 'adjusted_sea_surface_temperature', sst, routine)
call nc_get_attribute_from_variable(ncid, 'adjusted_sea_surface_temperature', &
                                    '_FillValue', missing_sst, routine)

! Incoming QC 
! 0: no_data, 1: bad_data, 2: worst_quality, 3: low_quality, 
! 4: acceptable_quality, 5: best_quality 
call nc_get_variable(ncid, 'quality_level', iqc, routine)

! Errors
call nc_get_variable(ncid, 'sses_standard_deviation', err, routine)

do ilat = 1, nlat
   do ilon = 1, nlon
      if (err(ilon, ilat) /= err(ilon, ilat)) err(ilon, ilat) = OBS_ERROR_SD_MIN            

      ! Clamp the errors
      if (sst(ilon, ilat) /= missing_sst .and. sst(ilon, ilat) == sst(ilon, ilat)) then 
         err(ilon, ilat) = max(err(ilon, ilat), OBS_ERROR_SD_MIN)
      endif
   enddo
enddo

! Time 
call nc_get_variable(ncid, 'time', time_since_init)
call get_time_units(ncid, datestr)

!call nc_get_attribute_from_variable(ncid, 'time', 'units', datestr, routine)

read(datestr, '(14x, i4, 5(1x,i2))', iostat=ios) year, month, day, hour, minute, second
if (ios /= 0) then
   write(string1, *) 'Unable to read time variable units. Error status was ', ios
   call error_handler(E_ERR, routine, string1, source)
endif
big_integer  = int(time_since_init, i8)

base_date    = set_date(year, month, day, hour, minute, second)
some_days    = big_integer / 86400
big_integer  = big_integer - (some_days * 86400)
some_seconds = int(big_integer, i4)

obs_time = base_date + set_time(some_seconds, some_days)
if (debug) call print_date(obs_time, str= '  * Reading SST; date ')

call nc_close_file(ncid, source)

! Start adding to the sequence
call get_time(obs_time, osec, oday)

do ilat = 1, nlat
   do ilon = 1, nlon
      if (sst(ilon, ilat) == missing_sst .or. & 
          sst(ilon, ilat) /= sst(ilon, ilat)) cycle 
      if (iqc(ilon, ilat) <= OBS_QC_LOW_QUALITY) cycle  

      if (debug) then 
         write(*, '(4X, 4(A, F10.6))') 'lon: ', lon(ilon), ', lat: ', lat(ilat), & 
                  ', sst: ', sst(ilon, ilat), ', QC: ', iqc(ilon, ilat)
      endif   

      call create_3d_obs(lat(ilat), lon(ilon), 0.0_r8, VERTISSURFACE, &
                         sst(ilon, ilat)-t_kelvin, SATELLITE_BLENDED_SST,      &  
                         err(ilon, ilat), oday, osec, obs_qc, obs)
      call add_obs_to_seq(obs_seq, obs, obs_time, prev_obs, prev_time, first_obs)
   enddo
enddo

end subroutine collect_data

!------------------------------------------------------
! Get the reference time. This could be using 
! NF90_STRING (and not NF90_CHAR). If that's the 
! case, just fall back to what we expect the reference 
! to be. 
subroutine get_time_units(ncid, datestr)

integer, intent(in)           :: ncid
character(len=*), intent(out) :: datestr

integer :: status, varid, xtype, attlen

! fallback default
datestr = 'seconds since 1970-01-01 00:00:00'

status = nf90_inq_varid(ncid, 'time', varid)
if (status /= nf90_noerr) return

status = nf90_inquire_attribute(ncid, varid, 'units', xtype=xtype, len=attlen)
if (status /= nf90_noerr) return

if (xtype == nf90_char) then
   status = nf90_get_att(ncid, varid, 'units', datestr)
   if (status /= nf90_noerr) then
      ! keep fallback
   endif
endif

end subroutine get_time_units


!------------------------------------------------------------
! All collected obs (if any) are written in the seq file. 
subroutine finish_obs()

integer :: obs_num        

obs_num = get_num_obs(obs_seq)

if (obs_num > 0) then 
   if (debug) write(*, '(/, A, i0, A)') '>>>> Ready to write ', obs_num, ' observations:'
   
   call write_obs_seq(obs_seq, file_out)
   call destroy_obs(obs)
   call destroy_obs_sequence(obs_seq)
else 
   string1 = 'No obs were converted.'
   call error_handler(E_MSG, source, string1)
endif

call error_handler(E_MSG, source, 'Finished successfully.')

end subroutine finish_obs

!------------------------------------------------------------
! Clear up memory 
subroutine cleanup()

deallocate(lon, lat, sst, iqc, err)

end subroutine cleanup

end program cmems_sst_to_obs
