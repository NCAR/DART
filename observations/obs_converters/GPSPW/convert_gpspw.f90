! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_gpspw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_gpspw - program that reads ground-based GPS 
!                   precipitable water from Suominet sites
!                   in netCDF and writes a DART
!                   obs_seq file using the DART library routines.
!
!     created Nov. 2014 Soyoung Ha, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8
use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              find_namelist_in_file, check_namelist_read 
use  netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file, &
                                  nc_get_global_attribute
use  time_manager_mod, only : time_type, set_calendar_type, set_date, set_time, &
                              increment_time, get_time, get_date, operator(-), GREGORIAN
use      location_mod, only : VERTISUNDEF 
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use      obs_kind_mod, only : GPS_PRECIPITABLE_WATER
use obs_utilities_mod, only : getvar_real, get_or_fill_QC, add_obs_to_seq, &
                              create_3d_obs, getvar_int, getdimlen, getvar_real_2d, &
                              getvar_int_2d, query_varname, set_missing_name

implicit none

character(len=20),  parameter :: gpspw_netcdf_file = 'gpspw_input.nc'
character(len=129), parameter :: gpspw_out_file    = 'obs_seq.gpspw'
character(len=129)            :: gpspw_outfile

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

integer  :: iunit, io
integer  :: ncid, nstn, n, i, oday, osec, nused, index, ntime, it
logical  :: file_exist, first_obs
real(r8) :: qc
real(r8) :: pwv_miss = -999.
!real(r8) :: pwverr_miss = -999.

!character(len=129), allocatable :: stationID(:)
real(r8), allocatable :: lat(:), lon(:), elev(:), toff(:)
!real(r8), allocatable :: latu(:), lonu(:), levu(:), tobu(:)
real(r8), allocatable :: pwv(:,:), pwv_err(:,:)
!real(r8), allocatable :: wdel(:), tdel(:), rain(:)
!real(r8), allocatable :: wdir(:), wspd(:), rh(:), psfc(:), tsfc(:)
real(r8)              :: obs_window

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

integer, parameter :: MAX_NAME = 256
character(len=MAX_NAME) :: varname(5)

! For the data resource, check Readme in data/ directory.
! Sumoinet data provides observation times only as offset [sec] from 00Z 
! at each gregorian day. Thus, the actual observation date is read from 
! global attributes available in the file.
! The Suominet data has two different data regions - either over CONUS or the
! whole globe. We add the region name to the output file name (as well as
! the observation time).
! In this converter, we can process either hourly conus data or daily global 
! data depending on the global_data flag.
! There are other surface variables available at the stations if one wants
! to read from this dataset. For now, we only get pwv.

integer  :: max_num_obs                 = 1000000
real(r8) :: obs_window_hr               = 6          ! read data only every 6 hr from 00Z.
logical  :: global_data                 = .true.

namelist /convert_gpspw_nml/ global_data, max_num_obs, obs_window_hr
!                             include_wet_delay, include_total_delay,                 &
!                             include_surface_pressure, include_relative_humidity,    &
!                             include_surface_temperature, include_rain_accumulation, &
!                             include_wind_speed, include_wind_direction

character(len=6) :: ftail 
character(len=19):: sdate
character(len= 8):: ymd            ! YYYYMMDD
character(len=10):: ymdh           ! YYYYMMDDHH
integer          :: iyear, iday, ihour
integer          :: iyr, imo, idy, ihr, imn, isc

!------------
! start of executable code
!------------

call initialize_utilities('convert_gpspw')

call find_namelist_in_file("input.nml", "convert_gpspw_nml", iunit)
read(iunit, nml = convert_gpspw_nml, iostat = io)
call check_namelist_read(iunit, io, "mpas_obs_preproc_nml")


first_obs = .true.

ncid = nc_open_file_readonly(gpspw_netcdf_file, 'convert_gpspw')
call nc_get_global_attribute(ncid, 'start_date', sdate)
read(sdate,'(i4,1x,i3,1x,i2)')  iyear, iday, ihour

call set_calendar_type(GREGORIAN)
comp_day0 = set_date(iyear,1,1,0,0,0)
call get_time(comp_day0, osec, oday)
!print *,'DATE for ',iyear,' 01-01_00:00 => (oday, osec): ',oday, osec
oday = oday + iday - 1
comp_day0 = set_time(ihour*3600,oday)
call get_date(comp_day0, iyr, imo, idy, ihr, imn, isc)

! Final output file name
if( global_data ) then
   ftail = '.globe'
   write(ymd,'(I4,2I2.2)') iyr, imo, idy
   gpspw_outfile = trim(gpspw_out_file) // ftail // '.' // trim(ymd)
   print *,'OBS_DATE: ',oday,' gregorian days which is ',ymd
else
   ftail = '.conus'
   write(ymdh,'(I4,3I2.2)') iyr, imo, idy, ihr
   gpspw_outfile = trim(gpspw_out_file) // ftail // '.' // trim(ymdh)
   print *,'OBS_DATE: ',oday,' gregorian days which is ',ymdh
endif
print *,'Output file: ',trim(gpspw_outfile)


call getdimlen(ncid, "station",      nstn)
call getdimlen(ncid, "time_offset", ntime)
call set_missing_name("missing_value")

obs_window = obs_window_hr * 3600.0_r8        ! obs frequency (from hours to seconds)

allocate( lat(nstn))         
allocate( lon(nstn))
allocate(elev(nstn))         
allocate(toff(ntime))
allocate( pwv(ntime,nstn))         ;  allocate(pwv_err(ntime,nstn))   ! for multiple obs times

! Set the DART data quality control.  Be consistent with NCEP codes;
! 0 is 'must use', 1 is good, no reason not to use it.
qc = 0.0_r8

! read in the data arrays

! we have gpspw data files which have different names for the 
! lat/lon/elev/obs arrays in the netcdf file.  there doesn't seem
! to be a global attr to say which one is in use, so for now try
! both options.  

varname(1) = 'lat'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real(ncid, varname(index),  lat            ) ! station latitude

varname(1) = 'lon'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real(ncid, varname(index),  lon            ) ! station longitude

varname(1) = 'height'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real(ncid, varname(index),  elev           ) ! station elevation

varname(1) = 'pwv'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real_2d(ncid, varname(index),  pwv       )   ! precipitable water

varname(1) = 'pwv_err'
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real_2d(ncid, varname(index),  pwv_err   )   ! pwv obs error

varname(1) = 'time_offset'    ! "PWV window midpoint time delta from start_time"
call query_varname(ncid, 1, varname, index, force=.true.)
call    getvar_real(ncid, varname(index),  toff        ) ! obs time offset in seconds

! if user says to use them, read in QCs if present

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=gpspw_outfile, exist=file_exist)

if ( file_exist ) then

  ! existing file found, append to it
   call read_obs_seq(gpspw_outfile, 0, 0, nstn, obs_seq)

else

  ! create a new one
  call init_obs_sequence(obs_seq, num_copies, num_qc, max_num_obs)
  do i = 1, num_copies
    call set_copy_meta_data(obs_seq, i, 'Ground GPS PWV observation')
  end do
  do i = 1, num_qc
    call set_qc_meta_data(obs_seq, i, 'Data QC')
  end do

endif


timloop: do it = 1, ntime

  if(obs_window_hr > 0) then
     if(mod(toff(it),obs_window).ne.0) cycle timloop
  endif

  ! compute time of observation
  time_obs = increment_time(comp_day0, nint(toff(it)))

  ! extract actual time of observation in file into oday, osec.
  call get_time(time_obs, osec, oday)
  print*,'time_obs for it = ',it,oday, osec

  nused = 0
  stnloop: do n = 1, nstn

  ! check the lat/lon values to see if they are ok
  if ( lat(n) >  90.0_r8 .or. lat(n) <  -90.0_r8 ) cycle stnloop
  if ( lon(n) > 180.0_r8 .or. lon(n) < -180.0_r8 ) cycle stnloop

  if ( lon(n) < 0.0_r8 )  lon(n) = lon(n) + 360.0_r8

  ! add pwv data to obs_seq
  if ( pwv(it,n) /= pwv_miss .and. pwv(it,n) > 0.0_r8 ) then

   !  Ha: pwv observation error should be at least 1.5 mm. 
   if ( pwv_err(it,n) < 1.5 ) pwv_err(it,n) = 1.5 

      call create_3d_obs(lat(n), lon(n), elev(n), VERTISUNDEF, pwv(it,n), &
                         GPS_PRECIPITABLE_WATER, pwv_err(it,n), oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
      nused = nused + 1

  endif

  end do stnloop
  print*, 'nused =',nused

end do timloop

deallocate(lat)
deallocate(lon)
deallocate(elev)
deallocate(toff)
deallocate(pwv)         
deallocate(pwv_err)

call nc_close_file(ncid, 'convert_gpspw')

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, gpspw_outfile)

! end of main program
call finalize_utilities()


end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
