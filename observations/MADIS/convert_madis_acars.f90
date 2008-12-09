!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_madis_acars - program that reads a MADIS netCDF ACARS 
!                         observation file and writes a text file of 
!                         observations within the analysis time.  The 
!                         text file can be used in other programs that 
!                         write obs_seq.out files.
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program convert_madis_acars

use        types_mod, only : r8
use       meteor_mod, only : pres_alt_to_pres, sat_vapor_pressure, & 
                             specific_humidity, wind_dirspd_to_uv
use ncep_obs_err_mod, only : ncep_acars_wind_error, ncep_acars_temp_error, &
                             ncep_acars_moist_error
use     location_mod, only : VERTISPRESSURE
use obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                             static_init_obs_sequence, init_obs, write_obs_seq, &
                             append_obs_to_seq, init_obs_sequence, get_num_obs, &
                             set_copy_meta_data, set_qc_meta_data
use time_manager_mod, only : time_type, set_calendar_type, set_date, &
                             increment_time, get_time, GREGORIAN, operator(-)
use     obs_kind_mod, only : ACARS_U_WIND_COMPONENT, ACARS_V_WIND_COMPONENT, &
                             ACARS_TEMPERATURE, ACARS_SPECIFIC_HUMIDITY
use           netcdf

implicit none

character(len=14),  parameter :: acars_netcdf_file = 'acars_input.nc'
character(len=129), parameter :: acars_out_file    = 'obs_seq.acars'

integer, parameter ::   num_copies = 1,   &   ! number of copies in sequence
                        num_qc     = 1        ! number of QC entries

character (len=129) :: meta_data
character (len=80)  :: name
character (len=19)  :: datime
integer :: rcode, ncid, varid, nobs, n, i, window_sec, dday, dsec, &
           oday, osec, nused, iyear, imonth, iday, ihour, imin, isec
logical :: file_exist
real(r8) :: palt_miss, tair_miss, relh_miss, wdir_miss, wspd_miss, uwnd, &
            vwnd, qobs, qsat, oerr, window_hours, pres, qc

integer,  allocatable :: tobs(:)
real(r8), allocatable :: lat(:), lon(:), palt(:), tair(:), relh(:), &
                         wdir(:), wspd(:), latu(:), lonu(:), palu(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs
type(time_type)         :: comp_day0, time_obs, time_anal

print*,'Enter the target assimilation time (yyyy-mm-dd_hh:mm:ss)'
read*,datime
print*,'Enter the observation window (hours)'
read*,window_hours

call set_calendar_type(GREGORIAN)
read(datime(1:4),   fmt='(i4)') iyear
read(datime(6:7),   fmt='(i2)') imonth
read(datime(9:10),  fmt='(i2)') iday
read(datime(12:13), fmt='(i2)') ihour
read(datime(15:16), fmt='(i2)') imin
read(datime(18:19), fmt='(i2)') isec
time_anal = set_date(iyear, imonth, iday, ihour, imin, isec)
call get_time(time_anal, osec, oday)
window_sec = nint(window_hours * 3600.0)
comp_day0 = set_date(1970, 1, 1, 0, 0, 0)

rcode = nf90_open(acars_netcdf_file, nf90_nowrite, ncid)

call check( nf90_inq_dimid(ncid, "recNum", varid) )
call check( nf90_inquire_dimension(ncid, varid, name, nobs) )

allocate( lat(nobs))  ;  allocate( lon(nobs))
allocate(palt(nobs))  ;  allocate(tobs(nobs))
allocate(tair(nobs))  ;  allocate(relh(nobs))
allocate(wdir(nobs))  ;  allocate(wspd(nobs))
allocate(latu(nobs))  ;  allocate(lonu(nobs))
allocate(palu(nobs))

! read the latitude array
call check( nf90_inq_varid(ncid, "latitude", varid) )
call check( nf90_get_var(ncid, varid, lat) )

! read the latitude array
call check( nf90_inq_varid(ncid, "longitude", varid) )
call check( nf90_get_var(ncid, varid, lon) )

! read the pressure altitude array
call check( nf90_inq_varid(ncid, "altitude", varid) )
call check( nf90_get_var(ncid, varid, palt) )
call check( nf90_get_att(ncid, varid, '_FillValue', palt_miss) )

! read the air temperature array
call check( nf90_inq_varid(ncid, "temperature", varid) )
call check( nf90_get_var(ncid, varid, tair) )
call check( nf90_get_att(ncid, varid, '_FillValue', tair_miss) )

! read the dew-point temperature array
call check( nf90_inq_varid(ncid, "downlinkedRH", varid) )
call check( nf90_get_var(ncid, varid, relh) )
call check( nf90_get_att(ncid, varid, '_FillValue', relh_miss) )

! read the wind direction array
call check( nf90_inq_varid(ncid, "windDir", varid) )
call check( nf90_get_var(ncid, varid, wdir) )
call check( nf90_get_att(ncid, varid, '_FillValue', wdir_miss) )

! read the wind speed array
call check( nf90_inq_varid(ncid, "windSpeed", varid) )
call check( nf90_get_var(ncid, varid, wspd) )
call check( nf90_get_att(ncid, varid, '_FillValue', wspd_miss) )

! read the observation time array
call check( nf90_inq_varid(ncid, "timeObs", varid) )
call check( nf90_get_var(ncid, varid, tobs) )

call check( nf90_close(ncid) )

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
inquire(file=acars_out_file, exist=file_exist)
if ( file_exist ) then

  call read_obs_seq(acars_out_file, 0, 0, 4*nobs, obs_seq)

else

  call init_obs_sequence(obs_seq, num_copies, num_qc, 4*nobs)
  do i = 1, num_copies
    meta_data = 'NCEP BUFR observation'
    call set_copy_meta_data(obs_seq, i, meta_data)
  end do
  do i = 1, num_qc
    meta_data = 'NCEP QC index'
    call set_qc_meta_data(obs_seq, i, meta_data)
  end do

end if

nused = 0
obsloop: do n = 1, nobs

  ! determine if the observation is within the window
  time_obs = increment_time(comp_day0, mod(tobs(n),86400), tobs(n) / 86400)
  call get_time((time_anal - time_obs), dsec, dday)
  if ( dsec > window_sec .or. dday > 0 ) cycle obsloop
  if ( lon(n) < 0.0_r8 )  lon(n) = lon(n) + 360.0_r8

  ! make sure this observation is not a repeat
  do i = 1, nused
    if ( lon(n)  == lonu(i) .and. lat(n) == latu(i) .and. & 
         palt(n) == palu(i) ) cycle obsloop
  end do
  qc = 1.0_r8
 
  if ( palt(n) == palt_miss ) cycle obsloop  
  pres = pres_alt_to_pres(palt(n))

  ! add wind component data to obs. sequence
  if ( wdir(n) /= wdir_miss .and. wspd(n) /= wspd_miss ) then

    call wind_dirspd_to_uv(wdir(n), wspd(n), uwnd, vwnd)
    if ( abs(uwnd) < 150.0_r8 .and. abs(vwnd) < 150.0_r8 ) then

      oerr = ncep_acars_wind_error
      call create_obs_type(lat(n), lon(n), pres, VERTISPRESSURE, uwnd, &
                           ACARS_U_WIND_COMPONENT, oerr, oday, osec, qc, obs)
      call append_obs_to_seq(obs_seq, obs)
      call create_obs_type(lat(n), lon(n), pres, VERTISPRESSURE, vwnd, &
                           ACARS_V_WIND_COMPONENT, oerr, oday, osec, qc, obs)
      call append_obs_to_seq(obs_seq, obs)

    end if

  end if

  ! add air temperature data to obs. sequence
  if ( tair(n) /= tair_miss ) then 
   
    oerr = ncep_acars_temp_error(pres * 0.01_r8) 
    call create_obs_type(lat(n), lon(n), pres, VERTISPRESSURE, tair(n), &
                         ACARS_TEMPERATURE, oerr, oday, osec, qc, obs)
    call append_obs_to_seq(obs_seq, obs)

  end if

  ! add relative humidity data to obs. sequence, but as specific humidity
  if ( tair(n) /= tair_miss .and. relh(n) /= relh_miss ) then

    qsat = specific_humidity(sat_vapor_pressure(tair(n)), pres)
    qobs = qsat * relh(n)
    oerr = max(ncep_acars_moist_error * qsat, 0.0001_r8)

    if ( abs(qobs) < 0.1_r8 ) then

      call create_obs_type(lat(n), lon(n), pres, VERTISPRESSURE, qobs, &
                         ACARS_SPECIFIC_HUMIDITY, oerr, oday, osec, qc, obs)
      call append_obs_to_seq(obs_seq, obs)

    end if

  end if

  nused = nused + 1
  latu(nused) = lat(n)
  lonu(nused) = lon(n)
  palu(nused) = palt(n)

end do obsloop

if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, acars_out_file)

stop
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   check - subroutine that checks the flag from a netCDF function.  If
!           there is an error, the message is displayed.
!
!    istatus - netCDF output flag
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check( istatus ) 

use netcdf

implicit none

integer, intent (in) :: istatus

if(istatus /= nf90_noerr) print*,'Netcdf error: ',trim(nf90_strerror(istatus))

end subroutine check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   create_obs_type - subroutine that is used to create an observation
!                     type from observation data.
!
!    lat   - latitude of observation
!    lon   - longitude of observation
!    pres  - pressure of observation
!    vcord - vertical coordinate
!    obsv  - observation value
!    okind - observation kind
!    oerr  - observation error
!    day   - gregorian day
!    sec   - gregorian second
!    qc    - quality control value
!    obs   - observation type
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_obs_type(lat, lon, pres, vcord, obsv, okind, oerr, day, sec, qc, obs)

use types_mod,        only : r8
use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location
use     location_mod, only : set_location
use time_manager_mod, only : set_time

implicit none

integer, intent(in)         :: okind, vcord, day, sec
real(r8), intent(in)        :: lat, lon, pres, obsv, oerr, qc
type(obs_type), intent(inout) :: obs

real(r8)              :: obs_val(1), qc_val(1)
type(obs_def_type)    :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, pres, vcord))
call set_obs_def_kind(obs_def, okind)
call set_obs_def_time(obs_def, set_time(sec, day))
call set_obs_def_error_variance(obs_def, oerr * oerr)
call set_obs_def(obs, obs_def)

obs_val(1) = obsv
call set_obs_values(obs, obs_val)
qc_val(1)  = qc
call set_qc(obs, qc_val)

return
end subroutine create_obs_type
