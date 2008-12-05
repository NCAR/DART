!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_madis_rawin - program that reads a netCDF file from the 
!                         MADIS database that contains rawinsonde data 
!                         and writes a text file that can be converted 
!                         into an obs_seq file.
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program convert_madis_rawin

use             types_mod, only : r8
use      time_manager_mod, only : time_type, set_calendar_type, set_date, &
                                  get_time, increment_time, GREGORIAN, operator(-)
use          location_mod, only : VERTISSURFACE, VERTISPRESSURE, VERTISHEIGHT
use      obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                                  static_init_obs_sequence, init_obs, write_obs_seq, &
                                  append_obs_to_seq, init_obs_sequence, get_num_obs, &
                                  set_copy_meta_data, set_qc_meta_data
use            meteor_mod, only : sat_vapor_pressure, specific_humidity, & 
                                  wind_dirspd_to_uv
use      ncep_obs_err_mod, only : ncep_rawin_temp_error, ncep_rawin_wind_error, &
                                  ncep_rawin_sfcpres_error, ncep_rawin_moist_error
use obs_def_altimeter_mod, only : compute_altimeter
use          obs_kind_mod, only : RADIOSONDE_U_WIND_COMPONENT,  & 
                                  RADIOSONDE_V_WIND_COMPONENT,  & 
                                  RADIOSONDE_TEMPERATURE,       & 
                                  RADIOSONDE_SPECIFIC_HUMIDITY, &
                                  RADIOSONDE_SURFACE_ALTIMETER 
use                netcdf

implicit none

character(len=19),  parameter :: rawin_in_file  = 'rawinsonde_input.nc'
character(len=129), parameter :: rawin_out_file = 'obs_seq.rawin'

integer, parameter ::   num_copies = 1,   &   ! number of copies in sequence
                        num_qc     = 1        ! number of QC entries

character (len=129) :: meta_data
character(len=80)   :: name
character(len=19)   :: datestr

integer :: oday, osec, iyear, imonth, iday, ihour, imin, isec, nman, nsig, nsound, & 
           nmaxml, nmaxsw, nmaxst, maxobs, k, n, fid, var_id, dsec, dday

integer, allocatable :: obscnt(:)

logical :: fexist, sigwnd, sigtmp

real(r8) :: obswindow, otime, lat, lon, elev, uwnd, vwnd, qobs, qsat, oerr, &
            pres_miss, wdir_miss, wspd_miss, tair_miss, tdew_miss, prespa, & 
            time_miss, qc, altim

real(r8), allocatable :: pres(:), wdir(:), wspd(:), tair(:), tdew(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs
type(time_type)         :: comp_day0, time_anal, time_obs

print*,'Enter the analysis time (yyyy-mm-dd_hh:mm:ss)'
read*,datestr
print*,'Enter the observation window (hours)'
read*,obswindow
print*,'Include significant level winds, temperature?'
read*,sigwnd, sigtmp

inquire(file=trim(rawin_in_file), exist=fexist) ! check for drop file 
if ( .NOT. fexist ) then
  print*,'Rawinsonde file ',rawin_in_file,' does not exist, exiting'
  stop
endif

call set_calendar_type(GREGORIAN)
read(datestr(1:4),   fmt='(i4)') iyear
read(datestr(6:7),   fmt='(i2)') imonth
read(datestr(9:10),  fmt='(i2)') iday
read(datestr(12:13), fmt='(i2)') ihour
read(datestr(15:16), fmt='(i2)') imin
read(datestr(18:19), fmt='(i2)') isec
time_anal = set_date(iyear, imonth, iday, ihour, imin, isec)
call get_time(time_anal, osec, oday)
obswindow = nint(obswindow * 3600.0_r8)
comp_day0 = set_date(1970, 1, 1, 0, 0, 0)

call check( nf90_open(rawin_in_file, nf90_nowrite, fid) )

call check( nf90_inq_dimid(fid, "recNum", var_id) )
call check( nf90_inquire_dimension(fid, var_id, name, nsound) )

allocate(obscnt(nsound))  
nmaxml = 0  ;  nmaxsw = 0  ;  nmaxst = 0
call check( nf90_inq_varid(fid, 'numMand', var_id) )
call check( nf90_get_var(fid, var_id, obscnt) )
do n = 1, nsound 
  if ( obscnt(n) > nmaxml .and. obscnt(n) < 25 )  nmaxml = obscnt(n)
end do

if ( sigwnd ) then
  call check( nf90_inq_varid(fid, 'numSigW', var_id) )
  call check( nf90_get_var(fid, var_id, obscnt) )
  do n = 1, nsound
    if ( obscnt(n) > nmaxsw .and. obscnt(n) < 150 )  nmaxsw = obscnt(n)
  end do
end if

if ( sigtmp ) then
  call check( nf90_inq_varid(fid, 'numSigT', var_id) )
  call check( nf90_get_var(fid, var_id, obscnt) )
  do n = 1, nsound
    if ( obscnt(n) > nmaxst .and. obscnt(n) < 150 )  nmaxst = obscnt(n)
  end do
end if

maxobs = nsound * (4 * nmaxml + 2 * nmaxsw + 2 * nmaxst + 1)
deallocate(obscnt)

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
inquire(file=rawin_out_file, exist=fexist)
if ( fexist ) then

  call read_obs_seq(rawin_out_file, 0, 0, maxobs, obs_seq)

else

  call init_obs_sequence(obs_seq, num_copies, num_qc, maxobs)
  do n = 1, num_copies
    meta_data = 'NCEP BUFR observation'
    call set_copy_meta_data(obs_seq, n, meta_data)
  end do
  do n = 1, num_qc
    meta_data = 'NCEP QC index'
    call set_qc_meta_data(obs_seq, n, meta_data)
  end do

end if

sondeloop : do n = 1, nsound !  loop over all soundings in the file 

  call check( nf90_inq_varid(fid, 'staLat', var_id) )
  call check( nf90_get_var(fid, var_id, lat,   start = (/ n /)) )
  call check( nf90_inq_varid(fid, 'staLon', var_id) )
  call check( nf90_get_var(fid, var_id, lon,   start = (/ n /)) )
  call check( nf90_inq_varid(fid, 'staElev', var_id) )
  call check( nf90_get_var(fid, var_id, elev,  start = (/ n /)) )
  call check( nf90_inq_varid(fid, 'numMand', var_id) )
  call check( nf90_get_var(fid, var_id, nman,  start = (/ n /)) ) 
  call check( nf90_inq_varid(fid, 'synTime', var_id) )
  call check( nf90_get_var(fid, var_id, otime, start = (/ n /)) )
  !call check( nf90_get_att(fid, var_id, '_FillValue', time_miss) )

  if ( otime < 0.0_r8 ) cycle sondeloop
  time_obs = increment_time(comp_day0,nint(mod(otime,86400.0_r8)),floor(otime/86400.0_r8))
  call get_time((time_anal - time_obs), dsec, dday)
  if ( dsec > obswindow .or. dday > 0 ) cycle sondeloop
  if ( nman <= 0 .or. nman > nmaxml   ) cycle sondeloop
  if ( lon < 0.0_r8 ) lon = lon + 360.0_r8
  qc = 1.0_r8

  allocate(pres(nman))  ;  allocate(tair(nman))  ;  allocate(tdew(nman))
  allocate(wdir(nman))  ;  allocate(wspd(nman))

  call check( nf90_inq_varid(fid, 'prMan', var_id) )
  call check( nf90_get_var(fid,var_id,pres,start=(/ 1, n /),count=(/ nman, 1 /)) )
  call check( nf90_get_att(fid, var_id, '_FillValue', pres_miss) )
  call check( nf90_inq_varid(fid, 'tpMan', var_id) )
  call check( nf90_get_var(fid,var_id,tair,start=(/ 1, n /),count=(/ nman, 1 /)) )
  call check( nf90_get_att(fid, var_id, '_FillValue', tair_miss) )
  call check( nf90_inq_varid(fid, 'tdMan', var_id) )
  call check( nf90_get_var(fid,var_id,tdew,start=(/ 1, n /),count=(/ nman, 1 /)) )
  call check( nf90_get_att(fid, var_id, '_FillValue', tdew_miss) )
  call check( nf90_inq_varid(fid, 'wdMan', var_id) )
  call check( nf90_get_var(fid,var_id,wdir,start=(/ 1, n /),count=(/ nman, 1 /)) )
  call check( nf90_get_att(fid, var_id, '_FillValue', wdir_miss) )
  call check( nf90_inq_varid(fid, 'wsMan', var_id) )
  call check( nf90_get_var(fid,var_id,wspd,start=(/ 1, n /),count=(/ nman, 1 /)) )
  call check( nf90_get_att(fid, var_id, '_FillValue', wspd_miss) )

  if ( pres(1) /= pres_miss ) then

    altim = compute_altimeter(pres(1), elev)
    call create_obs_type(lat, lon, elev, VERTISSURFACE, altim, &
                         RADIOSONDE_SURFACE_ALTIMETER, ncep_rawin_sfcpres_error, & 
                         oday, osec, qc, obs)
    call append_obs_to_seq(obs_seq, obs)

  end if

  do k = 2, nman   ! obtain the mandatory level data

    prespa = pres(k) * 100.0_r8

    if ( wdir(k) /= wdir_miss .and. wspd(k) /= wspd_miss ) then

      call wind_dirspd_to_uv(wdir(k), wspd(k), uwnd, vwnd)
      oerr = ncep_rawin_wind_error(pres(k))
      call create_obs_type(lat, lon, prespa, VERTISPRESSURE, uwnd, &
                           RADIOSONDE_U_WIND_COMPONENT, oerr, oday, osec, qc, obs)
      call append_obs_to_seq(obs_seq, obs)
      call create_obs_type(lat, lon, prespa, VERTISPRESSURE, vwnd, &
                           RADIOSONDE_V_WIND_COMPONENT, oerr, oday, osec, qc, obs)
      call append_obs_to_seq(obs_seq, obs)

    end if

    if ( tair(k) /= tair_miss ) then

      oerr = ncep_rawin_temp_error(pres(k))
      call create_obs_type(lat, lon, prespa, VERTISPRESSURE, tair(k), &
                           RADIOSONDE_TEMPERATURE, oerr, oday, osec, qc, obs)
      call append_obs_to_seq(obs_seq, obs)

    end if

    if ( tair(k) /= tair_miss .and. tdew(k) /= tdew_miss ) then

      qobs = tair(k) - tdew(k)
      qobs = specific_humidity(sat_vapor_pressure(qobs),    prespa)
      qsat = specific_humidity(sat_vapor_pressure(tair(k)), prespa)
      oerr = max(ncep_rawin_moist_error * qsat, 0.0001_r8)
      if ( qobs > 0.0_r8 .and. qobs < 0.070_r8 ) then
        call create_obs_type(lat, lon, prespa, VERTISPRESSURE, qobs, &
                             RADIOSONDE_SPECIFIC_HUMIDITY, oerr, oday, osec, qc, obs)
        call append_obs_to_seq(obs_seq, obs)
      end if

    end if

  end do
  deallocate(pres, wdir, wspd, tair, tdew)

  !  If desired, read the significant-level temperature data, write to obs. seq.
  call check( nf90_inq_varid(fid, 'numSigT', var_id) )
  call check( nf90_get_var(fid, var_id, nsig,  start = (/ n /)) )

  if ( sigtmp .and. nsig <= nmaxst ) then

    allocate(pres(nsig))  ;  allocate(tair(nsig))  ;  allocate(tdew(nsig))

    !  read significant level data
    call check( nf90_inq_varid(fid, 'prSigT', var_id) )
    call check( nf90_get_var(fid,var_id,pres,start=(/ 1, n /),count=(/ nsig, 1 /)) )
    call check( nf90_get_att(fid, var_id, '_FillValue', pres_miss) )
    call check( nf90_inq_varid(fid, 'tpSigT', var_id) )
    call check( nf90_get_var(fid,var_id,tair,start=(/ 1, n /),count=(/ nsig, 1 /)) )
    call check( nf90_get_att(fid, var_id, '_FillValue', tair_miss) )
    call check( nf90_inq_varid(fid, 'tdSigT', var_id) )
    call check( nf90_get_var(fid,var_id,tdew,start=(/ 1, n /),count=(/ nsig, 1 /)) )
    call check( nf90_get_att(fid, var_id, '_FillValue', tdew_miss) )

    do k = 1, nsig

      prespa = pres(k) * 100.0_r8

      if ( tair(k) /= tair_miss ) then

        oerr = ncep_rawin_temp_error(pres(k))
        call create_obs_type(lat, lon, prespa, VERTISPRESSURE, tair(k), &
                             RADIOSONDE_TEMPERATURE, oerr, oday, osec, qc, obs)
        call append_obs_to_seq(obs_seq, obs)

      end if

      if ( tair(k) /= tair_miss .and. tdew(k) /= tdew_miss ) then

        qobs = tair(k) - tdew(k)
        qobs = specific_humidity(sat_vapor_pressure(qobs),    prespa)
        qsat = specific_humidity(sat_vapor_pressure(tair(k)), prespa)
        oerr = max(ncep_rawin_moist_error * qsat, 0.0001_r8)
        if ( qobs > 0.0_r8 .and. qobs < 0.070_r8 ) then
          call create_obs_type(lat, lon, prespa, VERTISPRESSURE, qobs, &
                               RADIOSONDE_SPECIFIC_HUMIDITY, oerr, oday, osec, qc, obs)
          call append_obs_to_seq(obs_seq, obs)
        end if

      endif

    end do
    deallocate(pres, tair, tdew)

  end if

  !  If desired, read the significant-level wind data, write to obs. seq.
  call check( nf90_inq_varid(fid, 'numSigW', var_id) )
  call check( nf90_get_var(fid, var_id, nsig,  start = (/ n /)) )
  
  if ( sigwnd .and. nsig <= nmaxsw ) then

    allocate(pres(nsig))  ;  allocate(wdir(nsig))  ;  allocate(wspd(nsig))

    !  read significant level data
    call check( nf90_inq_varid(fid, 'htSigW', var_id) )
    call check( nf90_get_var(fid,var_id,pres,start=(/ 1, n /),count=(/ nsig, 1 /)) )
    call check( nf90_get_att(fid, var_id, '_FillValue', pres_miss) )
    call check( nf90_inq_varid(fid, 'wdSigW', var_id) )
    call check( nf90_get_var(fid,var_id,wdir,start=(/ 1, n /),count=(/ nsig, 1 /)) )
    call check( nf90_get_att(fid, var_id, '_FillValue', wdir_miss) )
    call check( nf90_inq_varid(fid, 'wsSigW', var_id) )
    call check( nf90_get_var(fid,var_id,wspd,start=(/ 1, n /),count=(/ nsig, 1 /)) )
    call check( nf90_get_att(fid, var_id, '_FillValue', wspd_miss) )

    do k = 1, nsig

      !  add data to the observation sequence here.
      if ( wdir(k) /= wdir_miss .and. wspd(k) /= wspd_miss ) then

        call wind_dirspd_to_uv(wdir(k), wspd(k), uwnd, vwnd)
        oerr = ncep_rawin_wind_error(500.0_r8)
        call create_obs_type(lat, lon, pres(k), VERTISHEIGHT, uwnd, &
                             RADIOSONDE_U_WIND_COMPONENT, oerr, oday, osec, qc, obs)
        call append_obs_to_seq(obs_seq, obs)
        call create_obs_type(lat, lon, pres(k), VERTISHEIGHT, vwnd, &
                             RADIOSONDE_V_WIND_COMPONENT, oerr, oday, osec, qc, obs)
        call append_obs_to_seq(obs_seq, obs)

      end if

    end do
    deallocate(pres, wdir, wspd)

  end if

enddo sondeloop

call check( nf90_close(fid) )
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, rawin_out_file)

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

if(istatus /= nf90_noerr) then 
  print*,'Netcdf error: ',trim(nf90_strerror(istatus))
  stop
end if

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
use     location_mod, only : location_type, set_location
use time_manager_mod, only : time_type, set_time

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
