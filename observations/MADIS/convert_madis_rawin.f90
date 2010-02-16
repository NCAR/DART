! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program convert_madis_rawin

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_madis_rawin - program that reads a netCDF file from the 
!                         MADIS database that contains rawinsonde data 
!                         and writes a DART obs_seq file using the DART 
!                         library routines.
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM
!     modified Dec. 2008 Soyoung Ha and David Dowell, NCAR/MMM
!     - added dewpoint as an output variable
!     - added relative humidity as an output variable
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use             types_mod, only : r8, missing_r8
use      time_manager_mod, only : time_type, set_calendar_type, set_date, &
                                  get_time, increment_time, GREGORIAN, operator(-)
use          location_mod, only : VERTISSURFACE, VERTISPRESSURE, VERTISHEIGHT
use      obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                                  static_init_obs_sequence, init_obs, write_obs_seq, &
                                  append_obs_to_seq, init_obs_sequence, get_num_obs, &
                                  set_copy_meta_data, set_qc_meta_data
use            meteor_mod, only : sat_vapor_pressure, specific_humidity, & 
                                  wind_dirspd_to_uv, pres_alt_to_pres, &
                                  temp_and_dewpoint_to_rh
use           obs_err_mod, only : rawin_temp_error, rawin_wind_error, &
                                  rawin_pres_error, rawin_rel_hum_error
use  dewpoint_obs_err_mod, only : dewpt_error_from_rh_and_temp, &
                                  rh_error_from_dewpt_and_temp
use obs_def_altimeter_mod, only : compute_altimeter
use          obs_kind_mod, only : RADIOSONDE_U_WIND_COMPONENT,      & 
                                  RADIOSONDE_V_WIND_COMPONENT,      & 
                                  RADIOSONDE_TEMPERATURE,           & 
                                  RADIOSONDE_SPECIFIC_HUMIDITY,     &
                                  RADIOSONDE_RELATIVE_HUMIDITY,     &
                                  RADIOSONDE_DEWPOINT, &
                                  RADIOSONDE_SURFACE_ALTIMETER 
use                netcdf

implicit none

character(len=19),  parameter :: rawin_in_file  = 'rawinsonde_input.nc'
character(len=129), parameter :: rawin_out_file = 'obs_seq.rawin'

! the following logical parameters control which water-vapor variables appear in the output file
! whether to use the NCEP error or Lin and Hubbard (2004) moisture error model, and if the
! input file has data quality control fields, whether to use or ignore them.
logical, parameter :: LH_err                    = .false.
logical, parameter :: include_specific_humidity = .true.
logical, parameter :: include_relative_humidity = .false.
logical, parameter :: include_dewpoint          = .false.
logical, parameter :: use_input_qc              = .true. 

integer, parameter ::   num_copies = 1,   &   ! number of copies in sequence
                        num_qc     = 1        ! number of QC entries

character (len=129) :: meta_data
character(len=80)   :: name
character(len=19)   :: datestr

integer :: oday, osec, iyear, imonth, iday, ihour, imin, isec, nman, nsig, nsound, & 
           nmaxml, nmaxsw, nmaxst, maxobs, nvars_man, nvars_sigt, k, n, ncid, varid, &
           dsec, dday, nfrc

integer, allocatable :: obscnt(:)

logical :: fexist, sigwnd, sigtmp, input_has_qc

real(r8) :: obswindow, otime, lat, lon, elev, uwnd, vwnd, qobs, qsat, dptk, oerr, &
            pres_miss, wdir_miss, wspd_miss, tair_miss, tdew_miss, prespa, & 
            time_miss, qc, altim, rh, qerr

real(r8), allocatable :: pres(:), wdir(:), wspd(:), tair(:), tdew(:)
integer,  allocatable :: qc_pres(:), qc_wdir(:), qc_wspd(:), qc_tair(:), qc_tdew(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs
type(time_type)         :: comp_day0, time_anal, time_obs

print*,'Enter the analysis time (yyyy-mm-dd_hh:mm:ss), window (hours): '
read*,datestr, obswindow
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

call check( nf90_open(rawin_in_file, nf90_nowrite, ncid) )

call check( nf90_inq_dimid(ncid, "recNum", varid) )
call check( nf90_inquire_dimension(ncid, varid, name, nsound) )

allocate(obscnt(nsound))  
nmaxml = 0  ;  nmaxsw = 0  ;  nmaxst = 0
call check( nf90_inq_varid(ncid, 'numMand', varid) )
call check( nf90_get_var(ncid, varid, obscnt) )
do n = 1, nsound 
  if ( obscnt(n) > nmaxml .and. obscnt(n) < 25 )  nmaxml = obscnt(n)
end do

if ( sigwnd ) then
  call check( nf90_inq_varid(ncid, 'numSigW', varid) )
  call check( nf90_get_var(ncid, varid, obscnt) )
  do n = 1, nsound
    if ( obscnt(n) > nmaxsw .and. obscnt(n) < 150 )  nmaxsw = obscnt(n)
  end do
end if

if ( sigtmp ) then
  call check( nf90_inq_varid(ncid, 'numSigT', varid) )
  call check( nf90_get_var(ncid, varid, obscnt) )
  do n = 1, nsound
    if ( obscnt(n) > nmaxst .and. obscnt(n) < 150 )  nmaxst = obscnt(n)
  end do
end if

nvars_man = 4
nvars_sigt = 1
if (include_specific_humidity) then
  nvars_man = nvars_man + 1
  nvars_sigt = nvars_sigt + 1
end if
if (include_relative_humidity) then
  nvars_man = nvars_man + 1
  nvars_sigt = nvars_sigt + 1
end if
if (include_dewpoint) then
  nvars_man = nvars_man + 1
  nvars_sigt = nvars_sigt + 1
end if

maxobs = nsound * (nvars_man * nmaxml + 2 * nmaxsw + nvars_sigt * nmaxst + 1)
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
    meta_data = 'MADIS observation'
    call set_copy_meta_data(obs_seq, n, meta_data)
  end do
  do n = 1, num_qc
    meta_data = 'Data QC'
    call set_qc_meta_data(obs_seq, n, meta_data)
  end do

end if

sondeloop : do n = 1, nsound !  loop over all soundings in the file 

  call check( nf90_inq_varid(ncid, 'staLat', varid) )
  call check( nf90_get_var(ncid, varid, lat,   start = (/ n /)) )
  call check( nf90_inq_varid(ncid, 'staLon', varid) )
  call check( nf90_get_var(ncid, varid, lon,   start = (/ n /)) )
  call check( nf90_inq_varid(ncid, 'staElev', varid) )
  call check( nf90_get_var(ncid, varid, elev,  start = (/ n /)) )
  call check( nf90_inq_varid(ncid, 'numMand', varid) )
  call check( nf90_get_var(ncid, varid, nman,  start = (/ n /)) ) 
  call check( nf90_inq_varid(ncid, 'synTime', varid) )
  call check( nf90_get_var(ncid, varid, otime, start = (/ n /)) )
  !call check( nf90_get_att(ncid, varid, '_FillValue', time_miss) )

  if ( otime < 0.0_r8 ) cycle sondeloop
  time_obs = increment_time(comp_day0,nint(mod(otime,86400.0_r8)),floor(otime/86400.0_r8))
  call get_time((time_anal - time_obs), dsec, dday)
  if ( dsec > obswindow .or. dday > 0 ) cycle sondeloop
  if ( nman <= 0 .or. nman > nmaxml   ) cycle sondeloop
  if ( lon < 0.0_r8 ) lon = lon + 360.0_r8
  qc = 1.0_r8

  allocate(pres(nman))  ;  allocate(tair(nman))  ;  allocate(tdew(nman))
  allocate(wdir(nman))  ;  allocate(wspd(nman))
  
  allocate(qc_pres(nman))  ;  allocate(qc_tair(nman))  ;  allocate(qc_tdew(nman))
  allocate(qc_wdir(nman))  ;  allocate(qc_wspd(nman))

  call check( nf90_inq_varid(ncid, 'prMan', varid) )
  call check( nf90_get_var(ncid,varid,pres,start=(/ 1, n /),count=(/ nman, 1 /)) )
  call check( nf90_get_att(ncid, varid, '_FillValue', pres_miss) )
  call check( nf90_inq_varid(ncid, 'tpMan', varid) )
  call check( nf90_get_var(ncid,varid,tair,start=(/ 1, n /),count=(/ nman, 1 /)) )
  call check( nf90_get_att(ncid, varid, '_FillValue', tair_miss) )
  call check( nf90_inq_varid(ncid, 'tdMan', varid) )
  call check( nf90_get_var(ncid,varid,tdew,start=(/ 1, n /),count=(/ nman, 1 /)) )
  call check( nf90_get_att(ncid, varid, '_FillValue', tdew_miss) )
  call check( nf90_inq_varid(ncid, 'wdMan', varid) )
  call check( nf90_get_var(ncid,varid,wdir,start=(/ 1, n /),count=(/ nman, 1 /)) )
  call check( nf90_get_att(ncid, varid, '_FillValue', wdir_miss) )
  call check( nf90_inq_varid(ncid, 'wsMan', varid) )
  call check( nf90_get_var(ncid,varid,wspd,start=(/ 1, n /),count=(/ nman, 1 /)) )
  call check( nf90_get_att(ncid, varid, '_FillValue', wspd_miss) )

  ! pick a random QC field and test for it.  if it's there, set
  ! the 'has qc' flag to true.  otherwise, set it to false.
  nfrc = nf90_inq_varid(ncid, "prManQCR", varid) 
  input_has_qc = (nfrc == nf90_noerr)

  ! read the QC check for each variable
  if (input_has_qc .and. use_input_qc) then
     call check( nf90_inq_varid(ncid, "prManQCR", varid) )
     call check( nf90_get_var(ncid, varid, qc_pres,start=(/ 1, n /),count=(/ nman, 1 /)) )
     call check( nf90_inq_varid(ncid, "tpManQCR", varid) )
     call check( nf90_get_var(ncid, varid, qc_tair,start=(/ 1, n /),count=(/ nman, 1 /)) )
     call check( nf90_inq_varid(ncid, "tdManQCR", varid) )
     call check( nf90_get_var(ncid, varid, qc_tdew,start=(/ 1, n /),count=(/ nman, 1 /)) )
     call check( nf90_inq_varid(ncid, "wdManQCR", varid) )
     call check( nf90_get_var(ncid, varid, qc_wdir,start=(/ 1, n /),count=(/ nman, 1 /)) )
     call check( nf90_inq_varid(ncid, "wsManQCR", varid) )
     call check( nf90_get_var(ncid, varid, qc_wspd,start=(/ 1, n /),count=(/ nman, 1 /)) )
  else
     ! if input contains no QCs, or user said skip them. assume all are ok.
     qc_pres = 0
     qc_tair = 0 ;  qc_tdew = 0
     qc_wdir = 0 ;  qc_wspd = 0
  endif

  if ( pres(1) /= pres_miss .and. qc_pres(1) == 0 ) then

    altim = compute_altimeter(pres(1), elev)
    oerr  = rawin_pres_error(pres_alt_to_pres(elev) * 0.01_r8)
    if ( altim >= 880.0_r8 .and. altim <= 1100.0_r8 .and. oerr /= missing_r8 ) then

      call create_obs_type(lat, lon, elev, VERTISSURFACE, altim, &
                           RADIOSONDE_SURFACE_ALTIMETER, oerr, oday, osec, qc, obs)
      call append_obs_to_seq(obs_seq, obs)

    end if

  end if

  do k = 2, nman   ! obtain the mandatory level data

    prespa = pres(k) * 100.0_r8

    if ( wdir(k) /= wdir_miss .and. wspd(k) /= wspd_miss .and. qc_wdir(k) == 0 .and. qc_wspd(k) == 0 ) then

      call wind_dirspd_to_uv(wdir(k), wspd(k), uwnd, vwnd)
      oerr = rawin_wind_error(pres(k))
      if ( abs(uwnd) <= 150.0_r8 .and. & 
           abs(vwnd) <= 150.0_r8 .and. oerr /= missing_r8 ) then

        call create_obs_type(lat, lon, prespa, VERTISPRESSURE, uwnd, &
                             RADIOSONDE_U_WIND_COMPONENT, oerr, oday, osec, qc, obs)
        call append_obs_to_seq(obs_seq, obs)
        call create_obs_type(lat, lon, prespa, VERTISPRESSURE, vwnd, &
                             RADIOSONDE_V_WIND_COMPONENT, oerr, oday, osec, qc, obs)
        call append_obs_to_seq(obs_seq, obs)

      end if

    end if

    if ( tair(k) /= tair_miss .and. qc_tair(k) == 0 ) then

      oerr = rawin_temp_error(pres(k))
      if ( tair(k) >= 180.0_r8 .and. tair(k) <= 330.0_r8 .and. oerr /= missing_r8 ) then

        call create_obs_type(lat, lon, prespa, VERTISPRESSURE, tair(k), &
                             RADIOSONDE_TEMPERATURE, oerr, oday, osec, qc, obs)
        call append_obs_to_seq(obs_seq, obs)

      end if

    end if

    if ( include_specific_humidity .and. tair(k) /= tair_miss .and. tdew(k) /= tdew_miss &
          .and. qc_tair(k) == 0 .and. qc_tdew(k) == 0 ) then

      ! tdew is the dewpoint depression
      dptk = tair(k) - tdew(k)
      qobs = specific_humidity(sat_vapor_pressure(dptk),    prespa)
      qsat = specific_humidity(sat_vapor_pressure(tair(k)), prespa)
      if (LH_err ) then
        qerr = rh_error_from_dewpt_and_temp(tair(k), dptk)
      else
        qerr = rawin_rel_hum_error(pres(k), tair(k), qobs / qsat)
      end if
      oerr = max(qerr * qsat, 0.0001_r8)

      if ( qobs > 0.0_r8 .and. qobs <= 0.070_r8 .and. qerr /= missing_r8 ) then
        call create_obs_type(lat, lon, prespa, VERTISPRESSURE, qobs, &
                             RADIOSONDE_SPECIFIC_HUMIDITY, oerr, oday, osec, qc, obs)
        call append_obs_to_seq(obs_seq, obs)
      end if

    end if

    if ( include_relative_humidity .and. tair(k) /= tair_miss .and. tdew(k) /= tdew_miss &
         .and. qc_tair(k) == 0 .and. qc_tdew(k) == 0 ) then

      ! tdew is the dewpoint depression
      dptk = tair(k) - tdew(k)
      rh = temp_and_dewpoint_to_rh(tair(k), dptk)
      if (LH_err ) then
        oerr = rh_error_from_dewpt_and_temp(tair(k), dptk)
      else
        oerr = rawin_rel_hum_error(pres(k), tair(k), rh)
      end if

      call create_obs_type(lat, lon, prespa, VERTISPRESSURE, rh, &
                           RADIOSONDE_RELATIVE_HUMIDITY, oerr, oday, osec, qc, obs)
      call append_obs_to_seq(obs_seq, obs)

    end if

    if ( include_dewpoint .and. tair(k) /= tair_miss .and. tdew(k) /= tdew_miss &
         .and. qc_tair(k) == 0 .and. qc_tdew(k) == 0 ) then

      ! tdew is the dewpoint depression
      dptk = tair(k) - tdew(k)
      rh = temp_and_dewpoint_to_rh(tair(k), dptk)
      oerr = dewpt_error_from_rh_and_temp(tair(k), rh)
      call create_obs_type(lat, lon, prespa, VERTISPRESSURE, dptk, &
                           RADIOSONDE_DEWPOINT, oerr, oday, osec, qc, obs)
      call append_obs_to_seq(obs_seq, obs)

    end if

  end do
  deallocate(pres, wdir, wspd, tair, tdew, qc_pres, qc_wdir, qc_wspd, qc_tair, qc_tdew)

  !  If desired, read the significant-level temperature data, write to obs. seq.
  call check( nf90_inq_varid(ncid, 'numSigT', varid) )
  call check( nf90_get_var(ncid, varid, nsig,  start = (/ n /)) )

  if ( sigtmp .and. nsig <= nmaxst ) then

    allocate(pres(nsig))  ;  allocate(tair(nsig))  ;  allocate(tdew(nsig))
    allocate(qc_pres(nsig))  ;  allocate(qc_tair(nsig))  ;  allocate(qc_tdew(nsig))

    !  read significant level data
    call check( nf90_inq_varid(ncid, 'prSigT', varid) )
    call check( nf90_get_var(ncid,varid,pres,start=(/ 1, n /),count=(/ nsig, 1 /)) )
    call check( nf90_get_att(ncid, varid, '_FillValue', pres_miss) )
    call check( nf90_inq_varid(ncid, 'tpSigT', varid) )
    call check( nf90_get_var(ncid,varid,tair,start=(/ 1, n /),count=(/ nsig, 1 /)) )
    call check( nf90_get_att(ncid, varid, '_FillValue', tair_miss) )
    call check( nf90_inq_varid(ncid, 'tdSigT', varid) )
    call check( nf90_get_var(ncid,varid,tdew,start=(/ 1, n /),count=(/ nsig, 1 /)) )
    call check( nf90_get_att(ncid, varid, '_FillValue', tdew_miss) )
    if (input_has_qc) then
       call check( nf90_inq_varid(ncid, "prSigTQCR", varid) )
       call check( nf90_get_var(ncid, varid, qc_pres,start=(/ 1, n /),count=(/ nsig, 1 /)) )
       call check( nf90_inq_varid(ncid, "tpSigTQCR", varid) )
       call check( nf90_get_var(ncid, varid, qc_tair,start=(/ 1, n /),count=(/ nsig, 1 /)) )
       call check( nf90_inq_varid(ncid, "tdSigTQCR", varid) )
       call check( nf90_get_var(ncid, varid, qc_tdew,start=(/ 1, n /),count=(/ nsig, 1 /)) )
    else
       qc_pres = 0
       qc_tair = 0
       qc_tdew = 0
    endif

    do k = 1, nsig

      prespa = pres(k) * 100.0_r8

      if ( tair(k) /= tair_miss .and. qc_tair(k) == 0 ) then

        oerr = rawin_temp_error(pres(k))
        if ( tair(k) >= 180.0_r8 .and. tair(k) <= 330.0_r8 .and. oerr /= missing_r8 ) then

          call create_obs_type(lat, lon, prespa, VERTISPRESSURE, tair(k), &
                               RADIOSONDE_TEMPERATURE, oerr, oday, osec, qc, obs)
          call append_obs_to_seq(obs_seq, obs)

        end if

      end if

      if ( include_specific_humidity .and. tair(k) /= tair_miss .and. tdew(k) /= tdew_miss &
           .and. qc_tair(k) == 0 .and. qc_tdew(k) == 0 ) then

        ! tdew is the dewpoint depression
        dptk = tair(k) - tdew(k)
        qobs = specific_humidity(sat_vapor_pressure(dptk),    prespa)
        qsat = specific_humidity(sat_vapor_pressure(tair(k)), prespa)
        if (LH_err ) then
          qerr = rh_error_from_dewpt_and_temp(tair(k), dptk)
        else
          qerr = rawin_rel_hum_error(pres(k), tair(k), qobs / qsat)
        end if
        oerr = max(qerr * qsat, 0.0001_r8)
        if ( qobs > 0.0_r8 .and. qobs <= 0.070_r8 .and. qerr /= missing_r8 ) then
          call create_obs_type(lat, lon, prespa, VERTISPRESSURE, qobs, &
                               RADIOSONDE_SPECIFIC_HUMIDITY, oerr, oday, osec, qc, obs)
          call append_obs_to_seq(obs_seq, obs)
        end if

      endif

      if ( include_relative_humidity .and. tair(k) /= tair_miss .and. tdew(k) /= tdew_miss &
           .and. qc_tair(k) == 0 .and. qc_tdew(k) == 0 ) then

      ! tdew is the dewpoint depression
        dptk = tair(k) - tdew(k)
        if (LH_err ) then
          oerr = rh_error_from_dewpt_and_temp(tair(k), dptk)
        else
          rh = temp_and_dewpoint_to_rh(tair(k), dptk)
          oerr = rawin_rel_hum_error(pres(k), tair(k), rh)
        end if
        call create_obs_type(lat, lon, prespa, VERTISPRESSURE, rh, &
                             RADIOSONDE_RELATIVE_HUMIDITY, oerr, oday, osec, qc, obs)
        call append_obs_to_seq(obs_seq, obs)

      endif

      if ( include_dewpoint .and. tair(k) /= tair_miss .and. tdew(k) /= tdew_miss &
            .and. qc_tair(k) == 0 .and. qc_tdew(k) == 0 ) then

      ! tdew is the dewpoint depression
        dptk = tair(k) - tdew(k)
        rh = temp_and_dewpoint_to_rh(tair(k), dptk)
        oerr = dewpt_error_from_rh_and_temp(tair(k), rh)
        call create_obs_type(lat, lon, prespa, VERTISPRESSURE, dptk, &
                             RADIOSONDE_DEWPOINT, oerr, oday, osec, qc, obs)
        call append_obs_to_seq(obs_seq, obs)

      endif

    end do
    deallocate(pres, tair, tdew, qc_pres, qc_tair, qc_tdew)

  end if

  !  If desired, read the significant-level wind data, write to obs. seq.
  call check( nf90_inq_varid(ncid, 'numSigW', varid) )
  call check( nf90_get_var(ncid, varid, nsig,  start = (/ n /)) )
  
  if ( sigwnd .and. nsig <= nmaxsw ) then

    allocate(pres(nsig))  ;  allocate(wdir(nsig))  ;  allocate(wspd(nsig))
    allocate(qc_pres(nsig))  ;  allocate(qc_wdir(nsig))  ;  allocate(qc_wspd(nsig))

    !  read significant level data
    call check( nf90_inq_varid(ncid, 'htSigW', varid) )
    call check( nf90_get_var(ncid,varid,pres,start=(/ 1, n /),count=(/ nsig, 1 /)) )
    call check( nf90_get_att(ncid, varid, '_FillValue', pres_miss) )
    call check( nf90_inq_varid(ncid, 'wdSigW', varid) )
    call check( nf90_get_var(ncid,varid,wdir,start=(/ 1, n /),count=(/ nsig, 1 /)) )
    call check( nf90_get_att(ncid, varid, '_FillValue', wdir_miss) )
    call check( nf90_inq_varid(ncid, 'wsSigW', varid) )
    call check( nf90_get_var(ncid,varid,wspd,start=(/ 1, n /),count=(/ nsig, 1 /)) )
    call check( nf90_get_att(ncid, varid, '_FillValue', wspd_miss) )
    if (input_has_qc) then
       call check( nf90_inq_varid(ncid, "wdSigTQCR", varid) )
       call check( nf90_get_var(ncid, varid, qc_wdir,start=(/ 1, n /),count=(/ nsig, 1 /)) )
       call check( nf90_inq_varid(ncid, "wsSigTQCR", varid) )
       call check( nf90_get_var(ncid, varid, qc_wspd,start=(/ 1, n /),count=(/ nsig, 1 /)) )
    else
       qc_wdir = 0
       qc_wspd = 0
    endif

    do k = 1, nsig

      !  add data to the observation sequence here.
      if ( wdir(k) /= wdir_miss .and. wspd(k) /= wspd_miss .and. qc_wdir(k) == 0  &
           .and. qc_wspd(k) == 0 ) then

        call wind_dirspd_to_uv(wdir(k), wspd(k), uwnd, vwnd)
        oerr = rawin_wind_error(pres_alt_to_pres(pres(k)) * 0.01_r8)
        if ( abs(uwnd) <= 150.0_r8 .and. & 
             abs(vwnd) <= 150.0_r8 .and. oerr /= missing_r8 ) then

          call create_obs_type(lat, lon, pres(k), VERTISHEIGHT, uwnd, &
                               RADIOSONDE_U_WIND_COMPONENT, oerr, oday, osec, qc, obs)
          call append_obs_to_seq(obs_seq, obs)
          call create_obs_type(lat, lon, pres(k), VERTISHEIGHT, vwnd, &
                               RADIOSONDE_V_WIND_COMPONENT, oerr, oday, osec, qc, obs)
          call append_obs_to_seq(obs_seq, obs)

        end if

      end if

    end do
    deallocate(pres, wdir, wspd, qc_wdir, qc_wspd)

  end if

enddo sondeloop

call check( nf90_close(ncid) )
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, rawin_out_file)

!end of main program

contains

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

end program convert_madis_rawin
