! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_madis_marine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_madis_marine - program that reads a MADIS netCDF marine
!                         surface observation file and writes a DART
!                         obs_seq file using the DART library routines.
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM
!     modified Dec. 2008 Soyoung Ha and David Dowell, NCAR/MMM
!     - added dewpoint as an output variable
!     - added relative humidity as an output variable
!
!
!     modified to include QC_flag check (Soyoung Ha, NCAR/MMM, 08-04-2009)
!
!     modified to use a common set of utilities, better netcdf error checks,
!     able to insert obs with any time correctly (not only monotonically
!     increasing times)    nancy collins,  ncar/image   11 march 2010
!     
!     keep original obs times, make source for all converters as similar
!     as possbile.   nancy collins,  ncar/image   26 march 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, missing_r8
use     utilities_mod, only : initialize_utilities, finalize_utilities
use  netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              increment_time, get_time, operator(-), GREGORIAN
use      location_mod, only : VERTISSURFACE
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use           meteor_mod, only : sat_vapor_pressure, specific_humidity, & 
                                 wind_dirspd_to_uv, pres_alt_to_pres, &
                                 temp_and_dewpoint_to_rh
use          obs_err_mod, only : fixed_marine_temp_error, fixed_marine_rel_hum_error, &
                                 fixed_marine_wind_error, fixed_marine_pres_error, &
                                 moving_marine_temp_error, moving_marine_rel_hum_error, &
                                 moving_marine_wind_error, moving_marine_pres_error
use dewpoint_obs_err_mod, only : dewpt_error_from_rh_and_temp, &
                                 rh_error_from_dewpt_and_temp
use         obs_kind_mod, only : MARINE_SFC_U_WIND_COMPONENT, MARINE_SFC_V_WIND_COMPONENT, &
                                 MARINE_SFC_TEMPERATURE, MARINE_SFC_SPECIFIC_HUMIDITY,     &
                                 MARINE_SFC_ALTIMETER, MARINE_SFC_DEWPOINT,   &
                                 MARINE_SFC_RELATIVE_HUMIDITY
use             sort_mod,  only : index_sort
use obs_def_altimeter_mod, only : compute_altimeter
use     obs_utilities_mod, only : getvar_real, get_or_fill_QC, add_obs_to_seq, &
                                  create_3d_obs, getvar_int, getdimlen, set_missing_name

implicit none

character(len=15),  parameter :: marine_netcdf_file = 'marine_input.nc'
character(len=129), parameter :: marine_out_file    = 'obs_seq.marine'

! the following logical parameters control which water-vapor variables appear in the output file,
! whether to use the NCEP error or Lin and Hubbard (2004) moisture error model, and if the
! input file has data quality control fields, whether to use or ignore them.
logical, parameter :: LH_err                    = .false.
logical, parameter :: include_specific_humidity = .true.
logical, parameter :: include_relative_humidity = .false.
logical, parameter :: include_dewpoint          = .false.
logical, parameter :: use_input_qc              = .true. 

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

real(r8), parameter :: def_elev = 0.0_r8

integer  :: ncid, nobs, nvars, n, i, oday, osec, nused
logical  :: file_exist, first_obs
real(r8) :: sfcp_miss, tair_miss, tdew_miss, wdir_miss, wspd_miss, uwnd, &
            vwnd, altim, palt, oerr, qobs, qerr, qsat, rh, slp_miss, elev_miss, qc

integer,  allocatable :: tobs(:), plid(:), tused(:), used(:), sorted_used(:)
real(r8), allocatable :: lat(:), lon(:), elev(:), sfcp(:), tair(:), slp(:), & 
                         tdew(:), wdir(:), wspd(:)
integer,  allocatable :: qc_sfcp(:), qc_slp(:), qc_tair(:), qc_tdew(:), qc_wdir(:), qc_wspd(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time


call initialize_utilities('convert_madis_marine')

! put the reference date into DART format
call set_calendar_type(GREGORIAN)
comp_day0 = set_date(1970, 1, 1, 0, 0, 0)

first_obs = .true.


ncid = nc_open_file_readonly(marine_netcdf_file, 'convert_madis_marine')

call getdimlen(ncid, "recNum", nobs)
call set_missing_name("missing_value")

allocate( lat(nobs))  ;  allocate( lon(nobs))
allocate(elev(nobs))  ;  allocate(plid(nobs))
allocate(sfcp(nobs))  ;  allocate(slp(nobs))
allocate(tair(nobs))  ;  allocate(tdew(nobs))
allocate(wdir(nobs))  ;  allocate(wspd(nobs))
allocate(tobs(nobs))  ;  allocate(tused(nobs))
allocate(used(nobs))  ;  allocate(sorted_used(nobs))

nvars = 4
if (include_specific_humidity) nvars = nvars + 1
if (include_relative_humidity) nvars = nvars + 1
if (include_dewpoint)          nvars = nvars + 1

allocate(qc_sfcp(nobs)) ;  allocate(qc_slp(nobs))
allocate(qc_tair(nobs)) ;  allocate(qc_tdew(nobs))
allocate(qc_wdir(nobs)) ;  allocate(qc_wspd(nobs))

! read in the data arrays
call getvar_real(ncid, "latitude",         lat            ) ! latitude
call getvar_real(ncid, "longitude",        lon            ) ! longitude
call getvar_real(ncid, "elevation",        elev, elev_miss) ! elevation
call getvar_int (ncid, "dataPlatformType", plid           ) ! platform type
call getvar_real(ncid, "stationPress",     sfcp, sfcp_miss) ! altimeter setting
call getvar_real(ncid, "seaLevelPress",    slp,  slp_miss ) ! sea-level pressure
call getvar_real(ncid, "temperature",      tair, tair_miss) ! air temperature
call getvar_real(ncid, "dewpoint",         tdew, tdew_miss) ! dew-point temperature
call getvar_real(ncid, "windDir",          wdir, wdir_miss) ! wind direction
call getvar_real(ncid, "windSpeed",        wspd, wspd_miss) ! wind speed
call getvar_int (ncid, "timeObs",          tobs           ) ! observation time

! if user says to use them, read in QCs if present
if (use_input_qc) then
   call get_or_fill_QC(ncid, "stationPressQCR",  qc_sfcp)
   call get_or_fill_QC(ncid, "seaLevelPressQCR", qc_slp )
   call get_or_fill_QC(ncid, "temperatureQCR",   qc_tair)
   call get_or_fill_QC(ncid, "dewpointQCR",      qc_tdew)
   call get_or_fill_QC(ncid, "windDirQCR",       qc_wdir)
   call get_or_fill_QC(ncid, "windSpeedQCR",     qc_wspd)
else
   qc_sfcp = 0 ;  qc_slp  = 0
   qc_tair = 0 ;  qc_tdew = 0
   qc_wdir = 0 ;  qc_wspd = 0
endif

call nc_close_file(ncid, 'convert_madis_marine')


!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=marine_out_file, exist=file_exist)

if ( file_exist ) then

  ! existing file found, append to it
  call read_obs_seq(marine_out_file, 0, 0, nvars*nobs, obs_seq)

else

  ! create a new one
  call init_obs_sequence(obs_seq, num_copies, num_qc, nvars*nobs)
  do i = 1, num_copies
    call set_copy_meta_data(obs_seq, i, 'MADIS observation')
  end do
  do i = 1, num_qc
    call set_qc_meta_data(obs_seq, i, 'Data QC')
  end do

endif

! Set the DART data quality control.  Be consistent with NCEP codes;
! 0 is 'must use', 1 is good, no reason not to use it.
qc = 1.0_r8

nused = 0
obsloop1: do n = 1, nobs

  ! check the lat/lon values to see if they are ok
  if ( lat(n) >  90.0_r8 .or. lat(n) <  -90.0_r8 ) cycle obsloop1
  if ( lon(n) > 180.0_r8 .or. lon(n) < -180.0_r8 ) cycle obsloop1

  ! change lon from -180 to 180 into 0-360
  if ( lon(n) < 0.0_r8 )  lon(n) = lon(n) + 360.0_r8

  ! Check for duplicate observations
  do i = 1, nused
    if ( lon(n) == lon(used(i)) .and. &
         lat(n) == lat(used(i)) .and. &
        tobs(n) == tobs(used(i)) ) cycle obsloop1
  end do

  nused = nused + 1
  used(nused) = n
  tused(nused) = tobs(n)

enddo obsloop1

! sort indices only by time
call index_sort(tused, sorted_used, nused)

obsloop2: do i = 1, nused

  ! get the next unique observation in sorted time order
  n = used(sorted_used(i))

  ! compute time of observation
  time_obs = increment_time(comp_day0, tobs(n))

  if ( elev(n) /= elev_miss ) then
    palt = pres_alt_to_pres(elev(n)) * 0.01_r8
  else
    palt = pres_alt_to_pres(def_elev) * 0.01_r8
  endif

  ! extract actual time of observation in file into oday, osec.
  call get_time(time_obs, osec, oday)

  ! add altimeter data to obs_seq
  if ( sfcp(n) /= sfcp_miss .and. qc_sfcp(n) == 0 .and. elev(n) /= elev_miss ) then

    altim = compute_altimeter(sfcp(n) * 0.01_r8, elev(n))
    if ( plid(n) == 0 ) then
      oerr = fixed_marine_pres_error(palt)
    else
      oerr = moving_marine_pres_error(palt)
    endif

    if ( altim >= 890.0_r8 .and. altim <= 1100.0_r8 .and. oerr /= missing_r8 ) then

      call create_3d_obs(lat(n), lon(n), elev(n), VERTISSURFACE, altim, &
                         MARINE_SFC_ALTIMETER, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

    endif

  !  if surface pressure and elevation do not exist, use SLP.
  else if ( slp(n) /= slp_miss .and. qc_slp(n) == 0 ) then

    altim = compute_altimeter(slp(n) * 0.01_r8, 0.0_r8)
    if ( plid(n) == 0 ) then
      oerr = fixed_marine_pres_error(palt)
    else
      oerr = moving_marine_pres_error(palt)
    endif

    if ( altim >= 890.0_r8 .and. altim <= 1100.0_r8 .and. oerr /= missing_r8 ) then

      call create_3d_obs(lat(n), lon(n), def_elev, VERTISSURFACE, altim, &
                         MARINE_SFC_ALTIMETER, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

    endif

  endif
  if ( elev(n) == elev_miss )  elev(n) = def_elev

  ! add wind component data to obs_seq
  if ( wdir(n) /= wdir_miss .and. qc_wdir(n) == 0 .and. &
       wspd(n) /= wspd_miss .and. qc_wspd(n) == 0  ) then

    call wind_dirspd_to_uv(wdir(n), wspd(n), uwnd, vwnd)
    if ( plid(n) == 0 ) then
      oerr = fixed_marine_wind_error(palt)
    else
      oerr = moving_marine_wind_error(palt)
    endif
    if ( abs(uwnd) < 150.0_r8 .and. abs(vwnd) < 150.0_r8 .and. oerr /= missing_r8 ) then

      call create_3d_obs(lat(n), lon(n), elev(n), VERTISSURFACE, uwnd, &
                         MARINE_SFC_U_WIND_COMPONENT, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

      call create_3d_obs(lat(n), lon(n), elev(n), VERTISSURFACE, vwnd, &
                         MARINE_SFC_V_WIND_COMPONENT, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

    endif

  endif

  ! add air temperature data to obs_seq
  if ( tair(n) /= tair_miss .and. qc_tair(n) == 0 ) then 

    if ( plid(n) == 0 ) then
      oerr = fixed_marine_temp_error(palt)
    else
      oerr = moving_marine_temp_error(palt)
    endif
    if ( tair(n) >= 200.0_r8 .and. tair(n) <= 335.0_r8 .and. oerr /= missing_r8 ) then

      call create_3d_obs(lat(n), lon(n), elev(n), VERTISSURFACE, tair(n), &
                         MARINE_SFC_TEMPERATURE, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

    endif

  endif

  ! if the air, dewpoint are ok, then see which of
  ! the three types of moisture obs to generate.
  if ( tair(n) /= tair_miss .and. qc_tair(n) == 0 .and. &
       tdew(n) /= tdew_miss .and. qc_tdew(n) == 0  ) then

       ! before we start computing things based on the dewpoint,
       ! make sure it isn't larger than the air temp.  if it is
       ! more than a degree larger, skip it completely.  if it is
       ! less, set them equal and continue.
       if (tdew(n) > tair(n)) then
          if (tdew(n) > tair(n) + 1.0_r8) goto 100
          tdew(n) = tair(n)
       endif

    ! add specific humidity to obs_seq
    if ( include_specific_humidity .and. &
         sfcp(n) /= sfcp_miss .and. qc_sfcp(n) == 0 ) then

      qobs = specific_humidity(sat_vapor_pressure(tdew(n)), sfcp(n))
      qsat = specific_humidity(sat_vapor_pressure(tair(n)), sfcp(n))
      if ( LH_err ) then
        qerr = rh_error_from_dewpt_and_temp(tair(n), tdew(n))
      else
        if ( plid(n) == 0 ) then
          !GSR - there was a bug here with the rh error assigned to oerr 
          !instead of qerr, leaving qerr undefined for calc below
          qerr = fixed_marine_rel_hum_error(palt, tair(n), qobs / qsat)
        else
          qerr = moving_marine_rel_hum_error(palt, tair(n), qobs / qsat)
        endif
      endif
      oerr = max(qerr * qsat, 0.0001_r8)
  
      if ( qobs >  0.0_r8  .and. &
           qobs <= 0.07_r8 .and. qerr /= missing_r8 ) then
  
        call create_3d_obs(lat(n), lon(n), def_elev, VERTISSURFACE, qobs, &
                           MARINE_SFC_SPECIFIC_HUMIDITY, oerr, oday, osec, qc, obs)
        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
      endif
  
    endif
  
    ! add relative humidity data to obs_seq
    if ( include_relative_humidity ) then
  
      rh = temp_and_dewpoint_to_rh(tair(n), tdew(n))
      if ( LH_err ) then
        oerr = rh_error_from_dewpt_and_temp(tair(n), tdew(n))
      else
        if ( plid(n) == 0 ) then
          oerr = fixed_marine_rel_hum_error(palt, tair(n), rh)
        else
          oerr = moving_marine_rel_hum_error(palt, tair(n), rh)
        endif
      endif
  
      if ( rh >  0.0_r8 .and. &
           rh <= 1.5_r8 .and. oerr /= missing_r8 ) then

        call create_3d_obs(lat(n), lon(n), def_elev, VERTISSURFACE, rh, &
                           MARINE_SFC_RELATIVE_HUMIDITY, oerr, oday, osec, qc, obs)
        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
      endif
    endif
  
    ! add dew-point temperature data to obs_seq
    if ( include_dewpoint ) then
  
      rh = temp_and_dewpoint_to_rh(tair(n), tdew(n))
      oerr = dewpt_error_from_rh_and_temp(tair(n), rh)
  
      if ( rh >  0.0_r8 .and. &
           rh <= 1.5_r8 .and. oerr /= missing_r8 ) then

        call create_3d_obs(lat(n), lon(n), def_elev, VERTISSURFACE, tdew(n), &
                           MARINE_SFC_DEWPOINT, oerr, oday, osec, qc, obs)
        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
      endif
    endif

  endif  ! quality control/missing check on tair, tdew

100 continue

end do obsloop2

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, marine_out_file)

! end of main program
call finalize_utilities()

end program convert_madis_marine

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
