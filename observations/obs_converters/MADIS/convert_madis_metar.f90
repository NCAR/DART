! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_madis_metar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_madis_metar - program that reads a MADIS netCDF land 
!                         surface observation file and writes a DART
!                         obs_seq file using the DART library routines.
!                         This version works on the standard METAR files.
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM
!     modified Dec. 2008 Soyoung Ha and David Dowell, NCAR/MMM
!     - added dewpoint as an output variable
!     - added relative humidity as an output variable
!
!     modified to include QC_flag check (Soyoung Ha, NCAR/MMM, 08-04-2009)
!     split from the mesonet version (Glen Romine, NCAR/MMM, Feb 2010)
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
use  netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file, nc_check
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                                  increment_time, get_time, operator(-), GREGORIAN
use      location_mod, only : VERTISSURFACE
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use        meteor_mod, only : sat_vapor_pressure, specific_humidity, & 
                              wind_dirspd_to_uv, invert_altimeter, pres_alt_to_pres, &
                              temp_and_dewpoint_to_rh
use       obs_err_mod, only : land_temp_error, land_wind_error, &
                              land_pres_error, land_rel_hum_error
use dewpoint_obs_err_mod, only : dewpt_error_from_rh_and_temp, &
                                 rh_error_from_dewpt_and_temp
use      obs_kind_mod, only : METAR_U_10_METER_WIND, METAR_V_10_METER_WIND, &
                              METAR_TEMPERATURE_2_METER, METAR_SPECIFIC_HUMIDITY_2_METER, & 
                              METAR_DEWPOINT_2_METER, METAR_RELATIVE_HUMIDITY_2_METER, &
                              METAR_ALTIMETER                
use             sort_mod,  only : index_sort
use  obs_utilities_mod, only : getvar_real, get_or_fill_QC, add_obs_to_seq, &
                               create_3d_obs, getvar_int, getdimlen, set_missing_name

use           netcdf

implicit none

! COMPILE TIME OPTION:  by default this converter only processes hourly
! observations.  if exclude_special is set to .false., then special obs
! will be converted as well as the normal obs.  but this is not on by default.
character(len=16),  parameter :: surface_netcdf_file = 'metar_input.nc'
character(len=129), parameter :: surface_out_file    = 'obs_seq.metar'
logical,            parameter :: exclude_special     = .true.

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

character (len=5)   :: rtype
integer  :: ncid, nobs, nvars, n, i, oday, osec, nused
logical  :: file_exist, first_obs
real(r8) :: alti_miss, tair_miss, tdew_miss, wdir_miss, wspd_miss, uwnd, &
            vwnd, palt, qobs, qsat, rh, oerr, pres, qerr, qc

integer,  allocatable :: tobs(:), tused(:), used(:), sorted_used(:)
real(r8), allocatable :: lat(:), lon(:), elev(:), alti(:), tair(:), & 
                         tdew(:), wdir(:), wspd(:)
integer,  allocatable :: qc_alti(:), qc_tair(:), qc_tdew(:), qc_wdir(:), qc_wspd(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time


call initialize_utilities('convert_madis_metar')

! put the reference date into DART format
call set_calendar_type(GREGORIAN)
comp_day0 = set_date(1970, 1, 1, 0, 0, 0)

first_obs = .true.


ncid = nc_open_file_readonly(surface_netcdf_file, 'convert_madis_metar')

call getdimlen(ncid, "recNum", nobs)
call set_missing_name("missing_value")

allocate( lat(nobs))  ;  allocate( lon(nobs))
allocate(elev(nobs))  ;  allocate(alti(nobs))
allocate(tair(nobs))  ;  allocate(tdew(nobs))
allocate(wdir(nobs))  ;  allocate(wspd(nobs))
allocate(tobs(nobs))  ;  allocate(tused(nobs))
allocate(used(nobs))  ;  allocate(sorted_used(nobs))

nvars = 4
if (include_specific_humidity) nvars = nvars + 1
if (include_relative_humidity) nvars = nvars + 1
if (include_dewpoint)          nvars = nvars + 1

allocate(qc_alti(nobs))
allocate(qc_tair(nobs)) ;  allocate(qc_tdew(nobs))
allocate(qc_wdir(nobs)) ;  allocate(qc_wspd(nobs))

! read in the data arrays
call getvar_real(ncid, "latitude",    lat            ) ! latitude
call getvar_real(ncid, "longitude",   lon            ) ! longitude
call getvar_real(ncid, "elevation",   elev           ) ! elevation
call getvar_real(ncid, "altimeter",   alti, alti_miss) ! altimeter setting
call getvar_real(ncid, "temperature", tair, tair_miss) ! air temperature
call getvar_real(ncid, "dewpoint",    tdew, tdew_miss) ! dew-point temperature
call getvar_real(ncid, "windDir",     wdir, wdir_miss) ! wind direction
call getvar_real(ncid, "windSpeed",   wspd, wspd_miss) ! wind speed
call getvar_int (ncid, "timeObs",     tobs           ) ! observation time

! if user says to use them, read in QCs if present
if (use_input_qc) then
   call get_or_fill_QC(ncid, "altimeterQCR",   qc_alti)
   call get_or_fill_QC(ncid, "temperatureQCR", qc_tair)
   call get_or_fill_QC(ncid, "dewpointQCR",    qc_tdew)
   call get_or_fill_QC(ncid, "windDirQCR",     qc_wdir)
   call get_or_fill_QC(ncid, "windSpeedQCR",   qc_wspd)
else
   qc_alti = 0
   qc_tair = 0 ;  qc_tdew = 0
   qc_wdir = 0 ;  qc_wspd = 0
endif

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=surface_out_file, exist=file_exist)

if ( file_exist ) then

  ! existing file found, append to it
  call read_obs_seq(surface_out_file, 0, 0, nvars*nobs, obs_seq)

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

  ! check to make sure this observation has not been used
  call getvar_char(ncid, "reportType", n, rtype)
  if ( rtype /= 'METAR' .and. exclude_special )  cycle obsloop1

  ! check the lat/lon values to see if they are ok
  if ( lat(n) >  90.0_r8 .or. lat(n) <  -90.0_r8 ) cycle obsloop1
  if ( lon(n) > 180.0_r8 .or. lon(n) < -180.0_r8 ) cycle obsloop1

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

! sort obs by time
call index_sort(tused, sorted_used, nused)

obsloop2: do i = 1, nused

  ! get the next unique observation in sorted time order
  n = used(sorted_used(i))

  ! compute time of observation
  time_obs = increment_time(comp_day0, tobs(n))

  palt = pres_alt_to_pres(elev(n)) * 0.01_r8

  ! extract actual time of observation in file into oday, osec.
  call get_time(time_obs, osec, oday)

  ! add altimeter data to obs_seq
  if ( alti(n) /= alti_miss .and. qc_alti(n) == 0 ) then

    pres = invert_altimeter(alti(n) * 0.01_r8, elev(n))
    oerr = land_pres_error(palt)
    if ( alti(n) >=  89000.0_r8 .and. &
         alti(n) <= 110000.0_r8 .and. oerr /= missing_r8 ) then

      call create_3d_obs(lat(n), lon(n), elev(n), VERTISSURFACE, alti(n) * 0.01_r8, & 
                         METAR_ALTIMETER, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

    endif

  endif

  ! add wind component data to obs_seq
  if ( wdir(n) /= wdir_miss .and. qc_wdir(n) == 0 .and. &
       wspd(n) /= wspd_miss .and. qc_wspd(n) == 0  ) then

    call wind_dirspd_to_uv(wdir(n), wspd(n), uwnd, vwnd)
    oerr = land_wind_error(palt)
    if ( abs(uwnd) < 150.0_r8 .and. abs(vwnd) < 150.0_r8 .and. oerr /= missing_r8 ) then

      call create_3d_obs(lat(n), lon(n), elev(n), VERTISSURFACE, uwnd, &
                         METAR_U_10_METER_WIND, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

      call create_3d_obs(lat(n), lon(n), elev(n), VERTISSURFACE, vwnd, &
                         METAR_V_10_METER_WIND, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

    endif

  endif

  ! add air temperature data to obs_seq
  if ( tair(n) /= tair_miss .and. qc_tair(n) == 0 ) then 

    oerr = land_temp_error(palt)
    if ( tair(n) >= 200.0_r8 .and. tair(n) <= 335.0_r8 .and. oerr /= missing_r8 ) then

      call create_3d_obs(lat(n), lon(n), elev(n), VERTISSURFACE, tair(n), &
                         METAR_TEMPERATURE_2_METER, oerr, oday, osec, qc, obs)
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
       if (tdew(n) > tair(n) + 1.0_r8) cycle obsloop2
       tdew(n) = tair(n)
    endif

    ! add specific humidity to obs_seq
    if ( include_specific_humidity .and. &
         alti(n) /= alti_miss .and. qc_alti(n) == 0 ) then

      qobs = specific_humidity(sat_vapor_pressure(tdew(n)), pres * 100.0_r8)
      qsat = specific_humidity(sat_vapor_pressure(tair(n)), pres * 100.0_r8)
      if ( LH_err ) then
        qerr = rh_error_from_dewpt_and_temp(tair(n), tdew(n))
      else
        qerr = land_rel_hum_error(pres, tair(n), qobs / qsat)
      endif
      oerr = max(qerr * qsat, 0.0001_r8)
  
      if ( qobs >  0.0_r8  .and. &
           qobs <= 0.07_r8 .and. qerr /= missing_r8 ) then
  
        call create_3d_obs(lat(n), lon(n), elev(n), VERTISSURFACE, qobs, &
                           METAR_SPECIFIC_HUMIDITY_2_METER, oerr, oday, osec, qc, obs)
        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
      endif
  
    endif
  
    ! add relative humidity data to obs_seq
    if ( include_relative_humidity ) then
  
      rh = temp_and_dewpoint_to_rh(tair(n), tdew(n))
      if ( LH_err ) then
        oerr = rh_error_from_dewpt_and_temp(tair(n), tdew(n))
      else
        oerr = land_rel_hum_error(pres, tair(n), rh)    
      endif
  
      if ( rh >  0.0_r8 .and. &
           rh <= 1.5_r8 .and. oerr /= missing_r8 ) then
  
        call create_3d_obs(lat(n), lon(n), elev(n), VERTISSURFACE, rh, &
                           METAR_RELATIVE_HUMIDITY_2_METER, oerr, oday, osec, qc, obs)
        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
      endif
  
    endif
  
    ! add dew-point temperature data to obs_seq
    if ( include_dewpoint ) then
  
      rh = temp_and_dewpoint_to_rh(tair(n), tdew(n))
      oerr = dewpt_error_from_rh_and_temp(tair(n), rh)
  
      if ( rh >  0.0_r8 .and. &
           rh <= 1.5_r8 .and. oerr /= missing_r8 ) then
  
        call create_3d_obs(lat(n), lon(n), elev(n), VERTISSURFACE, tdew(n), &
                           METAR_DEWPOINT_2_METER, oerr, oday, osec, qc, obs)
        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
      endif
  
  !    print*, 'temp (C), rh (%), oerr:  ', tair(n)-273.15_r8, rh*100.0_r8, oerr
  
    endif

  endif  ! quality control/missing check on tair, tdew

end do obsloop2

! need to wait to close file because in the loop it queries the
! report types.
call nc_close_file(ncid, 'convert_madis_metar')

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, surface_out_file)

! end of main program
call finalize_utilities()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   getvar_char - subroutine that inquires, gets the char variable
!     SPECIALIZED for this use - assumes start = (/ 1, n /)
!           so takes a scalar start, returns a character string
!
!      ncid - open netcdf file handle
!      varname - string name of netcdf variable
!      start - starting index in the 2d array.  integer
!      darray - output array.  character
!
!     created 11 Mar 2010,  nancy collins,  ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getvar_char(ncid, varname, start, darray)
 integer,            intent(in)   :: ncid
 character(len = *), intent(in)   :: varname
 integer,            intent(in)   :: start
 character(len = *), intent(out)  :: darray

integer :: varid

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'getvar_char', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, darray, start=(/1, start/) ), &
               'getvar_char', 'getting var '// trim(varname))

end subroutine getvar_char


end program convert_madis_metar

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
