! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_ok_mesonet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_ok_mesonet - program to convert Oklahoma Mesonet MDF files
!                        into DART observation sequence files.
!
!   The observation files can  be obtained from the Oklahoma Mesonet
!   archive using urls of the format:
!   http://www.mesonet.org/index.php/dataMdfMts/dataController/getFile/YYYYMMDDHHMM/mdf/TEXT/ 
!   where YYYYMMDDHHMM is the date and time of the desired set of observations.
!   files are available every 5 minutes.
!
!   Static fields are drawn from the station description file provided by
!   the OK Mesonet. Update the local file from:
!   http://www.mesonet.org/index.php/api/siteinfo/from_all_active_with_geo_fields/format/csv/ 
!
!   NOTE: you may want to consider using METAR surface ob errors
!   with Oklahoma Mesonet surface obs. See flag below.
!
!   Written by G. Romine using the madis converters as a template, Aug. 2013
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! may want to add a namelist - so leaving module components active for now
use         utilities_mod, only : get_unit, find_namelist_in_file, check_namelist_read, &
                                  do_nml_file, do_nml_term, logfileunit, nmlfileunit
use             types_mod, only : r8, missing_r8
use      time_manager_mod, only : time_type, set_calendar_type, set_date, &
                                  increment_time, get_time, operator(-), GREGORIAN
use          location_mod, only : location_type, set_location, get_location, &
                                  get_dist, VERTISSURFACE
use      obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                                  static_init_obs_sequence, init_obs, write_obs_seq, &
                                  init_obs_sequence, get_num_obs, &
                                  set_copy_meta_data, set_qc_meta_data
use           obs_def_mod, only : obs_def_type
use            meteor_mod, only : sat_vapor_pressure, specific_humidity, &
                                  wind_dirspd_to_uv, rh_and_temp_to_dewpoint
use           obs_err_mod, only : land_temp_error, land_wind_error, &
                                  land_pres_error, land_rel_hum_error, &
                                  metar_temp_error, metar_wind_error, &
                                  metar_pres_error, metar_rel_hum_error
use  dewpoint_obs_err_mod, only : dewpt_error_from_rh_and_temp, &
                                  rh_error_from_dewpt_and_temp
use          obs_kind_mod, only : LAND_SFC_U_WIND_COMPONENT, LAND_SFC_V_WIND_COMPONENT, &
                                  LAND_SFC_TEMPERATURE, LAND_SFC_SPECIFIC_HUMIDITY, &
                                  LAND_SFC_DEWPOINT, LAND_SFC_RELATIVE_HUMIDITY, &
                                  LAND_SFC_ALTIMETER
use     obs_utilities_mod, only : add_obs_to_seq, create_3d_obs
use obs_def_altimeter_mod, only : compute_altimeter

implicit none
! Mesonet observations are station lists at a common time with observations in
! the following order:
! stid (A4), stnm (I3), time (I4), relh (I4), tair (F6.1), wspd (F6.1), wvec
! (F6.1), wdir (I3), wdsd (F6.1), wssd (F6.1), wmax (F6.1), rain (F7.2), pres
! (F7.2), srad (I4), ta9m (F6.1), ws2m (F6.1), ts10 (F6.1), tb10 (F6.1), ts05
! (F6.1), tb05 (F6.1), ts30 (F6.1), tr05 (F6.1), tr25 (F6.1), tr60 (F6.1)
! first 3 lines are header, date info drawn from middle line

character(len=15),  parameter :: omeso_in_file  = 'okmeso_mdf.in'
character(len=129), parameter :: omeso_out_file = 'obs_seq.okmeso'

integer, parameter ::   nmax_out   = 4000, &   ! maximum number of reports for file
                        num_copies = 1,    &   ! number of copies in sequence
                        num_qc     = 1         ! number of QC entries

! Allow these as namelist entries? Default is convert it all
logical :: include_specific_humidity = .true.
logical :: include_relative_humidity = .true.
logical :: include_dewpoint          = .true.
logical :: LH_err                    = .true.
logical :: use_metar_ob_errors       = .true.

real(r8), parameter :: fmiss         = -999.0_r8 ! -998 is also used for some 
                                                 ! unused fields
logical  :: fexist, first_obs
integer :: iunit, obs_num, i
real(r8) :: qc
! fields in the mdf
real(r8)          :: tair, wspd, wvec, wdsd, wssd, wmax, rain, pres, ta9m, &
                     ws2m, ts10, tb10, ts05, tb05, ts30, tr05, tr25, tr60
integer           :: stnm, time, relh, wdir, srad
character(len=4)  :: stid
integer           :: y4, m2, d2, h2, n2, s2, tobs, osec, oday
! fields from the geo file
real(r8)          :: nlat, nlon, elev 
! local
real(r8)          :: uwnd, vwnd, alti, oerr, qerr, qv, dptk, qsat, tmpk
character(len=129) :: header, meta_data
! obs sequence vars
type(location_type)     :: obs_loc
type(obs_sequence_type) :: obs_seq
type(obs_def_type)      :: obs_def
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

obs_num   = 1
first_obs = .true.
qc        = 1.0_r8

inquire(file = omeso_in_file, exist = fexist)
if ( .NOT. fexist ) stop

iunit = get_unit()
open(unit=iunit, file = omeso_in_file, status = 'old')

!  either read existing obs_seq or create a new one
call set_calendar_type(GREGORIAN)
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
inquire(file=omeso_out_file, exist=fexist)

if ( fexist ) then
  call read_obs_seq(omeso_out_file, 0, 0, 4*nmax_out, obs_seq)
else
  call init_obs_sequence(obs_seq, num_copies, num_qc, 4*nmax_out)
  do i = 1, num_copies
    meta_data = 'OKMESO observation'
    call set_copy_meta_data(obs_seq, i, meta_data)
  end do
  do i = 1, num_qc
    meta_data = 'Data QC'
    call set_qc_meta_data(obs_seq, i, meta_data)
  end do
end if

 read(iunit,*) ! header
! get the date info from the file
 read(iunit,'(A50)',END=200) header
 read(header,12) y4, m2, d2
12 format(5X,I4,2(X,I2))
 write(*,*) 'Date ', y4, m2, d2
! the time is gathered below - minutes since 00 UTC
 read(iunit,*) ! header  need date and time from here

! loop through all of the obs until reaching the end of the file
obsloop: do
!  do while (1 .eq. 1)
    read(iunit,22,END=200) stid, stnm, time, relh, tair, wspd, wvec, wdir, &
                           wdsd, wssd, wmax, rain, pres

22 format(1X,A4,3X,I3,2x,I4,3X,I4,3(1X,F6.1),2X,I4,3(1X,F6.1), &
             1X,F7.2,2X,F7.2)
! Given the station id, get the lat, lon, and elevation
    call get_geo(stid,nlat,nlon,elev)
! Convert the wind speed and direction to u and v wind components, tair to K
    call wind_dirspd_to_uv(1.0_r8*wdir, wspd, uwnd, vwnd)
    tmpk = tair + 273.15_r8
! convert time to UTC
    h2 = time/60
    n2 = time - h2*60
    tobs = h2*100 + n2
    s2 = 0
    time_obs = set_date(y4, m2, d2, h2, n2, s2)
    call get_time(time_obs, osec, oday)
! NEEDED put in checks for missing obs before each type is added
! compute surface altimeter
    alti = compute_altimeter(pres, elev)   ! pres hPa and elev in m
    
!Debug    write(*,*) stid, nlat, nlon, elev, tobs, relh, tair, uwnd, vwnd, pres, alti

!    nlon = -nlon  
    if ( nlon < 0.0_r8 )  nlon = nlon + 360.0_r8

! Gen obs sequence
! altimeter
  if ( pres /= fmiss  ) then
    if (use_metar_ob_errors) then
      oerr = metar_pres_error(pres)
    else
      oerr = land_pres_error(pres)
    end if

    if ( alti >=  890.0_r8 .and. &
         alti <= 1100.0_r8 .and. oerr /= missing_r8 ) then
   
      call create_3d_obs(nlat, nlon, elev, VERTISSURFACE, alti, &
                         LAND_SFC_ALTIMETER, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
      obs_num = obs_num + 1
    end if
  end if

! winds
  if ( wspd /= fmiss .and. wdir /= fmiss ) then
    if (use_metar_ob_errors) then
      oerr = metar_wind_error(pres)
    else
      oerr = land_wind_error(pres)
    end if

    if ( abs(uwnd) < 150.0_r8 .and. abs(vwnd) < 150.0_r8 .and. oerr /= missing_r8 ) then

      call create_3d_obs(nlat, nlon, elev, VERTISSURFACE, uwnd, &
                         LAND_SFC_U_WIND_COMPONENT, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

      call create_3d_obs(nlat, nlon, elev, VERTISSURFACE, vwnd, &
                         LAND_SFC_V_WIND_COMPONENT, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
      obs_num = obs_num + 1

    endif
  endif

! temperature
  if ( tair /= fmiss  ) then
    if (use_metar_ob_errors) then
      oerr = metar_temp_error(pres)
    else
      oerr = land_temp_error(pres)
    endif
    if ( tmpk >= 200.0_r8 .and. tmpk <= 335.0_r8 .and. oerr /= missing_r8) then

      call create_3d_obs(nlat, nlon, elev, VERTISSURFACE, tmpk, &
                         LAND_SFC_TEMPERATURE, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
      obs_num = obs_num + 1

    endif
  endif

! moisture
  if ( tair /= fmiss .and. relh /= fmiss .and. pres /= fmiss ) then

    qsat = specific_humidity(sat_vapor_pressure(tmpk), pres * 100.0_r8)

    if ( include_specific_humidity ) then

      qv   = qsat * relh / 100.0_r8
      if (use_metar_ob_errors) then
        qerr = metar_rel_hum_error(pres, tmpk, qv / qsat)
      else
        qerr = land_rel_hum_error(pres, tmpk, qv / qsat)
      endif
      if ( LH_err ) then
        dptk = rh_and_temp_to_dewpoint(relh/100.0_r8, tmpk)  ! Kelvin temp, rh 0.0-1.0
        qerr = rh_error_from_dewpt_and_temp(tmpk, dptk)
      endif

      oerr = max(qerr * qsat, 0.0001_r8)

      if ( qv >= 0.0_r8 .and. oerr /= missing_r8) then

        call create_3d_obs(nlat, nlon, elev, VERTISSURFACE, qv, &
                           LAND_SFC_SPECIFIC_HUMIDITY, oerr, oday, osec, qc, obs)
        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
        obs_num = obs_num + 1

      end if

    end if

    if ( include_relative_humidity ) then

      if (use_metar_ob_errors) then
        oerr = metar_rel_hum_error(pres, tmpk, relh/100.0_r8)
      else
        oerr = land_rel_hum_error(pres, tmpk, relh/100.0_r8)
      endif
      if ( LH_err ) then
        dptk = rh_and_temp_to_dewpoint(relh/100.0_r8, tmpk)  ! Kelvin temp, rh 0.0-1.0
        oerr = rh_error_from_dewpt_and_temp(tmpk, dptk)
      endif

      if ( relh >=  0.0_r8 .and. relh <= 101.0_r8 .and. oerr /= missing_r8) then

        call create_3d_obs(nlat, nlon, elev, VERTISSURFACE, relh * 0.01_r8, &
                           LAND_SFC_RELATIVE_HUMIDITY, oerr, oday, osec, qc, obs)
        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
        obs_num = obs_num + 1

      end if

    end if

    if ( include_dewpoint ) then

      dptk = rh_and_temp_to_dewpoint(relh/100.0_r8, tmpk)  ! Kelvin temp, rh 0.0-1.0
      oerr = dewpt_error_from_rh_and_temp(tmpk, relh/100.0_r8)  ! This is the LH_err module

      if ( relh >=  0.0_r8 .and. relh <= 101.0_r8 .and. oerr /= missing_r8) then

        call create_3d_obs(nlat, nlon, elev, VERTISSURFACE, dptk, &
                           LAND_SFC_DEWPOINT, oerr, oday, osec, qc, obs)
        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
        obs_num = obs_num + 1

      end if

    end if

  end if
 
  end do obsloop
200 continue
close(iunit) 
! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, omeso_out_file)
   write (*,fmt='(A,I5,A)') 'Created ',obs_num, ' obs'
end program convert_ok_mesonet
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_geo(search_stid,nlat,nlon,elevf)

use         utilities_mod, only : get_unit, find_namelist_in_file, check_namelist_read, &
                                  do_nml_file, do_nml_term, logfileunit, nmlfileunit
use             types_mod, only : r8, missing_r8

implicit none

character(len=4), intent(in)  :: search_stid
real(r8), intent(out) :: nlat, nlon, elevf

! Given a station id, find the lat, lon and elevation from the geo file
! provided by the Oklahoma Mesonent website
!
! the geoinfo.csv file has data in the following order:
!
! stnm - station number
! stid - CHAR len=4
! name - CHAR of variable size
! city - CHAR of variable size
! rang - distance of station from city, float
! cdir - CHAR of variable size
! cnty - CHAR of variable size
! nlat - latitude, float
! nlon - longitude, float
! elev - elevation, float
! + other stuff we don't need
!

! temp vars
integer :: stnm, iunit, elev
character(len=4) :: stid
character(len=15) :: name, city, cdir, cnty
real :: rang
character(len=330) :: line
integer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, nend
logical :: stn_found


 stn_found = .false.
 iunit = get_unit()
 open(unit=iunit, file='geoinfo.csv', status='old')

 read(iunit,*) ! header

  do while (1 .eq. 1)
! read through each line and look for a match 
 read(iunit,'(A330)',END=200) line

! Ugh... a csv file with mixed var types is a 
! pain. Find the locations of delimiters and 
! assign vars as whatever is found between them. 
! Getting out to the 10th item gives us all we
! need. 
! Determine locations of the first 10 delimiters
n1  = index(line, ',')
nend = len_trim(line)
n2  = n1 + index(line(n1+1:nend), ',')
n3  = n2 + index(line(n2+1:nend), ',')
n4  = n3 + index(line(n3+1:nend), ',')
n5  = n4 + index(line(n4+1:nend), ',')
n6  = n5 + index(line(n5+1:nend), ',')
n7  = n6 + index(line(n6+1:nend), ',')
n8  = n7 + index(line(n7+1:nend), ',')
n9  = n8 + index(line(n8+1:nend), ',')
n10 = n9 + index(line(n9+1:nend), ',')

!
 read (line(1:n1-1),'(I3)') stnm
 read (line(n1+1:n2-1),'(A4)') stid
 read (line(n2+1:n3-1),'(A15)') name
 read (line(n3+1:n4-1),'(A15)') city
 read (line(n4+1:n5-1),'(f6.3)') rang
 read (line(n5+1:n6-1),'(A3)') cdir
 read (line(n6+1:n7-1),'(A15)') cnty
 read (line(n7+1:n8-1),'(f10.7)') nlat
 read (line(n8+1:n9-1),'(f10.7)') nlon
 read (line(n9+1:n10-1),'(I4)') elev
 elevf = 1.0_r8 * elev
  if (search_stid .eq. stid) then
    stn_found = .true.
    goto 200
  end if

  end do
200  continue
  if (.not. stn_found) then
    write(*,*) ' station ',stid,' not found. Update the geo file.'
  end if
  close(iunit)
end subroutine get_geo

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
