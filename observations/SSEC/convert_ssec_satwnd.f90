!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_ssec_satwnd - program that reads ASCII satellite wind data 
!                         from CIMSS/SSEC and writes an observation
!                         sequence file
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program convert_ssec_satwnd

use        types_mod, only : r8, missing_r8
use    utilities_mod, only : get_unit
use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                             increment_time, get_time, set_date, operator(-)
use     location_mod, only : VERTISPRESSURE
use obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                             static_init_obs_sequence, init_obs, write_obs_seq, &
                             append_obs_to_seq, init_obs_sequence, get_num_obs, &
                             set_copy_meta_data, set_qc_meta_data
use       meteor_mod, only : wind_dirspd_to_uv
use      obs_err_mod, only : sat_wind_error, sat_wv_wind_error
use     obs_kind_mod, only : SAT_U_WIND_COMPONENT, SAT_V_WIND_COMPONENT
use           netcdf

implicit none

character(len=16),  parameter :: ssec_sat_file = 'satwnd_input.txt'
character(len=129), parameter :: sat_wind_file = 'obs_seq.satwnd'

integer, parameter :: nmaxwnd = 50000,  &  ! maximum number of vectors
                      num_copies = 1,   &  ! number of copies in sequence
                      num_qc     = 1       ! number of QC entries

character (len=129) :: meta_data, header
character (len=19)  :: datestr
character (len=8)   :: datein
character (len=6)   :: sat
character (len=4)   :: band, hourin

logical :: iruse, visuse, wvuse, swiruse, file_exist, qcinfile

integer :: in_unit, i, days, secs, nused, iyear, imonth, iday, ihour, & 
           imin, isec, dsec, dday, dsecobs
real(r8) :: obs_window, minqc, lat, lon, pres, wdir, wspd, uwnd, vwnd, oerr, &
            latu(nmaxwnd), lonu(nmaxwnd), prsu(nmaxwnd), qc, qcin1, qcin2

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs
type(time_type)         :: time_anal

print*,'Enter the analysis time (yyyy-mm-dd_hh:mm:ss), window (hours), and min. qc'
read*, datestr, obs_window, minqc
print*,'Do you want to include IR, VISIBLE, WV data, and SW IR? (T/F, 4 times)'
read*, iruse, visuse, wvuse, swiruse 
dsecobs = nint(obs_window * 3600.0_r8)

call set_calendar_type(GREGORIAN)
read(datestr(1:4),   fmt='(i4)') iyear
read(datestr(6:7),   fmt='(i2)') imonth
read(datestr(9:10),  fmt='(i2)') iday
read(datestr(12:13), fmt='(i2)') ihour
read(datestr(15:16), fmt='(i2)') imin
read(datestr(18:19), fmt='(i2)') isec
time_anal = set_date(iyear, imonth, iday, ihour, imin, isec)
call get_time(time_anal, secs, days)

in_unit  = get_unit()  ;  qcinfile = .false.
open(unit=in_unit,  file = ssec_sat_file,  status='old')
read(in_unit,'(a100)') header
if ( header(75:76) == 'qi' ) qcinfile = .true.

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
inquire(file=sat_wind_file, exist=file_exist)
if ( file_exist ) then

  call read_obs_seq(sat_wind_file, 0, 0, 2*nmaxwnd, obs_seq)

else

  call init_obs_sequence(obs_seq, num_copies, num_qc, 2*nmaxwnd)
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
obsloop: do

  if ( qcinfile ) then
    read(in_unit,*,END=200) band, sat, datein, hourin, lat, lon, &
                            pres, wspd, wdir, qcin1, qcin2
    if ( qcin2 < minqc ) cycle obsloop
  else
    read(in_unit,*,END=200) band, sat, datein, hourin, lat, lon, pres, wspd, wdir
  end if

  read(datein(1:4), fmt='(i4)') iyear
  read(datein(5:6), fmt='(i2)') imonth
  read(datein(7:8), fmt='(i2)') iday
  read(hourin(1:2), fmt='(i2)') ihour
  read(hourin(3:4), fmt='(i2)') imin
  call get_time((time_anal-set_date(iyear, imonth, iday, ihour, imin, 0)), dsec, dday)
  if ( (dsec + dday * 86400) > dsecobs ) cycle obsloop

  !  check the satellite channel
  if ( trim(adjustl(band))      == 'IR'   .and. (.not. iruse)   ) cycle obsloop
  if ( trim(adjustl(band))      == 'SWIR' .and. (.not. swiruse) ) cycle obsloop
  if ( trim(adjustl(band))      == 'VIS'  .and. (.not. visuse)  ) cycle obsloop
  if ( trim(adjustl(band(1:2))) == 'WV'   .and. (.not. wvuse)   ) cycle obsloop
  lon = -lon  ;  if ( lon < 0.0_r8 )  lon = lon + 360.0_r8

  !  check to make sure this observation has not been used before
  do i = 1, nused
    if ( lon == lonu(i) .and. lat == latu(i) .and. pres == prsu(i) ) cycle obsloop
  end do
  qc = 1.0_r8

  if ( trim(adjustl(band(1:2))) == 'WV' ) then
    oerr = sat_wv_wind_error(pres)
  else
    oerr = sat_wind_error(pres)
  end if

  !  perform sanity checks on observation errors and values
  if ( oerr == missing_r8 .or. wdir < 0.0_r8 .or. wdir > 360.0_r8 .or. &
       wspd < 0.0_r8 .or. wspd > 120.0_r8 )  cycle obsloop

  call wind_dirspd_to_uv(wdir, wspd, uwnd, vwnd)

  call create_obs_type(lat, lon, pres * 100.0_r8, VERTISPRESSURE, uwnd, &
                       SAT_U_WIND_COMPONENT, oerr, days, secs, qc, obs)
  call append_obs_to_seq(obs_seq, obs)
  call create_obs_type(lat, lon, pres * 100.0_r8, VERTISPRESSURE, vwnd, &
                       SAT_V_WIND_COMPONENT, oerr, days, secs, qc, obs)
  call append_obs_to_seq(obs_seq, obs)

  nused = nused + 1
  latu(nused) = lat
  lonu(nused) = lon
  prsu(nused) = pres

end do obsloop
200 continue

close( in_unit)

if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, sat_wind_file)

! end of main program

contains

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

integer, intent(in)           :: okind, vcord, day, sec
real(r8), intent(in)          :: lat, lon, pres, obsv, oerr, qc
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

end program
