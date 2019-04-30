! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_ssec_satwnd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_ssec_satwnd - program that reads ASCII satellite wind data 
!                         from CIMSS/SSEC and writes an observation
!                         sequence file
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use          types_mod, only : r8, missing_r8
use      utilities_mod, only : get_unit
use   time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                               increment_time, get_time, set_date, operator(-)
use       location_mod, only : VERTISPRESSURE
use   obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                               static_init_obs_sequence, init_obs, write_obs_seq, &
                               append_obs_to_seq, init_obs_sequence, get_num_obs, &
                               set_copy_meta_data, set_qc_meta_data
use         meteor_mod, only : wind_dirspd_to_uv
use        obs_err_mod, only : sat_wind_error, sat_wv_wind_error
use       obs_kind_mod, only : SAT_U_WIND_COMPONENT, SAT_V_WIND_COMPONENT
use  obs_utilities_mod, only : create_3d_obs, add_obs_to_seq

implicit none

character(len=16),  parameter :: ssec_sat_file = 'satwnd_input.txt'
character(len=129), parameter :: sat_wind_file = 'obs_seq.satwnd'

integer, parameter :: nmaxwnd = 50000,  &  ! maximum number of vectors
                      num_copies = 1,   &  ! number of copies in sequence
                      num_qc     = 1       ! number of QC entries

character(len=129) :: header
character(len=8)   :: datein
character(len=6)   :: sat
character(len=4)   :: band, hourin

logical :: iruse, visuse, wvuse, swiruse, file_exist, qifile, eefile, &
           userfqc, useqiqc, useeeqc, first_obs

integer :: in_unit, i, oday, osec, nused, iyear, imonth, iday, ihour, & 
           imin, qctype
real(r8) :: lat, lon, pres, wdir, wspd, uwnd, vwnd, oerr, latu(nmaxwnd), &
            lonu(nmaxwnd), prsu(nmaxwnd), qc, qcthresh, rfqc, qiqc, eeqc

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

print*,'Do you want to include IR, VISIBLE, WV data, and SW IR? (T/F, 4 times)'
read*, iruse, visuse, wvuse, swiruse 

call set_calendar_type(GREGORIAN)
first_obs = .true.

qifile  = .false.  ;  eefile  = .false.
userfqc = .false.  ;  useqiqc = .false.  ;  useeeqc = .false.  

in_unit  = get_unit()
open(unit=in_unit,  file = ssec_sat_file,  status='old')
read(in_unit,'(a100)') header

if ( header(81:82) == 'ee' ) then

  print*,'Enter the QC type (0=none, 1=rf, 2=qi, 3=ee) and threshhold'
  read*, qctype, qcthresh

  select case ( qctype )
    case (1)  !  Regression Factor QC
      userfqc = .true.
    case (2)  !  QI format QC
      useqiqc = .true.
    case (3)
      useeeqc = .true.
  end select
  eefile = .true.

else if ( header(75:76) == 'qi' ) then
  
  print*,'Enter the QC type (0=none, 1=rf, 2=qi) and threshhold' 
  read*, qctype, qcthresh

  select case ( qctype )
    case (1)  !  Regression Factor QC
      userfqc = .true.
    case (2)  !  QI format QC
      useqiqc = .true.
  end select
  qifile = .true.

end if

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
inquire(file=sat_wind_file, exist=file_exist)
if ( file_exist ) then

  call read_obs_seq(sat_wind_file, 0, 0, 2*nmaxwnd, obs_seq)

else

  call init_obs_sequence(obs_seq, num_copies, num_qc, 2*nmaxwnd)
  do i = 1, num_copies
    call set_copy_meta_data(obs_seq, i, 'SSEC observation')
  end do
  do i = 1, num_qc
    call set_qc_meta_data(obs_seq, i, 'Data QC')
  end do

end if

nused = 0
obsloop: do

  if ( qifile ) then  !  QI format file

    read(in_unit,*,END=200) band, sat, datein, hourin, lat, lon, &
                            pres, wspd, wdir, rfqc, qiqc
    if ( userfqc .and. rfqc > qcthresh )  cycle obsloop
    if ( useqiqc .and. qiqc < qcthresh )  cycle obsloop

  else if ( eefile ) then  !  EE format file
  
    read(in_unit,*,END=200) band, sat, datein, hourin, lat, lon, &
                            pres, wspd, wdir, rfqc, qiqc, eeqc
    if ( userfqc .and. (rfqc > qcthresh) )  cycle obsloop
    if ( useqiqc .and. (qiqc < qcthresh) )  cycle obsloop
    if ( useeeqc .and. (eeqc > qcthresh) )  cycle obsloop

  else  !  file without QC information

    read(in_unit,*,END=200) band, sat, datein, hourin, lat, lon, pres, wspd, wdir

  end if

  read(datein(1:4), fmt='(i4)') iyear
  read(datein(5:6), fmt='(i2)') imonth
  read(datein(7:8), fmt='(i2)') iday
  read(hourin(1:2), fmt='(i2)') ihour
  read(hourin(3:4), fmt='(i2)') imin
  time_obs = set_date(iyear, imonth, iday, ihour, imin, 0)
  call get_time(time_obs, osec, oday)

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

!  if ( trim(adjustl(band(1:2))) == 'WV' ) then
!    oerr = sat_wv_wind_error(pres)
!  else
    oerr = sat_wind_error(pres)
!  end if

  !  perform sanity checks on observation errors and values
  if ( oerr == missing_r8 .or. wdir < 0.0_r8 .or. wdir > 360.0_r8 .or. &
       wspd < 0.0_r8 .or. wspd > 120.0_r8 )  cycle obsloop

  call wind_dirspd_to_uv(wdir, wspd, uwnd, vwnd)

  call create_3d_obs(lat, lon, pres * 100.0_r8, VERTISPRESSURE, uwnd, &
                     SAT_U_WIND_COMPONENT, oerr, oday, osec, qc, obs)
  call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

  call create_3d_obs(lat, lon, pres * 100.0_r8, VERTISPRESSURE, vwnd, &
                     SAT_V_WIND_COMPONENT, oerr, oday, osec, qc, obs)
  call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

  nused = nused + 1
  latu(nused) = lat
  lonu(nused) = lon
  prsu(nused) = pres

end do obsloop
200 continue

close( in_unit)

if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, sat_wind_file)

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
