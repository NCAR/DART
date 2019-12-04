! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$


!! THIS PROGRAM IS NOT FINISHED!  there are observations available in
!! netcdf format, but there were some missing fields that made completing
!! this program impossible.  We believe that now the missing data might
!! be added to the netcdf format files, so this could be finished.
!! It's in the release as-is so if someone wants to try to complete it
!! they have something to start from.

program convert_pb_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_pb_netcdf - program that reads a prep_bufr file that has
!                       been converted to netcdf format by MET version 5.
!                       (Model Evaluation Tools) Developmental Testbed Center
!
! adapted from MADIS converters by Soyoung Ha and Nancy Collins,  3/27/2018
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, missing_r8
use     utilities_mod, only : initialize_utilities, finalize_utilities
use  netcdf_utilities_mod, only : nc_check
use  time_manager_mod, only : time_type, set_calendar_type, set_date, set_time, &
                              increment_time, get_time, operator(-), GREGORIAN, &
                              set_time_missing, print_date
use      location_mod, only : VERTISSURFACE, VERTISPRESSURE
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use           meteor_mod, only : sat_vapor_pressure, specific_humidity, & 
                                 wind_dirspd_to_uv, pres_alt_to_pres, &
                                 temp_and_dewpoint_to_rh
use          obs_err_mod  ! all
use dewpoint_obs_err_mod, only : dewpt_error_from_rh_and_temp, &
                                 rh_error_from_dewpt_and_temp
use         obs_kind_mod  ! FIXME, all for now
use             sort_mod,  only : index_sort
use obs_def_altimeter_mod, only : compute_altimeter
use     obs_utilities_mod, only : getvar_real_2d, get_or_fill_QC, add_obs_to_seq, &
                                  create_3d_obs, getvar_int, getdimlen, set_missing_name, &
                                  getvar_char

use netcdf

implicit none

character(len=128), parameter :: pb_netcdf_file = 'prep_bufr_input.nc'
character(len=128), parameter :: pb_out_file    = 'obs_seq.prep_bufr'

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

integer, parameter :: MXSTRLEN = 16   ! must match what's in input files.


real(r8), parameter :: def_elev = 0.0_r8

integer  :: ncid, nobs, nhdr, mxstr, nvars, n, i, oday, osec, nused, sec1, sec2
integer  :: hdr, hdr_arr_len, obs_arr_len
logical  :: file_exist, first_obs
real(r8) :: hdr_miss, obs_miss

real(r8) :: sfcp_miss, tair_miss, tdew_miss, wdir_miss, wspd_miss, uwnd, &
            vwnd, altim, palt, oerr, qobs, qerr, qsat, rh, slp_miss, elev_miss

! in the input file you should find:
character(len=MXSTRLEN), allocatable :: obs_qty(:), hdr_typ(:), hdr_sid(:), hdr_vld(:)
real(r8), allocatable :: obs_arr(:,:), hdr_arr(:,:)

! we are going to compute these:
real(r8), allocatable :: lat(:), lon(:), pres(:), qc(:), stat_elev(:)
real(r8), allocatable :: obs_val(:), obs_err(:), obs_elev(:)
integer,  allocatable :: obs_typ(:)

! FIXME: leftovers - may not be needed
integer,  allocatable :: tobs(:), plid(:), tused(:), used(:), sorted_used(:)
real(r8), allocatable :: sfcp(:), tair(:), slp(:), tdew(:), wdir(:), wspd(:)
integer,  allocatable :: qc_sfcp(:), qc_slp(:), qc_tair(:), qc_tdew(:), qc_wdir(:), qc_wspd(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time, time_ref


! start of executable code
call initialize_utilities('convert_pb_netcdf')

! use the normal calendar
call set_calendar_type(GREGORIAN)

first_obs = .true.

call nc_check( nf90_open(pb_netcdf_file, nf90_nowrite, ncid), &
             'convert_pb_netcdf', 'opening file '//trim(pb_netcdf_file) )

call getdimlen(ncid, "nobs", nobs)
call getdimlen(ncid, "nhdr", nhdr)
call getdimlen(ncid, "mxstr", mxstr)
call getdimlen(ncid, "hdr_arr_len", hdr_arr_len)
call getdimlen(ncid, "obs_arr_len", obs_arr_len)

print *, 'getting dimensions: ', mxstr, hdr_arr_len, obs_arr_len, nobs, nhdr

! FIXME: is there a better way here?
if (mxstr /= MXSTRLEN) then
   print *, 'program assumes string lengths are ', MXSTRLEN, ' while file has ', mxstr
   stop
endif

! data from the input file
allocate(obs_qty(nobs), hdr_typ(nhdr), hdr_sid(nhdr), hdr_vld(nhdr))
allocate(obs_arr(obs_arr_len, nobs), hdr_arr(hdr_arr_len, nhdr))

! things we are going to extract and use
allocate(lat(nobs), lon(nobs), pres(nobs), qc(nobs), obs_typ(nobs))
allocate(stat_elev(nobs), obs_elev(nobs))
allocate(obs_err(nobs), obs_val(nobs))
allocate(tobs(nobs))  ;  allocate(tused(nobs))
allocate(used(nobs))  ;  allocate(sorted_used(nobs))

! read in the data arrays
call getvar_char   (ncid, "obs_qty",   obs_qty)             ! quality flag
call getvar_real_2d(ncid, "obs_arr",   obs_arr,  obs_miss)  ! obs values plus vert
call getvar_char   (ncid, "hdr_typ",   hdr_typ)             ! message type (eg ADPUPA)
call getvar_char   (ncid, "hdr_sid",   hdr_sid)             ! station id
call getvar_char   (ncid, "hdr_vld",   hdr_vld)             ! time: YYYYMMDD_HHMMSS
call getvar_real_2d(ncid, "hdr_arr",   hdr_arr,  hdr_miss)  ! location

call nc_check( nf90_close(ncid), &
             'convert_pb_netcdf', 'closing file '//trim(pb_netcdf_file) )

! we want a reference time that is earlier than any time
! that exists in this file.  for now take the first time
! and subtract 10 days.  these should be 6 hour files so
! that's overkill.
time_ref = string_to_time(hdr_vld(1))
time_ref = time_ref - set_time(0, 10)

inloop1: do n = 1, nobs
   ! FIXME: check obs_miss, hdr_miss and loop if missing

   ! extract header number and use it to index:
   !   lat, lon, time, elevation if surface obs
   !   message type and observation value
print *, n, obs_arr(:, n)
   hdr = obs_arr(1, n)
print *, 'f90 hdr = ', hdr+1, ' c hdr = ', hdr

   hdr = hdr + 1

   call get_time(string_to_time(hdr_vld(hdr)) - time_ref, tobs(n))

   qc(n) = string_to_qc(obs_qty(n))

   lat(n)       = hdr_arr(1, hdr)
   lon(n)       = hdr_arr(2, hdr)
   stat_elev(n) = hdr_arr(3, hdr)

   pres(n)     = obs_arr(3, n)
   obs_elev(n) = obs_arr(4, n)
   obs_val(n)  = obs_arr(5, n)

   ! we think:
   ! obs_typ is a combination of messsage type and grib code
   ! once we know the type we can call to get the obs error

print *, n, hdr_typ(n), obs_arr(2, n)
  
print *, lat(n), lon(n), stat_elev(n), qc(n), pres(n), obs_elev(n), obs_val(n)
call print_date(string_to_time(hdr_vld(hdr)))


enddo inloop1

stop

nused = 0
obsloop1: do n = 1, nobs

  ! check the lat/lon values to see if they are ok
  if ( lat(n) >  90.0_r8 .or. lat(n) <  -90.0_r8 ) cycle obsloop1
  if ( lon(n) > 180.0_r8 .or. lon(n) < -180.0_r8 ) cycle obsloop1

  ! change lon from -180 to 180 into 0-360
  if ( lon(n) < 0.0_r8 )  lon(n) = lon(n) + 360.0_r8

!  ! Check for duplicate observations
!  do i = 1, nused
!    if ( lon(n) == lon(used(i)) .and. &
!         lat(n) == lat(used(i)) .and. &
!        tobs(n) == tobs(used(i)) ) cycle obsloop1
!  end do

  nused = nused + 1
  used(nused) = n
  tused(nused) = tobs(n)

enddo obsloop1

! sort indices only by time
call index_sort(tused, sorted_used, nused)

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=pb_out_file, exist=file_exist)

if ( file_exist ) then

  ! existing file found, append to it
  call read_obs_seq(pb_out_file, 0, 0, nvars*nobs, obs_seq)

else

  ! create a new one
  call init_obs_sequence(obs_seq, num_copies, num_qc, nvars*nobs)
  do i = 1, num_copies
    call set_copy_meta_data(obs_seq, i, 'Observation')
  end do
  do i = 1, num_qc
    call set_qc_meta_data(obs_seq, i, 'Data QC')
  end do

endif

obsloop2: do i = 1, nused

  ! get the next unique observation in sorted time order
  n = used(sorted_used(i))

  ! convert offset in seconds to a real time type again
  time_obs = increment_time(time_ref, tobs(n))

!  if ( elev(n) /= elev_miss ) then
!    palt = pres_alt_to_pres(elev(n)) * 0.01_r8
!  else
!    palt = pres_alt_to_pres(def_elev) * 0.01_r8
!  endif

  ! extract actual time of observation in file into oday, osec.
  call get_time(time_obs, osec, oday)

!  ! add altimeter data to obs_seq
!  if ( sfcp(n) /= sfcp_miss .and. qc_sfcp(n) == 0 .and. elev(n) /= elev_miss ) then
!
!    altim = compute_altimeter(sfcp(n) * 0.01_r8, elev(n))
!    if ( plid(n) == 0 ) then
!      oerr = fixed_marine_pres_error(palt)
!    else
!      oerr = moving_marine_pres_error(palt)
!    endif

!    if ( altim >= 890.0_r8 .and. altim <= 1100.0_r8 .and. oerr /= missing_r8 ) then
!
!      call create_3d_obs(lat(n), lon(n), elev(n), VERTISSURFACE, altim, &
!                         MARINE_SFC_ALTIMETER, oerr, oday, osec, qc(n), obs)
!      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
!
!    endif
!
!  !  if surface pressure and elevation do not exist, use SLP.
!  else if ( slp(n) /= slp_miss .and. qc_slp(n) == 0 ) then
!
!    altim = compute_altimeter(slp(n) * 0.01_r8, 0.0_r8)
!    if ( plid(n) == 0 ) then
!      oerr = fixed_marine_pres_error(palt)
!    else
!      oerr = moving_marine_pres_error(palt)
!    endif
!

      call create_3d_obs(lat(n), lon(n), pres(n), VERTISPRESSURE, uwnd, &
                         RADIOSONDE_U_WIND_COMPONENT, oerr, oday, osec, qc(n), obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

      call create_3d_obs(lat(n), lon(n), pres(n), VERTISPRESSURE, vwnd, &
                         RADIOSONDE_V_WIND_COMPONENT, oerr, oday, osec, qc(n), obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)


100 continue

end do obsloop2

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, pb_out_file)

! end of main program
call finalize_utilities()

contains

function string_to_time(str_time)
 character(len=*), intent(in) :: str_time
 type(time_type) :: string_to_time

integer :: yr, mo, dy, hr, mn, sc
integer :: rc

read(str_time, "(I4,I2,I2,X,I2,I2,I2)", iostat=rc) yr, mo, dy, hr, mn, sc
if (rc /= 0) then
   print *, 'error converting string to time: '//trim(str_time)
   string_to_time = set_time_missing()
   return
endif

string_to_time = set_date(yr, mo, dy, hr, mn, sc)

end function string_to_time


function string_to_qc(str_qc)
 character(len=*), intent(in) :: str_qc
 real(r8) :: string_to_qc

integer :: qc_val, rc

read(str_qc, "(I2)", iostat=rc) qc_val
if (rc /= 0) then
   print *, 'error converting string to qc: '//trim(str_qc)
   string_to_qc = 99
   return
endif

string_to_qc = qc_val

end function string_to_qc

subroutine look_up_grib_code(grib_type, var_type)
! Find the variable type from 
! http://www.nco.ncep.noaa.gov/pmb/docs/on388/table2.html
 integer,          intent(in)  :: grib_type
 character(len=*), intent(out) :: var_type

select case (grib_type)
!    case(2) 
!      var_type = 'MSLP'
    case(7)
      var_type = 'GEOPOTENTIAL_HGT'
    case(11)
      var_type = 'TEMPERATURE'
    case(17)
      var_type = 'DEWPOINT'
!    case(32)
!      var_type = 'WIND_SPEED'
    case(33)
      var_type = 'U_WIND_COMPONENT'
    case(34)
      var_type = 'V_WIND_COMPONENT'
    case(51)
      var_type = 'SPECIFIC_HUMIDITY'
    case(52)
      var_type = 'RELATIVE_HUMIDITY'
!    case(53)
!      var_type = 'HUMIDITY MIXING RATIO'
    case(54)
      var_type = 'PRECIPITABLE_WATER'
    case default
      var_type = 'NOT_SUPPORTED'
end select

end subroutine look_up_grib_code

subroutine look_up_message_type(msg_type, obstype)
! Find observation type from 
! http://www.emc.ncep.noaa.gov/mmb/data_processing/prepbufr.doc/table_4.htm
 character(len=6), intent(in)  :: msg_type
 character(len=*), intent(out) :: obstype 

select case (msg_type)
  case ('ADPUPA')
    obstype = 'RADIOSONDE'   ! FIXME: DROPSONDE, too.
  case ('AIRCAR')
    obstype = 'ACARS'
  case ('AIRCFT')
    obstype = 'AIRCRAFT'
  case ('SATWND')
    obstype = 'SAT'
  case ('PROFLR')
    obstype = 'PROFILER'
  case ('ADPSFC')
    obstype = 'METAR'        ! FIXME: SYNOP, too.
  case ('MSONET')
    obstype = 'LAND_SFC'
  case ('SFCSHP')
    obstype = 'MARINE'
  case ('GPSIPW')
    obstype = 'GPS'
  case ('QKSWND')
    obstype = 'QKSWND'
  case default
    obstype = 'UNKNOWN'
end select

end subroutine look_up_message_type

subroutine look_up_obs_type(msg_type, grib_type, dart_obs_type)
 character(len=*), intent(in)  :: msg_type
 integer,          intent(in)  :: grib_type
 integer,          intent(out) :: dart_obs_type

select case (msg_type)
  case ('ADPUPA')
     ! t, q, p, u, v, sst
     select case(grib_type)
        case (11)
           dart_obs_type = RADIOSONDE_TEMPERATURE
        case (12)
           dart_obs_type = RADIOSONDE_TEMPERATURE
        case default
           dart_obs_type = -1
      end select
  case ('AIRCAR')
  case ('AIRCFT')
  case ('SATWND')
  case ('PROFLR')
  case ('VADWND')
  case ('SATEMP')
  case ('ADPSFC')
  case ('SFCSHP')
  case ('SFCBOG')
  case ('SPSSMI')
  case ('SYNDAT')
  case ('ERS1DA')
  case ('GOESND')
  case ('QKSWND')
  case ('MSONET')
  case ('GPSIPW')
  case ('RASSDA')
  case ('WDSATR')
  case ('ASCATW')
  case default
     print *, 'unknown message type: '//trim(msg_type)
end select

end subroutine look_up_obs_type

end program convert_pb_netcdf

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
