! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_madis_satwnd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_madis_satwnd - program that reads a MADIS netCDF satellite
!                          wind observation file and writes a DART
!                          obs_seq file using the DART library routines.
!                          This version works on the standard METAR files.
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
!     adapted May 2011, nancy collins, ncar/image
!     - split from the metar version
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, missing_r8
use     utilities_mod, only : initialize_utilities, finalize_utilities
use  netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              increment_time, get_time, operator(-), GREGORIAN
use      location_mod, only : VERTISPRESSURE
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use        meteor_mod, only : wind_dirspd_to_uv
use       obs_err_mod, only : sat_wind_error, sat_wv_wind_error
use      obs_kind_mod, only : SAT_U_WIND_COMPONENT, SAT_V_WIND_COMPONENT
use          sort_mod, only : index_sort
use obs_utilities_mod, only : getvar_real, get_or_fill_QC, add_obs_to_seq, &
                              create_3d_obs, getvar_int, getdimlen

use           netcdf

implicit none

character(len=16),  parameter :: satwnd_netcdf_file = 'satwnd_input.nc'
character(len=129), parameter :: satwnd_out_file    = 'obs_seq.satwnd'
logical,            parameter :: exclude_special     = .true.

logical, parameter :: use_input_qc              = .true. 

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries


logical :: iruse, visuse, wvuse, allbands

integer  :: ncid, nobs, n, i, oday, osec, nused
logical  :: file_exist, first_obs
real(r8) :: wdir_miss, wspd_miss
integer  :: band_miss

! FIXME:from ssec version
!logical :: iruse, visuse, wvuse, file_exist, qifile, eefile

real(r8) :: uwnd, vwnd, oerr, qc
! end FIXME


real(r8), allocatable :: lat(:), lon(:), pres(:),  &
                         tobs(:), wdir(:), wspd(:)
integer,  allocatable :: band(:), tused(:)
integer,  allocatable :: qc_wdir(:), qc_wspd(:)
integer,  allocatable :: used(:), sorted_used(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time


!------------
! start of executable code
!------------

call initialize_utilities('convert_madis_satwnd')

! band selections
print*,'Do you want to include IR, VISIBLE, and WV data? (T/F, 3 times)'
read*, iruse, visuse, wvuse

! set a flag to tell us if we are going to select by band or accept all
if (.not. iruse  .or.  .not. visuse  .or.  .not. wvuse) then
   allbands = .false.
else
   allbands = .true.
endif

! put the reference date into DART format
call set_calendar_type(GREGORIAN)
comp_day0 = set_date(1970, 1, 1, 0, 0, 0)

first_obs = .true.


ncid = nc_open_file_readonly(satwnd_netcdf_file, 'convert_madis_satwnd')

call getdimlen(ncid, "recNum", nobs)

allocate(lat(nobs), lon(nobs), pres(nobs), tobs(nobs))
allocate(band(nobs), wdir(nobs), wspd(nobs))
allocate(tused(nobs), used(nobs), sorted_used(nobs))
allocate(qc_wdir(nobs), qc_wspd(nobs))

! read in the data arrays
call getvar_real(ncid, "obLat",       lat            ) ! latitude
call getvar_real(ncid, "obLon",       lon            ) ! longitude
call getvar_real(ncid, "pressure",    pres           ) ! pressure
call getvar_real(ncid, "windDir",     wdir, wdir_miss) ! wind direction
call getvar_real(ncid, "windSpd",     wspd, wspd_miss) ! wind speed
call getvar_real(ncid, "validTime",   tobs           ) ! observation time
call getvar_int (ncid, "satelliteWindMethod", band, band_miss) ! band

! if user says to use them, read in QCs if present
if (use_input_qc) then
   call getvar_int(ncid, "windDirQCR",   qc_wdir ) ! wind direction qc
   call getvar_int(ncid, "windSpdQCR",   qc_wspd ) ! wind speed qc
else
   qc_wdir = 0 ;  qc_wspd = 0
endif

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=satwnd_out_file, exist=file_exist)

if ( file_exist ) then

  ! existing file found, append to it
  call read_obs_seq(satwnd_out_file, 0, 0, 2*nobs, obs_seq)

else

  ! create a new one
  call init_obs_sequence(obs_seq, num_copies, num_qc, 2*nobs)
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
    if ( lon(n) ==  lon(used(i)) .and. &
         lat(n) ==  lat(used(i)) .and. &
        pres(n) == pres(used(i)) .and. &
        band(n) == band(used(i)) .and. &
        tobs(n) == tobs(used(i)) ) cycle obsloop1
  end do

  ! if selecting only certain bands, cycle if not wanted
  if (.not. allbands) then
     if (.not. iruse  .and. band(n) == 1) cycle obsloop1
     if (.not. visuse .and. band(n) == 2) cycle obsloop1
     if (.not. wvuse  .and. &
         (band(n) == 3  .or.  band(n) == 5  .or. band(n) == 7)) cycle obsloop1
  endif

  ! the 'used' array are the index numbers of used obs
  ! the 'tused' array are the times of those obs so we can
  ! sort them later by time.
  nused = nused + 1
  used(nused) = n
  tused(nused) = tobs(n)

end do obsloop1

! sort by time
call index_sort(tused, sorted_used, nused)

obsloop2: do i = 1, nused

  ! get the next unique observation in sorted time order
  n = used(sorted_used(i))

  ! compute time of observation
  time_obs = increment_time(comp_day0, nint(tobs(n)))

  ! extract actual time of observation in file into oday, osec.
  call get_time(time_obs, osec, oday)

  ! add wind component data to obs_seq
  if ( wdir(n) /= wdir_miss .and. qc_wdir(n) == 0 .and. &
       wspd(n) /= wspd_miss .and. qc_wspd(n) == 0 ) then

    call wind_dirspd_to_uv(wdir(n), wspd(n), uwnd, vwnd)
    oerr = sat_wind_error(pres(n)/100.0_r8)  ! comes in as pascals already

   !  perform sanity checks on observation errors and values
   if ( oerr == missing_r8 .or. wdir(n) < 0.0_r8 .or. wdir(n) > 360.0_r8 .or. &
      wspd(n) < 0.0_r8 .or. wspd(n) > 120.0_r8 )  cycle obsloop2

      call create_3d_obs(lat(n), lon(n), pres(n), VERTISPRESSURE, uwnd, &
                         SAT_U_WIND_COMPONENT, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

      call create_3d_obs(lat(n), lon(n), pres(n), VERTISPRESSURE, vwnd, &
                         SAT_V_WIND_COMPONENT, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

  endif

end do obsloop2

! need to wait to close file because in the loop it queries the
! report types.
call nc_close_file(ncid, 'convert_madis_satwnd')

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, satwnd_out_file)

! end of main program
call finalize_utilities()

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
