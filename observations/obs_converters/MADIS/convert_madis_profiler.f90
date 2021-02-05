! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_madis_profiler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_madis_profiler - program that reads a MADIS netCDF profiler
!                          wind observation file and writes a DART
!                          obs_seq file using the DART library routines.
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM
!     modified Dec. 2008 Soyoung Ha and David Dowell, NCAR/MMM
!     modified to include QC_flag check (Soyoung Ha, NCAR/MMM, 08-04-2009)
!
!     modified to use a common set of utilities, better netcdf error checks,
!     able to insert obs with any time correctly (not only monotonically
!     increasing times)    nancy collins,  ncar/image   11 march 2010
!     
!     keep original obs times, make source for all converters as similar
!     as possible.   nancy collins,  ncar/image   26 march 2010
!
!     adapted June 2011, nancy collins, ncar/image, Glen Romine ncar/mmm
!     - split from the satwind version
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, missing_r8
use     utilities_mod, only : initialize_utilities, finalize_utilities
use  netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              increment_time, get_time, operator(-), GREGORIAN
use      location_mod, only : VERTISHEIGHT
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use       obs_err_mod, only : prof_wind_error
use      obs_kind_mod, only : PROFILER_U_WIND_COMPONENT, PROFILER_V_WIND_COMPONENT
use          sort_mod, only : index_sort
use obs_utilities_mod, only : getvar_real, get_or_fill_QC, add_obs_to_seq, &
                              create_3d_obs, getvar_int, getdimlen, getvar_real_2d, &
                              getvar_int_2d, query_varname

use           netcdf

implicit none

character(len=20),  parameter :: profiler_netcdf_file = 'profiler_input.nc'
character(len=129), parameter :: profiler_out_file    = 'obs_seq.profiler'

logical, parameter :: use_input_qc              = .true. 

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries


integer  :: ncid, nsta, nlev, n, i, oday, osec, nused, k, index
logical  :: file_exist, first_obs
real(r8) :: uwnd_miss, vwnd_miss

real(r8) :: oerr, qc

real(r8), allocatable :: lat(:), lon(:), latu(:), lonu(:), &
                         levs(:,:), elev(:), tobs(:), tobu(:), &
                         uwnd(:,:), vwnd(:,:), levu(:)
integer,  allocatable :: qc_uwnd(:,:), qc_vwnd(:,:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

character(len=NF90_MAX_NAME) :: namelist(5)

!------------
! start of executable code
!------------

call initialize_utilities('convert_madis_profiler')


! put the reference date into DART format
call set_calendar_type(GREGORIAN)
comp_day0 = set_date(1970, 1, 1, 0, 0, 0)

first_obs = .true.


ncid = nc_open_file_readonly(profiler_netcdf_file, 'convert_madis_profiler')

call getdimlen(ncid, "recNum", nsta)
call getdimlen(ncid, "level" , nlev)

allocate( lat(nsta))         ;  allocate( lon(nsta))
allocate(latu(nsta*nlev))    ;  allocate(lonu(nsta*nlev))
allocate(tobs(nsta))         ;  allocate(tobu(nsta*nlev))
allocate(elev(nsta))         ;  allocate(levu(nsta*nlev))
! 2d variables - levels at each station
allocate(uwnd(nlev,nsta))    ;  allocate(vwnd(nlev,nsta))
allocate(qc_uwnd(nlev,nsta)) ;  allocate(qc_vwnd(nlev,nsta))
allocate(levs(nlev,nsta))  

! read in the data arrays

! we have profiler data files which have different names for the 
! lat/lon/elev/obs arrays in the netcdf file.  there doesn't seem
! to be a global attr to say which one is in use, so for now try
! both options.  

namelist(1) = 'staLat'
namelist(2) = 'latitude'
call query_varname(ncid, 2, namelist, index, force=.true.)
call    getvar_real(ncid, namelist(index),  lat            ) ! station latitude

namelist(1) = 'staLon'
namelist(2) = 'longitude'
call query_varname(ncid, 2, namelist, index, force=.true.)
call    getvar_real(ncid, namelist(index),  lon            ) ! station longitude

namelist(1) = 'staElev'
namelist(2) = 'elevation'
call query_varname(ncid, 2, namelist, index, force=.true.)
call    getvar_real(ncid, namelist(index),  elev           ) ! station elevation

namelist(1) = 'timeObs'
namelist(2) = 'observationTime'
call query_varname(ncid, 2, namelist, index, force=.true.)
call    getvar_real(ncid, namelist(index),  tobs           ) ! observation time

call getvar_real_2d(ncid, "levels",      levs           ) ! height above station in meters
call getvar_real_2d(ncid, "uComponent",  uwnd, uwnd_miss) ! e-w component
call getvar_real_2d(ncid, "vComponent",  vwnd, vwnd_miss) ! n-s component

! if user says to use them, read in QCs if present
if (use_input_qc) then
   call getvar_int_2d(ncid, "uComponentQCR",   qc_uwnd ) ! wind direction qc
   call getvar_int_2d(ncid, "vComponentQCR",   qc_vwnd ) ! wind speed qc
else
   qc_uwnd = 0 ;  qc_vwnd = 0
endif

! levels is height above station.  we want actual elevation which means
! adding on the station elevation.  add it in here so when we use levs
! below it is the actual height above MSL.
do n = 1, nsta
   levs(:,n) = levs(:,n) + elev(n)
enddo

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=profiler_out_file, exist=file_exist)

if ( file_exist ) then

  ! existing file found, append to it
  call read_obs_seq(profiler_out_file, 0, 0, 2*nsta*nlev, obs_seq)

else

  ! create a new one
  call init_obs_sequence(obs_seq, num_copies, num_qc, 2*nsta*nlev)
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
staloop: do n = 1, nsta

  ! compute time of observation
  time_obs = increment_time(comp_day0, nint(tobs(n)))

  ! check the lat/lon values to see if they are ok
  if ( lat(n) >  90.0_r8 .or. lat(n) <  -90.0_r8 ) cycle staloop
  if ( lon(n) > 180.0_r8 .or. lon(n) < -180.0_r8 ) cycle staloop

  if ( lon(n) < 0.0_r8 )  lon(n) = lon(n) + 360.0_r8
  ! extract actual time of observation in file into oday, osec.
  call get_time(time_obs, osec, oday)

levloop: do k = 1, nlev
    ! Check for duplicate observations
    do i = 1, nused
      if ( lon(n)    == lonu(i) .and. &
           lat(n)    == latu(i) .and. &
           levs(k,n) == levu(i) .and. &
           tobs(n)   == tobu(i) ) cycle levloop
    end do


  ! add wind component data to obs_seq
  if ( uwnd(k,n) /= uwnd_miss .and. qc_uwnd(k,n) == 0 .and. &
       vwnd(k,n) /= vwnd_miss .and. qc_vwnd(k,n) == 0 ) then

     ! FIXME: we do not have pressure in these files, only height
     ! need pres(n) = convert_std_atm(height(n))
     !oerr = prof_wind_error(pres(n))  ! pressure based table
     oerr = 2.000_r8  ! this works mostly for surface to 100mb

   !  perform sanity checks on observation errors and values
   if ( oerr == missing_r8 .or. &
      abs(uwnd(k,n)) > 150.0_r8 .or. &
      abs(vwnd(k,n)) > 150.0_r8 )  cycle levloop

      call create_3d_obs(lat(n), lon(n), levs(k,n), VERTISHEIGHT, uwnd(k,n), &
                         PROFILER_U_WIND_COMPONENT, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

      call create_3d_obs(lat(n), lon(n), levs(k,n), VERTISHEIGHT, vwnd(k,n), &
                         PROFILER_V_WIND_COMPONENT, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

  endif

  nused = nused + 1
  latu(nused) = lat(n)
  lonu(nused) = lon(n)
  levu(nused) = levs(k,n)
  tobu(nused) = tobs(n)

end do levloop
end do staloop

! need to wait to close file because in the loop it queries the
! report types.
call nc_close_file(ncid, 'convert_madis_profiler')

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, profiler_out_file)

! end of main program
call finalize_utilities()


end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
