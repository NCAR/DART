! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_f16_edr_dsk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! convert_f16_edr_dsk - program that reads a SSUSI netCDF profiler
!                       wind observation file and writes a DART
!                       obs_seq file using the DART library routines.
!
! http://ssusi.jhuapl.edu/data_products
!
! The ON2_UNCERTAINTY variable in the netcdf files have IEEE NaN values,
! but none of the required metadata to interpret them correctly.
! These 2 lines will add the required attributes so that NaNs are replaced with
! a fill value that can be queried and checked for.
! Since the ON2_UNCERTAINTY is a standard deviation, it is enough to make it negative
!
! ncatted -a _FillValue,ON2_UNCERTAINTY,o,f,NaN        input_file.nc
! ncatted -a _FillValue,ON2_UNCERTAINTY,m,f,-1.0       input_file.nc
!
! These commands exist in the shell_scripts/netcdf_manip.csh script.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


use         types_mod, only : r8, MISSING_R8, digits12

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              find_namelist_in_file, check_namelist_read, do_output, &
                              error_handler, E_ERR, E_MSG

use  netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file

use  time_manager_mod, only : time_type, set_calendar_type, set_date, GREGORIAN, &
                              get_time, set_time, print_time, print_date

use      location_mod, only : VERTISUNDEF

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data

use      obs_kind_mod, only : SSUSI_O_N2_RATIO

use obs_utilities_mod, only : getdimlen, getvar_real, getvar_real_2d, &
                              getvar_int, getvar_int_2d, add_obs_to_seq, create_3d_obs

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   '$URL$'
character(len=*), parameter :: revision = '$Revision$'
character(len=*), parameter :: revdate  = '$Date$'


integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

integer  :: i, io, iunit, ncid, nused
integer  :: n_pix_along_day, n_pix_across_day, i_across, i_along
integer  :: seconds, days
logical  :: file_exist, first_obs

real(r8) :: oerr, qc, lon, lat, obs_value, org_missing

real(r8),       allocatable :: latitude(:,:), longitude(:,:)
integer,        allocatable :: year(:), doy(:)
real(digits12), allocatable :: time(:)
real(r8),       allocatable :: ON2(:,:), ON2_uncertainty(:,:)
integer,        allocatable :: ON2_quality(:,:)

! PIERCEPOINT_DAY_LATITUDE  &
!   "Dayside latitude of the pierce point; rebinned to the new grid." ;
! PIERCEPOINT_DAY_LONGITUDE &
!   "Dayside longitude of the pierce point; rebinned to the new grid." ;
! ON2                       &
!   "Ratio of the O to N2 vertical column densities on the disk" ;
! ON2_UNCERTAINTY           &
!   "Standard deviation of the ratio of O to N2 vertical column densities on the disk" ;
! DATA_QUALITY_DISK         &
!   "Disk data quality flags - see EDR docs for details" ;
 
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

character(len=512) :: string1, string2

!-------------------------------------------------------------------------------
! namelist with default values

character(len=256) :: input_netcdf_file  = 'input_netcdf.nc'
character(len=256) :: output_obs_file    = 'obs_seq.out'
integer            :: debug = 0

namelist /convert_f16_edr_dsk_nml/ input_netcdf_file, output_obs_file, debug
 
!-------------------------------------------------------------------------------
! start of executable code
!-------------------------------------------------------------------------------

call initialize_utilities('convert_f16_edr_dsk')

! Read the namelist entry for model_mod from input.nml
call find_namelist_in_file('input.nml', 'convert_f16_edr_dsk_nml', iunit)
read(iunit, nml = convert_f16_edr_dsk_nml, iostat = io)
call check_namelist_read(iunit, io, 'convert_f16_edr_dsk_nml')

ncid = nc_open_file_readonly(input_netcdf_file, 'convert_f16_edr_dsk')

call getdimlen(ncid,"N_PIX_ALONG_DAY",  n_pix_along_day)
call getdimlen(ncid,"N_PIX_ACROSS_DAY", n_pix_across_day)

allocate(       latitude(n_pix_across_day,n_pix_along_day))
allocate(      longitude(n_pix_across_day,n_pix_along_day))
allocate(            ON2(n_pix_across_day,n_pix_along_day))
allocate(ON2_uncertainty(n_pix_across_day,n_pix_along_day))
allocate(    ON2_quality(n_pix_across_day,n_pix_along_day))

allocate(year(n_pix_along_day))
allocate( doy(n_pix_along_day))
allocate(time(n_pix_along_day))

call getvar_int( ncid, 'YEAR', year ) 
call getvar_int( ncid, 'DOY' , doy  )
call getvar_real(ncid, 'TIME', time )

call getvar_real_2d(ncid, 'PIERCEPOINT_DAY_LATITUDE' , latitude )
call getvar_real_2d(ncid, 'PIERCEPOINT_DAY_LONGITUDE', longitude )
call getvar_real_2d(ncid, 'ON2'                      , ON2 )
call getvar_real_2d(ncid, 'ON2_UNCERTAINTY'          , ON2_uncertainty, org_missing)
call getvar_int_2d(ncid, 'DATA_QUALITY_DISK'        , ON2_quality )

! This test does not prevent the NaN failure unfortunately.
write(string1,*,iostat=io)'count of missing ON2_uncertainty values',count(ON2_uncertainty < 0.0_r8)
if (io /= 0) then
   write(string1,*)'['//trim(input_netcdf_file)//'] appears to have NaN values.'
   write(string2,*)'Remove them with shell_scripts/netcdf_manip.csh and retry.'
   call error_handler(E_ERR, 'convert_f16_edr_dsk', string1, &
              source, revision, revdate, text2=string2)
else
   call error_handler(E_MSG, 'convert_f16_edr_dsk:', string1)
endif

where (ON2_uncertainty == org_missing) ON2_uncertainty = MISSING_R8

call nc_close_file(ncid, 'convert_f16_edr_dsk')

! Ensure latitudes are within [-90,90]
where( latitude < -90.0_r8 ) latitude = -90.0_r8
where( latitude >  90.0_r8 ) latitude =  90.0_r8

! Ensure longitudes are within [0, 360.0)
where( longitude <   0.0_r8 ) longitude = longitude + 360.0_r8
where( longitude > 360.0_r8 ) longitude = longitude - 360.0_r8

!  either read existing obs_seq or create a new one

first_obs = .true.
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=output_obs_file, exist=file_exist)

if ( file_exist ) then ! existing file found, append to it

  call read_obs_seq(output_obs_file, 0, 0, n_pix_along_day*n_pix_across_day, obs_seq)

else ! create a new one

  call init_obs_sequence(obs_seq, num_copies, num_qc, n_pix_along_day*n_pix_across_day)
  do i = 1, num_copies
    call set_copy_meta_data(obs_seq, i, 'observation')
  end do
  do i = 1, num_qc
    call set_qc_meta_data(obs_seq, i, 'Quality Control')
  end do

endif

! Set the DART data quality control.  Be consistent with NCEP codes;
! 0 is 'must use', 1 is good, no reason not to use it.

nused = 0

call set_calendar_type(GREGORIAN)

along_day: do i_along = 1, n_pix_along_day

   ! calculate observation time and convert to DART format

   comp_day0 = set_date(year(i_along), 1, 1, 0, 0, 0)
   call get_time(comp_day0, seconds, days)

   seconds  = seconds + nint(time(i_along))
   days     = days + doy(i_along) - 1
   time_obs = set_time(seconds,days)
   call get_time(time_obs, seconds, days)

   if (do_output() .and. (debug > 0)) then   
      write(*,*)'year doy time',year(i_along),doy(i_along),time(i_along)
      call print_time(time_obs,'time of observation')
      call print_date(time_obs,'date of observation')
   endif

   across_day: do i_across = 1, n_pix_across_day
 
      oerr = ON2_uncertainty(i_across,i_along)
      if ( oerr == MISSING_R8 ) cycle across_day
   
      lat =   latitude(i_across,i_along)
      lon =  longitude(i_across,i_along)
      obs_value =  ON2(i_across,i_along)
      qc = ON2_quality(i_across,i_along)

      ! create_3d_obs takes the standard deviation of the observation error
      ! fortunately, this is what is reputed to be in ON2_UNCERTAINTY.

      call create_3d_obs(lat, lon, 0.0_r8, VERTISUNDEF, obs_value, &
                            SSUSI_O_N2_RATIO, oerr, days, seconds, qc, obs)

      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
      nused = nused + 1
   
      if (obs_value >= 0.0_r8 .and. obs_value <= 1.0_r8) then
         if (do_output() .and. (debug > 0)) then   
            write(*,*)'observation ',nused,' out of ',n_pix_along_day*n_pix_across_day
            write(*,*)'  values',lat,lon,qc,obs_value,oerr
         endif
      endif

   enddo across_day
enddo along_day

if (do_output()) then   
   write(*,*)'There were ',nused,' observations converted out of ', &
             n_pix_along_day*n_pix_across_day,' possible.'
endif

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, output_obs_file)

deallocate(latitude, longitude, ON2, ON2_uncertainty, ON2_quality)
deallocate(year, doy, time)

! end of main program
call finalize_utilities()

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
