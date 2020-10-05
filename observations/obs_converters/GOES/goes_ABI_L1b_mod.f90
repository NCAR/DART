! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module goes_ABI_L1b_mod

use types_mod,     only : r4, r8, deg2rad, rad2deg, PI, digits12, &
                          MISSING_R8

use utilities_mod, only : error_handler, E_MSG, E_ERR, &
                          is_longitude_between, register_module

use time_manager_mod, only : time_type, get_date, set_date,            &
                             get_time, set_time, set_calendar_type,    &
                             GREGORIAN, print_date, print_time,        &
                             operator(+)

use obs_sequence_mod, only : init_obs_sequence, init_obs, insert_obs_in_seq, &
                             set_obs_values, obs_sequence_type,              &
                             obs_type, set_copy_meta_data, set_qc_meta_data

use     location_mod, only : location_type, set_location, VERTISUNDEF, &
                             VERTISPRESSURE, get_location

use     obs_kind_mod,  only : get_index_for_type_of_obs

use obs_utilities_mod, only : add_obs_to_seq, create_3d_obs

use obs_def_rttov_mod, only : set_visir_metadata, &
                              get_rttov_option_logical

use netcdf_utilities_mod, only : nc_check

use netcdf

implicit none

private

public :: goes_abi_map_type, goes_load_abi_map, make_obs_sequence

!>@todo FIXME ... we should be using a lot more of the netcdf_utilities interfaces.
!>                they greatly simplify the code.

type goes_abi_map_type
   character(:), allocatable :: filename
   integer :: channel
   integer :: nx
   integer :: ny

   real(r8),    allocatable :: x(:)     
   integer,     allocatable :: x_raw(:) ! GOES fixed grid projection x-coordinate
   real(r4)                 :: x_scale
   real(r4)                 :: x_offset
   real(r8),    allocatable :: y(:) 
   integer,     allocatable :: y_raw(:) ! GOES fixed grid projection y-coordinate
   real(r4)                 :: y_scale
   real(r4)                 :: y_offset
   real(r8),    allocatable :: rad(:,:)
   integer,     allocatable :: rad_raw(:,:) ! radiances, unit (after scale and offset) is mW m-2 sr-1 (cm-1)-1
   real(r4)                 :: rad_scale
   real(r4)                 :: rad_offset
   integer(1),  allocatable :: dqf(:,:) ! data quality flag. 0 = good, 1 is conditionally okay, 2-4 is bad.
   real(r4)                 :: lon_origin
   real(r4)                 :: r_eq  ! radius of earth at equator
   real(r4)                 :: r_pol ! radius of earth at poles
   real(r4)                 :: H ! the total height of satellite 
   real(digits12)           :: t ! mid-point of scan. number of seconds since 2000 01/01 12:00
   real(r4)                 :: sc_lat ! space-craft latitude
   real(r4)                 :: sc_lon ! space-craft longitude
   real(r8),    allocatable :: lat(:,:) 
   real(r8),    allocatable :: lon(:,:)

end type

character(len=512) :: msgstring

integer :: fid, varid

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'goes_ABI_L1b_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical, save :: module_initialized = .false.

contains

subroutine initialize_module

call register_module(source, revision, revdate)

call set_calendar_type(GREGORIAN)

module_initialized = .true.

end subroutine initialize_module

! ----------------------------------------------------------------------------

subroutine goes_load_abi_map(l1b_file, map)

character(len=*),    intent(in)    :: l1b_file
type(goes_abi_map_type), intent(inout) :: map

logical,          parameter :: little_endian = .true.

character(len=*), parameter :: routine = 'goes_load_abi_map'

character(len=512) :: string1

integer :: i, j
real(digits12) :: a, b, c, r_fac, xv, yv, r_s
real(digits12) :: sx, sy, sz
real(digits12) :: time_r8

real(digits12), parameter :: r2d = 180.0_digits12/(atan(1.0_digits12)*4.0_digits12)

if ( .not. module_initialized ) call initialize_module

write(string1,*) 'Now loading L1C data from file ',trim(l1b_file)
call error_handler(E_MSG, routine, string1, source, revision, revdate)

allocate(character(len=len(l1b_file)) :: map%filename)
map%filename = l1b_file

! read x_raw from the file and convert it to x
call nc_check( nf90_open(l1b_file, nf90_nowrite, fid), 'file open', l1b_file)
call nc_check( nf90_inq_dimid(fid, "x", varid), 'inq dimid x', l1b_file)
call nc_check( nf90_inquire_dimension(fid, varid, len=map%nx), 'inq dim x', l1b_file)

allocate(map%x_raw(map%nx))
allocate(map%x(map%nx))

call nc_check( nf90_inq_varid(fid, "x", varid) , 'inq varid x', l1b_file)
call nc_check( nf90_get_var(fid, varid, map%x_raw), 'get var   x', l1b_file)
call nc_check( nf90_get_att(fid, varid, 'scale_factor', map%x_scale) ,'get_att x scale', l1b_file)
call nc_check( nf90_get_att(fid, varid, 'add_offset', map%x_offset) ,'get_att x offset', l1b_file)
map%x = real(map%x_raw)*map%x_scale + map%x_offset

! read y_raw from the file and convert it to y
call nc_check( nf90_inq_dimid(fid, "y", varid), 'inq dimid y', l1b_file)
call nc_check( nf90_inquire_dimension(fid, varid, len=map%ny), 'inq dim y', l1b_file)

allocate(map%y_raw(map%ny))
allocate(map%y(map%ny))

call nc_check( nf90_inq_varid(fid, "y", varid) , 'inq varid y', l1b_file)
call nc_check( nf90_get_var(fid, varid, map%y_raw), 'get var   y', l1b_file)
call nc_check( nf90_get_att(fid, varid, 'scale_factor', map%y_scale) ,'get_att y scale', l1b_file)
call nc_check( nf90_get_att(fid, varid, 'add_offset', map%y_offset) ,'get_att y offset', l1b_file)
map%y = real(map%y_raw)*map%y_scale + map%y_offset

! read rad_raw from the file and convert it to rad
allocate(map%rad_raw(map%nx,map%ny))
allocate(map%rad(map%nx,map%ny))

call nc_check( nf90_inq_varid(fid, "Rad", varid) , 'inq varid Rad', l1b_file)
call nc_check( nf90_get_var(fid, varid, map%Rad_raw), 'get var   Rad', l1b_file)
call nc_check( nf90_get_att(fid, varid, 'scale_factor', map%Rad_scale) ,'get_att Rad scale', l1b_file)
call nc_check( nf90_get_att(fid, varid, 'add_offset', map%Rad_offset) ,'get_att Rad offset', l1b_file)
map%Rad = real(map%Rad_raw)*map%Rad_scale + map%Rad_offset

allocate(map%dqf(map%nx,map%ny))

! read data quality flag (DQF)
call nc_check( nf90_inq_varid(fid, "DQF", varid), 'inq varid DQF', l1b_file)
call nc_check( nf90_get_var(fid, varid, map%DQF), 'get var   DQF', l1b_file)

! read the projection information
call nc_check( nf90_inq_varid(fid, "goes_imager_projection", varid) , 'inq varid goes_imager_projection', l1b_file)
call nc_check( nf90_get_att(fid, varid, 'longitude_of_projection_origin', map%lon_origin) ,'get_att lon_origin', l1b_file)
call nc_check( nf90_get_att(fid, varid, 'semi_major_axis', map%r_eq) ,'get_att semi_major_axis', l1b_file)
call nc_check( nf90_get_att(fid, varid, 'perspective_point_height', map%H) ,'get_att perspective_point_height', l1b_file)
call nc_check( nf90_get_att(fid, varid, 'semi_minor_axis', map%r_pol) ,'get_att semi_minor_axis', l1b_file)

! read the time at the mid-point (currently only this time is used for all observations)
call nc_check( nf90_inq_varid(fid, "t", varid) , 'inq varid t', l1b_file)
call nc_check( nf90_get_var(fid, varid, time_r8), 'get var  t', l1b_file)
map%t = time_r8

! read the channel 
call nc_check( nf90_inq_varid(fid, "band_id", varid) , 'inq varid band_id', l1b_file)
call nc_check( nf90_get_var(fid, varid, map%channel), 'get var  band_id', l1b_file)

! read the satellite latitude
call nc_check( nf90_inq_varid(fid, "nominal_satellite_subpoint_lat", varid) , 'inq varid sat_lat', l1b_file)
call nc_check( nf90_get_var(fid, varid, map%sc_lat), 'get var  sat_lat', l1b_file)

! read the satellite longitude
call nc_check( nf90_inq_varid(fid, "nominal_satellite_subpoint_lon", varid) , 'inq varid sat_lon', l1b_file)
call nc_check( nf90_get_var(fid, varid, map%sc_lon), 'get var  sat_lon', l1b_file)

! close the file as we are now done reading it
call nc_check( nf90_close(fid) , 'close file', l1b_file)

! convert lon_origin from degrees to radians
map%lon_origin = map%lon_origin/r2d

! convert H to total height to earth center
map%H = map%H + map%r_eq

! allocate lat and lon
allocate(map%lat(map%nx,map%ny))
allocate(map%lon(map%nx,map%ny))

map%lat = MISSING_R8
map%lon = MISSING_R8

r_fac = (map%r_eq**2)/(map%r_pol**2)

! now convert the GOES x/y projection to lat/lon pairs
do i=1,map%nx
    xv = map%x(i)
    do j=1,map%ny
        yv = map%y(j)
        ! convert up to digits12 precision here as b**2 - 4*a*c is a bit numerically nasty
        a = sin(xv)**2 + cos(xv)**2 * ( cos(yv)**2 + r_fac * sin(yv)**2 )
        b = -2.0_digits12 * map%H * cos(xv) * cos(yv)
        c = map%H**2 - map%r_eq**2

        if (b**2 >= 4.0_digits12*a*c) then
            r_s = (-b - sqrt(b**2-4.0_digits12*a*c))/(2.0_digits12*a)
            sx  =  r_s * cos(xv)*cos(yv)
            sy  = -r_s * sin(xv)
            sz  =  r_s * cos(xv)*sin(yv)

            map%lat(i,j) = r2d * atan(r_fac*(sz/sqrt((map%H-sx)**2+sy**2)))
            map%lon(i,j) = r2d * (map%lon_origin - atan(sy/(map%H-sx)))
        end if
    end do
end do

end subroutine goes_load_ABI_map

!------------------------------------------------------------------------------
!  extract the ABI channel observations from the map type
!  and convert to DART observation format.  allow caller to specify
!  a bounding box and only extract data within that region.

subroutine make_obs_sequence (seq, map, lon1, lon2, lat1, lat2, &
                              x_thin, y_thin, goes_num, reject_dqf_1, &
                              obs_err_spec, vloc_pres_hPa)

type(obs_sequence_type),    intent(inout) :: seq
type(goes_abi_map_type),    intent(in)    :: map
real(r8),                   intent(in)    :: lon1, lon2, lat1, lat2
integer,                    intent(in)    :: x_thin, y_thin
integer,                    intent(in)    :: goes_num
logical,                    intent(in)    :: reject_dqf_1
real(r8),                   intent(in)    :: obs_err_spec
real(r8),                   intent(in)    :: vloc_pres_hPa

type(obs_type)          :: obs, prev_obs

integer :: ix, iy
integer :: days, seconds
integer :: obs_num, key
integer :: which_vert 

real(r8) :: olon, olat, vloc
real(r8) :: obs_value, obs_err
real(r8) :: rqc
real(r8) :: latd, lond, beta

real(r8) :: lam1, lam2, phi1, phi2, r_fac

real(digits12) :: rdays, remainder

real(r8) :: sat_az, sat_ze, sun_az, sun_ze, specularity
integer :: platform_id, sat_id, sensor_id

type(time_type) :: obs_time, start_time

integer :: robstype
integer :: goes_channel

integer :: num_copies, num_qc
! max possible obs from this one map. in practice if the
! real number of processed channels is very much smaller, make
! another parameter so we don't allocate all these unused obs
! (takes time & space) and then delete them at the end.
integer :: max_num

logical :: is_first_obs
type(time_type) :: pre_time

character(len=512) :: obs_type_name

character(len=*), parameter :: routine = 'make_obs_sequence'

if ( .not. module_initialized ) call initialize_module

! one observation data value and one quality control value
! per obs.  if you change these you have to set additional
! metadata for them below.
num_copies  = 1
num_qc      = 1

! Initialize an obs_sequence
max_num = map%nx*map%ny
call init_obs_sequence(seq, num_copies, num_qc, max_num)

! set meta data of obs_seq
call set_copy_meta_data(seq, 1, 'observation')
call set_qc_meta_data(seq, 1, 'GOES QC')

! Initialize the obs variables
call init_obs(     obs, 1, 1)
call init_obs(prev_obs, 1, 1)

is_first_obs = .true.

! things known (and constant) about the input data and rttov

platform_id = 4   ! GOES

select case(goes_num)
    case (16)
        sat_id = goes_num
        obs_type_name = 'GOES_16_ABI_RADIANCE'
    case (17)
        sat_id = goes_num
        obs_type_name = 'GOES_17_ABI_RADIANCE'
    case (18)
        sat_id = goes_num
        obs_type_name = 'GOES_18_ABI_RADIANCE'
    case (19)
        sat_id = goes_num
        obs_type_name = 'GOES_19_ABI_RADIANCE'
    case default
        write(msgstring,*) 'Unknown GOES number ',goes_num,' should be between 16 and 19'
        call error_handler(E_ERR,routine,msgstring,source,revision,revdate)
end select        

sensor_id   = 44  ! ABI

!------------------------------------------------------------------------------
!  loop over all observations within the file

obs_num = 1
if (vloc_pres_hPa < 0.0_r8) then
    which_vert = VERTISUNDEF
else
    which_vert = VERTISPRESSURE
end if

! assign each observation the correct observation type
robstype = get_index_for_type_of_obs(obs_type_name)
if (robstype < 1) then
   write(msgstring,*) 'unknown observation type ',trim(obs_type_name)
   call error_handler(E_ERR,routine,msgstring,source,revision,revdate)
endif

r_fac = (map%r_eq**2)/(map%r_pol**2)

xloop: do ix=1,map%nx
    ! if we're going to subset x, we will cycle here
    if (x_thin > 0) then
        if (modulo(ix, x_thin) /= 0) cycle xloop
    endif

   ! columns are across-track, varying faster than rows.
   yloop:  do iy=1,map%ny

      ! if we're going to subset y, ditto
      if (y_thin > 0) then
         if (modulo(iy, y_thin) /= 0) cycle yloop
      endif

      ! check channel quality control
      rqc = map%DQF(ix, iy)
      ! reject bad scans here, depending on whether conditional DQF=1 is okay
      if (reject_dqf_1) then
        if (rqc /= 0) cycle yloop 
      else
        if (rqc > 1) cycle yloop
      end if

      ! observation lat, lon:
      olat  = map%lat (ix,iy) ! valid range [ -90.00,  90.00]
      olon  = map%lon (ix,iy) ! valid range [-180.00, 180.00]

      ! verify the location is not outside valid limits.  AIRS  uses -180/180
      if((olon > 180.0_r8) .or. (olon < -180.0_r8) .or.  &
         (olat >  90.0_r8) .or. (olat <  -90.0_r8)) then
         write(*,*)'WARNING : invalid location.  x,y,lon,lat = ', ix,iy,olon,olat
         cycle yloop
      endif

      ! reject observations outside the bounding box (allowing wrapping)
      if(( olat < lat1) .or. ( olat > lat2 ) .or. &
         (.not. is_longitude_between(olon, lon1, lon2))) cycle yloop

      ! make sure lon is between 0 and 360
      if (olon < 0.0_r8) olon = olon + 360.0_r8

      ! set the zenith angle (aka earth incidence angle)
      ! see https://svn.ssec.wisc.edu/repos/cloud_team_cr/trunk/viewing_geometry_module.f90

      latd = (map%lat(ix,iy) - map%sc_lat)*deg2rad
      lond = (map%lon(ix,iy) - map%sc_lon)*deg2rad
      beta = acos(cos(latd)*cos(lond))

      sat_ze = map%H * sin(beta) / sqrt( r_fac + map%H**2 - 2.0_r8*map%r_eq*map%H*cos(beta) )
      sat_ze = max(-1.0_r8,min(1.0_r8, sat_ze))
      sat_ze = asin(sat_ze) / deg2rad

      lam1 = deg2rad*map%lon(ix,iy)*deg2rad
      lam2 = deg2rad*map%sc_lon*deg2rad
      phi1 = deg2rad*map%lat(ix,iy)*deg2rad
      phi2 = deg2rad*map%sc_lat*deg2rad

      ! calculate the bearing between the obs lat/lon and the SClat/lon
      sat_az = atan2(sin(lam2-lam1)*cos(phi2),&
                     cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(lam2-lam1))/deg2rad

      rdays = floor(map%t/86400.0_digits12)
      remainder = map%t - rdays*86400.0_digits12

      !>@todo FIXME : the dart_time_io_mod.f90 or the netcdf_utilities_mod.f90 should
      !> have a routne to read the attribute and set the start_time ....
      ! start_time comes from the 't' variable 'units' attribute
      ! double t ;
      ! t:long_name = "J2000 epoch mid-point between the start and end image scan in seconds" ;
      ! t:standard_name = "time" ;
      ! t:units = "seconds since 2000-01-01 12:00:00" ;
      ! t:axis = "T" ;
      ! t:bounds = "time_bounds" ;

      start_time = set_date(2000,1,1,12,0,0)
      obs_time = set_time(nint(remainder),nint(rdays)) + start_time
      call get_time(obs_time, seconds, days)

      ! create the radiance obs for this observation, add to sequence

      obs_value = map%Rad(ix, iy)
      if (obs_value < 0.0_r8) cycle yloop

      ! TODO: specify this as a function of channel
      obs_err = obs_err_spec

      if (vloc_pres_hPa < 0.0_r8) then
          ! disable pressure location, use VERTISUNDEF, so no single vertical location
          vloc = 0.0_r8
      else
          ! assign the integrated column a vertical location
          ! can correspond to the height of the peak of the weighting function
          vloc = vloc_pres_hPa*100._r8 ! convert from hPa to Pa
      end if

      ! We don't yet have specularity data to add to the observations.
      if (get_rttov_option_logical('do_lambertian')) then
         write(msgstring,*) 'GOES observations do not yet support specularity or Lambertian'
         call error_handler(E_ERR,routine,msgstring,source,revision,revdate)
      else if (get_rttov_option_logical('addsolar')) then
         write(msgstring,*) 'GOES observations do not yet support solar calculations'
         call error_handler(E_ERR,routine,msgstring,source,revision,revdate)
      else
         sun_az = MISSING_R8
         sun_ze = MISSING_R8
         specularity = MISSING_R8
      end if

      ! the RTTOV ABI coefficients start from channel 7
      goes_channel = map%channel-6

      ! add additional metadata for this obs type.  returns key to use in create call
      call set_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, platform_id, sat_id, sensor_id, & 
         goes_channel, specularity)

      call create_3d_obs(olat, olon, vloc, which_vert, obs_value, robstype, &
                             obs_err, days, seconds, rqc, obs, key)

      call add_obs_to_seq(seq, obs, obs_time, prev_obs, pre_time, is_first_obs)

      obs_num = obs_num + 1
   enddo yloop
enddo xloop

!! Print a little summary
!call print_obs_seq(seq, '')

write(msgstring,*) 'Finished loading ',obs_num-1,' of ',key,'total GOES ABI observations for map ' // &
   trim(map%filename)
call error_handler(E_MSG, routine, msgstring, source, revision, revdate)

end subroutine make_obs_sequence

end module
