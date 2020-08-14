! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module gmi_L1_mod

use types_mod,     only : r8, deg2rad, rad2deg, PI, &
                          MISSING_R8
use utilities_mod, only : error_handler, E_MSG, E_ERR, &
                          is_longitude_between, register_module

use time_manager_mod, only : time_type, get_date, set_date,            &
                             get_time, set_time, set_calendar_type,    &
                             GREGORIAN, print_date, print_time

use obs_sequence_mod, only : init_obs_sequence, init_obs, insert_obs_in_seq, &
                             set_obs_values, obs_sequence_type,              &
                             obs_type, set_copy_meta_data, set_qc_meta_data, &
                             print_obs_seq_summary

use     location_mod, only : location_type, set_location, VERTISUNDEF, &
                             get_location

use     obs_kind_mod, only : get_index_for_type_of_obs
use obs_utilities_mod, only : add_obs_to_seq, create_3d_obs
use obs_def_rttov_mod, only : set_mw_metadata, &
                              get_rttov_option_logical

use hdf5

implicit none

private

public :: gmi_swath_type, gmi_l1c_granule_type, &
          gmi_load_l1c_granule, make_obs_sequence, MAXCHANS

type gmi_swath_type
   integer :: nchannels
   integer :: nscans
   integer :: npixels
   integer :: offset     ! offset for total number of channels. S1 = 0, S2 = 9
   character(3) :: dset_prefix

   real(4),    allocatable :: noiseLevel(:)
   real(4),    allocatable :: incidenceAngle(:,:,:)
   integer(1), allocatable :: incidenceAngleIndex(:,:)
   real(4),    allocatable :: Latitude(:,:)
   real(4),    allocatable :: Longitude(:,:)
   integer(1), allocatable :: Quality(:,:)
   integer(1), allocatable :: DayOfMonth(:)
   integer(2), allocatable :: DayOfYear(:)
   integer(1), allocatable :: Hour(:)
   integer(2), allocatable :: MilliSecond(:)
   integer(1), allocatable :: Minute(:)
   integer(1), allocatable :: Month(:)
   integer(1), allocatable :: Second(:)
   real(8),    allocatable :: SecondOfDay(:)
   integer(2), allocatable :: Year(:)
   real(8),    allocatable :: FractionalGranuleNumber(:)
   real(4),    allocatable :: SCaltitude(:)
   real(4),    allocatable :: SClatitude(:)
   real(4),    allocatable :: SClongitude(:)
   integer(2), allocatable :: SCorientation(:)
   integer(1), allocatable :: sunGlintAngle(:,:,:)
   real(4),    allocatable :: Tc(:,:,:)
end type

type gmi_l1c_granule_type
   character(len=:), allocatable :: filename
   type(gmi_swath_type) :: s1
   type(gmi_swath_type) :: s2
end type

integer, parameter :: MAXCHANS = 13

character(len=512) :: msgstring

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'gmi_L1_mod.f90'
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

subroutine gmi_load_l1c_granule(l1c_file, granule)

character(len=*),           intent(in)    :: l1c_file
type(gmi_l1c_granule_type), intent(inout) :: granule

integer :: hdferr

integer(HID_T) :: file_id       ! File identifier

logical, parameter :: little_endian = .true.

character(len=*), parameter :: routine = 'gmi_load_1lc_granule'

character(len=512) :: string1

if ( .not. module_initialized ) call initialize_module

write(string1,*) 'Now loading L1C data from file ',trim(l1c_file)
call error_handler(E_MSG, routine, string1, source, revision, revdate)

allocate(character(len=len(l1c_file)) :: granule%filename)
granule%filename = l1c_file

! initialize hdf if necessary
call h5open_f(hdferr); call h5check(hdferr,'h5open')

! open the file with the given access group
call h5fopen_f(l1c_file, H5F_ACC_RDONLY_F, file_id, hdferr); call h5check(hdferr,'h5fopen')

! set the H5 data set prefixes of the swaths
granule%s1%dset_prefix = '/S1'
granule%s2%dset_prefix = '/S2'

write(string1,*) 'Now loading S1 swath'
call error_handler(E_MSG, routine, string1, source, revision, revdate)
call gmi_load_swath(file_id,granule%s1,little_endian)

write(string1,*) 'Now loading S2 swath'
call error_handler(E_MSG, routine, string1, source, revision, revdate)
call gmi_load_swath(file_id,granule%s2,little_endian)

call h5fclose_f(file_id, hdferr); call h5check(hdferr, 'h5fclose_f')

call h5close_f(hdferr); call h5check(hdferr,'h5close_f')

allocate(granule%s1%noiseLevel(9))

! According to https://www.star.nesdis.noaa.gov/mirs/gpmgmi.php
granule%s1%noiseLevel(1) = 0.96 ! 10.65 V
granule%s1%noiseLevel(2) = 0.96 ! 10.65 H
granule%s1%noiseLevel(3) = 0.84 ! 18.7  V
granule%s1%noiseLevel(4) = 0.84 ! 18.7  H
granule%s1%noiseLevel(5) = 1.05 ! 23.8  V
granule%s1%noiseLevel(6) = 0.65 ! 36.5  V
granule%s1%noiseLevel(7) = 0.65 ! 36.5  H
granule%s1%noiseLevel(8) = 0.57 ! 89.0  V
granule%s1%noiseLevel(9) = 0.57 ! 89.0  H

allocate(granule%s2%noiseLevel(4))

granule%s2%noiseLevel(1) = 1.5 ! 166.0  V
granule%s2%noiseLevel(2) = 1.5 ! 166.0  H
granule%s2%noiseLevel(3) = 1.5 ! 183.31±3  V
granule%s2%noiseLevel(4) = 1.5 ! 183.31±7  V

granule%s1%offset = 0
granule%s2%offset = 9

end subroutine gmi_load_l1c_granule

! ----------------------------------------------------------------------------

subroutine gmi_load_swath(file_id,swath,little_endian)

integer(HID_T),       intent(in)    :: file_id       ! File identifier
type(gmi_swath_type), intent(inout) :: swath
logical,              intent(in)    :: little_endian

integer(HSIZE_T), allocatable :: dims(:)
integer(HSIZE_T) :: dims_ia(3)
integer(HSIZE_T) :: dims_cs(2)
integer(HSIZE_T) :: dims_ps(2)
integer(HSIZE_T) :: dims_s(1)

character(len=3) :: dset_prefix

character(len=*), parameter :: routine = 'gmi_load_swath'

character(len=512) :: string1

dset_prefix = swath%dset_prefix

if ( .not. module_initialized ) call initialize_module

! load the dimensions for the calibrated brightness temperatures
call load_h5_dims(file_id, dset_prefix // '/Tc', dims)

! set the swath metadata
swath%nchannels = dims(1)
swath%npixels   = dims(2)
swath%nscans    = dims(3)

write(string1,*) 'Dimensions are: nchan',swath%nchannels,'npixels:',swath%npixels,&
   'nscans:',swath%nscans
call error_handler(E_MSG, routine, string1, source, revision, revdate)

! load the actual data for Tc
call error_handler(E_MSG, routine, 'Now loading TC', source, revision, revdate)
call load_h5_3F_array(file_id, dset_prefix // '/Tc', dims, swath%Tc, little_endian)

! the incidence angle and sun glint angle have strange dimensions
dims_ia(1) = 1
dims_ia(2) = swath%npixels
dims_ia(3) = swath%nscans

! load the incidence angle and sun glint angle
call error_handler(E_MSG, routine, 'Now loading incidence/sun glint angle', source, revision, revdate)
call load_h5_3F_array(file_id, dset_prefix // '/incidenceAngle', dims_ia, swath%incidenceAngle, little_endian)
call load_h5_3B_array(file_id, dset_prefix // '/sunGlintAngle',  dims_ia, swath%sunGlintAngle,  little_endian)

! incidence angle index is channels x scans
dims_cs(1) = swath%nchannels
dims_cs(2) = swath%nscans

! load the incidence angle index
call error_handler(E_MSG, routine, 'Now loading incidence angle index and quality', source, revision, revdate)
call load_h5_2B_array(file_id, dset_prefix // '/incidenceAngleIndex', dims_cs, swath%incidenceAngleIndex, little_endian)

! lat/lon are pixels x scans
dims_ps(1) = swath%npixels
dims_ps(2) = swath%nscans

! load the lat/lon arrays
call error_handler(E_MSG, routine, 'Now loading lat/lon', source, revision, revdate)
call load_h5_2F_array(file_id, dset_prefix // '/Latitude',  dims_ps, swath%Latitude,  little_endian)
call load_h5_2F_array(file_id, dset_prefix // '/Longitude', dims_ps, swath%Longitude, little_endian)

! load the quality
call error_handler(E_MSG, routine, 'Now loading quality flag', source, revision, revdate)
call load_h5_2B_array(file_id, dset_prefix // '/Quality',   dims_ps, swath%Quality,   little_endian)

! load the scan time variables, which are 1d nscans
dims_s(1) = swath%nscans

call error_handler(E_MSG, routine, 'Now loading scan time values', source, revision, revdate)
call load_h5_1B_array(file_id, dset_prefix // '/ScanTime/DayOfMonth',  dims_s, swath%DayOfMonth,  little_endian)
call load_h5_1S_array(file_id, dset_prefix // '/ScanTime/DayOfYear',   dims_s, swath%DayOfYear,   little_endian)
call load_h5_1B_array(file_id, dset_prefix // '/ScanTime/Hour',        dims_s, swath%Hour,        little_endian)
call load_h5_1S_array(file_id, dset_prefix // '/ScanTime/MilliSecond', dims_s, swath%MilliSecond, little_endian)
call load_h5_1B_array(file_id, dset_prefix // '/ScanTime/Minute',      dims_s, swath%Minute,      little_endian)
call load_h5_1B_array(file_id, dset_prefix // '/ScanTime/Month',       dims_s, swath%Month,       little_endian)
call load_h5_1B_array(file_id, dset_prefix // '/ScanTime/Second',      dims_s, swath%Second,      little_endian)
call load_h5_1D_array(file_id, dset_prefix // '/ScanTime/SecondOfDay', dims_s, swath%SecondOfDay, little_endian)
call load_h5_1S_array(file_id, dset_prefix // '/ScanTime/Year',        dims_s, swath%Year,        little_endian)

! load the space-craft status arrays, which are 1d nscans
call error_handler(E_MSG, routine, 'Now loading SCstatus values', source, revision, revdate)

call load_h5_1D_array(file_id, dset_prefix // '/SCstatus/FractionalGranuleNumber',  dims_s, swath%FractionalGranuleNumber,  little_endian)
call load_h5_1F_array(file_id, dset_prefix // '/SCstatus/SCaltitude',               dims_s, swath%SCaltitude,               little_endian)
call load_h5_1F_array(file_id, dset_prefix // '/SCstatus/SClatitude',               dims_s, swath%SClatitude,               little_endian)
call load_h5_1F_array(file_id, dset_prefix // '/SCstatus/SClongitude',              dims_s, swath%SClongitude,              little_endian)
call load_h5_1S_array(file_id, dset_prefix // '/SCstatus/SCorientation',            dims_s, swath%SCorientation,            little_endian)

end subroutine gmi_load_swath

!------------------------------------------------------------------------------
!  extract the requested GMI channel observations from a granule
!  and convert to DART observation format.  allow caller to specify
!  a bounding box and only extract data within that region.

subroutine make_obs_sequence (seq, granule, lon1, lon2, lat1, lat2, &
                             use_channels, scan_thin, pix_thin)

type(obs_sequence_type),    intent(inout) :: seq
type(gmi_l1c_granule_type), intent(in)    :: granule
real(r8),                   intent(in)    :: lon1, lon2, lat1, lat2
logical,                    intent(in)    :: use_channels(:)
integer,                    intent(in)    :: scan_thin, pix_thin

type(obs_type)          :: obs, prev_obs

integer :: num_copies, num_qc
! max possible obs from this one granule. in practice if the
! real number of processed channels is very much smaller, make
! another parameter so we don't allocate all these unused obs
! (takes time & space) and then delete them at the end.
integer :: max_num

logical :: is_first_obs
type(time_type) :: pre_time

character(len=*), parameter :: routine = 'make_obs_sequence'

if ( .not. module_initialized ) call initialize_module

! one observation data value and one quality control value
! per obs.  if you change these you have to set additional
! metadata for them below.
num_copies  = 1
num_qc      = 1

! Initialize an obs_sequence
max_num = MAXCHANS * granule%s1%npixels * granule%s1%nscans
call init_obs_sequence(seq, num_copies, num_qc, max_num)

! set meta data of obs_seq
call set_copy_meta_data(seq, 1, 'observation')
call set_qc_meta_data(seq, 1, 'GMI QC')

! Initialize the obs variables
call init_obs(     obs, 1, 1)
call init_obs(prev_obs, 1, 1)

is_first_obs = .true.

call add_swath_observations(seq,lon1,lon2,lat1,lat2,use_channels,&
                            scan_thin, pix_thin, obs, prev_obs,  &
                            is_first_obs,pre_time,granule%s1)
call add_swath_observations(seq,lon1,lon2,lat1,lat2,use_channels,&
                            scan_thin, pix_thin, obs, prev_obs,  &
                            is_first_obs,pre_time,granule%s2)

end subroutine make_obs_sequence

!-----------------------------------------------------------------------------

subroutine add_swath_observations(seq, lon1, lon2, lat1, lat2, use_channels,  &
                            scan_thin, pix_thin, obs, prev_obs, is_first_obs, &
                            pre_time, swath)

type(obs_sequence_type),    intent(inout) :: seq
real(r8),                   intent(in)    :: lon1, lon2, lat1, lat2
logical,                    intent(in)    :: use_channels(:)
integer,                    intent(in)    :: scan_thin, pix_thin
type(obs_type),             intent(inout) :: obs, prev_obs
logical,                    intent(inout) :: is_first_obs
type(time_type),            intent(inout) :: pre_time
type(gmi_swath_type),       intent(in)    :: swath

integer :: ipix, iscan, ichan
integer :: days, seconds
integer :: obs_num, key
integer :: which_vert 

real(r8) :: olon, olat, vloc
real(r8) :: obs_value, obs_err
real(r8) :: rqc

real(r8) :: lam1, lam2, phi1, phi2, x, y

real(8) :: sat_az, sat_ze
integer :: platform_id, sat_id, sensor_id
real(8) :: mag_field, cosbk
real(8) :: fastem_p1
real(8) :: fastem_p2
real(8) :: fastem_p3
real(8) :: fastem_p4
real(8) :: fastem_p5

type(time_type) :: obs_time
character(len=*), parameter :: routine = 'add_swath_observations:'

integer :: robstype

! things known (and constant) about the input data and rttov

platform_id = 37  ! GPM
sat_id      = 1   ! Only only satellite
sensor_id   = 71  ! GMI 

!------------------------------------------------------------------------------
!  loop over all observations within the file

obs_num = 1
which_vert = VERTISUNDEF

! assign each observation the correct observation type
robstype = get_index_for_type_of_obs('GPM_1_GMI_TB')
if (robstype < 1) then
   msgstring = 'unknown observation type GPM_1_GMI_TB'
   call error_handler(E_ERR,routine,msgstring,source,revision,revdate)
endif

! rows are along-track, stepping in the direction the satellite is moving
scanloop:  do iscan=1,swath%nscans

   ! if we're going to subset rows, we will cycle here
   if (scan_thin > 0) then
      if (modulo(iscan, scan_thin) /= 0) cycle scanloop
   endif

   ! columns are across-track, varying faster than rows.
   pixloop:  do ipix=1,swath%npixels

      ! if we're going to subset columns, ditto
      if (pix_thin > 0) then
         if (modulo(ipix, pix_thin) /= 0) cycle pixloop
      endif

      ! check channel quality control
      rqc = swath%Quality(ipix, iscan)
      ! reject bad scans here - cycle pixloop
      if (rqc /= 0) cycle pixloop

      ! observation lat, lon:
      olat  = swath%Latitude (ipix,iscan) ! valid range [ -90.00,  90.00]
      olon  = swath%Longitude(ipix,iscan) ! valid range [-180.00, 180.00]

      ! verify the location is not outside valid limits.  AIRS  uses -180/180
      if((olon > 180.0_r8) .or. (olon < -180.0_r8) .or.  &
         (olat >  90.0_r8) .or. (olat <  -90.0_r8)) then
         write(*,*)'WARNING : invalid location.  col,row,lon,lat = ', ipix,iscan,olon,olat
         cycle pixloop
      endif

      ! reject observations outside the bounding box (allowing wrapping)
      if(( olat < lat1) .or. ( olat > lat2 ) .or. &
         (.not. is_longitude_between(olon, lon1, lon2))) cycle pixloop

      ! make sure lon is between 0 and 360
      if (olon < 0.0_r8) olon = olon + 360.0_r8

      ! set the zenith angle (aka earth incidence angle)
      sat_ze = swath%incidenceAngle(1,ipix,iscan)

      lam1 = deg2rad*swath%longitude(ipix,iscan)
      lam2 = deg2rad*swath%SClongitude(iscan)
      phi1 = deg2rad*swath%latitude(ipix,iscan)
      phi2 = deg2rad*swath%SClatitude(iscan)

      ! calculate the bearing between the obs lat/lon and the SClat/lon
      y = sin(lam2-lam1)*cos(phi2)
      x = cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(lam2-lam1)
      sat_az = rad2deg*atan2(y,x)

      obs_time = set_date(int(swath%year(iscan)),int(swath%month(iscan)), &
         int(swath%DayOfMonth(iscan)), int(swath%Hour(iscan)),            &
         int(swath%Minute(iscan)), int(swath%Second(iscan)))

      call get_time(obs_time, seconds, days)

      channel_loop: do ichan=1, swath%nchannels

         if (.not. use_channels(ichan+swath%offset)) cycle channel_loop

         ! create the radiance obs for this observation, add to sequence

         ! apparently -9999 is missing data, outside of qc mechanism
         obs_value = swath%Tc(ichan, ipix, iscan)
         if (obs_value < 0.0_r8) cycle channel_loop

         obs_err = swath%noiseLevel(ichan)

         ! column integrated value, so no vertical location
         vloc = 0.0_r8

         ! We don't yet have specularity data to add to the observations.
         if (get_rttov_option_logical('use_zeeman')) then
            write(msgstring,*) 'GMI observations do not yet support the Zeeman effect'
            call error_handler(E_ERR,routine,msgstring,source,revision,revdate)
         else
            mag_field = MISSING_R8
            cosbk = MISSING_R8
         end if

         ! FIXME: consider adding an atlas for looking up FASTEM parameters
         ! From the RTTOV User guide:
         ! Surface type
         ! Typical RTTOV default for land: 3.0, 5.0, 15.0, 0.1, 0.3
         ! Summer land surface:
         !        Forest: 1.7, 1.0, 163.0, 0.0, 0.5
         !        Open grass: 2.2, 1.3, 138.0, 0.0, 0.42
         !        Bare soil: 2.3, 1.9, 21.8, 0.0, 0.5
         ! Winter surface type:
         !        Forest and snow: 2.9, 3.4, 27.0, 0.0, 0.0
         !        Deep dry snow: 3.0, 24.0, 60.0, 0.1, 0.15
         !        Frozen soil: 117.8, 2.0, 0.19, 0.2 ,0.35
         ! Sea ice
         !        Grease ice: 23.7, 7.7, 17.3, 0.0, 0.15
         !        Baltic nilas: 1.6, 3.3, 2.2, 0.0, 0.0
         !        New ice (no snow): 2.9, 3.4, 27.0, 0.0, 0.0
         !        New ice (snow): 2.2, 3.7, 122.0, 0.0, 0.15
         !        Brash ice: 3.0, 5.5, 183.0, 0.0, 0.0
         !        Compact pack ice: 2.0, 1700000.0, 49000000.0, 0.0, 0.0
         !        Fast ice: 1.5, 77.8, 703.0, 0.1, 0.35
         !        Lake ice + snow: 1.8, 67.1, 534.0, 0.1, 0.15
         !        Multi-year ice: 1.5, 85000.0, 4700000.0, 0.0, 0.0
         fastem_p1 = 3.0d0
         fastem_p2 = 5.0d0
         fastem_p3 = 15.0d0
         fastem_p4 = 0.1d0
         fastem_p5 = 0.3d0

         ! add additional metadata for this obs type.  returns key to use in create call
         call set_mw_metadata(key, sat_az, sat_ze, platform_id, sat_id, sensor_id, & 
            ichan+swath%offset, mag_field, cosbk, fastem_p1, fastem_p2, fastem_p3, &
            fastem_p4, fastem_p5)

         call create_3d_obs(olat, olon, vloc, which_vert, obs_value, robstype, &
                            obs_err, days, seconds, rqc, obs, key)

         call add_obs_to_seq(seq, obs, obs_time, prev_obs, pre_time, is_first_obs)

         obs_num = obs_num + 1

      enddo channel_loop
   enddo pixloop
enddo scanloop

! Print a little summary
call print_obs_seq_summary(seq)

write(msgstring,*) 'Converted ',obs_num-1,' obs for swath ',swath%dset_prefix, &
                   '; total GMI obs = ',key
call error_handler(E_MSG, routine, msgstring, source, revision, revdate)

end subroutine add_swath_observations

! ----------------------------------------------------------------------------

subroutine load_h5_dims(file_id, location_in_file, dims)

integer(HID_T),                intent(in)  :: file_id       ! File identifier
character(len=*),              intent(in)  :: location_in_file
integer(HSIZE_T), allocatable, intent(out) :: dims(:)

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: dspace_id     ! Dataspace identifier
integer        :: ndims

integer :: hdferr

integer(HSIZE_T), allocatable :: maxdims(:)

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! open the dataspace for the data set
call h5dget_space_f(dset_id,dspace_id,hdferr); call h5check(hdferr,'h5dget_space_f')

! get the number of dimensions for the data space
call h5sget_simple_extent_ndims_f(dspace_id, ndims, hdferr)
call h5check(hdferr,'h5sget_simple_extent_ndims')

! allocate a place to store the dimensions
allocate(dims(ndims),maxdims(ndims))

! get the actual dimensions
call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr)
call h5check(hdferr,'h5sget_simple_extent_dims_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr)

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_3D_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),       intent(in)  :: file_id          ! File identifier
character(len=*),     intent(in)  :: location_in_file
integer(HSIZE_T),     intent(in)  :: dims(3)
real(8), allocatable, intent(out) :: array(:,:,:)
logical,              intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1),dims(2),dims(3)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_IEEE_F64LE
else
   mem_type = H5T_IEEE_F64BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_3F_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),       intent(in)  :: file_id          ! File identifier
character(len=*),     intent(in)  :: location_in_file
integer(HSIZE_T),     intent(in)  :: dims(3)
real(4), allocatable, intent(out) :: array(:,:,:)
logical,              intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1),dims(2),dims(3)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_IEEE_F32LE
else
   mem_type = H5T_IEEE_F32BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_3L_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),          intent(in)  :: file_id          ! File identifier
character(len=*),        intent(in)  :: location_in_file
integer(HSIZE_T),        intent(in)  :: dims(3)
integer(8), allocatable, intent(out) :: array(:,:,:)
logical,                 intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1),dims(2),dims(3)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_STD_I64LE
else
   mem_type = H5T_STD_I64BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_3I_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),          intent(in)  :: file_id          ! File identifier
character(len=*),        intent(in)  :: location_in_file
integer(HSIZE_T),        intent(in)  :: dims(3)
integer(4), allocatable, intent(out) :: array(:,:,:)
logical,                 intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1),dims(2),dims(3)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_STD_I32LE
else
   mem_type = H5T_STD_I32BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_3S_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),          intent(in)  :: file_id          ! File identifier
character(len=*),        intent(in)  :: location_in_file
integer(HSIZE_T),        intent(in)  :: dims(3)
integer(2), allocatable, intent(out) :: array(:,:,:)
logical,                 intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1),dims(2),dims(3)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_STD_I16LE
else
   mem_type = H5T_STD_I16BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_3B_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),          intent(in)  :: file_id          ! File identifier
character(len=*),        intent(in)  :: location_in_file
integer(HSIZE_T),        intent(in)  :: dims(3)
integer(1), allocatable, intent(out) :: array(:,:,:)
logical,                 intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1),dims(2),dims(3)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_STD_I8LE
else
   mem_type = H5T_STD_I8BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_2D_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),       intent(in)  :: file_id          ! File identifier
character(len=*),     intent(in)  :: location_in_file
integer(HSIZE_T),     intent(in)  :: dims(2)
real(8), allocatable, intent(out) :: array(:,:)
logical,              intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1),dims(2)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_IEEE_F64LE
else
   mem_type = H5T_IEEE_F64BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_2F_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),       intent(in)  :: file_id          ! File identifier
character(len=*),     intent(in)  :: location_in_file
integer(HSIZE_T),     intent(in)  :: dims(2)
real(4), allocatable, intent(out) :: array(:,:)
logical,              intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1),dims(2)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_IEEE_F32LE
else
   mem_type = H5T_IEEE_F32BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_2L_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),          intent(in)  :: file_id          ! File identifier
character(len=*),        intent(in)  :: location_in_file
integer(HSIZE_T),        intent(in)  :: dims(2)
integer(8), allocatable, intent(out) :: array(:,:)
logical,                 intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1),dims(2)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_STD_I64LE
else
   mem_type = H5T_STD_I64BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_2I_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),          intent(in)  :: file_id          ! File identifier
character(len=*),        intent(in)  :: location_in_file
integer(HSIZE_T),        intent(in)  :: dims(2)
integer(4), allocatable, intent(out) :: array(:,:)
logical,                 intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1),dims(2)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_STD_I32LE
else
   mem_type = H5T_STD_I32BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_2S_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),          intent(in)  :: file_id          ! File identifier
character(len=*),        intent(in)  :: location_in_file
integer(HSIZE_T),        intent(in)  :: dims(2)
integer(2), allocatable, intent(out) :: array(:,:)
logical,                 intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1),dims(2)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_STD_I16LE
else
   mem_type = H5T_STD_I16BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_2B_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),          intent(in)  :: file_id          ! File identifier
character(len=*),        intent(in)  :: location_in_file
integer(HSIZE_T),        intent(in)  :: dims(2)
integer(1), allocatable, intent(out) :: array(:,:)
logical,                 intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1),dims(2)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_STD_I8LE
else
   mem_type = H5T_STD_I8BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_1D_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),       intent(in)  :: file_id          ! File identifier
character(len=*),     intent(in)  :: location_in_file
integer(HSIZE_T),     intent(in)  :: dims(1)
real(8), allocatable, intent(out) :: array(:)
logical,              intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_IEEE_F64LE
else
   mem_type = H5T_IEEE_F64BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_1F_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),       intent(in)  :: file_id          ! File identifier
character(len=*),     intent(in)  :: location_in_file
integer(HSIZE_T),     intent(in)  :: dims(1)
real(4), allocatable, intent(out) :: array(:)
logical,              intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_IEEE_F32LE
else
   mem_type = H5T_IEEE_F32BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_1L_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),          intent(in)  :: file_id          ! File identifier
character(len=*),        intent(in)  :: location_in_file
integer(HSIZE_T),        intent(in)  :: dims(1)
integer(8), allocatable, intent(out) :: array(:)
logical,                 intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_STD_I64LE
else
   mem_type = H5T_STD_I64BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_1I_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),          intent(in)  :: file_id          ! File identifier
character(len=*),        intent(in)  :: location_in_file
integer(HSIZE_T),        intent(in)  :: dims(1)
integer(4), allocatable, intent(out) :: array(:)
logical,                 intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_STD_I32LE
else
   mem_type = H5T_STD_I32BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_1S_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),          intent(in)  :: file_id          ! File identifier
character(len=*),        intent(in)  :: location_in_file
integer(HSIZE_T),        intent(in)  :: dims(1)
integer(2), allocatable, intent(out) :: array(:)
logical,                 intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_STD_I16LE
else
   mem_type = H5T_STD_I16BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine load_h5_1B_array(file_id, location_in_file, dims, array, little_endian)

integer(HID_T),          intent(in)  :: file_id          ! File identifier
character(len=*),        intent(in)  :: location_in_file
integer(HSIZE_T),        intent(in)  :: dims(1)
integer(1), allocatable, intent(out) :: array(:)
logical,                 intent(in)  :: little_endian

integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: mem_type      ! Memory type (int, float, etc. and little/big endian)

integer :: hdferr

! allocate the array
allocate(array(dims(1)))

! open the data set
call h5dopen_f(file_id, location_in_file, dset_id, hdferr); call h5check(hdferr,'h5dopen_f')

! choose the memory type based on endianness
if (little_endian) then
   mem_type = H5T_STD_I8LE
else
   mem_type = H5T_STD_I8BE
end if

! read the data into the array
call h5dread_f(dset_id, mem_type, array, dims, hdferr); call h5check(hdferr, 'h5dread_f')

! close the data set
call h5dclose_f(dset_id, hdferr); call h5check(hdferr, 'h5dclose_f')

end subroutine

! ----------------------------------------------------------------------------

subroutine h5check(hdferr, message)

implicit none

integer,                    intent(in) :: hdferr
character(len=*), optional, intent(in) :: message

character(len=*), parameter :: routine = 'h5check'

character(len=512) :: string1

if (hdferr .lt. 0) then
   if (present(message)) then
      write(string1,*) trim(message),', hdf error code: ',hdferr
      call error_handler(E_ERR, routine, string1, source, revision, revdate)
   else
      write(string1,*) 'HDF Error code:',hdferr
      call error_handler(E_ERR, routine, string1, source, revision, revdate)
   end if
end if

end subroutine h5check

! ----------------------------------------------------------------------------


end module
