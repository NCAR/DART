module rttov_interface_mod

use     location_mod,  only : location_type, set_location, get_location, VERTISUNDEF, &
                              VERTISHEIGHT, VERTISLEVEL, set_location_missing, &
                              is_vertical, &
                              VERTISHEIGHT

use    utilities_mod,  only : register_module, error_handler, E_ERR, E_MSG, E_WARN, &
                              nmlfileunit, check_namelist_read,      &
                              find_namelist_in_file, do_nml_file, do_nml_term, &
                              ascii_file_format

! This code contains the DART to RTTOV interface. It uses a data structure of
! platforms/satellites/sensors to store runtime information for each sensor.
! This runtime information is initialized on the fly as needed to compute the
! RTTOV forward operator. 
!
! Copyright:
!    RTTOV was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    RTTOV is copyright 2017, EUMETSAT, All Rights Reserved.

  use types_mod, only : r8

  ! RTTOV module containing useful RTTOV constants
  use rttov_const, only :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name,           &
         pmax,                &
         pmin

  ! RTTOV module rttov_types contains definitions of all RTTOV data types
  use rttov_types, only :     &
         rttov_options,       &
         rttov_options_scatt, &
         rttov_scatt_coef,    &
         rttov_coefs,         &
         rttov_profile,       &
         rttov_profile_cloud, &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  use parkind1, only : jpim, jprb, jplm
  
  implicit none

! include the interface files as per the RTTOV standard
#include "rttov_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_read_scattcoeffs.interface"
#include "rttov_dealloc_scattcoeffs.interface"
#include "rttov_scatt_setupindex.interface"
#include "rttov_scatt.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_alloc_scatt_prof.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_opts_scatt.interface"
#include "rttov_print_profile.interface"
#include "rttov_print_cld_profile.interface"
#include "rttov_skipcommentline.interface"

!--------------------------
!
integer(kind=jpim), parameter :: ioout = 0    ! stdout for now

! Whether to output (copious amounts of) debug information
logical :: debug = .false.

! Expose publically the DART/RTTOV types and method calls
public :: visir_metadata_type,               &
          mw_metadata_type,                  &
          atmos_profile_type,                &
          aerosol_profile_type,              &
          cloud_profile_type,                &
          rttov_sensor_type,                 &
          rttov_sensor_runtime_type,         &
          rttov_satellite_type,              &
          rttov_platform_type,               &
          get_rttov_sensor,                  &
          read_sensor_db_file,               &
          sensor_runtime_setup,              &
          do_forward_model,                  &
          sensor_runtime_takedown

! Metadata for rttov observations.

type visir_metadata_type
   real(r8) :: sat_az      ! azimuth of satellite position (degrees)
   real(r8) :: sat_ze      ! zenith of satellite position (degrees)
   real(r8) :: sun_az      ! azimuth of solar position (degrees, only used with add_solar)
   real(r8) :: sun_ze      ! zenith of solar position  (degrees, only used with add_solar)
   integer  :: platform_id ! see rttov user guide, table 2
   integer  :: sat_id      ! see rttov user guide, table 2
   integer  :: sensor_id   ! see rttov user guide, table 3
   integer  :: channel     ! each channel is a different obs
   real(r8) :: specularity ! specularity (0-1, only used with do_lambertian)
end type visir_metadata_type

type mw_metadata_type
   real(r8) :: sat_az      ! azimuth of satellite position (degrees)
   real(r8) :: sat_ze      ! zenith of satellite position (degrees)
   integer  :: platform_id ! see rttov user guide, table 2
   integer  :: sat_id      ! see rttov user guide, table 2
   integer  :: sensor_id   ! see rttov user guide, table 3
   integer  :: channel     !  each channel is a different obs
   real(r8) :: mag_field   ! strength of mag_field (Gauss, )
   real(r8) :: cosbk       ! cosine of angle between mag field and viewing angle
   real(r8) :: fastem_land1 ! FASTEM land/sea ice parameter 1
   real(r8) :: fastem_land2 ! FASTEM land/sea ice parameter 2
   real(r8) :: fastem_land3 ! FASTEM land/sea ice parameter 3
   real(r8) :: fastem_land4 ! FASTEM land/sea ice parameter 4
   real(r8) :: fastem_land5 ! FASTEM land/sea ice parameter 5
end type mw_metadata_type

! DART container type to hold the essential atmosphere and surface fields.
! For 2D fields, the order is (ens_size, numlevels), while 1D fields are 
! (ens_size).
type atmos_profile_type
   real(r8), allocatable :: temperature(:,:)    ! mandatory, level temperature (K)
   real(r8), allocatable :: pressure(:,:)       ! mandatory, level pressure (hPa)
   real(r8), allocatable :: moisture(:,:)       ! mandatory, level water vapor (kg/kg)
   real(r8), allocatable :: sfc_p(:)            ! mandatory, surface pressure (hPa)
   real(r8), allocatable :: s2m_t(:)            ! mandatory, 2 meter temp (K)
   real(r8), allocatable :: skin_temp(:)        ! mandatory, surface skin temp (K)
   real(r8), allocatable :: sfc_elev(:)         ! mandatory, surface elevation (km)
   real(r8), allocatable :: surftype(:)         ! mandatory, surface type (land=0, sea=1, seaice = 2) 
   real(r8), allocatable :: s2m_q(:)            ! optional, 2 meter wator vapor (kg/kg) (used if add_q2m)
   real(r8), allocatable :: s10m_u(:)           ! optional, 10 meter u wind (m/s) (used if add_uv10m)
   real(r8), allocatable :: s10m_v(:)           ! optional, 10 meter v wind (m/s) (used if add_uv10m)
   real(r8), allocatable :: wfetch(:)           ! optional, wind fetch (m) (used if use_wfetch)
   real(r8), allocatable :: water_type(:)       ! optional, water type (fresh=0, ocean=1) (used if use_water_type)
   real(r8), allocatable :: sfc_salinity(:)     ! optional, ocean salinity (practial salinity unit) (used if use_salinity)
   real(r8), allocatable :: sfc_foam_frac(:)    ! optional, foam fraction (0-1) (used if supply_foam_fraction)
   real(r8), allocatable :: sfc_snow_frac(:)    ! optional, snow cover (0-1) (used if use_sfc_snow_frac)
end type atmos_profile_type

! container type for trace gasses. Note these values are on levels, so the 
! size of the arrays are (ens_size, numlevels)
type trace_gas_profile_type
   real(r8), allocatable :: ozone(:,:)         ! ozone concentration (kg/kg) (used if add_ozone)
   real(r8), allocatable :: co2(:,:)           ! CO2 concentration   (kg/kg) (used if add_co2)
   real(r8), allocatable :: n2o(:,:)           ! N2O concentration   (kg/kg) (used if add_n2o)
   real(r8), allocatable :: ch4(:,:)           ! CH4 concentration   (kg/kg) (used if add_ch4)
   real(r8), allocatable :: co(:,:)            ! CO concentration    (kg/kg) (used if add_co)
   real(r8), allocatable :: so2(:,:)           ! SO2 concentration   (kg/kg) (used if add_so2)
end type trace_gas_profile_type

! Container type for aerosols. Note these values are on layers, so the 
! size of the arrays are (ens_size, numlevels-1)
type aerosol_profile_type
   real(r8), allocatable :: insoluble(:,:)                  ! INSO (kg/kg), OPAC only
   real(r8), allocatable :: water_soluble(:,:)             ! WASO, OPAC
   real(r8), allocatable :: soot(:,:)                       ! SOOT, OPAC
   real(r8), allocatable :: sea_salt_accum(:,:)             ! SSAM, OPAC
   real(r8), allocatable :: sea_salt_coarse(:,:)            ! SSCM, OPAC
   real(r8), allocatable :: mineral_nucleus(:,:)            ! MINM, OPAC
   real(r8), allocatable :: mineral_accum(:,:)              ! MIAM, OPAC
   real(r8), allocatable :: mineral_coarse(:,:)             ! MICM, OPAC
   real(r8), allocatable :: mineral_transport(:,:)          ! MITR, OPAC
   real(r8), allocatable :: sulphated_droplets(:,:)         ! SUSO, OPAC
   real(r8), allocatable :: volcanic_ash(:,:)               ! VOLA, OPAC
   real(r8), allocatable :: new_volcanic_ash(:,:)           ! VAPO, OPAC
   real(r8), allocatable :: asian_dust(:,:)                 ! ASDU, OPAC
   real(r8), allocatable :: black_carbon(:,:)               ! BCAR, CAMS
   real(r8), allocatable :: dust_bin1(:,:)                  ! DUS1, CAMS
   real(r8), allocatable :: dust_bin2(:,:)                  ! DUS2, CAMS
   real(r8), allocatable :: dust_bin3(:,:)                  ! DUS3, CAMS
   real(r8), allocatable :: ammonium_sulphate(:,:)          ! SULP, CAMS
   real(r8), allocatable :: sea_salt_bin1(:,:)              ! SSA1, CAMS
   real(r8), allocatable :: sea_salt_bin2(:,:)              ! SSA2, CAMS
   real(r8), allocatable :: sea_salt_bin3(:,:)              ! SSA3, CAMS
   real(r8), allocatable :: hydrophilic_organic_matter(:,:) ! OMAT, CAMS
end type aerosol_profile_type

! container type for clouds - note RTTOV uses different cloud fields depending on scheme, frequency, type, etc. Note these values are on layers, not levels, so the 
! size of the arrays are (ens_size, numlevels-1)
type cloud_profile_type
   real(r8), allocatable :: cfrac(:,:)         ! (VIS/IR/MW) cloud fractional cover (0-1) 
   real(r8), allocatable :: simple_cfrac(:)    ! (VIS/IR)    cloud fraction for simple cloud (0-1) 
   real(r8), allocatable :: ctp(:)             ! (VIS/IR)    cloud top pressure for simple cloud (hPa) 
   real(r8), allocatable :: w(:,:)             ! (VIS/IR)    vertical velocity (used for classification of cumulus vs. stratus)
   real(r8), allocatable :: clw(:,:)           ! (VIS/IR/MW) cloud non-precipitating water
   real(r8), allocatable :: rain(:,:)          ! (VIS/IR/MW) cloud precipitating water
   real(r8), allocatable :: ciw(:,:)           ! (VIS/IR/MW) cloud non-precipitating ice concentration
   real(r8), allocatable :: snow(:,:)          ! (VIS/IR/MW) cloud precipitating ice (fluffy)
   real(r8), allocatable :: graupel(:,:)       ! (VIS/IR/MW) cloud precipitating ice (soft hail / snow pellets)
   real(r8), allocatable :: hail(:,:)          ! (VIS/IR/MW) cloud precipitating ice (hard hail)
   real(r8), allocatable :: clwde(:,:)         ! (VIS/IR)    cloud liquid effective diameter if clw_scheme = 2
   real(r8), allocatable :: icede(:,:)         ! (VIS/IR)    cloud liquid effective diameter if ice_scheme = 1
end type cloud_profile_type

! Container for the DART/rttov run-time structures to be used (per sensor)
! This is setup per sensor in accordance with Figure 1 of the RTTOV user guide,
! which indicates the allocation is per coefficient file, which is
! per sensor. RTTOV_alloc_direct also takes in the coefficient file. It's
! not clear if some of these structures could be reused across different
! sensors - to be safe, they are gathered together here. There will be 
! one runtime per sensor, but the runtime will only be allocated when the 
! sensor is used by a forward operator.
type rttov_sensor_runtime_type
   type(rttov_options),       pointer :: opts               => null() ! Options for RTTOV-DIRECT 
   type(rttov_options_scatt), pointer :: opts_scatt         => null() ! Options for RTTOV-SCATT
   type(rttov_coefs)                  :: coefs                        ! Sensor coefficients structure
   type(rttov_scatt_coef),    pointer :: coefs_scatt        => null() ! Sensor coefficients structure
   type(rttov_chanprof),      pointer :: chanprof(:)        => null() ! Input channel/profile list
   integer(kind=jpim),        pointer :: frequencies(:)     => null() ! Indexes into Mietable lookup
   integer(kind=jpim),        pointer :: frequencies_all(:) => null() ! Indexes into Mietable lookup
   logical(kind=jplm),        pointer :: calcemis(:)        => null() ! Flag to indicate calculation of emissivity within RTTOV
   type(rttov_emissivity),    pointer :: emissivity(:)      => null() ! Input/output surface emissivity
   logical(kind=jplm),        pointer :: calcrefl(:)        => null() ! Flag to indicate calculation of BRDF within RTTOV
   type(rttov_reflectance),   pointer :: reflectance(:)     => null() ! Input/output surface BRDF
   type(rttov_profile),       pointer :: profiles(:)        => null() ! Input profiles
   type(rttov_profile_cloud), pointer :: cld_profiles(:)    => null() ! Input profiles
   type(rttov_transmission)           :: transmission                 ! Output transmittances
   type(rttov_radiance)               :: radiance                     ! Output radiances
end type rttov_sensor_runtime_type

! a simple fortran type containing the generic sensor information 
type rttov_sensor_type
   integer              :: platform_id
   integer              :: satellite_id
   integer              :: sensor_id
   character(len=512)   :: obs_type
   character(len=3)     :: sensor_type_name ! mw, ir, or vis
   integer              :: sensor_type_num  ! vis = 1, ir = 2, mw = 3
   character(len=512)   :: coefficient_file
   integer, allocatable :: channels(:)      ! list of channels to use

   ! Runtime per-sensor RTTOV structures used to execute the forward operator
   type(rttov_sensor_runtime_type), pointer :: runtime => null()

   ! Linked-list of rttov_sensor_types, sorted (asc) by sensor_id 
   type(rttov_sensor_type),         pointer :: next_sensor  => null()
end type rttov_sensor_type

type rttov_satellite_type
   ! this is a singly-linked sorted (asc by sensor_id) list of sensors on the satellite
   type(rttov_sensor_type), pointer :: head
end type rttov_satellite_type

type rttov_platform_type
   ! this is a trie - each array index corresponds to the satellite number
   type(rttov_satellite_type), pointer :: sats(:)
end type rttov_platform_type

! this is a trie - each array index corresponds to the platform number
type(rttov_platform_type), pointer :: platforms(:)

integer, parameter :: NUM_PLATFORMS_INITIAL  = 60 ! first guess of # of platforms
integer, parameter :: NUM_SATELLITES_INITIAL = 25 ! first guess of # of satellites per platform

! arrays of length nlevels, initialized in do_forward_model
integer,  allocatable :: lvlidx(:)
integer,  allocatable :: ly1idx(:)
integer,  allocatable :: ly2idx(:)
real(r8), allocatable :: totalwater(:)
real(r8), allocatable :: totalice(:)

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

contains

!----------------------------------------------------------------------
! Add a sensor to the platform/sat/sensor data structures with the 
! given information. Note that this will create the platform, satellite,
! and/or sensors as necessary.
! 
subroutine add_sensor(platform_id, satellite_id, sensor_id, sensor_type_name, obs_type, &
   coef_file, channels)

   integer,          intent(in) :: platform_id        ! the RTTOV platform id
   integer,          intent(in) :: satellite_id       ! the RTTOV satellite id
   integer,          intent(in) :: sensor_id          ! the RTTOV sensor id
   character(len=*), intent(in) :: sensor_type_name   ! the type of sensor (vis, ir, mw)
   character(len=*), intent(in) :: obs_type           ! the DART observation type
   character(len=*), intent(in) :: coef_file          ! the RTTOV coefficient file
   integer,          intent(in) :: channels(:)        ! the channels to simulate (or 0-length for all)

   integer :: i, new_size

   type(rttov_platform_type),  pointer :: tmp_platforms(:)   ! temporary platform array, for dynamic growing
   type(rttov_satellite_type), pointer :: tmp_satellites(:)  ! temporary satellite array, for dynamic growing

   type(rttov_platform_type),  pointer :: platform          
   type(rttov_satellite_type), pointer :: satellite         
   type(rttov_sensor_type),    pointer :: sensor           
   type(rttov_sensor_type),    pointer :: new_sensor
   type(rttov_sensor_type),    pointer :: last_sensor

   character(len=*), parameter :: routine = 'add_sensor'

   character(len=512) :: string1
   character(len=512) :: string2

   string2 = 'name: ' // trim(obs_type)

   ! check 
   if (platform_id <= 0 .or. satellite_id < 0 .or. sensor_id < 0) then
      write(string1,*)'Invalid platform/satellite/sensor id:',&
         platform_id, ',',satellite_id,',',sensor_id
      call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
   end if

   ! allocate the platforms trie if necessary
   if (.not. associated(platforms)) then
      allocate(platforms(NUM_PLATFORMS_INITIAL))

      do i=1,size(platforms)
         ! by default, leave the contained sats trie null
         platforms(i) % sats => null()
      end do

      if (debug) then
         write(string1,*)'allocated platforms:',size(platforms)
         call error_handler(E_MSG, routine, string1, source, revision, revdate)
      end if
   end if

   ! check if we need to expand platforms array
   ! can't combine this with the above as platform_id may be too large at different times
   if (platform_id > size(platforms)) then

      new_size = 2*size(platforms)

      if (platform_id > new_size) then
         write(string1,*) 'Error: platform ID is too large: ',platform_id,'vs.',new_size
         call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
      else
         write(string1,*)'Resizing platform size from ',size(platforms),'to',&
            new_size
         call error_handler(E_WARN, routine, string1, source, revision, revdate)

         ! increase the size of platforms array
         allocate(tmp_platforms(new_size))

         do i=1,size(platforms)
            ! copy the pointer to the sats array
            tmp_platforms(i) % sats => platforms(i) % sats 
         end do

         do i=size(platforms)+1,new_size
            tmp_platforms(i) % sats => null()
         end do

         ! make the switch. deallocate the platform array
         deallocate(platforms)
         platforms => tmp_platforms
      end if 
   end if

   platform => platforms(platform_id)

   ! allocate the sats trie if necessary
   if (.not. associated(platform % sats)) then
      ! there are satellites with id 0, so start from 0 index
      allocate(platform % sats(0:NUM_SATELLITES_INITIAL-1))

      ! initialize the sensor linked list to null
      do i=0,size(platform % sats)-1
         platform % sats(i) % head => null()
      end do

      if (debug) then
         write(string1,*) 'allocated platform sats:',platform_id,size(platform % sats)
         call error_handler(E_MSG, routine, string1, source, revision, revdate)
      end if
   end if

   ! check if we would need to expand the sats array
   if (satellite_id > size(platform % sats)-1) then

      new_size = 2*size(platform % sats)

      if (satellite_id > new_size-1) then
         write(string1,*) 'Error: satellite ID is too large: ',satellite_id,'vs.',new_size
         call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
      else
         write(string1,*) 'Resizing platform % sats size from ',size(platform % sats),'to',&
            new_size
         call error_handler(E_WARN, routine, string1, source, revision, revdate)

         allocate(tmp_satellites(0:new_size))

         ! copy the pointers for the linked list over
         do i=0,size(platform % sats)-1
            tmp_satellites(i) % head => platform % sats(i) % head
         end do

         do i=size(platform % sats),new_size-1
            tmp_satellites(i) % head => null()
         end do

         ! make the switch
         deallocate(platform % sats)
         platform % sats => tmp_satellites
      end if 
   end if

   ! get the satellite pointer from the sats trie
   satellite => platform % sats(satellite_id)

   ! setup our new sensor structure to be added
   allocate(new_sensor)
   new_sensor % platform_id  = platform_id
   new_sensor % satellite_id = satellite_id
   new_sensor % sensor_id    = sensor_id
   new_sensor % obs_type     = obs_type
   new_sensor % sensor_type_name = sensor_type_name

   if (trim(sensor_type_name) == 'vis') then
      new_sensor % sensor_type_num = 1
   elseif (trim(sensor_type_name) == 'ir') then
      new_sensor % sensor_type_num = 2
   elseif (trim(sensor_type_name) == 'mw') then
      new_sensor % sensor_type_num = 3
   elseif (trim(sensor_type_name) == 'po') then
      new_sensor % sensor_type_num = 4
   elseif (trim(sensor_type_name) == 'hi') then
      new_sensor % sensor_type_num = 5
   else
      write(string1,*) 'Error: unknown sensor_type_name: ', trim(sensor_type_name),':',&
         platform_id,satellite_id,sensor_id
      call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
   end if

   new_sensor % coefficient_file = coef_file
   allocate(new_sensor % channels(size(channels)))
   new_sensor % channels(:) = channels(:)
   new_sensor % next_sensor => null()

   ! traverse the sensor list. Insert the node so the list is sorted.
   last_sensor => null()
   sensor => satellite % head

   linkedlist: do
      if (.not. associated(sensor)) then
         ! reached the end of the list. Add new_sensor to the end
         if (.not. associated(last_sensor)) then
            ! update head => new_sensor
            satellite % head => new_sensor
         else
            ! update last_sensor % next => sensor
            last_sensor % next_sensor => new_sensor
         end if

         exit ! the linked list loop
      end if

      ! check for duplicates
      if (sensor_id == sensor % sensor_id) then
         write(string1,*) 'Error: tried to add the same sensor to a satellite twice:',&
            platform_id,satellite_id,sensor_id
         call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
      end if

      ! if the sensor id to add is greater than the current one, we found our spot:
      ! we need to put the new sensor in between sensor and last_sensor
      if (sensor_id > sensor % sensor_id) then

         ! update new_sensor % next => sensor
         new_sensor % next_sensor => sensor

         if (associated(last_sensor)) then
            ! update last_sensor % next => new_sensor % next => sensor
            last_sensor % next_sensor => new_sensor
         else
            ! update head => new_sensor % next => sensor
            satellite % head => new_sensor
         end if

         exit ! the linked list loop
      else
         ! traverse forward in the list and keep looking
         last_sensor => sensor
         sensor => sensor % next_sensor
      end if
   end do linkedlist
end subroutine add_sensor

!----------------------------------------------------------------------
! Read the sensor DB file and setup the relevant data structures.
! This function will read through the DB file twice: once to 
! discover how many lines are in the file and the maximum line
! length, and a second time to actually parse the contents of the
! file. The contents of the file are assumed to be an "unformatted"
! (aka text) CSV file with the following contents:
! 
! <obs_type>,<platform_id>,<satellite_id>,<sensor_id>,<sensor_type>,<coef_file>,[channel list]
!
! where: 
!     obs_type     is the DART observation quantity,
!     platform_id  is the RTTOV platform id
!     satellite_id is the RTTOV satellite id
!     sensor_id    is the RTTOV sensor id
!     sensor_type  is vis, ir, or mw (visible, infrared, or microwave)
!     coef_file    is the RTTOV coefficient file
!     channel_list is an optional list of channels (can be zero-length, meaning all available channels should be used) 
! 
subroutine read_sensor_db_file(sensor_db_file)

character(len=*), intent(in) :: sensor_db_file

! a buffer to use for reading the file piece by piece
integer, parameter :: buffer_len = 80
character(len=buffer_len) :: buffer
character(len=2048) :: next_line

integer :: nlines, io, size_read, line_len, maxline_len
integer :: i
character(len=100) :: imsg

! a random fortran io unit
integer, parameter :: dbUnit = 193

character(len=256) :: obs_type
integer :: platform_id
integer :: satellite_id
integer :: sensor_id
integer :: ntok, toknum, lastind
integer, allocatable :: channels(:)
character(len=16)   :: sensor_type
character(len=512)  :: coef_file
character(len=512)  :: token
! assume the file is delimited by commas
character(len=*), parameter :: delimiter = ','

character(len=*), parameter :: routine = 'read_sensor_db_file'

character(len=512) :: string1

if (debug) then
   write(string1,*)'Now reading sensor_db_file: ',adjustl(trim(sensor_db_file))
   call error_handler(E_MSG, routine, string1, source, revision, revdate)
end if

open(unit=dbUnit, file=sensor_db_file)

! throw away the header row
read(dbUnit,'(A)') buffer

! Now count the number of lines in the file.
! While we are at it, find the size of a string large enough to hold
! the longest line for error checking
nlines = 0
maxline_len = -1

testlineloop: do 
   line_len = 0
   ! the next loop finds the length of each line
   bufferloop: do 
      read(dbUnit,'(A)',iostat=io,iomsg=imsg,advance='no',size=size_read) buffer

      if (io == 0) then
         ! 0 = no error, read the full buffer, continue
         line_len = line_len + buffer_len
      elseif (io == -1) then
         ! -1 = end of file reached. Exit both loops.
         exit ! the buffer loop
      elseif (io <= 0) then
         ! Different compilers have different end of record iostat values
         ! All return a negative value not equal to -1, however.
         line_len = line_len + size_read
         nlines = nlines + 1
         maxline_len = max(maxline_len,line_len)
         exit ! the buffer loop
      else
         ! for any other error, quit
         write(string1,*)'fatal error:',io,imsg
         call error_handler(E_ERR, routine, string1, source, revision, revdate)
      end if
   end do bufferloop

   if (io == -1) then
      ! all compilers have -1 = end of file
      exit ! the line loop
   end if
end do testlineloop

if (debug) then
   write(string1,*)'found nlines:',nlines,'max len:',maxline_len
   call error_handler(E_MSG, routine, string1, source, revision, revdate)
end if

if (maxline_len > len(next_line)) then
   write(string1,*)'The length of a line in the ' // trim(sensor_db_file) &
      //' file exceeds the maximum: ',maxline_len,'vs.',len(next_line)
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
end if

close(dbUnit)

open(unit=dbUnit, file=sensor_db_file)
! throw away the header again
read(dbUnit,'(A)') buffer

lineloop: do
   read(dbUnit,'(A)', iostat=io) next_line
   if (io == 0) then
      if (debug) then
         write(string1,*)'read the line:',trim(next_line)
         call error_handler(E_MSG, routine, string1, source, revision, revdate)
      end if
   else
      exit ! the line loop
   end if

   toknum = 0
   lastind = 1

   ntok = 0

   ! first pass: count the number of tokens
   do i=1,len(next_line)
      if (next_line(i:i) == delimiter .or. i==len(next_line)) then
         ntok = ntok + 1
      end if
   end do

   if (debug) then
      write(string1,*)'found ntokens:',ntok
      call error_handler(E_MSG, routine, string1, source, revision, revdate)
   end if

   allocate(channels(ntok-6))

   ! second pass: parse the tokens
   do i=1,len(next_line)
      if (next_line(i:i) == delimiter .or. i == len(next_line)) then
         toknum = toknum + 1
         token = adjustl(trim(next_line(lastind:i-1)))
         if (debug) then
            write(string1,*)'found the token :',toknum,trim(token)
            call error_handler(E_MSG, routine, string1, source, revision, revdate)
         end if
         select case (toknum)
            case (1)
               obs_type = token
            case (2)
               platform_id = str2int(token)
            case (3)
               satellite_id = str2int(token)
            case (4)
               sensor_id = str2int(token)
            case (5)
               sensor_type = token
            case (6)
               coef_file = token
            case default
               channels(toknum-6) = str2int(token)
         end select

         lastind = i + 1
      end if
   end do

   if (debug) then
      write(string1,*)'values:',trim(obs_type),platform_id,satellite_id,sensor_id,&
         trim(sensor_type),trim(coef_file),channels
      call error_handler(E_MSG, routine, string1, source, revision, revdate)
   end if

   call add_sensor(platform_id, satellite_id, sensor_id, sensor_type, obs_type, &
      coef_file, channels)

   deallocate(channels)
end do lineloop

end subroutine read_sensor_db_file

!----------------------------------------------------------------------
! Function to return the sensor associated with the given platform/sat/
! sensor id. If no sensor is found or if invalid id combinations are 
! passed in, an error will be generated.
! 
function get_rttov_sensor(instrument_id) result(sensor)

integer, intent(in) :: instrument_id(3)

integer :: platform_id
integer :: satellite_id
integer :: sensor_id

type(rttov_platform_type),  pointer :: platform
type(rttov_satellite_type), pointer :: satellite
type(rttov_sensor_type),    pointer :: sensor

character(len=*), parameter :: routine = 'get_rttov_sensor'

character(len=512) :: string1

platform_id  = instrument_id(1)
satellite_id = instrument_id(2)
sensor_id    = instrument_id(3)

if (platform_id <= 0 .or. platform_id > size(platforms) .or. &
    satellite_id < 0 .or. sensor_id < 0) then
   write(string1,*)'Invalid platform/satellite/sensor id:', platform_id, ',',&
      satellite_id,',',sensor_id
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
end if

platform => platforms(platform_id)

if (satellite_id > size(platform % sats)-1) then
   write(string1,*)'Invalid satellite id:', satellite_id,size(platform % sats)-1
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
end if

satellite => platform % sats(satellite_id)

sensor => satellite % head

linkedlist: do
   ! reached the end of the list, but didn't find what we were looking for.
   if (.not. associated(sensor)) then
      write(string1,*)'Could not find the sensor with the platform/satellite/sensor id:',&
         platform_id,satellite_id,sensor_id
      call error_handler(E_ERR, routine, string1, source, revision, revdate)
   end if

   ! check for our id
   if (sensor_id == sensor % sensor_id) then
      ! found it - return without error (function result is sensor)
      return
   end if

   ! if the sensor id is greater than the current one, we are out of luck
   if (sensor_id > sensor % sensor_id) then
      write(string1,*)'Could not find the sensor with the platform/satellite/sensor id:',&
         platform_id,satellite_id,sensor_id
      call error_handler(E_ERR, routine, string1, source, revision, revdate)
   else
      ! moving on to the next sensor in the list
      sensor => sensor % next_sensor 
   end if
end do linkedlist 

end function get_rttov_sensor

!----------------------------------------------------------------------
! Clean up (deallocate and nullify) the platforms/satellite/sensors data 
! structures.
! 
subroutine clean_up_sensors()

type(rttov_platform_type),  pointer :: platform
type(rttov_satellite_type), pointer :: satellite
type(rttov_sensor_type),    pointer :: sensor
type(rttov_sensor_type),    pointer :: next_sensor

integer :: i, j, error_status

do i=1,size(platforms)
   platform => platforms(i)

   do j=1,size(platform % sats)
      satellite => platform % sats(j)

      sensor => satellite % head

      do while (associated(sensor))
         next_sensor => sensor % next_sensor

         if (allocated(sensor % channels)) then
            deallocate(sensor % channels)
         end if

         if (associated(sensor % runtime)) then
            call sensor_runtime_takedown(sensor % runtime, error_status) 
            deallocate(sensor % runtime)
         end if

         deallocate(sensor)

         sensor => next_sensor
      end do

      nullify(satellite % head)
   end do

   deallocate(platform % sats)
   nullify(platform % sats)
end do

deallocate(platforms)
nullify(platforms)
end subroutine clean_up_sensors

!----------------------------------------------------------------------
! A helper function to convert a character(len=*) into an integer
!
function str2int(str) result(intv)

character(len=*), intent(in) :: str

integer :: strlength
integer :: intv

! the format for reading in the integer
character(len=4) :: strfmt

strlength = len(adjustl(trim(str)))
if (strlength > 9) then
   print *,'Error: integer string length is greater than 9 digits long.'
   stop
end if 

! write the string format, e.g. I6 for a 6 digit string
write(strfmt,'(A,I1,A)') '(I', strlength,')'

read(str,strfmt) intv

end function str2int

!----------------------------------------------------------------------
! Setup the sensor runtime structures for use with RTTOV in DART. Note
! that this loads the RTTOV coefficient file, allocates memory, sets the
! RTTOV options, and provides some basic-level of error checking.
! 
subroutine sensor_runtime_setup(sensor, ens_size, nlevs, opts, opts_scatt, &
   use_totalice)

type(rttov_sensor_type),             intent(inout)  :: sensor
integer,                             intent(in)     :: ens_size
integer,                             intent(in)     :: nlevs
type(rttov_options),                 pointer        :: opts
type(rttov_options_scatt), optional, pointer        :: opts_scatt
logical,                   optional, intent(in)     :: use_totalice ! only relevant if opts_scatt is present

character(len=*), parameter :: routine = 'sensor_runtime_setup'

character(len=512) :: string1
character(len=512) :: string2 ! used to hold the sensor % obs_type, don't overwrite

type(rttov_sensor_runtime_type), pointer :: runtime

integer :: nchannels

! Return error status of RTTOV subroutine calls
integer(kind=jpim)               :: errorstatus 

logical :: do_totalice

integer :: instrument(3) ! platform_id, satellite_id, sensor_id
integer :: alloc_status
integer :: i, ich, nch

logical(kind=jplm), pointer :: use_chan(:,:)  => null() ! flags to specify channels to simulate

! used below, don't overwrite
string2 = 'name: ' // trim(sensor % obs_type)

instrument(1) = sensor % platform_id
instrument(2) = sensor % satellite_id
instrument(3) = sensor % sensor_id

if (associated(sensor % runtime)) then
   write(string1,*)'Runtime information already initialized for platform/sat/sensor id combination:',&
      instrument
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
end if

if (debug) then
   write(string1,*) 'Now initializing the platform/sat/sensor id combination:',&
      instrument
   call error_handler(E_MSG, routine, string1, source, revision, revdate, text2=string2)
end if

allocate(runtime)
sensor % runtime => runtime

! --------------------------------------------------------------------------
! 1. Set options 
! --------------------------------------------------------------------------

runtime % opts => opts

if (present(opts_scatt)) then
   runtime % opts_scatt => opts_scatt
end if

if (present(use_totalice)) then
   do_totalice = use_totalice
else
   do_totalice = .false.
end if

! --------------------------------------------------------------------------
! 2. Read coefficients
! --------------------------------------------------------------------------

if (debug) then
   write(string1,*)'The coefficient file is:',trim(sensor % coefficient_file)
   call error_handler(E_MSG, routine, string1, source, revision, revdate)
end if

call rttov_read_coefs(errorstatus, runtime % coefs, runtime % opts, instrument=instrument)

if (errorstatus /= errorstatus_success) then
   write(string1,*)'Error reading coefs for platform/sat/sensor id combination:',&
      instrument
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
endif

if (present(opts_scatt)) then
   allocate(runtime % coefs_scatt)
   call rttov_read_scattcoeffs(errorstatus, runtime % opts_scatt, runtime % coefs, &
      runtime % coefs_scatt)
   if (errorstatus /= errorstatus_success) then
      write(string1,*) 'Error reading scatt coefs for platform/sat/sensor id combination:',&
         instrument
      call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
   endif
end if

if (size(sensor % channels) /= 0) then
   nchannels = size(sensor % channels)
else
   nchannels = runtime % coefs % coef % fmv_chn
end if

! Ensure input number of channels is not higher than number stored in coefficient file
IF (nchannels > runtime % coefs % coef % fmv_chn) THEN
  write(string1,*) 'Number of requested channels too large: ',nchannels,&
      ' vs. ',runtime % coefs % coef % fmv_chn
  call error_handler(E_ERR, routine, string1, source, revision, revdate)
endif

! Ensure the options and coefficients are consistent
if (runtime % opts % config % do_checkinput) then
   if (debug) then
      write(string1,*) 'Now running rttov_user_options_checkinput'
      call error_handler(E_MSG, routine, string1, source, revision, revdate)
   end if

   call rttov_user_options_checkinput(errorstatus, runtime % opts, runtime % coefs)

   if (errorstatus /= errorstatus_success) then
     write(string1,*) 'error in rttov_user_options_checkinput'
     call error_handler(E_ERR, routine, string1, source, revision, revdate)
   end if
end if

! --------------------------------------------------------------------------
! 3. Allocate RTTOV input and output structures
! --------------------------------------------------------------------------

! In DART RTTOV_Direct, we only simulate a single channel at a time for all ens members
allocate(runtime % profiles(ens_size))

! Allocate structures for rttov_direct
CALL rttov_alloc_direct(                  &
      errorstatus,                        &
      1_jpim,                             &  ! 1 => allocate
      ens_size,                           &
      ens_size,                           &
      nlevs,                              &
      runtime % chanprof,                 &
      runtime % opts,                     &
      runtime % profiles,                 &
      runtime % coefs,                    &
      runtime % transmission,             &
      runtime % radiance,                 &
      calcemis    =runtime % calcemis,    &
      emissivity  =runtime % emissivity,  &
      calcrefl    =runtime % calcrefl,    &
      reflectance =runtime % reflectance, &
      init=.TRUE._jplm)

IF (errorstatus /= errorstatus_success) THEN
   write(string1,*) 'allocation error for rttov_direct structures'
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
endif

if (debug) then
   write(string1,*) 'Successfully initialized rttov_direct platform/sat/sensor id combination:',&
      instrument
   call error_handler(E_MSG, routine, string1, source, revision, revdate, text2=string2)
end if

if (present(opts_scatt)) then
   ! Allocate cld_profile for rttov_scatt
   allocate(runtime % cld_profiles(ens_size), stat=alloc_status)
   if (alloc_status /= 0) then
      write(string1,*) 'allocation error for cld_profiles array'
      call error_handler(E_ERR, routine, string1, source, revision, revdate)
   end if

   CALL rttov_alloc_scatt_prof(               &
         err=errorstatus,                     &
         nprof=ens_size,                      &
         cld_profiles=runtime % cld_profiles, &
         nlev=nlevs,                          &
         use_totalice=do_totalice,            &
         asw=1,                               & ! 1 = allocate
         init=.TRUE._jplm,                    &
         mmr_snowrain=.TRUE._jplm)              ! true = kg/kg units for clouds

   if (errorstatus /= errorstatus_success) then
      write(string1,*) 'allocation error for rttov_direct structures'
      call error_handler(E_ERR, routine, string1, source, revision, revdate)
   endif

   allocate(use_chan(ens_size,runtime % coefs % coef % fmv_chn))
   allocate(runtime % frequencies(ens_size))

   ! only the channels to simulate will be set
   allocate(runtime % frequencies_all(runtime % coefs % coef % fmv_chn))

   if (size(sensor % channels) /= 0) then
      nch = size(sensor % channels)
   else
      nch = runtime % coefs % coef % fmv_chn
   end if

   do i=1,ens_size
      runtime % chanprof(i) % prof = i
   end do

   do i=1,nch
      use_chan(:,:) = .FALSE._jplm

      if (size(sensor % channels) /= 0) then
         ich = sensor % channels(i)
      else
         ich = i
      end if

      ! Set use_chan to .TRUE. only for the ith required channel
      use_chan(:,ich) = .TRUE._jplm

      runtime % chanprof(:) % chan = ich

      ! Populate chanprof and frequencies arrays for this one channel
      call rttov_scatt_setupindex (                    &
            nprofiles=ens_size,                        &
            n_chan=runtime % coefs % coef % fmv_chn,   &
            coef_rttov=runtime % coefs,                &
            nchannels=ens_size,                        &
            chanprof=runtime % chanprof,               &
            frequencies=runtime % frequencies,         &
            lchannel_subset=use_chan )

      runtime % frequencies_all(ich) = runtime % frequencies(1)
   end do

   if (debug) then
      write(string1,*) 'Successfully initialized RTTOV-scatt cloud profiles for platform/sat/sensor id combination:',&
         instrument
      call error_handler(E_MSG, routine, string1, source, revision, revdate, text2=string2)
   end if
end if

end subroutine sensor_runtime_setup

subroutine atmos_profile_setup(atmos, ens_size, numlevels, use_q2m, &
      use_uv10m, use_wfetch, use_water_type, use_salinity,          &
      supply_foam_fraction,  use_sfc_snow_frac)

type(atmos_profile_type), intent(inout) :: atmos
integer,                  intent(in)    :: ens_size
integer,                  intent(in)    :: numlevels
logical,                  intent(in)    :: use_q2m
logical,                  intent(in)    :: use_uv10m
logical,                  intent(in)    :: use_wfetch
logical,                  intent(in)    :: use_water_type
logical,                  intent(in)    :: use_salinity
logical,                  intent(in)    :: supply_foam_fraction
logical,                  intent(in)    :: use_sfc_snow_frac

allocate(atmos%temperature(ens_size, numlevels), &
         atmos%   moisture(ens_size, numlevels), &
         atmos%   pressure(ens_size, numlevels), &
         atmos%      sfc_p(ens_size),          &
         atmos%      s2m_t(ens_size),          &
         atmos%  skin_temp(ens_size),          &
         atmos%   sfc_elev(ens_size),          &
         atmos%   surftype(ens_size))

! zero the arrays as well
atmos%temperature = 0.d0
atmos%   moisture = 0.d0
atmos%   pressure = 0.d0
atmos%      sfc_p = 0.d0
atmos%      s2m_t = 0.d0
atmos%  skin_temp = 0.d0
atmos%   sfc_elev = 0.d0
atmos%   surftype = 0.d0

if (use_q2m) then
   allocate(atmos%s2m_q(ens_size))
   atmos%s2m_q = 0.d0
end if

if (use_uv10m) then
   allocate(atmos%s10m_u(ens_size))
   allocate(atmos%s10m_v(ens_size))

   atmos%s10m_u = 0.d0
   atmos%s10m_v = 0.d0
end if

if (use_wfetch) then
   allocate(atmos%wfetch(ens_size))

   atmos%wfetch = 0.d0
end if

if (use_water_type) then
   allocate(atmos%water_type(ens_size))

   atmos%water_type = -1
end if

if (use_salinity) then
   allocate(atmos%sfc_salinity(ens_size))

   atmos%sfc_salinity = 0.d0
end if

if (supply_foam_fraction) then
   allocate(atmos%sfc_foam_frac(ens_size))
   atmos%sfc_foam_frac = 0.d0
end if

if (use_sfc_snow_frac) then
   allocate(atmos%sfc_snow_frac(ens_size))
   atmos%sfc_snow_frac = 0.d0
end if

end subroutine atmos_profile_setup

subroutine trace_gas_profile_setup(trace_gas, ens_size, numlevels,  &
      ozone_data, co2_data, n2o_data, ch4_data, co_data)

type(trace_gas_profile_type), intent(inout) :: trace_gas
integer,                      intent(in)    :: ens_size
integer,                      intent(in)    :: numlevels
logical,                      intent(in)    :: ozone_data
logical,                      intent(in)    :: co2_data
logical,                      intent(in)    :: n2o_data
logical,                      intent(in)    :: ch4_data
logical,                      intent(in)    :: co_data

if (ozone_data) then
   allocate(trace_gas%ozone(ens_size, numlevels))
   trace_gas%ozone = 0.d0
end if

if (co2_data) then
   allocate(trace_gas%co2(ens_size, numlevels))
   trace_gas%co2 = 0.d0
end if

if (n2o_data) then
   allocate(trace_gas%n2o(ens_size, numlevels))
   trace_gas%n2o = 0.d0
end if

if (ch4_data) then
   allocate(trace_gas%ch4(ens_size, numlevels))
   trace_gas%ch4 = 0.d0
end if

if (co_data) then
   allocate(trace_gas%co(ens_size, numlevels))
   trace_gas%co = 0.d0
end if

end subroutine trace_gas_profile_setup


subroutine aerosol_profile_setup(aerosols, ens_size, numlevels,  &
      aerosl_type)

type(aerosol_profile_type), intent(inout) :: aerosols
integer,                    intent(in)    :: ens_size
integer,                    intent(in)    :: numlevels
integer,                    intent(in)    :: aerosl_type

character(len=512) :: string1
character(len=*), parameter :: routine = 'aerosol_profile_setup'

if (aerosl_type == 1) then
   ! OPAC
   allocate(aerosols%insoluble(ens_size,numlevels)) 
   allocate(aerosols%water_soluble(ens_size,numlevels)) 
   allocate(aerosols%soot(ens_size,numlevels)) 
   allocate(aerosols%sea_salt_accum(ens_size,numlevels)) 
   allocate(aerosols%sea_salt_coarse(ens_size,numlevels)) 
   allocate(aerosols%mineral_nucleus(ens_size,numlevels)) 
   allocate(aerosols%mineral_accum(ens_size,numlevels)) 
   allocate(aerosols%mineral_coarse(ens_size,numlevels)) 
   allocate(aerosols%mineral_transport(ens_size,numlevels)) 
   allocate(aerosols%sulphated_droplets(ens_size,numlevels)) 
   allocate(aerosols%volcanic_ash(ens_size,numlevels)) 
   allocate(aerosols%new_volcanic_ash(ens_size,numlevels)) 
   allocate(aerosols%asian_dust(ens_size,numlevels)) 

   aerosols%insoluble = 0.d0
   aerosols%water_soluble = 0.d0
   aerosols%soot = 0.d0
   aerosols%sea_salt_accum = 0.d0
   aerosols%sea_salt_coarse = 0.d0
   aerosols%mineral_nucleus = 0.d0
   aerosols%mineral_accum = 0.d0
   aerosols%mineral_coarse = 0.d0
   aerosols%mineral_transport = 0.d0
   aerosols%sulphated_droplets = 0.d0
   aerosols%volcanic_ash = 0.d0
   aerosols%new_volcanic_ash = 0.d0
   aerosols%asian_dust = 0.d0
elseif (aerosl_type == 2) then
   ! CAMS
   allocate(aerosols%black_carbon(ens_size,numlevels)) 
   allocate(aerosols%dust_bin1(ens_size,numlevels)) 
   allocate(aerosols%dust_bin2(ens_size,numlevels)) 
   allocate(aerosols%dust_bin3(ens_size,numlevels)) 
   allocate(aerosols%ammonium_sulphate(ens_size,numlevels)) 
   allocate(aerosols%sea_salt_bin1(ens_size,numlevels)) 
   allocate(aerosols%sea_salt_bin2(ens_size,numlevels)) 
   allocate(aerosols%sea_salt_bin3(ens_size,numlevels)) 
   allocate(aerosols%hydrophilic_organic_matter(ens_size,numlevels)) 

   aerosols%black_carbon = 0.d0
   aerosols%dust_bin1 = 0.d0
   aerosols%dust_bin2 = 0.d0
   aerosols%dust_bin3 = 0.d0
   aerosols%ammonium_sulphate = 0.d0
   aerosols%sea_salt_bin1 = 0.d0
   aerosols%sea_salt_bin2 = 0.d0
   aerosols%sea_salt_bin3 = 0.d0
   aerosols%hydrophilic_organic_matter = 0.d0
else
   ! error
   write(string1,*)"Unknown aerosl_type:",aerosl_type
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
end if

end subroutine aerosol_profile_setup

subroutine cloud_profile_setup(clouds, ens_size, numlevels,  &
      cfrac_data, clw_data, clw_scheme, rain_data, ciw_data, &
      ice_scheme, use_icede, snow_data, graupel_data,        &
      hail_data, w_data, htfrtc_simple_cloud)

type(cloud_profile_type), intent(inout) :: clouds
integer,                  intent(in)    :: ens_size
integer,                  intent(in)    :: numlevels
logical,                  intent(in)    :: cfrac_data
logical,                  intent(in)    :: clw_data
integer,                  intent(in)    :: clw_scheme
logical,                  intent(in)    :: rain_data
logical,                  intent(in)    :: ciw_data
integer,                  intent(in)    :: ice_scheme
logical,                  intent(in)    :: use_icede
logical,                  intent(in)    :: snow_data
logical,                  intent(in)    :: graupel_data
logical,                  intent(in)    :: hail_data
logical,                  intent(in)    :: w_data
logical,                  intent(in)    :: htfrtc_simple_cloud

! RTTOV wants layers, but models probably prefer levels
if (cfrac_data) then
   allocate(clouds%cfrac(ens_size, numlevels))
   clouds%cfrac = 0.d0
end if

if (clw_data) then
   allocate(clouds%clw(ens_size, numlevels))
   clouds%clw = 0.d0

   if (clw_scheme == 2) then
      allocate(clouds%clwde(ens_size, numlevels)) 
      clouds%clwde = 0.d0
   end if
end if

if (rain_data) then
   allocate(clouds%rain(ens_size, numlevels))
   clouds%rain = 0.d0
end if

if (ciw_data) then
   allocate(clouds%ciw(ens_size, numlevels))
   clouds%ciw = 0.d0
   if (ice_scheme == 1 .and. use_icede) then
      allocate(clouds%icede(ens_size, numlevels))
      clouds%icede = 0.d0
   end if
end if

if (snow_data) then
   allocate(clouds%snow(ens_size, numlevels))
   clouds%snow = 0.d0
end if

if (graupel_data) then
   allocate(clouds%graupel(ens_size, numlevels))
   clouds%graupel = 0.d0
end if

if (hail_data) then
   allocate(clouds%hail(ens_size, numlevels))
   clouds%hail = 0.d0
end if

if (w_data) then
   allocate(clouds%w(ens_size, numlevels))
   clouds%w = 0.d0
end if

if (htfrtc_simple_cloud) then
   allocate(clouds%simple_cfrac(ens_size))
   allocate(clouds%ctp(ens_size))
   clouds%simple_cfrac = 0.d0
   clouds%ctp = 0.d0
end if

end subroutine cloud_profile_setup

!----------------------------------------------------------------------
! Run the forward model for the given sensor, preallocating arrays and
! initializing the sensor if necessary (which will be on the first call
! to this operator and the first time using the sensor, respectively).

subroutine do_forward_model(ens_size, nlevels, location,   &
   atmos, trace_gas, clouds, aerosols, sensor, channel,    &
   first_lvl_is_sfc, mw_clear_sky_only, clw_scheme,        &
   ice_scheme, idg_scheme, aerosl_type, do_lambertian,     &
   use_totalice, use_zeeman, use_fastem_params, radiances, &
   error_status, visir_md, mw_md)

integer,                                intent(in)  :: ens_size
integer,                                intent(in)  :: nlevels
type(location_type),                    intent(in)  :: location
type(atmos_profile_type),               intent(in)  :: atmos
type(trace_gas_profile_type),           intent(in)  :: trace_gas
type(cloud_profile_type),               intent(in)  :: clouds
type(aerosol_profile_type),             intent(in)  :: aerosols
type(rttov_sensor_type),                pointer     :: sensor
integer,                                intent(in)  :: channel
logical,                                intent(in)  :: first_lvl_is_sfc
logical,                                intent(in)  :: mw_clear_sky_only
integer,                                intent(in)  :: clw_scheme
integer,                                intent(in)  :: ice_scheme
integer,                                intent(in)  :: idg_scheme
integer,                                intent(in)  :: aerosl_type
logical,                                intent(in)  :: do_lambertian
logical,                                intent(in)  :: use_totalice
logical,                                intent(in)  :: use_zeeman
logical,                                intent(in)  :: use_fastem_params
real(r8),                               intent(out) :: radiances(ens_size)
integer,                                intent(out) :: error_status(ens_size)
type(visir_metadata_type),     pointer, intent(in)  :: visir_md
type(mw_metadata_type),        pointer, intent(in)  :: mw_md

integer :: ilvl, imem 

! observation location variables
real(r8) :: lon, lat, obsloc(3)

integer(kind=jpim) :: j, nch

! Return error status of RTTOV subroutine calls
integer(kind=jpim)               :: errorstatus 

type(rttov_sensor_runtime_type), pointer :: runtime

character(len=512) :: string1
character(len=512) :: string2

character(len=*), parameter :: routine = 'do_forward_model'

real(r8) :: maxw

logical :: is_visir
logical :: is_mw
logical :: is_cumulus
integer :: instrument(3)
integer :: surftype

if (.not. associated(sensor)) then
   write(string1,*)'Passed an unassociated sensor'
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
end if

error_status(:) = 0 ! 0 is success

string2 = 'name: ' // trim(sensor % obs_type)

instrument(1) = sensor % platform_id
instrument(2) = sensor % satellite_id
instrument(3) = sensor % sensor_id

is_visir = associated(visir_md)
is_mw    = associated(mw_md)

if (.not. is_visir .and. .not. is_mw) then
   write(string1,*)'Neither vis/ir nor mw metadata were present for platform/sat/sensor id combination:',&
      sensor % platform_id,sensor % satellite_id,sensor % sensor_id
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
end if

if (is_visir .and. is_mw) then
   write(string1,*)'Both vis/ir and mw metadata were present (only one can be specified) for platform/sat/sensor id combination:',&
      sensor % platform_id,sensor % satellite_id,sensor % sensor_id
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
end if

runtime => sensor % runtime

if (.not. associated(runtime)) then
   write(string1,*)'Runtime information not initialized for platform/sat/sensor id combination:',&
      instrument
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
end if

obsloc   = get_location(location)

lon      = obsloc(1) ! degree: 0 to 360
lat      = obsloc(2) ! degree: -90 to 90
! ignore the 3rd dimension of the observation location for now
!height   = obsloc(3) ! (m)

! --------------------------------------------------------------------------
! 4. Build the list of profile/channel indices in chanprof
! --------------------------------------------------------------------------

nch = 0_jpim
do j = 1, ens_size
    nch = nch + 1_jpim
    runtime % chanprof(nch) % prof = j
    runtime % chanprof(nch) % chan = channel
end do

! We would like a level index array to allow either surface first or surface last order

! One would assume the number of levels would not change between calls, but check
if (allocated(lvlidx) .and. size(lvlidx) /= nlevels) then
   write(string1,*)'Number of levels changed from ',size(lvlidx),' to ',nlevels,&
      ' for platform/sat/sensor id combination:',instrument
   call error_handler(E_WARN, routine, string1, source, revision, revdate, text2=string2)
   
   ! deallocate so we can try again here with the right number of levels
   deallocate(lvlidx)
   deallocate(ly1idx)
   deallocate(ly2idx)
   deallocate(totalwater)
   deallocate(totalice)
end if

! (re)allocate if necessary
if (.not. allocated(lvlidx)) then
   allocate(lvlidx(nlevels))
   allocate(ly1idx(nlevels-1))
   allocate(ly2idx(nlevels-1))
   allocate(totalwater(nlevels))
   allocate(totalice(nlevels))
end if

! finally set the array to the correct order
if (first_lvl_is_sfc) then
   do ilvl=1,nlevels
      lvlidx(ilvl) = nlevels - ilvl + 1
   end do
else
   do ilvl=nlevels,1,-1
      lvlidx(ilvl) = ilvl
   end do
end if

! used for averaging levels to layers
ly1idx = lvlidx(1:nlevels-1)
ly2idx = lvlidx(2:nlevels)

! Loop over all of the ensemble members
! There is one profile per ensemble member
DO imem = 1, ens_size

   surftype   = nint(atmos%surftype(imem))

   runtime % profiles(imem) % nlevels = nlevels
   runtime % profiles(imem) % nlayers = nlevels - 1
   runtime % profiles(imem) % gas_units = 1 ! 1 = kg/kg, 2 = ppmv
   runtime % profiles(imem) % mmr_cldaer = .true. ! kg/kg
   runtime % profiles(imem) % clw_scheme = clw_scheme
   runtime % profiles(imem) % ice_scheme = ice_scheme
   runtime % profiles(imem) % idg        = idg_scheme

   ! set the required profile variables 
   runtime % profiles(imem) % p(:) = atmos % pressure(imem,lvlidx)/100.d0  ! Pa -> hPa
   runtime % profiles(imem) % t(:) = atmos % temperature(imem,lvlidx) 
   runtime % profiles(imem) % q(:) = max(atmos % moisture(imem,lvlidx),1d-8) 

   ! set trace gases if opts present and individual flags are used
   ! the arrays are assumed to have been allocated - this was already checked in obs_def
   if (is_visir .or. mw_clear_sky_only) then
      if (runtime % opts % rt_ir % ozone_data) then
         runtime % profiles(imem) % o3(:) = trace_gas % ozone(imem,lvlidx)
      end if

      if (runtime % opts % rt_ir % co2_data) then
         runtime % profiles(imem) % co2(:) = trace_gas % co2(imem,lvlidx)
      end if

      if (runtime % opts % rt_ir % n2o_data) then
         runtime % profiles(imem) % n2o(:) = trace_gas % n2o(imem,lvlidx)
      end if

      if (runtime % opts % rt_ir % ch4_data) then
         runtime % profiles(imem) % ch4(:) = trace_gas % ch4(imem,lvlidx)
      end if

      if (runtime % opts % rt_ir % co_data) then
         runtime % profiles(imem) % co(:)  = trace_gas % co(imem,lvlidx)
      end if

      if (runtime % opts % rt_ir % so2_data) then
         runtime % profiles(imem) % so2(:) = trace_gas % so2(imem,lvlidx)
      end if

      ! "clear-sky" cloud-liquid water data for microwave (not RTTOV-SCATT)
      if (is_mw .and. runtime % opts % rt_mw % clw_data) then
         runtime % profiles(imem) % clw(:) = clouds % clw(imem,lvlidx)
      end if

      ! add aerosols
      if (runtime % opts % rt_ir % addaerosl) then
         if (aerosl_type == 1) then 
            ! OPAC, convert from levels to layers
            runtime % profiles(imem) % aerosols(1,:) = 0.5d0*(aerosols % insoluble(imem,ly1idx) +           &
                                                              aerosols % insoluble(imem,ly2idx))
            runtime % profiles(imem) % aerosols(2,:) = 0.5d0*(aerosols % water_soluble(imem,ly1idx) +      &
                                                              aerosols % water_soluble(imem,ly2idx))
            runtime % profiles(imem) % aerosols(3,:) = 0.5d0*(aerosols % soot(imem,ly1idx) +                &
                                                              aerosols % soot(imem,ly2idx))
            runtime % profiles(imem) % aerosols(4,:) = 0.5d0*(aerosols % sea_salt_accum(imem,ly1idx) +      &
                                                              aerosols % sea_salt_accum(imem,ly2idx))
            runtime % profiles(imem) % aerosols(5,:) = 0.5d0*(aerosols % sea_salt_coarse(imem,ly1idx) +     &
                                                              aerosols % sea_salt_coarse(imem,ly2idx))
            runtime % profiles(imem) % aerosols(6,:) = 0.5d0*(aerosols % mineral_nucleus(imem,ly1idx) +     &
                                                              aerosols % mineral_nucleus(imem,ly2idx))
            runtime % profiles(imem) % aerosols(7,:) = 0.5d0*(aerosols % mineral_accum(imem,ly1idx) +       &
                                                              aerosols % mineral_accum(imem,ly2idx))
            runtime % profiles(imem) % aerosols(8,:) = 0.5d0*(aerosols % mineral_coarse(imem,ly1idx) +      &
                                                              aerosols % mineral_coarse(imem,ly2idx))
            runtime % profiles(imem) % aerosols(9,:) = 0.5d0*(aerosols % mineral_transport(imem,ly1idx) +   &
                                                              aerosols % mineral_transport(imem,ly2idx))
            runtime % profiles(imem) % aerosols(10,:) = 0.5d0*(aerosols % sulphated_droplets(imem,ly1idx) + &
                                                               aerosols % sulphated_droplets(imem,ly2idx))
            runtime % profiles(imem) % aerosols(11,:) = 0.5d0*(aerosols % volcanic_ash(imem,ly1idx) +       &
                                                               aerosols % volcanic_ash(imem,ly2idx))
            runtime % profiles(imem) % aerosols(12,:) = 0.5d0*(aerosols % new_volcanic_ash(imem,ly1idx) +   &
                                                               aerosols % new_volcanic_ash(imem,ly2idx))
            runtime % profiles(imem) % aerosols(13,:) = 0.5d0*(aerosols % asian_dust(imem,ly1idx) +         &
                                                               aerosols % asian_dust(imem,ly2idx))
         elseif (aerosl_type == 2) then
            ! CAMS, convert from levels to levels
            runtime % profiles(imem) % aerosols(1,:) = 0.5d0*(aerosols % black_carbon(imem,ly1idx) +               &
                                                              aerosols % black_carbon(imem,ly2idx))
            runtime % profiles(imem) % aerosols(2,:) = 0.5d0*(aerosols % dust_bin1(imem,ly1idx) +                  &
                                                              aerosols % dust_bin1(imem,ly2idx))
            runtime % profiles(imem) % aerosols(3,:) = 0.5d0*(aerosols % dust_bin2(imem,ly1idx) +                  &
                                                              aerosols % dust_bin2(imem,ly2idx))
            runtime % profiles(imem) % aerosols(4,:) = 0.5d0*(aerosols % dust_bin3(imem,ly1idx) +                  &
                                                              aerosols % dust_bin3(imem,ly2idx))
            runtime % profiles(imem) % aerosols(5,:) = 0.5d0*(aerosols % ammonium_sulphate(imem,ly1idx) +          &
                                                              aerosols % ammonium_sulphate(imem,ly2idx))
            runtime % profiles(imem) % aerosols(6,:) = 0.5d0*(aerosols % sea_salt_bin1(imem,ly1idx) +              &
                                                              aerosols % sea_salt_bin1(imem,ly2idx))
            runtime % profiles(imem) % aerosols(7,:) = 0.5d0*(aerosols % sea_salt_bin2(imem,ly1idx) +              &
                                                              aerosols % sea_salt_bin2(imem,ly2idx))
            runtime % profiles(imem) % aerosols(8,:) = 0.5d0*(aerosols % sea_salt_bin3(imem,ly1idx) +              &
                                                              aerosols % sea_salt_bin3(imem,ly2idx))
            runtime % profiles(imem) % aerosols(9,:) = 0.5d0*(aerosols % hydrophilic_organic_matter(imem,ly1idx) + &
                                                              aerosols % hydrophilic_organic_matter(imem,ly2idx))
         else
            write(string1,*)'Unknown aerosol type ',aerosl_type,&
               ' for platform/sat/sensor id combination:',instrument
            call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
         end if
      end if ! add aerosols

      ! add IR clouds
      if (runtime % opts % rt_ir % addclouds) then
         ! first find the total water and ice components

         totalwater(:) = 0.d0

         if (allocated(clouds % clw)) then
            totalwater(:) = totalwater(:) + clouds % clw(imem,lvlidx)
         end if

         if (allocated(clouds % rain)) then
            totalwater(:) = totalwater(:) + clouds % clw(imem,lvlidx)
         end if

         totalice(:) = 0.d0

         if (allocated(clouds % ciw)) then
            totalice(:) = totalice(:) + clouds % ciw(imem,lvlidx)
         end if 

         if (allocated(clouds % snow)) then
            totalice(:) = totalice(:) + clouds % snow(imem,lvlidx)
         end if 

         if (allocated(clouds % graupel)) then
            totalice(:) = totalice(:) + clouds % graupel(imem,lvlidx)
         end if 

         if (allocated(clouds % hail)) then
            totalice(:) = totalice(:) + clouds % hail(imem,lvlidx)
         end if 

         ! Classify liquid water cloud type. If the maximum absolute w is > 0.5 m/s, classify as a cumulus cloud.
         ! FIXME: consider adding 0.5d0 as a namelist parameter
         if (allocated(clouds % w)) then
            maxw = maxval(abs(clouds % w(imem,:)))
            is_cumulus = maxw > 0.5d0 ! m/s
         else
            ! assume cumulus if w not provided
            is_cumulus = .true. 
         end if

         ! depending on the vertical velocity and land type, classify clouds the way RTTOV wants 
         if (.not. is_cumulus) then
            ! stratus
            if (surftype == 0) then
               ! 1, Stratus Continental, STCO (kg/kg), convert levels to layers
               runtime % profiles(imem) % cloud(1,:) = 0.5d0*(totalwater(ly1idx) + &
                                                              totalwater(ly2idx))
            else
               ! 2, Stratus Maritime, STMA (kg/kg), convert levels to layers
               runtime % profiles(imem) % cloud(2,:) = 0.5d0*(totalwater(ly1idx) + &
                                                              totalwater(ly2idx))
            end if
         else
            ! cumulus
            if (surftype == 0) then
               ! 3, Cumulus Continental Clean, CUCC (kg/kg), convert levels to layers
               runtime % profiles(imem) % cloud(3,:) = 0.5d0*(totalwater(ly1idx) + &
                                                              totalwater(ly2idx))
               ! FIXME: give user (or obs) control over CUCC versus CUCP
               ! 4, Cumulus Continental Polluted, CUCP (kg/kg), convert levels to layers
               ! runtime % profiles(imem) % cloud(4,:) = 0.5d0*(totalwater(ly1idx) + totalwater(ly2idx))
            else
               ! 5, Cumulus Maritime, CUMA (kg/kg), convert levels to layers
               runtime % profiles(imem) % cloud(5,:) = 0.5d0*(totalwater(ly1idx) + &
                                                              totalwater(ly2idx))
            end if
         end if

         ! allow specification of cloud water effective diameter
         if (allocated(clouds % clwde)) then
            ! Liquid water effective diameter (microns), convert levels to layers
            runtime % profiles(imem)% clwde(:) = 0.5d0*(clouds % clwde(imem,ly1idx) + &
                                                        clouds % clwde(imem,ly2idx))
         end if

         ! all types of ice cloud (kg/kg) go in 6, Ice Cloud (CIRR, but not just cirrus)
         ! convert levels to layers
         runtime % profiles(imem) % cloud(6,:) = 0.5d0*(totalice(ly1idx) + &
                                                        totalice(ly2idx))

         ! allow specification of ice effective diameter
         if (allocated(clouds % icede)) then
            ! Ice effective diameter (microns), convert levels to layers
            runtime % profiles(imem) % icede(:) = 0.5d0*(clouds % icede(imem,ly1idx) + &
                                                         clouds % icede(imem,ly2idx))
         end if

         if (allocated(clouds % cfrac)) then
            ! Cloud fraction (0-1), convert levels to layers
            runtime % profiles(imem) % cfrac(:) = 0.5d0*(clouds % cfrac(imem,ly1idx) + &
                                                         clouds % cfrac(imem,ly2idx))
         end if
      end if ! add IR clouds
   else if (is_mw) then
      ! RTTOV-SCATT, add MW clouds
      if (allocated(clouds % cfrac)) then
         runtime % cld_profiles(imem) % cc(:) = clouds % cfrac(imem, lvlidx)
      else
         runtime % cld_profiles(imem) % cc(:) = 0.d0
      end if


      if (allocated(clouds % clw)) then
         runtime % cld_profiles(imem) % clw(:) = clouds % clw(imem,lvlidx)
      else
         runtime % cld_profiles(imem) % clw = 0.d0
      end if

      if (allocated(clouds % rain)) then
         runtime % cld_profiles(imem) % rain(:) = clouds % rain(imem,lvlidx)
      else
         runtime % cld_profiles(imem) % rain = 0.d0
      end if

      if (use_totalice) then
         totalice(:) = 0.d0

         if (allocated(clouds % ciw)) then
            totalice(:) = totalice(:) + clouds % ciw(imem,lvlidx)
         end if 

         if (allocated(clouds % snow)) then
            totalice(:) = totalice(:) + clouds % snow(imem,lvlidx)
         end if 

         if (allocated(clouds % graupel)) then
            totalice(:) = totalice(:) + clouds % graupel(imem,lvlidx)
         end if 

         if (allocated(clouds % hail)) then
            totalice(:) = totalice(:) + clouds % hail(imem,lvlidx)
         end if 

         runtime % cld_profiles(imem) % totalice(:) = totalice(:)
      else
         if (allocated(clouds % ciw)) then
            runtime % cld_profiles(imem) % ciw(:) = clouds % ciw(imem,lvlidx)
         end if 

         ! "totalice" here is being used only for solid precipitation (no ciw)
         totalice(:) = 0.d0

         if (allocated(clouds % snow)) then
            totalice(:) = totalice(:) + clouds % snow(imem,lvlidx)
         end if 

         if (allocated(clouds % graupel)) then
            totalice(:) = totalice(:) + clouds % graupel(imem,lvlidx)
         end if 

         if (allocated(clouds % hail)) then
            totalice(:) = totalice(:) + clouds % hail(imem,lvlidx)
         end if 

         runtime % cld_profiles(imem) % sp(:) = totalice(:)
      end if

      ! also add "half-level pressures" as requested by RTTOV-Scatt
      runtime % cld_profiles(imem) % ph(2:nlevels) = 0.5d0*(atmos % pressure(imem,ly1idx)+atmos % pressure(imem,ly2idx))/100.d0
      runtime % cld_profiles(imem) % ph(nlevels+1) = atmos % sfc_p(imem)/100.d0
      runtime % cld_profiles(imem) % ph(1) = 0.d0
   else
      write(string1,*)'Neither opts or opts_scatter were available for platform/sat/sensor id combination:',&
         instrument
      call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
   end if ! has runtime % opts or runtime % opts_scatt

   ! 2m above surface pressure (hPa)
   runtime % profiles(imem) % s2m % p  = atmos % sfc_p(imem)/100.d0 ! convert from Pa to hPa
   ! 2m above surface temperature (K)
   runtime % profiles(imem) % s2m % t  = atmos % s2m_t(imem)

   ! set the optional variables
   if (allocated(atmos % s2m_q)) then
      ! 2m above surface water vapor (kg/kg)
      runtime % profiles(imem) % s2m % q  = atmos % s2m_q(imem)
   end if
 
   if (allocated(atmos % s10m_u) .and. allocated(atmos % s10m_v)) then
      ! 10 m above surface u and v wind (m/s)
      runtime % profiles(imem) % s2m % u  = atmos % s10m_u(imem)
      runtime % profiles(imem) % s2m % v  = atmos % s10m_v(imem)
   end if
 
   if (allocated(atmos % wfetch)) then
      ! Wind fetch over the ocean (m)
      runtime % profiles(imem) % s2m % wfetc = atmos % wfetch(imem)  
   end if
   
   ! Surface type (0=land, 1=sea, 2=sea-ice)
   runtime % profiles(imem) % skin % surftype  = surftype

   if (allocated(atmos % water_type)) then 
      ! Water type (0=fresh, 1=ocean)
      runtime % profiles(imem) % skin % watertype = nint(atmos % water_type(imem))
   end if

   ! Surface skin temperature (K) 
   runtime % profiles(imem) % skin % t = atmos % skin_temp(imem)

   if (allocated(atmos % sfc_salinity)) then
      ! Surface ocean salinity (practical salinity unit)  
      runtime % profiles(imem) % skin % salinity = atmos % sfc_salinity(imem)
   end if

   if (allocated(atmos % sfc_foam_frac)) then
      ! Surface foam fraction (0-1)
      runtime % profiles(imem) % skin % foam_fraction = atmos % sfc_foam_frac(imem)
   end if

   if (allocated(atmos % sfc_snow_frac)) then
      ! Surface snow fraction (0-1)
      runtime % profiles(imem) % skin % snow_fraction = atmos % sfc_snow_frac(imem)
   end if

   if (is_mw .and. use_fastem_params) then
      ! FASTEM parameters, see RTTOV user guide e.g. Table 21
      runtime % profiles(imem) % skin % fastem(1)   = mw_md%fastem_land1
      runtime % profiles(imem) % skin % fastem(2)   = mw_md%fastem_land2
      runtime % profiles(imem) % skin % fastem(3)   = mw_md%fastem_land3
      runtime % profiles(imem) % skin % fastem(4)   = mw_md%fastem_land4
      runtime % profiles(imem) % skin % fastem(5)   = mw_md%fastem_land5
   end if

   if (is_visir .and. do_lambertian) then
      ! Surface specularity (0-1)
      runtime % profiles(imem) % skin % specularity = visir_md%specularity
   end if

   if (allocated(clouds % simple_cfrac)) then
      ! Simple (column) cloud fraction, 0-1
      runtime % profiles(imem) % cfraction = clouds % simple_cfrac(imem)
   end if

   if (allocated(clouds % ctp)) then
      ! Simple (column) cloud top pressure, hPa
      runtime % profiles(imem) % ctp = clouds % ctp(imem)
   end if

   if (is_visir) then
      ! Sat. zenith and azimuth angles (degrees)
      runtime % profiles(imem) % zenangle    = visir_md % sat_ze
      runtime % profiles(imem) % azangle     = visir_md % sat_az

      ! Solar zenith and azimuth angles (degrees), only relevant if use_solar
      runtime % profiles(imem) % sunzenangle = visir_md % sun_az
      runtime % profiles(imem) % sunazangle  = visir_md % sun_ze
   else
      ! Sat. zenith and azimuth angles (degrees)
      runtime % profiles(imem) % zenangle    = mw_md % sat_ze
      runtime % profiles(imem) % azangle     = mw_md % sat_az
   end if

   ! Elevation (km), latitude and longitude (degrees)
   runtime % profiles(imem) % latitude  = lat
   runtime % profiles(imem) % longitude = lon
   runtime % profiles(imem) % elevation = atmos % sfc_elev(imem)/1000.d0 ! m -> km

   if (use_zeeman .and. is_mw) then
      runtime % profiles(imem) % Be    = mw_md % mag_field ! mag field strength, Gauss
      runtime % profiles(imem) % cosbk = mw_md % cosbk     ! cosine of angle between mag field and viewing angle
   end if
end do ! profile data

! --------------------------------------------------------------------------
! 6. Specify surface emissivity and reflectance
! --------------------------------------------------------------------------

! In this example we have no values for input emissivities
runtime % emissivity(:) % emis_in = 0._jprb

! Calculate emissivity within RTTOV where the input emissivity value is
! zero or less (all channels in this case)
runtime % calcemis(:) = (runtime % emissivity(:) % emis_in <= 0._jprb)

! In this example we have no values for input reflectances
runtime % reflectance(:) % refl_in = 0._jprb

! Calculate BRDF within RTTOV where the input BRDF value is zero or less
! (all channels in this case)
runtime % calcrefl(:) = (runtime % reflectance(:) % refl_in <= 0._jprb)

! Use default cloud top BRDF for simple cloud in VIS/NIR channels
runtime % reflectance(:) % refl_cloud_top = 0._jprb

if (debug) then
   call rttov_print_opts(runtime % opts, lu=ioout)
   call rttov_print_profile(runtime % profiles(1), lu=ioout)
end if

if (is_visir .or. mw_clear_sky_only) then
   ! Call RTTOV forward model
   call rttov_direct (                         &
           errorstatus,                        & ! out   error flag
           runtime % chanprof,                 & ! in    channel and profile index structure
           runtime % opts,                     & ! in    options structure
           runtime % profiles,                 & ! in    profile array
           runtime % coefs,                    & ! in    coefficients structure
           runtime % transmission,             & ! inout computed transmittances
           runtime % radiance,                 & ! inout computed radiances
           calcemis    = runtime % calcemis,   & ! in    flag for internal emissivity calcs
           emissivity  = runtime % emissivity, & ! inout input/output emissivities per channel
           calcrefl    = runtime % calcrefl,   & ! in    flag for internal BRDF calcs
           reflectance = runtime % reflectance ) ! inout input/output BRDFs per channel

   if (errorstatus /= errorstatus_success) then
      write(string1,*)'RTTOV direct error for platform/satellite/sensor id:',&
         instrument
      call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
   end if
else
   ! set the frequencies based on which channel is being called
   runtime % frequencies(:) = runtime % frequencies_all(channel)

   if (debug) then
      call rttov_print_opts_scatt(runtime % opts_scatt, lu=ioout)
      call rttov_print_cld_profile(runtime % cld_profiles(1), lu=ioout)
   end if

   call rttov_scatt (          &
      errorstatus,            & ! out   error flag
      runtime % opts_scatt,   & ! in    RTTOV-SCATT options structure
      nlevels,                & ! in    number of profile levels
      runtime % chanprof,     & ! in    channel and profile index structure
      runtime % frequencies,  & ! in    channel indexes for Mietable lookup
      runtime % profiles,     & ! in    profile array
      runtime % cld_profiles, & ! in    cloud/hydrometeor profile array
      runtime % coefs,        & ! in    coefficients structure
      runtime % coefs_scatt,  & ! in    Mietable structure
      runtime % calcemis,     & ! in    flag for internal emissivity calcs
      runtime % emissivity,   & ! inout input/output emissivities per channel
      runtime % radiance )      ! inout computed radiances
   
   if (errorstatus /= errorstatus_success) then
      write(string1,*)'RTTOV scatt error for platform/satellite/sensor id:',&
         instrument
      call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=string2)
   end if
end if

if (is_visir) then
   do imem = 1, ens_size
      radiances(imem) = runtime % radiance % total(imem)
   end do
   if (debug) then
      print*, 'RADIANCE % TOTAL = ', runtime % radiance % total(:)
   end if
else
   do imem = 1, ens_size
      radiances(imem) = runtime % radiance % bt(imem)
   end do
   if (debug) then
      print*, 'RADIANCE % BT = ', runtime % radiance % bt(:)
   end if
end if

IF (errorstatus /= errorstatus_success) THEN
  WRITE (*,*) 'rttov_direct error'
  !CALL rttov_exit(errorstatus)
ENDIF

end subroutine do_forward_model

!--------------------------------------------------------------------------

subroutine sensor_runtime_takedown(runtime, error_status)

type(rttov_sensor_runtime_type), pointer :: runtime
integer, intent(out) :: error_status

integer :: nchanprof
integer :: nlevels

! Return error status of RTTOV subroutine calls
integer(kind=jpim)               :: errorstatus 

!FIXME - in the forward operator code we won't be deallocating
! this structure.

! initialize error status tp success
error_status = errorstatus_success

! Deallocate all RTTOV arrays and structures

nchanprof = size(runtime % chanprof)
nlevels = size(runtime % profiles(1) % p)

call rttov_alloc_direct(               &
      errorstatus,                     &
      0_jpim,                          &  ! 0 => deallocate
      nchanprof,                       &  ! DART only does one channel at a time, so nchanprof = nprof
      nchanprof,                       &
      nlevels,                         &
      runtime % chanprof,              &
      runtime % opts,                  &
      runtime % profiles,              &
      runtime % coefs,                 &
      runtime % transmission,          &
      runtime % radiance,              &
      calcemis=runtime % calcemis,     &
      emissivity=runtime % emissivity, &
      calcrefl=runtime % calcrefl,     &
      reflectance=runtime % reflectance)

IF (errorstatus /= errorstatus_success) THEN
  WRITE(*,*) 'deallocation error for rttov_direct structures'
  error_status = errorstatus
  return
ENDIF

CALL rttov_dealloc_coefs(errorstatus, runtime % coefs)
IF (errorstatus /= errorstatus_success) THEN
  WRITE(*,*) 'coefs deallocation error'
ENDIF

end subroutine sensor_runtime_takedown

!--------------------------------------------------------------------------

end module rttov_interface_mod
