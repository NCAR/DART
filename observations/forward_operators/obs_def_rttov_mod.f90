! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!----------------------------------------------------------------------
! This module provides support for computing forward operators for
! radiance observations by calling the rttov radiative transfer model.
!
! the additional metadata in each observation includes:
!
!  OBS            X
! rttov
!    sat az/el
!    sun az/el
!    platform
!    instrument
!    channel
!    <anything else useful>
!----------------------------------------------------------------------

! BEGIN DART PREPROCESS KIND LIST
! AIRS_AMSU_RADIANCE,    QTY_RADIANCE
! END DART PREPROCESS KIND LIST


! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_rttov_mod, only : read_rttov_metadata, &
!                                write_rttov_metadata, &
!                          interactive_rttov_metadata, &
!                                get_expected_radiance
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(AIRS_AMSU_RADIANCE)
!         call get_expected_radiance(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF


! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(AIRS_AMSU_RADIANCE)
!      call read_rttov_metadata(obs_def%key, key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(AIRS_AMSU_RADIANCE)
!      call write_rttov_metadata(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(AIRS_AMSU_RADIANCE)
!      call interactive_rttov_metadata(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_rttov_mod

use        types_mod, only : r8, PI, metadatalength, MISSING_R8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_WARN, E_MSG, &
                             logfileunit, get_unit, open_file, close_file, nc_check, &
                             file_exist, ascii_file_format
use     location_mod, only : location_type, set_location, get_location, &
                             VERTISHEIGHT, VERTISLEVEL, set_location_missing
use     obs_kind_mod, only : QTY_GEOPOTENTIAL_HEIGHT, QTY_SOIL_MOISTURE
use  assim_model_mod, only : interpolate

use obs_def_utilities_mod, only : track_status
use ensemble_manager_mod,  only : ensemble_type

use typesizes
use netcdf

implicit none
private

public ::            set_rttov_metadata, &
                     get_rttov_metadata, &
                    read_rttov_metadata, &
                   write_rttov_metadata, &
             interactive_rttov_metadata, &
                  get_expected_radiance

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=256) :: string1, string2
logical, save      :: module_initialized = .false.

! Metadata for rttov observations.
! There are soil parameters for each site that must be added to each
! observation in the sequence. Also COSMIC parameters ...

!FIXME
type site_metadata
   private
   real(r8)            :: bd         ! Dry Soil Bulk Density [  g / cm^3]
   real(r8)            :: lattwat    ! Lattice Water Content [M^3 /  M^3]
   real(r8)            :: N          ! High Energy Neutron Intensity
   real(r8)            :: alpha      ! Ratio of Fast Neutron Creation Factor
   real(r8)            :: L1         ! High Energy   Soil Attenuation Length
   real(r8)            :: L2         ! High Energy  Water Attenuation Length
   real(r8)            :: L3         ! Fast Neutron  Soil Attenuation Length
   real(r8)            :: L4         ! Fast Neutron Water Attenuation Length
end type site_metadata

type(site_metadata), allocatable, dimension(:) :: observation_metadata
type(site_metadata) :: missing_metadata
character(len=6), parameter :: RTTOVSTRING = 'rttov'

logical :: debug = .FALSE.
integer :: MAXrttovkey = 24*366  !FIXME
integer ::    rttovkey = 0       ! useful length of metadata arrays

contains


!----------------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)

module_initialized = .true.

missing_metadata%bd       = MISSING_R8
missing_metadata%lattwat  = MISSING_R8
missing_metadata%N        = MISSING_R8
missing_metadata%alpha    = MISSING_R8
missing_metadata%L1       = MISSING_R8
missing_metadata%L2       = MISSING_R8
missing_metadata%L3       = MISSING_R8
missing_metadata%L4       = MISSING_R8

allocate(observation_metadata(MAXrttovkey))

observation_metadata(:) = missing_metadata

end subroutine initialize_module



 subroutine set_rttov_metadata(key, bd, lattwat, N, alpha, L1, L2, L3, L4)
!----------------------------------------------------------------------
!subroutine set_rttov_metadata(key, bd, lattwat, N, alpha, L1, L2, L3, L4)
!
! Fill the module storage metadata for a particular observation.

integer,  intent(out) :: key
real(r8), intent(in)  :: bd, lattwat, N, alpha, L1, L2, L3, L4

if ( .not. module_initialized ) call initialize_module

rttovkey = rttovkey + 1  ! increase module storage used counter

! Make sure the new key is within the length of the metadata arrays.
call grow_metadata(rttovkey,'set_rttov_metadata')

key = rttovkey ! now that we know its legal

observation_metadata(key)%bd      = bd
observation_metadata(key)%lattwat = lattwat
observation_metadata(key)%N       = N
observation_metadata(key)%alpha   = alpha
observation_metadata(key)%L1      = L1
observation_metadata(key)%L2      = L2
observation_metadata(key)%L3      = L3
observation_metadata(key)%L4      = L4

end subroutine set_rttov_metadata


!----------------------------------------------------------------------
! Query the metadata in module storage for a particular observation.
! This can be useful for post-processing routines, etc.

subroutine get_rttov_metadata(key, bd, lattwat, N, alpha, L1, L2, L3, L4)
integer,  intent(in)  :: key
real(r8), intent(out) :: bd, lattwat, N, alpha, L1, L2, L3, L4

if ( .not. module_initialized ) call initialize_module

! Make sure the desired key is within the length of the metadata arrays.
call key_within_range(key,'get_rttov_metadata')

bd       = observation_metadata(key)%bd
lattwat  = observation_metadata(key)%lattwat
N        = observation_metadata(key)%N
alpha    = observation_metadata(key)%alpha
L1       = observation_metadata(key)%L1
L2       = observation_metadata(key)%L2
L3       = observation_metadata(key)%L3
L4       = observation_metadata(key)%L4

end subroutine get_rttov_metadata


!----------------------------------------------------------------------
! This routine reads the metadata for radiance obs

subroutine read_rttov_metadata(key, obsID, ifile, fform)
integer,          intent(out)          :: key    ! index into local metadata
integer,          intent(in)           :: obsID
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

! temp variables
logical           :: is_asciifile
integer           :: ierr
character(len=6)  :: header
integer           :: oldkey
real(r8)          :: bd, lattwat, N, alpha, L1, L2, L3, L4

if ( .not. module_initialized ) call initialize_module

is_asciifile = ascii_file_format(fform)

write(string2,*)'observation #',obsID

if ( is_asciifile ) then
   read(ifile, *, iostat=ierr) header
   call check_iostat(ierr,'read_rttov_metadata','header',string2)
   if (trim(header) /= trim(rttovSTRING)) then
       write(string1,*)"Expected radiance header ["//RTTOVSTRING//"] in input file, got ["//header//"]"
       call error_handler(E_ERR, 'read_rttov_metadata', string1, source, revision, revdate, text2=string2)
   endif
   read(ifile, *, iostat=ierr) bd, lattwat, N, alpha
   call check_iostat(ierr,'read_rttov_metadata','bd -> alpha',string2)
   read(ifile, *, iostat=ierr) L1, L2, L3, L4
   call check_iostat(ierr,'read_rttov_metadata','L1 -> L4',string2)
   read(ifile, *, iostat=ierr) oldkey
   call check_iostat(ierr,'read_rttov_metadata','oldkey',string2)
else
   read(ifile, iostat=ierr) header
   call  check_iostat(ierr,'read_rttov_metadata','header',string2)
   if (trim(header) /= trim(rttovSTRING)) then
       write(string1,*)"Expected radiance header ["//RTTOVSTRING//"] in input file, got ["//header//"]"
       call error_handler(E_ERR, 'read_rttov_metadata', string1, source, revision, revdate, text2=string2)
   endif
   read(ifile, iostat=ierr) bd, lattwat, N, alpha
   call  check_iostat(ierr,'read_rttov_metadata','bd -> alpha',string2)
   read(ifile, iostat=ierr) L1, L2, L3, L4
   call  check_iostat(ierr,'read_rttov_metadata','L1 -> L4',string2)
   read(ifile, iostat=ierr) oldkey
   call  check_iostat(ierr,'read_rttov_metadata','oldkey',string2)
endif

! The oldkey is thrown away.

! Store the metadata in module storage and record the new length of the metadata arrays.
call set_rttov_metadata(key, bd, lattwat, N, alpha, L1, L2, L3, L4)

! The new 'key' is returned.

end subroutine read_rttov_metadata


!----------------------------------------------------------------------
! writes the metadata for radiance observations.

subroutine write_rttov_metadata(key, ifile, fform)

integer,           intent(in)           :: key
integer,           intent(in)           :: ifile
character(len=*),  intent(in), optional :: fform

logical  :: is_asciifile
real(r8) :: bd, lattwat, N, alpha, L1, L2, L3, L4

if ( .not. module_initialized ) call initialize_module

! given the index into the local metadata arrays - retrieve
! the metadata for this particular observation.

call get_rttov_metadata(key, bd, lattwat, N, alpha, L1, L2, L3, L4)

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
   write(ifile, *) trim(rttovSTRING)
   write(ifile, *) bd, lattwat, N, alpha
   write(ifile, *) L1, L2, L3, L4
   write(ifile, *) key
else
   write(ifile   ) trim(rttovSTRING)
   write(ifile   ) bd, lattwat, N, alpha
   write(ifile   ) L1, L2, L3, L4
   write(ifile   ) key
endif

end subroutine write_rttov_metadata


!----------------------------------------------------------------------

subroutine interactive_rttov_metadata(key)
integer, intent(out) :: key

real(r8) :: bd, lattwat, N, alpha, L1, L2, L3, L4

if ( .not. module_initialized ) call initialize_module

! Prompt for input for the required metadata

bd      = interactive('"bd"      dry soil bulk density [g/cm^3]'                ,minvalue = 0.0_r8)
lattwat = interactive('"lattwat" lattice water content [m^3/m^3]'               ,minvalue = 0.0_r8)
N       = interactive('"N"       high energy neutron intensity [count]'         ,minvalue = 0.0_r8)
alpha   = interactive('"alpha"   ratio of fast neutron creation factor [Soil to Water]',minvalue=0.0_r8)
L1      = interactive('"L1"      high energy soil attenuation length [g/cm^2]'  ,minvalue = 0.0_r8)
L2      = interactive('"L2"      high energy water attenuation length [g/cm^2]' ,minvalue = 0.0_r8)
L3      = interactive('"L3"      fast neutron soil attenuation length [g/cm^2]' ,minvalue = 0.0_r8)
L4      = interactive('"L4"      fast neutron water attenuation length [g/cm^2]',minvalue = 0.0_r8)

call set_rttov_metadata(key, bd, lattwat, N, alpha, L1, L2, L3, L4)

end subroutine interactive_rttov_metadata


!----------------------------------------------------------------------

function interactive(str1,minvalue,maxvalue)
real(r8)                       :: interactive
character(len=*),   intent(in) :: str1
real(r8), optional, intent(in) :: minvalue
real(r8), optional, intent(in) :: maxvalue

integer :: i

interactive = MISSING_R8

! Prompt with a minimum amount of error checking

if     (present(minvalue) .and. present(maxvalue)) then

   interactive = minvalue - 1.0_r8
   MINMAXLOOP : do i = 1,10
      if ((interactive >= minvalue) .and. (interactive <= maxvalue)) exit MINMAXLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive
   end do MINMAXLOOP

elseif (present(minvalue)) then

   interactive = minvalue - 1.0_r8
   MINLOOP : do i=1,10
      if (interactive >= minvalue) exit MINLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive
   end do MINLOOP

elseif (present(maxvalue)) then

   interactive = maxvalue + 1.0_r8
   MAXLOOP : do i=1,10
      if (interactive <= maxvalue) exit MAXLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive
   end do MAXLOOP

else ! anything goes ... cannot check
      write(*, *) 'Enter '//str1
      read( *, *) interactive
endif

end function interactive


!----------------------------------------------------------------------

 subroutine get_expected_radiance(state_handle, ens_size, location, key, val, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location          ! location of obs
integer,             intent(in)  :: key               ! key into module metadata
real(r8),            intent(out) :: val(ens_size)     ! value of obs
integer,             intent(out) :: istatus(ens_size) ! status of the calculation

!FIXME - this all gets replaced by code from the example program

!========================================================================================
! COsmic-ray Soil Moisture Interaction Code (COSMIC) - Version 1.5
!
! W. James Shuttleworth and Rafael Rosolem - January/2012
! Fortran code developed by Rafael Rosolem
!========================================================================================
! Updates:
! 01/20/2012 - Version 1.0: * Original version based on SPAM
! 01/28/2012 - Version 1.1: * Some parameters are re-defined for better physical realism
! 02/17/2012 - Version 1.2: * After contribution from all angles are taken, need to
!                             multiply fastflux by 2/PI
!                           * Angle increments can be specified here (ideg)
!                             resolution
! 02/29/2012 - Version 1.3: * Reduced number of parameters based on relationship of ns
!                             and nw (now given as alpha = ns/nw)
! 04/03/2012 - Version 1.4: * Soil thickness (i.e., input.dat file) needs to be specified
!                             at finer resolution (i.e., 0.1 cm)
! 04/04/2012 - Version 1.5  * Now the contributions to soil and water densities/mass are
!                             taken to be at the center of a given soil layer
!=================================================================================
! COSMIC: Variables list
!=================================================================================

!rr: Use 1 mm layer increments (down to 3 meters) to compute the weighting function
!rr: (as originally done for COSMIC), and then compute the cumulative weights for
!rr: individual soil layers

integer, parameter :: nlyr=3000 ! Total number of soil layers - each 1mm for 3m

real(r8) :: bd     = 0.0_r8 ! Dry soil bulk density (g/m3)
real(r8) :: vwclat = 0.0_r8 ! Volumetric "lattice" water content (m3/m3)
real(r8) :: N      = 0.0_r8 ! High energy neutron flux (-)
real(r8) :: alpha  = 0.0_r8 ! Ratio of Fast Neutron Creation Factor (Soil to Water), alpha (-)
real(r8) :: L1     = 0.0_r8 ! High Energy Soil   Attenuation Length (g/cm2)
real(r8) :: L2     = 0.0_r8 ! High Energy Water  Attenuation Length (g/cm2)
real(r8) :: L3     = 0.0_r8 ! Fast Neutron Soil  Attenuation Length (g/cm2)
real(r8) :: L4     = 0.0_r8 ! Fast Neutron Water Attenuation Length (g/cm2)
real(r8) :: zdeg
real(r8) :: zrad
real(r8) :: ideg
real(r8) :: costheta
real(r8) :: dtheta
real(r8), dimension(ens_size) :: totflux     ! Total flux of above-ground fast neutrons

real(r8), dimension(:,:), allocatable :: dz          ! Soil layers (cm)
real(r8), dimension(:,:), allocatable :: zthick      ! Soil layer thickness (cm)
real(r8), dimension(:,:), allocatable :: vwc         ! Volumetric Water Content (m3/m3)
real(r8), dimension(:,:), allocatable :: isoimass    ! Integrated dry soil mass above layer (g)
real(r8), dimension(:,:), allocatable :: iwatmass    ! Integrated water mass above layer (g)
real(r8), dimension(:,:), allocatable :: hiflux      ! High energy neutron flux
real(r8), dimension(:,:), allocatable :: fastpot     ! Fast neutron source strength of layer
real(r8), dimension(:,:), allocatable :: h2oeffdens  ! "Effective" density of water in layer (g/cm3)
real(r8), dimension(:,:), allocatable :: idegrad     ! Integrated neutron degradation factor (-)
real(r8), dimension(:,:), allocatable :: fastflux    ! Contribution to above-ground neutron flux

!rr: Not needed for DART
!rr: real(r8), dimension(:), allocatable :: normfast ! Normalized contribution to neutron flux (-) [weighting factors]

real(r8), parameter   :: h2odens = 1000.0_r8 ! Density of water (g/cm3)

integer,  parameter   :: maxlayers = 1000   ! more than the maximum # of model soil layers
real(r8), allocatable :: layerz(:,:)        ! original soil layer depths
real(r8), allocatable :: soil_moisture(:,:) ! original soil layer moistures

integer  :: angle, angledz, maxangle  ! loop indices for an integration interval
integer  :: i, zi, nlevels
real(r8) :: loc_array(3)
real(r8) :: loc_lon, loc_lat 
real(r8) :: loc_value(ens_size)
type(location_type) :: loc
integer :: imem
integer :: loc_value_istatus(ens_size), loc_istatus(ens_size), layerz_istatus(ens_size)
logical :: return_now
character(len=*), parameter :: routine = 'get_expected_radiance'

!=================================================================================

if ( .not. module_initialized ) call initialize_module

val = 0.0_r8 ! set return value early

!=================================================================================
! COSMIC: Site specific-parameters come from the observation metadata
!=================================================================================

! Make sure the desired key is within the length of the metadata arrays.
call key_within_range(key, routine)

bd     = observation_metadata(key)%bd
vwclat = observation_metadata(key)%lattwat
N      = observation_metadata(key)%N
alpha  = observation_metadata(key)%alpha
L1     = observation_metadata(key)%L1
L2     = observation_metadata(key)%L2
L3     = observation_metadata(key)%L3
L4     = observation_metadata(key)%L4

!=================================================================================
! Determine the number of soil layers and their depths
! using only the standard DART interfaces to model
! set the locations for each of the model levels
!=================================================================================

loc_array = get_location(location) ! loc is in DEGREES
loc_lon   = loc_array(1)
loc_lat   = loc_array(2)

nlevels = 0
COUNTLEVELS : do i = 1,maxlayers
   loc = set_location(loc_lon, loc_lat, real(i,r8), VERTISLEVEL)
   call interpolate(state_handle, ens_size, loc, QTY_GEOPOTENTIAL_HEIGHT, loc_value, loc_istatus)
   if ( any(loc_istatus /= 0 ) ) exit COUNTLEVELS
   nlevels = nlevels + 1
enddo COUNTLEVELS

if ((nlevels == maxlayers) .or. (nlevels == 0)) then
   write(string1,*) 'FAILED to determine number of soil layers in model.'
   if (debug) call error_handler(E_MSG,routine,string1,source,revision,revdate)
   istatus = 1
   val     = MISSING_R8
   return 
else
!   write(*,*)routine // 'we have ',nlevels,' model levels'
endif

! Now actually find the depths at each level
! While we're at it, might as well get the soil moisture there, too.

allocate(layerz(ens_size, nlevels),&
         soil_moisture(ens_size, nlevels))

! Set all of the istatuses back to zero for track_status
istatus = 0

FINDLEVELS : do i = 1,nlevels
   loc = set_location(loc_lon, loc_lat, real(i,r8), VERTISLEVEL)
   call interpolate(state_handle, ens_size, loc, QTY_GEOPOTENTIAL_HEIGHT, layerz(:,i), layerz_istatus)
   call track_status(ens_size, layerz_istatus, val, istatus, return_now)
   if (return_now) return

   loc = set_location(loc_lon, loc_lat, layerz(1,i), VERTISHEIGHT)
   call interpolate(state_handle, ens_size, loc, QTY_SOIL_MOISTURE, loc_value, loc_value_istatus)
   call track_status(ens_size, loc_value_istatus, val, istatus, return_now)
   if ( any(loc_value_istatus /=0) ) then
      write(string1,*) 'FAILED to determine soil moisture for layer',i
      where (loc_value_istatus /= 0) istatus = 2
      if (return_now) return
   endif

   soil_moisture(:,i) = loc_value

enddo FINDLEVELS

!rr: DART soil layers are negative and given in meters while COSMIC needs them to be
!rr: positive and in centimeters

if ( all(layerz < 0.0_r8) ) then
   layerz = -1.0_r8 * 100.0_r8 * layerz
elseif ( any(layerz < 0.0_r8) ) then
   write(string1,*) 'unusual values of soil layers in model.'
   call error_handler(E_ERR,routine, string1,source,revision,revdate)
endif

!=================================================================================
! COSMIC: Allocate arrays and initialize variables
!=================================================================================

allocate( dz(ens_size,nlyr),          &
          zthick(ens_size,nlyr),      &
          vwc(ens_size,nlyr),         &
          hiflux(ens_size,nlyr),      &
          fastpot(ens_size,nlyr),     &
          h2oeffdens(ens_size,nlyr),  &
          idegrad(ens_size,nlyr),     &
          fastflux(ens_size,nlyr),    &
          isoimass(ens_size,nlyr),    &
          iwatmass(ens_size,nlyr) )

totflux(:) = 0.0_r8

do i = 1,nlyr

   dz(:,i)         = real(i,r8)/10.0_r8 ! 0.1 cm intervals ... for now.
   zthick(:,i)     = 0.0_r8
   h2oeffdens(:,i) = 0.0_r8
   vwc(:,i)        = 0.0_r8
   isoimass(:,i)   = 0.0_r8
   iwatmass(:,i)   = 0.0_r8
   hiflux(:,i)     = 0.0_r8
   fastpot(:,i)    = 0.0_r8
   idegrad(:,i)    = 0.0_r8
   fastflux(:,i)   = 0.0_r8

    !rr: Not needed for DART
    !rr: normfast(i)    = 0.0_r8

enddo

!=================================================================================
! Get soil moisture from individual model layers and assign them to
! 1 mm intervals (down to 3 meters)

do i = 1,nlyr
   !>@todo JH can we some how fix this logic to not include an imem
   !> loop?  can not use WHERE () when there is an exit.
   do imem = 1,ens_size
     SOIL : do zi = 1,nlevels
        if (dz(imem,i) >= layerz(imem, nlevels)) then
           vwc(imem, i) = soil_moisture(imem, nlevels)
           exit SOIL
        elseif (dz(imem, i) <= layerz(imem,  zi)) then
           vwc(imem, i) = soil_moisture(imem,  zi) ! soil moisture (m3/m3)
           exit SOIL
        endif
     enddo SOIL
   enddo
enddo

deallocate(layerz ,soil_moisture)

!=================================================================================
! COSMIC: Neutron flux calculation
!=================================================================================

! At some point, you might want to tinker around with non-uniform
! soil layer thicknesses.
zthick(:, 1) = dz(:, 1) - 0.0_r8 ! Surface layer
do i = 2,nlyr
   zthick(:, i) = dz(:, i) - dz(:, i-1) ! Remaining layers
enddo

! Angle distribution parameters (HARDWIRED)
!rr: Using 0.5 deg angle intervals appears to be sufficient
!rr: (smaller angles increase the computing time for COSMIC)
ideg     = 0.5_r8                   ! ideg ultimately controls the number of trips through
angledz  = nint(ideg*10.0_r8)       ! the ANGLE loop. Make sure the 10.0_r8 is enough
maxangle = 900 - angledz            ! to create integers with no remainder
dtheta   = ideg*(PI/180.0_r8)

if ( real(angledz,r8) /= ideg*10.0_r8 ) then
   write(string1,*) 'ideg*10.0 must result in an integer - it results in ',ideg*10.0_r8
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
endif

do i = 1,nlyr

   !> @todo put these computations in a routine and call it
   !> only if istatus is ok.

   ! High energy neutron downward flux
   ! The integration is now performed at the node of each layer (i.e., center of the layer)
   h2oeffdens(:,i) = ((vwc(:,i)+vwclat)*h2odens)/1000.0_r8

   if(i > 1) then
      ! Assuming an area of 1 cm2
      isoimass(:,i) = isoimass(:,i-1) + bd*(0.5_r8*zthick(:,i-1))*1.0_r8 + &
                                        bd*(0.5_r8*zthick(:,i  ))*1.0_r8
      ! Assuming an area of 1 cm2
      iwatmass(:,i) = iwatmass(:,i-1) + h2oeffdens(:, i-1)*(0.5_r8*zthick(:, i-1))*1.0_r8 + &
                                        h2oeffdens(:, i  )*(0.5_r8*zthick(:, i  ))*1.0_r8
   else
      isoimass(:, i) =               bd*(0.5_r8*zthick(:, i))*1.0_r8 ! Assuming an area of 1 cm2
      iwatmass(:, i) = h2oeffdens(:, i)*(0.5_r8*zthick(:, i))*1.0_r8 ! Assuming an area of 1 cm2
   endif

   !>@todo FIXME JH: probably will not do pointwise exponentiation
   hiflux( :, i) = N*exp(-(isoimass(:, i)/L1 + iwatmass(:, i)/L2) )
   fastpot(:, i) = zthick(:, i)*hiflux(:, i)*(alpha*bd + h2oeffdens(:, i))

   ! This second loop needs to be done for the distribution of angles for fast neutron release
   ! the intent is to loop from 0 to 89.5 by 0.5 degrees - or similar.
   ! Because Fortran loop indices are integers, we have to divide the indices by 10 - you get the idea.

   do angle=0,maxangle,angledz
      zdeg     = real(angle,r8)/10.0_r8   ! 0.0  0.5  1.0  1.5 ...
      zrad     = (zdeg*PI)/180.0_r8
      costheta = cos(zrad)

      ! Angle-dependent low energy (fast) neutron upward flux
      !>@todo FIXME JH: probably will not do pointwise exponentiation
      fastflux(:, i) = fastflux(:, i) + fastpot(:, i)*exp(-(isoimass(:, i)/L3 + iwatmass(:, i)/L4)/costheta)*dtheta
   enddo

   ! After contribution from all directions are taken into account,
   ! need to multiply fastflux by 2/PI

   fastflux(:, i) = (2.0_r8/PI)*fastflux(:, i)

   ! Low energy (fast) neutron upward flux
   totflux(:) = totflux(:) + fastflux(:, i)
enddo

deallocate(dz, zthick, vwc, hiflux, fastpot, &
           h2oeffdens, idegrad, fastflux, isoimass, iwatmass)

!=================================================================================
! ... and finally set the return the neutron intensity (i.e. totflux)

where (istatus == 0) val = totflux
where (istatus /= 0) val = missing_r8

return
end subroutine get_expected_radiance


!----------------------------------------------------------------------
! this is a fatal error routine.  use with care - if a forward
! operator fails it should return a bad error code, not die.

subroutine check_iostat(istat, routine, varname, msgstring)
integer,          intent(in) :: istat
character(len=*), intent(in) :: routine
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: msgstring

if ( istat /= 0 ) then
   write(string1,*)'istat should be 0 but is ',istat,' for '//varname
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=msgstring)
end if

end subroutine check_iostat


!----------------------------------------------------------------------
! Make sure we are addressing within the metadata arrays

subroutine key_within_range(key, routine)

integer,          intent(in) :: key
character(len=*), intent(in) :: routine

! fine -- no problem.
if ((key > 0) .and. (key <= rttovkey)) return

! Bad news. Tell the user.
write(string1, *) 'key (',key,') not within known range ( 1,', rttovkey,')'
call error_handler(E_ERR,routine,string1,source,revision,revdate)

end subroutine key_within_range


!----------------------------------------------------------------------
! If the allocatable metadata arrays are not big enough ... try again

subroutine grow_metadata(key, routine)

integer,          intent(in) :: key
character(len=*), intent(in) :: routine

integer :: orglength
type(site_metadata), allocatable, dimension(:) :: safe_metadata

! fine -- no problem.
if ((key > 0) .and. (key <= MAXrttovkey)) return

! Check for some error conditions.
if (key < 1) then
   write(string1, *) 'key (',key,') must be >= 1'
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
elseif (key >= 2*MAXrttovkey) then
   write(string1, *) 'key (',key,') really unexpected.'
   write(string2, *) 'doubling storage will not help.'
   call error_handler(E_ERR,routine,string1,source,revision,revdate, &
                      text2=string2)
endif

orglength   = MAXrttovkey
MAXrttovkey = 2 * orglength

! News. Tell the user we are increasing storage.
write(string1, *) 'key (',key,') exceeds Nmax_radiance_obs (',orglength,')'
write(string2, *) 'Increasing Nmax_radiance_obs to ',MAXrttovkey
call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)

allocate(safe_metadata(orglength))
safe_metadata(:) = observation_metadata(:)

deallocate(observation_metadata)
  allocate(observation_metadata(MAXrttovkey))

observation_metadata(1:orglength)              = safe_metadata(:)
observation_metadata(orglength+1:MAXrttovkey) = missing_metadata

deallocate(safe_metadata)

end subroutine grow_metadata

!----------------------------------------------------------------------

end module obs_def_rttov_mod

! END DART PREPROCESS MODULE CODE

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
