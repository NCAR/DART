! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_def_gps_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! BEGIN DART PREPROCESS KIND LIST
! GPSRO_REFRACTIVITY,  KIND_GPSRO
! END DART PREPROCESS KIND LIST


! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_gps_mod, only : get_expected_gpsro_ref, read_gpsro_ref, write_gpsro_ref
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(GPSRO_REFRACTIVITY)
!            call get_expected_gpsro_ref(state, location, obs_def%key, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF


! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(GPSRO_REFRACTIVITY)
!         call read_gpsro_ref(obs_def%key, ifile, fileformat)
!         continue
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(GPSRO_REFRACTIVITY)
!         call write_gpsro_ref(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(GPSRO_REFRACTIVITY)
!       ! call interactive_gpsro_ref(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF


use        types_mod, only : r8, missing_r8, RAD2DEG, DEG2RAD, PI
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, file_exist, &
                             open_file, logfileunit, close_file
use     location_mod, only : location_type, set_location, get_location , write_location, &
                             read_location
use time_manager_mod, only : time_type, read_time, write_time, set_time, set_time_missing, &
                             interactive_time
use  assim_model_mod, only : interpolate

use     obs_kind_mod, only : KIND_U_WIND_COMPONENT, &
                             KIND_V_WIND_COMPONENT, KIND_SURFACE_PRESSURE, &
                             KIND_TEMPERATURE, KIND_SPECIFIC_HUMIDITY, KIND_PRESSURE, &
                             KIND_GPSRO

implicit none
private

public :: set_gpsro_ref, write_gpsro_ref, read_gpsro_ref, get_expected_gpsro_ref

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical, save :: module_initialized = .false.

! Storage for the special information required for GPS RO observations
!
integer, parameter :: max_gpsro_obs = 100000
character(len=6) :: gpsro_ref_form(max_gpsro_obs)
real(r8) ::       ray_direction(3, max_gpsro_obs)
real(r8) ::                  rfict(max_gpsro_obs) ! local curvature radius
real(r8) ::              step_size(max_gpsro_obs)
real(r8) ::                ray_top(max_gpsro_obs)
integer  :: ii
integer  :: keycount = 0

contains

!---------------------------------------------------------------------------------


  subroutine initialize_module
!---------------------------------------------------------------------------------
! subroutine initialize_module
!
! pretty simple for this module.

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module



 subroutine set_gpsro_ref(key, nx, ny, nz, rfict0, ds, htop, subset0)
!---------------------------------------------------------------------------------
!subroutine set_gpsro_ref(key, nx, ny, nz, rfict0, ds, htop, subset0)
!

integer,          intent(in) :: key
real(r8),         intent(in) :: nx, ny, nz, rfict0, ds, htop
character(len=6), intent(in) :: subset0

if ( .not. module_initialized ) call initialize_module

ray_direction(1, key) = nx
ray_direction(2, key) = ny
ray_direction(3, key) = nz
gpsro_ref_form(key)   = subset0

    rfict(key) = rfict0
step_size(key) = ds
  ray_top(key) = htop

end subroutine set_gpsro_ref



 subroutine write_gpsro_ref(key, ifile, fform)
!---------------------------------------------------------------------------------
!subroutine write_gpsro_ref(key, ifile, fform)
!

integer,          intent(in)           :: key, ifile
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Write the 5 character identifier for verbose formatted output
! Write out the obs_def key for this observation
SELECT CASE (fileformat)

   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
    write(ifile) key
    write(ifile) rfict(key), step_size(key), ray_top(key), &
                (ray_direction(ii, key), ii=1, 3), gpsro_ref_form(key)
    continue


   CASE DEFAULT
    write(ifile,11)  key
    write(ifile, *) rfict(key), step_size(key), ray_top(key), &
                   (ray_direction(ii, key), ii=1, 3), gpsro_ref_form(key)
END SELECT
11  format('gpsroref', i8)

end subroutine write_gpsro_ref



 subroutine read_gpsro_ref(key, ifile, fform)
!---------------------------------------------------------------------------------
!subroutine read_gpsro_ref(key, ifile, fform)
!
! Every GPS observation has its own (metadata) key.
! When you read two gps observation sequence files, it is necessary to track the
! total number of metadata keys read ... not just the number in the current file.
! 

integer,          intent(out)          :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

integer             :: keyin    ! the metadata key in the current obs sequence

character(len=8)    :: header
character(len=32)   :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

keycount = keycount + 1    ! the total metadata key count from all sequences
key = keycount             ! copied to the output variable

if(key > max_gpsro_obs) then
   write(*, *) 'key (',key,') exceeds max_gpsro_obs (',max_gpsro_obs,')'
   call error_handler(E_ERR,'read_gpsro_ref', &
        'Increase max_gpsro_obs.', source, revision, revdate)
endif

! Read the character identifier for verbose formatted output
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")

    read(ifile) keyin          ! read and throw away
    read(ifile) rfict(key), step_size(key), ray_top(key), &
               (ray_direction(ii, key), ii=1, 3), gpsro_ref_form(key)
    continue

   CASE DEFAULT
      read(ifile, FMT='(a8, i8)') header, keyin    ! throw away keyin
      if(header /= 'gpsroref') then
       call error_handler(E_ERR,'read_gpsro_ref', &
            'Expected header "gpsroref" in input file', source, revision, revdate)
      endif
    read(ifile, *) rfict(key), step_size(key), ray_top(key), &
                  (ray_direction(ii, key), ii=1, 3), gpsro_ref_form(key)

END SELECT

end subroutine read_gpsro_ref



 subroutine get_expected_gpsro_ref(state_vector, location, key, ro_ref, istatus)
!---------------------------------------------------------------------------------
!subroutine get_expected_gpsro_ref(state_vector, location, key, ro_ref, istatus)
!
! Purpose: Calculate GPS RO local refractivity or non_local (integrated) 
!          refractivity (excess phase, Sergey Sokolovskiy et al., 2005)
!---------------------------------------------------------------------------------
!
! inputs:
!    state_vector:    DART state vector
!
! output parameters:
!    ro_ref: modeled local refractivity (N-1)*1.0e6 or non_local refractivity (excess phase, m)
!            (according to the input data parameter subset)
!    istatus:  =0 normal; =1 outside of domain.
!---------------------------------------------------------------------------------
!  Author: Hui Liu 
!  Version 1.1: June 15, 2004: Initial version CAM
!
!  Version 1.2: July 29, 2005: revised for new obs_def and WRF
!---------------------------------------------------------------------------------
implicit none

real(r8),            intent(in)  :: state_vector(:)
integer,             intent(in)  :: key
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: ro_ref
integer,             intent(out) :: istatus

! local variables

real(r8) :: nx, ny, nz       ! unit tangent direction of ray at perigee
real(r8) :: xo, yo, zo       ! perigee location in Cartesian coordinate

real(r8) :: ref_perigee, ref00, ref1, ref2, dist_to_perigee
real(r8) :: phase
real(r8) :: xx, yy, zz, height1, lat1, lon1, delta_phase1, delta_phase2

integer  :: iter, istatus0
real(r8) :: lon, lat, height, obsloc(3)

if ( .not. module_initialized ) call initialize_module

obsloc   = get_location(location)

lon      = obsloc(1)                       ! degree: 0 to 360
lat      = obsloc(2)                       ! degree: -90 to 90
height   = obsloc(3)                       ! (m)

! calculate refractivity at perigee

call ref_local(state_vector, location, height, lat, lon, ref_perigee, istatus0)
istatus = istatus0

choose: if(gpsro_ref_form(key) == 'GPSREF') then ! use local refractivity

    ro_ref = ref_perigee * 1.0e6      ! in (N-1)*1.0e6, same with obs

    if(istatus > 1) then
      istatus = 1
      ro_ref = missing_r8
    endif

    return

else  ! gpsro_ref_form(key) == 'GPSEXC'

    ! otherwise, use non_local refractivity(excess phase delay)

    ! Initilization
    phase = 0.0_r8  
    dist_to_perigee =  0.0_r8   ! distance to perigee from a point of the ray

    nx = ray_direction(1, key)
    ny = ray_direction(2, key)
    nz = ray_direction(3, key)

    ! convert location of the perigee from geodetic to Cartesian coordinate

    call geo2carte (height, lat, lon, xo, yo, zo, rfict(key) )

    ! currently, use a straight line passing the perigee point as ray model.
    ! later, more sophisticated ray models can be used.
    !
    ! Start the horizontal integrate of the model refractivity along a 
    ! straight line path in cartesian coordinate
    !
    ! (x-xo)/a = (y-yo)/b = (z-zo)/c,  (a,b,c) is the line direction

    ref1 = ref_perigee
    ref2 = ref_perigee

    ! step_size = 5.0_r8      ! (m)

    iter = 0
100 continue

    iter = iter + 1
    dist_to_perigee = dist_to_perigee + step_size(key)

    !  integrate to one direction of the ray for one step
    xx = xo + dist_to_perigee * nx
    yy = yo + dist_to_perigee * ny
    zz = zo + dist_to_perigee * nz
   
    ! convert the location of the point to geodetic coordinates 
    ! height(m), lat, lon(deg)

    call carte2geo(xx, yy, zz, height1, lat1, lon1, rfict(key) )  

    ! get the refractivity at this ray point(ref00)
    call ref_local(state_vector, location, height1, lat1, lon1, ref00, istatus0)
    istatus = istatus + istatus0

    ! get the excess phase due to this ray interval
    delta_phase1 = (ref1 + ref00) * step_size(key) * 0.5_r8

    ! save the refractivity for integration of next ray interval
    ref1 = ref00

    ! integrate to the other direction of the ray
    xx = xo - dist_to_perigee * nx
    yy = yo - dist_to_perigee * ny 
    zz = zo - dist_to_perigee * nz
   
    call carte2geo (xx, yy, zz, height1, lat1, lon1, rfict(key) )  

    ! get the refractivity at this ray point(ref00)
    call ref_local(state_vector, location, height1, lat1, lon1, ref00, istatus0)
    istatus = istatus + istatus0

    ! get the excess phase due to this ray interval
    delta_phase2 = (ref2 + ref00) * step_size(key) * 0.5_r8

    ! save the refractivity for integration of next ray interval
    ref2 = ref00

    phase = phase + delta_phase1 + delta_phase2
    ! print*, 'phase= ',  phase, delta_phase1, delta_phase2

    if(height1 .lt. ray_top(key) ) go to 100             !! do one more step

    ! finish the integration of the excess phase along the ray

    ro_ref = phase    ! in m

    ! when any point of the ray is problematic
    if(istatus > 1) then
      istatus = 1
      ro_ref = missing_r8
    endif

    ! print*, 'xx = ', lon, lat, height, ro_ref

endif choose

end subroutine get_expected_gpsro_ref



 subroutine ref_local(state_vector, location, height, lat, lon, ref00, istatus0)
!---------------------------------------------------------------------------------
!subroutine ref_local(state_vector, location, height, lat, lon, ref00, istatus0)
!
!   Calculate local refractivity at any GPS ray point (height, lat, lon)
!
!  inputs:
!         height, lat, lon:  GPS observation location (units: m, degree)
!
!  output:
!         ref00: modeled local refractivity at ray point(unit: N-1, ~1.0e-4 to e-6)
!
!---------------------------------------------------------------------------------
implicit none

real(r8), intent(in) :: state_vector(:)
real(r8), intent(in) :: lon, lat, height

real(r8), intent(out) :: ref00
integer,  intent(out) :: istatus0

real(r8), parameter::  rd = 287.05_r8, rv = 461.51_r8, c1 = 77.6d-6 , &
                       c2 = 3.73d-1,  rdorv = rd/rv
real(r8) :: lon2, t, q, p, tv, ew
type(location_type) :: location, location2
integer :: which_vert

if ( .not. module_initialized ) call initialize_module

! for integration of GPS ray path beyond the wraparound point
lon2 = lon
if(lon > 360.0_r8 ) lon2 = lon - 360.0_r8
if(lon <   0.0_r8 ) lon2 = lon + 360.0_r8

which_vert = 3
location2 = set_location(lon2, lat, height,  which_vert)

istatus0 = 0
call interpolate(state_vector, location2,  KIND_TEMPERATURE,       t, istatus0)
call interpolate(state_vector, location2,  KIND_SPECIFIC_HUMIDITY, q, istatus0)
call interpolate(state_vector, location2,  KIND_PRESSURE,          p, istatus0)

!  required variable units for calculation of GPS refractivity
!   t :  Kelvin, from top to bottom
!   q :  kg/kg, from top to bottom
!   p :  mb

p     = p * 0.01_r8      ! to mb

tv    = t * (1.0_r8+(rv/rd - 1.0_r8)*q)         ! virtual temperature
ew    = q * p/(rdorv + (1.0_r8-rdorv)*q )
ref00 = c1*p/t + c2*ew/(t**2)              ! (N-1)

end subroutine ref_local


 subroutine geo2carte (s1, s2, s3, x1, x2, x3, rfict0) 
!---------------------------------------------------------------------------------
!subroutine geo2carte (s1, s2, s3, x1, x2, x3, rfict0) 
!
!  Converts geodetical coordinates to cartesian with a reference sphere
!---------------------------------------------------------------------------------
!  input parameters:
!   s - geodetical coordinates
!        (height (m), latitude (degree), longitude (degree))
!                     -90 to 90           0 to 360
!  output parameters:
!   x - cartesian coordinates (m) connected with the earth(x, y, z-coordinate)
!---------------------------------------------------------------------------------
implicit none
real(r8), intent(in)  :: s1, s2, s3, rfict0    ! units: m
real(r8), intent(out) ::   x1, x2 ,x3
real(r8) :: g3, g4

if ( .not. module_initialized ) call initialize_module

g3 = s1 + rfict0
g4 = g3 * cos(s2*DEG2RAD) 
x1 = g4 * cos(s3*DEG2RAD)
x2 = g4 * sin(s3*DEG2RAD)
x3 = g3 * sin(s2*DEG2RAD)

end subroutine geo2carte


 subroutine carte2geo (x1, x2, x3, s1, s2, s3, rfict0)
!---------------------------------------------------------------------------------
!subroutine carte2geo (x1, x2, x3, s1, s2, s3, rfict0)
!
!  Converts cartesian coordinates to geodetical.
!
!   input parameters:
!        x - cartesian coordinates (x, y, z-coordinate, unit: m)
!
!   output parameters:
!        s - geodetical coordinates
!            (height (m), latitude (deg), longitude (deg))
!                          -90 to 90         0 to 360
!---------------------------------------------------------------------------------
implicit none
real(r8), intent(in)  :: x1, x2, x3, rfict0
real(r8), intent(out) :: s1, s2, s3

real(r8), parameter :: crcl  = 2.0_r8 * PI, &
                       crcl2 = 4.0_r8 * PI

real(r8) :: rho, sphi, azmth

if ( .not. module_initialized ) call initialize_module

rho   = sqrt (x1**2 + x2**2 + x3**2 ) 
sphi  = x3/rho
s1    = rho - rfict0
s2    = asin (sphi) 
azmth = atan2 (x2, x1)
s3    = mod((azmth + crcl2), crcl)

s2    = s2 * RAD2DEG
s3    = s3 * RAD2DEG

end  subroutine carte2geo

end module obs_def_gps_mod
