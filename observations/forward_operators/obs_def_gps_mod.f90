! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! Note:  This version has a namelist item for the max number of
! gps observations that can be read in, but it is currently commented out.
! Search below for the 'NAMELIST' string to find where to comment the
! code back in.  There are 2 places.  This maximum must be the total of
! all gps obs in all input files, which if you are reading multiple obs_seq
! files (e.g. for the obs_diag program) might be a larger number than 100K.

!>@todo we should have a local vs nonlocal forward operator for GPS RO,
!>so we don't have to add the metadata for the local operator.  big space
!>and time savings.  also, we should add GPSRO_BENDING_ANGLE if someone
!>can contribute a forward operator for it.

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! TEMPERATURE,             QTY_TEMPERATURE,        COMMON_CODE
! SPECIFIC_HUMIDITY,       QTY_SPECIFIC_HUMIDITY,  COMMON_CODE
! PRESSURE,                QTY_PRESSURE,           COMMON_CODE
! GPSRO_REFRACTIVITY,      QTY_GPSRO
! END DART PREPROCESS TYPE DEFINITIONS


! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_gps_mod, only : get_expected_gpsro_ref, interactive_gpsro_ref, &
!                              read_gpsro_ref, write_gpsro_ref
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(GPSRO_REFRACTIVITY)
!            call get_expected_gpsro_ref(state_handle, ens_size,  location, obs_def%key, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF


! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(GPSRO_REFRACTIVITY)
!            call read_gpsro_ref(obs_def%key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(GPSRO_REFRACTIVITY)
!            call write_gpsro_ref(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(GPSRO_REFRACTIVITY)
!            call interactive_gpsro_ref(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_gps_mod

use        types_mod, only : r8, missing_r8, RAD2DEG, DEG2RAD, PI
use    utilities_mod, only : register_module, error_handler, E_ERR, &
                             nmlfileunit, check_namelist_read,      &
                             find_namelist_in_file, do_nml_file, do_nml_term, &
                             ascii_file_format
use     location_mod, only : location_type, set_location, get_location, &
                             is_vertical, &
                             VERTISHEIGHT
use  assim_model_mod, only : interpolate

use     obs_kind_mod, only : QTY_TEMPERATURE, QTY_SPECIFIC_HUMIDITY, &
                             QTY_PRESSURE

use  ensemble_manager_mod, only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private

public :: set_gpsro_ref, get_gpsro_ref, write_gpsro_ref, read_gpsro_ref, &
          get_expected_gpsro_ref, interactive_gpsro_ref

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

! Storage for the special information required for GPS RO observations
!

! Because we are currently only generating one observation type
! (GPSRO_REFRACTIVITY), there must be enough of these to cover all gps
! obs in all obs_seq files that are read in (e.g. for obs_diag if you
! cover multiple days or weeks, you must have enough room for all of them.)
! the local operator needs none of this additional info; the best approach
! would be to keep a single QTY_GPSRO, but make 2 observation types.
! the local has no additional metadata; the nonlocal needs one of these
! allocated and filled in.
integer :: max_gpsro_obs = 100000

type gps_nonlocal_type
   private
   character(len=6) :: gpsro_ref_form
   real(r8)         :: ray_direction(3)
   real(r8)         :: rfict
   real(r8)         :: step_size
   real(r8)         :: ray_top
end type gps_nonlocal_type

type(gps_nonlocal_type), allocatable :: gps_data(:)

namelist /obs_def_gps_nml/ max_gpsro_obs

character(len=129) :: string1, string2
integer  :: ii
integer  :: keycount 

contains

!------------------------------------------------------------------------------


  subroutine initialize_module
!------------------------------------------------------------------------------
!
! initialize global gps private key number and allocate space for obs data
integer :: rc, iunit


call register_module(source, revision, revdate)
module_initialized = .true.

! global count of all gps observations from any input file
keycount = 0

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_def_gps_nml", iunit)
read(iunit, nml = obs_def_gps_nml, iostat = rc)
call check_namelist_read(iunit, rc, "obs_def_gps_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_def_gps_nml)
if (do_nml_term()) write(     *     , nml=obs_def_gps_nml)

! find max number of gps obs which can be stored, and initialize type
allocate(gps_data(max_gpsro_obs), stat = rc)
if (rc /= 0) then
   write(string1, *) 'initial allocation failed for gps observation data,', &
                       'itemcount = ', max_gpsro_obs
   call error_handler(E_ERR,'initialize_module', string1, &
                      source, revision, revdate)
endif

end subroutine initialize_module



 subroutine set_gpsro_ref(gpskey, nx, ny, nz, rfict0, ds, htop, subset0)
!------------------------------------------------------------------------------
!
! increment key and set all private data for this observation

integer,          intent(out) :: gpskey
real(r8),         intent(in)  :: nx, ny, nz, rfict0, ds, htop
character(len=6), intent(in)  :: subset0

if ( .not. module_initialized ) call initialize_module

keycount = keycount + 1
gpskey = keycount

if(gpskey > max_gpsro_obs) then
   write(string1, *) 'key (',gpskey,') exceeds max_gpsro_obs (',max_gpsro_obs,')'
   string2 = 'Increase max_gpsro_obs in input.nml &obs_def_gps_nml namelist.'
   call error_handler(E_ERR,'read_gpsro_ref', string1, &
                      source, revision, revdate, text2=string2)
endif

gps_data(gpskey)%ray_direction(1) = nx
gps_data(gpskey)%ray_direction(2) = ny
gps_data(gpskey)%ray_direction(3) = nz
gps_data(gpskey)%gpsro_ref_form   = subset0

gps_data(gpskey)%rfict     = rfict0
gps_data(gpskey)%step_size = ds
gps_data(gpskey)%ray_top   = htop

end subroutine set_gpsro_ref


 subroutine get_gpsro_ref(gpskey, nx, ny, nz, rfict0, ds, htop, subset0)
!------------------------------------------------------------------------------
!
! return all private data for this observation

integer,          intent(in)  :: gpskey
real(r8),         intent(out) :: nx, ny, nz, rfict0, ds, htop
character(len=6), intent(out) :: subset0

if ( .not. module_initialized ) call initialize_module

if (gpskey < 1 .or. gpskey > keycount) then
   write(string1, *) 'key (',gpskey,') out of valid range (1<=key<=',keycount,')'
   call error_handler(E_ERR,'get_gpsro_ref', string1, &
                      source, revision, revdate)
endif

nx = gps_data(gpskey)%ray_direction(1)
ny = gps_data(gpskey)%ray_direction(2)
nz = gps_data(gpskey)%ray_direction(3)
subset0 = gps_data(gpskey)%gpsro_ref_form

rfict0 = gps_data(gpskey)%rfict
ds     = gps_data(gpskey)%step_size
htop   = gps_data(gpskey)%ray_top

end subroutine get_gpsro_ref



 subroutine write_gpsro_ref(gpskey, ifile, fform)
!------------------------------------------------------------------------------
!

integer,          intent(in)           :: gpskey, ifile
character(len=*), intent(in), optional :: fform


if ( .not. module_initialized ) call initialize_module

! Write the 5 character identifier for verbose formatted output
! Write out the obs_def key for this observation
if (ascii_file_format(fform)) then
   write(ifile,11) gpskey
   write(ifile, *) gps_data(gpskey)%rfict, gps_data(gpskey)%step_size, &
                   gps_data(gpskey)%ray_top, &
                  (gps_data(gpskey)%ray_direction(ii), ii=1, 3), &
                   gps_data(gpskey)%gpsro_ref_form
11  format('gpsroref', i8)
else
   write(ifile) gpskey
   write(ifile) gps_data(gpskey)%rfict, gps_data(gpskey)%step_size, &
                gps_data(gpskey)%ray_top, &
               (gps_data(gpskey)%ray_direction(ii), ii=1, 3), &
                gps_data(gpskey)%gpsro_ref_form
endif

end subroutine write_gpsro_ref



 subroutine read_gpsro_ref(gpskey, ifile, fform)
!------------------------------------------------------------------------------
!
! Every GPS observation has its own (metadata) gpskey.
! When you read multiple gps observation sequence files, it is necessary 
! to track the total number of metadata gpskeys read, not just the number 
! in the current file.
! 

integer,          intent(out)          :: gpskey
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

integer :: keyin    ! the metadata key in the current obs sequence

real(r8) :: nx, ny, nz, rfict0, ds, htop
character(len=6) :: subset0
character(len=8) :: header

if ( .not. module_initialized ) call initialize_module

if (ascii_file_format(fform)) then
   read(ifile, FMT='(a8, i8)') header, keyin    ! throw away keyin
   if(header /= 'gpsroref') then
       call error_handler(E_ERR,'read_gpsro_ref', &
       'Expected header "gpsroref" in input file', source, revision, revdate)
   endif
   read(ifile, *) rfict0, ds, htop, nx, ny, nz, subset0
else
   read(ifile) keyin          ! read and throw away
   read(ifile) rfict0, ds, htop, nx, ny, nz, subset0
endif


! increment key and set all private data for this observation
call set_gpsro_ref(gpskey, nx, ny, nz, rfict0, ds, htop, subset0)

end subroutine read_gpsro_ref


subroutine interactive_gpsro_ref(gpskey)
!----------------------------------------------------------------------
!
! Interactively prompt for the info needed to create a gps refractivity 
! observation.  Increments the key number and returns it.

integer, intent(out) :: gpskey

real(r8) :: nx, ny, nz, rfict0, ds, htop
character(len=6) :: subset0
integer :: gpstype


if ( .not. module_initialized ) call initialize_module

!Now interactively obtain reflectivity type information
! valid choices are local or non-local

write(*, *)
write(*, *) 'Beginning to inquire information on reflectivity type.'
write(*, *)

100 continue
write(*, *) 'Enter 1 for local refractivity (GPSREF)'
write(*, *) 'Enter 2 for non-local refractivity/excess phase delay (GPSEXC)'
write(*, *)

read(*,*) gpstype

select case (gpstype)
   case (1)
      subset0 = 'GPSREF'
   case (2)
      subset0 = 'GPSEXC'
   case default
      write(*,*) 'Bad value, must enter 1 or 2'
      goto 100
end select

if (gpstype == 2) then
    ! FIXME:  i have no idea what valid values are for any
   !  of the following items, so i cannot add any error checking or
   !  guidance for the user.

   write(*, *)
   write(*, *) 'Enter X, Y, Z value for ray direction'
   write(*, *)
   read(*,*) nx, ny, nz

   write(*, *)
   write(*, *) 'Enter local curvature radius'
   write(*, *)
   read(*,*) rfict0

   write(*, *)
   write(*, *) 'Enter step size'
   write(*, *)
   read(*,*) ds

   write(*, *)
   write(*, *) 'Enter ray top'
   write(*, *)
   read(*,*) htop
else
   nx = 0.0
   ny = 0.0
   nz = 0.0
   rfict0 = 0.0
   ds = 0.0
   htop = 0.0
endif

! increment key and set all private data for this observation
call set_gpsro_ref(gpskey, nx, ny, nz, rfict0, ds, htop, subset0)

write(*, *)
write(*, *) 'End of specialized section for gps observation data.'
write(*, *) 'You will now have to enter the regular obs information.'
write(*, *)

end subroutine interactive_gpsro_ref

!> Distributed version of get_expected_gpsro_ref
 subroutine get_expected_gpsro_ref(state_handle, ens_size,  location, gpskey, ro_ref, istatus)
!------------------------------------------------------------------------------
!
! Purpose: Calculate GPS RO local refractivity or non_local (integrated) 
!          refractivity (excess phase, Sergey Sokolovskiy et al., 2005)
!------------------------------------------------------------------------------
!
! inputs:
!    state_vector:    DART state vector
!
! output parameters:
!    ro_ref: modeled local refractivity (N-1)*1.0e6 or non_local 
!            refractivity (excess phase, m)
!            (according to the input data parameter subset)
!    istatus:  =0 normal; =1 outside of domain.
!------------------------------------------------------------------------------
!  Author: Hui Liu 
!  Version 1.1: June 15, 2004: Initial version CAM
!
!  Version 1.2: July 29, 2005: revised for new obs_def and WRF
!------------------------------------------------------------------------------

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: gpskey
real(r8),            intent(out) :: ro_ref(ens_size)
integer,             intent(out) :: istatus(ens_size)

! local variables.  first set is per ensemble member.
! second set is independent of a particular ensemble member.

real(r8) :: phase(ens_size)
real(r8) :: delta_phase1(ens_size), delta_phase2(ens_size)
real(r8) :: ref_perigee(ens_size), ref1(ens_size), ref2(ens_size)
real(r8) :: ref00(ens_size)
integer  :: this_istatus(ens_size)

real(r8) :: nx, ny, nz       ! unit tangent direction of ray at perigee
real(r8) :: xo, yo, zo       ! perigee location in Cartesian coordinate
real(r8) :: height1, lat1, lon1
real(r8) :: xx, yy, zz, dist_to_perigee
real(r8) :: lon, lat, height, obsloc(3)
integer  :: iter
logical  :: return_now

if ( .not. module_initialized ) call initialize_module

if ( .not. is_vertical(location, "HEIGHT")) then
   write(string1, *) 'vertical location must be height; gps obs key ', gpskey
   call error_handler(E_ERR,'get_expected_gpsro_ref', string1, &
                      source, revision, revdate)
endif

obsloc   = get_location(location)


lon      = obsloc(1)                       ! degree: 0 to 360
lat      = obsloc(2)                       ! degree: -90 to 90
height   = obsloc(3)                       ! (m)

! to use track_status() start out with istatus all success.
istatus = 0

! calculate refractivity at perigee

call ref_local(state_handle, ens_size,  lat, lon, height, ref_perigee, this_istatus)
call track_status(ens_size, this_istatus, ro_ref, istatus, return_now)
if (return_now) return

!>@todo we need to add a different obs type and kind for the local operator,
!> so we can avoid allocating metadata because it isn't needed *at all*
!> for the local operator.  nsc 30oct2015
choosetype: if(gps_data(gpskey)%gpsro_ref_form == 'GPSREF') then
    ! use local refractivity
    where (istatus == 0) ro_ref = ref_perigee * 1.0e6      ! in (N-1)*1.0e6, same with obs

else  ! gps_data(gpskey)%gpsro_ref_form == 'GPSEXC'

    ! otherwise, use non_local refractivity(excess phase delay)

    ! Initialization
    phase = 0.0_r8  
    dist_to_perigee =  0.0_r8   ! distance to perigee from a point of the ray

    nx = gps_data(gpskey)%ray_direction(1)
    ny = gps_data(gpskey)%ray_direction(2)
    nz = gps_data(gpskey)%ray_direction(3)

    ! convert location of the perigee from geodetic to Cartesian coordinate

    call geo2carte (height, lat, lon, xo, yo, zo, gps_data(gpskey)%rfict )

    ! currently, use a straight line passing the perigee point as ray model.
    ! later, more sophisticated ray models can be used.
    !
    ! Start the horizontal integrate of the model refractivity along a 
    ! straight line path in cartesian coordinate
    !
    ! (x-xo)/a = (y-yo)/b = (z-zo)/c,  (a,b,c) is the line direction

    ref1 = ref_perigee
    ref2 = ref_perigee

    iter = 0
    INTEGRATE: do 

       iter = iter + 1
       dist_to_perigee = dist_to_perigee + gps_data(gpskey)%step_size
   
       !  integrate to one direction of the ray for one step
       ! HK These are now different for each ensemble member
       xx = xo + dist_to_perigee * nx
       yy = yo + dist_to_perigee * ny
       zz = zo + dist_to_perigee * nz
      
       ! convert the location of the point to geodetic coordinates 
       ! height(m), lat, lon(deg)
   
       call carte2geo(xx, yy, zz, height1, lat1, lon1, gps_data(gpskey)%rfict )  
       if (height1 >= gps_data(gpskey)%ray_top) exit INTEGRATE
   
       ! get the refractivity at this ray point(ref00)
       call ref_local(state_handle, ens_size, lat1, lon1, height1, ref00, this_istatus)
       call track_status(ens_size, this_istatus, ro_ref, istatus, return_now)
       if (return_now) return
   
       ! get the excess phase due to this ray interval
       where(istatus == 0) delta_phase1 = (ref1 + ref00) * gps_data(gpskey)%step_size * 0.5_r8
   
       ! save the refractivity for integration of next ray interval
       ref1 = ref00
   
       ! integrate to the other direction of the ray
       xx = xo - dist_to_perigee * nx
       yy = yo - dist_to_perigee * ny 
       zz = zo - dist_to_perigee * nz
      
       call carte2geo (xx, yy, zz, height1, lat1, lon1, gps_data(gpskey)%rfict )  
   
       ! get the refractivity at this ray point(ref00)
       call ref_local(state_handle, ens_size, lat1, lon1, height1, ref00, this_istatus)
       call track_status(ens_size, this_istatus, ro_ref, istatus, return_now)
       if (return_now) return
   
       ! get the excess phase due to this ray interval
       where(istatus == 0) delta_phase2 = (ref2 + ref00) * gps_data(gpskey)%step_size * 0.5_r8
   
       ! save the refractivity for integration of next ray interval
       ref2 = ref00
   
       where(istatus == 0) phase = phase + delta_phase1 + delta_phase2
       ! print*, 'phase= ',  phase, delta_phase1, delta_phase2

    end do INTEGRATE

    ! finish the integration of the excess phase along the ray

    where(istatus == 0) ro_ref = phase    ! in m

    ! print*, 'xx = ', lon, lat, height, ro_ref

endif choosetype

! make sure return is missing_r8 if failure.
!>@todo is the first line necessary?  i believe the second one
!> is a necessary test.
where (istatus /= 0) ro_ref = missing_r8
where (istatus == 0 .and. ro_ref < 0.0_r8)
   istatus = 5
   ro_ref = missing_r8
endwhere

end subroutine get_expected_gpsro_ref

!> Distributed version.  i removed the unused 'location' argument.
! i changed this to be lat, lon, height - which is a more logical layout of
! the location information.  it used to be height, lat, lon - be careful if
! anyone else calls this routine.  it isn't public so no code outside this
! module can be calling it.   nsc 30 oct 2015

subroutine ref_local(state_handle, ens_size, lat, lon, height, ref00, istatus0)
!------------------------------------------------------------------------------
!
! Calculate local refractivity at any GPS ray point (lat, lon, height)
!
! inputs:
!    lat, lon, height:  GPS observation location (units: degrees, degrees, meters)
!
! output:
!    ref00: modeled local refractivity at ray point(unit: N-1, ~1.0e-4 to e-6)
!
!------------------------------------------------------------------------------

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
real(r8), intent(in) :: lon, lat, height
real(r8),            intent(out) :: ref00(ens_size)
integer,             intent(out) :: istatus0(ens_size)

real(r8), parameter::  rd = 287.05_r8, rv = 461.51_r8, c1 = 77.6d-6 , &
                       c2 = 3.73d-1,  rdorv = rd/rv

real(r8) :: lon2
real(r8) :: t(ens_size), q(ens_size), p(ens_size), ew(ens_size)
integer  :: this_istatus(ens_size)
logical  :: return_now
type(location_type) :: location

if ( .not. module_initialized ) call initialize_module

! for integration of GPS ray path beyond the wraparound point
lon2 = lon
if(lon > 360.0_r8 ) lon2 = lon - 360.0_r8
if(lon <   0.0_r8 ) lon2 = lon + 360.0_r8

location = set_location(lon2, lat, height,  VERTISHEIGHT)

! for sequential sets of calls to interpolate, start with a 0 istatus
! and each call to track_status() will set new failures.
istatus0 = 0
ref00 = missing_r8

call interpolate(state_handle, ens_size, location,  QTY_TEMPERATURE, t, this_istatus)
call track_status(ens_size, this_istatus, ref00, istatus0, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_SPECIFIC_HUMIDITY, q, this_istatus)
call track_status(ens_size, this_istatus, ref00, istatus0, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location,  QTY_PRESSURE, p, this_istatus)
call track_status(ens_size, this_istatus, ref00, istatus0, return_now)
if (return_now) return

!  required variable units for calculation of GPS refractivity
!   t :  Kelvin, from top to bottom
!   q :  kg/kg, from top to bottom
!   p :  mb

where (istatus0 == 0) 
p     = p * 0.01_r8      ! to mb

ew    = q * p/(rdorv + (1.0_r8-rdorv)*q )
ref00 = c1*p/t + c2*ew/(t**2)              ! (N-1)
endwhere

end subroutine ref_local


 subroutine geo2carte (s1, s2, s3, x1, x2, x3, rfict0) 
!------------------------------------------------------------------------------
!
!  Converts geodetical coordinates to cartesian with a reference sphere
!------------------------------------------------------------------------------
!  input parameters:
!   s - geodetical coordinates
!        (height (m), latitude (degree), longitude (degree))
!                     -90 to 90           0 to 360
!  output parameters:
!   x - cartesian coordinates (m) connected with the earth(x, y, z-coordinate)
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
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

! END DART PREPROCESS MODULE CODE

