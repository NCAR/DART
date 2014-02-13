! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! Note:  This version has a namelist item for the max number of
! gps observations that can be read in, but it is currently commented out.
! Search below for the 'NAMELIST' string to find where to comment the
! code back in.  There are 2 places.  This maximum must be the total of
! all gps obs in all input files, which if you are reading multiple obs_seq
! files (e.g. for the obs_diag program) might be a larger number than 100K.

! BEGIN DART PREPROCESS KIND LIST
! TEMPERATURE,        KIND_TEMPERATURE,        COMMON_CODE
! SPECIFIC_HUMIDITY,  KIND_SPECIFIC_HUMIDITY,  COMMON_CODE
! PRESSURE,           KIND_PRESSURE,           COMMON_CODE
! GPSRO_REFRACTIVITY, KIND_GPSRO
! END DART PREPROCESS KIND LIST


! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_gps_mod, only : get_expected_gpsro_ref_distrib, interactive_gpsro_ref, &
!                              read_gpsro_ref, write_gpsro_ref
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(GPSRO_REFRACTIVITY)
!            call get_expected_gpsro_ref_distrib(state_ens_handle,  location, obs_def%key, expected_obs, istatus)
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
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             file_exist, open_file, close_file, nmlfileunit, &
                             check_namelist_read, find_namelist_in_file, &
                             do_output, do_nml_file, do_nml_term, &
                             ascii_file_format
use     location_mod, only : location_type, set_location, get_location, &
                             write_location, read_location, vert_is_height, &
                             VERTISHEIGHT
use time_manager_mod, only : time_type, read_time, write_time, &
                             set_time, set_time_missing, interactive_time
use  assim_model_mod, only : interpolate_distrib

use     obs_kind_mod, only : KIND_U_WIND_COMPONENT, &
                             KIND_V_WIND_COMPONENT, KIND_SURFACE_PRESSURE, &
                             KIND_TEMPERATURE, KIND_SPECIFIC_HUMIDITY, &
                             KIND_PRESSURE, KIND_GPSRO

use data_structure_mod, only : ensemble_type

implicit none
private

public :: set_gpsro_ref, get_gpsro_ref, write_gpsro_ref, read_gpsro_ref, &
          get_expected_gpsro_ref_distrib, interactive_gpsro_ref

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
! would be to keep a single KIND_GPSRO, but make 2 observation types.
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

!> Distributed version of get_expected_gpsro_ref_distrib
 subroutine get_expected_gpsro_ref_distrib(state_ens_handle,  location, gpskey, ro_ref, istatus)
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

!> @todo this is completely broken
type(ensemble_type), intent(in)  :: state_ens_handle
type(location_type), intent(in)  :: location
integer,             intent(in)  :: gpskey
real(r8),            intent(out) :: ro_ref(:)
integer,             intent(out) :: istatus(:)

! local variables

real(r8) :: nx, ny, nz       ! unit tangent direction of ray at perigee
real(r8) :: xo, yo, zo       ! perigee location in Cartesian coordinate

real(r8) :: height1, lat1, lon1

real(r8), allocatable :: phase(:)
real(r8), allocatable :: delta_phase1(:), delta_phase2(:)
real(r8), allocatable :: ref_perigee(:), ref1(:), ref2(:), dist_to_perigee(:)
real(r8), allocatable :: xx(:), yy(:), zz(:)
integer               :: iter
integer,  allocatable :: istatus0(:), track_status(:)
real(r8), allocatable :: ref00(:)
real(r8)              :: lon, lat, height, obsloc(3)
integer               :: e, ens_size

if ( .not. module_initialized ) call initialize_module

ens_size = state_ens_handle%num_copies -5
allocate(istatus0(ens_size), track_status(ens_size))
allocate(ref00(ens_size), ref_perigee(ens_size), ref1(ens_size), ref2(ens_size))
allocate(xx(ens_size), yy(ens_size), zz(ens_size))
allocate(delta_phase1(ens_size), delta_phase2(ens_size), phase(ens_size))

if ( .not. vert_is_height(location)) then
   write(string1, *) 'vertical location must be height; gps obs key ', gpskey
   call error_handler(E_ERR,'get_expected_gpsro_ref_distrib', string1, &
                      source, revision, revdate)
endif

obsloc   = get_location(location)

lon      = obsloc(1)                       ! degree: 0 to 360
lat      = obsloc(2)                       ! degree: -90 to 90
height   = obsloc(3)                       ! (m)

! calculate refractivity at perigee

call ref_local_distrib(state_ens_handle,  location, height, lat, lon, ref_perigee, istatus0)
! if istatus > 0, the interpolation failed and we should return failure now.

if(all(istatus0 > 0)) then
   istatus = istatus0
   ro_ref = missing_r8
   return
endif

track_status = istatus

choose: if(gps_data(gpskey)%gpsro_ref_form == 'GPSREF') then
    ! use local refractivity
    ro_ref = ref_perigee * 1.0e6      ! in (N-1)*1.0e6, same with obs

else  ! gps_data(gpskey)%gpsro_ref_form == 'GPSEXC'

    ! otherwise, use non_local refractivity(excess phase delay)
    !> @todo I have broken this.

    ! Initialization
    phase = 0.0_r8
    dist_to_perigee =  0.0_r8   ! distance to perigee from a point of the ray

    nx = gps_data(gpskey)%ray_direction(1)
    ny = gps_data(gpskey)%ray_direction(2)
    nz = gps_data(gpskey)%ray_direction(3)

    ! convert location of the perigee from geodetic to Cartesian coordinate
    call geo2carte(height, lat, lon, xo, yo, zo, gps_data(gpskey)%rfict )

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
    do 

       iter = iter + 1
       dist_to_perigee = dist_to_perigee + gps_data(gpskey)%step_size

       !  integrate to one direction of the ray for one step
       ! HK These are now different for each ensemble member
       xx = xo + dist_to_perigee * nx
       yy = yo + dist_to_perigee * ny
       zz = zo + dist_to_perigee * nz

       ! convert the location of the point to geodetic coordinates
       ! height(m), lat, lon(deg)

       do e = 1, ens_size
          !call carte2geo( xx(e), yy(e), zz(e), height1, lat1, lon1, gps_data(gpskey)%rfict )
          !if (height1(e) >= gps_data(gpskey)%ray_top) exit !HK Is height1 the same for all ensemble members? NO, of course not
         call error_handler(E_ERR, 'get_expected_gpsro_ref_distrib', 'broken Helen obs_def_mod.f90')
       enddo

       ! get the refractivity at this ray point(ref00)
       call ref_local_distrib(state_ens_handle,  location, height1, lat1, lon1, ref00, istatus0)
       do e = 1, ens_size
          ! when any point of the ray is problematic, return failure
          if(istatus0(e) > 0) then
            istatus(e) = istatus0(e)
            ro_ref(e) = missing_r8
          endif
       enddo

       track_status = istatus

       ! get the excess phase due to this ray interval
       delta_phase1 = (ref1 + ref00) * gps_data(gpskey)%step_size * 0.5_r8

       ! save the refractivity for integration of next ray interval
       ref1 = ref00

       ! integrate to the other direction of the ray
       xx = xo - dist_to_perigee * nx
       yy = yo - dist_to_perigee * ny
       zz = zo - dist_to_perigee * nz

       call carte2geo( xx(e), yy(e), zz(e), height1, lat1, lon1, gps_data(gpskey)%rfict )

       ! get the refractivity at this ray point(ref00)
       call ref_local_distrib(state_ens_handle,  location, height1, lat1, lon1, ref00, istatus0)
       ! when any point of the ray is problematic, return failure
       do e = 1, ens_size
          if(istatus0(e) > 0) then
            track_status(e) = istatus0(e)
            ro_ref(e) = missing_r8
          endif
       enddo

       ! get the excess phase due to this ray interval
       delta_phase2 = (ref2 + ref00) * gps_data(gpskey)%step_size * 0.5_r8

       ! save the refractivity for integration of next ray interval
       ref2 = ref00

       phase = phase + delta_phase1 + delta_phase2
       ! print*, 'phase= ',  phase, delta_phase1, delta_phase2

    end do

    ! finish the integration of the excess phase along the ray

    ro_ref = phase    ! in m

    ! print*, 'xx = ', lon, lat, height, ro_ref

endif choose

istatus = track_status

! Putting missing_r8 back in for failed ensemble members
do e = 1, ens_size
   if( istatus(e) /= 0 ) ro_ref(e) = missing_r8
   ! if the original height was too high, for example.  do not return a
   ! negative or 0 excess phase or refractivity.
   if (ro_ref(e) == missing_r8 .or. ro_ref(e) <= 0.0_r8) then
      istatus(e) = 5
      ro_ref(e) = missing_r8
   endif
enddo

! ended ok, return local refractivity or non-local excess phase accumulated value
!istatus = 0 !HK This is annoying 
!> @todo check the istatus return 

end subroutine get_expected_gpsro_ref_distrib

!> Distributed version 
subroutine ref_local_distrib(state_ens_handle,  location, height, lat, lon, ref00, istatus0)
!------------------------------------------------------------------------------
!
! Calculate local refractivity at any GPS ray point (height, lat, lon)
!
! inputs:
!    height, lat, lon:  GPS observation location (units: m, degree)
!
! output:
!    ref00: modeled local refractivity at ray point(unit: N-1, ~1.0e-4 to e-6)
!
!------------------------------------------------------------------------------

type(ensemble_type), intent(in)  :: state_ens_handle
real(r8),            intent(in)  :: lon, lat, height
real(r8),            intent(out) :: ref00(:)
integer,             intent(out) :: istatus0(:)

real(r8), parameter::  rd = 287.05_r8, rv = 461.51_r8, c1 = 77.6d-6 , &
                       c2 = 3.73d-1,  rdorv = rd/rv

integer,  allocatable :: track_status(:)
integer               :: e, ens_size
real(r8)              :: lon2
real(r8), allocatable :: t(:), q(:), p(:), tv(:), ew(:)
type(location_type)   :: location, location2
integer               :: which_vert

if ( .not. module_initialized ) call initialize_module

ens_size = state_ens_handle%num_copies -5 ! This could be n
allocate(t(ens_size), q(ens_size), p(ens_size), tv(ens_size), ew(ens_size), track_status(ens_size))

! for integration of GPS ray path beyond the wraparound point
lon2 = lon
if(lon > 360.0_r8 ) lon2 = lon - 360.0_r8
if(lon <   0.0_r8 ) lon2 = lon + 360.0_r8

which_vert = VERTISHEIGHT
location2 = set_location(lon2, lat, height,  which_vert)

! set return values assuming failure, so we can simply return if any
! of the interpolation calls below fail.
istatus0 = 3
ref00 = missing_r8

call interpolate_distrib(location2,  KIND_TEMPERATURE, istatus0,  t,  state_ens_handle)
if (all(istatus0 > 0)) return
track_status = istatus0

call interpolate_distrib(location2, KIND_SPECIFIC_HUMIDITY, istatus0, q, state_ens_handle)
if (all(istatus0 > 0)) return
do e = 1, ens_size
   if ( istatus0(e) /= 0 ) track_status(e) = istatus0(e)
enddo

call interpolate_distrib(location2,  KIND_PRESSURE, istatus0, p, state_ens_handle)
if (all(istatus0 > 0)) return
do e = 1, ens_size
   if ( istatus0(e) /= 0 ) track_status(e) = istatus0(e)
enddo

!  required variable units for calculation of GPS refractivity
!   t :  Kelvin, from top to bottom
!   q :  kg/kg, from top to bottom
!   p :  mb

p     = p * 0.01_r8      ! to mb

tv    = t * (1.0_r8+(rv/rd - 1.0_r8)*q)         ! virtual temperature
ew    = q * p/(rdorv + (1.0_r8-rdorv)*q )
ref00 = c1*p/t + c2*ew/(t**2)              ! (N-1)

istatus0 = track_status

! put failed copies back in
do e = 1, ens_size
   if( istatus0(e) /= 0 ) ref00(e) = missing_r8
enddo

end subroutine ref_local_distrib


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

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
