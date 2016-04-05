! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! Forward operator to compute total precipitable water in a column,
! in centimeters, over the ocean.   Can be used as an example of a
! forward operator that loops over either fixed pressure levels or
! over model levels.

! This code is only correct for TPW over the ocean where the surface
! is at 0m elevation.  For TWP over land you must be able to compute
! the surface elevation for any given lat/lon location.  This code
! assumes the model can return the surface pressure for any given
! lat/lon location, and the specific humidity for a given location
! where the vertical is either pressure or model levels.

! keep in mind that fortran allows only 31 characters in parameter
! definitions (which is what this string is going to be used for).
! if the platform name gets longer than 5 chars, consider going
! to something like xxx_TOTAL_PRECIP_WATER to give you room to
! put in more descriptive platform names.

! BEGIN DART PREPROCESS KIND LIST
!  AQUA_TOTAL_PRECIPITABLE_WATER, KIND_TOTAL_PRECIPITABLE_WATER
! TERRA_TOTAL_PRECIPITABLE_WATER, KIND_TOTAL_PRECIPITABLE_WATER
!  AMSR_TOTAL_PRECIPITABLE_WATER, KIND_TOTAL_PRECIPITABLE_WATER
! MODIS_TOTAL_PRECIPITABLE_WATER, KIND_TOTAL_PRECIPITABLE_WATER
!   GPS_TOTAL_PRECIPITABLE_WATER, KIND_TOTAL_PRECIPITABLE_WATER
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_tpw_mod, only : get_expected_tpw_distrib
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(AQUA_TOTAL_PRECIPITABLE_WATER,TERRA_TOTAL_PRECIPITABLE_WATER,GPS_TOTAL_PRECIPITABLE_WATER)
!            call get_expected_tpw_distrib(state_ens_handle, location, expected_obs, istatus)
!         case(AMSR_TOTAL_PRECIPITABLE_WATER,MODIS_TOTAL_PRECIPITABLE_WATER)
!            call get_expected_tpw_distrib(state_ens_handle, location, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(AQUA_TOTAL_PRECIPITABLE_WATER,TERRA_TOTAL_PRECIPITABLE_WATER,GPS_TOTAL_PRECIPITABLE_WATER)
!           continue
!         case(AMSR_TOTAL_PRECIPITABLE_WATER,MODIS_TOTAL_PRECIPITABLE_WATER)
!           continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(AQUA_TOTAL_PRECIPITABLE_WATER,TERRA_TOTAL_PRECIPITABLE_WATER,GPS_TOTAL_PRECIPITABLE_WATER)
!           continue
!         case(AMSR_TOTAL_PRECIPITABLE_WATER,MODIS_TOTAL_PRECIPITABLE_WATER)
!           continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(AQUA_TOTAL_PRECIPITABLE_WATER,TERRA_TOTAL_PRECIPITABLE_WATER,GPS_TOTAL_PRECIPITABLE_WATER)
!           continue
!         case(AMSR_TOTAL_PRECIPITABLE_WATER,MODIS_TOTAL_PRECIPITABLE_WATER)
!           continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_tpw_mod

use        types_mod, only : r8, missing_r8, RAD2DEG, DEG2RAD, PI
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             nmlfileunit, do_nml_file, do_nml_term, &
                             check_namelist_read, find_namelist_in_file
use     location_mod, only : location_type, set_location, get_location, &
                             write_location, read_location, &
                             VERTISLEVEL, VERTISPRESSURE, VERTISSURFACE
use time_manager_mod, only : time_type, read_time, write_time, &
                             set_time, set_time_missing
use  assim_model_mod, only : interpolate_distrib
use     obs_kind_mod, only : KIND_SURFACE_PRESSURE, KIND_SPECIFIC_HUMIDITY, &
                             KIND_TOTAL_PRECIPITABLE_WATER, KIND_PRESSURE
use  data_structure_mod, only : ensemble_type, copies_in_window


implicit none
private

public ::  get_expected_tpw_distrib

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

character(len=129) :: msgstring

real(r8), parameter :: gravity = 9.81_r8     ! gravitational acceleration (m s^-2)
real(r8), parameter :: density = 1000.0_r8   ! water density in kg/m^3

integer :: max_pressure_intervals = 1000   ! increase as needed

! default samples the atmosphere between the surface and 200 hPa 
! at the model level numbers.  if model_levels is set false,
! then the default samples at 40 heights, evenly divided in
! linear steps in pressure between the surface and top.

logical  :: model_levels = .true.        ! if true, use model levels, ignores num_pres_int
real(r8) :: pressure_top = 20000.0       ! top pressure in pascals
logical  :: separate_surface_level = .true.  ! false: level 1 of 3d grid is sfc
                                             ! true: sfc is separate from 3d grid
integer  :: num_pressure_intervals = 40  ! number of intervals if model_levels is F


namelist /obs_def_tpw_nml/ model_levels, pressure_top,  &
                           separate_surface_level, num_pressure_intervals

contains

!------------------------------------------------------------------------------
subroutine initialize_module()

! should be called once by the other routines in this file to be sure
! the namelist has been read and any initialization code is run.

integer :: iunit, rc

if (module_initialized) return

call register_module(source, revision, revdate)
module_initialized = .true.

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_def_tpw_nml", iunit)
read(iunit, nml = obs_def_tpw_nml, iostat = rc)
call check_namelist_read(iunit, rc, "obs_def_tpw_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_def_tpw_nml)
if (do_nml_term()) write(     *     , nml=obs_def_tpw_nml)

end subroutine initialize_module


!------------------------------------------------------------------------------
subroutine get_expected_tpw_distrib(state_ens_handle, location, tpw, istatus)

!------------------------------------------------------------------------------
! Purpose:  To calculate total precipitable water in a column over oceans.
! inputs:
!    state_vector:    DART state vector
!    location:        Observation location
!
! output parameters:
!    tpw:     total amount of liquid water (in cm) if all atmospheric water
!               vapor in the column was condensed.
!    istatus: 0 if ok, a positive value for error
!------------------------------------------------------------------------------
!  Author: Hui Liu ,  Version 1.1: May 25, 2011 for WRF
!  updated by n. collins  14 june 2012
!------------------------------------------------------------------------------

type(ensemble_type), intent(in)  :: state_ens_handle
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: tpw(:)
integer,             intent(out) :: istatus(:)

! local variables
real(r8) :: lon, lat, height, obsloc(3)
type(location_type) :: location2

! we'll compute the midpoint value for each pressure range, so allocate one
! more than the number of expected values.
!real(r8) :: pressure(max_pressure_intervals+1), qv(max_pressure_intervals+1)
real(r8), allocatable :: pressure(:, :), qv(:, :)
real(r8), allocatable :: pressure_interval(:), psfc(:), temp(:) ! ens size

integer  :: which_vert, k, lastk, first_non_surface_level
integer, allocatable :: track_status(:), lastk_ens(:)
logical, allocatable :: still_in(:)
integer  :: e, ens_size

if ( .not. module_initialized ) call initialize_module

ens_size = copies_in_window(state_ens_handle)

allocate(pressure(ens_size, max_pressure_intervals+1), qv(ens_size, max_pressure_intervals+1))
allocate(pressure_interval(ens_size), psfc(ens_size), lastk_ens(ens_size), track_status(ens_size))
allocate(temp(ens_size), still_in(ens_size))

! start assuming an error.  if we return early these will already
! be set.  replace the 99 error code with more specific values,
! and then 99 will only be returned in case of an uncaught error.

tpw(:) = missing_r8
istatus(:) = 99
track_status(:) = 99

! check for bad values in the namelist which exceed the size of the
! arrays that are hardcoded here.

if (num_pressure_intervals > max_pressure_intervals) then
   call error_handler(E_ERR, 'get_expected_tpw', &
                      'num_pressure_intervals greater than max allowed', &
                      source, revision, revdate, &
                      text2='increase max_pressure_intervals in obs_def_tpw_mod.f90', &
                      text3='and recompile.');
endif

! location is the lat/lon where we need to compute the column quantity
obsloc   = get_location(location)
lon      = obsloc(1)                       ! degree: 0 to 360
lat      = obsloc(2)                       ! degree: -90 to 90

! get the pressure at the surface first.

! This assumes the column is over an ocean where the surface
! is at 0m elevation.  If you are going to use this over land, the
! height below must be the elevation (in m) of the surface at this location.
which_vert = VERTISSURFACE
height = 0.0
location2 = set_location(lon, lat, height,  which_vert)

! interpolate the surface pressure and specific humidity at the desired location
! assumes the values returned from the interpolation will be in these units:
!   surface pressure :  Pa
!   moisture         :  kg/kg
call interpolate_distrib(location2, KIND_SURFACE_PRESSURE, istatus, pressure(:, 1), state_ens_handle)
do e = 1, ens_size
   if (istatus(e) /= 0) then
      tpw(e) = missing_r8
   endif
enddo
track_status = istatus
if (all(track_status /=0)) then
   return
endif


call interpolate_distrib(location2, KIND_SPECIFIC_HUMIDITY, istatus, qv(:, 1), state_ens_handle)
do e = 1, ens_size
   if (istatus(e) /= 0) then
      tpw(e) = missing_r8
      track_status(e) = istatus(e)
   endif
enddo

! save this for use below
psfc(:) = pressure(:, 1)

! there are two options for constructing the column of values.  if 'model_levels'
! is true, we query the model by vertical level number.  the 'separate_surface_level'
! flag should be set to indicate if the lowest level of the 3d grid is the
! surface or if the surface values are a separate quantity below the 3d grid.

if (model_levels) then

   ! some models have a 3d grid of values and the lowest level contains
   ! the surface quantities.  others have a separate field for the
   ! surface values and the 3d grid starts at some given elevation.
   ! if the namelist value 'separate_surface_level'  is true, we will
   ! ask to interpolate a surface pressure first and then work up the
   ! 3d column starting at level 1.  if it is false, we assume level 1
   ! was the surface pressure and we start here at level 2.

   if (separate_surface_level) then
      first_non_surface_level = 1
   else
      first_non_surface_level = 2
   endif

   ! construct a pressure column on model levels

   ! call the model until the interpolation call fails (above the top level)
   ! (this is not a fatal error unless the first call fails).
   ! also exit the loop if the pressure is above the namelist-specified pressure top

   lastk_ens(:) = 2
   still_in(:) = .true.

   lastk = 2 ! Last k could be different for each ensmeble member
   LEVELS: do k=first_non_surface_level, 10000   ! something unreasonably large

      ! call the model_mod to get the pressure and specific humidity at each level
      ! from the model and fill out the pressure and qv arrays.  the model must
      ! support a vertical type of level number.

      which_vert = VERTISLEVEL
      location2 = set_location(lon, lat, real(k, r8),  which_vert)

      call interpolate_distrib(location2, KIND_PRESSURE, istatus, temp, state_ens_handle)
      ! load up ensemble members still in
      do e = 1, ens_size
         if(still_in(e)) then
            pressure(e, lastk) = temp(e)
            if (istatus(e) /= 0) still_in(e) = .false.
            if (pressure(e, lastk) < pressure_top) still_in = .false.
         endif
      enddo

      call interpolate_distrib(location2, KIND_SPECIFIC_HUMIDITY, istatus, temp, state_ens_handle)
      do e = 1, ens_size
         if (still_in(e)) then
            qv(e, lastk) = temp(e)
            lastk_ens(e) = lastk
          endif
      enddo

      ! check to see if all ensemble members have failed
      if (all(still_in .eqv. .false.)) then
         exit LEVELS
      endif

      lastk = lastk + 1
   enddo LEVELS

   ! if we got no valid values, set istatus and return here.
   ! 'tpw' return value is already set to missing_r8
   do e = 1, ens_size
      if (lastk_ens(e) == 1) then ! ensemble member failed
         track_status(e) = 3
         !return ! can't return yet
      endif
   enddo

else

   ! construct an explicit pressure column and get qv at each pressure level.
   ! each column will be at a different set of heights because it
   ! divides the surface pressure and 'pressure_top' evenly into
   ! 'num_pressure_intervals'.  so each column will have the same
   ! number of samples, but each may be at different pressure values
   ! depending on the surface pressure.
   pressure_interval(:) = (psfc(:) - pressure_top)/num_pressure_intervals
   lastk = num_pressure_intervals + 1  ! This is the same for each ensemble member?
   lastk_ens(:) = lastk

   ! construct a pressure column at fixed pressure intervals
   ! pressure(1) is always the surface pressure.
   do k=2, lastk
      pressure(:, k) =  pressure(:, 1) - pressure_interval(:)* (k-1)
   end do

   ! call the model_mod to get the specific humidity at each location from the model
   ! and fill out the qv array.
   do e = 1, ens_size  ! The location is set by pressure so there is a loop around the enesemble
                       ! Would be better to not loop around
      do k=2, lastk

         which_vert = VERTISPRESSURE
         location2 = set_location(lon, lat, pressure(e, k),  which_vert)

         call interpolate_distrib(location2,  KIND_SPECIFIC_HUMIDITY, istatus, temp, state_ens_handle)
         if(istatus(e) == 0 ) then ! fill out the qv array for this ensemble member
            qv(e, k) = temp(e)
         endif

         if(istatus(e) /= 0 ) then ! fail
            track_status(e) = istatus(e)
         endif
      enddo

   enddo

endif

! whichever way the column was made (pressure levels or model levels),
! sum the values in the column, computing the area under the curve.
! pressure is in pascals (not hPa or mb), and moisture is in kg/kg.

tpw = 0.0
do e = 1, ens_size
   do k=1, lastk_ens(e) - 1
      if (track_status(e) == 0) then
         tpw(e) = tpw(e) + 0.5 * (qv(e, k) + qv(e, k+1) ) * (pressure(e, k) - pressure(e, k+1) )
      endif
   enddo
enddo

do e = 1, ens_size
   if (track_status(e) == 0) then
      tpw(e) = 100.0 * tpw(e) /(density*gravity)   ! -> cm
   else
      tpw(e) = missing_r8
   endif
enddo

istatus(:) = track_status(:)


deallocate(pressure, qv)
deallocate(pressure_interval, psfc, lastk_ens, track_status)
deallocate(temp, still_in)

end subroutine get_expected_tpw_distrib


end module obs_def_tpw_mod

! END DART PREPROCESS MODULE CODE

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
