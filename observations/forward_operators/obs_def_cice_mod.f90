! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

!>@todo FIXME: check to see if obs are of volume or thickness - for now we
! will assume volume.

! BEGIN DART PREPROCESS TYPE DEFINITIONS
!SYN_SEAICE_CONCENTR,             QTY_SEAICE_CONCENTR,           COMMON_CODE
!SAT_U_SEAICE_COMPONENT,          QTY_U_SEAICE_COMPONENT,        COMMON_CODE
!SAT_V_SEAICE_COMPONENT,          QTY_V_SEAICE_COMPONENT,        COMMON_CODE
!SAT_SEAICE_CONCENTR,             QTY_SEAICE_CONCENTR,           COMMON_CODE
!SAT_SEAICE_VOLUME,               QTY_SEAICE_VOLUME,             COMMON_CODE
!SAT_SEAICE_SNOWVOLUME,           QTY_SEAICE_SNOWVOLUME,         COMMON_CODE
!SAT_SEAICE_SURFACETEMP,          QTY_SEAICE_SURFACETEMP,        COMMON_CODE
!SAT_SEAICE_FY,                   QTY_SEAICE_FY,                 COMMON_CODE
!SAT_SEAICE_AGREG_FY,             QTY_SEAICE_AGREG_FY
!SAT_SEAICE_AGREG_SURFACETEMP,    QTY_SEAICE_AGREG_SURFACETEMP
!SAT_SEAICE_AGREG_FREEBOARD,      QTY_SEAICE_AGREG_FREEBOARD
!SAT_SEAICE_AGREG_CONCENTR,       QTY_SEAICE_AGREG_CONCENTR
!SAT_SEAICE_AGREG_VOLUME,         QTY_SEAICE_AGREG_VOLUME
!SAT_SEAICE_AGREG_SNOWVOLUME,     QTY_SEAICE_AGREG_SNOWVOLUME
!SAT_SEAICE_AGREG_THICKNESS,      QTY_SEAICE_AGREG_THICKNESS
!SAT_SEAICE_AGREG_SNOWDEPTH,      QTY_SEAICE_AGREG_SNOWDEPTH
! END DART PREPROCESS TYPE DEFINITIONS

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_cice_mod, only : get_expected_agreg_freeboard, &
!                               get_expected_agreg_over_grid, &
!                               get_expected_agreg_over_ice,  &
!                               get_expected_agreg_thickness
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!   case(SAT_SEAICE_AGREG_FREEBOARD)
!      call get_expected_agreg_freeboard(state_handle, ens_size, location, &
!               QTY_SEAICE_CONCENTR, QTY_SEAICE_VOLUME, QTY_SEAICE_SNOWVOLUME, &
!               expected_obs, istatus)
!   case(SAT_SEAICE_AGREG_THICKNESS)
!      call get_expected_agreg_thickness(state_handle, ens_size, location, &
!               QTY_SEAICE_VOLUME, expected_obs, istatus)
!   case(SAT_SEAICE_AGREG_SNOWDEPTH)
!      call get_expected_agreg_thickness(state_handle, ens_size, location, &
!               QTY_SEAICE_SNOWVOLUME, expected_obs, istatus)
!   case(SAT_SEAICE_AGREG_CONCENTR)
!      call get_expected_agreg_over_grid(state_handle, ens_size, location, &
!               QTY_SEAICE_CONCENTR, expected_obs,istatus)
!   case(SAT_SEAICE_AGREG_VOLUME)
!      call get_expected_agreg_over_grid(state_handle, ens_size, location, &
!               QTY_SEAICE_VOLUME, expected_obs, istatus)
!   case(SAT_SEAICE_AGREG_SNOWVOLUME)
!      call get_expected_agreg_over_grid(state_handle, ens_size, location, &
!               QTY_SEAICE_SNOWVOLUME, expected_obs, istatus)
!   case(SAT_SEAICE_AGREG_SURFACETEMP)
!      call get_expected_agreg_over_ice(state_handle, ens_size, location, &
!               QTY_SEAICE_SURFACETEMP, expected_obs, istatus)
!   case(SAT_SEAICE_AGREG_FY)
!      call get_expected_agreg_over_ice(state_handle, ens_size, location, &
!               QTY_SEAICE_FY, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!    case(SAT_SEAICE_AGREG_FREEBOARD,   &
!         SAT_SEAICE_AGREG_THICKNESS,   &
!         SAT_SEAICE_AGREG_SNOWDEPTH,   &
!         SAT_SEAICE_AGREG_CONCENTR,    &
!         SAT_SEAICE_AGREG_VOLUME,      &
!         SAT_SEAICE_AGREG_SNOWVOLUME,  &
!         SAT_SEAICE_AGREG_SURFACETEMP, &
!         SAT_SEAICE_AGREG_FY)
!       continue
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!    case(SAT_SEAICE_AGREG_FREEBOARD,   &
!         SAT_SEAICE_AGREG_THICKNESS,   &
!         SAT_SEAICE_AGREG_SNOWDEPTH,   &
!         SAT_SEAICE_AGREG_CONCENTR,    &
!         SAT_SEAICE_AGREG_VOLUME,      &
!         SAT_SEAICE_AGREG_SNOWVOLUME,  &
!         SAT_SEAICE_AGREG_SURFACETEMP, &
!         SAT_SEAICE_AGREG_FY)
!       continue
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!    case(SAT_SEAICE_AGREG_FREEBOARD,   &
!         SAT_SEAICE_AGREG_THICKNESS,   &
!         SAT_SEAICE_AGREG_SNOWDEPTH,   &
!         SAT_SEAICE_AGREG_CONCENTR,    &
!         SAT_SEAICE_AGREG_VOLUME,      &
!         SAT_SEAICE_AGREG_SNOWVOLUME,  &
!         SAT_SEAICE_AGREG_SURFACETEMP, &
!         SAT_SEAICE_AGREG_FY)
!       continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_cice_mod

use        types_mod, only : r8, missing_r8, PI, deg2rad

use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG,    &
                             check_namelist_read, find_namelist_in_file,      &
                             nmlfileunit, do_output, do_nml_file,do_nml_term, &
                             ascii_file_format

use     location_mod, only : location_type, set_location, get_location, &
                             VERTISLEVEL

use  assim_model_mod, only : interpolate

use     obs_kind_mod, only : QTY_SEAICE_VOLUME,      &
                             QTY_SEAICE_CONCENTR,    &
                             QTY_SEAICE_SNOWVOLUME,  &
                             QTY_SEAICE_FY,          &
                             QTY_SEAICE_SURFACETEMP, &
                             QTY_SEAICE_CATEGORY

use  ensemble_manager_mod, only : ensemble_type

use obs_def_utilities_mod, only : track_status

implicit none

public :: get_expected_agreg_freeboard, &
          get_expected_agreg_thickness, &
          get_expected_agreg_over_ice,  &
          get_expected_agreg_over_grid

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.
integer       :: Ncat = 0

character(len=512) :: string1, string2, string3

!>@todo FIXME ... instead of using multiple istatus? arrays, could use
!>                the track_status() routine

contains

!-----------------------------------------------------------------------------
!> The number of ice categories is needed for all the forward operators.
!> The number of categories will be determined ONCE.

subroutine initialize_module(state_handle, ens_size)
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size

if (module_initialized) return
call register_module(source, revision, revdate)
module_initialized = .true.

Ncat = get_number_of_ice_categories(state_handle, ens_size)

end subroutine initialize_module


!-----------------------------------------------------------------------------
!> 

subroutine get_expected_agreg_freeboard(state_handle, ens_size, location,  &
                                       var_sic, var_siv, var_snv, &
                                       agreg_fb, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: var_sic, var_siv, var_snv
real(r8),            intent(out) :: agreg_fb(ens_size)
integer,             intent(out) :: istatus(ens_size)

!fb_volume is over grid, agreg_fb is over sea ice
!for real observations we have fb
real(r8) :: fb_volume(ens_size)
real(r8) :: ice_volume(ens_size)
real(r8) :: snow_volume(ens_size)
real(r8) :: agreg_sic(ens_size)

real(r8) :: loc_array(3),llat,llon
integer  :: icat
integer  :: istatus1(ens_size)
integer  :: istatus2(ens_size)
integer  :: istatus3(ens_size)

type(location_type) :: location_fake

real(r8), parameter :: ice_dens   =  917.0_r8, &
                       snow_dens  =  330.0_r8, &
                       water_dens = 1026.0_r8

if (.not.module_initialized) call initialize_module(state_handle, ens_size)

loc_array = get_location(location)
llat = loc_array(1)
llon = loc_array(2)

istatus   = 0
fb_volume = 0.0_r8
agreg_fb  = 0.0_r8

call get_expected_agreg_over_grid(state_handle, ens_size, location, &
         QTY_SEAICE_CONCENTR, agreg_sic, istatus3)

! model_interpolate interpolates from the model grid to a specific location
! the 3rd variable location%z contains cat_index. The observation itself may 
! not contain cat_index, so we need create one.

do icat = 1, Ncat

   location_fake = set_location(llat,llon,real(icat,r8),VERTISLEVEL)

   call interpolate(state_handle, ens_size, location_fake, &
                    QTY_SEAICE_VOLUME, ice_volume, istatus1)

   call interpolate(state_handle, ens_size, location_fake, &
                    QTY_SEAICE_SNOWVOLUME, snow_volume, istatus2)

   where(istatus1==0 .and. istatus2==0) &
     fb_volume = fb_volume + ice_volume*(1 - ice_dens/water_dens) - &
                             snow_volume*snow_dens/water_dens
end do

where(istatus3==0 .and. agreg_sic>1e-6) agreg_fb = fb_volume/agreg_sic
where(istatus1/=0) istatus = istatus1
where(istatus2/=0) istatus = istatus2
where(istatus3/=0) istatus = istatus3

end subroutine get_expected_agreg_freeboard


!-----------------------------------------------------------------------------
!> The forward operator is to sum the model state over all the categories, 
!> unit per grid area. Nothing else is applied.
!> Good for aggregate sea ice volume and snow volume

subroutine get_expected_agreg_over_grid(state_handle, ens_size, location, &
               obstype, agreg_grid, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obstype
real(r8),            intent(out) :: agreg_grid(ens_size)
integer,             intent(out) :: istatus(ens_size)

type(location_type) :: location_fake

real(r8)  :: grid(ens_size)
integer   :: icat
real(r8)  :: loc_array(3),llat,llon

if (.not.module_initialized) call initialize_module(state_handle, ens_size)

istatus    = 0
agreg_grid = 0.0_r8

loc_array = get_location(location)
llat = loc_array(1)
llon = loc_array(2)

! model_interpolate interpolates from the model grid to a specific location
! the 3rd variable location%z contains cat_index. The observation itself may 
! not contain cat_index, so we need create one.

do icat = 1, Ncat
   location_fake = set_location(llat,llon,real(icat,r8),VERTISLEVEL)
   call interpolate(state_handle, ens_size, location_fake, obstype, grid ,istatus)
   where(istatus == 0) agreg_grid = agreg_grid + grid
end do

end subroutine get_expected_agreg_over_grid


!-----------------------------------------------------------------------------
!> This is to calculate the agregate value over ice only
!> Unit per ice area

subroutine get_expected_agreg_over_ice(state_handle, ens_size, location, &
               obstype, agreg_ice, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obstype
real(r8),            intent(out) :: agreg_ice(ens_size)
integer,             intent(out) :: istatus(ens_size)

real(r8)  :: sic(ens_size)
real(r8)  :: ice(ens_size)
real(r8)  :: agreg_sic(ens_size)
real(r8)  :: temp_ice(ens_size)
integer   :: istatus1(ens_size)
integer   :: istatus2(ens_size)
integer   :: istatus3(ens_size)
integer   :: icat

real(r8)  :: loc_array(3),llat,llon
type(location_type) :: location_fake

if (.not.module_initialized) call initialize_module(state_handle, ens_size)

loc_array = get_location(location)
llat = loc_array(1)
llon = loc_array(2)
istatus = 0

temp_ice = 0

call get_expected_agreg_over_grid(state_handle, ens_size, location, &
         QTY_SEAICE_CONCENTR, agreg_sic, istatus3)

! model_interpolate interpolates from the model grid to a specific location
! the 3rd variable location%z contains cat_index. The observation itself may 
! not contain cat_index, so we need create one.

do icat =1, Ncat
   location_fake = set_location(llat,llon,real(icat,r8),VERTISLEVEL)
   call interpolate(state_handle, ens_size, location_fake, &
                    QTY_SEAICE_CONCENTR, sic, istatus1)
   call interpolate(state_handle, ens_size, location_fake, obstype, ice, istatus2)

   where(istatus1 /= 0) istatus = istatus1
   where(istatus2 /= 0) istatus = istatus2
   where(istatus  == 0) temp_ice = temp_ice + sic * ice
end do

!where(agreg_sic>1e-6) agreg_ice = temp_ice/agreg_sic
agreg_ice = temp_ice/max(agreg_sic,1e-8)

end subroutine get_expected_agreg_over_ice


!-----------------------------------------------------------------------------
!> For thickness, it's sum(volume*concentr)
!> Good for sea ice thickness and snow thickness

subroutine get_expected_agreg_thickness(state_handle, ens_size, location, &
               obstype, agreg_thickness, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obstype
real(r8),            intent(out) :: agreg_thickness(ens_size)
integer,             intent(out) :: istatus(ens_size)

real(r8)  :: agreg_sic(ens_size)
real(r8)  :: agreg_volume(ens_size)
integer   :: istatus1(ens_size), istatus2(ens_size)
real(r8)  :: loc_array(3),llat,llon

if (.not.module_initialized) call initialize_module(state_handle, ens_size)

loc_array = get_location(location)
llat = loc_array(1)
llon = loc_array(2)
istatus = 0

agreg_thickness = 0

call get_expected_agreg_over_grid(state_handle, ens_size, location, &
         QTY_SEAICE_CONCENTR, agreg_sic, istatus1)

call get_expected_agreg_over_grid(state_handle, ens_size, location,&
         obstype, agreg_volume, istatus2)

where(istatus1 /= 0) istatus = istatus1
where(istatus2 /= 0) istatus = istatus2
!where(istatus == 0.and.agreg_sic>1e-6) agreg_thickness = agreg_volume/agreg_sic
where(istatus  == 0) agreg_thickness = agreg_volume/max(agreg_sic,1e-8)

end subroutine get_expected_agreg_thickness


!-----------------------------------------------------------------------------
!> This determines the number of sea ice categories.

function get_number_of_ice_categories(state_handle, ens_size)

integer :: get_number_of_ice_categories
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size

integer, parameter :: MAX_CATEGORIES = 100

integer  :: istatus(ens_size)
real(r8) :: category(ens_size)
real(r8) :: llat = 0.0_r8
real(r8) :: llon = 0.0_r8
integer  :: icat
type(location_type) :: location_fake

get_number_of_ice_categories = 0

CATEGORIES : do icat =1, MAX_CATEGORIES

   location_fake = set_location(llat, llon, real(icat,r8), VERTISLEVEL)

   call interpolate(state_handle, ens_size, location_fake, &
                    QTY_SEAICE_CATEGORY, category, istatus)

   if (any(istatus /= 0)) exit CATEGORIES

   get_number_of_ice_categories = get_number_of_ice_categories + 1

enddo CATEGORIES

if (get_number_of_ice_categories == MAX_CATEGORIES) then
   write(string1,*)'at capacity of ice categories.'
   write(string2,*)'If you have more than ',MAX_CATEGORIES,' ice categories'
   write(string3,*)'modify "MAX_CATEGORIES", recompile and try again.'
   call error_handler(E_ERR,'get_number_of_ice_categories',string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

if (get_number_of_ice_categories == 0) then
   write(string1,*)'Could not determine the number of ice categories.'
   call error_handler(E_ERR,'get_number_of_ice_categories',string1, &
              source, revision, revdate)
endif

end function get_number_of_ice_categories


end module obs_def_cice_mod

! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------

