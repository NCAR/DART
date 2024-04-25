! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! This is for the GRACE total water storage observation operator.

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! WATER_UNCONFINED_AQUIFER,   QTY_LIQUID_WATER,        COMMON_CODE
! GRACE_TOTAL_WATER_STORAGE,  QTY_TOTAL_WATER_STORAGE
! END DART PREPROCESS TYPE DEFINITIONS

!-------------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
! use obs_def_GRACE_mod, only: get_expected_TWS
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
! case(GRACE_TOTAL_WATER_STORAGE)
! call get_expected_TWS(state_handle, ens_size, location, &
!                       obs_def%key, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!    case(GRACE_TOTAL_WATER_STORAGE)
!       continue
! END DART PREPROCESS READ_OBS_DEF
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!    case(GRACE_TOTAL_WATER_STORAGE)
!       continue
! END DART PREPROCESS WRITE_OBS_DEF
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!    case(GRACE_TOTAL_WATER_STORAGE)
!       continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_GRACE_mod

use        types_mod, only : r4, r8, MISSING_R8, PI, DEG2RAD, RAD2DEG

use     obs_kind_mod, only : QTY_SOIL_MOISTURE, QTY_SNOW_WATER, QTY_AQUIFER_WATER

use          location_mod, only : location_type, get_close_type, get_dist, &
                                  get_close_obs, get_close_state, &
                                  convert_vertical_obs, convert_vertical_state, &
                                  set_location, set_location_missing, &
                                  query_location, write_location, &
                                  get_location, is_vertical,  &
                                  VERTISSURFACE, VERTISHEIGHT

use    utilities_mod, only : register_module, E_ERR, E_MSG, error_handler, &
                             check_namelist_read, find_namelist_in_file,   &
                             nmlfileunit, do_output, do_nml_file, do_nml_term, &
                             file_exist, is_longitude_between

use  assim_model_mod, only : interpolate

use  ensemble_manager_mod, only : ensemble_type
use obs_def_utilities_mod, only : track_status
use typesizes
use netcdf

implicit none
private

public ::  get_expected_TWS

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_GRACE_mod.f90'
character(len=*), parameter :: revision = '$Revision$'
character(len=*), parameter :: revdate  = '$Date$'

logical               :: module_initialized = .false.
logical               :: debug              = .false.
integer               :: max_grace_obs      = 100000
real(r8), allocatable :: grace_data(:)

namelist /obs_def_grace_nml/ max_grace_obs, debug

character(len=512) :: string1, string2, string3

contains

!-----------------------------------------------------------------------
!>! initialize key number and allocate space for obs data

subroutine initialize_module

integer :: iunit, rc

module_initialized = .true.
! Log the version of this source file.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file('input.nml', 'obs_def_grace_nml', iunit)
read(iunit, nml = obs_def_grace_nml, iostat = rc)
call check_namelist_read(iunit, rc, 'obs_def_grace_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_def_grace_nml)
if (do_nml_term()) write(     *     , nml=obs_def_grace_nml)

! find max number of grace obs which can be stored, and initialize type
allocate(grace_data(max_grace_obs), stat = rc)
if (rc /= 0) then
   write(string1, *) 'initial allocation failed for grace observation data,', &
                       'itemcount = ', max_grace_obs
   call error_handler(E_ERR,'initialize_module', string1, &
                      source, revision, revdate)
endif

end subroutine initialize_module


!-------------------------------------------------------------------------------
!> The forward operator for Total Water Storage.
!> Get the soil moisture at several layers, add in the snow amounts ...

subroutine get_expected_TWS(state_handle, ens_size, location, &
                            velkey, expected_obs, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
type(location_type)              :: location_soil
integer,             intent(in)  :: velkey
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

character(len=*), parameter :: routine = 'get_expected_TWS'

real(r8):: obs_h2osno(ens_size)
real(r8):: h2osoi_1(ens_size)
real(r8):: h2osoi_2(ens_size)
real(r8):: h2osoi_3(ens_size)
real(r8):: h2osoi_4(ens_size)
real(r8):: obs_h2osoi(ens_size)
real(r8):: obs_wa(ens_size)

integer :: istatus1(ens_size)
integer :: istatus2(ens_size)
integer :: istatus3(ens_size)
logical :: return_now

integer :: num
real(r8):: loc(3)

istatus      = 0
expected_obs = MISSING_R8

if ( .not. module_initialized ) call initialize_module

call error_handler(E_ERR,routine,'lots of work to do',source, revision, revdate)

!>@todo Determine number of layers and their thicknesses automatically

! Simple error check on key number before accessing the array
!call velkey_out_of_range(velkey,'get_expected_TWS')

loc = get_location(location)

!-------------------------------------------------------------------------------
! Get the soil moisture at 4 different levels

loc(3) = 0.05_r8
location_soil=set_location(loc(1),loc(2),loc(3),VERTISHEIGHT)
call interpolate(state_handle, ens_size, location_soil, QTY_SOIL_MOISTURE, &
                 h2osoi_1, istatus1)
call track_status(ens_size, istatus1, h2osoi_1, istatus, return_now)
if (return_now) return

loc(3) = 0.25_r8
location_soil=set_location(loc(1),loc(2),loc(3),VERTISHEIGHT)
call interpolate(state_handle, ens_size, location_soil, QTY_SOIL_MOISTURE, &
                 h2osoi_2, istatus1)
call track_status(ens_size, istatus1, h2osoi_2, istatus, return_now)
if (return_now) return

loc(3) = 0.7_r8
location_soil=set_location(loc(1),loc(2),loc(3),VERTISHEIGHT)
call interpolate(state_handle, ens_size, location_soil, QTY_SOIL_MOISTURE, &
                 h2osoi_3, istatus1)
call track_status(ens_size, istatus1, h2osoi_3, istatus, return_now)
if (return_now) return

loc(3) = 1.5_r8
location_soil=set_location(loc(1),loc(2),loc(3),VERTISHEIGHT)
call interpolate(state_handle, ens_size, location_soil, QTY_SOIL_MOISTURE, &
                 h2osoi_4, istatus1)
call track_status(ens_size, istatus1, h2osoi_4, istatus, return_now)
if (return_now) return

!>@todo replace the 100, 300, 600, 1000 with the actual layer thicknesses.
obs_h2osoi = h2osoi_1*100 + h2osoi_2*300 + h2osoi_3*600 + h2osoi_4*1000

if (debug) then
   write(string1,*)'soil moisture components:', &
        h2osoi_1(1),h2osoi_2(1),h2osoi_3(1),h2osoi_4(1),obs_h2osoi(1)
   call error_handler(E_MSG, 'get_expected_TWS', string1, &
                     source, revision, revdate)
endif

! Get the amount of water in the snow

call interpolate(state_handle, ens_size, location, QTY_SNOW_WATER, &
                 obs_h2osno, istatus2)
call track_status(ens_size, istatus2, obs_h2osno, istatus, return_now)
if (return_now) return

!-------------------------------------------------------------------------------
! Get the amount of water in the aquifer ...
!>@todo some models may not have aquifer ... should this be an error ...

call interpolate(state_handle, ens_size, location, QTY_AQUIFER_WATER, &
                 obs_wa, istatus3)
call track_status(ens_size, istatus3, obs_wa, istatus, return_now)
if (return_now) return

!-------------------------------------------------------------------------------
! Add up the water amounts from the different sources
! warn if things dont make sense

where (istatus == 0 ) expected_obs = obs_h2osno + obs_h2osoi + obs_wa

do num = 1, ens_size
   if (expected_obs(num) < 0.0_r8 .and. istatus(num) == 0) then

      call write_location(0,location,charstring=string3)
      write(string1,*)'TWS expected_obs for member ',num,' at ', &
                      trim(string3),' is unphysical ',expected_obs(num)
      write(string2,*)'obs_h2osno, obs_h2osoi, obs_wa', &
                     obs_h2osno(num), obs_h2osoi(num), obs_wa(num)
      write(string3,*)'istatus1, istatus2, istatus3', &
                     istatus1(num), istatus2(num), istatus3(num)  
      call error_handler(E_MSG, 'get_expected_TWS', string1, &
                     source, revision, revdate, text2=string2, text3=string3)
   
      expected_obs(num) = MISSING_R8
      istatus = 1
   
   endif
enddo

end subroutine get_expected_TWS


end module obs_def_GRACE_mod

! END DART PREPROCESS MODULE CODE
!-------------------------------------------------------------------------------
