! Copyright 2019 University Corporation for Atmospheric Research and 
! Colorado Department of Public Health and Environment.
!
! Licensed under the Apache License, Version 2.0 (the "License"); you may not use 
! this file except in compliance with the License. You may obtain a copy of the 
! License at      http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
! CONDITIONS OF ANY KIND, either express or implied. See the License for the 
! specific language governing permissions and limitations under the License.
!
! Development of this code utilized the RMACC Summit supercomputer, which is 
! supported by the National Science Foundation (awards ACI-1532235 and ACI-1532236),
! the University of Colorado Boulder, and Colorado State University.
! The Summit supercomputer is a joint effort of the University of Colorado Boulder
! and Colorado State University.

! BEGIN DART PREPROCESS KIND LIST
! MEXICO_AQS_PM10,         QTY_PM10
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_MEXICO_AQS_PM10_mod, only : get_expected_mexico_aqs_pm10
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!   case(MEXICO_AQS_PM10)
!        call get_expected_mexico_aqs_pm10(state_handle, ens_size, location, obs_def%key, expected_obs, istatus) 
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(MEXICO_AQS_PM10)
!     continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(MEXICO_AQS_PM10)
!     continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(MEXICO_AQS_PM10)
!     continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_MEXICO_AQS_PM10_mod

   use types_mod, only : r8, missing_r8
   use utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             nmlfileunit, check_namelist_read, &
                             find_namelist_in_file, do_nml_file, do_nml_term, &
                             ascii_file_format
   use location_mod, only : location_type, set_location, get_location, VERTISPRESSURE, VERTISLEVEL, VERTISSURFACE, VERTISUNDEF

   use assim_model_mod, only : interpolate
   use obs_kind_mod
   use ensemble_manager_mod,  only : ensemble_type
   use obs_def_utilities_mod, only : track_status

   implicit none
   private

   public :: get_expected_mexico_aqs_pm10

! version controlled file description for error handling, do not edit
   character(len=*), parameter :: source   = 'obs_def_MEXICO_AQS_PM10_mod.f90'
   character(len=*), parameter :: revision = ''
   character(len=*), parameter :: revdate  = ''

   character(len=512) :: string1, string2

   logical, save :: module_initialized = .false.
   logical       :: use_log_pm10
   namelist /obs_def_MEXICO_AQS_PM10_nml/ use_log_pm10

contains

! ---------------------------------------------------

subroutine initialize_module

   integer :: iunit, rc

! Prevent multiple calls frm executing this code more than once.
   if (module_initialized) return

   call register_module(source, revision, revdate)
   module_initialized = .true.

! Read the namelist entry.
   use_log_pm10=.false.
   call find_namelist_in_file("input.nml", "obs_def_MEXICO_AQS_PM10_nml", iunit)
   read(iunit, nml = obs_def_MEXICO_AQS_PM10_nml, iostat = rc)
   call check_namelist_read(iunit, rc, "obs_def_MEXICO_AQS_PM10_nml")

! Record the namelist values used for the run ... 
   if (do_nml_file()) write(nmlfileunit, nml=obs_def_MEXICO_AQS_PM10_nml)
   if (do_nml_term()) write(     *     , nml=obs_def_MEXICO_AQS_PM10_nml)
end subroutine initialize_module

!------------------------------------------------------------------------

subroutine get_expected_mexico_aqs_pm10(state_handle, ens_size, location, key, val, istatus)  
   type(ensemble_type), intent(in)  :: state_handle
   integer,             intent(in)  :: ens_size
   type(location_type), intent(in)  :: location
   integer,             intent(in)  :: key
   real(r8),            intent(out) :: val(:) ! (ug/m^3 dry air)
   integer,             intent(out) :: istatus(:)
   integer                          :: zstatus(ens_size)
   logical                          :: return_now
   character(len=*), parameter :: routine = 'get_expected_mexico_aqs_pm10'
!
! Forward operator for PM10.  The argument list to this routine
! must match the call in the GET_EXPECTED_OBS_FROM_DEF section above.

   type(location_type)  :: nloc
   real(r8)             :: mloc(3)
   real(r8)             :: level
   real(r8) :: prs(ens_size)
   real(r8) :: tmp(ens_size)
   real(r8) :: rho_d(ens_size)
!
   real(r8) :: p25(ens_size)      ! p25 
   real(r8) :: so4(ens_size)      ! Sulfate 
   real(r8) :: bc1(ens_size)      ! Hydrophobic Black Carbon 
   real(r8) :: bc2(ens_size)      ! Hydrophilic Black Carbon 
   real(r8) :: oc1(ens_size)      ! Hydrophobic Organic Carbon
   real(r8) :: oc2(ens_size)      ! Hydrophilic Organic Carbon 
   real(r8) :: dust1(ens_size)    ! Dust 1
   real(r8) :: dust2(ens_size)    ! Dust 2
   real(r8) :: dust3(ens_size)    ! Dust 3
   real(r8) :: dust4(ens_size)    ! Dust 4
   real(r8) :: dust5(ens_size)    ! Dust 5
   real(r8) :: ss1(ens_size)      ! Sea Salt 1
   real(r8) :: ss2(ens_size)      ! Sea Salt 2
   real(r8) :: ss3(ens_size)      ! Sea Salt 3
   real(r8) :: ss4(ens_size)      ! Sea Salt 4
   real(r8) :: eps
   real(r8) :: grav
   real(r8) :: Rd                 ! gas constant dry air
   real(r8) :: Ru                 ! universal dry air
   real(r8) :: Cp                 ! heat capacity dry air
   real(r8) :: Pa_to_torr         ! convert Pa to torr
!
! Initialize DART
   if ( .not. module_initialized ) call initialize_module

   eps = 0.61_r8         ! m/s^2
   grav = 9.8_r8         ! m/s^2
   Rd = 287.05_r8        ! J/kg/K
   Ru = 8.316_r8         ! kg/mole
   Cp = 1006.0_r8        ! J/kg/K
   Pa_to_torr = 133.322
   val(:) = 0.
!
   istatus(:)=0.
   return_now=.false.
!   
   mloc = get_location(location)
   if (mloc(2)>90.0_r8) then
      mloc(2)=90.0_r8
   elseif (mloc(2)<-90.0_r8) then
      mloc(2)=-90.0_r8
   endif
!
   level=1.0_r8   
   nloc = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
!
! pressure (Pa) at this location (midpoint) - this calls the model_mod code.
   zstatus(:)=0
   call interpolate(state_handle, ens_size, nloc, QTY_PRESSURE, prs, zstatus)
!
! temperature (K) at this location - this calls the model_mod code.
   zstatus(:)=0
   call interpolate(state_handle, ens_size, nloc, QTY_TEMPERATURE, tmp, zstatus)
!
! sulfate at this location - this calls the model_mod code.
   zstatus(:)=0.
   call interpolate(state_handle, ens_size, nloc, QTY_SO4, so4, zstatus)
   so4(:) = so4(:) * 1.e3
!
! dust at this location - this calls the model_mod code.
   zstatus(:)=0.
   call interpolate(state_handle, ens_size, nloc, QTY_DST01, dust1, zstatus)
!
! dust at this location - this calls the model_mod code.
   zstatus(:)=0.
   call interpolate(state_handle, ens_size, nloc, QTY_DST02, dust2, zstatus)
!
! dust at this location - this calls the model_mod code.
   zstatus(:)=0.
   call interpolate(state_handle, ens_size, nloc, QTY_DST03, dust3, zstatus)
!
! dust at this location - this calls the model_mod code.
   zstatus(:)=0.
   call interpolate(state_handle, ens_size, nloc, QTY_DST04, dust4, zstatus)
!
! dust at this location - this calls the model_mod code.
   zstatus(:)=0.
   call interpolate(state_handle, ens_size, nloc, QTY_DST05, dust5, zstatus)
!
! hydrophilic black carbon at this location - this calls the model_mod code.
   zstatus(:)=0.
   call interpolate(state_handle, ens_size, nloc, QTY_BC1, bc1, zstatus)
!
! hydrophilic black carbon at this location - this calls the model_mod code.
   zstatus(:)=0.
   call interpolate(state_handle, ens_size, nloc, QTY_BC2, bc2, zstatus)
!
! hydrophilic organic carbon at this location - this calls the model_mod code.
   zstatus(:)=0.
   call interpolate(state_handle, ens_size, nloc, QTY_OC1, oc1, zstatus)
!
! hydrophilic organic carbon at this location - this calls the model_mod code.
   zstatus(:)=0.
   call interpolate(state_handle, ens_size, nloc, QTY_OC2, oc2, zstatus)
!
! sea salt at this location - this calls the model_mod code.
   zstatus(:)=0.
   call interpolate(state_handle, ens_size, nloc, QTY_SSLT01, ss1, zstatus)
!
! sea salt at this location - this calls the model_mod code.
   zstatus(:)=0.
   call interpolate(state_handle, ens_size, nloc, QTY_SSLT02, ss2, zstatus)
!
! sea salt at this location - this calls the model_mod code.
   zstatus(:)=0.
   call interpolate(state_handle, ens_size, nloc, QTY_SSLT03, ss3, zstatus)
!
! sea salt at this location - this calls the model_mod code.
   zstatus(:)=0.
   call interpolate(state_handle, ens_size, nloc, QTY_SSLT04, ss4, zstatus)
!
! p25 at this location - this calls the model_mod code.
   zstatus(:)=0.
   call interpolate(state_handle, ens_size, nloc, QTY_P25, p25, zstatus)
!
! Check for missing values
   if(any(p25.lt.0.) .or. any(so4.lt.0.) .or. any(bc1.lt.0) .or. any(bc2.lt.0) .or. &
   any(oc1.lt.0) .or. any(oc2.lt.0) .or. any(dust1.lt.0) .or. any(dust2.lt.0) .or. &
   any(dust3.lt.0) .or. any(dust4.lt.0) .or. any(ss1.lt.0) .or. any(ss2.lt.0) .or. &
   any(ss3.lt.0)) then   
      zstatus(:)=20
      val(:)=missing_r8
      write(string1, *) 'APM: PM10 FO has missing values '
      call error_handler(E_MSG, routine, string1, source)
      call track_status(ens_size, zstatus, val, istatus, return_now)
      return
   endif
!
! The actual forward operator computation.  This is the value that
! will be returned.  istatus (the return code) of 0 is good,
! return any value > 0 for error.  (values < 0 reserved for
! system use.)
!
! Expected pm10 (ug/kg)
   if (use_log_pm10) then
      val(:) = (exp(p25(:)) + 1.375*exp(so4(:)) + exp(bc1(:)) + exp(bc2(:)) + &
      1.8*(exp(oc1(:))+exp(oc2(:))) + exp(dust1(:)) + exp(dust2(:)) + &
      exp(dust3(:)) + .87*exp(dust4(:)) + exp(ss1(:)) + exp(ss2(:)) +exp(ss3(:)))
   else
      val(:) = (p25(:) + 1.375*so4(:) + bc1(:) + bc2(:) + 1.8*(oc1(:)+oc2(:)) + &
      dust1(:) + dust2(:) + dust3(:) + .87*dust4(:) + ss1(:) + ss2(:) +ss3(:))
   endif
!
! Expected pm10 (ug/m^3) (LC)
   rho_d(:)=prs(:)/(Rd*tmp(:)) 
   val(:) = val(:)*rho_d
!
! Expected pm10 (STP) Should be this
   val(:) = val(:) * tmp(:)/298.15 * 760./prs(:)*Pa_to_torr
!
! Check expected observation
   if(any(val.lt.0)) then
      zstatus(:)=20
      val(:)=missing_r8
      write(string1, *) 'APM: PM10 expected observation is negative '
      call error_handler(E_MSG, routine, string1, source)
      call track_status(ens_size, zstatus, val, istatus, return_now)
      return
   endif
!  
end subroutine get_expected_mexico_aqs_pm10

!------------------------------------------------------------------------

end module obs_def_MEXICO_AQS_PM10_mod
! END DART PREPROCESS MODULE CODE

