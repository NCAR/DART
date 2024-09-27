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

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! MODIS_AOD_TOTAL_COL, QTY_AOD
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_MODIS_AOD_TOTAL_COL_mod, only : get_expected_modis_aod_total_col
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(MODIS_AOD_TOTAL_COL)                                                           
!            call get_expected_modis_aod_total_col(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(MODIS_AOD_TOTAL_COL)
!         continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(MODIS_AOD_TOTAL_COL)
!         continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(MODIS_AOD_TOTAL_COL)
!         continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_MODIS_AOD_TOTAL_COL_mod

   use types_mod, only : r8, missing_r8
   use utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             nmlfileunit, check_namelist_read, &
                             find_namelist_in_file, do_nml_file, do_nml_term, &
                             ascii_file_format, &
                             read_int_scalar, &
                             write_int_scalar, &       
                             read_r8_scalar, &
                             write_r8_scalar, &
                             read_r8_array, &
                             write_r8_array
                             
   use location_mod, only : location_type, set_location, get_location, VERTISPRESSURE, VERTISLEVEL, VERTISSURFACE, VERTISUNDEF

   use assim_model_mod, only : interpolate
   use obs_kind_mod
   use ensemble_manager_mod,  only : ensemble_type
   use obs_def_utilities_mod, only : track_status

   implicit none
   private

   public :: get_expected_modis_aod_total_col

! version controlled file description for error handling, do not edit
   character(len=*), parameter :: source   = 'obs_def_MODIS_AOD_TOTAL_COL_mod.f90'
   character(len=*), parameter :: revision = ''
   character(len=*), parameter :: revdate  = ''

   character(len=512) :: string1, string2
   character(len=200) :: upper_data_file
   character(len=200) :: upper_data_model
   character(len=200) :: model
   integer            :: ls_chem_dx, ls_chem_dy, ls_chem_dz, ls_chem_dt

   logical, save :: module_initialized = .false.
   logical       :: use_log_aod
   real          :: ptop
   integer       :: mdl_nlay
   namelist /obs_def_MODIS_AOD_nml/ upper_data_file, use_log_aod, ptop, mdl_nlay, &
   ls_chem_dx, ls_chem_dy, ls_chem_dz, ls_chem_dt, upper_data_model

contains

!----------------------------------------------------------------------

subroutine initialize_module

   integer :: iunit, rc

! Prevent multiple calls from executing this code more than once.
   if (module_initialized) return

   call register_module(source, revision, revdate)
   module_initialized = .true.

! Read the namelist entry.
   use_log_aod=.false.
   call find_namelist_in_file("input.nml", "obs_def_MODIS_AOD_nml", iunit)
   read(iunit, nml = obs_def_MODIS_AOD_nml, iostat = rc)
   call check_namelist_read(iunit, rc, "obs_def_MODIS_AOD_nml")

! Record the namelist values used for the run ... 
   if (do_nml_file()) write(nmlfileunit, nml=obs_def_MODIS_AOD_nml)
   if (do_nml_term()) write(     *     , nml=obs_def_MODIS_AOD_nml)
end subroutine initialize_module
!
subroutine get_expected_modis_aod_total_col(state_handle, ens_size, location, key, val, istatus)
!----------------------------------------------------------------------
!subroutine get_expected_modis_aod_total_col(state, ens_size, location, key,  val, istatus)
   type(ensemble_type), intent(in)  :: state_handle
   integer,             intent(in)  :: ens_size
   type(location_type), intent(in)  :: location
   integer,             intent(in)  :: key
   real(r8),            intent(out) :: val(:)
   integer,             intent(out) :: istatus(:)
   integer                          :: zstatus(ens_size)
   logical                          :: return_now
   character(len=*), parameter :: routine = 'get_expected_modis_aod_total_col'

! Forward operator for MODIS AOD.  The argument list to this routine
! must match the call in the GET_EXPECTED_OBS_FROM_DEF section above.

   type(location_type)  :: nloc,nploc
   integer  :: ilay
   real(r8) :: fac(ens_size)
   real(r8) :: qmr(ens_size)         ! Water Vapor Mixing Ratio
   real(r8) :: tmp(ens_size)         ! Temperature
   real(r8) :: tmp_vir(ens_size)     ! Virtual Temperature
   real(r8) :: psfc(ens_size)        ! Surface Pressure
   real(r8) :: prs(ens_size)         ! Pressure (current layer)
   real(r8) :: prs_p(ens_size)       ! Pressure (next layer)
   real(r8) :: ln_prs_up(ens_size)   ! Log Pressure (lower level)
   real(r8) :: ln_prs_dw(ens_size)   ! Log Pressure (upper level)
   real(r8) :: prs_up(ens_size)      ! Pressure (lower level)
   real(r8) :: prs_dw(ens_size)      ! Pressure (upper level)
   real(r8) :: thick(ens_size)       ! Layer Thickness 
!
   real(r8) :: so4(ens_size)         ! Sulfate
   real(r8) :: bc1(ens_size)         ! Hydrophobic Black Carbon 
   real(r8) :: bc2(ens_size)         ! Hydrophilic Black Carbon 
   real(r8) :: oc1(ens_size)         ! Hydrophobic Organic Carbon
   real(r8) :: oc2(ens_size)         ! Hydrophilic Organic Carbon 
   real(r8) :: dust1(ens_size)       ! Dust 1
   real(r8) :: dust2(ens_size)       ! Dust 2
   real(r8) :: dust3(ens_size)       ! Dust 3
   real(r8) :: dust4(ens_size)       ! Dust 4
   real(r8) :: dust5(ens_size)       ! Dust 5
   real(r8) :: ss1(ens_size)         ! Sea Salt 1
   real(r8) :: ss2(ens_size)         ! Sea Salt 2
   real(r8) :: ss3(ens_size)         ! Sea Salt 3
   real(r8) :: ss4(ens_size)         ! Sea Salt 4
   real(r8) :: eps
   real(r8) :: Rd                    ! gas constant dry air
   real(r8) :: Ru                    ! universal gas constant
   real(r8) :: mass_p_mole_dryair
   real(r8) :: Cp                    ! heat capacity dry air
   real(r8) :: Pa_to_torr            ! convert Pa to torr
   real(r8) :: grav                  ! gravity
   real(r8) :: ext_so4               ! extinction coefficient
   real(r8) :: ext_oc1               ! extinction coefficient
   real(r8) :: ext_oc2               ! extinction coefficient
   real(r8) :: ext_bc1               ! extinction coefficient
   real(r8) :: ext_bc2               ! extinction coefficient
   real(r8) :: ext_dst1              ! extinction coefficient
   real(r8) :: ext_dst2              ! extinction coefficient
   real(r8) :: ext_dst3              ! extinction coefficient
   real(r8) :: ext_dst4              ! extinction coefficient
   real(r8) :: ext_dst5              ! extinction coefficient
   real(r8) :: ext_ssl1              ! extinction coefficient
   real(r8) :: ext_ssl2              ! extinction coefficient
   real(r8) :: ext_ssl3              ! extinction coefficient
   real(r8) :: ext_ssl4              ! extinction coefficient
   real(r8) :: mloc(3)
   real(r8) :: level
!
! Initialize DART
   if ( .not. module_initialized ) call initialize_module
!
   eps = 0.61_r8         ! m/s^2
   grav = 9.8_r8         ! m/s^2
   Rd = 287.05_r8        ! J/kg/K
   Ru = 8.316_r8         ! kg/mole
   Cp = 1006.0_r8        ! J/kg/K
   Pa_to_torr = 133.322
   mass_p_mole_dryair=.02897  ! kg/mole 
   val(:) = 0.           ! dimensionless
   ptop = 100.           ! MODIS goes to top of atmosphere (not top of model)
!
   istatus(:)=0
   return_now = .false.
!
! From Liu et al. (2011)   
   ext_so4 = 3.13           ! (m^2)/g
   ext_oc1 = 2.65           ! (m^2)/g
   ext_oc2 = 2.65           ! (m^2)/g
   ext_bc1 = 9.16           ! (m^2)/g
   ext_bc2 = 9.16           ! (m^2)/g
   ext_ssl1 = 2.59          ! (m^2)/g
   ext_ssl2 = 0.90          ! (m^2)/g
   ext_ssl3 = 0.24          ! (m^2)/g
   ext_ssl4 = 0.097         ! (m^2)/g
   ext_dst1 = 1.61          ! (m^2)/g
   ext_dst2 = 0.51          ! (m^2)/g
   ext_dst3 = 0.27          ! (m^2)/g
   ext_dst4 = 0.14          ! (m^2)/g
   ext_dst5 = 0.076         ! (m^2)/g
!
   if(use_log_aod) then
      write(string1, *) 'APM: MODIS AOD TOTAL COL Forward Operator not set up for log_aod = T '
      call error_handler(E_MSG, 'obs_def_MODIS_AOD_TOTAL_COL', string1, source)
      call exit_all(-77)
   endif
!   
   mloc = get_location(location)
   if (mloc(2)>90.0_r8) then
      mloc(2)=90.0_r8
   elseif (mloc(2)<-90.0_r8) then
      mloc(2)=-90.0_r8
   endif
!
! surface pressure (Pa) at this location - this calls the model_mod code.
   level=0.0_r8   
   nloc = set_location(mloc(1), mloc(2), level, VERTISSURFACE)
   zstatus(:)=0
   call interpolate(state_handle, ens_size, nloc, QTY_SURFACE_PRESSURE, psfc, zstatus)
   ln_prs_dw(:)=log(psfc(:))
!
   do ilay=1,mdl_nlay-1
      level=real(ilay)
      nloc = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
!
      if(ilay.lt.mdl_nlay-1) then
         level=real(ilay+1)
         nploc = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
!
! pressure (Pa) at this location (midpoint + 1) - this calls the model_mod code.
         zstatus(:)=0
         call interpolate(state_handle, ens_size, nploc, QTY_PRESSURE, prs_p, zstatus)
      elseif(ilay.eq.mdl_nlay-1) then
         prs_p(:)=ptop
      endif
!
! pressure (Pa) at this location (midpoint) - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_PRESSURE, prs, zstatus)
!
! temperature (K) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_TEMPERATURE, tmp, zstatus)
!
! water vapor mixing ratio (kg/kg) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_VAPOR_MIXING_RATIO, qmr, zstatus)
!
! thickness (m)
      ln_prs_up(:)=(log(prs_p(:))+log(prs(:)))/2.
      tmp_vir(:)=(1.0_r8 + eps*qmr(:)) * tmp(:)
      prs_dw(:)=exp(ln_prs_dw(:))
      prs_up(:)=exp(ln_prs_up(:))
      thick(:)=Rd * tmp_vir(:) / grav * (ln_prs_dw(:) - ln_prs_up(:))
      ln_prs_dw(:)=ln_prs_up(:)
      fac(:)=thick(:) * prs(:) / (Rd * tmp(:))
!
      if(any(fac(:) == missing_r8)) then
         val(:)=missing_r8
         zstatus=20
         write(string1, *) &
         'APM: MODIS AOD TOTAL COL fac is missing '
         call error_handler(E_MSG, 'obs_def_MODIS_AOD_TOTAL_COL', string1, source)
         call track_status(ens_size, zstatus, val, istatus, return_now)
         return
      endif
!
      level=real(ilay)
      nloc = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
!
! so4 (ppmv) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_SO4, so4(:), zstatus)
!
! Convert from ppmv to kg/m^2
      so4(:) = so4(:) * 1.e-6 * fac(:)
!
! dust (ug/kg - dry air) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_DST01, dust1(:), zstatus)
!
! Convert from ug/kg to kg/m^2
      dust1(:) = dust1(:) * 1.e-9 * fac(:)
!
! dust (ug/kg - dry air) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_DST02, dust2(:), zstatus)
!
! Convert from ug/kg to kg/m^2
      dust2(:) = dust2(:) * 1.e-9 * fac(:)
!
! dust (ug/kg - dry air) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_DST03, dust3(:), zstatus)
!
! Convert from ug/kg to kg/m^2
      dust3(:) = dust3(:) * 1.e-9 * fac(:)
!
! dust (ug/kg - dry air) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_DST04, dust4(:), zstatus)
!
! Convert from ug/kg to kg/m^2
      dust4(:) = dust4(:) * 1.e-9 * fac(:)
!
! dust (ug/kg - dry air) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_DST05, dust5(:), zstatus)
!
! Convert from ug/kg to kg/m^2
      dust5(:) = dust5(:) * 1.e-9 * fac(:)
!
! hydrophilic black carbon (ug/kg - dry air) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_BC1, bc1(:), zstatus)
!
! Convert from ug/kg to kg/m^2
      bc1(:) = bc1(:) * 1.e-9 * fac(:)
!
! hydrophilic black carbon (ug/kg - dry air) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_BC2, bc2(:), zstatus)
!
! Convert from ug/kg to kg/m^2
      bc2(:) = bc2(:) * 1.e-9 * fac(:)
!
! hydrophilic organic carbon (ug/kg - dry air) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_OC1, oc1(:), zstatus)
!
! Convert from ug/kg to kg/m^2
      oc1(:) = oc1(:) * 1.e-9 * fac(:)
!
! hydrophilic organic carbon (ug/kg - dry air) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_OC2, oc2(:), zstatus)
!
! Convert from ug/kg to kg/m^2
      oc2(:) = oc2(:) * 1.e-9 * fac(:)
!
! sea salt (ug/kg - dry air) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_SSLT01, ss1(:), zstatus)
!
! Convert from ug/kg to kg/m^2
      ss1(:) = ss1(:) * 1.e-9 * fac(:)
!
! sea salt (ug/kg - dry air) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_SSLT02, ss2(:), zstatus)
!
! Convert from ug/kg to kg/m^2
      ss2(:) = ss2(:) * 1.e-9 * fac(:)
!
! sea salt (ug/kg - dry air) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_SSLT03, ss3(:), zstatus)
!
! Convert from ug/kg to kg/m^2
      ss3(:) = ss3(:) * 1.e-9 * fac(:)
!
! sea salt (ug/kg - dry air) at this location - this calls the model_mod code.
      zstatus(:)=0
      call interpolate(state_handle, ens_size, nloc, QTY_SSLT04, ss4(:), zstatus)
!
! Convert from ug/kg to kg/m^2
      ss4(:) = ss4(:) * 1.e-9 * fac(:)
!
! Check for missing values
   if(any(so4.lt.0.) .or. any(bc1.lt.0) .or. any(bc2.lt.0) .or. &
   any(oc1.lt.0) .or. any(oc2.lt.0) .or. any(dust1.lt.0) .or. any(dust2.lt.0) .or. &
   any(dust3.lt.0) .or. any(dust4.lt.0) .or. any(dust5.lt.0) .or. any(ss1.lt.0) .or. &
   any(ss2.lt.0) .or. any(ss3.lt.0) .or. any(ss4.lt.0)) then   
      zstatus(:)=20
      val(:)=missing_r8
!      write(string1, *) 'APM: MODIS AOD TOTAL CO FO has missing values '
!      call error_handler(E_MSG, routine, string1, source)
      call track_status(ens_size, zstatus, val, istatus, return_now)
      return
   endif
!
! Expected Modis AOD and convert from kg to g
      istatus=0      
      val(:) = val(:) + (ext_so4*so4(:) + ext_oc1*oc1(:) + ext_oc2*oc2(:) + &
      ext_bc1*bc1(:) + ext_bc2*bc2(:) + ext_ssl1*ss1(:) + ext_ssl2*ss2(:) + &
      ext_ssl3*ss3(:) + ext_ssl4*ss4(:) + ext_dst1*dust1(:) + ext_dst2*dust2(:) + &
      ext_dst3*dust3(:) + ext_dst4*dust4(:) + ext_dst5*dust5(:)) * 1.e3_r8
   enddo
!
! Check expected observation
   if(any(val.lt.0)) then
      zstatus(:)=20
      val(:)=missing_r8
!      write(string1, *) 'APM: MODIS AOD TOTAL COL expected value is negative '
!      call error_handler(E_MSG, 'obs_def_MODIS_AOD_TOTAL_COL', string1, source)
      call track_status(ens_size, zstatus, val, istatus, return_now)
      return
   endif
!
end subroutine get_expected_modis_aod_total_col
!
!----------------------------------------------------------------------
!
end module obs_def_MODIS_AOD_TOTAL_COL_mod
! END DART PREPROCESS MODULE CODE
