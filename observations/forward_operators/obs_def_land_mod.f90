! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! The obs_def_land_mod.f90 is intended to provide the forward operators
! necessary to convert land surface states to the expected values of the
! observations. Since the meaning of land model variables varies across
! land models, the forward operators here may try multiple approaches.


! The following observation kinds are currently created by the converters
! in the NASA_Earthdata directory. As NASA provides them
!	soil_moisture_x:long_name = "Volumetric Soil Moisture from X-band" ;
!	soil_moisture_x:units = "percent" ;
! CLM soil moisture in the restart files is another beast:
!       H2OSOI_LIQ:long_name = "liquid water" ;
!       H2OSOI_LIQ:units = "kg/m2" ;
! CLM soil moisture in the history files is another beast:
!       H2OSOI:long_name = "liquid water" ;
!       H2OSOI:units = "mm3/mm3" ;

!SMOS_A_SOIL_MOISTURE,         volumetric soil moisture  percent
!SMOS_D_SOIL_MOISTURE,         volumetric soil moisture  percent
!SMAP_A_SOIL_MOISTURE,         volumetric soil moisture  percent
!SMAP_D_SOIL_MOISTURE,         volumetric soil moisture  percent
!SSMI_A_SOIL_MOISTURE,         volumetric soil moisture  percent
!SSMI_D_SOIL_MOISTURE,         volumetric soil moisture  percent
!AMSRE_A_SOIL_MOISTURE_X,      volumetric soil moisture  percent (cvrts to a fraction)
!AMSRE_D_SOIL_MOISTURE_X,      volumetric soil moisture  percent (cvrts to a fraction)
!AMSRE_A_SOIL_MOISTURE_C,      volumetric soil moisture  percent (cvrts to a fraction)
!AMSRE_D_SOIL_MOISTURE_C,      volumetric soil moisture  percent (cvrts to a fraction)
!TRMM_SOIL_MOISTURE,           volumetric soil moisture  percent
!WINDSAT_SOIL_MOISTURE_X,      volumetric soil moisture  percent
!WINDSAT_SOIL_MOISTURE_C,      volumetric soil moisture  percent

! BEGIN DART PREPROCESS KIND LIST
!SOIL_TEMPERATURE,               QTY_SOIL_TEMPERATURE,           COMMON_CODE
!LPRM_SOIL_MOISTURE,             QTY_SOIL_MOISTURE,              COMMON_CODE
!SMOS_A_SOIL_MOISTURE,           QTY_SOIL_MOISTURE,              COMMON_CODE
!SMOS_D_SOIL_MOISTURE,           QTY_SOIL_MOISTURE,              COMMON_CODE
!SMAP_A_SOIL_MOISTURE,           QTY_SOIL_MOISTURE,              COMMON_CODE
!SMAP_D_SOIL_MOISTURE,           QTY_SOIL_MOISTURE,              COMMON_CODE
!SSMI_A_SOIL_MOISTURE,           QTY_SOIL_MOISTURE,              COMMON_CODE
!SSMI_D_SOIL_MOISTURE,           QTY_SOIL_MOISTURE,              COMMON_CODE
!AMSRE_A_SOIL_MOISTURE_X,        QTY_SOIL_MOISTURE,              COMMON_CODE
!AMSRE_D_SOIL_MOISTURE_X,        QTY_SOIL_MOISTURE,              COMMON_CODE
!AMSRE_A_SOIL_MOISTURE_C,        QTY_SOIL_MOISTURE,              COMMON_CODE
!AMSRE_D_SOIL_MOISTURE_C,        QTY_SOIL_MOISTURE,              COMMON_CODE
!TRMM_SOIL_MOISTURE,             QTY_SOIL_MOISTURE,              COMMON_CODE
!WINDSAT_SOIL_MOISTURE_X,        QTY_SOIL_MOISTURE,              COMMON_CODE
!WINDSAT_SOIL_MOISTURE_C,        QTY_SOIL_MOISTURE,              COMMON_CODE
!WATER_TABLE_DEPTH,              QTY_WATER_TABLE_DEPTH,          COMMON_CODE
!SOIL_TEMPERATURE,               QTY_SOIL_TEMPERATURE,           COMMON_CODE
!SOIL_MOISTURE,                  QTY_SOIL_MOISTURE,              COMMON_CODE
!LAYER_LIQUID_WATER,             QTY_SOIL_LIQUID_WATER,          COMMON_CODE
!LAYER_ICE,                      QTY_SOIL_ICE,                   COMMON_CODE
!SNOW_THICKNESS,                 QTY_SNOW_THICKNESS,             COMMON_CODE
!SNOW_WATER,                     QTY_SNOW_WATER,                 COMMON_CODE
!MODIS_SNOWCOVER_FRAC,           QTY_SNOWCOVER_FRAC,             COMMON_CODE
!MODIS_LEAF_AREA_INDEX,          QTY_LEAF_AREA_INDEX,            COMMON_CODE
!GIMMS_LEAF_AREA_INDEX,          QTY_LEAF_AREA_INDEX,            COMMON_CODE
!SP_LEAF_AREA_INDEX,             QTY_LEAF_AREA_INDEX,            COMMON_CODE
!MODIS_FPAR,                     QTY_FRACTION_ABSORBED_PAR
!BIOMASS,                        QTY_BIOMASS
!LEAF_CARBON,                    QTY_LEAF_CARBON,                COMMON_CODE
!LIVE_STEM_CARBON,               QTY_LIVE_STEM_CARBON,           COMMON_CODE
!DEAD_STEM_CARBON,               QTY_DEAD_STEM_CARBON,           COMMON_CODE
!LEAF_AREA_INDEX,                QTY_LEAF_AREA_INDEX,            COMMON_CODE
!LEAF_NITROGEN,                  QTY_LEAF_NITROGEN,              COMMON_CODE
!TOWER_AIR_TEMPERATURE,          QTY_TEMPERATURE,                COMMON_CODE
!TREE_RING_NPP_FLUX,             QTY_NET_PRIMARY_PROD_FLUX,      COMMON_CODE
!TOWER_SOIL_TEMPERATURE,         QTY_TEMPERATURE,                COMMON_CODE
!TOWER_U_WIND_COMPONENT,         QTY_U_WIND_COMPONENT,           COMMON_CODE
!TOWER_V_WIND_COMPONENT,         QTY_V_WIND_COMPONENT,           COMMON_CODE
!TOWER_GLOBAL_RADIATION,         QTY_RADIATION,                  COMMON_CODE
!TOWER_NET_CARBON_FLUX,          QTY_NET_CARBON_FLUX,            COMMON_CODE
!SURFACE_ALBEDO,                 QTY_SURFACE_ALBEDO
!OCO2_SIF,                       QTY_SOLAR_INDUCED_FLUORESCENCE, COMMON_CODE
!ECOSTRESS_ET,                   QTY_LATENT_HEAT_FLUX,           COMMON_CODE
! END DART PREPROCESS KIND LIST

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_land_mod, only : calculate_albedo, &
!                               calculate_biomass, &
!                               calculate_fpar
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!  case(SURFACE_ALBEDO)
!     call calculate_albedo(state_handle, ens_size, location, expected_obs, istatus)
!  case(BIOMASS)
!     call calculate_biomass(state_handle, ens_size, location, expected_obs, istatus)
!  case(MODIS_FPAR)
!     call calculate_fpar(state_handle, ens_size, location, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!    case(SURFACE_ALBEDO, &
!         BIOMASS, &
!         MODIS_FPAR)
!       continue
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!    case(SURFACE_ALBEDO, &
!         BIOMASS, &
!         MODIS_FPAR)
!       continue
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!    case(SURFACE_ALBEDO, &
!         BIOMASS, &
!         MODIS_FPAR)
!       continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_land_mod

! This is the basic module for the forward (observation) operators commonly
! used with land models.

use            types_mod, only : r8, MISSING_R8

use         location_mod, only : location_type, &
                                 write_location

use        utilities_mod, only : register_module, &
                                 error_handler, &
                                 E_ERR, E_MSG, &
                                 do_output

use      assim_model_mod, only : interpolate

use ensemble_manager_mod, only : ensemble_type

use         obs_kind_mod, only : QTY_RADIATION_VISIBLE_DOWN, &
                                 QTY_RADIATION_NEAR_IR_DOWN, &
                                 QTY_RADIATION_VISIBLE_UP, &
                                 QTY_RADIATION_NEAR_IR_UP, &
                                 QTY_LIVE_STEM_CARBON, &
                                 QTY_DEAD_STEM_CARBON, &
                                 QTY_LEAF_CARBON, &
                                 QTY_FRACTION_ABSORBED_PAR, &
                                 QTY_PAR_DIRECT, &
                                 QTY_PAR_DIFFUSE, &
                                 QTY_ABSORBED_PAR

implicit none
private

public :: calculate_albedo, &
          calculate_biomass, &
          calculate_fpar

character(len=*), parameter :: source   = 'obs_def_land_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical :: module_initialized = .false.

character(len=512) :: string1, string2, string3

! This might be useful, but not enough to warrant a namelist ... yet
logical :: debug = .false.


!===============================================================================
contains
!===============================================================================


subroutine initialize_module()

if (module_initialized) return

module_initialized = .true.

! Log the version of this source file.
call register_module(source, revision, revdate)

end subroutine initialize_module


!===============================================================================


subroutine calculate_albedo(state_handle, ens_size, location, obs_val, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: obs_val(ens_size)
integer,             intent(out) :: istatus(ens_size)

real(r8) ::  visible_in(ens_size)
real(r8) :: visible_out(ens_size)
real(r8) ::      nir_in(ens_size)
real(r8) ::     nir_out(ens_size)
real(r8) ::       numer(ens_size)
real(r8) ::       denom(ens_size)

integer :: stat(ens_size,4)
integer :: imem

istatus = 1           ! 0 == success, anything else is a failure
obs_val = MISSING_R8

call error_handler(E_ERR,'calculate_albedo','routine untested - stopping.', &
           source, revision, revdate)

if ( .not. module_initialized ) call initialize_module()

! Intentionally try to compute all required components before failing.
! The desire is to inform about ALL failed components instead of failing
! one-by-one.

call interpolate(state_handle, ens_size, location, QTY_RADIATION_VISIBLE_DOWN, &
        visible_in,  stat(:,1))
call interpolate(state_handle, ens_size, location, QTY_RADIATION_NEAR_IR_DOWN, &
        nir_in,      stat(:,2))
call interpolate(state_handle, ens_size, location, QTY_RADIATION_VISIBLE_UP, &
        visible_out, stat(:,3))
call interpolate(state_handle, ens_size, location, QTY_RADIATION_NEAR_IR_UP, &
        nir_out,     stat(:,4))

if (any(stat /= 0)) then
   istatus = stat(:,1)*1000 + stat(:,2)*100 + stat(:,3)*10 + stat(:,4)
   return
endif

numer = visible_out + nir_out
denom = visible_in  + nir_in

if (any(denom <= tiny(0.0_r8))) then
   call write_location(42,location,fform='FORMATTED',charstring=string1)
   call error_handler(E_MSG,'calculate_albedo','no incoming radiation',text2=string1)
   return
endif

obs_val = numer / denom

istatus = 0   ! success

if (debug .and. do_output()) then
   do imem = 1,ens_size
      write(string1,*)'incoming (visible nir sum status) ', &
         visible_in(imem), nir_in(imem), denom(imem), stat(imem,1)*1000 + stat(imem,2)*100
      write(string2,*)'outgoing (visible nir sum status) ', &
         visible_out(imem), nir_out(imem), numer(imem), stat(imem,3)*10 + stat(imem,4)
      write(string3,*)'albedo ', obs_val(imem)
      call error_handler(E_MSG,'calculate_albedo:',string1,text2=string2,text3=string3)
   enddo
endif

end subroutine calculate_albedo


!===============================================================================


subroutine calculate_biomass(state_handle, ens_size, location, obs_val, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: obs_val(ens_size)
integer,             intent(out) :: istatus(ens_size)

real(r8) ::      leaf_carbon(ens_size)
real(r8) :: live_stem_carbon(ens_size)
real(r8) :: dead_stem_carbon(ens_size)
integer  :: stat(ens_size,3)
integer  :: imem

istatus = 1           ! 0 == success, anything else is a failure
obs_val = MISSING_R8

if ( .not. module_initialized ) call initialize_module()

! Intentionally try to compute all required components before failing.
! The desire is to inform about ALL failed components instead of failing
! one-by-one.

call interpolate(state_handle, ens_size, location, QTY_LEAF_CARBON, &
                 leaf_carbon,      stat(:,1))
call interpolate(state_handle, ens_size, location, QTY_LIVE_STEM_CARBON, &
                 live_stem_carbon, stat(:,2))
call interpolate(state_handle, ens_size, location, QTY_DEAD_STEM_CARBON, &
                 dead_stem_carbon, stat(:,3))

if (any(stat /= 0)) then
   istatus = stat(:,1)*100 + stat(:,2)*10 + stat(:,3)
   return
endif

obs_val = leaf_carbon + live_stem_carbon + dead_stem_carbon

istatus = 0   ! success

if (debug .and. do_output()) then
   do imem = 1,ens_size
      write(string1,*)'biomass for ensemble member ', imem
      write(string2,*)'carbon: leaf,live_stem,dead_stem', &
         leaf_carbon(imem), live_stem_carbon(imem), dead_stem_carbon(imem)
      write(string3,*)'status: leaf,live_stem,dead_stem', stat(imem,:)
      call error_handler(E_MSG,'calculate_biomass:',string1,text2=string2,text3=string3)
   enddo
endif

end subroutine calculate_biomass


!===============================================================================
!> calculate the fraction of photosynthetically available radiation
!> The MODIS website https://modis.gsfc.nasa.gov/data/dataprod/mod15.php
!> states: 'FPAR is the fraction of photosynthetically active radiation
!>          (400-700 nm) absorbed by green vegetation.'


subroutine calculate_fpar(state_handle, ens_size, location, obs_val, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: obs_val(ens_size)
integer,             intent(out) :: istatus(ens_size)

real(r8) :: diffuse(ens_size)
real(r8) :: active(ens_size)
real(r8) :: absorbed(ens_size)
integer  :: stat(ens_size,3)
integer  :: imem
real(r8) :: denom

istatus = 1           ! 0 == success, anything else is a failure
obs_val = MISSING_R8

if ( .not. module_initialized ) call initialize_module()

! If the model state has it directly, this is simple.
! If it does not, we try to calculate it from what is available and what is absorbed

call interpolate(state_handle, ens_size, location, QTY_FRACTION_ABSORBED_PAR, &
                 obs_val, istatus)

if (all(istatus == 0)) return

! Intentionally try to compute all required components before failing.
! This is the part that needs scientific direction ...

call interpolate(state_handle, ens_size, location, QTY_PAR_DIRECT,   &
                 active,   stat(:,1))
call interpolate(state_handle, ens_size, location, QTY_PAR_DIFFUSE,   &
                 diffuse,  stat(:,2))
call interpolate(state_handle, ens_size, location, QTY_ABSORBED_PAR, &
                 absorbed, stat(:,3))

if (any(stat /= 0)) then
   istatus = stat(:,1)*1000 + stat(:,2)*100 + stat(:,3)
   return
endif

do imem = 1,ens_size

   ! If any of them are missing it is cause for failure and an early return
   if (absorbed(imem) == MISSING_R8 .or.  &
         active(imem) == MISSING_R8 .or.  &
        diffuse(imem) == MISSING_R8 ) then
      write(string1,*)'member',imem,'absorbed',absorbed(imem),'status',istatus(imem)
      write(string2,*)'values: active, absorbed',active(imem),diffuse(imem)
      call error_handler(E_MSG,'calculate_fpar:MISSING',string1,text2=string2)
      istatus(imem) = imem
      return
   endif 

   denom = active(imem) + diffuse(imem)

   if (absorbed(imem) < tiny(denom)) then
      obs_val(imem) = 0.0_r8
      istatus(imem) = 0
   elseif (denom <= tiny(denom)) then ! avoid dividing by zeroish
      write(string1,*)'member ',imem,' denom ',denom
      write(string2,*)'values: active, diffuse ',active(imem),diffuse(imem)
      call error_handler(E_MSG,'calculate_fpar:ZERO',string1,text2=string2)
      istatus(imem) = imem
   else
      obs_val(imem) = min(absorbed(imem) / denom, 1.0_r8)
      istatus(imem) = 0
   endif 

enddo

if (debug .and. do_output()) then
   do imem = 1,ens_size
      write(string1,*)'member ',imem,' fpar ',obs_val(imem),' status ',istatus(imem)
      write(string2,*)'values: active, diffuse, absorbed ', &
                      active(imem), diffuse(imem), absorbed(imem)
      write(string3,*)'status: active, diffuse, absorbed ',stat(imem,:)
      call error_handler(E_MSG,'calculate_fpar:',string1,text2=string2,text3=string3)
   enddo
endif

end subroutine calculate_fpar

end module obs_def_land_mod

! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------
