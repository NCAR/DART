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
! MONITOR_SO2, QTY_SO2
! MONITOR_NO2, QTY_NO2
! MONITOR_PM10, QTY_PM10
! MONITOR_CO, QTY_CO
! MONITOR_O3, QTY_O3
! MONITOR_PM25, QTY_PM25
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_monitor_mod, only : write_monitor_so2, read_monitor_so2, &
!                                  interactive_monitor_so2, get_expected_monitor_so2
!                                  
!   use obs_def_monitor_mod, only : write_monitor_no2, read_monitor_no2, &
!                                  interactive_monitor_no2, get_expected_monitor_no2
!                                  
!   use obs_def_monitor_mod, only : write_monitor_pm10, read_monitor_pm10, &
!                                  interactive_monitor_pm10, get_expected_monitor_pm10
!                                  
!   use obs_def_monitor_mod, only : write_monitor_co, read_monitor_co, &
!                                  interactive_monitor_co, get_expected_monitor_co
!                                  
!   use obs_def_monitor_mod, only : write_monitor_o3, read_monitor_o3, &
!                                  interactive_monitor_o3, get_expected_monitor_o3
!                                  
!   use obs_def_monitor_mod, only : write_monitor_pm25, read_monitor_pm25, &
!                                  interactive_monitor_pm25, get_expected_monitor_pm25
!                                  
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(MONITOR_SO2)                                                           
!            call get_expected_monitor_so2(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
!         case(MONITOR_NO2)                                                           
!            call get_expected_monitor_no2(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
!         case(MONITOR_PM10)                                                           
!            call get_expected_monitor_pm10(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
!         case(MONITOR_CO)                                                           
!            call get_expected_monitor_co(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
!         case(MONITOR_O3)                                                           
!            call get_expected_monitor_o3(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
!         case(MONITOR_PM25)                                                           
!            call get_expected_monitor_pm25(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(MONITOR_SO2)
!         call read_monitor_so2(obs_def%key, ifile, fileformat)
!      case(MONITOR_NO2)
!         call read_monitor_no2(obs_def%key, ifile, fileformat)
!      case(MONITOR_PM10)
!         call read_monitor_pm10(obs_def%key, ifile, fileformat)
!      case(MONITOR_CO)
!         call read_monitor_co(obs_def%key, ifile, fileformat)
!      case(MONITOR_O3)
!         call read_monitor_o3(obs_def%key, ifile, fileformat)
!      case(MONITOR_PM25)
!         call read_monitor_pm25(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(MONITOR_SO2)
!         call write_monitor_so2(obs_def%key, ifile, fileformat)
!      case(MONITOR_NO2)
!         call write_monitor_no2(obs_def%key, ifile, fileformat)
!      case(MONITOR_PM10)
!         call write_monitor_pm10(obs_def%key, ifile, fileformat)
!      case(MONITOR_CO)
!         call write_monitor_co(obs_def%key, ifile, fileformat)
!      case(MONITOR_O3)
!         call write_monitor_o3(obs_def%key, ifile, fileformat)
!      case(MONITOR_PM25)
!         call write_monitor_pm25(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(MONITOR_SO2)
!         call interactive_monitor_so2(obs_def%key)
!      case(MONITOR_NO2)
!         call interactive_monitor_no2(obs_def%key)
!      case(MONITOR_PM10)
!         call interactive_monitor_pm10(obs_def%key)
!      case(MONITOR_CO)
!         call interactive_monitor_co(obs_def%key)
!      case(MONITOR_O3)
!         call interactive_monitor_o3(obs_def%key)
!      case(MONITOR_PM25)
!         call interactive_monitor_pm25(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_monitor_mod

use        types_mod, only : r8, missing_r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             nmlfileunit, check_namelist_read, &
                             find_namelist_in_file, do_nml_file, do_nml_term, &
                             ascii_file_format
use     location_mod, only : location_type, set_location, get_location , write_location, &
                             read_location

use  assim_model_mod, only : interpolate
use    obs_kind_mod,  only : QTY_SO2, QTY_NO2, QTY_CO, QTY_O3, QTY_BC1, QTY_BC2, QTY_OC1, QTY_OC2, &
                             QTY_DST01, QTY_DST02, QTY_DST03, QTY_DST04, &
                             QTY_DST05, QTY_SO4, QTY_SSLT01, QTY_SSLT02, QTY_SSLT03, &
                             QTY_SSLT04, QTY_PM25, QTY_PM10,QTY_PRESSURE, QTY_TEMPERATURE 
use  ensemble_manager_mod, only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private

public :: write_monitor_so2,  read_monitor_so2, &
          write_monitor_no2,  read_monitor_no2, &
          write_monitor_co,   read_monitor_co, &
          write_monitor_o3,   read_monitor_o3, &
          write_monitor_pm10, read_monitor_pm10, &
          write_monitor_pm25, read_monitor_pm25, &
          interactive_monitor_so2,  get_expected_monitor_so2, &
          interactive_monitor_no2,  get_expected_monitor_no2, &
          interactive_monitor_co,   get_expected_monitor_co, &
          interactive_monitor_o3,   get_expected_monitor_o3, &
          interactive_monitor_pm10, get_expected_monitor_pm10, &
          interactive_monitor_pm25, get_expected_monitor_pm25

logical, parameter :: use_diag = .false.

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_monitor_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical, save :: module_initialized = .false.

contains

!----------------------------------------------------------------------
!>

subroutine initialize_module

! Prevent multiple calls from executing this code more than once.
if (module_initialized) return

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module


!----------------------------------------------------------------------
!>


subroutine read_monitor_so2(key, ifile, fform)

integer, intent(out)            :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional    :: fform

continue

end subroutine read_monitor_so2


!----------------------------------------------------------------------
!>


subroutine write_monitor_so2(key, ifile, fform)

integer, intent(in)             :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional :: fform

continue

end subroutine write_monitor_so2


!----------------------------------------------------------------------
!>


subroutine interactive_monitor_so2(key)

integer, intent(out) :: key

continue

end subroutine interactive_monitor_so2


!----------------------------------------------------------------------
!>


subroutine get_expected_monitor_so2(state_handle, ens_size, location, key, so2, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: key
real(r8),            intent(out) :: so2(ens_size)     ! so2 concentration(ug/m3)
integer,             intent(out) :: istatus(ens_size)

real(r8), PARAMETER :: mso2=64   ! molecular weight
                                 !   to avoid problems near zero in Bolton's equation
real(r8) :: p_Pa(ens_size)       ! pressure (Pa)
real(r8) :: T_k(ens_size)        ! temperature(K)
integer :: p_istatus(ens_size)
integer :: t_istatus(ens_size)
integer :: so2_istatus(ens_size)
logical  :: return_now

if ( .not. module_initialized ) call initialize_module

call interpolate(state_handle, ens_size, location, QTY_PRESSURE, p_Pa, p_istatus)
call track_status(ens_size, p_istatus, p_Pa, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, T_k, istatus)
call track_status(ens_size, t_istatus, T_k, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_SO2, so2, so2_istatus)! so2(ppmv)
call track_status(ens_size, so2_istatus, so2, istatus, return_now)
if (return_now) return

so2=so2*(mso2*p_pa)/(8.314*T_k) ! so2(ug/m3)


end subroutine get_expected_monitor_so2
 

!----------------------------------------------------------------------
!>


subroutine read_monitor_no2(key, ifile, fform)

integer, intent(out)            :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional    :: fform

continue

end subroutine read_monitor_no2


!----------------------------------------------------------------------
!>


subroutine write_monitor_no2(key, ifile, fform)


integer, intent(in)             :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional :: fform

continue

end subroutine write_monitor_no2


!----------------------------------------------------------------------
!>


subroutine interactive_monitor_no2(key)

integer, intent(out) :: key

continue

end subroutine interactive_monitor_no2


!----------------------------------------------------------------------
!>


subroutine get_expected_monitor_no2(state_handle, ens_size, location, key, no2, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: key
real(r8),            intent(out) :: no2(ens_size)              ! no2 concentration(ug/m3)
integer,             intent(out) :: istatus(ens_size)

real(r8), PARAMETER :: mno2=46   ! molecular weight
                                          !   to avoid problems near zero in Bolton's equation
real(r8) :: p_Pa(ens_size)                ! pressure (Pa)
real(r8) :: T_k(ens_size)                  ! tempreture(K)
logical :: return_now
integer :: p_istatus(ens_size)
integer :: t_istatus(ens_size)
integer :: no2_istatus(ens_size)

if ( .not. module_initialized ) call initialize_module

call interpolate(state_handle, ens_size, location, QTY_PRESSURE, p_Pa, p_istatus)
call track_status(ens_size, p_istatus, p_Pa, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, T_k, t_istatus)
call track_status(ens_size, t_istatus, T_k, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_NO2, no2, no2_istatus)! no2(ppmv)
call track_status(ens_size, no2_istatus, no2,istatus, return_now)
if (return_now) return

no2=no2*(mno2*p_pa)/(8.314*T_k) ! no2(ug/m3)

end subroutine get_expected_monitor_no2 


!----------------------------------------------------------------------
!>


subroutine read_monitor_o3(key, ifile, fform)

integer, intent(out)            :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional    :: fform

continue

end subroutine read_monitor_o3


!----------------------------------------------------------------------
!>


subroutine write_monitor_o3(key, ifile, fform)

integer, intent(in)             :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional :: fform

continue

end subroutine write_monitor_o3


!----------------------------------------------------------------------
!>


subroutine interactive_monitor_o3(key)

integer, intent(out) :: key

continue

end subroutine interactive_monitor_o3


!----------------------------------------------------------------------
!>


subroutine get_expected_monitor_o3(state_handle, ens_size, location, key, o3, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: key
real(r8),            intent(out) :: o3(ens_size)               ! o3 concentration(ug/m3)
integer,             intent(out) :: istatus(ens_size)

real(r8), PARAMETER :: mo3=48   ! molecular weight
                                          !   to avoid problems near zero in Bolton's equation
real(r8) :: p_Pa(ens_size)                ! pressure (Pa)
real(r8) :: T_k(ens_size)                  ! tempreture(K)
logical :: return_now
integer :: p_istatus(ens_size)
integer :: t_istatus(ens_size)
integer :: o3_istatus(ens_size)

if ( .not. module_initialized ) call initialize_module

call interpolate(state_handle, ens_size, location, QTY_PRESSURE, p_Pa, p_istatus)
call track_status(ens_size, p_istatus, p_Pa, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, T_k, t_istatus)
call track_status(ens_size, t_istatus, T_k, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_O3, o3, o3_istatus)! o3(ppmv)
call track_status(ens_size, o3_istatus, o3, istatus, return_now)
if (return_now) return

o3=o3*(mo3*p_pa)/(8.314*T_k) ! o3(ug/m3)


end subroutine get_expected_monitor_o3


!----------------------------------------------------------------------
!>


subroutine read_monitor_co(key, ifile, fform)

integer, intent(out)            :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional    :: fform

continue

end subroutine read_monitor_co


!----------------------------------------------------------------------
!>


subroutine write_monitor_co(key, ifile, fform)

integer, intent(in)             :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional :: fform

continue

end subroutine write_monitor_co


!----------------------------------------------------------------------
!>


subroutine interactive_monitor_co(key)

integer, intent(out) :: key

continue

end subroutine interactive_monitor_co


!----------------------------------------------------------------------
!>


subroutine get_expected_monitor_co(state_handle, ens_size, location, key, co, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: key
real(r8),            intent(out) :: co(ens_size)              ! co concentration(mg/m3)
integer,             intent(out) :: istatus(ens_size)

real(r8), PARAMETER :: mco=28   ! molecular weight
                                          !   to avoid problems near zero in Bolton's equation
real(r8) :: p_Pa(ens_size)                ! pressure (Pa)
real(r8) :: T_k(ens_size)                 ! tempreture(K)
logical :: return_now
integer :: p_istatus(ens_size)
integer :: t_istatus(ens_size)
integer :: co_istatus(ens_size)

if ( .not. module_initialized ) call initialize_module

call interpolate(state_handle, ens_size, location, QTY_PRESSURE, p_Pa, p_istatus)
call track_status(ens_size, p_istatus, p_Pa, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, T_k, t_istatus)
call track_status(ens_size, t_istatus, T_k, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_CO, co, co_istatus)! co(ppmv)
call track_status(ens_size, co_istatus, co, istatus, return_now)
if (return_now) return

co=co*(mco*p_pa)/(1000*8.314*T_k) ! co(mg/m3)

end subroutine get_expected_monitor_co


!----------------------------------------------------------------------
!>


subroutine read_monitor_pm25(key, ifile, fform)

integer, intent(out)            :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional    :: fform

continue

end subroutine read_monitor_pm25


!----------------------------------------------------------------------
!>


subroutine write_monitor_pm25(key, ifile, fform)

integer, intent(in)             :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional :: fform

continue

end subroutine write_monitor_pm25


!----------------------------------------------------------------------
!>


subroutine interactive_monitor_pm25(key)

integer, intent(out) :: key

continue

end subroutine interactive_monitor_pm25


!----------------------------------------------------------------------
!>


subroutine get_expected_monitor_pm25(state_handle, ens_size, location, key, ppm25, istatus)
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: key
real(r8),            intent(out) :: ppm25(ens_size)               ! pm2.5 concentration(ug/3)
integer,             intent(out) :: istatus(ens_size)

real(r8), PARAMETER :: m_dry_air=29   ! molecular weight
real(r8), PARAMETER :: mso4=96   ! molecular weight
                                          !   to avoid problems near zero in Bolton's equation
real(r8) :: p_Pa(ens_size)                ! pressure (Pa)
real(r8) :: T_k(ens_size)                  ! tempreture(K)
real(r8), dimension(ens_size) :: PM25,BC1,BC2,DST01,DST02,SSLT01,SSLT02,SO4,OC1,OC2 !components
real(r8) :: alt_temp(ens_size)
logical :: return_now
integer :: p_istatus(ens_size)
integer :: t_istatus(ens_size)
integer :: pm25_istatus(ens_size)
integer :: bc1_istatus(ens_size)
integer :: bc2_istatus(ens_size)
integer :: dst01_istatus(ens_size)
integer :: dst02_istatus(ens_size)
integer :: sslt01_istatus(ens_size)
integer :: sslt02_istatus(ens_size)
integer :: so4_istatus(ens_size)
integer :: oc1_istatus(ens_size)
integer :: oc2_istatus(ens_size)

if ( .not. module_initialized ) call initialize_module

istatus(:)=0
ppm25=0.
call interpolate(state_handle, ens_size, location, QTY_PRESSURE, p_Pa, p_istatus)
call track_status(ens_size, p_istatus, ppm25, istatus, return_now)
if (return_now) return

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, T_k, t_istatus)
call track_status(ens_size, t_istatus, ppm25, istatus, return_now)
if (return_now) return

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_PM25, PM25, pm25_istatus)
call track_status(ens_size, pm25_istatus, ppm25, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm25 = missing_r8
elsewhere
   ppm25=ppm25+PM25
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_BC1, BC1, bc1_istatus)
call track_status(ens_size, bc1_istatus, ppm25, istatus, return_now)
if (return_now) return
where (istatus /= 0) 
   ppm25 = missing_r8
elsewhere
   ppm25=ppm25+BC1
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_BC2, BC2, bc2_istatus)
call track_status(ens_size, bc2_istatus, ppm25, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm25 = missing_r8
elsewhere
   ppm25=ppm25+BC2
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_DST01, DST01, dst01_istatus)
call track_status(ens_size, dst01_istatus, ppm25, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm25 = missing_r8
elsewhere
   ppm25=ppm25+DST01
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_DST02, DST02, dst02_istatus)
call track_status(ens_size, dst02_istatus, ppm25, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm25 = missing_r8
elsewhere
   ppm25=ppm25+DST02*0.286
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_SSLT01, SSLT01, sslt01_istatus)
call track_status(ens_size, sslt01_istatus, ppm25, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm25 = missing_r8
elsewhere
   ppm25=ppm25+SSLT01
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_SSLT02, SSLT02, sslt02_istatus)
call track_status(ens_size, sslt02_istatus, ppm25, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm25 = missing_r8
elsewhere
   ppm25=ppm25+SSLT02*0.942
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_SO4, SO4, so4_istatus)
call track_status(ens_size, so4_istatus, ppm25, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm25 = missing_r8
elsewhere
   ppm25=ppm25+SO4*(mso4/m_dry_air)*1000.0*1.375
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_OC1, OC1, oc1_istatus)
call track_status(ens_size, oc1_istatus, ppm25, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm25 = missing_r8
elsewhere
   ppm25=ppm25+OC1*1.8
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_OC2, OC2, oc2_istatus)
call track_status(ens_size, oc2_istatus, ppm25, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm25 = missing_r8
elsewhere
   ppm25=ppm25+OC2*1.8 ! ug/kg
end where
alt_temp=1.0/(m_dry_air*p_Pa/(1000.0*8.314*T_k))
ppm25=ppm25/alt_temp ! pm25(ug/m3)


end subroutine get_expected_monitor_pm25


!----------------------------------------------------------------------
!>


subroutine read_monitor_pm10(key, ifile, fform)

integer, intent(out)            :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional    :: fform

continue

end subroutine read_monitor_pm10


!----------------------------------------------------------------------
!>


subroutine write_monitor_pm10(key, ifile, fform)

integer, intent(in)             :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional :: fform

continue

end subroutine write_monitor_pm10


!----------------------------------------------------------------------
!>


subroutine interactive_monitor_pm10(key)

integer, intent(out) :: key

continue

end subroutine interactive_monitor_pm10


!----------------------------------------------------------------------
!>


subroutine get_expected_monitor_pm10(state_handle, ens_size, location, key, ppm10, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: key
real(r8),            intent(out) :: ppm10(ens_size)              ! pm10 concentration(ug/3)
integer,             intent(out) :: istatus(ens_size)

real(r8), PARAMETER :: m_dry_air=29   ! molecular weight
real(r8), PARAMETER :: mso4=96   ! molecular weight
                                          !   to avoid problems near zero in Bolton's equation
real(r8) :: p_Pa(ens_size)                ! pressure (Pa)
real(r8) :: T_k(ens_size)                  ! tempreture(K)
real(r8), dimension(ens_size) :: PM25,BC1,BC2,DST01,DST02,DST03,DST04,SSLT01,SSLT02,SSLT03,SO4,OC1,OC2,PM10 !components
real(r8) :: alt_temp(ens_size)
logical :: return_now
integer :: p_istatus(ens_size)
integer :: t_istatus(ens_size)
integer :: pm10_istatus(ens_size)
integer :: pm25_istatus(ens_size)
integer :: bc1_istatus(ens_size)
integer :: bc2_istatus(ens_size)
integer :: dst01_istatus(ens_size)
integer :: dst02_istatus(ens_size)
integer :: dst03_istatus(ens_size)
integer :: dst04_istatus(ens_size)
integer :: sslt01_istatus(ens_size)
integer :: sslt02_istatus(ens_size)
integer :: sslt03_istatus(ens_size)
integer :: so4_istatus(ens_size)
integer :: oc1_istatus(ens_size)
integer :: oc2_istatus(ens_size)

if ( .not. module_initialized ) call initialize_module

istatus(:)=0
ppm10=0.
call interpolate(state_handle, ens_size, location, QTY_PRESSURE, p_Pa, p_istatus)
call track_status(ens_size, p_istatus, pm10, istatus, return_now)
if (return_now) return

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, T_k, t_istatus)
call track_status(ens_size, t_istatus, pm10, istatus, return_now)
if (return_now) return

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_PM25, PM25, pm25_istatus)
call track_status(ens_size, pm25_istatus, pm10, istatus, return_now)
if (return_now) return
where ( istatus /= 0) 
   ppm10 = missing_r8
elsewhere
   ppm10=ppm10+PM25
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_BC1, BC1, bc1_istatus)
call track_status(ens_size, bc1_istatus, pm10, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm10 = missing_r8
elsewhere
   ppm10=ppm10+BC1
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_BC2, BC2, bc2_istatus)
call track_status(ens_size, bc2_istatus, pm10, istatus, return_now)
if (return_now) return
where (istatus /= 0) 
   ppm10 = missing_r8
elsewhere
   ppm10=ppm10+BC2
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_DST01, DST01, dst01_istatus)
call track_status(ens_size, dst01_istatus, pm10, istatus, return_now)
if (return_now) return
where (istatus /= 0) 
   ppm10 = missing_r8
elsewhere
   ppm10=ppm10+DST01
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_DST02, DST02, dst02_istatus)
call track_status(ens_size, dst02_istatus, pm10, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm10 = missing_r8
elsewhere
   ppm10=ppm10+DST02
end where 

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_DST03, DST03, dst03_istatus)
call track_status(ens_size, dst03_istatus, pm10, istatus, return_now)
if (return_now) return
where (istatus /= 0) 
   ppm10 = missing_r8
elsewhere
   ppm10=ppm10+DST03
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_DST04, DST04, dst04_istatus)
call track_status(ens_size, dst04_istatus, pm10, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm10 = missing_r8
elsewhere
   ppm10=ppm10+DST04*0.87
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_SSLT01, SSLT01, sslt01_istatus)
call track_status(ens_size, sslt01_istatus, pm10, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm10 = missing_r8
elsewhere
   ppm10=ppm10+SSLT01
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_SSLT02, SSLT02, sslt02_istatus)
call track_status(ens_size, sslt02_istatus, pm10, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm10 = missing_r8
elsewhere
   ppm10=ppm10+SSLT02
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_SSLT03, SSLT03, sslt03_istatus)
call track_status(ens_size, sslt03_istatus, pm10, istatus, return_now)
if (return_now) return
where (istatus /= 0) 
   ppm10 = missing_r8
elsewhere
   ppm10=ppm10+SSLT03
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_SO4, SO4, so4_istatus)
call track_status(ens_size, so4_istatus, pm10, istatus, return_now)
if (return_now) return
where (istatus /= 0) 
   ppm10 = missing_r8
elsewhere
   ppm10=ppm10+SO4*(mso4/m_dry_air)*1000.0*1.375
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_OC1, OC1, oc1_istatus)
call track_status(ens_size, oc1_istatus, pm10, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm10 = missing_r8
elsewhere
   ppm10=ppm10+OC1*1.8
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_OC2, OC2, oc2_istatus)
call track_status(ens_size, oc2_istatus, pm10, istatus, return_now)
if (return_now) return
where (istatus /= 0) 
   ppm10 = missing_r8
elsewhere
   ppm10=ppm10+OC2*1.8 
end where

istatus(:)=0
call interpolate(state_handle, ens_size, location, QTY_PM10, PM10, pm10_istatus)
call track_status(ens_size, pm10_istatus, pm10, istatus, return_now)
if (return_now) return
where (istatus /= 0)
   ppm10 = missing_r8
elsewhere
   ppm10=ppm10+PM10 ! ug/kg
end where

alt_temp=1.0/(m_dry_air*p_Pa/(1000.0*8.314*T_k))
ppm10=ppm10/alt_temp ! pm10(ug/m3)


end subroutine get_expected_monitor_pm10


end module obs_def_monitor_mod
! END DART PREPROCESS MODULE CODE

