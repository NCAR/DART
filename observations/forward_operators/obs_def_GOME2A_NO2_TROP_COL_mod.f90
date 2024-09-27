 Copyright 2019 University Corporation for Atmospheric Research and 
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
!
! BEGIN DART PREPROCESS TYPE DEFINITIONS
! GOME2A_NO2_TROP_COL, QTY_NO2
! END DART PREPROCESS TYPE DEFINITIONS
!
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_gome2a_no2_trop_col_mod, only : get_expected_gome2a_no2_trop_col, &
!                                  read_gome2a_no2_trop_col, &
!                                  write_gome2a_no2_trop_col, &
!                                  interactive_gome2a_no2_trop_col
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(GOME2A_NO2_TROP_COL)                                                           
!         call get_expected_gome2a_no2_trop_col(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!
! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(GOME2A_NO2_TROP_COL)
!         call read_gome2a_no2_trop_col(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(GOME2A_NO2_TROP_COL)
!         call write_gome2a_no2_trop_col(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF
!
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(GOME2A_NO2_TROP_COL)
!         call interactive_gome2a_no2_trop_col(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_gome2a_no2_trop_col_mod

use             types_mod, only : r8, MISSING_R8

use         utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                                  nmlfileunit, check_namelist_read, &
                                  find_namelist_in_file, do_nml_file, do_nml_term, &
                                  ascii_file_format, &
                                  read_int_scalar, &
                                  write_int_scalar, &       
                                  read_r8_scalar, &
                                  write_r8_scalar, &
                                  read_r8_array, &
                                  write_r8_array

use          location_mod, only : location_type, set_location, get_location, &
                                  VERTISPRESSURE, VERTISSURFACE, VERTISLEVEL, &
                                  VERTISUNDEF

use       assim_model_mod, only : interpolate

use          obs_kind_mod, only : QTY_NO2, QTY_TEMPERATURE, QTY_SURFACE_PRESSURE, &
                                  QTY_PRESSURE, QTY_VAPOR_MIXING_RATIO

use  ensemble_manager_mod, only : ensemble_type

use obs_def_utilities_mod, only : track_status

implicit none
private

public :: write_gome2a_no2_trop_col, &
          read_gome2a_no2_trop_col, &
          interactive_gome2a_no2_trop_col, &
          get_expected_gome2a_no2_trop_col, &
          set_obs_def_gome2a_no2_trop_col

! Storage for the special information required for observations of this type
integer, parameter    :: max_gome2a_no2_obs = 10000000
integer               :: num_gome2a_no2_obs = 0
integer,  allocatable :: nlayer(:)
integer,  allocatable :: trop_indx(:)
real(r8), allocatable :: amf_trop(:)
real(r8), allocatable :: pressure(:,:)
real(r8), allocatable :: scat_wts(:,:)

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_gome2a_no2_trop_col_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

character(len=512) :: string1, string2
character(len=200) :: upper_data_file
character(len=200) :: upper_data_model
character(len=200) :: model
integer            :: ls_chem_dx, ls_chem_dy, ls_chem_dz, ls_chem_dt

logical, save :: module_initialized = .false.

! Namelist with default values
logical :: use_log_no2   = .false.
integer :: nlayer_model = -9999
integer :: nlayer_gome2a = -9999
integer :: nlayer_gome2a_no2_total_col = -9999
integer :: nlayer_gome2a_no2_trop_col = -9999

namelist /obs_def_GOME2A_NO2_nml/ upper_data_file, use_log_no2, nlayer_model, &
nlayer_gome2a_no2_total_col, nlayer_gome2a_no2_trop_col, &
ls_chem_dx, ls_chem_dy, ls_chem_dz, ls_chem_dt, upper_data_model

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

subroutine initialize_module

integer :: iunit, rc

! Prevent multiple calls from executing this code more than once.

if (module_initialized) return

call register_module(source, revision, revdate)
module_initialized = .true.

! Read namelist values
call find_namelist_in_file("input.nml", "obs_def_GOME2A_NO2_nml", iunit)
read(iunit, nml = obs_def_GOME2A_NO2_nml, iostat = rc)
call check_namelist_read(iunit, rc, "obs_def_GOME2A_NO2_nml")

! Record the namelist values
if (do_nml_file()) write(nmlfileunit, nml=obs_def_GOME2A_NO2_nml)
if (do_nml_term()) write(     *     , nml=obs_def_GOME2A_NO2_nml)
nlayer_gome2a=nlayer_gome2a_no2_trop_col

! Check for valid values

if (nlayer_model < 1) then
   write(string1,*)'obs_def_GOME2A_NO2_nml:nlayer_model must be > 0, it is ',nlayer_model
   call error_handler(E_ERR,'initialize_module',string1,source)
endif

if (nlayer_gome2a < 1) then
   write(string1,*)'obs_def_GOME2A_NO2_nml:nlayer_gome2a must be > 0, it is ',nlayer_gome2a
   call error_handler(E_ERR,'initialize_module',string1,source)
endif

allocate(    nlayer(max_gome2a_no2_obs))
allocate(    trop_indx(max_gome2a_no2_obs))
allocate(    amf_trop(max_gome2a_no2_obs))
allocate(  pressure(max_gome2a_no2_obs,nlayer_gome2a+1))
allocate(scat_wts(max_gome2a_no2_obs,nlayer_gome2a))

end subroutine initialize_module

!-------------------------------------------------------------------------------

subroutine read_gome2a_no2_trop_col(key, ifile, fform)

integer,          intent(out)          :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

! temporary arrays to hold buffer till we decide if we have enough room

integer               :: keyin
integer               :: nlayer_1
integer               :: trop_indx_1
real(r8)              :: amf_trop_1
real(r8), allocatable :: pressure_1(:)
real(r8), allocatable :: scat_wts_1(:)
character(len=32)     :: fileformat

integer, SAVE :: counts1 = 0

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii" 
if(present(fform)) fileformat = adjustl(fform)

! Need to know how many layers for this one
nlayer_1 = read_int_scalar( ifile, fileformat, 'nlayer_1')
trop_indx_1 = read_int_scalar( ifile, fileformat, 'trop_indx_1')
amf_trop_1 = read_r8_scalar( ifile, fileformat, 'amf_trop_1')

allocate(  pressure_1(nlayer_1+1))
allocate(scat_wts_1(nlayer_1))

call read_r8_array(ifile, nlayer_1+1, pressure_1,   fileformat, 'pressure_1')
call read_r8_array(ifile, nlayer_1,   scat_wts_1, fileformat, 'scat_wts_1')
keyin = read_int_scalar(ifile, fileformat, 'keyin')

counts1 = counts1 + 1
key     = counts1

if(counts1 > max_gome2a_no2_obs) then
   write(string1, *)'Not enough space for gome2a no2 trop col obs.'
   write(string2, *)'Can only have max_gome2a_no2_obs (currently ',max_gome2a_no2_obs,')'
   call error_handler(E_ERR,'read_gome2a_no2_trop_col',string1,source,text2=string2)
endif

call set_obs_def_gome2a_no2_trop_col(key, pressure_1, scat_wts_1, amf_trop_1, trop_indx_1, nlayer_1)

deallocate(pressure_1, scat_wts_1)

end subroutine read_gome2a_no2_trop_col

!-------------------------------------------------------------------------------

subroutine write_gome2a_no2_trop_col(key, ifile, fform)

integer,          intent(in)           :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"
if(present(fform)) fileformat = adjustl(fform)

! nlayer, amf_trop, trop_indx, pressure, and scatering weights are all scoped in this module
! you can come extend the context strings to include the key if desired.

call write_int_scalar(ifile,                     nlayer(key), fileformat,'nlayer')
call write_int_scalar(ifile,                     trop_indx(key), fileformat,'trop_indx')
call write_r8_scalar( ifile,                     amf_trop(key), fileformat,'amf_trop')
call write_r8_array(  ifile, nlayer(key)+1,      pressure(key,:), fileformat,'pressure')
call write_r8_array(  ifile, nlayer(key),        scat_wts(key,:), fileformat,'scat_wts')
call write_int_scalar(ifile,                     key, fileformat,'key')

end subroutine write_gome2a_no2_trop_col

!-------------------------------------------------------------------------------

subroutine interactive_gome2a_no2_trop_col(key)

integer, intent(out) :: key

if ( .not. module_initialized ) call initialize_module

! STOP because routine is not finished.
write(string1,*)'interactive_gome2a_no2_trop_col not yet working.'
call error_handler(E_ERR, 'interactive_gome2a_no2_trop_col', string1, source)

if(num_gome2a_no2_obs >= max_gome2a_no2_obs) then
   write(string1, *)'Not enough space for an gome2a no2 trop col obs.'
   write(string2, *)'Can only have max_gome2a_no2_obs (currently ',max_gome2a_no2_obs,')'
   call error_handler(E_ERR, 'interactive_gome2a_no2_trop_col', string1, &
              source, text2=string2)
endif

! Increment the index
num_gome2a_no2_obs = num_gome2a_no2_obs + 1
key            = num_gome2a_no2_obs

! Otherwise, prompt for input for the three required beasts

write(*, *) 'Creating an interactive_gome2a_no2_trop_col observation'
write(*, *) 'This featue is not setup '

end subroutine interactive_gome2a_no2_trop_col

!-------------------------------------------------------------------------------

subroutine get_expected_gome2a_no2_trop_col(state_handle, ens_size, location, key, expct_val, istatus)

   type(ensemble_type), intent(in)  :: state_handle
   type(location_type), intent(in)  :: location
   integer,             intent(in)  :: ens_size
   integer,             intent(in)  :: key
   integer,             intent(out) :: istatus(:)
   real(r8),            intent(out) :: expct_val(:)

   character(len=*), parameter :: routine = 'get_expected_gome2a_no2_trop_col'
   type(location_type) :: loc2

   integer :: layer_gome2a,level_gome2a
   integer :: layer_mdl,level_mdl
   integer :: k,imem,kend_gome2a
   integer :: interp_new
   integer, dimension(ens_size) :: zstatus

   real(r8) :: eps, AvogN, Rd, Ru, Cp, grav, msq2cmsq
   real(r8) :: missing,no2_min,tmp_max
   real(r8) :: level,del_prs
   real(r8) :: tmp_vir_k, tmp_vir_kp
   real(r8) :: mloc(3)
   real(r8) :: no2_val_conv
   real(r8) :: up_wt,dw_wt,tl_wt,lnpr_mid
   real(r8), dimension(ens_size) :: no2_mdl_1, tmp_mdl_1, qmr_mdl_1, prs_mdl_1
   real(r8), dimension(ens_size) :: no2_mdl_1p, tmp_mdl_1p, qmr_mdl_1p, prs_mdl_1p
   real(r8), dimension(ens_size) :: no2_mdl_n, tmp_mdl_n, qmr_mdl_n, prs_mdl_n
   real(r8), dimension(ens_size) :: no2_mdl_nm, tmp_mdl_nm, qmr_mdl_nm, prs_mdl_nm
   real(r8), dimension(ens_size) :: no2_temp, tmp_temp, qmr_temp, prs_sfc

   real(r8), allocatable, dimension(:)   :: thick, prs_gome2a, prs_gome2a_mem
   real(r8), allocatable, dimension(:,:) :: no2_val, tmp_val, qmr_val
   logical  :: return_now,no2_return_now,tmp_return_now,qmr_return_now

   if ( .not. module_initialized ) call initialize_module

   eps      =  0.61_r8
   Rd       = 287.05_r8     ! J/kg
   Ru       = 8.316_r8      ! J/kg
   Cp       = 1006.0        ! J/kg/K
   grav     =   9.8_r8
   no2_min  = 1.e-6_r8
   msq2cmsq = 1.e4_r8
   AvogN    = 6.02214e23_r8
   missing  = -888888_r8
   tmp_max  = 600.
   del_prs  = 5000.

   if(use_log_no2) then
      no2_min = log(no2_min)
   endif

! Assign vertical grid information

   layer_gome2a = nlayer(key)
   level_gome2a = nlayer(key)+1
   kend_gome2a  = trop_indx(key)
   layer_mdl = nlayer_model
   level_mdl = nlayer_model+1

   allocate(prs_gome2a(level_gome2a))
   allocate(prs_gome2a_mem(level_gome2a))
   prs_gome2a(1:level_gome2a)=pressure(key,1:level_gome2a)

! Get location infomation

   mloc = get_location(location)

   if (mloc(2) >  90.0_r8) then
      mloc(2) =  90.0_r8
   elseif (mloc(2) < -90.0_r8) then
      mloc(2) = -90.0_r8
   endif

! You could set a unique error code for each condition and then just return
! without having to issue a warning message. The error codes would then
! show up in the report from 'output_forward_op_errors'

   istatus(:) = 0  ! set this once at the beginning
   return_now=.false.

! pressure at model surface (Pa)

   zstatus=0
   level=0.0_r8
   loc2 = set_location(mloc(1), mloc(2), level, VERTISSURFACE)
   call interpolate(state_handle, ens_size, loc2, QTY_SURFACE_PRESSURE, prs_sfc, zstatus) 

   no2_mdl_1(:)=missing_r8
   tmp_mdl_1(:)=missing_r8
   qmr_mdl_1(:)=missing_r8
   prs_mdl_1(:)=missing_r8

   do k=1,layer_mdl
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_NO2, no2_mdl_1, zstatus) 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_1, &
      zstatus) 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_1, &
      zstatus) 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_1, &
      zstatus) 
!
      interp_new=0
      do imem=1,ens_size
         if(no2_mdl_1(imem).eq.missing_r8 .or. tmp_mdl_1(imem).eq.missing_r8 .or. &
         qmr_mdl_1(imem).eq.missing_r8 .or. prs_mdl_1(imem).eq.missing_r8) then
            interp_new=1
            exit
         endif
      enddo
      if(interp_new.eq.0) then
         exit
      endif    
   enddo

!   write(string1, *)'APM: no2 lower bound ',no2_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: tmp lower bound ',tmp_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: qmr lower bound ',qmr_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: prs lower bound ',prs_mdl_1
!   call error_handler(E_MSG, routine, string1, source)

   no2_mdl_n(:)=missing_r8
   tmp_mdl_n(:)=missing_r8
   qmr_mdl_n(:)=missing_r8
   prs_mdl_n(:)=missing_r8

   do k=layer_mdl,1,-1
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
      call interpolate(state_handle, ens_size, loc2, QTY_NO2, no2_mdl_n, zstatus) ! ppmv
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_n, zstatus) ! K 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_n, zstatus) ! kg / kg 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_n, zstatus) ! Pa
!
      interp_new=0
      do imem=1,ens_size
         if(no2_mdl_n(imem).eq.missing_r8 .or. tmp_mdl_n(imem).eq.missing_r8 .or. &
         qmr_mdl_n(imem).eq.missing_r8 .or. prs_mdl_n(imem).eq.missing_r8) then
            interp_new=1
            exit
         endif
      enddo
      if(interp_new.eq.0) then
         exit
      endif    
   enddo
   
!   write(string1, *)'APM: no2 upper bound ',no2_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: tmp upper bound ',tmp_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: qmr upper bound ',qmr_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: prs upper bound ',prs_mdl_n
!   call error_handler(E_MSG, routine, string1, source)

! Get profiles at GOME2A pressure levels

   allocate(no2_val(ens_size,level_gome2a))
   allocate(tmp_val(ens_size,level_gome2a))
   allocate(qmr_val(ens_size,level_gome2a))

   do k=1,level_gome2a
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), prs_gome2a(k), VERTISPRESSURE)
      call interpolate(state_handle, ens_size, loc2, QTY_NO2, no2_val(:,k), zstatus)  
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_val(:,k), zstatus)  
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_val(:,k), zstatus)  
!
! Correcting for expected failures near the surface
      do imem=1,ens_size
         if (prs_gome2a(k).ge.prs_mdl_1(imem)) then
            no2_val(imem,k) = no2_mdl_1(imem)
            tmp_val(imem,k) = tmp_mdl_1(imem)
            qmr_val(imem,k) = qmr_mdl_1(imem)
            cycle
         endif

! Correcting for expected failures near the top
         if (prs_gome2a(k).le.prs_mdl_n(imem)) then
            no2_val(imem,k) = no2_mdl_n(imem)
            tmp_val(imem,k) = tmp_mdl_n(imem)
            qmr_val(imem,k) = qmr_mdl_n(imem)
            cycle
         endif
      enddo
!
!      write(string1, *)'APM: no2 ',k,no2_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: tmp ',k,tmp_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: qmr ',k,qmr_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!
! Check data for missing values      
      do imem=1,ens_size
         if(no2_val(imem,k).eq.missing_r8 .or. tmp_val(imem,k).eq.missing_r8 .or. &
         qmr_val(imem,k).eq.missing_r8) then
            zstatus(:)=20
            expct_val(:)=missing_r8
            write(string1, *) &
            'APM: Input data has missing values '
            call error_handler(E_MSG, routine, string1, source)
            call track_status(ens_size, zstatus, expct_val, istatus, return_now)
            return
         endif
      enddo
!      
! Convert units for no2 from ppmv
      no2_val(:,k) = no2_val(:,k) * 1.e-6_r8
   enddo

   istatus(:)=0.
   zstatus(:)=0.
   expct_val(:)=0.0
   allocate(thick(layer_gome2a))
!
   do imem=1,ens_size
! Adjust the GOME2A pressure for WRF-Chem lower/upper boudary pressure
! (GOME2A NO2 vertical grid is bottom to top)
      prs_gome2a_mem(:)=prs_gome2a(:)
      if (prs_sfc(imem).gt.prs_gome2a_mem(1)) then
         prs_gome2a_mem(1)=prs_sfc(imem)
      endif   

! Calculate the thicknesses

      do k=1,kend_gome2a
         lnpr_mid=(log(prs_gome2a_mem(k+1))+log(prs_gome2a_mem(k)))/2.
         up_wt=log(prs_gome2a_mem(k))-lnpr_mid
         dw_wt=lnpr_mid-log(prs_gome2a_mem(k+1))
         tl_wt=up_wt+dw_wt
      
         tmp_vir_k  = (1.0_r8 + eps*qmr_val(imem,k))*tmp_val(imem,k)
         tmp_vir_kp = (1.0_r8 + eps*qmr_val(imem,k+1))*tmp_val(imem,k+1)
         thick(k)   = Rd*(dw_wt*tmp_vir_k + up_wt*tmp_vir_kp)/tl_wt/grav* &
                   log(prs_gome2a_mem(k)/prs_gome2a_mem(k+1))
      enddo

! Process the vertical summation

      do k=1,kend_gome2a
         lnpr_mid=(log(prs_gome2a_mem(k+1))+log(prs_gome2a_mem(k)))/2.
         up_wt=log(prs_gome2a_mem(k))-lnpr_mid
         dw_wt=lnpr_mid-log(prs_gome2a_mem(k+1))
         tl_wt=up_wt+dw_wt

! Convert from VMR to molar density (mol/m^3)
         if(use_log_no2) then
            no2_val_conv = (dw_wt*exp(no2_val(imem,k))+up_wt*exp(no2_val(imem,k+1)))/tl_wt * &
            (dw_wt*prs_gome2a_mem(k)+up_wt*prs_gome2a_mem(k+1)) / &
            (Ru*(dw_wt*tmp_val(imem,k)+up_wt*tmp_val(imem,k+1)))
         else
            no2_val_conv = (dw_wt*no2_val(imem,k)+up_wt*no2_val(imem,k+1))/tl_wt * &
            (dw_wt*prs_gome2a_mem(k)+up_wt*prs_gome2a_mem(k+1)) / &
            (Ru*(dw_wt*tmp_val(imem,k)+up_wt*tmp_val(imem,k+1)))
         endif
 
! Get expected observation

         expct_val(imem) = expct_val(imem) + thick(k) * no2_val_conv * &
         AvogN / msq2cmsq * scat_wts(key,k)
!         write(string1, *) 'APM: key,k,imem,expc_val,thick,no2_conv,scwt ', &
!         key,k,imem,expct_val(imem),thick(k),no2_val_conv,scat_wts(key,k)
!         call error_handler(E_MSG, routine, string1, source)
      enddo

      if(expct_val(imem).lt.0) then
         zstatus(imem)=20
         expct_val(:)=missing_r8
         write(string1, *) &
         'APM NOTICE: TROPOMI NO2 expected value is negative '
         call error_handler(E_MSG, routine, string1, source)
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
   enddo

! Clean up and return
   deallocate(no2_val, tmp_val, qmr_val)
   deallocate(thick)
   deallocate(prs_gome2a, prs_gome2a_mem)

end subroutine get_expected_gome2a_no2_trop_col

!-------------------------------------------------------------------------------

subroutine set_obs_def_gome2a_no2_trop_col(key, no2_pressure, no2_scat_wts, no2_amf_trop, no2_trop_indx, no2_nlayer)

integer,                           intent(in)   :: key, no2_nlayer
integer,                           intent(in)   :: no2_trop_indx
real(r8),                          intent(in)   :: no2_amf_trop
real(r8), dimension(no2_nlayer+1),  intent(in)   :: no2_pressure
real(r8), dimension(no2_nlayer),    intent(in)   :: no2_scat_wts

if ( .not. module_initialized ) call initialize_module

if(num_gome2a_no2_obs >= max_gome2a_no2_obs) then
   write(string1, *)'Not enough space for gome2a no2 trop col obs.'
   write(string2, *)'Can only have max_gome2a_no2_obs (currently ',max_gome2a_no2_obs,')'
   call error_handler(E_ERR,'set_obs_def_gome2a_no2_trop_col',string1,source,revision, &
   revdate,text2=string2)
endif

nlayer(key) = no2_nlayer
trop_indx(key) = no2_trop_indx
amf_trop(key) = no2_amf_trop
pressure(key,1:no2_nlayer+1) = no2_pressure(1:no2_nlayer+1)
scat_wts(key,1:no2_nlayer) = no2_scat_wts(1:no2_nlayer)

end subroutine set_obs_def_gome2a_no2_trop_col

end module obs_def_gome2a_no2_trop_col_mod

! END DART PREPROCESS MODULE CODE
