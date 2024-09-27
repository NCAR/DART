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
!
! BEGIN DART PREPROCESS KIND LIST
! IASI_CO_TOTAL_COL, QTY_CO
! END DART PREPROCESS KIND LIST
!
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_iasi_co_total_col_mod, only : get_expected_iasi_co_total_col, &
!                                  read_iasi_co_total_col, &
!                                  write_iasi_co_total_col, &
!                                  interactive_iasi_co_total_col
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(IASI_CO_TOTAL_COL)                                                           
!         call get_expected_iasi_co_total_col(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!
! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(IASI_CO_TOTAL_COL)
!         call read_iasi_co_total_col(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(IASI_CO_TOTAL_COL)
!         call write_iasi_co_total_col(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF
!
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(IASI_CO_TOTAL_COL)
!         call interactive_iasi_co_total_col(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_iasi_co_total_col_mod

use         apm_upper_bdy_mod, only :get_upper_bdy_fld, &
                                     get_MOZART_INT_DATA, &
                                     get_MOZART_REAL_DATA, &
                                     wrf_dart_ubval_interp, &
                                     apm_get_exo_coldens, &
                                     apm_get_upvals, &
                                     apm_interpolate

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

use          obs_kind_mod, only : QTY_CO, QTY_TEMPERATURE, QTY_PRESSURE, &
                                  QTY_VAPOR_MIXING_RATIO

use  ensemble_manager_mod, only : ensemble_type

use obs_def_utilities_mod, only : track_status

implicit none
private

public :: write_iasi_co_total_col, &
          read_iasi_co_total_col, &
          interactive_iasi_co_total_col, &
          get_expected_iasi_co_total_col, &
          set_obs_def_iasi_co_total_col

! Storage for the special information required for observations of this type
integer, parameter    :: max_iasi_co_obs = 10000000
integer               :: num_iasi_co_obs = 0
integer,  allocatable :: nlayer(:)
real(r8), allocatable :: avg_kernel(:,:)
real(r8), allocatable :: prior(:,:)

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_iasi_co_total_col_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

character(len=512) :: string1, string2
character(len=200) :: upper_data_file
character(len=200) :: upper_data_model
character(len=200) :: model
integer            :: ls_chem_dx, ls_chem_dy, ls_chem_dz, ls_chem_dt

logical, save :: module_initialized = .false.

! Namelist with default values
character(len=129)  :: IASI_CO_retrieval_type
logical :: use_log_co   = .false.
integer :: nlayer_model = -9999
integer :: nlayer_iasi = -9999
integer :: nlayer_iasi_co_total_col = -9999
integer :: nlayer_iasi_co_trop_col = -9999
integer :: nlayer_iasi_co_profile = -9999

   namelist /obs_def_IASI_CO_nml/upper_data_file, use_log_co, nlayer_model, &
   nlayer_iasi_co_total_col, nlayer_iasi_co_trop_col, nlayer_iasi_co_profile, &
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
call find_namelist_in_file("input.nml", "obs_def_IASI_CO_nml", iunit)
read(iunit, nml = obs_def_IASI_CO_nml, iostat = rc)
call check_namelist_read(iunit, rc, "obs_def_IASI_CO_nml")

! Record the namelist values
if (do_nml_file()) write(nmlfileunit, nml=obs_def_IASI_CO_nml)
if (do_nml_term()) write(     *     , nml=obs_def_IASI_CO_nml)
nlayer_iasi=nlayer_iasi_co_total_col

! Check for valid values

if (nlayer_model < 1) then
   write(string1,*)'obs_def_IASI_CO_nml:nlayer_model must be > 0, it is ',nlayer_model
   call error_handler(E_ERR,'initialize_module',string1,source)
endif

if (nlayer_iasi < 1) then
   write(string1,*)'obs_def_IASI_CO_nml:nlayer_iasi must be > 0, it is ',nlayer_iasi
   call error_handler(E_ERR,'initialize_module',string1,source)
endif

allocate(    nlayer(max_iasi_co_obs))
allocate(avg_kernel(max_iasi_co_obs,nlayer_iasi))
allocate(     prior(max_iasi_co_obs,nlayer_iasi))

end subroutine initialize_module

!-------------------------------------------------------------------------------

subroutine read_iasi_co_total_col(key, ifile, fform)

integer,          intent(out)          :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

! temporary arrays to hold buffer till we decide if we have enough room

integer               :: keyin
integer               :: nlayer_1
real(r8), allocatable :: avg_kernel_1(:)
real(r8), allocatable :: prior_1(:)
character(len=32)     :: fileformat

integer, SAVE :: counts1 = 0

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii" 
if(present(fform)) fileformat = adjustl(fform)

! Need to know how many layers for this one
nlayer_1 = read_int_scalar( ifile, fileformat, 'nlayer_1')

allocate(avg_kernel_1(nlayer_1))
allocate(     prior_1(nlayer_1))

call read_r8_array(ifile, nlayer_1,   avg_kernel_1, fileformat, 'avg_kernel_1')
call read_r8_array(ifile, nlayer_1,   prior_1,      fileformat, 'prior_1')
keyin = read_int_scalar(ifile, fileformat, 'nlayer_1')

counts1 = counts1 + 1
key     = counts1

if(counts1 > max_iasi_co_obs) then
   write(string1, *)'Not enough space for iasi co  col obs.'
   write(string2, *)'Can only have max_iasi_co_obs (currently ',max_iasi_co_obs,')'
   call error_handler(E_ERR,'read_iasi_co_tot_col',string1,source,text2=string2)
endif

call set_obs_def_iasi_co_total_col(key, avg_kernel_1,prior_1, nlayer_1)

deallocate(avg_kernel_1, prior_1)

end subroutine read_iasi_co_total_col

!-------------------------------------------------------------------------------

subroutine write_iasi_co_total_col(key, ifile, fform)

integer,          intent(in)           :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"
if(present(fform)) fileformat = adjustl(fform)

! nlayer, avg_kernel, and prior are all scoped in this module
! you can come extend the context strings to include the key if desired.

call write_int_scalar(ifile,                     nlayer(key), fileformat,'nlayer')
call write_r8_array(  ifile, nlayer(key),  avg_kernel(key,:), fileformat,'avg_kernel')
call write_r8_array(  ifile, nlayer(key),       prior(key,:), fileformat,'prior')
call write_int_scalar(ifile,                             key, fileformat,'key')

end subroutine write_iasi_co_total_col

!-------------------------------------------------------------------------------

subroutine interactive_iasi_co_total_col(key)

integer, intent(out) :: key

if ( .not. module_initialized ) call initialize_module

! STOP because routine is not finished.
write(string1,*)'interactive_iasi_co_total_col not yet working.'
call error_handler(E_ERR, 'interactive_iasi_co_total_col', string1, source)

if(num_iasi_co_obs >= max_iasi_co_obs) then
   write(string1, *)'Not enough space for an iasi co col obs.'
   write(string2, *)'Can only have max_iasi_co_obs (currently ',max_iasi_co_obs,')'
   call error_handler(E_ERR, 'interactive_iasi_co_total_col', string1, &
              source, text2=string2)
endif

! Increment the index
num_iasi_co_obs = num_iasi_co_obs + 1
key            = num_iasi_co_obs

! Otherwise, prompt for input for the three required beasts

write(*, *) 'Creating an interactive_iasi_co_total_col observation'
write(*, *) 'This featue is not setup '

end subroutine interactive_iasi_co_total_col

!-------------------------------------------------------------------------------

subroutine get_expected_iasi_co_total_col(state_handle, ens_size, location, key, expct_val, istatus)

type(ensemble_type), intent(in)  :: state_handle
type(location_type), intent(in)  :: location
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: key
integer,             intent(out) :: istatus(:)
real(r8),            intent(out) :: expct_val(:)

character(len=*), parameter :: routine = 'get_expected_iasi_co_total_col'
type(location_type) :: loc2

integer :: layer_iasi,level_iasi
integer :: layer_mdl,level_mdl
integer :: k,imem
integer :: co_istatus(ens_size)
integer, dimension(ens_size) :: tmp_istatus, qmr_istatus, prs_istatus

real(r8) :: eps, AvogN, Rd, Ru, grav, msq2cmsq
real(r8) :: co_min
real(r8) :: level
real(r8) :: tmp_vir_k, tmp_vir_kp
real(r8) :: mloc(3)
real(r8) :: co_val_conv
real(r8) :: up_wt,dw_wt,tl_wt,lnpr_mid
real(r8), dimension(ens_size) :: co_mdl_1, tmp_mdl_1, qmr_mdl_1, prs_mdl_1
real(r8), dimension(ens_size) :: co_mdl_n, tmp_mdl_n, qmr_mdl_n, prs_mdl_n
real(r8), dimension(ens_size) :: co_temp, tmp_temp, qmr_temp, prs_sfc

real(r8), allocatable, dimension(:)   :: thick, prs_iasi, prs_iasi_mem
real(r8), allocatable, dimension(:,:) :: co_val, tmp_val, qmr_val
logical  :: return_now

if ( .not. module_initialized ) call initialize_module

eps    =  0.61_r8
Rd     = 286.9_r8
Ru     = 8.316_r8
grav   =   9.8_r8
co_min = 1.e-6_r8
msq2cmsq = 1.e4_r8
AvogN = 6.02214e23_r8

if(use_log_co) then
   co_min = log(co_min)
endif

! Assign vertical grid information

layer_iasi = nlayer(key)
level_iasi = nlayer(key)+1
layer_mdl=nlayer_model
level_mdl=nlayer_model+1
allocate(prs_iasi(level_iasi))
allocate(prs_iasi_mem(level_iasi))

! for iasi need to figure pressure based on height
 prs_iasi(1)=1013.*100.
 prs_iasi(2)=1000.*100.
 prs_iasi(3)=950.*100.
 prs_iasi(4)=900.*100.
 prs_iasi(5)=850.*100.
 prs_iasi(6)=800.*100.
 prs_iasi(7)=750.*100.
 prs_iasi(8)=700.*100.
 prs_iasi(9)=650.*100.
 prs_iasi(10)=600.*100.
 prs_iasi(11)=550.*100.
 prs_iasi(12)=500.*100.
 prs_iasi(13)=450.*100.
 prs_iasi(14)=400.*100.
 prs_iasi(15)=350.*100.
 prs_iasi(16)=300.*100.
 prs_iasi(17)=250.*100.
 prs_iasi(18)=200.*100.
 prs_iasi(19)=150.*100.
 prs_iasi(20)=100.*100.
!
! Get location infomation

mloc = get_location(location)

if (    mloc(2) >  90.0_r8) then
        mloc(2) =  90.0_r8
elseif (mloc(2) < -90.0_r8) then
        mloc(2) = -90.0_r8
endif

! You could set a unique error code for each condition and then just return
! without having to issue a warning message. The error codes would then
! show up in the report from 'output_forward_op_errors'

istatus(:) = 0  ! set this once at the beginning
co_istatus = 0
tmp_istatus = 0
qmr_istatus = 0

! pressure at model surface (Pa)

level=1.0_r8
loc2 = set_location(mloc(1), mloc(2), level, VERTISSURFACE)
call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_sfc, prs_istatus) 
if(any(prs_istatus /= 0)) then
   write(string1, *)'APM NOTICE: MDL prs_sfc is bad '
   call error_handler(E_ERR, routine, string1, source)
endif
call track_status(ens_size, prs_istatus, prs_sfc, istatus, return_now)
if(return_now) return

! ozone at first model layer (ppmv)

level = 1.0_r8
loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
call interpolate(state_handle, ens_size, loc2, QTY_CO, co_mdl_1, co_istatus) 
if(any(co_istatus /= 0)) then
   write(string1, *)'APM NOTICE: MDL co_mdl_1 is bad '
   call error_handler(E_ERR, routine, string1, source)
endif
call track_status(ens_size, co_istatus, co_mdl_1, istatus, return_now)
if(return_now) return

! temperature at first model layer (K)

level=1.0_r8
loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_1, tmp_istatus) 
if(any(tmp_istatus /= 0)) then
   write(string1, *)'APM NOTICE: MDL tmp_mdl_1 is bad '
   call error_handler(E_ERR, routine, string1, source)
endif
call track_status(ens_size, tmp_istatus, tmp_mdl_1, istatus, return_now)
if(return_now) return

! vapor mixing ratio at first model layer (Kg/Kg)

level=1.0_r8
loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_1, qmr_istatus) 
if(any(qmr_istatus /= 0)) then
   write(string1, *)'APM NOTICE: MDL qmr_mdl_1 is bad '
   call error_handler(E_ERR, routine, string1, source)
endif
call track_status(ens_size, qmr_istatus, qmr_mdl_1, istatus, return_now)
if(return_now) return

! pressure at first model layer (Pa)

level=1.0_r8
loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_1, prs_istatus) 
if(any(prs_istatus /= 0)) then
   write(string1, *)'APM NOTICE: MDL prs_mdl_1 is bad '
   call error_handler(E_ERR, routine, string1, source)
endif
call track_status(ens_size, prs_istatus, prs_mdl_1, istatus, return_now)
if(return_now) return

! ozone at last model layer (ppmv)

level=real(layer_mdl-1)
loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
call interpolate(state_handle, ens_size, loc2, QTY_CO, co_mdl_n, co_istatus) 
if(any(co_istatus /= 0)) then
   write(string1, *)'APM NOTICE: MDL co_mdl_n is bad '
   call error_handler(E_ERR, routine, string1, source)
endif
call track_status(ens_size, co_istatus, co_mdl_n, istatus, return_now)
if(return_now) return

! temperature at last layer (K)

level=real(layer_mdl-1)
loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_n, tmp_istatus) 
if(any(tmp_istatus /= 0)) then
   write(string1, *)'APM NOTICE: MDL tmp_mdl_n is bad '
   call error_handler(E_ERR, routine, string1, source)
endif
call track_status(ens_size, tmp_istatus, tmp_mdl_n, istatus, return_now)
if(return_now) return

! vapor mixing ratio at last model layer (Kg/Kg)

level=real(layer_mdl-1)
loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_n, qmr_istatus) 
if(any(qmr_istatus /= 0)) then
   write(string1, *)'APM NOTICE: MDL qmr_mdl_n is bad '
   call error_handler(E_ERR, routine, string1, source)
endif
call track_status(ens_size, qmr_istatus, qmr_mdl_n, istatus, return_now)
if(return_now) return

! pressure at last model layer (Pa)

level=real(layer_mdl-1)
loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, &
prs_mdl_n, prs_istatus) 
if(any(prs_istatus /= 0)) then
   write(string1, *)'APM NOTICE: MDL prs_mdl_n is bad '
   call error_handler(E_ERR, routine, string1, source)
endif
call track_status(ens_size, prs_istatus, prs_mdl_n, istatus, return_now)
if(return_now) return

! Get profiles at IASI levels

allocate( co_val(ens_size,level_iasi))
allocate(tmp_val(ens_size,level_iasi))
allocate(qmr_val(ens_size,level_iasi))

do k=1,level_iasi
   co_istatus = 0
   tmp_istatus = 0
   qmr_istatus = 0

   loc2 = set_location(mloc(1), mloc(2), prs_iasi(k), VERTISPRESSURE)

   ! taking a different approach here ... interpolate all the required pieces
   ! for this level and then account for known special cases before determining
   ! if there is an error or not
   call interpolate(state_handle, ens_size, loc2, QTY_CO, co_temp, co_istatus)  
   call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_temp, tmp_istatus)  
   call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_temp, qmr_istatus)  

   ! Correcting for expected failures near the surface
   where(prs_iasi(k).ge.prs_mdl_1)
      co_istatus  = 0
      tmp_istatus = 0
      qmr_istatus = 0
      co_temp     = co_mdl_1
      tmp_temp    = tmp_mdl_1
      qmr_temp    = qmr_mdl_1
   endwhere

   ! Correcting for expected failures near the top
   where(prs_iasi(k).le.prs_mdl_n) 
      co_istatus  = 0
      tmp_istatus = 0
      qmr_istatus = 0
      co_temp     = co_mdl_n
      tmp_temp    = tmp_mdl_n
      qmr_temp    = qmr_mdl_n
   endwhere

   ! Report all issue before returning (when E_MSG is being used)
   if(any(co_istatus /= 0)) then
      write(string1,*) &
      'APM NOTICE: model CO obs value on IASI grid has a bad value '
      call error_handler(E_MSG, routine, string1, source)
      call track_status(ens_size, co_istatus, co_temp, istatus, return_now)
   endif
   
   if(any(tmp_istatus/=0)) then
      write(string1, *) &
      'APM NOTICE: model TMP obs value on IASI grid has a bad value'
      call error_handler(E_MSG, routine, string1, source)
      call track_status(ens_size, tmp_istatus, tmp_temp, istatus, return_now)
   endif
  
   if(any(tmp_istatus/=0)) then
      write(string1, *) &
      'APM NOTICE: model QMR obs value on IASI grid has a bad value '
      call error_handler(E_MSG, routine, string1, source)
      call track_status(ens_size, qmr_istatus, qmr_temp, istatus, return_now)
   endif
   if(return_now) return

   co_val(:,k) = co_temp(:)  
   tmp_val(:,k) = tmp_temp(:)  
   qmr_val(:,k) = qmr_temp(:)

   ! Convert units for co from ppmv
   co_val(:,k) = co_val(:,k) * 1.e-6_r8

enddo

! Adjust the IASI pressure for WRF-Chem lower/upper boudary pressure
! (IASI bottom to top)

expct_val(:)=0.0
allocate(thick(layer_iasi))
do imem=1,ens_size
   prs_iasi_mem(:)=prs_iasi(:)
   if (prs_sfc(imem).gt.prs_iasi_mem(1)) then
      prs_iasi_mem(1)=prs_sfc(imem)
   endif   

   ! Calculate the thicknesses

   thick(:)=0.
   do k=1,layer_iasi
      lnpr_mid=(log(prs_iasi_mem(k+1))+log(prs_iasi_mem(k)))/2.
      up_wt=log(prs_iasi_mem(k))-lnpr_mid
      dw_wt=log(lnpr_mid)-log(prs_iasi_mem(k+1))
      tl_wt=up_wt+dw_wt
      
      tmp_vir_k  = (1.0_r8 + eps*qmr_val(imem,k))*tmp_val(imem,k)
      tmp_vir_kp = (1.0_r8 + eps*qmr_val(imem,k+1))*tmp_val(imem,k+1)
      thick(k)   = Rd*(dw_wt*tmp_vir_k + up_wt*tmp_vir_kp)/tl_wt/grav* &
                   log(prs_iasi_mem(k)/prs_iasi_mem(k+1))
   enddo

   ! Process the vertical summation

   expct_val(imem)=0.0_r8
   do k=1,layer_iasi

   ! Convert from VMR to molar density (mol/m^3)
      if(use_log_co) then
         co_val_conv = (dw_wt*exp(co_val(imem,k))+up_wt*exp(co_val(imem,k+1)))/tl_wt * &
                        (dw_wt*prs_iasi_mem(k)+up_wt*prs_iasi_mem(k+1)) / &
                        (Ru*(dw_wt*tmp_val(imem,k)+up_wt*tmp_val(imem,k+1)))
      else
         co_val_conv = (dw_wt*co_val(imem,k)+up_wt*co_val(imem,k+1))/tl_wt * &
                        (dw_wt*prs_iasi_mem(k)+up_wt*prs_iasi_mem(k+1)) / &
                        (Ru*(dw_wt*tmp_val(imem,k)+up_wt*tmp_val(imem,k+1)))
      endif
 
   ! Get expected observation

      expct_val(imem) = expct_val(imem) + thick(k) * co_val_conv * &
                        avg_kernel(key,k) + (1.0_r8 - avg_kernel(key,k)) * &
                        prior(key,k)
   enddo
enddo

! Clean up and return
deallocate(co_val, tmp_val, qmr_val)
deallocate(thick)
deallocate(prs_iasi)

end subroutine get_expected_iasi_co_total_col

!-------------------------------------------------------------------------------

subroutine set_obs_def_iasi_co_total_col(key, co_avg_kernel, co_prior, co_nlayer)

integer,                           intent(in)   :: key, co_nlayer
real(r8), dimension(co_nlayer),    intent(in)   :: co_avg_kernel
real(r8), dimension(co_nlayer),    intent(in)   :: co_prior

if ( .not. module_initialized ) call initialize_module

if(num_iasi_co_obs >= max_iasi_co_obs) then
   write(string1, *)'Not enough space for iasi co col obs.'
   write(string2, *)'Can only have max_iasi_co_obs (currently ',max_iasi_co_obs,')'
   call error_handler(E_ERR,'set_obs_def_iasi_co_total_col',string1,source,revision, &
   revdate,text2=string2)
endif

nlayer(key) = co_nlayer
avg_kernel(key,:) = co_avg_kernel(:)
prior(key,:) = co_prior(:)

end subroutine set_obs_def_iasi_co_total_col

end module obs_def_iasi_co_total_col_mod

! END DART PREPROCESS MODULE CODE
