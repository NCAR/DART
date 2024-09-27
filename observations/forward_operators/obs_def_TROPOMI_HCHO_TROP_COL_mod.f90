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
! BEGIN DART PREPROCESS TYPE DEFINITIONS
! TROPOMI_HCHO_TROP_COL, QTY_HCHO
! END DART PREPROCESS TYPE DEFINITIONS
!
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_tropomi_hcho_trop_col_mod, only : get_expected_tropomi_hcho_trop_col, &
!                                  read_tropomi_hcho_trop_col, &
!                                  write_tropomi_hcho_trop_col, &
!                                  interactive_tropomi_hcho_trop_col
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(TROPOMI_HCHO_TROP_COL)                                                           
!         call get_expected_tropomi_hcho_trop_col(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!
! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(TROPOMI_HCHO_TROP_COL)
!         call read_tropomi_hcho_trop_col(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(TROPOMI_HCHO_TROP_COL)
!         call write_tropomi_hcho_trop_col(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF
!
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(TROPOMI_HCHO_TROP_COL)
!         call interactive_tropomi_hcho_trop_col(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_tropomi_hcho_trop_col_mod

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

use          obs_kind_mod, only : QTY_HCHO, QTY_TEMPERATURE, QTY_SURFACE_PRESSURE, &
                                  QTY_PRESSURE, QTY_VAPOR_MIXING_RATIO

use  ensemble_manager_mod, only : ensemble_type

use obs_def_utilities_mod, only : track_status

implicit none
private

public :: write_tropomi_hcho_trop_col, &
          read_tropomi_hcho_trop_col, &
          interactive_tropomi_hcho_trop_col, &
          get_expected_tropomi_hcho_trop_col, &
          set_obs_def_tropomi_hcho_trop_col

! Storage for the special information required for observations of this type
integer, parameter    :: max_tropomi_hcho_obs = 10000000
integer               :: num_tropomi_hcho_obs = 0
integer,  allocatable :: nlayer(:)
integer,  allocatable :: trop_indx(:)
real(r8), allocatable :: amf_trop(:)
real(r8), allocatable :: pressure(:,:)
real(r8), allocatable :: avg_kernel(:,:)

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_tropomi_hcho_trop_col_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

character(len=512) :: string1, string2
character(len=200) :: upper_data_file
character(len=200) :: upper_data_model
character(len=200) :: model
integer            :: ls_chem_dx, ls_chem_dy, ls_chem_dz, ls_chem_dt

logical, save :: module_initialized = .false.

! Namelist with default values
logical :: use_log_hcho   = .false.
integer :: nlayer_model = -9999
integer :: nlayer_tropomi = -9999
integer :: nlayer_tropomi_hcho_total_col = -9999
integer :: nlayer_tropomi_hcho_trop_col = -9999

namelist /obs_def_TROPOMI_HCHO_nml/ upper_data_file, use_log_hcho, nlayer_model, &
nlayer_tropomi_hcho_total_col, nlayer_tropomi_hcho_trop_col, &
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
call find_namelist_in_file("input.nml", "obs_def_TROPOMI_HCHO_nml", iunit)
read(iunit, nml = obs_def_TROPOMI_HCHO_nml, iostat = rc)
call check_namelist_read(iunit, rc, "obs_def_TROPOMI_HCHO_nml")

! Record the namelist values
if (do_nml_file()) write(nmlfileunit, nml=obs_def_TROPOMI_HCHO_nml)
if (do_nml_term()) write(     *     , nml=obs_def_TROPOMI_HCHO_nml)
nlayer_tropomi=nlayer_tropomi_hcho_trop_col

! Check for valid values

if (nlayer_model < 1) then
   write(string1,*)'obs_def_TROPOMI_HCHO_nml:nlayer_model must be > 0, it is ',nlayer_model
   call error_handler(E_ERR,'initialize_module',string1,source)
endif

if (nlayer_tropomi < 1) then
   write(string1,*)'obs_def_TROPOMI_HCHO_nml:nlayer_tropomi must be > 0, it is ',nlayer_tropomi
   call error_handler(E_ERR,'initialize_module',string1,source)
endif

allocate(    nlayer(max_tropomi_hcho_obs))
allocate(    trop_indx(max_tropomi_hcho_obs))
allocate(    amf_trop(max_tropomi_hcho_obs))
allocate(  pressure(max_tropomi_hcho_obs,nlayer_tropomi+1))
allocate(avg_kernel(max_tropomi_hcho_obs,nlayer_tropomi))

end subroutine initialize_module

!-------------------------------------------------------------------------------

subroutine read_tropomi_hcho_trop_col(key, ifile, fform)

integer,          intent(out)          :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

! temporary arrays to hold buffer till we decide if we have enough room

integer               :: keyin
integer               :: nlayer_1
integer               :: trop_indx_1
real(r8)              :: amf_trop_1
real(r8), allocatable :: pressure_1(:)
real(r8), allocatable :: avg_kernel_1(:)
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
allocate(avg_kernel_1(nlayer_1))

call read_r8_array(ifile, nlayer_1+1, pressure_1,   fileformat, 'pressure_1')
call read_r8_array(ifile, nlayer_1,   avg_kernel_1, fileformat, 'avg_kernel_1')
keyin = read_int_scalar(ifile, fileformat, 'keyin')

counts1 = counts1 + 1
key     = counts1

if(counts1 > max_tropomi_hcho_obs) then
   write(string1, *)'Not enough space for tropomi hcho total col obs.'
   write(string2, *)'Can only have max_tropomi_hcho_obs (currently ',max_tropomi_hcho_obs,')'
   call error_handler(E_ERR,'read_tropomi_hcho_trop_col',string1,source,text2=string2)
endif

call set_obs_def_tropomi_hcho_trop_col(key, pressure_1, avg_kernel_1, amf_trop_1, trop_indx_1, nlayer_1)

deallocate(pressure_1, avg_kernel_1)

end subroutine read_tropomi_hcho_trop_col

!-------------------------------------------------------------------------------

subroutine write_tropomi_hcho_trop_col(key, ifile, fform)

integer,          intent(in)           :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"
if(present(fform)) fileformat = adjustl(fform)

! nlayer, amf_trop, trop_indx, pressure, and avg_kernel are all scoped in this module
! you can come extend the context strings to include the key if desired.

call write_int_scalar(ifile,                     nlayer(key), fileformat,'nlayer')
call write_int_scalar(ifile,                     trop_indx(key), fileformat,'trop_indx')
call write_r8_scalar( ifile,                     amf_trop(key), fileformat,'amf_trop')
call write_r8_array(  ifile, nlayer(key)+1,      pressure(key,:), fileformat,'pressure')
call write_r8_array(  ifile, nlayer(key),        avg_kernel(key,:), fileformat,'avg_kernel')
call write_int_scalar(ifile,                     key, fileformat,'key')

end subroutine write_tropomi_hcho_trop_col

!-------------------------------------------------------------------------------

subroutine interactive_tropomi_hcho_trop_col(key)

integer, intent(out) :: key

if ( .not. module_initialized ) call initialize_module

! STOP because routine is not finished.
write(string1,*)'interactive_tropomi_hcho_trop_col not yet working.'
call error_handler(E_ERR, 'interactive_tropomi_hcho_trop_col', string1, source)

if(num_tropomi_hcho_obs >= max_tropomi_hcho_obs) then
   write(string1, *)'Not enough space for an tropomi hcho total colobs.'
   write(string2, *)'Can only have max_tropomi_hcho_obs (currently ',max_tropomi_hcho_obs,')'
   call error_handler(E_ERR, 'interactive_tropomi_hcho_trop_col', string1, &
              source, text2=string2)
endif

! Increment the index
num_tropomi_hcho_obs = num_tropomi_hcho_obs + 1
key            = num_tropomi_hcho_obs

! Otherwise, prompt for input for the three required beasts

write(*, *) 'Creating an interactive_tropomi_hcho_trop_col observation'
write(*, *) 'This featue is not setup '

end subroutine interactive_tropomi_hcho_trop_col

!-------------------------------------------------------------------------------

subroutine get_expected_tropomi_hcho_trop_col(state_handle, ens_size, location, key, expct_val, istatus)

   type(ensemble_type), intent(in)  :: state_handle
   type(location_type), intent(in)  :: location
   integer,             intent(in)  :: ens_size
   integer,             intent(in)  :: key
   integer,             intent(out) :: istatus(:)
   real(r8),            intent(out) :: expct_val(:)
   
   character(len=*), parameter :: routine = 'get_expected_tropomi_hcho_trop_col'
   type(location_type) :: loc2
   
   integer :: layer_tropomi,level_tropomi
   integer :: layer_mdl,level_mdl
   integer :: k,imem,kend_tropomi
   integer :: interp_new
   integer, dimension(ens_size) :: zstatus
   
   real(r8) :: eps, AvogN, Rd, Ru, Cp, grav, msq2cmsq
   real(r8) :: missing,hcho_min, tmp_max
   real(r8) :: level,del_prs
   real(r8) :: tmp_vir_k, tmp_vir_kp
   real(r8) :: mloc(3)
   real(r8) :: hcho_val_conv,amf
   real(r8) :: up_wt,dw_wt,tl_wt,lnpr_mid
   real(r8), dimension(ens_size) :: hcho_mdl_1, tmp_mdl_1, qmr_mdl_1, prs_mdl_1
   real(r8), dimension(ens_size) :: hcho_mdl_n, tmp_mdl_n, qmr_mdl_n, prs_mdl_n
   real(r8), dimension(ens_size) :: prs_sfc
   
   real(r8), allocatable, dimension(:)   :: thick, prs_tropomi, prs_tropomi_mem
   real(r8), allocatable, dimension(:,:) :: hcho_val, tmp_val, qmr_val
   logical  :: return_now,hcho_return_now,tmp_return_now,qmr_return_now
   
   if ( .not. module_initialized ) call initialize_module
   
   eps      =  0.61_r8
   Rd       = 287.05_r8     ! J/kg
   Ru       = 8.316_r8      ! J/kg
   Cp       = 1006.0        ! J/kg/K
   grav     =   9.8_r8
   hcho_min  = 1.e-6_r8
   msq2cmsq = 1.e4_r8
   AvogN    = 6.02214e23_r8
   missing  = -888888_r8
   tmp_max  = 600.
   del_prs  = 5000.
   
   if(use_log_hcho) then
      hcho_min = log(hcho_min)
   endif

! Assign vertical grid information

   layer_tropomi = nlayer(key)
   level_tropomi = nlayer(key)+1
   kend_tropomi  = trop_indx(key)
   layer_mdl = nlayer_model
   level_mdl = nlayer_model+1
   amf       = amf_trop(key)
   
   allocate(prs_tropomi(level_tropomi))
   allocate(prs_tropomi_mem(level_tropomi))
   prs_tropomi(1:layer_tropomi)=pressure(key,1:layer_tropomi+1)

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

   hcho_mdl_1(:)=missing_r8
   tmp_mdl_1(:)=missing_r8
   qmr_mdl_1(:)=missing_r8
   prs_mdl_1(:)=missing_r8

   do k=1,layer_mdl
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
      call interpolate(state_handle, ens_size, loc2, QTY_HCHO, hcho_mdl_1, zstatus) ! ppmv 
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_1, zstatus) ! K 
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_1, zstatus) ! kg / kg 
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_1, zstatus) ! Pa
!
      interp_new=0
      do imem=1,ens_size
         if(hcho_mdl_1(imem).eq.missing_r8 .or. tmp_mdl_1(imem).eq.missing_r8 .or. &
         qmr_mdl_1(imem).eq.missing_r8 .or. prs_mdl_1(imem).eq.missing_r8) then
            interp_new=1
            exit
         endif
      enddo
      if(interp_new.eq.0) then
         exit
      endif    
   enddo

!   write(string1, *)'APM: hcho lower bound 1 ',hcho_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: tmp lower bound 1 ',tmp_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: qmr lower bound 1 ',qmr_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: prs lower bound 1 ',prs_mdl_1
!   call error_handler(E_MSG, routine, string1, source)

   hcho_mdl_n(:)=missing_r8
   tmp_mdl_n(:)=missing_r8
   qmr_mdl_n(:)=missing_r8
   prs_mdl_n(:)=missing_r8

   do k=layer_mdl,1,-1
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)  
      call interpolate(state_handle, ens_size, loc2, QTY_HCHO, hcho_mdl_n, zstatus) ! ppmv
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_n, zstatus) ! K 
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_n, zstatus) ! kg / kg 
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_n, zstatus) ! Pa
!
      interp_new=0
      do imem=1,ens_size
         if(hcho_mdl_n(imem).eq.missing_r8 .or. tmp_mdl_n(imem).eq.missing_r8 .or. &
         qmr_mdl_n(imem).eq.missing_r8 .or. prs_mdl_n(imem).eq.missing_r8) then
            interp_new=1
            exit
         endif
      enddo
      if(interp_new.eq.0) then
         exit
      endif    
   enddo

!   write(string1, *)'APM: hcho upper bound 1 ',hcho_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: tmp upper bound 1 ',tmp_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: qmr upper bound 1 ',qmr_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: prs upper bound 1 ',prs_mdl_n
!   call error_handler(E_MSG, routine, string1, source)

! Get profiles at TROPOMI pressure levels

   allocate(hcho_val(ens_size,level_tropomi))
   allocate(tmp_val(ens_size,level_tropomi))
   allocate(qmr_val(ens_size,level_tropomi))

   do k=1,level_tropomi
      zstatus=0
      loc2 = set_location(mloc(1), mloc(2), prs_tropomi(k), VERTISPRESSURE)
      call interpolate(state_handle, ens_size, loc2, QTY_HCHO, hcho_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_val(:,k), zstatus)  
!
! Correcting for expected failures near the surface
      do imem=1,ens_size
         if (prs_tropomi(k).ge.prs_mdl_1(imem)) then
            hcho_val(imem,k) = hcho_mdl_1(imem)
            tmp_val(imem,k) = tmp_mdl_1(imem)
            qmr_val(imem,k) = qmr_mdl_1(imem)
            cycle
         endif

! Correcting for expected failures near the top
         if (prs_tropomi(k).le.prs_mdl_n(imem)) then
            hcho_val(imem,k) = hcho_mdl_n(imem)
            tmp_val(imem,k) = tmp_mdl_n(imem)
            qmr_val(imem,k) = qmr_mdl_n(imem)
            cycle
         endif
      enddo
!
!      write(string1, *)'APM: hcho ',k,hcho_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: tmp ',k,tmp_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: qmr ',k,qmr_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!
! Check data for missing values      
      do imem=1,ens_size
         if(hcho_val(imem,k).eq.missing_r8 .or. tmp_val(imem,k).eq.missing_r8 .or. &
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
! Convert units for hcho from ppmv
      hcho_val(:,k) = hcho_val(:,k) * 1.e-6_r8
   enddo

   istatus=0
   zstatus(:)=0
   expct_val(:)=0.0
   allocate(thick(layer_tropomi))

   do imem=1,ens_size
! Adjust the TROPOMI pressure for WRF-Chem lower/upper boudary pressure
! (TROPOMI HCHO vertical grid is bottom to top)
      prs_tropomi_mem(:)=prs_tropomi(:)
      if (prs_sfc(imem).gt.prs_tropomi_mem(1)) then
         prs_tropomi_mem(1)=prs_sfc(imem)
      endif   

! Calculate the thicknesses

      do k=1,kend_tropomi
         lnpr_mid=(log(prs_tropomi_mem(k+1))+log(prs_tropomi_mem(k)))/2.
         up_wt=log(prs_tropomi_mem(k))-lnpr_mid
         dw_wt=lnpr_mid-log(prs_tropomi_mem(k+1))
         tl_wt=up_wt+dw_wt      
         tmp_vir_k  = (1.0_r8 + eps*qmr_val(imem,k))*tmp_val(imem,k)
         tmp_vir_kp = (1.0_r8 + eps*qmr_val(imem,k+1))*tmp_val(imem,k+1)
         thick(k)   = Rd*(dw_wt*tmp_vir_k + up_wt*tmp_vir_kp)/tl_wt/grav* &
         log(prs_tropomi_mem(k)/prs_tropomi_mem(k+1))
      enddo

! Process the vertical summation

      do k=1,kend_tropomi
         lnpr_mid=(log(prs_tropomi_mem(k+1))+log(prs_tropomi_mem(k)))/2.
         up_wt=log(prs_tropomi_mem(k))-lnpr_mid
         dw_wt=lnpr_mid-log(prs_tropomi_mem(k+1))
         tl_wt=up_wt+dw_wt

! Convert from VMR to molar density (mol/m^3)
         if(use_log_hcho) then
            hcho_val_conv = (dw_wt*exp(hcho_val(imem,k))+up_wt*exp(hcho_val(imem,k+1)))/tl_wt * &
            (dw_wt*prs_tropomi_mem(k)+up_wt*prs_tropomi_mem(k+1)) / &
            (Ru*(dw_wt*tmp_val(imem,k)+up_wt*tmp_val(imem,k+1)))
         else
            hcho_val_conv = (dw_wt*hcho_val(imem,k)+up_wt*hcho_val(imem,k+1))/tl_wt * &
            (dw_wt*prs_tropomi_mem(k)+up_wt*prs_tropomi_mem(k+1)) / &
            (Ru*(dw_wt*tmp_val(imem,k)+up_wt*tmp_val(imem,k+1)))
         endif
 
! Get expected observation

         expct_val(imem) = expct_val(imem) + thick(k) * hcho_val_conv * &
         avg_kernel(key,k) * amf
      enddo

      if(expct_val(imem).lt.0) then
         zstatus(imem)=20
         expct_val(:)=missing_r8
         write(string1, *) &
         'APM NOTICE: TROPOMI HCHO expected value is negative '
         call error_handler(E_MSG, routine, string1, source)
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
   enddo

! Clean up and return
   deallocate(hcho_val, tmp_val, qmr_val)
   deallocate(thick)
   deallocate(prs_tropomi, prs_tropomi_mem)

end subroutine get_expected_tropomi_hcho_trop_col

!-------------------------------------------------------------------------------

subroutine set_obs_def_tropomi_hcho_trop_col(key, hcho_pressure, hcho_avg_kernel, hcho_amf_trop, hcho_trop_indx, hcho_nlayer)

integer,                           intent(in)   :: key, hcho_nlayer
integer,                           intent(in)   :: hcho_trop_indx
real(r8),                          intent(in)   :: hcho_amf_trop
real(r8), dimension(hcho_nlayer+1), intent(in)   :: hcho_pressure
real(r8), dimension(hcho_nlayer),    intent(in)   :: hcho_avg_kernel

if ( .not. module_initialized ) call initialize_module

if(num_tropomi_hcho_obs >= max_tropomi_hcho_obs) then
   write(string1, *)'Not enough space for tropomi hcho total col obs.'
   write(string2, *)'Can only have max_tropomi_hcho_obs (currently ',max_tropomi_hcho_obs,')'
   call error_handler(E_ERR,'set_obs_def_tropomi_hcho_trop_col',string1,source,revision, &
   revdate,text2=string2)
endif

nlayer(key) = hcho_nlayer
trop_indx(key) = hcho_trop_indx
amf_trop(key) = hcho_amf_trop
pressure(key,1:hcho_nlayer+1) = hcho_pressure(1:hcho_nlayer+1)
avg_kernel(key,1:hcho_nlayer) = hcho_avg_kernel(1:hcho_nlayer)

end subroutine set_obs_def_tropomi_hcho_trop_col

end module obs_def_tropomi_hcho_trop_col_mod

! END DART PREPROCESS MODULE CODE
