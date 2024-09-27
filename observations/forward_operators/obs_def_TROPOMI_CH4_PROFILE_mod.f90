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
! TROPOMI_CH4_PROFILE, QTY_CH4
! END DART PREPROCESS TYPE DEFINITIONS
!
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_tropomi_ch4_profile_mod, only : get_expected_tropomi_ch4_profile, &
!                                  read_tropomi_ch4_profile, &
!                                  write_tropomi_ch4_profile, &
!                                  interactive_tropomi_ch4_profile
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(TROPOMI_CH4_PROFILE)                                                           
!         call get_expected_tropomi_ch4_profile(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!
! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(TROPOMI_CH4_PROFILE)
!         call read_tropomi_ch4_profile(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(TROPOMI_CH4_PROFILE)
!         call write_tropomi_ch4_profile(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF
!
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(TROPOMI_CH4_PROFILE)
!         call interactive_tropomi_ch4_profile(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_tropomi_ch4_profile_mod

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

use          obs_kind_mod, only : QTY_CH4, QTY_TEMPERATURE, QTY_SURFACE_PRESSURE, &
                                  QTY_PRESSURE, QTY_VAPOR_MIXING_RATIO

use  ensemble_manager_mod, only : ensemble_type

use obs_def_utilities_mod, only : track_status

implicit none
private

public :: write_tropomi_ch4_profile, &
          read_tropomi_ch4_profile, &
          interactive_tropomi_ch4_profile, &
          get_expected_tropomi_ch4_profile, &
          set_obs_def_tropomi_ch4_profile

! Storage for the special information required for observations of this type
integer, parameter    :: max_tropomi_ch4_obs = 10000000
integer               :: num_tropomi_ch4_obs = 0
integer,  allocatable :: nlayer(:)
real(r8), allocatable :: altitude(:,:)
real(r8), allocatable :: avg_kernel(:,:)
real(r8), allocatable :: prior(:,:)

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_tropomi_ch4_profile_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

character(len=512) :: string1, string2
character(len=200) :: upper_data_file
character(len=200) :: upper_data_model
character(len=200) :: model
integer            :: ls_chem_dx, ls_chem_dy, ls_chem_dz, ls_chem_dt

logical, save :: module_initialized = .false.

! Namelist with default values
logical :: use_log_ch4   = .false.
integer :: nlayer_model = -9999
integer :: nlayer_tropomi = -9999
integer :: nlayer_tropomi_ch4_total_col = -9999
integer :: nlayer_tropomi_ch4_trop_col = -9999
integer :: nlayer_tropomi_ch4_profile = -9999

namelist /obs_def_TROPOMI_CH4_nml/ upper_data_file, use_log_ch4, nlayer_model, &
nlayer_tropomi_ch4_total_col, nlayer_tropomi_ch4_trop_col, nlayer_tropomi_ch4_profile, &
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
call find_namelist_in_file("input.nml", "obs_def_TROPOMI_CH4_nml", iunit)
read(iunit, nml = obs_def_TROPOMI_CH4_nml, iostat = rc)
call check_namelist_read(iunit, rc, "obs_def_TROPOMI_CH4_nml")

! Record the namelist values
if (do_nml_file()) write(nmlfileunit, nml=obs_def_TROPOMI_CH4_nml)
if (do_nml_term()) write(     *     , nml=obs_def_TROPOMI_CH4_nml)
nlayer_tropomi=nlayer_tropomi_ch4_profile

! Check for valid values

if (nlayer_model < 1) then
   write(string1,*)'obs_def_TROPOMI_CH4_nml:nlayer_model must be > 0, it is ',nlayer_model
   call error_handler(E_ERR,'initialize_module',string1,source)
endif

if (nlayer_tropomi < 1) then
   write(string1,*)'obs_def_TROPOMI_CH4_nml:nlayer_tropomi must be > 0, it is ',nlayer_tropomi
   call error_handler(E_ERR,'initialize_module',string1,source)
endif

allocate(    nlayer(max_tropomi_ch4_obs))
allocate(  altitude(max_tropomi_ch4_obs,nlayer_tropomi+1))
allocate(avg_kernel(max_tropomi_ch4_obs,nlayer_tropomi))
allocate(     prior(max_tropomi_ch4_obs,nlayer_tropomi))

end subroutine initialize_module

!-------------------------------------------------------------------------------

subroutine read_tropomi_ch4_profile(key, ifile, fform)

integer,          intent(out)          :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

! temporary arrays to hold buffer till we decide if we have enough room

integer               :: keyin
integer               :: nlayer_1
real(r8), allocatable :: altitude_1(:)
real(r8), allocatable :: avg_kernel_1(:)
real(r8), allocatable :: prior_1(:)
character(len=32)     :: fileformat

integer, SAVE :: counts1 = 0

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii" 
if(present(fform)) fileformat = adjustl(fform)

! Need to know how many layers for this one
nlayer_1 = read_int_scalar( ifile, fileformat, 'nlayer_1')

allocate(  altitude_1(nlayer_1+1))
allocate(avg_kernel_1(nlayer_1))
allocate(     prior_1(nlayer_1))

call read_r8_array(ifile, nlayer_1+1, altitude_1,   fileformat, 'altitude_1')
call read_r8_array(ifile, nlayer_1,   avg_kernel_1, fileformat, 'avg_kernel_1')
call read_r8_array(ifile, nlayer_1,   prior_1, fileformat, 'prior_1')
keyin = read_int_scalar(ifile, fileformat, 'keyin')

counts1 = counts1 + 1
key     = counts1

if(counts1 > max_tropomi_ch4_obs) then
   write(string1, *)'Not enough space for tropomi ch4 profile obs.'
   write(string2, *)'Can only have max_tropomi_ch4_obs (currently ',max_tropomi_ch4_obs,')'
   call error_handler(E_ERR,'read_tropomi_ch4_profile',string1,source,text2=string2)
endif

call set_obs_def_tropomi_ch4_profile(key, altitude_1, avg_kernel_1, prior_1, nlayer_1)

deallocate(altitude_1, avg_kernel_1, prior_1)

end subroutine read_tropomi_ch4_profile

!-------------------------------------------------------------------------------

subroutine write_tropomi_ch4_profile(key, ifile, fform)

integer,          intent(in)           :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"
if(present(fform)) fileformat = adjustl(fform)

call write_int_scalar(ifile,                     nlayer(key), fileformat,'nlayer')
call write_r8_array(  ifile, nlayer(key)+1,      altitude(key,:), fileformat,'altitude')
call write_r8_array(  ifile, nlayer(key),        avg_kernel(key,:), fileformat,'avg_kernel')
call write_r8_array(  ifile, nlayer(key),        prior(key,:), fileformat,'prior')
call write_int_scalar(ifile,                     key, fileformat,'key')

end subroutine write_tropomi_ch4_profile

!-------------------------------------------------------------------------------

subroutine interactive_tropomi_ch4_profile(key)

integer, intent(out) :: key

if ( .not. module_initialized ) call initialize_module

! STOP because routine is not finished.
write(string1,*)'interactive_tropomi_ch4_profile not yet working.'
call error_handler(E_ERR, 'interactive_tropomi_ch4_profile', string1, source)

if(num_tropomi_ch4_obs >= max_tropomi_ch4_obs) then
   write(string1, *)'Not enough space for an tropomi ch4 profile obs.'
   write(string2, *)'Can only have max_tropomi_ch4_obs (currently ',max_tropomi_ch4_obs,')'
   call error_handler(E_ERR, 'interactive_tropomi_ch4_profile', string1, &
   source, text2=string2)
endif

! Increment the index
num_tropomi_ch4_obs = num_tropomi_ch4_obs + 1
key            = num_tropomi_ch4_obs

! Otherwise, prompt for input for the three required beasts

write(*, *) 'Creating an interactive_tropomi_ch4_profile observation'
write(*, *) 'This featue is not setup '

end subroutine interactive_tropomi_ch4_profile

!-------------------------------------------------------------------------------

subroutine get_expected_tropomi_ch4_profile(state_handle, ens_size, location, key, expct_val, istatus)

   type(ensemble_type), intent(in)  :: state_handle
   type(location_type), intent(in)  :: location
   integer,             intent(in)  :: ens_size
   integer,             intent(in)  :: key
   integer,             intent(out) :: istatus(:)
   real(r8),            intent(out) :: expct_val(:)
   
   character(len=*), parameter :: routine = 'get_expected_tropomi_ch4_profile'
   type(location_type) :: loc2
   
   integer :: layer_tropomi,level_tropomi
   integer :: layer_mdl,level_mdl
   integer :: k,imem,kend_tropomi
   integer :: interp_new
   integer, dimension(ens_size) :: zstatus
   
   real(r8) :: eps, AvogN, Rd, Ru, Cp, grav, msq2cmsq
   real(r8) :: missing,ch4_min, tmp_max
   real(r8) :: level,del_prs
   real(r8) :: tmp_vir_k, tmp_vir_kp
   real(r8) :: mloc(3)
   real(r8) :: ch4_val_conv,amf
   real(r8) :: up_wt,dw_wt,tl_wt,lnpr_mid
   real(r8), dimension(ens_size) :: ch4_mdl_1, tmp_mdl_1, qmr_mdl_1, prs_mdl_1
   real(r8), dimension(ens_size) :: ch4_mdl_n, tmp_mdl_n, qmr_mdl_n, prs_mdl_n
   real(r8), dimension(ens_size) :: prs_sfc
   
   real(r8), allocatable, dimension(:)   :: thick, prs_tropomi, prs_tropomi_mem
   real(r8), allocatable, dimension(:,:) :: ch4_val, tmp_val, qmr_val
   logical  :: return_now,ch4_return_now,tmp_return_now,qmr_return_now
   
   if ( .not. module_initialized ) call initialize_module
   
   eps      =  0.61_r8
   Rd       = 287.05_r8     ! J/kg
   Ru       = 8.316_r8      ! J/kg
   Cp       = 1006.0        ! J/kg/K
   grav     =   9.8_r8
   ch4_min  = 1.e-6_r8
   msq2cmsq = 1.e4_r8
   AvogN    = 6.02214e23_r8
   missing  = -888888_r8
   tmp_max  = 600.
   del_prs  = 5000.
   
   if(use_log_ch4) then
      ch4_min = log(ch4_min)
   endif

! Assign vertical grid information

   layer_tropomi = nlayer(key)
   level_tropomi = nlayer(key)+1
   kend_tropomi  = nlayer(key)
   layer_mdl = nlayer_model
   level_mdl = nlayer_model+1
   
   allocate(prs_tropomi(level_tropomi))
   allocate(prs_tropomi_mem(level_tropomi))
   prs_tropomi(1:level_tropomi)=altitude(key,1:level_tropomi)

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

   ch4_mdl_1(:)=missing_r8
   tmp_mdl_1(:)=missing_r8
   qmr_mdl_1(:)=missing_r8
   prs_mdl_1(:)=missing_r8

   do k=1,layer_mdl
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
      call interpolate(state_handle, ens_size, loc2, QTY_CH4, ch4_mdl_1, zstatus) ! ppmv 
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_1, zstatus) ! K 
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_1, zstatus) ! kg / kg 
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_1, zstatus) ! Pa
!
      interp_new=0
      do imem=1,ens_size
         if(ch4_mdl_1(imem).eq.missing_r8 .or. tmp_mdl_1(imem).eq.missing_r8 .or. &
         qmr_mdl_1(imem).eq.missing_r8 .or. prs_mdl_1(imem).eq.missing_r8) then
            interp_new=1
            exit
         endif
      enddo
      if(interp_new.eq.0) then
         exit
      endif    
   enddo

!   write(string1, *)'APM: ch4 lower bound 1 ',ch4_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: tmp lower bound 1 ',tmp_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: qmr lower bound 1 ',qmr_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: prs lower bound 1 ',prs_mdl_1
!   call error_handler(E_MSG, routine, string1, source)

   ch4_mdl_n(:)=missing_r8
   tmp_mdl_n(:)=missing_r8
   qmr_mdl_n(:)=missing_r8
   prs_mdl_n(:)=missing_r8

   do k=layer_mdl,1,-1
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)  
      call interpolate(state_handle, ens_size, loc2, QTY_CH4, ch4_mdl_n, zstatus) ! ppmv
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_n, zstatus) ! K 
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_n, zstatus) ! kg / kg 
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_n, zstatus) ! Pa
!
      interp_new=0
      do imem=1,ens_size
         if(ch4_mdl_n(imem).eq.missing_r8 .or. tmp_mdl_n(imem).eq.missing_r8 .or. &
         qmr_mdl_n(imem).eq.missing_r8 .or. prs_mdl_n(imem).eq.missing_r8) then
            interp_new=1
            exit
         endif
      enddo
      if(interp_new.eq.0) then
         exit
      endif    
   enddo

!   write(string1, *)'APM: ch4 upper bound 1 ',ch4_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: tmp upper bound 1 ',tmp_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: qmr upper bound 1 ',qmr_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: prs upper bound 1 ',prs_mdl_n
!   call error_handler(E_MSG, routine, string1, source)

! Get profiles at TROPOMI pressure levels

   allocate(ch4_val(ens_size,level_tropomi))
   allocate(tmp_val(ens_size,level_tropomi))
   allocate(qmr_val(ens_size,level_tropomi))

   do k=1,level_tropomi
      zstatus=0
      loc2 = set_location(mloc(1), mloc(2), prs_tropomi(k), VERTISPRESSURE)
      call interpolate(state_handle, ens_size, loc2, QTY_CH4, ch4_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_val(:,k), zstatus)  
!
! Correcting for expected failures near the surface
      do imem=1,ens_size
         if (prs_tropomi(k).ge.prs_mdl_1(imem)) then
            ch4_val(imem,k) = ch4_mdl_1(imem)
            tmp_val(imem,k) = tmp_mdl_1(imem)
            qmr_val(imem,k) = qmr_mdl_1(imem)
            cycle
         endif

! Correcting for expected failures near the top
         if (prs_tropomi(k).le.prs_mdl_n(imem)) then
            ch4_val(imem,k) = ch4_mdl_n(imem)
            tmp_val(imem,k) = tmp_mdl_n(imem)
            qmr_val(imem,k) = qmr_mdl_n(imem)
            cycle
         endif
      enddo
!
!      write(string1, *)'APM: ch4 ',k,ch4_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: tmp ',k,tmp_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: qmr ',k,qmr_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!
! Check data for missing values      
      do imem=1,ens_size
         if(ch4_val(imem,k).eq.missing_r8 .or. tmp_val(imem,k).eq.missing_r8 .or. &
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
! Convert units for ch4 from ppmv
      ch4_val(:,k) = ch4_val(:,k) * 1.e-6_r8
   enddo

   istatus=0
   zstatus(:)=0
   expct_val(:)=0.0
   allocate(thick(layer_tropomi))

   do imem=1,ens_size
! Adjust the TROPOMI pressure for WRF-Chem lower/upper boudary pressure
! (TROPOMI CH4 vertical grid is bottom to top)
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
         if(use_log_ch4) then
            ch4_val_conv = (dw_wt*exp(ch4_val(imem,k))+up_wt*exp(ch4_val(imem,k+1)))/tl_wt * &
            (dw_wt*prs_tropomi_mem(k)+up_wt*prs_tropomi_mem(k+1)) / &
            (Ru*(dw_wt*tmp_val(imem,k)+up_wt*tmp_val(imem,k+1)))
         else
            ch4_val_conv = (dw_wt*ch4_val(imem,k)+up_wt*ch4_val(imem,k+1))/tl_wt * &
            (dw_wt*prs_tropomi_mem(k)+up_wt*prs_tropomi_mem(k+1)) / &
            (Ru*(dw_wt*tmp_val(imem,k)+up_wt*tmp_val(imem,k+1)))
         endif
 
! Get expected observation

         expct_val(imem) = expct_val(imem) + thick(k) * ch4_val_conv * &
         avg_kernel(key,k)
      enddo

      if(expct_val(imem).lt.0) then
         zstatus(imem)=20
         expct_val(:)=missing_r8
         write(string1, *) &
         'APM NOTICE: TROPOMI CH4 expected value is negative '
         call error_handler(E_MSG, routine, string1, source)
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
   enddo

! Clean up and return
   deallocate(ch4_val, tmp_val, qmr_val)
   deallocate(thick)
   deallocate(prs_tropomi, prs_tropomi_mem)

end subroutine get_expected_tropomi_ch4_profile

!-------------------------------------------------------------------------------

subroutine set_obs_def_tropomi_ch4_profile(key, ch4_altitude, ch4_avg_kernel, ch4_prior, ch4_nlayer)

integer,                            intent(in)  :: key, ch4_nlayer
real(r8), dimension(ch4_nlayer+1),  intent(in)  :: ch4_altitude
real(r8), dimension(ch4_nlayer),    intent(in)  :: ch4_avg_kernel
real(r8), dimension(ch4_nlayer),    intent(in)  :: ch4_prior

if ( .not. module_initialized ) call initialize_module

if(num_tropomi_ch4_obs >= max_tropomi_ch4_obs) then
   write(string1, *)'Not enough space for tropomi ch4 trop col obs.'
   write(string2, *)'Can only have max_tropomi_ch4_obs (currently ',max_tropomi_ch4_obs,')'
   call error_handler(E_ERR,'set_obs_def_tropomi_ch4_profile',string1,source,revision, &
   revdate,text2=string2)
endif

nlayer(key) = ch4_nlayer
altitude(key,1:ch4_nlayer+1) = ch4_altitude(1:ch4_nlayer+1)
avg_kernel(key,1:ch4_nlayer) = ch4_avg_kernel(1:ch4_nlayer)
prior(key,1:ch4_nlayer)      = ch4_prior(1:ch4_nlayer)

end subroutine set_obs_def_tropomi_ch4_profile

end module obs_def_tropomi_ch4_profile_mod

! END DART PREPROCESS MODULE CODE
