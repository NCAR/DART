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
! MLS_HNO3_PROFILE, QTY_HNO3
! END DART PREPROCESS TYPE DEFINITIONS
!
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_mls_hno3_profile_mod, only : get_expected_mls_hno3_profile, &
!                                  read_mls_hno3_profile, &
!                                  write_mls_hno3_profile, &
!                                  interactive_mls_hno3_profile
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(MLS_HNO3_PROFILE)                                                           
!         call get_expected_mls_hno3_profile(state_handle, ens_size, location, obs_def%key, obs_time, expected_obs, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!
! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(MLS_HNO3_PROFILE)
!         call read_mls_hno3_profile(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(MLS_HNO3_PROFILE)
!         call write_mls_hno3_profile(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF
!
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(MLS_HNO3_PROFILE)
!         call interactive_mls_hno3_profile(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_mls_hno3_profile_mod

   use         apm_upper_bdy_mod, only :get_upper_bdy_fld, &
                                        get_MOZART_INT_DATA, &
                                        get_MOZART_REAL_DATA, &
                                        wrf_dart_ubval_interp, &
                                        apm_get_exo_coldens, &
                                        apm_get_upvals, &
                                        apm_interpolate

   use             types_mod, only : r8, MISSING_R8
   
   use         utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, E_ALLMSG, &
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
   
   use          obs_kind_mod, only : QTY_HNO3, QTY_TEMPERATURE, QTY_SURFACE_PRESSURE, &
                                     QTY_PRESSURE, QTY_VAPOR_MIXING_RATIO
   
   use  ensemble_manager_mod, only : ensemble_type
   
   use obs_def_utilities_mod, only : track_status
   
   use mpi_utilities_mod,     only : my_task_id

   use      time_manager_mod, only : time_type, get_date, set_date, get_time, set_time
! get_date gets year, month, day, hour, minute, second from time_type
! get_time gets julian day and seconds from time_type
! set_date sets time_type from year, month, day, hour, minute, second
! set_time sets time_type from julian day and seconds
   implicit none
   private

   public :: write_mls_hno3_profile, &
             read_mls_hno3_profile, &
             interactive_mls_hno3_profile, &
             get_expected_mls_hno3_profile, &
             set_obs_def_mls_hno3_profile

! Storage for the special information required for observations of this type
   integer, parameter    :: max_mls_hno3_obs = 10000000
   integer               :: num_mls_hno3_obs = 0
   integer,  allocatable :: nlayer(:)
   integer,  allocatable :: klev(:)
   integer,  allocatable :: kend(:)
   real(r8), allocatable :: pressure(:,:)
   real(r8), allocatable :: avg_kernel(:,:)
   real(r8), allocatable :: prior(:,:)

! version controlled file description for error handling, do not edit
   character(len=*), parameter :: source   = 'obs_def_mls_hno3_profile_mod.f90'
   character(len=*), parameter :: revision = ''
   character(len=*), parameter :: revdate  = ''
   
   character(len=512) :: string1, string2
   character(len=200) :: upper_data_file
   character(len=200) :: upper_data_model
   character(len=200) :: model
   integer            :: ls_chem_dx, ls_chem_dy, ls_chem_dz, ls_chem_dt
   
   logical, save :: module_initialized = .false.

! Namelist with default values
   logical :: use_log_hno3   = .false.
   integer :: nlayer_model = -9999
   integer :: nlayer_mls = -9999
   integer :: nlayer_mls_hno3_total_col = -9999
   integer :: nlayer_mls_hno3_trop_col = -9999
   integer :: nlayer_mls_hno3_profile = -9999
   
   namelist /obs_def_MLS_HNO3_nml/ upper_data_file, use_log_hno3, nlayer_model, &
   nlayer_mls_hno3_total_col, nlayer_mls_hno3_trop_col, nlayer_mls_hno3_profile, &
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
   call find_namelist_in_file("input.nml", "obs_def_MLS_HNO3_nml", iunit)
   read(iunit, nml = obs_def_MLS_HNO3_nml, iostat = rc)
   call check_namelist_read(iunit, rc, "obs_def_MLS_HNO3_nml")

! Record the namelist values
   if (do_nml_file()) write(nmlfileunit, nml=obs_def_MLS_HNO3_nml)
   if (do_nml_term()) write(     *     , nml=obs_def_MLS_HNO3_nml)

! Check for valid values
   nlayer_mls=nlayer_mls_hno3_profile

   if (nlayer_model < 1) then
      write(string1,*)'obs_def_MLS_HNO3_nml:nlayer_model must be > 0, it is ',nlayer_model
      call error_handler(E_ERR,'initialize_module',string1,source)
   endif
   
   if (nlayer_mls < 1) then
      write(string1,*)'obs_def_MLS_HNO3_nml:nlayer_mls must be > 0, it is ',nlayer_mls
      call error_handler(E_ERR,'initialize_module',string1,source)
   endif
   
   allocate(    nlayer(max_mls_hno3_obs))
   allocate(    klev(max_mls_hno3_obs))
   allocate(    kend(max_mls_hno3_obs))
   allocate(  pressure(max_mls_hno3_obs,nlayer_mls))
   allocate(avg_kernel(max_mls_hno3_obs,nlayer_mls))
   allocate(     prior(max_mls_hno3_obs,nlayer_mls))
   
end subroutine initialize_module

!-------------------------------------------------------------------------------

subroutine read_mls_hno3_profile(key, ifile, fform)

   integer,          intent(out)          :: key
   integer,          intent(in)           :: ifile
   character(len=*), intent(in), optional :: fform

! temporary arrays to hold buffer till we decide if we have enough room

   integer               :: keyin
   integer               :: nlayer_1
   integer               :: klev_1
   integer               :: kend_1
   real(r8), allocatable :: pressure_1(:)
   real(r8), allocatable :: avg_kernel_1(:)
   real(r8), allocatable :: prior_1(:)
   character(len=32)     :: fileformat
   
   integer, SAVE :: counts1 = 0
   
   if ( .not. module_initialized ) call initialize_module
   
   fileformat = "ascii" 
   if(present(fform)) fileformat = adjustl(fform)
   
! Need to know how many layers for this one
   nlayer_1 = read_int_scalar( ifile, fileformat, 'nlayer_1')
   klev_1 = read_int_scalar( ifile, fileformat, 'klev_1')
   kend_1 = read_int_scalar( ifile, fileformat, 'kend_1')
   
   allocate(  pressure_1(nlayer_1))
   allocate(avg_kernel_1(nlayer_1))
   allocate(prior_1(nlayer_1))
   
   call read_r8_array(ifile, nlayer_1, pressure_1,   fileformat, 'pressure_1')
   call read_r8_array(ifile, nlayer_1,   avg_kernel_1, fileformat, 'avg_kernel_1')
   call read_r8_array(ifile, nlayer_1,   prior_1,      fileformat, 'prior_1')
   keyin = read_int_scalar(ifile, fileformat, 'keyin')
   
   counts1 = counts1 + 1
   key     = counts1
   
   if(counts1 > max_mls_hno3_obs) then
      write(string1, *)'Not enough space for mls hno3 obs.'
      write(string2, *)'Can only have max_mls_hno3_obs (currently ',max_mls_hno3_obs,')'
      call error_handler(E_ERR,'read_mls_hno3_profile',string1,source,text2=string2)
   endif
   
   call set_obs_def_mls_hno3_profile(key, pressure_1(1:nlayer_1), avg_kernel_1(1:nlayer_1), &
   prior_1(1:nlayer_1), nlayer_1, klev_1, kend_1)
   
   deallocate(pressure_1, avg_kernel_1, prior_1)
   
end subroutine read_mls_hno3_profile

!-------------------------------------------------------------------------------

subroutine write_mls_hno3_profile(key, ifile, fform)

   integer,          intent(in)           :: key
   integer,          intent(in)           :: ifile
   character(len=*), intent(in), optional :: fform
   
   character(len=32) :: fileformat
   
   if ( .not. module_initialized ) call initialize_module
   
   fileformat = "ascii"
   if(present(fform)) fileformat = adjustl(fform)
   
! nlayer, pressure, avg_kernel, and prior are all scoped in this module
! you can come extend the context strings to include the key if desired.

   call write_int_scalar(ifile,                     nlayer(key), fileformat,'nlayer')
   call write_int_scalar(ifile,                     klev(key), fileformat,'klev')
   call write_int_scalar(ifile,                     kend(key), fileformat,'kend')
   call write_r8_array(  ifile, nlayer(key),  pressure(key,:), fileformat,'pressure')
   call write_r8_array(  ifile, nlayer(key),  avg_kernel(key,:), fileformat,'avg_kernel')
   call write_r8_array(  ifile, nlayer(key),       prior(key,:), fileformat,'prior')
   call write_int_scalar(ifile,                             key, fileformat,'key')
   
end subroutine write_mls_hno3_profile

!-------------------------------------------------------------------------------

subroutine interactive_mls_hno3_profile(key)

   integer, intent(out) :: key
   
   if ( .not. module_initialized ) call initialize_module

! STOP because routine is not finished.
   write(string1,*)'interactive_mls_hno3_profile not yet working.'
   call error_handler(E_ERR, 'interactive_mls_hno3_profile', string1, source)
   
   if(num_mls_hno3_obs >= max_mls_hno3_obs) then
      write(string1, *)'Not enough space for an mls o3 obs.'
      write(string2, *)'Can only have max_mls_hno3_obs (currently ',max_mls_hno3_obs,')'
      call error_handler(E_ERR, 'interactive_mls_hno3_profile', string1, &
                 source, text2=string2)
   endif
   
! Increment the index
   num_mls_hno3_obs = num_mls_hno3_obs + 1
   key              = num_mls_hno3_obs

! Otherwise, prompt for input for the three required beasts

   write(*, *) 'Creating an interactive_mls_hno3_profile observation'
   write(*, *) 'This featue is not setup '

end subroutine interactive_mls_hno3_profile

!-------------------------------------------------------------------------------

subroutine get_expected_mls_hno3_profile(state_handle, ens_size, location, key, obs_time, expct_val, istatus)

   type(ensemble_type), intent(in)  :: state_handle
   type(location_type), intent(in)  :: location
   integer,             intent(in)  :: ens_size
   integer,             intent(in)  :: key
   type(time_type),     intent(in)  :: obs_time
   integer,             intent(out) :: istatus(:)
   real(r8),            intent(out) :: expct_val(:)
   
   character(len=*), parameter :: routine = 'get_expected_mls_hno3_profile'
   character(len=120)          :: data_file
   character(len=*),parameter  :: fld = 'HNO3_VMR_inst'
   type(location_type) :: loc2
   
   integer :: layer_mls,level_mls, klev_mls, kend_mls
   integer :: layer_mdl,level_mdl
   integer :: k,kk,imem,imemm,flg
   integer :: interp_new
   integer :: icnt,ncnt,kstart
   integer :: date_obs,datesec_obs
   integer, dimension(ens_size) :: zstatus,kbnd_1,kbnd_n
   
   real(r8) :: eps, AvogN, Rd, Ru, Cp, grav, msq2cmsq
   real(r8) :: missing,hno3_min,tmp_max
   real(r8) :: level,del_prs,prior_term
   real(r8) :: tmp_vir_k, tmp_vir_kp
   real(r8) :: mloc(3),obs_prs
   real(r8) :: hno3_val_conv, VMR_conv
   real(r8) :: up_wt,dw_wt,tl_wt,lnpr_mid
   real(r8) :: lon_obs,lat_obs,pi,rad2deg

   real(r8), dimension(ens_size) :: hno3_mdl_1, tmp_mdl_1, qmr_mdl_1, prs_mdl_1
   real(r8), dimension(ens_size) :: hno3_mdl_n, tmp_mdl_n, qmr_mdl_n, prs_mdl_n
   real(r8), dimension(ens_size) :: prs_sfc
   
   real(r8), allocatable, dimension(:)   :: thick, prs_mls, prs_mls_mem
   real(r8), allocatable, dimension(:,:) :: hno3_val, tmp_val, qmr_val
   logical  :: return_now,hno3_return_now,tmp_return_now,qmr_return_now
!
! Upper BC variables
   real     :: prs_del,delta,bdy_coef
   real     :: hno3_bot,prs_bot,tmp_bot,qmr_bot
   real     :: hno3_top,prs_top,tmp_top,qmr_top
   real(r8), allocatable, dimension(:)   :: hno3_prf_mdl,tmp_prf_mdl,qmr_prf_mdl
   real(r8), allocatable, dimension(:)   :: prs_mls_top
   
   if ( .not. module_initialized ) call initialize_module
   
   pi       = 4.*atan(1.)
   rad2deg  = 360./(2.*pi)
   eps      = 0.61_r8
   Rd       = 287.05_r8     ! J/(mole-kg)
   Ru       = 8.316_r8      ! J/(mole-kg)
   Cp       = 1006.0        ! J/kg/K
   grav     = 9.8_r8        ! m/s^2
   hno3_min = 1.e-6_r8
   msq2cmsq = 1.e4_r8
   AvogN    = 6.02214e23_r8
   missing  = -888888_r8
   tmp_max  = 600.
   del_prs  = 5000.
   VMR_conv = 28.9644/47.9982
   bdy_coef = 0.95
   prs_del  = 1000.         ! Pa  
! 
! WACCM - MMR
! WRFChem - VMR ppmv
! MLS HNO3 - DU   
!
! to convert from mass mixing ratio (MMR) to volume mixing ratio (VMR) multiply by
! the molar mass of dry air (28.9644 g) and divide by the molar mass of the constituent
! O3 - 47.9982 g
! CO - 28.0101 g
! NO2 - 46.0055 g
! SO2 - 64.0638 g
! CO2 - 44.0096 g
! CH4 - 16.0425 g
!
! to get VMR in ppb multiply by 1e9
! to get VMR in ppm multiply by 1e6   
!
   if(use_log_hno3) then
      hno3_min = log(hno3_min)
   endif
   
! Assign vertical grid information (MLS HNO3 grid is bottom to top)

   layer_mls = nlayer(key)
   level_mls = nlayer(key)+1
   klev_mls  = klev(key)
   kend_mls  = kend(key)
   layer_mdl   = nlayer_model
   level_mdl   = nlayer_model+1

   allocate(prs_mls(layer_mls))
   allocate(prs_mls_mem(layer_mls))
   prs_mls(1:layer_mls)=pressure(key,1:layer_mls)

! Get location infomation

   mloc = get_location(location)
   
   if (mloc(2) >  90.0_r8) then
      mloc(2) =  90.0_r8
   elseif (mloc(2) < -90.0_r8) then
      mloc(2) = -90.0_r8
   endif
   obs_prs=mloc(3)

! You could set a unique error code for each condition and then just return
! without having to issue a warning message. The error codes would then
! show up in the report from 'output_forward_op_errors'

   istatus(:) = 0  ! set this once at the beginning
   return_now=.false.
   
   zstatus=0
   level=0.0_r8
   loc2 = set_location(mloc(1), mloc(2), level, VERTISSURFACE)
   call interpolate(state_handle, ens_size, loc2, QTY_SURFACE_PRESSURE, prs_sfc, zstatus) 

   hno3_mdl_1(:)=missing_r8
   tmp_mdl_1(:)=missing_r8
   qmr_mdl_1(:)=missing_r8
   prs_mdl_1(:)=missing_r8

   do k=1,layer_mdl
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
      call interpolate(state_handle, ens_size, loc2, QTY_HNO3, hno3_mdl_1, zstatus) ! ppmv 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_1, zstatus) ! K 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_1, zstatus) ! kg / kg 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_1, zstatus) ! Pa
!
      interp_new=0
      do imem=1,ens_size
         if(hno3_mdl_1(imem).eq.missing_r8 .or. tmp_mdl_1(imem).eq.missing_r8 .or. &
         qmr_mdl_1(imem).eq.missing_r8 .or. prs_mdl_1(imem).eq.missing_r8) then
            interp_new=1
            exit
         endif
      enddo
      if(interp_new.eq.0) then
         exit
      endif
   enddo

!   write(string1, *) 'APM: hno3 lower bound ',key,hno3_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: tmp lower bound ',key,tmp_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: qmr lower bound ',key,qmr_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: prs lower bound ',key,prs_mdl_1
!   call error_handler(E_MSG, routine, string1, source)

   hno3_mdl_n(:)=missing_r8
   tmp_mdl_n(:)=missing_r8
   qmr_mdl_n(:)=missing_r8
   prs_mdl_n(:)=missing_r8

   do k=layer_mdl,1,-1
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
      call interpolate(state_handle, ens_size, loc2, QTY_HNO3, hno3_mdl_n, zstatus) ! ppmv
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_n, zstatus) ! K
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_n, zstatus) ! kg/kg
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_n, zstatus) ! Pa
!
      interp_new=0
      do imem=1,ens_size
         if(hno3_mdl_n(imem).eq.missing_r8 .or. tmp_mdl_n(imem).eq.missing_r8 .or. &
         qmr_mdl_n(imem).eq.missing_r8 .or. prs_mdl_n(imem).eq.missing_r8) then
            interp_new=1
            exit
         endif
      enddo
      if(interp_new.eq.0) then
         exit
      endif
   enddo

!   write(string1, *) 'APM: hno3 upper bound ',key,hno3_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: tmp upper bound ',key,tmp_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: qmr upper bound ',key,qmr_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: prs upper bound ',key,prs_mdl_n
!   call error_handler(E_MSG, routine, string1, source)

! Get profiles at MLS pressure levels (bottom to top)

   allocate(hno3_val(ens_size,layer_mls))
   allocate(tmp_val(ens_size,layer_mls))
   allocate(qmr_val(ens_size,layer_mls))

   do k=1,layer_mls
      zstatus=0
      loc2 = set_location(mloc(1), mloc(2), prs_mls(k), VERTISPRESSURE)
      call interpolate(state_handle, ens_size, loc2, QTY_HNO3, hno3_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_val(:,k), zstatus)  
!
! Correcting for expected failures near the surface
      do imem=1,ens_size
         if (prs_mls(k).ge.prs_mdl_1(imem)) then
            hno3_val(imem,k) = hno3_mdl_1(imem)
            tmp_val(imem,k) = tmp_mdl_1(imem)
            qmr_val(imem,k) = qmr_mdl_1(imem)
         endif
!
! Correcting for expected failures near the top
         if (prs_mls(k).le.prs_mdl_n(imem)) then
            hno3_val(imem,k) = hno3_mdl_n(imem)
            tmp_val(imem,k) = tmp_mdl_n(imem)
            qmr_val(imem,k) = qmr_mdl_n(imem)
         endif
      enddo

!      write(string1, *)'APM: hno3 ',key,k,hno3_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: tmp ',key,k,tmp_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: qmr ',key,k,qmr_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!
! Convert units for hno3 from ppmv
      hno3_val(:,k) = hno3_val(:,k) * 1.e-6_r8
   enddo
   hno3_mdl_1(:) = hno3_mdl_1(:) * 1.e-6_r8
   hno3_mdl_n(:) = hno3_mdl_n(:) * 1.e-6_r8
!
! Use large scale hno3 data above the regional model top
! MLS vertical grid is from bottom to top
   kstart=-1
   do imem=1,ens_size
      if (prs_mls(layer_mls).lt.prs_mdl_n(imem)) then
         do k=1,layer_mls
            if (prs_mls(k).le.prs_mdl_n(imem)) then
               kstart=k-1
               exit
            endif
         enddo
         ncnt=layer_mls-kstart+1
         allocate(prs_mls_top(ncnt))
         allocate(hno3_prf_mdl(ncnt),tmp_prf_mdl(ncnt),qmr_prf_mdl(ncnt))
         do k=kstart,layer_mls
            kk=k-kstart+1
            prs_mls_top(kk)=prs_mls(k)
         enddo
         prs_mls_top(:)=prs_mls_top(:)/100.
!
         lon_obs=mloc(1)/rad2deg
         lat_obs=mloc(2)/rad2deg
         call get_time(obs_time,datesec_obs,date_obs)
!
         data_file=trim(upper_data_file)
         model=trim(upper_data_model)
         
!         write(string1, *) &
!         'APM: Model pressure 1, n', kstart, prs_mdl_1(imem), prs_mdl_n(imem)
!         call error_handler(E_ALLMSG, routine, string1, source) 
!         do k=1,layer_mls
!            write(string1, *) &
!            'APM: MLS pressure ', k, prs_mls(k)
!            call error_handler(E_ALLMSG, routine, string1, source) 
!         enddo
         call get_upper_bdy_fld(fld,model,data_file,ls_chem_dx,ls_chem_dy, &
         ls_chem_dz,ls_chem_dt,lon_obs,lat_obs,prs_mls_top, &
         ncnt,hno3_prf_mdl,tmp_prf_mdl,qmr_prf_mdl,date_obs,datesec_obs)
!         do k=1,ncnt
!            write(string1, *) &
!            'APM: Uppr Bdy values ', k, hno3_prf_mdl(k), tmp_prf_mdl(k), qmr_prf_mdl(k)
!            call error_handler(E_ALLMSG, routine, string1, source) 
!         enddo
         
!
! Impose ensemble perturbations from level kstart-1
         do k=kstart,layer_mls
            kk=k-kstart+1
            hno3_val(imem,k)=hno3_prf_mdl(kk)*hno3_val(imem,kstart-1)/ &
            (sum(hno3_val(:,kstart-1))/real(ens_size))
            tmp_val(imem,k)=tmp_prf_mdl(kk)*tmp_val(imem,kstart-1)/ &
            (sum(tmp_val(:,kstart-1))/real(ens_size))
            qmr_val(imem,k)=qmr_prf_mdl(kk)*qmr_val(imem,kstart-1)/ &
            (sum(qmr_val(:,kstart-1))/real(ens_size))
         enddo
         deallocate(prs_mls_top)
         deallocate(hno3_prf_mdl,tmp_prf_mdl,qmr_prf_mdl)
      endif             
   enddo
!
! Print full profile examples
!   do imem=1,1
!      do k=1,layer_mls
!         write(string1, *) &
!         'APM: prs,hno3,tmp,qmr ',k,prs_mls(k)/100.,hno3_val(imem,k), &
!         tmp_val(imem,k),qmr_val(imem,k)
!         call error_handler(E_MSG, routine, string1, source)
!      enddo
!   enddo
!
! Check full profile for negative values
   do imem=1,ens_size
      flg=0
      do k=1,layer_mls
         if(hno3_val(imem,k).lt.0. .or. tmp_val(imem,k).lt.0. .or. &
         qmr_val(imem,k).lt.0.) then
            flg=1   
            write(string1, *) &
            'APM: Recentered full profile has negative values for key,imem ',key,imem
            call error_handler(E_ALLMSG, routine, string1, source)
         endif
      enddo
      if(flg.eq.1) then
         zstatus(imem)=20
         expct_val(:)=missing_r8
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
   enddo
!
! Calculate the expected retrievals
   istatus(:)=0
   zstatus(:)=0
   expct_val(:)=0.0
   allocate(thick(layer_mls))

   do imem=1,ens_size
! Define upper and lower values for layer grid
! (MLS HNO3 grid is bottom to top) prs is in Pa
      prs_mls_mem(:)=prs_mls(:)      
! Definitions for k=1 or k=layer_mls
      prs_bot=prs_sfc(imem)
      if (prs_bot.le.prs_mls_mem(1)) then
         prs_bot=prs_mls_mem(1)+prs_del
      endif   
! Bottom terms
      hno3_bot=hno3_mdl_1(imem)
      tmp_bot=tmp_mdl_1(imem)
      qmr_bot=qmr_mdl_1(imem)
! Top terms
      prs_top=prs_mls(layer_mls)+(prs_mls(layer_mls)-prs_mls(layer_mls-1))/2.
      if(prs_top.le.0.) prs_top=bdy_coef*prs_mls(layer_mls)
! hno3
      delta=(hno3_val(imem,layer_mls)-hno3_val(imem,layer_mls-1))/ &
      (prs_mls(layer_mls)-prs_mls(layer_mls-1))
      hno3_top=hno3_val(imem,layer_mls) + delta*(prs_top-prs_mls(layer_mls))
      if(hno3_top.le.0.) then
         if(delta.le.0.) hno3_top=bdy_coef*hno3_val(imem,layer_mls)
         if(delta.gt.0.) hno3_top=(2.-bdy_coef)*hno3_val(imem,layer_mls)
      endif
! tmp
      delta=(tmp_val(imem,layer_mls)-tmp_val(imem,layer_mls-1))/ &
      (prs_mls(layer_mls)-prs_mls(layer_mls-1))
      tmp_top=tmp_val(imem,layer_mls) + delta*(prs_top-prs_mls(layer_mls))
      if(tmp_top.le.0.) then
         if(delta.le.0.) tmp_top=bdy_coef*tmp_val(imem,layer_mls)
         if(delta.gt.0.) tmp_top=(2.-bdy_coef)*tmp_val(imem,layer_mls)
      endif
! qmr
      delta=(qmr_val(imem,layer_mls)-qmr_val(imem,layer_mls-1))/ &
      (prs_mls(layer_mls)-prs_mls(layer_mls-1))
      qmr_top=qmr_val(imem,layer_mls) + delta*(prs_top-prs_mls(layer_mls))
      if(qmr_top.le.0.) then
         if(delta.le.0.) qmr_top=bdy_coef*qmr_val(imem,layer_mls)
         if(delta.gt.0.) qmr_top=(2.-bdy_coef)*qmr_val(imem,layer_mls)
      endif
!
! VERTICAL SUMMATION
! k=1 term      
      k=1
! hno3 term (Units are VMR, calculate layer average)
         lnpr_mid=(log(prs_mls_mem(k+1))+log(prs_bot))/2.
         up_wt=log(prs_bot)-lnpr_mid
         dw_wt=lnpr_mid-log(prs_mls_mem(k+1))
         tl_wt=up_wt+dw_wt
         if(use_log_hno3) then
            hno3_val_conv = (dw_wt*exp(hno3_bot)+up_wt*exp(hno3_val(imem,k+1)))/tl_wt
         else
            hno3_val_conv = (dw_wt*hno3_bot+up_wt*hno3_val(imem,k+1))/tl_wt
         endif
         prior_term=-avg_kernel(key,k)
         if(k.eq.klev_mls) prior_term=1.-avg_kernel(key,k)
! expected retrieval sum
         expct_val(imem) = expct_val(imem) + hno3_val_conv * &
         avg_kernel(key,k) + prior_term*prior(key,k)

!         write(string1, *)'APM: expected retr ',k,expct_val(imem), &
!         avg_kernel(key,k), prior(key,k)
!         call error_handler(E_MSG, routine, string1, source)
!
! k=layer_mls term
      k=layer_mls
         lnpr_mid=(log(prs_top)+log(prs_mls_mem(k)))/2.
         up_wt=log(prs_mls_mem(k))-lnpr_mid
         dw_wt=lnpr_mid-log(prs_top)
         tl_wt=up_wt+dw_wt
! HNO3 term (Units are VMR, calculate layer average)
         if(use_log_hno3) then
            hno3_val_conv = (dw_wt*exp(hno3_val(imem,k))+up_wt*exp(hno3_top))/tl_wt
         else
            hno3_val_conv = (dw_wt*hno3_val(imem,k)+up_wt*hno3_top)/tl_wt
         endif
         prior_term=-avg_kernel(key,k)
         if(k.eq.klev_mls) prior_term=1.-avg_kernel(key,k)
! expected retrieval sum
         expct_val(imem) = expct_val(imem) + hno3_val_conv * &
         avg_kernel(key,k) + prior_term*prior(key,k)

!         write(string1, *)'APM: expected retr ',k,expct_val(imem), &
!         avg_kernel(key,k), prior(key,k)
!         call error_handler(E_MSG, routine, string1, source)
!
! remaining terms
      do k=2,layer_mls-1
         prs_bot=(prs_mls_mem(k-1)+prs_mls_mem(k))/2.
         prs_top=(prs_mls_mem(k)+prs_mls_mem(k+1))/2.
         hno3_bot=(hno3_val(imem,k-1)+hno3_val(imem,k))/2.
         hno3_top=(hno3_val(imem,k)+hno3_val(imem,k+1))/2.
         tmp_bot=(tmp_val(imem,k-1)+tmp_val(imem,k))/2.
         tmp_top=(tmp_val(imem,k)+tmp_val(imem,k+1))/2.
         qmr_bot=(qmr_val(imem,k-1)+qmr_val(imem,k))/2.
         qmr_top=(qmr_val(imem,k)+qmr_val(imem,k+1))/2.
         lnpr_mid=(log(prs_top)+log(prs_mls_mem(k)))/2.
         up_wt=log(prs_bot)-lnpr_mid
         dw_wt=lnpr_mid-log(prs_mls_mem(k+1))
         tl_wt=up_wt+dw_wt
! hno3 term (Units are VMR, calculate layer average)
         if(use_log_hno3) then
            hno3_val_conv = (dw_wt*exp(hno3_bot)+up_wt*exp(hno3_top))/tl_wt
         else
            hno3_val_conv = (dw_wt*hno3_bot+up_wt*hno3_top)/tl_wt
         endif
         prior_term=-avg_kernel(key,k)
         if(k.eq.klev_mls) prior_term=1.-avg_kernel(key,k)
! expected retrieval
         expct_val(imem) = expct_val(imem) + hno3_val_conv * &
         avg_kernel(key,k) + prior_term*prior(key,k)

!         write(string1, *)'APM: expected retr ',k,expct_val(imem), &
!         avg_kernel(key,k), prior(key,k)
!         call error_handler(E_MSG, routine, string1, source)
      enddo
!      write(string1, *)'APM: FINAL EXPECTED VALUE ',expct_val(imem)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'  '
!      call error_handler(E_MSG, routine, string1, source)
      if(expct_val(imem).lt.0) then
         zstatus(imem)=20
         expct_val(:)=missing_r8
!         write(string1, *) 'APM NOTICE: MLS HNO3 expected value is negative '
!         call error_handler(E_MSG, routine, string1, source)
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
   enddo

! Clean up and return
   deallocate(hno3_val, tmp_val, qmr_val)
   deallocate(thick)
   deallocate(prs_mls, prs_mls_mem)

end subroutine get_expected_mls_hno3_profile

!-------------------------------------------------------------------------------

subroutine set_obs_def_mls_hno3_profile(key, hno3_pressure, hno3_avg_kernel, hno3_prior, &
hno3_nlayer, hno3_klev, hno3_kend)

   integer,                           intent(in)   :: key, hno3_nlayer, hno3_klev, hno3_kend
   real(r8), dimension(hno3_nlayer),  intent(in)   :: hno3_pressure
   real(r8), dimension(hno3_nlayer),  intent(in)   :: hno3_avg_kernel
   real(r8), dimension(hno3_nlayer),  intent(in)   :: hno3_prior
   
   if ( .not. module_initialized ) call initialize_module
   
   if(num_mls_hno3_obs >= max_mls_hno3_obs) then
      write(string1, *)'Not enough space for mls hno3 obs.'
      write(string2, *)'Can only have max_mls_hno3_obs (currently ',max_mls_hno3_obs,')'
      call error_handler(E_ERR,'set_obs_def_mls_hno3_profile',string1,source,revision, &
      revdate,text2=string2)
   endif
   
   nlayer(key) = hno3_nlayer
   klev(key) = hno3_klev
   kend(key) = hno3_kend
   pressure(key,1:hno3_nlayer)   = hno3_pressure(1:hno3_nlayer)
   avg_kernel(key,1:hno3_nlayer) = hno3_avg_kernel(1:hno3_nlayer)
   prior(key,1:hno3_nlayer)      = hno3_prior(1:hno3_nlayer)
   
end subroutine set_obs_def_mls_hno3_profile

end module obs_def_mls_hno3_profile_mod

! END DART PREPROCESS MODULE CODE
