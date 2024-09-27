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
! CRIS_NH3_CPSR, QTY_NH3
! END DART PREPROCESS TYPE DEFINITIONS
!
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_cris_nh3_cpsr_mod, only : get_expected_cris_nh3_cpsr, &
!                                  read_cris_nh3_cpsr, &
!                                  write_cris_nh3_cpsr, &
!                                  interactive_cris_nh3_cpsr
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(CRIS_NH3_CPSR)                                                           
!         call get_expected_cris_nh3_cpsr(state_handle, ens_size, location, obs_def%key, obs_time, expected_obs, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!
! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(CRIS_NH3_CPSR)
!         call read_cris_nh3_cpsr(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(CRIS_NH3_CPSR)
!         call write_cris_nh3_cpsr(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF
!
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(CRIS_NH3_CPSR)
!         call interactive_cris_nh3_cpsr(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_cris_nh3_cpsr_mod

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
   
   use          obs_kind_mod, only : QTY_NH3, QTY_TEMPERATURE, QTY_SURFACE_PRESSURE, &
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

   public :: write_cris_nh3_cpsr, &
             read_cris_nh3_cpsr, &
          interactive_cris_nh3_cpsr, &
          get_expected_cris_nh3_cpsr, &
          set_obs_def_cris_nh3_cpsr

! Storage for the special information required for observations of this type
   integer, parameter    :: max_cris_nh3_obs = 10000000
   integer               :: num_cris_nh3_obs = 0
   integer,  allocatable :: nlayer(:)
   integer,  allocatable :: klev(:)
   integer,  allocatable :: kend(:)
   real(r8), allocatable :: pressure(:,:)
   real(r8), allocatable :: avg_kernel(:,:)
   real(r8), allocatable :: prior(:,:)

! version controlled file description for error handling, do not edit
   character(len=*), parameter :: source   = 'obs_def_cris_nh3_cpsr_mod.f90'
   character(len=*), parameter :: revision = ''
   character(len=*), parameter :: revdate  = ''
   
   character(len=512) :: string1, string2
   character(len=200) :: upper_data_file
   character(len=200) :: upper_data_model
   character(len=200) :: model
   integer            :: ls_chem_dx, ls_chem_dy, ls_chem_dz, ls_chem_dt
   
   logical, save :: module_initialized = .false.

! Namelist with default values
   logical :: use_log_nh3   = .false.
   integer :: nlayer_model = -9999
   integer :: nlayer_cris = -9999
   integer :: nlayer_cris_nh3_total_col = -9999
   integer :: nlayer_cris_nh3_trop_col = -9999
   integer :: nlayer_cris_nh3_profile = -9999
   
   namelist /obs_def_CRIS_NH3_nml/ upper_data_file, use_log_nh3, nlayer_model, &
   nlayer_cris_nh3_total_col, nlayer_cris_nh3_trop_col, nlayer_cris_nh3_profile, &
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
   call find_namelist_in_file("input.nml", "obs_def_CRIS_NH3_nml", iunit)
   read(iunit, nml = obs_def_CRIS_NH3_nml, iostat = rc)
   call check_namelist_read(iunit, rc, "obs_def_CRIS_NH3_nml")

! Record the namelist values
   if (do_nml_file()) write(nmlfileunit, nml=obs_def_CRIS_NH3_nml)
   if (do_nml_term()) write(     *     , nml=obs_def_CRIS_NH3_nml)

! Check for valid values
   nlayer_cris=nlayer_cris_nh3_profile

   if (nlayer_model < 1) then
      write(string1,*)'obs_def_CRIS_NH3_nml:nlayer_model must be > 0, it is ',nlayer_model
      call error_handler(E_ERR,'initialize_module',string1,source)
   endif
   
   if (nlayer_cris < 1) then
      write(string1,*)'obs_def_CRIS_NH3_nml:nlayer_cris must be > 0, it is ',nlayer_cris
      call error_handler(E_ERR,'initialize_module',string1,source)
   endif
   
   allocate(    nlayer(max_cris_nh3_obs))
   allocate(    klev(max_cris_nh3_obs))
   allocate(    kend(max_cris_nh3_obs))
   allocate(  pressure(max_cris_nh3_obs,nlayer_cris+1))
   allocate(avg_kernel(max_cris_nh3_obs,nlayer_cris))
   allocate(     prior(max_cris_nh3_obs,nlayer_cris))
   
end subroutine initialize_module

!-------------------------------------------------------------------------------

subroutine read_cris_nh3_cpsr(key, ifile, fform)

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
   
   allocate(  pressure_1(nlayer_1+1))
   allocate(avg_kernel_1(nlayer_1))
   allocate(prior_1(nlayer_1))
   
   call read_r8_array(ifile, nlayer_1+1, pressure_1,   fileformat, 'pressure_1')
   call read_r8_array(ifile, nlayer_1,   avg_kernel_1, fileformat, 'avg_kernel_1')
   call read_r8_array(ifile, nlayer_1,   prior_1,      fileformat, 'prior_1')
   keyin = read_int_scalar(ifile, fileformat, 'keyin')
   
   counts1 = counts1 + 1
   key     = counts1
   
   if(counts1 > max_cris_nh3_obs) then
      write(string1, *)'Not enough space for cris nh3 obs.'
      write(string2, *)'Can only have max_cris_nh3_obs (currently ',max_cris_nh3_obs,')'
      call error_handler(E_ERR,'read_cris_nh3_cpsr',string1,source,text2=string2)
   endif
   
   call set_obs_def_cris_nh3_cpsr(key, pressure_1(1:nlayer_1+1), avg_kernel_1(1:nlayer_1), &
   prior_1(1:nlayer_1), nlayer_1, klev_1, kend_1)
   
   deallocate(pressure_1, avg_kernel_1, prior_1)
   
end subroutine read_cris_nh3_cpsr

!-------------------------------------------------------------------------------

subroutine write_cris_nh3_cpsr(key, ifile, fform)

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
   call write_r8_array(  ifile, nlayer(key)+1,  pressure(key,:), fileformat,'pressure')
   call write_r8_array(  ifile, nlayer(key),  avg_kernel(key,:), fileformat,'avg_kernel')
   call write_r8_array(  ifile, nlayer(key),       prior(key,:), fileformat,'prior')
   call write_int_scalar(ifile,                             key, fileformat,'key')
   
end subroutine write_cris_nh3_cpsr

!-------------------------------------------------------------------------------

subroutine interactive_cris_nh3_cpsr(key)

   integer, intent(out) :: key
   
   if ( .not. module_initialized ) call initialize_module

! STOP because routine is not finished.
   write(string1,*)'interactive_cris_nh3_cpsr not yet working.'
   call error_handler(E_ERR, 'interactive_cris_nh3_cpsr', string1, source)
   
   if(num_cris_nh3_obs >= max_cris_nh3_obs) then
      write(string1, *)'Not enough space for an cris nh3 obs.'
      write(string2, *)'Can only have max_cris_nh3_obs (currently ',max_cris_nh3_obs,')'
      call error_handler(E_ERR, 'interactive_cris_nh3_cpsr', string1, &
                 source, text2=string2)
   endif
   
! Increment the index
   num_cris_nh3_obs = num_cris_nh3_obs + 1
   key            = num_cris_nh3_obs

! Otherwise, prompt for input for the three required beasts

   write(*, *) 'Creating an interactive_cris_nh3_cpsr observation'
   write(*, *) 'This featue is not setup '

end subroutine interactive_cris_nh3_cpsr

!-------------------------------------------------------------------------------

subroutine get_expected_cris_nh3_cpsr(state_handle, ens_size, location, key, obs_time, expct_val, istatus)

   type(ensemble_type), intent(in)  :: state_handle
   type(location_type), intent(in)  :: location
   integer,             intent(in)  :: ens_size
   integer,             intent(in)  :: key
   type(time_type),     intent(in)  :: obs_time
   integer,             intent(out) :: istatus(:)
   real(r8),            intent(out) :: expct_val(:)
   
   character(len=*), parameter :: routine = 'get_expected_cris_nh3_cpsr'
   character(len=120)          :: data_file
   character(len=*),parameter  :: fld = 'NH3_VMR_inst'
   type(location_type) :: loc2
   
   integer :: layer_cris,level_cris, klev_cris, kend_cris
   integer :: layer_mdl,level_mdl
   integer :: k,kk,imem,imemm,flg
   integer :: interp_new
   integer :: icnt,ncnt
   integer :: date_obs,datesec_obs
   integer, dimension(ens_size) :: zstatus,kbnd_1,kbnd_n,kstart
   
   real(r8) :: eps, AvogN, Rd, Ru, Cp, grav, msq2cmsq
   real(r8) :: missing,nh3_min,tmp_max
   real(r8) :: level,del_prs,prior_term
   real(r8) :: tmp_vir_k, tmp_vir_kp
   real(r8) :: mloc(3),obs_prs
   real(r8) :: nh3_val_conv, VMR_conv
   real(r8) :: up_wt,dw_wt,tl_wt,lnpr_mid
   real(r8) :: lon_obs,lat_obs,pi,rad2deg

   real(r8), dimension(ens_size) :: nh3_mdl_1, tmp_mdl_1, qmr_mdl_1, prs_mdl_1
   real(r8), dimension(ens_size) :: nh3_mdl_n, tmp_mdl_n, qmr_mdl_n, prs_mdl_n
   real(r8), dimension(ens_size) :: nh3_mdl_tmp, tmp_mdl_tmp, qmr_mdl_tmp, prs_mdl_tmp
   real(r8), dimension(ens_size) :: prs_sfc,rec_nh3_val,rec_tmp_val,rec_qmr_val
   
   real(r8), allocatable, dimension(:)   :: thick, prs_cris, prs_cris_mem
   real(r8), allocatable, dimension(:,:) :: nh3_val, tmp_val, qmr_val
   real(r8), allocatable, dimension(:)   :: nh3_prf_mdl,tmp_prf_mdl,qmr_prf_mdl
   real(r8), allocatable, dimension(:)   :: prs_cris_top   
   logical  :: return_now,nh3_return_now,tmp_return_now,qmr_return_now
   
   if ( .not. module_initialized ) call initialize_module
   
   pi       = 4.*atan(1.)
   rad2deg  = 360./(2.*pi)
   eps      = 0.61_r8
   Rd       = 287.05_r8     ! J/(mole-kg)
   Ru       = 8.316_r8      ! J/(mole-kg)
   Cp       = 1006.0        ! J/kg/K
   grav     = 9.8_r8        ! m/s^2
   nh3_min   = 1.e-6_r8
   msq2cmsq = 1.e4_r8
   AvogN    = 6.02214e23_r8
   missing  = -888888_r8
   tmp_max  = 600.
   del_prs  = 5000.
   VMR_conv = 28.9644/47.9982
! 
! WACCM - MMR
! WRFChem - VMR ppmv
! CRIS NH3 - DU   
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
   if(use_log_nh3) then
      nh3_min = log(nh3_min)
   endif
   
! Assign vertical grid information (CRIS NH3 grid is top to bottom)

   layer_cris = nlayer(key)
   level_cris = nlayer(key)+1
   klev_cris  = klev(key)
   kend_cris  = kend(key)
   layer_mdl   = nlayer_model
   level_mdl   = nlayer_model+1

   allocate(prs_cris(level_cris))
   allocate(prs_cris_mem(level_cris))
   prs_cris(1:level_cris)=pressure(key,1:level_cris)

! Get location infomation

   mloc = get_location(location)
   
   if (mloc(2) >  90.0_r8) then
      mloc(2) =  90.0_r8
   elseif (mloc(2) < -90.0_r8) then
      mloc(2) = -90.0_r8
   endif
   obs_prs=mloc(3)
!   write(string1, *) 'APM: observation ',key, ' lon ',mloc(1),' lat ',mloc(2)
!   call error_handler(E_MSG, routine, string1, source)
!
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

   nh3_mdl_tmp(:)=missing_r8
   tmp_mdl_tmp(:)=missing_r8
   qmr_mdl_tmp(:)=missing_r8
   prs_mdl_tmp(:)=missing_r8

   kbnd_1(:)=1
   do k=1,layer_mdl
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
      call interpolate(state_handle, ens_size, loc2, QTY_NH3, nh3_mdl_tmp, zstatus) ! ppmv 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_tmp, zstatus) ! K 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_tmp, zstatus) ! kg / kg 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_tmp, zstatus) ! Pa
!
      interp_new=0
      do imem=1,ens_size
         if(nh3_mdl_tmp(imem).eq.missing_r8 .or. tmp_mdl_tmp(imem).eq.missing_r8 .or. &
         qmr_mdl_tmp(imem).eq.missing_r8 .or. prs_mdl_tmp(imem).eq.missing_r8) then
            interp_new=1
         else
            kbnd_1(imem)=k
            nh3_mdl_1(imem)=nh3_mdl_tmp(imem)
            tmp_mdl_1(imem)=tmp_mdl_tmp(imem)
            qmr_mdl_1(imem)=qmr_mdl_tmp(imem)
            prs_mdl_1(imem)=prs_mdl_tmp(imem)
         endif
      enddo
      if(interp_new.eq.0) exit
   enddo
!
! Sometimes the WRF-Chem surface pressure is greater than the
! first model level pressure. This fixes that problem.   
   do imem=1,ens_size
      if(prs_sfc(imem).lt.prs_mdl_1(imem)) prs_mdl_1(imem)=prs_sfc(imem)
   enddo

!   write(string1, *) 'APM: nh3 lower bound ',key,nh3_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: tmp lower bound ',key,tmp_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: qmr lower bound ',key,qmr_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: prs lower bound ',key,prs_mdl_1
!   call error_handler(E_MSG, routine, string1, source)

! pressure at model top (Pa)

   nh3_mdl_tmp(:)=missing_r8
   tmp_mdl_tmp(:)=missing_r8
   qmr_mdl_tmp(:)=missing_r8
   prs_mdl_tmp(:)=missing_r8

   kbnd_n(:)=layer_mdl
   do k=layer_mdl,1,-1
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
      call interpolate(state_handle, ens_size, loc2, QTY_NH3, nh3_mdl_tmp, zstatus) 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_tmp, &
      zstatus) 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_tmp, &
      zstatus) 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_tmp, &
      zstatus) 
!
      interp_new=0
      do imem=1,ens_size
         if(nh3_mdl_tmp(imem).eq.missing_r8 .or. tmp_mdl_tmp(imem).eq.missing_r8 .or. &
         qmr_mdl_tmp(imem).eq.missing_r8 .or. prs_mdl_tmp(imem).eq.missing_r8) then
            interp_new=1
         else
            kbnd_n(imem)=k
            nh3_mdl_n(imem)=nh3_mdl_tmp(imem)
            tmp_mdl_n(imem)=tmp_mdl_tmp(imem)
            qmr_mdl_n(imem)=qmr_mdl_tmp(imem)
            prs_mdl_n(imem)=prs_mdl_tmp(imem)
         endif
      enddo
      if(interp_new.eq.0) exit
   enddo

!   write(string1, *) 'APM: nh3 upper bound ',key,nh3_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: tmp upper bound ',key,tmp_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: qmr upper bound ',key,qmr_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: prs upper bound ',key,prs_mdl_n
!   call error_handler(E_MSG, routine, string1, source)

! Get cpsrs at CRIS pressure levels

   allocate(nh3_val(ens_size,level_cris))
   allocate(tmp_val(ens_size,level_cris))
   allocate(qmr_val(ens_size,level_cris))

   do k=1,level_cris
      zstatus=0
      loc2 = set_location(mloc(1), mloc(2), prs_cris(k), VERTISPRESSURE)
      call interpolate(state_handle, ens_size, loc2, QTY_NH3, nh3_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_val(:,k), zstatus)  
!
! Correcting for expected failures near the surface
      do imem=1,ens_size
         if (prs_cris(k).ge.prs_mdl_1(imem)) then
            nh3_val(imem,k) = nh3_mdl_1(imem)
            tmp_val(imem,k) = tmp_mdl_1(imem)
            qmr_val(imem,k) = qmr_mdl_1(imem)
         endif
!
! Correcting for expected failures near the top
         if (prs_cris(k).le.prs_mdl_n(imem)) then
            nh3_val(imem,k) = nh3_mdl_n(imem)
            tmp_val(imem,k) = tmp_mdl_n(imem)
            qmr_val(imem,k) = qmr_mdl_n(imem)
         endif
      enddo

!      write(string1, *)'APM: nh3 ',key,k,nh3_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: tmp ',key,k,tmp_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: qmr ',key,k,qmr_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)

! Check data for missing values      
      do imem=1,ens_size
         if(nh3_val(imem,k).eq.missing_r8 .or. tmp_val(imem,k).eq.missing_r8 .or. &
         qmr_val(imem,k).eq.missing_r8) then
            zstatus(:)=20
            expct_val(:)=missing_r8
            write(string1, *) 'APM: Model cpsr data has missing values for obs, level ',key,k
            call error_handler(E_ALLMSG, routine, string1, source)
            call track_status(ens_size, zstatus, expct_val, istatus, return_now)
            do imemm=1,ens_size
               write(string1, *) &
               'APM: Model cpsr values: nh3,tmp,qmr',key,imem,k,nh3_val(imemm,k), &
               tmp_val(imemm,k),qmr_val(imemm,k)     
               call error_handler(E_ALLMSG, routine, string1, source)
            enddo
            return
         endif
      enddo
!
! Convert units for nh3 from ppmv
      nh3_val(:,k) = nh3_val(:,k) * 1.e-6_r8
   enddo
!
! Use large scale ozone data above the regional model top
! CRIS vertical is from top to bottom   
   kstart(:)=-1
   do imem=1,ens_size
      if (prs_cris(1).lt.prs_mdl_n(imem)) then
         do k=1,level_cris
            if (prs_cris(k).ge.prs_mdl_n(imem)) then
               kstart(imem)=k-1
               exit
            endif
         enddo
         ncnt=kstart(imem)
         allocate(prs_cris_top(ncnt))
         allocate(nh3_prf_mdl(ncnt),tmp_prf_mdl(ncnt),qmr_prf_mdl(ncnt))
         do k=1,kstart(imem)
            prs_cris_top(k)=prs_cris(k)
         enddo
         prs_cris_top(:)=prs_cris_top(:)/100.
!
         lon_obs=mloc(1)/rad2deg
         lat_obs=mloc(2)/rad2deg
         call get_time(obs_time,datesec_obs,date_obs)
!
         data_file=trim(upper_data_file)
         model=trim(upper_data_model)
         call get_upper_bdy_fld(fld,model,data_file,ls_chem_dx,ls_chem_dy, &
         ls_chem_dz,ls_chem_dt,lon_obs,lat_obs,prs_cris_top, &
         ncnt,nh3_prf_mdl,tmp_prf_mdl,qmr_prf_mdl,date_obs,datesec_obs)
!
! Impose ensemble perturbations from level kstart(imem)+1      
         nh3_prf_mdl(:)=nh3_prf_mdl(:)*VMR_conv
         do k=1,kstart(imem) 
            nh3_val(imem,k)=nh3_prf_mdl(k)*nh3_val(imem,kstart(imem)+1)/ &
            (sum(nh3_val(:,kstart(imem)+1))/real(ens_size))
            tmp_val(imem,k)=tmp_prf_mdl(k)*tmp_val(imem,kstart(imem)+1)/ &
            (sum(tmp_val(:,kstart(imem)+1))/real(ens_size))
            qmr_val(imem,k)=qmr_prf_mdl(k)*qmr_val(imem,kstart(imem)+1)/ &
            (sum(qmr_val(:,kstart(imem)+1))/real(ens_size))
         enddo
         deallocate(prs_cris_top)
         deallocate(nh3_prf_mdl,tmp_prf_mdl,qmr_prf_mdl)
      endif             
   enddo
!
! Check full cpsr for negative values
   do imem=1,ens_size
      flg=0
      do k=1,level_cris   
         if(nh3_val(imem,k).lt.0. .or. tmp_val(imem,k).lt.0. .or. &
         qmr_val(imem,k).lt.0.) then
            flg=1   
            write(string1, *) &
            'APM: Recentered full cpsr has negative values for key,imem ',key,imem
            call error_handler(E_ALLMSG, routine, string1, source)
            call track_status(ens_size, zstatus, expct_val, istatus, return_now)
            if(k.le.kstart(imem)) then
               do kk=1,level_cris
                  write(string1, *) &
                  'APM: prs, nh3, tmp, qmr',key,imem,prs_cris(kk),nh3_val(imem,kk), &
                  tmp_val(imem,kk),qmr_val(imem,kk)
                  call error_handler(E_ALLMSG, routine, string1, source)
               enddo
               exit
            endif
         endif
      enddo
      if(flg.eq.1) exit
   enddo

   istatus(:)=0
   zstatus(:)=0.
   expct_val(:)=0.0
   allocate(thick(layer_cris))

   do imem=1,ens_size
! Adjust the CRIS pressure for WRF-Chem lower/upper boudary pressure
! (CRIS NH3 vertical grid is top to bottom)
      prs_cris_mem(:)=prs_cris(:)
      if (prs_sfc(imem).gt.prs_cris_mem(level_cris)) then
         prs_cris_mem(level_cris)=prs_sfc(imem)
      endif   

! Calculate the thicknesses

      thick(:)=0.
      do k=1,layer_cris
         lnpr_mid=(log(prs_cris_mem(k+1))+log(prs_cris_mem(k)))/2.
         up_wt=log(prs_cris_mem(k+1))-lnpr_mid
         dw_wt=lnpr_mid-log(prs_cris_mem(k))
         tl_wt=up_wt+dw_wt
         tmp_vir_k  = (1.0_r8 + eps*qmr_val(imem,k))*tmp_val(imem,k)
         tmp_vir_kp = (1.0_r8 + eps*qmr_val(imem,k+1))*tmp_val(imem,k+1)
         thick(k)   = Rd*(dw_wt*tmp_vir_kp + up_wt*tmp_vir_k)/tl_wt/grav* &
         log(prs_cris_mem(k+1)/prs_cris_mem(k))
      enddo

! Process the vertical summation

      do k=1,layer_cris
         if(prior(key,k).lt.0.) then
            write(string1, *) &
            'APM: CRIS Prior is negative. Level may be below surface. Key,Layer: ',key,k
            call error_handler(E_MSG, routine, string1, source)
            write(string1, *) &
            'APM: Key ',key,' Prior: ',prior(key,k),' Avgk ',avg_kernel(key,k)
            call error_handler(E_MSG, routine, string1, source)
            cycle
         endif
!         
         lnpr_mid=(log(prs_cris_mem(k+1))+log(prs_cris_mem(k)))/2.
         up_wt=log(prs_cris_mem(k+1))-lnpr_mid
         dw_wt=lnpr_mid-log(prs_cris_mem(k))
         tl_wt=up_wt+dw_wt
   
! Convert from VMR to molar density (mol/m^3)
         if(use_log_nh3) then
            nh3_val_conv = (dw_wt*exp(nh3_val(imem,k+1))+up_wt*exp(nh3_val(imem,k)))/tl_wt * &
            (dw_wt*prs_cris_mem(k+1)+up_wt*prs_cris_mem(k)) / &
            (Ru*(dw_wt*tmp_val(imem,k+1)+up_wt*tmp_val(imem,k)))
         else
            nh3_val_conv = (dw_wt*nh3_val(imem,k+1)+up_wt*nh3_val(imem,k))/tl_wt * &
            (dw_wt*prs_cris_mem(k+1)+up_wt*prs_cris_mem(k)) / &
            (Ru*(dw_wt*tmp_val(imem,k+1)+up_wt*tmp_val(imem,k)))
         endif
 
! Get expected observation

         prior_term=-1.*avg_kernel(key,k)
         if(k.eq.klev_cris) prior_term=(1.0_r8 - avg_kernel(key,k)) 

         expct_val(imem) = expct_val(imem) + thick(k) * nh3_val_conv * &
         avg_kernel(key,k) + prior_term * prior(key,k)
         
!         write(string1, *) &
!         'APM: Mem ',imem,' Key ',key,' Expct Val Terms: prs ',k, &
!         (prs_cris_mem(k)+prs_cris_mem(k+1))/2.,' expct val ',expct_val(imem), &
!         'avgk*thick*nh3_conv ', avg_kernel(key,k)*thick(k)*nh3_val_conv, &
!         'prior_term*prior ', prior_term*prior(key,k)
!         call error_handler(E_MSG, routine, string1, source)
      enddo

! call exit_all(-77)

      if(expct_val(imem).lt.0.) then
         write(string1, *) &
         'APM: Member ',imem,'Key, Final Value ',key,expct_val(imem)
         call error_handler(E_ALLMSG, routine, string1, source)
      endif
!      
      if(isnan(expct_val(imem))) then
         zstatus(imem)=20
         expct_val(:)=missing_r8
         write(string1, *) &
         'APM NOTICE: CRIS NH3 expected value is NaN '
         call error_handler(E_MSG, routine, string1, source)
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
!
      if(expct_val(imem).lt.0) then
         zstatus(imem)=20
         expct_val(:)=missing_r8
         write(string1, *) &
         'APM NOTICE: CRIS NH3 expected value is negative '
         call error_handler(E_MSG, routine, string1, source)
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
   enddo

! Clean up and return
   deallocate(nh3_val, tmp_val, qmr_val)
   deallocate(thick)
   deallocate(prs_cris, prs_cris_mem)

end subroutine get_expected_cris_nh3_cpsr

!-------------------------------------------------------------------------------

subroutine set_obs_def_cris_nh3_cpsr(key, nh3_pressure, nh3_avg_kernel, nh3_prior, &
nh3_nlayer, nh3_klev, nh3_kend)

   integer,                           intent(in)   :: key, nh3_nlayer, nh3_klev, nh3_kend
   real(r8), dimension(nh3_nlayer+1),  intent(in)   :: nh3_pressure
   real(r8), dimension(nh3_nlayer),    intent(in)   :: nh3_avg_kernel
   real(r8), dimension(nh3_nlayer),    intent(in)   :: nh3_prior
   
   if ( .not. module_initialized ) call initialize_module
   
   if(num_cris_nh3_obs >= max_cris_nh3_obs) then
      write(string1, *)'Not enough space for cris nh3 obs.'
      write(string2, *)'Can only have max_cris_nh3_obs (currently ',max_cris_nh3_obs,')'
      call error_handler(E_ERR,'set_obs_def_cris_nh3_cpsr',string1,source,revision, &
      revdate,text2=string2)
   endif
   
   nlayer(key) = nh3_nlayer
   klev(key) = nh3_klev
   kend(key) = nh3_kend
   pressure(key,1:nh3_nlayer+1) = nh3_pressure(1:nh3_nlayer+1)
   avg_kernel(key,1:nh3_nlayer) = nh3_avg_kernel(1:nh3_nlayer)
   prior(key,1:nh3_nlayer)      = nh3_prior(1:nh3_nlayer)
   
end subroutine set_obs_def_cris_nh3_cpsr

end module obs_def_cris_nh3_cpsr_mod

! END DART PREPROCESS MODULE CODE
