
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
! TES_CH4_PROFILE, QTY_CH4
! END DART PREPROCESS TYPE DEFINITIONS
!
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_tes_ch4_profile_mod, only : get_expected_tes_ch4_profile, &
!                                  read_tes_ch4_profile, &
!                                  write_tes_ch4_profile, &
!                                  interactive_tes_ch4_profile
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(TES_CH4_PROFILE)                                                           
!         call get_expected_tes_ch4_profile(state_handle, ens_size, location, obs_def%key, obs_time, expected_obs, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!
! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(TES_CH4_PROFILE)
!         call read_tes_ch4_profile(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(TES_CH4_PROFILE)
!         call write_tes_ch4_profile(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF
!
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(TES_CH4_PROFILE)
!         call interactive_tes_ch4_profile(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_tes_ch4_profile_mod

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
   
   use          obs_kind_mod, only : QTY_CH4, QTY_TEMPERATURE, QTY_SURFACE_PRESSURE, &
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

   public :: write_tes_ch4_profile, &
             read_tes_ch4_profile, &
          interactive_tes_ch4_profile, &
          get_expected_tes_ch4_profile, &
          set_obs_def_tes_ch4_profile

! Storage for the special information required for observations of this type
   integer, parameter    :: max_tes_ch4_obs = 10000000
   integer               :: num_tes_ch4_obs = 0
   integer,  allocatable :: nlayer(:)
   integer,  allocatable :: klev(:)
   integer,  allocatable :: kend(:)
   real(r8), allocatable :: pressure(:,:)
   real(r8), allocatable :: avg_kernel(:,:)
   real(r8), allocatable :: prior(:,:)

! version controlled file description for error handling, do not edit
   character(len=*), parameter :: source   = 'obs_def_tes_ch4_profile_mod.f90'
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
   integer :: nlayer_tes = -9999
   integer :: nlayer_tes_ch4_total_col = -9999
   integer :: nlayer_tes_ch4_trop_col = -9999
   integer :: nlayer_tes_ch4_profile = -9999
   
   namelist /obs_def_TES_CH4_nml/ upper_data_file, use_log_ch4, nlayer_model, &
   nlayer_tes_ch4_total_col, nlayer_tes_ch4_trop_col, nlayer_tes_ch4_profile, &
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
   call find_namelist_in_file("input.nml", "obs_def_TES_CH4_nml", iunit)
   read(iunit, nml = obs_def_TES_CH4_nml, iostat = rc)
   call check_namelist_read(iunit, rc, "obs_def_TES_CH4_nml")

! Record the namelist values
   if (do_nml_file()) write(nmlfileunit, nml=obs_def_TES_CH4_nml)
   if (do_nml_term()) write(     *     , nml=obs_def_TES_CH4_nml)

! Check for valid values
   nlayer_tes=nlayer_tes_ch4_profile

   if (nlayer_model < 1) then
      write(string1,*)'obs_def_TES_CH4_nml:nlayer_model must be > 0, it is ',nlayer_model
      call error_handler(E_ERR,'initialize_module',string1,source)
   endif
   
   if (nlayer_tes < 1) then
      write(string1,*)'obs_def_TES_CH4_nml:nlayer_tes must be > 0, it is ',nlayer_tes
      call error_handler(E_ERR,'initialize_module',string1,source)
   endif
   
   allocate(    nlayer(max_tes_ch4_obs))
   allocate(    klev(max_tes_ch4_obs))
   allocate(    kend(max_tes_ch4_obs))
   allocate(  pressure(max_tes_ch4_obs,nlayer_tes))
   allocate(avg_kernel(max_tes_ch4_obs,nlayer_tes))
   allocate(     prior(max_tes_ch4_obs,nlayer_tes))
   
end subroutine initialize_module

!-------------------------------------------------------------------------------

subroutine read_tes_ch4_profile(key, ifile, fform)

   integer,          intent(out)          :: key
   integer,          intent(in)           :: ifile
   character(len=*), intent(in), optional :: fform

! tesrary arrays to hold buffer till we decide if we have enough room

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
   
   if(counts1 > max_tes_ch4_obs) then
      write(string1, *)'Not enough space for tes ch4 obs.'
      write(string2, *)'Can only have max_tes_ch4_obs (currently ',max_tes_ch4_obs,')'
      call error_handler(E_ERR,'read_tes_ch4_profile',string1,source,text2=string2)
   endif
   
   call set_obs_def_tes_ch4_profile(key, pressure_1(1:nlayer_1), avg_kernel_1(1:nlayer_1), &
   prior_1(1:nlayer_1), nlayer_1, klev_1, kend_1)
   
   deallocate(pressure_1, avg_kernel_1, prior_1)
   
end subroutine read_tes_ch4_profile

!-------------------------------------------------------------------------------

subroutine write_tes_ch4_profile(key, ifile, fform)

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
   
end subroutine write_tes_ch4_profile

!-------------------------------------------------------------------------------

subroutine interactive_tes_ch4_profile(key)

   integer, intent(out) :: key
   
   if ( .not. module_initialized ) call initialize_module

! STOP because routine is not finished.
   write(string1,*)'interactive_tes_ch4_profile not yet working.'
   call error_handler(E_ERR, 'interactive_tes_ch4_profile', string1, source)
   
   if(num_tes_ch4_obs >= max_tes_ch4_obs) then
      write(string1, *)'Not enough space for an tes ch4 obs.'
      write(string2, *)'Can only have max_tes_ch4_obs (currently ',max_tes_ch4_obs,')'
      call error_handler(E_ERR, 'interactive_tes_ch4_profile', string1, &
                 source, text2=string2)
   endif
   
! Increment the index
   num_tes_ch4_obs = num_tes_ch4_obs + 1
   key            = num_tes_ch4_obs

! Otherwise, prompt for input for the three required beasts

   write(*, *) 'Creating an interactive_tes_ch4_profile observation'
   write(*, *) 'This featue is not setup '

end subroutine interactive_tes_ch4_profile

!-------------------------------------------------------------------------------

subroutine get_expected_tes_ch4_profile(state_handle, ens_size, location, key, obs_time, expct_val, istatus)

   type(ensemble_type), intent(in)  :: state_handle
   type(location_type), intent(in)  :: location
   integer,             intent(in)  :: ens_size
   integer,             intent(in)  :: key
   type(time_type),     intent(in)  :: obs_time
   integer,             intent(out) :: istatus(:)
   real(r8),            intent(out) :: expct_val(:)
   
   character(len=*), parameter :: routine = 'get_expected_tes_ch4_profile'
   character(len=120)          :: data_file
   character(len=*),parameter  :: fld = 'CH4_VMR_inst'
   type(location_type) :: loc2
   
   integer :: layer_tes,level_tes
   integer :: layer_mdl,level_mdl
   integer :: k,kk,imem,imemm,flg
   integer :: interp_new
   integer :: icnt,ncnt,kstart,klev_tes
   integer :: date_obs,datesec_obs
   integer, dimension(ens_size) :: zstatus,kbnd_1,kbnd_n
   
   real(r8) :: eps, AvogN, Rd, Ru, Cp, grav, msq2cmsq
   real(r8) :: missing,ch4_min,tmp_max
   real(r8) :: level,del_prs,prior_term
   real(r8) :: tmp_vir_k, tmp_vir_kp
   real(r8) :: mloc(3),obs_prs
   real(r8) :: ch4_val_conv, VMR_conv
   real(r8) :: up_wt,dw_wt,tl_wt,lnpr_mid
   real(r8) :: lon_obs,lat_obs,pi,rad2deg

   real(r8), dimension(ens_size) :: ch4_mdl_1, tmp_mdl_1, qmr_mdl_1, prs_mdl_1
   real(r8), dimension(ens_size) :: ch4_mdl_n, tmp_mdl_n, qmr_mdl_n, prs_mdl_n
   real(r8), dimension(ens_size) :: prs_sfc
   
   real(r8), allocatable, dimension(:)   :: thick, prs_tes, prs_tes_mem
   real(r8), allocatable, dimension(:,:) :: ch4_val, tmp_val, qmr_val
   logical  :: return_now,ch4_return_now,tmp_return_now,qmr_return_now
!
! Upper BC variables
   real(r8), allocatable, dimension(:)   :: ch4_prf_mdl,tmp_prf_mdl,qmr_prf_mdl
   real(r8), allocatable, dimension(:)   :: prs_tes_top   
   
   if ( .not. module_initialized ) call initialize_module
   
   pi       = 4.*atan(1.)
   rad2deg  = 360./(2.*pi)
   eps      = 0.61_r8
   Rd       = 287.05_r8     ! J/(mole-kg)
   Ru       = 8.316_r8      ! J/(mole-kg)
   Cp       = 1006.0        ! J/kg/K
   grav     = 9.8_r8        ! m/s^2
   ch4_min   = 1.e-6_r8
   msq2cmsq = 1.e4_r8
   AvogN    = 6.02214e23_r8
   missing  = -888888_r8
   tmp_max  = 600.
   del_prs  = 5000.
   VMR_conv = 28.9644/47.9982
! 
! WACCM - MMR
! WRFChem - VMR ppmv
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
   if(use_log_ch4) then
      ch4_min = log(ch4_min)
   endif
   
! Assign vertical grid information (TES CH4 grid is top to bottom)

   layer_tes = nlayer(key)
   level_tes = nlayer(key)+1
   klev_tes  = klev(key)
   layer_mdl   = nlayer_model
   level_mdl   = nlayer_model+1

   allocate(prs_tes(level_tes))
   prs_tes(1:layer_tes)=pressure(key,1:layer_tes)

! Get location infomation

   mloc = get_location(location)
   
   if (mloc(2) >  90.0_r8) then
      mloc(2) =  90.0_r8
   elseif (mloc(2) < -90.0_r8) then
      mloc(2) = -90.0_r8
   endif
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

   ch4_mdl_1(:)=missing_r8
   tmp_mdl_1(:)=missing_r8
   qmr_mdl_1(:)=missing_r8
   prs_mdl_1(:)=missing_r8

   do k=1,layer_mdl
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
      call interpolate(state_handle, ens_size, loc2, QTY_CH4, ch4_mdl_1, zstatus) ! ppmv 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_1, zstatus) ! K 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_1, zstatus) ! kg / kg 
      zstatus(:)=0
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
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_n, zstatus) ! K 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_n, zstatus) ! kg / kg 
      zstatus(:)=0
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

! Get profiles at TES pressure levels

   allocate(ch4_val(ens_size,layer_tes))
   allocate(tmp_val(ens_size,layer_tes))
   allocate(qmr_val(ens_size,layer_tes))

   do k=1,layer_tes
      zstatus=0
      loc2 = set_location(mloc(1), mloc(2), prs_tes(k), VERTISPRESSURE)
      call interpolate(state_handle, ens_size, loc2, QTY_CH4, ch4_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_val(:,k), zstatus)  
!
! Correcting for expected failures near the surface
      do imem=1,ens_size
         if (prs_tes(k).ge.prs_mdl_1(imem)) then
            ch4_val(imem,k) = ch4_mdl_1(imem)
            tmp_val(imem,k) = tmp_mdl_1(imem)
            qmr_val(imem,k) = qmr_mdl_1(imem)
         endif
!
! Correcting for expected failures near the top
         if (prs_tes(k).le.prs_mdl_n(imem)) then
            ch4_val(imem,k) = ch4_mdl_n(imem)
            tmp_val(imem,k) = tmp_mdl_n(imem)
            qmr_val(imem,k) = qmr_mdl_n(imem)
         endif
      enddo

!      write(string1, *)'APM: ch4 ',key,k,ch4_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: tmp ',key,k,tmp_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: qmr ',key,k,qmr_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!
! Convert units for ch4 from ppmv
      ch4_val(:,k) = ch4_val(:,k) * 1.e-6_r8
   enddo
   ch4_mdl_1(:) = ch4_mdl_1(:) * 1.e-6_r8
   ch4_mdl_n(:) = ch4_mdl_n(:) * 1.e-6_r8
!
! Use large scale ch4 data above the regional model top
! TES vertical is from bottom to top   
   kstart=-1
   do imem=1,ens_size
      if (prs_tes(layer_tes).lt.prs_mdl_n(imem)) then
         do k=1,layer_tes
            if (prs_tes(k).le.prs_mdl_n(imem)) then
               kstart=k
               exit
            endif
         enddo
         ncnt=layer_tes-kstart+1
         allocate(prs_tes_top(ncnt))
         allocate(ch4_prf_mdl(ncnt),tmp_prf_mdl(ncnt),qmr_prf_mdl(ncnt))
         do k=kstart,layer_tes
            prs_tes_top(k-kstart+1)=prs_tes(k)
         enddo
         prs_tes_top(:)=prs_tes_top(:)/100.
!
         lon_obs=mloc(1)/rad2deg
         lat_obs=mloc(2)/rad2deg
         call get_time(obs_time,datesec_obs,date_obs)
!
         data_file=trim(upper_data_file)
         model=trim(upper_data_model)
         call get_upper_bdy_fld(fld,model,data_file,ls_chem_dx,ls_chem_dy, &
         ls_chem_dz,ls_chem_dt,lon_obs,lat_obs,prs_tes_top, &
         ncnt,ch4_prf_mdl,tmp_prf_mdl,qmr_prf_mdl,date_obs,datesec_obs)
!
! Impose ensemble perturbations from level kstart(imem)-1      
         do k=kstart,layer_tes
            kk=k-kstart+1
            ch4_val(imem,k)=ch4_prf_mdl(kk)*ch4_val(imem,kstart-1)/ &
            (sum(ch4_val(:,kstart-1))/real(ens_size))
            tmp_val(imem,k)=tmp_prf_mdl(kk)*tmp_val(imem,kstart-1)/ &
            (sum(tmp_val(:,kstart-1))/real(ens_size))
            qmr_val(imem,k)=qmr_prf_mdl(kk)*qmr_val(imem,kstart-1)/ &
            (sum(qmr_val(:,kstart-1))/real(ens_size))
         enddo
         deallocate(prs_tes_top)
         deallocate(ch4_prf_mdl,tmp_prf_mdl,qmr_prf_mdl)
      endif             
   enddo
!
! Check full profile for negative values

!   do imem=1,1
!      do k=1,layer_tes
!         write(string1, *) &
!         'APM: prs, ch4, tmp, qmr ',k,prs_tes(k),ch4_val(imem,k), &
!         tmp_val(imem,k),qmr_val(imem,k)
!         call error_handler(E_MSG, routine, string1, source)
!      enddo
!   enddo      
!
   do imem=1,ens_size
      flg=0
      do k=1,layer_tes   
         if(ch4_val(imem,k).lt.0. .or. tmp_val(imem,k).lt.0. .or. &
         qmr_val(imem,k).lt.0.) then
            flg=1   
            write(string1, *) &
            'APM: Recentered full profile has negative values for key,imem ',key,imem
            call error_handler(E_ALLMSG, routine, string1, source)
         endif
      enddo
      if(flg.eq.1) then
         zstatus(:)=20
         expct_val(:)=missing_r8
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
   enddo
!
! Calculate the expected retrievals
   istatus(:)=0
   zstatus(:)=0.
   expct_val(:)=0.0
!      
! Process the vertical summation
   do imem=1,ens_size
      do k=1,layer_tes
         if(prior(key,k).lt.0.) then
!            write(string1, *) &
!            'APM: TES Prior is negative. Level may be below surface. Key,Layer: ',key,k
!            call error_handler(E_MSG, routine, string1, source)
            cycle
         endif
!
! Get expected observation
         prior_term=-1.*avg_kernel(key,k)
         if(k.eq.klev_tes) prior_term=(1.0_r8 - avg_kernel(key,k)) 

         expct_val(imem) = expct_val(imem) + log(ch4_val(imem,k)) * &
         avg_kernel(key,k) + prior_term * log(prior(key,k))

!         write(string1, *) 'APM: exp_val, ch4, avgk, prior_trm, prior',imem,k, &
!         expct_val(imem),ch4_val(imem,k),avg_kernel(key,k),prior_term,prior(key,k)
!         call error_handler(E_MSG, routine, string1, source)
      enddo

      expct_val(imem)=exp(expct_val(imem))      
!      write(string1, *) 'APM: Finished vertical summation loop ',key,imem
!      call error_handler(E_MSG, routine, string1, source)
      
      if(isnan(expct_val(imem))) then
         zstatus(imem)=20
         expct_val(:)=missing_r8
!         write(string1, *) &
!         'APM NOTICE: TES CH4 expected value is NaN'
!         call error_handler(E_ALLMSG, routine, string1, source)
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
   enddo

! Clean up and return
   deallocate(ch4_val, tmp_val, qmr_val)
   deallocate(prs_tes)

end subroutine get_expected_tes_ch4_profile

!-------------------------------------------------------------------------------

subroutine set_obs_def_tes_ch4_profile(key, ch4_pressure, ch4_avg_kernel, ch4_prior, &
ch4_nlayer, ch4_klev, ch4_kend)

   integer,                           intent(in)   :: key, ch4_nlayer, ch4_klev, ch4_kend
   real(r8), dimension(ch4_nlayer),  intent(in)   :: ch4_pressure
   real(r8), dimension(ch4_nlayer),    intent(in)   :: ch4_avg_kernel
   real(r8), dimension(ch4_nlayer),    intent(in)   :: ch4_prior
   
   if ( .not. module_initialized ) call initialize_module
   
   if(num_tes_ch4_obs >= max_tes_ch4_obs) then
      write(string1, *)'Not enough space for tes ch4 obs.'
      write(string2, *)'Can only have max_tes_ch4_obs (currently ',max_tes_ch4_obs,')'
      call error_handler(E_ERR,'set_obs_def_tes_ch4_profile',string1,source,revision, &
      revdate,text2=string2)
   endif
   
   nlayer(key) = ch4_nlayer
   klev(key) = ch4_klev
   kend(key) = ch4_kend
   pressure(key,1:ch4_nlayer) = ch4_pressure(1:ch4_nlayer)
   avg_kernel(key,1:ch4_nlayer) = ch4_avg_kernel(1:ch4_nlayer)
   prior(key,1:ch4_nlayer)      = ch4_prior(1:ch4_nlayer)
   
end subroutine set_obs_def_tes_ch4_profile

end module obs_def_tes_ch4_profile_mod

! END DART PREPROCESS MODULE CODE
