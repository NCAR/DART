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
! TROPOMI_O3_CPSR, QTY_O3
! END DART PREPROCESS TYPE DEFINITIONS
!
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_tropomi_o3_cpsr_mod, only : get_expected_tropomi_o3_cpsr, &
!                                  read_tropomi_o3_cpsr, &
!                                  write_tropomi_o3_cpsr, &
!                                  interactive_tropomi_o3_cpsr
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(TROPOMI_O3_CPSR)                                                           
!         call get_expected_tropomi_o3_cpsr(state_handle, ens_size, location, obs_def%key, obs_time, expected_obs, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!
! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(TROPOMI_O3_CPSR)
!         call read_tropomi_o3_cpsr(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(TROPOMI_O3_CPSR)
!         call write_tropomi_o3_cpsr(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF
!
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(TROPOMI_O3_CPSR)
!         call interactive_tropomi_o3_cpsr(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_tropomi_o3_cpsr_mod

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
   
   use          obs_kind_mod, only : QTY_O3, QTY_TEMPERATURE, QTY_SURFACE_PRESSURE, &
                                     QTY_PRESSURE, QTY_VAPOR_MIXING_RATIO
   
   use  ensemble_manager_mod, only : ensemble_type
   
   use obs_def_utilities_mod, only : track_status
   
   use      time_manager_mod, only : time_type, get_date, set_date, get_time, set_time
! get_date gets year, month, day, hour, minute, second from time_type
! get_time gets julian day and seconds from time_type
! set_date sets time_type from year, month, day, hour, minute, second
! set_time sets time_type from julian day and seconds
   implicit none
   private

   public :: write_tropomi_o3_cpsr, &
             read_tropomi_o3_cpsr, &
          interactive_tropomi_o3_cpsr, &
          get_expected_tropomi_o3_cpsr, &
          set_obs_def_tropomi_o3_cpsr

! Storage for the special information required for observations of this type
   integer, parameter    :: max_tropomi_o3_obs = 10000000
   integer               :: num_tropomi_o3_obs = 0
   integer,  allocatable :: nlayer(:)
   real(r8), allocatable :: pressure(:,:)
   real(r8), allocatable :: avg_kernel(:,:)
   real(r8), allocatable :: prior(:)

! version controlled file description for error handling, do not edit
   character(len=*), parameter :: source   = 'obs_def_tropomi_o3_cpsr_mod.f90'
   character(len=*), parameter :: revision = ''
   character(len=*), parameter :: revdate  = ''
   
   character(len=512) :: string1, string2
   character(len=200) :: upper_data_file
   character(len=200) :: upper_data_model
   character(len=200) :: model
   integer            :: ls_chem_dx, ls_chem_dy, ls_chem_dz, ls_chem_dt
   
   logical, save :: module_initialized = .false.

! Namelist with default values
   logical :: use_log_o3   = .false.
   integer :: nlayer_model = -9999
   integer :: nlayer_tropomi = -9999
   integer :: nlayer_tropomi_o3_total_col = -9999
   integer :: nlayer_tropomi_o3_trop_col = -9999
   integer :: nlayer_tropomi_o3_profile = -9999
   
   namelist /obs_def_TROPOMI_O3_nml/ upper_data_file, use_log_o3, nlayer_model, &
   nlayer_tropomi_o3_total_col, nlayer_tropomi_o3_trop_col, nlayer_tropomi_o3_profile, &
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
   call find_namelist_in_file("input.nml", "obs_def_TROPOMI_O3_nml", iunit)
   read(iunit, nml = obs_def_TROPOMI_O3_nml, iostat = rc)
   call check_namelist_read(iunit, rc, "obs_def_TROPOMI_O3_nml")

! Record the namelist values
   if (do_nml_file()) write(nmlfileunit, nml=obs_def_TROPOMI_O3_nml)
   if (do_nml_term()) write(     *     , nml=obs_def_TROPOMI_O3_nml)
   nlayer_tropomi=nlayer_tropomi_o3_profile
   
! Check for valid values

   if (nlayer_model < 1) then
      write(string1,*)'obs_def_TROPOMI_O3_nml:nlayer_model must be > 0, it is ',nlayer_model
      call error_handler(E_ERR,'initialize_module',string1,source)
   endif
   
   if (nlayer_tropomi < 1) then
      write(string1,*)'obs_def_TROPOMI_O3_nml:nlayer_tropomi must be > 0, it is ',nlayer_tropomi
      call error_handler(E_ERR,'initialize_module',string1,source)
   endif
   
   allocate(    nlayer(max_tropomi_o3_obs))
   allocate(  pressure(max_tropomi_o3_obs,nlayer_tropomi+1))
   allocate(avg_kernel(max_tropomi_o3_obs,nlayer_tropomi))
   allocate(     prior(max_tropomi_o3_obs))
   
end subroutine initialize_module

!-------------------------------------------------------------------------------

subroutine read_tropomi_o3_cpsr(key, ifile, fform)

   integer,          intent(out)          :: key
   integer,          intent(in)           :: ifile
   character(len=*), intent(in), optional :: fform

! tropomirary arrays to hold buffer till we decide if we have enough room

   integer               :: keyin
   integer               :: nlayer_1
   real(r8)              :: prior_1
   real(r8), allocatable :: pressure_1(:)
   real(r8), allocatable :: avg_kernel_1(:)
   character(len=32)     :: fileformat
   
   integer, SAVE :: counts1 = 0
   
   if ( .not. module_initialized ) call initialize_module
   
   fileformat = "ascii" 
   if(present(fform)) fileformat = adjustl(fform)
   
! Need to know how many layers for this one
   nlayer_1 = read_int_scalar( ifile, fileformat, 'nlayer_1')
   prior_1 = read_r8_scalar( ifile, fileformat, 'prior_1')
   
   allocate(  pressure_1(nlayer_1+1))
   allocate(avg_kernel_1(nlayer_1))   

   call read_r8_array(ifile, nlayer_1+1, pressure_1,   fileformat, 'pressure_1')
   call read_r8_array(ifile, nlayer_1,   avg_kernel_1, fileformat, 'avg_kernel_1')
   keyin = read_int_scalar(ifile, fileformat, 'nlayer_1')
   
   counts1 = counts1 + 1
   key     = counts1
   
   if(counts1 > max_tropomi_o3_obs) then
      write(string1, *)'Not enough space for tropomi o3 obs.'
      write(string2, *)'Can only have max_tropomi_o3_obs (currently ',max_tropomi_o3_obs,')'
      call error_handler(E_ERR,'read_tropomi_o3_cpsr',string1,source,text2=string2)
   endif
   
   call set_obs_def_tropomi_o3_cpsr(key, pressure_1, avg_kernel_1, prior_1, nlayer_1)
   
   deallocate(pressure_1, avg_kernel_1)
   
end subroutine read_tropomi_o3_cpsr

!-------------------------------------------------------------------------------

subroutine write_tropomi_o3_cpsr(key, ifile, fform)

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
   call write_r8_scalar(ifile,                     prior(key), fileformat,'prior')
   call write_r8_array(  ifile, nlayer(key)+1,  pressure(key,:), fileformat,'pressure')
   call write_r8_array(  ifile, nlayer(key),  avg_kernel(key,:), fileformat,'avg_kernel')
   call write_int_scalar(ifile,                             key, fileformat,'key')
   
end subroutine write_tropomi_o3_cpsr

!-------------------------------------------------------------------------------

subroutine interactive_tropomi_o3_cpsr(key)

   integer, intent(out) :: key
   
   if ( .not. module_initialized ) call initialize_module

! STOP because routine is not finished.
   write(string1,*)'interactive_tropomi_o3_cpsr not yet working.'
   call error_handler(E_ERR, 'interactive_tropomi_o3_cpsr', string1, source)
   
   if(num_tropomi_o3_obs >= max_tropomi_o3_obs) then
      write(string1, *)'Not enough space for an tropomi o3 obs.'
      write(string2, *)'Can only have max_tropomi_o3_obs (currently ',max_tropomi_o3_obs,')'
      call error_handler(E_ERR, 'interactive_tropomi_o3_cpsr', string1, &
                 source, text2=string2)
   endif
   
! Increment the index
   num_tropomi_o3_obs = num_tropomi_o3_obs + 1
   key            = num_tropomi_o3_obs

! Otherwise, prompt for input for the three required beasts

   write(*, *) 'Creating an interactive_tropomi_o3_cpsr observation'
   write(*, *) 'This featue is not setup '

end subroutine interactive_tropomi_o3_cpsr

!-------------------------------------------------------------------------------

subroutine get_expected_tropomi_o3_cpsr(state_handle, ens_size, location, key, obs_time, expct_val, istatus)

   type(ensemble_type), intent(in)  :: state_handle
   type(location_type), intent(in)  :: location
   integer,             intent(in)  :: ens_size
   integer,             intent(in)  :: key
   type(time_type),     intent(in)  :: obs_time
   integer,             intent(out) :: istatus(:)
   real(r8),            intent(out) :: expct_val(:)
   
   character(len=*), parameter :: routine = 'get_expected_tropomi_o3_cpsr'
   character(len=120)          :: data_file
   character(len=*),parameter  :: fld = 'O3'
   type(location_type) :: loc2
   
   integer :: layer_tropomi,level_tropomi, kend_tropomi
   integer :: layer_mdl,level_mdl
   integer :: k,kk,imem
   integer :: interp_new
   integer :: icnt,ncnt
   integer :: date_obs,datesec_obs
   integer, dimension(ens_size) :: zstatus, kstart
   
   real(r8) :: eps, AvogN, Rd, Ru, Cp, grav, msq2cmsq, molec2du
   real(r8) :: missing,o3_min,tmp_max
   real(r8) :: level,del_prs
   real(r8) :: tmp_vir_k, tmp_vir_kp
   real(r8) :: mloc(3)
   real(r8) :: o3_val_conv, VMR_conv
   real(r8) :: up_wt,dw_wt,tl_wt,lnpr_mid
   real(r8) :: lon_obs,lat_obs,pi,rad2deg
   real(r8) :: ensavg_o3,ensavg_tmp,ensavg_qmr
   real(r8) :: fac_o3,fac_tmp,fac_qmr

   real(r8), dimension(ens_size) :: o3_mdl_1, tmp_mdl_1, qmr_mdl_1, prs_mdl_1
   real(r8), dimension(ens_size) :: o3_mdl_n, tmp_mdl_n, qmr_mdl_n, prs_mdl_n
   real(r8), dimension(ens_size) :: prs_sfc
   
   real(r8), allocatable, dimension(:)   :: thick, prs_tropomi, prs_tropomi_mem
   real(r8), allocatable, dimension(:,:) :: o3_val, tmp_val, qmr_val
   real(r8), allocatable, dimension(:)   :: o3_prf_mdl,tmp_prf_mdl,qmr_prf_mdl
   real(r8), allocatable, dimension(:)   :: prs_tropomi_top   
   logical  :: return_now,o3_return_now,tmp_return_now,qmr_return_now
   
   if ( .not. module_initialized ) call initialize_module
   
   pi       = 4.*atan(1.)
   rad2deg  = 360./(2.*pi)
   eps      =  0.61_r8
   Rd       = 287.05_r8     ! J/kg
   Ru       = 8.316_r8      ! J/kg
   Cp       = 1006.0        ! J/kg/K
   grav     =   9.8_r8
   o3_min   = 1.e-6_r8
   msq2cmsq = 1.e4_r8
   AvogN    = 6.02214e23_r8
   missing  =-888888_r8
   tmp_max  = 600.
   del_prs  = 5000.
   VMR_conv = 28.9644/47.9982
   molec2du = 1. / 2.6867e20
! 
! WACCM - MMR
! WRFChem - VMR ppmv
! TROPOMI O3 - DU   
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
   if(use_log_o3) then
      o3_min = log(o3_min)
   endif
   
! Assign vertical grid information

   layer_tropomi = nlayer(key)
   level_tropomi = nlayer(key)+1
!   kend_tropomi  = kend(key)
   layer_mdl=nlayer_model
   level_mdl=nlayer_model+1

   allocate(prs_tropomi(level_tropomi))
   allocate(prs_tropomi_mem(level_tropomi))
   prs_tropomi(1:level_tropomi)=pressure(key,1:level_tropomi)

! Get location infomation

   mloc = get_location(location)
   
   if (    mloc(2) >  90.0_r8) then
           mloc(2) =  90.0_r8
   elseif (mloc(2) < -90.0_r8) then
           mloc(2) = -90.0_r8
   endif
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

   o3_mdl_1(:)=missing_r8
   tmp_mdl_1(:)=missing_r8
   qmr_mdl_1(:)=missing_r8
   prs_mdl_1(:)=missing_r8

   do k=1,layer_mdl
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
      call interpolate(state_handle, ens_size, loc2, QTY_O3, o3_mdl_1, zstatus) ! ppmv 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_1, zstatus) ! K 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_1, zstatus) ! kg / kg 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_1, zstatus) ! Pa
!
      interp_new=0
      do imem=1,ens_size
         if(o3_mdl_1(imem).eq.missing_r8 .or. tmp_mdl_1(imem).eq.missing_r8 .or. &
         qmr_mdl_1(imem).eq.missing_r8 .or. prs_mdl_1(imem).eq.missing_r8) then
            interp_new=1
            exit
         endif
      enddo
      if(interp_new.eq.0) then
         exit
      endif    
   enddo
      
!   write(string1, *) 'APM: o3 lower bound ',o3_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: tmp lower bound ',tmp_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: qmr lower bound ',qmr_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: prs lower bound ',prs_mdl_1
!   call error_handler(E_MSG, routine, string1, source)

   o3_mdl_n(:)=missing_r8
   tmp_mdl_n(:)=missing_r8
   qmr_mdl_n(:)=missing_r8
   prs_mdl_n(:)=missing_r8

   do k=layer_mdl,1,-1
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
      call interpolate(state_handle, ens_size, loc2, QTY_O3, o3_mdl_n, zstatus) 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_n, &
      zstatus) 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_n, &
      zstatus) 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_n, &
      zstatus) 
!
      interp_new=0
      do imem=1,ens_size
         if(o3_mdl_n(imem).eq.missing_r8 .or. tmp_mdl_n(imem).eq.missing_r8 .or. &
         qmr_mdl_n(imem).eq.missing_r8 .or. prs_mdl_n(imem).eq.missing_r8) then
            interp_new=1
            exit
         endif
      enddo
      if(interp_new.eq.0) then
         exit
      endif    
   enddo

!   write(string1, *) 'APM: o3 upper bound ',o3_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: tmp upper bound ',tmp_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: qmr upper bound ',qmr_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: prs upper bound ',prs_mdl_n
!   call error_handler(E_MSG, routine, string1, source)

! Get cpsrs at TROPOMI pressure levels

   allocate( o3_val(ens_size,level_tropomi))
   allocate(tmp_val(ens_size,level_tropomi))
   allocate(qmr_val(ens_size,level_tropomi))

   do k=1,level_tropomi
      zstatus=0
      loc2 = set_location(mloc(1), mloc(2), prs_tropomi(k), VERTISPRESSURE)
      call interpolate(state_handle, ens_size, loc2, QTY_O3, o3_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_val(:,k), zstatus)  
!
! Correcting for expected failures near the surface
      do imem=1,ens_size
         if (prs_tropomi(k).ge.prs_mdl_1(imem)) then
            o3_val(imem,k) = o3_mdl_1(imem)
            tmp_val(imem,k) = tmp_mdl_1(imem)
            qmr_val(imem,k) = qmr_mdl_1(imem)
            cycle
         endif
!
! Correcting for expected failures near the top
         if (prs_tropomi(k).le.prs_mdl_n(imem)) then
            o3_val(imem,k) = o3_mdl_n(imem)
            tmp_val(imem,k) = tmp_mdl_n(imem)
            qmr_val(imem,k) = qmr_mdl_n(imem)
            cycle
         endif
      enddo

!      write(string1, *)'APM: o3 ',k,o3_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: tmp ',k,tmp_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: qmr ',k,qmr_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)

! Check data for missing values      
      do imem=1,ens_size
         if(o3_val(imem,k).eq.missing_r8 .or. tmp_val(imem,k).eq.missing_r8 .or. &
         qmr_val(imem,k).eq.missing_r8) then
            zstatus(:)=20
            expct_val(:)=missing_r8
            write(string1, *) 'APM: Input data has missing values ',o3_val(imem,k), &
            tmp_val(imem,k),qmr_val(imem,k)
            call error_handler(E_MSG, routine, string1, source)
            call track_status(ens_size, zstatus, expct_val, istatus, return_now)
            return
         endif
      enddo
!
! Convert units for o3 from ppmv
      o3_val(:,k) = o3_val(:,k) * 1.e-6_r8
   enddo
!
! Use large scale ozone data above the regional model top
   kstart(:)=-1
   do imem=1,ens_size
      do k=1,level_tropomi
         if (prs_tropomi(k).ge.prs_mdl_n(imem)) then
            kstart(imem)=k-1
!            write(string1, *) 'APM: imem,k-1,prs,mdl_n ',imem,k-1,prs_tropomi(k-1), &
!            prs_mdl_n(imem),prs_tropomi(k)
!            call error_handler(E_MSG, routine, string1, source)
            exit
         endif
      enddo
      if(kstart(imem).lt.0.) then
         write(string1, *) 'APM: Member ',imem,' kstart less than zero'
         call error_handler(E_MSG, routine, string1, source)
      endif   
      ncnt=kstart(imem)
      allocate(prs_tropomi_top(ncnt))
      allocate(o3_prf_mdl(ncnt),tmp_prf_mdl(ncnt),qmr_prf_mdl(ncnt))
      do k=1,kstart(imem)
         prs_tropomi_top(k)=prs_tropomi(k)
      enddo
      prs_tropomi_top(:)=prs_tropomi_top(:)/100.
!
      lon_obs=mloc(1)/rad2deg
      lat_obs=mloc(2)/rad2deg
      call get_time(obs_time,datesec_obs,date_obs)
!
      data_file=trim(upper_data_file)
      model=trim(upper_data_model)
      call get_upper_bdy_fld(fld,model,data_file,ls_chem_dx,ls_chem_dy, &
      ls_chem_dz,ls_chem_dt,lon_obs,lat_obs,prs_tropomi_top, &
      ncnt,o3_prf_mdl,tmp_prf_mdl,qmr_prf_mdl,date_obs,datesec_obs)

      o3_prf_mdl(:)=o3_prf_mdl(:)*VMR_conv
!
! Save upper BC data in the cpsrs   
      do k=1,kstart(imem)
         o3_val(imem,k)=o3_prf_mdl(k)
         tmp_val(imem,k)=tmp_prf_mdl(k)
         qmr_val(imem,k)=qmr_prf_mdl(k)
      enddo
      deallocate(prs_tropomi_top)
      deallocate(o3_prf_mdl,tmp_prf_mdl,qmr_prf_mdl)
!      do k=1,level_tropomi
!         write(string1, *) 'APM: imem,k,prs,o3,tmp,qmr ',imem,k,prs_tropomi(k), &
!         o3_val(imem,k),tmp_val(imem,k),qmr_val(imem,k)
!         call error_handler(E_MSG, routine, string1, source)
!      enddo
   enddo
!
! Impose ensemble perturbations from level kstart+1   
   do imem=1,ens_size
      ensavg_o3=0.
      ensavg_tmp=0.
      ensavg_qmr=0.
      do k=1,ens_size
         ensavg_o3=ensavg_o3+o3_val(k,kstart(imem)+1)/real(ens_size)
         ensavg_tmp=ensavg_tmp+tmp_val(k,kstart(imem)+1)/real(ens_size)
         ensavg_qmr=ensavg_qmr+qmr_val(k,kstart(imem)+1)/real(ens_size)
!         write(string1, *) 'APM: imem,k,kstart,ensavg,o3_val ', &
!         imem,k,kstart(imem),ensavg_o3,o3_val(k,kstart(imem)-1)
!         call error_handler(E_MSG, routine, string1, source)
      enddo
!      write(string1, *) 'APM: o3, tmp, qmr ',imem,ensavg_o3,ensavg_tmp,ensavg_qmr 
!      call error_handler(E_MSG, routine, string1, source)
!
      fac_o3=o3_val(imem,kstart(imem)+1)/ensavg_o3
      fac_tmp=tmp_val(imem,kstart(imem)+1)/ensavg_tmp
      fac_qmr=qmr_val(imem,kstart(imem)+1)/ensavg_qmr      
      do k=1,kstart(imem)
         o3_val(imem,k)=o3_val(imem,k)*fac_o3
         tmp_val(imem,k)=tmp_val(imem,k)*fac_tmp
         qmr_val(imem,k)=qmr_val(imem,k)*fac_qmr
      enddo
   enddo
!   do k=1,level_tropomi
!      write(string1, *) 'APM: o3 ',k,o3_val(1,k),o3_val(int(ens_size/2),k), &
!      o3_val(ens_size,k)
!      call error_handler(E_MSG, routine, string1, source)
!   enddo
   istatus=0
   zstatus(:)=0.
   expct_val(:)=0.0
   allocate(thick(layer_tropomi))

   do imem=1,ens_size
! Adjust the TROPOMI pressure for WRF-Chem lower/upper boudary pressure
! (TROPOMI O3 vertical grid is top to bottom)
      prs_tropomi_mem(:)=prs_tropomi(:)
      if (prs_sfc(imem).gt.prs_tropomi_mem(level_tropomi)) then
         prs_tropomi_mem(level_tropomi)=prs_sfc(imem)
      endif   

! Calculate the thicknesses

      thick(:)=0.
      do k=1,layer_tropomi
         lnpr_mid=(log(prs_tropomi_mem(k))+log(prs_tropomi_mem(k+1)))/2.
         up_wt=log(prs_tropomi_mem(k+1))-lnpr_mid
         dw_wt=log(lnpr_mid)-log(prs_tropomi_mem(k))
         tl_wt=up_wt+dw_wt
         tmp_vir_k  = (1.0_r8 + eps*qmr_val(imem,k))*tmp_val(imem,k)
         tmp_vir_kp = (1.0_r8 + eps*qmr_val(imem,k+1))*tmp_val(imem,k+1)
         thick(k)   = Rd*(dw_wt*tmp_vir_kp + up_wt*tmp_vir_k)/tl_wt/grav* &
         log(prs_tropomi_mem(k+1)/prs_tropomi_mem(k))
      enddo

! Process the vertical summation
   
      do k=1,layer_tropomi
         lnpr_mid=(log(prs_tropomi_mem(k))+log(prs_tropomi_mem(k+1)))/2.
         up_wt=log(prs_tropomi_mem(k+1))-lnpr_mid
         dw_wt=log(lnpr_mid)-log(prs_tropomi_mem(k))
         tl_wt=up_wt+dw_wt
   
! Convert from VMR to molar density (mol/m^3)
         if(use_log_o3) then
            o3_val_conv = (dw_wt*exp(o3_val(imem,k+1))+up_wt*exp(o3_val(imem,k)))/tl_wt * &
            (dw_wt*prs_tropomi_mem(k+1)+up_wt*prs_tropomi_mem(k)) / &
            (Ru*(dw_wt*tmp_val(imem,k+1)+up_wt*tmp_val(imem,k)))
         else
            o3_val_conv = (dw_wt*o3_val(imem,k+1)+up_wt*o3_val(imem,k))/tl_wt * &
            (dw_wt*prs_tropomi_mem(k+1)+up_wt*prs_tropomi_mem(k)) / &
            (Ru*(dw_wt*tmp_val(imem,k+1)+up_wt*tmp_val(imem,k)))
         endif
!
! Convert from mol/m^2 to DU 
         o3_val_conv=o3_val_conv*AvogN*molec2du
 
! Get expected observation
         expct_val(imem) = expct_val(imem) + thick(k) * o3_val_conv * &
         avg_kernel(key,k)
         if(imem.eq.1) then         
            write(string1, *) 'APM Summation : ', &
            k,expct_val(imem),(prs_tropomi_mem(k)+prs_tropomi_mem(k+1))/2., &
            thick(k), o3_val_conv, avg_kernel(key,k)
            call error_handler(E_MSG, routine, string1, source)
         endif
         expct_val(imem)=expct_val(imem)
      enddo
!
      if(imem.eq.1) then         
         write(string1, *) &
         'APM: Member ',imem,'Expected Value ',expct_val(imem)
         call error_handler(E_MSG, routine, string1, source)
      endif
!      
      if(isnan(expct_val(imem))) then
         zstatus(imem)=20
         expct_val(:)=missing_r8
         write(string1, *) &
         'APM NOTICE: TROPOMI O3 expected value is NaN '
         call error_handler(E_MSG, routine, string1, source)
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
!
      if(expct_val(imem).lt.0) then
         zstatus(imem)=20
         expct_val(:)=missing_r8
         write(string1, *) &
         'APM NOTICE: TROPOMI O3 expected value is negative '
         call error_handler(E_MSG, routine, string1, source)
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
   enddo

! Clean up and return
   deallocate(o3_val, tmp_val, qmr_val)
   deallocate(thick)
   deallocate(prs_tropomi, prs_tropomi_mem)

end subroutine get_expected_tropomi_o3_cpsr

!-------------------------------------------------------------------------------

subroutine set_obs_def_tropomi_o3_cpsr(key, o3_pressure, o3_avg_kernel, o3_prior, o3_nlayer)

   integer,                           intent(in)   :: key, o3_nlayer
   real(r8),                          intent(in)   :: o3_prior
   real(r8), dimension(o3_nlayer+1),  intent(in)   :: o3_pressure
   real(r8), dimension(o3_nlayer),    intent(in)   :: o3_avg_kernel
   
   if ( .not. module_initialized ) call initialize_module
   
   if(num_tropomi_o3_obs >= max_tropomi_o3_obs) then
      write(string1, *)'Not enough space for tropomi o3 obs.'
      write(string2, *)'Can only have max_tropomi_o3_obs (currently ',max_tropomi_o3_obs,')'
      call error_handler(E_ERR,'set_obs_def_tropomi_o3_cpsr',string1,source,revision, &
      revdate,text2=string2)
   endif
   
   nlayer(key) = o3_nlayer
   prior(key)  = o3_prior
   pressure(key,1:o3_nlayer+1) = o3_pressure(1:o3_nlayer+1)
   avg_kernel(key,1:o3_nlayer) = o3_avg_kernel(1:o3_nlayer)
   
end subroutine set_obs_def_tropomi_o3_cpsr

end module obs_def_tropomi_o3_cpsr_mod

! END DART PREPROCESS MODULE CODE
