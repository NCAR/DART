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
! OMI_SO2_PBL_COL, QTY_SO2
! END DART PREPROCESS TYPE DEFINITIONS
!
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_omi_so2_pbl_col_mod, only : get_expected_omi_so2_pbl_col, &
!                                   read_omi_so2_pbl_col, &
!                                   write_omi_so2_pbl_col, &
!                                   interactive_omi_so2_pbl_col
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(OMI_SO2_PBL_COL)                                                           
!         call get_expected_omi_so2_pbl_col(state_handle, ens_size, location, obs_def%key, obs_time, expected_obs, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!
! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(OMI_SO2_PBL_COL)
!         call read_omi_so2_pbl_col(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(OMI_SO2_PBL_COL)
!         call write_omi_so2_pbl_col(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF
!
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(OMI_SO2_PBL_COL)
!         call interactive_omi_so2_pbl_col(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_omi_so2_pbl_col_mod

   use         apm_upper_bdy_mod, only :get_upper_bdy_fld, &
                                        get_MOZART_INT_DATA, &
                                        get_MOZART_REAL_DATA, &
                                        wrf_dart_ubval_interp, &
                                        apm_get_exo_coldens, &
                                        apm_get_upvals, &
                                        apm_interpolate

   use         types_mod, only : r8, MISSING_R8

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

   use          obs_kind_mod, only : QTY_SO2, QTY_TEMPERATURE, QTY_SURFACE_PRESSURE, &
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

   public :: write_omi_so2_pbl_col, &
          read_omi_so2_pbl_col, &
          interactive_omi_so2_pbl_col, &
          get_expected_omi_so2_pbl_col, &
          set_obs_def_omi_so2_pbl_col

! Storage for the special information required for observations of this type
   integer, parameter    :: max_omi_so2_obs = 10000000
   integer               :: num_omi_so2_obs = 0
   integer,  allocatable :: nlayer(:)
   integer,  allocatable :: kend(:)
   real(r8), allocatable :: pressure(:,:)
   real(r8), allocatable :: scat_wt(:,:)

! version controlled file description for error handling, do not edit
   character(len=*), parameter :: source   = 'obs_def_omi_so2_pbl_col_mod.f90'
   character(len=*), parameter :: revision = ''
   character(len=*), parameter :: revdate  = ''
   
   character(len=512) :: string1, string2
   character(len=200) :: upper_data_file
   character(len=200) :: upper_data_model
   character(len=200) :: model
   integer            :: ls_chem_dx, ls_chem_dy, ls_chem_dz, ls_chem_dt
   
   logical, save :: module_initialized = .false.

! Namelist with default values
   logical :: use_log_so2   = .false.
   integer :: nlayer_model = -9999
   integer :: nlayer_omi = -9999
   integer :: nlayer_omi_so2_total_col = -9999
   integer :: nlayer_omi_so2_pbl_col = -9999
   
   namelist /obs_def_OMI_SO2_nml/ upper_data_file, use_log_so2, &
   nlayer_model, nlayer_omi_so2_total_col, nlayer_omi_so2_pbl_col, &
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
   call find_namelist_in_file("input.nml", "obs_def_OMI_SO2_nml", iunit)
   read(iunit, nml = obs_def_OMI_SO2_nml, iostat = rc)
   call check_namelist_read(iunit, rc, "obs_def_OMI_SO2_nml")

! Record the namelist values
   if (do_nml_file()) write(nmlfileunit, nml=obs_def_OMI_SO2_nml)
   if (do_nml_term()) write(     *     , nml=obs_def_OMI_SO2_nml)
   nlayer_omi=nlayer_omi_so2_pbl_col

! Check for valid values

   if (nlayer_model < 1) then
      write(string1,*)'obs_def_OMI_SO2_nml:nlayer_model must be > 0, it is ',nlayer_model
      call error_handler(E_ERR,'initialize_module',string1,source)
   endif
   
   if (nlayer_omi < 1) then
      write(string1,*)'obs_def_OMI_SO2_nml:nlayer_omi must be > 0, it is ',nlayer_omi
      call error_handler(E_ERR,'initialize_module',string1,source)
   endif
   
   allocate(   nlayer(max_omi_so2_obs))
   allocate(   kend(max_omi_so2_obs))
   allocate( pressure(max_omi_so2_obs,nlayer_omi+1))
   allocate(  scat_wt(max_omi_so2_obs,nlayer_omi))

end subroutine initialize_module

!-------------------------------------------------------------------------------

subroutine read_omi_so2_pbl_col(key, ifile, fform)

   integer,          intent(out)          :: key
   integer,          intent(in)           :: ifile
   character(len=*), intent(in), optional :: fform
   
! temporary arrays to hold buffer till we decide if we have enough room

   integer               :: keyin
   integer               :: nlayer_1
   integer               :: kend_1
   real(r8), allocatable :: pressure_1(:)
   real(r8), allocatable :: scat_wt_1(:)
   character(len=32)     :: fileformat
   
   integer, SAVE :: counts1 = 0
   
   if ( .not. module_initialized ) call initialize_module
   
   fileformat = "ascii" 
   if(present(fform)) fileformat = adjustl(fform)

! Need to know how many layers for this one
   nlayer_1   = read_int_scalar( ifile, fileformat, 'nlayer_1')
   kend_1   = read_int_scalar( ifile, fileformat, 'kend_1')
   
   allocate( pressure_1(nlayer_1+1))
   allocate(  scat_wt_1(nlayer_1))
   
   call read_r8_array(ifile, nlayer_1+1, pressure_1,   fileformat, 'pressure_1')
   call read_r8_array(ifile, nlayer_1,   scat_wt_1, fileformat, 'scat_wt_1')
   keyin = read_int_scalar(ifile, fileformat, 'keyin')
   
   counts1 = counts1 + 1
   key     = counts1
   
   if(counts1 > max_omi_so2_obs) then
      write(string1, *)'Not enough space for omi so2 obs.'
      write(string2, *)'Can only have max_omi_so2_obs (currently ',max_omi_so2_obs,')'
      call error_handler(E_ERR,'read_omi_so2_pbl_col',string1,source,text2=string2)
   endif
   
   call set_obs_def_omi_so2_pbl_col(key, pressure_1, scat_wt_1, nlayer_1, kend_1)
   
   deallocate(pressure_1, scat_wt_1)

end subroutine read_omi_so2_pbl_col

!-------------------------------------------------------------------------------

subroutine write_omi_so2_pbl_col(key, ifile, fform)

   integer,          intent(in)           :: key
   integer,          intent(in)           :: ifile
   character(len=*), intent(in), optional :: fform
   
   character(len=32) :: fileformat
   
   if ( .not. module_initialized ) call initialize_module
   
   fileformat = "ascii"
   if(present(fform)) fileformat = adjustl(fform)

! nlayer, pressure, and scat_wt are all scoped in this module
! you can come extend the context strings to include the key if desired.

   call write_int_scalar(ifile,                     nlayer(key), fileformat,'nlayer')
   call write_int_scalar(ifile,                     kend(key), fileformat,'kend')
   call write_r8_array(  ifile, nlayer(key)+1,  pressure(key,:), fileformat,'pressure')
   call write_r8_array(  ifile, nlayer(key),  scat_wt(key,:), fileformat,'scat_wt')
   call write_int_scalar(ifile,                             key, fileformat,'key')

end subroutine write_omi_so2_pbl_col

!-------------------------------------------------------------------------------

subroutine interactive_omi_so2_pbl_col(key)

   integer, intent(out) :: key
   
   if ( .not. module_initialized ) call initialize_module

! STOP because routine is not finished.
   write(string1,*)'interactive_omi_so2_pbl_col not yet working.'
   call error_handler(E_ERR, 'interactive_omi_so2_pbl_col', string1, source)
   
   if(num_omi_so2_obs >= max_omi_so2_obs) then
      write(string1, *)'Not enough space for an omi so2 obs.'
      write(string2, *)'Can only have max_omi_so2_obs (currently ',max_omi_so2_obs,')'
      call error_handler(E_ERR, 'interactive_omi_so2_pbl_col', string1, &
                 source, text2=string2)
   endif

! Increment the index
   num_omi_so2_obs = num_omi_so2_obs + 1
   key            = num_omi_so2_obs

! Otherwise, prompt for input for the three required beasts

   write(*, *) 'Creating an interactive_omi_so2_pbl_col observation'
   write(*, *) 'This featue is not setup '

end subroutine interactive_omi_so2_pbl_col

!-------------------------------------------------------------------------------

subroutine get_expected_omi_so2_pbl_col(state_handle, ens_size, location, key, obs_time, expct_val, istatus)

   type(ensemble_type), intent(in)  :: state_handle
   type(location_type), intent(in)  :: location
   integer,             intent(in)  :: ens_size
   integer,             intent(in)  :: key
   type(time_type),     intent(in)  :: obs_time
   integer,             intent(out) :: istatus(:)
   real(r8),            intent(out) :: expct_val(:)

   character(len=*), parameter :: routine = 'get_expected_omi_so2_pbl_col'
   character(len=120)          :: data_file
   character(len=*),parameter  :: fld = 'SO2_VMR_inst'
   type(location_type) :: loc2

   integer :: layer_omi,level_omi, klev_omi, kend_omi
   integer :: layer_mdl,level_mdl
   integer :: k,kk,imem,imemm,flg
   integer :: interp_new
   integer :: icnt,ncnt,kstart
   integer :: date_obs,datesec_obs,pbl_index
   integer, dimension(ens_size) :: zstatus,kbnd_1,kbnd_n
   
   real(r8) :: eps, AvogN, Rd, Ru, Cp, grav, msq2cmsq
   real(r8) :: missing,so2_min,tmp_max
   real(r8) :: level,del_prs,prior_term,pbl_sum
   real(r8) :: tmp_vir_k,tmp_vir_kp
   real(r8) :: mloc(3),obs_prs
   real(r8) :: so2_val_conv, VMR_conv
   real(r8) :: up_wt,dw_wt,tl_wt,lnpr_mid
   real(r8) :: lon_obs,lat_obs,pi,rad2deg

   real(r8), dimension(ens_size) :: so2_mdl_tmp, tmp_mdl_tmp, qmr_mdl_tmp, prs_mdl_tmp
   real(r8), dimension(ens_size) :: so2_mdl_1, tmp_mdl_1, qmr_mdl_1, prs_mdl_1
   real(r8), dimension(ens_size) :: so2_mdl_n, tmp_mdl_n, qmr_mdl_n, prs_mdl_n
   real(r8), dimension(ens_size) :: prs_sfc,rec_so2_val,rec_tmp_val,rec_qmr_val

   real(r8), allocatable, dimension(:)   :: thick, prs_omi, prs_omi_mem
   real(r8), allocatable, dimension(:,:) :: so2_val, tmp_val, qmr_val
   real(r8), allocatable, dimension(:)   :: so2_prf_mdl,tmp_prf_mdl,qmr_prf_mdl
   real(r8), allocatable, dimension(:)   :: prs_omi_top
   logical  :: return_now,so2_return_now,tmp_return_now,qmr_return_now
   
   if ( .not. module_initialized ) call initialize_module

   pi       = 4.*atan(1.)
   rad2deg  = 360./(2.*pi)
   eps      =  0.61_r8
   Rd       = 287.05_r8     ! J/kg
   Ru       = 8.316_r8      ! J/kg
   Cp       = 1006.0        ! J/kg/K
   grav     = 9.8_r8
   so2_min  = 1.e-6_r8
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
   if(use_log_so2) then
      so2_min = log(so2_min)
   endif

! Assign vertical grid information (OMI SO2 grid is top to bottom)

   layer_omi = nlayer(key)
   level_omi = nlayer(key)+1
!  klev_omi  = klev(key)
   kend_omi  = kend(key)
   layer_mdl = nlayer_model
   level_mdl = nlayer_model+1
   
   allocate(prs_omi(level_omi))
   allocate(prs_omi_mem(level_omi))
   prs_omi(1:level_omi)=pressure(key,1:level_omi)*100.

! Get location infomation

   mloc = get_location(location)

   if (    mloc(2) >  90.0_r8) then
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

   so2_mdl_1(:)=missing_r8
   tmp_mdl_1(:)=missing_r8
   qmr_mdl_1(:)=missing_r8
   prs_mdl_1(:)=missing_r8

   kbnd_1(:)=1
   do k=1,layer_mdl
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
      call interpolate(state_handle, ens_size, loc2, QTY_SO2, so2_mdl_1, zstatus) ! ppmv 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_1, zstatus) ! K 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_1, zstatus) ! kg / kg 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_1, zstatus) ! Pa

      interp_new=0
      do imem=1,ens_size
         if(so2_mdl_1(imem).eq.missing_r8 .or. tmp_mdl_1(imem).eq.missing_r8 .or. &
         qmr_mdl_1(imem).eq.missing_r8 .or. prs_mdl_1(imem).eq.missing_r8) then
            interp_new=1
            exit
         endif
      enddo
      if(interp_new.eq.0) then
         exit
      endif
   enddo

!   write(string1, *) 'APM: so2 lower bound ',key,so2_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: tmp lower bound ',key,tmp_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: qmr lower bound ',key,qmr_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: prs lower bound ',key,prs_mdl_1
!   call error_handler(E_MSG, routine, string1, source)

   so2_mdl_n(:)=missing_r8
   tmp_mdl_n(:)=missing_r8
   qmr_mdl_n(:)=missing_r8
   prs_mdl_n(:)=missing_r8

   kbnd_n(:)=layer_mdl
   do k=layer_mdl,1,-1
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
      call interpolate(state_handle, ens_size, loc2, QTY_SO2, so2_mdl_n, zstatus) 
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
         if(so2_mdl_n(imem).eq.missing_r8 .or. tmp_mdl_n(imem).eq.missing_r8 .or. &
         qmr_mdl_n(imem).eq.missing_r8 .or. prs_mdl_n(imem).eq.missing_r8) then
            interp_new=1
            exit
         endif
      enddo
      if(interp_new.eq.0) then
         exit
      endif
   enddo

!   write(string1, *) 'APM: so2 upper bound ',key,so2_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: tmp upper bound ',key,tmp_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: qmr upper bound ',key,qmr_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *) 'APM: prs upper bound ',key,prs_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
   
! Get values at OMI pressure levels (Pa)

   allocate(so2_val(ens_size,level_omi))
   allocate(tmp_val(ens_size,level_omi))
   allocate(qmr_val(ens_size,level_omi))

   do k=1,level_omi
      zstatus=0
      loc2 = set_location(mloc(1), mloc(2), prs_omi(k), VERTISPRESSURE)
      call interpolate(state_handle, ens_size, loc2, QTY_SO2, so2_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_val(:,k), zstatus)  
!
! Correcting for expected failures near the surface
      do imem=1,ens_size
         if (prs_omi(k).ge.prs_mdl_1(imem)) then
            so2_val(imem,k) = so2_mdl_1(imem)
            tmp_val(imem,k) = tmp_mdl_1(imem)
            qmr_val(imem,k) = qmr_mdl_1(imem)
         endif
!
! Correcting for expected failures near the top
         if (prs_omi(k).le.prs_mdl_n(imem)) then
            so2_val(imem,k) = so2_mdl_n(imem)
            tmp_val(imem,k) = tmp_mdl_n(imem)
            qmr_val(imem,k) = qmr_mdl_n(imem)
         endif
      enddo
!
!      write(string1, *)'APM: so2 ',key,k,so2_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: tmp ',key,k,tmp_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!      write(string1, *)'APM: qmr ',key,k,qmr_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!
! Check data for missing values      
      do imem=1,ens_size
         if(so2_val(imem,k).eq.missing_r8 .or. tmp_val(imem,k).eq.missing_r8 .or. &
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
! Convert units for so2 from ppmv
      so2_val(:,k) = so2_val(:,k) * 1.e-6_r8
   enddo
   so2_mdl_1(:)=so2_mdl_1(:) * 1.e-6_r8
   so2_mdl_n(:)=so2_mdl_n(:) * 1.e-6_r8
!   imem=1
!   do k=1,level_omi
!      write(string1, *)'APM: ',k,prs_omi(k), &
!      so2_val(imem,k)*1.e6_8,tmp_val(imem,k),qmr_val(imem,k)
!      call error_handler(E_MSG, routine, string1, source)
!   enddo                                                                                                    
!
! Use large scale so2 data above the regional model top
! OMI vertical is from top to bottom   
   kstart=-1
   do imem=1,ens_size
      if (prs_omi(1).lt.prs_mdl_n(imem)) then
         do k=1,level_omi
            if (prs_omi(k).ge.prs_mdl_n(imem)) then
               kstart=k
               exit
            endif
         enddo
         ncnt=kstart
         allocate(prs_omi_top(ncnt))
         allocate(so2_prf_mdl(ncnt),tmp_prf_mdl(ncnt),qmr_prf_mdl(ncnt))
         do k=1,kstart
            prs_omi_top(k)=prs_omi(k)
         enddo
         prs_omi_top(:)=prs_omi_top(:)/100.
!
         lon_obs=mloc(1)/rad2deg
         lat_obs=mloc(2)/rad2deg
         call get_time(obs_time,datesec_obs,date_obs)
!
         data_file=trim(upper_data_file)
         model=trim(upper_data_model)
         call get_upper_bdy_fld(fld,model,data_file,ls_chem_dx,ls_chem_dy, &
         ls_chem_dz,ls_chem_dt,lon_obs,lat_obs,prs_omi_top, &
         ncnt,so2_prf_mdl,tmp_prf_mdl,qmr_prf_mdl,date_obs,datesec_obs)
!
! Impose ensemble perturbations from level kstart+1      
         do k=1,kstart 
            so2_val(imem,k)=so2_prf_mdl(k)*so2_val(imem,kstart+1)/ &
            (sum(so2_val(:,kstart+1))/real(ens_size))
            tmp_val(imem,k)=tmp_prf_mdl(k)*tmp_val(imem,kstart+1)/ &
            (sum(tmp_val(:,kstart+1))/real(ens_size))
            qmr_val(imem,k)=qmr_prf_mdl(k)*qmr_val(imem,kstart+1)/ &
            (sum(qmr_val(:,kstart+1))/real(ens_size))
         enddo
         deallocate(prs_omi_top)
         deallocate(so2_prf_mdl,tmp_prf_mdl,qmr_prf_mdl)
      endif
   enddo
!
! Print full profile examples
!   do imem=1,1
!      do k=1,level_omi
!         write(string1, *) &
!         'APM: prs,so2,tmp,qmr ',k,prs_omi(k)/100.,so2_val(imem,k), &
!         tmp_val(imem,k),qmr_val(imem,k)
!         call error_handler(E_MSG, routine, string1, source)
!      enddo
!   enddo
!
! Check full profile for negative values
   do imem=1,ens_size
      flg=0
      do k=1,level_omi   
         if(so2_val(imem,k).lt.0. .or. tmp_val(imem,k).lt.0. .or. &
         qmr_val(imem,k).lt.0.) then
            flg=1   
            write(string1, *) &
            'APM: Recentered full profile has negative values for key,imem ',key,imem
            call error_handler(E_MSG, routine, string1, source)
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
! Calculate the expected retrieval   
   istatus(:)=0
   zstatus(:)=0.
   expct_val(:)=0.0
   allocate(thick(layer_omi))

   prs_omi_mem(:)=prs_omi(:)
   do imem=1,ens_size

! Calculate the thicknesses (grid is top to bottom)

      thick(:)=0.
      do k=1,layer_omi
         lnpr_mid=(log(prs_omi_mem(k+1))+log(prs_omi_mem(k)))/2.
         up_wt=log(prs_omi_mem(k+1))-lnpr_mid
         dw_wt=lnpr_mid-log(prs_omi_mem(k))
         tl_wt=up_wt+dw_wt 
         tmp_vir_k  = (1.0_r8 + eps*qmr_val(imem,k))*tmp_val(imem,k)
         tmp_vir_kp = (1.0_r8 + eps*qmr_val(imem,k+1))*tmp_val(imem,k+1)
         thick(k)   = Rd*(dw_wt*tmp_vir_kp + up_wt*tmp_vir_k)/tl_wt/grav* &
         log(prs_omi_mem(k+1)/prs_omi_mem(k))

!         write(string1, *) &
!         'APM: Key, Thickness calcs ', k, thick(k), up_wt, dw_wt, tmp_vir_k,tmp_vir_kp, &
!         prs_omi_mem(k), prs_omi_mem(k+1)     
!         call error_handler(E_MSG, routine, string1, source)
         
      enddo
!
! Find PBL index      
      pbl_sum=0.
      pbl_index=1
      do k=1,layer_omi
         kk=layer_omi-k+1
         pbl_sum=pbl_sum+thick(kk)
         if(pbl_sum.ge.1000.) then
            pbl_index=kk
            exit
         endif
      enddo
      if(pbl_index.lt.layer_omi/2 .or. pbl_index.ge.layer_omi) then
         zstatus(imem)=20
         expct_val(:)=missing_r8
         write(string1, *) &
         'APM NOTICE: OMI SO2 PBL COL - pbl_index is unreasonable ',pbl_index
         call error_handler(E_ALLMSG, routine, string1, source)
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
      
! Process the vertical summation (OMI SO2 units are mole per m^2)

      do k=pbl_index,layer_omi
         lnpr_mid=(log(prs_omi_mem(k+1))+log(prs_omi_mem(k)))/2.
         up_wt=log(prs_omi_mem(k+1))-lnpr_mid
         dw_wt=lnpr_mid-log(prs_omi_mem(k))
         tl_wt=up_wt+dw_wt

! Convert from VMR to molar density (mol/m^3)
         if(use_log_so2) then
            so2_val_conv = (up_wt*exp(so2_val(imem,k))+dw_wt*exp(so2_val(imem,k+1)))/tl_wt * &
            (up_wt*prs_omi_mem(k)+dw_wt*prs_omi_mem(k+1)) / &
            (Ru*(up_wt*tmp_val(imem,k)+dw_wt*tmp_val(imem,k+1)))
         else
            so2_val_conv = (up_wt*so2_val(imem,k)+dw_wt*so2_val(imem,k+1))/tl_wt * &
            (up_wt*prs_omi_mem(k)+dw_wt*prs_omi_mem(k+1)) / &
            (Ru*(up_wt*tmp_val(imem,k)+dw_wt*tmp_val(imem,k+1)))
         endif
 
! Get expected observation (molec/cm^2)

         expct_val(imem) = expct_val(imem) + thick(k) * so2_val_conv * &
         AvogN/msq2cmsq * scat_wt(key,k)

!         write(string1, *) &
!         'APM: Key, Expected Value Terms ',k, expct_val(imem), thick(k), &
!         thick(k)*so2_val_conv, so2_val_conv, scat_wt(key,k)
!         call error_handler(E_MSG, routine, string1, source)

      enddo
!      write(string1, *) &
!      'APM: Expected Value is ',imem, expct_val(imem)
!      call error_handler(E_MSG, routine, string1, source)

      if(isnan(expct_val(imem))) then
         zstatus(imem)=20
         expct_val(:)=missing_r8
         write(string1, *) &
         'APM NOTICE: OMI SO2 expected value is NaN '
         call error_handler(E_MSG, routine, string1, source)
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
!
      if(expct_val(imem).le.0) then
         zstatus(imem)=20
         expct_val(:)=missing_r8
         write(string1, *) &
         'APM NOTICE: OMI SO2 expected value is negative '
         call error_handler(E_MSG, routine, string1, source)
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
   enddo

! Clean up and return
   deallocate(so2_val, tmp_val, qmr_val)
   deallocate(thick)
   deallocate(prs_omi, prs_omi_mem)

end subroutine get_expected_omi_so2_pbl_col

!-------------------------------------------------------------------------------

subroutine set_obs_def_omi_so2_pbl_col(key, so2_pressure, so2_scat_wt, so2_nlayer, so2_kend)

   integer,                            intent(in)   :: key, so2_nlayer, so2_kend
   real(r8), dimension(so2_nlayer+1),  intent(in)   :: so2_pressure
   real(r8), dimension(so2_nlayer),    intent(in)   :: so2_scat_wt
   
   if ( .not. module_initialized ) call initialize_module
   
   if(num_omi_so2_obs >= max_omi_so2_obs) then
      write(string1, *)'Not enough space for omi so2 obs.'
      write(string2, *)'Can only have max_omi_so2_obs (currently ',max_omi_so2_obs,')'
      call error_handler(E_ERR,'set_obs_def_omi_so2_pbl_col',string1,source,revision, &
      revdate,text2=string2)
   endif
   
   nlayer(key) = so2_nlayer
   kend(key) = so2_kend
   pressure(key,1:so2_nlayer+1) = so2_pressure(1:so2_nlayer+1)
   scat_wt(key,1:so2_nlayer) = so2_scat_wt(1:so2_nlayer)

end subroutine set_obs_def_omi_so2_pbl_col

end module obs_def_omi_so2_pbl_col_mod

! END DART PREPROCESS MODULE CODE
