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
! OMI_HCHO_TOTAL_COL, QTY_HCHO
! END DART PREPROCESS TYPE DEFINITIONS
!
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_omi_hcho_total_col_mod, only : get_expected_omi_hcho_total_col, &
!                                   read_omi_hcho_total_col, &
!                                   write_omi_hcho_total_col, &
!                                   interactive_omi_hcho_total_col
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(OMI_HCHO_TOTAL_COL)                                                           
!         call get_expected_omi_hcho_total_col(state_handle, ens_size, location, obs_def%key, obs_time, expected_obs, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!
! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(OMI_HCHO_TOTAL_COL)
!         call read_omi_hcho_total_col(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(OMI_HCHO_TOTAL_COL)
!         call write_omi_hcho_total_col(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF
!
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(OMI_HCHO_TOTAL_COL)
!         call interactive_omi_hcho_total_col(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_omi_hcho_total_col_mod

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

use          obs_kind_mod, only : QTY_HCHO, QTY_TEMPERATURE, QTY_SURFACE_PRESSURE, &
                                  QTY_PRESSURE, QTY_VAPOR_MIXING_RATIO

use  ensemble_manager_mod, only : ensemble_type

use obs_def_utilities_mod, only : track_status

use mpi_utilities_mod,     only : my_task_id

use      time_manager_mod, only : time_type, get_date, set_date, get_time, set_time

implicit none
private

public :: write_omi_hcho_total_col, &
          read_omi_hcho_total_col, &
          interactive_omi_hcho_total_col, &
          get_expected_omi_hcho_total_col, &
          set_obs_def_omi_hcho_total_col

! Storage for the special information required for observations of this type
integer, parameter    :: max_omi_hcho_obs = 10000000
integer               :: num_omi_hcho_obs = 0
integer,  allocatable :: nlayer(:)
integer,  allocatable :: kend(:)
real(r8), allocatable :: pressure(:,:)
real(r8), allocatable :: scat_wt(:,:)
! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_omi_hcho_total_col_mod.f90'
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
integer :: nlayer_omi = -9999
integer :: nlayer_omi_hcho_total_col = -9999
integer :: nlayer_omi_hcho_trop_col = -9999

namelist /obs_def_OMI_HCHO_nml/ upper_data_file, use_log_hcho, &
nlayer_model, nlayer_omi_hcho_total_col, nlayer_omi_hcho_trop_col, &
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
call find_namelist_in_file("input.nml", "obs_def_OMI_HCHO_nml", iunit)
read(iunit, nml = obs_def_OMI_HCHO_nml, iostat = rc)
call check_namelist_read(iunit, rc, "obs_def_OMI_HCHO_nml")

! Record the namelist values
if (do_nml_file()) write(nmlfileunit, nml=obs_def_OMI_HCHO_nml)
if (do_nml_term()) write(     *     , nml=obs_def_OMI_HCHO_nml)
nlayer_omi=nlayer_omi_hcho_total_col

! Check for valid values

if (nlayer_model < 1) then
   write(string1,*)'obs_def_OMI_HCHO_nml:nlayer_model must be > 0, it is ',nlayer_model
   call error_handler(E_ERR,'initialize_module',string1,source)
endif

if (nlayer_omi < 1) then
   write(string1,*)'obs_def_OMI_HCHO_nml:nlayer_omi must be > 0, it is ',nlayer_omi
   call error_handler(E_ERR,'initialize_module',string1,source)
endif

allocate(   nlayer(max_omi_hcho_obs))
allocate(   kend(max_omi_hcho_obs))
allocate( pressure(max_omi_hcho_obs,nlayer_omi))
allocate(  scat_wt(max_omi_hcho_obs,nlayer_omi))

end subroutine initialize_module

!-------------------------------------------------------------------------------

subroutine read_omi_hcho_total_col(key, ifile, fform)

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
kend_1     = read_int_scalar( ifile, fileformat, 'kend_1')

allocate( pressure_1(nlayer_1))
allocate(  scat_wt_1(nlayer_1))

call read_r8_array(ifile, nlayer_1, pressure_1,   fileformat, 'pressure_1')
call read_r8_array(ifile, nlayer_1,   scat_wt_1, fileformat, 'scat_wt_1')
keyin = read_int_scalar(ifile, fileformat, 'keyin')

counts1 = counts1 + 1
key     = counts1

if(counts1 > max_omi_hcho_obs) then
   write(string1, *)'Not enough space for omi hcho obs.'
   write(string2, *)'Can only have max_omi_hcho_obs (currently ',max_omi_hcho_obs,')'
   call error_handler(E_ERR,'read_omi_hcho_total_col',string1,source,text2=string2)
endif

call set_obs_def_omi_hcho_total_col(key, pressure_1, scat_wt_1, kend_1, nlayer_1)

deallocate(pressure_1, scat_wt_1)

end subroutine read_omi_hcho_total_col

!-------------------------------------------------------------------------------

subroutine write_omi_hcho_total_col(key, ifile, fform)

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
call write_r8_array(  ifile, nlayer(key),  pressure(key,:), fileformat,'pressure')
call write_r8_array(  ifile, nlayer(key),  scat_wt(key,:), fileformat,'scat_wt')
call write_int_scalar(ifile,                             key, fileformat,'key')

end subroutine write_omi_hcho_total_col

!-------------------------------------------------------------------------------

subroutine interactive_omi_hcho_total_col(key)

integer, intent(out) :: key

if ( .not. module_initialized ) call initialize_module

! STOP because routine is not finished.
write(string1,*)'interactive_omi_hcho_total_col not yet working.'
call error_handler(E_ERR, 'interactive_omi_hcho_total_col', string1, source)

if(num_omi_hcho_obs >= max_omi_hcho_obs) then
   write(string1, *)'Not enough space for an omi hcho obs.'
   write(string2, *)'Can only have max_omi_hcho_obs (currently ',max_omi_hcho_obs,')'
   call error_handler(E_ERR, 'interactive_omi_hcho_total_col', string1, &
              source, text2=string2)
endif

! Increment the index
num_omi_hcho_obs = num_omi_hcho_obs + 1
key            = num_omi_hcho_obs

! Otherwise, prompt for input for the three required beasts

write(*, *) 'Creating an interactive_omi_hcho_total_col observation'
write(*, *) 'This featue is not setup '

end subroutine interactive_omi_hcho_total_col

!-------------------------------------------------------------------------------

subroutine get_expected_omi_hcho_total_col(state_handle, ens_size, location, key, obs_time, expct_val, istatus)

   type(ensemble_type), intent(in)  :: state_handle
   type(location_type), intent(in)  :: location
   type(time_type),     intent(in)  :: obs_time
   integer,             intent(in)  :: ens_size
   integer,             intent(in)  :: key
   integer,             intent(out) :: istatus(:)
   real(r8),            intent(out) :: expct_val(:)
   
   character(len=*), parameter :: routine = 'get_expected_omi_hcho_total_col'
   character(len=120)          :: data_file
   character(len=*),parameter  :: fld = 'HCHO_VMR_inst'
   type(location_type) :: loc2
   
   integer :: layer_omi,level_omi
   integer :: layer_mdl,level_mdl
   integer :: k,kend_omi,imem,imemm
   integer :: interp_new
   integer :: date_obs,datesec_obs,pbl_index
   integer, dimension(ens_size) :: zstatus
   
   real(r8) :: eps, AvogN, Rd, Ru, Cp, grav, msq2cmsq
   real(r8) :: missing,hcho_min,tmp_max
   real(r8) :: level
   real(r8) :: tmp_vir_k, tmp_vir_kp
   real(r8) :: mloc(3)
   real(r8) :: hcho_val_conv
   real(r8) :: up_wt,dw_wt,tl_wt,lnpr_mid
   real(r8), dimension(ens_size) :: hcho_mdl_1, tmp_mdl_1, qmr_mdl_1, prs_mdl_1
   real(r8), dimension(ens_size) :: hcho_mdl_n, tmp_mdl_n, qmr_mdl_n, prs_mdl_n
   real(r8), dimension(ens_size) :: prs_sfc
   
   real(r8), allocatable, dimension(:)   :: thick, prs_omi, prs_omi_mem
   real(r8), allocatable, dimension(:,:) :: hcho_val, tmp_val, qmr_val
   logical  :: return_now,hcho_return_now,tmp_return_now,qmr_return_now
!
! Upper BC variables
   integer  :: ncnt,kstart,flg
   real(r8) :: lon_obs,lat_obs,pi,rad2deg
   real     :: prs_del,delta,bdy_coef
   real     :: hcho_bot,prs_bot,tmp_bot,qmr_bot
   real     :: hcho_top,prs_top,tmp_top,qmr_top
   real(r8), allocatable, dimension(:)   :: hcho_prf_mdl,tmp_prf_mdl,qmr_prf_mdl
   real(r8), allocatable, dimension(:)   :: prs_omi_top
   
   if ( .not. module_initialized ) call initialize_module
   
   pi       = 4.*atan(1.)
   rad2deg  = 360./(2.*pi)
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
   bdy_coef = 0.95
   prs_del  = 1000.         ! Pa  
   if(use_log_hcho) then
      hcho_min = log(hcho_min)
   endif

! Assign vertical grid information

   layer_omi = nlayer(key)
   level_omi = nlayer(key)+1
   kend_omi  = kend(key)
   layer_mdl = nlayer_model
   level_mdl = nlayer_model+1
   
   allocate(prs_omi(layer_omi))
   allocate(prs_omi_mem(layer_omi))
   prs_omi(1:layer_omi)=pressure(key,1:layer_omi)*100.
   
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
   zstatus(:)=0
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
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_1, zstatus) ! K 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_1, zstatus) ! kg / kg 
      zstatus(:)=0
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
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_n, zstatus) ! K 
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_n, zstatus) ! kg / kg 
      zstatus(:)=0
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

! Get values at OMI pressure levels

   allocate(hcho_val(ens_size,layer_omi))
   allocate(tmp_val(ens_size,layer_omi))
   allocate(qmr_val(ens_size,layer_omi))

   do k=1,layer_omi
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), prs_omi(k), VERTISPRESSURE)
      call interpolate(state_handle, ens_size, loc2, QTY_HCHO, hcho_val(:,k), zstatus)  
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_val(:,k), zstatus)  
      zstatus(:)=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_val(:,k), zstatus)  

! Correcting for expected failures near the surface
      do imem=1,ens_size
         if (prs_omi(k).ge.prs_mdl_1(imem)) then
            hcho_val(imem,k) = hcho_mdl_1(imem)
            tmp_val(imem,k) = tmp_mdl_1(imem)
            qmr_val(imem,k) = qmr_mdl_1(imem)
            cycle
         endif

! Correcting for expected failures near the top
         if (prs_omi(k).le.prs_mdl_n(imem)) then
            hcho_val(imem,k) = hcho_mdl_n(imem)
            tmp_val(imem,k) = tmp_mdl_n(imem)
            qmr_val(imem,k) = qmr_mdl_n(imem)
            cycle
         endif
      enddo
!
!      write(string1, *)'APM: hcho ',k,layer_omi,prs_omi(k),hcho_val(1,k)
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
            write(string1, *) 'APM: Input data has missing values '
            call error_handler(E_MSG, routine, string1, source)
            call track_status(ens_size, zstatus, expct_val, istatus, return_now)
!            do imemm=1,ens_size
!               write(string1, *) &
!               'APM: Model profile values: hcho,tmp,qmr',key,imem,k,hcho_val(imemm,k), &
!               tmp_val(imemm,k),qmr_val(imemm,k)     
!               call error_handler(E_ALLMSG, routine, string1, source)
!            enddo
            return
         endif
      enddo
!
! Convert units for hcho from ppmv
      hcho_val(:,k) = hcho_val(:,k) * 1.e-6_r8
      hcho_mdl_1(:)=hcho_mdl_1(:) * 1.e-6_r8
      hcho_mdl_n(:)=hcho_mdl_n(:) * 1.e-6_r8
   enddo
!
! Use large scale hcho data above the regional model top
! OMI vertical is from bottom to top   
   kstart=-1
   do imem=1,ens_size
      if (prs_omi(layer_omi).lt.prs_mdl_n(imem)) then
         do k=1,layer_omi
            if (prs_omi(k).le.prs_mdl_n(imem)) then
               kstart=k-1
               exit
            endif
         enddo
         do k=kstart,layer_omi
             hcho_val(imem,k)=hcho_val(imem,kstart-1)
         enddo
!
!         ncnt=layer_omi-kstart+1
!         allocate(prs_omi_top(ncnt))
!         allocate(hcho_prf_mdl(ncnt),tmp_prf_mdl(ncnt),qmr_prf_mdl(ncnt))
!         do k=kstart,layer_omi
!            prs_omi_top(k-kstart+1)=prs_omi(k)
!         enddo
!         prs_omi_top(:)=prs_omi_top(:)/100.
!
!         lon_obs=mloc(1)/rad2deg
!         lat_obs=mloc(2)/rad2deg
!         call get_time(obs_time,datesec_obs,date_obs)
!
! Mozart does not have HCHO          
!         data_file=trim(upper_data_file)
!         model=trim(upper_data_model)
!         call get_upper_bdy_fld(fld,model,data_file,ls_chem_dx,ls_chem_dy, &
!         ls_chem_dz,ls_chem_dt,lon_obs,lat_obs,prs_omi_top, &
!         ncnt,hcho_prf_mdl,tmp_prf_mdl,qmr_prf_mdl,date_obs,datesec_obs)
!
! Impose ensemble perturbations from level kstart-1      
!         do k=kstart,layer_omi 
!            hcho_val(imem,k)=hcho_prf_mdl(k-kstart+1)*hcho_val(imem,kstart-1)/ &
!            (sum(hcho_val(:,kstart-1))/real(ens_size))
!            tmp_val(imem,k)=tmp_prf_mdl(k-kstart+1)*tmp_val(imem,kstart-1)/ &
!            (sum(tmp_val(:,kstart-1))/real(ens_size))
!            qmr_val(imem,k)=qmr_prf_mdl(k-kstart+1)*qmr_val(imem,kstart-1)/ &
!            (sum(qmr_val(:,kstart-1))/real(ens_size))
!         enddo
!         deallocate(prs_omi_top)
!         deallocate(hcho_prf_mdl,tmp_prf_mdl,qmr_prf_mdl)
      endif
   enddo
!
! Check full profile for negative values
   do imem=1,ens_size
      flg=0
      do k=1,layer_omi    
!         write(string1, *) &
!         'APM: hcho values: imem,k,hcho ',imem,k,hcho_val(imem,k)
!         call error_handler(E_MSG, routine, string1, source)
         if(hcho_val(imem,k).lt.0. .or. tmp_val(imem,k).lt.0. .or. &
         qmr_val(imem,k).lt.0.) then
            flg=1   
            write(string1, *) &
            'APM: Recentered full profile has negative values for key,imem ',key,imem
            call error_handler(E_ALLMSG, routine, string1, source)
            call track_status(ens_size, zstatus, expct_val, istatus, return_now)
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
   allocate(thick(layer_omi))

   do imem=1,ens_size
! Define upper and lower values for layer grid
! (OMI HCHO grid is bottom to top) prs is in Pa
      prs_omi_mem(:)=prs_omi(:)      
! Definitions for k=1 or k=layer_omi
      prs_bot=prs_sfc(imem)
      if (prs_bot.le.prs_omi_mem(1)) then
         prs_bot=prs_omi_mem(1)+prs_del
      endif   
! Bottom terms
      hcho_bot=hcho_mdl_1(imem)
      tmp_bot=tmp_mdl_1(imem)
      qmr_bot=qmr_mdl_1(imem)
! Top terms
      prs_top=prs_omi(layer_omi)+(prs_omi(layer_omi)-prs_omi(layer_omi-1))/2.
      if(prs_top.le.0.) prs_top=bdy_coef*prs_omi(layer_omi)
! hcho
      delta=(hcho_val(imem,layer_omi)-hcho_val(imem,layer_omi-1))/ &
      (prs_omi(layer_omi)-prs_omi(layer_omi-1))
      hcho_top=hcho_val(imem,layer_omi) + delta*(prs_top-prs_omi(layer_omi))
      if(hcho_top.le.0.) then
         if(delta.le.0.) hcho_top=bdy_coef*hcho_val(imem,layer_omi)
         if(delta.gt.0.) hcho_top=(2.-bdy_coef)*hcho_val(imem,layer_omi)
      endif
! tmp
      delta=(tmp_val(imem,layer_omi)-tmp_val(imem,layer_omi-1))/ &
      (prs_omi(layer_omi)-prs_omi(layer_omi-1))
      tmp_top=tmp_val(imem,layer_omi) + delta*(prs_top-prs_omi(layer_omi))
      if(tmp_top.le.0.) then
         if(delta.le.0.) tmp_top=bdy_coef*tmp_val(imem,layer_omi)
         if(delta.gt.0.) tmp_top=(2.-bdy_coef)*tmp_val(imem,layer_omi)
      endif
! qmr
      delta=(qmr_val(imem,layer_omi)-qmr_val(imem,layer_omi-1))/ &
      (prs_omi(layer_omi)-prs_omi(layer_omi-1))
      qmr_top=qmr_val(imem,layer_omi) + delta*(prs_top-prs_omi(layer_omi))
      if(qmr_top.le.0.) then
         if(delta.le.0.) qmr_top=bdy_coef*qmr_val(imem,layer_omi)
         if(delta.gt.0.) qmr_top=(2.-bdy_coef)*qmr_val(imem,layer_omi)
      endif
!
! VERTICAL SUMMATION      
! k=1 term
      k=1
         lnpr_mid=(log(prs_omi_mem(k+1))+log(prs_bot))/2.
         up_wt=log(prs_bot)-lnpr_mid
         dw_wt=lnpr_mid-log(prs_omi_mem(k+1))
         tl_wt=up_wt+dw_wt
         tmp_vir_k  = (1.0_r8 + eps*qmr_bot)*tmp_bot
         tmp_vir_kp = (1.0_r8 + eps*qmr_val(imem,k+1))*tmp_val(imem,k+1)
! thickness
         thick(k)   = Rd*(dw_wt*tmp_vir_k + up_wt*tmp_vir_kp)/tl_wt/grav* &
         log(prs_bot/prs_omi_mem(k+1))
! hcho term (convert VMR to molar density (mol/m^3))
         if(use_log_hcho) then
            hcho_val_conv = (dw_wt*exp(hcho_bot)+up_wt*exp(hcho_val(imem,k+1)))/tl_wt * &
            (dw_wt*prs_bot+up_wt*prs_omi(k+1)) / &
            (Ru*(dw_wt*tmp_bot+up_wt*tmp_val(imem,k+1)))
         else
            hcho_val_conv = (dw_wt*hcho_bot+up_wt*hcho_val(imem,k+1))/tl_wt * &
            (dw_wt*prs_bot+up_wt*prs_omi(k+1)) / &
            (Ru*(dw_wt*tmp_bot+up_wt*tmp_val(imem,k+1)))
         endif
! expected retrieval sum
         expct_val(imem) = expct_val(imem) + thick(k) * hcho_val_conv * &
         AvogN/msq2cmsq * scat_wt(key,k)
!         if(key.eq.1 .and. imem.eq.1) then
!            write(string1, *)'APM: expected retr ',k,expct_val(imem),thick(k), &
!            hcho_val_conv, scat_wt(key,k)
!            call error_handler(E_MSG, routine, string1, source)
!         endif
!
! k=layer_omi term
      k=layer_omi
         lnpr_mid=(log(prs_top)+log(prs_omi_mem(k)))/2.
         up_wt=log(prs_omi_mem(k))-lnpr_mid
         dw_wt=lnpr_mid-log(prs_top)
         tl_wt=up_wt+dw_wt
         tmp_vir_k  = (1.0_r8 + eps*qmr_val(imem,k))*tmp_val(imem,k)
         tmp_vir_kp = (1.0_r8 + eps*qmr_top)*tmp_top
! thickness
         thick(k)   = Rd*(dw_wt*tmp_vir_k + up_wt*tmp_vir_kp)/tl_wt/grav* &
         log(prs_omi_mem(k)/prs_top)
! hcho term (convert VMR to molar density (mol/m^3))
         if(use_log_hcho) then
            hcho_val_conv = (dw_wt*exp(hcho_val(imem,k))+up_wt*exp(hcho_top))/tl_wt * &
            (dw_wt*prs_omi(k)+up_wt*prs_top) / &
            (Ru*(dw_wt*tmp_val(imem,k)+up_wt*tmp_top))
         else
            hcho_val_conv = (dw_wt*hcho_val(imem,k)+up_wt*hcho_top)/tl_wt * &
            (dw_wt*prs_omi(k)+up_wt*prs_top) / &
            (Ru*(dw_wt*tmp_val(imem,k)+up_wt*tmp_top))
         endif
! expected retrieval sum
         expct_val(imem) = expct_val(imem) + thick(k) * hcho_val_conv * &
         AvogN/msq2cmsq * scat_wt(key,k)
!         if(key.eq.1 .and. imem.eq.1) then
!            write(string1, *)'APM: expected retr ',k,expct_val(imem),thick(k), &
!            hcho_val_conv, scat_wt(key,k)
!            call error_handler(E_MSG, routine, string1, source)
!         endif 
!
! remaining terms
      do k=2,layer_omi-1
         prs_bot=(prs_omi_mem(k-1)+prs_omi_mem(k))/2.
         prs_top=(prs_omi_mem(k)+prs_omi_mem(k+1))/2.
         hcho_bot=(hcho_val(imem,k-1)+hcho_val(imem,k))/2.
         hcho_top=(hcho_val(imem,k)+hcho_val(imem,k+1))/2.
         tmp_bot=(tmp_val(imem,k-1)+tmp_val(imem,k))/2.
         tmp_top=(tmp_val(imem,k)+tmp_val(imem,k+1))/2.
         qmr_bot=(qmr_val(imem,k-1)+qmr_val(imem,k))/2.
         qmr_top=(qmr_val(imem,k)+qmr_val(imem,k+1))/2.
         lnpr_mid=(log(prs_top)+log(prs_omi_mem(k)))/2.
         up_wt=log(prs_bot)-lnpr_mid
         dw_wt=lnpr_mid-log(prs_omi_mem(k+1))
         tl_wt=up_wt+dw_wt
         tmp_vir_k  = (1.0_r8 + eps*qmr_bot)*tmp_bot
         tmp_vir_kp = (1.0_r8 + eps*qmr_top)*tmp_top
! thickness
         thick(k)   = Rd*(dw_wt*tmp_vir_k + up_wt*tmp_vir_kp)/tl_wt/grav* &
         log(prs_bot/prs_top)
! hcho term (convert VMR to molar density (mol/m^3))
         if(use_log_hcho) then
            hcho_val_conv = (dw_wt*exp(hcho_bot)+up_wt*exp(hcho_top))/tl_wt * &
            (dw_wt*prs_bot+up_wt*prs_top) / (Ru*(dw_wt*tmp_bot+up_wt*tmp_top))
         else
            hcho_val_conv = (dw_wt*hcho_bot+up_wt*hcho_top)/tl_wt * &
            (dw_wt*prs_bot+up_wt*prs_top) / (Ru*(dw_wt*tmp_bot+up_wt*tmp_top))
         endif
! expected retrieval
         expct_val(imem) = expct_val(imem) + thick(k) * hcho_val_conv * &
         AvogN/msq2cmsq * scat_wt(key,k)
!         if(key.eq.1 .and. imem.eq.1) then
!            write(string1, *)'APM: expected retr ',k,expct_val(imem),thick(k), &
!            hcho_val_conv, scat_wt(key,k)
!            call error_handler(E_MSG, routine, string1, source)
!         endif
      enddo       
!      if(key.eq.1 .and. imem.eq.1) then
!         write(string1, *)'APM: FINAL EXPECTED VALUE ',expct_val(imem)
!         call error_handler(E_MSG, routine, string1, source)
!         write(string1, *)'  '
!         call error_handler(E_MSG, routine, string1, source)
!      endif

      if(expct_val(imem).lt.0) then
         zstatus(imem)=20
         expct_val(:)=missing_r8
         write(string1, *) &
         'APM NOTICE: OMI HCHO expected value is negative '
         call error_handler(E_MSG, routine, string1, source)
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
   enddo

! Clean up and return
   deallocate(hcho_val, tmp_val, qmr_val)
   deallocate(thick)
   deallocate(prs_omi, prs_omi_mem)

end subroutine get_expected_omi_hcho_total_col

!-------------------------------------------------------------------------------

subroutine set_obs_def_omi_hcho_total_col(key, hcho_pressure, hcho_scat_wt, hcho_kend, hcho_nlayer)

integer,                            intent(in)   :: key, hcho_nlayer
integer,                            intent(in)   :: hcho_kend
real(r8), dimension(hcho_nlayer),  intent(in)   :: hcho_pressure
real(r8), dimension(hcho_nlayer),    intent(in)   :: hcho_scat_wt

if ( .not. module_initialized ) call initialize_module

if(num_omi_hcho_obs >= max_omi_hcho_obs) then
   write(string1, *)'Not enough space for omi hcho total col obs.'
   write(string2, *)'Can only have max_omi_hcho_obs (currently ',max_omi_hcho_obs,')'
   call error_handler(E_ERR,'set_obs_def_omi_hcho_total_col',string1,source,revision, &
   revdate,text2=string2)
endif

nlayer(key) = hcho_nlayer
kend(key) = hcho_kend
pressure(key,1:hcho_nlayer) = hcho_pressure(1:hcho_nlayer)
scat_wt(key,1:hcho_nlayer) = hcho_scat_wt(1:hcho_nlayer)

end subroutine set_obs_def_omi_hcho_total_col

end module obs_def_omi_hcho_total_col_mod

! END DART PREPROCESS MODULE CODE
