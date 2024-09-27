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
! MAXDOAS_NO2_COLUMN, QTY_NO2
! END DART PREPROCESS KIND LIST
!
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_maxdoas_no2_mod, only : get_expected_maxdoas_no2, &
!                                  read_maxdoas_no2, &
!                                  write_maxdoas_no2, &
!                                  interactive_maxdoas_no2
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(MAXDOAS_NO2_COLUMN)                                                           
!         call get_expected_maxdoas_no2(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!
! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(MAXDOAS_NO2_COLUMN)
!         call read_maxdoas_no2(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(MAXDOAS_NO2_COLUMN)
!         call write_maxdoas_no2(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF
!
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(MAXDOAS_NO2_COLUMN)
!         call interactive_maxdoas_no2(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_maxdoas_no2_mod

use             types_mod, only : r8, MISSING_R8

use         utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                                  nmlfileunit, check_namelist_read, &
                                  find_namelist_in_file, do_nml_file, do_nml_term, &
                                  ascii_file_format

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

public :: write_maxdoas_no2, &
          read_maxdoas_no2, &
          interactive_maxdoas_no2, &
          get_expected_maxdoas_no2, &
          set_obs_def_maxdoas_no2

! Storage for the special information required for observations of this type
integer, parameter    :: max_maxdoas_no2_obs = 10000000
integer               :: num_maxdoas_no2_obs = 0
integer,  allocatable :: nlayer(:)
real(r8), allocatable :: prior(:)
real(r8), allocatable :: pressure(:,:)
real(r8), allocatable :: avg_kernel(:,:)

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_maxdoas_no2_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

character(len=512) :: string1, string2

logical, save :: module_initialized = .false.

! Namelist with default values
logical :: use_log_no2   = .false.
integer :: nlayer_model = -9999
integer :: nlayer_maxdoas_no2 = -9999

namelist /obs_def_MAXDOAS_NO2_nml/ use_log_no2, nlayer_model, nlayer_maxdoas_no2

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
call find_namelist_in_file("input.nml", "obs_def_MAXDOAS_NO2_nml", iunit)
read(iunit, nml = obs_def_MAXDOAS_NO2_nml, iostat = rc)
call check_namelist_read(iunit, rc, "obs_def_MAXDOAS_NO2_nml")

! Record the namelist values
if (do_nml_file()) write(nmlfileunit, nml=obs_def_MAXDOAS_NO2_nml)
if (do_nml_term()) write(     *     , nml=obs_def_MAXDOAS_NO2_nml)

! Check for valid values

if (nlayer_model < 1) then
   write(string1,*)'obs_def_MAXDOAS_NO2_nml:nlayer_model must be > 0, it is ',nlayer_model
   call error_handler(E_ERR,'initialize_module',string1,source)
endif

if (nlayer_maxdoas_no2 < 1) then
   write(string1,*)'obs_def_MAXDOAS_NO2_nml:nlayer_maxdoas_no2 must be > 0, it is ',nlayer_maxdoas_no2
   call error_handler(E_ERR,'initialize_module',string1,source)
endif

allocate(    nlayer(max_maxdoas_no2_obs))
allocate(    prior(max_maxdoas_no2_obs))
allocate(  pressure(max_maxdoas_no2_obs,nlayer_maxdoas_no2+1))
allocate(avg_kernel(max_maxdoas_no2_obs,nlayer_maxdoas_no2))

end subroutine initialize_module

!-------------------------------------------------------------------------------

subroutine read_maxdoas_no2(key, ifile, fform)

integer,          intent(out)          :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

! temporary arrays to hold buffer till we decide if we have enough room

integer               :: keyin
integer               :: nlayer_1
real(r8), allocatable :: prior_1
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
keyin = read_int_scalar(ifile, fileformat, 'keyin')

counts1 = counts1 + 1
key     = counts1

if(counts1 > max_maxdoas_no2_obs) then
   write(string1, *)'Not enough space for maxdoas no2 obs.'
   write(string2, *)'Can only have max_maxdoas_no2_obs (currently ',max_maxdoas_no2_obs,')'
   call error_handler(E_ERR,'read_maxdoas_no2',string1,source,text2=string2)
endif

call set_obs_def_maxdoas_no2(key, pressure_1, avg_kernel_1, prior_1, nlayer_1)

deallocate(pressure_1, avg_kernel_1)

end subroutine read_maxdoas_no2

!-------------------------------------------------------------------------------

subroutine write_maxdoas_no2(key, ifile, fform)

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
call write_r8_scalar(ifile,                     prior(key), fileformat,'prior')
call write_r8_array(  ifile, nlayer(key)+1,      pressure(key,:), fileformat,'pressure')
call write_r8_array(  ifile, nlayer(key),        avg_kernel(key,:), fileformat,'avg_kernel')
call write_int_scalar(ifile,                     key, fileformat,'key')

end subroutine write_maxdoas_no2

!-------------------------------------------------------------------------------

subroutine interactive_maxdoas_no2(key)

integer, intent(out) :: key

if ( .not. module_initialized ) call initialize_module

! STOP because routine is not finished.
write(string1,*)'interactive_maxdoas_no2 not yet working.'
call error_handler(E_ERR, 'interactive_maxdoas_no2', string1, source)

if(num_maxdoas_no2_obs >= max_maxdoas_no2_obs) then
   write(string1, *)'Not enough space for an maxdoas no2 obs.'
   write(string2, *)'Can only have max_maxdoas_no2_obs (currently ',max_maxdoas_no2_obs,')'
   call error_handler(E_ERR, 'interactive_maxdoas_no2', string1, &
              source, text2=string2)
endif

! Increment the index
num_maxdoas_no2_obs = num_maxdoas_no2_obs + 1
key            = num_maxdoas_no2_obs

! Otherwise, prompt for input for the three required beasts

write(*, *) 'Creating an interactive_maxdoas_no2 observation'
write(*, *) 'This featue is not setup '

end subroutine interactive_maxdoas_no2

!-------------------------------------------------------------------------------

subroutine get_expected_maxdoas_no2(state_handle, ens_size, location, key, expct_val, istatus)

   type(ensemble_type), intent(in)  :: state_handle
   type(location_type), intent(in)  :: location
   integer,             intent(in)  :: ens_size
   integer,             intent(in)  :: key
   integer,             intent(out) :: istatus(:)
   real(r8),            intent(out) :: expct_val(:)
   
   character(len=*), parameter :: routine = 'get_expected_maxdoas_no2'
   type(location_type) :: loc2
   
   integer :: layer_maxdoas,level_maxdoas
   integer :: layer_mdl,level_mdl
   integer :: k,imem,kend_maxdoas
   integer :: interp_new
   integer, dimension(ens_size) :: zstatus
   
   real(r8) :: eps, AvogN, Rd, Ru, Cp, grav, msq2cmsq
   real(r8) :: missing,no2_min, tmp_max
   real(r8) :: level,del_prs
   real(r8) :: tmp_vir_k, tmp_vir_kp
   real(r8) :: mloc(3)
   real(r8) :: no2_val_conv,amf
   real(r8) :: up_wt,dw_wt,tl_wt,lnpr_mid
   real(r8), dimension(ens_size) :: no2_mdl_1, tmp_mdl_1, qmr_mdl_1, prs_mdl_1
   real(r8), dimension(ens_size) :: no2_mdl_n, tmp_mdl_n, qmr_mdl_n, prs_mdl_n
   real(r8), dimension(ens_size) :: prs_sfc
   
   real(r8), allocatable, dimension(:)   :: thick, prs_maxdoas, prs_maxdoas_mem
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

   layer_maxdoas = nlayer(key)
   level_maxdoas = nlayer(key)+1
   layer_mdl = nlayer_model
   level_mdl = nlayer_model+1
   
   allocate(prs_maxdoas(level_maxdoas))
   allocate(prs_maxdoas_mem(level_maxdoas))
   prs_maxdoas(1:level_maxdoas)=pressure(key,1:level_maxdoas)

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
      call interpolate(state_handle, ens_size, loc2, QTY_NO2, no2_mdl_1, zstatus) ! ppmv 
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_1, zstatus) ! K 
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_1, zstatus) ! kg / kg 
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_1, zstatus) ! Pa
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

!   write(string1, *)'APM: no2 lower bound 1 ',no2_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: tmp lower bound 1 ',tmp_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: qmr lower bound 1 ',qmr_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: prs lower bound 1 ',prs_mdl_1
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
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_mdl_n, zstatus) ! K 
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_mdl_n, zstatus) ! kg / kg 
      zstatus=0
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

!   write(string1, *)'APM: no2 upper bound 1 ',no2_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: tmp upper bound 1 ',tmp_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: qmr upper bound 1 ',qmr_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: prs upper bound 1 ',prs_mdl_n
!   call error_handler(E_MSG, routine, string1, source)

! Get profiles at MAXDOAS pressure levels

   allocate(no2_val(ens_size,level_maxdoas))
   allocate(tmp_val(ens_size,level_maxdoas))
   allocate(qmr_val(ens_size,level_maxdoas))

   do k=1,level_maxdoas
      zstatus=0
      loc2 = set_location(mloc(1), mloc(2), prs_maxdoas(k), VERTISPRESSURE)
      call interpolate(state_handle, ens_size, loc2, QTY_NO2, no2_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_TEMPERATURE, tmp_val(:,k), zstatus)  
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_VAPOR_MIXING_RATIO, qmr_val(:,k), zstatus)  
!
! Correcting for expected failures near the surface
      do imem=1,ens_size
         if (prs_maxdoas(k).ge.prs_mdl_1(imem)) then
            no2_val(imem,k) = no2_mdl_1(imem)
            tmp_val(imem,k) = tmp_mdl_1(imem)
            qmr_val(imem,k) = qmr_mdl_1(imem)
            cycle
         endif

! Correcting for expected failures near the top
         if (prs_maxdoas(k).le.prs_mdl_n(imem)) then
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

   istatus=0
   zstatus(:)=0
   expct_val(:)=0.0
   allocate(thick(layer_maxdoas))

   do imem=1,ens_size
! Adjust the MAXDOAS pressure for WRF-Chem lower/upper boudary pressure
! (MAXDOAS NO2 vertical grid is bottom to top)
      prs_maxdoas_mem(:)=prs_maxdoas(:)
      if (prs_sfc(imem).gt.prs_maxdoas_mem(1)) then
         prs_maxdoas_mem(1)=prs_sfc(imem)
      endif   

! Calculate the thicknesses

      do k=1,kend_maxdoas
         lnpr_mid=(log(prs_maxdoas_mem(k+1))+log(prs_maxdoas_mem(k)))/2.
         up_wt=log(prs_maxdoas_mem(k))-lnpr_mid
         dw_wt=lnpr_mid-log(prs_maxdoas_mem(k+1))
         tl_wt=up_wt+dw_wt      
         tmp_vir_k  = (1.0_r8 + eps*qmr_val(imem,k))*tmp_val(imem,k)
         tmp_vir_kp = (1.0_r8 + eps*qmr_val(imem,k+1))*tmp_val(imem,k+1)
         thick(k)   = Rd*(dw_wt*tmp_vir_k + up_wt*tmp_vir_kp)/tl_wt/grav* &
         log(prs_maxdoas_mem(k)/prs_maxdoas_mem(k+1))
      enddo

! Process the vertical summation

      do k=1,kend_maxdoas
         lnpr_mid=(log(prs_maxdoas_mem(k+1))+log(prs_maxdoas_mem(k)))/2.
         up_wt=log(prs_maxdoas_mem(k))-lnpr_mid
         dw_wt=lnpr_mid-log(prs_maxdoas_mem(k+1))
         tl_wt=up_wt+dw_wt

! Convert from VMR to molar density (mol/m^3)
         if(use_log_no2) then
            no2_val_conv = (dw_wt*exp(no2_val(imem,k))+up_wt*exp(no2_val(imem,k+1)))/tl_wt * &
            (dw_wt*prs_maxdoas_mem(k)+up_wt*prs_maxdoas_mem(k+1)) / &
            (Ru*(dw_wt*tmp_val(imem,k)+up_wt*tmp_val(imem,k+1)))
         else
            no2_val_conv = (dw_wt*no2_val(imem,k)+up_wt*no2_val(imem,k+1))/tl_wt * &
            (dw_wt*prs_maxdoas_mem(k)+up_wt*prs_maxdoas_mem(k+1)) / &
            (Ru*(dw_wt*tmp_val(imem,k)+up_wt*tmp_val(imem,k+1)))
         endif
 
! Get expected observation

         expct_val(imem) = expct_val(imem) + thick(k) * no2_val_conv * &
         avg_kernel(key,k)
      enddo

      if(expct_val(imem).lt.0) then
         zstatus(imem)=20
         expct_val(:)=missing_r8
         write(string1, *) &
         'APM NOTICE: MAXDOAS NO2 expected value is negative '
         call error_handler(E_MSG, routine, string1, source)
         call track_status(ens_size, zstatus, expct_val, istatus, return_now)
         return
      endif
   enddo

! Clean up and return
   deallocate(no2_val, tmp_val, qmr_val)
   deallocate(thick)
   deallocate(prs_maxdoas, prs_maxdoas_mem)

end subroutine get_expected_maxdoas_no2

!-------------------------------------------------------------------------------

subroutine set_obs_def_maxdoas_no2(key, no2_pressure, no2_avg_kernel, no2_prior,  no2_nlayer)

integer,                           intent(in)   :: key, no2_nlayer
real(r8),                          intent(in)   :: no2_prior
real(r8), dimension(no2_nlayer+1),  intent(in)   :: no2_pressure
real(r8), dimension(no2_nlayer),    intent(in)   :: no2_avg_kernel

if ( .not. module_initialized ) call initialize_module

if(num_maxdoas_no2_obs >= max_maxdoas_no2_obs) then
   write(string1, *)'Not enough space for maxdoas no2 obs.'
   write(string2, *)'Can only have max_maxdoas_no2_obs (currently ',max_maxdoas_no2_obs,')'
   call error_handler(E_ERR,'set_obs_def_maxdoas_no2',string1,source,revision, &
   revdate,text2=string2)
endif

nlayer(key) = no2_nlayer
prior(key) = no2_prior
pressure(key,1:no2_nlayer+1) = no2_pressure(1:no2_nlayer+1)
avg_kernel(key,1:no2_nlayer) = no2_avg_kernel(1:no2_nlayer)

end subroutine set_obs_def_maxdoas_no2

!-------------------------------------------------------------------------------

function read_int_scalar(ifile, fform, context)

integer                      :: read_int_scalar
integer,          intent(in) :: ifile
character(len=*), intent(in) :: fform
character(len=*), intent(in) :: context

integer :: io

if (ascii_file_format(fform)) then
   read(ifile, *, iostat = io) read_int_scalar
else
   read(ifile, iostat = io) read_int_scalar
endif
if ( io /= 0 ) then
   call error_handler(E_ERR,'read_int_scalar', context, source)
endif

end function read_int_scalar

!-------------------------------------------------------------------------------

subroutine write_int_scalar(ifile, my_scalar, fform, context)

integer,          intent(in) :: ifile
integer,          intent(in) :: my_scalar
character(len=*), intent(in) :: fform
character(len=*), intent(in) :: context

integer :: io

if (ascii_file_format(fform)) then
   write(ifile, *, iostat=io) my_scalar
else
   write(ifile, iostat=io) my_scalar
endif
if ( io /= 0 ) then
   call error_handler(E_ERR, 'write_int_scalar', context, source)
endif

end subroutine write_int_scalar

!-------------------------------------------------------------------------------

function read_r8_scalar(ifile, fform, context)

real(r8)                     :: read_r8_scalar
integer,          intent(in) :: ifile
character(len=*), intent(in) :: fform
character(len=*), intent(in) :: context

integer :: io

if (ascii_file_format(fform)) then
   read(ifile, *, iostat = io) read_r8_scalar
else
   read(ifile, iostat = io) read_r8_scalar
endif
if ( io /= 0 ) then
   call error_handler(E_ERR,'read_r8_scalar', context, source)
endif

end function read_r8_scalar

!-------------------------------------------------------------------------------

subroutine write_r8_scalar(ifile, my_scalar, fform, context)

integer,          intent(in) :: ifile
real(r8),         intent(in) :: my_scalar
character(len=*), intent(in) :: fform
character(len=*), intent(in) :: context

integer :: io

if (ascii_file_format(fform)) then
   write(ifile, *, iostat=io) my_scalar
else
   write(ifile, iostat=io) my_scalar
endif
if ( io /= 0 ) then
   call error_handler(E_ERR, 'write_r8_scalar', context, source)
endif

end subroutine write_r8_scalar

!-------------------------------------------------------------------------------

subroutine read_r8_array(ifile, num_items, r8_array, fform, context)

integer,          intent(in)  :: ifile, num_items
real(r8),         intent(out) :: r8_array(:)
character(len=*), intent(in)  :: fform
character(len=*), intent(in)  :: context

integer :: io

if (ascii_file_format(fform)) then
   read(ifile, *, iostat = io) r8_array(1:num_items)
else
   read(ifile, iostat = io) r8_array(1:num_items)
endif
if ( io /= 0 ) then
   call error_handler(E_ERR, 'read_r8_array', context, source)
endif

end subroutine read_r8_array

!-------------------------------------------------------------------------------

subroutine write_r8_array(ifile, num_items, array, fform, context)

integer,          intent(in) :: ifile, num_items
real(r8),         intent(in) :: array(:)
character(len=*), intent(in) :: fform
character(len=*), intent(in) :: context

integer :: io

if (ascii_file_format(fform)) then
   write(ifile, *, iostat = io) array(1:num_items)
else
   write(ifile, iostat = io) array(1:num_items)
endif
if ( io /= 0 ) then
   call error_handler(E_ERR, 'write_r8_array', context, source)
endif

end subroutine write_r8_array

!-------------------------------------------------------------------------------




end module obs_def_maxdoas_no2_mod

! END DART PREPROCESS MODULE CODE
