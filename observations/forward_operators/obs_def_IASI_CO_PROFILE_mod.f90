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

! BEGIN DART PREPROCESS KIND LIST
! IASI_CO_PROFILE, QTY_CO
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_IASI_CO_PROFILE_mod, only : write_iasi_co_profile, read_iasi_co_profile, &
                                    interactive_iasi_co_profile, get_expected_iasi_co_profile, &
                                    set_obs_def_iasi_co_profile
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(IASI_CO_PROFILE)                                                           
!          call get_expected_iasi_co_profile(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(IASI_CO_PROFILE)
!          call read_iasi_co_profile(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(IASI_CO_PROFILE)
!          call write_iasi_co_profile(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(IASI_CO_PROFILE)
!          call interactive_iasi_co_profile(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!
! BEGIN DART PREPROCESS SET_OBS_DEF_IASI_CO_PROFILE
!      case(IASI_CO_PROFILE)
!          call set_obs_def_iasi_co_profile(obs_def%key)
! END DART PREPROCESS SET_OBS_DEF_IASI_CO_PROFILE
!
! BEGIN DART PREPROCESS MODULE CODE
module obs_def_IASI_CO_PROFILE_mod

   use        types_mod, only : r8, missing_r8
   use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             nmlfileunit, check_namelist_read, &
                             find_namelist_in_file, do_nml_file, do_nml_term, &
                             ascii_file_format, &
                             read_int_scalar, &
                             write_int_scalar, &       
                             read_r8_scalar, &
                             write_r8_scalar, &
                             read_r8_array, &
                             write_r8_array

   use     location_mod, only : location_type, set_location, get_location, VERTISPRESSURE, VERTISLEVEL, VERTISSURFACE, VERTISUNDEF

   use  assim_model_mod, only : interpolate
   use    obs_kind_mod, only  : QTY_CO, QTY_SURFACE_PRESSURE, QTY_PRESSURE, QTY_LANDMASK
   use ensemble_manager_mod,  only : ensemble_type
   use obs_def_utilities_mod, only : track_status

   implicit none
   private

   public :: write_iasi_co_profile, &
          read_iasi_co_profile, &
          interactive_iasi_co_profile, &
          get_expected_iasi_co_profile, &
          set_obs_def_iasi_co_profile

! Storage for the special information required for observations of this type
   integer, parameter               :: MAX_IASI_CO_OBS = 10000000
   integer, parameter               :: IASI_DIM = 19
   integer, parameter               :: IASI_DIMP = 20
   integer                          :: num_iasi_co_obs = 0
!
   real(r8), allocatable, dimension(:,:) :: avg_kernel
   real(r8), allocatable, dimension(:,:) :: pressure
   real(r8), allocatable, dimension(:) :: iasi_prior
   real(r8), allocatable, dimension(:) :: iasi_psurf
   integer,  allocatable, dimension(:) :: iasi_nlayers
   integer,  allocatable, dimension(:) :: iasi_nlevels

! version controlled file description for error handling, do not edit
   character(len=*), parameter :: source   = 'obs_def_IASI_CO_PROFILE_mod.f90'
   character(len=*), parameter :: revision = ''
   character(len=*), parameter :: revdate  = ''

   character(len=512) :: string1, string2
   character(len=200) :: upper_data_file
   character(len=200) :: upper_data_model
   character(len=200) :: model
   integer            :: ls_chem_dx, ls_chem_dy, ls_chem_dz, ls_chem_dt

   logical, save :: module_initialized = .false.
   integer  :: counts1 = 0

   character(len=129)  :: IASI_CO_retrieval_type
   logical             :: use_log_co
   integer :: nlayer_model = -9999
   integer :: nlayer_iasi = -9999
   integer :: nlayer_iasi_co_total_col = -9999
   integer :: nlayer_iasi_co_trop_col = -9999
   integer :: nlayer_iasi_co_profile = -9999
!
! IASI_CO_retrieval_type:
!     RAWR - retrievals in format from supplier
!     RETR - retrievals in retrieval (ppbv) format
!     QOR  - quasi-optimal retrievals
!     CPSR - compact phase space retrievals

   namelist /obs_def_IASI_CO_nml/upper_data_file, use_log_co, nlayer_model, &
   nlayer_iasi_co_total_col, nlayer_iasi_co_trop_col, nlayer_iasi_co_profile, &
   ls_chem_dx, ls_chem_dy, ls_chem_dz, ls_chem_dt, upper_data_model
   
contains

!----------------------------------------------------------------------

subroutine initialize_module

   integer :: iunit, rc

! Prevent multiple calls from executing this code more than once.
   if (module_initialized) return

   call register_module(source, revision, revdate)
   module_initialized = .true.

   allocate (avg_kernel(   MAX_IASI_CO_OBS,IASI_DIM))
   allocate (pressure(     MAX_IASI_CO_OBS,IASI_DIMP))
   allocate (iasi_prior(   MAX_IASI_CO_OBS))
   allocate (iasi_psurf(   MAX_IASI_CO_OBS))
   allocate (iasi_nlayers( MAX_IASI_CO_OBS))
   allocate (iasi_nlevels(MAX_IASI_CO_OBS))

! Read the namelist entry.
   IASI_CO_retrieval_type='RAWR'
   use_log_co=.false.
   call find_namelist_in_file("input.nml", "obs_def_IASI_CO_nml", iunit)
   read(iunit, nml = obs_def_IASI_CO_nml, iostat = rc)
   call check_namelist_read(iunit, rc, "obs_def_IASI_CO_nml")

! Record the namelist values used for the run ... 
   if (do_nml_file()) write(nmlfileunit, nml=obs_def_IASI_CO_nml)
   if (do_nml_term()) write(     *     , nml=obs_def_IASI_CO_nml)
   nlayer_iasi=nlayer_iasi_co_profile

!
end subroutine initialize_module
!
subroutine read_iasi_co_profile(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine read_iasi_co_profile(key, ifile, fform)

integer,          intent(out)          :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=32)              :: fileformat

integer                        :: iasi_nlayers_1
integer                        :: iasi_nlevels_1
real(r8)                       :: iasi_prior_1
real(r8)                       :: iasi_psurf_1
real(r8), dimension(IASI_DIM)  :: avg_kernels_1
real(r8), dimension(IASI_DIMP) :: pressure_1
integer                        :: keyin

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading
avg_kernels_1(:) = 0.0_r8
pressure_1(:) = 0.0_r8
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   iasi_nlayers_1 = read_iasi_nlayers(ifile, fileformat)
   iasi_nlevels_1 = iasi_nlayers_1+1
   iasi_prior_1 = read_iasi_prior(ifile, fileformat)
   iasi_psurf_1 = read_iasi_psurf(ifile, fileformat)
   avg_kernels_1(:) = read_iasi_avg_kernels(ifile, iasi_nlayers_1, fileformat)
   pressure_1(:) = read_iasi_pressure(ifile, iasi_nlevels_1, fileformat)
   read(ifile) keyin
   CASE DEFAULT
   iasi_nlayers_1 = read_iasi_nlayers(ifile, fileformat)
   iasi_nlevels_1 = iasi_nlayers_1+1
   iasi_prior_1 = read_iasi_prior(ifile, fileformat)
   iasi_psurf_1 = read_iasi_psurf(ifile, fileformat)
   avg_kernels_1(:) = read_iasi_avg_kernels(ifile, iasi_nlayers_1, fileformat)
   pressure_1(:) = read_iasi_pressure(ifile, iasi_nlevels_1, fileformat)
   read(ifile, *) keyin
END SELECT

counts1 = counts1 + 1
key = counts1
call set_obs_def_iasi_co_profile(key, avg_kernels_1, pressure_1, iasi_prior_1, iasi_psurf_1, &
                           iasi_nlayers_1, iasi_nlevels_1)
end subroutine read_iasi_co_profile
!
subroutine write_iasi_co_profile(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine write_iasi_co_profile(key, ifile, fform)

integer,          intent(in)           :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat
real(r8), dimension(IASI_DIM)   :: avg_kernels_temp
real(r8), dimension(IASI_DIMP)  :: pressure_temp
if ( .not. module_initialized ) call initialize_module
fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))
!
! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading
avg_kernels_temp(:)=avg_kernel(key,:)
pressure_temp(:)=pressure(key,:)
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   call write_iasi_nlayers(ifile, iasi_nlayers(key), fileformat)
   call write_iasi_prior(ifile, iasi_prior(key), fileformat)
   call write_iasi_psurf(ifile, iasi_psurf(key), fileformat)
   call write_iasi_avg_kernels(ifile, avg_kernels_temp, iasi_nlayers(key), fileformat)
   call write_iasi_pressure(ifile, pressure_temp, iasi_nlevels(key), fileformat)
   write(ifile) key
   CASE DEFAULT
   call write_iasi_nlayers(ifile, iasi_nlayers(key), fileformat)
   call write_iasi_prior(ifile, iasi_prior(key), fileformat)
   call write_iasi_psurf(ifile, iasi_psurf(key), fileformat)
   call write_iasi_avg_kernels(ifile, avg_kernels_temp, iasi_nlayers(key), fileformat)
   call write_iasi_pressure(ifile, pressure_temp, iasi_nlevels(key), fileformat)
   write(ifile, *) key
END SELECT 
end subroutine write_iasi_co_profile
!
subroutine interactive_iasi_co_profile(key)
!----------------------------------------------------------------------
!subroutine interactive_iasi_co_profile(key)
!
! Initializes the specialized part of a IASI observation
! Passes back up the key for this one
integer, intent(out) :: key

if ( .not. module_initialized ) call initialize_module
!
! Make sure there's enough space, if not die for now (clean later)
if(num_iasi_co_obs >= MAX_IASI_CO_OBS) then
   write(string1, *)'Not enough space for a iasi CO obs.'
   write(string2, *)'Can only have max_iasi_co obs (currently ',MAX_IASI_CO_OBS,')'
   call error_handler(E_ERR,'interactive_iasi_co_profile',string1,source,revision,revdate,text2=string2)
endif
!
! Increment the index
num_iasi_co_obs = num_iasi_co_obs + 1
key = num_iasi_co_obs
!
! Otherwise, prompt for input for the three required beasts
write(*, *) 'Creating an interactive_iasi_co_profile observation'
write(*, *) 'Input the IASI Prior '
read(*, *) iasi_prior
write(*, *) 'Input IASI Surface Pressure '
read(*, *) iasi_psurf(num_iasi_co_obs)
write(*, *) 'Input the 19 Averaging Kernel Weights '
read(*, *) avg_kernel(num_iasi_co_obs,:)
write(*, *) 'Input the 20 Averaging Pressure Levels '
read(*, *) pressure(num_iasi_co_obs,:)
end subroutine interactive_iasi_co_profile
!
subroutine get_expected_iasi_co_profile(state_handle, ens_size, location, key, val, istatus)
!----------------------------------------------------------------------
!subroutine get_expected_iasi_co_profile(state_handle, ens_size, location, key, val, istatus)
   type(ensemble_type), intent(in)  :: state_handle
   integer,             intent(in)  :: ens_size
   type(location_type), intent(in)  :: location
   integer,             intent(in)  :: key
   real(r8),            intent(out) :: val(ens_size)
   integer,             intent(out) :: istatus(ens_size)
!
   type(location_type) :: loc2

   integer             :: i, k, imem, kstr
   integer             :: interp_new
   integer             :: zstatus(ens_size)
   integer             :: layer_iasi,level_iasi
   integer             :: layer_mdl,level_mdl

   real(r8)            :: prs_iasi_sfc
   real(r8)            :: co_min, co_min_log, level
   real(r8)            :: mloc(3)
   real(r8),dimension(ens_size)  :: prs_mdl_sfc
   real(r8),dimension(ens_size)  :: prs_mdl_1, co_mdl_1
   real(r8),dimension(ens_size)  :: prs_mdl_n, co_mdl_n
   real(r8),dimension(:),allocatable   :: prs_iasi
   real(r8),dimension(:,:),allocatable :: co_val

   character(len=*), parameter :: routine = 'get_expected_iasi_co_profile'
   logical :: return_now
!
! Initialize DART
   istatus(:) = 0
   if ( .not. module_initialized ) call initialize_module
!
! Initialize variables (IASI is ppbv; WRF CO is ppmv)
   co_min      = 1.e-2
   co_min_log  = log(co_min)

   layer_iasi = iasi_nlayers(key)
   level_iasi = iasi_nlevels(key)
   layer_mdl = nlayer_model
   level_mdl = layer_mdl+1

   allocate(prs_iasi(layer_iasi))
   do k=1,layer_iasi
      prs_iasi(k)=(pressure(key,k)+pressure(key,k+1))/2.
   enddo
   
   if ( use_log_co ) then
      co_min=co_min_log
   endif
!
! Get location infomation
   mloc = get_location(location)
   if (mloc(2)>90.0_r8) then
      mloc(2)=90.0_r8
   elseif (mloc(2)<-90.0_r8) then
      mloc(2)=-90.0_r8
   endif
!
! pressure at model surface (Pa)
   level=0.0_r8
   loc2 = set_location(mloc(1), mloc(2), level, VERTISSURFACE)
   zstatus(:)=0
   call interpolate(state_handle, ens_size, loc2, QTY_SURFACE_PRESSURE, prs_mdl_sfc, zstatus)

   co_mdl_1(:)=missing_r8
   prs_mdl_1(:)=missing_r8

   do k=1,layer_mdl
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)
      call interpolate(state_handle, ens_size, loc2, QTY_CO, co_mdl_1, zstatus) ! ppmv 
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_1, zstatus) ! Pa
!
      interp_new=0
      do imem=1,ens_size
         if(co_mdl_1(imem).eq.missing_r8 .or. prs_mdl_1(imem).eq.missing_r8) then
            interp_new=1
            exit
         endif
      enddo
      if(interp_new.eq.0) then
         exit
      endif    
   enddo

!   write(string1, *)'APM: co lower bound ',co_mdl_1
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: prs lower bound ',prs_mdl_1
!   call error_handler(E_MSG, routine, string1, source)

   co_mdl_n(:)=missing_r8
   prs_mdl_n(:)=missing_r8

   do k=layer_mdl,1,-1
      level=real(k)
      zstatus(:)=0
      loc2 = set_location(mloc(1), mloc(2), level, VERTISLEVEL)  
      call interpolate(state_handle, ens_size, loc2, QTY_CO, co_mdl_n, zstatus) ! ppmv
      zstatus=0
      call interpolate(state_handle, ens_size, loc2, QTY_PRESSURE, prs_mdl_n, zstatus) ! Pa
!
      interp_new=0
      do imem=1,ens_size
         if(co_mdl_n(imem).eq.missing_r8 .or. prs_mdl_n(imem).eq.missing_r8) then
            interp_new=1
            exit
         endif
      enddo
      if(interp_new.eq.0) then
         exit
      endif    
   enddo

!   write(string1, *)'APM: co upper bound ',co_mdl_n
!   call error_handler(E_MSG, routine, string1, source)
!   write(string1, *)'APM: prs upper bound ',prs_mdl_n
!   call error_handler(E_MSG, routine, string1, source)

! Get profiles at IASI pressure levels
!    if (prs_mdl_sfc(imem).gt.iasi_psurf(key)) then
!       prs_iasi(1)=(prs_mdl_sfc(imem)+pressure(key,2))/2.
!    endif   
!
   allocate(co_val(ens_size,layer_iasi))
   do k=1,layer_iasi
      zstatus=0
      loc2 = set_location(mloc(1), mloc(2), prs_iasi(k), VERTISPRESSURE)
      call interpolate(state_handle, ens_size, loc2, QTY_CO, co_val(:,k), zstatus)  
!
! Correcting for expected failures near the surface
      do imem=1,ens_size
         if (prs_iasi(k).ge.prs_mdl_1(imem)) then
            co_val(imem,k) = co_mdl_1(imem)
            cycle
         endif

! Correcting for expected failures near the top
         if (prs_iasi(k).le.prs_mdl_n(imem)) then
            co_val(imem,k) = co_mdl_n(imem)
            cycle
         endif
      enddo
!
!      write(string1, *)'APM: co ',k,co_val(1,k)
!      call error_handler(E_MSG, routine, string1, source)
!
! Check data for missing values      
      do imem=1,ens_size
         if(co_val(imem,k).eq.missing_r8) then
            zstatus(:)=20
            val(:)=missing_r8
            write(string1, *) &
            'APM: Input data has missing values '
            call error_handler(E_MSG, routine, string1, source)
            call track_status(ens_size, zstatus, val, istatus, return_now)
            return
         endif
!
! Convert units for co to ppbv for IASI
         if(use_log_co) then
            co_val(imem,k) = log(exp(co_val(imem,k)) * 1.e3_r8)
         else
            co_val(imem,k) = co_val(imem,k) * 1.e3_r8
         endif
      enddo
   enddo
!
! Apply IASI Averaging kernel A and IASI Prior (I-A)xa
! x = Axm + (I-A)xa , where x is a 19 element vector 
!
! loop through IASI layers
   val(:) = 0.0_r8
   do imem=1,ens_size
      do k=1,layer_iasi
         if( use_log_co ) then
            val(imem) = val(imem) + avg_kernel(key,k) * exp(co_val(imem,k))  
         else
            val(imem) = val(imem) + avg_kernel(key,k) * co_val(imem,k)
!            if(imem.eq.1) then
!               write(string1, *) 'APM: imem,k,val,avgk,co_val ',imem,k,val(imem), &
!               avg_kernel(key,k),co_val(imem,k)
!               call error_handler(E_MSG, routine, string1, source)
!            endif
         endif         
      enddo
!
! NOTE: For the following the iasi_prior is zero due to the QOR subtraction
      if (trim(IASI_CO_retrieval_type).eq.'RAWR' .or. trim(IASI_CO_retrieval_type).eq.'QOR' &
      .or. trim(IASI_CO_retrieval_type).eq.'CPSR') then
         val(imem) = val(imem) + iasi_prior(key)
      elseif (trim(IASI_CO_retrieval_type).eq.'RETR') then
         val(imem) = val(imem) + iasi_prior(key)
      endif

      if(val(imem).lt.0) then
         zstatus(imem)=20
         val(:)=missing_r8
         write(string1, *) &
         'APM NOTICE: IASI CO expected value is negative '
         call error_handler(E_MSG, routine, string1, source)
         call track_status(ens_size, zstatus, val, istatus, return_now)
         return
      endif
   enddo

! Clean up and return
   deallocate(co_val)
   deallocate(prs_iasi)
!
end subroutine get_expected_iasi_co_profile
!
subroutine set_obs_def_iasi_co_profile(key, co_avgker, co_press, co_prior, co_psurf, co_nlayers, co_nlevels)
!----------------------------------------------------------------------
! subroutine set_obs_def_iasi_co_profile(key, co_avgker, co_press, co_prior, co_psurf, co_nlayers, co_nlevels)

! Allows passing of obs_def special information 

integer,                 intent(in) :: key, co_nlayers, co_nlevels
real(r8), dimension(IASI_DIM),  intent(in) :: co_avgker
real(r8), dimension(IASI_DIMP), intent(in) :: co_press
real(r8),                intent(in) :: co_prior
real(r8),                intent(in) :: co_psurf

character(len=*), parameter :: routine = 'set_obs_def_iasi_co_profile'

if ( .not. module_initialized ) call initialize_module

if(num_iasi_co_obs >= MAX_IASI_CO_OBS) then
   write(string1, *)'Not enough space for a iasi CO obs.'
   write(string2, *)'Can only have MAX_IASI_CO_OBS (currently ',MAX_IASI_CO_OBS,')'
   call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2)
endif

avg_kernel(   key,1:co_nlayers)  = co_avgker(1:co_nlayers)
pressure(     key,1:co_nlevels) = co_press(1:co_nlevels)
iasi_prior(   key)   = co_prior
iasi_psurf(   key)   = co_psurf
iasi_nlayers(key)    = co_nlayers
iasi_nlevels(key)    = co_nlevels

end subroutine set_obs_def_iasi_co_profile
!
function read_iasi_prior(ifile, fform)
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform
real(r8)                               :: read_iasi_prior

character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_prior
   CASE DEFAULT
      read(ifile, *) read_iasi_prior
END SELECT
end function read_iasi_prior
!
function read_iasi_nlayers(ifile, fform)
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform
integer                                :: read_iasi_nlayers

character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_nlayers
   CASE DEFAULT
      read(ifile, *) read_iasi_nlayers
END SELECT
end function read_iasi_nlayers
!
function read_iasi_nlevels(ifile, fform)
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform
integer                                :: read_iasi_nlevels

character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_nlevels
   CASE DEFAULT
      read(ifile, *) read_iasi_nlevels
END SELECT
end function read_iasi_nlevels
!
subroutine write_iasi_prior(ifile, iasi_prior_temp, fform)
integer,          intent(in) :: ifile
real(r8),         intent(in) :: iasi_prior_temp
character(len=*), intent(in) :: fform

character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) iasi_prior_temp
   CASE DEFAULT
      write(ifile, *) iasi_prior_temp
END SELECT
end subroutine write_iasi_prior
!
subroutine write_iasi_nlayers(ifile, iasi_nlayers_temp, fform)
integer,          intent(in) :: ifile
integer,          intent(in) :: iasi_nlayers_temp
character(len=*), intent(in) :: fform

character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) iasi_nlayers_temp
   CASE DEFAULT
      write(ifile, *) iasi_nlayers_temp
END SELECT
end subroutine write_iasi_nlayers
!
subroutine write_iasi_nlevels(ifile, iasi_nlevels_temp, fform)
integer,          intent(in) :: ifile
integer,          intent(in) :: iasi_nlevels_temp
character(len=*), intent(in) :: fform

character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) iasi_nlevels_temp
   CASE DEFAULT
      write(ifile, *) iasi_nlevels_temp
END SELECT
end subroutine write_iasi_nlevels
!
function read_iasi_psurf(ifile, fform)
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform
real(r8)                               :: read_iasi_psurf

character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_psurf
   CASE DEFAULT
      read(ifile, *) read_iasi_psurf
END SELECT
end function read_iasi_psurf
!
subroutine write_iasi_psurf(ifile, iasi_psurf_temp, fform)
integer,          intent(in) :: ifile
real(r8),         intent(in) :: iasi_psurf_temp
character(len=*), intent(in) :: fform

character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) iasi_psurf_temp
   CASE DEFAULT
      write(ifile, *) iasi_psurf_temp
END SELECT
end subroutine write_iasi_psurf
!
function read_iasi_avg_kernels(ifile, nlayers, fform)
integer,          intent(in)           :: ifile, nlayers
character(len=*), intent(in), optional :: fform
real(r8), dimension(IASI_DIM)           :: read_iasi_avg_kernels

character(len=32)  :: fileformat
read_iasi_avg_kernels(:) = 0.0_r8
if ( .not. module_initialized ) call initialize_module
fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_avg_kernels(1:nlayers)
   CASE DEFAULT
      read(ifile, *) read_iasi_avg_kernels(1:nlayers)
END SELECT
end function read_iasi_avg_kernels
!
subroutine write_iasi_avg_kernels(ifile, avg_kernels_temp, nlayers_temp, fform)
integer,                 intent(in) :: ifile, nlayers_temp
real(r8), dimension(:), intent(in) :: avg_kernels_temp
character(len=*),        intent(in) :: fform

character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) avg_kernels_temp(1:nlayers_temp)
   CASE DEFAULT
      write(ifile, *) avg_kernels_temp(1:nlayers_temp)
END SELECT
end subroutine write_iasi_avg_kernels
!
function read_iasi_pressure(ifile, nlevels, fform)
integer,          intent(in)           :: ifile, nlevels
character(len=*), intent(in), optional :: fform
real(r8), dimension(IASI_DIMP)         :: read_iasi_pressure

character(len=32) :: fileformat
read_iasi_pressure(:) = 0.0_r8
if ( .not. module_initialized ) call initialize_module
fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_pressure(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_iasi_pressure(1:nlevels)
END SELECT
end function read_iasi_pressure
!
subroutine write_iasi_pressure(ifile, pressure_temp, nlevels_temp, fform)
integer,                 intent(in) :: ifile, nlevels_temp
real(r8), dimension(IASI_DIMP), intent(in)  :: pressure_temp
character(len=32),       intent(in) :: fform

character(len=32)  :: fileformat
if ( .not. module_initialized ) call initialize_module
fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) pressure_temp(1:nlevels_temp)
   CASE DEFAULT
      write(ifile, *) pressure_temp(1:nlevels_temp)
END SELECT
end subroutine write_iasi_pressure
!
!----------------------------------------------------------------------
!
end module obs_def_IASI_CO_PROFILE_mod
! END DART PREPROCESS MODULE CODE
