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

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! MOPITT_CO_TOTAL_COL, QTY_CO
! END DART PREPROCESS TYPE DEFINITIONS

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_mopitt_co_total_col_mod, only : write_mopitt_co_total_col, read_mopitt_co_total_col, &
!                                  interactive_mopitt_co_total_col, get_expected_mopitt_co_total_col, &
!                                  set_obs_def_mopitt_co_total_col
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(MOPITT_CO_TOTAL_COL)                                                           
!            call get_expected_mopitt_co_total_col(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(MOPITT_CO_TOTAL_COL)
!         call read_mopitt_co_total_col(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(MOPITT_CO_TOTAL_COL)
!         call write_mopitt_co_total_col(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(MOPITT_CO_TOTAL_COL)
!         call interactive_mopitt_co_total_col(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS SET_OBS_DEF_MOPITT_CO_TOTAL_COL
!      case(MOPITT_CO_TOTAL_COL)
!         call set_obs_def_mopitt_co_total_col(obs_def%key)
! END DART PREPROCESS SET_OBS_DEF_MOPITT_CO_TOTAL_COL

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_mopitt_co_total_col_mod

   use         apm_upper_bdy_mod, only :get_upper_bdy_fld, &
                                        get_MOZART_INT_DATA, &
                                        get_MOZART_REAL_DATA, &
                                        wrf_dart_ubval_interp, &
                                        apm_get_exo_coldens, &
                                        apm_get_upvals, &
                                        apm_interpolate

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

   public :: write_mopitt_co_total_col, &
             read_mopitt_co_total_col, &
             interactive_mopitt_co_total_col, &
             get_expected_mopitt_co_total_col, &
             set_obs_def_mopitt_co_total_col

! Storage for the special information required for observations of this type
   integer, parameter               :: MAX_MOPITT_CO_OBS = 10000000
   integer, parameter               :: MOPITT_DIM = 10
   integer                          :: num_mopitt_co_obs = 0
!
! MOPITT level pressures
   real(r8)   :: mopitt_prs(MOPITT_DIM) =(/ &
   90000.,80000.,70000.,60000.,50000.,40000.,30000.,20000.,10000.,5000. /)
!
! MOPITT layer pressures (layer 1 value is place holder for sfc_prs/90000 midpoint)
   real(r8)   :: mopitt_prs_mid(MOPITT_DIM) =(/ &
   100000.,85000.,75000.,65000.,55000.,45000.,35000.,25000.,15000.,7500. /)

   real(r8), allocatable, dimension(:,:) :: avg_kernel
   real(r8), allocatable, dimension(:)   :: mopitt_prior
   real(r8), allocatable, dimension(:)   :: mopitt_psurf
   integer,  allocatable, dimension(:)   :: mopitt_nlayers

! version controlled file description for error handling, do not edit
   character(len=*), parameter :: source   = 'obs_def_MOPITT_CO_TOTAL_COL_mod.f90'
   character(len=*), parameter :: revision = ''
   character(len=*), parameter :: revdate  = ''

   character(len=512) :: string1, string2
   character(len=200) :: upper_data_file
   character(len=200) :: upper_data_model
   character(len=200) :: model
   integer            :: ls_chem_dx, ls_chem_dy, ls_chem_dz, ls_chem_dt

   logical, save :: module_initialized = .false.
   integer  :: counts1 = 0

   character(len=129)  :: MOPITT_CO_retrieval_type
   logical             :: use_log_co
   integer             :: nlayer_model
   integer             :: nlayer_mopitt
   integer             :: nlayer_mopitt_co_total_col
   integer             :: nlayer_mopitt_co_trop_col
   integer             :: nlayer_mopitt_co_profile
!
! MOPITT_CO_retrieval_type:
!     RAWR - retrievals in format from supplier
!     RETR - retrievals in retrieval (ppbv) format
!     QOR  - quasi-optimal retrievals
!     CPSR - compact phase space retrievals

   namelist /obs_def_MOPITT_CO_nml/ upper_data_file, use_log_co, nlayer_model, &
   nlayer_mopitt_co_total_col, nlayer_mopitt_co_trop_col, nlayer_mopitt_co_profile, &
   ls_chem_dx, ls_chem_dy, ls_chem_dz, ls_chem_dt, upper_data_model

contains

!----------------------------------------------------------------------

subroutine initialize_module

   integer :: iunit, rc

! Prevent multiple calls from executing this code more than once.
   if (module_initialized) return

   call register_module(source, revision, revdate)
   module_initialized = .true.

   allocate (avg_kernel(    MAX_MOPITT_CO_OBS,MOPITT_DIM))
   allocate (mopitt_prior(  MAX_MOPITT_CO_OBS))
   allocate (mopitt_psurf(  MAX_MOPITT_CO_OBS))
   allocate (mopitt_nlayers(MAX_MOPITT_CO_OBS))

! Read the namelist entry.
   MOPITT_CO_retrieval_type='RETR'
   use_log_co=.false.
   call find_namelist_in_file("input.nml", "obs_def_MOPITT_CO_nml", iunit)
   read(iunit, nml = obs_def_MOPITT_CO_nml, iostat = rc)
   call check_namelist_read(iunit, rc, "obs_def_MOPITT_CO_nml")
   nlayer_mopitt=nlayer_mopitt_co_total_col

! Record the namelist values used for the run ... 
   if (do_nml_file()) write(nmlfileunit, nml=obs_def_MOPITT_CO_nml)
   if (do_nml_term()) write(     *     , nml=obs_def_MOPITT_CO_nml)

end subroutine initialize_module

subroutine read_mopitt_co_total_col(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine read_mopitt_co_total_col(key, ifile, fform)

integer,          intent(out)          :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=32)               :: fileformat

integer                         :: mopitt_nlayers_1
real(r8)                        :: mopitt_prior_1
real(r8)                        :: mopitt_psurf_1
real(r8), dimension(MOPITT_DIM) :: avg_kernels_1
integer                         :: keyin

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading
avg_kernels_1(:) = 0.0_r8
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   mopitt_nlayers_1 = read_mopitt_nlayers(ifile, fileformat)
   mopitt_prior_1 = read_mopitt_prior(ifile, fileformat)
   mopitt_psurf_1 = read_mopitt_psurf(ifile, fileformat)
   avg_kernels_1(1:mopitt_nlayers_1)  = read_mopitt_avg_kernels(ifile, mopitt_nlayers_1, fileformat)
   read(ifile) keyin
   CASE DEFAULT
   mopitt_nlayers_1 = read_mopitt_nlayers(ifile, fileformat)
   mopitt_prior_1 = read_mopitt_prior(ifile, fileformat)
   mopitt_psurf_1 = read_mopitt_psurf(ifile, fileformat)
   avg_kernels_1(1:mopitt_nlayers_1)  = read_mopitt_avg_kernels(ifile, mopitt_nlayers_1, fileformat)
   read(ifile, *) keyin
END SELECT

counts1 = counts1 + 1
key = counts1
call set_obs_def_mopitt_co_total_col(key, avg_kernels_1, mopitt_prior_1, &
     mopitt_psurf_1, mopitt_nlayers_1)

end subroutine read_mopitt_co_total_col

 subroutine write_mopitt_co_total_col(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine write_mopitt_co_total_col(key, ifile, fform)

integer,          intent(in)           :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat
real(r8), dimension(MOPITT_DIM) :: avg_kernels_temp

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading
   
avg_kernels_temp=avg_kernel(key,:)

SELECT CASE (fileformat)
   
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   call write_mopitt_nlayers(ifile, mopitt_nlayers(key), fileformat)
   call write_mopitt_prior(ifile, mopitt_prior(key), fileformat)
   call write_mopitt_psurf(ifile, mopitt_psurf(key), fileformat)
   call write_mopitt_avg_kernels(ifile, avg_kernels_temp, mopitt_nlayers(key), fileformat)
   write(ifile) key

   CASE DEFAULT
   call write_mopitt_nlayers(ifile, mopitt_nlayers(key), fileformat)
   call write_mopitt_prior(ifile, mopitt_prior(key), fileformat)
   call write_mopitt_psurf(ifile, mopitt_psurf(key), fileformat)
   call write_mopitt_avg_kernels(ifile, avg_kernels_temp, mopitt_nlayers(key), fileformat)
   write(ifile, *) key
END SELECT 
end subroutine write_mopitt_co_total_col
!
subroutine interactive_mopitt_co_total_col(key)
!----------------------------------------------------------------------
!subroutine interactive_mopitt_co_total_col(key)
!
! Initializes the specialized part of a MOPITT observation
! Passes back up the key for this one

integer, intent(out) :: key

if ( .not. module_initialized ) call initialize_module

! Make sure there's enough space, if not die for now (clean later)
if(num_mopitt_co_obs >= MAX_MOPITT_CO_OBS) then
   write(string1, *)'Not enough space for a mopitt CO obs.'
   write(string2, *)'Can only have MAX_MOPITT_CO_OBS (currently ',MAX_MOPITT_CO_OBS,')'
   call error_handler(E_ERR,'interactive_mopitt_co',string1,source,revision,revdate, text2=string2)
endif

! Increment the index
num_mopitt_co_obs = num_mopitt_co_obs + 1
key = num_mopitt_co_obs

! Otherwise, prompt for input for the three required beasts
write(*, *) 'Creating an interactive_mopitt_co observation'
write(*, *) 'Input the MOPITT Prior '
read(*, *) mopitt_prior
write(*, *) 'Input MOPITT Surface Pressure '
read(*, *) mopitt_psurf(num_mopitt_co_obs)
write(*, *) 'Input the 10 Averaging Kernel Weights '
read(*, *) avg_kernel(num_mopitt_co_obs,:)
end subroutine interactive_mopitt_co_total_col
!
subroutine get_expected_mopitt_co_total_col(state_handle, ens_size, location, key, val, istatus)
!----------------------------------------------------------------------
!subroutine get_expected_mopitt_co_total_col(state_handle, ens_size, location, key, val, istatus)
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
   integer             :: layer_mopitt,level_mopitt
   integer             :: layer_mdl,level_mdl

   real(r8)            :: prs_mopitt_sfc, temp_sfc
   real(r8)            :: co_min, co_min_log, level
   real(r8)            :: mloc(3)
   real(r8),dimension(ens_size) :: prs_mdl_sfc
   real(r8),dimension(ens_size) :: prs_mdl_1, co_mdl_1
   real(r8),dimension(ens_size) :: prs_mdl_n, co_mdl_n
   real(r8),dimension(:),allocatable   :: prs_mopitt
   real(r8),dimension(:,:),allocatable :: co_val

   character(len=*), parameter :: routine = 'get_expected_mopitt_co_total_col'
   logical :: return_now
!
! Initialize DART
   
   istatus(:) = 0
   if ( .not. module_initialized ) call initialize_module
!
! Initialize variables (MOPITT is ppbv; WRF CO is ppmv)
   co_min      = 1.e-2
   co_min_log  = log(co_min)
   kstr        = nlayer_mopitt-mopitt_nlayers(key)+1
   val(:)      = 0.0
   layer_mopitt = mopitt_nlayers(key)
   level_mopitt = layer_mopitt+1
   layer_mdl = nlayer_model
   level_mdl = layer_mdl+1
   
   allocate(prs_mopitt(layer_mopitt))
   prs_mopitt(1:layer_mopitt)=mopitt_prs_mid(kstr:nlayer_mopitt)
   prs_mopitt(1)=(mopitt_psurf(key)+mopitt_prs(kstr))/2.

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

! Get profiles at MOPITT pressure levels
   do imem=1,ens_size   
      if (prs_mdl_sfc(imem).gt.mopitt_psurf(key)) then
         temp_sfc=(prs_mdl_sfc(imem)+mopitt_prs(kstr))/2.
         if (temp_sfc.gt.prs_mopitt(1)) then
            prs_mopitt(1)=temp_sfc
         endif
      endif
   enddo
!
   allocate(co_val(ens_size,layer_mopitt))
   do k=1,layer_mopitt
      zstatus=0
      loc2 = set_location(mloc(1), mloc(2), prs_mopitt(k), VERTISPRESSURE)
      call interpolate(state_handle, ens_size, loc2, QTY_CO, co_val(:,k), zstatus)  
!
! Correcting for expected failures near the surface
      do imem=1,ens_size
         if (prs_mopitt(k).ge.prs_mdl_1(imem)) then
            co_val(imem,k) = co_mdl_1(imem)
            cycle
         endif

! Correcting for expected failures near the top
         if (prs_mopitt(k).le.prs_mdl_n(imem)) then
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
! Convert units for co to ppbv for MOPITT
         if(use_log_co) then
            co_val(imem,k) = log(exp(co_val(imem,k)) * 1.e3_r8)
         else
            co_val(imem,k) = co_val(imem,k) * 1.e3_r8
         endif
      enddo
   enddo
!
! Apply MOPITT Averaging kernel A and MOPITT Prior (I-A)xa
! x = Axm + (I-A)xa , where x is a 10 element vector 
!
! loop through MOPITT layers
   val(:) = 0.0_r8
   do imem=1,ens_size
      do k = 1, layer_mopitt
         if(use_log_co) then
            val(imem) = val(imem) + avg_kernel(key,k) * log10(exp(co_val(imem,k)))  
         else
            val(imem) = val(imem) + avg_kernel(key,k) * log10(co_val(imem,k))  
         endif
      enddo
!
! NOTE: For the following the mopitt_prior is zero due to the QOR subtraction
      if (trim(MOPITT_CO_retrieval_type).eq.'RAWR' .or. trim(MOPITT_CO_retrieval_type).eq.'QOR' &
      .or. trim(MOPITT_CO_retrieval_type).eq.'CPSR') then
         val(imem) = val(imem) + mopitt_prior(key)
      elseif (trim(MOPITT_CO_retrieval_type).eq.'RETR') then
         val(imem) = val(imem) + mopitt_prior(key)
         val(imem) = (10.**val(imem)) * 1.e-3
      endif

      if(val(imem).lt.0) then
         zstatus(imem)=20
         val(:)=missing_r8
         write(string1, *) &
         'APM NOTICE: MOPITT CO expected value is negative '
         call error_handler(E_MSG, routine, string1, source)
         call track_status(ens_size, zstatus, val, istatus, return_now)
         return
      endif
   enddo

! Clean up and return
   deallocate(co_val)
   deallocate(prs_mopitt)
!
end subroutine get_expected_mopitt_co_total_col
!
!----------------------------------------------------------------------

 subroutine set_obs_def_mopitt_co_total_col(key, co_avgker, co_prior, co_psurf, co_nlayers)
!----------------------------------------------------------------------
! Allows passing of obs_def special information 

integer,                 intent(in) :: key, co_nlayers
real(r8), dimension(10), intent(in) :: co_avgker
real(r8),                intent(in) :: co_prior
real(r8),                intent(in) :: co_psurf

if ( .not. module_initialized ) call initialize_module

if(num_mopitt_co_obs >= MAX_MOPITT_CO_OBS) then
   write(string1, *)'Not enough space for a mopitt CO obs.'
   call error_handler(E_MSG,'set_obs_def_mopitt_co',string1,source,revision,revdate)
   write(string1, *)'Can only have MAX_MOPITT_CO_OBS (currently ',MAX_MOPITT_CO_OBS,')'
   call error_handler(E_ERR,'set_obs_def_mopitt_co',string1,source,revision,revdate)
endif

avg_kernel(key,:)   = co_avgker(:)
mopitt_prior(key)   = co_prior
mopitt_psurf(key)   = co_psurf
mopitt_nlayers(key) = co_nlayers

end subroutine set_obs_def_mopitt_co_total_col

function read_mopitt_prior(ifile, fform)

integer,                    intent(in) :: ifile
real(r8)                               :: read_mopitt_prior
character(len=*), intent(in), optional :: fform

character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_mopitt_prior
   CASE DEFAULT
      read(ifile, *) read_mopitt_prior
END SELECT

end function read_mopitt_prior

function read_mopitt_nlayers(ifile, fform)

integer,                    intent(in) :: ifile
integer                               :: read_mopitt_nlayers
character(len=*), intent(in), optional :: fform

character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_mopitt_nlayers
   CASE DEFAULT
      read(ifile, *) read_mopitt_nlayers
END SELECT

end function read_mopitt_nlayers


subroutine write_mopitt_prior(ifile, mopitt_prior_temp, fform)

integer,           intent(in) :: ifile
real(r8),          intent(in) :: mopitt_prior_temp
character(len=32), intent(in) :: fform

character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) mopitt_prior_temp
   CASE DEFAULT
      write(ifile, *) mopitt_prior_temp
END SELECT

end subroutine write_mopitt_prior

subroutine write_mopitt_nlayers(ifile, mopitt_nlayers_temp, fform)

integer,                    intent(in) :: ifile
integer,                    intent(in) :: mopitt_nlayers_temp
character(len=32),          intent(in) :: fform

character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) mopitt_nlayers_temp
   CASE DEFAULT
      write(ifile, *) mopitt_nlayers_temp
END SELECT

end subroutine write_mopitt_nlayers


function read_mopitt_psurf(ifile, fform)

integer,                    intent(in) :: ifile
real(r8)                               :: read_mopitt_psurf
character(len=*), intent(in), optional :: fform

character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_mopitt_psurf
   CASE DEFAULT
      read(ifile, *) read_mopitt_psurf
END SELECT

end function read_mopitt_psurf

subroutine write_mopitt_psurf(ifile, mopitt_psurf_temp, fform)

integer,           intent(in) :: ifile
real(r8),          intent(in) :: mopitt_psurf_temp
character(len=32), intent(in) :: fform

character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) mopitt_psurf_temp
   CASE DEFAULT
      write(ifile, *) mopitt_psurf_temp
END SELECT

end subroutine write_mopitt_psurf

function read_mopitt_avg_kernels(ifile, nlayers, fform)

integer,                    intent(in) :: ifile, nlayers
real(r8), dimension(10)        :: read_mopitt_avg_kernels
character(len=*), intent(in), optional :: fform

character(len=32)  :: fileformat

read_mopitt_avg_kernels(:) = 0.0_r8

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_mopitt_avg_kernels(1:nlayers)
   CASE DEFAULT
      read(ifile, *) read_mopitt_avg_kernels(1:nlayers)
END SELECT

end function read_mopitt_avg_kernels

subroutine write_mopitt_avg_kernels(ifile, avg_kernels_temp, nlayers_temp, fform)

integer,                    intent(in) :: ifile, nlayers_temp
real(r8), dimension(10), intent(in)  :: avg_kernels_temp
character(len=32),          intent(in) :: fform

character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) avg_kernels_temp(1:nlayers_temp)
   CASE DEFAULT
      write(ifile, *) avg_kernels_temp(1:nlayers_temp)
END SELECT

end subroutine write_mopitt_avg_kernels

end module obs_def_mopitt_co_total_col_mod
! END DART PREPROCESS MODULE CODE

