! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: obs_def_CO_Nadir_mod.f90 11289 2017-03-10 21:56:06Z hendric@ucar.edu $
!!! observation operator for CO profiles
!!! March 2021


! BEGIN DART PREPROCESS KIND LIST
! MOPITT_CO_RETRIEVAL, QTY_CO
! MOPITT_CO_RETRIEVAL_TIR, QTY_CO
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_profile_mopitt_mod, only : write_mopitt_co, read_mopitt_co, interactive_mopitt_co, &
!                                          get_expected_mopitt_co, set_obs_def_mopitt_co
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(MOPITT_CO_RETRIEVAL, MOPITT_CO_RETRIEVAL_TIR)                                                           
!            call get_expected_mopitt_co(state_handle, ens_size, location, obs_def%key, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(MOPITT_CO_RETRIEVAL, MOPITT_CO_RETRIEVAL_TIR)
!         call read_mopitt_co(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(MOPITT_CO_RETRIEVAL, MOPITT_CO_RETRIEVAL_TIR)
!         call write_mopitt_co(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(MOPITT_CO_RETRIEVAL, MOPITT_CO_RETRIEVAL_TIR)
!         call interactive_mopitt_co(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS SET_OBS_DEF_MOPITT_CO
!      case(MOPITT_CO_RETRIEVAL, MOPITT_CO_RETRIEVAL_TIR)
!         call set_obs_def_mopitt_co(obs_def%key)
! END DART PREPROCESS SET_OBS_DEF_MOPITT_CO


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_profile_mopitt_mod

use typeSizes
use        types_mod, only : i8, r8, MISSING_R8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, set_location, get_location, query_location, &
                             VERTISPRESSURE, VERTISLEVEL, VERTISSURFACE 

use  assim_model_mod, only : interpolate, get_state_meta_data
use    obs_kind_mod, only  : QTY_CO, QTY_PRESSURE, QTY_SURFACE_PRESSURE
use ensemble_manager_mod,  only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none

public :: write_mopitt_co, read_mopitt_co, interactive_mopitt_co, &
          get_expected_mopitt_co, set_obs_def_mopitt_co

! Storage for the special information required for observations of this type
integer, parameter               :: max_mopitt_co_obs = 10000000
integer, parameter               :: mopitt_dim = 10
integer, parameter               :: max_model_levs = 33   
integer                          :: num_mopitt_co_obs = 0
! KDR replace 10 with mopitt_dim?
! No it is specific to the 10th level
! real(r8), dimension(max_mopitt_co_obs,10) :: avg_kernel
real(r8), dimension(max_mopitt_co_obs,mopitt_dim) :: avg_kernel
real(r8), dimension(max_mopitt_co_obs)            :: mopitt_prior
! Hardcoded pressure levels for MOPITT CO obs, it shouldn't be...
real(r8)   :: mopitt_pressure(mopitt_dim) =(/ &
                              95000.,90000.,80000.,70000.,60000.,50000.,40000.,30000.,20000.,10000. /)
real(r8), dimension(max_mopitt_co_obs)   :: mopitt_psurf
integer,  dimension(max_mopitt_co_obs)   :: mopitt_nlevels

! For now, read in all info on first read call, write all info on first write call
logical :: already_read = .false., already_written = .false.

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/observations/forward_operators/obs_def_CO_Nadir_mod.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 11289 $"
character(len=128), parameter :: revdate  = "$Date: 2017-03-10 14:56:06 -0700 (Fri, 10 Mar 2017) $"


logical, save :: module_initialized = .false.
integer  :: counts1 = 0

contains

!----------------------------------------------------------------------

subroutine initialize_module
!----------------------------------------------------------------------------
! subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module



subroutine read_mopitt_co(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine read_mopitt_co(key, ifile, fform)

integer, intent(out)            :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional    :: fform
character(len=32)               :: fileformat
integer                         :: mopitt_nlevels_1
real(r8)                        :: mopitt_prior_1
real(r8)                        :: mopitt_psurf_1
real(r8), dimension(mopitt_dim) :: avg_kernels_1
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
   mopitt_nlevels_1 = read_mopitt_nlevels(ifile, fileformat)
   mopitt_prior_1 = read_mopitt_prior(ifile, fileformat)
   mopitt_psurf_1 = read_mopitt_psurf(ifile, fileformat)
   avg_kernels_1(1:mopitt_nlevels_1) = read_mopitt_avg_kernels(ifile, mopitt_nlevels_1, fileformat)
   read(ifile) keyin

   CASE DEFAULT
   mopitt_nlevels_1 = read_mopitt_nlevels(ifile, fileformat)
   mopitt_prior_1 = read_mopitt_prior(ifile, fileformat)
   mopitt_psurf_1 = read_mopitt_psurf(ifile, fileformat)
   avg_kernels_1(1:mopitt_nlevels_1)  = read_mopitt_avg_kernels(ifile, mopitt_nlevels_1, fileformat)
   read(ifile, *) keyin
END SELECT

counts1 = counts1 + 1
key = counts1
call set_obs_def_mopitt_co(key, avg_kernels_1, mopitt_prior_1, mopitt_psurf_1, &
                           mopitt_nlevels_1)

end subroutine read_mopitt_co

subroutine write_mopitt_co(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine write_mopitt_co(key, ifile, fform)

integer, intent(in)             :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional :: fform

character(len=32)               :: fileformat
real(r8), dimension(mopitt_dim) :: avg_kernels_temp

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading
   
avg_kernels_temp=avg_kernel(key,:)

SELECT CASE (fileformat)
   
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   call write_mopitt_nlevels(ifile, mopitt_nlevels(key), fileformat)
   call write_mopitt_prior(ifile, mopitt_prior(key), fileformat)
   call write_mopitt_psurf(ifile, mopitt_psurf(key), fileformat)
   call write_mopitt_avg_kernels(ifile, avg_kernels_temp, mopitt_nlevels(key), fileformat)
   write(ifile) key

   CASE DEFAULT
   call write_mopitt_nlevels(ifile, mopitt_nlevels(key), fileformat)
   call write_mopitt_prior(ifile, mopitt_prior(key), fileformat)
   call write_mopitt_psurf(ifile, mopitt_psurf(key), fileformat)
   call write_mopitt_avg_kernels(ifile, avg_kernels_temp, mopitt_nlevels(key), fileformat)
   write(ifile, *) key
END SELECT 

end subroutine write_mopitt_co


subroutine interactive_mopitt_co(key)
!----------------------------------------------------------------------
!subroutine interactive_mopitt_co(key)
!
! Initializes the specialized part of a MOPITT observation
! Passes back up the key for this one

integer, intent(out) :: key

character(len=129) :: msgstring

if ( .not. module_initialized ) call initialize_module

! Make sure there's enough space, if not die for now (clean later)
if(num_mopitt_co_obs >= max_mopitt_co_obs) then
   ! PUT IN ERROR HANDLER CALL
   write(msgstring, *)'Not enough space for a mopitt CO obs.'
   call error_handler(E_MSG,'interactive_mopitt_co',msgstring,source,revision,revdate)
   write(msgstring, *)'Can only have max_mopitt_co_obs (currently ',max_mopitt_co_obs,')'
   call error_handler(E_ERR,'interactive_mopitt_co',msgstring,source,revision,revdate)
endif

! Increment the index
num_mopitt_co_obs = num_mopitt_co_obs + 1
key = num_mopitt_co_obs

! Otherwise, prompt for input for the three required beasts
write(*, *) 'Creating an interactive_mopitt_co observation'
write(*, *) 'Input the MOPITT Prior '
read(*, *) mopitt_prior(key)
write(*, *) 'Input MOPITT Surface Pressure '
read(*, *) mopitt_psurf(key)
write(*, *) 'Input the 10 Averaging Kernel Weights '
read(*, *) avg_kernel(key,:)

end subroutine interactive_mopitt_co

subroutine get_expected_mopitt_co(state_handle, ens_size, location, key, val, istatus)
!----------------------------------------------------------------------
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: key
real(r8),            intent(out) :: val(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer :: i,j
type(location_type)     :: loc0,loc1,loc2,loc3,loc3p,loc3m,locS
real(r8)                :: mloc(3), mloc1(3), mloc2(3)
real(r8)                :: obs_val(ens_size), obs_val_int(ens_size)

integer                 :: nlevels, start_i, end_i
integer                 :: iens ! BG

real(r8)                :: top_pres(ens_size), bot_pres(ens_size), coef(ens_size), mop_layer_wght(ens_size)
real(r8)                :: i_top_pres(ens_size), i_bot_pres(ens_size), i_pres
integer                 :: lev, n_model_lev, n_obs_lev
real(r8)                :: p0(ens_size), p1(ens_size)
real(r8), allocatable   :: p_col(:, :), co_col(:, :)
real(r8), allocatable   :: mopitt_pres_local(:, :)
integer                 :: p_col_istatus(ens_size), co_col_istatus(ens_size), obs_val_int_istatus(ens_size)
logical                 :: return_now, islog
integer                 :: imem

if ( .not. module_initialized ) call initialize_module
mloc = get_location(location)

val = 0.0_r8

! prelude - Check if position is valid
if (mloc(2)>90.0_r8) then
    mloc(2)=90.0_r8
elseif (mloc(2)<-90.0_r8) then
    mloc(2)=-90.0_r8
endif

! part 1 - get model quatities
! model profile of pressure at interfaces
! model profile of concentration

! "silly" way to get model levels but no other way to do this...
istatus = 0
p_col = MISSING_R8
lev = 1
model_levels: do
   locS = set_location(mloc(1),mloc(2),real(lev,r8),VERTISLEVEL)
   call interpolate(state_handle, ens_size, loc0, QTY_PRESSURE, p0, istatus)
   if (any(istatus /= 2)) then
      if (any(istatus /= 0)) then
         n_model_lev = lev
         exit model_levels
      endif
   endif
   lev = lev + 1
enddo model_levels
write(*,*) 'number of levels'
write(*,*) n_model_lev

allocate(p_col(ens_size, n_model_lev+1), co_col(ens_size, n_model_lev))
! get profile concentration and profile pressure at interfaces (approximate).
do i = 1, n_model_lev

   ! get locations
   loc0 = set_location(mloc(1),mloc(2),real(i,r8),VERTISLEVEL)
   loc1 = set_location(mloc(1),mloc(2),real(i+1,r8),VERTISLEVEL)

   ! get pressure at interfaces
   call interpolate(state_handle, ens_size, loc0, QTY_PRESSURE, p0, istatus)
   call interpolate(state_handle, ens_size, loc1, QTY_PRESSURE, p1, istatus)
   if (i == 1) then
      p_col(:, i) = 0.5_r8 * p0 ! approximation since there not such capability in model mods to get pressure at INTERFACES
   else
      p_col(:, i) = 0.5_r8 * (p0 + p1) ! approximation since there not such capability in model mods to get pressure at INTERFACES
   endif

   !get concentration profile at grid cell center
   call interpolate(state_handle, ens_size, loc0, QTY_CO, co_col(:, lev), istatus)
enddo

!get surface pressure at the bottom
call interpolate(state_handle, ens_size, locS, QTY_SURFACE_PRESSURE, p_col(:, n_model_lev+1), istatus)

locS = set_location(mloc(1),mloc(2),0.0_r8, VERTISSURFACE)

call track_status(ens_size, p_col_istatus, p_col(:,n_model_lev+1), istatus, return_now)
if (return_now) return

! part 2 - get retrieval info
! avg_kernel
! mopitt_prior
! mopitt_pressure at interfaces

! to be changed later once the obs seq are updated: convention is that we want to work with increasing pressure coordinates as in model
! therefore we indexing like (n:1:-1)
n_obs_lev = mopitt_nlevels(key)
avg_kernel(key,:) = avg_kernel(key,n_obs_lev:1:-1)

allocate(mopitt_pres_local(ens_size, n_obs_lev+1))
do imem = 1, ens_size
   mopitt_pres_local(imem, 2:n_obs_lev+1) = mopitt_pressure(n_obs_lev:1:-1)
   mopitt_pres_local(imem, n_obs_lev+1) = p_col(imem,n_model_lev+1) ! here we want to make sure havew have the same surface pressure as in model
   mopitt_pres_local(imem, 1) = p_col(imem, 1) ! mopitt top interface pressure is TOA but we want to have the same top interface pressure as model top
enddo

! part 3 - call the column operator function and then deallocate

do imem = 1, ens_size
   call simulate_column_ob(n_obs_lev, n_model_lev, avg_kernel, &
   mopitt_pres_local(imem, :), p_col(imem, :), co_col(imem, :), obs_val(imem), islog=.true.)
enddo

deallocate(p_col, co_col, mopitt_pres_local)

end subroutine get_expected_mopitt_co

!!!! BEGIN: THIS SOULD GO IN A SEPARATE MODULE

subroutine simulate_column_ob(nlayers_obs, nlayers_model, avgkernel_obs, &
                              prsi_obs, prsi_model, profile_model, hofx, islog)

integer, intent(in   ) :: nlayers_obs, nlayers_model
real(r8), intent(in   ), dimension(nlayers_obs) :: avgkernel_obs
real(r8), intent(in   ), dimension(nlayers_obs+1) :: prsi_obs
real(r8), intent(in   ), dimension(nlayers_model+1) :: prsi_model
real(r8), intent(in   ), dimension(nlayers_model) :: profile_model
real(r8), intent(  out) :: hofx
real(r8) :: wf_a, wf_b, avgkernel, grav, M_dryair
real(r8), dimension(nlayers_obs) :: profile_obslayers
real(r8), dimension(nlayers_obs+1) :: pobs
real(r8), dimension(nlayers_model+1) :: pmod
integer, parameter :: max_string=800
character(len=max_string) :: err_msg
integer :: k, j, wi_a, wi_b
logical :: islog

hofx = 0.0_r8
profile_obslayers = 0.0_r8
grav = 9.80665_r8
M_dryair = 0.0289645_r8 ! kg/mol, dry air molar mass
avgkernel = 1.0_r8
do k=1,nlayers_obs
   ! get obs layer bound model indexes and weights for staggered
   ! obs and geoval levels
   call vert_interp_weights(nlayers_model+1, pobs(k), pmod, wi_a, wf_a)
   call vert_interp_weights(nlayers_model+1, pobs(k+1), pmod, wi_b, wf_b)

   !check if pmod is monotonic and decreasing
   if ((pmod(wi_a+1) < pmod(wi_a)) .or. (pmod(wi_b+1) < pmod(wi_b))) then
     write(*, *) "Error: inverted pressure coordinate in geovals, &
             &convention: top->bottom, decreasing pressures"
     call abor1_ftn(err_msg)
   end if

   ! when multiple mopdel levels are in a obs layer
   if ( wi_a < wi_b ) then
      profile_obslayers(k) = profile_obslayers(k) + profile_model(wi_a) * &
           (pmod(wi_a+1)-pmod(wi_a)) * wf_a / (M_dryair*grav)
      do j=wi_a+1,wi_b-1
         profile_obslayers(k) = profile_obslayers(k) + profile_model(j) * &
              (pmod(j+1)-pmod(j)) / (M_dryair*grav)
      enddo
      profile_obslayers(k) = profile_obslayers(k) + profile_model(wi_b) * &
           (pmod(wi_b+1)-pmod(wi_b)) * (1.0_r8-wf_b) / (M_dryair*grav)

   ! when multiple obs layers are in a model level
   else if ( wi_a == wi_b ) then
      profile_obslayers(k) = profile_obslayers(k) + profile_model(wi_a) * &
           (pmod(wi_a+1)-pmod(wi_a)) * (wf_a-wf_b) / (M_dryair*grav)

   ! if pressures coordinates are inverted return exception
   else if ( wi_a > wi_b ) then
      write(*, *) "Error: inverted pressure coordinate in obs, &
              &convention: top->bottom, decreasing pressures"
      call abor1_ftn(err_msg)
   end if

   if (islog) then
      hofx = hofx + (avgkernel * log10(profile_obslayers(k)))
   else if (.not. islog) then
      ! if not log, then just multiply by the avg kernel
      hofx = hofx + (avgkernel * profile_obslayers(k))
   else
      write(*, *) "Error: islog must be .true. or .false."
      call abor1_ftn(err_msg)
   end if
end do

if (islog) then
   hofx = 10.0**hofx
endif

end subroutine simulate_column_ob

subroutine vert_interp_weights(nlev,obl,vec,wi,wf)

   implicit none
   integer,         intent(in ) :: nlev       !Number of model levels
   real(r8), intent(in ) :: obl        !Observation location
   real(r8), intent(in ) :: vec(nlev)  !Structured vector of grid points
   integer,         intent(out) :: wi         !Index for interpolation
   real(r8), intent(out) :: wf         !Weight for interpolation
   
   integer         :: k
   
   if (vec(1) < vec(nlev)) then !Pressure increases with index
     if (obl < vec(1)) then
        wi = 1
        wf = 1.0
     elseif (obl > vec(nlev)) then
        wi = nlev - 1
        wf = 0.0
     else
        do k = 1,nlev-1
           if (obl >= vec(k) .and. obl <= vec(k+1)) then
              wi = k
           endif
        enddo
        wf = (vec(wi+1) - obl)/(vec(wi+1) - vec(wi))
     endif
   else !Pressure decreases with index
     if (obl > vec(1)) then
        wi = 1
        wf = 1.0
     elseif (obl < vec(nlev)) then
        wi = nlev - 1
        wf = 0.0
     else
        do k = 1,nlev-1
           if (obl >= vec(k+1) .and. obl <= vec(k)) then
              wi = k
           endif
        enddo
        wf = (vec(wi+1) - obl)/(vec(wi+1) - vec(wi))
     endif
   endif
   
end subroutine vert_interp_weights

!!!! THIS SOULD GO IN A SEPARATE MODULE: END

subroutine set_obs_def_mopitt_co(key, co_avgker, co_prior, co_psurf, co_nlevels)
!----------------------------------------------------------------------
! Allows passing of obs_def special information 

integer,                        intent(in) :: key, co_nlevels
real(r8),dimension(mopitt_dim), intent(in) :: co_avgker
real(r8),                       intent(in) :: co_prior
real(r8),                       intent(in) :: co_psurf
character(len=129)                         :: msgstring

if ( .not. module_initialized ) call initialize_module

if(num_mopitt_co_obs >= max_mopitt_co_obs) then
   ! PUT IN ERROR HANDLER CALL
   write(msgstring, *)'Not enough space for a mopitt CO obs.'
   call error_handler(E_MSG,'set_obs_def_mopitt_co',msgstring,source,revision,revdate)
   write(msgstring, *)'Can only have max_mopitt_co_obs (currently ',max_mopitt_co_obs,')'
   call error_handler(E_ERR,'set_obs_def_mopitt_co',msgstring,source,revision,revdate)
endif

avg_kernel(key,:)       = co_avgker(:)
mopitt_prior(key)       = co_prior
mopitt_psurf(key)       = co_psurf
mopitt_nlevels(key)     = co_nlevels

end subroutine set_obs_def_mopitt_co

function read_mopitt_prior(ifile, fform)

integer,                    intent(in) :: ifile
real(r8)                               :: read_mopitt_prior
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
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

function read_mopitt_nlevels(ifile, fform)

integer,                    intent(in) :: ifile
integer                               :: read_mopitt_nlevels
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_mopitt_nlevels
   CASE DEFAULT
      read(ifile, *) read_mopitt_nlevels
END SELECT

end function read_mopitt_nlevels


subroutine write_mopitt_prior(ifile, mopitt_prior_temp, fform)

integer,                    intent(in) :: ifile
real(r8),                   intent(in) :: mopitt_prior_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
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

subroutine write_mopitt_nlevels(ifile, mopitt_nlevels_temp, fform)

integer,                    intent(in) :: ifile
integer,                    intent(in) :: mopitt_nlevels_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) mopitt_nlevels_temp
   CASE DEFAULT
      write(ifile, *) mopitt_nlevels_temp
END SELECT

end subroutine write_mopitt_nlevels


function read_mopitt_psurf(ifile, fform)

integer,                    intent(in) :: ifile
real(r8)                               :: read_mopitt_psurf
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
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

integer,                    intent(in) :: ifile
real(r8),                   intent(in) :: mopitt_psurf_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
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

function read_mopitt_avg_kernels(ifile, nlevels, fform)

integer,                    intent(in) :: ifile, nlevels
real(r8), dimension(mopitt_dim)        :: read_mopitt_avg_kernels
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

read_mopitt_avg_kernels(:) = 0.0_r8

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_mopitt_avg_kernels(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_mopitt_avg_kernels(1:nlevels)
END SELECT

end function read_mopitt_avg_kernels

subroutine write_mopitt_avg_kernels(ifile, avg_kernels_temp, nlevels_temp, fform)

integer,                    intent(in) :: ifile, nlevels_temp
real(r8), dimension(mopitt_dim), intent(in)  :: avg_kernels_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) avg_kernels_temp(1:nlevels_temp)
   CASE DEFAULT
      write(ifile, *) avg_kernels_temp(1:nlevels_temp)
END SELECT

end subroutine write_mopitt_avg_kernels


end module obs_def_profile_mopitt_mod
! END DART PREPROCESS MODULE CODE

! <next few lines under version control, do not edit>
! $URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/observations/forward_operators/obs_def_CO_Nadir_mod.f90 $
! $Id: obs_def_CO_Nadir_mod.f90 11289 2017-03-10 21:56:06Z hendric@ucar.edu $
! $Revision: 11289 $
! $Date: 2017-03-10 14:56:06 -0700 (Fri, 10 Mar 2017) $
