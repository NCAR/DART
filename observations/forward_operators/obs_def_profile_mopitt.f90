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
use        types_mod, only : r8, MISSING_R8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, set_location, get_location, VERTISPRESSURE, &
                             VERTISLEVEL, VERTISSURFACE

use  assim_model_mod, only : interpolate
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
type(location_type) :: loc1,loc2,loc3,loc3p,loc3m,locS
real(r8)            :: mloc(3), mloc1(3), mloc2(3)
real(r8)            :: obs_val(ens_size), obs_val_int(ens_size)

integer             :: nlevels, start_i, end_i
integer             :: iens ! BG

real(r8)            :: top_pres(ens_size), bot_pres(ens_size), coef(ens_size), mop_layer_wght(ens_size)
real(r8)            :: i_top_pres(ens_size), i_bot_pres(ens_size), i_pres
integer             :: num_levs, lev
real(r8)            :: p_col(ens_size, max_model_levs)
real(r8)            :: mopitt_pres_local(ens_size, mopitt_dim)
integer,  allocatable :: dim_sizes(:)
integer             :: p_col_istatus(ens_size), obs_val_int_istatus(ens_size)
logical             :: return_now
integer             :: imem

if ( .not. module_initialized ) call initialize_module
mloc = get_location(location)

! Apply MOPITT Column Averaging kernel vector a_t and MOPITT Prior xa- a_t xa
! x = a_t xm + xa- a_t xa , where x is a 10 element vector 
! xm is the 10-element partial columnis at MOPITT grid
! also here, the averaging kernel is operated in log VMR.
! i.e. the state vector CO is in log VMR 
! edit dart_to_cam and cam_to_dart for the transformation
! also, it turns our xa - a_txa is a scalar quantity  
! KDR Why not initialize to mopitt_prior(key), instead of adding it after 
!     the end of the i-loop?
val = 0.0_r8

if (mloc(2)>90.0_r8) then
    mloc(2)=90.0_r8
elseif (mloc(2)<-90.0_r8) then
    mloc(2)=-90.0_r8
endif

do imem = 1, ens_size
   mopitt_pres_local(imem, :) = mopitt_pressure
enddo

nlevels = mopitt_nlevels(key)
! Modify AFAJ (072513)
! KDR start_i >= 1
start_i = mopitt_dim-nlevels+1
! KDR Ack! redefining one element of a global array.
!     Various elements will be redefined during various calls to this subroutine,
!     because the array is initialized in the specification statement.
mopitt_pres_local(:, start_i)=mopitt_psurf(key)
end_i = mopitt_dim

!write(*,*) 'BG in get_expected_mopitt_co'

! Find the number of model levels and the pressures on them.
istatus = 0
p_col = MISSING_R8
lev = 1
model_levels: do
   locS = set_location(mloc(1),mloc(2),real(lev,r8),VERTISLEVEL)

   !write(*,*) 'call interpolate'
   call interpolate(state_handle, ens_size, locS, QTY_PRESSURE, p_col(:, lev), p_col_istatus)
   !write(*,*) 'After call interpolate'

   !write(*,*) p_col_istatus
   !write(*,*) p_col(:, lev)
   if (any(p_col_istatus /= 2)) then
      if (any(p_col_istatus /= 0)) then
         !write(*,*) 'bb'
         p_col(:, lev) = MISSING_R8
         num_levs = lev - 1
         !write(*,*) 'BG in loop loop BG EXITING'
         !write(*,*) num_levs
         !write(*,*) 'EXIT loop loop BG EXITING'
         exit model_levels
       endif
   endif
   !write(*,*) 'cc'

   lev = lev + 1
   !write(*,*) 'lev = lev + 1'
   !write(*,*) lev

enddo model_levels

!write(*,*) 'after model_levels loop BG'
!write(*,*) lev
!write(*,*) num_levs

!Barre: here p_col is the pressure levels at the grid points, we will need to have in
!the future the pressure values at the mid-points (need to create a new function
!plevs_cam in model_mod using the the mid-levels hybrid coefs)
! KDR This comment looks confused.  plevs_cam returns pressures on what CESM refers to
!     as layer midpoints.  Maybe MOPITT needs pressures on the model interfaces,
!     in which case a new subroutine using hyai and hybi would replace hyam and hybm.
!     Anyway, p_col(30) is being redefined as the CAM surface pressure,
!     but the rest of the p_col is not redefined.
!     This all may be fine, since val is an accumulation of contributions 
!     from pressure sub-layers.

locS = set_location(mloc(1),mloc(2),0.0_r8, VERTISSURFACE)

!write(*,*) 'after set_location BG ## VERTISSURFACE'

call interpolate(state_handle, ens_size, locS, QTY_SURFACE_PRESSURE, p_col(:, num_levs), p_col_istatus)
call track_status(ens_size, p_col_istatus, p_col(:,num_levs), istatus, return_now)
if (return_now) return

! KDR Ack! redefining one element of a global array.
!> @todo not global array
!     This may be the 2nd redefinition of (1), or the first, if start_i /= 1.
mopitt_pres_local(:,1) = p_col(:,num_levs)

! KDR; Algorithm; 
!      Work through the MOPITT pressure level layers. (i loop)
!      Find CAM levels that are in each layer. (j loop)
!      Calculate CO from CAM state using level as the vertical location.
!      Accumulate contributions of CO from the CAM layers,
!         weighted by the the thickness of the CAM layer that lies within the MOPITT layer.
! Work through the MOPITT pressure level layers. (i loop)
! Loop from some number to 10; the mopitt pressure layers above the ground at this location.
! BG January 10 2019
! let's start over with V8

do i=start_i, end_i
   obs_val=0.0_r8

   bot_pres(:)=mopitt_pres_local(:, i)
   if (i == mopitt_dim) then
      top_pres(:)=mopitt_pres_local(:, i)/2.0_r8 
   else
      top_pres(:)=mopitt_pres_local(:, i+1)
   endif

   ! This is used in all 'coef's, below.
   ! first simplification  I agree it is used in all cases
   mop_layer_wght = 1.0_r8/abs(bot_pres-top_pres)

   !do j=10,num_levs
   !   i_bot_pres(:)=p_col(:,j)
   !   write(*,*) 'level j, i_bot_pres', j, i_bot_pres
   !enddo
   !write(*,*) '****************************************************** MOPITT obs '
   !write(*,*) 'MOPITT bot_pres :: ', bot_pres
   !write(*,*) 'MOPITT top_pres :: ', top_pres


   !  KDR Search through CAM layers for the ones which overlap the current mopitt layer.
   ! Yes
   do j=15,num_levs
      i_bot_pres(:)=p_col(:,j)
      if (j == 1) then
         i_top_pres(:)=0.0_r8
      else
         i_top_pres(:)=p_col(:,j-1)
      endif

      coef=0.0_r8
      obs_val_int=0.0_r8
      !write(*,*) 'level j, i_bot_pres ', j, i_bot_pres(1)
      !write(*,*) '         i_top_pres ', i_top_pres(1)
      ! purpose of the loop
      ! find the level pressures [i_bot_pres(iens) and i_top_pres(iens)]  that match MOPITT layers [bot_pres(iens) and top_pres(iens)] 
      ! and weight the respective CO accordingly
      if ( bot_pres(1) > top_pres(1) ) then
         if (i_bot_pres(1) <= bot_pres(1) .and. i_bot_pres(1) > top_pres(1)) then
            if (i_top_pres(1) <= top_pres(1)) then
               loc3 = set_location(mloc(1),mloc(2),real(j,r8), VERTISLEVEL)

               coef=abs(i_bot_pres(1)-top_pres(1))/abs(bot_pres(1)-top_pres(1))
               call interpolate(state_handle, ens_size, loc3, QTY_CO, obs_val_int, obs_val_int_istatus)

!               if (istatus /= 0) then
!                  if ( obs_val_int > 0. ) then
!                       istatus=0
!                  endif
!                       write(*,*) 'Case 1 obs_val_int :', obs_val_int, coef
!                       write(*,*) 'Case 1:', i_bot_pres, i_top_pres
!                       write(*,*) 'Case 1: / istatus', bot_pres, top_pres, istatus
!               endif
            endif
            if (i_top_pres(1) > top_pres(1)) then
               loc3 = set_location(mloc(1),mloc(2),real(j,r8), VERTISLEVEL)
               coef=abs(i_bot_pres(1)-i_top_pres(1))/abs(bot_pres(1)-top_pres(1))
               call interpolate(state_handle, ens_size, loc3, QTY_CO, obs_val_int, obs_val_int_istatus)

!               if (istatus /= 0) then
!                  if ( obs_val_int > 0. ) then
!                       istatus=0
!                  endif
!                 write(*,*) 'Case 2 obs_val_int :', obs_val_int, coef
!                 write(*,*) 'Case 2:', i_bot_pres(1), i_top_pres(1)
!                 write(*,*) 'Case 2: / istatus', bot_pres(1), top_pres(1), istatus
!               endif
            endif
         endif
         if (i_bot_pres(1) > bot_pres(1) .and. i_top_pres(1) < bot_pres(1)) then
            if (i_top_pres(1) <= top_pres(1)) then
               loc3 = set_location(mloc(1),mloc(2),real(j,r8), VERTISLEVEL)
               coef=abs(bot_pres(1)-top_pres(1))/abs(bot_pres(1)-top_pres(1))
               call interpolate(state_handle, ens_size, loc3, QTY_CO, obs_val_int, obs_val_int_istatus)
!              if (istatus /= 0) then
!                 if ( obs_val_int > 0. ) then
!                      istatus=0
!                 endif
!              write(*,*) 'Case 3 obs_val_int :', obs_val_int, coef
!              write(*,*) 'Case 3:', i_bot_pres(1), i_top_pres(1)
!              write(*,*) 'Case 3: / istatus ', bot_pres(1), top_pres(1), istatus
!              endif
            endif
            if (i_top_pres(1) > top_pres(1)) then
               loc3 = set_location(mloc(1),mloc(2),real(j,r8), VERTISLEVEL)
               coef=abs(bot_pres(1)-i_top_pres(1))/abs(bot_pres(1)-top_pres(1))
               call interpolate(state_handle, ens_size, loc3, QTY_CO, obs_val_int, obs_val_int_istatus)
!               if (istatus /= 0) then
!                  if ( obs_val_int > 0. ) then
!                       istatus=0
!                  endif
!              write(*,*) 'Case 4 obs_val_int :', obs_val_int, coef
!              write(*,*) 'Case 4:', i_bot_pres(1), i_top_pres(1)
!              write(*,*) 'Case 4: / istatus', bot_pres(1), top_pres(1), istatus
!               endif
            endif
         endif
      endif

      obs_val=obs_val+obs_val_int*coef
   enddo 
   where (obs_val > 0.0_r8 .and. avg_kernel(key,i) > -700) val = val + avg_kernel(key,i) * log10(obs_val)
enddo

where (istatus == 0) val = val + mopitt_prior(key)
where (istatus == 0) val=10.0**val

where (val < 0.0_r8)
   val=MISSING_R8
   istatus = 6
end where

end subroutine get_expected_mopitt_co


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
