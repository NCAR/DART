! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

! BEGIN DART PREPROCESS KIND LIST
! IASI_O3_RETRIEVAL, KIND_O3
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_iasi_mod, only : write_iasi_o3, read_iasi_o3, &
!                                  interactive_iasi_o3, get_expected_iasi_o3, &
!                                  set_obs_def_iasi_o3
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(IASI_O3_RETRIEVAL)                                                           
!            call get_expected_iasi_o3(state, location, obs_def%key, obs_val, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(IASI_O3_RETRIEVAL)
!         call read_iasi_o3(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(IASI_O3_RETRIEVAL)
!         call write_iasi_o3(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(IASI_O3_RETRIEVAL)
!         call interactive_iasi_o3(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS SET_OBS_DEF_IASI_O3
!      case(IASI_O3_RETRIEVAL)
!         call set_obs_def_iasi_o3(obs_def%key)
! END DART PREPROCESS SET_OBS_DEF_IASI_O3


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_iasi_mod

use        types_mod, only : r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, set_location, get_location, VERTISHEIGHT,&
                             VERTISPRESSURE, VERTISLEVEL
use  assim_model_mod, only : interpolate
use    obs_kind_mod, only  : KIND_O3, KIND_O3_COLUMN
!AFAJ
use mpi_utilities_mod, only : my_task_id  
implicit none

public :: write_iasi_o3, read_iasi_o3, interactive_iasi_o3, &
          get_expected_iasi_o3, set_obs_def_iasi_o3

! Storage for the special information required for observations of this type
integer, parameter               :: max_iasi_o3_obs = 1000000
integer                          :: num_iasi_o3_obs = 0

! nominal iasi number of levels
integer, parameter                        :: iasi_dim = 40

! number of levels used
integer, dimension(max_iasi_o3_obs)     :: iasi_nlevels

! mopitt averaging kernel --note that this is transformed and scaled
real(r8), dimension(max_iasi_o3_obs,iasi_dim) :: avg_kernel
! prior term of x=Ax + (I-A)xa + Gey
real(r8), dimension(max_iasi_o3_obs)    :: iasi_o3_prior
!real(r8), dimension(max_iasi_o3_obs)    :: iasi_air_column

! nominal iasi height levels in m
real(r8)                                  :: iasi_altitude(iasi_dim) =(/ &
                                             500.,1500.,2500.,3500.,4500., &
                                             5500.,6500.,7500.,8500.,9500., &
                                             10500.,11500.,12500.,13500.,14500., &
                                             15500.,16500.,17500.,18500.,19500., &
                                             20500.,21500.,22500.,23500.,24500., &
                                             25500.,26500.,27500.,28500.,29500., &
                                             30500.,31500.,32500.,33500.,34500., &
                                             35500.,36500.,37500.,38500.,39500. /) 
! iasi retrieval  heights
real(r8), dimension(max_iasi_o3_obs,iasi_dim) :: iasi_heights

! iasi retrieval  prior profile
real(r8), dimension(max_iasi_o3_obs,iasi_dim) :: iasi_prior_prof

! iasi retrieval  air column profile
real(r8), dimension(max_iasi_o3_obs,iasi_dim) :: iasi_air_column


! For now, read in all info on first read call, write all info on first write call
logical :: already_read = .false., already_written = .false.

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source: /home/thoar/CVS.REPOS/DART/obs_def/obs_def_iasi_mod.f90,v $", &
revision = "$Revision: 1.1 $", &
revdate  = "$Date: 2005/10/05 15:19:28 $"

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



 subroutine read_iasi_o3(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine read_iasi_o3(key, ifile, fform)

integer, intent(out)            :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional    :: fform
character(len=32) 		:: fileformat

! temp variables
real(r8)                        :: prior_1
integer                         :: nlevel_1
real(r8),  dimension(iasi_dim)  :: avg_kernel_1
real(r8),  dimension(iasi_dim)  :: altitude_1
real(r8),  dimension(iasi_dim)  :: prior_prof_1
real(r8),  dimension(iasi_dim)  :: aircol_1
integer 			:: keyin

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading

nlevel_1 = read_iasi_num_levels(ifile, fileformat)
prior_1  = read_iasi_prior_column(ifile, fileformat)
aircol_1  = read_iasi_air_column(ifile, nlevel_1,fileformat)
altitude_1 = read_iasi_heights(ifile, nlevel_1, fileformat)
avg_kernel_1 = read_iasi_avg_kernels(ifile, nlevel_1, fileformat) 
prior_prof_1 = read_iasi_prior_prof(ifile, nlevel_1, fileformat) 

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   read(ifile) keyin

   CASE DEFAULT
   read(ifile, *) keyin
END SELECT

counts1 = counts1 + 1
key = counts1
call set_obs_def_iasi_o3(key,avg_kernel_1(1:nlevel_1),prior_1,altitude_1(1:nlevel_1), &
                         aircol_1, nlevel_1, prior_prof_1)

end subroutine read_iasi_o3


 subroutine write_iasi_o3(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine write_iasi_o3(key, ifile, fform)

integer, intent(in)             :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional 	:: fform
character(len=32) 		:: fileformat

! dummy variables
real(r8), dimension(iasi_dim)   :: avg_kernel_1
real(r8), dimension(iasi_dim)   :: altitude_1
real(r8), dimension(iasi_dim)   :: prior_prof_1
real(r8), dimension(iasi_dim)   :: iasi_air_column_1
integer                         :: nlevel_1

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading
nlevel_1 = iasi_nlevels(key)
avg_kernel_1(1:nlevel_1) = avg_kernel(key,:)
altitude_1(1:nlevel_1) = iasi_heights(key,:)
prior_prof_1(1:nlevel_1) = iasi_prior_prof(key,:)
iasi_air_column_1(1:nlevel_1) = iasi_air_column(key,:)

call write_iasi_num_levels(ifile, nlevel_1, fileformat)
call write_iasi_prior_column(ifile, iasi_o3_prior(key), fileformat)
call write_iasi_air_column(ifile, iasi_air_column_1(1:nlevel_1), nlevel_1, fileformat)
call write_iasi_heights(ifile, altitude_1(1:nlevel_1), nlevel_1, fileformat)
call write_iasi_avg_kernels(ifile, avg_kernel_1(1:nlevel_1), nlevel_1, fileformat)
call write_iasi_prior_prof(ifile, prior_prof_1(1:nlevel_1), nlevel_1, fileformat)

   
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   write(ifile) key

   CASE DEFAULT
   write(ifile, *) key
END SELECT 

end subroutine write_iasi_o3


 subroutine interactive_iasi_o3(key)
!----------------------------------------------------------------------
!subroutine interactive_iasi_o3(key)
!
! Initializes the specialized part of a IASI observation
! Passes back up the key for this one

integer, intent(out) :: key

character(len=129) :: msgstring

if ( .not. module_initialized ) call initialize_module

! Make sure there's enough space, if not die for now (clean later)
if(num_iasi_o3_obs >= max_iasi_o3_obs) then
   ! PUT IN ERROR HANDLER CALL
   write(msgstring, *)'Not enough space for a iasi O3 obs.'
   call error_handler(E_MSG,'interactive_iasi_o3',msgstring,source,revision,revdate)
   write(msgstring, *)'Can only have max_iasi_o3_obs (currently ',max_iasi_o3_obs,')'
   call error_handler(E_ERR,'interactive_iasi_o3',msgstring,source,revision,revdate)
endif

! Increment the index
num_iasi_o3_obs = num_iasi_o3_obs + 1
key = num_iasi_o3_obs

! Otherwise, prompt for input for the three required beasts
write(*, *) 'Creating an interactive_iasi_o3 observation'

end subroutine interactive_iasi_o3



 subroutine get_expected_iasi_o3(state, location, key, val, istatus)
!----------------------------------------------------------------------
!subroutine get_expected_iasi_o3(state, location, key, val, istatus)

real(r8), intent(in)            :: state(:)
type(location_type), intent(in) :: location
integer, intent(in)             :: key
real(r8), intent(out)           :: val
integer, intent(out)            :: istatus

integer :: i, ilev, valid_lev
type(location_type) :: loc2
real(r8)            :: mloc(3)
real(r8)	    :: obs_val, tmp_val, level

integer :: nlevels
if ( .not. module_initialized ) call initialize_module

mloc = get_location(location)
 
val = 0.0_r8
if (mloc(2)>90.0_r8) then
    mloc(2)=90.0_r8
elseif (mloc(2)<-90.0_r8) then
    mloc(2)=-90.0_r8
endif

! Apply IASI Averaging kernel A and IASI Prior (I-A)xa
! x = Axm + (I-A)xa , where x is a 40-element vector 

nlevels = iasi_nlevels(key)
obs_val = 0.0_r8
istatus  = 0
valid_lev  = 0
level = 1.0_r8

! loop through all levels
do ilev = 1, nlevels

    !get location of obs
    if (ilev == 1) then
       loc2 = set_location(mloc(1),mloc(2),level, VERTISLEVEL)
    else
       loc2 = set_location(mloc(1),mloc(2),iasi_heights(key,ilev), VERTISHEIGHT)
    endif

    !initialize every level
    obs_val = 0.0_r8
    istatus = 0

    ! interpolate to obs location
    call interpolate(state, loc2, KIND_O3, obs_val, istatus)  

    ! check for problems 
    if (istatus /= 0) then
        if (istatus == 2) then
            ! it looks like it's a problem with vertical level interpolation
            ! use the prior iasi sub column instead
            tmp_val = iasi_prior_prof(key,ilev)
            istatus = 0
        else
            ! interpolation failed (outside of domain?)
            tmp_val = 0.0_r8
            istatus = istatus
            return
        endif
    else
        ! assign interpolated subcolumn (multiply by air column and 1000)
        if ( obs_val > 0.0_r8 ) then
            tmp_val = obs_val * iasi_air_column(key,ilev) * 1000.0_r8
            istatus = 0
        else
            ! the interpolated value cannot be negative
            tmp_val = 0.0_r8
            istatus = 25
            return
        endif 
    endif
    ! update expected obs with the averaging kernel
    val = val + avg_kernel(key,ilev) * tmp_val
    !if (my_task_id() == 0 ) then
    !endif
enddo
val = val + iasi_o3_prior(key)  

end subroutine get_expected_iasi_o3
!----------------------------------------------------------------------

 subroutine set_obs_def_iasi_o3(key, akcol, apcol_val, altretlev, aircol_val, nlev_use, prior_prof)
!----------------------------------------------------------------------
! Allows passing of obs_def special information 

integer,	 	intent(in)	:: key
character(len=129) 			:: msgstring

real(r8),               intent(in)      :: apcol_val
real*8, dimension(40),  intent(in)      :: akcol
real*8, dimension(40),  intent(in)      :: aircol_val
real*8, dimension(40),  intent(in)      :: altretlev
real*8, dimension(40),  intent(in)      :: prior_prof
integer,                intent(in)      :: nlev_use

if ( .not. module_initialized ) call initialize_module

if(num_iasi_o3_obs >= max_iasi_o3_obs) then
   ! PUT IN ERROR HANDLER CALL
   write(msgstring, *)'Not enough space for a iasi O3 obs.'
   call error_handler(E_MSG,'set_obs_def_iasi_o3',msgstring,source,revision,revdate)
   write(msgstring, *)'Can only have max_iasi_o3_obs (currently ',max_iasi_o3_obs,')'
   call error_handler(E_ERR,'set_obs_def_iasi_o3',msgstring,source,revision,revdate)
endif

avg_kernel(key, 1:nlev_use)   		= akcol(1:nlev_use)
iasi_prior_prof(key, 1:nlev_use)	= prior_prof(1:nlev_use)
iasi_o3_prior(key)			= apcol_val
iasi_air_column(key,1:nlev_use)		= aircol_val(1:nlev_use)
iasi_heights(key,1:nlev_use)		= altretlev
iasi_nlevels(key)			= nlev_use

end subroutine set_obs_def_iasi_o3


!=================================
! other functions and subroutines
!=================================
function read_iasi_prior_column(ifile, fform)

integer,                    intent(in) :: ifile
real(r8)                               :: read_iasi_prior_column
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_prior_column
   CASE DEFAULT
      read(ifile, *) read_iasi_prior_column
END SELECT

end function read_iasi_prior_column

subroutine write_iasi_prior_column(ifile, iasi_prior_temp, fform)

integer,                    intent(in) :: ifile
real(r8),                   intent(in) :: iasi_prior_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) iasi_prior_temp
   CASE DEFAULT
      write(ifile, *) iasi_prior_temp
END SELECT

end subroutine write_iasi_prior_column

function read_iasi_num_levels(ifile, fform)

integer,                    intent(in) :: ifile
integer                                :: read_iasi_num_levels
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_num_levels
   CASE DEFAULT
      read(ifile, *) read_iasi_num_levels
END SELECT

end function read_iasi_num_levels

subroutine write_iasi_num_levels(ifile, number_of_levels_temp, fform)

integer,                    intent(in) :: ifile
integer,                    intent(in) :: number_of_levels_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) number_of_levels_temp
   CASE DEFAULT
      write(ifile, *) number_of_levels_temp
END SELECT

end subroutine write_iasi_num_levels

function read_iasi_avg_kernels(ifile, nlevels,fform)

integer,                    intent(in) :: ifile, nlevels
real(r8), dimension(40)        :: read_iasi_avg_kernels
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat
   
if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_avg_kernels(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_iasi_avg_kernels(1:nlevels)
END SELECT

end function read_iasi_avg_kernels

function read_iasi_heights(ifile, nlevels,fform)

integer,                    intent(in) :: ifile, nlevels
real(r8), dimension(40)        :: read_iasi_heights
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_heights(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_iasi_heights(1:nlevels)
END SELECT

end function read_iasi_heights

function read_iasi_prior_prof(ifile, nlevels,fform)

integer,                    intent(in) :: ifile, nlevels
real(r8), dimension(40)        :: read_iasi_prior_prof
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_prior_prof(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_iasi_prior_prof(1:nlevels)
END SELECT

end function read_iasi_prior_prof

function read_iasi_air_column(ifile, nlevels,fform)

integer,                    intent(in) :: ifile, nlevels
real(r8), dimension(40)        :: read_iasi_air_column
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_air_column(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_iasi_air_column(1:nlevels)
END SELECT

end function read_iasi_air_column




subroutine write_iasi_avg_kernels(ifile, avg_kernels_temp, nlevels, fform)

integer,                    intent(in) :: ifile, nlevels
real(r8), dimension(40), intent(in)  :: avg_kernels_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) avg_kernels_temp(1:nlevels)
   CASE DEFAULT
      write(ifile, *) avg_kernels_temp(1:nlevels)
END SELECT

end subroutine write_iasi_avg_kernels

subroutine write_iasi_heights(ifile, height_temp, nlevels, fform)

integer,                    intent(in) :: ifile, nlevels
real(r8), dimension(40), intent(in)  :: height_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) height_temp(1:nlevels)
   CASE DEFAULT
      write(ifile, *) height_temp(1:nlevels)
END SELECT

end subroutine write_iasi_heights

subroutine write_iasi_prior_prof(ifile, prior_prof_temp, nlevels, fform)

integer,                    intent(in) :: ifile, nlevels
real(r8), dimension(40), intent(in)  :: prior_prof_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) prior_prof_temp(1:nlevels)
   CASE DEFAULT
      write(ifile, *) prior_prof_temp(1:nlevels)
END SELECT

end subroutine write_iasi_prior_prof

subroutine write_iasi_air_column(ifile, aircol_prof_temp, nlevels, fform)

integer,                    intent(in) :: ifile, nlevels
real(r8), dimension(40), intent(in)  :: aircol_prof_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) aircol_prof_temp(1:nlevels)
   CASE DEFAULT
      write(ifile, *) aircol_prof_temp(1:nlevels)
END SELECT

end subroutine write_iasi_air_column








end module obs_def_iasi_mod
! END DART PREPROCESS MODULE CODE
