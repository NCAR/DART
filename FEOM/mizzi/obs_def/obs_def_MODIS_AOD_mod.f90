! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

! BEGIN DART PREPROCESS KIND LIST
! MODIS_AOD_RETRIEVAL, KIND_AOD
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_modis_mod, only : write_modis_aod, read_modis_aod, &
!                                  interactive_modis_aod, get_expected_modis_aod, &
!                                  set_obs_def_modis_aod
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(MODIS_AOD_RETRIEVAL)                                                           
!            call get_expected_modis_aod(state, location, obs_def%key, obs_val, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(MODIS_AOD_RETRIEVAL)
!         call read_modis_aod(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(MODIS_AOD_RETRIEVAL)
!         call write_modis_aod(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(MODIS_AOD_RETRIEVAL)
!         call interactive_modis_aod(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS SET_OBS_DEF_MODIS_AOD
!      case(MODIS_AOD_RETRIEVAL)
!         call set_obs_def_modis_aod(obs_def%key)
! END DART PREPROCESS SET_OBS_DEF_MODIS_AOD


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_modis_mod

use        types_mod, only : r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, set_location, get_location, VERTISPRESSURE, VERTISUNDEF

use  assim_model_mod, only : interpolate
use    obs_kind_mod, only  : KIND_AOD

implicit none

public :: write_modis_aod, read_modis_aod, interactive_modis_aod, &
          get_expected_modis_aod, set_obs_def_modis_aod

! Storage for the special information required for observations of this type
integer, parameter               :: max_modis_aod_obs = 10000000
integer                          :: num_modis_aod_obs = 0


! For now, read in all info on first read call, write all info on first write call
logical :: already_read = .false., already_written = .false.

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source: /home/thoar/CVS.REPOS/DART/obs_def/obs_def_modis_mod.f90,v $", &
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



 subroutine read_modis_aod(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine read_modis_aod(key, ifile, fform)

integer, intent(out)            :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional    :: fform
character(len=32) 		:: fileformat

integer 			:: keyin

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   read(ifile) keyin

   CASE DEFAULT
   read(ifile, *) keyin
END SELECT

counts1 = counts1 + 1
key = counts1
call set_obs_def_modis_aod(key)

end subroutine read_modis_aod


 subroutine write_modis_aod(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine write_modis_aod(key, ifile, fform)

integer, intent(in)             :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional 	:: fform

character(len=32) 		:: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading
   
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   write(ifile) key

   CASE DEFAULT
   write(ifile, *) key
END SELECT 

end subroutine write_modis_aod


 subroutine interactive_modis_aod(key)
!----------------------------------------------------------------------
!subroutine interactive_modis_aod(key)
!
! Initializes the specialized part of a MODIS observation
! Passes back up the key for this one

integer, intent(out) :: key

character(len=129) :: msgstring

if ( .not. module_initialized ) call initialize_module

! Make sure there's enough space, if not die for now (clean later)
if(num_modis_aod_obs >= max_modis_aod_obs) then
   ! PUT IN ERROR HANDLER CALL
   write(msgstring, *)'Not enough space for a modis AOD obs.'
   call error_handler(E_MSG,'interactive_modis_aod',msgstring,source,revision,revdate)
   write(msgstring, *)'Can only have max_modis_aod_obs (currently ',max_modis_aod_obs,')'
   call error_handler(E_ERR,'interactive_modis_aod',msgstring,source,revision,revdate)
endif

! Increment the index
num_modis_aod_obs = num_modis_aod_obs + 1
key = num_modis_aod_obs

! Otherwise, prompt for input for the three required beasts
write(*, *) 'Creating an interactive_modis_aod observation'

end subroutine interactive_modis_aod



 subroutine get_expected_modis_aod(state, location, key, val, istatus)
!----------------------------------------------------------------------
!subroutine get_expected_modis_aod(state, location, key, val, istatus)

real(r8), intent(in)            :: state(:)
type(location_type), intent(in) :: location
integer, intent(in)             :: key
real(r8), intent(out)           :: val
integer, intent(out)            :: istatus

integer :: i
type(location_type) :: loc2
real(r8)            :: mloc(3)
real(r8)	    :: obs_val

if ( .not. module_initialized ) call initialize_module

mloc = get_location(location)
 
val = 0.0_r8
if (mloc(2)>90.0_r8) then
    mloc(2)=90.0_r8
elseif (mloc(2)<-90.0_r8) then
    mloc(2)=-90.0_r8
endif
loc2 = set_location(mloc(1),mloc(2),mloc(3), VERTISUNDEF)
call interpolate(state, loc2, KIND_AOD, obs_val, istatus)  
if (istatus /= 0) then
    val = 0
    return
endif
val = obs_val  
!print *, 'AFAJ DEBUG AOD VAL ', val, istatus

end subroutine get_expected_modis_aod
!----------------------------------------------------------------------

 subroutine set_obs_def_modis_aod(key)
!----------------------------------------------------------------------
! Allows passing of obs_def special information 

integer,	 	intent(in)	:: key
character(len=129) 			:: msgstring

if ( .not. module_initialized ) call initialize_module

if(num_modis_aod_obs >= max_modis_aod_obs) then
   ! PUT IN ERROR HANDLER CALL
   write(msgstring, *)'Not enough space for a modis AOD obs.'
   call error_handler(E_MSG,'set_obs_def_modis_aod',msgstring,source,revision,revdate)
   write(msgstring, *)'Can only have max_modis_aod_obs (currently ',max_modis_aod_obs,')'
   call error_handler(E_ERR,'set_obs_def_modis_aod',msgstring,source,revision,revdate)
endif

end subroutine set_obs_def_modis_aod


end module obs_def_modis_mod
! END DART PREPROCESS MODULE CODE
