! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_def_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

use        types_mod, only : r8
use   utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     obs_kind_mod, only : obs_kind_type, read_kind, write_kind, interactive_kind
use     location_mod, only : location_type, read_location, write_location, interactive_location
use time_manager_mod, only : time_type, read_time, write_time, set_time

implicit none
private

interface assignment(=)
   module procedure copy_obs_def
end interface

public init_obs_def, get_obs_def_location, get_obs_def_kind, get_obs_def_time, &
   get_obs_def_error_variance, set_obs_def_location, set_obs_def_kind, set_obs_def_time, &
   set_obs_def_error_variance, interactive_obs_def, write_obs_def, read_obs_def, obs_def_type , &
   destroy_obs_def, copy_obs_def, assignment(=)
!!!   get_obs_platform, get_obs_platform_qc, get_obs_aperture

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type obs_def_type
! In revision, obs_kind module is responsibe for taking care of identity obs kinds, too
   private
   type(location_type) :: location      ! center of mass, so to speak
   type(obs_kind_type) :: kind          ! keyword, BUFR values for now
   type(time_type)     :: time
   real(r8)            :: error_variance
!  type(platform_type) :: platform
!  type(instrument_type) :: instrument
!  type(instrument_name) :: instrument_name
!  real(r8), pointer :: platform_qc(:)
!  type(aperture_type) :: aperture
end type obs_def_type

logical, save :: module_initialized = .false.

contains


  subroutine initialize_module
!----------------------------------------------------------------------------
! subroutine initialize_module

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module


!----------------------------------------------------------------------------

subroutine init_obs_def(obs_def, location, kind, time, error_variance)
! Need to add additional component arguments as optionals as needed

! Constructor for an obs_def 

type(obs_def_type), intent(out) :: obs_def
type(location_type), intent(in) :: location
type(obs_kind_type), intent(in) :: kind
type(time_type), intent(in) :: time
real(r8), intent(in) :: error_variance

if ( .not. module_initialized ) call initialize_module

obs_def%location = location
obs_def%kind = kind
obs_def%time = time
obs_def%error_variance = error_variance

end subroutine init_obs_def

!---------------------------------------------------------------------

subroutine copy_obs_def(obs_def1, obs_def2)

! Copy function to be overloaded with '='

type(obs_def_type), intent(out) :: obs_def1
type(obs_def_type), intent(in) :: obs_def2

if ( .not. module_initialized ) call initialize_module

obs_def1%location = obs_def2%location
obs_def1%kind = obs_def2%kind
obs_def1%time = obs_def2%time
obs_def1%error_variance = obs_def2%error_variance
!obs_def1%platform = obs_def2% platform
!deallocate(obs_def1%platform_qc)
!allocate(obs_def1%platform_qc(size(obs_def2%platform_qc))
! Should this be pointer assignment or regular
!obs_def1%platform_qc >= or == obs_def2%platform_qc
!obs_def1%aperture = obs_def2%aperture

end subroutine copy_obs_def

!----------------------------------------------------------------------------

function get_obs_def_error_variance(obs_def)

type(obs_def_type), intent(in) :: obs_def
real(r8) :: get_obs_def_error_variance

if ( .not. module_initialized ) call initialize_module

get_obs_def_error_variance = obs_def%error_variance

end function get_obs_def_error_variance

!----------------------------------------------------------------------------

function get_obs_def_location(obs_def)

! Returns observation location.

type(location_type) :: get_obs_def_location
type(obs_def_type), intent(in) :: obs_def

if ( .not. module_initialized ) call initialize_module

get_obs_def_location = obs_def%location

end function get_obs_def_location

!----------------------------------------------------------------------------

function get_obs_def_kind(obs_def)

! Returns observation kind

type(obs_kind_type) :: get_obs_def_kind
type(obs_def_type), intent(in) :: obs_def

if ( .not. module_initialized ) call initialize_module

get_obs_def_kind = obs_def%kind

end function get_obs_def_kind

!----------------------------------------------------------------------------

function get_obs_def_time(obs_def)

! Returns observation time

type(time_type) :: get_obs_def_time
type(obs_def_type), intent(in) :: obs_def

if ( .not. module_initialized ) call initialize_module

get_obs_def_time = obs_def%time

end function get_obs_def_time

!----------------------------------------------------------------------------

subroutine set_obs_def_location(obs_def, location)

! Sets the location of an obs_def

type(obs_def_type), intent(inout) :: obs_def
type(location_type), intent(in) :: location

if ( .not. module_initialized ) call initialize_module

obs_def%location = location

end subroutine set_obs_def_location

!----------------------------------------------------------------------------

subroutine set_obs_def_error_variance(obs_def, error_variance)

! Sets the error variance of an obs_def

type(obs_def_type), intent(inout) :: obs_def
real(r8), intent(in) :: error_variance

if ( .not. module_initialized ) call initialize_module

obs_def%error_variance = error_variance

end subroutine set_obs_def_error_variance

!----------------------------------------------------------------------------

subroutine set_obs_def_kind(obs_def, kind)

! Sets the kind of an obs_def

type(obs_def_type), intent(inout) :: obs_def
type(obs_kind_type), intent(in) :: kind

if ( .not. module_initialized ) call initialize_module

obs_def%kind = kind

end subroutine set_obs_def_kind

!----------------------------------------------------------------------------

subroutine set_obs_def_time(obs_def, time)

! Sets the time of an obs_def

type(obs_def_type), intent(inout) :: obs_def
type(time_type), intent(in) :: time

if ( .not. module_initialized ) call initialize_module

obs_def%time = time

end subroutine set_obs_def_time

!----------------------------------------------------------------------------

subroutine read_obs_def(file, obs_def, fform)

! Reads an obs_def from file which is just an integer unit number in the 
! current preliminary implementation.

type(obs_def_type), intent(inout) :: obs_def
integer, intent(in) :: file
character(len=*), intent(in), optional :: fform

character(len=5) :: header
character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Begin by reading five character ascii header, then location, kind, error variance, index

! Need to add additional error checks on read
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      continue
   CASE DEFAULT
      read(file, 11) header
11    Format(a5)
      if(header /= 'obdef') then
   call error_handler(E_ERR,'read_obs_def', 'Expected location header "obdef" in input file', &
                      source, revision, revdate)
      endif
END SELECT

! Read the location, kind and error variance
obs_def%location = read_location(file, fform)
obs_def%kind = read_kind(file, fform)
obs_def%time = read_time(file, fform)

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      read(file) obs_def%error_variance
   CASE DEFAULT
      read(file, *) obs_def%error_variance
END SELECT

end subroutine read_obs_def

!----------------------------------------------------------------------------

subroutine write_obs_def(file, obs_def, fform)

! Writes an obs_def to file.

integer, intent(in) :: file
type(obs_def_type), intent(in) :: obs_def
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Write the 5 character identifier for verbose formatted output
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      continue
   CASE DEFAULT
      write(file, 11)
11    format('obdef')
END SELECT

! Write out the location, kind and error variance
call write_location(file, obs_def%location, fform)
call write_kind(file, obs_def%kind, fform)
call write_time(file, obs_def%time, fform)

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      write(file) obs_def%error_variance
   CASE DEFAULT
      write(file, *) obs_def%error_variance
END SELECT

end subroutine write_obs_def


subroutine interactive_obs_def(obs_def)
!---------------------------------------------------------------------------
!
! Allows interactive creation of an observation

type(obs_def_type), intent(inout) :: obs_def

integer :: ind
real(r8) :: error_variance
integer :: seconds, days

if ( .not. module_initialized ) call initialize_module

! Get the location
call interactive_location(obs_def%location)

! Get the obsrvation kind
call interactive_kind(obs_def%kind)

! Get the time
!!!call interactive_time(obs_def%time)
! Eventually this should be done in time manager to allow for calendar type use
! but this could be done after ESMF adoption?
write(*, *) 'input time in seconds and days '
read(*, *) seconds, days
obs_def%time = set_time(seconds, days)

write(*, *) 'Input error variance for this observation definition '
read(*, *) obs_def%error_variance

end subroutine interactive_obs_def

!----------------------------------------------------------------

subroutine destroy_obs_def(obs_def)

type(obs_def_type), intent(in) :: obs_def

if ( .not. module_initialized ) call initialize_module

! TECHNICALLY NEED TO CALL DESTRUCTORS FOR ALL SUBCOMPONENTS, NO ALLOCATED STORAGE YET

end subroutine destroy_obs_def

!-------------------------------------------------------

end module obs_def_mod
