! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module platform_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, read_location, write_location, interactive_location

implicit none
private

interface assignment(=)
   module procedure copy_platform
end interface

public platform_type, set_platform_location, copy_platform, write_platform, assignment(=)

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type platform_type
   private
   type(location_type)   :: location
end type platform_type

logical, save :: module_initialized = .false.

contains


  subroutine initialize_module
!----------------------------------------------------------------------------

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module


!---------------------------------------------------------------------

subroutine copy_platform(platform1, platform2)

! Copy function to be overloaded with '='

type(platform_type), intent(out) :: platform1
type(platform_type), intent(in) :: platform2

if ( .not. module_initialized ) call initialize_module

platform1%location = platform2%location

end subroutine copy_platform

!----------------------------------------------------------------------------

subroutine set_platform_location(platform, location)

! Sets the location of a platform

type(platform_type), intent(inout) :: platform
type(location_type), intent(in)    :: location

if ( .not. module_initialized ) call initialize_module

platform%location = location

end subroutine set_platform_location

!-------------------------------------------------------

!----------------------------------------------------------------------------

subroutine write_platform(ifile, platform, fform)

implicit none

integer,                    intent(in) :: ifile
type(platform_type),        intent(in) :: platform
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
      write(ifile, 11)
11    format('platform')
END SELECT

call write_location(ifile, platform%location, fform)

end subroutine write_platform

end module platform_mod
