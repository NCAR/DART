! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module platform_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

use     types_mod, only : r8
use utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use  location_mod, only : location_type, read_location, write_location, interactive_location

implicit none
private

interface assignment(=)
   module procedure copy_platform
end interface

public :: platform_type, read_platform, copy_platform, write_platform, assignment(=), &
       set_platform_location, get_platform_location, &
       set_platform_orientation, get_platform_orientation

public :: write_orientation, read_orientation

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type platform_type
   private
   type(location_type) :: location
   real(r8)            :: dir(3)
end type platform_type

logical, save :: module_initialized = .false.

contains

!---------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!---------------------------------------------------------------------

subroutine copy_platform(platform1, platform2)

! Copy function to be overloaded with '='

type(platform_type), intent(out) :: platform1
type(platform_type), intent(in)  :: platform2

if ( .not. module_initialized ) call initialize_module

platform1%location = platform2%location

platform1%dir(1) = platform2%dir(1)
platform1%dir(2) = platform2%dir(2)
platform1%dir(3) = platform2%dir(3)

end subroutine copy_platform

!---------------------------------------------------------------------

subroutine set_platform_location(platform, location)

! Sets the location of a platform

type(platform_type), intent(inout) :: platform
type(location_type), intent(in)    :: location

if ( .not. module_initialized ) call initialize_module

platform%location = location

end subroutine set_platform_location

!---------------------------------------------------------------------

function get_platform_location(platform)

! Sets the location of a platform

type(platform_type), intent(in) :: platform
type(location_type)             :: get_platform_location

if ( .not. module_initialized ) call initialize_module

get_platform_location = platform%location

end function get_platform_location

!---------------------------------------------------------------------

subroutine set_platform_orientation(platform, dir)

! Sets the location of a platform

type(platform_type),    intent(inout) :: platform
real(r8),               intent(in)    :: dir(3)

if ( .not. module_initialized ) call initialize_module

platform%dir(1) = dir(1)
platform%dir(2) = dir(2)
platform%dir(3) = dir(3)

end subroutine set_platform_orientation

!---------------------------------------------------------------------

function get_platform_orientation(platform)

! Sets the orientation of a platform

type(platform_type), intent(in) :: platform
real(r8)                        :: get_platform_orientation(3)

if ( .not. module_initialized ) call initialize_module

get_platform_orientation(1) = platform%dir(1)
get_platform_orientation(2) = platform%dir(2)
get_platform_orientation(3) = platform%dir(3)

end function get_platform_orientation

!---------------------------------------------------------------------

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

call write_location(ifile, platform%location, fileformat)
call write_orientation(ifile, platform%dir, fileformat)

end subroutine write_platform

!---------------------------------------------------------------------

subroutine read_platform(ifile, platform, fform)

implicit none

integer,                    intent(in)  :: ifile
type(platform_type),        intent(out) :: platform
character(len=*), intent(in), optional  :: fform

character(len=8)  :: header
character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Read the character identifier for verbose formatted output
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      continue
   CASE DEFAULT
      read(ifile, FMT='(a8)') header
      if(header /= 'platform') then
         call error_handler(E_ERR,'read_platform', &
              'Expected location header "platform" in input file', &
              source, revision, revdate)
      endif
END SELECT

platform%location = read_location(ifile, fileformat)
platform%dir = read_orientation(ifile, fileformat)

end subroutine read_platform

subroutine write_orientation(ifile, dir, fform)
!----------------------------------------------------------------------------
!

implicit none

integer,                    intent(in) :: ifile
real(r8),                   intent(in) :: dir(3)
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) dir(1), dir(2), dir(3)
   CASE DEFAULT
      write(ifile, '(''dir3d'')' ) 
      write(ifile, *) dir(1), dir(2), dir(3)
END SELECT

end subroutine write_orientation


function read_orientation(ifile, fform)
!----------------------------------------------------------------------------
!
! Reads orientation from file that was written by write_orientation.
! See write_orientation for additional discussion.

implicit none

integer,                    intent(in) :: ifile
real(r8)                               :: read_orientation(3)
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_orientation(1), read_orientation(2), read_orientation(3)
   CASE DEFAULT
      read(ifile, '(a5)' ) header

      if(header /= 'dir3d') then
         write(errstring,*)'Expected orientation header "dir3d" in input file, got ', header
         call error_handler(E_ERR, 'read_orientation', errstring, source, revision, revdate)
      endif
! Now read the orientation data value
      read(ifile, *) read_orientation(1), read_orientation(2), read_orientation(3)
END SELECT


end function read_orientation


end module platform_mod
