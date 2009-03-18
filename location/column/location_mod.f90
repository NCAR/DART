! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module location_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! Implements location interfaces for a one dimensional column domain

use      types_mod, only : r8, MISSING_R8, MISSING_I
use  utilities_mod, only : register_module, error_handler, E_ERR

implicit none
private

public :: location_type, get_dist, get_location, set_location, &
          set_location2, set_location_missing, is_location_in_region, &
          write_location, read_location, interactive_location, vert_is_undef, &
          vert_is_surface, vert_is_pressure, vert_is_level, vert_is_height, &
          query_location, LocationDims, LocationName, LocationLName, &
          alloc_get_close_obs, get_close_obs, &
          operator(==), operator(/=)

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

! The possible values for the location_type%which_vert component.
! These are intended to be PRIVATE to this module. Do not make public.

integer, parameter :: VERTISUNDEF    = -2 ! has no vertical location (undefined)
integer, parameter :: VERTISSURFACE  = -1 ! surface value
integer, parameter :: VERTISLEVEL    =  1 ! by level
integer, parameter :: VERTISPRESSURE =  2 ! by pressure
integer, parameter :: VERTISHEIGHT   =  3 ! by height

type location_type
   private
   real(r8) :: vloc
   integer  :: which_vert
end type location_type

logical, save :: module_initialized = .false.

integer,              parameter :: LocationDims = 1
character(len = 129), parameter :: LocationName = "loc1d"
character(len = 129), parameter :: LocationLName = "one-dimensional location"

interface operator(==); module procedure loc_eq; end interface
interface operator(/=); module procedure loc_ne; end interface

contains

  subroutine initialize_location
!----------------------------------------------------------------------------
! subroutine initialize_location
!
! pretty simple for this module.

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_location



function get_dist(loc1, loc2)
!----------------------------------------------------------------------------

implicit none

type(location_type), intent(in) :: loc1, loc2
real(r8) :: get_dist

if ( .not. module_initialized ) call initialize_location

! Reentrant domain, if distance is greater than half wraparound the other way.
get_dist = abs(loc1%vloc - loc2%vloc)

end function get_dist



function get_location(loc)
!---------------------------------------------------------------------------
!
! Given a location type, return the z coordinate

implicit none

type(location_type), intent(in) :: loc
real(r8) :: get_location

if ( .not. module_initialized ) call initialize_location

get_location = loc%vloc

end function get_location



function set_location(vert_loc, which_vert)
!----------------------------------------------------------------------------
!
! Given a location type and a double precision value between 0 and 1
! puts this value into the location.

implicit none

type (location_type) :: set_location
real(r8), intent(in) :: vert_loc
integer,  intent(in) :: which_vert

if ( .not. module_initialized ) call initialize_location

set_location%vloc = vert_loc
set_location%which_vert = which_vert

end function set_location


function set_location2(list)
!----------------------------------------------------------------------------
!
! location semi-independent interface routine
! given 2 float numbers, call the underlying set_location routine

implicit none

type (location_type) :: set_location2
real(r8), intent(in) :: list(:)

character(len=129) :: errstring

if ( .not. module_initialized ) call initialize_module

if (size(list) /= 2) then
   write(errstring,*)'requires 2 input values'
   call error_handler(E_ERR, 'set_location2', errstring, source, revision, revdate)
endif

set_location2 = set_location(list(1), nint(list(2)))

end function set_location2



function vert_is_undef(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is undefined, else false

logical :: vert_is_undef
type(location_type), intent(in) :: loc

if ( .not. module_initialized ) call initialize_location

if(loc%which_vert == VERTISUNDEF) then
   vert_is_undef = .true.
else
   vert_is_undef = .false.
endif

end function vert_is_undef



function vert_is_surface(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is surface, else false

logical :: vert_is_surface
type(location_type), intent(in) :: loc

if ( .not. module_initialized ) call initialize_location

if(loc%which_vert == VERTISSURFACE) then
   vert_is_surface = .true.
else
   vert_is_surface = .false.
endif

end function vert_is_surface



function vert_is_pressure(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is pressure, else false

logical :: vert_is_pressure
type(location_type), intent(in) :: loc

if ( .not. module_initialized ) call initialize_location

if(loc%which_vert == VERTISPRESSURE) then
   vert_is_pressure = .true.
else
   vert_is_pressure = .false.
endif

end function vert_is_pressure



function vert_is_height(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is height, else false

logical :: vert_is_height
type(location_type), intent(in) :: loc

if ( .not. module_initialized ) call initialize_location

if(loc%which_vert == VERTISHEIGHT) then
   vert_is_height = .true.
else
   vert_is_height = .false.
endif

end function vert_is_height



function vert_is_level(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is level, else false

logical :: vert_is_level
type(location_type), intent(in) :: loc

if ( .not. module_initialized ) call initialize_location

if(loc%which_vert == VERTISLEVEL) then
   vert_is_level = .true.
else
   vert_is_level = .false.
endif

end function vert_is_level



function loc_eq(loc1,loc2)
!---------------------------------------------------------------------------
!
! interface operator used to compare two locations.
! Returns true only if all components are 'the same' to within machine
! precision.

implicit none

type(location_type), intent(in) :: loc1, loc2
logical :: loc_eq

if ( .not. module_initialized ) call initialize_location

loc_eq = .false.

if ( abs(loc1%vloc - loc2%vloc) > epsilon(loc1%vloc) ) return

loc_eq = .true.

end function loc_eq



function loc_ne(loc1,loc2)
!---------------------------------------------------------------------------
!
! interface operator used to compare two locations.
! Returns true if locations are not identical to machine precision.

implicit none

type(location_type), intent(in) :: loc1, loc2
logical :: loc_ne

if ( .not. module_initialized ) call initialize_location

loc_ne = (.not. loc_eq(loc1,loc2))

end function loc_ne


function set_location_missing()
!----------------------------------------------------------------------------
!
implicit none

type (location_type) :: set_location_missing

if ( .not. module_initialized ) call initialize_location

set_location_missing%vloc = MISSING_R8
set_location_missing%which_vert = missing_i

end function set_location_missing



function query_location(loc,attr) result(fval)
!---------------------------------------------------------------------------
!
! Returns the value of the attribute
!

implicit none

type(location_type),        intent(in) :: loc
character(len=*), optional, intent(in) :: attr
real(r8)                               :: fval

character(len=16) :: attribute

if ( .not. module_initialized ) call initialize_location

attribute = 'which_vert'
if (present(attr)) attribute = attr
selectcase(adjustl(attribute))
 case ('which_vert','WHICH_VERT')
   fval = loc%which_vert
 case ('vloc','VLOC')
   fval = loc%vloc
 case default
   call error_handler(E_ERR, 'query_location; column', &
         'Only vloc and which_vert are legal attributes to request from location', source, revision, revdate)
end select

end function query_location



subroutine write_location(locfile, loc, fform)
!----------------------------------------------------------------------------
!
! Writes a oned location to the file. Implemented as a subroutine but  could
! rewrite as a function with error control info returned. For initial implementation,
! file is just an integer file unit number. Probably want to replace this with file
! as a file_type allowing more flexibility for IO at later point. file_type and 
! associated operations would have to be supported. The mpp_io intefaces are a good
! place to head with this, perhaps, when we need to extend to supporting parallel
! platforms. 

implicit none

integer, intent(in) :: locfile
type(location_type), intent(in) :: loc
character(len = *), intent(in), optional :: fform

character(len = 32) :: fileformat

if ( .not. module_initialized ) call initialize_location

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! For now, output a character tag followed by the r8 value. 

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      write(locfile) loc%vloc, loc%which_vert
   CASE DEFAULT
      write(locfile, '(''loc1d'')' ) 
      write(locfile, *) loc%vloc, loc%which_vert
END SELECT

end subroutine write_location



function read_location(locfile, fform)
!----------------------------------------------------------------------------
!
! Reads a oned location from file that was written by write_location. 
! See write_location for additional discussion.

implicit none

integer, intent(in) :: locfile
type(location_type) :: read_location
character(len = *), intent(in), optional :: fform

character(len=5) :: header
character(len = 32) :: fileformat

if ( .not. module_initialized ) call initialize_location

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      read(locfile) read_location%vloc, read_location%which_vert
   CASE DEFAULT
      read(locfile, '(a5)' ) header
      if(header /= 'loc1d') call error_handler(E_ERR, 'read_location', &
          'Expected location header "loc1d" in input file', source, revision, revdate)
! Now read the location data value
      read(locfile, *) read_location%vloc, read_location%which_vert
END SELECT

end function read_location



subroutine interactive_location(location, set_to_default)
!--------------------------------------------------------------------------
!
! Allows for interactive input of a location.

implicit none

type(location_type), intent(out) :: location
logical, intent(in), optional :: set_to_default

real(r8) :: x

if ( .not. module_initialized ) call initialize_location

! If set_to_default is true, then just zero out and return
if(present(set_to_default)) then
   if(set_to_default) then
      location%vloc = 0.0_r8
      location%which_vert = 0
      return
   endif
endif

write(*, *)'Vertical co-ordinate options'
write(*, *)VERTISUNDEF,' --> vertical coordinate undefined'
write(*, *)VERTISSURFACE,' --> surface'
write(*, *)VERTISLEVEL,' --> model level'
write(*, *)VERTISPRESSURE,' --> pressure'
write(*, *)VERTISHEIGHT,' --> height'

100   read(*, *) location%which_vert
if(location%which_vert == VERTISLEVEL ) then
   write(*, *) 'Vertical co-ordinate model level'
   read(*, *) location%vloc
else if(location%which_vert == VERTISPRESSURE ) then
   write(*, *) 'Vertical co-ordinate Pressure (in hPa)'
   read(*, *) location%vloc
   location%vloc = 100.0_r8 * location%vloc
else if(location%which_vert == VERTISHEIGHT ) then
   write(*, *) 'Vertical co-ordinate height (in gpm)'
   read(*, *) location%vloc
else if(location%which_vert == VERTISSURFACE ) then
   write(*, *) 'Vertical co-ordinate surface pressure (in hPa)'
   read(*, *) location%vloc
   location%vloc = 100.0_r8 * location%vloc
else
   write(*, *) 'Wrong choice of which_vert try again between ',VERTISUNDEF, &
               ' and ',VERTISHEIGHT
   go to 100
end if

! TJH is there the possibility the person can enter an illegal value?
! A nonsensical value?
! i.e. should there be some error checking?
! I removed the 'random' component, and the associated modules.

end subroutine interactive_location



 subroutine alloc_get_close_obs(num, obs, cutoff, obs_box)
!----------------------------------------------------------------------------
!subroutine alloc_get_close_obs(num, obs, cutoff, obs_box)
!
! This does pre-computing for close obs; no function needed in one dimension

implicit none

integer,             intent(in)  :: num
type(location_type), intent(in)  :: obs(num)
real(r8),            intent(in)  :: cutoff
integer,             intent(out) :: obs_box(num)

return

end subroutine alloc_get_close_obs



 subroutine get_close_obs(base_ob, num, obs, cutoff, obs_box, num_close, close_ind, dist)
!----------------------------------------------------------------------------
!subroutine get_close_obs(base_ob, num, obs, cutoff, obs_box, num_close, close_ind, dist)
!
! Default version with no smarts; no need to be smart in 1D

implicit none

integer,             intent(in)  :: base_ob, num
type(location_type), intent(in)  :: obs(num)
real(r8),            intent(in)  :: cutoff
integer,             intent(in)  :: obs_box(num)
integer,             intent(out) :: num_close, close_ind(num)
real(r8),            intent(out) :: dist(num)

integer :: i
real(r8) :: this_dist

! Return list of obs that are within cutoff and their distances
num_close = 0
do i = 1, num
   this_dist = get_dist(obs(base_ob), obs(i))
   if(this_dist <= cutoff) then
      ! Add this ob to the list
      num_close = num_close + 1
      close_ind(num_close) = i
      dist(num_close) = this_dist
   endif
end do

end subroutine get_close_obs

function is_location_in_region(loc, minl, maxl)
!----------------------------------------------------------------------------
!
! Returns true if the first arg is between the other two.

implicit none

logical                          :: is_location_in_region
type(location_type), intent(in)  :: loc, minl, maxl

character(len=129) :: errstring

if ( .not. module_initialized ) call initialize_module

if ((minl%which_vert /= maxl%which_vert) .or. &
    (minl%which_vert /= loc%which_vert)) then
   write(errstring,*)'which_vert (',loc%which_vert,') must be same in all args'
   call error_handler(E_ERR, 'is_location_in_region', errstring, source, revision, revdate)
endif

! assume failure and return as soon as we are confirmed right.
! set to success only at the bottom after all tests have passed.
is_location_in_region = .false.

if ((loc%vloc < minl%vloc) .or. (loc%vloc > maxl%vloc)) return

is_location_in_region = .true.

end function is_location_in_region


!----------------------------------------------------------------------------
! end of location/column/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod
