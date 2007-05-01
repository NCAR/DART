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

! Implements location interfaces for a one dimensional periodic domain. Initial 
! implementation has domain 'longitude' running from 0 to 1. May want to investigate
! allowing an arbitrary real domain size at some point.

use      types_mod, only : r8, MISSING_R8
use  utilities_mod, only : register_module, error_handler, E_ERR
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none
private

public :: location_type, get_location, set_location, set_location_missing, &
          write_location, read_location, interactive_location, query_location, &
          LocationDims, LocationName, LocationLName, get_close_obs, &
          get_close_maxdist_init, get_close_obs_init, get_close_type, &
          operator(==), operator(/=), get_dist, get_close_obs_destroy

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type location_type
   private
   real(r8) :: x
end type location_type

! Needed as stub but not used in this low-order model
type get_close_type
   integer  :: num
   real(r8) :: maxdist
end type get_close_type

type(random_seq_type) :: ran_seq
logical :: ran_seq_init = .false.
logical, save :: module_initialized = .false.

integer,              parameter :: LocationDims = 1
character(len = 129), parameter :: LocationName = "loc1d"
character(len = 129), parameter :: LocationLName = "location on unit circle"

interface operator(==); module procedure loc_eq; end interface
interface operator(/=); module procedure loc_ne; end interface

contains

  subroutine initialize_module
!----------------------------------------------------------------------------
! subroutine initialize_module
!
! pretty simple for this module.

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module



function get_dist(loc1, loc2, kind1, kind2)
!----------------------------------------------------------------------------

implicit none

type(location_type), intent(in) :: loc1, loc2
integer,             intent(in) :: kind1, kind2
real(r8)                        :: get_dist

if ( .not. module_initialized ) call initialize_module

! Reentrant domain, if distance is greater than half wraparound the other way.
get_dist = abs(loc1%x - loc2%x)
if(get_dist > 0.5_r8) get_dist = 1.0_r8 - get_dist

end function get_dist



function loc_eq(loc1,loc2)
!---------------------------------------------------------------------------
!
! interface operator used to compare two locations.
! Returns true only if all components are 'the same' to within machine
! precision.

implicit none

type(location_type), intent(in) :: loc1, loc2
logical :: loc_eq

if ( .not. module_initialized ) call initialize_module

loc_eq = .false.

if ( abs(loc1%x  - loc2%x ) > epsilon(loc1%x ) ) return

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

if ( .not. module_initialized ) call initialize_module

loc_ne = (.not. loc_eq(loc1,loc2))

end function loc_ne



function get_location(loc)
!---------------------------------------------------------------------------
!
! Given a location type, return the x coordinate

implicit none

type(location_type), intent(in) :: loc
real(r8) :: get_location

if ( .not. module_initialized ) call initialize_module

get_location = loc%x

end function get_location



function set_location(x)
!----------------------------------------------------------------------------
!
! Given a location type and a double precision value between 0 and 1
! puts this value into the location.

implicit none

type (location_type) :: set_location
real(r8), intent(in) :: x

if ( .not. module_initialized ) call initialize_module

if(x < 0.0_r8 .or. x > 1.0_r8) call error_handler(E_ERR, 'set_location', &
         'Value of x is out of 0->1 range', source, revision, revdate)

set_location%x = x

end function set_location



function set_location_missing()
!----------------------------------------------------------------------------
!
implicit none

type (location_type) :: set_location_missing

if ( .not. module_initialized ) call initialize_module

set_location_missing%x = MISSING_R8

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

if ( .not. module_initialized ) call initialize_module

attribute = 'x'
if (present(attr)) attribute = attr
selectcase(adjustl(attribute))
 case ('x','X')
   fval = loc%x
 case default
   call error_handler(E_ERR, 'query_location; oned', &
         'Only x is legal attribute to request from location', source, revision, revdate)
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

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! For now, output a character tag followed by the r8 value. 

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      write(locfile) loc%x
   CASE DEFAULT
      write(locfile, '(''loc1d'')' ) 
      write(locfile, *) loc%x
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

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      read(locfile) read_location%x
   CASE DEFAULT
      read(locfile, '(a5)' ) header
      if(header /= 'loc1d') call error_handler(E_ERR, 'read_location', &
          'Expected location header "loc1d" in input file', source, revision, revdate)
! Now read the location data value
      read(locfile, *) read_location%x
END SELECT

end function read_location



subroutine interactive_location(location, set_to_default)
!--------------------------------------------------------------------------
!
! Allows for interactive input of a location. Also gives option of selecting
! a uniformly distributed random location (what the heck).

implicit none

type(location_type), intent(out) :: location
logical, intent(in), optional :: set_to_default

real(r8) :: x

if ( .not. module_initialized ) call initialize_module

! If set_to_default is true, then just zero out and return
if(present(set_to_default)) then
   if(set_to_default) then
      location%x = 0
      return
   endif
endif

write(*, *) 'Input location for this obs: value 0 to 1 or a negative number for '
write(*, *) 'Uniformly distributed random location'
read(*, *) x

do while(x > 1.0_r8)
   write(*, *) 'Input value greater than 1.0 is illegal, please try again'
   read(*, *) x
end do

if(x < 0.0_r8) then

   ! Need to make sure random sequence is initialized

   if(.not. ran_seq_init) then
      call init_random_seq(ran_seq)
      ran_seq_init = .TRUE.
   endif

   ! Uniform location from 0 to 1 for this location type

   location%x = random_uniform(ran_seq)
   write(*, *) 'random location is ', location%x

else
   location%x = x
end if

end subroutine interactive_location

!----------------------------------------------------------------------------

subroutine get_close_obs_init(gc, num, obs)

! Initializes part of get_close accelerator that depends on the particular obs
! Currently not doing much in this oned location 

implicit none

type(get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
type(location_type),  intent(in)    :: obs(num)

! Set the value of num_obs in the structure
gc%num = num

end subroutine get_close_obs_init

!----------------------------------------------------------------------------

subroutine get_close_obs_destroy(gc)

type(get_close_type), intent(inout) :: gc

end subroutine get_close_obs_destroy
!----------------------------------------------------------------------------

subroutine get_close_maxdist_init(gc, maxdist)

implicit none

type(get_close_type), intent(inout) :: gc
real(r8),             intent(in)    :: maxdist

! Set the maximum distance in the structure
gc%maxdist = maxdist

end subroutine get_close_maxdist_init

!----------------------------------------------------------------------------

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs, obs_kind, &
   num_close, close_ind, dist)

! Default version with no smarts; no need to be smart in 1D
! Kinds are available here if one wanted to do more refined distances.

implicit none

type(get_close_type), intent(in)  :: gc
type(location_type),  intent(in)  :: base_obs_loc, obs(:)
integer,              intent(in)  :: base_obs_kind, obs_kind(:)
integer,              intent(out) :: num_close, close_ind(:)
real(r8), optional,   intent(out) :: dist(:)

integer :: i
real(r8) :: this_dist

! Return list of obs that are within maxdist and their distances
num_close = 0
do i = 1, gc%num
   this_dist = get_dist(base_obs_loc, obs(i), base_obs_kind, obs_kind(i))
   if(this_dist <= gc%maxdist) then
      ! Add this ob to the list
      num_close = num_close + 1
      close_ind(num_close) = i
      if (present(dist)) dist(num_close) = this_dist 
   endif
end do

end subroutine get_close_obs

!----------------------------------------------------------------------------
! end of location/oned/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod
