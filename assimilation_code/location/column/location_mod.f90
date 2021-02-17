! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module location_mod

! Implements location interfaces for a one dimensional column domain

use            types_mod, only : i8, r8, MISSING_R8, MISSING_I
use ensemble_manager_mod, only : ensemble_type
use        utilities_mod, only : error_handler, E_ERR, ascii_file_format

implicit none
private

public :: location_type, get_location, set_location, &
          set_location_missing, is_location_in_region, get_maxdist, &
          write_location, read_location, interactive_location, query_location, &
          LocationDims, LocationName, LocationLName, LocationStorageOrder, LocationUnits, &
          get_close_type, get_close_init, get_close_obs, get_close_state, get_close_destroy, &
          operator(==), operator(/=), get_dist, has_vertical_choice, vertical_localization_on, &
          set_vertical, is_vertical, get_vertical_localization_coord, get_close, &
          set_vertical_localization_coord, convert_vertical_obs, convert_vertical_state, &
          VERTISUNDEF, VERTISSURFACE, VERTISLEVEL, VERTISPRESSURE, &
          VERTISHEIGHT, VERTISSCALEHEIGHT

character(len=*), parameter :: source = 'column/location_mod.f90'

! The possible values for the location_type%which_vert component.
! These are intended to be PRIVATE to this module. Do not make public.

integer, parameter :: VERTISUNDEF       = -2  ! has no specific vertical location (undefined)
integer, parameter :: VERTISSURFACE     = -1  ! surface value (value is surface elevation in m)
integer, parameter :: VERTISLEVEL       =  1  ! by level
integer, parameter :: VERTISPRESSURE    =  2  ! by pressure (in pascals)
integer, parameter :: VERTISHEIGHT      =  3  ! by height (in meters)
integer, parameter :: VERTISSCALEHEIGHT =  4  ! by scale height (unitless)

type location_type
   private
   real(r8) :: vloc
   integer  :: which_vert
end type location_type

! Needed as stub but not used in this low-order model
type get_close_type
   private
   integer  :: num
   real(r8) :: maxdist
end type get_close_type
 
integer :: location_vertical_localization_coord = 0

logical, save :: module_initialized = .false.

integer,             parameter :: LocationDims = 1
character(len = 64), parameter :: LocationName = "loc1Dcolumn"
character(len = 64), parameter :: LocationLName = "one-dimensional column"
character(len = 64), parameter :: LocationStorageOrder = "Vertical"
character(len = 64), parameter :: LocationUnits = "which_vert"


character(len = 129) :: errstring

interface operator(==); module procedure loc_eq; end interface
interface operator(/=); module procedure loc_ne; end interface

interface set_location
   module procedure set_location_single
   module procedure set_location_array
end interface set_location

contains

!----------------------------------------------------------------------------

subroutine initialize_module
 
if (module_initialized) return

module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------------

function get_dist(loc1, loc2, type1, kind2)

! Return distance between 2 locations.  The vertical type must
! be the same for this to return a value.

type(location_type), intent(in) :: loc1, loc2
integer, optional,   intent(in) :: type1, kind2
real(r8)                        :: get_dist

if ( .not. module_initialized ) call initialize_module

! if not same, return very long (but positive) distance: -1*missing_r8
if (loc1%which_vert /= loc2%which_vert) then
   get_dist = -1.0 * MISSING_R8
   return
endif

get_dist = abs(loc1%vloc - loc2%vloc)

end function get_dist

!---------------------------------------------------------------------------

function loc_eq(loc1,loc2)
 
! Interface operator used to compare two locations.
! Returns true only if all components are 'the same' to within machine
! precision.  vert type must also match.

type(location_type), intent(in) :: loc1, loc2
logical                         :: loc_eq

if ( .not. module_initialized ) call initialize_module

loc_eq = .false.

if ( loc1%which_vert /= loc2%which_vert ) return
if ( abs(loc1%vloc - loc2%vloc) > epsilon(loc1%vloc) ) return

loc_eq = .true.

end function loc_eq

!---------------------------------------------------------------------------

function loc_ne(loc1,loc2)
 
! Interface operator used to compare two locations.
! Returns true if locations are not identical to machine precision.

type(location_type), intent(in) :: loc1, loc2
logical                         :: loc_ne

if ( .not. module_initialized ) call initialize_module

loc_ne = (.not. loc_eq(loc1,loc2))

end function loc_ne

!---------------------------------------------------------------------------

function get_location(loc)
 
! Given a location type, return the z coordinate

type(location_type), intent(in) :: loc
real(r8)                        :: get_location

if ( .not. module_initialized ) call initialize_module

get_location = loc%vloc

end function get_location

!----------------------------------------------------------------------------

function set_location_single(vert_loc, which_vert)
 
! Store a vertical coordinate and vertical type in a location type.

real(r8), intent(in) :: vert_loc
integer,  intent(in) :: which_vert
type(location_type)  :: set_location_single

if ( .not. module_initialized ) call initialize_module

if (which_vert /= VERTISUNDEF    .and. &
    which_vert /= VERTISSURFACE  .and. &
    which_vert /= VERTISLEVEL    .and. &
    which_vert /= VERTISPRESSURE .and. &
    which_vert /= VERTISHEIGHT) then
   write(errstring,*) 'Illegal value for "which_vert", must be -2, -1, 1, 2, or 3'
   call error_handler(E_ERR, 'set_location', errstring, source)
endif
set_location_single%which_vert = which_vert
set_location_single%vloc = vert_loc

end function set_location_single

!----------------------------------------------------------------------------

function set_location_array(list)
 
! location semi-independent interface routine
! given 2 float numbers, call the underlying set_location routine

real(r8), intent(in) :: list(:)
type(location_type)  :: set_location_array

if ( .not. module_initialized ) call initialize_module

if (size(list) < 2) then
   write(errstring,*) 'requires 2 input values'
   call error_handler(E_ERR, 'set_location', errstring, source)
endif

set_location_array = set_location_single(list(1), nint(list(2)))

end function set_location_array

!----------------------------------------------------------------------------

function set_location_missing()

! Initialize a location type to indicate the contents are unset.

type (location_type) :: set_location_missing

if ( .not. module_initialized ) call initialize_module

set_location_missing%vloc = MISSING_R8
set_location_missing%which_vert = MISSING_I

end function set_location_missing

!---------------------------------------------------------------------------

function query_location(loc,attr)
 
! Returns the value of the attribute
! returns the requested part of the location type

type(location_type),        intent(in) :: loc
character(len=*), optional, intent(in) :: attr
real(r8)                               :: query_location

if ( .not. module_initialized ) call initialize_module

! see the long comment in this routine in the threed_sphere
! module for warnings about compiler bugs before you change
! this code.

query_location = loc%which_vert

if (.not. present(attr)) return

select case(attr)
 case ('which_vert','WHICH_VERT')
   query_location = loc%which_vert
 case ('vloc','VLOC')
   query_location = loc%vloc
 case default
   call error_handler(E_ERR, 'query_location; column', &
         'Only vloc and which_vert are legal attributes to request from location', source)
end select

end function query_location

!----------------------------------------------------------------------------

subroutine write_location(locfile, loc, fform, charstring)
 
! Writes a location to either a file, or if specified, fills in the character
! string instead.

! additional functionality: if optional argument charstring is specified,
! it must be long enough to hold the string, and the location information is
! written into it instead of to a file.  fform must be ascii (which is the
! default if not specified) to use this option.

integer, intent(in)                        :: locfile
type(location_type), intent(in)            :: loc
character(len = *),  intent(in),  optional :: fform
character(len = *),  intent(out), optional :: charstring

integer             :: charlength
logical             :: writebuf

! 10 format(1x,(f22.14,1x),i4)  ! old
10 format(1X,(G25.16,1X),I2) 

if ( .not. module_initialized ) call initialize_module

! writing to a file (normal use) or to a character buffer?
writebuf = present(charstring)

! output file; test for ascii or binary, write what's asked, and return
if (.not. writebuf) then
   if (ascii_file_format(fform)) then
      write(locfile, '(''loc1c'')' ) 
      write(locfile, 10) loc%vloc, loc%which_vert
   else
      write(locfile) loc%vloc, loc%which_vert
   endif
   return
endif

! you only get here if you're writing to a buffer and not
! to a file, and you can't have binary format set.
if (.not. ascii_file_format(fform)) then
   call error_handler(E_ERR, 'write_location', &
      'Cannot use string buffer with binary format', source)
endif

! format the location to be more human-friendly; meaning
! kilometers for height, hectopascals instead of pascals 
! for pressure, etc.


! this must be the sum of the longest of the formats below.
charlength = 32

if (len(charstring) < charlength) then
   write(errstring, *) 'charstring buffer must be at least ', charlength, ' chars long'
   call error_handler(E_ERR, 'write_location', errstring, source)
endif

select case  (loc%which_vert)
   case (VERTISUNDEF)
      write(charstring, '(A,A)')       'Vert: ', '               Undefined'
   case (VERTISSURFACE)
      write(charstring, '(A,F12.5,A)') 'Vert: ', loc%vloc, '   surface (m)'
   case (VERTISLEVEL)
      write(charstring, '(A,F5.0,A)')  'Vert: ', loc%vloc, '          level'
   case (VERTISPRESSURE)
      write(charstring, '(A,F12.7,A)') 'Vert: ', loc%vloc / 100.0_r8, '   hPa'
   case (VERTISHEIGHT)
      write(charstring, '(A,F12.7,A)') 'Vert: ', loc%vloc / 1000.0_r8, '   km'
   case (VERTISSCALEHEIGHT)
      write(charstring, '(A,F13.7,A)') 'Vert: ', loc%vloc, ' scale ht'
   case default
      write(errstring, *) 'unrecognized key for vertical type: ', loc%which_vert
      call error_handler(E_ERR, 'write_location', errstring, source)
end select


end subroutine write_location

!----------------------------------------------------------------------------

function read_location(locfile, fform)
 
! Reads a location from file that was written by write_location. 
! See write_location for additional discussion.

integer, intent(in)                      :: locfile
type(location_type)                      :: read_location
character(len = *), intent(in), optional :: fform

character(len=5) :: header

if ( .not. module_initialized ) call initialize_module

if (ascii_file_format(fform)) then
   read(locfile, '(a5)' ) header
   if(header /= 'loc1c') then
      write(errstring,*)'Expected location header "loc1c" in input file, got ', header 
      call error_handler(E_ERR, 'read_location', errstring, source)
   endif
   ! Now read the location data value
   read(locfile, *) read_location%vloc, read_location%which_vert
else
   read(locfile) read_location%vloc, read_location%which_vert
endif

end function read_location

!--------------------------------------------------------------------------

subroutine interactive_location(location, set_to_default)
 
! Allows for interactive input of a location.

type(location_type), intent(out) :: location
logical, intent(in), optional    :: set_to_default

if ( .not. module_initialized ) call initialize_module

! If set_to_default is true, then just zero out and return
if(present(set_to_default)) then
   if(set_to_default) then
      location%vloc = 0.0_r8
      location%which_vert = 0
      return
   endif
endif

write(*, *)'Vertical coordinate options'
write(*, *)VERTISUNDEF,' --> vertical coordinate undefined'
write(*, *)VERTISSURFACE,' --> surface'
write(*, *)VERTISLEVEL,' --> model level'
write(*, *)VERTISPRESSURE,' --> pressure'
write(*, *)VERTISHEIGHT,' --> height'
write(*, *)VERTISSCALEHEIGHT,' --> scale height'

100   read(*, *) location%which_vert
if(location%which_vert == VERTISLEVEL ) then
   write(*, *) 'Vertical coordinate model level'
   read(*, *) location%vloc
else if(location%which_vert == VERTISPRESSURE ) then
   write(*, *) 'Vertical coordinate Pressure (in hPa)'
   read(*, *) location%vloc
   location%vloc = 100.0_r8 * location%vloc
else if(location%which_vert == VERTISHEIGHT ) then
   write(*, *) 'Vertical coordinate height (in meters)'
   read(*, *) location%vloc
else if(location%which_vert == VERTISSCALEHEIGHT ) then
   write(*, *) 'Vertical coordinate scale height (-ln(p/ps))'
   read(*, *) location%vloc
else if(location%which_vert == VERTISSURFACE ) then
   write(*, *) 'Vertical coordinate surface height'
   read(*, *) location%vloc
else if(location%which_vert == VERTISUNDEF ) then
   location%vloc = 0.0_r8
else
   write(*, *) 'Wrong choice of which_vert try again between ',VERTISUNDEF, &
               ' and ',VERTISHEIGHT
   go to 100
end if

! FIXME: the other location mods have a random location option.

end subroutine interactive_location

!----------------------------------------------------------------------------

subroutine get_close_destroy(gc)

type(get_close_type), intent(inout) :: gc

end subroutine get_close_destroy

!----------------------------------------------------------------------------

subroutine get_close_init(gc, num, maxdist, locs, maxdist_list)

type(get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
real(r8),             intent(in)    :: maxdist
type(location_type),  intent(in)    :: locs(:)
real(r8), intent(in), optional      :: maxdist_list(:)

! Set the count and maximum distance in the structure
gc%num = num
gc%maxdist = maxdist

end subroutine get_close_init

!----------------------------------------------------------------------------

subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

! The specific type of the base observation, plus the generic kinds list
! for either the state or obs lists are available if a more sophisticated
! distance computation is needed.

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(inout) :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:), loc_types(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ens_handle

call get_close(gc, base_loc, base_type, locs, loc_qtys, &
               num_close, close_ind, dist, ens_handle)

end subroutine get_close_obs

!----------------------------------------------------------------------------

subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)

! The specific type of the base observation, plus the generic kinds list
! for either the state or obs lists are available if a more sophisticated
! distance computation is needed.

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:)
integer(i8),                   intent(in)  :: loc_indx(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ens_handle

call get_close(gc, base_loc, base_type, locs, loc_qtys, &
               num_close, close_ind, dist, ens_handle)

end subroutine get_close_state

!----------------------------------------------------------------------------

subroutine get_close(gc, base_loc, base_type, locs, loc_qtys, &
                     num_close, close_ind, dist, ens_handle)

! The specific type of the base observation, plus the generic kinds list
! for either the state or obs lists are available if a more sophisticated
! distance computation is needed.

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ens_handle

! Default version with no smarts; no need to be smart in 1D
! Quantities are available here if one wanted to do more refined distances.

integer :: i
real(r8) :: this_dist

! Return list of locs that are within maxdist and their distances
num_close = 0
do i = 1, gc%num
   this_dist = get_dist(base_loc, locs(i), base_type, loc_qtys(i))
   if(this_dist <= gc%maxdist) then
      ! Add this ob to the list
      num_close = num_close + 1
      close_ind(num_close) = i
      if (present(dist)) dist(num_close) = this_dist 
   endif
end do

end subroutine get_close

!---------------------------------------------------------------------------

function get_maxdist(gc, obs_type)
type(get_close_type), intent(in) :: gc
integer, optional,    intent(in) :: obs_type
real(r8) :: get_maxdist

get_maxdist = gc%maxdist

end function get_maxdist

!----------------------------------------------------------------------------

function is_location_in_region(loc, minl, maxl)
 
! Returns true if the given location is between the other two.

logical                          :: is_location_in_region
type(location_type), intent(in)  :: loc, minl, maxl

if ( .not. module_initialized ) call initialize_module

if ((minl%which_vert /= maxl%which_vert) .or. &
    (minl%which_vert /= loc%which_vert)) then
   write(errstring,*)'which_vert (',loc%which_vert,') must be same in all args'
   call error_handler(E_ERR, 'is_location_in_region', errstring, source)
endif

! assume failure and return as soon as we are confirmed right.
! set to success only at the bottom after all tests have passed.
is_location_in_region = .false.

if ((loc%vloc < minl%vloc) .or. (loc%vloc > maxl%vloc)) return
 
is_location_in_region = .true.

end function is_location_in_region

!----------------------------------------------------------------------------

function has_vertical_choice()

logical :: has_vertical_choice

if ( .not. module_initialized ) call initialize_module

has_vertical_choice = .true.

end function has_vertical_choice

!----------------------------------------------------------------------------

function get_vertical_localization_coord()

integer :: get_vertical_localization_coord

if ( .not. module_initialized ) call initialize_module

get_vertical_localization_coord = location_vertical_localization_coord

end function get_vertical_localization_coord

!----------------------------------------------------------------------------

subroutine set_vertical_localization_coord(which_vert)

integer, intent(in) :: which_vert

if ( .not. module_initialized ) call initialize_module

location_vertical_localization_coord = which_vert

end subroutine set_vertical_localization_coord

!---------------------------------------------------------------------------

function vertical_localization_on()

logical :: vertical_localization_on

if ( .not. module_initialized ) call initialize_module

!>@todo FIXME
!> Always returns false since this type of location doesn't support
!> vertical localization (no namelist choice).  but it could since
!> it makes sense for this type.
vertical_localization_on = .false.

end function vertical_localization_on

!----------------------------------------------------------------------------
!> use a string so caller doesn't have to have access to VERTISxxx values

function is_vertical(loc, which_vert)

logical                          :: is_vertical
type(location_type), intent(in)  :: loc
character(len=*),    intent(in)  :: which_vert

select case  (which_vert)
   case ("UNDEFINED")
      is_vertical = (VERTISUNDEF == loc%which_vert)
   case ("SURFACE")
      is_vertical = (VERTISSURFACE == loc%which_vert)
   case ("LEVEL")
      is_vertical = (VERTISLEVEL == loc%which_vert) 
   case ("PRESSURE")
      is_vertical = (VERTISPRESSURE == loc%which_vert)
   case ("HEIGHT") 
      is_vertical = (VERTISHEIGHT == loc%which_vert)
   case ("SCALE_HEIGHT")
      is_vertical = (VERTISSCALEHEIGHT == loc%which_vert)
   case default
      write(errstring, *) 'unrecognized key for vertical type: ', which_vert
      call error_handler(E_ERR, 'is_vertical', errstring, source)
end select

end function is_vertical
   
!--------------------------------------------------------------------

subroutine set_vertical(loc, vloc, which_vert)

type(location_type), intent(inout) :: loc
real(r8), optional,  intent(in)    :: vloc
integer,  optional,  intent(in)    :: which_vert

if (present(vloc)) loc%vloc = vloc
if (present(which_vert)) loc%which_vert = which_vert

end subroutine set_vertical

!--------------------------------------------------------------------

subroutine convert_vertical_obs(ens_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, status)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:), loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: status(:)

status(:) = 0

end subroutine convert_vertical_obs

!--------------------------------------------------------------------

subroutine convert_vertical_state(ens_handle, num, locs, loc_qtys, loc_indx, &
                                  which_vert, istatus)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer(i8),         intent(in)    :: loc_indx(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus

istatus = 0

end subroutine convert_vertical_state


!----------------------------------------------------------------------------
! end of location/column/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod

