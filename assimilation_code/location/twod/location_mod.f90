! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module location_mod

! Implements location interfaces for a two dimensional cyclic region.
! The internal representation of the location is currently implemented
! as (x, y) from 0.0 to 1.0 in both dimensions.

use      types_mod, only : r8, MISSING_R8, i8
use  utilities_mod, only : error_handler, E_ERR, ascii_file_format
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform
use default_location_mod, only : has_vertical_choice, vertical_localization_on, &
                                 get_vertical_localization_coord, &
                                 set_vertical_localization_coord
use ensemble_manager_mod, only : ensemble_type

implicit none
private

public :: location_type, get_location, set_location, &
          set_location_missing, is_location_in_region, get_maxdist, &
          write_location, read_location, interactive_location, query_location, &
          LocationDims, LocationName, LocationLName, LocationStorageOrder, LocationUnits, &
          get_close_obs, get_close_maxdist_init, get_close_init, get_close_type, &
          operator(==), operator(/=), get_dist, get_close_destroy, &
          is_vertical, set_vertical, vertical_localization_on, &
          set_vert, get_vert, set_which_vert, has_vertical_choice, &
          get_close_state, &
          convert_vertical_obs, convert_vertical_state, get_vertical_localization_coord, &
          set_vertical_localization_coord

character(len=*), parameter :: source = 'twod/location_mod.f90'

type location_type
   private
   real(r8) :: x, y
end type location_type

! Needed as stub but not used in this low-order model
type get_close_type
   private
   integer  :: num
   real(r8) :: maxdist
end type get_close_type

type(random_seq_type) :: ran_seq
logical :: ran_seq_init = .false.
logical, save :: module_initialized = .false.

integer,              parameter :: LocationDims = 2
character(len = 129), parameter :: LocationName = "loc2D"
character(len = 129), parameter :: LocationLName = "twod cyclic locations: x, y"
character(len = 129), parameter :: LocationStorageOrder = "x, y"
character(len = 129), parameter :: LocationUnits = " "

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

! Return the distance between 2 locations.  Since this is a periodic
! domain, the shortest distance may wrap around.

type(location_type), intent(in) :: loc1, loc2
integer, optional,   intent(in) :: type1, kind2
real(r8)                        :: get_dist

real(r8) :: x_dif, y_dif

if ( .not. module_initialized ) call initialize_module

! Periodic domain, if distance is greater than half wraparound the other way.
x_dif = abs(loc1%x - loc2%x)
if (x_dif > 0.5_r8) x_dif = 1.0_r8 - x_dif
y_dif = abs(loc1%y - loc2%y)
if (y_dif > 0.5_r8) y_dif = 1.0_r8 - y_dif

get_dist = sqrt ( x_dif * x_dif + y_dif * y_dif)

end function get_dist

!---------------------------------------------------------------------------

function loc_eq(loc1,loc2)
 
! interface operator used to compare two locations.
! Returns true only if all components are 'the same' to within machine
! precision.

type(location_type), intent(in) :: loc1, loc2
logical                         :: loc_eq

if ( .not. module_initialized ) call initialize_module

loc_eq = .false.

if ( abs(loc1%x  - loc2%x ) > epsilon(loc1%x ) ) return
if ( abs(loc1%y  - loc2%y ) > epsilon(loc1%y ) ) return

loc_eq = .true.

end function loc_eq

!---------------------------------------------------------------------------

function loc_ne(loc1,loc2)
 
! interface operator used to compare two locations.
! Returns true if locations are not identical to machine precision.

type(location_type), intent(in) :: loc1, loc2
logical                         :: loc_ne

if ( .not. module_initialized ) call initialize_module

loc_ne = (.not. loc_eq(loc1,loc2))

end function loc_ne

!---------------------------------------------------------------------------

function get_location(loc)
 
! Given a location type, return the x, y

type(location_type), intent(in) :: loc
real(r8), dimension(2) :: get_location

if ( .not. module_initialized ) call initialize_module

get_location(1) = loc%x 
get_location(2) = loc%y

end function get_location

!----------------------------------------------------------------------------

function set_location_single(x, y)
 
! Given an x, y pair, put the values into the location.

real(r8), intent(in) :: x, y
type (location_type) :: set_location_single

if ( .not. module_initialized ) call initialize_module

if(x < 0.0_r8 .or. x > 1.0_r8) then
   write(errstring,*)'x (',x,') is not within range [0,1]'
   call error_handler(E_ERR, 'set_location', errstring, source)
endif

if(y < 0.0_r8 .or. y > 1.0_r8) then
   write(errstring,*)'y (',y,') is not within range [0,1]'
   call error_handler(E_ERR, 'set_location', errstring, source)
endif

set_location_single%x = x
set_location_single%y = y

end function set_location_single

!----------------------------------------------------------------------------

function set_location_array(list)
 
! location semi-independent interface routine
! given 2 float numbers, call the underlying set_location routine

real(r8), intent(in) :: list(:)
type (location_type) :: set_location_array

if ( .not. module_initialized ) call initialize_module

if (size(list) < 2) then
   write(errstring,*)'requires 2 input values'
   call error_handler(E_ERR, 'set_location', errstring, source)
endif

set_location_array = set_location_single(list(1), list(2))

end function set_location_array

!----------------------------------------------------------------------------

function set_location_missing()

! fill in the contents to a known value.

type (location_type) :: set_location_missing

if ( .not. module_initialized ) call initialize_module

set_location_missing%x = MISSING_R8
set_location_missing%y = MISSING_R8

end function set_location_missing

!---------------------------------------------------------------------------

function query_location(loc,attr)
 
! Returns the value of the attribute
!

type(location_type),        intent(in) :: loc
character(len=*), optional, intent(in) :: attr
real(r8)                               :: query_location

if ( .not. module_initialized ) call initialize_module

! see the long comment in this routine in the threed_sphere
! module for warnings about compiler bugs before you change
! this code.

query_location = loc%x

if (.not. present(attr)) return

select case(attr)
 case ('x','X')
   query_location = loc%x
 case ('y','Y')
   query_location = loc%y
 case default
   call error_handler(E_ERR, 'query_location; twod', &
         'Only x or y are legal attributes to request from location', source)
end select

end function query_location

!----------------------------------------------------------------------------

subroutine write_location(locfile, loc, fform, charstring)
 
! Writes a 2D location to the file.
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

! 10 format(1x,2(f22.14,1x))   ! old
10 format(1X,2(F20.16,1X)) 

if ( .not. module_initialized ) call initialize_module

! writing to a file (normal use) or to a character buffer?
writebuf = present(charstring)

! output file; test for ascii or binary, write what's asked, and return
if (.not. writebuf) then
   if (ascii_file_format(fform)) then
      write(locfile, '(''loc2D'')' ) 
      write(locfile, 10) loc%x, loc%y
   else
      write(locfile) loc%x, loc%y
   endif
   return
endif

! you only get here if you're writing to a buffer and not
! to a file, and you can't have binary format set.
if (.not. ascii_file_format(fform)) then
   call error_handler(E_ERR, 'write_location', &
      'Cannot use string buffer with binary format', source)
endif

! format the location to be more human-friendly; which in
! this case doesn't change the value.

! this must be the sum of the formats below.
charlength = 25

if (len(charstring) < charlength) then
   write(errstring, *) 'charstring buffer must be at least ', charlength, ' chars long'
   call error_handler(E_ERR, 'write_location', errstring, source)
endif

write(charstring, '(A,F9.7,2X,F9.7)') 'X/Y: ',  loc%x, loc%y


end subroutine write_location

!----------------------------------------------------------------------------

function read_location(locfile, fform)
 
! Reads a 2D location from locfile that was written by write_location. 
! See write_location for additional discussion.

integer, intent(in)                      :: locfile
type(location_type)                      :: read_location
character(len = *), intent(in), optional :: fform

character(len=5) :: header

if ( .not. module_initialized ) call initialize_module

if (ascii_file_format(fform)) then
   read(locfile, '(a5)' ) header
   if(header /= 'loc2D') then
      write(errstring,*)'Expected location header "loc2D" in input file, got ', header 
      call error_handler(E_ERR, 'read_location', errstring, source)
   endif
   ! Now read the location data value
   read(locfile, *) read_location%x, read_location%y
else
   read(locfile) read_location%x, read_location%y
endif

end function read_location

!--------------------------------------------------------------------------

subroutine interactive_location(location, set_to_default)
 
! Allows for interactive input of a location. Also gives option of selecting
! a uniformly distributed random location.

type(location_type), intent(out) :: location
logical, intent(in), optional    :: set_to_default

real(r8) :: x, y

if ( .not. module_initialized ) call initialize_module

! If set_to_default is true, then just zero out and return
if(present(set_to_default)) then
   if(set_to_default) then
      location%x = 0.0
      location%y = 0.0
      return
   endif
endif

write(*, *) 'Input X location for this obs: value 0 to 1 or a negative number for '
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
   write(*, *) 'random X location is ', location%x

else
   location%x = x
endif

write(*, *) 'Input Y location for this obs: value 0 to 1 or a negative number for '
write(*, *) 'Uniformly distributed random location'
read(*, *) y

do while(y > 1.0_r8)
   write(*, *) 'Input value greater than 1.0 is illegal, please try again'
   read(*, *) y
end do

if(y < 0.0_r8) then

   ! Need to make sure random sequence is initialized

   if(.not. ran_seq_init) then
      call init_random_seq(ran_seq)
      ran_seq_init = .TRUE.
   endif

   ! Uniform location from 0 to 1 for this location type

   location%y = random_uniform(ran_seq)
   write(*, *) 'random Y location is ', location%y

else
   location%y = y
endif

end subroutine interactive_location

!----------------------------------------------------------------------------
! Initializes get_close accelerator - unused in this location module

subroutine get_close_init(gc, num, maxdist, locs, maxdist_list)

type(get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
real(r8),             intent(in)    :: maxdist
type(location_type),  intent(in)    :: locs(:)
real(r8), intent(in), optional      :: maxdist_list(:)

! Set the maximum localization distance
gc%maxdist = maxdist

! Save the value of num_locs
gc%num = num

if (present(maxdist_list)) then
   write(errstring,*)'twod locations does not support different cutoff distances by type'
   call error_handler(E_ERR, 'get_close_init', errstring, source)
endif

end subroutine get_close_init

!----------------------------------------------------------------------------

subroutine get_close_destroy(gc)

type(get_close_type), intent(inout) :: gc

end subroutine get_close_destroy

!----------------------------------------------------------------------------

subroutine get_close_maxdist_init(gc, maxdist, maxdist_list)

type(get_close_type), intent(inout) :: gc
real(r8),             intent(in)    :: maxdist
real(r8), intent(in), optional      :: maxdist_list(:)

! Set the maximum distance in the structure
gc%maxdist = maxdist

end subroutine get_close_maxdist_init

!----------------------------------------------------------------------------
subroutine get_close(gc, base_obs_loc, base_obs_type, obs, obs_kind, &
   num_close, close_ind, dist)

! Default version with no smarts; no need to be smart in 1D
! Kinds are available here if one wanted to do more refined distances.

type(get_close_type), intent(in)  :: gc
type(location_type),  intent(in)  :: base_obs_loc, obs(:)
integer,              intent(in)  :: base_obs_type, obs_kind(:)
integer,              intent(out) :: num_close, close_ind(:)
real(r8), optional,   intent(out) :: dist(:)

integer :: i
real(r8) :: this_dist

! the list of locations in the obs() argument must be the same
! as the list of locations passed into get_close_obs_init(), so
! gc%num and size(obs) better be the same.   if the list changes,
! you have to destroy the old gc and init a new one.
if (size(obs) /= gc%num) then
   write(errstring,*)'obs() array must match one passed to get_close_obs_init()'
   call error_handler(E_ERR, 'get_close_obs', errstring, source)
endif

! Return list of obs that are within maxdist and their distances
num_close = 0
do i = 1, gc%num
   this_dist = get_dist(base_obs_loc, obs(i), base_obs_type, obs_kind(i))
   if(this_dist <= gc%maxdist) then
      ! Add this ob to the list
      num_close = num_close + 1
      close_ind(num_close) = i
      if (present(dist)) dist(num_close) = this_dist 
   endif
end do

end subroutine get_close

!---------------------------------------------------------------------------
subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)


type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(in)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:), loc_types(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ens_handle

call get_close(gc, base_loc, base_type, locs, loc_qtys, &
               num_close, close_ind, dist)

end subroutine get_close_obs

!---------------------------------------------------------------------------
subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)


type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(in)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:)
integer(i8),                   intent(in)  :: loc_indx(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ens_handle

call get_close(gc, base_loc, base_type, locs, loc_qtys, &
               num_close, close_ind, dist)

end subroutine get_close_state

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

! assume failure and return as soon as we are confirmed right.
! set to success only at the bottom after all tests have passed.
is_location_in_region = .false.

! FIXME: this is a doubly cyclic domain.  check if min
! limit > max; if so, then wrap around.
!if (minl%x <= maxl%x) .and.  ...
if ((loc%x < minl%x) .or. (loc%x > maxl%x)) return
if ((loc%y < minl%y) .or. (loc%y > maxl%y)) return
 
is_location_in_region = .true.

end function is_location_in_region

!----------------------------------------------------------------------------
! stubs - always say no, but allow this code to be compiled with
!         common code that sometimes needs vertical info.
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
! stubs - here only because they have a location type as one of the arguments
!----------------------------------------------------------------------------

function is_vertical(loc, which_vert)

logical                          :: is_vertical
type(location_type), intent(in)  :: loc
character(len=*),    intent(in)  :: which_vert

is_vertical = .false.

end function is_vertical

!--------------------------------------------------------------------

subroutine set_vertical(loc, vloc, which_vert)

type(location_type), intent(inout) :: loc
real(r8), optional,  intent(in)    :: vloc
integer,  optional,  intent(in)    :: which_vert


end subroutine set_vertical

!--------------------------------------------------------------------

function has_vertical_localization()
 
! Always returns false since this type of location doesn't support
! vertical localization.

logical :: has_vertical_localization

if ( .not. module_initialized ) call initialize_module

has_vertical_localization = .false.

end function has_vertical_localization

!--------------------------------------------------------------------
!> dummy routine for models that don't have a vertical location
function get_vert(loc)

type(location_type), intent(in) :: loc
real(r8) :: get_vert

get_vert = 1 ! any old value

end function get_vert

!--------------------------------------------------------------------
!> dummy routine for models that don't have a vertical location
subroutine set_vert(loc, vloc)

type(location_type), intent(inout) :: loc
real(r8), intent(in) :: vloc


end subroutine set_vert

!----------------------------------------------------------------------------
!> set the which vert
subroutine set_which_vert(loc, which_vert)

type(location_type), intent(inout) :: loc
integer,                intent(in) :: which_vert !< vertical coordinate type


end subroutine set_which_vert

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
! end of location/twod/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod

