! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module location_mod

! Implements location interfaces for a three dimensional cyclic region.
! The internal representation of the location is currently implemented
! as (x, y) from 0.0 to 1.0 in both dimensions.
!
! If you are looking for a geophysical locations module, look at either
! the threed_sphere or threed_cartesian versions of this file.

use            types_mod, only : r8, MISSING_R8, i8
use        utilities_mod, only : error_handler, E_ERR, ascii_file_format
use       random_seq_mod, only : random_seq_type, init_random_seq, random_uniform
use ensemble_manager_mod, only : ensemble_type
use default_location_mod, only : has_vertical_choice, vertical_localization_on, &
                                 get_vertical_localization_coord, &
                                 set_vertical_localization_coord

implicit none
private

public :: location_type, get_location, set_location, &
          set_location_missing, is_location_in_region, get_maxdist, &
          write_location, read_location, interactive_location, query_location, &
          LocationDims, LocationName, LocationLName, LocationStorageOrder, LocationUnits, &
          get_close_type, get_close_init, get_close_obs, get_close_state, get_close_destroy, &
          operator(==), operator(/=), get_dist, has_vertical_choice, vertical_localization_on, &
          set_vertical, is_vertical, get_vertical_localization_coord, &
          set_vertical_localization_coord, convert_vertical_obs, convert_vertical_state

character(len=*), parameter :: source = 'threed/location_mod.f90'

type location_type
   private
   real(r8) :: x, y, z 
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

integer,              parameter :: LocationDims = 3
character(len = 129), parameter :: LocationName = "loc3D"
character(len = 129), parameter :: LocationLName = "threed cyclic locations: x, y, z"
character(len = 129), parameter :: LocationStorageOrder = "X Y Z"
character(len = 129), parameter :: LocationUnits = "none none none"


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

real(r8) :: x_dif, y_dif, z_dif

if ( .not. module_initialized ) call initialize_module

! Periodic domain, if distance is greater than half wraparound the other way.
x_dif = abs(loc1%x - loc2%x)
if (x_dif > 0.5_r8) x_dif = 1.0_r8 - x_dif
y_dif = abs(loc1%y - loc2%y)
if (y_dif > 0.5_r8) y_dif = 1.0_r8 - y_dif
z_dif = abs(loc1%z - loc2%z)
if (z_dif > 0.5_r8) z_dif = 1.0_r8 - z_dif

get_dist = sqrt ( x_dif * x_dif + y_dif * y_dif + z_dif * z_dif)

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
if ( abs(loc1%z  - loc2%z ) > epsilon(loc1%z ) ) return

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
 
! Given a location type, return the x, y, z

type(location_type), intent(in) :: loc
real(r8), dimension(3) :: get_location

if ( .not. module_initialized ) call initialize_module

get_location(1) = loc%x 
get_location(2) = loc%y
get_location(3) = loc%z

end function get_location

!----------------------------------------------------------------------------

function set_location_single(x, y, z)
 
! Given an x, y, z triplet, put the values into the location.

real(r8), intent(in) :: x, y, z 
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

if(z < 0.0_r8 .or. z > 1.0_r8) then
   write(errstring,*)'z (',z,') is not within range [0,1]'
   call error_handler(E_ERR, 'set_location', errstring, source)
endif

set_location_single%x = x
set_location_single%y = y
set_location_single%z = z

end function set_location_single

!----------------------------------------------------------------------------

function set_location_array(list)
 
! location semi-independent interface routine
! given 3 float numbers, call the underlying set_location routine

real(r8), intent(in) :: list(:)
type (location_type) :: set_location_array

if ( .not. module_initialized ) call initialize_module

if (size(list) < 3) then
   write(errstring,*)'requires 3 input values'
   call error_handler(E_ERR, 'set_location', errstring, source)
endif

set_location_array = set_location_single(list(1), list(2), list(3))

end function set_location_array

!----------------------------------------------------------------------------

function set_location_missing()

! fill in the contents to a known value.

type (location_type) :: set_location_missing

if ( .not. module_initialized ) call initialize_module

set_location_missing%x = MISSING_R8
set_location_missing%y = MISSING_R8
set_location_missing%z = MISSING_R8

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
 case ('z','Z')
   query_location = loc%z
 case default
   call error_handler(E_ERR, 'query_location; threed', &
         'Only x, y, or z are legal attributes to request from location', source)
end select

end function query_location

!----------------------------------------------------------------------------

subroutine write_location(locfile, loc, fform, charstring)
 
! Writes a 3D location to the file.
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
10 format(1X,3(F20.16,1X)) 

if ( .not. module_initialized ) call initialize_module

! writing to a file (normal use) or to a character buffer?
writebuf = present(charstring)

! output file; test for ascii or binary, write what's asked, and return
if (.not. writebuf) then
   if (ascii_file_format(fform)) then
      write(locfile, '(''loc3D'')' ) 
      write(locfile, 10) loc%x, loc%y, loc%z
   else
      write(locfile) loc%x, loc%y, loc%z
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
charlength = 38

if (len(charstring) < charlength) then
   write(errstring, *) 'charstring buffer must be at least ', charlength, ' chars long'
   call error_handler(E_ERR, 'write_location', errstring, source)
endif

write(charstring, '(A,F9.7,2(2X,F9.7))') 'X/Y/Z: ',  loc%x, loc%y, loc%z


end subroutine write_location

!----------------------------------------------------------------------------

function read_location(locfile, fform)
 
! Reads a 3D location from locfile that was written by write_location. 
! See write_location for additional discussion.

integer, intent(in)                      :: locfile
type(location_type)                      :: read_location
character(len = *), intent(in), optional :: fform

character(len=5) :: header

if ( .not. module_initialized ) call initialize_module

if (ascii_file_format(fform)) then
   read(locfile, '(a5)' ) header
   if(header /= 'loc3D') then
      write(errstring,*)'Expected location header "loc3D" in input file, got ', header 
      call error_handler(E_ERR, 'read_location', errstring, source)
   endif
   ! Now read the location data value
   read(locfile, *) read_location%x, read_location%y, read_location%z
else
   read(locfile) read_location%x, read_location%y, read_location%z
endif

end function read_location

!--------------------------------------------------------------------------

subroutine interactive_location(location, set_to_default)
 
! Allows for interactive input of a location. Also gives option of selecting
! a uniformly distributed random location.

type(location_type), intent(out) :: location
logical, intent(in), optional    :: set_to_default

real(r8) :: v(3)
character(len=1) :: l(3)
integer :: i

if ( .not. module_initialized ) call initialize_module

! If set_to_default is true, then just zero out and return
if(present(set_to_default)) then
   if(set_to_default) then
      location%x = 0.0
      location%y = 0.0
      location%z = 0.0
      return
   endif
endif

l(1) = 'X'
l(2) = 'Y'
l(3) = 'Z'

do i=1, 3
   write(*, *) 'Input ', l(i), ' location for this obs: value 0 to 1 or a negative number for '
   write(*, *) 'Uniformly distributed random location'
   read(*, *) v(i)
   
   do while(v(i) > 1.0_r8)
      write(*, *) 'Input value greater than 1.0 is illegal, please try again'
      read(*, *) v(i)
   end do
   
   if(v(i) < 0.0_r8) then
   
      ! Need to make sure random sequence is initialized
   
      if(.not. ran_seq_init) then
         call init_random_seq(ran_seq)
         ran_seq_init = .TRUE.
      endif
   
      ! Uniform location from 0 to 1 for this location type
   
      v(i) = random_uniform(ran_seq)
      write(*, *) 'random ',l(i),' location is ', v(i)
   
   endif
enddo

location%x = v(1)
location%y = v(2)
location%z = v(3)

end subroutine interactive_location

!----------------------------------------------------------------------------

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
   write(errstring,*)'threed locations does not support different cutoff distances by type'
   call error_handler(E_ERR, 'get_close_init', errstring, source)
endif

end subroutine get_close_init

!----------------------------------------------------------------------------

subroutine get_close_destroy(gc)

type(get_close_type), intent(inout) :: gc

end subroutine get_close_destroy

!----------------------------------------------------------------------------

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
               num_close, close_ind, dist, ens_handle)

end subroutine get_close_state

!----------------------------------------------------------------------------
subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(in)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:), loc_types(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ens_handle

call get_close(gc, base_loc, base_type, locs, loc_qtys, &
               num_close, close_ind, dist, ens_handle)

end subroutine get_close_obs


!----------------------------------------------------------------------------
subroutine get_close(gc, base_loc, base_type, locs, loc_qtys, &
                     num_close, close_ind, dist, ensemble_handle)

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(in)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ensemble_handle

integer :: i
real(r8) :: this_dist

! the list of locations in the locs() argument must be the same
! as the list of locations passed into get_close_init(), so
! gc%num and size(locs) better be the same.   if the list changes,
! you have to destroy the old gc and init a new one.
if (size(locs) /= gc%num) then
   write(errstring,*)'locs() array must match one passed to get_close_init()'
   call error_handler(E_ERR, 'get_close', errstring, source)
endif

! Return list of obs that are within maxdist and their distances
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

! assume failure and return as soon as we are confirmed right.
! set to success only at the bottom after all tests have passed.
is_location_in_region = .false.

! FIXME: this is a triply cyclic domain.  check if min
! limit > max; if so, then wrap around.
!if (minl%x <= maxl%x) .and.  ...
if ((loc%x < minl%x) .or. (loc%x > maxl%x)) return
if ((loc%y < minl%y) .or. (loc%y > maxl%y)) return
if ((loc%z < minl%z) .or. (loc%z > maxl%z)) return
 
is_location_in_region = .true.

end function is_location_in_region

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

subroutine convert_vertical_obs(ens_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, status)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer,             intent(in)    :: loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: status(:)

status(:) = 0

end subroutine convert_vertical_obs

!--------------------------------------------------------------------

subroutine convert_vertical_state(ens_handle, num, locs, loc_qtys, loc_indx, &
                                  which_vert, status)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer(i8),         intent(in)    :: loc_indx(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: status

status = 0

end subroutine convert_vertical_state

!----------------------------------------------------------------------------
! end of location/threed/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod

