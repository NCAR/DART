! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!>@todo the channel location_mod.html needs to be written.

module location_mod

! Implements location interfaces for a 3d channel in X,Y,Z where X is periodic,
! Y has walls (limited domain), and Z is infinite

use      types_mod, only : r8, i8, MISSING_R8, MISSING_I, PI, RAD2DEG, DEG2RAD
use  utilities_mod, only : error_handler, E_ERR, ascii_file_format, &
                           E_MSG, open_file, close_file, set_output,                 &
                           logfileunit, nmlfileunit, find_namelist_in_file,          &
                           check_namelist_read, do_output, do_nml_file,              &
                           do_nml_term, is_longitude_between
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform
use mpi_utilities_mod, only : my_task_id, task_count
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
          set_vertical_localization_coord, convert_vertical_obs, convert_vertical_state, &
          print_get_close_type, find_nearest

character(len=*), parameter :: source = 'channel/location_mod.f90'

type location_type
   private
   real(r8) :: x, y, z
end type location_type

type get_close_type
   private
   integer           :: num
   real(r8)          :: maxdist
   integer, pointer  :: loc_box(:)           ! (nloc); List of loc indices in boxes
   integer, pointer  :: count(:, :, :)       ! (nx, ny, nz); # of locs in each box
   integer, pointer  :: start(:, :, :)       ! (nx, ny, nz); Start of list of locs in this box
   real(r8)          :: bot_x, top_x         ! extents in x, y, z
   real(r8)          :: bot_y, top_y 
   real(r8)          :: bot_z, top_z 
   real(r8)          :: x_width, y_width, z_width    ! widths of boxes in x,y,z
   real(r8)          :: nboxes_x, nboxes_y, nboxes_z ! based on maxdist how far to search
end type get_close_type

type(random_seq_type) :: ran_seq
logical               :: ran_seq_init = .false.
logical, save         :: module_initialized = .false.

integer,              parameter :: LocationDims = 3
character(len = 129), parameter :: LocationName = "loc3Dchan"
character(len = 129), parameter :: LocationLName = &
                                   "threed channel locations: x, y, z"
character(len = 129), parameter :: LocationStorageOrder = "X Y Z"
character(len = 129), parameter :: LocationUnits = "none none none"

character(len = 512) :: errstring

integer :: nx               = 10
integer :: ny               = 10
integer :: nz               = 10

!-----------------------------------------------------------------
! Namelist with default values

logical :: output_box_info  = .false.
integer :: print_box_level  = 0
! tuning options
integer :: nboxes           = 1000 ! compute nx/ny/nz based on cube root of this #

! Option for verification using exhaustive search
! .true. means recompute everything - expensive
logical :: compare_to_correct = .false.

namelist /location_nml/ nboxes, compare_to_correct, &
                        output_box_info, print_box_level

!-----------------------------------------------------------------

interface operator(==); module procedure loc_eq; end interface
interface operator(/=); module procedure loc_ne; end interface

interface set_location
   module procedure set_location_single
   module procedure set_location_array
end interface set_location

contains

!----------------------------------------------------------------------------

subroutine initialize_module
 
! things which need doing exactly once.

integer :: iunit, io

if (module_initialized) return

module_initialized = .true.

! Read the namelist entry
call find_namelist_in_file("input.nml", "location_nml", iunit)
read(iunit, nml = location_nml, iostat = io)
call check_namelist_read(iunit, io, "location_nml")

! Write the namelist values to the log file

if(do_nml_file()) write(nmlfileunit, nml=location_nml)
if(do_nml_term()) write(     *     , nml=location_nml)

end subroutine initialize_module

!----------------------------------------------------------------------------

function get_dist(loc1, loc2, type1, kind2)

! returns the distance between 2 locations 

! In spite of the names, the 3rd and 4th argument are actually specific types
! (e.g. RADIOSONDE_TEMPERATURE, AIRCRAFT_TEMPERATURE).  The types are part of
! the interface in case user-code wants to do a more sophisticated distance
! calculation based on the base or target types.  In the usual case this
! code still doesn't use the types, but there's an undocumented feature that
! allows you to maintain the original vertical normalization even when
! changing the cutoff distance in the horizontal.  For that to work we
! do need to know the type, and we use the type of loc1 to control it.
! 

type(location_type), intent(in) :: loc1, loc2
integer, optional,   intent(in) :: type1, kind2
real(r8)                        :: get_dist

real(r8) :: x_dif, y_dif, z_dif

if ( .not. module_initialized ) call initialize_module

x_dif = loc1%x - loc2%x
y_dif = loc1%y - loc2%y
z_dif = loc1%z - loc2%z

get_dist = sqrt(x_dif * x_dif + y_dif * y_dif + z_dif * z_dif)

end function get_dist

!---------------------------------------------------------------------------

function loc_eq(loc1,loc2)
 
! Interface operator used to compare two locations.
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
 
! Interface operator used to compare two locations.
! Returns true if locations are not identical to machine precision.

type(location_type), intent(in) :: loc1, loc2
logical                         :: loc_ne

if ( .not. module_initialized ) call initialize_module

loc_ne = (.not. loc_eq(loc1,loc2))

end function loc_ne

!---------------------------------------------------------------------------

function get_location(loc)
 
! Given a location type return the x,y,z coordinates

type(location_type), intent(in) :: loc
real(r8), dimension(3) :: get_location

if ( .not. module_initialized ) call initialize_module

get_location(1) = loc%x
get_location(2) = loc%y
get_location(3) = loc%z

end function get_location

!---------------------------------------------------------------------------

function set_location_single(x, y, z)
 
! Puts the x, y, z into a location datatype.

real(r8), intent(in) :: x, y, z
type (location_type) :: set_location_single

if ( .not. module_initialized ) call initialize_module

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

! Initialize a location type to indicate the contents are unset.

type (location_type) :: set_location_missing

if ( .not. module_initialized ) call initialize_module

set_location_missing%x = MISSING_R8
set_location_missing%y = MISSING_R8
set_location_missing%z = MISSING_R8

end function set_location_missing

!---------------------------------------------------------------------------

function query_location(loc, attr)
 
! Returns the value of the attribute

type(location_type),        intent(in) :: loc
character(len=*), optional, intent(in) :: attr
real(r8)                               :: query_location

if ( .not. module_initialized ) call initialize_module

! before you change any of the code in this subroutine,
! check out the extensive comments in the threed_sphere
! version of this routine.  then proceed with caution.

if (.not. present(attr)) then
   query_location = loc%x
   return
endif

select case(attr)
   case ('x','X')
      query_location = loc%x
   case ('y','Y')
      query_location = loc%y
   case ('z','Z')
      query_location = loc%z
   case default
      call error_handler(E_ERR, 'query_location:', &
         'Only "X","Y","Z" are legal attributes to request from location', source)
end select

end function query_location

!----------------------------------------------------------------------------

subroutine write_location(locfile, loc, fform, charstring)
 
! Writes a location to a file.
! most recent change: adding the optional charstring option.  if present,
! locfile is ignored, and a pretty-print formatting is done into charstring.

integer, intent(in)                        :: locfile
type(location_type), intent(in)            :: loc
character(len = *),  intent(in),  optional :: fform
character(len = *),  intent(out), optional :: charstring

integer             :: charlength
logical             :: writebuf

10 format(1X,3(G25.16,1X))

if ( .not. module_initialized ) call initialize_module

! writing to a file (normal use) or to a character buffer?
writebuf = present(charstring)

! output file; test for ascii or binary, write what's asked, and return
if (.not. writebuf) then
   if (ascii_file_format(fform)) then
      write(locfile, '(''loc3Dchan'')' ) 
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

! format the location to be human-friendly

! the output can be no longer than this
charlength = 70

if (len(charstring) < charlength) then
   write(errstring, *) 'charstring buffer must be at least ', charlength, ' chars long'
   call error_handler(E_ERR, 'write_location', errstring, source)
endif

! format into the outout string
write(charstring, '(A,3(G20.8,1X))') 'X/Y/Z: ', loc%x, loc%y, loc%z

end subroutine write_location

!----------------------------------------------------------------------------

function read_location(locfile, fform)
 
! Reads a location from a file that was written by write_location. 
! See write_location for additional discussion.

integer, intent(in)                      :: locfile
character(len = *), intent(in), optional :: fform
type(location_type)                      :: read_location

character(len=9) :: header

if ( .not. module_initialized ) call initialize_module

if (ascii_file_format(fform)) then
   read(locfile, '(A9)' ) header
   if(header /= 'loc3Dchan') then
         write(errstring,*)'Expected location header "loc3Dchan" in input file, got ', header
      call error_handler(E_ERR, 'read_location', errstring, source)
   endif
   ! Now read the location data value
   read(locfile, *)read_location%x, read_location%y, read_location%z
else
   read(locfile)read_location%x, read_location%y, read_location%z
endif

end function read_location

!--------------------------------------------------------------------------

subroutine interactive_location(location, set_to_default)
 
! Allows for interactive input of a location. Also gives option of selecting
! a uniformly distributed random location.

type(location_type), intent(out) :: location
logical, intent(in), optional    :: set_to_default

real(r8) :: v(3), minv, maxv
character(len=1) :: l(3)
integer :: i, r

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

! prompt for an explicit location or a random one.
! if random, generate all 3 x/y/z values randomly.
! if you want to make some combination of x/y/z random
! and specify others, you would have to move this read 
! into the loop.

r = 1
do while (r > 0) 

   write(*, *) 'Input 0 to specify a value for the location, or'
   write(*, *) '-1 for a uniformly distributed random location'
   read(*, *) r

   if (r > 0) write(*, *) 'Please input 0 or -1 for selection'
enddo

do i = 1, 3
   if (r == 0) then
      write(*, *) 'Input value for ', l(i)
      read (*,*) v(i)

   else 
      ! Need to make sure random sequence is initialized
   
      if(.not. ran_seq_init) then
         call init_random_seq(ran_seq)
         ran_seq_init = .TRUE.
      endif
   
      write(*, *) 'Input minimum ', l(i), ' value '
      read(*, *) minv
   
      write(*, *) 'Input maximum ', l(i), ' value '
      read(*, *) maxv
   
      v(i) = random_uniform(ran_seq) * (maxv-minv) + minv
   
      write(*, *) 'random location is ', v(i)
   
   endif
   
enddo

location%x = v(1)
location%y = v(2)
location%z = v(3)

end subroutine interactive_location

!----------------------------------------------------------------------------

subroutine get_close_init(gc, num, maxdist, locs, maxdist_list)
 
! Initializes part of get_close accelerator that depends on the particular obs

type(get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
real(r8),             intent(in)    :: maxdist
type(location_type),  intent(in)    :: locs(:)
real(r8), intent(in), optional      :: maxdist_list(:)

integer :: i, j, k, cum_start, l
integer :: x_box(num), y_box(num), z_box(num)
integer :: tstart(nx, ny, nz)

if ( .not. module_initialized ) call initialize_module

! Set the maximum localization distance
gc%maxdist = maxdist

! Allocate storage for number dependent part
allocate(gc%loc_box(num))
gc%loc_box(:) = -1

! FIXME: compute nx, ny, nz from nboxes?  or put in namelist
nx = nint(real(nboxes, r8)**0.33333)   ! roughly cube root
ny = nint(real(nboxes, r8)**0.33333)   ! roughly cube root
nz = nint(real(nboxes, r8) / real(nx * ny, r8))  ! whatever is left

! Allocate the storage for the grid dependent boxes
allocate(gc%count(nx,ny,nz), gc%start(nx,ny,nz))
gc%count  = -1
gc%start  = -1

! Set the value of num_locs in the structure
gc%num = num

! If num == 0, no point in going any further.
if (num == 0) return


! Allocate storage for obs number dependent part
allocate(gc%loc_box(num))
gc%loc_box(:) = -1

! Determine where the boxes should be for this set of locs and maxdist
call find_box_ranges(gc, locs, num)

! Begin by computing the number of locations in each box in x,y,z
gc%count = 0
do i = 1, num

!print *, i, locs(i)%x, locs(i)%y, locs(i)%z
   x_box(i) = floor((locs(i)%x - gc%bot_x) / gc%x_width) + 1
   if(x_box(i) > nx) x_box(i) = nx
   if(x_box(i) < 1)  x_box(i) = 1

   y_box(i) = floor((locs(i)%y - gc%bot_y) / gc%y_width) + 1
   if(y_box(i) > ny) y_box(i) = ny
   if(y_box(i) < 1)  y_box(i) = 1

   z_box(i) = floor((locs(i)%z - gc%bot_z) / gc%z_width) + 1
   if(z_box(i) > nz) z_box(i) = nz
   if(z_box(i) < 1)  z_box(i) = 1

   gc%count(x_box(i), y_box(i), z_box(i)) = gc%count(x_box(i), y_box(i), z_box(i)) + 1
!print *, 'adding count to box ', x_box(i), y_box(i), z_box(i), &
!                                 gc%count(x_box(i), y_box(i), z_box(i))
end do

! Figure out where storage for each boxes members should begin
cum_start = 1
do i = 1, nx
   do j = 1, ny
      do k = 1, nz
         gc%start(i, j, k) = cum_start
         cum_start = cum_start + gc%count(i, j, k)
      end do
   end do
end do

! Now we know how many are in each box, get a list of which are in each box
tstart = gc%start
do i = 1, num
   gc%loc_box(tstart(x_box(i), y_box(i), z_box(i))) = i
   tstart(x_box(i), y_box(i), z_box(i)) = tstart(x_box(i), y_box(i), z_box(i)) + 1
end do

do i = 1, nx
   do j = 1, ny
      do k = 1, nz
if (gc%count(i,j,k) > 0) print *, i,j,k, gc%count(i,j,k), gc%start(i,j,k)
         do l=1, gc%count(i,j,k)
!print *, l, gc%loc_box(l)
         enddo
      end do
   end do
end do


! info on how well the boxes are working.  by default print nothing.
! set print_box_level to higher values to get more and more detail.
! user info should be level 1; 2 and 3 should be for debug only.
! special for grid-decomposition debugging; set print level to -8.
if (output_box_info) then
   ! if this task normally prints, call the print routine.
   ! if print level > 2, set all tasks to print and call print.
   ! then reset the status to off again.
   if (do_output()) then
      call print_get_close_type(gc, print_box_level)
   else if (print_box_level >= 2 .or. print_box_level < 0) then
      ! print status was false, but turn on temporarily
      ! to output box info from all tasks.
      call set_output(.true.)
      call print_get_close_type(gc, print_box_level)
      call set_output(.false.)
   endif
endif

end subroutine get_close_init

!----------------------------------------------------------------------------

subroutine get_close_destroy(gc)

type(get_close_type), intent(inout) :: gc

deallocate(gc%loc_box, gc%count, gc%start)

end subroutine get_close_destroy

!----------------------------------------------------------------------------

subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

! The specific type of the base observation, plus the generic kinds list
! for either the state or obs lists are available if a more sophisticated
! distance computation is needed.

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

subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)

! The specific type of the base observation, plus the generic kinds list
! for either the state or obs lists are available if a more sophisticated
! distance computation is needed.

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

subroutine get_close(gc, base_loc, base_type, locs, loc_qtys, &
                     num_close, close_ind, dist, ens_handle)

type(get_close_type), intent(in)  :: gc
type(location_type),  intent(in)  :: base_loc,  locs(:)
integer,              intent(in)  :: base_type, loc_qtys(:)
integer,              intent(out) :: num_close, close_ind(:)
real(r8), optional,   intent(out) :: dist(:)
type(ensemble_type),  intent(in)  :: ens_handle

! If dist is NOT present, just find everybody in a box, put them in the list,
! but don't compute any distances

integer :: x_box, y_box, z_box, i, j, k, l
integer :: start_x, end_x, start_y, end_y, start_z, end_z
integer ::  n_in_box, st, t_ind
real(r8) :: this_dist, this_maxdist

! Variables needed for comparing against correct case.
! these could be large - make them allocatable
! and only allocate them if needed.
integer :: cnum_close
integer, allocatable :: cclose_ind(:)
real(r8), allocatable :: cdist(:)

! First, set the intent out arguments to a missing value
num_close = 0
close_ind = -99
if(present(dist)) dist = -1e38_r8  ! big but negative
this_dist = 1e38_r8                ! something big and positive.

! the list of locations in the locs() argument must be the same
! as the list of locations passed into get_close_init(), so
! gc%num and size(locs) better be the same.   if the list changes,
! you have to destroy the old gc and init a new one.
if (size(locs) /= gc%num) then
   write(errstring,*)'locs() array must match one passed to get_close_init()'
   call error_handler(E_ERR, 'get_close', errstring, source)
endif

! If num == 0, no point in going any further. 
if (gc%num == 0) return

this_maxdist = gc%maxdist


! For validation, it is useful to be able to compare against exact
! exhaustive search
if(compare_to_correct) then
   call exhaustive_collect(gc, base_loc, locs, &
                           cnum_close, cclose_ind, cdist)
endif


! Begin by figuring out which box the base locs is in
x_box = floor((base_loc%x - gc%bot_x) / gc%x_width) + 1
y_box = floor((base_loc%y - gc%bot_y) / gc%y_width) + 1
z_box = floor((base_loc%z - gc%bot_z) / gc%z_width) + 1

! If it is not in any box, then it is more than the maxdist away from everybody
if(x_box > nx .or. x_box < 1 .or. x_box < 0) return
if(y_box > ny .or. y_box < 1 .or. y_box < 0) return
if(z_box > nz .or. z_box < 1 .or. z_box < 0) return

! figure out how many boxes need searching
! FIXME: if we support a variable maxdist, nboxes_X will need to
! be computed on the fly here instead of precomputed at init time.
!print *, 'nboxes x, y, z = ', gc%nboxes_x, gc%nboxes_y, gc%nboxes_z
!print *, 'base_loc in box ', x_box, y_box, z_box

start_x = x_box - gc%nboxes_x
if (start_x < 1) start_x = 1
end_x = x_box + gc%nboxes_x
if (end_x > nx) end_x = nx

start_y = y_box - gc%nboxes_y
if (start_y < 1) start_y = 1
end_y = y_box + gc%nboxes_y
if (end_y > ny) end_y = ny

start_z = z_box - gc%nboxes_z
if (start_z < 1) start_z = 1
end_z = z_box + gc%nboxes_z
if (end_z > nz) end_z = nz

!print *, 'looping from '
!print *, 'x: ', start_x, end_x
!print *, 'y: ', start_y, end_y
!print *, 'z: ', start_z, end_z

! Next, loop through each box that is close to this box
do i = start_x, end_x
   do j = start_y, end_y
      do k = start_z, end_z

         ! Box to search is i,j,k
         n_in_box = gc%count(i, j, k)
         st = gc%start(i,j,k)


         ! Loop to check how close all loc in the box are; add those that are close
         do l = 1, n_in_box

            t_ind = gc%loc_box(st - 1 + l)
!print *, 'l, t_ind = ', l, t_ind

            ! Only compute distance if dist is present
            if(present(dist)) then
               this_dist = get_dist(base_loc, locs(t_ind))
!print *, 'this_dist = ', this_dist
               ! If this loc's distance is less than cutoff, add it in list
               if(this_dist <= this_maxdist) then
                  num_close = num_close + 1
                  close_ind(num_close) = t_ind
                  dist(num_close) = this_dist
               endif
            else
               ! Dist isn't present; add this ob to list without computing distance
               num_close = num_close + 1
               close_ind(num_close) = t_ind
            endif

         end do
      end do
   end do
end do


! Verify by comparing to exhaustive search
if(compare_to_correct) then
   call exhaustive_report(cnum_close, num_close, cclose_ind, close_ind, cdist, dist)
endif


end subroutine get_close

!--------------------------------------------------------------------------

subroutine find_box_ranges(gc, locs, num)
 
! Finds boundaries for x,y,z boxes.  
! FIXME: ways boxes could be divided:
!  - evenly along each axis
!  - octree-like, divide each axis so roughly half the points are
!     on each side of the dividing plane.
!  - about 100 other schemes
  
type(get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
type(location_type),  intent(in)    :: locs(num)


! FIXME: this space could be very sparse

gc%bot_x = minval(locs(:)%x)
gc%bot_y = minval(locs(:)%y)
gc%bot_z = minval(locs(:)%z)

gc%top_x = maxval(locs(:)%x)
gc%top_y = maxval(locs(:)%y)
gc%top_z = maxval(locs(:)%z)

gc%x_width = (gc%top_x - gc%bot_x) / nx
gc%y_width = (gc%top_y - gc%bot_y) / ny 
gc%z_width = (gc%top_z - gc%bot_z) / nz

! FIXME:  compute a sphere of radius maxdist and see how
! many boxes in x, y, z that would include.
gc%nboxes_x = aint((gc%maxdist + (gc%x_width-1)) / gc%x_width) 
gc%nboxes_y = aint((gc%maxdist + (gc%y_width-1)) / gc%y_width) 
gc%nboxes_z = aint((gc%maxdist + (gc%z_width-1)) / gc%z_width) 

end subroutine find_box_ranges

!----------------------------------------------------------------------------

subroutine find_nearest(gc, base_loc, loc_list, nearest, rc)
 type(get_close_type), intent(in), target  :: gc
 type(location_type),  intent(in)  :: base_loc
 type(location_type),  intent(in)  :: loc_list(:)
 integer,              intent(out) :: nearest
 integer,              intent(out) :: rc

! find the nearest point in the get close list to the specified point

integer :: x_box, y_box, z_box, i, j, k, l
integer :: start_x, end_x, start_y, end_y, start_z, end_z
integer ::  n_in_box, st, t_ind
real(r8) :: this_dist, dist

! First, set the intent out arguments to a missing value
nearest = -99
rc = -1
dist = 1e38_r8                ! something big and positive.

! the list of locations in the loc() argument must be the same
! as the list of locations passed into get_close_init(), so
! gc%num and size(loc) better be the same.   if the list changes,
! you have to destroy the old gc and init a new one.
if (size(loc_list) /= gc%num) then
   write(errstring,*)'loc() array must match one passed to get_close_init()'
   call error_handler(E_ERR, 'find_nearest_boxes', errstring, source)
endif

! If num == 0, no point in going any further. 
if (gc%num == 0) return

!--------------------------------------------------------------

! Begin by figuring out which box the base loc is in
x_box = floor((base_loc%x - gc%bot_x) / gc%x_width) + 1
y_box = floor((base_loc%y - gc%bot_y) / gc%y_width) + 1
z_box = floor((base_loc%z - gc%bot_z) / gc%z_width) + 1

! If it is not in any box, then it is more than the maxdist away from everybody
if(x_box > nx .or. x_box < 1 .or. x_box < 0) return
if(y_box > ny .or. y_box < 1 .or. y_box < 0) return
if(z_box > nz .or. z_box < 1 .or. z_box < 0) return

! figure out how many boxes need searching
! FIXME: if we support a variable maxdist, nboxes_X will need to
! be computed on the fly here instead of precomputed at init time.
!print *, 'nboxes x, y, z = ', gc%nboxes_x, gc%nboxes_y, gc%nboxes_z
!print *, 'base_loc in box ', x_box, y_box, z_box

start_x = x_box - gc%nboxes_x
if (start_x < 1) start_x = 1
end_x = x_box + gc%nboxes_x
if (end_x > nx) end_x = nx

start_y = y_box - gc%nboxes_y
if (start_y < 1) start_y = 1
end_y = y_box + gc%nboxes_y
if (end_y > ny) end_y = ny

start_z = z_box - gc%nboxes_z
if (start_z < 1) start_z = 1
end_z = z_box + gc%nboxes_z
if (end_z > nz) end_z = nz

!print *, 'looping from '
!print *, 'x: ', start_x, end_x
!print *, 'y: ', start_y, end_y
!print *, 'z: ', start_z, end_z

! Next, loop through each box that is close to this box
do i = start_x, end_x
   do j = start_y, end_y
      do k = start_z, end_z

         ! Box to search is i,j,k
         n_in_box = gc%count(i, j, k)
         st = gc%start(i,j,k)


         ! Loop to check how close all loc in the box are; add those that are close
         do l = 1, n_in_box

            t_ind = gc%loc_box(st - 1 + l)
!print *, 'l, t_ind = ', l, t_ind

            this_dist = get_dist(base_loc, loc_list(t_ind))
!print *, 'this_dist = ', this_dist
            ! If this loc's distance is less than current nearest, it's new nearest
            if(this_dist <= dist) then
               nearest = t_ind
               dist = this_dist
               if (rc < 0) rc = 0
            endif
         end do
      end do
   end do
end do

end subroutine find_nearest

!---------------------------------------------------------------------------

function get_maxdist(gc, obs_type)
type(get_close_type), intent(in) :: gc
integer, optional,    intent(in) :: obs_type
real(r8) :: get_maxdist

get_maxdist = gc%maxdist

end function get_maxdist

!----------------------------------------------------------------------------

subroutine print_get_close_type(gc, amount)
 
! print out debugging statistics, or optionally print out a full
! dump from all mpi tasks in a format that can be plotted with matlab.

type(get_close_type), intent(in), target :: gc
integer, intent(in), optional            :: amount

integer :: i, j, k, l, first, index, mytask, alltasks
integer :: sample, nfull, nempty, howmuch, total, maxcount, maxi, maxj, maxk
logical :: tickmark(gc%num), iam0
real(r8) :: x_cen, y_cen, z_cen

logical, save :: write_now = .true.
integer, save :: been_called = 0
integer :: funit
character(len=64) :: fname

! cumulative times through this routine
been_called = been_called + 1

! second arg is now an int, not logical, and means:
! 0 = very terse, only box summary (default).  
! 1 = structs and first part of arrays.
! 2 = all parts of all arrays.
! -8 = special for grid-decomposition debugging

! by default do not print all the loc_box or start contents (it can
! be very long).  but give the option to print more info or even an
! entire contents dump.  'sample' is the number to print for the
! short version.  (this value prints about 5-6 lines of data.)
! to get a full dump, change print_box_level to 2 or more in the namelist.
howmuch = 0
sample = 10
mytask = my_task_id() 
alltasks = task_count()
iam0 = (mytask == 0)

if (present(amount)) then
   howmuch = amount
endif

if (howmuch == -8) then
   if (.not. write_now) howmuch = 0
endif

!! SPECIAL - debugging
! if you enable debugging, maybe you want to turn it off for really
! large counts?  often it's easy to construct a case that has a lot of
! locations from the state vector in one set of boxes, but just a few
! locations from the observations in another.  this lets you turn off
! the debugging level for the large set and leave it on for the small.
!if (gc%num > 100) howmuch = 0

! print the get_close_type derived type values

if (howmuch /= 0 .and. iam0) then
   write(errstring,*) 'get_close_type values:'
   call error_handler(E_MSG, 'loc', errstring)

   write(errstring,*) ' num = ', gc%num
   call error_handler(E_MSG, 'loc', errstring)

   write(errstring,*) ' nx, ny, nz = ', nx, ny, nz
   call error_handler(E_MSG, 'loc', errstring)

   write(errstring,"(A,F12.6)") ' maxdist = ', gc%maxdist
   call error_handler(E_MSG, 'loc', errstring)
   write(errstring, "(A,3(F12.6))") ' x_box: bot, top, width = ', gc%bot_x, gc%top_x, gc%x_width
   call error_handler(E_MSG, 'loc', errstring)
   write(errstring, "(A,3(F12.6))") ' y_box: bot, top, width = ', gc%bot_y, gc%top_y, gc%y_width
   call error_handler(E_MSG, 'loc', errstring)
   write(errstring, "(A,3(F12.6))") ' z_box: bot, top, width = ', gc%bot_z, gc%top_z, gc%z_width
   call error_handler(E_MSG, 'loc', errstring)

endif

! this one can be very large.   print only the first nth unless
! instructed otherwise.  (print n+1 because 1 more value fits on
! the line because it prints ( i ) and not ( i, j ) like the others.)
if (associated(gc%loc_box)) then
   i = size(gc%loc_box,1)
   if (i/= gc%num) then
      write(errstring,*) ' warning: size of loc_box incorrect, nlocs, i =', gc%num, i
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
   if (howmuch > 1) then
      ! DEBUG
      write(errstring,"(A,I8,A,36(I8,1X))") ' loc_box(',i,') =', gc%loc_box(1:min(i,36))  ! (nlocs)
      !write(errstring,*) ' loc_box(',i,') =', gc%loc_box    ! (nlocs)
      call error_handler(E_MSG, 'locations_mod', errstring)
   else if(howmuch > 0) then
      write(errstring,*) ' loc_box(',i,') =', gc%loc_box(1:min(i,sample+1))
      call error_handler(E_MSG, 'locations_mod', errstring)
      write(errstring,*) '  <rest of loc_box omitted>'
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
else
   if (howmuch > 0) then
      write(errstring,*) ' loc_box unallocated'
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
endif

! like loc_box, this one can be very large.   print only the first nth unless
! instructed otherwise
if (associated(gc%start)) then
   i = size(gc%start,1)
   j = size(gc%start,2)
   k = size(gc%start,3)
   if ((i /= nx) .or. (j /= ny) .or. (k /= nz)) then
      write(errstring,*) ' warning: size of start incorrect, nx, ny, nz, i, j, k =', nx, ny, nz, i, j, k
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
   if (howmuch > 1) then
      write(errstring,*) ' start(',i,j,k,') ='              ! (nx, ny, nz)
      call error_handler(E_MSG, 'locations_mod', errstring)
      do l=1, j
         write(errstring,"(36(I8,1X))") gc%start(1:min(i,36), l, 1)
         call error_handler(E_MSG, 'locations_mod', errstring)
      enddo
   else if (howmuch > 0) then
      write(errstring,*) ' start(',i,j,k,') =', gc%start(1:min(i,sample), 1, 1)
      call error_handler(E_MSG, 'locations_mod', errstring)
      write(errstring,*) '  <rest of start omitted>'
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
else
   if (howmuch > 0) then
      write(errstring,*) ' start unallocated'
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
endif

! as above, print only first n unless second arg is .true.
if (associated(gc%count)) then
   i = size(gc%count,1)
   j = size(gc%count,2)
   k = size(gc%count,3)
   if ((i /= nx) .or. (j /= ny) .or. (k /= nz)) then
      write(errstring,*) ' warning: size of count incorrect, nx, ny, nz, i, j, k =', nx, ny, nz, i, j, k
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
   if (howmuch > 1) then
      write(errstring,*) ' count(',i,j,k,') ='              ! (nx, ny, nz)
      call error_handler(E_MSG, 'locations_mod', errstring)
      do l=1, j
         write(errstring,"(36(I8,1X))") gc%count(1:min(i,36), l, 1) 
         call error_handler(E_MSG, 'locations_mod', errstring)
      enddo
   else if (howmuch > 0) then
      write(errstring,*) ' count(',i,j,k,') =', gc%count(1:min(i,sample), 1, 1)
      call error_handler(E_MSG, 'locations_mod', errstring)
      write(errstring,*) '  <rest of count omitted>'
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
else
   if (howmuch > 0) then
      write(errstring,*) ' count unallocated'
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
endif


! end of print.  strictly speaking the following code is validation code,
! not a simple print.


! initialize all ticks to false.  turn them true as they are found
! in the loc_box list, and complain about duplicates or misses.
tickmark = .FALSE.

do i=1, nx
   do j=1, ny
      do k=1, nz
         first = gc%start(i, j, k)
         do l=1, gc%count(i, j, k)
            index = first + l - 1
            if ((index < 1) .or. (index > gc%num)) then
               write(errstring, *) 'exiting at first bad value; could be more'
               call error_handler(E_MSG, 'locations_mod', errstring)
               write(errstring, *) 'bad locs list index, in box: ', index, i, j, k
               call error_handler(E_ERR, 'locations_mod', errstring)
            endif
            if (tickmark(index)) then
               write(errstring, *) 'exiting at first bad value; could be more'
               call error_handler(E_MSG, 'locations_mod', errstring)
               write(errstring, *) 'error: loc found in more than one box list.  index, box: ', &
                            index, i, j, k
               call error_handler(E_ERR, 'locations_mod', errstring)
            endif
            tickmark(index) = .TRUE.
         enddo
      enddo
   enddo
enddo

do i=1, gc%num
  if (.not. tickmark(i)) then
     write(errstring, *) 'exiting at first bad value; could be more'
     call error_handler(E_MSG, 'locations_mod', errstring)
     write(errstring,*) 'loc not found in any box list: ', i
     call error_handler(E_ERR, 'locations_mod', errstring)
  endif
enddo

! print out some hopefully useful stats
nfull = 0
nempty = 0
total = 0
maxcount = 0
maxi = 0
maxj = 0

if (howmuch == -8) then
   if (iam0) then
      fname = 'loc_dump_header.m'
      funit = open_file(fname, action='write')
      write(funit,'(A,I2,A,I4,A)') 'xlocs = zeros(', nx, ',', alltasks, ');'
      write(funit,'(A,I2,A,I4,A)') 'ylocs = zeros(', ny, ',', alltasks, ');'
      write(funit,'(A,I2,A,I4,A)') 'zlocs = zeros(', nz, ',', alltasks, ');'
      write(funit,'(3(A,I2),A,I4,A)') 'boxes = zeros(', nx, ',', ny, ',', nz, ',', alltasks, ');'
      call close_file(funit)
   endif
   write(fname, '(A,I3.3,A)')  'loc_dump_', mytask, '.m'
   funit = open_file(fname, action='write')
endif

do i=1, nx
   if (howmuch == -8) then
      x_cen = gc%bot_x + ((i-1)*gc%x_width) + (gc%x_width/2.0)
      write(funit, '(A,I2,A,I4,A,F12.9,A)') 'xlocs(', i, ',', mytask+1, ') = ',  x_cen, ';'
   endif
   do j=1, ny
      if (howmuch == -8 .and. i==1) then
         y_cen = gc%bot_y + ((j-1)*gc%y_width) + (gc%y_width/2.0)
         write(funit, '(A,I2,A,I4,A,F12.9,A)') 'ylocs(', j, ',', mytask+1, ') = ',  y_cen, ';'
      endif
      do k=1, nz
         if (howmuch == -8 .and. i==1) then
            z_cen = gc%bot_z + ((j-1)*gc%z_width) + (gc%z_width/2.0)
            write(funit, '(A,I2,A,I4,A,F12.9,A)') 'zlocs(', k, ',', mytask+1, ') = ',  z_cen, ';'
         endif
         if (gc%count(i, j, k) > 0) then
            nfull = nfull + 1
            total = total + gc%count(i, j, k)
            if (gc%count(i, j, k) > maxcount) then
               maxcount = gc%count(i, j, k)
               maxi = i
               maxj = j
               maxk = k
            endif
         else
            nempty = nempty + 1
         endif
         ! output for grid boxes; in matlab-friendly format
         if (howmuch == -8) then
            write(funit, '(3(A,I2),A,I4,A,I8,A)') 'boxes(', i, ', ', j, ', ', k, &
                                   ',', mytask+1, ') = ', gc%count(i, j, k), ';'
         endif
      enddo
   enddo
enddo

if (howmuch == -8) then
   call close_file(funit)
   write_now = .false.
endif

! these print out always - make sure they are useful to end users.
write(errstring, '(a)') "Location module statistics:"
call error_handler(E_MSG, 'locations_mod', errstring)
write(errstring, '(a,i9)') " Total boxes (nx * ny * nz): ", nfull + nempty
call error_handler(E_MSG, 'locations_mod', errstring)
write(errstring, '(a,i9)') " Total items to put in boxes: ", gc%num
call error_handler(E_MSG, 'locations_mod', errstring)
if (howmuch > 0) then
   write(errstring, '(a,i9)') " Total boxes with 1+ items: ", nfull
   call error_handler(E_MSG, 'locations_mod', errstring)
   write(errstring, '(a,i9)') " Total boxes empty: ", nempty
   call error_handler(E_MSG, 'locations_mod', errstring)
endif
if (nfull > 0) then
   write(errstring, '(a,f7.2)') " Percent boxes with 1+ items: ", nfull / real(nfull + nempty, r8) * 100.
   call error_handler(E_MSG, 'locations_mod', errstring)
   write(errstring, '(a,f12.2)') " Average #items per non-empty box: ", real(total, r8) / nfull
   call error_handler(E_MSG, 'locations_mod', errstring)
endif
if (maxcount > 0) then
   write(errstring, '(a,i9)') " Largest #items in one box: ", maxcount
   call error_handler(E_MSG, 'locations_mod', errstring)
! leave this out for now.  one, if there are multiple boxes with
! the same maxcount this is just the last one found.  two, the
! index numbers do not seem very helpful.
!   if (howmuch > 0) then
!      write(errstring, '(a,i9,i9)') " That box index: ", maxi, maxj, maxz
!      call error_handler(E_MSG, 'locations_mod', errstring)
!   endif
endif

end subroutine print_get_close_type

!----------------------------------------------------------------------------

function is_location_in_region(loc, minl, maxl)
 
! Returns true if the given location is inside the rectangular
! region defined by minl as the lower left, maxl the upper right.
! test is inclusive; values on the edges are considered inside.

logical                          :: is_location_in_region
type(location_type), intent(in)  :: loc, minl, maxl

if ( .not. module_initialized ) call initialize_module

! assume failure and return as soon as we are confirmed right.
! set to success only at the bottom after all tests have passed.
is_location_in_region = .false.

if ((loc%x < minl%x) .or. (loc%x > maxl%x)) return
if ((loc%y < minl%y) .or. (loc%y > maxl%y)) return 
if ((loc%z < minl%z) .or. (loc%z > maxl%z)) return
 
is_location_in_region = .true.

end function is_location_in_region

!---------------------------------------------------------------------------

subroutine exhaustive_collect(gc, base_loc, loc_list, num_close, close_ind, close_dist)

! For validation, it is useful to be able to compare against exact
! exhaustive search

type(get_close_type),  intent(in)  :: gc
type(location_type),   intent(in)  :: base_loc, loc_list(:)
integer,               intent(out) :: num_close
integer,  allocatable, intent(out) :: close_ind(:)
real(r8), allocatable, intent(out) :: close_dist(:)

real(r8) :: this_dist
integer :: i

allocate(close_ind(size(loc_list)), close_dist(size(loc_list)))
num_close = 0
do i = 1, gc%num 
   this_dist = get_dist(base_loc, loc_list(i))
   if(this_dist <= gc%maxdist) then
      ! Add this loc to correct list
      num_close = num_close + 1
      close_ind(num_close) = i
      close_dist(num_close) = this_dist
   endif
end do

end subroutine exhaustive_collect

!---------------------------------------------------------------------------

subroutine exhaustive_report(cnum_close, num_close, cclose_ind, close_ind, cclose_dist, close_dist)

! For validation, it is useful to be able to compare against exact
! exhaustive search

integer,                         intent(in)    :: cnum_close, num_close
integer,  allocatable,           intent(inout) :: cclose_ind(:)
integer,                         intent(in)    :: close_ind(:)
real(r8), allocatable,           intent(inout) :: cclose_dist(:)
real(r8),              optional, intent(in)    :: close_dist(:)

! Do comparisons against full search
if((num_close /= cnum_close) .and. present(close_dist)) then
   write(errstring, *) 'get_close (', num_close, ') should equal exhaustive search (', cnum_close, ')'
   call error_handler(E_ERR, 'get_close', errstring, source, &
                      text2='optional arg "dist" is present; we are computing exact distances', &
                      text3='the exhaustive search should find an identical number of locations')
else if (num_close < cnum_close) then
   write(errstring, *) 'get_close (', num_close, ') should not be smaller than exhaustive search (', cnum_close, ')'
   call error_handler(E_ERR, 'get_close', errstring, source, &
                      text2='optional arg "dist" not present; we are returning a superset of close locations', &
                      text3='the exhaustive search should find an equal or lesser number of locations')
endif

! if they do not compare, we have the exhaustive lists here and can print out
! exactly which items in the list differ.

deallocate(cclose_ind, cclose_dist)

end subroutine exhaustive_report

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

subroutine convert_vertical_obs(ens_handle, num, locs, loc_kinds, loc_types, &
                                which_vert, status)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_kinds(:), loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: status(:)

status(:) = 0

end subroutine convert_vertical_obs

!--------------------------------------------------------------------

subroutine convert_vertical_state(ens_handle, num, locs, loc_kinds, loc_indx, &
                                  which_vert, istatus)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_kinds(:)
integer(i8),         intent(in)    :: loc_indx(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus

istatus = 0

end subroutine convert_vertical_state


!----------------------------------------------------------------------------
! end of location/channel/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod

