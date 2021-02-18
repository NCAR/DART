! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module xyz_location_mod

! A subset of the 3d cartesian locations module which only does 'nearest point'
! computations.  It has a limited set of interfaces to initialize the search struct
! and to return the single nearest point.  Does nothing else that a normal
! locations module should do.  

use      types_mod, only : r8, MISSING_R8, MISSING_I, PI, RAD2DEG, DEG2RAD
use  utilities_mod, only : error_handler, E_ERR,           &
                           E_MSG, open_file, close_file, set_output,        &
                           logfileunit, nmlfileunit, find_namelist_in_file, &
                           check_namelist_read, do_nml_file, do_nml_term
use mpi_utilities_mod, only : my_task_id, task_count

implicit none
private

public :: xyz_location_type,         &
          xyz_get_location,          &
          xyz_set_location,          &
          xyz_get_close_type,        &
          xyz_get_close_init,        &
          xyz_get_close_destroy,     &
          xyz_find_nearest,          &
          xyz_find_nearest_N,        &
          xyz_use_great_circle_dist, &
          xyz_get_dist,              &
          xyz_get_ll_location,       &
          xyz_write_location

character(len=*), parameter :: source = 'threed_cartesian/xyz_location_mod.f90'

type xyz_location_type
   private
   real(r8) :: x, y, z
end type xyz_location_type

! This version supports only regularly spaced boxes. it originally had code to
! use an octree division of the space, but finding nearest boxes on each side
! (where there might be multiples) was too complex.  i still think its a good
! idea but need a better data structure to capture the relationships of leaving
! any box via a face and listing all the other boxes that adjoin that face.
!
! the basic idea for octrees was to divide each dim in half until N numbers of filled 
! boxes, or octree reaches some max depth.  give some threshold where you don't
! divide a box with less than N points in it.

! contrast with kD-trees (divide along dimensions, not points), and there are
! two types of octrees - PR (point region) where the regions split at an
! explicit point, vs MX tree where the split is defined to be at the center
! of the region.

! if the underlying geometry is spherical, there will be many many empty boxes 
! if we uniformly divide up space, and worse, existing locations will be 
! clustered in a few boxes.


type box_type
   private
   integer, pointer  :: loc_box(:)           ! (nloc); List of loc indices in boxes
   integer, pointer  :: count(:, :, :)       ! (nx, ny, nz); # of loc in each box
   integer, pointer  :: start(:, :, :)       ! (nx, ny, nz); Start of list of loc in this box
   real(r8)          :: bot_x, top_x         ! extents in x, y, z
   real(r8)          :: bot_y, top_y 
   real(r8)          :: bot_z, top_z 
   real(r8)          :: x_width, y_width, z_width    ! widths of boxes in x,y,z
   real(r8)          :: nboxes_x, nboxes_y, nboxes_z ! based on maxdist how far to search - unused FIXME
end type box_type

! Type to facilitate efficient computation of locations close to a given location
type xyz_get_close_type
   private
   integer           :: num
   real(r8)          :: maxdist
   type(box_type)    :: box
end type xyz_get_close_type

logical, save         :: module_initialized = .false.

real(r8) :: sphere_radius          = -1.0_r8
logical  :: line_of_sight_distance = .true.    ! alternative: great circle

character(len = 512) :: errstring

!-----------------------------------------------------------------
! Namelist with default values

! count of boxes (for box option) in each dim.
integer :: nx         = 20
integer :: ny         = 20
integer :: nz         = 20

namelist /xyz_location_nml/ &
   nx, ny, nz

!-----------------------------------------------------------------

interface xyz_set_location
   module procedure set_location_single
!   module procedure set_location_array
   module procedure set_location_lonlat
end interface xyz_set_location

contains

!----------------------------------------------------------------------------

subroutine initialize_module
 
! things which need doing exactly once.

integer :: iunit, io

if (module_initialized) return

module_initialized = .true.

! Read the namelist entry
call find_namelist_in_file("input.nml", "xyz_location_nml", iunit)
read(iunit, nml = xyz_location_nml, iostat = io)
call check_namelist_read(iunit, io, "xyz_location_nml")

! Write the namelist values to the log file

if(do_nml_file()) write(nmlfileunit, nml=xyz_location_nml)
if(do_nml_term()) write(     *     , nml=xyz_location_nml)

end subroutine initialize_module

!----------------------------------------------------------------------------

subroutine xyz_use_great_circle_dist(radius)
 real(r8), intent(in) :: radius

sphere_radius = radius
line_of_sight_distance = .false.

end subroutine xyz_use_great_circle_dist

!----------------------------------------------------------------------------

function xyz_get_dist(loc1, loc2, type1, kind2)

! returns the distance between 2 locations 

! the 3rd name is a specific type, the 4th is a generic kind.
! These types/kinds are part of
! the interface in case user-code wants to do a more sophisticated distance
! calculation based on the base type or target kind.  In the usual case this
! code still doesn't use the types, but there's an undocumented feature that
! allows you to maintain the original vertical normalization even when
! changing the cutoff distance in the horizontal.  For that to work we
! do need to know the type, and we use the type of loc1 to control it.
! 

type(xyz_location_type), intent(in) :: loc1, loc2
integer, optional,       intent(in) :: type1, kind2
real(r8)                            :: xyz_get_dist

real(r8) :: x_dif, y_dif, z_dif
real(r8) :: mag1, mag2, norm

if ( .not. module_initialized ) call initialize_module

if (line_of_sight_distance) then

   ! straight line in 3d space
x_dif = loc1%x - loc2%x
y_dif = loc1%y - loc2%y
z_dif = loc1%z - loc2%z

xyz_get_dist = sqrt(x_dif * x_dif + y_dif * y_dif + z_dif * z_dif)

   !print *, 'los dist: ', xyz_get_dist

else

   ! great circle distance
   mag1 = sqrt(loc1%x * loc1%x + loc1%y * loc1%y + loc1%z * loc1%z)
   mag2 = sqrt(loc2%x * loc2%x + loc2%y * loc2%y + loc2%z * loc2%z)
   norm = (loc1%x*loc2%x + loc1%y*loc2%y + loc1%z*loc2%z) / (mag1 * mag2)

   xyz_get_dist = sphere_radius * acos(norm)

   !print *, 'mag1, 2, norm, acos, radius: ', mag1, mag2, norm, acos(norm), sphere_radius
   !print *, ' gc dist: ', xyz_get_dist

endif

end function xyz_get_dist

!---------------------------------------------------------------------------

function xyz_get_location(loc)
 
! Given a location type return the x,y,z coordinates

type(xyz_location_type), intent(in) :: loc
real(r8), dimension(3) :: xyz_get_location

if ( .not. module_initialized ) call initialize_module

xyz_get_location(1) = loc%x
xyz_get_location(2) = loc%y
xyz_get_location(3) = loc%z

end function xyz_get_location

!------------------------------------------------------------

subroutine xyz_get_ll_location(loc, radius, lon, lat)

! Given a cartesian x, y, z coordinate relative to an origin at
! the center of the sphere, using a fixed radius specified by the 
! caller, return the corresponding lon, lat location in degrees.

type(xyz_location_type) :: loc
real(r8), intent(in)  :: radius
real(r8), intent(out) :: lon, lat

real(r8) :: rlat, rlon

! don't call this in performance-critical code.
! the arc versions of routines are expensive.

rlat = PI/2.0_r8 - acos(loc%z/radius)
rlon = atan2(loc%y,loc%x)
if (rlon < 0) rlon = rlon + PI*2

lat = rlat * rad2deg
lon = rlon * rad2deg

end subroutine xyz_get_ll_location

!---------------------------------------------------------------------------

function set_location_single(x, y, z)
 
! Puts the x, y, z into a location datatype.

real(r8), intent(in) :: x, y, z
type (xyz_location_type) :: set_location_single

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
type (xyz_location_type) :: set_location_array

if ( .not. module_initialized ) call initialize_module

if (size(list) < 3) then
   write(errstring,*)'requires 3 input values'
   call error_handler(E_ERR, 'set_location', errstring, source)
endif

set_location_array = set_location_single(list(1), list(2), list(3))

end function set_location_array

!----------------------------------------------------------------------------

function set_location_lonlat(lon, lat, height, radius)
 
! given a lon, lat, and radius, compute X,Y,Z and set location

real(r8), intent(in) :: lon, lat, height, radius
type (xyz_location_type) :: set_location_lonlat

real(r8) :: x, y, z
real(r8) :: rlat, rlon, rr

if ( .not. module_initialized ) call initialize_module


rlat = lat * deg2rad
rlon = lon * deg2rad

rr = radius + height

x = rr * cos(rlon) * cos(rlat)
y = rr * sin(rlon) * cos(rlat)
z = rr * sin(rlat)

set_location_lonlat%x = x
set_location_lonlat%y = y
set_location_lonlat%z = z

end function set_location_lonlat

!----------------------------------------------------------------------------

subroutine xyz_get_close_init(gc, maxdist, num, locs)

type(xyz_get_close_type), intent(inout) :: gc
real(r8),             intent(in)    :: maxdist
integer,              intent(in)    :: num
type(xyz_location_type),  intent(in)    :: locs(num)
integer :: i, j, k, cum_start, l
integer :: x_box(num), y_box(num), z_box(num)
integer :: tstart(nx, ny, nz)

if ( .not. module_initialized ) call initialize_module

! set the default value.
gc%maxdist = maxdist
!print *, 'setting maxdist to ', maxdist

! Allocate the storage for the grid dependent boxes
allocate(gc%box%count(nx,ny,nz), gc%box%start(nx,ny,nz))
gc%box%count  = -1
gc%box%start  = -1

! Allocate storage for loc number dependent part
allocate(gc%box%loc_box(num))
gc%box%loc_box(:) = -1

! Set the value of num_loc in the structure
gc%num = num

! If num == 0, no point in going any further.
if (num == 0) return

! Determine where the boxes should be for this set of loc and maxdist
call find_box_ranges(gc, locs, num)

! Begin by computing the number of locations in each box in x,y,z
gc%box%count = 0
do i = 1, num

!print *, i, locs(i)%x, locs(i)%y, locs(i)%z
   x_box(i) = floor((locs(i)%x - gc%box%bot_x) / gc%box%x_width) + 1
   if(x_box(i) > nx) x_box(i) = nx
   if(x_box(i) < 1)  x_box(i) = 1

   y_box(i) = floor((locs(i)%y - gc%box%bot_y) / gc%box%y_width) + 1
   if(y_box(i) > ny) y_box(i) = ny
   if(y_box(i) < 1)  y_box(i) = 1

   z_box(i) = floor((locs(i)%z - gc%box%bot_z) / gc%box%z_width) + 1
   if(z_box(i) > nz) z_box(i) = nz
   if(z_box(i) < 1)  z_box(i) = 1

   gc%box%count(x_box(i), y_box(i), z_box(i)) = gc%box%count(x_box(i), y_box(i), z_box(i)) + 1
!print *, 'adding count to box ', x_box(i), y_box(i), z_box(i), &
!                                 gc%box%count(x_box(i), y_box(i), z_box(i))
end do

! Figure out where storage for each boxes members should begin
cum_start = 1
do i = 1, nx
   do j = 1, ny
      do k = 1, nz
         gc%box%start(i, j, k) = cum_start
         cum_start = cum_start + gc%box%count(i, j, k)
      end do
   end do
end do

! Now we know how many are in each box, get a list of which are in each box
tstart = gc%box%start
do i = 1, num
   gc%box%loc_box(tstart(x_box(i), y_box(i), z_box(i))) = i
   tstart(x_box(i), y_box(i), z_box(i)) = tstart(x_box(i), y_box(i), z_box(i)) + 1
end do

do i = 1, nx
   do j = 1, ny
      do k = 1, nz
!if (gc%box%count(i,j,k) > 0) print *, i,j,k, gc%box%count(i,j,k), gc%box%start(i,j,k)
         do l=1, gc%box%count(i,j,k)
!print *, l, gc%box%loc_box(l)
         enddo
      end do
   end do
end do


end subroutine xyz_get_close_init

!----------------------------------------------------------------------------

subroutine xyz_get_close_destroy(gc)

type(xyz_get_close_type), intent(inout) :: gc

deallocate(gc%box%loc_box, gc%box%count, gc%box%start)

end subroutine xyz_get_close_destroy

!----------------------------------------------------------------------------
!> Finds boundaries for x,y,z boxes.  

subroutine find_box_ranges(gc, locs, num)
 
type(xyz_get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
type(xyz_location_type),  intent(in)    :: locs(num)

logical :: old_out
integer :: i


! FIXME: this space could be very sparse

! find the rectangular prism that encloses all the given locations
! with maxdist boundaries on all sides
gc%box%bot_x = minval(locs(:)%x) - gc%maxdist
gc%box%bot_y = minval(locs(:)%y) - gc%maxdist
gc%box%bot_z = minval(locs(:)%z) - gc%maxdist
  
gc%box%top_x = maxval(locs(:)%x) + gc%maxdist
gc%box%top_y = maxval(locs(:)%y) + gc%maxdist
gc%box%top_z = maxval(locs(:)%z) + gc%maxdist

!! for now, add 1% margin around the space.  
!! should this be a namelist item as a percentage, or an actual distance?
!gc%box%bot_x = gc%box%bot_x - (0.01_r8 * gc%box%bot_x)
!gc%box%bot_y = gc%box%bot_y - (0.01_r8 * gc%box%bot_y)
!gc%box%bot_z = gc%box%bot_z - (0.01_r8 * gc%box%bot_z)
!   
!gc%box%top_x = gc%box%top_x + (0.01_r8 * gc%box%top_x)
!gc%box%top_y = gc%box%top_y + (0.01_r8 * gc%box%top_y)
!gc%box%top_z = gc%box%top_z + (0.01_r8 * gc%box%top_z)

! divide the space into boxes
gc%box%x_width = (gc%box%top_x - gc%box%bot_x) / nx
gc%box%y_width = (gc%box%top_y - gc%box%bot_y) / ny 
gc%box%z_width = (gc%box%top_z - gc%box%bot_z) / nz

! FIXME:  compute a sphere of radius maxdist and see how
! many boxes in x, y, z that would include.  unused at present.
!gc%box%nboxes_x = aint((gc%maxdist + (gc%box%x_width-1)) / gc%box%x_width) 
!gc%box%nboxes_y = aint((gc%maxdist + (gc%box%y_width-1)) / gc%box%y_width) 
!gc%box%nboxes_z = aint((gc%maxdist + (gc%box%z_width-1)) / gc%box%z_width) 

!print *, 'min xyz = ', gc%box%bot_x, gc%box%bot_y, gc%box%bot_z
!print *, 'max xyz = ', gc%box%top_x, gc%box%top_y, gc%box%top_z
!print *, 'wid xyz = ', gc%box%x_width, gc%box%y_width, gc%box%z_width
!print *, 'nbx xyz = ', nx, ny, nz
!!print *, 'nbx xyz = ', gc%box%nboxes_x, gc%box%nboxes_y, gc%box%nboxes_z

end subroutine find_box_ranges

!----------------------------------------------------------------------------
!> find the nearest point in the get close list to the specified point
!> optionally return the exact distance since we have to compute it anyway

subroutine xyz_find_nearest(gc, base_loc, loc_list, nearest, rc, dist)
 type(xyz_get_close_type), intent(in), target  :: gc
 type(xyz_location_type),  intent(in)  :: base_loc
 type(xyz_location_type),  intent(in)  :: loc_list(:)
 integer,              intent(out) :: nearest
 integer,              intent(out) :: rc
 real(r8), optional,   intent(out) :: dist


call find_nearest_boxes(gc, base_loc, loc_list, nearest, rc, dist)

end subroutine xyz_find_nearest

!----------------------------------------------------------------------------
!> find the nearest N points in the get close list to the specified location

subroutine xyz_find_nearest_N(gc, base_loc, loc_list, n_wanted, nearest, n_found, rc, dist)
 type(xyz_get_close_type), intent(in), target  :: gc
 type(xyz_location_type),  intent(in)  :: base_loc
 type(xyz_location_type),  intent(in)  :: loc_list(:)
 integer,                  intent(in)  :: n_wanted
 integer,              intent(out) :: nearest(:)
 integer,              intent(out) :: n_found
 integer,              intent(out) :: rc
 real(r8), optional,   intent(out) :: dist(:)

call find_nearest_n_boxes(gc, base_loc, loc_list, n_wanted, nearest, n_found, rc, dist)

end subroutine xyz_find_nearest_N

!----------------------------------------------------------------------------
!> search all boxes which are potentially close enough

subroutine find_nearest_boxes(gc, base_loc, loc_list, nearest, rc, dist)
 type(xyz_get_close_type), intent(in), target  :: gc
 type(xyz_location_type),  intent(in)  :: base_loc
 type(xyz_location_type),  intent(in)  :: loc_list(:)
 integer,              intent(out) :: nearest
 integer,              intent(out) :: rc
 real(r8), optional,   intent(out) :: dist


integer :: x_box, y_box, z_box, i, j, k, l
integer :: start_x, end_x, start_y, end_y, start_z, end_z
integer ::  n_in_box, st, t_ind(1), ghost, n_found
real(r8) :: this_dist, mindist(1)

! First, set the intent out arguments to a missing value
nearest = -99
t_ind = -99
rc = -1
n_found = 0
mindist = 1e38_r8                  ! something big and positive.
if (present(dist)) dist = 1e38_r8  ! ditto

! the list of locations in the loc_list() argument must be the same
! as the list of locations passed into get_close_init(), so
! gc%num and size(loc_list) better be the same.   if the list changes
! you have to destroy the old gc and init a new one.
if (size(loc_list) /= gc%num) then
   write(errstring,*)'loc_list() array must match one passed to xyz_get_close_init()'
   call error_handler(E_ERR, 'find_nearest_boxes', errstring, source)
endif

! If num == 0, no point in going any further. 
if (gc%num == 0) return

!--------------------------------------------------------------

! Begin by figuring out which box the base loc is in
x_box = floor((base_loc%x - gc%box%bot_x) / gc%box%x_width) + 1
if(x_box > nx) x_box = nx
if(x_box < 1)  x_box = 1
y_box = floor((base_loc%y - gc%box%bot_y) / gc%box%y_width) + 1
if(y_box > ny) y_box = ny
if(y_box < 1)  y_box = 1
z_box = floor((base_loc%z - gc%box%bot_z) / gc%box%z_width) + 1
if(z_box > nz) z_box = nz
if(z_box < 1)  z_box = 1

!print *, 'base_loc box ', x_box, y_box, z_box
!print *, 'nx, ny, nz = ', nx, ny, nz

! If it is not in any box, then it is more than the maxdist away from everybody
if(x_box > nx .or. x_box < 1) return
if(y_box > ny .or. y_box < 1) return
if(z_box > nz .or. z_box < 1) return

!print *, 'good box'


!  First, search all points in the box that contains the base loc
call do_this_box(gc, x_box, y_box, z_box, base_loc, loc_list, 1, t_ind, mindist, n_found)

! if box small enough that no points match, expand search.
! and for now, it is quite possible that points in a neighboring
! box are closer than all points in the current box.  so search
! the neighbors until we are far enough away that we know we
! have the nearest ones.
ghost = 0

10 continue
if (ghost == 0 .or. n_found < 1)  then
   ghost = ghost + 1

   start_x = x_box - ghost
   if (start_x < 1) start_x = 1
   end_x = x_box + ghost
   if (end_x > nx) end_x = nx

   start_y = y_box - ghost
   if (start_y < 1) start_y = 1
   end_y = y_box + ghost
   if (end_y > ny) end_y = ny

   start_z = z_box - ghost
   if (start_z < 1) start_z = 1
   end_z = z_box + ghost
   if (end_z > nz) end_z = nz

   !print *, 'looping from '
   !print *, 'x: ', start_x, end_x
   !print *, 'y: ', start_y, end_y
   !print *, 'z: ', start_z, end_z
   
   ! Next, loop through each box that is close to this box
   do i = start_x, end_x
      do j = start_y, end_y
         do k = start_z, end_z
   
            call do_this_box(gc, i, j, k, base_loc, loc_list, 1, t_ind, mindist, n_found)

         end do
      end do
end do

   if (n_found < 1) then
      ! if we have searched the entire space, punt.
      if (start_x == 1 .and. end_x == nx .and. &
          start_y == 1 .and. end_y == ny .and. &
          start_z == 1 .and. end_z == nz) return
        
      ! repeat search with larger radius of boxes
      goto 10
   endif
endif

nearest = t_ind(1)
rc = 0

! if the caller asked for the actual distance, return it
if (present(dist)) dist = mindist(1)

end subroutine find_nearest_boxes

!----------------------------------------------------------------------------

subroutine find_nearest_n_boxes(gc, base_loc, loc_list, n_wanted, nearest, n_found, rc, dist)
 type(xyz_get_close_type), intent(in), target  :: gc
 type(xyz_location_type),  intent(in)  :: base_loc
 type(xyz_location_type),  intent(in)  :: loc_list(:)
 integer,                  intent(in)  :: n_wanted
 integer,              intent(out) :: nearest(:)
 integer,              intent(out) :: n_found
 integer,              intent(out) :: rc
 real(r8), optional,   intent(out) :: dist(:)

! find the nearest N points in the get close list to the specified point

integer :: x_box, y_box, z_box, i, j, k, l
integer :: start_x, end_x, start_y, end_y, start_z, end_z
integer ::  n_in_box, st, t_ind, ghost
real(r8) :: this_dist, mindist(n_wanted), largest_dist
logical :: been_searched(nx, ny, nz)

! First, set the intent out arguments to a missing value
nearest = -99
n_found = 0
rc = -1
mindist = 1e38_r8                     ! something big and positive.
largest_dist = 1e38_r8                ! ditto
if (present(dist)) dist(:) = 1e38_r8  ! ditto

! make sure they want something - else return or fail?
if (n_wanted < 1) then
   write(errstring,*)'n_wanted must be 1 or larger'
   call error_handler(E_ERR, 'find_nearest_n_boxes', errstring, source)
endif

if (present(dist)) then
  if (size(dist) < n_wanted) then
    write(errstring,*)'dist(:) array size must be equal to or larger than n_wanted'
    call error_handler(E_ERR, 'find_nearest_n_boxes', errstring, source)
  endif
endif

! the list of locations in the loc_list argument must be the same
! as the list of locations passed into get_close_init(), so
! gc%num and size(loc_list better be the same.   if the list changes
! you have to destroy the old gc and init a new one.
if (size(loc_list) /= gc%num) then
   write(errstring,*)'loc_list() array must match one passed to xyz_get_close_init()'
   call error_handler(E_ERR, 'find_nearest_n_boxes', errstring, source)
endif

! If num == 0, no point in going any further. 
if (gc%num == 0) return

!--------------------------------------------------------------

! Begin by figuring out which box the base loc is in
x_box = floor((base_loc%x - gc%box%bot_x) / gc%box%x_width) + 1
if(x_box > nx) x_box = nx
if(x_box < 1)  x_box = 1
y_box = floor((base_loc%y - gc%box%bot_y) / gc%box%y_width) + 1
if(y_box > ny) y_box = ny
if(y_box < 1)  y_box = 1
z_box = floor((base_loc%z - gc%box%bot_z) / gc%box%z_width) + 1
if(z_box > nz) z_box = nz
if(z_box < 1)  z_box = 1

!print *, 'base_loc box ', x_box, y_box, z_box
!print *, 'nx, ny, nz = ', nx, ny, nz

! If it is not in any box, then is this an error?  what do we do about maxdist? FIXME
if(x_box > nx .or. x_box < 1) return
if(y_box > ny .or. y_box < 1) return
if(z_box > nz .or. z_box < 1) return

!print *, 'good box'

been_searched(:,:,:) = .false.

! find the closest N in the box that contains the base loc
call do_this_box(gc, x_box, y_box, z_box, base_loc, loc_list, n_wanted, nearest, mindist, n_found)
been_searched(x_box, y_box, z_box) = .true.

! if box small enough that not enough points match, expand search.
! also, it is quite possible that points in a neighboring
! box are closer than some points in the current box.  so search
! the neighbors until we are far enough away that we know we
! have the nearest ones.  FIXME: be smarter about this.
! have it set a flag for boxes it has already checked?
ghost = 0

10 continue
if (ghost == 0 .or. n_found < n_wanted)  then
   ghost = ghost + 1

   start_x = x_box - ghost
   if (start_x < 1) start_x = 1
   end_x = x_box + ghost
   if (end_x > nx) end_x = nx

   start_y = y_box - ghost
   if (start_y < 1) start_y = 1
   end_y = y_box + ghost
   if (end_y > ny) end_y = ny

   start_z = z_box - ghost
   if (start_z < 1) start_z = 1
   end_z = z_box + ghost
   if (end_z > nz) end_z = nz

   !print *, 'looping from '
   !print *, 'x: ', start_x, end_x
   !print *, 'y: ', start_y, end_y
   !print *, 'z: ', start_z, end_z
   
   ! Next, loop through each box that is close to this box
   do i = start_x, end_x
      do j = start_y, end_y
         do k = start_z, end_z
   
            if (been_searched(i,j,k)) cycle
   
            call do_this_box(gc, i, j, k, base_loc, loc_list, n_wanted, nearest, mindist, n_found)
            been_searched(i,j,k) = .true.
   
         end do
      end do
   end do

   if (n_found < n_wanted) then
      ! if we have searched the entire space, punt.
      if (start_x == 1 .and. end_x == nx .and. &
          start_y == 1 .and. end_y == ny .and. &
          start_z == 1 .and. end_z == nz) return
        
      ! repeat search with larger radius of boxes
      goto 10
   endif
endif

! if they asked for the explicit distances, return them
if (present(dist)) dist(:) = mindist(:)

end subroutine find_nearest_n_boxes

!----------------------------------------------------------------------------

subroutine do_this_box(gc, i, j, k, base_loc, loc_list, n_wanted, nearest, dist, n_found)
 type(xyz_get_close_type), intent(in), target  :: gc
 integer,                  intent(in)          :: i, j, k
 type(xyz_location_type),  intent(in)          :: base_loc
 type(xyz_location_type),  intent(in)          :: loc_list(:)
 integer,                  intent(in)          :: n_wanted
 integer,                  intent(inout)       :: nearest(:)
 real(r8),                 intent(inout)       :: dist(:)
 integer,                  intent(inout)       :: n_found    !< how many we already have
 

integer :: n_in_box, st, l, m, n, t_ind
real(r8) :: this_dist
logical :: this_one_is_wanted

! Box to search is i,j,k
n_in_box = gc%box%count(i, j, k)
st = gc%box%start(i,j,k)
   
! Loop to check how close all locs in the box are; add those that are close
do l = 1, n_in_box
   
   t_ind = gc%box%loc_box(st - 1 + l)
   ! print *, 'l, t_ind = ', l, t_ind
   
   this_dist = xyz_get_dist(base_loc, loc_list(t_ind))
   ! print *, 'this_dist = ', this_dist

   ! if we haven't filled up the number of wanted near points,
   ! or this new distance is smaller than the largest one, we are
   ! going to be adding it.  make these separate tests since if
   ! n_found is 0, you'll get an out-of-bounds on dist(n_found)
   this_one_is_wanted = .false.
   if(n_found < n_wanted) then
      this_one_is_wanted = .true. 
   else if (this_dist < dist(n_found)) then
      this_one_is_wanted = .true.
   endif
   
   if (this_one_is_wanted) then
      if (n_found < n_wanted) n_found = n_found + 1
      do m=1, n_wanted     ! updated value
         if (this_dist >= dist(m)) cycle
         if (m == n_found) then
            nearest(m) = t_ind
            dist(m) = this_dist
            exit
         else 
            do n=n_found, m+1, -1
               nearest(n) = nearest(n-1)
               dist(n) = dist(n-1)
            enddo
            nearest(m) = t_ind
            dist(m) = this_dist
            exit
         endif
      enddo
   endif
end do ! n_in_box

end subroutine do_this_box

!----------------------------------------------------------------------------

subroutine xyz_write_location(loc, buf)
  type(xyz_location_type), intent(in)  :: loc
  character(len=*),        intent(out) :: buf

  if (len(buf) < 18*3) print *, 'buffer too short in xyz_write_location'
  write(buf, '(3F18.6)') loc%x, loc%y, loc%z

end subroutine

!----------------------------------------------------------------------------
! end of location/threed_cartesian/xyz_location_mod.f90
!----------------------------------------------------------------------------

end module xyz_location_mod

