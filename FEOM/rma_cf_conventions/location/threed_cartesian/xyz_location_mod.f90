! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module xyz_location_mod

! A subset of the 3d cartesian locations module which only does 'nearest point'
! computations.  It has a limited set of interfaces to initialize the search struct
! and to return the single nearest point.  Does nothing else that a normal
! locations module should do.  

use      types_mod, only : r8, MISSING_R8, MISSING_I, PI, RAD2DEG, DEG2RAD
use  utilities_mod, only : register_module, error_handler, E_ERR,           &
                           E_MSG, open_file, close_file, set_output,        &
                           logfileunit, nmlfileunit, find_namelist_in_file, &
                           check_namelist_read, do_nml_file, do_nml_term
use mpi_utilities_mod, only : my_task_id, task_count

implicit none
private

public :: xyz_location_type, xyz_get_location, xyz_set_location, &
          xyz_get_close_maxdist_init, xyz_get_close_obs_init, xyz_get_close_type, &
          xyz_find_nearest, xyz_get_close_obs_destroy, xyz_get_dist, xyz_get_ll_location

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type xyz_location_type
   private
   real(r8) :: x, y, z
end type xyz_location_type

! This version supports both regularly spaced boxes, and octree division
! of the space.  for octrees, divide each dim in half until N numbers of filled 
! boxes, or octree reaches some depth?  give some threshold where you don't
! divide a box with less than N points in it?

! contrast with kD-trees (divide along dimensions, not points), and there are
! two types of octrees - PR (point region) where the regions split at an
! explicit point, vs MX tree where the split is defined to be at the center
! of the region.

! if the underlying geometry is spherical, there will be many many empty boxes 
! if we uniformly divide up space, and worse, existing locations will be 
! clustered in a few boxes.


type box_type
   private
   integer, pointer  :: obs_box(:)           ! (nobs); List of obs indices in boxes
   integer, pointer  :: count(:, :, :)       ! (nx, ny, nz); # of obs in each box
   integer, pointer  :: start(:, :, :)       ! (nx, ny, nz); Start of list of obs in this box
   real(r8)          :: bot_x, top_x         ! extents in x, y, z
   real(r8)          :: bot_y, top_y 
   real(r8)          :: bot_z, top_z 
   real(r8)          :: x_width, y_width, z_width    ! widths of boxes in x,y,z
   real(r8)          :: nboxes_x, nboxes_y, nboxes_z ! based on maxdist how far to search
end type box_type

! Type to facilitate efficient computation of observations close to a given location
type xyz_get_close_type
   private
   integer           :: num
   real(r8)          :: maxdist
   type(box_type)    :: box
end type xyz_get_close_type

logical, save         :: module_initialized = .false.

character(len = 512) :: errstring

!real(r8) :: radius     ! used only for converting points on a sphere into x,y,z and back

!-----------------------------------------------------------------
! Namelist with default values

! count of boxes (for box option) in each dim.
integer :: nx               = 20
integer :: ny               = 20
integer :: nz               = 20

! tuning options
integer :: filled           = 10   ! threshold at which you quit splitting
logical :: use_octree       = .false.  ! if false, use regular boxes

! extensible options - these may be useful for tuning the octree 
! integer :: nboxes           = 1000 ! suggestion for max number of nodes
! integer :: maxdepth         = 4    ! suggestion for max tree depth

namelist /xyz_location_nml/ &
   filled, use_octree, &
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

call register_module(source, revision, revdate)
module_initialized = .true.

! Read the namelist entry
call find_namelist_in_file("input.nml", "xyz_location_nml", iunit)
read(iunit, nml = xyz_location_nml, iostat = io)
call check_namelist_read(iunit, io, "xyz_location_nml")

! Write the namelist values to the log file

if(do_nml_file()) write(nmlfileunit, nml=xyz_location_nml)
if(do_nml_term()) write(     *     , nml=xyz_location_nml)

if (filled < 1) then
   write(errstring,*)'filled sets limit for number of points per box.  must be >= 1'
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
endif

! FIXME:
use_octree   = .false.  ! if false, use regular boxes

end subroutine initialize_module

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

if ( .not. module_initialized ) call initialize_module

x_dif = loc1%x - loc2%x
y_dif = loc1%y - loc2%y
z_dif = loc1%z - loc2%z

xyz_get_dist = sqrt(x_dif * x_dif + y_dif * y_dif + z_dif * z_dif)

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

! Given a cartesian x, y, z coordinate relative to the origin
! at the center of the earth, using a fixed radius specified
! by MPAS (in the grid generation step), return the corresponding
! lon, lat location in degrees.

type(xyz_location_type) :: loc
real(r8), intent(in)  :: radius
real(r8), intent(out) :: lon, lat

real(r8) :: rlat, rlon

! right now this is only needed for debugging messages.
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
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
endif

set_location_array = set_location_single(list(1), list(2), list(3))

end function set_location_array

!----------------------------------------------------------------------------

function set_location_lonlat(lon, lat, height, radius)
 
! location semi-independent interface routine
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

!---------------------------------------------------------------------------

subroutine xyz_get_close_obs_init(gc, num, obs)
 
! Initializes part of get_close accelerator that depends on the particular obs

type(xyz_get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
type(xyz_location_type),  intent(in)    :: obs(num)

   call get_close_init_boxes(gc, num, obs)

end subroutine xyz_get_close_obs_init

!----------------------------------------------------------------------------

subroutine get_close_init_boxes(gc, num, obs)
 
! Initializes part of get_close accelerator that depends on the particular obs

type(xyz_get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
type(xyz_location_type),  intent(in)    :: obs(num)

integer :: i, j, k, cum_start, l
integer :: x_box(num), y_box(num), z_box(num)
integer :: tstart(nx, ny, nz)

if ( .not. module_initialized ) call initialize_module

! Allocate storage for obs number dependent part
allocate(gc%box%obs_box(num))
gc%box%obs_box(:) = -1

! Set the value of num_obs in the structure
gc%num = num

! If num == 0, no point in going any further.
if (num == 0) return

! Determine where the boxes should be for this set of obs and maxdist
call find_box_ranges(gc, obs, num)

! Begin by computing the number of observations in each box in x,y,z
gc%box%count = 0
do i = 1, num

!print *, i, obs(i)%x, obs(i)%y, obs(i)%z
   x_box(i) = floor((obs(i)%x - gc%box%bot_x) / gc%box%x_width) + 1
   if(x_box(i) > nx) x_box(i) = nx
   if(x_box(i) < 1)  x_box(i) = 1

   y_box(i) = floor((obs(i)%y - gc%box%bot_y) / gc%box%y_width) + 1
   if(y_box(i) > ny) y_box(i) = ny
   if(y_box(i) < 1)  y_box(i) = 1

   z_box(i) = floor((obs(i)%z - gc%box%bot_z) / gc%box%z_width) + 1
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
   gc%box%obs_box(tstart(x_box(i), y_box(i), z_box(i))) = i
   tstart(x_box(i), y_box(i), z_box(i)) = tstart(x_box(i), y_box(i), z_box(i)) + 1
end do

do i = 1, nx
   do j = 1, ny
      do k = 1, nz
!if (gc%box%count(i,j,k) > 0) print *, i,j,k, gc%box%count(i,j,k), gc%box%start(i,j,k)
         do l=1, gc%box%count(i,j,k)
!print *, l, gc%box%obs_box(l)
         enddo
      end do
   end do
end do


end subroutine get_close_init_boxes

!----------------------------------------------------------------------------

subroutine xyz_get_close_obs_destroy(gc)

type(xyz_get_close_type), intent(inout) :: gc

   call get_close_destroy_boxes(gc)

end subroutine xyz_get_close_obs_destroy

!----------------------------------------------------------------------------

subroutine get_close_destroy_boxes(gc)

type(xyz_get_close_type), intent(inout) :: gc

deallocate(gc%box%obs_box, gc%box%count, gc%box%start)

end subroutine get_close_destroy_boxes

!----------------------------------------------------------------------------

subroutine xyz_get_close_maxdist_init(gc, maxdist)

type(xyz_get_close_type), intent(inout) :: gc
real(r8),             intent(in)    :: maxdist

! set the default value.
gc%maxdist = maxdist
!print *, 'setting maxdist to ', maxdist

if (.not. use_octree) then
   ! Allocate the storage for the grid dependent boxes
   allocate(gc%box%count(nx,ny,nz), gc%box%start(nx,ny,nz))
   gc%box%count  = -1
   gc%box%start  = -1
endif

end subroutine xyz_get_close_maxdist_init

!----------------------------------------------------------------------------

subroutine find_box_ranges(gc, locs, num)
 
! Finds boundaries for x,y,z boxes.  
! FIXME: ways boxes could be divided:
!  - evenly along each axis
!  - octree-like, divide each axis so roughly half the points are
!     on each side of the dividing plane.
!  - about 100 other schemes
  
type(xyz_get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
type(xyz_location_type),  intent(in)    :: locs(num)

! integer :: i

! FIXME: this space could be very sparse

gc%box%bot_x = minval(locs(:)%x)
gc%box%bot_y = minval(locs(:)%y)
gc%box%bot_z = minval(locs(:)%z)
  
gc%box%top_x = maxval(locs(:)%x)
gc%box%top_y = maxval(locs(:)%y)
gc%box%top_z = maxval(locs(:)%z)

!gc%box%bot_x = locs(1)%x
!gc%box%bot_y = locs(1)%y
!gc%box%bot_z = locs(1)%z
!  
!gc%box%top_x = locs(1)%x
!gc%box%top_y = locs(1)%y
!gc%box%top_z = locs(1)%z
!
!do i=2, num
!   gc%box%bot_x = min(gc%box%bot_x, locs(i)%x)
!   gc%box%bot_y = min(gc%box%bot_y, locs(i)%y)
!   gc%box%bot_z = min(gc%box%bot_z, locs(i)%z)
!   
!   gc%box%top_x = max(gc%box%top_x, locs(i)%x)
!   gc%box%top_y = max(gc%box%top_y, locs(i)%y)
!   gc%box%top_z = max(gc%box%top_z, locs(i)%z)
!enddo

gc%box%x_width = (gc%box%top_x - gc%box%bot_x) / nx
gc%box%y_width = (gc%box%top_y - gc%box%bot_y) / ny 
gc%box%z_width = (gc%box%top_z - gc%box%bot_z) / nz

! FIXME:  compute a sphere of radius maxdist and see how
! many boxes in x, y, z that would include.
gc%box%nboxes_x = aint((gc%maxdist + (gc%box%x_width-1)) / gc%box%x_width) 
gc%box%nboxes_y = aint((gc%maxdist + (gc%box%y_width-1)) / gc%box%y_width) 
gc%box%nboxes_z = aint((gc%maxdist + (gc%box%z_width-1)) / gc%box%z_width) 

!print *, 'min xyz = ', gc%box%bot_x, gc%box%bot_y, gc%box%bot_z
!print *, 'max xyz = ', gc%box%top_x, gc%box%top_y, gc%box%top_z
!print *, 'wid xyz = ', gc%box%x_width, gc%box%y_width, gc%box%z_width
!print *, 'nbx xyz = ', nx, ny, nz
!print *, 'nbx xyz = ', gc%box%nboxes_x, gc%box%nboxes_y, gc%box%nboxes_z

end subroutine find_box_ranges

!----------------------------------------------------------------------------

subroutine xyz_find_nearest(gc, base_loc, loc_list, nearest, rc)
 type(xyz_get_close_type), intent(in), target  :: gc
 type(xyz_location_type),  intent(in)  :: base_loc
 type(xyz_location_type),  intent(in)  :: loc_list(:)
 integer,              intent(out) :: nearest
 integer,              intent(out) :: rc

! find the nearest point in the get close list to the specified point

   call find_nearest_boxes(gc, base_loc, loc_list, nearest, rc)

end subroutine xyz_find_nearest

!----------------------------------------------------------------------------

subroutine find_nearest_boxes(gc, base_loc, loc_list, nearest, rc)
 type(xyz_get_close_type), intent(in), target  :: gc
 type(xyz_location_type),  intent(in)  :: base_loc
 type(xyz_location_type),  intent(in)  :: loc_list(:)
 integer,              intent(out) :: nearest
 integer,              intent(out) :: rc

! find the nearest point in the get close list to the specified point

integer :: x_box, y_box, z_box, i, j, k, l
integer :: start_x, end_x, start_y, end_y, start_z, end_z
integer ::  n_in_box, st, t_ind, ghost
real(r8) :: this_dist, dist

! First, set the intent out arguments to a missing value
nearest = -99
rc = -1
dist = 1e38_r8                ! something big and positive.

! the list of locations in the obs() argument must be the same
! as the list of locations passed into get_close_obs_init(), so
! gc%num and size(obs) better be the same.   if the list changes,
! you have to destroy the old gc and init a new one.
if (size(loc_list) /= gc%num) then
   write(errstring,*)'obs() array must match one passed to get_close_obs_init()'
   call error_handler(E_ERR, 'get_close_obs', errstring, source, revision, revdate)
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


! First, search all points in this box.

! Box to search is x_box,y_box,z_box
n_in_box = gc%box%count(x_box,y_box,z_box)
st = gc%box%start(x_box,y_box,z_box)

! find the closest one in this box
do l = 1, n_in_box

   t_ind = gc%box%obs_box(st - 1 + l)
!print *, 'l, t_ind = ', l, t_ind

   this_dist = xyz_get_dist(base_loc, loc_list(t_ind))
!print *, 'this_dist = ', this_dist
   ! If this obs' distance is less than current nearest, it's new nearest
   if(this_dist <= dist) then
      nearest = t_ind
      dist = this_dist
      if (rc < 0) rc = 0
   endif
end do

! if box small enough that no points match, expand search
ghost = 0

10 continue
if (nearest < 0 .or. ghost == 0)  then
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
   
            ! Box to search is i,j,k
            n_in_box = gc%box%count(i, j, k)
            st = gc%box%start(i,j,k)
   
   
            ! Loop to check how close all obs in the box are; add those that are close
            do l = 1, n_in_box
   
               t_ind = gc%box%obs_box(st - 1 + l)
  ! print *, 'l, t_ind = ', l, t_ind
   
               this_dist = xyz_get_dist(base_loc, loc_list(t_ind))
  ! print *, 'this_dist = ', this_dist
               ! If this obs' distance is less than current nearest, it's new nearest
               if(this_dist <= dist) then
                  nearest = t_ind
                  dist = this_dist
                  if (rc < 0) rc = 0
               endif
            end do
         end do
      end do
   end do

   if (nearest < 0) then
      ! if we have searched the entire space, punt.
      if (start_x == 1 .and. end_x == nx .and. &
          start_y == 1 .and. end_y == ny .and. &
          start_z == 1 .and. end_z == nz) return
        
      ! repeat search with larger radius of boxes
      goto 10
   endif
endif

end subroutine find_nearest_boxes

!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
! end of location/threed_cartesian/xyz_location_mod.f90
!----------------------------------------------------------------------------

end module xyz_location_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
