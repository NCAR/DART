! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module location_mod

! Implements location interfaces for a 3d channel in X,Y,Z where X is periodic,
! Y has walls (limited domain), and Z is infinite

use      types_mod, only : r8, MISSING_R8, MISSING_I, PI, RAD2DEG, DEG2RAD
use  utilities_mod, only : register_module, error_handler, E_ERR, ascii_file_format, &
                           nc_check, E_MSG, open_file, close_file, set_output,       &
                           logfileunit, nmlfileunit, find_namelist_in_file,          &
                           check_namelist_read, do_output, do_nml_file,              &
                           do_nml_term, is_longitude_between
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform
use   obs_kind_mod, only : get_num_obs_kinds, get_obs_kind_name
use mpi_utilities_mod, only : my_task_id, task_count

implicit none
private

public :: location_type, get_location, set_location, &
          set_location_missing, is_location_in_region, &
          write_location, read_location, interactive_location, query_location, &
          LocationDims, LocationName, LocationLName, get_close_obs, &
          get_close_maxdist_init, get_close_obs_init, get_close_type, &
          operator(==), operator(/=), get_dist, get_close_obs_destroy, &
          nc_write_location_atts, nc_get_location_varids, nc_write_location, &
          vert_is_height, vert_is_pressure, vert_is_undef, vert_is_level, &
          vert_is_surface, vert_is_scale_height, has_vertical_localization, &
          print_get_close_type, find_nearest, &
          set_vert, get_vert, set_which_vert


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type location_type
   private
   real(r8) :: x, y, z
end type location_type

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


! fortran doesn't let you make arrays of pointers, but you can make a
! derived type containing a pointer, and then make arrays of that derived type.
! i'm sure if i think about this hard enough i'll figure out why this is so,
! but for now i'll just believe the great google which tells me it's this way.
type octree_ptr
   private
   type(octree_type), pointer :: p
end type octree_ptr

type octree_type
   private
   integer             :: count    ! count in this cube, -1 for non-terminal cube
   integer, pointer    :: index(:) ! list of indices in this cube, count long
   type(octree_ptr), allocatable :: children(:,:,:)  ! subcubes
   type(octree_type), pointer    :: parent           ! who made you
   type(location_type) :: llb      ! xyz of lower left bottom
   type(location_type) :: split    ! xyz of split point
   type(location_type) :: urt      ! xyz of upper right top
end type octree_type

type box_type
   private
   integer, pointer  :: loc_box(:)           ! (nloc); List of loc indices in boxes
   integer, pointer  :: count(:, :, :)       ! (nx, ny, nz); # of locs in each box
   integer, pointer  :: start(:, :, :)       ! (nx, ny, nz); Start of list of locs in this box
   real(r8)          :: bot_x, top_x         ! extents in x, y, z
   real(r8)          :: bot_y, top_y 
   real(r8)          :: bot_z, top_z 
   real(r8)          :: x_width, y_width, z_width    ! widths of boxes in x,y,z
   real(r8)          :: nboxes_x, nboxes_y, nboxes_z ! based on maxdist how far to search
end type box_type

! Type to facilitate efficient computation of observations close to a given location
type get_close_type
   private
   integer           :: num
   real(r8)          :: maxdist
   type(box_type)    :: box
   type(octree_type) :: root
end type get_close_type

type(random_seq_type) :: ran_seq
logical               :: ran_seq_init = .false.
logical, save         :: module_initialized = .false.

integer,              parameter :: LocationDims = 3
character(len = 129), parameter :: LocationName = "loc3Dchan"
character(len = 129), parameter :: LocationLName = &
                                   "threed channel locations: x, y, z"

character(len = 512) :: errstring

real(r8) :: radius     ! used only for converting points on a sphere into x,y,z and back

! If maxdist stays the same, don't need to do box distance calculations
integer :: last_maxdist = -1.0

integer :: nx               = 10
integer :: ny               = 10
integer :: nz               = 10

!-----------------------------------------------------------------
! Namelist with default values

logical :: output_box_info  = .false.
integer :: print_box_level  = 0
! tuning options
integer :: nboxes           = 1000 ! suggestion for max number of nodes
integer :: maxdepth         = 4    ! suggestion for max tree depth
integer :: filled           = 10   ! threshold at which you quit splitting
logical :: use_octree       = .false.  ! if false, use regular boxes, true = octree

! Option for verification using exhaustive search
logical :: compare_to_correct = .true.    ! normally false

namelist /location_nml/ &
   filled, nboxes, maxdepth, use_octree, &
   compare_to_correct, output_box_info, print_box_level

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

integer :: iunit, io, i
character(len=129) :: str1

if (module_initialized) return

call register_module(source, revision, revdate)
module_initialized = .true.

! Read the namelist entry
call find_namelist_in_file("input.nml", "location_nml", iunit)
read(iunit, nml = location_nml, iostat = io)
call check_namelist_read(iunit, io, "location_nml")

! Write the namelist values to the log file

if(do_nml_file()) write(nmlfileunit, nml=location_nml)
if(do_nml_term()) write(     *     , nml=location_nml)

if (filled < 1) then
   write(errstring,*)'filled sets limit for number of points per box.  must be >= 1'
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
endif

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
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
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
         'Only "X","Y","Z" are legal attributes to request from location', &
          source, revision, revdate)
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
character(len=129)  :: string1

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
      'Cannot use string buffer with binary format', &
       source, revision, revdate)
endif

! format the location to be human-friendly

! the output can be no longer than this
charlength = 70

if (len(charstring) < charlength) then
   write(errstring, *) 'charstring buffer must be at least ', charlength, ' chars long'
   call error_handler(E_ERR, 'write_location', errstring, source, revision, revdate)
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

character(len=8) :: header

if ( .not. module_initialized ) call initialize_module

if (ascii_file_format(fform)) then
   read(locfile, '(A8)' ) header
   if(header /= 'loc3Dchan') then
         write(errstring,*)'Expected location header "loc3Dchan" in input file, got ', header
      call error_handler(E_ERR, 'read_location', errstring, source, revision, revdate)
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

function nc_write_location_atts( ncFileID, fname, ObsNumDimID ) result (ierr)
 
! Writes the "location module" -specific attributes to a netCDF file.

use typeSizes
use netcdf

integer,          intent(in) :: ncFileID     ! handle to the netcdf file
character(len=*), intent(in) :: fname        ! file name (for printing purposes)
integer,          intent(in) :: ObsNumDimID  ! handle to the dimension that grows
integer                      :: ierr

integer :: LocDimID
integer :: VarID

if ( .not. module_initialized ) call initialize_module

ierr = -1 ! assume things will fail ...

! define the rank/dimension of the location information
call nc_check(nf90_def_dim(ncid=ncFileID, name='location', len=LocationDims, &
       dimid = LocDimID), 'nc_write_location_atts', 'def_dim:location '//trim(fname))

! Define the location variable and attributes

call nc_check(nf90_def_var(ncid=ncFileID, name='location', xtype=nf90_double, &
          dimids=(/ LocDimID, ObsNumDimID /), varid=VarID), &
            'nc_write_location_atts', 'location:def_var')

call nc_check(nf90_put_att(ncFileID, VarID, 'description', &
        'location coordinates'), 'nc_write_location_atts', 'location:description')
call nc_check(nf90_put_att(ncFileID, VarID, 'location_type', &
        trim(LocationName)), 'nc_write_location_atts', 'location:location_type')
call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', &
        trim(LocationLName)), 'nc_write_location_atts', 'location:long_name')
call nc_check(nf90_put_att(ncFileID, VarID, 'storage_order',     &
        'X Y Z'), 'nc_write_location_atts', 'location:storage_order')
call nc_check(nf90_put_att(ncFileID, VarID, 'units',     &
        'X Y Z'), 'nc_write_location_atts', 'location:units')

ierr = 0

end function nc_write_location_atts

!----------------------------------------------------------------------------

subroutine nc_get_location_varids( ncFileID, fname, LocationVarID, WhichVertVarID )

! Return the LocationVarID and WhichVertVarID variables from a given netCDF file.
!
! ncFileId         the netcdf file descriptor
! fname            the name of the netcdf file (for error messages only)
! LocationVarID    the integer ID of the 'location' variable in the netCDF file
! WhichVertVarID   the integer ID of the 'which_vert' variable in the netCDF file

use typeSizes
use netcdf

integer,          intent(in)  :: ncFileID   ! handle to the netcdf file
character(len=*), intent(in)  :: fname      ! file name (for printing purposes)
integer,          intent(out) :: LocationVarID, WhichVertVarID

if ( .not. module_initialized ) call initialize_module

call nc_check(nf90_inq_varid(ncFileID, 'location', varid=LocationVarID), &
          'nc_get_location_varids', 'inq_varid:location '//trim(fname))

WhichVertVarID = -1

end subroutine nc_get_location_varids

!----------------------------------------------------------------------------

subroutine nc_write_location(ncFileID, LocationVarID, loc, obsindex, WhichVertVarID)
 
! Writes a SINGLE location to the specified netCDF variable and file.
! The LocationVarID and WhichVertVarID must be the values returned from
! the nc_get_location_varids call.

use typeSizes
use netcdf

integer,             intent(in) :: ncFileID, LocationVarID
type(location_type), intent(in) :: loc
integer,             intent(in) :: obsindex
integer,             intent(in) :: WhichVertVarID

real(r8), dimension(LocationDims) :: locations
integer,  dimension(1) :: intval

if ( .not. module_initialized ) call initialize_module

locations = get_location( loc )

call nc_check(nf90_put_var(ncFileID, LocationVarId, locations, &
          start=(/ 1, obsindex /), count=(/ LocationDims, 1 /) ), &
            'nc_write_location', 'put_var:location')

end subroutine nc_write_location

!----------------------------------------------------------------------------

subroutine get_close_obs_init(gc, num, obs)
 
! Initializes part of get_close accelerator that depends on the particular obs

type(get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
type(location_type),  intent(in)    :: obs(num)

if (use_octree) then
   call get_close_init_otree(gc, num, obs)
else
   call get_close_init_boxes(gc, num, obs)
endif

end subroutine get_close_obs_init

!----------------------------------------------------------------------------

subroutine get_close_init_boxes(gc, num, locs)
 
! Initializes part of get_close accelerator that depends on the particular loc

type(get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
type(location_type),  intent(in)    :: locs(num)

integer :: i, j, k, cum_start, l
integer :: x_box(num), y_box(num), z_box(num)
integer :: tstart(nx, ny, nz)

if ( .not. module_initialized ) call initialize_module

! Allocate storage for obs number dependent part
allocate(gc%box%loc_box(num))
gc%box%loc_box(:) = -1

! Set the value of num_locs in the structure
gc%num = num

! If num == 0, no point in going any further.
if (num == 0) return

! FIXME: compute nx, ny, nz from nboxes?  or put in namelist
nx = nint(real(nboxes, r8)**0.33333)   ! roughly cube root
ny = nint(real(nboxes, r8)**0.33333)   ! roughly cube root
nz = nint(real(nboxes, r8) / real(nx * ny, r8))  ! whatever is left

! Determine where the boxes should be for this set of locs and maxdist
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
if (gc%box%count(i,j,k) > 0) print *, i,j,k, gc%box%count(i,j,k), gc%box%start(i,j,k)
         do l=1, gc%box%count(i,j,k)
!print *, l, gc%box%loc_box(l)
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

end subroutine get_close_init_boxes

!----------------------------------------------------------------------------

subroutine get_close_init_otree(gc, num, locs)
 
! Octree version

! Initializes part of get_close accelerator that depends on the particular locs

type(get_close_type), intent(inout), target :: gc
integer,              intent(in)    :: num
type(location_type),  intent(in)    :: locs(num)

integer :: i, j, k, n

type(octree_type), pointer :: r, c
real(r8) :: xl, xu, yl, yu, zl, zu

if ( .not. module_initialized ) call initialize_module

r => gc%root

! Set the value of num_locs in the structure
gc%num = num

! If num == 0, no point in going any further.
if (num == 0) return

r%count = num

! need to include space outside the limits of the initialization set,
! so points outside boundary but closer than maxdist will still match.
r%llb = set_location(minval(locs(:)%x)-gc%maxdist, &
                     minval(locs(:)%y)-gc%maxdist, &
                     minval(locs(:)%z)-gc%maxdist)
r%urt = set_location(maxval(locs(:)%x)+gc%maxdist, &
                     maxval(locs(:)%y)+gc%maxdist, &
                     maxval(locs(:)%z)+gc%maxdist)
! for now, split is midpoint in each dim.  this does NOT have to be the case.
r%split = set_location((r%urt%x + r%llb%x) / 2.0_r8, &
                       (r%urt%y + r%llb%y) / 2.0_r8, &
                       (r%urt%z + r%llb%z) / 2.0_r8)

! initially everyone is in the original list
allocate(r%index(num))
do i=1,num
   r%index(i) = i
enddo

! recursion starts here
if (r%count > filled) then
   call split_tree(r)
   call move_to_children(r, locs)
endif

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

end subroutine get_close_init_otree

!----------------------------------------------------------------------------

subroutine get_close_obs_destroy(gc)

type(get_close_type), intent(inout) :: gc

if (use_octree) then
   call get_close_destroy_otree(gc)
else
   call get_close_destroy_boxes(gc)
endif

end subroutine get_close_obs_destroy

!----------------------------------------------------------------------------

subroutine get_close_destroy_boxes(gc)

type(get_close_type), intent(inout) :: gc

deallocate(gc%box%loc_box, gc%box%count, gc%box%start)

end subroutine get_close_destroy_boxes

!----------------------------------------------------------------------------

subroutine get_close_destroy_otree(gc)

type(get_close_type), intent(inout) :: gc

! FIXME: do a depth-first search and deallocate
! the index() arrays, then deallocate the octree structs
! one by one.

! gc%root -> children until unallocated

end subroutine get_close_destroy_otree

!----------------------------------------------------------------------------

subroutine get_close_maxdist_init(gc, maxdist)

type(get_close_type), intent(inout) :: gc
real(r8),             intent(in)    :: maxdist

character(len=129) :: str1
integer :: i

! set the default value.
gc%maxdist = maxdist
!print *, 'setting maxdist to ', maxdist

if (.not. use_octree) then
   ! Allocate the storage for the grid dependent boxes
   allocate(gc%box%count(nx,ny,nz), gc%box%start(nx,ny,nz))
   gc%box%count  = -1
   gc%box%start  = -1
endif

end subroutine get_close_maxdist_init

!----------------------------------------------------------------------------

subroutine get_close_obs(gc, base_obs_loc, base_obs_type, obs, obs_kind, &
   num_close, close_ind, dist)

! FIXME: these work on any locations. the names of these args should be:  
!   gc, base_loc, base_type, locs, locs_kinds, ... 

type(get_close_type), intent(in)  :: gc
type(location_type),  intent(in)  :: base_obs_loc,  obs(:)
integer,              intent(in)  :: base_obs_type, obs_kind(:)
integer,              intent(out) :: num_close, close_ind(:)
real(r8), optional,   intent(out) :: dist(:)

if (use_octree) then
   call get_close_otree(gc, base_obs_loc, base_obs_type, obs, obs_kind, &
                        num_close, close_ind, dist)
else
   call get_close_boxes(gc, base_obs_loc, base_obs_type, obs, obs_kind, &
                        num_close, close_ind, dist)
endif

end subroutine get_close_obs

!----------------------------------------------------------------------------

subroutine get_close_boxes(gc, base_loc, base_type, loc, loc_kind, &
   num_close, close_ind, dist)

type(get_close_type), intent(in)  :: gc
type(location_type),  intent(in)  :: base_loc,  loc(:)
integer,              intent(in)  :: base_type, loc_kind(:)
integer,              intent(out) :: num_close, close_ind(:)
real(r8), optional,   intent(out) :: dist(:)

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

! the list of locations in the loc() argument must be the same
! as the list of locations passed into get_close_obs_init(), so
! gc%num and size(loc) better be the same.   if the list changes,
! you have to destroy the old gc and init a new one.
if (size(loc) /= gc%num) then
   write(errstring,*)'loc() array must match one passed to get_close_obs_init()'
   call error_handler(E_ERR, 'get_close_boxes', errstring, source, revision, revdate)
endif

! If num == 0, no point in going any further. 
if (gc%num == 0) return

this_maxdist = gc%maxdist


! For validation, it is useful to be able to compare against exact
! exhaustive search
if(compare_to_correct) then
   call exhaustive_collect(gc, base_loc, loc, &
                           cnum_close, cclose_ind, cdist)
endif


! Begin by figuring out which box the base loc is in
x_box = floor((base_loc%x - gc%box%bot_x) / gc%box%x_width) + 1
y_box = floor((base_loc%y - gc%box%bot_y) / gc%box%y_width) + 1
z_box = floor((base_loc%z - gc%box%bot_z) / gc%box%z_width) + 1

! If it is not in any box, then it is more than the maxdist away from everybody
if(x_box > nx .or. x_box < 1 .or. x_box < 0) return
if(y_box > ny .or. y_box < 1 .or. y_box < 0) return
if(z_box > nz .or. z_box < 1 .or. z_box < 0) return

! figure out how many boxes need searching
! FIXME: if we support a variable maxdist, nboxes_X will need to
! be computed on the fly here instead of precomputed at init time.
!print *, 'nboxes x, y, z = ', gc%box%nboxes_x, gc%box%nboxes_y, gc%box%nboxes_z
!print *, 'base_loc in box ', x_box, y_box, z_box

start_x = x_box - gc%box%nboxes_x
if (start_x < 1) start_x = 1
end_x = x_box + gc%box%nboxes_x
if (end_x > nx) end_x = nx

start_y = y_box - gc%box%nboxes_y
if (start_y < 1) start_y = 1
end_y = y_box + gc%box%nboxes_y
if (end_y > ny) end_y = ny

start_z = z_box - gc%box%nboxes_z
if (start_z < 1) start_z = 1
end_z = z_box + gc%box%nboxes_z
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


         ! Loop to check how close all loc in the box are; add those that are close
         do l = 1, n_in_box

            t_ind = gc%box%loc_box(st - 1 + l)
!print *, 'l, t_ind = ', l, t_ind

            ! Only compute distance if dist is present
            if(present(dist)) then
               this_dist = get_dist(base_loc, loc(t_ind))
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


end subroutine get_close_boxes

!----------------------------------------------------------------------------

subroutine get_close_otree(gc, base_loc, base_type, loc, loc_kind, &
   num_close, close_ind, dist)

type(get_close_type), intent(in), target  :: gc
type(location_type),  intent(in)  :: base_loc,  loc(:)
integer,              intent(in)  :: base_type, loc_kind(:)
integer,              intent(out) :: num_close, close_ind(:)
real(r8), optional,   intent(out) :: dist(:)

! If dist is NOT present, just find everybody in a box, put them in the list,
! but don't compute any distances

type(octree_type), pointer :: r, c

integer :: x_box, y_box, z_box, i, j, k, l, n
integer :: start_x, end_x, start_y, end_y, start_z, end_z
integer ::  n_in_box, st, t_ind
real(r8) :: this_dist

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

! the list of locations in the loc() argument must be the same
! as the list of locations passed into get_close_obs_init(), so
! gc%num and size(loc) better be the same.   if the list changes,
! you have to destroy the old gc and init a new one.
if (size(loc) /= gc%num) then
   write(errstring,*)'loc() array must match one passed to get_close_obs_init()'
   call error_handler(E_ERR, 'get_close_otree', errstring, source, revision, revdate)
endif

! If num == 0, no point in going any further. 
if (gc%num == 0) return


! For validation, it is useful to be able to compare against exact
! exhaustive search
if(compare_to_correct) then
   call exhaustive_collect(gc, base_loc, loc, &
                           cnum_close, cclose_ind, cdist)
endif


! revised plan:
! find min/max extents in each dim, +/- maxdist from base
! and only descend lower in a branch if the region intersects
! the range that this branch covers.  maxdist doesn't have to
! be a constant in this case - it can be specified per search.

call collect_nearby(gc%root, base_loc, gc%maxdist, loc, base_type, loc_kind, &
                    num_close, close_ind, dist)


! Verify by comparing to exhaustive search
if(compare_to_correct) then
   call exhaustive_report(cnum_close, num_close, cclose_ind, close_ind, cdist, dist)
endif


end subroutine get_close_otree

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

logical :: old_out


! FIXME: this space could be very sparse

gc%box%bot_x = minval(locs(:)%x)
gc%box%bot_y = minval(locs(:)%y)
gc%box%bot_z = minval(locs(:)%z)

gc%box%top_x = maxval(locs(:)%x)
gc%box%top_y = maxval(locs(:)%y)
gc%box%top_z = maxval(locs(:)%z)

gc%box%x_width = (gc%box%top_x - gc%box%bot_x) / nx
gc%box%y_width = (gc%box%top_y - gc%box%bot_y) / ny 
gc%box%z_width = (gc%box%top_z - gc%box%bot_z) / nz

! FIXME:  compute a sphere of radius maxdist and see how
! many boxes in x, y, z that would include.
gc%box%nboxes_x = aint((gc%maxdist + (gc%box%x_width-1)) / gc%box%x_width) 
gc%box%nboxes_y = aint((gc%maxdist + (gc%box%y_width-1)) / gc%box%y_width) 
gc%box%nboxes_z = aint((gc%maxdist + (gc%box%z_width-1)) / gc%box%z_width) 


!if(compare_to_correct) then
!   old_out = do_output()
!   call set_output(.true.)
!   write(errstring, *) 'x bot, top, width, nboxes ', gc%box%bot_x, gc%box%top_x, gc%box%x_width, gc%box%nboxes_x
!   call error_handler(E_MSG, 'find_box_ranges', errstring)
!   write(errstring, *) 'y bot, top, width, nboxes ', gc%box%bot_y, gc%box%top_y, gc%box%y_width, gc%box%nboxes_y
!   call error_handler(E_MSG, 'find_box_ranges', errstring)
!   write(errstring, *) 'z bot, top, width, nboxes ', gc%box%bot_z, gc%box%top_z, gc%box%z_width, gc%box%nboxes_z
!   call error_handler(E_MSG, 'find_box_ranges', errstring)
!   call set_output(old_out)
!endif

end subroutine find_box_ranges

!----------------------------------------------------------------------------

recursive subroutine split_tree(r)
 type(octree_type), pointer :: r

integer :: i, j, k
real(r8) :: xl, xu, yl, yu, zl, zu
type(octree_type), pointer :: c


allocate(r%children(2,2,2))

do i=1,2
   if (i == 1) then
      xl = r%llb%x
      xu = r%split%x
   else
      xl = r%split%x
      xu = r%urt%x
   endif
   do j=1,2
      if (j == 1) then
         yl = r%llb%y
         yu = r%split%y
      else
         yl = r%split%y
         yu = r%urt%y
      endif
      do k=1,2
         if (k == 1) then
            zl = r%llb%z
            zu = r%split%z
         else
            zl = r%split%z
            zu = r%urt%z
         endif

         allocate(r%children(i,j,k)%p)
         c => r%children(i,j,k)%p
         c%count = 0
         nullify(c%index)
         c%llb = set_location(xl, yl, zl)
         c%urt = set_location(xu, yu, zu)
         c%split = set_location((c%urt%x + c%llb%x) / 2.0_r8, &
                                (c%urt%y + c%llb%y) / 2.0_r8, &
                                (c%urt%z + c%llb%z) / 2.0_r8)

      enddo
   enddo
enddo

end subroutine split_tree

!----------------------------------------------------------------------------

recursive subroutine move_to_children(r, locs)
 type(octree_type),   pointer    :: r
 type(location_type), intent(in) :: locs(:)

integer :: n, i, j, k
type(octree_type), pointer :: c

! child counts need to be allocated to a reasonable size
do i=1,2
   do j=1,2
      do k=1,2
         c => r%children(i,j,k)%p
         ! FIXME: this is too big now
         allocate(c%index(r%count))
      enddo
   enddo
enddo

do n=1, r%count
   if (locs(r%index(n))%x < r%split%x) then
      i = 1
   else
      i = 2
   endif

   if (locs(r%index(n))%y < r%split%y) then
      j = 1
   else
      j = 2
   endif

   if (locs(r%index(n))%z < r%split%z) then
      k = 1
   else
      k = 2
   endif

   c => r%children(i,j,k)%p
   c%count = c%count + 1
   c%index(c%count) = r%index(n)

enddo

r%count = -1
deallocate(r%index)
nullify(r%index)

do i=1,2
   do j=1,2
      do k=1,2
         c => r%children(i,j,k)%p
         c%parent => r
         if (c%count > filled) then
            call split_tree(c)
            call move_to_children(c, locs)
         endif
      enddo
   enddo
enddo


end subroutine move_to_children

!----------------------------------------------------------------------------

subroutine find_nearest(gc, base_loc, loc_list, nearest, rc)
 type(get_close_type), intent(in), target  :: gc
 type(location_type),  intent(in)  :: base_loc
 type(location_type),  intent(in)  :: loc_list(:)
 integer,              intent(out) :: nearest
 integer,              intent(out) :: rc

! find the nearest point in the get close list to the specified point

if (use_octree) then
   call find_nearest_otree(gc, base_loc, loc_list, nearest, rc)
else
   call find_nearest_boxes(gc, base_loc, loc_list, nearest, rc)
endif

end subroutine find_nearest

!----------------------------------------------------------------------------

subroutine find_nearest_otree(gc, base_loc, loc_list, nearest, rc)
 type(get_close_type), intent(in), target  :: gc
 type(location_type),  intent(in)  :: base_loc
 type(location_type),  intent(in)  :: loc_list(:)
 integer,              intent(out) :: nearest
 integer,              intent(out) :: rc

! find the nearest point in the get close list to the specified point

integer :: n, i, j, k, nindex
real(r8) :: ndist, dist
type(octree_type), pointer :: r, c

! until children is null, descend based on the side of the split point
! the base loc is on

! arrgghh - what do we do about points outside the original cube?
! there is still a closest point defined somewhere...  do we use
! maxdist here at all?  (i'm thinking no.)

! FIXME: we cannot reject points outside the cube if we're just
! looking for any nearest point.  if we say nearest has to be
! within maxdist, and if we define the root cube to have +/- maxdist
! around all edges, then this is ok.

! FIXME: we must add maxdist and search that region for the
! closest point - the closest point might be in the neighboring block.
! (maybe maxdist isn't the value to use?)

r => gc%root
if (.not. is_location_in_region(base_loc, r%llb, r%urt)) then
   nearest = -1
   rc = -1
   return
endif

r = find_child(r, base_loc)

ndist = 1e38_r8
nindex = -1
do i=1, r%count
   n = r%index(i)
   dist =  get_dist(base_loc, loc_list(n))
   if (dist < ndist) then
      nindex = n
      ndist = dist
   endif
enddo

if (nindex < 0) then 
   nearest = -1
   rc = -1
   return
endif
   
nearest = nindex
rc = 0

end subroutine find_nearest_otree

!----------------------------------------------------------------------------

subroutine find_nearest_boxes(gc, base_loc, loc_list, nearest, rc)
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
! as the list of locations passed into get_close_obs_init(), so
! gc%num and size(loc) better be the same.   if the list changes,
! you have to destroy the old gc and init a new one.
if (size(loc_list) /= gc%num) then
   write(errstring,*)'loc() array must match one passed to get_close_obs_init()'
   call error_handler(E_ERR, 'find_nearest_boxes', errstring, source, revision, revdate)
endif

! If num == 0, no point in going any further. 
if (gc%num == 0) return

!--------------------------------------------------------------

! Begin by figuring out which box the base loc is in
x_box = floor((base_loc%x - gc%box%bot_x) / gc%box%x_width) + 1
y_box = floor((base_loc%y - gc%box%bot_y) / gc%box%y_width) + 1
z_box = floor((base_loc%z - gc%box%bot_z) / gc%box%z_width) + 1

! If it is not in any box, then it is more than the maxdist away from everybody
if(x_box > nx .or. x_box < 1 .or. x_box < 0) return
if(y_box > ny .or. y_box < 1 .or. y_box < 0) return
if(z_box > nz .or. z_box < 1 .or. z_box < 0) return

! figure out how many boxes need searching
! FIXME: if we support a variable maxdist, nboxes_X will need to
! be computed on the fly here instead of precomputed at init time.
!print *, 'nboxes x, y, z = ', gc%box%nboxes_x, gc%box%nboxes_y, gc%box%nboxes_z
!print *, 'base_loc in box ', x_box, y_box, z_box

start_x = x_box - gc%box%nboxes_x
if (start_x < 1) start_x = 1
end_x = x_box + gc%box%nboxes_x
if (end_x > nx) end_x = nx

start_y = y_box - gc%box%nboxes_y
if (start_y < 1) start_y = 1
end_y = y_box + gc%box%nboxes_y
if (end_y > ny) end_y = ny

start_z = z_box - gc%box%nboxes_z
if (start_z < 1) start_z = 1
end_z = z_box + gc%box%nboxes_z
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


         ! Loop to check how close all loc in the box are; add those that are close
         do l = 1, n_in_box

            t_ind = gc%box%loc_box(st - 1 + l)
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

end subroutine find_nearest_boxes

!----------------------------------------------------------------------------

function find_child(r, base_loc)
 type(octree_type), pointer :: r
 type(location_type)        :: base_loc
 type(octree_type), pointer :: find_child

integer :: i, j, k

! until we get to the leaves, descend.
do while(r%count < 0)
! print *, 'rcount = ', r%count
! print *, 'llb = ', r%llb%x, r%llb%y, r%llb%z
! print *, 'spl = ', r%split%x, r%split%y, r%split%z
! print *, 'urt = ', r%urt%x, r%urt%y, r%urt%z
! print *, 'loc = ', base_loc%x, base_loc%y, base_loc%z
   if (base_loc%x < r%split%x) then
      i = 1;
   else
      i = 2;
   endif
   if (base_loc%y < r%split%y) then
      j = 1;
   else
      j = 2;
   endif
   if (base_loc%z < r%split%z) then
      k = 1;
   else
      k = 2;
   endif
! print *, 'ijk= ', i,j,k

   if (.not. allocated(r%children)) then
 print *, 'children array not allocated, stop'
 stop
   endif
   r => r%children(i,j,k)%p
enddo

find_child => r

end function find_child

!----------------------------------------------------------------------------

subroutine print_get_close_type(gc, amount)
 
! print out debugging statistics, or optionally print out a full
! dump from all mpi tasks in a format that can be plotted with matlab.

type(get_close_type), intent(in), target :: gc
integer, intent(in), optional            :: amount

if (use_octree) then
   call print_get_close_otree(gc, amount)
else
   call print_get_close_boxes(gc, amount)
endif

end subroutine print_get_close_type

!----------------------------------------------------------------------------

subroutine print_get_close_otree(gc, amount)
 
! print out debugging statistics, or optionally print out a full
! dump from all mpi tasks in a format that can be plotted with matlab.

type(get_close_type), intent(in), target :: gc
integer, intent(in), optional            :: amount

type(octree_type), pointer :: r

r => gc%root
call print_tree(r, 0)

end subroutine print_get_close_otree

!----------------------------------------------------------------------------

subroutine print_get_close_boxes(gc, amount)
 
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
   write(errstring, "(A,3(F12.6))") ' x_box: bot, top, width = ', gc%box%bot_x, gc%box%top_x, gc%box%x_width
   call error_handler(E_MSG, 'loc', errstring)
   write(errstring, "(A,3(F12.6))") ' y_box: bot, top, width = ', gc%box%bot_y, gc%box%top_y, gc%box%y_width
   call error_handler(E_MSG, 'loc', errstring)
   write(errstring, "(A,3(F12.6))") ' z_box: bot, top, width = ', gc%box%bot_z, gc%box%top_z, gc%box%z_width
   call error_handler(E_MSG, 'loc', errstring)

endif

! this one can be very large.   print only the first nth unless
! instructed otherwise.  (print n+1 because 1 more value fits on
! the line because it prints ( i ) and not ( i, j ) like the others.)
if (associated(gc%box%loc_box)) then
   i = size(gc%box%loc_box,1)
   if (i/= gc%num) then
      write(errstring,*) ' warning: size of loc_box incorrect, nlocs, i =', gc%num, i
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
   if (howmuch > 1) then
      ! DEBUG
      write(errstring,"(A,I8,A,36(I8,1X))") ' loc_box(',i,') =', gc%box%loc_box(1:min(i,36))  ! (nlocs)
      !write(errstring,*) ' loc_box(',i,') =', gc%box%loc_box    ! (nlocs)
      call error_handler(E_MSG, 'locations_mod', errstring)
   else if(howmuch > 0) then
      write(errstring,*) ' loc_box(',i,') =', gc%box%loc_box(1:min(i,sample+1))
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
if (associated(gc%box%start)) then
   i = size(gc%box%start,1)
   j = size(gc%box%start,2)
   k = size(gc%box%start,3)
   if ((i /= nx) .or. (j /= ny) .or. (k /= nz)) then
      write(errstring,*) ' warning: size of start incorrect, nx, ny, nz, i, j, k =', nx, ny, nz, i, j, k
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
   if (howmuch > 1) then
      write(errstring,*) ' start(',i,j,k,') ='              ! (nx, ny, nz)
      call error_handler(E_MSG, 'locations_mod', errstring)
      do l=1, j
         write(errstring,"(36(I8,1X))") gc%box%start(1:min(i,36), l, 1)
         call error_handler(E_MSG, 'locations_mod', errstring)
      enddo
   else if (howmuch > 0) then
      write(errstring,*) ' start(',i,j,k,') =', gc%box%start(1:min(i,sample), 1, 1)
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
if (associated(gc%box%count)) then
   i = size(gc%box%count,1)
   j = size(gc%box%count,2)
   k = size(gc%box%count,3)
   if ((i /= nx) .or. (j /= ny) .or. (k /= nz)) then
      write(errstring,*) ' warning: size of count incorrect, nx, ny, nz, i, j, k =', nx, ny, nz, i, j, k
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
   if (howmuch > 1) then
      write(errstring,*) ' count(',i,j,k,') ='              ! (nx, ny, nz)
      call error_handler(E_MSG, 'locations_mod', errstring)
      do l=1, j
         write(errstring,"(36(I8,1X))") gc%box%count(1:min(i,36), l, 1) 
         call error_handler(E_MSG, 'locations_mod', errstring)
      enddo
   else if (howmuch > 0) then
      write(errstring,*) ' count(',i,j,k,') =', gc%box%count(1:min(i,sample), 1, 1)
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
         first = gc%box%start(i, j, k)
         do l=1, gc%box%count(i, j, k)
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
      x_cen = gc%box%bot_x + ((i-1)*gc%box%x_width) + (gc%box%x_width/2.0)
      write(funit, '(A,I2,A,I4,A,F12.9,A)') 'xlocs(', i, ',', mytask+1, ') = ',  x_cen, ';'
   endif
   do j=1, ny
      if (howmuch == -8 .and. i==1) then
         y_cen = gc%box%bot_y + ((j-1)*gc%box%y_width) + (gc%box%y_width/2.0)
         write(funit, '(A,I2,A,I4,A,F12.9,A)') 'ylocs(', j, ',', mytask+1, ') = ',  y_cen, ';'
      endif
      do k=1, nz
         if (howmuch == -8 .and. i==1) then
            z_cen = gc%box%bot_z + ((j-1)*gc%box%z_width) + (gc%box%z_width/2.0)
            write(funit, '(A,I2,A,I4,A,F12.9,A)') 'zlocs(', k, ',', mytask+1, ') = ',  z_cen, ';'
         endif
         if (gc%box%count(i, j, k) > 0) then
            nfull = nfull + 1
            total = total + gc%box%count(i, j, k)
            if (gc%box%count(i, j, k) > maxcount) then
               maxcount = gc%box%count(i, j, k)
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
                                   ',', mytask+1, ') = ', gc%box%count(i, j, k), ';'
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

end subroutine print_get_close_boxes

!----------------------------------------------------------------------------

recursive subroutine print_tree(r, l)
 type(octree_type), pointer :: r
 integer, intent(in) :: l

integer :: i, j, k
   
   print *, '--start of block--'
   print *, 'level = ', l
   print *, 'count = ', r%count
   print *, 'llb   = ', r%llb%x, r%llb%y, r%llb%z
   print *, 'split = ', r%split%x, r%split%y, r%split%z
   print *, 'urt   = ', r%urt%x, r%urt%y, r%urt%z
   if (associated(r%index)) then
      do i=1, r%count
         print *, i, r%index(i)
      enddo
   endif
   print *, 'parent associated = ', associated(r%parent)
   print *, 'child allocated = ', allocated(r%children)
   if (allocated(r%children)) then
      do i=1,2
         do j=1,2
            do k=1,2
               print *, i,j,k
               call print_tree(r%children(i,j,k)%p, l+1)
            enddo
         enddo
      enddo
   endif
   print *, '--end of block--'

end subroutine print_tree

!----------------------------------------------------------------------------

function region_intersects_block(loc, maxdist, r)

! Returns true if the base location +/- maxdist in each dime
! intersects any part of the region defined by block r.

logical                          :: region_intersects_block
type(location_type), intent(in)  :: loc
real(r8),            intent(in)  :: maxdist
type(octree_type),   intent(in)  :: r

type(location_type) :: llb, urt

region_intersects_block = .false.

llb = set_location(loc%x-maxdist, loc%y-maxdist, loc%z-maxdist)
urt = set_location(loc%x+maxdist, loc%y+maxdist, loc%z+maxdist)

if (r%urt%x < llb%x) return
if (r%llb%x > urt%x) return

if (r%urt%y < llb%y) return
if (r%llb%y > urt%y) return

if (r%urt%z < llb%z) return
if (r%llb%z > urt%z) return

region_intersects_block = .true.

end function region_intersects_block

!----------------------------------------------------------------------------

recursive subroutine collect_nearby(r, base_loc, maxdist, loc_list, base_type, kind_list, &
                                   num_close, close_ind, dist)

type(octree_type), intent(in) :: r
type(location_type), intent(in)    :: base_loc, loc_list(:)
integer,             intent(in)    :: base_type, kind_list(:)
real(r8),            intent(in)    :: maxdist
integer,             intent(inout) :: num_close, close_ind(:)
real(r8), optional,  intent(inout) :: dist(:)

integer  :: i, j, k, t_ind
real(r8) :: this_dist

! revised plan:
! find min/max extents in each dim, +/- maxdist from base
! and only descend lower in a branch if the region intersects
! the range that this branch covers.  maxdist doesn't have to
! be a constant in this case - it can be specified per search.

! if this block does not intersect the region around the point,
! head back up.
if (.not. region_intersects_block(base_loc, maxdist, r)) return
!print *, 'location in/near region'

! FIXME: cannot go all the way down to child here in one routine.
! have to make a recursive routine to test each child one by one 
! and descend if possible intersection.

! if child node isn't a terminal node, test and either descend
! or skip.  if terminal node, search for nearby points and add
! them to the list if 'close enough'.

if (r%count < 0) then

   do i=1,2
    do j=1,2
     do k=1,2      
        call collect_nearby(r%children(i,j,k)%p, base_loc, maxdist, loc_list, base_type, kind_list, &
                                   num_close, close_ind, dist)
     enddo
    enddo
   enddo

else

   do i=1, r%count
      t_ind = r%index(i)
!print *, i, 'next index = ', t_ind

      ! Only compute distance if dist is present
      if(present(dist)) then
         this_dist = get_dist(base_loc, loc_list(t_ind), base_type, kind_list(t_ind))
!print *, 'this_dist = ', this_dist
         ! If dist is present and this loc's distance is less than cutoff, add it in list
         if(this_dist <= maxdist) then
            num_close = num_close + 1
            close_ind(num_close) = t_ind
            dist(num_close) = this_dist
!print *, 'adding loc to list'
         else
!print *, 'not adding loc to list, dist > maxdist', this_dist, maxdist
         endif
      else
         ! Dist isn't present; add this ob to list without computing distance
         num_close = num_close + 1
         close_ind(num_close) = t_ind
!print *, 'adding loc to list, dist not present'
      endif

   enddo

endif

end subroutine collect_nearby

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
   call error_handler(E_ERR, 'get_close_obs', errstring, source, revision, revdate, &
                      text2='optional arg "dist" is present; we are computing exact distances', &
                      text3='the exhaustive search should find an identical number of locations')
else if (num_close < cnum_close) then
   write(errstring, *) 'get_close (', num_close, ') should not be smaller than exhaustive search (', cnum_close, ')'
   call error_handler(E_ERR, 'get_close_obs', errstring, source, revision, revdate, &
                      text2='optional arg "dist" not present; we are returning a superset of close locations', &
                      text3='the exhaustive search should find an equal or lesser number of locations')
endif

! if they do not compare, we have the exhaustive lists here and can print out
! exactly which items in the list differ.

deallocate(cclose_ind, cclose_dist)

end subroutine exhaustive_report

!---------------------------------------------------------------------------

function vert_is_undef(loc)
 
! Given a location, return true if vertical coordinate is undefined, else false

logical                          :: vert_is_undef
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

vert_is_undef = .false.

end function vert_is_undef

!---------------------------------------------------------------------------

function vert_is_surface(loc)
 
! Given a location, return true if vertical coordinate is surface, else false

logical                          :: vert_is_surface
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

vert_is_surface = .false.

end function vert_is_surface

!---------------------------------------------------------------------------

function vert_is_pressure(loc)
 
! Given a location, return true if vertical coordinate is pressure, else false

logical                          :: vert_is_pressure
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

vert_is_pressure = .false.

end function vert_is_pressure

!---------------------------------------------------------------------------

function vert_is_height(loc)
 
! Given a location, return true if vertical coordinate is height, else false

logical                          :: vert_is_height
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

vert_is_height = .false.

end function vert_is_height

!---------------------------------------------------------------------------

function vert_is_level(loc)
 
! Given a location, return true if vertical coordinate is level, else false

logical                          :: vert_is_level
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

vert_is_level = .false.

end function vert_is_level

!---------------------------------------------------------------------------

function vert_is_scale_height(loc)

! Given a location, return true if vertical coordinate is scale height, else false

logical                          :: vert_is_scale_height
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

vert_is_scale_height = .false.

end function vert_is_scale_height

!---------------------------------------------------------------------------

function has_vertical_localization()
 
! Return the (opposite) namelist setting for horiz_dist_only.

logical :: has_vertical_localization

if ( .not. module_initialized ) call initialize_module

has_vertical_localization = .true.

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


!----------------------------------------------------------------------------
! end of location/channel/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
