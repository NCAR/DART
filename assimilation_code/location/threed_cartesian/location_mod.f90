! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module location_mod

! Implements location interfaces for full 3d cartesian space.
! Has interfaces to convert spherical lat/lon coords in degrees,
! plus a radius, into cartesian coords.

use      types_mod, only : r8, i8, MISSING_R8, MISSING_I, PI, RAD2DEG, DEG2RAD, OBSTYPELENGTH
use  utilities_mod, only : error_handler, E_ERR, ascii_file_format, &
                           E_MSG, open_file, close_file, set_output,       &
                           logfileunit, nmlfileunit, find_namelist_in_file,          &
                           check_namelist_read, do_output, do_nml_file,              &
                           do_nml_term, is_longitude_between
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform
use mpi_utilities_mod, only : my_task_id, task_count
use         obs_kind_mod, only : get_num_types_of_obs, get_name_for_type_of_obs, get_index_for_type_of_obs
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
          print_get_close_type, find_nearest, set_periodic

character(len=*), parameter :: source = 'threed_cartesian/location_mod.f90'

integer,              parameter :: LocationDims = 3
character(len = 129), parameter :: LocationName = "loc3Dcartesian"
character(len = 129), parameter :: LocationLName = &
                                   "threed cartesian locations: x, y, z"
character(len = 129), parameter :: LocationStorageOrder = "X Y Z"
character(len = 129), parameter :: LocationUnits = "none none none"

! The numeric value for the vertical component
! Used with CM1 when using normalized vertical localization
integer, parameter :: VERTISHEIGHT      =  1  ! by height (in meters)

type location_type
   private
   real(r8) :: x, y, z
end type location_type

!> @todo: is this how to approach it?
!> there are 8 different offsets depending on which side
!> of each axis midline the original point is on.
!> (the values are the same; the sign changes.  i can
!> add if() tests and do add/subtracts, or i can set the
!> offsets differently and do a single if() on the midline sign
!> or box number tests?  require number of boxes to be even,
!> and do integer tests on middle box values.
type periodic_info_type
   private
   real(r8) :: midline(3)
   integer  :: midbox(3)
   real(r8) :: offset(3)
end type periodic_info_type

! This version supports only regularly spaced boxes (the non-working octree code
! was removed.) 

type box_type
   private
   integer           :: num                      
   real(r8)          :: maxdist              ! furthest seperation between "close" locations  
   integer, allocatable  :: loc_box(:)           ! (nloc); List of loc indices in boxes
   integer, allocatable  :: count(:, :, :)       ! (nx, ny, nz); # of locs in each box
   integer, allocatable  :: start(:, :, :)       ! (nx, ny, nz); Start of list of locs in this box
   real(r8)          :: bot_x, top_x         ! extents in x, y, z
   real(r8)          :: bot_y, top_y
   real(r8)          :: bot_z, top_z
   real(r8)          :: x_width, y_width, z_width    ! widths of boxes in x,y,z
   real(r8)          :: nboxes_x, nboxes_y, nboxes_z ! based on maxdist how far to search
end type box_type




! Support more than a single cutoff distance.  nt is the count of
! distinct cutoffs, which are selected per specific observation type.
! The map associates the incoming location type with the
! corresponding gtt index.  There are only as many close_types
! as there are distinct cutoff distances; if more than one specific
! type has the same cutoff distances they share the type.
type get_close_type
   private
   integer                     :: nt                      ! The number of distinct cutoffs
   integer, allocatable        :: type_to_cutoff_map(:)   ! mapping of types to index
   type(box_type),allocatable  :: box(:)                  ! Array of box types 
end type get_close_type


! Some calls include information about the type or kind of the location.
! If the location refers to an observation it is possible to have both
! a specific type and a generic kind.  If the location refers to a
! state vector item, it only has a generic kind.  Some routines have
! a specific type and a generic kind as arguments; look carefully at
! the argument names before using. - Message stollen from 3d_Sphere
 
type(random_seq_type) :: ran_seq
logical               :: ran_seq_init = .false.
logical, save         :: module_initialized = .false.

character(len = 512) :: errstring
character(len = 512) :: msgstring, msgstring1, msgstring2

! Horizontal localization/cutoff values are passed in by the caller.
! The Vertical normalization values are globals; are specified by namelist
! here, and apply to all specific types unless a 'special list' is also specified
! that overrides the default values.

logical :: has_special_vertical_norms = .false.
integer :: num_special_vert_norms = 0

! Global storage for vertical distance normalization factors
integer, parameter    :: VERT_TYPE_COUNT = 1 ! only height supported
real(r8)              :: vert_normalization(VERT_TYPE_COUNT)
real(r8), allocatable :: per_type_vert_norm(:,:)  ! if doing per-type


! for sanity when i'm using arrays of length(3):
integer, parameter :: IX = 1
integer, parameter :: IY = 2
integer, parameter :: IZ = 3

! for quicker tests to see whether we do the easy version or descend
! into periodic hell, and to which level of hell we go into.
logical :: any_periodic = .false.
logical :: x_periodic = .false.
logical :: y_periodic = .false.
logical :: z_periodic = .false.
logical :: xy_periodic = .false.
logical :: xz_periodic = .false.
logical :: yz_periodic = .false.
logical :: xyz_periodic = .false.

! periodic option.  need axis, and min/max values

logical  :: x_is_periodic = .false.
real(r8) :: min_x_for_periodic = MISSING_R8
real(r8) :: max_x_for_periodic = MISSING_R8

logical  :: y_is_periodic = .false.
real(r8) :: min_y_for_periodic = MISSING_R8
real(r8) :: max_y_for_periodic = MISSING_R8

logical  :: z_is_periodic = .false.
real(r8) :: min_z_for_periodic = MISSING_R8
real(r8) :: max_z_for_periodic = MISSING_R8

type(periodic_info_type) :: loopy

!-----------------------------------------------------------------
! Namelist with default values

! Option for verification using exhaustive search, and debugging
! compare_to_correct = .true. is expensive
logical :: compare_to_correct = .false.

logical :: output_box_info  = .false.
integer :: print_box_level  = 0
integer :: debug  = 0

! for boxes
!integer :: nboxes           = 10000   ! currently unused
integer :: nx               = 10   ! box counts in each dimension
integer :: ny               = 10
integer :: nz               = 10


! Additional Items for Vertical Localization
!------------------------------------------------------------------------
! Namelist with default values
!
! vert_normalization_height            -> Number meters that give a distance equivalent
!                                         to one radian in horizontal
! special_vert_normalization_obs_types -> Which obs types to modify the default vert
!                                         normalization values
! special_vert_normalization_heights   -> value for each obs type
real(r8) :: vert_normalization_height       = 1.0_r8
integer, parameter :: MAX_ITEMS = 500
character(len=OBSTYPELENGTH) :: special_vert_normalization_obs_types(MAX_ITEMS)
real(r8) :: special_vert_normalization_heights(MAX_ITEMS)

namelist /location_nml/ &
   compare_to_correct, output_box_info, print_box_level,    &
   x_is_periodic, min_x_for_periodic, max_x_for_periodic,   &
   y_is_periodic, min_y_for_periodic, max_y_for_periodic,   &
   z_is_periodic, min_z_for_periodic, max_z_for_periodic,   &
   nx, ny, nz, vert_normalization_height,                   &
   special_vert_normalization_obs_types,                    &
   special_vert_normalization_heights, debug


!> @todo:
!>  for any point and combination of periodic axes, there are 2^N
!>  phantom points where we must compute distance.  if 'o' is
!>  the original location of the point, here are the possibilities:
!>    nothing periodic, 1 point:
!>         o
!>    single axis periodic, 2 points, e.g. x:
!>         o
!>         o + x_offset
!>    double axis periodic, 4 points, e.g. x,y:
!>         o
!>         o + x_offset
!>         o + y_offset
!>         o + x_offset + y_offset
!>    triple axis periodic, 8 points:
!>         o
!>         o + x_offset
!>         o + y_offset
!>         o + z_offset
!>         o + x_offset + y_offset
!>         o + x_offset + z_offset
!>         o + y_offset + z_offset
!>         o + x_offset + y_offset + z_offset
!>
!>  offsets:
!>  points below x_mid have a positive offset, points equal to
!>  and above x_mid have a negative offset.  doesn't matter if
!>  we add or subtract for points exactly == the midline as long
!>  as we do one or the other but not both.
!>
!-----------------------------------------------------------------

interface operator(==); module procedure loc_eq; end interface
interface operator(/=); module procedure loc_ne; end interface

interface set_location
   module procedure set_location_single
   module procedure set_location_array
   module procedure set_location_lonlat
end interface set_location

contains

!----------------------------------------------------------------------------

subroutine initialize_module

! things which need doing exactly once.

integer :: iunit, io, i, k, typecount, type_index

if (module_initialized) return

module_initialized = .true.

! give these initial values before reading them from the namelist.
special_vert_normalization_obs_types(:) = 'null'
special_vert_normalization_heights(:)       = missing_r8

! Read the namelist entry
call find_namelist_in_file("input.nml", "location_nml", iunit)
read(iunit, nml = location_nml, iostat = io)
call check_namelist_read(iunit, io, "location_nml")

! Write the namelist values to the log file

if(do_nml_file()) write(nmlfileunit, nml=location_nml)
if(do_nml_term()) write(     *     , nml=location_nml)

! Copy the normalization factors in the vertical into an array
vert_normalization(VERTISHEIGHT)      = vert_normalization_height

! If the user has set a namelist to make different vertical normalization factors
! when computing the localization distances, allocate and set that up here.
! It overrides the defaults above.

if (special_vert_normalization_obs_types(1) /= 'null' .or. &
    special_vert_normalization_heights(1)   /= missing_r8) then

   ! FIXME: add code to check for mismatched length lists.

   typecount = get_num_types_of_obs()  ! ignore function name, this is specific type count
   allocate(per_type_vert_norm(VERT_TYPE_COUNT, typecount))

   ! Set the defaults for all specific types not listed in the special list
   per_type_vert_norm(VERTISHEIGHT, :)      = vert_normalization_height

   ! Go through special-treatment observation kinds, if any.
   num_special_vert_norms = 0
   k = 0
   do i = 1, MAX_ITEMS
      if(special_vert_normalization_obs_types(i) == 'null') exit
      k = k + 1
   enddo
   num_special_vert_norms = k

   if (num_special_vert_norms > 0) has_special_vertical_norms = .true.

   do i = 1, num_special_vert_norms
      type_index = get_index_for_type_of_obs(special_vert_normalization_obs_types(i))
      if (type_index < 0) then
         write(msgstring, *) 'unrecognized TYPE_ in the special vertical normalization namelist:'
         call error_handler(E_ERR,'location_mod:', msgstring, source, &
                            text2=trim(special_vert_normalization_obs_types(i)))
      endif

      per_type_vert_norm(VERTISHEIGHT,      type_index) = special_vert_normalization_heights(i)

   enddo

   if (any(per_type_vert_norm == missing_r8)) then
      write(msgstring, *) 'one or more special vertical normalization values is uninitialized.'
      call error_handler(E_ERR,'location_mod:', &
                        'special vert normalization value namelist requires all 4 values per type', &
                        source, text2=trim(msgstring))
   endif

endif
!if (horiz_dist_only) then
!      call error_handler(E_MSG,'location_mod:', &
!      'Ignoring vertical separation when computing distances; horizontal distances only', &
!      source)
!else
call error_handler(E_MSG,'location_mod:', &
   'Including vertical separation when computing distances:', source)
write(msgstring,'(A,f17.5)') '        # meters ~ 1 horiz radian: ', vert_normalization_height
call error_handler(E_MSG,'location_mod:',msgstring,source)

if (allocated(per_type_vert_norm)) then
   typecount = get_num_types_of_obs()  ! ignore function name, this is specific type count
   do i = 1, typecount
      if    (per_type_vert_norm(VERTISHEIGHT,      i) /= vert_normalization_height) then
         write(msgstring,'(2A)') 'Altering default vertical normalization for type ', trim(get_name_for_type_of_obs(i))
         call error_handler(E_MSG,'location_mod:',msgstring,source)
         if (per_type_vert_norm(VERTISHEIGHT,      i) /= vert_normalization_height) then
            write(msgstring,'(A,f17.5)') '        # meters ~ 1 horiz meter: ', &
                  per_type_vert_norm(VERTISHEIGHT, i)
            call error_handler(E_MSG,'location_mod:',msgstring,source)
         endif
      endif
   enddo
endif


end subroutine initialize_module

!----------------------------------------------------------------------------

function get_dist(loc1, loc2, type1, kind2)

! @todo JDL I don't think I need to include the no_vert variable, but keep this
! in mind in case problems arise later on.  Probbably also don't need kind2 since 
! we are just looking at heights

! returns the distance between 2 locations

! the names are correct here - the first location gets the corresponding
! specific type; the second location gets a generic kind.
! these are included in case user-code wants to do a more sophisticated
! distance computation based on the kinds or types. The DART lib code
! doesn't use either of these.
!

type(location_type), intent(in) :: loc1, loc2
integer, optional,   intent(in) :: type1, kind2
real(r8)                        :: get_dist

real(r8) :: diff(3), next_diff(3)
logical  :: below_L1(3), below_L2(3)
real(r8) :: square_dist, this_dist


if ( .not. module_initialized ) call initialize_module

if (debug > 0) write(0,*)  'get_dist() called for 2 locations'

diff(IX) = loc1%x - loc2%x
diff(IY) = loc1%y - loc2%y
diff(IZ) = loc1%z - loc2%z

! the simple non-periodic case
if (.not. any_periodic) then
   if (debug > 0) write(0,*)  'non-periodic distance: '
   get_dist = dist_3d(diff,type1=type1)
   return
endif

! everything below deals with the periodic options
if (debug > 0) write(0,*)  'periodic distance: '
square_dist = dist_3d_sq(diff,type1=type1)
if (debug > 0) write(0,*)  'non-periodic 0 sq_dist, dist: ', square_dist, sqrt(square_dist)

! ok, not simple.
!> @todo: can we separate this?

!> while searching for the closest distance skip doing the
!> square root.  do it only at the end before returning.
!>
!> for all combinations of location 1 and location 2 where
!> they aren't both above or both below the midpoint add the
!> offset, compute the distance, and keep the minimium
below_L1 = find_my_quadrant(loc1)
if (debug > 0) write(0,*) 'below_L1, find my quadrant: ', below_L1
below_L2 = find_my_quadrant(loc2)
if (debug > 0) write(0,*) 'below_L2, find my quadrant: ', below_L2

if (xyz_periodic) then

   if (debug > 0) write(0,*)  'in the xyz_periodic case'

   if ((below_L1(IX) .neqv. below_L2(IX)) .or. &
       (below_L1(IY) .neqv. below_L2(IY)) .or. &
       (below_L1(IZ) .neqv. below_L2(IZ))) then
      next_diff = diff
      if (below_L1(IX) .neqv. below_L2(IX)) then ! Adjust x
         if (loc1%x > loopy%midline(IX)) then
            next_diff(IX) = loopy%offset(IX) - next_diff(IX)
         else
            next_diff(IX) = loopy%offset(IX) + next_diff(IX)
         endif
      endif
      if (below_L1(IY) .neqv. below_L2(IY)) then ! Adjust y
         if (loc1%y > loopy%midline(IY)) then
            next_diff(IY) = loopy%offset(IY) - next_diff(IY)
         else
            next_diff(IY) = loopy%offset(IY) + next_diff(IY)
         endif
      endif
      if (below_L1(IZ) .neqv. below_L2(IZ)) then ! Adjust z
         if (loc1%z > loopy%midline(IZ)) then
            next_diff(IZ) = loopy%offset(IZ) - next_diff(IZ)
         else
            next_diff(IZ) = loopy%offset(IZ) + next_diff(IZ)
         endif
      endif
      this_dist = dist_3d_sq(next_diff,type1=type1)
  if (debug > 0) write(0,*)  'periodic XYZ dist: ', square_dist, this_dist
      if(this_dist < square_dist) then
         square_dist = this_dist
      endif
   else
  if (debug > 0) write(0,*)  'XYZ not on diff sides of midline'
   endif

endif

if (xy_periodic) then

   if (debug > 0) write(0,*)  'in the xy_periodic case'

   if ((below_L1(IX) .neqv. below_L2(IX)) .or. &
       (below_L1(IY) .neqv. below_L2(IY))) then
      next_diff = diff
      if (below_L1(IX) .neqv. below_L2(IX)) then
         if (loc1%x > loopy%midline(IX)) then
            next_diff(IX) = loopy%offset(IX) - next_diff(IX)
         else
            next_diff(IX) = loopy%offset(IX) + next_diff(IX)
         endif
      endif
      if (below_L1(IY) .neqv. below_L2(IY)) then
         if (loc1%y > loopy%midline(IY)) then
            next_diff(IY) = loopy%offset(IY) - next_diff(IY)
         else
            next_diff(IY) = loopy%offset(IY) + next_diff(IY)
         endif
      endif

      this_dist = dist_3d_sq(next_diff,type1=type1)
  if (debug > 0) write(0,*)  'periodic XY dist: ', square_dist, this_dist
      if(this_dist < square_dist) then
         square_dist = this_dist
      endif
   else
  if (debug > 0) write(0,*)  'XY not on diff sides of midline'
   endif

endif


if (xz_periodic) then

   if (debug > 0) write(0,*)  'in the xz_periodic case'

   if ((below_L1(IX) .neqv. below_L2(IX)) .or. &
       (below_L1(IZ) .neqv. below_L2(IZ))) then
      next_diff = diff
      if (below_L1(IX) .neqv. below_L2(IX)) then ! Adjust x
         if (loc1%x > loopy%midline(IX)) then
            next_diff(IX) = loopy%offset(IX) - next_diff(IX)
         else
            next_diff(IX) = loopy%offset(IX) + next_diff(IX)
         endif
      endif
      if (below_L1(IZ) .neqv. below_L2(IZ)) then ! Adjust z
         if (loc1%z > loopy%midline(IZ)) then
            next_diff(IZ) = loopy%offset(IZ) - next_diff(IZ)
         else
            next_diff(IZ) = loopy%offset(IZ) + next_diff(IZ)
         endif
      endif
      this_dist = dist_3d_sq(next_diff,type1=type1)
  if (debug > 0) write(0,*)  'periodic XZ dist: ', square_dist, this_dist
      if(this_dist < square_dist) then
         square_dist = this_dist
      endif
   else
  if (debug > 0) write(0,*)  'XZ not on diff sides of midline'
   endif

endif

if (yz_periodic) then

   if (debug > 0) write(0,*)  'in the yz_periodic case'

   if ((below_L1(IY) .neqv. below_L2(IY)) .or. &
       (below_L1(IZ) .neqv. below_L2(IZ))) then
      next_diff = diff
      if (below_L1(IY) .neqv. below_L2(IY)) then ! Adjust y
         if (loc1%y > loopy%midline(IY)) then
            next_diff(IY) = loopy%offset(IY) - next_diff(IY)
         else
            next_diff(IY) = loopy%offset(IY) + next_diff(IY)
         endif
      endif
      if (below_L1(IZ) .neqv. below_L2(IZ)) then ! Adjust z
         if (loc1%z > loopy%midline(IZ)) then
            next_diff(IZ) = loopy%offset(IZ) - next_diff(IZ)
         else
            next_diff(IZ) = loopy%offset(IZ) + next_diff(IZ)
         endif
      endif
      this_dist = dist_3d_sq(next_diff,type1=type1)
  if (debug > 0) write(0,*)  'periodic YZ dist: ', square_dist, this_dist
      if(this_dist < square_dist) then
         square_dist = this_dist
      endif
   else
  if (debug > 0) write(0,*)  'YZ not on diff sides of midline'
   endif

endif

if (x_periodic) then

   if (debug > 0) write(0,*)  'in the x_periodic case'

   if (below_L1(IX) .neqv. below_L2(IX)) then
      next_diff = diff
      if (loc1%x > loopy%midline(IX)) then
         next_diff(IX) = loopy%offset(IX) - next_diff(IX)
      else
         next_diff(IX) = loopy%offset(IX) + next_diff(IX)
      endif
      this_dist = dist_3d_sq(next_diff,type1=type1)

  if (debug > 0) write(0,*)  'periodic X dist: ', square_dist, this_dist
      if(this_dist < square_dist) then
         square_dist = this_dist
      endif
   else
  if (debug > 0) write(0,*)  'X not on diff sides of midline'
   endif

endif

if (y_periodic) then
   if (debug > 0) write(0,*)  'in the y_periodic case'

   if (below_L1(IY) .neqv. below_L2(IY)) then
      next_diff = diff
      if (loc1%y >  loopy%midline(IY)) then
         next_diff(IY) = loopy%offset(IY) - next_diff(IY) 
      else
         next_diff(IY) = loopy%offset(IY) + next_diff(IY)    
      endif

      this_dist = dist_3d_sq(next_diff,type1=type1)

  if (debug > 0) write(0,*)  'periodic Y dist: ', square_dist, this_dist
      if(this_dist < square_dist) then
         square_dist = this_dist
      endif
   else
  if (debug > 0) write(0,*)  'Y not on diff sides of midline'
   endif

endif

if (z_periodic) then

   if (debug > 0) write(0,*)  'in the z_periodic case'

   if (below_L1(IZ) .neqv. below_L2(IZ)) then
      next_diff = diff
      if (loc1%z >  loopy%midline(IZ)) then
         next_diff(IZ) = loopy%offset(IZ) - next_diff(IZ)
      else
         next_diff(IZ) = loopy%offset(IZ) + next_diff(IZ)
      endif
      this_dist = dist_3d_sq(next_diff,type1=type1)
  if (debug > 0) write(0,*)  'periodic Y dist: ', square_dist, this_dist
      if(this_dist < square_dist) then
         square_dist = this_dist
      endif
   else
  if (debug > 0) write(0,*)  'Z not on diff sides of midline'
   endif

endif

  if (debug > 0) write(0,*)  'periodic dist before sqrt: ', square_dist
if (square_dist < 0.0_r8) then
  if (debug > 0) write(0,*)  'trouble!  square_dist is negative: ', square_dist
  stop
endif

get_dist = sqrt(square_dist)
if (debug > 0) write(0,*) 'get_dist is returning: ', get_dist

end function get_dist

!---------------------------------------------------------------------------

! return the 3d distance given the separation along each axis

function dist_3d(separation,type1) result(val)

real(r8), intent(in) :: separation(3)
real(r8) :: val, vert_normal
integer, optional,   intent(in) :: type1 ! JDL Addition

if (allocated(per_type_vert_norm)) then
   if (.not. present(type1)) then
      write(msgstring, *) 'obs type required in get_dist`() if doing per-type vertical normalization 3d'
      call error_handler(E_MSG, 'get_dist', msgstring, source)
   endif 
   vert_normal = separation(IZ)/per_type_vert_norm(VERTISHEIGHT, type1)
else
   vert_normal = separation(IZ)/vert_normalization(VERTISHEIGHT)

endif

val = sqrt(separation(IX)*separation(IX) + &
           separation(IY)*separation(IY) + &
           vert_normal*vert_normal )

if (debug > 0) write(0,*)  'dist_3d called, distance computed: ', val
if (debug > 0) write(0,*)  'XYZ separations: ', separation
end function dist_3d

!---------------------------------------------------------------------------

! return the square of the 3d distance given the separation along each axis
! (saves doing a square root)

function dist_3d_sq(separation,type1) result(val)

real(r8), intent(in) :: separation(3)
real(r8) :: val, vert_normal
integer, optional,   intent(in) :: type1 ! JDL Addition

if (allocated(per_type_vert_norm)) then
   if (.not. present(type1)) then
      write(msgstring, *) 'obs type required in get_dist`() if doing per-type vertical normalization 3d sq'
      call error_handler(E_MSG, 'get_dist', msgstring, source)
   endif 
   vert_normal = separation(IZ) / per_type_vert_norm(VERTISHEIGHT, type1)
else
   vert_normal = separation(IZ) / vert_normalization(VERTISHEIGHT)
endif

val = separation(IX)*separation(IX) + &
      separation(IY)*separation(IY) + &
      vert_normal * vert_normal

if (debug > 0) write(0,*)  'dist_3d_sq called, distance computed: ', val
if (debug > 0) write(0,*)  'XYZ separations: ', separation
end function dist_3d_sq

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

get_location(IX) = loc%x
get_location(IY) = loc%y
get_location(IZ) = loc%z

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

set_location_array = set_location_single(list(IX), list(IY), list(IZ))

end function set_location_array

!----------------------------------------------------------------------------

function set_location_lonlat(lon, lat, height, radius)

! location semi-independent interface routine
! given a lon, lat, height and radius, compute X,Y,Z and set location

real(r8), intent(in) :: lon, lat, height, radius
type (location_type) :: set_location_lonlat

real(r8) :: x, y, z
real(r8) :: rlat, rlon

if ( .not. module_initialized ) call initialize_module


rlat = lat * deg2rad
rlon = lon * deg2rad

x = radius * cos(rlon) * cos(rlat)
y = radius * sin(rlon) * cos(rlat)
z = radius * sin(rlat)

set_location_lonlat%x = x
set_location_lonlat%y = y
set_location_lonlat%z = z

end function set_location_lonlat

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
      write(locfile, '(''loc3Dxyz'')' )
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

character(len=8) :: header

if ( .not. module_initialized ) call initialize_module

if (ascii_file_format(fform)) then
   read(locfile, '(A8)' ) header
   if(header /= 'loc3Dxyz') then
         write(errstring,*)'Expected location header "loc3Dxyz" in input file, got ', header
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

l(IX) = 'X'
l(IY) = 'Y'
l(IZ) = 'Z'

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
      write(*,*) 'min/max was: ', minv, maxv

   endif

enddo

location%x = v(IX)
location%y = v(IY)
location%z = v(IZ)

end subroutine interactive_location

!----------------------------------------------------------------------------
! Initializes get_close accelerator 

subroutine get_close_init(gc, num, maxdist, locs, maxdist_list)
type(get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
real(r8),             intent(in)    :: maxdist
type(location_type),  intent(in)    :: locs(:)
real(r8), intent(in), optional      :: maxdist_list(:)

integer :: i, j, k, cum_start, l, n
integer :: x_box(num), y_box(num), z_box(num)
integer :: tstart(nx, ny, nz)


integer :: typecount, distcount
real(r8), allocatable :: distlist(:)


if ( .not. module_initialized ) call initialize_module

typecount = get_num_types_of_obs()
allocate(gc%type_to_cutoff_map(typecount)) 

if (present(maxdist_list)) then
   if (size(maxdist_list) .ne. typecount) then
      write(msgstring,'(A,I8,A,I8)')'maxdist_list len must equal number of specific types, ', &
                                    size(maxdist_list), ' /= ', typecount
      call error_handler(E_ERR, 'get_close_init', msgstring, source)
   endif

   allocate(distlist(typecount))
   call distinct_values(maxdist_list, distcount, distlist, gc%type_to_cutoff_map)
   gc%nt = distcount
   if (gc%nt <= 0) then
      write(msgstring,'(A)')'error getting count of distinct cutoff dists; should not happen'
      call error_handler(E_ERR, 'get_close_init', msgstring, source)
   endif
else
   gc%nt = 1
   gc%type_to_cutoff_map(:) = 1
endif

allocate(gc%box(gc%nt))

if (present(maxdist_list)) then
   do i=1, gc%nt
      gc%box(i)%maxdist = distlist(i)
   enddo
else
   ! no per-type settings, everyone uses same distance
   gc%box(1)%maxdist = maxdist
endif

if (present(maxdist_list)) deallocate(distlist)

! Allocate the storage for the grid dependent boxes
do i=1, gc%nt
   allocate(gc%box(i)%count(nx, ny, nz), gc%box(i)%start(nx, ny, nz))
   gc%box(i)%count      = -1
   gc%box(i)%start      = -1
enddo

! store the location counts in all derived types
do i=1, gc%nt
   gc%box(i)%num = num
enddo


! If there are no locs to operate on, no point in going any further.
if (num == 0) return

! Set the maximum localization distance
!gc%maxdist = maxdist

! Allocate storage for number dependent part
!allocate(gc%box%loc_box(num))
!gc%box%loc_box(:) = -1

! Allocate the storage for the grid dependent boxes
!allocate(gc%box%count(nx,ny,nz), gc%box%start(nx,ny,nz))
!gc%box%count  = -1
!gc%box%start  = -1

do i=1, gc%nt
   ! Allocate storage for locs number dependent part
   allocate(gc%box(i)%loc_box(num))
   gc%box(i)%loc_box(:) = -1
enddo

! Set the value of num_locs in the structure
!gc%num = num

! If num == 0, no point in going any further.
!if (num == 0) return

!! FIXME: compute nx, ny, nz from nboxes?  or put in namelist
!nx = nint(real(nboxes, r8)**0.33333)   ! roughly cube root
!ny = nint(real(nboxes, r8)**0.33333)   ! roughly cube root
!nz = nint(real(nboxes, r8) / real(nx * ny, r8))  ! whatever is left

!> TODO:  plan --
!> - if periodic in a dimension, translate all lower half points
!> into a halo box beyond the upper half, and translate all upper half
!> points into a halo box below the lower half.
!> - subtract or add the original min or max to translate the locations
!> - nboxes along that dim must be even.
!> - use the original location indices in all box lists?
!> - must have a way to have or compute both the original x/y/z
!> values for a location, but also the translated location(s) when
!> computing distance.
!> - mark the periodic boxes somehow - with offsets?  original boxes
!> have 0 offsets, halo boxes have them.
!> - need a box offset (which is nboxes/2 for that axis) when
!> computing box number
!> - negative box numbers for < 1?  allocate(-10:30) for 20 boxes?
!> or have an offset and 0 is the halo, nboxes/2 is the real start
!> and nboxes*2 is the end of the halo on the upper side.

! Determine where the boxes should be for this set of locs and maxdist
do n=1, gc%nt 
    ! Determine where the boxes should be for this set of locs and maxdist
    call find_box_ranges(gc%box(n), locs, num) 

    ! Begin by computing the number of locations in each box in x,y,z
    gc%box(n)%count = 0
    do i = 1, num

    !write(0,*)  i, locs(i)%x, locs(i)%y, locs(i)%z
       x_box(i) = floor((locs(i)%x - gc%box(n)%bot_x) / gc%box(n)%x_width) + 1
       if(x_box(i) > nx) x_box(i) = nx
       if(x_box(i) < 1)  x_box(i) = 1

       y_box(i) = floor((locs(i)%y - gc%box(n)%bot_y) / gc%box(n)%y_width) + 1
       if(y_box(i) > ny) y_box(i) = ny
       if(y_box(i) < 1)  y_box(i) = 1

       z_box(i) = floor((locs(i)%z - gc%box(n)%bot_z) / gc%box(n)%z_width) + 1
       if(z_box(i) > nz) z_box(i) = nz
       if(z_box(i) < 1)  z_box(i) = 1

       gc%box(n)%count(x_box(i), y_box(i), z_box(i)) = gc%box(n)%count(x_box(i), y_box(i), z_box(i)) + 1
    !write(0,*)  'adding count to box ', x_box(i), y_box(i), z_box(i), &
    !                                 gc%box%count(x_box(i), y_box(i), z_box(i))
    end do

    ! Figure out where storage for each boxes members should begin
    cum_start = 1
    do i = 1, nx
      do j = 1, ny
         do k = 1, nz
            gc%box(n)%start(i, j, k) = cum_start
            cum_start = cum_start + gc%box(n)%count(i, j, k)
         end do
      end do
   end do

   ! Now we know how many are in each box, get a list of which are in each box
   tstart = gc%box(n)%start
   do i = 1, num
      gc%box(n)%loc_box(tstart(x_box(i), y_box(i), z_box(i))) = i
      tstart(x_box(i), y_box(i), z_box(i)) = tstart(x_box(i), y_box(i), z_box(i)) + 1
   end do
   do i = 1, nx
      do j = 1, ny
         do k = 1, nz
   !if (gc%box%count(i,j,k) > 0) write(0,*)  i,j,k, gc%box%count(i,j,k), gc%box%start(i,j,k)
            do l=1, gc%box(n)%count(i,j,k)
   !write(0,*)  l, gc%box%loc_box(l)
            enddo
         end do
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
      call print_get_close_type(gc, 1, print_box_level)
   else if (print_box_level >= 2 .or. print_box_level < 0) then
      ! print status was false, but turn on temporarily
      ! to output box info from all tasks.
      call set_output(.true.)
      call print_get_close_type(gc, 1, print_box_level)
      call set_output(.false.)
   endif
endif

end subroutine get_close_init

!----------------------------------------------------------------------------

subroutine get_close_destroy(gc)
type(get_close_type), intent(inout) :: gc
integer :: i

do i = 1, gc%nt 
  if (allocated(gc%box(i)%loc_box)) deallocate(gc%box(i)%loc_box)
  deallocate(gc%box(i)%count, gc%box(i)%start)
enddo
deallocate(gc%type_to_cutoff_map)
deallocate(gc%box)
!deallocate(gc%box%loc_box, gc%box%count, gc%box%start) ! ORIG

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

integer :: x_box, y_box, z_box, i, j, k, l, bt
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


! handle identity obs correctly -- only support them if there are no
! per-obs-type cutoffs.  identity obs don't have a specific type, they
! only have a generic kind based on what index into the state vector is
! specified.  if this is absolutely needed by someone, i can rejigger the
! code to try to use the default cutoff for identity obs, but it's tricky
! to identify which obs type is using the default.  in theory it's possible
! to specify a default cutoff and then specify a per-type cutoff for every
! other type in the system.  in that case i don't have a map entry for the
! default, and it's a pain to construct one.  so i'm punting for now.
if (base_type < 0) then
   if (gc%nt > 1) then
      write(msgstring,  '(A)') 'no support for identity obs if per-obs-type cutoffs are specified'
      write(msgstring1, '(A)') 'contact dart support if you have a need for this combination'
      call error_handler(E_ERR, 'get_close', msgstring, source, text2=msgstring1)
   endif
   bt = 1
else
   ! map from type index to gtt index
   if (base_type < 1 .or. base_type > size(gc%type_to_cutoff_map)) then
      write(msgstring,'(A,I8)')'base_type out of range, is ', base_type
      write(msgstring1,'(A,2I8)')'must be between ', 1, size(gc%type_to_cutoff_map)
      call write_location (0, base_loc, charstring=msgstring2)
      call error_handler(E_ERR, 'get_close', msgstring, source, &
                         text2=msgstring1, text3=msgstring2)
   endif
   bt = gc%type_to_cutoff_map(base_type)
   if (bt < 1 .or. bt > gc%nt) then
      write(msgstring,'(A,I8)')'mapped type index out of range, is ', bt
      write(msgstring1,'(A,2I8)')'must be between ', 1, gc%nt
      write(msgstring2, '(A)')'internal error, should not happen.  Contact DART Support'
      call error_handler(E_ERR, 'get_close', msgstring, source, &
                         text2=msgstring1, text3=msgstring2)
   endif
endif



! the list of locations in the loc() argument must be the same
! as the list of locations passed into get_close_init(), so
! gc%num and size(loc) better be the same.   if the list changes,
! you have to destroy the old gc and init a new one.
if (size(locs) /= gc%box(bt)%num) then
   !write(errstring,*)'locs() array must match one passed to get_close_init()'
   !call error_handler(E_ERR, 'get_close_boxes', errstring, source)
   write(msgstring,'(A)')'locs() array must match one passed to get_close_init()'
   write(msgstring1,'(A,2I8)')'init size, current list size: ', gc%box(bt)%num, size(locs)
   write(msgstring2,'(A,I8)')'bt = ', bt
   call error_handler(E_ERR, 'get_close', msgstring, source, &
                      text2=msgstring1, text3=msgstring2)
endif

! If num == 0, no point in going any further.
if (gc%box(bt)%num == 0) return

this_maxdist = gc%box(bt)%maxdist

!> @todo this is doing an exhaustive search each time.  expensive
!> but should give the right answer.
if(.true.) then
   if (present(dist)) then
      call exhaustive_collect(gc%box(bt), base_loc, base_type, locs, &
                              num_close, close_ind, dist)
   else
      allocate(cdist(size(locs)))
      call exhaustive_collect(gc%box(bt), base_loc, base_type, locs, &
                              num_close, close_ind, cdist)
      deallocate(cdist)
   endif
   return
endif

! For validation, it is useful to be able to compare against exact
! exhaustive search
if(compare_to_correct) then
   allocate(cclose_ind(size(locs)), cdist(size(locs)))
   call exhaustive_collect(gc%box(bt), base_loc, base_type, locs, &
                           cnum_close, cclose_ind, cdist)
endif


! Begin by figuring out which box the base loc is in
x_box = floor((base_loc%x - gc%box(bt)%bot_x) / gc%box(bt)%x_width) + 1
y_box = floor((base_loc%y - gc%box(bt)%bot_y) / gc%box(bt)%y_width) + 1
z_box = floor((base_loc%z - gc%box(bt)%bot_z) / gc%box(bt)%z_width) + 1

! If it is not in any box, then it is more than the maxdist away from everybody
if(x_box > nx .or. x_box < 1 .or. x_box < 0) return
if(y_box > ny .or. y_box < 1 .or. y_box < 0) return
if(z_box > nz .or. z_box < 1 .or. z_box < 0) return

! figure out how many boxes need searching
! FIXME: if we support a variable maxdist, nboxes_X will need to
! be computed on the fly here instead of precomputed at init time.
!write(0,*)  'nboxes x, y, z = ', gc%box%nboxes_x, gc%box%nboxes_y, gc%box%nboxes_z
!write(0,*)  'base_loc in box ', x_box, y_box, z_box

start_x = x_box - gc%box(bt)%nboxes_x
if (start_x < 1) start_x = 1
end_x = x_box + gc%box(bt)%nboxes_x
if (end_x > nx) end_x = nx

start_y = y_box - gc%box(bt)%nboxes_y
if (start_y < 1) start_y = 1
end_y = y_box + gc%box(bt)%nboxes_y
if (end_y > ny) end_y = ny

start_z = z_box - gc%box(bt)%nboxes_z
if (start_z < 1) start_z = 1
end_z = z_box + gc%box(bt)%nboxes_z
if (end_z > nz) end_z = nz

!write(0,*)  'looping from '
!write(0,*)  'x: ', start_x, end_x
!write(0,*)  'y: ', start_y, end_y
!write(0,*)  'z: ', start_z, end_z

! Next, loop through each box that is close to this box
do i = start_x, end_x
   do j = start_y, end_y
      do k = start_z, end_z

         ! Box to search is i,j,k
         n_in_box = gc%box(bt)%count(i, j, k)
         st = gc%box(bt)%start(i,j,k)


         ! Loop to check how close all loc in the box are; add those that are close
         do l = 1, n_in_box

            t_ind = gc%box(bt)%loc_box(st - 1 + l)
!write(0,*)  'l, t_ind = ', l, t_ind

            ! Only compute distance if dist is present
            if(present(dist)) then
               this_dist = get_dist(base_loc, locs(t_ind), type1=base_type)
!write(0,*)  'this_dist = ', this_dist
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
   deallocate(cclose_ind, cdist)
endif


end subroutine get_close

!--------------------------------------------------------------------------

subroutine find_box_ranges(box, locs, num)

! Finds boundaries for x,y,z boxes.
! FIXME: ways boxes could be divided:
!  - evenly along each axis
!  - octree-like, divide each axis so roughly half the points are
!     on each side of the dividing plane.
!  - about 100 other schemes

type(box_type), intent(inout)       :: box
integer,              intent(in)    :: num
type(location_type),  intent(in)    :: locs(num)

!logical :: old_out

if (x_is_periodic) then
   box%bot_x = min_x_for_periodic
   box%top_x = max_x_for_periodic
else
   box%bot_x = minval(locs(:)%x - box%maxdist)
   box%top_x = maxval(locs(:)%x + box%maxdist)
endif

if (y_is_periodic) then
   box%bot_y = min_x_for_periodic
   box%top_y = max_x_for_periodic
else
   box%bot_y = minval(locs(:)%y - box%maxdist)
   box%top_y = maxval(locs(:)%y + box%maxdist)
endif

if (z_is_periodic) then
   box%bot_z = min_x_for_periodic
   box%top_z = max_x_for_periodic
else
   box%bot_z = minval(locs(:)%z - box%maxdist)
   box%top_z = maxval(locs(:)%z + box%maxdist)
endif

if (debug > 0) print *, 'nx/ny/nz: ', nx, ny, nz
if (debug > 0) print *, 'bots: ', box%bot_x, box%bot_y, box%bot_z
if (debug > 0) print *, 'tops: ', box%top_x, box%top_y, box%top_z

box%x_width = max(1.0_r8, (box%top_x - box%bot_x) / nx)
box%y_width = max(1.0_r8, (box%top_y - box%bot_y) / ny)
box%z_width = max(1.0_r8, (box%top_z - box%bot_z) / nz)

if (debug > 0) print *, 'widths = ', box%x_width, box%y_width, box%z_width

! FIXME:  compute a sphere of radius maxdist and see how
! many boxes in x, y, z that would include.
if (box%x_width <= 0.0_r8) &
   call error_handler(E_ERR, 'find_box_ranges', 'x_width <= 0', source)
if (box%y_width <= 0.0_r8) &
   call error_handler(E_ERR, 'find_boy_ranges', 'y_width <= 0', source)
if (box%z_width <= 0.0_r8) &
   call error_handler(E_ERR, 'find_boz_ranges', 'z_width <= 0', source)

box%nboxes_x = aint((box%maxdist + (box%x_width-1)) / box%x_width)
box%nboxes_y = aint((box%maxdist + (box%y_width-1)) / box%y_width)
box%nboxes_z = aint((box%maxdist + (box%z_width-1)) / box%z_width)

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

!subroutine find_nearest(box, base_loc, loc_list, nearest, rc)
 !type(box_type), intent(in), target  :: box
subroutine find_nearest(gc,base_loc,loc_list,nearest, rc)
 type(get_close_type), intent(in), target :: gc
 type(location_type),  intent(in)  :: base_loc
 type(location_type),  intent(in)  :: loc_list(:)
 integer,              intent(out) :: nearest
 integer,              intent(out) :: rc

! find the nearest point in the get close list to the specified point

integer :: x_box, y_box, z_box, i, j, k, l
integer :: start_x, end_x, start_y, end_y, start_z, end_z
integer :: n_in_box, st, t_ind
integer :: box_num
real(r8) :: this_dist, dist

! First, set the intent out arguments to a missing value
nearest = -99
rc = -1
dist = 1e38_r8                ! something big and positive.
box_num = 1                   ! JDL - HARD CODE TO JUST GET THE FIRST BOX OBJECT
! the list of locations in the loc() argument must be the same
! as the list of locations passed into get_close_init(), so
! gc%num and size(loc) better be the same.   if the list changes,
! you have to destroy the old gc and init a new one.
if (size(loc_list) /= gc%box(box_num)%num) then
   write(errstring,*)'loc() array must match one passed to get_close_init()'
   call error_handler(E_ERR, 'find_nearest_boxes', errstring, source)
endif

! If num == 0, no point in going any further.
if (gc%box(box_num)%num == 0) return

!--------------------------------------------------------------

! Begin by figuring out which box the base loc is in
x_box = floor((base_loc%x - gc%box(box_num)%bot_x) / gc%box(box_num)%x_width) + 1
y_box = floor((base_loc%y - gc%box(box_num)%bot_y) / gc%box(box_num)%y_width) + 1
z_box = floor((base_loc%z - gc%box(box_num)%bot_z) / gc%box(box_num)%z_width) + 1

! FIXME: this should figure out if it's > n or < 0 and
! set to n or 0 and always return something.
! If it is not in any box, then it is more than the maxdist away from everybody
if(x_box > nx .or. x_box < 1 .or. x_box < 0) return
if(y_box > ny .or. y_box < 1 .or. y_box < 0) return
if(z_box > nz .or. z_box < 1 .or. z_box < 0) return

! figure out how many boxes need searching
! FIXME: if we support a variable maxdist, nboxes_X will need to
! be computed on the fly here instead of precomputed at init time.
!write(0,*)  'nboxes x, y, z = ', gc%box%nboxes_x, gc%box%nboxes_y, gc%box%nboxes_z
!write(0,*)  'base_loc in box ', x_box, y_box, z_box

start_x = x_box - gc%box(box_num)%nboxes_x
if (start_x < 1) start_x = 1
end_x = x_box + gc%box(box_num)%nboxes_x
if (end_x > nx) end_x = nx

start_y = y_box - gc%box(box_num)%nboxes_y
if (start_y < 1) start_y = 1
end_y = y_box + gc%box(box_num)%nboxes_y
if (end_y > ny) end_y = ny

start_z = z_box - gc%box(box_num)%nboxes_z
if (start_z < 1) start_z = 1
end_z = z_box + gc%box(box_num)%nboxes_z
if (end_z > nz) end_z = nz

!write(0,*)  'looping from '
!write(0,*)  'x: ', start_x, end_x
!write(0,*)  'y: ', start_y, end_y
!write(0,*)  'z: ', start_z, end_z
! JDL WARNING STATEMENT IN CASE YOU HAPPEN TO CALL GET_DIST() FROM THIS SUBROUTINE
! YOU SHOULDN'T BE USING THIS SUBROUTINE
print*,'WARNING - IF YOUR HERE MAKE SURE TO WORK ON FIND_NEAREST SUBROUTINE'

! Next, loop through each box that is close to this box
do i = start_x, end_x
   do j = start_y, end_y
      do k = start_z, end_z

         ! Box to search is i,j,k
         n_in_box = gc%box(box_num)%count(i, j, k)
         st = gc%box(box_num)%start(i,j,k)


         ! Loop to check how close all loc in the box are; add those that are close
         do l = 1, n_in_box

            t_ind = gc%box(box_num)%loc_box(st - 1 + l)
!write(0,*)  'l, t_ind = ', l, t_ind

            this_dist = get_dist(base_loc, loc_list(t_ind)) !JDL Probably do not need to include base type here because this is never called during DA (OR IT SHOULDNT BE) 
!write(0,*)  'this_dist = ', this_dist
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

!----------------------------------------------------------------------------

! return an array of 3 logicals to say whether the given
! location is below the midpoint along each axis.  points
! exactly on the midline are not below.

function find_my_quadrant(loc) result(val)

type(location_type), intent(in) :: loc
logical, dimension(3) :: val

if (debug > 0) write(0,*) 'testing X against midline: ', loc%x, loopy%midline(IX)
if (debug > 0) write(0,*) 'testing Y against midline: ', loc%y, loopy%midline(IY)
if (debug > 0) write(0,*) 'testing Z against midline: ', loc%z, loopy%midline(IZ)

val(IX) = (loc%x < loopy%midline(IX))
val(IY) = (loc%y < loopy%midline(IY))
val(IZ) = (loc%z < loopy%midline(IZ))

if (debug > 0) write(0,*) 'find_my_quadrant returns: ', val

end function find_my_quadrant

!----------------------------------------------------------------------------

! set an axis to be periodic, and give limits.
! specify the axis with a single letter:  X, Y, or Z
! the min and max values will be used to compute the
! size and location of the domain along that axis.

subroutine set_periodic(axis, minlimit, maxlimit)

character(len=1), intent(in) :: axis
real(r8),         intent(in) :: minlimit
real(r8),         intent(in) :: maxlimit

select case (axis)
   case ('X', 'x')
      x_is_periodic = .true.
      min_x_for_periodic = minlimit
      max_x_for_periodic = maxlimit

   case ('Y', 'y')
      y_is_periodic = .true.
      min_y_for_periodic = minlimit
      max_y_for_periodic = maxlimit

   case ('Z', 'z')
      z_is_periodic = .true.
      min_z_for_periodic = minlimit
      max_z_for_periodic = maxlimit

   case default
      call error_handler(E_ERR, 'set_periodic', 'unrecognized axis name', source)
end select

call recompute_periodic()

end subroutine set_periodic

!----------------------------------------------------------------------------

! internal routine - reset the state of the axes.
! and precompute some stuff to make this simpler.

subroutine recompute_periodic()

! Shortcut to see if any axis is periodic
any_periodic = (x_is_periodic .or. y_is_periodic .or. z_is_periodic)

! For now, if any periodic force all nx,ny,nz to be even
if (any_periodic) then
   if ((((nx/2.0_r8)*2.0_r8) /= nx) .or. &
       (((ny/2.0_r8)*2.0_r8) /= ny) .or. &
       (((nz/2.0_r8)*2.0_r8) /= nz)) then
      call error_handler(E_ERR, 'locations XYZ', 'If any axes are periodic, all box counts must be even', &
                         source)
   endif

   if (x_is_periodic) x_periodic = .true.
   if (y_is_periodic) y_periodic = .true.
   if (z_is_periodic) z_periodic = .true.
   if (x_is_periodic .and. y_is_periodic) xy_periodic = .true.
   if (x_is_periodic .and. z_is_periodic) xz_periodic = .true.
   if (y_is_periodic .and. z_is_periodic) yz_periodic = .true.
   if (x_is_periodic .and. y_is_periodic .and. z_is_periodic) xyz_periodic = .true.
endif


! fill this in in either case.  if nothing periodic, fill with
! values that are likely to cause crashes if used unintentionally
if (any_periodic) then
   loopy%midline(IX) = (min_x_for_periodic + max_x_for_periodic ) / 2.0_r8
   loopy%midline(IY) = (min_y_for_periodic + max_y_for_periodic ) / 2.0_r8
   loopy%midline(IZ) = (min_z_for_periodic + max_z_for_periodic ) / 2.0_r8

   if (debug > 0) write(0,*) 'setting midlines: ', loopy%midline
   ! update these in the get_close init routines
   loopy%midbox(IX) = -1
   loopy%midbox(IY) = -1
   loopy%midbox(IZ) = -1

   loopy%offset(IX) = max_x_for_periodic - min_x_for_periodic
   loopy%offset(IY) = max_y_for_periodic - min_y_for_periodic
   loopy%offset(IZ) = max_z_for_periodic - min_z_for_periodic
   if (debug > 0) write(0,*) 'setting offsets: ', loopy%offset
else
   loopy%midline(:) = MISSING_R8
   loopy%midbox(:) = MISSING_I
   loopy%offset(:) = MISSING_R8
endif

if (debug > 0) write(0,*)  'namelist read, module initialized, loopy filled in'
if (debug > 0) write(0,*)  'any_periodic = ', any_periodic

end subroutine recompute_periodic


!---------------------------------------------------------------------------
subroutine distinct_values(in_list, mycount, values, map)
!> parse an input list of values and return:
!>  1) the count of distinct values
!>  2) the list of unique values
!>  3) the mapping of the input list to the value list
!> the values and map list should already be allocated, and be the same
!> length as the incoming list length.

real(r8), intent(in)    :: in_list(:)   !< incoming list of all values
integer,  intent(out)   :: mycount      !< count of distinct values
real(r8), intent(inout) :: values(:)    !< list of distinct values
integer,  intent(inout) :: map(:)       !< mapping of in_list to values

integer :: i, j, listsize, nextslot
logical :: foundnew
real(r8) :: newval

! set return values now; if we error out then we can
! just return.
mycount = 0
values(:) = -1.0_r8
map(:) = -1

listsize = size(in_list)
if (listsize <= 0) return
if (size(values) /= size(in_list)) return
if (size(map) /= size(in_list)) return

! set values() with only the unique distances.
! when done, the valid values are only 'count' long,
! not 'listsize'.
OUTER: do i=1, listsize
  newval = in_list(i)
  foundnew = .true.
  INNER: do j=1, listsize
     nextslot = j
     if (values(j) < 0.0_r8) exit INNER
     if (abs(values(j) - newval) <= epsilon(newval)) then
        foundnew = .false.
        map(i) = j
        exit INNER
     endif
  enddo INNER
  if (foundnew) then
     values(nextslot) = newval
     map(i) = nextslot
     mycount = nextslot
  endif
enddo OUTER

end subroutine distinct_values


!---------------------------------------------------------------------------

function get_maxdist(gc, obs_type)
type(get_close_type), intent(in) :: gc
integer, optional,    intent(in) :: obs_type
real(r8) :: get_maxdist

integer :: bt

bt = gc%type_to_cutoff_map(obs_type)
get_maxdist = gc%box(bt)%maxdist

!get_maxdist = gc%maxdist

end function get_maxdist

!----------------------------------------------------------------------------

subroutine print_get_close_type(gc, tt, amount)

! print out debugging statistics, or optionally print out a full
! dump from all mpi tasks in a format that can be plotted with matlab.
type(get_close_type), intent(in) :: gc
integer, intent(in), optional    :: tt
integer, intent(in), optional    :: amount

integer :: i, j, k, l, first, index, mytask, alltasks, whichtt
integer :: sample, nfull, nempty, howmuch, total, maxcount, maxi, maxj, maxk
logical :: tickmark(gc%box(1)%num), iam0
real(r8) :: x_cen, y_cen, z_cen

logical, save :: write_now = .true.
integer, save :: been_called = 0
integer :: funit
character(len=64) :: fname

! cumulative times through this routine
been_called = been_called + 1

! second arg is optional, defaults to 1, and selects which
! of the cutoff structs to print
if (present(tt)) then
   whichtt = tt
else
   whichtt = 1
endif

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
!if (gc%box(whichtt)%num > 100) howmuch = 0

! print the get_close_type derived type values

if (howmuch /= 0 .and. iam0) then
   write(errstring,*) 'get_close_type values:'
   call error_handler(E_MSG, 'loc', errstring)

   write(errstring,*) ' num = ', gc%box(whichtt)%num
   call error_handler(E_MSG, 'loc', errstring)

   write(errstring,*) ' nx, ny, nz = ', nx, ny, nz
   call error_handler(E_MSG, 'loc', errstring)

   write(errstring,"(A,F12.6)") ' maxdist = ', gc%box(whichtt)%maxdist
   call error_handler(E_MSG, 'loc', errstring)
   write(errstring, "(A,3(F12.6))") ' x_box: bot, top, width = ', gc%box(whichtt)%bot_x, gc%box(whichtt)%top_x, gc%box(whichtt)%x_width
   call error_handler(E_MSG, 'loc', errstring)
   write(errstring, "(A,3(F12.6))") ' y_box: bot, top, width = ', gc%box(whichtt)%bot_y, gc%box(whichtt)%top_y, gc%box(whichtt)%y_width
   call error_handler(E_MSG, 'loc', errstring)
   write(errstring, "(A,3(F12.6))") ' z_box: bot, top, width = ', gc%box(whichtt)%bot_z, gc%box(whichtt)%top_z, gc%box(whichtt)%z_width
   call error_handler(E_MSG, 'loc', errstring)

endif

! this one can be very large.   print only the first nth unless
! instructed otherwise.  (print n+1 because 1 more value fits on
! the line because it prints ( i ) and not ( i, j ) like the others.)
if (allocated(gc%box(whichtt)%loc_box)) then
   i = size(gc%box(whichtt)%loc_box,1)
   if (i/= gc%box(whichtt)%num) then
      write(errstring,*) ' warning: size of loc_box incorrect, nlocs, i =', gc%box(whichtt)%num, i
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
   if (howmuch > 1) then
      ! DEBUG
      write(errstring,"(A,I8,A,36(I8,1X))") ' loc_box(',i,') =', gc%box(whichtt)%loc_box(1:min(i,36))  ! (nlocs)
      !write(errstring,*) ' loc_box(',i,') =', gc%box%loc_box    ! (nlocs)
      call error_handler(E_MSG, 'locations_mod', errstring)
   else if(howmuch > 0) then
      write(errstring,*) ' loc_box(',i,') =', gc%box(whichtt)%loc_box(1:min(i,sample+1))
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
if (allocated(gc%box(whichtt)%start)) then
   i = size(gc%box(whichtt)%start,1)
   j = size(gc%box(whichtt)%start,2)
   k = size(gc%box(whichtt)%start,3)
   if ((i /= nx) .or. (j /= ny) .or. (k /= nz)) then
      write(errstring,*) ' warning: size of start incorrect, nx, ny, nz, i, j, k =', nx, ny, nz, i, j, k
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
   if (howmuch > 1) then
      write(errstring,*) ' start(',i,j,k,') ='              ! (nx, ny, nz)
      call error_handler(E_MSG, 'locations_mod', errstring)
      do l=1, j
         write(errstring,"(36(I8,1X))") gc%box(whichtt)%start(1:min(i,36), l, 1)
         call error_handler(E_MSG, 'locations_mod', errstring)
      enddo
   else if (howmuch > 0) then
      write(errstring,*) ' start(',i,j,k,') =', gc%box(whichtt)%start(1:min(i,sample), 1, 1)
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
if (allocated(gc%box(whichtt)%count)) then
   i = size(gc%box(whichtt)%count,1)
   j = size(gc%box(whichtt)%count,2)
   k = size(gc%box(whichtt)%count,3)
   if ((i /= nx) .or. (j /= ny) .or. (k /= nz)) then
      write(errstring,*) ' warning: size of count incorrect, nx, ny, nz, i, j, k =', nx, ny, nz, i, j, k
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
   if (howmuch > 1) then
      write(errstring,*) ' count(',i,j,k,') ='              ! (nx, ny, nz)
      call error_handler(E_MSG, 'locations_mod', errstring)
      do l=1, j
         write(errstring,"(36(I8,1X))") gc%box(whichtt)%count(1:min(i,36), l, 1)
         call error_handler(E_MSG, 'locations_mod', errstring)
      enddo
   else if (howmuch > 0) then
      write(errstring,*) ' count(',i,j,k,') =', gc%box(whichtt)%count(1:min(i,sample), 1, 1)
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
         first = gc%box(whichtt)%start(i, j, k)
         do l=1, gc%box(whichtt)%count(i, j, k)
            index = first + l - 1
            if ((index < 1) .or. (index > gc%box(whichtt)%num)) then
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

do i=1, gc%box(whichtt)%num
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
      x_cen = gc%box(whichtt)%bot_x + ((i-1)*gc%box(whichtt)%x_width) + (gc%box(whichtt)%x_width/2.0)
      write(funit, '(A,I2,A,I4,A,F12.9,A)') 'xlocs(', i, ',', mytask+1, ') = ',  x_cen, ';'
   endif
   do j=1, ny
      if (howmuch == -8 .and. i==1) then
         y_cen = gc%box(whichtt)%bot_y + ((j-1)*gc%box(whichtt)%y_width) + (gc%box(whichtt)%y_width/2.0)
         write(funit, '(A,I2,A,I4,A,F12.9,A)') 'ylocs(', j, ',', mytask+1, ') = ',  y_cen, ';'
      endif
      do k=1, nz
         if (howmuch == -8 .and. i==1) then
            z_cen = gc%box(whichtt)%bot_z + ((j-1)*gc%box(whichtt)%z_width) + (gc%box(whichtt)%z_width/2.0)
            write(funit, '(A,I2,A,I4,A,F12.9,A)') 'zlocs(', k, ',', mytask+1, ') = ',  z_cen, ';'
         endif
         if (gc%box(whichtt)%count(i, j, k) > 0) then
            nfull = nfull + 1
            total = total + gc%box(whichtt)%count(i, j, k)
            if (gc%box(whichtt)%count(i, j, k) > maxcount) then
               maxcount = gc%box(whichtt)%count(i, j, k)
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
                                   ',', mytask+1, ') = ', gc%box(whichtt)%count(i, j, k), ';'
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
write(errstring, '(a,i9)') " Total items to put in boxes: ", gc%box(whichtt)%num
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

subroutine exhaustive_collect(box, base_loc, base_type, loc_list, num_close, close_ind, close_dist)

! For validation, it is useful to be able to compare against exact
! exhaustive search

type(box_type),        intent(in)  :: box
type(location_type),   intent(in)  :: base_loc, loc_list(:)
integer,               intent(in)  :: base_type ! JDL Addition
integer,               intent(out) :: num_close
integer,               intent(out) :: close_ind(:)
real(r8),              intent(out) :: close_dist(:)

real(r8) :: this_dist
integer :: i

num_close = 0
do i = 1, box%num
   this_dist = get_dist(base_loc, loc_list(i), type1=base_type)
   if(this_dist <= box%maxdist) then
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

integer,            intent(in)    :: cnum_close, num_close
integer,            intent(inout) :: cclose_ind(:)
integer,            intent(in)    :: close_ind(:)
real(r8),           intent(inout) :: cclose_dist(:)
real(r8), optional, intent(in)    :: close_dist(:)

! Do comparisons against full search
if((num_close /= cnum_close) .and. present(close_dist)) then
   write(errstring, *) 'get_close (', num_close, ') should equal exhaustive search (', cnum_close, ')'
   call error_handler(E_ERR, 'get_close_obs', errstring, source, &
                      text2='optional arg "dist" is present; we are computing exact distances', &
                      text3='the exhaustive search should find an identical number of locations')
else if (num_close < cnum_close) then
   write(errstring, *) 'get_close (', num_close, ') should not be smaller than exhaustive search (', cnum_close, ')'
   call error_handler(E_ERR, 'get_close_obs', errstring, source, &
                      text2='optional arg "dist" not present; we are returning a superset of close locations', &
                      text3='the exhaustive search should find an equal or lesser number of locations')
endif

! if they do not compare, we have the exhaustive lists here and can print out
! exactly which items in the list differ.

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
! end of location/threed_cartesian/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod

