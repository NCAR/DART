! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

!> Interpolation routines for longitude/latitude grids which are logically 
!> rectangular and either fully regular, partially regular or fully deformed.
!>
!> This module includes initialization, search, interpolation, and finalization 
!> routines. 
!> 
!> The size of the grid is specified by a count of the longitudes and latitudes.
!>
!> The actual coordinates of the grid can be specified in one of 3 ways:
!>
!>    fully regular grid, evenly spaced and fully orthogonal:
!>     (origin, delta) each for lon and lat.
!>
!>    fully orthogonal but possibly irregularly spaced along the axes:
!>     1D array(counts) each for lon and lat.
!>
!>    logically rectangular but corners of quads are fully irregular:
!>     2D array(lon counts, lat counts) each for lon and lat.
!>
!> To search for a given location and return the (i,j) indices of the
!> corners of the enclosing quad:
!>
!>    for the fully regular grid the enclosing quad is found by computing
!>    the corresponding (i,j) indices.
!>
!>    for the irregularly spaced grid, the enclosing quad is found by searching
!>    each of the 1D arrays for the enclosing (i,j) indices.
!>
!>    For the fully irregular grid, the enclosing quad is found by this process:
!>     At setup time:
!>     1) create a coarse fully regular grid.
!>     2) compute the intersection of each regular grid box with the target grid quads
!>     3) keep a list of which target grid quads overlap for each regular grid box.
!>     At search time:
!>     4) find the location in the fully regular grid (which can be done quickly)
!>     5) do an exhaustive search of the target grid quads which overlap that 
!>        regular grid box, returning when one of the target grid quads encloses 
!>        the given location.
!>

module quad_utils_mod

use        types_mod, only : r8, i8, MISSING_R8, PI, deg2rad

use     location_mod, only : location_type, get_location

use    utilities_mod, only : register_module, error_handler,         &
                             E_ERR, E_WARN, E_MSG, nmlfileunit,      &
                             do_output, do_nml_file, do_nml_term,    &
                             find_namelist_in_file, check_namelist_read, &
                             log_it, array_dump

use      options_mod, only : get_missing_ok_status

implicit none
private

public :: quad_interp_handle,              & ! derived type which holds the grid and option info
          init_quad_interp,                & ! pass in grid type and counts here
          finalize_quad_interp,            & ! release storage 
          set_quad_coords,                 & ! overload these 3: set_reg_xx, set_1d_xx, set_2d_xx
          quad_lon_lat_locate,             & ! given lat,lon return above and below Is and Js
          quad_lon_lat_evaluate,           & ! given i,j and all 4 corner values, return interp val
          GRID_QUAD_FULLY_REGULAR,         &
          GRID_QUAD_IRREG_SPACED_REGULAR,  &
          GRID_QUAD_FULLY_IRREGULAR,       &
          GRID_QUAD_UNKNOWN_TYPE,          &
          QUAD_LOCATED_UNKNOWN,            &
          QUAD_LOCATED_CELL_CENTERS,       &
          QUAD_LOCATED_LON_EDGES,          &
          QUAD_LOCATED_LAT_EDGES,          &
          QUAD_LOCATED_CELL_CORNERS,       &
          get_quad_grid_size,              &
          get_quad_global,                 &
          print_quad_handle                ! debug


! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

! message strings
character(len=512) :: string1, string2, string3

logical, save :: module_initialized = .false.

logical :: missing_ok_in_state

integer  :: debug = 0               ! turn up for more and more debug messages
integer  :: interpolation_type = 1  ! add cases for different strategies
logical  :: do_rotate = .false.     ! rotate edge from pts 1,2 to horizontal before interp

namelist /quad_interpolate_nml/ do_rotate, debug

!> @todo FIXME internal routines could use h for the handle; externally callable
!> routines should use interp_handle for clarity in the interface.

! the grid must always be logically rectangular, so knowing the i and j of a
! quad corner, the next quad starts at index i+1 and j+1.

! 2d grid types:
!  can lats/lons be defined by giving only start, delta?  type 1
!  are they each 1d arrays and the grid is defined by the cross product?  type 2
!  are lat and lon both full 2d arrays, so completely irregular?  type 3
integer, parameter :: GRID_QUAD_FULLY_REGULAR        =  1
integer, parameter :: GRID_QUAD_IRREG_SPACED_REGULAR =  2
integer, parameter :: GRID_QUAD_FULLY_IRREGULAR      =  3
integer, parameter :: GRID_QUAD_UNKNOWN_TYPE         = -1

! where the locations are relative to each grid cell
integer, parameter :: QUAD_LOCATED_UNKNOWN       =  -1
integer, parameter :: QUAD_LOCATED_CELL_CENTERS  =   1
integer, parameter :: QUAD_LOCATED_LON_EDGES     =   2
integer, parameter :: QUAD_LOCATED_LAT_EDGES     =   3
integer, parameter :: QUAD_LOCATED_CELL_CORNERS  =   4

! options that control how the interpolation handles specific cases
type quad_grid_options
   private

   ! if any corner of the cell is masked to false, fail to find the location.
   logical :: uses_mask = .false.
   logical, allocatable :: grid_mask(:,:)

   ! global grid?  if yes, both spans_lon_zero and pole_wrap are true as well.
   ! if not, either or both could still be true but we can't assume.
   logical :: global_grid = .false.

   ! separate out the single logical 'spans' flag into two cases?
   ! case 1:  the longitude grid is cyclic; either a global grid or a band between 
   !  two latitude lines that circles the globe.  all longitude values are valid.
   ! case 2: a regional grid that crosses the prime meridian.
   !  this will contain a discontinuity around 360 -> 0 that should be a valid region.
   ! can a single flag handle both of these cases?  

   ! are there valid values between array(N) and array(1) that should be interpolated?
   ! always true for global grids. 
   logical :: lon_cyclic = .false.  

   ! are there valid values between array(X) and array(X+1) in the interior of
   ! the longitude array where a(X) > a(X+1) (e.g. around where a(X) is near 360 
   ! and a(X+1) is near 0) that should be interpolated?

   ! always true for global grids. 
   ! for partially regular grids (1D lon and 1D lat arrays) this can be
   ! true if a regional grid crosses the prime meridian.
   ! for irregular grids (2D lons and 2D lats) this is true if 
   ! any lon(:,1) > lon(:,nlons)
   logical :: spans_lon_zero = .false.   

   ! do we handle wrap over the poles?  
   ! (complicated for interpolating vector values.)
   logical :: pole_wrap = .false.

   ! are the latitudes specified from smallest (south) to largest (north)
   ! or largest (north) to smallest (south)?
   logical :: north_to_south = .false.

   ! i don't want to know this, but apparently we might
   ! have to know if the points given are cell-centered or
   ! on the edges.  (god forbid we need to know U stagger vs V stagger!)
   ! (for detecting if we are at the poles, for example)
   integer :: cell_relative = QUAD_LOCATED_CELL_CENTERS
   ! or should default be UNKNOWN, or CORNERS ?

end type quad_grid_options


! fully regular grid, evenly spaced and fully orthogonal
type quad_reg_grid_coords
   private

   real(r8) :: lat_start = MISSING_R8
   real(r8) :: lat_delta = MISSING_R8
   real(r8) :: lon_start = MISSING_R8
   real(r8) :: lon_delta = MISSING_R8

end type quad_reg_grid_coords

! fully orthorgonal but irregularly spaced along the axis
type quad_irregspaced_grid_coords
   private

   real(r8), allocatable :: lats_1D(:)
   real(r8), allocatable :: lons_1D(:)

end type quad_irregspaced_grid_coords

! logically rectangular but corners of quads are fully irregular
type quad_irreg_grid_coords
   private

   ! These arrays store the longitude and latitude of the lower left corner of
   ! each of the quadrilaterals and must be set by the user.

   real(r8), allocatable :: lats_2D(:,:)
   real(r8), allocatable :: lons_2D(:,:)

   ! the sizes of these depend on the grid size.  these are good defaults for ?
   integer  :: num_reg_x = 180
   integer  :: num_reg_y = 180
   integer  :: max_reg_list_num = 800
   real(r8) :: min_lon =     0.0_r8
   real(r8) :: max_lon =   360.0_r8
   real(r8) :: lon_width = 360.0_r8
   real(r8) :: min_lat =   -90.0_r8
   real(r8) :: max_lat =    90.0_r8
   real(r8) :: lat_width = 180.0_r8

   ! these next 2 should be allocated num_reg_x, num_reg_y
   integer, allocatable :: grid_start(:,:)
   integer, allocatable :: grid_num  (:,:)
   integer, allocatable :: grid_lon_list(:)
   integer, allocatable :: grid_lat_list(:)
   ! these replace u_dipole_xxx and t_dipole_xxx

end type quad_irreg_grid_coords


! public type -- derived type to hold search and interp info
type quad_interp_handle
   private

   integer :: grid_type = GRID_QUAD_UNKNOWN_TYPE

   integer :: nlat = -1   ! grid sizes in each dim
   integer :: nlon = -1

   ! which ones of these are allocated/used depends on the grid type

   type(quad_reg_grid_coords)         :: rr
   type(quad_irregspaced_grid_coords) :: ir
   type(quad_irreg_grid_coords)       :: ii

   type(quad_grid_options) :: opt

end type quad_interp_handle

interface quad_lon_lat_locate
   module procedure quad_lon_lat_locate_ii
   module procedure quad_lon_lat_locate_ir    ! handles both ir and rr
end interface

interface quad_lon_lat_evaluate
   module procedure quad_lon_lat_evaluate_ii_single
   module procedure quad_lon_lat_evaluate_ii_array
   module procedure quad_lon_lat_evaluate_ir_single  ! handles both ir and rr
   module procedure quad_lon_lat_evaluate_ir_array   ! handles both ir and rr
end interface

interface set_quad_coords
   module procedure set_reg_quad_coords
   module procedure set_irregspaced_quad_coords
   module procedure set_irreg_quad_coords
end interface


!------------------------------------------------

! NOTE (dipole/tripole grids): since both of the dipole and tripole
! grids are logically rectangular we can use the same interpolation
! scheme originally implemented for the dipole grid. Here we can
! interchange dipole and tripole when reading the code.

! The regular grid used for dipole interpolation divides the sphere into
! a set of regularly spaced lon-lat boxes. The number of boxes in
! longitude and latitude are set by num_reg_x and num_reg_y. Making the
! number of regular boxes smaller decreases the computation required for
! doing each interpolation but increases the static storage requirements
! and the initialization computation (which seems to be pretty small).
! FIX ME: to account for various grid sizes we should dynamically
! allocate these numbers.  To keep max_reg_list_num < 100 we can use:
!    tx0.1v2 num_reg_x = num_reg_y = 900
!    tx0.5v1 num_reg_x = num_reg_y = 180
!    gx1v6   num_reg_x = num_reg_y = 90
! Larger num_reg_(x,y) values require more temporary storage in
! ureg_list_lon, ureg_list_lat, treg_list_lon, treg_list_lat. For now
! we can use num_reg_(x,y) = 180 and max_reg_list_num = 800 to account
! for all of the currently implemented grid types.
!integer, parameter :: num_reg_x = 180, num_reg_y = 180

! The max_reg_list_num controls the size of temporary storage used for
! initializing the regular grid. Four arrays
! of size num_reg_x*num_reg_y*max_reg_list_num are needed. The initialization
! fails and returns an error if max_reg_list_num is too small. With 180 regular
! lat lon boxes a value of 30 is sufficient for the gx3 POP grid, 80 for the
! gx1 grid, 180 for the tx0.5 grid and 800 for the tx0.1 grid.
! FIX ME: we should declare this at runtime depending on the grid size.
!integer, parameter :: max_reg_list_num = 800

! The dipole interpolation keeps a list of how many and which dipole quads
! overlap each regular lon-lat box. The number for the u and t grids are
! stored in u_dipole_num and t_dipole_num. The allocatable arrays
! u_dipole_lon(lat)_list and t_dipole_lon(lat)_list list the longitude
! and latitude indices for the overlapping dipole quads. The entry in
! u_dipole_start and t_dipole_start for a given regular lon-lat box indicates
! where the list of dipole quads begins in the u_dipole_lon(lat)_list and
! t_dipole_lon(lat)_list arrays.

! Need to check for pole quads: for now we are not interpolating in them
integer :: pole_x, t_pole_y, u_pole_y

! Have a global variable saying whether this is dipole or regular lon-lat grid
! This should be initialized static_init_model. Code to do this is below.

! FIXME: remove this - it becomes part of the interp handle
!logical :: dipole_grid


contains

!------------------------------------------------------------------
!------------------------------------------------------------------

!> Called to do one time initialization of the module.

subroutine initialize_module()

integer :: iunit, io

if (module_initialized) return

module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'quad_interpolate_nml', iunit)
read(iunit, nml = quad_interpolate_nml, iostat = io)
call check_namelist_read(iunit, io, 'quad_interpolate_nml')

if (do_nml_file()) write(nmlfileunit, nml=quad_interpolate_nml)
if (do_nml_term()) write(     *     , nml=quad_interpolate_nml)

! are MISSING_R8 values possible in the model state?
missing_ok_in_state = get_missing_ok_status()

end subroutine initialize_module

!------------------------------------------------------------------

!> initialize a grid handle with the grid type and the sizes of each axis.
!> after this other routines are called based on the grid type to specify
!> the actual lat/lon arrays.

!>@todo FIXME cell_relative matters for poles and longitude wrap - should
!> it be optional since regional areas not near 0 lon don't care?  maybe it
!> is better to be explicit and just say CELL_CENTERS is a good default.

subroutine init_quad_interp(grid_type, num_lons, num_lats, cell_relative, &
                            global, spans_lon_zero, pole_wrap, interp_handle)

integer,                  intent(in)  :: grid_type
integer,                  intent(in)  :: num_lons
integer,                  intent(in)  :: num_lats
integer,                  intent(in)  :: cell_relative
logical,                  intent(in)  :: global, spans_lon_zero, pole_wrap
type(quad_interp_handle), intent(out) :: interp_handle

if (.not. module_initialized) call initialize_module()

interp_handle%nlat = num_lats
interp_handle%nlon = num_lons

!>@todo  : add sanity checking between global, spans_lon_zero, pole_wrap settings

if (global) then
   interp_handle%opt%global_grid = .true.
   interp_handle%opt%spans_lon_zero = .true.
   interp_handle%opt%pole_wrap = .true.
else
   if (spans_lon_zero) interp_handle%opt%spans_lon_zero = .true.
   if (pole_wrap) interp_handle%opt%pole_wrap = .true.
endif

interp_handle%grid_type = grid_type

select case (grid_type)
   case(GRID_QUAD_FULLY_REGULAR)
      ! nothing to do for this case

   case(GRID_QUAD_IRREG_SPACED_REGULAR)
      allocate(interp_handle%ir%lats_1D(num_lats), &
               interp_handle%ir%lons_1D(num_lons))
      interp_handle%ir%lats_1D(num_lats) = MISSING_R8
      interp_handle%ir%lons_1D(num_lons) = MISSING_R8
      if (debug > 10) then
         write(string1, *) 'nlats, nlons: ', num_lats, num_lons
         call log_it(string1)
         write(string1, *) 'first lat: ', interp_handle%ir%lats_1D(1)
         call log_it(string1)
         write(string1, *) 'first lon: ', interp_handle%ir%lons_1D(1)
         call log_it(string1)
      endif

   case(GRID_QUAD_FULLY_IRREGULAR)
      allocate(interp_handle%ii%lats_2D(num_lons, num_lats), &
               interp_handle%ii%lons_2D(num_lons, num_lats))
      interp_handle%ii%lats_2D(num_lons, num_lats) = MISSING_R8
      interp_handle%ii%lons_2D(num_lons, num_lats) = MISSING_R8

      ! tx0.1v2 num_reg_x = num_reg_y = 900
      ! tx0.5v1 num_reg_x = num_reg_y = 180
      ! gx1v6   num_reg_x = num_reg_y = 90
      ! max_reg_list_num = 800

      ! adjust num_regs here based on numlons, numlats
      if (num_lats * num_lons > 6 * 1000 * 1000) then  ! ~1/10th degree
         interp_handle%ii%num_reg_x = 900
         interp_handle%ii%num_reg_y = 900
         interp_handle%ii%max_reg_list_num = 800   !todo  what is good val?
         if(debug > 10) then
            write(string1, *) 'case 1: ', interp_handle%ii%num_reg_x, interp_handle%ii%num_reg_y, &
                               interp_handle%ii%max_reg_list_num
            call log_it(string1)
         endif

      else if (num_lats * num_lons > 250 * 1000) then  ! ~1/2th degree
         interp_handle%ii%num_reg_x = 180
         interp_handle%ii%num_reg_y = 180
         interp_handle%ii%max_reg_list_num = 800
         if(debug > 10) then
            write(string1, *) 'case 2: ', interp_handle%ii%num_reg_x, interp_handle%ii%num_reg_y, &
                               interp_handle%ii%max_reg_list_num
            call log_it(string1)
         endif

      else
         interp_handle%ii%num_reg_x = 90
         interp_handle%ii%num_reg_y = 90
         interp_handle%ii%max_reg_list_num = 800
         if(debug > 10) then
            write(string1, *) 'case 3: ', interp_handle%ii%num_reg_x, interp_handle%ii%num_reg_y, &
                               interp_handle%ii%max_reg_list_num
            call log_it(string1)
         endif

      endif

      allocate(interp_handle%ii%grid_start(interp_handle%ii%num_reg_x, &
                                           interp_handle%ii%num_reg_y))
      allocate(interp_handle%ii%grid_num(  interp_handle%ii%num_reg_x, &
                                           interp_handle%ii%num_reg_y))

      interp_handle%ii%grid_num = 0

   case default
      write(string1, *) 'unrecognized grid type: ', grid_type
      write(string2, *) 'should be one of: GRID_QUAD_FULLY_REGULAR, ', &
                        'GRID_QUAD_IRREG_SPACED_REGULAR, GRID_QUAD_FULLY_IRREGULAR'
      call error_handler(E_ERR, 'init_quad_interp', string1, &
                         source, revision, revdate, text2=string2)

end select

select case (cell_relative)
   case (QUAD_LOCATED_CELL_CENTERS, QUAD_LOCATED_LON_EDGES, QUAD_LOCATED_LAT_EDGES, &
         QUAD_LOCATED_CELL_CORNERS)
      interp_handle%opt%cell_relative = cell_relative

   case default
      write(string1, *) 'unrecognized cell_relative type: ', cell_relative
      write(string2, *) 'should be one of: QUAD_LOCATED_CELL_CENTERS, ', &
            'QUAD_LOCATED_LON_EDGES, QUAD_LOCATED_LAT_EDGES, QUAD_LOCATED_CELL_CORNERS'
      write(string3, *) 'important if handling poles and/or longitude wrap across prime meridian'
      call error_handler(E_ERR, 'init_quad_interp', string1, &
                         source, revision, revdate, text2=string2, text3=string3)

end select

if (debug > 2) then
   write(string1, *) 'calling init for nlons/nlats/type = ', num_lons, num_lats, grid_type
   call log_it(string1)
endif

end subroutine init_quad_interp

!------------------------------------------------------------------

subroutine print_quad_handle(interp_handle)
type(quad_interp_handle), intent(in) :: interp_handle

if (debug > 10) then
   write(string1, *) 'nlat, nlon, grid type: ', interp_handle%nlat, interp_handle%nlon, &
                      interp_handle%grid_type
   call log_it(string1)
endif

select case (interp_handle%grid_type)
   case(GRID_QUAD_FULLY_REGULAR)
      call log_it('fully regular quad grid')
      write(string1, *) 'lon start, delta, count: ', interp_handle%rr%lon_start, &
                         interp_handle%rr%lon_delta, interp_handle%nlon
      call log_it(string1)
      write(string1, *) 'lat start, delta, count: ', interp_handle%rr%lat_start, &
                         interp_handle%rr%lat_delta, interp_handle%nlat
      call log_it(string1)

   case(GRID_QUAD_IRREG_SPACED_REGULAR)
      call log_it('irregularly spaced but orthogonal quad grid')
      write(string1, *) 'nlons: ', interp_handle%nlon
      call log_it(string1)
      call array_dump(interp_handle%ir%lons_1D, label = 'lon values')
      write(string1, *) 'nlats: ', interp_handle%nlat
      call log_it(string1)
      call array_dump(interp_handle%ir%lats_1D, label = 'lat values')

   case(GRID_QUAD_FULLY_IRREGULAR)
      call log_it('fully irregular quad grid')
      write(string1, *) 'nlons: ', interp_handle%nlon
      call log_it(string1)
      call array_dump(interp_handle%ii%lons_2D, label = 'lon values')
      write(string1, *) 'nlats: ', interp_handle%nlat
      call log_it(string1)
      call array_dump(interp_handle%ii%lats_2D, label = 'lat values')

   case default
      write(string1, *) 'unrecognized grid type: ', interp_handle%grid_type
      write(string2, *) 'should be one of: GRID_QUAD_FULLY_REGULAR, '&
                        &'GRID_QUAD_IRREG_SPACED_REGULAR, GRID_QUAD_FULLY_IRREGULAR'
      call error_handler(E_ERR, 'print_quad_handle', string1, &
                         source, revision, revdate, text2=string2)

end select

if (debug > 10) then
   write(string1, *) 'cell relative flag: ', interp_handle%opt%cell_relative
   call log_it(string1)
endif

end subroutine print_quad_handle

!------------------------------------------------------------------

subroutine finalize_quad_interp(interp_handle)

type(quad_interp_handle), intent(inout) :: interp_handle

! reset vals and deallocate storage

!>@todo FIXME: make this call individual subtype destructors?

interp_handle%nlat = -1
interp_handle%nlon = -1

interp_handle%grid_type = GRID_QUAD_UNKNOWN_TYPE

if (allocated(interp_handle%ir%lats_1D)) deallocate(interp_handle%ir%lats_1D)
if (allocated(interp_handle%ir%lons_1D)) deallocate(interp_handle%ir%lons_1D)

if (allocated(interp_handle%ii%lats_2D)) deallocate(interp_handle%ii%lats_2D)
if (allocated(interp_handle%ii%lons_2D)) deallocate(interp_handle%ii%lons_2D)
if (allocated(interp_handle%ii%grid_start)) deallocate(interp_handle%ii%grid_start)
if (allocated(interp_handle%ii%grid_num)) deallocate(interp_handle%ii%grid_num)
if (allocated(interp_handle%ii%grid_lon_list)) deallocate(interp_handle%ii%grid_lon_list)
if (allocated(interp_handle%ii%grid_lat_list)) deallocate(interp_handle%ii%grid_lat_list)

interp_handle%opt%cell_relative = QUAD_LOCATED_UNKNOWN

end subroutine finalize_quad_interp

!------------------------------------------------------------

subroutine set_reg_quad_coords(interp_handle, lon_start, lon_delta, &
                                              lat_start, lat_delta)

type(quad_interp_handle), intent(inout) :: interp_handle
real(r8),                 intent(in)    :: lon_start, lon_delta
real(r8),                 intent(in)    :: lat_start, lat_delta

interp_handle%rr%lon_start = lon_start
interp_handle%rr%lon_delta = lon_delta
interp_handle%rr%lat_start = lat_start
interp_handle%rr%lat_delta = lat_delta

if (lon_delta == 0.0_r8 .or. lat_delta == 0.0_r8) then
   write(string1, *) 'neither lon_delta nor lat_delta can equal 0'
   write(string2, *) 'lon_delta: ', lon_delta, ' lat_delta: ', lat_delta
   call error_handler(E_ERR, 'set_quad_coords', string1, &
                      source, revision, revdate, text2=string2)
endif

end subroutine set_reg_quad_coords

!------------------------------------------------------------

subroutine set_irregspaced_quad_coords(interp_handle, lons, lats)

type(quad_interp_handle), intent(inout) :: interp_handle
real(r8),                 intent(in)    :: lons(:)
real(r8),                 intent(in)    :: lats(:)

if (size(lons) /= interp_handle%nlon) then
   write(string1, *) 'longitude count in handle: ', interp_handle%nlon, &
                       ' must match length of 1D lons array: ', size(lons)
   call error_handler(E_ERR, 'set_irregspaced_quad_coords', string1, &
                      source, revision, revdate)
endif

if (size(lats) /= interp_handle%nlat) then
   write(string1, *) 'latitude count in handle: ', interp_handle%nlat, &
                       ' must match length of 1D lats array: ', size(lats)
   call error_handler(E_ERR, 'set_irregspaced_quad_coords', string1, &
                      source, revision, revdate)
endif

! lons and lats are declared intent(in).  any code that just needs to
! test values can use them.  any code that is going to modify values
! has to use the xxx_1D arrays in the structures.

interp_handle%ir%lons_1D(:) = lons
interp_handle%ir%lats_1D(:) = lats

! inverted order for latitudes?
if (lats(1) > lats(interp_handle%nlat)) then
   interp_handle%opt%north_to_south = .true.
endif

! -180 to 180 instead of 0 to 360?  add 360, which makes all
! the values valid (we require longitudes to be between 0 and 360)
! but it also makes a partially regular grid start at 180, go up
! to 360, then back to 0, then up to 180.  set the 'spans 0' flag.
if (any(lons < 0.0_r8)) then
  where(interp_handle%ir%lons_1D < 0.0_r8) &
     interp_handle%ir%lons_1D = interp_handle%ir%lons_1D + 360.0_r8
     interp_handle%opt%spans_lon_zero = .true. 
endif

! validate ranges 
if (any(interp_handle%ir%lats_1D < -90.0_r8) .or. any(interp_handle%ir%lats_1D > 90.0_r8)) then
   write(string1, *) 'latitude values must be between -90 and 90.', &
                       ' out of range values found in latitude array. '
   write(string2, *) 'min, max values: ', minval(interp_handle%ir%lats_1D), &
                                          maxval(interp_handle%ir%lats_1D)
   call error_handler(E_ERR, 'set_irregspaced_quad_coords', string1, &
                      source, revision, revdate, text2=string2)
endif

if (any(interp_handle%ir%lons_1D < 0.0_r8) .or. any(interp_handle%ir%lons_1D > 360.0_r8)) then
   write(string1, *) 'longitude values must be between 0 and 360, or -180 and 180.', &
                       ' out of range values found in longitude array. '
   write(string2, *) 'min, max values: ', minval(interp_handle%ir%lons_1D), &
                                          maxval(interp_handle%ir%lons_1D)
   call error_handler(E_ERR, 'set_irregspaced_quad_coords', string1, &
                      source, revision, revdate, text2=string2)
endif

!>@todo FIXME i would like to put something like this to check
!>for degenerate grids, but i don't know how to avoid throwing
!>an error at the poles, for example.  i'm leaving this here
!>but commented out to remind me to try to add some way of 
!>catching bad values at init time.
!do i=1, nlons-1
!   lon_delta = interp_handle%ir%lons_1d(i+1) - interp_handle%ir%lons_1d(i)
!   if (lon_delta == 0.0_r8) then
!      write(string1, *) 'no lon_deltas can equal 0'
!      write(string2, *) 'i, lons_1d(i), lons_1d(i+1): ', i, lons_1d(i), lons_1d(i+1)
!      call error_handler(E_ERR, 'set_quad_coords', string1, &
!                         source, revision, revdate, text2=string2)
!   endif
!enddo
!do j=1, nlats-1
!   lat_delta = interp_handle%ir%lats_1d(j+1) - interp_handle%ir%lats_1d(j)
!   if (lat_delta == 0.0_r8) then
!      write(string1, *) 'no lat_deltas can equal 0'
!      wrjte(string2, *) 'j, lats_1d(j), lats_1d(j+1): ', j, lats_1d(j), lats_1d(j+1)
!      call error_handler(E_ERR, 'set_quad_coords', string1, &
!                         source, revision, revdate, text2=string2)
!   endif
!enddo

end subroutine set_irregspaced_quad_coords

!------------------------------------------------------------

subroutine set_irreg_quad_coords(interp_handle, lons, lats, mask)

type(quad_interp_handle), intent(inout) :: interp_handle
real(r8),                 intent(in)    :: lons(:,:)
real(r8),                 intent(in)    :: lats(:,:)
logical, optional,        intent(in)    :: mask(:,:)

integer :: gridsize(2)

gridsize = shape(lons)
call shapecheck(interp_handle, gridsize, 'longitude')

gridsize = shape(lats)
call shapecheck(interp_handle, gridsize, 'latitude')

interp_handle%ii%lons_2D(:,:) = lons(:,:)
interp_handle%ii%lats_2D(:,:) = lats(:,:)

if (present(mask)) then
   interp_handle%opt%uses_mask = .true.
   gridsize = shape(mask)
   call shapecheck(interp_handle, gridsize, 'mask')

   allocate(interp_handle%opt%grid_mask(gridsize(1), gridsize(2)))
   interp_handle%opt%grid_mask(:,:) = mask(:,:)
endif

! Initialize the interpolation routines
call init_irreg_interp(interp_handle)

end subroutine set_irreg_quad_coords

!------------------------------------------------------------

subroutine shapecheck(h, gridsize, name)

type(quad_interp_handle), intent(in) :: h
integer,                  intent(in) :: gridsize(2)
character(len=*),         intent(in) :: name

if (gridsize(1) /= h%nlon .or. gridsize(2) /= h%nlat) then
   write(string1, *) 'longitude/latitude counts in handle: ', h%nlon, h%nlat, &
                       ' must match shape of 2D '//trim(name)//' array: ', gridsize
   call error_handler(E_ERR, 'shapecheck', string1, source, revision, revdate)
endif

end subroutine shapecheck

!------------------------------------------------------------
!> Build the data structure for interpolation for an irregular quad grid

subroutine init_irreg_interp(h)

type(quad_interp_handle), intent(inout) :: h

character(len=*), parameter :: routine = 'init_irreg_interp'

! Need a temporary data structure to build this.
! These arrays keep a list of the x and y indices of dipole quads
! that potentially overlap the regular boxes.
integer, allocatable :: reg_list_lon(:,:,:)
integer, allocatable :: reg_list_lat(:,:,:)

real(r8) :: u_c_lons(4), u_c_lats(4), pole_row_lon
integer  :: i, j, k, pindex, nx, ny, nrx, nry, istatus
integer  :: reg_lon_ind(2), reg_lat_ind(2), u_total, u_index
logical  :: cyclic, pole
integer  :: xlim

allocate(reg_list_lon(h%ii%num_reg_x, h%ii%num_reg_y, h%ii%max_reg_list_num))
allocate(reg_list_lat(h%ii%num_reg_x, h%ii%num_reg_y, h%ii%max_reg_list_num))

! poles?  span?
cyclic = h%opt%spans_lon_zero
pole   = h%opt%pole_wrap
nx  = h%nlon
ny  = h%nlat
nrx = h%ii%num_reg_x
nry = h%ii%num_reg_y

reg_list_lon(:, :, :) = 0
reg_list_lat(:, :, :) = 0

! for a global grid, the initial values have already been set in
! the derived type.  otherwise, find the min/max of lons and lats.
if (.not. h%opt%global_grid) then
   h%ii%min_lon = minval(h%ii%lons_2d)
   h%ii%max_lon = maxval(h%ii%lons_2d)
   h%ii%lon_width = h%ii%max_lon - h%ii%min_lon  ! FIXME: wrap?

   if (h%ii%lon_width < 0) then
      if(h%opt%spans_lon_zero) then
         h%ii%lon_width = h%ii%lon_width + 360.0_r8
      else
         write(string1,*)'min_lon, max_lon, lon_width, spans_lon_zero: ', &
                    h%ii%min_lon, h%ii%max_lon, h%ii%lon_width, h%opt%spans_lon_zero
         call error_handler(E_ERR,routine,'regional grid with bad longitudes', &
                          source, revision, revdate, text2=string1)
      endif
   endif

   h%ii%min_lat = minval(h%ii%lats_2d)
   h%ii%max_lat = maxval(h%ii%lats_2d)
   h%ii%lat_width = h%ii%max_lat - h%ii%min_lat
endif

if (cyclic) then
   ! Begin by finding the quad that contains the pole for the dipole t_grid.
   ! To do this locate the u quad with the pole on its right boundary. This is on
   ! the row that is opposite the shifted pole and exactly follows a lon circle.
   pole_x = nx / 2;
   ! Search for the row at which the longitude flips over the pole
   pole_row_lon = h%ii%lons_2d(pole_x, 1);
   do i = 1, ny
      pindex = i
      if(h%ii%lons_2d(pole_x, i) /= pole_row_lon) exit
   enddo

   ! Pole boxes for u have indices pole_x or pole_x-1 and index - 1;
   ! (it's right before the flip).
   u_pole_y = pindex - 1;

   ! Locate the T dipole quad that contains the pole.
   ! We know it is in either the same lat quad as the u pole edge or one higher.
   ! Figure out if the pole is more or less than halfway along
   ! the u quad edge to find the right one.
   if(h%ii%lats_2d(pole_x, u_pole_y) > h%ii%lats_2d(pole_x, u_pole_y + 1)) then
      t_pole_y = u_pole_y;
   else
      t_pole_y = u_pole_y + 1;
   endif
endif

! Loop through each of the dipole grid quads
if (cyclic) then
   xlim = nx
else
   xlim = nx - 1
endif

do i = 1, xlim
   ! There's no wraparound in y, one box less than grid boundaries
   do j = 1, ny - 1

      if( all_corners_valid(h%opt, i,j, nx) ) then

 !>@todo is istatus /= 0 a failure condition

         ! Set up array of lons and lats for the corners of these u quads
         call get_quad_corners(h%ii%lons_2d, i, j, cyclic, pole, nx, ny, u_c_lons, istatus)
         if (istatus /= 0) print *, 'get_quad_corners for lons returns failure'

         call get_quad_corners(h%ii%lats_2d, i, j, cyclic, pole, nx, ny, u_c_lats, istatus)
         if (istatus /= 0) print *, 'get_quad_corners for lats returns failure'

         !print *, 'get_quad_corners returns ', u_c_lons, u_c_lats, ' for ', &
         !          h%ii%lons_2d(i,j), h%ii%lats_2d(i,j), ' index ', i, j

         ! Get list of regular boxes that cover this u dipole quad
         ! false indicates that for the u grid there's nothing special about pole
         call reg_box_overlap(h, u_c_lons, u_c_lats, .false., reg_lon_ind, reg_lat_ind)
         ! Update the temporary data structures for the u quad
         call update_reg_list(h%ii%grid_num, reg_list_lon, reg_list_lat, &
                    reg_lon_ind, reg_lat_ind, nrx, nry, h%ii%max_reg_list_num, i, j)
      endif

   enddo
enddo

write(string1,*)'to determine (minimum) max_reg_list_num values for new grids ...'
write(string2,*)'interp_handle%ii%grid_num is ',maxval(h%ii%grid_num)
call error_handler(E_MSG, routine, string1, text2=string2)

! Invert the temporary data structure. The total number of entries will be
! the sum of the number of dipole cells for each regular cell.
u_total = sum(h%ii%grid_num)

! Allocate storage for the final structures in module storage
allocate(h%ii%grid_lon_list(u_total), h%ii%grid_lat_list(u_total))

! Fill up the long list by traversing the temporary structure. Need indices
! to keep track of where to put the next entry.
u_index = 1
! Loop through each regular grid box
do i = 1, h%ii%num_reg_x
   do j = 1, h%ii%num_reg_y

      ! The list for this regular box starts at the current indices.
      h%ii%grid_start(i, j) = u_index

      ! Copy all the close dipole quads for regular u box(i, j)
      do k = 1, h%ii%grid_num(i, j)
         h%ii%grid_lon_list(u_index) = reg_list_lon(i, j, k)
         h%ii%grid_lat_list(u_index) = reg_list_lat(i, j, k)
         u_index = u_index + 1
      enddo

   enddo
enddo

! Confirm that the indices come out okay as debug
if(u_index /= u_total + 1) then
   string1 = 'Storage indices did not balance for U grid: : contact DART developers'
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
endif

deallocate(reg_list_lon, reg_list_lat)

end subroutine init_irreg_interp

!------------------------------------------------------------

!> @todo FIXME: this is the original code, for reference.
!> the init_irreg_interp() routine above should replace it.

!%! subroutine init_dipole_interp()
!%!
!%! ! Build the data structure for interpolation for a dipole grid.
!%!
!%! ! Need a temporary data structure to build this.
!%! ! These arrays keep a list of the x and y indices of dipole quads
!%! ! that potentially overlap the regular boxes. Need one for the u
!%! ! and one for the t grid.
!%! integer, allocatable :: ureg_list_lon(:,:,:)
!%! integer, allocatable :: ureg_list_lat(:,:,:)
!%! integer, allocatable :: treg_list_lon(:,:,:)
!%! integer, allocatable :: treg_list_lat(:,:,:)
!%!
!%! real(r8) :: u_c_lons(4), u_c_lats(4), t_c_lons(4), t_c_lats(4), pole_row_lon
!%! integer  :: i, j, k, pindex
!%! integer  :: reg_lon_ind(2), reg_lat_ind(2), u_total, t_total, u_index, t_index
!%! logical  :: is_pole
!%! integer  :: surf_index
!%!
!%! allocate(ureg_list_lon(num_reg_x, num_reg_y, max_reg_list_num))
!%! allocate(ureg_list_lat(num_reg_x, num_reg_y, max_reg_list_num))
!%! allocate(treg_list_lon(num_reg_x, num_reg_y, max_reg_list_num))
!%! allocate(treg_list_lat(num_reg_x, num_reg_y, max_reg_list_num))
!%!
!%! ! this is the level threshold for deciding whether we are over land
!%! ! or water.  to be valid all 4 corners of the quad must have a level
!%! ! number greater than this index.  (so 0 excludes all land points.)
!%! ! if you wanted to assimilate only in regions where the water depth is
!%! ! deeper than some threshold, set this index to N and only quads where
!%! ! all the level numbers are N+1 or deeper will be used.
!%! surf_index = 1
!%!
!%! ! Begin by finding the quad that contains the pole for the dipole t_grid.
!%! ! To do this locate the u quad with the pole on its right boundary. This is on
!%! ! the row that is opposite the shifted pole and exactly follows a lon circle.
!%! pole_x = nx / 2;
!%! ! Search for the row at which the longitude flips over the pole
!%! pole_row_lon = ulon(pole_x, 1);
!%! do i = 1, ny
!%!    pindex = i
!%!    if(ulon(pole_x, i) /= pole_row_lon) exit
!%! enddo
!%!
!%! ! Pole boxes for u have indices pole_x or pole_x-1 and index - 1;
!%! ! (it's right before the flip).
!%! u_pole_y = pindex - 1;
!%!
!%! ! Locate the T dipole quad that contains the pole.
!%! ! We know it is in either the same lat quad as the u pole edge or one higher.
!%! ! Figure out if the pole is more or less than halfway along
!%! ! the u quad edge to find the right one.
!%! if(ulat(pole_x, u_pole_y) > ulat(pole_x, u_pole_y + 1)) then
!%!    t_pole_y = u_pole_y;
!%! else
!%!    t_pole_y = u_pole_y + 1;
!%! endif
!%!
!%! ! Loop through each of the dipole grid quads
!%! do i = 1, nx
!%!    ! There's no wraparound in y, one box less than grid boundaries
!%!    do j = 1, ny - 1
!%!
!%!       if( all_corners_valid(i,j) ) then
!%!          ! Set up array of lons and lats for the corners of these u quads
!%!          call get_quad_corners(ulon, i, j, cyclic, u_c_lons)
!%!          call get_quad_corners(ulat, i, j, cyclic, u_c_lats)
!%!
!%!          ! Get list of regular boxes that cover this u dipole quad
!%!          ! false indicates that for the u grid there's nothing special about pole
!%!          call reg_box_overlap(u_c_lons, u_c_lats, .false., reg_lon_ind, reg_lat_ind)
!%!          ! Update the temporary data structures for the u quad
!%!          call update_reg_list(u_dipole_num, ureg_list_lon, &
!%!             ureg_list_lat, reg_lon_ind, reg_lat_ind, i, j)
!%!       endif
!%!
!%!       ! Repeat for t dipole quads.
!%!       ! Only update regular boxes that contain all valid corners
!%!       if( all_corners_valid(i,j) ) then
!%!          ! Set up array of lons and lats for the corners of these t quads
!%!          call get_quad_corners(tlon, i, j, cyclic, t_c_lons)
!%!          call get_quad_corners(tlat, i, j, cyclic, t_c_lats)
!%!
!%!          ! Is this the pole quad for the T grid?
!%!          is_pole = (i == pole_x .and. j == t_pole_y)
!%!
!%!          call reg_box_overlap(t_c_lons, t_c_lats, is_pole, reg_lon_ind, reg_lat_ind)
!%!          call update_reg_list(t_dipole_num, treg_list_lon, &
!%!             treg_list_lat, reg_lon_ind, reg_lat_ind, i, j)
!%!       endif
!%!    enddo
!%! enddo
!%!
!%! if (do_output()) write(*,*)'to determine (minimum) max_reg_list_num values for new grids ...'
!%! if (do_output()) write(*,*)'u_dipole_num is ',maxval(u_dipole_num)
!%! if (do_output()) write(*,*)'t_dipole_num is ',maxval(t_dipole_num)
!%!
!%! ! Invert the temporary data structure. The total number of entries will be
!%! ! the sum of the number of dipole cells for each regular cell.
!%! u_total = sum(u_dipole_num)
!%! t_total = sum(t_dipole_num)
!%!
!%! ! Allocate storage for the final structures in module storage
!%! allocate(u_dipole_lon_list(u_total), u_dipole_lat_list(u_total))
!%! allocate(t_dipole_lon_list(t_total), t_dipole_lat_list(t_total))
!%!
!%! ! Fill up the long list by traversing the temporary structure. Need indices
!%! ! to keep track of where to put the next entry.
!%! u_index = 1
!%! t_index = 1
!%!
!%! ! Loop through each regular grid box
!%! do i = 1, num_reg_x
!%!    do j = 1, num_reg_y
!%!
!%!       ! The list for this regular box starts at the current indices.
!%!       u_dipole_start(i, j) = u_index
!%!       t_dipole_start(i, j) = t_index
!%!
!%!       ! Copy all the close dipole quads for regular u box(i, j)
!%!       do k = 1, u_dipole_num(i, j)
!%!          u_dipole_lon_list(u_index) = ureg_list_lon(i, j, k)
!%!          u_dipole_lat_list(u_index) = ureg_list_lat(i, j, k)
!%!          u_index = u_index + 1
!%!       enddo
!%!
!%!       ! Copy all the close dipoles for regular t box (i, j)
!%!       do k = 1, t_dipole_num(i, j)
!%!          t_dipole_lon_list(t_index) = treg_list_lon(i, j, k)
!%!          t_dipole_lat_list(t_index) = treg_list_lat(i, j, k)
!%!          t_index = t_index + 1
!%!       enddo
!%!
!%!    enddo
!%! enddo
!%!
!%! ! Confirm that the indices come out okay as debug
!%! if(u_index /= u_total + 1) then
!%!    string1 = 'Storage indices did not balance for U grid: : contact DART developers'
!%!    call error_handler(E_ERR, 'init_dipole_interp', string1, source, revision, revdate)
!%! endif
!%! if(t_index /= t_total + 1) then
!%!    string1 = 'Storage indices did not balance for T grid: : contact DART developers'
!%!    call error_handler(E_ERR, 'init_dipole_interp', string1, source, revision, revdate)
!%! endif
!%!
!%! end subroutine init_dipole_interp

!------------------------------------------------------------
!>@todo FIXME if i'm doing this right, we shouldn't have
!> to have this routine or the next one.

subroutine get_quad_grid_size(interp_handle, nlon, nlat)

type(quad_interp_handle), intent(in)  :: interp_handle
integer,                  intent(out) :: nlon, nlat

nlon = interp_handle%nlon
nlat = interp_handle%nlat

end subroutine get_quad_grid_size

!------------------------------------------------------------

function get_quad_global(interp_handle)

type(quad_interp_handle), intent(in)  :: interp_handle
logical                               :: get_quad_global

get_quad_global = interp_handle%opt%global_grid

end function get_quad_global

!------------------------------------------------------------
!> Given a longitude and latitude in degrees returns the index of the regular
!> lon-lat box that contains the point.  if this a global grid it cannot fail.
!> if this is a regional grid, the given (lon,lat) might be outside of the region.
!> if istatus=0, good return.  istatus=1 bad box numbers.

subroutine get_reg_box_indices(h, lon, lat, x_ind, y_ind, istatus)

type(quad_interp_handle), intent(in)  :: h
real(r8),                 intent(in)  :: lon, lat
integer,                  intent(out) :: x_ind, y_ind
integer,                  intent(out) :: istatus

istatus = 0

call get_reg_lon_box(h, lon, x_ind)
call get_reg_lat_box(h, lat, y_ind)

if ( (.not. get_quad_global(h)) .and. &
     (x_ind < 1 .or. x_ind > h%ii%num_reg_x .or. &
      y_ind < 1 .or. y_ind > h%ii%num_reg_y)) then
      istatus = 1
endif

end subroutine get_reg_box_indices

!------------------------------------------------------------
!> Determine which regular longitude box contains the longitude of interest

subroutine get_reg_lon_box(h, lon, x_ind)

type(quad_interp_handle), intent(in)  :: h
real(r8),                 intent(in)  :: lon
integer,                  intent(out) :: x_ind

x_ind = int(h%ii%num_reg_x * (lon - h%ii%min_lon) / h%ii%lon_width) + 1
!print *, 'get_reg_lon_box: ', h%ii%num_reg_x, lon, h%ii%min_lon, h%ii%lon_width, ((lon - h%ii%min_lon) / h%ii%lon_width), x_ind

! Watch out for exactly at top; assume all lats and lons in legal range
if(lon == h%ii%max_lon) x_ind = h%ii%num_reg_x

end subroutine get_reg_lon_box

!------------------------------------------------------------
!> Determine which regular latitude box contains the latitude of interest

subroutine get_reg_lat_box(h, lat, y_ind)

type(quad_interp_handle), intent(in)  :: h
real(r8),                 intent(in)  :: lat
integer,                  intent(out) :: y_ind

y_ind = int(h%ii%num_reg_y * (lat - h%ii%min_lat) / h%ii%lat_width) + 1
!print *, 'get_reg_lat_box: ', h%ii%num_reg_y, lat, h%ii%min_lat, h%ii%lat_width, ((lat - h%ii%min_lat) / h%ii%lat_width), y_ind

! Watch out for exactly at top; assume all lats and lons in legal range
if(lat == h%ii%max_lat)  y_ind = h%ii%num_reg_y

end subroutine get_reg_lat_box

!------------------------------------------------------------
!> Find a set of regular lat lon boxes that covers all of the area covered by
!> a distorted grid quad whose corners are given by the dimension four x_corners
!> and y_corners arrays.

subroutine reg_box_overlap(h, x_corners, y_corners, is_pole, &
                           reg_lon_ind, reg_lat_ind)

type(quad_interp_handle), intent(in) :: h
real(r8),                 intent(in)  :: x_corners(4), y_corners(4)
logical,                  intent(in)  :: is_pole
integer,                  intent(out) :: reg_lon_ind(2), reg_lat_ind(2)

! The two dimensional arrays reg_lon_ind and reg_lat_ind
! return the first and last indices of the regular boxes in latitude and
! longitude respectively. These indices may wraparound for reg_lon_ind.
! A special computation is needed for a dipole quad that has the true north
! pole in its interior. The logical is_pole is set to true if this is the case.
! This can only happen for the t grid.  If the longitude boxes overlap 0
! degrees, the indices returned are adjusted by adding the total number of
! boxes to the second index (e.g. the indices might be 88 and 93 for a case
! with 90 longitude boxes).

real(r8) :: lat_min, lat_max, lon_min, lon_max
integer  :: i, nrx, nry

nrx = h%ii%num_reg_x
nry = h%ii%num_reg_y

!  A quad containing the pole is fundamentally different
if(is_pole) then
   ! Need all longitude boxes
   reg_lon_ind(1) = 1
   reg_lon_ind(2) = nrx
   ! Need to cover from lowest latitude to top box
   lat_min = minval(y_corners)
   reg_lat_ind(1) = int(nry * (lat_min + 90.0_r8) / 180.0_r8) + 1
   call get_reg_lat_box(h, lat_min, reg_lat_ind(1))
   reg_lat_ind(2) = nry
else
   ! All other quads do not contain pole (pole could be on edge but no problem)
   ! This is specific to the dipole POP grids that do not go to the south pole
   ! Finding the range of latitudes is cake
   lat_min = minval(y_corners)
   lat_max = maxval(y_corners)

   ! Figure out the indices of the regular boxes for min and max lats
   call get_reg_lat_box(h, lat_min, reg_lat_ind(1))
   call get_reg_lat_box(h, lat_max, reg_lat_ind(2))

   ! Lons are much trickier. Need to make sure to wraparound the
   ! right way. There is no guarantee on direction of lons in the
   ! high latitude dipole rows.
   ! All longitudes for non-pole rows have to be within 180 degrees
   ! of one another.
   lon_min = minval(x_corners)
   lon_max = maxval(x_corners)
   if((lon_max - lon_min) > 180.0_r8) then
      ! If the max longitude value is more than 180
      ! degrees larger than the min, then there must be wraparound.
      ! Then, find the smallest value > 180 and the largest < 180 to get range.
      lon_min = 360.0_r8
      lon_max =   0.0_r8
      do i=1, 4
         if(x_corners(i) > 180.0_r8 .and. x_corners(i) < lon_min) lon_min = x_corners(i)
         if(x_corners(i) < 180.0_r8 .and. x_corners(i) > lon_max) lon_max = x_corners(i)
      enddo
   endif

   ! Get the indices for the extreme longitudes
   call get_reg_lon_box(h, lon_min, reg_lon_ind(1))
   call get_reg_lon_box(h, lon_max, reg_lon_ind(2))

   ! Watch for wraparound again; make sure that second index is greater than first
   if(reg_lon_ind(2) < reg_lon_ind(1)) reg_lon_ind(2) = reg_lon_ind(2) + nrx
endif

end subroutine reg_box_overlap

!------------------------------------------------------------
!> Grabs the corners for a given quadrilateral from the global array of lower
!> right corners. Note that corners go counterclockwise around the quad.

subroutine get_quad_corners(x, i, j, cyclic, pole, nx, ny, corners, istatus)

real(r8), intent(in)  :: x(:, :)
integer,  intent(in)  :: i, j
logical,  intent(in)  :: cyclic, pole
integer,  intent(in)  :: nx, ny
real(r8), intent(out) :: corners(4)
integer,  intent(out) :: istatus

integer :: ip1, jp1

! for global grids have to worry about wrapping.
! for regional grids you might be at the grid edge.
istatus = 1
corners(:) = MISSING_R8

! find the other indices for this quad, respecting cyclic and pole settings.
call quad_index_neighbors(i, j, nx, ny, cyclic, pole, ip1, jp1)
if (ip1 < 0 .or. jp1 < 0) then
   istatus = -1
   return
endif

corners(1) = x(i,   j  )
corners(2) = x(ip1, j  )
corners(3) = x(ip1, jp1)
corners(4) = x(i,   jp1)

istatus = 0

end subroutine get_quad_corners

!-----------------------------------------------------------------------
!> given lon/lat indices, add one to lon and lat.
!> check for wraparound in lon, and pole for lat.

subroutine quad_index_neighbors(lon_index, lat_index, nx, ny, cyclic, pole, &
                                next_lon, next_lat)
integer, intent(in)  :: lon_index
integer, intent(in)  :: lat_index
integer, intent(in)  :: nx, ny
logical, intent(in)  :: cyclic, pole
integer, intent(out) :: next_lon
integer, intent(out) :: next_lat

! if cyclic, wrap back to 1
next_lon = lon_index+1
if (next_lon > nx) then
   if (cyclic) then
      next_lon = 1
   else
      ! FIXME: is this expected for a regional grid?
      ! or have we already bounds checked the edges?
      write(string1, *) ' next_lon > nx and not cyclic', &
                          lon_index, lat_index, cyclic, nx, next_lon
      call error_handler(E_MSG, 'quad_index_neighbors', string1)
      next_lon = -1
   endif
endif

! if over poles, subtract one to continue down the other side
next_lat = lat_index+1
if (next_lat > ny) then
   if (pole) then
      next_lat = ny - 1
   else
      ! FIXME: is this expected for a regional grid?
      ! or have we already bounds checked the edges?
      write(string1, *) ' next_lat > ny and not polar', &
                          lon_index, lat_index, pole, ny, next_lat
      call error_handler(E_MSG, 'quad_index_neighbors', string1)
      next_lat = -1
   endif
endif

end subroutine quad_index_neighbors


!------------------------------------------------------------
!> Updates the data structure listing dipole quads that are in a given regular box

subroutine update_reg_list(reg_list_num, reg_list_lon, reg_list_lat, &
                           reg_lon_ind, reg_lat_ind, nrx, nry, maxlist, &
                           grid_lon_index, grid_lat_index)

integer, intent(inout) :: reg_list_num(:, :)
integer, intent(inout) :: reg_list_lon(:, :, :)
integer, intent(inout) :: reg_list_lat(:, :, :)
integer, intent(inout) :: reg_lon_ind(2)
integer, intent(inout) :: reg_lat_ind(2)
integer, intent(in)    :: nrx, nry, maxlist
integer, intent(in)    :: grid_lon_index, grid_lat_index

integer :: ind_x, index_x, ind_y, index_y

!print *, 'update_reg_list called for ', grid_lon_index, grid_lat_index
!print *, 'update_reg_list bins: ', reg_lon_ind(1), reg_lon_ind(2), reg_lat_ind(1), reg_lat_ind(2)

! Loop through indices for each possible regular cell
! Have to watch for wraparound in longitude
if(reg_lon_ind(2) < reg_lon_ind(1)) reg_lon_ind(2) = reg_lon_ind(2) + nrx

do ind_x = reg_lon_ind(1), reg_lon_ind(2)
   ! Inside loop, need to go back to wraparound indices to find right box
   index_x = ind_x
   if(index_x > nrx) index_x = index_x - nrx

   do ind_y = reg_lat_ind(1), reg_lat_ind(2)
      index_y = ind_y
      if(index_y > nry) index_y = index_y - nry

      if ((index_x < 1 .or. index_x > nrx) .or. (index_y < 1 .or. index_y > nry)) then
         string1 = 'unable to find right box'
         write(string2,*) 'index_x may be out-of-range: ', 1, index_x, nrx
         write(string3,*) 'index_y may be out-of-range: ', 1, index_y, nry
         call error_handler(E_ERR,'update_reg_list',string1, &
                 source, revision, revdate, text2=string2, text3=string3)
      endif

      ! Make sure the list storage isn't full
!print *, 'reg_list_num, x, y = ', reg_list_num, index_x, index_y
      if(reg_list_num(index_x, index_y) >= maxlist) then
         write(string1,*) 'max_reg_list_num (',maxlist,') is too small ... increase'
         write(string2,*) 'adding 1 to bin ', index_x, index_y
         write(string3,*) 'bins: ', reg_lon_ind(1), reg_lon_ind(2), &
                                    reg_lat_ind(1), reg_lat_ind(2)
         call error_handler(E_ERR, 'update_reg_list', string1, &
                            source, revision, revdate, text2=string2, text3=string3)
      endif

      ! Increment the count
      reg_list_num(index_x, index_y) = reg_list_num(index_x, index_y) + 1
      ! Store this quad in the list for this regular box
      reg_list_lon(index_x, index_y, reg_list_num(index_x, index_y)) = grid_lon_index
      reg_list_lat(index_x, index_y, reg_list_num(index_x, index_y)) = grid_lat_index
      !print *, 'adding 1 to bin ', index_x, index_y, ' for ', grid_lon_index, grid_lat_index, &
      !          ' now entries = ', reg_list_num(index_x, index_y)
   enddo
enddo

end subroutine update_reg_list

!------------------------------------------------------------------
!> Subroutine to locate the given lon lat location and return the
!> longitude/latitude indices for the 4 corners.  This routine is
!> for irregular grids only.
!>
!> Two different types of grids are used here. The irregular grid
!> regions are referred to as quads (short for quadrilateral).
!> The lon/lat regular grid regions are called boxes.
!> All grids are referenced by the index of the lower left corner
!> of the quad or box.
!>

subroutine quad_lon_lat_locate_ii(interp_handle, lon, lat, &
                                  four_lon_indices, four_lat_indices, istatus)
type(quad_interp_handle), intent(in)  :: interp_handle
real(r8),                 intent(in)  :: lon, lat
integer,                  intent(out) :: four_lon_indices(4), four_lat_indices(4)
integer,                  intent(out) :: istatus

! Local storage
integer  :: num_inds, start_ind
integer  :: x_ind, y_ind, nx, ny
integer  :: lon_bot, lat_bot, lon_top, lat_top
logical  :: cyclic, pole
real(r8) :: x_corners(4), y_corners(4)

character(len=*), parameter :: routine = 'quad_lon_lat_locate:quad_lon_lat_locate_ii'

! Succesful return has istatus of 0
istatus = 0

! shorter var names for ease in reading the code
nx = interp_handle%nlon
ny = interp_handle%nlat
cyclic = interp_handle%opt%spans_lon_zero
pole   = interp_handle%opt%pole_wrap

select case (interp_handle%grid_type)

 case (GRID_QUAD_FULLY_IRREGULAR)
   ! Figure out which of the regular grid boxes this is in
   call get_reg_box_indices(interp_handle, lon, lat, x_ind, y_ind, istatus)
   if (istatus /= 0) return

   num_inds =  interp_handle%ii%grid_num  (x_ind, y_ind)
   start_ind = interp_handle%ii%grid_start(x_ind, y_ind)

   ! If there are no quads overlapping, can't do interpolation
   if(num_inds == 0) then
      istatus = 1
      return
   endif

   ! Search the list of quads to see if (lon, lat) is in one
   call get_grid_quad(lon, lat, interp_handle%ii%lons_2d, interp_handle%ii%lats_2d, &
                      num_inds, start_ind, interp_handle%ii%grid_lon_list, &
                      interp_handle%ii%grid_lat_list, cyclic, pole, nx, ny, &
                      lon_bot, lat_bot, istatus)
   if (debug > 10) print *, 'get_grid_quad returns lon/lat bot: ', lon_bot, lat_bot
   if (istatus /= 0) return

   ! Getting corners for accurate interpolation
   call get_quad_corners(interp_handle%ii%lons_2d, lon_bot, lat_bot, cyclic, pole, &
                         nx, ny, x_corners, istatus)
   if (debug > 10) print *, 'get_quad_corners returns x_corners: ', x_corners
   if (istatus /= 0) return
   call get_quad_corners(interp_handle%ii%lats_2d, lon_bot, lat_bot, cyclic, pole, &
                         nx, ny, y_corners, istatus)
   if (debug > 10) print *, 'get_quad_corners returns y_corners: ', y_corners
   if (istatus /= 0) return

   ! this test shouldn't be needed
   if ( .not. all_corners_valid(interp_handle%opt, lon_bot, lat_bot, nx)) then
      string1 = 'got into a quad where at least one of the corners is not valid. should not happen'
      write(string2,*) 'lon/lat bot, nx, lon/lat', lon_bot, lat_bot, nx, lon, lat
      call error_handler(E_ERR, routine, string1, &
                         source, revision, revdate, text2=string2)
   endif

   ! Fail if point is in one of the U boxes that go through the
   ! pole (this could be fixed up if necessary)
   if (lat_bot == u_pole_y .and. &
      (lon_bot == pole_x -1 .or. lon_bot == pole_x)) then
      istatus = 4
      return
   endif

 case (GRID_QUAD_IRREG_SPACED_REGULAR)
 case (GRID_QUAD_FULLY_REGULAR)
   string1 = 'this version of the call only work on fully irregular grids'
   call error_handler(E_ERR, routine, string1, source, revision, revdate)

 case default
   call error_handler(E_ERR, routine, 'unrecognized grid type', &
                      source, revision, revdate)
end select

call quad_index_neighbors(lon_bot, lat_bot, nx, ny, cyclic, pole, &
                          lon_top, lat_top)

if (lon_top < 0 .or. lat_top < 0) then
   istatus = 2
   return
endif

!! Find the indices to get the values for interpolating
!lat_top = lat_bot + 1
!if(lat_top > ny) then
!   istatus = 2
!   return
!endif
!
!! Watch for wraparound in longitude
!lon_top = lon_bot + 1
!if(lon_top > nx) then
!   if (cyclic) then
!      lon_top = 1
!   else
!      istatus = 2
!      return
!   endif
!endif

! the 4 return values set here are:  lon_bot, lat_bot, lon_top, lat_top
four_lon_indices(1) = lon_bot
four_lon_indices(2) = lon_top
four_lon_indices(3) = lon_top
four_lon_indices(4) = lon_bot

four_lat_indices(1) = lat_bot
four_lat_indices(2) = lat_bot
four_lat_indices(3) = lat_top
four_lat_indices(4) = lat_top

end subroutine quad_lon_lat_locate_ii

!------------------------------------------------------------------
!> Subroutine to interpolate to a lon lat location given the state vector
!> for that level, x. This works just on one horizontal slice.
!> This routine works for either the dipole or a regular lat-lon grid.
!> Successful interpolation returns istatus=0.

!>@todo FIXME should this still return lon_fract, lat_fract?
!>(thinking yes)
subroutine quad_lon_lat_locate_ir(interp_handle, lon, lat, &
                four_lons, four_lats, lon_fract, lat_fract, istatus)

type(quad_interp_handle), intent(in)  :: interp_handle
real(r8),                 intent(in)  :: lon, lat
integer,                  intent(out) :: four_lons(4), four_lats(4)
real(r8),                 intent(out) :: lon_fract, lat_fract
integer,                  intent(out) :: istatus

! Local storage
integer  :: nx, ny
logical  :: cyclic
integer  :: lon_bot, lat_bot, lon_top, lat_top

character(len=*), parameter :: routine = 'quad_lon_lat_locate:quad_lon_lat_locate_ir'

! Succesful return has istatus of 0
istatus = 0

! shorter var names for ease in reading the code
nx = interp_handle%nlon
ny = interp_handle%nlat
cyclic = interp_handle%opt%spans_lon_zero

select case (interp_handle%grid_type)

 case (GRID_QUAD_FULLY_IRREGULAR)
   string1 = 'this version of the call only work on partially or fully regular grids'
   call error_handler(E_ERR, routine, string1, source, revision, revdate)

 case (GRID_QUAD_IRREG_SPACED_REGULAR)
   ! This is an irregular grid (irregular == spacing; still completely orthogonal)
   call get_semireg_box(lon, lat, nx, ny, &
         interp_handle%ir%lons_1d, interp_handle%ir%lats_1d, &
         interp_handle%opt%north_to_south, &
         lon_bot, lat_bot, lon_fract, lat_fract, istatus)

 case (GRID_QUAD_FULLY_REGULAR)
   ! evenly spaced and orthogonal
   call get_reg_box(lon, lat, nx, ny, &
         interp_handle%rr%lon_start, interp_handle%rr%lat_start, &
         interp_handle%rr%lon_delta, interp_handle%rr%lat_delta, &
         lon_bot, lat_bot, lon_fract, lat_fract, istatus)

 case default
   call error_handler(E_ERR, routine, 'unrecognized grid type', &
                      source, revision, revdate)
end select

if (istatus /= 0) return

! Find the indices to get the values for interpolating
lat_top = lat_bot + 1
if(lat_top > ny) then
   istatus = 2
   return
endif

! Watch for wraparound in longitude
lon_top = lon_bot + 1
if(lon_top > nx) then
   if (cyclic) then
      lon_top = 1
   else
      istatus = 2
      return
   endif
endif

! the 6 values set so far in this routine are:
!    lon_bot, lat_bot, lon_top, lat_top, lon_fract, lat_fract
!
! now fill arrays so they are easy to process in the calling
! code in the right order, which is counterclockwise around the quad:
!
!  (lon_bot, lat_bot), (lon_top, lat_bot), (lon_top, lat_top), (lon_bot, lat_top)

four_lons(1) = lon_bot
four_lons(2) = lon_top
four_lons(3) = lon_top
four_lons(4) = lon_bot

four_lats(1) = lat_bot
four_lats(2) = lat_bot
four_lats(3) = lat_top
four_lats(4) = lat_top

end subroutine quad_lon_lat_locate_ir

!------------------------------------------------------------
!> Given a longitude and latitude of a point (lon and lat) and the
!> longitudes and latitudes of the lower left corner of the regular grid
!> boxes, gets the indices of the grid box that contains the point and
!> the fractions along each direction for interpolation.

subroutine get_irreg_box(lon, lat, nx, ny, lon_array, lat_array, cyclic, &
                         found_x, found_y, lon_fract, lat_fract, istatus)

real(r8),  intent(in)  :: lon, lat
integer,   intent(in)  :: nx, ny
real(r8),  intent(in)  :: lon_array(nx, ny), lat_array(nx, ny)
logical,   intent(in)  :: cyclic
real(r8),  intent(out) :: lon_fract, lat_fract
integer,   intent(out) :: found_x, found_y, istatus

! Local storage
integer  :: lat_status, lon_top, lat_top

! Succesful return has istatus of 0
istatus = 0

! Get latitude box boundaries and fraction
! lat_bounds() is called for both the fully irregular case and for
! the partially regular case.  for irreg, we don't need the options
! for covering the poles or inverted latitude arrays, so hardcode
! those to false.
call lat_bounds(lat, ny, lat_array(1,:), .false., .false., &
                found_y, lat_top, lat_fract, lat_status)

! Check for error on the latitude interpolation
if(lat_status /= 0) then
   istatus = 1
   return
endif

call lon_bounds(lon, nx, lon_array(1,:), cyclic, &
                found_x, lon_top, lon_fract, istatus)

end subroutine get_irreg_box

!------------------------------------------------------------
!> Given a longitude and latitude of a point (lon and lat) and the
!> start and deltas of the lons and lats, get the lower left indicies
!> of the grid box that contains the point and
!> the fractions along each direction for interpolation.

subroutine get_reg_box(lon, lat, nx, ny, lon_min, lat_min, lon_del, lat_del, &
                       found_x, found_y, lon_fract, lat_fract, istatus)

real(r8),   intent(in) :: lon, lat
integer,    intent(in) :: nx, ny
real(r8),   intent(in) :: lon_min, lat_min, lon_del, lat_del
integer,   intent(out) :: found_x, found_y
real(r8),  intent(out) :: lon_fract, lat_fract
integer,   intent(out) :: istatus

! Local storage
integer  :: lat_status, lon_top, lat_top, i
real(r8) :: lon_array(nx), lat_array(ny)

! Succesful return has istatus of 0
istatus = 0

!>@todo  FIXME: hack to get code running.  don't expand arrays - slow.
! search directly in a loop w/ deltas.
do i=1, nx
   lon_array(i) = lon_min + (i-1)*lon_del
enddo
do i=1, ny
   lat_array(i) = lat_min + (i-1)*lat_del
enddo

! Get latitude box boundaries
call lat_bounds(lat, ny, lat_array, .false., .false., &
                found_y, lat_top, lat_fract, lat_status)

! Check for error on the latitude interpolation
if(lat_status /= 0) then
   istatus = 1
   return
endif

! Find out what longitude box and fraction - FIXME: cyclic flag
call lon_bounds(lon, nx, lon_array, .true., found_x, lon_top, lon_fract, istatus)

end subroutine get_reg_box


!------------------------------------------------------------
!> Given a longitude and latitude array for irregular spaced,
!> orthogonal grids, get the lower left indices of the grid box
!> that contains the point and the fractions along each direction
!> for interpolation.

subroutine get_semireg_box(lon, lat, nx, ny, lon_array, lat_array, invert_lat, &
                           found_x, found_y, lon_fract, lat_fract, istatus)

real(r8),  intent(in)  :: lon, lat
integer,   intent(in)  :: nx, ny
real(r8),  intent(in)  :: lon_array(:), lat_array(:)
logical,   intent(in)  :: invert_lat
integer,   intent(out) :: found_x, found_y
real(r8),  intent(out) :: lon_fract, lat_fract
integer,   intent(out) :: istatus

! Local storage
integer  :: lat_status, lon_top, lat_top, i

! Succesful return has istatus of 0
istatus = 0

! Get latitude box boundaries
!>@todo FIXME check on the pole wrap and cyclic flags
call lat_bounds(lat, ny, lat_array, .false., invert_lat, &
                found_y, lat_top, lat_fract, lat_status)

! Check for error on the latitude interpolation
if(lat_status /= 0) then
   istatus = 1
   return
endif

! Find out what longitude box and fraction - FIXME: cyclic flag
call lon_bounds(lon, nx, lon_array, .true., found_x, lon_top, lon_fract, istatus)

end subroutine get_semireg_box

!------------------------------------------------------------
!> assumes longitudes can be described by a single 1D array.
!> Given a longitude lon, the array of longitudes for grid boundaries, and the
!> number of longitudes in the grid, returns the indices of the longitude
!> below and above the location longitude and the fraction of the distance
!> between. if 'cyclic=.true' the longitude wraps around for a global grid.
!> Algorithm fails for a silly grid that has only two longitudes separated by 180 degrees.

subroutine lon_bounds(lon, nlons, lon_array, cyclic, bot, top, fract, istatus)

real(r8),   intent(in)  :: lon
integer,    intent(in)  :: nlons
real(r8),   intent(in)  :: lon_array(nlons)
logical,    intent(in)  :: cyclic
integer,    intent(out) :: bot, top
real(r8),   intent(out) :: fract
integer,    intent(out) :: istatus

! Local storage
integer  :: i
real(r8) :: dist_bot, dist_top
logical  :: span ! FIXME: unneeded?

! Success should return 0, failure a positive number.
istatus = 0

! If not cyclic, check for too far west or east
! span is true if the longitude array crosses the prime meridian
if (.not. cyclic) then
   if(lon < lon_array(1)) then
      istatus = 1
      return
   else if(lon > lon_array(nlons)) then
      istatus = 2
      return
   endif
   span = (lon_array(1) > lon_array(nlons))
else
   span = .true.
endif


! search through middle
do i = 2, nlons
   dist_bot = lon_dist(lon, lon_array(i - 1))
   dist_top = lon_dist(lon, lon_array(i))
   if (debug > 3) print *, 'lon: i, lon_array(i-1), lon_array(i), bot, top: ', &
                            i, lon_array(i-1), lon_array(i), dist_bot, dist_top
   if(dist_bot <= 0.0_r8 .and. dist_top > 0.0_r8) then
      bot = i - 1
      top = i
      if ((abs(dist_bot) + dist_top) == 0.0_r8) then
         istatus = 2
         return
      endif
      fract = abs(dist_bot) / (abs(dist_bot) + dist_top)
      if (debug > 3) print *, 'lon: returning bot, top, fract', bot, top, fract
      return
   endif
enddo


! Falling off the end means it's in between; wraparound
if (cyclic) then
   bot = nlons
   top = 1
   dist_bot = lon_dist(lon, lon_array(bot))
   dist_top = lon_dist(lon, lon_array(top))
   if ((abs(dist_bot) + dist_top) == 0.0_r8) then
      istatus = 2
      return
   endif
   fract = abs(dist_bot) / (abs(dist_bot) + dist_top)
else
   string1 = 'end reached. internal error, should not happen'
   write(string2,*)'lon of interest is ',lon
   call error_handler(E_ERR, 'lon_bounds', string1, &
                      source, revision, revdate, text2=string2)
endif


end subroutine lon_bounds

!-------------------------------------------------------------
!> assumes latitudes can be described by a single 1D array.
!> Given a latitude lat, the array of latitudes for grid boundaries, and the
!> number of latitudes in the grid, returns the indices of the latitude
!> below and above the location latitude and the fraction of the distance
!> between. istatus is returned as 0 unless the location latitude is
!> south of the southernmost grid point (1 returned) or north of the
!> northernmost (2 returned). If one really had lots of polar obs would
!> want to worry about interpolating around poles.

subroutine lat_bounds(lat, nlats, lat_array, polar, invert_lat, bot, top, fract, istatus)

real(r8),  intent(in)  :: lat
integer,   intent(in)  :: nlats
real(r8),  intent(in)  :: lat_array(nlats)
logical,   intent(in)  :: polar
logical,   intent(in)  :: invert_lat
integer,   intent(out) :: bot, top
real(r8),  intent(out) :: fract
integer,   intent(out) :: istatus

! Local storage
integer :: i, north, south, start, end, increment

! Success should return 0, failure a positive number.
istatus = 0

! FIXME: polar is for future expansion, ignored for now

! normally grids start at -90 and go to 90 (south to north).
! but some grids start at 90 and go to -90 (north to south).
! try to handle both with a minimum of replicated code.
if (invert_lat) then
  north = lat_array(1)
  south = lat_array(nlats)
  start = nlats - 1
  end = 1
  increment = -1
else
  north = lat_array(nlats)
  south = lat_array(1)
  start = 2
  end = nlats
  increment = 1
endif

! Check for too far south or north
if(lat < south) then
   istatus = 1
   return
else if(lat > north) then
   istatus = 2
   return
endif

! In the middle, search through
do i = start, end, increment
   if (debug > 3) print *, 'lat: i, lat, lat(i): ', i, lat, lat_array(i)
   if(lat <= lat_array(i)) then
      if (invert_lat) then
         bot = i
         top = i + 1
      else
         bot = i - 1
         top = i
      endif
      if (lat_array(top) - lat_array(bot) == 0.0_r8) then
         istatus = 2
         return
      endif
      fract = (lat - lat_array(bot)) / (lat_array(top) - lat_array(bot))
      if (debug > 3) print *, 'lat: returning bot, top, fract', bot, top, fract
      return
   endif
enddo

string1 = 'end reached. internal error, should not happen'
write(string2,*)'lat of interest is ',lat
call error_handler(E_ERR, 'lat_bounds', string1, &
                   source, revision, revdate, text2=string2)

end subroutine lat_bounds

!------------------------------------------------------------------
!> Returns the smallest signed distance between lon1 and lon2 on the sphere

function lon_dist(lon1, lon2)

real(r8), intent(in) :: lon1, lon2
real(r8)             :: lon_dist

! If lon1 is less than 180 degrees east of lon2 the distance is negative
! If lon1 is less than 180 degrees west of lon2 the distance is positive

lon_dist = lon2 - lon1
if(lon_dist >= -180.0_r8 .and. lon_dist <= 180.0_r8) then
   return
else if(lon_dist< -180.0_r8) then
   lon_dist = lon_dist + 360.0_r8
else
   lon_dist = lon_dist - 360.0_r8
endif

end function lon_dist

!------------------------------------------------------------
!> Given the lon and lat of a point, and a list of the
!> indices of the quads that might contain a point at (lon, lat), determines
!> which quad contains the point.  istatus is returned as 0 if all went
!> well and 1 if the point was not found to be in any of the quads.

subroutine get_grid_quad(lon, lat, qlons, qlats, num_inds, start_ind, &
                         x_inds, y_inds, cyclic, pole, nx, ny, found_x, found_y, istatus)

real(r8), intent(in)  :: lon, lat, qlons(:, :), qlats(:, :)
integer,  intent(in)  :: num_inds, start_ind, x_inds(:), y_inds(:)
logical,  intent(in)  :: cyclic, pole
integer,  intent(in)  :: nx, ny
integer,  intent(out) :: found_x, found_y, istatus


integer :: i, my_index
real(r8) :: x_corners(4), y_corners(4)

! Loop through all the quads and see if the point is inside
do i = 1, num_inds
   my_index = start_ind + i - 1
   call get_quad_corners(qlons, x_inds(my_index), y_inds(my_index), &
                         cyclic, pole, nx, ny, x_corners, istatus)
   if (istatus /= 0) return

   call get_quad_corners(qlats, x_inds(my_index), y_inds(my_index), &
                         cyclic, pole, nx, ny, y_corners, istatus)
   if (istatus /= 0) return

   ! Ssearch in this individual quad
   if(in_quad(lon, lat, x_corners, y_corners)) then
      found_x = x_inds(my_index)
      found_y = y_inds(my_index)
      return
   endif
enddo

! Falling off the end means search failed, return istatus 1
istatus = 1

end subroutine get_grid_quad

!------------------------------------------------------------
!> Return in_quad true if the point (lon, lat) is in the quad with
!> the given corners.

function in_quad(lon, lat, x_corners, y_corners)

real(r8), intent(in)  :: lon
real(r8), intent(in)  :: lat
real(r8), intent(in)  :: x_corners(4)
real(r8), intent(in)  :: y_corners(4)
logical               :: in_quad

! Do this by line tracing in latitude for now. For non-pole point, want a vertical
! line from the lon, lat point to intersect a side of the quad both above
! and below the point.

real(r8) :: x(2), y(2)
logical  :: cant_be_in_box, in_box
integer  :: intercepts_above(4), intercepts_below(4), i
integer  :: num_above, num_below

! Default answer is point is not in quad
in_quad = .false.

! Loop through the sides and compute intercept (if any) with a vertical line
! from the point. This line has equation x=lon.
do i = 1, 4
   ! Load up the sides endpoints
   if(i <= 3) then
      x(1:2) = x_corners(i:i+1)
      y(1:2) = y_corners(i:i+1)
   else
      x(1) = x_corners(4)
      x(2) = x_corners(1)
      y(1) = y_corners(4)
      y(2) = y_corners(1)
   endif

   ! Check to see how a vertical line from the point is related to this side
   call line_intercept(x, y, lon, lat, cant_be_in_box, in_box, &
                       intercepts_above(i), intercepts_below(i))

   ! If cant_be_in_box is true, can return right away
   if(cant_be_in_box) then
      in_quad = .false.
      return
   ! Return true if on a side
   else if(in_box) then
      in_quad = .true.
      return
   endif

enddo

! See if the line intercepted a side of the quad both above and below
num_above = sum(intercepts_above)
num_below = sum(intercepts_below)

if(num_above > 0 .and. num_below > 0) then
   in_quad = .true.
endif

end function in_quad

!------------------------------------------------------------
!> Find the intercept of a vertical line from point (x_point, y_point) and
!> a line segment with endpoints side_x and side_y.

subroutine line_intercept(side_x_in, side_y, x_point_in, y_point, &
                          cant_be_in_box, in_box, intercept_above, intercept_below)

real(r8), intent(in)  :: side_x_in(2)
real(r8), intent(in)  :: side_y(2)
real(r8), intent(in)  :: x_point_in
real(r8), intent(in)  :: y_point
logical,  intent(out) :: cant_be_in_box
logical,  intent(out) :: in_box
integer,  intent(out) :: intercept_above
integer,  intent(out) :: intercept_below

! For a given side have endpoints (side_x1, side_y1) and (side_x2, side_y2)
! so equation of segment is y = side_y1 + m(x-side_x1) for y
! between side_y1 and side_y2.
! Intersection of vertical line and line containing side
! occurs at y = side_y1 + m(x_point - side_x1); need this
! y to be between side_y1 and side_y2.
! If the vertical line is colinear with the side but the point is not on the side, return
! cant_be_in_box as true. If the point is on the side, return in_box true.
! If the intersection of the vertical line and the side occurs at a point above
! the given point, return 1 for intercept_above. If the intersection occurs
! below, return 1 for intercept_below. If the vertical line does not intersect
! the segment, return false and 0 for all intent out arguments.

! WARNING: CERTAINLY A PROBLEM FOR THE POLE BOX!!! POLE BOX COULD
! HAVE SIDES THAT ARE LONGER THAN 180. For now pole boxes are excluded.

! This can probably be made much cleaner and more efficient.

real(r8) :: slope, y_intercept, side_x(2), x_point

! May have to adjust the longitude intent in values, so copy
side_x = side_x_in
x_point = x_point_in

! See if the side wraps around in longitude
if(maxval(side_x) - minval(side_x) > 180.0_r8) then
   if(side_x(1) < 180.0_r8)  side_x(1) = side_x(1) + 360.0_r8
   if(side_x(2) < 180.0_r8)  side_x(2) = side_x(2) + 360.0_r8
   if(x_point   < 180.0_r8)  x_point   = x_point   + 360.0_r8
endif

! Initialize the default returns
cant_be_in_box  = .false.
in_box          = .false.
intercept_above = 0
intercept_below = 0

! First easy check, if x_point is not between endpoints of segment doesn't intersect
if(x_point < minval(side_x) .or. x_point > maxval(side_x)) return

! Otherwise line must intersect the segment

! First subblock, slope is undefined
if(side_x(2) == side_x(1)) then
   ! The line is colinear with the side
   ! If y_point is between endpoints then point is on this side
   if(y_point <= maxval(side_y) .and. y_point >= minval(side_y)) then
      in_box = .true.
      return
   ! If not on side but colinear with side, point cant be in quad
   else
      cant_be_in_box = .true.
      return
   endif

else

   ! Second possibility; slope is defined
   ! FIXME: watch out for numerical instability.
   ! near-zero x's and large side_y's may cause overflow
   slope = (side_y(2) - side_y(1)) / (side_x(2) - side_x(1))

   ! Intercept of vertical line through is at x_point and...
   y_intercept = side_y(1) + slope * (x_point - side_x(1))

   ! Intersects the segment, is it above, below, or at the point
   if(y_intercept == y_point) then
      in_box = .true.
      return
   else if(y_intercept > y_point) then
      intercept_above = 1
      return
   else
      intercept_below = 1
      return
   endif
endif

end subroutine line_intercept

!------------------------------------------------------------
! Given a longitude and latitude (lon_in, lat_in), the longitude and
! latitude of the 4 corners of a quadrilateral and the values at the
! four corners, interpolates to (lon_in, lat) which is assumed to
! be in the quad. This is done by bilinear interpolation, fitting
! a function of the form a + bx + cy + dxy to the four points and
! then evaluating this function at (lon, lat). The fit is done by
! solving the 4x4 system of equations for a, b, c, and d. The system
! is reduced to a 3x3 by eliminating a from the first three equations
! and then solving the 3x3 before back substituting. There is concern
! about the numerical stability of this implementation. Implementation
! checks showed accuracy to seven decimal places on all tests.

subroutine quad_bilinear_interp(lon_in, lat_in, x_corners_in, y_corners_in, cyclic, &
                                p, expected_obs)

real(r8),  intent(in) :: lon_in, lat_in, x_corners_in(4), y_corners_in(4), p(4)
logical,   intent(in) :: cyclic
real(r8), intent(out) :: expected_obs

integer :: i
real(r8) :: m(3, 3), v(3), r(3), a, b(2), c(2), d
real(r8) :: x_corners(4), lon, y_corners(4), lat
real(r8) :: lon_mean, lat_mean, interp_val, angle

! Watch out for wraparound on x_corners.
lon = lon_in
x_corners = x_corners_in
lat = lat_in
y_corners = y_corners_in

if (debug > 10) write(*,'(A,4F12.3)') 'corner data values: ', p
if (debug > 10) write(*,'(A,4F12.3)') 'original x_corners: ', x_corners
if (debug > 10) write(*,'(A,4F12.3)') 'original y_corners: ', y_corners

!> @todo FIXME does this depend on cyclic or span flag???

! See if the side wraps around in longitude. If the corners longitudes
! wrap around 360, then the corners and the point to interpolate to
! must be adjusted to be in the range from 180 to 540 degrees.
if(maxval(x_corners) - minval(x_corners) > 180.0_r8) then
   if(lon < 180.0_r8) lon = lon + 360.0_r8
   do i = 1, 4
      if(x_corners(i) < 180.0_r8) x_corners(i) = x_corners(i) + 360.0_r8
   enddo
endif

!>@todo FIXME here is where can select and test various interpolation types

!*******
! Problems with extremes in polar cell interpolation can be reduced
! by this block, but it is not clear that it is needed for actual
! ocean grid data
!! Find the mean longitude of corners and remove
!lon_mean = sum(x_corners) / 4.0_r8
!lat_mean = sum(y_corners) / 4.0_r8
!
!x_corners = x_corners - lon_mean
!lon = lon - lon_mean
!! Multiply everybody by the cos of the latitude - why?
!do i = 1, 4
!   !x_corners(i) = x_corners(i) * cos(y_corners(i) * deg2rad)
!enddo
!!lon = lon * cos(lat * deg2rad)
!!lon_mean = lon_mean * cos(lat * deg2rad)
!
!y_corners = y_corners - lat_mean
!lat = lat - lat_mean

! try something else.  compute offsets from lower left,
! rotate so line segment 1-2 is horizontal, and then
! compute values.

if (do_rotate) then
   !print *, 'rotating quads before interp'
   !do i=1, 4
   !   print *,  'before', i, x_corners(i), y_corners(i)
   !enddo
   !print *, lat, lon
   do i = 2, 4
      x_corners(i) = x_corners(i) - x_corners(1)
      y_corners(i) = y_corners(i) - y_corners(1)
   enddo
   lon = lon - x_corners(1)
   lat = lat - y_corners(1)
   x_corners(1) = 0.0_r8
   y_corners(1) = 0.0_r8

   !do i=1, 4
   !   print *,  'xform ', i, x_corners(i), y_corners(i)
   !enddo
   !print *, lat, lon

   b(1) = x_corners(2)
   b(2) = y_corners(2)
   ! avoid degenerate cases where grid rotated
   ! exactly +/- 90 degrees.
   if (abs(x_corners(2)) > 0.001_r8) then
      c(1) = x_corners(2)
      c(2) = 0.0_r8
   else
      c(1) = 0.0_r8
      c(2) = y_corners(2)
   endif

!print *, b, c
   angle = angle2(b, c)
   !print *, 'angle = ', angle

   if (abs(angle) > 0.001_r8) then
   do i = 2, 4
     b(1) = x_corners(i)
     b(2) = y_corners(i)
     b = rotate2(b, angle)
     x_corners(i) = b(1)
     y_corners(i) = b(2)
   enddo
   b(1) = lon
   b(2) = lat
   b = rotate2(b, angle)
   lon = b(1)
   lat = b(2)
   endif
else
   !print *, 'NOT rotating quads before interp'
endif

! now everything is in degrees relative to the lower left and rotated.

if (debug > 10) write(*,'(A,5F15.5)') 'xformed x_corners, lon: ', x_corners, lon
if (debug > 10) write(*,'(A,5F15.5)') 'xformed y_corners, lat: ', y_corners, lat

!*******

! Fit a surface and interpolate; solve for 3x3 matrix
do i = 1, 3
   ! Eliminate a from the first 3 equations
   m(i, 1) = x_corners(i) - x_corners(i + 1)
   m(i, 2) = y_corners(i) - y_corners(i + 1)
   m(i, 3) = x_corners(i)*y_corners(i) - x_corners(i + 1)*y_corners(i + 1)
   v(i) = p(i) - p(i + 1)
if (debug > 10) write(*,'(A,I3,7F12.3)') 'i, m(3), p(2), v: ', i, m(i,:), p(i), p(i+1), v(i)
enddo

! look for degenerate matrix and rotate if needed
! compute deter of m
!d = deter3(m)

! Solve the matrix for b, c and d
call mat3x3(m, v, r)
if (debug > 10) print *, 'r ', r
if (debug > 10) print *, 'p ', p


! r contains b, c, and d; solve for a
a = p(4) - r(1) * x_corners(4) - &
           r(2) * y_corners(4) - &
           r(3) * x_corners(4)*y_corners(4)


!----------------- Implementation test block
! When interpolating on dipole x3 never exceeded 1e-9 error in this test
if (debug > 10)  write(*,'(A,8F12.3)') 'test corners: a, r(1), r(2), r(3)', a, r(1), r(2), r(3)
do i = 1, 4
   interp_val = a + r(1)*x_corners(i) + r(2)*y_corners(i)+ r(3)*x_corners(i)*y_corners(i)

   if(abs(interp_val - p(i)) > 1e-9) &
      write(*, *) 'large interp residual ', i, interp_val, p(i), interp_val - p(i)
if (debug > 10)  write(*,'(A,I3,8F12.5)') 'test corner: i, interp_val, x_corn, y_corn: ',  &
                                                        i, interp_val, x_corners(i), y_corners(i)
enddo

!----------------- Implementation test block


! Now do the interpolation

expected_obs = a + r(1)*lon + r(2)*lat + r(3)*lon*lat

if (debug > 10)  write(*,'(A,8F15.5)') 'poly: expected,     lon, lat, a,  r(1)*lon,  r(2)*lat,  r(3)*lon*lat: ', &
                                              expected_obs, lon, lat, a,  r(1)*lon,  r(2)*lat,  r(3)*lon*lat


!********
! Avoid exceeding maxima or minima as stopgap for poles problem
! When doing bilinear interpolation in quadrangle, can get interpolated
! values that are outside the range of the corner values
if(expected_obs > maxval(p)) then
!   expected_obs = maxval(p)
if (debug > 10)  write(*,'(A,3F12.3)') 'expected obs > maxval (diff): ', expected_obs, maxval(p), abs(expected_obs - maxval(p))
else if(expected_obs < minval(p)) then
!   expected_obs = minval(p)
if (debug > 10)  write(*,'(A,3F12.3)') 'expected obs < minval (diff): ', expected_obs, minval(p), abs(expected_obs - minval(p))
endif
!********

end subroutine quad_bilinear_interp

!------------------------------------------------------------
!> Solves rank 3 linear system mr = v for r using Cramer's rule.

subroutine mat3x3(m, v, r)

real(r8),  intent(in) :: m(3, 3), v(3)
real(r8), intent(out) :: r(3)

! Cramer's rule isn't the best choice
! for speed or numerical stability so might want to replace
! this at some point.

real(r8) :: m_sub(3, 3), numer, denom
integer  :: i

! Compute the denominator, det(m)
denom = deter3(m)

! Loop to compute the numerator for each component of r
do i = 1, 3
   m_sub = m
   m_sub(:, i) = v
   numer = deter3(m_sub)
   r(i) = numer / denom
if (debug > 10) write(*,'(A,I3,7F12.3)') 'mat: i, numer, denom, r: ', i, numer, denom, r(i)
enddo

end subroutine mat3x3

!------------------------------------------------------------
!> Computes determinant of 3x3 matrix m

function deter3(m)

real(r8), intent(in) :: m(3, 3)
real(r8)             :: deter3

deter3 = m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) + &
         m(1,3)*m(2,1)*m(3,2) - m(3,1)*m(2,2)*m(1,3) - &
         m(1,1)*m(2,3)*m(3,2) - m(3,3)*m(2,1)*m(1,2)

end function deter3

!------------------------------------------------------------
! Computes dot product of two 2-vectors

function dot2(a, b)
 real(r8), intent(in) :: a(2), b(2)
 real(r8)             :: dot2

dot2 = a(1)*b(1) + a(2)*b(2)

end function dot2

!------------------------------------------------------------
! compute the magnitude of a 2-vector

function mag2(a)
 real(r8), intent(in) :: a(2)
 real(r8)             :: mag2

mag2 = sqrt(a(1)*a(1) + a(2)*a(2))

end function mag2

!------------------------------------------------------------
! compute the angle between two 2-vectors

function angle2(a, b)
 real(r8), intent(in) :: a(2), b(2)
 real(r8)             :: angle2

angle2 = acos(dot2(a,b) / (mag2(a) * mag2(b)))

end function angle2

!------------------------------------------------------------
! rotate vector a counterclockwise by angle theta (in radians)

function rotate2(a, theta)
 real(r8), intent(in) :: a(2)
 real(r8), intent(in) :: theta
 real(r8)             :: rotate2(2)

real(r8) :: r(2,2)

r(1,1) = cos(theta)
r(1,2) = sin(theta)
r(2,1) = sin(-theta)
r(2,2) = cos(theta)

rotate2(1) = r(1,1)*a(1) + r(1,2)*a(2)
rotate2(2) = r(2,1)*a(1) + r(2,2)*a(2)

end function rotate2

!------------------------------------------------------------------

!> masked locations are invalid and cannot be the corner of any quad

function is_masked(opt, lon_index, lat_index)

type(quad_grid_options), intent(in) :: opt
integer,                 intent(in) :: lon_index, lat_index
logical                             :: is_masked

if (.not. opt%uses_mask) then
   is_masked = .false.
else
   is_masked = opt%grid_mask(lon_index, lat_index)
endif

end function is_masked

!------------------------------------------------------------------

function all_corners_valid(opt, lon_ind, lat_ind, nx)

type(quad_grid_options), intent(in) :: opt
integer,                 intent(in) :: lon_ind, lat_ind
integer,                 intent(in) :: nx
logical                             :: all_corners_valid

integer :: lon_ind_p1

! set to fail so we can return early.
all_corners_valid = .false.

! Might have to worry about wrapping in longitude but not in latitude
lon_ind_p1 = lon_ind + 1
if (opt%spans_lon_zero .and. lon_ind_p1 > nx) lon_ind_p1 = 1

if (is_masked(opt, lon_ind,    lat_ind  )) return
if (is_masked(opt, lon_ind_p1, lat_ind  )) return
if (is_masked(opt, lon_ind_p1, lat_ind+1)) return
if (is_masked(opt, lon_ind,    lat_ind+1)) return

all_corners_valid = .true.

end function all_corners_valid

!------------------------------------------------------------

! single item wrapper

subroutine quad_lon_lat_evaluate_ii_single(interp_handle, lon, lat, &
                      four_lons, four_lats, invals, outval, istatus)

type(quad_interp_handle), intent(in)  :: interp_handle
real(r8),                 intent(in)  :: lon, lat
integer,                  intent(in)  :: four_lons(4), four_lats(4)
real(r8),                 intent(in)  :: invals(4)
real(r8),                 intent(out) :: outval
integer,                  intent(out) :: istatus

real(r8) :: in_array(4, 1), out_array(1)

in_array(:, 1) = invals
call quad_lon_lat_evaluate_ii_array(interp_handle, lon, lat, four_lons, four_lats, &
                                    1, in_array, out_array, istatus)
outval = out_array(1)
istatus = 0

end subroutine quad_lon_lat_evaluate_ii_single

!------------------------------------------------------------

!> This is a different interface because you don't need fractions for the
!> irregular case.

subroutine quad_lon_lat_evaluate_ii_array(interp_handle, lon, lat, &
             four_lons, four_lats, nitems, invals, outvals, istatus)

type(quad_interp_handle), intent(in)  :: interp_handle
real(r8),                 intent(in)  :: lon, lat
integer,                  intent(in)  :: four_lons(4), four_lats(4)
integer,                  intent(in)  :: nitems
real(r8),                 intent(in)  :: invals(4, nitems)
real(r8),                 intent(out) :: outvals(nitems)
integer,                  intent(out) :: istatus

real(r8) :: x_corners(4), y_corners(4)
integer  :: e

character(len=*), parameter :: routine = 'quad_lon_lat_evaluate:quad_lon_lat_evaluate_ii_array'

! Full bilinear interpolation for quads
if(interp_handle%grid_type == GRID_QUAD_FULLY_IRREGULAR) then

   ! lons and lats are integer indices.  x_corners and y_corners are the real*8 locations.
   ! Get corner grid locations for accurate interpolation
   call get_quad_corners(interp_handle%ii%lons_2D, four_lons(1), four_lats(1), &
                         interp_handle%opt%spans_lon_zero, interp_handle%opt%pole_wrap, &
                         interp_handle%nlon, interp_handle%nlat, x_corners, istatus)
   if (istatus /= 0) return

   call get_quad_corners(interp_handle%ii%lats_2D, four_lons(1), four_lats(1), &
                         interp_handle%opt%spans_lon_zero, interp_handle%opt%pole_wrap, &
                         interp_handle%nlon, interp_handle%nlat, y_corners, istatus)
   if (istatus /= 0) return

   if (debug > 10) write(*,'(A,8F12.3)') 'evaluate: x_corners = ', x_corners
   if (debug > 10) write(*,'(A,8F12.3)') 'evaluate: y_corners = ', y_corners

   if (debug > 10) write(*,'(A,8F12.3)') 'evaluate: invals ens1 = ', invals(:, 1)
   do e = 1, nitems
      call quad_bilinear_interp(lon, lat, x_corners, y_corners, &
                        interp_handle%opt%spans_lon_zero, invals(:,e), outvals(e))
   enddo
   if (debug > 10) write(*,'(A,8F12.3)') 'evaluate: outvals ens1 = ', outvals(1)
else
   string1 = 'wrong interface for this grid'
   write(string2,*)'grid type is ',interp_handle%grid_type
   write(string3,*)'expected     ',GRID_QUAD_FULLY_IRREGULAR
   call error_handler(E_ERR, routine, string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

istatus = 0

end subroutine quad_lon_lat_evaluate_ii_array

!------------------------------------------------------------------
!> single item wrapper

subroutine quad_lon_lat_evaluate_ir_single(interp_handle, lon_fract, lat_fract, &
                                           invals, outval, istatus)

type(quad_interp_handle), intent(in)  :: interp_handle
real(r8),                 intent(in)  :: lon_fract, lat_fract
real(r8),                 intent(in)  :: invals(4)
real(r8),                 intent(out) :: outval
integer,                  intent(out) :: istatus

real(r8) :: in_array(4, 1), out_array(1)
integer  :: stat(1)

in_array(:, 1) = invals
call quad_lon_lat_evaluate_ir_array(interp_handle, lon_fract, lat_fract, &
                                    1, in_array, out_array, stat)
outval = out_array(1)
istatus = stat(1)

end subroutine quad_lon_lat_evaluate_ir_single

!------------------------------------------------------------

!> In the regular, orthogonal case you only need the fractions
!> across the quad at this point.

subroutine quad_lon_lat_evaluate_ir_array(interp_handle, lon_fract, lat_fract, &
                                          nitems, invals, outvals, istatus)

type(quad_interp_handle), intent(in)  :: interp_handle
real(r8),                 intent(in)  :: lon_fract, lat_fract
integer,                  intent(in)  :: nitems
real(r8),                 intent(in)  :: invals(4, nitems)
real(r8),                 intent(out) :: outvals(nitems)
integer,                  intent(out) :: istatus(nitems)

real(r8) :: xbot(nitems), xtop(nitems)
real(r8) :: x_corners(4), y_corners(4)
integer  :: i

character(len=*), parameter :: routine = 'quad_lon_lat_evaluate:quad_lon_lat_evaluate_ir_array'

! Full bilinear interpolation for quads
if(interp_handle%grid_type == GRID_QUAD_FULLY_IRREGULAR) then

   string1 = 'wrong interface for this grid'
   write(string2,*)'grid type is ',interp_handle%grid_type
   write(string3,*)'cannot be    ',GRID_QUAD_FULLY_IRREGULAR
   call error_handler(E_ERR, routine, string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

! Rectangular bilinear interpolation
!>@todo FIXME should this code check invals(:) for MISSING_R8?
!> it costs time and for grids that don't have missing data it is
!> not needed.  should it call allow_missing_in_state() on init and
!> key off that?  (i think yes.)

if (missing_ok_in_state) then

   ! have to do the items individually because some items might
   ! have missing and others not.
   do i=1, nitems
      if (any(invals(:, i) == MISSING_R8)) then
         outvals(i) = MISSING_R8
         istatus(i) = 1
      else
         xbot(1) = invals(1, i) + lon_fract * (invals(2, i) - invals(1, i))
         xtop(1) = invals(4, i) + lon_fract * (invals(3, i) - invals(4, i))
         outvals(i) = xbot(1) + lat_fract * (xtop(1) - xbot(1))
         istatus(i) = 0
      endif
   enddo
   return

else

   ! can use array syntax and do them all at once.  no missing vals.
   xbot = invals(1, :) + lon_fract * (invals(2, :) - invals(1, :))
   xtop = invals(4, :) + lon_fract * (invals(3, :) - invals(4, :))

   outvals(:) = xbot + lat_fract * (xtop - xbot)
   istatus(:) = 0

endif

end subroutine quad_lon_lat_evaluate_ir_array

!------------------------------------------------------------------

end module quad_utils_mod


! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
