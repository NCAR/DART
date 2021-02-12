! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Implements location interfaces for a three dimensional spherical shell 
!> with a choice of vertical coordinates.
!> Horizontal coordinates are always latitude and longitude.
!> Vertical coordinate choices include pressure, height, model level,
!> scale height, surface, and non-specific (column-integrated values, or 
!> with no logically defined vertical location, e.g. hurricane vortex center)
!> The internal representation of the location is stored as
!> radians from 0 to 2 PI for longitude and -PI/2 to PI/2 for latitude to
!> minimize computational cost for distances. However, the external 
!> representation is longitude in degrees from 0 to 360 and latitude 
!> from -90 to 90 for consistency with most applications in the field.
!>
!> This version supports multiple cutoff distances in an efficient manner.
!> Smaller cutoff values will do less searching than larger ones.  (This was
!> not true in earlier implementations of this code.)
!>
module location_mod

use      types_mod, only : r8, MISSING_R8, MISSING_I, PI, RAD2DEG, DEG2RAD, OBSTYPELENGTH, i8
use  utilities_mod, only : error_handler, E_ERR, ascii_file_format, &
                           E_MSG, open_file, close_file, set_output,                 &
                           logfileunit, nmlfileunit, find_namelist_in_file,          &
                           check_namelist_read, do_output, do_nml_file,              &
                           do_nml_term, is_longitude_between
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform
use   obs_kind_mod, only : get_num_types_of_obs, get_name_for_type_of_obs, get_index_for_type_of_obs
use mpi_utilities_mod, only : my_task_id, task_count
use ensemble_manager_mod, only : ensemble_type

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
          VERTISHEIGHT, VERTISSCALEHEIGHT, print_get_close_type


character(len=*), parameter :: source = 'threed_sphere/location_mod.f90'

integer,          parameter :: LocationDims = 3
character(len=*), parameter :: LocationName = "loc3Dsphere"
character(len=*), parameter :: LocationLName = &
                                   "threed sphere locations: lon, lat, vertical"
character(len=*), parameter :: LocationStorageOrder = "Lon Lat Vertical"
character(len=*), parameter :: LocationUnits = "degrees degrees which_vert"


! The possible numeric values for the location_type%which_vert component.
! The numeric values are PRIVATE to this module. The parameter names are PUBLIC.

integer, parameter :: VERTISUNDEF       = -2  ! has no specific vertical location (undefined)
integer, parameter :: VERTISSURFACE     = -1  ! surface value (value is surface elevation in m)
integer, parameter :: VERTISLEVEL       =  1  ! by level
integer, parameter :: VERTISPRESSURE    =  2  ! by pressure (in pascals)
integer, parameter :: VERTISHEIGHT      =  3  ! by height (in meters)
integer, parameter :: VERTISSCALEHEIGHT =  4  ! by scale height (unitless)

type location_type
   private
   real(r8) :: lon, lat        ! lon, lat are stored in radians
   real(r8) :: vloc            ! units vary based on value of which_vert
   integer  :: which_vert      ! determines if vert is level, height, pressure, ...
end type location_type

! Derived type to facilitate efficient computation of locations close to a given observation.

type get_close_type_by_type
   private
   integer               :: num
   real(r8)              :: maxdist              ! furthest separation between "close" locations
   integer, allocatable  :: lon_offset(:, :)     ! (nlat, nlat), lon box indices searched (nlat 2x IS correct)
   integer, allocatable  :: loc_box(:)           ! (nloc), List of loc indices in boxes
   integer, allocatable  :: count(:, :)          ! (nlon, nlat), # of loc in each box
   integer, allocatable  :: start(:, :)          ! (nlon, nlat), Start of list of loc in this box
   real(r8)              :: bot_lat, top_lat     ! Bottom and top latitudes of latitude boxes
   real(r8)              :: bot_lon, top_lon     ! Bottom and top longitudes of longitude boxes
   real(r8)              :: lon_width, lat_width ! Width of boxes in lon and lat
   logical               :: lon_cyclic           ! Do boxes wraparound in longitude?
end type get_close_type_by_type

! Support more than a single cutoff distance.  nt is the count of
! distinct cutoffs, which are selected per specific observation type.
! The map associates the incoming location type with the 
! corresponding gtt index.  There are only as many close_types
! as there are distinct cutoff distances; if more than one specific
! type has the same cutoff distances they share the type.

type get_close_type
   integer              :: nt                             ! number of distinct cutoffs
   integer, allocatable :: type_to_cutoff_map(:)          ! mapping of types to index
   type(get_close_type_by_type), allocatable :: gtt(:)    ! array of close_types by type
end type get_close_type

! Horizontal localization/cutoff values are passed in by the caller.
! The Vertical normalization values are globals; are specified by namelist
! here, and apply to all specific types unless a 'special list' is also specified
! that overrides the default values.

logical :: has_special_vertical_norms = .false.
integer :: num_special_vert_norms = 0

integer :: location_vertical_localization_coord = 0

! Some calls include information about the type or kind of the location. 
! If the location refers to an observation it is possible to have both
! a specific type and a generic kind.  If the location refers to a
! state vector item, it only has a generic kind.  Some routines have
! a specific type and a generic kind as arguments; look carefully at
! the argument names before using.

type(random_seq_type) :: ran_seq
logical               :: ran_seq_init = .false.
logical, save         :: module_initialized = .false.

character(len = 512) :: msgstring, msgstring1, msgstring2

! Global storage for vertical distance normalization factors
! The 4 below is for the 4 vertical units (pressure, level, height,
! scale height).  undefined and surface don't need vert distances.
! NOTE: Code that uses VERT_TYPE_COUNT depends on pressure, level,
! height, and scale height having actual values between 1 and 4, or
! this code will break.
integer, parameter    :: VERT_TYPE_COUNT = 4
real(r8)              :: vert_normalization(VERT_TYPE_COUNT)
real(r8), allocatable :: per_type_vert_norm(:,:)  ! if doing per-type

! Global storage for fast approximate sin and cosine lookups
! PAR For efficiency for small cases might want to fill tables as needed
! Also these could be larger for more accurate computations, if needed.
! 630 is 2 * PI rounded up, times 100.
real(r8), parameter :: SINCOS_DELTA = 100.0_r8
integer,  parameter :: SINCOS_LIMIT = 630
real(r8), parameter :: ACOS_DELTA = 1000.0_r8
integer,  parameter :: ACOS_LIMIT = 1000
real(r8) :: my_sin(-SINCOS_LIMIT:SINCOS_LIMIT)
real(r8) :: my_cos(-SINCOS_LIMIT:SINCOS_LIMIT)
real(r8) :: my_acos(-ACOS_LIMIT:ACOS_LIMIT)

! Tolerance for the top latitude boundary test.  All locations which
! are located on a box boundary are added to the bin on the larger
! side of the boundary (e.g. for north-south latitudes, it rounds
! towards the north).  But when a value falls exactly on the edge
! of the last box, technically it is inside the region but would be
! rounded up and outside the region unless handled specially.
! this tolerance below is used to determine if a value is within
! the range of the last box boundary and if so, the location is
! included in the last box.  In particular, for global grids this
! preserves locations which are at exactly 90.0 degrees latitude.
real(r8), parameter :: EDGE_TOLERANCE = 100.0_r8 * epsilon(0.0_r8)

! Option for verification using exhaustive search
logical :: COMPARE_TO_CORRECT = .false.    ! normally false

!-----------------------------------------------------------------
! Namelist with default values
! horiz_dist_only == .true.       -> Only the great circle horizontal distance is
!                                    computed in get_dist.
! horiz_dist_only == .false.      -> Square root of sum of squared horizontal and
!                                    normalized vertical dist computed in get_dist
! vert_normalization_pressure     -> Number pascals that give a distance equivalent
!                                    to one radian in horizontal
! vert_normalization_height       -> Number meters that give a distance equivalent 
!                                    to one radian in horizontal
! vert_normalization_level        -> Number levels that give a distance equivalent
!                                    to one radian in horizontal
! vert_normalization_scale_height -> Number scale heights that give a distance 
!                                    equivalent to one radian in horizontal
! approximate_distance            -> Use a faster table lookup for the trig math.
!                                    Works well for global models and large areas,
!                                    and improves performance.  For smaller regions
!                                    might introduce banding, so leave .false.
! nlon                            -> Number longitude boxes for get_close
!                                    nlon MUST BE ODD
! nlat                            -> Number latitude boxes for get_close
! output_box_info                 -> Useful for debugging performance problems.
! print_box_level                 -> How much data to print out.
! special_vert_normalization_obs_types -> Which obs types to modify the default vert
!                                    normalization values
! special_vert_normalization_pressure       -> must give all 4 values for each type listed
! special_vert_normalization_heights
! special_vert_normalization_levels
! special_vert_normalization_scale_heights


logical  :: horiz_dist_only                 = .true.
real(r8) :: vert_normalization_pressure     = 100000.0_r8
real(r8) :: vert_normalization_height       = 10000.0_r8
real(r8) :: vert_normalization_level        = 20.0_r8
real(r8) :: vert_normalization_scale_height = 5.0_r8
logical  :: approximate_distance            = .false.
integer  :: nlon                            = 71
integer  :: nlat                            = 36
logical  :: output_box_info                 = .false.
integer  :: print_box_level                 = 0

!>@todo make these allocatable instead of fixed length
integer, parameter :: MAX_ITEMS = 500
character(len=OBSTYPELENGTH) :: special_vert_normalization_obs_types(MAX_ITEMS)
real(r8) :: special_vert_normalization_pressures(MAX_ITEMS)
real(r8) :: special_vert_normalization_heights(MAX_ITEMS)
real(r8) :: special_vert_normalization_levels(MAX_ITEMS)
real(r8) :: special_vert_normalization_scale_heights(MAX_ITEMS)


namelist /location_nml/ horiz_dist_only, vert_normalization_pressure, &
   vert_normalization_height, vert_normalization_level,               &
   vert_normalization_scale_height, approximate_distance, nlon, nlat, &
   output_box_info, print_box_level, &
   special_vert_normalization_obs_types, special_vert_normalization_pressures, &
   special_vert_normalization_heights, special_vert_normalization_levels, &
   special_vert_normalization_scale_heights


!-----------------------------------------------------------------

interface operator(==); module procedure loc_eq; end interface
interface operator(/=); module procedure loc_ne; end interface

interface set_location
   module procedure set_location_single
   module procedure set_location_array
end interface set_location

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! basic location routines
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

! things which need doing exactly once.

subroutine initialize_module()

integer :: iunit, io, i, k, typecount, type_index


if (module_initialized) return

module_initialized = .true.

! give these initial values before reading them from the namelist.
special_vert_normalization_obs_types(:) = 'null'
special_vert_normalization_pressures(:)     = missing_r8
special_vert_normalization_heights(:)       = missing_r8
special_vert_normalization_levels(:)        = missing_r8
special_vert_normalization_scale_heights(:) = missing_r8


! Read the namelist entry
call find_namelist_in_file("input.nml", "location_nml", iunit)
read(iunit, nml = location_nml, iostat = io)
call check_namelist_read(iunit, io, "location_nml")

! Write the namelist values to the log file

if(do_nml_file()) write(nmlfileunit, nml=location_nml)
if(do_nml_term()) write(     *     , nml=location_nml)

call error_handler(E_MSG, 'location_mod:', 'using code with optimized cutoffs', &
                   source)

! Simple error checking for nlon, nlat to be sure we can use them for allocations.
if (nlon < 1 .or. nlat < 1) then
   write(msgstring, '(A,I5,A,I5)') 'from the namelist, nlon is ', nlon, '; nlat is ', nlat
   write(msgstring2, '(A)') 'both must be >= 1, and nlon must be odd'
   call error_handler(E_ERR, 'location_mod', msgstring, source, text2=msgstring2)
endif

! One additional restriction; number of longitudes, nlon, for get_close
if((nlon / 2) * 2 == nlon) then
   write(msgstring, '(A,I5,A)') 'nlon is ', nlon, '. Must be odd'
   call error_handler(E_ERR, 'location_mod', msgstring, source)
endif

! Copy the normalization factors in the vertical into an array
vert_normalization(VERTISLEVEL)       = vert_normalization_level
vert_normalization(VERTISPRESSURE)    = vert_normalization_pressure
vert_normalization(VERTISHEIGHT)      = vert_normalization_height
vert_normalization(VERTISSCALEHEIGHT) = vert_normalization_scale_height

! If the user has set a namelist to make different vertical normalization factors
! when computing the localization distances, allocate and set that up here.
! It overrides the defaults above.

if (special_vert_normalization_obs_types(1) /= 'null' .or. &
    special_vert_normalization_pressures(1)     /= missing_r8 .or. &
    special_vert_normalization_heights(1)       /= missing_r8 .or. &
    special_vert_normalization_levels(1)        /= missing_r8 .or. &
    special_vert_normalization_scale_heights(1) /= missing_r8) then

   ! FIXME: add code to check for mismatched length lists.  are we going to force
   ! users to specify all 4 values for any obs type that is not using the defaults?

   typecount = get_num_types_of_obs()  ! ignore function name, this is specific type count
   allocate(per_type_vert_norm(VERT_TYPE_COUNT, typecount))
  
   ! Set the defaults for all specific types not listed in the special list
   per_type_vert_norm(VERTISLEVEL, :)       = vert_normalization_level
   per_type_vert_norm(VERTISPRESSURE, :)    = vert_normalization_pressure
   per_type_vert_norm(VERTISHEIGHT, :)      = vert_normalization_height
   per_type_vert_norm(VERTISSCALEHEIGHT, :) = vert_normalization_scale_height

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

      per_type_vert_norm(VERTISLEVEL,       type_index) = special_vert_normalization_levels(i)
      per_type_vert_norm(VERTISPRESSURE,    type_index) = special_vert_normalization_pressures(i)
      per_type_vert_norm(VERTISHEIGHT,      type_index) = special_vert_normalization_heights(i)
      per_type_vert_norm(VERTISSCALEHEIGHT, type_index) = special_vert_normalization_scale_heights(i)
      
   enddo

   if (any(per_type_vert_norm == missing_r8)) then
      write(msgstring, *) 'one or more special vertical normalization values is uninitialized.'
      call error_handler(E_ERR,'location_mod:', &
                        'special vert normalization value namelist requires all 4 values per type', &
                        source, text2=trim(msgstring))
   endif

endif

! if the namelist control says we are only basing the localization on
! distances in the horizontal, log that in the dart log file.
if (horiz_dist_only) then
   call error_handler(E_MSG,'location_mod:', &
      'Ignoring vertical separation when computing distances; horizontal distances only', &
      source)
else
   call error_handler(E_MSG,'location_mod:', &
      'Including vertical separation when computing distances:', source)
   write(msgstring,'(A,f17.5)') '       # pascals ~ 1 horiz radian: ', vert_normalization_pressure
   call error_handler(E_MSG,'location_mod:',msgstring,source)
   write(msgstring,'(A,f17.5)') '        # meters ~ 1 horiz radian: ', vert_normalization_height
   call error_handler(E_MSG,'location_mod:',msgstring,source)
   write(msgstring,'(A,f17.5)') '  # model levels ~ 1 horiz radian: ', vert_normalization_level
   call error_handler(E_MSG,'location_mod:',msgstring,source)
   write(msgstring,'(A,f17.5)') ' # scale heights ~ 1 horiz radian: ', vert_normalization_scale_height
   call error_handler(E_MSG,'location_mod:',msgstring,source)

   if (allocated(per_type_vert_norm)) then
      typecount = get_num_types_of_obs()  ! ignore function name, this is specific type count
      do i = 1, typecount
         if ((per_type_vert_norm(VERTISLEVEL,       i) /= vert_normalization_level) .or. &
             (per_type_vert_norm(VERTISPRESSURE,    i) /= vert_normalization_pressure) .or. &
             (per_type_vert_norm(VERTISHEIGHT,      i) /= vert_normalization_height) .or. &
             (per_type_vert_norm(VERTISSCALEHEIGHT, i) /= vert_normalization_scale_height)) then
 
            write(msgstring,'(2A)') 'Altering default vertical normalization for type ', trim(get_name_for_type_of_obs(i))
            call error_handler(E_MSG,'location_mod:',msgstring,source)
            if (per_type_vert_norm(VERTISPRESSURE,    i) /= vert_normalization_pressure) then
               write(msgstring,'(A,f17.5)') '       # pascals ~ 1 horiz radian: ', &
                     per_type_vert_norm(VERTISPRESSURE, i)
               call error_handler(E_MSG,'location_mod:',msgstring,source)
            endif
            if (per_type_vert_norm(VERTISHEIGHT,      i) /= vert_normalization_height) then
               write(msgstring,'(A,f17.5)') '        # meters ~ 1 horiz radian: ', &
                     per_type_vert_norm(VERTISHEIGHT, i)
               call error_handler(E_MSG,'location_mod:',msgstring,source)
            endif
            if (per_type_vert_norm(VERTISLEVEL,       i) /= vert_normalization_level) then
               write(msgstring,'(A,f17.5)') '  # model levels ~ 1 horiz radian: ', &
                     per_type_vert_norm(VERTISLEVEL, i)
               call error_handler(E_MSG,'location_mod:',msgstring,source)
            endif
            if (per_type_vert_norm(VERTISSCALEHEIGHT, i) /= vert_normalization_scale_height) then
               write(msgstring,'(A,f17.5)') ' # scale heights ~ 1 horiz radian: ', &
                     per_type_vert_norm(VERTISSCALEHEIGHT, i)
               call error_handler(E_MSG,'location_mod:',msgstring,source)
            endif
         endif
      enddo
   endif
endif

! Set up a lookup table for cos and sin for approximate but fast distances
! Don't worry about rounding errors as long as one gives more weight
! Really only need tables half this size, too (sin from -pi/2 to pi/2, cos only +)
if(approximate_distance) then
   do i = -SINCOS_LIMIT, SINCOS_LIMIT
      my_cos(i) = cos(i / SINCOS_DELTA)
      my_sin(i) = sin(i / SINCOS_DELTA)
   end do
   do i = -ACOS_LIMIT, ACOS_LIMIT
      my_acos(i) = acos(i / ACOS_DELTA)
   end do
   call error_handler(E_MSG,'location_mod:', &
      'Using table-lookup approximation for distance computations', source)
endif

end subroutine initialize_module

!----------------------------------------------------------------------------

function get_dist(loc1, loc2, type1, kind2, no_vert)

! returns the distance between 2 locations in units of radians.
! Distance depends on only horizontal distance or both horizontal and 
! vertical distance. 
! The choice is determined by horiz_dist_only and the which_vert of loc1.
! May want to allow return of some joint distance in the long run? 
! Or just a distance that is a function of all 3 things.
! The namelist controls whether default computations use just horizontal distance.
! However, this behavior can be over-ridden by the no_vert optional argument.
! If set to false, this will always do full 3d distance if possible. If set to
! true it will never do the full 3d distance. At present asking to do a vertical
! distance computation for incompatible vertical location units results 
! in a fatal error unless one of the vertical units is UNDEFINED.

! CHANGE from previous versions:  the 3rd argument is now a specific type
! (e.g. RADIOSONDE_TEMPERATURE, AIRCRAFT_SPECIFIC_HUMIDITY) associated
! with loc1, while the 4th argument is a generic kind (QTY_TEMPERATURE, 
! QTY_U_WIND_COMPONENT) associated with loc2.
! The type and kind are part of the interface in case user-code wants to do 
! a more sophisticated distance calculation based on the base type or target
! kind. In the usual case this code still doesn't use the kind/type, but 
! it does require at least type1 if using per-type vertical normalization.

type(location_type), intent(in) :: loc1, loc2
integer, optional,   intent(in) :: type1, kind2
logical, optional,   intent(in) :: no_vert
real(r8)                        :: get_dist

real(r8) :: lon_dif, vert_dist, rtemp
integer  :: lat1_ind, lat2_ind, lon_ind, temp  ! indexes into lookup tables
logical  :: comp_h_only


if ( .not. module_initialized ) call initialize_module()

! Begin with the horizontal distance
! Compute great circle path shortest route between two points
lon_dif = loc1%lon - loc2%lon

if(approximate_distance) then
   ! Option 1: Use table lookup; faster but less accurate
   lat1_ind = int(loc1%lat*SINCOS_DELTA)
   lat2_ind = int(loc2%lat*SINCOS_DELTA)
   lon_ind  = int(lon_dif *SINCOS_DELTA)
   temp     = int(ACOS_DELTA * (my_sin(lat2_ind) * my_sin(lat1_ind) + &
                  my_cos(lat2_ind) * my_cos(lat1_ind) * my_cos(lon_ind)))
   get_dist = my_acos(temp)
else
   ! Option 2: Use pre-defined trig functions: accurate but slow
   ! First 2 ifs avoids round-off error that can kill acos;
   if(abs(loc1%lat) >= PI/2.0_r8 .or. abs(loc2%lat) >= PI/2.0_r8 .or. &
      lon_dif == 0.0_r8) then
      get_dist = abs(loc2%lat - loc1%lat)
   else
      ! This test is for apparent roundoff error which may be a result of
      ! running r8 == r4. 
      rtemp = sin(loc2%lat) * sin(loc1%lat) + &
              cos(loc2%lat) * cos(loc1%lat) * cos(lon_dif)
      if (rtemp < -1.0_r8) then
         get_dist = PI
      else if (rtemp > 1.0_r8) then
         get_dist = 0.0_r8
      else
         get_dist = acos(rtemp)
      endif
   endif
endif

! Now compute a vertical distance if requested.  Highest priority is
! the optional no_vert argument, so test it first.
if(present(no_vert)) then
   comp_h_only = no_vert
! Namelist horizontal only has second highest priority
else 
   comp_h_only = horiz_dist_only
endif

! If which_vert has no vertical definition for either location do only horizontal
if(loc1%which_vert == VERTISUNDEF .or. loc2%which_vert == VERTISUNDEF) comp_h_only = .true.
! If both verts are surface, do only horizontal
if(loc1%which_vert == VERTISSURFACE .and. loc2%which_vert == VERTISSURFACE) comp_h_only = .true.

! Add in vertical component if required
if(.not. comp_h_only) then
   ! Vert distance can only be done for like vertical locations units
   if(loc1%which_vert /= loc2%which_vert) then
      write(msgstring,*)'loc1%which_vert (',loc1%which_vert, &
                   ') /= loc2%which_vert (',loc2%which_vert,')'
      call error_handler(E_MSG, 'get_dist', msgstring, source)
      call write_location(logfileunit,loc1)
      call write_location(logfileunit,loc2)
      call write_location(0,loc1)
      call write_location(0,loc2)
      call error_handler(E_ERR, 'get_dist', &
         'Dont know how to compute vertical distance for unlike vertical coordinates', &
         source)
   endif

   ! Compute the difference and divide by the appropriate normalization factor
   ! Normalization factor computes relative distance in vertical compared to one radian
   ! This is new - if per-type localization distances given, use the specific type of loc1
   ! to determine the vertical mapping distance.  it defaults to the 4 standard ones,
   ! but can be specified separately if desired.

   ! note that per-type vertical conversion factors are a global here, and were set
   ! by namelist values.  they apply to all calls to get_dist() based on the obs type.
   if (allocated(per_type_vert_norm)) then 
      if (.not. present(type1)) then
         write(msgstring, *) 'obs type required in get_dist() if doing per-type vertical normalization'
         call error_handler(E_MSG, 'get_dist', msgstring, source)
      endif 
      vert_dist = abs(loc1%vloc - loc2%vloc) / per_type_vert_norm(loc1%which_vert, type1)
   else
      vert_dist = abs(loc1%vloc - loc2%vloc) / vert_normalization(loc1%which_vert)
   endif

   ! Spherical distance shape is computed here, other flavors can be computed
   get_dist = sqrt(get_dist**2 + vert_dist**2)
endif

end function get_dist

!---------------------------------------------------------------------------

function loc_eq(loc1,loc2)
 
! Interface operator used to compare two locations.
! Returns true only if all components are 'the same' to within machine
! precision. There is some debate as to whether or not the vertical
! locations need to be identical if 'VERTISUNDEF' ... hard to peruse
! the code tree to find where this may be affected. 

type(location_type), intent(in) :: loc1, loc2
logical                         :: loc_eq

if ( .not. module_initialized ) call initialize_module()

loc_eq = .false.

! if ( loc1%which_vert /= loc2%which_vert ) return
if ( abs(loc1%lon  - loc2%lon) > epsilon(loc1%lon) ) return
if ( abs(loc1%lat  - loc2%lat) > epsilon(loc1%lat) ) return

!if ( loc1%which_vert /= VERTISUNDEF ) then
   if ( abs(loc1%vloc - loc2%vloc) > epsilon(loc1%vloc) ) return
!endif

loc_eq = .true.

end function loc_eq

!---------------------------------------------------------------------------

function loc_ne(loc1,loc2)
 
! Interface operator used to compare two locations.
! Returns true if locations are not identical to machine precision.

type(location_type), intent(in) :: loc1, loc2
logical                         :: loc_ne

if ( .not. module_initialized ) call initialize_module()

loc_ne = (.not. loc_eq(loc1,loc2))

end function loc_ne

!---------------------------------------------------------------------------

function get_location(loc)
 
! Given a location type (in radians), 
! return the longitude, latitude (in degrees) and vertical value

type(location_type), intent(in) :: loc
real(r8), dimension(3) :: get_location

if ( .not. module_initialized ) call initialize_module()

get_location(1) = loc%lon * RAD2DEG                 
get_location(2) = loc%lat * RAD2DEG                 
get_location(3) = loc%vloc     

end function get_location

!---------------------------------------------------------------------------

function set_location_single(lon, lat, vert_loc,  which_vert)
 
! Puts the given longitude, latitude, and vertical location
! into a location datatype.  Arguments to this function are in degrees,
! but the values are stored as radians.

real(r8), intent(in) :: lon, lat
real(r8), intent(in) :: vert_loc
integer,  intent(in) :: which_vert
type (location_type) :: set_location_single

if ( .not. module_initialized ) call initialize_module()

if(lon < 0.0_r8 .or. lon > 360.0_r8) then
   write(msgstring,*)'longitude (',lon,') is not within range [0,360]'
   call error_handler(E_ERR, 'set_location', msgstring, source)
endif

if(lat < -90.0_r8 .or. lat > 90.0_r8) then
   write(msgstring,*)'latitude (',lat,') is not within range [-90,90]'
   call error_handler(E_ERR, 'set_location', msgstring, source)
endif

set_location_single%lon = lon * DEG2RAD
set_location_single%lat = lat * DEG2RAD

if(which_vert < VERTISUNDEF .or. which_vert == 0 .or. which_vert > VERTISSCALEHEIGHT) then
   write(msgstring,*)'which_vert (',which_vert,') must be one of -2, -1, 1, 2, 3, or 4'
   call error_handler(E_ERR,'set_location', msgstring, source)
endif

set_location_single%vloc = vert_loc
set_location_single%which_vert = which_vert

end function set_location_single

!----------------------------------------------------------------------------

function set_location_array(list)
 
! location semi-independent interface routine
! given 4 float numbers, call the underlying set_location routine

real(r8), intent(in) :: list(:)
type (location_type) :: set_location_array

if ( .not. module_initialized ) call initialize_module()

if (size(list) < 4) then
   write(msgstring,*)'requires 4 input values'
   call error_handler(E_ERR, 'set_location', msgstring, source)
endif

set_location_array = set_location_single(list(1), list(2), list(3), nint(list(4)))

end function set_location_array

!----------------------------------------------------------------------------

function set_location_missing()

! Initialize a location type to indicate the contents are unset.

type (location_type) :: set_location_missing

if ( .not. module_initialized ) call initialize_module()

set_location_missing%lon        = MISSING_R8
set_location_missing%lat        = MISSING_R8
set_location_missing%vloc       = MISSING_R8
set_location_missing%which_vert = MISSING_I

end function set_location_missing

!---------------------------------------------------------------------------

function query_location(loc, attr)
 
! Returns the value of the attribute

type(location_type),        intent(in) :: loc
character(len=*), optional, intent(in) :: attr
real(r8)                               :: query_location

if ( .not. module_initialized ) call initialize_module()

! Workaround for apparent bug in mac osx intel 10.x fortran compiler.
! Previous code had a 16 byte local character variable which was
! apparently not getting deallocated after this function returned.

! Workaround for a second compiler bug.  The PGI 6.1.x compiler
! refused to compile the select statement below when it was
! selectcase(adjustl(attr)).  i'm rearranging the code again because
! this particular routine has been so troublesome.  i'm removing
! the result(fval) construct, and setting the default to which_vert
! and then overwriting it the case we recognize another quantity.
! if we fall through to the end of the routine, the return value
! is which_vert.  (This routine is rapidly becoming my favorite
! problem child.  How can such a short piece of code be so troublesome?)

! set the default here, and then only overwrite it if we
! recognize one of the other valid queries.

query_location = real(loc%which_vert, r8)  ! this is really an int

if (.not. present(attr)) return

select case(attr)
   case ('lat','LAT')
      query_location = loc%lat
   case ('lon','LON')
      query_location = loc%lon
   case ('vloc','VLOC')
      query_location = loc%vloc
   case ('which_vert','WHICH_VERT')
      ! already set
   case default
      call error_handler(E_ERR, 'query_location:', &
         'Only "lon","lat","vloc","which_vert" are legal attributes to request from location', &
          source)
end select

end function query_location

!----------------------------------------------------------------------------

subroutine write_location(locfile, loc, fform, charstring)
 
! Writes a location to a file.
! most recent change: adding the optional charstring option.  if present,
! locfile is ignored, and a pretty-print formatting is done into charstring.
! the locations are converted to lat/lon, and the vert is put into more
! common units (e.g. hPa, km, etc)

integer, intent(in)                        :: locfile
type(location_type), intent(in)            :: loc
character(len = *),  intent(in),  optional :: fform
character(len = *),  intent(out), optional :: charstring

integer             :: charlength
logical             :: writebuf
character(len=129)  :: string1

! 10 format(1x,3(f22.14,1x),i4)  ! old
10 format(1X,3(G25.16,1X),I2)

if ( .not. module_initialized ) call initialize_module()

! writing to a file (normal use) or to a character buffer?
writebuf = present(charstring)

! output file; test for ascii or binary, write what's asked, and return
if (.not. writebuf) then
   if (ascii_file_format(fform)) then
      write(locfile, '(''loc3d'')' ) 
      write(locfile, 10) loc%lon, loc%lat, loc%vloc, loc%which_vert
   else
      write(locfile) loc%lon, loc%lat, loc%vloc, loc%which_vert
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
! degrees instead of radians, and kilometers for height,
! hectopascals instead of pascals for pressure, etc.

! this must be the sum of the longest of the formats below.
charlength = 72

if (len(charstring) < charlength) then
   write(msgstring, *) 'charstring buffer must be at least ', charlength, ' chars long'
   call error_handler(E_ERR, 'write_location', msgstring, source)
endif

! format the horizontal into a temp string
write(string1, '(A,F12.8,1X,F12.8,1X,A)') 'Lon/Lat(deg): ',  loc%lon*RAD2DEG, &
   loc%lat*RAD2DEG, ' Vert:'

! then pretty up the vertical choices, trying to get them to line up in
! case the caller is listing out locations with different vert units.
! concatinate the vertical on the end of the horizontal and put it all
! into the return string. 
select case  (loc%which_vert)
   case (VERTISUNDEF)
      write(charstring, '(A,A)')       trim(string1), '              Undefined'
   case (VERTISSURFACE)
      write(charstring, '(A,F13.5,A)') trim(string1), loc%vloc, ' surface (m)'
   case (VERTISLEVEL)
      write(charstring, '(A,F13.6,A)') trim(string1), loc%vloc, '  level'
   case (VERTISPRESSURE)
      write(charstring, '(A,F13.7,A)') trim(string1), loc%vloc / 100.0_r8, ' hPa'
   case (VERTISHEIGHT)
      write(charstring, '(A,F13.7,A)') trim(string1), loc%vloc / 1000.0_r8, ' km'
   case (VERTISSCALEHEIGHT)
      write(charstring, '(A,F13.7,A)') trim(string1), loc%vloc, ' scale ht'
   case default
      write(msgstring, *) 'unrecognized key for vertical type: ', loc%which_vert
      call error_handler(E_ERR, 'write_location', msgstring, source)
end select


end subroutine write_location

!----------------------------------------------------------------------------

function read_location(locfile, fform)
 
! Reads a location from a file that was written by write_location. 
! See write_location for additional discussion.

integer, intent(in)                      :: locfile
character(len = *), intent(in), optional :: fform
type(location_type)                      :: read_location

character(len=5) :: header

if ( .not. module_initialized ) call initialize_module()

if (ascii_file_format(fform)) then
   read(locfile, '(a5)' ) header
   if(header /= 'loc3d') then
         write(msgstring,*)'Expected location header "loc3d" in input file, got ', header
      call error_handler(E_ERR, 'read_location', msgstring, source)
   endif
   ! Now read the location data value
   read(locfile, *)read_location%lon, read_location%lat, &
                   read_location%vloc, read_location%which_vert
else
   read(locfile)read_location%lon, read_location%lat, &
                read_location%vloc, read_location%which_vert
endif

end function read_location

!--------------------------------------------------------------------------

subroutine interactive_location(location, set_to_default)
 
! Allows for interactive input of a location. Also gives option of selecting
! a uniformly distributed random location.

type(location_type), intent(out) :: location
logical, intent(in), optional    :: set_to_default

real(r8) :: lon, lat, minlon, maxlon, minlat, maxlat

if ( .not. module_initialized ) call initialize_module()

! If set_to_default is true, then just zero out and return
if(present(set_to_default)) then
   if(set_to_default) then
      location%lon = 0.0
      location%lat = 0.0
      location%vloc = 0.0
      location%which_vert = 0   ! note that 0 isn't a valid vert type
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
   location%vloc = 100.0 * location%vloc
else if(location%which_vert == VERTISHEIGHT ) then
   write(*, *) 'Vertical coordinate height (in meters)'
   read(*, *) location%vloc
else if(location%which_vert == VERTISSURFACE ) then
   ! most 3d sphere users want height in meters, not pressure.
   ! original code asked for pressure:
   !write(*, *) 'Vertical coordinate surface pressure (in hPa)'
   !location%vloc = 100.0 * location%vloc  ! only applies to pressure
   write(*, *) 'Vertical coordinate height'
   read(*, *) location%vloc
else if(location%which_vert == VERTISSCALEHEIGHT ) then
   write(*, *) 'Vertical coordinate scale height (-ln(p/ps))'
   read(*, *) location%vloc
else if(location%which_vert == VERTISUNDEF ) then
   ! a valid floating point value, but should be unused.
   location%vloc = MISSING_R8
else
   write(*, *) 'Wrong choice of which_vert try again between ',VERTISUNDEF, &
               ' and ',VERTISHEIGHT
   go to 100
end if

write(*, *) 'Input longitude: value 0 to 360.0 or a negative number for '
write(*, *) 'Uniformly distributed random location in the horizontal'
read(*, *) lon

do while(lon > 360.0_r8)
   write(*, *) 'Input value greater than 360.0 is illegal, please try again'
   read(*, *) lon
end do

if(lon < 0.0_r8) then

   ! Need to make sure random sequence is initialized

   if(.not. ran_seq_init) then
      call init_random_seq(ran_seq)
      ran_seq_init = .TRUE.
   endif

   write(*, *) 'Input minimum longitude (0 to 360.0)'
   read(*, *) minlon
   minlon = minlon * DEG2RAD

   write(*, *) 'Input maximum longitude (0 to 360.0)'
   read(*, *) maxlon
   maxlon = maxlon * DEG2RAD

   ! Longitude is random from minlon to maxlon
   location%lon = random_uniform(ran_seq) * (maxlon-minlon) + minlon

   write(*, *) 'Input minimum latitude (-90.0 to 90.0)'
   read(*, *) minlat
   minlat = sin(minlat * DEG2RAD)

   write(*, *) 'Input maximum latitude (-90.0 to 90.0)'
   read(*, *) maxlat
   maxlat = sin(maxlat * DEG2RAD)

   ! Latitude must be area weighted
   location%lat = asin(random_uniform(ran_seq) * (maxlat-minlat) + minlat)

   write(*, *) 'random location is ', location%lon / DEG2RAD, &
                                      location%lat / DEG2RAD

else

   write(*, *) 'Input latitude: value -90.0 to 90.0'
   read(*, *) lat

   do while(lat < -90.0_r8 .or. lat > 90.0_r8)
      write(*, *) 'Input value < -90.0 or > 90.0 is illegal, please try again'
      read(*, *) lat
   end do

   location%lat = lat*DEG2RAD
   location%lon = lon*DEG2RAD

end if

end subroutine interactive_location

!----------------------------------------------------------------------------

function is_location_in_region(loc, minl, maxl)
 
! Returns true if the given location is inside the rectangular
! region defined by minl as the lower left, maxl the upper right.
! test is inclusive; values on the edges are considered inside.
! Periodic in longitude (box can cross the 2PI -> 0 line)

logical                          :: is_location_in_region
type(location_type), intent(in)  :: loc, minl, maxl

if ( .not. module_initialized ) call initialize_module()

! maybe could use VERTISUNDEF in the minl and maxl args to indicate
! we want to test only in horizontal?  and if not, vtypes must match?
!if ( (minl%which_vert /= maxl%which_vert) .or. &
! ((minl%which_vert /= loc%which_vert).and.(minl%which_vert /= VERTISUNDEF))) then
!   write(msgstring,*)'which_vert (',loc%which_vert,') must be same in all args'
!   call error_handler(E_ERR, 'is_location_in_region', msgstring, source)
!endif

! assume failure and return as soon as we are confirmed right.
! set to success only at the bottom after all tests have passed.
is_location_in_region = .false.

! latitude: we do not allow wrap of rectangular regions over the poles.
if ((loc%lat < minl%lat) .or. (loc%lat > maxl%lat)) return

! use common routine in utilities module to do all the wrapping
if (.not. is_longitude_between(loc%lon, minl%lon, maxl%lon, doradians=.TRUE.)) return

! once we decide what to do about diff vert units, this is the test.
!if ((minl%which_vert .ne. VERTISUNDEF) .and. 
!    (loc%vloc < minl%vloc) .or. (loc%vloc > maxl%vloc)) return
 
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

vertical_localization_on = .not. horiz_dist_only

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
      write(msgstring, *) 'unrecognized key for vertical type: ', which_vert
      call error_handler(E_ERR, 'is_vertical', msgstring, source)
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
!----------------------------------------------------------------------------
! get close routines
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine get_close_init(gc, num, maxdist, locs, maxdist_list)
 
! Initializes part of get_close accelerator dependent on the particular location

type(get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
real(r8),             intent(in)    :: maxdist
type(location_type),  intent(in)    :: locs(:)
real(r8), intent(in), optional      :: maxdist_list(:)

integer :: blat_ind, tlat_ind
integer :: bot_tlat_ind, top_tlat_ind
real(r8) :: base_lat(2), target_lat(2), del_lon, cos_del_lon
real(r8) :: max_del_lon
integer :: i, j, n, cum_start, lon_box(num), lat_box(num), tstart(nlon, nlat), tj, bj
integer :: typecount, distcount
real(r8), allocatable :: distlist(:)

if ( .not. module_initialized ) call initialize_module()

! Support per-loc-type localization more efficiently.
typecount = get_num_types_of_obs()  ! ignore function name, this is specific type count
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

! make a gtt type array for each different cutoff distance
! (a single type is the most common case)
allocate(gc%gtt(gc%nt))

if (present(maxdist_list)) then
   do i=1, gc%nt
      gc%gtt(i)%maxdist = distlist(i)
   enddo
else
   ! no per-type settings, everyone uses same distance
   gc%gtt(1)%maxdist = maxdist
endif

if (present(maxdist_list)) deallocate(distlist)   ! temp storage

! Allocate the storage for the grid dependent boxes
do i=1, gc%nt
   allocate(gc%gtt(i)%count(nlon, nlat), gc%gtt(i)%start(nlon, nlat))
   allocate(gc%gtt(i)%lon_offset(nlat, nlat))
   gc%gtt(i)%lon_offset = -1
   gc%gtt(i)%count      = -1
   gc%gtt(i)%start      = -1
enddo

! store the location counts in all derived types
do i=1, gc%nt
   gc%gtt(i)%num = num
enddo

! If there are no locs to operate on, no point in going any further.
if (num == 0) return

do i=1, gc%nt
   ! Allocate storage for locs number dependent part 
   allocate(gc%gtt(i)%loc_box(num))
   gc%gtt(i)%loc_box(:) = -1
enddo

do n=1, gc%nt
   ! Determine where the boxes should be for this set of locs and maxdist
   call find_box_ranges(gc%gtt(n), num, locs)

   ! Figure out which boxes are close to a box on a given latitude circle
   ! MIGHT AVOID DOING THIS WITH A COPY ROUTINE: HAVE SAME BOXES OFTEN
   do blat_ind = 1, nlat
      ! Search from east side of base block
      ! Start searching out, have to look for closest point in box being checked
      ! Only have to search latitude boxes that are within maximum distance
      bot_tlat_ind = blat_ind - floor(gc%gtt(n)%maxdist / gc%gtt(n)%lat_width) - 1
      if(bot_tlat_ind < 1) bot_tlat_ind = 1
      top_tlat_ind = blat_ind + floor(gc%gtt(n)%maxdist / gc%gtt(n)%lat_width) + 1
      if(top_tlat_ind > nlat) top_tlat_ind = nlat
      do tlat_ind = bot_tlat_ind, top_tlat_ind
         ! Spherical geometry can be tricky
         ! We want to find the MINIMUM distance between two lat/lon boxes
         ! This distance is known to be corners of the boxes. It is also known
         ! to be between the corners such that the longitude difference between
         ! the corners is a minimum. HOWEVER, determining whether it is between
         ! the closest latitudes or not is a non-trivial computation. Hence,
         ! since this isn't done much, we just check all four possible combinations
         ! of latitude and pick the one that gives the closest distance.
         do j = 1, 2
            base_lat(j)   = gc%gtt(n)%bot_lat + (blat_ind - 2 + j) * gc%gtt(n)%lat_width
            target_lat(j) = gc%gtt(n)%bot_lat + (tlat_ind - 2 + j) * gc%gtt(n)%lat_width
         end do
   
         ! If the max distance > PI, then everybody is close.
         ! Do a test for something slightly less than PI to avoid round-off error.
         ! Set max_del_lon to something much larger than PI since it doesn't matter.
         if(gc%gtt(n)%maxdist > PI - 0.0001_r8) then
            max_del_lon = 2.0 * PI
         else
            ! Find the maximum longitude offset for the different possible latitudes edges
            max_del_lon = 0.0_r8
            do tj = 1, 2
               do bj = 1, 2
                  ! Compute the lon offset directly by inverting distance
                  cos_del_lon = (cos(gc%gtt(n)%maxdist) - sin(base_lat(bj)) * sin(target_lat(tj))) / &
                     (cos(base_lat(bj)) * cos(target_lat(tj)))
                  if(cos_del_lon < -1.0_r8) then
                     del_lon = PI
                  else if(cos_del_lon > 1.0_r8) then
                     del_lon = 0.0_r8
                  else
                     del_lon = acos(cos_del_lon)
                  endif
                  if(del_lon > max_del_lon) max_del_lon = del_lon
               
               end do
            end do
         endif
         
         ! Compute the number of boxes to search in longitude for maximum del_lon
         gc%gtt(n)%lon_offset(blat_ind, tlat_ind) = floor(max_del_lon / gc%gtt(n)%lon_width) + 1
         ! Watch for roundoff leading to a search of more offsets than exist
         if(gc%gtt(n)%lon_offset(blat_ind, tlat_ind) > nlon / 2) &
            gc%gtt(n)%lon_offset(blat_ind, tlat_ind) = nlon / 2
   
      end do
   end do

   ! Begin by computing the number of locations in each box in lat/lon
   gc%gtt(n)%count = 0
   do i = 1, num
      lon_box(i) = get_lon_box(gc%gtt(n), locs(i)%lon)
      if(lon_box(i) < 0 .or. lon_box(i) > nlon) then
         write(msgstring, *) 'Contact Dart Developers: this error should not happen'
         call error_handler(E_MSG, 'get_close_init', msgstring, source)
         write(msgstring, *) 'location outside grid boxes, index value:',  lon_box(i)
         call error_handler(E_ERR, 'get_close_init', msgstring, source)
      endif
   
      lat_box(i) = floor((locs(i)%lat - gc%gtt(n)%bot_lat) / gc%gtt(n)%lat_width) + 1
      if(lat_box(i) > nlat) lat_box(i) = nlat
      if(lat_box(i) < 1) lat_box(i) = 1
   
      gc%gtt(n)%count(lon_box(i), lat_box(i)) = gc%gtt(n)%count(lon_box(i), lat_box(i)) + 1
   end do
   
   ! Figure out where storage for each boxes members should begin
   cum_start = 1
   do i = 1, nlon
      do j = 1, nlat
         gc%gtt(n)%start(i, j) = cum_start
         cum_start = cum_start + gc%gtt(n)%count(i, j)
      end do
   end do
   
   ! Now we know how many are in each box, get a list of which are in each box
   tstart = gc%gtt(n)%start
   do i = 1, num
      gc%gtt(n)%loc_box(tstart(lon_box(i), lat_box(i))) = i
      tstart(lon_box(i), lat_box(i)) = tstart(lon_box(i), lat_box(i)) + 1
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

! order matters here.  do all arrays inside the gtt types
! before releasing the array of types
do i=1, gc%nt
  if (allocated(gc%gtt(i)%loc_box)) deallocate(gc%gtt(i)%loc_box)
  deallocate(gc%gtt(i)%lon_offset, gc%gtt(i)%count, gc%gtt(i)%start)
enddo
deallocate(gc%type_to_cutoff_map)
deallocate(gc%gtt)

end subroutine get_close_destroy

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

! If dist is NOT present, just find everybody in a box, put them in the list,
! but don't compute any distances

integer :: lon_box, lat_box, i, j, k, n_lon, lon_ind, n_in_box, st, t_ind, bt
real(r8) :: this_dist, this_maxdist

! Variables needed for comparing against correct case
integer :: cnum_close, cclose_ind(size(locs))
real(r8) :: cdist(size(locs))


! First, set the intent out arguments to a missing value
num_close = 0
close_ind = -99
if(present(dist)) dist = -99.0_r8
this_dist = 999999.0_r8   ! something big.

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

! the list of locations in the locs() argument must be the same
! as the list of locations passed into get_close_init(), so
! gc%gtt(bt)%num and size(locs) better be the same.   if the list 
! changes you have to destroy the old gc and init a new one.
if (size(locs) /= gc%gtt(bt)%num) then
   write(msgstring,'(A)')'locs() array must match one passed to get_close_init()'
   write(msgstring1,'(A,2I8)')'init size, current list size: ', gc%gtt(bt)%num, size(locs)
   write(msgstring2,'(A,I8)')'bt = ', bt
   call error_handler(E_ERR, 'get_close', msgstring, source, &
                      text2=msgstring1, text3=msgstring2)
endif

! If num == 0, no point in going any further. 
if (gc%gtt(bt)%num == 0) return


! local variable for what the maxdist is in this particular case.
this_maxdist = gc%gtt(bt)%maxdist

!--------------------------------------------------------------
! For validation, it is useful to be able to compare against exact
! exhaustive search
if(COMPARE_TO_CORRECT) then
   cnum_close = 0
   do i = 1, gc%gtt(bt)%num 
      if (locs(i)%which_vert /= base_loc%which_vert) then
         this_dist = get_dist(base_loc, locs(i), base_type, loc_qtys(i), &
                     no_vert = .true.)
      else
         this_dist = get_dist(base_loc, locs(i), base_type, loc_qtys(i))
      endif
      if(this_dist <= this_maxdist) then
         ! Add this location to correct list
         cnum_close = cnum_close + 1
         cclose_ind(cnum_close) = i
         cdist(cnum_close) = this_dist
      endif
   end do
endif

!--------------------------------------------------------------

! Begin by figuring out which box the base_ob is in.
! Note that the boxes will not cover the entire sphere for sets of locations 
! that cover less than a hemisphere. In that case you expect to get out of range
! box numbers if the location is outside the covered region. But you do have to 
! handle the case of a location which is exactly on the top latitude boundary.  
! Include it in the top box if it is in the next box but still within N*epsilon 
! of the lower edge.
! FIXME: the longitude tests don't try to be smart about this, probably
! because of the complication of wrapping around the greenwich line.
! it just bins all values in box nlon+1 into box nlon.  we could do 
! the same thing here.
lon_box = get_lon_box(gc%gtt(bt), base_loc%lon)
lat_box = floor((base_loc%lat - gc%gtt(bt)%bot_lat) / gc%gtt(bt)%lat_width) + 1
! FIXME: cheaper to test now or later?
!if ((lat_box == nlat+1) .and. &
!    (base_loc%lat <= gc%gtt(bt)%top_lat + EDGE_TOLERANCE)) then
!!DEBUG write(*,'(A,4(G25.16,1X))') 'add to top lat_box', base_loc%lat, gc%gtt(bt)%top_lat, &
!!DEBUG                             epsilon(0.0_r8), (base_loc%lat - gc%gtt(bt)%top_lat)/epsilon(0.0_r8)
!   lat_box = nlat
!endif
! consistent with the current longitude code in get_lon_box():
if (lat_box == nlat+1) lat_box = nlat

!  If it is not in any box, then it is more than the maxdist away from everybody
if(lat_box > nlat .or. lat_box < 1 .or. lon_box < 0) return

! Next, loop through to find each box that is close to this box
do j = 1, nlat
   n_lon = gc%gtt(bt)%lon_offset(lat_box, j)
   if(n_lon >= 0) then
      LON_OFFSET: do i = -1 * n_lon, n_lon
         lon_ind = lon_box + i
         ! Search a box at this latitude j offset in longitude by i
         ! If domain is cyclic, need to wrap around
         if(gc%gtt(bt)%lon_cyclic) then
            if(lon_ind > nlon) lon_ind = lon_ind - nlon
            if(lon_ind < 1) lon_ind = lon_ind + nlon
         else
            ! Domain is not cyclic, don't search if outside of range
            if(lon_ind > nlon .or. lon_ind < 1) cycle LON_OFFSET
         endif
         ! Box to search is lon_ind, j
         n_in_box = gc%gtt(bt)%count(lon_ind, j)
         st = gc%gtt(bt)%start(lon_ind, j)
         ! Loop to check how close all locs in the box are; add those that are close
         do k = 1, n_in_box

            ! SHOULD ADD IN OPTIONAL ARGUMENT FOR DOING THIS!!!
            ! Could avoid adding any that have nums lower than base_ob???
            t_ind = gc%gtt(bt)%loc_box(st - 1 + k)

            if(.not. present(dist)) then
               ! Dist isn't present; add this ob to list without computing distance
               num_close = num_close + 1
               close_ind(num_close) = t_ind
            else
               if(base_loc%which_vert == locs(t_ind)%which_vert) then
                  ! Can compute total distance here if verts are the same
                  this_dist = get_dist(base_loc, locs(t_ind), base_type, loc_qtys(t_ind))
               else 
                  ! Otherwise can just get horizontal distance
                  this_dist = get_dist(base_loc, locs(t_ind), base_type, loc_qtys(t_ind), &
                     no_vert = .true.)
               endif

               ! If this locations distance is less than cutoff, add it to the list
               if(this_dist <= this_maxdist) then
                  num_close = num_close + 1
                  close_ind(num_close) = t_ind
                  dist(num_close) = this_dist
               endif
            endif
         end do
      end do LON_OFFSET
   endif
end do

!------------------------ Verify by comparing to exhaustive search --------------
if(COMPARE_TO_CORRECT) then
   ! Do comparisons against full search
   if((num_close /= cnum_close) .and. present(dist)) then
      write(msgstring, *) 'get_close (', num_close, ') should equal exhaustive search (', cnum_close, ')'
      call error_handler(E_ERR, 'get_close', msgstring, source, &
           text2='optional arg "dist" is present; we are computing exact distances', &
           text3='the exhaustive search should find an identical number of locations')
   else if (num_close < cnum_close) then
      write(msgstring, *) 'get_close (', num_close, ') should not be smaller than exhaustive search (', cnum_close, ')'
      call error_handler(E_ERR, 'get_close', msgstring, source, &
           text2='optional arg "dist" not present; we are returning a superset of close locations', &
           text3='the exhaustive search should find an equal or lesser number of locations')
   endif
endif
!--------------------End of verify by comparing to exhaustive search --------------

end subroutine get_close

!--------------------------------------------------------------------------

subroutine find_box_ranges(gtt, num, locs)
 
! Finds boundaries for boxes in N/S direction. If data is localized in N/S
! tries to find boxes that only span the range of the data.
  
type(get_close_type_by_type), intent(inout) :: gtt
integer,                      intent(in)    :: num
type(location_type),          intent(in)    :: locs(num)

real(r8) :: min_lat, max_lat, beg_box_lon, end_box_lon, first_loc_lon, last_loc_lon
real(r8) :: longitude_range, degrees, lon_dist
integer  :: i, indx, gap_start, gap_end, gap_length
logical  :: lon_box_full(360)

! Initialize boxes used to see where locations are.
! Assume regional until we prove that we have to use the
! entire 360 in longitude.
lon_box_full = .false.
gtt%lon_cyclic = .false.

! Figure out domain over which an additional locs MIGHT be close to one in this set
! If any points within maxdist of the poles, our boxes have to cover all 360 of
! longitude - no point in trying to restrict boxes to a region of the globe.
min_lat = minval(locs(:)%lat) - gtt%maxdist
max_lat = maxval(locs(:)%lat) + gtt%maxdist
if(min_lat <= -PI / 2.0_r8) then
   min_lat = -PI / 2.0_r8
   gtt%lon_cyclic = .true.
endif 
if(max_lat >= PI / 2.0_r8) then
   max_lat = PI / 2.0_r8
   gtt%lon_cyclic = .true.
endif

! Put this into storage for this get_close_type
gtt%bot_lat = min_lat
gtt%top_lat = max_lat
gtt%lat_width = (max_lat - min_lat) / nlat
! don't have to do all this work if we already know it has to be cyclic
if (.not. gtt%lon_cyclic) then

   ! Finding the longitude range is tricky because of cyclic nature
   ! Want to find minimum range spanned by locs even if they wrap-around Greenwich
   ! Would like to do this without sorting if possible at low-cost
   ! First, partition into 360 1-degree boxes and find the biggest gap
   do i = 1, num
      degrees = locs(i)%lon * 180.0_r8 / PI
      ! If the value of the longitude is very close to an integer number of degrees
      ! a roundoff can occur that leads to an assignment in the wrong box.  We avoid this
      ! by first testing to see if this is possible and then setting both boxes to full.
      ! If this is not the case, then we fill the box the location is in.
      if (abs(degrees - nint(degrees)) < 0.00001_r8) then
         indx = nint(degrees)
         if(indx <   1) indx = 360
         lon_box_full(indx) = .true.
   
         indx = nint(degrees) + 1
         if(indx > 360) indx = 1
         lon_box_full(indx) = .true.
      else
         indx = floor(degrees) + 1
         lon_box_full(indx) = .true.
      endif
   end do
   
   ! Find the longest sequence of empty boxes
   call find_longest_gap(lon_box_full, 360, gap_start, gap_end, gap_length)
   
   if (gap_length > 0) then

      ! There is a gap; figure out locs that are closest to ends of non-gap
      beg_box_lon = (gap_end / 180.0_r8) * PI
      end_box_lon = ((gap_start -1) / 180.0_r8) * PI
      first_loc_lon = find_closest_to_start(beg_box_lon, num, locs)
      last_loc_lon  = find_closest_to_end  (end_box_lon, num, locs)
      ! Determine the final longitude range
      longitude_range = last_loc_lon - first_loc_lon
      if(longitude_range <= 0.0_r8) longitude_range = longitude_range + 2.0_r8 * PI
      
      ! Add on the extra distance needed for the boxes
      ! To avoid any hard thinking about wraparound with sub-domain boxes
      ! Must span less than 180 degrees to get smaller boxes
      ! If addition of halos for close location fills more than half of space 
      ! things go 0 to 2PI

      ! other places we are computing in radians.  here we are computing in
      ! lat/lon, and you can't just add maxdist to the edges - that doesn't
      ! take into account the great-circle distance.  the separation in longitude
      ! varies with latitude.  compute the delta longitude based on the most
      ! poleward latitude and add that onto the edges of both boxes.  that
      ! overestimates for points closer to the equator, but that's better
      ! than underestimating and excluding points that are within maxdist.
      lon_dist = find_del_lon(minval(locs(:)%lat), maxval(locs(:)%lat), gtt%maxdist)

      if(longitude_range + 2.0_r8 * lon_dist > PI) then
         gtt%lon_cyclic = .true.
      else
         first_loc_lon = first_loc_lon - lon_dist
         if(first_loc_lon < 0.0_r8    ) first_loc_lon = first_loc_lon + 2.0_r8 * PI
         last_loc_lon  = last_loc_lon + lon_dist
         if(last_loc_lon > 2.0_r8 * PI) last_loc_lon  = last_loc_lon  - 2.0_r8 * PI
         gtt%lon_cyclic = .false.
      endif
   else
      ! No gap was found: all 360 boxes had an location in them
      gtt%lon_cyclic = .true.
   endif
endif

if (gtt%lon_cyclic) then
   first_loc_lon = 0.0_r8
   last_loc_lon  = 2.0_r8 * PI
endif

! Put in storage for structure
gtt%bot_lon = first_loc_lon
gtt%top_lon = last_loc_lon
longitude_range = last_loc_lon - first_loc_lon
if(longitude_range <= 0.0_r8) longitude_range = longitude_range + 2.0_r8 * PI
gtt%lon_width = longitude_range / nlon

end subroutine find_box_ranges

!----------------------------------------------------------------------------

subroutine find_longest_gap(lon_box_full, num_boxes, gap_start, gap_end, gap_length)
 
! Find the longest gap in the boxes (take the first one if there's a tie)

integer, intent(in) :: num_boxes
integer, intent(out) :: gap_start, gap_end, gap_length
logical, intent(in)  :: lon_box_full(num_boxes)

integer :: g_start, g_end, g_length, next_box, i, full_count
logical :: all_done

! Initialize these to known values.
gap_start  = -1
gap_end    = -1
gap_length = -1

! If more than half of the boxes are full, then assume that there is no
! meaningful gap and just use the whole domain for get_close
full_count = 0
do i = 1, num_boxes
   if(lon_box_full(i)) full_count = full_count + 1
end do
if(full_count >= num_boxes / 2) return
 
! More than half of the boxes were empty, try to hone in on potentially
! local locations
next_box = 1
all_done = .false.
! Loop long enough to be sure we go around
do i = 1, num_boxes
   call find_next_gap(next_box, lon_box_full, num_boxes, g_start, g_end, g_length) 
   next_box = g_end
   ! Easy way to terminate at cost of some small additional computation
   if(g_start == gap_start .and. g_end == gap_end) all_done = .true.
   if(g_length > gap_length) then
      gap_start  = g_start
      gap_end    = g_end
      gap_length = g_length
   endif
   ! Exit if longest has been found twice
   if(all_done) return
end do

end subroutine find_longest_gap

!----------------------------------------------------------------------------

subroutine find_next_gap(start_box, lon_box_full, num_boxes, gap_start, gap_end, gap_length)
integer, intent(in)  :: start_box, num_boxes
integer, intent(out) :: gap_start, gap_end, gap_length
logical, intent(in)  :: lon_box_full(num_boxes)

integer :: next_full

! This never gets called unless at least half of the boxes are empty
! No need to error check for all boxes empty (means num is 0) or for
! all boxes full.

! Finds the next gap of empty boxes in the cyclic set
! First, find the next full box from the start
next_full = next_full_box(start_box, num_boxes, lon_box_full)
! Find the next empty after that, make it the start of the gap
gap_start = next_empty_box(next_full, num_boxes, lon_box_full)
! Find the next full, box before that is the end of the gap
gap_end = next_full_box(gap_start, num_boxes, lon_box_full) - 1
if(gap_end < 1) gap_end = gap_end + num_boxes
! Carefully compute gap length
if(gap_end >= gap_start) then
   gap_length = gap_end - gap_start + 1
else
   gap_length = gap_end - gap_start + 360 + 1
endif

end subroutine find_next_gap


!----------------------------------------------------------------------------

function next_full_box(start_box, num_boxes, lon_box_full)
 
integer, intent(in) :: start_box
integer, intent(in) :: num_boxes
logical, intent(in) :: lon_box_full(num_boxes)
integer             :: next_full_box

integer :: i, indx

do i = 0, num_boxes
   indx = start_box + i
   if(indx > num_boxes) indx = indx - num_boxes
   if(lon_box_full(indx)) then
      next_full_box = indx
      return
   endif
end do

! Should never fall off the end since all boxes should not be empty
! Fatal error if this happens
call error_handler(E_ERR, 'next_full_box', 'All boxes empty: should not happen', source)

end function next_full_box

!----------------------------------------------------------------------------

function next_empty_box(start_box, num_boxes, lon_box_full)
 
integer, intent(in) :: start_box
integer, intent(in) :: num_boxes
logical, intent(in) :: lon_box_full(num_boxes)
integer             :: next_empty_box

integer :: i, indx

do i = 0, num_boxes
   indx = start_box + i
   if(indx > num_boxes) indx = indx - num_boxes
   if(.not. lon_box_full(indx)) then
      next_empty_box = indx
      return
   endif
end do

! Should never fall off the end since all boxes should not be full
! Fatal error if this happens
call error_handler(E_ERR, 'next_empty_box', 'All boxes full: should not happen', source)

end function next_empty_box

!----------------------------------------------------------------------------

function find_closest_to_start(beg_box_lon, num, locs)
 
real(r8),            intent(in) :: beg_box_lon
integer,             intent(in) :: num
type(location_type), intent(in) :: locs(num)
real(r8)                        :: find_closest_to_start

real(r8) :: least_dist, dist
integer  :: i

! Start with large value
least_dist = 2.0_r8 * PI

do i = 1, num
   dist = locs(i)%lon - beg_box_lon
   if(dist < 0.0_r8) dist = dist + 2.0_r8 * PI
   if(dist < least_dist) then
      least_dist = dist
      find_closest_to_start = locs(i)%lon
   endif 
end do

end function find_closest_to_start

!----------------------------------------------------------------------------

function find_closest_to_end(end_box_lon, num, locs)
 
real(r8),            intent(in) :: end_box_lon
integer,             intent(in) :: num
type(location_type), intent(in) :: locs(num)
real(r8)                        :: find_closest_to_end

real(r8) :: least_dist, dist
integer  :: i

! Start with large value
least_dist = 2.0_r8 * PI

do i = 1, num
   dist = end_box_lon - locs(i)%lon
   if(dist < 0.0_r8) dist = dist + 2.0_r8 * PI
   if(dist < least_dist) then
      least_dist = dist
      find_closest_to_end = locs(i)%lon
   endif 
end do

end function find_closest_to_end


!----------------------------------------------------------------------------

function get_lon_box(gtt, lon)

type(get_close_type_by_type), intent(in) :: gtt
real(r8),                     intent(in) :: lon
integer                                  :: get_lon_box
 
real(r8) :: del_lon

del_lon = lon - gtt%bot_lon
if(del_lon < 0.0_r8) del_lon = del_lon + 2.0_r8 * PI
get_lon_box = floor(del_lon / gtt%lon_width) + 1
! On wraparound, correct for truncation
! If not wraparound, then we're not in one of the boxes
if(get_lon_box > nlon) then
   if(gtt%lon_cyclic) then
      get_lon_box = 1
   else
      ! technically this should have a tolerance and only
      ! bin items which are on the boundary of boxes nlon and
      ! nlon+1.  this bins points which won't be close enough.
      ! FIXME: evaluate which is cheaper; 1) binning and computing
      ! distances on points later which can be excluded now, or
      ! 2) computing now whether the points are on the boundary
      ! and only keeping those ones.
      if(get_lon_box == nlon+1) then
         get_lon_box = nlon
      else
         get_lon_box = -1
      endif
   endif
endif

end function get_lon_box


!----------------------------------------------------------------------------

function find_del_lon(minlat, maxlat, maxdist)
 
! for the given latitudes, find the furthest longitude that is still
! within maxdist away.  this will be at a different latitude at any
! location other than the equator.  all values specified in radians.
! distance returned in radians.  if either lat is closer to the pole
! than maxdist, it returns 2*PI.

real(r8), intent(in) :: minlat, maxlat, maxdist
real(r8)             :: find_del_lon

real(r8) :: a, b, c
real(r8) :: latval, poleward_lat

! find the most poleward of the two latitudes
poleward_lat = max(abs(minlat), abs(maxlat))

! if either latitude is within maxdist of either pole, return 2 PI
! because you are now covering all possible longitudes.
if (poleward_lat + maxdist > (PI / 2.0_r8)) then
   find_del_lon = 2.0_r8 * PI
   return
endif

! compute some values we will reuse a couple times
a = cos(maxdist)
b = sin(poleward_lat)
c = cos(poleward_lat)

! lat at which max offset is found
latval = asin(b/a)

! distance to furthest lon, at latval
find_del_lon = acos((a - (b*sin(latval))) / (c*cos(latval)))

end function find_del_lon

!---------------------------------------------------------------------------
!> returns the maximum distance for the cutoff specified for the 
!> observation type of interest.
!> May be useful in custom 'get_close' applications.

function get_maxdist(gc, obs_type)
type(get_close_type), intent(in) :: gc
integer, optional,    intent(in) :: obs_type
real(r8) :: get_maxdist

integer :: bt

bt = gc%type_to_cutoff_map(obs_type)
get_maxdist = gc%gtt(bt)%maxdist

end function get_maxdist

!----------------------------------------------------------------------------

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

!----------------------------------------------------------------------------

subroutine print_get_close_type(gc, tt, amount)
 
! FIXME:  this is very useful for debugging but the current state
! of this code is ATROCIOUS!!  fix it.   possibly take the val=8
! details out into a separate subroutine so the simple case isn't so
! hard to follow.
!
! print out debugging statistics, or optionally print out a full
! dump from all mpi tasks in a format that can be plotted with matlab.

type(get_close_type), intent(in) :: gc
integer, intent(in), optional    :: tt
integer, intent(in), optional    :: amount

integer :: i, j, k, first, myindex, mytask, alltasks, whichtt
integer :: sample, nfull, nempty, howmuch, total, maxcount, maxi, maxj
logical :: tickmark(gc%gtt(1)%num), iam0
real(r8) :: lon_cen, lat_cen

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

! third arg is now an int, not logical, and means:
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
! locations from the locations in another.  this lets you turn off
! the debugging level for the large set and leave it on for the small.
!if (gc%gtt(whichtt)%num > 100) howmuch = 0

! print the get_close_type derived type values

if (howmuch /= 0 .and. iam0) then
   write(msgstring,*) 'get_close_type values:'
   call error_handler(E_MSG, 'loc', msgstring)

   write(msgstring,*) ' num = ', gc%gtt(whichtt)%num
   call error_handler(E_MSG, 'loc', msgstring)

   write(msgstring,*) ' nlon, nlat = ', nlon, nlat
   call error_handler(E_MSG, 'loc', msgstring)

   write(msgstring, "(A,F12.6)") ' maxdist = ', gc%gtt(whichtt)%maxdist
   call error_handler(E_MSG, 'loc', msgstring)
   write(msgstring, "(A,3(F12.6))") ' latbox: bot, top, width = ', gc%gtt(whichtt)%bot_lat, gc%gtt(whichtt)%top_lat, gc%gtt(whichtt)%lat_width
   call error_handler(E_MSG, 'loc', msgstring)
   write(msgstring, "(A,3(F12.6))") ' lonbox: bot, top, width = ', gc%gtt(whichtt)%bot_lon, gc%gtt(whichtt)%top_lon, gc%gtt(whichtt)%lon_width
   call error_handler(E_MSG, 'loc', msgstring)

   write(msgstring, "(A,F12.6)") ' maxdist = ', RAD2DEG*gc%gtt(whichtt)%maxdist
   call error_handler(E_MSG, 'loc', msgstring)
   write(msgstring, "(A,3(F12.6))") ' latbox: bot, top, width = ', RAD2DEG*gc%gtt(whichtt)%bot_lat, RAD2DEG*gc%gtt(whichtt)%top_lat, RAD2DEG*gc%gtt(whichtt)%lat_width
   call error_handler(E_MSG, 'loc', msgstring)
   write(msgstring, "(A,3(F12.6))") ' lonbox: bot, top, width = ', RAD2DEG*gc%gtt(whichtt)%bot_lon, RAD2DEG*gc%gtt(whichtt)%top_lon, RAD2DEG*gc%gtt(whichtt)%lon_width
   call error_handler(E_MSG, 'loc', msgstring)

   write(msgstring,*) ' lon_cyclic = ', gc%gtt(whichtt)%lon_cyclic
   call error_handler(E_MSG, 'loc', msgstring)
endif

! this one can be very large.   print only the first nth unless
! instructed otherwise.  (print n+1 because 1 more value fits on
! the line because it prints ( i ) and not ( i, j ) like the others.)
if (allocated(gc%gtt(whichtt)%loc_box)) then
   i = size(gc%gtt(whichtt)%loc_box,1)
   if (i/= gc%gtt(whichtt)%num) then
      write(msgstring,*) ' warning: size of loc_box incorrect, nlocs, i =', gc%gtt(whichtt)%num, i
      call error_handler(E_MSG, 'locations_mod', msgstring)
   endif
   if (howmuch > 1) then
      ! DEBUG
      write(msgstring,"(A,I8,A,36(I8,1X))") ' loc_box(',i,') =', gc%gtt(whichtt)%loc_box(1:min(i,36))  ! (nlocs)
      !write(msgstring,*) ' loc_box(',i,') =', gc%gtt(whichtt)%loc_box    ! (nlocs)
      call error_handler(E_MSG, 'locations_mod', msgstring)
   else if(howmuch > 0) then
      write(msgstring,*) ' loc_box(',i,') =', gc%gtt(whichtt)%loc_box(1:min(i,sample+1))
      call error_handler(E_MSG, 'locations_mod', msgstring)
      write(msgstring,*) '  <rest of loc_box omitted>'
      call error_handler(E_MSG, 'locations_mod', msgstring)
   endif
else
   if (howmuch > 0) then
      write(msgstring,*) ' loc_box unallocated'
      call error_handler(E_MSG, 'locations_mod', msgstring)
   endif
endif

! like loc_box, this one can be very large.   print only the first nth unless
! instructed otherwise
if (allocated(gc%gtt(whichtt)%start)) then
   i = size(gc%gtt(whichtt)%start,1)
   j = size(gc%gtt(whichtt)%start,2)
   if ((i /= nlon) .or. (j /= nlat)) then
      write(msgstring,*) ' warning: size of start incorrect, nlon, nlat, i, j =', nlon, nlat, i, j
      call error_handler(E_MSG, 'locations_mod', msgstring)
   endif
   if (howmuch > 1) then
      write(msgstring,*) ' start(',i,j,') ='              ! (nlon, nlat)
      call error_handler(E_MSG, 'locations_mod', msgstring)
      do k=1, j
         write(msgstring,"(36(I8,1X))") gc%gtt(whichtt)%start(1:min(i,36), k)
         call error_handler(E_MSG, 'locations_mod', msgstring)
      enddo
   else if (howmuch > 0) then
      write(msgstring,*) ' start(',i,j,') =', gc%gtt(whichtt)%start(1:min(i,sample), 1)
      call error_handler(E_MSG, 'locations_mod', msgstring)
      write(msgstring,*) '  <rest of start omitted>'
      call error_handler(E_MSG, 'locations_mod', msgstring)
   endif
else
   if (howmuch > 0) then
      write(msgstring,*) ' start unallocated'
      call error_handler(E_MSG, 'locations_mod', msgstring)
   endif
endif

! as above, print only first n unless second arg is .true.
if (allocated(gc%gtt(whichtt)%lon_offset)) then
   i =  size(gc%gtt(whichtt)%lon_offset,1)
   j =  size(gc%gtt(whichtt)%lon_offset,2)
   if ((i /= nlat) .or. (j /= nlat)) then
      write(msgstring,*) ' warning: size of lon_offset incorrect, nlat, i, j =', nlat, i, j
      call error_handler(E_MSG, 'locations_mod', msgstring)
   endif
   if (howmuch > 1) then
      write(msgstring,*) ' lon_offset(',i,j,') ='                 ! (nlat, nlat)
      call error_handler(E_MSG, 'locations_mod', msgstring)
      do k=1, j
         write(msgstring,"(36(I8,1X))") gc%gtt(whichtt)%lon_offset(1:min(i,36), k) 
         call error_handler(E_MSG, 'locations_mod', msgstring)
      enddo
   else if (howmuch > 0) then
      write(msgstring,*) ' lon_offset(',i,j,') =', gc%gtt(whichtt)%lon_offset(1:min(i,sample), 1)
      call error_handler(E_MSG, 'locations_mod', msgstring)
      write(msgstring,*) '  <rest of lon_offset omitted>'
      call error_handler(E_MSG, 'locations_mod', msgstring)
   endif
else
   if (howmuch > 0) then
      write(msgstring,*) ' lon_offset unallocated'
      call error_handler(E_MSG, 'locations_mod', msgstring)
   endif
endif

! as above, print only first n unless second arg is .true.
if (allocated(gc%gtt(whichtt)%count)) then
   i = size(gc%gtt(whichtt)%count,1)
   j = size(gc%gtt(whichtt)%count,2)
   if ((i /= nlon) .or. (j /= nlat)) then
      write(msgstring,*) ' warning: size of count incorrect, nlon, nlat, i, j =', &
                      nlon, nlat, i, j
      call error_handler(E_MSG, 'locations_mod', msgstring)
   endif
   if (howmuch > 1) then
      write(msgstring,*) ' count(',i,j,') ='              ! (nlon, nlat)
      call error_handler(E_MSG, 'locations_mod', msgstring)
      do k=1, j
         write(msgstring,"(36(I8,1X))") gc%gtt(whichtt)%count(1:min(i,36), k) 
         call error_handler(E_MSG, 'locations_mod', msgstring)
      enddo
   else if (howmuch > 0) then
      write(msgstring,*) ' count(',i,j,') =', gc%gtt(whichtt)%count(1:min(i,sample), 1)
      call error_handler(E_MSG, 'locations_mod', msgstring)
      write(msgstring,*) '  <rest of count omitted>'
      call error_handler(E_MSG, 'locations_mod', msgstring)
   endif
else
   if (howmuch > 0) then
      write(msgstring,*) ' count unallocated'
      call error_handler(E_MSG, 'locations_mod', msgstring)
   endif
endif


! end of print.  strictly speaking the following code is validation code,
! not a simple print.


! initialize all ticks to false.  turn them true as they are found
! in the loc_box list, and complain about duplicates or misses.
tickmark = .FALSE.

do i=1, nlon
   do j=1, nlat
      first = gc%gtt(whichtt)%start(i, j)
      do k=1, gc%gtt(whichtt)%count(i, j)
         myindex = first + k - 1
         if ((myindex < 1) .or. (myindex > gc%gtt(whichtt)%num)) then
            write(msgstring, *) 'exiting at first bad value; could be more'
            call error_handler(E_MSG, 'locations_mod', msgstring)
            write(msgstring, *) 'bad loc list index, in box: ', myindex, i, j
            call error_handler(E_ERR, 'locations_mod', msgstring)
         endif
         if (tickmark(myindex)) then
            write(msgstring, *) 'exiting at first bad value; could be more'
            call error_handler(E_MSG, 'locations_mod', msgstring)
            write(msgstring, *) 'error: locs found in more than one box list.  index, box: ', &
                         myindex, i, j
            call error_handler(E_ERR, 'locations_mod', msgstring)
         endif
         tickmark(myindex) = .TRUE.
      enddo
   enddo
enddo

do i=1, gc%gtt(whichtt)%num
  if (.not. tickmark(i)) then
     write(msgstring, *) 'exiting at first bad value; could be more'
     call error_handler(E_MSG, 'locations_mod', msgstring)
     write(msgstring,*) 'locs not found in any box list: ', i
     call error_handler(E_ERR, 'locations_mod', msgstring)
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
      write(funit,'(A,I2,A,I4,A)') 'xlocs = zeros(', nlon, ',', alltasks, ');'
      write(funit,'(A,I2,A,I4,A)') 'ylocs = zeros(', nlat, ',', alltasks, ');'
      write(funit,'(A,I2,A,I2,A,I4,A)') 'boxes = zeros(', nlon, ',', nlat, ',', alltasks, ');'
      call close_file(funit)
   endif
   write(fname, '(A,I3.3,A)')  'loc_dump_', mytask, '.m'
   funit = open_file(fname, action='write')
endif

do i=1, nlon
   if (howmuch == -8) then
      lon_cen = gc%gtt(whichtt)%bot_lon + ((i-1)*gc%gtt(whichtt)%lon_width) + (gc%gtt(whichtt)%lon_width/2.0)
      write(funit, '(A,I2,A,I4,A,F12.9,A)') 'xlocs(', i, ',', mytask+1, ') = ',  lon_cen, ';'
   endif
   do j=1, nlat
      if (howmuch == -8 .and. i==1) then
         lat_cen = gc%gtt(whichtt)%bot_lat + ((j-1)*gc%gtt(whichtt)%lat_width) + (gc%gtt(whichtt)%lat_width/2.0)
         write(funit, '(A,I2,A,I4,A,F12.9,A)') 'ylocs(', j, ',', mytask+1, ') = ',  lat_cen, ';'
      endif
      if (gc%gtt(whichtt)%count(i, j) > 0) then
         nfull = nfull + 1
         total = total + gc%gtt(whichtt)%count(i, j)
         if (gc%gtt(whichtt)%count(i, j) > maxcount) then
            maxcount = gc%gtt(whichtt)%count(i, j)
            maxi = i
            maxj = j
         endif
      else
         nempty = nempty + 1
      endif
      ! output for grid boxes; in matlab-friendly format
      if (howmuch == -8) then
         write(funit, '(A,I2,A,I2,A,I4,A,I8,A)') 'boxes(', i, ', ', j, &
                                ',', mytask+1, ') = ', gc%gtt(whichtt)%count(i, j), ';'
      endif
   enddo
enddo

if (howmuch == -8) then
   call close_file(funit)
   write_now = .false.
endif

! these print out always - make sure they are useful to end users.
write(msgstring, '(a)') "Location module statistics:"
call error_handler(E_MSG, 'locations_mod', msgstring)
write(msgstring, '(a,i9)') " Total boxes (nlon * nlat): ", nfull + nempty
call error_handler(E_MSG, 'locations_mod', msgstring)
write(msgstring, '(a,i9)') " Total items to put in boxes: ", gc%gtt(tt)%num
call error_handler(E_MSG, 'locations_mod', msgstring)
if (howmuch > 0) then
   write(msgstring, '(a,i9)') " Total boxes with 1+ items: ", nfull
   call error_handler(E_MSG, 'locations_mod', msgstring)
   write(msgstring, '(a,i9)') " Total boxes empty: ", nempty
   call error_handler(E_MSG, 'locations_mod', msgstring)
endif
if (nfull > 0) then
   write(msgstring, '(a,f7.2)') " Percent boxes with 1+ items: ", nfull / real(nfull + nempty, r8) * 100.
   call error_handler(E_MSG, 'locations_mod', msgstring)
   write(msgstring, '(a,f12.2)') " Average #items per non-empty box: ", real(total, r8) / nfull
   call error_handler(E_MSG, 'locations_mod', msgstring)
endif
if (maxcount > 0) then
   write(msgstring, '(a,i9)') " Largest #items in one box: ", maxcount
   call error_handler(E_MSG, 'locations_mod', msgstring)
! leave this out for now.  one, if there are multiple boxes with
! the same maxcount this is just the last one found.  two, the
! index numbers do not seem very helpful.
!   if (howmuch > 0) then
!      write(msgstring, '(a,i9,i9)') " That box index: ", maxi, maxj
!      call error_handler(E_MSG, 'locations_mod', msgstring)
!   endif
endif


end subroutine print_get_close_type

!----------------------------------------------------------------------------
! end of location/threed_sphere/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod

