! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module location_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! Implements location interfaces for a three dimensional spherical shell 
! with a pressure vertical coordinate plus
! a vertical coordinate based on the models native set of
! discrete levels. In the long run, it would be nice to separate the 
! location detail for the vertical and horizontal when possible.
! The internal representation of the location is currently implemented
! as radians from 0 to 2 PI for longitude and -PI/2 to PI/2 for latitude to
! minimize computational cost for distances. However, the external 
! representation is longitude in degrees from 0 to 360 and latitude 
! from -90 to 90 for consistency with most applications in the field.
! Note that for now, lev = -1 represents a surface quantity independent
! of vertical discretization as required for Bgrid surface pressure.

use      types_mod, only : r8, PI, RAD2DEG, DEG2RAD, MISSING_R8, MISSING_I
use  utilities_mod, only : register_module, error_handler, E_ERR, E_MSG,    &
                           logfileunit, nmlfileunit, find_namelist_in_file, &
                           check_namelist_read, do_output, do_nml_file,     &
                           do_nml_term, is_longitude_between
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none
private

public :: location_type, get_location, set_location,                          &
          set_location2, set_location_missing, is_location_in_region,         &
          write_location, read_location, interactive_location, vert_is_undef, &
          vert_is_surface, vert_is_pressure, vert_is_level, vert_is_height,   &
          query_location, LocationDims, LocationName, LocationLName,          &
          horiz_dist_only, get_close_obs, get_close_type,                     &
          get_close_maxdist_init, get_close_obs_init, get_close_obs_destroy,  &
          operator(==), operator(/=), VERTISUNDEF, VERTISSURFACE,             &
          VERTISLEVEL, VERTISPRESSURE, VERTISHEIGHT, get_dist,                &
          print_get_close_type

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

! The possible values for the location_type%which_vert component.
! These are intended to be PRIVATE to this module. Do not make public.

integer, parameter :: VERTISUNDEF    = -2 ! has no vertical location (undefined)
integer, parameter :: VERTISSURFACE  = -1 ! surface value
integer, parameter :: VERTISLEVEL    =  1 ! by level
integer, parameter :: VERTISPRESSURE =  2 ! by pressure
integer, parameter :: VERTISHEIGHT   =  3 ! by height

type location_type
   private 
   real(r8) :: lon, lat, vloc ! lon, lat stored in radians
   integer  :: which_vert     ! determines if by level, height, pressure, ...
end type location_type

! Type to facilitate efficient compuation of observations close to a given location
type get_close_type
   integer          :: num
   real(r8)         :: maxdist 
   integer, pointer :: lon_offset(:, :)     ! (nlat, nlat); 
   integer, pointer :: obs_box(:)           ! (nobs); List of obs indices in boxes
   integer, pointer :: count(:, :)          ! (nlon, nlat); # of obs in each box
   integer, pointer :: start(:, :)          ! (nlon, nlat); Start of list of obs in this box
   real(r8)         :: bot_lat, top_lat     ! Bottom and top latitudes of latitude boxes
   real(r8)         :: bot_lon, top_lon     ! Bottom and top longitudes of longitude boxes
   real(r8)         :: lon_width, lat_width ! Width of boxes in lon and lat
   logical          :: lon_cyclic           ! Do boxes wraparound in longitude?
end type get_close_type

type(random_seq_type) :: ran_seq
logical :: ran_seq_init = .false.
logical,save :: module_initialized = .false.

! Global storage for vertical distance normalization factors
real(r8) :: vert_normalization(3)

integer,              parameter :: LocationDims = 3
character(len = 129), parameter :: LocationName = "loc3Dsphere"
character(len = 129), parameter :: LocationLName = &
                                   "threed sphere locations: lon, lat, lev or pressure"

! Global storage for fast approximate sin and cosine lookups
! PAR For efficiency for small cases might want to fill tables as needed
real(r8) :: my_sin(-630:630), my_cos(-630:630), my_acos(-1000:1000)

! If maxdist stays the same, don't need to do box distance calculations
integer :: last_maxdist = -1.0

character(len=129) :: errstring

! Option for verification using exhaustive search
logical :: COMPARE_TO_CORRECT = .false.

!-----------------------------------------------------------------
! Namelist with default values
! horiz_dist_only == .true.  ->  Only the great circle horizontal distance is
!                                computed in get_dist.
! horiz_dist_only == .false. ->  Square root of sum of squared horizontal and
!                                normalized vertical distance computed in get_dist
! vert_normalization_pressure    Number of pascals that give a distance equivalent
!                                to one radian in horizontal
! vert_normalization_height ->   Number of meters that give a distance equivalent 
!                                to one radian in horizontal
! vert_normalization_level ->    Number of levels that give a distance equivalent
!                                to one radian in horizontal
! nlon ->                        Number of longitude boxes for get_close_obs MUST BE ODD
! nlat ->                        Number of latitude boxes for get_close_obs

logical  :: horiz_dist_only             = .true.
real(r8) :: vert_normalization_pressure = 100000.0_r8
real(r8) :: vert_normalization_height   = 10000.0_r8
real(r8) :: vert_normalization_level    = 20.0_r8 
logical  :: approximate_distance        = .false.
integer  :: nlon                        = 71
integer  :: nlat                        = 36
logical  :: output_box_info             = .false.
! should be a namelist item at some point
integer  :: print_box_level             = 0

namelist /location_nml/ horiz_dist_only, vert_normalization_pressure,         &
   vert_normalization_height, vert_normalization_level, approximate_distance, &
   nlon, nlat, output_box_info
!-----------------------------------------------------------------

interface operator(==); module procedure loc_eq; end interface
interface operator(/=); module procedure loc_ne; end interface

contains



subroutine initialize_module
!----------------------------------------------------------------------------
! subroutine initialize_module

integer :: iunit, io, i
character(len=129) :: str1,str2

call register_module(source, revision, revdate)
module_initialized = .true.

! Read the namelist entry
call find_namelist_in_file("input.nml", "location_nml", iunit)
read(iunit, nml = location_nml, iostat = io)
call check_namelist_read(iunit, io, "location_nml")

! Write the namelist values to the log file

if(do_nml_file()) write(nmlfileunit, nml=location_nml)
if(do_nml_term()) write(     *     , nml=location_nml)

! Make sure that the number of longitudes, nlon, for get_close_obs is odd
if(nlon / 2 * 2 == nlon) then
   call error_handler(E_ERR, 'initialize_module', 'nlon must be odd', &
      source, revision, revdate)
endif

! Copy the normalization factors in the vertical into an array
vert_normalization(1) = vert_normalization_level
vert_normalization(2) = vert_normalization_pressure
vert_normalization(3) = vert_normalization_height

write(str1,*) 'vert-normalization ', vert_normalization
write(str2,*) 'horizontal only ', horiz_dist_only
call error_handler(E_MSG,'location_mod:initialize_module',str1,source,revision,revdate)
call error_handler(E_MSG,'location_mod:initialize_module',str2,source,revision,revdate)

! Set up a lookup table for cos and sin for approximate but fast distances
! Don't worry about rounding errors as long as one gives more weight
! Really only need tables half this size, too (sin from -pi/2 to pi/2, cos only +)
if(approximate_distance) then
   do i = -630, 630
      my_cos(i) = cos(i / 100.0_r8)
      my_sin(i) = sin(i / 100.0_r8)
   end do
   do i = -1000, 1000
      my_acos(i) = acos(i / 1000.0_r8)
   end do
endif

end subroutine initialize_module



function get_dist(loc1, loc2, kind1, kind2, no_vert)
!----------------------------------------------------------------------------

implicit none

! Distance depends on only horizontal distance or both horizontal and vertical distance. 
! The choice is determined by horiz_dist_only and the which_vert of loc1.
! May want to allow return of some joint distance in the long run? 
! Or just a distance that is a function of all 3 things.
! The namelist controls whether default computations use just horizontal distance.
! However, this behavior can be over-ridden by the no_vert optional argument.
! If set to true, this will always do full 3d distance if possible. If set to
! false it will never do the full 3d distance. At present asking to do a vertical
! distance computation for incompatible vertical location types results in a fatal
! error.

! The kinds are available if a more sophisticated distance computation is required

type(location_type), intent(in) :: loc1, loc2
integer,             intent(in) :: kind1, kind2
real(r8)                        :: get_dist
logical, optional,   intent(in)  :: no_vert

real(r8) :: lon_dif, vert_dist, rtemp
integer  :: lat1_ind, lat2_ind, lon_ind, temp  ! indexes into lookup tables
logical  :: comp_h_only


if ( .not. module_initialized ) call initialize_module

! Begin with the horizontal distance
! Compute great circle path shortest route between two points
lon_dif = loc1%lon - loc2%lon

if(approximate_distance) then
   ! Option 1: Use table lookup; faster but less accurate
   lat1_ind = int(loc1%lat*100.0_r8)
   lat2_ind = int(loc2%lat*100.0_r8)
   lon_ind = int(lon_dif*100.0_r8)
   temp     = int(1000.0_r8 * (my_sin(lat2_ind) * my_sin(lat1_ind) + &
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

! Now compute a vertical distance if requested, 
! Highest priority is optional no_vert argument, test it first
if(present(no_vert)) then
   comp_h_only = no_vert
! Namelist horizontal only has second highest priority
else 
   comp_h_only = horiz_dist_only
endif
! Finally, if which_vert has no vertical definition for either location do only horizontal
if(loc1%which_vert == VERTISUNDEF .or. loc2%which_vert == VERTISUNDEF) comp_h_only = .true.

! Add in vertical component if required
if(.not. comp_h_only) then
   ! Vert distance can only be done for like vertical locations types
   if(loc1%which_vert /= loc2%which_vert) then
      write(errstring,*)'loc1%which_vert (',loc1%which_vert, &
                   ') /= loc2%which_vert (',loc2%which_vert,')'
      call error_handler(E_MSG, 'get_dist', errstring, source, revision, revdate)
      call write_location(logfileunit,loc1)
      call write_location(logfileunit,loc2)
      call write_location(0,loc1)
      call write_location(0,loc2)
      call error_handler(E_ERR, 'get_dist', &
         'Dont know how to compute vertical distance for unlike vertical coordinates', &
         source, revision, revdate)
   endif

   ! Compute the difference and divide by the appropriate normalization factor
   ! Normalization factor computes relative distance in vertical compared to one radian
   vert_dist = abs(loc1%vloc - loc2%vloc) / vert_normalization(loc1%which_vert)

   ! Spherical distance shape is computed here, other flavors can be computed
   get_dist = sqrt(get_dist**2 + vert_dist**2)
endif

end function get_dist




function get_location(loc)
!---------------------------------------------------------------------------
!
! Given a location type (in radians), 
! return the longitude, latitude (in degrees) and level 

implicit none

type(location_type), intent(in) :: loc
real(r8), dimension(3) :: get_location

   if ( .not. module_initialized ) call initialize_module

   get_location(1) = loc%lon * RAD2DEG                 
   get_location(2) = loc%lat * RAD2DEG                 
   get_location(3) = loc%vloc     

end function get_location



function vert_is_undef(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is undefined, else false

logical :: vert_is_undef
type(location_type), intent(in) :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == VERTISUNDEF) then
   vert_is_undef = .true.
else
   vert_is_undef = .false.
endif

end function vert_is_undef



function vert_is_surface(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is surface, else false

logical :: vert_is_surface
type(location_type), intent(in) :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == VERTISSURFACE) then
   vert_is_surface = .true.
else
   vert_is_surface = .false.
endif

end function vert_is_surface



function vert_is_pressure(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is pressure, else false

logical :: vert_is_pressure
type(location_type), intent(in) :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == VERTISPRESSURE) then
   vert_is_pressure = .true.
else
   vert_is_pressure = .false.
endif

end function vert_is_pressure



function vert_is_height(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is height, else false

logical :: vert_is_height
type(location_type), intent(in) :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == VERTISHEIGHT ) then
   vert_is_height = .true.
else
   vert_is_height = .false.
endif

end function vert_is_height



function vert_is_level(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is level, else false

logical :: vert_is_level
type(location_type), intent(in) :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == VERTISLEVEL) then
   vert_is_level = .true.
else
   vert_is_level = .false.
endif

end function vert_is_level



function loc_eq(loc1,loc2)
!---------------------------------------------------------------------------
!
! interface operator used to compare two locations.
! Returns true only if all components are 'the same' to within machine
! precision. There is some debate as to whether or not the vertical
! locations need to be identical if 'VERTISUNDEF' ... hard to peruse
! the code tree to find where this may be affected. 

implicit none

type(location_type), intent(in) :: loc1, loc2
logical :: loc_eq

if ( .not. module_initialized ) call initialize_module

loc_eq = .false.

! if ( loc1%which_vert /= loc2%which_vert ) return
if ( abs(loc1%lon  - loc2%lon ) > epsilon(loc1%lon ) ) return
if ( abs(loc1%lat  - loc2%lat ) > epsilon(loc1%lat ) ) return

!if ( loc1%which_vert /= VERTISUNDEF ) then
   if ( abs(loc1%vloc - loc2%vloc) > epsilon(loc1%vloc) ) return
!endif

loc_eq = .true.

end function loc_eq



function loc_ne(loc1,loc2)
!---------------------------------------------------------------------------
!
! interface operator used to compare two locations.
! Returns true if locations are not identical to machine precision.

implicit none

type(location_type), intent(in) :: loc1, loc2
logical :: loc_ne

if ( .not. module_initialized ) call initialize_module

loc_ne = (.not. loc_eq(loc1,loc2))

end function loc_ne



function get_location_lon(loc)
!---------------------------------------------------------------------------
!
! Given a location type, return the longitude.  Values stored in radians but
! returned in degrees.

implicit none

type(location_type), intent(in) :: loc
real(r8) :: get_location_lon

if ( .not. module_initialized ) call initialize_module

get_location_lon = loc%lon * RAD2DEG    

end function get_location_lon

function get_location_lat(loc)
!---------------------------------------------------------------------------
!
! Given a location type, return the latitude.  Values stored in radians but
! returned in degrees.

implicit none

type(location_type), intent(in) :: loc
real(r8) :: get_location_lat

if ( .not. module_initialized ) call initialize_module

get_location_lat = loc%lat * RAD2DEG      

end function get_location_lat



function set_location(lon, lat, vert_loc,  which_vert)
!----------------------------------------------------------------------------
!
! Puts the given longitude, latitude, and vertical location
! into a location datatype.  Arguments to this function are in degrees,
! but the values are stored as radians.
!

implicit none

type (location_type) :: set_location
real(r8), intent(in) :: lon, lat
real(r8), intent(in) :: vert_loc
integer,  intent(in) :: which_vert

if ( .not. module_initialized ) call initialize_module

if(lon < 0.0_r8 .or. lon > 360.0_r8) then
   write(errstring,*)'longitude (',lon,') is not within range [0,360]'
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
endif

if(lat < -90.0_r8 .or. lat > 90.0_r8) then
   write(errstring,*)'latitude (',lat,') is not within range [-90,90]'
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
endif

set_location%lon = lon * DEG2RAD
set_location%lat = lat * DEG2RAD

if(which_vert < VERTISUNDEF .or. which_vert == 0 .or. which_vert > VERTISHEIGHT) then
   write(errstring,*)'which_vert (',which_vert,') must be one of -2, -1, 1, 2, or 3'
   call error_handler(E_ERR,'set_location', errstring, source, revision, revdate)
endif

set_location%which_vert = which_vert
set_location%vloc = vert_loc

end function set_location


function set_location2(list)
!----------------------------------------------------------------------------
!
! location semi-independent interface routine
! given 4 float numbers, call the underlying set_location routine

implicit none

type (location_type) :: set_location2
real(r8), intent(in) :: list(:)

if ( .not. module_initialized ) call initialize_module

if (size(list) /= 4) then
   write(errstring,*)'requires 4 input values'
   call error_handler(E_ERR, 'set_location2', errstring, source, revision, revdate)
endif

set_location2 = set_location(list(1), list(2), list(3), nint(list(4)))

end function set_location2



function set_location_missing()
!----------------------------------------------------------------------------
!

type (location_type) :: set_location_missing

set_location_missing%lon        = missing_r8
set_location_missing%lat        = missing_r8
set_location_missing%vloc       = missing_r8
set_location_missing%which_vert = missing_i

end function set_location_missing



function query_location(loc,attr)
!---------------------------------------------------------------------------
!
! Returns the value of the requested location quantity
!

implicit none

type(location_type),        intent(in) :: loc
character(len=*), optional, intent(in) :: attr
real(r8)                               :: query_location

if ( .not. module_initialized ) call initialize_module

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

query_location = loc%which_vert

if (.not. present(attr)) return

selectcase(attr)
   case ('lat','LAT')
      query_location = loc%lat
   case ('lon','LON')
      query_location = loc%lon
   case ('vloc','VLOC')
      query_location = loc%vloc
end select

end function query_location


subroutine write_location(ifile, loc, fform)
!----------------------------------------------------------------------------
!
! Writes location to the file. Implemented as a subroutine but  could
! rewrite as a function with error control info returned. For initial implementation,
! file is just an integer file unit number. Probably want to replace this with file
! as a file_type allowing more flexibility for IO at later point. file_type and 
! associated operations would have to be supported. The mpp_io interfaces are a good
! place to head with this, perhaps, when we need to extend to supporting parallel
! platforms.

implicit none

integer,                    intent(in) :: ifile
type(location_type),        intent(in) :: loc
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile)loc%lon, loc%lat, loc%vloc, loc%which_vert
   CASE DEFAULT
      write(ifile, '(''loc3d'')' ) 
! Write out pressure or level along with integer tag
!     write(ifile, *)loc%lon, loc%lat, loc%vloc, loc%which_vert
! TJH -- would like to use a fixed-format write instead.:
!     the free-format write may/may not "wrap" depending on the precision
!     you compile with -- makes reading unpredictable ...
      write(ifile, '(1x,3(f22.14,1x),i8)')loc%lon, loc%lat, loc%vloc, loc%which_vert

END SELECT

end subroutine write_location


function read_location(ifile, fform)
!----------------------------------------------------------------------------
!
! Reads location from file that was written by write_location.
! See write_location for additional discussion.

implicit none

integer, intent(in) :: ifile
type(location_type) :: read_location
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile)read_location%lon, read_location%lat, &
         read_location%vloc, read_location%which_vert
   CASE DEFAULT
      read(ifile, '(a5)' ) header

      if(header /= 'loc3d') then
         write(errstring,*)'Expected location header "loc3d" in input file, got ', header
         call error_handler(E_ERR, 'read_location', errstring, source, revision, revdate)
      endif
! Now read the location data value
      read(ifile, *)read_location%lon, read_location%lat, &
         read_location%vloc, read_location%which_vert
END SELECT


end function read_location

subroutine interactive_location(location, set_to_default)
!--------------------------------------------------------------------------
!
! Allows for interactive input of a location. Also gives option of selecting
! a uniformly distributed random location (what the heck).

implicit none

type(location_type), intent(out) :: location
logical, intent(in), optional :: set_to_default

real(r8) :: lon, lat, minlon, maxlon, minlat, maxlat

if ( .not. module_initialized ) call initialize_module

! If set_to_default is true, then just zero out and return
if(present(set_to_default)) then
   if(set_to_default) then
      location%lon = 0.0
      location%lat = 0.0
      location%vloc = 0.0
      location%which_vert = 0
      return
   endif
endif

write(*, *)'Vertical co-ordinate options'
write(*, *)VERTISUNDEF,' --> vertical coordinate undefined'
write(*, *)VERTISSURFACE,' --> surface'
write(*, *)VERTISLEVEL,' --> model level'
write(*, *)VERTISPRESSURE,' --> pressure'
write(*, *)VERTISHEIGHT,' --> height'

100   read(*, *) location%which_vert
if(location%which_vert == VERTISLEVEL ) then
   write(*, *) 'Vertical co-ordinate model level'
   read(*, *) location%vloc
else if(location%which_vert == VERTISPRESSURE ) then
   write(*, *) 'Vertical co-ordinate Pressure (in hPa)'
   read(*, *) location%vloc
   location%vloc = 100.0 * location%vloc
else if(location%which_vert == VERTISHEIGHT ) then
   write(*, *) 'Vertical co-ordinate height (in gpm)'
   read(*, *) location%vloc
else if(location%which_vert == VERTISSURFACE ) then
   ! my understanding is that most of our users want height as
   ! the actual data value if the type of vert is surface.  but
   ! the original code was asking for pressure.  until i get a
   ! consensus, i will leave this alone. but it could be easily
   ! changed by commenting out the first and fourth lines below,
   ! and commenting in the second line (third line remains):
   write(*, *) 'Vertical co-ordinate surface pressure (in hPa)'
   !write(*, *) 'Vertical co-ordinate height (in gpm)'
   read(*, *) location%vloc
   location%vloc = 100.0 * location%vloc  ! only applies to pressure
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
!   location%lon = random_uniform(ran_seq) * 2.0_r8 * PI
   location%lon = random_uniform(ran_seq) * (maxlon-minlon) + minlon

   write(*, *) 'Input minimum latitude (-90.0 to 90.0)'
   read(*, *) minlat
   minlat = sin(minlat * DEG2RAD)

   write(*, *) 'Input maximum latitude (-90.0 to 90.0)'
   read(*, *) maxlat
   maxlat = sin(maxlat * DEG2RAD)

   ! Latitude must be area weighted
!   location%lat = asin(random_uniform(ran_seq) * 2.0_r8 - 1.0_r8)
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



subroutine nc_write_location(ncFileID, LocationVarID, loc, start)
!----------------------------------------------------------------------------
!
! Writes a SINGLE location to the specified netCDF variable and file.
!

use typeSizes
use netcdf

implicit none

integer, intent(in)             :: ncFileID, LocationVarID
type(location_type), intent(in) :: loc
integer, intent(in)             :: start

if ( .not. module_initialized ) call initialize_module

call check(nf90_put_var(ncFileID, LocationVarID, loc%lon,  (/ start, 1 /) ))
call check(nf90_put_var(ncFileID, LocationVarID, loc%lat,  (/ start, 2 /) ))
call check(nf90_put_var(ncFileID, LocationVarID, loc%vloc, (/ start, 3 /) ))

contains
  
  ! Internal subroutine - checks error status after each netcdf, prints
  !                       text message each time an error code is returned.
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) call error_handler(E_ERR,'nc_write_location', &
         trim(nf90_strerror(status)), source, revision, revdate )

  end subroutine check

end subroutine nc_write_location

!----------------------------------------------------------------------------

subroutine get_close_obs_init(gc, num, obs)

! Initializes part of get_close accelerator that depends on the particular obs

implicit none

type(get_close_type), intent(inout) :: gc            
integer,              intent(in)    :: num
type(location_type),  intent(in)    :: obs(num)

integer :: blat_ind, tlat_ind
integer :: bot_tlat_ind, top_tlat_ind
real(r8) :: base_lat(2), target_lat(2), del_lon, cos_del_lon
real(r8) :: max_del_lon
integer :: i, j, cum_start, lon_box(num), lat_box(num), tstart(nlon, nlat), tj, bj

if ( .not. module_initialized ) call initialize_module

! Allocate storage for obs number dependent part
allocate(gc%obs_box(num))
gc%obs_box(:) = -1

! Set the value of num_obs in the structure
gc%num = num

! If num == 0, no point in going any further.
if (num == 0) return

! Determine where the boxes should be for this set of obs and maxdist
call find_box_ranges(gc, obs, num)

! Figure out which boxes are close to a box on a given latitude circle
! MIGHT AVOID DOING THIS WITH A COPY ROUTINE: HAVE SAME BOXES OFTEN
do blat_ind = 1, nlat
   ! Search from east side of base block
   ! Start searching out, have to look for closest point in box being checked
   ! Only have to search latitude boxes that are within maximum distance
   bot_tlat_ind = blat_ind - floor(gc%maxdist / gc%lat_width ) - 1
   if(bot_tlat_ind < 1) bot_tlat_ind = 1
   top_tlat_ind = blat_ind + floor(gc%maxdist / gc%lat_width) + 1
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
         base_lat(j) = gc%bot_lat + (blat_ind - 2 + j) * gc%lat_width
         target_lat(j) = gc%bot_lat + (tlat_ind - 2 + j) * gc%lat_width
      end do

      ! If the max distance > PI, then everybody is close.
      ! Do a test for something slightly less than PI to avoid round-off error.
      ! Set max_del_lon to something much larger than PI since it doesn't matter.
      if(gc%maxdist > PI - 0.0001_r8) then
         max_del_lon = 2.0 * PI
      else
         ! Find the maximum longitude offset for the different possible latitudes edges
         max_del_lon = 0.0_r8
         do tj = 1, 2
            do bj = 1, 2
               ! Compute the lon offset directly by inverting distance
               cos_del_lon = (cos(gc%maxdist) - sin(base_lat(bj)) * sin(target_lat(tj))) / &
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
      gc%lon_offset(blat_ind, tlat_ind) = floor(max_del_lon / gc%lon_width) + 1
      ! Watch for roundoff leading to a search of more offsets than exist
      if(gc%lon_offset(blat_ind, tlat_ind) > nlon / 2) &
         gc%lon_offset(blat_ind, tlat_ind) = nlon / 2

   end do
end do

! Begin by computing the number of observations in each box in lat/lon
gc%count = 0
do i = 1, num
   lon_box(i) = get_lon_box(gc, obs(i)%lon)
   if(lon_box(i) < 0 .or. lon_box(i) > nlon) then
      write(errstring, *) 'Contact Dart Developers: this error should not happen'
      call error_handler(E_MSG, 'get_close_obs_init', errstring, source, revision, revdate)
      write(errstring, *) 'obs outside grid boxes, index value:',  lon_box(i)
      call error_handler(E_ERR, 'get_close_obs_init', errstring, source, revision, revdate)
   endif

   lat_box(i) = floor((obs(i)%lat - gc%bot_lat) / gc%lat_width) + 1
   if(lat_box(i) > nlat) lat_box(i) = nlat
   if(lat_box(i) < 1) lat_box(i) = 1

   gc%count(lon_box(i), lat_box(i)) = gc%count(lon_box(i), lat_box(i)) + 1
end do

! Figure out where storage for each boxes members should begin
cum_start = 1
do i = 1, nlon
   do j = 1, nlat
      gc%start(i, j) = cum_start
      cum_start = cum_start + gc%count(i, j)
   end do
end do

! Now we know how many are in each box, get a list of which are in each box
tstart = gc%start
do i = 1, num
   gc%obs_box(tstart(lon_box(i), lat_box(i))) = i
   tstart(lon_box(i), lat_box(i)) = tstart(lon_box(i), lat_box(i)) + 1
end do

! info on how well the boxes are working.  by default print nothing.
! set print_box_level to higher values to get more and more detail.
! user info should be level 1; 2 and 3 should be for debug only.
if (output_box_info .and. do_output()) call print_get_close_type(gc, print_box_level)

end subroutine get_close_obs_init



!----------------------------------------------------------------------------

subroutine get_close_obs_destroy(gc)

type(get_close_type), intent(inout) :: gc

deallocate(gc%obs_box, gc%lon_offset, gc%count, gc%start)

end subroutine get_close_obs_destroy

!----------------------------------------------------------------------------

subroutine get_close_maxdist_init(gc, maxdist)

implicit none

type(get_close_type), intent(inout) :: gc
real(r8),             intent(in)    :: maxdist

! Allocate the storage for the grid dependent boxes
allocate(gc%lon_offset(nlat, nlat), gc%count(nlon, nlat), gc%start(nlon, nlat))
gc%lon_offset = -1
gc%count      = -1
gc%start      = -1

! Set the maximum distance in the structure
gc%maxdist = maxdist

end subroutine get_close_maxdist_init

!----------------------------------------------------------------------------

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs, obs_kind, &
   num_close, close_ind, dist)

!!!ADD IN SOMETHING TO USE EFFICIENTLY IF IT"S AT SAME LOCATION AS PREVIOUS OB!!!

! The kinds are available to do more sophisticated distance computations if needed

implicit none

type(get_close_type),             intent(in)  :: gc
type(location_type),              intent(in)  :: base_obs_loc, obs(:)
integer,                          intent(in)  :: base_obs_kind, obs_kind(:)
integer,                          intent(out) :: num_close, close_ind(:)
real(r8),             optional,   intent(out) :: dist(:)

! If dist is NOT present, just find everybody in a box, put them in the list,
! but don't compute any distances

integer :: lon_box, lat_box, i, j, k, n_lon, lon_ind, n_in_box, st, t_ind
real(r8) :: this_dist

! Variables needed for comparing against correct case
integer :: cnum_close, cclose_ind(size(obs))
real(r8) :: cdist(size(obs))

! First, set the intent out arguments to a missing value
num_close = 0
close_ind = -99
if(present(dist)) dist = -99.0_r8
this_dist = 999999.0_r8   ! something big.

! If num == 0, no point in going any further. 
if (gc%num == 0) return


!--------------------------------------------------------------
! For validation, it is useful to be able to compare against exact
! exhaustive search
if(COMPARE_TO_CORRECT) then
   cnum_close = 0
   do i = 1, gc%num
      this_dist = get_dist(base_obs_loc, obs(i), base_obs_kind, obs_kind(i))
      if(this_dist < gc%maxdist) then
         ! Add this obs to correct list
         cnum_close = cnum_close + 1
         cclose_ind(cnum_close) = i
         cdist(cnum_close) = this_dist
      endif
   end do
endif

!--------------------------------------------------------------

! Begin by figuring out which box the base_ob is in
lon_box = get_lon_box(gc, base_obs_loc%lon)
lat_box = floor((base_obs_loc%lat - gc%bot_lat) / gc%lat_width) + 1

! If it is not in any box, then it is more than the maxdist away from everybody
if(lat_box > nlat .or. lat_box < 1 .or. lon_box < 0) return

! Next, loop through to find each box that is close to this box
do j = 1, nlat
   n_lon = gc%lon_offset(lat_box, j)
   if(n_lon >= 0) then
      LON_OFFSET: do i = -1 * n_lon, n_lon
         lon_ind = lon_box + i
         ! Search a box at this latitude j offset in longitude by i
         ! If domain is cyclic, need to wrap around
         if(gc%lon_cyclic) then
            if(lon_ind > nlon) lon_ind = lon_ind - nlon
            if(lon_ind < 1) lon_ind = lon_ind + nlon
         else
            ! Domain is not cyclic, don't search if outside of range
            if(lon_ind > nlon .or. lon_ind < 1) cycle LON_OFFSET
         endif
         ! Box to search is lon_ind, j
         n_in_box = gc%count(lon_ind, j)
         st = gc%start(lon_ind, j)
         ! Loop to check how close all obs in the box are; add those that are close
         do k = 1, n_in_box

            ! SHOULD ADD IN OPTIONAL ARGUMENT FOR DOING THIS!!!
            ! Could avoid adding any that have nums lower than base_ob???
            t_ind = gc%obs_box(st - 1 + k)
            ! Can compute total distance here if verts are the same
            ! Only compute distance if dist is present
            if(present(dist)) then
               if(base_obs_loc%which_vert == obs(t_ind)%which_vert) then
                  this_dist = get_dist(base_obs_loc, obs(t_ind), base_obs_kind, obs_kind(t_ind))
               else
               ! Otherwise can just get horizontal distance
                  this_dist = get_dist(base_obs_loc, obs(t_ind), base_obs_kind, obs_kind(t_ind), &
                     no_vert = .true.)
               endif
            else
               ! Dist isn't present; add this ob to list without computing distance
               num_close = num_close + 1
               close_ind(num_close) = t_ind
            endif

            ! If dist is present and this obs' distance is less than cutoff, add it in list
            if(present(dist) .and. this_dist <= gc%maxdist) then
               num_close = num_close + 1
               close_ind(num_close) = t_ind
               dist(num_close) = this_dist
            endif
         end do
      end do LON_OFFSET
   endif
end do

!------------------------ Verify by comparing to exhaustive search --------------
if(COMPARE_TO_CORRECT) then
   ! Do comparisons against full search
   if(num_close /= cnum_close) then
      write(*, *) 'ERROR: num_close, cnum_close', num_close, cnum_close
      stop
   endif
endif
!--------------------End of verify by comparing to exhaustive search --------------


end subroutine get_close_obs



subroutine find_box_ranges(gc, obs, num)
!--------------------------------------------------------------------------
!
! Finds boundaries for boxes in N/S direction. If data is localized in N/S
! tries to find boxes that only span the range of the data.
! 
type(get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
type(location_type),  intent(in)    :: obs(num)

real(r8) :: min_lat, max_lat, beg_box_lon, end_box_lon, first_obs_lon, last_obs_lon
real(r8) :: longitude_range, degrees
integer  :: i, indx, gap_start, gap_end, gap_length
logical  :: lon_box_full(360)

! Initialize boxes used to see where observations are
lon_box_full = .false.

! Figure out domain over which an additional obs MIGHT be close to one in this set
min_lat = minval(obs(:)%lat) - gc%maxdist
max_lat = maxval(obs(:)%lat) + gc%maxdist
if(min_lat < -PI / 2.0_r8) min_lat = -PI / 2.0_r8
if(max_lat > PI / 2.0_r8) max_lat = PI / 2.0_r8

! Put this into storage for this get_close_type
gc%bot_lat = min_lat
gc%top_lat  = max_lat
gc%lat_width = (max_lat - min_lat) / nlat
if(COMPARE_TO_CORRECT) write(*, *) 'min and max lat and width', gc%bot_lat, gc%top_lat, gc%lat_width

! Finding the longitude range is tricky because of cyclic nature
! Want to find minimum range spanned by obs even if they wrap-around Greenwich
! Would like to do this without sorting if possible at low-cost
! First, partition into 360 1-degree boxes and find the biggest gap
do i = 1, num
   degrees = obs(i)%lon * 180.0_r8 / PI
   ! If the value of the longitude is very close to an integer number of degrees
   ! a roundoff can occur that leads to an assignment in the wrong box.  We avoid this
   ! by first testing to see if this is possible and then setting both boxes to full.
   ! If this is not the case, then we fill the box the observation is in.
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
if(gap_length > 0) then
   ! There is a gap; figure out obs that are closest to ends of non-gap
   beg_box_lon = (gap_end / 180.0_r8) * PI
   end_box_lon = ((gap_start -1) / 180.0_r8) * PI
   first_obs_lon = find_closest_to_start(beg_box_lon, obs, num)
   last_obs_lon  = find_closest_to_end  (end_box_lon, obs, num)
   ! Determine the final longitude range
   longitude_range = last_obs_lon - first_obs_lon
   if(longitude_range <= 0.0_r8) longitude_range = longitude_range + 2.0_r8 * PI
   
   ! Add on the extra distance needed for the boxes

   ! To avoid any hard thinking about wraparound with sub-domain boxes
   ! Must span less than 180 degrees to get smaller boxes
   ! If addition of halos for close obs fills more than half of space things go 0 to 2PI
   if(longitude_range + 2.0_r8 * gc%maxdist > PI) then
      first_obs_lon = 0.0_r8
      last_obs_lon  = 2.0_r8 * PI
      gc%lon_cyclic = .true.
   else
      first_obs_lon = first_obs_lon - gc%maxdist 
      if(first_obs_lon < 0.0_r8) first_obs_lon = first_obs_lon + 2.0_r8 * PI
      last_obs_lon   = last_obs_lon + gc%maxdist 
      if(last_obs_lon > 2.0_r8 * PI) last_obs_lon = last_obs_lon - 2.0_r8 * PI
      gc%lon_cyclic = .false.
   endif
else
   ! No gap was found: all 360 boxes had an observation in them
   first_obs_lon = 0.0_r8
   last_obs_lon  = 2.0_r8 * PI
   gc%lon_cyclic = .true.
endif

! Put in storage for structure
gc%bot_lon = first_obs_lon
gc%top_lon = last_obs_lon
longitude_range = last_obs_lon - first_obs_lon
if(longitude_range <= 0.0_r8) longitude_range = longitude_range + 2.0_r8 * PI
gc%lon_width = longitude_range / nlon
if(COMPARE_TO_CORRECT) write(*, *) 'lon bot, top, width ', gc%bot_lon, gc%top_lon, gc%lon_width

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
next_full = next_full_box(start_box, lon_box_full, num_boxes)
! Find the next empty after that, make it the start of the gap
gap_start = next_empty_box(next_full, lon_box_full, num_boxes)
! Find the next full, box before that is the end of the gap
gap_end = next_full_box(gap_start, lon_box_full, num_boxes) - 1
if(gap_end < 1) gap_end = gap_end + num_boxes
! Carefully compute gap length
if(gap_end >= gap_start) then
   gap_length = gap_end - gap_start + 1
else
   gap_length = gap_end - gap_start + 360 + 1
endif

end subroutine find_next_gap


!----------------------------------------------------------------------------

function next_full_box(start_box, lon_box_full, num_boxes)

integer             :: next_full_box
integer, intent(in) :: start_box, num_boxes
logical, intent(in) :: lon_box_full(num_boxes)

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
call error_handler(E_ERR, 'next_full_box', 'All boxes empty:should not happen', &
   source, revision, revdate)

end function next_full_box

!----------------------------------------------------------------------------

function next_empty_box(start_box, lon_box_full, num_boxes)

integer             :: next_empty_box
integer, intent(in) :: start_box, num_boxes
logical, intent(in) :: lon_box_full(num_boxes)

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
call error_handler(E_ERR, 'next_empty_box', 'All boxes full:should not happen', &
   source, revision, revdate)

end function next_empty_box

!----------------------------------------------------------------------------

function find_closest_to_start(beg_box_lon, obs, num)

real(r8)                        :: find_closest_to_start
real(r8),            intent(in) :: beg_box_lon
integer,             intent(in) :: num
type(location_type), intent(in) :: obs(num)

real(r8) :: least_dist, dist
integer  :: i

! Start with large value
least_dist = 2.0_r8 * PI

do i = 1, num
   dist = obs(i)%lon - beg_box_lon
   if(dist < 0.0_r8) dist = dist + 2.0_r8 * PI
   if(dist < least_dist) then
      least_dist = dist
      find_closest_to_start = obs(i)%lon
   endif 
end do

end function find_closest_to_start

!----------------------------------------------------------------------------

function find_closest_to_end(end_box_lon, obs, num)

real(r8)                        :: find_closest_to_end
real(r8),            intent(in) :: end_box_lon
integer,             intent(in) :: num
type(location_type), intent(in) :: obs(num)

real(r8) :: least_dist, dist
integer  :: i

! Start with large value
least_dist = 2.0_r8 * PI

do i = 1, num
   dist = end_box_lon - obs(i)%lon
   if(dist < 0.0_r8) dist = dist + 2.0_r8 * PI
   if(dist < least_dist) then
      least_dist = dist
      find_closest_to_end = obs(i)%lon
   endif 
end do

end function find_closest_to_end


!----------------------------------------------------------------------------

function get_lon_box(gc, lon)

integer                          :: get_lon_box
type(get_close_type), intent(in) :: gc
real(r8),             intent(in) :: lon

real(r8) :: del_lon

del_lon = lon - gc%bot_lon
if(del_lon < 0.0_r8) del_lon = del_lon + 2.0_r8 * PI
get_lon_box = floor(del_lon / gc%lon_width) + 1
! On wraparound, correct for truncation
! If not wraparound, then we're not in one of the boxes
if(get_lon_box > nlon) then
   if(gc%lon_cyclic) then
      get_lon_box = 1
   else
      if(get_lon_box == nlon+1) then
         get_lon_box = nlon
      else
         get_lon_box = -1
      endif
   endif
endif

end function get_lon_box


!----------------------------------------------------------------------------

subroutine print_get_close_type(gc, amount)

type(get_close_type), intent(in) :: gc
integer, intent(in), optional    :: amount

integer :: i, j, k, first, index
integer :: sample, nfull, nempty, howmuch, total, maxcount, maxi, maxj
logical :: tickmark(gc%num)

! second arg is now an int, not logical, and means:
! 0 = very terse, only box summary (default).  
! 1 = structs and first part of arrays.
! 2 = all parts of all arrays.

! by default do not print all the obs_box or start contents (it can
! be very long).  but give an option, which if true, forces an
! entire contents dump.  'sample' is the number to print for
! the short version.  (this value prints about 5-6 lines of data.)
! to get a full dump, change howmuch to 2 below.
howmuch = 0
sample = 39

if (present(amount)) then
   howmuch = amount
endif

! print the get_close_type derived type values

if (howmuch > 0) then
   write(*,*) 'get_close_type values are:'

   write(*,*) ' num = ', gc%num
   write(*,*) ' maxdist = ', gc%maxdist
   write(*,*) ' bot_lat, top_lat = ', gc%bot_lat, gc%top_lat
   write(*,*) ' bot_lon, top_lon = ', gc%bot_lon, gc%top_lon
   write(*,*) ' lon_width, lat_width = ', gc%lon_width, gc%lat_width
   write(*,*) ' lon_cyclic = ', gc%lon_cyclic
endif

! this one can be very large.   print only the first nth unless
! instructed otherwise.  (print n+1 because 1 more value fits on
! the line because it prints ( i ) and not ( i, j ) like the others.)
if (associated(gc%obs_box)) then
   i = size(gc%obs_box,1)
   if (i/= gc%num) then
      write(*,*) ' warning: size of obs_box incorrect, nobs, i =', gc%num, i
   endif
   if (howmuch > 1) then
      write(*,*) ' obs_box(',i,') =', gc%obs_box    ! (nobs)
   else if(howmuch > 0) then
      write(*,*) ' obs_box(',i,') =', gc%obs_box(1:min(i,sample+1))    ! (nobs)
      write(*,*) '  <rest of obs_box omitted>'
   endif
else
   if (howmuch > 0) write(*,*) ' obs_box unallocated'
endif

! like obs_box, this one can be very large.   print only the first nth unless
! instructed otherwise
if (associated(gc%start)) then
   i = size(gc%start,1)
   j = size(gc%start,2)
   if ((i /= nlon) .or. (j /= nlat)) then
      write(*,*) ' warning: size of start incorrect, nlon, nlat, i, j =', nlon, nlat, i, j
   endif
   if (howmuch > 1) then
      write(*,*) ' start(',i,j,') =', gc%start    ! (nlon, nlat)
   else if (howmuch > 0) then
      write(*,*) ' start(',i,j,') =', gc%start(1:min(i,sample), 1)    ! (nlon, nlat)
      write(*,*) '  <rest of start omitted>'
   endif
else
   if (howmuch > 0) write(*,*) ' start unallocated'
endif

! as above, print only first n unless second arg is .true.
if (associated(gc%lon_offset)) then
   i = size(gc%lon_offset,1)
   j = size(gc%lon_offset,2)
   if ((i /= nlat) .or. (j /= nlat)) then
      write(*,*) ' warning: size of lon_offset incorrect, nlat, i, j =', nlat, i, j
   endif
   if (howmuch > 1) then
      write(*,*) ' lon_offset(',i,j,') =', gc%lon_offset    ! (nlon, nlat)
   else if (howmuch > 0) then
      write(*,*) ' lon_offset(',i,j,') =', gc%lon_offset(1:min(i,sample), 1)    ! (nlon, nlat)
      write(*,*) '  <rest of lon_offset omitted>'
   endif
else
   if (howmuch > 0) write(*,*) ' lon_offset unallocated'
endif

! as above, print only first n unless second arg is .true.
if (associated(gc%count)) then
   i = size(gc%count,1)
   j = size(gc%count,2)
   if ((i /= nlon) .or. (j /= nlat)) then
      write(*,*) ' warning: size of count incorrect, nlon, nlat, i, j =', &
                      nlon, nlat, i, j
   endif
   if (howmuch > 1) then
      write(*,*) ' count(',i,j,') =', gc%count    ! (nlon, nlat)
   else if (howmuch > 0) then
      write(*,*) ' count(',i,j,') =', gc%count(1:min(i,sample), 1)    ! (nlon, nlat)
      write(*,*) '  <rest of count omitted>'
   endif
else
   if (howmuch > 0) write(*,*) ' count unallocated'
endif


! end of print.  strictly speaking the following code is validation code,
! not a simple print.


! initialize all ticks to false.  turn them true as they are found
! in the obs_box list, and complain about duplicates or misses.
tickmark = .FALSE.

do i=1, nlon
   do j=1, nlat
      first = gc%start(i, j)
      do k=1, gc%count(i, j)
         index = first + k - 1
         if ((index < 1) .or. (index > gc%num)) then
            write(*, *) 'error: bad obs list index, in box: ', index, i, j
            write(*, *) 'exiting without checking further'
            exit
         endif
         if (tickmark(index)) then
            write(*, *) 'error: obs found in more than one box list.  index, box: ', &
                         index, i, j
            write(*, *) 'exiting without checking further'
            exit
         endif
         tickmark(index) = .TRUE.
      enddo
   enddo
enddo

do i=1, gc%num
  if (.not. tickmark(i)) then
     write(*,*) 'error: obs found in no box list: ', i
     write(*, *) 'exiting without checking further'
     exit
  endif
enddo

! print out some hopefully useful stats
nfull = 0
nempty = 0
total = 0
maxcount = 0
maxi = 0
maxj = 0

do i=1, nlon
   do j=1, nlat
      if (gc%count(i, j) > 0) then
         nfull = nfull + 1
         total = total + gc%count(i, j)
         if (gc%count(i, j) > maxcount) then
            maxcount = gc%count(i, j)
            maxi = i
            maxj = j
         endif
      else
         nempty = nempty + 1
      endif
   enddo
enddo

! these print out always - make sure they are useful to end users.
write(*, '(a)') "Location module statistics:"
write(*, '(a,i9)') " Total boxes (nlon * nlat): ", nfull + nempty
write(*, '(a,i9)') " Total items to put in boxes: ", gc%num
if (howmuch > 0) then
   write(*, '(a,i9)') " Total boxes with 1+ items: ", nfull
write(*, '(a,i9)') " Total boxes empty: ", nempty
endif
if (nfull > 0) then
   write(*, '(a,f7.2)') " Percent boxes with 1+ items: ", nfull / real(nfull + nempty, r8) * 100.
   write(*, '(a,f12.2)') " Average #items per non-empty box: ", real(total, r8) / nfull
endif
if (maxcount > 0) then
   write(*, '(a,i9)') " Largest #items in one box: ", maxcount
! leave this out for now.  one, if there are multiple boxes with
! the same maxcount this is just the last one found.  two, the
! index numbers do not seem very helpful.
!   if (howmuch > 0) write(*, '(a,i9,i9)') " That box index: ", maxi, maxj
endif


end subroutine print_get_close_type


function is_location_in_region(loc, minl, maxl)
!----------------------------------------------------------------------------
!
! Returns true if the given location is inside the rectangular
! region defined by minl as the lower left, maxl the upper right.
! test is inclusive; values on the edges are considered inside.
! Periodic in longitude (box can cross the 2PI -> 0 line)

implicit none

logical                          :: is_location_in_region
type(location_type), intent(in)  :: loc, minl, maxl

if ( .not. module_initialized ) call initialize_module

! maybe could use VERTISUNDEF in the minl and maxl args to indicate
! we want to test only in horizontal?  and if not, vtypes must match?
!if ( (minl%which_vert /= maxl%which_vert) .or. &
! ((minl%which_vert /= loc%which_vert).and.(minl%which_vert /= VERTISUNDEF))) then
!   write(errstring,*)'which_vert (',loc%which_vert,') must be same in all args'
!   call error_handler(E_ERR, 'is_location_in_region', errstring, source, revision, revdate)
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
! end of location/threed_sphere/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod
