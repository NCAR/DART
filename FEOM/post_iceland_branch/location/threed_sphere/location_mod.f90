! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module location_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 
!

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
use  utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, logfileunit, &
                           find_namelist_in_file, check_namelist_read
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none
private

public :: location_type, get_dist, get_location, set_location, set_location_missing, &
          write_location, read_location, interactive_location, vert_is_undef, &
          vert_is_surface, vert_is_pressure, vert_is_level, vert_is_height, &
          query_location, LocationDims, LocationName, LocationLName, &
          horiz_dist_only, alloc_get_close_obs, get_close_obs, &
          operator(==), operator(/=), VERTISUNDEF, VERTISSURFACE, &
          VERTISLEVEL, VERTISPRESSURE, VERTISHEIGHT

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
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
   real(r8) :: lon, lat, vloc
   integer  :: which_vert     ! determines if by level, height, pressure, ...
end type location_type

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
real(r8) :: my_sin(-630:630), my_cos(-630:630), my_acos(-1000:1000)

! Global storage for efficient get_close_obs search
! WARNING: NLON MUST BE ODD
integer, parameter :: nlon = 71, nlat = 36
!integer, parameter :: nlon = 9, nlat = 10

! Permanent global storage for efficient get_close_obs
integer :: count(nlon, nlat), lon_offset(nlat, nlat)
integer :: cum_start, start(nlon, nlat)

! If cutoff stays the same, don't need to do box distance calculations
integer :: last_cutoff = -1.0


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

logical  :: horiz_dist_only             = .true.
real(r8) :: vert_normalization_pressure = 100000.0_r8
real(r8) :: vert_normalization_height   = 10000.0_r8
real(r8) :: vert_normalization_level    = 20.0_r8 
logical  :: approximate_distance        = .false.

namelist /location_nml/ horiz_dist_only, vert_normalization_pressure, &
   vert_normalization_height, vert_normalization_level, approximate_distance
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

call error_handler(E_MSG,'initialize_module','location_nml values',' ',' ',' ')
write(logfileunit, nml=location_nml)
write(     *     , nml=location_nml)

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



function get_dist(loc1, loc2, no_vert)
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

type(location_type), intent(in) :: loc1, loc2
real(r8) :: get_dist
logical, optional :: no_vert

real(r8) :: lon_dif, vert_dist
integer  :: lat1_ind, lat2_ind, lon_ind, temp  ! indexes into lookup tables
logical  :: comp_h_only

character(len=129) :: errstring

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
      get_dist = acos(sin(loc2%lat) * sin(loc1%lat) + &
         cos(loc2%lat) * cos(loc1%lat) * cos(lon_dif))
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
! precision.

implicit none

type(location_type), intent(in) :: loc1, loc2
logical :: loc_eq

if ( .not. module_initialized ) call initialize_module

loc_eq = .false.

if ( abs(loc1%lon  - loc2%lon ) > epsilon(loc1%lon ) ) return
if ( abs(loc1%lat  - loc2%lat ) > epsilon(loc1%lat ) ) return
if ( abs(loc1%vloc - loc2%vloc) > epsilon(loc1%vloc) ) return

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
! Given a location type, return the longitude

implicit none

type(location_type), intent(in) :: loc
real(r8) :: get_location_lon

if ( .not. module_initialized ) call initialize_module

get_location_lon = loc%lon * RAD2DEG    

end function get_location_lon

function get_location_lat(loc)
!---------------------------------------------------------------------------
!
! Given a location type, return the latitude

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
! into a location datatype.
!

implicit none

type (location_type) :: set_location
real(r8), intent(in) :: lon, lat
real(r8), intent(in) :: vert_loc
integer,  intent(in) :: which_vert

character(len=129) :: errstring

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



function set_location_missing()
!----------------------------------------------------------------------------
!

type (location_type) :: set_location_missing

set_location_missing%lon        = missing_r8
set_location_missing%lat        = missing_r8
set_location_missing%vloc       = missing_r8
set_location_missing%which_vert = missing_i

end function set_location_missing



function query_location(loc,attr) result(fval)
!---------------------------------------------------------------------------
!
! Returns the value of the attribute
!

implicit none

type(location_type),        intent(in) :: loc
character(len=*), optional, intent(in) :: attr
real(r8)                               :: fval

character(len=16) :: attribute

if ( .not. module_initialized ) call initialize_module

attribute = 'which_vert'

if (present(attr)) attribute = attr
selectcase(adjustl(attribute))
 case ('which_vert','WHICH_VERT')
   fval = loc%which_vert
 case ('lat','LAT')
   fval = loc%lat
 case ('lon','LON')
   fval = loc%lon
 case ('vloc','VLOC')
   fval = loc%vloc
 case default
   fval = loc%which_vert
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
character(len=129) :: errstring
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
   write(*, *) 'Vertical co-ordinate surface pressure (in hPa)'
   read(*, *) location%vloc
   location%vloc = 100.0 * location%vloc
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

subroutine alloc_get_close_obs(num, obs, cutoff, obs_box)

implicit none

integer, intent(in) :: num
type(location_type), intent(in) :: obs(num)
real(r8), intent(in) :: cutoff
integer, intent(out) :: obs_box(num)

integer :: i, j, blat_ind, tlat_ind
integer :: lon_box(num), lat_box(num), tstart(nlon, nlat)
integer :: bot_tlat_ind, top_tlat_ind
real(r8) :: base_lon, base_lat, target_lat, del_lon, cos_del_lon

! Begin by computing the number of observations in each box in lat/lon
count = 0
do i = 1, num
   lon_box(i) = int(nlon * obs(i)%lon / (2.0_r8 * PI)) + 1
   if(lon_box(i) > nlon) lon_box(i) = 1
   lat_box(i) = int(nlat * (obs(i)%lat + PI / 2.0_r8) / PI) + 1
   if(lat_box(i) > nlat) lat_box(i) = nlat
   if(lat_box(i) < 1) lat_box(i) = 1
   count(lon_box(i), lat_box(i)) = count(lon_box(i), lat_box(i)) + 1
end do

! Figure out where storage for each boxes members should begin
cum_start = 1
do i = 1, nlon
   do j = 1, nlat
      start(i, j) = cum_start
      cum_start = cum_start + count(i, j)
   end do
end do

tstart = start
! Now we know how many are in each box, get a list of which are in each box
do i = 1, num
   obs_box(tstart(lon_box(i), lat_box(i))) = i
   tstart(lon_box(i), lat_box(i)) = tstart(lon_box(i), lat_box(i)) + 1
end do

! Don't need to do the close box computation again if cutoff hasn't changed
if(cutoff == last_cutoff) then
   return
else
   last_cutoff = cutoff
endif

! Figure out which boxes are close to a box on a given latitude circle
do blat_ind = 1, nlat
   ! Search from east side of base block
   base_lon = 360.0_r8 /nlon
   ! Start searching out, have to look for closest point in box being checked
   ! Only have to search latitude boxes that are within cutoff distance
   bot_tlat_ind = blat_ind - int(cutoff * nlat / PI) - 1
   if(bot_tlat_ind < 1) bot_tlat_ind = 1
   top_tlat_ind = blat_ind + int(cutoff * nlat / PI) + 1
   if(top_tlat_ind > nlat) top_tlat_ind = nlat
   do tlat_ind = bot_tlat_ind, top_tlat_ind
      ! For lats north of me, search from NE corner as base to SW corner as target
      if(tlat_ind > blat_ind) then
         base_lat = -90.0_r8 + blat_ind * 180.0_r8 / nlat
         target_lat = -90.0_r8 + (tlat_ind - 1) * 180.0_r8 / nlat
      else if(tlat_ind < blat_ind) then
         ! Otherwise, want to search from SE corner to NW corner
         base_lat = -90.0_r8 + (blat_ind - 1) * 180.0_r8 / nlat
         target_lat = -90.0_r8 + tlat_ind * 180.0_r8 / nlat
      else
         ! When same latitude, do both from southern corner
         base_lat = -90.0_r8 + (blat_ind - 1) * 180.0_r8 / nlat
         target_lat = base_lat
      endif

      ! Compute the lon offset directly by inverting distance
      cos_del_lon = (cos(cutoff) - sin(base_lat*DEG2RAD) * sin(target_lat*DEG2RAD)) / &
         (cos(base_lat*DEG2RAD) * cos(target_lat*DEG2RAD))
      if(cos_del_lon < -1) then
         del_lon = PI
      else if(cos_del_lon > 1) then
         del_lon = 0
      else
         del_lon = acos(cos_del_lon)
      endif
      lon_offset(blat_ind, tlat_ind) = int(del_lon * nlon / (2.0 * PI)) + 1

   end do
end do

end subroutine alloc_get_close_obs


!----------------------------------------------------------------------------

subroutine get_close_obs(base_ob, num, obs, cutoff, obs_box, num_close, close_ind, dist)

!!!ADD IN SOMETHING TO USE EFFICIENTLY IF IT"S AT SAME LOCATION AS PREVIOUS OB!!!

implicit none

integer, intent(in) :: base_ob, num
type(location_type), intent(in) :: obs(num)
real(r8), intent(in) :: cutoff
integer, intent(in) :: obs_box(num)
integer, intent(out) :: num_close, close_ind(num)
real(r8), intent(out) :: dist(num)

integer :: lon_box, lat_box, i, j, k, n_lon, lon_ind, n_in_box, st, t_ind
real(r8) :: this_dist

! Begin by figuring out which box the base_ob is in
lon_box = int(nlon * obs(base_ob)%lon / (2.0_r8 * PI)) + 1
if(lon_box > nlon) lon_box = 1
lat_box = int(nlat * (obs(base_ob)%lat + PI / 2.0_r8) / PI) + 1
if(lat_box > nlat) lat_box = nlat

num_close = 0

! Next, loop through to find each box that is close to this box
do j = 1, nlat
   n_lon = lon_offset(lat_box, j)
   if(n_lon > 0) then
      do i = -1 * (n_lon -1), n_lon - 1
         ! Search a box at this latitude j offset in longitude by i
         lon_ind = lon_box + i
         if(lon_ind > nlon) lon_ind = lon_ind - nlon
         if(lon_ind < 1) lon_ind = lon_ind + nlon
         ! Box to search is lon_ind, j
         n_in_box = count(lon_ind, j)
         st = start(lon_ind, j)
         ! Loop to check how close all obs in the box are; add those that are close
         do k = 1, n_in_box
            ! Could avoid adding any that have nums lower than base_ob???
            t_ind = obs_box(st - 1 + k)
            ! Can compute total distance here if verts are the same
            if(obs(base_ob)%which_vert == obs(t_ind)%which_vert) then
               this_dist = get_dist(obs(base_ob), obs(t_ind))
            else
            ! Otherwise can just get horizontal distance
               this_dist = get_dist(obs(base_ob), obs(t_ind), no_vert = .true.)
            endif
            if(this_dist <= cutoff) then
               ! Add this ob to the list
               num_close = num_close + 1
               close_ind(num_close) = t_ind
               dist(num_close) = this_dist
            endif
         end do
      end do
   endif
end do

end subroutine get_close_obs


!----------------------------------------------------------------------------
! end of location/threed_sphere/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod
