! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module location_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
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
use  utilities_mod, only : file_exist, open_file, close_file, check_nml_error, &
                           register_module, error_handler, E_ERR, E_MSG, logfileunit
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none
private

public location_type, get_dist, get_location, set_location, set_location_missing, &
       write_location, read_location, interactive_location, &
       vert_is_pressure, vert_is_level, vert_is_height, query_location, &
       LocationDims, LocationName, LocationLName

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type location_type
   private 
   real(r8) :: lon, lat, vloc
   integer  :: which_vert
   ! which_vert determines if the location is by level or by height/pressure
   ! 1 ===> obs is by level
   ! 2 ===> obs is by pressure
   ! 3 ===> obs is by height

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
! vert_normalization_level ->    Number of levels that give a distnace equivalent
!                                to one radian in horizontal

logical  :: horiz_dist_only = .true.
real(r8) :: vert_normalization_pressure = 100000.0_r8
real(r8) :: vert_normalization_height   = 10000.0_r8
real(r8) :: vert_normalization_level    = 20.0_r8 

namelist /location_nml/ horiz_dist_only, vert_normalization_pressure, &
   vert_normalization_height, vert_normalization_level
!-----------------------------------------------------------------

contains



subroutine initialize_module
!----------------------------------------------------------------------------
! subroutine initialize_module

integer :: iunit, ierr, io

call register_module(source, revision, revdate)
module_initialized = .true.

! Read namelist for run time control

if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   ierr = 1

   READBLOCK: do while(ierr /= 0)
      read(iunit, nml = location_nml, iostat = io)
      if ( io < 0 ) exit READBLOCK          ! end-of-file
      ierr = check_nml_error(io, 'location_nml')
   enddo READBLOCK

   call close_file(iunit)
endif

! Write the namelist values to the log file

call error_handler(E_MSG,'initialize_module','location namelist values',' ',' ',' ')
write(logfileunit, nml=location_nml)

! Copy the normalization factors in the vertical into an array
vert_normalization(1) = vert_normalization_level
vert_normalization(2) = vert_normalization_pressure
vert_normalization(3) = vert_normalization_height

write(*, *) 'vert-normalization ', vert_normalization
write(*, *) 'horizontal only ', horiz_dist_only

end subroutine initialize_module



function get_dist(loc1, loc2)
!----------------------------------------------------------------------------

implicit none

! For now, distance depends only on horizontal distance. May want to allow
! return of some joint distance in the long run? Or just a distance that
! is a function of all 3 things.

type(location_type), intent(in) :: loc1, loc2
real(r8) :: get_dist

real(r8) :: lon_dif, horiz_dist, vert_dist

if ( .not. module_initialized ) call initialize_module

! Returns horizontal distance in radians (independent of diameter of sphere)

! Begin with the horizontal distance
! Compute great circle path shortest route between two points
lon_dif = abs(loc1%lon - loc2%lon)
if(lon_dif > PI) lon_dif = 2.0_r8 * PI - lon_dif

if(cos(loc1%lat) == 0) then
   horiz_dist = abs(loc2%lat - loc1%lat)
else
   horiz_dist = acos(sin(loc2%lat) * sin(loc1%lat) + &
      cos(loc2%lat) * cos(loc1%lat) * cos(lon_dif))
endif

! Now compute a vertical distance if requested, 
if(horiz_dist_only) then
   get_dist = horiz_dist
else
   ! Vert distance can only be done for like vertical locations types
   if(loc1%which_vert /= loc2%which_vert) then
      call error_handler(E_ERR, 'get_dist', &
         'Dont know how to compute vertical distance for unlike vertical coordinates', &
         source, revision, revdate)
   endif

   ! Compute the difference and divide by the appropriate normalization factor
   ! Normalization factor computes relative distance in vertical compared to one radian
   vert_dist = abs(loc1%vloc - loc2%vloc) / vert_normalization(loc1%which_vert)

   ! Spherical distance shape is computed here, other flavors can be computed
   get_dist = sqrt(horiz_dist**2 + vert_dist**2)
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



function vert_is_pressure(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is pressure, else false

logical :: vert_is_pressure
type(location_type), intent(in) :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == 2) then
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

if(loc%which_vert == 3 ) then
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

if(loc%which_vert == 1) then
   vert_is_level = .true.
else
   vert_is_level = .false.
endif

end function vert_is_level

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

if(which_vert < 1 .or. which_vert > 3  ) then
   write(errstring,*)'which_vert (',which_vert,') must be one of 1, 2, or 3'
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
      write(ifile, *)loc%lon, loc%lat, loc%vloc, loc%which_vert
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
write(*, *)'-1 -> surface, 1 -> model level, 2 -> pressure, 3 -> height'
100   read(*, *) location%which_vert
if(location%which_vert == 1 ) then
   write(*, *) 'Vertical co-ordinate model level'
   read(*, *) location%vloc
else if(location%which_vert == 2 ) then
   write(*, *) 'Vertical co-ordinate Pressure (in hPa)'
   read(*, *) location%vloc
   location%vloc = 100.0 * location%vloc
else if(location%which_vert == 3 ) then
   write(*, *) 'Vertical co-ordinate height (in gpm)'
   read(*, *) location%vloc
else if(location%which_vert == -1 ) then
   write(*, *) 'Vertical co-ordinate surface pressure (in hPa)'
   read(*, *) location%vloc
   location%vloc = 100.0 * location%vloc
else
   write(*, *) 'Wrong choice of which_vert try again between -1 and 3'
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
! end of location/threed_sphere/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod
