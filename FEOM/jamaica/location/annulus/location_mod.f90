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

! Implements location interfaces for a three dimensional annulus
! with a vertical coordinate based on the models native set of
! discrete levels. The internal representation of the location is 
! currently implemented as radians from 0 to 2 PI for the azimuthal
! direction (longitude-like).  The radial distance is latitude-like,
! and the vertical coordinate is zero at the bottom of the annulus.
!
! Note: in an effort to maintain consistency with the other 
! location_mods, will use lon for azimuthal direction, lat for
! radial direction, and lev for depth.

use      types_mod, only : r8, PI, RAD2DEG, DEG2RAD, MISSING_R8, MISSING_I
use  utilities_mod, only : register_module, error_handler, E_ERR
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none
private

public :: location_type, get_dist, get_location, set_location, set_location_missing, &
          write_location, read_location, interactive_location, &
          vert_is_pressure, vert_is_level, vert_is_height, query_location, &
          LocationDims, LocationName, LocationLName, get_close_obs, alloc_get_close_obs, &
          operator(==), operator(/=)

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type location_type
   private 
   real(r8) :: lon, lat, vloc
   integer  :: which_vert
   ! which_vert determines if the location is by level or by height/pressure
   ! 1 ===> obs is by level
   ! 2 ===> obs is by pressure (Not supported in the annulus)
   ! 3 ===> obs is by height
end type location_type

type(random_seq_type) :: ran_seq
logical :: ran_seq_init = .false.
logical,save :: module_initialized = .false.

integer,              parameter :: LocationDims = 3
character(len = 129), parameter :: LocationName = "loc_annulus"
character(len = 129), parameter :: LocationLName = &
                                   "Annulus location: azimuthal angle, radius, and level"

interface operator(==); module procedure loc_eq; end interface
interface operator(/=); module procedure loc_ne; end interface

contains



  subroutine initialize_module
!----------------------------------------------------------------------------
! subroutine initialize_module

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module



function get_dist(loc1, loc2)
!----------------------------------------------------------------------------

implicit none

! To comply with current (Guam) configuration, distance depends only on 
! horizontal distance.
!
! Note that for compatibility with other location modules, lat represents
! the radial direction, and lon represents the azimuthal angle.

type(location_type), intent(in) :: loc1, loc2
real(r8) :: get_dist

real(r8) :: x1, y1, x2, y2

if ( .not. module_initialized ) call initialize_module

! convert from cylindrical to cartesian coordinates
x1 = loc1%lat*cos(loc1%lon)
y1 = loc1%lat*sin(loc1%lon)
x2 = loc2%lat*cos(loc2%lon)
y2 = loc2%lat*sin(loc2%lon)

! Returns distance in m
get_dist = sqrt((x1 - x2)**2 + (y1 - y2)**2)

end function get_dist



function get_location(loc)
!---------------------------------------------------------------------------
!
! Given a location type (where the azimuthal angle is in radians), this 
! routine return the azimuthal angle in degrees, the radius, and the level 

implicit none

type(location_type), intent(in) :: loc
real(r8), dimension(3) :: get_location

   if ( .not. module_initialized ) call initialize_module

   get_location(1) = loc%lon * RAD2DEG                 
   get_location(2) = loc%lat
   get_location(3) = loc%vloc     

end function get_location



function vert_is_pressure(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is pressure, else false.
!
! Always returns false, as vertical coordinate is never pressure for the annulus.

logical :: vert_is_pressure
type(location_type), intent(in) :: loc

if ( .not. module_initialized ) call initialize_module

vert_is_pressure = .false.

end function vert_is_pressure



function vert_is_height(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is height, else false.

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
! Given a location, return true if vertical coordinate is pressure, else false.

logical :: vert_is_level
type(location_type), intent(in) :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == 1) then
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
! Given a longitude location, return the azimuthal angle in degrees

implicit none

type(location_type), intent(in) :: loc
real(r8) :: get_location_lon

if ( .not. module_initialized ) call initialize_module

get_location_lon = loc%lon * RAD2DEG    

end function get_location_lon



function get_location_lat(loc)
!---------------------------------------------------------------------------
!
! Given a latitude location, return the radius

implicit none

type(location_type), intent(in) :: loc
real(r8) :: get_location_lat

if ( .not. module_initialized ) call initialize_module

get_location_lat = loc%lat

end function get_location_lat



function set_location(lon, lat, vert_loc, which_vert)
!----------------------------------------------------------------------------
!
! Puts the given azimuthal angle (in degrees), radius, and vertical 
! location into a location datatype.

implicit none

type (location_type) :: set_location
real(r8), intent(in) :: lon, lat
real(r8), intent(in) :: vert_loc
integer,  intent(in) :: which_vert

character(len=129) :: errstring

if ( .not. module_initialized ) call initialize_module

if(lon < 0.0_r8 .or. lon > 360.0_r8) then
   write(errstring,*)'azimuthal angle (',lon,') is not within range [0,360]'
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
endif

set_location%lon = lon * DEG2RAD
set_location%lat = lat 

if(which_vert < 1 .or. which_vert == 2 .or. which_vert > 3  ) then
   write(errstring,*)'which_vert (',which_vert,') must be either 1 or 3'
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
! Writes location to the file. Implemented as a subroutine but could
! rewrite as a function with error control info returned. For initial 
! implementation, file is just an integer file unit number. 

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

integer, intent(in)                    :: ifile
type(location_type)                    :: read_location
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

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
      location%lon        = 0.0_r8
      location%lat        = 0.0_r8
      location%vloc       = 0.0_r8
      location%which_vert = 0
      return
   endif
endif

write(*, *)'Vertical co-ordinate options'
write(*, *)'-1 -> surface, 1 -> model level, 3 -> depth'

100   read(*, *) location%which_vert
if(location%which_vert == 1 ) then
   write(*, *) 'Vertical co-ordinate model level'
   read(*, *) location%vloc
else if(location%which_vert == 3 ) then
   write(*, *) 'Vertical co-ordinate depth (in negative m)'
   read(*, *) location%vloc
   do while (location%vloc > 0)
      write(*, *) 'Depth must be negative (zero at top of fluid), please try again'
      read(*, *) location%vloc
   end do
else if(location%which_vert == -1 ) then
   write(*, *) 'Vertical co-ordinate surface pressure (in hPa)'
   read(*, *) location%vloc
   location%vloc = 100.0 * location%vloc
else
   write(*, *) 'Wrong choice of which_vert try again between -1, 1, and 3'
   go to 100
end if

write(*, *) 'Input azimuthal angle: value 0 to 360.0 or a negative number '
write(*, *) 'for uniformly distributed random location in the horizontal.'
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

   write(*, *) 'Input minimum azimuthal angle (0 to 360.0)'
   read(*, *) minlon
   minlon = minlon * DEG2RAD

   write(*, *) 'Input maximum azimuthal angle(0 to 360.0)'
   read(*, *) maxlon
   maxlon = maxlon * DEG2RAD

   ! Longitude is random from minlon to maxlon
   location%lon = random_uniform(ran_seq) * (maxlon-minlon) + minlon

   write(*, *) 'Input minimum radius '
   read(*, *) minlat

   write(*, *) 'Input maximum radius '
   read(*, *) maxlat

   ! Latitude must be area weighted to obtain proper random realizations
   location%lat = sqrt(random_uniform(ran_seq) * (maxlat-minlat) + minlat)

   write(*, *) 'random location is ', location%lon / DEG2RAD, &
                                      location%lat 

else

   write(*, *) 'Input radius '
   read(*, *) lat

   do while(lat < 0)
      write(*, *) 'Input value is illegal, please try again'
      read(*, *) lat
   end do

   location%lat = lat
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

! This does pre-computing for close obs; no function needed in one dimension

return

end subroutine alloc_get_close_obs


!----------------------------------------------------------------------------

subroutine get_close_obs(base_ob, num, obs, cutoff, obs_box, num_close, close_ind, dist)

! Default version with no smarts; no need to be smart in 1D

implicit none

integer, intent(in) :: base_ob, num
type(location_type), intent(in) :: obs(num)
real(r8), intent(in) :: cutoff
integer, intent(in) :: obs_box(num)
integer, intent(out) :: num_close, close_ind(num)
real(r8), intent(out) :: dist(num)

integer :: i
real(r8) :: this_dist

! Return list of obs that are within cutoff and their distances
num_close = 0
do i = 1, num
   this_dist = get_dist(obs(base_ob), obs(i))
   if(this_dist <= cutoff) then
      ! Add this ob to the list
      num_close = num_close + 1
      close_ind(num_close) = i
      dist(num_close) = this_dist 
   endif
end do

end subroutine get_close_obs



!----------------------------------------------------------------------------
! end of location/annulus/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod
