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

! Implements location interfaces for a two dimensional spherical shell 
! plus a naive vertical coordinate based on the models native set of
! discrete levels. In the long run, it would be nice to separate the 
! location detail for the vertical and horizontal when possible.
! The internal representation of the location is currently implemented
! as radians from 0 to 2 PI for longitude and -PI/2 to PI/2 for latitude to
! minimize computational cost for distances. However, the external 
! representation is longitude in degrees from 0 to 360 and latitude 
! from -90 to 90 for consistency with most applications in the field.
! Note that for now, lev = -1 represents a surface quantity independent
! of vertical discretization as required for Bgrid surface pressure.

use      types_mod, only : r8, PI
use  utilities_mod, only : register_module, error_handler, E_ERR
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none
private

public location_type, get_dist, get_location, set_location, &
       write_location, read_location, interactive_location, &
       vert_is_level, &
       LocationDims, LocationName, LocationLName, &
       read_ncep_obs_location

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type location_type
   private
   real(r8) :: lon, lat, lev, pressure
   integer :: which_vert
! which vert determines if location has a pressure or level location in vertical
! 1 implies levels for now, 2 implies pressure, could add other options
end type location_type

type(random_seq_type) :: ran_seq
logical :: ran_seq_init = .false.
logical, save :: module_initialized = .false.

integer,              parameter :: LocationDims = 3
character(len = 129), parameter :: LocationName = "loc3Dsphere"
character(len = 129), parameter :: LocationLName = &
                                   "simple threed sphere locations: lon, lat, lev"


contains


  subroutine initialize_module
!----------------------------------------------------------------------------
! subroutine initialize_module

   call register_module(source,revision,revdate)
   module_initialized = .true.

end subroutine initialize_module



function get_dist(loc1, loc2)
!----------------------------------------------------------------------------

implicit none

! For now, distance depends only on horizontal distance. May want to allow
! return of some joint distance in the long run?

type(location_type), intent(in) :: loc1, loc2
real(r8) :: get_dist

real(r8) :: lon_dif

if ( .not. module_initialized ) call initialize_module

! Returns distance in radians (independent of diameter of sphere)

! Compute great circle path shortest route between two points
lon_dif = abs(loc1%lon - loc2%lon)
if(lon_dif > PI) lon_dif = 2.0_r8 * PI - lon_dif

if(cos(loc1%lat) == 0) then
   get_dist = abs(loc2%lat - loc1%lat)
else
   get_dist = acos(sin(loc2%lat) * sin(loc1%lat) + &
                   cos(loc2%lat) * cos(loc1%lat) * cos(lon_dif))
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

get_location(1) = loc%lon * 180.0_r8 / PI
get_location(2) = loc%lat * 180.0_r8 / PI
get_location(3) = loc%lev

end function get_location



function vert_is_level(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is pressure, else false

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

get_location_lon = loc%lon * 180.0_r8 / PI

end function get_location_lon



function get_location_lat(loc)
!---------------------------------------------------------------------------
!
! Given a location type, return the latitude

implicit none

type(location_type), intent(in) :: loc
real(r8) :: get_location_lat

if ( .not. module_initialized ) call initialize_module

get_location_lat = loc%lat * 180.0_r8 / PI

end function get_location_lat



function get_location_lev(loc)
!---------------------------------------------------------------------------
!
! Given a location type, return the level

implicit none

type(location_type), intent(in) :: loc
real(r8) :: get_location_lev

if ( .not. module_initialized ) call initialize_module

get_location_lev = loc%lev

end function get_location_lev



function set_location(lon, lat, lev)
!----------------------------------------------------------------------------
!
! Given a longitude and latitude
! puts this value into the location.

implicit none

type (location_type) :: set_location
real(r8), intent(in) :: lon, lat, lev

if ( .not. module_initialized ) call initialize_module

if(lon < 0.0_r8 .or. lon > 360.0_r8) call error_handler(E_ERR, 'set_location', &
        'Longitude is out of [0,360] range', source, revision, revdate)

if(lat < -90.0_r8 .or. lat > 90.0_r8) call error_handler(E_ERR, 'set_location', &
        'Latitude is out of [-90,90] range', source, revision, revdate)

set_location%lon = lon * PI / 180.0_r8
set_location%lat = lat * PI / 180.0_r8
set_location%lev = lev

end function set_location



subroutine write_location(ifile, loc)
!----------------------------------------------------------------------------
!
! Writes a oned location to the file. Implemented as a subroutine but  could
! rewrite as a function with error control info returned. For initial implementation,
! file is just an integer file unit number. Probably want to replace this with file
! as a file_type allowing more flexibility for IO at later point. file_type and 
! associated operations would have to be supported. The mpp_io intefaces are a good
! place to head with this, perhaps, when we need to extend to supporting parallel
! platforms. 

implicit none

integer, intent(in) :: ifile
type(location_type), intent(in) :: loc

if ( .not. module_initialized ) call initialize_module

! For now, output a character tag followed by the r8 value. 

write(ifile, '(''loc2s'')' ) 
write(ifile, *) loc%lon, loc%lat, loc%lev

end subroutine write_location



function read_location(ifile)
!----------------------------------------------------------------------------
!
! Reads a oned location from ifile that was written by write_location. 
! See write_location for additional discussion.

implicit none

integer, intent(in) :: ifile
type(location_type) :: read_location

character(len=5) :: header

if ( .not. module_initialized ) call initialize_module

! Will want to add additional error checks on the read
read(ifile, '(a5)' ) header

if(header /= 'loc2s') call error_handler(E_ERR, 'read_location', &
    'Expected location header "loc2s" in input file', source, revision, revdate)

! Now read the location data value
read(ifile, *) read_location%lon, read_location%lat, read_location%lev

end function read_location



subroutine interactive_location(location)
!--------------------------------------------------------------------------
!
! Allows for interactive input of a location. Also gives option of selecting
! a uniformly distributed random location (what the heck).

implicit none

type(location_type), intent(out) :: location

real(r8) :: lon, lat, lev

if ( .not. module_initialized ) call initialize_module

write(*, *) 'Input level for this observation: -1 for surface '
read(*, *) lev

write(*, *) 'Input longitude for this obs: value 0 to 360.0 or a negative number for '
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

   ! Longitude is random from 0 to 2 PI
   location%lon = random_uniform(ran_seq) * 2.0_r8 * PI

   ! Latitude must be area weightedA
   location%lat = asin(random_uniform(ran_seq) * 2.0_r8 - 1.0_r8)
   location%lev = lev

   write(*, *) 'random location is ', location%lon * 180.0 / PI, &
                                      location%lat * 180.0 / PI

else

   write(*, *) 'Input latitude for this obs: value -90.0 to 90.0'
   read(*, *) lat
   do while(lat < -90.0_r8 .or. lat > 90.0_r8)
      write(*, *) 'Input value < -90.0 or > 90.0 is illegal, please try again'
      read(*, *) lat
   end do
   location = set_location(lon, lat, lev)
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

call check(nf90_put_var(ncFileID, LocationVarID, loc%lon, (/ start, 1 /) ))
call check(nf90_put_var(ncFileID, LocationVarID, loc%lat, (/ start, 2 /) ))
call check(nf90_put_var(ncFileID, LocationVarID, loc%lev, (/ start, 3 /) ))

contains
  
  ! Internal subroutine - checks error status after each netcdf, prints
  !                       text message each time an error code is returned.
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) call error_handler(E_ERR, 'nc_write_locations', &
         trim(nf90_strerror(status)),source, revision, revdate)
    
  end subroutine check

end subroutine nc_write_location



!=============================================================
  subroutine read_ncep_obs_location(location, obsunit, obsindex, var)
!=============================================================
! read location (lon,lat,pressure) from NCEP observation files
!  Input units are: radians and hPa.

implicit none

type(location_type) :: location

integer :: obs_prof
integer, intent(in) :: obsunit
integer, intent(out) :: obsindex
real (r8), intent(out) :: var
real (r8) :: lon,lat,lev,zob, dummy,cnt,rtime,rtype

if ( .not. module_initialized ) call initialize_module

! Read location, kind and error variance of NCEP data

    read(obsunit, 880) var, lon, lat, lev, zob, dummy,cnt,rtime,rtype
  880 format(f4.2, 2f7.3, f7.1, f7.2, f7.2, f9.0, f7.3, f5.0)

    location%lon = lon     ! in radian
    location%lat = lat     ! in radian

!   set up observation kind
    obs_prof = cnt/1000000

    if(obs_prof == 2) obsindex = 1
    if(obs_prof == 9) obsindex = 2
    if(obs_prof == 3) obsindex = 3
    if(obs_prof == 1) obsindex = 4

    if (obsindex .ne. 3) then      ! for u,v,t
    location%pressure = lev*100.0   ! (transfer from mb to Pascal)
    location%which_vert = 2
    else
    location%lev = -1      ! for Ps
    location%which_vert = 1
    var = var*100.0     ! convert to Pascal
    endif

end subroutine read_ncep_obs_location

!----------------------------------------------------------------------------
! end of location/simple_threed_sphere/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod
