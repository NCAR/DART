module location_mod
!
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

use      types_mod
!use  utilities_mod, only : output_err, E_ERR
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

private

public location_type, get_dist, get_location, set_location, &
       write_location, read_location, interactive_location, &
       vert_is_pressure, vert_is_level, &
       LocationDims, LocationName, LocationLName, &
       read_ncep_obs_location

! let CVS fill strings ... DO NOT EDIT ...
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

! CVS Generated file description for error handling, do not edit
character(len = 129), parameter :: &
   e_src = "$Source$", &
   e_rev = "$Revision$", &
   e_dat = "$Date$", &
   e_aut = "$Author$"

! There needs to be some sort of public metadata for the location module.
! The number of dimensions, the name of each dimension, that sort of thing.
! TJH Sept. 16, 2002
!
!type location_meta
!   integer             :: ndims = 3
!   character (len=129) :: name = "loc3Dsphere"
!   character (len=129) :: longname = "threed sphere locations: lon, lat, lev"
!   character (len=129) :: rev = "$Revision$"
!   character (len=129) :: dat = "$Date$"
!   character (len=129) :: aut = "$Author$"
!   character (len=129) :: src = &
!   "$Source$"
!end type location_meta

integer,              parameter :: LocationDims = 3
character(len = 129), parameter :: LocationName = "loc3Dsphere"
character(len = 129), parameter :: LocationLName = &
                                   "threed sphere locations: lon, lat, lev or pressure"


contains


function get_dist(loc1, loc2)
!----------------------------------------------------------------------------

implicit none

! For now, distance depends only on horizontal distance. May want to allow
! return of some joint distance in the long run? Or just a distance that
! is a function of all 3 things.

type(location_type), intent(in) :: loc1, loc2
real(r8) :: get_dist

real(r8) :: lon_dif

! Returns distance in radians (independent of diameter of sphere)

! Compute great circle path shortest route between two points
lon_dif = abs(loc1%lon - loc2%lon)
if(lon_dif > pi) lon_dif = 2.0_r8 * pi - lon_dif

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

get_location(1) = loc%lon * 360.0_r8 / (2.0_r8 * pi)
get_location(2) = loc%lat * 180.0_r8 / pi
! Return level or pressure as indicated
if(loc%which_vert == 1) then
   get_location(3) = loc%lev
else
   get_location(3) = loc%pressure
end if

end function get_location



function vert_is_pressure(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is pressure, else false

logical :: vert_is_pressure
type(location_type), intent(in) :: loc

if(loc%which_vert == 2) then
   vert_is_pressure = .true.
else
   vert_is_pressure = .false.
endif

end function vert_is_pressure



function vert_is_level(loc)
!---------------------------------------------------------------------------
!
! Given a location, return true if vertical coordinate is pressure, else false

logical :: vert_is_level
type(location_type), intent(in) :: loc

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

get_location_lon = loc%lon * 360.0_r8 / (2.0_r8 * pi)

end function get_location_lon



function get_location_lat(loc)
!---------------------------------------------------------------------------
!
! Given a location type, return the latitude

implicit none

type(location_type), intent(in) :: loc
real(r8) :: get_location_lat

get_location_lat = loc%lat * 180.0_r8 / pi

end function get_location_lat



function set_location(lon, lat, lev, pressure)
!----------------------------------------------------------------------------
!
! Given a longitude and latitude
! puts this value into the location.

implicit none

type (location_type) :: set_location
real(r8), intent(in) :: lon, lat
real(r8), intent(in), optional :: lev, pressure

!if(lon < 0.0_r8 .or. lon > 360.0_r8) call output_err(E_ERR, e_src, e_rev, e_dat, e_aut, &
!   'set_location', 'Longitude is out of 0->360 range')
if(lon < 0.0_r8 .or. lon > 360.0_r8) then
   write(*, *) 'set location: Longitude is out of 0->360 range'
   stop
endif

!if(lat < -90.0_r8 .or. lat > 90.0_r8) call output_err(E_ERR, e_src, e_rev, e_dat, e_aut, &
!   'set_location', 'Latitude is out of -90->90 range')
if(lat < -90.0_r8 .or. lat > 90.0_r8) then
   write(*, *) 'set_location  Latitude is out of -90->90 range'
   stop
endif


set_location%lon = lon * 2.0_r8 * pi / 360.0_r8
set_location%lat = lat * pi / 180.0_r8

! Set the vertical coordinate to pressure or level as required; watch for out of range
if(present(lev)) then
   if(lev < -1) then
      write(*, *) 'set_location: Level is less than -1'
      stop
   endif
   set_location%lev = lev
   set_location%which_vert = 1
else if(present(pressure)) then
   if(pressure < 0.0) then
      write(*, *) 'set_location: Pressure is less than 0'
      stop
   endif
   set_location%pressure = pressure
   set_location%which_vert = 2
endif

end function set_location



subroutine write_location(file, loc)
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

integer, intent(in) :: file
type(location_type), intent(in) :: loc

! For now, output a character tag followed by the r8 value. Is this written at 
! machine precision ???

write(file, 11) 
11 format('loc3d')       ! TJH request

! Write out pressure or level along with integer tag
if(loc%which_vert == 1) then
   write(file, *) loc%lon, loc%lat, loc%lev, loc%which_vert
else
   write(file, *) loc%lon, loc%lat, loc%pressure, loc%which_vert
endif

! I need to learn how to do formatted IO with the types package
!21 format(r8)

end subroutine write_location



function read_location(file)
!----------------------------------------------------------------------------
!
! Reads a oned location from file that was written by write_location. See write_location
! for additional discussion.

implicit none

integer, intent(in) :: file
type(location_type) :: read_location
real(r8) :: lev_or_press

character(len=5) :: header

! Will want to add additional error checks on the read
read(file, 11) header
11 format(a5)
!if(header /= 'loc2s') call output_err(E_ERR, e_src, e_rev, e_dat, e_aut, &
!   'read_location', 'Expected location header "loc1d" in input file')
if(header /= 'loc3d') then
   write(*, *) 'read_location Expected location header "loc3d" in input file'
   stop
endif

! Now read the location data value
read(file, *) read_location%lon, read_location%lat, lev_or_press, read_location%which_vert
if(read_location%which_vert == 1) then
   read_location%lev = lev_or_press
else
   read_location%pressure = lev_or_press
endif

end function read_location



subroutine interactive_location(location)
!--------------------------------------------------------------------------
!
! Allows for interactive input of a location. Also gives option of selecting
! a uniformly distributed random location (what the heck).

implicit none

type(location_type), intent(out) :: location

real(r8) :: lon, lat, lev, pressure

write(*, *) 'Input model level with -1 for surface or a -99 to input pressure'

read(*, *) lev

if(lev < -98.0) then
! Convert from Millibars to Pascals as required for model internals
   write(*, *) 'Input pressure for this observation in PASCALS'
   read(*, *) location%pressure
   location%which_vert = 2
else
   location%lev = lev
   location%which_vert = 1
endif

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

   write(*, *) 'random location is ', location%lon * 180.0 / PI, &
                                      location%lat * 180.0 / PI

else
   write(*, *) 'Input latitude for this obs: value -90.0 to 90.0'
   read(*, *) lat
   do while(lat < -90.0_r8 .or. lat > 90.0_r8)
      write(*, *) 'Input value < -90.0 or > 90.0 is illegal, please try again'
      read(*, *) lat
   end do
   location%lon = lon
   location%lat = lat
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

call check(nf90_put_var(ncFileID, LocationVarID, loc%lon, (/ start, 1 /) ))
call check(nf90_put_var(ncFileID, LocationVarID, loc%lat, (/ start, 2 /) ))
! ??? Need output for pressure locations, too ???
call check(nf90_put_var(ncFileID, LocationVarID, loc%lev, (/ start, 3 /) ))

contains
  
  ! Internal subroutine - checks error status after each netcdf, prints
  !                       text message each time an error code is returned.
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      print *,'location_mod:nc_write_location'
      stop
    end if
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
real (r8) :: lon,lat,lev,zob, dummy,count,time,type                                         
                                                                                            
! Read location, kind and error variance of NCEP data 

    read(obsunit, 880) var, lon, lat, lev, zob, dummy,count,time,type               
  880 format(f4.2, 2f7.3, f7.1, f7.2, f7.2, f9.0, f7.3, f5.0)

    location%lon = lon     ! in radian
    location%lat = lat     ! in radian
                                                                                            
!   set up observation kind
    obs_prof = count/1000000

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
!
!----------------------------------------------------------------------------
! end of location/threed_sphere/location_mod.f90
!----------------------------------------------------------------------------
!
end module location_mod
