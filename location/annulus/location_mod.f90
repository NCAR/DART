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
! (Note2: not clear this is a good thing - what consistency do we buy?
!  might want to use real names for the values for clarity.)

use      types_mod, only : r8, PI, RAD2DEG, DEG2RAD, MISSING_R8, MISSING_I
use  utilities_mod, only : register_module, error_handler, E_ERR, ascii_file_format
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none
private

public :: location_type, get_dist, get_location, set_location, &
          set_location2, set_location_missing, is_location_in_region, &
          write_location, read_location, interactive_location, &
          vert_is_pressure, vert_is_level, vert_is_height, query_location, &
          LocationDims, LocationName, LocationLName, &
          get_close_obs, alloc_get_close_obs, get_close_type, &
          get_close_maxdist_init, get_close_obs_init, get_close_obs_destroy, &
          operator(==), operator(/=)

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type location_type
   private 
   real(r8) :: azm, rad, vloc
   integer  :: which_vert
   ! which_vert determines if the location is by level or by height/pressure
   ! 1 ===> obs is by level
   ! 2 ===> obs is by pressure (Not supported in the annulus)
   ! 3 ===> obs is by height
end type location_type

type get_close_type
   private
   integer :: maxdist
end type get_close_type

type(random_seq_type) :: ran_seq
logical :: ran_seq_init = .false.
logical,save :: module_initialized = .false.

character(len=129) :: errstring

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

! To comply with current (Guam) configuration, distance depends only on 
! horizontal distance.
!

type(location_type), intent(in) :: loc1, loc2
real(r8) :: get_dist

real(r8) :: x1, y1, x2, y2

if ( .not. module_initialized ) call initialize_module

! convert from cylindrical to cartesian coordinates
x1 = loc1%rad*cos(loc1%azm)
y1 = loc1%rad*sin(loc1%azm)
x2 = loc2%rad*cos(loc2%azm)
y2 = loc2%rad*sin(loc2%azm)

! Returns distance in m
get_dist = sqrt((x1 - x2)**2 + (y1 - y2)**2)

end function get_dist



function get_location(loc)
!---------------------------------------------------------------------------
!
! Given a location type (where the azimuthal angle is in radians), this 
! routine return the azimuthal angle in degrees, the radius, and the level 

type(location_type), intent(in) :: loc
real(r8), dimension(3) :: get_location

   if ( .not. module_initialized ) call initialize_module

   get_location(1) = loc%azm * RAD2DEG                 
   get_location(2) = loc%rad
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

type(location_type), intent(in) :: loc1, loc2
logical :: loc_eq

if ( .not. module_initialized ) call initialize_module

loc_eq = .false.

if ( abs(loc1%azm  - loc2%azm ) > epsilon(loc1%azm ) ) return
if ( abs(loc1%rad  - loc2%rad ) > epsilon(loc1%rad ) ) return
if ( abs(loc1%vloc - loc2%vloc) > epsilon(loc1%vloc) ) return

loc_eq = .true.

end function loc_eq



function loc_ne(loc1,loc2)
!---------------------------------------------------------------------------
!
! interface operator used to compare two locations.
! Returns true if locations are not identical to machine precision.

type(location_type), intent(in) :: loc1, loc2
logical :: loc_ne

if ( .not. module_initialized ) call initialize_module

loc_ne = (.not. loc_eq(loc1,loc2))

end function loc_ne



function get_location_lon(loc)
!---------------------------------------------------------------------------
!
! Given a longitude location, return the azimuthal angle in degrees

type(location_type), intent(in) :: loc
real(r8) :: get_location_lon

if ( .not. module_initialized ) call initialize_module

get_location_lon = loc%azm * RAD2DEG    

end function get_location_lon



function get_location_lat(loc)
!---------------------------------------------------------------------------
!
! Given a latitude location, return the radius

type(location_type), intent(in) :: loc
real(r8) :: get_location_lat

if ( .not. module_initialized ) call initialize_module

get_location_lat = loc%rad

end function get_location_lat



function set_location(lon, lat, vert_loc, which_vert)
!----------------------------------------------------------------------------
!
! Puts the given azimuthal angle (in degrees), radius, and vertical 
! location into a location datatype.

type (location_type) :: set_location
real(r8), intent(in) :: lon, lat
real(r8), intent(in) :: vert_loc
integer,  intent(in) :: which_vert

if ( .not. module_initialized ) call initialize_module

if(lon < 0.0_r8 .or. lon > 360.0_r8) then
   write(errstring,*)'azimuthal angle (',lon,') is not within range [0,360]'
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
endif

set_location%azm = lon * DEG2RAD
set_location%rad = lat 

if(which_vert /= -1 .and. which_vert /= 1 .and. which_vert /= 3  ) then
   write(errstring,*)'which_vert (',which_vert,') must be -1, 1 or 3'
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

set_location_missing%azm        = missing_r8
set_location_missing%rad        = missing_r8
set_location_missing%vloc       = missing_r8
set_location_missing%which_vert = missing_i

end function set_location_missing



function query_location(loc,attr) result(fval)
!---------------------------------------------------------------------------
!
! Returns the value of the attribute
!

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
   fval = loc%rad
 case ('lon','LON')
   fval = loc%azm
 case ('vloc','VLOC')
   fval = loc%vloc
 case default
   fval = loc%which_vert
end select



end function query_location

subroutine write_location(ifile, loc, fform, charstring)
!----------------------------------------------------------------------------
!
! Writes location to the file. Implemented as a subroutine but could
! rewrite as a function with error control info returned. For initial 
! implementation, file is just an integer file unit number. 

integer,                     intent(in) :: ifile
type(location_type),         intent(in) :: loc
character(len=*), intent(in),  optional :: fform
character(len=*), intent(out), optional :: charstring

logical            :: writebuf
integer            :: charlength

if ( .not. module_initialized ) call initialize_module

! incoming char string must be long enough for format
writebuf = present(charstring)
if (writebuf) then
   charlength = 82
   if (len(charstring) < charlength) then
      write(errstring, *) 'charstring buffer must be at least ', charlength, ' chars long'
      call error_handler(E_ERR, 'write_location', errstring, source, revision, revdate)
   endif
endif


if (ascii_file_format(fform)) then
   if (writebuf) then
      write(charstring,'(3(f22.14,1x),i4)') &
            loc%azm, loc%rad, loc%vloc, loc%which_vert
   else
      ! Write out pressure or level along with integer tag
      write(ifile, '(''loc3d'')' ) 
      write(ifile, '(1x,3(f22.14,1x),i4)')loc%azm, loc%rad, loc%vloc, loc%which_vert
   endif
else
   if (writebuf) then
      call error_handler(E_ERR, 'write_location', &
           'Cannot use string buffer with binary format', &
            source, revision, revdate)
   endif

   write(ifile)loc%azm, loc%rad, loc%vloc, loc%which_vert
endif

end subroutine write_location



function read_location(ifile, fform)
!----------------------------------------------------------------------------
!
! Reads location from file that was written by write_location.
! See write_location for additional discussion.

integer, intent(in)                    :: ifile
type(location_type)                    :: read_location
character(len=*), intent(in), optional :: fform

character(len=5)   :: header

if ( .not. module_initialized ) call initialize_module

if (ascii_file_format(fform)) then
   read(ifile, '(a5)' ) header

   if(header /= 'loc3d') then
      write(errstring,*)'Expected location header "loc3d" in input file, got ', header
      call error_handler(E_ERR, 'read_location', errstring, source, revision, revdate)
   endif
   ! Now read the location data value
   read(ifile, *)read_location%azm, read_location%rad, &
                 read_location%vloc, read_location%which_vert
else
   read(ifile)read_location%azm, read_location%rad, &
              read_location%vloc, read_location%which_vert
endif

end function read_location



subroutine interactive_location(location, set_to_default)
!--------------------------------------------------------------------------
!
! Allows for interactive input of a location. Also gives option of selecting
! a uniformly distributed random location (what the heck).

type(location_type), intent(out) :: location
logical, intent(in), optional :: set_to_default

real(r8) :: azm, rad, minazm, maxazm, minrad, maxrad

if ( .not. module_initialized ) call initialize_module

! If set_to_default is true, then just zero out and return
if(present(set_to_default)) then
   if(set_to_default) then
      location%azm        = 0.0_r8
      location%rad        = 0.0_r8
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
read(*, *) azm

do while(azm > 360.0_r8)
   write(*, *) 'Input value greater than 360.0 is illegal, please try again'
   read(*, *) azm
end do

if(azm < 0.0_r8) then

   ! Need to make sure random sequence is initialized

   if(.not. ran_seq_init) then
      call init_random_seq(ran_seq)
      ran_seq_init = .TRUE.
   endif

   write(*, *) 'Input minimum azimuthal angle (0 to 360.0)'
   read(*, *) minazm
   minazm = minazm * DEG2RAD

   write(*, *) 'Input maximum azimuthal angle(0 to 360.0)'
   read(*, *) maxazm
   maxazm = maxazm * DEG2RAD

   ! Longitude is random from minlon to maxlon
   location%azm = random_uniform(ran_seq) * (maxazm-minazm) + minazm

   write(*, *) 'Input minimum radius '
   read(*, *) minrad

   write(*, *) 'Input maximum radius '
   read(*, *) maxrad

   ! Latitude must be area weighted to obtain proper random realizations
   location%rad = sqrt(random_uniform(ran_seq) * (maxrad-minrad) + minrad)

   write(*, *) 'random location is ', location%azm / DEG2RAD, &
                                      location%rad 

else

   write(*, *) 'Input radius '
   read(*, *) rad

   do while(rad < 0)
      write(*, *) 'Input value is illegal, please try again'
      read(*, *) rad
   end do

   location%rad = rad
   location%azm = azm*DEG2RAD

end if

end subroutine interactive_location



subroutine nc_write_location(ncFileID, LocationVarID, loc, start)
!----------------------------------------------------------------------------
!
! Writes a SINGLE location to the specified netCDF variable and file.
!

use typeSizes
use netcdf

integer, intent(in)             :: ncFileID, LocationVarID
type(location_type), intent(in) :: loc
integer, intent(in)             :: start

if ( .not. module_initialized ) call initialize_module

call check(nf90_put_var(ncFileID, LocationVarID, loc%azm,  (/ start, 1 /) ))
call check(nf90_put_var(ncFileID, LocationVarID, loc%rad,  (/ start, 2 /) ))
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

subroutine get_close_obs_destroy(gc)

type(get_close_type), intent(inout) :: gc


end subroutine get_close_obs_destroy

!----------------------------------------------------------------------------

subroutine get_close_maxdist_init(gc, maxdist)

type(get_close_type), intent(inout) :: gc
real(r8),             intent(in)    :: maxdist

gc%maxdist = maxdist

end subroutine get_close_maxdist_init

!----------------------------------------------------------------------------

subroutine get_close_obs_init(gc, num, obs)

! Initializes part of get_close accelerator that depends on the particular obs

type(get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
type(location_type),  intent(in)    :: obs(num)

end subroutine get_close_obs_init


!----------------------------------------------------------------------------

subroutine alloc_get_close_obs(num, obs, cutoff, obs_box)
! CURRENTLY UNUSED, I BELIEVE.

integer, intent(in) :: num
type(location_type), intent(in) :: obs(num)
real(r8), intent(in) :: cutoff
integer, intent(out) :: obs_box(num)

! There might not need to be code here but if the get_close_obs() call
! gets too slow, precomputing can be done here.

! set this to satisfy the intent(out) directive.
obs_box(:) = 0
return

end subroutine alloc_get_close_obs


!----------------------------------------------------------------------------

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs, obs_kind, &
   num_close, close_ind, dist)

! Return how many locations are closer than cutoff distance, along with a
! count, and the actual distances if requested.  The kinds are available if
! more sophisticated distance computations are wanted.

type(get_close_type),             intent(in)  :: gc
type(location_type),              intent(in)  :: base_obs_loc, obs(:)
integer,                          intent(in)  :: base_obs_kind, obs_kind(:)
integer,                          intent(out) :: num_close, close_ind(:)
real(r8),             optional,   intent(out) :: dist(:)

integer :: i
real(r8) :: this_dist

! Return list of obs that are within cutoff and their distances
num_close = 0
do i = 1, size(obs)
   this_dist = get_dist(base_obs_loc, obs(i))
   if(this_dist < gc%maxdist) then
      ! Add this ob to the list
      num_close = num_close + 1
      close_ind(num_close) = i
      if (present(dist)) dist(num_close) = this_dist 
   endif
end do

end subroutine get_close_obs



function is_location_in_region(loc, minl, maxl)
!----------------------------------------------------------------------------
!
! Returns true if the given location is between the other two.

logical                          :: is_location_in_region
type(location_type), intent(in)  :: loc, minl, maxl


if ( .not. module_initialized ) call initialize_module

if ((minl%which_vert /= maxl%which_vert) .or. &
    (minl%which_vert /= loc%which_vert)) then
   write(errstring,*)'which_vert (',loc%which_vert,') must be same in all args'
   call error_handler(E_ERR, 'is_location_in_region', errstring, source, revision, revdate)
endif

! assume failure and return as soon as we are confirmed right.
! set to success only at the bottom after all tests have passed.
is_location_in_region = .false.

if ((loc%azm < minl%azm) .or. (loc%azm > maxl%azm)) return
if ((loc%rad < minl%rad) .or. (loc%rad > maxl%rad)) return
if ((loc%vloc < minl%vloc) .or. (loc%vloc > maxl%vloc)) return
 
is_location_in_region = .true.

end function is_location_in_region


!----------------------------------------------------------------------------
! end of location/annulus/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod
