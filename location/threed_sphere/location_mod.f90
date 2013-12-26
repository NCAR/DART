! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module location_mod

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

use      types_mod, only : r8, MISSING_R8, MISSING_I, PI, RAD2DEG, DEG2RAD
use  utilities_mod, only : register_module, error_handler, E_ERR, ascii_file_format, &
                           nc_check, E_MSG, open_file, close_file, set_output,       &
                           logfileunit, nmlfileunit, find_namelist_in_file,          &
                           check_namelist_read, do_output, do_nml_file,              &
                           do_nml_term, is_longitude_between
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform
use   obs_kind_mod, only : get_num_obs_kinds, get_obs_kind_name
use mpi_utilities_mod, only : my_task_id, task_count

implicit none
private

public :: location_type, get_location, set_location, &
          set_location_missing, is_location_in_region, &
          write_location, read_location, interactive_location, query_location, &
          LocationDims, LocationName, LocationLName, get_close_obs, &
          get_close_maxdist_init, get_close_obs_init, get_close_type, &
          operator(==), operator(/=), get_dist, get_close_obs_destroy, &
          nc_write_location_atts, nc_get_location_varids, nc_write_location, &
          vert_is_height, vert_is_pressure, vert_is_undef, vert_is_level, &
          vert_is_surface, vert_is_scale_height, has_vertical_localization, &
          VERTISUNDEF, VERTISSURFACE, VERTISLEVEL, VERTISPRESSURE, &
          VERTISHEIGHT, VERTISSCALEHEIGHT, print_get_close_type, horiz_dist_only, &
          location_type_distrib,copy_location_type, get_seq_distrib_vloc, &
          get_seq_distrib_lon, get_seq_distrib_lat, set_seq_distrib_vloc


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! The possible values for the location_type%which_vert component.
! These are intended to be PRIVATE to this module. Do not make public.

integer, parameter :: VERTISUNDEF       = -2 ! has no vertical location (undefined)
integer, parameter :: VERTISSURFACE     = -1 ! surface value
integer, parameter :: VERTISLEVEL       =  1 ! by level
integer, parameter :: VERTISPRESSURE    =  2 ! by pressure
integer, parameter :: VERTISHEIGHT      =  3 ! by height
integer, parameter :: VERTISSCALEHEIGHT =  4 ! by scale height

type location_type
   !private
   real(r8) :: lon, lat, vloc ! lon, lat stored in radians
   integer  :: which_vert     ! determines if by level, height, pressure, ...
end type location_type

type location_type_distrib !HK
   private
   real(r8) :: lon, lat, vloc ! lon, lat stored in radians
   integer  :: which_vert     ! determines if by level, height, pressure, ...
   integer  :: which_vert_localize ! for localization, not sure if we need this
   real(r8) :: localize_vloc  ! height for localization
end type


! Type to facilitate efficient computation of observations close to a given location
type get_close_type
   private
   integer           :: num
   real(r8)          :: maxdist
   integer, pointer  :: lon_offset(:, :)     ! (nlat, nlat); 
   integer, pointer  :: obs_box(:)           ! (nobs); List of obs indices in boxes
   integer, pointer  :: count(:, :)          ! (nlon, nlat); # of obs in each box
   integer, pointer  :: start(:, :)          ! (nlon, nlat); Start of list of obs in this box
   real(r8)          :: bot_lat, top_lat     ! Bottom and top latitudes of latitude boxes
   real(r8)          :: bot_lon, top_lon     ! Bottom and top longitudes of longitude boxes
   real(r8)          :: lon_width, lat_width ! Width of boxes in lon and lat
   logical           :: lon_cyclic           ! Do boxes wraparound in longitude?
   integer           :: num_types            ! if > 0, cutoffs per type
   real(r8), pointer :: special_maxdist(:)
   real(r8)          :: original_maxdist
end type get_close_type

type(random_seq_type) :: ran_seq
logical               :: ran_seq_init = .false.
logical, save         :: module_initialized = .false.

integer,              parameter :: LocationDims = 3
character(len = 129), parameter :: LocationName = "loc3Dsphere"
character(len = 129), parameter :: LocationLName = &
                                   "threed sphere locations: lon, lat, vertical"

character(len = 512) :: errstring

! Global storage for vertical distance normalization factors
real(r8)              :: vert_normalization(4)
real(r8), allocatable :: special_vert_norm(:,:)  ! if doing per-type

! Global storage for fast approximate sin and cosine lookups
! PAR For efficiency for small cases might want to fill tables as needed
real(r8) :: my_sin(-630:630), my_cos(-630:630), my_acos(-1000:1000)

! If maxdist stays the same, don't need to do box distance calculations
integer :: last_maxdist = -1.0

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
! nlon                            -> Number longitude boxes for get_close_obs 
!                                    nlon MUST BE ODD
! nlat                            -> Number latitude boxes for get_close_obs

logical  :: horiz_dist_only                 = .true.
real(r8) :: vert_normalization_pressure     = 100000.0_r8
real(r8) :: vert_normalization_height       = 10000.0_r8
real(r8) :: vert_normalization_level        = 20.0_r8
real(r8) :: vert_normalization_scale_height = 5.0_r8
logical  :: maintain_original_vert          = .false. 
logical  :: approximate_distance            = .false.
integer  :: nlon                            = 71
integer  :: nlat                            = 36
logical  :: output_box_info                 = .false.
integer  :: print_box_level                 = 0
! obsolete now - code fixed.  leave in for backwards compatibility
! but does nothing now.
logical  :: num_tasks_insensitive           = .false.

namelist /location_nml/ horiz_dist_only, vert_normalization_pressure, &
   vert_normalization_height, vert_normalization_level,               &
   vert_normalization_scale_height, approximate_distance, nlon, nlat, &
   output_box_info, print_box_level, num_tasks_insensitive,           &
   maintain_original_vert

!-----------------------------------------------------------------

interface operator(==); module procedure loc_eq; end interface
interface operator(/=); module procedure loc_ne; end interface

interface set_location
   module procedure set_location_single
   module procedure set_location_array
end interface set_location

contains

!----------------------------------------------------------------------------

subroutine initialize_module
 
! things which need doing exactly once.

integer :: iunit, io, i
character(len=129) :: str1

if (module_initialized) return

call register_module(source, revision, revdate)
module_initialized = .true.

! Read the namelist entry
call find_namelist_in_file("input.nml", "location_nml", iunit)
read(iunit, nml = location_nml, iostat = io)
call check_namelist_read(iunit, io, "location_nml")

! Write the namelist values to the log file

if(do_nml_file()) write(nmlfileunit, nml=location_nml)
if(do_nml_term()) write(     *     , nml=location_nml)

! deprecated namelist item
if (num_tasks_insensitive) then
  call error_handler(E_MSG, 'location_mod:', &
                    'WARNING: namelist item "num_tasks_insensitive" is deprecated and will be removed in the future')
endif

! Make sure that the number of longitudes, nlon, for get_close_obs is odd
if(nlon / 2 * 2 == nlon) then
   call error_handler(E_ERR, 'initialize_module', 'nlon must be odd', &
      source, revision, revdate)
endif

! Copy the normalization factors in the vertical into an array
vert_normalization(1) = vert_normalization_level
vert_normalization(2) = vert_normalization_pressure
vert_normalization(3) = vert_normalization_height
vert_normalization(4) = vert_normalization_scale_height

if (horiz_dist_only) then
   call error_handler(E_MSG,'location_mod:', &
      'Ignoring vertical when computing distances; horizontal only', &
      source,revision,revdate)
else
   call error_handler(E_MSG,'location_mod:', &
      'Including vertical separation when computing distances:', &
      source,revision,revdate)
   write(str1,'(A,f17.5)') '       # pascals ~ 1 horiz radian: ', vert_normalization_pressure
   call error_handler(E_MSG,'location_mod:',str1,source,revision,revdate)
   write(str1,'(A,f17.5)') '        # meters ~ 1 horiz radian: ', vert_normalization_height
   call error_handler(E_MSG,'location_mod:',str1,source,revision,revdate)
   write(str1,'(A,f17.5)') '  # model levels ~ 1 horiz radian: ', vert_normalization_level
   call error_handler(E_MSG,'location_mod:',str1,source,revision,revdate)
   write(str1,'(A,f17.5)') ' # scale heights ~ 1 horiz radian: ', vert_normalization_scale_height
   call error_handler(E_MSG,'location_mod:',str1,source,revision,revdate)
endif

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
   call error_handler(E_MSG,'location_mod:', &
      'Using table-lookup approximation for distance computations', &
      source,revision,revdate)
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
! distance computation for incompatible vertical location types results 
! in a fatal error unless one of the vertical types is UNDEFINED.

! The 3rd argument is actually a specific type, (e.g. RADIOSONDE_TEMPERATURE, 
! AIRCRAFT_TEMPERATURE).  The 4th argument is a generic kind (e.g. KIND_TEMPERATURE).
! The types/kinds are part of
! the interface in case user-code wants to do a more sophisticated distance
! calculation based on the base or target types.  In the usual case this
! code still doesn't use the types, but there's an undocumented feature that
! allows you to maintain the original vertical normalization even when
! changing the cutoff distance in the horizontal.  For that to work we
! do need to know the type, and we use the type of loc1 to control it.
! 

type(location_type), intent(in) :: loc1, loc2
integer, optional,   intent(in) :: type1, kind2
logical, optional,   intent(in) :: no_vert
real(r8)                        :: get_dist

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
   ! This is new - if per-type localization distances given, use the kind of loc1
   ! to determine the vertical mapping distance.  it defaults to the 4 standard ones,
   ! but can be specified separately if desired.

   ! note that per-type vertical conversion factors are a global here, and were set
   ! by the last call to gc_init that specified per/type cutoffs.
   if (allocated(special_vert_norm)) then 
      vert_dist = abs(loc1%vloc - loc2%vloc) / special_vert_norm(loc1%which_vert, type1)
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
 
! Given a location type (in radians), 
! return the longitude, latitude (in degrees) and level 

type(location_type), intent(in) :: loc
real(r8), dimension(3) :: get_location

if ( .not. module_initialized ) call initialize_module

get_location(1) = loc%lon * RAD2DEG                 
get_location(2) = loc%lat * RAD2DEG                 
get_location(3) = loc%vloc     

end function get_location

!---------------------------------------------------------------------------

!FIXME: this routine does not exist in any other locations module, and i cannot
! find any code which uses it.  I propose we obsolete this code.  get_location()
! returns the contents of a locations type in a location-independent routine.
function get_location_lon(loc)
 
! Given a location type, return the longitude.  Values stored in radians but
! returned in degrees.

type(location_type), intent(in) :: loc
real(r8) :: get_location_lon

if ( .not. module_initialized ) call initialize_module

get_location_lon = loc%lon * RAD2DEG    

end function get_location_lon

!---------------------------------------------------------------------------

!FIXME: ditto for comment above.  this should be deprecated.
function get_location_lat(loc)
 
! Given a location type, return the latitude.  Values stored in radians but
! returned in degrees.

type(location_type), intent(in) :: loc
real(r8) :: get_location_lat

if ( .not. module_initialized ) call initialize_module

get_location_lat = loc%lat * RAD2DEG      

end function get_location_lat

!----------------------------------------------------------------------------

function set_location_single(lon, lat, vert_loc,  which_vert)
 
! Puts the given longitude, latitude, and vertical location
! into a location datatype.  Arguments to this function are in degrees,
! but the values are stored as radians.

real(r8), intent(in) :: lon, lat
real(r8), intent(in) :: vert_loc
integer,  intent(in) :: which_vert
type (location_type) :: set_location_single

if ( .not. module_initialized ) call initialize_module

if(lon < 0.0_r8 .or. lon > 360.0_r8) then
   write(errstring,*)'longitude (',lon,') is not within range [0,360]'
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
endif

if(lat < -90.0_r8 .or. lat > 90.0_r8) then
   write(errstring,*)'latitude (',lat,') is not within range [-90,90]'
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
endif

set_location_single%lon = lon * DEG2RAD
set_location_single%lat = lat * DEG2RAD

if(which_vert < VERTISUNDEF .or. which_vert == 0 .or. which_vert > VERTISSCALEHEIGHT) then
   write(errstring,*)'which_vert (',which_vert,') must be one of -2, -1, 1, 2, 3, or 4'
   call error_handler(E_ERR,'set_location', errstring, source, revision, revdate)
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

if ( .not. module_initialized ) call initialize_module

if (size(list) < 4) then
   write(errstring,*)'requires 4 input values'
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
endif

set_location_array = set_location_single(list(1), list(2), list(3), nint(list(4)))

end function set_location_array

!----------------------------------------------------------------------------

function set_location_missing()

! Initialize a location type to indicate the contents are unset.

type (location_type) :: set_location_missing

if ( .not. module_initialized ) call initialize_module

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
          source, revision, revdate)
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

if ( .not. module_initialized ) call initialize_module

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
      'Cannot use string buffer with binary format', &
       source, revision, revdate)
endif

! format the location to be more human-friendly; meaning
! degrees instead of radians, and kilometers for height,
! hectopascals instead of pascals for pressure, etc.

! this must be the sum of the longest of the formats below.
charlength = 72

if (len(charstring) < charlength) then
   write(errstring, *) 'charstring buffer must be at least ', charlength, ' chars long'
   call error_handler(E_ERR, 'write_location', errstring, source, revision, revdate)
endif

! format the horizontal into a temp string
write(string1, '(A,F12.8,1X,F12.8,1X,A)') 'Lon/Lat(deg): ',  loc%lon*RAD2DEG, &
   loc%lat*RAD2DEG, ' Vert:'

! then pretty up the vertical choices, trying to get them to line up in
! case the caller is listing out locations with different vert types.
! concatinate the vertical on the end of the horizontal and put it all
! into the return string.
select case  (loc%which_vert)
   case (VERTISUNDEF)
      write(charstring, '(A,1X,A)')       trim(string1), '              Undefined'
   case (VERTISSURFACE)
      write(charstring, '(A,1X,F12.5,A)') trim(string1), loc%vloc, '  surface (m)'
   case (VERTISLEVEL)
      write(charstring, '(A,1X,F5.0,A)')  trim(string1), loc%vloc, '         level'
   case (VERTISPRESSURE)
      write(charstring, '(A,1X,F12.7,A)') trim(string1), loc%vloc / 100.0_r8, '  hPa'
   case (VERTISHEIGHT)
      write(charstring, '(A,1X,F12.7,A)') trim(string1), loc%vloc / 1000.0_r8, '  km'
   case (VERTISSCALEHEIGHT)
      write(charstring, '(A,1X,F12.7,A)') trim(string1), loc%vloc, '  scale ht'
   case default
      write(errstring, *) 'unrecognized key for vertical type: ', loc%which_vert
      call error_handler(E_ERR, 'write_location', errstring, source, revision, revdate)
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

if ( .not. module_initialized ) call initialize_module

if (ascii_file_format(fform)) then
   read(locfile, '(a5)' ) header
   if(header /= 'loc3d') then
         write(errstring,*)'Expected location header "loc3d" in input file, got ', header
      call error_handler(E_ERR, 'read_location', errstring, source, revision, revdate)
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

if ( .not. module_initialized ) call initialize_module

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

!----------------------------------------------------------------------------

function nc_write_location_atts( ncFileID, fname, ObsNumDimID ) result (ierr)
 
! Writes the "location module" -specific attributes to a netCDF file.

use typeSizes
use netcdf

integer,          intent(in) :: ncFileID     ! handle to the netcdf file
character(len=*), intent(in) :: fname        ! file name (for printing purposes)
integer,          intent(in) :: ObsNumDimID  ! handle to the dimension that grows
integer                      :: ierr

integer :: LocDimID
integer :: VarID

if ( .not. module_initialized ) call initialize_module

ierr = -1 ! assume things will fail ...

! define the rank/dimension of the location information
call nc_check(nf90_def_dim(ncid=ncFileID, name='location', len=LocationDims, &
       dimid = LocDimID), 'nc_write_location_atts', 'def_dim:location '//trim(fname))

! Define the observation location variable and attributes

call nc_check(nf90_def_var(ncid=ncFileID, name='location', xtype=nf90_double, &
          dimids=(/ LocDimID, ObsNumDimID /), varid=VarID), &
            'nc_write_location_atts', 'location:def_var')

call nc_check(nf90_put_att(ncFileID, VarID, 'description', &
        'location coordinates'), 'nc_write_location_atts', 'location:description')
call nc_check(nf90_put_att(ncFileID, VarID, 'location_type', &
        trim(LocationName)), 'nc_write_location_atts', 'location:location_type')
call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', &
        trim(LocationLName)), 'nc_write_location_atts', 'location:long_name')
call nc_check(nf90_put_att(ncFileID, VarID, 'storage_order',     &
        'Lon Lat Vertical'), 'nc_write_location_atts', 'location:storage_order')
call nc_check(nf90_put_att(ncFileID, VarID, 'units',     &
        'degrees degrees which_vert'), 'nc_write_location_atts', 'location:units')

! Define the ancillary vertical array and attributes

call nc_check(nf90_def_var(ncid=ncFileID, name='which_vert', xtype=nf90_int, &
          dimids=(/ ObsNumDimID /), varid=VarID), &
            'nc_write_location_atts', 'which_vert:def_var')

call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', 'vertical coordinate system code'), &
           'nc_write_location_atts', 'which_vert:long_name')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISUNDEF', VERTISUNDEF), &
           'nc_write_location_atts', 'which_vert:VERTISUNDEF')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISSURFACE', VERTISSURFACE), &
           'nc_write_location_atts', 'which_vert:VERTISSURFACE')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISLEVEL', VERTISLEVEL), &
           'nc_write_location_atts', 'which_vert:VERTISLEVEL')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISPRESSURE', VERTISPRESSURE), &
           'nc_write_location_atts', 'which_vert:VERTISPRESSURE')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISHEIGHT', VERTISHEIGHT), &
           'nc_write_location_atts', 'which_vert:VERTISHEIGHT')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISSCALEHEIGHT', VERTISSCALEHEIGHT), &
           'nc_write_location_atts', 'which_vert:VERTISSCALEHEIGHT')

! If we made it to here without error-ing out ... we're good.

ierr = 0

end function nc_write_location_atts

!----------------------------------------------------------------------------

subroutine nc_get_location_varids( ncFileID, fname, LocationVarID, WhichVertVarID )

! Return the LocationVarID and WhichVertVarID variables from a given netCDF file.
!
! ncFileId         the netcdf file descriptor
! fname            the name of the netcdf file (for error messages only)
! LocationVarID    the integer ID of the 'location' variable in the netCDF file
! WhichVertVarID   the integer ID of the 'which_vert' variable in the netCDF file

use typeSizes
use netcdf

integer,          intent(in)  :: ncFileID   ! handle to the netcdf file
character(len=*), intent(in)  :: fname      ! file name (for printing purposes)
integer,          intent(out) :: LocationVarID, WhichVertVarID

if ( .not. module_initialized ) call initialize_module

call nc_check(nf90_inq_varid(ncFileID, 'location', varid=LocationVarID), &
          'nc_get_location_varids', 'inq_varid:location '//trim(fname))

call nc_check(nf90_inq_varid(ncFileID, 'which_vert', varid=WhichVertVarID), &
          'nc_get_location_varids', 'inq_varid:which_vert '//trim(fname))

end subroutine nc_get_location_varids

!----------------------------------------------------------------------------

subroutine nc_write_location(ncFileID, LocationVarID, loc, obsindex, WhichVertVarID)
 
! Writes a SINGLE location to the specified netCDF variable and file.
! The LocationVarID and WhichVertVarID must be the values returned from
! the nc_get_location_varids call.

use typeSizes
use netcdf

integer,             intent(in) :: ncFileID, LocationVarID
type(location_type), intent(in) :: loc
integer,             intent(in) :: obsindex
integer,             intent(in) :: WhichVertVarID

real(r8), dimension(LocationDims) :: locations
integer,  dimension(1) :: intval

if ( .not. module_initialized ) call initialize_module

locations = get_location( loc ) ! converts from radians to degrees, btw

call nc_check(nf90_put_var(ncFileID, LocationVarId, locations, &
          start=(/ 1, obsindex /), count=(/ LocationDims, 1 /) ), &
            'nc_write_location', 'put_var:location')

intval = loc%which_vert
call nc_check(nf90_put_var(ncFileID, WhichVertVarID, intval, &
          start=(/ obsindex /), count=(/ 1 /) ), &
            'nc_write_location','put_var:vert' )

end subroutine nc_write_location

!----------------------------------------------------------------------------

subroutine get_close_obs_init(gc, num, obs)
 
! Initializes part of get_close accelerator that depends on the particular obs

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
! special for grid-decomposition debugging; set print level to -8.
if (output_box_info) then
   ! if this task normally prints, call the print routine.
   ! if print level > 2, set all tasks to print and call print.
   ! then reset the status to off again.
   if (do_output()) then
      call print_get_close_type(gc, print_box_level)
   else if (print_box_level >= 2 .or. print_box_level < 0) then
      ! print status was false, but turn on temporarily
      ! to output box info from all tasks.
      call set_output(.true.)
      call print_get_close_type(gc, print_box_level)
      call set_output(.false.)
   endif
endif

end subroutine get_close_obs_init

!----------------------------------------------------------------------------

subroutine get_close_obs_destroy(gc)

type(get_close_type), intent(inout) :: gc

deallocate(gc%obs_box, gc%lon_offset, gc%count, gc%start)
if (gc%num_types > 0) deallocate(gc%special_maxdist)
! since this is a global, keep it around.  it is always allocated
! the same size and can be reused for any gc.
!if (allocated(special_vert_norm)) deallocate(special_vert_norm)

end subroutine get_close_obs_destroy

!----------------------------------------------------------------------------

subroutine get_close_maxdist_init(gc, maxdist, maxdist_list)

type(get_close_type), intent(inout) :: gc
real(r8),             intent(in)    :: maxdist
real(r8), intent(in), optional      :: maxdist_list(:)

character(len=129) :: str1
integer :: i
! a bit of a hack - it only prints the first time this routine is
! called.  in filter it's called twice with the same args and this
! code is counting on that.  
logical, save :: firsttime = .true.

! Allocate the storage for the grid dependent boxes
allocate(gc%lon_offset(nlat, nlat), gc%count(nlon, nlat), gc%start(nlon, nlat))
gc%lon_offset = -1
gc%count      = -1
gc%start      = -1

! set the default value.  if there is a list and any distance in the list 
! is larger, use that instead.  (the boxes need to be calculated based on 
! the largest possible distance).
gc%maxdist = maxdist
gc%original_maxdist = maxdist
gc%num_types = 0

if (present(maxdist_list)) then
   gc%num_types = get_num_obs_kinds()
   if (size(maxdist_list) .ne. gc%num_types) then
      write(errstring,'(A,I8,A,I8)')'maxdist_list len must equal number of kinds, ', &
                                    size(maxdist_list), ' /= ', gc%num_types
      call error_handler(E_ERR, 'get_close_maxdist_init', errstring, source, revision, revdate)
   endif
   allocate(gc%special_maxdist(gc%num_types))
   if (.not.allocated(special_vert_norm)) &
      allocate(special_vert_norm(4, gc%num_types))  ! fill from namelist here, or assim_tools?
  
   gc%special_maxdist(:) = maxdist_list(:)
   gc%maxdist = maxval(gc%special_maxdist)

   ! by default, the vertical changes along with the horizontal to keep 
   ! the aspect ratio of the elipsoid constant.
   special_vert_norm(1, :) = vert_normalization_level
   special_vert_norm(2, :) = vert_normalization_pressure
   special_vert_norm(3, :) = vert_normalization_height
   special_vert_norm(4, :) = vert_normalization_scale_height

   ! keep the original vertical distance constant even as you change
   ! the horizontal localization distance.
   if (maintain_original_vert) then
      special_vert_norm(1, :) = vert_normalization_level        * &
                                   (gc%original_maxdist / maxdist_list(:))
      special_vert_norm(2, :) = vert_normalization_pressure     * &
                                   (gc%original_maxdist / maxdist_list(:))
      special_vert_norm(3, :) = vert_normalization_height       * & 
                                   (gc%original_maxdist / maxdist_list(:))
      special_vert_norm(4, :) = vert_normalization_scale_height * & 
                                   (gc%original_maxdist / maxdist_list(:))
      if (firsttime) then
         do i = 1, gc%num_types
            if (abs(gc%original_maxdist - maxdist_list(i)) < 1.0e-14_r8) cycle 
    
            write(str1,'(2A)') 'Altering vertical normalization for type ', get_obs_kind_name(i)
            call error_handler(E_MSG,'location_mod:',str1,source,revision,revdate)
            write(str1,'(A,f17.5)') '       # pascals ~ 1 horiz radian: ', special_vert_norm(2, i)
            call error_handler(E_MSG,'location_mod:',str1,source,revision,revdate)
            write(str1,'(A,f17.5)') '        # meters ~ 1 horiz radian: ', special_vert_norm(3, i)
            call error_handler(E_MSG,'location_mod:',str1,source,revision,revdate)
            write(str1,'(A,f17.5)') '  # model levels ~ 1 horiz radian: ', special_vert_norm(1, i)
            call error_handler(E_MSG,'location_mod:',str1,source,revision,revdate)
            write(str1,'(A,f17.5)') ' # scale heights ~ 1 horiz radian: ', special_vert_norm(4, i)
            call error_handler(E_MSG,'location_mod:',str1,source,revision,revdate)
         enddo
         firsttime = .false.
      endif
   endif
endif

end subroutine get_close_maxdist_init

!----------------------------------------------------------------------------

subroutine get_close_obs(gc, base_obs_loc, base_obs_type, obs, obs_kind, &
   num_close, close_ind, dist)

! In spite of the names, the specific types are available to do a more
! sophisticated distance computation if needed.

type(get_close_type), intent(in)  :: gc
type(location_type),  intent(in)  :: base_obs_loc, obs(:)
integer,              intent(in)  :: base_obs_type, obs_kind(:)
integer,              intent(out) :: num_close, close_ind(:)
real(r8), optional,   intent(out) :: dist(:)

! If dist is NOT present, just find everybody in a box, put them in the list,
! but don't compute any distances

integer :: lon_box, lat_box, i, j, k, n_lon, lon_ind, n_in_box, st, t_ind
real(r8) :: this_dist, this_maxdist

! Variables needed for comparing against correct case
integer :: cnum_close, cclose_ind(size(obs))
real(r8) :: cdist(size(obs))

! First, set the intent out arguments to a missing value
num_close = 0
close_ind = -99
if(present(dist)) dist = -99.0_r8
this_dist = 999999.0_r8   ! something big.

! the list of locations in the obs() argument must be the same
! as the list of locations passed into get_close_obs_init(), so
! gc%num and size(obs) better be the same.   if the list changes,
! you have to destroy the old gc and init a new one.
if (size(obs) /= gc%num) then
   write(errstring,*)'obs() array must match one passed to get_close_obs_init()'
   call error_handler(E_ERR, 'get_close_obs', errstring, source, revision, revdate)
endif

! If num == 0, no point in going any further. 
if (gc%num == 0) return


! set a local variable for what the maxdist is in this particular case.
! if per-type distances are set, use those.  otherwise, use the global val.
if (gc%num_types > 0) then
   this_maxdist = gc%special_maxdist(base_obs_type)
else
   this_maxdist = gc%maxdist
endif

!--------------------------------------------------------------
! For validation, it is useful to be able to compare against exact
! exhaustive search
if(COMPARE_TO_CORRECT) then
   cnum_close = 0
   do i = 1, gc%num 
   this_dist = get_dist(base_obs_loc, obs(i), base_obs_type, obs_kind(i))
      if(this_dist <= this_maxdist) then
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
                  this_dist = get_dist(base_obs_loc, obs(t_ind), base_obs_type, obs_kind(t_ind))
               else
               ! Otherwise can just get horizontal distance
                  this_dist = get_dist(base_obs_loc, obs(t_ind), base_obs_type, obs_kind(t_ind), &
                     no_vert = .true.)
               endif
            else
               ! Dist isn't present; add this ob to list without computing distance
               num_close = num_close + 1
               close_ind(num_close) = t_ind
            endif

            ! If dist is present and this obs' distance is less than cutoff, add it in list
            if(present(dist) .and. this_dist <= this_maxdist) then
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
   if((num_close /= cnum_close) .and. present(dist)) then
      write(errstring, *) 'get_close (', num_close, ') should equal exhaustive search (', cnum_close, ')'
      call error_handler(E_ERR, 'get_close_obs', errstring, source, revision, revdate, &
                         text2='optional arg "dist" is present; we are computing exact distances', &
                         text3='the exhaustive search should find an identical number of locations')
   else if (num_close < cnum_close) then
      write(errstring, *) 'get_close (', num_close, ') should not be smaller than exhaustive search (', cnum_close, ')'
      call error_handler(E_ERR, 'get_close_obs', errstring, source, revision, revdate, &
                         text2='optional arg "dist" not present; we are returning a superset of close locations', &
                         text3='the exhaustive search should find an equal or lesser number of locations')
   endif
endif
!--------------------End of verify by comparing to exhaustive search --------------


end subroutine get_close_obs

!--------------------------------------------------------------------------

subroutine find_box_ranges(gc, obs, num)
 
! Finds boundaries for boxes in N/S direction. If data is localized in N/S
! tries to find boxes that only span the range of the data.
  
type(get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
type(location_type),  intent(in)    :: obs(num)

real(r8) :: min_lat, max_lat, beg_box_lon, end_box_lon, first_obs_lon, last_obs_lon
real(r8) :: longitude_range, degrees, lon_dist
integer  :: i, indx, gap_start, gap_end, gap_length
logical  :: lon_box_full(360), old_out

! Initialize boxes used to see where observations are.
! Assume regional until we prove that we have to use the
! entire 360 in longitude.
lon_box_full = .false.
gc%lon_cyclic = .false.

! Figure out domain over which an additional obs MIGHT be close to one in this set
! If any points within maxdist of the poles, our boxes have to cover all 360 of
! longitude - no point in trying to restrict boxes to a region of the globe.
min_lat = minval(obs(:)%lat) - gc%maxdist
max_lat = maxval(obs(:)%lat) + gc%maxdist
if(min_lat <= -PI / 2.0_r8) then
   min_lat = -PI / 2.0_r8
   gc%lon_cyclic = .true.
endif 
if(max_lat >= PI / 2.0_r8) then
   max_lat = PI / 2.0_r8
   gc%lon_cyclic = .true.
endif

! Put this into storage for this get_close_type
gc%bot_lat = min_lat
gc%top_lat = max_lat
gc%lat_width = (max_lat - min_lat) / nlat
! don't have to do all this work if we already know it has to be cyclic
if (.not. gc%lon_cyclic) then

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
   
   if (gap_length > 0) then

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
      ! If addition of halos for close obs fills more than half of space 
      ! things go 0 to 2PI

      ! other places we are computing in radians.  here we are computing in
      ! lat/lon, and you can't just add maxdist to the edges - that doesn't
      ! take into account the great-circle distance.  the separation in longitude
      ! varies with latitude.  compute the delta longitude based on the most
      ! poleward latitude and add that onto the edges of both boxes.  that
      ! overestimates for points closer to the equator, but that's better
      ! than underestimating and excluding points that are within maxdist.
      lon_dist = find_del_lon(minval(obs(:)%lat), maxval(obs(:)%lat), gc%maxdist)

      if(longitude_range + 2.0_r8 * lon_dist > PI) then
         gc%lon_cyclic = .true.
      else
         first_obs_lon = first_obs_lon - lon_dist
         if(first_obs_lon < 0.0_r8    ) first_obs_lon = first_obs_lon + 2.0_r8 * PI
         last_obs_lon  = last_obs_lon + lon_dist
         if(last_obs_lon > 2.0_r8 * PI) last_obs_lon  = last_obs_lon  - 2.0_r8 * PI
         gc%lon_cyclic = .false.
      endif
   else
      ! No gap was found: all 360 boxes had an observation in them
      gc%lon_cyclic = .true.
   endif
endif

if (gc%lon_cyclic) then
   first_obs_lon = 0.0_r8
   last_obs_lon  = 2.0_r8 * PI
endif

! Put in storage for structure
gc%bot_lon = first_obs_lon
gc%top_lon = last_obs_lon
longitude_range = last_obs_lon - first_obs_lon
if(longitude_range <= 0.0_r8) longitude_range = longitude_range + 2.0_r8 * PI
gc%lon_width = longitude_range / nlon

if(COMPARE_TO_CORRECT) then
   old_out = do_output()
   call set_output(.true.)
   write(errstring, *) 'lat bot, top, width ', gc%bot_lat, gc%top_lat, gc%lat_width
   call error_handler(E_MSG, 'find_box_ranges', errstring)
   write(errstring, *) 'lon bot, top, width ', gc%bot_lon, gc%top_lon, gc%lon_width
   call error_handler(E_MSG, 'find_box_ranges', errstring)
   call set_output(old_out)
endif

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

!----------------------------------------------------------------------------

subroutine print_get_close_type(gc, amount)
 
! print out debugging statistics, or optionally print out a full
! dump from all mpi tasks in a format that can be plotted with matlab.

type(get_close_type), intent(in) :: gc
integer, intent(in), optional    :: amount

integer :: i, j, k, first, index, mytask, alltasks
integer :: sample, nfull, nempty, howmuch, total, maxcount, maxi, maxj
logical :: tickmark(gc%num), iam0
real(r8) :: lon_cen, lat_cen

logical, save :: write_now = .true.
integer, save :: been_called = 0
integer :: funit
character(len=64) :: fname

! cumulative times through this routine
been_called = been_called + 1

! second arg is now an int, not logical, and means:
! 0 = very terse, only box summary (default).  
! 1 = structs and first part of arrays.
! 2 = all parts of all arrays.
! -8 = special for grid-decomposition debugging

! by default do not print all the obs_box or start contents (it can
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
!if (gc%num > 100) howmuch = 0

! print the get_close_type derived type values

if (howmuch /= 0 .and. iam0) then
   write(errstring,*) 'get_close_type values:'
   call error_handler(E_MSG, 'loc', errstring)

   write(errstring,*) ' num = ', gc%num
   call error_handler(E_MSG, 'loc', errstring)

   write(errstring,*) ' nlon, nlat = ', nlon, nlat
   call error_handler(E_MSG, 'loc', errstring)

   write(errstring,"(A,F12.6)") ' maxdist = ', gc%maxdist
   call error_handler(E_MSG, 'loc', errstring)
   write(errstring, "(A,3(F12.6))") ' latbox: bot, top, width = ', gc%bot_lat, gc%top_lat, gc%lat_width
   call error_handler(E_MSG, 'loc', errstring)
   write(errstring, "(A,3(F12.6))") ' lonbox: bot, top, width = ', gc%bot_lon, gc%top_lon, gc%lon_width
   call error_handler(E_MSG, 'loc', errstring)

   write(errstring,"(A,F12.6)") ' maxdist = ', RAD2DEG*gc%maxdist
   call error_handler(E_MSG, 'loc', errstring)
   write(errstring, "(A,3(F12.6))") ' latbox: bot, top, width = ', RAD2DEG*gc%bot_lat, RAD2DEG*gc%top_lat, RAD2DEG*gc%lat_width
   call error_handler(E_MSG, 'loc', errstring)
   write(errstring, "(A,3(F12.6))") ' lonbox: bot, top, width = ', RAD2DEG*gc%bot_lon, RAD2DEG*gc%top_lon, RAD2DEG*gc%lon_width
   call error_handler(E_MSG, 'loc', errstring)

   write(errstring,*) ' lon_cyclic = ', gc%lon_cyclic
   call error_handler(E_MSG, 'loc', errstring)
endif

! this one can be very large.   print only the first nth unless
! instructed otherwise.  (print n+1 because 1 more value fits on
! the line because it prints ( i ) and not ( i, j ) like the others.)
if (associated(gc%obs_box)) then
   i = size(gc%obs_box,1)
   if (i/= gc%num) then
      write(errstring,*) ' warning: size of obs_box incorrect, nobs, i =', gc%num, i
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
   if (howmuch > 1) then
      ! DEBUG
      write(errstring,"(A,I8,A,36(I8,1X))") ' obs_box(',i,') =', gc%obs_box(1:min(i,36))  ! (nobs)
      !write(errstring,*) ' obs_box(',i,') =', gc%obs_box    ! (nobs)
      call error_handler(E_MSG, 'locations_mod', errstring)
   else if(howmuch > 0) then
      write(errstring,*) ' obs_box(',i,') =', gc%obs_box(1:min(i,sample+1))
      call error_handler(E_MSG, 'locations_mod', errstring)
      write(errstring,*) '  <rest of obs_box omitted>'
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
else
   if (howmuch > 0) then
      write(errstring,*) ' obs_box unallocated'
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
endif

! like obs_box, this one can be very large.   print only the first nth unless
! instructed otherwise
if (associated(gc%start)) then
   i = size(gc%start,1)
   j = size(gc%start,2)
   if ((i /= nlon) .or. (j /= nlat)) then
      write(errstring,*) ' warning: size of start incorrect, nlon, nlat, i, j =', nlon, nlat, i, j
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
   if (howmuch > 1) then
      write(errstring,*) ' start(',i,j,') ='              ! (nlon, nlat)
      call error_handler(E_MSG, 'locations_mod', errstring)
      do k=1, j
         write(errstring,"(36(I8,1X))") gc%start(1:min(i,36), k)
         call error_handler(E_MSG, 'locations_mod', errstring)
      enddo
   else if (howmuch > 0) then
      write(errstring,*) ' start(',i,j,') =', gc%start(1:min(i,sample), 1)
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
if (associated(gc%lon_offset)) then
   i = size(gc%lon_offset,1)
   j = size(gc%lon_offset,2)
   if ((i /= nlat) .or. (j /= nlat)) then
      write(errstring,*) ' warning: size of lon_offset incorrect, nlat, i, j =', nlat, i, j
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
   if (howmuch > 1) then
      write(errstring,*) ' lon_offset(',i,j,') ='                   ! (nlat, nlat)
      call error_handler(E_MSG, 'locations_mod', errstring)
      do k=1, j
         write(errstring,"(36(I8,1X))") gc%lon_offset(1:min(i,36), k) 
         call error_handler(E_MSG, 'locations_mod', errstring)
      enddo
   else if (howmuch > 0) then
      write(errstring,*) ' lon_offset(',i,j,') =', gc%lon_offset(1:min(i,sample), 1)
      call error_handler(E_MSG, 'locations_mod', errstring)
      write(errstring,*) '  <rest of lon_offset omitted>'
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
else
   if (howmuch > 0) then
      write(errstring,*) ' lon_offset unallocated'
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
endif

! as above, print only first n unless second arg is .true.
if (associated(gc%count)) then
   i = size(gc%count,1)
   j = size(gc%count,2)
   if ((i /= nlon) .or. (j /= nlat)) then
      write(errstring,*) ' warning: size of count incorrect, nlon, nlat, i, j =', &
                      nlon, nlat, i, j
      call error_handler(E_MSG, 'locations_mod', errstring)
   endif
   if (howmuch > 1) then
      write(errstring,*) ' count(',i,j,') ='              ! (nlon, nlat)
      call error_handler(E_MSG, 'locations_mod', errstring)
      do k=1, j
         write(errstring,"(36(I8,1X))") gc%count(1:min(i,36), k) 
         call error_handler(E_MSG, 'locations_mod', errstring)
      enddo
   else if (howmuch > 0) then
      write(errstring,*) ' count(',i,j,') =', gc%count(1:min(i,sample), 1)
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
! in the obs_box list, and complain about duplicates or misses.
tickmark = .FALSE.

do i=1, nlon
   do j=1, nlat
      first = gc%start(i, j)
      do k=1, gc%count(i, j)
         index = first + k - 1
         if ((index < 1) .or. (index > gc%num)) then
            write(errstring, *) 'exiting at first bad value; could be more'
            call error_handler(E_MSG, 'locations_mod', errstring)
            write(errstring, *) 'bad obs list index, in box: ', index, i, j
            call error_handler(E_ERR, 'locations_mod', errstring)
         endif
         if (tickmark(index)) then
            write(errstring, *) 'exiting at first bad value; could be more'
            call error_handler(E_MSG, 'locations_mod', errstring)
            write(errstring, *) 'error: obs found in more than one box list.  index, box: ', &
                         index, i, j
            call error_handler(E_ERR, 'locations_mod', errstring)
         endif
         tickmark(index) = .TRUE.
      enddo
   enddo
enddo

do i=1, gc%num
  if (.not. tickmark(i)) then
     write(errstring, *) 'exiting at first bad value; could be more'
     call error_handler(E_MSG, 'locations_mod', errstring)
     write(errstring,*) 'obs not found in any box list: ', i
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
      lon_cen = gc%bot_lon + ((i-1)*gc%lon_width) + (gc%lon_width/2.0)
      write(funit, '(A,I2,A,I4,A,F12.9,A)') 'xlocs(', i, ',', mytask+1, ') = ',  lon_cen, ';'
   endif
   do j=1, nlat
      if (howmuch == -8 .and. i==1) then
         lat_cen = gc%bot_lat + ((j-1)*gc%lat_width) + (gc%lat_width/2.0)
         write(funit, '(A,I2,A,I4,A,F12.9,A)') 'ylocs(', j, ',', mytask+1, ') = ',  lat_cen, ';'
      endif
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
      ! output for grid boxes; in matlab-friendly format
      if (howmuch == -8) then
         write(funit, '(A,I2,A,I2,A,I4,A,I8,A)') 'boxes(', i, ', ', j, &
                                ',', mytask+1, ') = ', gc%count(i, j), ';'
      endif
   enddo
enddo

if (howmuch == -8) then
   call close_file(funit)
   write_now = .false.
endif

! these print out always - make sure they are useful to end users.
write(errstring, '(a)') "Location module statistics:"
call error_handler(E_MSG, 'locations_mod', errstring)
write(errstring, '(a,i9)') " Total boxes (nlon * nlat): ", nfull + nempty
call error_handler(E_MSG, 'locations_mod', errstring)
write(errstring, '(a,i9)') " Total items to put in boxes: ", gc%num
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
!      write(errstring, '(a,i9,i9)') " That box index: ", maxi, maxj
!      call error_handler(E_MSG, 'locations_mod', errstring)
!   endif
endif


end subroutine print_get_close_type

!----------------------------------------------------------------------------

function is_location_in_region(loc, minl, maxl)
 
! Returns true if the given location is inside the rectangular
! region defined by minl as the lower left, maxl the upper right.
! test is inclusive; values on the edges are considered inside.
! Periodic in longitude (box can cross the 2PI -> 0 line)

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

!---------------------------------------------------------------------------

function vert_is_undef(loc)
 
! Given a location, return true if vertical coordinate is undefined, else false

logical                          :: vert_is_undef
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == VERTISUNDEF) then
   vert_is_undef = .true.
else
   vert_is_undef = .false.
endif

end function vert_is_undef

!---------------------------------------------------------------------------

function vert_is_surface(loc)
 
! Given a location, return true if vertical coordinate is surface, else false

logical                          :: vert_is_surface
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == VERTISSURFACE) then
   vert_is_surface = .true.
else
   vert_is_surface = .false.
endif

end function vert_is_surface

!---------------------------------------------------------------------------

function vert_is_pressure(loc)
 
! Given a location, return true if vertical coordinate is pressure, else false

logical                          :: vert_is_pressure
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == VERTISPRESSURE) then
   vert_is_pressure = .true.
else
   vert_is_pressure = .false.
endif

end function vert_is_pressure

!---------------------------------------------------------------------------

function vert_is_height(loc)
 
! Given a location, return true if vertical coordinate is height, else false

logical                          :: vert_is_height
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == VERTISHEIGHT ) then
   vert_is_height = .true.
else
   vert_is_height = .false.
endif

end function vert_is_height

!---------------------------------------------------------------------------

function vert_is_level(loc)
 
! Given a location, return true if vertical coordinate is level, else false

logical                          :: vert_is_level
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == VERTISLEVEL) then
   vert_is_level = .true.
else
   vert_is_level = .false.
endif

end function vert_is_level

!---------------------------------------------------------------------------

function vert_is_scale_height(loc)

! Given a location, return true if vertical coordinate is scale height, else false

logical                          :: vert_is_scale_height
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == VERTISSCALEHEIGHT ) then
   vert_is_scale_height = .true.
else
   vert_is_scale_height = .false.
endif

end function vert_is_scale_height

!---------------------------------------------------------------------------

function has_vertical_localization()
 
! Return the (opposite) namelist setting for horiz_dist_only.

logical :: has_vertical_localization

if ( .not. module_initialized ) call initialize_module

has_vertical_localization = .not. horiz_dist_only

end function has_vertical_localization

!> trying to reduce the commication in convert vert
!---------------------------------------------------------------------------
subroutine copy_location_type(loc_min, loc_distrib, n)
! This is doing more than copy, it is initiallizing the distributed version

integer,                     intent(in)    :: n
type(location_type),         intent(in)    :: loc_min(n)
type(location_type_distrib), intent(inout) :: loc_distrib(n)

integer :: i

do i = 1, n
   loc_distrib(i)%lon = loc_min(i)%lon
   loc_distrib(i)%lat = loc_min(i)%lat
   loc_distrib(i)%which_vert = loc_min(i)%which_vert
enddo

loc_distrib(:)%which_vert_localize = VERTISPRESSURE ! This assumes CAM watch out.
loc_distrib(:)%localize_vloc = MISSING_R8

end subroutine copy_location_type

!----------------------------------------------------------------------------
function get_seq_distrib_vloc(seq)
! Returns the vertical cooridinate used in localization

type(location_type_distrib) :: seq
real(r8) :: get_seq_distrib_vloc

if ( .not. module_initialized ) call initialize_module

get_seq_distrib_vloc = seq%localize_vloc

end function

!----------------------------------------------------------------------------
function get_seq_distrib_lat(seq)
! Returns the vertical cooridinate used in localization

type(location_type_distrib) :: seq
real(r8) :: get_seq_distrib_lat

if ( .not. module_initialized ) call initialize_module

get_seq_distrib_lat = seq%lat

end function

!----------------------------------------------------------------------------
function get_seq_distrib_lon(seq)
! Returns the vertical cooridinate used in localization

type(location_type_distrib) :: seq
real(r8) :: get_seq_distrib_lon

if ( .not. module_initialized ) call initialize_module

get_seq_distrib_lon = seq%lon

end function

!----------------------------------------------------------------------------
subroutine set_seq_distrib_vloc(seq, vloc)

type(location_type_distrib) :: seq
real(r8) :: vloc

if ( .not. module_initialized ) call initialize_module

seq%localize_vloc = vloc

end subroutine

!----------------------------------------------------------------------------
! end of location/threed_sphere/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
