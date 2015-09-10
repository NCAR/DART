! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module location_mod

! Implements location interfaces for a one dimensional column domain

use      types_mod, only : r8, MISSING_R8, MISSING_I
use  utilities_mod, only : register_module, error_handler, E_ERR, ascii_file_format, &
                           nc_check

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
          vert_is_surface, has_vertical_localization, VERTISUNDEF, VERTISSURFACE, &
          VERTISLEVEL, VERTISPRESSURE, VERTISHEIGHT, &
          set_vert, get_vert, set_which_vert


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! The possible values for the location_type%which_vert component.
! These are intended to be PRIVATE to this module. Do not make public.

integer, parameter :: VERTISUNDEF    = -2 ! has no vertical location (undefined)
integer, parameter :: VERTISSURFACE  = -1 ! surface value
integer, parameter :: VERTISLEVEL    =  1 ! by level
integer, parameter :: VERTISPRESSURE =  2 ! by pressure
integer, parameter :: VERTISHEIGHT   =  3 ! by height

type location_type
   private
   real(r8) :: vloc
   integer  :: which_vert
end type location_type

! Needed as stub but not used in this low-order model
type get_close_type
   private
   integer  :: num
   real(r8) :: maxdist
end type get_close_type

logical, save :: module_initialized = .false.

integer,              parameter :: LocationDims = 1
character(len = 129), parameter :: LocationName = "loc1Dcolumn"
character(len = 129), parameter :: LocationLName = "one-dimensional column"

character(len = 129) :: errstring

interface operator(==); module procedure loc_eq; end interface
interface operator(/=); module procedure loc_ne; end interface

interface set_location
   module procedure set_location_single
   module procedure set_location_array
end interface set_location

contains

!----------------------------------------------------------------------------

subroutine initialize_module
 
if (module_initialized) return

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------------

function get_dist(loc1, loc2, type1, kind2)

! Return distance between 2 locations.  The vertical type must
! be the same for this to return a value.

type(location_type), intent(in) :: loc1, loc2
integer, optional,   intent(in) :: type1, kind2
real(r8)                        :: get_dist

if ( .not. module_initialized ) call initialize_module

! if not same, return very long (but positive) distance: -1*missing_r8
if (loc1%which_vert /= loc2%which_vert) then
   get_dist = -1.0 * MISSING_R8
   return
endif

get_dist = abs(loc1%vloc - loc2%vloc)

end function get_dist

!---------------------------------------------------------------------------

function loc_eq(loc1,loc2)
 
! Interface operator used to compare two locations.
! Returns true only if all components are 'the same' to within machine
! precision.  vert type must also match.

type(location_type), intent(in) :: loc1, loc2
logical                         :: loc_eq

if ( .not. module_initialized ) call initialize_module

loc_eq = .false.

if ( loc1%which_vert /= loc2%which_vert ) return
if ( abs(loc1%vloc - loc2%vloc) > epsilon(loc1%vloc) ) return

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
 
! Given a location type, return the z coordinate

type(location_type), intent(in) :: loc
real(r8)                        :: get_location

if ( .not. module_initialized ) call initialize_module

get_location = loc%vloc

end function get_location

!----------------------------------------------------------------------------

function set_location_single(vert_loc, which_vert)
 
! Store a vertical coordinate and vertical type in a location type.

real(r8), intent(in) :: vert_loc
integer,  intent(in) :: which_vert
type(location_type)  :: set_location_single

if ( .not. module_initialized ) call initialize_module

if (which_vert /= VERTISUNDEF    .and. &
    which_vert /= VERTISSURFACE  .and. &
    which_vert /= VERTISLEVEL    .and. &
    which_vert /= VERTISPRESSURE .and. &
    which_vert /= VERTISHEIGHT) then
   write(errstring,*) 'Illegal value for "which_vert", must be -2, -1, 1, 2, or 3'
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
endif
set_location_single%which_vert = which_vert
set_location_single%vloc = vert_loc

end function set_location_single

!----------------------------------------------------------------------------

function set_location_array(list)
 
! location semi-independent interface routine
! given 2 float numbers, call the underlying set_location routine

real(r8), intent(in) :: list(:)
type(location_type)  :: set_location_array

if ( .not. module_initialized ) call initialize_module

if (size(list) < 2) then
   write(errstring,*) 'requires 2 input values'
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
endif

set_location_array = set_location_single(list(1), nint(list(2)))

end function set_location_array

!----------------------------------------------------------------------------

function set_location_missing()

! Initialize a location type to indicate the contents are unset.

type (location_type) :: set_location_missing

if ( .not. module_initialized ) call initialize_module

set_location_missing%vloc = MISSING_R8
set_location_missing%which_vert = MISSING_I

end function set_location_missing

!---------------------------------------------------------------------------

function query_location(loc,attr)
 
! Returns the value of the attribute
! returns the requested part of the location type

type(location_type),        intent(in) :: loc
character(len=*), optional, intent(in) :: attr
real(r8)                               :: query_location

if ( .not. module_initialized ) call initialize_module

! see the long comment in this routine in the threed_sphere
! module for warnings about compiler bugs before you change
! this code.

query_location = loc%which_vert

if (.not. present(attr)) return

select case(attr)
 case ('which_vert','WHICH_VERT')
   query_location = loc%which_vert
 case ('vloc','VLOC')
   query_location = loc%vloc
 case default
   call error_handler(E_ERR, 'query_location; column', &
         'Only vloc and which_vert are legal attributes to request from location', source, revision, revdate)
end select

end function query_location

!----------------------------------------------------------------------------

subroutine write_location(locfile, loc, fform, charstring)
 
! Writes a location to either a file, or if specified, fills in the character
! string instead.

! additional functionality: if optional argument charstring is specified,
! it must be long enough to hold the string, and the location information is
! written into it instead of to a file.  fform must be ascii (which is the
! default if not specified) to use this option.

integer, intent(in)                        :: locfile
type(location_type), intent(in)            :: loc
character(len = *),  intent(in),  optional :: fform
character(len = *),  intent(out), optional :: charstring

integer             :: charlength
logical             :: writebuf

! 10 format(1x,(f22.14,1x),i4)  ! old
10 format(1X,(G25.16,1X),I2) 

if ( .not. module_initialized ) call initialize_module

! writing to a file (normal use) or to a character buffer?
writebuf = present(charstring)

! output file; test for ascii or binary, write what's asked, and return
if (.not. writebuf) then
   if (ascii_file_format(fform)) then
      write(locfile, '(''loc1c'')' ) 
      write(locfile, 10) loc%vloc, loc%which_vert
   else
      write(locfile) loc%vloc, loc%which_vert
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
! kilometers for height, hectopascals instead of pascals 
! for pressure, etc.


! this must be the sum of the longest of the formats below.
charlength = 32

if (len(charstring) < charlength) then
   write(errstring, *) 'charstring buffer must be at least ', charlength, ' chars long'
   call error_handler(E_ERR, 'write_location', errstring, source, revision, revdate)
endif

select case  (loc%which_vert)
   case (VERTISUNDEF)
      write(charstring, '(A,A)')       'Vert: ', '               Undefined'
   case (VERTISSURFACE)
      write(charstring, '(A,F12.5,A)') 'Vert: ', loc%vloc, '   surface (m)'
   case (VERTISLEVEL)
      write(charstring, '(A,F5.0,A)')  'Vert: ', loc%vloc, '          level'
   case (VERTISPRESSURE)
      write(charstring, '(A,F12.7,A)') 'Vert: ', loc%vloc / 100.0_r8, '   hPa'
   case (VERTISHEIGHT)
      write(charstring, '(A,F12.7,A)') 'Vert: ', loc%vloc / 1000.0_r8, '   km'
   case default
      write(errstring, *) 'unrecognized key for vertical type: ', loc%which_vert
      call error_handler(E_ERR, 'write_location', errstring, source, revision, revdate)
end select


end subroutine write_location

!----------------------------------------------------------------------------

function read_location(locfile, fform)
 
! Reads a location from file that was written by write_location. 
! See write_location for additional discussion.

integer, intent(in)                      :: locfile
type(location_type)                      :: read_location
character(len = *), intent(in), optional :: fform

character(len=5) :: header

if ( .not. module_initialized ) call initialize_module

if (ascii_file_format(fform)) then
   read(locfile, '(a5)' ) header
   if(header /= 'loc1c') then
      write(errstring,*)'Expected location header "loc1c" in input file, got ', header 
      call error_handler(E_ERR, 'read_location', errstring, source, revision, revdate)
   endif
   ! Now read the location data value
   read(locfile, *) read_location%vloc, read_location%which_vert
else
   read(locfile) read_location%vloc, read_location%which_vert
endif

end function read_location

!--------------------------------------------------------------------------

subroutine interactive_location(location, set_to_default)
 
! Allows for interactive input of a location.

type(location_type), intent(out) :: location
logical, intent(in), optional    :: set_to_default

real(r8) :: x

if ( .not. module_initialized ) call initialize_module

! If set_to_default is true, then just zero out and return
if(present(set_to_default)) then
   if(set_to_default) then
      location%vloc = 0.0_r8
      location%which_vert = 0
      return
   endif
endif

write(*, *)'Vertical coordinate options'
write(*, *)VERTISUNDEF,' --> vertical coordinate undefined'
write(*, *)VERTISSURFACE,' --> surface'
write(*, *)VERTISLEVEL,' --> model level'
write(*, *)VERTISPRESSURE,' --> pressure'
write(*, *)VERTISHEIGHT,' --> height'

100   read(*, *) location%which_vert
if(location%which_vert == VERTISLEVEL ) then
   write(*, *) 'Vertical coordinate model level'
   read(*, *) location%vloc
else if(location%which_vert == VERTISPRESSURE ) then
   write(*, *) 'Vertical coordinate Pressure (in hPa)'
   read(*, *) location%vloc
   location%vloc = 100.0_r8 * location%vloc
else if(location%which_vert == VERTISHEIGHT ) then
   write(*, *) 'Vertical coordinate height (in meters)'
   read(*, *) location%vloc
else if(location%which_vert == VERTISSURFACE ) then
   write(*, *) 'Vertical coordinate surface height'
   read(*, *) location%vloc
else if(location%which_vert == VERTISUNDEF ) then
   location%vloc = 0.0_r8
else
   write(*, *) 'Wrong choice of which_vert try again between ',VERTISUNDEF, &
               ' and ',VERTISHEIGHT
   go to 100
end if

! FIXME: the other location mods have a random location option.

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
        'Vertical'), 'nc_write_location_atts', 'location:storage_order')
call nc_check(nf90_put_att(ncFileID, VarID, 'units',     &
        'which_vert'), 'nc_write_location_atts', 'location:units')

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

locations = get_location( loc )

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

! Set the value of num_obs in the structure
gc%num = num

end subroutine get_close_obs_init

!----------------------------------------------------------------------------

subroutine get_close_obs_destroy(gc)

type(get_close_type), intent(inout) :: gc

end subroutine get_close_obs_destroy

!----------------------------------------------------------------------------

subroutine get_close_maxdist_init(gc, maxdist, maxdist_list)

type(get_close_type), intent(inout) :: gc
real(r8),             intent(in)    :: maxdist
real(r8), intent(in), optional      :: maxdist_list(:)

! Set the maximum distance in the structure
gc%maxdist = maxdist

end subroutine get_close_maxdist_init

!----------------------------------------------------------------------------

subroutine get_close_obs(gc, base_obs_loc, base_obs_type, obs, obs_kind, &
   num_close, close_ind, dist)

! Default version with no smarts; no need to be smart in 1D
! Kinds are available here if one wanted to do more refined distances.

type(get_close_type), intent(in)  :: gc
type(location_type),  intent(in)  :: base_obs_loc, obs(:)
integer,              intent(in)  :: base_obs_type, obs_kind(:)
integer,              intent(out) :: num_close, close_ind(:)
real(r8), optional,   intent(out) :: dist(:)

integer :: i
real(r8) :: this_dist

! the list of locations in the obs() argument must be the same
! as the list of locations passed into get_close_obs_init(), so
! gc%num and size(obs) better be the same.   if the list changes,
! you have to destroy the old gc and init a new one.
if (size(obs) /= gc%num) then
   write(errstring,*)'obs() array must match one passed to get_close_obs_init()'
   call error_handler(E_ERR, 'get_close_obs', errstring, source, revision, revdate)
endif

! Return list of obs that are within maxdist and their distances
num_close = 0
do i = 1, gc%num
   this_dist = get_dist(base_obs_loc, obs(i), base_obs_type, obs_kind(i))
   if(this_dist <= gc%maxdist) then
      ! Add this ob to the list
      num_close = num_close + 1
      close_ind(num_close) = i
      if (present(dist)) dist(num_close) = this_dist 
   endif
end do

end subroutine get_close_obs

!----------------------------------------------------------------------------

function is_location_in_region(loc, minl, maxl)
 
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

if ((loc%vloc < minl%vloc) .or. (loc%vloc > maxl%vloc)) return
 
is_location_in_region = .true.

end function is_location_in_region

!----------------------------------------------------------------------------

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

!----------------------------------------------------------------------------

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

!----------------------------------------------------------------------------

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

!----------------------------------------------------------------------------

function vert_is_height(loc)
 
! Given a location, return true if vertical coordinate is height, else false

logical                          :: vert_is_height
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == VERTISHEIGHT) then
   vert_is_height = .true.
else
   vert_is_height = .false.
endif

end function vert_is_height

!----------------------------------------------------------------------------

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

function has_vertical_localization()
 
! Always returns false since this type of location doesn't support
! vertical localization. (but it could - it makes sense for this type.)

logical :: has_vertical_localization

if ( .not. module_initialized ) call initialize_module

has_vertical_localization = .false.

end function has_vertical_localization

!----------------------------------------------------------------------------
!> return the vertical location 
function get_vert(loc)

type(location_type), intent(in) :: loc
real(r8) :: get_vert

if ( .not. module_initialized ) call initialize_module

get_vert = loc%vloc

end function get_vert

!----------------------------------------------------------------------------
!> set the vertical location
subroutine set_vert(loc, vloc)

type(location_type), intent(inout) :: loc
real(r8),            intent(in) :: vloc !< vertical location

if ( .not. module_initialized ) call initialize_module

loc%vloc = vloc

end subroutine set_vert

!----------------------------------------------------------------------------
!> set the which vert
subroutine set_which_vert(loc, which_vert)

type(location_type), intent(inout) :: loc
integer,                intent(in) :: which_vert !< vertical coordinate type

if ( .not. module_initialized ) call initialize_module

loc%which_vert = which_vert

end subroutine set_which_vert


!----------------------------------------------------------------------------
! end of location/column/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
