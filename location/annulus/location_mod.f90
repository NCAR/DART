! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module location_mod

! Implements location interfaces for a three dimensional annulus
! with a vertical coordinate based on the models native set of
! discrete levels. The internal representation of the location is 
! currently implemented as radians from 0 to 2 PI for the azimuthal
! direction (longitude-like).  The radial distance is latitude-like,
! and the vertical coordinate is zero at the bottom of the annulus.
!

use      types_mod, only : r8, PI, RAD2DEG, DEG2RAD, MISSING_R8, MISSING_I
use  utilities_mod, only : register_module, error_handler, E_ERR, ascii_file_format, &
                           nc_check, find_namelist_in_file, check_namelist_read, &
                           do_output, do_nml_file, do_nml_term, nmlfileunit, &
                           open_file, close_file, nc_check, is_longitude_between
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

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
          vert_is_surface, has_vertical_localization, VERTISSURFACE, &
          VERTISLEVEL, VERTISHEIGHT, &
          set_vert, get_vert, set_which_vert


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type location_type
   private
   real(r8) :: azm, rad, vloc
   integer  :: which_vert
   ! which_vert determines if the location is by level or by height
   ! -1 ==> obs is on surface
   ! 1 ===> obs is by level
   ! 3 ===> obs is by height
end type location_type

integer, parameter :: VERTISSURFACE  = -1 ! surface
integer, parameter :: VERTISLEVEL    =  1 ! by level
integer, parameter :: VERTISHEIGHT   =  3 ! by height

type get_close_type
   private
   integer  :: num
   real(r8) :: maxdist
end type get_close_type

type(random_seq_type) :: ran_seq
logical :: ran_seq_init = .false.
logical, save :: module_initialized = .false.

integer,              parameter :: LocationDims = 3
character(len = 129), parameter :: LocationName = "loc_annulus"
character(len = 129), parameter :: LocationLName = &
                                   "Annulus location: azimuthal angle, radius, and height"

! really just a placeholder.  there was a comment that this code
! needs a namelist with a min & max limit on the radius, but 
! the code no longer has one.  
real(r8) :: min_radius = 0.0_r8
real(r8) :: max_radius = 100000.0_r8

! FIXME: if we set mins & maxes, then we have to enforce the values
! in various places in the code.  none of that is there at this point.
namelist /location_nml/ min_radius, max_radius

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
 
! read namelist, set up anything that needs one-time initialization

integer :: iunit, io

! only do this code once
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

! copy code from threed sphere module for handing the
! options in the vertical?  e.g. distances?

end subroutine initialize_module

!----------------------------------------------------------------------------

function get_dist(loc1, loc2, type1, kind2)

! Compute distance between 2 locations.  Right now the distance only
! depends on the horizontal.  A namelist option might need to be added
! that supports computing a true 3d distance.

type(location_type), intent(in) :: loc1, loc2
integer, optional,   intent(in) :: type1, kind2
real(r8)                        :: get_dist

real(r8) :: x1, y1, x2, y2

if ( .not. module_initialized ) call initialize_module

! FIXME: this does not take into account any vertical separation
! convert from cylindrical to cartesian coordinates
x1 = loc1%rad * cos(loc1%azm)
y1 = loc1%rad * sin(loc1%azm)
x2 = loc2%rad * cos(loc2%azm)
y2 = loc2%rad * sin(loc2%azm)

! Returns distance in m
get_dist = sqrt((x1 - x2)**2 + (y1 - y2)**2)

end function get_dist

!---------------------------------------------------------------------------

function loc_eq(loc1,loc2)
 
! Interface operator used to compare two locations.
! Returns true only if all components are 'the same' to within machine
! precision.

type(location_type), intent(in) :: loc1, loc2
logical                         :: loc_eq

if ( .not. module_initialized ) call initialize_module

loc_eq = .false.

if ( loc1%which_vert /= loc2%which_vert ) return
if ( abs(loc1%azm  - loc2%azm ) > epsilon(loc1%azm ) ) return
if ( abs(loc1%rad  - loc2%rad ) > epsilon(loc1%rad ) ) return
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
 
! Given a location type (where the azimuthal angle is in radians), this 
! routine return the azimuthal angle in degrees, the radius, and the vert 

type(location_type), intent(in) :: loc
real(r8), dimension(3) :: get_location

if ( .not. module_initialized ) call initialize_module

get_location(1) = loc%azm * RAD2DEG                 
get_location(2) = loc%rad
get_location(3) = loc%vloc     

end function get_location

!----------------------------------------------------------------------------

function set_location_single(azm, rad, vert_loc, which_vert)
 
! Puts the given azimuthal angle (in degrees), radius, and vertical 
! location into a location datatype.

real(r8), intent(in) :: azm, rad
real(r8), intent(in) :: vert_loc
integer,  intent(in) :: which_vert
type (location_type) :: set_location_single

if ( .not. module_initialized ) call initialize_module

if(azm < 0.0_r8 .or. azm > 360.0_r8) then
   write(errstring,*)'azimuthal angle (',azm,') is not within range [0,360]'
   call error_handler(E_ERR, 'set_location', errstring, source, revision, revdate)
endif

set_location_single%azm = azm * DEG2RAD
set_location_single%rad = rad 

if(which_vert /= VERTISSURFACE .and. &
   which_vert /= VERTISLEVEL   .and. &
   which_vert /= VERTISHEIGHT  ) then
   write(errstring,*)'which_vert (',which_vert,') must be -1, 1 or 3'
   call error_handler(E_ERR,'set_location', errstring, source, revision, revdate)
endif

set_location_single%which_vert = which_vert
set_location_single%vloc = vert_loc

end function set_location_single

!----------------------------------------------------------------------------

function set_location_array(list)
 
! Location semi-independent interface routine
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

set_location_missing%azm        = MISSING_R8
set_location_missing%rad        = MISSING_R8
set_location_missing%vloc       = MISSING_R8
set_location_missing%which_vert = MISSING_I

end function set_location_missing

!---------------------------------------------------------------------------

function query_location(loc,attr)
 
! Returns the value of the attribute
!

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
   case ('rad','RAD','radius','RADIUS')
      query_location = loc%rad
   case ('azm','AZM','azimuth','AZIMUTH')
      query_location = loc%azm
   case ('vloc','VLOC')
      query_location = loc%vloc
   case default
       call error_handler(E_ERR, 'query_location:', &
          'Only "azm","rad","vloc","which_vert" are legal attributes to request from location', &
          source, revision, revdate)
end select

end function query_location

!----------------------------------------------------------------------------

subroutine write_location(locfile, loc, fform, charstring)
 
! Writes a location to the file.
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
character(len = 128) :: string1

! 10 format(1x,3(f22.14,1x),i4)   ! old
10 format(1X,F21.16,2(1X,G25.16),1X,I2)

if ( .not. module_initialized ) call initialize_module

! writing to a file (normal use) or to a character buffer?
writebuf = present(charstring)

! output file; test for ascii or binary, write what's asked, and return
if (.not. writebuf) then
   if (ascii_file_format(fform)) then
      ! Write out pressure or level along with integer tag
      ! we know azm is between 0, 360, and which_vert is a single digit.
      write(locfile, '(''loc3a'')' ) 
      write(locfile, 10) loc%azm, loc%rad, loc%vloc, loc%which_vert
   else
      write(locfile) loc%azm, loc%rad, loc%vloc, loc%which_vert
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
charlength = 85

if (len(charstring) < charlength) then
   write(errstring, *) 'charstring buffer must be at least ', charlength, ' chars long'
   call error_handler(E_ERR, 'write_location', errstring, source, revision, revdate)
endif

write(string1, '(A,F12.8,1X,G15.8,A)') 'Azm(deg)/Radius(m): ',  &
   loc%azm*RAD2DEG, loc%rad, ' Depth:'

! i am attempting to make these line up so if you have a list of mixed
! vertical units, they all take the same number of columns.  thus the extra
! white space around some of the labels below.
select case  (loc%which_vert)
   case (VERTISSURFACE)
      write(charstring, '(A,1X,G15.6,A)') trim(string1), loc%vloc, ' surface (hPa)'
   case (VERTISLEVEL)
      write(charstring, '(A,1X,F6.0,A)')  trim(string1), loc%vloc, '          level'
   case (VERTISHEIGHT)
      write(charstring, '(A,1X,G15.6,A)') trim(string1), loc%vloc / 1000.0_r8, ' km'
   case default
      write(errstring, *) 'unrecognized key for vertical type: ', loc%which_vert
      call error_handler(E_ERR, 'write_location', errstring, source, revision, revdate)
end select

end subroutine write_location

!----------------------------------------------------------------------------

function read_location(locfile, fform)
 
! Reads location from file that was written by write_location.
! See write_location for additional discussion.

integer, intent(in)                      :: locfile
type(location_type)                      :: read_location
character(len = *), intent(in), optional :: fform

character(len=5) :: header

if ( .not. module_initialized ) call initialize_module

if (ascii_file_format(fform)) then
   read(locfile, '(a5)' ) header

   if(header /= 'loc3a') then
      write(errstring,*)'Expected location header "loc3a" in input file, got ', header
      call error_handler(E_ERR, 'read_location', errstring, source, revision, revdate)
   endif
   ! Now read the location data value
   read(locfile, *) read_location%azm, read_location%rad, &
                    read_location%vloc, read_location%which_vert
else
   read(locfile) read_location%azm, read_location%rad, &
                 read_location%vloc, read_location%which_vert
endif

end function read_location

!--------------------------------------------------------------------------

subroutine interactive_location(location, set_to_default)
 
! Allows for interactive input of a location. Also gives option of selecting
! a uniformly distributed random location.

type(location_type), intent(out) :: location
logical, intent(in), optional    :: set_to_default

real(r8) :: azm, rad, minazm, maxazm, minrad, maxrad

if ( .not. module_initialized ) call initialize_module

! If set_to_default is true, then just zero out and return
if(present(set_to_default)) then
   if(set_to_default) then
      location%azm        = 0.0_r8
      location%rad        = 0.0_r8
      location%vloc       = 0.0_r8
      location%which_vert = 0    ! zero is an invalid vert type
      return
   endif
endif

write(*, *)'Vertical coordinate options'
write(*, *)'-1 -> surface, 1 -> model level, 3 -> depth'

100   read(*, *) location%which_vert
if(location%which_vert == VERTISLEVEL ) then
   write(*, *) 'Vertical coordinate model level'
   read(*, *) location%vloc
else if(location%which_vert == VERTISHEIGHT ) then
   write(*, *) 'Vertical coordinate depth (in negative m)'
   read(*, *) location%vloc
   do while (location%vloc > 0)
      write(*, *) 'Depth must be negative (zero at top of fluid), please try again'
      read(*, *) location%vloc
   end do
else if(location%which_vert == VERTISSURFACE ) then
   write(*, *) 'Vertical coordinate surface pressure (in hPa)'
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

   minazm = -1.0
   do while (minazm < 0.0 .or. minazm > 360.0) 
      write(*, *) 'Input minimum azimuthal angle (0 to 360.0)'
      read(*, *) minazm
      if (minazm < 0.0 .or. minazm > 360.0) then
         write(*, *) 'Angle must be between 0 to 360.0'
      endif
   enddo
   minazm = minazm * DEG2RAD

   maxazm = -1.0
   do while (maxazm < 0.0 .or. maxazm > 360.0) 
      write(*, *) 'Input maximum azimuthal angle (0 to 360.0)'
      read(*, *) maxazm
      if (maxazm < 0.0 .or. maxazm > 360.0) then
         write(*, *) 'Angle must be between 0 to 360.0'
      endif
   enddo
   maxazm = maxazm * DEG2RAD

   ! Azimuth is random from minazm to maxazm, handle wrap around 360.0
   if (minazm > maxazm) maxazm = maxazm + 2.0_r8 * PI
   location%azm = random_uniform(ran_seq) * (maxazm-minazm) + minazm
   if (location%azm > 2.0_r8 * PI) location%azm = location%azm - 2.0_r8 * PI

   minrad = -1.0
   do while (minrad < min_radius) 
      write(*, *) 'Input minimum radius '
      read(*, *) minrad
      if (minrad < min_radius) then
         write(*, *) 'Radius must be larger or equal to ', min_radius
      endif
   enddo

   maxrad = -1.0
   do while (maxrad > max_radius .or. maxrad <= minrad) 
      write(*, *) 'Input maximum radius '
      read(*, *) maxrad
      if (maxrad > max_radius .or. maxrad <= minrad) then
         write(*, *) 'Radius must be greater than minrad and less than ', max_radius
      endif
   enddo

   ! Radius must be area weighted to obtain proper random realizations
   location%rad = sqrt(random_uniform(ran_seq)) * (maxrad-minrad) + minrad

   write(*, *) 'random location is ', location%azm / DEG2RAD, &
                                      location%rad 

else

   rad = -1.0
   do while (rad < min_radius .or. rad > max_radius) 
      write(*, *) 'Input radius '
      read(*, *) rad
      if (rad < min_radius .or. rad > max_radius) then
         write(*, *) 'Radius must be between ', min_radius, ' and ', max_radius
      endif
   enddo

   location%rad = rad
   location%azm = azm*DEG2RAD

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
        'Azimuth Radius Vertical'), 'nc_write_location_atts', 'location:storage_order')
call nc_check(nf90_put_att(ncFileID, VarID, 'units',     &
        'degrees meters which_vert'), 'nc_write_location_atts', 'location:units')

! Define the ancillary vertical array and attributes

call nc_check(nf90_def_var(ncid=ncFileID, name='which_vert', xtype=nf90_int, &
          dimids=(/ ObsNumDimID /), varid=VarID), &
            'nc_write_location_atts', 'which_vert:def_var')

call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', 'vertical coordinate system code'), &
           'nc_write_location_atts', 'which_vert:long_name')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISSURFACE', VERTISSURFACE), &
           'nc_write_location_atts', 'which_vert:VERTISSURFACE')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISLEVEL', VERTISLEVEL), &
           'nc_write_location_atts', 'which_vert:VERTISLEVEL')
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

! Return how many locations are closer than cutoff distance, along with a
! count, and the actual distances if requested.  The kinds are available if
! more sophisticated distance computations are wanted.

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

! use the code in the utils module that knows how to wrap longitude/radians.
if (.not. is_longitude_between(loc%azm, minl%azm, maxl%azm, doradians=.true.)) return
if ((loc%rad < minl%rad) .or. (loc%rad > maxl%rad)) return
if ((loc%vloc < minl%vloc) .or. (loc%vloc > maxl%vloc)) return
 
is_location_in_region = .true.

end function is_location_in_region

!----------------------------------------------------------------------------

function vert_is_undef(loc)
 
! Stub, always returns false.

logical                          :: vert_is_undef
type(location_type), intent(in)  :: loc

vert_is_undef = .false.

end function vert_is_undef

!----------------------------------------------------------------------------

function vert_is_surface(loc)
 
! Given a location, return true if vertical coordinate is surface, else false.

logical                          :: vert_is_surface
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == VERTISSURFACE ) then
   vert_is_surface = .true.
else
   vert_is_surface = .false.
endif

end function vert_is_surface

!----------------------------------------------------------------------------

function vert_is_pressure(loc)
 
! Always returns false, as vertical coordinate is never pressure for the annulus.

logical                          :: vert_is_pressure
type(location_type), intent(in)  :: loc

vert_is_pressure = .false.

end function vert_is_pressure

!----------------------------------------------------------------------------

function vert_is_height(loc)
 
! Given a location, return true if vertical coordinate is height, else false.

logical                          :: vert_is_height
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == VERTISHEIGHT ) then
   vert_is_height = .true.
else
   vert_is_height = .false.
endif

end function vert_is_height

!----------------------------------------------------------------------------

function vert_is_level(loc)
 
! Given a location, return true if vertical coordinate is level, else false.

logical                          :: vert_is_level
type(location_type), intent(in)  :: loc

if ( .not. module_initialized ) call initialize_module

if(loc%which_vert == VERTISLEVEL ) then
   vert_is_level = .true.
else
   vert_is_level = .false.
endif

end function vert_is_level

!---------------------------------------------------------------------------

function has_vertical_localization()
 
! Always returns false since this type of location doesn't support
! vertical localization.

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
! end of location/annulus/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
