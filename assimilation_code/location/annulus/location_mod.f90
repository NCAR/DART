! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module location_mod

! Implements location interfaces for a three dimensional annulus
! with a vertical coordinate based on the models native set of
! discrete levels. The internal representation of the location is 
! currently implemented as radians from 0 to 2 PI for the azimuthal
! direction (longitude-like).  The radial distance is latitude-like,
! and the vertical coordinate is zero at the bottom of the annulus.
!

use      types_mod, only : r8, PI, RAD2DEG, DEG2RAD, MISSING_R8, MISSING_I
use  utilities_mod, only : error_handler, E_ERR, ascii_file_format, &
                           find_namelist_in_file, check_namelist_read, &
                           do_output, do_nml_file, do_nml_term, nmlfileunit, &
                           open_file, close_file, is_longitude_between
use ensemble_manager_mod, only : ensemble_type
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none
private

public :: location_type, get_location, set_location, &
          set_location_missing, is_location_in_region, get_maxdist, &
          write_location, read_location, interactive_location, query_location, &
          LocationDims, LocationName, LocationLName, LocationStorageOrder, LocationUnits, &
          get_close_type, get_close_init, get_close, get_close_destroy, &
          operator(==), operator(/=), get_dist, &
          vertical_localization_on, is_vertical, set_vertical, &
          get_vertical_localization_coordinate, set_vertical_localization_coordinate, &
          VERTISSURFACE, VERTISLEVEL, VERTISHEIGHT

character(len=*), parameter :: source = 'annulus/location_mod.f90'

type location_type
   private
   real(r8) :: azm, rad, vloc
   integer  :: which_vert
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
character(len=256) :: msgstring

integer,          parameter :: LocationDims = 3
character(len=*), parameter :: LocationName = "loc_annulus"
character(len=*), parameter :: LocationLName = &
                                   "Annulus location: azimuthal angle, radius, and height"
character(len=*), parameter :: LocationStorageOrder = 'Azimuth Radius Vertical'
character(len=*), parameter :: LocationUnits = 'degrees meters which_vert'


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
   call error_handler(E_ERR, 'set_location', errstring, source)
endif

set_location_single%azm = azm * DEG2RAD
set_location_single%rad = rad 

if(which_vert /= VERTISSURFACE .and. &
   which_vert /= VERTISLEVEL   .and. &
   which_vert /= VERTISHEIGHT  ) then
   write(errstring,*)'which_vert (',which_vert,') must be -1, 1 or 3'
   call error_handler(E_ERR,'set_location', errstring, source)
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
   call error_handler(E_ERR, 'set_location', errstring, source)
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
          source)
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
      'Cannot use string buffer with binary format', source)
endif

! format the location to be more human-friendly; meaning
! degrees instead of radians, and kilometers for height,
! hectopascals instead of pascals for pressure, etc.

! this must be the sum of the longest of the formats below.
charlength = 85

if (len(charstring) < charlength) then
   write(errstring, *) 'charstring buffer must be at least ', charlength, ' chars long'
   call error_handler(E_ERR, 'write_location', errstring, source)
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
      call error_handler(E_ERR, 'write_location', errstring, source)
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
      call error_handler(E_ERR, 'read_location', errstring, source)
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

subroutine get_close_init(gc, num, maxdist, locs, maxdist_list)
 
! Initializes part of get_close accelerator that depends on the particular obs

type(get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
real(r8),             intent(in)    :: maxdist
type(location_type),  intent(in)    :: locs(:)
real(r8), intent(in), optional      :: maxdist_list(:)

! Set the maximum distance in the structure
gc%maxdist = maxdist

! Set the value of num_locs in the structure
gc%num = num

end subroutine get_close_init

!----------------------------------------------------------------------------

subroutine get_close_destroy(gc)

type(get_close_type), intent(inout) :: gc

end subroutine get_close_destroy

!----------------------------------------------------------------------------
!> Return how many locations are closer than cutoff distance, along with a
!> count, and the actual distances if requested.  The base type, the
!> loc_quantity and the ensemble_handle are unused here but are available
!> for higher level code to use if needed.

subroutine get_close(gc, base_loc, base_type, locs, loc_quantities, &
                     num_close, close_ind, dist, ensemble_handle)

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(in)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_quantities(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in) :: ensemble_handle

integer :: i
real(r8) :: this_dist

! the list of locations in the locs() argument must be the same
! as the list of locations passed into get_close_init(), so
! gc%num and size(locs) better be the same.   if the list changes,
! you have to destroy the old gc and init a new one.
if (size(locs) /= gc%num) then
   write(errstring,*)'locs() array must match one passed to get_close_init()'
   call error_handler(E_ERR, 'get_close', errstring, source)
endif

! Return list of obs that are within maxdist and their distances
num_close = 0
do i = 1, gc%num
   this_dist = get_dist(base_loc, locs(i), base_type, loc_quantities(i))
   if(this_dist <= gc%maxdist) then
      ! Add this ob to the list
      num_close = num_close + 1
      close_ind(num_close) = i
      if (present(dist)) dist(num_close) = this_dist 
   endif
end do

end subroutine get_close

!---------------------------------------------------------------------------

function get_maxdist(gc, obs_type)
type(get_close_type), intent(in) :: gc
integer, optional,    intent(in) :: obs_type
real(r8) :: get_maxdist

get_maxdist = gc%maxdist

end function get_maxdist

!----------------------------------------------------------------------------
!> Returns true if the given location is between the other two.

function is_location_in_region(loc, minl, maxl)

type(location_type), intent(in)  :: loc, minl, maxl
logical                          :: is_location_in_region

if ( .not. module_initialized ) call initialize_module

if ((minl%which_vert /= maxl%which_vert) .or. &
    (minl%which_vert /= loc%which_vert)) then
   write(errstring,*)'which_vert (',loc%which_vert,') must be same in all args'
   call error_handler(E_ERR, 'is_location_in_region', errstring, source)
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

!---------------------------------------------------------------------------
!> true, this location type has more than one vertical coordinate

function has_vertical_choice()

logical :: has_vertical_choice

if ( .not. module_initialized ) call initialize_module

has_vertical_choice = .false.

end function has_vertical_choice

!----------------------------------------------------------------------------
!> use a string so caller doesn't have to have access to VERTISxxx values

function is_vertical(loc, which_vert)

logical                          :: is_vertical
type(location_type), intent(in)  :: loc
character(len=*),    intent(in)  :: which_vert

select case  (which_vert)
   case ("SURFACE")
      is_vertical = (VERTISSURFACE == loc%which_vert)
   case ("LEVEL")
      is_vertical = (VERTISLEVEL == loc%which_vert)
   case ("HEIGHT")
      is_vertical = (VERTISHEIGHT == loc%which_vert)
   case default
      write(msgstring, *) 'unrecognized key for vertical type: ', which_vert
      call error_handler(E_ERR, 'is_vertical', msgstring, source)
end select

end function is_vertical

!---------------------------------------------------------------------------
!> Always returns false since this type of location doesn't 
!> currently support vertical localization.

function vertical_localization_on()
 
logical :: vertical_localization_on

if ( .not. module_initialized ) call initialize_module

vertical_localization_on = .false.

end function vertical_localization_on

!----------------------------------------------------------------------------
!> set the vertical location and type
!> ( get it with query_location(loc, 'VLOC' or 'WHICH_VERT') )

subroutine set_vertical(loc, vloc, which_vert)

type(location_type), intent(inout) :: loc
real(r8), optional,  intent(in)    :: vloc       !< vertical location
real(r8), optional,  intent(in)    :: which_vert !< vertical type

if ( .not. module_initialized ) call initialize_module

if (present(vloc)) loc%vloc = vloc
if (present(which_vert)) loc%which_vert = which_vert

end subroutine set_vertical

!---------------------------------------------------------------------------

subroutine set_vertical_localization_coordinate()
 
! does nothing in this module

end subroutine set_vertical_localization_coordinate

!---------------------------------------------------------------------------

function get_vertical_localization_coordinate()
 
integer :: get_vertical_localization_coordinate

! does nothing in this module - error to call?
! set it to an undefined value.
get_vertical_localization_coordinate = 0

end function get_vertical_localization_coordinate


!----------------------------------------------------------------------------
! end of location/annulus/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod

