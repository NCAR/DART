! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module location_mod

! Implements location interfaces for a two dimensional annulus
! discrete levels. The internal representation of the azimuth is 
! implemented as radians and is namelist selectable to be valid 
! from -PI/2 to PI/2 (like latitudes) or from 0 to 2 PI (like longitudes)
! The radial distance is in meters and the inner and outer boundary
! distances are namelist settable.
!

use      types_mod, only : r8, PI, RAD2DEG, DEG2RAD, MISSING_R8, MISSING_I, i8
use  utilities_mod, only : error_handler, E_ERR, ascii_file_format, &
                           find_namelist_in_file, check_namelist_read, &
                           do_output, do_nml_file, do_nml_term, nmlfileunit, &
                           open_file, close_file, is_longitude_between
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform
use ensemble_manager_mod, only : ensemble_type
use default_location_mod, only : has_vertical_choice, vertical_localization_on, &
                                 get_vertical_localization_coord, &
                                 set_vertical_localization_coord, &
                                 is_vertical, set_vertical, &
                                 convert_vertical_obs, convert_vertical_state

implicit none
private

public :: location_type, get_location, set_location, &
          set_location_missing, is_location_in_region, get_maxdist,  &
          write_location, read_location, interactive_location, query_location, &
          LocationDims, LocationName, LocationLName, LocationStorageOrder, LocationUnits, &
          get_close_type, get_close_init, get_close_obs, get_close_state, get_close_destroy, &
          operator(==), operator(/=), get_dist, has_vertical_choice, vertical_localization_on, &
          set_vertical, is_vertical, get_vertical_localization_coord, &
          set_vertical_localization_coord, convert_vertical_obs, convert_vertical_state

character(len=*), parameter :: source   = 'twod_annulus/location_mod.f90'

type location_type
   private
   real(r8) :: azm, rad
   integer  :: which_azm
   ! which_azm determines the valid boundaries of the azimuth
   ! -1 ==> no limits
   !  0 ==> -PI/2 to PI/2 (like latitudes)
   !  1 ==> 0 to 2PI (like longitudes)
end type location_type

integer, parameter :: AZMISUNBOUND  = -1 ! unbounded
integer, parameter :: AZMISLAT      =  0 ! like lats
integer, parameter :: AZMISLON      =  1 ! like lons

type get_close_type
   private
   integer  :: num
   real(r8) :: maxdist
end type get_close_type

type(random_seq_type) :: ran_seq
logical :: ran_seq_init = .false.
logical, save :: module_initialized = .false.

integer,          parameter :: LocationDims = 2
character(len=*), parameter :: LocationName = "loc_2d_annulus"
character(len=*), parameter :: LocationLName = &
                                "2D Annulus location: azimuthal angle, radius"
character(len=*), parameter :: LocationStorageOrder = "X Y"
character(len=*), parameter :: LocationUnits = "degrees none"

! limits on the radius (sets inner and outer radius values)
real(r8) :: min_radius = 0.0_r8
real(r8) :: max_radius = 100000.0_r8

! FIXME: if we set mins & maxes, then we have to enforce the values
! in various places in the code.  none of that is there at this point.
namelist /location_nml/ min_radius, max_radius

character(len = 512) :: errstring

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
! distances?

end subroutine initialize_module

!----------------------------------------------------------------------------

function get_dist(loc1, loc2, type1, kind2)

! Compute distance between 2 locations.

type(location_type), intent(in) :: loc1, loc2
integer, optional,   intent(in) :: type1, kind2
real(r8)                        :: get_dist

real(r8) :: x1, y1, x2, y2

if ( .not. module_initialized ) call initialize_module

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

if ( loc1%which_azm /= loc2%which_azm ) return
if ( abs(loc1%azm  - loc2%azm ) > epsilon(loc1%azm ) ) return
if ( abs(loc1%rad  - loc2%rad ) > epsilon(loc1%rad ) ) return

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
! routine return the azimuthal angle in degrees, and the radius.

type(location_type), intent(in) :: loc
real(r8), dimension(2) :: get_location

if ( .not. module_initialized ) call initialize_module

get_location(1) = loc%azm * RAD2DEG                 
get_location(2) = loc%rad

end function get_location

!----------------------------------------------------------------------------

function set_location_single(azm, rad, which_azm)
 
! Puts the given azimuthal angle (in degrees), and radius
! location into a location datatype.

real(r8), intent(in) :: azm, rad
integer,  intent(in) :: which_azm
type (location_type) :: set_location_single

if ( .not. module_initialized ) call initialize_module

! FIXME: test range based on which_azm
if(azm < 0.0_r8 .or. azm > 360.0_r8) then
   write(errstring,*)'azimuthal angle (',azm,') is not within range [0,360]'
   call error_handler(E_ERR, 'set_location', errstring, source)
endif

set_location_single%azm = azm * DEG2RAD
set_location_single%rad = rad 

set_location_single%which_azm = which_azm

end function set_location_single

!----------------------------------------------------------------------------

function set_location_array(list)
 
! Location semi-independent interface routine
! given 3 float numbers, call the underlying set_location routine

real(r8), intent(in) :: list(:)
type (location_type) :: set_location_array

if ( .not. module_initialized ) call initialize_module

if (size(list) < 3) then
   write(errstring,*)'requires 3 input values'
   call error_handler(E_ERR, 'set_location', errstring, source)
endif

set_location_array = set_location_single(list(1), list(2), nint(list(3)))

end function set_location_array

!----------------------------------------------------------------------------

function set_location_missing()

! Initialize a location type to indicate the contents are unset.

type (location_type) :: set_location_missing

if ( .not. module_initialized ) call initialize_module

set_location_missing%azm        = MISSING_R8
set_location_missing%rad        = MISSING_R8
set_location_missing%which_azm  = MISSING_I

end function set_location_missing

!---------------------------------------------------------------------------

function query_location(loc,attr)
 
! Returns the value of the attribute

type(location_type),        intent(in) :: loc
character(len=*), optional, intent(in) :: attr
real(r8)                               :: query_location

if ( .not. module_initialized ) call initialize_module

! see the long comment in this routine in the threed_sphere
! module for warnings about compiler bugs before you change
! this code.

query_location = loc%which_azm

if (.not. present(attr)) return

select case(attr)
   case ('which_azm','WHICH_AZM')
      query_location = loc%which_azm
   case ('rad','RAD','radius','RADIUS')
      query_location = loc%rad
   case ('azm','AZM','azimuth','AZIMUTH')
      query_location = loc%azm
   case default
       call error_handler(E_ERR, 'query_location:', &
          'Only "azm","rad","which_azm" are legal attributes to request from location', &
          source)
end select

end function query_location

!----------------------------------------------------------------------------

subroutine write_location(locfile, loc, fform, charstring)
 
! Writes a location to a file.
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

10 format(1X,F21.16,1X,G25.16,1X,I2)

if ( .not. module_initialized ) call initialize_module

! writing to a file (normal use) or to a character buffer?
writebuf = present(charstring)

! output file; test for ascii or binary, write what's asked, and return
if (.not. writebuf) then
   if (ascii_file_format(fform)) then
      ! azm is between -90, 360, and which_vert is a single digit.
      write(locfile, '(''loc2a'')' ) 
      write(locfile, 10) loc%azm, loc%rad, loc%which_azm
   else
      write(locfile) loc%azm, loc%rad, loc%which_azm
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
! degrees instead of radians

! this must be the sum of the longest of the formats below.
charlength = 48

if (len(charstring) < charlength) then
   write(errstring, *) 'charstring buffer must be at least ', charlength, ' chars long'
   call error_handler(E_ERR, 'write_location', errstring, source)
endif

write(charstring, '(A,F12.8,1X,G15.8)') 'Azm(deg)/Radius(m): ',  &
   loc%azm*RAD2DEG, loc%rad

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

   if(header /= 'loc2a') then
      write(errstring,*)'Expected location header "loc2a" in input file, got ', header
      call error_handler(E_ERR, 'read_location', errstring, source)
   endif
   ! Now read the location data value
   read(locfile, *) read_location%azm, read_location%rad, &
                    read_location%which_azm
else
   read(locfile) read_location%azm, read_location%rad, &
                 read_location%which_azm
endif

end function read_location

!--------------------------------------------------------------------------

subroutine interactive_location(location, set_to_default)
 
! Allows for interactive input of a location. Also gives option of selecting
! a uniformly distributed random location.

type(location_type), intent(out) :: location
logical, intent(in), optional    :: set_to_default

real(r8) :: azm, rad, minazm, maxazm, minrad, maxrad, minv, maxv
integer :: r

if ( .not. module_initialized ) call initialize_module

! If set_to_default is true, then just zero out and return
if(present(set_to_default)) then
   if(set_to_default) then
      location%azm        = 0.0_r8
      location%rad        = min_radius
      location%which_azm = -1   ! unlimited angle
      return
   endif
endif

write(*, *)'Azimuth coordinate limit options'
write(*, *)'-1 -> none, 0 -> -90, 90, 1 -> 0, 360'

100   read(*, *) location%which_azm
if(location%which_azm == AZMISUNBOUND ) then
   minazm = -HUGE(r8)
   maxazm =  HUGE(r8)
else if(location%which_azm == AZMISLAT ) then
   minazm = -PI/2
   maxazm =  PI/2
else if(location%which_azm == AZMISLON ) then
   minazm = 0.0_r8
   maxazm = 2.0_r8*PI
else
   write(*, *) 'Wrong choice of which_azm try again, valid: -1, 0, and 1'
   go to 100
end if

r = 1
do while (r > 0)

   write(*, *) 'Input 0 to specify a value for the location, or'
   write(*, *) '-1 for a uniformly distributed random location'
   read(*, *) r

   if (r > 0) write(*, *) 'Please input 0 or -1 for selection'
enddo

if (r == 0) then
101 continue
   if (location%which_azm == AZMISUNBOUND) then
      write(*, *) 'Input value for aziumth in degrees '
      read(*,*) azm
   else if (location%which_azm == AZMISLAT) then
      write(*, *) 'Input value for aziumth in degrees (-90 to 90) '
      read(*,*) azm
      if (azm < -90.0_r8 .or. azm > 90.0_r8) then
         write(*,*) 'Illegal value; must be between -90 and 90'
         goto 101
      endif
   else if (location%which_azm == AZMISLON) then
      write(*, *) 'Input value for aziumth in degrees (0 to 360)'
      read(*,*) azm
      if (azm < 0.0_r8 .or. azm > 360.0_r8) then
         write(*,*) 'Illegal value; must be between 0 and 360'
         goto 101
      endif
   endif
   location%azm = azm * DEG2RAD

else
   ! Need to make sure random sequence is initialized

   if(.not. ran_seq_init) then
      call init_random_seq(ran_seq)
      ran_seq_init = .TRUE.
   endif

102 continue
   write(*, *) 'Input minimum azimuth value in degrees '
   read(*, *) minv
   if (location%which_azm == AZMISLAT) then
      if (minv < -90.0_r8) then
         write(*,*) 'Illegal value; minimum must be >= -90'
         goto 102
      endif   
   else if (location%which_azm == AZMISLON) then
      if (minv < 0.0_r8) then
         write(*,*) 'Illegal value; minimum must be >= 0'
         goto 102
      endif   
   endif

103 continue
   write(*, *) 'Input maximum azimuth value in degrees '
   read(*, *) maxv
   if (location%which_azm == AZMISLAT) then
      if (maxv > 90.0_r8) then
         write(*,*) 'Illegal value; maximum must be <= 90'
         goto 103
      endif   
   else if (location%which_azm == AZMISLON) then
      if (maxv > 360.0_r8) then
         write(*,*) 'Illegal value; maximum must be <= 360'
         goto 103
      endif   
   endif

   minv = minv * DEG2RAD
   maxv = maxv * DEG2RAD

   ! Azimuth is random from minazm to maxazm, handle wrap around 360.0
   if (location%which_azm == AZMISLON) then
      if (minv > maxv) maxv = maxv + 2.0_r8 * PI
   endif

   location%azm = random_uniform(ran_seq) * (maxv-minv) + minv

   if (location%which_azm == AZMISLON) then
      if (location%azm > 2.0_r8 * PI) location%azm = location%azm - 2.0_r8 * PI
   endif

endif

if (r == 0) then
   rad = -1.0
   do while (rad < min_radius .or. rad > max_radius) 
      write(*, *) 'Input radius '
      read(*, *) rad
      if (rad < min_radius .or. rad > max_radius) then
         write(*, *) 'Radius must be between ', min_radius, ' and ', max_radius
      endif
   enddo

   location%rad = rad

else

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

   write(*, *) 'random location is ', location%azm * RAD2DEG, &
                                      location%rad 
endif

end subroutine interactive_location

!----------------------------------------------------------------------------
! Initializes get_close accelerator - unused in this location module

subroutine get_close_init(gc, num, maxdist, locs, maxdist_list)

type(get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
real(r8),             intent(in)    :: maxdist
type(location_type),  intent(in)    :: locs(:)
real(r8), intent(in), optional      :: maxdist_list(:)

! Set the maximum localization distance
gc%maxdist = maxdist

! Save the value of num_locs
gc%num = num

if (present(maxdist_list)) then
   write(errstring,*)'twod locations does not support different cutoff distances by type'
   call error_handler(E_ERR, 'get_close_init', errstring, source)
endif

end subroutine get_close_init

!----------------------------------------------------------------------------

subroutine get_close_destroy(gc)

type(get_close_type), intent(inout) :: gc

end subroutine get_close_destroy

!----------------------------------------------------------------------------

subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ensemble_handle)

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(in)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:), loc_types(:) 
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ensemble_handle

call get_close(gc, base_loc, base_type, locs, loc_qtys, &
               num_close, close_ind, dist, ensemble_handle)

end subroutine get_close_obs

!----------------------------------------------------------------------------

subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ensemble_handle)

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(in)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:)
integer(i8),                   intent(in)  :: loc_indx(:) 
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ensemble_handle

call get_close(gc, base_loc, base_type, locs, loc_qtys, &
               num_close, close_ind, dist, ensemble_handle)

end subroutine get_close_state

!----------------------------------------------------------------------------

subroutine get_close(gc, base_loc, base_type, locs, loc_qtys, &
                     num_close, close_ind, dist, ensemble_handle)

type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(in)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:)
integer,                       intent(out) :: num_close, close_ind(:)
real(r8),            optional, intent(out) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ensemble_handle

integer :: i
real(r8) :: this_dist

! the list of locations in the locs() argument must be the same
! as the list of locations passed into get_close_init(), so
! gc%num and size(locs) better be the same.   if the list changes
! you have to destroy the old gc and init a new one.
if (size(locs) /= gc%num) then
   write(errstring,*)'locs() array must match one passed to get_close_init()'
   call error_handler(E_ERR, 'get_close', errstring, source)
endif

! Return list of obs that are within maxdist and their distances
num_close = 0
do i = 1, gc%num
   this_dist = get_dist(base_loc, locs(i), base_type, loc_qtys(i))
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

function is_location_in_region(loc, minl, maxl)
 
! Returns true if the given location is between the other two.

logical                          :: is_location_in_region
type(location_type), intent(in)  :: loc, minl, maxl

if ( .not. module_initialized ) call initialize_module

if ((minl%which_azm /= maxl%which_azm) .or. &
    (minl%which_azm /= loc%which_azm)) then
   write(errstring,*)'which_azm (',loc%which_azm,') must be same in all args'
   call error_handler(E_ERR, 'is_location_in_region', errstring, source)
endif

! assume failure and return as soon as we are confirmed right.
! set to success only at the bottom after all tests have passed.
is_location_in_region = .false.

if (loc%which_azm /= AZMISLON) then
   if ((loc%azm < minl%azm) .or. (loc%azm > maxl%azm)) return
else
   ! use the code in the utils module that knows how to wrap longitude/radians.
   if (.not. is_longitude_between(loc%azm, minl%azm, maxl%azm, doradians=.true.)) return
endif
if ((loc%rad < minl%rad) .or. (loc%rad > maxl%rad)) return
 
is_location_in_region = .true.

end function is_location_in_region


!----------------------------------------------------------------------------
! end of location/twod_annulus/location_mod.f90
!----------------------------------------------------------------------------

end module location_mod

