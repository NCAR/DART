! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_def_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

use       types_mod, only : r8, pi
use   utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use    obs_kind_mod, only : obs_kind_type, read_kind, write_kind, set_obs_kind, &
                            IDENTITY_OBSERVATION, set_ncep_obs_kind, get_obs_kind
use    location_mod, only : location_type, read_location, write_location, &
                            get_location, vert_is_level, read_ncep_obs_location
use   obs_model_mod, only : take_obs, interactive_def
use assim_model_mod, only : am_get_close_states=>get_close_states, &
                            am_get_num_close_states=>get_num_close_states, &
                            get_state_meta_data

implicit none
private

public init_obs_def, get_expected_obs, get_error_variance, get_obs_location, &
   get_obs_def_kind, get_num_close_states, get_close_states, &
   set_obs_def_kind, read_obs_def, write_obs_def, obs_def_type, &
   interactive_obs_def, &
   read_ncep_obs_def, get_seq_loc, get_obs_location4, get_obs_kind4

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

! Need overloaded interface for init_obs_def
interface init_obs_def
   module procedure init_obs_def1
   module procedure init_obs_def2
end interface

type obs_def_type
   private
   type(location_type) :: location
   type(obs_kind_type) :: kind
   real(r8)            :: error_variance
   integer :: model_state_index ! indexes into model state for identity obs
end type obs_def_type

logical, save :: module_initialized = .false.


contains


  subroutine initialize_module
!----------------------------------------------------------------------------
! subroutine initialize_module

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module




  function init_obs_def1(location, knd, error_variance)
!----------------------------------------------------------------------------
! function init_obs_def1(location, knd, error_variance)
!
! Constructor for an obs_def that is not identity observation.

implicit none

type(obs_def_type) :: init_obs_def1
type(location_type), intent(in) :: location
type(obs_kind_type), intent(in) :: knd
real(r8), intent(in) :: error_variance

if ( .not. module_initialized ) call initialize_module

init_obs_def1%location       = location
init_obs_def1%kind           = knd
init_obs_def1%error_variance = error_variance
init_obs_def1%model_state_index = -1       ! This is not an identity observation

end function init_obs_def1


!----------------------------------------------------------------------------

function init_obs_def2(ind, error_variance)

! Constructor for an obs_def that is an identity observation.

implicit none

type(obs_def_type) :: init_obs_def2
integer, intent(in) :: ind
real(r8), intent(in) :: error_variance

if ( .not. module_initialized ) call initialize_module

! Having a non-zero index indicates this is an identity observation
init_obs_def2%model_state_index = ind
init_obs_def2%error_variance = error_variance

! Location comes from model
call get_state_meta_data(ind, init_obs_def2%location)

! Define kind as identity
init_obs_def2%kind = set_obs_kind(IDENTITY_OBSERVATION)

end function init_obs_def2

!----------------------------------------------------------------------------

function get_expected_obs(obs_def, state_vector)

! Given an obs_def and a state vector from a model returns the expected value
! of this observation (forward operator, h).

implicit none

real(r8) :: get_expected_obs
type(obs_def_type), intent(in) :: obs_def
real(r8), intent(in) :: state_vector(:)

if ( .not. module_initialized ) call initialize_module

! If this is identity obs, can just return from state vector now
if(obs_def%model_state_index > 0) then
   get_expected_obs = state_vector(obs_def%model_state_index)
else
   ! Need to figure out exactly how to expand this
   get_expected_obs = take_obs(state_vector, obs_def%location, obs_def%kind)
endif

end function get_expected_obs

!----------------------------------------------------------------------------

function get_error_variance(obs_def)

implicit none

real(r8) :: get_error_variance
type(obs_def_type), intent(in) :: obs_def

if ( .not. module_initialized ) call initialize_module

get_error_variance = obs_def%error_variance

end function get_error_variance

!----------------------------------------------------------------------------

function get_obs_location(obs_def)

! Returns observation location.

implicit none

type(location_type) :: get_obs_location
type(obs_def_type), intent(in) :: obs_def

if ( .not. module_initialized ) call initialize_module

get_obs_location = obs_def%location

end function get_obs_location

!----------------------------------------------------------------------------

function get_obs_def_kind(obs_def)

! Returns observation kind

implicit none

type(obs_kind_type) :: get_obs_def_kind
type(obs_def_type), intent(in) :: obs_def

if ( .not. module_initialized ) call initialize_module

get_obs_def_kind = obs_def%kind

end function get_obs_def_kind

!----------------------------------------------------------------------------

function get_num_close_states(obs_def, radius)

! Returns the number of state variables that are within distance radius of this
! obs_def location. This is a function of the class data for the state, not a 
! particular state so no state vector argument is needed. This limits things to
! one model per executable which might need to be generalized far in the future.
! F90 limitations make this difficult.

implicit none

integer :: get_num_close_states
type(obs_def_type), intent(in) :: obs_def
real(r8), intent(in) :: radius

if ( .not. module_initialized ) call initialize_module

! Call to assim_model level which knows how to work with locations
get_num_close_states = am_get_num_close_states(obs_def%location, radius)

end function get_num_close_states

!----------------------------------------------------------------------------

subroutine get_close_states(obs_def, radius, number, close_state_list, dist)

! Returns the indices of those state variables that are within distance radius
! of the location of the obs_def along with the number of these states. In the 
! initial quick implementation, close_state_list is a fixed size real array
! and an error is returned by setting number to -number if the number of close states
! is larger than the array. Eventually may want to clean this up and make it 
! more efficient by allowing a dynamic storage allocation return.

implicit none

type(obs_def_type), intent(in) :: obs_def
real(r8), intent(in) :: radius
integer, intent(out) :: number
integer, intent(out) :: close_state_list(:)
real(r8), intent(out) :: dist(:)

if ( .not. module_initialized ) call initialize_module

! For now, do this in inefficient redundant way; need to make more efficient soon
! NOTE: Could do the error checking on storage in assim_model if desired, probably
! have to do it there anyway.

number = get_num_close_states(obs_def, radius)

! Check for insufficient storage
if(number > size(close_state_list)) then
   number = -1 * number
   close_state_list = -1
else
   call am_get_close_states(obs_def%location, radius, number, &
      close_state_list, dist)
endif

end subroutine get_close_states

!----------------------------------------------------------------------------

function set_obs_location(obs_def, location)

! Sets the location of an obs_def

implicit none

type(obs_def_type) :: set_obs_location
type(obs_def_type), intent(in) :: obs_def
type(location_type), intent(in) :: location

if ( .not. module_initialized ) call initialize_module

set_obs_location = obs_def
set_obs_location%location = location

end function set_obs_location

!----------------------------------------------------------------------------

function set_error_variance(obs_def, error_variance)

! Sets the error variance of an obs_def

implicit none

type(obs_def_type) :: set_error_variance
type(obs_def_type), intent(in) :: obs_def
real(r8), intent(in) :: error_variance

if ( .not. module_initialized ) call initialize_module

set_error_variance = obs_def
set_error_variance%error_variance = error_variance

end function set_error_variance

!----------------------------------------------------------------------------


function set_obs_def_kind(obs_def, kind)

! Sets the kind of an obs_def

implicit none

type(obs_def_type) :: set_obs_def_kind
type(obs_def_type), intent(in) :: obs_def
type(obs_kind_type), intent(in) :: kind

if ( .not. module_initialized ) call initialize_module

set_obs_def_kind = obs_def
set_obs_def_kind%kind = kind

end function set_obs_def_kind

!----------------------------------------------------------------------------

function read_obs_def(file)

! Reads an obs_def from file which is just an integer unit number in the 
! current preliminary implementation.

implicit none

type(obs_def_type) :: read_obs_def
integer, intent(in) :: file

character(len=5) :: header

if ( .not. module_initialized ) call initialize_module

! Begin by reading five character ascii header, then location, kind, error variance, index

! Need to add additional error checks on read
read(file, 11) header
11 format(a5)
if(header /= 'obdef') then
   call error_handler(E_ERR,'read_obs_def', 'Expected location header "obdef" in input file', &
                      source, revision, revdate)
endif

! Read the location, kind and error variance
read_obs_def%location = read_location(file)
read_obs_def%kind = read_kind(file)
read(file, *) read_obs_def%error_variance
read(file, *) read_obs_def%model_state_index

end function read_obs_def

!----------------------------------------------------------------------------

subroutine write_obs_def(file, obs_def)

! Writes an obs_def to file.

implicit none

integer, intent(in) :: file
type(obs_def_type), intent(in) :: obs_def

if ( .not. module_initialized ) call initialize_module

! Write the 5 character identifier
write(file, 11)
11 format('obdef')

! Write out the location, kind and error variance
call write_location(file, obs_def%location)
call write_kind(file, obs_def%kind)
write(file, *) obs_def%error_variance
write(file, *) obs_def%model_state_index

end subroutine write_obs_def



function interactive_obs_def()
!---------------------------------------------------------------------------
!
! Allows interactive creation of an observation

implicit none

type(obs_def_type) :: interactive_obs_def

integer :: ind
real(r8) :: error_variance

if ( .not. module_initialized ) call initialize_module

write(*, *) 'Input error variance for this observation definition '
read(*, *) error_variance

write(*, *) 'Input an integer index if this is identity observation, else -1'
read(*, *) ind

if(ind > 0) then
   interactive_obs_def = init_obs_def2(ind, error_variance) 
else

! Call obs_model interfaces to set up kind and location
   call interactive_def(interactive_obs_def%location, interactive_obs_def%kind)
   interactive_obs_def%error_variance = error_variance
   interactive_obs_def%model_state_index = -1
endif

end function interactive_obs_def


!----------------------------------
function read_ncep_obs_def(obsunit)
!----------------------------------
implicit none

type(obs_def_type) :: read_ncep_obs_def

integer, intent(in) :: obsunit
integer :: obsindex
real (r8) :: var

   if ( .not. module_initialized ) call initialize_module

   ! Read in NCEP observation location, kind and error variance
   call read_ncep_obs_location(read_ncep_obs_def%location, obsunit, obsindex, var)

   read_ncep_obs_def%kind              = set_ncep_obs_kind(obsindex)
   read_ncep_obs_def%error_variance    = var**2
   read_ncep_obs_def%model_state_index = -1

end function read_ncep_obs_def

!======================================
subroutine get_seq_loc(obs_def, obsloc0)
!======================================
implicit none

real(r8), intent(out) :: obsloc0(3)
type(obs_def_type), intent(in) :: obs_def

type(location_type) :: location                           
real(r8) :: lon, lat, level, lon_lat_lev(3), pressure                     

   if ( .not. module_initialized ) call initialize_module

   lon_lat_lev = get_location(obs_def%location)                                       
   lon = lon_lat_lev(1); lat = lon_lat_lev(2);   
                                                    
   if(vert_is_level(obs_def%location)) then                                           
      level = lon_lat_lev(3)                    
   else                      
     pressure = lon_lat_lev(3)       
   endif                                                                                        
   obsloc0(1) = lon       ! degree
   obsloc0(2) = lat       ! degree
   obsloc0(3) = pressure  ! Pascal

end subroutine get_seq_loc

!======================================
subroutine get_obs_location4(obs_def, obsloc0)

implicit none
real(r8), intent(out) :: obsloc0(3)
type(obs_def_type), intent(in) :: obs_def

type(location_type) :: location                           
real(r8) :: lon, lat, level, lon_lat_lev(3), pressure       

   if ( .not. module_initialized ) call initialize_module

   lon_lat_lev = get_location(obs_def%location)                                       
   lon = lon_lat_lev(1); lat = lon_lat_lev(2);   
                                                    
   if(vert_is_level(obs_def%location)) then                                           
      level = lon_lat_lev(3)                    
   else                      
      pressure = lon_lat_lev(3)       
   endif                                                                                        
   obsloc0(1) = lon*pi/180.0_r8    ! degree
   obsloc0(2) = lat*pi/180.0_r8    ! degree
   obsloc0(3) = pressure           ! Pascal
                                       
end subroutine get_obs_location4


subroutine get_obs_kind4(obs_def, obskind0)
implicit none
real(r8), intent(out) :: obskind0
type(obs_def_type), intent(in) :: obs_def
type(obs_kind_type) :: kind                           

   if ( .not. module_initialized ) call initialize_module

   obskind0 = get_obs_kind(obs_def%kind)                                       
                                       
end subroutine get_obs_kind4


end module obs_def_mod
