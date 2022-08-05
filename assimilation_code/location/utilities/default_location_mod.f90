! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module default_location_mod

! For location modules which do NOT have multiple choices for
! the vertical coordinate, this routine can be 'use'd by other
! location modules to pass through routines which aren't asking
! useful questions.  this code shouldn't be called directly by
! anything outside of the location modules.
!
! because this module is called by the location module, it cannot
! use the location type.  drat fortran and being unable to support
! circular references.  see the vertical module for those routines.

use            types_mod, only : r8, i8
use ensemble_manager_mod, only : ensemble_type
use    location_type_mod, only : location_type

implicit none
private

public :: has_vertical_choice, vertical_localization_on, &
          get_vertical_localization_coord, set_vertical_localization_coord, &
          is_vertical, set_vertical, convert_vertical_obs, convert_vertical_state

integer :: location_vertical_localization_coord = 0

logical, save :: module_initialized = .false.

contains

!----------------------------------------------------------------------------

subroutine initialize_module
 
if (module_initialized) return

module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------------

function has_vertical_choice()

logical :: has_vertical_choice
 
if ( .not. module_initialized ) call initialize_module

has_vertical_choice = .false.

end function has_vertical_choice

!----------------------------------------------------------------------------

function get_vertical_localization_coord()

integer :: get_vertical_localization_coord

if ( .not. module_initialized ) call initialize_module

get_vertical_localization_coord = location_vertical_localization_coord

end function get_vertical_localization_coord

!----------------------------------------------------------------------------

subroutine set_vertical_localization_coord(which_vert)
 
integer, intent(in) :: which_vert

if ( .not. module_initialized ) call initialize_module

location_vertical_localization_coord = which_vert

end subroutine set_vertical_localization_coord

!---------------------------------------------------------------------------

function vertical_localization_on()
 
logical :: vertical_localization_on

if ( .not. module_initialized ) call initialize_module

vertical_localization_on = .false.

end function vertical_localization_on

!----------------------------------------------------------------------------

function is_vertical(loc, which_vert)

logical                          :: is_vertical
type(location_type), intent(in)  :: loc
character(len=*),    intent(in)  :: which_vert

is_vertical = .false.

end function is_vertical

!--------------------------------------------------------------------

subroutine set_vertical(loc, vloc, which_vert)

type(location_type), intent(inout) :: loc
real(r8), optional,  intent(in)    :: vloc
integer,  optional,  intent(in)    :: which_vert

end subroutine set_vertical

!--------------------------------------------------------------------

subroutine convert_vertical_obs(ens_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, my_status)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer,             intent(in)    :: loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: my_status(:)

my_status(:) = 0

end subroutine convert_vertical_obs

!--------------------------------------------------------------------

subroutine convert_vertical_state(ens_handle, num, locs, loc_qtys, loc_indx, &
                                  which_vert, istatus)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer(i8),         intent(in)    :: loc_indx(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus

istatus = 0

end subroutine convert_vertical_state

!--------------------------------------------------------------------

end module default_location_mod

