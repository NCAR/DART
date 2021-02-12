! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module default_location_mod

! For location modules which do NOT have multiple choices for
! the vertical coordinate, this routine can be 'use'd by other
! location modules to pass through routines which aren't asking
! useful questions.  this code shouldn't be called directly by
! anything outside of the location modules.

use            types_mod, only : r8
use ensemble_manager_mod, only : ensemble_type

implicit none
private

public :: has_vertical_choice, vertical_localization_on, &
          get_vertical_localization_coord, set_vertical_localization_coord

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

end module default_location_mod

