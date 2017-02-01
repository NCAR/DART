! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module null_vert_convert_mod

use ensemble_manager_mod, only : ensemble_type
use location_mod,       only : location_type, get_close_obs, get_close_type
use types_mod,          only : r8

implicit none

private

public :: query_vert_localization_coord, vert_convert_distrib, get_close_obs_distrib

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

contains

!--------------------------------------------------------------------
!> pass the vertical localization coordinate to assim_tools_mod
function query_vert_localization_coord()

integer :: query_vert_localization_coord

query_vert_localization_coord = 1 ! any old value

end function query_vert_localization_coord

!--------------------------------------------------------------------
!> This is used in the filter_assim. The vertical conversion is done using the 
!> mean state.
!> Calling this is a waste of time
subroutine vert_convert_distrib(state_ens_handle, location, obs_kind, istatus)

type(ensemble_type), intent(in)  :: state_ens_handle
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_kind
integer,             intent(out) :: istatus

istatus = 0

end subroutine vert_convert_distrib

!--------------------------------------------------------------------
subroutine get_close_obs_distrib(gc, base_obs_loc, base_obs_kind, obs_loc, &
                                 obs_kind, num_close, close_ind, dist, state_ens_handle)

type(ensemble_type),         intent(inout)  :: state_ens_handle
type(get_close_type),        intent(in)     :: gc
type(location_type),         intent(inout)  :: base_obs_loc, obs_loc(:)
integer,                     intent(in)     :: base_obs_kind, obs_kind(:)
integer,                     intent(out)    :: num_close, close_ind(:)
real(r8),                    intent(out)    :: dist(:)


call get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                          num_close, close_ind, dist)

end subroutine get_close_obs_distrib


!--------------------------------------------------------------------

end module null_vert_convert_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
