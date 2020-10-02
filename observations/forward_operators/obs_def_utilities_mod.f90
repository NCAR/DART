! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

!----------------------------------------------------------------------
!> Utilities that can be used by any of the obs_def_xxx_mods.
!> In a separate module to obs_def_mod to avoid a circular dependency.
module obs_def_utilities_mod

use types_mod, only : r8, MISSING_R8

implicit none

private

public :: track_status, set_debug_fwd_op

logical :: debug = .false. ! Flag for debugging

contains

!----------------------------------------------------------------------
!> This routine can be used to set the debugging level of track_status.
!> It is called from filter_mod.
subroutine set_debug_fwd_op(setme)

logical, intent(in) :: setme

debug = setme

end subroutine set_debug_fwd_op

!----------------------------------------------------------------------
!> track_status can be used to keep track of the status of
!> each ensemble member during multiple calls to model_interpolate
!> for a given obs_def.
!> It assumes that you are starting with istatus(:) = 0
!> If debugging, return_now is only set to true if all istatuses are non-zero
!> If not debugging, return_now is set to true if any istatues are non-zero
!>  and any remaining zero istatues are set to 1.
subroutine track_status(ens_size, val_istatus, val_data, istatus, return_now)

integer,  intent(in)    :: ens_size
integer,  intent(in)    :: val_istatus(ens_size)
real(r8), intent(inout) :: val_data(ens_size) !! expected_obs for obs_def
integer,  intent(inout) :: istatus(ens_size) !! istatus for obs_def
logical,  intent(out)   :: return_now

where (istatus == 0) istatus = val_istatus
where (istatus /= 0) val_data = MISSING_R8

return_now = .false.
if (debug) then
   if( all(istatus /= 0))then
      return_now = .true.
      val_data(:) = missing_r8
   endif
else
   if( any(istatus /= 0)) then
      return_now = .true.
      val_data(:) = missing_r8
      where (istatus == 0) istatus = 1
   endif
endif

end subroutine track_status
!----------------------------------------------------------------------

end module obs_def_utilities_mod
