! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module arestart_pnetcdf_mod

use time_manager_mod, only : time_type
use types_mod,        only : r8

implicit none
private

public :: aread_state_restart_parallel, aread_state_restart_giant, awrite_state_restart_parallel, awrite_state_restart_giant

contains

!----------------------------------------------------------------------
!> null version - does nothing
subroutine aread_state_restart_parallel(stateId, timeId, nblocks, model_time, model_state, funit, target_time)

integer                                   :: stateId 
integer                                   :: timeId 
integer,         intent(in)               :: nblocks
type(time_type), intent(out)              :: model_time
real(r8),        intent(out)              :: model_state(:) 
integer,         intent(in)               :: funit 
type(time_type), optional, intent(out)    :: target_time


end subroutine aread_state_restart_parallel

!----------------------------------------------------------------------
!> null version - does nothing
subroutine aread_state_restart_giant(state_length, stateId, timeId, nblocksVars, nblocksCopies, model_time, model_state, transpose, funit, target_time)

integer,                    intent(in) :: state_length !< just to make the overloading work?
integer,                    intent(in) :: stateId !< Id for state_vector
integer,                    intent(in) :: timeId !< Id for time
integer,                    intent(in) :: nblocksVars, nblocksCopies !< num vars
type(time_type),         intent(inout) :: model_time(:)
real(r8),                intent(inout) :: model_state(:, :) !< see inout in ensemble manager call
logical,                    intent(in) :: transpose
integer,                    intent(in) :: funit !< restart file
type(time_type), optional,  intent(in) :: target_time

end subroutine aread_state_restart_giant

!----------------------------------------------------------------------
!> null version - does nothing
subroutine awrite_state_restart_parallel(timeId, stateId, nblocks, model_time, model_state, state_length, funit, target_time)
!----------------------------------------------------------------------

integer,                   intent(in) :: stateId !< Id for state_vector
integer,                   intent(in) :: timeId !< Id for time
integer,                   intent(in) :: nblocks !< num vars
type(time_type),           intent(in) :: model_time
real(r8),                  intent(inout) :: model_state(:) !< see inout in ensemble manager call
integer,                   intent(in) :: state_length !< num vars
integer,                   intent(in) :: funit !< restart file
type(time_type), optional, intent(in) :: target_time


end subroutine awrite_state_restart_parallel

!----------------------------------------------------------------------
!> null version - does nothing
subroutine awrite_state_restart_giant(timeId, stateId, nblocksVars, nblocksCopies, model_time, model_state, transpose, funit, target_time)
!----------------------------------------------------------------------

integer,                        intent(in) :: stateId !< Id for state_vector
integer,                        intent(in) :: timeId !< Id for time
integer,                        intent(in) :: nblocksVars, nblocksCopies !< num vars
type(time_type),                intent(in) :: model_time !< one time?
real(r8),                    intent(inout) :: model_state(:, :) !< see inout in ensemble manager call
logical,                        intent(in) :: transpose !< pnetcdf stride and count depends on (copies, vars), (vars, copies)
integer,                        intent(in) :: funit !< restart file
type(time_type), optional,      intent(in) :: target_time


end subroutine awrite_state_restart_giant


!----------------------------------------------------------------------

end module arestart_pnetcdf_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
