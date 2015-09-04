!> Temporary module to clean up diagnostic file writing.
!> Needs to:
!>   * write existing diagnstic files (using model_mod interface)
!>        - this is temporary
!>
!> Removed making skeleton diagnostic files.  Do we need to make skeleton
!> diagnostic files?
!> Can we write diagnostic files using state_structure_mod?
!>   * The advantage of this is that writing diagnostic files is then
!>     not model specific - can do parallel writes or anything we want
!>   * nc_write_model_atts and nc_write_model_vars are model specific for
!>     the diagnostic files we have now.5
!> For single timestep, do we even want diagnostic files?
!>   * You don't know that you have single step until you are in AdvanceTime
!>     do loop
!>
!> Does the smoother write number of lags diagnostic files?
!>
!> Why is part of the diagnostic code in assim_model_mod? 

module state_space_diag_mod

use        types_mod,     only : r8, i8
use time_manager_mod,     only : time_type, get_time
use ensemble_manager_mod, only : ensemble_type, map_task_to_pe, get_copy, &
                                 all_copies_to_all_vars
use assim_model_mod,      only : netcdf_file_type, aoutput_diagnostics, &
                                 netcdf_file_type
use adaptive_inflate_mod, only : adaptive_inflate_type
use mpi_utilities_mod,    only : my_task_id
use adaptive_inflate_mod, only : do_varying_ss_inflate, do_single_ss_inflate

implicit none

private

public :: filter_state_space_diagnostics

contains
!-----------------------------------------------------------
!> Not sure if anyone will ever want a skeleton diagnostic file.
!> 
!> Skeleton version just to write the time to the diagnostic file
!> This needs to add to the netcdf file (restarts) for multi-step assimilation
!> using aoutput_diagnostics -> nc_get_tindex
subroutine skeleton_filter_state_space_diagnostics(curr_ens_time, out_unit, ens_handle, model_size, &
            num_output_state_members, output_state_mean_index, output_state_spread_index, &
            ENS_MEAN_COPY, ENS_SD_COPY, inflate, INF_COPY, INF_SD_COPY)

type(time_type),             intent(in)    :: curr_ens_time
type(netcdf_file_type),      intent(inout) :: out_unit
type(ensemble_type),         intent(inout) :: ens_handle
integer(i8),                 intent(in)    :: model_size
integer,                     intent(in)    :: num_output_state_members
integer,                     intent(in)    :: output_state_mean_index, output_state_spread_index
! temp_ens is passed from above to avoid extra storage
!real(r8),                    intent(out)   :: temp_ens(:)
type(adaptive_inflate_type), intent(in)    :: inflate
integer,                     intent(in)    :: ENS_MEAN_COPY, ENS_SD_COPY, INF_COPY, INF_SD_COPY

type(time_type) :: temp_time
integer         :: ens_offset, j
real(r8)        :: temp_ens(1) ! junk value

! Assumes that mean and spread have already been computed

! just to write the time
if(my_task_id() == 0) call aoutput_diagnostics(out_unit, curr_ens_time, temp_ens, output_state_mean_index)

end subroutine skeleton_filter_state_space_diagnostics

!-----------------------------------------------------------
!> Full version - this has a model specific call
subroutine filter_state_space_diagnostics(curr_ens_time, out_unit, ens_handle, model_size, &
            num_output_state_members, output_state_mean_index, output_state_spread_index, &
           output_inflation, ENS_MEAN_COPY, ENS_SD_COPY, inflate, INF_COPY, INF_SD_COPY)

type(time_type),             intent(in)    :: curr_ens_time
type(netcdf_file_type),      intent(inout) :: out_unit
type(ensemble_type),         intent(inout) :: ens_handle
integer(i8),                 intent(in)    :: model_size
integer,                     intent(in)    :: num_output_state_members
integer,                     intent(in)    :: output_state_mean_index, output_state_spread_index
! temp_ens is passed from above to avoid extra storage - not any more
! You only do this if there is an ens_mean
!real(r8),                    intent(out)   :: temp_ens(model_size)
type(adaptive_inflate_type), intent(in)    :: inflate
integer,                     intent(in)    :: ENS_MEAN_COPY, ENS_SD_COPY, INF_COPY, INF_SD_COPY
logical,                     intent(in)    :: output_inflation

type(time_type)        :: temp_time
integer                :: ens_offset, j
real(r8), allocatable  :: temp_ens(:) ! junk value

! Assumes that mean and spread have already been computed
allocate(ens_handle%vars(ens_handle%num_vars, ens_handle%my_num_copies))
call all_copies_to_all_vars(ens_handle)

! HK note that for single source we could pass the ens_mean storage
! from above as the trunk does. The trunk passes ens_mean to temp_ens

! task 0 needs some space
if (my_task_id() == 0) then
   allocate(temp_ens(model_size))
else
   allocate(temp_ens(1))
endif

! Output ensemble mean
call get_copy(map_task_to_pe(ens_handle, 0), ens_handle, ENS_MEAN_COPY, temp_ens)
if(my_task_id() == 0) call aoutput_diagnostics(out_unit, curr_ens_time, temp_ens,  &
   output_state_mean_index)

! Output ensemble spread
call get_copy(map_task_to_pe(ens_handle, 0), ens_handle, ENS_SD_COPY, temp_ens) 
if(my_task_id() == 0) call aoutput_diagnostics(out_unit, curr_ens_time, temp_ens, &
   output_state_spread_index)

! Compute the offset for copies of the ensemble
ens_offset = 2

! Output state diagnostics as required: NOTE: Prior has been inflated
do j = 1, num_output_state_members
   ! Get this state copy to task 0; then output it
   call get_copy(map_task_to_pe(ens_handle, 0), ens_handle, j, temp_ens, temp_time)
   if(my_task_id() == 0) call aoutput_diagnostics( out_unit, temp_time, temp_ens, ens_offset + j)
end do

! Unless specifically asked not to, output inflation
if (output_inflation) then
   ! Output the spatially varying inflation if used
   if(do_varying_ss_inflate(inflate) .or. do_single_ss_inflate(inflate)) then
      call get_copy(map_task_to_pe(ens_handle, 0), ens_handle, INF_COPY, temp_ens)
   else
      ! Output inflation value as 1 if not in use (no inflation)
      temp_ens = 1.0_r8
   endif

   if(my_task_id() == 0) call aoutput_diagnostics(out_unit,  curr_ens_time, temp_ens, &
     ens_offset + num_output_state_members + 1)  


   if(do_varying_ss_inflate(inflate) .or. do_single_ss_inflate(inflate)) then
      call get_copy(map_task_to_pe(ens_handle, 0), ens_handle, INF_SD_COPY, temp_ens)
   else
      ! Output inflation sd as 0 if not in use
      temp_ens = 0.0_r8
   endif

   if(my_task_id() == 0) call aoutput_diagnostics(out_unit, curr_ens_time, temp_ens, &
      ens_offset + num_output_state_members + 2) 

endif

deallocate(ens_handle%vars, temp_ens)

end subroutine filter_state_space_diagnostics

!-----------------------------------------------------------

end module state_space_diag_mod