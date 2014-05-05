! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: ensemble_manager_mod.f90 6323 2013-07-26 19:32:37Z hkershaw $

!> @brief Aim: to abstract the data structure from the ensemeble manager.
!>
!> The reason for doing this is to enable a model_mod to have knowledge of the
!> location of state_ensemble handle
!> Currently this is not possible because of the dependencies
!> ensemble_manager_mod -> assim_model_mod -> model_mod
!> so you can not have
!>   model_mod -> ensemble_manager
!>
!> I would like to put the indicies into ensemble storage in this module.
!> e.g. ENS_MEAN_COPY etc. 
module data_structure_mod

use types_mod,         only : r8
use time_manager_mod,  only : time_type
use mpi_utilities_mod, only : task_count, datasize, my_task_id

implicit none
private

public ensemble_type, map_pe_to_task, get_var_owner_index, copies_in_window, mean_row

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma/ensemble_manager/ensemble_manager_mod.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 6323 $"
character(len=128), parameter :: revdate  = "$Date: 2013-07-26 13:32:37 -0600 (Fri, 26 Jul 2013) $"

type ensemble_type
   !DIRECT ACCESS INTO STORAGE IS USED TO REDUCE COPYING: BE CAREFUL
   !!!private
   integer                  :: num_copies, num_vars, my_num_copies, my_num_vars
   integer,         pointer :: my_copies(:), my_vars(:)
   ! Storage in next line is to be used when each pe has all copies of subset of vars
   real(r8),        pointer :: copies(:, :)         ! Dimensioned (num_copies, my_num_vars)
   ! Storage on next line is used when each pe has subset of copies of all vars
   real(r8),        pointer :: vars(:, :)           ! Dimensioned (num_vars, my_num_copies)
   ! Time is only related to var complete
   type(time_type), pointer :: time(:)
   integer                  :: distribution_type
   integer                  :: valid     ! copies modified last, vars modified last, both same
   integer                  :: id_num
   integer, allocatable     :: task_to_pe_list(:), pe_to_task_list(:) ! List of tasks
   ! Flexible my_pe, layout_type which allows different task layouts for different ensemble handles
   integer                  :: my_pe    
   integer                  :: layout_type           

end type ensemble_type

!! Module storage for pe information for this process avoids recomputation
! Would need to have an initialize module call to assign this
!integer              :: num_pes

contains

!-----------------------------------------------------------------
!> Given the var number, returns which PE stores it when copy complete
!> and its index in that pes local storage. Depends on distribution_type
!> with only option 1 currently implemented.
!> Assumes that all tasks are used in the ensemble
subroutine get_var_owner_index(var_number, owner, owners_index)

integer, intent(in)  :: var_number !> index into state vector
integer, intent(out) :: owner !> pe who owns the state element
integer, intent(out) :: owners_index !> local index on the owner

integer :: div
integer :: num_pes

num_pes = task_count() !> @todo Fudge task_count()

! Asummes distribution type 1: rount robin
div = (var_number - 1) / num_pes
owner = var_number - div * num_pes - 1
owners_index = div + 1

end subroutine get_var_owner_index

!--------------------------------------------------------------------------------
!> Return the physical task for my_pe
function map_pe_to_task(ens_handle, p)

type(ensemble_type), intent(in) :: ens_handle
integer,             intent(in) :: p  !> pe number
integer                         :: map_pe_to_task !> physical task number

map_pe_to_task = ens_handle%pe_to_task_list(p + 1)

end function map_pe_to_task

!--------------------------------------------------------------------------------
!> return the numbers of copies in the window
function copies_in_window(state_ens_handle)

type(ensemble_type), intent(in) :: state_ens_handle
integer                         :: copies_in_window

copies_in_window = state_ens_handle%num_copies -6

end function copies_in_window

!--------------------------------------------------------------------------------
!> return the index of the mean row
!> mean row is the row in state_ens_handle%copies(:,:) which is the mean. Typically
!> has been state_ens_handle%copies -6 ( just the regular ensemble members
function mean_row(state_ens_handle)

type(ensemble_type), intent(in) :: state_ens_handle
integer                         :: mean_row

mean_row = state_ens_handle%num_copies -5

end function mean_row


end module data_structure_mod
