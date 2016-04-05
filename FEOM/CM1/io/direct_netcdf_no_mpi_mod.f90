!> At the moment this is to have a null version of read_transpose
!> and transpose_write for programs like closest_member_tool that
!> read in an state vector (old school) but don't use mpi.  
!> These programs need to be able to use state_vector_io_mod but 
!> without mpi.
!> Not sure if this is the best way to organize the code.
!> Should this module reaad from netcdf and put into a copies array on a single task?
module direct_netcdf_mod

use ensemble_manager_mod, only : ensemble_type

implicit none

private

public :: read_transpose, transpose_write

contains

!-------------------------------------------------
subroutine read_transpose(state_ens_handle, domain, dart_index, read_limit_mem, read_limit_procs)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: domain
integer,             intent(inout) :: dart_index !< This is for mulitple domains
integer,             intent(in)    :: read_limit_mem !< How many state elements you can read at once
integer,             intent(in)    :: read_limit_procs !< How many processors are involved in the transpose



end subroutine read_transpose
!-------------------------------------------------

subroutine transpose_write(state_ens_handle, num_extras, domain, dart_index, isprior, write_limit_mem, write_limit_procs, write_single_precision)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: num_extras ! non restart copies
integer,             intent(in)    :: domain
integer,             intent(inout) :: dart_index
logical,             intent(in)    :: isprior
integer,             intent(in)    :: write_limit_mem !< How many state elements you can read at once
integer,             intent(in)    :: write_limit_procs !< How many processors are involved in the transpose
logical,             intent(in)    :: write_single_precision



end subroutine transpose_write
!-------------------------------------------------

end module direct_netcdf_mod