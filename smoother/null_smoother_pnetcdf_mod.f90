!> Aim: provide null module so you can compile without pnetcdf
module smoother_pnetcdf_mod

use ensemble_manager_mod,    only : ensemble_type

implicit none
private

public :: filter_state_space_diagnostics_parallel, query_pnetcdf

contains

!-----------------------------------------------------------
!> Tell filter whether or not the code has been compiled with pnetcdf.
function query_pnetcdf()

logical :: query_pnetcdf
query_pnetcdf = .false.

end function

!-----------------------------------------------------------
!> null version - empty
subroutine filter_state_space_diagnostics_parallel(state_ens_handle, start_copy, end_copy, diag_filename)

type(ensemble_type), intent(inout) :: state_ens_handle 
integer,             intent(in)    :: start_copy
integer,             intent(in)    :: end_copy 
character(len=*),    intent(in)    :: diag_filename

! pnetcdf variables
integer :: ret !< pnetcdf return code
integer :: ncfile !< pnetcdf file identifier
integer :: ndims, dimIds(2), stateId
integer :: num_copies, num_vars, my_num_vars
integer :: copies_dim, vars_dim
integer :: start(2), count(2), stride(2) !< for state copies
integer :: bufcount !< my_num_vars * output num_copies
! timing variables

!-----------------------------------------------------------

end subroutine filter_state_space_diagnostics_parallel



end module smoother_pnetcdf_mod