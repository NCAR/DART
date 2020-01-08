! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> Aim: to be an addition to smoother_mod that uses pnetcdf to write the diagnostic files.
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
query_pnetcdf = .true.

end function

!-----------------------------------------------------------
!> Use pnetcdf to write out the copies array from each processor to a [vars x copies]
!> diagnostic file.
!> It is easy to create a netcdf variable that is bigger than the allowable size of 2GB
!> so, you need to split up the state array into variables less than 2GB. No this is not true,
!> the last variable can be huge.  You do need to switch on large file support if you are 
!> creating the file with netcdf.

subroutine filter_state_space_diagnostics_parallel(state_ens_handle, start_copy, end_copy, diag_filename)

use pnetcdf_utilities_mod, only : pnet_check
use mpi_utilities_mod,     only : my_task_id, task_count, datasize

use mpi
use pnetcdf

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: start_copy !! copy to start from. This is going to be annoying for the posterior output.
integer,             intent(in)    :: end_copy !! copy to output up to
character(len=*),    intent(in)    :: diag_filename

! pnetcdf variables
integer                       :: ret !! pnetcdf return code
integer                       :: ncfile !! pnetcdf file identifier
integer                       :: ndims, dimIds(2), stateId
integer(KIND=MPI_OFFSET_KIND) :: num_copies, num_vars, my_num_vars
integer                       :: copies_dim, vars_dim
integer(KIND=MPI_OFFSET_KIND) :: start(2), count(2), stride(2) !! for state copies
integer(KIND=MPI_OFFSET_KIND) :: bufcount !! my_num_vars * output num_copies
! timing variables
double precision :: start_at_time

start_at_time = MPI_WTIME()

! open diagnostic file
ret = nfmpi_open(mpi_comm_world, diag_filename, NF_WRITE, mpi_info_null, ncfile)
call pnet_check(ret, 'filter_state_space_diagnostics_parallel', 'opening file')

! define variables
ret = nfmpi_redef(ncfile)
call pnet_check(ret, 'filter_state_space_diagnostics_parallel', 'going into define mode')

! define state(copies, vars)
ndims = 2 ! two dimensional state

!>@todo only output what you need to?
num_copies = end_copy - start_copy + 1 

num_vars = state_ens_handle%num_vars
my_num_vars = state_ens_handle%my_num_vars

ret = nfmpi_def_dim(ncfile, 'vars', num_vars, vars_dim) ! length of state vector
call pnet_check(ret, 'filter_state_space_diagnostics_parallel', 'defining vars')

ret = nfmpi_def_dim(ncfile, 'copies', num_copies, copies_dim) ! number of copies
call pnet_check(ret, 'filter_state_space_diagnostics_parallel', 'defining copies')

dimIds = (/copies_dim, vars_dim/)

if (datasize == MPI_REAL4) then ! single precsion state

   ret = nfmpi_def_var(ncfile, 'state', NF_FLOAT, ndims, dimIds, stateId)
   call pnet_check(ret, 'filter_state_space_diagnostics_parallel', 'defining variables')

else ! double precision

   ret = nfmpi_def_var(ncfile, 'state', NF_DOUBLE, ndims, dimIds, stateId)
   call pnet_check(ret, 'filter_state_space_diagnostics_parallel', 'defining variables')

endif

ret = nfmpi_enddef(ncfile) ! metadata IO occurs in this
call pnet_check(ret, 'filter_state_space_diagnostics_parallel', 'metadata IO')

! write state copies to diagnostic file
start = (/1 , my_task_id() + 1/)
count = (/num_copies, my_num_vars/) ! my_num_vars
stride = (/1, task_count()/)
bufcount = state_ens_handle%my_num_vars * num_copies

if (datasize == MPI_REAL4) then ! single precsion state

   ret = nfmpi_put_vars_all(ncfile, stateId, start, count, stride, state_ens_handle%copies(start_copy:end_copy, :), bufcount, MPI_REAL4)
   call pnet_check(ret, 'filter_state_space_diagnostics_parallel', 'writing copies')

else ! double precision

   ret = nfmpi_put_vars_all(ncfile, stateId, start, count, stride, state_ens_handle%copies(start_copy:end_copy, :), bufcount, MPI_REAl8)
   call pnet_check(ret, 'filter_state_space_diagnostics_parallel', 'writing copies')

endif

! close netcdf file
ret = nfmpi_close(ncfile)
call pnet_check(ret, 'filter_state_space_diagnostics_parallel', 'closing file')

if (my_task_id() == 0) print*, 'parallel diagnostic time :', MPI_WTIME() - start_at_time

end subroutine filter_state_space_diagnostics_parallel

!-----------------------------------------------------------




end module smoother_pnetcdf_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
