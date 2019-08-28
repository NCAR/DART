! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module arestart_pnetcdf_mod

use pnetcdf_utilities_mod, only : pnet_check
use types_mod,        only : r8
use time_manager_mod, only : time_type, set_time, get_time
use mpi_utilities_mod,     only : task_count, my_task_id, datasize

use pnetcdf
use mpi

implicit none
private

public :: aread_state_restart_parallel, &
          aread_state_restart_giant, &
          awrite_state_restart_parallel, &
          awrite_state_restart_giant

contains

!----------------------------------------------------------------------
!> reads a netcdf restart file in parallel. No complete state.
subroutine aread_state_restart_parallel(stateId, timeId, nblocks, model_time, model_state, funit, target_time)

implicit none

integer,         intent(in)               :: nblocks !! number of blocks - my_num_vars
type(time_type), intent(out)              :: model_time
real(r8),        intent(out)              :: model_state(:) !! now distributed
integer,         intent(in)               :: funit !! pnetcdf file identifier
type(time_type), optional, intent(out)    :: target_time

character(len = 16) :: open_format
integer :: ios, int1, int2

! Parallel netcdf variables
integer(KIND=MPI_OFFSET_KIND) :: start(1) !! state is one dimensional
integer(KIND=MPI_OFFSET_KIND) :: count(1)
integer(KIND=MPI_OFFSET_KIND) :: stride(1)
integer(KIND=MPI_OFFSET_KIND) :: bufcount !! my num vars
integer                       :: stateId !! Id of state variable
integer                       :: timeId !! Id of time variable
integer                       :: ret !! return code
integer                       :: time(2) !! for reading from netcdf file

ios = 0

!! Figure out whether the file is opened FORMATTED or UNFORMATTED
!inquire(funit, FORM=open_format) !HK the file is open already? Yes

start(1) = my_task_id() + 1
count(1) = nblocks
stride(1) = task_count()
bufcount = nblocks

! Strided read of the state the values from the netCDF variable

if ( datasize == MPI_REAL4 ) then ! trying to catch single precision

   ret = nfmpi_get_vars_all(funit, stateId, start, count, stride, model_state, bufcount, MPI_REAL4)
   call pnet_check(ret, 'aread_state_restart', 'reading state')

else ! double precision

   ret = nfmpi_get_vars_all(funit, stateId, start, count, stride, model_state, bufcount, MPI_REAL8)
   call pnet_check(ret, 'aread_state_restart', 'reading state')

endif

start(1) = 1
count(1) = 2
stride(1) = 1

! Read time from the netCDF variable

!>@todo does this need to be a collective call?
ret = nfmpi_get_vars_int_all(funit, timeId, start, count, stride, time)

call pnet_check(ret, 'aread_state_restart', 'reading time')
model_time = set_time(time(1), time(2)) ! seconds, days

end subroutine aread_state_restart_parallel

!----------------------------------------------------------------------
!> giant restart files
subroutine aread_state_restart_giant(state_length, stateId, timeId, nblocksVars, nblocksCopies, model_time, model_state, transpose, funit, target_time)

use mpi
use pnetcdf

integer(KIND=MPI_OFFSET_KIND),  intent(in) :: state_length !! just to make the overloading work?
integer,                        intent(in) :: stateId !! Id for state_vector
integer,                        intent(in) :: timeId !! Id for time
integer(KIND=MPI_OFFSET_KIND),  intent(in) :: nblocksVars, nblocksCopies !! num vars
type(time_type),             intent(inout) :: model_time(:)
real(r8),                    intent(inout) :: model_state(:, :) !! see inout in ensemble manager call
logical,                        intent(in) :: transpose
integer,                        intent(in) :: funit !! restart file
type(time_type), optional,      intent(in) :: target_time

integer :: time(2), seconds, days

! pnetcdf variables
integer                       :: ret !! return code for pnetcdf calls
integer(KIND=MPI_OFFSET_KIND) :: start(1) !! time is one dimensional
integer(KIND=MPI_OFFSET_KIND) :: count(1)
integer(KIND=MPI_OFFSET_KIND) :: stride(1)
integer(KIND=MPI_OFFSET_KIND) :: state_start(2) !! state is two dimensional
integer(KIND=MPI_OFFSET_KIND) :: state_count(2)
integer(KIND=MPI_OFFSET_KIND) :: state_stride(2)
integer :: time_length

time_length = 2

! read the state vector
! the transpose depends on whether your giant restart file is (copies, vars) or (vars, copies)

if (transpose) then

   state_start = (/ my_task_id() + 1, 1 /)
   state_count = (/ nblocksVars,  nblocksCopies  /)
   state_stride = (/ task_count(), 1 /)

else

   state_start = (/ 1, my_task_id() + 1 /)
   state_count = (/ nblocksCopies, nblocksVars /)
   state_stride = (/ 1, task_count() /)

endif

if ( datasize == MPI_REAL4 ) then ! single precision

   ret = nfmpi_get_vars_all(funit, stateId, state_start, state_count, state_stride, model_state, nblocksCopies*nblocksVars, MPI_REAL4)
   call pnet_check(ret, 'aread_state_restart', 'reading state')

else ! double precision

   ret = nfmpi_get_vars_all(funit, stateId, state_start, state_count, state_stride, model_state, nblocksCopies*nblocksVars, MPI_REAL8)
   call pnet_check(ret, 'aread_state_restart', 'reading state')

endif

! time - does this need to be a collective call? Are the times ever different
start(1) = 1
count(1) = 2
stride(1) = 1

ret = nfmpi_get_vars_int_all(funit, timeId, start, count, stride, time)
call pnet_check(ret, 'aread_state_restart', 'reading time')
model_time(1:nblocksCopies) = set_time(time(1), time(2)) ! seconds, days - one time?

end subroutine aread_state_restart_giant

!----------------------------------------------------------------------
!> parallel write of restarts using pnetcdf
subroutine awrite_state_restart_parallel(timeId, stateId, nblocks, model_time, model_state, state_length, funit, target_time)
!----------------------------------------------------------------------
!
!HK what does this comment mean?
! Write a restart file given a model extended state and a unit number 
! opened to the restart file. (Need to reconsider what is passed to 
! identify file or if file can even be opened within this routine).

use mpi_utilities_mod,     only : task_count, my_task_id, datasize
use pnetcdf_utilities_mod, only :pnet_check
use mpi
use pnetcdf

integer,                        intent(in) :: stateId !! Id for state_vector
integer,                        intent(in) :: timeId !! Id for time
integer(KIND=MPI_OFFSET_KIND),  intent(in) :: nblocks !! num vars
type(time_type),                intent(in) :: model_time
real(r8),                    intent(inout) :: model_state(:) !! see inout in ensemble manager call
integer(KIND=MPI_OFFSET_KIND),  intent(in) :: state_length !! num vars
integer,                        intent(in) :: funit !! restart file
type(time_type), optional,      intent(in) :: target_time

integer :: i, io, rc
character(len = 16) :: open_format
character(len=128) :: filename
logical :: is_named

integer :: time(2), seconds, days

! pnetcdf variables
integer                       :: ret !! return code for pnetcdf calls
integer(KIND=MPI_OFFSET_KIND) :: start(1) !! state and time are one dimensional
integer(KIND=MPI_OFFSET_KIND) :: count(1)
integer(KIND=MPI_OFFSET_KIND) :: stride(1)
integer :: time_length

time_length = 2

! Write the state vector

! state
start(1) = my_task_id() + 1
count(1) = nblocks ! my_num_vars
stride(1) = task_count()

if ( datasize == MPI_REAL4 ) then ! trying to catch single precision

   ret = nfmpi_put_vars_all(funit, stateId, start, count, stride, model_state, nblocks, MPI_REAL4)
   call pnet_check(ret, 'awrite_state_restart', 'writing state')

else ! double precision

   ret = nfmpi_put_vars_all(funit, stateId, start, count, stride, model_state, nblocks, MPI_REAL8)
   call pnet_check(ret, 'awrite_state_restart', 'writing state')

endif

! time - does this need to be a collective call? Are the times ever different
start(1) = 1
count(1) = 2
stride(1) = 1

call get_time(model_time, seconds, days)
time = (/seconds, days/)
ret = nfmpi_put_vars_int_all(funit, timeId, start, count, stride, time)
call pnet_check(ret, 'awrite_state_restart', 'writing time')

end subroutine awrite_state_restart_parallel

!----------------------------------------------------------------------
!> for giant restart file
subroutine awrite_state_restart_giant(timeId, stateId, nblocksVars, nblocksCopies, model_time, model_state, transpose, funit, target_time)
!----------------------------------------------------------------------
!
!HK what does this comment mean?
! Write a restart file given a model extended state and a unit number 
! opened to the restart file. (Need to reconsider what is passed to 
! identify file or if file can even be opened within this routine).

integer,                        intent(in) :: stateId !! Id for state_vector
integer,                        intent(in) :: timeId !! Id for time
integer(KIND=MPI_OFFSET_KIND),  intent(in) :: nblocksVars, nblocksCopies !! num vars
type(time_type),                intent(in) :: model_time !! one time?
real(r8),                    intent(inout) :: model_state(:, :) !! see inout in ensemble manager call
logical,                        intent(in) :: transpose !! pnetcdf stride and count depends on (copies, vars), (vars, copies)
integer,                        intent(in) :: funit !! restart file
type(time_type), optional,      intent(in) :: target_time

integer :: i, io, rc
character(len = 16) :: open_format
character(len=128) :: filename
logical :: is_named

integer :: time(2), seconds, days

! pnetcdf variables
integer                       :: ret !! return code for pnetcdf calls
integer(KIND=MPI_OFFSET_KIND) :: start(1) ! state and time are one dimensional
integer(KIND=MPI_OFFSET_KIND) :: count(1)
integer(KIND=MPI_OFFSET_KIND) :: stride(1)
integer(KIND=MPI_OFFSET_KIND) :: state_start(2) ! state and time are one dimensional
integer(KIND=MPI_OFFSET_KIND) :: state_count(2)
integer(KIND=MPI_OFFSET_KIND) :: state_stride(2)

! time - does this need to be a collective call? Are the times ever different
start(1) = 1
count(1) = 2
stride(1) = 1

call get_time(model_time, seconds, days)
time = (/seconds, days/)
ret = nfmpi_put_vars_int_all(funit, timeId, start, count, stride, time)
call pnet_check(ret, 'awrite_state_restart', 'writing time')


! Write the state vector
! the transpose depends on whether your giant restart file is (copies, vars) or (vars, copies)

if (transpose) then

   state_start = (/ my_task_id() + 1, 1 /)
   state_count = (/ nblocksVars,  nblocksCopies  /)
   state_stride = (/ task_count(), 1 /)

else

   state_start = (/ 1, my_task_id() + 1 /)
   state_count = (/ nblocksCopies, nblocksVars /) ! my_num_vars
   state_stride = (/ 1, task_count() /)

endif

if ( datasize == MPI_REAL4 ) then ! trying to catch single precision

   ret = nfmpi_put_vars_all(funit, stateId, state_start, state_count, state_stride, model_state, nblocksCopies*nblocksVars, MPI_REAL4)
   call pnet_check(ret, 'awrite_state_restart', 'writing state')

else ! double precision

   ret = nfmpi_put_vars_all(funit, stateId, state_start, state_count, state_stride, model_state, nblocksCopies*nblocksVars, MPI_REAL8)
   call pnet_check(ret, 'awrite_state_restart', 'writing state')

endif

end subroutine awrite_state_restart_giant

!----------------------------------------------------------------------

end module arestart_pnetcdf_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
