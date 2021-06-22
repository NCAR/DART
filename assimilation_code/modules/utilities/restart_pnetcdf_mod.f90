! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module restart_pnetcdf_mod

use assim_model_mod,       only : awrite_state_restart, aread_state_restart
use time_manager_mod,      only : time_type
use types_mod,             only : r8
use ensemble_manager_mod,  only : ensemble_type
use pnetcdf_utilities_mod, only : pnet_check
use mpi_utilities_mod,     only : datasize, my_task_id

use mpi
use pnetcdf

use assim_model_mod,   only : aread_state_restart

implicit none
private

public :: read_ensemble_restart_parallel, write_ensemble_restart_parallel, query_pnetcdf

contains

!-----------------------------------------------------------
!> Tell ensemble manager whether or not the code has been compiled with pnetcdf.
function query_pnetcdf()

logical :: query_pnetcdf
query_pnetcdf = .true.

end function

!------------------------------------------------------
!> parallel read of netcdf restart file using pnetcdf
subroutine read_ensemble_restart_parallel(ens_handle, file_name, start_copy, end_copy, giant_restart, transpose_giant, init_time)

type(ensemble_type),  intent(inout) :: ens_handle
character(len = *),   intent(in)    :: file_name
integer,              intent(in)    :: start_copy, end_copy
logical,              intent(in)    :: giant_restart
logical,              intent(in)    :: transpose_giant
type(time_type),      intent(in),    optional :: init_time

character(len = 4)       :: extension !< for restart file number
integer                  :: i !< for loops
character(len=256)       :: this_file_name !< giant restart filename
real(r8), allocatable    :: model_state(:,:) !< to hold transposed copies
integer                  :: my_num_vars

! Parallel netcdf variables
integer                       :: ret !< return code for pnetcdf calls
integer                       :: ncfile !< ncfile
integer                       :: stateId !< Id for state_vector
integer                       :: timeId !< Id for time
integer                       :: stateDimId !< Id of the dimension of state
integer                       :: timeDimId !< Id of the dimension of time
integer                       :: varsDimId !< for giant read
integer                       :: copiesDimId !< for giant read
integer(KIND=MPI_OFFSET_KIND) :: start(1) ! state is one dimensional
integer(KIND=MPI_OFFSET_KIND) :: count(1)
integer(KIND=MPI_OFFSET_KIND) :: stride(1)
integer(KIND=MPI_OFFSET_KIND) :: state_length, time_length
integer(KIND=MPI_OFFSET_KIND) :: num_blocks
integer(KIND=MPI_OFFSET_KIND) :: num_copies


! timing variables
double precision :: time_at_start

num_blocks = ens_handle%my_num_vars ! Just me

time_at_start = MPI_WTIME()

if (giant_restart) then

   this_file_name = 'filter_ic_giant.nc'

   ! Open netcdf file
   ret = nfmpi_open(mpi_comm_world, this_file_name, NF_NOWRITE, mpi_info_null, ncfile)
   call pnet_check(ret, 'read_ensemble_restart', 'open file')

   ! get id for vars_dim - the length of the state
   ret = nfmpi_inq_dimid(ncfile, 'vars', varsDimId)
   call pnet_check(ret, 'read_ensemble_restart', 'cannot get vars dim')

   ! get num vars
   ret = nfmpi_inq_dimlen(ncfile, varsDimId, state_length)
   call pnet_check(ret, 'read_ensemble_restart', 'get state_length')

   ! get id for copies_dim - the number of copies
   ret = nfmpi_inq_dimid(ncfile, 'copies', copiesDimId)
   call pnet_check(ret, 'read_ensemble_restart', 'cannot get state dim')

   ! get num_copies
   ret = nfmpi_inq_dimlen(ncfile, copiesDimId, num_copies)
   call pnet_check(ret, 'read_ensemble_restart', 'get num_copies')

   ! get id for time_dim - the length of the time (should be 2)
    ret = nfmpi_inq_dimid(ncfile, 'time', timeDimId)
    call pnet_check(ret, 'read_ensemble_restart', 'cannot get time dim')

   ! Get status of varible
   ret = nfmpi_inq_varid(ncfile, 'state_vector', stateId)
   call pnet_check(ret, 'read_ensemble_restart', 'cannot get state_vector id')

   ! Get status of varible
   ret = nfmpi_inq_varid(ncfile, 'time', timeId)
   call pnet_check(ret, 'read_ensemble_restart', 'cannot get time id')

   if (transpose_giant) then

      ! need to read transposed state
      allocate(model_state(num_blocks, num_copies))

      call aread_state_restart(state_length, stateId, timeId, num_blocks, num_copies, ens_handle%time, model_state, transpose_giant, ncfile)

      ens_handle%copies(1:num_copies,:) = transpose(model_state)

      deallocate(model_state)

   else

      call aread_state_restart(state_length, stateId, timeId, num_blocks, num_copies, ens_handle%time, ens_handle%copies(1:num_copies,:),transpose_giant, ncfile)
   endif

   ret = nfmpi_close(ncfile)
   call pnet_check(ret, 'read_ensemble_restart', 'closing file')

   if(present(init_time)) ens_handle%time(1:num_copies) = init_time

   if (my_task_id() == 0) print*, 'giant read time ', MPI_WTIME() - time_at_start, 'copies = ', end_copy - start_copy + 1
else

   ! Loop to read in all ensemble members
   PNETCDF_RESTARTS: do i = start_copy, end_copy

      ! Get global index for my ith ensemble
      write(extension, '(i4.4)') i
      this_file_name = trim(file_name) // '.' // extension // '.nc'
      ! Open netcdf file
      ret = nfmpi_open(mpi_comm_world, this_file_name, NF_NOWRITE, mpi_info_null, ncfile)
      call pnet_check(ret, 'read_ensemble_restart', 'open file')

      ! get id for state_dim - the length of the state
      ret = nfmpi_inq_dimid(ncfile, 'state', stateDimId)
      call pnet_check(ret, 'read_ensemble_restart', 'get state_dim')

      ! get id for time_dim - the length of the time (should be 2)
      ret = nfmpi_inq_dimid(ncfile, 'time', timeDimId)
      call pnet_check(ret, 'read_ensemble_restart', 'get time_dim')

      ! get state length ( a dimension called state )
      ret = nfmpi_inq_dimlen(ncfile, stateDimId, state_length)
      call pnet_check(ret, 'read_ensemble_restart', 'get state length')

      ! get time length ( a dimension called time ) should be 2
      ret = nfmpi_inq_dimlen(ncfile, timeDimId, time_length)
      call pnet_check(ret, 'read_ensemble_restart', 'get time length')

      ! I think the model_mod gets the model size (state_length) separately

      ! Get status of varible - Do you need to do this to get the variable id?
      ret = nfmpi_inq_varid(ncfile, 'state_vector', stateId)
      call pnet_check(ret, 'read_ensemble_restart', 'get state id')

      ! Get status of varible - Do you need to do this to get the variable id?
      ret = nfmpi_inq_varid(ncfile, 'time', timeId)
      call pnet_check(ret, 'read_ensemble_restart', 'get time id')

      ! Read the file directly into storage
      my_num_vars = ens_handle%my_num_vars
      call aread_state_restart(stateId, timeId, my_num_vars, ens_handle%time(i), ens_handle%copies(i, :), ncfile)
      if(present(init_time)) ens_handle%time(i) = init_time

      ! Close the restart file
      ret = nfmpi_close(ncfile)
      call pnet_check(ret, 'read_ensemble_restart', 'close restart file')

   end do PNETCDF_RESTARTS

   if (my_task_id() == 0) print*, 'Read time : ', MPI_WTIME() - time_at_start, 'average :', (MPI_WTIME() - time_at_start) / (end_copy - start_copy + 1)
 
endif

end subroutine read_ensemble_restart_parallel

!------------------------------------------------------
!> write restart files in parallel using pnetcdf
subroutine write_ensemble_restart_parallel(ens_handle, file_name, start_copy, end_copy, giant_restart, transpose_giant, init_time)

type(ensemble_type),  intent(inout) :: ens_handle
character(len = *),   intent(in)    :: file_name
integer,              intent(in)    :: start_copy, end_copy
logical,              intent(in)    :: giant_restart
logical,              intent(in)    :: transpose_giant
type(time_type),      intent(in),    optional :: init_time

character(len = 4)                  :: extension !< for restart file number
integer                             :: i !< for loops
character(len=256)                  :: this_file_name !< giant restart filename


! pnetcdf variables
integer                       :: ret !< return code for pnetcdf calls
integer                       :: ncfile !< ncfile
integer                       :: stateId !< Id for state_vector
integer                       :: timeId !< Id for time
integer                       :: stateDimId !< Id of the dimension of state
integer                       :: timeDimId !< Id of the dimension of time
integer                       :: varsDimId
integer                       :: copiesDimId
integer                       :: stateDim(1) !< Array needed to hold state dimension
integer                       :: timeDim(1) !< Array needed to hold state dimension
integer                       :: varsCopiesDim(2) !< Array needed to hold state dimension
integer                       :: num_dims
integer(KIND=MPI_OFFSET_KIND) :: time_length
integer(KIND=MPI_OFFSET_KIND) :: state_length
integer(KIND=MPI_OFFSET_KIND) :: num_blocks
integer(KIND=MPI_OFFSET_KIND) :: num_copies

real(r8), allocatable :: model_state(:,:)

! timing variables
double precision :: start_at_time

time_length = 2
state_length = ens_handle%num_vars ! Whole state
num_blocks = ens_handle%my_num_vars ! Just me
num_copies = ens_handle%num_copies -6

start_at_time = MPI_WTIME()

if (giant_restart) then

   this_file_name = 'filter_restart_giant.nc'

   ! create netcdf file
   ret = nfmpi_create(mpi_comm_world, this_file_name, NF_CLOBBER, mpi_info_null, ncfile)
   call pnet_check(ret, 'write_ensemble_restart', 'cannot create netcdf file')

   ! load up state and time in netcdf metadata
   ! Define the dimensions
   ret = nfmpi_def_dim(ncfile, 'time', time_length, timeDimId) !time
   call pnet_check(ret, 'write_ensemble_restart', 'defining time dim')

   ret = nfmpi_def_dim(ncfile, 'vars', state_length, varsDimId) !vars
   call pnet_check(ret, 'write_ensemble_restart', 'defining state dim')

   ret = nfmpi_def_dim(ncfile, 'copies', num_copies, copiesDimId) !copies
   call pnet_check(ret, 'write_ensemble_restart', 'defining state dim')

   ! The dimids array is used to pass the IDs of the dimensions of the variables
   if ( transpose_giant) then
         varsCopiesDim = (/varsDimId, copiesDimId/)
   else
      varsCopiesDim = (/copiesDimId, varsDimId/) !- these are backwards?
   endif
   timeDim = (/timeDimId/)

   ! Define the variables - state needs to go last so it can be unlimited size
   num_dims = 1
   ret = nfmpi_def_var(ncfile, "time", NF_INT, num_dims, timeDim, timeId) ! time
      call pnet_check(ret, 'write_ensemble_restart', 'defining time var')

   num_dims = 2
   if (datasize == MPI_REAL4) then ! single precision state

      ret = nfmpi_def_var(ncfile, "state_vector", NF_FLOAT, num_dims, varsCopiesDim, stateId) ! state
      call pnet_check(ret, 'write_ensemble_restart', 'defining state var')

   else ! double precision

      ret = nfmpi_def_var(ncfile, "state_vector", NF_DOUBLE, num_dims, varsCopiesDim, stateId) ! state
   call pnet_check(ret, 'write_ensemble_restart', 'defining state var')

   endif

   ret = nfmpi_enddef(ncfile) ! metadata IO occurs in this
   call pnet_check(ret, 'write_ensemble_restart', 'metadata IO')

   if (transpose_giant) then

      allocate(model_state(num_blocks, num_copies))
      model_state = transpose(ens_handle%copies(1:num_copies,:))

      call awrite_state_restart(timeId, stateId, num_blocks, num_copies, ens_handle%time(1), model_state, transpose_giant, ncfile)

      deallocate(model_state)

   else

      call awrite_state_restart(timeId, stateId, num_blocks, num_copies, ens_handle%time(1), ens_handle%copies(1:num_copies,:), transpose_giant, ncfile)

   endif

   ret = nfmpi_close(ncfile)
   call pnet_check(ret, 'write_ensemble_restart', 'closing file')

   if (my_task_id() == 0) print*, 'giant write time ', MPI_WTIME() - start_at_time, 'copies = ', end_copy - start_copy + 1

else

   PNETCDF_RESTARTS: do i = start_copy, end_copy
      write(extension, '(i4.4)') i
      this_file_name = trim(file_name) // '.' // extension // '.nc'

      ! create netcdf file
      ret = nfmpi_create(mpi_comm_world, this_file_name, NF_CLOBBER, mpi_info_null, ncfile)
      call pnet_check(ret, 'write_ensemble_restart', 'cannot create netcdf file')

      ! load up state and time in netcdf metadata
      ! Define the dimensions
      ret = nfmpi_def_dim(ncfile, 'time', time_length, timeDimId) !time
      call pnet_check(ret, 'write_ensemble_restart', 'defining time dim')

      ret = nfmpi_def_dim(ncfile, 'state', state_length, stateDimId) !state
      call pnet_check(ret, 'write_ensemble_restart', 'defining state dim')

      ! The dimids array is used to pass the IDs of the dimensions of the variables
      stateDim = (/stateDimId/)
      timeDim = (/timeDimId/)

      ! Define the variables - state needs to go last because it is huge

      num_dims = 1
      ret = nfmpi_def_var(ncfile, "time", NF_INT, num_dims, timeDim, timeId) ! time
      call pnet_check(ret, 'write_ensemble_restart', 'defining time var')

      if (datasize == MPI_REAL4) then ! single precision state

         ret = nfmpi_def_var(ncfile, "state_vector", NF_FLOAT, num_dims, stateDim, stateId) ! state
         call pnet_check(ret, 'write_ensemble_restart', 'defining state var')

      else ! double precision

         ret = nfmpi_def_var(ncfile, "state_vector", NF_DOUBLE, num_dims, stateDim, stateId) ! state
         call pnet_check(ret, 'write_ensemble_restart', 'defining state var')

      endif

      ret = nfmpi_enddef(ncfile) ! metadata IO occurs in this
      call pnet_check(ret, 'write_ensemble_restart', 'metadata IO')

      call awrite_state_restart(timeId, stateId, num_blocks, ens_handle%time(i), ens_handle%copies(i, :), state_length,ncfile)

      ret = nfmpi_close(ncfile)
      call pnet_check(ret, 'write_ensemble_restart', 'closing file')

   end do PNETCDF_RESTARTS

   if (my_task_id() == 0) print*, 'write time: ', MPI_WTIME() - start_at_time, 'average :', (MPI_WTIME() - start_at_time) / (end_copy - start_copy + 1)

endif

end subroutine write_ensemble_restart_parallel

!------------------------------------------------------

end module restart_pnetcdf_mod

