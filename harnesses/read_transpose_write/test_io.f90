!> Test harness for a limited transpose code
program test_io_read_transpose

use filter_mod,           only : filter_read_restart_direct, filter_set_initial_time, &
                                 filter_write_restart_direct
use types_mod,            only : r8, i8
use time_manager_mod,     only : time_type
use mpi_utilities_mod,    only : get_dart_mpi_comm, initialize_mpi_utilities, &
                                 finalize_mpi_utilities, datasize, my_task_id, task_sync
use ensemble_manager_mod, only : ensemble_type, init_ensemble_manager, end_ensemble_manager
use model_mod,            only : get_model_size
use assim_model_mod,      only : static_init_assim_model
use utilities_mod,        only : find_namelist_in_file, check_namelist_read
use io_filenames_mod,     only : set_filenames, io_filenames_init

use state_vector_io_mod
use state_structure_mod,     only : static_init_state_type
use state_structure_mod,     only : static_init_state_type

use pio_transpose_mod

use mpi

implicit none

type(ensemble_type) :: state_ens_handle
type(time_type)     :: time1

integer(i8)         :: model_size ! size of the whole state vector
integer             :: ens_size
logical             :: pio_read_restart  = .false.
logical             :: pio_write_restart = .false.
! timing variables
real(r8)            ::read_time, write_time, start 

character(len = 129) :: inflation_in(2), inflation_out(2)

character(len = 129) :: filewrite

integer             :: num_domains
integer             :: icopy

integer :: io, iunit

namelist /test_io_nml/ ens_size, pio_read_restart, pio_write_restart

! initialize utilities
call initialize_mpi_utilities('test_io')
call io_filenames_init()

! read test_io namelist for ensemble size
call find_namelist_in_file('input.nml', 'test_io_nml', iunit)
read(iunit, nml = test_io_nml, iostat = io)
call check_namelist_read(iunit, io, 'test_io_nml')

! initalize state type for netcdf variables
call static_init_state_type()
  
! intialize model mod to get number of domains and model size
call static_init_assim_model()
call state_vector_io_init()

! Allocate model size storage and ens_size storage for metadata for outputting ensembles
model_size = get_model_size()

! make space for state ensemble
call init_ensemble_manager(state_ens_handle, ens_size, model_size)

call set_filenames(state_ens_handle, ens_size, inflation_in, inflation_out)

! read all of the files
start = MPI_WTIME()
if (pio_read_restart) then
   call pio_read_transpose(state_ens_handle)
else
   call setup_read_write(ens_size)
   call turn_read_copy_on(1,ens_size)
   call filter_set_initial_time(time1)
   call filter_read_restart_direct(state_ens_handle, time1, ens_size)
endif

call task_sync()
read_time = MPI_WTIME()-start
if(my_task_id() == 0) print*, 'read time: ', read_time

if (pio_write_restart) then
  call pio_transpose_write(state_ens_handle)
else
   ! set up arrays for which copies to read/write
   call setup_read_write(ens_size)
   
   ! write restarts
   call turn_write_copy_on(1, ens_size)
   call filter_write_restart_direct(state_ens_handle,0, .true.)
endif

call task_sync()
write_time = MPI_WTIME()-read_time
if(my_task_id() == 0) print*, 'write time: ', write_time

! clean up ensemble and utilities
call end_ensemble_manager(state_ens_handle)
call finalize_mpi_utilities(async=0)

end program test_io_read_transpose
