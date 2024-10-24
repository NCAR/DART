program test_window

use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities, my_task_id, task_count, task_sync

use ensemble_manager_mod, only : init_ensemble_manager, end_ensemble_manager, ensemble_type, set_num_extra_copies 
use distributed_state_mod, only : create_state_window, free_state_window, get_state
use types_mod, only : i8, r8

use test ! fortran-testanything

implicit none

integer :: num_copies = 10
integer :: real_ens_members = 3
real(r8) :: res(3)
integer(i8) :: num_vars = 201
type(ensemble_type) :: ens_handle

call initialize_mpi_utilities('test_window')

if (my_task_id() == 0 ) then
  call plan(3*task_count())
endif 

call init_ensemble_manager(ens_handle, num_copies, num_vars)
call set_num_extra_copies(ens_handle, num_copies - real_ens_members)

ens_handle%copies(1:real_ens_members,:) = my_task_id() 
ens_handle%copies(real_ens_members+1:num_copies,:) = -100

call create_state_window(ens_handle)

! result should be index-1 mod task_count() for round robin distribution
res = get_state(1_i8, ens_handle)
call ok(res(1) ==  mod(1-1, task_count()))

res = get_state(27_i8, ens_handle)
call ok(res(1) ==  mod(27-1, task_count()))

res = get_state(198_i8, ens_handle)
call ok(res(1) ==  mod(198-1, task_count()))

call free_state_window(ens_handle)

call end_ensemble_manager(ens_handle)

call finalize_mpi_utilities()


end program test_window