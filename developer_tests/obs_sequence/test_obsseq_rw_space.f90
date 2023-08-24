program test_obsseq_rw_space

use        types_mod,  only : r8, i8, missing_r8, metadatalength

use    utilities_mod,  only : register_module, initialize_utilities,        &
                              find_namelist_in_file, check_namelist_read,   &
                              error_handler, E_ERR, E_MSG, nmlfileunit,     &
                              do_nml_file, do_nml_term, get_next_filename,  &
                              open_file, close_file, finalize_utilities

use      location_mod, only : location_type, set_location, write_location

use       obs_def_mod, only : obs_def_type, get_obs_def_time,               &
                              get_obs_def_type_of_obs

use      obs_kind_mod, only : max_defined_types_of_obs, get_name_for_type_of_obs

use mpi_utilities_mod, only : initialize_mpi_utilities, my_task_id,         &
                              finalize_mpi_utilities, task_sync

use obs_sequence_mod,  only : obs_sequence_type, obs_type, write_obs_seq,   &
                              read_obs_seq, read_obs_seq_header, init_obs,  &
                              init_obs_sequence, static_init_obs_sequence,  &
                              validate_obs_seq_time, replace_obs_values,    &
                              replace_qc


use ensemble_manager_mod, only : init_ensemble_manager, allocate_single_copy, &
                                 ensemble_type, get_copy, map_pe_to_task, &
                                 all_copies_to_all_vars

use           filter_mod, only : filter_setup_obs_sequence

implicit none

type(obs_sequence_type) :: seq_in
integer                 :: num_copies_in, num_qc_in
integer                 :: num_obs_in, max_num_obs
integer                 :: file_id
character(len=128)      :: read_format
logical                 :: pre_I_format, cal
character(len=256)      :: msgstring, msgstring1, msgstring2

! variable declarations for for obs_space_diagnostics
type(ensemble_type)     :: obs_fwd_op_ens_handle, qc_ens_handle
integer,    allocatable :: keys(:)
integer                 :: last_used_key, key_bounds(2)
integer                 :: in_obs_copy, obs_val_index
integer                 :: ens_size, num_groups
integer                 :: total_obs_copies, num_obs_in_set
integer                 :: num_output_obs_members, members_index
integer                 :: DART_qc_index, input_qc_index
integer                 :: prior_obs_mean_index, prior_obs_spread_index
logical                 :: do_post



integer                 :: OBS_KEY_COPY, OBS_MEAN_START, OBS_VAR_START
integer                 :: OBS_GLOBAL_QC_COPY, OBS_VAL_COPY
integer                 :: OBS_ERR_VAR_COPY

integer :: i

integer :: iunit, io

! integer :: TOTAL_OBS_COPIES = ens_size + 5 + 2*num_groups
! test hard code num of obs
!num_obs_in_set = 1091025
!total_obs_copies = ens_size + 5 + 2*num_groups


!=======================================================
! namelist input default values

character(len=256) :: file_in = ''
character(len=256) :: file_out = ''

namelist /test_obsseq_rw_space_nml/ &
         file_in, file_out

ens_size = 20
num_groups = 1
num_obs_in_set = 1091025 ! use num_obs_in from read_obs_seq
total_obs_copies = ens_size + 5 + 2*num_groups
members_index = 1

num_output_obs_members = ens_size
prior_obs_mean_index   = ens_size + 1
prior_obs_spread_index = ens_size + 2
do_post = .false.

OBS_KEY_COPY         = ens_size + 3
OBS_MEAN_START       = ens_size + 6
OBS_GLOBAL_QC_COPY   = ens_size + 4
OBS_VAL_COPY         = ens_size + 2
OBS_ERR_VAR_COPY     = ens_size + 1


!=======================================================
! main executable

call initialize_mpi_utilities('test_obsseq_rw_space')

! read namelist entries
call find_namelist_in_file("input.nml", "test_obsseq_rw_space_nml", iunit)
read(iunit, nml = test_obsseq_rw_space_nml, iostat = io)
call check_namelist_read(iunit, io, "test_obsseq_rw_space_nml")

! record namelist used
!if (do_nml_file()) write(nmlfileunit, nml=test_obsseq_rw_nml)
!if (do_nml_file()) write(     *     , nml=test_obsseq_rw_nml)

! write(*,*) 'hello from task ', my_task_id(), 'file name: ', trim(file_in)

call read_obs_seq_header(file_in, num_copies_in, num_qc_in, &
      num_obs_in, max_num_obs, file_id, read_format, pre_I_format, &
      close_the_file = .true.)

call read_obs_seq(file_in, 0, 0, 0, seq_in)


call init_ensemble_manager(obs_fwd_op_ens_handle, TOTAL_OBS_COPIES, &
                           int(num_obs_in_set, i8), 1, transpose_type_in = 2)

call init_ensemble_manager(qc_ens_handle, ens_size, &
                           int(num_obs_in_set, i8), 1, transpose_type_in = 2)  
! allocate and populate keys
! test obs_seq has all keys starting from 1 - num. of obs.
allocate(keys(num_obs_in))

do i = 1, num_obs_in
   keys(i) = i
end do

! call filter_setup_obs_sequence from filter_mod
call filter_setup_obs_sequence(seq_in, in_obs_copy, obs_val_index, &
                               input_qc_index, DART_qc_index, do_post)


! call obs_space_diagnostics
call obs_space_diagnostics(obs_fwd_op_ens_handle, qc_ens_handle, ens_size, &
           seq_in, keys, num_output_obs_members, in_obs_copy+1, &
           obs_val_index, OBS_KEY_COPY, &
           prior_obs_mean_index, prior_obs_spread_index, num_obs_in_set, &
           OBS_MEAN_START, OBS_VAR_START, OBS_GLOBAL_QC_COPY, &
           OBS_VAL_COPY, OBS_ERR_VAR_COPY, DART_qc_index, do_post)



if (my_task_id() == 0) then
   call write_obs_seq(seq_in, file_out)
endif

call finalize_mpi_utilities()

! end program test_obsseq_rw_space

!=======================================================
contains

subroutine obs_space_diagnostics(obs_fwd_op_ens_handle, qc_ens_handle, ens_size, &
   seq, keys, num_output_members, members_index, &
   obs_val_index, OBS_KEY_COPY, &
   ens_mean_index, ens_spread_index, num_obs_in_set, &
   OBS_MEAN_START, OBS_VAR_START, OBS_GLOBAL_QC_COPY, OBS_VAL_COPY, &
   OBS_ERR_VAR_COPY, DART_qc_index, do_post)


! Do observation space diagnostics on the set of obs corresponding to keys

type(ensemble_type),     intent(inout) :: obs_fwd_op_ens_handle, qc_ens_handle
integer,                 intent(in)    :: ens_size
integer,                 intent(in)    :: num_obs_in_set
integer,                 intent(in)    :: keys(num_obs_in_set)
integer,                 intent(in)    :: num_output_members, members_index
integer,                 intent(in)    :: obs_val_index
integer,                 intent(in)    :: OBS_KEY_COPY
integer,                 intent(in)    :: ens_mean_index, ens_spread_index
type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(in)    :: OBS_MEAN_START, OBS_VAR_START
integer,                 intent(in)    :: OBS_GLOBAL_QC_COPY, OBS_VAL_COPY
integer,                 intent(in)    :: OBS_ERR_VAR_COPY, DART_qc_index
logical,                 intent(in)    :: do_post

integer               :: j, k, ens_offset, copy_factor
integer               :: ivalue, io_task, my_task
real(r8), allocatable :: obs_temp(:)
real(r8)              :: rvalue(1)

! Do verbose forward operator output if requested
! if(output_forward_op_errors) call verbose_forward_op_output(qc_ens_handle, prior_post, ens_size, keys)

! this is a query routine to return which task has 
! logical processing element 0 in this ensemble.

io_task = map_pe_to_task(obs_fwd_op_ens_handle, 0)
my_task = my_task_id()

! single value per member if no posterior, else 2
if (do_post) then
   copy_factor = 2
else
   copy_factor = 1
endif

! Make var complete for get_copy() calls below.
! Optimize: Could we use a gather instead of a transpose and get copy?
call all_copies_to_all_vars(obs_fwd_op_ens_handle)

! allocate temp space for sending data only on the task that will
! write the obs_seq.final file
if (my_task == io_task) then
   allocate(obs_temp(num_obs_in_set))
else ! TJH: this change became necessary when using Intel 19.0.5 ...
   allocate(obs_temp(1))
endif

! Update the ensemble mean
call get_copy(io_task, obs_fwd_op_ens_handle, OBS_MEAN_START, obs_temp)
if(my_task == io_task) then
   do j = 1, obs_fwd_op_ens_handle%num_vars
      rvalue(1) = obs_temp(j)
      call replace_obs_values(seq, keys(j), rvalue, ens_mean_index)
     end do
  endif

! Update the ensemble spread
call get_copy(io_task, obs_fwd_op_ens_handle, OBS_VAR_START, obs_temp)
if(my_task == io_task) then
   do j = 1, obs_fwd_op_ens_handle%num_vars
      if (obs_temp(j) /= missing_r8) then
         rvalue(1) = sqrt(obs_temp(j))
      else
         rvalue(1) = obs_temp(j)
      endif
      call replace_obs_values(seq, keys(j), rvalue, ens_spread_index)
   end do
endif

! Update any requested ensemble members
ens_offset = members_index + 2*copy_factor
do k = 1, num_output_members
   call get_copy(io_task, obs_fwd_op_ens_handle, k, obs_temp)
   if(my_task == io_task) then
      ivalue = ens_offset + copy_factor * (k - 1)
      do j = 1, obs_fwd_op_ens_handle%num_vars
         rvalue(1) = obs_temp(j)
         call replace_obs_values(seq, keys(j), rvalue, ivalue)
      end do
   endif
end do

! Update the qc global value
call get_copy(io_task, obs_fwd_op_ens_handle, OBS_GLOBAL_QC_COPY, obs_temp)
if(my_task == io_task) then
   do j = 1, obs_fwd_op_ens_handle%num_vars
      rvalue(1) = obs_temp(j)
      call replace_qc(seq, keys(j), rvalue, DART_qc_index)
   end do
endif

deallocate(obs_temp)

end subroutine obs_space_diagnostics

end program test_obsseq_rw_space
