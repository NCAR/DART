program test_obsseq_rw

use        types_mod,  only : r8, missing_r8, metadatalength

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
                              validate_obs_seq_time


implicit none

type(obs_sequence_type) :: seq_in
integer                 :: num_copies_in, num_qc_in
integer                 :: num_obs_in, max_num_obs
integer                 :: file_id
character(len=128)      :: read_format
logical                 :: pre_I_format, cal
character(len=256)      :: msgstring, msgstring1, msgstring2

!=======================================================
! namelist input default values

character(len=256) :: file_in = ''
character(len=256) :: file_out = ''

namelist /test_obsseq_rw_nml/ &
         file_in, file_out

! io variables

integer :: iunit, io

!=======================================================
! main executable

call initialize_mpi_utilities('test_obsseq_rw')

! read namelist entries
call find_namelist_in_file("input.nml", "test_obsseq_rw_nml", iunit)
read(iunit, nml = test_obsseq_rw_nml, iostat = io)
call check_namelist_read(iunit, io, "test_obsseq_rw_nml")

! record namelist used
!if (do_nml_file()) write(nmlfileunit, nml=test_obsseq_rw_nml)
!if (do_nml_file()) write(     *     , nml=test_obsseq_rw_nml)

! write(*,*) 'hello from task ', my_task_id(), 'file name: ', trim(file_in)


call read_obs_seq_header(file_in, num_copies_in, num_qc_in, &
      num_obs_in, max_num_obs, file_id, read_format, pre_I_format, &
      close_the_file = .true.)


call read_obs_seq(file_in, 0, 0, 0, seq_in)

call validate_obs_seq_time(seq_in, file_in)

call task_sync()

if (my_task_id() == 0) then
   call write_obs_seq(seq_in, file_out)
endif

call finalize_mpi_utilities()

!=======================================================


end program test_obsseq_rw
