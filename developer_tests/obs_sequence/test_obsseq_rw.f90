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
                              validate_obs_seq_time, add_qc, set_qc_meta_data

use filter_mod,        only : get_blank_qc_index, get_obs_qc_index, get_obs_dartqc_index, &
                              filter_generate_copy_meta_data              

implicit none

type(obs_sequence_type) :: seq
integer                 :: num_copies_in, tnum_qc
integer                 :: num_obs_in, max_num_obs
integer                 :: file_id
character(len=128)      :: read_format
logical                 :: pre_I_format, cal
character(len=256)      :: msgstring, msgstring1, msgstring2
integer                 :: copies_num_inc, qc_num_inc, DART_qc_index, input_qc_index
integer                 :: io_task
character(len=metadatalength) :: no_qc_meta_data = 'No incoming data QC'
character(len=metadatalength) :: dqc_meta_data   = 'DART quality control'
integer :: prior_obs_mean_index, posterior_obs_mean_index, &
           prior_obs_spread_index, posterior_obs_spread_index


!=======================================================
! namelist input default values

character(len=256) :: file_in = ''
character(len=256) :: file_out = ''
logical :: do_post = .true.
integer :: num_output_obs_members = 0

namelist /test_obsseq_rw_nml/ &
         file_in, file_out, do_post, num_output_obs_members

! io variables

integer :: iunit, io

!=======================================================
! main executable
! compare to filter_setup_obs_sequence

call initialize_mpi_utilities('test_obsseq_rw')

! read namelist entries
call find_namelist_in_file("input.nml", "test_obsseq_rw_nml", iunit)
read(iunit, nml = test_obsseq_rw_nml, iostat = io)
call check_namelist_read(iunit, io, "test_obsseq_rw_nml")

call read_obs_seq_header(file_in, num_copies_in, tnum_qc, &
      num_obs_in, max_num_obs, file_id, read_format, pre_I_format, &
      close_the_file = .true.)

io_task = 0

! only the task writing the obs_seq.final file needs space for the
! additional copies/qcs.  for large numbers of individual members
! in the final file this takes quite a bit of memory. 

if (my_task_id() == io_task) then
   ! Determine the number of output obs space fields
   if (do_post) then
      ! 4 is for prior/posterior mean and spread, plus
      ! prior/posterior values for all requested members
      copies_num_inc = 4 + (2 * num_output_obs_members)
   else
      ! 2 is for prior mean and spread, plus
      ! prior values for all requested members
      copies_num_inc = 2 + (1 * num_output_obs_members)
   endif
else
   copies_num_inc = 0
endif


! if there are less than 2 incoming qc fields, we will need
! to make at least 2 (one for the dummy data qc and one for
! the dart qc) on task 0.  other tasks just need 1 for incoming qc.
if (tnum_qc < 2) then
   if (my_task_id() == io_task) then
      qc_num_inc = 2 - tnum_qc
   else
      qc_num_inc = 1 - tnum_qc
   endif
else
   qc_num_inc = 0
endif


! Read in with enough space for diagnostic output values and add'l qc field(s)
! ONLY ADD SPACE ON TASK 0.  everyone else just read in the original obs_seq file.
call read_obs_seq(file_in, copies_num_inc, qc_num_inc, 0, seq)

call filter_generate_copy_meta_data(seq, num_copies_in, &
        prior_obs_mean_index, posterior_obs_mean_index, &
        prior_obs_spread_index, posterior_obs_spread_index, &
        do_post)

! check to be sure that we have an incoming qc field.  if not, look for
! a blank qc field
input_qc_index = get_obs_qc_index(seq)
if (input_qc_index < 0) then
   input_qc_index = get_blank_qc_index(seq)
   if (input_qc_index < 0) then
      ! Need 1 new qc field for dummy incoming qc
      call add_qc(seq, 1)
      input_qc_index = get_blank_qc_index(seq)
      if (input_qc_index < 0) then
         call error_handler(E_ERR,'filter_setup_obs_sequence', &
           'error adding blank qc field to sequence; should not happen')
      endif
   endif
   ! Since we are constructing a dummy QC, label it as such
   call set_qc_meta_data(seq, input_qc_index, no_qc_meta_data)
endif

! check to be sure we either find an existing dart qc field and
! reuse it, or we add a new one. only on task 0.
DART_qc_index = get_obs_dartqc_index(seq)
if (DART_qc_index < 0 .and. my_task_id() == io_task) then
   DART_qc_index = get_blank_qc_index(seq)
   if (DART_qc_index < 0) then
      ! Need 1 new qc field for the DART quality control
      call add_qc(seq, 1)
      DART_qc_index = get_blank_qc_index(seq)
      if (DART_qc_index < 0) then
         call error_handler(E_ERR,'filter_setup_obs_sequence', &
           'error adding blank qc field to sequence; should not happen')
      endif
   endif
   call set_qc_meta_data(seq, DART_qc_index, dqc_meta_data)
endif


call task_sync()

if (my_task_id() == 0) then
   call write_obs_seq(seq, file_out)
endif

call finalize_mpi_utilities()

!=======================================================


end program test_obsseq_rw
