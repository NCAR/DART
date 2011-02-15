! innov_to_obs_seq
! --------------------
! 
program innov_to_obs_seq

   use navdas_innov_mod, only : seq,                         &
                                num_copies,                  &
                                num_qc,                      &
                                read_ngt_file,               &
                                open_ngt_file,               &
                                open_innov_file,             &
                                read_innov_header,           &
                                read_innov_data,             &
                                get_max_obs,                 &
                                init_navdas_innov_mod,       &
                                terminate_navdas_innov_mod,  &
                                obs_seq_in_name

   use obs_sequence_mod, only : obs_sequence_type,           &
                                init_obs_sequence,           &
                                set_copy_meta_data,          &
                                set_qc_meta_data,            &
                                write_obs_seq

  implicit none

  character(len=128) :: &
       source   = " ",  &
       revision = " ",  &
       revdate  = " "

  character(len=*), parameter :: routine = 'innov_to_obs_seq'

  integer :: max_obs

  call init_navdas_innov_mod()

  call open_innov_file()
  call read_innov_header()

  call open_ngt_file()

  max_obs = get_max_obs()
  
  call init_obs_sequence(seq, num_copies, num_qc, max_obs)
  call set_copy_meta_data(seq, num_copies, 'observation')
  call set_qc_meta_data(seq, num_qc, 'navdas qc')

  call read_ngt_file()

  call read_innov_data()

  call write_obs_seq(seq, obs_seq_in_name)

  call terminate_navdas_innov_mod()

end program innov_to_obs_seq
