! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

!> \dir coamps_nest Main program contained here
!> \file innov_to_obs_seq.f90 Main program

program innov_to_obs_seq

!> \mainage innov_to_obs_seq Program to convert navdas-format 'innovations' file to DART format
!> @{ \brief convert innovations to observation format
!>

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

use    utilities_mod, only : initialize_utilities,        &
                             finalize_utilities

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"


character(len=*), parameter :: routine = 'innov_to_obs_seq'
integer :: max_obs


call initialize_utilities(routine)
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
call finalize_utilities(routine)

!> @}

end program innov_to_obs_seq

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
