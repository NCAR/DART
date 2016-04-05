! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: test_read_write_restarts.f90 9553 2016-01-20 17:26:41Z hendric $

program test_read_write_restarts

use        types_mod, only : r8, i8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG,       &
                             open_file, close_file, nc_check, get_next_filename, &
                             find_namelist_in_file, check_namelist_read,         &
                             do_nml_file, do_nml_term, nmlfileunit,              &
                             initialize_utilities, finalize_utilities
use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities,       &
                              task_sync, my_task_id
use assim_model_mod,      only : static_init_assim_model, get_model_size
use state_vector_io_mod,  only : read_state, write_state
use ensemble_manager_mod, only : init_ensemble_manager, ensemble_type
use io_filenames_mod,     only : io_filenames_init, file_info_type
use time_manager_mod,     only : time_type
use filter_mod,           only : filter_set_initial_time

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_test_read_write_restarts_dir/io/test/test_read_write_restarts.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 9553 $"
character(len=128), parameter :: revdate  = "$Date: 2016-01-20 10:26:41 -0700 (Wed, 20 Jan 2016) $"

logical, save :: module_initialized = .false.

integer :: iunit, io

type(ensemble_type)  :: ens_handle
type(file_info_type) :: file_info
type(time_type)      :: time1

logical :: single_restart_file_in = .false.
logical :: single_restart_file_out = .false.
logical :: use_restart_list = .false.
logical :: output_restart = .true.
logical :: output_restart_mean = .false.
logical :: add_domain_extension = .true.
logical :: overwrite_state_input = .false.
logical :: read_time_from_file
logical :: perturb_from_single_instance = .false.

character(len=512) :: restart_list_file(10) = 'null'
character(len=129) :: inf_in_file_name(2) = 'not_initialized'
character(len=129) :: inf_out_file_name(2) = 'not_initialized'

integer(i8) :: model_size
integer     :: ens_size = 3

character(len=129) :: restart_in_file_name
character(len=129) :: restart_out_file_name
logical :: direct_netcdf_read = .false.
logical :: direct_netcdf_write = .false.

! namelist items we are going to create/overwrite
namelist /test_read_write_restarts_nml/  restart_in_file_name, restart_out_file_name, &
direct_netcdf_read, direct_netcdf_write, restart_list_file, use_restart_list

! main code here
 
! initialize the dart libs
call initialize_module()

! Read back the namelist entry
call find_namelist_in_file("input.nml", "test_read_write_restarts_nml", iunit)
read(iunit, nml = test_read_write_restarts_nml, iostat = io)
call check_namelist_read(iunit, io, "test_read_write_restarts_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=test_read_write_restarts_nml)
if (do_nml_term()) write(     *     , nml=test_read_write_restarts_nml)

model_size = get_model_size()

call init_ensemble_manager(ens_handle, ens_size, model_size)

call filter_set_initial_time(-1,-1,time1,read_time_from_file)

file_info = io_filenames_init(ens_handle, single_restart_file_in, single_restart_file_out, &
              restart_in_file_name, restart_out_file_name, output_restart, direct_netcdf_read, &
              direct_netcdf_write, output_restart_mean, add_domain_extension, use_restart_list, &
              restart_list_file, overwrite_state_input, inf_in_file_name, inf_out_file_name)

call read_state(ens_handle, file_info, read_time_from_file, time1)

call write_state(ens_handle, file_info)

! finalize test_read_write_restarts
call error_handler(E_MSG,'test_read_write_restarts','Finished successfully.',source,revision,revdate)
call finalize_mpi_utilities()

! end of main code


contains

!----------------------------------------------------------------------

subroutine initialize_module

  call initialize_mpi_utilities('test_read_write_restarts')
  call register_module(source, revision, revdate)
  call static_init_assim_model()
  module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------

end program

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_test_read_write_restarts_dir/io/test/test_read_write_restarts.f90 $
! $Id: test_read_write_restarts.f90 9553 2016-01-20 17:26:41Z hendric $
! $Revision: 9553 $
! $Date: 2016-01-20 10:26:41 -0700 (Wed, 20 Jan 2016) $
