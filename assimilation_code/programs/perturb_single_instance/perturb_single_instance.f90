! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download


!> This is a utility program that computes an ensemble of restarts
!> using the model_mod pert_model_copies

program perturb_single_instance

use types_mod,            only : r8, i8, obstypelength, MAX_NUM_DOMS, MAX_FILES

use time_manager_mod,     only : time_type, set_time_missing, operator(/=), &
                                 print_time
 
use utilities_mod,        only : find_namelist_in_file,        &
                                 error_handler, nmlfileunit, E_MSG, E_ERR,      &
                                 check_namelist_read, do_nml_file, do_nml_term, &
                                 open_file, close_file, set_multiple_filename_lists

use  location_mod,        only : location_type

use  obs_kind_mod,        only : get_num_quantities, get_index_for_quantity, &
                                 get_name_for_quantity

use  sort_mod,            only : index_sort

use assim_model_mod,      only : static_init_assim_model, get_model_size, &
                                 get_state_meta_data, pert_model_copies

use state_vector_io_mod,  only : read_state, write_state

use io_filenames_mod,     only : file_info_type, io_filenames_init,        &
                                 set_io_copy_flag, set_file_metadata,      &
                                 set_member_file_metadata, file_info_dump, &
                                 stage_metadata_type, get_stage_metadata,  &
                                 get_restart_filename, READ_COPY, WRITE_COPY

use state_structure_mod,  only : get_num_domains

use mpi_utilities_mod,    only : initialize_mpi_utilities, task_count, &
                                 finalize_mpi_utilities, my_task_id,   &
                                 send_sum_to

use ensemble_manager_mod, only : ensemble_type, init_ensemble_manager, compute_copy_mean, &
                                 get_my_num_vars, end_ensemble_manager

implicit none

character(len=*), parameter :: source = 'perturb_single_instance.f90'

!----------------------------------------------------------------
! These variables are namelist-controllable.
!
integer               :: ens_size                       = 1
character(len=256)    :: input_files(MAX_FILES)         = ''
character(len=256)    :: output_file_list(MAX_NUM_DOMS) = ''
character(len=256)    :: output_files(MAX_FILES)        = '' 
real(r8)              :: perturbation_amplitude         = 0.0
logical               :: single_restart_file_in         = .false.

namelist /perturb_single_instance_nml/  &
   ens_size,                &
   input_files,             &
   output_files,            &
   output_file_list,        &
   perturbation_amplitude,  &
   single_restart_file_in

!----------------------------------------------------------------
! Additional global variables
!
type(ensemble_type)             :: ens_handle
character(len=256), allocatable :: file_array_input(:,:)
character(len=256), allocatable :: file_array_output(:,:)
character(len=512)              :: msgstring, msgstring1
character(len=256)              :: my_base, my_desc
integer                         :: idom, imem, iunit, io, i
integer                         :: ndomains
logical                         :: interf_provided
integer(i8)                     :: model_size
type(time_type)                 :: member_time
type(file_info_type)            :: file_info_input, file_info_output
type(stage_metadata_type)       :: input_restart_files
type(stage_metadata_type)       :: output_restart_files

!----------------------------------------------------------------
! program start
!----------------------------------------------------------------

call initialize_mpi_utilities('perturb_single_instance')


! Read the namelist entry and print it
call find_namelist_in_file("input.nml", "perturb_single_instance_nml", iunit)
read(iunit, nml = perturb_single_instance_nml, iostat = io)
call check_namelist_read(iunit, io, "perturb_single_instance_nml")

if (do_nml_file()) write(nmlfileunit, nml=perturb_single_instance_nml)
if (do_nml_term()) write(     *     , nml=perturb_single_instance_nml)

if (single_restart_file_in) then
   write(msgstring,  *) 'single_restart_file_in is not supported.'
   write(msgstring1, *) 'Please contact DART if you would like to use this capability.'
   call error_handler(E_ERR,msgstring,msgstring1)
endif

!----------------------------------------------------------------------
! Calling static_init_assim_model() is required, which also calls
! static_init_model(), so there is no need to explicitly call it.
!----------------------------------------------------------------------

call static_init_assim_model()

!----------------------------------------------------------------------
! initialization code, model size
!----------------------------------------------------------------------

model_size = get_model_size()

!----------------------------------------------------------------------
! Make space that is ensemble handle
!----------------------------------------------------------------------
call init_ensemble_manager(ens_handle, ens_size, model_size)

!----------------------------------------------------------------------
! Allocate space for file arrays.  
! Contains a matrix of files (ncopies x ndomains)
! If perturbing from a single instance the number of 
! input files does not have to be ens_size but rather 
! a single file (or multiple files if more than one domain)
!----------------------------------------------------------------------

ndomains = get_num_domains()

!----------------------------------------------------------------------
! can be ens_size but rather a single file 
! (or multiple files if more than one domain)
!----------------------------------------------------------------------
allocate(file_array_input(ens_size, ndomains))

file_array_input = RESHAPE(input_files,  (/1,  ndomains/))

!----------------------------------------------------------------------
! read in a single ensemble member
!----------------------------------------------------------------------
call io_filenames_init(file_info_input, &
                       ncopies       = 1, &
                       cycling       = single_restart_file_in, &
                       single_file   = single_restart_file_in, &
                       restart_files = file_array_input)

!----------------------------------------------------------------------
! Read the template file to get the shape of netCDF file
! and its variables. It is possible to have multiple domains
! but only require one member.
!----------------------------------------------------------------------
write(my_base,'(A)') 'template'
write(my_desc,'(A)') 'template file'
call set_file_metadata(file_info_input,                  &
                       cnum     = 1,                     &
                       fnames   = file_array_input(1,:), &
                       basename = my_base,               &
                       desc     = my_desc)

call set_io_copy_flag(file_info_input, &
                      cnum    = 1,     &
                      io_flag = READ_COPY)

input_restart_files = get_stage_metadata(file_info_input)

imem = 1
do idom = 1, ndomains
   write(msgstring1, *) '- Reading File : ', &
                     trim(get_restart_filename(input_restart_files, &
                                               copy   = imem,       &
                                               domain = idom))
   call error_handler(E_MSG, 'perturb_single_instance: ', msgstring1, source)
enddo

!----------------------------------------------------------------------
! Read the ensemble from files
!----------------------------------------------------------------------
member_time = set_time_missing()
call read_state(ens_handle, file_info_input, read_time_from_file=.true., model_time=member_time)

!----------------------------------------------------------------------
! Copy from ensemble member 1 to the other copies
!----------------------------------------------------------------------
do i = 1, get_my_num_vars(ens_handle)
   ens_handle%copies(2:ens_size, i) = ens_handle%copies(1, i)
enddo

call pert_model_copies(ens_handle, ens_size, perturbation_amplitude, interf_provided)
if (.not. interf_provided) then
   call error_handler(E_ERR, 'model_mod::pert_model_copies interface required', source)
endif 


!----------------------------------------------------------------------
! can be ens_size but rather a single file 
! (or multiple files if more than one domain)
!----------------------------------------------------------------------
allocate(file_array_output(ens_size, ndomains))

!----------------------------------------------------------------------
! Given either a vector of in/output_files or a text file containing
! a list of files, return a vector of files containing the filenames.
!----------------------------------------------------------------------
call set_multiple_filename_lists(output_files(:), &
                                 output_file_list(:), &
                                 ndomains, &
                                 ens_size, &
                                 'perturb_single_instance', &
                                 'output_files', &
                                 'output_file_list')

file_array_output = RESHAPE(output_files,  (/ens_size,  ndomains/))

!----------------------------------------------------------------------
! output ens_size perturbed restarts
!----------------------------------------------------------------------
call io_filenames_init(file_info_output, &
                       ncopies       = ens_size, &
                       cycling       = single_restart_file_in, &
                       single_file   = single_restart_file_in, &
                       restart_files = file_array_output)

do imem = 1, ens_size
   write(my_base,'(A,I0.2)') 'output_',                 imem
   write(my_desc,'(A,I0.2)') 'output ensemble member ', imem
   call set_file_metadata(file_info_output,                    &
                          cnum     = imem,                     &
                          fnames   = file_array_output(imem,:), &
                          basename = my_base,                  &
                          desc     = my_desc)

   call set_io_copy_flag(file_info_output,      &
                         cnum    = imem,     &
                         io_flag = WRITE_COPY)
enddo

output_restart_files = get_stage_metadata(file_info_output)
do imem = 1, ens_size
   do idom = 1, ndomains
      write(msgstring1, *) '- Writing File : ', imem, idom, &
                        trim(get_restart_filename(output_restart_files, &
                                                  copy   = imem,       &
                                                  domain = idom))
      call error_handler(E_MSG, 'perturb_single_instance: ', msgstring1, source)
   enddo
enddo

call write_state(ens_handle, file_info_output)

!----------------------------------------------------------------------
! clean up allocated memory
!----------------------------------------------------------------------
call end_ensemble_manager(ens_handle)
deallocate(file_array_output, file_array_input)

call finalize_mpi_utilities()

!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------

end program perturb_single_instance

