! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
      
module filter_io_diag_mod
     
!------------------------------------------------------------------------------

use types_mod,             only : r8, missing_r8

use ensemble_manager_mod,  only : ensemble_type, compute_copy_mean_sd

use random_seq_mod,        only : random_seq_type, init_random_seq, random_gaussian

use io_filenames_mod,      only : file_info_type, set_member_file_metadata, &
                                  set_file_metadata, COPY_NOT_PRESENT, &
                                  io_filenames_init, set_io_copy_flag, &
                                  READ_COPY, WRITE_COPY

use assim_model_mod,       only : pert_model_copies

use state_vector_io_mod,   only : get_stage_to_write, write_state

use adaptive_inflate_mod,  only : do_rtps_inflate, adaptive_inflate_type

use utilities_mod,         only : set_multiple_filename_lists

use state_structure_mod,   only : get_num_domains

!------------------------------------------------------------------------------

! Ensemble copy numbers; Initialized to not be output
integer :: ENS_START_COPY                  = COPY_NOT_PRESENT
integer :: ENS_END_COPY                    = COPY_NOT_PRESENT
integer :: ENS_MEAN_COPY                   = COPY_NOT_PRESENT
integer :: ENS_SD_COPY                     = COPY_NOT_PRESENT
integer :: PRIOR_INF_COPY                  = COPY_NOT_PRESENT
integer :: PRIOR_INF_SD_COPY               = COPY_NOT_PRESENT
integer :: POST_INF_COPY                   = COPY_NOT_PRESENT
integer :: POST_INF_SD_COPY                = COPY_NOT_PRESENT
integer :: RTPS_PRIOR_SPREAD_COPY          = COPY_NOT_PRESENT

! Identifier for different copies for diagnostic files
integer, parameter :: ENS_START     = 1
integer, parameter :: ENS_END       = 2
integer, parameter :: ENS_MEAN      = 3
integer, parameter :: ENS_SD        = 4
integer, parameter :: PRIOR_INF     = 5
integer, parameter :: PRIOR_INF_SD  = 6
integer, parameter :: POST_INF      = 7
integer, parameter :: POST_INF_SD   = 8

! Data structure for different diag and output stages
integer, parameter :: NUM_SCOPIES    = 8
integer :: DIAG_FILE_COPIES( NUM_SCOPIES )     = COPY_NOT_PRESENT
                                  
!------------------------------------------------------------------------------

public ::  create_ensemble_from_single_file, do_stage_output, &
   count_state_ens_copies, init_input_file_info, &
   initialize_file_information,    &
   ENS_MEAN_COPY, ENS_SD_COPY, PRIOR_INF_COPY, &
   PRIOR_INF_SD_COPY, POST_INF_COPY, POST_INF_SD_COPY, RTPS_PRIOR_SPREAD_COPY

!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!> Produces an ensemble by copying my_vars of the 1st ensemble member
!> and then perturbing the copies array.
!> Mimicks the behaviour of pert_model_state:
!> pert_model_copies is called:
!>   if no model perturb is provided, perturb_copies_task_bitwise is called.
!> Note: Not enforcing a model_mod to produce a
!> pert_model_copies that is bitwise across any number of
!> tasks, although there is enough information in the
!> ens_handle to do this.
!>
!> Some models allow missing_r8 in the state vector.  If missing_r8 is
!> allowed the locations of missing_r8s are stored before the perturb,
!> then the missing_r8s are put back in after the perturb.

subroutine create_ensemble_from_single_file(ens_handle, ens_size, &
   perturbation_amplitude, missing_ok)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: ens_size
real(r8),            intent(in)    :: perturbation_amplitude
logical,             intent(in)    :: missing_ok

integer               :: i
logical               :: interf_provided ! model does the perturbing
logical, allocatable  :: miss_me(:)
integer               :: partial_state_on_my_task ! the number of elements ON THIS TASK

! Copy from ensemble member 1 to the other copies
do i = 1, ens_handle%my_num_vars
   ens_handle%copies(2:ens_size, i) = ens_handle%copies(1, i)  ! How slow is this?
enddo

! If the state allows missing values, we have to record their locations
! and restore them in all the new perturbed copies.

if (missing_ok) then
   partial_state_on_my_task = size(ens_handle%copies,2)
   allocate(miss_me(partial_state_on_my_task))
   miss_me = .false.
   where(ens_handle%copies(1, :) == missing_r8) miss_me = .true.
endif

call pert_model_copies(ens_handle, ens_size, perturbation_amplitude, interf_provided)
if (.not. interf_provided) then
   call perturb_copies_task_bitwise(ens_handle, ens_size, perturbation_amplitude)
endif

! Restore the missing_r8
if (missing_ok) then
   do i = 1, ens_size
      where(miss_me) ens_handle%copies(i, :) = missing_r8
   enddo
   deallocate(miss_me)
endif

end subroutine create_ensemble_from_single_file

!------------------------------------------------------------------------------
! Perturb the copies array in a way that is bitwise reproducible
! no matter how many task you run on.

subroutine perturb_copies_task_bitwise(ens_handle, ens_size, perturbation_amplitude)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: ens_size
real(r8),            intent(in)    :: perturbation_amplitude

integer               :: i, j ! loop variables
type(random_seq_type) :: r(ens_size)
real(r8)              :: random_array(ens_size) ! array of random numbers
integer               :: local_index

! Need ens_size random number sequences.
do i = 1, ens_size
   call init_random_seq(r(i), i)
enddo

local_index = 1 ! same across the ensemble

! Only one task is going to update per i.  This will not scale at all.
do i = 1, ens_handle%num_vars

   do j = 1, ens_size
     ! Can use %copies here because the random number
     ! is only relevant to the task than owns element i.
     random_array(j)  =  random_gaussian(r(j), ens_handle%copies(j, local_index), perturbation_amplitude)
   enddo

   if (ens_handle%my_vars(local_index) == i) then
      ens_handle%copies(1:ens_size, local_index) = random_array(:)
      local_index = local_index + 1 ! task is ready for the next random number
      local_index = min(local_index, ens_handle%my_num_vars)
   endif

enddo

end subroutine perturb_copies_task_bitwise

!------------------------------------------------------------------------------

subroutine do_stage_output(stage_name, output_interval, time_step_number, &
   state_ens_handle, file_info, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

character(len = *),   intent(in)    :: stage_name
integer,              intent(in)    :: output_interval, time_step_number
type(ensemble_type),  intent(inout) :: state_ens_handle
type(file_info_type), intent(inout) :: file_info
integer,              intent(in)    :: ens_size, ENS_MEAN_COPY, ENS_SD_COPY

if(get_stage_to_write(stage_name)) then
   if((output_interval > 0) .and. &
      (time_step_number / output_interval * output_interval == time_step_number)) then

      ! Compute the ensemble mean and standard deviation
      ! For efficiency could make this optional in call
      call compute_copy_mean_sd(state_ens_handle, 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

      call write_state(state_ens_handle, file_info)
   endif
endif

end subroutine do_stage_output

!------------------------------------------------------------------------------

! Assign the indices to the storage in the ensemble manager
! Initialize indices for where diagnostic copies are found

subroutine count_state_ens_copies(ens_size, output_mean, output_sd,          &
   do_prior_inflate, do_posterior_inflate, prior_inflate, posterior_inflate, &
   num_copies, num_extras)

integer,                     intent(in)  :: ens_size
logical,                     intent(in)  :: output_mean, output_sd
logical,                     intent(in)  :: do_prior_inflate, do_posterior_inflate
type(adaptive_inflate_type), intent(in)  :: prior_inflate, posterior_inflate
integer,                     intent(out) :: num_copies, num_extras

! First ens_size entries are the actual ensemble members
ENS_START_COPY = 1
ENS_END_COPY   = ens_size

! Filter Extra Copies For Assimilation
ENS_MEAN_COPY     = ens_size + 1
ENS_SD_COPY       = ens_size + 2
PRIOR_INF_COPY    = ens_size + 3
PRIOR_INF_SD_COPY = ens_size + 4
POST_INF_COPY     = ens_size + 5
POST_INF_SD_COPY  = ens_size + 6
num_copies = ens_size + 6

! If Whitaker/Hamill (2012) relaxation-to-prior-spread (RTPS) inflation
! then we need an extra copy to hold (save) the prior ensemble spread
if ( do_rtps_inflate(posterior_inflate) ) then
   RTPS_PRIOR_SPREAD_COPY = ens_size + 7
   num_copies = ens_size + 7
endif                  
   
! Copies in excess of ensemble
num_extras = num_copies - ens_size

! Specifying copies to be output for state diagnostics files
DIAG_FILE_COPIES(ENS_START) = ENS_START_COPY
DIAG_FILE_COPIES(ENS_END) = ENS_END_COPY
if(output_mean) then
   DIAG_FILE_COPIES(ENS_MEAN) = ENS_MEAN_COPY
else
   DIAG_FILE_COPIES(ENS_MEAN) = COPY_NOT_PRESENT
endif
if(output_sd) then
   DIAG_FILE_COPIES(ENS_SD) = ENS_SD_COPY
else
   DIAG_FILE_COPIES(ENS_SD) = COPY_NOT_PRESENT
endif
if(do_prior_inflate) then
   DIAG_FILE_COPIES(PRIOR_INF) = PRIOR_INF_COPY
   DIAG_FILE_COPIES(PRIOR_INF_SD) = PRIOR_INF_SD_COPY
else
   DIAG_FILE_COPIES(PRIOR_INF) = COPY_NOT_PRESENT
   DIAG_FILE_COPIES(PRIOR_INF_SD) = COPY_NOT_PRESENT
endif

if(do_posterior_inflate) then
   DIAG_FILE_COPIES(POST_INF) = POST_INF_COPY
   DIAG_FILE_COPIES(POST_INF_SD) = POST_INF_SD_COPY
else
   DIAG_FILE_COPIES(POST_INF) = COPY_NOT_PRESENT
   DIAG_FILE_COPIES(POST_INF_SD) = COPY_NOT_PRESENT
endif

end subroutine count_state_ens_copies

!------------------------------------------------------------------------------

!> Set file name information.  For members restarts can be read from
!> an input_state_file_list or constructed using a stage name and
!> num_ens.  The file_info handle knows whether or not there is an
!> associated input_state_file_list. If no list is provided member
!> filenames are written as :
!>    stage_member_####.nc (ex. preassim_member_0001.nc)
!> extra copies are stored as :
!>    stage_basename.nc (ex. preassim_mean.nc)

subroutine set_filename_info(file_info, stage, num_ens, STAGE_COPIES)

type(file_info_type), intent(inout) :: file_info
character(len=*),     intent(in)    :: stage
integer,              intent(in)    :: num_ens
integer,              intent(inout) :: STAGE_COPIES(NUM_SCOPIES)

call set_member_file_metadata(file_info, num_ens, STAGE_COPIES(ENS_START))


STAGE_COPIES(ENS_END) = STAGE_COPIES(ENS_START) + num_ens - 1

call set_file_metadata(file_info, STAGE_COPIES(ENS_MEAN),      stage, 'mean',          'ensemble mean')
call set_file_metadata(file_info, STAGE_COPIES(ENS_SD),        stage, 'sd',            'ensemble sd')
call set_file_metadata(file_info, STAGE_COPIES(PRIOR_INF),     stage, 'priorinf_mean', 'prior inflation mean')
call set_file_metadata(file_info, STAGE_COPIES(PRIOR_INF_SD),   stage, 'priorinf_sd',   'prior inflation sd')
call set_file_metadata(file_info, STAGE_COPIES(POST_INF),  stage, 'postinf_mean',  'posterior inflation mean')
call set_file_metadata(file_info, STAGE_COPIES(POST_INF_SD),    stage, 'postinf_sd',    'posterior inflation sd')

end subroutine set_filename_info

!------------------------------------------------------------------------------

subroutine init_input_file_info(ncopies, ens_size, single_file_in, &
   perturb_from_single_instance, has_cycling, input_state_files,   &
   input_state_file_list, do_prior_inflate, prior_inflate_from_restart, &
   do_posterior_inflate, posterior_inflate_from_restart, file_info_input)
                                 
integer,              intent(in)  :: ncopies
integer,              intent(in)  :: ens_size
type(file_info_type), intent(out) :: file_info_input
logical,              intent(in)  :: single_file_in
logical,              intent(in)  :: perturb_from_single_instance
logical,              intent(in)  :: has_cycling
character(len=256),   intent(inout)  :: input_state_files(:)
logical,              intent(in)  :: do_prior_inflate, prior_inflate_from_restart
logical,              intent(in)  :: do_posterior_inflate
logical,              intent(in)  :: posterior_inflate_from_restart
character(len=256),   intent(in)  :: input_state_file_list(:)
                                 
character(len=256), allocatable :: file_array_input(:,:)
integer :: ninput_files

! This should have its own variant of DIAG_FILE_COPIES

! Determine number of files
if (single_file_in .or. perturb_from_single_instance)  then
   ninput_files = 1               
else
   ninput_files    = ens_size ! number of incomming ensemble members
endif

! Allocate space for file arrays.  contains a matrix of files (num_ens x num_domains)
! If perturbing from a single instance the number of input files does not have to
! be ens_size but rather a single file (or multiple files if more than one domain)
allocate(file_array_input(ninput_files, get_num_domains()))    
file_array_input  = RESHAPE(input_state_files,  (/ninput_files,  get_num_domains()/))
                       
! Given vector of input_state_files or a text file containing  
! a list of files, return a vector of files containing the filenames.
call set_multiple_filename_lists(input_state_files(:), input_state_file_list(:), &
   get_num_domains(), ninput_files, 'filter', 'input_state_files', 'input_state_file_list')

! Allocate space for the filename handles
call io_filenames_init(file_info_input, ncopies, has_cycling, single_file_in, &
   file_array_input, 'input')

! Set filename metadata information
call set_filename_info(file_info_input, 'input', ens_size, DIAG_FILE_COPIES )
   
! Set file IO information
call set_input_file_info( file_info_input, ens_size, perturb_from_single_instance, &
   do_prior_inflate, prior_inflate_from_restart, do_posterior_inflate,             &
   posterior_inflate_from_restart, DIAG_FILE_COPIES ) 

end subroutine init_input_file_info

!------------------------------------------------------------------------------

subroutine set_input_file_info( file_info, num_ens, perturb_from_single_instance, &
   do_prior_inflate, prior_inflate_from_restart, &
   do_posterior_inflate, posterior_inflate_from_restart, STAGE_COPIES )
                                  
type(file_info_type), intent(inout) :: file_info
integer,              intent(in)    :: num_ens
logical,              intent(in)    :: perturb_from_single_instance
logical,              intent(in)    :: do_prior_inflate, prior_inflate_from_restart
logical,              intent(in)    :: do_posterior_inflate
logical,              intent(in)    ::  posterior_inflate_from_restart
integer,              intent(in)    :: STAGE_COPIES(NUM_SCOPIES)                              
                                  
if ( perturb_from_single_instance ) then
   call set_io_copy_flag(file_info, STAGE_COPIES(ENS_START), READ_COPY)                       
   !>@todo know whether we are perturbing or not
   !#! call set_perturb_members(file_info, ENS_START, num_ens)
else                              
   call set_io_copy_flag(file_info, STAGE_COPIES(ENS_START), STAGE_COPIES(ENS_START)+num_ens-1, READ_COPY)
endif
                                  
if ( do_prior_inflate .and. prior_inflate_from_restart) then
   call set_io_copy_flag(file_info, STAGE_COPIES(PRIOR_INF), READ_COPY, inherit_units=.false.)
   call set_io_copy_flag(file_info, STAGE_COPIES(PRIOR_INF_SD), READ_COPY, inherit_units=.false.)
endif                             

if ( do_posterior_inflate .and. posterior_inflate_from_restart) then
   call set_io_copy_flag(file_info, STAGE_COPIES(POST_INF), READ_COPY, inherit_units=.false.)
   call set_io_copy_flag(file_info, STAGE_COPIES(POST_INF_SD), READ_COPY, inherit_units=.false.)
endif

end subroutine set_input_file_info

!------------------------------------------------------------------------------
   
subroutine set_output_file_info( file_info, num_ens, STAGE_COPIES, &
   output_mean, output_sd, do_prior_inflate, do_posterior_inflate, &
   output_members, do_clamping, force_copy)

type(file_info_type), intent(inout) :: file_info
integer,              intent(in)    :: num_ens
integer,              intent(in)    :: STAGE_COPIES(NUM_SCOPIES)
logical,              intent(in)    :: output_mean, output_sd
logical,              intent(in)    :: do_prior_inflate, do_posterior_inflate
logical,              intent(in)    :: output_members
logical,              intent(in)    :: do_clamping
logical,              intent(in)    :: force_copy
                           
!>@todo revisit if we should be clamping mean copy for file_info_output
if ( num_ens > 0 .and. output_members ) then
   call set_io_copy_flag(file_info, STAGE_COPIES(ENS_START), &
      STAGE_COPIES(ENS_START)+num_ens-1, WRITE_COPY, num_ens, &
      do_clamping, force_copy)
endif                      
                           
if(output_mean) call set_io_copy_flag(file_info, STAGE_COPIES(ENS_MEAN), WRITE_COPY, &
                   inherit_units=.true., clamp_vars=do_clamping, force_copy_back=force_copy)
if(output_sd) call set_io_copy_flag(file_info, STAGE_COPIES(ENS_SD), WRITE_COPY, &
                   inherit_units=.true., force_copy_back=force_copy)
if(do_prior_inflate) call set_io_copy_flag(file_info, STAGE_COPIES(PRIOR_INF), WRITE_COPY, &
                        inherit_units=.false., force_copy_back=force_copy)
if(do_prior_inflate) call set_io_copy_flag(file_info, STAGE_COPIES(PRIOR_INF_SD), &
                        WRITE_COPY, inherit_units=.false., force_copy_back=force_copy)
if (do_posterior_inflate) call set_io_copy_flag(file_info, STAGE_COPIES(POST_INF), &
                             WRITE_COPY, inherit_units=.false., force_copy_back=force_copy)
if (do_posterior_inflate) call set_io_copy_flag(file_info, STAGE_COPIES(POST_INF_SD), &
                             WRITE_COPY, inherit_units=.false., force_copy_back=force_copy)

end subroutine set_output_file_info

!------------------------------------------------------------------------------
!> initialize file names and which copies should be read and or written

subroutine init_output_file_info(output_file_name, ncopies, ens_size, file_info, &
    output_state_files, output_state_file_list, single_file_out,   &
    has_cycling, do_prior_inflate, do_posterior_inflate)
                                  
character(len=*),     intent(in)  :: output_file_name                                  
integer,              intent(in)  :: ncopies
integer,              intent(in)  :: ens_size
type(file_info_type), intent(out) :: file_info
character(len=256),   intent(inout)  :: output_state_files(:)
character(len=256),   intent(in)  :: output_state_file_list(:)
logical,              intent(in)  :: single_file_out
logical,              intent(in)  :: has_cycling
logical,              intent(in)  :: do_prior_inflate, do_posterior_inflate
   
integer :: noutput_members, noutput_files, ndomains
character(len=256), allocatable :: file_array_output(:,:)
                
! local variable to shorten the name for function input
ndomains        = get_num_domains()
noutput_files   = ens_size ! number of incomming ensemble members
   
! Assign the correct number of output files.
if (single_file_out) noutput_files = 1
   
! Given vector of output_state_files or a text file containing
! a list of files, return a vector of files containing the filenames.
call set_multiple_filename_lists(output_state_files(:), output_state_file_list(:), &
   ndomains, noutput_files, 'filter', 'output_state_files', 'output_state_file_list')
                                  
! Allocate space for file arrays.  contains a matrix of files (num_ens x num_domains)
allocate(file_array_output(noutput_files, ndomains))
file_array_output = RESHAPE(output_state_files, (/noutput_files, ndomains/))

! Write restart from output_state_file_list if provided
call io_filenames_init(file_info, ncopies, has_cycling, single_file_out, &
   file_array_output, output_file_name, check_output_compatibility = .true.)

!   Output Files
call set_filename_info(file_info, output_file_name, ens_size, DIAG_FILE_COPIES )

!   Output Files
call set_output_file_info(file_info, ens_size, DIAG_FILE_COPIES, &
   .true., .true., do_prior_inflate, do_posterior_inflate, &
   output_members = .true., do_clamping = .true., force_copy = .false. )

end subroutine init_output_file_info

!------------------------------------------------------------------------------
!> initialize file names and which copies should be read and or written

subroutine init_diag_file_info(diag_file_name, ncopies, ens_size, file_info, &
    single_file_out,   &
    has_cycling, output_mean, output_sd,             &
    output_members, do_prior_inflate, do_posterior_inflate)

character(len=*),     intent(in)  :: diag_file_name                                  
integer,              intent(in)  :: ncopies
integer,              intent(in)  :: ens_size
type(file_info_type), intent(out) :: file_info
logical,              intent(in)  :: single_file_out
logical,              intent(in)  :: has_cycling
logical,              intent(in)  :: output_mean, output_sd, output_members
logical,              intent(in)  :: do_prior_inflate, do_posterior_inflate
   
integer :: noutput_members, noutput_files, ndomains
character(len=256), allocatable :: file_array_output(:,:)
                
! local variable to shorten the name for function input
noutput_members = num_output_state_members 
ndomains        = get_num_domains()
noutput_files   = ens_size ! number of incomming ensemble members
   
! Assign the correct number of output files.
if (single_file_out) noutput_files = 1
   
! Output Files (we construct the filenames)
call io_filenames_init(file_info, ncopies, has_cycling, single_file_out, &
   root_name = diag_file_name)

!   Output Files
call set_filename_info(file_info, diag_file_name,  noutput_members,  DIAG_FILE_COPIES )

!   Output Files
call set_output_file_info(file_info, noutput_members, DIAG_FILE_COPIES,  &
   output_mean, output_sd, do_prior_inflate, do_posterior_inflate, output_members, &
   do_clamping  = .false., force_copy = .true.)

end subroutine init_diag_file_info

!------------------------------------------------------------------------------

end module filter_io_diag_mod
