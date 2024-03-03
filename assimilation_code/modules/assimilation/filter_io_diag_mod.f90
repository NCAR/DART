! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
      
module filter_io_diag_mod
     
!------------------------------------------------------------------------------

use types_mod,             only : r8, missing_r8

use ensemble_manager_mod,  only : ensemble_type, compute_copy_mean_sd

use time_manager_mod,      only : time_type

use random_seq_mod,        only : random_seq_type, init_random_seq, random_gaussian

use io_filenames_mod,      only : file_info_type, set_member_file_metadata, &
                                  set_file_metadata, COPY_NOT_PRESENT, &
                                  io_filenames_init, set_io_copy_flag, &
                                  check_file_info_variable_shape, READ_COPY, WRITE_COPY

use assim_model_mod,       only : pert_model_copies

use state_vector_io_mod,   only : write_state, read_state

use adaptive_inflate_mod,  only : do_rtps_inflate, adaptive_inflate_type

use utilities_mod,         only : set_multiple_filename_lists

use state_structure_mod,   only : get_num_domains

!------------------------------------------------------------------------------

implicit none
private

public ::  create_ensemble_from_single_file, ens_copies_type, &
   count_state_ens_copies, read_state_and_inflation, output_diagnostics, init_output_file_info

!------------------------------------------------------------------------------

type ens_copies_type
   integer           :: ens_size 
   integer           :: num_state_ens_copies
   integer           :: num_extras
   integer           :: num_output_state_members

   ! A list of all possible copies supported here 
   integer           :: ENS_START_COPY         = COPY_NOT_PRESENT
   integer           :: ENS_END_COPY           = COPY_NOT_PRESENT
   integer           :: ENS_MEAN_COPY          = COPY_NOT_PRESENT
   integer           :: ENS_SD_COPY            = COPY_NOT_PRESENT
   integer           :: PRIOR_INF_COPY         = COPY_NOT_PRESENT
   integer           :: PRIOR_INF_SD_COPY      = COPY_NOT_PRESENT
   integer           :: POST_INF_COPY          = COPY_NOT_PRESENT
   integer           :: POST_INF_SD_COPY       = COPY_NOT_PRESENT
   integer           :: RTPS_PRIOR_SPREAD_COPY = COPY_NOT_PRESENT

   ! Pointers to where different copies are stored in array
   integer            :: ENS_START    = 1
   integer            :: ENS_END      = 2
   integer            :: ENS_MEAN     = 3
   integer            :: ENS_SD       = 4
   integer            :: PRIOR_INF    = 5
   integer            :: PRIOR_INF_SD = 6
   integer            :: POST_INF     = 7 
   integer            :: POST_INF_SD  = 8 

   ! Data structure for different diag and output stages
   integer            :: DIAG_FILE_COPIES(1:8)     = COPY_NOT_PRESENT
end type ens_copies_type
                                  
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

! Assign the indices to the storage in the ensemble manager
! Initialize indices for where diagnostic copies are found

subroutine count_state_ens_copies(ens_copies, ens_size, output_mean, output_sd, &
   do_prior_inflate, do_posterior_inflate, posterior_inflate, num_output_state_members)
   

type(ens_copies_type),       intent(inout) :: ens_copies
integer,                     intent(in)  :: ens_size
logical,                     intent(in)  :: output_mean, output_sd
logical,                     intent(in)  :: do_prior_inflate, do_posterior_inflate
type(adaptive_inflate_type), intent(in)  :: posterior_inflate
integer,                     intent(in) :: num_output_state_members

! Initialize the fixed size information
ens_copies%ens_size = ens_size

ens_copies%num_output_state_members = min(num_output_state_members, ens_size)

! First ens_size entries are the actual ensemble members
ens_copies%ENS_START_COPY = 1
ens_copies%ENS_END_COPY    = ens_size

! Filter Extra Copies For Assimilation
ens_copies%ENS_MEAN_COPY        = ens_size + 1
ens_copies%ENS_SD_COPY          = ens_size + 2
ens_copies%PRIOR_INF_COPY       = ens_size + 3
ens_copies%PRIOR_INF_SD_COPY    = ens_size + 4
ens_copies%POST_INF_COPY        = ens_size + 5
ens_copies%POST_INF_SD_COPY     = ens_size + 6
ens_copies%num_state_ens_copies = ens_size + 6

! If Whitaker/Hamill (2012) relaxation-to-prior-spread (RTPS) inflation
! then we need an extra copy to hold (save) the prior ensemble spread
if ( do_rtps_inflate(posterior_inflate) ) then
   ens_copies%RTPS_PRIOR_SPREAD_COPY = ens_size + 7
   ens_copies%num_state_ens_copies = ens_size + 7
endif                  
   
! Copies in excess of ensemble
ens_copies%num_extras = ens_copies%num_state_ens_copies - ens_size

! Specifying copies to be output for state diagnostics files
ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_START) = ens_copies%ENS_START_COPY
ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_END)   = ens_copies%ENS_END_COPY
if(output_mean) then
   ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_MEAN) = ens_copies%ENS_MEAN_COPY
else
   ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_MEAN) = COPY_NOT_PRESENT
endif
if(output_sd) then
   ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_SD) = ens_copies%ENS_SD_COPY
else
   ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_SD) = COPY_NOT_PRESENT
endif
if(do_prior_inflate) then
   ens_copies%DIAG_FILE_COPIES(ens_copies%PRIOR_INF)    = ens_copies%PRIOR_INF_COPY
   ens_copies%DIAG_FILE_COPIES(ens_copies%PRIOR_INF_SD) = ens_copies%PRIOR_INF_SD_COPY
else
   ens_copies%DIAG_FILE_COPIES(ens_copies%PRIOR_INF)    = COPY_NOT_PRESENT
   ens_copies%DIAG_FILE_COPIES(ens_copies%PRIOR_INF_SD) = COPY_NOT_PRESENT
endif

if(do_posterior_inflate) then
   ens_copies%DIAG_FILE_COPIES(ens_copies%POST_INF)    = ens_copies%POST_INF_COPY
   ens_copies%DIAG_FILE_COPIES(ens_copies%POST_INF_SD) = ens_copies%POST_INF_SD_COPY
else
   ens_copies%DIAG_FILE_COPIES(ens_copies%POST_INF)    = COPY_NOT_PRESENT
   ens_copies%DIAG_FILE_COPIES(ens_copies%POST_INF_SD) = COPY_NOT_PRESENT
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

subroutine set_filename_info(file_info, stage, ens_copies, num_ens)

type(file_info_type), intent(inout) :: file_info
character(len=*),     intent(in)    :: stage
type(ens_copies_type), intent(inout) :: ens_copies
integer,              intent(in)    :: num_ens

call set_member_file_metadata(file_info, num_ens, ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_START))

ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_END) = &
   ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_START) + num_ens - 1

call set_file_metadata(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_MEAN), &
   stage, 'mean', 'ensemble mean')
call set_file_metadata(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_SD), &
   stage, 'sd', 'ensemble sd')
call set_file_metadata(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%PRIOR_INF), &
   stage, 'priorinf_mean', 'prior inflation mean')
call set_file_metadata(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%PRIOR_INF_SD), &
   stage, 'priorinf_sd', 'prior inflation sd')
call set_file_metadata(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%POST_INF), &
   stage, 'postinf_mean',  'posterior inflation mean')
call set_file_metadata(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%POST_INF_SD), &
   stage, 'postinf_sd',    'posterior inflation sd')

end subroutine set_filename_info

!------------------------------------------------------------------------------

subroutine read_state_and_inflation(ens_copies, state_ens_handle, single_file_in, &
   perturb_from_single_instance, has_cycling, input_state_files,   &
   input_state_file_list, prior_inflate, do_prior_inflate, prior_inflate_from_restart, &
   posterior_inflate, do_posterior_inflate, posterior_inflate_from_restart, file_info_input, &
   read_time_from_file, time1)
                                
type(ens_copies_type),intent(inout)  :: ens_copies 
type(ensemble_type),  intent(inout) :: state_ens_handle
type(file_info_type), intent(out) :: file_info_input
logical,              intent(in)  :: single_file_in
logical,              intent(in)  :: perturb_from_single_instance
logical,              intent(in)  :: has_cycling
character(len=256),   intent(inout)  :: input_state_files(:)
type(adaptive_inflate_type), intent(in) :: prior_inflate
logical,              intent(in)  :: do_prior_inflate, prior_inflate_from_restart
type(adaptive_inflate_type), intent(in) :: posterior_inflate
logical,              intent(in)  :: do_posterior_inflate
logical,              intent(in)  :: posterior_inflate_from_restart
character(len=256),   intent(in)  :: input_state_file_list(:)
logical,              intent(in)  :: read_time_from_file
type(time_type),      intent(inout) :: time1

                                 
character(len=256), allocatable :: file_array_input(:,:)
integer :: ninput_files

! Best to fail if already initialized
!!!if(.not. file_info_input%initialized) then

! Determine number of files
if (single_file_in .or. perturb_from_single_instance)  then
   ninput_files = 1               
else
   ninput_files    = ens_copies%ens_size ! number of incomming ensemble members
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
call io_filenames_init(file_info_input, ens_copies%num_state_ens_copies, has_cycling, single_file_in, &
   file_array_input, 'input')

! Set filename metadata information
call set_filename_info(file_info_input, 'input', ens_copies, ens_copies%ens_size)
   
! Set file IO information
call set_input_file_info( file_info_input, ens_copies, perturb_from_single_instance, &
   do_prior_inflate, prior_inflate_from_restart, do_posterior_inflate,             &
   posterior_inflate_from_restart)

! Do the reading
call read_state(state_ens_handle, file_info_input, read_time_from_file, time1,      &
                ens_copies%PRIOR_INF_COPY, ens_copies%PRIOR_INF_SD_COPY, &
                ens_copies%POST_INF_COPY, ens_copies%POST_INF_SD_COPY, &
                prior_inflate, posterior_inflate,                                   &
                prior_inflate_from_restart, posterior_inflate_from_restart,         &         
                perturb_from_single_instance)

end subroutine read_state_and_inflation

!------------------------------------------------------------------------------

subroutine set_input_file_info( file_info, ens_copies, perturb_from_single_instance, &
   do_prior_inflate, prior_inflate_from_restart, &
   do_posterior_inflate, posterior_inflate_from_restart)
                                  
type(file_info_type), intent(inout) :: file_info
type(ens_copies_type), intent(inout) :: ens_copies
logical,              intent(in)    :: perturb_from_single_instance
logical,              intent(in)    :: do_prior_inflate, prior_inflate_from_restart
logical,              intent(in)    :: do_posterior_inflate
logical,              intent(in)    ::  posterior_inflate_from_restart
     
if ( perturb_from_single_instance ) then
   ! Only reading a single ensemble member, the first one                             
   call set_io_copy_flag(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_START), READ_COPY)                       
else                              
   ! Will read all the ensemble members
   call set_io_copy_flag(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_START), &
      ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_START) + ens_copies%ens_size - 1, READ_COPY)
endif

! Reading prior inflation mean and sd from restart                                  
if ( do_prior_inflate .and. prior_inflate_from_restart) then
   call set_io_copy_flag(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%PRIOR_INF), &
       READ_COPY, inherit_units=.false.)
   call set_io_copy_flag(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%PRIOR_INF_SD), &
      READ_COPY, inherit_units=.false.)
endif                             

! Reading posterior inflation mean and sd from restart                                  
if ( do_posterior_inflate .and. posterior_inflate_from_restart) then
   call set_io_copy_flag(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%POST_INF), &
      READ_COPY, inherit_units=.false.)
   call set_io_copy_flag(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%POST_INF_SD), &
      READ_COPY, inherit_units=.false.)
endif

end subroutine set_input_file_info

!------------------------------------------------------------------------------
   
subroutine set_output_file_info( file_info, ens_copies, num_ens, &
   output_mean, output_sd, do_prior_inflate, do_posterior_inflate, &
   output_members, do_clamping, force_copy)

type(file_info_type), intent(inout) :: file_info
type(ens_copies_type), intent(inout) :: ens_copies
integer,              intent(in)    :: num_ens
logical,              intent(in)    :: output_mean, output_sd
logical,              intent(in)    :: do_prior_inflate, do_posterior_inflate
logical,              intent(in)    :: output_members
logical,              intent(in)    :: do_clamping
logical,              intent(in)    :: force_copy

! Output file ensemble members                            
if ( num_ens > 0 .and. output_members ) then
   call set_io_copy_flag(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_START), &
      ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_START)+num_ens-1, WRITE_COPY, num_ens, &
      do_clamping, force_copy)
endif                      

! Output file copies for the mean and sd if requested                           
if(output_mean) call set_io_copy_flag(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_MEAN), &
   WRITE_COPY, inherit_units=.true., clamp_vars=do_clamping, force_copy_back=force_copy)
if(output_sd) call set_io_copy_flag(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%ENS_SD), &
   WRITE_COPY, inherit_units=.true., force_copy_back=force_copy)

! Output file copies for the prior inflation and sd
if(do_prior_inflate) then
   call set_io_copy_flag(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%PRIOR_INF), &
      WRITE_COPY, inherit_units=.false., force_copy_back=force_copy)
   call set_io_copy_flag(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%PRIOR_INF_SD), &
      WRITE_COPY, inherit_units=.false., force_copy_back=force_copy)
endif

! Output copies for the posterior inflation and sd
if (do_posterior_inflate) then
   call set_io_copy_flag(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%POST_INF), &
      WRITE_COPY, inherit_units=.false., force_copy_back=force_copy)
   call set_io_copy_flag(file_info, ens_copies%DIAG_FILE_COPIES(ens_copies%POST_INF_SD), &
      WRITE_COPY, inherit_units=.false., force_copy_back=force_copy)
endif

end subroutine set_output_file_info

!------------------------------------------------------------------------------
!> initialize file names and which copies should be read and or written

subroutine init_output_file_info(output_file_name, state_ens_handle, ens_copies, &
   file_info, output_state_files, output_state_file_list, single_file_out,   &
   has_cycling, do_prior_inflate, do_posterior_inflate)
                                  
character(len=*),     intent(in)  :: output_file_name                                  
type(ensemble_type),  intent(in)  :: state_ens_handle
type(ens_copies_type),intent(inout)  :: ens_copies
type(file_info_type), intent(out) :: file_info
character(len=256),   intent(inout)  :: output_state_files(:)
character(len=256),   intent(in)  :: output_state_file_list(:)
logical,              intent(in)  :: single_file_out
logical,              intent(in)  :: has_cycling
logical,              intent(in)  :: do_prior_inflate, do_posterior_inflate
   
integer :: noutput_files, ndomains
character(len=256), allocatable :: file_array_output(:,:)
                
! local variable to shorten the name for function input
ndomains        = get_num_domains()
noutput_files   = ens_copies%ens_size ! number of incomming ensemble members
   
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
call io_filenames_init(file_info, ens_copies%num_state_ens_copies, has_cycling, single_file_out, &
   file_array_output, output_file_name, check_output_compatibility = .true.)

!   Output Files
call set_filename_info(file_info, output_file_name, ens_copies, ens_copies%ens_size)

!   Output Files
call set_output_file_info(file_info, ens_copies, ens_copies%ens_size, &
   .true., .true., do_prior_inflate, do_posterior_inflate, &
   output_members = .true., do_clamping = .true., force_copy = .false. )

! Make sure the size of the variables in the file and the ensemble storage are consistent
call check_file_info_variable_shape(file_info, state_ens_handle)

end subroutine init_output_file_info

!------------------------------------------------------------------------------
!> initialize diagnostic file if needed and output current diagnostics from state ensemble

subroutine output_diagnostics(diag_file_name, state_ens_handle, ens_copies, &
   file_info, single_file_out, has_cycling, output_mean, output_sd, output_members, &
   do_prior_inflate, do_posterior_inflate)

character(len=*),     intent(in)  :: diag_file_name                                  
type(ensemble_type),  intent(inout) :: state_ens_handle
type(ens_copies_type), intent(inout) :: ens_copies
type(file_info_type), intent(inout) :: file_info
logical,              intent(in)  :: single_file_out
logical,              intent(in)  :: has_cycling
logical,              intent(in)  :: output_mean, output_sd, output_members
logical,              intent(in)  :: do_prior_inflate, do_posterior_inflate
   
integer :: noutput_files, ndomains
           
! Don't need to initialize if already done
if(.not. file_info%initialized) then
     
   ! local variable to shorten the name for function input
   ndomains        = get_num_domains()
   noutput_files   = ens_copies%ens_size ! number of incomming ensemble members
   
   ! Assign the correct number of output files.
   if (single_file_out) noutput_files = 1
   
   ! Output Files (we construct the filenames)
   call io_filenames_init(file_info, ens_copies%num_state_ens_copies, has_cycling, single_file_out, &
      root_name = diag_file_name)

   !   Output Files
   call set_filename_info(file_info, diag_file_name,  ens_copies, ens_copies%num_output_state_members)

   !   Output Files
   call set_output_file_info(file_info, ens_copies, ens_copies%num_output_state_members, &
      output_mean, output_sd, do_prior_inflate, do_posterior_inflate, output_members, &
      do_clamping  = .false., force_copy = .true.)

endif

! Compute the ensemble mean and standard deviation; for efficiency could make this optional
call compute_copy_mean_sd(state_ens_handle, 1, ens_copies%ens_size, &
   ens_copies%ENS_MEAN_COPY, ens_copies%ENS_SD_COPY)

! Write state diagnostics to the netcdf file
call write_state(state_ens_handle, file_info)

end subroutine output_diagnostics

!------------------------------------------------------------------------------

end module filter_io_diag_mod
