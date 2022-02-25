! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Generate initial inflation files from namelist values.
!> This way an experiment can always start from a restart file 
!> without having to alter the namelist between cycles 1 and 2.
!>
!> an alternative to running this program is to use the nco utilities thus:
!>
!> Here is an example using version 4.4.2 or later of the NCO tools:
!>   ncap2 -s "T=1.0;U=1.0;V=1.0" wrfinput_d01 prior_inflation_mean.nc
!>   ncap2 -s "T=0.6;U=0.6;V=0.6" wrfinput_d01 prior_inflation_sd.nc'

program fill_inflation_restart

use             types_mod, only : r8, i8, missing_r8

use         utilities_mod, only : error_handler, E_MSG, E_ERR, &
                                  find_namelist_in_file, check_namelist_read,   &
                                  nmlfileunit, do_nml_file, do_nml_term

use       assim_model_mod, only : static_init_assim_model

use      time_manager_mod, only : time_type

use  ensemble_manager_mod, only : init_ensemble_manager, ensemble_type, &
                                  end_ensemble_manager

use   state_vector_io_mod, only : state_vector_io_init, read_state, write_state

use   state_structure_mod, only : get_num_domains, set_update_list, get_num_variables

use      io_filenames_mod, only : io_filenames_init, file_info_type,       &
                                  stage_metadata_type, get_stage_metadata, &
                                  get_restart_filename, set_file_metadata, &
                                  set_io_copy_flag, READ_COPY, WRITE_COPY

use             model_mod, only : get_model_size

use     mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities

implicit none

character(len=*), parameter :: source = 'fill_inflation_restart.f90'

! MAX_FILES is max number of domains
integer, parameter :: MAX_FILES = 10
integer, parameter :: ss_inflate_index    = 1 
integer, parameter :: ss_inflate_sd_index = 2

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

logical            :: write_prior_inf = .FALSE. 
real(r8)           :: prior_inf_mean = MISSING_R8
real(r8)           :: prior_inf_sd   = MISSING_R8
logical            :: write_post_inf = .FALSE.
real(r8)           :: post_inf_mean  = MISSING_R8
real(r8)           :: post_inf_sd    = MISSING_R8
logical            :: single_file    = .FALSE.
character(len=256) :: input_state_files(MAX_FILES)  = ''
logical            :: verbose        = .FALSE.

namelist /fill_inflation_restart_nml/              &
             prior_inf_mean, prior_inf_sd,         &
             post_inf_mean, post_inf_sd, verbose,  &
             write_prior_inf, write_post_inf,      &
             input_state_files, single_file

! io variables
integer                   :: iunit, io
type(file_info_type)      :: file_info_input
type(file_info_type)      :: file_info_output
type(stage_metadata_type) :: input_restart_files
type(stage_metadata_type) :: output_restart_files


! model state variables
type(ensemble_type)   :: ens_handle
type(time_type)       :: model_time
integer(i8)           :: model_size
integer               :: ncopies = 2    ! {prior,posterior}_inf_{mean,sd}
logical               :: read_time_from_file = .true.

! counter variables
integer :: idom, imem, ndomains, ivars

! message strings
character(len=512) :: my_base, my_desc, my_stage
character(len=512) :: string1, string2, string3

character(len=256), allocatable :: file_array_input(:,:)

!======================================================================
! start of executable code
!======================================================================

call initialize_modules_used()

call find_namelist_in_file("input.nml", "fill_inflation_restart_nml", iunit)
read(iunit, nml = fill_inflation_restart_nml, iostat = io)
call check_namelist_read(iunit, io, "fill_inflation_restart_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=fill_inflation_restart_nml)
if (do_nml_term()) write(     *     , nml=fill_inflation_restart_nml)

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
! read/write restart files
!----------------------------------------------------------------------

! Set up the ensemble storage and read in the restart file
call init_ensemble_manager(ens_handle, ncopies, model_size)

! Allocate space for file arrays.  
! Contains a matrix of files (ncopies x ndomains)
! If perturbing from a single instance the number of 
! input files does not have to be ens_size but rather 
! a single file (or multiple files if more than one domain)

ndomains = get_num_domains()

allocate(file_array_input(ncopies, ndomains)) 

file_array_input  = RESHAPE(input_state_files,  (/1,  ndomains/))

if(single_file) then
  call error_handler(E_ERR, 'fill_inflation_restart: ', &
                     'single_file not yet supported, please contact DART', source)
endif

call io_filenames_init(file_info_input,             &
                       ncopies      = ncopies,      &
                       cycling      = single_file,  &
                       single_file  = single_file,  &
                       restart_files = file_array_input(:,:))

! Read the template file to get the shape of netCDF file
! and its variables. It is possible to have multiple domains
! but only require one member.
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
   write(string1, *) '- Reading File : ', &
                     trim(get_restart_filename(input_restart_files, &
                                               copy   = imem,       &
                                               domain = idom))
   call error_handler(E_MSG, 'fill_inflation_restart: ', string1, source)
enddo

call read_state(ens_handle, file_info_input, read_time_from_file, model_time)

! if users have variables in the dart state which are set to 'no update'
! for assimilation purposes, they still have to have corresponding fields 
! in the inflation files. (science question is what those values should be.)
!
! the normal write routine skips 'no update' fields, so the output of
! this program would be missing those fields and cause an error when
! filter tries to read the inflation files.  the following code 
! overwrites the update flag values so all variables are always written.
!
! also tell users to script their workflows to copy the input inflation
! file to the output name before running filter, so there is an existing
! output inflation file for filter to update. it will already include the
! 'no update' fields, whose inflation values will remain unchanged. 
! this prevents confusing errors at the next execution of filter.
! [ed note: tim is brilliant.]

do idom = 1, ndomains
   ivars = get_num_variables(idom)   
   call set_update_list(idom, ivars, (/ ( .TRUE., imem=1,ivars ) /) )
enddo

! Initialize file output to default inflation file names
call io_filenames_init(file_info_output,           &
                       ncopies      = ncopies,     &
                       cycling      = single_file, &
                       single_file  = single_file)

if(write_prior_inf) then  
   call fill_inflation_files(prior_inf_mean, prior_inf_sd, 'prior')
endif

if(write_post_inf) then  
   call fill_inflation_files(post_inf_mean,  post_inf_sd,   'post')
endif

deallocate(file_array_input)

call finalize_modules_used()

!======================================================================
contains
!======================================================================

!> fill inflation values in a separate file

subroutine fill_inflation_files(inf_mean, inf_sd, stage)
real(r8),         intent(in) :: inf_mean
real(r8),         intent(in) :: inf_sd
character(len=*), intent(in) :: stage

if (inf_mean == MISSING_R8 .or. inf_sd == MISSING_R8) then
   write(string1,*) 'you must specify both inflation mean and inflation standard deviation values'
   write(string2,*) 'you have "',trim(stage),'_inf_mean = ', inf_mean,'" and '
   write(string3,*) '         "',trim(stage),'_inf_sd   = ', inf_sd,  '"     '
   call error_handler(E_MSG, 'fill_inflation_restart: ', string1,  &
                      source, text2=string2, text3=string3)
   return
endif
ens_handle%copies(ss_inflate_index   , :) = inf_mean
ens_handle%copies(ss_inflate_sd_index, :) = inf_sd

write(my_stage,'(3A)') 'input_', stage, 'inf'
write(my_base, '(A)')  'mean'
write(my_desc, '(2A)') stage, ' inflation mean'
call set_file_metadata(file_info_output,    &
                       cnum     = ss_inflate_index, &
                       stage    = my_stage, &
                       basename = my_base,  &
                       desc     = my_desc)

call set_io_copy_flag(file_info_output,    &
                      cnum    = 1,         &
                      io_flag = WRITE_COPY)

write(my_base, '(A)') 'sd'
write(my_desc, '(2A)') stage, ' inflation sd'
call set_file_metadata(file_info_output,    &
                       cnum     = ss_inflate_sd_index, &
                       stage    = my_stage, &
                       basename = my_base,  &
                       desc     = my_desc)

call set_io_copy_flag(file_info_output,    &
                      cnum    = 2,         &
                      io_flag = WRITE_COPY)

output_restart_files = get_stage_metadata(file_info_output)

do idom = 1, ndomains
   ! write inflation mean
   write(string1, *) '- Writing ',trim(stage), ' Mean Inflation File : "',   &
                     trim(get_restart_filename(output_restart_files,         &
                                               copy   = ss_inflate_index,    &
                                               domain = idom)), '"'
   write(string2, *) ' with value = ', inf_mean
   call error_handler(E_MSG, 'fill_inflation_restart: ', string1, source, text2=string2)

   ! write inflation stadard deviation
   write(string1, *) '- Writing ',trim(stage), ' SD   Inflation File : "',   &
                     trim(get_restart_filename(output_restart_files,         &
                                               copy   = ss_inflate_sd_index, &
                                               domain = idom)), '"'
   write(string2, *) ' with value = ', inf_sd
   call error_handler(E_MSG, 'fill_inflation_restart: ', string1, source, text2=string2)
enddo


call write_state(ens_handle, file_info_output)

end subroutine fill_inflation_files

!----------------------------------------------------------------------
!> initialize modules that need it

subroutine initialize_modules_used()

call initialize_mpi_utilities('fill_inflation_restart')


call state_vector_io_init()

end subroutine initialize_modules_used

!----------------------------------------------------------------------
!> clean up before exiting

subroutine finalize_modules_used()

call end_ensemble_manager(ens_handle)

! this must be last, and you can't print/write anything
! after this is called.
call finalize_mpi_utilities()

end subroutine finalize_modules_used

end program fill_inflation_restart

