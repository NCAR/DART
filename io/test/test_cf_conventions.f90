! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: test_cf_conventions.f90 9553 2016-01-20 17:26:41Z hendric $

program test_cf_conventions

use             types_mod, only : r4, r8, i8, metadatalength, MISSING_R8 
use         utilities_mod, only : register_module, error_handler, E_MSG
use  adaptive_inflate_mod, only : adaptive_inflate_end, adaptive_inflate_init, &
                                  adaptive_inflate_type
use     mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities
use       assim_model_mod, only : static_init_assim_model, get_model_size
use   state_vector_io_mod, only : read_state, write_state
use  ensemble_manager_mod, only : init_ensemble_manager, ensemble_type, set_num_extra_copies
use      io_filenames_mod, only : io_filenames_init, end_io_filenames, file_info_type, &
                                  get_output_file
use     copies_on_off_mod, only : ENS_MEAN_COPY, ENS_SD_COPY, &
                                  PRIOR_INF_COPY, PRIOR_INF_SD_COPY, &
                                  POST_INF_COPY, POST_INF_SD_COPY, &
                                  SPARE_PRIOR_MEAN, SPARE_PRIOR_SPREAD, &
                                  SPARE_PRIOR_INF_MEAN, SPARE_PRIOR_INF_SPREAD, &
                                  SPARE_POST_INF_MEAN, SPARE_POST_INF_SPREAD
use state_space_diag_mod,  only : filter_state_space_diagnostics, netcdf_file_type, &
                                  init_diag_output, finalize_diag_output,           &
                                  skip_diag_files
use   state_structure_mod, only : get_xtype,             &
                                  get_units,             &
                                  get_long_name,         &
                                  get_short_name,        &
                                  get_has_missing_value, &
                                  get_FillValue,         &
                                  get_missing_value,     &
                                  get_add_offset,        &
                                  get_scale_factor
use      time_manager_mod, only : time_type, set_time
use            filter_mod, only : filter_set_initial_time
use            assert_mod, only : assert_equal

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_test_read_write_restarts_dir/io/test/test_cf_conventions.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 9553 $"
character(len=128), parameter :: revdate  = "$Date: 2016-01-20 10:26:41 -0700 (Wed, 20 Jan 2016) $"

logical, save :: module_initialized = .false.

type(ensemble_type)         :: ens_handle
type(file_info_type)        :: file_info
type(time_type)             :: time1
type(adaptive_inflate_type) :: prior_inflate_handle, post_inflate_handle
type(netcdf_file_type)      :: PriorStateUnit_handle, PosteriorStateUnit_handle
type(time_type)             :: curr_ens_time

integer(i8) :: model_size
integer     :: ens_size, num_extras, num_copies

logical :: read_time_from_file

integer :: num_output_state_members = 3
integer :: output_state_mean_index, output_state_spread_index

logical :: output_inflation  = .true. ! This is for the diagnostic files, no separate option for prior and posterior

integer :: domid = 1 ! only one domain
integer :: var_xtype
character(len=NF90_MAX_NAME) :: var_units, blank_string, var_att_name

integer  :: missINT
real(r4) :: missR4
real(r8) :: missR8, var_offset, var_scale_factor

blank_string = ' '

! main code here
 
! initialize the dart libs
call initialize_module()

model_size = get_model_size()

write(*,*) " model size : ", model_size

ens_size   = 3
num_extras = 12

num_copies = ens_size + num_extras

! initalize routines needed for read_state and write_state
call init_ensemble_manager(ens_handle, num_copies, model_size)
call set_num_extra_copies(ens_handle, num_extras)
call filter_set_initial_time(0,0,time1,read_time_from_file)
call initialize_copy_numbers(ens_size)
call initialize_adaptive_inflate(ens_handle, prior_inflate_handle, post_inflate_handle)
call initialize_diagnostics(PriorStateUnit_handle, PosteriorStateUnit_handle)
file_info = initialize_filenames(ens_handle, overwrite_state_input=.false.)

curr_ens_time = set_time(0, 0)

! read in restarts
call read_state(ens_handle, file_info, read_time_from_file, time1)

! If needed, store copies(mean, sd, inf_mean, inf_sd) that would have
! gone in Prior_Diag.nc and write them at the end.
call store_prior(ens_handle)

call filter_state_space_diagnostics(file_info, curr_ens_time, PriorStateUnit_handle, ens_handle, &
      model_size, num_output_state_members, &
      output_state_mean_index, output_state_spread_index, output_inflation,&
      ENS_MEAN_COPY, ENS_SD_COPY, &
      prior_inflate_handle, PRIOR_INF_COPY, PRIOR_INF_SD_COPY)

! write out all files possible including restarts, mean, sd, and prior/posterior inflation files
call write_state(ens_handle, file_info, prior_inflate_handle, post_inflate_handle, skip_diag_files())

write(*,*)' '
write(*,*)'======================================================================'
write(*,*)' Unit Test for CF-Conventions'
write(*,*)'======================================================================'
write(*,*)' '


write(*,*)'Testing get_xtype'

var_xtype = get_xtype(domid,1)
call assert_equal(var_xtype, NF90_INT,    'variable1:get_xtype')

var_xtype = get_xtype(domid,2)
call assert_equal(var_xtype, NF90_FLOAT,  'variable2:get_xtype')

var_xtype = get_xtype(domid,3)
call assert_equal(var_xtype, NF90_DOUBLE, 'variable3:get_xtype')

var_xtype = get_xtype(domid,3)
call assert_equal(var_xtype, NF90_DOUBLE, 'variable4:get_xtype')

! test units

write(*,*)'Testing get_units'

var_units = get_units(domid,1)
call assert_equal(var_units, 'units A',     'variable1:get_units')

var_units = get_units(domid,2)
call assert_equal(var_units, 'units B',     'variable2:get_units')

var_units = get_units(domid,3)
call assert_equal(var_units, 'units C',     'variable3:get_units')

var_units = get_units(domid,4)
call assert_equal(var_units, blank_string,  'variable4:get_units')


write(*,*)'Testing get_long_name'

var_att_name = get_long_name(domid,1)
call assert_equal(var_att_name, 'variable A',  'variable1:get_long_name')

var_att_name = get_long_name(domid,2)
call assert_equal(var_att_name, 'variable B',  'variable2:get_long_name')

var_att_name = get_long_name(domid,3)
call assert_equal(var_att_name, 'variable C',  'variable3:get_long_name')

! var_att_name = get_long_name(domid,4)
! call assert_equal(var_att_name, 'D',  'variable4:get_long_name')


write(*,*)'Testing get_short_name'

var_att_name = get_short_name(domid,1)
call assert_equal(var_att_name, 'short A',     'variable1:get_short_name')

var_att_name = get_short_name(domid,2)
call assert_equal(var_att_name, 'short B',     'variable2:get_short_name')

var_att_name = get_short_name(domid,3)
call assert_equal(var_att_name, 'short C',     'variable3:get_short_name')

var_att_name = get_short_name(domid,4)

call assert_equal(var_att_name, blank_string,  'variable4:get_short_name')


write(*,*)'Testing get_missing_value'

call get_missing_value(domid,1,missINT)
call assert_equal(missINT, -77,                   'variable1:get_missing_value')

call get_missing_value(domid,2,missR4)
call assert_equal(missR4, -777.77_r4,             'variable2:get_missing_value')

call get_missing_value(domid,3,missR8)
call assert_equal(missR8, -88888.88888_r8,        'variable3:get_missing_value')

! ! variable 4 has no missing value this will fail
! call get_missing_value(domid,4,missR8)
! call assert_equal(missR8, -88888.88888_r8,        'variable4:get_missing_value')

write(*,*)'Testing get_FillValue'

call get_FillValue(domid,1,missINT)
call assert_equal(missINT, -77,                   'variable1:get_FillValue')

call get_FillValue(domid,2,missR4)
call assert_equal(missR4, -777.77_r4,             'variable2:get_FillValue')

call get_FillValue(domid,3,missR8)
call assert_equal(missR8, -88888.88888_r8,        'variable3:get_FillValue')

! ! variable 4 has no _FillValue this will fail
! call get_FillValue(domid,4,missR8)
! call assert_equal(missR8, -88888.88888_r8,        'variable4:get_FillValue')

write(*,*)'Testing offset and scale factor'

! this is only for r8 at the moment
! since it is not being used within DART

var_offset = get_add_offset(domid,3)
call assert_equal(var_offset, 2.0_r8,            'variable3:get_var_offset')

var_offset = get_add_offset(domid,4)
call assert_equal(var_offset, MISSING_R8 ,       'variable4:get_var_offset')


var_scale_factor = get_scale_factor(domid,3)
call assert_equal(var_scale_factor, 0.2_r8 ,     'variable3:get_scale_factor')

var_scale_factor = get_scale_factor(domid,4)
call assert_equal(var_scale_factor, MISSING_R8 , 'variable3:get_scale_factor')

write(*,*)' '
write(*,*)'======================================================================'
write(*,*)' Finished Unit Test'
write(*,*)'======================================================================'
write(*,*)' '

!>@todo FIXME : do we want to go through and do a unit test for the files that are being output?
!>              I think that this may be more of a time hole than it is worth.  I have checked
!>              the file attributes by hand for both creating files from scratch and appending
!>              clamping data to files that have global attributes for clamping.
! call exit(0)
! 
! call end_io_filenames(file_info)
! file_info = initialize_filenames(ens_handle, overwrite_state_input=.true.)
! 
! ! write out all files possible including restarts, mean, sd, and prior/posterior inflation files
! call write_state(ens_handle, file_info, prior_inflate_handle, post_inflate_handle, skip_diag_files())

! finalize test_cf_conventions
call error_handler(E_MSG,'test_cf_conventions','Finished successfully.',source,revision,revdate)
call finalize_mpi_utilities()

! end of main code


contains

!----------------------------------------------------------------------

subroutine initialize_module

  call initialize_mpi_utilities('test_cf_conventions')
  call register_module(source, revision, revdate)
  call static_init_assim_model()
  module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------

function initialize_filenames(ensemble_handle, overwrite_state_input) result(file_handle)
type(ensemble_type), intent(inout) :: ensemble_handle
logical,             intent(in)    :: overwrite_state_input
type(file_info_type) :: file_handle

logical :: single_restart_file_in  = .false.
logical :: single_restart_file_out = .false.
logical :: use_restart_list        = .false.
logical :: output_restart          = .true.
logical :: output_restart_mean     = .true.
logical :: add_domain_extension    = .false.
! logical :: perturb_from_single_instance = .false.

character(len=512) :: restart_list_file(10) = 'null'
character(len=129) :: inf_in_file_name(2)   = 'inf_in'
character(len=129) :: inf_out_file_name(2)  = 'inf_out'

character(len=129) :: restart_in_file_name  = "cf_test"
character(len=129) :: restart_out_file_name = "cf_test_out"

logical :: direct_netcdf_read  = .true.
logical :: direct_netcdf_write = .true.

file_handle = io_filenames_init(ensemble_handle, single_restart_file_in, single_restart_file_out, &
                restart_in_file_name, restart_out_file_name, output_restart, direct_netcdf_read, &
                direct_netcdf_write, output_restart_mean, add_domain_extension, use_restart_list, &
                restart_list_file, overwrite_state_input, inf_in_file_name, inf_out_file_name)

end function initialize_filenames

!----------------------------------------------------------------------

subroutine initialize_adaptive_inflate(ensemble_handle, prior_inflate, post_inflate)
type(ensemble_type),         intent(inout) :: ensemble_handle
type(adaptive_inflate_type), intent(inout) :: prior_inflate, post_inflate

! Inflation namelist entries follow, first entry for prior, second for posterior
! inf_flavor is 0:none, 1:obs space, 2: varying state space, 3: fixed state_space
integer              :: inf_flavor(2)             = 2
logical              :: inf_initial_from_restart(2)    = .false.
logical              :: inf_sd_initial_from_restart(2) = .false.

! old way
logical              :: inf_output_restart(2)     = .true.
! new way
!logical              :: inf_output_prior(2) = .false. ! mean sd
!logical              :: inf_output_post(2)  = .false. ! mean sd

logical              :: inf_deterministic(2)      = .true.
character(len = 129) :: inf_in_file_name(2)       = 'inf_in',    &
                        inf_out_file_name(2)      = 'inf_out',   &
                        inf_diag_file_name(2)     = 'inf_diag'
real(r8)             :: inf_initial(2)            = 1.0_r8
real(r8)             :: inf_sd_initial(2)         = 0.0_r8
!real(r8)             :: inf_damping(2)            = 1.0_r8
real(r8)             :: inf_lower_bound(2)        = 1.0_r8
real(r8)             :: inf_upper_bound(2)        = 1000000.0_r8
real(r8)             :: inf_sd_lower_bound(2)     = 0.0_r8

logical              :: allow_missing             = .false.

inf_out_file_name(1) = 'inf_out1'
inf_out_file_name(2) = 'inf_out2'

! Initialize the adaptive inflation module
call adaptive_inflate_init(prior_inflate, inf_flavor(1), inf_initial_from_restart(1), &
   inf_sd_initial_from_restart(1), inf_output_restart(1), inf_deterministic(1),       &
   inf_in_file_name(1), inf_out_file_name(1), inf_diag_file_name(1), inf_initial(1),  &
   inf_sd_initial(1), inf_lower_bound(1), inf_upper_bound(1), inf_sd_lower_bound(1),  &
   ensemble_handle, PRIOR_INF_COPY, PRIOR_INF_SD_COPY, allow_missing, 'Prior')

call adaptive_inflate_init(post_inflate, inf_flavor(2), inf_initial_from_restart(2),  &
   inf_sd_initial_from_restart(2), inf_output_restart(2), inf_deterministic(2),       &
   inf_in_file_name(2), inf_out_file_name(2), inf_diag_file_name(2), inf_initial(2),  &
   inf_sd_initial(2), inf_lower_bound(2), inf_upper_bound(2), inf_sd_lower_bound(2),  &
   ensemble_handle, POST_INF_COPY, POST_INF_SD_COPY, allow_missing, 'Posterior')

call adaptive_inflate_end(prior_inflate, ensemble_handle, PRIOR_INF_COPY, PRIOR_INF_SD_COPY)
call adaptive_inflate_end(post_inflate,  ensemble_handle, POST_INF_COPY, POST_INF_SD_COPY)

end subroutine initialize_adaptive_inflate

!----------------------------------------------------------------------

subroutine initialize_copy_numbers(ens_size)
integer, intent(in) :: ens_size

ENS_MEAN_COPY           = ens_size + 1
ENS_SD_COPY             = ens_size + 2
PRIOR_INF_COPY          = ens_size + 3
PRIOR_INF_SD_COPY       = ens_size + 4
POST_INF_COPY           = ens_size + 5
POST_INF_SD_COPY        = ens_size + 6
SPARE_PRIOR_MEAN        = ens_size + 7
SPARE_PRIOR_SPREAD      = ens_size + 8
SPARE_PRIOR_INF_MEAN    = ens_size + 9
SPARE_PRIOR_INF_SPREAD  = ens_size + 10
SPARE_POST_INF_MEAN     = ens_size + 11
SPARE_POST_INF_SPREAD   = ens_size + 12

end subroutine initialize_copy_numbers

!----------------------------------------------------------------------

subroutine initialize_diagnostics(PriorStateUnit, PosteriorStateUnit )
type(netcdf_file_type), intent(inout) :: PriorStateUnit, PosteriorStateUnit

integer :: i, ensemble_offset, num_state_copies
character(len=metadatalength) :: state_meta(num_output_state_members + 4)

! Section for state variables + other generated data stored with them.

! Ensemble mean goes first 
num_state_copies = num_output_state_members + 2
output_state_mean_index = 1
state_meta(output_state_mean_index) = 'ensemble mean'

! Ensemble spread goes second
output_state_spread_index = 2
state_meta(output_state_spread_index) = 'ensemble spread'
! Compute starting point for ensemble member output
ensemble_offset = 2

! Set up the metadata for the output state diagnostic files
do i = 1, ens_size
   write(state_meta(i + ensemble_offset), '(a15, 1x, i6)') 'ensemble member', i
end do

! Next two slots are for inflation mean and sd metadata
! To avoid writing out inflation values to the Prior and Posterior netcdf files,
! set output_inflation to false in the filter section of input.nml 
if(output_inflation) then
   num_state_copies = num_state_copies + 2
   state_meta(num_state_copies-1) = 'inflation mean'
   state_meta(num_state_copies)   = 'inflation sd'
endif

! Set up diagnostic output for model state
! All task call init and finalize diag_output.  The choice can then be made 
! in state_space_diag_mod to use a collective call (e.g. pnetcdf) or not.
PriorStateUnit     = init_diag_output('Prior_Diag', &
                     'prior ensemble state', num_state_copies, state_meta)
PosteriorStateUnit = init_diag_output('Posterior_Diag', &
                     'posterior ensemble state', num_state_copies, state_meta)

end subroutine initialize_diagnostics

!------------------------------------------------------------------
!> Copy the current mean, sd, inf_mean, inf_sd to spare copies
!> Assuming that if the spare copy is there you should fill it
subroutine store_prior(ens_handle)

type(ensemble_type), intent(inout) :: ens_handle

ens_handle%copies(SPARE_PRIOR_MEAN, :)       = ens_handle%copies(ENS_MEAN_COPY, :)
ens_handle%copies(SPARE_PRIOR_SPREAD, :)     = ens_handle%copies(ENS_SD_COPY, :)
ens_handle%copies(SPARE_PRIOR_INF_MEAN, :)   = ens_handle%copies(PRIOR_INF_COPY, :)
ens_handle%copies(SPARE_PRIOR_INF_SPREAD, :) = ens_handle%copies(PRIOR_INF_SD_COPY, :)

end subroutine store_prior


end program

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_test_read_write_restarts_dir/io/test/test_cf_conventions.f90 $
! $Id: test_cf_conventions.f90 9553 2016-01-20 17:26:41Z hendric $
! $Revision: 9553 $
! $Date: 2016-01-20 10:26:41 -0700 (Wed, 20 Jan 2016) $
