! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!>  A variety of operations required by assimilation.
module assim_tools_mod

!> \defgroup assim_tools assim_tools_mod
!>
!> @{
use      types_mod,       only : r8, i8, PI, missing_r8

use    options_mod,       only : get_missing_ok_status

use  utilities_mod,       only : file_exist, get_unit, check_namelist_read, do_output,    &
                                 find_namelist_in_file, error_handler,   &
                                 E_ERR, E_MSG, nmlfileunit, do_nml_file, do_nml_term,     &
                                 open_file, close_file, timestamp
use       sort_mod,       only : index_sort
use random_seq_mod,       only : random_seq_type, random_gaussian, init_random_seq,       &
                                 random_uniform

use obs_sequence_mod,     only : obs_sequence_type, obs_type, get_num_copies, get_num_qc, &
                                 init_obs, get_obs_from_key, get_obs_def, get_obs_values, &
                                 destroy_obs

use          obs_def_mod, only : obs_def_type, get_obs_def_location, get_obs_def_time,    &
                                 get_obs_def_error_variance, get_obs_def_type_of_obs

use         obs_kind_mod, only : get_num_types_of_obs, get_index_for_type_of_obs,                   &
                                 get_quantity_for_type_of_obs, assimilate_this_type_of_obs

use       cov_cutoff_mod, only : comp_cov_factor

use       reg_factor_mod, only : comp_reg_factor

use       obs_impact_mod, only : allocate_impact_table, read_impact_table, free_impact_table

use sampling_error_correction_mod, only : get_sampling_error_table_size, &
                                          read_sampling_error_correction

use         location_mod, only : location_type, get_close_type, query_location,           &
                                 operator(==), set_location_missing, write_location,      &
                                 LocationDims, is_vertical, vertical_localization_on,     &
                                 set_vertical, has_vertical_choice, get_close_init,       &
                                 get_vertical_localization_coord, get_close_destroy,      &
                                 set_vertical_localization_coord

use ensemble_manager_mod, only : ensemble_type, get_my_num_vars, get_my_vars,             &
                                 compute_copy_mean_var, get_var_owner_index,              &
                                 map_pe_to_task

use mpi_utilities_mod,    only : my_task_id, broadcast_send, broadcast_recv,              &
                                 sum_across_tasks, task_count, start_mpi_timer,           &
                                 read_mpi_timer, task_sync

use adaptive_inflate_mod, only : do_obs_inflate,  do_single_ss_inflate, do_ss_inflate,    &
                                 do_varying_ss_inflate,                                   &
                                 update_inflation, update_single_state_space_inflation,   &
                                 update_varying_state_space_inflation,                    &
                                 inflate_ens, adaptive_inflate_type,                      &
                                 deterministic_inflate, solve_quadratic

use time_manager_mod,     only : time_type, get_time

use assim_model_mod,      only : get_state_meta_data,                                     &
                                 get_close_obs,         get_close_state,                  &
                                 convert_vertical_obs,  convert_vertical_state

use distributed_state_mod, only : create_mean_window, free_mean_window

use quality_control_mod, only : good_dart_qc, DARTQC_FAILED_VERT_CONVERT

use probit_transform_mod, only : transform_to_probit, transform_from_probit, &
                                   transform_all_from_probit

use normal_distribution_mod, only : normal_cdf, inv_weighted_normal_cdf

use algorithm_info_mod, only : probit_dist_info, obs_inc_info, EAKF, ENKF, &
                               BOUNDED_NORMAL_RHF, UNBOUNDED_RHF, GAMMA_FILTER, &
                               KERNEL, OBS_PARTICLE

use gamma_distribution_mod, only : gamma_cdf, inv_gamma_cdf, gamma_mn_var_to_shape_scale, &
                                   gamma_gamma_prod

use bnrh_distribution_mod, only   :  inv_bnrh_cdf, bnrh_cdf, inv_bnrh_cdf_like

use distribution_params_mod, only : distribution_params_type, deallocate_distribution_params
                               

implicit none
private

public :: filter_assim, &
          set_assim_tools_trace, &
          test_state_copies, &
          update_ens_from_weights

! Indicates if module initialization subroutine has been called yet
logical :: module_initialized = .false.

integer :: print_timestamps    = 0
integer :: print_trace_details = 0

! True if random sequence needs to be initialized
logical                :: first_inc_ran_call = .true.
type (random_seq_type) :: inc_ran_seq

integer                :: num_types = 0
real(r8), allocatable  :: cutoff_list(:)
logical                :: has_special_cutoffs
logical                :: close_obs_caching = .true.

! true if we have multiple vert choices and we're doing vertical localization
! (make it a local variable so we don't keep making subroutine calls)
logical                :: is_doing_vertical_conversion = .false.

character(len=512)     :: msgstring, msgstring2, msgstring3

! Need to read in table for off-line based sampling correction and store it
integer                :: sec_table_size
real(r8), allocatable  :: exp_true_correl(:), alpha(:)

! if adjust_obs_impact is true, read in triplets from the ascii file
! and fill this 2d impact table.
real(r8), allocatable  :: obs_impact_table(:,:)

character(len=*), parameter :: source = 'assim_tools_mod.f90'

!============================================================================

!---- namelist with default values

real(r8) :: cutoff                          = 0.2_r8
logical  :: sort_obs_inc                    = .true.
logical  :: spread_restoration              = .false.
logical  :: sampling_error_correction       = .false.
integer  :: adaptive_localization_threshold = -1
real(r8) :: adaptive_cutoff_floor           = 0.0_r8
integer  :: print_every_nth_obs             = 0

! since this is in the namelist, it has to have a fixed size.
integer, parameter   :: MAX_ITEMS = 300
character(len = 129) :: special_localization_obs_types(MAX_ITEMS)
real(r8)             :: special_localization_cutoffs(MAX_ITEMS)

logical              :: output_localization_diagnostics = .false.
character(len = 129) :: localization_diagnostics_file = "localization_diagnostics"

! Following only relevant for filter_kind = UNBOUNDED_RHF
logical  :: rectangular_quadrature          = .true.
logical  :: gaussian_likelihood_tails       = .false.

! False by default; if true, expect to read in an ascii table
! to adjust the impact of obs on other state vector and obs values.
logical            :: adjust_obs_impact  = .false.
character(len=256) :: obs_impact_filename = ''
logical            :: allow_any_impact_values = .false.

! These next two only affect models with multiple options
! for vertical localization:
!
! "convert_state" is false by default; it depends on the model
! what is faster - do the entire state up front and possibly
! do unneeded work, or do the conversion during the assimilation
! loop. we think this depends heavily on how much of the state
! is going to be adjusted by the obs.  for a global model
! we think false may be better; for a regional model with
! a lot of obs and full coverage true may be better.
!
! "convert_obs" is true by default; in general it seems to
! be better for each task to convert the obs vertical before
! going into the loop but again this depends on how many
! obs per task and whether the mean is distributed or
! replicated on each task.
logical :: convert_all_state_verticals_first = .false.
logical :: convert_all_obs_verticals_first   = .true.

! Not in the namelist; this var disables the experimental
! linear and spherical case code in the adaptive localization
! sections.  to try out the alternatives, set this to .false.
logical  :: only_area_adapt  = .true.

! Option to distribute the mean.  If 'false' each task will have a full
! copy of the ensemble mean, which speeds models doing vertical conversion.
! If 'true' the mean will be spread across all tasks which reduces the
! memory needed per task but requires communication if the mean is used
! for vertical conversion.  We have changed the default to be .false.
! compared to previous versions of this namelist item.
logical  :: distribute_mean  = .false.

namelist / assim_tools_nml / cutoff, sort_obs_inc,                         &
   spread_restoration, sampling_error_correction,                          &
   adaptive_localization_threshold, adaptive_cutoff_floor,                 &
   print_every_nth_obs, rectangular_quadrature, gaussian_likelihood_tails, &
   output_localization_diagnostics, localization_diagnostics_file,         &
   special_localization_obs_types, special_localization_cutoffs,           &
   distribute_mean, close_obs_caching,                                     &
   adjust_obs_impact, obs_impact_filename, allow_any_impact_values,        &
   convert_all_state_verticals_first, convert_all_obs_verticals_first

!============================================================================

contains

!-------------------------------------------------------------

subroutine assim_tools_init()

integer :: iunit, io, i, j
integer :: num_special_cutoff, type_index
logical :: cache_override = .false.


! do this up front
module_initialized = .true.

! give these guys initial values at runtime *before* we read
! in the namelist.  this is to help detect how many items are
! actually given in the namelist.
special_localization_obs_types(:)  = 'null'
special_localization_cutoffs(:)    =  missing_r8

! Read the namelist entry
call find_namelist_in_file("input.nml", "assim_tools_nml", iunit)
read(iunit, nml = assim_tools_nml, iostat = io)
call check_namelist_read(iunit, io, "assim_tools_nml")

! Write the namelist values to the log file
if (do_nml_file()) write(nmlfileunit, nml=assim_tools_nml)
if (do_nml_term()) write(     *     , nml=assim_tools_nml)

! Forcing distributed_mean for single processor.
! Note null_win_mod.f90 ignores distibute_mean.
if (task_count() == 1) distribute_mean = .true.

if(spread_restoration) then
   write(msgstring, *) 'The spread_restoration option is not supported in this version of ', &
                       'DART. Contact the DAReS team if this option is needed '
   call error_handler(E_ERR,'assim_tools_init:', msgstring, source)
endif

! allocate a list in all cases - even the ones where there is only
! a single cutoff value.  note that in spite of the name these
! are specific types (e.g. RADIOSONDE_TEMPERATURE, AIRCRAFT_TEMPERATURE)
! because that's what get_close() is passed.   and because i've confused
! myself several times -- we define generic kinds starting at 0, but
! the specific types are autogenerated and always start at 1.  so the
! cutoff list is never (0:num_types); it is always (num_types).
num_types = get_num_types_of_obs()
allocate(cutoff_list(num_types))
cutoff_list(:) = cutoff
has_special_cutoffs = .false.

! Go through special-treatment observation kinds, if any.
num_special_cutoff = 0
j = 0
do i = 1, MAX_ITEMS
   if(special_localization_obs_types(i) == 'null') exit
   if(special_localization_cutoffs(i) == MISSING_R8) then
      write(msgstring, *) 'cutoff value', i, ' is uninitialized.'
      call error_handler(E_ERR,'assim_tools_init:', &
                         'special cutoff namelist for types and distances do not match', &
                         source, &
                         text2='kind = '//trim(special_localization_obs_types(i)), &
                         text3=trim(msgstring))
   endif
   j = j + 1
enddo
num_special_cutoff = j

if (num_special_cutoff > 0) has_special_cutoffs = .true.

do i = 1, num_special_cutoff
   type_index = get_index_for_type_of_obs(special_localization_obs_types(i))
   if (type_index < 0) then
      write(msgstring, *) 'unrecognized TYPE_ in the special localization namelist:'
      call error_handler(E_ERR,'assim_tools_init:', msgstring, source, &
                         text2=trim(special_localization_obs_types(i)))
   endif
   cutoff_list(type_index) = special_localization_cutoffs(i)
end do

! cannot cache previous obs location if different obs types have different
! localization radii.  change it to false, and warn user why.
if (has_special_cutoffs .and. close_obs_caching) then
   cache_override = .true.
   close_obs_caching = .false.
endif

if(sampling_error_correction) then
   sec_table_size = get_sampling_error_table_size()
   allocate(exp_true_correl(sec_table_size), alpha(sec_table_size))
   ! we can't read the table here because we don't have access to the ens_size
endif

is_doing_vertical_conversion = (has_vertical_choice() .and. vertical_localization_on())

call log_namelist_selections(num_special_cutoff, cache_override)

end subroutine assim_tools_init

!-------------------------------------------------------------

subroutine filter_assim(ens_handle, obs_ens_handle, obs_seq, keys,           &
   ens_size, num_groups, obs_val_index, inflate, ENS_MEAN_COPY, ENS_SD_COPY, &
   ENS_INF_COPY, ENS_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY,          &
   OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END, OBS_PRIOR_VAR_START,            &
   OBS_PRIOR_VAR_END, inflate_only)

type(ensemble_type),         intent(inout) :: ens_handle, obs_ens_handle
type(obs_sequence_type),     intent(in)    :: obs_seq
integer,                     intent(in)    :: keys(:)
integer,                     intent(in)    :: ens_size, num_groups, obs_val_index
! JLA: At present, this only needs to be inout because of the possible use of
! non-determinstic obs_space adaptive inflation that is not currently supported.
! Implementing that would require communication of the info about the inflation
! values as each observation updated them.
type(adaptive_inflate_type), intent(inout) :: inflate
integer,                     intent(in)    :: ENS_MEAN_COPY, ENS_SD_COPY, ENS_INF_COPY
integer,                     intent(in)    :: ENS_INF_SD_COPY
integer,                     intent(in)    :: OBS_KEY_COPY, OBS_GLOBAL_QC_COPY
integer,                     intent(in)    :: OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END
integer,                     intent(in)    :: OBS_PRIOR_VAR_START, OBS_PRIOR_VAR_END
logical,                     intent(in)    :: inflate_only

! changed the ensemble sized things here to allocatable

real(r8) :: obs_prior(ens_size), obs_inc(ens_size)
real(r8) :: obs_post(ens_size), probit_obs_prior(ens_size), probit_obs_post(ens_size)
real(r8) :: final_factor
real(r8) :: net_a(num_groups), correl(num_groups)
real(r8) :: obs(1), obs_err_var, my_inflate, my_inflate_sd
real(r8) :: obs_qc, cutoff_rev, cutoff_orig
real(r8) :: orig_obs_prior_mean(num_groups), orig_obs_prior_var(num_groups)
real(r8) :: obs_prior_mean(num_groups), obs_prior_var(num_groups)
real(r8) :: vertvalue_obs_in_localization_coord, whichvert_real
real(r8), allocatable :: close_obs_dist(:)
real(r8), allocatable :: close_state_dist(:)

integer(i8) :: state_index
integer(i8), allocatable :: my_state_indx(:)
integer(i8), allocatable :: my_obs_indx(:)

integer :: my_num_obs, i, j, owner, owners_index, my_num_state
integer :: obs_mean_index, obs_var_index
integer :: grp_beg(num_groups), grp_end(num_groups), grp_size, grp_bot, grp_top, group
integer :: num_close_obs, obs_index, num_close_states
integer :: last_num_close_obs, last_num_close_states
integer :: base_obs_kind, base_obs_type, nth_obs
integer :: num_close_obs_cached, num_close_states_cached
integer :: num_close_obs_calls_made, num_close_states_calls_made
integer :: whichvert_obs_in_localization_coord
integer :: istatus, localization_unit
integer, allocatable :: close_obs_ind(:)
integer, allocatable :: close_state_ind(:)
integer, allocatable :: my_obs_kind(:)
integer, allocatable :: my_obs_type(:)
integer, allocatable :: my_state_kind(:)
integer, allocatable :: vstatus(:)

type(location_type)  :: base_obs_loc, last_base_obs_loc, last_base_states_loc
type(location_type)  :: dummyloc
type(location_type), allocatable :: my_obs_loc(:)
type(location_type), allocatable :: my_state_loc(:)

type(get_close_type) :: gc_obs, gc_state
type(obs_type)       :: observation
type(obs_def_type)   :: obs_def
type(time_type)      :: obs_time

logical :: allow_missing_in_state
logical :: local_single_ss_inflate
logical :: local_varying_ss_inflate
logical :: local_ss_inflate
logical :: local_obs_inflate

! Storage for normal probit conversion, keeps prior mean and sd for all state ensemble members
type(distribution_params_type) :: state_dist_params(ens_handle%my_num_vars)
type(distribution_params_type) :: obs_dist_params(obs_ens_handle%my_num_vars)
integer :: dist_for_state, dist_for_obs
type(distribution_params_type) :: temp_dist_params
logical  :: bounded_below, bounded_above
real(r8) :: lower_bound,   upper_bound
real(r8) :: probit_ens(ens_size)

! allocate rather than dump all this on the stack
allocate(close_obs_dist(     obs_ens_handle%my_num_vars), &
         close_obs_ind(      obs_ens_handle%my_num_vars), &
         vstatus(            obs_ens_handle%my_num_vars), &
         my_obs_indx(        obs_ens_handle%my_num_vars), &
         my_obs_kind(        obs_ens_handle%my_num_vars), &
         my_obs_type(        obs_ens_handle%my_num_vars), &
         my_obs_loc(         obs_ens_handle%my_num_vars))

allocate(close_state_dist(     ens_handle%my_num_vars), &
         close_state_ind(      ens_handle%my_num_vars), &
         my_state_indx(        ens_handle%my_num_vars), &
         my_state_kind(        ens_handle%my_num_vars), &
         my_state_loc(         ens_handle%my_num_vars))
! end alloc

! Initialize assim_tools_module if needed
if (.not. module_initialized) call assim_tools_init()

!HK make window for mpi one-sided communication
! used for vertical conversion in get_close_obs
! Need to give create_mean_window the mean copy
call create_mean_window(ens_handle, ENS_MEAN_COPY, distribute_mean)

! Open the localization diagnostics file
if(output_localization_diagnostics .and. my_task_id() == 0) &
  localization_unit = open_file(localization_diagnostics_file, action = 'append')

! For performance, make local copies of these settings which
! are really in the inflate derived type.
local_single_ss_inflate  = do_single_ss_inflate(inflate)
local_varying_ss_inflate = do_varying_ss_inflate(inflate)
local_ss_inflate         = do_ss_inflate(inflate)
local_obs_inflate        = do_obs_inflate(inflate)

! Default to printing nothing
nth_obs = -1

! Divide ensemble into num_groups groups.
! make sure the number of groups and ensemble size result in
! at least 2 members in each group (to avoid divide by 0) and
! that the groups all have the same number of members.
grp_size = ens_size / num_groups
if ((grp_size * num_groups) /= ens_size) then
   write(msgstring,  *) 'The number of ensemble members must divide into the number of groups evenly.'
   write(msgstring2, *) 'Ensemble size = ', ens_size, '  Number of groups = ', num_groups
   write(msgstring3, *) 'Change number of groups or ensemble size to avoid remainders.'
   call error_handler(E_ERR,'filter_assim:', msgstring, source, &
                         text2=msgstring2,text3=msgstring3)
endif
if (grp_size < 2) then
   write(msgstring,  *) 'There must be at least 2 ensemble members in each group.'
   write(msgstring2, *) 'Ensemble size = ', ens_size, '  Number of groups = ', num_groups
   write(msgstring3, *) 'results in < 2 members/group.  Decrease number of groups or increase ensemble size'
   call error_handler(E_ERR,'filter_assim:', msgstring, source, &
                         text2=msgstring2,text3=msgstring3)
endif
do group = 1, num_groups
   grp_beg(group) = (group - 1) * grp_size + 1
   grp_end(group) = grp_beg(group) + grp_size - 1
enddo

! Put initial value of state space inflation in copy normally used for SD
! This is to avoid weird storage footprint in filter
ens_handle%copies(ENS_SD_COPY, :) = ens_handle%copies(ENS_INF_COPY, :)

! For single state or obs space inflation, the inflation is like a token
! Gets passed from the processor with a given obs on to the next
if(local_single_ss_inflate) then
   my_inflate    = ens_handle%copies(ENS_INF_COPY,    1)
   my_inflate_sd = ens_handle%copies(ENS_INF_SD_COPY, 1)
end if

! Get info on my number and indices for obs
my_num_obs = get_my_num_vars(obs_ens_handle)
call get_my_vars(obs_ens_handle, my_obs_indx)

! Construct an observation temporary
call init_obs(observation, get_num_copies(obs_seq), get_num_qc(obs_seq))

! Get the locations for all of my observations
! HK I would like to move this to before the calculation of the forward operator so you could
! overwrite the vertical location with the required localization vertical coordinate when you
! do the forward operator calculation
call get_my_obs_loc(obs_ens_handle, obs_seq, keys, my_obs_loc, my_obs_kind, my_obs_type, obs_time)

if (convert_all_obs_verticals_first .and. is_doing_vertical_conversion) then
   ! convert the vertical of all my observations to the localization coordinate
   if (obs_ens_handle%my_num_vars > 0) then
      call convert_vertical_obs(ens_handle, obs_ens_handle%my_num_vars, my_obs_loc, &
                                my_obs_kind, my_obs_type, get_vertical_localization_coord(), vstatus)
      do i = 1, obs_ens_handle%my_num_vars
         if (good_dart_qc(nint(obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, i)))) then
            !> @todo Can I just use the OBS_GLOBAL_QC_COPY? Is it ok to skip the loop?
            if (vstatus(i) /= 0) obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, i) = DARTQC_FAILED_VERT_CONVERT
         endif
      enddo
   endif 
endif

! Get info on my number and indices for state
my_num_state = get_my_num_vars(ens_handle)
call get_my_vars(ens_handle, my_state_indx)

! Get the location and kind of all my state variables
do i = 1, ens_handle%my_num_vars
   call get_state_meta_data(my_state_indx(i), my_state_loc(i), my_state_kind(i))

   ! Need to specify what kind of prior to use for each
   call probit_dist_info(my_state_kind(i), .true., .false., dist_for_state, &
      bounded_below, bounded_above, lower_bound, upper_bound)

   ! Convert all my state variables to appropriate probit space
   call transform_to_probit(ens_size, ens_handle%copies(1:ens_size, i), dist_for_state, &
      state_dist_params(i), probit_ens, .false., &
      bounded_below, bounded_above, lower_bound, upper_bound)
   ens_handle%copies(1:ens_size, i) = probit_ens
end do

!> optionally convert all state location verticals
if (convert_all_state_verticals_first .and. is_doing_vertical_conversion) then
   if (ens_handle%my_num_vars > 0) then
      call convert_vertical_state(ens_handle, ens_handle%my_num_vars, my_state_loc, my_state_kind,  &
                                  my_state_indx, get_vertical_localization_coord(), istatus)
   endif
endif

! Get mean and variance of each group's observation priors for adaptive inflation
! Important that these be from before any observations have been used
if(local_ss_inflate) then
   do group = 1, num_groups
      obs_mean_index = OBS_PRIOR_MEAN_START + group - 1
      obs_var_index  = OBS_PRIOR_VAR_START  + group - 1
         call compute_copy_mean_var(obs_ens_handle, grp_beg(group), grp_end(group), &
           obs_mean_index, obs_var_index)
   end do
endif

! Have gotten the mean and variance from original ensembles, can convert all my obs to probit
! CAN WE DO THE ADAPTIVE INFLATION ENTIRELY IN PROBIT SPACE TO MAKE IT DISTRIBUTION INDEPENDENT????
! WOULD NEED AN OBSERVATION ERROR VARIANCE IN PROBIT SPACE SOMEHOW. IS THAT POSSIBLE???

do i = 1, my_num_obs
   obs_qc = obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, i)
   ! Only do conversion of qc if forward operator is good
   if(nint(obs_qc) == 0) then
      ! Need to specify what kind of prior to use for each
      call probit_dist_info(my_obs_kind(i), .false., .false., dist_for_obs, &
         bounded_below, bounded_above, lower_bound, upper_bound)
   
      ! Convert all my obs (extended state) variables to appropriate probit space
      call transform_to_probit(ens_size, obs_ens_handle%copies(1:ens_size, i), dist_for_obs, &
         obs_dist_params(i), probit_ens, .false., &
         bounded_below, bounded_above, lower_bound, upper_bound)
      obs_ens_handle%copies(1:ens_size, i) = probit_ens
   endif
end do

! Initialize the method for getting state variables close to a given ob on my process
if (has_special_cutoffs) then
   call get_close_init(gc_state, my_num_state, 2.0_r8*cutoff, my_state_loc, 2.0_r8*cutoff_list)
else
   call get_close_init(gc_state, my_num_state, 2.0_r8*cutoff, my_state_loc)
endif

! Initialize the method for getting obs close to a given ob on my process
if (has_special_cutoffs) then
   call get_close_init(gc_obs, my_num_obs, 2.0_r8*cutoff, my_obs_loc, 2.0_r8*cutoff_list)
else
   call get_close_init(gc_obs, my_num_obs, 2.0_r8*cutoff, my_obs_loc)
endif

if (close_obs_caching) then
   ! Initialize last obs and state get_close lookups, to take advantage below
   ! of sequential observations at the same location (e.g. U,V, possibly T,Q)
   ! (this is getting long enough it probably should go into a subroutine. nsc.)
   last_base_obs_loc           = set_location_missing()
   last_base_states_loc        = set_location_missing()
   last_num_close_obs          = -1
   last_num_close_states       = -1
   num_close_obs_cached        = 0
   num_close_states_cached     = 0
   num_close_obs_calls_made    = 0
   num_close_states_calls_made = 0
endif

allow_missing_in_state = get_missing_ok_status()

! Loop through all the (global) observations sequentially
SEQUENTIAL_OBS: do i = 1, obs_ens_handle%num_vars
   ! Some compilers do not like mod by 0, so test first.
   if (print_every_nth_obs > 0) nth_obs = mod(i, print_every_nth_obs)

   ! If requested, print out a message every Nth observation
   ! to indicate progress is being made and to allow estimates
   ! of how long the assim will take.
   if (nth_obs == 0) then
      write(msgstring, '(2(A,I8))') 'Processing observation ', i, &
                                         ' of ', obs_ens_handle%num_vars
      if (print_timestamps == 0) then
         call error_handler(E_MSG,'filter_assim',msgstring)
      else
         call timestamp(trim(msgstring), pos="brief")
      endif
   endif

   ! Every pe has information about the global obs sequence
   call get_obs_from_key(obs_seq, keys(i), observation)
   call get_obs_def(observation, obs_def)
   base_obs_loc = get_obs_def_location(obs_def)
   obs_err_var = get_obs_def_error_variance(obs_def)
   base_obs_type = get_obs_def_type_of_obs(obs_def)
   if (base_obs_type > 0) then
      base_obs_kind = get_quantity_for_type_of_obs(base_obs_type)
   else
      call get_state_meta_data(-1 * int(base_obs_type,i8), dummyloc, base_obs_kind)  ! identity obs
   endif
   ! Get the value of the observation
   call get_obs_values(observation, obs, obs_val_index)

   ! Find out who has this observation and where it is
   call get_var_owner_index(ens_handle, int(i,i8), owner, owners_index)

   ! Following block is done only by the owner of this observation
   !-----------------------------------------------------------------------
   if(ens_handle%my_pe == owner) then
      ! each task has its own subset of all obs.  if they were converted in the
      ! vertical up above, then we need to broadcast the new values to all the other
      ! tasks so they're computing the right distances when applying the increments.
      if (is_doing_vertical_conversion) then
         vertvalue_obs_in_localization_coord = query_location(my_obs_loc(owners_index), "VLOC")
         whichvert_obs_in_localization_coord = query_location(my_obs_loc(owners_index), "WHICH_VERT")
      else
         vertvalue_obs_in_localization_coord = 0.0_r8
         whichvert_obs_in_localization_coord = 0
      endif

      obs_qc = obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, owners_index)

      ! Only value of 0 for DART QC field should be assimilated
      IF_QC_IS_OKAY: if(nint(obs_qc) ==0) then
         ! Note that these are before DA starts, so can be different from current obs_prior
         orig_obs_prior_mean = obs_ens_handle%copies(OBS_PRIOR_MEAN_START: &
            OBS_PRIOR_MEAN_END, owners_index)
         orig_obs_prior_var  = obs_ens_handle%copies(OBS_PRIOR_VAR_START:  &
            OBS_PRIOR_VAR_END, owners_index)

         ! If QC is okay, convert this observation ensemble from probit to regular space
         call transform_from_probit(ens_size, obs_ens_handle%copies(1:ens_size, owners_index) , &
            obs_dist_params(owners_index), obs_ens_handle%copies(1:ens_size, owners_index))

         obs_prior = obs_ens_handle%copies(1:ens_size, owners_index)
      endif IF_QC_IS_OKAY

      !Broadcast the info from this obs to all other processes
      ! orig_obs_prior_mean and orig_obs_prior_var only used with adaptive inflation
      ! my_inflate and my_inflate_sd only used with single state space inflation
      ! vertvalue_obs_in_localization_coord and whichvert_real only used for vertical
      ! coordinate transformation
      whichvert_real = real(whichvert_obs_in_localization_coord, r8)
      call broadcast_send(map_pe_to_task(ens_handle, owner), obs_prior,    &
         orig_obs_prior_mean, orig_obs_prior_var,                          &
         scalar1=obs_qc, scalar2=vertvalue_obs_in_localization_coord,      &
         scalar3=whichvert_real, scalar4=my_inflate, scalar5=my_inflate_sd)

   ! Next block is done by processes that do NOT own this observation
   !-----------------------------------------------------------------------
   else
      call broadcast_recv(map_pe_to_task(ens_handle, owner), obs_prior,    &
         orig_obs_prior_mean, orig_obs_prior_var,                          & 
         scalar1=obs_qc, scalar2=vertvalue_obs_in_localization_coord,      &
         scalar3=whichvert_real, scalar4=my_inflate, scalar5=my_inflate_sd)
      whichvert_obs_in_localization_coord = nint(whichvert_real)

   endif
   !-----------------------------------------------------------------------

   ! Everybody is doing this section, cycle if qc is bad
   if(nint(obs_qc) /= 0) cycle SEQUENTIAL_OBS

   !> all tasks must set the converted vertical values into the 'base' version of this loc
   !> because that's what we pass into the get_close_xxx() routines below.
   if (is_doing_vertical_conversion) &
      call set_vertical(base_obs_loc, vertvalue_obs_in_localization_coord, whichvert_obs_in_localization_coord)

   ! Compute observation space increments for each group
   do group = 1, num_groups
      grp_bot = grp_beg(group); grp_top = grp_end(group)
      call obs_increment(obs_prior(grp_bot:grp_top), grp_size, obs(1), &
         obs_err_var, base_obs_kind, obs_inc(grp_bot:grp_top), inflate, my_inflate,   &
         my_inflate_sd, net_a(group))
      obs_post(grp_bot:grp_top) = obs_prior(grp_bot:grp_top) + obs_inc(grp_bot:grp_top)

      ! Convert both the prior and posterior to probit space (efficiency for prior???)
      ! Running probit space with groups needs to be studied more carefully
      ! EFFICIENCY NOTE: FOR RHF, THE OBS_INCREMENT HAS TO DO A SORT
      ! THE POSTERIOR WOULD HAVE THE SAME RANK STATISTICS, SO THIS SORT WOULD BE THE SAME
      ! THE SECOND CONVERT_TO_PROBIT CAN BE MUCH MORE EFFICIENT USING A SORT
      ! SHOULD FIGURE OUT A WAY TO PASS THE SORT ORDER
      ! NOTE 2: THIS CONVERSION IS USING THE INFO FROM THE CURRENT (UPDATED) PRIOR ENSEMBLE. THIS
      ! IS GENERALLY GOING TO BE A DIFFERENT PROBIT TRANSFORMED ENSEMBLE THAN THE ONE THAT WAS JUST
      ! CONVERTED FROM PROBIT SPACE BY THE PROCESS THAT OWNS THIS OBSERVATION. 

      ! Need to specify what kind of prior to use for obs being assimilated
      call probit_dist_info(base_obs_kind, .false., .false., dist_for_obs, &
         bounded_below, bounded_above, lower_bound, upper_bound)

      ! Convert the prior and posterior for this observation to probit space
      call transform_to_probit(grp_size, obs_prior(grp_bot:grp_top), dist_for_obs, &
         temp_dist_params, probit_obs_prior(grp_bot:grp_top), .false., &
         bounded_below, bounded_above, lower_bound, upper_bound)
      call transform_to_probit(grp_size, obs_post(grp_bot:grp_top), dist_for_obs, &
         temp_dist_params, probit_obs_post(grp_bot:grp_top), .true., &
         bounded_below, bounded_above, lower_bound, upper_bound)
      ! Free up the storage used for this transform
      call deallocate_distribution_params(temp_dist_params)

      ! Copy back into original storage
      obs_prior(grp_bot:grp_top) = probit_obs_prior(grp_bot:grp_top)
      obs_post(grp_bot:grp_top) = probit_obs_post(grp_bot:grp_top)
      ! Recompute obs_inc in probit space
      obs_inc(grp_bot:grp_top) = obs_post(grp_bot:grp_top) - obs_prior(grp_bot:grp_top)

      ! Also compute prior mean and variance of obs for efficiency here
      obs_prior_mean(group) = sum(obs_prior(grp_bot:grp_top)) / grp_size
      obs_prior_var(group) = sum((obs_prior(grp_bot:grp_top) - obs_prior_mean(group))**2) / &
         (grp_size - 1)
      if (obs_prior_var(group) < 0.0_r8) obs_prior_var(group) = 0.0_r8
   end do

   ! Compute updated values for single state space inflation
   if(local_single_ss_inflate) then
      ! Update for each group separately
      do group = 1, num_groups
         call update_single_state_space_inflation(inflate, my_inflate, my_inflate_sd, &
            ens_handle%copies(ENS_SD_COPY, 1), orig_obs_prior_mean(group), &
            orig_obs_prior_var(group), obs(1), obs_err_var, grp_size, inflate_only)
      end do
   endif
  
   ! Adaptive localization needs number of other observations within localization radius.
   ! Do get_close_obs first, even though state space increments are computed before obs increments.
   call  get_close_obs_cached(gc_obs, base_obs_loc, base_obs_type,      &
      my_obs_loc, my_obs_kind, my_obs_type, num_close_obs, close_obs_ind, close_obs_dist,  &
      ens_handle, last_base_obs_loc, last_num_close_obs, num_close_obs_cached,             &
      num_close_obs_calls_made)

   ! set the cutoff default, keep a copy of the original value, and avoid
   ! looking up the cutoff in a list if the incoming obs is an identity ob
   ! (and therefore has a negative kind).  specific types can never be 0;
   ! generic kinds (not used here) start their numbering at 0 instead of 1.
   if (base_obs_type > 0) then
      cutoff_orig = cutoff_list(base_obs_type)
   else
      cutoff_orig = cutoff
   endif

   ! JLA, could also cache for adaptive_localization which may be expensive?
   call adaptive_localization_and_diags(cutoff_orig, cutoff_rev, adaptive_localization_threshold, &
      adaptive_cutoff_floor, num_close_obs, close_obs_ind, close_obs_dist, my_obs_type, &
      i, base_obs_loc, obs_def, localization_unit)

   ! Find state variables on my process that are close to observation being assimilated
   call  get_close_state_cached(gc_state, base_obs_loc, base_obs_type,      &
      my_state_loc, my_state_kind, my_state_indx, num_close_states, close_state_ind, close_state_dist,  &
      ens_handle, last_base_states_loc, last_num_close_states, num_close_states_cached,              &
      num_close_states_calls_made)
   !call test_close_obs_dist(close_state_dist, num_close_states, i)

   ! Loop through to update each of my state variables that is potentially close
   STATE_UPDATE: do j = 1, num_close_states
      state_index = close_state_ind(j)

      if ( allow_missing_in_state ) then
         ! Don't allow update of state ensemble with any missing values
         if (any(ens_handle%copies(1:ens_size, state_index) == MISSING_R8)) cycle STATE_UPDATE
      endif

      ! Compute the covariance localization and adjust_obs_impact factors (module storage)
      final_factor = cov_and_impact_factors(base_obs_loc, base_obs_type, my_state_loc(state_index), &
         my_state_kind(state_index), close_state_dist(j), cutoff_rev)

      if(final_factor <= 0.0_r8) cycle STATE_UPDATE
      
      call obs_updates_ens(ens_size, num_groups, ens_handle%copies(1:ens_size, state_index), &
         my_state_loc(state_index), my_state_kind(state_index), obs_prior, obs_inc, &
         obs_prior_mean, obs_prior_var, base_obs_loc, base_obs_type, obs_time, &
         net_a, grp_size, grp_beg, grp_end, i, &
         my_state_indx(state_index), final_factor, correl, local_varying_ss_inflate, inflate_only)

      ! Compute spatially-varying state space inflation
      if(local_varying_ss_inflate) then
         do group = 1, num_groups
            call update_varying_state_space_inflation(inflate,                     &
               ens_handle%copies(ENS_INF_COPY, state_index),                       &
               ens_handle%copies(ENS_INF_SD_COPY, state_index),                    &
               ens_handle%copies(ENS_SD_COPY, state_index),                        &
               orig_obs_prior_mean(group), orig_obs_prior_var(group), obs(1),      &
               obs_err_var, grp_size, final_factor, correl(group), inflate_only)
         end do
      endif
   end do STATE_UPDATE

   if(.not. inflate_only) then
      ! Now everybody updates their obs priors (only ones after this one)
      OBS_UPDATE: do j = 1, num_close_obs
         obs_index = close_obs_ind(j)

         ! Only have to update obs that have not yet been used
         if(my_obs_indx(obs_index) > i) then

            ! If forward observation operator failed, no need to update unassimilated observations
            if (any(obs_ens_handle%copies(1:ens_size, obs_index) == MISSING_R8)) cycle OBS_UPDATE

         ! Compute the covariance localization and adjust_obs_impact factors (module storage)
            final_factor = cov_and_impact_factors(base_obs_loc, base_obs_type, my_obs_loc(obs_index), &
            my_obs_kind(obs_index), close_obs_dist(j), cutoff_rev)

            if(final_factor <= 0.0_r8) cycle OBS_UPDATE

            call obs_updates_ens(ens_size, num_groups, obs_ens_handle%copies(1:ens_size, obs_index), &
               my_obs_loc(obs_index), my_obs_kind(obs_index), obs_prior, obs_inc, &
               obs_prior_mean, obs_prior_var, base_obs_loc, base_obs_type, obs_time, &
               net_a, grp_size, grp_beg, grp_end, i, &
               -1*my_obs_indx(obs_index), final_factor, correl, .false., inflate_only)
         endif
      end do OBS_UPDATE
   endif
end do SEQUENTIAL_OBS

! Do the inverse probit transform for state variables
call transform_all_from_probit(ens_size, ens_handle%my_num_vars, ens_handle%copies, &
   state_dist_params, ens_handle%copies)

! Every pe needs to get the current my_inflate and my_inflate_sd back
if(local_single_ss_inflate) then
   ens_handle%copies(ENS_INF_COPY, :) = my_inflate
   ens_handle%copies(ENS_INF_SD_COPY, :) = my_inflate_sd
end if

! Free up the storage
call destroy_obs(observation)
call get_close_destroy(gc_state)
call get_close_destroy(gc_obs)

! do some stats - being aware that unless we do a reduce() operation
! this is going to be per-task.  so only print if something interesting
! shows up in the stats?  maybe it would be worth a reduce() call here?

! Assure user we have done something
if (print_trace_details >= 0) then
   write(msgstring, '(A,I8,A)') 'Processed', obs_ens_handle%num_vars, ' total observations'
   call error_handler(E_MSG,'filter_assim:',msgstring)
endif

! diagnostics for stats on saving calls by remembering obs at the same location.
! change .true. to .false. in the line below to remove the output completely.
if (close_obs_caching) then
   if (num_close_obs_cached > 0 .and. do_output()) then
      print *, "Total number of calls made    to get_close_obs for obs/states:    ", &
                num_close_obs_calls_made + num_close_states_calls_made
      print *, "Total number of calls avoided to get_close_obs for obs/states:    ", &
                num_close_obs_cached + num_close_states_cached
      if (num_close_obs_cached+num_close_obs_calls_made+ &
          num_close_states_cached+num_close_states_calls_made > 0) then
         print *, "Percent saved: ", 100.0_r8 * &
                   (real(num_close_obs_cached+num_close_states_cached, r8) /  &
                   (num_close_obs_calls_made+num_close_obs_cached +           &
                    num_close_states_calls_made+num_close_states_cached))
      endif
   endif
endif

!call test_state_copies(ens_handle, 'end')

! Close the localization diagnostics file
if(output_localization_diagnostics .and. my_task_id() == 0) call close_file(localization_unit)

! get rid of mpi window
call free_mean_window()

! deallocate space
deallocate(close_obs_dist,      &
           my_obs_indx,         &
           my_obs_kind,         &
           my_obs_type,         &
           close_obs_ind,       &
           vstatus,             &
           my_obs_loc)

deallocate(close_state_dist,      &
           my_state_indx,         &
           close_state_ind,       &
           my_state_kind,         &
           my_state_loc)

! end dealloc

end subroutine filter_assim

!-------------------------------------------------------------

subroutine obs_increment(ens_in, ens_size, obs, obs_var, obs_kind, obs_inc, &
   inflate, my_cov_inflate, my_cov_inflate_sd, net_a)

! Given the ensemble prior for an observation, the observation, and
! the observation error variance, computes increments and adjusts
! observation space inflation values

integer,                     intent(in)    :: ens_size
real(r8),                    intent(in)    :: ens_in(ens_size), obs, obs_var
integer,                     intent(in)    :: obs_kind
real(r8),                    intent(out)   :: obs_inc(ens_size)
type(adaptive_inflate_type), intent(inout) :: inflate
real(r8),                    intent(inout) :: my_cov_inflate, my_cov_inflate_sd
real(r8),                    intent(out)   :: net_a

real(r8) :: ens(ens_size), inflate_inc(ens_size)
real(r8) :: prior_mean, prior_var, new_val(ens_size)
integer  :: i, ens_index(ens_size), new_index(ens_size)

real(r8) :: rel_weights(ens_size)

integer  :: filter_kind
logical  :: bounded_below, bounded_above
real(r8) :: lower_bound,   upper_bound

! Declarations for bounded rank histogram filter
real(r8) :: likelihood(ens_size), like_sum

! Copy the input ensemble to something that can be modified
ens = ens_in

! Null value of net spread change factor is 1.0
net_a = 0.0_r8

! Compute prior variance and mean from sample
prior_mean = sum(ens) / ens_size
prior_var  = sum((ens - prior_mean)**2) / (ens_size - 1)

! If observation space inflation is being done, compute the initial
! increments and update the inflation factor and its standard deviation
! as needed. my_cov_inflate < 0 means don't do any of this.
if(do_obs_inflate(inflate)) then
   ! If my_cov_inflate_sd is <= 0, just retain current my_cov_inflate setting
   if(my_cov_inflate_sd > 0.0_r8) &
      ! Gamma set to 1.0 because no distance for observation space
      call update_inflation(inflate, my_cov_inflate, my_cov_inflate_sd, prior_mean, &
         prior_var, ens_size, obs, obs_var, gamma_corr = 1.0_r8)

   ! Now inflate the ensemble and compute a preliminary inflation increment
   call inflate_ens(inflate, ens, prior_mean, my_cov_inflate, prior_var)
   ! Keep the increment due to inflation alone
   inflate_inc = ens - ens_in

   ! Need to recompute variance if non-deterministic inflation (mean is unchanged)
   if(.not. deterministic_inflate(inflate)) &
      prior_var  = sum((ens - prior_mean)**2) / (ens_size - 1)
endif

! Gets information about the increment method (filter_kind, bounds) for the current observation.
call obs_inc_info(obs_kind, filter_kind, bounded_below, bounded_above, &
                  lower_bound, upper_bound)

! The first three options in the next if block of code may be inappropriate for 
! some more general filters; need to revisit
! If obs_var == 0, delta function.  The mean becomes obs value with no spread.
! If prior_var == 0, obs has no effect.  The increments are 0.
! If both obs_var and prior_var == 0 there is no right thing to do, so Stop.
if ((obs_var == 0.0_r8) .and. (prior_var == 0.0_r8)) then

   ! fail if both obs variance and prior spreads are 0.
   write(msgstring,  *) 'Observation value is ', obs, ' ensemble mean value is ', prior_mean
   write(msgstring2, *) 'The observation has 0.0 error variance, and the ensemble members have 0.0 spread.'
   write(msgstring3, *) 'These require inconsistent actions and the algorithm cannot continue.'
   call error_handler(E_ERR, 'obs_increment', msgstring, &
           source, text2=msgstring2, text3=msgstring3)

else if (obs_var == 0.0_r8) then

   ! new mean is obs value, so increments are differences between obs
   ! value and current value.  after applying obs, all state will equal obs.
   obs_inc(:) = obs - ens

else if (prior_var == 0.0_r8) then

   ! if all state values are the same, nothing changes.
   obs_inc(:) = 0.0_r8

else

   ! Call the appropriate filter option to compute increments for ensemble
   ! note that at this point we've taken care of the cases where either the
   ! obs_var or the prior_var is 0, so the individual routines no longer need
   ! to have code to test for those cases.
   if(filter_kind == EAKF) then
      call obs_increment_eakf(ens, ens_size, prior_mean, prior_var, &
         obs, obs_var, obs_inc, net_a)
   else if(filter_kind == ENKF) then
      call obs_increment_enkf(ens, ens_size, prior_var, obs, obs_var, obs_inc)
   else if(filter_kind == KERNEL) then
      call obs_increment_kernel(ens, ens_size, obs, obs_var, obs_inc)
   else if(filter_kind == OBS_PARTICLE) then
      call obs_increment_particle(ens, ens_size, obs, obs_var, obs_inc)
   else if(filter_kind == UNBOUNDED_RHF) then
      call obs_increment_rank_histogram(ens, ens_size, prior_var, obs, obs_var, obs_inc)
   else if(filter_kind == GAMMA_FILTER) then
      call obs_increment_gamma(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
   !--------------------------------------------------------------------------
   else if(filter_kind == BOUNDED_NORMAL_RHF) then

      ! Use bounded normal likelihood; Could use an arbitrary likelihood
      do i = 1, ens_size
         likelihood(i) = get_truncated_normal_like(ens(i), obs, obs_var, &
            bounded_below, bounded_above, lower_bound, upper_bound)
      end do

      ! Normalize the likelihood here
      like_sum = sum(likelihood)
      ! If likelihood underflow, assume flat likelihood, so no increments
      if(like_sum <= 0.0_r8) then
         obs_inc = 0.0_r8
         return
      else
         likelihood = likelihood / like_sum
      endif

      call obs_increment_bounded_norm_rhf(ens, likelihood, ens_size, prior_var, &
         obs_inc, bounded_below, bounded_above, lower_bound, upper_bound)

      ! Do test of inversion for an uninformative likelihood
      !!!t_likelihood = 1.0
      !!!t_likelihood = t_likelihood / sum(t_likelihood)
      !!!call obs_increment_bounded_norm_rhf(ens, t_likelihood, ens_size, prior_var, &
         !!!obs_inc_temp, bounded_below, bounded_above, lower_bound, upper_bound)
      !!!if(maxval(abs(obs_inc_temp)) > 1e-3_r8) then
         !!!write(msgstring, *) 'Null increment tests exceed the threshold ', maxval(abs(obs_inc_temp))
         !!!call error_handler(E_ERR, 'obs_increment', msgstring, source)
      !!!endif

   !--------------------------------------------------------------------------
   else
      call error_handler(E_ERR,'obs_increment', &
              'Illegal value of filter_kind', source)
   endif
endif

! Add in the extra increments if doing observation space covariance inflation
if(do_obs_inflate(inflate)) obs_inc = obs_inc + inflate_inc

! Get the net change in spread if obs space inflation was used
if(do_obs_inflate(inflate)) net_a = net_a * sqrt(my_cov_inflate)

end subroutine obs_increment



subroutine obs_increment_gamma(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
!========================================================================
!
! Gamma version of obs increment. This demonstrates the updat

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), prior_mean, prior_var, obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: prior_shape, prior_scale, like_shape, like_scale, post_shape, post_scale
real(r8) :: q(ens_size), post(ens_size)
integer :: i

! Compute the prior quantiles of each ensemble member in the prior gamma distribution
call gamma_mn_var_to_shape_scale(prior_mean, prior_var, prior_shape, prior_scale)
do i = 1, ens_size
   q(i) = gamma_cdf(ens(i), prior_shape, prior_scale, .true., .false., 0.0_r8, missing_r8) 
end do

! Compute the statistics of the continous posterior distribution
call gamma_mn_var_to_shape_scale(obs, obs_var, like_shape, like_scale)
call gamma_gamma_prod(prior_shape, prior_scale, like_shape, like_scale, &
   post_shape, post_scale)

! Check for illegal values. This can occur if the distributions are getting too
! concentrated towards the bound
if(post_shape <= 0.0_r8) then
   write(msgstring, *) 'Posterior gamma shape is negative ', post_shape
   call error_handler(E_ERR, 'obs_increment_gamma', msgstring, source)
endif

! Now invert the quantiles with the posterior distribution
do i = 1, ens_size
   post(i) = inv_gamma_cdf(q(i), post_shape, post_scale, .true., .false., 0.0_r8, missing_r8)
end do

obs_inc = post - ens

end subroutine obs_increment_gamma



subroutine obs_increment_eakf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc, a)
!========================================================================
!
! EAKF version of obs increment

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), prior_mean, prior_var, obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
real(r8), intent(out) :: a

real(r8) :: new_mean, var_ratio

! Compute the new mean
var_ratio = obs_var / (prior_var + obs_var)
new_mean  = var_ratio * (prior_mean  + prior_var*obs / obs_var)

! Compute sd ratio and shift ensemble
a = sqrt(var_ratio)
obs_inc = a * (ens - prior_mean) + new_mean - ens

end subroutine obs_increment_eakf


subroutine obs_increment_bounded_norm_rhf(ens, ens_like, ens_size, prior_var, &
   obs_inc, bounded_below, bounded_above, lower_bound, upper_bound)
!------------------------------------------------------------------------
integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size)
real(r8), intent(in)  :: ens_like(ens_size)
real(r8), intent(in)  :: prior_var
real(r8), intent(out) :: obs_inc(ens_size)
logical,  intent(in)  :: bounded_below, bounded_above
real(r8), intent(in)  :: lower_bound,   upper_bound

! Does bounded RHF assuming that the prior in outer regions is part of a normal. 

real(r8) :: sort_ens(ens_size), sort_ens_like(ens_size)
real(r8) :: post(ens_size), sort_post(ens_size), q(ens_size)
real(r8) :: tail_amp_left,  tail_mean_left,  tail_sd_left
real(r8) :: tail_amp_right, tail_mean_right, tail_sd_right
logical  :: do_uniform_tail_left, do_uniform_tail_right
integer  :: i, sort_ind(ens_size)

! If all ensemble members are identical, this algorithm becomes undefined, so fail
if(prior_var <= 0.0_r8) then
      msgstring = 'Ensemble variance <= 0 '
      call error_handler(E_ERR, 'obs_increment_bounded_norm_rhf', msgstring, source)
endif

! Do an index sort of the ensemble members; Use prior info for efficiency in the future
call index_sort(ens, sort_ind, ens_size)

! Get the sorted ensemble
sort_ens = ens(sort_ind)

! Get the sorted likelihood
sort_ens_like = ens_like(sort_ind)

! Generate the prior information for a BNRH for this ensemble
call bnrh_cdf(ens, ens_size, bounded_below, bounded_above, lower_bound, upper_bound, &
   sort_ens, q, tail_amp_left,  tail_mean_left,  tail_sd_left,  do_uniform_tail_left,  &
   tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right)

! Invert the bnrh cdf after it is multiplied by the likelihood
call inv_bnrh_cdf_like(q, ens_size, sort_ens, &
   bounded_below, bounded_above, lower_bound, upper_bound, &
   tail_amp_left,  tail_mean_left,  tail_sd_left,  do_uniform_tail_left, &
   tail_amp_right, tail_mean_right, tail_sd_right, do_uniform_tail_right, post, &
   sort_ens_like)
! The posterior needs to be sorted to get increments; this can be done more efficiently
sort_post = post(sort_ind)

! These are increments for sorted ensemble; convert to increments for unsorted
do i = 1, ens_size
   obs_inc(sort_ind(i)) = sort_post(i) - ens(sort_ind(i))
   ! It may be possible, although apparently exceedingly unusual, to generate an increment
   ! here that when added back onto the prior leads to a posterior that is greater than 
   ! the bounds. Unclear if there is any direct way to fix this given that increments are
   ! being passed out.
end do

end subroutine obs_increment_bounded_norm_rhf




! Computes a normal or truncated normal (above and/or below) likelihood.
function get_truncated_normal_like(x, obs, obs_var, &
   bounded_below, bounded_above, lower_bound, upper_bound)
!------------------------------------------------------------------------
real(r8)             :: get_truncated_normal_like
real(r8), intent(in) :: x
real(r8), intent(in) :: obs, obs_var
logical,  intent(in) :: bounded_below, bounded_above
real(r8), intent(in) :: lower_bound,   upper_bound

real(r8) :: cdf(2), obs_sd, weight

! A zero observation error variance is a degenerate case
if(obs_var <= 0.0_r8) then
   if(x == obs) then
      get_truncated_normal_like = 1.0_r8
   else
      get_truncated_normal_like = 0.0_r8
   endif
   return
endif

obs_sd = sqrt(obs_var)

! If the truth were at point x, what is the weight of the truncated normal obs error dist?
! If no bounds, the whole cdf is possible
cdf(1) = 0.0_r8
cdf(2) = 1.0_r8

! Compute the cdf's at the bounds if they exist
if(bounded_below) cdf(1) = normal_cdf(lower_bound, x, obs_sd)
if(bounded_above) cdf(2) = normal_cdf(upper_bound, x, obs_sd)

! The weight is the reciprocal of the fraction of the cdf that is in legal range
weight = 1.0_r8 / (cdf(2) - cdf(1))

get_truncated_normal_like = weight * exp(-1.0_r8 * (x - obs)**2 / (2.0_r8 * obs_var))

end function get_truncated_normal_like



subroutine obs_increment_particle(ens, ens_size, obs, obs_var, obs_inc)
!------------------------------------------------------------------------
!
! A observation space only particle filter implementation for a
! two step sequential update filter. Second version, 2 October, 2003.

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: weight(ens_size), rel_weight(ens_size), cum_weight(0:ens_size)
real(r8) :: base, frac, new_val(ens_size), weight_sum
integer  :: i, j, indx(ens_size)

! Begin by computing a weight for each of the prior ensemble members
do i = 1, ens_size
   weight(i) = exp(-1.0_r8 * (ens(i) - obs)**2 / (2.0_r8 * obs_var))
end do

! Compute relative weight for each ensemble member
weight_sum = sum(weight)
do i = 1, ens_size
   rel_weight(i) = weight(i) / weight_sum
end do

! Compute cumulative weights at boundaries
cum_weight(0) = 0.0_r8
do i = 1, ens_size
   cum_weight(i) = cum_weight(i - 1) + rel_weight(i)
!   write(*,'(1x,i3,3(e10.4,1x))') i, weight(i), rel_weight(i), cum_weight(i)
end do
! Fix up for round-off error if any
cum_weight(ens_size) = 1.0_r8

! Do a deterministic implementation: just divide interval into ens_size parts and see
! which interval this is in (careful to offset; not start at 0)
base = 1.0_r8 / (ens_size * 2.0_r8)

do i = 1, ens_size

   frac = base + (i - 1.0_r8) / ens_size

   ! Now search in the cumulative range to see where this frac falls
   ! Can make this search more efficient by limiting base
   do j = 1, ens_size
      if(cum_weight(j - 1) < frac .and. frac < cum_weight(j)) then
         indx(i) = j
!         write(*, *) i, frac, 'gets index ', j
         goto 111
      end if
   end do

111 continue

end do

! Set the new values for the ensemble members
do i = 1, ens_size
   new_val(i) = ens(indx(i))
!   write(*, *) 'new_val ', i, new_val(i)
end do

! Generate increments
obs_inc = new_val - ens

end subroutine obs_increment_particle



subroutine obs_increment_enkf(ens, ens_size, prior_var, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment_enkf(ens, ens_size, obs, obs_var, obs_inc)
!

! ENKF version of obs increment

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), prior_var, obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: obs_var_inv, prior_var_inv, new_var, new_mean(ens_size)
! real(r8) :: sx, s_x2
real(r8) :: temp_mean, temp_obs(ens_size)
real(r8) :: new_val(ens_size)
integer  :: i, ens_index(ens_size), new_index(ens_size)

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0_r8 / obs_var
prior_var_inv = 1.0_r8 / prior_var

new_var       = 1.0_r8 / (prior_var_inv + obs_var_inv)

! If this is first time through, need to initialize the random sequence.
! This will reproduce exactly for multiple runs with the same task count,
! but WILL NOT reproduce for a different number of MPI tasks.
! To make it independent of the number of MPI tasks, it would need to
! use the global ensemble number or something else that remains constant
! as the processor count changes.  this is not currently an argument to
! this function and so we are not trying to make it task-count invariant.
if(first_inc_ran_call) then
   call init_random_seq(inc_ran_seq, my_task_id() + 1)
   first_inc_ran_call = .false.
endif

! Generate perturbed obs
do i = 1, ens_size
    temp_obs(i) = random_gaussian(inc_ran_seq, obs, sqrt(obs_var))
end do

! Move this so that it has original obs mean
temp_mean = sum(temp_obs) / ens_size
temp_obs(:) = temp_obs(:) - temp_mean + obs

! Loop through pairs of priors and obs and compute new mean
do i = 1, ens_size
   new_mean(i) = new_var * (prior_var_inv * ens(i) + temp_obs(i) / obs_var)
   obs_inc(i)  = new_mean(i) - ens(i)
end do

! To minimize regression errors, may want to sort to minimize increments
! This makes sense for any of the non-deterministic algorithms
! By doing it here, can take care of both standard non-deterministic updates
! plus non-deterministic obs space covariance inflation. This is expensive, so
! don't use it if it's not needed.
if (sort_obs_inc) then 
   new_val = ens + obs_inc
   ! Sorting to make increments as small as possible
   call index_sort(ens, ens_index, ens_size)
   call index_sort(new_val, new_index, ens_size)
   do i = 1, ens_size
      obs_inc(ens_index(i)) = new_val(new_index(i)) - ens(ens_index(i))
   end do
endif

! Can also adjust mean (and) variance of final sample; works fine
!sx         = sum(new_mean)
!s_x2       = sum(new_mean * new_mean)
!temp_mean = sx / ens_size
!temp_var  = (s_x2 - sx**2 / ens_size) / (ens_size - 1)
!new_mean = (new_mean - temp_mean) * sqrt(new_var / temp_var) + updated_mean
!obs_inc = new_mean - ens


end subroutine obs_increment_enkf



subroutine obs_increment_kernel(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment_kernel(ens, ens_size, obs, obs_var, obs_inc)
!

! Kernel version of obs increment

integer, intent(in)             :: ens_size
real(r8), intent(in)            :: ens(ens_size), obs, obs_var
real(r8), intent(out)           :: obs_inc(ens_size)

real(r8) :: obs_var_inv
real(r8) :: prior_mean, prior_cov_inv, new_cov, prior_cov
real(r8) :: sx
real(r8) :: weight(ens_size), new_mean(ens_size)
real(r8) :: cum_weight, total_weight, cum_frac(ens_size)
real(r8) :: unif, norm, new_member(ens_size)

integer :: i, j, kernel

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0_r8 / obs_var

! Compute prior mean and covariance
sx         = sum(ens)
prior_mean = sx / ens_size
prior_cov  = sum((ens - prior_mean)**2) / (ens_size - 1)

prior_cov     = prior_cov / 10.0_r8     ! For kernels, scale the prior covariance
prior_cov_inv = 1.0_r8 / prior_cov

! Compute new covariance once for these kernels
new_cov = 1.0_r8 / (prior_cov_inv + obs_var_inv)

! New mean is computed ens_size times as is weight
do i = 1, ens_size
   new_mean(i) = new_cov*(prior_cov_inv * ens(i) + obs / obs_var)
   weight(i) =  2.71828_r8 ** (-0.5_r8 * (ens(i)**2 * prior_cov_inv + &
      obs**2 * obs_var_inv - new_mean(i)**2 / new_cov))
end do

! Compute total weight
total_weight = sum(weight)
cum_weight   = 0.0_r8
do i = 1, ens_size
   cum_weight  = cum_weight + weight(i)
   cum_frac(i) = cum_weight / total_weight
end do

! If this is first time through, need to initialize the random sequence.
! This will reproduce exactly for multiple runs with the same task count,
! but WILL NOT reproduce for a different number of MPI tasks.
! To make it independent of the number of MPI tasks, it would need to
! use the global ensemble number or something else that remains constant
! as the processor count changes.  this is not currently an argument to
! this function and so we are not trying to make it task-count invariant.
if(first_inc_ran_call) then
   call init_random_seq(inc_ran_seq, my_task_id() + 1)
   first_inc_ran_call = .false.
endif

! Generate a uniform random number and a Gaussian for each new member
do i = 1, ens_size
   unif = random_uniform(inc_ran_seq)
   ! Figure out which kernel it's in
   whichk: do j = 1, ens_size
      if(unif < cum_frac(j)) then
         kernel = j
         exit whichk
      end if
   end do whichk

   ! Next calculate a unit normal in this kernel
   norm = random_gaussian(inc_ran_seq, 0.0_r8, sqrt(new_cov))
   ! Now generate the new ensemble member
   new_member(i) = new_mean(kernel) + norm
end do

! Generate the increments
obs_inc = new_member - ens

end subroutine obs_increment_kernel



subroutine update_from_obs_inc(obs, obs_prior_mean, obs_prior_var, obs_inc, &
               state, ens_size, state_inc, reg_coef, net_a_in, correl_out)
!========================================================================

! Does linear regression of a state variable onto an observation and
! computes state variable increments from observation increments

integer,            intent(in)    :: ens_size
real(r8),           intent(in)    :: obs(ens_size), obs_inc(ens_size)
real(r8),           intent(in)    :: obs_prior_mean, obs_prior_var
real(r8),           intent(in)    :: state(ens_size)
real(r8),           intent(out)   :: state_inc(ens_size), reg_coef
real(r8),           intent(in) :: net_a_in
real(r8), optional, intent(inout) :: correl_out

real(r8) :: obs_state_cov, intermed
real(r8) :: restoration_inc(ens_size), state_mean, state_var, correl
real(r8) :: factor, exp_true_correl, mean_factor, net_a


! For efficiency, just compute regression coefficient here unless correl is needed

state_mean = sum(state) / ens_size
obs_state_cov = sum( (state - state_mean) * (obs - obs_prior_mean) ) / (ens_size - 1)

if (obs_prior_var > 0.0_r8) then
   reg_coef = obs_state_cov/obs_prior_var
else
   reg_coef = 0.0_r8
endif

! If correl_out is present, need correl for adaptive inflation
! Also needed for file correction below.

! WARNING: we have had several different numerical problems in this
! section, especially with users running in single precision floating point.
! Be very cautious if changing any code in this section, taking into
! account underflow and overflow for 32 bit floats.

if(present(correl_out) .or. sampling_error_correction) then
   if (obs_state_cov == 0.0_r8 .or. obs_prior_var <= 0.0_r8) then
      correl = 0.0_r8
   else
      state_var = sum((state - state_mean)**2) / (ens_size - 1)
      if (state_var <= 0.0_r8) then
         correl = 0.0_r8
      else
         intermed = sqrt(obs_prior_var) * sqrt(state_var)
         if (intermed <= 0.0_r8) then
            correl = 0.0_r8
         else
            correl = obs_state_cov / intermed
         endif
      endif
   endif
   if(correl >  1.0_r8) correl =  1.0_r8
   if(correl < -1.0_r8) correl = -1.0_r8
endif
if(present(correl_out)) correl_out = correl


! Get the expected actual correlation and the regression weight reduction factor
if(sampling_error_correction) then
   call get_correction_from_table(correl, mean_factor, exp_true_correl, ens_size)
   ! Watch out for division by zero; if correl is really small regression is safely 0
   if(abs(correl) > 0.001_r8) then
      reg_coef = reg_coef * (exp_true_correl / correl) * mean_factor
   else
      reg_coef = 0.0_r8
   endif
   correl = exp_true_correl
endif



! Then compute the increment as product of reg_coef and observation space increment
state_inc = reg_coef * obs_inc

!! NOTE: if requested to be returned, correl_out is set further up in the
!! code, before the sampling error correction, if enabled, is applied.
!! this means it's returning a different larger value than the correl
!! being returned here.  it's used by the adaptive inflation and so the
!! inflation will see a slightly different correlation value.  it isn't
!! clear that this is a bad thing; it means the inflation might be a bit
!! larger than it would otherwise.  before we move any code this would
!! need to be studied to see what the real impact would be.

end subroutine update_from_obs_inc


!------------------------------------------------------------------------

subroutine get_correction_from_table(scorrel, mean_factor, expected_true_correl, ens_size)

real(r8),  intent(in) :: scorrel
real(r8), intent(out) :: mean_factor, expected_true_correl
integer,  intent(in)  :: ens_size

! Uses interpolation to get correction factor into the table

integer             :: low_indx, high_indx
real(r8)            :: correl, fract, low_correl, low_exp_correl, low_alpha
real(r8)            :: high_correl, high_exp_correl, high_alpha

logical, save :: first_time = .true.

if (first_time) then
   call read_sampling_error_correction(ens_size, exp_true_correl, alpha)
   first_time = .false.
endif

! Interpolate to get values of expected correlation and mean_factor
if(scorrel < -1.0_r8) then
   correl = -1.0_r8
   mean_factor = 1.0_r8
else if(scorrel > 1.0_r8) then
   correl = 1.0_r8
   mean_factor = 1.0_r8
else if(scorrel <= -0.995_r8) then
   fract = (scorrel + 1.0_r8) / 0.005_r8
   correl = (exp_true_correl(1) + 1.0_r8) * fract - 1.0_r8
   mean_factor = (alpha(1) - 1.0_r8) * fract + 1.0_r8
else if(scorrel >= 0.995_r8) then
   fract = (scorrel - 0.995_r8) / 0.005_r8
   correl = (1.0_r8 - exp_true_correl(sec_table_size)) * fract + exp_true_correl(sec_table_size)
   mean_factor = (1.0_r8 - alpha(sec_table_size)) * fract + alpha(sec_table_size)
else
   ! given the ifs above, the floor() computation below for low_indx
   ! should always result in a value in the range 1 to 199.  but if this
   ! code is compiled with r8=r4 (single precision reals) it turns out
   ! to be possible to get values a few bits below 0 which results in
   ! a very large negative integer.  the limit tests below ensure the
   ! index stays in a legal range.
   low_indx = floor((scorrel + 0.995_r8) / 0.01_r8 + 1.0_r8)
   if (low_indx <   1) low_indx =   1
   if (low_indx > 199) low_indx = 199
   low_correl = -0.995_r8 + (low_indx - 1) * 0.01_r8
   low_exp_correl = exp_true_correl(low_indx)
   low_alpha = alpha(low_indx)
   high_indx = low_indx + 1
   high_correl = low_correl + 0.01_r8
   high_exp_correl = exp_true_correl(high_indx)
   high_alpha = alpha(high_indx)
   fract = (scorrel - low_correl) / (high_correl - low_correl)
   correl = (high_exp_correl - low_exp_correl) * fract + low_exp_correl
   mean_factor = (high_alpha - low_alpha) * fract + low_alpha
endif

expected_true_correl = correl

! Don't want Monte Carlo interpolation problems to put us outside of a
! ratio between 0 and 1 for expected_true_correl / sample_correl
! If they have different signs, expected should just be 0
if(expected_true_correl * scorrel <= 0.0_r8) then
   expected_true_correl = 0.0_r8
else if(abs(expected_true_correl) > abs(scorrel)) then
   ! If same sign, expected should not be bigger in absolute value
   expected_true_correl = scorrel
endif

end subroutine get_correction_from_table


subroutine obs_increment_rank_histogram(ens, ens_size, prior_var, &
   obs, obs_var, obs_inc)
!------------------------------------------------------------------------
!
! Does observation space update by approximating the prior distribution by
! a rank histogram. Prior and posterior are assumed to have 1/(n+1) probability
! mass between each ensemble member. The tails are assumed to be gaussian with
! a variance equal to sample variance of the entire ensemble and a mean
! selected so that 1/(n+1) of the mass is in each tail.
!
! The likelihood between the extreme ensemble members is approximated by
! quadrature. Two options are available and controlled by the namelist entry
! rectangular_quadrature. If this namelist is true than the likelihood between
! a pair of ensemble members is assumed to be uniform with the average of
! the likelihood computed at the two ensemble members. If it is false then
! the likelihood between two ensemble members is approximated by a line
! connecting the values of the likelihood computed at each of the ensemble
! members (trapezoidal quadrature).
!
! Two options are available for approximating the likelihood on the tails.
! If gaussian_likelihood_tails is true that the likelihood is assumed to
! be N(obs, obs_var) on the tails. If this is false, then the likelihood
! on the tails is taken to be uniform (to infinity) with the value at the
! outermost ensemble members.
!
! A product of the approximate prior and approximate posterior is taken
! and new ensemble members are located so that 1/(n+1) of the mass is between
! each member and on the tails.

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), prior_var, obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

integer  :: i, e_ind(ens_size), lowest_box, j
real(r8) :: prior_sd, var_ratio, umass, left_amp, right_amp, norm_const
real(r8) :: left_mean, right_mean
real(r8) :: mass(ens_size + 1), like(ens_size), cumul_mass(0:ens_size + 1)
real(r8) :: nmass(ens_size + 1)
real(r8) :: new_mean_left, new_mean_right, prod_weight_left, prod_weight_right
real(r8) :: new_var_left, new_var_right, new_sd_left, new_sd_right
real(r8) :: new_ens(ens_size), mass_sum
real(r8) :: x(ens_size)
real(r8) :: like_dense(2:ens_size), height(2:ens_size)
real(r8) :: dist_for_unit_sd
real(r8) :: a, b, c, hright, hleft, r1, r2, adj_r1, adj_r2

! Do an index sort of the ensemble members; Will want to do this very efficiently
call index_sort(ens, e_ind, ens_size)
x = ens(e_ind)

! Define normal PDF constant term
norm_const = 1.0_r8 / sqrt(2.0_r8 * PI * obs_var)
! Compute likelihood for each ensemble member; just evaluate the gaussian
do i = 1, ens_size
   like(i) = norm_const * exp(-1.0_r8 * (x(i) - obs)**2 / (2.0_r8 * obs_var))
end do

! Compute approx likelihood each interior bin (average of bounding likelihoods)
do i = 2, ens_size
   like_dense(i) = ((like(i - 1) + like(i)) / 2.0_r8)
end do

! Compute the s.d. of the ensemble for getting the gaussian tails
prior_sd = sqrt(prior_var)

! For unit normal, find distance from mean to where cdf is 1/(n+1)
! Lots of this can be done once in first call and then saved
dist_for_unit_sd = inv_weighted_normal_cdf(1.0_r8, 0.0_r8, 1.0_r8, &
   1.0_r8 / (ens_size + 1.0_r8))
dist_for_unit_sd = -1.0_r8 * dist_for_unit_sd

! Have variance of tails just be sample prior variance
! Mean is adjusted so that 1/(n+1) is outside
left_mean = x(1) + dist_for_unit_sd * prior_sd
! Same for right tail
right_mean = x(ens_size) - dist_for_unit_sd * prior_sd

if(gaussian_likelihood_tails) then
   !*************** Block to do Gaussian-Gaussian on tail **************
   ! Compute the product of the obs likelihood gaussian with the priors
   ! Left tail gaussian first
   var_ratio = obs_var / (prior_var + obs_var)
   new_var_left = var_ratio * prior_var
   new_sd_left = sqrt(new_var_left)
   new_mean_left  = var_ratio * (left_mean  + prior_var*obs / obs_var)
   ! REMEMBER, this product has an associated weight which must be taken into account
   ! See Anderson and Anderson for this weight term (or tutorial kernel filter)
   prod_weight_left =  exp(-0.5_r8 * (left_mean**2 / prior_var + &
         obs**2 / obs_var - new_mean_left**2 / new_var_left)) / &
         sqrt(prior_var + obs_var) / sqrt(2.0_r8 * PI)
   ! Determine how much mass is in the updated tails by computing gaussian cdf
   mass(1) = normal_cdf(x(1), new_mean_left, new_sd_left) * prod_weight_left

   ! Same for the right tail
   var_ratio = obs_var / (prior_var + obs_var)
   new_var_right = var_ratio * prior_var
   new_sd_right = sqrt(new_var_right)
   new_mean_right  = var_ratio * (right_mean  + prior_var*obs / obs_var)
   prod_weight_right =  exp(-0.5_r8 * (right_mean**2 / prior_var + &
         obs**2 / obs_var - new_mean_right**2 / new_var_right)) / &
         sqrt(prior_var + obs_var) / sqrt(2.0_r8 * PI)
   ! Determine how much mass is in the updated tails by computing gaussian cdf
   mass(ens_size + 1) = (1.0_r8 - normal_cdf(x(ens_size), new_mean_right, new_sd_right)) * &
      prod_weight_right

   !************ End Block to do Gaussian-Gaussian on tail **************
else
   !*************** Block to do flat tail for likelihood ****************
   ! Flat tails: THIS REMOVES ASSUMPTIONS ABOUT LIKELIHOOD AND CUTS COST
   new_sd_left = prior_sd
   new_mean_left = left_mean
   prod_weight_left = like(1)
   mass(1) = like(1) / (ens_size + 1.0_r8)

   ! Same for right tail
   new_sd_right = prior_sd
   new_mean_right = right_mean
   prod_weight_right = like(ens_size)
   mass(ens_size + 1) = like(ens_size) / (ens_size + 1.0_r8)
   !*************** End block to do flat tail for likelihood ****************
endif

! The mass in each interior box is the height times the width
! The height of the likelihood is like_dense
! For the prior, mass is 1/(n+1),   and mass = height x width so...
! The height of the prior is 1 / ((n+1) width);   multiplying by width leaves 1/(n+1)

! In prior, have 1/(n+1) mass in each bin, multiply by mean likelihood density
! to get approximate mass in updated bin
do i = 2, ens_size
   mass(i) = like_dense(i) / (ens_size + 1.0_r8)
   ! Height of prior in this bin is mass/width; Only needed for trapezoidal
   ! If two ensemble members are the same, set height to -1 as flag
   if(x(i) == x(i - 1)) then
      height(i) = -1.0_r8
   else
      height(i) = 1.0_r8 / ((ens_size + 1.0_r8) * (x(i) - x(i-1)))
   endif
end do

! Now normalize the mass in the different bins to get a pdf
mass_sum = sum(mass)
nmass = mass / mass_sum

! Get the weight for the final normalized tail gaussians
! This is the same as left_amp=(ens_size + 1)*nmass(1)
left_amp = prod_weight_left / mass_sum
! This is the same as right_amp=(ens_size + 1)*nmass(ens_size + 1)
right_amp = prod_weight_right / mass_sum

! Find cumulative mass at each box boundary and middle boundary
cumul_mass(0) = 0.0_r8
do i = 1, ens_size + 1
   cumul_mass(i) = cumul_mass(i - 1) + nmass(i)
end do

! Begin intenal box search at bottom of lowest box, update for efficiency
lowest_box = 1

! Find each new ensemble members location
do i = 1, ens_size
   ! Each update ensemble member has 1/(n+1) mass before it
   umass = (1.0_r8 * i) / (ens_size + 1.0_r8)

   ! If it is in the inner or outer range have to use normal
   if(umass < cumul_mass(1)) then
      ! It's in the left tail
      ! Get position of x in weighted gaussian where the cdf has value umass
      new_ens(i) = inv_weighted_normal_cdf(left_amp, new_mean_left, new_sd_left, &
         umass)
   else if(umass > cumul_mass(ens_size)) then
      ! It's in the right tail
      ! Get position of x in weighted gaussian where the cdf has value umass
      new_ens(i) = inv_weighted_normal_cdf(right_amp, new_mean_right, new_sd_right, &
         1.0_r8 - umass)
      ! Coming in from the right, use symmetry after pretending its on left
      new_ens(i) = new_mean_right + (new_mean_right - new_ens(i))
   else
      ! In one of the inner uniform boxes.
      FIND_BOX:do j = lowest_box, ens_size - 1
         ! Find the box that this mass is in
         if(umass >= cumul_mass(j) .and. umass <= cumul_mass(j + 1)) then

            if(rectangular_quadrature) then
               !********* Block for rectangular quadrature *******************
               ! Linearly interpolate in mass
               new_ens(i) = x(j) + ((umass - cumul_mass(j)) / &
                  (cumul_mass(j+1) - cumul_mass(j))) * (x(j + 1) - x(j))
               !********* End block for rectangular quadrature *******************

            else

               !********* Block for trapezoidal interpolation *******************
               ! Assume that mass has linear profile, quadratic interpolation
               ! If two ensemble members are the same, just keep that value
               if(height(j + 1) < 0) then
                  new_ens(i) = x(j)
               else
                  ! Height on left side and right side
                  hleft = height(j + 1) * like(j) / mass_sum
                  hright = height(j + 1) * like(j + 1) / mass_sum
                  ! Will solve a quadratic for desired x-x(j)
                  ! a is 0.5(hright - hleft) / (x(j+1) - x(j))
                  a = 0.5_r8 * (hright - hleft) / (x(j+1) - x(j))
                  ! b is hleft
                  b = hleft
                  ! c is cumul_mass(j) - umass
                  c = cumul_mass(j) - umass
                  ! Use stable quadratic solver
                  call solve_quadratic(a, b, c, r1, r2)
                  adj_r1 = r1 + x(j)
                  adj_r2 = r2 + x(j)
                  if(adj_r1 >= x(j) .and. adj_r1 <= x(j+1)) then
                     new_ens(i) = adj_r1
                  elseif (adj_r2 >= x(j) .and. adj_r2 <= x(j+1)) then
                     new_ens(i) = adj_r2
                  else
                     msgstring = 'Did not get a satisfactory quadratic root'
                     call error_handler(E_ERR, 'obs_increment_rank_histogram', msgstring, &
                        source)
                  endif
               endif
               !********* End block for quadratic interpolation *******************

            endif

            ! Don't need to search lower boxes again
            lowest_box = j
            exit FIND_BOX
         end if
      end do FIND_BOX
   endif
end do

! Convert to increments for unsorted
do i = 1, ens_size
   obs_inc(e_ind(i)) = new_ens(i) - x(i)
end do

end subroutine obs_increment_rank_histogram




subroutine update_ens_from_weights(ens, ens_size, rel_weight, ens_inc)
!------------------------------------------------------------------------
! Given relative weights for an ensemble, compute increments for the
! ensemble members. Assumes that prior distributon is equal uniform mass
! between each ensemble member. On the edges, have a normal with the
! sample mean and s.d. BUT normalized by a factor alpha so that only
! 1/(2*ens_size) of the total mass lies on each flank.

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), rel_weight(ens_size)
real(r8), intent(out) :: ens_inc(ens_size)

integer  :: i, j, lowest_box
integer  :: e_ind(ens_size)
real(r8) :: x(1:2*ens_size - 1), cumul_mass(1:2*ens_size - 1), new_ens(ens_size)
real(r8) :: sort_inc(ens_size), updated_mass(2 * ens_size)
real(r8) :: sx, prior_mean, prior_var, prior_sd, mass
real(r8) :: total_mass_left, total_mass_right, alpha(2)

! Initialize assim_tools_module if needed
if (.not. module_initialized) call assim_tools_init()

call error_handler(E_ERR,'update_ens_from_weight','Routine needs testing.', &
           source, text2='Talk to Jeff before using.')

! Do an index sort of the ensemble members
call index_sort(ens, e_ind, ens_size)

! Have half boxes between all ensembles in the interior
! Total number of mass boxes is 2*ens_size

! Compute the points that bound all the updated mass boxes; start with ensemble
do i = 1, ens_size
   x(2*i - 1) = ens(e_ind(i))
end do
! Compute the mid-point interior boundaries; these are halfway between ensembles
do i = 2, 2*ens_size - 2, 2
   x(i) = (x(i - 1) + x(i + 1)) / 2.0_r8
end do

! Compute the mean and s.d. of the prior ensemble to handle wings
sx         = sum(ens)
prior_mean = sx / ens_size
prior_var  = sum((ens - prior_mean)**2) / (ens_size - 1)
prior_sd = sqrt(prior_var)

! Need to normalize the wings so they have 1/(2*ens_size) mass outside
! Use cdf to find out how much mass is left of 1st member, right of last
total_mass_left = normal_cdf(ens(e_ind(1)), prior_mean, prior_sd)
total_mass_right = 1.0_r8 - normal_cdf(ens(e_ind(ens_size)), prior_mean, prior_sd)

! Find the mass in each division given the initial equal partition and the weights
updated_mass(1) = rel_weight(e_ind(1)) / (2.0_r8 * ens_size)
updated_mass(2 * ens_size) = rel_weight(e_ind(ens_size)) / (2.0_r8 * ens_size)
do i = 2, 2*ens_size - 2, 2
   updated_mass(i) = rel_weight(e_ind(i / 2)) / (2.0_r8 * ens_size)
end do
do i = 3, 2*ens_size - 1, 2
   updated_mass(i) = rel_weight(e_ind((i+1) / 2)) / (2.0_r8 * ens_size)
end do

! Normalize the mass; (COULD IT EVER BE 0 necessitating error check?)
updated_mass = updated_mass / sum(updated_mass)

! Find a normalization factor to get tail mass right
if(total_mass_left > 0.0_r8) then
   alpha(1) = updated_mass(1) / total_mass_left
else
   alpha(1) = 0.0_r8
endif
if(total_mass_right > 0.0_r8) then
   alpha(2) = updated_mass(2 * ens_size) / total_mass_right
else
   alpha(2) = 0.0_r8
endif

! Find cumulative mass at each box boundary and middle boundary
cumul_mass(1) = updated_mass(1)
do i = 2, 2*ens_size - 1
   cumul_mass(i) = cumul_mass(i - 1) + updated_mass(i)
end do

! Get resampled position an inefficient way
! Need 1/ens_size between each EXCEPT for outers which get half of this
mass = 1.0_r8 / (2.0_r8 * ens_size)

do i = 1, ens_size
   ! If it's in the inner or outer range have to use normal
   if(mass < cumul_mass(1)) then
      ! In the first normal box
      new_ens(i) = inv_weighted_normal_cdf(alpha(1), prior_mean, prior_sd, mass)
   else if(mass > cumul_mass(2*ens_size - 1)) then
      ! In the last normal box; Come in from the outside
      new_ens(i) = inv_weighted_normal_cdf(alpha(2), prior_mean, prior_sd, 1.0_r8 - mass)
      new_ens(i) = prior_mean + (prior_mean - new_ens(i))
   else
      ! In one of the inner uniform boxes. Make this much more efficient search?
      lowest_box = 1
      FIND_BOX:do j = lowest_box, 2 * ens_size - 2
         ! Find the box that this mass is in
         if(mass >= cumul_mass(j) .and. mass <= cumul_mass(j + 1)) then
            new_ens(i) = x(j) + ((mass - cumul_mass(j)) / (cumul_mass(j+1) - cumul_mass(j))) * &
               (x(j + 1) - x(j))
            ! Don't need to search lower boxes again
            lowest_box = j
            exit FIND_BOX
         end if
      end do FIND_BOX
   endif
   ! Want equally partitioned mass in update with exception that outermost boxes have half
   mass = mass + 1.0_r8 / ens_size
end do

! Can now compute sorted increments
do i = 1, ens_size
   sort_inc(i) = new_ens(i) - ens(e_ind(i))
end do

! Now, need to convert to increments for unsorted
do i = 1, ens_size
   ens_inc(e_ind(i)) = sort_inc(i)
end do

end subroutine update_ens_from_weights
!---------------------------------------------------------------

subroutine obs_updates_ens(ens_size, num_groups, ens, ens_loc, ens_kind, &
   obs_prior, obs_inc, obs_prior_mean, obs_prior_var, obs_loc, obs_type, obs_time,    &
   net_a, grp_size, grp_beg, grp_end, reg_factor_obs_index,         &
   reg_factor_ens_index, final_factor, correl, correl_needed, inflate_only)

integer,             intent(in)  :: ens_size
integer,             intent(in)  :: num_groups
real(r8),            intent(inout)  :: ens(ens_size)
type(location_type), intent(in)  :: ens_loc
integer,             intent(in)  :: ens_kind
real(r8),            intent(in)  :: obs_prior(ens_size)
real(r8),            intent(in)  :: obs_inc(ens_size)
real(r8),            intent(in)  :: obs_prior_mean(num_groups)
real(r8),            intent(in)  :: obs_prior_var(num_groups)
type(location_type), intent(in)  :: obs_loc
integer,             intent(in)  :: obs_type
type(time_type),     intent(in)  :: obs_time
real(r8),            intent(in)  :: net_a(num_groups)
integer,             intent(in)  :: grp_size
integer,             intent(in)  :: grp_beg(num_groups)
integer,             intent(in)  :: grp_end(num_groups)
integer,             intent(in)  :: reg_factor_obs_index
integer(i8),         intent(in)  :: reg_factor_ens_index
real(r8),            intent(inout) :: final_factor
real(r8),            intent(out) :: correl(num_groups)
logical,             intent(in)  :: correl_needed
logical,             intent(in)  :: inflate_only

real(r8) :: reg_coef(num_groups), increment(ens_size)
real(r8) :: reg_factor
integer  :: group, grp_bot, grp_top

! Loop through groups to update the state variable ensemble members
do group = 1, num_groups
   grp_bot = grp_beg(group); grp_top = grp_end(group)
   ! Do update of state, correl only needed for varying ss inflate
   if(correl_needed) then
      call update_from_obs_inc(obs_prior(grp_bot:grp_top), obs_prior_mean(group), &
         obs_prior_var(group), obs_inc(grp_bot:grp_top), ens(grp_bot:grp_top), grp_size, &
         increment(grp_bot:grp_top), reg_coef(group), net_a(group), correl(group))
   else
      call update_from_obs_inc(obs_prior(grp_bot:grp_top), obs_prior_mean(group), &
         obs_prior_var(group), obs_inc(grp_bot:grp_top), ens(grp_bot:grp_top), grp_size, &
         increment(grp_bot:grp_top), reg_coef(group), net_a(group))
   endif
end do

if(num_groups > 1) then
   reg_factor = comp_reg_factor(num_groups, reg_coef, obs_time, &
      reg_factor_obs_index, reg_factor_ens_index)
   final_factor = min(final_factor, reg_factor)
endif

! Get the updated ensemble
if(.not. inflate_only) ens = ens + final_factor * increment

end subroutine obs_updates_ens

!-------------------------------------------------------------

function cov_and_impact_factors(base_obs_loc, base_obs_type, state_loc, state_kind, &
dist, cutoff_rev)

! Computes the cov_factor and multiplies by obs_impact_factor if selected

real(r8) :: cov_and_impact_factors
type(location_type), intent(in) :: base_obs_loc
integer, intent(in) :: base_obs_type
type(location_type), intent(in) :: state_loc
integer, intent(in) :: state_kind
real(r8), intent(in) :: dist
real(r8), intent(in) :: cutoff_rev

real(r8) :: impact_factor, cov_factor

! Get external impact factors, cycle if impact of this ob on this state is zero
if (adjust_obs_impact) then
   ! Get the impact factor from the table if requested
   impact_factor = obs_impact_table(base_obs_type, state_kind)
   if(impact_factor <= 0.0_r8) then
      ! Avoid the cost of computing cov_factor if impact is 0
      cov_and_impact_factors = 0.0_r8
      return
   endif
else
   impact_factor = 1.0_r8
endif

! Compute the covariance factor
cov_factor = comp_cov_factor(dist, cutoff_rev, &
   base_obs_loc, base_obs_type, state_loc, state_kind)

! Combine the impact_factor and the cov_factor
cov_and_impact_factors = cov_factor * impact_factor

end function cov_and_impact_factors


!------------------------------------------------------------------------

subroutine set_assim_tools_trace(execution_level, timestamp_level)
 integer, intent(in) :: execution_level
 integer, intent(in) :: timestamp_level

! set module local vars from the calling code to indicate how much
! output we should generate from this code.  execution level is
! intended to make it easier to figure out where in the code a crash
! is happening; timestamp level is intended to help with gross levels
! of overall performance profiling.  eventually, a level of 1 will
! print out only basic info; level 2 will be more detailed.
! (right now, only > 0 prints anything and it doesn't matter how
! large the value is.)

! Initialize assim_tools_module if needed
if (.not. module_initialized) call assim_tools_init()

print_trace_details = execution_level
print_timestamps    = timestamp_level

end subroutine set_assim_tools_trace

!--------------------------------------------------------------------

function revised_distance(orig_dist, newcount, oldcount, base, cutfloor)
 real(r8),            intent(in) :: orig_dist
 integer,             intent(in) :: newcount, oldcount
 type(location_type), intent(in) :: base
 real(r8),            intent(in) :: cutfloor

 real(r8)                        :: revised_distance

! take the ratio of the old and new counts, and revise the
! original cutoff distance to match.

! for now, only allow the code to do a 2d area adaption.
! to experiment with other schemes, set this local variable
! to .false. at the top of the file and recompile.

if (only_area_adapt) then

   revised_distance = orig_dist * sqrt(real(newcount, r8) / oldcount)

   ! allow user to set a minimum cutoff, so even if there are very dense
   ! observations the cutoff distance won't go below this floor.
   if (revised_distance < cutfloor) revised_distance = cutfloor
   return

endif

! alternatives for different dimensionalities and schemes

! Change the cutoff radius to get the appropriate number
if (LocationDims == 1) then
   ! linear (be careful of cyclic domains; if > domain, this is
   ! not going to be right)
   revised_distance = orig_dist * real(newcount, r8) / oldcount

else if (LocationDims == 2) then
   ! do an area scaling
   revised_distance = orig_dist * sqrt(real(newcount, r8) / oldcount)

else if (LocationDims == 3) then
   ! do either a volume or area scaling (depending on whether we are
   ! localizing in the vertical or not.)   if surface obs, assume a hemisphere
   ! and shrink more.

   if (vertical_localization_on()) then
      ! cube root for volume
      revised_distance = orig_dist * ((real(newcount, r8) / oldcount) &
                                      ** 0.33333333333333333333_r8)

      ! Cut the adaptive localization threshold in half again for 'surface' obs
      if (is_vertical(base, "SURFACE")) then
         revised_distance = revised_distance * (0.5_r8 ** 0.33333333333333333333_r8)
      endif
   else
      ! do an area scaling, even if 3d obs
      revised_distance = orig_dist * sqrt(real(newcount, r8) / oldcount)

      ! original code was:
      !cutoff_rev =  sqrt((2.0_r8*cutoff)**2 * adaptive_localization_threshold / &
      !   total_num_close_obs) / 2.0_r8

      ! original comment
      ! Need to get thinning out of assim_tools and into something about locations
   endif
else
   call error_handler(E_ERR, 'revised_distance', 'unknown locations dimension, not 1, 2 or 3', &
      source)
endif

! allow user to set a minimum cutoff, so even if there are very dense
! observations the cutoff distance won't go below this floor.
if (revised_distance < cutfloor) revised_distance = cutfloor

end function revised_distance

!--------------------------------------------------------------------

function count_close(num_close, index_list, my_types, dist, maxdist)
 integer, intent(in)  :: num_close, index_list(:), my_types(:)
 real(r8), intent(in) :: dist(:), maxdist
 integer :: count_close

! return the total number of items from the index_list which
! are types which are going to be assimilated, and within distance.
! this excludes items on the eval list only, not listed, or
! items too far away.   this routine does a global communication
! so if any MPI tasks make this call, all must.

integer :: k, thistype, local_count

local_count = 0
do k=1, num_close

   ! only accept items closer than limit
   if (dist(k) > maxdist) cycle

   ! include identity obs, plus types on assim list.
   ! you have to do the if tests separately because fortran allows
   ! both parts of an if(a .or. b) test to be eval'd at the same time.
   ! you'd be using a negative index if it was an identity obs.
   thistype = my_types(index_list(k))
   if (thistype < 0) then
      local_count = local_count + 1
   else if (assimilate_this_type_of_obs(thistype)) then
      local_count = local_count + 1
   endif
end do

! broadcast sums from all tasks to compute new total
call sum_across_tasks(local_count, count_close)

end function count_close

!----------------------------------------------------------------------
! Revise the cutoff for this observation if adaptive localization is required
! Output diagnostics for localization if requested

subroutine adaptive_localization_and_diags(cutoff_orig, cutoff_rev, adaptive_localization_threshold, &
   adaptive_cutoff_floor, num_close_obs, close_obs_ind, close_obs_dist, my_obs_type, &
   base_obs_index, base_obs_loc, obs_def, out_unit)

real(r8),            intent(in)  :: cutoff_orig
real(r8),            intent(out) :: cutoff_rev
integer,             intent(in)  :: adaptive_localization_threshold
real(r8),            intent(in)  :: adaptive_cutoff_floor
integer,             intent(in)  :: num_close_obs
integer,             intent(in)  :: close_obs_ind(:)
real(r8),            intent(in)  :: close_obs_dist(:)
integer,             intent(in)  :: my_obs_type(:)
integer,             intent(in)  :: base_obs_index
type(location_type), intent(in)  :: base_obs_loc
type(obs_def_type),  intent(in)  :: obs_def
integer,             intent(in)  :: out_unit

integer :: total_num_close_obs, rev_num_close_obs, secs, days
type(time_type) :: this_obs_time
character(len = 200) :: base_loc_text   ! longer than longest location formatting possible

! Default is that cutoff is not revised
cutoff_rev = cutoff_orig

! For adaptive localization, need number of other obs close to the chosen observation
if(adaptive_localization_threshold > 0) then
   ! this does a cross-task sum, so all tasks must make this call.
   total_num_close_obs = count_close(num_close_obs, close_obs_ind, my_obs_type, &
                                     close_obs_dist, cutoff_rev*2.0_r8)

   ! Want expected number of close observations to be reduced to some threshold;
   ! accomplish this by cutting the size of the cutoff distance.
   if(total_num_close_obs > adaptive_localization_threshold) then
      cutoff_rev = revised_distance(cutoff_rev*2.0_r8, adaptive_localization_threshold, &
                                    total_num_close_obs, base_obs_loc, &
                                    adaptive_cutoff_floor*2.0_r8) / 2.0_r8
   endif
endif

if ( output_localization_diagnostics ) then
   ! Warning, this can be costly and generate large output
   ! This is referred to as revised in case adaptive localization was done
   rev_num_close_obs = count_close(num_close_obs, close_obs_ind, my_obs_type, &
                                     close_obs_dist, cutoff_rev*2.0_r8)

   ! Output diagnostic information about the number of close obs
   if (my_task_id() == 0) then
      this_obs_time = get_obs_def_time(obs_def)
      call get_time(this_obs_time,secs,days)
      call write_location(-1, base_obs_loc, charstring=base_loc_text)

      ! If adaptive localization did something, output info about what it did
      ! Probably would be more consistent to just output for all observations
      if(adaptive_localization_threshold > 0 .and. &
         total_num_close_obs > adaptive_localization_threshold) then
         write(out_unit,'(i12,1x,i5,1x,i8,1x,A,2(f14.5,1x,i12))') base_obs_index, &
            secs, days, trim(base_loc_text), cutoff_orig, total_num_close_obs, cutoff_rev, &
            rev_num_close_obs
      else
         write(out_unit,'(i12,1x,i5,1x,i8,1x,A,f14.5,1x,i12)') base_obs_index, &
            secs, days, trim(base_loc_text), cutoff_rev, rev_num_close_obs
      endif
   endif
endif

end subroutine adaptive_localization_and_diags

!----------------------------------------------------------------------
!> gets the location of of all my observations
subroutine get_my_obs_loc(obs_ens_handle, obs_seq, keys, my_obs_loc, my_obs_kind, my_obs_type, my_obs_time)

type(ensemble_type),      intent(in)  :: obs_ens_handle
type(obs_sequence_type),  intent(in)  :: obs_seq
integer,                  intent(in)  :: keys(:)
type(location_type),      intent(out) :: my_obs_loc(:)
integer,                  intent(out) :: my_obs_type(:), my_obs_kind(:)
type(time_type),          intent(out) :: my_obs_time

type(obs_type) :: observation
type(obs_def_type)   :: obs_def
integer :: this_obs_key
integer i
type(location_type) :: dummyloc

Get_Obs_Locations: do i = 1, obs_ens_handle%my_num_vars

   this_obs_key = keys(obs_ens_handle%my_vars(i)) ! if keys becomes a local array, this will need changing
   call get_obs_from_key(obs_seq, this_obs_key, observation)
   call get_obs_def(observation, obs_def)
   my_obs_loc(i)  = get_obs_def_location(obs_def)
   my_obs_type(i) = get_obs_def_type_of_obs(obs_def)
   if (my_obs_type(i) > 0) then
         my_obs_kind(i) = get_quantity_for_type_of_obs(my_obs_type(i))
   else
      call get_state_meta_data(-1 * int(my_obs_type(i),i8), dummyloc, my_obs_kind(i))
   endif
end do Get_Obs_Locations

! Need the time for regression diagnostics potentially; get from first observation
my_obs_time = get_obs_def_time(obs_def)

end subroutine get_my_obs_loc

!--------------------------------------------------------------------
!> Get close obs from cache if appropriate. Cache new get_close_obs info
!> if requested.

subroutine get_close_obs_cached(gc_obs, base_obs_loc, base_obs_type, &
   my_obs_loc, my_obs_kind, my_obs_type, num_close_obs, close_obs_ind, close_obs_dist,  &
   ens_handle, last_base_obs_loc, last_num_close_obs, num_close_obs_cached,               &
   num_close_obs_calls_made)

type(get_close_type),          intent(in)  :: gc_obs
type(location_type),           intent(inout) :: base_obs_loc, my_obs_loc(:)
integer,                       intent(in)  :: base_obs_type, my_obs_kind(:), my_obs_type(:)
integer,                       intent(out) :: num_close_obs
integer,                       intent(inout) :: close_obs_ind(:)
real(r8),                      intent(inout) :: close_obs_dist(:)
type(ensemble_type),           intent(in)  :: ens_handle
type(location_type), intent(inout) :: last_base_obs_loc
integer, intent(inout) :: last_num_close_obs
integer, intent(inout) :: num_close_obs_cached, num_close_obs_calls_made

! This logic could be arranged to make code less redundant
if (.not. close_obs_caching) then
   call get_close_obs(gc_obs, base_obs_loc, base_obs_type, &
                      my_obs_loc, my_obs_kind, my_obs_type, &
                      num_close_obs, close_obs_ind, close_obs_dist, ens_handle)
else
   if (base_obs_loc == last_base_obs_loc) then
      num_close_obs     = last_num_close_obs
      num_close_obs_cached = num_close_obs_cached + 1
   else
      call get_close_obs(gc_obs, base_obs_loc, base_obs_type, &
                         my_obs_loc, my_obs_kind, my_obs_type, &
                         num_close_obs, close_obs_ind, close_obs_dist, ens_handle)

      last_base_obs_loc      = base_obs_loc
      last_num_close_obs     = num_close_obs
      num_close_obs_calls_made = num_close_obs_calls_made +1
   endif
endif

end subroutine get_close_obs_cached

!--------------------------------------------------------------------
!> Get close state from cache if appropriate. Cache new get_close_state info
!> if requested.

subroutine get_close_state_cached(gc_state, base_obs_loc, base_obs_type, &
   my_state_loc, my_state_kind, my_state_indx, num_close_states, close_state_ind, close_state_dist,  &
   ens_handle, last_base_states_loc, last_num_close_states, num_close_states_cached,               &
   num_close_states_calls_made)

type(get_close_type),          intent(in)    :: gc_state
type(location_type),           intent(inout) :: base_obs_loc, my_state_loc(:)
integer,                       intent(in)    :: base_obs_type, my_state_kind(:)
integer(i8),                   intent(in)    :: my_state_indx(:)
integer,                       intent(out)   :: num_close_states
integer,                       intent(inout) :: close_state_ind(:)
real(r8),                      intent(inout) :: close_state_dist(:)
type(ensemble_type),           intent(in)    :: ens_handle
type(location_type), intent(inout) :: last_base_states_loc
integer, intent(inout) :: last_num_close_states
integer, intent(inout) :: num_close_states_cached, num_close_states_calls_made

! This logic could be arranged to make code less redundant
if (.not. close_obs_caching) then
   call get_close_state(gc_state, base_obs_loc, base_obs_type, &
                      my_state_loc, my_state_kind, my_state_indx, &
                      num_close_states, close_state_ind, close_state_dist, ens_handle)
else
   if (base_obs_loc == last_base_states_loc) then
      num_close_states     = last_num_close_states
      num_close_states_cached = num_close_states_cached + 1
   else
      call get_close_state(gc_state, base_obs_loc, base_obs_type, &
                         my_state_loc, my_state_kind, my_state_indx, &
                         num_close_states, close_state_ind, close_state_dist, ens_handle)

      last_base_states_loc      = base_obs_loc
      last_num_close_states     = num_close_states
      num_close_states_calls_made = num_close_states_calls_made +1
   endif
endif

end subroutine get_close_state_cached

!--------------------------------------------------------------------
!> log what the user has selected via the namelist choices

subroutine log_namelist_selections(num_special_cutoff, cache_override)

integer, intent(in) :: num_special_cutoff
logical, intent(in) :: cache_override

integer :: i

if (adjust_obs_impact) then
   call allocate_impact_table(obs_impact_table)
   call read_impact_table(obs_impact_filename, obs_impact_table, allow_any_impact_values, "allow_any_impact_values")
   call error_handler(E_MSG, 'assim_tools_init:', &
                      'Using observation impact table from file "'//trim(obs_impact_filename)//'"')
endif

write(msgstring,  '(A,F18.6)') 'The cutoff namelist value is ', cutoff
write(msgstring2, '(A)') 'cutoff is the localization half-width parameter,'
write(msgstring3, '(A,F18.6)') 'so the effective localization radius is ', cutoff*2.0_r8
call error_handler(E_MSG,'assim_tools_init:', msgstring, text2=msgstring2, text3=msgstring3)

if (has_special_cutoffs) then
   call error_handler(E_MSG, '', '')
   call error_handler(E_MSG,'assim_tools_init:','Observations with special localization treatment:')
   call error_handler(E_MSG,'assim_tools_init:','(type name, specified cutoff distance, effective localization radius)')

   do i = 1, num_special_cutoff
      write(msgstring, '(A32,F18.6,F18.6)') special_localization_obs_types(i), &
            special_localization_cutoffs(i), special_localization_cutoffs(i)*2.0_r8
      call error_handler(E_MSG,'assim_tools_init:', msgstring)
   end do
   call error_handler(E_MSG,'assim_tools_init:','all other observation types will use the default cutoff distance')
   call error_handler(E_MSG, '', '')
endif

if (cache_override) then
   call error_handler(E_MSG,'assim_tools_init:','Disabling the close obs caching because specialized localization')
   call error_handler(E_MSG,'assim_tools_init:','distances are enabled. ')
endif

if(adaptive_localization_threshold > 0) then
   write(msgstring, '(A,I10,A)') 'Using adaptive localization, threshold ', &
                                  adaptive_localization_threshold, ' obs'
   call error_handler(E_MSG,'assim_tools_init:', msgstring)
   if(adaptive_cutoff_floor > 0.0_r8) then
      write(msgstring, '(A,F18.6)') 'Minimum cutoff will not go below ', &
                                     adaptive_cutoff_floor
      call error_handler(E_MSG,'assim_tools_init:', 'Using adaptive localization cutoff floor.', &
                         text2=msgstring)
   endif
endif

if(output_localization_diagnostics) then
   call error_handler(E_MSG,'assim_tools_init:', 'Writing localization diagnostics to file:')
   call error_handler(E_MSG,'assim_tools_init:', trim(localization_diagnostics_file))
endif

if(sampling_error_correction) then
   call error_handler(E_MSG,'assim_tools_init:', 'Using Sampling Error Correction')
endif

if (task_count() > 1) then
    if(distribute_mean) then
       msgstring  = 'Distributing one copy of the ensemble mean across all tasks'
       msgstring2 = 'uses less memory per task but may run slower if doing vertical '
    else
       msgstring  = 'Replicating a copy of the ensemble mean on every task'
       msgstring2 = 'uses more memory per task but may run faster if doing vertical '
    endif
    call error_handler(E_MSG,'assim_tools_init:', msgstring, text2=msgstring2, &
                       text3='coordinate conversion; controlled by namelist item "distribute_mean"')
endif

if (has_vertical_choice()) then
   if (.not. vertical_localization_on()) then
      msgstring = 'Not doing vertical localization, no vertical coordinate conversion required'
      call error_handler(E_MSG,'assim_tools_init:', msgstring)
   else
      msgstring = 'Doing vertical localization, vertical coordinate conversion may be required'
      if (convert_all_state_verticals_first) then
         msgstring2 = 'Converting all state vector verticals to localization coordinate first.'
      else
         msgstring2 = 'Converting all state vector verticals only as needed.'
      endif
      if (convert_all_obs_verticals_first) then
         msgstring3 = 'Converting all observation verticals to localization coordinate first.'
      else
         msgstring3 = 'Converting all observation verticals only as needed.'
      endif
      call error_handler(E_MSG,'assim_tools_init:', msgstring, text2=msgstring2, text3=msgstring3)
   endif
endif

end subroutine log_namelist_selections

!===========================================================
! TEST FUNCTIONS BELOW THIS POINT
!-----------------------------------------------------------
!> test get_state_meta_data
!> Write out the resutls of get_state_meta_data for each task
!> They should be the same as the Trunk version
subroutine test_get_state_meta_data(locations, num_vars)

type(location_type), intent(in) :: locations(:)
integer,             intent(in) :: num_vars

character*20  :: task_str !< string to hold the task number
character*129 :: file_meta !< output file name
character(len=128) :: locinfo
integer :: i

write(task_str, '(i10)') my_task_id()
file_meta = TRIM('test_get_state_meta_data' // TRIM(ADJUSTL(task_str)))

open(15, file=file_meta, status = 'unknown')

do i = 1, num_vars
   call write_location(-1, locations(i), charstring=locinfo)
   write(15,*) trim(locinfo)
enddo

close(15)


end subroutine test_get_state_meta_data

!--------------------------------------------------------
!> dump out the copies array for the state ens handle
subroutine test_state_copies(state_ens_handle, information)

type(ensemble_type), intent(in) :: state_ens_handle
character(len=*),        intent(in) :: information

character*20  :: task_str !< string to hold the task number
character*129 :: file_copies !< output file name
integer :: i

write(task_str, '(i10)') state_ens_handle%my_pe
file_copies = TRIM('statecopies_'  // TRIM(ADJUSTL(information)) // '.' // TRIM(ADJUSTL(task_str)))
open(15, file=file_copies, status ='unknown')

do i = 1, state_ens_handle%num_copies - state_ens_handle%num_extras
   write(15, *) state_ens_handle%copies(i,:)
enddo

close(15)

end subroutine test_state_copies

!--------------------------------------------------------
!> dump out the distances calculated in get_close_obs
subroutine test_close_obs_dist(distances, num_close, ob)

real(r8), intent(in) :: distances(:) !< array of distances calculated in get_close
integer,  intent(in) :: num_close !< number of close obs
integer,  intent(in) :: ob

character*20  :: task_str !< string to hold the task number
character*20  :: ob_str !< string to hold ob number
character*129 :: file_dist !< output file name
integer :: i

write(task_str, '(i10)') my_task_id()
write(ob_str, '(i20)') ob
file_dist = TRIM('distances'   // TRIM(ADJUSTL(task_str)) // '.' // TRIM(ADJUSTL(ob_str)))
open(15, file=file_dist, status ='unknown')

write(15, *) num_close

do i = 1, num_close
   write(15, *) distances(i)
enddo

close(15)

end subroutine test_close_obs_dist

!> @}

!========================================================================
! end module assim_tools_mod
!========================================================================

end module assim_tools_mod

