! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
!
!
!!!!!!!!!!!!!!!!!!!!!!  PLEASE READ BEFORE USING  !!!!!!!!!!!!!!!!!!!!!!
!
! This version of assim_tools_mod supports ongoing developments toward a
! localized particle filter (Poterjoy 2016; Poterjoy and Anderson 2016;
! Poterjoy et al. 2019; Poterjoy 2021). It should NOT be used with 'filter 
! kind' other than 9 or with any sort of adaptive prior/posterior inflation.
! A breakdown of namelist options specific to the local PF are described
! below. They include options to use a hybrid PF-EAKF update, which performs
! a partial assimilation of observations with the PF followed by a second
! update using the EAKF (Poterjoy 2021). Please contact the developer
! (Jon Poterjoy, poterjoy@umd.edu) before using.
!
!>  A variety of operations required by assimilation.
module assim_tools_mod

!> \defgroup assim_tools assim_tools_mod
!>
!> @{
use      types_mod,       only : r8, i8, digits12, PI, missing_r8

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
                                 read_mpi_timer, task_sync, get_global_max

use adaptive_inflate_mod, only : do_obs_inflate,  do_single_ss_inflate, do_ss_inflate,    &
                                 do_varying_ss_inflate,                                   &
                                 update_inflation,                                        &
                                 inflate_ens, adaptive_inflate_type,                      &
                                 deterministic_inflate, solve_quadratic

use time_manager_mod,     only : time_type, get_time

use assim_model_mod,      only : get_state_meta_data,                                     &
                                 get_close_obs,         get_close_state,                  &
                                 convert_vertical_obs,  convert_vertical_state

use distributed_state_mod, only : create_mean_window, free_mean_window

use quality_control_mod, only : good_dart_qc, DARTQC_FAILED_VERT_CONVERT

implicit none
private

public :: filter_assim, &
          set_assim_tools_trace, &
          test_state_copies, &
          update_ens_from_weights  ! Jeff thinks this routine is in the wild.

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
real(r8), parameter    :: small = epsilon(1.0_r8)   ! threshold for avoiding NaNs/Inf

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

character(len=*), parameter :: source = 'assim_tools_mod.pf.f90'

!============================================================================

!---- namelist with default values

! Filter kind selects type of observation space filter
!      1 = EAKF filter
!      2 = ENKF
!      3 = Kernel filter
!      4 = particle filter
!      5 = random draw from posterior
!      6 = deterministic draw from posterior with fixed kurtosis
!      8 = Rank Histogram Filter (see Anderson 2011)
!      9 = Localized particle filter (Poterjoy Nov. 2014)
!
!  special_localization_obs_types -> Special treatment for the specified observation types
!  special_localization_cutoffs   -> Different cutoff value for each specified obs type
!
!
! PF namelist variables: 
!      frac_neff -> Effective ensemble size for regularization 
!                   (set as fraction of ensemble size)
!      pf_alpha -> Mixing coefficient for PF update step
!      pf_kddm -> Flag for turning on probability mapping step
!      sampling_weighted_prior -> Flag for turning on resampling from weighted prior
!      pf_enkf_hybrid -> Flag for turning on hybrid PF-EnKF update
!      min_residual -> Min residual for PF iterations (also used as hybrid parameter)
!      pf_maxiter -> Max number of PF iterations for tempering
!      pf_kf_rtps_coeff -> Relaxation to prior spread factor used for EnKF part of hybrid

integer  :: filter_kind                     = 9
real(r8) :: cutoff                          = 0.3_r8
logical  :: sort_obs_inc                    = .false.
logical  :: spread_restoration              = .false.
logical  :: sampling_error_correction       = .false.
integer  :: adaptive_localization_threshold = -1
real(r8) :: adaptive_cutoff_floor           = 0.0_r8
integer  :: print_every_nth_obs             = 0
real(r8) :: frac_neff                       = 0.20_r8
real(r8) :: pf_alpha                        = 0.30_r8
integer  :: pf_kddm                         = 0
logical  :: sampling_weighted_prior         = .true.
logical  :: pf_enkf_hybrid                  = .false.
real(r8) :: min_residual                    = 0.5_r8
integer  :: pf_maxiter                      = 3
real(r8) :: pf_kf_rtps_coeff                = 0.0_r8

! since this is in the namelist, it has to have a fixed size.
integer, parameter   :: MAX_ITEMS = 300
character(len = 129) :: special_localization_obs_types(MAX_ITEMS)
real(r8)             :: special_localization_cutoffs(MAX_ITEMS)

logical              :: output_localization_diagnostics = .false.
character(len = 129) :: localization_diagnostics_file = "localization_diagnostics"

! Following only relevant for filter_kind = 8
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

! New namelist variables added for local PF
namelist / assim_tools_nml / filter_kind, cutoff, sort_obs_inc, &
   spread_restoration, sampling_error_correction,                          &
   adaptive_localization_threshold, adaptive_cutoff_floor,                 &
   print_every_nth_obs, rectangular_quadrature, gaussian_likelihood_tails, &
   pf_kddm, frac_neff, pf_alpha, sampling_weighted_prior, pf_enkf_hybrid,  &
   min_residual, pf_maxiter, pf_kf_rtps_coeff,                             &
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

! FOR NOW, can only do spread restoration with filter option 1 (need to extend this)
if(spread_restoration .and. .not. filter_kind == 1) then
   write(msgstring, *) 'cannot combine spread_restoration and filter_kind ', filter_kind
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
type(adaptive_inflate_type), intent(inout) :: inflate
integer,                     intent(in)    :: ENS_MEAN_COPY, ENS_SD_COPY, ENS_INF_COPY
integer,                     intent(in)    :: ENS_INF_SD_COPY
integer,                     intent(in)    :: OBS_KEY_COPY, OBS_GLOBAL_QC_COPY
integer,                     intent(in)    :: OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END
integer,                     intent(in)    :: OBS_PRIOR_VAR_START, OBS_PRIOR_VAR_END
logical,                     intent(in)    :: inflate_only

! changed the ensemble sized things here to allocatable

real(r8) :: obs_prior(ens_size), obs_inc(ens_size), increment(ens_size)
real(r8) :: reg_factor, impact_factor
real(r8) :: net_a(num_groups), reg_coef(num_groups), correl(num_groups)
real(r8) :: cov_factor, obs(1), obs_err_var, my_inflate, my_inflate_sd
real(r8) :: varying_ss_inflate, varying_ss_inflate_sd
real(r8) :: ss_inflate_base, obs_qc, cutoff_rev, cutoff_orig
real(r8) :: gamma, ens_obs_mean, ens_obs_var, ens_var_deflate
real(r8) :: r_mean, r_var
real(r8) :: orig_obs_prior_mean(num_groups), orig_obs_prior_var(num_groups)
real(r8) :: obs_prior_mean(num_groups), obs_prior_var(num_groups)
real(r8) :: diff_sd, outlier_ratio
real(r8) :: vertvalue_obs_in_localization_coord, whichvert_real
real(r8), allocatable :: close_obs_dist(:)
real(r8), allocatable :: close_state_dist(:)
real(r8), allocatable :: last_close_obs_dist(:)
real(r8), allocatable :: last_close_state_dist(:)

integer(i8) :: state_index
integer(i8), allocatable :: my_state_indx(:)
integer(i8), allocatable :: my_obs_indx(:)

character(8)  :: date
character(10) :: time
real(r8), parameter :: beta_max = 1E10_r8

real(r8) :: orig_obs_prior(ens_size), min_res, max_res
real(r8) :: ens_init(ens_size,ens_handle%my_num_vars), pf_infl(obs_ens_handle%my_num_vars)
real(r8) :: obs_ens_init(ens_size,obs_ens_handle%my_num_vars), obs_err_infl, hw(ens_size), wo(ens_size)
real(r8) :: lw(ens_size,ens_handle%my_num_vars), lhw(ens_size,obs_ens_handle%my_num_vars), wp(ens_size,obs_ens_handle%my_num_vars)
real(r8) :: beta(ens_handle%my_num_vars), beta_y(obs_ens_handle%my_num_vars)
real(r8) :: w(ens_size), ens_mean, ens_var, ws, d(ens_size), wt(ens_size), temp1, temp2, beta_hw, neff
integer  :: indx(ens_size), filter_kind_orig
real(r8) :: res(ens_handle%my_num_vars), res_y(obs_ens_handle%my_num_vars)
integer  :: iter, maxiter

integer :: my_num_obs, i, j, n, owner, owners_index, my_num_state
integer :: this_obs_key, obs_mean_index, obs_var_index
integer :: grp_beg(num_groups), grp_end(num_groups), grp_size, grp_bot, grp_top, group
integer :: num_close_obs, obs_index, num_close_states
integer :: total_num_close_obs, last_num_close_obs, last_num_close_states
integer :: base_obs_kind, base_obs_type, nth_obs
integer :: num_close_obs_cached, num_close_states_cached
integer :: num_close_obs_calls_made, num_close_states_calls_made
integer :: localization_unit, secs, days, rev_num_close_obs
integer :: whichvert_obs_in_localization_coord
integer :: istatus
integer, allocatable :: close_obs_ind(:)
integer, allocatable :: close_state_ind(:)
integer, allocatable :: last_close_obs_ind(:)
integer, allocatable :: last_close_state_ind(:)
integer, allocatable :: my_obs_kind(:)
integer, allocatable :: my_obs_type(:)
integer, allocatable :: my_state_kind(:)
integer, allocatable :: vstatus(:)

character(len = 200) :: base_loc_text   ! longer than longest location formatting possible

type(location_type)  :: base_obs_loc, last_base_obs_loc, last_base_states_loc
type(location_type)  :: dummyloc
type(location_type), allocatable :: my_obs_loc(:)
type(location_type), allocatable :: my_state_loc(:)

type(get_close_type) :: gc_obs, gc_state
type(obs_type)       :: observation
type(obs_def_type)   :: obs_def
type(time_type)      :: obs_time, this_obs_time

logical :: do_adapt_inf_update
logical :: allow_missing_in_state
logical :: local_single_ss_inflate
logical :: local_varying_ss_inflate
logical :: local_ss_inflate
logical :: local_obs_inflate

type(location_type) :: lc(1)
integer             :: kd(1)

! timing related vars:
! set timing(N) true to collect and print timing info
integer, parameter :: Ntimers = 5
integer, parameter :: MLOOP  = 1  ! main assimilation loop
integer, parameter :: LG_GRN = 2  ! large section timings
integer, parameter :: SM_GRN = 3  ! inner loops - use carefully!
integer, parameter :: GC     = 4  ! get_close() related loops
logical        :: timing(Ntimers)   ! enable or disable w/ this
real(digits12) :: t_base(Ntimers)   ! storage for time info
integer(i8)    :: t_items(Ntimers)  ! count of number of calls
integer(i8)    :: t_limit(Ntimers)  ! limit on number printed
real(digits12), allocatable :: elapse_array(:)

integer, allocatable :: n_close_state_items(:), n_close_obs_items(:)

! timing disabled by default
timing(:)  = .false.
t_base(:)  = 0.0_r8
t_items(:) = 0_i8
t_limit(:) = 0_i8

! how about this?  look for imbalances in the tasks
allocate(n_close_state_items(obs_ens_handle%num_vars), &
         n_close_obs_items(  obs_ens_handle%num_vars))

! turn these on carefully - they can generate a lot of output!
! also, to be readable - at least with ifort:
!  setenv FORT_FMT_RECL 1024
! so output lines don't wrap.

!timing(MLOOP)  = .true.
!timing(LG_GRN) = .true.

if (timing(MLOOP)) allocate(elapse_array(obs_ens_handle%num_vars))

! use maxitems limit here or drown in output.
!timing(SM_GRN) = .false.
!t_limit(SM_GRN) = 4_i8

!timing(GC) = .true.
!t_limit(GC) = 4_i8


! allocate rather than dump all this on the stack
allocate(close_obs_dist(     obs_ens_handle%my_num_vars), &
         last_close_obs_dist(obs_ens_handle%my_num_vars), &
         close_obs_ind(      obs_ens_handle%my_num_vars), &
         last_close_obs_ind( obs_ens_handle%my_num_vars), &
         vstatus(            obs_ens_handle%my_num_vars), &
         my_obs_indx(        obs_ens_handle%my_num_vars), &
         my_obs_kind(        obs_ens_handle%my_num_vars), &
         my_obs_type(        obs_ens_handle%my_num_vars), &
         my_obs_loc(         obs_ens_handle%my_num_vars))

allocate(close_state_dist(     ens_handle%my_num_vars), &
         last_close_state_dist(ens_handle%my_num_vars), &
         close_state_ind(      ens_handle%my_num_vars), &
         last_close_state_ind( ens_handle%my_num_vars), &
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

! filter kinds 1 and 8 return sorted increments, however non-deterministic
! inflation can scramble these. the sort is expensive, so help users get better
! performance by rejecting namelist combinations that do unneeded work.
if (sort_obs_inc) then
   if(deterministic_inflate(inflate) .and. ((filter_kind == 1) .or. (filter_kind == 8))) then
      write(msgstring,  *) 'With a deterministic filter [assim_tools_nml:filter_kind = ',filter_kind,']'
      write(msgstring2, *) 'and deterministic inflation [filter_nml:inf_deterministic = .TRUE.]'
      write(msgstring3, *) 'assim_tools_nml:sort_obs_inc = .TRUE. is not needed and is expensive.'
      call error_handler(E_MSG,'', '')  ! whitespace
      call error_handler(E_MSG,'WARNING filter_assim:', msgstring, source, &
                         text2=msgstring2,text3=msgstring3)
      call error_handler(E_MSG,'', '')  ! whitespace
      sort_obs_inc = .FALSE.
   endif
endif

!GSR open the dignostics file
if(output_localization_diagnostics .and. my_task_id() == 0) then
  localization_unit = open_file(localization_diagnostics_file, action = 'append')
endif

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
   if (timing(LG_GRN)) call start_timer(t_base(LG_GRN))
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
   if (timing(LG_GRN)) call read_timer(t_base(LG_GRN), 'convert_vertical_obs')
endif

! Get info on my number and indices for state
my_num_state = get_my_num_vars(ens_handle)
call get_my_vars(ens_handle, my_state_indx)

! Get the location and kind of all my state variables
if (timing(LG_GRN)) call start_timer(t_base(LG_GRN))
do i = 1, ens_handle%my_num_vars
   call get_state_meta_data(my_state_indx(i), my_state_loc(i), my_state_kind(i))
end do
if (timing(LG_GRN)) call read_timer(t_base(LG_GRN), 'get_state_meta_data')

!call test_get_state_meta_data(my_state_loc, ens_handle%my_num_vars)

!> optionally convert all state location verticals
if (convert_all_state_verticals_first .and. is_doing_vertical_conversion) then
   if (timing(LG_GRN)) call start_timer(t_base(LG_GRN))
   if (ens_handle%my_num_vars > 0) then
      call convert_vertical_state(ens_handle, ens_handle%my_num_vars, my_state_loc, my_state_kind,  &
                                  my_state_indx, get_vertical_localization_coord(), istatus)
   endif
   if (timing(LG_GRN)) call read_timer(t_base(LG_GRN), 'convert_vertical_state')
endif

! PAR: MIGHT BE BETTER TO HAVE ONE PE DEDICATED TO COMPUTING
! INCREMENTS. OWNING PE WOULD SHIP IT'S PRIOR TO THIS ONE
! BEFORE EACH INCREMENT.

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

! The computations in the two get_close_maxdist_init are redundant

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

! Set parameters for iterative PF updates
filter_kind_orig = filter_kind

if (filter_kind == 9) then

  if (pf_enkf_hybrid) then
     min_res = min_residual
  else
     min_res = 0.0_r8
  end if

! Set number of iterations for PF step
  res = 1.0_r8 - min_res
  res_y = 1.0_r8 - min_res
  pf_infl = 1.0_r8

  ! Add extra iteration for hybrid
  maxiter = pf_maxiter + 1

else
   maxiter = 1
end if

if (close_obs_caching) then
   ! Initialize last obs and state get_close lookups, to take advantage below
   ! of sequential observations at the same location (e.g. U,V, possibly T,Q)
   ! (this is getting long enough it probably should go into a subroutine. nsc.)
   last_base_obs_loc           = set_location_missing()
   last_base_states_loc        = set_location_missing()
   last_num_close_obs          = -1
   last_num_close_states       = -1
   last_close_obs_ind(:)       = -1
   last_close_state_ind(:)     = -1
   last_close_obs_dist(:)      = 888888.0_r8   ! something big, not small
   last_close_state_dist(:)    = 888888.0_r8   ! ditto
   num_close_obs_cached        = 0
   num_close_states_cached     = 0
   num_close_obs_calls_made    = 0
   num_close_states_calls_made = 0
endif

allow_missing_in_state = get_missing_ok_status()

! use MLOOP for the overall outer loop times; LG_GRN is for
! sections inside the overall loop, including the total time
! for the state_update and obs_update loops.  use SM_GRN for
! sections inside those last 2 loops and be careful - they will
! be called nobs * nstate * ntasks.


! Main loop for iterative PF updates
ITERATIONS: do iter = 1,maxiter

  ens_init = ens_handle%copies(1:ens_size,:)
  obs_ens_init = obs_ens_handle%copies(1:ens_size,:)

  ! Initiate local PF arrays
  if (filter_kind == 9) then

   ! Need to store initial model space and observation space ensemble and initialize 
   ! weights before observation loop. For now, store weights for each state variable. 
   ! In the future, weights are needed only for each model grid point.
   lw  = 0.0_r8
   lhw = 0.0_r8
   wp = 1.0_r8 / ens_size
   beta = beta_max
   beta_y = beta_max

   ! Regularization strategy requires calculating the -log() of localized
   ! obs-and model-space weights prior to DA step
   if (timing(LG_GRN)) call start_timer(t_base(LG_GRN))

   REGULARIZATION: do i = 1, obs_ens_handle%num_vars

      ! Every pe has information about the global obs sequence
      call get_obs_from_key(obs_seq, keys(i), observation)
      call get_obs_def(observation, obs_def)
      base_obs_loc = get_obs_def_location(obs_def)
      obs_err_var = get_obs_def_error_variance(obs_def)
      base_obs_type = get_obs_def_type_of_obs(obs_def)
      if (base_obs_type > 0) then
         base_obs_kind = get_quantity_for_type_of_obs(base_obs_type)
      else
         call get_state_meta_data(-1 * int(base_obs_type,i8), dummyloc, base_obs_kind)
      endif
   
      ! Get the value of the observation
      call get_obs_values(observation, obs, obs_val_index)
   
      ! Find out who has this observation and where it is
      call get_var_owner_index(ens_handle, int(i,i8), owner, owners_index)

      ! Owner calculates weights for current ob
      if(ens_handle%my_pe == owner) then
   
        obs_qc = obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, owners_index)

        if(nint(obs_qc) == 0) then

          ! Likelihood calculations
          orig_obs_prior = obs_ens_init(1:ens_size, owners_index)
          d = (obs(1) - orig_obs_prior)**2 / (2.0_r8*obs_err_var)
          d = d - minval(d)

          ! Determine whether to skip ob
          hw = exp( -d )
          hw = hw / sum(hw)

          if (1.0_r8 > ens_size * 0.98_r8 *sum(hw**2) ) then
            ! write(*,*) 'Skipping with Neff =',1.0_r8 / sum(hw**2)
            obs_qc = 1
            obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, owners_index) = 1
          else

            temp1 = max( 1.0_r8, maxval(d)/80.0_r8 )
            call pf_regularization_minw(d/temp1, 1E-10_r8 / ens_size, ens_size, obs_err_infl)

            pf_infl(owners_index) = obs_err_infl*temp1

            hw = exp( -d/pf_infl(owners_index) )
            hw = hw / sum(hw)

            ! Save inflation and weights
            wp(1:ens_size,owners_index) = hw

          end if

        end if

        call broadcast_send(map_pe_to_task(ens_handle, owner), hw, scalar1=obs_qc)

      else

         call broadcast_recv(map_pe_to_task(ens_handle, owner), hw, scalar1=obs_qc)

      end if

      ! Skip ob if flagged by qc
      if(nint(obs_qc) /= 0) cycle REGULARIZATION

      ! Get localization coefficients for nearby obs- and model-space priors
      if (.not. close_obs_caching) then
         call get_close_obs(gc_obs, base_obs_loc, base_obs_type, &
                            my_obs_loc, my_obs_kind, my_obs_type, &
                            num_close_obs, close_obs_ind, close_obs_dist, ens_handle)
      else
 
         if (base_obs_loc == last_base_obs_loc) then
            num_close_obs     = last_num_close_obs
            close_obs_ind(:)  = last_close_obs_ind(:)
            close_obs_dist(:) = last_close_obs_dist(:)
            num_close_obs_cached = num_close_obs_cached + 1
         else
            call get_close_obs(gc_obs, base_obs_loc, base_obs_type, &
                               my_obs_loc, my_obs_kind, my_obs_type, &
                               num_close_obs, close_obs_ind, close_obs_dist, ens_handle)
            last_base_obs_loc      = base_obs_loc
            last_num_close_obs     = num_close_obs
            last_close_obs_ind(:)  = close_obs_ind(:)
            last_close_obs_dist(:) = close_obs_dist(:)
            num_close_obs_calls_made = num_close_obs_calls_made +1
         endif
      endif

      if (base_obs_type > 0) then
         cutoff_orig = cutoff_list(base_obs_type)
      else
         cutoff_orig = cutoff
      endif

      cutoff_rev = cutoff_orig

      if (.not. close_obs_caching) then
         call get_close_state(gc_state, base_obs_loc, base_obs_type, &
                              my_state_loc, my_state_kind, my_state_indx, &
                              num_close_states, close_state_ind, close_state_dist, ens_handle)
      else
         if (base_obs_loc == last_base_states_loc) then
            num_close_states    = last_num_close_states
            close_state_ind(:)  = last_close_state_ind(:)
            close_state_dist(:) = last_close_state_dist(:)
            num_close_states_cached = num_close_states_cached + 1
         else
            call get_close_state(gc_state, base_obs_loc, base_obs_type, &
                                 my_state_loc, my_state_kind, my_state_indx, &
                                 num_close_states, close_state_ind, close_state_dist, ens_handle)
            last_base_states_loc     = base_obs_loc
            last_num_close_states    = num_close_states
            last_close_state_ind(:)  = close_state_ind(:)
            last_close_state_dist(:) = close_state_dist(:)
            num_close_states_calls_made = num_close_states_calls_made + 1 
         endif
      endif

      ! Redundant part of weight calculation
      wt = ens_size*hw - 1.0_r8

      ! Calculate logw for prior states that are potentially close
      STATE_WEIGHTS: do j = 1, num_close_states
         state_index = close_state_ind(j)

         ! the "any" is an expensive test when you do it for every ob.  don't test
         ! if we know there aren't going to be missing values in the state.
         if ( allow_missing_in_state ) then
            ! Some models can take evasive action if one or more of the ensembles have
            ! a missing value. Generally means 'do nothing' (as opposed to DIE)
            if (any(ens_handle%copies(1:ens_size, state_index) == MISSING_R8)) cycle STATE_WEIGHTS
         endif
     
         ! Compute the distance and covariance factor 
         cov_factor = comp_cov_factor(close_state_dist(j), cutoff_rev, &
            base_obs_loc, base_obs_type, my_state_loc(state_index), my_state_kind(state_index))

         ! Update only when prior variance is not zero
         if ( maxval( ens_handle%copies(1:ens_size, state_index)) /= &
              minval( ens_handle%copies(1:ens_size, state_index)) ) then

            ! Take running sum of -log() of weights
            if (cov_factor == 1.0_r8) then
              d = log(ens_size*hw)
            else
              d = wt*cov_factor
              do n = 1,ens_size
                if (abs(d(n)) > 0.1_r8) then
                  d(n) = log( d(n) + 1.0_r8 )
                end if
              end do
            end if
            lw(1:ens_size,state_index) = lw(1:ens_size,state_index) - d
            lw(1:ens_size,state_index) = lw(1:ens_size,state_index) - minval(lw(1:ens_size,state_index))

         end if 

      end do STATE_WEIGHTS

      ! Calculate logw for prior obs-space states that are potentially close
      OBS_WEIGHTS: do j = 1, num_close_obs
         obs_index = close_obs_ind(j)

         ! If the forward observation operator failed, no need to 
         ! update the unassimilated observations 
         if (any(obs_ens_handle%copies(1:ens_size, obs_index) == MISSING_R8)) cycle OBS_WEIGHTS

         ! Compute the distance and the covar_factor
         cov_factor = comp_cov_factor(close_obs_dist(j), cutoff_rev, &
            base_obs_loc, base_obs_type, my_obs_loc(obs_index), my_obs_kind(obs_index))

         ! Update only when prior variance is not zero
         if ( maxval( obs_ens_handle%copies(1:ens_size, obs_index)) /= &
              minval( obs_ens_handle%copies(1:ens_size, obs_index)) ) then

            ! Take running sum of -log() of weights
            if (cov_factor == 1.0_r8) then
              d = log(ens_size*hw)
            else
              d = wt*cov_factor
              do n = 1,ens_size
                if (abs(d(n)) > 0.1_r8) then
                  d(n) = log( d(n) + 1.0_r8 )
                end if
              end do
            end if
            lhw(1:ens_size,obs_index) = lhw(1:ens_size,obs_index) - d
            lhw(1:ens_size,obs_index) = lhw(1:ens_size,obs_index) - minval(lhw(1:ens_size,obs_index))
         end if  

      end do OBS_WEIGHTS

   end do REGULARIZATION


   ! Calculate model-space regularization coefficients 
   max_res = 0.0_r8
   do i = 1,ens_handle%my_num_vars

      if ( (res(i) > 0.0_r8) ) then

         call pf_regularization(lw(1:ens_size,i),ens_size,frac_neff*ens_size,beta(i),beta_max)

         ! Fix beta if its inverse exceeds residual
         if (res(i) <= 1.0_r8/beta(i)) then
            beta(i) = 1.0_r8/res(i)
            res(i) = 0.0_r8
         else
            res(i) = res(i) - 1.0_r8/beta(i)
         end if

         ! Store residuals
         beta(i) = min(beta(i),beta_max)
         max_res = max(res(i),max_res)

      else

         beta(i) = beta_max

      end if

   end do

   ! Calculate obs-space regularization coefficients 
   do i = 1,obs_ens_handle%my_num_vars

      obs_qc = obs_ens_handle%copies(OBS_GLOBAL_QC_COPY,i)

      if ( (res_y(i) > 0.0_r8) .and. (nint(obs_qc) == 0) ) then

         call pf_regularization(lhw(1:ens_size,i),ens_size,frac_neff*ens_size,beta_y(i),beta_max)

         ! Fix beta if its inverse exceeds residual
         if (res_y(i) <= 1.0_r8/beta_y(i)) then
            beta_y(i) = 1.0_r8/res_y(i)
            res_y(i) = 0.0_r8
         else
            res_y(i) = res_y(i) - 1.0_r8/beta_y(i)
         end if

         ! Store residuals
         beta_y(i) = min(beta_y(i),beta_max)

      else

         beta_y(i) = beta_max

      end if

   end do

   call get_global_max(max_res)

   if (timing(LG_GRN)) call read_timer(t_base(LG_GRN), 'regularization')

   ! Reset -log() weights to zero for sequential udpate step
   lw  = 0.0_r8
   lhw = 0.0_r8

  end if  ! local PF


  ! Loop through all the (global) observations sequentially
  SEQUENTIAL_OBS: do i = 1, obs_ens_handle%num_vars

   if (timing(MLOOP))  call start_timer(t_base(MLOOP))
   if (timing(LG_GRN)) call start_timer(t_base(LG_GRN))

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

      ! Inflate obs error variance for EnKF step in hybrid
      if (filter_kind == 1 .and. iter > 1) then
         obs_err_infl = 1.0_r8 / min_res
         obs_qc = obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, owners_index)
         obs_err_var = obs_err_var*obs_err_infl
         if ( obs_err_infl > 500 ) obs_qc = 1
      end if

      ! Only value of 0 for DART QC field should be assimilated
      IF_QC_IS_OKAY: if(nint(obs_qc) == 0) then

        obs_prior = obs_ens_handle%copies(1:ens_size, owners_index)

        ! Get current PF obs-space weights
        if (filter_kind == 9) then

          ! Prior for current ob
          if (sampling_weighted_prior) then
            obs_prior = obs_ens_init(1:ens_size, owners_index)
          end if

          ! Contribution of current ob for vector weight calculations
          hw = wp(1:ens_size, owners_index)

          ! Obs-space weights used for sampling
          if (sampling_weighted_prior) then

            w = lhw(1:ens_size,owners_index) - log( ens_size*hw )
            w = w - minval(w)
            w = exp( - w / beta_y(owners_index) )
            w = w / sum(w)

          else

            ! Obs error inflation 
            obs_err_infl = pf_infl(owners_index)*beta_y(owners_index)

            ! Scalar weights for resampling particles
            d = (obs(1) - obs_prior)**2 / (2.0_r8*obs_err_var*obs_err_infl)
            w = exp( -d )
            if ( sum(w) == 0.0_r8 ) then
              n = minloc(d, 1, mask=d.gt.0)
              w(n) = 1.0_r8
            end if
            w = w / sum(w)

          end if

          if (sampling_weighted_prior) then
            ws = 1.0_r8 / 2.0_r8
            w = w**ws
            w = w / sum(w)
          end if

          ! Compute obs space prior information for adaptive inflation
          if(local_varying_ss_inflate) then
             orig_obs_prior_mean = obs_ens_handle%copies(OBS_PRIOR_MEAN_START: &
                OBS_PRIOR_MEAN_END, owners_index)
             orig_obs_prior_var  = obs_ens_handle%copies(OBS_PRIOR_VAR_START:  &
                OBS_PRIOR_VAR_END, owners_index)
          endif

        else

          ! Compute the prior mean and variance for this observation
          orig_obs_prior_mean = obs_ens_handle%copies(OBS_PRIOR_MEAN_START: &
             OBS_PRIOR_MEAN_END, owners_index)

          orig_obs_prior_var  = obs_ens_handle%copies(OBS_PRIOR_VAR_START:  &
             OBS_PRIOR_VAR_END, owners_index)

          ! Compute observation space increments for each group
          do group = 1, num_groups
             grp_bot = grp_beg(group)
             grp_top = grp_end(group)
             call obs_increment(obs_prior(grp_bot:grp_top), grp_size, obs(1), &
                obs_err_var, obs_inc(grp_bot:grp_top), inflate, my_inflate,   &
                my_inflate_sd, net_a(group))
          end do

        endif

        ! Compute updated values for single state space inflation
        SINGLE_SS_INFLATE: if(local_single_ss_inflate) then
           ss_inflate_base = ens_handle%copies(ENS_SD_COPY, 1)
           ! Update for each group separately
           do group = 1, num_groups
             ! If either inflation or sd is not positive, not really doing inflation
             if(my_inflate > 0.0_r8 .and. my_inflate_sd > 0.0_r8) then
                 ! For case with single spatial inflation, use gamma = 1.0_r8
                  ! See adaptive inflation module for details
                  gamma = 1.0_r8
                  ! Deflate the inflated variance; required for efficient single pass
                  ! This is one of many places that assumes linear state/obs relation
                  ! over range of ensemble; Essentially, we are removing the inflation
                  ! which has already been applied in filter to see what inflation should
                  ! have been needed.
                  ens_obs_mean = orig_obs_prior_mean(group)
                  ens_obs_var = orig_obs_prior_var(group)
                  ! gamma is hardcoded as 1.0, so no test is needed here.
                  ens_var_deflate = ens_obs_var / &
                     (1.0_r8 + gamma*(sqrt(ss_inflate_base) - 1.0_r8))**2

                  ! If this is inflate_only (i.e. posterior) remove impact of this obs.
                  ! This is simulating independent observation by removing its impact.
                  if(inflate_only .and. &
                        ens_var_deflate               > small .and. &
                        obs_err_var                   > small .and. &
                        obs_err_var - ens_var_deflate > small ) then
                     r_var = 1.0_r8 / (1.0_r8 / ens_var_deflate - 1.0_r8 / obs_err_var)
                     r_mean = r_var *(ens_obs_mean / ens_var_deflate - obs(1) / obs_err_var)
                  else
                     r_var = ens_var_deflate
                     r_mean = ens_obs_mean
                  endif

                  if (timing(SM_GRN)) call start_timer(t_base(SM_GRN), t_items(SM_GRN), t_limit(SM_GRN), do_sync=.false.)
                  ! Update the inflation value
                  call update_inflation(inflate, my_inflate, my_inflate_sd, &
                     r_mean, r_var, grp_size, obs(1), obs_err_var, gamma)
                  if (timing(SM_GRN)) call read_timer(t_base(SM_GRN), 'update_inflation_C', &
                                                      t_items(SM_GRN), t_limit(SM_GRN), do_sync=.false.)
               endif
            end do
         endif SINGLE_SS_INFLATE

      endif IF_QC_IS_OKAY

      !Broadcast the info from this obs to all other processes
      ! What gets broadcast depends on what kind of inflation is being done
      !>@todo it should also depend on if vertical is being converted.  the last
      !>two values aren't needed unless vertical conversion is happening.
      !>@todo FIXME: this is messy, but should we have 6 different broadcasts,
      !>the three below and three more which omit the 2 localization values?
      !>how much does this cost in time? time this and see.
      whichvert_real = real(whichvert_obs_in_localization_coord, r8)

      if (filter_kind == 9) then

         if(local_varying_ss_inflate) then
            call broadcast_send(map_pe_to_task(ens_handle, owner), obs_prior, &
               orig_obs_prior_mean, orig_obs_prior_var, w, hw, scalar1=obs_qc, &
               scalar2=vertvalue_obs_in_localization_coord, scalar3=whichvert_real)
         else
            call broadcast_send(map_pe_to_task(ens_handle, owner), obs_prior, w, hw, scalar1=obs_qc, &
               scalar2=vertvalue_obs_in_localization_coord, scalar3=whichvert_real)
         endif

      else

         if(local_varying_ss_inflate) then
            call broadcast_send(map_pe_to_task(ens_handle, owner), obs_prior, obs_inc, &
               orig_obs_prior_mean, orig_obs_prior_var, net_a, scalar1=obs_qc, &
               scalar2=vertvalue_obs_in_localization_coord, scalar3=whichvert_real)
         else if(local_single_ss_inflate .or. local_obs_inflate) then
            call broadcast_send(map_pe_to_task(ens_handle, owner), obs_prior, obs_inc, &
              net_a, scalar1=my_inflate, scalar2=my_inflate_sd, scalar3=obs_qc, &
              scalar4=vertvalue_obs_in_localization_coord, scalar5=whichvert_real)
         else
            call broadcast_send(map_pe_to_task(ens_handle, owner), obs_prior, obs_inc, &
              net_a, scalar1=obs_qc, &
              scalar2=vertvalue_obs_in_localization_coord, scalar3=whichvert_real)
         endif

      endif

   ! Next block is done by processes that do NOT own this observation
   !-----------------------------------------------------------------------
   else
      ! I don't store this obs; receive the obs prior and increment from broadcast
      ! Also get qc and inflation information if needed
      ! also a converted vertical coordinate if needed
      !>@todo FIXME see the comment in the broadcast_send() section about
      !>the cost of sending unneeded values

      ! PF needs to broadcast different variables than other filters 
      if (filter_kind == 9) then

         if(local_varying_ss_inflate) then
            call broadcast_recv(map_pe_to_task(ens_handle, owner), obs_prior, &
               orig_obs_prior_mean, orig_obs_prior_var, w, hw, scalar1=obs_qc,  &
               scalar2=vertvalue_obs_in_localization_coord, scalar3=whichvert_real)
         else
            call broadcast_recv(map_pe_to_task(ens_handle, owner), obs_prior, w, hw, scalar1=obs_qc, &
               scalar2=vertvalue_obs_in_localization_coord, scalar3=whichvert_real)
         endif
      else
         if(local_varying_ss_inflate) then
            call broadcast_recv(map_pe_to_task(ens_handle, owner), obs_prior, obs_inc, &
               orig_obs_prior_mean, orig_obs_prior_var, net_a, scalar1=obs_qc, &
               scalar2=vertvalue_obs_in_localization_coord, scalar3=whichvert_real)
         else if(local_single_ss_inflate .or. local_obs_inflate) then
            call broadcast_recv(map_pe_to_task(ens_handle, owner), obs_prior, obs_inc, &
               net_a, scalar1=my_inflate, scalar2=my_inflate_sd, scalar3=obs_qc, &
               scalar4=vertvalue_obs_in_localization_coord, scalar5=whichvert_real)
         else
            call broadcast_recv(map_pe_to_task(ens_handle, owner), obs_prior, obs_inc, &
              net_a, scalar1=obs_qc, &
              scalar2=vertvalue_obs_in_localization_coord, scalar3=whichvert_real)
         endif
      endif
      whichvert_obs_in_localization_coord = nint(whichvert_real)

   endif
   !-----------------------------------------------------------------------

   ! Everybody is doing this section, cycle if qc is bad
   if(nint(obs_qc) /= 0) then
      if (timing(MLOOP)) then
         write(msgstring, '(A32,I7)') 'sequential obs cycl: obs', keys(i)
         call read_timer(t_base(MLOOP), msgstring, elapsed = elapse_array(i))
      endif
      cycle SEQUENTIAL_OBS
   endif
   
   !> all tasks must set the converted vertical values into the 'base' version of this loc
   !> because that's what we pass into the get_close_xxx() routines below.
   if (is_doing_vertical_conversion) &
      call set_vertical(base_obs_loc, vertvalue_obs_in_localization_coord, whichvert_obs_in_localization_coord)
   
   ! Skip for PF
   if (filter_kind /= 9) then

      ! Can compute prior mean and variance of obs for each group just once here
      do group = 1, num_groups
         grp_bot = grp_beg(group)
         grp_top = grp_end(group)
         obs_prior_mean(group) = sum(obs_prior(grp_bot:grp_top)) / grp_size
         obs_prior_var(group) = sum((obs_prior(grp_bot:grp_top) - obs_prior_mean(group))**2) / &
            (grp_size - 1)
         if (obs_prior_var(group) < 0.0_r8) obs_prior_var(group) = 0.0_r8
      end do

   endif

   ! If we are doing adaptive localization then we need to know the number of
   ! other observations that are within the localization radius.  We may need
   ! to shrink it, and so we need to know this before doing get_close() for the
   ! state space (even though the state space increments will be computed and
   ! applied first).

   !******************************************


   if (.not. close_obs_caching) then
      if (timing(GC)) call start_timer(t_base(GC), t_items(GC), t_limit(GC), do_sync=.false.)
      call get_close_obs(gc_obs, base_obs_loc, base_obs_type, &
                         my_obs_loc, my_obs_kind, my_obs_type, &
                         num_close_obs, close_obs_ind, close_obs_dist, ens_handle)
      if (timing(GC)) then
         write(msgstring, '(A32,3I7)') 'gc_ob_NC:nobs,tot,obs# ', num_close_obs, obs_ens_handle%my_num_vars, keys(i)
         call read_timer(t_base(GC), msgstring, t_items(GC), t_limit(GC), do_sync=.false.)
      endif

   else

      if (base_obs_loc == last_base_obs_loc) then
         num_close_obs     = last_num_close_obs
         close_obs_ind(:)  = last_close_obs_ind(:)
         close_obs_dist(:) = last_close_obs_dist(:)
         num_close_obs_cached = num_close_obs_cached + 1
      else
         if (timing(GC)) call start_timer(t_base(GC), t_items(GC), t_limit(GC), do_sync=.false.)
         call get_close_obs(gc_obs, base_obs_loc, base_obs_type, &
                            my_obs_loc, my_obs_kind, my_obs_type, &
                            num_close_obs, close_obs_ind, close_obs_dist, ens_handle)
         if (timing(GC)) then
            write(msgstring, '(A32,3I7)') 'gc_ob_C: nobs,tot,obs# ', num_close_obs, obs_ens_handle%my_num_vars, keys(i)
            call read_timer(t_base(GC), msgstring, t_items(GC), t_limit(GC), do_sync=.false.)
         endif

         last_base_obs_loc      = base_obs_loc
         last_num_close_obs     = num_close_obs
         last_close_obs_ind(:)  = close_obs_ind(:)
         last_close_obs_dist(:) = close_obs_dist(:)
         num_close_obs_calls_made = num_close_obs_calls_made +1
      endif
   endif

   n_close_obs_items(i) = num_close_obs
    !print*, 'base_obs _oc', base_obs_loc, 'rank ', my_task_id()
    !call test_close_obs_dist(close_obs_dist, num_close_obs, i)
    !print*, 'num close ', num_close_obs

   ! set the cutoff default, keep a copy of the original value, and avoid
   ! looking up the cutoff in a list if the incoming obs is an identity ob
   ! (and therefore has a negative kind).  specific types can never be 0;
   ! generic kinds (not used here) start their numbering at 0 instead of 1.
   if (base_obs_type > 0) then
      cutoff_orig = cutoff_list(base_obs_type)
   else
      cutoff_orig = cutoff
   endif

   cutoff_rev = cutoff_orig

   ! For adaptive localization, need number of other obs close to the chosen observation
   if(adaptive_localization_threshold > 0) then

      if (timing(GC)) call start_timer(t_base(GC), t_items(GC), t_limit(GC), do_sync=.false.)

      ! this does a cross-task sum, so all tasks must make this call.
      total_num_close_obs = count_close(num_close_obs, close_obs_ind, my_obs_type, &
                                        close_obs_dist, cutoff_rev*2.0_r8)
      if (timing(GC)) call read_timer(t_base(GC), 'count_close', t_items(GC), t_limit(GC), do_sync=.false.)


      ! Want expected number of close observations to be reduced to some threshold;
      ! accomplish this by cutting the size of the cutoff distance.
      if(total_num_close_obs > adaptive_localization_threshold) then

         cutoff_rev = revised_distance(cutoff_rev*2.0_r8, adaptive_localization_threshold, &
                                       total_num_close_obs, base_obs_loc, &
                                       adaptive_cutoff_floor*2.0_r8) / 2.0_r8

         if ( output_localization_diagnostics ) then

            ! to really know how many obs are left now, you have to
            ! loop over all the obs, again, count how many kinds are
            ! going to be assim, and explicitly check the distance and
            ! see if it's closer than the new cutoff ( times 2 ), and
            ! then do a global sum to get the total.  since this costs,
            ! do it only when diagnostics are requested.

            ! this does a cross-task sum, so all tasks must make this call.
            rev_num_close_obs = count_close(num_close_obs, close_obs_ind, my_obs_type, &
                                              close_obs_dist, cutoff_rev*2.0_r8)


            ! GSR output the new cutoff
            ! Here is what we might want:
            ! time, ob index #, ob location, new cutoff, the assimilate obs count, owner (which process has this ob)
            ! obs_time, obs_val_index, base_obs_loc, cutoff_rev, total_num_close_obs, owner
            ! break up the time into secs and days, and break up the location into lat, lon and height
            ! nsc - the min info here that can't be extracted from the obs key is:
            !  key (obs#), total_num_close_obs (close w/ original cutoff), revised cutoff & new count
            if (my_task_id() == 0) then
               call get_obs_def(observation, obs_def)
               this_obs_time = get_obs_def_time(obs_def)
               call get_time(this_obs_time,secs,days)
               call write_location(-1, base_obs_loc, charstring=base_loc_text)

               write(localization_unit,'(i12,1x,i5,1x,i8,1x,A,2(f14.5,1x,i12))') i, secs, days, &
                     trim(base_loc_text), cutoff_orig, total_num_close_obs, cutoff_rev, rev_num_close_obs
            endif
         endif

      endif

   else if (output_localization_diagnostics) then

      ! if you aren't adapting but you still want to know how many obs are within the
      ! localization radius, set the diag output.  this could be large, use carefully.

      ! this does a cross-task sum, so all tasks must make this call.
      total_num_close_obs = count_close(num_close_obs, close_obs_ind, my_obs_type, &
                                        close_obs_dist, cutoff_rev*2.0_r8)

      if (my_task_id() == 0) then
         call get_obs_def(observation, obs_def)
         this_obs_time = get_obs_def_time(obs_def)
         call get_time(this_obs_time,secs,days)
         call write_location(-1, base_obs_loc, charstring=base_loc_text)

         write(localization_unit,'(i12,1x,i5,1x,i8,1x,A,f14.5,1x,i12)') i, secs, days, &
               trim(base_loc_text), cutoff_rev, total_num_close_obs
      endif
   endif

   ! Now everybody updates their close states
   ! Find state variables on my process that are close to observation being assimilated
   if (.not. close_obs_caching) then
      if (timing(GC)) call start_timer(t_base(GC), t_items(GC), t_limit(GC), do_sync=.false.)
      call get_close_state(gc_state, base_obs_loc, base_obs_type, &
                           my_state_loc, my_state_kind, my_state_indx, &
                           num_close_states, close_state_ind, close_state_dist, ens_handle)
      if (timing(GC)) then
         write(msgstring, '(A32,3I7)') 'gc_st_NC:nsts,tot,obs# ', num_close_states, ens_handle%my_num_vars, keys(i)
         call read_timer(t_base(GC), msgstring, t_items(GC), t_limit(GC), do_sync=.false.)
      endif
   else
      if (base_obs_loc == last_base_states_loc) then
         num_close_states    = last_num_close_states
         close_state_ind(:)  = last_close_state_ind(:)
         close_state_dist(:) = last_close_state_dist(:)
         num_close_states_cached = num_close_states_cached + 1
      else
         if (timing(GC)) call start_timer(t_base(GC), t_items(GC), t_limit(GC), do_sync=.false.)
         call get_close_state(gc_state, base_obs_loc, base_obs_type, &
                              my_state_loc, my_state_kind, my_state_indx, &
                              num_close_states, close_state_ind, close_state_dist, ens_handle)
         if (timing(GC)) then
            write(msgstring, '(A32,3I7)') 'gc_st_C: nsts,tot,obs# ', num_close_states, ens_handle%my_num_vars, keys(i)
            call read_timer(t_base(GC), msgstring, t_items(GC), t_limit(GC), do_sync=.false.)
         endif

         last_base_states_loc     = base_obs_loc
         last_num_close_states    = num_close_states
         last_close_state_ind(:)  = close_state_ind(:)
         last_close_state_dist(:) = close_state_dist(:)
         num_close_states_calls_made = num_close_states_calls_made + 1
      endif
   endif

   n_close_state_items(i) = num_close_states
   !print*, 'num close state', num_close_states
   !call test_close_obs_dist(close_state_dist, num_close_states, i)
   !call test_state_copies(ens_handle, 'beforeupdates')

   if (timing(LG_GRN)) then
      write(msgstring, '(A32,I7)') 'before_state_update: obs', keys(i)
      call read_timer(t_base(LG_GRN), msgstring)
   endif


   if (filter_kind == 9) then

     ! Skip update step if weights at ob location are uniform
     if (1.0_r8 > 0.98_r8*ens_size*sum(w**2) ) then
       cycle SEQUENTIAL_OBS
     end if

     ! Get sampling indices
     call pf_sample(obs_prior, w(1:ens_size), ens_size, indx(1:ens_size))

     ! Redundant part of weight calculation   
     wt = ens_size*hw - 1.0_r8

   endif

   ! Loop through to update each of my state variables that is potentially close
   if (timing(LG_GRN)) call start_timer(t_base(LG_GRN))
   STATE_UPDATE: do j = 1, num_close_states
      state_index = close_state_ind(j)

      ! Skip if regularization reaches threashold value 
      if (filter_kind == 9) then
         if ( beta(state_index) == beta_max ) cycle STATE_UPDATE
      end if

      ! the "any" is an expensive test when you do it for every ob.  don't test
      ! if we know there aren't going to be missing values in the state.
      if ( allow_missing_in_state ) then
         ! Some models can take evasive action if one or more of the ensembles have
         ! a missing value. Generally means 'do nothing' (as opposed to DIE)
         if (any(ens_handle%copies(1:ens_size, state_index) == MISSING_R8)) cycle STATE_UPDATE
      endif

      ! Get the initial values of inflation for this variable if state varying inflation
      if(local_varying_ss_inflate .and. filter_kind /= 9) then
         varying_ss_inflate    = ens_handle%copies(ENS_INF_COPY,    state_index)
         varying_ss_inflate_sd = ens_handle%copies(ENS_INF_SD_COPY, state_index)
      else
         varying_ss_inflate    = 0.0_r8
         varying_ss_inflate_sd = 0.0_r8
      endif

      ! Compute the distance and covariance factor
      cov_factor = comp_cov_factor(close_state_dist(j), cutoff_rev, &
         base_obs_loc, base_obs_type, my_state_loc(state_index), my_state_kind(state_index))

      ! if external impact factors supplied, factor them in here
      ! FIXME: this would execute faster for 0.0 impact factors if
      ! we check for that before calling comp_cov_factor.  but it makes
      ! the logic more complicated - this is simpler if we do it after.
      if (adjust_obs_impact) then
         impact_factor = obs_impact_table(base_obs_type, my_state_kind(state_index))
         cov_factor = cov_factor * impact_factor
      endif

      ! If no weight is indicated, no more to do with this state variable
      if (filter_kind == 9) then ! note other filter kinds will need to cycle
         if ( cov_factor == 0.0_r8 ) cycle STATE_UPDATE
      end if

      ! Update state for PF
      if (filter_kind == 9) then

         ! Update only when prior variance is not zero
         if ( maxval( ens_handle%copies(1:ens_size, state_index)) /= &
              minval( ens_handle%copies(1:ens_size, state_index)) ) then

            ! Take running sum of -log() of weights
            if (cov_factor == 1.0_r8) then
              d = log(ens_size*hw)
            else
              d = wt*cov_factor
              do n = 1,ens_size
                if (abs(d(n)) > 0.1_r8) then
                  d(n) = log( d(n) + 1.0_r8 )
                end if
              end do
            end if
            lw(1:ens_size,state_index) = lw(1:ens_size,state_index) - d
            lw(1:ens_size,state_index) = lw(1:ens_size,state_index) - minval(lw(1:ens_size,state_index))

            ! Get state-space weights
            wo = exp(-lw(1:ens_size,state_index)/beta(state_index))
            wo = wo/sum(wo)

            ! Use weights to calculate posterior mean and variance
            ens_mean = sum( wo * ens_init(1:ens_size,state_index) )
           
            if (sum(ens_handle%copies(1:ens_size,state_index))/ens_size .ne. ens_mean) then

               ens_var = sum( wo * ( ens_init(1:ens_size,state_index) - ens_mean )**2 ) &
                         / ( 1.0_r8 - sum(wo**2) )

               ! Perform sampling from weighted prior or unweighted posterior 
               if (sampling_weighted_prior) then
                 d = ens_init(indx,state_index)
               else
                 d = ens_handle%copies(indx,state_index)
               end if

               ! Combine newly sampled particles with prior particles
               call pf_update(ens_handle%copies(1:ens_size,state_index), ens_mean, ens_var, &
                              increment(1:ens_size), ens_size, cov_factor, d, pf_alpha)

            else
   
               increment(1:ens_size) = 0.0_r8

            end if
  
            else

           increment(1:ens_size) = 0.0_r8

         endif

         ! Set reg_factor to 1
         reg_factor = 1.0_r8

         ! Need correl for ss inflation
         if(local_varying_ss_inflate .and. varying_ss_inflate > 0.0_r8 .and. &
            varying_ss_inflate_sd > 0.0_r8 .and. filter_kind /= 9) then

            call pf_calc_correl(obs_prior(1:ens_size),ens_handle%copies(1:ens_size, state_index), &
                      ens_size, correl(1))

            ! Include localization in correl
            correl(1) = correl(1) * cov_factor
         endif

      else

         ! All other filters

         ! Loop through groups to update the state variable ensemble members
         do group = 1, num_groups
            grp_bot = grp_beg(group)
            grp_top = grp_end(group)
            ! Do update of state, correl only needed for varying ss inflate
            if(local_varying_ss_inflate .and. varying_ss_inflate > 0.0_r8 .and. &
               varying_ss_inflate_sd > 0.0_r8 .and. filter_kind /= 9) then
               call update_from_obs_inc(obs_prior(grp_bot:grp_top), obs_prior_mean(group), &
                  obs_prior_var(group), obs_inc(grp_bot:grp_top), &
                  ens_handle%copies(grp_bot:grp_top, state_index), grp_size, &
                  increment(grp_bot:grp_top), reg_coef(group), net_a(group), correl(group))
            else
               call update_from_obs_inc(obs_prior(grp_bot:grp_top), obs_prior_mean(group), &
                  obs_prior_var(group), obs_inc(grp_bot:grp_top), &
                  ens_handle%copies(grp_bot:grp_top, state_index), grp_size, &
                  increment(grp_bot:grp_top), reg_coef(group), net_a(group))
            endif
         end do
         if (timing(SM_GRN)) call read_timer(t_base(SM_GRN), 'update_from_obs_inc_S', &
                                          t_items(SM_GRN), t_limit(SM_GRN), do_sync=.false.)

         ! Compute an information factor for impact of this observation on this state
         if(num_groups == 1) then
             reg_factor = 1.0_r8
         else
            ! Pass the time along with the index for possible diagnostic output
            ! Compute regression factor for this obs-state pair
            reg_factor = comp_reg_factor(num_groups, reg_coef, obs_time, i, my_state_indx(state_index))
         endif

         ! The final factor is the minimum of group regression factor and localization cov_factor
         reg_factor = min(reg_factor, cov_factor)

      endif ! Filter other than PF

!PAR NEED TO TURN STUFF OFF MORE EFFICEINTLY
      ! If doing full assimilation, update the state variable ensemble with weighted increments
      if(.not. inflate_only) then
         ens_handle%copies(1:ens_size, state_index) = &
            ens_handle%copies(1:ens_size, state_index) + reg_factor * increment
      endif

      ! Compute spatially-varying state space inflation
      if(local_varying_ss_inflate .and. filter_kind /= 9) then
         ! base is the initial inflate value for this state variable
         ss_inflate_base = ens_handle%copies(ENS_SD_COPY, state_index)
         ! Loop through each group to update inflation estimate
         GroupInflate: do group = 1, num_groups
            if(varying_ss_inflate > 0.0_r8 .and. varying_ss_inflate_sd > 0.0_r8) then
               ! Gamma is less than 1 for varying ss, see adaptive inflate module
               gamma = reg_factor * abs(correl(group))
               ! Deflate the inflated variance using the INITIAL state inflate
               ! value (before these obs started gumming it up).
               ens_obs_mean = orig_obs_prior_mean(group)
               ens_obs_var =  orig_obs_prior_var(group)

               ! Remove the impact of inflation to allow efficient single pass with assim.
               if ( abs(gamma) > small ) then
                  ens_var_deflate = ens_obs_var / &
                     (1.0_r8 + gamma*(sqrt(ss_inflate_base) - 1.0_r8))**2
               else
                  ens_var_deflate = ens_obs_var
               endif

               ! If this is inflate only (i.e. posterior) remove impact of this obs.
               if(inflate_only .and. &
                     ens_var_deflate               > small .and. &
                     obs_err_var                   > small .and. &
                     obs_err_var - ens_var_deflate > small ) then
                  r_var  = 1.0_r8 / (1.0_r8 / ens_var_deflate - 1.0_r8 / obs_err_var)
                  r_mean = r_var *(ens_obs_mean / ens_var_deflate - obs(1) / obs_err_var)
               else
                  r_var = ens_var_deflate
                  r_mean = ens_obs_mean
               endif

               ! IS A TABLE LOOKUP POSSIBLE TO ACCELERATE THIS?
               ! Update the inflation values
               if (timing(SM_GRN)) call start_timer(t_base(SM_GRN), t_items(SM_GRN), t_limit(SM_GRN), do_sync=.false.)
               call update_inflation(inflate, varying_ss_inflate, varying_ss_inflate_sd, &
                  r_mean, r_var, grp_size, obs(1), obs_err_var, gamma)
               if (timing(SM_GRN)) call read_timer(t_base(SM_GRN), 'update_inflation_V', &
                                                   t_items(SM_GRN), t_limit(SM_GRN), do_sync=.false.)
            else
               ! if we don't go into the previous if block, make sure these
               ! have good values going out for the block below
               r_mean = orig_obs_prior_mean(group)
               r_var =  orig_obs_prior_var(group)
            endif

            ! Update adaptive values if posterior outlier_ratio test doesn't fail.
            ! Match code in obs_space_diags() in filter.f90
            do_adapt_inf_update = .true.
            if (inflate_only) then
               diff_sd = sqrt(obs_err_var + r_var)
               if (diff_sd > 0.0_r8) then
                  outlier_ratio = abs(obs(1) - r_mean) / diff_sd
                  do_adapt_inf_update = (outlier_ratio <= 3.0_r8)
               endif
            endif
            if (do_adapt_inf_update) then
               ens_handle%copies(ENS_INF_COPY, state_index) = varying_ss_inflate
               ens_handle%copies(ENS_INF_SD_COPY, state_index) = varying_ss_inflate_sd
            endif
         end do GroupInflate
      endif

   end do STATE_UPDATE
   if (timing(LG_GRN)) then
      write(msgstring, '(A32,I7)') 'state_update: obs', keys(i)
      call read_timer(t_base(LG_GRN), msgstring)
   endif

   !call test_state_copies(ens_handle, 'after_state_updates')

   !------------------------------------------------------

   ! Now everybody updates their obs priors (only ones after this one)
   if (timing(LG_GRN)) call start_timer(t_base(LG_GRN))
   OBS_UPDATE: do j = 1, num_close_obs
      obs_index = close_obs_ind(j)

      ! The local PF now iterates multiple times over observations, so it still needs 
      ! to update obs-space priors for measurements that are already assimilated.
      ! For other filter_kinds, this line should be put back in to cycle for obs which 
      ! have already been processed. 
      ! if (my_obs_indx(obs_index) <= i) cycle OBS_UPDATE

      ! Skip if regularization reaches threashold value 
      if (filter_kind == 9) then
         if ( beta_y(obs_index) == beta_max ) cycle OBS_UPDATE
      end if

      ! If the forward observation operator failed, no need to 
      ! update the unassimilated observations 
      if (any(obs_ens_handle%copies(1:ens_size, obs_index) == MISSING_R8)) cycle OBS_UPDATE

      ! Compute the distance and the covar_factor
      cov_factor = comp_cov_factor(close_obs_dist(j), cutoff_rev, &
         base_obs_loc, base_obs_type, my_obs_loc(obs_index), my_obs_kind(obs_index))

      ! if external impact factors supplied, factor them in here
      ! FIXME: this would execute faster for 0.0 impact factors if
      ! we check for that before calling comp_cov_factor.  but it makes
      ! the logic more complicated - this is simpler if we do it after.
      if (adjust_obs_impact) then
         impact_factor = obs_impact_table(base_obs_type, my_obs_kind(obs_index))
         cov_factor = cov_factor * impact_factor
      endif
 
      if (filter_kind == 9) then ! note other filter kinds will also need to cycle
         if ( cov_factor == 0.0_r8 ) cycle OBS_UPDATE
      end if

      ! Update obs prior for PF
      if (filter_kind == 9) then

         ! Update only when prior variance is not zero
         if ( maxval( obs_ens_handle%copies(1:ens_size, obs_index)) /= &
              minval( obs_ens_handle%copies(1:ens_size, obs_index)) ) then

            ! Take running sum of -log() of weights
            if (cov_factor == 1.0_r8) then
              d = log(ens_size*hw)
            else
              d = wt*cov_factor
              do n = 1,ens_size
                if (abs(d(n)) > 0.1_r8) then
                  d(n) = log( d(n) + 1.0_r8 )
                end if
              end do
            end if
            lhw(1:ens_size,obs_index) = lhw(1:ens_size,obs_index) - d
            lhw(1:ens_size,obs_index) = lhw(1:ens_size,obs_index) - minval(lhw(1:ens_size,obs_index))

            ! Get obs-space weights
            wo = exp(-lhw(1:ens_size,obs_index)/beta_y(obs_index))
            wo = wo/sum(wo)

            ! Use weights to calculate posterior mean and variance
            ens_mean = sum( wo * obs_ens_init(1:ens_size,obs_index) )

            if (sum(obs_ens_handle%copies(1:ens_size,obs_index))/ens_size .ne. ens_mean) then

               ens_var = sum( wo * ( obs_ens_init(1:ens_size,obs_index) - ens_mean )**2 ) &
                         / ( 1.0_r8 - sum(wo**2) )

               ! Perform sampling from weighted prior or unweighted posterior 
               if (sampling_weighted_prior) then
                 d = obs_ens_init(indx,obs_index)
               else
                 d = obs_ens_handle%copies(indx, obs_index)
               end if

               ! Combine newly sampled particles with prior particles
               call pf_update(obs_ens_handle%copies(1:ens_size,obs_index), ens_mean, ens_var, &
                              increment(1:ens_size), ens_size, cov_factor, d, pf_alpha)


            else
   
               increment(1:ens_size) = 0.0_r8
   
            end if

   
         else

            increment(1:ens_size) = 0.0_r8

         endif

         ! Set reg_factor to 1
         reg_factor = 1.0_r8

      else

         ! Loop through and update ensemble members in each group
         do group = 1, num_groups
            grp_bot = grp_beg(group)
            grp_top = grp_end(group)
            call update_from_obs_inc(obs_prior(grp_bot:grp_top), obs_prior_mean(group), &
               obs_prior_var(group), obs_inc(grp_bot:grp_top), &
                obs_ens_handle%copies(grp_bot:grp_top, obs_index), grp_size, &
                increment(grp_bot:grp_top), reg_coef(group), net_a(group))
         end do
         if (timing(SM_GRN)) call read_timer(t_base(SM_GRN), 'update_from_obs_inc_O', &
                                             t_items(SM_GRN), t_limit(SM_GRN), do_sync=.false.)

         ! FIXME: could we move the if test for inflate only to here?

         ! Compute an information factor for impact of this observation on this state
         if(num_groups == 1) then
             reg_factor = 1.0_r8
         else
            ! Pass the time along with the index for possible diagnostic output
            ! Compute regression factor for this obs-state pair
            ! Negative indicates that this is an observation index
            reg_factor = comp_reg_factor(num_groups, reg_coef, obs_time, i, -1*my_obs_indx(obs_index))
         endif

         ! Final weight is min of group and localization factors
         reg_factor = min(reg_factor, cov_factor)

      endif ! Filter other than PF

      ! Only update state if indicated (otherwise just getting inflation)
      if(.not. inflate_only) then
         obs_ens_handle%copies(1:ens_size, obs_index) = &
           obs_ens_handle%copies(1:ens_size, obs_index) + reg_factor * increment
      endif

   end do OBS_UPDATE
   if (timing(LG_GRN)) then
      write(msgstring, '(A32,I7)') 'obs_update: obs', keys(i)
      call read_timer(t_base(LG_GRN), msgstring)
   endif

   !call test_state_copies(ens_handle, 'after_obs_updates')

   if (timing(MLOOP)) then
      write(msgstring, '(A32,I7)') 'outer sequential obs loop: obs', keys(i)
      call read_timer(t_base(MLOOP), msgstring, elapsed = elapse_array(i))
   endif

end do SEQUENTIAL_OBS

! Additional corrections to state variables using KDDM
if (filter_kind == 9 .and. pf_kddm > 0 ) then

   if (my_task_id() == 0) then
      call date_and_time( date, time )
      write(msgstring,*) 'KDDM begin time: ',time(1:2),':',time(3:4),':',time(5:6)
      call error_handler(E_MSG,'',msgstring)
   endif

   do i = 1,ens_handle%my_num_vars

     ! Update only when posterior ensemble is outside span of prior ensemble
     if ( ( maxval(ens_handle%copies(1:ens_size, i)) > maxval(ens_init(1:ens_size, i) ) ) .or. & 
          ( minval(ens_handle%copies(1:ens_size, i)) < minval(ens_init(1:ens_size, i) ) ) ) then

       wt = exp(-lw(1:ens_size,i)/beta(i))
       wt = wt/sum(wt)

       call pf_kddm_update(ens_handle%copies(1:ens_size,i), ens_init(1:ens_size,i), &
           wt, ens_size, increment)

       ens_handle%copies(1:ens_size, i) = ens_handle%copies(1:ens_size, i) + increment

     endif

   enddo

   do i = 1,obs_ens_handle%my_num_vars

     ! Update only when posterior ensemble is outside span of prior ensemble
     if ( maxval(obs_ens_handle%copies(1:ens_size, i)) > maxval(obs_ens_init(1:ens_size, i) ) .or. & 
        minval(obs_ens_handle%copies(1:ens_size, i)) < minval(obs_ens_init(1:ens_size, i) ) ) then

       wt = exp(-lhw(1:ens_size,i)/beta_y(i))
       wt = wt/sum(wt)

       call pf_kddm_update(obs_ens_handle%copies(1:ens_size,i), obs_ens_init(1:ens_size,i), &
               wt, ens_size, increment)

       obs_ens_handle%copies(1:ens_size, i) = obs_ens_handle%copies(1:ens_size, i) + increment
  
     endif

   enddo

   if (my_task_id() == 0) then
      call date_and_time( date, time )
      write(msgstring,*) 'KDDM end time:   ',time(1:2),':',time(3:4),':',time(5:6)
      call error_handler(E_MSG,'',msgstring)
   endif

endif ! KDDM

! Perform RTPS step for hybrid. For hybrid PF-EAKF the last iteration switches to EAKF (filter_kind=1)
if (filter_kind == 1 .and. .not. local_varying_ss_inflate) then

   if (my_task_id() == 0) then

      write(msgstring, *) 'Performing RTPS with alpha of',pf_kf_rtps_coeff
      call error_handler(E_MSG,'filter_assim:',msgstring)

   end if

   do i = 1,ens_handle%my_num_vars

      wt = ens_init(1:ens_size, i) - sum(ens_init(1:ens_size, i)) / ens_size
      temp1 = sqrt( sum(wt**2) / (ens_size - 1.0_r8) )
      ens_mean = sum(ens_handle%copies(1:ens_size, i)) / ens_size
      hw = ens_handle%copies(1:ens_size, i) - ens_mean
      temp2 = sqrt( sum(hw**2) / (ens_size - 1.0_r8) )
      if (temp2 > 0.0_r8) then
         ens_handle%copies(1:ens_size, i) = ens_mean + hw*( pf_kf_rtps_coeff*(temp1 - temp2)/temp2 + 1.0_r8)
      end if

   end do

   do i = 1,obs_ens_handle%my_num_vars

      wt = obs_ens_init(1:ens_size, i) - sum(obs_ens_init(1:ens_size, i)) / ens_size
      temp1 = sqrt( sum(wt**2) / (ens_size - 1.0_r8) )
      ens_mean = sum(obs_ens_handle%copies(1:ens_size, i)) / ens_size
      hw = obs_ens_handle%copies(1:ens_size, i) - ens_mean
      temp2 = sqrt( sum(hw**2) / (ens_size - 1.0_r8) )
      if (temp2 > 0.0_r8) then
         obs_ens_handle%copies(1:ens_size, i) = ens_mean + hw*( pf_kf_rtps_coeff*(temp1 - temp2)/temp2 + 1.0_r8)
      end if

   end do


   if (my_task_id() == 0) then

      write(msgstring, *) 'Completed iteration with EAKF'
      call error_handler(E_MSG,'filter_assim:',msgstring)

   end if

   ! USE RTPS FOR EAKF ONLY
   if (iter > 1) then
      exit ITERATIONS
   end if

end if

! Exit PF loop by setting max_res to zero once maxiter is reached
if (iter == maxiter-1) then
  max_res = 0.0_r8
end if

! Exit PF iteration loop when sum of regularization coefficients surpasses threshold
if (filter_kind == 9 .and. max_res == 0.0_r8) then

  if (my_task_id() == 0) then

     write(msgstring, *) 'Number of PF iterations: ',iter
     call error_handler(E_MSG,'filter_assim:',msgstring)

  end if

  ! Exit if EAKF is not needed
  if (.not. pf_enkf_hybrid) then
     exit ITERATIONS
  end if

  ! Switch to EAKF for last iteration 
  filter_kind = 1

end if

end do ITERATIONS

! Switch back to PF if needed
if (filter_kind_orig == 9) then
  filter_kind = 9
end if

! Every pe needs to get the current my_inflate and my_inflate_sd back
if(local_single_ss_inflate) then
   ens_handle%copies(ENS_INF_COPY, :) = my_inflate
   ens_handle%copies(ENS_INF_SD_COPY, :) = my_inflate_sd
end if

! Free up the storage
call destroy_obs(observation)
call get_close_destroy(gc_state)
call get_close_destroy(gc_obs)

! print some stats about the assimilation
! (if interesting, could print exactly which obs # was fastest and slowest)
if (my_task_id() == 0 .and. timing(MLOOP)) then
   write(msgstring, *) 'average assim time: ', sum(elapse_array) / size(elapse_array)
   call error_handler(E_MSG,'filter_assim:',msgstring)

   write(msgstring, *) 'minimum assim time: ', minval(elapse_array)
   call error_handler(E_MSG,'filter_assim:',msgstring)

   write(msgstring, *) 'maximum assim time: ', maxval(elapse_array)
   call error_handler(E_MSG,'filter_assim:',msgstring)
endif

if (timing(MLOOP)) deallocate(elapse_array)

! do some stats - being aware that unless we do a reduce() operation
! this is going to be per-task.  so only print if something interesting
! shows up in the stats?  maybe it would be worth a reduce() call here?

!>@todo FIXME:  
!  we have n_close_obs_items and n_close_state_items for each assimilated
!  observation.  what we really want to know is across the tasks is there
!  a big difference in counts?  so that means communication.  maybe just
!  the largest value?  and the number of 0 values?  and if the largest val
!  is way off compared to the other tasks, warn the user?
!  we don't have space or time to do all the obs * tasks but could we
!  send enough info to make a histogram?  compute N bin counts and then
!  reduce that across all the tasks and have task 0 print out?
! still thinking on this idea.
!   write(msgstring, *) 'max state items per observation: ', maxval(n_close_state_items)
!   call error_handler(E_MSG, 'filter_assim:', msgstring)
! if i come up with something i like, can we use the same idea
! for the threed_sphere locations boxes?

! Assure user we have done something
if (print_trace_details >= 0) then
write(msgstring, '(A,I8,A)') &
   'Processed', obs_ens_handle%num_vars, ' total observations'
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

!GSR close the localization diagnostics file
if(output_localization_diagnostics .and. my_task_id() == 0) then
  call close_file(localization_unit)
end if

! get rid of mpi window
call free_mean_window()

! deallocate space
deallocate(close_obs_dist,      &
           last_close_obs_dist, &
           my_obs_indx,         &
           my_obs_kind,         &
           my_obs_type,         &
           close_obs_ind,       &
           last_close_obs_ind,  &
           vstatus,             &
           my_obs_loc)

deallocate(close_state_dist,      &
           last_close_state_dist, &
           my_state_indx,         &
           close_state_ind,       &
           last_close_state_ind,  &
           my_state_kind,         &
           my_state_loc)

deallocate(n_close_state_items, &
           n_close_obs_items)
! end dealloc

end subroutine filter_assim

!-------------------------------------------------------------

subroutine obs_increment(ens_in, ens_size, obs, obs_var, obs_inc, &
   inflate, my_cov_inflate, my_cov_inflate_sd, net_a)

! Given the ensemble prior for an observation, the observation, and
! the observation error variance, computes increments and adjusts
! observation space inflation values

integer,                     intent(in)    :: ens_size
real(r8),                    intent(in)    :: ens_in(ens_size), obs, obs_var
real(r8),                    intent(out)   :: obs_inc(ens_size)
type(adaptive_inflate_type), intent(inout) :: inflate
real(r8),                    intent(inout) :: my_cov_inflate, my_cov_inflate_sd
real(r8),                    intent(out)   :: net_a

real(r8) :: ens(ens_size), inflate_inc(ens_size)
real(r8) :: prior_mean, prior_var, new_val(ens_size)
integer  :: i, ens_index(ens_size), new_index(ens_size)

real(r8) :: rel_weights(ens_size)

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
   if(filter_kind == 1) then
      call obs_increment_eakf(ens, ens_size, prior_mean, prior_var, &
         obs, obs_var, obs_inc, net_a)
   else if(filter_kind == 2) then
      call obs_increment_enkf(ens, ens_size, prior_var, obs, obs_var, obs_inc)
   else if(filter_kind == 3) then
      call obs_increment_kernel(ens, ens_size, obs, obs_var, obs_inc)
   else if(filter_kind == 4) then
      call obs_increment_particle(ens, ens_size, obs, obs_var, obs_inc)
   else if(filter_kind == 5) then
      call obs_increment_ran_kf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
   else if(filter_kind == 6) then
      call obs_increment_det_kf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
   else if(filter_kind == 7) then
      call obs_increment_boxcar(ens, ens_size, obs, obs_var, obs_inc, rel_weights)
   else if(filter_kind == 8) then
      call obs_increment_rank_histogram(ens, ens_size, prior_var, obs, obs_var, obs_inc)
   else
      call error_handler(E_ERR,'obs_increment', &
              'Illegal value of filter_kind in assim_tools namelist [1-8 OK]', source)
   endif
endif

! Add in the extra increments if doing observation space covariance inflation
if(do_obs_inflate(inflate)) obs_inc = obs_inc + inflate_inc

! To minimize regression errors, may want to sort to minimize increments
! This makes sense for any of the non-deterministic algorithms
! By doing it here, can take care of both standard non-deterministic updates
! plus non-deterministic obs space covariance inflation. This is expensive, so
! don't use it if it's not needed.
if (sort_obs_inc) then
   new_val = ens_in + obs_inc
   ! Sorting to make increments as small as possible
   call index_sort(ens_in, ens_index, ens_size)
   call index_sort(new_val, new_index, ens_size)
   do i = 1, ens_size
      obs_inc(ens_index(i)) = new_val(new_index(i)) - ens_in(ens_index(i))
   end do
endif

! Get the net change in spread if obs space inflation was used
if(do_obs_inflate(inflate)) net_a = net_a * sqrt(my_cov_inflate)

end subroutine obs_increment



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


subroutine obs_increment_ran_kf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
!========================================================================
!
! Forms a random sample of the Gaussian from the update equations.
! This is very close to what a true 'ENSEMBLE' Kalman Filter would
! look like. Note that outliers, multimodality, etc., get tossed.

integer,   intent(in)  :: ens_size
real(r8),  intent(in)  :: prior_mean, prior_var
real(r8),  intent(in)  :: ens(ens_size), obs, obs_var
real(r8),  intent(out) :: obs_inc(ens_size)

real(r8) :: new_mean, var_ratio
real(r8) :: temp_mean, temp_var, new_ens(ens_size), new_var
integer  :: i

var_ratio = obs_var / (prior_var + obs_var)
new_var = var_ratio * prior_var
new_mean  = var_ratio * (prior_mean  + prior_var*obs / obs_var)

! This will reproduce exactly for multiple runs with the same task count,
! but WILL NOT reproduce for a different number of MPI tasks.
! To make it independent of the number of MPI tasks, it would need to
! use the global ensemble number or something else that remains constant
! as the processor count changes.  this is not currently an argument to
! this function and so we are not trying to make it task-count invariant.

! Form a random sample from the updated distribution
! Then adjust the mean (what about adjusting the variance?)!
! Definitely need to sort with this; sort is done in main obs_increment
if(first_inc_ran_call) then
   call init_random_seq(inc_ran_seq, my_task_id() + 1)
   first_inc_ran_call = .false.
endif

do i = 1, ens_size
   new_ens(i) = random_gaussian(inc_ran_seq, new_mean, sqrt(prior_var*var_ratio))
end do

! Adjust the mean of the new ensemble
temp_mean = sum(new_ens) / ens_size
new_ens(:) = new_ens(:) - temp_mean + new_mean

! Compute prior variance and mean from sample
temp_var  = sum((new_ens - new_mean)**2) / (ens_size - 1)
! Adjust the variance, also
new_ens = (new_ens - new_mean) * sqrt(new_var / temp_var) + new_mean

! Get the increments
obs_inc = new_ens - ens

end subroutine obs_increment_ran_kf



subroutine obs_increment_det_kf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
!========================================================================
!
! Does a deterministic ensemble layout for the updated Gaussian.
! Note that all outliers, multimodal behavior, etc. get tossed.

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: prior_mean, prior_var
real(r8), intent(in)  :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: new_mean, var_ratio, temp_var, new_ens(ens_size), new_var
integer :: i

var_ratio = obs_var / (prior_var + obs_var)
new_var = var_ratio * prior_var
new_mean = var_ratio * (prior_mean  + prior_var*obs / obs_var)

! Want a symmetric distribution with kurtosis 3 and variance new_var and mean new_mean
if(ens_size /= 20) then
   write(*, *) 'EXPERIMENTAL version obs_increment_det_kf only works for ens_size 20 now'
   stop
endif

! This has kurtosis of 3.0, verify again from initial uniform
!new_ens(1) = -2.146750_r8
!new_ens(2) = -1.601447_r8
!new_ens(3) = -1.151582_r8
!new_ens(4) = -0.7898650_r8
!new_ens(5) = -0.5086292_r8
!new_ens(6) = -0.2997678_r8
!new_ens(7) = -0.1546035_r8
!new_ens(8) = -6.371084E-02_r8
!new_ens(9) = -1.658448E-02_r8
!new_ens(10) = -9.175255E-04_r8

! This has kurtosis of 3.0, verify again from initial inverse gaussian
!new_ens(1) = -2.188401_r8
!new_ens(2) = -1.502174_r8
!new_ens(3) = -1.094422_r8
!new_ens(4) = -0.8052422_r8
!new_ens(5) = -0.5840152_r8
!new_ens(6) = -0.4084518_r8
!new_ens(7) = -0.2672727_r8
!new_ens(8) = -0.1547534_r8
!new_ens(9) = -6.894587E-02_r8
!new_ens(10) = -1.243549E-02_r8

! This has kurtosis of 2.0, verify again
new_ens(1) = -1.789296_r8
new_ens(2) = -1.523611_r8
new_ens(3) = -1.271505_r8
new_ens(4) = -1.033960_r8
new_ens(5) = -0.8121864_r8
new_ens(6) = -0.6077276_r8
new_ens(7) = -0.4226459_r8
new_ens(8) = -0.2598947_r8
new_ens(9) = -0.1242189_r8
new_ens(10) = -2.539018E-02_r8

! This has kurtosis of 1.7, verify again
!new_ens(1) = -1.648638_r8
!new_ens(2) = -1.459415_r8
!new_ens(3) = -1.272322_r8
!new_ens(4) = -1.087619_r8
!new_ens(5) = -0.9056374_r8
!new_ens(6) = -0.7268229_r8
!new_ens(7) = -0.5518176_r8
!new_ens(8) = -0.3816142_r8
!new_ens(9) = -0.2179997_r8
!new_ens(10) = -6.538583E-02_r8
do i = 11, 20
   new_ens(i) = -1.0_r8 * new_ens(20 + 1 - i)
end do

! Right now, this ensemble has mean 0 and some variance
! Compute prior variance and mean from sample
temp_var  = sum((new_ens)**2) / (ens_size - 1)

! Adjust the variance of this ensemble to match requirements and add in the mean
new_ens = new_ens * sqrt(new_var / temp_var) + new_mean

! Get the increments
obs_inc = new_ens - ens

end subroutine obs_increment_det_kf




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
integer  :: i

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



subroutine pf_regularization(lw, ens_size, Neff, beta, beta_max)
!------------------------------------------------------------------------
!
!  Calculate regularization factors for particle filter: J. Poterjoy Jun. 2019
! 

integer,  intent(in)    :: ens_size
real(r8), intent(in)    :: Neff, beta_max, lw(ens_size)
real(r8), intent(out)   :: beta

real(r8) :: Neff_init, Neff_final, ke, km, ks, ws
real(r8) :: tol, fke, fkm, fks, w(ens_size), beta_base
integer  :: i, tot

! Initial weights and Neff
w = exp(-lw)
ws = sum(w)
w = w/ws
Neff_init = 1.0_r8 / sum( w**2 )

! Inflate if effective ensemble size is smaller than threshold
if ( ( Neff_init < Neff ) .or. ( ws == 0.0_r8 ) ) then

   ! Initial start and end bounds
   ks = 1
!   ke = max(10.0_r8,maxval(lw)**4)
   ke = max(10.0_r8,maxval(lw)**2)

   ! Apply bisection method to solve for k
   tol = 1E-3_r8
   do i = 1,1000
 
      ! Mid point
      km = (ke + ks) / 2.0_r8
  
      ! Evaluate function at end points
      w = exp(-lw/ks)
      if (sum(w) == 0.0_r8) then
         fks = Neff - 1.0_r8
      else
         w = w / sum(w)
         fks = Neff - 1.0_r8 / sum(w**2)
      end if

      w = exp(-lw/ke)
      w = w / sum(w)
      fke = Neff - 1.0_r8 / sum(w**2)

      ! Evaluate function at mid points
      w = exp(-lw/km)
      if (sum(w) == 0.0_r8) then
         fkm = Neff - 1.0_r8
      else
         w = w / sum(w)
         fkm = Neff - 1.0_r8 / sum(w**2)
      end if

      ! Exit critera
      if ( abs(ke-ks) < tol ) exit
 
      ! New end points 
      if ( fkm * fks > 0.0_r8 ) then
        ks = km
      else
        ke = km
      end if

   end do

   beta = km

   ! Underflow errors can still lead to wrong result
   w = exp( -lw/beta)
   w = w / sum(w)
   Neff_final = 1.0_r8 / sum(w**2)
   ! Target Neff is not always obtainable when multiple members have zero weights.
   ! Set beta to max value when min value of 2 is not reached.
   if (Neff_final < Neff - 0.1_r8) then
      beta = beta_max
      write(*,*) 'Warning: setting beta to beta_max'
      write(*,*) 'Neff:',Neff_final
      write(*,*) 'min w:',minval(w)
      write(*,*) 'neff iter:',minval(w)
  end if

  ! Sanity check
  !write(*,*) ' Starting Neff: ',Neff_init,' Target Neff: ',Neff,'New Neff: ',Neff_final,'beta: ',beta,'iterations: ',i

else

   beta = 1.0_r8

end if


end subroutine pf_regularization




subroutine pf_regularization_minw(lw, minwt, ens_size, beta)
!------------------------------------------------------------------------
!
!  Calculate regularization factors for particle filter: J. Poterjoy Jun. 2019
! 

integer,  intent(in)    :: ens_size
real(r8), intent(in)    :: lw(ens_size), minwt
real(r8), intent(out)   :: beta

real(r8) :: minw, ke, km, ks, ws
real(r8) :: tol, fke, fkm, fks, w(ens_size)
integer  :: i

! Initial weights and minw
w = exp(-lw)
ws = sum(w)
w = w/ws

! Inflate if min normalized weight is smaller than threshold
if ( ( minval(w) < minwt ) .or. ( ws == 0.0_r8 ) ) then

   ! Initial start and end bounds
   ks = 1
!   ke = maxval(lw/30.0_r8)
   ke = maxval(lw)
  
   ! Apply bisection method to solve for k
   tol = 0.001_r8

   do i = 1,1000
 
      ! Mid point
      km = (ke + ks) / 2.0_r8
  
      ! Evaluate function at end points
      w = exp(-lw/ks)
      if (sum(w) == 0.0_r8) then
         fks = minwt
      else
         w = w / sum(w)
         fks = minwt - minval(w)
      end if

      w = exp(-lw/ke)
      w = w / sum(w)
      fke = minwt - minval(w)

      ! Evaluate function at mid points
      w = exp(-lw/km)
      w = w / sum(w)
      fkm = minwt - minval(w)

      ! Exit critera
      if ( abs(ke-ks)/2.0_r8 < tol ) exit
 
      ! New end points 
      if ( fkm * fks > 0.0_r8 ) then
        ks = km
      else
        ke = km
      end if

   end do

   beta = km

   w = exp(-lw/beta)
   w = w / sum(w)

  ! Sanity check
  !write(*,*) ' Target min w: ',minwt,'New min w: ',minval(w),'beta: ',beta,'iterations: ',i

else

   beta = 1.0_r8

end if


end subroutine pf_regularization_minw



subroutine pf_calc_correl(obs, state, ens_size, correl)
               
!========================================================================

! Compute correl for adaptive inflation

integer,            intent(in)    :: ens_size
real(r8),           intent(in)    :: obs(ens_size)
real(r8),           intent(in)    :: state(ens_size)
real(r8),           intent(out)   :: correl

real(r8) :: obs_state_cov, intermed, state_mean, state_var, obs_prior_mean, obs_prior_var

obs_prior_mean = sum(obs(1:ens_size)) / ens_size
obs_prior_var  = sum((obs(1:ens_size) - obs_prior_mean)**2) / (ens_size - 1)

state_mean = sum(state) / ens_size
obs_state_cov = sum( (state - state_mean) * (obs - obs_prior_mean) ) / (ens_size - 1)

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

end subroutine pf_calc_correl



subroutine pf_sample(ens, w, ens_size, indx2)
!------------------------------------------------------------------------
!
!  Perform sampling step of particle filter: J. Poterjoy Nov. 2014
! 

integer,  intent(in)    :: ens_size
real(r8), intent(in)    :: w(ens_size), ens(ens_size)
integer,  intent(out)   :: indx2(ens_size)

real(r8) :: cw(0:ens_size), base, frac, dum
integer  :: i, j, indx0(ens_size), indx1(ens_size), m, ind(ens_size)

! Find sorting indices and sort weights
call index_sort(ens, ind, ens_size)

! Perform deterministic resampling
cw(0) = 0.0_r8
do i = 1, ens_size
   cw(i) = cw(i - 1) + w(ind(i))
end do

! Divide interval into ens_size parts and choose new particles
! based on the interval they accumulate in

base = 1.0_r8 / ens_size / 2.0_r8

j = 1
do i = 1, ens_size

   frac = base + (i - 1.0_r8) / ens_size

   ! Search in the cumulative range to see where frac falls
   m = 0
   do while (m == 0)
      if(cw(j - 1) < frac .and. frac <= cw(j)) then
         indx1(i) = j
         m = 1
      else
         j = j + 1
      end if
   end do

end do

! Unsort indices
indx1 = ind(indx1)
indx0 = indx1

! If a particle is removed, it is replaced by a duplicated
! particle. This is accomplished by looping through indx1
! and flagging replicated indices with a zero, and
! indicating their location in indx2

! Locate the removed indices in indx1
do i = 1, ens_size

   ! Locate first occurance of index i in indx1
   m = minloc(indx1, 1, mask=indx1.eq.i)

   if ( m == 0 ) then
      ! If i is not in indx1, flag the index with a zero in indx2
      indx2(i) = 0
   else
      ! If i is in indx1, indicate value in indx2
      indx2(i) = i
      ! Flag value in indx1 with a zero to show it was removed
      indx1(m) = 0
   endif

end do

! TEMP: Uncomment/comment this loop when commenting/uncommenting 
!       the next one
! Replace the removed indices with duplicated ones
do i = 1, ens_size
  if (indx2(i) == 0) then
     do m = 1, ens_size
       if (indx1(m) /= 0) exit
     end do
     indx2(i) = indx1(m)
     indx1(m) = 0
  endif
end do

!! TEMP: Maximize covariance between sampled and removed particles
!do while (sum(indx1) > 0)
!  i = minloc(ens(indx0), 1, mask=indx1.gt.0)
!  m = minloc(ens, 1, mask=indx2.eq.0)
!  indx2(m) = indx1(i)
!  indx1(i) = 0
!end do



! FOR TESTING PURPOSE
!do i = 1,ens_size
!  if ( ens(i) .ne. ens(indx2(i)) ) then
!    write(*,*) 'replacing',ens(i),'with',ens(indx2(i))
!  end if
!end do


end subroutine pf_sample



subroutine pf_update(ens, ens_mean, ens_var, incr, ens_size, loc, ens_s, pf_alpha)
!------------------------------------------------------------------------
!
!  Perform update of particles from weights: J. Poterjoy Nov. 2014
! 

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), loc, pf_alpha, ens_s(ens_size)
real(r8), intent(in)  :: ens_mean, ens_var
real(r8), intent(out) :: incr(ens_size)

real(r8) :: r1, r2, c, c2, em, rho, T1, T2, T3, m1, m2, v1, v2, v3
real(r8) :: ens_post(ens_size), pf_alpha2, smin, smax, alpha
integer  :: i, k

! Calculate weights for updating
c = (1.0_r8-loc)/loc

! r1 and r2 determine coefficients for updating ensemble
v1 = 0.0_r8
v2 = 0.0_r8
v3 = 0.0_r8
do i = 1, ens_size
  v1 = v1 + ( ens_s(i) - ens_mean )**2
  v2 = v2 + ( ens(i) - ens_mean )**2
  v3 = v3 + ( ens(i) - ens_mean )*( ens_s(i) - ens_mean )
end do

c2 = c*c
r1 = v1 + v2*c2 + 2.0_r8*v3*c
r2 = c2/r1

! The coeffiecient pf_alpha reduces part of the update to maintain particle diversity 
! near observation. While pf_alpha is specified, pf_alpha2 is derived to maintain the
! correct amount of spread

alpha = pf_alpha

r1 = alpha*sqrt((ens_size-1.0_r8)*ens_var/r1)
r2 = sqrt((ens_size-1.0_r8)*ens_var*r2)

m1 = sum(ens_s - ens_mean)/ens_size
m2 = sum(ens - ens_mean)/ens_size
v1 = v1 - ens_size*m1**2
v2 = v2 - ens_size*m2**2
v3 = v3 - ens_size*m1*m2
T1 = v2
T2 = 2.0_r8*( r1*v3 + r2*v2 )
T3 = v1*r1**2 + v2*r2**2 + 2.0_r8*v3*r1*r2 - (ens_size-1.0_r8)*ens_var

pf_alpha2 = ( - T2 + sqrt( T2**2 - 4.0_r8*T1*T3 ) ) / (2.0_r8*T1)

r2 = r2 + pf_alpha2

! Update ensemble using mix of prior members and sampled members
do i = 1, ens_size
  ens_post(i) = ens_mean + r1*( ens_s(i) - ens_mean ) + r2*( ens(i) - ens_mean )
end do

! Adjust posterior mean and variance to correct for sampling errors
em = sum( ens_post ) / ens_size
ens_post = ens_mean + (ens_post - em)
incr = ens_post - ens

end subroutine pf_update


subroutine pf_kddm_update(ens1, ens2, w, ens_size, incr)

!------------------------------------------------------------------------
!
!  Apply kernel density distribution mapping method proposed by Seth McGinnis
!  to map a sample of particles into posterior particles:  J. Poterjoy Jan. 2015
! 

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens1(ens_size), ens2(ens_size)
real(r8), intent(inout) :: w(ens_size)
real(r8), intent(out) :: incr(ens_size)
integer,  parameter   :: npoints = 3000
integer               :: i, m, ind(ens_size)
real(r8)              :: xd(npoints), cda(npoints), qf(ens_size), x(ens_size)
real(r8)              :: w2, w1, r(ens_size)

! Note: The specified npoints should depend on range and bandwidth of domain

! Use kernels to approximate quantiles and posterior cdf
call pf_get_q_cda(ens1,ens2,ens_size,npoints,w,xd,qf,cda)

! Correct prior quantiles that land outside span of ensemble
if ( minval(qf) < minval(cda) ) then

   r(1) = 1.0_r8
   r(ens_size) = 0.0_r8
   do i = 2,ens_size-1
     r(i) = r(i-1) - 1.0_r8 / (ens_size - 1.0_r8)
   end do

   ! Sorting indices for quantiles
   call index_sort(ens1, ind, ens_size)

   ! Perform correction
   qf(ind) = qf(ind) + r*( minval(cda) - minval(qf) )

end if

if ( maxval(qf) > maxval(cda) ) then

   r(1) = 0.0_r8
   r(ens_size) = 1.0_r8
   do i = 2,ens_size-1
     r(i) = r(i-1) + 1.0_r8 / (ens_size - 1.0_r8)
   end do

   ! Sorting indices for quantiles
   call index_sort(ens1, ind, ens_size)

   ! Perform correction
   qf(ind) = qf(ind) + r*( maxval(cda) - maxval(qf) )

end if

! Invert posterior cdf to find values at prior quantiles
do i = 1,ens_size

   if ( qf(i) >= maxval(cda) ) then 
      x(i) = maxval(xd)
   else if ( qf(i) <= minval(cda) ) then 
      x(i) = minval(xd)
   else

      m = minloc(cda, 1, mask=cda.gt.qf(i))
  
      if ( (qf(i) == cda(m)) .or. (m == 1) ) then
         x(i) = xd(m)
      else

         if (cda(m) > qf(i)) m = m - 1

         if ( cda(m+1) - cda(m) < 1E-20_r8 ) then
            w1 = ( cda(m+1) - qf(i) ) / ( cda(m+1) - cda(m) )
            w2 = ( qf(i) - cda(m) ) / ( cda(m+1) - cda(m) )
            x(i) = w1 * xd(m) + w2 * xd(m+1)
         else
            x(i) = xd(m)
         end if

      end if

   end if

end do

incr = x - ens1

end subroutine pf_kddm_update



subroutine pf_get_q_cda(ens1,ens2,ens_size,npoints,w,x,qf,cda)
!------------------------------------------------------------------------
!
! Gaussian kernel density estimation:  J. Poterjoy Jan. 2015
! 
! This subroutine returns prior quantiles and posterior cdf
! estimated using Gaussian kernels.

integer,            intent(in)    :: ens_size, npoints
real(r8),           intent(in)    :: ens1(ens_size), ens2(ens_size), w(ens_size)
real(r8),           intent(out)   :: x(npoints), qf(ens_size), cda(npoints)
integer                           :: i, m
real(r8)                          :: bw, xmin, xmax, v2, range

! Bandwidth is set to sample standard deviation
!v2 = sum(w*ens1)
!v2 = sum(w*( ens1 - v2 )**2 ) / (1.0_r8 - sum(w**2))
!v2 =  sum(ens2)/ens_size
!v2 = sum( (ens2 - v2)**2 ) /(ens_size - 1.0_r8)
!bw = sqrt(v2)/10.0_r8

v2 = sum(ens2)/ens_size
v2 = sum( ( ens2 - v2 )**2 ) / (ens_size - 1.0_r8)
bw = sqrt(v2)/4.0_r8

! Domain for calculating posterior cdf
!xmin = min(minval(ens1),minval(ens2)) - 2*bw
!xmax = max(maxval(ens1),maxval(ens2)) + 2*bw
xmin = minval(ens2) - 2*bw
xmax = maxval(ens2) + 2*bw
range = xmax-xmin
do i=1,npoints
  x(i) = xmin + (i-1.0_r8)*range/(npoints-1.0_r8)
end do

! Estimate quantiles and cdfs by taking sum over Gaussian cdfs
qf  = 0.0_r8
cda = 0.0_r8
do i = 1,ens_size

   ! Prior quantiles
   qf = qf + ( 1.0_r8 + erf( (ens1 - ens1(i) )/(sqrt(2.0_r8)*bw) ) )/(2.0_r8*ens_size)

   ! Posterior cdf
   cda = cda + w(i) * ( 1.0_r8 + erf( (x - ens2(i) )/(sqrt(2.0_r8)*bw) ) )/2.0_r8

end do



!write(*,*) minval(cda),maxval(cda)

end subroutine pf_get_q_cda


subroutine update_from_obs_inc(obs, obs_prior_mean, obs_prior_var, obs_inc, &
               state, ens_size, state_inc, reg_coef, net_a, correl_out)
!========================================================================

! Does linear regression of a state variable onto an observation and
! computes state variable increments from observation increments

integer,            intent(in)    :: ens_size
real(r8),           intent(in)    :: obs(ens_size), obs_inc(ens_size)
real(r8),           intent(in)    :: obs_prior_mean, obs_prior_var
real(r8),           intent(in)    :: state(ens_size)
real(r8),           intent(out)   :: state_inc(ens_size), reg_coef
real(r8),           intent(inout) :: net_a
real(r8), optional, intent(inout) :: correl_out

real(r8) :: obs_state_cov, intermed
real(r8) :: restoration_inc(ens_size), state_mean, state_var, correl
real(r8) :: factor, exp_true_correl, mean_factor


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

!
! FIXME: craig schwartz has a degenerate case involving externally computed
! forward operators in which the obs prior variance is in fact exactly 0.
! adding this test allowed him to continue to  use spread restoration
! without numerical problems.  we don't know if this is sufficient;
! for now we'll leave the original code but it needs to be revisited.
!
! Spread restoration algorithm option.
!if(spread_restoration .and. obs_prior_var > 0.0_r8) then
!

! Spread restoration algorithm option.
if(spread_restoration) then
   ! Don't use this to reduce spread at present (should revisit this line)
   if(net_a > 1.0_r8) net_a = 1.0_r8

   ! Default restoration increment is 0.0
   restoration_inc = 0.0_r8

   ! Compute the factor by which to inflate
   ! These come from correl_error.f90 in system_simulation and the files ens??_pairs and
   ! ens_pairs_0.5 in work under system_simulation. Assume a linear reduction from 1
   ! as a function of the net_a. Assume that the slope of this reduction is a function of
   ! the reciprocal of the ensemble_size (slope = 0.80 / ens_size). These are empirical
   ! for now. See also README in spread_restoration_paper documentation.
   !!!factor = 1.0_r8 / (1.0_r8 + (net_a - 1.0_r8) * (0.8_r8 / ens_size)) - 1.0_r8
   factor = 1.0_r8 / (1.0_r8 + (net_a - 1.0_r8) / (-2.4711_r8 + 1.6386_r8 * ens_size)) - 1.0_r8
   !!!factor = 1.0_r8 / (1.0_r8 + (net_a**2 - 1.0_r8) * (-0.0111_r8 + .8585_r8 / ens_size)) - 1.0_r8

   ! Variance restoration
   state_mean = sum(state) / ens_size
   restoration_inc = factor * (state - state_mean)
   state_inc = state_inc + restoration_inc
endif

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



subroutine obs_increment_boxcar(ens, ens_size, obs, obs_var, obs_inc, rel_weight)
!------------------------------------------------------------------------
!
! An observation space update that uses a set of boxcar kernels plus two
! half-gaussians on the wings to represent the prior distribution. If N is
! the ensemble size, 1/(N+1) of the mass is placed between each ensemble
! member. This is reminiscent of the ranked historgram approach for
! evaluating ensembles. The prior distribution on the wings is
! represented by a half gaussian with mean being the outermost ensemble
! member (left or right) and variance being somewhat arbitrarily chosen
! as half the total ensemble sample variance. A particle
! filter like algorithm is then used for the update. The weight associated
! with each prior ensemble member is computed by evaluating the likelihood.
! For the interior, the domain for each boxcar is divided in half and each
! half is associated with the nearest ensemble member. The updated mass in
! each half box is the product of the prior mass and the ensemble weight.
! In the wings, the observation likelihood gaussian is convolved with the
! prior gaussian to get an updated weighted gaussian that is assumed to
! represent the posterior outside of the outermost ensemble members. The
! updated ensemble members are chosen so that 1/(N+1) of the updated
! mass is between each member and also on the left and right wings. This
! algorithm is able to deal well with outliers, bimodality and other
! non-gaussian behavior in observation space. It could also be modified to
! deal with non-gaussian likelihoods in the future.

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
real(r8), intent(out) :: rel_weight(ens_size)

integer  :: i, e_ind(ens_size), lowest_box, j
real(r8) :: sx, prior_mean, prior_var, prior_var_d2
real(r8) :: var_ratio, new_var, new_sd, umass, left_weight, right_weight
real(r8) :: mass(2*ens_size), weight(ens_size), cumul_mass(0:2*ens_size)
real(r8) :: new_mean_left, new_mean_right, prod_weight_left, prod_weight_right
real(r8) :: new_ens(ens_size), mass_sum, const_term
real(r8) :: x(1:2*ens_size - 1), sort_inc(ens_size)

! The factor a is not defined for this filter for now (could it be???)

! The relative weights could be used for a multi-dimensional particle-type
! update using update_ens_from_weights. There are algorithmic challenges
! with outliers so this is not currently a supported option. For now,
! rel_weight is simply set to 0 and is unused elsewhere.
rel_weight = 0.0_r8

! Do an index sort of the ensemble members; Need sorted ensemble
call index_sort(ens, e_ind, ens_size)

! Prior distribution is boxcar in the central bins with 1/(n+1) density
! in each intermediate bin. BUT, distribution on the wings is a normal with
! 1/(n + 1) of the mass on each side.

! Begin by computing a weight for each of the prior ensemble membersA
! This is just evaluating the gaussian likelihood
const_term = 1.0_r8 / (sqrt(2.0_r8 * PI) * sqrt(obs_var))
do i = 1, ens_size
   weight(i) = const_term * exp(-1.0_r8 * (ens(i) - obs)**2 / (2.0_r8 * obs_var))
end do

! Compute the points that bound all the updated mass boxes; start with ensemble
do i = 1, ens_size
   x(2*i - 1) = ens(e_ind(i))
end do
! Compute the mid-point interior boundaries; these are halfway between ensembles
do i = 2, 2*ens_size - 2, 2
   x(i) = (x(i - 1) + x(i + 1)) / 2.0_r8
end do

! Compute the s.d. of the ensemble for getting the gaussian wings
sx         = sum(ens)
prior_mean = sx / ens_size
prior_var  = sum((ens - prior_mean)**2) / (ens_size - 1)

! Need to normalize the wings so they have 1/(ens_size + 1) mass outside
! Since 1/2 of a normal is outside, need to multiply by 2 / (ens_size + 1)

! Need some sort of width for the boundary kernel, try 1/2 the VAR for now
prior_var_d2 = prior_var / 2.0_r8

! Compute the product of the obs error gaussian with the prior gaussian (EAKF)
! Left wing first
var_ratio = obs_var / (prior_var_d2 + obs_var)
new_var = var_ratio * prior_var_d2
new_sd = sqrt(new_var)
new_mean_left  = var_ratio * (ens(e_ind(1))  + prior_var_d2*obs / obs_var)
new_mean_right  = var_ratio * (ens(e_ind(ens_size))  + prior_var_d2*obs / obs_var)
! REMEMBER, this product has an associated weight which must be taken into account
! See Anderson and Anderson for this weight term (or tutorial kernel filter)
prod_weight_left =  2.71828_r8 ** (-0.5_r8 * (ens(e_ind(1))**2 / prior_var_d2 + &
      obs**2 / obs_var - new_mean_left**2 / new_var)) / sqrt(2.0_r8 * PI)

prod_weight_right =  2.71828_r8 ** (-0.5_r8 * (ens(e_ind(ens_size))**2 / prior_var_d2 + &
      obs**2 / obs_var - new_mean_right**2 / new_var)) / sqrt(2.0_r8 * PI)

! Split into 2*ens_size domains; mass in each is computed
! Start by computing mass in the outermost (gaussian) regions
mass(1) = norm_cdf(ens(e_ind(1)), new_mean_left, new_sd) * &
   prod_weight_left * (2.0_r8 / (ens_size + 1.0_r8))
mass(2*ens_size) = (1.0_r8 - norm_cdf(ens(e_ind(ens_size)), new_mean_right, &
   new_sd)) * prod_weight_right * (2.0_r8 / (ens_size + 1.0_r8))

! Compute mass in the inner half boxes that have ensemble point on the left
do i = 2, 2*ens_size - 2, 2
   mass(i) = (1.0_r8 / (2.0_r8 * (ens_size + 1.0_r8))) * weight(e_ind(i/2))
end do

! Now right inner half boxes
do i = 3, 2*ens_size - 1, 2
   mass(i) = (1.0_r8 / (2.0_r8 * (ens_size + 1.0_r8))) * weight(e_ind(i/2 + 1))
end do

! Now normalize the mass in the different bins
mass_sum = sum(mass)
mass = mass / mass_sum

! Find cumulative mass at each box boundary and middle boundary
cumul_mass(0) = 0.0_r8
do i = 1, 2*ens_size
   cumul_mass(i) = cumul_mass(i - 1) + mass(i)
end do

! Get resampled ensemble, Need 1/(ens_size + 1) between each
umass = 1.0_r8 / (ens_size + 1.0_r8)

! Begin search at bottom of lowest box, but then update for efficiency
lowest_box = 1

! Find each new ensemble members location
do i = 1, ens_size
   ! If it's in the inner or outer range have to use normal
   if(umass < cumul_mass(1)) then
      ! In the first normal box
      left_weight = (1.0_r8 / mass_sum) * prod_weight_left * (2.0_r8 / (ens_size + 1.0_r8))
      call weighted_norm_inv(left_weight, new_mean_left, new_sd, umass, new_ens(i))
   else if(umass > cumul_mass(2*ens_size - 1)) then
      ! In the last normal box; Come in from the outside
      right_weight = (1.0_r8 / mass_sum) * prod_weight_right * (2.0_r8 / (ens_size + 1.0_r8))
      call weighted_norm_inv(right_weight, new_mean_right, new_sd, 1.0_r8 - umass, new_ens(i))
      new_ens(i) = new_mean_right + (new_mean_right - new_ens(i))
   else
      ! In one of the inner uniform boxes.
      FIND_BOX:do j = lowest_box, 2 * ens_size - 2
         ! Find the box that this mass is in
         if(umass >= cumul_mass(j) .and. umass <= cumul_mass(j + 1)) then
            new_ens(i) = x(j) + ((umass - cumul_mass(j)) / (cumul_mass(j+1) - cumul_mass(j))) * &
               (x(j + 1) - x(j))
            ! Don't need to search lower boxes again
            lowest_box = j
            exit FIND_BOX
         end if
      end do FIND_BOX
   endif
   ! Want equally partitioned mass in update with exception that outermost boxes have half
   umass = umass + 1.0_r8 / (ens_size + 1.0_r8)
end do

! Can now compute sorted increments
do i = 1, ens_size
   sort_inc(i) = new_ens(i) - ens(e_ind(i))
end do

! Now, need to convert to increments for unsorted
do i = 1, ens_size
   obs_inc(e_ind(i)) = sort_inc(i)
end do

end subroutine obs_increment_boxcar



subroutine obs_increment_rank_histogram(ens, ens_size, prior_var, &
   obs, obs_var, obs_inc)
!------------------------------------------------------------------------
!
! Revised 14 November 2008
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

! This code is still under development. Please contact Jeff Anderson at
! jla@ucar.edu if you are interested in trying it.

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), prior_var, obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

integer  :: i, e_ind(ens_size), lowest_box, j
real(r8) :: prior_sd, var_ratio, umass, left_amp, right_amp
real(r8) :: left_sd, left_var, right_sd, right_var, left_mean, right_mean
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

do i = 1, ens_size
   ! The boundaries of the interior bins are just the sorted ensemble members
   x(i) = ens(e_ind(i))
   ! Compute likelihood for each ensemble member; just evaluate the gaussian
   ! No need to compute the constant term since relative likelihood is what matters
   like(i) = exp(-1.0_r8 * (x(i) - obs)**2 / (2.0_r8 * obs_var))
end do

! Prior distribution is boxcar in the central bins with 1/(n+1) density
! in each intermediate bin. BUT, distribution on the tails is a normal with
! 1/(n + 1) of the mass on each side.

! Can now compute the mean likelihood density in each interior bin
do i = 2, ens_size
   like_dense(i) = ((like(i - 1) + like(i)) / 2.0_r8)
end do

! Compute the s.d. of the ensemble for getting the gaussian tails
prior_sd = sqrt(prior_var)

! For unit normal, find distance from mean to where cdf is 1/(n+1)
! Lots of this can be done once in first call and then saved
call weighted_norm_inv(1.0_r8, 0.0_r8, 1.0_r8, &
   1.0_r8 / (ens_size + 1.0_r8), dist_for_unit_sd)
dist_for_unit_sd = -1.0_r8 * dist_for_unit_sd

! Have variance of tails just be sample prior variance
! Mean is adjusted so that 1/(n+1) is outside
left_mean = x(1) + dist_for_unit_sd * prior_sd
left_var = prior_var
left_sd = prior_sd
! Same for right tail
right_mean = x(ens_size) - dist_for_unit_sd * prior_sd
right_var = prior_var
right_sd = prior_sd

if(gaussian_likelihood_tails) then
   !*************** Block to do Gaussian-Gaussian on tail **************
   ! Compute the product of the obs likelihood gaussian with the priors
   ! Left tail gaussian first
   var_ratio = obs_var / (left_var + obs_var)
   new_var_left = var_ratio * left_var
   new_sd_left = sqrt(new_var_left)
   new_mean_left  = var_ratio * (left_mean  + left_var*obs / obs_var)
   ! REMEMBER, this product has an associated weight which must be taken into account
   ! See Anderson and Anderson for this weight term (or tutorial kernel filter)
   ! NOTE: The constant term has been left off the likelihood so we don't have
   ! to divide by sqrt(2 PI) in this expression
   prod_weight_left =  exp(-0.5_r8 * (left_mean**2 / left_var + &
         obs**2 / obs_var - new_mean_left**2 / new_var_left)) / &
         sqrt(left_var + obs_var)
   ! Determine how much mass is in the updated tails by computing gaussian cdf
   mass(1) = norm_cdf(x(1), new_mean_left, new_sd_left) * prod_weight_left

   ! Same for the right tail
   var_ratio = obs_var / (right_var + obs_var)
   new_var_right = var_ratio * right_var
   new_sd_right = sqrt(new_var_right)
   new_mean_right  = var_ratio * (right_mean  + right_var*obs / obs_var)
   ! NOTE: The constant term has been left off the likelihood so we don't have
   ! to divide by sqrt(2 PI) in this expression
   prod_weight_right =  exp(-0.5_r8 * (right_mean**2 / right_var + &
         obs**2 / obs_var - new_mean_right**2 / new_var_right)) / &
         sqrt(right_var + obs_var)
   ! Determine how much mass is in the updated tails by computing gaussian cdf
   mass(ens_size + 1) = (1.0_r8 - norm_cdf(x(ens_size), new_mean_right, &
      new_sd_right)) * prod_weight_right
   !************ End Block to do Gaussian-Gaussian on tail **************
else
   !*************** Block to do flat tail for likelihood ****************
   ! Flat tails: THIS REMOVES ASSUMPTIONS ABOUT LIKELIHOOD AND CUTS COST
   new_var_left = left_var
   new_sd_left = left_sd
   new_mean_left = left_mean
   prod_weight_left = like(1)
   mass(1) = like(1) / (ens_size + 1.0_r8)

   ! Same for right tail
   new_var_right = right_var
   new_sd_right = right_sd
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
      call weighted_norm_inv(left_amp, new_mean_left, new_sd_left, &
         umass, new_ens(i))
   else if(umass > cumul_mass(ens_size)) then
      ! It's in the right tail
      ! Get position of x in weighted gaussian where the cdf has value umass
      call weighted_norm_inv(right_amp, new_mean_right, new_sd_right, &
         1.0_r8 - umass, new_ens(i))
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
total_mass_left = norm_cdf(ens(e_ind(1)), prior_mean, prior_sd)
total_mass_right = 1.0_r8 - norm_cdf(ens(e_ind(ens_size)), prior_mean, prior_sd)

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
      call weighted_norm_inv(alpha(1), prior_mean, prior_sd, mass, new_ens(i))
   else if(mass > cumul_mass(2*ens_size - 1)) then
      ! In the last normal box; Come in from the outside
      call weighted_norm_inv(alpha(2), prior_mean, prior_sd, 1.0_r8 - mass, new_ens(i))
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


!------------------------------------------------------------------------

function norm_cdf(x_in, mean, sd)

! Approximate cumulative distribution function for normal
! with mean and sd evaluated at point x_in
! Only works for x>= 0.

real(r8)             :: norm_cdf
real(r8), intent(in) :: x_in, mean, sd

real(digits12) :: x, p, b1, b2, b3, b4, b5, t, density, nx

! Convert to a standard normal
nx = (x_in - mean) / sd

x = abs(nx)


! Use formula from Abramowitz and Stegun to approximate
p = 0.2316419_digits12
b1 = 0.319381530_digits12
b2 = -0.356563782_digits12
b3 = 1.781477937_digits12
b4 = -1.821255978_digits12
b5 = 1.330274429_digits12

t = 1.0_digits12 / (1.0_digits12 + p * x)

density = (1.0_digits12 / sqrt(2.0_digits12 * PI)) * exp(-x*x / 2.0_digits12)

norm_cdf = 1.0_digits12 - density * &
   ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t

if(nx < 0.0_digits12) norm_cdf = 1.0_digits12 - norm_cdf

!write(*, *) 'cdf is ', norm_cdf

end function norm_cdf


!------------------------------------------------------------------------

subroutine weighted_norm_inv(alpha, mean, sd, p, x)

! Find the value of x for which the cdf of a N(mean, sd) multiplied times
! alpha has value p.

real(r8), intent(in)  :: alpha, mean, sd, p
real(r8), intent(out) :: x

real(r8) :: np

! Can search in a standard normal, then multiply by sd at end and add mean
! Divide p by alpha to get the right place for weighted normal
np = p / alpha

! Find spot in standard normal
call norm_inv(np, x)

! Add in the mean and normalize by sd
x = mean + x * sd

end subroutine weighted_norm_inv


!------------------------------------------------------------------------

subroutine norm_inv(p, x)

real(r8), intent(in)  :: p
real(r8), intent(out) :: x

! normal inverse
! translate from http://home.online.no/~pjacklam/notes/invnorm
! a routine written by john herrero

real(r8) :: p_low,p_high
real(r8) :: a1,a2,a3,a4,a5,a6
real(r8) :: b1,b2,b3,b4,b5
real(r8) :: c1,c2,c3,c4,c5,c6
real(r8) :: d1,d2,d3,d4
real(r8) :: q,r
a1 = -39.69683028665376_digits12
a2 =  220.9460984245205_digits12
a3 = -275.9285104469687_digits12
a4 =  138.357751867269_digits12
a5 = -30.66479806614716_digits12
a6 =  2.506628277459239_digits12
b1 = -54.4760987982241_digits12
b2 =  161.5858368580409_digits12
b3 = -155.6989798598866_digits12
b4 =  66.80131188771972_digits12
b5 = -13.28068155288572_digits12
c1 = -0.007784894002430293_digits12
c2 = -0.3223964580411365_digits12
c3 = -2.400758277161838_digits12
c4 = -2.549732539343734_digits12
c5 =  4.374664141464968_digits12
c6 =  2.938163982698783_digits12
d1 =  0.007784695709041462_digits12
d2 =  0.3224671290700398_digits12
d3 =  2.445134137142996_digits12
d4 =  3.754408661907416_digits12
p_low  = 0.02425_digits12
p_high = 1_digits12 - p_low
! Split into an inner and two outer regions which have separate fits
if(p < p_low) then
   q = sqrt(-2.0_digits12 * log(p))
   x = (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
      ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_digits12)
else if(p > p_high) then
   q = sqrt(-2.0_digits12 * log(1.0_digits12 - p))
   x = -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
      ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_digits12)
else
   q = p - 0.5_digits12
   r = q*q
   x = (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6)*q / &
      (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1.0_digits12)
endif

end subroutine norm_inv

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
!> wrappers for timers
!>
!> t_space is where we store the time information 
!> itemcount is the running count of how many times we've been called
!> maxitems is a limit on the number of times we want this timer to print.
!> do_sync overrides the default for whether we want to do a task sync
!> or not.  right now this code defaults to yes, sync before getting
!> the time.  for very large processor counts this increases overhead.

subroutine start_timer(t_space, itemcount, maxitems, do_sync)
 real(digits12), intent(out) :: t_space
 integer(i8),    intent(inout), optional :: itemcount
 integer(i8),    intent(in),    optional :: maxitems
 logical,        intent(in),    optional :: do_sync

logical :: sync_me

if (present(itemcount) .and. present(maxitems)) then
  itemcount = itemcount + 1
  if (itemcount > maxitems) then
     ! if called enough, this can roll over the integer limit.  
     ! set itemcount to maxitems+1 here to avoid this.
     ! also, go ahead and set the time because there is
     ! an option to print large time values even if over the
     ! number of calls limit in read_timer()

     itemcount = maxitems + 1
     call start_mpi_timer(t_space)
     return
  endif
endif

sync_me = .true.
if (present(do_sync)) sync_me = do_sync

if (sync_me) call task_sync()
call start_mpi_timer(t_space)

end subroutine start_timer

!--------------------------------------------------------------------
!>
!> t_space is where we store the time information 
!> label is the string to print out with the time.  limited to ~60 chars.
!> itemcount is the running count of how many times we've been called
!> maxitems is a limit on the number of times we want this timer to print.
!> do_sync overrides the default for whether we want to do a task sync
!> or not.  right now this code defaults to yes, sync before getting
!> the time.  for very large processor counts this increases overhead.
!> elapsed is an optional return of the value instead of only printing here

subroutine read_timer(t_space, label, itemcount, maxitems, do_sync, elapsed)
 real(digits12),   intent(in) :: t_space
 character(len=*), intent(in) :: label
 integer(i8),      intent(inout), optional :: itemcount
 integer(i8),      intent(in),    optional :: maxitems
 logical,          intent(in),    optional :: do_sync
 real(digits12),   intent(out),   optional :: elapsed

real(digits12) :: interval
logical :: sync_me
character(len=132) :: buffer

! if interval time (in seconds) is > this, go ahead and
! print even if item count is over limit.
integer(i8), parameter :: T_ALWAYS_PRINT = 1.0_r8

! get the time first and then figure out what we're doing

interval = read_mpi_timer(t_space)

sync_me = .true.
if (present(do_sync)) sync_me = do_sync

! if there's a limit on number of prints don't allow sync 
! because if we return early we could hang everyone else.
if (present(maxitems)) sync_me = .false.

! print out large values no matter what
! (large is defined above locally in this routine)
if (present(itemcount) .and. present(maxitems)) then
  if (interval < T_ALWAYS_PRINT .and. itemcount > maxitems) return
endif

! if syncing, wait and read the timer again
if (sync_me) then
   call task_sync()
   interval = read_mpi_timer(t_space)
endif

if (sync_me) then
   if (my_task_id() == 0) then
      write(buffer,'(A75,F15.8)') "timer: "//trim(label)//" time ", interval
      write(*,*) buffer
   endif
else
   write(buffer,'(A75,F15.8,A6,I7)') "timer: "//trim(label)//" time ", interval, " rank ", my_task_id()
   write(*,*) buffer
endif

if (present(elapsed)) elapsed = interval

end subroutine read_timer

!--------------------------------------------------------------------
!> log what the user has selected via the namelist choices

subroutine log_namelist_selections(num_special_cutoff, cache_override)

integer, intent(in) :: num_special_cutoff
logical, intent(in) :: cache_override

integer :: i

select case (filter_kind)
 case (1)
   msgstring = 'Ensemble Adjustment Kalman Filter (EAKF)'
 case (2)
   msgstring = 'Ensemble Kalman Filter (ENKF)'
 case (3)
   msgstring = 'Kernel filter'
 case (4)
   msgstring = 'observation space particle filter'
 case (5)
   msgstring = 'random draw from posterior'
 case (6)
   msgstring = 'deterministic draw from posterior with fixed kurtosis'
 case (7)
   msgstring = 'Boxcar'
 case (8)
   msgstring = 'Rank Histogram Filter'
 case (9) 
   msgstring = 'Local Particle Filter (Poterjoy)'
 case default
   call error_handler(E_ERR, 'assim_tools_init:', 'illegal filter_kind value, valid values are 1-9', &
                      source)
end select
call error_handler(E_MSG, 'assim_tools_init:', 'Selected filter type is '//trim(msgstring))

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
