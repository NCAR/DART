! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module assim_tools_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
! A variety of operations required by assimilation.

use      types_mod, only : r8
use  utilities_mod, only : file_exist, open_file, close_file, check_nml_error, get_unit, &
                           register_module, error_handler, E_ERR, E_MSG, logfileunit
use       sort_mod, only : index_sort 
use random_seq_mod, only : random_seq_type, random_gaussian, &
                           init_random_seq, random_uniform

use obs_sequence_mod, only : obs_sequence_type, obs_type, get_num_copies, get_num_qc, &
   init_obs, get_obs_from_key, get_obs_def, get_obs_values, get_qc, set_qc, &
   set_obs, get_copy_meta_data, read_obs_seq, destroy_obs_sequence, get_num_obs, &
   write_obs_seq
   
use obs_def_mod, only      : obs_def_type, get_obs_def_error_variance, get_obs_def_location
use cov_cutoff_mod, only   : comp_cov_factor
use obs_model_mod, only    : get_expected_obs, get_close_states
use reg_factor_mod, only   : comp_reg_factor
use location_mod, only     : location_type, get_dist
use time_manager_mod, only : time_type
use ensemble_manager_mod, only : get_ensemble_region, put_ensemble_region, ensemble_type

implicit none
private

public :: assim_tools_init, obs_increment, update_from_obs_inc, look_for_bias, &
   filter_assim, async_assim_region

type (random_seq_type) :: inc_ran_seq
logical :: first_inc_ran_call = .true.

! CVS Generated file description for error handling, do not edit
character(len=128), parameter :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

!============================================================================

!---- namelist with default values

logical :: prior_spread_correction = .false.
! Filter kind selects type of observation space filter
!      1 = EAKF filter
!      2 = ENKF
!      3 = Kernel filter
!      4 = particle filter
integer  :: filter_kind     = 1
real(r8) :: slope_threshold = 1.0_r8
logical :: do_parallel = .false.
integer :: num_domains = 1
character(len=129) :: parallel_command = './assim_filter.csh'

namelist / assim_tools_nml / prior_spread_correction, filter_kind, slope_threshold, &
   do_parallel, num_domains, parallel_command

!============================================================================

contains


subroutine assim_tools_init()
!============================================================================
! subroutine assim_tools_init()
!

integer :: iunit, ierr, io

call register_module(source, revision, revdate)

! Read namelist for run time control

if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   ierr = 1

   READBLOCK: do while(ierr /= 0)
      read(iunit, nml = assim_tools_nml, iostat = io)
      if ( io < 0 ) exit READBLOCK          ! end-of-file
      ierr = check_nml_error(io, 'assim_tools_nml')
   enddo READBLOCK
 
   call close_file(iunit)
endif

! Write the namelist values to the log file

call error_handler(E_MSG,'assim_tools_init','assim_tools namelist values',' ',' ',' ')
write(logfileunit, nml=assim_tools_nml)

end subroutine assim_tools_init

!-------------------------------------------------------------

subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc, &
                         slope, a, bias_ratio_out)

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
real(r8), intent(in)  :: slope
real(r8), intent(out) :: a
real(r8), intent(out), optional :: bias_ratio_out

if(filter_kind == 1) then
   call     obs_increment_eakf(ens, ens_size, obs, obs_var, obs_inc, slope, a, bias_ratio_out)
else if(filter_kind == 2) then
   call     obs_increment_enkf(ens, ens_size, obs, obs_var, obs_inc, slope, a, bias_ratio_out)
else if(filter_kind == 3) then
   call   obs_increment_kernel(ens, ens_size, obs, obs_var, obs_inc, slope, a, bias_ratio_out)
else if(filter_kind == 4) then
   call obs_increment_particle(ens, ens_size, obs, obs_var, obs_inc, slope, a, bias_ratio_out)
else 
   call error_handler(E_ERR,'obs_increment', &
              'Illegal value of filter_kind in assim_tools namelist [1-4 OK]', &
              source, revision, revdate)
endif

end subroutine obs_increment



subroutine obs_increment_eakf(ens, ens_size, obs, obs_var, obs_inc, &
                              slope, a, bias_ratio_out)
!========================================================================
!
! EAKF version of obs increment

integer, intent(in)   :: ens_size
real(r8), intent(in)  :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
real(r8), intent(in)  :: slope
real(r8), intent(out) :: a
real(r8), intent(out), optional :: bias_ratio_out

real(r8) :: prior_mean, new_mean, prior_var, var_ratio, sum_x
real(r8) :: error, diff_sd, ratio, factor, b, c, inf_obs_var
real(r8) :: inf_obs_var_inv, inf_prior_var, inf_prior_var_inv
real(r8) :: new_var, inf_ens(ens_size)

! Compute prior variance and mean from sample
sum_x      = sum(ens)
prior_mean = sum_x / ens_size
prior_var  = (sum(ens * ens) - sum_x**2 / ens_size) / (ens_size - 1)

if (obs_var /= 0.0_r8) then
   var_ratio = obs_var / (prior_var + obs_var)
   new_mean  = var_ratio * (prior_mean  + prior_var*obs / obs_var)
else
   if (prior_var /= 0.0_r8) then
      var_ratio = 0.0_r8
      new_mean  = obs
   else
      call error_handler(E_ERR,'obs_increment_eakf', &
           'Both obs_var and prior_var are zero. This is inconsistent', &
           source, revision, revdate)
!!$      if (prior_mean == obs) then
!!$         var_ratio = 1.0_r8
!!$         new_mean  = prior_mean
!!$      else
!!$         call error_handler(E_ERR,'obs_increment_eakf','Both obs_var and prior_var are zero and prior_mean is different from the obs. This is inconsistent',source, revision, revdate)
!!$      endif
   endif
endif

! THIS POINT CAN EVALUATE INCONSISTENCY if needed
if(abs(slope) > 0.0000001_r8 .or. present(bias_ratio_out)) then
   error = prior_mean - obs
   diff_sd = sqrt(obs_var + prior_var)
   ratio   = abs(error / diff_sd)
endif
if(present(bias_ratio_out)) bias_ratio_out = ratio

if(abs(slope) < 0.0000001_r8) then
   a = sqrt(var_ratio)
   obs_inc = a * (ens - prior_mean) + new_mean - ens
else

   ! Do confidence slope systematic error adjustment if requested
   ! Only modify if the ratio exceeds the threshold value
   if(ratio > slope_threshold) then
      b = (1.0_r8 / slope) ** (1.0_r8 / (slope - 1.0_r8))
      c = -1.0_r8 * b**slope
      factor = ratio / ((ratio - slope_threshold + b)**slope + c + slope_threshold)
   else
      factor = 1.0_r8
   endif

   ! Can now inflate by this ratio and then do adjustment
   inf_obs_var     = factor**2 * obs_var
   inf_obs_var_inv = 1.0_r8 / inf_obs_var

   ! Form inflated ensemble
   inf_ens       = prior_mean + factor * (ens - prior_mean)
   inf_prior_var = factor**2 * prior_var

   inf_prior_var_inv = 1.0_r8 / inf_prior_var
   new_var           = 1.0_r8 / (inf_prior_var_inv + inf_obs_var_inv)

   a       = sqrt(new_var * inf_prior_var_inv)
   obs_inc = a * (inf_ens - prior_mean) + new_mean - ens

endif

end subroutine obs_increment_eakf



subroutine obs_increment_particle(ens, ens_size, obs, obs_var, obs_inc, &
                    slope, a, bias_ratio_out)
!------------------------------------------------------------------------
!
! A observation space only particle filter implementation for a
! two step sequential update filter. Second version, 2 October, 2003.

integer, intent(in)             :: ens_size
real(r8), intent(in)            :: ens(ens_size), obs, obs_var
real(r8), intent(out)           :: obs_inc(ens_size)
real(r8), intent(in)            :: slope
real(r8), intent(out)           :: a
real(r8), intent(out), optional :: bias_ratio_out

real(r8) :: weight(ens_size), rel_weight(ens_size), cum_weight(0:ens_size)
real(r8) :: base, frac, new_val(ens_size), weight_sum
integer  :: i, j, indx(ens_size), ens_index(ens_size), new_index(ens_size)

! The factor a is not defined for particle filters
a = -1.0_r8

! Slope correction not currently implemented with particle filter
if( abs(slope) > 0.0000001_r8 ) call error_handler(E_ERR,'obs_increment_particle', &
                'Confidence slope bias correction is not implemented.', &
                 source, revision, revdate)

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

! Try sorting to make increments as small as possible
call index_sort(ens, ens_index, ens_size)
call index_sort(new_val, new_index, ens_size)
do i = 1, ens_size
   obs_inc(ens_index(i)) = new_val(new_index(i)) - ens(ens_index(i))
end do

end subroutine obs_increment_particle



subroutine obs_increment_enkf(ens, ens_size, obs, obs_var, obs_inc, &
                              slope, a, bias_ratio_out)
!========================================================================
! subroutine obs_increment_enkf(ens, ens_size, obs, obs_var, obs_inc)
!

! ENKF version of obs increment

integer, intent(in)             :: ens_size
real(r8), intent(in)            :: ens(ens_size), obs, obs_var
real(r8), intent(out)           :: obs_inc(ens_size)
real(r8), intent(in)            :: slope
real(r8), intent(out)           :: a
real(r8), intent(out), optional :: bias_ratio_out

real(r8) :: obs_var_inv
real(r8) :: prior_mean, prior_cov_inv, new_cov, new_mean(ens_size)
real(r8) :: sx, s_x2, prior_cov
real(r8) :: temp_mean, temp_obs(ens_size)
integer  :: ens_index(ens_size), new_index(ens_size)

integer  :: i

! The factor a is not defined for kernel filters
a = -1.0_r8

! Slope correction not currently implemented
if(abs(slope) > 0.0000001_r8) call error_handler(E_ERR,'obs_increment_enkf', &
          'Confidence slope bias correction is not implemented.', &
          source, revision, revdate)

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0_r8 / obs_var

! Compute prior mean and covariance
sx         = sum(ens)
s_x2       = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov  = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

prior_cov_inv = 1.0_r8 / prior_cov
new_cov       = 1.0_r8 / (prior_cov_inv + obs_var_inv)

! If this is first time through, need to initialize the random sequence
if(first_inc_ran_call) then
   call init_random_seq(inc_ran_seq)
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
   new_mean(i) = new_cov * (prior_cov_inv * ens(i) + temp_obs(i) / obs_var)
   obs_inc(i)  = new_mean(i) - ens(i)
end do

! Try sorting to make increments as small as possible
call index_sort(ens, ens_index, ens_size)
call index_sort(new_mean, new_index, ens_size)
do i = 1, ens_size
   obs_inc(ens_index(i)) = new_mean(new_index(i)) - ens(ens_index(i))
end do

end subroutine obs_increment_enkf



subroutine obs_increment_kernel(ens, ens_size, obs, obs_var, obs_inc, &
                                slope, a, bias_ratio_out)
!========================================================================
! subroutine obs_increment_kernel(ens, ens_size, obs, obs_var, obs_inc)
!

! Kernel version of obs increment

integer, intent(in)             :: ens_size
real(r8), intent(in)            :: ens(ens_size), obs, obs_var
real(r8), intent(out)           :: obs_inc(ens_size)
real(r8), intent(in)            :: slope
real(r8), intent(out)           :: a
real(r8), intent(out), optional :: bias_ratio_out

real(r8) :: obs_var_inv
real(r8) :: prior_mean, prior_cov_inv, new_cov, prior_cov
real(r8) :: sx, s_x2
real(r8) :: weight(ens_size), new_mean(ens_size)
real(r8) :: cum_weight, total_weight, cum_frac(ens_size)
real(r8) :: unif, norm, new_member(ens_size)

integer :: i, j, kernel, ens_index(ens_size), new_index(ens_size)

! The factor a is not defined for kernel filters
a = -1.0_r8

! Slope correction not currently implemented with kernel filter
if( abs(slope) > 0.0000001_r8 ) call error_handler(E_ERR,'obs_increment_kernel', &
              'Confidence slope bias correction is not implemented.', &
               source, revision, revdate)

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0_r8 / obs_var

! Compute prior mean and covariance
sx         = sum(ens)
s_x2       = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov  = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

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

! If this is first time through, need to initialize the random sequence
if(first_inc_ran_call) then
   call init_random_seq(inc_ran_seq)
   first_inc_ran_call = .false.
endif

! Generate a uniform random number and a Gaussian for each new member
do i = 1, ens_size
   unif = random_uniform(inc_ran_seq)
   ! Figure out which kernel it's in
   do j = 1, ens_size
      if(unif < cum_frac(j)) then
         kernel = j
         goto 10
      end if
   end do
10 continue

   ! Next calculate a unit normal in this kernel
   norm = random_gaussian(inc_ran_seq, 0.0_r8, sqrt(new_cov))
   ! Now generate the new ensemble member
   new_member(i) = new_mean(kernel) + norm
end do

! Try sorting to make increments as small as possible
call index_sort(ens, ens_index, ens_size)
call index_sort(new_member, new_index, ens_size)

do i = 1, ens_size
   obs_inc(ens_index(i)) = new_member(new_index(i)) - ens(ens_index(i))
end do

end subroutine obs_increment_kernel



subroutine update_from_obs_inc(obs, obs_inc, state, ens_size, &
               a, state_inc, reg_coef, correl_out)
!========================================================================

! Does linear regression of a state variable onto an observation and
! computes state variable increments from observation increments

integer, intent(in)             :: ens_size
real(r8), intent(in)            :: obs(ens_size), obs_inc(ens_size)
real(r8), intent(in)            :: state(ens_size)
real(r8), intent(inout)         :: a
real(r8), intent(out)           :: state_inc(ens_size), reg_coef
real(r8), intent(out), optional :: correl_out

real(r8) :: sum_x, t(ens_size), sum_t2, sum_ty
real(r8) :: inf_state(ens_size), correl

real(r8) :: sum_y, sum_y2
real(r8) :: factor, state_var_norm


! For efficiency, just compute regression coefficient here
sum_x  = sum(obs)
t      = obs - sum_x/ens_size
sum_t2 = sum(t * t)
sum_ty = sum(t * state)

if (sum_t2 /= 0.0_r8) then
   reg_coef = sum_ty/sum_t2
else
   reg_coef = 0.0_r8
endif

! Compute the sample correlation
if(present(correl_out) .or. prior_spread_correction) then
   sum_y          = sum(state)
   sum_y2         = sum(state*state)
   state_var_norm = sum_y2 - sum_y**2 / ens_size
   correl         = reg_coef * sqrt(sum_t2 / state_var_norm)
endif
if(present(correl_out)) correl_out = correl

! Then compute the increment as product of reg_coef and observation space increment
state_inc = reg_coef * obs_inc

if(.not. prior_spread_correction) return

! Won't work with non-EAKF filters
if( a < 0.0_r8 ) call error_handler(E_ERR,'update_from_obs_inc', & 
   'The prior_spread_correction algorithm only works with EAKF filters. Do not select prior_spread_correction = .true. with other filters.', &
              source, revision, revdate)

! Add in a slope factor for continuity
! Following flat triad works fairly well for base L96 cases; Use it for ens_size >= 20
if(ens_size >= 20) then
   if(abs(correl) < 0.3_r8) then
      factor = 1.0_r8 / (1.0_r8 - 1.0_r8 * (1.0_r8 - a) / 30.0_r8) 
   else if(abs(correl) < 0.6_r8) then
      factor = 1.0_r8 / (1.0_r8 - 1.0_r8 * (1.0_r8 - a) / 35.0_r8) 
   else
      factor = 1.0_r8
   endif
else
   !FOLLOWING ONE WORKS PRETTY WELL FOR 10 MEMBER ENSEMBLES IN L96!
   if(abs(correl) < 0.40_r8) then
      factor = 1.0_r8 / (1.0_r8 - 1.0_r8 * (1.0_r8 - a) / 9.0_r8)
   else if(abs(correl) < 0.7_r8) then
      factor = 1.0_r8 / (1.0_r8 - 1.0_r8 * (1.0_r8 - a) / 35.0_r8)
   else
      factor = 1.0_r8
   endif
endif

! Inflate the state estimate as required
inf_state = factor * (state - sum_y / ens_size) + (sum_y / ens_size)
state_inc = state_inc + inf_state - state

end subroutine update_from_obs_inc



!========================================================================


subroutine look_for_bias(ens, n, obs, obs_var, var_ratio)

integer, intent(in)   :: n
real(r8), intent(in)  :: ens(n), obs, obs_var
real(r8), intent(out) :: var_ratio

real(r8) :: sx, s_x2, prior_mean, prior_var, sq_err, tot_var

! Compute variance of the ensemble prior for this obs
sx         = sum(ens)
s_x2       = sum(ens * ens)
prior_mean = sx / n
prior_var  = (s_x2 - sx**2 / n) / (n - 1)

! Variance of difference between obs and mean should be sum of variances
sq_err    = (obs - prior_mean)**2
!!!sq_err = sum((obs - ens)**2) / n
tot_var   = obs_var + prior_var
var_ratio = sq_err / tot_var

end subroutine look_for_bias

!========================================================================

! DONT USE THIS RIGHT NOW
subroutine seq_filter_assim(ens, ens_obs, compute_obs, ens_size, &
     num_obs_in_set, num_groups, seq, keys, obs_val_index, &
     confidence_slope, cutoff, save_reg_series, reg_series_unit)

integer, intent(in) :: ens_size, num_groups, num_obs_in_set, keys(num_obs_in_set)
real(r8), intent(inout) :: ens(:, :), ens_obs(ens_size, num_obs_in_set)
logical, intent(in) :: compute_obs(num_obs_in_set)
integer, intent(in) :: obs_val_index, reg_series_unit
type(obs_sequence_type), intent(inout) :: seq
real(r8), intent(in) :: confidence_slope, cutoff
logical, intent(in) :: save_reg_series

integer :: i, j, jjj, k, kkk, istatus, ind, grp_size, group, grp_bot, grp_top
integer :: num_close_ptr(1), indx
type(obs_type) :: observation
type(obs_def_type) :: obs_def
real(r8) :: obs(num_obs_in_set), obs_err_var(num_obs_in_set), cov_factor, regress(num_groups)
real(r8) :: qc(1), obs_inc(ens_size), a_returned(num_groups), obs_dist
real(r8) :: ens_inc(ens_size)
integer, parameter :: first_num_close = 100000
integer :: order(num_obs_in_set)
integer, allocatable :: close_ptr(:, :)
real(r8), allocatable :: dist_ptr(:, :)
real(r8) :: swath(ens_size), reg_factor
type(location_type) :: obs_loc(num_obs_in_set)

! Set an initial size for the close state pointers
allocate(close_ptr(1, first_num_close), dist_ptr(1, first_num_close))

!write(*, *) 'initializing an observation temporary', get_num_copies(seq), get_num_qc(seq)
! Construnct an observation temporary
call init_obs(observation, get_num_copies(seq), get_num_qc(seq))

! Get the locations for all of the observations
do i = 1, num_obs_in_set
   call get_obs_from_key(seq, keys(i), observation)
   call get_obs_def(observation, obs_def)
   obs_loc(i) = get_obs_def_location(obs_def)
end do

! Reorder the observations to get rid of the ones that can't be recomputed first
indx = 1
do i = 1, num_obs_in_set
   if(.not. compute_obs(i)) then
      order(indx) = i
      indx = indx + 1
   endif
end do

! Then put in all the ones that can be recomputed
do i = 1, num_obs_in_set
   if(compute_obs(i)) then
      order(indx) = i
      indx = indx + 1
   endif
end do

do i = 1, num_obs_in_set
   write(*, *) 'order ', i, order(i)
end do

!---------------------- Sequential filter section --------------------------
! Loop through each observation in the set
!Observations : do j = 1, num_obs_in_set
Observations : do jjj = 1, num_obs_in_set
   j = jjj
   
   ! Get this observation in a temporary and get its qc value for possible update
   call get_obs_from_key(seq, keys(j), observation)

   ! Get the observation value and the error variance
   call get_obs_def(observation, obs_def)
   ! Get the value of the observation and the error variance
   call get_obs_values(observation, obs(j:j), obs_val_index)
   obs_err_var(j) = get_obs_def_error_variance(obs_def)
   ! Get the qc value set so far
   call get_qc(observation, qc, 1)
   ! Just skip observations that have failed prior qc ( >4 for now)
   if(qc(1) > 4.01_r8) cycle Observations
   
   ! Compute the ensemble prior for this ob
   if(compute_obs(j)) then
      do k = 1, ens_size
         ! Only compute forward operator if allowed
         call get_expected_obs(seq, keys(j:j), ens(k, :), ens_obs(k:k, j), istatus)
         ! Inability to compute forward operator implies skip this observation
         if (istatus > 0) then
            qc(1) = qc(1) + 4000.0_r8
            goto 333
         endif
      end do
   endif
   ! Divide ensemble into num_groups groups
   grp_size = ens_size / num_groups
   Group1: do group = 1, num_groups
      grp_bot = (group - 1) * grp_size + 1
      grp_top = grp_bot + grp_size - 1

      ! Call obs_increment to do observation space
      call obs_increment(ens_obs(grp_bot:grp_top, j), ens_size/num_groups, obs(j), &
         obs_err_var(j), obs_inc(grp_bot:grp_top), confidence_slope, &
         a_returned(group))
   end do Group1

   ! Getting close states for each scalar observation for now
   ! Need model state for some distance computations in sigma, have
   ! to pick some state, only current ones are ensemble, just pass 1
222   call get_close_states(seq, keys(j), 2.0_r8*cutoff, num_close_ptr(1), &
      close_ptr(1, :), dist_ptr(1, :), ens(1, :))
   if(num_close_ptr(1) < 0) then
      deallocate(close_ptr, dist_ptr)
      allocate(close_ptr(1, -1 * num_close_ptr(1)), dist_ptr(1, -1 * num_close_ptr(1)))
      goto 222
   endif

   ! Now loop through each close state variable for this observation
   do k = 1, num_close_ptr(1)
      ind = close_ptr(1, k)
      ! Compute distance dependent envelope
      cov_factor = comp_cov_factor(dist_ptr(1, k), cutoff)
      ! WARNING: TEST FOR REPRODUCING with next line
      if(cov_factor == 0) cycle

      ! Get the ensemble elements for this state variable and do regression
       swath = ens(:, ind)

      ! Loop through the groups
      Group2: do group = 1, num_groups
         grp_bot = (group - 1) * grp_size + 1
         grp_top = grp_bot + grp_size - 1
         call update_from_obs_inc(ens_obs(grp_bot:grp_top, j), &
            obs_inc(grp_bot:grp_top), swath(grp_bot:grp_top), ens_size/num_groups, &
            a_returned(group), ens_inc(grp_bot:grp_top), regress(group))
      end do Group2

      ! Compute an information factor for impact of this observation on this state
      if(num_groups == 0) then
          reg_factor = 1.0_r8
      else
! PROBLEM WITH TIME INFO ON REGRESSION COEF
!!!         reg_factor = comp_reg_factor(num_groups, regress, i, j, ind)
         reg_factor = comp_reg_factor(num_groups, regress, 1, j, ind)
      endif
      if(save_reg_series) write(reg_series_unit, *) j, k, reg_factor

      ! COV_FACTOR NOW COMES IN HERE FOR USE WITH GROUPS; JLA 11/19/03
      reg_factor = min(reg_factor, cov_factor)

      ! Do the final update for this state variable
      ens(:, ind) = ens(:, ind) + reg_factor * ens_inc(:)
   end do


   ! Also need to update all obs variable priors that cannot be recomputed
   ! Obs already used will never be used again, be careful not to do the jth one here!
   do kkk = jjj + 1, num_obs_in_set
      k = order(kkk)

      ! Don't need to do this if observation can be recomputed
      if(compute_obs(k)) cycle
      ! Get location of the two observations (can be done once for efficiency)
      obs_dist = get_dist(obs_loc(k), obs_loc(j))
      ! Compute distance dependent envelope
      cov_factor = comp_cov_factor(obs_dist, cutoff)
      if(cov_factor == 0.0) cycle
      ! Get the ensemble elements for this state variable and do regression
       swath = ens_obs(:, k)

      ! Loop through the groups
      Group3: do group = 1, num_groups
         grp_bot = (group - 1) * grp_size + 1
         grp_top = grp_bot + grp_size - 1
         call update_from_obs_inc(ens_obs(grp_bot:grp_top, j), &
            obs_inc(grp_bot:grp_top), swath(grp_bot:grp_top), ens_size/num_groups, &
            a_returned(group), ens_inc(grp_bot:grp_top), regress(group))
      end do Group3

      ! Compute an information factor for impact of this observation on this state
      if(num_groups == 0) then
          reg_factor = 1.0_r8
      else
! PROBLEM WITH TIME INFO ON REGRESSION COEF
!!!         reg_factor = comp_reg_factor(num_groups, regress, i, j, ind)
         reg_factor = comp_reg_factor(num_groups, regress, 1, j, ind)
      endif
      if(save_reg_series) write(reg_series_unit, *) j, k, reg_factor

      ! COV_FACTOR NOW COMES IN HERE FOR USE WITH GROUPS; JLA 11/19/03
      reg_factor = min(reg_factor, cov_factor)

      ! Do the final update for this state variable
      ens_obs(:, k) = ens_obs(:, k) + reg_factor * ens_inc(:)

   end do

!------------------

   ! Write out the updated quality control for this observation
   333 call set_qc(observation, qc, 1)
   call set_obs(seq, observation, keys(j))

end do Observations

deallocate(close_ptr, dist_ptr)

end subroutine seq_filter_assim

!===========================================================================

subroutine filter_assim_region(my_domain, domain_size, ens_size, model_size, &
   num_groups, num_obs_in_set, obs_val_index, confidence_slope, cutoff, &
   save_reg_series, reg_series_unit, ens, ens_obs_in, compute_obs, seq, keys, my_state)

integer, intent(in) :: my_domain, domain_size
integer, intent(in) :: ens_size, model_size, num_groups, num_obs_in_set, keys(num_obs_in_set)
real(r8), intent(inout) :: ens(ens_size, domain_size)
real(r8), intent(in) :: ens_obs_in(ens_size, num_obs_in_set)
logical, intent(in) :: compute_obs(num_obs_in_set)
integer, intent(in) :: obs_val_index, reg_series_unit
type(obs_sequence_type), intent(inout) :: seq
real(r8), intent(in) :: confidence_slope, cutoff
logical, intent(in) :: save_reg_series, my_state(model_size)

integer :: i, j, jjj, k, kkk, istatus, ind, grp_size, group, grp_bot, grp_top
integer :: num_close_ptr(1), indx
type(obs_type) :: observation
type(obs_def_type) :: obs_def
real(r8) :: obs(num_obs_in_set), obs_err_var(num_obs_in_set), cov_factor, regress(num_groups)
real(r8) :: qc(1), obs_inc(ens_size), a_returned(num_groups), obs_dist
real(r8) :: ens_inc(ens_size), ens_obs(ens_size, num_obs_in_set)
integer, parameter :: first_num_close = 100000
integer :: order(num_obs_in_set)
integer, allocatable :: close_ptr(:, :)
real(r8), allocatable :: dist_ptr(:, :)
real(r8) :: swath(ens_size), reg_factor
type(location_type) :: obs_loc(num_obs_in_set)
logical :: close_to_any

integer :: inv_indices(model_size)
integer, allocatable :: indices(:)

! Generate an array with the indices of my state variables
allocate(indices(domain_size))
inv_indices = 0
indx = 1
do i = 1, model_size
   if(my_state(i)) then
      indices(indx) = i
      inv_indices(i) = indx
      indx = indx + 1
   endif
end do 

do i = 1, model_size
!   write(*, *) 'my_state, inv_index ', i, my_state(i), inv_indices(i)
end do
do i = 1, domain_size
!   write(*, *) 'indices ', i, indices(i)
end do

! Need copy of obs to be modified
ens_obs = ens_obs_in

! Set an initial size for the close state pointers
allocate(close_ptr(1, first_num_close), dist_ptr(1, first_num_close))

!write(*, *) 'initializing an observation temporary', get_num_copies(seq), get_num_qc(seq)
! Construnct an observation temporary
call init_obs(observation, get_num_copies(seq), get_num_qc(seq))

! Get the locations for all of the observations
do i = 1, num_obs_in_set
   call get_obs_from_key(seq, keys(i), observation)
   call get_obs_def(observation, obs_def)
   obs_loc(i) = get_obs_def_location(obs_def)
end do

! Reorder the observations to get rid of the ones that can't be recomputed first
indx = 1
do i = 1, num_obs_in_set
   if(.not. compute_obs(i)) then
      order(indx) = i
      indx = indx + 1
   endif
end do

! Then put in all the ones that can be recomputed
do i = 1, num_obs_in_set
   if(compute_obs(i)) then
      order(indx) = i
      indx = indx + 1
   endif
end do

do i = 1, num_obs_in_set
!!!   write(*, *) 'order ', i, order(i)
end do

!---------------------- Sequential filter section --------------------------
! Loop through each observation in the set
!Observations : do j = 1, num_obs_in_set
Observations : do jjj = 1, num_obs_in_set
   j = jjj
   
   ! Get this observation in a temporary and get its qc value for possible update
   call get_obs_from_key(seq, keys(j), observation)

   ! Get the observation value and the error variance
   call get_obs_def(observation, obs_def)
   ! Get the value of the observation and the error variance
   call get_obs_values(observation, obs(j:j), obs_val_index)
   obs_err_var(j) = get_obs_def_error_variance(obs_def)
   ! Get the qc value set so far
   call get_qc(observation, qc, 1)
   ! Just skip observations that have failed prior qc ( >4 for now)
   if(qc(1) > 4.01_r8) cycle Observations
   
   ! Compute the ensemble prior for this ob
   if(compute_obs(j)) then
      do k = 1, ens_size
         ! Only compute forward operator if allowed
         if(1 == 1) then
            write(*, *) 'OOOPS Computing get_expected_obs'
            stop
         endif
         call get_expected_obs(seq, keys(j:j), ens(k, :), ens_obs(k:k, j), istatus)
         ! Inability to compute forward operator implies skip this observation
         if (istatus > 0) then
            qc(1) = qc(1) + 4000.0_r8
            goto 333
         endif
      end do
   endif
   ! Divide ensemble into num_groups groups
   grp_size = ens_size / num_groups
   Group1: do group = 1, num_groups
      grp_bot = (group - 1) * grp_size + 1
      grp_top = grp_bot + grp_size - 1

      ! Call obs_increment to do observation space
      call obs_increment(ens_obs(grp_bot:grp_top, j), ens_size/num_groups, obs(j), &
         obs_err_var(j), obs_inc(grp_bot:grp_top), confidence_slope, &
         a_returned(group))
   end do Group1

   ! Getting close states for each scalar observation for now
   ! Need model state for some distance computations in sigma, have
   ! to pick some state, only current ones are ensemble, just pass 1
222   call get_close_states(seq, keys(j), 2.0_r8*cutoff, num_close_ptr(1), &
      close_ptr(1, :), dist_ptr(1, :), ens(1, :))
   if(num_close_ptr(1) < 0) then
      deallocate(close_ptr, dist_ptr)
      allocate(close_ptr(1, -1 * num_close_ptr(1)), dist_ptr(1, -1 * num_close_ptr(1)))
      goto 222
   endif

   ! Determine if this observation is close to ANY state variable for this process
   close_to_any = .false.
   ! Now loop through each close state variable for this observation
   CLOSE_STATE: do k = 1, num_close_ptr(1)
      ind = close_ptr(1, k)

      ! If this state variable is not in the list to be updated for this PE, skip
      if(.not. my_state(ind)) cycle CLOSE_STATE
      close_to_any = .true.


      ! Compute distance dependent envelope
      cov_factor = comp_cov_factor(dist_ptr(1, k), cutoff)
      ! WARNING: TEST FOR REPRODUCING with next line
      if(cov_factor <= 0) cycle CLOSE_STATE

      ! Get the ensemble elements for this state variable and do regression
       swath = ens(:, inv_indices(ind))

      ! Loop through the groups
      Group2: do group = 1, num_groups
         grp_bot = (group - 1) * grp_size + 1
         grp_top = grp_bot + grp_size - 1
         call update_from_obs_inc(ens_obs(grp_bot:grp_top, j), &
            obs_inc(grp_bot:grp_top), swath(grp_bot:grp_top), ens_size/num_groups, &
            a_returned(group), ens_inc(grp_bot:grp_top), regress(group))
      end do Group2

      ! Compute an information factor for impact of this observation on this state
      if(num_groups == 0) then
          reg_factor = 1.0_r8
      else
! PROBLEM WITH TIME INFO ON REGRESSION COEF
!!!         reg_factor = comp_reg_factor(num_groups, regress, i, j, ind)
         reg_factor = comp_reg_factor(num_groups, regress, 1, j, ind)
      endif
      if(save_reg_series) write(reg_series_unit, *) j, k, reg_factor

      ! COV_FACTOR NOW COMES IN HERE FOR USE WITH GROUPS; JLA 11/19/03
      reg_factor = min(reg_factor, cov_factor)

      ! Do the final update for this state variable
      ens(:, inv_indices(ind)) = ens(:, inv_indices(ind)) + reg_factor * ens_inc(:)
   end do CLOSE_STATE


   ! If this observation is not close to ANY state variable, can skip obs update
   ! This will NOT bitwise reproduce. Notion is that we could reorder the use
   ! of observations so that this one would come after all meaningful ones.
   if(.not. close_to_any) goto 333

   ! Also need to update all obs variable priors that cannot be recomputed
   ! Obs already used will never be used again, be careful not to do the jth one here!
   do kkk = jjj + 1, num_obs_in_set
      k = order(kkk)

      ! Don't need to do this if observation can be recomputed
      if(compute_obs(k)) cycle
      ! Get location of the two observations (can be done once for efficiency)
      obs_dist = get_dist(obs_loc(k), obs_loc(j))
      ! Compute distance dependent envelope
      cov_factor = comp_cov_factor(obs_dist, cutoff)
      if(cov_factor <= 0.0) cycle
      ! Get the ensemble elements for this state variable and do regression
       swath = ens_obs(:, k)

      ! Loop through the groups
      Group3: do group = 1, num_groups
         grp_bot = (group - 1) * grp_size + 1
         grp_top = grp_bot + grp_size - 1
         call update_from_obs_inc(ens_obs(grp_bot:grp_top, j), &
            obs_inc(grp_bot:grp_top), swath(grp_bot:grp_top), ens_size/num_groups, &
            a_returned(group), ens_inc(grp_bot:grp_top), regress(group))
      end do Group3

      ! Compute an information factor for impact of this observation on this state
      if(num_groups == 0) then
          reg_factor = 1.0_r8
      else
! PROBLEM WITH TIME INFO ON REGRESSION COEF
!!!         reg_factor = comp_reg_factor(num_groups, regress, i, j, ind)
         reg_factor = comp_reg_factor(num_groups, regress, 1, j, ind)
      endif
      !!!if(save_reg_series) write(reg_series_unit, *) j, k, reg_factor

      ! COV_FACTOR NOW COMES IN HERE FOR USE WITH GROUPS; JLA 11/19/03
      reg_factor = min(reg_factor, cov_factor)

      ! Do the final update for this state variable
      ens_obs(:, k) = ens_obs(:, k) + reg_factor * ens_inc(:)

   end do

!------------------

   ! Write out the updated quality control for this observation
   333 call set_qc(observation, qc, 1)
   call set_obs(seq, observation, keys(j))

end do Observations

deallocate(close_ptr, dist_ptr)

end subroutine filter_assim_region


!---------------------------------------------------------------------------

subroutine filter_assim(ens_handle, ens_obs, compute_obs_in, ens_size, model_size, num_obs_in_set, &
   num_groups, seq, keys, confidence_slope, cutoff, save_reg_series, reg_series_unit, &
   obs_sequence_in_name)

integer, intent(in) :: ens_size, model_size, num_groups, num_obs_in_set, keys(num_obs_in_set)
type(ensemble_type), intent(in) :: ens_handle
real(r8), intent(in) :: ens_obs(ens_size, num_obs_in_set)
logical, intent(in) :: compute_obs_in(num_obs_in_set)
integer, intent(in) :: reg_series_unit
type(obs_sequence_type), intent(inout) :: seq
real(r8), intent(in) :: confidence_slope, cutoff
logical, intent(in) :: save_reg_series
character(len=*), intent(in) :: obs_sequence_in_name

integer :: i, j, domain_size, indx, obs_val_index, iunit, control_unit
integer :: base_size, remainder, domain_top, domain_bottom
logical :: compute_obs(num_obs_in_set), my_state(model_size)
real(r8), allocatable :: ens(:, :)
integer, allocatable :: indices(:)
type(time_type) :: ens_time(ens_size)
character(len=28) :: in_file_name(num_domains), out_file_name(num_domains)
character(len=129) :: temp_obs_seq_file, errstring

compute_obs = compute_obs_in

! Set compute obs to false for all obs for easy parallel testing
compute_obs = .false.

! Determine which copy has actual obs
do j = 1, get_num_copies(seq)
   obs_val_index = j
   ! Need to look for 'observation'
   if(index(get_copy_meta_data(seq, j), 'observation') > 0) goto 444
end do
! Falling of end means 'observations' not found; die
call error_handler(E_ERR,'filter_assim', &
        'Did not find observation copy with metadata "observation"', &
         source, revision, revdate)



! DOING A SINGLE DOMAIN AND ALLOWING RECOMPUTATION OF ALL OBS GIVES TRADITIONAL
! SEQUENTIAL ANSWER
444 continue
      
! Write out the most up-to-date obs_sequence file, too, which includes qc
if(do_parallel) then
   temp_obs_seq_file = 'filter_assim_obs_seq'
   call write_obs_seq(seq, temp_obs_seq_file)
endif

! Having more domains than state variables is silly
if(num_domains > model_size) then
call error_handler(E_ERR,'filter_assim', &
        'Having more domains than state variables does not make sense', &
         source, revision, revdate)
endif

base_size = model_size / num_domains
remainder = model_size - base_size * num_domains
domain_top = 0
do j = 1, num_domains
   domain_bottom = domain_top + 1

   ! Find out which state variables are in my domain
   !call whats_in_my_domain(num_domains, j, my_state, domain_size)
   ! Temporarily just do the domain stuff here
   domain_size = base_size
   if(j <= remainder) domain_size = domain_size + 1
   domain_top = domain_bottom + domain_size - 1
   my_state = .false.
   my_state(domain_bottom : domain_top) = .true.

   ! Generate an array with the indices of my state variables
   allocate(ens(ens_size, domain_size), indices(domain_size))

   indx = 1
   do i = 1, model_size
      if(my_state(i)) then
         indices(indx) = i
         indx = indx + 1
      endif
   end do 

   ! Get the ensemble (whole thing for now) from storage
   call get_ensemble_region(ens_handle, ens, ens_time, state_vars_in = indices)

! Do ensemble filter update for this region
   if(.not. do_parallel) then
      call filter_assim_region(j, domain_size, ens_size, model_size, num_groups, &
         num_obs_in_set, obs_val_index, confidence_slope, cutoff, save_reg_series, &
         reg_series_unit, ens, ens_obs, compute_obs, seq, keys, my_state)
   else
   ! Loop to write out arguments for each region
      if(j < 10) then
         write(in_file_name(j), 11) 'filter_assim_region__in', j
         write(out_file_name(j), 11) 'filter_assim_region_out', j
      else if(j < 100) then
         write(in_file_name(j), 21) 'filter_assim_region__in', j
         write(out_file_name(j), 21) 'filter_assim_region_out', j
      else if(j < 1000) then
         write(in_file_name(j), 31) 'filter_assim_region__in', j
         write(out_file_name(j), 31) 'filter_assim_region_out', j
      else if(j < 10000) then
         write(in_file_name(j), 41) 'filter_assim_region__in', j
         write(out_file_name(j), 41) 'filter_assim_region_out', j
      else
         write(errstring,*)'Trying to use ',ens_size,' model states -- too many.'
         call error_handler(E_MSG,'filter_assim',errstring,source,revision,revdate)
         call error_handler(E_ERR,'filter_assim','Use less than 10000 model states.',source,revision,revdate)
      endif

 11   format(a23, i1)
 21   format(a23, i2)
 31   format(a23, i3)
 41   format(a23, i4)

! Development test of passing this interface by file in preparation for separate executables
      iunit = get_unit()
      open(unit = iunit, file = in_file_name(j), action = 'write', form = 'formatted', &
         status = 'replace')
      write(iunit, *) num_domains, j, domain_size, ens_size, model_size, num_groups, &
         num_obs_in_set, obs_val_index, confidence_slope, cutoff, save_reg_series, &
         reg_series_unit, obs_sequence_in_name
      write(iunit, *) ens, ens_obs, compute_obs, keys, my_state
      close(iunit)

! Following line allows single processor test of communications
!!!!      call async_assim_region(in_file_name(j), out_file_name(j))

   endif

   ! Put this region into storage for single executable which is now finished
   if(.not. do_parallel) call put_ensemble_region(ens_handle, ens, ens_time, state_vars_in = indices)

   deallocate(ens, indices)

end do

! All the files are out there, now need to spawn off jobs if parallel
if(do_parallel) then

   ! Write out the file names to a control file
   control_unit = get_unit()
   open(unit = control_unit, file = 'assim_region_control')
   write(control_unit, *) num_domains
   do i = 1, num_domains
      write(control_unit, '(a28)' ) in_file_name(i)
      write(control_unit, '(a28)' ) out_file_name(i)
   end do
   close(control_unit)

! Spawn the job to do parallel assims
!   call system('echo go > batchflag; ./assim_filter.csh; sleep 1')
   call system('echo go > batchflag; '//parallel_command//' ; sleep 1')

   do
      if(file_exist('batchflag')) then
         call system('sleep 10')
      else
         exit
      endif
   end do

   base_size = model_size / num_domains
   remainder = model_size - base_size * num_domains
   domain_top = 0
   do j = 1, num_domains
   ! Find out which state variables are in my domain
   !call whats_in_my_domain(num_domains, j, my_state, domain_size)
   ! Temporarily just do the domain stuff here
      domain_bottom = domain_top + 1

      domain_size = base_size
      if(j <= remainder) domain_size = domain_size + 1
      domain_top = domain_bottom + domain_size - 1
      my_state = .false.
      my_state(domain_bottom : domain_top) = .true.

      ! Generate an array with the indices of my state variables
      allocate(ens(ens_size, domain_size), indices(domain_size))

      indx = 1
      do i = 1, model_size
         if(my_state(i)) then
            indices(indx) = i
            indx = indx + 1
         endif
      end do

      ! Read in the update ensemble region (also need the qc eventually)
      iunit = get_unit()
      open(unit = iunit, file = out_file_name(j), action = 'read', form = 'formatted')
      read(iunit, *) ens
      close(iunit)
      call put_ensemble_region(ens_handle, ens, ens_time, state_vars_in = indices)
write(*, *) 'preparing to deallocate ens and indices', j
      deallocate(ens, indices)
write(*, *) 'finished with deallocate ens and indices', j
   end do

endif
write(*, *) 'done with subroutine filter_assim'

end subroutine filter_assim


!-----------------------------------------------------------------------

subroutine async_assim_region(in_file_name, out_file_name)

character(len = *), intent(in) :: in_file_name, out_file_name

! DONT FORGET THE QC!!!
! AND IN GENERAL NEED TO WORK ON REGSERIES STUFF!!!

integer :: n_domains, my_domain, domain_size, ens_size, model_size, num_groups
integer :: num_obs_in_set, obs_val_index, reg_series_unit
real(r8) :: confidence_slope, cutoff
character(len = 129) :: obs_sequence_in_name
logical :: save_reg_series
real(r8), allocatable :: ens(:, :), ens_obs(:, :)
integer, allocatable :: keys(:)
logical, allocatable :: my_state(:), compute_obs(:)

type(obs_sequence_type) :: seq2
integer :: iunit
type(obs_type) :: obs1, obs2
real(r8) :: val1(10), val2(10)

! Read in all the arguments from the file
iunit = get_unit()
open(unit = iunit, file = in_file_name, action = 'read', form = 'formatted')
read(iunit, *) n_domains, my_domain, domain_size, ens_size, model_size, num_groups, &
   num_obs_in_set, obs_val_index, confidence_slope, cutoff, save_reg_series, &
   reg_series_unit, obs_sequence_in_name

write(*, *) 'first try, n_domains '

! Allocate storage
allocate(ens(ens_size, domain_size), ens_obs(ens_size, num_obs_in_set), &
   compute_obs(num_obs_in_set), keys(num_obs_in_set), my_state(model_size))
read(iunit, *) ens, ens_obs, compute_obs, keys, my_state

close(iunit)

! Now need to open the obs sequence file to proceed with call to region assim
!val1 = 0.0; val2 = 0.0
obs_sequence_in_name = 'filter_assim_obs_seq'
call read_obs_seq(obs_sequence_in_name, 0, 0, 0, seq2)

! Do ensemble filter update for this region
write(*, *) 'calling filter_assim_region in async_filter_assim, REALLY'
write(*, *) 'n_domains is ', n_domains
write(*, *) n_domains, my_domain, domain_size, ens_size, model_size, num_groups, num_obs_in_set, obs_val_index, confidence_slope, cutoff, save_reg_series, reg_series_unit
call filter_assim_region(my_domain, domain_size, ens_size, model_size, num_groups, &
   num_obs_in_set, obs_val_index, confidence_slope, cutoff, save_reg_series, &
   reg_series_unit, ens, ens_obs, compute_obs, seq2, keys, my_state)
write(*, *) 'back from filter_assim_region in async_filter_assim'

! Write out the updated ensemble region, need to do QC also
iunit = get_unit()
open(unit = iunit, file = out_file_name, action = 'write', form = 'formatted', &
   status = 'replace')
write(iunit, *) ens
close(iunit)

call destroy_obs_sequence(seq2)

write(*, *) 'done with destroy_obs_seq and async_assim_region'

end subroutine async_assim_region

!========================================================================
! end module assim_tools_mod
!========================================================================

end module assim_tools_mod
