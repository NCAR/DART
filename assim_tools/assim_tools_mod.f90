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

use      types_mod, only : r8, missing_r8
use  utilities_mod, only : file_exist, open_file, close_file, check_nml_error, get_unit, &
                           register_module, error_handler, E_ERR, E_MSG, logfileunit
use       sort_mod, only : index_sort 
use random_seq_mod, only : random_seq_type, random_gaussian, &
                           init_random_seq, random_uniform

use obs_sequence_mod, only : obs_sequence_type, obs_type, get_num_copies, get_num_qc, &
   init_obs, get_obs_from_key, get_obs_def, get_obs_values, get_qc, set_qc, &
   set_obs, get_copy_meta_data, read_obs_seq, destroy_obs_sequence, get_num_obs, &
   write_obs_seq, destroy_obs
   
use obs_def_mod, only      : obs_def_type, get_obs_def_error_variance, get_obs_def_location
use cov_cutoff_mod, only   : comp_cov_factor
use obs_model_mod, only    : get_expected_obs, get_close_states
use reg_factor_mod, only   : comp_reg_factor
use location_mod, only     : location_type, get_dist
use time_manager_mod, only : time_type
use ensemble_manager_mod, only : get_ensemble_region, put_ensemble_region, ensemble_type, &
                                 transpose_ens_to_regions, transpose_regions_to_ens, &
                                 get_ensemble_member, put_region_by_number, get_region_by_number

implicit none
private

public :: assim_tools_init, obs_increment, update_from_obs_inc, look_for_bias, &
   filter_assim, async_assim_region

type (random_seq_type) :: inc_ran_seq
logical :: first_inc_ran_call = .true.
character(len = 129) :: errstring

! CVS Generated file description for error handling, do not edit
character(len=128), parameter :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

!============================================================================

!---- namelist with default values

! Filter kind selects type of observation space filter
!      1 = EAKF filter
!      2 = ENKF
!      3 = Kernel filter
!      4 = particle filter
integer  :: filter_kind     = 1
logical  :: sort_obs_inc = .false.
real(r8) :: cov_inflate = -1.0
real(r8) :: cov_inflate_sd = 0.05
real(r8) :: sd_lower_bound = 0.05
logical  :: deterministic_cov_inflate = .true.
logical  :: start_from_assim_restart = .false.
character(len = 129) :: assim_restart_in_file_name = 'assim_ics'
character(len = 129) :: assim_restart_out_file_name = 'assim_restart'
integer :: do_parallel = 0
integer :: num_domains = 1
character(len=129) :: parallel_command = './assim_filter.csh'

namelist / assim_tools_nml / filter_kind, sort_obs_inc, cov_inflate, &
   cov_inflate_sd, sd_lower_bound, deterministic_cov_inflate, &
   start_from_assim_restart, assim_restart_in_file_name, &
   assim_restart_out_file_name, do_parallel, num_domains, parallel_command

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
write(     *     , nml=assim_tools_nml)
! If requested, read cov_inflate and cov_inflate_sd from a restart file
!if(start_from_assim_restart) then
!   restart_unit = get_unit()
!   open(unit = restart_unit, file = assim_restart_in_file_name)
!   read(unit, *) cov_inflate, cov_inflate_sd
!   close(unit)
!endif

write(*, *) 'sd controls are ', cov_inflate, cov_inflate_sd, sd_lower_bound

!!! PROBLEMS WITH REGIONS AND STATE!!!

! Look for an error in the do_parallel options
if(do_parallel /= 0 .and. do_parallel /= 2 .and. do_parallel /= 3) then
   write(errstring, *) 'do_parallel option ', do_parallel, ' not supported'
   call error_handler(E_ERR,'assim_tools_init', errstring, source, revision, revdate)
endif

end subroutine assim_tools_init

!-------------------------------------------------------------

subroutine obs_increment(ens_in, ens_size, obs, obs_var, obs_inc)
                         

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens_in(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens(ens_size), inflate_inc(ens_size)
real(r8) :: sum_x, prior_mean, prior_var, new_cov_inflate, new_cov_inflate_sd
real(r8) :: rand_sd, new_val(ens_size), enhanced_inflate
integer :: i, ens_index(ens_size), new_index(ens_size)

! If observation space inflation is being done, compute the initial 
! increments and update the inflation factor and its standard deviation
! as needed. cov_inflate < 0 means don't do any of this.
if(cov_inflate > 0.0) then
   ! Compute prior variance and mean from sample
   sum_x      = sum(ens_in)
   prior_mean = sum_x / ens_size
   prior_var  = (sum(ens_in * ens_in) - sum_x**2 / ens_size) / (ens_size - 1)

   ! If cov_inflate_sd is <= 0, just retain current cov_inflate setting
   if(cov_inflate_sd > 0.0) then
      call bayes_cov_inflate(prior_mean, sqrt(prior_var), obs, sqrt(obs_var), &
         cov_inflate, cov_inflate_sd, new_cov_inflate, new_cov_inflate_sd)

      ! Keep the covariance inflation > 1.0 for most model situations
      if(new_cov_inflate < 1.0) new_cov_inflate = 1.0
      !write(*, *) 'old, new cov_inflate', cov_inflate, new_cov_inflate

      cov_inflate = new_cov_inflate

      !!!if(new_cov_inflate_sd < cov_inflate_sd) cov_inflate_sd = new_cov_inflate_sd
      cov_inflate_sd = new_cov_inflate_sd
      if(cov_inflate_sd < sd_lower_bound) cov_inflate_sd = sd_lower_bound

      ! Ad hoc mechanism for reducing cov_inflate_sd
      !!!cov_inflate_sd = cov_inflate_sd - (cov_inflate_sd - sd_lower_bound) / 10000.0
   endif

   ! Now inflate the ensemble and compute a preliminary inflation increment
   if(deterministic_cov_inflate) then
      ens = (ens_in - prior_mean) * sqrt(cov_inflate) + prior_mean
      inflate_inc = ens - ens_in

   else
      ! Another option at this point would be to add in random noise designed
      ! To increase the variance of the prior ensemble to the appropriate level
      ! Would probably want to keep the mean fixed by shifting AND do a sort
      ! on the final prior/posterior increment pairs to avoid large regression
      ! error as per stocahstic filter algorithms. This might help to avoid 
      ! problems with generating gravity waves in the Bgrid model, for instance.
      ! If this is first time through, need to initialize the random sequence
      if(first_inc_ran_call) then
         call init_random_seq(inc_ran_seq)
         first_inc_ran_call = .false.
      endif

      ! Figure out required sd for random noise being added
      if(cov_inflate > 1.0) then
         ! Don't allow covariance deflation in this version
         rand_sd = sqrt(cov_inflate*prior_var - prior_var)
         !write(*, *) 'rand_sd increment needed is ', sqrt(prior_var), rand_sd
         ! Add random sample from this noise into the ensemble
         do i = 1, ens_size
            ens(i) = random_gaussian(inc_ran_seq, ens_in(i), rand_sd)
         end do
         ! Adjust the mean
         ens = ens - (sum(ens) / ens_size - prior_mean)
      else
         ens = ens_in
      endif
      inflate_inc = ens - ens_in

   endif

else
   ! No covariance inflation is being done, just copy initial ensemble
   ens = ens_in
endif

! Call the appropriate filter option to compute increments for ensemble
if(filter_kind == 1) then
   call     obs_increment_eakf(ens, ens_size, obs, obs_var, obs_inc)
else if(filter_kind == 2) then
   call     obs_increment_enkf(ens, ens_size, obs, obs_var, obs_inc)
else if(filter_kind == 3) then
   call   obs_increment_kernel(ens, ens_size, obs, obs_var, obs_inc)
else if(filter_kind == 4) then
   call obs_increment_particle(ens, ens_size, obs, obs_var, obs_inc)
else 
   call error_handler(E_ERR,'obs_increment', &
              'Illegal value of filter_kind in assim_tools namelist [1-4 OK]', &
              source, revision, revdate)
endif

! Add in the extra increments if doing observation space covariance inflation
if(cov_inflate > 0.0) then
   obs_inc = obs_inc + inflate_inc
endif


! For random algorithm, may need to sort to have minimum increments
! To minimize regression errors, may want to sort to minimize increments
! This makes sense for any of the non-deterministic algorithms
! By doing it here, can take care of both standard non-deterministic updates
! plus non-deterministic obs space covariance inflation. This is expensive, so
! don't use it if it's not needed.
if(sort_obs_inc) then
   new_val = ens_in + obs_inc
   ! Sorting to make increments as small as possible
   call index_sort(ens_in, ens_index, ens_size)
   call index_sort(new_val, new_index, ens_size)
   do i = 1, ens_size
      obs_inc(ens_index(i)) = new_val(new_index(i)) - ens_in(ens_index(i))
! The following erroneous line can provide improved results; understand why
      !!!obs_inc(ens_index(i)) = new_val(new_index(i)) - ens(ens_index(i))
   end do
end if

end subroutine obs_increment



subroutine obs_increment_eakf(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
!
! EAKF version of obs increment

integer, intent(in)   :: ens_size
real(r8), intent(in)  :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: a, prior_mean, new_mean, prior_var, var_ratio, sum_x

! Compute prior variance and mean from sample
sum_x      = sum(ens)
prior_mean = sum_x / ens_size
prior_var  = sum((ens - prior_mean)**2) / (ens_size - 1)

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

a = sqrt(var_ratio)

   obs_inc = a * (ens - prior_mean) + new_mean - ens

end subroutine obs_increment_eakf



subroutine obs_increment_particle(ens, ens_size, obs, obs_var, obs_inc)
!------------------------------------------------------------------------
!
! A observation space only particle filter implementation for a
! two step sequential update filter. Second version, 2 October, 2003.

integer, intent(in)             :: ens_size
real(r8), intent(in)            :: ens(ens_size), obs, obs_var
real(r8), intent(out)           :: obs_inc(ens_size)

real(r8) :: a, weight(ens_size), rel_weight(ens_size), cum_weight(0:ens_size)
real(r8) :: base, frac, new_val(ens_size), weight_sum
integer  :: i, j, indx(ens_size), ens_index(ens_size), new_index(ens_size)

! The factor a is not defined for particle filters
a = -1.0_r8

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

end subroutine obs_increment_particle



subroutine obs_increment_enkf(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment_enkf(ens, ens_size, obs, obs_var, obs_inc)
!

! ENKF version of obs increment

integer, intent(in)             :: ens_size
real(r8), intent(in)            :: ens(ens_size), obs, obs_var
real(r8), intent(out)           :: obs_inc(ens_size)

real(r8) :: a, obs_var_inv
real(r8) :: prior_mean, prior_cov_inv, new_cov, new_mean(ens_size)
real(r8) :: sx, s_x2, prior_cov
real(r8) :: temp_mean, temp_obs(ens_size)
integer  :: ens_index(ens_size), new_index(ens_size)

integer  :: i

! The factor a is not defined for kernel filters
a = -1.0_r8

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
real(r8) :: sx, s_x2
real(r8) :: weight(ens_size), new_mean(ens_size)
real(r8) :: cum_weight, total_weight, cum_frac(ens_size)
real(r8) :: unif, norm, new_member(ens_size)

integer :: i, j, kernel, ens_index(ens_size), new_index(ens_size)

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

end subroutine obs_increment_kernel



subroutine update_from_obs_inc(obs, obs_inc, state, ens_size, &
               state_inc, reg_coef, correl_out)
!========================================================================

! Does linear regression of a state variable onto an observation and
! computes state variable increments from observation increments

integer, intent(in)             :: ens_size
real(r8), intent(in)            :: obs(ens_size), obs_inc(ens_size)
real(r8), intent(in)            :: state(ens_size)
real(r8), intent(out)           :: state_inc(ens_size), reg_coef
real(r8), intent(out), optional :: correl_out

real(r8) :: sum_x, t(ens_size), sum_t2, sum_ty, correl
real(r8) :: sum_y, sum_y2, state_var_norm

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
if(present(correl_out)) then
   sum_y          = sum(state)
   sum_y2         = sum(state*state)
   state_var_norm = sum_y2 - sum_y**2 / ens_size
   correl         = reg_coef * sqrt(sum_t2 / state_var_norm)
endif
if(present(correl_out)) correl_out = correl

! Then compute the increment as product of reg_coef and observation space increment
state_inc = reg_coef * obs_inc

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
real(r8) :: qc(num_obs_in_set), obs_inc(ens_size), a_returned(num_groups), obs_dist
real(r8) :: ens_inc(ens_size), ens_obs(ens_size, num_obs_in_set)
integer, parameter :: first_num_close = 100000
integer :: order(num_obs_in_set)
integer, allocatable :: close_ptr(:, :)
real(r8), allocatable :: dist_ptr(:, :)
real(r8) :: swath(ens_size), reg_factor
type(location_type) :: obs_loc(num_obs_in_set)
logical :: close_to_any, local_close_state(num_obs_in_set)

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

! Need copy of obs to be modified
ens_obs = ens_obs_in

! Set an initial size for the close state pointers
allocate(close_ptr(1, first_num_close), dist_ptr(1, first_num_close))

! Construct an observation temporary
call init_obs(observation, get_num_copies(seq), get_num_qc(seq))

! Initialize search for close states in my domain
if(num_domains == 1) then 
   local_close_state = .true.
else
   local_close_state = .false.
endif

! Get the locations for all of the observations
Locations: do i = 1, num_obs_in_set
   call get_obs_from_key(seq, keys(i), observation)
   call get_obs_def(observation, obs_def)
   obs_loc(i) = get_obs_def_location(obs_def)

   ! Get the value of the observation
   call get_obs_values(observation, obs(i:i), obs_val_index)
   ! Get the qc value set so far; if qc is bad, don't be close to any state
   call get_qc(observation, qc(i:i), 1)
   if(qc(i) > 4.01_r8 .or. obs(i) == missing_r8) cycle Locations

   ! Determine if this observation is close to any state variables in my domain
   if(num_domains > 1) then
555   call get_close_states(seq, keys(i), 2.0_r8*cutoff, num_close_ptr(1), &
      close_ptr(1, :), dist_ptr(1, :), ens(1, :))
      if(num_close_ptr(1) < 0) then
         deallocate(close_ptr, dist_ptr)
         allocate(close_ptr(1, -1 * num_close_ptr(1)), dist_ptr(1, -1 * num_close_ptr(1)))
         goto 555
      endif
      do j = 1, num_close_ptr(1)
         if(my_state(close_ptr(1, j))) then
            local_close_state(i) = .true.
            cycle Locations
         endif
      end do
   endif
end do Locations

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

!---------------------- Sequential filter section --------------------------
! Loop through each observation in the set
Observations : do jjj = 1, num_obs_in_set
   j = order(jjj)
   if(.not. local_close_state(j)) cycle Observations
   
   ! Get this observation in a temporary and get its qc value for possible update
   call get_obs_from_key(seq, keys(j), observation)

   ! Get the observation value and the error variance
   call get_obs_def(observation, obs_def)
   ! Get the value of the observation and the error variance
!!$   call get_obs_values(observation, obs(j:j), obs_val_index)
   obs_err_var(j) = get_obs_def_error_variance(obs_def)
   ! Just skip observations that have failed prior qc ( >4 for now)
   if(qc(j) > 4.01_r8 .or. obs(j) == missing_r8) cycle Observations
   
   ! Compute the ensemble prior for this ob
   if(compute_obs(j)) then
      do k = 1, ens_size
         ! Only compute forward operator if allowed
         call get_expected_obs(seq, keys(j:j), ens(k, :), ens_obs(k:k, j), istatus)
         ! Inability to compute forward operator implies skip this observation
         if (istatus > 0) then
            qc(j) = qc(j) + 4000.0_r8
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
         obs_err_var(j), obs_inc(grp_bot:grp_top))
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
            ens_inc(grp_bot:grp_top), regress(group))
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

      ! Don't need to update if not close to any state variables
      if(.not. local_close_state(k)) cycle

      ! Don't need to do this if qc is bad
      if (qc(k) > 4.01_r8) cycle 

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
            ens_inc(grp_bot:grp_top), regress(group))
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
   333 call set_qc(observation, qc(j:j), 1)
   call set_obs(seq, observation, keys(j))

end do Observations

! Free up storage
deallocate(close_ptr, dist_ptr)
call destroy_obs(observation)

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

integer :: which_domain(model_size), region_size(num_domains)
real(r8) :: new_ens_member(ens_size, model_size), old_ens_member(ens_size, model_size)
type(time_type) :: temp_time

write(*, *) 'STARTING filter_assim'
compute_obs = compute_obs_in

! Set compute obs to false for more than one region (can be modified later)
if(num_domains /= 1) compute_obs = .false.

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

444 continue

! Write out the most up-to-date obs_sequence file, too, which includes qc
if(do_parallel == 2 .or. do_parallel == 3) then
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
! TEMPORARY TEST OF ENSEMBLE TRANSPOSING
do j = 1, num_domains
   domain_bottom = domain_top + 1
   domain_size = base_size
   if(j <= remainder) domain_size = domain_size + 1
   region_size(j) = domain_size
   domain_top = domain_bottom + domain_size - 1
   which_domain(domain_bottom : domain_top) = j
end do

! TRANSPOSE TO REGIONS
call transpose_ens_to_regions(ens_handle, num_domains, which_domain, region_size)



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
   !!!call get_ensemble_region(ens_handle, ens, ens_time, state_vars_in = indices)
   call get_region_by_number(ens_handle, j, domain_size, ens, indices)

! Do ensemble filter update for this region
   if(do_parallel == 0) then
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
         call error_handler(E_ERR,'filter_assim','Use less than 10000 domains.',source,revision,revdate)
      endif

 11   format(a23, i1)
 21   format(a23, i2)
 31   format(a23, i3)
 41   format(a23, i4)

! Development test of passing this interface by file in preparation for separate executables
      iunit = get_unit()
      open(unit = iunit, file = in_file_name(j), action = 'write', form = 'unformatted', &
         status = 'replace')
     ! write(iunit, *) num_domains, j, domain_size, ens_size, model_size, num_groups, &
      write(iunit) num_domains, j, domain_size, ens_size, model_size, num_groups, &
         num_obs_in_set, obs_val_index, confidence_slope, cutoff, save_reg_series, &
         reg_series_unit, obs_sequence_in_name
      !write(iunit, *) ens, ens_obs, compute_obs, keys, my_state
      write(iunit) ens, ens_obs, compute_obs, keys, my_state
      close(iunit)
   endif

   ! Put this region into storage for single executable which is now finished
   !!!if(do_parallel == 0) call put_ensemble_region(ens_handle, ens, ens_time, state_vars_in = indices)
   if(do_parallel == 0) call put_region_by_number(ens_handle, j, domain_size, ens, indices)

   deallocate(ens, indices)

end do

! All the files are out there, now need to spawn off jobs if parallel
if(do_parallel == 2 .or. do_parallel == 3) then

   ! Write out the file names to a control file
   control_unit = get_unit()
   open(unit = control_unit, file = 'assim_region_control')
   write(control_unit, *) num_domains
   do i = 1, num_domains
      write(control_unit, '(a28)' ) in_file_name(i)
      write(control_unit, '(a28)' ) out_file_name(i)
   end do
   close(control_unit)


! Spawn the job to do parallel assims for do_parallel 2
   if(do_parallel == 2) then
   call system('echo go > batchflag; '//parallel_command//' ; sleep 1')

   do
      if(file_exist('batchflag')) then
         call system('sleep 10')
      else
         exit
      endif
   end do
! Do the continuously running script case
   else if(do_parallel == 3) then
      call system('echo a > go_assim_regions')
      do
         if(file_exist('go_assim_regions')) then
            call system('sleep 1')
         else
            exit
         endif
      enddo

   end if

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
      open(unit = iunit, file = out_file_name(j), action = 'read', form = 'unformatted')
      !read(iunit, *) ens
      read(iunit) ens
      close(iunit)
      !!!call put_ensemble_region(ens_handle, ens, ens_time, state_vars_in = indices)
      call put_region_by_number(ens_handle, j, domain_size, ens, indices)
      deallocate(ens, indices)
   end do

endif


! UNTRANSPOSE TO REGIONS
call transpose_regions_to_ens(ens_handle, num_domains, which_domain, region_size)


write(*, *) 'done with subroutine filter_assim cov_inflate is ', cov_inflate

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

! Read in all the arguments from the file
iunit = get_unit()
open(unit = iunit, file = in_file_name, action = 'read', form = 'unformatted')
read(iunit) n_domains, my_domain, domain_size, ens_size, model_size, num_groups, &
   num_obs_in_set, obs_val_index, confidence_slope, cutoff, save_reg_series, &
   reg_series_unit, obs_sequence_in_name

! Allocate storage
allocate(ens(ens_size, domain_size), ens_obs(ens_size, num_obs_in_set), &
   compute_obs(num_obs_in_set), keys(num_obs_in_set), my_state(model_size))
read(iunit) ens, ens_obs, compute_obs, keys, my_state

close(iunit)

! Now need to open the obs sequence file to proceed with call to region assim
obs_sequence_in_name = 'filter_assim_obs_seq'
call read_obs_seq(obs_sequence_in_name, 0, 0, 0, seq2)

! Do ensemble filter update for this region
call filter_assim_region(my_domain, domain_size, ens_size, model_size, num_groups, &
   num_obs_in_set, obs_val_index, confidence_slope, cutoff, save_reg_series, &
   reg_series_unit, ens, ens_obs, compute_obs, seq2, keys, my_state)

! Write out the updated ensemble region, need to do QC also
iunit = get_unit()
open(unit = iunit, file = out_file_name, action = 'write', form = 'unformatted', &
   status = 'replace')
write(iunit) ens
close(iunit)

call destroy_obs_sequence(seq2)

end subroutine async_assim_region




!-----------------------------------------------------------------------

subroutine bayes_cov_inflate(x_p, sigma_p, y_o, sigma_o, l_mean, l_sd, &
   new_cov_inflate, new_cov_inflate_sd)

real(r8), intent(in) :: x_p, sigma_p, y_o, sigma_o, l_mean, l_sd
real(r8), intent(out) :: new_cov_inflate, new_cov_inflate_sd

integer, parameter :: max_num_iters = 100
integer, parameter :: sd_range = 4
integer, parameter :: sd_intervals = 4
real(r8), parameter :: tol = 0.00000001
integer :: i, j

real(r8) :: d(3), lam(3), d_new, lam_new, new_max, e_minus_half, ratio, x_dist
real(r8) :: cov_sd_sum, cov_sd_sample(-sd_intervals:sd_intervals)

! Compute the maximum value of the updated probability
lam(1) = l_mean - 2.0 * l_sd
! Can't go below zero, would get negative undefined inflation
if(lam(1) <= 0.0) lam(1) = 0.0
lam(2) = l_mean
lam(3) = l_mean + 2.0 * l_sd
do i = 1, 3
   d(i) = compute_new_density(x_p, sigma_p, y_o, sigma_o, l_mean, l_sd, lam(i))
end do

if(d(2) < d(1) .or. d(2) < d(3)) then
   write(*, *)'mid-point is not greater than endpoints in bayes_cov_inflate'
   write(*, *) 'need to broaden the search range'
   do i = 1, 3
      write(*, *) 'd ', i, lam(i), d(i)
   end do
new_cov_inflate = l_mean
   !!!stop
endif

! Do a simple 1D optimization by bracketing
do i = 1, max_num_iters
   ! Bring in the right side
   lam_new = lam(1) + (lam(2) - lam(1)) / 2.0
   d_new = compute_new_density(x_p, sigma_p, y_o, sigma_o, l_mean, l_sd, lam_new)
   if(d_new > d(2)) then
      d(3) = d(2); lam(3) = lam(2)
      d(2) = d_new; lam(2) = lam_new
   else
      d(1) = d_new; lam(1) = lam_new
   endif

   ! Now bring in the left side (this can be coded more glamorously, but I'm rushed)
   lam_new = lam(2) + (lam(3) - lam(2)) / 2.0
   d_new = compute_new_density(x_p, sigma_p, y_o, sigma_o, l_mean, l_sd, lam_new)
   if(d_new > d(2)) then
      d(1) = d(2); lam(1) = lam(2)
      d(2) = d_new; lam(2) = lam_new
   else
      d(3) = d_new; lam(3) = lam_new
   endif

   ! Exit if tolerance gets down to 0.000001
   if(abs(lam(1) - lam(2)) < tol) then
      !write(*, *) 'needed ', i, 'iterations'
      exit
   endif
end do

new_cov_inflate = lam(2)
!write(*, *) 'New maximum value is at lambda = ', lam(2)


! Temporarily bail out to save cost when lower bound is reached'
if(l_sd <= sd_lower_bound) then
   new_cov_inflate_sd = l_sd
   return
endif

! Try at a more economical method for computing updated SD
! First compute the new_max value for normalization purposes
new_max = compute_new_density(x_p, sigma_p, y_o, sigma_o, l_mean, l_sd, new_cov_inflate)

! Look at a series of points to deal with skewness???
cov_sd_sum = 0.0
do i = -sd_intervals, sd_intervals
   x_dist = i * sd_range * l_sd / sd_intervals
   ! Next compute the value one old sd out
   lam(1) = new_cov_inflate + x_dist
   if(lam(1) < 0.0) then
      lam(1) = 0.0
      x_dist = new_cov_inflate
   endif
   d(1) = compute_new_density(x_p, sigma_p, y_o, sigma_o, l_mean, l_sd, lam(1))
   ! Compute the ratio of this value to the value at the mean (the max value)
   ratio = d(1) / new_max 
   
   ! Can now compute the standard deviation consistent with this as
   ! sigma = sqrt(-x^2 / (2 ln(r))  where r is ratio and x is l_sd (distance from mean)
   new_cov_inflate_sd = sqrt( -1.0 * x_dist**2 / (2.0 * log(ratio)))
   write(*, *) 'at x new_sd ', lam(1), new_cov_inflate_sd
   if(i /= 0) cov_sd_sum = cov_sd_sum + new_cov_inflate_sd
   cov_sd_sample(i) = new_cov_inflate_sd
   if(i == 0) cov_sd_sample(i) = 1e10
end do

new_cov_inflate_sd = minval(cov_sd_sample)
write(*, *) 'old, new min sd      ', l_sd, new_cov_inflate_sd


end subroutine bayes_cov_inflate

!------------------------------------------------------------------

function compute_new_density(x_p, sigma_p, y_o, sigma_o, l_mean, l_sd, lambda)

real :: compute_new_density
real, intent(in) :: x_p, sigma_p, y_o, sigma_o, l_mean, l_sd, lambda

real :: dist, exponent, l_prob, var, sd, prob

! Compute distance between prior mean and observed value
dist = abs(x_p - y_o)

! Compute probability of this lambda being correct
exponent = (-1.0 * (lambda - l_mean)**2) / (2.0 * l_sd**2)
l_prob = 1.0 / (sqrt(2.0 * 3.14159) * l_sd) * exp(exponent)

! Compute probability that observation would have been observed given this lambda
var = (sigma_p * sqrt(lambda))**2 + sigma_o**2
sd = sqrt(var)

exponent = (-1.0 * dist**2) / (2.0 * var)
prob = 1.0 / (sqrt(2.0 * 3.14159) * sd) * exp(exponent)

! Compute the updated probability density for lambda
compute_new_density = l_prob * prob

end function compute_new_density

!========================================================================
! end module assim_tools_mod
!========================================================================

end module assim_tools_mod

