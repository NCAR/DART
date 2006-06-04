! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section, 
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module assim_tools_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 
!
! A variety of operations required by assimilation.

use      types_mod, only : r8, missing_r8, PI
use  utilities_mod, only : file_exist, get_unit, check_namelist_read, find_namelist_in_file, &
                           register_module, error_handler, E_ERR, E_MSG, logfileunit
use       sort_mod, only : index_sort 
use random_seq_mod, only : random_seq_type, random_gaussian, &
                           init_random_seq, random_uniform

use obs_sequence_mod, only : obs_sequence_type, obs_type, get_num_copies, get_num_qc, &
   init_obs, get_obs_from_key, get_obs_def, get_obs_values, get_qc, set_qc, &
   set_obs, get_copy_meta_data, read_obs_seq, destroy_obs_sequence, get_num_obs, &
   write_obs_seq, destroy_obs, get_expected_obs
   
use          obs_def_mod, only : obs_def_type, get_obs_def_error_variance, &
                                 get_obs_def_location, get_obs_def_time
use       cov_cutoff_mod, only : comp_cov_factor
use        obs_model_mod, only : get_close_states
use       reg_factor_mod, only : comp_reg_factor
use         location_mod, only : location_type, get_dist, alloc_get_close_obs, &
                                 get_close_obs
use ensemble_manager_mod, only : ensemble_type, transpose_ens_to_regions, &
                                 transpose_regions_to_ens, put_region_by_number, &
                                 get_region_by_number, is_ens_in_core
use adaptive_inflate_mod, only : adaptive_inflate_obs_init, &
                                 deterministic_inflate, obs_inflate, obs_inflate_sd, &
                                 obs_sd_lower_bound, do_obs_inflate, do_single_ss_inflate, &
                                 do_varying_ss_inflate, &
                                 obs_inf_upper_bound, update_inflation, ss_inf_upper_bound, &
                                 inflate_ens

implicit none
private

public :: assim_tools_init, filter_assim, async_assim_region

type (random_seq_type) :: inc_ran_seq
logical :: first_inc_ran_call = .true.
character(len = 129) :: errstring

logical :: first_get_correction = .true.
real(r8) :: exp_true_correl(200), alpha(200)                                                                      
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
real(r8) :: cutoff = 0.2_r8
logical  :: sort_obs_inc = .false.
integer :: do_parallel = 0
integer :: num_domains = 1
character(len=129) :: parallel_command = './assim_filter.csh'
logical :: spread_restoration = .false.
real(r8) :: internal_outlier_threshold = -1.0_r8
logical :: sampling_error_correction = .false.
integer :: num_close_threshold = 10000000

namelist / assim_tools_nml / filter_kind, cutoff, sort_obs_inc, &
   do_parallel, num_domains, parallel_command, spread_restoration, &
   internal_outlier_threshold, sampling_error_correction, num_close_threshold

!============================================================================

contains

subroutine assim_tools_init(dont_read_restart)
!============================================================================
! subroutine assim_tools_init()
!

logical, intent(in), optional :: dont_read_restart

integer :: iunit, io

call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "assim_tools_nml", iunit)
read(iunit, nml = assim_tools_nml, iostat = io)
call check_namelist_read(iunit, io, "assim_tools_nml")

! Write the namelist values to the log file

call error_handler(E_MSG,'assim_tools_init','assim_tools namelist values',' ',' ',' ')
write(logfileunit, nml=assim_tools_nml)
write(     *     , nml=assim_tools_nml)

! Check for illegal combination of parallel with single region (doesn't make sense)
if(do_parallel /= 0 .and. num_domains == 1) then
   write(errstring, *) 'do_parallel not 0 and num_domains 1 is not supported'
   call error_handler(E_ERR,'assim_tools_init', errstring, source, revision, revdate)
endif

! Initialize the obs_space inflation storage in case it's needed
call adaptive_inflate_obs_init(num_domains, dont_read_restart)

! Look for an error in the do_parallel options
if(do_parallel /= 0 .and. do_parallel /= 2 .and. do_parallel /= 3) then
   write(errstring, *) 'do_parallel option ', do_parallel, ' not supported'
   call error_handler(E_ERR,'assim_tools_init', errstring, source, revision, revdate)
endif

! FOR NOW, can only do spread restoration with filter option 1 (need to extend this)
if(spread_restoration .and. .not. filter_kind == 1) then
   write(errstring, *) 'cant combine spread_restoration and filter_kind ', filter_kind
   call error_handler(E_ERR,'assim_tools_init', errstring, source, revision, revdate)
endif

! AT PRESENT, internal_outlier_threshold only works if adaptive cov_inflate is in use
! Need to modify this in an efficient way in the future
if(internal_outlier_threshold > 0.0_r8 .and. .not. do_obs_inflate) then
   write(errstring, *) 'internal_outlier_threshold checking only works do_obs_inflate true at present'
   call error_handler(E_ERR,'assim_tools_init', errstring, source, revision, revdate)
endif

end subroutine assim_tools_init

!-------------------------------------------------------------

subroutine obs_increment(ens_in, ens_size, obs, obs_var, obs_inc, &
   my_cov_inflate, my_cov_inflate_sd, net_a)


integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens_in(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
real(r8), intent(inout) :: my_cov_inflate, my_cov_inflate_sd
real(r8), intent(out) :: net_a

real(r8) :: ens(ens_size), inflate_inc(ens_size)
real(r8) :: sum_x, prior_mean, prior_var
real(r8) :: new_val(ens_size)
integer :: i, ens_index(ens_size), new_index(ens_size)
real(r8) :: a

! Copy the input ensemble to something that can be modified
ens = ens_in

! If observation space inflation is being done, compute the initial 
! increments and update the inflation factor and its standard deviation
! as needed. my_cov_inflate < 0 means don't do any of this.
if(do_obs_inflate) then
   ! Compute prior variance and mean from sample
   sum_x      = sum(ens)
   prior_mean = sum_x / ens_size
   prior_var = sum((ens - prior_mean)**2) / (ens_size - 1)


   ! If my_cov_inflate_sd is <= 0, just retain current my_cov_inflate setting
   if(my_cov_inflate_sd > 0.0_r8) & 
      ! Gamma set to 1.0 because no distance for observation space
      call update_inflation(my_cov_inflate, my_cov_inflate_sd, prior_mean, prior_var, &
         obs, obs_var, 1.0_r8, obs_sd_lower_bound, obs_inf_upper_bound)

   ! If internal_outlier_threshold is exceeded, don't use this observation
   if(internal_outlier_threshold > 0.0_r8) then
      if(abs(prior_mean - obs) / sqrt(my_cov_inflate * prior_var + obs_var) > &
         internal_outlier_threshold) then
         !!!write(*, *) 'QZ BOUND EXCEEDED: returning'
         obs_inc = 0.0_r8
         net_a = 1.0_r8
         return
      endif
   endif

   ! Now inflate the ensemble and compute a preliminary inflation increment
   call inflate_ens(ens, prior_mean, my_cov_inflate, prior_var)

   ! Keep the increment due to inflation alone 
   inflate_inc = ens - ens_in
endif

! Call the appropriate filter option to compute increments for ensemble
if(filter_kind == 1) then
   call     obs_increment_eakf(ens, ens_size, obs, obs_var, obs_inc, a)
else if(filter_kind == 2) then
   call     obs_increment_enkf(ens, ens_size, obs, obs_var, obs_inc)
else if(filter_kind == 3) then
   call   obs_increment_kernel(ens, ens_size, obs, obs_var, obs_inc)
else if(filter_kind == 4) then
   call obs_increment_particle(ens, ens_size, obs, obs_var, obs_inc)
else if(filter_kind == 5) then
   call obs_increment_ran_kf(ens, ens_size, obs, obs_var, obs_inc)
else if(filter_kind == 6) then
   call obs_increment_det_kf(ens, ens_size, obs, obs_var, obs_inc)
else 
   call error_handler(E_ERR,'obs_increment', &
              'Illegal value of filter_kind in assim_tools namelist [1-6 OK]', &
              source, revision, revdate)
endif

! Add in the extra increments if doing observation space covariance inflation
if(do_obs_inflate) then
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

if(do_obs_inflate) then
   net_a = a * sqrt(my_cov_inflate)
else
   net_a = a
endif

end subroutine obs_increment



subroutine obs_increment_eakf(ens, ens_size, obs, obs_var, obs_inc, a)
!========================================================================
!
! EAKF version of obs increment

integer, intent(in)   :: ens_size
real(r8), intent(in)  :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
real(r8), intent(out) :: a

real(r8) :: prior_mean, new_mean, prior_var, var_ratio, sum_x

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
   endif
endif

a = sqrt(var_ratio)

obs_inc = a * (ens - prior_mean) + new_mean - ens

end subroutine obs_increment_eakf


subroutine obs_increment_ran_kf(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
!
! Forms a random sample of the Gaussian from the update equations.
! This is very close to what a true 'ENSEMBLE' Kalman Filter would 
! look like. Note that outliers, multimodality, etc., get tossed.

integer, intent(in)   :: ens_size
real(r8), intent(in)  :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: prior_mean, new_mean, prior_var, var_ratio, sum_x
real(r8) :: temp_mean, temp_var, new_ens(ens_size), new_var
integer :: i

! Compute prior variance and mean from sample
sum_x      = sum(ens)
prior_mean = sum_x / ens_size
prior_var  = sum((ens - prior_mean)**2) / (ens_size - 1)

if (obs_var /= 0.0_r8) then
   var_ratio = obs_var / (prior_var + obs_var)
   new_var = var_ratio * prior_var
   new_mean  = var_ratio * (prior_mean  + prior_var*obs / obs_var)
else
   if (prior_var /= 0.0_r8) then
      var_ratio = 0.0_r8
      new_var = var_ratio * prior_var
      new_mean  = obs
   else
      call error_handler(E_ERR,'obs_increment_ran_kf', &
           'Both obs_var and prior_var are zero. This is inconsistent', &
           source, revision, revdate)
   endif
endif

! Now, just from a random sample from the updated distribution
! Then adjust the mean (what about adjusting the variance?)!
! Definitely need to sort with this; sort is done in main obs_increment
if(first_inc_ran_call) then
   call init_random_seq(inc_ran_seq)
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



subroutine obs_increment_det_kf(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
!
! Does a deterministic ensemble layout for the updated Gaussian.
! Note that all outliers, multimodal behavior, etc. get tossed.

integer, intent(in)   :: ens_size
real(r8), intent(in)  :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: prior_mean, new_mean, prior_var, var_ratio, sum_x
real(r8) :: temp_var, new_ens(ens_size), new_var
integer :: i

! Compute prior variance and mean from sample
sum_x      = sum(ens)
prior_mean = sum_x / ens_size
prior_var  = sum((ens - prior_mean)**2) / (ens_size - 1)

if (obs_var /= 0.0_r8) then
   var_ratio = obs_var / (prior_var + obs_var)
   new_var = var_ratio * prior_var
   new_mean  = var_ratio * (prior_mean  + prior_var*obs / obs_var)
else
   if (prior_var /= 0.0_r8) then
      var_ratio = 0.0_r8
      new_var = var_ratio * prior_var
      new_mean  = obs
   else
      call error_handler(E_ERR,'obs_increment_det_kf', &
           'Both obs_var and prior_var are zero. This is inconsistent', &
           source, revision, revdate)
   endif
endif

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

integer, intent(in)             :: ens_size
real(r8), intent(in)            :: ens(ens_size), obs, obs_var
real(r8), intent(out)           :: obs_inc(ens_size)

real(r8) :: a, weight(ens_size), rel_weight(ens_size), cum_weight(0:ens_size)
real(r8) :: base, frac, new_val(ens_size), weight_sum
integer  :: i, j, indx(ens_size)

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

! Generate increments
obs_inc = new_val - ens

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

integer  :: i

! The factor a is not defined for kernel filters
a = -1.0_r8

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0_r8 / obs_var

! Compute prior mean and covariance
sx         = sum(ens)
s_x2       = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov  = sum((ens - prior_mean)**2) / (ens_size - 1)

prior_cov_inv = 1.0_r8 / prior_cov
new_cov       = 1.0_r8 / (prior_cov_inv + obs_var_inv)

! Temporary for adjustment test
!   var_ratio = obs_var / (prior_cov + obs_var)
!   updated_mean  = var_ratio * (prior_mean  + prior_cov*obs / obs_var)

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

! Can also adjust mean (and) variance of final sample; works fine
!sx         = sum(new_mean)
!s_x2       = sum(new_mean * new_mean)
!temp_mean = sx / ens_size
!temp_cov  = (s_x2 - sx**2 / ens_size) / (ens_size - 1)
!new_mean = (new_mean - temp_mean) * sqrt(new_cov / temp_cov) + updated_mean
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
real(r8) :: sx, s_x2
real(r8) :: weight(ens_size), new_mean(ens_size)
real(r8) :: cum_weight, total_weight, cum_frac(ens_size)
real(r8) :: unif, norm, new_member(ens_size)

integer :: i, j, kernel

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0_r8 / obs_var

! Compute prior mean and covariance
sx         = sum(ens)
s_x2       = sum(ens * ens)
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

! Generate the increments
obs_inc = new_member - ens

end subroutine obs_increment_kernel



subroutine update_from_obs_inc(obs, obs_inc, state, ens_size, &
               state_inc, reg_coef, correl, net_a)
!========================================================================

! Does linear regression of a state variable onto an observation and
! computes state variable increments from observation increments

integer, intent(in)             :: ens_size
real(r8), intent(in)            :: obs(ens_size), obs_inc(ens_size)
real(r8), intent(in)            :: state(ens_size)
real(r8), intent(out)           :: state_inc(ens_size), reg_coef, correl
real(r8), intent(inout)         :: net_a

real(r8) :: sum_x, t(ens_size), sum_t2, sum_ty, pair(2, ens_size)
real(r8) :: restoration_inc(ens_size), state_mean
real(r8) :: factor
real(r8) :: exp_true_correl, mean_factor

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

! Temporarily need correlation returned for use in adaptive state space
! Computation of correlation                                                         
pair(1, :) = state                                                                   
pair(2, :) = obs                                                                     
call comp_correl(pair, ens_size, correl)        


! BEGIN TEST OF CORRECTION FROM FILE +++++++++++++++++++++++++++++++++++++++++++++++++
! Get the expected actual correlation and the regression weight reduction factor
if(sampling_error_correction) then
   call get_correction_from_file(ens_size, correl, mean_factor, exp_true_correl)
   !write(*, *) correl, exp_true_correl, mean_factor
   ! Watch out for division by zero; if correl is really small regression is safely 0
   if(abs(correl) > 0.001) then
      reg_coef = reg_coef * (exp_true_correl / correl) * mean_factor
   else
      reg_coef = 0.0_r8
   endif
   correl = exp_true_correl
endif

! END TEST OF CORRECTION FROM FILE +++++++++++++++++++++++++++++++++++++++++++++++++



! Then compute the increment as product of reg_coef and observation space increment
state_inc = reg_coef * obs_inc








! Spread restoration algorithm option
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
   ! for now.
   factor = 1.0_r8 / (1.0_r8 - (1.0_r8 - net_a) * (0.8_r8 / ens_size)) - 1.0_r8

   ! Variance restoration
   state_mean = sum(state) / ens_size
   restoration_inc = factor * (state - state_mean)
   state_inc = state_inc + restoration_inc
endif

end subroutine update_from_obs_inc

!========================================================================

subroutine filter_assim_region(domain_size, ens_size, model_size, &
   num_groups, num_obs_in_set, obs_val_index, &
   ens, ens_obs_in, compute_obs, &
   ss_inflate, ss_inflate_sd, ss_sd_lower_bound, seq, keys, my_state, &
   my_cov_inflate, my_cov_inflate_sd)

integer, intent(in) :: domain_size
integer, intent(in) :: ens_size, model_size, num_groups, num_obs_in_set, keys(num_obs_in_set)
real(r8), intent(inout) :: ens(ens_size, domain_size)
real(r8), intent(in) :: ens_obs_in(ens_size, num_obs_in_set)
real(r8), intent(inout) :: ss_inflate(:), ss_inflate_sd(:)
real(r8), intent(in)    ::  ss_sd_lower_bound
logical, intent(in) :: compute_obs(num_obs_in_set)
integer, intent(in) :: obs_val_index
type(obs_sequence_type), intent(inout) :: seq
logical, intent(in) :: my_state(model_size)
real(r8), intent(inout) :: my_cov_inflate, my_cov_inflate_sd

!!!character(len=129) :: msgstring
integer :: i, j, jjj, k, kkk, istatus, ind, grp_size, group, grp_bot, grp_top
integer :: grp_beg(num_groups), grp_end(num_groups)
integer :: num_close_ptr(1), indx
type(obs_type) :: observation
type(obs_def_type) :: obs_def
real(r8) :: obs(num_obs_in_set), obs_err_var(num_obs_in_set), cov_factor, regress(num_groups)
real(r8) :: correl(num_groups)
real(r8) :: qc(num_obs_in_set), obs_inc(ens_size)
real(r8) :: ens_inc(ens_size), ens_obs(ens_size, num_obs_in_set)
integer :: first_num_close
integer :: order(num_obs_in_set)
integer, allocatable :: close_ptr(:, :)
real(r8), allocatable :: dist_ptr(:, :)
real(r8) :: swath(ens_size), reg_factor
type(location_type) :: obs_loc(num_obs_in_set)
real(r8) :: ens_obs_mean(num_obs_in_set, num_groups), ens_obs_var(num_obs_in_set, num_groups)
real(r8) :: gamma, ens_var_deflate

! Added for get_close obs
integer :: obs_box(num_obs_in_set), close_ind(num_obs_in_set), num_close
integer :: num_cant_recompute
real(r8) :: close_dist(num_obs_in_set)
integer :: inv_indices(model_size), inv_index
real(r8) :: net_a(num_groups)
logical :: evaluate_this_ob, assimilate_this_ob
real(r8) :: cutoff_rev

! Specify initial storage size for number of close states
first_num_close = min(domain_size, 400000)

! Divide ensemble into num_groups groups
grp_size = ens_size / num_groups

do group = 1, num_groups
   grp_beg(group) = (group - 1) * grp_size + 1
   grp_end(group) = grp_beg(group) + grp_size - 1
enddo

! Generate an array with the indices of my state variables
inv_indices = 0
indx = 1
do i = 1, model_size
   if(my_state(i)) then
      inv_indices(i) = indx
      indx = indx + 1
   endif
end do

! Need copy of obs to be modified
ens_obs = ens_obs_in

! Need to compute mean and variance of prior obs before doing any sequential assimilation
! For use with adaptive state space inflation algorithms. Should be able to avoid this
! computation if adaptive state space is not in use.
if(do_varying_ss_inflate .or. do_single_ss_inflate) then
   do i = 1, num_obs_in_set
      do j = 1, num_groups
         ens_obs_mean(i, j) = sum(ens_obs(grp_beg(j):grp_end(j), i)) / grp_size
         ens_obs_var(i, j) = &
            sum((ens_obs(grp_beg(j):grp_end(j), i) - ens_obs_mean(i, j))**2) / (grp_size - 1)
      end do
   end do
endif

! Set an initial size for the close state pointers
allocate(close_ptr(1, first_num_close), dist_ptr(1, first_num_close))

! Construct an observation temporary
call init_obs(observation, get_num_copies(seq), get_num_qc(seq))

! Get the locations for all of the observations
Locations: do i = 1, num_obs_in_set
   call get_obs_from_key(seq, keys(i), observation)
   call get_obs_def(observation, obs_def)
   obs_loc(i) = get_obs_def_location(obs_def)
   ! Get the observation error variance
   obs_err_var(i) = get_obs_def_error_variance(obs_def)
   ! Get the value of the observation
   call get_obs_values(observation, obs(i:i), obs_val_index)
   ! Get the qc value set so far
   call get_qc(observation, qc(i:i), 1)
end do Locations

! Reorder the observations to get rid of the ones that can't be recomputed first
indx = 1
do i = 1, num_obs_in_set
   if(.not. compute_obs(i)) then
      order(indx) = i
      indx = indx + 1
   endif
end do

num_cant_recompute = indx - 1

! Then put in all the ones that can be recomputed
do i = 1, num_obs_in_set
   if(compute_obs(i)) then
      order(indx) = i
      indx = indx + 1
   endif
end do

! If there are some that can't be recomputed, need to setup get close states
! ONLY HAVE TO DO THIS FOR THOSE THAT CAN"T BE RECOMPUTED
! ALWAYS NEED THIS FOR OBSERVATIONAL DENSITY THINNING
!!!if(num_cant_recompute > 0) &
   call alloc_get_close_obs(num_obs_in_set, obs_loc, 2.0_r8*cutoff, obs_box)


!---------------------- Sequential filter section --------------------------
! Loop through each observation in the set
Observations : do jjj = 1, num_obs_in_set
   j = order(jjj)

   ! Just skip observations that have failed prior qc ( >3 for now)
   if(qc(j) > 3.01_r8 .or. obs(j) == missing_r8) cycle Observations

   ! Observational density thinning of observations
   ! Get a list of obs that are close to the ob being processed
   call get_close_obs(j, num_obs_in_set, obs_loc, 2.0_r8*cutoff, obs_box, &
      num_close, close_ind, close_dist)

   ! Want expected number of close observations to be reduced to some threshold
   if(num_close > num_close_threshold) then
      ! Change the cutoff radius to get the appropriate number in the circle
      ! This is specific to models on sphere's
      ! Need to get thinning out of assim_tools and into something about locations
      ! Compute a new radius if the num_close is greater than the desired as
      ! 2*cutoff_rev = sqrt(2*cutoff * num_close_threshold / num_close)
      cutoff_rev =  sqrt((2.0_r8*cutoff)**2 * num_close_threshold / num_close) / 2.0_r8
   else 
      cutoff_rev = cutoff
   endif

   
   ! Compute the ensemble prior for this ob if it can be computed here
   if(compute_obs(j)) then
      do k = 1, ens_size
         ! Only compute forward operator if allowed
         call get_expected_obs(seq, keys(j:j), ens(k, :), ens_obs(k:k, j), istatus, &
              assimilate_this_ob, evaluate_this_ob)
         ! If evaluating but not assimilating, just skip but don't change qc
         if(evaluate_this_ob .and. .not. assimilate_this_ob) cycle Observations
         ! Inability to compute forward operator implies skip this observation
         ! Also skip on qc if not being used at all (NOT CLEAR THIS COULD EVER HAPPEN)
         if (istatus > 0 .or. (.not. assimilate_this_ob .and. .not. evaluate_this_ob)) then
            !!!write(msgstring, FMT='(a,i7,a,i7)') 'Skipping obs ',j
            !!!call error_handler(E_MSG,'filter',msgstring,source,revision,revdate)
            qc(j) = qc(j) + 200.0_r8
            call get_obs_from_key(seq, keys(j), observation)
            ! Write out the updated quality control for this observation
            ! At present: THIS DOES NOT GO BACK TO FILE FOR MULTI-DOMAIN RUNS
            call set_qc(observation, qc(j:j), 1)
            call set_obs(seq, observation, keys(j))
            cycle Observations
         endif
      end do
   endif

   ! Getting close states for each scalar observation for now
   ! Need model state for some distance computations in sigma, have
   ! to pick some state, only current ones are ensemble, just pass 1
222 call get_close_states(seq, keys(j), 2.0_r8*cutoff_rev, num_close_ptr(1), &
         close_ptr(1, :), dist_ptr(1, :), ens(1, :))
   if(num_close_ptr(1) < 0) then
      deallocate(close_ptr, dist_ptr)
      allocate(close_ptr(1, -1 * num_close_ptr(1)), dist_ptr(1, -1 * num_close_ptr(1)))
      goto 222
   endif
   
   !!! These calls are useful for WRF but generate too much output for low-order models. 
   !!!write(msgstring, FMT='(a,i7,a,i7)') 'Variables updated for obs ',j,' : ',num_close_ptr(1)
   !!!call error_handler(E_MSG,'filter',msgstring,source,revision,revdate)

   ! Determine if this observation is close to ANY state variable for this process
   CLOSE_STATE1: do k = 1, num_close_ptr(1)
      ind = close_ptr(1, k)
      ! If this state variable is in the list to be updated for this PE, proceed
      if(my_state(ind)) then
         goto 777
      endif
   end do CLOSE_STATE1
   ! Falling off end means no close states, skip this observation
   cycle Observations

   777 continue

   Group1: do group = 1, num_groups
      grp_bot = grp_beg(group)
      grp_top = grp_end(group)

      ! Call obs_increment to do observation space
      call obs_increment(ens_obs(grp_bot:grp_top, j), grp_size, obs(j), &
         obs_err_var(j), obs_inc(grp_bot:grp_top), my_cov_inflate, my_cov_inflate_sd, net_a(group))
   end do Group1
      


   ! Spatially-varying state space inflation
   if(do_single_ss_inflate) then
      do group = 1, num_groups
         if(ss_inflate(1) > 0.0_r8 .and. ss_inflate_sd(1) > 0.0_r8) &
            ! For case with single spatial inflation, use gamma = 1.0_r8
            gamma = 1.0_r8
            ! Deflate the inflated variance
            ! NOTE: Technically, the initial value of the inflation from the start
            ! of the routine should be used instead of ss_inflate(1) which may have
            ! been modified. In generaly, the ss_inflate value changes so little
            ! during one observation set that this doesn't appear worth the cost.
            ens_var_deflate = ens_obs_var(j, group) / &
               (1.0_r8 + gamma*(sqrt(ss_inflate(1)) - 1.0_r8))**2
            call update_inflation(ss_inflate(1), ss_inflate_sd(1), &
               ens_obs_mean(j, group), ens_var_deflate, obs(j), obs_err_var(j), &
               gamma, ss_sd_lower_bound, ss_inf_upper_bound)
      end do
   endif




   ! Now loop through each close state variable for this observation
   CLOSE_STATE: do k = 1, num_close_ptr(1)
      ind = close_ptr(1, k)

      ! If this state variable is not in the list to be updated for this PE, skip
      if(.not. my_state(ind)) cycle CLOSE_STATE

      ! Compute distance dependent envelope
      cov_factor = comp_cov_factor(dist_ptr(1, k), cutoff_rev)
      if(cov_factor <= 0.0_r8) cycle CLOSE_STATE

      ! Get the ensemble elements for this state variable and do regression
      inv_index = inv_indices(ind)
      swath = ens(:, inv_index)

      ! Loop through the groups
      Group2: do group = 1, num_groups
         grp_bot = grp_beg(group)
         grp_top = grp_end(group)
         call update_from_obs_inc(ens_obs(grp_bot:grp_top, j), &
            obs_inc(grp_bot:grp_top), swath(grp_bot:grp_top), grp_size, &
            ens_inc(grp_bot:grp_top), regress(group), correl(group), net_a(group))
      end do Group2

      ! Compute an information factor for impact of this observation on this state
      if(num_groups == 1) then
          reg_factor = 1.0_r8
      else

         ! Pass the time along with the index for possible diagnostic output
         call get_obs_def(observation, obs_def)
         ! Compute regression factor for this obs-state pair
         reg_factor = comp_reg_factor(num_groups, regress, &
            get_obs_def_time(obs_def), j, ind)
      endif

      ! COV_FACTOR NOW COMES IN HERE FOR USE WITH GROUPS; JLA 11/19/03
      reg_factor = min(reg_factor, cov_factor)

      ! Do the final update for this state variable
      ens(:, inv_index) = ens(:, inv_index) + reg_factor * ens_inc(:)

      ! Spatially-varying state space inflation
      if(do_varying_ss_inflate) then
         GroupInflate: do group = 1, num_groups
            if(ss_inflate(inv_index) > 0.0_r8 .and. ss_inflate_sd(inv_index) > 0.0_r8) &
               gamma = reg_factor * abs(correl(group))
               ! Deflate the inflated variance
               ! NOTE: Technically, the initial value of the inflation from the start
               ! of the routine should be used instead of ss_inflate(1) which may have
               ! been modified. In generaly, the ss_inflate value changes so little
               ! during one observation set that this doesn't appear worth the cost.
               ens_var_deflate = ens_obs_var(j, group) / &
                  (1.0_r8 + gamma*(sqrt(ss_inflate(inv_index)) - 1.0_r8))**2
               call update_inflation(ss_inflate(inv_index), ss_inflate_sd(inv_index), &
                  ens_obs_mean(j, group), ens_var_deflate, obs(j), obs_err_var(j), &
                  gamma, ss_sd_lower_bound, ss_inf_upper_bound)
         end do GroupInflate
      endif

   end do CLOSE_STATE

   ! Also need to update all obs variable priors that cannot be recomputed
   ! If there are no obs that can't be recomputed can skip obs update section
   if(num_cant_recompute == 0) cycle Observations

   ! GET_CLOSE_OBS; get a list of obs that are close to the ob being processed
   ! Done up front for observational density thinning
   !!!call get_close_obs(j, num_obs_in_set, obs_loc, 2.0_r8*cutoff, obs_box, &
   !!!   num_close, close_ind, close_dist)

   ! Loop through the list of obs close to this one
   UPDATE_OBS: do kkk = 1, num_close 
      ! Get the index of the next close observation
      k = order(close_ind(kkk))

      ! Obs already used will never be used again, be careful not to do the jth one here!
      if(k <= j) cycle UPDATE_OBS

      ! Don't need to do this if qc is bad
      if (qc(k) > 3.01_r8) cycle UPDATE_OBS

      ! Don't need to do this if observation can be recomputed
      if(compute_obs(k)) cycle UPDATE_OBS
      ! Compute distance dependent envelope
      cov_factor = comp_cov_factor(close_dist(kkk), cutoff_rev)
      if(cov_factor <= 0.0_r8) cycle UPDATE_OBS

      ! Get the ensemble elements for this state variable and do regression
      swath = ens_obs(:, k)

      ! Loop through the groups
      Group3: do group = 1, num_groups
         grp_bot = grp_beg(group)
         grp_top = grp_end(group)
         call update_from_obs_inc(ens_obs(grp_bot:grp_top, j), &
            obs_inc(grp_bot:grp_top), swath(grp_bot:grp_top), grp_size, &
            ens_inc(grp_bot:grp_top), regress(group), correl(group), net_a(group))
      end do Group3

      ! Compute an information factor for impact of this observation on this state
      if(num_groups == 1) then
          reg_factor = 1.0_r8
      else
         ! Pass the time along with the index for possible diagnostic output
         call get_obs_def(observation, obs_def)
         ! Compute regression factor for this obs-state pair
         ! Negative final argument indicates that this is an obs-obs regression computation.
         reg_factor = comp_reg_factor(num_groups, regress, &
            get_obs_def_time(obs_def), j, -1*k)
      endif

      ! COV_FACTOR NOW COMES IN HERE FOR USE WITH GROUPS; JLA 11/19/03
      reg_factor = min(reg_factor, cov_factor)

      ! Do the final update for this obs. variable
      ens_obs(:, k) = ens_obs(:, k) + reg_factor * ens_inc(:)

   end do UPDATE_OBS

!------------------

end do Observations

! Free up storage
deallocate(close_ptr, dist_ptr)
call destroy_obs(observation)

end subroutine filter_assim_region


!---------------------------------------------------------------------------

subroutine filter_assim(ens_handle, ens_obs, compute_obs_in, &
   ens_size, model_size, num_obs_in_set, num_groups, ss_inflate, ss_inflate_sd, &
   ss_sd_lower_bound, seq, keys, obs_sequence_in_name)

integer,                 intent(in) :: ens_size, model_size, num_groups, num_obs_in_set
integer,                 intent(in) :: keys(num_obs_in_set)
type(ensemble_type),     intent(inout) :: ens_handle
real(r8),                intent(in) :: ens_obs(ens_size, num_obs_in_set)
logical,                 intent(in) :: compute_obs_in(num_obs_in_set)
real(r8),                intent(inout) :: ss_inflate(:), ss_inflate_sd(:)
real(r8),                intent(in) :: ss_sd_lower_bound
type(obs_sequence_type), intent(inout) :: seq
character(len=*),        intent(in) :: obs_sequence_in_name

integer :: i, j, indx, obs_val_index, iunit, control_unit
integer :: base_size, remainder, domain_top, domain_bottom
logical :: compute_obs(num_obs_in_set), my_state(model_size)

real(r8), allocatable :: ens(:,:), my_ss_inf(:), my_ss_inf_sd(:)
integer,  allocatable :: indices(:)

character(len=28) :: in_file_name(num_domains), out_file_name(num_domains)
character(len=129) :: temp_obs_seq_file, errstring

integer :: which_domain(model_size), region_size(num_domains)

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
   region_size(j) = base_size
   if(j <= remainder) region_size(j) = region_size(j) + 1
   domain_top = domain_bottom + region_size(j) - 1
   which_domain(domain_bottom : domain_top) = j
end do

! TRANSPOSE TO REGIONS
call transpose_ens_to_regions(ens_handle, num_domains, which_domain, region_size)



domain_top = 0
do j = 1, num_domains
   domain_bottom = domain_top + 1

   ! Find out which state variables are in my domain
   !call whats_in_my_domain(num_domains, j, my_state, region_size(j))
   ! Temporarily just do the domain stuff here
   domain_top = domain_bottom + region_size(j) - 1
   my_state = .false.
   my_state(domain_bottom : domain_top) = .true.

   if(num_domains > 1 .or. .not. is_ens_in_core(ens_handle)) then
      ! Generate an array with the indices of my state variables
      ! Only need full size array for my_ss_inf's IF spatially varying ss inflate
      if(do_varying_ss_inflate) then
         allocate(ens(ens_size, region_size(j)), indices(region_size(j)), &
            my_ss_inf(region_size(j)), my_ss_inf_sd(region_size(j)))
      else 
         allocate(ens(ens_size, region_size(j)), indices(region_size(j)), &
            my_ss_inf(1), my_ss_inf_sd(1))
      endif

      ! For single inflate, copy the single value of the inflate and sd
      if(do_single_ss_inflate) then
         my_ss_inf(1) = ss_inflate(1)
         my_ss_inf_sd(1) = ss_inflate_sd(1)
      endif

      indx = 1
      do i = 1, model_size
         if(my_state(i)) then
            indices(indx) = i
            ! If spatially varying inflate, need to copy over all values for domain
            if(do_varying_ss_inflate) then
               my_ss_inf(indx) = ss_inflate(i)
               my_ss_inf_sd(indx) = ss_inflate_sd(i)
            endif
            indx = indx + 1
         endif
      end do

      ! Get the ensemble (whole thing for now) from storage
      call get_region_by_number(ens_handle, j, region_size(j), ens, indices)
   endif

   ! Do ensemble filter update for this region
   if(do_parallel == 0) then
      if(num_domains > 1 .or. .not. is_ens_in_core(ens_handle)) then
         call filter_assim_region(region_size(j), ens_size, model_size, num_groups, &
            num_obs_in_set, obs_val_index, &
            ens, ens_obs, compute_obs, &
            my_ss_inf, my_ss_inf_sd, ss_sd_lower_bound, seq, keys, my_state, &
            obs_inflate(j), obs_inflate_sd(j))
         ! Copy the inflation update back to standard vector
         ! Only copy for last domain for single state space
         if(do_single_ss_inflate .and. j == num_domains) then
            ss_inflate(1) = my_ss_inf(1)
            ss_inflate_sd(1) = my_ss_inf_sd(1)
         endif

         do i = 1, region_size(j)
            indx = indices(i)
            if(do_varying_ss_inflate) then
               ss_inflate(indx) = my_ss_inf(i)
               ss_inflate_sd(indx) = my_ss_inf_sd(i)
               !!!??? Efficiency if single domain...
            endif
         end do

      else
! Single domain, in core for direct access to ensemble storage
         call filter_assim_region(region_size(j), ens_size, model_size, num_groups, &
            num_obs_in_set, obs_val_index, &
            ens_handle%ens, ens_obs, compute_obs, &
            ss_inflate, ss_inflate_sd, ss_sd_lower_bound, seq, keys, my_state, &
            obs_inflate(j), obs_inflate_sd(j))
      endif
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

! Passing this interface by file in preparation for separate executables
      iunit = get_unit()
      open(unit = iunit, file = in_file_name(j), action = 'write', form = 'unformatted', &
      !!!open(unit = iunit, file = in_file_name(j), action = 'write', form = 'formatted', &
         status = 'replace')
      write(iunit) num_domains, region_size(j), ens_size, model_size, num_groups, &
         num_obs_in_set, obs_val_index, obs_sequence_in_name
      write(iunit) ens, ens_obs, compute_obs, my_ss_inf, my_ss_inf_sd, ss_sd_lower_bound
      write(iunit) keys, my_state, obs_inflate(j), obs_inflate_sd(j)
      close(iunit)
   endif

   ! Put this region into storage for single executable which is now finished
   if(do_parallel == 0 .and. (num_domains > 1 .or. .not. is_ens_in_core(ens_handle))) & 
      call put_region_by_number(ens_handle, j, region_size(j), ens, indices)
   ! Free up the storage allocated for the region
   if(num_domains > 1 .or. .not. is_ens_in_core(ens_handle)) &
      deallocate(ens, indices, my_ss_inf, my_ss_inf_sd)

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

   domain_top = 0
   do j = 1, num_domains
   ! Find out which state variables are in my domain
   !call whats_in_my_domain(num_domains, j, my_state, region_size(j))
   ! Temporarily just do the domain stuff here
      domain_bottom = domain_top + 1

      domain_top = domain_bottom + region_size(j) - 1
      my_state = .false.
      my_state(domain_bottom : domain_top) = .true.

      ! Generate an array with the indices of my state variables
      ! Only need full size array for my_ss_inf's IF spatially varying ss inflate
      if(do_varying_ss_inflate) then
         allocate(ens(ens_size, region_size(j)), indices(region_size(j)), &
            my_ss_inf(region_size(j)), my_ss_inf_sd(region_size(j)))
      else 
         allocate(ens(ens_size, region_size(j)), indices(region_size(j)), &
            my_ss_inf(1), my_ss_inf_sd(1))
      endif


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
      !!!open(unit = iunit, file = out_file_name(j), action = 'read', form = 'formatted')
      read(iunit) ens, my_ss_inf, my_ss_inf_sd, obs_inflate(j), obs_inflate_sd(j)
      close(iunit)
  
      ! Only one value to copy back if single state space inflation
      ! Only copy single state space if this is the last domain
      if(do_single_ss_inflate .and. j == num_domains) then
         ss_inflate(1) = my_ss_inf(1)
         ss_inflate_sd(1) = my_ss_inf_sd(1)
      endif
      

      ! Copy the inflation update back to standard vector for varying state space inflate
      if(do_varying_ss_inflate) then
         do i = 1, region_size(j)
            indx = indices(i)
            ss_inflate(indx) = my_ss_inf(i)
            ss_inflate_sd(indx) = my_ss_inf_sd(i)
         end do
      endif

      call put_region_by_number(ens_handle, j, region_size(j), ens, indices)
      deallocate(ens, indices, my_ss_inf, my_ss_inf_sd)
   end do

endif


! UNTRANSPOSE TO REGIONS
call transpose_regions_to_ens(ens_handle, num_domains, which_domain, region_size)

end subroutine filter_assim



!-----------------------------------------------------------------------

subroutine async_assim_region(in_file_name, out_file_name)

character(len = *), intent(in) :: in_file_name, out_file_name

! DONT FORGET THE QC!!!
! AND IN GENERAL NEED TO WORK ON REGSERIES STUFF!!!

integer :: n_domains, domain_size, ens_size, model_size, num_groups
integer :: num_obs_in_set, obs_val_index
character(len = 129) :: obs_sequence_in_name
real(r8), allocatable :: ens(:, :), ens_obs(:, :), my_ss_inf(:), my_ss_inf_sd(:)
integer, allocatable :: keys(:)
logical, allocatable :: my_state(:), compute_obs(:)
real(r8) :: my_cov_inflate, my_cov_inflate_sd, ss_sd_lower_bound

type(obs_sequence_type) :: seq2
integer :: iunit

! Read in all the arguments from the file
iunit = get_unit()
open(unit = iunit, file = in_file_name, action = 'read', form = 'unformatted')
!!!open(unit = iunit, file = in_file_name, action = 'read', form = 'formatted')
read(iunit) n_domains, domain_size, ens_size, model_size, num_groups, &
   num_obs_in_set, obs_val_index, obs_sequence_in_name

! Allocate storage; size of state space inflate is model_size only if spatially varying
if(do_varying_ss_inflate) then
   allocate(ens(ens_size, domain_size), ens_obs(ens_size, num_obs_in_set), &
      compute_obs(num_obs_in_set), keys(num_obs_in_set), my_state(model_size), &
      my_ss_inf(domain_size), my_ss_inf_sd(domain_size))
else
   allocate(ens(ens_size, domain_size), ens_obs(ens_size, num_obs_in_set), &
      compute_obs(num_obs_in_set), keys(num_obs_in_set), my_state(model_size), &
      my_ss_inf(1), my_ss_inf_sd(1))
endif

read(iunit) ens, ens_obs, compute_obs, my_ss_inf, my_ss_inf_sd, ss_sd_lower_bound
read(iunit) keys, my_state, my_cov_inflate, my_cov_inflate_sd

close(iunit)

! Now need to open the obs sequence file to proceed with call to region assim
obs_sequence_in_name = 'filter_assim_obs_seq'
call read_obs_seq(obs_sequence_in_name, 0, 0, 0, seq2)

! Do ensemble filter update for this region
call filter_assim_region(domain_size, ens_size, model_size, num_groups, &
   num_obs_in_set, obs_val_index, &
   ens, ens_obs, compute_obs, my_ss_inf, my_ss_inf_sd, ss_sd_lower_bound, &
   seq2, keys, my_state, my_cov_inflate, my_cov_inflate_sd)

! Write out the updated ensemble region, need to do QC also
iunit = get_unit()
open(unit = iunit, file = out_file_name, action = 'write', form = 'unformatted', &
!!!open(unit = iunit, file = out_file_name, action = 'write', form = 'formatted', &
   status = 'replace')
write(iunit) ens, my_ss_inf, my_ss_inf_sd, my_cov_inflate, my_cov_inflate_sd
close(iunit)

call destroy_obs_sequence(seq2)

end subroutine async_assim_region


!------------------------------------------------------------------


subroutine comp_correl(ens, n, correl)

implicit none

integer, intent(in) :: n
real(r8), intent(in) :: ens(2, n)
real(r8), intent(out) :: correl
real(r8) :: sum_x, sum_y, sum_xy, sum_x2, sum_y2, denom_2

! There is an efficiency issue with the regression being computed
! independently from this. Should be addressed.

sum_x = sum(ens(2, :))
sum_y = sum(ens(1, :))
sum_xy = sum(ens(2, :) * ens(1, :))
sum_x2 = sum(ens(2, :) * ens(2, :))

! Computation of correlation
sum_y2 = sum(ens(1, :) * ens(1, :))

! Avoid overflow for no variance and illegal values from round-off
denom_2 = (n * sum_x2 - sum_x**2) * (n * sum_y2 - sum_y**2)
if(denom_2 <= 0.0_r8) then
   correl = 0.0_r8
else
   correl = (n * sum_xy - sum_x * sum_y) / sqrt(denom_2)
   if(correl >  1.0_r8) correl =  1.0_r8
   if(correl < -1.0_r8) correl = -1.0_r8
endif

end subroutine comp_correl


!------------------------------------------------------------------------

subroutine get_correction_from_file(ens_size, scorrel, mean_factor, expected_true_correl)

integer, intent(in) :: ens_size
real(r8), intent(in) :: scorrel
real(r8), intent(out) :: mean_factor, expected_true_correl

! Reads in a regression error file for a give ensemble_size and uses interpolation
! to get correction factor into the file

integer :: iunit, i
real(r8) :: temp, temp2, correl
real(r8) :: fract, low_correl, low_exp_correl, low_alpha
real(r8) :: high_correl, high_exp_correl, high_alpha
integer :: low_indx, high_indx
character(len = 20) :: correction_file_name
character(len = 129) :: errstring

if(first_get_correction) then
   ! Compute the file name for this ensemble size
   if(ens_size < 10) then
      write(correction_file_name, 11) 'final_full.', ens_size
   else if(ens_size < 100) then
      write(correction_file_name, 21) 'final_full.', ens_size
   else if(ens_size < 1000) then
      write(correction_file_name, 31) 'final_full.', ens_size
   else if(ens_size < 10000) then
      write(correction_file_name, 41) 'final_full.', ens_size
   else
      write(errstring,*)'Trying to use ',ens_size,' model states -- too many.'
      call error_handler(E_MSG,'get_correction_from_file',errstring,source,revision,revdate)
      call error_handler(E_ERR,'get_correction_from_file','Use less than 10000 ensemble.',source,revision,revdate)

    11   format(a11, i1)
    21   format(a11, i2)
    31   format(a11, i3)
    41   format(a11, i4)
   endif
 
   ! Make sure that the correction file exists, else an error
   if(.not. file_exist(correction_file_name)) then
      write(errstring,*) 'Correction file ', correction_file_name, ' does not exist'
      call error_handler(E_ERR,'get_correction_from_file',errstring,source,revision,revdate)
   endif

   ! Read in file to get the expected value of the true correlation given the sample
   iunit = get_unit()
   open(unit = iunit, file = correction_file_name)
   do i = 1, 200
      read(iunit, *) temp, temp2, exp_true_correl(i), alpha(i)
   end do
   close(iunit)

   first_get_correction = .false.
endif


! First quick modification of correl to expected true correl for test (should interp)
if(scorrel < -1.0_r8) then
   correl = -1.0_r8
   mean_factor = 1.0_r8
else if(scorrel > 1.0_r8) then
   correl = 1.0_r8
   mean_factor = 1.0_r8
else if(scorrel <= -0.995_r8) then
   fract = (scorrel + 1.0_r8) / 0.05_r8
   correl = (exp_true_correl(1) + 1.0_r8) * fract - 1.0_r8
   mean_factor = (alpha(1) - 1.0_r8) * fract + 1.0_r8
else if(scorrel >= 0.995_r8) then
   fract = (scorrel - 0.995_r8) / 0.05_r8
   correl = (1.0_r8 - exp_true_correl(200)) * fract + exp_true_correl(200)
   mean_factor = (1.0_r8 - alpha(200)) * fract + alpha(200)
else
   low_indx = (scorrel + 0.995_r8) / 0.01_r8 + 1.0_r8
   low_correl = -0.995_r8 + (low_indx - 1) * 0.01
   low_exp_correl = exp_true_correl(low_indx)
   low_alpha = alpha(low_indx)
   high_indx = low_indx + 1
   high_correl = low_correl + 0.01
   high_exp_correl = exp_true_correl(high_indx)
   high_alpha = alpha(low_indx)
   fract = (scorrel - low_correl) / (high_correl - low_correl)
   correl = (high_exp_correl - low_exp_correl) * fract + low_exp_correl
   mean_factor = (high_alpha - low_alpha) * fract + low_alpha
endif

expected_true_correl = correl

end subroutine get_correction_from_file

!------------------------------------------------------------------------




!========================================================================
! end module assim_tools_mod
!========================================================================

end module assim_tools_mod
