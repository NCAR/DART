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
use  utilities_mod, only : file_exist, open_file, close_file, check_nml_error, &
                           register_module, error_handler, E_ERR, E_MSG, logfileunit
use       sort_mod, only : index_sort 
use random_seq_mod, only : random_seq_type, random_gaussian, &
                           init_random_seq, random_uniform

implicit none
private

public :: assim_tools_init, obs_increment, update_from_obs_inc, look_for_bias

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

namelist / assim_tools_nml / prior_spread_correction, filter_kind, slope_threshold

!============================================================================

contains


subroutine assim_tools_init()
!============================================================================
! subroutine assim_tools_init()
!

implicit none

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

implicit none

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

implicit none

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

var_ratio  = obs_var / (prior_var + obs_var)
new_mean   = var_ratio * (prior_mean  + prior_var*obs / obs_var)

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

implicit none

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

implicit none

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

implicit none

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

implicit none

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

! The compute the increment as product of reg_coef and observation space increment
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

implicit none

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
! end module assim_tools_mod
!========================================================================

end module assim_tools_mod
