! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$ 

!> Operations and storage required for various adaptive inflation algorithms

module adaptive_inflate_mod

!> \defgroup adaptive_inflate adaptive_inflate_mod
!> @{

use types_mod,            only : r8, PI, missing_r8
use time_manager_mod,     only : time_type, get_time
use utilities_mod,        only : register_module, open_file, close_file, &
                                 error_handler, E_ERR, E_MSG
use random_seq_mod,       only : random_seq_type, random_gaussian, init_random_seq
use ensemble_manager_mod, only : ensemble_type,   &
                                 get_copy_owner_index, prepare_to_write_to_vars,                &
                                 prepare_to_read_from_vars, prepare_to_update_vars, &
                                 map_pe_to_task

use mpi_utilities_mod,    only : my_task_id, send_to, receive_from, send_minmax_to

implicit none
private

public :: update_inflation,                                          do_obs_inflate,     &
          do_varying_ss_inflate,      do_single_ss_inflate,          inflate_ens,        &
          adaptive_inflate_init,      adaptive_inflate_type,                             &
                                      deterministic_inflate,         solve_quadratic,    &
          log_inflation_info,         get_minmax_task_zero,          mean_from_restart,  &
          sd_from_restart,                                                               &
          output_inf_restart,         get_inflate_mean,              get_inflate_sd,     &
          get_is_prior,               get_is_posterior,              do_ss_inflate,      &
          set_inflation_mean_copy,    set_inflation_sd_copy,         get_inflation_mean_copy, &
          get_inflation_sd_copy,      do_rtps_inflate


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! Manages both observation space and state space inflation
! Handles initial values and restarts, diagnostic output, and computations
! Algorithm options at present include a single fixed observation space,
! a single fixed state space adaptive inflation,
! and a spatially-varying state space inflation that carries
! a mean and variance for the state space inflation at each point. 

! Type to keep track of information for inflation
type adaptive_inflate_type
   private
   ! Flavor can be 0:none, 1:obs_inflate, 2:varying_ss_inflate, 3:single_ss_inflate
   ! Deprecating 1:obs_inflate, there is concerns how the observation space inflation
   ! is happening. JPH.
   integer               :: inflation_flavor
   logical               :: output_restart = .false.
   logical               :: deterministic
   real(r8)              :: inflate, sd, sd_lower_bound, inf_lower_bound, inf_upper_bound
   ! Include a random sequence type in case non-deterministic inflation is used
   type(random_seq_type) :: ran_seq
   logical               :: allow_missing_in_clm
   real(r8)              :: minmax_mean(2), minmax_sd(2)
   logical               :: mean_from_restart
   logical               :: sd_from_restart
   logical               :: prior = .false.
   logical               :: posterior = .false.
   integer               :: input_mean_copy = -1 !>todo NO_COPY_PRESENT
   integer               :: input_sd_copy   = -1
end type adaptive_inflate_type

! Module storage for writing error messages
character(len=512) :: msgstring, msgstring2

! Flag indicating whether module has been initialized
logical :: initialized = .false.

!============================================================================

contains

!------------------------------------------------------------------
! Accessor functions for adaptive inflate type
function mean_from_restart(inflation)

type(adaptive_inflate_type) :: inflation
logical :: mean_from_restart

mean_from_restart = inflation%mean_from_restart

end function mean_from_restart

!------------------------------------------------------------------

function sd_from_restart(inflation)

type(adaptive_inflate_type) :: inflation
logical :: sd_from_restart

sd_from_restart = inflation%sd_from_restart

end function sd_from_restart

!------------------------------------------------------------------

function output_inf_restart(inflation)

type(adaptive_inflate_type) :: inflation
logical :: output_inf_restart

output_inf_restart = inflation%output_restart

end function

!------------------------------------------------------------------
function get_inflate_mean(inflation)

type(adaptive_inflate_type) :: inflation
real(r8)  :: get_inflate_mean

get_inflate_mean = inflation%inflate

end function

!------------------------------------------------------------------
function get_inflate_sd(inflation)

type(adaptive_inflate_type) :: inflation
real(r8)  :: get_inflate_sd

get_inflate_sd = inflation%sd

end function

!------------------------------------------------------------------
function get_is_prior(inflation)

type(adaptive_inflate_type) :: inflation
logical :: get_is_prior

get_is_prior = inflation%prior

end function get_is_prior

!------------------------------------------------------------------
function get_is_posterior(inflation)

type(adaptive_inflate_type) :: inflation
logical :: get_is_posterior

get_is_posterior = inflation%posterior

end function get_is_posterior

!------------------------------------------------------------------
! Check for state_space_inflation
function do_ss_inflate(inflation)

type(adaptive_inflate_type), intent(in) :: inflation
logical :: do_ss_inflate

if (do_single_ss_inflate(inflation) .or. do_varying_ss_inflate(inflation) .or. &
    do_rtps_inflate(inflation)) then
   do_ss_inflate = .true.
else
   do_ss_inflate = .false.
endif

end function do_ss_inflate
!------------------------------------------------------------------

subroutine adaptive_inflate_init(inflate_handle, inf_flavor, mean_from_restart, &
   sd_from_restart, output_inflation, deterministic, & 
   inf_initial, sd_initial, inf_lower_bound, inf_upper_bound, &
   sd_lower_bound, ens_handle, missing_ok, label)

! Initializes an adaptive_inflate_type 

type(adaptive_inflate_type), intent(inout) :: inflate_handle
integer,                     intent(in)    :: inf_flavor
logical,                     intent(in)    :: mean_from_restart
logical,                     intent(in)    :: sd_from_restart
logical,                     intent(in)    :: output_inflation
logical,                     intent(in)    :: deterministic
real(r8),                    intent(in)    :: inf_initial, sd_initial
real(r8),                    intent(in)    :: inf_lower_bound, inf_upper_bound
real(r8),                    intent(in)    :: sd_lower_bound
type(ensemble_type),         intent(inout) :: ens_handle
logical,                     intent(in)    :: missing_ok
character(len = *),          intent(in)    :: label

! Record the module version if this is first initialize call
if(.not. initialized) then
   initialized = .true.
   call register_module(source, revision, revdate)
endif

! If non-deterministic inflation is being done, need to initialize random sequence.
! use the task id number (plus 1 since they start at 0) to set the initial seed.
! NOTE: non-deterministic inflation does NOT reproduce as process count is varied!
if(.not. deterministic) then
   call init_random_seq(inflate_handle%ran_seq, my_task_id()+1)
endif

! Load up the structure first to keep track of all details of this inflation type
inflate_handle%inflation_flavor   = inf_flavor
inflate_handle%output_restart     = output_inflation
inflate_handle%deterministic      = deterministic
inflate_handle%inflate            = inf_initial
inflate_handle%sd                 = sd_initial
inflate_handle%inf_lower_bound    = inf_lower_bound
inflate_handle%inf_upper_bound    = inf_upper_bound
inflate_handle%sd_lower_bound     = sd_lower_bound
inflate_handle%allow_missing_in_clm = missing_ok
inflate_handle%mean_from_restart  = mean_from_restart
inflate_handle%sd_from_restart    = sd_from_restart
! Prior and posterior are intialized to false
if (trim(label)=='Prior') inflate_handle%prior = .true.
if (trim(label)=='Posterior') inflate_handle%posterior = .true.

! Cannot support non-determistic inflation and an inf_lower_bound < 1
if(.not. deterministic .and. inf_lower_bound < 1.0_r8) then
   write(msgstring, *) 'Cannot have non-deterministic inflation and inf_lower_bound < 1'
   call error_handler(E_ERR, 'adaptive_inflate_init', msgstring, source, revision, revdate)
endif

! give these distinctive values; if inflation is being used
! (e.g. inf_flavor > 0) then they should be set in all cases.
inflate_handle%minmax_mean(:) = missing_r8
inflate_handle%minmax_sd(:)   = missing_r8

! State space inflation is read in the IO routine read_state.

! Read type 1 (observation space inflation)
if(inf_flavor == 1) then

   write(msgstring,  *) 'No longer supporting observation space inflation ', &
                        '(i.e. inf_flavor = 1).'
   write(msgstring2, *) 'Please contact dart@ucar.edu if you would like to use ', &
                        'observation space inflation'
   call error_handler(E_ERR, 'adaptive_inflate_init', &
      msgstring, source, revision, revdate, text2=msgstring2)

endif

end subroutine adaptive_inflate_init

!------------------------------------------------------------------

function do_obs_inflate(inflate_handle)

! Returns true if this inflation type indicates observation space inflation

logical                                 :: do_obs_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

do_obs_inflate = (inflate_handle%inflation_flavor == 1)

if (do_obs_inflate) then
  write(msgstring,  *) 'observation space inflation not suppported (i.e. inf_flavor = 1)'
  write(msgstring2, *) 'please contact dart if you would like to use this functionality'
  call error_handler(E_ERR, 'do_obs_inflate', &
     msgstring, source, revision, revdate, text2=msgstring2)
endif

end function do_obs_inflate

!------------------------------------------------------------------

function do_varying_ss_inflate(inflate_handle)

! Returns true if this inflation type indicates varying state space inflation

logical                                 :: do_varying_ss_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

do_varying_ss_inflate = (inflate_handle%inflation_flavor == 2)

end function do_varying_ss_inflate

!------------------------------------------------------------------

function do_single_ss_inflate(inflate_handle)

! Returns true if this inflation type indicates fixed state space inflation

logical                                 :: do_single_ss_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

do_single_ss_inflate = (inflate_handle%inflation_flavor == 3)

end function do_single_ss_inflate

!------------------------------------------------------------------

function do_rtps_inflate(inflate_handle)

! Returns true if this inflation type indicates posterior relaxion-to-prior-spread
! (whitaker & Hamill, 2012)

logical                                 :: do_rtps_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

do_rtps_inflate = (inflate_handle%inflation_flavor == 4)

end function do_rtps_inflate

!------------------------------------------------------------------

function deterministic_inflate(inflate_handle)

! Returns true if deterministic inflation is indicated

logical :: deterministic_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

deterministic_inflate = inflate_handle%deterministic

end function deterministic_inflate

!------------------------------------------------------------------

subroutine inflate_ens(inflate_handle, ens, mean, inflate, var_in, fsprd, asprd)

! Inflates subset of ensemble members given mean and inflate
! Selects between deterministic and stochastic inflation

type(adaptive_inflate_type), intent(inout) :: inflate_handle
real(r8),                    intent(inout) :: ens(:)
real(r8),                    intent(in)    :: mean, inflate
real(r8), optional,          intent(in)    :: var_in
real(r8), optional,          intent(in)    :: fsprd, asprd

integer  :: i, ens_size
real(r8) :: rand_sd, var, sd_inflate

! it's possible to have MISSING_R8s in the state
! vector now.  so we need to be able to avoid changing
! MISSING_R8 values by inflation here.
if (inflate_handle%allow_missing_in_clm) then
   if (any(ens == MISSING_R8)) return
endif

if(inflate_handle%deterministic) then

   if ( do_rtps_inflate(inflate_handle)) then
      if ( .not. present(fsprd) .or. .not. present(asprd)) then 
         write(msgstring, *) 'missing arguments for RTPS inflation, should not happen'
         call error_handler(E_ERR,'inflate_ens',msgstring,source,revision,revdate) 
      endif 
      ! only inflate if spreads are > 0
      if ( asprd .gt. 0.0_r8 .and. fsprd .gt. 0.0_r8) &
          ens = mean + (ens-mean) * ( inflate*((fsprd-asprd)/asprd) + 1.0_r8 )
   else 

      ! Spread the ensemble out linearly for deterministic
      ! Following line can lead to inflation of 1.0 changing ens on some compilers
      !!! ens = (ens - mean) * sqrt(inflate) + mean
      ! Following gives 1.0 inflation having no impact on known compilers
      sd_inflate = sqrt(inflate) 
      ens = ens * sd_inflate + mean * (1.0_r8 - sd_inflate)

   endif

else
   ! Use a stochastic algorithm to spread out.
   ens_size = size(ens)

   ! If var is not present, go ahead and compute it here.
   if(.not. present(var_in)) then
      var = sum((ens - mean)**2) / (ens_size - 1)
   else
      var = var_in
   endif
   
   ! To increase the variance of the prior ensemble to the appropriate level
   ! probably want to keep the mean fixed by shifting AND do a sort
   ! on the final prior/posterior increment pairs to avoid large regression
   ! error as per stochastic filter algorithms. This might help to avoid
   ! problems with generating gravity waves in the Bgrid model, for instance.

   ! The following code does not do the sort.
   
   ! Figure out required sd for random noise being added
   ! Don't allow covariance deflation in this version
   if(inflate > 1.0_r8) then
      rand_sd = sqrt(inflate*var - var)
      ! Add random sample from this noise into the ensemble
      do i = 1, ens_size
         ens(i) = random_gaussian(inflate_handle%ran_seq, ens(i), rand_sd)
      end do
      ! Adjust the mean back to the original value
      ens = ens - (sum(ens) / ens_size - mean)
   endif
endif

end subroutine inflate_ens

!------------------------------------------------------------------

subroutine update_inflation(inflate_handle, inflate, inflate_sd, prior_mean, prior_var, &
   obs, obs_var, gamma)

! Given information from an inflate type, scalar values for inflate and inflate_sd,
! the ensemble prior_mean and prior_var for an observation, and the obsered value
! and observational error variance, computes updated values for the inflate and
! inflate_sd values using the algorithms documented on the DART website.
! The gamma paramter gives the localized prior correlation times the localization
! which is computed in the assim_tools routine filter_assim. For single state
! space inflation it is 1.0.

type(adaptive_inflate_type), intent(in)    :: inflate_handle
real(r8),                    intent(inout) :: inflate, inflate_sd
real(r8),                    intent(in)    :: prior_mean, prior_var, obs, obs_var, gamma

real(r8) :: new_inflate, new_inflate_sd

! If the inflate_sd not positive, keep everything the same
if(inflate_sd <= 0.0_r8) return

! A lower bound on the updated inflation sd and an upper bound
! on the inflation itself are provided in the inflate_handle. 

! Use bayes theorem to update
call bayes_cov_inflate(prior_mean, prior_var, obs, obs_var, inflate, &
   inflate_sd, gamma, new_inflate, new_inflate_sd, inflate_handle%sd_lower_bound)

! Make sure inflate satisfies constraints
inflate = new_inflate
if(inflate < inflate_handle%inf_lower_bound) inflate = inflate_handle%inf_lower_bound
if(inflate > inflate_handle%inf_upper_bound) inflate = inflate_handle%inf_upper_bound

! Make sure sd satisfies constraints
inflate_sd = new_inflate_sd
if(inflate_sd < inflate_handle%sd_lower_bound) inflate_sd = inflate_handle%sd_lower_bound

end subroutine update_inflation

!------------------------------------------------------------------

subroutine bayes_cov_inflate(x_p, sigma_p_2, y_o, sigma_o_2, lambda_mean, lambda_sd, &
   gamma, new_cov_inflate, new_cov_inflate_sd, sd_lower_bound_in)

! Uses algorithms in references on DART web site to update the distribution of inflation.

real(r8), intent(in)  :: x_p, sigma_p_2, y_o, sigma_o_2, lambda_mean, lambda_sd, gamma
real(r8), intent(in)  :: sd_lower_bound_in
real(r8), intent(out) :: new_cov_inflate, new_cov_inflate_sd

integer  :: i, mlambda_index(1)

real(r8) :: new_1_sd, new_max, ratio, lambda_sd_2
real(r8) :: dist_2, b, c, d, Q, R, disc, alpha, beta, cube_root_alpha, cube_root_beta, x
real(r8) :: rrr, cube_root_rrr, angle, mx(3), sep(3), mlambda(3)


! If gamma is 0, nothing happens
if(gamma <= 0.0_r8) then
   new_cov_inflate = lambda_mean
   new_cov_inflate_sd = lambda_sd
   return
endif

! Computation saver
lambda_sd_2 = lambda_sd**2
dist_2 = (x_p - y_o)**2

! Use ONLY the linear approximation, cubic solution below can be numerically
! unstable for extreme cases. Should look at this later.
!!!if(gamma > 0.99_r8) then
if(gamma > 1.01_r8) then

! FIXME: rectify the comment in the following line (about gamma
! and 1.0) and the test above, which clearly used to be true for 
! 1.0 but is not anymore.

! The solution of the cubic below only works if gamma is 1.0
! Can analytically find the maximum of the product: d/dlambda is a
! cubic polynomial in lambda**2; solve using cubic formula for real root
! Can write so that coefficient of x**3 is 1, other coefficients are:
   b = -1.0_r8 * (sigma_o_2 + sigma_p_2 * lambda_mean)
   c = lambda_sd_2 * sigma_p_2**2 / 2.0_r8
   d = -1.0_r8 * (lambda_sd_2 * sigma_p_2**2 * dist_2) / 2.0_r8

   Q = c - b**2 / 3
   R = d + (2 * b**3) / 27 - (b * c) / 3

   ! Compute discriminant, if this is negative have 3 real roots, else 1 real root
   disc = R**2 / 4 + Q**3 / 27

   if(disc < 0.0_r8) then
      rrr = sqrt(-1.0 * Q**3 / 27)
      ! Note that rrr is positive so no problem for cube root
      cube_root_rrr = rrr ** (1.0 / 3.0)
      angle = acos(-0.5 * R / rrr)
      do i = 0, 2
         mx(i+1) = 2.0_r8 * cube_root_rrr * cos((angle + i * 2.0_r8 * PI) / 3.0_r8) - b / 3.0_r8
         mlambda(i + 1) = (mx(i + 1) - sigma_o_2) / sigma_p_2
         sep(i+1) = abs(mlambda(i + 1) - lambda_mean)
      end do
      ! Root closest to initial peak is appropriate
      mlambda_index = minloc(sep)
      new_cov_inflate = mlambda(mlambda_index(1))

   else
      ! Only one real root here, find it.

      ! Compute the two primary terms
      alpha = -R/2 + sqrt(disc)
      beta = R/2 + sqrt(disc)

      cube_root_alpha = abs(alpha) ** (1.0 / 3.0) * abs(alpha) / alpha
      cube_root_beta = abs(beta) ** (1.0 / 3.0) * abs(beta) / beta
   
      x = cube_root_alpha - cube_root_beta - b / 3.0

      ! This root is the value of x = theta**2
      new_cov_inflate = (x - sigma_o_2) / sigma_p_2

   endif

   ! Put in code to approximate the mode (new_cov_inflate)
   !write(*, *) 'old, orig mode is ', lambda_mean, new_cov_inflate
else
   ! Approximate with Taylor series for likelihood term
   call linear_bayes(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd_2, gamma, &
      new_cov_inflate)
endif

! Bail out to save cost when lower bound is reached on lambda standard deviation
if(lambda_sd <= sd_lower_bound_in) then
   new_cov_inflate_sd = lambda_sd
else
   ! Compute by forcing a Gaussian fit at one positive SD
! First compute the new_max value for normalization purposes
   new_max = compute_new_density(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd, &
      gamma, new_cov_inflate)

! Find value at a point one OLD sd above new mean value
   new_1_sd = compute_new_density(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd, gamma, &
      new_cov_inflate + lambda_sd)

   ! If either the numerator or denominator of the following computation 
   ! of 'ratio' is going to be zero (or almost so), return the original incoming
   ! inflation value.  The computation would have resulted in either Inf or NaN.
   if (abs(new_max) <= TINY(0.0_r8) .or. abs(new_1_sd) <= TINY(0.0_r8)) then
      new_cov_inflate_sd = lambda_sd
      return
   endif

   ratio = new_1_sd / new_max 

   ! Another error for numerical issues; if ratio is larger than 0.99, bail out
   if(ratio > 0.99) then
      new_cov_inflate_sd = lambda_sd
      return
   endif

   ! Can now compute the standard deviation consistent with this as
      ! sigma = sqrt(-x^2 / (2 ln(r))  where r is ratio and x is lambda_sd (distance from mean)
   new_cov_inflate_sd = sqrt( -1.0_r8 * lambda_sd_2 / (2.0_r8 * log(ratio)))

   ! Prevent an increase in the sd of lambda???
   ! For now, this is mostly countering numerical errors in this computation
   if(new_cov_inflate_sd > lambda_sd) new_cov_inflate_sd = lambda_sd

endif

end subroutine bayes_cov_inflate

!------------------------------------------------------------------

function compute_new_density(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd, gamma, lambda)

! Used to update density by taking approximate gaussian product

real(r8)             :: compute_new_density
real(r8), intent(in) :: dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd, gamma, lambda

real(r8) :: theta_2, theta
real(r8) :: exponent_prior, exponent_likelihood


! Compute probability of this lambda being correct
exponent_prior = (lambda - lambda_mean)**2 / (-2.0_r8 * lambda_sd**2)

! Compute probability that observation would have been observed given this lambda
theta_2 = (1.0_r8 + gamma * (sqrt(lambda) - 1.0_r8))**2 * sigma_p_2 + sigma_o_2
theta = sqrt(theta_2)

exponent_likelihood = dist_2 / ( -2.0_r8 * theta_2)

! Compute the updated probability density for lambda
! Have 1 / sqrt(2 PI) twice, so product is 1 / (2 PI)
compute_new_density = exp(exponent_likelihood + exponent_prior) / &
   (2.0_r8 * PI * lambda_sd * theta)

end function compute_new_density


!---------------------------------------------------------------------

subroutine linear_bayes(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd_2, gamma, &
   new_cov_inflate)

real(r8), intent(in)    :: dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd_2
real(r8), intent(in)    :: gamma
real(r8), intent(inout) :: new_cov_inflate

real(r8) :: theta_bar_2, u_bar, like_exp_bar, v_bar, like_bar, like_prime, theta_bar
real(r8) :: a, b, c, plus_root, minus_root, dtheta_dlambda
   
! Compute value of theta at current lambda_mean
theta_bar_2 = (1.0_r8 + gamma * (sqrt(lambda_mean) - 1.0_r8))**2 * sigma_p_2 + sigma_o_2
theta_bar = sqrt(theta_bar_2)
! Compute constant coefficient for likelihood at lambda_bar
u_bar = 1.0_r8 / (sqrt(2.0_r8 * PI) * theta_bar)
! Compute exponent of likelihood at lambda_bar
like_exp_bar = dist_2 / (-2.0_r8 * theta_bar_2)
! Compute exponential part of likelihood at lambda_bar
v_bar = exp(like_exp_bar)
! Compute value of likelihood at current lambda_bar value
like_bar = u_bar * v_bar

! If like_bar goes to 0, can't do anything, so just keep current values
if(like_bar <= 0.0_r8) then
   new_cov_inflate = lambda_mean
   return
endif

! Next compute derivative of likelihood at this point

! First compute d/dlambda of theta evaluated at lambda_mean
! Verified correct by finite difference, 1 January, 2006
dtheta_dlambda = 0.5_r8 * sigma_p_2 * gamma *(1.0_r8 - gamma + gamma*sqrt(lambda_mean)) / &
   (theta_bar * sqrt(lambda_mean))
like_prime = (u_bar * v_bar * dtheta_dlambda / theta_bar) * (dist_2 / theta_bar_2 - 1.0_r8)

! If like_prime goes to 0, can't do anything, so just keep current values
if(like_prime == 0.0_r8) then
   new_cov_inflate = lambda_mean
   return
endif

a = 1.0_r8
b = like_bar / like_prime - 2.0_r8 * lambda_mean
c = lambda_mean**2 -lambda_sd_2 - like_bar * lambda_mean / like_prime

! Use nice scaled quadratic solver to avoid precision issues
call solve_quadratic(a, b, c, plus_root, minus_root)

! Do a check to pick closest root
if(abs(minus_root - lambda_mean) < abs(plus_root - lambda_mean)) then
   new_cov_inflate = minus_root
else
   new_cov_inflate = plus_root
endif

end subroutine linear_bayes

!------------------------------------------------------------------------

subroutine solve_quadratic(a, b, c, r1, r2)

real(r8), intent(in)  :: a, b, c
real(r8), intent(out) :: r1, r2

real(r8) :: scaling, as, bs, cs, disc

! Scale the coefficients to get better round-off tolerance
scaling = max(abs(a), abs(b), abs(c))
as = a / scaling
bs = b / scaling
cs = c / scaling

! Get discriminant of scaled equation
disc = sqrt(bs**2 - 4.0_r8 * as * cs)

if(bs > 0.0_r8) then
   r1 = (-bs - disc) / (2 * as)
else
   r1 = (-bs + disc) / (2 * as)
endif

! Compute the second root given the larger one
r2 = (cs / as) / r1

end subroutine solve_quadratic


!------------------------------------------------------------------------
!> Write to log file what kind of inflation is being used.  
subroutine log_inflation_info(inflation_handle, mype, label)

type(adaptive_inflate_type), intent(in) :: inflation_handle
integer,                     intent(in) :: mype
character(len = *),          intent(in) :: label

character(len = 128) :: det, tadapt, sadapt, akind

if(inflation_handle%deterministic) then
  det = 'deterministic,'
else
  det = 'random-noise,'
endif
if (inflation_handle%minmax_sd(2) > inflation_handle%sd_lower_bound) then
   det = trim(det) // ' variance adaptive,'
endif
if (inflation_handle%inf_lower_bound < 1.0_r8) then
   det = trim(det) // ' deflation permitted,'
endif
if (inflation_handle%minmax_sd(2) > 0.0_r8) then
  tadapt = ' time-adaptive,'
else
  tadapt = ' time-constant,'
endif

select case(inflation_handle%inflation_flavor)
   case (0)
      det = ''
      tadapt = ''
      sadapt = ''
      akind = 'None '
   case (1)
      sadapt = ''
      akind = ' observation-space'
   case (2)
      sadapt = ' spatially-varying,'
      akind = ' state-space '
   case (3)
      sadapt = ' spatially-constant,'
      akind = ' state-space'
   case (4)
      sadapt = ' spatially-varying relaxation-to-prior-spread,'
      akind = ' state-space'
   case default
      write(msgstring, *) 'Illegal inflation value for ', label
      call error_handler(E_ERR, 'adaptive_inflate_init', msgstring, source, revision, revdate)
end select

! say in plain english what kind of inflation was selected.
write(msgstring, '(4A)') &
   trim(det), trim(tadapt), trim(sadapt), trim(akind)
call error_handler(E_MSG, trim(label) // ' inflation:', msgstring, source, revision, revdate)

! if inflation was set from a namelist, that message has already been
! printed (further up in this routine).  for values set from a restart
! file, if the inflation flavor is 2, the values printed are the min and
! max from the entire array.  for flavors 1 and 3 there is only a single
! value to print out. 
if (inflation_handle%inflation_flavor > 0) then
   if (inflation_handle%mean_from_restart .and. mype == 0) then
      if (inflation_handle%inflation_flavor == 2) then
         write(msgstring, '(A, F8.3, A, F8.3)') &
            'inf mean   from restart file: min value: ', inflation_handle%minmax_mean(1), ' max value: ', inflation_handle%minmax_mean(2)
      else
         write(msgstring, '(A, F8.3, A, F8.3)') &
            'inf mean   from restart file: value: ', inflation_handle%minmax_mean(1)
      endif
      call error_handler(E_MSG, trim(label) // ' inflation:', msgstring, source, revision, revdate)
   endif
   if (inflation_handle%sd_from_restart .and. mype == 0) then
      if (inflation_handle%inflation_flavor == 2) then
         write(msgstring, '(A, F8.3, A, F8.3)') &
            'inf stddev from restart file: min value: ', inflation_handle%minmax_sd(1), ' max value: ', inflation_handle%minmax_sd(2)
      else
         write(msgstring, '(A, F8.3, A, F8.3)') &
            'inf stddev from restart file: value: ', inflation_handle%minmax_sd(1)
      endif
      call error_handler(E_MSG, trim(label) // ' inflation:', msgstring, source, revision, revdate)
   endif
endif

end subroutine log_inflation_info


!-----------------------------------------------------------------------
! Collect the min and max of inflation on task 0
! this block figures out what the min/max value of the mean/sd is
! if we are reading in the values from a restart file.  it is used
! in diagnostic output so it needs to get to PE0.  we also could check it
! against the limits in the namelist to be sure the file values aren't
! already outside the requested limits. (not sure that's necessary -
! depends on when the code that changes the values imposes the limits.)
subroutine get_minmax_task_zero(inflation_handle, ens_handle, ss_inflate_index, ss_inflate_sd_index)

type(adaptive_inflate_type), intent(inout) :: inflation_handle
type(ensemble_type),         intent(inout) :: ens_handle
integer,                     intent(in)    :: ss_inflate_index
integer,                     intent(in)    :: ss_inflate_sd_index

real(r8) :: minmax_mean(2), minmax_sd(2), global_val(2)

if (inflation_handle%inflation_flavor >= 0) then
if (inflation_handle%mean_from_restart) then

   ! find min and max on each processor
   minmax_mean(1) = minval(ens_handle%copies(ss_inflate_index, :))
   minmax_mean(2) = maxval(ens_handle%copies(ss_inflate_index, :))

   ! collect on pe 0
   call send_minmax_to(minmax_mean, map_pe_to_task(ens_handle, 0), global_val)
   if (ens_handle%my_pe == 0) inflation_handle%minmax_mean = global_val

endif

if (inflation_handle%sd_from_restart) then

   ! find min and max on each processor
   minmax_sd(1) = minval(ens_handle%copies(ss_inflate_sd_index, :))
   minmax_sd(2) = maxval(ens_handle%copies(ss_inflate_sd_index, :))

   ! collect on task 0
   call send_minmax_to(minmax_sd, map_pe_to_task(ens_handle, 0), global_val)
   if (ens_handle%my_pe == 0) inflation_handle%minmax_sd = minmax_sd

endif
endif

end subroutine get_minmax_task_zero

!-----------------------------------------------------------------------

subroutine set_inflation_mean_copy(inflation_handle, c)
type(adaptive_inflate_type), intent(inout) :: inflation_handle
integer,                     intent(in)    :: c

inflation_handle%input_mean_copy = c

end subroutine set_inflation_mean_copy

!-----------------------------------------------------------------------

subroutine set_inflation_sd_copy(inflation_handle, c)
type(adaptive_inflate_type), intent(inout) :: inflation_handle
integer,                     intent(in)    :: c

inflation_handle%input_sd_copy = c

end subroutine set_inflation_sd_copy

!-----------------------------------------------------------------------

function get_inflation_mean_copy(inflation_handle) result (c)
type(adaptive_inflate_type), intent(in) :: inflation_handle
integer :: c

c = inflation_handle%input_mean_copy

end function get_inflation_mean_copy

!-----------------------------------------------------------------------

function get_inflation_sd_copy(inflation_handle) result (c)
type(adaptive_inflate_type), intent(in) :: inflation_handle
integer :: c

c = inflation_handle%input_sd_copy

end function get_inflation_sd_copy


!========================================================================
! end module adaptive_inflate_mod
!========================================================================

!> @}

end module adaptive_inflate_mod

! <next few lines under version control, do not edit>
! $URL$ 
! $Id$ 
! $Revision$ 
! $Date$ 
