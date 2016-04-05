! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module adaptive_inflate_mod

! <next few lines under version control, do not edit>
! $URL$ 
! $Id$ 
! $Revision$ 
! $Date$ 
!
! Operations and storage required for various adaptive inflation algorithms

use types_mod,            only : r8, PI
use time_manager_mod,     only : time_type, get_time, set_time
use utilities_mod,        only : file_exist, get_unit, register_module, &
                                 error_handler, E_ERR, E_MSG, logfileunit
use random_seq_mod,       only : random_seq_type, random_gaussian, init_random_seq, &
                                 random_uniform
use ensemble_manager_mod, only : ensemble_type, all_vars_to_all_copies, all_copies_to_all_vars, &
                                 read_ensemble_restart, write_ensemble_restart

implicit none
private

public :: update_inflation,           adaptive_inflate_end,          do_obs_inflate,     &
          do_varying_ss_inflate,      do_single_ss_inflate,          inflate_ens,        &
          adaptive_inflate_init,      adaptive_inflate_type,         get_inflate,        &
          get_sd,                     set_inflate,                   set_sd,             &
          output_inflate_diagnostics, deterministic_inflate


! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

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
   integer               :: inflation_flavor, obs_diag_unit
   logical               :: start_from_restart, output_restart, deterministic
   character(len = 129)  :: in_file_name, out_file_name, diag_file_name
   real(r8)              :: inflate, sd, sd_lower_bound, inf_lower_bound, inf_upper_bound
   ! Include a random sequence type in case non-deterministic inflation is used
   type(random_seq_type) :: ran_seq
end type adaptive_inflate_type

! Module storage for writing error messages
character(len = 129) :: errstring

! Flag indicating whether module has been initialized
logical :: initialized = .false.

!============================================================================

contains

!------------------------------------------------------------------

subroutine adaptive_inflate_init(inflate_handle, inf_flavor, start_from_restart, &
   output_restart, deterministic, in_file_name, out_file_name, diag_file_name, &
   inf_initial, sd_initial, inf_lower_bound, inf_upper_bound, sd_lower_bound, &
   ens_handle, ss_inflate_index, ss_inflate_sd_index)

! Initializes an adaptive_inflate_type 

type(adaptive_inflate_type), intent(inout) :: inflate_handle
integer,                     intent(in)    :: inf_flavor
logical,                     intent(in)    :: start_from_restart
logical,                     intent(in)    :: output_restart
logical,                     intent(in)    :: deterministic
character(len = *),          intent(in)    :: in_file_name
character(len = *),          intent(in)    :: out_file_name
character(len = *),          intent(in)    :: diag_file_name
real(r8),                    intent(in)    :: inf_initial, sd_initial
real(r8),                    intent(in)    :: inf_lower_bound, inf_upper_bound
real(r8),                    intent(in)    ::  sd_lower_bound
type(ensemble_type),         intent(inout) :: ens_handle
integer,                     intent(in)    :: ss_inflate_index, ss_inflate_sd_index

integer :: restart_unit

! Load up the structure first to keep track of all details of this inflation type
inflate_handle%inflation_flavor   = inf_flavor
inflate_handle%start_from_restart = start_from_restart
inflate_handle%output_restart     = output_restart
inflate_handle%deterministic      = deterministic
inflate_handle%in_file_name       = in_file_name
inflate_handle%out_file_name      = out_file_name
inflate_handle%diag_file_name     = diag_file_name
inflate_handle%inflate            = inf_initial
inflate_handle%sd                 = sd_initial
inflate_handle%inf_lower_bound    = inf_lower_bound
inflate_handle%inf_upper_bound    = inf_upper_bound
inflate_handle%sd_lower_bound     = sd_lower_bound

! Set obs_diag unit to -1 indicating it has not been opened yet
inflate_handle%obs_diag_unit = -1

! Record the module version if this is first initialize call
if(.not. initialized) then
   initialized = .true.
   call register_module(source, revision, revdate)
endif

! If non-deterministic inflation is being done, need to initialize random sequence
! NOTE: non-deterministic inflation does NOT reproduce as process count is varied!
if(.not. deterministic) call init_random_seq(inflate_handle%ran_seq)

! Cannot support non-determistic inflation and an inf_lower_bound < 1
if(.not. deterministic .and. inf_lower_bound < 1.0_r8) then
   write(errstring, *) 'Cannot have non-deterministic inflation and inf_lower_bound < 1'
   call error_handler(E_ERR, 'adaptive_inflate_init', errstring, source, revision, revdate)
endif

!------ Block for state space inflation initialization ------

! Types 2 and 3 are state space inflation types
if(inf_flavor >= 2) then
   ! Initialize state space inflation, copies in ensemble are given
   ! by the inflate and inflate_sd indices. These should be contiguous.

   ! Verify that indices are contiguous
   if(ss_inflate_sd_index /= ss_inflate_index + 1) then
      write(errstring, *) 'ss_inflate_index = ', ss_inflate_index, &
         ' and ss_inflate_sd_index = ', ss_inflate_sd_index, ' must be continguous'
      call error_handler(E_ERR, 'adaptive_inflate_init', &
         errstring, source, revision, revdate)
   endif

   ! Read in initial values from file OR get from subroutine arguments
   if(start_from_restart) then
      call read_ensemble_restart(ens_handle, ss_inflate_index, ss_inflate_sd_index, &
         start_from_restart, in_file_name, force_single_file = .true.)
   else
      ! Get initial values from higher level; requires pe's to have all copies of some vars
      call all_vars_to_all_copies(ens_handle)
      ens_handle%copies(ss_inflate_index, :)    = inf_initial
      ens_handle%copies(ss_inflate_sd_index, :) =  sd_initial
      call all_copies_to_all_vars(ens_handle)
   endif

!------ Block for obs. space inflation initialization ------

! Type 1 is observation space inflation
else if(inf_flavor == 1) then

   ! Initialize observation space inflation values from restart files
   ! Only values are inflation, inflation_sd
   if(start_from_restart) then
      ! Open the file
      restart_unit = get_unit()
      open(unit = restart_unit, file = in_file_name, action = 'read', form = 'formatted')
      read(restart_unit, *) inflate_handle%inflate, inflate_handle%sd
      close(restart_unit)
   endif

endif

end subroutine adaptive_inflate_init

!------------------------------------------------------------------

subroutine adaptive_inflate_end(inflate_handle, ens_handle, ss_inflate_index, &
   ss_inflate_sd_index)

type(adaptive_inflate_type), intent(in)    :: inflate_handle
type(ensemble_type),         intent(inout) :: ens_handle
integer,                     intent(in)    :: ss_inflate_index, ss_inflate_sd_index

integer :: restart_unit

if(inflate_handle%output_restart) then
   ! Use the ensemble manager to output restart for state space (flavors 2 or 3)
   if(do_varying_ss_inflate(inflate_handle) .or. do_single_ss_inflate(inflate_handle)) then
      ! Verify that indices are contiguous
      if(ss_inflate_sd_index /= ss_inflate_index + 1) then
         write(errstring, *) 'ss_inflate_index = ', ss_inflate_index, &
            ' and ss_inflate_sd_index = ', ss_inflate_sd_index, ' must be continguous'
         call error_handler(E_ERR, 'adaptive_inflate_end', &
            errstring, source, revision, revdate)
      endif

      ! Write the inflate and inflate_sd as two copies for a restart
      call write_ensemble_restart(ens_handle, inflate_handle%out_file_name, &
         ss_inflate_index, ss_inflate_sd_index, force_single_file = .true.)

   ! Flavor 1 is observation space, write its restart directly
   else if(do_obs_inflate(inflate_handle)) then
      ! Open the restart file
      restart_unit = get_unit()
      open(unit = restart_unit, file = inflate_handle%out_file_name, &
         action = 'write', form = 'formatted')
      write(restart_unit, *) inflate_handle%inflate, inflate_handle%sd
      close(unit = restart_unit)
   endif
endif

! Need to close diagnostic files for observation space if in use
if(inflate_handle%obs_diag_unit > -1) close(inflate_handle%obs_diag_unit)
   
end subroutine adaptive_inflate_end


!------------------------------------------------------------------

function do_obs_inflate(inflate_handle)

! Returns true if this inflation type indicates observation space inflation

logical                                 :: do_obs_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

do_obs_inflate = (inflate_handle%inflation_flavor == 1)

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

function deterministic_inflate(inflate_handle)

! Returns true if deterministic inflation is indicated

logical :: deterministic_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

deterministic_inflate = inflate_handle%deterministic

end function deterministic_inflate

!------------------------------------------------------------------

function get_inflate(inflate_handle)

! The single real value inflate contains the obs_space inflation value
! when obs_space inflation is in use and this retrieves it.

real(r8)                                :: get_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

get_inflate = inflate_handle%inflate

end function get_inflate

!------------------------------------------------------------------

function get_sd(inflate_handle)

! The single real value inflate_sd contains the obs_space inflate_sd value
! when obs_space inflation is in use and this retrieves it.

real(r8)                                :: get_sd
type(adaptive_inflate_type), intent(in) :: inflate_handle

get_sd = inflate_handle%sd

end function get_sd

!------------------------------------------------------------------

subroutine set_inflate(inflate_handle, inflate)

! Sets the single inflation value in the type

type(adaptive_inflate_type), intent(inout) :: inflate_handle
real(r8),                    intent(in)    :: inflate

inflate_handle%inflate = inflate

end subroutine set_inflate

!------------------------------------------------------------------

subroutine set_sd(inflate_handle, sd)

! Sets the single sd value in the type

type(adaptive_inflate_type), intent(inout) :: inflate_handle
real(r8),                    intent(in)    :: sd

inflate_handle%sd = sd

end subroutine set_sd

!------------------------------------------------------------------

subroutine output_inflate_diagnostics(inflate_handle, time)

type(adaptive_inflate_type), intent(inout) :: inflate_handle
type(time_type),             intent(in)    :: time

integer :: days, seconds

! Diagnostics for state space inflate are done by the filter on the state-space
! netcdf diagnostic files. Here, need to do initial naive ascii dump for obs space.
! Values can come from storage in this module directly.

! Only need to do something if obs_space
if(do_obs_inflate(inflate_handle)) then
   ! If unit is -1, it hasn't been opened yet, do it.
   if(inflate_handle%obs_diag_unit == -1) then
      ! Open the file
      inflate_handle%obs_diag_unit = get_unit()
      open(unit = inflate_handle%obs_diag_unit, file = inflate_handle%diag_file_name, &
         action = 'write', form = 'formatted')
   endif

   ! Get the time in days and seconds
   call get_time(time, seconds, days)
   ! Write out the time followed by the values
   write(inflate_handle%obs_diag_unit, *) days, seconds, inflate_handle%inflate, &
      inflate_handle%sd
endif

end subroutine output_inflate_diagnostics

!------------------------------------------------------------------

subroutine inflate_ens(inflate_handle, ens, mean, inflate, var_in)

! Inflates subset of ensemble members given mean and inflate
! Selects between deterministic and stochastic inflation

type(adaptive_inflate_type), intent(inout) :: inflate_handle
real(r8),                    intent(inout) :: ens(:)
real(r8),                    intent(in)    :: mean, inflate
real(r8), optional,          intent(in)    :: var_in

integer  :: i, ens_size
real(r8) :: rand_sd, var, sd_inflate

if(inflate_handle%deterministic) then
   ! Just spread the ensemble out linearly for deterministic
   ! Following line can lead to inflation of 1.0 changing ens on some compilers
   !!! ens = (ens - mean) * sqrt(inflate) + mean
   ! Following gives 1.0 inflation having no impact on known compilers
   sd_inflate = sqrt(inflate) 
   ens = ens * sd_inflate + mean * (1.0_r8 - sd_inflate)

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

! If the inflate_sd is negative, just keep everything the same
if(inflate_sd < 0.0_r8) return

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

subroutine comp_likelihood(dist_2, sigma_p_2, sigma_o_2, lambda_mean, gamma, like_bar)
real(r8), intent(in)  :: dist_2, sigma_p_2, sigma_o_2, lambda_mean, gamma
real(r8), intent(out) :: like_bar

real :: theta_bar_2, u_bar, like_exp_bar, v_bar

! Compute value of theta at current lambda_bar
theta_bar_2 = (1.0_r8 + gamma * (sqrt(lambda_mean) - 1.0_r8))**2 * sigma_p_2 + sigma_o_2
! Compute constanc coefficient for likelihood at lambda_bar
u_bar = 1.0_r8 / (sqrt(2.0_r8 * PI) * sqrt(theta_bar_2))
! Compute exponent of likelihood at lambda_bar
like_exp_bar = dist_2 / (-2.0_r8 * theta_bar_2)
! Compute exponential part of likelihood at lambda_bar
v_bar = exp(like_exp_bar)
! Compute value of likelihood at current lambda_bar value
like_bar = u_bar * v_bar

end subroutine comp_likelihood

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


!========================================================================
! end module adaptive_inflate_mod
!========================================================================

end module adaptive_inflate_mod
