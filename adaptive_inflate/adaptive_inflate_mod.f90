! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section, 
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module adaptive_inflate_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 
!
! Operations and storage required for various adaptive inflation algorithms

use        types_mod, only : r8, PI
use time_manager_mod, only : time_type, get_time
use    utilities_mod, only : file_exist, get_unit, check_namelist_read, find_namelist_in_file, &
                             register_module, error_handler, E_ERR, E_MSG, logfileunit
use   random_seq_mod, only : random_seq_type, random_gaussian, &
                             init_random_seq, random_uniform


implicit none
private

public :: update_inflation, adaptive_inflate_ss_init, ss_inflate, ss_inflate_sd, &
   ss_sd_lower_bound, adaptive_inflate_end, do_obs_inflate, do_varying_ss_inflate, &
   do_single_ss_inflate, adaptive_inflate_obs_init, deterministic_inflate, &
   obs_inflate, obs_inflate_sd, obs_sd_lower_bound, obs_inf_upper_bound, &
   ss_inf_upper_bound, inflate_ens, output_inflate_diagnostics

character(len = 129) :: errstring

! CVS Generated file description for error handling, do not edit
character(len=128), parameter :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

! Manages both observation space and state space inflation
! Handles initial values and restarts, diagnostic output, and computations
! Algorithm options at present include a single fixed observation space
! adaptive inflation and a spatially-varying state space inflation that carries
! a mean and variance for the state space inflation at each point. 

! Storage for fixed observation space inflation; publically accessible; one per domain/region
real(r8), allocatable :: obs_inflate(:), obs_inflate_sd(:)

! Storage for spatially varying state space inflation; publically accessible
real(r8), allocatable :: ss_inflate(:), ss_inflate_sd(:)

! Flag indicating whether module has been initialized
logical :: initialized = .false.

! Controls for a random sequence for non-deterministic inflation
type (random_seq_type) :: inc_ran_seq
logical :: first_inc_ran_call = .true.

! Controls opening file for inflate diagnostic output
logical :: first_diag_call = .true.
integer :: diag_unit

!============================================================================

!---- namelist with default values

logical              :: do_obs_inflate             = .false.
logical              :: do_varying_ss_inflate      = .false.
logical              :: do_single_ss_inflate       = .false.
logical              :: start_from_inflate_restart = .false.
logical              :: output_restart             = .true.
logical              :: deterministic_inflate      = .true.
character(len = 129) :: inflate_in_file_name       = "inflate_ics"
character(len = 129) :: inflate_out_file_name      = "inflate_restart"
real(r8)             :: obs_inf_initial            = 1.0_r8
real(r8)             :: obs_inf_sd_initial         = 0.0_r8
real(r8)             :: obs_sd_lower_bound         = 0.0_r8
real(r8)             :: obs_inf_upper_bound        = 1000000_r8
real(r8)             :: ss_inf_initial             = 1.0_r8
real(r8)             :: ss_inf_sd_initial          = 0.0_r8
real(r8)             :: ss_sd_lower_bound          = 0.0_r8
real(r8)             :: ss_inf_upper_bound         = 1000000_r8
character(len = 129) :: diagnostic_file_name       = "inflate_diag"

namelist / adaptive_inflate_nml / do_obs_inflate, do_varying_ss_inflate, &
   do_single_ss_inflate, start_from_inflate_restart, inflate_in_file_name, &
   output_restart, deterministic_inflate, inflate_out_file_name, obs_inf_initial, &
   obs_inf_sd_initial, obs_sd_lower_bound, obs_inf_upper_bound, ss_inf_initial, &
   ss_inf_upper_bound, ss_inf_sd_initial, ss_sd_lower_bound, diagnostic_file_name

!============================================================================

contains


!------------------------------------------------------------------

subroutine adaptive_inflate_init()

integer :: iunit, io, true_count

if(initialized) return

call register_module(source, revision, revdate)                                         

! Read the namelist entry
call find_namelist_in_file("input.nml", "adaptive_inflate_nml", iunit)
read(iunit, nml = adaptive_inflate_nml, iostat = io)                                         
call check_namelist_read(iunit, io, "adaptive_inflate_nml")                                  
                                                                                        
! Write the namelist values to the log file                                             
call error_handler(E_MSG,'adaptive_inflate','adaptive_inflate namelist values',' ',' ',' ')  
write(logfileunit, nml=adaptive_inflate_nml)                                                 
write(     *     , nml=adaptive_inflate_nml)    

! Can't do more than one type of inflation
true_count = 0
if(do_obs_inflate) true_count = true_count + 1
if(do_varying_ss_inflate) true_count = true_count + 1
if(do_single_ss_inflate) true_count = true_count + 1

if(true_count > 1) then
   write(errstring, *) 'Can only use one type of inflation at a time'
   call error_handler(E_ERR, 'adaptive_inflate_init', errstring, source, revision, revdate)
endif

end subroutine adaptive_inflate_init

!------------------------------------------------------------------

subroutine adaptive_inflate_ss_init(model_size)

integer, intent(in) :: model_size

integer :: restart_unit, restart_model_size, required_size

! Make sure basic part of module has been initialized already
call adaptive_inflate_init()

! If spatially varying state space, allocate storage
if(do_varying_ss_inflate) then
   allocate(ss_inflate(model_size), ss_inflate_sd(model_size))
else if(do_single_ss_inflate) then
   allocate(ss_inflate(1), ss_inflate_sd(1))
else
   ! either do_obs_inflation or no inflation 
   return
endif

! Read in initial values from file OR get from namelist as requested
! Initially, file is binary only
if(start_from_inflate_restart) then
   ! Determine the required size of the restart file
   if(do_varying_ss_inflate) then
      required_size = model_size
   else if(do_single_ss_inflate) then
      required_size = 1
   endif
   ! Open the file
   restart_unit = get_unit()
   open(unit = restart_unit, file = inflate_in_file_name, action = 'read', form = 'formatted')
   read(restart_unit, *) restart_model_size
   if(restart_model_size /= required_size) then
      write(errstring, *) 'Size of state space restart file is incorrect'
      call error_handler(E_ERR, 'adaptive_inflate_ss_init', &
         errstring, source, revision, revdate)
   endif
   read(restart_unit, *) ss_inflate, ss_inflate_sd, ss_sd_lower_bound 
   close(restart_unit)
else
   ! Get initial values from namelist
   ss_inflate = ss_inf_initial
   ss_inflate_sd = ss_inf_sd_initial
endif

end subroutine adaptive_inflate_ss_init


!------------------------------------------------------------------

subroutine adaptive_inflate_obs_init(num_domains, dont_read_restart)

integer, intent(in) :: num_domains
logical, intent(in) :: dont_read_restart

! Don't read restart when doing parallel runs driven by assim_region

integer :: restart_unit, restart_num_domains

! Make sure basic part of module has been initialized already
call adaptive_inflate_init()

! It is convenient to allocate this storage no matter what
! Avoids a lot of messy logic in assim_tools
allocate(obs_inflate(num_domains), obs_inflate_sd(num_domains))

! Need to initialize the values
! If obs_space inflation is not being done, set everything to -1
if(.not. do_obs_inflate) then
   obs_inflate = -1.0_r8
   obs_inflate_sd = -1.0_r8
else 
   ! Doing obs space inflation, follow rest of namelist params
   if(start_from_inflate_restart .and. .not. dont_read_restart) then
      ! Open the file
      restart_unit = get_unit()
      open(unit = restart_unit, file = inflate_in_file_name, action = 'read', &
         form = 'formatted')
      read(restart_unit, *) restart_num_domains
      ! Number of domains in restart file must match
      if(restart_num_domains /= num_domains) then
         write(errstring, *) 'Number of domains in restart not same as number of domains'
         call error_handler(E_ERR, 'adaptive_inflate_obs_init', &
            errstring, source, revision, revdate)
      endif
      read(restart_unit, *) obs_inflate, obs_inflate_sd, obs_sd_lower_bound
      close(restart_unit)
   else
      ! Get initial values from namelist
      obs_inflate = obs_inf_initial
      obs_inflate_sd = obs_inf_sd_initial
   endif
endif

end subroutine adaptive_inflate_obs_init


!------------------------------------------------------------------


subroutine adaptive_inflate_end

integer :: restart_unit

! if no inflation, no restart file
if(.not.do_obs_inflate .and. .not.do_varying_ss_inflate .and. &
   .not.do_single_ss_inflate) return

! Write restart files if requested
if(output_restart) then
   ! Open the file
   restart_unit = get_unit()
   open(unit = restart_unit, file = inflate_out_file_name, action = 'write', form = 'formatted')

   if(do_obs_inflate) then
      ! Write the size to allow for error check on input
      write(restart_unit, *) size(obs_inflate)
      write(restart_unit, *) obs_inflate, obs_inflate_sd, obs_sd_lower_bound
   else if(do_varying_ss_inflate .or. do_single_ss_inflate) then
      ! Write the size to allow for error check on input
      write(restart_unit, *) size(ss_inflate)
      write(restart_unit, *) ss_inflate, ss_inflate_sd, ss_sd_lower_bound
   endif
endif

if(do_varying_ss_inflate .or. do_single_ss_inflate) then
   deallocate(ss_inflate, ss_inflate_sd)
else if(do_obs_inflate) then
   deallocate(obs_inflate, obs_inflate_sd)
endif

end subroutine adaptive_inflate_end


!------------------------------------------------------------------

subroutine output_inflate_diagnostics(time)

type(time_type), intent(in) :: time

integer :: days, seconds, i

! Diagnostics for spatially-varying state space inflate are done
! By the filter on the state-space netcdf diagnostic files.
! Here, need to do initial naive ascii dump for obs space or
! fixed state space. Assume that values can come from storage in
! this module directly.

! Only need to do something if obs_space or single state space inflate
if(do_obs_inflate .or. do_single_ss_inflate) then
   ! Open the file if this is first time through
   ! Just make it ascii flat for now
   ! Both prior and posterior are going in here for now
   if(first_diag_call) then
      first_diag_call = .false.
      ! Open the file
      diag_unit = get_unit()
      open(unit = diag_unit, file = diagnostic_file_name, action = 'write', form = 'formatted')
   endif

   ! Get the time in days and seconds
   call get_time(time, seconds, days)
   ! Write out the time followed by the values
   if(do_single_ss_inflate) then
      write(diag_unit, *) days, seconds, ss_inflate, ss_inflate_sd
   else if(do_obs_inflate) then
      do i = 1, size(obs_inflate)
         write(diag_unit, *) days, seconds, obs_inflate(i), obs_inflate_sd(i)
      end do
   endif
endif

end subroutine output_inflate_diagnostics


!------------------------------------------------------------------

subroutine inflate_ens(ens, mean, inflate, var_in)

! Inflates subset of ensemble members given mean and inflate
! Should select between deterministic and stochastic inflation

real(r8), intent(inout) :: ens(:)
real(r8), intent(in) :: mean, inflate
real(r8), intent(in), optional :: var_in

integer :: i, ens_size
real(r8) :: rand_sd, var

if(deterministic_inflate) then
   ! Just spread the ensemble out linearly
   ens = (ens - mean) * sqrt(inflate) + mean
else
   ens_size = size(ens)

   ! Use a stochastic algorithm to spread out.
   ! WARNING: This requires that the WHOLE ensemble is available here.
   ! There is no way to confirm that at this level.  For now, (Jan 06)
   ! filter has an error check for the case of ensemble stored out of 
   ! core that flags the use of non-deterministic as illegal.
   ! If var is not present, go ahead and compute it here.
   if(.not. present(var_in)) then
      var = sum((ens - mean)**2) / (ens_size - 1)
   else
      var = var_in
   endif
   
   ! Another option at this point would be to add in random noise designed
   ! To increase the variance of the prior ensemble to the appropriate level
   ! Would probably want to keep the mean fixed by shifting AND do a sort
   ! on the final prior/posterior increment pairs to avoid large regression
   ! error as per stochastic filter algorithms. This might help to avoid
   ! problems with generating gravity waves in the Bgrid model, for instance.
   ! If this is first time through, need to initialize the random sequence
   if(first_inc_ran_call) then
      call init_random_seq(inc_ran_seq)
      first_inc_ran_call = .false.
   endif
   
    ! Figure out required sd for random noise being added
    if(inflate > 1.0_r8) then
       ! Don't allow covariance deflation in this version
       rand_sd = sqrt(inflate*var - var)
       !write(*, *) 'rand_sd increment needed is ', sqrt(prior_var), rand_sd
       ! Add random sample from this noise into the ensemble
       do i = 1, ens_size
          ens(i) = random_gaussian(inc_ran_seq, ens(i), rand_sd)
       end do
       ! Adjust the mean
       ens = ens - (sum(ens) / ens_size - mean)
   endif
endif

end subroutine inflate_ens

!------------------------------------------------------------------

subroutine update_inflation(inflate, inflate_sd, prior_mean, prior_var, &
   obs, obs_var, gamma, sd_lower_bound, inf_upper_bound)

real(r8), intent(inout) :: inflate, inflate_sd
real(r8), intent(in)    :: prior_mean, prior_var, obs, obs_var, gamma
real(r8), intent(in)    :: sd_lower_bound, inf_upper_bound

real(r8) :: new_inflate, new_inflate_sd

! If the inflate_sd is negative, just keep everything the same
if(inflate_sd < 0.0_r8) return
   

! Updates the distribution of inflation given prior mean and sd
! plus the prior mean and variance and obs value and variance.
! gamma is the 'correlation * localization' limiting factor.
! A lower bound on the updated inflation sd and an upper bound
! on the inflation itself are provided. 

call bayes_cov_inflate(prior_mean, prior_var, obs, obs_var, inflate, &
   inflate_sd, gamma, new_inflate, new_inflate_sd, sd_lower_bound)

! Make sure inflate satisfies constraints
inflate = new_inflate
! Should we continue to enforce > 1.0 on inflate?
if(inflate < 1.0_r8) inflate = 1.0_r8
if(inflate > inf_upper_bound) inflate = inf_upper_bound

! Make sure sd satisfies constraints
inflate_sd = new_inflate_sd
if(inflate_sd < sd_lower_bound) inflate_sd = sd_lower_bound

end subroutine update_inflation

!------------------------------------------------------------------


subroutine bayes_cov_inflate(x_p, sigma_p_2, y_o, sigma_o_2, lambda_mean, lambda_sd, &
   gamma, new_cov_inflate, new_cov_inflate_sd, sd_lower_bound_in)

real(r8), intent(in) :: x_p, sigma_p_2, y_o, sigma_o_2, lambda_mean, lambda_sd, gamma
real(r8), intent(in) :: sd_lower_bound_in
real(r8), intent(out) :: new_cov_inflate, new_cov_inflate_sd

integer :: i, mlambda_index(1)

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

if(gamma > 0.99_r8) then
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
   ! If gamma is non-zero, have to approximate with Taylor series for likelihood term
   call linear_bayes(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd_2, gamma, &
      new_cov_inflate, new_cov_inflate_sd, sd_lower_bound_in)
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

real(r8) :: compute_new_density
real(r8), intent(in) :: dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd, gamma, lambda

real(r8) :: theta_2, theta
real(r8) :: exponent_prior, exponent_likelihood


! Compute probability of this lambda being correct
exponent_prior = (lambda - lambda_mean)**2 / (-2.0_r8 * lambda_sd**2)

! Compute probability that observation would have been observed given this lambda
!!! Legacy code for no gamma: theta_2 = lambda * sigma_p_2 + sigma_o_2
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
   new_cov_inflate, new_cov_inflate_sd, sd_lower_bound_in)

real(r8), intent(in) :: dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd_2
real(r8), intent(in) :: gamma, sd_lower_bound_in
real(r8), intent(inout) :: new_cov_inflate, new_cov_inflate_sd

real(r8) :: theta_bar_2, u_bar, like_exp_bar, v_bar, like_bar, like_prime, theta_bar
!!!real(r8) :: like_plus, like_minus, delta_lambda
real(r8) :: a, b, c, disc, plus_root, minus_root, dtheta_dlambda

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

! Next compute derivative of likelihood at this point
!!! Case for no gamma was checked by finite difference, 26 Dec. 2005
!!!like_prime = (u_bar * v_bar * sigma_p_2 / 2.0_r8) * (dist_2 / theta_bar_2**2 - 1.0_r8 / theta_bar_2)

! First compute d/dlambda of theta evaluated at lambda_mean
! Verified correct by finite difference, 1 January, 2006
dtheta_dlambda = 0.5_r8 * sigma_p_2 * gamma *(1.0_r8 - gamma + gamma*sqrt(lambda_mean)) / &
   (theta_bar * sqrt(lambda_mean))
like_prime = (u_bar * v_bar * dtheta_dlambda / theta_bar) * (dist_2 / theta_bar_2 - 1.0_r8)

! FOR NOW, TO AVOID PROBLEMS WITH OVERFLOW/UNDERFLOW, if like_prime is small, assume it is 0
! SHOULD REALLY FIX THIS BETTER
if(abs(like_prime) < 1e-8_r8) then
   new_cov_inflate = lambda_mean
   new_cov_inflate_sd = sqrt(lambda_sd_2)
   return
endif

!!!write(*, *) 'like_bar, like_prime is ', like_bar, like_prime

! Check by computing finite difference like_prime
!delta_lambda = 0.00001_r8
!call comp_likelihood(dist_2, sigma_p_2, sigma_o_2, lambda_mean + delta_lambda, gamma, like_plus)
!call comp_likelihood(dist_2, sigma_p_2, sigma_o_2, lambda_mean - delta_lambda, gamma, like_minus)
!like_prime = (like_plus - like_minus) / (2.0_r8 * delta_lambda)
!write(*, *) 'finite_diff like_prime ', like_prime

! Given like_bar and like_prime can find mode of product of linear likelihood and prior
! Get a quadratic for lambda to find estimated mode
!!!a = -1.0_r8 * like_prime / lambda_sd_2
!!!b = 2.0_r8 * like_prime * lambda_mean / lambda_sd_2  - like_bar / lambda_sd_2
!!!c = like_prime + like_bar * lambda_mean / lambda_sd_2 - like_prime * lambda_mean**2 / lambda_sd_2

a = 1.0_r8
b = like_bar / like_prime - 2.0_r8 * lambda_mean
c = lambda_mean**2 -lambda_sd_2 - like_bar * lambda_mean / like_prime

! Find the roots using quadratic formula
!!!disc = b**2 - 4.0_r8 * a * c
disc = b**2 - 4.0_r8 * c
if(disc < 0.0_r8) then
   write(*, *) 'disc is negative in linear_bayes: An algorithmic failure', disc
   stop
endif   

!!!plus_root = (-1.0*b + sqrt(disc)) / (2.0_r8 * a)
!!!minus_root = (-1.0*b - sqrt(disc)) / (2.0_r8 * a)
plus_root = (-1.0*b + sqrt(disc)) / 2.0_r8
minus_root = (-1.0*b - sqrt(disc)) / 2.0_r8

! Always pick the minus root???
! Do a check to pick closest root
if(abs(minus_root - lambda_mean) < abs(plus_root - lambda_mean)) then
   new_cov_inflate = minus_root
   !write(*, *) 'minus root selected in linear_bayes: This is unexpected but probably ok'
   !stop
else
   new_cov_inflate = plus_root
   !write(*, *) 'plus root selected in linear_bayes: This is unexpected but probably ok'
   !stop
endif


end subroutine linear_bayes


!------------------------------------------------------------------------

subroutine comp_likelihood(dist_2, sigma_p_2, sigma_o_2, lambda_mean, gamma, like_bar)
real(r8), intent(in) :: dist_2, sigma_p_2, sigma_o_2, lambda_mean, gamma
real(r8), intent(out) :: like_bar

real :: theta_bar_2, u_bar, like_exp_bar, v_bar

! Compute value of theta at current lambda_bar
!!! Non-gamma legacy code: theta_bar_2 = lambda_mean * sigma_p_2 + sigma_o_2
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




!========================================================================
! end module adaptive_inflate_mod
!========================================================================

end module adaptive_inflate_mod
