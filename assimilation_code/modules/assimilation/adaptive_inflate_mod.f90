! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Implements adaptive inflation algorithms. Defines a type that keeps trace of 
! algorithm parameter settings. Implements constant state space inflation, single
! value adaptive state space inflation, spatially-varying state space inflation
! via the legacy algorithm of Anderson and the enhanced algorithm of Gharamti. 

module adaptive_inflate_mod

!> \defgroup adaptive_inflate adaptive_inflate_mod
!> @{

use types_mod,            only : r8, PI, MISSING_R8
use utilities_mod,        only : error_handler, E_ERR, E_MSG, &
                                 nmlfileunit, do_nml_file, do_nml_term,              &
                                 check_namelist_read, find_namelist_in_file
use random_seq_mod,       only : random_seq_type, random_gaussian, init_random_seq
use mpi_utilities_mod,    only : my_task_id

implicit none
private

public :: adaptive_inflate_type, adaptive_inflate_init, log_inflation_info,         &
          update_inflation, update_varying_state_space_inflation,                   &
          inflate_ens, solve_quadratic,                                             & 
          set_inflate_flavor, get_inflate_initial_mean, get_inflate_initial_sd,     &
          do_varying_ss_inflate, do_single_ss_inflate, do_obs_inflate,              &
          do_deterministic_inflate, do_ss_inflate, do_rtps_inflate,                 &
          NO_INFLATION

character(len=*), parameter :: source = 'adaptive_inflate_mod.f90'

! Encode the different inflation options
! OBS_INFLATION is currently deprecated.

integer, parameter :: NO_INFLATION               = 0
integer, parameter :: OBS_INFLATION              = 1
integer, parameter :: VARYING_SS_INFLATION       = 2
integer, parameter :: SINGLE_SS_INFLATION        = 3
integer, parameter :: RELAXATION_TO_PRIOR_SPREAD = 4
integer, parameter :: ENHANCED_SS_INFLATION      = 5

! Type that defines the parameter setting for an application of inflation
type adaptive_inflate_type
   private
   integer               :: flavor
   logical               :: deterministic
   real(r8)              :: initial_mean
   real(r8)              :: initial_sd 
   real(r8)              :: mean_lower_bound 
   real(r8)              :: mean_upper_bound
   real(r8)              :: sd_lower_bound 
   real(r8)              :: sd_max_change
   real(r8)              :: damping
   real(r8)              :: rtps_relaxation
   ! Include a random sequence type in case non-deterministic inflation is used
   type(random_seq_type) :: ran_seq
end type adaptive_inflate_type

! Module storage for writing error messages
character(len=512) :: string1, string2

! Flag indicating whether module has been initialized
logical :: initialized = .false.

! Used for precision tests in inflation update routines
real(r8), parameter    :: small = epsilon(1.0_r8)   ! threshold for avoiding NaNs/Inf


!----------------------------------------------------------------
! Namelist input with default values
!
integer  :: flavor                       = NO_INFLATION
logical  :: deterministic                = .true.
real(r8) :: initial_mean                 = 1.0_r8
real(r8) :: initial_sd                   = 0.0_r8
real(r8) :: mean_lower_bound             = 1.0_r8
real(r8) :: mean_upper_bound             = 1000000.0_r8
real(r8) :: sd_lower_bound               = 0.0_r8
real(r8) :: sd_max_change                = 1.05_r8
real(r8) :: damping                      = 1.0_r8
real(r8) :: rtps_relaxation              = 1.0_r8

namelist /adaptive_inflate_nml/ flavor, &
   sd_max_change,                       &
   deterministic,                       &
   damping,                             &
   initial_mean,                        &
   initial_sd,                          &
   mean_lower_bound,                    &
   mean_upper_bound,                    &
   sd_lower_bound,                      &
   rtps_relaxation


!===============================================================================

contains

!-------------------------------------------------------------------------------
!>

subroutine set_inflate_flavor(inflation_handle, flavor)

type(adaptive_inflate_type) :: inflation_handle
integer, intent(in)         :: flavor

! No error checking done on value
inflation_handle%flavor = flavor

end subroutine set_inflate_flavor

!-------------------------------------------------------------------------------
!>

function get_inflate_initial_mean(inflation)

type(adaptive_inflate_type) :: inflation
real(r8)  :: get_inflate_initial_mean

get_inflate_initial_mean = inflation%initial_mean

end function


!-------------------------------------------------------------------------------
!>

function get_inflate_initial_sd(inflation)

type(adaptive_inflate_type) :: inflation
real(r8)  :: get_inflate_initial_sd

get_inflate_initial_sd = inflation%initial_sd

end function


!-------------------------------------------------------------------------------
!>

! Returns true if any of the adaptive state space inflations are in use
function do_ss_inflate(inflation)

type(adaptive_inflate_type), intent(in) :: inflation
logical :: do_ss_inflate

if (do_single_ss_inflate(inflation) .or. &
    do_orig_varying_ss_inflate(inflation) .or. &
    do_enhanced_varying_ss_inflate(inflation) then
   do_ss_inflate = .true.
else
   do_ss_inflate = .false.
endif

end function do_ss_inflate


!-------------------------------------------------------------------------------
!>

function do_varying_ss_inflate(inflation)

! Returns true for any of the spatially varying inflations
! ALL OF THESE NEED TO BE REMOVED FROM THIS MODULE

type(adaptive_inflate_type), intent(in) :: inflation
logical :: do_varying_ss_inflate

if (do_enhanced_varying_ss_inflate(inflation) .or. &
    do_orig_varying_ss_inflate(inflation)) then
   do_varying_ss_inflate = .true.
else
   do_varying_ss_inflate = .false.
endif

end function do_varying_ss_inflate


!-------------------------------------------------------------------------------
!> Initializes an adaptive_inflate_type 

subroutine adaptive_inflate_init(inflate_handle)

type(adaptive_inflate_type), intent(inout) :: inflate_handle

! random value
integer, save :: salt = 139
integer       :: iunit, io

! Record the module version if this is first initialize call
if(.not. initialized) then

   ! Read the namelist settings for prior and posterior inflation
   call find_namelist_in_file("input.nml", "adaptive_inflate_nml", iunit)
   read(iunit, nml = adaptive_inflate_nml, iostat = io)
   call check_namelist_read(iunit, io, "adpative_inflate_nml")
   
   ! Record the namelist values used for the run ...
   if (do_nml_file()) write(nmlfileunit, nml=adaptive_inflate_nml)
   if (do_nml_term()) write(     *     , nml=adaptive_inflate_nml)


   initialized = .true.
endif

! Load up the structure first to keep track of all details of this inflation type
inflate_handle%flavor               = flavor
inflate_handle%deterministic        = deterministic
inflate_handle%initial_mean         = initial_mean
inflate_handle%initial_sd           = initial_sd
inflate_handle%mean_lower_bound     = mean_lower_bound
inflate_handle%mean_upper_bound     = mean_upper_bound
inflate_handle%sd_lower_bound       = sd_lower_bound
inflate_handle%sd_max_change        = sd_max_change
inflate_handle%damping              = damping
inflate_handle%rtps_relaxation      = rtps_relaxation

! Make sure selected options are okay
call validate_inflate_options(inflate_handle)

! If non-deterministic inflation is being done, need to initialize random sequence.
! use the task id number (plus 1 since they start at 0) to set the initial seed.
! NOTE: non-deterministic inflation does NOT reproduce as process count is varied!
! The use of the saved variable salt means that multiple calls to the initialization
! will use different seeds so there is no danger of unexpected correlations.
if(.not. inflate_handle%deterministic) then
   call init_random_seq(inflate_handle%ran_seq, my_task_id()+1 + salt)
   salt = salt + 1000 
endif

! Not currently supporting observation space inflation
if(inflate_handle%flavor == OBS_INFLATION) then

   write(string1,  *) 'No longer supporting observation space inflation ', &
                        '(i.e. flavor = 1).'
   write(string2, *) 'Please contact dart@ucar.edu if you would like to use ', &
                        'observation space inflation'
   call error_handler(E_ERR, 'adaptive_inflate_init', string1, source, text2=string2)

endif

end subroutine adaptive_inflate_init


!-------------------------------------------------------------------------------
!> Returns true if this inflation type indicates observation space inflation

function do_obs_inflate(inflate_handle)

logical                                 :: do_obs_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

do_obs_inflate = (inflate_handle%flavor == OBS_INFLATION)

if (do_obs_inflate) then
  write(string1,  *) 'observation space inflation not suppported (i.e. flavor = 1)'
  write(string2, *) 'please contact dart if you would like to use this functionality'
  call error_handler(E_ERR, 'do_obs_inflate', string1, source, text2=string2)
endif

end function do_obs_inflate


!-------------------------------------------------------------------------------
!> Returns true if this inflation type indicates varying state space inflation

function do_orig_varying_ss_inflate(inflate_handle)

logical                                 :: do_orig_varying_ss_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

do_orig_varying_ss_inflate = (inflate_handle%flavor == VARYING_SS_INFLATION)

end function do_orig_varying_ss_inflate


!-------------------------------------------------------------------------------
!> Returns true if this inflation type indicates fixed state space inflation

function do_single_ss_inflate(inflate_handle)

logical                                 :: do_single_ss_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

do_single_ss_inflate = (inflate_handle%flavor == SINGLE_SS_INFLATION)

end function do_single_ss_inflate


!-------------------------------------------------------------------------------
!> Returns true if this inflation type indicates posterior relaxion-to-prior-spread
!> (whitaker & Hamill, 2012)

function do_rtps_inflate(inflate_handle)

logical                                 :: do_rtps_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

do_rtps_inflate = (inflate_handle%flavor == RELAXATION_TO_PRIOR_SPREAD)

end function do_rtps_inflate


!-------------------------------------------------------------------------------
!>
!> Returns true if this inflation sub type indicates enhanced state space inflation
!> Moha Gharamti, 2017

function do_enhanced_varying_ss_inflate(inflate_handle)

logical                                 :: do_enhanced_varying_ss_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

do_enhanced_varying_ss_inflate = (inflate_handle%flavor == ENHANCED_SS_INFLATION)

end function do_enhanced_varying_ss_inflate


!-------------------------------------------------------------------------------
!> Returns true if deterministic inflation is indicated

function do_deterministic_inflate(inflate_handle)

logical :: do_deterministic_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

do_deterministic_inflate = inflate_handle%deterministic

end function do_deterministic_inflate


!-------------------------------------------------------------------------------
!> Make sure the combination of inflation options are legal

subroutine validate_inflate_options(inflation_handle)

type(adaptive_inflate_type), intent(inout) :: inflation_handle

character(len=32) :: string

if(inflation_handle%flavor < NO_INFLATION .or. inflation_handle%flavor > ENHANCED_SS_INFLATION) then
   write(string1, *) 'flavor=', inflation_handle%flavor, ' Must be 0, 1, 2, 3, 4, or 5 '
   call error_handler(E_ERR,'validate_inflate_options', string1, source, &
                          text2='Inflation type for '//string)
endif

if(inflation_handle%damping < 0.0_r8 .or. inflation_handle%damping > 1.0_r8) then
   write(string1, *) 'damping=', inflation_handle%damping, ' Must be 0.0 <= d <= 1.0'
   call error_handler(E_ERR,'validate_inflate_options', string1, source, &
                             text2='Inflation damping for '//string)
endif

! Observation space inflation not currently supported
if(inflation_handle%flavor == OBS_INFLATION)             &
   call error_handler(E_ERR, 'validate_inflate_options', &
   'observation space inflation (type 1) not currently supported', source, &
   text2='contact DART developers if you are interested in using it.')

! enhanced inflation checks - this is before we set the subflavor in the structure.
if (inflation_handle%flavor     == ENHANCED_SS_INFLATION) then

   ! check sd_max_change() for valid range
   if (inflation_handle%sd_max_change < 1.0_r8 .or. inflation_handle%sd_max_change > 2.0_r8) then
      write(string1, *) 'sd_max_change=', inflation_handle%sd_max_change, &
                        ' Must be 1.0 <= X <= 2.0'
      call error_handler(E_ERR,'validate_inflate_options', string1, source, &
                                text2='Inflation stddev max change for '//string)
   endif
endif

! Cannot support non-determistic inflation and a mean_lower_bound < 1
if(.not. inflation_handle%deterministic .and. inflation_handle%mean_lower_bound < 1.0_r8) then
   write(string1, *) 'Cannot have non-deterministic inflation and mean_lower_bound < 1'
   call error_handler(E_ERR, 'validate_inflate_options', string1, source)
endif


end subroutine validate_inflate_options


!-------------------------------------------------------------------------------
!> Inflates subset of ensemble members given mean and inflate
!> Selects between deterministic and stochastic inflation

subroutine inflate_ens(inflate_handle, ens, mean, inflate, var_in, fsprd, asprd)

type(adaptive_inflate_type), intent(inout) :: inflate_handle
real(r8),                    intent(inout) :: ens(:)
real(r8),                    intent(in)    :: mean 
real(r8),                    intent(inout) :: inflate
real(r8), optional,          intent(in)    :: var_in
real(r8), optional,          intent(in)    :: fsprd, asprd

integer  :: i, ens_size
real(r8) :: rand_sd, var, sd_inflate

! This used to be restricted to case where clm supports missing
! For now, if there is a missing_r8 in any ensemble, don't do anything and return
if (any(ens == MISSING_R8)) return

! Damp the inflation if requested; 
if(inflate_handle%damping /= 1.0_r8) then
   inflate = 1.0_r8 + inflate_handle%damping * (inflate - 1.0_r8)     
endif 

if ( do_rtps_inflate(inflate_handle)) then
   if ( .not. present(fsprd) .or. .not. present(asprd)) then 
      write(string1, *) 'missing arguments for RTPS inflation'
      call error_handler(E_ERR,'inflate_ens',string1,source) 
   endif 
   ! only inflate if spreads are > 0
   if ( asprd .gt. 0.0_r8 .and. fsprd .gt. 0.0_r8) then
       ! Saving inflate for possible diagnostics
       inflate = 1.0_r8 + inflate_handle%rtps_relaxation * ((fsprd-asprd) / asprd) 
   else
       inflate = 1.0_r8
   endif
   ens = ens * inflate + mean * (1.0_r8 - inflate)
else 
   if(inflate_handle%deterministic) then

      ! Spread the ensemble out linearly for deterministic
      ! Following line can lead to inflation of 1.0 changing ens on some compilers
      !!! ens = (ens - mean) * sqrt(inflate) + mean
      ! Following gives 1.0 inflation having no impact on known compilers
      sd_inflate = sqrt(inflate) 
      ens = ens * sd_inflate + mean * (1.0_r8 - sd_inflate)

   else
      ! Use a stochastic algorithm to spread out.

      ! If var is not present, go ahead and compute it here.
      if(.not. present(var_in)) then
         ens_size = size(ens)
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
endif
end subroutine inflate_ens


!-------------------------------------------------------------------------------
!> This routine is given information from an inflate type, scalar values for inflate 
!> and inflate_sd, the ensemble prior_mean and prior_var for an observation, the
!> ensemble (or group) size, the observed value, observational error variance, 
!> and gamma_corr (see below for description).  It computes updated values for the 
!> inflate and inflate_sd items using algorithms described in filter.html (including
!> references to papers).
!>
!> The gamma_corr parameter gives the localized prior correlation times the localization
!> which is computed in the assim_tools routine filter_assim. For single state
!> space inflation it is 1.0.

subroutine update_inflation(inflate_handle, inflate, inflate_sd, prior_mean, prior_var, &
   ens_size, obs, obs_var, gamma_corr)

type(adaptive_inflate_type), intent(in)    :: inflate_handle
real(r8),                    intent(inout) :: inflate, inflate_sd
real(r8),                    intent(in)    :: prior_mean, prior_var
integer,                     intent(in)    :: ens_size
real(r8),                    intent(in)    :: obs, obs_var, gamma_corr

real(r8) :: new_inflate, new_inflate_sd
integer :: inf_type

! If the inflate_sd not positive, keep everything the same
if(inflate_sd <= 0.0_r8) return

! A lower bound on the updated inflation sd and an upper bound
! on the inflation itself are provided in the inflate_handle. 

! Use bayes theorem to update
call bayes_cov_inflate(ens_size, inf_type, prior_mean, prior_var, obs, obs_var, inflate, &
   inflate_sd, gamma_corr, inflate_handle, &
   new_inflate, new_inflate_sd)

! Make sure inflate satisfies constraints
inflate = new_inflate
if(inflate < inflate_handle%mean_lower_bound) inflate = inflate_handle%mean_lower_bound
if(inflate > inflate_handle%mean_upper_bound) inflate = inflate_handle%mean_upper_bound

! Make sure sd satisfies constraints
inflate_sd = new_inflate_sd
if(inflate_sd < inflate_handle%sd_lower_bound) inflate_sd = inflate_handle%sd_lower_bound

end subroutine update_inflation

!-------------------------------------------------------------------------------
!> Computes updated inflation mean and inflation sd for varying state space inflation

subroutine update_varying_state_space_inflation(inflate, inflate_mean, inflate_sd, &
   ss_inflate_base, orig_obs_prior_mean, orig_obs_prior_var, obs, obs_err_var, &
   ens_size, reg_factor, correl, inflate_only)

type(adaptive_inflate_type), intent(in)    :: inflate
real(r8),                    intent(inout) :: inflate_mean
real(r8),                    intent(inout) :: inflate_sd
real(r8),                    intent(in)    :: ss_inflate_base
real(r8),                    intent(in)    :: orig_obs_prior_mean
real(r8),                    intent(in)    :: orig_obs_prior_var
real(r8),                    intent(in)    :: obs
real(r8),                    intent(in)    :: obs_err_var
integer,                     intent(in)    :: ens_size
real(r8),                    intent(in)    :: reg_factor
real(r8),                    intent(in)    :: correl
logical,                     intent(in)    :: inflate_only

real(r8) :: gamma, ens_var_deflate, r_var, r_mean

if(inflate_mean <= 0.0_r8 .or. inflate_sd <= 0.0_r8) return

! Gamma is less than 1 for varying ss, see adaptive inflate module
gamma = reg_factor * abs(correl)

! Remove the impact of inflation to allow efficient single pass with assim.
if ( abs(gamma) > small ) then
   ens_var_deflate = orig_obs_prior_var / &
      (1.0_r8 + gamma*(sqrt(ss_inflate_base) - 1.0_r8))**2
else
   ens_var_deflate = orig_obs_prior_var
endif

! If this is inflate only (i.e. posterior) remove impact of this obs.
if(inflate_only .and. &
      ens_var_deflate               > small .and. &
      obs_err_var                   > small .and. &
      obs_err_var - ens_var_deflate > small ) then
   r_var  = 1.0_r8 / (1.0_r8 / ens_var_deflate - 1.0_r8 / obs_err_var)
   r_mean = r_var *(orig_obs_prior_mean / ens_var_deflate - obs / obs_err_var)
else
   r_var = ens_var_deflate
   r_mean = orig_obs_prior_mean
endif

! IS A TABLE LOOKUP POSSIBLE TO ACCELERATE THIS?
! Update the inflation values
call update_inflation(inflate, inflate_mean, inflate_sd, &
   r_mean, r_var, ens_size, obs, obs_err_var, gamma)

end subroutine update_varying_state_space_inflation


!-------------------------------------------------------------------------------
!> Uses one of 2 algorithms in references on DART web site to update the 
!> distribution of inflation:  Anderson 2007, 2009 or Gharamti 2017

subroutine bayes_cov_inflate(ens_size, inf_type, x_p, sigma_p_2, y_o, sigma_o_2, &
                 lambda_mean, lambda_sd, gamma_corr, inflate_handle, &
                 new_cov_inflate, new_cov_inflate_sd)

integer , intent(in)  :: ens_size, inf_type
real(r8), intent(in)  :: x_p, sigma_p_2, y_o, sigma_o_2, lambda_mean, lambda_sd
real(r8), intent(in)  :: gamma_corr
type(adaptive_inflate_type), intent(in) :: inflate_handle
real(r8), intent(out) :: new_cov_inflate, new_cov_inflate_sd

real(r8) :: dist_2, rate, shape_old, shape_new, rate_new
real(r8) :: lambda_sd_2, density_1, density_2, omega, ratio
real(r8) :: new_1_sd, new_max

! If gamma is 0, nothing changes
if(gamma_corr <= 0.0_r8) then
   new_cov_inflate = lambda_mean
   new_cov_inflate_sd = lambda_sd
   return
endif

! Inflation variance
lambda_sd_2 = lambda_sd**2

! Squared Innovation
dist_2 = (y_o - x_p)**2
   
! this block of code no longer being used.  it's here for historical purposes.

!integer  :: i, mlambda_index(1)
!real(r8) :: b, c, d, Q, R, disc, alpha, beta, cube_root_alpha, cube_root_beta, x
!real(r8) :: rrr, cube_root_rrr, angle, mx(3), sep(3), mlambda(3)

!   ! Use ONLY the linear approximation, cubic solution below can be numerically
!   ! unstable for extreme cases. Should look at this later.
!   if(gamma_corr > 0.99_r8) then
!   
!   ! The solution of the cubic below only works if gamma is 1.0
!   ! Can analytically find the maximum of the product: d/dlambda is a
!   ! cubic polynomial in lambda**2; solve using cubic formula for real root
!   ! Can write so that coefficient of x**3 is 1, other coefficients are:
!      b = -1.0_r8 * (sigma_o_2 + sigma_p_2 * lambda_mean)
!      c = lambda_sd_2 * sigma_p_2**2 / 2.0_r8
!      d = -1.0_r8 * (lambda_sd_2 * sigma_p_2**2 * dist_2) / 2.0_r8
!   
!      Q = c - b**2 / 3
!      R = d + (2 * b**3) / 27 - (b * c) / 3
!   
!      ! Compute discriminant, if this is negative have 3 real roots, else 1 real root
!      disc = R**2 / 4 + Q**3 / 27
!   
!      if(disc < 0.0_r8) then
!         rrr = sqrt(-1.0 * Q**3 / 27)
!         ! Note that rrr is positive so no problem for cube root
!         cube_root_rrr = rrr ** (1.0 / 3.0)
!         angle = acos(-0.5 * R / rrr)
!         do i = 0, 2
!            mx(i+1) = 2.0_r8 * cube_root_rrr * cos((angle + i * 2.0_r8 * PI) / 3.0_r8) - b / 3.0_r8
!            mlambda(i + 1) = (mx(i + 1) - sigma_o_2) / sigma_p_2
!            sep(i+1) = abs(mlambda(i + 1) - lambda_mean)
!         end do
!         ! Root closest to initial peak is appropriate
!         mlambda_index = minloc(sep)
!         new_cov_inflate = mlambda(mlambda_index(1))
!   
!      else
!         ! Only one real root here, find it.
!   
!         ! Compute the two primary terms
!         alpha = -R/2 + sqrt(disc)
!         beta = R/2 + sqrt(disc)
!   
!         cube_root_alpha = abs(alpha) ** (1.0 / 3.0) * abs(alpha) / alpha
!         cube_root_beta = abs(beta) ** (1.0 / 3.0) * abs(beta) / beta
!   
!         x = cube_root_alpha - cube_root_beta - b / 3.0
!   
!         ! This root is the value of x = theta**2
!         new_cov_inflate = (x - sigma_o_2) / sigma_p_2
!   
!      endif
!   
!      ! Put in code to approximate the mode (new_cov_inflate)
!      !write(*, *) 'old, orig mode is ', lambda_mean, new_cov_inflate
!   endif

! NOTE THAT THE SINGLE CASE CAN NOW COME THROUGH HERE
if (do_orig_varying_ss_inflate(inflate_handle) .or. &
    do_single_ss_inflate(inflate_handle)) then

   ! Approximate with Taylor series for likelihood term
   call linear_bayes(dist_2, sigma_p_2, sigma_o_2, lambda_mean, lambda_sd_2, gamma_corr, &
      new_cov_inflate)

   ! Bail out to save cost when lower bound is reached on lambda standard deviation
   ! The original test to see if lambda_sd was less than the lower bound
   ! would sometimes return false because of roundoff error and the computation
   ! would go through the expensive part of the code when the minimum was
   ! really reached.  (sd_lower_bound comes from a namelist and precision
   ! errors may be because of the conversion between ascii, single precision
   ! and double precision.)  In any case, the test was changed to return if
   ! the value is within TINY of the limit.
   if(abs(lambda_sd - inflate_handle%sd_lower_bound) <= TINY(0.0_r8)) then
      new_cov_inflate_sd = lambda_sd
      return
   else
      ! Compute by forcing a Gaussian fit at one positive SD
      ! First compute the new_max value for normalization purposes
      new_max = compute_new_density(dist_2, sigma_p_2, sigma_o_2, &
                  lambda_mean, lambda_sd, gamma_corr, new_cov_inflate)
   
      ! Find value at a point one OLD sd above new mean value
      new_1_sd = compute_new_density(dist_2, sigma_p_2, sigma_o_2, &
                 lambda_mean, lambda_sd, gamma_corr, new_cov_inflate + lambda_sd)
   
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
      ! sigma = sqrt(-x^2 / (2 ln(r))
      ! where r is ratio and x is lambda_sd (distance from mean)
      new_cov_inflate_sd = sqrt( -1.0_r8 * lambda_sd_2 / (2.0_r8 * log(ratio)))
   
      ! Prevent an increase in the sd of lambda???
      ! For now, this is mostly countering numerical errors in this computation
      if(new_cov_inflate_sd > lambda_sd) then
         new_cov_inflate_sd = lambda_sd
         return
      endif
   
   endif

else if (do_enhanced_varying_ss_inflate(inflate_handle)) then

   ! Transform Gaussian prior to Inverse Gamma
   call change_GA_IG(lambda_mean, lambda_sd_2, rate)

   ! Approximate with Taylor series for likelihood term
   call enh_linear_bayes(dist_2, sigma_p_2, sigma_o_2,lambda_mean, &
                    gamma_corr, ens_size, rate, new_cov_inflate)

   ! Bail out to save cost when lower bound is reached on lambda standard deviation
   ! See comment in Anderson case for why we use abs and TINY for this comparison.
   if(abs(lambda_sd - inflate_handle%sd_lower_bound) <= TINY(0.0_r8)) then
      new_cov_inflate_sd = lambda_sd
      return 
   else
      ! Compute the shape parameter of the prior IG
      ! This comes from the assumption that the mode of the IG is the mean/mode 
      ! of the input Gaussian
      shape_old = rate / lambda_mean - 1.0_r8
      if (shape_old <= 2.0_r8) then
         new_cov_inflate_sd = lambda_sd
         return
      endif
   
      ! Evaluate exact IG posterior at p1: \lambda_u+\sigma_{\lambda_b} & p2: \lambda_u
      density_1 = enh_compute_new_density(dist_2, ens_size, sigma_p_2, sigma_o_2, &
                         shape_old, rate, gamma_corr, new_cov_inflate+lambda_sd)
      density_2 = enh_compute_new_density(dist_2, ens_size, sigma_p_2, sigma_o_2, &
                         shape_old, rate, gamma_corr, new_cov_inflate)
   
      ! Computational errors check (small numbers + NaNs)
      if (abs(density_1) <= TINY(0.0_r8) .OR. &
          abs(density_2) <= TINY(0.0_r8) .OR. &
          density_1 /= density_1 .OR. density_1 /= density_1 .OR. &
          density_2 /= density_2 .OR. density_2 /= density_2) then
         new_cov_inflate_sd = lambda_sd
         return
      endif
   
      ! Now, compute omega and the new distribution parameters
      ratio     = density_1 / density_2
      omega     = log(new_cov_inflate          )/new_cov_inflate + 1.0_r8/new_cov_inflate - &
                  log(new_cov_inflate+lambda_sd)/new_cov_inflate - 1.0_r8/(new_cov_inflate+lambda_sd)
      rate_new  = log(ratio) / omega
      shape_new = rate_new / new_cov_inflate - 1.0_r8
   
      ! Finally, get the sd of the IG posterior
      if (shape_new <= 2.0_r8) then
         new_cov_inflate_sd = lambda_sd
         return
      endif
      new_cov_inflate_sd = sqrt(rate_new**2 / ( (shape_new-1.0_r8)**2 * (shape_new-2.0_r8) ))
   
      ! If the updated variance is more than xx% the prior variance, keep the prior unchanged 
      ! for stability reasons. Also, if the updated variance is NaN (not sure why this
      ! can happen; never did when developing this code), keep the prior variance unchanged. 
      if ( new_cov_inflate_sd > inflate_handle%sd_max_change*lambda_sd .OR. &
           new_cov_inflate_sd /= new_cov_inflate_sd) then
         new_cov_inflate_sd = lambda_sd
         return
      endif
   endif
   
else
   write(string1, *) 'Internal error, should not happen.  Illegal value for bayes type.'
   call error_handler(E_ERR, 'bayes_cov_inflate', string1, source)
endif

end subroutine bayes_cov_inflate


!-------------------------------------------------------------------------------
!> Used to update density by taking approximate gaussian product
!> original routine.

function compute_new_density(dist_2, sigma_p_2, sigma_o_2, &
                             lambda_mean, lambda_sd, gamma, lambda)

real(r8), intent(in) :: dist_2, sigma_p_2, sigma_o_2
real(r8), intent(in) :: lambda_mean, lambda_sd, gamma, lambda
real(r8)             :: compute_new_density

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


!-------------------------------------------------------------------------------
!>

function enh_compute_new_density(dist_2, ens_size, sigma_p_2, sigma_o_2, &
                                 alpha, beta, gamma_corr, lambda)

! Used to update density by taking approximate gaussian product
real(r8), intent(in) :: dist_2
integer , intent(in) :: ens_size
real(r8), intent(in) :: sigma_p_2, sigma_o_2, gamma_corr, lambda
real(r8), intent(in) :: alpha, beta
real(r8)             :: enh_compute_new_density

real(r8) :: theta, fac1, fac2
real(r8) :: exp_prior, exp_like

! Compute probability of this lambda being correct
exp_prior = - beta / lambda

! Compute probability that observation would have been observed given this lambda
fac1 = (1.0_r8 + gamma_corr * (sqrt(lambda) - 1.0_r8))**2
fac2 = -1.0_r8 / ens_size
if ( fac1 < abs(fac2) ) fac2 = 0.0_r8

theta    = sqrt( (fac1+fac2) * sigma_p_2 + sigma_o_2 )
exp_like = - 0.5_r8 * dist_2 / theta**2

! Modern versions of compilers have intrinsic gamma functions.
! If you have an unresolved external for the gamma function, you should
! either try a newer compiler or code your own gamma function here.
! We know pg compiler versions before pgf90 15.1 do not contain the
! gamma function. If you are never going to use ENHANCED_SS_INFLATION
! you could also just comment out the computation for enh_compute_new_density
! and uncomment the following code block. This will ensure that if you ever
! did try to use ENHANCED_SS_INFLATION, it would appropriately fail.

! write(string1,*)'gamma function not available'
! write(string2,*)'when available uncomment block below and recompile'
! call error_handler(E_ERR, 'enh_compute_new_density', string1, source, text2=string2)

! Compute the updated probability density for lambda
enh_compute_new_density = beta**alpha / gamma(alpha)  * &
                         lambda**(- alpha - 1.0_r8)  / &
                         (sqrt(2.0_r8 * PI) * theta) * &
                         exp(exp_like + exp_prior)

end function enh_compute_new_density


!-------------------------------------------------------------------------------
!> original linear_bayes routine

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


!-------------------------------------------------------------------------------
!> enhanced linear bayes

subroutine enh_linear_bayes(dist_2, sigma_p_2, sigma_o_2, &
           lambda_mean, gamma_corr, ens_size, beta, new_cov_inflate)

real(r8), intent(in)    :: dist_2, sigma_p_2, sigma_o_2, lambda_mean
real(r8), intent(in)    :: gamma_corr, beta
real(r8), intent(inout) :: new_cov_inflate

integer  :: ens_size
real(r8) :: theta_bar_2, like_bar, like_prime, theta_bar
real(r8) :: a, b, c, plus_root, minus_root, deriv_theta
real(r8) :: fac1, fac2, like_ratio

! Scaling factors
fac1 = (1.0_r8 + gamma_corr * (sqrt(lambda_mean) - 1.0_r8))**2
fac2 = -1.0_r8 / ens_size

! Compute value of theta at current lambda_mean
if ( fac1 < abs(fac2) ) fac2 = 0.0_r8
theta_bar_2 = (fac1+fac2) * sigma_p_2 + sigma_o_2
theta_bar   = sqrt(theta_bar_2)

! Compute constant coefficient for likelihood at lambda_bar
like_bar = exp(- 0.5_r8 * dist_2 / theta_bar_2) / (sqrt(2.0_r8 * PI) * theta_bar)

! If like_bar goes to 0, can't do anything, so just keep current values
! Density at current inflation value must be positive
if(like_bar <= 0.0_r8) then
   new_cov_inflate = lambda_mean
   return
endif

! Next compute derivative of likelihood at this point
deriv_theta = 0.5_r8 * sigma_p_2 * gamma_corr * ( 1.0_r8 - gamma_corr + &
              gamma_corr * sqrt(lambda_mean) ) / ( theta_bar * sqrt(lambda_mean) )
like_prime  = like_bar * deriv_theta * (dist_2 / theta_bar_2 - 1.0_r8) / theta_bar

! If like_prime goes to 0, can't do anything, so just keep current values
! We're dividing by the derivative in the quadratic equation, so this
! term better non-zero!
if(     like_prime == 0.0_r8 .OR. &
   abs(like_bar)   <= TINY(0.0_r8) .OR. &
   abs(like_prime) <= TINY(0.0_r8) ) then
   new_cov_inflate = lambda_mean
   return
endif
like_ratio = like_bar / like_prime

a = 1.0_r8 - lambda_mean / beta
b = like_ratio - 2.0_r8 * lambda_mean
c = lambda_mean**2 - like_ratio * lambda_mean

! Use nice scaled quadratic solver to avoid precision issues
call solve_quadratic(a, b, c, plus_root, minus_root)

! Do a check to pick closest root
if(abs(minus_root - lambda_mean) < abs(plus_root - lambda_mean)) then
   new_cov_inflate = minus_root
else
   new_cov_inflate = plus_root
endif

! Do a final check on the sign of the updated factor
! Sometimes the factor can be very small (almost zero) 
! From the selection process above it can be negative
! if the positive root is far away from it. 
! As such, keep the current factor value
if(new_cov_inflate <= 0.0_r8 .OR. new_cov_inflate /= new_cov_inflate) &
   new_cov_inflate = lambda_mean

end subroutine enh_linear_bayes


!-------------------------------------------------------------------------------
!>

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


!-------------------------------------------------------------------------------
!> Routine to change the Gaussian prior into an inverse gamma (IG).
!> The Gaussian prior is represented by a mode (:= mean) and a variance; var 

subroutine change_GA_IG(mode, var, beta)

real(r8), intent(in)  :: mode, var
real(r8), intent(out) :: beta

integer :: i
real(r8) :: var_p(3), mode_p(9)   ! var and mode to the Nth power
real(r8) :: AA, BB, CC, DD, EE

! Computation savers - powers are computationally expensive
var_p(1) = var
do i=2, 3
   var_p(i) = var_p(i-1)*var
enddo

mode_p(1) = mode
do i = 2, 9
  mode_p(i) = mode_p(i-1)*mode
enddo

! Calculate the rate parameter for IG distribution.
! It's a function of both the prior mean and variannce, 
! obtained as a "real" solution to a cubic polynomial.
AA = mode_p(4) * sqrt((var_p(2) + 47.0_r8*var*mode_p(2) + 3.0_r8*mode_p(4)) / var_p(3))
BB = 75.0_r8*var_p(2)*mode_p(5)
CC = 21.0_r8*var*mode_p(7)
DD = var_p(3)*mode_p(3)
EE = (CC + BB + DD + mode_p(9) + 6.0_r8*sqrt(3.0_r8)*AA*var_p(3)) / var_p(3)

beta = (7.0_r8*var*mode + mode_p(3))/(3.0_r8*var)                               + &
       EE**(1.0_r8/3.0_r8)/3.0_r8 + mode_p(2)*(var_p(2) + 14.0_r8*var*mode_p(2) + &
       mode_p(4)) / (3.0_r8*var_p(2)*EE**(1.0_r8/3.0_r8))

end subroutine change_GA_IG


!-------------------------------------------------------------------------------
!> Write to log file what kind of inflation is being used.  

subroutine log_inflation_info(inflation_handle)

type(adaptive_inflate_type), intent(in) :: inflation_handle

character(len = 128) :: det, tadapt, sadapt, akind
write(*, *) 'flaor ', inflation_handle%flavor
select case(inflation_handle%flavor)
   case(NO_INFLATION)
      call error_handler(E_MSG, 'log_inflation_info', 'No inflation is being applied')
      return
   case (OBS_INFLATION)
      call error_handler(E_MSG, 'log_inflation_info', 'Observation space inflation was selected but is deprecated')
      return
   case (VARYING_SS_INFLATION)
      call error_handler(E_MSG, 'log_inflation_info', 'Original Anderson spatially-varying adaptive inflation selected')
   case (ENHANCED_SS_INFLATION)
      call error_handler(E_MSG, 'log_inflation_info', 'Gharamti enhanced spatially-varying adaptive inflation selected')
   case (SINGLE_SS_INFLATION)
      call error_handler(E_MSG, 'log_inflation_info', 'Single value time adaptive inflation selected')
   case (RELAXATION_TO_PRIOR_SPREAD)
      call error_handler(E_MSG, 'log_inflation_info', 'Relatation to prior spread inflation selected')
      write(string1, *) 'Relaxation coefficient is ', rtps_relaxation
      call error_handler(E_MSG, 'log_inflation_info', string1)
      return
   case default
      write(string1, *) 'Illegal inflation value for '
      call error_handler(E_ERR, 'log_inflation_info', string1, source)
end select

! For spatially varying, add other details
if(inflation_handle%deterministic) then
   call error_handler(E_MSG, 'log_inflation_info', 'Inflation is deterministic')
else
   call error_handler(E_MSG, 'log_inflation_info', 'Inflation is not deterministic')
endif
write(string1, *) 'Inflation mean lower bound is ', inflation_handle%mean_lower_bound
call error_handler(E_MSG, 'log_inflation_info', string1)

write(string1, *) 'Inflation mean upper bound is ', inflation_handle%mean_upper_bound
call error_handler(E_MSG, 'log_inflation_info', string1)

write(string1, *) 'Inflation sd lower bound is ', inflation_handle%sd_lower_bound
call error_handler(E_MSG, 'log_inflation_info', string1)

if(inflation_handle%flavor == ENHANCED_SS_INFLATION) then
   write(string1, *) 'Inflation sd max change is ', inflation_handle%sd_max_change
   call error_handler(E_MSG, 'log_inflation_info', string1)
endif

write(string1, *) 'Inflation_damping is ', inflation_handle%damping
call error_handler(E_MSG, 'log_inflation_info', string1)

end subroutine log_inflation_info

!===============================================================================
! end module adaptive_inflate_mod
!===============================================================================

!> @}

end module adaptive_inflate_mod
