! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Operations and storage required for various adaptive inflation algorithms

module adaptive_inflate_mod

!> \defgroup adaptive_inflate adaptive_inflate_mod
!> @{

use types_mod,            only : r8, PI, MISSING_R8
use time_manager_mod,     only : time_type, get_time
use utilities_mod,        only : open_file, close_file, error_handler, E_ERR, E_MSG
use random_seq_mod,       only : random_seq_type, random_gaussian, init_random_seq
use ensemble_manager_mod, only : ensemble_type, map_pe_to_task
use mpi_utilities_mod,    only : my_task_id, send_to, receive_from, send_minmax_to

implicit none
private

public :: update_inflation,                                 do_obs_inflate,           &
          update_single_state_space_inflation, update_varying_state_space_inflation,  &
          do_varying_ss_inflate,    do_single_ss_inflate,   inflate_ens,              &
          adaptive_inflate_init,    adaptive_inflate_type,                            &
                                    deterministic_inflate,  solve_quadratic,          &
          log_inflation_info,       get_minmax_task_zero,   mean_from_restart,        &
          sd_from_restart,                                                            &
          output_inf_restart,       get_inflate_mean,       get_inflate_sd,           &
          get_is_prior,             get_is_posterior,       do_ss_inflate,            &
          set_inflation_mean_copy,  set_inflation_sd_copy,  get_inflation_mean_copy,  &
          get_inflation_sd_copy,    do_rtps_inflate,        validate_inflate_options, &
          print_inflation_restart_filename,                                           &
          PRIOR_INF, POSTERIOR_INF, NO_INFLATION, OBS_INFLATION, VARYING_SS_INFLATION, &
          SINGLE_SS_INFLATION, RELAXATION_TO_PRIOR_SPREAD, ENHANCED_SS_INFLATION

character(len=*), parameter :: source = 'adaptive_inflate_mod.f90'

! Manages both observation space and state space inflation
! Handles initial values and restarts, diagnostic output, and computations
! Algorithm options at present include a single fixed observation space,
! a single fixed state space adaptive inflation,
! and a spatially-varying state space inflation that carries
! a mean and variance for the state space inflation at each point. 

! There are many length 2 variables declared to hold items. Some of them 
! refer to the concept of prior and posterior. Those that do can be referenced
! with the PRIOR_INF and POSTERIOR_INF - which MUST have values 1 and 2, do not change.

integer, parameter :: PRIOR_INF = 1
integer, parameter :: POSTERIOR_INF = 2

!>@todo Eventually the namelist options corresponding to inflation should move from 
!> filter into this module and then possibly become two different namelists so we
!> don't have these arrays of length (2).

! Encode the different inflation options
! OBS_INFLATION is currently deprecated.
! ENHANCED_SS_INFLATION is an extension of VARYING_SS_INFLATION

integer, parameter :: NO_INFLATION               = 0
integer, parameter :: OBS_INFLATION              = 1
integer, parameter :: VARYING_SS_INFLATION       = 2
integer, parameter :: SINGLE_SS_INFLATION        = 3
integer, parameter :: RELAXATION_TO_PRIOR_SPREAD = 4
integer, parameter :: ENHANCED_SS_INFLATION      = 5

! Type to keep track of information for inflation
type adaptive_inflate_type
   private
   integer               :: inflation_flavor
   integer               :: inflation_sub_flavor
   logical               :: output_restart = .false.
   logical               :: deterministic
   real(r8)              :: inflate, sd, sd_lower_bound, inf_lower_bound, inf_upper_bound
   real(r8)              :: sd_max_change
   ! Include a random sequence type in case non-deterministic inflation is used
   type(random_seq_type) :: ran_seq
   logical               :: allow_missing_in_clm
   real(r8)              :: minmax_mean(2), minmax_sd(2)
   logical               :: mean_from_restart
   logical               :: sd_from_restart
   logical               :: prior = .false.
   logical               :: posterior = .false.
   integer               :: input_mean_copy = -1 !!todo NO_COPY_PRESENT
   integer               :: input_sd_copy   = -1
end type adaptive_inflate_type

! types for updating the inflation
integer, parameter :: GHA2017 = 1
integer, parameter :: AND2009 = 2

! Module storage for writing error messages
character(len=512) :: string1, string2

! Flag indicating whether module has been initialized
logical :: initialized = .false.

! Used for precision tests in inflation update routines
real(r8), parameter    :: small = epsilon(1.0_r8)   ! threshold for avoiding NaNs/Inf


!===============================================================================

contains

!-------------------------------------------------------------------------------
!> Accessor functions for adaptive inflate type

function mean_from_restart(inflation)

type(adaptive_inflate_type) :: inflation
logical :: mean_from_restart

mean_from_restart = inflation%mean_from_restart

end function mean_from_restart


!-------------------------------------------------------------------------------
!>

function sd_from_restart(inflation)

type(adaptive_inflate_type) :: inflation
logical :: sd_from_restart

sd_from_restart = inflation%sd_from_restart

end function sd_from_restart


!-------------------------------------------------------------------------------
!>

function output_inf_restart(inflation)

type(adaptive_inflate_type) :: inflation
logical :: output_inf_restart

output_inf_restart = inflation%output_restart

end function


!-------------------------------------------------------------------------------
!>

function get_inflate_mean(inflation)

type(adaptive_inflate_type) :: inflation
real(r8)  :: get_inflate_mean

get_inflate_mean = inflation%inflate

end function


!-------------------------------------------------------------------------------
!>

function get_inflate_sd(inflation)

type(adaptive_inflate_type) :: inflation
real(r8)  :: get_inflate_sd

get_inflate_sd = inflation%sd

end function


!-------------------------------------------------------------------------------
!>

function get_is_prior(inflation)

type(adaptive_inflate_type) :: inflation
logical :: get_is_prior

get_is_prior = inflation%prior

end function get_is_prior


!-------------------------------------------------------------------------------
!>

function get_is_posterior(inflation)

type(adaptive_inflate_type) :: inflation
logical :: get_is_posterior

get_is_posterior = inflation%posterior

end function get_is_posterior


!-------------------------------------------------------------------------------
!>

function do_ss_inflate(inflation)

type(adaptive_inflate_type), intent(in) :: inflation
logical :: do_ss_inflate

if (do_single_ss_inflate(inflation) .or. &
    do_varying_ss_inflate(inflation) .or. &
    do_rtps_inflate(inflation)) then
   do_ss_inflate = .true.
else
   do_ss_inflate = .false.
endif

end function do_ss_inflate


!-------------------------------------------------------------------------------
!> Initializes an adaptive_inflate_type 

subroutine adaptive_inflate_init(inflate_handle, inf_flavor, mean_from_restart, &
   sd_from_restart, output_inflation, deterministic, & 
   inf_initial, sd_initial, inf_lower_bound, inf_upper_bound, &
   sd_lower_bound, sd_max_change, missing_ok, label)

type(adaptive_inflate_type), intent(inout) :: inflate_handle
integer,                     intent(in)    :: inf_flavor
logical,                     intent(in)    :: mean_from_restart
logical,                     intent(in)    :: sd_from_restart
logical,                     intent(in)    :: output_inflation
logical,                     intent(in)    :: deterministic
real(r8),                    intent(in)    :: inf_initial, sd_initial
real(r8),                    intent(in)    :: inf_lower_bound, inf_upper_bound
real(r8),                    intent(in)    :: sd_lower_bound
real(r8),                    intent(in)    :: sd_max_change
logical,                     intent(in)    :: missing_ok
character(len = *),          intent(in)    :: label

! random value
integer, save :: salt = 139

! Record the module version if this is first initialize call
if(.not. initialized) then
   initialized = .true.
endif

! If non-deterministic inflation is being done, need to initialize random sequence.
! use the task id number (plus 1 since they start at 0) to set the initial seed.
! NOTE: non-deterministic inflation does NOT reproduce as process count is varied!
!> this used to set the same seed for prior & posterior - add a constant value
!> in case it's called a second time.
if(.not. deterministic) then
   call init_random_seq(inflate_handle%ran_seq, my_task_id()+1 + salt)
   salt = salt + 1000 
endif

! Load up the structure first to keep track of all details of this inflation type
inflate_handle%inflation_flavor     = inf_flavor
inflate_handle%inflation_sub_flavor = inf_flavor    ! see code below
inflate_handle%output_restart       = output_inflation
inflate_handle%deterministic        = deterministic
inflate_handle%inflate              = inf_initial
inflate_handle%sd                   = sd_initial
inflate_handle%inf_lower_bound      = inf_lower_bound
inflate_handle%inf_upper_bound      = inf_upper_bound
inflate_handle%sd_lower_bound       = sd_lower_bound
inflate_handle%sd_max_change        = sd_max_change
inflate_handle%allow_missing_in_clm = missing_ok
inflate_handle%mean_from_restart    = mean_from_restart
inflate_handle%sd_from_restart      = sd_from_restart

!Overwriting the initial value of inflation with 1.0
!as in other inflation flavors (usually start from 1). 
!This is required for RTPS because inf_initial(2)
!is not really the inflation factor but rather the weighting
!parameter, say alpha: 
!RTPS: lambda = alpha * (sd_b - sd_a) / sd_a + 1
!where; sd_b (sd_a): prior (posteriro) spread
if(inf_flavor == RELAXATION_TO_PRIOR_SPREAD) inflate_handle%inflate = 1.0_r8

! Prior and posterior are intialized to false
if (trim(label)=='Prior') inflate_handle%prior = .true.
if (trim(label)=='Posterior') inflate_handle%posterior = .true.

! ENHANCED_SS_INFLATION is a subset of VARYING_SS_INFLATION. modify the main flavor here.
if (inf_flavor == ENHANCED_SS_INFLATION) then
   inflate_handle%inflation_flavor = VARYING_SS_INFLATION
endif

! Cannot support non-determistic inflation and an inf_lower_bound < 1
if(.not. deterministic .and. inf_lower_bound < 1.0_r8) then
   write(string1, *) 'Cannot have non-deterministic inflation and inf_lower_bound < 1'
   call error_handler(E_ERR, 'adaptive_inflate_init', string1, source)
endif

! give these distinctive values; if inflation is being used
! (e.g. inf_flavor > 0) then they should be set in all cases.
inflate_handle%minmax_mean(:) = MISSING_R8
inflate_handle%minmax_sd(:)   = MISSING_R8

! State space inflation is read in the IO routine read_state.

! Read type 1 (observation space inflation)
if(inf_flavor == OBS_INFLATION) then

   write(string1,  *) 'No longer supporting observation space inflation ', &
                        '(i.e. inf_flavor = 1).'
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

do_obs_inflate = (inflate_handle%inflation_flavor == OBS_INFLATION)

!>@todo I am not sure you can get here given the check in adaptive_inflate_init()
if (do_obs_inflate) then
  write(string1,  *) 'observation space inflation not suppported (i.e. inf_flavor = 1)'
  write(string2, *) 'please contact dart if you would like to use this functionality'
  call error_handler(E_ERR, 'do_obs_inflate', string1, source, text2=string2)
endif

end function do_obs_inflate


!-------------------------------------------------------------------------------
!> Returns true if this inflation type indicates varying state space inflation

function do_varying_ss_inflate(inflate_handle)

logical                                 :: do_varying_ss_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

do_varying_ss_inflate = (inflate_handle%inflation_flavor == VARYING_SS_INFLATION)

end function do_varying_ss_inflate


!-------------------------------------------------------------------------------
!> Returns true if this inflation type indicates fixed state space inflation

function do_single_ss_inflate(inflate_handle)

logical                                 :: do_single_ss_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

do_single_ss_inflate = (inflate_handle%inflation_flavor == SINGLE_SS_INFLATION)

end function do_single_ss_inflate


!-------------------------------------------------------------------------------
!> Returns true if this inflation type indicates posterior relaxion-to-prior-spread
!> (whitaker & Hamill, 2012)

function do_rtps_inflate(inflate_handle)

logical                                 :: do_rtps_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

do_rtps_inflate = (inflate_handle%inflation_flavor == RELAXATION_TO_PRIOR_SPREAD)

end function do_rtps_inflate


!-------------------------------------------------------------------------------
!> *private* accessor routine for the subtype
!>
!> Returns true if this inflation sub type indicates enhanced state space inflation
!> Moha Gharamti, 2017

function do_enhanced_ss_inflate(inflate_handle)

logical                                 :: do_enhanced_ss_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

do_enhanced_ss_inflate = ((inflate_handle%inflation_flavor == VARYING_SS_INFLATION) .and. &
                          (inflate_handle%inflation_sub_flavor == ENHANCED_SS_INFLATION))

end function do_enhanced_ss_inflate


!-------------------------------------------------------------------------------
!> Returns true if deterministic inflation is indicated

function deterministic_inflate(inflate_handle)

logical :: deterministic_inflate
type(adaptive_inflate_type), intent(in) :: inflate_handle

deterministic_inflate = inflate_handle%deterministic

end function deterministic_inflate


!-------------------------------------------------------------------------------
!> Make sure the combination of inflation options are legal

subroutine validate_inflate_options(inf_flavor, inf_damping, inf_initial_from_restart, &
                    inf_sd_initial_from_restart, inf_deterministic, inf_sd_max_change, &
                    do_prior_inflate, do_posterior_inflate, output_inflation, &
                    compute_posterior)

integer,  intent(in)    :: inf_flavor(2)
real(r8), intent(inout) :: inf_damping(2)
logical,  intent(inout) :: inf_initial_from_restart(2)
logical,  intent(inout) :: inf_sd_initial_from_restart(2)
logical,  intent(inout) :: inf_deterministic(2)
real(r8), intent(in)    :: inf_sd_max_change(2)
logical,  intent(out)   :: do_prior_inflate
logical,  intent(out)   :: do_posterior_inflate
logical,  intent(out)   :: output_inflation 
logical,  intent(in)    :: compute_posterior

integer :: i
character(len=32) :: string(2)

do_prior_inflate     = .false.
do_posterior_inflate = .false.
output_inflation     = .false.

! for error messages
string(PRIOR_INF)     = 'Prior'
string(POSTERIOR_INF) = 'Posterior'

do i = PRIOR_INF, POSTERIOR_INF
   if(inf_flavor(i) < NO_INFLATION .or. inf_flavor(i) > ENHANCED_SS_INFLATION) then
      write(string1, *) 'inf_flavor=', inf_flavor(i), ' Must be 0, 1, 2, 3, 4, or 5 '
      call error_handler(E_ERR,'validate_inflate_options', string1, source, &
                                text2='Inflation type for '//string(i))
   endif

   if(inf_damping(i) < 0.0_r8 .or. inf_damping(i) > 1.0_r8) then
      write(string1, *) 'inf_damping=', inf_damping(i), ' Must be 0.0 <= d <= 1.0'
      call error_handler(E_ERR,'validate_inflate_options', string1, source, &
                                text2='Inflation damping for '//string(i))
   endif
end do

! Check to see if state space inflation is turned on
if (inf_flavor(PRIOR_INF)     /= NO_INFLATION .and. &
    inf_flavor(PRIOR_INF)     /= OBS_INFLATION)  do_prior_inflate     = .true.
if (inf_flavor(POSTERIOR_INF) /= NO_INFLATION .and. &
    inf_flavor(POSTERIOR_INF) /= OBS_INFLATION)  do_posterior_inflate = .true.

if (do_prior_inflate .or. do_posterior_inflate) output_inflation = .true.

! Observation space inflation not currently supported
if(inf_flavor(PRIOR_INF) == OBS_INFLATION .or. &
   inf_flavor(POSTERIOR_INF) == OBS_INFLATION) &
   call error_handler(E_ERR, 'validate_inflate_options', &
   'observation space inflation (type 1) not currently supported', source, &
   text2='contact DART developers if you are interested in using it.')

! Relaxation-to-prior-spread (RTPS) is only an option for posterior inflation
if(inf_flavor(PRIOR_INF) == RELAXATION_TO_PRIOR_SPREAD) &
   call error_handler(E_ERR, 'validate_inflate_options', &
   'RTPS inflation (type 4) only supported for Posterior inflation', source)

! Cannot select posterior options if not computing posterior
if(.not. compute_posterior .and. inf_flavor(POSTERIOR_INF) /= NO_INFLATION) then
   write(string1, *) 'cannot enable posterior inflation if not computing posterior values'
   call error_handler(E_ERR,'validate_inflate_options', string1, source, &
           text2='"compute_posterior" is false; posterior inflation flavor must be 0')
endif

! RTPS needs a single parameter from namelist: inf_initial(2).  
! Do not read in any files.  Also, no damping.  
! Warn the user if they try to set different values in the namelist.
if (inf_flavor(POSTERIOR_INF) == RELAXATION_TO_PRIOR_SPREAD) then
   if (inf_initial_from_restart(POSTERIOR_INF) .or. inf_sd_initial_from_restart(POSTERIOR_INF)) &
      call error_handler(E_MSG, 'validate_inflate_options:', &
         'RTPS inflation (type 4) overrides posterior inflation restart file with value in namelist', &
         text2='posterior inflation standard deviation value not used in RTPS')

   ! Get parameter from namelist inf_initial(2), not from file
   ! inf_sd not used in this algorithm
      inf_initial_from_restart(POSTERIOR_INF) = .false.
   inf_sd_initial_from_restart(POSTERIOR_INF) = .false.

   if (.not. inf_deterministic(POSTERIOR_INF)) &
      call error_handler(E_MSG, 'validate_inflate_options:', &
              'RTPS inflation (type 4) overrides posterior inf_deterministic with .true.')
   inf_deterministic(POSTERIOR_INF) = .true.  ! this algorithm is deterministic

   if (inf_damping(POSTERIOR_INF) /= 1.0_r8) &
      call error_handler(E_MSG, 'validate_inflate_options:', &
                        'RTPS inflation (type 4) disables posterior inf_damping')
   inf_damping(POSTERIOR_INF) = 1.0_r8  ! no damping
endif

! enhanced inflation checks - this is before we set the subflavor in the structure.
if (inf_flavor(PRIOR_INF)     == ENHANCED_SS_INFLATION .or. &
    inf_flavor(POSTERIOR_INF) == ENHANCED_SS_INFLATION) then

   ! check inf_sd_max_change() for valid range
   do i=PRIOR_INF, POSTERIOR_INF
      if (inf_sd_max_change(i) < 1.0_r8 .or. inf_sd_max_change(i) > 2.0_r8) then
         write(string1, *) 'inf_sd_max_change=', inf_sd_max_change(i), &
                           ' Must be 1.0 <= X <= 2.0'
         call error_handler(E_ERR,'validate_inflate_options', string1, source, &
                                   text2='Inflation stddev max change for '//string(i))
      endif
   enddo
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

! it's possible to have MISSING_R8s in the state vector now.  
! so we need to be able to avoid changing MISSING_R8 values by inflation here.
if (inflate_handle%allow_missing_in_clm) then
   if (any(ens == MISSING_R8)) return
endif

if(inflate_handle%deterministic) then

   if ( do_rtps_inflate(inflate_handle)) then
      if ( .not. present(fsprd) .or. .not. present(asprd)) then 
         write(string1, *) 'missing arguments for RTPS inflation, should not happen'
         call error_handler(E_ERR,'inflate_ens',string1,source) 
      endif 
      ! only inflate if spreads are > 0
      if ( asprd .gt. 0.0_r8 .and. fsprd .gt. 0.0_r8) &
          inflate = 1.0_r8 + inflate * ((fsprd-asprd) / asprd) 
          ens = ens * inflate + mean * (1.0_r8 - inflate)
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

! select which method to update with
if (do_enhanced_ss_inflate(inflate_handle)) then
   inf_type = GHA2017
else
   inf_type = AND2009
endif

! Use bayes theorem to update
call bayes_cov_inflate(ens_size, inf_type, prior_mean, prior_var, obs, obs_var, inflate, &
   inflate_sd, gamma_corr, inflate_handle%sd_lower_bound, inflate_handle%sd_max_change, &
   new_inflate, new_inflate_sd)

! Make sure inflate satisfies constraints
inflate = new_inflate
if(inflate < inflate_handle%inf_lower_bound) inflate = inflate_handle%inf_lower_bound
if(inflate > inflate_handle%inf_upper_bound) inflate = inflate_handle%inf_upper_bound

! Make sure sd satisfies constraints
inflate_sd = new_inflate_sd
if(inflate_sd < inflate_handle%sd_lower_bound) inflate_sd = inflate_handle%sd_lower_bound

end subroutine update_inflation

!-------------------------------------------------------------------------------
!> Computes updated inflation mean and inflation sd for single state space inflation

subroutine update_single_state_space_inflation(inflate, inflate_mean, inflate_sd, &
   ss_inflate_base, orig_obs_prior_mean, orig_obs_prior_var, obs, obs_err_var, &
   ens_size, inflate_only)

type(adaptive_inflate_type), intent(in)    :: inflate
real(r8),                   intent(inout)  :: inflate_mean
real(r8),                   intent(inout)  :: inflate_sd
real(r8),                   intent(in)     :: ss_inflate_base
real(r8),                   intent(in)     :: orig_obs_prior_mean
real(r8),                   intent(in)     :: orig_obs_prior_var
real(r8),                   intent(in)     :: obs
real(r8),                   intent(in)     :: obs_err_var
integer,                    intent(in)     :: ens_size
logical,                    intent(in)     :: inflate_only

real(r8) :: gamma, ens_var_deflate, r_var, r_mean

! If either inflation or sd is not positive, not really doing inflation
if(inflate_mean <= 0.0_r8 .or. inflate_sd <= 0.0_r8) return

! For case with single spatial inflation, use gamma = 1.0_r8
gamma = 1.0_r8
! Deflate the inflated variance; required for efficient single pass
! This is one of many places that assumes linear state/obs relation
! over range of ensemble; Essentially, we are removing the inflation
! which has already been applied in filter to see what inflation should
! have been needed.
ens_var_deflate = orig_obs_prior_var / &
   (1.0_r8 + gamma*(sqrt(ss_inflate_base) - 1.0_r8))**2

! If this is inflate_only (i.e. posterior) remove impact of this obs.
! This is simulating independent observation by removing its impact.
if(inflate_only .and. &
      ens_var_deflate               > small .and. &
      obs_err_var                   > small .and. &
      obs_err_var - ens_var_deflate > small ) then
   r_var = 1.0_r8 / (1.0_r8 / ens_var_deflate - 1.0_r8 / obs_err_var)
   r_mean = r_var *(orig_obs_prior_mean / ens_var_deflate - obs / obs_err_var)
else
   r_var = ens_var_deflate
   r_mean = orig_obs_prior_mean
endif

! Update the inflation mean value and standard deviation
call update_inflation(inflate, inflate_mean, inflate_sd, &
   r_mean, r_var, ens_size, obs, obs_err_var, gamma)

end subroutine update_single_state_space_inflation

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
                 lambda_mean, lambda_sd, gamma_corr, sd_lower_bound_in, &
                 sd_max_change_in, new_cov_inflate, new_cov_inflate_sd)

integer , intent(in)  :: ens_size, inf_type
real(r8), intent(in)  :: x_p, sigma_p_2, y_o, sigma_o_2, lambda_mean, lambda_sd
real(r8), intent(in)  :: gamma_corr, sd_lower_bound_in, sd_max_change_in
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

if (inf_type == AND2009) then

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
   if(abs(lambda_sd - sd_lower_bound_in) <= TINY(0.0_r8)) then
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

else if (inf_type == GHA2017) then

   ! Transform Gaussian prior to Inverse Gamma
   call change_GA_IG(lambda_mean, lambda_sd_2, rate)

   ! Approximate with Taylor series for likelihood term
   call enh_linear_bayes(dist_2, sigma_p_2, sigma_o_2,lambda_mean, &
                    gamma_corr, ens_size, rate, new_cov_inflate)

   ! Bail out to save cost when lower bound is reached on lambda standard deviation
   ! See comment in Anderson case for why we use abs and TINY for this comparison.
   if(abs(lambda_sd - sd_lower_bound_in) <= TINY(0.0_r8)) then
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
      if ( new_cov_inflate_sd > sd_max_change_in*lambda_sd .OR. &
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

subroutine log_inflation_info(inflation_handle, mype, label, single_file)

type(adaptive_inflate_type), intent(in) :: inflation_handle
integer,                     intent(in) :: mype
character(len = *),          intent(in) :: label
logical,                     intent(in) :: single_file

character(len = 128) :: det, tadapt, sadapt, akind, from

! nothing to do if not task 0
if (mype /= 0) return

! if inflation is off, say so and return now
if (inflation_handle%inflation_flavor <= NO_INFLATION) then
   call error_handler(E_MSG, trim(label) // ' inflation:', 'None', source)
   return
endif

! construct english language version of our complicated combinations
! of inflation-related parameters.
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
   if (inflation_handle%sd_lower_bound < inflation_handle%minmax_sd(2) .or. &
       inflation_handle%inflation_sub_flavor == ENHANCED_SS_INFLATION) then
      tadapt = trim(tadapt) // ' time-rate adaptive,'
   endif
else
  tadapt = ' time-constant,'
endif
if (inflation_handle%inflation_sub_flavor == ENHANCED_SS_INFLATION) then
  tadapt = ' enhanced' //trim(tadapt)
endif

select case(inflation_handle%inflation_flavor)
   case (OBS_INFLATION)
      sadapt = ' (deprecated),'
      akind = ' observation-space'
   case (VARYING_SS_INFLATION)
      sadapt = ' spatially-varying,'
      akind = ' state-space '
   case (SINGLE_SS_INFLATION)
      sadapt = ' spatially-constant,'
      akind = ' state-space'
   case (RELAXATION_TO_PRIOR_SPREAD)
      tadapt = ' time-adaptive,'    ! IS THIS TRUE??
      sadapt = ' spatially-varying relaxation-to-prior-spread,'
      akind = ' state-space'
   case default
      write(string1, *) 'Illegal inflation value for ', label
      call error_handler(E_ERR, 'adaptive_inflate_init', string1, source)
end select

! say what basic kind of inflation was selected.
write(string1, '(4A)') trim(det), trim(tadapt), trim(sadapt), trim(akind)
call error_handler(E_MSG, trim(label) // ' inflation:', string1, source)

! print out details about the type of inflation selected.
! inflation flavor 2 has min/max values if from restart,
! flavor 2 from namelist and flavor 3 only a single value

! do this twice - for mean and sd

! combination file, individual file or namelist
call set_from_string(inflation_handle%mean_from_restart, single_file, from)
if (nvalues_to_log(inflation_handle, from) == 1) then
   write(string1,  '(A, F8.3)') &
         'inf mean   '//trim(from)//', value: ', inflation_handle%minmax_mean(1)
else
   write(string1,  '(A, 2F8.3)') &
         'inf mean   '//trim(from)//', min/max values: ', inflation_handle%minmax_mean
endif
call error_handler(E_MSG, trim(label) // ' inflation:', string1,  source)

call set_from_string(inflation_handle%sd_from_restart, single_file, from)
if (nvalues_to_log(inflation_handle, from) == 1) then
   write(string1,  '(A, F8.3)') &
         'inf stddev '//trim(from)//', value: ', inflation_handle%minmax_sd(1)
else
   write(string1,  '(A, 2F8.3)') &
         'inf stddev '//trim(from)//', min/max values: ', inflation_handle%minmax_sd
endif
call error_handler(E_MSG, trim(label) // ' inflation:', string1,  source)

if (inflation_handle%inflation_sub_flavor == ENHANCED_SS_INFLATION) then
   write(string1, '(A, F8.3)') &
            'inf stddev max change: ', inflation_handle%sd_max_change
   call error_handler(E_MSG, trim(label) // ' inflation:', string1, source)
endif

end subroutine log_inflation_info


!-------------------------------------------------------------------------------
!>

subroutine set_from_string(from_restart, single_file, from_string)
logical,          intent(in)  :: from_restart
logical,          intent(in)  :: single_file
character(len=*), intent(out) :: from_string

if (from_restart) then
   if (single_file) then
      from_string = 'variable from input file'
   else
      from_string = 'restart file'
   endif
else
   from_string = 'from namelist'
endif

end subroutine set_from_string


!-------------------------------------------------------------------------------
!>

function nvalues_to_log(inflation_handle, from_string)
type(adaptive_inflate_type), intent(in) :: inflation_handle
character(len=*),            intent(in) :: from_string
integer :: nvalues_to_log

if ((inflation_handle%inflation_flavor == SINGLE_SS_INFLATION) .or. &
    (inflation_handle%inflation_flavor == VARYING_SS_INFLATION .and. &
     from_string == 'from namelist')) then
   nvalues_to_log = 1
else
   nvalues_to_log = 2
endif

end function nvalues_to_log


!-------------------------------------------------------------------------------
!>

subroutine print_inflation_restart_filename(inflation_handle, fname, which)
type(adaptive_inflate_type), intent(in) :: inflation_handle
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: which

character(len=32) :: label

if (inflation_handle%prior) then
   label = "Prior"
elseif (inflation_handle%posterior) then
   label = "Posterior"
else
   label = ""
endif

write(string1,*) trim(which)//' read from restart file: ' // trim(fname)
call error_handler(E_MSG, trim(label) // ' inflation:', trim(string1), source)

end subroutine print_inflation_restart_filename


!-------------------------------------------------------------------------------
!> Collect the min and max of inflation on task 0
!> this block handles communicating the min/max local values to PE 0
!> if running with MPI, or just sets the min/max directly if reading
!> from a namelist.

subroutine get_minmax_task_zero(inflation_handle, ens_handle, &
                       ss_inflate_index, ss_inflate_sd_index)

type(adaptive_inflate_type), intent(inout) :: inflation_handle
type(ensemble_type),         intent(in)    :: ens_handle
integer,                     intent(in)    :: ss_inflate_index
integer,                     intent(in)    :: ss_inflate_sd_index

real(r8) :: minmax_mean(2), minmax_sd(2), global_val(2)

! if not using inflation, return now
if (inflation_handle%inflation_flavor <= NO_INFLATION) return

if (inflation_handle%mean_from_restart) then

   ! find min and max on each processor
   minmax_mean(1) = minval(ens_handle%copies(ss_inflate_index, :))
   minmax_mean(2) = maxval(ens_handle%copies(ss_inflate_index, :))

   ! collect on pe 0
   call send_minmax_to(minmax_mean, map_pe_to_task(ens_handle, 0), global_val)
   if (ens_handle%my_pe == 0) inflation_handle%minmax_mean = global_val

else 
   inflation_handle%minmax_mean = inflation_handle%inflate
endif

if (inflation_handle%sd_from_restart) then

   ! find min and max on each processor
   minmax_sd(1) = minval(ens_handle%copies(ss_inflate_sd_index, :))
   minmax_sd(2) = maxval(ens_handle%copies(ss_inflate_sd_index, :))

   ! collect on pe 0
   call send_minmax_to(minmax_sd, map_pe_to_task(ens_handle, 0), global_val)
   if (ens_handle%my_pe == 0) inflation_handle%minmax_sd = global_val
else
   inflation_handle%minmax_sd = inflation_handle%sd 
endif

end subroutine get_minmax_task_zero


!-------------------------------------------------------------------------------
!>

subroutine set_inflation_mean_copy(inflation_handle, c)
type(adaptive_inflate_type), intent(inout) :: inflation_handle
integer,                     intent(in)    :: c

inflation_handle%input_mean_copy = c

end subroutine set_inflation_mean_copy


!-------------------------------------------------------------------------------
!>

subroutine set_inflation_sd_copy(inflation_handle, c)
type(adaptive_inflate_type), intent(inout) :: inflation_handle
integer,                     intent(in)    :: c

inflation_handle%input_sd_copy = c

end subroutine set_inflation_sd_copy


!-------------------------------------------------------------------------------
!>

function get_inflation_mean_copy(inflation_handle) result (c)
type(adaptive_inflate_type), intent(in) :: inflation_handle
integer :: c

c = inflation_handle%input_mean_copy

end function get_inflation_mean_copy


!-------------------------------------------------------------------------------
!>

function get_inflation_sd_copy(inflation_handle) result (c)
type(adaptive_inflate_type), intent(in) :: inflation_handle
integer :: c

c = inflation_handle%input_sd_copy

end function get_inflation_sd_copy


!===============================================================================
! end module adaptive_inflate_mod
!===============================================================================

!> @}

end module adaptive_inflate_mod
