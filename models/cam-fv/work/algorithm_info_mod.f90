! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module algorithm_info_mod

use types_mod, only : r8, i8, missing_r8

use obs_def_mod, only : obs_def_type, get_obs_def_type_of_obs, get_obs_def_error_variance
use obs_kind_mod, only : get_quantity_for_type_of_obs

! Get the QTY definitions that are needed (aka kind)
use obs_kind_mod, only : QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT, QTY_SURFACE_PRESSURE, &
                         QTY_TEMPERATURE, QTY_SPECIFIC_HUMIDITY, QTY_CLOUD_LIQUID_WATER, &
                         QTY_CLOUD_ICE, QTY_GPSRO

use assim_model_mod, only : get_state_meta_data
use location_mod, only    : location_type

use distribution_params_mod, only : NORMAL_DISTRIBUTION, BOUNDED_NORMAL_RH_DISTRIBUTION, &
   GAMMA_DISTRIBUTION, BETA_DISTRIBUTION, LOG_NORMAL_DISTRIBUTION, UNIFORM_DISTRIBUTION,  &
   PARTICLE_FILTER_DISTRIBUTION

implicit none
private

! Defining parameter strings for different observation space filters
! For now, retaining backwards compatibility in assim_tools_mod requires using
! these specific integer values and there is no point in using these in assim_tools.
! That will change if backwards compatibility is removed in the future.
integer, parameter :: EAKF               = 1
integer, parameter :: ENKF               = 2
integer, parameter :: UNBOUNDED_RHF      = 8
integer, parameter :: GAMMA_FILTER       = 11
integer, parameter :: BOUNDED_NORMAL_RHF = 101 

public :: obs_error_info, probit_dist_info, obs_inc_info, &
          EAKF, ENKF, BOUNDED_NORMAL_RHF, UNBOUNDED_RHF, GAMMA_FILTER

! Provides routines that give information about details of algorithms for 
! observation error sampling, observation increments, and the transformations
! for regression and inflation in probit space. 
! For now, it is convenient to have these in a single module since several
! users will be developing their own problem specific versions of these
! subroutines. This will avoid constant merge conflicts as other parts of the
! assimilation code are updated.

contains

!-------------------------------------------------------------------------
subroutine obs_error_info(obs_def, error_variance, &
   bounded_below, bounded_above, lower_bound, upper_bound)

! Computes information needed to compute error sample for this observation
! This is called by perfect_model_obs when generating noisy obs
type(obs_def_type), intent(in)  :: obs_def
real(r8),           intent(out) :: error_variance
logical,            intent(out) :: bounded_below, bounded_above
real(r8),           intent(out) :: lower_bound, upper_bound

integer     :: obs_type, obs_kind
integer(i8) :: state_var_index
type(location_type) :: temp_loc

! Get the kind of the observation
obs_type = get_obs_def_type_of_obs(obs_def)
! If it is negative, it is an identity obs
if(obs_type < 0) then
   state_var_index = -1 * obs_type
   call get_state_meta_data(state_var_index, temp_loc, obs_kind)
else
   obs_kind = get_quantity_for_type_of_obs(obs_type)
endif

! Get the default error variance
error_variance = get_obs_def_error_variance(obs_def)

! Set the observation error details for each type of quantity
bounded_below  = .false.;    bounded_above  = .false.
lower_bound = missing_r8;    upper_bound = missing_r8

end subroutine obs_error_info


!-------------------------------------------------------------------------


subroutine probit_dist_info(kind, is_state, is_inflation, dist_type, &
   bounded_below, bounded_above, lower_bound, upper_bound)

! Computes the details of the probit transform for initial experiments
! with Molly 

integer,  intent(in)  :: kind
logical,  intent(in)  :: is_state      ! True for state variable, false for obs
logical,  intent(in)  :: is_inflation  ! True for inflation transform
integer,  intent(out) :: dist_type
logical,  intent(out) :: bounded_below, bounded_above
real(r8), intent(out) :: lower_bound,   upper_bound

! Have input information about the kind of the state or observation being transformed
! along with additional logical info that indicates whether this is an observation
! or state variable and about whether the transformation is being done for inflation
! or for regress. 
! Need to select the appropriate transform. At present, options are NORMAL_PRIOR
! which does nothing or BOUNDED_NORMAL_RH_PRIOR. 
! If the BNRH is selected then information about the bounds must also be set.
! The two dimensional logical array 'bounded' is set to false for no bounds and true
! for bounded. the first element of the array is for the lower bound, the second for the upper.
! If bounded is chosen, the corresponding bound value(s) must be set in the two dimensional 
! real array 'bounds'.
! For example, if my_state_kind corresponds to a sea ice fraction then an appropriate choice
! would be:
! bounded_below = .true.;  bounded_above = .true.
! lower_bound  = 0.0_r8;   upper_bounds  = 1.0_r8

! In the long run, may not have to have separate controls for each of the input possibilities
! However, for now these are things that need to be explored for science understanding

select case(kind)
   case(QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT, QTY_SURFACE_PRESSURE, QTY_TEMPERATURE)
   dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
!   dist_type = NORMAL_PRIOR
   bounded_below = .false.;       bounded_above = .false.
   lower_bound = missing_r8;      upper_bound = missing_r8

!--------------
   case(QTY_SPECIFIC_HUMIDITY)
   dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
!   dist_type = NORMAL_PRIOR
!   bounded_below = .false.;     bounded_above = .false.
   bounded_below = .true.;     bounded_above = .true.
   lower_bound = 0.0_r8;       upper_bound = 1.0_r8

!--------------
   case(QTY_CLOUD_LIQUID_WATER, QTY_CLOUD_ICE)
   dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
!   dist_type = NORMAL_PRIOR
!   bound_below = .false.;     bounded_above = .false.
   bounded_below = .true.;     bounded_above = .false.
   lower_bound = 0.0_r8;       upper_bound = missing_r8

!--------------
   case(QTY_GPSRO)
   dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
!   dist_type = NORMAL_PRIOR
!   bounded_below = .false.;     bounded_above = .false.
   bounded_below = .true.;     bounded_above = .false.
   lower_bound = 0.0_r8;       upper_bound = missing_r8

!--------------
   case DEFAULT
      write(*, *) 'Unexpected QTY in algorithm_info_mod ', kind
      stop
end select

end subroutine probit_dist_info

!------------------------------------------------------------------------


subroutine obs_inc_info(obs_kind, filter_kind, rectangular_quadrature, gaussian_likelihood_tails, &
   sort_obs_inc, spread_restoration, bounded_below, bounded_above, lower_bound, upper_bound)

integer,  intent(in)  :: obs_kind
integer,  intent(inout) :: filter_kind
logical,  intent(inout) :: rectangular_quadrature, gaussian_likelihood_tails
logical,  intent(inout) :: sort_obs_inc
logical,  intent(inout) :: spread_restoration
logical,  intent(inout) :: bounded_below, bounded_above
real(r8), intent(inout) :: lower_bound,  upper_bound

! The information arguments are all intent (inout). This means that if they are not set
! here, they retain the default values from the assim_tools_mod namelist. Bounds don't exist 
! in that namelist, so default values are set in assim_tools_mod just before the call to here.

! Temporary approach for setting the details of how to assimilate this observation
! This example is designed to reproduce the squared forward operator results from paper

select case(obs_kind)
   case(QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT, QTY_SURFACE_PRESSURE, QTY_TEMPERATURE)
      ! Set the observation increment details for each type of quantity
      filter_kind = BOUNDED_NORMAL_RHF
      bounded_below = .false.;     bounded_above = .false.
      lower_bound = missing_r8;    upper_bound = missing_r8

   case(QTY_GPSRO)
      filter_kind = BOUNDED_NORMAL_RHF
      bounded_below = .true.;     bounded_above = .false.
!      bounded_below = .false.;     bounded_above = .false.
      lower_bound = 0.0_r8;       upper_bound = missing_r8

   case(QTY_SPECIFIC_HUMIDITY)
      filter_kind = BOUNDED_NORMAL_RHF
      bounded_below = .true.;     bounded_above = .true.
!      bounded_below = .false.;     bounded_above = .false.
      lower_bound = 0.0_r8;       upper_bound = 1.0_r8

   case DEFAULT
      write(*, *) 'Unexpected QTY in algorithm_info_mod ', obs_kind
      stop
end select

! Default settings for now for Icepack and tracer model tests
sort_obs_inc = .false.
spread_restoration = .false.

! Only need to set these two for options the original RHF implementation
!!!rectangular_quadrature = .true.
!!!gaussian_likelihood_tails = .false.

end subroutine obs_inc_info

!------------------------------------------------------------------------

end module algorithm_info_mod
