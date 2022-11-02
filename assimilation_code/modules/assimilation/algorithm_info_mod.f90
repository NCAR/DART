! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module algorithm_info_mod

use types_mod, only : r8

use obs_def_mod, only : obs_def_type, get_obs_def_type_of_obs, get_obs_def_error_variance
use obs_kind_mod, only : get_quantity_for_type_of_obs

implicit none
private

integer, parameter :: NORMAL_PRIOR = 1
integer, parameter :: BOUNDED_NORMAL_RH_PRIOR = 2

public :: obs_error_info, probit_dist_info, obs_inc_info, &
          NORMAL_PRIOR, BOUNDED_NORMAL_RH_PRIOR

! Provides routines that give information about details of algorithms for 
! observation error sampling, observation increments, and the transformations
! for regression and inflation in probit space. 
! For now, it is convenient to have these in a single module since several
! users will be developing their own problem specific versions of these
! subroutines. This will avoid constant merge conflicts as other parts of the
! assimilation code are updated.

contains

!-------------------------------------------------------------------------
subroutine obs_error_info(obs_def, error_variance, bounded, bounds)

! Computes information needed to compute error sample for this observation
! This is called by perfect_model_obs when generating noisy obs
type(obs_def_type), intent(in)  :: obs_def
real(r8),           intent(out) :: error_variance
logical,            intent(out) :: bounded(2)
real(r8),           intent(out) :: bounds(2)

integer :: obs_type, obs_kind

! Get the kind of the observation
obs_type = get_obs_def_type_of_obs(obs_def)
obs_kind = get_quantity_for_type_of_obs(obs_type)

write(*, *) 'obs_error_info type and kind ', obs_type, obs_kind

! Get the default error variance
error_variance = get_obs_def_error_variance(obs_def)

! Do some computation with the obs_kind
bounded(1) = .false.;     bounded(2) = .false.
bounds(1) = 0.0_r8;

end subroutine obs_error_info


!-------------------------------------------------------------------------


subroutine probit_dist_info(kind, is_state, is_inflation, dist_type, &
   bounded, bounds)

! Computes the details of the probit transform for initial experiments
! with Molly 

integer,  intent(in)  :: kind
logical,  intent(in)  :: is_state      ! True for state variable, false for obs
logical,  intent(in)  :: is_inflation  ! True for inflation transform
integer,  intent(out) :: dist_type
logical,  intent(out) :: bounded(2)
real(r8), intent(out) :: bounds(2)

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
! bounded(1) = .true.;  bounded(2) = .true.
! bounds(1)  = 0.0_r8;  bounds(2)  = 1.0_r8

! This logic is consistent with Quantile Regression paper square experiments
if(is_inflation) then
   dist_type = BOUNDED_NORMAL_RH_PRIOR
   bounded = .false.
elseif(is_state) then
   dist_type = BOUNDED_NORMAL_RH_PRIOR
   bounded = .false.
else
   dist_type = BOUNDED_NORMAL_RH_PRIOR
   bounded(1) = .true.;    bounded(2) = .false.
   bounds(1)  = 0.0_r8
endif

end subroutine probit_dist_info

!------------------------------------------------------------------------


subroutine obs_inc_info(obs_kind, filter_kind, rectangular_quadrature, gaussian_likelihood_tails, &
   sort_obs_inc, spread_restoration, bounded, bounds)

integer,  intent(in)  :: obs_kind
integer,  intent(out) :: filter_kind
logical,  intent(out) :: rectangular_quadrature, gaussian_likelihood_tails
logical,  intent(out) :: sort_obs_inc
logical,  intent(out) :: spread_restoration
logical,  intent(out) :: bounded(2)
real(r8), intent(out) :: bounds(2)

! Temporary approach for setting the details of how to assimilate this observation
! This example is designed to reproduce the squared forward operator results from paper

filter_kind = 101
sort_obs_inc = .false.
spread_restoration = .false.
bounded(1) = .true.;   bounded(2) = .false.
bounds(1) = 0.0_r8;

! Only need to set these two for options on old RHF implementation
! rectangular_quadrature = .true.
! gaussian_likelihood_tails = .false.

end subroutine obs_inc_info

!------------------------------------------------------------------------

end module algorithm_info_mod
