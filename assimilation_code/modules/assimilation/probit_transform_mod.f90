! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! A variety of PDFs, CDFs, quantile functions and other tools for working with distributions
! to implement quantile conserving filters in observation space and regression in quantile space.

module probit_transform_mod

use types_mod, only : r8, missing_r8

use sort_mod,  only : index_sort

use utilities_mod, only : E_ERR, error_handler, do_nml_file, do_nml_term, nmlfileunit, &
                          find_namelist_in_file, check_namelist_read

use distribution_params_mod, only : distribution_params_type, deallocate_distribution_params, &
                                    NORMAL_DISTRIBUTION, BOUNDED_NORMAL_RH_DISTRIBUTION, &
                                    GAMMA_DISTRIBUTION, BETA_DISTRIBUTION,               &
                                    LOG_NORMAL_DISTRIBUTION, UNIFORM_DISTRIBUTION,       &
                                    PARTICLE_FILTER_DISTRIBUTION

use normal_distribution_mod, only : normal_cdf, inv_std_normal_cdf

use gamma_distribution_mod, only : gamma_cdf_params, inv_gamma_cdf_params, &
                                   set_gamma_params_from_ens 

use beta_distribution_mod,  only : beta_cdf_params, inv_beta_cdf_params, &
                                   set_beta_params_from_ens

use bnrh_distribution_mod,    only : bnrh_cdf_initialized_vector, bnrh_cdf_params, &
                                     inv_bnrh_cdf_params, get_bnrh_sd

implicit none
private

public :: transform_to_probit, transform_from_probit, transform_all_to_probit, &
   transform_all_from_probit

character(len=512)     :: errstring
character(len=*), parameter :: source = 'probit_transform_mod.f90'

! Global to indicate module has been initialized
logical :: module_initialized = .false.

! Namelist with default value
! Logical to fix bounds violations for bounded_normal_rh
logical :: fix_bound_violations = .false.
! Should we use a logit transform instead of the default probit transform
logical :: use_logit_instead_of_probit = .false.
! Set to true to do a check of the probit to/from transforms for inverse accuracy
logical :: do_inverse_check = .false.

namelist /probit_transform_nml/ fix_bound_violations, &
          use_logit_instead_of_probit, do_inverse_check

contains

!------------------------------------------------------------------------

subroutine transform_all_to_probit(ens_size, num_vars, state_ens, distribution_type, &
   p, probit_ens, use_input_p, bounded_below, bounded_above, lower_bound, upper_bound)

integer, intent(in)                           :: ens_size
integer, intent(in)                           :: num_vars
real(r8), intent(in)                          :: state_ens(:, :)
integer, intent(in)                           :: distribution_type(num_vars)
type(distribution_params_type), intent(inout) :: p(num_vars)
real(r8), intent(out)                         :: probit_ens(:, :)
logical, intent(in)                           :: use_input_p
logical, intent(in)                           :: bounded_below, bounded_above
real(r8), intent(in)                          :: lower_bound,   upper_bound


! NOTE THAT WILL MAKE HELEN CRAZY: THIS WORKS WITH THE INPUT CALLING ARGUMENTS FOR STATE_ENS AND
! PROBIT_ENS BEING THE SAME. A TEMP IS USED TO AVOID OVERWRITING ISSUES. IS THIS YUCKY?

! Note that the input and output arrays may have extra copies (first subscript). Passing sections of a
! leading index could be inefficient for time and storage, so avoiding that for now.

! Assumes that the bounds are the same for any variables that are BNRH for now
! The bounds variables are not used for the normal case or the case where the input p is used

integer  :: i
real(r8) :: temp_ens(ens_size)

do i = 1, num_vars
   call transform_to_probit(ens_size, state_ens(1:ens_size, i), distribution_type(i), &
      p(i), temp_ens, use_input_p, bounded_below, bounded_above, lower_bound, upper_bound)
   probit_ens(1:ens_size, i) = temp_ens
end do

end subroutine transform_all_to_probit

!------------------------------------------------------------------------

subroutine transform_to_probit(ens_size, state_ens_in, distribution_type, p, &
   probit_ens, use_input_p, bounded_below, bounded_above, lower_bound, upper_bound)

integer, intent(in)                           :: ens_size
real(r8), intent(in)                          :: state_ens_in(ens_size)
integer, intent(in)                           :: distribution_type
type(distribution_params_type), intent(inout) :: p
real(r8), intent(out)                         :: probit_ens(ens_size)
logical, intent(in)                           :: use_input_p
logical, intent(in)                           :: bounded_below, bounded_above
real(r8), intent(in)                          :: lower_bound,   upper_bound

real(r8) :: state_ens(ens_size)
real(r8) :: probit_ens_temp(ens_size), state_ens_temp(ens_size), diff(ens_size)
type(distribution_params_type) :: p_temp
integer :: i

! If not initialized, read in the namelist
if(.not. module_initialized) call initialize_probit_transform

! Fix bounds violations if requested
if(fix_bound_violations) then
   do i = 1, ens_size
      state_ens(i) = fix_bounds(state_ens_in(i), bounded_below, bounded_above, &
         lower_bound, upper_bound) 
   end do
else
   state_ens = state_ens_in
endif

! Set the type of the distribution in the parameters defined type
p%distribution_type = distribution_type

if(p%distribution_type == NORMAL_DISTRIBUTION) then 
   ! No transformation is done for a normal
   probit_ens = state_ens
elseif(p%distribution_type == LOG_NORMAL_DISTRIBUTION) then 
   call to_probit_log_normal(ens_size, state_ens, probit_ens)
elseif(p%distribution_type == UNIFORM_DISTRIBUTION) then 
   call to_probit_uniform(ens_size, state_ens, p, probit_ens, use_input_p, lower_bound, upper_bound)
elseif(p%distribution_type == GAMMA_DISTRIBUTION) then 
   call to_probit_gamma(ens_size, state_ens, p, probit_ens, use_input_p, &
      bounded_below, bounded_above, lower_bound, upper_bound)
elseif(p%distribution_type == BETA_DISTRIBUTION) then 
   call to_probit_beta(ens_size, state_ens, p, probit_ens, use_input_p, &
      lower_bound, upper_bound)
elseif(p%distribution_type == BOUNDED_NORMAL_RH_DISTRIBUTION) then
   call to_probit_bounded_normal_rh(ens_size, state_ens, p, probit_ens, &
      use_input_p, bounded_below, bounded_above, lower_bound, upper_bound)

!----------------------------------------------------------------------------------
! The following code block tests that the to/from probit calls are nearly inverse
! for all of the calls made during an assimilation
   if(do_inverse_check) then
      if(.not. use_input_p) then
         call to_probit_bounded_normal_rh(ens_size, state_ens, p_temp, probit_ens_temp, &
            use_input_p, bounded_below, bounded_above, lower_bound, upper_bound)
         call from_probit_bounded_normal_rh(ens_size, probit_ens_temp, p_temp, state_ens_temp)
         diff = state_ens - state_ens_temp
         if(abs(maxval(diff)) > 1.0e-8_r8) then
            write(*, *) 'Maximum allowed value of probit to/from difference exceeded'
            write(*, *) 'Location of minimum ensemble member ', minloc(state_ens)
            write(*, *) 'Location of maximum ensemble member ', maxloc(state_ens)
            do i = 1, ens_size
               write(*, *) i, state_ens(i), state_ens_temp(i), diff(i)
            enddo
            stop
         endif
      endif
   
      if(use_input_p) then
         call to_probit_bounded_normal_rh(ens_size, state_ens, p, probit_ens_temp, &
            use_input_p, bounded_below, bounded_above, lower_bound, upper_bound)
         call from_probit_bounded_normal_rh(ens_size, probit_ens_temp, p, state_ens_temp)
         diff = state_ens - state_ens_temp
         if(abs(maxval(diff)) > 1.0e-8_r8) then
            write(*, *) 'Maximum allowed value of probit to/from difference for input p exceeded'
            write(*, *) 'Location of minimum ensemble member ', minloc(state_ens)
            write(*, *) 'Location of maximum ensemble member ', maxloc(state_ens)
            do i = 1, ens_size
               write(*, *) i, state_ens(i), state_ens_temp(i), diff(i)
            enddo
            stop
         endif
      
      endif
   endif
!----------------------------------------------------------------------------------


!!!elseif(p%distribution_type == PARTICLE_FILTER_DISTRIBUTION) then
   !!!call to_probit_particle(ens_size, state_ens, p, probit_ens, use_input_p, &
       !!!bounded_below, bounded_above, lower_bound, upper_bound)
else
   write(errstring, *) 'Illegal distribution type', p%distribution_type
   call error_handler(E_ERR, 'transform_to_probit', errstring, source)
endif

end subroutine transform_to_probit

!------------------------------------------------------------------------

subroutine to_probit_log_normal(ens_size, state_ens, probit_ens)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: state_ens(ens_size)
real(r8), intent(out)                :: probit_ens(ens_size)

! Taking the logarithm leads directly to a normal distribution
! This normal may not be standard normal, but needs no further adjustment like 
! the regular normal
probit_ens = log(state_ens)

end subroutine to_probit_log_normal

!------------------------------------------------------------------------

subroutine to_probit_uniform(ens_size, state_ens, p, probit_ens, use_input_p, &
   lower_bound_in, upper_bound_in)

integer, intent(in)                           :: ens_size
real(r8), intent(in)                          :: state_ens(ens_size)
type(distribution_params_type), intent(inout) :: p
real(r8), intent(out)                         :: probit_ens(ens_size)
logical, intent(in)                           :: use_input_p
real(r8), intent(in)                          :: lower_bound_in, upper_bound_in

real(r8) :: lower_bound, upper_bound, d_range, quantile
integer :: i

! There is no distribution_mod for uniform at the moment so params setup is done here
if(use_input_p) then
   lower_bound = p%lower_bound
   upper_bound = p%upper_bound
else
   lower_bound = lower_bound_in
   upper_bound = upper_bound_in
   ! Save the bounds in the distribution_params_type
   p%lower_bound = lower_bound
   p%upper_bound = upper_bound
endif

d_range = upper_bound - lower_bound
do i = 1, ens_size
   ! Transform to quantile; U(lower_bound, upper_bound) to U(0, 1)
   quantile = (state_ens(i) - lower_bound) / d_range
   ! Transform to probit/logit space 
   probit_ens(i) = probit_or_logit_transform(quantile)
end do

end subroutine to_probit_uniform

!------------------------------------------------------------------------

subroutine to_probit_gamma(ens_size, state_ens, p, probit_ens, use_input_p, &
   bounded_below, bounded_above, lower_bound, upper_bound)

integer,  intent(in)                          :: ens_size
real(r8), intent(in)                          :: state_ens(ens_size)
type(distribution_params_type), intent(inout) :: p
real(r8), intent(out)                         :: probit_ens(ens_size)
logical,  intent(in)                          :: use_input_p
logical,  intent(in)                          :: bounded_below, bounded_above
real(r8), intent(in)                          :: lower_bound,   upper_bound

! Probit transform for gamma.
real(r8) :: quantile
integer  :: i

! Bounds other than a lower bound at 0 not yet implemented for gamma distribution

! Get the parameters for this distribution if not already available
if(.not. use_input_p) then 

   ! In full generality, gamma must be bounded either below or above
   if(.not. (bounded_below .neqv. bounded_above)) then
      errstring = 'Gamma distribution requires either bounded above or below to be true'
      call error_handler(E_ERR, 'to_probit_gamma', errstring, source)
   endif

   call set_gamma_params_from_ens(state_ens, ens_size, bounded_below, bounded_above, &
      lower_bound, upper_bound, p)
endif

do i = 1, ens_size
   ! First, get the quantile for this ensemble member
   quantile = gamma_cdf_params(state_ens(i), p)
   ! Transform to probit space 
   probit_ens(i) = probit_or_logit_transform(quantile)
end do

end subroutine to_probit_gamma

!------------------------------------------------------------------------

subroutine to_probit_beta(ens_size, state_ens, p, probit_ens, use_input_p, &
   lower_bound, upper_bound)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: state_ens(ens_size)
type(distribution_params_type), intent(inout) :: p
real(r8), intent(out)                :: probit_ens(ens_size)
logical, intent(in)                  :: use_input_p
real(r8), intent(in)                 :: lower_bound, upper_bound

! Probit transform for beta.
real(r8) :: quantile
integer  :: i

! Get the parameters for this distribution if not already available
if(.not. use_input_p) then
   call set_beta_params_from_ens(state_ens, ens_size, lower_bound, upper_bound, p)
endif

do i = 1, ens_size
   ! First, get the quantile for this ensemble member
   quantile = beta_cdf_params(state_ens(i), p)
   ! Transform to probit/logit space 
   probit_ens(i) = probit_or_logit_transform(quantile)
end do

end subroutine to_probit_beta

!------------------------------------------------------------------------

subroutine to_probit_bounded_normal_rh(ens_size, state_ens, p, probit_ens, &
   use_input_p, bounded_below, bounded_above, lower_bound, upper_bound)

! Note that this is just for transforming back and forth, not for doing the RHF observation update
! This means that we know a prior that the quantiles associated with the initial ensemble are
! uniformly spaced which can be used to simplify transforming.

! How to handle identical ensemble members is an open question for now. This is also a problem
! for ensemble members that are identical to one of the bounds. 

integer,  intent(in)                          :: ens_size
real(r8), intent(in)                          :: state_ens(ens_size)
type(distribution_params_type), intent(inout) :: p
real(r8), intent(out)                         :: probit_ens(ens_size)
logical,  intent(in)                          :: use_input_p
logical,  intent(in)                          :: bounded_below, bounded_above
real(r8), intent(in)                          :: lower_bound,   upper_bound

! Probit transform for bounded normal rh.
integer  :: i
real(r8) :: quantile(ens_size)

if(use_input_p) then
   ! Do not know what to do if sd of original ensemble is 0 (or small, work on this later)
   if(get_bnrh_sd(p) <= 0.0_r8) then
      ! Just return the original ensemble
      probit_ens = state_ens 
      return
   endif

   ! Get the quantiles for each of the ensemble members in a BNRH distribution
   call bnrh_cdf_initialized_vector(state_ens, ens_size, p, quantile)

else
   ! Get all the info about the rank histogram cdf
   call bnrh_cdf_params(state_ens, ens_size, bounded_below, bounded_above, &
      lower_bound, upper_bound, p, quantile)

   ! Do not know what to do if sd is 0 (or small, work on this later)
   if(get_bnrh_sd(p) <= 0.0_r8) then
      ! Just return the original ensemble
      probit_ens = state_ens 
      return
   endif

endif

! Transform the quantiles to probit space
do i = 1, ens_size
   probit_ens(i) = probit_or_logit_transform(quantile(i)) 
end do

end subroutine to_probit_bounded_normal_rh

!------------------------------------------------------------------------

!!!subroutine to_probit_particle(ens_size, state_ens, p, probit_ens, &
   !!!use_input_p, bounded_below_in, bounded_above_in, lower_bound_in, upper_bound_in)
!!!
!!!! Doing a particle filter. Quantiles are (2i-1) / 2n 
!!!
!!!integer, intent(in)                  :: ens_size
!!!real(r8), intent(in)                 :: state_ens(ens_size)
!!!type(distribution_params_type), intent(inout) :: p
!!!real(r8), intent(out)                :: probit_ens(ens_size)
!!!logical, intent(in)                  :: use_input_p
!!!logical, intent(in)                  :: bounded_below_in, bounded_above_in
!!!real(r8), intent(in)                 :: lower_bound_in,   upper_bound_in
!!!
!!!integer  :: i, j, indx
!!!integer  :: ens_index(ens_size)
!!!real(r8) :: quantile
!!!
!!!! This should fail if any of the input states are not the same as one of the 
!!!! original ensemble states when use_input_p is false. 
!!!if(use_input_p) then
   !!!! The particles are available from a previous call
   !!!! The input member gets the same quantile as the corresponding member from the previous call
   !!!! This can be done vastly more efficiently with either binary searches or by first sorting the
   !!!! incoming state_ens so that the lower bound for starting the search is updated with each ensemble member
   !!! 
   !!!do i = 1, ens_size
      !!!! Loop through the previous ensemble members
      !!!quantile = -99_r8
      !!!do j = 1, ens_size
         !!!! Is exact equivalence a problem here?
         !!!if(state_ens(i) == p%params(j)) then
            !!!quantile = 2*(j-1) / (2*ens_size)
            !!!exit
         !!!endif
         !!!! Test failed to find a match
         !!!if(quantile < 0.0_r8) then
            !!!write(errstring, *) 'Unable to find prior for use_input_p', state_ens(i)
            !!!call error_handler(E_ERR, 'to_probit_particle', errstring, source)
         !!!endif
         !!!! Do probit/logit transform
         !!!probit_ens(i) = probit_or_logit_transform(quantile)
      !!!end do
   !!!end do
  !!! 
!!!else
   !!!! Not using a pre-existing distribution
   !!!! Take care of space for the transform data structure, just need to know sorted prior members
   !!!if(allocated(p%params)) deallocate(p%params)
   !!!allocate(p%params(ens_size))
!!!
   !!!! For particle filter, the required data for inversion is the original ensemble values
   !!!! Having them in sorted order is useful for subsequent inversion
   !!!call index_sort(state_ens, ens_index, ens_size)
   !!!p%params(1:ens_size) = state_ens(ens_index)
!!!
   !!!! Get the quantiles for each of the ensemble members
   !!!do i = 1, ens_size
      !!!indx = ens_index(i)
      !!!! The quantiles for a particle filter are just 2(i-1) / 2n
      !!!quantile = 2*(indx - 1) / (2 * ens_size) 
!!!
      !!!! Transform the quantiles to probit/logit space
      !!!probit_ens(indx) = probit_or_logit_transform(quantile)
   !!!end do 
!!!
!!!endif
!!!
!!!end subroutine to_probit_particle
!!!
!------------------------------------------------------------------------

subroutine transform_all_from_probit(ens_size, num_vars, probit_ens, p, state_ens)

integer, intent(in)                  :: ens_size
integer, intent(in)                  :: num_vars
real(r8), intent(in)                 :: probit_ens(:, :)
type(distribution_params_type), intent(inout) :: p(num_vars)
real(r8), intent(out)                :: state_ens(:, :)

! Transform back to the original space
integer  :: i
real(r8) :: temp_ens(ens_size)

do i = 1, num_vars
   call transform_from_probit(ens_size, probit_ens(1:ens_size, i), p(i), temp_ens)
   state_ens(1:ens_size, i) = temp_ens
end do

end subroutine transform_all_from_probit

!------------------------------------------------------------------------

subroutine transform_from_probit(ens_size, probit_ens, p, state_ens)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: probit_ens(ens_size)
type(distribution_params_type), intent(inout) :: p
real(r8), intent(out)                :: state_ens(ens_size)

! If not initialized, read in the namelist
if(.not. module_initialized) call initialize_probit_transform

! Transform back to the original space
if(p%distribution_type == NORMAL_DISTRIBUTION) then
   ! No need to do any transformation for a normal
   state_ens = probit_ens
elseif(p%distribution_type == LOG_NORMAL_DISTRIBUTION) then
   call from_probit_log_normal(ens_size, probit_ens, state_ens)
elseif(p%distribution_type == UNIFORM_DISTRIBUTION) then
   call from_probit_uniform(ens_size, probit_ens, p, state_ens)
elseif(p%distribution_type == GAMMA_DISTRIBUTION) then
   call from_probit_gamma(ens_size, probit_ens, p, state_ens)
elseif(p%distribution_type == BETA_DISTRIBUTION) then
   call from_probit_beta(ens_size, probit_ens, p, state_ens)
elseif(p%distribution_type == BOUNDED_NORMAL_RH_DISTRIBUTION) then
   call from_probit_bounded_normal_rh(ens_size, probit_ens, p, state_ens)
!!!elseif(p%distribution_type == PARTICLE_FILTER_DISTRIBUTION) then
   !!!call from_probit_particle(ens_size, probit_ens, p, state_ens)
else
   write(errstring, *) 'Illegal distribution type', p%distribution_type
   call error_handler(E_ERR, 'transform_from_probit', errstring, source)
   stop
endif

! Deallocate any allocatable storage that was used for this distribution
call deallocate_distribution_params(p)

end subroutine transform_from_probit

!------------------------------------------------------------------------

subroutine from_probit_log_normal(ens_size, probit_ens, state_ens)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: probit_ens(ens_size)
real(r8), intent(out)                :: state_ens(ens_size)

! Take the inverse of the log to get back to original space
state_ens = exp(probit_ens)

end subroutine from_probit_log_normal

!------------------------------------------------------------------------

subroutine from_probit_uniform(ens_size, probit_ens, p, state_ens)

integer, intent(in)                           :: ens_size
real(r8), intent(in)                          :: probit_ens(ens_size)
type(distribution_params_type), intent(inout) :: p
real(r8), intent(out)                         :: state_ens(ens_size)

real(r8) :: quantile
integer :: i

do i = 1, ens_size
   ! First, invert the probit to get a quantile
   quantile = inv_probit_or_logit_transform(probit_ens(i))
   ! Transform from U(0, 1) to U(lower_bound, upper_bound)
   state_ens(i) = p%lower_bound + quantile * (p%upper_bound - p%lower_bound)
end do

end subroutine from_probit_uniform

!------------------------------------------------------------------------

subroutine from_probit_gamma(ens_size, probit_ens, p, state_ens)

integer,  intent(in)                          :: ens_size
real(r8), intent(in)                          :: probit_ens(ens_size)
type(distribution_params_type), intent(inout) :: p
real(r8), intent(out)                         :: state_ens(ens_size)

! Transform back to the original space
real(r8) :: quantile
integer  :: i

do i = 1, ens_size
   ! First, invert the probit/logit to get a quantile
   quantile = inv_probit_or_logit_transform(probit_ens(i))
   ! Invert the gamma quantiles to get physical space
   state_ens(i) = inv_gamma_cdf_params(quantile, p)
end do

end subroutine from_probit_gamma

!------------------------------------------------------------------------

subroutine from_probit_beta(ens_size, probit_ens, p, state_ens)

integer, intent(in)                  :: ens_size
real(r8), intent(in)                 :: probit_ens(ens_size)
type(distribution_params_type), intent(inout) :: p
real(r8), intent(out)                :: state_ens(ens_size)

! Transform back to the original space
real(r8) :: quantile
integer  :: i

do i = 1, ens_size
   ! First, invert the probit/logit to get a quantile
   quantile = inv_probit_or_logit_transform(probit_ens(i))
   ! Invert the beta quantiles to get scaled physical space
   state_ens(i) = inv_beta_cdf_params(quantile, p)
end do

end subroutine from_probit_beta

!------------------------------------------------------------------------

subroutine from_probit_bounded_normal_rh(ens_size, probit_ens, p, state_ens)

integer, intent(in)                          :: ens_size
real(r8), intent(in)                         :: probit_ens(ens_size)
type(distribution_params_type), intent(inout) :: p
real(r8), intent(out)                        :: state_ens(ens_size)

integer :: i
real(r8) :: quantiles(ens_size)

! Do not know what to do if original ensemble had all members the same (or nearly so???)
if(get_bnrh_sd(p) <= 0.0_r8) then
   state_ens = probit_ens
else

   ! Transform each probit ensemble member back to physical space
   do i = 1, ens_size
      ! First, invert the probit/logit to get quantiles
      quantiles(i) = inv_probit_or_logit_transform(probit_ens(i))
   end do
   
   ! Invert the rank histogram CDF to get the physical space ensemble
   call inv_bnrh_cdf_params(quantiles, ens_size, p, state_ens)
endif

end subroutine from_probit_bounded_normal_rh

!------------------------------------------------------------------------

!!!subroutine from_probit_particle(ens_size, probit_ens, p, state_ens)
!!!
!!!integer, intent(in)                  :: ens_size
!!!real(r8), intent(in)                 :: probit_ens(ens_size)
!!!type(distribution_params_type), intent(inout) :: p
!!!real(r8), intent(out)                :: state_ens(ens_size)
!!!
!!!integer :: i, indx
!!!real(r8) :: quantile
!!!
!!!do i = 1, ens_size
   !!!! First invert the probit/logit transform to tg
   !!!quantile = inv_probit_or_logit_transform(probit_ens(i))
!!!
   !!!! Invert the quantile for a particle prior
   !!!! There is a prior ensemble member associated with each 1/ens_size fraction of the quantile 
   !!!! range
   !!!indx = floor(quantile * ens_size) + 1
   !!!if(indx <= 0) indx = 1
   !!!state_ens(i) = p%more_params(indx)
!!!end do
!!!
!!!! Probably do this explicitly 
!!!! Free the storage
!!!deallocate(p%more_params)
!!!
!!!end subroutine from_probit_particle

!------------------------------------------------------------------------

function probit_or_logit_transform(quantile)

real(r8)             :: probit_or_logit_transform
real(r8), intent(in) :: quantile

! Transform the quantile 
if(use_logit_instead_of_probit) then
   probit_or_logit_transform =  log(quantile / (1.0_r8 - quantile))
else
   probit_or_logit_transform = inv_std_normal_cdf(quantile)
endif

end function probit_or_logit_transform

!------------------------------------------------------------------------

function inv_probit_or_logit_transform(p)

real(r8)             :: inv_probit_or_logit_transform
real(r8), intent(in) :: p 

! Transform back to get a quantile
if(use_logit_instead_of_probit) then
   inv_probit_or_logit_transform = 1.0_r8 / (1.0_r8 + exp(-p))
else
   inv_probit_or_logit_transform = normal_cdf(p, 0.0_r8, 1.0_r8)
endif

end function inv_probit_or_logit_transform

!------------------------------------------------------------------------
subroutine initialize_probit_transform()

integer :: iunit, io

module_initialized = .true.

! Read the namelist entry
call find_namelist_in_file("input.nml", "probit_transform_nml", iunit)
read(iunit, nml = probit_transform_nml, iostat = io)
call check_namelist_read(iunit, io, "probit_transform_nml")

if (do_nml_file()) write(nmlfileunit,nml=probit_transform_nml)
if (do_nml_term()) write(     *     ,nml=probit_transform_nml)

end subroutine initialize_probit_transform

!------------------------------------------------------------------------
function fix_bounds(x, bounded_below, bounded_above,  lower_bound, upper_bound)

real(r8)             :: fix_bounds
real(r8), intent(in) :: x
logical,  intent(in) :: bounded_below, bounded_above
real(r8), intent(in) :: lower_bound,   upper_bound

! A variety of round off errors can lead to small violations of the bounds for state and
! observation quantities. This function corrects the violations if they are small. If 
! they are bigger than the egregious bound set here, then execution is terminated.
      
real(r8), parameter :: egregious_bound_threshold = 1.0e-12_r8

! Default behavior is to leave x unchanged
fix_bounds = x

! Fail here on egregious violations; this could be removed 
if(bounded_below) then
   if(lower_bound - x > egregious_bound_threshold) then
      write(errstring, *) 'Egregious lower bound violation (see code)', x, lower_bound
      call error_handler(E_ERR, 'fix_bounds', errstring, source)
   else
      fix_bounds = max(x, lower_bound)
   endif
endif

if(bounded_above) then
   if(x - upper_bound > egregious_bound_threshold) then
      write(errstring, *) 'Egregious upper bound violoation first check(see code)', x, upper_bound
      call error_handler(E_ERR, 'fix_bounds', errstring, source)
   else
      fix_bounds = min(x, upper_bound)
   endif
endif

end function fix_bounds

!------------------------------------------------------------------------

end module probit_transform_mod
