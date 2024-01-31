module distribution_params_mod

! Provides data structure and tools to represent probability distribution families for DART

use types_mod, only : r8

implicit none
private

type distribution_params_type
   integer               :: distribution_type
   logical               :: bounded_below, bounded_above
   real(r8)              :: lower_bound,   upper_bound
   real(r8)              :: params(2)
   integer               :: ens_size
   real(r8), allocatable :: ens(:)
   real(r8), allocatable :: more_params(:)
end type

! Defining parameter strings for different prior distributions that can be used for probit transform
integer, parameter :: NORMAL_DISTRIBUTION            = 1
integer, parameter :: BOUNDED_NORMAL_RH_DISTRIBUTION = 2
integer, parameter :: GAMMA_DISTRIBUTION             = 3
integer, parameter :: BETA_DISTRIBUTION              = 4
integer, parameter :: LOG_NORMAL_DISTRIBUTION        = 5
integer, parameter :: UNIFORM_DISTRIBUTION           = 6
integer, parameter :: PARTICLE_FILTER_DISTRIBUTION   = 7

public :: distribution_params_type, deallocate_distribution_params, &
   NORMAL_DISTRIBUTION, BOUNDED_NORMAL_RH_DISTRIBUTION, GAMMA_DISTRIBUTION, BETA_DISTRIBUTION, &
      LOG_NORMAL_DISTRIBUTION, UNIFORM_DISTRIBUTION, PARTICLE_FILTER_DISTRIBUTION

contains

!----------------------------------------------------------------------

subroutine deallocate_distribution_params(p)

type(distribution_params_type), intent(inout) :: p

! Free up the allocatable storage
if(allocated(p%ens)) deallocate(p%ens)
if(allocated(p%more_params)) deallocate(p%more_params)

end subroutine deallocate_distribution_params

!----------------------------------------------------------------------

end module distribution_params_mod
