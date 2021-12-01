! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Operations and storage required for the adaptive hybrid scheme

module hybrid_mod

!> \defgroup hybrid_mod
!> @{

use types_mod,            only : r8, PI, MISSING_R8
use time_manager_mod,     only : time_type, get_time
use utilities_mod,        only : open_file, close_file, error_handler, E_ERR, E_MSG
use random_seq_mod,       only : random_seq_type, random_gaussian, init_random_seq
use ensemble_manager_mod, only : ensemble_type, map_pe_to_task
use mpi_utilities_mod,    only : my_task_id, send_to, receive_from, send_minmax_to
use adaptive_inflate_mod, only : set_from_string, SINGLE_SS_INFLATION

implicit none
private

public :: NO_HYBRID, CONSTANT_HYBRID, SINGLE_SS_HYBRID, VARYING_SS_HYBRID,  & 
          hybrid_type, validate_hybrid_options, adaptive_hybrid_init,       &
          hyb_mean_from_restart, hyb_sd_from_restart, set_hybrid_mean_copy, &
          set_hybrid_sd_copy, do_single_ss_hybrid, get_hybrid_mean_copy,    &
          get_hybrid_sd_copy, do_ss_hybrid, get_hybrid_mean, get_hybrid_sd, &
          do_varying_ss_hybrid, hyb_minmax_task_zero, log_hybrid_info,      & 
          update_hybrid, get_hybrid_sample_size

character(len=*), parameter :: source = 'hybrid_mod.f90'

integer, parameter :: NO_HYBRID         = 0
integer, parameter :: CONSTANT_HYBRID   = 1
integer, parameter :: SINGLE_SS_HYBRID  = 2
integer, parameter :: VARYING_SS_HYBRID = 3

! Type to keep track of information for hybrid 
type hybrid_type
   private
   integer               :: flavor
   integer               :: sample_size
   logical               :: output_restart = .false.
   real(r8)              :: weight, sd
   real(r8)              :: minmax_mean(2), minmax_sd(2)
   logical               :: mean_from_restart
   logical               :: sd_from_restart
   integer               :: input_mean_copy = -1 !!todo NO_COPY_PRESENT
   integer               :: input_sd_copy   = -1
end type hybrid_type

! Module storage for writing error messages
character(len=512) :: string1

! Flag indicating whether module has been initialized
logical :: initialized = .false.

! Climatology sample size
integer, parameter     :: clim_size_accept = 1000

! Bounds for the weighting coefficients
real(r8) :: weight_upper_bound = 1.0_r8
real(r8) :: weight_lower_bound = 0.0_r8

!===============================================================================

contains


!-------------------------------------------------------------------------------
!> Make sure the combination of hybrid options are legal

subroutine validate_hybrid_options(inf_flavor, hyb_flavor, hyb_ens_size, & 
                  hyb_initial_from_restart, hyb_sd_initial_from_restart, & 
                  do_hybrid, output_hybrid)

integer, intent(in)    :: inf_flavor(2)
integer, intent(in)    :: hyb_flavor
integer, intent(in)    :: hyb_ens_size
logical, intent(inout) :: hyb_initial_from_restart
logical, intent(inout) :: hyb_sd_initial_from_restart
logical, intent(out)   :: do_hybrid
logical, intent(out)   :: output_hybrid

do_hybrid     = .false.
output_hybrid = .false.

if(hyb_flavor < NO_HYBRID .or. hyb_flavor > VARYING_SS_HYBRID) then
   write(string1, *) 'hyb_flavor=', hyb_flavor, ' Must be 0, 1, 2, 3 '

   call error_handler(E_ERR,'validate_hybrid_options', string1, source)
endif

! Check to see if state space hybrid is turned on
if (hyb_flavor /= NO_HYBRID) do_hybrid = .true.
   
! Adaptive hybrid scheme is only spatially-constant for now, 
! Inflation is relied on to spread the information properly in space
! Only support hybrid with inf-flavors 0, 2, 4, 5 (1 is deprecated)
if (do_hybrid .and. hyb_flavor >= SINGLE_SS_HYBRID) then 
   if (inf_flavor(1) == SINGLE_SS_INFLATION .or. inf_flavor(2) == SINGLE_SS_INFLATION) then 
      write(string1, '(A, i2, A)') 'Adaptive hybrid scheme does not support inf_flavor:', SINGLE_SS_HYBRID, ' (i.e., spatially-uniform inflation)' 
      call error_handler(E_ERR, 'validate_hybrid_options:', string1, source)
   endif
endif
 

! Check the size of the climatological ensemble
! We need this to be fairly large ~O(1000)
if (do_hybrid .and. hyb_ens_size < clim_size_accept) then 
   write(string1, '(A, i5, 2A, i5, A)') 'Climatology ensemble is of size ', hyb_ens_size, '. The recommended sample size ', &
                       'should be at least ', clim_size_accept, ' for statistical robustness.' 
   call error_handler(E_MSG, 'validate_hybrid_options:', string1, source)
endif

if (do_hybrid) output_hybrid = .true.

! Observation space inflation not currently supported
if(hyb_flavor == VARYING_SS_HYBRID) then 
   call error_handler(E_ERR, 'validate_hybrid_options', &
                      'Varying state-space hybridization (type 3) not currently supported', source, &
                      text2 = 'Look out for this type in future DART releases.')
endif

end subroutine validate_hybrid_options


!-------------------------------------------------------------------------------
!> Initializes a hybrid_type 

subroutine adaptive_hybrid_init(hybrid_handle, hyb_flavor, ens_size, mean_from_restart, sd_from_restart, & 
                                output_hybrid, hyb_weight_initial, hyb_weight_sd_initial)

type(hybrid_type), intent(inout) :: hybrid_handle
integer,           intent(in)    :: hyb_flavor
integer,           intent(in)    :: ens_size
logical,           intent(in)    :: mean_from_restart
logical,           intent(in)    :: sd_from_restart
logical,           intent(in)    :: output_hybrid
real(r8),          intent(in)    :: hyb_weight_initial, hyb_weight_sd_initial

! Record the module version if this is first initialize call
if(.not. initialized) then
   initialized = .true.
endif

! Load up the structure first to keep track of all details of this inflation type
hybrid_handle%flavor            = hyb_flavor
hybrid_handle%sample_size       = ens_size
hybrid_handle%output_restart    = output_hybrid
hybrid_handle%weight            = hyb_weight_initial
hybrid_handle%sd                = hyb_weight_sd_initial
hybrid_handle%mean_from_restart = mean_from_restart
hybrid_handle%sd_from_restart   = sd_from_restart

! give these distinctive values; if inflation is being used
hybrid_handle%minmax_mean(:) = MISSING_R8
hybrid_handle%minmax_sd(:)   = MISSING_R8

end subroutine adaptive_hybrid_init

!-------------------------------------------------------------------------------
!> Accessor functions for adaptive hybrid type

function hyb_mean_from_restart(hybrid_handle)

type(hybrid_type) :: hybrid_handle
logical           :: hyb_mean_from_restart

hyb_mean_from_restart = hybrid_handle%mean_from_restart

end function hyb_mean_from_restart

!-------------------------------------------------------------------------------
!>

function hyb_sd_from_restart(hybrid_handle)

type(hybrid_type) :: hybrid_handle
logical           :: hyb_sd_from_restart

hyb_sd_from_restart = hybrid_handle%sd_from_restart

end function hyb_sd_from_restart

!-------------------------------------------------------------------------------
!>

function do_ss_hybrid(hybrid_handle)

type(hybrid_type), intent(in) :: hybrid_handle
logical                       :: do_ss_hybrid

if (do_single_ss_hybrid(hybrid_handle) .or. do_varying_ss_hybrid(hybrid_handle) ) then
   do_ss_hybrid = .true.
else
   do_ss_hybrid = .false.
endif

end function do_ss_hybrid

!-------------------------------------------------------------------------------
!> Returns true if varying state space hybrid scheme

function do_varying_ss_hybrid(hybrid_handle)

logical                       :: do_varying_ss_hybrid
type(hybrid_type), intent(in) :: hybrid_handle

do_varying_ss_hybrid = (hybrid_handle%flavor == VARYING_SS_HYBRID)

end function do_varying_ss_hybrid


!-------------------------------------------------------------------------------
!> Returns true if fixed state space hybrid scheme

function do_single_ss_hybrid(hybrid_handle)

logical                       :: do_single_ss_hybrid
type(hybrid_type), intent(in) :: hybrid_handle

do_single_ss_hybrid = (hybrid_handle%flavor == SINGLE_SS_HYBRID)

end function do_single_ss_hybrid


!-------------------------------------------------------------------------------
!> Given prior mean and sd values for the hybrid weighting coefficient, 
!> this routine finds the updated (posterior) values of both the mean and sd.
!> The algorithm assumes the weighting coefficient to follow a Gaussian 
!> distribution. The in-coming data (also Gaussian) is used to update the 
!> prior distribution. More details can be found here: 
!> 
!> El Gharamti, M. (2021). Hybrid Ensemble-Variational Filter: A Spatially and 
!> Temporally Varying Adaptive Algorithm to Estimate Relative Weighting. 
!> Monthly Weather Review, 149(1), 65-76.  
!>
!> May need to explore using a beta distribution for the weighting factor
!> instead of a normal one. Beta is nicely bounded between 0 and 1. 
!>
!> "r" parameter gives the localized prior correlation times the localization
!> which is computed in the assim_tools routine filter_assim. For single state
!> space hybrid it is 1.0.

subroutine update_hybrid(weight, weight_sd, se2, ss2, so2, x, obs, rho)

real(r8),          intent(inout) :: weight, weight_sd
real(r8),          intent(in)    :: x, se2, ss2 
real(r8),          intent(in)    :: obs, so2 
real(r8),          intent(in)    :: rho 

real(r8), parameter :: small_diff = 1.0e-8

integer  :: k, find_index(1)
real(r8) :: m, v, d2, ss, Y, Z, Y2, Z2
real(r8) :: a, b, c, a0, a1, a2, a3
real(r8) :: Q, R, G, H, theta 
real(r8) :: disc, sol(3), abssep(3)
real(r8) :: mulfac, addfac, PIfac
real(r8) :: new_weight, new_weight_sd

! If the weight_sd not positive, keep everything the same
if(weight_sd <= 0.0_r8) return

m  = weight
v  = weight_sd**2
d2 = (obs - x)**2

ss = se2 - ss2 

!print *, 'd2: ', d2
!print *, 'ss: ', ss

if (ss == 0.0_r8 .or. abs(ss) <= small_diff) then
   ! If the ensemble and static variances are very close 
   ! Then just keep the current weighting. 

   write(string1, *) 'Difference between ensmeble and climatology variance for this obs ', &
                     'is negligible. Skipping the update of the weighting factor.'
   call error_handler(E_MSG, 'update hybrid:', string1, source)

   return !can't do anything because we need to divide by that no. 
endif

! Simplify coefficients
Y  = so2+ss2
Z  = rho*ss
Y2 = Y**2
Z2 = Z**2

! Polynormial coefficients:
! a0x^3 + a1x^2 + a2x + a3 = 0
a0 = -2*Z2
a1 =  2*m*Z2  - 4*Y*Z
a2 =  4*m*Y*Z - 2*Y2        - rho*v*Z*ss
a3 =  2*m*Y2  + d2*rho*v*ss - rho*v*Y*ss

! Rearrage the polynomial such that a == 1
! x^3 + ax^2 + bx + c = 0
a = a1 / a0
b = a2 / a0
c = a3 / a0

! Find discriminant: 
! Q and R are real because a, b, c are real
Q = ( a**2 - 3.0_r8 * b ) / 9.0_r8
R = ( 2.0_r8 * a**3 - 9.0_r8 * a * b + 27.0_r8 * c ) / 54.0_r8

disc = R**2 - Q**3

!print *, 'disc: ', disc
!print *, 'Q: ', Q, ' R: ', R

! Find cubic roots
if (disc < 0.0_r8) then
   ! 3 distict real roots
 
   theta = acos( R / sqrt(Q**3) ) / 3.0_r8  
   !print *, 'theta: ', theta

   mulfac = - 2.0_r8 * sqrt(Q)
   addfac = - a / 3.0_r8
   PIfac  = 2.0_r8 * PI / 3.0_r8
   
   sol(1) = mulfac * cos(theta        ) + addfac
   sol(2) = mulfac * cos(theta + PIfac) + addfac
   sol(3) = mulfac * cos(theta - PIfac) + addfac 
  
   abssep = abs(sol - m)

   !print *, 'sol: ', sol

   ! Closest root
   find_index = minloc(abssep)
   new_weight = sol(find_index(1))

else
   ! only one real

   G = - abs(R) / R * (abs(R) + sqrt(disc))**(1.0/3.0)
   if (G == 0.0_r8) then 
      H = 0.0_r8
   else
      H = Q / G
   endif
   
   !print *, 'G: ', G, 'H: ', H
   
   new_weight = G + H - a / 3.0_r8

endif

!print *, 'new_weight = ', new_weight, ' d2 = ', d2, ' se2 = ', se2, ' ss2 = ', ss2, ' so2 = ', so2 

! Make sure weight satisfies constraints
weight = new_weight
if(weight < weight_lower_bound) weight = weight_lower_bound
if(weight > weight_upper_bound) weight = weight_upper_bound

! For now, weight sd is assumed fixed in time. 

end subroutine update_hybrid

!-------------------------------------------------------------------------------
!> Write to log file what kind of hybridization is being used.  

subroutine log_hybrid_info(hybrid_handle, mype, single_file)

type(hybrid_type), intent(in) :: hybrid_handle
integer,           intent(in) :: mype
logical,           intent(in) :: single_file

character(len = 128) :: det, tadapt, sadapt, akind, from

! nothing to do if not task 0
if (mype /= 0) return

! if inflation is off, say so and return now
if (hybrid_handle%flavor <= NO_HYBRID) then
   call error_handler(E_MSG, 'Hybrid scheme:', 'None', source)
   return
endif

write(string1, '(A, i5)') 'Sample size is ', hybrid_handle%sample_size
call error_handler(E_MSG, 'Hybrid scheme:', string1, source)

select case(hybrid_handle%flavor)
   case (CONSTANT_HYBRID)
      tadapt = 'time-invariant,'
      sadapt = ' spatially-constant'
   case (SINGLE_SS_HYBRID)
      tadapt = 'time-varying,'
      sadapt = ' spatially-constant '
   case (VARYING_SS_HYBRID)
      tadapt = 'time-varying,'
      sadapt = ' spatially-varying'
   case default
      write(string1, *) 'Illegal hybrid flavor '
      call error_handler(E_ERR, 'log_hybrid_info', string1, source)
end select

write(string1, '(3A)') trim(tadapt), trim(sadapt)
call error_handler(E_MSG, 'Hybrid scheme:', string1, source)

! combination file, individual file or namelist
call set_from_string(hybrid_handle%mean_from_restart, single_file, from)

if (nval_to_log(hybrid_handle, from) == 1) then
   write(string1,  '(A, F8.3)') &
         'hybrid weight mean   '//trim(from)//', value: ', hybrid_handle%minmax_mean(1)
else
   write(string1,  '(A, 2F8.3)') &
         'hybrid weight mean   '//trim(from)//', min/max values: ', hybrid_handle%minmax_mean
endif
call error_handler(E_MSG, 'Hybrid scheme:', string1,  source)

call set_from_string(hybrid_handle%sd_from_restart, single_file, from)

if (nval_to_log(hybrid_handle, from) == 1) then
   write(string1,  '(A, F8.3)') &
         'hybrid weight stddev '//trim(from)//', value: ', hybrid_handle%minmax_sd(1)
else
   write(string1,  '(A, 2F8.3)') &
         'hybrid weight stddev '//trim(from)//', min/max values: ', hybrid_handle%minmax_sd
endif
call error_handler(E_MSG, 'Hybrid scheme:', string1,  source)

end subroutine log_hybrid_info


!-------------------------------------------------------------------------------
!>

function nval_to_log(hybrid_handle, from_string)

type(hybrid_type), intent(in) :: hybrid_handle
character(len=*),  intent(in) :: from_string
integer :: nval_to_log

if ((hybrid_handle%flavor == SINGLE_SS_HYBRID) .or. &
    (hybrid_handle%flavor == VARYING_SS_HYBRID .and. &
     from_string == 'from namelist')) then
   nval_to_log = 1
else
   nval_to_log = 2
endif

end function nval_to_log


!-------------------------------------------------------------------------------
!> Collect the min and max of hybrid weights on task 0
!> this block handles communicating the min/max local values to PE 0
!> if running with MPI, or just sets the min/max directly if reading
!> from a namelist.

subroutine hyb_minmax_task_zero(hybrid_handle, ens_handle, &
                    ss_weight_index, ss_weight_sd_index)

type(hybrid_type),   intent(inout) :: hybrid_handle
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ss_weight_index
integer,             intent(in)    :: ss_weight_sd_index

real(r8) :: minmax_mean(2), minmax_sd(2), global_val(2)

! if not using inflation, return now
if (hybrid_handle%flavor <= NO_HYBRID) return

if (hybrid_handle%mean_from_restart) then

   ! find min and max on each processor
   minmax_mean(1) = minval(ens_handle%copies(ss_weight_index, :))
   minmax_mean(2) = maxval(ens_handle%copies(ss_weight_index, :))

   ! collect on pe 0
   call send_minmax_to(minmax_mean, map_pe_to_task(ens_handle, 0), global_val)
   if (ens_handle%my_pe == 0) hybrid_handle%minmax_mean = global_val

else
 
   hybrid_handle%minmax_mean = hybrid_handle%weight

endif

if (hybrid_handle%sd_from_restart) then

   ! find min and max on each processor
   minmax_sd(1) = minval(ens_handle%copies(ss_weight_sd_index, :))
   minmax_sd(2) = maxval(ens_handle%copies(ss_weight_sd_index, :))

   ! collect on pe 0
   call send_minmax_to(minmax_sd, map_pe_to_task(ens_handle, 0), global_val)
   if (ens_handle%my_pe == 0) hybrid_handle%minmax_sd = global_val
else

   hybrid_handle%minmax_sd = hybrid_handle%sd 

endif

end subroutine hyb_minmax_task_zero


!-------------------------------------------------------------------------------
!>

subroutine set_hybrid_mean_copy(hybrid_handle, c)
type(hybrid_type), intent(inout) :: hybrid_handle
integer,           intent(in)    :: c

hybrid_handle%input_mean_copy = c

end subroutine set_hybrid_mean_copy


!-------------------------------------------------------------------------------
!>

subroutine set_hybrid_sd_copy(hybrid_handle, c)
type(hybrid_type), intent(inout) :: hybrid_handle
integer,           intent(in)    :: c

hybrid_handle%input_sd_copy = c

end subroutine set_hybrid_sd_copy


!-------------------------------------------------------------------------------
!>

function get_hybrid_mean_copy(hybrid_handle) result (c)
type(hybrid_type), intent(in) :: hybrid_handle
integer :: c

c = hybrid_handle%input_mean_copy

end function get_hybrid_mean_copy


!-------------------------------------------------------------------------------
!>

function get_hybrid_sd_copy(hybrid_handle) result (c)
type(hybrid_type), intent(in) :: hybrid_handle
integer :: c

c = hybrid_handle%input_sd_copy

end function get_hybrid_sd_copy

!-------------------------------------------------------------------------------
!>

function get_hybrid_mean(hybrid_handle)

type(hybrid_type) :: hybrid_handle
real(r8)          :: get_hybrid_mean

get_hybrid_mean = hybrid_handle%weight

end function


!-------------------------------------------------------------------------------
!>

function get_hybrid_sd(hybrid_handle)

type(hybrid_type) :: hybrid_handle
real(r8)          :: get_hybrid_sd

get_hybrid_sd = hybrid_handle%sd

end function

!-------------------------------------------------------------------------------
!>

function get_hybrid_sample_size(hybrid_handle) result (c)
type(hybrid_type), intent(in) :: hybrid_handle
integer :: c

c = hybrid_handle%sample_size

end function get_hybrid_sample_size

!===============================================================================
! end module hybrid_mod
!===============================================================================

!> @}

end module hybrid_mod
