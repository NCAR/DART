module forecast_eval_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
! Computes ensemble forecasts and statistics for evaluating their quality

use types_mod
use assim_diag_mod, only : diag_type, output_diagnostics, save_diagnostics, &
                           assim_diag_init
use model_mod,      only : advance

implicit none

private
public forecast_advance, forecast_out

! Maximum number of observation periods for which forecasts should extend

integer, parameter :: max_leads = 20
integer :: last_time_step = -1
integer :: ind = -1, num_avail

! Storage to keep forecasts, rotates
real(r8), allocatable :: f(:, :, :) 

! Diagnostic type to keep statistics for forecasts of each lead time
type (diag_type) :: diag_stuff(max_leads)

contains



  subroutine forecast_advance(truth, ens, time, var_list, num_var, add_on)
!=====================================================================
! subroutine forecast_advance(truth, ens, time, var_list, num_var, add_on)
! 
! Given the truth at a particular time, the assimilated ensemble, and the
! time step, evaluates current forecast sets of available leads and starts
! the newest forecast set integration.
!
! truth(model_size):		True state of system
! ens(model_size, ens_size):	Ensemble of assimilated states
! time:				Integer, time step corresponding to truth
!-------------------------------------------------------------------------

implicit none

real(r8), intent(in) :: truth(:), ens(:, :)
integer,  intent(in) :: time, num_var, var_list(num_var)
logical :: add_on

integer  :: ens_size, model_size, num_steps, i, j
real(r8) :: ens_mean(size(ens, 1))

ens_size   = size(ens, 2)
model_size = size(truth)

! Is this first call? If so initialize storage and pointers

if( last_time_step == -1 ) then

   last_time_step = time
   allocate(f(model_size, ens_size, max_leads))
   f(:, :, 1) = ens

   ! Set pointer into rotating storage for next new ensemble

   ind = 2
   num_avail = 1

   ! Also need to initialize assim_diag structures (added 12 Feb. 2001)

   write(*, *) 'initializing forecast_advance ens_size ', ens_size
   do i = 1, max_leads
      call assim_diag_init(diag_stuff(i), model_size, ens_size, var_list, num_var, add_on)
   end do
   return

endif

! Not the first call; advance each current forecast requisite number of steps

num_steps = time - last_time_step
do i = 1, num_avail
   do j = 1, ens_size
      call advance(f(:, j, i), num_steps, f(:, j, i))
   end do   
end do

! Do evaluation if all leads are available;
! Compute ensemble mean of current ensemble for validation

ens_mean = sum(ens, dim=2) / real(ens_size)

if( num_avail == max_leads ) then
   do i = 1, max_leads

      j = ind - i        ! Compute index into rotating forecast storage
      if(j < 1) j = j + max_leads

      ! Instead of obs, use mean of current analyzed state as control

      ens_mean = sum(f(:, :, j), dim = 2) / real(ens_size)
!     call save_diagnostics(f(:, :, j), ens_mean, truth, diag_stuff(i), add_on)
      call save_diagnostics(f(:, :, j), ens_mean, truth, diag_stuff(i))
   end do
endif

! Now load in the new forecast at the ind slot

f(:, :, ind) = ens
ind          = ind + 1
if(ind == max_leads + 1) ind = 1
if(num_avail < max_leads) num_avail = num_avail + 1

last_time_step = time          ! Update last time

end subroutine forecast_advance



  subroutine forecast_out(chan, diag_output_index)
!======================================================================
! subroutine forecast_out(chan, diag_output_index)
!
! Outputs statistics for each forecast lead time

implicit none

integer, intent(in) :: chan
integer, intent(in) :: diag_output_index(:)

integer :: i

do i = 1, max_leads
   write(*, *) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   write(*, *) 'output for forecast lead ', i
   write(chan, *) 'output for forecast lead ', i
!   call output_diagnostics(diag_stuff(i), chan, diag_output_index)
   call output_diagnostics(diag_stuff(i), chan)
end do

end subroutine forecast_out

!=======================================================================
! End of forecast_eval_mod.f90
!=======================================================================

end module forecast_eval_mod
