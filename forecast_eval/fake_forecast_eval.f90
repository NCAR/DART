module forecast_eval_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
! Used to do forecast evaluation at leads less frequent than the assimilation
! frequency for the PE 9var paper. This should NOT be used in any applications
! where it is not a perfect model because it advances the truth directly.

! Computes ensemble forecasts and statistics for evaluating their quality

use assim_diag_mod, only : diag_type, output_diagnostics, save_diagnostics, assim_diag_init

use model_mod, only : advance

implicit none

private
public forecast_advance, forecast_out

! Maximum number of observation periods for which forecasts should extend
integer, parameter :: max_leads = 100, eval_freq = 10
integer :: last_time_step = -1
integer :: ind = -1, num_avail

! Diagnostic type to keep statistics for forecasts of each lead time
type (diag_type) :: diag_stuff(max_leads)

contains

!=====================================================================

subroutine forecast_advance(truth_in, ens_in, time, var_list, num_var, add_on)

implicit none

!-------------------------------------------------------------------------
! Given the truth at a particular time, the assimilated ensemble, and the
! time step, evaluates current forecast sets of available leads and starts
! the newest forecast set integration.
!
! truth(model_size):		True state of system
! ens(model_size, ens_size):	Ensemble of assimilated states
! time:				Integer, time step corresponding to truth
!-------------------------------------------------------------------------

double precision, intent(in) :: truth_in(:), ens_in(:, :)
integer, intent(in) :: time, num_var, var_list(num_var)
logical :: add_on

integer :: ens_size, model_size, num_steps, i, j
double precision :: ens_mean(size(ens_in, 1))
double precision :: truth(size(truth_in)), ens(size(ens_in, 1), size(ens_in, 2))

ens_size = size(ens, 2)
model_size = size(truth)

! Is this first call? If so initialize storage and pointers
if(last_time_step == -1) then
   last_time_step = 1
! Also need to initialize assim_diag structures (added 12 Feb. 2001)
   write(*, *) 'initializeing in forecast_advance ens_size ', ens_size
   do i = 1, max_leads
      call assim_diag_init(diag_stuff(i), model_size, ens_size, var_list, num_var, add_on)
   end do
endif

! All calls; advance each forecast and truth requisite number of steps
ens = ens_in
truth = truth_in
do i = 1, max_leads
   do j = 1, ens_size
      call advance(ens(:, j), eval_freq, ens(:, j))
   end do   
   call advance(truth, eval_freq, truth)

! Compute ensemble mean of current ensemble for validation
   ens_mean = sum(ens, dim=2) / ens_size

!      call save_diagnostics(ens(:, :), ens_mean, truth, diag_stuff(i), add_on)
   call save_diagnostics(ens(:, :), ens_mean, truth, diag_stuff(i))
end do

end subroutine forecast_advance

!======================================================================

subroutine forecast_out(chan, diag_output_index)

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

end module forecast_eval_mod
