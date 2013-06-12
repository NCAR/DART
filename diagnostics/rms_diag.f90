! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program rms_diag

! Sample diagnostic program that computes rms error as function 
! of time between a truth found in one file, priors found in 
! another, and posterior found in a third.

use       types_mod, only : r8, metadatalength
use    location_mod, only : location_type
use assim_model_mod, only : init_diag_input, get_diag_input_copy_meta_data, &
   get_model_state_vector, input_diagnostics, assim_model_type, &
   static_init_assim_model, init_assim_model

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer :: prior_unit, posterior_unit, truth_unit, ens_size, model_size, copies
type(location_type), allocatable :: location(:)
character(len = metadatalength), allocatable :: meta_data(:), global_meta_data
character(len = 129) :: prior_file, posterior_file, truth_file
integer :: i, j, num_samples
real(r8), allocatable :: prior(:, :), posterior(:, :), truth(:)
type(assim_model_type) :: truth_state, prior_state, posterior_state
real(r8), allocatable :: prior_ens_mean(:), posterior_ens_mean(:)
real(r8) :: prior_rms_err, posterior_rms_err, sum_prior_rms = 0, sum_posterior_rms = 0
real(r8) :: sum_prior_mean_rms, sum_posterior_mean_rms, prior_mean_rms, posterior_mean_rms
real(r8) :: desired


! OBVIOUSLY, a mature facility would have more error checks and make more
! direct use of available meta data. This is just a demo.

! Initialize the assim_model
call static_init_assim_model()
call init_assim_model(truth_state)
call init_assim_model(prior_state)
call init_assim_model(posterior_state)

write(*, *) 'input file name for prior state ensembles'
read(*, *) prior_file
prior_unit = init_diag_input(prior_file, global_meta_data, model_size, ens_size)
allocate(meta_data(ens_size), location(model_size), truth(model_size), &
   prior(model_size, ens_size), posterior(model_size, ens_size), &
   prior_ens_mean(model_size), posterior_ens_mean(model_size))
call get_diag_input_copy_meta_data(prior_unit, model_size, ens_size, location, meta_data)

write(*, *) 'model size and ensemble size ', model_size, ens_size

write(*, *) 'input file name for posterior state ensembles'
read(*, *) posterior_file
posterior_unit = init_diag_input(posterior_file, global_meta_data, model_size, ens_size)
call get_diag_input_copy_meta_data(posterior_unit, model_size, ens_size, location, meta_data)

write(*, *) 'input file name for truth state'
read(*, *) truth_file
truth_unit = init_diag_input(truth_file, global_meta_data, model_size, copies)
call get_diag_input_copy_meta_data(truth_unit, model_size, copies, location, &
   meta_data(1:copies))

! LAME, need to fix this up so end of file is detected
write(*, *) 'Input number of samples in files '
read(*, *) num_samples


! Loop through to read each entry, assume that times are all synched up for now
! CHANGE TO A CHECKING WHILE LOOP
do i = 1, num_samples
   do j = 1, ens_size
      call input_diagnostics(prior_unit, prior_state, copies)
      prior(:, j) = get_model_state_vector(prior_state)
      call input_diagnostics(posterior_unit, posterior_state, copies)
      posterior(:, j) = get_model_state_vector(posterior_state)
   end do
   call input_diagnostics(truth_unit, truth_state, copies)
   truth = get_model_state_vector(truth_state)

   ! Do some very basic diagnotics 
   prior_ens_mean = sum(prior, dim=2) / ens_size
   posterior_ens_mean = sum(posterior, dim=2) / ens_size
   prior_rms_err = sqrt(sum((prior_ens_mean - truth) * (prior_ens_mean - truth)))
   posterior_rms_err = sqrt(sum((posterior_ens_mean - truth) * &
      (posterior_ens_mean - truth)))
   sum_prior_rms = sum_prior_rms + prior_rms_err
   sum_posterior_rms = sum_posterior_rms + posterior_rms_err
   write(*, *) i, prior_rms_err, posterior_rms_err

   prior_mean_rms = 0.0
   posterior_mean_rms = 0.0
   do j = 1, ens_size
      prior_mean_rms = prior_mean_rms + sqrt(sum((prior(:, j) - truth) * (prior(:, j) - truth)))
      posterior_mean_rms = posterior_mean_rms +  sqrt(sum((posterior(:, j) - truth) * &
         (posterior(:, j) - truth)))
   end do
   sum_prior_mean_rms = sum_prior_mean_rms + prior_mean_rms / ens_size
   sum_posterior_mean_rms = sum_posterior_mean_rms + posterior_mean_rms / ens_size
end do

! Verify the ratio 
desired = sqrt((ens_size + 1) / (2. * ens_size))

write(*, *) 'mean prior rms is     ', sum_prior_rms / num_samples
write(*, *) 'prior mean rms is     ', sum_prior_mean_rms / num_samples
write(*, *) 'normalize ratio is ', (sum_prior_rms / sum_prior_mean_rms) / desired

write(*, *) '-----------------'

write(*, *) 'mean posterior rms is ', sum_posterior_rms / num_samples
write(*, *) 'posterior mean rms is ', sum_posterior_mean_rms / num_samples
write(*, *) 'normalized ratio is ', (sum_posterior_rms / sum_posterior_mean_rms) / desired

end program rms_diag

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
