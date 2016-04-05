! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program leeuwen_particle_filter

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

! Created 21 August, 2002. Based directly on resampling filter in draft paper
! by Leeuwen but consistent with methods in Anderson and Anderson, standard
! statistical literature, and in earlier repository versions.


use types_mod
use obs_sequence_mod, only : obs_sequence_type, write_obs_sequence, &
   read_obs_sequence, get_num_obs_sets, get_obs_sequence_time, &
   get_num_obs_in_set, get_expected_obs, get_diag_obs_err_cov, &
   get_obs_values, &
   obs_sequence_def_copy, inc_num_obs_copies, set_obs_values, &
   set_single_obs_value, get_obs_def_index
use time_manager_mod, only : time_type, set_time, print_time, operator(/=)
use utilities_mod,    only :  get_unit, open_file, close_file, check_nml_error, &
   file_exist
use assim_model_mod,  only : assim_model_type, static_init_assim_model, &
   get_model_size, get_initial_condition, get_closest_state_time_to, &
   advance_state, set_model_time, get_model_time, init_diag_output, &
   output_diagnostics, init_assim_model, get_state_vector_ptr, &
   write_state_restart, read_state_restart

use random_seq_mod,   only : random_seq_type, init_random_seq, random_gaussian, &
   random_uniform

implicit none

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

! Define a type for doing direct access to ensemble state vectors
type model_state_ptr_type
   real(r8), pointer :: state(:)
end type model_state_ptr_type

type(obs_sequence_type) :: seq, prior_seq, posterior_seq
type(time_type)         :: time, time2
type(random_seq_type)   :: random_seq


integer :: i, ii, j, jj, k, ind, unit, prior_obs_unit, posterior_obs_unit, io
integer :: prior_state_unit, posterior_state_unit, num_obs_in_set, ierr

! Need to set up namelists for controlling all of this mess, too!
integer, parameter :: ens_size = 50

! ADDED NOISE NEEDED FOR PARTICLE FILTER WITH DETERMINISTIC MODELS
real(r8), parameter :: added_noise_sd = 0.01
integer            :: model_size, num_obs_sets

! Storage for direct access to ensemble state vectors
type(model_state_ptr_type) :: ens_ptr(ens_size), x_ptr

type(assim_model_type) :: x, ens(ens_size)
real(r8)               :: obs_inc(ens_size), ens_inc(ens_size)
real(r8)               :: swath(ens_size), ens_mean
real(r8), allocatable  :: obs_err_cov(:), obs(:), ens_obs(:, :)
real(r8)               :: cov_factor
character(len = 129)   :: ens_copy_meta_data(ens_size)

!----------------------------------------------------------------
! Namelist input with default values
!
real(r8) :: cutoff = 1.0, cov_inflate = 1.0_r8
integer :: cache_size = 10
logical :: start_from_restart = .false., output_restart = .false.
! ens_size needs to be added as parameter at some point, but requires
! massive changes to allocate storage
!integer :: ens_size = 20

character(len = 129) :: obs_sequence_file_name = "obs_sequence", &
                        restart_in_file_name = 'filter_restart_in', &
                        restart_out_file_name = 'filter_restart_out'

namelist /filter_nml/ cutoff, cov_inflate, cache_size,  &
   start_from_restart, output_restart, &
   obs_sequence_file_name, restart_in_file_name, restart_out_file_name
!----------------------------------------------------------------


! Begin by reading the namelist input
if(file_exist('input.nml')) then
   unit = open_file('input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(unit, nml = filter_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'filter_nml')
   enddo
 11 continue
   call close_file(unit)
endif

! Input the obs_sequence
unit = get_unit()
write(*, *) 'obs_sequence_file_name is ', obs_sequence_file_name
open(unit = unit, file = obs_sequence_file_name)
seq = read_obs_sequence(unit)
close(unit)
! Count of number of sets in the sequence
num_obs_sets = get_num_obs_sets(seq)


! Copy just the definitions part of the sequence to the two output obs sequences
call obs_sequence_def_copy(prior_seq, seq)
call obs_sequence_def_copy(posterior_seq, seq)
! Set up the metadata for the output ensemble observations
do i = 1, ens_size
   write(ens_copy_meta_data(i), *) 'ensemble ', i
end do

! For now output all ensemble members for prior and posterior; add space
call inc_num_obs_copies(prior_seq, ens_size, ens_copy_meta_data)
call inc_num_obs_copies(posterior_seq, ens_size, ens_copy_meta_data)

! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()
model_size = get_model_size()

! Set up diagnostic output for model state

prior_state_unit = init_diag_output('prior_diag', &
   'prior ensemble state', ens_size, ens_copy_meta_data)
posterior_state_unit = init_diag_output('posterior_diag', &
   'posterior ensemble state', ens_size, ens_copy_meta_data)

!------------------- Read restart if requested ----------------------
if(start_from_restart) then
   unit = get_unit()
   open(unit = unit, file = restart_in_file_name)
   do i = 1, ens_size
      write(*, *) 'trying to read restart ', i
      call read_state_restart(ens(i), unit)
   end do
   close(unit)
   write(*, *) 'successfully read state restart'
!-----------------  Restart read in --------------------------------

else

!-----  Block to do cold start initialization of ensembles ----------
! Initialize the control and ensemble states and set up direct pointers
   call init_assim_model(x)
   x_ptr%state => get_state_vector_ptr(x)
   do i = 1, ens_size
      call init_assim_model(ens(i))
      ens_ptr(i)%state => get_state_vector_ptr(ens(i))
   end do

! Get the initial condition
   call get_initial_condition(x)

! Advance for a long time (5 days) to get things started?
! This should all be parameterized and controlled
!   time = set_time(0, 5)
!   write(*, *) 'calling advance state for x 5 days'
!   call advance_state(x, time)
!   write(*, *) 'back from advance state for x 5 days'


! Initialize a repeatable random sequence for perturbations
! Where should the magnitude of the perturbations come from here???
   call init_random_seq(random_seq)
! Perturb for ensembles;
   do i = 1, ens_size
      do j = 1, model_size
         ens_ptr(i)%state(j) = random_gaussian(random_seq, x_ptr%state(j), 1.0_r8)
      end do
   end do
!-------------------- End of cold start ensemble initialization block ------
endif

! Advance the model and ensemble to the closest time to the next
! available observations (need to think hard about these model time interfaces).

AdvanceTime : do i = 1, num_obs_sets

   time = get_obs_sequence_time(seq, i)
   write(*, *) 'time of obs set ', i
   call print_time(time)

! For now, set initial time of ensembles to time of first obs set
! Want more flexibility than this through namelist at some point
   if(i == 1) then
      do j = 1, ens_size
         call set_model_time(ens(j), time)
      end do
   endif

   time2 = get_closest_state_time_to(ens(1), time)
   write(*, *) 'advancing to time2 '
   call  print_time(time2)
   do j = 1, ens_size      ! Advance the ensembles to this time
      write(*, *) 'advancing ensemble member ', j
! Advancing to same time causes problem with B-grid diag calls
      if(time2 /= get_model_time(ens(j))) call advance_state(ens(j), time2)

! FOR PARTICLE FILTER, NEED NON-DETERMINISTIC MODEL; JUST ADD GAUSIAN NOISE FOR NOW
      do ii = 1, ens_size
         do jj = 1, model_size
           ens_ptr(ii)%state(jj) = &
              random_gaussian(random_seq, ens_ptr(ii)%state(jj), added_noise_sd)
         end do
      end do
      
      call output_diagnostics(prior_state_unit, ens(j), j)
   end do

   ! Do a covariance inflation for now?
   ! Inflate the ensemble state estimates
   do k = 1, model_size
      ens_mean = get_ens_mean(ens_ptr, ens_size, k)
      do j = 1, ens_size
         ens_ptr(j)%state(k) = ens_mean + (ens_ptr(j)%state(k) - &
            ens_mean) * sqrt(cov_inflate)
      end do
   end do

   ! How many observations in this set
   num_obs_in_set = get_num_obs_in_set(seq, i)
      write(*, *) 'num_obs_in_set is ', num_obs_in_set

   ! Allocate storage for the ensemble priors for this number of observations
   allocate(obs_err_cov(num_obs_in_set), obs(num_obs_in_set), &
      ens_obs(num_obs_in_set, ens_size))

   ! Get the observational error covariance (diagonal at present)
   call get_diag_obs_err_cov(seq, i, obs_err_cov)

   ! Get the observations; from copy 1 for now
   call get_obs_values(seq, i, obs, 1)

! Get the expected value of the obs for each ensemble member
   do k = 1, ens_size
      call get_expected_obs(seq, i, ens_ptr(k)%state, ens_obs(:, k))
   end do

! Do the update
   call get_new_ens(ens_ptr, ens_obs, obs, obs_err_cov, model_size, ens_size, num_obs_in_set) 

! Output the ensemble prior and posterior to diagnostic files NOT IMPLEMENTED HERE
      !do k = 1, ens_size
      !   call set_single_obs_value(prior_seq, i, j, ens_obs(k), k)
      !   call set_single_obs_value(posterior_seq, i, j, ens_obs(k) + obs_inc(k), k)
      !end do


! Add in posterior state diagnostics
   do j = 1, ens_size
        call output_diagnostics(posterior_state_unit, ens(j), j)
   end do

   ! Deallocate the ens_obs storage for this obs set
   deallocate(obs_err_cov, obs, ens_obs)

end do AdvanceTime


! Initialize the model state space diagnostic output files

! Output the observation space diagnostic files

prior_obs_unit = get_unit()
open(unit = prior_obs_unit, file = 'prior_obs_diagnostics')
call write_obs_sequence(prior_obs_unit, prior_seq)
close(prior_obs_unit)

posterior_obs_unit = get_unit()
open(unit = posterior_obs_unit, file = 'posterior_obs_diagnostics')
call write_obs_sequence(posterior_obs_unit, posterior_seq)
close(posterior_obs_unit)

! Output a restart file if requested
if(output_restart) then
   unit = get_unit()
   open(unit = unit, file = restart_out_file_name)
   do i = 1, ens_size
      call write_state_restart(ens(i), unit)
   end do
   close(unit)
endif



!===========================================================

contains


subroutine get_new_ens(ens_ptr, ens_obs, obs, obs_var, model_size, ens_size, num_obs)
!------------------------------------------------------------------------
! 
! Does the weighted particle filter update with resampling for a set of
! observations with obs_var variance (diagonal Gaussian)

implicit none

integer, intent(in) :: model_size, ens_size, num_obs
type(model_state_ptr_type), intent(inout) :: ens_ptr(ens_size)
real(r8), intent(in) :: ens_obs(num_obs, ens_size), obs(num_obs), obs_var(num_obs)

real(r8) :: weight(ens_size)
integer :: i, j

! This version does obs group together making use of independence of errors
weight = 1.0
do j = 1, num_obs
   do i = 1, ens_size
      weight(i) = weight(i) * exp(-1.0 * (ens_obs(j, i) - obs(j))**2 / &
         (2.0 * obs_var(j)))
   end do
end do

call update_ens_given_weight(ens_ptr, model_size, ens_size, weight)

end subroutine get_new_ens




subroutine update_ens_given_weight(ens_ptr, model_size, ens_size, weight)
!-------------------------------------------------------------------
!
! Given the weights for each ensemble member, does resampling for ensemble

implicit none

integer, intent(in) :: model_size, ens_size
type(model_state_ptr_type), intent(inout) :: ens_ptr(ens_size)
real(r8) :: weight(ens_size)

integer :: i, j, ind, index(ens_size)
real(r8) :: weight_sum, rel_weight(ens_size), new_ens(model_size, ens_size)
real(r8) :: resid_weight(ens_size), cum_weight(0:ens_size), y

! Compute relative weight for each ensemble member
weight_sum = sum(weight)
do i = 1, ens_size
   rel_weight(i) = weight(i) / weight_sum
end do

! First deterministically select integer number of copies
ind = 0
do i = 1, ens_size
   do j = 1, ens_size * rel_weight(i)
      ind = ind + 1
      new_ens(:, ind) = ens_ptr(i)%state
   end do
! Adjust the relative weight
   resid_weight(i) = rel_weight(i) - int(ens_size * rel_weight(i)) / (1.0 * ens_size)
end do

! Get total residual weight
weight_sum = sum(resid_weight)

! Renormalize the residual weights to 1
do i = 1, ens_size
   rel_weight(i) = resid_weight(i) / weight_sum
end do

! Compute cumulative weights at boundaries
cum_weight(0) = 0
do i = 1, ens_size
   cum_weight(i) = cum_weight(i - 1) + rel_weight(i)
end do
! Fix up for round-off error if any
cum_weight(ens_size) = 1.0

! Loop through rest of indices, generate random uniform [0,1] and pick appropriate member
do i = ind + 1, ens_size
   y = random_uniform(random_seq)
   do j = 1, ens_size
      if(y > cum_weight(j - 1) .and. y < cum_weight(j)) then
         new_ens(:, i) = ens_ptr(j)%state
         goto 20
      end if
   end do
20 continue
end do

! Copy the new ensemble back into the ptr array
do i = 1, ens_size
   ens_ptr(i)%state = new_ens(:, i)
end do

end subroutine update_ens_given_weight





function get_ens_mean(ens_ptr, ens_size, index)
!----------------------------------------------------------
!
! Computes the ensemble mean of the index element of the
! state vector.

implicit none

integer,                    intent(in) :: ens_size, index
type(model_state_ptr_type), intent(in) :: ens_ptr(ens_size)
real(r8)                               :: get_ens_mean

integer :: i

get_ens_mean = ens_ptr(1)%state(index)
do i = 2, ens_size
   get_ens_mean = get_ens_mean + ens_ptr(i)%state(index)
end do

get_ens_mean = get_ens_mean / ens_size

end function get_ens_mean

end program leeuwen_particle_filter
