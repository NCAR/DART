! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program filter

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$
!

use        types_mod, only : r8
use obs_sequence_mod, only : obs_sequence_type, write_obs_sequence, &
       read_obs_sequence, get_num_obs_sets, get_obs_sequence_time, &
       get_num_obs_in_set, get_expected_obs, get_diag_obs_err_cov, &
       get_obs_values, obs_sequence_def_copy, inc_num_obs_copies, &
       set_single_obs_value, get_num_close_states, get_close_states
use time_manager_mod, only : time_type, set_time, print_time, operator(/=), &
                             operator(>)
use    utilities_mod, only : get_unit, open_file, close_file, &
                             check_nml_error, file_exist
use  assim_model_mod, only : assim_model_type, static_init_assim_model, &
       get_model_size, get_closest_state_time_to, &
       advance_state, set_model_time, get_model_time, init_diag_output, &
       output_diagnostics, init_assim_model, get_state_vector_ptr, &
       write_state_restart, read_state_restart, get_state_meta_data
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use  assim_tools_mod, only : obs_increment, update_from_obs_inc, &
       look_for_bias, obs_increment17
use   cov_cutoff_mod, only : comp_cov_factor
use     location_mod, only : location_type

use typesizes         ! From the F90 netCDF interface
use netcdf            ! From the F90 netCDF interface

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

! Define a type for doing direct access to ensemble state vectors
type model_state_ptr_type
   private
   real(r8), pointer :: state(:)
end type model_state_ptr_type

type(obs_sequence_type) :: seq, prior_seq, posterior_seq
type(time_type)         :: time, time2
type(random_seq_type)   :: random_seq


integer :: i, j, k, ind, iunit, prior_obs_unit, posterior_obs_unit, io
integer :: prior_state_unit, posterior_state_unit, num_obs_in_set, ierr
integer :: PriorStateUnit, PosteriorStateUnit
integer :: lji, meta_data_size
integer :: output_ens_mean_index, output_ens_spread_index
integer :: model_size, num_obs_sets

! Storage for direct access to ensemble state vectors
type(model_state_ptr_type), allocatable :: ens_ptr(:)
type(model_state_ptr_type) :: x_ptr, ens_mean_ptr, ens_spread_ptr

type(assim_model_type), allocatable :: ens(:)
type(assim_model_type) :: x, ens_mean, ens_spread
real(r8), allocatable  :: obs_inc(:), ens_inc(:), ens_obs(:), swath(:)
real(r8), allocatable  :: obs_err_cov(:), obs(:)
real(r8)               :: cov_factor, mean_inc, sd_ratio
character(len = 129), allocatable   :: ens_copy_meta_data(:)

integer,  allocatable :: num_close_ptr(:), close_ptr(:, :), num_close(:)
real(r8), allocatable :: dist_ptr(:, :)

! Test storage for variance ratio
real(r8) :: var_ratio_sum, var_ratio

! Temporary storage to allow decent CAM initial ensemble perturbations
integer :: var_type
type(location_type) :: location

! Temporary storage to test adaptive error capability; should be moved to assim_tools
integer  :: slope_index = 0
real(r8) :: slope       = 0.0_r8

!----------------------------------------------------------------
! Namelist input with default values
!
integer  :: async       = 0, ens_size = 20
real(r8) :: cutoff      = 200.0_r8
real(r8) :: cov_inflate = 1.0_r8
integer  :: cache_size  = 10
logical  :: start_from_restart = .false., output_restart = .false.
! if init_time_days and seconds are negative initial time is 0, 0
! for no restart or comes from restart if restart exists
integer  :: init_time_days    = -1
integer  :: init_time_seconds = -1
! Control diagnostic output for state variables
logical  :: output_state_ens_mean = .true., output_state_ens_spread = .true.
integer  :: num_output_ens_members = 0
integer  :: output_interval = 1

character(len = 129) :: obs_sequence_file_name = "obs_sequence", &
                        restart_in_file_name = 'filter_restart_in', &
                        restart_out_file_name = 'filter_restart_out'

namelist /filter_nml/async, ens_size, cutoff, cov_inflate, cache_size, &
   start_from_restart, output_restart, &
   obs_sequence_file_name, restart_in_file_name, restart_out_file_name, &
   init_time_days, init_time_seconds, output_state_ens_mean, &
   output_state_ens_spread, num_output_ens_members, output_interval
!----------------------------------------------------------------

! Begin by reading the namelist input
if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(iunit, nml = filter_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'filter_nml')
   enddo
 11 continue
   call close_file(iunit)
endif

! Now know the ensemble size; allocate all the storage
write(*, *) 'the ensemble size is ', ens_size
allocate(ens_ptr(ens_size), ens(ens_size), obs_inc(ens_size), &
   ens_inc(ens_size), ens_obs(ens_size), swath(ens_size), &
   ens_copy_meta_data(ens_size + 2))

! Input the obs_sequence
iunit = get_unit()
open(unit = iunit, file = obs_sequence_file_name)
seq = read_obs_sequence(iunit)
close(iunit)
! Count of number of sets in the sequence
num_obs_sets = get_num_obs_sets(seq)

! Copy just the definitions part of the sequence to the two output obs sequences
!!!call obs_sequence_def_copy(prior_seq, seq)
!!!call obs_sequence_def_copy(posterior_seq, seq)
! Set up the metadata for the output ensemble observations
do i = 1, num_output_ens_members
   write(ens_copy_meta_data(i), *) 'ensemble ', i
end do
meta_data_size = num_output_ens_members
if(output_state_ens_mean) then
   meta_data_size = meta_data_size + 1
   ens_copy_meta_data(meta_data_size) = 'ensemble mean'
   output_ens_mean_index = meta_data_size
endif
if(output_state_ens_spread) then
   meta_data_size = meta_data_size + 1
   ens_copy_meta_data(meta_data_size) = 'ensemble spread'
   output_ens_spread_index = meta_data_size
endif


! For now output all ensemble members for prior and posterior; add space
!!!call inc_num_obs_copies(prior_seq, ens_size, ens_copy_meta_data)
!!!call inc_num_obs_copies(posterior_seq, ens_size, ens_copy_meta_data)

! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()
model_size = get_model_size()

! Allocate storage for ensemble mean (and spread if needed)
call init_assim_model(ens_mean)
ens_mean_ptr%state => get_state_vector_ptr(ens_mean)
if(output_state_ens_spread) then
   call init_assim_model(ens_spread)
   ens_spread_ptr%state => get_state_vector_ptr(ens_spread)
endif

! Set up diagnostic output for model state, if output is desired

if(  output_state_ens_spread .or. output_state_ens_mean .or. &
    ( num_output_ens_members > 0 ) ) then
   PriorStateUnit     = init_diag_output('Prior_Diag', &
                           'prior ensemble state', meta_data_size, ens_copy_meta_data)
   PosteriorStateUnit = init_diag_output('Posterior_Diag', &
                           'posterior ensemble state', meta_data_size, ens_copy_meta_data)
endif

! Set a time type for initial time if namelist inputs are not negative
if(init_time_days >= 0) then
   time = set_time(init_time_seconds, init_time_days)
else
   time = set_time(0, 0)
endif

!------------------- Read restart if requested ----------------------
if(start_from_restart) then
   iunit = get_unit()
   open(unit = iunit, file = restart_in_file_name)
   call init_assim_model(x)
   x_ptr%state => get_state_vector_ptr(x)
   do i = 1, ens_size
      write(*, *) 'trying to read restart ', i
      call init_assim_model(ens(i))
      ens_ptr(i)%state => get_state_vector_ptr(ens(i))

      call read_state_restart(ens(i), iunit)
      ! If init_time_days an init_time_seconds are not < 0, set time to them
      if(init_time_days >= 0) call set_model_time(ens(i) , time)
   end do
   close(iunit)
   !-----------------  Restart read in --------------------------------
else
   !-----  Block to do cold start initialization of ensembles ----------
   ! Initialize the control and ensemble states and set up direct pointers

   ! WARNING: THIS IS COUNTERINTUITIVE: IF START FROM RESTART IS FALSE,
   ! STILL USE A RESTART FILE TO GET SINGLE CONTROL RUN TO PERTURB AROUND.
   iunit = get_unit()
   open(unit = iunit, file = restart_in_file_name)
   call init_assim_model(x)
   x_ptr%state => get_state_vector_ptr(x)

   do i = 1, ens_size
      call init_assim_model(ens(i))
      ens_ptr(i)%state => get_state_vector_ptr(ens(i))
   end do

   ! Get the initial condition
   call read_state_restart(x, iunit)
   close(iunit)

   ! Initialize a repeatable random sequence for perturbations
   ! Where should the magnitude of the perturbations come from here???
   call init_random_seq(random_seq)
   ! Perturb for ensembles; 
   do i = 1, ens_size
      do j = 1, model_size
         ! TEMPORARY KLUGE FOR GETTING CAM ROLLING: NEED A PERTURB_MODEL_ENS interface
         call get_state_meta_data(j, location, var_type)
         if(var_type < 4) then 
            ens_ptr(i)%state(j) = random_gaussian(random_seq, x_ptr%state(j), 2.0_r8) 
         else
            ens_ptr(i)%state(j) = random_gaussian(random_seq, x_ptr%state(j), 0.000001_r8) 
         endif
         
      end do
      ! Set time to 0, 0 if none specified, otherwise to specified
      call set_model_time(ens(i), time)
   end do
   !-------------------- End of cold start ensemble initialization block ------
endif

! Temporary print of initial model time
write(*, *) 'initial model time is '
call print_time(get_model_time(ens(1)))


! Advance the model and ensemble to the closest time to the next
! available observations (need to think hard about these model time interfaces).

AdvanceTime : do i = 1, num_obs_sets

   call get_obs_sequence_time(seq, i, time)
   write(*, *) ' '
   write(*, *) 'time of obs set ', i
   call print_time(time)

   ! If the model time is past the obs set time, just need to skip???
   if(get_model_time(ens(1)) > time) cycle AdvanceTime

   time2 = get_closest_state_time_to(ens(1), time)
   !   write(*, *) 'advancing to time2 '
   call  print_time(time2)

   ! Advance all the ensembles (to the time of the first ensemble)
   if(time2 /= get_model_time(ens(1))) call advance_state(ens, ens_size, time2, async)

   ! Tag the ensemble mean and spread with the current time
   ens_mean%time   = get_model_time(ens(1))
   ens_spread%time = get_model_time(ens(1))

   ! Do a covariance inflation for now? 
   ! Inflate the ensemble state estimates
   do k = 1, model_size
      ens_mean_ptr%state(k) = get_ens_mean(ens_ptr, ens_size, k)
      do j = 1, ens_size
         ens_ptr(j)%state(k) = ens_mean_ptr%state(k) + (ens_ptr(j)%state(k) - &
            ens_mean_ptr%state(k)) * sqrt(cov_inflate)
      end do
   end do

   ! Output state diagnostics as required: NOTE: Prior has been inflated
   if(i / output_interval * output_interval == i) then
      do j = 1, num_output_ens_members
         ! TJH debugging block.
         print *,'Going one time, i,noem= ',j,num_output_ens_members
         call output_diagnostics(     PriorStateUnit, ens(j), j)
      end do
   end if

   ! Output ensemble mean if requested
   if(output_state_ens_mean .and. i / output_interval * output_interval == i) &
      call output_diagnostics(PriorStateUnit, ens_mean, output_ens_mean_index)

   ! Compute and output ensemble spread if requested
   if(output_state_ens_spread  .and. i / output_interval * output_interval == i) then
      do k = 1, model_size
         ens_spread_ptr%state(k) = get_ens_spread(ens_ptr, &
            ens_mean_ptr%state(k), ens_size, k)
      end do
      call output_diagnostics(PriorStateUnit, ens_spread, output_ens_spread_index)
   endif

   ! How many observations in this set
   num_obs_in_set = get_num_obs_in_set(seq, i)

   ! Allocate storage for the ensemble priors for this number of observations
   allocate(obs_err_cov(num_obs_in_set), obs(num_obs_in_set), &
      num_close(num_obs_in_set), num_close_ptr(num_obs_in_set))

   ! Get the observational error covariance (diagonal at present)
   call get_diag_obs_err_cov(seq, i, obs_err_cov)

   ! Get the observations; from copy 1 for now
   call get_obs_values(seq, i, obs, 1)


   ! Following block replaces cache block -----------------------------
   write(*, *) 'calling get_close_states'
   call get_num_close_states(seq, i, 2.0*cutoff, num_close)

   ! Allocate storage for num_close, close, and distance
   allocate(close_ptr(num_obs_in_set, maxval(num_close)), &
      dist_ptr(num_obs_in_set, maxval(num_close)))

   ! Now get the values
   call get_close_states(seq, i, 2.0*cutoff, num_close_ptr, close_ptr, dist_ptr)
   write(*, *) 'back form get_close_states'

!--------------------------------------------------------------------

! A preliminary search for bias???
   var_ratio_sum = 0.0
!!!   do j = 1, num_obs_in_set
! Get all the prior estimates
!!!      do k = 1, ens_size
!!!         call get_expected_obs(seq, i, ens_ptr(k)%state, ens_obs(k:k), j)
!!!      end do
! Call a subroutine to evaluate how many S.D.s away the ob value is
!!!      call look_for_bias(ens_obs, ens_size, obs(j), obs_err_cov(j), var_ratio)
!!!      var_ratio_sum = var_ratio_sum + var_ratio
!!!   end do
  
!!!   write(*, *) 'mean var_ratio is ', var_ratio_sum / num_obs_in_set
   
! Adjust the slope for error_correcting filter given var_ratio
!!!   if(var_ratio_sum / num_obs_in_set > 1.0 .and. slope < 1.0) then
!!!      slope_index = slope_index + 1
!!!   endif

!!!   if(var_ratio_sum / num_obs_in_set < 1.00) then
!!!      slope_index = slope_index - 1
!!!   endif
!!! !!!   slope = (slope_index * 0.02) ** 2
!!!   if(slope_index > 3) then 
!!!      slope = ((slope_index - 3) * 0.04) ** 2
!!!   else if(slope_index < -3) then
!!!      slope = (slope_index + 3) * 0.01
!!!   else
!!!      slope = 0.0
!!!   endif


!!! LOOKING AT TUNING FOR PS ONLY
!!!slope = 0.05
!!!   write(*, *) 'UPDATED SLOPE IS ', slope


   ! Loop through each observation in the set
   Observations : do j = 1, num_obs_in_set
      ! Compute the ensemble prior for this ob
      do k = 1, ens_size
         call get_expected_obs(seq, i, ens_ptr(k)%state, ens_obs(k:k), j)
      end do

      call obs_increment(ens_obs, ens_size, obs(j), obs_err_cov(j), obs_inc)
!!!      call obs_increment17(ens_obs, ens_size, obs(j), obs_err_cov(j), obs_inc, slope)

      ! Output the ensemble prior and posterior to diagnostic files
      do k = 1, ens_size
!!!         call set_single_obs_value(prior_seq, i, j, ens_obs(k), k)
!!!         call set_single_obs_value(posterior_seq, i, j, ens_obs(k) + obs_inc(k), k)
      end do

      ! Now loop through each close state variable for this observation
      do k = 1, num_close_ptr(j)
         ind = close_ptr(j, k)
         ! Compute distance dependent envelope
         cov_factor = comp_cov_factor(dist_ptr(j, k), cutoff)

         ! Get the ensemble elements for this state variable and do regression
         swath = get_ens_swath(ens_ptr, ens_size, ind)

         call update_from_obs_inc(ens_obs, obs_inc, &
                   swath, ens_size, ens_inc, cov_factor)

         call inc_ens_swath(ens_ptr, ens_size, ind, ens_inc)
      end do

   end do Observations

   ! Free up the storage for close
   deallocate(close_ptr, dist_ptr)
   deallocate(obs_err_cov, obs, num_close, num_close_ptr)

   ! Output posterior diagnostics
   ! Output state diagnostics as requested
   if(i / output_interval * output_interval == i) then
      do j = 1, num_output_ens_members
          call output_diagnostics(     PosteriorStateUnit, ens(j), j)
      end do
   end if

   ! Compute ensemble mean if either mean or spread to be output
   if(output_state_ens_mean .or. output_state_ens_spread  &
      .and. i / output_interval * output_interval == i) then
      do k = 1, model_size
         ens_mean_ptr%state(k) = get_ens_mean(ens_ptr, ens_size, k)
      end do
   endif

   ! Output an ensemble mean if requested
   if(output_state_ens_mean  .and. i / output_interval * output_interval == i) &
      call output_diagnostics(PosteriorStateUnit, ens_mean, output_ens_mean_index)

   ! Compute and output state_ens_spread if requested
   if(output_state_ens_spread  .and. i / output_interval * output_interval == i) then
      do k = 1, model_size
         ens_spread_ptr%state(k) = get_ens_spread(ens_ptr, &
            ens_mean_ptr%state(k), ens_size, k)
      end do
      call output_diagnostics(PosteriorStateUnit, ens_spread, output_ens_spread_index)
   endif

   ! ierr = NF90_sync(PosteriorStateUnit)   ! just for good measure -- TJH 

end do AdvanceTime

! properly dispose of the diagnostics files

ierr = NF90_close(PriorStateUnit)
ierr = NF90_close(PosteriorStateUnit)

! Initialize the model state space diagnostic output files

! Output the observation space diagnostic files

!!!prior_obs_unit = get_unit()
!!!open(unit = prior_obs_unit, file = 'prior_obs_diagnostics')
!!!call write_obs_sequence(prior_obs_unit, prior_seq)
!!!close(prior_obs_unit)

!!!posterior_obs_unit = get_unit()
!!!open(unit = posterior_obs_unit, file = 'posterior_obs_diagnostics')
!!!call write_obs_sequence(posterior_obs_unit, posterior_seq)
!!!close(posterior_obs_unit)

! Output a restart file if requested
if(output_restart) then
   iunit = get_unit()
   open(unit = iunit, file = restart_out_file_name)
   do i = 1, ens_size
      call write_state_restart(ens(i), iunit)
   end do
   close(iunit)
endif



!===========================================================

contains


function get_ens_swath(ens_ptr, ens_size, stateindex)
!----------------------------------------------------------
!
! Returns the ensemble for state variable index

implicit none

integer,                    intent(in) :: ens_size, stateindex
type(model_state_ptr_type), intent(in) :: ens_ptr(ens_size)
real(r8)                               :: get_ens_swath(ens_size)

integer :: i

do i = 1, ens_size
   get_ens_swath(i) = ens_ptr(i)%state(stateindex)
end do

end function get_ens_swath

         

subroutine inc_ens_swath(ens_ptr, ens_size, stateindex, ens_inc)
!-----------------------------------------------------------
!
! Increments all ensemble members for a given state variable
! index.

implicit none

integer,                    intent(in)    :: ens_size, stateindex
type(model_state_ptr_type), intent(inout) :: ens_ptr(ens_size)
real(r8),                   intent(in)    :: ens_inc(ens_size)

integer :: i

do i = 1, ens_size
   ens_ptr(i)%state(stateindex) = ens_ptr(i)%state(stateindex) + ens_inc(i)
end do

end subroutine inc_ens_swath




function get_ens_mean(ens_ptr, ens_size, stateindex)
!----------------------------------------------------------
!
! Computes the ensemble mean of the index element of the
! state vector.

implicit none

integer,                    intent(in) :: ens_size, stateindex
type(model_state_ptr_type), intent(in) :: ens_ptr(ens_size)
real(r8)                               :: get_ens_mean

integer :: i

get_ens_mean = ens_ptr(1)%state(stateindex)
do i = 2, ens_size
   get_ens_mean = get_ens_mean + ens_ptr(i)%state(stateindex)
end do

get_ens_mean = get_ens_mean / ens_size

end function get_ens_mean




function get_ens_spread(ens_ptr, ens_mean, ens_size, stateindex)
!----------------------------------------------------------
!
! Computes the ensemble mean of the index element of the
! state vector.

implicit none

integer,                    intent(in) :: ens_size, stateindex
real(r8),                   intent(in) :: ens_mean
type(model_state_ptr_type), intent(in) :: ens_ptr(ens_size)
real(r8)                               :: get_ens_spread

integer :: i

get_ens_spread = (ens_ptr(1)%state(stateindex) - ens_mean)**2
do i = 2, ens_size
   get_ens_spread = get_ens_spread + (ens_ptr(i)%state(stateindex) - ens_mean)**2
end do

get_ens_spread = sqrt(get_ens_spread / (ens_size - 1))

end function get_ens_spread


end program filter
