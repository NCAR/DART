! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program filter

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!

use        types_mod, only : r8, missing_r
use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
   get_obs_from_key, set_copy_meta_data, get_copy_meta_data, get_obs_def, get_obs_time_range, &
   get_time_range_keys, set_obs_values, set_obs, write_obs_seq, get_num_obs, &
   get_next_obs, get_num_times, get_obs_values, init_obs, assignment(=), &
   get_num_copies, static_init_obs_sequence, get_qc
use obs_def_mod, only : obs_def_type, get_obs_def_error_variance, get_obs_def_time
use time_manager_mod, only : time_type, set_time, print_time, operator(/=), &
   operator(>)
use    utilities_mod, only :  get_unit, open_file, close_file, register_module, &
                              check_nml_error, file_exist, error_handler, E_ERR, &
                              logfileunit, initialize_utilities, finalize_utilities
use  assim_model_mod, only : assim_model_type, static_init_assim_model, &
   get_model_size, get_closest_state_time_to, &
   advance_state, set_model_time, get_model_time, init_diag_output, &
   output_diagnostics, finalize_diag_output, init_assim_model, get_state_vector_ptr, &
   write_state_restart, read_state_restart, get_state_meta_data, &
   binary_restart_files, aoutput_diagnostics, aread_state_restart, &
   aget_closest_state_time_to, awrite_state_restart, Aadvance_state, pert_model_state
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use  assim_tools_mod, only : obs_increment, update_from_obs_inc, assim_tools_init
use   cov_cutoff_mod, only : comp_cov_factor
use   reg_factor_mod, only : comp_reg_factor
use         sort_mod, only : sort
use    obs_model_mod, only : get_close_states, get_expected_obs

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: seq
type(obs_type)          :: observation, next_obs
type(obs_def_type)      :: obs_def
type(time_type)         :: time1, time2, next_time
type(random_seq_type)   :: random_seq


integer :: i, j, k, ind, iunit, io, istatus
integer :: num_obs_in_set, ierr
integer :: PriorStateUnit, PosteriorStateUnit
integer :: model_size, num_obs_sets
integer :: grp_size, grp_bot, grp_top, group
real(r8) :: reg_factor
real(r8), allocatable :: sum_reg_factor(:, :), reg_factor_series(:, :, :)
real(r8), allocatable :: regress(:), a_returned(:), obs_vals(:)

integer, allocatable :: keys(:)
integer :: key_bounds(2)
integer :: num_state_copies, num_obs_copies
integer :: output_state_mean_index, output_state_spread_index
integer :: prior_obs_mean_index, posterior_obs_mean_index
integer :: prior_obs_spread_index, posterior_obs_spread_index

! Storage for direct access to ensemble state vectors
real(r8),        allocatable :: ens(:, :), ens_mean(:), ens_spread(:), x(:)
type(time_type), allocatable :: ens_time(:)
type(time_type)              :: ens_mean_time, ens_spread_time, x_time

real(r8), allocatable  :: obs_inc(:), ens_inc(:), ens_obs(:), swath(:), rstatus(:,:)
real(r8), allocatable  :: obs_err_cov(:), obs(:), qc(:)
real(r8)               :: cov_factor, obs_mean(1), obs_spread(1)
character(len = 129), allocatable   :: prior_copy_meta_data(:), posterior_copy_meta_data(:)

logical :: interf_provided, out_of_range, is_there_one, is_this_last

! Set a reasonable upper bound on number of close states, will be increased if needed
integer, parameter    :: first_num_close = 100000
integer               :: num_close_ptr(1)
integer,  allocatable :: close_ptr(:, :)         ! First element size should be 1
real(r8), allocatable ::  dist_ptr(:, :)         ! First element size should be 1

!----------------------------------------------------------------
! Namelist input with default values
!
integer  :: async = 0, ens_size = 20
real(r8) :: cutoff      = 0.2_r8
real(r8) :: cov_inflate = 1.0_r8
logical  :: start_from_restart = .false., output_restart = .false.
! if init_time_days and seconds are negative initial time is 0, 0
! for no restart or comes from restart if restart exists
integer  :: init_time_days    = 0
integer  :: init_time_seconds = 0
! Control diagnostic output for state variables
logical  :: output_state_ens_mean = .true., output_state_ens_spread = .true.
logical  :: output_obs_ens_mean   = .true., output_obs_ens_spread   = .true.
integer  :: num_output_state_members = 0
integer  :: num_output_obs_members   = 0
integer  :: output_interval = 1
integer  :: num_groups = 1
real(r8) :: confidence_slope = 0.0_r8
logical  :: get_mean_reg = .false., get_median_reg = .false.

character(len = 129) :: obs_sequence_in_name  = "obs_seq.out",    &
                        obs_sequence_out_name = "obs_seq.final",  &
                        restart_in_file_name  = 'filter_ics',     &
                        restart_out_file_name = 'filter_restart', &
                        adv_ens_command       = './advance_ens.csh'

! adv_ens_command  == 'qsub advance_ens.csh' -> system call advances ensemble by
!                                               qsub submission of a batch job
!                                               -l nodes=# can be inserted after qsub
!                  == './advance_ens.csh'    -> advance ensemble using a script which
!                                               explicitly distributes ensemble among nodes
! advance_ens.csh is currently written to handle both batch submissions (qsub) and
!                 non-batch executions.

namelist /filter_nml/async, adv_ens_command, ens_size, cutoff, cov_inflate, &
   start_from_restart, output_restart, &
   obs_sequence_in_name, obs_sequence_out_name, restart_in_file_name, restart_out_file_name, &
   init_time_days, init_time_seconds, output_state_ens_mean, &
   output_state_ens_spread, output_obs_ens_mean, output_obs_ens_spread, &
   num_output_state_members, num_output_obs_members, output_interval, &
   num_groups, confidence_slope, get_mean_reg, get_median_reg

!----------------------------------------------------------------
! Start of the routine
!----------------------------------------------------------------

call initialize_utilities
call register_module(source,revision,revdate)
write(logfileunit,*)'STARTING filter ...'
call assim_tools_init()

! Initialize the obs sequence module
call static_init_obs_sequence()

! Initialize the observation type variables
call init_obs(observation, 0, 0) 
call init_obs(next_obs, 0, 0)

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
write(logfileunit, nml=filter_nml)

! Now know the ensemble size; allocate all the storage
write(*, *) 'the ensemble size is ', ens_size
allocate(obs_inc(ens_size), ens_inc(ens_size), ens_obs(ens_size), swath(ens_size), rstatus(ens_size,1), &
   prior_copy_meta_data(ens_size + 2), posterior_copy_meta_data(ens_size + 1), &
   regress(num_groups), a_returned(num_groups), obs_vals(ens_size))

! Set an initial size for the close state pointers
allocate(close_ptr(1, first_num_close), dist_ptr(1, first_num_close))

! Determine the number of output obs space fields
num_obs_copies = 2 * num_output_obs_members
if(output_obs_ens_mean) num_obs_copies = num_obs_copies + 1
if(output_obs_ens_spread) num_obs_copies = num_obs_copies + 1

! Read in with enough space for diagnostic output values
call read_obs_seq(obs_sequence_in_name, num_obs_copies, 0, 0, seq)

! Count of number of sets in the sequence
num_obs_sets = get_num_times(seq)

! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()
model_size = get_model_size()
allocate(ens(ens_size, model_size), ens_time(ens_size), ens_mean(model_size))

! Allocate storage for ensemble spread if needed
if(output_state_ens_spread) allocate(ens_spread(model_size))

! Initialize the output sequences and state files and set their meta data
call generate_copy_meta_data(output_state_ens_mean, output_state_ens_spread, &
   num_output_state_members, output_obs_ens_mean, output_obs_ens_spread, &
   num_output_obs_members, num_state_copies, num_obs_copies, output_state_mean_index, &
   output_state_spread_index, prior_obs_mean_index, posterior_obs_mean_index, &
   prior_obs_spread_index, posterior_obs_spread_index, &
   PriorStateUnit, PosteriorStateUnit, seq)

! Set a time type for initial time if namelist inputs are not negative
if(init_time_days >= 0) then
   time1 = set_time(init_time_seconds, init_time_days)
else
   time1 = set_time(0, 0)
endif

!------------------- Read restart if requested ----------------------

if(start_from_restart) then
   iunit = get_unit()
   if (binary_restart_files ) then
      open(unit = iunit, file = restart_in_file_name, form = "unformatted")
   else
      open(unit = iunit, file = restart_in_file_name)
   endif

   do i = 1, ens_size
      write(*, *) 'trying to read restart ', i
      if (binary_restart_files ) then
         call aread_state_restart(ens_time(i), ens(i, :), iunit, "unformatted")
      else
         call aread_state_restart(ens_time(i), ens(i, :), iunit)
      endif

      ! If init_time_days and init_time_seconds are not < 0, set time to them
      if(init_time_days >= 0) ens_time(i) = time1
   end do
   close(iunit)

   !-----------------  Restart read in --------------------------------
else
   !-----  Block to do cold start initialization of ensembles ----------
   ! Initialize the control and ensemble states and set up direct pointers

   ! WARNING: THIS IS COUNTERINTUITIVE: IF START FROM RESTART IS FALSE,
   ! STILL USE A RESTART FILE TO GET SINGLE CONTROL RUN TO PERTURB AROUND.
   allocate(x(model_size))
   iunit = get_unit()
   if (binary_restart_files ) then
      open(unit = iunit, file = restart_in_file_name, form = "unformatted")
   else
      open(unit = iunit, file = restart_in_file_name)
   endif

   ! Get the initial condition
   if (binary_restart_files ) then
      call aread_state_restart(x_time, x, iunit, "unformatted")
   else
      call aread_state_restart(x_time, x, iunit)
   endif
   close(iunit)

   ! Initialize a repeatable random sequence for perturbations
   call init_random_seq(random_seq)

   ! Perturb for ensembles; 
   do i = 1, ens_size
      call pert_model_state(x, ens(i, :), interf_provided)
      ! If model does not provide a perturbing interface, do it here with uniform 0.002
      if(.not. interf_provided) then
         do j = 1, model_size
            ens(i, j) = random_gaussian(random_seq, x(j), 0.002_r8) 
         end do
      endif
      ! Set time to 0, 0 if none specified, otherwise to specified
      ens_time(i) = time1
   end do
   deallocate(x)
   !-------------------- End of cold start ensemble initialization block ------
endif

! Temporary print of initial model time
write(*, *) 'initial model time is '
call print_time(ens_time(1))

! Get the time of the first observation in the sequence
is_there_one = get_first_obs(seq, observation)
! Test for no data at all here?
call get_obs_def(observation, obs_def)
next_time = get_obs_def_time(obs_def)

! Advance the model and ensemble to the closest time to the next
! available observations (need to think hard about these model time interfaces).

AdvanceTime : do i = 1, num_obs_sets

! WANT TO GIVE IT A HINT WHERE TO START WITH OPTIONAL ARGUMENTS FOR BIG OBS SEQUENCES
   call get_obs_time_range(seq, next_time, next_time, key_bounds, num_obs_in_set, out_of_range, observation)
   allocate(keys(num_obs_in_set))
   call get_time_range_keys(seq, key_bounds, num_obs_in_set, keys)
   write(*, *) 'time of obs set ', i
   call print_time(next_time)

   ! If the model time is past the obs set time, just need to skip???
   if(ens_time(1) > next_time) cycle AdvanceTime

   time2 = aget_closest_state_time_to(ens_time(1), next_time)

   ! Advance all the ensembles (to the time of the first ensemble)
   if(time2 /= ens_time(1)) call Aadvance_state(ens_time, ens, ens_size, time2, async, adv_ens_command)

   ! Tag the ensemble mean and spread with the current time
   ens_mean_time   = ens_time(1)
   ens_spread_time = ens_time(1)

   ! Inflate each group separately;  Divide ensemble into num_groups groups
   grp_size = ens_size / num_groups
   Group_inflate: do group = 1, num_groups
      grp_bot = (group - 1) * grp_size + 1
      grp_top = grp_bot + grp_size - 1
      do k = 1, model_size
         ens_mean(k) = sum(ens(grp_bot:grp_top, k)) / grp_size
         do j = grp_bot, grp_top
            ens(j, k) = ens_mean(k) + (ens(j, k) - ens_mean(k)) * sqrt(cov_inflate)
         end do
      end do
   end do Group_inflate

   ! Need global ensemble mean for diagnostics after this
   do k = 1, model_size
      ens_mean(k) = sum(ens(:, k)) / ens_size
   end do

   ! Output state diagnostics as required: NOTE: Prior has been inflated
   if(i / output_interval * output_interval == i) then
      do j = 1, num_output_state_members
         call aoutput_diagnostics( PriorStateUnit, ens_time(j), ens(j, :), j)
      end do
   end if

   ! Output ensemble mean if requested
   if(output_state_ens_mean .and. i / output_interval * output_interval == i) then
      call aoutput_diagnostics(PriorStateUnit, ens_mean_time, ens_mean, output_state_mean_index)
   endif

   ! Compute and output ensemble spread if requested
   if(output_state_ens_spread  .and. i / output_interval * output_interval == i) then
      do k = 1, model_size
         ens_spread(k) = get_ens_spread(ens, ens_mean(k), ens_size, k)
      end do
      call aoutput_diagnostics(PriorStateUnit, ens_spread_time, ens_spread, output_state_spread_index)
   endif


   ! Storage for diagnosing fixed obs set reg_factor
   if(i == 1) then
      if(get_mean_reg) then
         allocate(sum_reg_factor(num_obs_in_set, model_size))
         sum_reg_factor = 0.0_r8
      endif
      if(get_median_reg) then
         allocate(reg_factor_series(num_obs_in_set, model_size, num_obs_sets))
         reg_factor_series = 0.0_r8
      endif
   endif


   ! Allocate storage for the ensemble priors for this number of observations
   allocate(obs_err_cov(num_obs_in_set), obs(num_obs_in_set), qc(num_obs_in_set)) 

   ! Get the observational error covariance (diagonal at present) and the obs values
   do j = 1, num_obs_in_set
      call get_obs_from_key(seq, keys(j), observation)
      call get_obs_def(observation, obs_def)
! ???? Here's where the copy should not be hard coded to 1
      call get_obs_values(observation, obs(j:j), 1)
      call get_qc(observation, qc(j:j), 1)
      obs_err_cov(j) = get_obs_def_error_variance(obs_def)

! Each ensemble member
      do k = 1, ens_size
         call get_expected_obs(seq, keys(j:j),  ens(k, :), obs_vals(k:k))
      enddo

! Compute observation prior ensemble mean
      obs_mean(1) = sum(obs_vals) / ens_size
      obs_spread(1) = sqrt(sum((obs_vals - obs_mean(1))**2) / (ens_size - 1))

! Output all of these that are required to the sequence file
      do k = 1, num_output_obs_members
         call set_obs_values(observation, obs_vals(k:k), 2*k - 1)
      end do
! If requested output the ensemble mean
      if(output_obs_ens_mean) & 
         call set_obs_values(observation, obs_mean, prior_obs_mean_index)

! If requested output the ensemble spread  
      if(output_obs_ens_spread) &
         call set_obs_values(observation, obs_spread, prior_obs_spread_index)

! Finally store the prior observation space stuff for this observation into sequence
      call set_obs(seq, observation, keys(j))
   end do

!------------------------------------------------

   ! Loop through each observation in the set
   Observations : do j = 1, num_obs_in_set
      ! Compute the ensemble prior for this ob
      do k = 1, ens_size
         call get_expected_obs(seq, keys(j:j), ens(k, :), ens_obs(k:k), istatus, rstatus(k:k,1:1))
      end do

      ! Divide ensemble into num_groups groups
      grp_size = ens_size / num_groups
      Group1: do group = 1, num_groups
         grp_bot = (group - 1) * grp_size + 1
         grp_top = grp_bot + grp_size - 1

         ! Call obs_increment to do observation space
         call obs_increment(ens_obs(grp_bot:grp_top), ens_size/num_groups, obs(j), &
            obs_err_cov(j), obs_inc(grp_bot:grp_top), confidence_slope, a_returned(group))
      end do Group1

      ! Getting close states for each scalar observation for now
222   call get_close_states(seq, keys(j), 2.0*cutoff, num_close_ptr(1), close_ptr(1, :), dist_ptr(1, :))
      if(num_close_ptr(1) < 0) then
         deallocate(close_ptr, dist_ptr)
         allocate(close_ptr(1, -1 * num_close_ptr(1)), dist_ptr(1, -1 * num_close_ptr(1)))
         goto 222
      endif

      do k = 1, ens_size
         if (rstatus(k,1) /= 0.0_r8) num_close_ptr(1) = 0
      end do
      if (qc(j) /= 0.0_r8) num_close_ptr(1) = 0

      write(*, *) 'Variables updated for obs ',j,' : ',num_close_ptr(1)

      ! Now loop through each close state variable for this observation
      do k = 1, num_close_ptr(1)
         ind = close_ptr(1, k)
         ! Compute distance dependent envelope
         cov_factor = comp_cov_factor(dist_ptr(1, k), cutoff)

         ! Get the ensemble elements for this state variable and do regression
          swath = ens(:, ind)

         ! Loop through the groups
         Group2: do group = 1, num_groups
            grp_bot = (group - 1) * grp_size + 1
            grp_top = grp_bot + grp_size - 1
            call update_from_obs_inc(ens_obs(grp_bot:grp_top), &
               obs_inc(grp_bot:grp_top), swath(grp_bot:grp_top), ens_size/num_groups, &
               a_returned(group), ens_inc(grp_bot:grp_top), regress(group))
         end do Group2

         ! Compute an information factor for impact of this observation on this state
         if(num_groups == 0) then
             reg_factor = 1.0_r8
         else
            reg_factor = comp_reg_factor(num_groups, regress, i, j, ind)
            if(get_mean_reg) sum_reg_factor(j, k) = sum_reg_factor(j, k) + reg_factor
            if(get_median_reg) reg_factor_series(j, k, i) = reg_factor
         endif

         ! COV_FACTOR NOW COMES IN HERE FOR USE WITH GROUPS; JLA 11/19/03
         reg_factor = min(reg_factor, cov_factor)

         ! Do the final update for this state variable
         ens(:, ind) = ens(:, ind) + reg_factor * ens_inc(:)
      end do

   end do Observations

   ! Free up the storage for this obs set
   deallocate(obs_err_cov, obs)

   ! Output posterior diagnostics

   ! Output state diagnostics as requested
   if(i / output_interval * output_interval == i) then
      do j = 1, num_output_state_members
          call aoutput_diagnostics(     PosteriorStateUnit, ens_time(j), ens(j, :), j)
      end do
   end if

   ! Compute ensemble mean if either mean or spread to be output
   if(output_state_ens_mean .or. output_state_ens_spread  &
      .and. i / output_interval * output_interval == i) then
      do k = 1, model_size
         ens_mean(k) = sum(ens(:, k)) / ens_size
      end do
   endif

!--------------------------------------------------------------------
! Do posterior observation space diagnostics
   do j = 1, num_obs_in_set
      call get_obs_from_key(seq, keys(j), observation)

! Each ensemble member
      do k = 1, ens_size
         call get_expected_obs(seq, keys(j:j),  ens(k, :), obs_vals(k:k))
      enddo

! Compute observation posterior ensemble mean
      obs_mean(1) = sum(obs_vals) / ens_size
      obs_spread(1) = sqrt(sum((obs_vals - obs_mean(1))**2) / (ens_size - 1))

! Output all of these that are required to the sequence file
      do k = 1, num_output_obs_members
         call set_obs_values(observation, obs_vals(k:k), 2*k)
      end do
! If requested output the ensemble mean
      if(output_obs_ens_mean) &
         call set_obs_values(observation, obs_mean, posterior_obs_mean_index)

! If requested output the ensemble spread
      if(output_obs_ens_spread) &
         call set_obs_values(observation, obs_spread, posterior_obs_spread_index)

! Finally store the posterior observation space stuff for this observation into sequence
      call set_obs(seq, observation, keys(j))
   end do

!--------------------------------------------------------------------

   ! Output an ensemble mean in state space if requested
   if(output_state_ens_mean  .and. i / output_interval * output_interval == i) &
      call aoutput_diagnostics(PosteriorStateUnit, ens_mean_time, ens_mean, output_state_mean_index)

   ! Compute and output state_ens_spread if requested
   if(output_state_ens_spread  .and. i / output_interval * output_interval == i) then
      do k = 1, model_size
         ens_spread(k) = get_ens_spread(ens, ens_mean(k), ens_size, k)
      end do
      call aoutput_diagnostics(PosteriorStateUnit, ens_spread_time, ens_spread, output_state_spread_index)
   endif

! Get the next time (if any) in the sequence; ASSUMES THAT OBSERVATION IS THE LAST ONE
   call get_next_obs(seq, observation, next_obs, is_this_last)
   if(is_this_last) exit
   call get_obs_def(next_obs, obs_def)
   next_time = get_obs_def_time(obs_def)

! Deallocate storage used for each set
   deallocate(keys)

end do AdvanceTime

! properly dispose of the diagnostics files

ierr = finalize_diag_output(PriorStateUnit)
ierr = finalize_diag_output(PosteriorStateUnit)

! Output the observation space diagnostic file
call write_obs_seq(seq, obs_sequence_out_name)

! Output a restart file if requested
if(output_restart) then
   iunit = get_unit()
   if (binary_restart_files ) then
      open(unit = iunit, file = restart_out_file_name, form = "unformatted")
      do i = 1, ens_size
         call awrite_state_restart(ens_time(i), ens(i, :), iunit, "unformatted")
      end do
   else
      open(unit = iunit, file = restart_out_file_name)
      do i = 1, ens_size
         call awrite_state_restart(ens_time(i), ens(i, :), iunit)
      end do
   endif
   close(iunit)
endif

! Output the regression factor means if requested
if(get_mean_reg) then
   iunit = get_unit()
   open(unit = iunit, file = 'time_mean_reg')
   write(iunit, *) num_obs_in_set, model_size
   do j = 1, num_obs_in_set
      do i = 1, model_size
         write(iunit, *) j, i, sum_reg_factor(j, i) / num_obs_sets
      end do
   end do
   close(iunit)
endif

! Output the regression factor medians if requested
if(get_median_reg) then
   iunit = get_unit()
   open(unit = iunit, file = 'time_median_reg')
   write(iunit, *) num_obs_in_set, model_size
   do j = 1, num_obs_in_set
      do i = 1, model_size
         reg_factor_series(j, i, :) = sort(reg_factor_series(j, i, :))
         write(iunit, *) j, i, reg_factor_series(j, i, num_obs_sets / 2)
      end do
   end do
   close(iunit)
endif

!===========================================================

write(logfileunit,*)'FINISHED filter.'
write(logfileunit,*)

call finalize_utilities ! closes the log file.

contains



function get_ens_spread(ens, ens_mean, ens_size, index)
!----------------------------------------------------------
!
! Computes the ensemble mean of the index element of the
! state vector.

implicit none

integer,  intent(in) :: ens_size, index
real(r8), intent(in) :: ens_mean
real(r8), intent(in) :: ens(:, :)
real(r8)             :: get_ens_spread

get_ens_spread = sum((ens(:, index) - ens_mean)**2)
get_ens_spread = sqrt(get_ens_spread / (ens_size - 1))

end function get_ens_spread


!-----------------------------------------------------------

subroutine generate_copy_meta_data(output_state_ens_mean, output_state_ens_spread, &
   num_output_state_members, output_obs_ens_mean, output_obs_ens_spread, &
   num_output_obs_members, num_state_copies, num_obs_copies, output_state_mean_index, &
   output_state_spread_index, prior_obs_mean_index, posterior_obs_mean_index, &
   prior_obs_spread_index, posterior_obs_spread_index, &
   PriorStateUnit, PosteriorStateUnit, seq)

implicit none

! Figures out the strings describing the output copies for the three output files.
! THese are the prior and posterior state output files and the observation sequence
! output file which contains both prior and posterior data.

logical, intent(in) :: output_state_ens_mean, output_state_ens_spread
logical, intent(in) :: output_obs_ens_mean, output_obs_ens_spread
integer, intent(in) :: num_output_state_members, num_output_obs_members
integer, intent(out) :: output_state_mean_index, output_state_spread_index
integer, intent(out) :: prior_obs_mean_index, posterior_obs_mean_index
integer, intent(out) :: prior_obs_spread_index, posterior_obs_spread_index
integer, intent(out) :: num_state_copies, num_obs_copies
integer, intent(out) :: PriorStateUnit, PosteriorStateUnit
type(obs_sequence_type), intent(inout) :: seq


character(len=129) :: prior_meta_data, posterior_meta_data
character(len=129) :: state_meta(num_output_state_members + 2)
integer :: i


! Set up the metadata for the output state diagnostic files
do i = 1, num_output_state_members
   if(i < 10000) then
      write(state_meta(i), '(a15, 1x, i6)') 'ensemble member', i
   else
      write(*, *) 'output metadata in filter needs ensemble size < 10000'
      stop
   endif
end do

num_state_copies = num_output_state_members
if(output_state_ens_mean) then
   num_state_copies = num_state_copies + 1
   state_meta(num_state_copies) = 'ensemble mean'
   output_state_mean_index = num_state_copies
endif
if(output_state_ens_spread) then
   num_state_copies = num_state_copies + 1
   state_meta(num_state_copies) = 'ensemble spread'
   output_state_spread_index = num_state_copies
endif


! Set up diagnostic output for model state, if output is desired
if(  output_state_ens_spread .or. output_state_ens_mean .or. &
    ( num_output_state_members > 0 ) ) then
   PriorStateUnit     = init_diag_output('Prior_Diag', &
                           'prior ensemble state', num_state_copies, state_meta)
   PosteriorStateUnit = init_diag_output('Posterior_Diag', &
                           'posterior ensemble state', num_state_copies, state_meta)
endif




! Set up the metadata for the output ensemble observations space file
do i = 1, num_output_obs_members
   if(i < 10000) then
      write(prior_meta_data, '(a21, 1x, i6)') 'prior ensemble member', i
      write(posterior_meta_data, '(a25, 1x, i6)') 'posterior ensemble member', i
   else
      write(*, *) 'output metadata in filter needs ensemble size < 10000'
      stop
   endif
   call set_copy_meta_data(seq, 2*i - 1, prior_meta_data)
   call set_copy_meta_data(seq, 2*i, posterior_meta_data)
end do

num_obs_copies = 2 * num_output_obs_members
if(output_obs_ens_mean) then
   num_obs_copies = num_obs_copies + 1
   prior_meta_data = 'prior ensemble mean'
   call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
   prior_obs_mean_index = num_obs_copies
   num_obs_copies = num_obs_copies + 1
   posterior_meta_data = 'posterior ensemble mean'
   call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
   posterior_obs_mean_index = num_obs_copies 
endif

if(output_obs_ens_spread) then
   num_obs_copies = num_obs_copies + 1
   prior_meta_data = 'prior ensemble spread'
   call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
   prior_obs_spread_index = num_obs_copies
   num_obs_copies = num_obs_copies + 1
   posterior_meta_data = 'posterior ensemble spread'
   call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
   posterior_obs_spread_index = num_obs_copies
endif

end subroutine generate_copy_meta_data

!-------------------------------------------------------------------------



end program filter
