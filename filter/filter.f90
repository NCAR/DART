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

!-----------------------------------------------------------------------------------------
use        types_mod, only : r8
use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
   get_obs_from_key, set_copy_meta_data, get_copy_meta_data, get_obs_def, get_obs_time_range, &
   get_time_range_keys, set_obs_values, set_obs, write_obs_seq, get_num_obs, &
   get_next_obs, get_num_times, get_obs_values, init_obs, assignment(=), &
   get_num_copies, get_qc, get_num_qc, set_qc, static_init_obs_sequence, destroy_obs, &
   read_obs_seq_header, set_qc_meta_data
use obs_def_mod, only : obs_def_type, get_obs_def_error_variance, get_obs_def_time
use time_manager_mod, only : time_type, get_time, set_time, print_time, &
                             operator(/=), operator(>)
use    utilities_mod, only :  get_unit, open_file, close_file, register_module, &
                              check_nml_error, file_exist, error_handler, &
                              E_ERR, E_WARN, E_MSG, E_DBG, initialize_utilities, &
                              logfileunit, finalize_utilities, &
                              timestamp
use  assim_model_mod, only : static_init_assim_model, get_model_size, &
   netcdf_file_type, init_diag_output, output_diagnostics, finalize_diag_output, & 
   binary_restart_files, aoutput_diagnostics, aread_state_restart, &
   awrite_state_restart, pert_model_state
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use  assim_tools_mod, only : obs_increment, update_from_obs_inc, assim_tools_init
use   cov_cutoff_mod, only : comp_cov_factor
use   reg_factor_mod, only : comp_reg_factor
use    obs_model_mod, only : get_close_states, get_expected_obs, move_ahead
!-----------------------------------------------------------------------------------------

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: seq
type(obs_type)          :: observation
type(obs_def_type)      :: obs_def
type(time_type)         :: time1
type(random_seq_type)   :: random_seq

character(len=129) :: msgstring
integer :: i, j, k, ind, iunit, io, istatus, days, secs, reg_series_unit
integer :: num_obs_in_set, ierr, num_qc, last_key_used, model_size
type(netcdf_file_type) :: PriorStateUnit, PosteriorStateUnit
integer :: grp_size, grp_bot, grp_top, group
real(r8) :: reg_factor
real(r8), allocatable :: regress(:), a_returned(:)

integer, allocatable :: keys(:)
integer :: key_bounds(2)
integer :: num_state_copies, num_obs_copies, in_obs_copy
integer :: output_state_mean_index, output_state_spread_index
integer :: prior_obs_mean_index, posterior_obs_mean_index
integer :: prior_obs_spread_index, posterior_obs_spread_index

! Storage for direct access to ensemble state vectors
real(r8),        allocatable :: ens(:, :), ens_mean(:), ens_spread(:)
type(time_type), allocatable :: ens_time(:)
type(time_type)              :: ens_mean_time, ens_spread_time

real(r8), allocatable  :: obs_inc(:), ens_inc(:), ens_obs(:), swath(:)
real(r8), allocatable  :: obs_err_var(:), obs(:)
real(r8)               :: cov_factor, obs_mean(1), obs_spread(1), qc(1)
character(len = 129), allocatable   :: prior_copy_meta_data(:), posterior_copy_meta_data(:)

logical :: interf_provided

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
real(r8) :: outlier_threshold = -1.0_r8
logical  :: save_reg_series = .false.

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
   num_groups, confidence_slope, outlier_threshold, save_reg_series

!----------------------------------------------------------------
! Start of the routine
!----------------------------------------------------------------

call filter_initialize_modules_used()

call filter_read_namelist()

call filter_alloc_ens_size_storage()

call filter_setup_obs_sequence()

call filter_alloc_model_size_storage()

! Initialize the output sequences and state files and set their meta data
call filter_generate_copy_meta_data()

! Set a time type for initial time if namelist inputs are not negative
call filter_set_initial_time()

call filter_read_restart()

! Open an output file for the regression series if requested
if(save_reg_series) then
   reg_series_unit = get_unit()
   open(unit = reg_series_unit, file = 'reg_time_series')
endif

! Start out with no previously used observations
last_key_used = -99

AdvanceTime : do 

   ! Get the model to a good time to use a next set of observations
   call move_ahead(ens, ens_time, ens_size, model_size, seq, last_key_used, &
      key_bounds, num_obs_in_set, async, adv_ens_command)
   if(key_bounds(1) < 0) exit AdvanceTime

   ! Write the divider for the regression series if requested
   if(save_reg_series) write(reg_series_unit, *) -99, -99, -99.0

   ! Allocate storage for the ensemble priors for this number of observations
   allocate(keys(num_obs_in_set), obs_err_var(num_obs_in_set), obs(num_obs_in_set)) 

   ! Get all the keys associated with this set of observations
   call get_time_range_keys(seq, key_bounds, num_obs_in_set, keys)

   ! Tag the ensemble mean and spread with the current time
   ens_mean_time = ens_time(1); ens_spread_time = ens_time(1)

   ! Inflate each of the groups
   call filter_ensemble_inflate()

   ! Do prior state space diagnostic output as required
   if(i / output_interval * output_interval == i)then 
      call filter_state_space_diagnostics(PriorStateUnit)
   endif

   ! Get the observational values, error covariance, and input qc value
   call filter_get_obs_info()

   ! Do prior state space diagnostics and associated quality control
   call obs_space_diagnostics(ens, ens_size, model_size, seq, keys, &
      num_obs_in_set, obs, obs_err_var, outlier_threshold, .true., 0, &
      num_output_obs_members, in_obs_copy + 1, output_obs_ens_mean, &
      prior_obs_mean_index, output_obs_ens_spread, prior_obs_spread_index)

!---------------------- Sequential filter section --------------------------
   ! Loop through each observation in the set
   Observations : do j = 1, num_obs_in_set
      
      ! Get this observation in a temporary and get its qc value for possible update
      call get_obs_from_key(seq, keys(j), observation)
      ! Get the qc value set so far
      call get_qc(observation, qc, 1)
      ! Just skip observations that have failed prior qc ( >4 for now)
      if(qc(1) > 4) cycle Observations

      
      ! Compute the ensemble prior for this ob
      do k = 1, ens_size
         call get_expected_obs(seq, keys(j:j), ens(k, :), ens_obs(k:k), istatus)
         ! Inability to compute forward operator implies skip this observation
         if (istatus > 0) then
            qc(1) = qc(1) + 4000.0_r8
            goto 333
         endif
      end do
      ! Divide ensemble into num_groups groups
      grp_size = ens_size / num_groups
      Group1: do group = 1, num_groups
         grp_bot = (group - 1) * grp_size + 1
         grp_top = grp_bot + grp_size - 1

         ! Call obs_increment to do observation space
         call obs_increment(ens_obs(grp_bot:grp_top), ens_size/num_groups, obs(j), &
            obs_err_var(j), obs_inc(grp_bot:grp_top), confidence_slope, &
            a_returned(group))
      end do Group1

      ! Getting close states for each scalar observation for now
      ! Need model state for some distance computations in sigma, have
      ! to pick some state, only current ones are ensemble, just pass 1
222   call get_close_states(seq, keys(j), 2.0*cutoff, num_close_ptr(1), &
         close_ptr(1, :), dist_ptr(1, :), ens(1, :))
      if(num_close_ptr(1) < 0) then
         deallocate(close_ptr, dist_ptr)
         allocate(close_ptr(1, -1 * num_close_ptr(1)), dist_ptr(1, -1 * num_close_ptr(1)))
         goto 222
      endif

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
         endif
         if(save_reg_series) write(reg_series_unit, *) j, k, reg_factor

         ! COV_FACTOR NOW COMES IN HERE FOR USE WITH GROUPS; JLA 11/19/03
         reg_factor = min(reg_factor, cov_factor)

         ! Do the final update for this state variable
         ens(:, ind) = ens(:, ind) + reg_factor * ens_inc(:)
      end do

!--------------------- End sequential filter section --------------
   ! Write out the updated quality control for this observation
      333 call set_qc(observation, qc, 1)
      call set_obs(seq, observation, keys(j))

   end do Observations

   ! Do prior state space diagnostic output as required
   if(i / output_interval * output_interval == i) &
      call filter_state_space_diagnostics(PosteriorStateUnit)

! Do posterior observation space diagnostics
   call obs_space_diagnostics(ens, ens_size, model_size, seq, keys, &
      num_obs_in_set, obs, obs_err_var, outlier_threshold, .false., 2, &
      num_output_obs_members, in_obs_copy + 2, output_obs_ens_mean, &
      posterior_obs_mean_index, output_obs_ens_spread, posterior_obs_spread_index)

! Deallocate storage used for each set
   deallocate(keys, obs_err_var, obs)

! The last key used is updated to move forward in the observation sequence
   last_key_used = key_bounds(2)

end do AdvanceTime

! properly dispose of the diagnostics files
ierr = finalize_diag_output(PriorStateUnit)
ierr = finalize_diag_output(PosteriorStateUnit)

! Output the observation space diagnostic file
call write_obs_seq(seq, obs_sequence_out_name)

! Output a restart file if requested
call filter_output_restart()

! Close regrssion time series file if needed
if(save_reg_series) close(reg_series_unit)

write(logfileunit,*)'FINISHED filter.'
write(logfileunit,*)

call timestamp(source,revision,revdate,'end') ! closes the log file.

contains

! WARNING: THERE IS SOME DANGER IN USING THESE SCOPED SUBROUTINES
!==========================================================================
!==========================================================================



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

subroutine filter_generate_copy_meta_data()

! Figures out the strings describing the output copies for the three output files.
! THese are the prior and posterior state output files and the observation sequence
! output file which contains both prior and posterior data.

character(len=129) :: prior_meta_data, posterior_meta_data, msgstring
character(len=129) :: state_meta(num_output_state_members + 2)
integer :: i

! Set up the metadata for the output state diagnostic files
do i = 1, num_output_state_members
   if(i < 10000) then
      write(state_meta(i), '(a15, 1x, i6)') 'ensemble member', i
   else
      write(msgstring, *)'output metadata in filter needs state ensemble size < 10000, not ', &
                         num_output_state_members
      call error_handler(E_ERR,'filter_generate_copy_meta_data',msgstring,source,revision,revdate)
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
      write(msgstring, *)'output metadata in filter needs obs ensemble size < 10000, not ',&
                         num_output_obs_members
      call error_handler(E_ERR,'filter_generate_copy_meta_data',msgstring,source,revision,revdate)
   endif
   call set_copy_meta_data(seq, in_obs_copy + 2*i - 1, prior_meta_data)
   call set_copy_meta_data(seq, in_obs_copy + 2*i, posterior_meta_data)
end do

num_obs_copies = in_obs_copy + 2 * num_output_obs_members
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

end subroutine filter_generate_copy_meta_data

!-------------------------------------------------------------------------

subroutine filter_initialize_modules_used()

! Initialize modules used that require it

call initialize_utilities
call register_module(source,revision,revdate)
call error_handler(E_MSG,'filter','STARTING',source,revision,revdate)

call assim_tools_init()

! Initialize the obs sequence module
call static_init_obs_sequence()

! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()

end subroutine filter_initialize_modules_used

!-------------------------------------------------------------------------

subroutine filter_read_namelist()

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

end subroutine filter_read_namelist

!-------------------------------------------------------------------------

subroutine filter_alloc_ens_size_storage()

! Now know the ensemble size; allocate all the storage
write(msgstring, *) 'the ensemble size is ', ens_size
call error_handler(E_MSG,'filter',msgstring,source,revision,revdate)
allocate(obs_inc(ens_size), ens_inc(ens_size), ens_obs(ens_size), swath(ens_size), &
   prior_copy_meta_data(ens_size + 2), &
   posterior_copy_meta_data(ens_size + 1), &
   regress(num_groups), a_returned(num_groups))

! Set an initial size for the close state pointers
allocate(close_ptr(1, first_num_close), dist_ptr(1, first_num_close))

end subroutine filter_alloc_ens_size_storage

!-------------------------------------------------------------------------

subroutine filter_setup_obs_sequence()

integer :: tnum_copies, tnum_qc, tnum_obs, tmax_num_obs, qc_num_inc, num_obs
real(r8) :: qc(1)
character(len = 129) :: qc_meta_data = 'quality control'

! Determine the number of output obs space fields
num_obs_copies = 2 * num_output_obs_members
if(output_obs_ens_mean) num_obs_copies = num_obs_copies + 2
if(output_obs_ens_spread) num_obs_copies = num_obs_copies + 2

! For now, want to have a single qc field, increment if one doesn't exist
call read_obs_seq_header(obs_sequence_in_name, tnum_copies, tnum_qc, tnum_obs, tmax_num_obs)
if(tnum_qc == 0) then
   qc_num_inc = 1
else if(tnum_qc == 1) then
   qc_num_inc = 0
else
   write(*, *) 'Error: filter is only prepared to ingest one qc field for now'
   stop
endif

write(*, *) 'increasing number of qc fields by ', qc_num_inc 

! Read in with enough space for diagnostic output values and add'l qc field
call read_obs_seq(obs_sequence_in_name, num_obs_copies, qc_num_inc, 0, seq)

! Get num of obs copies and num_qc
num_qc = get_num_qc(seq)
in_obs_copy = get_num_copies(seq) - num_obs_copies

! Create an observation type temporaries for use in filter
call init_obs(observation, get_num_copies(seq), num_qc)

! If no qc existed in input file, need to set qc to 0 for all observations
if(qc_num_inc == 1) then
   qc(1) = 0.0_r8
   call set_qc_meta_data(seq, 1, qc_meta_data)
   do i = 1, get_num_obs(seq)
      call get_obs_from_key(seq, i, observation)
      call set_qc(observation, qc, 1)
      call set_obs(seq, observation, i)
   end do
endif


end subroutine filter_setup_obs_sequence

!-------------------------------------------------------------------------

subroutine filter_alloc_model_size_storage()

model_size = get_model_size()
allocate(ens(ens_size, model_size), ens_time(ens_size), ens_mean(model_size))

! Allocate storage for ensemble spread if needed
if(output_state_ens_spread) allocate(ens_spread(model_size))

end subroutine filter_alloc_model_size_storage

!-------------------------------------------------------------------------

subroutine filter_set_initial_time

if(init_time_days >= 0) then
   time1 = set_time(init_time_seconds, init_time_days)
else
   time1 = set_time(0, 0)
endif

end subroutine filter_set_initial_time

!-------------------------------------------------------------------------

subroutine filter_read_restart()

if(start_from_restart) then
   iunit = get_unit()
   if (binary_restart_files ) then
      open(unit = iunit, file = restart_in_file_name, form = "unformatted")
   else
      open(unit = iunit, file = restart_in_file_name)
   endif

   do i = 1, ens_size
      write(msgstring, *) 'trying to read restart ', i
      call error_handler(E_DBG,'filter',msgstring,source,revision,revdate)
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
   iunit = get_unit()
   if (binary_restart_files ) then
      open(unit = iunit, file = restart_in_file_name, form = "unformatted")
   else
      open(unit = iunit, file = restart_in_file_name)
   endif

   ! Get the initial condition
   if (binary_restart_files ) then
   ! Read the basic state into ens_mean to conserve storage
      call aread_state_restart(ens_mean_time, ens_mean, iunit, "unformatted")
   else
      call aread_state_restart(ens_mean_time, ens_mean, iunit)
   endif
   close(iunit)

   ! Initialize a repeatable random sequence for perturbations
   call init_random_seq(random_seq)

   ! Perturb for ensembles; 
   do i = 1, ens_size
      call pert_model_state(ens_mean, ens(i, :), interf_provided)
      ! If model does not provide a perturbing interface, do it here with uniform 0.002
      if(.not. interf_provided) then
         do j = 1, model_size
            ens(i, j) = random_gaussian(random_seq, ens_mean(j), 0.002_r8) 
         end do
      endif
      ! Set time to 0, 0 if none specified, otherwise to specified
      ens_time(i) = time1
   end do
   !-------------------- End of cold start ensemble initialization block ------
endif

! Temporary print of initial model time
call get_time(ens_time(1),secs,days)
write(msgstring, *) 'initial model time of first ensemble member (days,seconds) ',days,secs
call error_handler(E_DBG,'filter',msgstring,source,revision,revdate)

end subroutine filter_read_restart

!-------------------------------------------------------------------------

subroutine filter_ensemble_inflate()

! Inflate each group separately;  Divide ensemble into num_groups groups
grp_size = ens_size / num_groups

do group = 1, num_groups
   grp_bot = (group - 1) * grp_size + 1
   grp_top = grp_bot + grp_size - 1
   do k = 1, model_size
      ens_mean(k) = sum(ens(grp_bot:grp_top, k)) / grp_size
      do j = grp_bot, grp_top
         ens(j, k) = ens_mean(k) + (ens(j, k) - ens_mean(k)) * sqrt(cov_inflate)
      end do
   end do
end do

end subroutine filter_ensemble_inflate

!-------------------------------------------------------------------------

subroutine get_ens_mean_spread(ens, mean, spread, model_size, ens_size, do_spread)

implicit none

real(r8), intent(in) :: ens(ens_size, model_size)
real(r8), intent(out) :: mean(model_size), spread(:)
integer, intent(in) :: model_size, ens_size
logical, intent(in) :: do_spread

integer :: k

do k = 1, model_size
   mean(k) = sum(ens(:, k)) / ens_size
end do

if(do_spread) then
   do k = 1, model_size
      spread(k) = sqrt(sum((ens(:, k) - mean(k))**2) / (ens_size - 1))
   end do
endif


end subroutine get_ens_mean_spread

!-------------------------------------------------------------------------

subroutine filter_state_space_diagnostics(out_unit)

implicit none

type(netcdf_file_type), intent(inout) :: out_unit

! Compute ensemble mean and spread if needed for output
if(output_state_ens_mean .or. output_state_ens_spread) call get_ens_mean_spread(ens, &
   ens_mean, ens_spread, model_size, ens_size, output_state_ens_spread)

! Output state diagnostics as required: NOTE: Prior has been inflated
do j = 1, num_output_state_members
   call aoutput_diagnostics( out_unit, ens_time(j), ens(j, :), j)
end do

! Output ensemble mean if requested
call get_time(ens_time(1),secs,days)
if(output_state_ens_mean) call aoutput_diagnostics(out_unit, &
   ens_mean_time, ens_mean, output_state_mean_index)

! Output ensemble spread if requested
if(output_state_ens_spread) call aoutput_diagnostics(out_unit, &
   ens_spread_time, ens_spread, output_state_spread_index)

end subroutine filter_state_space_diagnostics

!-------------------------------------------------------------------------

subroutine filter_get_obs_info()

integer :: obs_val_index

! Want the obs value to come from field with appropriate meta-data
do j = 1, get_num_copies(seq)
   obs_val_index = j
!   write(*, *) 'meta data copy ', j, trim(get_copy_meta_data(seq, j))
!!!   if(trim(get_copy_meta_data(seq, j)) == 'observations') goto 333
! Need to look for 'observations' or 'NCEP BUFR observation' at present
      if(index(get_copy_meta_data(seq, j), 'observation') > 0) goto 333
end do
! Falling off end means 'observations' not found; die
write(*, *) 'Did not find observation copy with metadata "observations"'
stop

! Get the observational values, error covariance, and input qc value
333 continue
do j = 1, num_obs_in_set
   call get_obs_from_key(seq, keys(j), observation)
   call get_obs_def(observation, obs_def)
! Get the value associated with the observation copy
   call get_obs_values(observation, obs(j:j), obs_val_index)
  obs_err_var(j) = get_obs_def_error_variance(obs_def)
end do

end subroutine filter_get_obs_info

!-------------------------------------------------------------------------

subroutine obs_space_diagnostics(ens, ens_size, model_size, seq, keys, &
   num_obs_in_set, obs, obs_err_var, outlier_threshold, do_qc, prior_post, &
   num_output_members, members_index, &
   output_ens_mean, ens_mean_index, output_ens_spread, ens_spread_index)

! Do prior observation space diagnostics on the set of obs corresponding to keys

implicit none

real(r8), intent(in) :: ens(ens_size, model_size), outlier_threshold
real(r8), intent(in) :: obs(num_obs_in_set), obs_err_var(num_obs_in_set)
type(obs_sequence_type), intent(inout) :: seq
integer, intent(in) :: ens_size, model_size
logical, intent(in) :: do_qc
integer, intent(in) :: num_obs_in_set, keys(num_obs_in_set), prior_post
integer, intent(in) :: num_output_members, members_index, ens_mean_index, ens_spread_index
logical, intent(in) :: output_ens_mean, output_ens_spread

integer :: j, k, istatus
real(r8) :: obs_vals(ens_size), qc(1), obs_mean(1), obs_spread(1)
real(r8) :: error, diff_sd, ratio
type(obs_type) :: observation

! Construnct an observation temporary
call init_obs(observation, get_num_copies(seq), get_num_qc(seq))

do j = 1, num_obs_in_set
   call get_obs_from_key(seq, keys(j), observation)
   obs_vals = 0.0 
   ! Get the qc value set so far
   call get_qc(observation, qc, 1)
   do k = 1, ens_size
      call get_expected_obs(seq, keys(j:j), ens(k, :), obs_vals(k:k), istatus)
      if(istatus > 0) then 
         qc(1) = qc(1) + 2**prior_post * 1000
         exit
      endif
   end do
   ! Compute ensemble mean and spread, zero if qc problem occurred
   obs_mean(1) = sum(obs_vals) / ens_size
   obs_spread(1) = sqrt(sum((obs_vals - obs_mean(1))**2) / (ens_size - 1))

   ! This is efficient place to do observation space quality control
   ! For now just looking for outliers from prior
   ! Need to get the observation value for this
   if(outlier_threshold > 0.0 .and. do_qc) then
      error = obs_mean(1) - obs(j)
      diff_sd = sqrt(obs_spread(1)**2 + obs_err_var(j))
      ratio = abs(error / diff_sd)
      if(ratio > outlier_threshold) qc(1) = qc(1) + 2**prior_post * 10000
   endif

   ! Output all of these ensemble priors that are required to sequence file
   do k = 1, num_output_members
      call set_obs_values(observation, obs_vals(k:k), members_index + 2 * (k - 1))
   end do

   ! If requested output the ensemble mean
   if(output_ens_mean) call set_obs_values(observation, obs_mean, ens_mean_index)
   ! If requested output the ensemble spread
   if(output_ens_spread) call set_obs_values(observation, obs_spread, ens_spread_index)
   
   ! Set the qc value, too
   call set_qc(observation, qc, 1)

   ! Store the observation into the sequence
   call set_obs(seq, observation, keys(j))
end do

call destroy_obs(observation)

end subroutine obs_space_diagnostics

!-------------------------------------------------------------------------

subroutine filter_output_restart()

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

end subroutine filter_output_restart

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

end program filter
