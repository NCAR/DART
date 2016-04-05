! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program  smoother

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

!-----------------------------------------------------------------------------------------
use        types_mod, only : r8, missing_r8
use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, &
   get_obs_from_key, set_copy_meta_data, get_copy_meta_data, get_obs_def, &
   get_time_range_keys, set_obs_values, set_obs, write_obs_seq, get_num_obs, &
   get_obs_values, init_obs, assignment(=), &
   get_num_copies, get_qc, get_num_qc, set_qc, static_init_obs_sequence, destroy_obs, &
   read_obs_seq_header, set_qc_meta_data, get_expected_obs
use obs_def_mod, only : obs_def_type, get_obs_def_error_variance
use time_manager_mod, only : time_type, get_time, set_time, operator(/=), operator(>)
use    utilities_mod, only :  get_unit, open_file, close_file, register_module, &
                              file_exist, error_handler, &
                              E_ERR, E_MSG, E_DBG, initialize_utilities, &
                              logfileunit, timestamp, find_namelist_in_file, check_namelist_read
use  assim_model_mod, only : static_init_assim_model, get_model_size, &
   netcdf_file_type, init_diag_output, finalize_diag_output, & 
   aoutput_diagnostics, aread_state_restart, &
   pert_model_state, open_restart_read, close_restart
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use  assim_tools_smoother_mod, only : assim_tools_init, filter_assim, assim_tools_end
use    obs_model_mod, only : move_ahead
use ensemble_manager_mod, only : init_ensemble_manager, get_ensemble_member, &
   put_ensemble_member, update_ens_mean_spread, end_ensemble_manager, &
   ensemble_type, is_ens_in_core, get_ensemble_time, copy_ens


!-----------------------------------------------------------------------------------------

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(ensemble_type), allocatable     :: ens_handle(:)
type(obs_sequence_type) :: seq
type(obs_type)          :: observation
type(obs_def_type)      :: obs_def
type(time_type)         :: time1
type(random_seq_type)   :: random_seq

character(len=129) :: msgstring
integer :: i, j, iunit, io, days, secs, si
integer :: time_step_number
integer :: num_obs_in_set, ierr, num_qc, last_key_used, model_size
type(netcdf_file_type) :: PriorStateUnit, PosteriorStateUnit
type(netcdf_file_type), allocatable :: SmootherStateUnit(:)
integer :: grp_size, grp_bot, grp_top, group

integer, allocatable :: keys(:)
integer :: key_bounds(2)
integer :: num_state_copies, num_obs_copies, in_obs_copy
integer :: output_state_mean_index, output_state_spread_index
integer :: prior_obs_mean_index, posterior_obs_mean_index
integer :: prior_obs_spread_index, posterior_obs_spread_index

! Storage for direct access to ensemble state vectors
real(r8),        allocatable :: ens_mean(:), temp_ens(:)
type(time_type)              :: ens_mean_time, temp_time

! Storage for use with parallelizable efficient filter
real(r8), allocatable  :: ens_obs(:, :)
real(r8), allocatable  :: obs_err_var(:), obs(:)
character(len = 129), allocatable :: prior_copy_meta_data(:), posterior_copy_meta_data(:)
character(len = 129), allocatable :: smoother_restart_file(:)

logical :: interf_provided
logical, allocatable :: compute_obs(:)

! an integer for writing a file to keep track of step number
integer :: sfile_unit_number 

!----------------------------------------------------------------
! Namelist input with default values
!
integer  :: async = 0, ens_size = 20
real(r8) :: cov_inflate = 1.0_r8
logical  :: start_from_smoother_restart = .false.
logical  :: output_smoother_restart = .true.
logical  :: start_from_restart = .false.
! if init_time_days and seconds are negative -- the time comes from the
! restart file. If you are not starting from a smoother restart file, the initial time is then 0, 0
integer  :: init_time_days    = 0
integer  :: init_time_seconds = 0
! Control diagnostic output for state variables
logical  :: output_state_ens_mean = .true., output_state_ens_spread = .true.
logical  :: output_obs_ens_mean   = .true., output_obs_ens_spread   = .true.
integer  :: num_output_state_members = 0
integer  :: num_output_obs_members   = 0
integer  :: output_interval = 1
integer  :: num_groups = 1
real(r8) :: outlier_threshold = -1.0_r8
integer  :: lag = 0   ! for smoother set to something greater than 0

character(len = 129) :: obs_sequence_in_name  = "obs_seq.out",    &
                        obs_sequence_out_name = "obs_seq.final",  &
                        restart_in_file_name  = 'filter_ics',     &
                        smoother_restart_in_base  = 'smoother_ics',     &
                        smoother_restart_out_base = 'smoother_restart', &
                        adv_ens_command       = './advance_ens.csh'

! adv_ens_command  == 'qsub advance_ens.csh' -> system call advances ensemble by
!                                               qsub submission of a batch job
!                                               -l nodes=# can be inserted after qsub
!                  == './advance_ens.csh'    -> advance ensemble using a script which
!                                               explicitly distributes ensemble among nodes
! advance_ens.csh is currently written to handle both batch submissions (qsub) and
!                 non-batch executions.

namelist /smoother_nml/async, adv_ens_command, ens_size, cov_inflate, &
   start_from_smoother_restart, output_smoother_restart, &
   smoother_restart_in_base, smoother_restart_out_base, &
   start_from_restart, &
   restart_in_file_name, &  
   obs_sequence_in_name, obs_sequence_out_name, &    
   init_time_days, init_time_seconds, output_state_ens_mean, &
   output_state_ens_spread, output_obs_ens_mean, output_obs_ens_spread, &
   num_output_state_members, num_output_obs_members, output_interval, &
   num_groups, outlier_threshold, lag

!----------------------------------------------------------------
! Start of the routine
!----------------------------------------------------------------

! Delete the semaphore files that are used for parallel version 3
call system('rm -f go_advance_model go_end_filter go_assim_regions')

call filter_initialize_modules_used()

! Read the namelist entry
call find_namelist_in_file("input.nml", "smoother_nml", iunit)
read(iunit, nml = smoother_nml, iostat = io)
call check_namelist_read(iunit, io, "smoother_nml")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'smoother_nml','smoother_nml values are',' ',' ',' ')
write(logfileunit, nml=smoother_nml)
write(     *     , nml=smoother_nml)

! Can't output more ensemble members than exist
if(num_output_state_members > ens_size) num_output_state_members = ens_size
if(num_output_obs_members   > ens_size) num_output_obs_members   = ens_size

call filter_alloc_ens_size_storage()

call filter_setup_obs_sequence()

! Allocate model size storage
model_size = get_model_size()
allocate(ens_mean(model_size), temp_ens(model_size))

! Allocate a collection of ensemble objects 
! element 0 is 'lag 0' -- the current ensemble
! element 1 is 'lag 1' -- the lag-1 smoothed posterior
! element 2 is 'lag 2' -- the lag-2 smoothed posterior, etc.
allocate(ens_handle(0:lag), SmootherStateUnit(lag), smoother_restart_file(0:lag))


! Initialize the output sequences and state files and set their meta data
call filter_generate_copy_meta_data()

! Set a time type for initial time if namelist inputs are not negative
call filter_set_initial_time()

! Read restart and initialize storage for 
! Read in the restart for the whole array of ensemble handles
!smoother_restart_file(0)  = "s000"
!smoother_restart_file(9)  = "s009"
if (start_from_smoother_restart) then
   do si = 0, lag
      write(smoother_restart_file(si),'(a,i3.3)')   &
                     trim(adjustl(smoother_restart_in_base)),si
      write(*,*)'smoother restart name is ',trim(adjustl(smoother_restart_file(si)))
      call smoother_read_restart(si)
   enddo
else
   ! do normal filter restart - may or may not have cold initialization
   ! depending on namelist parameter
   ! if namelist start_from_restart = .false. then, have cold restart
   do si = 0, lag
      call filter_read_restart(si)
   enddo 
endif

! Start out with no previously used observations
last_key_used = -99

! Time step number is used to do periodic diagnostic output
time_step_number = 0


AdvanceTime : do

   write(*, *) 'starting advance time loop;'
   time_step_number = time_step_number + 1



   ! Get the model to a good time to use a next set of observations
   call move_ahead(ens_handle(0), ens_size, model_size, seq, last_key_used, &
      key_bounds, num_obs_in_set, async, adv_ens_command)
   if(key_bounds(1) < 0) exit AdvanceTime
   
   if (start_from_smoother_restart) then
   else  
      ! not starting from smoother restart
      ! and we're on the first step, therefore need to 
      ! copy the forecast ensemble
      if (time_step_number == 1) then
         do i = 1, lag
            call copy_ens(ens_handle(0), ens_handle(i), .true.)
         enddo
      endif
   endif


   ! Allocate storage for the ensemble priors for this number of observations
   allocate(keys(num_obs_in_set), obs_err_var(num_obs_in_set), obs(num_obs_in_set), &
            ens_obs(ens_size, num_obs_in_set), compute_obs(num_obs_in_set)) 

   ! For starters allow all obs to be computed as before
   compute_obs = .true.

   ! Get all the keys associated with this set of observations
   call get_time_range_keys(seq, key_bounds, num_obs_in_set, keys)

   ! Inflate each of the groups if required
   ! Uses ens_handle(0) which implies we're only inflating the filtering
   ! solution - in future applications of smoothers it may be prudent
   ! to inflate smoother states ie. in model error cases.
   if(cov_inflate > 1.0_r8) call filter_ensemble_inflate()


   ! Do prior state space diagnostic output as required
   if(time_step_number / output_interval * output_interval == time_step_number) then 
      call filter_state_space_diagnostics(PriorStateUnit)
   endif

   ! Get the observational values, error covariance, and input qc value
   call filter_get_obs_info()

   ! Do prior observation space diagnostics and associated quality control
   ! Uses ens_handle(0)
   call obs_space_diagnostics(ens_obs, ens_size, seq, keys, &
      num_obs_in_set, obs, obs_err_var, outlier_threshold, .true., 0, &
      num_output_obs_members, in_obs_copy + 1, output_obs_ens_mean, &
      prior_obs_mean_index, output_obs_ens_spread, prior_obs_spread_index)

   call filter_assim(ens_handle, ens_obs, compute_obs, ens_size, model_size, num_obs_in_set, &
         num_groups, seq, keys, obs_sequence_in_name, lag)

   ! Do posterior state space diagnostic output as required
   if(time_step_number / output_interval * output_interval == time_step_number) &
      call filter_state_space_diagnostics(PosteriorStateUnit)
      ! above call still generates the usual filter state in same netcdf format

   ! Do the smoother copies if we have non-zero lag
   if (lag > 0) then
   
      if (start_from_smoother_restart) then
         do si = 1, lag
            call filter_state_space_diagnostics(SmootherStateUnit(si), si)
            write(*,*)'Finished outputting lag ',si 
         enddo
      else
         if (time_step_number > lag) then
            do si = 1, lag
               call filter_state_space_diagnostics(SmootherStateUnit(si), si)
               write(*,*)'Finished outputting lag ',si 
            enddo
         endif      
      endif
      ! copies to get appropriate Prior array of ensemble handles
      do si = lag, 1, -1
         ! Time only gets advanced under model propagation 
         write(*,*)' copying  lag ',si-1, ' to ens_handle ',si
         call copy_ens(ens_handle(si-1), ens_handle(si), .true.)
      enddo     

   endif 

   ! Do posterior observation space diagnostics
   call obs_space_diagnostics(ens_obs, ens_size, seq, keys, &
      num_obs_in_set, obs, obs_err_var, outlier_threshold, .false., 2, &
      num_output_obs_members, in_obs_copy + 2, output_obs_ens_mean, &
      posterior_obs_mean_index, output_obs_ens_spread, posterior_obs_spread_index)

   ! Deallocate storage used for each set
   deallocate(keys, obs_err_var, obs, ens_obs, compute_obs)

   ! The last key used is updated to move forward in the observation sequence
   last_key_used = key_bounds(2)

enddo AdvanceTime

! properly dispose of the diagnostics files
ierr = finalize_diag_output(PriorStateUnit)
ierr = finalize_diag_output(PosteriorStateUnit)
do si = 1,lag
ierr = finalize_diag_output(SmootherStateUnit(si))
enddo

! Output the observation space diagnostic file
call write_obs_seq(seq, obs_sequence_out_name)

! Output a restart file if requested
! Note, if run w/ smoother-restart - also have `traditional smoother-restart'
! the '000' restart file is the traditional 'filter_restart'.
do si = 0, lag
   write(smoother_restart_file(si),'(a,i3.3)')   &
                     trim(adjustl(smoother_restart_out_base)),si
   call smoother_output_restart(si)
enddo


! Output the restart for the assim_tools inflation parameters
call assim_tools_end()

write(logfileunit,*)'FINISHED smoother.'
write(logfileunit,*)

deallocate(ens_mean, temp_ens)
deallocate(ens_handle, SmootherStateUnit)

! Send a message to the asynchronous version 3 that all is done
! Must be done always because this also terminates option 3 for assim_tools!
call system('echo a > go_end_filter')

call timestamp(source,revision,revdate,'end') ! That closes the log file, too.

contains

! WARNING: THERE IS SOME DANGER IN USING THESE SCOPED SUBROUTINES
!==========================================================================
!==========================================================================


!-----------------------------------------------------------

subroutine filter_generate_copy_meta_data()

! Figures out the strings describing the output copies for the three output files.
! THese are the prior and posterior state output files and the observation sequence
! output file which contains both prior and posterior data.

character(len=129) :: prior_meta_data, posterior_meta_data, msgstring
character(len=129) :: state_meta(num_output_state_members + 2)
character(len= 18) :: smoothname

integer :: i, ensemble_offset

! Ensemble mean goes first 
num_state_copies = num_output_state_members
if(output_state_ens_mean) then
   num_state_copies = num_state_copies + 1
   output_state_mean_index = 1
   state_meta(output_state_mean_index) = 'ensemble mean'
endif

! Ensemble spread goes second
if(output_state_ens_spread) then
   num_state_copies = num_state_copies + 1
   if(output_state_ens_mean) then
      output_state_spread_index = 2
   else
      output_state_spread_index = 1
   endif
   state_meta(output_state_spread_index) = 'ensemble spread'
endif

! Check for too many output ensemble members
if(num_output_state_members > 10000) then
   write(msgstring, *)'output metadata in filter needs state ensemble size < 10000, not ', &
                      num_output_state_members
   call error_handler(E_ERR,'filter_generate_copy_meta_data',msgstring,source,revision,revdate)
endif

! Compute starting point for ensemble member output
ensemble_offset = 0
if(output_state_ens_mean)   ensemble_offset = ensemble_offset + 1
if(output_state_ens_spread) ensemble_offset = ensemble_offset + 1

! Set up the metadata for the output state diagnostic files
do i = 1, num_output_state_members
   write(state_meta(i + ensemble_offset), '(a15, 1x, i6)') 'ensemble member', i
end do

! Set up diagnostic output for model state, if output is desired
if(  output_state_ens_spread .or. output_state_ens_mean .or. &
    ( num_output_state_members > 0 ) ) then
   PriorStateUnit     = init_diag_output('Prior_Diag', &
                           'prior ensemble state', num_state_copies, state_meta)
   PosteriorStateUnit = init_diag_output('Posterior_Diag', &
                           'posterior ensemble state', num_state_copies, state_meta)

   write(*,*)'Starting lag ',si, ' of ', lag 
   write(logfileunit,*)'Starting lag ',si, ' of ', lag 

   do si = 1,lag
      write(smoothname,'(''Posterior_Lag_'',i4.4)')si
   SmootherStateUnit(si)  = init_diag_output(smoothname, &
                           'smoothed ensemble state', num_state_copies, state_meta, si)

            write(*,*)'Finished initializing lag ',si 
            write(logfileunit,*)'Finished initializing lag ',si 

   enddo
endif


! Set up the metadata for the output ensemble observations space file
! Set up ensemble mean if requested
num_obs_copies = in_obs_copy
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

! Set up ensemble spread if requested
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

! Make sure there are not too many copies requested
if(num_output_obs_members > 10000) then
   write(msgstring, *)'output metadata in filter needs obs ensemble size < 10000, not ',&
                      num_output_obs_members
   call error_handler(E_ERR,'filter_generate_copy_meta_data',msgstring,source,revision,revdate)
endif

! Set up ensemble members as requested
do i = 1, num_output_obs_members
   write(prior_meta_data, '(a21, 1x, i6)') 'prior ensemble member', i
   write(posterior_meta_data, '(a25, 1x, i6)') 'posterior ensemble member', i
   num_obs_copies = num_obs_copies + 1
   call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
   num_obs_copies = num_obs_copies + 1
   call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
end do


end subroutine filter_generate_copy_meta_data

!-------------------------------------------------------------------------

subroutine filter_initialize_modules_used()

! Initialize modules used that require it

call initialize_utilities('Smoother')
call register_module(source,revision,revdate)
call assim_tools_init()

! Initialize the obs sequence module
call static_init_obs_sequence()

! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()

end subroutine filter_initialize_modules_used

!-------------------------------------------------------------------------

subroutine filter_alloc_ens_size_storage()

! Now know the ensemble size; allocate all the storage
write(msgstring, *) 'the ensemble size is ', ens_size
call error_handler(E_MSG,'filter_alloc_ens_size_storage',msgstring,source,revision,revdate)
allocate( prior_copy_meta_data(ens_size + 2), posterior_copy_meta_data(ens_size + 1))

end subroutine filter_alloc_ens_size_storage

!-------------------------------------------------------------------------

subroutine filter_setup_obs_sequence()

integer :: tnum_copies, tnum_qc, tnum_obs, tmax_num_obs, qc_num_inc
real(r8) :: qc(1)
character(len = 129) :: qc_meta_data = 'quality control'
logical :: pre_I_format
character(len = 129) :: obs_seq_read_format
integer :: obs_seq_file_id

! Determine the number of output obs space fields
num_obs_copies = 2 * num_output_obs_members
if(output_obs_ens_mean) num_obs_copies = num_obs_copies + 2
if(output_obs_ens_spread) num_obs_copies = num_obs_copies + 2

! For now, want to have a single qc field, increment if one doesn't exist
call read_obs_seq_header(obs_sequence_in_name, tnum_copies, tnum_qc, tnum_obs, tmax_num_obs, &
   obs_seq_file_id, obs_seq_read_format, pre_I_format, close_the_file = .true.)
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

subroutine filter_set_initial_time

if(init_time_days >= 0) then
   time1 = set_time(init_time_seconds, init_time_days)
else
   time1 = set_time(0, 0)
endif

end subroutine filter_set_initial_time

!-------------------------------------------------------------------------

subroutine filter_read_restart(si)

integer, intent(in) :: si

if(start_from_restart) then
   if(init_time_days >= 0) then
      call init_ensemble_manager(ens_handle(si), ens_size, model_size, restart_in_file_name, &
         time1)
   else
      call init_ensemble_manager(ens_handle(si), ens_size, model_size, restart_in_file_name)
   endif


   !-----------------  Restart read in --------------------------------
else
   !-----  Block to do cold start initialization of ensembles ----------
   ! Initialize the control and ensemble states and set up direct pointers

   ! WARNING: THIS IS COUNTERINTUITIVE: IF START FROM RESTART IS FALSE,
   ! STILL USE A RESTART FILE TO GET SINGLE CONTROL RUN TO PERTURB AROUND.
   call init_ensemble_manager(ens_handle(si), ens_size, model_size)
   iunit = open_restart_read(restart_in_file_name)

   ! Get the initial condition
   ! Read the basic state into ens_mean to conserve storage
   call aread_state_restart(ens_mean_time, ens_mean, iunit)
   call close_restart(iunit)

   ! Initialize a repeatable random sequence for perturbations
   call init_random_seq(random_seq)

   ! Perturb for ensembles; 
   do i = 1, ens_size
      call pert_model_state(ens_mean, temp_ens, interf_provided)
      ! If model does not provide a perturbing interface, do it here with uniform 0.002
      if(.not. interf_provided) then
         do j = 1, model_size
            temp_ens(j) = random_gaussian(random_seq, ens_mean(j), 0.002_r8) 
         end do
      endif
      ! Set this ensemble member 
      call put_ensemble_member(ens_handle(si), i, temp_ens, time1)
   end do
   !-------------------- End of cold start ensemble initialization block ------
endif

! Temporary print of initial model time
call get_time(time1,secs,days)
write(msgstring, *) 'initial model time of first ensemble member (days,seconds) ',days,secs
call error_handler(E_DBG,'filter_read_restart',msgstring,source,revision,revdate)

end subroutine filter_read_restart



 subroutine filter_ensemble_inflate()
!-------------------------------------------------------------------------
!subroutine filter_ensemble_inflate()


! Inflate each group separately;  Divide ensemble into num_groups groups
grp_size = ens_size / num_groups

do group = 1, num_groups
   grp_bot = (group - 1) * grp_size + 1
   grp_top = grp_bot + grp_size - 1
   ens_mean = 0.0_r8
   do j = grp_bot, grp_top
   if(is_ens_in_core(ens_handle(0))) then
      !ens_mean = ens_mean + ens_direct(j, :)
      ens_mean = ens_mean + ens_handle(0)%ens(j, :)
   else
      call get_ensemble_member(ens_handle(0), j, temp_ens, temp_time)
      ens_mean = ens_mean + temp_ens
   endif
   end do
   ens_mean = ens_mean / grp_size
   do j = grp_bot, grp_top
   if(is_ens_in_core(ens_handle(0))) then
      ens_handle(0)%ens(j, :) = ens_mean + sqrt(cov_inflate) * (ens_handle(0)%ens(j, :) - ens_mean)
   else
      call get_ensemble_member(ens_handle(0), j, temp_ens, temp_time)
      temp_ens = ens_mean + sqrt(cov_inflate) * (temp_ens - ens_mean)
      call put_ensemble_member(ens_handle(0), j, temp_ens, temp_time)
   endif
   end do
end do

end subroutine filter_ensemble_inflate



 subroutine filter_state_space_diagnostics(out_unit, lagindx )
!-------------------------------------------------------------------------
!subroutine filter_state_space_diagnostics(out_unit, lagindx )
!
! can write out ANY part of the ens_handle
!
! out_unit
! lagindx    the index to the array of ensemble handles i.e. 'lag'
!
! Right now the time array does not reflect the 'lagged' time.
! That is the intended purpose of the lagindx variable.

implicit none

type(netcdf_file_type), intent(inout) :: out_unit
integer, OPTIONAL,      intent(in)    :: lagindx

integer :: ens_offset
integer :: ind_write


if (present(lagindx)) then
   ind_write = lagindx
else
   ind_write = 0
endif

! Compute ensemble mean and spread if needed for output
if(output_state_ens_mean .or. output_state_ens_spread) &
   call update_ens_mean_spread(ens_handle(ind_write))

! Output ensemble mean if requested
if(output_state_ens_mean) then
   if(is_ens_in_core(ens_handle(ind_write))) then
      call get_ensemble_time(ens_handle(ind_write), 0, temp_time)
      call aoutput_diagnostics(out_unit, temp_time, ens_handle(ind_write)%mean, output_state_mean_index)
   else
      call get_ensemble_member(ens_handle(ind_write), 0, temp_ens, temp_time)
      call aoutput_diagnostics(out_unit, temp_time, temp_ens, output_state_mean_index)
   endif
endif

! Output ensemble spread if requested
if(output_state_ens_spread) then
   if(is_ens_in_core(ens_handle(ind_write))) then
      call get_ensemble_time(ens_handle(ind_write), -1, temp_time)
      call aoutput_diagnostics(out_unit, temp_time, ens_handle(ind_write)%spread, output_state_spread_index)
   else
      call get_ensemble_member(ens_handle(ind_write), -1, temp_ens, temp_time)
      call aoutput_diagnostics(out_unit, temp_time, temp_ens, output_state_spread_index)
   endif
endif

! Compute the offset for copies of the ensemble
ens_offset = 0
if(output_state_ens_mean)   ens_offset = ens_offset + 1
if(output_state_ens_spread) ens_offset = ens_offset + 1

! Output state diagnostics as required: NOTE: Prior has been inflated
do j = 1, num_output_state_members
   if(is_ens_in_core(ens_handle(ind_write))) then
      call get_ensemble_time(ens_handle(ind_write), j, temp_time)
      call aoutput_diagnostics( out_unit, temp_time, ens_handle(ind_write)%ens(j, :), ens_offset + j)
   else
      call get_ensemble_member(ens_handle(ind_write), j, temp_ens, temp_time)
      call aoutput_diagnostics( out_unit, temp_time, temp_ens, ens_offset + j)
   endif
end do

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
call error_handler(E_ERR, 'filter_get_obs_info', &
   'Did not find observation copy with metadata "observations"', &
   source, revision, revdate)

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

subroutine obs_space_diagnostics(ens_obs, ens_size, seq, keys, &
   num_obs_in_set, obs, obs_err_var, outlier_threshold, do_qc, prior_post, &
   num_output_members, members_index, &
   output_ens_mean, ens_mean_index, output_ens_spread, ens_spread_index)

! Do prior observation space diagnostics on the set of obs corresponding to keys

implicit none

integer,  intent(in) :: ens_size
integer,  intent(in) :: num_obs_in_set, keys(num_obs_in_set), prior_post
real(r8), intent(inout) :: ens_obs(ens_size, num_obs_in_set)
integer,  intent(in) :: num_output_members, members_index, ens_mean_index, ens_spread_index
real(r8), intent(in) :: outlier_threshold
real(r8), intent(in) :: obs(num_obs_in_set), obs_err_var(num_obs_in_set)
type(obs_sequence_type), intent(inout) :: seq
logical,  intent(in) :: do_qc
logical,  intent(in) :: output_ens_mean, output_ens_spread

integer :: j, k, istatus, ens_offset
real(r8) :: qc(num_obs_in_set)
real(r8) :: obs_mean(1), obs_spread(1)
real(r8) :: error, diff_sd, ratio
type(obs_type) :: observation
logical :: evaluate_this_ob, assimilate_this_ob

! Construct an observation temporary
call init_obs(observation, get_num_copies(seq), get_num_qc(seq))

ens_obs = 0.0_r8

do k = 1, ens_size
   if(.not. is_ens_in_core(ens_handle(0))) then
      call get_ensemble_member(ens_handle(0), k, temp_ens, temp_time)
   endif
   do j = 1, num_obs_in_set
      call get_obs_from_key(seq, keys(j), observation)
      ! Get the qc value set so far
      if(k == 1) call get_qc(observation, qc(j:j), 1)
      if(is_ens_in_core(ens_handle(0))) then
         call get_expected_obs(seq, keys(j:j), ens_handle(0)%ens(k, :), ens_obs(k, j:j), istatus, &
            assimilate_this_ob, evaluate_this_ob)
      else
         call get_expected_obs(seq, keys(j:j), temp_ens, ens_obs(k, j:j), istatus, &
            assimilate_this_ob, evaluate_this_ob)
      endif


      ! If the observation is not being used at all, also set qc to not able to compute
      if(istatus > 0 .or. (.not. assimilate_this_ob .and. .not. evaluate_this_ob)) then
         if (prior_post == 0) then
            qc(j) = qc(j) + 1000.0_r8       ! Prior
         elseif (prior_post == 2) then
            qc(j) = qc(j) + 1000000.0_r8    ! Posterior
         else
            call error_handler(E_ERR, 'obs_space_diagnostics', &
                 'prior_post value not known.', &
                 source, revision, revdate)
         endif

         ! TEST LINE for doing quality control pass through
         call set_qc(observation, qc(j:j), 1)
         !!! exit
      endif
   end do
end do

do j = 1, num_obs_in_set
   call get_obs_from_key(seq, keys(j), observation)
   ! Compute ensemble mean and spread, zero if qc problem occurred
   if(qc(j) < 1000.0_r8) then
      obs_mean(1) = sum(ens_obs(:, j)) / ens_size
      obs_spread(1) = sqrt(sum((ens_obs(:, j) - obs_mean(1))**2) / (ens_size - 1))
   else
      obs_mean(1) = missing_r8
      obs_spread(1) = 0.0_r8
   endif

   ! This is efficient place to do observation space quality control
   ! For now just looking for outliers from prior
   ! Need to get the observation value for this
   if((outlier_threshold > 0.0_r8) .and. do_qc .and. (qc(j) < 1000.0_r8)) then
      error = obs_mean(1) - obs(j)
      diff_sd = sqrt(obs_spread(1)**2 + obs_err_var(j))
      ratio = abs(error / diff_sd)
      if(ratio > outlier_threshold) qc(j) = qc(j) + 2**prior_post * 100
   endif

   ! If requested output the ensemble mean
   if(output_ens_mean) call set_obs_values(observation, obs_mean, ens_mean_index)
   ! If requested output the ensemble spread
   if(output_ens_spread) call set_obs_values(observation, obs_spread, ens_spread_index)

   ! Compute base location for output of ensemble members
   ens_offset = members_index
   if(output_ens_mean) ens_offset = ens_offset + 2
   if(output_ens_spread) ens_offset = ens_offset + 2
   ! Output all of these ensemble priors that are required to sequence file
   do k = 1, num_output_members
      call set_obs_values(observation, ens_obs(k:k, j), ens_offset + 2 * (k - 1))
   end do

   ! Set the qc value, too
   call set_qc(observation, qc(j:j), 1)

   ! Store the observation into the sequence
   call set_obs(seq, observation, keys(j))
end do

call destroy_obs(observation)

end subroutine obs_space_diagnostics
!-------------------------------------------------------------------------



!-------------------------------------------------------------------------
subroutine smoother_output_restart(si)

implicit none

integer, intent(in) :: si

! Output a restart if requested
if(output_smoother_restart) then
   call end_ensemble_manager(ens_handle(si), smoother_restart_file(si))
else
   call end_ensemble_manager(ens_handle(si))
end if

end subroutine smoother_output_restart
!-------------------------------------------------------------------------


 subroutine smoother_read_restart(si)
!-------------------------------------------------------------------------
!subroutine smoother_read_restart(si)
!
! If the namelist parameters defining the time for the restart are
! 'supplied', then use those as starting values for all the restart files.
! If they are not, (i.e. they are negative) get the restart times from
! the restart files. Those times represent the 'valid time'.

implicit none

integer, intent(in) :: si

type(time_type) :: ensemble_valid_time

if(init_time_days >= 0) then
   call init_ensemble_manager(ens_handle(si), ens_size, model_size, &
                              smoother_restart_file(si), time1)
else
   call init_ensemble_manager(ens_handle(si), ens_size, model_size, &
                              smoother_restart_file(si))
endif

! Print of initial model time
call get_ensemble_time(ens_handle(si), 0, ensemble_valid_time)
call get_time(ensemble_valid_time,secs,days)
write(msgstring,100) si,days,secs
call error_handler(E_MSG,'smoother_read_restart',msgstring,source,revision,revdate)

100 format('ensemble member ',i4,' (days,seconds) = ',i4,1x,i4)

end subroutine smoother_read_restart
!-------------------------------------------------------------------------





end program smoother
