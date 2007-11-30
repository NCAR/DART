! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module smoother_mod 

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
!
! Tools for turning the filter into a fixed lag smoother for the full state vector.

use      types_mod,       only : r8
use  mpi_utilities_mod,   only : my_task_id
use  utilities_mod,       only : file_exist, get_unit, check_namelist_read, do_output, &
                                 find_namelist_in_file, register_module, error_handler, &
                                 E_ERR, E_MSG, logfileunit
use ensemble_manager_mod, only : ensemble_type, init_ensemble_manager, read_ensemble_restart, &
                                 write_ensemble_restart, all_vars_to_all_copies, &
                                 duplicate_ens, compute_copy_mean, compute_copy_mean_sd, &
                                 all_copies_to_all_vars, get_copy
use time_manager_mod,     only : time_type
use assim_model_mod,      only : static_init_assim_model, get_model_size,                    &
                                 netcdf_file_type, init_diag_output, finalize_diag_output,   &
                                 aoutput_diagnostics
use assim_tools_mod,      only : filter_assim
use obs_sequence_mod,     only : obs_sequence_type
use adaptive_inflate_mod, only : adaptive_inflate_type, adaptive_inflate_init, &
                                 do_varying_ss_inflate, do_single_ss_inflate
implicit none
private

public :: smoother_read_restart, advance_smoother,                          &
   smoother_gen_copy_meta_data, smoother_write_restart, init_smoother, &
   do_smoothing, smoother_mean_spread, smoother_assim,                      &
   filter_state_space_diagnostics, smoother_ss_diagnostics,        &
   smoother_end, smoother_inc_lags


! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical :: module_initialized = .false.

character(len = 129) :: errstring

! State for the smoother
integer :: smoother_head = 1, num_current_lags = 0
integer :: smoother_state_mean_index, smoother_state_spread_index
type(ensemble_type), allocatable :: lag_handle(:)
type(netcdf_file_type), allocatable :: SmootherStateUnit(:)
type(adaptive_inflate_type)         :: lag_inflate

!============================================================================

!---- namelist with default values

integer  :: num_lags = 0
logical  :: start_from_restart = .false.
logical  :: output_restart     = .false.

character(len = 129) :: restart_in_file_name  = 'smoother_ics', &
                        restart_out_file_name = 'smoother_restart'

namelist /smoother_nml/ num_lags, start_from_restart, &
                        output_restart, restart_in_file_name, &
                        restart_out_file_name

!-------------------------------------------------------------------------
contains


subroutine static_init_smoother()

integer :: iunit, io

! First call to init_smoother must initialize module and read namelist
if ( .not. module_initialized ) then
   ! Initialize the module with utilities
   call register_module(source, revision, revdate)
   module_initialized = .true.

   ! Read the namelist entry
   call find_namelist_in_file("input.nml", "smoother_nml", iunit)
   read(iunit, nml = smoother_nml, iostat = io)
   call check_namelist_read(iunit, io, "smoother_nml")

   call error_handler(E_MSG,'init_smoother','smoother_nml values are',' ',' ',' ')
   if (do_output()) write(logfileunit, nml=smoother_nml)
   if (do_output()) write(     *     , nml=smoother_nml)

   ! Allocate space for lag_size storage
   if(num_lags > 0) &
      allocate(lag_handle(num_lags), SmootherStateUnit(num_lags))

endif

end subroutine static_init_smoother

!-------------------------------------------------------------------------

subroutine init_smoother(ens_handle, POST_INF_COPY, POST_INF_SD_COPY)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in) :: POST_INF_COPY, POST_INF_SD_COPY

! static_init_smoother initializes module and read namelist
if ( .not. module_initialized ) call static_init_smoother()

! Initialize a null adaptive_inflate type since inflation is not done at lags
! NOTE: Using ens_handle here (not lag_handle) so it doesn't die for 0 lag choice
if(num_lags > 0) call adaptive_inflate_init(lag_inflate, 0, .false., .false., .false., &
   .true., 'no_lag_inflate', 'no_lag_inflate', 'no_lag_inflate', 1.0_r8, 0.0_r8,       &
   1.0_r8, 1.0_r8, 0.0_r8, ens_handle, POST_INF_COPY, POST_INF_SD_COPY, "Lag")

end subroutine init_smoother

!-------------------------------------------------------------------------

subroutine smoother_read_restart(ens_handle, ens_size, model_size, time1, init_time_days)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: model_size, ens_size
type(time_type),     intent(inout) :: time1
integer,             intent(in)    :: init_time_days

! Can either get restart by copying ensemble_ic into every lag or by reading from
! restart files. If reading from restart, may want to override the times as indicated
! by the filter namelist.
integer             :: i
character(len = 13) :: file_name

! must have called init_smoother() before using this routine
if ( .not. module_initialized ) then
   write(errstring, *)'cannot be called before init_smoother() called'
   call error_handler(E_ERR,'smoother_read_restart',errstring,source,revision,revdate)
endif

! Initialize the storage space for the lag ensembles
do i = 1, num_lags
   ! Only need 2 extra copies but have 6 for now to comply with filter's needs
   call init_ensemble_manager(lag_handle(i), ens_size + 6, model_size, 1)
end do

! If starting from restart, read these in
if(start_from_restart) then
   do i = 1, num_lags
      write(file_name, '("lag_", i5.5, "_ics")') i
      if(init_time_days >= 0) then
         call read_ensemble_restart(lag_handle(i), 1, ens_size, &
            start_from_restart, file_name, time1)
      else
         call read_ensemble_restart(lag_handle(i), 1, ens_size, &
            start_from_restart, file_name)
      endif
   end do
   ! All lags in restart are assumed to be current
   num_current_lags = num_lags
else
   ! If not starting from restart, just copy the filter ics to all lags
   do i = 1, num_lags
      call duplicate_ens(ens_handle, lag_handle(i), .true.)
   end do
   ! None of these are current to start with
   num_current_lags = 0
endif

end subroutine smoother_read_restart

!-------------------------------------------------------------------------

subroutine advance_smoother(ens_handle)

type(ensemble_type), intent(in) :: ens_handle

integer         :: smoother_tail, j

! must have called init_smoother() before using this routine
if ( .not. module_initialized ) then
   write(errstring, *)'cannot be called before init_smoother() called'
   call error_handler(E_ERR,'advance_smoother',errstring,source,revision,revdate)
endif

! Copy the newest state from the ensemble over the oldest state
! Storage is cyclic
smoother_tail = smoother_head - 1
if(smoother_tail == 0) smoother_tail = num_lags
call duplicate_ens(ens_handle, lag_handle(smoother_tail), .true.)

! Now make the head point to the latest
smoother_head = smoother_tail
 
! Update time of lag1 to be the updated filter time;
! Rest of times propagated back in the advance_smoother step
if(num_current_lags > 0) then
   do j = 1, ens_handle%my_num_copies
      lag_handle(smoother_head)%time(j) = ens_handle%time(j)
   end do
endif

end subroutine advance_smoother

!-------------------------------------------------------------------------

subroutine smoother_gen_copy_meta_data(num_output_state_members, output_inflation)

integer, intent(in) :: num_output_state_members
logical, intent(in) :: output_inflation

! Figures out the strings describing the output copies for the smoother state output files.
! These are the prior and posterior state output files. 

! The 4 is for ensemble mean and spread plus inflation mean and spread
character(len = 129) :: state_meta(num_output_state_members + 4)
character(len = 14)  :: file_name
character(len = 15)  :: meta_data_string
integer              :: i, ensemble_offset, num_state_copies

! must have called init_smoother() before using this routine
if ( .not. module_initialized ) then
   write(errstring, *)'cannot be called before init_smoother() called'
   call error_handler(E_ERR,'smoother_gen_copy_meta_data',errstring,source,revision,revdate)
endif

! Ensemble mean goes first 
num_state_copies = num_output_state_members + 2
smoother_state_mean_index = 1
state_meta(smoother_state_mean_index) = 'ensemble mean'

! Ensemble spread goes second
smoother_state_spread_index = 2
state_meta(smoother_state_spread_index) = 'ensemble spread'

! Check for too many output ensemble members
if(num_output_state_members > 10000) then
   write(errstring, *)'output metadata in smoother needs state ensemble size < 10000, not ', &
                      num_output_state_members
   call error_handler(E_ERR,'smoother_gen_copy_meta_data',errstring,source,revision,revdate)
endif

! Compute starting point for ensemble member output
ensemble_offset = 2

! Set up the metadata for the output state diagnostic files
do i = 1, num_output_state_members
   write(state_meta(i + ensemble_offset), '(a15, 1x, i6)') 'ensemble member', i
end do

! Netcdf output diagnostics for inflation; inefficient for single spatial
! NOTE: this should be changed to handle prior, posterior, both, or 
! neither inflation.  it must match what is written out in 
! filter_state_space_diagnostics() below.  space for an array state_space size
! long is created unless turned off by the namelist 'output_inflation' value,
! regardless of the actual inflation settings.
!!!if(do_single_ss_inflate(prior_inflate) .or. do_varying_ss_inflate(prior_inflate)) then
if (output_inflation) then
   num_state_copies = num_state_copies + 2
   state_meta(num_state_copies -1) = 'inflation mean'
   state_meta(num_state_copies) = 'inflation sd'
endif

! Set up diagnostic output for model state, if output is desired
do i = 1, num_lags
   ! Generate file name and metadata for lag i output
   write(file_name, '("Lag_", i5.5, "_Diag")') i
   write(meta_data_string, '("lag ", i5.5, " state")') i
   SmootherStateUnit(i) = init_diag_output(file_name, meta_data_string, &
      num_state_copies, state_meta)
end do

end subroutine smoother_gen_copy_meta_data

!-----------------------------------------------------------

subroutine smoother_write_restart(start_copy, end_copy)

integer,             intent(in)    :: start_copy, end_copy

! Write out each smoother lag to a restart file as requested
! Single versus multiple file status is selected by ensemble manager
! Names for file are set here

character(len = 17) :: file_name
integer             :: i, index

! must have called init_smoother() before using this routine
if ( .not. module_initialized ) then
   write(errstring, *)'cannot be called before init_smoother() called'
   call error_handler(E_ERR,'smoother_write_restart',errstring,source,revision,revdate)
endif

! Write out restart to each lag in turn
! Storage is cyclic with lag 1 pointed to by head
do i = 1, num_lags
   index = smoother_head + i - 1
   if(index > num_lags) index = index - num_lags
   write(file_name, '("lag_", i5.5, "_restart")') index
   call write_ensemble_restart(lag_handle(index), file_name, start_copy, end_copy)
end do

end subroutine smoother_write_restart

!-----------------------------------------------------------

subroutine smoother_assim(obs_ens_handle, seq, keys, ens_size, num_groups, obs_val_index, &
   ENS_MEAN_COPY, ENS_SD_COPY, PRIOR_INF_COPY, PRIOR_INF_SD_COPY, OBS_KEY_COPY, &
   OBS_GLOBAL_QC_COPY, OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END, OBS_PRIOR_VAR_START, &
   OBS_PRIOR_VAR_END)

type(ensemble_type),         intent(inout) :: obs_ens_handle
type(obs_sequence_type),     intent(in) :: seq
integer,                     intent(in)    :: keys(:)
integer,                     intent(in)    :: ens_size, num_groups, obs_val_index
integer,                     intent(in)    :: ENS_MEAN_COPY, ENS_SD_COPY, PRIOR_INF_COPY
integer,                     intent(in)    :: PRIOR_INF_SD_COPY
integer,                     intent(in)    :: OBS_KEY_COPY, OBS_GLOBAL_QC_COPY
integer,                     intent(in)    :: OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END
integer,                     intent(in)    :: OBS_PRIOR_VAR_START, OBS_PRIOR_VAR_END

integer :: i

! must have called init_smoother() before using this routine
if ( .not. module_initialized ) then
   write(errstring, *)'cannot be called before init_smoother() called'
   call error_handler(E_ERR,'smoother_assim',errstring,source,revision,revdate)
endif

do i = 1, num_lags
   call all_vars_to_all_copies(lag_handle(i))
   ! NEED A LAG INFLATE TYPE THAT DOES NO INFLATION FOR NOW
   call filter_assim(lag_handle(i), obs_ens_handle, seq, keys, ens_size, num_groups, &
      obs_val_index, lag_inflate, ENS_MEAN_COPY, ENS_SD_COPY, &
      PRIOR_INF_COPY, PRIOR_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, &
      OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END, OBS_PRIOR_VAR_START, &
      OBS_PRIOR_VAR_END, inflate_only = .false.)
end do

end subroutine smoother_assim

!-----------------------------------------------------------

function do_smoothing()

logical :: do_smoothing

! static_init_smoother initializes module and read namelist
! (which sets num_lags which if > 0 enables the whole process)

if ( .not. module_initialized ) call static_init_smoother()

do_smoothing = num_lags > 0

end function do_smoothing

!-----------------------------------------------------------

subroutine smoother_mean_spread(ens_size, ENS_MEAN_COPY, ENS_SD_COPY)

integer, intent(in) :: ens_size, ENS_MEAN_COPY, ENS_SD_COPY

integer :: i

! must have called init_smoother() before using this routine
if ( .not. module_initialized ) then
   write(errstring, *)'cannot be called before init_smoother() called'
   call error_handler(E_ERR,'smoother_mean_spread',errstring,source,revision,revdate)
endif

do i = 1, num_lags
   call compute_copy_mean_sd(lag_handle(i), 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)
end do

! Now back to var complete for diagnostics
do i = 1, num_lags
   call all_copies_to_all_vars(lag_handle(i))
end do

end subroutine smoother_mean_spread

!-----------------------------------------------------------

subroutine filter_state_space_diagnostics(out_unit, ens_handle, model_size, &
   num_output_state_members, output_state_mean_index, output_state_spread_index, &
   output_inflation, temp_ens, ENS_MEAN_COPY, ENS_SD_COPY, inflate, INF_COPY, INF_SD_COPY)

type(netcdf_file_type),      intent(inout) :: out_unit
type(ensemble_type),         intent(inout) :: ens_handle
integer,                     intent(in)    :: model_size, num_output_state_members
integer,                     intent(in)    :: output_state_mean_index, output_state_spread_index
! temp_ens is passed from above to avoid extra storage
real(r8),                    intent(out)   :: temp_ens(model_size)
type(adaptive_inflate_type), intent(in)    :: inflate
integer,                     intent(in)    :: ENS_MEAN_COPY, ENS_SD_COPY, INF_COPY, INF_SD_COPY
logical,                     intent(in)    :: output_inflation

type(time_type) :: temp_time
integer         :: ens_offset, j

! Assumes that mean and spread have already been computed

! must have called init_smoother() before using this routine
if ( .not. module_initialized ) then
   write(errstring, *)'cannot be called before init_smoother() called'
   call error_handler(E_ERR,'smoother_state_space_diagnostics',errstring,source,revision,revdate)
endif

! Output ensemble mean
call get_copy(0, ens_handle, ENS_MEAN_COPY, temp_ens)
if(my_task_id() == 0) call aoutput_diagnostics(out_unit, ens_handle%time(1), temp_ens, output_state_mean_index)

! Output ensemble spread
call get_copy(0, ens_handle, ENS_SD_COPY, temp_ens)
if(my_task_id() == 0) call aoutput_diagnostics(out_unit, ens_handle%time(1), temp_ens, output_state_spread_index)

! Compute the offset for copies of the ensemble
ens_offset = 2

! Output state diagnostics as required: NOTE: Prior has been inflated
do j = 1, num_output_state_members
   ! Get this state copy to PE 0; then output it
   call get_copy(0, ens_handle, j, temp_ens, temp_time)
   if(my_task_id() == 0) call aoutput_diagnostics( out_unit, temp_time, temp_ens, ens_offset + j)
end do

! Unless specifically asked not to, output inflation
if (output_inflation) then
   ! Output the spatially varying inflation if used
   if(do_varying_ss_inflate(inflate) .or. do_single_ss_inflate(inflate)) then
      call get_copy(0, ens_handle, INF_COPY, temp_ens)
   else
      ! Output inflation value as 1 if not in use (no inflation)
      temp_ens = 1.0_r8
   endif
   if(my_task_id() == 0) call aoutput_diagnostics(out_unit, ens_handle%time(1), temp_ens, &
      ens_offset + num_output_state_members + 1)


   if(do_varying_ss_inflate(inflate) .or. do_single_ss_inflate(inflate)) then
      call get_copy(0, ens_handle, INF_SD_COPY, temp_ens)
   else
      ! Output inflation sd as 0 if not in use
      temp_ens = 0.0_r8
   endif
   if(my_task_id() == 0) call aoutput_diagnostics(out_unit, ens_handle%time(1), temp_ens, &
      ens_offset + num_output_state_members + 2)
endif

end subroutine filter_state_space_diagnostics


!-----------------------------------------------------------

subroutine smoother_ss_diagnostics(model_size, num_output_state_members, output_inflation, &
   temp_ens, ENS_MEAN_COPY, ENS_SD_COPY, POST_INF_COPY, POST_INF_SD_COPY)

   integer,  intent(in)  :: model_size, num_output_state_members
   logical,  intent(in)  :: output_inflation
   real(r8), intent(out) :: temp_ens(model_size)
   integer,  intent(in)  :: ENS_MEAN_COPY, ENS_SD_COPY, POST_INF_COPY, POST_INF_SD_COPY

integer :: smoother_index, i

! must have called init_smoother() before using this routine
if ( .not. module_initialized ) then
   write(errstring, *)'cannot be called before init_smoother() called'
   call error_handler(E_ERR,'smoother_ss_diagnostics',errstring,source,revision,revdate)
endif

do i = 1, num_current_lags
   smoother_index = smoother_head + i - 1
   if(smoother_index > num_lags) smoother_index = smoother_index - num_lags
   call filter_state_space_diagnostics(SmootherStateUnit(i), lag_handle(smoother_index), &
      model_size, num_output_state_members, &
      smoother_state_mean_index, smoother_state_spread_index, output_inflation, temp_ens, &
      ENS_MEAN_COPY, ENS_SD_COPY, lag_inflate, POST_INF_COPY, POST_INF_SD_COPY)
end do


end subroutine smoother_ss_diagnostics

!-----------------------------------------------------------

subroutine smoother_end()

! Free up smoother allocated storage

! Probably need to free up ensemble manager storage too

deallocate(lag_handle, SmootherStateUnit)

end subroutine smoother_end

!-----------------------------------------------------------

subroutine smoother_inc_lags()

! Increment the number of lags that are current for the smoother storage

if(num_current_lags < num_lags) num_current_lags = num_current_lags + 1

end subroutine smoother_inc_lags
!-----------------------------------------------------------

end module smoother_mod
