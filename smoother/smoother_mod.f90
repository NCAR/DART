! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module smoother_mod 

! Tools for turning the filter into a fixed lag smoother for the full state vector.

use      types_mod,       only : r8, metadatalength
use  mpi_utilities_mod,   only : my_task_id
use  utilities_mod,       only : file_exist, get_unit, check_namelist_read, do_output,  &
                                 find_namelist_in_file, register_module, error_handler, &
                                 E_ERR, E_MSG, nmlfileunit, logfileunit, timestamp,     &
                                 do_nml_file, do_nml_term
use ensemble_manager_mod, only : ensemble_type, init_ensemble_manager, read_ensemble_restart, &
                                 write_ensemble_restart, all_vars_to_all_copies,              &
                                 duplicate_ens, compute_copy_mean, compute_copy_mean_sd,      &
                                 all_copies_to_all_vars, get_copy, map_task_to_pe
use time_manager_mod,     only : time_type, operator(==), print_time
use assim_model_mod,      only : static_init_assim_model, get_model_size,                    &
                                 netcdf_file_type, init_diag_output, finalize_diag_output,   &
                                 aoutput_diagnostics
use assim_tools_mod,      only : filter_assim, get_missing_ok_status
use obs_sequence_mod,     only : obs_sequence_type
use adaptive_inflate_mod, only : adaptive_inflate_type, adaptive_inflate_init, &
                                 do_varying_ss_inflate, do_single_ss_inflate
use smoother_pnetcdf_mod ! this could be a null module if you don't want to use pnetcdf

implicit none
private

public :: smoother_read_restart, advance_smoother,                     &
   smoother_gen_copy_meta_data, smoother_write_restart, init_smoother, &
   do_smoothing, smoother_mean_spread, smoother_assim,                 &
   filter_state_space_diagnostics, smoother_ss_diagnostics,            &
   smoother_end, set_smoother_trace, query_pnetcdf


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical :: module_initialized = .false.
integer :: print_trace_details = 0
integer :: print_timestamps    = 0

character(len = 129) :: errstring

! State for the smoother
integer :: smoother_head = 1, num_current_lags = 0
integer :: smoother_state_mean_index, smoother_state_spread_index
type(ensemble_type), allocatable :: lag_handle(:)
type(netcdf_file_type), allocatable :: SmootherStateUnit(:)
type(adaptive_inflate_type)         :: lag_inflate

!============================================================================

!---- namelist with default values

integer  :: num_lags           = 0
logical  :: start_from_restart = .false.
logical  :: output_restart     = .false.

! These will be prepended with Lag_NNNNN_ (and optionally appended with
! ens number depending on namelist settings).
character(len = 129) :: restart_in_file_name  = 'ics', &
                        restart_out_file_name = 'restart'

namelist /smoother_nml/ num_lags, start_from_restart, &
                        output_restart, restart_in_file_name, &
                        restart_out_file_name

interface filter_state_space_diagnostics
   module procedure filter_state_space_diagnostics_complete, filter_state_space_diagnostics_parallel
end interface

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

   if (do_nml_file()) write(nmlfileunit, nml=smoother_nml)
   if (do_nml_term()) write(     *     , nml=smoother_nml)

   ! Allocate space for lag_size storage
   if(num_lags > 0) &
      allocate(lag_handle(num_lags), SmootherStateUnit(num_lags))

endif

end subroutine static_init_smoother

!-------------------------------------------------------------------------

subroutine init_smoother(ens_handle, POST_INF_COPY, POST_INF_SD_COPY)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in) :: POST_INF_COPY, POST_INF_SD_COPY

logical :: allow_missing

! static_init_smoother initializes module and read namelist
if ( .not. module_initialized ) call static_init_smoother()

! find out if it is ok to have missing values in the state vector
allow_missing = get_missing_ok_status()

! Initialize a null adaptive_inflate type since inflation is not done at lags
! NOTE: Using ens_handle here (not lag_handle) so it doesn't die for 0 lag choice
if(num_lags > 0) call adaptive_inflate_init(lag_inflate, 0, .false., .false., .false., &
   .true., 'no_lag_inflate', 'no_lag_inflate', 'no_lag_inflate', 1.0_r8, 0.0_r8,       &
   1.0_r8, 1.0_r8, 0.0_r8, ens_handle, POST_INF_COPY, POST_INF_SD_COPY, allow_missing, "Lag")

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
integer              :: i, j, smoother_index
character(len = 256) :: file_name, temp_name

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

! assume no lags are current, and increment this as the files are found.
! this allows a model which has only advanced N times to have > N lags.
num_current_lags = 0

! If starting from restart, read these in
if(start_from_restart) then
   READ_LAGS: do i = 1, num_lags
      smoother_index = next_index(i)
      write(file_name, '("Lag_", I5.5, "_", A)') i, trim(restart_in_file_name)
      write(temp_name, '(A, ".", I4.4)') trim(file_name), 1
      if (file_exist(file_name) .or. file_exist(temp_name)) then
         if(init_time_days >= 0) then
            call read_ensemble_restart(lag_handle(smoother_index), 1, ens_size, &
               start_from_restart, file_name, time1)
         else
            call read_ensemble_restart(lag_handle(smoother_index), 1, ens_size, &
               start_from_restart, file_name)
         endif
         num_current_lags = num_current_lags + 1
   !write(errstring, '(A,I4,A,I4)') 'reading restart file ', i, ' into cycle number', smoother_index
   !call error_handler(E_MSG, 'smoother_read_restart', errstring)
      else
         ! lag ic file does not exist yet, duplicate the filter ics 
         ! for the rest of the lags and break out of the i loop
         do j = i, num_lags
            smoother_index = next_index(j)
            call duplicate_ens(ens_handle, lag_handle(smoother_index), .true.)
   !write(errstring, '(A,I4,A,I4)') 'filling restart data ', i, ' into cycle number', smoother_index
   !call error_handler(E_MSG, 'smoother_read_restart', errstring)
         end do
         exit READ_LAGS
      endif
   end do READ_LAGS

   !write(errstring, '(i4, A)') num_current_lags, ' smoother restart files processed' 
   !call error_handler(E_MSG,'smoother_read_restart',errstring)

else
   ! If not starting from restart, just copy the filter ics to all lags
   do i = 1, num_lags
      call duplicate_ens(ens_handle, lag_handle(i), .true.)
   end do
   !write(errstring, '(A,I4,A)') 'no restart data found, filling all', num_lags, ' lags'
   !call error_handler(E_MSG, 'smoother_read_restart', errstring)
endif

!write(errstring, '(A,I4,A)') 'head =', smoother_head
!call error_handler(E_MSG, 'smoother_read_restart', errstring)

end subroutine smoother_read_restart

!-------------------------------------------------------------------------

subroutine advance_smoother(ens_handle)

type(ensemble_type), intent(in) :: ens_handle

integer         :: smoother_tail

! must have called init_smoother() before using this routine
if ( .not. module_initialized ) then
   write(errstring, *)'cannot be called before init_smoother() called'
   call error_handler(E_ERR,'advance_smoother',errstring,source,revision,revdate)
endif

!call error_handler(E_MSG, 'advance_smoother', 'start of routine')
!if (num_lags /= num_current_lags) then
!   write(errstring, '(A,I4)') 'num_current_lags starts =', num_current_lags
!   call error_handler(E_MSG, 'advance_smoother', errstring)
!endif

!write(errstring, '(A,I4)') 'apparently time has advanced, head starts =', smoother_head
!call error_handler(E_MSG, 'advance_smoother', errstring)

! Copy the newest state from the ensemble over the oldest state
! Storage is cyclic
smoother_tail = smoother_head - 1
if(smoother_tail == 0) smoother_tail = num_lags
call duplicate_ens(ens_handle, lag_handle(smoother_tail), .true.)
!write(errstring, '(A,I4)') 'copied current ens data into tail =', smoother_tail
!call error_handler(E_MSG, 'advance_smoother', errstring)


! Make the head point to the most recent data copy
smoother_head = smoother_tail

! Add one to smoother ens if necessary
call smoother_inc_lags()

!write(errstring, '(A,I4,A,I4)') 'head now =', smoother_head, ' num_current_lags =', num_current_lags
!call error_handler(E_MSG, 'advance_smoother', errstring)

! debug only
!call print_time(lag_handle(smoother_head)%time(1), ' advance_smoother: head time now', logfileunit)
!call print_time(lag_handle(smoother_head)%time(1), ' advance_smoother: head time now')
!smoother_tail = smoother_head - 1
!if(smoother_tail == 0) smoother_tail = num_lags
!call print_time(lag_handle(smoother_tail)%time(1), ' advance_smoother: tail time now', logfileunit)
!call print_time(lag_handle(smoother_tail)%time(1), ' advance_smoother: tail time now')

end subroutine advance_smoother

!-------------------------------------------------------------------------

subroutine smoother_gen_copy_meta_data(num_output_state_members, output_inflation)

integer, intent(in) :: num_output_state_members
logical, intent(in) :: output_inflation

! Figures out the strings describing the output copies for the smoother state output files.
! These are the prior and posterior state output files. 

! The 4 is for ensemble mean and spread plus inflation mean and spread
character(len = metadatalength) :: state_meta(num_output_state_members + 4)
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

character(len = 256) :: file_name
integer              :: i, smoother_index

! must have called init_smoother() before using this routine
if ( .not. module_initialized ) then
   write(errstring, *)'cannot be called before init_smoother() called'
   call error_handler(E_ERR,'smoother_write_restart',errstring,source,revision,revdate)
endif

! namelist controlled
if (.not. output_restart) return

! Write out restart to each lag in turn
! Storage is cyclic with oldest lag pointed to by head, and 
! head + 1 the most recent lag
do i = 1, num_current_lags
   smoother_index = next_index(i)
   write(file_name, '("Lag_", I5.5, "_", A)') i, trim(restart_out_file_name)
   call write_ensemble_restart(lag_handle(smoother_index), file_name, start_copy, end_copy)
   !write(errstring, '(A,I4,A,I4)') 'writing restart file ', i, ' from cycle number', smoother_index
   !call error_handler(E_MSG, 'smoother_write_restart', errstring)
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

integer :: smoother_index, i

! must have called init_smoother() before using this routine
if ( .not. module_initialized ) then
   write(errstring, *)'cannot be called before init_smoother() called'
   call error_handler(E_ERR,'smoother_assim',errstring,source,revision,revdate)
endif

do i = 1, num_current_lags
   smoother_index = next_index(i)
   call all_vars_to_all_copies(lag_handle(smoother_index))

   !write(errstring, '(A,I4,A,I4)') 'starting assimilate pass for lag', i, &
   !                                ' data, cycle index', smoother_index
   !call error_handler(E_MSG,'smoother_assim',errstring)

   write(errstring, '(A,I4,A)') 'Starting reassimilate pass for lag', i, ' data'
   if (print_trace_details >= 0) call error_handler(E_MSG,'smoother_assim:',errstring)
   if (print_trace_details >= 1) then
      call print_time(lag_handle(smoother_index)%time(1), &
         ' smoother_assim: Time of lagged data is: ', logfileunit)
      call print_time(lag_handle(smoother_index)%time(1), &
         ' smoother_assim: Time of lagged data is: ')
   endif

   ! NEED A LAG INFLATE TYPE THAT DOES NO INFLATION FOR NOW
   call filter_assim(lag_handle(smoother_index), obs_ens_handle, &
      seq, keys, ens_size, num_groups, &
      obs_val_index, lag_inflate, ENS_MEAN_COPY, ENS_SD_COPY, &
      PRIOR_INF_COPY, PRIOR_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, &
      OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END, OBS_PRIOR_VAR_START, &
      OBS_PRIOR_VAR_END, inflate_only = .false.)

   !write(errstring, '(A,I8,A,I4,A)') 'finished assimilating ', size(keys), &
   !                                  ' observations for lag', i, ' data'
   !call error_handler(E_MSG,'smoother_assim', errstring)
   !write(errstring, '(A,I4)') ' smoother_assim: smoother_index =', smoother_index
   !call error_handler(E_MSG,'smoother_assim', errstring)
   !call print_time(lag_handle(smoother_index)%time(1), errstring)
   write(errstring, '(A,I8,A,I4,A)') 'finished assimilating ', size(keys), &
                                     ' observations for lag', i, ' data'
   if (do_output() .and. print_timestamps > 0) call timestamp(errstring, pos='debug')
end do

!call error_handler(E_MSG,'smoother_assim', 'end of routine, returning')

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

integer :: smoother_index, i

! must have called init_smoother() before using this routine
if ( .not. module_initialized ) then
   write(errstring, *)'cannot be called before init_smoother() called'
   call error_handler(E_ERR,'smoother_mean_spread',errstring,source,revision,revdate)
endif

! Must be called when smoother handles are copy complete.
! Leaves the data var complete before it returns.
do i = 1, num_current_lags
   smoother_index = next_index(i)

   call compute_copy_mean_sd(lag_handle(smoother_index), 1, ens_size, ENS_MEAN_COPY, ENS_SD_COPY)
   call all_copies_to_all_vars(lag_handle(smoother_index))
end do

end subroutine smoother_mean_spread

!-----------------------------------------------------------

subroutine filter_state_space_diagnostics_complete(curr_ens_time, out_unit, ens_handle, model_size, &
            num_output_state_members, output_state_mean_index, output_state_spread_index, &
           output_inflation, temp_ens, ENS_MEAN_COPY, ENS_SD_COPY, inflate, INF_COPY, INF_SD_COPY)

type(time_type),             intent(in)    :: curr_ens_time
type(netcdf_file_type),      intent(inout) :: out_unit
type(ensemble_type),         intent(inout) :: ens_handle
integer,                     intent(in)    :: model_size, num_output_state_members
integer,                     intent(in)    :: output_state_mean_index, output_state_spread_index
! temp_ens is passed from above to avoid extra storage
real(r8),                    intent(out)   :: temp_ens(:)
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
call get_copy(map_task_to_pe(ens_handle, 0), ens_handle, ENS_MEAN_COPY, temp_ens)
if(my_task_id() == 0) call aoutput_diagnostics(out_unit, curr_ens_time, temp_ens,  &
   output_state_mean_index)

! Output ensemble spread
call get_copy(map_task_to_pe(ens_handle, 0), ens_handle, ENS_SD_COPY, temp_ens) 
if(my_task_id() == 0) call aoutput_diagnostics(out_unit, curr_ens_time, temp_ens, &
   output_state_spread_index)

! Compute the offset for copies of the ensemble
ens_offset = 2

! Output state diagnostics as required: NOTE: Prior has been inflated
do j = 1, num_output_state_members
   ! Get this state copy to task 0; then output it
   call get_copy(map_task_to_pe(ens_handle, 0), ens_handle, j, temp_ens, temp_time)
   if(my_task_id() == 0) call aoutput_diagnostics( out_unit, temp_time, temp_ens, ens_offset + j)
end do

! Unless specifically asked not to, output inflation
if (output_inflation) then
   ! Output the spatially varying inflation if used
   if(do_varying_ss_inflate(inflate) .or. do_single_ss_inflate(inflate)) then
      call get_copy(map_task_to_pe(ens_handle, 0), ens_handle, INF_COPY, temp_ens)
   else
      ! Output inflation value as 1 if not in use (no inflation)
      temp_ens = 1.0_r8
   endif

   if(my_task_id() == 0) call aoutput_diagnostics(out_unit,  curr_ens_time, temp_ens, &
     ens_offset + num_output_state_members + 1)  


   if(do_varying_ss_inflate(inflate) .or. do_single_ss_inflate(inflate)) then
      call get_copy(map_task_to_pe(ens_handle, 0), ens_handle, INF_SD_COPY, temp_ens)
   else
      ! Output inflation sd as 0 if not in use
      temp_ens = 0.0_r8
   endif

   if(my_task_id() == 0) call aoutput_diagnostics(out_unit, curr_ens_time, temp_ens, &
      ens_offset + num_output_state_members + 2) 

endif

end subroutine filter_state_space_diagnostics_complete

!-----------------------------------------------------------

subroutine smoother_ss_diagnostics(model_size, num_output_state_members, output_inflation, &
   temp_ens, ENS_MEAN_COPY, ENS_SD_COPY, POST_INF_COPY, POST_INF_SD_COPY)

use mpi

integer,         intent(in)  :: model_size, num_output_state_members
logical,         intent(in)  :: output_inflation
real(r8),        intent(out) :: temp_ens(model_size)
integer,         intent(in)  :: ENS_MEAN_COPY, ENS_SD_COPY, POST_INF_COPY, POST_INF_SD_COPY

integer :: smoother_index, i

! timing variables
double precision :: start_at_time

start_at_time = MPI_WTIME()

! must have called init_smoother() before using this routine
if ( .not. module_initialized ) then
   write(errstring, *)'cannot be called before init_smoother() called'
   call error_handler(E_ERR,'smoother_ss_diagnostics',errstring,source,revision,revdate)
endif

do i = 1, num_current_lags
   smoother_index = next_index(i)

   ! FIXME: is this still needed?
   call all_copies_to_all_vars(lag_handle(smoother_index))
   ! only ensemble copies have the time
   call filter_state_space_diagnostics(lag_handle(smoother_index)%time(1), SmootherStateUnit(i), &
      lag_handle(smoother_index), model_size, num_output_state_members, &
      smoother_state_mean_index, smoother_state_spread_index, output_inflation, temp_ens, &
      ENS_MEAN_COPY, ENS_SD_COPY, lag_inflate, POST_INF_COPY, POST_INF_SD_COPY)
end do

if (my_task_id() == 0) print*, 'serial diagnostic time :', MPI_WTIME() - start_at_time

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

! must have called init_smoother() before using this routine
if ( .not. module_initialized ) then
   write(errstring, *)'cannot be called before init_smoother() called'
   call error_handler(E_ERR,'smoother_inc_lags',errstring,source,revision,revdate)
endif


if(num_current_lags < num_lags) num_current_lags = num_current_lags + 1

!write(errstring, *)'smoother_inc_lags called, num_current_lags =', num_current_lags
!call error_handler(E_MSG,'smoother_inc_lags',errstring)


end subroutine smoother_inc_lags

!-----------------------------------------------------------

function next_index(i)
 integer, intent(in) :: i
 integer             :: next_index
! the index numbers for the collection of lag_handles is
! a circular list.  increment and do the wrap if needed.

if ( .not. module_initialized ) then
   write(errstring, *)'cannot be called before init_smoother() called'
   call error_handler(E_ERR,'smoother_inc_lags',errstring,source,revision,revdate)
endif


next_index = smoother_head + i - 1
if(next_index > num_lags) next_index = next_index - num_lags

!write(errstring, *)'next_index called, index =', next_index
!call error_handler(E_MSG,'next_index',errstring)

end function next_index

!-----------------------------------------------------------

subroutine set_smoother_trace(execution_level, timestamp_level)
 integer, intent(in) :: execution_level
 integer, intent(in) :: timestamp_level

! set module local vars from the calling code to indicate how much
! output we should generate from this code.  execution level is
! intended to make it easier to figure out where in the code a crash
! is happening; timestamp level is intended to help with gross levels
! of overall performance profiling.  eventually, a level of 1 will
! print out only basic info; level 2 will be more detailed.  
! (right now, only > 0 prints anything and it doesn't matter how 
! large the value is.)
 

print_trace_details = execution_level
print_timestamps    = timestamp_level

end subroutine set_smoother_trace

!-----------------------------------------------------------

end module smoother_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
