MODULE smoother_mod
===================

.. attention::

  The DART smoother works with versions of DART *before* Manhattan (9.x.x) and has yet to be updated. If you are interested in
  using ``smoother_mod`` with more recent versions of DART, contact DAReS staff to assess the feasibility of an update.
  Until that time, you should consider this documentation as out-of-date.

Overview
--------

| Implements a fixed lag ensemble smoother as part of the filter. For now, this is done inefficiently with a separate
  call to ``assim_tools_mod:filter_assim()`` for each lag.
| To enable the smoother, set the number of lags (num_lags) to something larger than 0 in the ``smoother_nml`` section
  of your ``input.nml`` file and run ``filter`` as before.

.. container:: routine

   ::

      &smoother_nml
         num_lags              = 10,
         start_from_restart    = .false.,
         output_restart        = .true.,
         restart_in_file_name  = "ics",
         restart_out_file_name = "restart"  /

| In the low order models, 10 is a plausible number.
| In addition to generating ``preassim.nc`` and ``analysis.nc`` files, files of the form ``Lag_NNNNN_Diag.nc`` will be
  generated. Each of these has N fewer timesteps than the lag=0 run, starting at the same time but ending N timesteps
  sooner. The ``obs_seq.final`` file and the ``preassim.nc`` and ``analysis.nc`` files will be the same as the
  non-lagged version; the new output will be in each of the ``Lag_NNNNN_Diag.nc`` files.

Example
-------

| If you have a ``true_state.nc`` file and want to use the ``plot_total_err`` matlab function to plot the error, you
  must do the following steps to generate analogs of lagged ``true_state.nc`` files to use as a comparison. (The logic
  is not currently implemented in the matlab scripts to be able to compare netCDF files with unequal time coordinates.)
| Make N separate versions of the true_state.nc with the last N timesteps removed. Using the netCDF NCO operator program
  'ncks' is one way. If the true_state.nc file has 1000 time steps, then this command removes the last one:

.. container:: unix

   ncks -d time,0,998 true_state.nc True_Lag01.nc

Note that the first time is at index 0, so the last timestep is index 999 in the full file, and 998 in the truncated
file. Repeat this step for all N lags. Here are NCO commands to generate 10 truth files for num_lags = 10, 1000 time
steps in true_state.nc:

.. container:: unix

   ncks -d time,0,998 true_state.nc True_Lag01.nc
   ncks -d time,0,997 true_state.nc True_Lag02.nc
   ncks -d time,0,996 true_state.nc True_Lag03.nc
   ncks -d time,0,995 true_state.nc True_Lag04.nc
   ncks -d time,0,994 true_state.nc True_Lag05.nc
   ncks -d time,0,993 true_state.nc True_Lag06.nc
   ncks -d time,0,992 true_state.nc True_Lag07.nc
   ncks -d time,0,991 true_state.nc True_Lag08.nc
   ncks -d time,0,990 true_state.nc True_Lag09.nc
   ncks -d time,0,989 true_state.nc True_Lag10.nc

Here is an example matlab session which plots the lag=0 results and then odd numbered lags from 1 to 9. It uses the
``plot_total_err`` function from the $DART/matlab directory:

::

   datadir    = '.';
   truth_file = fullfile(datadir,'true_state.nc');
   diagn_file = fullfile(datadir,'preassim.nc');
   plot_total_err
   reply = input('original data.  hit enter to continue ');

   truth_file = fullfile(datadir,'True_Lag01.nc');
   diagn_file = fullfile(datadir,'Lag_00001_Diag.nc');
   plot_total_err
   reply = input('Lag 01.  hit enter to continue ');

   truth_file = fullfile(datadir,'True_Lag03.nc');
   diagn_file = fullfile(datadir,'Lag_00003_Diag.nc');
   plot_total_err
   reply = input('Lag 03.  hit enter to continue ');

   truth_file = fullfile(datadir,'True_Lag05.nc');
   diagn_file = fullfile(datadir,'Lag_00005_Diag.nc');
   plot_total_err
   reply = input('Lag 05.  hit enter to continue ');

   truth_file = fullfile(datadir,'True_Lag07.nc');
   diagn_file = fullfile(datadir,'Lag_00007_Diag.nc');
   plot_total_err
   reply = input('Lag 07.  hit enter to continue ');

   truth_file = fullfile(datadir,'True_Lag09.nc');
   diagn_file = fullfile(datadir,'Lag_00009_Diag.nc');
   plot_total_err
   reply = input('Lag 09.  hit enter to continue ');

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &smoother_nml
      num_lags              = 0,
      start_from_restart    = .false.,
      output_restart        = .false.,
      restart_in_file_name  = 'ics',
      restart_out_file_name = 'restart'  
   /

| 

.. container::

   +-----------------------+--------------------+-----------------------------------------------------------------------+
   | Item                  | Type               | Description                                                           |
   +=======================+====================+=======================================================================+
   | num_lags              | integer            | Number of smoother lags; < 1 means no smoother.                       |
   +-----------------------+--------------------+-----------------------------------------------------------------------+
   | start_from_restart    | logical            | True if smoother states are to come from restart file(s). False if    |
   |                       |                    | they are to be spun up from scratch.                                  |
   +-----------------------+--------------------+-----------------------------------------------------------------------+
   | output_restart        | logical            | True if restart file(s) are to be written, else false.                |
   +-----------------------+--------------------+-----------------------------------------------------------------------+
   | restart_in_file_name  | character(len=129) | String used to construct the file name from which to read restart     |
   |                       |                    | data. ``Lag_NNNNN_`` will be prepended to the specified value to      |
   |                       |                    | create the actual filename. If each ensemble is to be read from a     |
   |                       |                    | separate file, the .NNNN ensemble number will also be appended. e.g.  |
   |                       |                    | specifying 'ics' here results in 'Lag_00001_ics' if all ensemble      |
   |                       |                    | members are read from a single file, 'Lag_00001_ics.0001',            |
   |                       |                    | 'Lag_00001_ics.0002', etc for multiples.                              |
   +-----------------------+--------------------+-----------------------------------------------------------------------+
   | restart_out_file_name | character(len=129) | String used to construct the file name to which to write restart      |
   |                       |                    | data. ``Lag_NNNNN_`` will be prepended to the specified value to      |
   |                       |                    | create the actual filename. If each ensemble is to be written to a    |
   |                       |                    | separate file, the .NNNN ensemble number will also be appended. e.g.  |
   |                       |                    | specifying 'restart' here results in 'Lag_00001_restart' if all       |
   |                       |                    | ensemble members are written to a single file,                        |
   |                       |                    | 'Lag_00001_restart.0001', 'Lag_00001_restart.0002', etc for           |
   |                       |                    | multiples.                                                            |
   +-----------------------+--------------------+-----------------------------------------------------------------------+

| 

Other modules used
------------------

::

   types_mod
   mpi_utilities_mod
   utilities_mod
   ensemble_manager_mod
   time_manager_mod
   assim_model_mod
   assim_tools_mod
   obs_sequence_mod
   adaptive_inflate_mod

Public interfaces
-----------------

========================== ==============================
*use smoother_mod, only :* smoother_read_restart
\                          advance_smoother
\                          smoother_gen_copy_meta_data
\                          smoother_write_restart
\                          init_smoother
\                          do_smoothing
\                          smoother_mean_spread
\                          smoother_assim
\                          filter_state_space_diagnostics
\                          smoother_ss_diagnostics
\                          smoother_end
========================== ==============================

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *call smoother_read_restart(ens_handle, ens_size, model_size, time1, init_time_days)*
   ::

      type(ensemble_type), intent(inout) :: ens_handle
      integer, intent(in)                :: ens_size
      integer, intent(in)                :: model_size
      type(time_type), intent(inout)     :: time1
      integer, intent(in)                :: init_time_days

.. container:: indent1

   Reads in ensemble of states for all lag estimates from a restart file.

   ================== =========================================================================================
   ``ens_handle``     Handle of ensemble manager structure of single state; copied into all lags for startup.
   ``ens_size``       Size of the ensemble.
   ``model_size``     Size of the model state vector.
   ``time1``          Overwrite the time in the restart file with this value if init_time_days is non-negative.
   ``init_time_days`` If non-negative, use time1 instead of time in restart file.
   ================== =========================================================================================

| 

.. container:: routine

   *call advance_smoother(ens_handle)*
   ::

      type(ensemble_type), intent(in) :: ens_handle

.. container:: indent1

   Advances smoother state estimates at all lags forward in time. This entails copying the most recent smoother state,
   contained in ens_handle, into the lag 1 smoother state and pushing back all other lags by 1 (i.e. lag 1 becomes lag
   2, etc.).

   ============== ================================================
   ``ens_handle`` Ensemble handle with most recent filtered state.
   ============== ================================================

| 

.. container:: routine

   *call smoother_gen_copy_meta_data(num_output_state_members, output_inflation)*
   ::

      integer, intent(in) :: num_output_state_members
      logical, intent(in) :: output_inflation

.. container:: indent1

   Initializes the metadata required for the smoother state space diagnostic files.

   ============================ ==========================================================================================
   ``num_output_state_members`` Number of copies of smoother state vector that should be in state space diagnostic output.
   ``output_inflation``         True if smoother state space output should include inflation values.
   ============================ ==========================================================================================

| 

.. container:: routine

   *call smoother_write_restart(start_copy, end_copy)*
   ::

      integer, intent(in) :: start_copy
      integer, intent(in) :: end_copy

.. container:: indent1

   Outputs restart files for all lags of smoother state. Integer arguments specify the start and end global indices of a
   continguous set of copies that contain the ensemble members.

   ============== ===================================================================================
   ``start_copy`` Global index of ensemble copy that starts the actual ensemble members for smoother.
   ``end_copy``   Global index of ensemble copy that ends the actual ensemble members for smoother.
   ============== ===================================================================================

| 

.. container:: routine

   *call init_smoother(ens_handle, POST_INF_COPY, POST_INF_SD_COPY)*
   ::

      type(ensemble_type), intent(inout) :: ens_handle
      integer, intent(in)                :: POST_INF_COPY
      integer, intent(in)                :: POST_INF_SD_COPY

.. container:: indent1

   Initializes the storage needed for a smoother. Also initializes an adaptive inflation type that does NO inflation
   (not currently supported for smoothers).

   ==================== ==========================================================================================
   ``ens_handle``       An ensemble handle for the filter that contains information about ensemble and model size.
   ``POST_INF_COPY``    Global index of ensemble copy that holds posterior state space inflation values.
   ``POST_INF_SD_COPY`` Global index of ensemble copy that holds posterior inflation standard deviation values.
   ==================== ==========================================================================================

| 

.. container:: routine

   *var = do_smoothing()*
   ::

      logical, intent(out) :: do_smoothing

.. container:: indent1

   Returns true if smoothing is to be done, else false.

   ================ ========================================
   ``do_smoothing`` Returns true if smoothing is to be done.
   ================ ========================================

| 

.. container:: routine

   *call smoother_mean_spread(ens_size,ENS_MEAN_COPY,ENS_SD_COPY, output_state_ens_mean,output_state_ens_spread)*
   ::

      integer, intent(in) :: ens_size
      integer, intent(in) :: ENS_MEAN_COPY
      integer, intent(in) :: ENS_SD_COPY
      logical, intent(in) :: output_state_ens_mean
      logical, intent(in) :: output_state_ens_spread

.. container:: indent1

   Computes the ensemble mean (and spread if required) of all state variables for all lagged ensembles. Spread is only
   computed if it is required for output.

   =========================== ===================================================================
   ``ens_size``                Size of ensemble.
   ``ENS_MEAN_COPY``           Global index of copy that stores ensemble mean.
   ``ENS_SD_COPY``             Global index of copy that stores ensemble spread.
   ``output_state_ens_mean``   True if the ensemble mean is to be output to state diagnostic file.
   ``output_state_ens_spread`` True if ensemble spread is to be output to state diagnostic file.
   =========================== ===================================================================

| 

.. container:: routine

   *call smoother_assim(obs_ens_handle, seq, keys, ens_size, num_groups, obs_val_index, ENS_MEAN_COPY, ENS_SD_COPY,
   PRIOR_INF_COPY, PRIOR_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY, OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END,
   OBS_PRIOR_VAR_START, OBS_PRIOR_VAR_END)*
   ::

      type(ensemble_type), intent(inout)  :: obs_ens_handle
      type(obs_sequence_type), intent(in) :: seq
      integer, dimension(:), intent(in)   :: keys
      integer, intent(in)                 :: ens_size
      integer, intent(in)                 :: num_groups
      integer, intent(in)                 :: obs_val_index
      integer, intent(in)                 :: ENS_MEAN_COPY
      integer, intent(in)                 :: ENS_SD_COPY
      integer, intent(in)                 :: PRIOR_INF_COPY
      integer, intent(in)                 :: PRIOR_INF_SD_COPY
      integer, intent(in)                 :: OBS_KEY_COPY
      integer, intent(in)                 :: OBS_GLOBAL_QC_COPY
      integer, intent(in)                 :: OBS_PRIOR_MEAN_START
      integer, intent(in)                 :: OBS_PRIOR_MEAN_END
      integer, intent(in)                 :: OBS_PRIOR_VAR_START
      integer, intent(in)                 :: OBS_PRIOR_VAR_END

.. container:: indent1

   Does assimilation of a set of observations for each smoother lag.

   +--------------------------+------------------------------------------------------------------------------------------+
   | ``obs_ens_handle``       | Handle for ensemble manager holding prior estimates of observations.                     |
   +--------------------------+------------------------------------------------------------------------------------------+
   | ``seq``                  | Observation sequence being assimilated.                                                  |
   +--------------------------+------------------------------------------------------------------------------------------+
   | ``keys``                 | A one dimensional array containing indices in seq of observations to as similate at      |
   |                          | current time.                                                                            |
   +--------------------------+------------------------------------------------------------------------------------------+
   | ``ens_size``             | Ensemble size.                                                                           |
   +--------------------------+------------------------------------------------------------------------------------------+
   | ``num_groups``           | Number of groups in filter.                                                              |
   +--------------------------+------------------------------------------------------------------------------------------+
   | ``obs_val_index``        | Integer index of copy of data in seq that contains the observed value from instruments.  |
   +--------------------------+------------------------------------------------------------------------------------------+
   | ``ENS_MEAN_COPY``        | Global index in smoother's state ensemble that holds ensemble mean.                      |
   +--------------------------+------------------------------------------------------------------------------------------+
   | ``ENS_SD_COPY``          | Global index in smoother's state ensemble that holds ensemble standard deviation.        |
   +--------------------------+------------------------------------------------------------------------------------------+
   | ``PRIOR_INF_COPY``       | Global index in obs_ens_handle that holds inflation values (not used for smoother).      |
   +--------------------------+------------------------------------------------------------------------------------------+
   | ``PRIOR_INF_SD_COPY``    | Global index in obs_ens_handle that holds inflation sd values (not used for smoother).   |
   +--------------------------+------------------------------------------------------------------------------------------+
   | ``OBS_KEY_COPY``         | Global index in obs_ens_handle that holds the key for the observation.                   |
   +--------------------------+------------------------------------------------------------------------------------------+
   | ``OBS_GLOBAL_QC_COPY``   | Global index in obs_ens_handle that holds the quality control value.                     |
   +--------------------------+------------------------------------------------------------------------------------------+
   | ``OBS_PRIOR_MEAN_START`` | Global index in obs_ens_handle that holds the first group's prior mean.                  |
   +--------------------------+------------------------------------------------------------------------------------------+
   | ``OBS_PRIOR_MEAN_END``   | Global index in obs_ens_handle that holds the last group's prior mean.                   |
   +--------------------------+------------------------------------------------------------------------------------------+
   | ``OBS_PRIOR_VAR_START``  | Global index in obs_ens_handle that holds the first group's prior variance.              |
   +--------------------------+------------------------------------------------------------------------------------------+
   | ``OBS_PRIOR_VAR_END``    | Global index in obs_ens_handle that holds the last group's prior variance.               |
   +--------------------------+------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call filter_state_space_diagnostics(out_unit, ens_handle, model_size, num_output_state_members,
   output_state_mean_index, output_state_spread_index, output_inflation, temp_ens, ENS_MEAN_COPY, ENS_SD_COPY, inflate,
   INF_COPY, INF_SD_COPY)*
   ::

      type(netcdf_file_type), intent(inout)   :: out_unit
      type(ensemble_type), intent(inout)      :: ens_handle
      integer, intent(in)                     :: model_size
      integer, intent(in)                     :: num_output_state_members
      integer, intent(in)                     :: output_state_mean_index
      integer, intent(in)                     :: output_state_spread_index
      logical, intent(in)                     :: output_inflation
      real(r8), intent(out)                   :: temp_ens(model_size)
      integer, intent(in)                     :: ENS_MEAN_COPY
      integer, intent(in)                     :: ENS_SD_COPY
      type(adaptive_inflate_type), intent(in) :: inflate
      integer, intent(in)                     :: INF_COPY
      integer, intent(in)                     :: INF_SD_COPY

.. container:: indent1

   Writes state space diagnostic values including ensemble members, mean and spread, and inflation mean and spread to a
   netcdf file.

   ============================= ==================================================================
   ``out_unit``                  Descriptor for the netcdf file being written.
   ``ens_handle``                Ensemble handle whose state space values are to be written.
   ``model_size``                Size of the model state vector.
   ``num_output_state_members``  Number of individual state members to be output.
   ``output_state_mean_index``   Index in netcdf file for ensemble mean.
   ``output_state_spread_index`` Index in netcdf file for ensemble spread.
   ``output_inflation``          True if the inflation values are to be output. Default is .TRUE.
   ``temp_ens``                  Storage passed in to avoid having to allocate extra space.
   ``ENS_MEAN_COPY``             Global index in ens_handle for ensemble mean.
   ``ENS_SD_COPY``               Global index in ens_handle for ensemble spread.
   ``inflate``                   Contains description and values of state space inflation.
   ``INF_COPY``                  Global index in ens_handle of inflation values.
   ``INF_SD_COPY``               Global index in ens_handle of inflation standard deviation values.
   ============================= ==================================================================

| 

.. container:: routine

   *call smoother_ss_diagnostics(model_size, num_output_state_members, output_inflation, temp_ens, ENS_MEAN_COPY,
   ENS_SD_COPY, POST_INF_COPY, POST_INF_SD_COPY)*
   ::

      integer, intent(in)   :: model_size
      integer, intent(in)   :: num_output_state_members
      logical, intent(in)   :: output_inflation
      real(r8), intent(out) :: temp_ens(model_size)
      integer, intent(in)   :: ENS_MEAN_COPY
      integer, intent(in)   :: ENS_SD_COPY
      integer, intent(in)   :: POST_INF_COPY
      integer, intent(in)   :: POST_INF_SD_COPY

.. container:: indent1

   Outputs state space diagnostics files for all smoother lags.

   +------------------------------+--------------------------------------------------------------------------------------+
   | ``model_size``               | Size of the model state vector.                                                      |
   +------------------------------+--------------------------------------------------------------------------------------+
   | ``num_output_state_members`` | Number of state copies to be output in the state space diagnostics file.             |
   +------------------------------+--------------------------------------------------------------------------------------+
   | ``output_inflation``         | True if the inflation values are to be output. Default is .TRUE.                     |
   +------------------------------+--------------------------------------------------------------------------------------+
   | ``temp_ens``                 | Storage passed in to avoid having to allocate extra space.                           |
   +------------------------------+--------------------------------------------------------------------------------------+
   | ``ENS_MEAN_COPY``            | Global index of the ensemble mean in the lag smoother ensemble handles.              |
   +------------------------------+--------------------------------------------------------------------------------------+
   | ``ENS_SD_COPY``              | Global index of the ensemble spread in the lag smoother ensemble handles.            |
   +------------------------------+--------------------------------------------------------------------------------------+
   | ``POST_INF_COPY``            | Global index of the inflation value in the lag smoother ensemble handles (not        |
   |                              | currently used).                                                                     |
   +------------------------------+--------------------------------------------------------------------------------------+
   | ``POST_INF_SD_COPY``         | Global index of the inflation spread in the lag smoother ensemble handles (not       |
   |                              | currently used).                                                                     |
   +------------------------------+--------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call smoother_end()*

.. container:: indent1

   Releases storage allocated for smoother.

| 

.. container:: routine

   *call smoother_inc_lags()*

.. container:: indent1

   Increments the number of lags that are in use for smoother. Used when a smoother is being started up and there have
   not been enough times to propagate the state to all requested lags.

| 

Files
-----

-  input.nml
-  smoother initial condition files
-  smoother restart files

References
----------

#. none

Private components
------------------

N/A
