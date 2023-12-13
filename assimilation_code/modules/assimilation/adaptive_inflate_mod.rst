MODULE adaptive_inflate_mod
===========================

Overview
--------

This module implements a variety of hierarchical Bayesian adaptive inflation algorithms for use with ensemble filters.
It can provide constant valued inflation in state or observation space, consistent with previous DART releases. It can
provide spatially-constant, time-varying adaptive inflation. It can provide spatially-varying, time-varying adaptive
inflation and it can provide temporally-varying observation space inflation. And finally, it can provide adaptive damped
inflation, which decreases inflation through time when observation density varies. Diagnostic output and restart files
are available. Several papers on the NSF NCAR `DART <https://dart.ucar.edu/publications/>`__ website document the
algorithms in detail. The ``DART/tutorial/section12`` chapter has more information.

Details on controlling the inflation options are contained in the documentation for the filter. The filter_nml controls
what inflation options are used.

Inflation flavor 3 (spatially-constant state space) reads and writes a restart file that is the full size of the state
vector, however it takes the first value in the array and replicates that throughout the array. This allows one to
switch between flavors 2 and 3. Going from inflation flavor 3 to 2 the initial value for all items in the state vector
will be a constant value and will then start to adapt. Going from inflation flavor 2 to 3 whatever value is in the array
at index 1 will be replicated and used for the entire rest of the state vector items.

Other modules used
------------------

::

   types_mod
   utilities_mod
   random_seq_mod
   time_manager_mod
   ensemble_manager_mod

Public interfaces
-----------------

================================== ==========================
*use adaptive_inflate_mod, only :* update_inflation
\                                  adaptive_inflate_end
\                                  inflate_ens
\                                  output_inflate_diagnostics
\                                  do_obs_inflate
\                                  do_single_ss_inflate
\                                  do_varying_ss_inflate
\                                  adaptive_inflate_init
\                                  adaptive_inflate_type
\                                  get_inflate
\                                  set_inflate
\                                  set_sd
\                                  set_sd
\                                  deterministic_inflate
================================== ==========================

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *call update_inflation(inflate_handle, inflate, inflate_sd, prior_mean, prior_var, obs, obs_var, gamma)*
   ::

      type(adaptive_inflate_type), intent(in)    :: inflate_handle
      real(r8),                    intent(inout) :: inflate
      real(r8),                    intent(inout) :: inflate_sd
      real(r8),                    intent(in)    :: prior_mean
      real(r8),                    intent(in)    :: prior_var
      real(r8),                    intent(in)    :: obs
      real(r8),                    intent(in)    :: obs_var
      real(r8),                    intent(in)    :: gamma

.. container:: indent1

   Updates the mean and standard deviation of an inflation distribution given the prior values, the prior observation
   ensemble mean and variance, and the observation and its error variance. The factor gamma is the expected impact (0 to
   1) of the state variable corresponding to the inflation on the observation and is the product of the ensemble
   correlation plus an additional localization factor or group regression factors.

   ================== ================================================================================
   ``inflate_handle`` Handle to object that describes the inflation type and values.
   ``inflate``        Prior mean value of the inflation distribution.
   ``inflate_sd``     Prior standard deviation of the inflation distribution.
   ``prior_mean``     Mean of the prior observation ensemble.
   ``prior_var``      Variance of the prior observation ensemble.
   ``obs``            The observed value.
   ``obs_var``        Observational error variance.
   ``gamma``          Expected impact factor, product of correlation, localization, regression factor.
   ================== ================================================================================

| 

.. container:: routine

   *call adaptive_inflate_end(inflate_handle, ens_handle, ss_inflate_index, ss_inflate_sd_index)*
   ::

      type(adaptive_inflate_type), intent(in)    :: inflate_handle
      type(ensemble_type),         intent(inout) :: ens_handle
      integer,                     intent(in)    :: ss_inflate_index
      integer,                     intent(in)    :: ss_inflate_sd_index

.. container:: indent1

   Outputs the values of inflation to restart files using the ensemble_manager for state space inflation and file output
   for observation space inflation. Releases allocated storage in inflate_handle.

   ======================= ==============================================================================
   ``inflate_handle``      Handle for the details of the inflation being performed.
   ``ens_handle``          Handle for ensemble storage that holds values of state space inflation.
   ``ss_inflate_index``    Index in ensemble storage copies for state space inflation.
   ``ss_inflate_sd_index`` Index in ensemble storage copies for state space inflation standard deviation.
   ======================= ==============================================================================

| 

.. container:: routine

   *call inflate_ens(inflate_handle, ens,mean, inflate [,var_in])*
   ::

      type(adaptive_inflate_type),               intent(in)  :: inflate_handle
      real(r8),                    dimension(:), intent(out) :: ens
      real(r8),                                  intent(in)  :: mean
      real(r8),                                  intent(in)  :: inflate
      real(r8),                    optional,     intent(in)  :: var_in

.. container:: indent1

   Given an ensemble, its mean and the covarance inflation factor, inflates the ensemble.

   ================== ========================================================
   ``inflate_handle`` Handle for the details of the inflation being performed.
   ``ens``            Values for the ensemble to be inflated
   ``mean``           The mean of the ensemble.
   ``inflate``        The covariance inflation factor.
   ``var_in``         The variance of the ensemble.
   ================== ========================================================

| 

.. container:: routine

   *call output_inflate_diagnostics(inflate_handle, time)*
   ::

      type(adaptive_inflate_type), intent(in) :: inflate_handle
      type(time_type),             intent(in) :: time

.. container:: indent1

   Outputs diagnostic record of inflation for the observation space of spatially constant state space inflation.
   Spatially varying state space diagnostics are in the Posterior and Prior Diagnostic netcdf files and are written with
   calls from filter.f90.

   ================== ========================================================
   ``inflate_handle`` Handle for the details of the inflation being performed.
   ``time``           Time of this diagnostic info.
   ================== ========================================================

| 

.. container:: routine

   *var = do_obs_inflate(inflate_handle)*
   ::

      logical,               intent(out) :: do_obs_inflate
      adaptive_inflate_type, intent(in)  :: inflate_handle

.. container:: indent1

   Returns true if observation space inflation is being done by this handle.

   ================== =========================================================
   ``do_obs_inflate`` True if obs space inflation is being done by this handle.
   ``inflate_handle`` Handle to inflation details.
   ================== =========================================================

| 

.. container:: routine

   *var = do_varying_ss_inflate(inflate_handle)*
   ::

      logical,               intent(out) :: do_varying_ss_inflate
      adaptive_inflate_type, intent(in)  :: inflate_handle

.. container:: indent1

   Returns true if spatially varying state space inflation is being done by this handle.

   ========================= =============================================================================
   ``do_varying_ss_inflate`` True if spatially varying state space inflation is being done by this handle.
   ``inflate_handle``        Handle to inflation details.
   ========================= =============================================================================

| 

.. container:: routine

   *var = do_single_ss_inflate(inflate_handle)*
   ::

      logical,               intent(out) :: do_single_ss_inflate
      adaptive_inflate_type, intent(in)  :: inflate_handle

.. container:: indent1

   Returns true if spatially fixed state space inflation is being done by this handle.

   ======================== ===========================================================================
   ``do_single_ss_inflate`` True if spatially fixed state space inflation is being done by this handle.
   ``inflate_handle``       Handle to inflation details.
   ======================== ===========================================================================

| 

.. container:: routine

   *call adaptive_inflate_init(inflate_handle, inf_flavor, mean_from_restart, sd_from_restart, output_restart,
   deterministic, in_file_name, out_file_name, diag_file_name, inf_initial, sd_initial, inf_lower_bound,
   inf_upper_bound, sd_lower_bound, ens_handle, ss_inflate_index, ss_inflate_sd_index, label)*
   ::

      type(adaptive_inflate_type), intent(inout) :: inflate_handle
      integer, intent(in)                        :: inf_flavor
      logical, intent(in)                        :: mean_from_restart
      logical, intent(in)                        :: sd_from_restart
      logical, intent(in)                        :: output_restart
      logical, intent(in)                        :: deterministic
      character(len=*), intent(in)               :: in_file_name
      character(len=*), intent(in)               :: out_file_name
      character(len=*), intent(in)               :: diag_file_name
      real(r8), intent(in)                       :: inf_initial
      real(r8), intent(in)                       :: sd_initial
      real(r8), intent(in)                       :: inf_lower_bound
      real(r8), intent(in)                       :: inf_upper_bound
      real(r8), intent(in)                       :: sd_lower_bound
      type(ensemble_type), intent(inout)         :: ens_handle
      integer, intent(in)                        :: ss_inflate_index
      integer, intent(in)                        :: ss_inflate_sd_index
      character(len=*), intent(in)               :: label

.. container:: indent1

   Initializes a descriptor of an inflation object.

   ======================= ============================================================================
   ``inflate_handle``      Handle for the inflation descriptor being initialized.
   ``inf_flavor``          Type of inflation, 1=obs_inflate, 2=varying_ss_inflate, 3=single_ss_inflate.
   ``mean_from_restart``   True if inflation mean values to be read from restart file.
   ``sd_from_restart``     True if inflation standard deviation values to be read from restart file.
   ``output_restart``      True if an inflation restart file is to be output.
   ``deterministic``       True if deterministic inflation is to be done.
   ``in_file_name``        File name from which to read restart.
   ``out_file_name``       File name to which to write restart.
   ``diag_file_name``      File name to which to write diagnostic output; obs space inflation only .
   ``inf_initial``         Initial value of inflation for start_from_restart=.false.
   ``sd_initial``          Initial value of inflation standard deviation for start_from_restart=.false.
   ``inf_lower_bound``     Lower bound on inflation value.
   ``inf_upper_bound``     Upper bound on inflation value.
   ``sd_lower_bound``      Lower bound on inflation standard deviation.
   ``ens_handle``          Ensemble handle with storage for state space inflation.
   ``ss_inflate_index``    Index op copy in ensemble storage for inflation value.
   ``ss_inflate_sd_index`` Index of copy in ensemble storage for inflation standard deviation.
   ``label``               Character label to be used in diagnostic output (e.g. 'Prior', 'Posterior').
   ======================= ============================================================================

| 

.. container:: routine

   *var = get_sd(inflate_handle)*
   ::

      real(r8), intent(out)                   :: get_sd
      type(adaptive_inflate_type), intent(in) :: inflate_handle

.. container:: indent1

   Returns value of observation space inflation standard deviation.

   ================== =================================================
   ``get_sd``         Returns the value of observation space inflation.
   ``inflate_handle`` Handle for inflation descriptor.
   ================== =================================================

| 

.. container:: routine

   *var = get_inflate(inflate_handle)*
   ::

      real(r8), intent(out)                   :: get_inflate
      type(adaptive_inflate_type), intent(in) :: inflate_handle

.. container:: indent1

   Returns value of observation space inflation.

   ================== =================================================
   ``get_inflate``    Returns the value of observation space inflation.
   ``inflate_handle`` Handle for inflation descriptor.
   ================== =================================================

| 

.. container:: routine

   *call set_inflate(inflate_handle, inflate)*
   ::

      type(adaptive_inflate_type), intent(inout) :: inflate_handle
      real(r8), intent(in)                       :: inflate

.. container:: indent1

   Set the value of observation space inflation.

   ================== ==============================================
   ``inflate_handle`` Handle for inflation descriptor.
   ``inflate``        Set observation space inflation to this value.
   ================== ==============================================

| 

.. container:: routine

   *call set_sd(inflate_handle, sd)*
   ::

      type(adaptive_inflate_type), intent(inout) :: inflate_handle
      real(r8), intent(in)                       :: sd

.. container:: indent1

   Set the value of observation space inflation standard deviation.

   ================== =================================================================
   ``inflate_handle`` Handle for inflation descriptor.
   ``sd``             Set observation space inflation standard deviation to this value.
   ================== =================================================================

| 

.. container:: routine

   *var = deterministic_inflate(inflate_handle)*
   ::

      logical, intent(out)                    :: deterministic_inflate
      type(adaptive_inflate_type), intent(in) :: inflate_handle

.. container:: indent1

   Returns true if deterministic inflation is being done.

   ========================= ======================================================
   ``deterministic_inflate`` Returns true if deterministic inflation is being done.
   ``inflate_handle``        Handle for inflation descriptor.
   ========================= ======================================================

| 

.. container:: type

   ::

      type adaptive_inflate_type
         private
         integer :: inflation_flavor
         integer :: obs_diag_unit
         logical :: start_from_restart
         logical :: output_restart
         logical :: deterministic
         character(len = 129) :: in_file_name
         character(len = 129) :: out_file_name
         character(len = 129) :: diag_file_name
         real(r8) :: inflate
         real(r8) :: sd
         real(r8) :: sd_lower_bound
         real(r8) :: inf_lower_bound
         real(r8) :: inf_upper_bound
         type(random_seq_type) :: ran_seq
      end type adaptive_inflate_type

.. container:: indent1

   Provides a handle for a descriptor of inflation. Includes type of inflation, values controlling it, input and output
   file names, an output file descriptor for observation space inflation diagnotics, and a random sequence for doing
   reproducible non-determinstic inflation. There are 2 instances of this type, one for Prior and one for Posterior
   inflation.

   ================== ================================================================================
   Component          Description
   ================== ================================================================================
   inflation_flavor   Type of inflation; 0=none, 1=obs. space, 2=spatially varying, 3=spatially-fixed.
   obs_diag_unit      Unit descriptor for output diagnostic file.
   start_from_restart True if initial inflation to be read from file.
   output_restart     True if final inflation values to be written to file.
   deterministic      True if inflation is to be done be deterministic algorithm.
   in_file_name       File name containing restart.
   out_file_name      File to contain output restart.
   diag_file_name     File to hold observation space diagnostics.
   inflate            Initial value of inflation for all types; current value for obs. space.
   sd                 Initial value of sd for all types; current value for obs. space.
   sd_lower_bound     Don't allow standard deviation to get smaller than this.
   inf_lower_bound    Don't let inflation get smaller than this.
   inf_upper_bound    Don't let inflation get larger than this.
   ran_seq            Handle to random number sequence to allow reproducing non-deterministic inflate.
   ================== ================================================================================

| 

Namelist
--------

The adaptive_inflate module no longer has a namelist. Control has been moved to
`&filter_nml <filter_mod.html#Namelist>`__ in filter.

Files
-----

Three files are opened from this module, but all names are passed in from the filter_nml now, and there are 2 values for
each name: one for the prior and one for the posterior inflation.

-  inf_in_file_name
   Mean and standard deviation values read in restart file format.
-  inf_out_file_name
   Mean and standard deviation values written in restart file format.
-  inf_diag_file_name
   Contains diagnostic history of inflation values for obs space and spatially-fixed state space inflation. Diagnostics
   for spatially-varying state space inflation are extra fields on the Posterior and Prior diagnostic netcdf files
   created in filter.f90.

References
----------

-  Anderson, J. L., 2007: An adaptive covariance inflation error correction algorithm for ensemble filters. Tellus A,
   59, 210-224.
   `doi: 10.1111/j.1600-0870.2006.00216.x <http://dx.doi.org/10.1111/j.1600-0870.2006.00216.x>`__
-  Anderson, J. L., 2009: Spatially and temporally varying adaptive covariance inflation for ensemble filters. Tellus A,
   61, 72-83.
   `doi: 10.1111/j.1600-0870.2008.00361.x <http://dx.doi.org/10.1111/j.1600-0870.2008.00361.x>`__

Private components
------------------

no discussion
