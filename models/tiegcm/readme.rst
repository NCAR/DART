TIEGCM
======

Overview
--------

| This is the DART interface to the Thermosphere Ionosphere Electrodynamic General Circulation Model
  (`TIEGCM <http://www.hao.ucar.edu/modeling/tgcm/tie.php>`__), which is a community model developed at the NCAR High
  Altitude Observatory. TIEGCM is widely used by the space physics and aeronomy community and is one of the most
  well-validated models of the Earth's upper atmosphere. DART/TIEGCM has been used to assimilate neutral mass density
  retrieved from satellite-borne accelerometers and electon density obtained from ground-based and space-based GNSS
  signals. Unlike other ionospheric data assimilation applications, this approach allows simultaneous assimilation of
  thermospheric and ionospheric parameters by taking advantage of the coupling of plasma and neutral constituents
  described in TIEGCM. DART/TIEGCM's demonstrated capability to infer under-observed thermospheric parameters from
  abundant electron density observations has important implications for the future of upper atmosphere research.
| DART is designed so that the TIEGCM source code can be used with no modifications, as DART runs TIEGCM as a completely
  separate executable. The TIEGCM source code and restart files are **not** included in DART, so you must obtain them
  from the NCAR High Altitude Observatory (`download website <http://www.hao.ucar.edu/modeling/tgcm/download.php>`__).
  It is **strongly** recommended that you become familiar with running TIEGCM **before** you try to run DART/TIEGCM (See
  the `TIEGCM User's Guide <http://www.hao.ucar.edu/modeling/tgcm/doc/userguide/html>`__). Some assumptions are made
  about the mannner in which TIEGCM is run: (1) There can only be 1 each of the TIEGCM primary (restart) and secondary
  NetCDF history files. The TIEGCM primary history files contain the prognostic variables necessary to restart the
  model, while the secondary history files contain diagnostic variables; (2) The last timestep in the restart file is
  the only timestep which is converted to a DART state vector, and only the last timestep in the TIEGCM primary file is
  ever modified by DART. The TIEGCM variables to be included in a DART state vector, and possibly updated by the
  assimilation, are specified in the DART namelist. (Some of the TIEGCM variables used to compute observation priors
  need not to be updated.) It is required to associate the TIEGCM variable name with a 'generic' DART counterpart (e.g.,
  ``NE`` is ``QTY_ELECTRON_DENSITY``). The composition of the DART state vector and which variables get updated in the
  TIEGCM primary file are under complete user control.
| In the course of a filtering experiment, it is necessary to make a short forecast with TIEGCM. DART writes out an
  ancillary file with the information necessary to advance TIEGCM to the required time. The DART script
  ``advance_model.csh`` reads this information and modifies the TIEGCM namelist ``tiegcm.nml`` such that TIEGCM runs
  upto the requested time when DART assimilates the next set of observations. The run scripts ``run_filter.csh`` and
  ``run_perfect_model_obs.csh`` are configured to run under the LSF queueing system. The scripting examples exploit an
  'embarassingly-simple' parallel paradigm in that each TIEGCM instance is a single-threaded executable and all ensemble
  members may run simultaneously. To use these run scripts, the TIECGM executable needs to be compiled with no MPI
  option. As such, there is an advantage to matching the ensemble size to the number of tasks. Requesting more tasks
  than the number of ensemble members may speed up the DART portion of an assimilation (i.e., ``filter``) but will not
  make the model advance faster. The ``filter`` may be compiled with MPI and can exploit all available tasks.

Quickstart guide to running
---------------------------

| It is important to understand basic DART nomenclature and mechanisms. Please take the time to read and run the `DART
  tutorial <../../tutorial/index.pdf>`__.
| Both ``run_filter.csh`` and ``run_perfect_model_obs.csh`` are heavily internally commented. Please read and understand
  the scripts. The overall process is to

#. Specify resources (wall-clock time, number of nodes, tasks that sort of thing).
#. Set shell variables to identify the location of the DART exectuables, the TIEGCM executables, initial ensemble, etc.
#. Establish a temporary working directory for the experiment.
#. Populate that directory with the initial ensemble and required namelists.
#. Convert each TIEGCM ensemble member to a DART initial conditions file.
#. Run either ``filter`` or ``run_perfect_model_obs.csh``.

#. ``perfect_model_obs`` will
#. Check for any desired observations at the current time of the model state and create the synthetic observations for
   all observation times in the specified assimilation window. If the model needs to be advanced, it then
#. creates a unique run-time directory for the model advance,
#. copies the required information into that directory,
#. conveys the desired forecast stopping time to TIEGCM via the ``tiegcm.nml`` and
#. runs a single executable of TIEGCM.
#. Steps 1-5 are repeated until the input DART observation sequence file has been exhausted.

#. ``filter`` will
#. Check for any desired observations at the current time of the model state and assimilates all the observations in the
   specified assimilation window. If the model needs to be advanced, it then
#. creates a set of run-time directories, one for each task. A single task may be responsible for advancing more than
   one TIEGCM instance. If so, each instance is done serially, one after another. See the documentation for
   :doc:`../../docs/html/filter_async_modes`.
#. Copy the required information into that directory.
#. Update the TIEGCM restart file with the most current DART-modified state and convey the desired forecast stopping
   time to TIEGCM via the unique ``tiegcm.nml`` for this ensemble member.
#. Runs a single executable of TIEGCM.
#. Steps 1-5 are repeated until the input DART observation sequence file

What to check when things go wrong
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The scripts are designed to send email to the user that contains the run-time output from the script. Check that first.
If that does not provide the information needed, go to the run directory (i.e. CENTRALDIR) and check the
``dart_log.out``. It usually provides the same information as the email, but sometimes it can help. If that does not
help, go to any of the CENTRALDIR/*advance_temp\ nnnn* directories and read the *log_advance.\ nnnn.txt* file.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &model_nml 
      output_state_vector         = .false.
      tiegcm_restart_file_name    = 'tiegcm_restart_p.nc'
      tiegcm_secondary_file_name  = 'tiegcm_s.nc'
      tiegcm_namelist_file_name   = 'tiegcm.nml'
      assimilation_period_seconds = 3600
      estimate_f10_7              = .false.
      debug                       = 1
      variables = 'NE',    'QTY_ELECTRON_DENSITY',          '1000.0',  'NA',      'restart',    'UPDATE'
                  'OP',    'QTY_DENSITY_ION_OP',            'NA',      'NA',      'restart',    'UPDATE',
                  'TI',    'QTY_TEMPERATURE_ION',           'NA',      'NA',      'restart',    'UPDATE',
                  'TE',    'QTY_TEMPERATURE_ELECTRON',      'NA',      'NA',      'restart',    'UPDATE',
                  'OP_NM', 'QTY_DENSITY_ION_OP',            'NA',      'NA',      'restart',    'UPDATE',
                  'O1',    'QTY_ATOMIC_OXYGEN_MIXING_RATIO','0.00001', '0.99999', 'secondary',  'NO_COPY_BACK',
                  'O2',    'QTY_MOLEC_OXYGEN_MIXING_RATIO', '0.00001', '0.99999', 'secondary',  'NO_COPY_BACK',
                  'TN',    'QTY_TEMPERATURE',               '0.0',     '6000.0',  'secondary',  'NO_COPY_BACK',
                  'ZG',    'QTY_GEOMETRIC_HEIGHT',          'NA',      'NA',      'secondary',  'NO_COPY_BACK',
                  'VTEC',  'QTY_VERTICAL_TEC',              'NA',      'NA',      'calculate',  'NO_COPY_BACK'
      /

.. container::

   +---------------------------------------+---------------------------------------+---------------------------------------+
   | Item                                  | Type                                  | Description                           |
   +=======================================+=======================================+=======================================+
   | output_state_vector                   | logical                               | If .true. write state vector as a 1D  |
   |                                       |                                       | array to the DART diagnostic output   |
   |                                       |                                       | files. If .false. break state vector  |
   |                                       |                                       | up into variables before writing to   |
   |                                       |                                       | the output files.                     |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | tiegcm_restart_file_name              | character(len=256)                    | The TIEGCM restart file name.         |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | tiegcm_secondary_file_name            | character(len=256)                    | The TIEGCM secondary file name.       |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | tiegcm_namelist_file_name             | character(len=256)                    | The TIEGCM namelist file name.        |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | assimilation_period_seconds           | integer                               | This specifies the width of the       |
   |                                       |                                       | assimilation window. The current      |
   |                                       |                                       | model time is used as the center time |
   |                                       |                                       | of the assimilation window. All       |
   |                                       |                                       | observations in the assimilation      |
   |                                       |                                       | window are assimilated. BEWARE: if    |
   |                                       |                                       | you put observations that occur       |
   |                                       |                                       | before the beginning of the           |
   |                                       |                                       | assimilation_period, DART will error  |
   |                                       |                                       | out because it cannot move the model  |
   |                                       |                                       | 'back in time' to process these       |
   |                                       |                                       | observations.                         |
   |                                       |                                       | ``assimilation_period_seconds`` must  |
   |                                       |                                       | be an integer number of TIEGCM        |
   |                                       |                                       | dynamical timesteps (as specified by  |
   |                                       |                                       | tiegcm.nml:STEP) AND be able to be    |
   |                                       |                                       | expressed by tiegcm.nml:STOP. Since   |
   |                                       |                                       | STOP has three components:            |
   |                                       |                                       | day-of-year, hour, and minute, the    |
   |                                       |                                       | ``assimilation_period_seconds`` must  |
   |                                       |                                       | be an integer number of minutes.      |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | estimate_f10_7                        | logical                               | Switch to specify that the f10.7      |
   |                                       |                                       | index should be estimated by          |
   |                                       |                                       | augmenting the DART state vector with |
   |                                       |                                       | a scalar. The location of the f10.7   |
   |                                       |                                       | index is taken to be longitude of     |
   |                                       |                                       | local noon and latitude zero.         |
   |                                       |                                       | WARNING: this is provided with no     |
   |                                       |                                       | guarantees. Please read the comments  |
   |                                       |                                       | in ``model_mod.f90`` and act          |
   |                                       |                                       | accordingly.                          |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | debug                                 | integer                               | Set to 0 (zero) for minimal output.   |
   |                                       |                                       | Successively larger values generate   |
   |                                       |                                       | successively more output.             |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | variables                             | character(:,6)                        | Strings that identify the TIEGCM      |
   |                                       |                                       | variables, their DART kind, the min & |
   |                                       |                                       | max values, what file to read from,   |
   |                                       |                                       | and whether or not the file should be |
   |                                       |                                       | updated after the assimilation. The   |
   |                                       |                                       | DART kind must be one found in the    |
   |                                       |                                       | ``DART/assimilation_code/mo           |
   |                                       |                                       | dules/observations/obs_kind_mod.f90`` |
   |                                       |                                       | AFTER it gets built by                |
   |                                       |                                       | ``preprocess``. Most of the upper     |
   |                                       |                                       | atmosphere observation kinds are      |
   |                                       |                                       | specified by                          |
   |                                       |                                       | ``DART/observations/forward_o         |
   |                                       |                                       | perators/obs_def_upper_atm_mod.f90``, |
   |                                       |                                       | so it should be specified in the      |
   |                                       |                                       | ``preprocess_nml``:``input_files``    |
   |                                       |                                       | variable. Since TIEGCM has an entire  |
   |                                       |                                       | class of variables (all the variables |
   |                                       |                                       | that end in ``_NM``) that are simply  |
   |                                       |                                       | 1 dynamical timestep behind the       |
   |                                       |                                       | variables at the output time, it is   |
   |                                       |                                       | **imperative** that these variables   |
   |                                       |                                       | be specified to occur AFTER their     |
   |                                       |                                       | counterparts in the DART namelist.    |
   |                                       |                                       | This will ensure that the most        |
   |                                       |                                       | current variables are used in the     |
   |                                       |                                       | calculation of the forward            |
   |                                       |                                       | observation operators.                |
   |                                       |                                       |                                       |
   |                                       |                                       | +----------------+----------------+   |
   |                                       |                                       | | ``va           | Specifies the  |   |
   |                                       |                                       | | riables(:,1)`` | TIEGCM         |   |
   |                                       |                                       | |                | variable name  |   |
   |                                       |                                       | |                | in the netCDF  |   |
   |                                       |                                       | |                | file.          |   |
   |                                       |                                       | +----------------+----------------+   |
   |                                       |                                       | | ``va           | Specifies the  |   |
   |                                       |                                       | | riables(:,2)`` | DART kind for  |   |
   |                                       |                                       | |                | that variable. |   |
   |                                       |                                       | +----------------+----------------+   |
   |                                       |                                       | | ``va           | Specifies a    |   |
   |                                       |                                       | | riables(:,3)`` | minimum bound  |   |
   |                                       |                                       | |                | (if any) for   |   |
   |                                       |                                       | |                | that variable. |   |
   |                                       |                                       | +----------------+----------------+   |
   |                                       |                                       | | ``va           | Specifies a    |   |
   |                                       |                                       | | riables(:,4)`` | maximum bound  |   |
   |                                       |                                       | |                | (if any) for   |   |
   |                                       |                                       | |                | that variable. |   |
   |                                       |                                       | +----------------+----------------+   |
   |                                       |                                       | | ``va           | Specifies what |   |
   |                                       |                                       | | riables(:,5)`` | file the       |   |
   |                                       |                                       | |                | variable       |   |
   |                                       |                                       | |                | should come    |   |
   |                                       |                                       | |                | from. The only |   |
   |                                       |                                       | |                | valid          |   |
   |                                       |                                       | |                | possibilies    |   |
   |                                       |                                       | |                | are "restart", |   |
   |                                       |                                       | |                | "secondary",   |   |
   |                                       |                                       | |                | or             |   |
   |                                       |                                       | |                | "calculate".   |   |
   |                                       |                                       | |                | "restart" will |   |
   |                                       |                                       | |                | read from      |   |
   |                                       |                                       | |                | whatever file  |   |
   |                                       |                                       | |                | is specified   |   |
   |                                       |                                       | |                | by             |   |
   |                                       |                                       | |                | `              |   |
   |                                       |                                       | |                | `tiegcm_restar |   |
   |                                       |                                       | |                | t_file_name``. |   |
   |                                       |                                       | |                | "secondary"    |   |
   |                                       |                                       | |                | will read from |   |
   |                                       |                                       | |                | whatever file  |   |
   |                                       |                                       | |                | is specified   |   |
   |                                       |                                       | |                | by             |   |
   |                                       |                                       | |                | ``t            |   |
   |                                       |                                       | |                | iegcm_secondar |   |
   |                                       |                                       | |                | y_file_name``. |   |
   |                                       |                                       | |                | "calculate"    |   |
   |                                       |                                       | |                | will call a    |   |
   |                                       |                                       | |                | vari           |   |
   |                                       |                                       | |                | able-dependent |   |
   |                                       |                                       | |                | function --    |   |
   |                                       |                                       | |                | see            |   |
   |                                       |                                       | |                | ``m            |   |
   |                                       |                                       | |                | odel_mod.f90`` |   |
   |                                       |                                       | |                | :``tiegcm_to_d |   |
   |                                       |                                       | |                | art_vector()`` |   |
   |                                       |                                       | |                | for the        |   |
   |                                       |                                       | |                | ``c            |   |
   |                                       |                                       | |                | reate_vtec()`` |   |
   |                                       |                                       | |                | example.       |   |
   |                                       |                                       | +----------------+----------------+   |
   |                                       |                                       | | ``va           | Specifies if   |   |
   |                                       |                                       | | riables(:,6)`` | the variable   |   |
   |                                       |                                       | |                | should be      |   |
   |                                       |                                       | |                | updated in the |   |
   |                                       |                                       | |                | TIEGCM restart |   |
   |                                       |                                       | |                | file. The      |   |
   |                                       |                                       | |                | value may be   |   |
   |                                       |                                       | |                | "UPDATE" or    |   |
   |                                       |                                       | |                | anything else. |   |
   |                                       |                                       | |                | If **and only  |   |
   |                                       |                                       | |                | if** the       |   |
   |                                       |                                       | |                | variable comes |   |
   |                                       |                                       | |                | from the       |   |
   |                                       |                                       | |                | restart file   |   |
   |                                       |                                       | |                | **and**        |   |
   |                                       |                                       | |                | ``va           |   |
   |                                       |                                       | |                | riables(:,6)`` |   |
   |                                       |                                       | |                | == "UPDATE"    |   |
   |                                       |                                       | |                | will the       |   |
   |                                       |                                       | |                | variable be    |   |
   |                                       |                                       | |                | modified in    |   |
   |                                       |                                       | |                | the TIEGCM     |   |
   |                                       |                                       | |                | restart file.  |   |
   |                                       |                                       | |                | No variables   |   |
   |                                       |                                       | |                | in the         |   |
   |                                       |                                       | |                | secondary file |   |
   |                                       |                                       | |                | are EVER       |   |
   |                                       |                                       | |                | modified.      |   |
   |                                       |                                       | +----------------+----------------+   |
   +---------------------------------------+---------------------------------------+---------------------------------------+

Other modules used
------------------

::

   adaptive_inflate_mod.f90
   assim_model_mod.f90
   assim_tools_mod.f90
   types_mod.f90
   cov_cutoff_mod.f90
   ensemble_manager_mod.f90
   filter.f90
   location/threed_sphere/location_mod.f90
   [null_,]mpi_utilities_mod.f90
   obs_def_mod.f90
   obs_kind_mod.f90
   obs_model_mod.f90
   obs_sequence_mod.f90
   random_seq_mod.f90
   reg_factor_mod.f90
   smoother_mod.f90
   sort_mod.f90
   time_manager_mod.f90
   utilities_mod.f90

Public interfaces - required
----------------------------

======================= ======================
*use model_mod, only :* get_model_size
\                       adv_1step
\                       get_state_meta_data
\                       model_interpolate
\                       get_model_time_step
\                       static_init_model
\                       end_model
\                       init_time
\                       init_conditions
\                       nc_write_model_atts
\                       nc_write_model_vars
\                       pert_model_state
\                       get_close_maxdist_init
\                       get_close_obs_init
\                       get_close_obs
\                       ens_mean_for_model
======================= ======================

Public interfaces - optional
----------------------------

======================= =====================
*use model_mod, only :* tiegcm_to_dart_vector
\                       dart_vector_to_tiegcm
\                       get_f107_value
\                       test_interpolate
======================= =====================

A namelist interface ``&model_nml`` is defined by the module, and is read from file ``input.nml``.

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *model_size = get_model_size( )*
   ::

      integer :: get_model_size

.. container:: indent1

   Returns the length of the model state vector. Required.

   ============== =====================================
   ``model_size`` The length of the model state vector.
   ============== =====================================

| 

.. container:: routine

   *call adv_1step(x, time)*
   ::

      real(r8), dimension(:), intent(inout) :: x
      type(time_type),        intent(in)    :: time

.. container:: indent1

   Since TIEGCM is not called as a subroutine, this is a NULL interface. TIEGCM is advanced as a separate executable -
   i.e. ``async == 2``. *adv_1step* only gets called if ``async == 0``. The subroutine must still exist, but contains no
   code and will not be called. An error message is issued if an unsupported value of
   ``filter,perfect_model_obs``:``async`` is used.

| 

.. container:: routine

   *call get_state_meta_data (index_in, location, [, var_kind] )*
   ::

      integer,             intent(in)  :: index_in
      type(location_type), intent(out) :: location
      integer, optional,   intent(out) ::  var_kind 

.. container:: indent1

   Given an integer index into the state vector structure, returns the associated location. A second intent(out)
   optional argument returns the generic kind of this item, e.g. QTY_MOLEC_OXYGEN_MIXING_RATIO, QTY_ELECTRON_DENSITY,
   ... This interface is required to be functional for all applications.

   ============ ===================================================================
   ``index_in`` Index of state vector element about which information is requested.
   ``location`` The location of state variable element.
   *var_kind*   The generic kind of the state variable element.
   ============ ===================================================================

| 

.. container:: routine

   *call model_interpolate(x, location, ikind, obs_val, istatus)*
   ::

      real(r8), dimension(:), intent(in)  :: x
      type(location_type),    intent(in)  :: location
      integer,                intent(in)  :: ikind
      real(r8),               intent(out) :: obs_val
      integer,                intent(out) :: istatus

.. container:: indent1

   Given a state vector, a location, and a model state variable kind interpolates the state variable field to that
   location and returns the value in obs_val. The istatus variable should be returned as 0 unless there is some problem
   in computing the interpolation in which case a positive value should be returned. The ikind variable is one of the
   KIND parameters defined in the :doc:`../../assimilation_code/modules/observations/obs_kind_mod` file and defines
   which generic kind of item is being interpolated.

   ============ ========================================================================================
   ``x``        A model state vector.
   ``location`` Location to which to interpolate.
   ``itype``    Kind of state field to be interpolated.
   ``obs_val``  The interpolated value from the model.
   ``istatus``  Integer value returning 0 for success. Other values can be defined for various failures.
   ============ ========================================================================================

| 

.. container:: routine

   *var = get_model_time_step()*
   ::

      type(time_type) :: get_model_time_step

.. container:: indent1

   Returns the smallest useful forecast length (time step) of the model. This is set by
   ``input.nml``:``assimilation_period_seconds`` and must be an integer number of TIEGCM dynamical timesteps (as
   specified by ``tiegcm.nml``:``STEP``) AND be able to be expressed by ``tiegcm.nml``:``STOP``. Since ``STOP`` has
   three components: day-of-year, hour, and minute, the ``assimilation_period_seconds`` must be an integer number of
   minutes.

   ======= ================================
   ``var`` Smallest forecast step of model.
   ======= ================================

| 

.. container:: routine

   *call static_init_model()*

.. container:: indent1

   Called to do one-time initialization of the model. There are no input arguments. ``static_init_model`` reads the DART
   and TIEGCM namelists and reads the grid geometry and constructs the shape of the DART vector given the TIEGCM
   variables specified in the DART namelist.

| 

.. container:: routine

   *call end_model()*

.. container:: indent1

   Does all required shutdown and clean-up needed.

| 

.. container:: routine

   *call init_time(time)*
   ::

      type(time_type), intent(out) :: time

.. container:: indent1

   This is a NULL INTERFACE for TIEGCM. If ``input.nml``:``start_from_restart == .FALSE.``, this routine is called and
   will generate a fatal error.

| 

.. container:: routine

   *call init_conditions(x)*
   ::

      real(r8), dimension(:), intent(out) :: x

.. container:: indent1

   This is a NULL INTERFACE for TIEGCM. If ``input.nml``:``start_from_restart == .FALSE.``, this routine is called and
   will generate a fatal error.

| 

.. container:: routine

   *ierr = nc_write_model_atts(ncFileID)*
   ::

      integer             :: nc_write_model_atts
      integer, intent(in) :: ncFileID

.. container:: indent1

   This routine writes the model-specific attributes to a netCDF file. This includes the coordinate variables and any
   metadata, but NOT the model state vector. We do have to allocate SPACE for the model state vector, but that variable
   gets filled as the model advances. If ``input.nml``:``model_nml:output_state_vector == .TRUE.``, the DART state
   vector is written as one long vector. If ``input.nml``:``model_nml:output_state_vector == .FALSE.``, the DART state
   vector is reshaped into the original TIEGCM variables and those variables are written.

   ============ =========================================================
   ``ncFileID`` Integer file descriptor to previously-opened netCDF file.
   ``ierr``     Returns a 0 for successful completion.
   ============ =========================================================

| 

.. container:: routine

   *ierr = nc_write_model_vars(ncFileID, statevec, copyindex, timeindex)*
   ::

      integer                            :: nc_write_model_vars
      integer,                intent(in) :: ncFileID
      real(r8), dimension(:), intent(in) :: statevec
      integer,                intent(in) :: copyindex
      integer,                intent(in) :: timeindex

.. container:: indent1

   This routine writes the DART state vector to a netCDF file. If
   ``input.nml``:``model_nml:output_state_vector == .TRUE.``, the DART state vector is written as one long vector. If
   ``input.nml``:``model_nml:output_state_vector == .FALSE.``, the DART state vector is reshaped into the original
   TIEGCM variables and those variables are written.

   ============= =================================================
   ``ncFileID``  file descriptor to previously-opened netCDF file.
   ``statevec``  A model state vector.
   ``copyindex`` Integer index of copy to be written.
   ``timeindex`` The timestep counter for the given state.
   ``ierr``      Returns 0 for normal completion.
   ============= =================================================

| 

.. container:: routine

   *call pert_model_state(state, pert_state, interf_provided)*
   ::

      real(r8), dimension(:), intent(in)  :: state
      real(r8), dimension(:), intent(out) :: pert_state
      logical,                intent(out) :: interf_provided

.. container:: indent1

   | ``pert_model_state`` is intended to take a single model state vector and perturbs it in some way to generate
     initial conditions for spinning up ensembles. TIEGCM does this is a manner that is different than most other
     models. The F10_7 parameter must be included in the DART state vector as a QTY_1D_PARAMETER and gaussian noise is
     added to it. That value must be conveyed to the tiegcm namelist and used to advance the model.
   | Most other models simply add noise with certain characteristics to the model state.

   +---------------------+-----------------------------------------------------------------------------------------------+
   | ``state``           | State vector to be perturbed.                                                                 |
   +---------------------+-----------------------------------------------------------------------------------------------+
   | ``pert_state``      | Perturbed state vector.                                                                       |
   +---------------------+-----------------------------------------------------------------------------------------------+
   | ``interf_provided`` | This is returned as .TRUE. since the routine exists. A value of .FALSE. would indicate that   |
   |                     | the default DART routine should just add noise to every element of state.                     |
   +---------------------+-----------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call get_close_maxdist_init(gc, maxdist)*
   ::

      type(get_close_type), intent(inout) :: gc
      real(r8),             intent(in)    :: maxdist

.. container:: indent1

   This is a PASS-THROUGH routine, the actual routine is the default one in ``location_mod``. In distance computations
   any two locations closer than the given ``maxdist`` will be considered close by the ``get_close_obs()`` routine.
   ``get_close_maxdist_init`` is listed on the ``use`` line for the locations_mod, and in the public list for this
   module, but has no subroutine declaration and no other code in this module.

| 

.. container:: routine

   *call get_close_obs_init(gc, num, obs)*
   ::

      type(get_close_type), intent(inout) :: gc
      integer,              intent(in)    :: num
      type(location_type),  intent(in)    :: obs(num)

.. container:: indent1

   This is a PASS-THROUGH routine. The default routine in the location module precomputes information to accelerate the
   distance computations done by ``get_close_obs()``. Like the other PASS-THROUGH ROUTINES it is listed on the use line
   for the locations_mod, and in the public list for this module, but has no subroutine declaration and no other code in
   this module:

| 

.. container:: routine

   *call get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, num_close, close_ind [, dist])*
   ::

      type(get_close_type), intent(in)  :: gc
      type(location_type),  intent(in)  :: base_obs_loc
      integer,              intent(in)  :: base_obs_kind
      type(location_type),  intent(in)  :: obs_loc(:)
      integer,              intent(in)  :: obs_kind(:)
      integer,              intent(out) :: num_close
      integer,              intent(out) :: close_ind(:)
      real(r8), optional,   intent(out) :: dist(:)

.. container:: indent1

   | Given a location and kind, compute the distances to all other locations in the ``obs_loc`` list. The return values
     are the number of items which are within maxdist of the base, the index numbers in the original obs_loc list, and
     optionally the distances. The ``gc`` contains precomputed information to speed the computations.
   | This is different than the default ``location_mod:get_close_obs()`` in that it is possible to modify the 'distance'
     based on the DART 'kind'. This allows one to apply specialized localizations.

   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``gc``            | The get_close_type which stores precomputed information about the locations to speed up         |
   |                   | searching                                                                                       |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``base_obs_loc``  | Reference location. The distances will be computed between this location and every other        |
   |                   | location in the obs list                                                                        |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``base_obs_kind`` | The kind of base_obs_loc                                                                        |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``obs_loc``       | Compute the distance between the base_obs_loc and each of the locations in this list            |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``obs_kind``      | The corresponding kind of each item in the obs list                                             |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``num_close``     | The number of items from the obs_loc list which are within maxdist of the base location         |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``close_ind``     | The list of index numbers from the obs_loc list which are within maxdist of the base location   |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | *dist*            | If present, return the distance between each entry in the close_ind list and the base location. |
   |                   | If not present, all items in the obs_loc list which are closer than maxdist will be added to    |
   |                   | the list but the overhead of computing the exact distances will be skipped.                     |
   +-------------------+-------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call ens_mean_for_model(ens_mean)*
   ::

      real(r8), dimension(:), intent(in) :: ens_mean

.. container:: indent1

   A model-size vector with the means of the ensembles for each of the state vector items. The model should save a local
   copy of this data if it needs to use it later to compute distances or other values. This routine is called after each
   model advance and contains the updated means.

   ============ ==========================================
   ``ens_mean`` State vector containing the ensemble mean.
   ============ ==========================================

TIEGCM public routines
~~~~~~~~~~~~~~~~~~~~~~

| 

.. container:: routine

   *call tiegcm_to_dart_vector(statevec, model_time)*
   ::

      real(r8), dimension(:), intent(out) :: statevec
      type(time_type),        intent(out) :: model_time

.. container:: indent1

   Read TIEGCM fields from the TIEGCM restart file and/or TIEGCM secondary file and pack them into a DART vector.

   ============== ================================================================
   ``statevec``   variable that contains the DART state vector
   ``model_time`` variable that contains the LAST TIME in the TIEGCM restart file.
   ============== ================================================================

| 

.. container:: routine

   *call dart_vector_to_tiegcm(statevec, dart_time)*
   ::

      real(r8), dimension(:), intent(in) :: statevec
      type(time_type),        intent(in) :: dart_time

.. container:: indent1

   | Unpacks a DART vector and updates the TIEGCM restart file variables. Only those variables designated as 'UPDATE'
     are put into the TIEGCM restart file. All variables are written to the DART diagnostic files **prior** to the
     application of any "clamping". The variables **are "clamped"** before being written to the TIEGCM restart file. The
     clamping limits are specified in columns 3 and 4 of ``&model_nml:variables``.
   | The time of the DART state is compared to the time in the restart file to ensure that we are not improperly
     updating a restart file.

   ============= ======================================================
   ``statevec``  Variable containing the DART state vector.
   ``dart_time`` Variable containing the time of the DART state vector.
   ============= ======================================================

| 

.. container:: routine

   *var = get_f107_value(x)*
   ::

      real(r8)                           :: get_f107_value
      real(r8), dimension(:), intent(in) :: x

.. container:: indent1

   If the F10_7 value is part of the DART state, return that value. If it is not part of the DART state, just return the
   F10_7 value from the TIEGCM namelist.

   ======= ==========================================
   ``x``   Variable containing the DART state vector.
   ``var`` The f10_7 value.
   ======= ==========================================

| 

.. container:: routine

   *call test_interpolate(x, locarray)*
   ::

      real(r8), dimension(:), intent(in) :: x
      real(r8), dimension(3), intent(in) :: locarray

.. container:: indent1

   This function is **only** used by
   :doc:`../../assimilation_code/programs/model_mod_check/model_mod_check.html%20models/POP/model_mod_check` and can be
   modified to suit your needs. ``test_interpolate()`` exercises ``model_interpolate()``, ``get_state_meta_data()``,
   ``static_init_model()`` and a host of supporting routines.

   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``x``                                                     | variable containing the DART state vector.                |
   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``locarray``                                              | variable containing the location of interest.             |
   |                                                           | locarray(1) is the longitude (in degrees East)            |
   |                                                           | locarray(2) is the latitude (in degrees North)            |
   |                                                           | locarray(3) is the height (in meters).                    |
   +-----------------------------------------------------------+-----------------------------------------------------------+

Files
-----

+--------------------------+------------------------------------------------------------------------------------------+
| ``filename``             | purpose                                                                                  |
+==========================+==========================================================================================+
| ``tiegcm.nml``           | TIEGCM control file modified to control starting and stopping.                           |
+--------------------------+------------------------------------------------------------------------------------------+
| ``input.nml``            | to read the model_mod namelist                                                           |
+--------------------------+------------------------------------------------------------------------------------------+
| ``tiegcm_restart_p.nc``  | both read and modified by the TIEGCM model_mod                                           |
+--------------------------+------------------------------------------------------------------------------------------+
| ``tiegcm_s.nc``          | read by the GCOM model_mod for metadata purposes.                                        |
+--------------------------+------------------------------------------------------------------------------------------+
| ``namelist_update``      | DART file containing information useful for starting and stopping TIEGCM.                |
|                          | ``advance_model.csh`` uses this to update the TIEGCM file ``tiegcm.nml``                 |
+--------------------------+------------------------------------------------------------------------------------------+
| ``dart_log.out``         | the run-time diagnostic output                                                           |
+--------------------------+------------------------------------------------------------------------------------------+
| ``dart_log.nml``         | the record of all the namelists (and their values) actually USED                         |
+--------------------------+------------------------------------------------------------------------------------------+
| *log_advance.\ nnnn.txt* | the run-time output of everything that happens in ``advance_model.csh``. This file will  |
|                          | be in the *advance_temp\ nnnn* directory.                                                |
+--------------------------+------------------------------------------------------------------------------------------+

References
----------

-  Matsuo, T., and E. A. Araujo-Pradere (2011),
   Role of thermosphere-ionosphere coupling in a global ionosphere specification,
   *Radio Science*, **46**, RS0D23, `doi:10.1029/2010RS004576 <http://dx.doi.org/doi:10.1029/2010RS004576>`__
-  
-  Lee, I. T., T, Matsuo, A. D. Richmond, J. Y. Liu, W. Wang, C. H. Lin, J. L. Anderson, and M. Q. Chen (2012),
   Assimilation of FORMOSAT-3/COSMIC electron density profiles into thermosphere/Ionosphere coupling model by using
   ensemble Kalman filter,
   *Journal of Geophysical Research*, **117**, A10318,
   `doi:10.1029/2012JA017700 <http://dx.doi.org/doi:10.1029/2012JA017700>`__
-  
-  Matsuo, T., I. T. Lee, and J. L. Anderson (2013),
   Thermospheric mass density specification using an ensemble Kalman filter,
   *Journal of Geophysical Research*, **118**, 1339-1350,
   `doi:10.1002/jgra.50162 <http://dx.doi.org/doi:10.1002/jgra.50162>`__
-  
-  Lee, I. T., H. F. Tsai, J. Y. Liu, Matsuo, T., and L. C. Chang (2013),
   Modeling impact of FORMOSAT-7/COSMIC-2 mission on ionospheric space weather monitoring,
   *Journal of Geophysical Research*, **118**, 6518-6523,
   `doi:10.1002/jgra.50538 <http://dx.doi.org/doi:10.1002/jgra.50538>`__
-  
-  Matsuo, T. (2014),
   Upper atmosphere data assimilation with an ensemble Kalman filter, in Modeling the Ionosphere-Thermosphere System,
   *Geophys. Monogr. Ser.*, vol. 201, edited by J. Huba, R. Schunk, and G. Khazanov, pp. 273-282, John Wiley & Sons,
   Ltd, Chichester, UK, `doi:10.1002/9781118704417 <http://dx.doi.org/doi:10.1002/9781118704417>`__
-  
-  Hsu, C.-H., T. Matsuo, W. Wang, and J. Y. Liu (2014),
   Effects of inferring unobserved thermospheric and ionospheric state variables by using an ensemble Kalman filter on
   global ionospheric specification and forecasting,
   *Journal of Geophysical Research*, **119**, 9256-9267,
   `doi:10.1002/2014JA020390 <http://dx.doi.org/doi:10.1002/2014JA020390>`__
-  
-  Chartier, A., T. Matsuo, J. L. Anderson, G. Lu, T. Hoar, N. Collins, A. Coster, C. Mitchell, L. Paxton, G. Bust
   (2015),
   Ionospheric Data Assimilation and Forecasting During Storms,
   *Journal of Geophysical Research*, under review
-  
