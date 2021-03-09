MITgcm_ocean
============

Overview
--------

| The `MIT ocean GCM <http://mitgcm.org/>`__ version 'checkpoint59a' is the foundation of this directory. It was
  modified by Ibrahim Hoteit of Scripps for his use, and so it differs from the original distribution.
| Since the model is highly parallelized, it can be compiled with a target number of processors in mind. From DART's
  perspective, the most logical strategy is to run ``filter`` or ``perfect_model_obs`` with **async=4**: advance the
  model in parallel ... one ensemble member after another. In this mode, the same set of processors are used for the
  data assimilation. The performance of the parallel assimilation algorithm has been tested up through 64 processors,
  and should scale well beyond that - but it remains to be quantified. The scaling for the ocean model is unknown to me,
  but Ibrahim routinely runs with many more than 64 processors.
| As for all DART experiments, the overall design for an experiment is this: the DART program ``filter`` will read the
  initial conditions file, the observation sequence file, and the DART namelist to decide whether or not to advance the
  ocean model. All of the control of the execution of the ocean model is done by DART directly. If the model needs to be
  advanced, ``filter`` makes a call to the shell to execute the script ``advance_model.csh``. ``advance_model.csh`` is
  ENTIRELY responsible for getting all the input files, data files, namelists, etc. into a temporary directory, running
  the model, and copying the results back to the parent directory (which we call CENTRALDIR). The whole process hinges
  on setting the ocean model namelist values such that it is doing a cold start for every model advance.

| 

Observations
^^^^^^^^^^^^

The observations for the ocean model were the first observations of oceanic quantities, so there is an
``observations/forward_operators/obs_def_MITgcm_ocean_mod.f90`` file containing the novel observation definitions like
*salinity, sea surface height, current components ...*. In keeping with the DART philosophy, there is a concept of
inheritance between platform-specific observations like *DRIFTER_U_CURRENT_COMPONENT* and the general
*U_CURRENT_COMPONENT*. Using the specific types when possible will allow flexibility specifying what kinds of
observations to assimilate. :doc:`./create_ocean_obs` is the program to create a DART observation sequence from a very
particular ASCII file.

| 

Converting between DART and the model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are a set of support programs:

+---------------------------+-----------------------------------------------------------------------------------------+
| :doc:`./trans_pv_sv`      | converts the ocean model snapshot files into a DART-compatible format                   |
+---------------------------+-----------------------------------------------------------------------------------------+
| :doc:`./trans_sv_pv`      | converts the DART output into snapshot files to be used as ocean model input datasets   |
|                           | (specified in ``data``\ ``&PARM05``); creates a new ``data`` namelist file              |
|                           | (``data.DART``) containing the correct ``&PARM03;startTime,endTime`` values to advance  |
|                           | the ocean model the expected amount; and creates a new ``data.cal`` namelist file       |
|                           | (``data.cal.DART``) containing the calendar information.                                |
+---------------------------+-----------------------------------------------------------------------------------------+
| :doc:`./create_ocean_obs` | create observation sequence files                                                       |
+---------------------------+-----------------------------------------------------------------------------------------+

The data assimilation period is controlled in the ``input.nml``\ ``&model_nml`` namelist. In combination with the ocean
model dynamics timestep ``data``\ ``&PARM03:deltaTClock`` this determines the amount of time the model will advance for
each assimilation cycle.

| 

Generating the initial ensemble
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| The MITgcm_ocean model cannot (as of Oct 2008) take one single model state and generate its own ensemble (typically
  done with pert_model_state). This means I don't really know how to perform a 'perfect model' experiment until I find a
  way to correctly perturb a single state to create an ensemble.
| The ensemble has to come from 'somewhere else'. I ran the model forward (outside the DART framework) for 14 days and
  output snapshot files ever 12 hours. One state vector can be generated from a set of snapshot files using
  ``trans_pv_sv``. I called this my 'initial ensemble' - it's better than nothing, but it is ENTIRELY unknown if this
  creates an intial ensemble with sufficient spread. Just for comparison, the initial ensemble for the atmospheric
  models is derived from 'climatological' values. If they need an 80-member ensemble for July 14, 2008; they use the
  July 1 estimates of the atmosphere from 1900 to 1979. By the time they assimilate (every 6 hours) for several days,
  things are on-track.
| There is a ``shell_scripts/MakeInitialEnsemble.csh`` script that was intended to automate this process - with modest
  success. It does illustrate the steps needed to convert each snapshot file to a DART initial conditions file and then
  run the `restart_file_utility <../../utilities/restart_file_utility.f90>`__ to overwrite the timestep in the header of
  the initial conditions file. After you have created all the initial conditions files, you can simply 'cat' them all
  together. Even if the script doesn't work *out-of-the-box*, it should be readable enough to be some help.

| 

Fortran direct-access big-endian data files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The MITgcm_ocean model uses Fortran direct-access big-endian data files. It is up to you to determine the proper
compiler flags to compile DART such that DART can read and write these files. Every compiler/architecture is different,
but we have put notes in each ``mkmf.template`` if we know how to achieve this.

| 

Controlling the model advances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| The assimilation period is specified by two namelist parameters in the ``input.nml``\ ``&model_nml`` namelist:
  ``assimilation_period_days`` and ``assimilation_period_seconds``. Normally, all observations within (+/-) HALF of the
  total assimilation period are used in the assimilation.
| The time of the initial conditions is specified by two namelist parameters in the ``input.nml``\ ``&model_nml``
  namelist: ``init_time_days`` and ``init_time_seconds``; depending on the settings of these parameters, the times may
  or may not come directly from the DART initial conditions files.
| The ocean model **MUST always** start from the input datasets defined in the ``data``\ ``&PARM05`` namelist.
  Apparently, this requires ``data``\ ``&PARM03:startTime`` to be **0.0**. One of the DART support routines
  (:doc:`./trans_sv_pv`) converts the DART state vector to the files used in ``data``\ ``&PARM05`` and creates new
  ``data.cal``\ ``&CAL_NML`` and ``data``\ ``&PARM03`` namelists with values appropriate to advance the model to the
  desired time.
| The ocean model then advances till ``data``\ ``&PARM03:endTime`` and writes out snapshot files. :doc:`./trans_pv_sv`
  converts the snapshot files to a DART-compatible file which is ingested by ``filter``. ``filter`` also reads the
  observation sequence file to determine which observations are within the assimilation window, assimilates them, and
  writes out a set of restart files, one for each ensemble member. ``filter`` then waits for each instance of the ocean
  model (one instance for each ensemble member) to advance to ``data``\ ``&PARM03:endTime``. The whole process repeats
  until 1) there are no more observations to assimilate (i.e. the observation sequence file is exhausted) or 2) the time
  specified by ``input.nml``\ ``&filter_nml:last_obs_days,last_obs_seconds`` has been reached.

| 

Getting started
^^^^^^^^^^^^^^^

I always like running something akin to a 'perfect model' experiment to start. Since I have not come up with a good way
to perturb a single model state to generate an ensemble, here's the next best thing. Please keep in mind that the
details for running each program are covered in their own documentation.

#. create a set of initial conditions for DART as described in Generating the intial ensemble and keep a copy of the
   'middle' snapshot - then use it as the initial condition for ``perfect_model_obs``.
#. create a TINY set of 'perfect' observations in the normal fashion:
   :doc:`../../assimilation_code/programs/create_obs_sequence/create_obs_sequence` and then
   :doc:`../../assimilation_code/programs/create_fixed_network_seq/create_fixed_network_seq` to create an empty
   observation sequence file (usually called ``obs_seq.in``)
#. modify ``data``, ``data.cal``, and ``input.nml`` to control the experiment and populate the observation sequence file
   by running :doc:`../../assimilation_code/programs/perfect_model_obs/perfect_model_obs`
#. Now use the full ensemble of initial conditions from Step 1 and run
   :doc:`../../assimilation_code/programs/filter/filter`

A perfectly sensible approach to get to know the system would be to try to

#. assimilate data for the first assimilation period and stop. Do not advance the model at all. The filter namelist can
   control all of this and you do not need to have a working ``advance_model.csh`` script, or even a working ocean model
   (as long as you have input data files).
#. advance the model first and then assimilate data for the first assimilation period and stop.
#. advance, assimilate and advance again. This tests the whole DART facility.

Exploring the output
^^^^^^^^^^^^^^^^^^^^

Is pretty much like any other model. The netCDF files have the model prognostic variables before and after the
assimilation. There are Matlab® scripts for perusing the netCDF files in the ``DART/matlab`` directory. There are
Matlab® scripts for exploring the performance of the assimilation in observation-space (after running
:doc:`../../assimilation_code/programs/obs_diag/threed_sphere/obs_diag` to explore the ``obs_seq.final`` file) - use the
scripts starting with ``'plot_'``, e.g. ``DART/diagnostics/matlab/plot_*.m``. As always, there are some model-specific
item you should know about in ``DART/models/MITgcm_ocean/matlab``, and ``DART/models/MITgcm_ocean/shell_scripts``.

Other modules used
------------------

::

   types_mod
   time_manager_mod
   threed_sphere/location_mod
   utilities_mod
   obs_kind_mod
   mpi_utilities_mod
   random_seq_mod

Public interfaces
-----------------

Only a select number of interfaces used are discussed here.

========================== ================================================================================
*use location_mod, only :* `location_type <../../location/threed_sphere/location_mod.html#location_type>`__
\                          `get_location <../../location/threed_sphere/location_mod.html#get_location>`__
\                          `set_location <../../location/threed_sphere/location_mod.html#set_location>`__
========================== ================================================================================

The ocean model namelists ``data``, and ``data.cal`` *MUST* be present. These namelists are needed to reconstruct the
valid time of the snapshot files created by the ocean model. Be aware that as DART advances the model, the ``data``
namelist gets modified to reflect the current time of the model output.

Required Interface Routines

*use model_mod, only :*

get_model_size

adv_1step

get_state_meta_data

model_interpolate

get_model_time_step

static_init_model

end_model

init_time

init_conditions

nc_write_model_atts

nc_write_model_vars

pert_model_state

get_close_maxdist_init

get_close_obs_init

get_close_obs

ens_mean_for_model

Unique Interface Routines

*use model_mod, only :*

MIT_meta_type

read_meta

write_meta

prog_var_to_vector

vector_to_prog_var

read_snapshot

write_snapshot

get_gridsize

snapshot_files_to_sv

sv_to_snapshot_files

timestep_to_DARTtime

DARTtime_to_MITtime

DARTtime_to_timestepindex

write_data_namelistfile

Ocean model namelist interfaces ``&PARM03``, ``&PARM04``, and ``&PARM04`` are read from file ``data``. Ocean model
namelist interface ``&CAL_NML``, is read from file ``data.cal``.

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

   ``adv_1step`` is not used for the MITgcm_ocean model. Advancing the model is done through the ``advance_model``
   script. This is a NULL_INTERFACE, provided only for compatibility with the DART requirements.

   ======== ==========================================
   ``x``    State vector of length model_size.
   ``time`` Specifies time of the initial model state.
   ======== ==========================================

| 

.. container:: routine

   *call get_state_meta_data (index_in, location, [, var_type] )*
   ::

      integer,             intent(in)  :: index_in
      type(location_type), intent(out) :: location
      integer, optional,   intent(out) ::  var_type 

.. container:: indent1

   ``get_state_meta_data`` returns metadata about a given element of the DART representation of the model state vector.
   Since the DART model state vector is a 1D array and the native model grid is multidimensional,
   ``get_state_meta_data`` returns information about the native model state vector representation. Things like the
   ``location``, or the type of the variable (for instance: salinity, temperature, u current component, ...). The
   integer values used to indicate different variable types in ``var_type`` are themselves defined as public interfaces
   to model_mod if required.

   +--------------+------------------------------------------------------------------------------------------------------+
   | ``index_in`` | Index of state vector element about which information is requested.                                  |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``location`` | Returns the 3D location of the indexed state variable. The ``location_ type`` comes from             |
   |              | ``DART/location/threed_sphere/location_mod.f90``. Note that the lat/lon are specified in degrees by  |
   |              | the user but are converted to radians internally.                                                    |
   +--------------+------------------------------------------------------------------------------------------------------+
   | *var_type*   | Returns the type of the indexed state variable as an optional argument. The type is one of the list  |
   |              | of supported observation types, found in the block of code starting                                  |
   |              | ``! Integer definitions for DART TYPES`` in                                                          |
   |              | ``DART/assimilation_code/modules/observations/obs_kind_mod.f90``                                     |
   +--------------+------------------------------------------------------------------------------------------------------+

   The list of supported variables in ``DART/assimilation_code/modules/observations/obs_kind_mod.f90`` is created by
   ``preprocess`` using the entries in ``input.nml``\ [``&preprocess_nml, &obs_kind_nml``], ``DEFAULT_obs_kin_mod.F90``
   and ``obs_def_MITgcm_ocean_mod.f90``.

| 

.. container:: routine

   *call model_interpolate(x, location, itype, obs_val, istatus)*
   ::

      real(r8), dimension(:), intent(in)  :: x
      type(location_type),    intent(in)  :: location
      integer,                intent(in)  :: itype
      real(r8),               intent(out) :: obs_val
      integer,                intent(out) :: istatus

.. container:: indent1

   | Given a model state, ``model_interpolate`` returns the value of the desired observation type (which could be a
     state variable) that would be observed at the desired location. The interpolation method is either completely
     specified by the model, or uses some standard 2D or 3D scalar interpolation routines. Put another way,
     ``model_interpolate`` will apply the forward operator **H** to the model state to create an observation at the
     desired location.
   | If the interpolation is valid, ``istatus = 0``. In the case where the observation operator is not defined at the
     given location (e.g. the observation is below the lowest model level, above the top level, or 'dry'), interp_val is
     returned as 0.0 and istatus = 1.

   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``x``                                                     | A model state vector.                                     |
   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``location``                                              | Location to which to interpolate.                         |
   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``itype``                                                 | Not used.                                                 |
   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``obs_val``                                               | The interpolated value from the model.                    |
   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``istatus``                                               | Integer flag indicating the success of the interpolation. |
   |                                                           | success == 0, failure == anything else                    |
   +-----------------------------------------------------------+-----------------------------------------------------------+

| 

.. container:: routine

   *var = get_model_time_step()*
   ::

      type(time_type) :: get_model_time_step

.. container:: indent1

   ``get_model_time_step`` returns the forecast length to be used as the "model base time step" in the filter. This is
   the minimum amount of time the model can be advanced by ``filter``. *This is also the assimilation window*. All
   observations within (+/-) one half of the forecast length are used for the assimilation. In the ``MITgcm_ocean``
   case, this is set from the namelist values for
   ``input.nml``\ ``&model_nml:assimilation_period_days, assimilation_period_seconds``, after ensuring the forecast
   length is a multiple of the ocean model dynamical timestep declared by ``data``\ ``&PARM03:deltaTClock``.

   ======= ============================
   ``var`` Smallest time step of model.
   ======= ============================

   Please read the note concerning Controlling the model advances

| 

.. container:: routine

   *call static_init_model()*

.. container:: indent1

   | ``static_init_model`` is called for runtime initialization of the model. The namelists are read to determine
     runtime configuration of the model, the calendar information, the grid coordinates, etc. There are no input
     arguments and no return values. The routine sets module-local private attributes that can then be queried by the
     public interface routines.
   | The namelists (all mandatory) are:
   | ``input.nml``\ ``&model_mod_nml``,
   | ``data.cal``\ ``&CAL_NML``,
   | ``data``\ ``&PARM03``,
   | ``data``\ ``&PARM04``, and
   | ``data``\ ``&PARM05``.

| 

.. container:: routine

   *call end_model()*

.. container:: indent1

   ``end_model`` is used to clean up storage for the model, etc. when the model is no longer needed. There are no
   arguments and no return values. This is required by DART but nothing needs to be done for the MITgcm_ocean model.

| 

.. container:: routine

   *call init_time(time)*
   ::

      type(time_type), intent(out) :: time

.. container:: indent1

   ``init_time`` returns the time at which the model will start if no input initial conditions are to be used. This is
   frequently used to spin-up models from rest, but is not meaningfully supported for the MITgcm_ocean model. The only
   time this routine would get called is if the ``input.nml``\ ``&perfect_model_obs_nml:start_from_restart`` is .false.,
   which is not supported in the MITgcm_ocean model.

   +----------+----------------------------------------------------------------------------------------------------------+
   | ``time`` | the starting time for the model if no initial conditions are to be supplied. As of Oct 2008, this is     |
   |          | hardwired to 0.0                                                                                         |
   +----------+----------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call init_conditions(x)*
   ::

      real(r8), dimension(:), intent(out) :: x

.. container:: indent1

   ``init_conditions`` returns default initial conditions for model; generally used for spinning up initial model
   states. For the MITgcm_ocean model it is just a stub because the initial state is always provided by the input files.

   ===== ==========================================================================
   ``x`` Model state vector. [default is 0.0 for every element of the state vector]
   ===== ==========================================================================

| 

.. container:: routine

   *ierr = nc_write_model_atts(ncFileID)*
   ::

      integer             :: nc_write_model_atts
      integer, intent(in) :: ncFileID

.. container:: indent1

   ``nc_write_model_atts`` writes model-specific attributes to an opened netCDF file: In the MITgcm_ocean case, this
   includes information like the coordinate variables (the grid arrays: XG, XC, YG, YC, ZG, ZC, ...), information from
   some of the namelists, and either the 1D state vector or the prognostic variables (S,T,U,V,Eta). All the required
   information (except for the netCDF file identifier) is obtained from the scope of the ``model_mod`` module.

   ============ =========================================================
   ``ncFileID`` Integer file descriptor to previously-opened netCDF file.
   ``ierr``     Returns a 0 for successful completion.
   ============ =========================================================

   ``nc_write_model_atts`` is responsible for the model-specific attributes in the following DART-output netCDF files:
   ``true_state.nc``, ``preassim.nc``, and ``analysis.nc``.

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

   ``nc_write_model_vars`` writes a copy of the state variables to a NetCDF file. Multiple copies of the state for a
   given time are supported, allowing, for instance, a single file to include multiple ensemble estimates of the state.
   Whether the state vector is parsed into prognostic variables (S,T,U,V,Eta) or simply written as a 1D array is
   controlled by ``input.nml``\ ``&model_mod_nml:output_state_vector``. If ``output_state_vector = .true.`` the state
   vector is written as a 1D array (the simplest case, but hard to explore with the diagnostics). If
   ``output_state_vector = .false.`` the state vector is parsed into prognostic variables before being written.

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

   | Given a model state, ``pert_model_state`` produces a perturbed model state. This is used to generate ensemble
     initial conditions perturbed around some control trajectory state when one is preparing to spin-up ensembles. Since
     the DART state vector for the MITgcm_ocean model contains both 'wet' and 'dry' cells, (the 'dry' cells having a
     value of a perfect 0.0 - not my choice) it is imperative to provide an interface to perturb **just** the wet cells
     (``interf_provided == .true.``).
   | At present (Oct 2008) the magnitude of the perturbation is wholly determined by
     ``input.nml``\ ``&model_mod_nml:model_perturbation_amplitude`` and **utterly, completely fails**. The resulting
     model states cause a fatal error when being read in by the ocean model - something like

   ::

      *** ERROR *** S/R INI_THETA: theta = 0 identically. 
      If this is intentional you will need to edit ini_theta.F to avoid this safety check

   A more robust perturbation mechanism is needed (see, for example this routine in the CAM model_mod.f90). Until then,
   you can avoid using this routine by using your own ensemble of initial conditions. This is determined by setting
   ``input.nml``\ ``&filter_nml:start_from_restart = .false.`` See also Generating the initial ensemble at the start of
   this document.

   +---------------------+-----------------------------------------------------------------------------------------------+
   | ``state``           | State vector to be perturbed.                                                                 |
   +---------------------+-----------------------------------------------------------------------------------------------+
   | ``pert_state``      | The perturbed state vector.                                                                   |
   +---------------------+-----------------------------------------------------------------------------------------------+
   | ``interf_provided`` | Because of the 'wet/dry' issue discussed above, this is always ``.true.``, indicating a       |
   |                     | model-specific perturbation is available.                                                     |
   +---------------------+-----------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call get_close_maxdist_init(gc, maxdist)*
   ::

      type(get_close_type), intent(inout) :: gc
      real(r8),             intent(in)    :: maxdist

.. container:: indent1

   Pass-through to the 3-D sphere locations module. See
   `get_close_maxdist_init() <../../location/threed_sphere/location_mod.html#get_close_maxdist_init>`__ for the
   documentation of this subroutine.

| 

.. container:: routine

   *call get_close_obs_init(gc, num, obs)*
   ::

      type(get_close_type), intent(inout) :: gc
      integer,              intent(in)    :: num
      type(location_type),  intent(in)    :: obs(num)

.. container:: indent1

   Pass-through to the 3-D sphere locations module. See
   `get_close_obs_init() <../../location/threed_sphere/location_mod.html#get_close_obs_init>`__ for the documentation of
   this subroutine.

| 

.. container:: routine

   *call get_close_obs(gc, base_obs_loc, base_obs_kind, obs, obs_kind, num_close, close_ind [, dist])*
   ::

      type(get_close_type), intent(in)  :: gc
      type(location_type),  intent(in)  :: base_obs_loc
      integer,              intent(in)  :: base_obs_kind
      type(location_type),  intent(in)  :: obs(:)
      integer,              intent(in)  :: obs_kind(:)
      integer,              intent(out) :: num_close
      integer,              intent(out) :: close_ind(:)
      real(r8), optional,   intent(out) :: dist(:)

.. container:: indent1

   Pass-through to the 3-D sphere locations module. See
   `get_close_obs() <../../location/threed_sphere/location_mod.html#get_close_obs>`__ for the documentation of this
   subroutine.

| 

.. container:: routine

   *call ens_mean_for_model(ens_mean)*
   ::

      real(r8), dimension(:), intent(in) :: ens_mean

.. container:: indent1

   ``ens_mean_for_model`` saves a copy of the ensemble mean to module-local storage. Sometimes the ensemble mean is
   needed rather than individual copy estimates. This is a NULL_INTERFACE for the MITgcm_ocean model. At present there
   is no application which requires module-local storage of the ensemble mean. No storage is allocated.

   ============ ==========================
   ``ens_mean`` Ensemble mean state vector
   ============ ==========================

| 

Unique interface routines
-------------------------

| 

.. container:: type

   ::

      type MIT_meta_type
         private
         integer           :: nDims
         integer           :: dimList(3)
         character(len=32) :: dataprec
         integer           :: reclen
         integer           :: nrecords
         integer           :: timeStepNumber
      end type MIT_meta_type

.. container:: indent1

   ``MIT_meta_type`` is a derived type used to codify the metadata associated with a snapshot file.

   +----------------+----------------------------------------------------------------------------------------------------+
   | Component      | Description                                                                                        |
   +================+====================================================================================================+
   | nDims          | the number of dimensions for the associated object. S,T,U,V all have nDims==3, Eta has nDims==2    |
   +----------------+----------------------------------------------------------------------------------------------------+
   | dimList        | the extent of each of the dimensions                                                               |
   +----------------+----------------------------------------------------------------------------------------------------+
   | dataprec       | a character string depicting the precision of the data storage. Commonly 'float32'                 |
   +----------------+----------------------------------------------------------------------------------------------------+
   | reclen         | the record length needed to correctly read using Fortran direct-access. This is tricky business.   |
   |                | Each vendor has their own units for record length. Sometimes it is bytes, sometimes words,         |
   |                | sometimes ???. See comments in code for ``item_size_direct_access``                                |
   +----------------+----------------------------------------------------------------------------------------------------+
   | nrecords       | the number of records (either 2D or 3D hyperslabs) in the snapshot file                            |
   +----------------+----------------------------------------------------------------------------------------------------+
   | timeStepNumber | the timestep number ... the snapshot filenames are constructed using the timestepcount as the      |
   |                | unique part of the filename. To determine the valid time of the snapshot, you must multiply the    |
   |                | timeStepNumber by the amount of time in each timestep and add the start time.                      |
   +----------------+----------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *metadata = read_meta(fbase [, vartype])*
   ::

      character(len=*),           intent(in)  ::  fbase 
      character(len=*), OPTIONAL, intent(in)  ::  vartype 
      type(MIT_meta_type),        intent(out) ::  metadata 

.. container:: indent1

   | ``read_meta`` reads the metadata file for a particular snapshot file. This routine is primarily bulletproofing,
     since the snapshot files tend to move around a lot. I don't want to use a snapshot file from a 70-level case in a
     40-level experiment; and without checking the metadata, you'd never know. The metadata for the file originally
     comes from the namelist values specifying the grid resolution, etc. If the metadata file exists, the metadata in
     the file is compared to the original specifications. If the metadata file does not exist, no comparison is done.
   | The filename is fundamentally comprised of three parts. Take 'U.0000000024.meta' for example. The first part of the
     name is the variable, the second part of the name is the timestepnumber, the last part is the file extension. For
     various reasons, sometimes it is convenient to call this function without the building the entire filename outside
     the function and then passing it in as an argument. Since the '.meta' extension seems to be fixed, we will only
     concern ourselves with building the 'base' part of the filename, i.e., the first two parts.

   +--------------+------------------------------------------------------------------------------------------------------+
   | ``fbase``    | If *vartype* is supplied, this is simply the timestepnumber converted to a character string of       |
   |              | length 10. For example, '0000000024'. If *vartype* is **not** supplied, it is the entire filename    |
   |              | without the extension; 'U.0000000024', for example.                                                  |
   +--------------+------------------------------------------------------------------------------------------------------+
   | *vartype*    | is an optional argument specifying the first part of the snapshot filename. Generally,               |
   |              | 'S','T','U','V', or 'Eta'.                                                                           |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``metadata`` | The return value of the function is the metadata for the file, packed into a user-derived variable   |
   |              | type specifically designed for the purpose.                                                          |
   +--------------+------------------------------------------------------------------------------------------------------+

   .. rubric:: Metadata example
      :name: metadata-example
      :class: indent1

   ::

      metadata = read_meta('U.0000000024')
       ... or ...
      metadata = read_meta('0000000024','U')

| 

.. container:: routine

   *call write_meta(metadata, filebase)*
   ::

      type(MIT_meta_type),        intent(in) ::  metadata 
      character(len=*),           intent(in) ::  filebase 

.. container:: indent1

   ``write_meta`` writes a metadata file. This routine is called by routines ``write_2d_snapshot``, and
   ``write_3d_snapshot`` to support converting the DART state vector to something the ocean model can ingest.

   ============ =======================================================================================================
   ``metadata`` The user-derived varible, filled with the metadata for the file.
   ``filebase`` the filename without the extension; 'U.0000000024', for example. (see the Description in ``read_meta``)
   ============ =======================================================================================================

| 

.. container:: routine

   *call prog_var_to_vector(s,t,u,v,eta,x)*
   ::

      real(r4), dimension(:,:,:), intent(in)  :: s,t,u,v
      real(r4), dimension(:,:),   intent(in)  :: eta
      real(r8), dimension(:),     intent(out) :: x

.. container:: indent1

   ``prog_var_to_vector`` packs the prognostic variables [S,T,U,V,Eta] read from the snapshot files into a DART vector.
   The DART vector is simply a 1D vector that includes all the 'dry' cells as well as the 'wet' ones. This routine is
   not presently used (since we never have [S,T,U,V,Eta] as such in memory). See snapshot_files_to_sv.

   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``s,t,u,v`` | The 3D arrays read from the individual snapshot files.                                                |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``eta``     | The 2D array read from its snapshot file.                                                             |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``x``       | the 1D array containing the concatenated s,t,u,v,eta variables. To save storage, it is possible to    |
   |             | modify the definition of ``r8`` in ``DART/common/types_mod.f90`` to be the same as that of ``r4``.    |
   +-------------+-------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call vector_to_prog_var(x,varindex,hyperslab)*
   ::

      real(r8), dimension(:),     intent(in)  :: x
      integer,                    intent(in)  :: varindex
      real(r4), dimension(:,:,:), intent(out) :: hyperslab -or-
      real(r4), dimension(:,:),   intent(out) :: hyperslab

.. container:: indent1

   ``vector_to_prog_var`` unpacks a prognostic variable [S,T,U,V,Eta] from the DART vector ``x``.

   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``x``                                                     | the 1D array containing the 1D DART state vector.         |
   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``varindex``                                              | an integer code specifying which variable to unpack. The  |
   |                                                           | following parameters are in module storage:               |
   |                                                           | ::                                                        |
   |                                                           |                                                           |
   |                                                           |    integer, parameter :: S_index   = 1                    |
   |                                                           |    integer, parameter :: T_index   = 2                    |
   |                                                           |    integer, parameter :: U_index   = 3                    |
   |                                                           |    integer, parameter :: V_index   = 4                    |
   |                                                           |    integer, parameter :: Eta_index = 5                    |
   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``hyperslab``                                             | The N-D array containing the prognostic variable. The     |
   |                                                           | function is overloaded to be able to return both 2D and   |
   |                                                           | 3D arrays.                                                |
   +-----------------------------------------------------------+-----------------------------------------------------------+

   .. rubric:: Vector_to_prog_var
      :name: vector_to_prog_var
      :class: indent1

   ::

      call vector_to_prog_var(statevec,V_index,data_3d)
       - or - 
      call vector_to_prog_var(statevec,Eta_index,data_2d)

| 

.. container:: routine

   *call read_snapshot(fbase, x, timestep, vartype)*
   ::

      character(len=*),           intent(in)  :: fbase
      real(r4), dimension(:,:,:), intent(out) :: x - or - 
      real(r4), dimension(:,:),   intent(out) :: x
      integer,                    intent(out) :: timestep
      character(len=*), optional, intent(in)  :: vartype

.. container:: indent1

   ``read_snapshot`` reads a snapshot file and returns a hyperslab that includes all the 'dry' cells as well as the
   'wet' ones. By design, the MITgcm_ocean model writes out Fortran direct-access big-endian binary files, independent
   of the platform. Since it is not guaranteed that the binary file we need to read is on the same architecture that
   created the file, getting the compiler settings in ``mkmf.template`` correct to read Fortran direct-access big-endian
   binary files is **imperative** to the process. Since each compiler issues its own error, there's no good way to even
   summarize the error messages you are likely to encounter by improperly reading the binary files. Read each template
   file for hints about the proper settings. See also the section Fortran direct-access big-endian datafiles in the
   "Discussion" of this document.

   +--------------+------------------------------------------------------------------------------------------------------+
   | ``fbase``    | The 'base' portion of the filename, i.e., without the [.meta, .data] extension. If *vartype* is      |
   |              | supplied, *vartype* is prepended to ``fbase`` to create the 'base' portion of the filename.          |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``x``        | The hyperslab containing what is read. The function is overloaded to be able to return a 2D or 3D    |
   |              | array. ``x`` must be allocated before the call to ``read_snapshot``.                                 |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``timestep`` | the timestepcount in the ``'fbase'``.meta file, if the .meta file exists. Provided for               |
   |              | bulletproofing.                                                                                      |
   +--------------+------------------------------------------------------------------------------------------------------+
   | *vartype*    | The character string representing the 'prognostic variable' portion of the snapshot filename.        |
   |              | Commonly 'S','T','U','V', or 'Eta'. If supplied, this is prepended to ``fbase`` to create the 'base' |
   |              | portion of the filename.                                                                             |
   +--------------+------------------------------------------------------------------------------------------------------+

   .. rubric:: Code snippet
      :name: code-snippet

   ::

      real(r4), allocatable :: data_2d_array(:,:), data_3d_array(:,:,:)
      ...
      allocate(data_2d_array(Nx,Ny), data_3d_array(Nx,Ny,Nz))
      ...
      call read_snapshot('S.0000000024', data_3d_array, timestepcount_out)
      call read_snapshot(  '0000000024', data_2d_array, timestepcount_out, 'Eta')
      call read_snapshot(  '0000000024', data_3d_array, timestepcount_out, 'T')
      ...

| 

.. container:: routine

   *call write_snapshot(x, fbase, timestepcount)*
   ::

      real(r4), dimension(:,:),   intent(in) :: x - or -
      real(r4), dimension(:,:,:), intent(in) :: x
      character(len=*),           intent(in) :: fbase
      integer, optional,          intent(in) :: timestepcount

.. container:: indent1

   ``write_snapshot`` writes a hyperslab of data to a snapshot file and corresponding metadata file. This routine is an
   integral part of sv_to_snapshot_files, the routine that is responsible for unpacking the DART state vector and
   writing out a set of snapshot files used as input to the ocean model.

   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``x``             | The hyperslab containing the prognostic variable data to be written. The function is overloaded |
   |                   | to be able to ingest a 2D or 3D array.                                                          |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``fbase``         | The 'base' portion of the filename, i.e., without the [.meta, .data] extension.                 |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``timestepcount`` | the timestepcount to be written into the ``'fbase'``.meta file. If none is supplied,            |
   |                   | ``timestepcount`` is 0. I'm not sure this is ever used, since the timestepcount can be gotten   |
   |                   | from ``fbase``.                                                                                 |
   +-------------------+-------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call get_gridsize( num_x, num_y, num_z)*
   ::

      integer, intent(out) :: num_x, num_y, num_z

.. container:: indent1

   ``get_gridsize`` returns the dimensions of the compute domain. The gridsize is determined from
   ``data``\ ``&PARM04:delY,delX``, and ``delZ`` when the namelist is read by ``static_init_model``. The MITgcm_ocean
   model is interesting in that it has a staggered grid but all grid variables are declared the same length.

   ========= ======================================
   ``num_x`` The number of longitudinal gridpoints.
   ``num_y`` The number of latitudinal gridpoints.
   ``num_z`` The number of vertical gridpoints.
   ========= ======================================

| 

.. container:: routine

   *call snapshot_files_to_sv(timestepcount, state_vector)*
   ::

      integer,  intent(in)    :: timestepcount
      real(r8), intent(inout) :: state_vector

.. container:: indent1

   ``snapshot_files_to_sv`` reads the snapshot files for a given timestepcount and concatenates them into a
   DART-compliant 1D array. All the snapshot filenames are constructed given the ``timestepcount`` - read the
   'Description' section of read_meta, particularly the second paragraph.

   ================= ============================================================================
   ``timestepcount`` The integer that corresponds to the middle portion of the snapshot filename.
   ``state_vector``  The 1D array of the DART state vector.
   ================= ============================================================================

   The files are read in this order [S,T,U,V,Eta] (almost alphabetical!) and the multidimensional arrays are unwrapped
   with the leftmost index being the fastest-varying. You shouldn't need to know this, but it is critical to the way
   ``prog_var_to_vector`` and ``vector_to_prog_var`` navigate the array.

   ::

      do k = 1, Nz   ! depth
      do j = 1, Ny   ! latitudes
      do i = 1, Nx   ! longitudes
         state_vector(indx) = data_3d_array(i, j, k)
         indx = indx + 1
      enddo
      enddo
      enddo

| 

.. container:: routine

   *call sv_to_snapshot_files(state_vector, date1, date2)*
   ::

      real(r8), intent(in)    :: state_vector
      integer,  intent(in)    :: date1, date2

.. container:: indent1

   ``sv_to_snapshot_files`` takes the DART state vector and creates a set of snapshot files. The filenames of these
   snapshot files is different than that of snapshot files created by the ocean model. See the 'Notes' section for an
   explanation.

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``state_vector`` | The DART 1D state vector.                                                                        |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``date1``        | The year/month/day of the valid time for the state vector, in YYYYMMDD format - an 8-digit       |
   |                  | integer. This is the same format as ``data.cal``\ ``&CAL_NML:startDate_1``                       |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``date2``        | The hour/min/sec of the valid time for the state vector, in HHMMSS format. This is the same      |
   |                  | format as ``data.cal``\ ``&CAL_NML:startDate_2``                                                 |
   +------------------+--------------------------------------------------------------------------------------------------+

   Since the snapshot files have the potential to move around a lot, I thought it best to have a more descriptive name
   than simply the snapshot number. DART creates snapshot files with names like ``S.19960718.060000.data`` to let you
   know it is a snapshot file for 06Z 18 July 1996. This is intended to make it easier to create initial conditions
   files and, should the assimilation fail, inform as to \_when\_ the assimilation failed. Since DART needs the ocean
   model to coldstart (``data``\ ``&PARM02:startTime = 0.0``) for every model advance, every snapshot file has the same
   timestamp. The ``advance_model.csh`` script actually has to rename the DART-written snapshot files to that declared
   by the ``data``\ ``&PARM05`` namelist, so the name is not really critical from that perspective. **However**, the
   components of the DART-derived snapshot files **are** used to create an appropriate ``data.cal``\ ``&CAL_NML`` for
   each successive model advance.

| 

.. container:: routine

   *mytime = timestep_to_DARTtime(TimeStepIndex)*
   ::

      integer,         intent(in)  :: TimeStepIndex
      type(time_type), intent(out) :: mytime

.. container:: indent1

   ``timestep_to_DARTtime`` combines the ``TimeStepIndex`` with the time per timestep (from ``data``\ ``&PARM03``) and
   the start date supplied by ``data.cal``\ ``&CAL_NML`` to form a Gregorian calendar date which is then converted to a
   DART time object. As of Oct 2008, this ``model_mod`` is forced to use the Gregorian calendar.

   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``TimeStepIndex`` | an integer referring to the ocean model timestep ... the middle part of the ocean-model-flavor  |
   |                   | snapshot filename.                                                                              |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``mytime``        | The DART representation of the time indicated by the ``TimeStepIndex``                          |
   +-------------------+-------------------------------------------------------------------------------------------------+

   The time per timestep is something I don't understand that well. The ``data``\ ``&PARM03`` namelist has three
   variables: ``deltaTmom``, ``deltaTtracer``, and ``deltaTClock``. Since I don't know which one is relavent, and every
   case I looked at had them set to be the same, I decided to require that they all be identical and then it wouldn't
   matter which one I used. The values are checked when the namelist is read.

   ::

      ! Time stepping parameters are in PARM03
      call find_namelist_in_file("data", "PARM03", iunit)
      read(iunit, nml = PARM03, iostat = io)
      call check_namelist_read(iunit, io, "PARM03")

      if ((deltaTmom   == deltaTtracer) .and. &
          (deltaTmom   == deltaTClock ) .and. &
          (deltaTClock == deltaTtracer)) then
         timestep       = deltaTmom                    ! need a time_type version
      else
         write(msgstring,*)"namelist PARM03 has deltaTmom /= deltaTtracer /= deltaTClock"
         call error_handler(E_MSG,"static_init_model", msgstring, source, revision, revdate)
         write(msgstring,*)"values were ",deltaTmom, deltaTtracer, deltaTClock
         call error_handler(E_MSG,"static_init_model", msgstring, source, revision, revdate)
         write(msgstring,*)"At present, DART only supports equal values."
         call error_handler(E_ERR,"static_init_model", msgstring, source, revision, revdate)
      endif

| 

.. container:: routine

   *call DARTtime_to_MITtime(darttime, date1, date2)*
   ::

      type(time_type), intent(in)  :: darttime
      integer,         intent(out) :: date1, date2

.. container:: indent1

   ``DARTtime_to_MITtime`` converts the DART time to a pair of integers that are compatible with the format used in
   ``data.cal``\ ``&CAL_NML``

   +--------------+------------------------------------------------------------------------------------------------------+
   | ``darttime`` | The DART time to be converted.                                                                       |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``date1``    | The year/month/day component of the time in YYYYMMDD format - an 8-digit integer. This is the same   |
   |              | format as ``data.cal``\ ``&CAL_NML:startDate_1``                                                     |
   +--------------+------------------------------------------------------------------------------------------------------+
   | ``date2``    | The hour/min/sec component of the time in HHMMSS format. This is the same format as                  |
   |              | ``data.cal``\ ``&CAL_NML:startDate_2``                                                               |
   +--------------+------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *timeindex = DARTtime_to_timestepindex(darttime)*
   ::

      type(time_type), intent(in)  :: darttime
      integer,         intent(out) :: timeindex

.. container:: indent1

   ``DARTtime_to_timestepindex`` converts the DART time to an integer representing the number of timesteps since the
   date in ``data.cal``\ ``&CAL_NML``, i.e., the start of the model run. The size of each timestep is determined as
   discussed in the timestep_to_DARTtime section.

   ============= =========================================================
   ``darttime``  The DART time to be converted.
   ``timeindex`` The number of timesteps corresponding to the DARTtime ...
   ============= =========================================================

| 

.. container:: routine

   *call write_data_namelistfile()*

.. container:: indent1

   | There are no input arguments to ``write_data_namelistfile``. ``write_data_namelistfile`` reads the ``data``
     namelist file and creates an almost-identical copy named ``data.DART`` that differs only in the namelist parameters
     that control the model advance.
   | (NOTE) ``advance_model.csh`` is designed to first run ``trans_sv_pv`` to create appropriate ``data.DART`` and
     ``data.cal.DART`` files. The script then renames them to that expected by the ocean model.

| 

Namelists
---------

We adhere to the F90 standard of starting a namelist with an ampersand '&' and terminating with a slash '/' for all our
namelist input. Consider yourself forewarned that character strings that contain a '/' must be enclosed in quotes to
prevent them from prematurely terminating the namelist.

.. container:: namelist

   ::

      namelist /model_nml/  assimilation_period_days, &
           assimilation_period_seconds, output_state_vector, model_perturbation_amplitude

.. container:: indent1

   This namelist is read in a file called ``input.nml``. This namelist provides control over the assimilation period for
   the model. All observations within (+/-) half of the assimilation period are assimilated. The assimilation period is
   the minimum amount of time the model can be advanced, and checks are performed to ensure that the assimilation window
   is a multiple of the ocean model dynamical timestep indicated by ``PARM03:deltaTClock``.

   +------------------------------+-----------------------------+-------------------------------------------------------+
   | Contents                     | Type                        | Description                                           |
   +==============================+=============================+=======================================================+
   | assimilation_period_days     | integer *[default: 7]*      | The number of days to advance the model for each      |
   |                              |                             | assimilation.                                         |
   +------------------------------+-----------------------------+-------------------------------------------------------+
   | assimilation_period_seconds  | integer *[default: 0]*      | In addition to ``assimilation_period_days``, the      |
   |                              |                             | number of seconds to advance the model for each       |
   |                              |                             | assimilation.                                         |
   +------------------------------+-----------------------------+-------------------------------------------------------+
   | output_state_vector          | logical *[default: .true.]* | The switch to determine the form of the state vector  |
   |                              |                             | in the output netcdf files. If ``.true.`` the state   |
   |                              |                             | vector will be output exactly as DART uses it ... one |
   |                              |                             | long array. If ``.false.``, the state vector is       |
   |                              |                             | parsed into prognostic variables and output that way  |
   |                              |                             | -- much easier to use with 'ncview', for example.     |
   +------------------------------+-----------------------------+-------------------------------------------------------+
   | model_perturbation_amplitude | real(r8) *[default: 0.2]*   | The amount of noise to add when trying to perturb a   |
   |                              |                             | single state vector to create an ensemble. Only       |
   |                              |                             | needed when                                           |
   |                              |                             | ``inpu                                                |
   |                              |                             | t.nml``\ ``&filter_nml:start_from_restart = .false.`` |
   |                              |                             | See also Generating the initial ensemble at the start |
   |                              |                             | of this document. units: standard deviation of a      |
   |                              |                             | gaussian distribution with the mean at the value of   |
   |                              |                             | the state vector element.                             |
   +------------------------------+-----------------------------+-------------------------------------------------------+

   .. rubric:: Model namelist
      :name: model-namelist

   ::

      &model_nml
         assimilation_period_days     = 1, 
         assimilation_period_seconds  = 0, 
         model_perturbation_amplitude = 0.2, 
         output_state_vector          = .false.  /

| 

.. container:: namelist

   ::

      namelist /CAL_NML/  TheCalendar, startDate_1, startDate_2, calendarDumps

.. container:: indent1

   | This namelist is read in a file called ``data.cal`` This namelist is the same one that is used by the ocean model.
     The values **must** correspond to the date at the start of an experiment. This is more important for
     ``create_ocean_obs, trans_pv_sv`` than for ``filter`` and :doc:`./trans_sv_pv` since ``trans_sv_pv`` takes the
     start time of the experiment from the DART initial conditions file and actually writes a new ``data.cal.DART`` and
     a new ``data.DART`` file. ``advance_model.csh`` renames ``data.DART`` and ``data.cal.DART`` to be used for the
     model advance.
   | Still, the files must exist before DART runs to avoid unnecessarily complex logic. If you are running the support
     programs in a standalone fashion (as you might if you are converting snapshot files into an intial ensemble), it is
     critical that the values in this namelist are correct to have accurate times in the headers of the restart files.
     You can always patch the times in the headers with ``restart_file_utility``.

| 

.. container:: namelist

   ::

      namelist /PARM03/  startTime, endTime, deltaTmom, &
                              deltaTtracer, deltaTClock, dumpFreq, taveFreq, ...

.. container:: indent1

   | This namelist is read in a file called ``data``. This namelist is the same one that is used by the ocean model.
     Only the variables listed here are used by the DART programs, there are more variables that are used only by the
     ocean model.
   | There are two scenarios of interest for this namelist.

   #. During an experiment, the ``advance_model.csh`` script is invoked by ``filter`` and the namelist is read by
      ``trans_sv_pv`` and REWRITTEN for use by the ocean model. Since this all happens in a local directory for the
      model advance, only a copy of the input ``data`` file is overwritten. The intent is that the ``data`` file is
      preserved 'perfectly' except for the values in ``&PARM03`` that pertain to controlling the model advance:
      ``endTime``, ``dumpFreq``, and ``taveFreq``.
   #. Outside the confines of ``trans_sv_pv``, this namelist is always simply read and is unchanged.

   +--------------------------------------+----------+----------------------------------------------------+
   | Contents                             | Type     | Description                                        |
   +======================================+==========+====================================================+
   | startTime                            | real(r8) | This **must** be 0.0 to tell the ocean model to    |
   |                                      |          | read from the input files named in                 |
   |                                      |          | ``data``\ ``&PARM05``.                             |
   +--------------------------------------+----------+----------------------------------------------------+
   | endTime                              | real(r8) | The number of seconds for one model advance.       |
   |                                      |          | (normally set by ``trans_sv_pv``)                  |
   +--------------------------------------+----------+----------------------------------------------------+
   | deltaTmom, deltaTtracer, deltaTClock | real(r8) | These are used when trying to interpret the        |
   |                                      |          | timestepcount in the snapshot files. They must all |
   |                                      |          | be identical unless someone can tell me which one  |
   |                                      |          | is used when the ocean model creates snapshot      |
   |                                      |          | filenames.                                         |
   +--------------------------------------+----------+----------------------------------------------------+
   | dumpFreq, taveFreq                   | real(r8) | Set to the same value value as ``endTime``. I have |
   |                                      |          | never run with different settings, my one concern  |
   |                                      |          | would be how this affects a crappy piece of logic  |
   |                                      |          | in ``advance_model.csh`` that requires there to be |
   |                                      |          | exactly ONE set of snapshot files - and that they  |
   |                                      |          | correspond to the completed model advance.         |
   +--------------------------------------+----------+----------------------------------------------------+

   This namelist is the same one that is used by the ocean model. Only some of the namelist variables are needed by
   DART; the rest are ignored by DART but could be needed by the ocean model. Here is a fragment for a daily
   assimilation timestep with the model dynamics having a much shorter timestep.

   .. rubric:: Parm03 namelist
      :name: parm03-namelist
      :class: indent1

   ::

      &PARM03
         startTime    =     0.,
           endTime    = 86400.,
         deltaTmom    =   900.,
         deltaTtracer =   900.,
         deltaTClock  =   900.,
         dumpFreq     = 86400.,
         taveFreq     = 86400.,
           ...

   This would result in snapshot files with names like ``[S,T,U,V,Eta].0000000096.data`` since 86400/900 = 96. These
   values remain fixed for the entire assimilation experiment, the only thing that changes from the ocean model's
   perspective is a new ``data.cal`` gets created for every new assimilation cycle. ``filter`` is responsible for
   starting and stopping the ocean model. The DART model state has a valid time associated with it, this information is
   used to create the new ``data.cal``.

| 

.. container:: namelist

   ::

      namelist /PARM04/  phiMin, thetaMin, delY, delX, delZ, ...

.. container:: indent1

   This namelist is read in a file called ``data``. This namelist is the same one that is used by the ocean model. Only
   the variables listed here are used by the DART programs, there are more variables that are used only by the ocean
   model.

   +----------+---------------------------+-----------------------------------------------------------------------------+
   | Contents | Type                      | Description                                                                 |
   +==========+===========================+=============================================================================+
   | phiMin   | real(r8)                  | The latitude of the southmost grid edge. In degrees.                        |
   +----------+---------------------------+-----------------------------------------------------------------------------+
   | thetaMin | real(r8)                  | The longitude of the leftmost grid edge. In degrees.                        |
   +----------+---------------------------+-----------------------------------------------------------------------------+
   | delY     | real(r8), dimension(1024) | The latitudinal distance between grid cell edges. In degrees. The array has |
   |          |                           | a default value of 0.0. The number of non-zero entries determines the       |
   |          |                           | number of latitudes. static_init_model() converts the namelist values to    |
   |          |                           | grid centroids and edges.                                                   |
   +----------+---------------------------+-----------------------------------------------------------------------------+
   | delX     | real(r8), dimension(1024) | The longitudinal distance between grid cell edges. In degrees. The array    |
   |          |                           | has a default value of 0.0. The number of non-zero entries determines the   |
   |          |                           | number of longitudes. static_init_model() converts the namelist values to   |
   |          |                           | grid centroids and edges.                                                   |
   +----------+---------------------------+-----------------------------------------------------------------------------+
   | delZ     | real(r8), dimension(512)  | The vertical distance between grid cell edges i.e., the thickness of the    |
   |          |                           | layer. In meters. The array has a default value of 0.0. The number of       |
   |          |                           | non-zero entries determines the number of depths. static_init_model()       |
   |          |                           | converts the namelist values to grid centroids and edges.                   |
   +----------+---------------------------+-----------------------------------------------------------------------------+

   This namelist is the same one that is used by the ocean model. Only some of the namelist variables are needed by
   DART; the rest are ignored by DART but could be needed by the ocean model. Here is a fragment for a (NY=225, NX=256,
   NZ=...) grid

   .. rubric:: Parm04 namelist
      :name: parm04-namelist

   ::

      &PARM04
         phiMin   =     8.4,
         thetaMin =   262.0,
         delY     = 225*0.1,
         delX     = 256*0.1,
         delZ     =  5.0037,
                     5.5860,
                     6.2725,
                     7.0817,
                     8.0350,
                     9.1575,
                    10.4786,
                    12.0322,
                    13.8579,
                    16.0012,
                      ...

   Note that the ``225*0.1`` construct exploits the Fortran repeat mechanism to achieve 225 evenly-spaced gridpoints
   without having to manually enter 225 identical values. No such construct exists for the unevenly-spaced vertical
   layer thicknesses, so each layer thickness is explicitly entered.

| 

.. container:: namelist

   ::

      namelist /PARM05/  bathyFile, hydrogSaltFile, hydrogThetaFile, &
                       uVelInitFile, vVelInitFile, pSurfInitFile

.. container:: indent1

   This namelist is read in a file called ``data``. The only DART component to use this namelist is the shell script
   responsible for advancing the model - ``advance_model.csh``.

   +-----------------+------------------+-------------------------------------------------------------------------------+
   | Contents        | Type             | Description                                                                   |
   +=================+==================+===============================================================================+
   | bathyFile       | character(len=*) | The Fortran direct-access big-endian binary file containing the bathymetry.   |
   +-----------------+------------------+-------------------------------------------------------------------------------+
   | hydrogSaltFile  | character(len=*) | The Fortran direct-access big-endian binary (snapshot) file containing the    |
   |                 |                  | salinity. ``S.0000000096.data``, for example. Units: psu                      |
   +-----------------+------------------+-------------------------------------------------------------------------------+
   | hydrogThetaFile | character(len=*) | The Fortran direct-access big-endian binary (snapshot) file containing the    |
   |                 |                  | temperatures. ``T.0000000096.data``, for example. Units: degrees C            |
   +-----------------+------------------+-------------------------------------------------------------------------------+
   | uVelInitFile    | character(len=*) | The Fortran direct-access big-endian binary (snapshot) file containing the U  |
   |                 |                  | current velocities. ``U.0000000096.data``, for example. Units: m/s            |
   +-----------------+------------------+-------------------------------------------------------------------------------+
   | vVelInitFile    | character(len=*) | The Fortran direct-access big-endian binary (snapshot) file containing the V  |
   |                 |                  | current velocities. ``V.0000000096.data``, for example. Units: m/s            |
   +-----------------+------------------+-------------------------------------------------------------------------------+
   | pSurfInitFile   | character(len=*) | The Fortran direct-access big-endian binary (snapshot) file containing the    |
   |                 |                  | sea surface heights. ``Eta.0000000096.data``, for example. Units: m           |
   +-----------------+------------------+-------------------------------------------------------------------------------+

   This namelist specifies the input files to the ocean model. DART must create these input files. ``advance_model.csh``
   has an ugly block of code that actually 'reads' this namelist and extracts the names of the input files expected by
   the ocean model. ``advance_model.csh`` then **renames** the snapshot files to be that expected by the ocean model.
   For this reason (and several others) a DART experiment occurrs in a separate directory we call CENTRALDIR, and each
   model advance happens in a run-time subdirectory. The data files copied to the run-time directory are deemed to be
   volatile, i.e., we can overwrite them and change them during the course of an experiment.

| 

Files
-----

-  input namelist files: ``data, data.cal, input.nml``
-  input data file: ``filter_ics, perfect_ics``
-  output data files: ``[S,T,U,V,Eta].YYYYMMDD.HHMMSS.[data,meta]``

Please note that there are **many** more files needed to advance the ocean model, none of which are discussed here.

References
----------

-  none

Private components
------------------

N/A
