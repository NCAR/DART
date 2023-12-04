WRF-Hydro
=========

Overview
--------

The Weather Research and Forecasting Hydrologic Model 
(`WRF-Hydro <http://www.ral.ucar.edu/projects/wrf_hydro/overview>`_)
is a community modeling system and framework for hydrologic modeling and model
coupling. WRF-Hydro is configured to use the Noah-MP Land Surface Model to 
simulate land surface processes. Combined with DART, the facility is called
*HydroDART*.

The development of HydroDART was a collaboration between **James McCreight**
of the Research Applications Laboratory of NSF NCAR and **Moha Gharamti** of
the Data Assimilation Research Section of NSF NCAR.

Streamflow assimilation is an active area of research and provides many
interesting research challenges. 

Description of this directory within the DART repository
--------------------------------------------------------

Contents of the ``$DARTROOT/models/wrf_hydro/``:

.. code-block::

   ├── ensemble_config_files/
   │      # Files which configure ensembles in wrfhydropy.
   ├── experiment_config_files/
   │      # File which configure hydro_dart_py experiments.
   ├── hydro_dart_py/
   │      # Python package/library for configuring and executing experiments.
   ├── python/
   │      # Python scripts for various purposes.
   ├── R/
   │      # R scripts for various purposes.
   ├── shell_scripts/
   │      # Shell scripts for various purposes.
   ├── templates/
   │      # Obsolete?
   ├── work/
   │      # Dart executables build directory and other testing.
   ├── model_mod.html
   │      # The model_mod documentation.
   ├── model_mod.nml
   │      # The model_mod namelist (subsumed by work/input.nml)
   ├── model_mod.f90
   │      # The model_mod code.
   ├── noah_hydro_mod.f90
   │      # Some model_mod interfaces more specific to Noah?
   ├── create_identity_streamflow_obs.f90
   │      # For creating identity streamflow obs for the NHDPlus-based
   │      # channel-network configuration of WRF-Hydro.
   ├── README.rst
          # This file.

To set up an experiment
-----------------------

To set up an experiment, consult the ``./python/experiment`` directory.

Description of external directories on GLADE
--------------------------------------------

The gridded version of the model has bits/bobs in these directories:

- ``/gpfs/fs1/work/jamesmcc/domains/public/croton_NY/Gridded/DOMAIN``
- ``/gpfs/fs1/work/jamesmcc/domains/public/croton_NY/Gridded/RESTART``

Only the gridcells with flow are retained in the ``qlink[1,2]``, ``hlink``
variables, so they must be unpacked in EXACTLY the same way as wrfHydo packs
them from the grid to their 'sparse' representation.

Namelist
--------

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand ``&`` and terminate with a slash ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist.

.. code-block:: fortran

   &model_nml
        assimilation_period_days       = 0
        assimilation_period_seconds    = 3600
        lsm_model_choice               = 'noahMP'
        model_perturbation_amplitude   = 0.5
        perturb_distribution           = 'lognormal'
        max_link_distance              = 2000.0
        streamflow_4_local_multipliers = 0.0001
        debug                          = 0
        domain_order                   = 'hydro'
        domain_shapefiles              = 'restart.hydro.nc'
        lsm_variables    = 'SH2O',              'QTY_SOIL_LIQUID_WATER', '0.0',   'NA', 'NOUPDATE',
                           'SUBSURFACE_FLUX',   'QTY_SUBSURFACE',        '0.0',   'NA', 'NOUPDATE',
                           'OVERLAND_FLUX',     'QTY_OVERLAND_FLOW',     '0.0',   'NA', 'NOUPDATE'
        hydro_variables  = 'qlink1',            'QTY_STREAM_FLOW',       '0.0',   'NA', 'UPDATE',
                           'z_gwsubbas',        'QTY_AQUIFER_WATER',     'NA',    'NA', 'UPDATE'
        parameters       = 'qBucketMult',       'QTY_BUCKET_MULTIPLIER', '0.001', '50', 'UPDATE',
                           'qSfcLatRunoffMult', 'QTY_RUNOFF_MULTIPLIER', '0.001', '50', 'UPDATE'
   /

This namelist is read from a file called ``input.nml``. This namelist provides
control over the assimilation period for the model. All observations within
(+/-) half of the assimilation period are assimilated. The assimilation period
is the minimum amount of time the model can be advanced, and checks are
performed to ensure that the assimilation window is a multiple of the NOAH
model dynamical timestep.

+-------------------------------------+-----------------------------------+------------------------------------------+
| Item                                | Type                              | Description                              |
+=====================================+===================================+==========================================+
| assimilation_period_days            | integer                           | The number of days to advance the model  |
|                                     |                                   | for each assimilation. [default: ``1``]  |
+-------------------------------------+-----------------------------------+------------------------------------------+
| assimilation_period_seconds         | integer                           | In addition to                           |
|                                     |                                   | ``assimilation_period_days``, the number |
|                                     |                                   | of seconds to advance the model for each |
|                                     |                                   | assimilation. [default: ``0``]           |
+-------------------------------------+-----------------------------------+------------------------------------------+
| lsm_model_choice                    | character(len=128)                | case-insensitive specification of the    |
|                                     |                                   | Land Surface model. Valid values are     |
|                                     |                                   | ``noahmp`` and ``noahmp_36``             |
+-------------------------------------+-----------------------------------+------------------------------------------+
| model_perturbation_amplitude        | real(r8)                          | The amount of noise to add when trying   |
|                                     |                                   | to perturb a single state vector to      |
|                                     |                                   | create an ensemble. Only used when       |
|                                     |                                   | ``input.nml`` is set with                |
|                                     |                                   | ``&filter_nml:start_from_restart =       |
|                                     |                                   | .false.``. See also                      |
|                                     |                                   | Generating the initial ensemble.         |
|                                     |                                   | units: standard deviation of the         |
|                                     |                                   | specified distribution the mean at the   |
|                                     |                                   | value of the state vector element.       |
+-------------------------------------+-----------------------------------+------------------------------------------+
| perturb_distribution                | character(len=256)                | The switch to determine the distribution |
|                                     |                                   | of the perturbations used to create an   |
|                                     |                                   | initial ensemble from a single model     |
|                                     |                                   | state. Valid values are :                |
|                                     |                                   | ``lognormal`` or ``gaussian``            |
+-------------------------------------+-----------------------------------+------------------------------------------+
| max_link_distance                   | real(r8)                          | The along-the-stream localization        |
|                                     |                                   | distance. In meters.                     |
+-------------------------------------+-----------------------------------+------------------------------------------+
| streamflow_4_local_multipliers      | real(r8)                          |                                          |
+-------------------------------------+-----------------------------------+------------------------------------------+
| debug                               | integer                           | The switch to specify the run-time       |
|                                     |                                   | verbosity.                               |
|                                     |                                   |                                          |
|                                     |                                   | - ``0`` is as quiet as it gets           |
|                                     |                                   | - ``> 1`` provides more run-time         |
|                                     |                                   |   messages                               |
|                                     |                                   | - ``> 5`` provides ALL run-time          |
|                                     |                                   |   messages                               |
|                                     |                                   |                                          |
|                                     |                                   | All values above 0 will also write a     |
|                                     |                                   | netCDF file of the grid information and  |
|                                     |                                   | perform a grid interpolation test.       |
|                                     |                                   | [default: ``0``]                         |
+-------------------------------------+-----------------------------------+------------------------------------------+
| domain_order                        | character(len=256)::              | There are three possible domains to      |
|                                     | dimension(3)                      | include in the HydroDART state:          |
|                                     |                                   | ``hydro``, ``parameters``, ``lsm``       |
|                                     |                                   | This variable specifies the ordering of  |
|                                     |                                   | the domains.                             |
+-------------------------------------+-----------------------------------+------------------------------------------+
| domain_shapefiles                   | character(len=256)::              | There are input files used to determine  |
|                                     | dimension(3)                      | the shape of the input variables and any |
|                                     |                                   | geographic metadata.                     |
|                                     |                                   | They must be specified in the same       |
|                                     |                                   | order as listed in  ``domain_order``     |
+-------------------------------------+-----------------------------------+------------------------------------------+
| lsm_variables                       | character(len=32)::               | The list of variable names in the NOAH   |
|                                     | dimension(5,40)                   | restart file to use to create the DART   |
|                                     |                                   | state vector and their corresponding     |
|                                     |                                   | DART QUANTITY. [see example below]       |
+-------------------------------------+-----------------------------------+------------------------------------------+
| hydro_variables                     | character(len=32)::               | The list of variable names in the channel|
|                                     | dimension(5,40)                   | model file to use to create the DART     |
|                                     |                                   | state vector and their corresponding     |
|                                     |                                   | DART QUANTITY. [see example below]       |
+-------------------------------------+-----------------------------------+------------------------------------------+
| parameters                          | character(len=32)::               | The list of variable names in the        |
|                                     | dimension(5,40)                   | parameter file to use to create the DART |
|                                     |                                   | state vector and their corresponding     |
|                                     |                                   | DART QUANTITY. [see example below]       |
+-------------------------------------+-----------------------------------+------------------------------------------+


The columns of ``lsm_variables``, ``hydro_variables``, and ``parameters`` needs 
some explanation. Starting with the column 5,
``UPDATE`` denotes whether or not to replace the variable with the Posterior (i.e.
assimilated) value. Columns 3 and 4 denote lower and upper bounds that should be
enforced when writing to the files used to restart the model. These limits are not
enforced for the DART diagnostic files. Column 2 specifies the relationship between
the netCDF variable name for the model and the corresponding DART QUANTITY.

Support for these QUANTITYs is provided by
running ``preprocess`` with the following namelist settings:

.. code-block::

   &preprocess_nml
              overwrite_output = .true.
       input_obs_kind_mod_file = '../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90'
      output_obs_kind_mod_file = '../../../assimilation_code/modules/observations/obs_kind_mod.f90'
        input_obs_def_mod_file = '../../../observations/forward_operators/DEFAULT_obs_def_mod.F90'
       output_obs_def_mod_file = '../../../observations/forward_operators/obs_def_mod.f90'
      input_files              = '../../../observations/forward_operators/obs_def_streamflow_mod.f90',
                                 '../../../observations/forward_operators/obs_def_land_mod.f90',
                                 '../../../observations/forward_operators/obs_def_COSMOS_mod.f90'
     /
