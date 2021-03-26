NOAH, NOAH-MP
=============

Overview
--------

The Manhattan-compliant version of the NOAH (technically NOAH-MP) supports 
NOAH-MP V3.6 and was largely updated in support of the data assimilation
efforts with **wrf_hydro**. Experiments to perform data assimilation strictly
with the NOAH-MP model have been run at the University of Texas at Austin by
**Jingjing Liang**. We know other people are using DART and NOAH-MP.
*however, we have not had the chance to update the documentation for the 
Manhattan release.* Consequently, we readily welcome any advice on how to
improve the documentation and heartily encourage participation.


The `NOAH <http://www.ral.ucar.edu/research/land/technology/lsm.php>`_ **Land
Surface Model** and **Data Assimilation Research Testbed (DART)** may now be
used for assimilation experiments.
The Classic or Lanai version should be considered an 'alpha' release -- 
the code has only been tested for a single column configuration of NOAH. 

Any of the variables in the NOAH restart file are available to be adjusted by
the assimilation. The list of variables is set though a simple namelist
interface. Since we are testing in a column configuration, there is no
practical reason not to include all the variables necessary for a bit-for-bit
restart: *SOIL_T*, *SOIL_M*, *SOIL_W*, *SKINTEMP*, *SNODEP*, *WEASD*,
*CANWAT*, and *QFX*. These variables are then adjusted to be consistent with
real observations and stuffed back into the same netCDF restart files. Since
DART is an ensemble algorithm there are multiple restart files for a single
restart time; one for each ensemble member. Creating the initial ensemble of
land surface states is an area of active research. At present, it may be
sufficient to use a climatological ensemble; e.g., using the restarts for '1
January 00Z' from 50 consecutive years.

There is reason to believe that the ensemble system will benefit from 
having unique atmospheric forcing for each ensemble member.
A reasonable ensemble size is 50 or 80 or so.

DART reads the NOAH namelist ``&NOAHLSM_OFFLINE`` from a file called
``namelist.hrldas`` for several pieces of information. DART is responsible for
starting/stopping NOAH; the restart information is conveyed through the NOAH
namelist. **Unpleasant Reality #1 :** managing the tremendous number of hourly
forcing files for every ensemble member is tedious. To facilitate matters, the
DART/NOAH system uses a *single* netCDF file for each ensemble member that
contains ALL of the forcing for that ensemble member.

+------------------------------------------+--------------------------------------------------------------------------+
| `dart_to_noah.f90 <dart_to_noah.html>`__ | **updates** some or all of a NOAH restart file with the posterior DART   |
|                                          | state vector. There is the ability to selectively avoid updating the     |
|                                          | NOAH variables. This allows one to include NOAH variables in the DART    |
|                                          | state vector to aid in the application of observation operators, etc.,   |
|                                          | without having to modify those variables in the NOAH restart file.       |
|                                          | `[dart_to_noah.html] <dart_to_noah.html>`__                              |
+------------------------------------------+--------------------------------------------------------------------------+

Running a "Perfect Model" experiment ... OSSE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The example requires a basic knowledge of running NOAH.
Four scripts are provided to demonstrate how to set up and run a perfect model
experiment for a single site - with one caveat. You must provide your own
initial ensemble for the experiment. The scripts are not intended to be black
boxes. You are expected to read them and modify them to your own purpose.

The scripts assume the directory containing the DART executables is
``${DARTDIR}/work``, and assume that the directory containing the NOAH
executables is ``${NOAHDIR}/Run``.

+----------------------------------------------------------+------------------------------------------------------------------------------+
| 1. ``shell_scripts/setup_pmo.csh``                       | This script stages the run of                                                |
|                                                          | :doc:`../../assimilation_code/programs/perfect_model_obs/perfect_model_obs`. |
|                                                          | The directory where you run the script is called ``CENTRALDIR`` and          |
|                                                          | will be the working directory for the experiment. The required input         |
|                                                          | observation sequence file must be created in the normal DART way. This       |
|                                                          | ``obs_seq.in`` file must exist before running this script. All the           |
|                                                          | necessary data files and exectuables for a perfect model experiment get      |
|                                                          | copied to CENTRALDIR so that you may run multiple experiments at the         |
|                                                          | same time - in separate ``CENTRALDIRs``.                                     |
+----------------------------------------------------------+------------------------------------------------------------------------------+
| 2. ``shell_scripts/run_pmo.csh``                         | very simply - it advances NOAH and applies the observation operator to       |
|                                                          | put the "perfect" observations in an observation sequence file that can      |
|                                                          | then be used for an assimilation.                                            |
+----------------------------------------------------------+------------------------------------------------------------------------------+
| 3. ``shell_scripts/setup_filter.csh``                    | builds upon the work of ``setup_pmo.csh`` and stages a PRE-EXISTING          |
|                                                          | initial ensemble.                                                            |
+----------------------------------------------------------+------------------------------------------------------------------------------+
| 4. ``shell_scripts/run_filter.csh``                      | Actually runs the filtering (assimilation) experiment.                       |
+----------------------------------------------------------+------------------------------------------------------------------------------+


Generating the initial ensemble
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Creating the initial ensemble of soil moisture states is an area of active
research. The ensemble must come from 'somewhere else'. At present, it may be
sufficient to use a climatological ensemble; e.g., using the NOAH restarts for
'1 January 00Z' from 50 consecutive years from a hindcast experiment. It may
also be sufficient to take a single model state, replicate it N times and
force each of the N instances with different atmospheric conditions for 'a
long time'.

By The Way
~~~~~~~~~~

Experience has shown that having a paired (unique) atmospheric forcing maintains
the ensemble spread during an assimilation better than simply forcing all the
ensemble members with one single atmospheric state.

DART has routines to perturb a single NOAH state and generate its own ensemble
(typically done with ``pert_model_state``), but this produces model states that
are incompatible with NOAH. We are interested in adopting/adapting strategies
to create sensible initial conditions for NOAH.

If you have an algorithm you believe will be useful, please contact us!

Observations
------------

Some novel observations come from the Cosmic-ray Soil Moisture Observing System:
`COSMOS <http://cosmos.hwr.arizona.edu/>`__ and are processed by DART routines
in the ``$DARTROOT/observations/COSMOS`` directory.

DART has a very object-oriented approach to observation support. All
observations that are intended to be supported must be preprocessed (see 
``$DARTROOT/preprocess/`` into a single ``obs_def_mod.f90`` and
``obs_kind_mod.f90`` in the standard DART way.

Exploring the Output
~~~~~~~~~~~~~~~~~~~~

There are Matlab® scripts for exploring the performance of the assimilation in
observation-space (after running ``obs_diag``). See ``$DARTROOT/diagnostics/threed_sphere/obs_diag.html``
to explore the *obs_seq.final* file) - use the scripts starting with ``plot_``,
i.e. ``$DARTROOT/diagnostics/matlab/plot_*.m*``. As always, there are some
model-specific items Matlab® will need to know about in
``$DARTROOT/models/NOAH/matlab``.

The ``Prior_Diag.nc`` and ``Posterior_Diag.nc`` (and possibly ``True_State.nc``)
netCDF files have the model prognostic variables before and after the
assimilation. The ``./matlab`` scripts for NOAH are under development.

It is also worthwhile to convert your ``obs_seq.final`` file to a netCDF format
obs_sequence file with ``obs_seq_to_netcdf``. See
``$DARTROOT/obs_sequence/obs_seq_to_netcdf.html`` and use any of the standard
plots. Be aware that the COSMOS site-specific metadata will not get conveyed to
the netCDF file.

Namelist
--------

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand ``&`` and terminate with a slash ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist. The standard values are shown below:

.. code-block:: fortran

   &model_nml
      lsm_model_choice             = 'noahMP_36'
      domain_shapefiles            = 'RESTART.2003051600_DOMAIN1_01'
      assimilation_period_days     =    0
      assimilation_period_seconds  = 3600
      model_perturbation_amplitude = 0.2
      perturb_distribution         = 'gaussian'
      debug                        = 0
      polar                        = .false.
      periodic_x                   = .false.
      periodic_y                   = .false.
      lsm_variables = 'SOIL_T',   'QTY_SOIL_TEMPERATURE',   '0.0',  'NA', 'UPDATE',
                      'SMC',      'QTY_SOIL_MOISTURE',      '0.0', '1.0', 'UPDATE',
                      'WA',       'QTY_AQUIFER_WATER',      '0.0',  'NA', 'UPDATE',
                      'SNEQV',    'QTY_SNOW_WATER',         '0.0',  'NA', 'UPDATE',
                      'FSNO',     'QTY_SNOWCOVER_FRAC',     '0.0', '1.0', 'UPDATE'
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
| lsm_model_choice                    | character(len=256)                | The version of the NOAH namelist to read |
+-------------------------------------+-----------------------------------+------------------------------------------+
| domain_shapefiles                   | an array of character(len=256)    | The name of the NOAH RESTART files to    |
|                                     |                                   | use to specify the shape of the variables|
|                                     |                                   | and geographic metadata. One per domain. |
+-------------------------------------+-----------------------------------+------------------------------------------+
| assimilation_period_days            | integer                           | The number of days to advance the model  |
|                                     |                                   | for each assimilation.                   |
+-------------------------------------+-----------------------------------+------------------------------------------+
| assimilation_period_seconds         | integer                           | In addition to                           |
|                                     |                                   | ``assimilation_period_days``, the number |
|                                     |                                   | of seconds to advance the model for each |
|                                     |                                   | assimilation.                            |
+-------------------------------------+-----------------------------------+------------------------------------------+
| model_perturbation_amplitude        | real(r8)                          | The amount of noise to add when trying   |
|                                     |                                   | to perturb a single state vector to      |
|                                     |                                   | create an ensemble. Only used when       |
|                                     |                                   | ``input.nml`` is set with                |
|                                     |                                   | ``&filter_nml:start_from_restart =       |
|                                     |                                   | .false.``. See also                      |
|                                     |                                   | `Generating the initial ensemble`_.      |
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
| periodic_x                          | logical                           | Switch to determine if the configuration |
|                                     |                                   | has periodicity in the X direction.      |
+-------------------------------------+-----------------------------------+------------------------------------------+
| periodic_y                          | logical                           | Switch to determine if the configuration |
|                                     |                                   | has periodicity in the Y direction.      |
+-------------------------------------+-----------------------------------+------------------------------------------+
| lsm_variables                       | character(len=32)::               | The list of variable names in the NOAH   |
|                                     | dimension(5,40)                   | restart file to use to create the DART   |
|                                     |                                   | state vector and their corresponding     |
|                                     |                                   | DART kind. [default: see example below]  |
+-------------------------------------+-----------------------------------+------------------------------------------+


The columns of ``lsm_variables`` needs some explanation. Starting with the column 5, 
``UPDATE`` denotes whether or not to replace the variable with the Posterior (i.e. 
assimilated) value. Columns 3 and 4 denote lower and upper bounds that should be 
enforced when writing to the files used to restart the model. These limits are not 
enforced for the DART diagnostic files. Column 2 specifies the relationship between 
the netCDF variable name for the model and the corresponding DART QUANTITY.

The DART 'QTY's match what the ``model_mod`` knows how to interpolate, so you can't
just add a new quantity and expect it to work. There is a complex interplay between
``obs_def_mod`` and ``preprocess``, and ``model_mod`` that defines what QUANTITIES
are supported. There is only a single QUANTITY that works with each variable and
the example shows the current QUANTITYs. Support for these QUANTITYs was provided by
running ``preprocess`` with the following namelist settings:

.. code-block::

   &preprocess_nml
       input_obs_kind_mod_file = '../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90'
      output_obs_kind_mod_file = '../../../assimilation_code/modules/observations/obs_kind_mod.f90'
        input_obs_def_mod_file = '../../../observations/forward_operators/DEFAULT_obs_def_mod.F90'
       output_obs_def_mod_file = '../../../observations/forward_operators/obs_def_mod.f90'
      input_files              = '../../../observations/forward_operators/obs_def_land_mod.f90',
                                 '../../../observations/forward_operators/obs_def_COSMOS_mod.f90',
                                 '../../../observations/forward_operators/obs_def_GRACE_mod.f90'
     /


NOAHLSM_OFFLINE NAMELIST
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: fortran

   namelist /NOAHLSM_OFFLINE/
      hrldas_constants_file, &
      indir, outdir,  &
      restart_filename_requested, &
      khour,  kday, &
      forcing_timestep, &
      noah_timestep,  &
      output_timestep, &
      restart_frequency_hours, &
      split_output_count, &
      nsoil, &
      zsoil

The remaining variables are not used by DART - but are used by NOAH. Since DART
verifies namelist accuracy, any namelist entry in NOAHLSM_OFFLINE that is not
in the following list will cause a FATAL DART ERROR.

.. code-block:: fortran

   zlvl, zlvl_wind, iz0tlnd, sfcdif_option, update_snow_from_forcing,
   start_year, start_month, start_day, start_hour, start_min,
   external_fpar_filename_template, external_lai_filename_template,
   subwindow_xstart, subwindow_xend, subwindow_ystart, subwindow_yend

This namelist is read from a file called ``namelist.hrldas``. This namelist is
the same one that is used by NOAH. The values are explained in full in the NOAH
documentation. Only the namelist variables of interest to DART are discussed.
All other namelist variables are ignored by DART - but mean something to NOAH.

+-------------------------------------+-----------------------------------+------------------------------------------+
| Item                                | Type                              | Description                              |
+=====================================+===================================+==========================================+
| hrldas_constants_file               | character(len=256)                | The name of the netCDF file containing   |
|                                     |                                   | the grid information. [default:          |
|                                     |                                   | ``wrfinput``]                            |
+-------------------------------------+-----------------------------------+------------------------------------------+
| indir                               | character(len=256)                | The DART/NOAH environment requires all   |
|                                     |                                   | the input files to be in the current     |
|                                     |                                   | working directory. [default: ``'.'``]    |
+-------------------------------------+-----------------------------------+------------------------------------------+
| outdir                              | character(len=256)                | The DART/NOAH environment requires all   |
|                                     |                                   | output files are in the current working  |
|                                     |                                   | directory. [default: ``'.'``]            |
+-------------------------------------+-----------------------------------+------------------------------------------+
| restart_filename_requested          | character(len=256)                | The name of the file containing the grid |
|                                     |                                   | information. The default value is        |
|                                     |                                   | implicitly used by the scripting         | 
|                                     |                                   | examples. Change at your own risk.       |
|                                     |                                   | [default: ``'restart.nc'``]              |
+-------------------------------------+-----------------------------------+------------------------------------------+
| khour                               | integer                           | The duration (in hours) of the model     |
|                                     |                                   | integration. [default: ``1``]            |
+-------------------------------------+-----------------------------------+------------------------------------------+
| kday                                | integer                           | The duration (in days) of the model      |
|                                     |                                   | integration. [default: ``0``]            |
+-------------------------------------+-----------------------------------+------------------------------------------+
| forcing_timestep                    | integer                           | The timestep (in seconds) of the         |
|                                     |                                   | atmospheric forcing. [default: ``3600``] |
+-------------------------------------+-----------------------------------+------------------------------------------+
| noah_timestep                       | integer                           | The internal (dynamical) timestep (in    |
|                                     |                                   | seconds). [default: ``3600``]            |
+-------------------------------------+-----------------------------------+------------------------------------------+
| output_timestep                     | integer                           | The output interval (in seconds).        |
|                                     |                                   | [default: ``3600``]                      |
+-------------------------------------+-----------------------------------+------------------------------------------+
| restart_frequency_hours             | integer                           | How often the NOAH restart files get     |
|                                     |                                   | written. [default: ``1``]                |
+-------------------------------------+-----------------------------------+------------------------------------------+
| split_output_count                  | integer                           | should be 1 or bad things happen.        |
|                                     |                                   | [default: ``1``]                         |
+-------------------------------------+-----------------------------------+------------------------------------------+
| nsoil                               | integer                           | The number of soil interfaces. As I      |
|                                     |                                   | understand it, NOAH requires this to be  |
|                                     |                                   | 4. [default: ``4``]                      |
+-------------------------------------+-----------------------------------+------------------------------------------+
| zsoil                               | integer(NSOLDX)                   | The depth (in meters) of the soil        |
|                                     |                                   | interfaces. [default: ``-0.1, -0.4,      |
|                                     |                                   | -1.0, -2.04``]                           |
+-------------------------------------+-----------------------------------+------------------------------------------+

Example
~~~~~~~

Note: the ``FORCING_FILE_DIRECTORY`` line is not required by NOAH but IS required
by DART - specifically in the *advance_model.csh* script.

.. code-block:: fortran

   ### THIS IS FOR DART ###
   FORCING_FILE_DIRECTORY = "/path/to/your/forcing/files"
   
   &NOAHLSM_OFFLINE
      HRLDAS_CONSTANTS_FILE = "wrfinput"
      INDIR  = "."
      OUTDIR = "."
      RESTART_FILENAME_REQUESTED = "restart.nc"
      KHOUR                   = 1
      FORCING_TIMESTEP        = 3600
      NOAH_TIMESTEP           = 3600
      OUTPUT_TIMESTEP         = 3600
      RESTART_FREQUENCY_HOURS = 1
      SPLIT_OUTPUT_COUNT      = 1
      NSOIL=4
      ZSOIL(1) = -0.10
      ZSOIL(2) = -0.40
      ZSOIL(3) = -1.00
      ZSOIL(4) = -2.00
   /


Input Files
-----------

+-----------------------------------+-----------------------------------------+
| filename                          | purpose                                 |
+===================================+=========================================+
| input.nml                         | to read the model_mod namelist          |
+-----------------------------------+-----------------------------------------+
| namelist.hrldas                   | to read the NOAHLSM_OFLINE namelist     |
+-----------------------------------+-----------------------------------------+
| wrfinput                          | provides NOAH grid information          |
+-----------------------------------+-----------------------------------------+
| *&model_nml:noah_netcdf_filename* | the RESTART file containing the NOAH    |
|                                   | model state.                            |
+-----------------------------------+-----------------------------------------+
