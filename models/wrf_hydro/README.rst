################
WRF-Hydro README
################

Contents
========

#. `Overview`_
#. `Description of this Directory within the DART Repository`_
#. `To Setup an Experiment`_
#. `Description of External Directories on GLADE`_
#. `For the Impatient - Skip the Fundamentals and Get Frustrated`_
#. `Observations`_
#. `Namelist`_
#. `Input Files`_
#. `Output Files`_
#. `Terms of Use`_

Overview
========





The `NOAH <http://www.ral.ucar.edu/research/land/technology/lsm.php>`_ **Land
Surface Model** and **Data Assimilation Research Testbed (DART)** may now be
used for assimilation experiments.

This should be considered an 'alpha' release -- the code has only been tested
for a single column configuration of NOAH. The offline 2D driver code:
`High-Resolution Land Data Assimilation System (HRLDAS) v3.3
<http://www.ral.ucar.edu/research/land/technology/lsm.php>`_ (April 2011) is the
version used to develop the interface. This must be downloaded and compiled
outside the framework of DART.

The ``simple_driver-v3.3`` was deemed suboptimal because of poor restart
characteristics. Since the offline 2D driver code can be run at a single
gridpoint and has excellent restart characteristics, there are no plans to
support the ``simple-driver-v3.3`` version of NOAH.

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

The ``offline hrldas-v3.3`` driver is designed to use atmospheric forcing files
that cover a VERY short period - about every hour. Consequently, the
directories that contain the forcing files get very cluttered for experiments
that may cover years, for example. There is reason to believe that the
ensemble system will benefit from having unique atmospheric forcing for each
ensemble member. A reasonable ensemble size is 50 or 80 or so. You do the
math.

DART reads the NOAH namelist ``&NOAHLSM_OFFLINE`` from a file called
``namelist.hrldas`` for several pieces of information. DART is responsible for
starting/stopping NOAH; the restart information is conveyed through the NOAH
namelist. **Unpleasant Reality #1 :** managing the tremendous number of hourly
forcing files for every ensemble member is tedious. To facilitate matters, the
DART/NOAH system uses a *single* netCDF file for each ensemble member that
contains ALL of the forcing for that ensemble member.

The DART ``advance_model.csh`` in the ``./shell_scripts/`` directory for NOAH
uses the ``ncks`` program to extract the required forcing timesteps into the
unique forcing files required by NOAH. The simplest way we could think of to
provide the location of the 'master' forcing files was by adding a *comment* in
the ``namelist.hrldas`` file of the form:

.. code-block:: fortran

  FORCING_FILE_DIRECTORY = "/path/to/your/forcing/files"

We felt it important to specify this information in the ``namelist.hrldas`` file
since it is a strong analogue of the traditional NOAH ``INDIR`` variable.

Description of this Directory within the DART Repository
========================================================

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

To Setup an Experiment
======================

Please consult the ``./python/experiment`` directory.

Description of External Directories on GLADE
============================================

The gridded version of the model has bits/bobs in these directories:

- ``/gpfs/fs1/work/jamesmcc/domains/public/croton_NY/Gridded/DOMAIN``
- ``/gpfs/fs1/work/jamesmcc/domains/public/croton_NY/Gridded/RESTART``

Only the gridcells with flow are retained in the ``qlink[1,2]``, ``hlink``
variables, so they must be unpacked in EXACTLY the same way as wrfHydo packs
them from the grid to their 'sparse' representation.

For the Impatient - Skip the Fundamentals and Get Frustrated
============================================================

There are several scripts in the *DART/models/noah/shell_scripts* directory that
assist in configuring and running experiments: *setup_pmo.csh*, *run_pmo.csh*,
*setup_filter.csh*, and *run_filter.csh* . The scripts are not intended to be
black boxes. **You are expected to read them and modify them to your own
purpose.**

The offline 2D HRLDAS-V3.3
--------------------------

The offline 2D HRLDAS-V3.3 and the development branch of DART was tested and run
on a MacBook Pro laptop running OS X 10.7.4 (Lion) with gfortran v 4.5.0. The
*Noah_hrldas_beta* source code was downloaded from `the Unified NOAH/LSM
website <http://www.ral.ucar.edu/research/land/technology/lsm.php>`__ and built
- outside the DART framework - to produce a standalone executable.

Trivial modifications to the distributed hrldas makefiles were necessary to
compile *the Noah_hrldas_beta* on a case-insensitive filesystem. All of the
original 14 *Makefile* files have to be changed:

.. code-block:: perl

   $(CPP) $(CPPFLAGS) $(CPPHRLDAS) $(*).F > $(*).f
   to
   $(CPP) $(CPPFLAGS) $(CPPHRLDAS) $(*).F > $(*).\ *f90*

In essense, all instances of ".f" must be changed to ".f90". If you do not have
a case-insensitive filesystem, the original Makefiles will not need these
modifications.

The DART components were built for debugging with the following settings within
``mkmf.template``:

.. code-block:: perl

   FC = gfortran
   LD = gfortran
   INCS = -I${NETCDF}/include
   LIBS = -L${NETCDF}/lib -lnetcdff -lnetcdf -lcurl -lhdf5_hl -lhdf5 -lz  -lm
   FFLAGS = -g -O0 -ffree-line-length-none -fbounds-check -frecord-marker=4 -ffpe-trap=invalid $(INCS)
   LDFLAGS = $(FFLAGS) $(LIBS)
      
A note about NOAH times and dates
---------------------------------

HRLDAS-V3.3 has some unusual conventions as pertains the contents of 'restart'
files. You should confirm the following to prove it to yourself - **outside** of
the DART framework. Your actual numbers will change, but the point should be
obvious.

#. Assume there is a NOAH restart file:

   - ``RESTART.2009010118_DOMAIN1``, containing:

     .. code-block::

        Times = "2009-01-01_18:00:00" ;
        SOIL_M = 0.100843, 0.1324335, 0.1104331, 0.16349 ;
        QFX = 12.03423 ;
             

#. Now, run/advance NOAH for a single (1 hour) timestep (i.e. from 18Z to 19Z).
   The following files are created:
   
   - ``2009010118.LDASOUT_DOMAIN1``, containing:

     .. code-block::
        
        Times = "2009-01-01_18:00:00" ;
        SOIL_M = 0.1007659, 0.132411, 0.1104217, 0.1634897 ;
        QFX = 14.90484 ;

   - ``RESTART.2009010119_DOMAIN1``, containing:

     .. code-block::

        Times = "2009-01-01_19:00:00" ;
        SOIL_M = 0.1007659, 0.132411, 0.1104217, 0.1634897 ;
        QFX = 14.90484 ;

   - ``2009010119.LDASOUT_DOMAIN1``, containing:

     .. code-block::
     
        Times = "2009-01-01_19:00:00" ;
        SOIL_M = 0.1006793, 0.1323851, 0.1104083, 0.1634894 ;
        QFX = 17.13207 ;
             

The contents of the ``2009010118.LDASOUT_DOMAIN1`` contain the **same** values
as the ``RESTART.2009010119_DOMAIN1`` -- **although the names and Times imply
they are from different times.** The *Times* in the LDASOUT files are
fundamentally "now"casts and reflect the valid time of the model state. This is
the time that DART requires.

.. important::

   DART reads and modifies the RESTART files. DART internally adjusts the times
   in the restart files to correspond to the time of the companion LDASOUT file.
   The *Times* array is unchanged. NOTE: The DART/NOAH interface uses the
   Gregorian calendar.

Converting between DART files and NOAH restart files
----------------------------------------------------

The information about how the NOAH variables are stored in the DART state vector
comes from the order in which the variables are specified in the ``input.nml``
file's ``&model_nml:noah_state_variables`` entry, as is the the name of the NOAH
restart file.

DART also needs to read the ``namelist.hrldas`` ``&NOAHLSM_OFFLINE`` namelist.

There are two programs - both use the ``model_mod`` module, and both have their
own documentation:


+------------------------------------------+--------------------------------------------------------------------------+
| `noah_to_dart.f90 <noah_to_dart.html>`__ | converts a NOAH restart file into a DART-compatible file normally called |
|                                          | *dart_ics* . We usually wind up linking the NOAH restart files to a      |
|                                          | static name (*restart.nc*). `[noah_to_dart.html] <noah_to_dart.html>`__  |
+------------------------------------------+--------------------------------------------------------------------------+
| `dart_to_noah.f90 <dart_to_noah.html>`__ | **updates** some or all of a NOAH restart file with the posterior DART   |
|                                          | state vector. There is the ability to selectively avoid updating the     |
|                                          | NOAH variables. This allows one to include NOAH variables in the DART    |
|                                          | state vector to aid in the application of observation operators, etc.,   |
|                                          | without having to modify those variables in the NOAH restart file.       |
|                                          | `[dart_to_noah.html] <dart_to_noah.html>`__                              |
+------------------------------------------+--------------------------------------------------------------------------+

Running a "Perfect Model" experiment ... OSSE
---------------------------------------------

The example requires a basic knowledge of running NOAH for a single column.
Since the single-column version requires such small netCDF input files, it is
possible to simply use *ncdump* to convert these files to text files, edit them,
and generate new netCDF files using *ncgen*. An appropriate ``wrfinput`` file
may be generated by editing the ``templates/wrfinput.template.ascii`` file and
using ``ncgen``, for example.

Four scripts are provided to demonstrate how to set up and run a perfect model
experiment for a single site - with one caveat. You must provide your own
initial ensemble for the experiment. The scripts are not intended to be black
boxes. You are expected to read them and modify them to your own purpose.

The scripts assume the directory containing the DART executables is
``${DARTDIR}/work``, and assume that the directory containing the NOAH
executables is ``${NOAHDIR}/Run``.

+----------------------------------------------------------+----------------------------------------------------------------+
| 1. ` setup_pmo.csh <shell_scripts/setup_pmo.csh>`_       | This script stages the run of                                  |
|                                                          | `perfect_model_obs                                             |
|                                                          |  <../../perfect_model_obs/perfect_model_obs.html>`__.          |
|                                                          | The directory where you run the script is called               |
|                                                          | ``CENTRALDIR`` and will be the working directory for the       |
|                                                          | experiment. The required input observation sequence file       |
|                                                          | must be created in the normal DART way (`one way is to         |
|                                                          | create synthetic observations                                  |
|                                                          | <https://dart.ucar.edu/pages/Observations.html#obs_seq_osse>`_ |
|                                                          | and must exist before running this script. All the             |
|                                                          | necessary data files and exectuables for a perfect model       |
|                                                          | experiment get copied to CENTRALDIR so that you may run        | 
|                                                          | multiple experiments at the same time - in separate            |
|                                                          | ``CENTRALDIRs``.                                               |
+----------------------------------------------------------+----------------------------------------------------------------+
| 2. ` run_pmo.csh <shell_scripts/run_pmo.csh>`_           | very simply - it advances NOAH and applies the observation     |
|                                                          | operator to put the "perfect" observations in an observation   |
|                                                          | sequence file that can then be used for an assimilation.       |
+----------------------------------------------------------+----------------------------------------------------------------+
| 3. `setup_filter.csh <shell_scripts/setup_filter.csh>`_  | builds upon the work of ``setup_pmo.csh`` and stages a         |
|                                                          | PRE-EXISTING initial ensemble.                                 |
+----------------------------------------------------------+----------------------------------------------------------------+
| 4. `run_filter.csh <shell_scripts/run_filter.csh>`_      | Actually runs the filtering (assimilation) experiment.         |
+----------------------------------------------------------+----------------------------------------------------------------+

Generating the initial ensemble
-------------------------------

Creating the initial ensemble of soil moisture states is an area of active
research. The ensemble must come from 'somewhere else'. At present, it may be
sufficient to use a climatological ensemble; e.g., using the NOAH restarts for
'1 January 00Z' from 50 consecutive years from a hindcast experiment. It may
also be sufficient to take a single model state, replicate it N times and
force each of the N instances with different atmospheric conditions for 'a
long time'.

By The Way
~~~~~~~~~
Experience has shown that having a paired (unique) atmospheric forcing maintains
the ensemble spread during an assimilation better than simply forcing all the
ensemble members with one single atmospheric state.

DART has routines to perturb a single NOAH state and generate its own ensemble
(typically done with ``pert_model_state``), but this produces model states that
are incompatible with NOAH. We are interested in adopting/adapting strategies
to create sensible initial conditions for NOAH.

If you have an algorithm you believe will be useful, please contact us!

Observations
============

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
========

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand ``&`` and terminate with a slash ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist.

.. code-block:: fortran

   &model_nml
      noah_netcdf_filename,
      assimilation_period_days,
      assimilation_period_seconds,
      model_perturbation_amplitude,
      output_state_vector,
      debug,
      noah_state_variables
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
| noah_netcdf_filename                | character(len=128)                | The name of the NOAH RESTART file to     |
|                                     |                                   | use to create the DART state vector.     |
|                                     |                                   | For convenience, the                     |
|                                     |                                   | ``advance_model.csh`` script usually     |
|                                     |                                   | links the most recent restart file to    |
|                                     |                                   | a static name. [default: ``restart.nc``] |
+-------------------------------------+-----------------------------------+------------------------------------------+
| assimilation_period_days            | integer                           | The number of days to advance the model  |
|                                     |                                   | for each assimilation. [default: ``1``]  |
+-------------------------------------+-----------------------------------+------------------------------------------+
| assimilation_period_seconds         | integer                           | In addition to                           |
|                                     |                                   | ``assimilation_period_days``, the number |
|                                     |                                   | of seconds to advance the model for each |
|                                     |                                   | assimilation. [default: ``0``]           |
+-------------------------------------+-----------------------------------+------------------------------------------+
| model_perturbation_amplitude        | real(r8)                          | The amount of noise to add when trying   |
|                                     |                                   | to perturb a single state vector to      |
|                                     |                                   | create an ensemble. Only used when       |
|                                     |                                   | ``input.nml`` is set with                |
|                                     |                                   | ``&filter_nml:start_from_restart =       |
|                                     |                                   | .false.``. See also                      |
|                                     |                                   | `Generating the initial ensemble`_.      |
|                                     |                                   | units: standard deviation of a Gaussian  |
|                                     |                                   | distribution with the mean at the value  |
|                                     |                                   | of the state vector element.             |
|                                     |                                   | [default: ``0.2``]                       |
+-------------------------------------+-----------------------------------+------------------------------------------+
| output_state_vector                 | logical                           | The switch to determine the form of the  |
|                                     |                                   | state vector in the output netCDF files. |
|                                     |                                   | If ``.true.`` the state vector will be   |
|                                     |                                   | output exactly as DART uses it, as one   |
|                                     |                                   | long array. If ``*.false.*``, the state  |
|                                     |                                   | vector is parsed into prognostic         |
|                                     |                                   | variables and output that way -- much    |
|                                     |                                   | easier to use with ``ncview``, for       |
|                                     |                                   | example. [default: ``.false.``]          |
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
| noah_state_variables                | character(len=32)::               | The list of variable names in the NOAH   |
|                                     | dimension(2,40)                   | restart file to use to create the DART   |
|                                     |                                   | state vector and their corresponding     |
|                                     |                                   | DART kind. [default: see example below]  |
+-------------------------------------+-----------------------------------+------------------------------------------+

Example
-------

.. code-block:: fortran

   &model_nml
      noah_netcdf_file             = 'restart.nc',
      assimilation_period_days     = 0,
      assimilation_period_seconds  = 3600,
      model_perturbation_amplitude = 0.0,
      output_state_vector          = .false.,
      debug                        = 0,
      noah_state_variables         = 'SOIL_T',   'KIND_SOIL_TEMPERATURE',
                                       'SOIL_M',   'KIND_SOIL_MOISTURE',
                                       'SOIL_W',   'KIND_SOIL_LIQUID_WATER',
                                       'SKINTEMP', 'KIND_SKIN_TEMPERATURE',
                                       'SNODEP',   'KIND_SNOW_THICKNESS',
                                       'WEASD',    'KIND_SNOW_WATER',
                                       'CANWAT',   'KIND_CANOPY_WATER',
                                       'QFX',      'KIND_LATENT_HEAT_FLUX',
                                       'HFX',      'KIND_SENSIBLE_HEAT_FLUX',
                                       'GRDFLX',   'KIND_GROUND_HEAT_FLUX'
   /

The second column of ``noah_state_variables`` needs some explanation. The DART
'KIND's match what the ``model_mod`` knows how to interpolate, so you can't
just add a new kind and expect it to work. There is a complex interplay between
``obs_def_mod`` and ``preprocess``, and ``model_mod`` that defines what KINDs
are supported. There is only a single KIND that works with each variable and
the example shows the current KINDs. Support for these KINDs was provided by
running ``preprocess`` with the following namelist settings:

.. code-block::

   &preprocess_nml
      input_obs_kind_mod_file = '../../../obs_kind/DEFAULT_obs_kind_mod.F90',
      output_obs_kind_mod_file = '../../../obs_kind/obs_kind_mod.f90',
      input_obs_def_mod_file = '../../../obs_def/DEFAULT_obs_def_mod.F90',
      output_obs_def_mod_file = '../../../obs_def/obs_def_mod.f90',
      input_files              = '../../../obs_def/obs_def_tower_mod.f90'
   /

NOAHLSM_OFFLINE NAMELIST
------------------------

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
|                                     |                                   |  [default: ``'restart.nc'``]             |
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
===========

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

Output Files
============

+-----------------------------------+-----------------------------------------+
| *&model_nml:noah_netcdf_filename* | the updated RESTART file containing the |
|                                   | NOAH model state.                       |
+-----------------------------------+-----------------------------------------+
| True_State.nc                     | the time-history of the "true" model    |
|                                   | state from an OSSE                      |
+-----------------------------------+-----------------------------------------+
| Prior_Diag.nc                     | the time-history of the model state(s)  |
|                                   | before assimilation                     |
+-----------------------------------+-----------------------------------------+
| Posterior_Diag.nc                 | the time-history of the model state(s)  |
|                                   | after assimilation                      |
+-----------------------------------+-----------------------------------------+
| dart_log.out [default name]       | the run-time diagnostic output          |
+-----------------------------------+-----------------------------------------+
| dart_log.nml [default name]       | the record of all the namelists         |
|                                   | actually USED - contains the default    |
|                                   | values                                  |
+-----------------------------------+-----------------------------------------+

Terms of Use
============

|Copyright| University Corporation for Atmospheric Research

Licensed under the `Apache License, Version 2.0
<http://www.apache.org/licenses/LICENSE-2.0>`__. Unless required by applicable
law or agreed to in writing, software distributed under this license is
distributed on an "as is" basis, without warranties or conditions of any kind,
either express or implied.

.. |Copyright| unicode:: 0xA9 .. copyright sign
