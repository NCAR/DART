GITM
====

.. attention::

   ``GITM`` works with versions of DART *before* Manhattan (9.x.x) and has yet to be updated.
   If you are interested in using ``GITM`` with more recent versions of DART, contact DAReS 
   staff to assess the feasibility of an update.
   Until that time, you should consider this documentation as out-of-date.

   A preliminary Manhattan-compliant interface to GITM exists and has been used for science 
   applications but has not yet been released for public use.

Overview
--------

The `Global Ionosphere Thermosphere Model
(GITM) <http://ccmc.gsfc.nasa.gov/models/modelinfo.php?model=GITM>`__ is a
3-dimensional spherical code that models the Earth's thermosphere and ionosphere
system using a stretched grid in latitude and altitude.

The **GITM** interface for **Data Assimilation Research Testbed (DART)** is
under development. If you wish to use GITM, you are urged to contact us. The
original scripts were configured to run on the University of Michigan machine
**NYX** using the Portable Batch System (**PBS**). We have attempted to extend
the scripts to work with both PBS and LSF and are only partway through the
process.

**DART does not come with the GITM code.** You need to get that on your own.
The normal procedure of building GITM creates some resource files that are
subsequently needed by DART - just to compile. These include:

#. ``models/gitm/GITM2/src/ModConstants.f90``
#. ``models/gitm/GITM2/src/ModEarth.f90``
#. ``models/gitm/GITM2/src/ModKind.f90``
#. ``models/gitm/GITM2/src/ModOrbital.f90``
#. ``models/gitm/GITM2/src/ModSize.f90``
#. ``models/gitm/GITM2/src/ModTime.f90``
#. ``models/gitm/GITM2/src/time_routines.f90``


GITM uses binary files for their restart mechanisms, so no metadata is available
to confirm the number and order of fields in the file. Care must be used to make
sure the namelist-controlled set of variables to be included in the DART state
vector is consistent with the restart files. Each variable must also correspond
to a DART "KIND"; required for the DART interpolate routines.

For example, this configuration of ``input.nml`` is nowhere close to being
correct:

.. code-block:: fortran

   &model_nml
      gitm_state_variables = 'Temperature',      'QTY_TEMPERATURE',
                             'eTemperature',     'QTY_TEMPERATURE_ELECTRON',
                             'ITemperature',     'QTY_TEMPERATURE_ION',
                             'iO_3P_NDensityS',  'QTY_DENSITY_NEUTRAL_O3P',
                             'iO2_NDensityS',    'QTY_DENSITY_NEUTRAL_O2',
                             'iN2_NDensityS',    'QTY_DENSITY_NEUTRAL_N2',
                                  ...                    ...
   /
      

These variables are then adjusted to be consistent with observations and stuffed
back into the same netCDF restart files. Since DART is an ensemble algorithm,
there are multiple restart files for a single restart time: one for each
ensemble member. Creating the initial ensemble of states is an area of active
research.

DART reads grid information for GITM from several sources. The ``UAM.in`` file
specifies the number of latitudes/longitudes per block, and the number of blocks
comes from the ``GITM2/src/ModSize.f90`` module. Internal to the DART code, the
following variables exist:

+-----------------+---------------+------------------------------------------+
| Item            | Type          | Description                              |
+=================+===============+==========================================+
| LON(:)          | real(r8)      | longitude array [0, 360)                 |
+-----------------+---------------+------------------------------------------+
| LAT(:)          | real(r8)      | latitude array (-90,90)                  |
+-----------------+---------------+------------------------------------------+
| ALT(:)          | real(r8)      | altitude array (0,~inf)                  |
+-----------------+---------------+------------------------------------------+
| NgridLon        | integer       | the length of the longitude array        |
+-----------------+---------------+------------------------------------------+
| NgridLat        | integer       | the length of the latitude array         |
+-----------------+---------------+------------------------------------------+
| NgridAlt        | integer       | the length of the altitude array         |
+-----------------+---------------+------------------------------------------+

Compiling
---------

GITM has been sucessfully tested with DART using the ``gfortran`` compiler,
version ``4.2.3``. The DART components were built with the following
``mkmf.template`` settings.

.. code-block::

   FC = gfortran
   LD = gfortran
   NETCDF = /Users/thoar/GNU
   INCS = -I${NETCDF}/include
   LIBS = -L${NETCDF}/lib -lnetcdf -lcurl -lhdf5_hl -lhdf5 -lz  -lm
   FFLAGS = -O0 -fbounds-check -frecord-marker=4 -ffpe-trap=invalid $(INCS)
   LDFLAGS = $(FFLAGS) $(LIBS)
         
Converting Between DART Files and GITM Restart Files
----------------------------------------------------

The binary GITM files contain no metadata, so care is needed when converting
between DART state variables and GITM files.

There are two programs - both require the list of GITM variables to use in the
DART state vector: the ``&model_nml:gitm_state_variables`` variable in the
``input.nml`` file.

+-------------------------+--------------------------------------+
| ``gitm_to_dart.f90``    | converts a set of GITM restart files |
|                         | (there is one restart file per       |
|                         | block) *bxxxx.rst* into a            |
|                         | DART-compatible file normally called |
|                         | *dart_ics* . We usually wind up      |
|                         | linking to this static filename.     |
+-------------------------+--------------------------------------+
| ``dart_to_gitm.f90``    | inserts the DART output into         |
|                         | existing GITM restart files. There   |
|                         | are two different types of DART      |
|                         | output files, so there is a namelist |
|                         | option to specify if the DART file   |
|                         | has two time records or just one. If |
|                         | there are two, the first one is the  |
|                         | 'advance_to' time, followed by the   |
|                         | 'valid_time' of the ensuing state.   |
|                         | If there is just one, it is the      |
|                         | 'valid_time' of the ensuing state.   |
|                         | *dart_to_gitm* determines the GITM   |
|                         | restart file name from the           |
|                         | *input.nml*                          |
|                         | *model_nml:gitm_restart_dirname*. If |
|                         | the DART file contains an            |
|                         | 'advance_to' time, *dart_to_gitm*    |
|                         | creates a                            |
|                         | *DART_GITM_time_control.txt* file    |
|                         | which can be used to control the     |
|                         | length of the GITM integration.      |
+-------------------------+--------------------------------------+

Simple Test
-----------

The simplest way to test the converter is to compile GITM and run a single
model state forward using ``work/clean.sh``. To build GITM ... download GITM
and unpack the code into ``DART/models/gitm/GITM2`` and run the following 
commands:

.. code-block:: bash

   $ cd models/gitm/GITM2
   $ ./Config.pl -install -compiler=ifortmpif90 -earth
   $ make
   $ cd ../work
   $ ./clean.sh 1 1 0 150.0 170.0 1.0

Namelist
--------

We adhere to the F90 standard of starting a namelist with an ampersand ``&``
and terminating with a slash ``/`` for all our namelist input. Character
strings that contain a ``/`` **must** be enclosed in quotes to prevent them
from prematurely terminating the namelist.

This namelist is read from a file called ``input.nml``. This namelist provides
control over the assimilation period for the model. All observations within
(+/-) half of the assimilation period are assimilated. The assimilation period
is the minimum amount of time the model can be advanced, and checks are
performed to ensure that the assimilation window is a multiple of the model
dynamical timestep.

Sample input.nml Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: fortran

   # The list of variables to put into the state vector is here:
   # The definitions for the DART kinds are in DART/observations/forward_operators/obs_def*f90
   # The order doesn't matter to DART. It may to you.

   &model_nml
      gitm_restart_dirname         = 'advance_temp_e1/UA/restartOUT',
      assimilation_period_days     = 0,
      assimilation_period_seconds  = 1800,
      model_perturbation_amplitude = 0.2,
      output_state_vector          = .false.,
      calendar                     = 'Gregorian',
      debug                        = 0,
      gitm_state_variables  = 'Temperature',            'QTY_TEMPERATURE',
                              'eTemperature',           'QTY_TEMPERATURE_ELECTRON',
                              'ITemperature',           'QTY_TEMPERATURE_ION',
                              'iO_3P_NDensityS',        'QTY_DENSITY_NEUTRAL_O3P',
                              'iO2_NDensityS',          'QTY_DENSITY_NEUTRAL_O2',
                              'iN2_NDensityS',          'QTY_DENSITY_NEUTRAL_N2',
                              'iN_4S_NDensityS',        'QTY_DENSITY_NEUTRAL_N4S',
                              'iNO_NDensityS',          'QTY_DENSITY_NEUTRAL_NO',
                              'iN_2D_NDensityS',        'QTY_DENSITY_NEUTRAL_N2D',
                              'iN_2P_NDensityS',        'QTY_DENSITY_NEUTRAL_N2P',
                              'iH_NDensityS',           'QTY_DENSITY_NEUTRAL_H',
                              'iHe_NDensityS',          'QTY_DENSITY_NEUTRAL_HE',
                              'iCO2_NDensityS',         'QTY_DENSITY_NEUTRAL_CO2',
                              'iO_1D_NDensityS',        'QTY_DENSITY_NEUTRAL_O1D',
                              'iO_4SP_IDensityS',       'QTY_DENSITY_ION_O4SP',
                              'iO2P_IDensityS',         'QTY_DENSITY_ION_O2P',
                              'iN2P_IDensityS',         'QTY_DENSITY_ION_N2P',
                              'iNP_IDensityS',          'QTY_DENSITY_ION_NP',
                              'iNOP_IDensityS',         'QTY_DENSITY_ION_NOP',
                              'iO_2DP_IDensityS',       'QTY_DENSITY_ION_O2DP',
                              'iO_2PP_IDensityS',       'QTY_DENSITY_ION_O2PP',
                              'iHP_IDensityS',          'QTY_DENSITY_ION_HP',
                              'iHeP_IDensityS',         'QTY_DENSITY_ION_HEP',
                              'ie_IDensityS',           'QTY_DENSITY_ION_E',
                              'U_Velocity_component',   'QTY_VELOCITY_U',
                              'V_Velocity_component',   'QTY_VELOCITY_V',
                              'W_Velocity_component',   'QTY_VELOCITY_W',
                              'U_IVelocity_component',  'QTY_VELOCITY_U_ION',
                              'V_IVelocity_component',  'QTY_VELOCITY_V_ION',
                              'W_IVelocity_component',  'QTY_VELOCITY_W_ION',
                              'iO_3P_VerticalVelocity', 'QTY_VELOCITY_VERTICAL_O3P',
                              'iO2_VerticalVelocity',   'QTY_VELOCITY_VERTICAL_O2',
                              'iN2_VerticalVelocity',   'QTY_VELOCITY_VERTICAL_N2',
                              'iN_4S_VerticalVelocity', 'QTY_VELOCITY_VERTICAL_N4S',
                              'iNO_VerticalVelocity',   'QTY_VELOCITY_VERTICAL_NO',
                              'f107',                   'QTY_1D_PARAMETER',
                              'Rho',                    'QTY_DENSITY',
         /

Description of Each Term in the Namelist
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+-------------------------------------+-----------------------------------+------------------------------------------+
| Item                                | Type                              | Description                              |
+=====================================+===================================+==========================================+
| gitm_restart_dirname                | character(len=256)                | The name of the directory containing the |
|                                     |                                   | GITM restart files and runtime control   |
|                                     |                                   | information.                             |
+-------------------------------------+-----------------------------------+------------------------------------------+
| assimilation_period_days            | integer                           | The number of days to advance the model  |
|                                     |                                   | for each assimilation.                   |
+-------------------------------------+-----------------------------------+------------------------------------------+
| assimilation_period_seconds         | integer                           | In addition to                           |
|                                     |                                   | ``assimilation_period_days`` the number  |
|                                     |                                   | of seconds to advance the model for each |
|                                     |                                   | each assimilation.                       |
+-------------------------------------+-----------------------------------+------------------------------------------+
| model_perturbation_amplitude        | real(r8)                          | Reserved for future use.                 |
+-------------------------------------+-----------------------------------+------------------------------------------+
| output_state_vector                 | logical                           | The switch to determine the form of the  |
|                                     |                                   | of the state vector in the output netCDF |
|                                     |                                   | files.                                   |
|                                     |                                   | If ``.true.`` the state vector will be   |
|                                     |                                   | output exactly as DART uses it ... one   |
|                                     |                                   | long array. If ``.false.``, the state    |
|                                     |                                   | vector is parsed into prognostic         |
|                                     |                                   | variables and output that way -- much    |
|                                     |                                   | easier to use with 'ncview', for         |
|                                     |                                   | example.                                 |
+-------------------------------------+-----------------------------------+------------------------------------------+
| calendar                            | character(len=32)                 | Character string specifying the calendar |
|                                     |                                   | being used by GITM.                      |
+-------------------------------------+-----------------------------------+------------------------------------------+
| debug                               | integer                           | The switch to specify the run-time       | 
|                                     |                                   | verbosity.                               |
|                                     |                                   |                                          |
|                                     |                                   | - ``0`` is as quiet as it gets           |
|                                     |                                   | - ``> 1`` provides more run-time         |
|                                     |                                   |   messages                               |
|                                     |                                   | - ``> 5`` provides ALL run-time messages |
+-------------------------------------+-----------------------------------+------------------------------------------+
| gitm_state_variables                | character                         | The table that relates the GITM          |
|                                     | (len=NF90_MAX_NAME)::             | variables to use to build the DART state |
|                                     | dimension(160)                    | vector, and the corresponding DART kinds |
|                                     |                                   | for those variables.                     |
+-------------------------------------+-----------------------------------+------------------------------------------+

Files
-----

+--------------------------------------+--------------------------------------+
| filename                             | purpose                              |
+======================================+======================================+
| input.nml                            | to read the model_mod namelist       |
+--------------------------------------+--------------------------------------+
| Several GITM source modules:         | provides grid dimensions, model      |
| ModConstants, ModSizeGitm, ModEarth  | state, and 'valid_time' of the model |
| ...                                  | state                                |
+--------------------------------------+--------------------------------------+
| header.rst, bNNNN.rst                | provides the 'valid_time' of the     |
|                                      | model state and the model state,     |
|                                      | respectively                         |
+--------------------------------------+--------------------------------------+
| true_state.nc                        | the time-history of the "true" model |
|                                      | state from an OSSE                   |
+--------------------------------------+--------------------------------------+
| preassim.nc                          | the time-history of the model state  |
|                                      | before assimilation                  |
+--------------------------------------+--------------------------------------+
| analysis.nc                          | the time-history of the model state  |
|                                      | after assimilation                   |
+--------------------------------------+--------------------------------------+
| dart_log.out [default name]          | the run-time diagnostic output       |
+--------------------------------------+--------------------------------------+
| dart_log.nml [default name]          | the record of all the namelists      |
|                                      | actually USED - contains the default |
|                                      | values                               |
+--------------------------------------+--------------------------------------+

References
----------

NASA's official *GITM* description can be found at their `Community Coordinated
Modeling Center website <http://ccmc.gsfc.nasa.gov/models/modelinfo.php?model=GITM>`_.
