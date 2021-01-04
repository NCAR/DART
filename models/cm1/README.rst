##########
CM1 README
##########

#. `Overview`_
#. `namelist.input`_
#. `About Testing CM1 and DART`_
#. `Namelist`_
#. `Terms of Use`_

Overview
========

Cloud Model 1 (CM1) version 18 (CM1r18) is compatible with the DART. CM1 is a
non-hydrostatic numerical model in Cartesian 3D coordinates designed for the
study of micro-to-mesoscale atmospheric phenomena in idealized to
semi-idealized simulations.

The CM1 model was developed and is maintained by George Bryan at the National
Center for Atmospheric Research (NCAR) Mesoscale and Microscale Meteorology
Laboratory (MMM).

The model code is freely available from the `CM1 website <http://www2.mmm.ucar.edu/people/bryan/cm1/>`_
and must be downloaded and compiled outside of DART.

This model interface and scripting support were created by Luke Madaus.
**Thanks Luke!**

namelist.input
==============

Several modifications to the CM1 namelist ``namelist.input`` are required to
produce model output compatible with DART. The values are described here and an
example is shown below.

The ``namelist.input`` file is partitioned into several distinct namelists.
These namelists are denoted ``&param0``, ``&param1``, ``&param2``, ... 
``&param13``.

These namelists start with an ampersand ``&`` and terminate with a slash ``/``.
Thus, character strings that contain a ``/`` must be enclosed in quotes to
prevent them from prematurely terminating the namelist.

Using CM1 output files as a prior ensemble state in DART requires each ensemble
member to produce a restart file in netCDF format (which requires setting
``restart_format=2`` in the ``&param9`` namelist) and these restart files must
only contain output at the analysis time (which requires setting
``restart_filetype=2`` in the ``&param9`` namelist).

Here is an example configuration of the ``&param9`` namelist in
``namelist.input``:

.. code-block:: fortran

   &param9 
     restart_format     = 2         restart needs to be netCDF
     restart_filetype   = 2         restart must be the analysis time - ONLY
     restart_file_theta = .true.    make sure theta is in restart file
     restart_use_theta  = .true.
   /

.. important::

   The only required state variable to be updated is potential temperature
   (``theta``). Thus two additional settings in the ``&param9`` namelist  --
   ``restart_file_theta = .true.`` and ``restart_use_theta = .true.`` must be
   set to ensure ``theta`` is output the CM1 restart files.

Additional state variables that have been tested within DART include:

``ua, va, wa, ppi, u0, v0, u10, v10, t2, th2, tsk, q2, psfc, qv, qc, qr, qi qs, & qg``.
  
At present, observation times are evaluated relative to the date and time
specified in the ``&param11`` namelist.

Observation locations are specified in meters relative to the domain origin as
defined the ``iorigin`` setting of ``&param2``.

About Testing CM1 and DART
==========================

There are two sets of scripts in the ``shell_scripts`` directory. Luke
contributed a set written in python, and the DART team had a set written in
csh. The csh scripts have not been tested in quite some time, so use with the
understanding that they will need work. Those csh scripts and some unfinished
python scripts reside in a ``shell_scripts/unfinished`` directory and should be
used with the understanding that they require effort on the part of the user
before the scripts will actually work.

Strategy and Instructions for Using the Python Scripts
------------------------------------------------------

A List of Prerequisites
~~~~~~~~~~~~~~~~~~~~~~~

#. CM1 is required to use netCDF restart files.
#. A collection of CM1 model states for initial conditions will be
   available.
#. There is a separate observation sequence file for each assimilation
   time.
#. The DART *input.nml* file has some required values as defined below.
#. Each time CM1 is advanced, it will start from the same filename, and
   the restart number in that filename will be 000001 - ALWAYS. That
   filename will be a link to the most current model state.

Testing a Cycling Experiment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The big picture: three scripts (``setup_filter.py``, ``run_filter.py``, and
``advance_ensemble.py``) are alternated to configure an experiment, perform an
assimilation on a set of restart files, and make the ensemble forecast. Time
management is controlled through command-line arguments.

It is required that you have generated the DART executables before you test.
The term ``{centraldir}`` refers to a filesystem and directory that will be
used to run the experiment, the working directory. ``{centraldir}`` should have
a lot of capacity, as ensemble data assimilation will require lots of disk. The
term ``{dart_dir}`` will refer to the location of the DART source code.

The data referenced in the directories (the initial ensemble, etc.) are
provided as a compressed tar file `cm1r18_3member_example_data.tar.gz
<http://www.image.ucar.edu/pub/DART/CM1/cm1r18_3member_example_data.tar.gz>`_.

You will have to download the tar file, uncompress it, and modify the scripts to
use these directories instead of the example directories in the scripts. You
will also have to compile your own cm1 executable.

#. Set some variables in both ``shell_scripts/setup_filter.py`` and
   ``shell_scripts/advance_ensemble.py`` as described below.
#. In the ``{dart_dir}/models/cm1/shell_scripts`` directory, run:
   
   .. code-block:: bash
   
      $ ./setup_filter.py -d YYYYmmDDHHMMSS -i
      
   where ``YYYYmmDDHHMMSS`` is the date and time of the first assimilation
   cycle (the ``-i`` option indicates this is the initial setup and extra work
   will be performed). This will create the working directory ``{centraldir}``,
   link in required executables, copy in the initial conditions for each member
   from some predetermined location, copy in the observation sequence file for
   this assimilation time from some predetermined location, modify namelists,
   and build a queue submission script in the ``{centraldir}``:
   ``run_filter.py``.
#. Change into ``{centraldir}`` and verify the contents of ``run_filter.py``.
   Ensure the assimilation settings in ``input.nml`` are correct. Once you
   are satisfied, submit ``run_filter.py`` to the queue to perform an
   assimilation.
#. After the assimilation job completes, check to be sure that the assimilation
   completed successfully, and the archived files requested in the
   ``setup_filter.py`` ``files_to_archive`` variable are in
   ``{centraldir}/archive/YYYYmmDDHHMMSS``.
#. Change into ``{dart_dir}/models/cm1/shell_scripts`` and advance the ensemble
   to the next assimilation time by running:
   
   .. code-block:: bash

      $ ./advance_ensemble.py -d YYYYmmDDHHMMSS -l nnnn
   
   where ``YYYYmmDDHHMMSS`` is the date of the COMPLETED analysis (the start
   time for the model) and ``nnnn`` is the length of model integration in
   seconds (the forecast length). (The forecast length option is
   specified by 'hypen ell' - the lowercase letter L, not the number one.)
   ``advance_ensemble.py`` will submit jobs to the queue to advance the
   ensemble.
#. After all ensemble members have successfully completed, run:

   .. code-block:: bash
   
      $ ./setup_filter.py -d YYYYmmDDHHMMSS
   
   where $YYYYmmDDHHMMSS$ is the **new** current analysis time. Note the $-i$
   flag is NOT used here, as we do not need to (should not need to!)
   re-initialize the entire directory structure.
#. Change into ``{centraldir}`` and run:

   .. code-block:: bash
   
      $ ``run_filter.py``
   
   to perform the assimilation.
#. Go back to step 4 and repeat steps 4-7 for each assimilation cycle
   until the end of the experiment.

Within the ``setup_filter.py`` and ``advance_ensemble.py`` scripts, the
following variables need to be set between the "BEGIN USER-DEFINED VARIABLES"
and "END USER-DEFINED VARIABLES" comment blocks:

``jobname``

   A name for this experiment, will be included in the working directory path.

``ens_size``

   Number of ensemble members.

``restart_filename``

   The filename for each ensemble member's restart. Highly recommended to leave
   this as ``cm1out_rst_000001.nc``

``window_mins``

   The assimilation window width (in minutes) for each assimilation cycle.

``copy``

   The copy command with desired flags for this system.

``link``

   The link command with desired flags for this system.

``remove``

   The remove command with desired flags for this system.

``files_to_archive``

   A list of DART output files to archive for each assimilation cycle. Note
   that any inflation files generated are automatically carried over.

``centraldir``

   Directory (which will be created if ``setup_filter.py`` is run in
   intialization mode) where the assimilation and model advances will take
   place. Should be on a system with enough space to allow for several
   assimilation cycles of archived output.

``dart_dir``

   Path to the cm1 subdirectory of DART.

``cm1_dir``

   Path to the cm1 model executable (*cm1.exe*)

``icdir``

   Path to the ensemble of initial conditions. It is assumed that within this
   directory, each ensemble member has a subdirectory (*m1*, *m2*, *m3*, ...)
   that contains:

   -  a restart file for cm1 at the desired start time and having the
      filename defined in ``restart_filename`` above
   -  a ``namelist.input`` file compatible with the generation of that
      restart file.

``obsdir``

   Path to a directory containing observation sequence files to be assimilated.
   It is assumed that the observation sequence files are named following the
   convention ``YYYYmmDDHHMMSS_obs_seq.prior``, where the date of the analysis
   time whose observations are contained in that file is the first part of the
   file name.

   ``setup_filter.py`` and ``advance_ensemble.py`` assume that mpi queue
   submissions are required to run ``cm1.exe`` and ``filter``. These variables
   control how that is handled.

``queue_system``

   The name of the queueing system

``mpi_run_command``

   The command used in a submitted script to execute an mpi task in the queue,
   including any required flags

``queue_sub_command``

   The command used to submit a script to the queue

``job_sub_info``

   A dictionary of all flags required to execute a job in the queue, with
   the key being the flag and the value being the variable. e.g. {'-P' :
   'PROJECT CODE HERE', '-W' : '00:20'}, etc.

.. note:

   ``{dart_dir}/work/input.nml`` should be modified with the desired
   assimilation settings. Some of the variables listed above will override the
   values in ``{dart_dir}/work/input.nml`` should be modified.


Namelist
========

The ``&model_nml`` namelist is read from the ``input.nml`` file. Again,
namelists start with an ampersand ``&`` and terminate with a slash ``/``.
Character strings that contain a ``/`` must be enclosed in quotes to prevent
them from prematurely terminating the namelist.

.. code-block:: fortran

   &model_nml 
      assimilation_period_days     = 0
      assimilation_period_seconds  = 21600
      model_perturbation_amplitude = 0.2
      cm1_template_file            = 'null'
      calendar                     = 'Gregorian'
      periodic_x                   = .true.
      periodic_y                   = .true.
      periodic_z                   = .false.
      debug                        = 0
      model_variables              = ' '
   /

Description of each namelist entry
----------------------------------

+------------------------------------+-----------------------+-------------------------------------------------------+
| Item                               | Type                  | Description                                           |
+====================================+=======================+=======================================================+
| assimilation_period_[days,seconds] | integer               | This specifies the width of the assimilation window.  |
|                                    |                       | The current model time is used as the center time of  |
|                                    |                       | the assimilation window. All observations in the      |
|                                    |                       | assimilation window are assimilated. BEWARE: if you   |
|                                    |                       | put observations that occur before the beginning of   |
|                                    |                       | the assimilation_period, DART will error out because  |
|                                    |                       | it cannot move the model 'back in time' to process    |
|                                    |                       | these observations.                                   |
+------------------------------------+-----------------------+-------------------------------------------------------+
| model_perturbation_amplitude       | real(r8)              | unsupported                                           |
+------------------------------------+-----------------------+-------------------------------------------------------+
| cm1_template_file                  | character(len=256)    | filename used to read the variable sizes, location    |
|                                    |                       | metadata, etc.                                        |
+------------------------------------+-----------------------+-------------------------------------------------------+
| calendar                           | character(len=256)    | Character string to specify the calendar in use.      |
|                                    |                       | Usually 'Gregorian' (since that is what the           |
|                                    |                       | observations use).                                    |
+------------------------------------+-----------------------+-------------------------------------------------------+
| model_variables                    | character(:,5)        | Strings that identify the CM1 variables, their DART   |
|                                    |                       | quantity, the minimum & maximum possible values, and  |
|                                    |                       | whether or not the posterior values should be written |
|                                    |                       | to the output file. The DART QUANTITY must be one     |
|                                    |                       | found in the *DART/obs_kind/obs_kind_mod.f90* AFTER   |
|                                    |                       | it gets built by *preprocess*.                        |
|                                    |                       |                                                       |
|                                    |                       | +--------------------------+----------------------+   |
|                                    |                       | | *model_variables(:,1)*   | Specifies the CM1    |   |
|                                    |                       | |                          | variable name in the |   |
|                                    |                       | |                          | netCDF file.         |   |
|                                    |                       | +--------------------------+----------------------+   |
|                                    |                       | | *model_variables(:,2)*   | Specifies the DART   |   |
|                                    |                       | |                          | quantity for that    |   |
|                                    |                       | |                          | variable.            |   |
|                                    |                       | +--------------------------+----------------------+   |
|                                    |                       | | *model_variables(:,3)*   | Specifies a minimum  |   |
|                                    |                       | |                          | bound (if any) for   |   |
|                                    |                       | |                          | that variable.       |   |
|                                    |                       | +--------------------------+----------------------+   |
|                                    |                       | | *model_variables(:,4)*   | Specifies a maximum  |   |
|                                    |                       | |                          | bound (if any) for   |   |
|                                    |                       | |                          | that variable.       |   |
|                                    |                       | +--------------------------+----------------------+   |
|                                    |                       | | *model_variables(:,5)*   | Specifies if the     |   |
|                                    |                       | |                          | variable should be   |   |
|                                    |                       | |                          | updated in the       |   |
|                                    |                       | |                          | restart file. The    |   |
|                                    |                       | |                          | value may be         |   |
|                                    |                       | |                          | "UPDATE" or anything |   |
|                                    |                       | |                          | else.                |   |
|                                    |                       | +--------------------------+----------------------+   |
|                                    |                       |                                                       |
+------------------------------------+-----------------------+-------------------------------------------------------+
| periodic_x                         | logical               | a value of *.true.* means the 'X' dimension is        |
|                                    |                       | periodic.                                             |
+------------------------------------+-----------------------+-------------------------------------------------------+
| periodic_y                         | logical               | a value of *.true.* means the 'Y' dimension is        |
|                                    |                       | periodic.                                             |
+------------------------------------+-----------------------+-------------------------------------------------------+
| periodic_z                         | logical               | unsupported                                           |
+------------------------------------+-----------------------+-------------------------------------------------------+
| debug                              | integer               | switch to control the amount of run-time output is    |
|                                    |                       | produced. Higher values produce more output. 0        |
|                                    |                       | produces the least.                                   |
+------------------------------------+-----------------------+-------------------------------------------------------+

.. note::

   The values above are the default values. A more realistic example is shown
   below and closely matches the values in the default ``input.nml``.

.. code-block:: fortran

   &model_nml 
      assimilation_period_days     = 0
      assimilation_period_seconds  = 60
      cm1_template_file            = 'cm1out_rst_000001.nc'
      calendar                     = 'Gregorian'
      periodic_x                   = .true.
      periodic_y                   = .true.
      periodic_z                   = .false.
      debug                        = 0
      model_variables = 'ua'   , 'QTY_U_WIND_COMPONENT'      , 'NULL', 'NULL', 'UPDATE',
                        'va'   , 'QTY_V_WIND_COMPONENT'      , 'NULL', 'NULL', 'UPDATE',
                        'wa'   , 'QTY_VERTICAL_VELOCITY'     , 'NULL', 'NULL', 'UPDATE',
                        'theta', 'QTY_POTENTIAL_TEMPERATURE' , 0.0000, 'NULL', 'UPDATE',
                        'ppi'  , 'QTY_PRESSURE'              , 'NULL', 'NULL', 'UPDATE',
                        'u10'  , 'QTY_10M_U_WIND_COMPONENT'  , 'NULL', 'NULL', 'UPDATE',
                        'v10'  , 'QTY_10M_V_WIND_COMPONENT'  , 'NULL', 'NULL', 'UPDATE',
                        't2'   , 'QTY_2M_TEMPERATURE'        , 0.0000, 'NULL', 'UPDATE',
                        'th2'  , 'QTY_POTENTIAL_TEMPERATURE' , 0.0000, 'NULL', 'UPDATE',
                        'tsk'  , 'QTY_SURFACE_TEMPERATURE'   , 0.0000, 'NULL', 'UPDATE',
                        'q2'   , 'QTY_SPECIFIC_HUMIDITY'     , 0.0000, 'NULL', 'UPDATE',
                        'psfc' , 'QTY_SURFACE_PRESSURE'      , 0.0000, 'NULL', 'UPDATE',
                        'qv'   , 'QTY_VAPOR_MIXING_RATIO'    , 0.0000, 'NULL', 'UPDATE',
                        'qc'   , 'QTY_CLOUD_LIQUID_WATER'    , 0.0000, 'NULL', 'UPDATE',
                        'qr'   , 'QTY_RAINWATER_MIXING_RATIO', 0.0000, 'NULL', 'UPDATE',
                        'qi'   , 'QTY_CLOUD_ICE'             , 0.0000, 'NULL', 'UPDATE',
                        'qs'   , 'QTY_SNOW_MIXING_RATIO'     , 0.0000, 'NULL', 'UPDATE',
                        'qg'   , 'QTY_GRAUPEL_MIXING_RATIO'  , 0.0000, 'NULL', 'UPDATE'
   /

Terms of Use
============

|Copyright| University Corporation for Atmospheric Research

Licensed under the `Apache License, Version 2.0
<http://www.apache.org/licenses/LICENSE-2.0>`__. Unless required by applicable
law or agreed to in writing, software distributed under this license is
distributed on an "as is" basis, without warranties or conditions of any kind,
either express or implied.

.. |Copyright| unicode:: 0xA9 .. copyright sign