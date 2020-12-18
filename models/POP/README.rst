##########
POP README
##########

Contents
========

#. `Overview`_
#. `Observations`_
#. `Namelist`_
#. `Files`_
#. `Terms of Use`_

Overview
========

The **Parallel Ocean Program (POP)** may be used with the **Data Assimilation
Research Testbed (DART)**. Both main variants - `LANL POP
<https://climatemodeling.science.energy.gov/projects/climate-ocean-and-sea-ice-modeling-cosim>`_
and `CESM2.0 POP2 <https://www.cesm.ucar.edu/models/cesm2/ocean/>`_ have been
tested under the Lanai framework, only use with **CESM** is supported under the 
Manhattan framework.

We are assimilating salinity and temperature observations from either the World
Ocean Database (2005) or (2009).

The following POP variables are extracted from the POP **netCDF** restart files
and are conveyed to DART: *SALT_CUR*, *TEMP_CUR*, *UVEL_CUR*, *VVEL_CUR*, and
*PSURF_CUR*. These variables are then adjusted to be consistent with real
observations and stuffed back into the same netCDF restart files. Since DART is
an ensemble algorithm, there are multiple restart files for a single restart
time: one for each ensemble member. Creating the initial ensemble of ocean
states is an area of active research. At present, it may be sufficient to use a
climatological ensemble; e.g., using the restarts for '1 January 00Z' from 50
consecutive years.

Experience has shown that having a paired (unique) atmospheric forcing maintains
the ensemble spread better than simply forcing all the ocean ensemble members
with one single atmospheric state.

DART reads the grid information for POP from the files specified in POP's
``&grid_nml``. When DART is responsible for starting/stopping POP, the 
information is conveyed through POP's ``&time_manager_nml``.


CESM1 POP2
----------

was tested and run in production on NCAR's bluefire computer from IBM. This
implementation is a significant departure from the DART 'business as usual'
model in that DART is not responsible for advancing the model - in this case,
ALL of CESM! Instead, the CESM Interactive Ensemble facility is used to manage
the ensemble and the Flux Coupler is responsible for stopping POP at the times
required to perform an assimilation. DART simply runs 'end-to-end' at every
assimilation time, while CESM runs continuously.

This is a complete role-reversal from the normal DART operation but was
relatively simple to implement because CESM had infrastructure to exploit.

Several modifications to CESM CASEROOT scripts will be required and will be
documented more later in this document. The
``DART/models/POP/shell_scripts/cesm1_x/assimilate.csh`` script is inserted into
the CESM run script. The Flux Coupler stops POP every midnight and all the
observations within +/- 12 hours are assimilated. The observation sequence files
have been parsed into 'daylong' chunks and have names derived from the date to
facilitate manipulation in the UNIX shell.

The DART components were built with the following settings in ``mkmf.template``:

.. code-block:: perl

   MPIFC = mpxlf95_r
   MPILD = mpxlf95_r
   FC = xlf90_r
   LD = xlf90_r
   INCS = -I/usr/local/lib64/r4i4 -I/usr/local/include
   LIBS = -L/usr/local/lib64/r4i4 -lnetcdf
   FFLAGS = -qsuffix=f=f90:cpp=F90 -q64 -b64 -qarch=auto -qmaxmem=-1 -O2 $(INCS)
   LDFLAGS = $(FFLAGS) $(LIBS)

LANL-POP
--------

Has not been tested with this version of DART and must be verified.

It is invoked the same way as any other high-order model.

.. important:: 

   This interface was tested with the LANL/POP 2_0_1 version of POP ... but
   STILL CANNOT BE USED for assimilation until the POP code is modified to do a
   forward euler timestep for an 'assimilation' restart. DART is invoked and POP
   is started/stopped multiple times. It was checked in the gx3v5 geometry and
   POP was built in the 'default' configuration: one which requires no forcing
   files, no boundary conditions, etc., so I have no idea what to expect when
   confronting this with real observations ...

.. code-block:: bash

   setenv ARCHDIR linux
   gmake OPTIMIZE=no COUPLED=no
         
Given the wide range of input files and modes for running POP - the DART scripts
will surely have to be modified to accomodate moving the boundary/forcing files
required for different usage patterns.

There are several scripts in the ``DART/models/POP/shell_scripts`` directory
that are employed when using DART to assimilate with a LANL/POP model:

- ``advance_model.csh``
- ``run_perfect_model_obs.batch``, and
- ``run_filter.batch``.

The DART components compile and run on our Intel-based cluster running SLES10
with the ifort 10.1 20090203 compiler with the following flags (the value of
NETCDF was appropriate for our system in ``mkmf.template``:

.. code-block:: perl

   MPIFC = mpif90
   MPILD = mpif90
   FC = ifort
   LD = ifort
   INCS = -I$(NETCDF)/include
   LIBS = -L$(NETCDF)/lib -lnetcdf -lmkl -lmkl_lapack -lguide -lpthread
   FFLAGS = -O0 -fpe0 -vec-report0 -assume byterecl $(INCS)
   LDFLAGS = $(FFLAGS) $(LIBS)
         

Intel-based machines are natively little-endian, so I like to append a ".le"
suffix on all binary files.

On our machine, with the openmpi framework, it is necessary to specify
``input.nml:&mpi_utilities_nml:reverse_task_layout = .true.,`` to be able to
simultaneously run (2) MPI programs on the same set of nodes.

Observations
============

The observations come from the `World Ocean Database 2005
<http://www.nodc.noaa.gov/OC5/WOD05/pr_wod05.html>`_ and are processed by DART
routines in the ``$DART/observations/obs_converters/`` directory.

Converting between DART files and POP restart files
---------------------------------------------------

**Is not needed.** DART natively reads and writes netCDF files.

Generating the initial ensemble
-------------------------------

Creating the initial ensemble of ocean states is an area of active research.
The POP model cannot take one single model state and generate its own ensemble.

The ensemble has to come from 'somewhere else'. At present, it may be sufficient
to use a climatological ensemble; e.g., using the POP restarts for '1 January
00Z' from 50 consecutive years from a hindcast experiment.

By The Way
~~~~~~~~~~

Experience has shown that having a paired (unique) atmospheric forcing maintains
the ensemble spread better than simply forcing all the ocean ensemble members
with one single atmospheric state.

Generating a set of observations for a 'perfect model' experiment using the LANL/POP executable and scripts
-----------------------------------------------------------------------------------------------------------

A perfectly sensible approach to get to know the system would be to try to

#. assimilate data for the first assimilation period and stop. Do not advance
   the model at all. The filter namelist can control all of this and you do
   not need to have a working *advance_model.csh* script, or even a working
   ocean model (as long as you have input data files).
#. advance the model first and then assimilate data for the first
   assimilation period and stop.
#. advance, assimilate and advance again. This tests the whole DART facility.

I always like running something akin to a 'perfect model' experiment to
start. Since I have not come up with a good way to perturb a single model
state to generate an ensemble, here's the next best thing. The details for
running each program are covered in their own documentation.

#. Create a set of initial conditions for DART by running one instance of POP
   for a very long time and saving restart files 'every so often'. Use one of
   these as the initial condition for ``perfect_model_obs`` and the rest as the
   ensemble for the assimilation experiment. Since no one in their right mind
   would use a high-resolution model for a proof-of-concept case (hint,
   hint), running a low-resolution model for a 'very long time' should not be
   a problem.
#. create a TINY (i.e. 1) set of 'perfect' observations in the normal
   fashion using ``create_obs_sequence`` and then use
   ``create_fixed_network_seq`` to create an empty observation sequence file
   (usually called ``obs_seq.in``). The programs will prompt you for all the
   information they require. Read their documentation if necessary.
#. break the ``pop_in`` namelist that comes with POP into two pieces - one
   called ``pop_in.part1``, that contains the ``&time_manager_nml`` and put the
   rest in ``pop_in.part2``. The ``&time_manager_nml`` will be repeatedly
   updated as the POP model is repeatedly called by ``advance_model.csh``.
#. modify ``POP/work/input.nml`` as needed.
#. modify ``DART/models/POP/shell_scriptsrun_perfect_model_obs.batch`` to
   reflect the location of your DART directory, the POP directory, and which
   POPFILE to use as the initial condition.
#. Run the experiment and populate the observation sequence file by
   executing/submitting the script
   ``DART/models/POP/shell_scripts/run_perfect_model_obs.batch``. The script
   may require some modification, but not much. Please let me know if I can
   improve the readability or comments. ``run_perfect_model_obs.batch`` runs
   ``perfect_model_obs``.
#. ``run_filter.batch`` runs ``filter`` in a similar fashion. I have not
   finished the documentation for this yet.

Exploring the Output
--------------------

Is pretty much like any other model. The netCDF files have the model prognostic
variables before and after the assimilation. There are Matlab® scripts for
perusing the netCDF files in the ``DART/matlab`` directory. There are Matlab®
scripts for exploring the performance of the assimilation in observation-space
(after running ``obs_diag``) to explore the ``obs_seq.final`` file) - use the
scripts starting with ``plot_``, i.e. ``DART/diagnostics/matlab/plot_*.m*``.

As always, there are some model-specific items you should know about in
``DART/models/POP/matlab``, and ``DART/models/POP/shell_scripts``.

It is also worthwhile to convert your ``obs_seq.final`` file to a netCDF format
obs_sequence file with ``obs_seq_to_netcdf`` in
``DART/assimilation_code/programs/obs_seq_to_netcdf/``.

Namelist
========

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand ``&`` and terminate with a slash ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist.

The variables and their default values are listed here:

.. code-block:: fortran

   &model_nml
      assimilation_period_days     = -1
      assimilation_period_seconds  = -1
      model_perturbation_amplitude = 0.2
      binary_grid_file_format      = 'big_endian'
      debug                        = 0,
      model_state_variables        = 'SALT_CUR ', 'QTY_SALINITY             ', 'UPDATE',
                                       'TEMP_CUR ', 'QTY_POTENTIAL_TEMPERATURE', 'UPDATE',
                                       'UVEL_CUR ', 'QTY_U_CURRENT_COMPONENT  ', 'UPDATE',
                                       'VVEL_CUR ', 'QTY_V_CURRENT_COMPONENT  ', 'UPDATE',
                                       'PSURF_CUR', 'QTY_SEA_SURFACE_PRESSURE ', 'UPDATE'
   /

This namelist provides control over the assimilation period for the model. All
observations within (+/-) half of the assimilation period are assimilated. The
assimilation period is the minimum amount of time the model can be advanced, and
checks are performed to ensure that the assimilation window is a multiple of the
ocean model dynamical timestep.

+-------------------------------------+-------------------+------------------------------------------------------------+
| Item                                | Type              | Description                                                |
+=====================================+===================+============================================================+
| ``assimilation_period_days``        | integer           | The number of days to advance the model for each           | 
|                                     |                   | assimilation. If both ``assimilation_period_days`` and     |
|                                     |                   | ``assimilation_period_seconds`` are ≤ 0; the value of the  | 
|                                     |                   | POP namelist variables ``restart_freq`` and                |
|                                     |                   | ``restart_freq_opt`` are used to determine the             |
|                                     |                   | assimilation period.                                       |
|                                     |                   |                                                            |
|                                     |                   | *WARNING:* in the CESM framework, the ``restart_freq`` is  |
|                                     |                   | set to a value that is not useful so DART defaults to 1    |
|                                     |                   | day - even if you are using POP in the LANL framework.     |
+-------------------------------------+-------------------+------------------------------------------------------------+
| ``assimilation_period_seconds``     | integer           | In addition to ``assimilation_period_days``, the number    |
|                                     |                   | of seconds to advance the model for each assimilation.     |
|                                     |                   | Make sure you read the description of                      |
|                                     |                   | ``assimilation_period_days*.                               |
+-------------------------------------+-------------------+------------------------------------------------------------+
| ``model_perturbation_amplitude``    | real(r8)          | Reserved for future use.                                   |
+-------------------------------------+-------------------+------------------------------------------------------------+
| ``binary_grid_file_format``         | character(len=32) | The POP grid files are in a binary format. Valid values    |
|                                     |                   | are ``native``, ``big_endian``, or ``little_endian``.      |
|                                     |                   | Modern versions of Fortran allow you to specify the        |
|                                     |                   | endianness of the file you wish to read when they are      |
|                                     |                   | opened as opposed to needing to set a compiler switch or   |
|                                     |                   | environment variable.                                      |
+-------------------------------------+-------------------+------------------------------------------------------------+
| ``debug``                           | integer           | The switch to specify the run-time verbosity.              |
|                                     |                   |                                                            |
|                                     |                   | - ``0`` is as quiet as it gets.                            |
|                                     |                   | - ``> 1`` provides more run-time messages.                 |
|                                     |                   | - ``> 5`` provides ALL run-time messages.                  |
|                                     |                   |                                                            |
|                                     |                   | All values above ``0`` will also write a netCDF file of    |
|                                     |                   | the grid information and perform a grid interpolation      |
|                                     |                   | test.                                                      |
+-------------------------------------+-------------------+------------------------------------------------------------+
| ``model_state_variables``           | character(:,3)    | Strings that associate POP variables with a DART kind and  |
|                                     |                   | whether or not to write the updated values to the restart  |
|                                     |                   | files. These variables will be read from the POP restart   |
|                                     |                   | file and modified by the assimilation. Some (perhaps all)  |
|                                     |                   | will be used by the forward observation operators. If the  |
|                                     |                   | 3rd column is 'UPDATE', the output files will have the     |
|                                     |                   | modified (assimilated,posterior) values. If the 3rd        |
|                                     |                   | column is 'NO_COPY_BACK', that variable will not be        |
|                                     |                   | written to the restart files. **The DART diagnostic files  |
|                                     |                   | will always have the (modified) posterior values.**        |
|                                     |                   | Diagnostic variables that are useful for the calculation   |
|                                     |                   | of the forward observation operator but have no impact on  |
|                                     |                   | the forecast trajectory of the model could have a value of |
|                                     |                   | ``NO_COPY_BACK``. The DART kind must be one found in the   |
|                                     |                   | ``obs_kind_mod.f90`` source code file kept in              |
|                                     |                   | ``DART/assimilation_code/modules/observations/`` **AFTER** |
|                                     |                   | it gets built by ``preprocess``. Most of the ocean         |
|                                     |                   | observation kinds are specified within the                 |
|                                     |                   | ``obs_def_ocean_mod.f90`` source code file kept in         |
|                                     |                   | ``DART/observations/forward_operators/``, so it should be  |
|                                     |                   | specified in the ``&preprocess_nml:input_files``           |
|                                     |                   | variable.                                                  |
+-------------------------------------+-------------------+------------------------------------------------------------+

.. code-block:: fortran

      namelist /time_manager_nml/  runid, stop_option, stop_count, &
             time_mix_opt, fit_freq, time_mix_freq, dt_option, dt_count, &
             impcor, laccel, accel_file, dtuxcel, allow_leapyear, &
             date_separator, iyear0, imonth0, iday0, ihour0, iminute0, isecond0

This namelist is read in a file called ``pop_in``. This namelist is the same
one that is used by the ocean model and is used to control the integration
length of POP.

It is unimportant for the CESM/POP experiments but is critically important for
the LANL/POP experiments. The values are explained in full in the POP
documentation. The DART code reads the namelist and simply overwrites several
values with the new time integration information. All the other values are
unchanged.

``dart_to_pop`` writes out a new ``&time_manager_nml`` in ``pop_in.DART`` if the
DART state being converted has the ``'advance_to_time'`` record in it. This is
the case during the middle of a DART experiment, but is not typically
encountered if one is working with DART 'initial conditions' or 'restart'
files. The ``pop_in.DART`` must be concatenated with the other namelists
needed by POP into a file called ``pop_in``. We have chosen to store the
other namelists (which contain static information) in a file called
``pop_in.part2``. Initially, the ``time_manager_nml`` is stored in a companion
file called ``pop_in.part1`` and the two files are concatenated into the
expected ``pop_in`` - then, during the course of an assimilation experiment,
DART keeps writing out a new ``time_manager_nml`` with new integration
information - which gets appended with the static information in
``pop_in.part2``.

If you are running the support programs in a standalone fashion (as you
might if you are converting restart files into an intial ensemble), the
'valid time' of the model state comes from the restart file - NOT - the
namelist. You can always patch the times in the headers with
``restart_file_utility``.

Only the namelist variables of interest to DART are discussed. All other
namelist variables are ignored by DART - but mean something to POP.

+-------------------------------------+-----------------------------------+------------------------------------------+
| Item                                | Type                              | Description                              |
+=====================================+===================================+==========================================+
| ``stop_option``                     | character [default: ``'nday'``]   | The units for ``stop_count``.            |
+-------------------------------------+-----------------------------------+------------------------------------------+
| ``stop_count``                      | integer [default: ``1``]          | The duration of the model integration.   |
|                                     |                                   | The units come from ``stop_option``.     |
+-------------------------------------+-----------------------------------+------------------------------------------+

Example Namelist
----------------

.. code-block:: fortran

   &time_manager_nml
      runid          = 'gx3v5'
      stop_option    = 'nday'
      stop_count     = 1
      time_mix_opt   = 'avgfit'
      fit_freq       = 1
      time_mix_freq  = 17
      dt_option      = 'auto_dt'
      dt_count       = 1
      impcor         = .true.
      laccel         = .false.
      accel_file     = 'unknown_accel_file'
      dtuxcel        = 1.0
      allow_leapyear = .true.
      iyear0         = 2000
      imonth0        = 1
      iday0          = 1
      ihour0         = 0
      iminute0       = 0
      isecond0       = 0
      date_separator = '-'
   /

Files
=====

+-------------------------------+---------------------------------------------+
| filename                      | purpose                                     |
+===============================+=============================================+
| input.nml                     | to read the model_mod namelist              |
+-------------------------------+---------------------------------------------+
| pop_in                        | to read the model_mod namelist              |
+-------------------------------+---------------------------------------------+
| pop.r.nc                      | provides grid dimensions and 'valid_time'   |
|                               | of the model state                          |
+-------------------------------+---------------------------------------------+
| *&grid_nml* "horiz_grid_file" | contains the values of the horizontal grid  |
+-------------------------------+---------------------------------------------+
| *&grid_nml* "vert_grid_file"  | contains the number and values of the       |
|                               | vertical levels                             |
+-------------------------------+---------------------------------------------+
| true_state.nc                 | the time-history of the "true" model state  |
|                               | from an OSSE                                |
+-------------------------------+---------------------------------------------+
| preassim.nc                   | the time-history of the model state before  |
|                               | assimilation                                |
+-------------------------------+---------------------------------------------+
| analysis.nc                   | the time-history of the model state after   |
|                               | assimilation                                |
+-------------------------------+---------------------------------------------+
| dart_log.out [default name]   | the run-time diagnostic output              |
+-------------------------------+---------------------------------------------+
| dart_log.nml [default name]   | the record of all the namelists actually    |
|                               | USED - contains the default values          |
+-------------------------------+---------------------------------------------+

Terms of Use
============

DART software - Copyright UCAR. This open source software is provided by UCAR,
"as is", without charge, subject to all terms of use at
http://www.image.ucar.edu/DAReS/DART/DART_download

.. |DART project logo| image:: ../../docs/images/Dartboard7.png
   :height: 70px
