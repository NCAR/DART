MITgcm_ocean
============


Overview
--------

The `MIT ocean GCM <http://mitgcm.org/>`__ version 'checkpoint59a' is the foundation of this directory. It was
modified by Ibrahim Hoteit of Scripps for his use, and so it differs from the original distribution.

Since the model is highly parallelized, it can be compiled with a target number of processors in mind. From DART's
perspective, the most logical strategy is to run ``filter`` or ``perfect_model_obs`` with **async=4**: advance the
model in parallel ... one ensemble member after another. In this mode, the same set of processors are used for the
data assimilation. The performance of the parallel assimilation algorithm has been tested up through 64 processors,
and should scale well beyond that - but it remains to be quantified. The scaling for the ocean model is unknown to me,
but Ibrahim routinely runs with many more than 64 processors.

As for all DART experiments, the overall design for an experiment is this: the DART program ``filter`` will read the
initial conditions file, the observation sequence file, and the DART namelist to decide whether or not to advance the
ocean model. All of the control of the execution of the ocean model is done by DART directly. If the model needs to be
advanced, ``filter`` makes a call to the shell to execute the script ``advance_model.csh``. ``advance_model.csh`` is
ENTIRELY responsible for getting all the input files, data files, namelists, etc. into a temporary directory, running
the model, and copying the results back to the parent directory (which we call CENTRALDIR). The whole process hinges
on setting the ocean model namelist values such that it is doing a cold start for every model advance.



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

+---------------------------+-----------------------------------------------------------------------------------------+
| trans_mitdart             | converts the ocean model snapshot files into a DART netcdf files and back               |
+---------------------------+-----------------------------------------------------------------------------------------+

The data assimilation period is controlled in the ``input.nml``\ ``&model_nml`` namelist. In combination with the ocean
model dynamics timestep ``data``\ ``&PARM03:deltaTClock`` this determines the amount of time the model will advance for
each assimilation cycle.


Fortran direct-access big-endian data files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The MITgcm_ocean model uses Fortran direct-access big-endian data files. It is up to you to determine the proper
compiler flags to compile DART such that DART can read and write these files. Every compiler/architecture is different,
but we have put notes in each ``mkmf.template`` if we know how to achieve this.


Controlling the model advances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| The assimilation period is specified by two namelist parameters in the ``input.nml``\ ``&model_nml`` namelist:
  ``assimilation_period_days`` and ``assimilation_period_seconds``. Normally, all observations within (+/-) HALF of the
  total assimilation period are used in the assimilation.
| The time of the initial conditions is specified by two namelist parameters in the ``input.nml``\ &model_nml``
  namelist: ``init_time_days`` and ``init_time_seconds``; depending on the settings of these parameters, the times may
  or may not come directly from the DART initial conditions files.
| The ocean model **MUST always** start from the input datasets defined in the ``data``\ ``&PARM05`` namelist.
  Apparently, this requires ``data``\ ``&PARM03:startTime`` to be **0.0**. One of the DART support routines
  ``trans_mitdart`` converts the DART state vector to the files used in ``data``\ ``&PARM05`` and creates new
  ``data.cal``\ ``&CAL_NML`` and ``data``\ ``&PARM03`` namelists with values appropriate to advance the model to the
  desired time.
| The ocean model then advances till ``data``\ ``&PARM03:endTime`` and writes out snapshot files. ``trans_mitdart``
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


Files
-----

-  input namelist files: ``data, data.cal, input.nml``
-  input data file: ``filter_ics, perfect_ics``
-  output data files: ``[S,T,U,V,Eta].YYYYMMDD.HHMMSS.[data,meta]``

Please note that there are **many** more files needed to advance the ocean model, none of which are discussed here.

