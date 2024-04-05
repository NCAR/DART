MITgcm_ocean
============


Overview
--------

The `MIT ocean GCM <http://mitgcm.org/>`__ version 'checkpoint59a' is the foundation of this directory. It was
modified by Ibrahim Hoteit of Scripps for his use, and so it differs from the original distribution. The whole
process hinges on setting the ocean model namelist values such that it is doing a cold start for every model advance.

The data assimilation period is controlled in the ``input.nml``\ ``&model_nml`` namelist. In combination with the ocean
model dynamics timestep ``data``\ ``&PARM03:deltaTClock`` this determines the amount of time the model will advance for
each assimilation cycle.


Observations
^^^^^^^^^^^^

The forward operators for ocean observations are in ``observations/forward_operators/obs_def_oceam_mod.f90``.
:doc:`./create_ocean_obs` is the program to create a DART observation sequence from a particular ASCII file.


Converting between DART and the model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The MITgcm_ocean model reads and writes Fortran direct-access big-endian data files. To convert MITgcm_ocean files
into netcdf files for input to DART use the program ``mit_to_dart``.  To convert the netcdf output from DART to
MITgcm_ocean files use the program ``dart_to_mit``. When converting from DART to MIT and back, the following options
can be set in the ``&trans_mitdart_nml`` namelist in ``input.nml``.

.. code-block:: fortran

   &trans_mitdart_nml
     do_bgc = .false.         ! change to .true. if doing bio-geo-chemistry
     log_transform = .false.  ! change to .true. if using log_transform
     compress = .false.       ! change to .true. to compress the state vector
   /

``compress = .true.`` can be used to generate netcdf files for use with DART which has missing values (land) removed.
For some datasets this reduces the state vector size significantly. For example, the state vector size is
reduced by approximately 90% for the Red Sea. The program ``expand_netcdf`` can be used to uncompress the netcdf
file to view the data in a convenient form.


.. Warning::

   The ``trans_mit_dart_mod`` module has hardcoded record lengths, ``recl3d`` and ``recl2d``.
   Be sure to check these are correct for the compiler you are using.

  .. code-block:: fortran
  
    ! set record lengths
    recl3d = Nx*Ny*Nz*4
    recl2d = Nx*Ny*4
  

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



Getting started
^^^^^^^^^^^^^^^


#. create an initial ensemble for DART.
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

There are Matlab® scripts for perusing netCDF files in the ``DART/matlab`` directory. There are
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

