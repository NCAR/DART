PROGRAM ``dart_to_cam``
=======================

Overview
--------

``dart_to_cam`` is the program that reads a DART restart or model advance file (e.g.
``perfect_ics, filter_ics, assim_model_state_id ...`` ). and overwrites the part of the CAM data in a single CAM restart
file (usually ``caminput.nc``) which is in the DART state vector. If you have multiple input files, you will need to
rename the output files as you create them.

The list of variables extracted from the DART state vector and exported to the CAM netCDF file is controlled by the set
of ``input.nml`` ``&model_nml:state_names_*`` variables.

If the input is a model advance file, containing 2 timestamps (the current model time and a future time the model should
run until), this program also creates a separate file named ``times`` that contains three lines: the advance-to time,
the current model time, and the number of hours to advance. These will need to be extracted and inserted in a CAM
namelist to indicate to CAM how long to run.

This program also updates the ``date`` and ``datesec`` variables in the CAM netcdf file. Generally these are identical
times since the assimilation doesn't change the time of the data, but in case the original file had a different time
that was overwritten in the state vector, it will update the time for consistency.

Conditions required for successful execution of ``dart_to_cam``:

-  a valid ``input.nml`` namelist file for DART
-  a CAM 'phis' netCDF file [default: ``cam_phis.nc``]
-  a DART restart file [default: ``dart_ics``] (read)
-  a CAM restart file [default: ``caminput.nc``] (read and written)

Since this program is called repeatedly for every ensemble member, we have found it convenient to link the DART input
and CAM restart files to the default filenames ``dart_ics`` and ``caminput.nc``). The output files may be moved or
relinked as necessary.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &dart_to_cam_nml
      dart_to_cam_input_file  = 'dart_ics',
      dart_to_cam_output_file = 'caminput.nc',
      advance_time_present    = .true.,
      /

| 

.. container::

   +-------------------------+--------------------+---------------------------------------------------------------------+
   | Item                    | Type               | Description                                                         |
   +=========================+====================+=====================================================================+
   | dart_to_cam_input_file  | character(len=128) | The name of the DART restart file containing the CAM state.         |
   +-------------------------+--------------------+---------------------------------------------------------------------+
   | dart_to_cam_output_file | character(len=128) | The name of the CAM restart netcdf file.                            |
   +-------------------------+--------------------+---------------------------------------------------------------------+
   | advance_time_present    | logical            | Set to .false. for DART initial condition and restart files. Use    |
   |                         |                    | the .true. setting for the files written by filter during a model   |
   |                         |                    | advance.                                                            |
   +-------------------------+--------------------+---------------------------------------------------------------------+

| 

Modules used
------------

::

   assim_model_mod.f90
   types_mod.f90
   threed_sphere/location_mod.f90
   model_mod.f90
   null_mpi_utilities_mod.f90
   obs_kind_mod.f90
   random_seq_mod.f90
   time_manager_mod.f90
   utilities_mod.f90

Files read
----------

-  DART namelist file; ``input.nml``
-  DART initial conditions/restart file; e.g. ``dart_ics`` (read)
-  CAM restart file; ``caminput.nc`` (read and written)
-  CAM "phis" file specified in ``&model_nml::cam_phis`` (normally ``cam_phis.nc``)

Files written
-------------

-  CAM restart file; ``caminput.nc`` (read and written)

References
----------

none
