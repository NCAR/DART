PROGRAM ``rad_3dvar_to_dart``
=============================

Overview
--------

Programs to convert MM5 3D-VAR 2.0 Radar data files into DART observation sequence files. The capability of the program
is limited to DOPPLER_RADIAL_VELOCITY and RADAR_REFLECTIVITY.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &rad_3dvar_to_dart_nml
      var_file = 'qc_radr_3dvar_2002083100.dat',
      obs_seq_out_file_name = 'obs_seq.out',
      calendar_type = 3  
   /

| 

.. container::

   ===================== ================== ==========================================================================
   Item                  Type               Description
   ===================== ================== ==========================================================================
   var_file              character(len=129) This is the name of the file containing MM5 3D-VAR 2.0 Radar observations.
   obs_seq_out_file_name character(len=129) File name for output observation sequence file.
   calendar_type         integer            Calendar type. We recommend using 3 (GREGORIAN).
   ===================== ================== ==========================================================================

| 

Modules directly used
---------------------

::

   types_mod
   obs_sequence_mod
   obs_def_mod
   obs_def/obs_def_radar_mod
   obs_kind_mod
   location/threed_sphere/location_mod
   time_manager_mod
   utilities_mod

Modules indirectly used
-----------------------

::

   assim_model_mod
   models/wrf/model_mod
   models/wrf/module_map_utils
   random_seq_mod

Files
-----

-  input namelist ; ``input.nml``
-  Input observation file; ``qc_radr_3dvar_2002083100.dat``
-  Output observation file; ``obs_seq.out``

File formats
~~~~~~~~~~~~

``input.nml`` and ``qc_radr_3dvar_2002083100.dat`` are ASCII files. ``obs_seq.out`` is either ASCII or binary, depending
on the logical write_binary_obs_sequence, which is the namelist entry for obs_sequence_mod.

References
----------

-  `3DVAR GROUP PAGE <https://www.mmm.ucar.edu/wrf-administration>`__
