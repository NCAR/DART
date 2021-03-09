PROGRAM ``cam_to_dart``
=======================

Overview
--------

| ``cam_to_dart`` is the program that reads a CAM restart file (usually ``caminput.nc``) and creates a single DART
  output/restart file (e.g. ``perfect_ics, filter_ics, ...`` ). If you have multiple input files, you will need to
  rename the output files as you create them.
| The list of variables extracted from the CAM netCDF file and conveyed to DART is controlled by the set of
  ``input.nml`` ``&model_nml:state_names_*`` variables. The ``date`` and ``datesec`` variables in the CAM netcdf file
  are used to specify the valid time of the state vector. The time may be changed with the
  :doc:`../../assimilation_code/programs/restart_file_tool/restart_file_tool` if desired.
| Some CAM restart files are from climatological runs and have a valid time that predates the use of the Gregorian
  calendar. In such instances, the year component of the original date is changed to be a valid Gregorian year (by
  adding 1601). A warning is issued to the screen and to the logfile. Please use the
  :doc:`../../assimilation_code/programs/restart_file_tool/restart_file_tool` to change this time.
| Conditions required for successful execution of ``cam_to_dart``:

-  a valid ``input.nml`` namelist file for DART
-  a CAM 'phis' netCDF file [default: ``cam_phis.nc``]
-  a CAM restart file [default: ``caminput.nc``].

Since this program is called repeatedly for every ensemble member, we have found it convenient to link the CAM restart
files to the default input filename (``caminput.nc``). The default DART output filename is ``dart_ics`` - this may be
moved or linked as necessary.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &cam_to_dart_nml
      cam_to_dart_input_file  = 'caminput.nc',
      cam_to_dart_output_file = 'dart_ics', 
      /

| 

.. container::

   +-------------------------+--------------------+---------------------------------------------------------------------+
   | Item                    | Type               | Description                                                         |
   +=========================+====================+=====================================================================+
   | cam_to_dart_input_file  | character(len=128) | The name of the DART file containing the CAM state.                 |
   +-------------------------+--------------------+---------------------------------------------------------------------+
   | cam_to_dart_output_file | character(len=128) | The name of the DART file containing the model state derived from   |
   |                         |                    | the CAM restart file.                                               |
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
-  CAM restart file; ``caminput.nc``
-  CAM "phis" file specified in ``&model_nml::cam_phis`` (normally ``cam_phis.nc``)

Files written
-------------

-  DART initial conditions/restart file; e.g. ``dart_ics``

References
----------

none
