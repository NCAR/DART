PROGRAM ``model_to_dart`` for MPAS OCN
======================================

.. attention::

   ``mpas_ocn`` was being developed with versions of DART *before* Manhattan (9.x.x) and has yet to be updated. If you are interested in
   using ``mpas_ocn`` with more recent versions of DART, contact DAReS staff to assess the feasibility of an update.
   Until that time, you should consider this documentation as out-of-date.


Overview
--------

| ``model_to_dart`` is the program that reads an MPAS OCN analysis file (nominally named ``mpas_restart.nc``) and
  creates a DART state vector file (e.g. ``perfect_ics, filter_ics, ...`` ). The MPAS analysis files have a **Time**
  UNLIMITED Dimension, which indicates there may (at some point) be more than one timestep in the file. The DART
  routines are currently designed to use the LAST timestep. If the Time dimension of length 3, we use the third
  timestep. A warning message is issued and indicates exactly the time being used.
| ``input.nml``\ ``&mpas_vars_nml`` defines the list of MPAS variables used to build the DART state vector. This
  namelist is more fully described in the :doc:`./readme` documentation. For example:

::

   &mpas_vars_nml
      mpas_state_variables = 'temperature',  'QTY_TEMPERATURE',
                             'salinity',     'QTY_SALINITY',
                             'rho',          'QTY_DENSITY',
                             'u',            'QTY_EDGE_NORMAL_SPEED',
                             'h',            'QTY_SEA_SURFACE_HEIGHT'
                             'tracer1',      'QTY_TRACER_CONCENTRATION'
      /

Conditions required for successful execution of ``model_to_dart`` are:

-  a valid ``input.nml`` namelist file for DART which contains
-  a MPAS OCN analysis file (nominally named ``mpas_analysis.nc``).

Since this program is called repeatedly for every ensemble member, we have found it convenient to link the MPAS OCN
analysis files to a static input filename (e.g. ``mpas_analysis.nc``). The default DART filename is ``dart_ics`` - this
may be moved or linked as necessary.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &model_to_dart_nml
      model_to_dart_output_file = 'dart_ics'
      /

| 

::

   &model_nml
      model_analysis_filename  = 'mpas_analysis.nc'
      /

   (partial namelist)

| 

::

   &mpas_vars_nml
      mpas_state_variables = '',
      mpas_state_bounds = '',
      /

| 

The model_to_dart namelist includes:

.. container::

   +---------------------------+--------------------+-------------------------------------------------------------------+
   | Item                      | Type               | Description                                                       |
   +===========================+====================+===================================================================+
   | model_to_dart_output_file | character(len=128) | The name of the DART file containing the model state derived from |
   |                           |                    | the MPAS analysis file.                                           |
   +---------------------------+--------------------+-------------------------------------------------------------------+

| 

Two more namelists need to be mentioned. The `model_nml <readme.html#Namelist>`__ namelist specifies the MPAS
analysis file to be used as the source. The `mpas_vars_nml <readme.html#mpas_vars_nml>`__ namelist specifies the MPAS
variables that will comprise the DART state vector.

For example:

::

   &mpas_vars_nml
      mpas_state_variables = 'temperature',  'QTY_TEMPERATURE',
                             'salinity',     'QTY_SALINITY',
                             'rho',          'QTY_DENSITY',
                             'u',            'QTY_EDGE_NORMAL_SPEED',
                             'h',            'QTY_SEA_SURFACE_HEIGHT'
                             'tracer1',      'QTY_TRACER_CONCENTRATION'
      /

| 

Modules used
------------

::

   assim_model_mod.f90
   types_mod.f90
   location_mod.f90
   model_to_dart.f90
   model_mod.f90
   null_mpi_utilities_mod.f90
   obs_kind_mod.f90
   random_seq_mod.f90
   time_manager_mod.f90
   utilities_mod.f90

Files read
----------

-  MPAS analysis file; ``mpas_analysis.nc``
-  DART namelist file; ``input.nml``

Files written
-------------

-  DART initial conditions/restart file; e.g. ``dart_ics``

References
----------

none
