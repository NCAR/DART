PROGRAM ``ncommas_to_dart``
===========================

.. attention::

   ``NCOMMAS`` works with versions of DART *before* Manhattan (9.x.x) and has yet to be updated. If you are interested in
   using ``NCOMMAS`` with more recent versions of DART, contact DAReS staff to assess the feasibility of an update.
   Until that time, you should consider this documentation as out-of-date.



| ``ncommas_to_dart`` is the program that reads a ncommas restart file (usually ``ncommas_restart.nc``) and creates a
  DART state vector file (e.g. ``perfect_ics, filter_ics, ...`` ).
| The list of variables used to create the DART state vector are specified in the ``ncommas_vars.nml`` file.
| Conditions required for successful execution of ``ncommas_to_dart``:

-  a valid ``input.nml`` namelist file for DART
-  a valid ``ncommas_vars.nml`` namelist file for ncommas
-  the ncommas restart file mentioned in the ``input.nml&model_nml:ncommas_restart_filename`` variable.

Since this program is called repeatedly for every ensemble member, we have found it convenient to link the ncommas
restart files to the default input filename (``ncommas_restart.nc``). The default DART state vector filename is
``dart_ics`` - this may be moved or linked as necessary.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &model_nml
      ncommas_restart_filename     = 'ncommas_restart.nc';
      assimilation_period_days     = 1,
      assimilation_period_seconds  = 0,
      model_perturbation_amplitude = 0.2,
      output_state_vector          = .true.,
      calendar                     = 'Gregorian',
      debug                        = 0
   /

| 

::

   &ncommas_to_dart_nml
      ncommas_to_dart_output_file = 'dart_ics'  
   /

| 

``ncommas_to_dart_nml`` and ``model_nml`` are always read from a file called ``input.nml``. The full description of the
``model_nml`` namelist is documented in the `NCOMMAS model_mod <model_mod.html#Namelist>`__.

.. container::

   +-----------------------------+--------------------+-----------------------------------------------------------------+
   | Item                        | Type               | Description                                                     |
   +=============================+====================+=================================================================+
   | ncommas_to_dart_output_file | character(len=128) | The name of the DART file which contains the updated model      |
   |                             |                    | state info that should be written into the NCOMMAS file.        |
   +-----------------------------+--------------------+-----------------------------------------------------------------+

| 

``ncommas_vars_nml`` is always read from a file called ``ncommas_vars.nml``.

.. container::

   +---------------------------------------+---------------------------------------+---------------------------------------+
   | Item                                  | Type                                  | Description                           |
   +=======================================+=======================================+=======================================+
   | ncommas_state_variables               | character(len=NF90_MAX_NAME) ::       | The list of variable names in the     |
   |                                       | dimension(160)                        | NCOMMAS restart file to use to create |
   |                                       |                                       | the DART state vector and their       |
   |                                       |                                       | corresponding DART kind.              |
   +---------------------------------------+---------------------------------------+---------------------------------------+

| 

::

   &ncommas_vars_nml
      ncommas_state_variables = 'U',   'QTY_U_WIND_COMPONENT',
                                'V',   'QTY_V_WIND_COMPONENT',
                                'W',   'QTY_VERTICAL_VELOCITY',
                                'TH',  'QTY_POTENTIAL_TEMPERATURE',
                                'DBZ', 'QTY_RADAR_REFLECTIVITY',
                                'WZ',  'QTY_VERTICAL_VORTICITY',
                                'PI',  'QTY_EXNER_FUNCTION',
                                'QV',  'QTY_VAPOR_MIXING_RATIO',
                                'QC',  'QTY_CLOUDWATER_MIXING_RATIO',
                                'QR',  'QTY_RAINWATER_MIXING_RATIO',
                                'QI',  'QTY_ICE_MIXING_RATIO',
                                'QS',  'QTY_SNOW_MIXING_RATIO',
                                'QH',  'QTY_GRAUPEL_MIXING_RATIO'
     /

| 

Modules used
------------

::

   assim_model_mod
   location_mod
   model_mod
   null_mpi_utilities_mod
   obs_kind_mod
   random_seq_mod
   time_manager_mod
   types_mod
   utilities_mod

Files read
----------

-  ncommas restart file; ``ncommas_restart.nc``
-  DART namelist files; ``input.nml`` and ``ncommas_vars.nml``

Files written
-------------

-  DART state vector file; e.g. ``dart_ics``

References
----------

none
