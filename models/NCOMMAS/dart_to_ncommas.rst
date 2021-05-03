PROGRAM ``dart_to_ncommas``
===========================

.. attention::

   ``NCOMMAS`` works with versions of DART *before* Manhattan (9.x.x) and has yet to be updated. If you are interested in
   using ``NCOMMAS`` with more recent versions of DART, contact DAReS staff to assess the feasibility of an update.
   Until that time, you should consider this documentation as out-of-date.



| ``dart_to_ncommas`` is the program that **updates** a ncommas netCDF-format restart file (usually
  ``ncommas_restart.nc``) with the state information contained in a DART output/restart file (e.g.
  ``perfect_ics, filter_ics, ...`` ). Only the CURRENT values in the ncommas restart file will be updated. The DART
  model time is compared to the time in the ncommas restart file. If the last time in the restart file does not match
  the DART model time, the program issues an error message and aborts.
| From the user perspective, most of the time ``dart_to_ncommas`` will be used on DART files that have a header
  containing one time stamp followed by the model state.
| The dart_to_ncommas_nml namelist allows ``dart_to_ncommas`` to read the ``assim_model_state_ic`` files that have *two*
  timestamps in the header. These files are temporarily generated when DART is used to advance the model. One timestamp
  is the 'advance_to' time, the other is the 'valid_time' of the model state. In this case, a namelist for ncommas
  (called ``ncommas_in.DART``) is written that contains the ``&time_manager_nml`` settings appropriate to advance
  ncommas to the time requested by DART. The repository version of the ``advance_model.csh`` script has a section to
  ensure the proper DART namelist settings for this case.
| Conditions required for successful execution of ``dart_to_ncommas``:

-  a valid ``input.nml`` namelist file for DART
-  a valid ``ncommas_vars.nml`` namelist file for ncommas - the same one used to create the DART state vector,
   naturally,
-  a DART file (typically ``filter_restart.xxxx`` or ``filter_ics.xxxx``)
-  a ncommas restart file (typically ``ncommas_restart.nc``).

Since this program is called repeatedly for every ensemble member, we have found it convenient to link the DART input
file to the default input filename (``dart_restart``). The same thing goes true for the ncommas output filename
``ncommas_restart.nc``.

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

   &dart_to_ncommas_nml
      dart_to_ncommas_input_file = 'dart_restart',
      advance_time_present   = .false.  
   /

| 

``dart_to_ncommas_nml`` and ``model_nml`` are always read from a file called ``input.nml``. The full description of the
``model_nml`` namelist is documented in the `NCOMMAS model_mod <model_mod.html#Namelist>`__.

.. container::

   +----------------------------+--------------------+------------------------------------------------------------------+
   | Item                       | Type               | Description                                                      |
   +============================+====================+==================================================================+
   | dart_to_ncommas_input_file | character(len=128) | The name of the DART file containing the model state to insert   |
   |                            |                    | into the ncommas restart file.                                   |
   +----------------------------+--------------------+------------------------------------------------------------------+
   | advance_time_present       | logical            | If you are converting a DART initial conditions or restart file  |
   |                            |                    | this should be ``.false.``; these files have a single timestamp  |
   |                            |                    | describing the valid time of the model state. If ``.true.`` TWO  |
   |                            |                    | timestamps are expected to be the DART file header. In this      |
   |                            |                    | case, a namelist for ncommas (called ``ncommas_in.DART``) is     |
   |                            |                    | created that contains the ``&time_manager_nml`` settings         |
   |                            |                    | appropriate to advance ncommas to the time requested by DART.    |
   +----------------------------+--------------------+------------------------------------------------------------------+

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

-  DART initial conditions/restart file; e.g. ``filter_ic``
-  DART namelist file; ``input.nml``
-  ncommas namelist file; ``ncommas_vars.nml``
-  ncommas restart file ``ncommas_restart.nc``

Files written
-------------

-  ncommas restart file; ``ncommas_restart.nc``
-  ncommas namelist file; ``ncommas_in.DART``

References
----------

none
