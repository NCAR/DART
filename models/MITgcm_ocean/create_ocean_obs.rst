PROGRAM ``create_ocean_obs``
============================


``create_ocean_obs`` is responsible for converting an interim ASCII file of ocean observations into a DART observation
sequence file. The interim ASCII file is a simple 'whitespace separated' table where each row is an observation and each
column is specific information about the observation.

+---------------+--------------------------+---------------------------------------------------------------------------------------------------------------------+
| column number | quantity                 | description                                                                                                         |
+===============+==========================+=====================================================================================================================+
| 1             | longitude (in degrees)   | longitude of the observation                                                                                        |
+---------------+--------------------------+---------------------------------------------------------------------------------------------------------------------+
| 2             | latitude (in degrees)    | latitude of the observation                                                                                         |
+---------------+--------------------------+---------------------------------------------------------------------------------------------------------------------+
| 3             | depth (in meters)        | depth of the observation                                                                                            |
+---------------+--------------------------+---------------------------------------------------------------------------------------------------------------------+
| 4             | observation value        | such as it is ...                                                                                                   |
+---------------+--------------------------+---------------------------------------------------------------------------------------------------------------------+
| 5             | vertical coordinate flag | There is a pathological difference between a surface observation and an observation with a depth of zero. See       |
|               |                          | `location_mod:location_type <../../assimilation_code/location/threed_sphere/location_mod.html#location_type>`__     |
|               |                          | for a full explanation. The short explanation is that *surface == -1*, and *depth == 3*                             |
+---------------+--------------------------+---------------------------------------------------------------------------------------------------------------------+
| 6             | observation variance     | good luck here ...                                                                                                  |
+---------------+--------------------------+---------------------------------------------------------------------------------------------------------------------+
| 7             | Quality Control flag     | integer value passed through to DART. There is a namelist parameter for                                             |
|               |                          | ``filter`` to ignore any observation with a QC value <=                                                             |
|               |                          | `input_qc_threshold <../../assimilation_code/programs/filter/filter.html#Namelist>`__                               |
+---------------+--------------------------+---------------------------------------------------------------------------------------------------------------------+
| 8             | obs_kind_name            | a character string that must match a string in                                                                      |
|               |                          | :doc:`../../observations/forward_operators/obs_def_ocean_mod`                                                       |
+---------------+--------------------------+---------------------------------------------------------------------------------------------------------------------+
| 9             | startDate_1              | the year-month-date of the observation (YYYYMMDD format)                                                            |
+---------------+--------------------------+---------------------------------------------------------------------------------------------------------------------+
| 10            | startDate_2              | the hour-minute-second of the observation (HHMMSS format)                                                           |
+---------------+--------------------------+---------------------------------------------------------------------------------------------------------------------+

For example:

::

   273.7500 21.3500 -2.5018 28.0441  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
   273.7500 21.4500 -2.5018 28.1524  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
   273.7500 21.5500 -2.5018 28.0808  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
   273.7500 21.6500 -2.5018 28.0143  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
   273.7500 21.7500 -2.5018 28.0242  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
   273.7500 21.8500 -2.5018 28.0160  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
   273.7500 21.9500 -2.5018 28.0077  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
   273.7500 22.0500 -2.5018 28.3399  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
   273.7500 22.1500 -2.5018 27.8852  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
   273.7500 22.2500 -2.5018 27.8145  3 0.0400  1  GLIDER_TEMPERATURE 19960101  10000
   ...

It is always possible to combine observation sequence files with the program
:doc:`../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool`, so it was simply convenient to generate a
separate file for each observation platform and type ('GLIDER' and 'TEMPERATURE'), however it is by no means required.

Modules used
------------

Some of these modules use modules ... **those** modules and namelists are not discussed here. probably should be ...

::

   types_mod
   utilities_mod
   dart_MITocean_mod
   obs_sequence_mod

Namelist
--------

This program has a namelist of its own, and some of the underlying modules require namelists. To avoid duplication and,
possibly, some inconsistency in the documentation; only a list of the required namelists is provided - with a hyperlink
to the full documentation for each namelist.

+----------------------------------------------------------+----------------------------------------------------------+
| Namelist                                                 | Primary Purpose                                          |
+==========================================================+==========================================================+
| `utilities_nml <../../assimilatio                        | set the termination level and file name for the run-time |
| n_code/modules/utilities/utilities_mod.html#Namelist>`__ | log                                                      |
+----------------------------------------------------------+----------------------------------------------------------+
| `obs_sequence_nml <../../assimilation_code               | write binary or ASCII observation sequence files         |
| /modules/observations/obs_sequence_mod.html#Namelist>`__ |                                                          |
+----------------------------------------------------------+----------------------------------------------------------+

We adhere to the F90 standard of starting a namelist with an ampersand '&' and terminating with a slash '/'. Consider
yourself forewarned that filenames that contain a '/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.

.. container:: namelist

   ::

      namelist /create_ocean_obs_nml/  year, month, day, &
               tot_days, max_num, fname, output_name, lon1, lon2, lat1, lat2

.. container:: indent1

   This namelist is read in a file called ``input.nml``

   +---------------------------------------+---------------------------------------+---------------------------------------+
   | Contents                              | Type                                  | Description                           |
   +=======================================+=======================================+=======================================+
   | year                                  | integer *[default: 1996]*             | The first year of interest.           |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | month                                 | integer *[default: 1]*                | The first month of interest.          |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | day                                   | integer *[default: 1]*                | The first day of interest.            |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | tot_days                              | integer *[default: 31]*               | Stop processing after this many days. |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | max_num                               | integer *[default: 800000]*           | The maximum number of observations to |
   |                                       |                                       | read/write.                           |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | fname                                 | character(len=129)                    | The name of the interim ASCII file of |
   |                                       | *[default: 'raw_ocean_obs.txt']*      | observations.                         |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | output_name                           | character(len=129)                    | The output file name.                 |
   |                                       | *[default: 'raw_ocean_obs_seq.out']*  |                                       |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | lon1                                  | real *[default: 0.0]*                 | The leftmost longitude of interest.   |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | lon2                                  | real *[default: 360.0]*               | The rightmost longitude of interest.  |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | lat1                                  | real *[default: -90.0]*               | The most southern latitude of         |
   |                                       |                                       | interest.                             |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | lat2                                  | real *[default: 90.0]*                | The most northern latitude of         |
   |                                       |                                       | interest.                             |
   +---------------------------------------+---------------------------------------+---------------------------------------+

| 

Files
-----

-  input namelist file: ``input.nml``
-  input data file: as listed by ``input.nml``\ ``&create_ocean_obs_nml:fname``
-  output data file: as listed by ``input.nml``\ ``&create_ocean_obs_nml:output_name``

References
----------

-  none
