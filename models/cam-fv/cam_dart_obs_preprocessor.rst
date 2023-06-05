program ``cam_dart_obs_preprocessor``
=====================================

Overview
--------

This program provides a way to modifiy the incoming observation stream before running filter
for the CAM model.  Currently runs with the CAM-FV variant.  It uses the CAM model_mod namelist
to select the CAM level above which to remove observations.  The output obs_seq file contains
all observations below that level and none above.

This program could be extended to superob dense obs, adjust obs error values near the surface,
or any other desired changes to the incoming observation stream.

|

Usage
-----

Set the input and output obs_seq filenames in the namelist and run.

|

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::


   &cam_dart_obs_preprocessor_nml
      filename_in  = 'obs_seq.out'
      filename_out = 'obs_seq.nohighobs'
      verbose      = .true.
      calendar     = 'Gregorian'
      print_every  = 5000
      /

| 

Items in this namelist set the input and output files.

.. container::

   +--------------+--------------------+--------------------------------------------------------------------------------+
   | Item         | Type               | Description                                                                    |
   +==============+====================+================================================================================+
   | filename_in  | character(len=256) | Observation sequence file to read                                              |
   +--------------+--------------------+--------------------------------------------------------------------------------+
   | filename_out | character(len=256) | Observation sequence file to create and write. If this file exists it will be  |
   |              |                    | overwritten.                                                                   |
   +--------------+--------------------+--------------------------------------------------------------------------------+
   | verbose      | logical            | If true, write more debugging/status information.                              |
   |              |                    |                                                                                |
   +--------------+--------------------+--------------------------------------------------------------------------------+
   | calendar     | character(len=32)  | The string name of a valid DART calendar type. (See the                        |
   |              |                    | :doc:`../../../assimilation_code/modules/utilities/time_manager_mod` document  |
   |              |                    | for a list of valid types.) The setting here does not change what is written to|
   |              |                    | the output file; it only changes how the date information is printed to the    |
   |              |                    | screen in the informational messages.                                          |
   +--------------+--------------------+--------------------------------------------------------------------------------+
   | print_every  | integer            | If set to a positive integer, print a status message each time after           |
   |              |                    | processing that number of observations.                                        |
   +--------------+--------------------+--------------------------------------------------------------------------------+

| 

Discussion
----------

See the documentation in the obs_kind and obs_def modules for things you can query about an observation, and how to set
(overwrite) existing values.

| 

Building
--------

This program will be built when ``models/cam-fv/work/quickbuild.sh`` is run.

The ``preprocess`` program must be built and run first to define what set of observation types will be supported. See
the :doc:`../../../assimilation_code/programs/preprocess/preprocess` for more details on how to define the list and run
it. The ``&preprocess_nml`` namelist in the ``input.nml`` file must contain files with definitions for the combined set
of all observation types which will be encountered over all input obs_seq files.

If you have observation types which are not part of the default list in the &preprocess_nml namelist, add them to the
input.nml file and then run quickbuild.sh.


Files
-----

.. container::

   +---------------+-------------------------------------------------------------+
   | Filename      | Purpose                                                     |
   +===============+=============================================================+
   | input.nml     | to read the &cam_dart_obs_preprocess namelist               |
   +---------------+-------------------------------------------------------------+
   | caminput.nml  | (or the name in model_mod_nml:cam_template_filename)        |
   +---------------+-------------------------------------------------------------+
   |               | File required by static_init_model to run the program.      |
   +---------------+-------------------------------------------------------------+
   | cam_phis.nc   | (or the name inmodel_mod_nml:cam_phis_filename)             |
   +---------------+-------------------------------------------------------------+
   |               | File required by static_init_model to run the program.      |
   +---------------+-------------------------------------------------------------+

|

References
----------

#. none

Error codes and conditions
--------------------------

.. container:: errors

   +---------------------------+------------------------------------+-------------------------------------------------------+
   | Routine                   | Message                            | Comment                                               |
   +===========================+====================================+=======================================================+
   | cam_obs_dart_preprocess   | No obs in input sequence file      | Number of reported obs is 0a                          |
   +---------------------------+------------------------------------+-------------------------------------------------------+
   | cam_obs_dart_preprocess   | no first observation in {file}     | Unable to find the first observation in the filea     |
   +---------------------------+------------------------------------+-------------------------------------------------------+
   | cam_obs_dart_preprocess   | No obs will be written to {file}   | All obs were excluded; file_out will not be writtena  |
   +---------------------------+------------------------------------+-------------------------------------------------------+

|

Future plans
------------

none
