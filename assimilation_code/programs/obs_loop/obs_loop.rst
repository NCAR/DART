program ``obs_loop``
====================

Overview
--------

This program is a template that is intended to be modified by the user to do any desired operations on an observation
sequence file.

Usage
-----

This program is intended to be used as a template to read in observations from one obs_seq file and write them,
optionally modified in some way, to another obs_seq file. It can be compiled and run as-is, but it simply makes an exact
copy of the input file.

There are comments in the code (search for ``MODIFY HERE`` ) where you can test values, types, times, error values, and
either modify them or skip copying that observation to the output.

There are build files in ``observations/utilities/oned`` and ``observations/utilities/threed_sphere`` to build the
``obs_loop`` program.

| 

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &obs_loop_nml
      filename_in  = ''
      filename_out = '' 
      print_only   = .false.
      calendar     = 'Gregorian'
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
   | print_only   | logical            | If .TRUE. then do the work but only print out information about what would be  |
   |              |                    | written as output without actually creating the output file.                   |
   +--------------+--------------------+--------------------------------------------------------------------------------+
   | calendar     | character(len=32)  | The string name of a valid DART calendar type. (See the                        |
   |              |                    | :doc:`../../modules/utilities/time_manager_mod` documentation for a list of    |
   |              |                    | valid types.) The setting here does not change what is written to the output   |
   |              |                    | file; it only changes how the date information is printed to the screen in the |
   |              |                    | informational messages.                                                        |
   +--------------+--------------------+--------------------------------------------------------------------------------+

| 

Discussion
----------

See the documentation in the obs_kind and obs_def modules for things you can query about an observation, and how to set
(overwrite) existing values.

| 

Building
--------

There are build files in ``observations/utilities/oned`` and ``observations/utilities/threed_sphere`` to build the
``obs_loop`` program.

The ``preprocess`` program must be built and run first to define what set of observation types will be supported. See
the :doc:`../../../assimilation_code/programs/preprocess/preprocess` for more details on how to define the list and run
it. The ``&preprocess_nml`` namelist in the ``input.nml`` file must contain files with definitions for the combined set
of all observation types which will be encountered over all input obs_seq files.

If you have observation types which are not part of the default list in the &preprocess_nml namelist, add them to the
input.nml file and then run quickbuild.sh.


Files
-----

========= ==================================
filename  purpose
========= ==================================
input.nml to read the &obs_loop_nml namelist
========= ==================================

References
----------

#. none

Error codes and conditions
--------------------------

.. container:: errors

   ======== ======= =======
   Routine  Message Comment
   ======== ======= =======
   obs_loop         
   obs_loop         
   ======== ======= =======

Future plans
------------

none
