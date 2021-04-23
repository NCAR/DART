program ``obs_keep_a_few``
==========================

Overview
--------

This program creates an output observation sequence (obs_seq) file that is shorter than the input obs_seq file. There
are two ways to restrict the number of observations copied to the output: the total number of observations regardless of
observation type, or up to N observations of each type. Observations in an obs_seq file are processed in time order so
the observations with the earliest timestamps will be copied.

Set either limit to -1 to disable it. If both the maximum count per type and maximum total count are given the copying
stops when the first limit is reached.

If you want to subset an obs_seq file starting at a later time see the :doc:`../obs_sequence_tool/obs_sequence_tool` for
subsetting by time and then use this tool on the output. That tool also allows you to subset by obs type, location, data
value, and a variety of other options.

The ``obs_keep_a_few`` program only subsets by numbers of observations. It is expected to be useful when prototyping
experiments so the run time is short, or for debugging or testing. Setting a limit per type ensures you have up to N of
each type of observation present in the output file.

Identity observations are all considered to be the same identity "observation type" by this tool.

Other modules used
------------------

::

   types_mod
   utilities_mod
   location_mod
   obs_def_mod
   obs_kind_mod
   time_manager_mod
   obs_sequence_mod

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &obs_keep_a_few_nml
      filename_in        = ''
      filename_out       = '' 
      max_count_per_type = 10
      max_total_count    = -1
      print_only         = .false.
      calendar           = 'Gregorian'
      /

| 

.. container::

   +--------------------+--------------------+--------------------------------------------------------------------------+
   | Item               | Type               | Description                                                              |
   +====================+====================+==========================================================================+
   | filename_in        | character(len=256) | Name of the observation sequence file to read.                           |
   +--------------------+--------------------+--------------------------------------------------------------------------+
   | filename_out       | character(len=256) | Name of the observation sequence file to create. An existing file will   |
   |                    |                    | be overwritten.                                                          |
   +--------------------+--------------------+--------------------------------------------------------------------------+
   | max_count_per_type | integer            | The first N observations of each different type will be copied to the    |
   |                    |                    | output file. Observation sequence files are processed in time order so   |
   |                    |                    | these will be the ones with the earliest time stamps relative to other   |
   |                    |                    | observations of this same type. Set to -1 to disable this limit.         |
   +--------------------+--------------------+--------------------------------------------------------------------------+
   | max_total_count    | integer            | If greater than 0, sets the upper limit on the total number of           |
   |                    |                    | observations to be copied to the output file regardless of type. The     |
   |                    |                    | program quits when either this limit is reached or when there are N of   |
   |                    |                    | each different obs type in the output. Set to -1 to disable.             |
   +--------------------+--------------------+--------------------------------------------------------------------------+
   | print_only         | logical            | If true, does all the work and prints out what the output file would     |
   |                    |                    | have in it (timestamps and counts of each obs type) but doesn't create   |
   |                    |                    | the output file.                                                         |
   +--------------------+--------------------+--------------------------------------------------------------------------+
   | calendar           | character(len=256) | Name of the DART calendar type to use. Generally 'Gregorian' or 'No      |
   |                    |                    | calendar'. See the DART time manager for more options. Only controls the |
   |                    |                    | formatting of how the times in the output summary messages are           |
   |                    |                    | displayed.                                                               |
   +--------------------+--------------------+--------------------------------------------------------------------------+

| 

Files
-----

-  filename_in is read.
-  filename_out is written.

References
----------

-  none
