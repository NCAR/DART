PROGRAM ``closest_member_tool``
===============================

Overview
--------

Utility program to compare the ensemble mean to an ensemble of restart files, which can now be run in parallel. The
program prints out a sorted order of which members are 'closest' to the mean, where the method used to determine 'close'
is selectable by namelist option. It also creates a file with a single number or character string in it, for ease in
scripting, which identifies the closest member.

The ensemble mean is computed from the input ensemble. The difference is computed point by point across the ensemble
members. There is an option to restrict the computation to just a subset of the entire state vector by listing one or
more generic quantities. In this case, only state vector items matching one of these quantities will contribute to the
total difference value.

Available methods are:

1 - simple absolute difference:
   The absolute value of the difference between each item in the mean vector and the corresponding item in each ensemble
   member, accumulated over the entire state vector.
2 - normalized absolute difference:
   The absolute value of the difference between each item in the mean vector and the corresponding item in each ensemble
   member normalized by the mean value, accumulated over the entire state vector.
3 - simple RMS difference:
   The square root of the accumulated sum of the square of the difference between each item in the mean vector and the
   corresponding item in each ensemble member.
4 - normalized RMS difference:
   The square root of the accumulated sum of the square of the normalized difference between each item in the mean
   vector and the corresponding item in each ensemble member.

This program could be used to select one or more ensemble members to run a free model forecast forward in time after the
assimilation is finished. Each member is an equally likely representation of the model state. Using the ensemble mean
may not be the best choice since the mean may not have self-consistent fine-scale structures in the data.

In addition to printing out data about all members to both the console and to the dart log file, this program creates a
single output file containing information about the closest member. If the input restart data is in a single file, the
output file 'closest_restart' contains a single number which is the ensemble member number. If the input restart data is
in separate files, the output file contains the full filename of the closest member, e.g. 'filter_restart.0004' if
member 4 is closest. For scripting the contents of this file can be used to copy the corresponding member data and
convert it to the model input format for a free forecast, for example.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &closest_member_tool_nml
      input_restart_files     = ''
      input_restart_file_list = ''      
      output_file_name        = 'closest_restart'
      ens_size                = 1
      difference_method       = 4      
      use_only_qtys           = ''
      single_restart_file_in  = .false.      
     /

.. container::

   +---------------------------------------+---------------------------------------+---------------------------------------+
   | Item                                  | Type                                  | Description                           |
   +=======================================+=======================================+=======================================+
   | input_restart_files                   | character(len=256),                   | An array of filenames each containing |
   |                                       | dimension(ens_size x num_domains)     | a list DART restart data.             |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | input_restart_file_list               | character(len=256),                   | A file containing a list of filenames |
   |                                       | dimension(num_domains)                | for DART restart data, one for each   |
   |                                       |                                       | domain.                               |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | output_file_name                      | character(len=256)                    | This is a file containing the member  |
   |                                       |                                       | number that is closest to the         |
   |                                       |                                       | ensemble mean.                        |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | ens_size                              | integer                               | Total number of ensemble members.     |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | difference_method                     | integer                               | Select which method is used to        |
   |                                       |                                       | compute 'distance' from mean:         |
   |                                       |                                       |                                       |
   |                                       |                                       | -  1 = simple absolute difference     |
   |                                       |                                       | -  2 = absolute difference normalized |
   |                                       |                                       |    by the mean                        |
   |                                       |                                       | -  3 = simple RMS difference          |
   |                                       |                                       | -  4 = RMS of the normalized          |
   |                                       |                                       |    difference                         |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | use_only_quantities                   | character(len=32)                     | If unspecified, all items in the      |
   |                                       |                                       | state vector contribute to the total  |
   |                                       |                                       | difference. If one or more quantities |
   |                                       |                                       | are listed here, only items in the    |
   |                                       |                                       | state vector of these quantities      |
   |                                       |                                       | contribute to the total difference.   |
   |                                       |                                       | These are the generic quantities,     |
   |                                       |                                       | such as QTY_TEMPERATURE,              |
   |                                       |                                       | QTY_U_WIND_COMPONENT, QTY_DENSITY,    |
   |                                       |                                       | etc. and not specific types like      |
   |                                       |                                       | RADIOSONDE_TEMPERATURE. Consult the   |
   |                                       |                                       | model interface code to determine     |
   |                                       |                                       | which possible quantities are         |
   |                                       |                                       | returned by the                       |
   |                                       |                                       | `get_state_me                         |
   |                                       |                                       | ta_data() <../../../models/template/m |
   |                                       |                                       | odel_mod.html#get_state_meta_data>`__ |
   |                                       |                                       | routine.                              |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | single_restart_file_in                | logical                               | **Not supported yet.** Contact        |
   |                                       |                                       | dart@ucar.edu if you are interested   |
   |                                       |                                       | in using this tool with files that    |
   |                                       |                                       | contain all ensemble members in a     |
   |                                       |                                       | single file.                          |
   +---------------------------------------+---------------------------------------+---------------------------------------+

Below is an example of a typical namelist for the closest_member_tool.

::

   &closest_member_tool_nml
      input_restart_files     = ''
      input_restart_file_list = 'restart_list.txt'      
      output_file_name        = 'closest_restart.txt'
      ens_size                = 3
      single_restart_file_in  = .false.      
      difference_method       = 4      
      use_only_qtys           = ''
     /

where ``restart_list.txt`` contains

::

   cam_restart_0001.nc
   cam_restart_0002.nc
   cam_restart_0003.nc

Currently ``single_restart_file_in`` is not supported. This is typically used for simpler models that have built in
model advances such as ``lorenz_96``.

| 

Files
-----

-  inputfile.####.nc (list of restarts to find closest member) -or-
-  ``restart_list.txt`` (a file containing a list of restart files) and,
-  ``input.nml``

References
----------

-  none
