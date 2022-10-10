PROGRAM ``tc_to_obs``
=====================

Tropical Cyclone ATCF File to DART Converter
============================================

Overview
--------

Tropical Cyclone data created by the 'Automated Tropical Cyclone Forecast (ATCF) System' can be converted into DART
observations of the storm center location, minimum sea level pressure, and maximum wind speed. Several of the options
can be customized at runtime by setting values in a Fortran namelist. See the namelist section below for more details.
In the current release of DART only the :doc:`../../../models/wrf/readme` has forward operator code to generate
expected obs values for these vortex observations.

`This webpage <http://www.ral.ucar.edu/hurricanes/realtime/index.php#about_atcf_data_files>`__ documents many things
about the ATCF system and the various file formats that are used for storm track data and other characteristics.

The converter in this directory is only configured to read the packed "b-deck" format (as described on the webpage
referenced above). There are sections in the fortran code which can be filled in to read other format variants. This
should mostly be a matter of changing the read format string to match the data in the file.

| 

Data sources
------------

A collection of past storm ATCF information can be found `here <http://www.ral.ucar.edu/hurricanes/repository>`__. For
each observation you will need a location, a data value, a type, a time, and some kind of error estimate. The error
estimates will need to be hardcoded or computed in the converter since they are not available in the input data. See
below for more details on selecting an appropriate error value.

| 

Programs
--------

The ``tc_to_obs.f90`` file is the source for the main converter program. Look at the source code where it reads the
example data file. Given the variety of formatting details in different files, you may quite possibly need to change the
"read" statement to match your data format. There is a 'select case' section which is intended to let you add more
formats and select them at runtime via namelist.

To compile and test, go into the work subdirectory and run the ``quickbuild.sh`` script to build the converter and a
couple of general purpose utilities. ``advance_time`` helps with calendar and time computations, and the
``obs_sequence_tool`` manipulates DART observation files once they have been created.

This converter creates observation types defined in the ``DART/observations/forward_operators/obs_def_vortex_mod.f90``
file. This file must be listed in the ``input.nml`` namelist file, in the ``&preprocess_nml`` namelist, in the
'input_files' variable, for any programs which are going to process these observations. If you have to change the
``&preprocess_nml`` namelist you will have to run ``quickbuild.sh`` again to build and execute the ``preprocess``
program before compiling other executables. It remakes the table of supported observation types before trying to
recompile other source code.

There is an example b-deck data file in the ``data`` directory. This format is what is supported in the code as
distributed. There are other variants of this format which have more spaces so the columns line up, and variants which
have many more fields than what is read here.

| 

Specifying expected error
-------------------------

The ATCF files DO NOT include any estimated error values. The source code currently has hardcoded values for location,
sea level pressure, and max wind errors. These may need to be adjusted as needed if they do not give the expected
results.

| 

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &tc_to_obs_nml
      input_atcf_file         = 'input.txt'
      fileformat              = 'b-deck'
      obs_out_file            = 'obs_seq.out'
      append_to_existing_file = .false.
      debug                   = .false.
    /

| 

.. container::

   +--------------------------+---------------------+---------------------------------------------------------------------------------+
   | Item                     | Type                | Description                                                                     |
   +==========================+=====================+=================================================================================+
   | input_atcf_file          | character(len=256)  | Name of the input ascii text file in ATCF format.                               |
   +--------------------------+---------------------+---------------------------------------------------------------------------------+
   | fileformat               | character(len=128)  | Currently only supports 'b-deck' but if other format strings                    |
   |                          |                     | are added, can switch at runtime between reading                                |
   |                          |                     | different varieties of ATCF file formats.                                       |
   +--------------------------+---------------------+---------------------------------------------------------------------------------+
   | obs_out_file             | character(len=256)  | Name of the output observation sequence file to create.                         |
   +--------------------------+---------------------+---------------------------------------------------------------------------------+
   | append_to_existing_file  | logical             | If .false., this program will overwrite an existing file. If .true.             |
   |                          |                     | and if a file already exists with the same name the newly converted             |
   |                          |                     | observations will be appended to that file. Useful if you have multiple         |
   |                          |                     | small input files that you want to concatenate into a single output             |
   |                          |                     | file. However, there is no code to check for duplicated observations. If        |
   |                          |                     | this is .true. and you run the converter twice you will get duplicate           |
   |                          |                     | observations in the file which is bad. (It will affect the quality of           |
   |                          |                     | your assimilation results.) Use with care.  You can concatenate multiple        |
   |                          |                     | obs sequence files as a postprocessing step with the                            |
   |                          |                     | :doc:`../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool`  |
   |                          |                     | which comes with DART and is built by the quickbuild.sh script in               |
   |                          |                     | the TC converter work directory.                                                |
   +--------------------------+---------------------+---------------------------------------------------------------------------------+
   | debug                    | logical             | Set to .true. to print out more details during the conversion process.          |
   +--------------------------+---------------------+---------------------------------------------------------------------------------+

| 
