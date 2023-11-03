PROGRAM ``snow_to_obs``
=======================

MODIS snowcover fraction observation converter
----------------------------------------------

Overview
^^^^^^^^

There are several satellite sources for snow observations. Generally the data is distributed in HDF-EOS format. The
converter code in this directory DOES NOT READ HDF FILES as input. It expects the files to have been preprocessed to
contain text, one line per observation, with northern hemisphere data only.

Data sources
------------

not sure.

Programs
--------

The ``snow_to_obs.f90`` file is the source for the main converter program.

To compile and test, go into the work subdirectory and run the ``quickbuild.sh`` script to build the converter and a
couple of general purpose utilities. ``advance_time`` helps with calendar and time computations, and the
``obs_sequence_tool`` manipulates DART observation files once they have been created.

This converter creates observations of the "MODIS_SNOWCOVER_FRAC" type.

There is another program in this directory called ``snow_to_obs_netcdf.f90`` which is a prototype for reading netcdf
files that contain some metadata and presumably have been converted from the original HDF. THIS HAS NOT BEEN TESTED but
if you have such data, please contact dart@ucar.edu for more assistance. If you write something that reads the HDF-EOS
MODIS files directly, please, please contact us! Thanks.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &snow_to_obs_nml
     longrid         = 360,
     latgrid         = 90, 
     year            = 2000, 
     doy             = 1,
     snow_input_file = 'snowdata.input', 
     missing_value   = -20.0, 
     debug           = .false.
   /

+-----------------+--------------------+-----------------------------------------------------------------------------+
| Item            | Type               | Description                                                                 |
+=================+====================+=============================================================================+
| longrid         | integer            | The number of divisions in the longitude dimension.                         |
+-----------------+--------------------+-----------------------------------------------------------------------------+
| latgrid         | integer            | The number of divisions in the latitude dimension. This converter assumes   |
|                 |                    | the data is for the northern hemisphere only. A namelist item could be      |
|                 |                    | added to select northern verses southern hemisphere if needed.              |
+-----------------+--------------------+-----------------------------------------------------------------------------+
| year            | integer            | The year number of the data.                                                |
+-----------------+--------------------+-----------------------------------------------------------------------------+
| doy             | integer            | The day number in the year. Valid range 1 to 365 in a non-leap year, 1 to   |
|                 |                    | 366 in a leap year.                                                         |
+-----------------+--------------------+-----------------------------------------------------------------------------+
| snow_input_file | character(len=128) | The name of the input file.                                                 |
+-----------------+--------------------+-----------------------------------------------------------------------------+
| missing_value   | real(r8)           | The value used to mark missing data.                                        |
+-----------------+--------------------+-----------------------------------------------------------------------------+
| debug           | logical            | If set to .true. the converter will print out more information as it does   |
|                 |                    | the conversion.                                                             |
+-----------------+--------------------+-----------------------------------------------------------------------------+


Known Bugs
^^^^^^^^^^
This program is hardcoded to read only northern hemisphere data. It should handle global values.


Future Plans
^^^^^^^^^^^^
This program should use the HDF-EOS libraries to read the native MODIS granule files. Right now the ascii intermediate files contain no metadata, so if the namelist values don't match the actual division of the globe, bad things will happen.



