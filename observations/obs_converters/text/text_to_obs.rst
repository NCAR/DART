PROGRAM ``text_to_obs``
=======================

Text file to DART converter
---------------------------

Overview
~~~~~~~~

If you have observations in spreadsheet or column format, in text, with a single line per observation, then the files
this directory are a template for how to convert these observations into a format suitable for DART use.

The workflow is usually:

-  read in the needed information about each observation - location, time, data value, observation type - from a data
   source (usually a file)
-  call a series of DART library routines to construct a derived type that contains all the information about a single
   observation
-  call another set of DART library routines to put it into a time-sorted series
-  repeat the last 2 steps until all observations are processed
-  finally, call a write subroutine that writes out the entire series to a file in a format that DART can read in

It is not recommended that you try to mimic the ascii file format by other means; the format is subject to change and
the library routines will continue to be supported even if the physical format changes.

If your input data is in some kind of format like netCDF or HDF, then one of the other converters (e.g. the MADIS ones
for netCDF) might be a better starting place for adapting code.

Data sources
------------

This part is up to you. For each observation you will need a location, a data value, a type, a time, and some kind of
error estimate. The error estimate can be hardcoded in the converter if they are not available in the input data. See
below for more details on selecting an appropriate error value.

Programs
--------

The ``text_to_obs.f90`` file is the source for the main converter program. Look at the source code where it reads the
example data file. You will almost certainly need to change the "read" statement to match your data format. The example
code reads each text line into a character buffer and then reads from that buffer to parse up the data items.

To compile and test, go into the work subdirectory and run the ``quickbuild.sh`` script to build the converter and a
couple of general purpose utilities. ``advance_time`` helps with calendar and time computations, and the
``obs_sequence_tool`` manipulates DART observation files once they have been created.

To change the observation types, look in the ``DART/observations/forward_operators`` directory. If you can find an "obs_def_XXX_mod.f90" file
with an appropriate set of observation types, change the 'use' lines in the converter source to include those types.
Then add that filename in the ``input.nml`` namelist file to the &preprocess_nml namelist, the 'input_files' variable.
Multiple files can be listed. Then run quickbuild.sh again. It remakes the table of supported observation types before
trying to recompile the source code.

An example script for converting batches of files is in the ``shell_scripts`` directory. A tiny example data file is in
the ``data`` directory. These are *NOT* intended to be turnkey scripts; they will certainly need to be customized for
your use. There are comments at the top of the script saying what options they include, and should be commented enough
to indicate where changes will be likely to need to be made.

Decisions you might need to make
--------------------------------

See the discussion in the :doc:`../../../guide/creating-obs-seq-real` page about what options are available
for the things you need to specify. These include setting a time, specifying an expected error, setting a location, and
an observation type.
