PROGRAM ``prepbufr``
====================

Overview
--------

Translating NCEP PREPBUFR files into DART obs_seq.out files (input file to filter) is a 2 stage process. The first stage
uses NCEP software to translate the PREPBUFR file into an intermediate text file. This is described in this document.
The second step is to translate the intermediate files into obs_seq.out files, which is done by create_real_obs, as
described in :doc:`../ascii_to_obs/create_real_obs` .

Instructions
------------

The prep_bufr package is free-standing and has not been completely assimilated into the DART architecture. It also
requires adaptation of the sources codes and scripts to the computing environment where it will be run. It is not so
robust that it can be controlled just with input parameters. It may not have the same levels of error detection and
warning that the rest of DART has, so the user should very careful about checking the end product for correctness.

Overview of what needs to be built and run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

More detailed instructions follow, but this section describes a quick overview of what programs you will be building and
running.

Building
^^^^^^^^

Running the install.sh script will build the library and main executable. You will probably have to edit this script to
set which fortran compiler is available on your system.

If you have raw unblocked PREPBUFR files you will need to convert them to blocked format (what prepbufr expects as
input). The blk/ublk section of the build script compiles the ``cword.x`` converter program.

If you are running on an Intel (little-endian) based machine you will need the ``grabbufr`` byte swapping program that
is also built by this script.

One-shot execution
^^^^^^^^^^^^^^^^^^

If you are converting a single obs file, or are walking through the process by hand for the first time, you can follow
the more detailed build instructions below, and then run the prep_bufr.x program by hand. This involves the following
steps:

-  building the executables.
-  running the blocker if needed (generally not if you have downloaded the blocked format PREPBUFR files).
-  running the binary format converter if you are on an Intel (little-endian) machine.
-  linking the input file to a fixed input filename
-  running prepbufr.x to convert the file
-  copying the fixed output filename to the desired output filename

Production mode
^^^^^^^^^^^^^^^

If you have multiple days (or months) of observations that you are intending to convert, there is a script in the work
subdirectory which is set up to run the converter on a sequence of raw data files, and concatenate the output files
together into one output file per day. Edit the work/prepbufr.csh script and set the necessary values in the 'USER SET
PARAMETERS' section near the top. This script can either be run from the command line, or it can be submitted to a batch
queue for a long series of conversion runs.

Installation of the ncep prepbufr decoding program
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This package is currently organized into files under the DART/observations/NCEP/prep_bufr directory:

::

   src           Source code of the NCEP PREPBUFR decoder
   lib           NCEP BUFR library source
   install.sh    A script to install the NCEP PREPBUFR decoder and the NCEP BUFR library.
   exe           Executables of the decoder and converter.
   data          Where the NCEP PREPBUFR files (prepqm****) could be loaded into
                 from the NSF NCAR Mass Store (the script assumes this is the default location).
   work          Where we run the script to do the decoding.
   convert_bufr  Source code (grabbufr) to convert the binary big-endian PREPBUFR files to 
                 little-endian files, and a script to compile the program.
   blk_ublk      Source code (cwordsh) to convert between blocked and unblocked format.
   docs          Some background information about NCEP PREPBUFR observations.

The decoding program: src/prepbufr.f
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The program prepbufr.f is used to decode the NCEP reanalysis PREPBUFR data into intermediate text files. This program
was originally developed by NCEP. It has been modified to output surface pressure, dry temperature, specific humidity,
and wind components (U/V) of conventional radiosonde, aircraft reports, and satellite cloud motion derived wind. There
are additional observation types on the PREPBUFR files, but using them they would require significant modifications of
prepbufr and require detailed knowledge of the NCEP PREPBUFR files. The NCEP quality control indexes for these
observations based on NCEP forecasts are also output and used in DART observation sequence files. The NCEP PREPBUFR
decoding program is written in Fortran 77 and has been successfully compiled on Linux computers using pgi90, SGI®
computers with f77, IBM® SP® systems with xlf, and Intel® based Mac® with gfortran.

If your operating system uses modules you may need to remove the default compiler and add the one desired for this
package. For example

-  which pgf90 (to see if pgf90 is available.)
-  module rm intel64 netcdf64 mpich64
-  module add pgi32

To compile the BUFR libraries and the decoding program, set the CPLAT variable in the install.sh script to match the
compilers available on your system. CPLAT = linux is the default. Execute the install.sh script to complete the
compilations for the main decoding program, the NCEP BUFR library, and the conversion utilities.

The executables (i.e., prepbufr.x, prepbufr_03Z.x) are placed in the ../exe directory.

Platforms tested:

-  Linux clusters with Intel, PGI, Pathscale, GNU Fortran,
-  Mac OS X with Intel, GNU Fortran,
-  SGI Altix with Intel
-  Cray with Intel, Cray Fortran.

The byte-swapping program convert_bufr/grabbufr.f
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For platforms with little-endian binary file format (e.g. Intel, AMD®, and non-MIPS SGI processors) the program
grabbufr.f is used to convert the big-endian format NCEP PREPBUFR data into little-endian format. The grabbufr.f code is
written in Fortran 90, and has been compiled can be compiled with the pgf90 compiler on a Linux system, with gfortran on
an Intel based Mac, and the ifort compiler on other Linux machines. More detailed instructions for building it can be
found in convert_bufr/README, but the base install script should build this by default. In case of problems, cd into the
convert_bufr subdirectory, edit convert_bufr.csh to set your compiler, and run it to compile the converter code
(grabbufr).

This program reads the whole PREPBUFR file into memory, and needs to know the size of the file (in bytes).
Unfortunately, the system call STAT() returns this size as one number in an array, and the index into that array differs
depending on the system and sometimes the word size (32 vs 64) of the compiler. To test that the program is using the
right offset into this array, you can compile and run the stat_test.f program. It takes a single filename argument and
prints out information about that file. One of the numbers will be the file size in bytes. Compare this to the size you
see with the 'ls -l' command for that same file. If the numbers do not agree, find the right index and edit the
grabbufr.f source file. Look for the INDEXVAL line near the first section of executable code.

If grabbufr.f does not compile because the getarg() or iargc() subroutines are not found or not available, then either
use the arg_test.f program to debug how to get command line arguments into a fortran program on your system, or simply
go into the grabbufr.f source and comment out the section which tries to parse command line arguments and comment in the
hardcoded input and output filenames. Now to run this program you must either rename the data files to these
predetermined filenames, or you can use links to temporarily give the files the names needed.

The blocking program blk_ublk/cword.x
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The prepbufr.x program expects to read a blocked input file, which is generally what is available for download. However,
if you have an unblocked file that you need to convert, there is a conversion program. The install.sh script will try to
build this by default, but in case of problems you can build it separately. Change directories into the blk_ublk
subdirectory and read the README_cwordsh file for more help. The cwordsh shell-script wrapper shows how to run the
executable cwordsh.x executable.

Note that if you can get the blocked file formats to begin with, this program is not needed.

Getting the ncep reanalysis prepbufr format data from NSF NCAR HPSS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The NCEP PREPBUFR files (prepqmYYMMDDHH) can be found within the NCEP reanalysis dataset, ds090.0, on NSF NCAR Mass Store
System (HPSS).

To find the files:

-  go to the `NSF NCAR/NCEP reanalysis archive. <http://rda.ucar.edu/datasets/ds090.0/>`__
-  Click on the "Inventories" tab.
-  Select the year you are interested in.
-  Search for files with the string "prepqm" in the name.
-  Depending on the year the format of the filenames change, but they should contain the year, usually as 2 digits, the
   month, and then either the start/stop day for weekly files, or the letters A and B for semi-monthly files.

Depending on the year you select, the prepqm files can be weekly, monthly, or semi-monthly. Each tar file has a unique
dataset number of the form "A#####". For example, for January of 2003, the 4 HPSS TAR files are: A21899, A21900, A21901,
A21902. After September 2003, these files include AIRCRAFT data (airplane readings taken at cruising elevation) but not
ACARS data (airplane readings taken during takeoff and landing). There are different datasets which include ACARS data
but their use is restricted and you must contact the RDA group to get access.

| If you are running on a machine with direct access to the NSF NCAR HPSS, then change directories into the prep_bufr/data
  subdirectory and run:
| *> hsi get /DSS/A##### rawfile*
| where ##### is the data set number you want.

| These files may be readable tar files, or they may require running the ``cosconvert`` program first. See if the
  ``tar`` command can read them:
| *> tar -tvf rawfile*
| If you get a good table of contents then simply rename the file and untar it:
| *> mv rawfile data.tar*
| *> tar -xvf data.tar*
| However, if you get an error from the tar command you will need to run the ``cosconvert`` program to convert the file
  into a readable tar file. On the NSF NCAR machine *yellowstone*, run:
| *> /glade/u/home/rdadata/bin/cosconvert -b rawfile data.tar*
| On other platforms, download the appropriate version from: http://rda.ucar.edu/libraries/io/cos_blocking/utils/ .
  Build and run the converter and then you should have a tar file you can unpack.

The output of tar should yield individual 6-hourly NCEP PREPBUFR data files for the observations in the +/- 3-hour time
windows of 00Z, 06Z, 12Z, and 18Z of each day. Note that DART obs_seq files are organized such that a 24 hour file with
4 observation times would contain observations from 3:01Z to 3:00Z of the next day, centered on 6Z, 12Z, 18Z and "24Z".
In addition, there are some observations at 3:00Z on the PREPBUFR file labelled with 06Z. Then, in order to make a full
day intermediate file incorporating all the required obs from the "next" day, you'll need the PREPBUFR files through 6Z
of the day after the last day of interest. For example, to generate the observation sequence for Jan 1, 2003, the
decoded NCEP PREPBUFR text files for Jan 1 and 2, 2003 are needed, and hence the PREPBUFR files

-  prepqm03010106
-  prepqm03010112
-  prepqm03010118
-  prepqm03010200
-  prepqm03010206

are needed.

Running the ncep prepbufr decoding program
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In prep_bufr/work/prepbufr.csh set the appropriate values of the year, month, first day, and last day of the period you
desire, and the variable "convert" to control conversion from big- to little-endian. Confirm that the raw PREPBUFR files
are in ../data, or that prepbufr.csh has been changed to find them. Execute prepbufr.csh in the work directory. It has
code for running in the LSF batch environment, but not PBS.

Currently, this script generates decoded PREPBUFR text data each 24 hours which contains the observations within the
time window of -3:01 hours to +3:00Z within each six-hour synoptic time. These daily output text files are named as
temp_obs.yyyymmdd. These text PREPBUFR data files can then be read by
DART/observations/NCEP/ascii_to_obs/work/:doc:`../ascii_to_obs/create_real_obs` to generate the DART daily observation
sequence files.

There is an alternate section in the script which creates a decoded PREPBUFR text data file each 6 hours (so they are
1-for-1 with the original PREPBUFR files). Edit the script prepbufr.csh and look for the commented out code which
outputs 4 individual files per day. Note that if you chose this option, you will have to make corresponding changes in
the create_obs_seq.csh script in step 2.

Other modules used
------------------

This is a piece of code that is intended to be 'close' to the original, as such, we have not modified it to use the DART
build mechanism. This code does not use any DART modules.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &prep_bufr_nml
      obs_window       = 1.5,
      obs_window_upa   = 1.5,
      obs_window_air   = 1.5,
      obs_window_sfc   = 0.8,
      obs_window_cw    = 1.5,
      land_temp_error  = 2.5,
      land_wind_error  = 3.5,
      land_moist_error = 0.2,
      otype_use        = missing,
      qctype_use       = missing,
   /

| 

.. container::

   +---------------------+--------------+-------------------------------------------------------------------------------+
   | Item                | Type         | Description                                                                   |
   +=====================+==============+===============================================================================+
   | obs_window          | real         | Window of time to include observations. If > 0, overrides all the other more  |
   |                     |              | specific window sizes. Set to -1.0 to use different time windows for          |
   |                     |              | different obs types. The window is +/- this number of hours, so the total     |
   |                     |              | window size is twice this value.                                              |
   +---------------------+--------------+-------------------------------------------------------------------------------+
   | obs_window_upa      | real         | Window of time to include sonde observations (+/- hours) if obs_window is <   |
   |                     |              | 0, otherwise ignored.                                                         |
   +---------------------+--------------+-------------------------------------------------------------------------------+
   | obs_window_air      | real         | Window of time to include aircraft observations (+/- hours) if obs_window is  |
   |                     |              | < 0, otherwise ignored.                                                       |
   +---------------------+--------------+-------------------------------------------------------------------------------+
   | obs_window_sfc      | real         | Window of time to include surface observations (+/- hours) if obs_window is < |
   |                     |              | 0, otherwise ignored.                                                         |
   +---------------------+--------------+-------------------------------------------------------------------------------+
   | obs_window_cw       | real         | Window of time to include cloud wind observations (+/- hours) if obs_window   |
   |                     |              | is < 0, otherwise ignored.                                                    |
   +---------------------+--------------+-------------------------------------------------------------------------------+
   | otype_use           | real(300)    | Report Types to extract from bufr file. If unspecified, all types will be     |
   |                     |              | converted.                                                                    |
   +---------------------+--------------+-------------------------------------------------------------------------------+
   | qctype_use          | integer(300) | QC types to include from the bufr file. If unspecified, all QC values will be |
   |                     |              | accepted.                                                                     |
   +---------------------+--------------+-------------------------------------------------------------------------------+
   | land_temp_error     | real         | observation error for land surface temperature observations when none is in   |
   |                     |              | the input file.                                                               |
   +---------------------+--------------+-------------------------------------------------------------------------------+
   | land_wind_error     | real         | observation error for land surface wind observations when none is in the      |
   |                     |              | input file.                                                                   |
   +---------------------+--------------+-------------------------------------------------------------------------------+
   | land_moisture_error | real         | observation error for land surface moisture observations when none is in the  |
   |                     |              | input file.                                                                   |
   +---------------------+--------------+-------------------------------------------------------------------------------+

| 

Files
-----

-  input file(s); NCEP PREPBUFR observation files named using ObsBase with the "yymmddhh" date tag on the end. Input to
   grabbufr if big- to little-endian is to be done. Input to prepbufr if not.
-  intermediate (binary) prepqm.little; output from grabbufr, input to prepbufr.
-  intermediate (text) file(s) "temp_obs.yyyymmddhh"; output from prepbufr, input to create_real_obs

References
----------

DART/observations/NCEP/prep_bufr/docs/\* (NCEP text files describing the PREPBUFR files)
