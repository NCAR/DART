PROGRAM ``prepbufr``
====================

General Overview
----------------

Converting NCEP PREPBUFR files into DART obs_seq.out files is a 2 step process. The first step
uses NCEP software to convert the PREPBUFR file into an intermediate text file. This is described in this document.
The second step converts the intermediate files into obs_seq.out files, which is done by create_real_obs, as
described in :doc:`../ascii_to_obs/create_real_obs` .  For details of the conversion process we refer users to the
the sections starting with **Prepbufr Overview** (on this page below) and the create_real_obs documentation page.  
For those with access to Derecho looking for a quick overview of all steps in the conversion process see the **Quickstart
Instructions** below.

Quickstart Instructions
-----------------------

This section provides a complete list of steps within the conversion process (prepbufr file to obs_seq.out) performed on Derecho. 
It is only recommended if you have a Derecho account, or prefer a comprehensive list of the required steps/commands. 
Otherwise please proceed to the **Prepbufr Overview** section for more details.


- The following Derecho modules were loaded during the compiling process. This exact setup
  is not mandatory, but given for reference.

   ::

     Currently Loaded Modules:
        1) ncarenv/23.09 (S)   3) intel/2023.2.1        5) cray-mpich/8.1.27  7) netcdf/4.9.2   9) ncl/6.6.2 
        2) craype/2.7.23       4) ncarcompilers/1.0.0   6) hdf5/1.12.2        8) nco/5.2.4 

-  Go to the DART prep_bufr observation converter directory and
   build the PREPBUFR utilities as follows:

   ::

      cd $DART_DIR/observations/obs_converters/NCEP/prep_bufr
      ./install.sh

   The ``install.sh`` script can be edited to use different compilers by setting the CCOMP and FCOMP 
   sections at the start of the script.

   Confirm the exe directory contains the executables ``prepbufr.x``, ``prepbufr_03Z.x``,
   ``grabbufr.x`` and ``cword.x``.

-  Go to the ``$DART_DIR/build_templates/`` directory and use the intel build template
   for compiling DART

   ::

      cd $DART_DIR/build_templates
      cp mkmf.template.intel.linux mkmf.template

-  Go to the ``$DART_DIR/observations/obs_converters/NCEP/prep_bufr/work/``
   directory and run ``quickbuild.sh`` to build the ``advance_time`` executable:

   ::

      cd $DART_DIR/observations/obs_converters/NCEP/prep_bufr/work
      ./quickbuild.sh


-  Make a ``prepqm`` subdirectory in the data directory:

   ::

      cd $DART_DIR/observations/obs_converters/NCEP/prep_bufr/data
      mkdir prepqm
 

-  Obtain the PREPBUFR observations for your desired time. In the following examples
   observations from April 27, 2017 will be used.  On Derecho you can go to the
   the campaign collections directory at:

   ::

      cd /glade/campaign/collections/rda/data/d090000/2017/
      cp A25328-201704prepqmB.tar $DART_DIR/observations/obs_converters/NCEP/prep_bufr/data/prepqm/


   Alternatively, download the files directly from the `NSF NCAR Research Data
   Archive <NCEP+NCAR_obs_>`_ page for the
   NCEP/NSF NCAR Global Reanalysis Products. Register on the site, click on
   the "Data Access" tab, and locate the "prepqm: BUFR observations files".
   Use a similar approach to  obtain the `RDA ds337 <NCEP_obs_>`_
   prepbufr files. Locate and download the GDAS PREPBUFR files of your choice.   
  

-  You may not be able to untar the downloaded files. In that case
   cosblocking command will unblock the file in place as:

   ::

       
      cd $DART_DIR/observations/obs_converters/NCEP/prep_bufr/data/prepqm/
      /glade/campaign/cisl/dares/Observations/cosconvert -b A25328-201704prepqmB.tar


-  Untar the data (organized into half months), then unzip the days of interest (6 hour chunks): 
   To process a full day of observations you will need both the day of interest, and the following day.

   ::

      tar -xvf A25328-201704prepqmB.tar
      cd 201704
      gzip -d prepqm17042*nr.gz


   Some tar files create a subdirectory when the archive is unpacked.  It may be easier to
   move all the files from the subdirectory to the main prepqm directory:

   ::
     
     cd $DART_DIR/observations/obs_converters/NCEP/prep_bufr/data/prepqm/201704
     mv prepqm* ..

   Now all the individual 6-hour prepqm* files are in the data/prepqm directory.
   (If you decide to leave the files in the 201704 subdirectory, edit the BUFR_in
   variable in the prepbufr.csh script to include the subdirectory name.)


- Confirm these files exist:  ``prepqm17042706.nr``, ``prepqm17042712.nr``, ``prepqm17042718.nr``,
  and ``prepqm17042800.nr`` and ``prepqm17042806.nr``.



-  In the ``$DART_DIR/observations/obs_converters/NCEP/prep_bufr/work``
   directory, edit the ``input.nml`` file. This file will control what
   observations will be used for your experiment, so the namelist
   options are worth investigating a bit here. 

   For example, for weather applications 
   where you want only atmospheric observations close to your assimilation time, you could
   use the following:

   ::

      &prep_bufr_nml
         obs_window    = 1.0
         obs_window_cw = 1.5
         otype_use     = 120.0, 130.0, 131.0, 132.0, 133.0, 180.0
                         181.0, 182.0, 220.0, 221.0, 230.0, 231.0
                         232.0, 233.0, 242.0, 243.0, 245.0, 246.0
                         252.0, 253.0, 255.0, 280.0, 281.0, 282.0
         qctype_use    = 0,1,2,3,15
         /

   This defines an observation time window (``obs_window``) of +/- 1.0 hours, while cloud
   motion vectors (``obs_window_cw``) will be used over a window of +/- 1.5 hours.  Observations
   outside of this time window are excluded from the analysis. 

   This will use observation types:

       - sounding temps (120),
       - aircraft temps (130,131),
       - dropsonde temps (132),
       - mdcars aircraft temps (133),
       - marine temp (180),
       - land humidity (181),
       - ship humidity (182),
       - rawinsonde U,V (220),
       - pibal U,V (221),
       - aircraft U,V (230,231,232),
       - cloudsat winds (242,243,245),
       - GOES water vapor (246),
       - sat winds (252,253,255), and
       - ship obs (280, 281,282)
 
   Additionally, it will include observations with specified qc
   types only. Skip to the prepfbufr **Namelist** section at the bottom
   of this page for more available namelist controls.

   For applications where you want to use all available observations 
   in the 6-hour assimilation window, you would set obs_window 
   to 3.0 or larger to avoid excluding any observations.

-  Within the ``$DART_DIR/observations/obs_converters/NCEP/prep_bufr/work``
   directory, edit the ``prepbufr.csh`` file and change the following
   variables to match the locations and format of the data you downloaded
   as shown below:
  
   ::

      set BUFR_dir  = ../data
      set BUFR_idir = ${BUFR_dir}/prepqm
      set BUFR_odir = ${BUFR_dir}/prepout
      set DART_exec_dir = ../exe
      ..
      set BUFR_in = ${BUFR_idir}/prepqm${sdtg}.nr
      
      

-  The ``daily`` variable in ``prepbufr.csh`` controls whether the
   output of the script is a single 24 hour file or 4, 6-hour files.  Depending 
   on what other observations you intend to merge later on,
   one option may be easier to use.  For this example, we are
   leaving ``daily = no``, and the script will make 4, 6-hour files.

-  Run the ``prepbufr.csh`` script for a single day:

   ::

      cd $DART_DIR/observations/obs_converters/NCEP/prep_bufr/work
      ./prepbufr.csh 2017 04 27



-  Your PREPBUFR files have now been converted to an intermediate ASCII
   format. Confirm that ``temp_obs.20170427*`` files within 
   ``~/data/prepout`` exist. Please note that the script can function
   with only the *06* prepqm input file, but will also need the
   *12*, *18*, and following day *00* files to run to completion. 


-  There is a separate executable to convert the
   observations from this intermediate text format into
   the native DART format. 
   Edit the ``input.nml`` namelist file in the
   ``DART_DIR/observations/obs_converters/NCEP/ascii_to_obs/work``.
   For this example:

   ::

      &ncepobs_nml
         year       = 2017,
         month      = 4,
         day        = 27,
         tot_days   = 1,
         max_num    = 800000,
         select_obs = 0,
         ObsBase = '../../prepbufr/data/prepout/temp_obs.',
         daily_file = .false.,
         lat1       = 15.0,
         lat2       = 60.0,
         lon1       = 200.0,
         lon2       = 330.0
         /

   Setting ``select_obs = 0`` will select all the observations in the
   ASCII file. Set ``ObsBase`` to the directory you output the intermediate files from
   during the last step. If you wish to choose specific observations
   from the ASCII intermediate file or control other program behavior,
   there are many namelist options documented on the
   :doc:`create_real_obs <../../../../observations/obs_converters/NCEP/ascii_to_obs/create_real_obs>`
   page.

-  It is now time to build *ascii_to_obs*. Run the following:

   ::

      cd $DART_DIR/observations/obs_converters/NCEP/ascii_to_obs/work
      ./quickbuild.sh

-  Run ``create_real_obs`` to create the DART observation
   sequence files:


   ::

      cd $DART_DIR/observations/obs_converters/NCEP/ascii_to_obs/work
      ./create_real_obs

-  The executable ``create_real_obs`` will create observation sequence files
   with one file for each six hour window. Confirm that the ``obs_seq20170427*``
   files have been generated.  For a cycled experiment, there
   are several options for naming the output observation files.  Each month
   could be in a separate directory, or each 6-hour file could be in its
   own directory.

   For a single obs file per directory, you can rename each as *obs_seq.out*,
   or you can include the timestamp in the filename, e.g. *obs_seq.2017042706.out*.

-  The observation types within the file should look like:

   ::

      obs_sequence
      obs_type_definitions
          19
           1 RADIOSONDE_U_WIND_COMPONENT
           2 RADIOSONDE_V_WIND_COMPONENT
           5 RADIOSONDE_TEMPERATURE
           6 RADIOSONDE_SPECIFIC_HUMIDITY
          12 AIRCRAFT_U_WIND_COMPONENT
          13 AIRCRAFT_V_WIND_COMPONENT
          14 AIRCRAFT_TEMPERATURE
          16 ACARS_U_WIND_COMPONENT
          17 ACARS_V_WIND_COMPONENT
          18 ACARS_TEMPERATURE
          20 MARINE_SFC_U_WIND_COMPONENT
          21 MARINE_SFC_V_WIND_COMPONENT
          22 MARINE_SFC_TEMPERATURE
          23 MARINE_SFC_SPECIFIC_HUMIDITY
          30 SAT_U_WIND_COMPONENT
          31 SAT_V_WIND_COMPONENT
          40 RADIOSONDE_SURFACE_ALTIMETER
          42 MARINE_SFC_ALTIMETER
          43 LAND_SFC_ALTIMETER
     num_copies:            1  num_qc:            1
     num_obs:        68107  max_num_obs:        68107
 
-  Some models include a preprocessing program to do additional processing
   of the observations, for example, limiting obs to a sub-domain, superobbing spatially dense obs,
   and increasing the error for near boundary observations.

   For example, the WRF weather model includes the
   :doc:`wrf_dart_obs_preprocess <../../../../models/wrf/WRF_DART_utilities/wrf_dart_obs_preprocess>`
   program.  There may be no further processing needed for some models, and the observation sequence
   file is ready to be used.


**You have completed the Quickstart Instructions**. See the following sections for more details of the 
prepbufr conversion package.



Prepbufr Overview
-----------------

The prep_bufr package is external NCEP code and has not been completely incorporated into the DART architecture. It
requires adaptation of the source codes and scripts to the computing environment where it will be run. It is not so
robust that it can be controlled just with input parameters. It may not have the same levels of error detection and
warning that the rest of DART has, so the user should very careful about checking the end product for correctness.


Install Prepbufr package
^^^^^^^^^^^^^^^^^^^^^^^^

Running the ``install.sh`` script located within the ``$DART_DIR$/observations/obs_converters/NCEP/prep_bufr`` directory will build the library
and main executable. You will probably have to edit this script to set the fortran compiler on your system.

If you have raw unblocked PREPBUFR files you will need to convert them to blocked format (what prepbufr expects as
input). The blk/ublk section of the build script compiles the ``cword.x`` converter program.

If you are running on an Intel (little-endian) based machine you will need the ``grabbufr`` byte swapping program that
is also built by this script.

One-shot mode
^^^^^^^^^^^^^

If you are converting a single obs file, or are walking through the process by hand for the first time, you can follow
the more detailed build instructions below, and then run the prep_bufr.x program by hand. This involves the following
steps:

-  Build the executables.
-  Run the blocker if needed (generally not if you have downloaded the blocked format PREPBUFR files).
-  Run the binary format converter if you are on an Intel (little-endian) machine.
-  Link the input file to a fixed input filename
-  Run prepbufr.x to convert the file
-  Copy the fixed output filename to the desired output filename

Production mode
^^^^^^^^^^^^^^^

If you have multiple days (or months) of observations to convert, there is a script in the work
subdirectory which is set up to run the converter on a sequence of raw data files, and concatenate the output files
together into one output file per day. Edit the ``work/prepbufr.csh`` script (as described in the Quickstart section) 
and set the necessary values in the 'USER SET PARAMETERS' section near the top. This script can either be run from 
the command line, or it can be submitted to a batch queue for a long series of conversion runs.

Overview of Prepbufr package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This package is currently organized into files under the ``DART/observations/NCEP/prep_bufr`` directory:

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

Decoding program: src/prepbufr.f
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The program ``prepbufr.f`` is used to decode the NCEP reanalysis PREPBUFR data into intermediate text files. This program
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

To compile the BUFR libraries and the decoding program, set the CCOMP and FCOMP variables in the install.sh script 
to match the compilers available on your system.  Execute the install.sh script to complete the
compilations for the main decoding program, the NCEP BUFR library, and the conversion utilities.

The executables (i.e., prepbufr.x, prepbufr_03Z.x) are placed in the ../exe directory.

Platforms tested:

-  Linux clusters with Intel, PGI, Pathscale, GNU Fortran,
-  Mac OS X with Intel, GNU Fortran,
-  SGI Altix with Intel
-  Cray with Intel, Cray Fortran.

Byte-swapping program: convert_bufr/grabbufr.f
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For platforms with little-endian binary file format (e.g. Intel, AMD®, and non-MIPS SGI processors) the program
``grabbufr.f`` is used to convert the big-endian format NCEP PREPBUFR data into little-endian format. The ``grabbufr.f`` code is
written in Fortran 90, and has been compiled can be compiled with the pgf90 compiler on a Linux system, with gfortran on
an Intel based Mac, and the ifort compiler on other Linux machines. The ``install.sh`` script should build this by default, 
however instructions are in ``convert_bufr/README``.  In case of problems, go to the ``convert_bufr`` subdirectory, 
edit ``convert_bufr.csh`` to set your compiler, and run it to compile the converter code (grabbufr).

This program reads the PREPBUFR file into memory, and needs to know the size of the file (in bytes).
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

Blocking program blk_ublk/cword.x
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``prepbufr.x`` program expects to read a blocked input file, which is generally what is available for download. However,
if you have an unblocked file that you need to convert, there is a conversion program. The ``install.sh`` script will try to
build this by default, but in case of problems you can build it separately. Change directories into the ``blk_ublk``
subdirectory and read the ``README_cwordsh`` file for more help. The cwordsh shell-script wrapper shows how to run the
executable ``cwordsh.x`` executable.

This program is not required for blocked file formats.

Downloading Prepbufr raw data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The NCEP PREPBUFR files (prepqmYYMMDDHH) can be found within the NCEP/NCAR Global Reanalysis Products dataset, d090000, 
on NSF NCAR Research Data Archive (RDA).  Operational observations can be found in the NCEP ADP Global Upper Air and 
Surface Weather Observations dataset, d337000.

To find the files:

-  go to the `NSF NCAR/NCEP reanalysis archive. <NCEP+NCAR_obs_>`_
-  Click on the "Data Access" tab.
-  Locate the **preqm: BUFR observation files**
-  Click on Complete File List link and Select the year you are interested in.
-  Depending on the year the format of the filenames change, but they should contain the year, usually as 2 digits, the
   month, and then either the start/stop day for weekly files, or the letters A and B for semi-monthly files.

Depending on the year you select, the prepqm files can be weekly, monthly, or semi-monthly. Each tar file has a unique
dataset number of the form "A#####". For example, for January of 2003, the 4 HPSS TAR files are: A21899, A21900, A21901,
A21902. After September 2003, these files include AIRCRAFT data (airplane readings taken at cruising elevation) but not
ACARS data (airplane readings taken during takeoff and landing). There are different datasets which include ACARS data
but their use is restricted and you must contact the RDA group to get access.

| If you are running on a machine with direct access to the NSF NCAR HPSS, then change directories into the prep_bufr/data
  subdirectory and obtain the prepqm rawfile from:
| *> cd /glade/campaign/collections/rda/data/d#####*
| where ##### is the data set number you want.

| These files may be readable tar files, or they may require running the ``cosconvert`` program first. See if the
  ``tar`` command can read them:
| *> tar -tvf rawfile*
| If you get a good table of contents then simply rename the file and untar it:
| *> mv rawfile data.tar*
| *> tar -xvf data.tar*
| However, if you get an error from the tar command, on the NSF NCAR machine Derecho run:
| *> /glade/campaign/cisl/dares/Observations/cosconvert -b data.tar*

The output of tar should yield individual 6-hourly NCEP PREPBUFR data files for the observations in the +/- 3-hour time
windows of 00Z, 06Z, 12Z, and 18Z of each day. Note that DART obs_seq files are organized such that a 24 hour file with
4 6-hour observation windows would contain observations from 3:01Z to 3:00Z of the next day, centered on 6Z, 12Z, 18Z and "24Z".
In addition, there are some observations at 3:00Z on the PREPBUFR file labelled with 06Z. Then, in order to make a full
day intermediate file incorporating all the required obs from the "next" day, you'll need the PREPBUFR files through 6Z
of the day after the last day of interest. For example, to generate the observation sequence for Jan 1, 2003, the
decoded NCEP PREPBUFR text files for Jan 1 and 2, 2003 are needed, and hence the following PREPBUFR files are needed:

-  prepqm03010106
-  prepqm03010112
-  prepqm03010118
-  prepqm03010200
-  prepqm03010206


Execution of Prepbufr
~~~~~~~~~~~~~~~~~~~~~

In ``prep_bufr/work/prepbufr.csh`` set the appropriate values of the year, month, first day, and last day of the period you
desire, and the variable "convert" to control conversion from big- to little-endian. Confirm that the raw PREPBUFR files
are in ../data, or that prepbufr.csh has been changed to find them. Execute ``prepbufr.csh`` in the work directory.

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
      obs_window       = 3.0,
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
