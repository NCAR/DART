PROGRAM ``sst_to_obs, oi_sst_to_obs``
=====================================

Overview
--------

There are two gridded SST observation converters in this directory, one for data from PODAAC, and one from NOAA/NCDC.
``sst_to_obs`` converts data from PODAAC and has been used by Romain Escudier for regional studies with ROMS.
``oi_sst_to_obs`` converts data from NOAA/NCDC and has been used by Fred Castruccio for global studies with POP.

sst_to_obs -- GHRSST to DART observation sequence converter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These routines are designed to convert the `GHRSST Level 4 AVHRR_OI Global Blended Sea Surface Temperature Analysis (GDS
version 2) from NCEI data <https://podaac.jpl.nasa.gov/dataset/AVHRR_OI-NCEI-L4-GLOB-v2.0>`__ distributed by the
`Physical Oceanography Distributed Active Archive Center <http://podaac.jpl.nasa.gov>`__. Please remember to cite the
data in your publications, `specific instructions from PODAAC are available
here. <https://podaac.jpl.nasa.gov/dataset/AVHRR_OI-NCEI-L4-GLOB-v2.0>`__ This is an example:

   National Centers for Environmental Information. 2016. GHRSST Level 4 AVHRR_OI Global Blended Sea Surface Temperature
   Analysis (GDS version 2) from NCEI. Ver. 2.0. PO.DAAC, CA, USA. Dataset accessed [YYYY-MM-DD] at
   http://dx.doi.org/10.5067/GHAAO-4BC02.

**Many thanks to Romain Escudier (then at Rutgers) who did the bulk of the work and graciously contributed his efforts
to the DART project.** Romain gave us scripts and source code to download the data from the PODAAC site, subset the
global files to a region of interest, and convert that subsetted file to a DART observation sequence file. Those scripts
and programs have been only lightly modified to work with the Manhattan version of DART and contain a bit more
documentation.

The workflow is usually:

#. compile the converters by running ``work/quickbuild.sh`` in the usual way.
#. customize the ``shell_scripts/parameters_SST`` resource file to specify variables used by the rest of the scripting.
#. run ``shell_scripts/get_sst_ftp.sh`` to download the data from PODAAC.
#. provide a mask for the desired study area.
#. run ``shell_scripts/Prepare_SST.sh`` to subset the PODAAC data and create the DART observation sequence files. Be
   aware that the ``Prepare_SST.sh`` modifies the ``shell_scripts/input.nml.template`` file and generates its own
   ``input.nml``. ``work/input.nml`` is not used.
#. combine all output files for the region and timeframe of interest into one file using the
   :doc:`../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool`

Example
~~~~~~~

It is worth describing a small example. If you configure ``get_sst_ftp.sh`` to download the last two days of 2010 and
then specify the mask to subset for the NorthWestAtlantic (NWA) and run ``Prepare_SST.sh`` your directory structure
should look like the following:

::


   0[1234] cheyenne6:/<6>obs_converters/SST
   .
   |-- ObsData
   |   `-- SST
   |       |-- ncfile
   |       |   `-- 2010
   |       |       |-- 20101230120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0.nc
   |       |       `-- 20101231120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0.nc
   |       `-- nwaSST
   |           `-- 2010
   |               |-- 20101230120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0_NWA.nc
   |               `-- 20101231120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0_NWA.nc
   |-- oi_sst_to_obs.f90
   |-- oi_sst_to_obs.nml
   |-- sst_to_obs.f90
   |-- sst_to_obs.nml
   |-- shell_scripts
   |   |-- Prepare_SST.sh
   |   |-- functions.sh
   |   |-- get_sst_ftp.sh
   |   |-- input.nml
   |   |-- input.nml.template
   |   |-- my_log.txt
   |   |-- parameters_SST
   |   `-- prepare_SST_file_NWA.sh
   |-- masks
   |   |-- Mask_NWA-NCDC-L4LRblend-GLOB-v01-fv02_0-AVHRR_OI.nc
   |   `-- Mask_NWA120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0.nc
   `-- work
       |-- Makefile
       |-- advance_time
       |-- input.nml
       |-- obs_sequence_tool
       |-- oi_sst_to_obs
       |-- preprocess
       |-- quickbuild.sh
       `-- sst_to_obs

The location of the DART observation sequence files is specified by ``parameter_SST``:``DIR_OUT_DART``. That directory
should contain the following two files:

::

   0[1236] cheyenne6:/<6>v2/Err30 > ls -l
   'total 7104
   -rw-r--r-- 1 thoar p86850054 3626065 Jan 10 11:08 obs_seq.sst.20101230
   -rw-r--r-- 1 thoar p86850054 3626065 Jan 10 11:08 obs_seq.sst.20101231

oi_sst_to_obs -- noaa/ncdc to DART observation sequence converter
-----------------------------------------------------------------

``oi_sst_to_obs`` is designed to convert the `NOAA High-resolution Blended Analysis: Daily Values using AVHRR
only <https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html>`__ data. The global metadata of a
typical file is shown here:

::


   :Conventions = "CF-1.5" ;
   :title = "NOAA High-resolution Blended Analysis: Daily Values using AVHRR only" ;
   :institution = "NOAA/NCDC" ;
   :source = "NOAA/NCDC  ftp://eclipse.ncdc.noaa.gov/pub/OI-daily-v2/" ;
   :comment = "Reynolds, et al., 2007:
        Daily High-Resolution-Blended Analyses for Sea Surface Temperature.
        J. Climate, 20, 5473-5496.
        Climatology is based on 1971-2000 OI.v2 SST, 
        Satellite data: Navy NOAA17 NOAA18 AVHRR, Ice data: NCEP ice." ;
   :history = "Thu Aug 24 13:46:51 2017: ncatted -O -a References,global,d,, sst.day.mean.2004.v2.nc\n",
           "Version 1.0" ;
   :references = "https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html" ;
   :dataset_title = "NOAA Daily Optimum Interpolation Sea Surface Temperature" ;

The workflow is usually:

#. compile the converters by running ``work/quickbuild.sh`` in the usual way.
#. `download the desired data. <https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html>`__
#. customize the ``work/input.nml`` file.
#. run ``work/oi_sst_to_obs`` to create a single DART observation sequence file.
#. combine all output files for the region and timeframe of interest into one file using the
   :doc:`../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool`

sst_to_obs namelist
-------------------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &sst_to_obs_nml
      sst_netcdf_file     = '1234567.nc'
      sst_netcdf_filelist = 'sst_to_obs_filelist'
      sst_out_file        = 'obs_seq.sst'
      subsample_intv      = 1
      sst_rep_error       = 0.3
      debug               = .false.
      /

.. container::

   +----------------------+---------------------+----------------------------------------------------------------------------------------------------------+
   | Contents             | Type                | Description                                                                                              |
   +======================+=====================+==========================================================================================================+
   | sst_netcdf_file      | character(len=256)  | Name of the (usually subsetted) netcdf data file. This may be a relative or absolute filename.           |
   |                      |                     | If you run the scripts 'as is', this will be something like:                                             |
   |                      |                     | ``../ObsData/SST/nwaSST/2010/20101231120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0_NWA.nc``  |
   +----------------------+---------------------+----------------------------------------------------------------------------------------------------------+
   | sst_netcdf_filelist  | character(len=256)  | Name of the file that contains a list of (usually subsetted) data files, one per line.                   |
   |                      |                     | **You may not specify both sst_netcdf_file AND sst_netcdf_filelist.** One of them must be empty.         |
   +----------------------+---------------------+----------------------------------------------------------------------------------------------------------+
   | sst_out_file         | character(len=256)  | Name of the output observation sequence file.                                                            |
   +----------------------+---------------------+----------------------------------------------------------------------------------------------------------+
   | subsample_intv       | integer             | It is possible to 'thin' the observations. ``subsample_intv`` allows one to take every Nth observation.  |
   +----------------------+---------------------+----------------------------------------------------------------------------------------------------------+
   | sst_rep_error        | real                | In DART the observation error variance can be thought of as having two components, an instrument error   |
   |                      |                     | and a representativeness error. In ``sst_to_obs`` the instrument error is specified in the netCDF        |
   |                      |                     | file by the variable ``analysis_error``. The representativeness error is specified by                    |
   |                      |                     | ``sst_rep_error``, which is specified as a standard deviation.  These two values are added together      |
   |                      |                     | and squared and used as the observation error variance. **Note:**  This algorithm maintains backwards    |
   |                      |                     | compatibility, but is technically not the right way to combine these two quantities. If they both        |
   |                      |                     | specified variance, adding them together and then taking the square root would correctly specify         |
   |                      |                     | a standard deviation. Variances add, standard deviations do not. Since the true observation error        |
   |                      |                     | variance (in general) is not known, we are content to live with an algorithm that produces useful        |
   |                      |                     | observation error variances. If your research comes to a more definitive conclusion,                     |
   |                      |                     | please let us know.                                                                                      |
   +----------------------+---------------------+----------------------------------------------------------------------------------------------------------+
   | debug                | logical             | Print extra information during the ``sst_to_obs`` execution.                                             |
   +----------------------+---------------------+----------------------------------------------------------------------------------------------------------+

oi_sst_to_obs namelist
----------------------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &oi_sst_to_obs_nml
      input_file       = '1234567.nc'
      output_file_base = 'obs_seq.sst'
      subsample_intv   = 1
      sst_error_std    = 0.3
      debug            = .false.
      /

.. container::

   +-------------------+---------------------+----------------------------------------------------------------------------------------------------------+
   | Contents          | Type                | Description                                                                                              |
   +===================+=====================+==========================================================================================================+
   | input_file        | character(len=256)  | Name of the input netcdf data file. This may be a relative or absolute                                   |
   |                   |                     | filename. If you run the scripts 'as is', this will be something like:                                   |
   |                   |                     | ``../ObsData/SST/nwaSST/2010/20101231120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0_NWA.nc``  |
   +-------------------+---------------------+----------------------------------------------------------------------------------------------------------+
   | output_file_base  | character(len=256)  | Partial filename for the output file.  The date and time are appended to ``output_file_base``            |
   |                   |                     | to construct a unique filename reflecting the time of the observations in the file.                      |
   +-------------------+---------------------+----------------------------------------------------------------------------------------------------------+
   | subsample_intv    | integer             | It is possible to 'thin' the observations. ``subsample_intv``                                            |
   |                   |                     | allows one to take every Nth observation.                                                                |
   +-------------------+---------------------+----------------------------------------------------------------------------------------------------------+
   | sst_error_std     | real                | This is the total observation error standard deviation.                                                  |
   +-------------------+---------------------+----------------------------------------------------------------------------------------------------------+
   | debug             | logical             | Print extra information during the ``oi_sst_to_obs`` execution.                                          |
   +-------------------+---------------------+----------------------------------------------------------------------------------------------------------+

Decisions you might need to make
--------------------------------

See the general discussion in the :doc:`../../../guide/creating-obs-seq-real` page about what options are
available for the things you need to specify. These include setting a time, specifying an expected error, setting a
location, and an observation type.


Known Bugs
----------

I do not believe ``sst_to_obs`` will work correctly
if given multiple files in ``sst_netcdf_filelist``.
The number of observation used to declare the length of the output 
observation sequence is based on a single file ... yet seems to be used 
by many. I have not tested this configuration, since the scripting does 
not use the ``sst_netcdf_filelist`` mechanism.
