PROGRAMS ``LPRM_L3_to_obs.f90`` ``AMSR_E_L2_to_obs.f90``
===========================================================================
This is a brief description of the converters and utilities in this directory
retrieved from the `NASA Earthdata portal <https://earthdata.nasa.gov/>`__
This is a front end for many (Distributed Active Archive Center) DAAC portals
and for the `Goddard Earth Sciences Data and Information Services Center <https://disc.gsfc.nasa.gov>`__

These directories contain satellite retrieval data for land surface soil moisture
and leaf area index (LAI). 


The general workflow for each of the observation converters
(described in more detail below) is usually:


   1. Download the data for the period in question
      (see DATA SOURCES below)
   2. Build the DART executables with support for the soil moisture observations.
      This is done by running preprocess with
      ``obs_def_land_mod.f90`` in the list of input_files
      for ``preprocess_nml``.
   3. Convert each data file individually (e.g. executing ``LPRM_L3_to_obs``)
   4. Combine or subset all output files for the region and timeframe of interest
      into one file using :doc:`../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool`

For some models (CLM, for example), it is required to reorganize the observation sequence
files into a series of files that contains ONLY the observations for each assimilation.
This can be achieved with the `~/models/clm/shell_scripts/makedaily.sh` script. Since
there are subtleties for each model, makedaily.sh is generally found in the shell_scripts
directory of the model.
 

Soil Moisture Observation Converters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Program  ``AMSR_E_L2_to_obs.f90``

Description
-----------

AMSR-E/Aqua surface soil moisture (LPRM) L2B V002 is a Level 2 (swath) data set. 
Its land surface parameters, surface soil moisture, land surface (skin) temperature, 
and vegetation water content are derived from passive microwave remote sensing data from 
the Advanced Microwave Scanning Radiometer-Earth Observing System (AMSR-E), using the 
Land Parameter Retrieval Model (LPRM). Each swath is packaged with associated geolocation fields. 
The data set covers the period from June 2002 to October 2011 (when the AMSR-E on the NASA EOS 
Aqua satellite stopped producing data due to a problem with the rotation of its antenna), at the 
spatial resolution (nominally 56 and 38 km, respectively) of AMSR-E's C and X bands 
(6.9 and 10.7 GHz, respectively).

NAMELIST
--------

This namelist is read from the file input.nml.
Namelists start with an ampersand
'&amp;' and terminate with a slash '/'.
Character strings that contain a '/' must be
enclosed in quotes to prevent them from
prematurely terminating the namelist.

::

  &AMSRE_L2_to_obs_nml
     input_file          = 'LPRM-AMSR_E_L2_D_SOILM2_V002_20030630025503.nc'
     obs_out_file        = 'obs_seq.out'
     max_rfi_code        = 2
     amsre_rep_error     = 2
     /

Description of namelist variables:

+--------------------+--------------------+---------------------------------------------------------------------------+
| Contents           | Type               | Description                                                               |
+====================+====================+===========================================================================+
| input_file         | character(len=256) | Name of the netcdf soil moisture data file.                               |
+--------------------+--------------------+---------------------------------------------------------------------------+
| obs_out_file       | character(len=256) | Name of the output observation sequence file.                             |
+--------------------+--------------------+---------------------------------------------------------------------------+
| max_rfi_code       | integer            | Maximum Radio Frequency Interference. Soil moisture values with a         |
|                    |                    | max rfi above this threshold are excluded from obs_out_file               |
+--------------------+--------------------+---------------------------------------------------------------------------+
| amsre_rep_error    | integer            | Representativeness Error (standard deviation units). This value is added  |
|                    |                    | to the instrument error provided from the input file to calculate the     |
|                    |                    | total observation error variance                                          |
+--------------------+--------------------+---------------------------------------------------------------------------+

.. Important::

  The total error (instrument error + representativeness error) is used to calculate the observation error
  variance in the ``obs_out_file``. The instrument error is taken from the data file variable (sm_c_error
  or sm_x_error), whereas the represenativeness error (amsre_rep_error) is a namelist variable input by the
  user at runtime.



Data Source
-----------

The dataset (LPRM_AMSRE_SOILM2: AMSR-E/Aqua surface soil moisture (LPRM) L2B V002) can be
found `here. <https://disc.gsfc.nasa.gov/datasets/LPRM_AMSRE_SOILM2_002/summary>`__



Program ``LPRM_L3_to_obs.f90``

Description
-----------

TMI/TRMM surface soil moisture (LPRM) L3 1 day 25 km x 25 km nighttime V001
is Level 3 (gridded) data set. Its land surface parameters, surface soil moisture,
land surface (skin) temperature, and vegetation water content, are derived from
passive microwave remote sensing data from the Tropical Rainfall Measuring Mission (TRMM)
Microwave Imager (TMI), using the Land Parameter Retrieval Model (LPRM). There are
two files per day, one daytime and one nighttime, archived as two different products.
This document is for the nighttime product. The data set covers the period from
December 1997 to April 2015 (when the instruments on the TRMM satellite were shut
down in preparation for its reentry into the earth's atmosphere).

The LPRM is based on a forward radiative transfer model to retrieve surface
soil moisture and vegetation optical depth. The land surface temperature is
derived separately from TMI's Ka-band (37 GHz). A unique feature of this method
is that it can be applied at any microwave frequency, making it very suitable to
exploit all the available passive microwave data from various satellites.


NAMELIST
--------

This namelist is read from the file input.nml.
Namelists start with an ampersand
'&amp;' and terminate with a slash '/'.
Character strings that contain a '/' must be
enclosed in quotes to prevent them from
prematurely terminating the namelist.



::

  &LPRM_L3_to_obs_nml
     input_file        = 'LPRM-TMI_L3_NT_SOILM3_V001-20120411T144345Z_20120407.nc'
     output_file       = 'obs_seq.out'  
     lon_bounds        = 0.0, 360.0
     lat_bounds        = -90.0,  90.0
     /

Description of namelist variables:

+--------------------+--------------------+---------------------------------------------------------------------------+
| Contents           | Type               | Description                                                               |
+====================+====================+===========================================================================+
| input_file         | character(len=256) | Name of the netcdf soil moisture data file.                               |
+--------------------+--------------------+---------------------------------------------------------------------------+
| output_file        | character(len=256) | Name of the output observation sequence file.                             |
+--------------------+--------------------+---------------------------------------------------------------------------+
| lon_bounds         | real(r8)           | Longitude bounds. Observations outside these bounds are excluded from     |
|                    |                    | the output_file                                                           |
+--------------------+--------------------+---------------------------------------------------------------------------+
| lat_bounds         | real(r8)           | Latitude bounds. Observations outside these bounds are excluded from      |
|                    |                    | the output_file                                                           |
+--------------------+--------------------+---------------------------------------------------------------------------+

.. Important::

  The total error (instrument error + representativeness error) is used to calculate the observation error
  variance in the ``output_file``. The instrument error is taken from the data file variable 
  (sm_x_error), whereas the representativeness error is set to 0.1 within the ``LPRM_L3_to_obs``.


Data Source
-----------

The dataset (LPRM_TMI_NT_SOILM3: TMI/TRMM surface soil moisture (LPRM) L3 1 day 25km x 25km nighttime V001) can be
found `here. <https://disc.gsfc.nasa.gov/datasets/LPRM_TMI_NT_SOILM3_001/summary>`__





Leaf Area Index Observation Converter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Program ``netCDF_to_obs.f90``

Description
-----------

This dataset provides a global 0.25 degree x 0.25 degree gridded monthly 
mean leaf area index (LAI) climatology as averaged over the period from 
August 1981 to August 2015. The data were derived from the Advanced Very
High Resolution Radiometer (AVHRR) Global Inventory Modeling and Mapping 
Studies (GIMMS) LAI3g version 2, a bi-weekly data product from 1981 to 2015
(GIMMS-LAI3g version 2). The LAI3g version 2 (raw) data were first regridded
from 1/12 x 1/12 degree to 0.25 x 0.25 degree resolution, then processed to 
remove missing and unreasonable values, scaled to obtain LAI values, and the
bi-weekly LAI values were averaged for every month. Finally, the monthly 
long-term mean LAI (1981-2015) was calculated.


The Global Monthly Mean Leaf Area Index Climatology, (1981-2015) dataset
may be converted with the ``netCDF_to_obs`` program.   Since these are monthly means,
each timestep is read and output as their own observation sequence file that has the 
date and time appended to the filename. 


NAMELIST
--------

::

  &netCDF_to_obs_nml
     input_file        = 'LAI_mean_monthly_1981-2015.nc4'
     output_file_base  = 'obs_seq.out'  
     lon_bounds        = 0.0, 360.0
     lat_bounds        = -90.0,  90.0
     debug               = .FALSE.
     observation_varname          = 'LAI'
     observation_type             = 'GIMMS_LEAF_AREA_INDEX'
     obs_error_standard_deviation = 0.2     
     /


Description of namelist variables:

+------------------------------+--------------------+---------------------------------------------------------------------------+
| Contents                     | Type               | Description                                                               |
+==============================+====================+===========================================================================+
| input_file                   | character(len=256) | Name of the netcdf LAI data file.                                         |
+------------------------------+--------------------+---------------------------------------------------------------------------+
| output_file_base             | character(len=256) | Name of the output observation sequence file.                             |
+------------------------------+--------------------+---------------------------------------------------------------------------+
| lon_bounds                   | real(r8)           | Longitude bounds. Observations outside these bounds are excluded from     |
|                              |                    | the output file                                                           |
+------------------------------+--------------------+---------------------------------------------------------------------------+
| lat_bounds                   | real(r8)           | Latitude bounds. Observations outside these bounds are excluded from      |
|                              |                    | the output file                                                           |
+------------------------------+--------------------+---------------------------------------------------------------------------+
| debug                        | logical            | If .TRUE. prints out extra information on data file characteristics       |
+------------------------------+--------------------+---------------------------------------------------------------------------+
| observation_varname          | character(len=256) | Name of of the leaf area variable within the netcdf data file             |
+------------------------------+--------------------+---------------------------------------------------------------------------+
| observation_type             | character(len=256) | Name of the DART observation type                                         |
+------------------------------+--------------------+---------------------------------------------------------------------------+
| obs_error_standard_deviation | character(len=256) | The observation error standard deviation (not provided within data file)  |
+------------------------------+--------------------+---------------------------------------------------------------------------+  


Data Source
-----------

The Global Monthly Mean Leaf Area Index Climatology, (1981-2015) data can be found
`here. <https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1653>`__



