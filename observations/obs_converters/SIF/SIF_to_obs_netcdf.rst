PROGRAM ``SIF_to_obs_netcdf``
=============================

Overview
--------

Harmonized SIF data product to DART observation sequence converter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This routine converts a harmonized satellite SIF product 
(Harmonized long-term SIF; SIF005) to DART ``obs_seq`` format.
The SIF product is described by
`JPL <https://climatesciences.jpl.nasa.gov/sif/download-data/level-3/>`__ 
and combines GOME-2 and SCIAMACHY SIF retrievals, along with MODIS data
to produce a single continuous, monthly, 0.05 degree SIF data set.  
See `Wen et al., 2020 RSE <https://doi.org/10.1016/j.rse.2020.111644>`__ 
for a more detailed description.  
The conversion script was designed and tested for version SIF005v2. 
Download instructions can be found in the `Data Sources`_ section below.

This SIF data product also comes with its own uncertainty value, and quality 
control flag described below.  Namelist options also include a wavelength option
(e.g. 740 nm or 755 nm) to specify the location the SIF value is centered upon. 


Standard workflow:

#. Download the Level 3 data for the months of interest. Years 2002-2018 are available
   as of 5/18/21.  (see `Data Sources`_ below)
#. Make note of the SIF wavelength the data is centered upon. This information is 
   included in the SIF variable of netcdf file ``SIF_740_daily_corr``  
#. Build the DART executables with support for land observations. This is done by running 
   ``quickbuild.sh`` with ``obs_def_land_mod.f90`` in the list of ``input_files`` for 
   ``preprocess_nml``.
#. Provide basic information via the ``SIF_to_obs_netcdf_nml`` (e.g. verbose, wavelength)
#. Convert single or multiple SIF netcdf data files using ``SIF_to_obs_netcdf``. Converting
   one file at a time results in better memory management, but this is unlikely to be an
   issue in most cases.
#. Combine all output files for the region and timeframe of interest into one file using
   :doc:`../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool`

For some models (CLM, for example), it is required to reorganize the observation sequence 
files into a series of files that contains ONLY the observations for each assimilation. 
This can be achieved with the `makedaily.sh` script which can be found in 
the `DART/models/clm/shell_scripts` directory.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' 
and terminate with a slash '/'.  Character strings that contain a '/' must be enclosed in
quotes to prevent them from prematurely terminating the namelist.

::

   &SIF_to_obs_netcdf_nml
      input_file_list = 'SIF.input.txt',
      verbose         = 0
      wavelength      = 740
      /


.. container::

   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | Contents        | Type               | Description                                                                 |
   +=================+====================+=============================================================================+
   | input_file_list | character(len=256) | Name of the Level 3 netcdf containing with SIF data. This may be a          |
   |                 |                    | relative or absolute filename.                                              |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | verbose         | integer            | Print more/less information during the ``SIF_to_obs_netcdf`` execution.     |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | wavelength      | integer            | The wavelength at which SIF irradiance is centered upon (e.g. 740 nm)       | 
   +-----------------+--------------------+-----------------------------------------------------------------------------+


Data Sources
------------

The datasets are available from the
`Cornell University Box service <https://cornell.app.box.com/s/gkp4moy4grvqsus1q5oz7u5lc30i7o41/folder/100438579357>`__,
and have names like:

``SIF005_YYYYMM.nc``, ``SIF005_200504.nc``, ``SIF005_201808.nc`` 

The Level 3 SIF product is provided within netcdf files for monthly average values
from 2002-2018. This ``SIF_obs_to_netcdf`` converter was tested using SIF005v2 files,
although older (SIF005v1) and newer (SIF005v2.2) versions exist with similar format.

The data product variables are provided in global, gridded (lat/lon) format:

+---------------+----------------------+------------------------------+--------------------------+--------------------------------+-------------+
| Units         | Variable             | Description                  | Observation TYPE         | DART QUANTITY                  | DART units  |
+===============+======================+==============================+==========================+================================+=============+
| mW/m^2/nm/sr  | SIF_740_daily_corr   | Solar Induced                | HARMONIZED_SIF           | QTY_SOLAR_INDUCED_FLUORESCENCE | mW/m^2/nm/sr|
|               |                      | Fluorescence Irradiance      |                          |                                |             |
+---------------+----------------------+------------------------------+--------------------------+--------------------------------+-------------+
| mW/m^2/nm/sr  | SIF_740_daily_corr_SD| Solar Induced Fluorescence   |   N/A                    |    N/A                         | mW/m^2/nm/sr|
|               |                      | Irradiance Standard Deviation|                          |                                |             |
+---------------+----------------------+------------------------------+--------------------------+--------------------------------+-------------+
| See below     | EVI_Quality          | MODIS EVI Quality Flag       |   N/A                    |    N/A                         | See below   |
+---------------+----------------------+------------------------------+--------------------------+--------------------------------+-------------+
| degrees       | lat                  | latitude                     |   N/A                    |    N/A                         | radians     |
+---------------+----------------------+------------------------------+--------------------------+--------------------------------+-------------+
| degrees       | lon                  | longitude                    |   N/A                    |    N/A                         | radians     |
+---------------+----------------------+------------------------------+--------------------------+--------------------------------+-------------+



The ``SIF_740_daily_corr`` value is the SIF satellite derived irradiance value. 
It is most closely related to the 'top of the vegetation canopy' emitted SIF as simulated
from land surface models.  This is distinct from 'leaf-level' emitted SIF.

The ``SIF_740_daily_corr_SD`` value is an algorithm based uncertainty estimate 
provided by the data product providers.  It is most closely related to instrument 
uncertainty inherent to the SIF retrievals and does not account for
representativeness error when compared to the simulated SIF from a land surface model.
We recommend this uncertainty value be used as a minimum baseline when performing
data assimilation.

The ``EVI_Quality`` is a data quality estimate for the ``SIF_740_daily_corr``.
The ``EVI_Quality`` is derived from the MODIS retrieval of EVI (enhanced vegetation index)
which is one of the explanatory variables used in the algorithm to calculate 
``SIF_740_daily_corr``.  The ``EVI_Quality`` is an integer (representing a 16 bit field)
that evaluates quality through 9 parameters that include VI (Vegetation Index) Quality, 
VI Usefulness, Aerosol Quantity, Adjacent Cloud Detection, Atmosphere BRDF correction, 
Mixed Clouds, Land/Water Mask, possible snow/ice, possible shadow.  See Table 5 of the
`MODIS Vegetation Index Users Guide <https://lpdaac.usgs.gov/documents/103/MOD13_User_Guide_V6.pdf>`__ 
for more information.  

The DART-compatible QC value assigned to the `obs_seq.out` uses the criteria from 
the MODIS EVI Quality and EVI Usefulness only.  The DART-compatible QC is based on
NCEP-like error codes and ``SIF_to_obs_netcdf`` assigns values as follows:

+--------------------------+
| 0  = best quality        |
+--------------------------+
| 1  = less quality        |
+--------------------------+
| "........."              |
+--------------------------+
| 17 = least quality       |
+--------------------------+     
| 50  = faulty, no utility |
+--------------------------+

The `input_qc_threshold` namelist value can be used to test whether or not lesser 
quality observations improve the result or not.  Thus, all observations (except those
that are defined as faulty/no utility) are included in `obs_seq.out` and the exclusion
of observations is left up to the user based upon the `input_qc_threshold`.

The qc value assignment is such where values given an EVI quality value of 
'good' (00), are assigned a QC from 1-7 based on the EVI Quality Usefulness
Parameter (see table below).  Values where the 'EVI is produced, but should be checked
with additional QA' (01) are assigned a QC from 10-17. Anything with an
EVI Quality Usefulness Parameter of '1101' or higher is given a QC of 50 (or more) and
is currently skipped and **not written** to the output observation sequence file.

+------------------------------------------------------+----+-----------------------------+----+---------------------------+
| EVI Quality Usefulness Parameter                     | QC | EVI Quality Value (00)      | QC | EVI Quality Value (01)    |
+======+===============================================+====+=============================+====+===========================+
| 0000 |  Highest quality                              | 0  | Highest quality             | 10 | Decreasing quality        |
+------+-----------------------------------------------+----+-----------------------------+----+---------------------------+
| 0001 | Lower quality                                 | 1  | Lower quality               | 11 | Decreasing quality        |
+------+-----------------------------------------------+----+-----------------------------+----+---------------------------+
| 0010 | Decreasing quality                            | 2  | Decreasing quality          | 12 | Decreasing quality        |
+------+-----------------------------------------------+----+-----------------------------+----+---------------------------+
| 0100 | Decreasing quality                            | 3  | Decreasing quality          | 13 | Decreasing quality        |
+------+-----------------------------------------------+----+-----------------------------+----+---------------------------+
| 1000 | Decreasing quality                            | 4  | Decreasing quality          | 14 | Decreasing quality        |
+------+-----------------------------------------------+----+-----------------------------+----+---------------------------+
| 1001 | Decreasing quality                            | 5  | Decreasing quality          | 15 | Decreasing quality        |
+------+-----------------------------------------------+----+-----------------------------+----+---------------------------+
| 1010 | Decreasing quality                            | 6  | Decreasing quality          | 16 | Decreasing quality        |
+------+-----------------------------------------------+----+-----------------------------+----+---------------------------+
| 1100 | Lowest     quality                            | 7  | Decreasing quality          | 17 | Least quality             |
+------+-----------------------------------------------+----+-----------------------------+----+---------------------------+
| 1101 | Quality so low that it is not useful          | 50 | Not used                    | 50 | Not used                  |
+------+-----------------------------------------------+----+-----------------------------+----+---------------------------+
| 1110 | L1B data faulty                               | 50 | Not used                    | 50 | Not used                  |
+------+-----------------------------------------------+----+-----------------------------+----+---------------------------+
| 1111 | Not useful for any other reason/not processed | 50 | Not used                    | 50 | Not used                  |
+------+-----------------------------------------------+----+-----------------------------+----+---------------------------+



Citation
--------

Wen, J., P. Köhler, G. Duveiller, N. C. Parazoo, T. S. Magney, G. Hooker, L. Yu, 
C. Y. Chang, and Y. Sun. "A framework for harmonizing multiple satellite instruments 
to generate a long-term global high spatial-resolution solar-induced chlorophyll 
fluorescence (SIF)." Remote Sensing of Environment 239 (2020): 
111644.https://doi.org/10.1016/j.rse.2020.111644

