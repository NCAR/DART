Program ``convert_airs_L2`` 
===========================

.. caution:: 

   Before you begin: Installing the libraries needed to read these files can be
   fairly troublesome. The NASA Earthdata Data Access Services website is the
   `download site <https://wiki.earthdata.nasa.gov/display/DAS/Toolkit+Downloads>`__
   for the necessary libraries. An example build script (`AIRS/Build_HDF-EOS.sh`)
   is intended to provide some guidance.


Overview
--------

The Atmospheric Infrared Sounder (AIRS) is a facility instrument aboard the second 
Earth Observing System (EOS) polar-orbiting platform, EOS Aqua. In combination with 
the Advanced Microwave Sounding Unit (AMSU) and the Humidity Sounder for Brazil (HSB),
AIRS constitutes an innovative atmospheric sounding group of visible, infrared, and 
microwave sensors. AIRS data will be generated continuously. Global coverage will 
be obtained twice daily (day and night) on a 1:30pm sun synchronous orbit from a 
705-km altitude.

The AIRS Standard Retrieval Product consists of retrieved estimates of cloud 
and surface properties, plus profiles of retrieved temperature, water vapor, 
ozone, carbon monoxide and methane. Estimates of the errors associated with these 
quantities will also be part of the Standard Product. The temperature profile 
vertical resolution is 28 levels total between 1100 mb and 0.1 mb, while moisture 
profile is reported at 14 atmospheric layers between 1100 mb and 50 mb. The 
horizontal resolution is 50 km. An AIRS granule has been set as 6 minutes of data, 
30 footprints cross track by 45 lines along track. The Shortname for this product 
is AIRX2RET. (AIRS2RET is the same product but without the AMSU data.)

Atmospheric Infrared Sounder (AIRS) Level 2 observations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Several types of AIRS data, with varying levels of processing, are available.
The following descriptions are taken from the
`V5_Data_Release_UG <http://disc.sci.gsfc.nasa.gov/AIRS/documentation/v5_docs/AIRS_V5_Release_User_Docs/V5_Data_Release_UG.pdf>`__
document:

   The L1B data product includes geolocated, calibrated observed microwave, 
   infrared and visible/near infrared radiances, as well as Quality Assessment 
   (QA) data. The radiances are well calibrated; however, not all QA data have 
   been validated. Each product granule contains 6 minutes of data. Thus there 
   are 240 granules of each L1B product produced every day.

   The L2 data product includes geolocated, calibrated cloud-cleared radiances and 
   2-dimensional and 3-dimensional retrieved physical quantities (e.g., surface 
   properties and temperature, moisture, ozone, carbon monoxide and methane profiles 
   throughout the atmosphere). Each product granule contains 6 minutes of data. 
   Thus there are 240 granules of each L2 product produced every day.

   The L3 data are created from the L2 data product by binning them in 1°x1° grids.
   There are three products: daily, 8-day and monthly. Each product provides separate 
   ascending (daytime) and descending (nighttime) binned data sets.

The converter in this directory processes level 2 (L2) data files, using data 
set ``AIRS_DP`` and data product ``AIRX2RET`` or ``AIRS2RET`` without ``HSB`` 
(the instrument measuring humidity which failed).

Getting the data currently means putting in a start/stop time at 
`this web page <http://mirador.gsfc.nasa.gov/cgi-bin/mirador/homepageAlt.pl?keyword=AIRX2RET>`__.
The keyword is ``AIRX2RET`` and put in the time range of interest and optionally a 
geographic region. Each file contains 6 minutes of data, is about 2.3 Megabytes, 
and globally there are 240 files/day (about 550 Megabytes/day). There are additional 
options for getting only particular variables of interest, but the current reader 
expects whole files to be present. Depending on your connection to the internet, 
there are various options for downloading. We have chosen to download a ``wget`` 
script which is created by the web page after adding the selected files to a 'cart' 
and 'checking out'. The script has a series of ``wget`` commands which downloads 
each file, one at a time, which is run on the machine where you want the data.

convert_airs_L2.f90
-------------------

The ``convert_airs_L2`` converter is for **temperature and moisture retrievals** from
the L2 data. The temperature observations are at the 
corresponding vertical pressure levels. However, the moisture obs are the mean for 
the layer, so the location in the vertical is the midpoint, in log space, of the 
current layer and the layer above it. There is an alternative computation for the 
moisture across the layer which may be more accurate, but requires a forward 
operator subroutine to be written and for the observation to contain metadata. 
The observation could be defined with a layer top, in pressure, and a number of 
points to use for the integration across the layer. Then the forward operator would 
query the model at each of the N points in the vertical for a given horizontal 
location, and compute the mean moisture value. This code has not been implemented 
yet, and would require a different QTY_xxx to distinguish it from the simple 
location/value moisture obs. See the GPS non-local operator code for an example 
of how this would need to be implemented.

The temperature observations are located on standard levels; there is a single array 
of heights in each file and all temperature data is located on one of these levels. 
The moisture observations, however, are an integrated quantity for the space between 
the levels; in their terminology the fixed heights are 'levels' and the space between 
them are 'layers'. The current converter locates the moisture obs at the midpoint, 
in log space, between the levels.

The hdf files need to be downloaded from the data server, in any manner you choose. 
The converter program reads each hdf granule and outputs a DART obs_seq file 
containing up to 56700 observations. Only those with a quality control of 0 (Best) 
are kept. The resulting obs_seq files can be merged with the 
:doc:`../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool` into 
larger time periods.

It is possible to restrict the output observation sequence to contain data from a 
region of interest throught the use of the namelist parameters. If you need a region 
that spans the Prime Meridian lon1 can be a larger number than lon2, for example, 
a region from 300 E to 40 E and 60 S to 30 S (some of the South Atlantic), 
would be *lon1 = 300, lon2 = 40, lat1 = -60, lat2 = -30*.

The ``DART/observations/obs_converters/AIRS/shell_scripts`` directory includes scripts
(``download_L2.sh`` and ``oneday_down.sh``) that make use of the fact that the AIRS data 
is also archived on the NCAR HPSS (tape library) in daily tar files. 
``oneday_down.sh`` has options to download a day of granule files, convert them, merge them 
into daily files, and remove the original data files and repeat the process for any 
specified time period.


Namelist
--------

This namelist is read in a file called ``input.nml``. We adhere to the F90 
standard of starting a namelist with an ampersand '&' and terminating with a 
slash '/' for all our namelist input. Character strings that contain a '/' must be
enclosed in quotes to prevent them from prematurely terminating the namelist.
The default values are shown below. More realistic values are provided in
``AIRS/work/input.nml``

::

   &convert_airs_L2_nml
      l2_files           = ''
      l2_file_list       = ''
      outputfile         = ''
      lon1               =   0.0
      lon2               = 360.0
      lat1               = -90.0
      lat2               =  90.0
      min_MMR_threshold  = 1.0e-30
      top_pressure_level = 0.0001
      cross_track_thin   = 0
      along_track_thin   = 0
      use_NCEP_errs      = .false.
      version            = 6
   /

| 

.. container::

   +--------------------+------------------------+--------------------------------------------------------------+
   | Contents           | Type                   | Description                                                  |
   +====================+========================+==============================================================+
   | l2_files           | character(len=256),    | A list of one or more names of the HDF file(s) to read,      |
   |                    | dimension(512)         | NOT including the directory. If multiple files are listed,   |
   |                    |                        | each will be read and the results will be placed in a        |
   |                    |                        | separate file with an output filename constructed based on   |
   |                    |                        | the input filename.                                          |
   +--------------------+------------------------+--------------------------------------------------------------+
   | l2_file_list       | character(len=256)     | The name of an ascii text file which contains one filename   |
   |                    |                        | per line, NOT including the directory. Each file will be     |
   |                    |                        | read and the observations converted into an output file      |
   |                    |                        | where the output filename is based on the input filename.    |
   |                    |                        | Only one of 'l2_files' and 'l2_file_list' can be             |
   |                    |                        | specified. The other must be ' ' (empty).                    |
   +--------------------+------------------------+--------------------------------------------------------------+
   | outputfile         | character(len=256)     | The name of the output observation sequence file.            |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lon1               | real(r8)               | the West-most longitude of interest in degrees. [0.0, 360]   |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lon2               | real(r8)               | the East-most longitude of interest in degrees. [0.0, 360]   |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lat1               | real(r8)               | the South-most latitude of interest in degrees. [-90.0,90.0] |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lat2               | real(r8)               | the North-most latitude of interest in degrees. [-90.0,90.0] |
   +--------------------+------------------------+--------------------------------------------------------------+
   | min_MMR_threshold  | real(r8)               | The data files contains 'Retrieved Water Vapor Mass Mixing   |
   |                    |                        | Ratio'. This is the minimum threshold, in gm/kg, that will   |
   |                    |                        | be converted into a specific humidity observation.           |
   +--------------------+------------------------+--------------------------------------------------------------+
   | top_pressure_level | real(r8)               | The highest pressure level of interest (in mb).              |
   +--------------------+------------------------+--------------------------------------------------------------+
   | cross_track_thin   | integer                | provides ability to thin the data by keeping every Nth data  |
   |                    |                        | value in the cross-track scan.   [0,30]                      |
   |                    |                        | e.g. 3 == keep every third value. 0 is no thinning.          |
   +--------------------+------------------------+--------------------------------------------------------------+
   | along_track_thin   | integer                | provides ability to thin the data by keeping every Nth data  |
   |                    |                        | value in the along-track scan.   [0,45]                      |
   |                    |                        | e.g. 4 == keep only every 4th row. 0 is no thinning.         |
   +--------------------+------------------------+--------------------------------------------------------------+
   | use_NCEP_errs      | logical                | if .true. use the maximum observation error from either the  |
   |                    |                        | granule or the NCEP equivalent (from ``obs_error_mod.f90``)  |
   +--------------------+------------------------+--------------------------------------------------------------+
   | version            | integer                | The AIRS file format version.                                |
   +--------------------+------------------------+--------------------------------------------------------------+


Dependencies
~~~~~~~~~~~~

See the :doc:`Dependencies Section<./README>` of the AIRS/README.

Known Bugs
~~~~~~~~~~

Earlier versions of this converter mistakenly put the moisture obs
at level heights, in the same location as the temperature observations.
The moisture observations are in fact an integrated value across the
distance between two levels.
This means the location was shifted 1/2 level in the vertical from 
the center of the layer.  The fixed converter outputs the location
at the center, in log space, of each layer.


Future Plans
~~~~~~~~~~~~
If a more accurate moisture observation was needed, the observation value
could be computed by actually integrating multiple values between the levels.
At this point it doesn't seem necessary.
 
