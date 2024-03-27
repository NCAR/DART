Program ``convert_airs_L2`` 
===========================

Overview
--------

The Atmospheric Infrared Sounder `(AIRS) <http://airs.jpl.nasa.gov/>`_ is a facility
instrument aboard the second Earth Observing System (EOS) polar-orbiting platform
`Aqua <http://aqua.nasa.gov>`_. Aqua is one of a group of satellites flying close
together in a polar orbit, collectively known as the “A-train”. In combination with
the Advanced Microwave Sounding Unit (AMSU) and the Humidity Sounder for Brazil (HSB),
AIRS constitutes an innovative atmospheric sounding group of visible, infrared, and 
microwave sensors. AIRS data will be generated continuously. 

The AIRS Standard Retrieval Product consists of retrieved estimates of cloud 
and surface properties, plus profiles of retrieved temperature, water vapor, 
ozone, carbon monoxide and methane. Estimates of the errors associated with these 
quantities will also be part of the Standard Product. The temperature profile 
vertical resolution is 28 levels total between 1100 and 0.1 hPa, while moisture 
profile is reported at 14 atmospheric layers between 1100 hPa and 50 hPa. The 
horizontal resolution is 50 km. An AIRS granule has been set as 6 minutes of data, 
There are 240 granules per day, with orbit repeat cycle of approximately 16 days.

Overview of L1-L3 Atmospheric Infrared Sounder (AIRS) Observations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``convert_airs_L2`` converter is designed specifically for 
**temperature and moisture retrievals for L2 observations** only. 
For reference, we provide a brief description of the L1-L3 AIRS data
products below. For more detailed information please see the 
`AIRS documentation page: <https://disc.gsfc.nasa.gov/information/documents?title=AIRS%20Documentation>`_


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


Downloading Atmospheric Infrared Sounder (AIRS) L2 Observations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several data file types and versions that contain L2
observations for temperature and moisture profiles.  **We recommend the use of
the AIRS2RET version 7 (AIRS2RETv7) data product.**  The ``AIRS2RET`` data (AIRS data only)
product is preferred to the ``AIRX2RET`` (AIRS/AMSU data) because the radiometric
noise in several AMSU channels increased (since June 2007) degrading the
``AIRX2RET`` product. Furthermore, the version 7 product is higher quality than version 6
because of an improved retrieval algorithm leading to significantly improved RMSE and bias statistics.
See the `AIRS2RETv7 documentation <https://disc.gsfc.nasa.gov/datasets/AIRS2RET_7.0/summary>`_ 
for more information.

Although we recommend ``AIRS2RETv7``, the  ``convert_airs_L2`` converter is compatible
with ``AIRS2RET`` and ``AIRX2RET`` versions 5-7. Version 5 is no longer available
within the GES DISC database. For more information on these data products see the
links below:

- `AIRS2RETv6 <https://disc.gsfc.nasa.gov/datasets/AIRS2RET_006/summary>`_
- `AIRX2RETv6 <https://disc.gsfc.nasa.gov/datasets/AIRX2RET_006/summary>`_
- `AIRX2RETv7 <https://disc.gsfc.nasa.gov/datasets/AIRX2RET_007/summary>`_

The AIRS data is located within the Goddard Earth Sciences Data and Information
Services Center (GES DISC) `located here <https://disc.gsfc.nasa.gov/>`_. You need
to create an Earthdata account before you can download data. As an example, to 
access the AIRS2RETv7 data, search on keyword ``AIRS2RET`` and locate
the AIRS2RET 7.0 data set within your search results. The full name is listed as
**Aqua/AIRS L2 Standard Physical Retrieval (AIRS-only) V7.0 (AIRS2RET)**. Next, click on the 
``Subset/Get Data`` link within the `Data Access` portion of the webpage. This will
bring up a separate window that allows you to refine your search results 
by 1) ``Refine range (time)`` and 2) ``Refine region (spatial)``. 

There are various options for downloading, however, the most straightforward approach
for macOS and Linux users is to use the ``wget`` command.  The ``download instructions``
provide the proper wget flags/options.  The ``Download Links List`` provides 
the AIRS file list based on your search results. 

convert_airs_L2.f90
-------------------

The ``convert_airs_L2`` converter is for **temperature and moisture retrievals** from
the L2 data. 
The vertical coordinate is pressure.
The temperature observations are defined at standard pressure levels (see Overview).
Those are defined in each file by the array 
'StdPressureLev:L2_Standard_atmospheric&surface_product'.
Between 2 levels is a "layer".
A moisture observation is an average across the layer
and is defined at the midpoint (in log(pressure)) of the layer.
This choice makes half of the mass of the layer above the midpoint and half below.
The midpoints are defined in 'H2OPressureLay:L2_Standard_atmospheric&surface_product'.

There is an alternative computation for the moisture across the layer
which may be more accurate, but requires a forward operator subroutine
to be written and for the observation converter to include additional metadata
to support this forward operator.
For more information see the Future Plans section below.

The converter program reads each AIRS hdf file granule and outputs a DART obs_seq file 
containing up to 56700 observations. Only those with a quality control of 0 (Best) 
are kept. The resulting obs_seq files can be merged with the 
:ref:`obs sequence tool` into
larger time periods.

During the excecution of the obs converter, It is possible to restrict the output
observation sequence to contain data from a region of interest throught the use of
the namelist parameters (described in Namelist section below). If you need a region
that spans the Prime Meridian, ``lon1`` can be a larger number than ``lon2``. 
For example, a region from 300 E to 40 E and 60 S to 30 S (some of the South Atlantic), 
would be ``lon1 = 300``, ``lon2 = 40``, ``lat1 = -60``, ``lat2 = -30``.


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
      version            = 7
   /

| 

.. container::

   +--------------------+------------------------+--------------------------------------------------------------+
   | Contents           | Type                   | Description                                                  |
   +====================+========================+==============================================================+
   | l2_files           | character(len=256),    | A list of one or more names of the HDF file(s) to read.      |
   |                    | dimension(512)         | If multiple files are listed, each will be read and          |
   |                    |                        | the results will be placed in a separate file with           |
   |                    |                        | an output filename constructed based on the input filename.  |
   +--------------------+------------------------+--------------------------------------------------------------+
   | l2_file_list       | character(len=256)     | The name of an ascii text file which contains one filename   |
   |                    |                        | per line.  Each file will be read and the observations       |
   |                    |                        | converted into an output file where the output filename      |
   |                    |                        | is based on the input filename.                              |
   |                    |                        | Only one of 'l2_files' and 'l2_file_list' can be  specified. |
   |                    |                        | The other must be ' ' (empty).                               |
   +--------------------+------------------------+--------------------------------------------------------------+
   | outputfile         | character(len=256)     | The name of the output observation sequence file.            |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lon1               | real(r8)               | The West-most longitude of interest in degrees. [0.0, 360]   |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lon2               | real(r8)               | The East-most longitude of interest in degrees. [0.0, 360]   |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lat1               | real(r8)               | The South-most latitude of interest in degrees. [-90.0,90.0] |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lat2               | real(r8)               | The North-most latitude of interest in degrees. [-90.0,90.0] |
   +--------------------+------------------------+--------------------------------------------------------------+
   | min_MMR_threshold  | real(r8)               | The data files contain 'Retrieved Water Vapor Mass Mixing    |
   |                    |                        | Ratio'. This is the minimum threshold (g/kg) that will       |
   |                    |                        | be converted into a specific humidity observation (kg/kg).   |
   +--------------------+------------------------+--------------------------------------------------------------+
   | top_pressure_level | real(r8)               | The highest pressure level of interest (in hPa).             |
   +--------------------+------------------------+--------------------------------------------------------------+
   | cross_track_thin   | integer                | Provides ability to thin the data by keeping every Nth data  |
   |                    |                        | value in the cross-track scan.   [0,30]                      |
   |                    |                        | e.g. 3 == keep every third value. 0 is no thinning.          |
   +--------------------+------------------------+--------------------------------------------------------------+
   | along_track_thin   | integer                | Provides ability to thin the data by keeping every Nth data  |
   |                    |                        | value in the along-track scan.   [0,45]                      |
   |                    |                        | e.g. 4 == keep only every 4th row. 0 is no thinning.         |
   +--------------------+------------------------+--------------------------------------------------------------+
   | use_NCEP_errs      | logical                | If .true. use the maximum observation error from either the  |
   |                    |                        | granule or the NCEP equivalent (from ``obs_error_mod.f90``)  |
   +--------------------+------------------------+--------------------------------------------------------------+
   | version            | integer                | The AIRS file format version. Version 7 is recommended, but  |
   |                    |                        | the converter is compatible with versions 5-7.               | 
   +--------------------+------------------------+--------------------------------------------------------------+

   | Included here are some example values for the l2_files namelist option.
   | Version 5 file: ``l2_files = '../data/AIRS.2007.11.01.001.L2.RetStd.v5.2.2.0.G08078150655.hdf'``
   | Version 6 file: ``l2_files = '../data/AIRS.2017.01.01.110.L2.RetStd_IR.v6.0.31.1.G19058124823.hdf'``
   | Version 7 file: ``l2_files = '../data/AIRS.2020.06.15.224.L2.RetStd_IR.v7.0.4.0.G20330033505.hdf'``

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
The observation could be defined with a layer top, in pressure, and a number of
points to use for the integration across the layer. Then the forward operator would
query the model at each of the N points in the vertical for a given horizontal
location, and compute the mean moisture value. This code has not been implemented
yet, and would require a different QTY_xxx to distinguish it from the simple
location/value moisture obs. The observation converter would also have to bring
in moisture observation metadata for this forward operator. See the 
GPS non-local operator code (:ref:`gps`)  for an example of how this
would need to be implemented.
 
