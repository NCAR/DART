Program ``convert_amsu_L1``
===========================

Overview
---------

There is a little bit of confusing history to be aware of for AMSU/A:

https://en.wikipedia.org/wiki/Advanced_microwave_sounding_unit#History

AMSU/A was flown on NOAA 15-17. It is also on the Aqua satellite (that
also houses AIRS) as well as the European MetOp. It has been replaced by
ATMS on NOAA-20.

The datset of interest is: “AIRS/Aqua L1B AMSU (A1/A2) geolocated and
calibrated brightness temperatures V005 (AIRABRAD) at GES DISC” The
*short name* for this dataset is ‘AIRABRAD’

The introductory paragraph for the dataset is:

   Version 5 is the current version of the data set.tmospheric Infrared
   Sounder (AIRS) is a grating spectrometer (R = 1200) aboard the second
   Earth Observing System (EOS) polar-orbiting platform, EOS Aqua. In
   combination with the Advanced Microwave Sounding Unit (AMSU) and the
   Humidity Sounder for Brazil (HSB), AIRS constitutes an innovative
   atmospheric sounding group of visible, infrared, and microwave
   sensors. The AMSU-A instrument is co-aligned with AIRS so that
   successive blocks of 3 x 3 AIRS footprints are contained within one
   AMSU-A footprint. AMSU-A is primarily a temperature sounder that
   provides atmospheric information in the presence of clouds, which can
   be used to correct the AIRS infrared measurements for the effects of
   clouds. This is possible because non-precipitating clouds are for the
   most part transparent to microwave radiation, in contrast to visible
   and infrared radiation which are strongly scattered and absorbed by
   clouds. AMSU-A1 has 13 channels from 50 - 90 GHz and AMSU-A2 has 2
   channels from 23 - 32 GHz. The AIRABRAD_005 products are stored in
   files (often referred to as “granules”) that contain 6 minutes of
   data, 30 footprints across track by 45 lines along track.

The citation information for this dataset is:

   Title: AIRS/Aqua L1B AMSU (A1/A2) geolocated and calibrated
   brightness temperatures V005 Version: 005 Creator: AIRS project
   Publisher: Goddard Earth Sciences Data and Information Services
   Center (GES DISC) Release Date: 2007-07-26T00:00:00.000Z Linkage:
   https://disc.gsfc.nasa.gov/datacollection/AIRABRAD_005.html

NASA provides a `README.AIRABRAD.pdf <https://docserver.gesdisc.eosdis.nasa.gov/repository/Mission/AIRS/3.3_ScienceDataProductDocumentation/3.3.4_ProductGenerationAlgorithms/README.AIRABRAD.pdf>`__
through the Goddard Earth Sciences Data and Information Services Center.

Advanced Microwave Sounding Unit (AMSU-A) L1B Brightness Temperatures
---------------------------------------------------------------------

Converting AMSU_L1 observations is a two-step process:

- convert the data from HDF to netCDF
- run ``convert_amsua_L1`` to convert the netCDF file to an obs_seq file.

.. note::

   The native HDF-EOS2 format files must be converted to netCDF before running
   ``convert_amsua_L1``.

To convert from HDF-EOS2 to netCDF use the
`h4tonccf_nc4 <http://hdfeos.org/software/h4cflib.php>` from the HDF-EOS
tools.


The netCDF files have two global
attributes that are exceedingly large and uninformative. If needed you can
remove these attributes, you can use the ``ncatted`` command from
`NCO <http://nco.sourceforge.net/nco.html>`_.

::

   h4tonccf_nc4 AIRS.2019.06.22.236.L1B.AMSU_Rad.v5.0.0.0.G19174110442.hdf bob.nc
   ncatted -a coremetadata,global,d,,, -a StructMetadata_0,global,d,,, bob.nc bill.nc




As you can imagine, you need to download each satellite’s data in a
different way. Also, just for your information, AMSU/B has been replaced
on newer satellites by MHS and HSB, but especially MHS is almost
identical.

Namelist
~~~~~~~~

``convert_amsua_L1`` makes use of :doc:`../../forward_operators/obs_def_rttov_mod`
Only two &obs_def_rttov_nml options are required when converting
the observations: *use_zeeman* and *rttov_sensor_db_file*.

Be aware that if the RTTOV namelist option ``use_zeeman = .true.``
certain metadata must be available in the observation. This is not fully
implemented in the AMSU-A observation converter. For more information,
please see GitHub Issue 99 “`AIRS AMSUA observation converter … Zeeman
coefficients and channels <https://github.com/NCAR/DART/issues/99>`__”

Namelists are read in a file called ``input.nml``. We adhere to the F90 
standard of starting a namelist with an ampersand '&' and terminating with a 
slash '/' for all our namelist input. Character strings that contain a '/' must be
enclosed in quotes to prevent them from prematurely terminating the namelist.
The default values are shown below. More realistic values are provided in
``AIRS/work/input.nml``

::

   &convert_amsua_L1_nml
      l1_files           = ''
      l1_file_list       = ''
      outputfile         = ''
      append_output      = .false.
      channel_list       = 'null'
      along_track_thin   = 0
      cross_track_thin   = 0
      lon1               =   0.0
      lon2               = 360.0
      lat1               = -90.0
      lat2               =  90.0
      verbose            = 0
   /



.. container::

   +--------------------+------------------------+--------------------------------------------------------------+
   | Contents           | Type                   | Description                                                  |
   +====================+========================+==============================================================+
   | l1_files           | character(len=256),    | A list of one or more names of the netCDF file(s) to read.   |
   |                    | dimension(512)         |                                                              |
   +--------------------+------------------------+--------------------------------------------------------------+
   | l1_file_list       | character(len=256)     | The name of an ascii text file which contains one filename   |
   |                    |                        | per line. Each file will be read and the observations        |
   |                    |                        | converted into a single output file.                         |
   |                    |                        | Only one of 'l1_files' and 'l1_file_list' can be             |
   |                    |                        | specified. The other must be ' ' (empty).                    |
   +--------------------+------------------------+--------------------------------------------------------------+
   | outputfile         | character(len=256)     | The name of the output observation sequence file.            |
   +--------------------+------------------------+--------------------------------------------------------------+
   | append_output      | logical                | If the output observation sequence file exists it is possible|
   |                    |                        | to add to it. The observations are added consistent with the |
   |                    |                        | paradigm that the observation linked list will be traversed  |
   |                    |                        | in temporally-ascending fashion, no matter the physical      |
   |                    |                        | location of the observation in the file. ``.true.`` adds the |
   |                    |                        | new observations to the existing file, ``.false.`` will      |
   |                    |                        | cause an existing output file to be overwritten.             |
   +--------------------+------------------------+--------------------------------------------------------------+
   | channel_list       | character(len=8),      | The AMSU channels desired.                                   |
   |                    | dimension(15)          | See the table below for valid input.                         |
   +--------------------+------------------------+--------------------------------------------------------------+
   | along_track_thin   | integer                | provides ability to thin the data by keeping every Nth data  |
   |                    |                        | value in the along-track scan.   [0,45]                      |
   |                    |                        | e.g. 4 == keep only every 4th row. 0 is no thinning.         |
   +--------------------+------------------------+--------------------------------------------------------------+
   | cross_track_thin   | integer                | provides ability to thin the data by keeping every Nth data  |
   |                    |                        | value in the cross-track scan.   [0,30]                      |
   |                    |                        | e.g. 3 == keep every third value. 0 is no thinning.          |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lon1               | real(r8)               | the West-most longitude of interest in degrees. [0.0, 360]   |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lon2               | real(r8)               | the East-most longitude of interest in degrees. [0.0, 360]   |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lat1               | real(r8)               | the South-most latitude of interest in degrees. [-90.0,90.0] |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lat2               | real(r8)               | the North-most latitude of interest in degrees. [-90.0,90.0] |
   +--------------------+------------------------+--------------------------------------------------------------+
   | verbose            | integer                | Controls the amount of run-time output.                      |
   |                    |                        | 0 == bare minimum. 3 is very verbose.                        |
   |                    |                        | Only use 3 if converting one or two files for testing.       |
   +--------------------+------------------------+--------------------------------------------------------------+


Channel Specification
~~~~~~~~~~~~~~~~~~~~~

   "AMSU-A primarily provides temperature soundings. It is a 15-channel microwave
   temperature sounder implemented as two independently operated modules. Module 1
   (AMSU-A1) has 12 channels in the 50-58 GHz oxygen absorption band which provide
   the primary temperature sounding capabilities and 1 channel at 89 GHz which provides
   surface and moisture information. Module 2 (AMSU-A2) has 2 channels: one at 23.8
   GHz and one at 31.4 GHz which provide surface and moisture information (total
   precipitable water and cloud liquid water)."


To facilitate the selection of channels, either the 'Integer' or 'String' values
may be used to specify ``channel_list``. The 'Documentation' and 'netCDF' values
are provided for reference only. The 'Documentation' values are from the 
`README.AIRABRAD.pdf <https://docserver.gesdisc.eosdis.nasa.gov/repository/Mission/AIRS/3.3_ScienceDataProductDocumentation/3.3.4_ProductGenerationAlgorithms/README.AIRABRAD.pdf>`__ document.


.. container::


   +---------+---------+---------------+---------------+
   |         |         | Documentation | netCDF        |
   | Integer | String  | Frequency     | `center_freq` |
   +=========+=========+===============+===============+
   | Module 2 - surface and moisture information       |
   +---------+---------+---------------+---------------+
   | 1       | 'A2-1'  | 23.8          | 23.8          |
   +---------+---------+---------------+---------------+
   | 2       | 'A2-2'  | 31.4          | 31.4          |
   +---------+---------+---------------+---------------+
   | Module 1 - primary temperature sounding capability|
   +---------+---------+---------------+---------------+
   | 3       | 'A1-1'  | 50.3          | 50.3          |
   +---------+---------+---------------+---------------+
   | 4       | 'A1-2'  | 52.8          | 52.8          |
   +---------+---------+---------------+---------------+
   | 5       | 'A1-3'  | 53.596        | 53.596        |
   +---------+---------+---------------+---------------+
   | 6       | 'A1-4'  | 54.4          | 54.4          |
   +---------+---------+---------------+---------------+
   | 7       | 'A1-5'  | 54.94         | 54.94         |
   +---------+---------+---------------+---------------+
   | 8       | 'A1-6'  | 55.5          | 55.5          |
   +---------+---------+---------------+---------------+
   | 9       | 'A1-7'  | 57.29034      | 57.29034      |
   +---------+---------+---------------+---------------+
   | 10      | 'A1-8'  |               | 57.29034      |
   +---------+---------+---------------+---------------+
   | 11      | 'A1-9'  |               | 57.29034      |
   +---------+---------+---------------+---------------+
   | 12      | 'A1-10' |               | 57.29034      |
   +---------+---------+---------------+---------------+
   | 13      | 'A1-11' |               | 57.29034      |
   +---------+---------+---------------+---------------+
   | 14      | 'A1-12' |               | 57.29034      |
   +---------+---------+---------------+---------------+
   | 15      | 'A1-13' | 89            | 89            |
   +---------+---------+---------------+---------------+


.. _instructions-to-download-the-airabrad-dataset-1:

Instructions to download the AIRABRAD dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Go to https://earthdata.nasa.gov
2. Log in (or create an account if necessary)
3. Search for AIRABRAD
4. Scroll down past datasets to “Matching results.”

-  Follow the link to “AIRS/Aqua L1B AMSU (A1/A2) geolocated and
   calibrated brightness temperatures V005 (AIRABRAD) at GES DISC”

5. You should now be at
   ‘https://cmr.earthdata.nasa.gov/search/concepts/C1243477366-GES_DISC.html’
   (unless they’ve changed the site).

-  Select the ‘Download data’ tab
-  Select ‘Earthdata search’
-  Select the AIRS link under ‘Matching datasets’ (I have not tested the
   NRT products)

6. You can now select ‘Granule filters’ to choose your start and end
   dates.
7. Select the granules you want, then click ‘download all’ and 
   'download data’
8. Click download access script
9. Follow the instructions on that page to download the data.


| Each granule is about 560K and has names like

::

   AIRS.2019.06.22.236.L1B.AMSU_Rad.v5.0.0.0.G19174110442.hdf


Actually converting to netCDF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




