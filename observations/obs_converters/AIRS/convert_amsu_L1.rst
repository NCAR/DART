Program ``convert_amsu_L1``
===========================

Overview
---------

The following is an excerpt from the AIRS L1B AMSU-A documentation.
The complete documentation provided by the Goddard Earth Sciences Data 
and Information Services Center `(GES DISC) <https://disc.gsfc.nasa.gov/>`_ 
can be within the Documentation->README Document `found here <https://disc.gsfc.nasa.gov/datasets/AIRABRAD_005/summary>`_.

The Atmospheric Infrared Sounder (AIRS) Version 5 Level 1B Advanced Microwave
Sounding Unit (AMSU)-A Products (AIRABRAD) contain calibrated and 
geolocated brightness temperatures in degrees Kelvin. AIRABRAD_NRT (Near Real Time)
products are also available within ~3 hours of observations globally and stay for
about 5 days from the time they are generated. This data set is generated from 
AMSU-A level 1A digital numbers (DN) and contains 15 microwave channels in the
50-90 GHz and 23-32 GHz regions of the spectrum. A day's worth of data is divided
into 240 scenes (granules), each of 6 minute duration. An AMSU-A scene contains 
30 cross-track footprints in each of 45 along-track scanlines, for a total of 
45 x 30 = 1350 footprints per scene. AMSU-A scans three times as slowly as AIRS 
(once per 8 seconds) and its footprints are approximately three times as large as
those of AIRS (45 km at nadir). This results in three AIRS scans per AMSU-A scans
and nine AIRS footprints per AMSU-A footprint.

For more details on the history of the AMSU/A satellite instrumentation
see the following `link <https://en.wikipedia.org/wiki/Advanced_microwave_sounding_unit#History>`_.

To summarize, AMSU/A was flown on satellites NOAA 15-17. Versions of AMSU-A also
fly on the Aqua satellite (that also houses AIRS) as well as the European MetOp
satellite. It has been replaced by the Advance Technology Microwave Sounder (ATMS)
on the satellite NOAA-20.

Instructions to download the AMSU-A L1B Version 5 (AIRABRAD) dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The AMSU-A data is located within the Goddard Earth Sciences Data and Information
Services Center (GES DISC) `located here <https://disc.gsfc.nasa.gov/>`_. You need
to create an Earthdata account before you can download data. To access the 
AMSU-A data, search on keyword ``AIRABRAD`` and locate
the **AIRS/Aqua L1B AMSU (A1/A2) geolocated and calibrated brightness temperatures V005
(AIRABRAD)** heading within your search results. 

Next, under the Data Access header, click on `Subset/Get Data`, then refine your
search results by 1) data range (time) and 2) spatial region.

There are various options for downloading, however, the most straightforward approach
for macOS and Linux users is to use the ``wget`` command.  The ``download instructions``
provide the proper wget flags/options.  The ``Download Links List`` provides
the AMSU-A file list based on your search results.


| Each granule is about 560K and has names like

::

   AIRS.2019.06.22.236.L1B.AMSU_Rad.v5.0.0.0.G19174110442.hdf

Advanced Microwave Sounding Unit (AMSU-A) L1B Brightness Temperatures
---------------------------------------------------------------------

Perform the following steps to convert the AMSU_L1 observations:

1. Download the `h4tonccf_nc4 tool <http://hdfeos.org/software/h4cflib.php>`_ provided 
   from the hdf-eos website. Options are provided for Mac, Linux and Windows platforms. 
   For example, the following command downloads the CentOS7 v1.3 executable that
   works for Derecho:
   ::
   
     wget  https://hdfeos.org/software/h4cflib/bin/linux/v1.3/CentOS7/h4tonccf_nc4

2. Convert the AMSU data file from HDF-EOS to netCDF format using the ``h4tonccf_nc4``
   exectuable as shown below. Be sure to provide execute permission first:
   ::
   
      chmod +x h4tonccf_nc4
      ./h4tonccf_nc4 AMSU.hdf
   
      Done with writing netcdf file AMSU.nc

2. b. Optional: The netCDF files have two global attributes that are exceedingly large and uninformative. If needed you can remove these attributes, you can use the 
   ``ncatted`` command from
   `NCO <http://nco.sourceforge.net/nco.html>`_ through the following command:
   ::

       module load nco
       ncatted -a coremetadata,global,d,,, -a StructMetadata_0,global,d,,, AMSU.nc AMSU_final.nc

3. Run ``convert_amsu_L1`` to convert the AMSU_final.nc file to the DART obs_seq format.
   **Important:** Be sure to configure your namelist settings (below) before running the 
   converter.  Also be sure you have compiled the ``convert_amsu_L1`` executable using
   the proper ~/DART/build_templates/mkmf.template that includes both RTTOV and HDF-EOS2
   libraries as described here: :doc:`./README`  
       
   ::
 
   ./convert_amsu_L1


Check the completed ``obs_seq``. It should include brightness temperatures for
the ``EOS_2_AMSUA_TB`` observation type.  The converter should also produce the
following metadata underneath the ``mw`` (microwave) header as shown in the table
below. For more information on the metadata see the
`RTTOV documentation <https://www.nwpsaf.eu/site/software/rttov/documentation/>`_

.. container::

   +-----------------------+------------------------+
   | Metadata variable Name| Description            | 
   +=======================+========================+
   | Sat_az                | Azimuth of satellite   |
   |                       | position (degrees)     |
   +-----------------------+------------------------+
   | Sat_ze                | Aenith of satellite    |
   |                       | position (degrees)     |
   +-----------------------+------------------------+
   | Platform_id           | EOS (9), RTTOV User    | 
   |                       | Guide, Table 2         |
   +-----------------------+------------------------+
   | Sat_id                | (2), RTTOV User        | 
   |                       | Guide, Table 2         | 
   +-----------------------+------------------------+       
   | Sensor_id             | AMSU-A (3), RTTOV User |                        
   |                       | Guide, Table 2         | 
   +-----------------------+------------------------+
   | Channel               | Microwave frequency    |
   |                       | channel (1-15)         | 
   +-----------------------+------------------------+
   | Mag_field             | Earth magnetic field   | 
   |                       | strength (Gauss)       | 
   +-----------------------+------------------------+
   | cosbk                 | Cosine of angle between|                 
   |                       | magnetic field and     | 
   |                       | viewing direction      |
   +-----------------------+------------------------+
   | Fastem_p(1-5)         | Land/sea-ice parameters|                                        
   |                       | 1-5 for FASTEM         | 
   |                       | emissivity model       |
   |                       | Table 21, RTTOV User   |
   |                       | Guide                  |
   +-----------------------+------------------------+



Namelist
~~~~~~~~

The ``convert_amsu_L1`` converter requires :ref:`obs_def_rttov_mod`.
Only two ``&obs_def_rttov_nml`` options are required when converting
the observations: ``use_zeeman`` and ``rttov_sensor_db_file``.

Be aware that if the RTTOV namelist option ``use_zeeman = .true.``
certain metadata must be available in the observation. This is not fully
implemented in the AMSU-A observation converter, so we recommend setting
``use_zeeman = .false.``. For more information,
please see GitHub Issue 99 “`AIRS AMSUA observation converter … Zeeman
coefficients and channels <https://github.com/NCAR/DART/issues/99>`__”

Namelists are read in a file called ``input.nml``. We adhere to the F90 
standard of starting a namelist with an ampersand '&' and terminating with a 
slash '/' for all our namelist input. Character strings that contain a '/' must be
enclosed in quotes to prevent them from prematurely terminating the namelist.
The default values are shown below. More realistic values are provided in
``AIRS/work/input.nml``

::

   &convert_amsu_L1_nml
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

::

  &obs_def_rttov_nml
   rttov_sensor_db_file   = '../../../forward_operators/rttov_sensor_db.csv'
   use_zeeman             = .false.
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
   | along_track_thin   | integer                | Provides ability to thin the data by keeping every Nth data  |
   |                    |                        | value in the along-track scan.   [0,45]                      |
   |                    |                        | e.g. 4 == keep only every 4th row. 0 is no thinning.         |
   +--------------------+------------------------+--------------------------------------------------------------+
   | cross_track_thin   | integer                | Provides ability to thin the data by keeping every Nth data  |
   |                    |                        | value in the cross-track scan.   [0,30]                      |
   |                    |                        | e.g. 3 == keep every third value. 0 is no thinning.          |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lon1               | real(r8)               | The West-most longitude of interest in degrees. [0.0, 360]   |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lon2               | real(r8)               | The East-most longitude of interest in degrees. [0.0, 360]   |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lat1               | real(r8)               | The South-most latitude of interest in degrees. [-90.0,90.0] |
   +--------------------+------------------------+--------------------------------------------------------------+
   | lat2               | real(r8)               | The North-most latitude of interest in degrees. [-90.0,90.0] |
   +--------------------+------------------------+--------------------------------------------------------------+
   | verbose            | integer                | Controls the amount of run-time output.                      |
   |                    |                        | 0 == bare minimum. 3 is very verbose.                        |
   |                    |                        | Only use 3 if converting one or two files for testing.       |
   +--------------------+------------------------+--------------------------------------------------------------+


Channel Specification
~~~~~~~~~~~~~~~~~~~~~

The following channel description is excerpted from the 
Documentation->README Document `found here <https://disc.gsfc.nasa.gov/datasets/AIRABRAD_005/summary>`_.


   "AMSU-A primarily provides temperature soundings. It is a 15-channel microwave
   temperature sounder implemented as two independently operated modules. Module 1
   (AMSU-A1) has 12 channels in the 50-58 GHz oxygen absorption band which provide
   the primary temperature sounding capabilities and 1 channel at 89 GHz which provides
   surface and moisture information. Module 2 (AMSU-A2) has 2 channels: one at 23.8
   GHz and one at 31.4 GHz which provide surface and moisture information (total
   precipitable water and cloud liquid water)."


To facilitate the selection of channels, either the ``Integer`` or ``String`` values
may be used to specify ``channel_list`` within ``&convert_amsu_L1_nml``. The 
`Documentation` and `netCDF` values are provided for reference only.

For example the following ``channel list`` settings are identical and
specify the AMSU channels centered on 50.3 and 89 GHz:

::

 channel_list       = 3,15
 channel_list       = 'A1-1','A1-13'

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




