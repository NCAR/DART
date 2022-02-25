Program ``convert_amsu_L1``
===========================

.. caution:: 

   Before you begin: Installing the libraries needed to read these files can be
   fairly troublesome. The NASA Earthdata Data Access Services website is the
   `download site <https://wiki.earthdata.nasa.gov/display/DAS/Toolkit+Downloads>`__
   for the necessary libraries. An example build script (`AIRS/Build_HDF-EOS.sh`)
   is intended to provide some guidance.

Overview
--------

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

convert_amsua_L1.f90
--------------------

``convert_amsua_L1`` converts the L1B AMSU-A Brightness
Temperatures in netCDF format to the DART observation sequence file format.
The native HDF-EOS2 format files must be converted to netCDF.
The conversion from HDF-EOS2 to netCDF is easily performed by the 
`h4tonccf_nc4 <http://hdfeos.org/software/h4cflib.php>`__ converter.

As you can imagine, you need to download each satellite’s data in a
different way. Also, just for your information, AMSU/B has been replaced
on newer satellites by MHS and HSB, but especially MHS is almost
identical.

Namelist
~~~~~~~~

DARTs design structure has the support for radiance observations (like brightness
temperatures) provided by the :doc:`../../forward_operators/obs_def_rttov_mod`
which depends on HDF5 libraries. Consequently, the ``obs_def_rttov_mod_nml`` namelist 
must appear in the ``input.nml``. However, only two options are used when converting
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


Known Bugs
~~~~~~~~~~

None.


Future Plans
~~~~~~~~~~~~

None.


----------


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


Build
^^^^^^

See the :doc:`Dependencies Section<./README>` of the AIRS/README.

Because the data are distributed in HDF-EOS format, and the RTTOV
libraries require HDF5 (incompatible with HDF-EOS) a two-step conversion
is necessary. The data must be converted from HDF to netCDF (which can
be done without HDF5) and then the netCDF files can be converted to DART
radiance observation format - which is the part that requires
``obs_def_rttov_mod.f90``, which is the part that requires HDF5.

The NASA Earthdata Data Access Services website is the `download
site <https://wiki.earthdata.nasa.gov/display/DAS/Toolkit+Downloads>`__,
at press time, the following packages were required to build HDF-EOS
Release v2.20:

-  hdf-4.2.13.tar.gz
-  HDF-EOS2.20v1.00.tar.Z
-  HDF-EOS2.20v1.00_TestDriver.tar.Z
-  HDF-EOS_REF.pdf
-  HDF-EOS_UG.pdf
-  jpegsrc.v9b.tar.gz
-  zlib-1.2.11.tar.gz

Similarly for HDF-EOS5 Release v5.1.16:

-  HDF-EOS5.1.16.tar.Z
-  HDF-EOS5.1.16_TESTDRIVERS.tar.Z
-  HDF-EOS5_REF.pdf
-  HDF-EOS5_UG.pdf
-  hdf5-1.8.19.tar.gz
-  szip-2.1.1.tar.gz

DART provides a script ``DART/observations/obs_converters/AIRS/BUILD_HDF-EOS.sh`` 
that may help provide support for these libraries. You *will* have to modify it for your
system, and you *probably will* have to iterate on that process. The
script takes the stance that if you have to build HDF4, HDF-EOS, HDF5 …
you might as well build HDF-EOS5 too. The HDF-EOS5 is entirely optional.
The HDF5 will be needed by RTTOV.

Converting from HDF4 to netCDF
------------------------------

There are multiple ways to convert from HDF4 to netCDF. The HDF-EOS
Tools and Information Center provides binaries for several common
platforms as well as source code should you need to build your own.

HDF4 CF CONVERSION TOOLKIT
~~~~~~~~~~~~~~~~~~~~~~~~~~

The HDF-EOS Tools and Information Center provides the `HDF4 CF
CONVERSION TOOLKIT <http://hdfeos.org/software/h4cflib.php>`__

   The HDF4 CF (H4CF) Conversion Toolkit can access various NASA HDF4
   external and HDF-EOS2 external files by following the CF conventions
   external. The toolkit includes a conversion library for application
   developers and a conversion utility for NetCDF users. We have
   translated the information obtained from various NASA HDF-EOS2 and
   HDF4 files and the corresponding product documents into the
   information required by CF into the conversion library. We also have
   implemented an HDF4-to-NetCDF (either NetCDF-3 or NetCDF-4 classic)
   conversion tool by using this conversion library. In this web page,
   we will first introduce how to build the conversion library and the
   tool from the source. Then, we will provide basic usage of the tool
   and the conversion library APIs. The information for the supported
   NASA HDF-EOS2 and HDF4 products and visualization screenshots of some
   converted NetCDF files will also be presented.

If you download a binary, it’s a good habit to verify the checksum.
The download page has a link
to a .pdf that has the known checksums. 
`Here’s how to generate the checksum <https://security.stackexchange.com/questions/189000/how-to-verify-the-checksum-of-a-downloaded-file-pgp-sha-etc>`__.
Be aware that when I downloaded the file (via Chrome or ‘wget’) on an
OSX system, the checksum did not match. When I downloaded the file on a
linux system, the checksum *did* match.

If you download the source, the tar file comes with a ``README`` and an 
``INSTALL``. Please become familiar with them. DART also has a build script:
``AIRS/shell_scripts/Build_HDF_to_netCDF.csh`` that you can customize
after you read the ``INSTALL`` document.

Actually converting to netCDF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While the converter creates very nice netCDF files, there are two global
attributes that are exceedingly large and uninformative. Should you want
to remove them, I suggest using the ``ncatted`` command from
`NCO <http://nco.sourceforge.net/nco.html>`__.

::

   h4tonccf_nc4 AIRS.2019.06.22.236.L1B.AMSU_Rad.v5.0.0.0.G19174110442.hdf bob.nc
   ncatted -a coremetadata,global,d,,, -a StructMetadata_0,global,d,,, bob.nc bill.nc

The DART ``L1_AMSUA_to_netcdf.f90`` program
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before I became aware of ``h4tonccf_nc4``, I was in the process of
writing my own converter ``L1_AMSUA_to_netcdf.f90``. *It is not
finished.* Furthermore, at this stage, I don’t know which variables are
needed to be a viable DART observation sequence file, and I don’t see
the point in converting EVERYTHING.
