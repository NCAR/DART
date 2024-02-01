AIRS and AMSU
=============

.. caution:: 

   Before you begin: Installing the libraries needed to read these files can be
   fairly troublesome. The NASA Earthdata Data Access Services website is the
   `download site <https://wiki.earthdata.nasa.gov/display/DAS/Toolkit+Downloads>`__
   for the necessary libraries. An example build script (`AIRS/Build_HDF-EOS.sh`)
   is intended to provide some guidance.

This directory covers two observation converters:

- :doc:`./convert_airs_L2` for temperature and moisture retrievals.

- :doc:`./convert_amsu_L1` for radiances.

Both converters are in the AIRS directory because of the complicated history
of the data used to create the AIRS L2 product (which includes some AMSU observations).
Since both datasets are HDF - it was believed that some of the routines could be
used by both converters. Alas, that has not proven to be the case.

Atmospheric Infrared Sounder (AIRS) Level 2 observations 
--------------------------------------------------------

The `AIRS <http://airs.jpl.nasa.gov/>`__ instrument is an Atmospheric
Infrared Sounder flying on the `Aqua <http://aqua.nasa.gov>`__
spacecraft. Aqua is one of a group of satellites flying close together
in a polar orbit, collectively known as the “A-train”. The programs in
this directory help to extract the data from the distribution files and
put them into DART observation sequence (obs_seq) file format.

AIRS data includes atmospheric temperature in the troposphere, derived
moisture profiles, land and ocean surface temperatures, surface
emissivity, cloud fraction, cloud top height, and ozone burden in the
atmosphere.


Advanced Microwave Sounding Unit (AMSU-A) L1B Brightness Temperatures 
---------------------------------------------------------------------

The *DART/observations/obs_converters/AIRS* directory contains the code 
to convert the L1B AMSU-A Brightness Temperatures in HDF-EOS2 format to 
the DART observation sequence file format.

There is a little bit of confusing history to be aware of for AMSU/A:

https://en.wikipedia.org/wiki/Advanced_microwave_sounding_unit#History

AMSU/A was flown on NOAA 15-17. It is also on the Aqua satellite (that
also houses AIRS) as well as the European MetOp. It has been replaced by
ATMS on NOAA-20.

Dependencies
------------

Both *convert_airs_L2* and *convert_amsu_L1* require the HDF-EOS libraries.
*convert_amsu_L1* also requires HDF5 support because of
the RTTOV libraries. HDF5 is incompatible with HDF-EOS, so a two-step 
conversion is necessary for the AMSU observations. 
The data must be converted from HDF to netCDF 
(which can be done without HDF5) and then the netCDF files can be 
converted to DART radiance observation format - which requires
``obs_def_rttov_mod.f90``, which depends on HDF5.  To simplify things,
An example build script (*DART/observations/obs_converters/AIRS/Build_HDF-EOS.sh*)
is supplied and may provide some guidance on downloading and building
the libraries required by NASA.

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

*BUILD_HDF-EOS.sh* may help you build these libraries. 
You *will* have to modify it for your
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
^^^^^^^^^^^^^^^^^^^^^^^^^^

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

If you download a binary, it’s a good habit to verify the checksum. The download page has a link
to a .pdf that has the known checksums. `Here’s how to generate the
checksum <https://security.stackexchange.com/questions/189000/how-to-verify-the-checksum-of-a-downloaded-file-pgp-sha-etc>`__.
Be aware that when I downloaded the file (via Chrome or ‘wget’) on an
OSX system, the checksum did not match. When I downloaded the file on a
linux system, the checksum *did* match.

If you download the source, the tar file comes with a ``README`` and an ``INSTALL``. Please become
familiar with them. DART also has a build script:
``AIRS/shell_scripts/Build_HDF_to_netCDF.csh`` that you can customize
after you read the ``INSTALL`` document.

