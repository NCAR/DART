
# Before you begin:

Installing the libraries needed to read these files can be fairly troublesome.
The NASA Earthdata Data Access Services website is the 
[download site](https://wiki.earthdata.nasa.gov/display/DAS/Toolkit+Downloads)
for the necessary libraries. 


 - [AIRS Level 2](#AIRS)
 - [AMSU-A](#AMSUA)
 - [Building the required libraries](#Build)


<a name="AIRS"></a>
# Atmospheric Infrared Sounder (AIRS)

Please open and read the `AIRS.html` document once you have the file local.

The [AIRS](http://airs.jpl.nasa.gov/)
instrument is an Atmospheric Infrared Sounder flying on the
[Aqua](http://aqua.nasa.gov) spacecraft.
Aqua is one of a group of satellites flying close together
in a polar orbit, collectively known as the "A-train".
The programs in this directory help to extract the
data from the distribution files and put them into
DART observation sequence (obs_seq) file format.

AIRS data includes atmospheric temperature in the troposphere,
derived moisture profiles, land and ocean surface temperatures,
surface emmissivity, cloud fraction, cloud top height,
and ozone burden in the atmosphere.

## AIRS DATA SOURCES

Access to the web pages where the AIRS data are stored is available by
[registering](https://airs.jpl.nasa.gov/data/registration)
as a data user.

There are two products this converter can be used on: AIRX2RET, which is the
L2 standard retrieval product using AIRS IR and AMSU (without-HSB); and
AIRS2RET, which is the L2 standard retrieval product using AIRS IR-only.

More detailed information on the
[AIRS2RET data product](http://disc.sci.gsfc.nasa.gov/AIRS/data-holdings/by-data-product-v5/airsL2_Std_AIRS_only.shtml)
and the 
[AIRX2RET data product](http://disc.sci.gsfc.nasa.gov/AIRS/data-holdings/by-data-product/airsL2_Std.shtml)
is available from the nasa web pages.

The data is distributed in
[HDF-4](http://www.hdfgroup.org) format, using
some additional conventions for metadata called
[HDF-EOS](http://hdfeos.org/software.php).
There is a basic library for accessing data in hdf files, and a variety of
[generic tools](http://www.hdfgroup.org/products/index.html)
that work with hdf files.
The specific library we use is the
[HDF-EOS2](http://hdfeos.org/software/library.php#HDF-EOS2)
The web page has a link to specific build
instructions.  Also, see [Build](#build) on this page
for very specific instructions for getting the required software
and building it.  If you find more recent instructions online, use those.
But in the absence of anything else, it's someplace to start.

Besides the programs in this directory, a variety of
[specific tools](http://disc.sci.gsfc.nasa.gov/AIRS/tools.shtml)
targeted at AIRS data are available to help read and browse the data.
General information on using hdf in the earth sciences is available
[here](http://eosweb.larc.nasa.gov/HBDOCS/hdf.html).

Specific information on the AIRS support in DART is further explained in `AIRS.html`.


<a name="AMSU-A L1B"></a>
# L1B AMSU-A Brightness Temperatures

This directory contains the code to convert the L1B AMSU-A Brightness Temperatures
in HDF-EOS2 format to the DART observation sequence file format.

There is a little bit of confusing history to be aware of for AMSU/A:

https://en.wikipedia.org/wiki/Advanced_microwave_sounding_unit#History

AMSU/A was flown on NOAA 15-17. It is also on the Aqua satellite
(that also houses AIRS) as well as the European MetOp.
It has been replaced by ATMS on NOAA-20.

As you can imagine, you need to download each satellite's data in a different way.
Also, just for your information, AMSU/B has been replaced on newer satellites
by MHS and HSB, but especially MHS is almost identical.

Be aware that if the RTTOV namelist option `use_zeeman = .true.` certain metadata must
be available in the observation. This is not fully implemented in the AMSU-A
observation converter. For more information, please see GitHub Issue 99 
"[AIRS AMSUA observation converter ... Zeeman coefficients and channels](https://github.com/NCAR/DART/issues/99)"

## Instructions to download the AIRABRAD dataset

The datset of interest is: "AIRS/Aqua L1B AMSU (A1/A2) geolocated and 
calibrated brightness temperatures V005 (AIRABRAD) at GES DISC"
The _short name_ for this dataset is 'AIRABRAD'

The introductory paragraph for the dataset is:

> Version 5 is the current version of the data set.tmospheric Infrared Sounder (AIRS)
> is a grating spectrometer (R = 1200) aboard the second Earth Observing System (EOS)
> polar-orbiting platform, EOS Aqua. In combination with the Advanced Microwave
> Sounding Unit (AMSU) and the Humidity Sounder for Brazil (HSB), AIRS constitutes
> an innovative atmospheric sounding group of visible, infrared, and microwave
> sensors. The AMSU-A instrument is co-aligned with AIRS so that successive blocks
> of 3 x 3 AIRS footprints are contained within one AMSU-A footprint. AMSU-A is
> primarily a temperature sounder that provides atmospheric information in the
> presence of clouds, which can be used to correct the AIRS infrared measurements
> for the effects of clouds. This is possible because non-precipitating clouds are
> for the most part transparent to microwave radiation, in contrast to visible and
> infrared radiation which are strongly scattered and absorbed by clouds. AMSU-A1
> has 13 channels from 50 - 90 GHz and AMSU-A2 has 2 channels from 23 - 32 GHz.
> The AIRABRAD_005 products are stored in files (often referred to as "granules")
> that contain 6 minutes of data, 30 footprints across track by 45 lines along track.

The citation information for this dataset is:

> Title: AIRS/Aqua L1B AMSU (A1/A2) geolocated and calibrated brightness temperatures V005
> Version: 005
> Creator: AIRS project
> Publisher: Goddard Earth Sciences Data and Information Services Center (GES DISC)
> Release Date: 2007-07-26T00:00:00.000Z
> Linkage: https://disc.gsfc.nasa.gov/datacollection/AIRABRAD_005.html

[README.AIRABRAD.pdf](https://docserver.gesdisc.eosdis.nasa.gov/repository/Mission/AIRS/3.3_ScienceDataProductDocumentation/3.3.4_ProductGenerationAlgorithms/README.AIRABRAD.pdf)

## Instructions to download the AIRABRAD dataset

1. Go to https://earthdata.nasa.gov
2. Log in (or create an account if necessary) 
2. Search for AIRABRAD 
3. Scroll down past datasets to "Matching results."
  - Follow the link to "AIRS/Aqua L1B AMSU (A1/A2) geolocated and calibrated brightness temperatures V005 (AIRABRAD) at GES DISC"
4. You should now be at 'https://cmr.earthdata.nasa.gov/search/concepts/C1243477366-GES_DISC.html' (unless they've changed the site).
  - Select the 'Download data' tab
  - Select 'Earthdata search'
  - Select the AIRS link under 'Matching datasets' (I have not tested the NRT products)
5. You can now select 'Granule filters' to choose your start and end dates.
6. Select the granules you want, then click 'download all' and 'download data'
7. Click download access script
8. Follow the instructions on that page to download the data. 

Each granule is about 560K and has names like:  
`AIRS.2019.06.22.236.L1B.AMSU_Rad.v5.0.0.0.G19174110442.hdf`

<a name="Build"></a>
# Build

Because the data are distributed in HDF-EOS format, and the RTTOV libraries 
require HDF5 (incompatible with HDF-EOS) a two-step conversion is necessary.
The data must be converted from HDF to netCDF (which can be done without HDF5)
and then the netCDF files can be converted to DART radiance observation format - 
which is the part that requires `obs_def_rttov_mod.f90`, which is the part that
requires HDF5.

The NASA Earthdata Data Access Services website is the 
[download site](https://wiki.earthdata.nasa.gov/display/DAS/Toolkit+Downloads),
at press time, the following packages were required to build HDF-EOS Release v2.20:

 - hdf-4.2.13.tar.gz
 - HDF-EOS2.20v1.00.tar.Z
 - HDF-EOS2.20v1.00_TestDriver.tar.Z
 - HDF-EOS_REF.pdf
 - HDF-EOS_UG.pdf
 - jpegsrc.v9b.tar.gz
 - zlib-1.2.11.tar.gz

Similarly for HDF-EOS5 Release v5.1.16:

 - HDF-EOS5.1.16.tar.Z
 - HDF-EOS5.1.16_TESTDRIVERS.tar.Z
 - HDF-EOS5_REF.pdf
 - HDF-EOS5_UG.pdf
 - hdf5-1.8.19.tar.gz
 - szip-2.1.1.tar.gz

DART provides a script `BUILD_HDF-EOS.sh`
that may help provide support for these libraries. You _will_ have to modify
it for your system, and you _probably will_ have to iterate on that process.
The script takes the stance that if you have to build HDF4, HDF-EOS, HDF5 ...
you might as well build HDF-EOS5 too. The HDF-EOS5 is entirely optional.
The HDF5 will be needed by RTTOV.

## Converting from HDF4 to netCDF

There are multiple ways to convert from HDF4 to netCDF.
The HDF-EOS Tools and Information Center provides binaries for several common 
platforms as well as source code should you need to build your own. 

### HDF4 CF CONVERSION TOOLKIT

The HDF-EOS Tools and Information Center provides the 
[HDF4 CF CONVERSION TOOLKIT](http://hdfeos.org/software/h4cflib.php)

> The HDF4 CF (H4CF) Conversion Toolkit can access various NASA HDF4 external and 
> HDF-EOS2 external files by following the CF conventions external. The toolkit includes 
> a conversion library for application developers and a conversion utility for NetCDF 
> users. We have translated the information obtained from various NASA HDF-EOS2 and 
> HDF4 files and the corresponding product documents into the information required by 
> CF into the conversion library. We also have implemented an HDF4-to-NetCDF (either 
> NetCDF-3 or NetCDF-4 classic) conversion tool by using this conversion library. 
> In this web page, we will first introduce how to build the conversion library and 
> the tool from the source. Then, we will provide basic usage of the tool and the 
> conversion library APIs. The information for the supported NASA HDF-EOS2 and 
> HDF4 products and visualization screenshots of some converted NetCDF files will 
> also be presented.

#### If you download a binary,

it's a good habit to verify the checksum. The download 
page has a link to a .pdf that has the known checksums. 
[Here's how to generate the checksum](https://security.stackexchange.com/questions/189000/how-to-verify-the-checksum-of-a-downloaded-file-pgp-sha-etc).
Be aware that when I downloaded the file (via Chrome or 'wget') on an OSX system, 
the checksum did not match. When I downloaded the file on a linux system, 
the checksum *did* match.

#### If you download the source,

the tar file comes with a `README` and an `INSTALL`. Please become familiar with them.
DART also has a build script: `AIRS/shell_scripts/Build_HDF_to_netCDF.csh` that you 
can customize after you read the `INSTALL` document.

#### Actually converting to netCDF

While the converter creates very nice netCDF files, there are two global attributes 
that are exceedingly large and uninformative. Should you want to remove them, I suggest
using the `ncatted` command from [NCO](http://nco.sourceforge.net/nco.html). 

```
h4tonccf_nc4 AIRS.2019.06.22.236.L1B.AMSU_Rad.v5.0.0.0.G19174110442.hdf bob.nc
ncatted -a coremetadata,global,d,,, -a StructMetadata_0,global,d,,, bob.nc bill.nc
```

### The DART `L1_AMSUA_to_netcdf.f90` program

Before I became aware of `h4tonccf_nc4`, I was in the process of writing my 
own converter `L1_AMSUA_to_netcdf.f90`.  _It is not finished._
Furthermore, at this stage, I don't know which variables are needed to be a viable
DART observation sequence file, and I don't see the point in converting EVERYTHING.

