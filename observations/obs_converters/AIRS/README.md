
This directory contains the code to convert the L1B AMSU-A Brightness Temperatures
in HDF-EOS format to the DART observation sequence file format.

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

Title: AIRS/Aqua L1B AMSU (A1/A2) geolocated and calibrated brightness temperatures V005
Version: 005
Creator: AIRS project
Publisher: Goddard Earth Sciences Data and Information Services Center (GES DISC)
Release Date: 2007-07-26T00:00:00.000Z
Linkage: https://disc.gsfc.nasa.gov/datacollection/AIRABRAD_005.html

https://docserver.gesdisc.eosdis.nasa.gov/repository/Mission/AIRS/3.3_ScienceDataProductDocumentation/3.3.4_ProductGenerationAlgorithms/README.AIRABRAD.pdf

## Instructions to download the AIRABRAD dataset for the AMSUA converter

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

Each granule is about 28M and has names like:  
`1C.GPM.GMI.XCAL2016-C.20160621-S001235-E014508.013137.V05A.HDF5`


TIM: see the GMI readme for more content - this file is not done!


https://wiki.earthdata.nasa.gov/display/DAS/HDF-EOS2+to+HDF-EOS5+Conversion+Tool
https://wiki.earthdata.nasa.gov/display/DAS/Toolkit+Downloads
I downloaded 
* hdf-4.2.13.tar.gz
* HDF-EOS2.20v1.00.tar.Z
* HDF-EOS2.20v1.00_TestDriver.tar.Z
* HDF-EOS_REF.pdf
* HDF-EOS_UG.pdf
* jpegsrc.v9b.tar.gz
* zlib-1.2.11.tar.gz


The HDF Support portal is https://portal.hdfgroup.org/display/support
and has links for documentation, downloads, HDF-EOS, etc...
There is a download for the 'H4H5Tools 2.2.5' that I will be using to
create a standalone converter so that the AIRS c



https://hdfeos.org/software/library.php#HDF-EOS2

Compatibility Library download location:
https://opensource.gsfc.nasa.gov/projects/HDF-EOS2/index.php
he2he5_lib.zip
he2to5.zip

https://wiki.earthdata.nasa.gov/display/DAS/Toolkit+Downloads  has a link to hdf-eos v2.20

