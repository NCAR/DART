AIRS Observations
=================

Overview
--------

| The `AIRS <http://airs.jpl.nasa.gov/>`__ instrument is an Atmospheric Infrared Sounder flying on the
  `Aqua <http://aqua.nasa.gov/>`__ spacecraft. Aqua is one of a group of satellites flying close together in a polar
  orbit, collectively known as the "A-train". The programs in this directory help to extract the data from the
  distribution files and put them into DART observation sequence (obs_seq) file format.
| AIRS data includes atmospheric temperature in the troposphere, derived moisture profiles, land and ocean surface
  temperatures, surface emmissivity, cloud fraction, cloud top height, and ozone burden in the atmosphere.

Data sources
------------

Access to the web pages where the AIRS data are stored is available by
`registering <https://airs.jpl.nasa.gov/data/registration>`__ as a data user.

There are two products this converter can be used on: AIRX2RET, which is the L2 standard retrieval product using AIRS IR
and AMSU (without-HSB); and AIRS2RET, which is the L2 standard retrieval product using AIRS IR-only. More detailed
information on the `AIRS2RET data
product <http://disc.sci.gsfc.nasa.gov/AIRS/data-holdings/by-data-product-v5/airsL2_Std_AIRS_only.shtml>`__ and the
`AIRX2RET data product <http://disc.sci.gsfc.nasa.gov/AIRS/data-holdings/by-data-product/airsL2_Std.shtml>`__ is
available from the nasa web pages.

The data is distributed in `HDF-4 <http://www.hdfgroup.org/>`__ format, using some additional conventions for metadata
called `HDF-EOS <http://hdfeos.org/software.php>`__. There is a basic library for accessing data in hdf files, and a
variety of `generic tools <http://www.hdfgroup.org/products/index.html>`__ that work with hdf files. The specific
libraries we use are the `HDF-EOS2 <http://hdfeos.org/software/library.php#HDF-EOS2>`__ library built on HDF4. The web
page has a link to specific build instructions. Also, see below on this web page for very specific instructions for
getting the required software and building it. If you find more recent instructions online, use those. But in the
absence of anything else, it's someplace to start.

Besides the programs in this directory, a variety of `specific tools <http://disc.sci.gsfc.nasa.gov/AIRS/tools.shtml>`__
targeted at AIRS data are available to help read and browse the data. General information on using hdf in the earth
sciences is available `here <http://eosweb.larc.nasa.gov/HBDOCS/hdf.html>`__.

Several types of AIRS data, with varying levels of processing, are available. The following descriptions are taken from
the
`V5_Data_Release_UG <http://disc.sci.gsfc.nasa.gov/AIRS/documentation/v5_docs/AIRS_V5_Release_User_Docs/V5_Data_Release_UG.pdf>`__
document:

   The L1B data product includes geolocated, calibrated observed microwave, infrared and visible/near infrared
   radiances, as well as Quality Assessment (QA) data. The radiances are well calibrated; however, not all QA data have
   been validated. Each product granule contains 6 minutes of data. Thus there are 240 granules of each L1B product
   produced every day.

   The L2 data product includes geolocated, calibrated cloud-cleared radiances and 2-dimensional and 3-dimensional
   retrieved physical quantities (e.g., surface properties and temperature, moisture, ozone, carbon monoxide and methane
   profiles throughout the atmosphere). Each product granule contains 6 minutes of data. Thus there are 240 granules of
   each L2 product produced every day.

   The L3 data are created from the L2 data product by binning them in 1°x1° grids. There are three products: daily,
   8-day and monthly. Each product provides separate ascending (daytime) and descending (nighttime) binned data sets.

The converter in this directory processes level 2 (L2) data files, using data set ``AIRS_DP`` and data product
``AIRX2RET`` or ``AIRS2RET`` without ``HSB`` (the instrument measuring humidity which failed).

The Atmospheric Infrared Sounder (AIRS) is a facility instrument aboard the second Earth Observing System (EOS)
polar-orbiting platform, EOS Aqua. In combination with the Advanced Microwave Sounding Unit (AMSU) and the Humidity
Sounder for Brazil (HSB), AIRS constitutes an innovative atmospheric sounding group of visible, infrared, and microwave
sensors. AIRS data will be generated continuously. Global coverage will be obtained twice daily (day and night) on a
1:30pm sun synchronous orbit from a 705-km altitude.

The AIRS Standard Retrieval Product consists of retrieved estimates of cloud and surface properties, plus profiles of
retrieved temperature, water vapor, ozone, carbon monoxide and methane. Estimates of the errors associated with these
quantities will also be part of the Standard Product. The temperature profile vertical resolution is 28 levels total
between 1100 mb and 0.1 mb, while moisture profile is reported at 14 atmospheric layers between 1100 mb and 50 mb. The
horizontal resolution is 50 km. An AIRS granule has been set as 6 minutes of data, 30 footprints cross track by 45 lines
along track. The Shortname for this product is AIRX2RET. (AIRS2RET is the same product but without the AMSU data.)

The converter outputs temperature observations at the corresponding vertical pressure levels. However, the moisture obs
are the mean for the layer, so the location in the vertical is the midpoint, in log space, of the current layer and the
layer above it. There is an alternative computation for the moisture across the layer which may be more accurate, but
requires a forward operator subroutine to be written and for the observation to contain metadata. The observation could
be defined with a layer top, in pressure, and a number of points to use for the integration across the layer. Then the
forward operator would query the model at each of the N points in the vertical for a given horizontal location, and
compute the mean moisture value. This code has not been implemented yet, and would require a different QTY_xxx to
distinguish it from the simple location/value moisture obs. See the GPS non-local operator code for an example of how
this would need to be implemented.

Getting the data currently means putting in a start/stop time at `this web
page <http://mirador.gsfc.nasa.gov/cgi-bin/mirador/homepageAlt.pl?keyword=AIRX2RET>`__. The keyword is ``AIRX2RET`` and
put in the time range of interest and optionally a geographic region. Each file contains 6 minutes of data, is about 2.3
Megabytes, and globally there are 240 files/day (about 550 Megabytes/day). There are additional options for getting only
particular variables of interest, but the current reader expects whole files to be present. Depending on your connection
to the internet, there are various options for downloading. We have chosen to download a ``wget`` script which is
created by the web page after adding the selected files to a 'cart' and 'checking out'. The script has a series of
``wget`` commands which downloads each file, one at a time, which is run on the machine where you want the data to end
up.

Programs
--------

The temperature observations are located on standard levels; there is a single array of heights in each file and all
temperature data is located on one of these levels. The moisture observations, however, are an integrated quantity for
the space between the levels; in their terminology the fixed heights are 'levels' and the space between them are
'layers'. The current converter locates the moisture obs at the midpoint, in log space, between the levels.

The hdf files need to be downloaded from the data server, in any manner you choose. The converter program reads each hdf
granule and outputs a DART obs_seq file containing up to 56700 observations. Only those with a quality control of 0
(Best) are kept. The resulting obs_seq files can be merged with the
:doc:`../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool` into larger time periods.

It is possible to restrict the output observation sequence to contain data from a region of interest throught the use of
the namelist parameters. If you need a region that spans the Prime Meridian lon1 can be a larger number than lon2, for
example, a region from 300 E to 40 E and 60 S to 30 S (some of the South Atlantic), would be *lon1 = 300, lon2 = 40,
lat1 = -60, lat2 = -30*.

The scripts directory here includes some shell scripts that make use of the fact that the AIRS data is also archived on
the NCAR HPSS (tape library) in daily tar files. The script has options to download a day of granule files, convert
them, merge them into daily files, and remove the original data files and repeat the process for any specified time
period. (See ``oneday_down.sh``)

Here is a very specific script I used to build the required libraries on a Linux cluster. If you find more up-to-date
instructions, use those. But in the absence of anything else, here's a place to start:

   ::

       
      wget https://observer.gsfc.nasa.gov/ftp/edhs/hdfeos/latest_release/*

      # NOTE: direct ftp does not work for me anymore

      ##ftp edhs1.gsfc.nasa.gov
      ### (log in as 'anonymous' and your email as the password)
      ##cd /edhs/hdfeos/latest_release
      ##mget *
      ##quit
       
      # mar 2013, the dir contents:
      # 
      # hdf-4.2.6.tar.gz
      # HDF-EOS2.18v1.00.tar.Z
      # HDF-EOS2.18v1.00_TestDriver.tar.Z
      # HDF_EOS_REF.pdf
      # HDF_EOS_UG.pdf
      # jpegsrc.v6b.tar.gz
      # zlib-1.2.5.tar.gz
      # 
      # (i skipped a 'windows' dir).
      # 
      # mar 2019 contents:
      #      HDF-EOS2.20v1.00.tar.Z  08-Jan-2018 15:21  7.3M  
      #      HDF-EOS2.20v1.00_Tes..> 08-Jan-2018 15:21  9.5M  
      #      HDF-EOS_REF.pdf         07-Nov-2018 13:45  695K  
      #      HDF-EOS_UG.pdf          08-Jan-2018 15:28  429K  
      #      hdf-4.2.13.tar.gz       08-Jan-2018 15:14  4.3M  
      #      jpegsrc.v9b.tar.gz      09-Jan-2018 13:44  1.0M  
      #      zlib-1.2.11.tar.gz      08-Jan-2018 15:22  593K  
      #
      for i in *.tar.gz
      do
        tar -zxvf $i
      done

      # 
      # start with smaller libs, work up to HDF-EOS.
      # 
      # 

      echo zlib:

      cd zlib-1.2.11
      ./configure --prefix=/glade/p/work/nancy
      make
      make test 
      make install

      echo jpeg:

      cd jpeg-9b
      ./configure --prefix=/glade/p/work/nancy
      make
      make test 
      mkdir /glade/p/work/nancy/{bin,man,man/man1} 
      make install

      # (make install wouldn't create the dirs if they didn't exist.
      # lib was there from the zlib install, the others weren't.)

      echo hdf:

      cd hdf-4.2.13
      ./configure --prefix=/glade/p/work/nancy
      # (it found zlib and jpeg, from the install prefix i guess)
      make
      # (there is apparently no 'make test')
      make install

      echo hdf-eos:

      cd hdfeos
      ./configure CC='/glade/p/work/nancy/bin/h4cc -Df2cFortran' --prefix=/glade/p/work/nancy
      # (the CC= is crucial)
      make
      # (i didn't build the test drivers so i didn't do make test)
      make install


      echo AIRS converter:

      cd $DART/observations/AIRS/work

      echo edit mkmf_convert_airs_L2 to have all the base paths
      echo be /glade/p/work/nancy instead of whatever.  make it look like:
      echo ' '
      echo 'set JPGDIR = /glade/work/nancy'
      echo 'set HDFDIR = /glade/work/nancy'
      echo 'set EOSDIR = /glade/work/nancy'
      echo ' '

      ./quickbuild.csh

      exit 0

Namelist
--------

This namelist is read in a file called ``input.nml``. We adhere to the F90 standard of starting a namelist with an
ampersand '&' and terminating with a slash '/' for all our namelist input. Character strings that contain a '/' must be
enclosed in quotes to prevent them from prematurely terminating the namelist.

::

   &convert_airs_L2_nml
      l2_files           = 'input.hdf',
      l2_file_list       = '',
      datadir            = '.',
      outputdir          = '.',
      lon1               =   0.0,
      lon2               = 360.0,
      lat1               = -90.0,
      lat2               =  90.0,
      min_MMR_threshold  = 1.0e-30,
      top_pressure_level = 0.0001,
      cross_track_thin   = 0,
      along_track_thin   = 0,
   /

| 

.. container::

   +-------------------+------------------------+------------------------------------------------------------+---------+
   | Contents          | Type                   | Description                                                | Default |
   +===================+========================+============================================================+=========+
   | l2_files          | character(len=128) (:) | A list of one or more names of the HDF file(s) to read,    |         |
   |                   |                        | NOT including the directory. If multiple files are listed, |         |
   |                   |                        | each will be read and the results will be placed in a      |         |
   |                   |                        | separate file with an output filename constructed based on |         |
   |                   |                        | the input filename.                                        |         |
   +-------------------+------------------------+------------------------------------------------------------+---------+
   | l2_file_list      | character(len=128)     | The name of an ascii text file which contains one filename |         |
   |                   |                        | per line, NOT including the directory. Each file will be   |         |
   |                   |                        | read and the observations converted into an output file    |         |
   |                   |                        | where the output filename is based on the input filename.  |         |
   |                   |                        | Only one of 'l2_files' and 'l2_file_list' can be           |         |
   |                   |                        | specified. The other must be ' ' (empty).                  |         |
   +-------------------+------------------------+------------------------------------------------------------+---------+
   | datadir           | character(len=128)     | The directory containing the HDF files                     |         |
   +-------------------+------------------------+------------------------------------------------------------+---------+
   | outputdir         | character(len=128)     | The directory for the output observation sequence files.   |         |
   +-------------------+------------------------+------------------------------------------------------------+---------+
   | lon1              | real(r4)               | the West-most longitude of interest in degrees. [0.0, 360] |         |
   +-------------------+------------------------+------------------------------------------------------------+---------+
   | lon2              | real(r4)               | the East-most longitude of interest in degrees. [0.0, 360] |         |
   +-------------------+------------------------+------------------------------------------------------------+---------+
   | lat1              | real(r4)               | the South-most latitude of interest in degrees. [-90.0,    |         |
   |                   |                        | 90.0]                                                      |         |
   +-------------------+------------------------+------------------------------------------------------------+---------+
   | lat2              | real(r8)               | the North-most latitude of interest in degrees. [-90.0,    |         |
   |                   |                        | 90.0]                                                      |         |
   +-------------------+------------------------+------------------------------------------------------------+---------+
   | min_MMR_threshold | real(r8)               | The data files contains 'Retrieved Water Vapor Mass Mixing |         |
   |                   |                        | Ratio'. This is the minimum threshold, in gm/kg, that will |         |
   |                   |                        | be converted into a specific humidity observation.         |         |
   +-------------------+------------------------+------------------------------------------------------------+---------+
   | cross_track_thin  | integer                | provides ability to thin the data by keeping only every    |         |
   |                   |                        | Nth data value in a particular row. e.g. 3 == keep every   |         |
   |                   |                        | third value.                                               |         |
   +-------------------+------------------------+------------------------------------------------------------+---------+
   | along_track_thin  | integer                | provides ability to thin the data by keeping only every    |         |
   |                   |                        | Nth row. e.g. 4 == keep only every 4th row.                |         |
   +-------------------+------------------------+------------------------------------------------------------+---------+

| 
