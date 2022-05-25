SSEC Data Center
================

Overview
--------

The program in this directory takes satellite wind data from the University of Wisconsin-Madison Space Science and
Engineering Center, and converts it into DART format observation sequence files, for use in assimilating with the DART
filter program.

Data sources
------------

The Space Science and Engineering Center (SSEC) at University of Wisconsin-Madison has an online `data
center <http://www.ssec.wisc.edu/data>`__ with both real-time and archival weather satellite data.

The last 2 day's worth of data is available from ftp://cyclone.ssec.wisc.edu/pub/fnoc.

There is a second satellite wind DART converter in the :doc:`../MADIS/MADIS` directory which converts wind observations
which originate from `NESDIS <http://www.nesdis.noaa.gov>`__. The data from this converter is processed at the SSEC and
the observations will be different from the ones distributed by MADIS.

Programs
--------

| Conversion program ``convert_ssec_satwnd`` converts the ascii data in the input files into a DART observation sequence
  file. Go into the ``work`` directory and run the ``quickbuild.sh`` script to compile the necessary files.
| The program reads standard input for the data time range, which types of observations to convert, and then, if quality
  control information is found in the input file, what type of quality control algorithm to use when deciding whether
  the observation is of good quality or not. See the references below.

References
----------

-  RF method: Velden, C. S., T. L. Olander, and S. Wanzong, 1998: The impact of multispectral GOES-8 wind information on
   Atlantic tropical cyclone track forecasts in 1995. Part I: Dataset methodology, description, and case analysis. Mon.
   Wea. Rev., 126, 1202-1218.
-  QI method: Holmlund, K., 1998: The utilization of statistical properties of satellite-derived atmospheric motion
   vectors to derive quality indicators. Wea. Forecasting, 13, 1093-1104.
-  Comparison of two methods: Holmlund, K., C.S. Velden, and M. Rohn, 2001: Enhanced Automated Quality Control Applied
   to High-Density Satellite-Derived Winds. Mon. Wea. Rev., 129, 517-529.
