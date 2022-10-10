MADIS Data Ingest System
========================

Overview
--------

The `MADIS <http://madis.noaa.gov/>`__ (Meteorological Assimilation Data Ingest System) service provides access to
real-time and archived data of a variety of types, with added Quality Control (QC) and integration of data from a
variety of sources.

To convert a series of MADIS data files (where different types of observations are distributed in separate files), one
high level view of the workflow is:

#. convert each madis file, by platform type, into an obs_seq file. one file in, one file out. no time changes. use the
   ``shell_scripts/madis_conv.csh`` script. there are script options for hourly output files, or a single daily output
   file.
#. if you aren't using the wrf preprocessing program, you're ready to go.
#. if you do want to do subsequent wrf preprocessing, you need to:

   #. decide on the windowing. each platform has a different convention and if you're going to put them into the wrf
      preprocessing you'll need to have the windowing match. use the ``shell_scripts/windowing.csh`` script.
   #. the wrf preprocessing takes a list of files and assumes they will all be assimilated at the same time, for
      superob'ing purposes, so it should match the expected assimilation window when running filter.

Data sources
------------

`http://madis.noaa.gov <http://madis.noaa.gov/>`__

There are two satellite wind converter programs; the one in this directory and one in the :doc:`../SSEC/SSEC` directory.
The observations distributed here come from `NESDIS <http://www.nesdis.noaa.gov>`__. The SSEC observations are processed
by SSEC itself and will differ from the observations converted here.

Programs
--------

The programs in the ``DART/observations/MADIS/`` directory extract data from the distribution files and create DART
observation sequence (obs_seq) files. Build them in the ``work`` directory by running the ``./quickbuild.sh`` script.
In addition to the converters, the ``advance_time`` and ``obs_sequence_tool`` utilities will be built.

There are currently converters for these data types:

=========================== ======================
ACARS aircraft T,U,V,Q data convert_madis_acars
Marine surface data         convert_madis_marine
Mesonet surface data        convert_madis_mesonet
Metar data                  convert_madis_metar
Wind Profiler data          convert_madis_profiler
Rawinsonde/Radiosonde data  convert_madis_rawin
Satellite Wind data         convert_madis_satwnd
=========================== ======================

Example data files are in the ``data`` directory. Example scripts for converting batches of these files are in the
``shell_scripts`` directory. These are *NOT* intended to be turnkey scripts; they will certainly need to be customized
for your use. There are comments at the top of the scripts saying what options they include, and should be commented
enough to indicate where changes will be likely to need to be made.

Several converters have compile-time choices for outputting various types of moist variables. Check the source code for
more details. Some converters also read multiple T/F strings from the console (standard input) to control at run-time
what types of observations to convert. Again, check the source code for more details.

Each converter has hard-coded input and output filenames:

======================= ================= ================
convert_madis_acars:    acars_input.nc    obs_seq.acars
convert_madis_marine:   marine_input.nc   obs_seq.marine
convert_madis_mesonet:  mesonet_input.nc  obs_seq.mesonet
convert_madis_metar:    metar_input.nc    obs_seq.metar
convert_madis_profiler: profiler_input.nc obs_seq.profiler
convert_madis_rawin:    rawin_input.nc    obs_seq.rawin
convert_madis_satwnd:   satwnd_input.nc   obs_seq.satwnd
======================= ================= ================

The expected usage pattern is that a script will copy, rename, or make a symbolic link from the actual input file (which
often contains a timestamp in the name) to the fixed input name before conversion, and move the output file to an
appropriate filename before the next invocation of the converter. If an existing observation sequence file of the same
output name is found when the converter is run again, it will open that file and append the next set of observations to
it.
