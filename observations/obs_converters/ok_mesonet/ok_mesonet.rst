Oklahoma Mesonet MDF Data
=========================

Overview
--------

| Program to convert Oklahoma Mesonet MDF files into DART observation sequence files.

Data sources
------------

The observation files can be obtained from the Oklahoma Mesonet archive using urls of the format:
``http://www.mesonet.org/index.php/dataMdfMts/dataController/getFile/YYYYMMDDHHMM/mdf/TEXT``
where YYYYMMDDHHMM is the date
and time of the desired set of observations. Files are available every 5 minutes.

If you are located outside of Oklahoma or are going to use this for a non-research purpose see this web page for
information about access: http://www.mesonet.org/index.php/site/about/data_access_and_pricing

Static fields are drawn from the station description file provided by the OK Mesonet. Update the local file from:
http://www.mesonet.org/index.php/api/siteinfo/from_all_active_with_geo_fields/format/csv

Programs
--------

The programs in the ``DART/observations/ok_mesonet/`` directory extract data from the distribution files and create DART
observation sequence (obs_seq) files. Build them in the ``work`` directory by running the ``./quickbuild.sh`` script.
In addition to the converters, the ``advance_time`` and ``obs_sequence_tool`` utilities will be built.

The converter is a preliminary version which has no namelist inputs. It has hard-coded input and output filenames. It
always reads a data file named ``okmeso_mdf.in`` and creates an output file named ``obs_seq.okmeso``. The converter also
requires a text file with the location of all the observating stations, called ``geoinfo.csv``.

The converter creates observations of the following types:

-  LAND_SFC_ALTIMETER
-  LAND_SFC_U_WIND_COMPONENT
-  LAND_SFC_V_WIND_COMPONENT
-  LAND_SFC_TEMPERATURE
-  LAND_SFC_SPECIFIC_HUMIDITY
-  LAND_SFC_DEWPOINT
-  LAND_SFC_RELATIVE_HUMIDITY

Example data files are in the ``data`` directory. Example scripts for converting batches of these files are in the
``shell_scripts`` directory. These are *NOT* intended to be turnkey scripts; they will certainly need to be customized
for your use. There are comments at the top of the scripts saying what options they include, and should be commented
enough to indicate where changes will be likely to need to be made.

The expected usage pattern is that a script will copy, rename, or make a symbolic link from the actual input file (which
often contains a timestamp in the name) to the fixed input name before conversion, and move the output file to an
appropriate filename before the next invocation of the converter. If an existing observation sequence file of the same
output name is found when the converter is run again, it will open that file and append the next set of observations to
it.
