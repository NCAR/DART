PROGRAM ``MOD15A2_to_obs``
==========================

MODIS land product subsets (collection 5) to DART observation sequence converter
--------------------------------------------------------------------------------

Overview
~~~~~~~~

This routine is designed to convert the `MODIS Land Product Subsets <http://daac.ornl.gov/MODIS/modis.shtml>`__ data of
Leaf Area Index (**LAI**) and Fraction of Photosynthetically Active Radiation (**FPAR**) 8 day composite [MOD15A2] to a
DART observation sequence file. According to the `MODIS LAI/FPAR Product User's
Guide <https://lpdaac.usgs.gov/sites/default/files/public/modis/docs/MODIS-LAI-FPAR-User-Guide.pdf>`__:

   Leaf area index (LAI; dimensionless) is defined as the one-sided green leaf area per unit ground area in broadleaf
   canopies and as one-half the total needle surface area per unit ground area in coniferous canopies.
   Fraction of Photosynthetically Active Radiation absorbed by vegetation (FPAR; dimensionless) is defined as the
   fraction of incident photosynthetically active radiation (400-700 nm) absorbed by the green elements of a vegetation
   canopy.

Specifically, the composites are comma-separated-values (.csv format) ASCII files where each line is a record. The input
``.csv`` files are directly from the Oak Ridge National Laboratory `DAAC <http://daac.ornl.gov>`__. There are two
streams to download the data formats we support, they differ only in the very first line of the file. One of the formats
has a header record, the other does not. Other than that, the file formats are identical. The format with the header
record is fully described in https://lpdaac.usgs.gov/dataset_discovery/modis. Please remember to cite the data in your
publications, `specific instructions from LP DAAC are available
here. <https://lpdaac.usgs.gov/about/citing_lp_daac_and_data>`__ This is an example:

   Data Citation: Oak Ridge National Laboratory Distributed Active Archive Center (ORNL DAAC). 2012. MODIS subsetted
   land products, Collection 5. Available on-line [http://daac.ornl.gov/MODIS/modis.html] from ORNL DAAC, Oak Ridge,
   Tennessee, U.S.A. Accessed *Month dd, yyyy*.

For more information on *downloading* the data, see DATA SOURCES below. The `MODIS Land Product
Subsets <http://daac.ornl.gov/MODIS/modis.shtml>`__ page indicates that the Collection 5 MODIS Subsets are available
three ways:

#. `Field Site and Flux tower <http://daac.ornl.gov/cgi-bin/MODIS/GR_col5_1/mod_viz.html>`__. Since the files are
   preprocessed, the download is immediate. The current state of the converter supports this format.
#. `Global Tool <http://daac.ornl.gov/cgi-bin/MODIS/GLBVIZ_1_Glb/modis_subset_order_global_col5.pl>`__. This requires
   exact knowledge of the location(s) of interest. Because some processing to fulfill the request is needed, a job is
   scheduled on the DAAC server and an email notification is sent with instuctions on how to retrieve the file(s) of
   interest. The converter **does not** currently support this format, but will soon. Worst case scenario is that you
   make your own header file and add your 'site' to the metadata file described below.
#. `Web Service <https://lpdaac.usgs.gov/tools/lp_daac_web_services>`__. I have not used the Web Service.

The DART workflow is usually:

#. download the MOD15A2 data for the sites and years in question (see DATA SOURCES below)
#. build the DART executables with support for ``MODIS_LEAF_AREA_INDEX`` and ``MODIS_FPAR`` observations. This is done
   by running ``preprocess`` with ``obs_def_land_mod.f90`` in the list of ``input_files`` for ``preprocess_nml`` and
   then building ``MOD15A2_to_obs`` in the usual DART way.
#. provide basic information via the ``input.nml``:``MOD15A2_to_obs_nml`` namelist
#. convert each MODIS data file individually using ``MOD15A2_to_obs``
#. combine all output files for the region and timeframe of interest into one file using
   :doc:`../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool`

For some models (CLM, for example), it is required to reorganize the observation sequence files into a series of files
that contains ONLY the observations for each assimilation. This can be achieved with the 
``DART/observations/obs_converters/MODIS/shell_scripts/makedaily.sh`` script.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &MOD15A2_to_obs_nml
      text_input_file = 'MOD15A2.fn_usbouldr.txt',
      metadata_file   = 'MOD15A2_site_metadata.txt',
      obs_out_file    = 'obs_seq.out',
      maxgoodqc       = 10,
      verbose         = .false.
      /

.. container::

   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | Contents        | Type               | Description                                                                 |
   +=================+====================+=============================================================================+
   | text_input_file | character(len=256) | Name of the MODIS file of comma-separated values. This may be a relative or |
   |                 |                    | absolute filename.                                                          |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | metadata_file   | character(len=256) | Name of the file that contains the location information for the specific    |
   |                 |                    | sites. This may be a relative or absolute filename. If this file does not   |
   |                 |                    | exist, it is **presumed** that the location information is part of the      |
   |                 |                    | 'site' column. If this is not true, the program will fail. For more         |
   |                 |                    | information see the section Presumed Format                                 |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | obs_out_file    | character(len=128) | Name of the output observation sequence file.                               |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | maxgoodqc       | real               | maximum value of any observation quality control flag to pass through to    |
   |                 |                    | the output observation sequence. Keep in mind that ``filter`` has the       |
   |                 |                    | ability to discriminate on the value, so there is really little to be       |
   |                 |                    | gained by rejecting them during the conversion. The QC value is passed      |
   |                 |                    | through in its native value, i.e. it is not converted to play nicely with   |
   |                 |                    | observations that have values 0,1,2,3,4,5 etc.                              |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | verbose         | logical            | Print extra information during the ``MOD15A2_to_obs`` execution.            |
   +-----------------+--------------------+-----------------------------------------------------------------------------+

Data sources
------------

Field site and flux tower
~~~~~~~~~~~~~~~~~~~~~~~~~

| The download site for the 'Field Site and Flux tower' data is
| http://daac.ornl.gov/cgi-bin/MODIS/GR_col5_1/mod_viz.html. Since the files are preprocessed, the download is
  immediate. This method results in files **with** the header record, and requires a small amount of additional work:

-  Download the metadata file containing the locations for the Field Sites
   ftp://daac.ornl.gov/data/modis_ascii_subsets/5_MODIS_SUBSETS_C5_&_FLUXNET.csv
-  I usually convert this to UNIX format with the UNIX utility ``dos2unix`` and rename it to
   ``MOD15A2_site_metadata.txt``

| The data files have names like ``MOD15A2.fn_uswiirpi.txt`` or ``MOD15A2.fn_dehambur.txt`` and have very long lines.
  The first line (i.e. record) of the file is a comma-separated list explaining the file format for all the remaining
  lines/records.
| These files contain records with 49 pixel values where each pixel represents the values for a 1km by 1km voxel. The
  center pixel is the only value converted to a DART observation value.

.. container:: unix

   ::

      MODIS_LAI % head -1 MOD15A2.fn_dehambur.txt
      HDFname,Product,Date,Site,ProcessDate,Band,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49

The format of the ``Site`` in these files is the predominant difference between the files from the download methods. The
``Site`` fields in these files have specified site names that must have a case-sensitive match to a site in the metadata
file specified by ``input.nml``:``metadata_file`` .

Global tool
~~~~~~~~~~~

| **This format is not supported yet.**
| The download site for the 'Global Tool' data is
| http://daac.ornl.gov/cgi-bin/MODIS/GLBVIZ_1_Glb/modis_subset_order_global_col5.pl. Because some processing to fulfill
  the request is needed, a job is scheduled on the DAAC server and an email notification is sent with instuctions on how
  to retrieve the file(s) of interest. **This method requires exact knowledge of the location(s) of interest.**
  ``MOD15A2_to_obs`` presumes prior knowledge of the file format and that the latitude and longitude are coded in the
  site name (which is the default behavior). **Do not change the format of the file.** Please follow the download
  instructions below - **exactly.** These instructions were accurate as of 11 April 2014.

#. go to the DAAC `download site for MODIS global
   data <http://daac.ornl.gov/cgi-bin/MODIS/GLBVIZ_1_Glb/modis_subset_order_global_col5.pl>`__.
#. Select either

   #. "Country" (it helps to FIRST clear out the values from the "lat/lon" boxes)
   #. or a specific latitude and longitude. Be precise. This will specify the center pixel location.

#. click "Continue"
#. Select the "[MOD15A2] Leaf Area Index (LAI) and Fraction of Photsyntetically Active Radiation (FPAR) 8 Day Composite"
   from the pull-down menu.
#. **Important:** Specify 3 **and only 3** kilometers to encompass the center location. This results in the 7 km by 7 km
   resolution required by ``MOD15A2_to_obs``.
#. click "Continue"
#. select the Starting Date and Ending Date from the list. You can convert the entire dataset into one long DART
   observation sequence file and then subset it later if need be.
#. **Important:** Make sure you check the button "Generate GeoTIFF and Reproject to Geographic Lat/long"
#. Supply your REAL email address
#. click "Continue"
#. Review the confirmation page. Make sure the requested resolution and area is correct. You should see something like
   "The Requested Data Area is Approximately 7 Kilometers Wide and 7 Kilometers High"
#. click "Continue"
#. At some point later (perhaps even days), you will get an email with the subject "ORNL DAAC MODIS MOD15A2 order",
   follow the instructions to complete the download.

The resulting ASCII files will have the same format as described below. The 'site name' column for these files is of the
form: ``Lat47.61666667Lon12.58333333Samp7Line7`` which provides the location information otherwise provided by the
``MOD15A2_site_metadata.txt`` file for the predefined sites.

Web service
~~~~~~~~~~~

I have not used the `Web Service <https://lpdaac.usgs.gov/tools/lp_daac_web_services>`__.

Format
------

| The data product "Leaf Area Index - Fraction of Photosynthetically Active Radiation 8-Day L4 Global 1km" (**MOD15A2**)
  is described in https://lpdaac.usgs.gov/products/modis_products_table/mod15a2 (**expand the 'Layers' tab**). The units
  and the QC values are described there. What I have not been able to determine is how to interpret the 'Date' ... if it
  is 2000049 ... It is day 49 of year 2000. Is that the start of the 8 day composite, the middle, the end? If you know
  the answer, please let me know.
| Taken (almost) directly from https://lpdaac.usgs.gov/tools/lp_daac_web_services and modified only slightly with
  examples more appropriate for the LAI/FPAR product.
| The MODIS MOD15A2 products in question are ASCII files of comma-separated values. If the file contains a header
  record/line, all columns are interpreted based on this header column. If the file does not contain a header, the
  following format is REQUIRED.

-  ASCII values are comma delimited
-  Row 1 is the header row (which may not exist for products generated by the Global Tool)
-  Data values start in row 2 if the header row is present.
-  Rows of QC data are interleaved with measurement data as indicated in Column 6.
-  Note that values may contain embedded periods, dashes, and underscores (".,-, \_").

+-------------+---------------------------------------------------------------------------------------------------------------+---------------------------------------------------------+
| **Column**  | **Column Description**                                                                                        | **Example Values**                                      |
+-------------+---------------------------------------------------------------------------------------------------------------+---------------------------------------------------------+
| 1           | Unique row identifier                                                                                         | MOD15A2.A2000049.fn_ruyakuts.005.2006268205917.Fpar_1km |
|             |                                                                                                               | MOD15A2.A2000049.fn_ruyakuts.005.2006268205917.Lai_1km  |
+-------------+---------------------------------------------------------------------------------------------------------------+---------------------------------------------------------+
| 2           | MODIS Land Product Code                                                                                       | MOD15A2                                                 |
+-------------+---------------------------------------------------------------------------------------------------------------+---------------------------------------------------------+
| 3           | MODIS Acquisition Date                                                                                        | A2000049 ( ?this is an 8 day average)                   |
|             | A(YYYYDDD)                                                                                                    | What does 49 indicate? start? middle? end?              |
+-------------+---------------------------------------------------------------------------------------------------------------+---------------------------------------------------------+
| 4           | SiteID                                                                                                        | fn_ustnwalk,                                            |
|             | Each site is assigned a unique ID.                                                                            |                                                         |
|             | To get the Site name information from SiteID,                                                                 | Lat47.61666667Lon12.58333333Samp7Line7                  |
|             | `click here <ftp://daac.ornl.gov/data/modis_ascii_subsets/MODIS_Subset_Sites_Information_Collection5.csv>`__  |                                                         |
+-------------+---------------------------------------------------------------------------------------------------------------+---------------------------------------------------------+
| 5           | MODIS Processing Date (YYYYDDDHHMMSS)                                                                         | 2006269073558                                           |
+-------------+---------------------------------------------------------------------------------------------------------------+---------------------------------------------------------+
| 6           | Product Scientific Data Set (Band):                                                                           | MOD15A2: FparExtra_QC, FparLai_QC,                      |
|             | Indicates type of values to follow.                                                                           | FparStdDev_1km, Fpar_1km,                               |
|             | Specific values vary by Product. Data                                                                         | LaiStdDev_1km, Lai_1km                                  |
|             | quality information are interleaved.                                                                          |                                                         |
+-------------+---------------------------------------------------------------------------------------------------------------+---------------------------------------------------------+
| 7 to N      | Data values of type as specified.                                                                             | QC: 00100001,01100001,01100001, ...                     |
|             | Number of data columns as given in                                                                            | Measurement:                                            |
|             | Column 4. Definition of QC component values vary by Scientific Data Set                                       | 2,2,1,1,1,1,1,0,0,0,1,1,0,0, to N                       |
+-------------+---------------------------------------------------------------------------------------------------------------+---------------------------------------------------------+

QC flags are binary-coded ascii strings e.g., 10011101 bits 5,6,7 (the last three) are decoded as follows:

-  000 ... Main(RT) method used, best result possible (no saturation)
-  001 ... Main(RT) method used with saturation, Good, very usable
-  010 ... Main(RT) method failed due to bad geometry, empirical algorithm used
-  011 ... Main(RT) method failed due to other problems
-  100 ... pixel not produced at all

Consequently, the last three digits are used by DART's data processing logic.

Programs
--------

| The ``MOD15A2_to_obs.f90`` file is the source for the main converter program. Look at the source code where it reads
  the example data file. You will almost certainly need to change the "read" statement to match your data format. The
  example code reads each text line into a character buffer and then reads from that buffer to parse up the data items.
| FIXME Explain the 10% for the obs error for FPAR and question the LAIStddev ...

To compile and test, go into the work subdirectory and run the ``quickbuild.sh`` script to build the converter and a
couple of general purpose utilities. ``advance_time`` helps with calendar and time computations, and the
``obs_sequence_tool`` manipulates DART observation files once they have been created.

To change the observation types, look in the ``DART/obs_def`` directory. If you can find an obs_def_XXX_mod.f90 file
with an appropriate set of observation types, change the 'use' lines in the converter source to include those types.
Then add that filename in the ``input.nml`` namelist file to the &preprocess_nml namelist, the 'input_files' variable.
Multiple files can be listed. Then run quickbuild.sh again. It remakes the table of supported observation types before
trying to recompile the source code.

An example script for converting batches of files is in the ``shell_scripts`` directory. A tiny example data file is in
the ``data`` directory. These are *NOT* intended to be turnkey scripts; they will certainly need to be customized for
your use. There are comments at the top of the script saying what options they include, and should be commented enough
to indicate where changes will be likely to need to be made.

Decisions you might need to make
--------------------------------

See the general discussion in the :doc:`../../../guide/creating-obs-seq-real` page about what options are
available for the things you need to specify. These include setting a time, specifying an expected error, setting a
location, and an observation type.


Future plans
------------

- Support for the data records without the header, as created by
  the Global Tool.
- The work that remains is to get the IGBP landcover code for the site and
  incorporate that into the observation metadata. I *almost* have
  everything I need. Once that happens, the forward observation operator
  can be made to be much more accurate by only using model landunits
  that have the right landcover class.


PROGRAM ``MOD15A2_to_obs``
==========================

MODIS land product subsets (collection 5) to DART observation sequence converter
