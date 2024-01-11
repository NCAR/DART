WOD Observations
================

Overview
--------

The World Ocean Database (WOD) is a collection of data from various sources,
combined into a single format with uniform treatment. WOD is created by the 
National Centers for Environmental Information (NCEI) of the National Oceanic
and Atmospheric Administration (NOAA).

An updated version of the dataset is released approximately every four years.
It was first produced in 1994 and has been released in 1998, 2001, 2005, 2009,
2013 and 2018.

The `WOD website <https://www.ncei.noaa.gov/products/world-ocean-atlas>`__ has
detailed information about the repository, observations, and datasets. The
programs in this directory convert from the packed ASCII files found in the
repository into DART observation sequence (obs_seq) file format.

There are two sets of available files: the raw observations and the
observations binned onto standard levels.

.. note::

   DAReS staff recommend using the datasets on standard levels for
   assimilation. The raw data can be very dense in the vertical and are not
   truly independent observations. The correlation between nearby observations
   leads to too much certainty in the updated values during the assimilation.

Data sources
------------

Use already existing obs_seq files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NSF NCAR staff have prepared datasets already converted to DART's obs_seq file 
format for the World Ocean Database 2013 (WOD13) and the World Ocean Database
2009 (WOD09).

.. Warning::  

   The WOD data is in PSU, the dart observation converter ``wod_to_obs`` converts the observations to MSU.
 
   | PSU = g/kg
   | MSU = PSU/1000 = kg/kg
   
   The WOD observation sequence files availiable from NSF NCAR's RDA are in MSU.

WOD13
~~~~~

The already-converted WOD13 dataset comprises data from 2005-01-01 to
2016-12-31 and was created by *Fred Castruccio*. Thanks Fred! The files are
stored in the following directory on GLADE:

.. code-block::

   /glade/p/cisl/dares/Observations/WOD13

The subdirectories are formatted in ``YYYYMM`` order and contain the
following observation types:

+--------------------------------------+--------------------------------------+
| FLOAT_SALINITY                       | FLOAT_TEMPERATURE                    |
+--------------------------------------+--------------------------------------+
| DRIFTER_SALINITY                     | DRIFTER_TEMPERATURE                  |
+--------------------------------------+--------------------------------------+
| GLIDER_SALINITY                      | GLIDER_TEMPERATURE                   |
+--------------------------------------+--------------------------------------+
| MOORING_SALINITY                     | MOORING_TEMPERATURE                  |
+--------------------------------------+--------------------------------------+
| BOTTLE_SALINITY                      | BOTTLE_TEMPERATURE                   |
+--------------------------------------+--------------------------------------+
| CTD_SALINITY                         | CTD_TEMPERATURE                      |
+--------------------------------------+--------------------------------------+
| XCTD_SALINITY                        | XCTD_TEMPERATURE                     |
+--------------------------------------+--------------------------------------+
| APB_SALINITY                         | APB_TEMPERATURE                      |
+--------------------------------------+--------------------------------------+
| XBT_TEMPERATURE                      |                                      |
+--------------------------------------+--------------------------------------+

If you use WOD13, please cite Boyer et al. (2013). [1]_

WOD09
~~~~~

The already-converted WOD09 dataset, which comprises data from 1960-01-01 to
2008-12-31, is stored in the following directory on GLADE:

.. code-block::

   /glade/p/cisl/dares/Observations/WOD09

If you use WOD09, please cite Johnson et al. (2009). [2]_ 

Download WOD from NCEI
^^^^^^^^^^^^^^^^^^^^^^

Data from each of the WOD releases can be downloaded interactively from the 
`WOD website <https://www.ncei.noaa.gov/products/world-ocean-atlas>`__.

Download WOD from NSF NCAR
^^^^^^^^^^^^^^^^^^^^^^^^^^

WOD09 can also be downloaded from NSF NCAR's `research data archive (RDA) dataset 
285.0 <https://rda.ucar.edu/datasets/ds285.0/>`__.

Programs
--------

The data is distributed in a specialized packed ASCII format. In this directory is a program called ``wodFOR.f`` which
is an example reader program to print out data values from the files. The program ``wod_to_obs`` converts these packed
ASCII files into DART obs_sequence files.

As with most other DART directories, the ``work`` directory contains a ``quickbuild.sh`` script to build all necessary
executables.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &wod_to_obs_nml
      wod_input_file       =  'XBTS2005',
      wod_input_filelist   =  '',
      wod_out_file         =  'obs_seq.wod',
      avg_obs_per_file     =  500000,
      debug                =  .false.,
      timedebug            =  .false.,
      print_qc_summary     =  .true.,
      max_casts            =  -1,
      no_output_file       =  .false.,
      print_every_nth_cast =  -1,
      temperature_error    =  0.5,
      salinity_error       =  0.5, 
    /
   ! temperature error is in degrees C, salinity error in g/kg.

| 

.. container::

   +----------------------+--------------------+------------------------------------------------------------------------+
   | Item                 | Type               | Description                                                            |
   +======================+====================+========================================================================+
   | wod_input_file       | character(len=128) | The input filename when converting a single file. Only one of the two  |
   |                      |                    | namelist items that specify input files can have a valid value, so to  |
   |                      |                    | use a single filename set the list name 'wod_input_filelist' to the    |
   |                      |                    | empty string (' ').                                                    |
   +----------------------+--------------------+------------------------------------------------------------------------+
   | wod_input_filelist   | character(len=128) | To convert one or more files in a single execution create a text file  |
   |                      |                    | which contains each input filename, in ascii, one filename per line.   |
   |                      |                    | Set this item to the name of that file, and set 'wod_input_file' to    |
   |                      |                    | the empty string (' ').                                                |
   +----------------------+--------------------+------------------------------------------------------------------------+
   | wod_out_file         | character(len=128) | The output file to be created. Note that unlike earlier versions of    |
   |                      |                    | some converters, this program will overwrite an existing output file   |
   |                      |                    | instead of appending to it. The risk of replicated observations, which |
   |                      |                    | are difficult to detect since most of the contents are floating point  |
   |                      |                    | numbers, outweighed the possible utility.                              |
   +----------------------+--------------------+------------------------------------------------------------------------+
   | avg_obs_per_file     | integer            | The code needs an upper limit on the number of observations generated  |
   |                      |                    | by this program. It can be larger than the actual number of            |
   |                      |                    | observations converted. The total number of obs is computed by         |
   |                      |                    | multiplying this number by the number of input files. If you get an    |
   |                      |                    | error because there is no more room to add observations to the output  |
   |                      |                    | file, increase this number. Do not make this an unreasonably huge      |
   |                      |                    | number, however, since the code does preallocate space and will be     |
   |                      |                    | slow if the number of obs becomes very large.                          |
   +----------------------+--------------------+------------------------------------------------------------------------+
   | print_every_nth_cast | integer            | If a value greater than 0, the program will print a message after      |
   |                      |                    | processing every N casts. This allows the user to monitor the progress |
   |                      |                    | of the conversion.                                                     |
   +----------------------+--------------------+------------------------------------------------------------------------+
   | print_qc_summary     | logical            | If .TRUE. the program will print out a summary of the number of casts  |
   |                      |                    | which had a non-zero quality control values (current files appear to   |
   |                      |                    | use values of 1-9).                                                    |
   +----------------------+--------------------+------------------------------------------------------------------------+
   | debug                | logical            | If .TRUE. the program will print out debugging information.            |
   +----------------------+--------------------+------------------------------------------------------------------------+
   | timedebug            | logical            | If .TRUE. the program will print out specialized time-related          |
   |                      |                    | debugging information.                                                 |
   +----------------------+--------------------+------------------------------------------------------------------------+
   | max_casts            | integer            | If a value greater than 0 the program will only convert at most this   |
   |                      |                    | number of casts from each input file. Generally only expected to be    |
   |                      |                    | useful for debugging. A negative value will convert all data from the  |
   |                      |                    | input file.                                                            |
   +----------------------+--------------------+------------------------------------------------------------------------+
   | no_output_file       | logical            | If .TRUE. the converter will do all the work needed to convert the     |
   |                      |                    | observations, count the number of each category of QC values, etc, but |
   |                      |                    | will not create the final obs_seq file. Can be useful if checking an   |
   |                      |                    | input file for problems, or for getting QC statistics without waiting  |
   |                      |                    | for a full output file to be constructed, which can be slow for large  |
   |                      |                    | numbers of obs. Only expected to be useful for debugging.              |
   +----------------------+--------------------+------------------------------------------------------------------------+
   | temperature_error    | real(r8)           | The combined expected error of temperature observations from all       |
   |                      |                    | sources, including instrument error, model bias, and                   |
   |                      |                    | representativeness error (e.g. larger or smaller grid box sizes        |
   |                      |                    | affecting expected accuracy), in degrees Centigrade. Values in output  |
   |                      |                    | file are error variance, which will be this value squared.             |
   +----------------------+--------------------+------------------------------------------------------------------------+
   | salinity_error       | real(r8)           | The combined expected error of salinity observations from all sources, |
   |                      |                    | including instrument error, model bias, and representativeness error   |
   |                      |                    | (e.g. larger or smaller grid box sizes affecting expected accuracy) in |
   |                      |                    | g/kg (psu). Values in output file are error variance, and use units of |
   |                      |                    | msu (kg/kg), so the numbers will be this value / 1000.0, squared.      |
   +----------------------+--------------------+------------------------------------------------------------------------+

| 

Modules used
------------

::

   types_mod
   time_manager_mod
   utilities_mod
   location_mod
   obs_sequence_mod
   obs_def_mod
   obs_def_ocean_mod
   obs_kind_mod

Errors and known bugs
---------------------

The code for setting observation error variances is using fixed values, and we are not certain if they are correct.
Incoming QC values larger than 0 are suspect, but it is not clear if they really signal unusable values or whether there
are some codes we should accept.

Future Plans
------------

- This converter is currently being used on WOD09 data, but the standard files generally stop with early 2009 data.
  There are subsequent additional new obs files available from the download site.

- The fractional-time field, and sometimes the day-of-month field in a small percentage of the obs have bad values. 
  The program currently discards these obs, but it may be possible to recover the original good day number and/or time of
  day. There is a subroutine at the end of the *wod_to_obs.f90* file which contains all the reject/accept/correction 
  information for the year, month, day, time fields. To accept or correct the times on more obs, edit this subroutine
  and make the necessary changes.

References
----------

.. [1] Boyer, T.P., J. I. Antonov, O. K. Baranova, C. Coleman, H. E. Garcia,
       A. Grodsky, D. R. Johnson, R. A. Locarnini, A. V. Mishonov, T.D.
       O'Brien, C.R. Paver, J.R. Reagan, D. Seidov, I. V. Smolyar, and M. M.
       Zweng, 2013: World Ocean Database 2013, NOAA Atlas NESDIS 72, S.
       Levitus, Ed., A. Mishonov, Technical Ed.; Silver Spring, MD, 209 pp., `doi:10.7289/V5NZ85MT <http://doi.org/10.7289/V5NZ85MT>`_.

.. [2] Johnson, D.R., T.P. Boyer, H.E. Garcia, R.A. Locarnini, O.K. Baranova,
       and M.M. Zweng,  2009. World Ocean Database 2009 Documentation. Edited
       by Sydney Levitus. NODC Internal Report 20, NOAA Printing Office, Silver
       Spring, MD, 175 pp., http://www.nodc.noaa.gov/OC5/WOD09/pr_wod09.html.
