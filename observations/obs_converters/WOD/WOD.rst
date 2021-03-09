WOD Observations
================

Overview
--------

| The WOD (World Ocean Database) data is a collection of data from various sources, combined into a single format with
  uniform treatment. The `WOD 2009 page <http://www.nodc.noaa.gov/OC5/WOD09/pr_wod09.html>`__ has detailed information
  about the repository, observations, and datasets. The programs in this directory convert from the packed ASCII files
  found in the repository into DART observation sequence (obs_seq) file format.
| There are 2 sets of available files - the raw observations, and the observations binned onto standard levels. The
  recommended datasets are the ones on standard levels. The raw data can be very dense in the vertical and are not truly
  independent observations. This leads to too much certainty in the updated values during the assimilation.

Data sources
------------

Data from the WOD09 can be downloaded interactively from links on `this
page <http://www.nodc.noaa.gov/OC5/WOD09/pr_wod09.html>`__. One suggestion is to pick the 'sorted by year link' and
download the files (you can select multiple years at once) for each data type for the years of interest. Make sure to
select the standard level versions of each dataset.

UCAR/NCAR users with access to the DSS data repository can download WOD09 files from
`here <http://dss.ucar.edu/datazone/dsszone/ds285.0/#WOD09>`__. A UCAR DSS userid is required to access this page. The
files to use are named "yearly_*_STD.tar".

Requested citation if you use this data:

::

   Johnson, D.R., T.P. Boyer, H.E. Garcia, R.A. Locarnini, O.K. Baranova, and M.M. Zweng, 
   2009. World Ocean Database 2009 Documentation. Edited by Sydney Levitus. NODC 
   Internal Report 20, NOAA Printing Office, Silver Spring, MD, 175 pp.  
   Available at http://www.nodc.noaa.gov/OC5/WOD09/pr_wod09.html. 

Programs
--------

The data is distributed in a specialized packed ASCII format. In this directory is a program called ``wodFOR.f`` which
is an example reader program to print out data values from the files. The program ``wod_to_obs`` converts these packed
ASCII files into DART obs_sequence files.

As with most other DART directories, the ``work`` directory contains a ``quickbuild.csh`` script to build all necessary
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
