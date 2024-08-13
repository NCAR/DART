GTSPP Observations
==================

Overview
--------

GTSPP (Global Temperature-Salinity Profile Program) data measures vertical profiles of ocean temperature and salinity.
The `GTPSS home page <https://www.ncei.noaa.gov/products/global-temperature-and-salinity-profile-programme>`__ 
has detailed information about the repository,
observations, and datasets. The programs in this directory convert from the netcdf files found in the repository into
DART observation sequence (obs_seq) file format.

Data sources
------------

Data from the GTSPP can be downloaded interactively from
`the GTSPP data server <http://www.nodc.noaa.gov/cgi-bin/gtspp/gtsppform01.cgi>`__. It is delivered in
`netCDF <http://www.unidata.ucar.edu/software/netcdf>`__ file format, one vertical profile per netCDF file.

Currently each vertical profile is stored in a separate file, so converting a months's worth of observations involves
downloading many individual files. The converter program can take a list of input files, so it is easy to collect a
month of observations together into a single output file with one execution of the converter program.

The units in the source file are degrees C for temperature, g/kg for salinity, and so far we have not found any error
information (not quality control, but observation instrument error values). There is probably instrument source
information encoded in these files, but so far we don't have the key. The quality control values are read and only those
with a QC of 1 are retained.

.. Note::  

   The GTSPP data is in PSU, the dart observation converter ``gtspp_to_obs`` converts the observations to MSU.
 
   | PSU = g/kg
   | MSU = PSU/1000 = kg/kg

Programs
--------

The data is distributed in `netCDF <http://www.unidata.ucar.edu/software/netcdf>`__ file format. DART requires all
observations to be in a proprietary format often called DART "obs_seq" format. The files in this directory, a
combination of C shell scripts and a Fortran source executable, do this data conversion.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &gtspp_to_obs_nml
      gtspp_netcdf_file     = '1234567.nc'
      gtspp_netcdf_filelist = 'gtspp_to_obs_filelist'
      gtspp_out_file        = 'obs_seq.gtspp'
      avg_obs_per_file      = 500
      debug                 = .false.
   /

| 

.. container::

   +-----------------------+--------------------+-----------------------------------------------------------------------+
   | Item                  | Type               | Description                                                           |
   +=======================+====================+=======================================================================+
   | gtspp_netcdf_file     | character(len=128) | The input filename when converting a single profile. Only one of the  |
   |                       |                    | two file or filelist items can have a valid value, so to use the      |
   |                       |                    | single filename set the list name 'gtspp_netcdf_filelist' to the      |
   |                       |                    | empty string (' ').                                                   |
   +-----------------------+--------------------+-----------------------------------------------------------------------+
   | gtspp_netcdf_filelist | character(len=128) | To convert a series of profiles in a single execution create a text   |
   |                       |                    | file which contains each input file, in ascii, one filename per line. |
   |                       |                    | Set this item to the name of that file, and set 'gtspp_netcdf_file'   |
   |                       |                    | to the empty string (' ').                                            |
   +-----------------------+--------------------+-----------------------------------------------------------------------+
   | gtspp_out_file        | character(len=128) | The output file to be created. To be compatible with earlier versions |
   |                       |                    | of this program, if this file already exists it will be read in and   |
   |                       |                    | the new data will be inserted into that file.                         |
   +-----------------------+--------------------+-----------------------------------------------------------------------+
   | avg_obs_per_file      | integer            | The code needs an upper limit on the number of observations generated |
   |                       |                    | by this program. It can be larger than the actual number of           |
   |                       |                    | observations converted. The total number of obs is computed by        |
   |                       |                    | multiplying this number by the number of input files. If you get an   |
   |                       |                    | error because there is no more room to add observations to the output |
   |                       |                    | file, increase this number.                                           |
   +-----------------------+--------------------+-----------------------------------------------------------------------+
   | debug                 | logical            | If true, output more debugging messages.                              |
   +-----------------------+--------------------+-----------------------------------------------------------------------+

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
   netcdf

Known Bugs
----------

Does not have correct code for setting observation error variance yet. Also, not sure if the incoming data qc is strict enough.

