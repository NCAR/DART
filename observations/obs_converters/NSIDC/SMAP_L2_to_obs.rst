PROGRAM ``SMAP_L2_to_obs``
==========================

Overview
--------

Soil Moisture Active Passive (SMAP) passive microwave radiometer 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
to DART Observation Sequence Converter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This routine is designed to convert the SMAP L2 Radiometer Half-Orbit
36 km EASE-Grid Soil Moisture, `Soil moisture product (Versions 4-8)
<https://nsidc.org/data/SPL2SMP>`__ 
to a DART observation sequence file.  

Quoting the `NSIDC: <https://nsidc.org>`__

`This Level-2 (L2) soil moisture product provides estimates of global land 
surface conditions retrieved by the Soil Moisture Active Passive (SMAP) 
passive microwave radiometer during 6:00 a.m. descending and 6:00 p.m. 
ascending half-orbit passes. SMAP L-band brightness temperatures are used  to
derive soil moisture data, which are then resampled to an Earth-fixed, global,
cylindrical 36 km Equal-Area Scalable Earth Grid, Version 2.0 (EASE-Grid 2.0).`

`Data Set ID: SPL2SMP
SMAP L2 Radiometer Half-Orbit 36 km EASE-Grid Soil Moisture`

`Surface soil moisture (0-5 cm) in m3/m3 derived from brightness temperatures
(TBs) is output on a fixed global 36 km EASE-Grid 2.0. Also included are 
brightness temperatures in kelvin representing the weighted average of 
Level-1B brightness temperatures whose boresights fall within a 36 km 
EASE-Grid 2.0 cell.`

.. Important::

  ``SMAP_L2_to_obs`` uses an
  observation error standard deviation of 0.01 m3/m3 or 20% of the soil moisture 
  value, whatever is higher. **These numbers have no scientific basis and 
  should be thoroughly explored.**  The data files I have explored have 
  a variable soil_moisture_error but these appear to be 
  empty - i.e. full of _FillValue values. A better way
  of specifying the observation error standard deviation is needed.

.. Important::

  SMAP_L2_to_obs has only
  been thoroughly tested with the half-orbit files - these files have names like
  'SMAP_L2_SM_P_02526_D_20150723T070211_R12170_001.h5'. 

.. Important::

  SMAP_L2_to_obs is not compatible with SMAP L3 files as they are formatted 
  differently



The workflow is usually: 


   1. Download the data for the period in question 
      (see DATA SOURCES below)
   2. Build the DART executables with support for the soil moisture observations.
      This is done by running preprocess with 
      ``obs_def_land_mod.f90`` in the list of input_files
      for ``preprocess_nml``.
   3. Convert each data file individually using ``SMAP_L2_to_obs``
   4. Combine or subset all output files for the region and timeframe of interest into one file 
      using :doc:`../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool`


For some models (CLM, for example), it is required to reorganize the observation sequence
files into a series of files that contains ONLY the observations for each assimilation.
This can be achieved with the `~/models/clm/shell_scripts/makedaily.sh`
script. Since there are subtleties for each model, 
makedaily.sh is generally found in the
shell_scripts directory of the model.

NAMELIST
--------

This namelist is read from the file input.nml.
Namelists start with an ampersand
'&amp;' and terminate with a slash '/'.
Character strings that contain a '/' must be
enclosed in quotes to prevent them from
prematurely terminating the namelist.

::

  &SMAP_L2_to_obs_nml
     input_file_list     = 'file_list.txt'
     obs_out_file        = 'obs_seq.out'
     max_num_input_files = 10
     verbose             = .false.
     /

Description of namelist variables:

+--------------------+--------------------+---------------------------------------------------------------------------+
| Contents           | Type               | Description                                                               |
+====================+====================+===========================================================================+
| input_file_list    | character(len=256) | Name of the file containing the list of input data files.                 |
|                    |                    | Each input data file name must be on a separate line. No blank lines      |
|                    |                    | are allowed. This may be a relative or absolute filename.                 |
+--------------------+--------------------+---------------------------------------------------------------------------+
| obs_out_file       | character(len=256) | Name of the output observation sequence file.                             |
+--------------------+--------------------+---------------------------------------------------------------------------+
| max_num_input_files| integer            | The maximum number of filenames in 'input_file_list' to convert.          |
|                    |                    | This should be a small number - converting all frequencies, all           |
|                    |                    | polarizations, both passes into one file is not recommended               |
+--------------------+--------------------+---------------------------------------------------------------------------+
| verbose            | logical            | Print extra information during the `SMAP_L2_to_obs` execution.            |
+--------------------+--------------------+---------------------------------------------------------------------------+



Data Sources
------------

The SMAP L2 Radiometer Half-Orbit 36 km EASE-Grid Soil Moisture, Version 4 
`data <https://nsidc.org/data/SPL2SMP/versions/4>`__

Data Citation
-------------

The following example shows how to cite the use of this data set in a publication.
For more information, see the `Use and Copyright web page <http://nsidc.org/about/use_copyright.html>`__

  O'Neill, P. E., S. Chan, E. G. Njoku, T. Jackson, and R. Bindlish. 2016. 
  SMAP L2 Radiometer Half-Orbit 36 km EASE-Grid Soil Moisture, Version 4. 
  [Indicate subset used]. 
  Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center.
  doi: http://dx.doi.org/10.5067/XPJTJT812XFY. [Date Accessed].


PROGRAMS
--------

The ``SMAP_L2_to_obs.f90`` file is the source
for the converter program.
To compile and test,
go into the work subdirectory and run the ``quickbuild.sh``
script to build the converter and a couple of general purpose utilities.
:doc:`../../../assimilation_code/programs/advance_time/advance_time`
helps with calendar and time computations, and the
:doc:`../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool`
manipulates DART observation files once they have been created.



DECISIONS YOU MIGHT NEED TO MAKE
--------------------------------

See the discussion in the
:doc:`../../../guide/creating-obs-seq-real/`
introduction page about what options are available for the things you need to
specify.  These include setting a time, specifying an expected error,
setting a location, and an observation type.




Terms of Use
------------

DART software - Copyright UCAR. This open source software is provided
by UCAR, "as is", without charge, subject to all terms of use at
`http://www.image.ucar.edu/DAReS/DART/DART_download
<http://www.image.ucar.edu/DAReS/DART/DART_download>`__
