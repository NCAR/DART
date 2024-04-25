PROGRAM ``level4_to_obs``
=========================

Overview
--------

AmeriFlux level 4 data to DART observation sequence converter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| This routine was designed to convert the flux tower Level 4 data from the `AmeriFlux <http://ameriflux.lbl.gov>`__
  network of observations from micrometeorological tower sites. The download link and flux data format for this code
  *has been discontinued* (e.g. ``USBar2004_L4_h.txt``).  Thus if you are using flux obs converters for the first time
  *PLEASE USE* the updated ``Fluxnetfull_to_obs.f90`` converter and follow the documentation there :doc:`./fluxnetfull_to_obs`
  We have kept this code available if you still use the older Ameriflux data format.  Also this code uses a general approach
  to calculating sensible, latent and net ecosystem exchange uncertainty, that may be helpful. 
| The AmeriFlux Level 4 data products are provided in the local time of the flux tower location. DART observation sequence
  files are provided in UTC, thus this routine includes a time conversion. 

.. warning::

   There was a pretty severe bug in the converter that swapped latent heat flux and sensible heat flux. The bug was
   present through revision 7200. It was corrected on 30 Dec 2016. 

The workflow is usually:

#. download the Level 4 data for the towers and years in question (see DATA SOURCES below)
#. record the TIME ZONE, latitude, longitude, and elevation for each tower
#. build the DART executables with support for the tower observations. This is done by running ``preprocess`` with
   ``obs_def_tower_mod.f90`` in the list of ``input_files`` for ``preprocess_nml``.
#. provide basic tower information via the ``level4_to_obs_nml`` namelist since this information is not contained in the
   Level 4 data file
#. convert each Level 4 data file individually using ``level4_to_obs``
#. combine all output files for the region and timeframe of interest into one file using
   :doc:`../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool`

For some models (CLM, for example), it is required to reorganize the observation sequence files into a series of files
that contains ONLY the observations for each assimilation. This can be achieved with the `makedaily.sh <makedaily.sh>`__
script.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &level4_to_obs_nml
      text_input_file = 'textdata.input',
      obs_out_file    = 'obs_seq.out',
      year            = -1,
      timezoneoffset  = -1,
      latitude        = -1.0,
      longitude       = -1.0,
      elevation       = -1.0,
      flux_height     = -1.0,
      maxgoodqc       = 3,
      verbose         = .false.
      /

.. container::

   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | Contents        | Type               | Description                                                                 |
   +=================+====================+=============================================================================+
   | text_input_file | character(len=128) | Name of the Level 4 ASCII file of comma-separated values. This may be a     |
   |                 |                    | relative or absolute filename.                                              |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | obs_out_file    | character(len=128) | Name of the output observation sequence file.                               |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | year            | integer            | The year of the observations in the Level 4 text file.                      |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | timezoneoffset  | real               | the time zone offset (in hours) of the station. The tower observation times |
   |                 |                    | are local time, we need to convert them to GMT.                             |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | latitude        | real               | Latitude (in degrees N) of the tower.                                       |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | longitude       | real               | Longitude (in degrees E) of the tower. For internal consistency, DART uses  |
   |                 |                    | longitudes in the range [0,360]. An input value of -90 will be converted to |
   |                 |                    | 270, for example.                                                           |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | elevation       | real               | surface elevation (in meters) of the tower.                                 |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | flux_height     | real               | height (in meters) of the flux instrument on the tower.                     |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | maxgoodqc       | real               | maximum value of any observation quality control flag to pass through to    |
   |                 |                    | the output observation sequence. Keep in mind that ``filter`` has the       |
   |                 |                    | ability to discriminate on the value, so there is really little to be       |
   |                 |                    | gained by rejecting them during the conversion.                             |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | verbose         | logical            | Print extra information during the ``level4_to_obs`` execution.             |
   +-----------------+--------------------+-----------------------------------------------------------------------------+

Data sources
------------

| The data was acquired from http://cdiac.ornl.gov/ftp/ameriflux/data/Level4/Sites_ByName
| and have names like
|  ``USBar2004_L4_h.txt, USHa12004_L4_h.txt, USNR12004_L4_h.txt``, 
|  ``USSP32004_L4_h.txt, USSRM2004_L4_h.txt, USWCr2004_L4_h.txt, USWrc2004_L4_h.txt, ...``
| The Level 4 products in question are ASCII files of comma-separated values taken every 30 minutes for an entire year.
  The first line is a comma-separated list of column descriptors, all subsequent lines are comma-separated numerical
  values. The converter presently searches for the columns pertaining to *NEE_or_fMDS*, *H_f*, *LE_f*, their
  corresponding quality control fields, and those columns pertaining to the time of the observation. These values are
  mapped as follows:

+---------------+------------------+--------------------------+--------------------------+---------------------------+------------+
| Level 4 units | Level 4 variable | description              | DART type                | DART kind                 | DART units |
+===============+==================+==========================+==========================+===========================+============+
| W/m^2         | LE_f             | Latent Heat Flux         | TOWER_LATENT_HEAT_FLUX   | QTY_LATENT_HEAT_FLUX      | W/m^2      |
+---------------+------------------+--------------------------+--------------------------+---------------------------+------------+
| [0-3]         | LE_fqc           | QC for LE_f              | N/A                      | N/A                       | same       |
+---------------+------------------+--------------------------+--------------------------+---------------------------+------------+
| W/m^2         | H_f              | Sensible Heat Flux       | TOWER_SENSIBLE_HEAT_FLUX | QTY_SENSIBLE_HEAT_FLUX    | W/m^2      |
+---------------+------------------+--------------------------+--------------------------+---------------------------+------------+
| [0-3]         | H_fqc            | QC for H_f               | N/A                      | N/A                       | same       |
+---------------+------------------+--------------------------+--------------------------+---------------------------+------------+
| umolCO2/m^2/s | NEE_or_fMDS      | Net Ecosystem Production | TOWER_NETC_ECO_EXCHANGE  | QTY_NET_CARBON_PRODUCTION | gC/m^2/s   |
+---------------+------------------+--------------------------+--------------------------+---------------------------+------------+
| [0-3]         | NEE_or_fMDSqc    | QC for NEE_or_fMDS       | N/A                      | N/A                       | same       |
+---------------+------------------+--------------------------+--------------------------+---------------------------+------------+




The ``LE_fqc``, ``H_fqc``, and ``NEE_or_fMDSqc`` variables use the following convention:

   0 = original, 1 = category A (most reliable), 2 = category B (medium), 3 = category C (least reliable). (Refer to
   Reichstein et al. 2005 Global Change Biology for more information)


I am repeating the AmeriFlux `Data Fair-Use Policy <http://ameriflux.lbl.gov/Data/Pages/DataUsagePolicy.aspx>`__ because
I believe it is important to be a good scientific citizen:

   "The AmeriFlux data provided on this site are freely available and were furnished by individual AmeriFlux scientists
   who encourage their use.
   Please kindly inform in writing (or e-mail) the appropriate AmeriFlux scientist(s) of how you intend to use the data
   and of any publication plans. It is also important to contact the AmeriFlux investigator to assure you are
   downloading the latest revision of the data and to prevent potential misuse or misinterpretation of the data.
   Please acknowledge the data source as a citation or in the acknowledgments if no citation is available. If the
   AmeriFlux Principal Investigators (PIs) feel that they should be acknowledged or offered participation as authors,
   they will let you know and we assume that an agreement on such matters will be reached before publishing and/or use
   of the data for publication.
   If your work directly competes with the PI's analysis they may ask that they have the opportunity to submit a
   manuscript before you submit one that uses unpublished data. In addition, when publishing please acknowledge the
   agency that supported the research.
   Lastly, we kindly request that those publishing papers using AmeriFlux data provide reprints to the PIs providing the
   data and to the AmeriFlux archive via ameriflux.lbl.gov."

Programs
--------

The ``level4_to_obs.f90`` file is the source for the main converter program. Look at the source code where it reads the
example data file. You will almost certainly need to change the "read" statement to match your data format. The example
code reads each text line into a character buffer and then reads from that buffer to parse up the data items.

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

See the discussion in the :doc:`../../../guide/creating-obs-seq-real` page about what options are available
for the things you need to specify. These include setting a time, specifying an expected error, setting a location, and
an observation type.
