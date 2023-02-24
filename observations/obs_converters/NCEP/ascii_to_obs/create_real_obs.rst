PROGRAM create_real_obs
=======================

Overview
--------

Translating NCEP BUFR files into DART obs_seq.out files (input file to filter) is a 2 stage process. The first stage
uses NCEP software to translate the BUFR file into an "intermediate" text file. This is described in
:doc:`../prep_bufr/prep_bufr`. The second step is to translate the intermediate files into an ``obs_seq.out`` files,
which is done by ``create_real_obs``, as described in this document.

This program provides a number of options to select several observation types (radiosonde, aircraft, and satellite data,
etc.) and the DART observation variables (U, V, T, Q, Ps) which are specified in its optional namelist interface
``&ncepobs_nml`` which may be read from file ``input.nml``.

Instructions
------------

-  Go to DART/observations/obs_converters/NCEP/ascii_to_obs/work
-  Use ``quickbuild.sh`` to compile all executable programs in the directory.
-  Make appropriate changes to the ``&ncep_obs_nml`` namelist in ``input.nml``, as follows.
-  run ``create_real_obs``.

The selection of any combinations of the specific observation fields (T, Q, U/V, and surface pressure) and types
(radiosonde, aircraft reports, or satellite wind, etc.) is made in the namelist ``&ncepobs_nml``. All the available
combinations of fields X types (i.e. ADPUPA and obs_U) will be written to the obs_seq file. (You will be able to select
which of those to use during an assimilation in another namelist (``assimilate_these_obs``, in ``&obs_kind_nml``), so be
sure to include all the fields and types you might want.) You should change ``Obsbase`` to the pathname of the decoded
PREPBUFR text data files. Be sure that ``daily_file`` is set to .TRUE. to create a single 24 hour file; .FALSE. converts
input files one-for-one with output files. The default action is to tag each observation with the exact time it was
taken and is the recommended setting. However, if you want to bin the observations in time, for example to do additional
post-processing, the time on all observations in the window can be overwritten and set to the nearest synoptic time
(e.g. 0Z, 6Z, 12Z, or 18Z), by setting ``obs_time`` to false.

Generally you will want to customize the namelist for your own use. For example, here is a sample namelist:

::

   &ncepobs_nml
     year = 2007, 
     month = 3,
     day = 1,
     tot_days = 31,
     max_num = 700000,
     ObsBase = '../prep_bufr/work/temp_obs.'
     select_obs  = 1,
     ADPUPA = .true., 
     AIRCAR = .false.,  
     AIRCFT = .true., 
     SATEMP = .false., 
     SFCSHP = .false.,
     ADPSFC = .false.,  
     SATWND = .true., 
     obs_U  = .true., 
     obs_V  = .true.,
     obs_T  = .true.,
     obs_PS = .false.,
     obs_QV = .false.,
     daily_file = .true.
     obs_time = .true.,
   /

   &obs_sequence_nml
     write_binary_obs_sequence = .false.  
   /

This will produce daily observation sequence files for the period of March 2007, which have the selected observation
types and fields; T, U, and V from radiosondes (ADPUPA) and aircraft (AIRCFT). No surface pressure or specific humidity
would appear in the obs_seq files, nor observations from ACARS, satellites, and surface stations. The output files look
like "obs_seq200703dd", with dd = 1,...,31.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &ncepobs_nml
      year       = 2003,
      month      = 1,
      day        = 1,
      tot_days   = 31,
      max_num    = 800000,
      select_obs = 0,
      ObsBase    = 'temp_obs.',
      ADPUPA     = .false., 
      AIRCAR     = .false., 
      AIRCFT     = .false., 
      SATEMP     = .false., 
      SFCSHP     = .false., 
      ADPSFC     = .false., 
      SATWND     = .false.,
      obs_U      = .false., 
      obs_V      = .false., 
      obs_T      = .false.,
      obs_PS     = .false.,
      obs_QV     = .false.,
      daily_file = .true.,
      obs_time   = .true.,
      lon1       =   0.0,
      lon2       = 360.0,
      lat1       = -90.0,
      lat2       =  90.0  
   /

| 

.. container::

   +------------------------------+-----------------------+---------------------------------------+
   | Item                         | Type                  | Description                           |
   +==============================+=======================+=======================================+
   | year, month, day             | integer               | Beginning year, month, day of the     |
   |                              |                       | observation period.                   |
   +------------------------------+-----------------------+---------------------------------------+
   | tot_days                     | integer               | Total days in the observation period. |
   |                              |                       | The converter cannot cross month      |
   |                              |                       | boundaries.                           |
   +------------------------------+-----------------------+---------------------------------------+
   | max_num                      | integer               | Maximum observation number for the    |
   |                              |                       | current one day files.                |
   +------------------------------+-----------------------+---------------------------------------+
   | select_obs                   | integer               | Controls whether to select a subset   |
   |                              |                       | of observations from the NCEP BUFR    |
   |                              |                       | decoded daily ascii files.            |
   |                              |                       |                                       |
   |                              |                       | -  0 = All observations are selected. |
   |                              |                       | -  1 = Select observations using the  |
   |                              |                       |    logical parameters below.          |
   +------------------------------+-----------------------+---------------------------------------+
   | daily_file                   | logical               | Controls timespan of observations in  |
   |                              |                       | each obs_seq file:                    |
   |                              |                       |                                       |
   |                              |                       | -  true = 24 hour spans (3:01Z to     |
   |                              |                       |    3:00Z of the next day). Filenames  |
   |                              |                       |    have the form obs_seqYYYYMMDD.     |
   |                              |                       | -  false = 6 hour spans (3:01Z to     |
   |                              |                       |    9:00Z, 9:01Z to 15:00Z, 15:01Z to  |
   |                              |                       |    21:00Z, and 21:01Z to 3:00Z of the |
   |                              |                       |    next day. Filenames have the form  |
   |                              |                       |    obs_seqYYYYMMDDHH, where HH is 06, |
   |                              |                       |    12, 18, and 24.                    |
   +------------------------------+-----------------------+---------------------------------------+
   | ObsBase                      | character(len=129)    | Path that contains the decoded NCEP   |
   |                              |                       | BUFR daily observation files. To work |
   |                              |                       | with the example scripts this should  |
   |                              |                       | be 'temp_obs.', or if it includes a   |
   |                              |                       | pathname then it should end with a    |
   |                              |                       | '/temp_obs.'                          |
   +------------------------------+-----------------------+---------------------------------------+
   | include_specific_humidity,   | logical               | Controls which moisture observations  |
   | include_relative_humidity,   |                       | are created. The default is to create |
   | include_dewpoint             |                       | only specific humidity obs, but any,  |
   |                              |                       | all, or none can be requested. Set to |
   |                              |                       | .TRUE. to output that obs type,       |
   |                              |                       | .FALSE. skips it.                     |
   +------------------------------+-----------------------+---------------------------------------+
   | ADPUPA                       | logical               | Select the NCEP type ADPUPA           |
   |                              |                       | observations which includes land and  |
   |                              |                       | ship launched radiosondes and pibals  |
   |                              |                       | as well as a few profile dropsonde.   |
   |                              |                       | This involves, at 00Z and 12Z, about  |
   |                              |                       | 650 - 1000 stations, and at 06Z and   |
   |                              |                       | 18Z (which are mostly pibals), about  |
   |                              |                       | 150 - 400 stations.                   |
   +------------------------------+-----------------------+---------------------------------------+
   | AIRCFT                       | logical               | Select the NCEP type AIRCFT           |
   |                              |                       | observations, which includes          |
   |                              |                       | commercial, some military and         |
   |                              |                       | reconnaissance reports. They are      |
   |                              |                       | flight level reports.                 |
   +------------------------------+-----------------------+---------------------------------------+
   | AIRCAR                       | logical               | Select the NCEP type AIRCAR           |
   |                              |                       | observations, which includes data     |
   |                              |                       | from aircraft takeoff and landings.   |
   |                              |                       | Sometimes referred to as ACARS obs.   |
   +------------------------------+-----------------------+---------------------------------------+
   | SATEMP                       | logical               | Select the NCEP type SATEMP           |
   |                              |                       | observations, which includes NESDIS   |
   |                              |                       | ATOVS virtual temperature soundings.  |
   +------------------------------+-----------------------+---------------------------------------+
   | SFCSHP                       | logical               | Select the NCEP type SFCSHP           |
   |                              |                       | observations, which includes surface  |
   |                              |                       | marine (ship, buoy, c-man) reports.   |
   +------------------------------+-----------------------+---------------------------------------+
   | ADPSFC                       | logical               | Select the NCEP type ADPSFC           |
   |                              |                       | observations, which includes surface  |
   |                              |                       | land synoptic station reports.        |
   +------------------------------+-----------------------+---------------------------------------+
   | SATWND                       | logical               | Select the NCEP type SATWND           |
   |                              |                       | observations, which includes winds    |
   |                              |                       | derived from satellite cloud drift    |
   |                              |                       | analysis.                             |
   +------------------------------+-----------------------+---------------------------------------+
   | obs_U                        | logical               | Select u-component of wind            |
   |                              |                       | observations.                         |
   +------------------------------+-----------------------+---------------------------------------+
   | obs_V                        | logical               | Select v-component of wind            |
   |                              |                       | observations.                         |
   +------------------------------+-----------------------+---------------------------------------+
   | obs_T                        | logical               | Select temperature observations.      |
   +------------------------------+-----------------------+---------------------------------------+
   | obs_PS                       | logical               | Select surface pressure observations. |
   +------------------------------+-----------------------+---------------------------------------+
   | obs_QV                       | logical               | Select specific humidity              |
   |                              |                       | observations.                         |
   +------------------------------+-----------------------+---------------------------------------+
   | lon1                         | real                  | Western longitude bound of            |
   |                              |                       | observations to keep.                 |
   +------------------------------+-----------------------+---------------------------------------+
   | lon2                         | real                  | Eastern longitude bound of            |
   |                              |                       | observations to keep. Can be less     |
   |                              |                       | than lon1 if region crosses prime     |
   |                              |                       | meridian.                             |
   +------------------------------+-----------------------+---------------------------------------+
   | lat1                         | real                  | Lower latitude bound of observations  |
   |                              |                       | to keep.                              |
   +------------------------------+-----------------------+---------------------------------------+
   | lat2                         | real                  | upper latitude bound of observations  |
   |                              |                       | to keep.                              |
   +------------------------------+-----------------------+---------------------------------------+
   | obs_time                     | logical               | If .true. use the full time in the    |
   |                              |                       | input data. To force all observation  |
   |                              |                       | times in the output to the synoptic   |
   |                              |                       | time (e.g. 0Z, 6Z, 12Z, or 18Z) set   |
   |                              |                       | this to .false. (not recommended).    |
   +------------------------------+-----------------------+---------------------------------------+

| 

Modules used
------------

::

   types_mod
   utilities_mod
   obs_utilities_mod
   obs_sequence_mod
   obs_kind_mod
   obs_def_mod
   assim_model_mod
   model_mod
   cov_cutoff_mod
   location_mod
   random_seq_mod
   time_manager_mod
   null_mpi_utilities_mod
   real_obs_mod

Files
-----

-  temp_obs.yyyymmdd; (input) NCEP BUFR (decoded/intermediate) observation file(s) Each one has 00Z of the next day on
   it.
-  input.nml; the namelist file used by create_real_obs.
-  obs_seqYYYYMMDD[HH]; (output) the obs_seq files used by DART.

References
----------

-  DART/observations/obs_converters/NCEP/prep_bufr/docs (NCEP text files describing the BUFR files)
