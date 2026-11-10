.. _cmems_ssh_to_obs:

===============================
CMEMS SSH Along-Track Converter
===============================

.. contents:: 
   :depth: 3
   :local:

Overview
--------
This converter ingests **Copernicus Marine Environment Monitoring Service (CMEMS)**
along-track sea surface height (SSH) observations and converts them into a DART
``obs_seq`` file.

The supported product is:

- *Global Ocean Along Track L3 Sea Surface Heights NRT*
  (Product ID: ``SEALEVEL_GLO_PHY_L3_NRT_008_044``)

The converter reads CSV files exported from CMEMS and extracts **filtered sea level
anomaly (SLA)** observations for assimilation into ocean models such as ROMS.

The resulting observations are assigned the DART type: ``SATELLITE_SSH``. The 
observations in the sequence file are in units of meters. 

Data Description
----------------
The converter uses the following variable: ``SLA_FILTERED``

This variable represents **low-pass filtered sea level anomaly**, including corrections for 
dynamic atmospheric correction (DAC), ocean tides, internal tides, and long-wavelength errors.
Compared to ``SLA_UNFILTERED``, this field is smoother, less noisy, and more 
representative of mesoscale variability. This makes it more appropriate for 
comparison with model sea surface height (e.g., ROMS ``zeta``).

Observation Error
-----------------
The converter assigns a uniform observation error standard deviation to all SSH
observations. The error value is specified by the user through the converter
namelist and is written to the output ``obs_seq`` file as observation error
variance.

Optional Post-Processing
^^^^^^^^^^^^^^^^^^^^^^^^
Additional filtering and depth-dependent error adjustment can be performed
after the converter is run using the utility script:

::

   models/ROMS_Rutgers/preprocess_ocean_obs.py

This script is intended for model-specific preprocessing of ocean observations
prior to assimilation. It can:

- remove observations located in shallow-water regions,
- apply depth-dependent observation error models,
- process only selected observation types,
- preserve all other observation types unchanged.

The script currently uses ROMS bathymetry and land masks from a ROMS restart
or history file.

For SSH observations, a commonly used workflow is to remove observations in
water shallower than 200 m and increase the observation error in coastal
regions where representativeness errors and unresolved processes are larger.

The depth-dependent error model is:

::

   sigma = sigma_min + (sigma_max - sigma_min) * exp(-h / h0)

where:

- ``sigma_min``: deep-ocean observation error standard deviation,
- ``sigma_max``: shallow-water observation error standard deviation,
- ``h``: local bathymetric depth (m),
- ``h0``: transition depth scale (m).

Example usage:

::

   python preprocess_ocean_obs.py obs_seq.all obs_seq.trim \
       --roms-file roms_hist.nc \
       --obs-type SATELLITE_SSH

Typical output:

::

   Reading obs_seq:  obs_seq.all
   Reading ROMS grid: roms_hist.nc

   Obs types in file: ['BOTTLE_TEMPERATURE', 'CTD_TEMPERATURE',
   'FLOAT_TEMPERATURE', 'MOORING_TEMPERATURE',
   'SATELLITE_SSH', 'XBT_TEMPERATURE']

   Input obs (SATELLITE_SSH): 3042
   Looking up bathymetric depths ...

   Summary (SATELLITE_SSH)
     Kept:    2305  (depth > 200 m)
     Removed: 737
     Depth range kept: 205.0 – 5000.0 m
     Error model: sigma = 0.04 + (0.08 - 0.04) * exp(-h / 500.0)

     Other obs types kept unchanged: 571

     Wrote: obs_seq.trim


Time Handling
-------------
Observation times are read from ISO-formatted strings, for example:

::

   2026-04-16T08:43:54.000Z

These are converted to DART ``time_type`` using the Gregorian calendar.

Usage
-----
1. Prepare a list of input CSV files and place them in a text file, for instance ``ssh_list.txt`` 

2. Configure the namelist in ``input.nml``:

::

   &cmems_ssh_to_obs_nml
      file_list        = 'ssh_file_list.txt'
      file_out         = 'obs_seq.ssh'
      avg_obs_per_file = 500000    
      obs_error_sd     = 0.04
      debug            = .true.
   /

+----------------------+-----------------------+---------------------------+--------------------------------------------------------------+
| Namelist item        | Type                  | Default Values            | Description                                                  |
+======================+=======================+===========================+==============================================================+
| ``file_list``        | character(len=256)    | ``''`` (empty string)     | Path to a text file containing a list of input CSV files,    |
|                      |                       |                           | one filename per line.                                       |
+----------------------+-----------------------+---------------------------+--------------------------------------------------------------+
| ``file_out``         | character(len=256)    | ``'obs_seq.ssh'``         | Name of the DART observation sequence output file. If the    |
|                      |                       |                           | file already exists, it will be replaced.                    |
+----------------------+-----------------------+---------------------------+--------------------------------------------------------------+
| ``avg_obs_per_file`` | integer               | ``500000``                | Estimated average number of observations per input file.     |   
|                      |                       |                           | Used to pre-allocate sequence memory efficiently.            |  
+----------------------+-----------------------+---------------------------+--------------------------------------------------------------+
| ``obs_error_sd``     | real(r8)              | ``0.04``                  | Uniform observation error standard deviation (meters)        |
|                      |                       |                           | assigned to all SSH observations.                            |
+----------------------+-----------------------+---------------------------+--------------------------------------------------------------+
| ``debug``            | logical               | ``.true.``                | If true, prints diagnostic output including sample           |
|                      |                       |                           | observations and processing information.                     |
+----------------------+-----------------------+---------------------------+--------------------------------------------------------------+

1. Run the converter: ``./cmems_ssh_to_obs``

2. Output: The resulting observation sequence file is stored in ``obs_seq.ssh``

When ``debug`` is turned on, the progress of the converter can be monitored 
which would look similar to the following:

::

  Input file: #2 ../data/obs_ssh_2.csv

   CSV Header Fields (10 columns):
     1  parameter
     2  platformId
     3  platformType
     4  time
     5  longitude
     6  latitude
     7  depth
     8  pressure
     9  value
    10  valueQc

    * obs #  1, lat: -14.2924, lon: 127.8070, SSH: 0.2940, QC: 0, date: 2026 Apr 14 19:26:36
    * obs #  2, lat: -14.2351, lon: 127.8010, SSH: 0.2820, QC: 0, date: 2026 Apr 14 19:26:37
    * obs #  3, lat: -14.1778, lon: 127.7949, SSH: 0.2690, QC: 0, date: 2026 Apr 14 19:26:38
    * obs #  4, lat: -14.1206, lon: 127.7888, SSH: 0.2560, QC: 0, date: 2026 Apr 14 19:26:39
    * obs #  5, lat: -14.0633, lon: 127.7827, SSH: 0.2420, QC: 0, date: 2026 Apr 14 19:26:40
    * obs #  6, lat: -14.0060, lon: 127.7766, SSH: 0.2280, QC: 0, date: 2026 Apr 14 19:26:41
    * obs #  7, lat: -13.9487, lon: 127.7705, SSH: 0.2130, QC: 0, date: 2026 Apr 14 19:26:42
    * obs #  8, lat: -13.8915, lon: 127.7645, SSH: 0.2000, QC: 0, date: 2026 Apr 14 19:26:43
    * obs #  9, lat: -13.8342, lon: 127.7584, SSH: 0.1910, QC: 0, date: 2026 Apr 14 19:26:43
    * obs # 10, lat: -13.7769, lon: 127.7523, SSH: 0.1850, QC: 0, date: 2026 Apr 14 19:26:44

 > Ready to write 3042 observations:
   write_obs_seq  opening formatted observation sequence file "obs_seq.ssh"
   write_obs_seq  closed observation sequence file "obs_seq.ssh"
   cmems_ssh_to_obs Finished successfully.

References
----------
- Copernicus Marine Service:
  https://marine.copernicus.eu

- Product documentation:
  https://data.marine.copernicus.eu/product/SEALEVEL_GLO_PHY_L3_NRT_008_044
