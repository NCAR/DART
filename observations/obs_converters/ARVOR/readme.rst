ARVOR (Profiling Floats)
========================

This utility converts in-situ profiling float (T/S profile) ASCII files 
into a DART observation sequence (`obs_seq`) file.

.. contents::
   :depth: 2
   :local:

Overview
--------
Profiling floats such as ARVOR-C and ARVOR-I (ice-capable) record seawater temperature and 
salinity at multiple depths during dives and ascents.  
This converter ingests ASCII (.csv format) files of these profiles and 
writes them into a DART obs_seq format for assimilation with DART. The 
csv files are expected to have on the first line a header with 
different fields names followed by data entries. 

The converter reads the following from the raw data: 
 * ``dat`` -- timestamp
 * ``lat`` -- latitude in degrees
 * ``lon`` -- longitude in degrees
 * ``pres`` -- pressure in dbar converted to depth in meters
 * ``sea_water_temperature`` -- temperature in Kelvin, converted to C.
 * ``salinity`` -- salinity in psu 

.. note::

   **Time handling:** The converter reads the ``dat`` column, parses the timestamp, and
   constructs a DART ``time_type`` using the Gregorian calendar. The 
   expected time format is ``YYYY-MM-DDThh:mm:ssZ``

Build & Run
-----------
Build as usual with DART converters. Then run, ``./arvor_to_hf``, with an 
`input.nml` that includes a namelist section that looks like:

.. code-block:: fortran

   &arvor_to_obs_nml
      file_in           = 'ARVOR_20251006_300534063313500.csv',
      file_list         = '',               ! path to list of files (or '')
      file_out          = 'obs_seq.arvor',
      obs_error_temp    = 0.02,             ! temperature error standard deviation (C)
      obs_error_sal     = 0.02,             ! salinity error standard deviation (PSU)
      avg_obs_per_file  = 500000,           ! pre-allocation limit
      debug             = .true.
      /

* Exactly one of `file_in` or `file_list` must be non-empty.
* Temperature values in the raw file are in Kelvin and are converted to degrees Celsius.
* Pressure values (dbar) are converted to depth in meters (+ down).
* Any missing values or bad rows (`lat`, `lon`, `pres`, `temp`, or `saln` = missing) are skipped.
* Observation kinds: `FLOAT_TEMPERATURE` and `FLOAT_SALINITY`. 

Namelist Summary
****************
.. list-table::
   :header-rows: 1
   :widths: 20 10 25 45

   * - **Namelist item**
     - **Type**
     - **Default / Allowed values**
     - **Description**

   * - ``file_in``
     - character(len=256)
     - ``''`` (empty)
     - Path to single ASCII profile file.

   * - ``file_list``
     - character(len=256)
     - ``''`` (empty)
     - Text file listing profile files.

   * - ``file_out``
     - character(len=256)
     - ``'obs_seq.arvor'``
     - Output obs_seq file name.

   * - ``obs_error_temp``
     - real(r8)
     - ``0.02``
     - Standard deviation error for temperature.

   * - ``obs_error_sal``
     - real(r8)
     - ``0.02``
     - Standard deviation error for salinity.

   * - ``avg_obs_per_file``
     - integer
     - ``500000``
     - Estimate of valid obs per file. Used for pre-allocation. Number of files times this number must be larger than the total number of output observations.

   * - ``debug``
     - logical
     - ``.true.``
     - If true, prints detailed conversion log.

Output
------
A new observation sequence file named according to `file_out`.  
Each observation record includes location (lat, lon, depth), time, and 
one of the two kinds (temperature or salinity), with error variance 
set to the square of `obs_error_temp` or `obs_error_sal`.

If ``debug = .true.``, a summary will be printed (similar to the one below):

.. code-block:: text

   Input file: #2 T_HCSV00_C_BMKG_20251006061000_ARVORC_300534063012970.csv
     >> Total number of obs: 25
     >> Valid observations count: 24
        * obs #  1, lat: -5.7368, lon: 113.8890, dep: 65.0270, T: 28.6840, S: 33.9690, date: 2025 Oct 06 06:10:00
        * obs #  2, lat: -5.7368, lon: 113.8890, dep: 62.2036, T: 28.6840, S: 33.9680, date: 2025 Oct 06 06:11:05
        * obs #  3, lat: -5.7368, lon: 113.8890, dep: 59.3703, T: 28.6840, S: 33.9680, date: 2025 Oct 06 06:12:10
        * obs #  4, lat: -5.7368, lon: 113.8890, dep: 56.5468, T: 28.6840, S: 33.9680, date: 2025 Oct 06 06:13:15
        * obs #  5, lat: -5.7368, lon: 113.8890, dep: 53.7233, T: 28.6840, S: 33.9680, date: 2025 Oct 06 06:14:20
        * obs #  6, lat: -5.7368, lon: 113.8890, dep: 50.8898, T: 28.6840, S: 33.9680, date: 2025 Oct 06 06:15:25
        * obs #  7, lat: -5.7368, lon: 113.8890, dep: 48.0663, T: 28.6840, S: 33.9680, date: 2025 Oct 06 06:16:30
        ... 

   > Ready to write 286 observations:
     write_obs_seq  opening formatted observation sequence file "obs_seq.arvor"
     write_obs_seq  closed observation sequence file "obs_seq.arvor"
     arvor_to_obs Finished successfully.

Further reading
---------------
For background on the observing platform supported by this converter:

* `Argo Program: <https://argo.ucsd.edu/how-do-floats-work/float-types/>`_ *"Float types"* -- Official overview of the various Argo float families, including ARVOR, PROVOR, and DEEP-Arvor.  
* `Le Reste et al. (2016): <https://archimer.ifremer.fr/doc/00318/42969/42473.pdf>`_ *"Deep-Arvor: A new profiling float to extend the Argo observations to 4000 m"* -- Detailed design and performance study of the deep variant of the ARVOR float.  
* `NKE Instrumentation: <https://nke-instrumentation.com/standard-profiling-floats/>`_ *"Profiling floats"* -- Manufacturer specifications for the ARVOR and PROVOR families (pressure range, sensors, telemetry).
