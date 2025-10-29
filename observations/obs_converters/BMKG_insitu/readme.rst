================================
BMKG In-situ Floats and Drifters
================================
This utility converts Indonesian BMKG/Baron Weather in-situ ocean observations
into a DART observation sequence (``obs_seq``) file. It currently supports:

* **SVP** (Surface Velocity Program) drifters: sea surface temperature (SST), U, V
* **ARVOR-C** profiling floats: temperature/salinity profiles
* **ARVOR-I** (ice-capable) profiling floats: temperature/salinity profiles

.. contents::
   :depth: 3
   :local:

Overview
--------
Data arrive as ASCII text with a one-line header of field names followed by rows
of observations. Files may be delimited by commas (``','``) or semicolons (``';'``);
the converter auto-detects the delimiter and preserves empty cells so column
indices remain consistent between header and data.

The converter:
  * identifies the product family (SVP/ARVOR-C/ARVOR-I) using the ``process`` column,
  * parses timestamps, locations, and variables of interest,
  * converts temperatures (K to C) and pressure (dbar) to depth (m, positive downward), and 
  * writes DART observations with configurable error standard deviations.

SVP observations are written at the surface (``VERTISSURFACE``), ARVOR profiles use 
positive-down depth with ``VERTISHEIGHT`` as used by DART’s ocean conventions.

Forward operators
-----------------
* ``DRIFTER_TEMPERATURE`` uses observed SVP ``sea_temperature`` (converted to C)
* ``DRIFTER_U_CURRENT_COMPONENT`` uses observed SVP ``u`` (m/s)
* ``DRIFTER_V_CURRENT_COMPONENT`` uses observed SVP ``v`` (m/s)
* ``FLOAT_TEMPERATURE`` uses observed ARVOR ``sea_water_temperature`` (converted to C)
* ``FLOAT_SALINITY`` uses observed ARVOR ``salinity`` (PSU)

Input data expectations
-----------------------
**Common header fields**

* ``dat``  – ISO time stamp, e.g. ``YYYY-MM-DDThh:mm:ssZ``
* ``lat``  – latitude (degrees)
* ``lon``  – longitude (degrees, will be mapped to [0, 360))
* ``process`` – string label: ``SVP`` or ``DRIFTER``, ``ARVORC``, ``ARVORI``

**SVP files (drifters)**

* Variables: ``sea_temperature``, ``u`` (eastward), ``v`` (northward)
* Each row is a surface sample (depth set to 0 m; ``VERTISSURFACE``).

**ARVOR files (profiles)**

* Variables: ``sea_water_pressure`` (dbar) converted to depth, 
  ``sea_water_temperature``, ``salinity``
* Each row is a level in a vertical profile at the given time/position
  (written with ``VERTISHEIGHT``).

Time handling
^^^^^^^^^^^^^
The converter parses timestamps like ``2025-10-06T06:10:00Z`` into a DART
``time_type``. If your source provides alternative ISO forms with a trailing
``'Z'`` or timezone-less UTC, they will parse as well.

Build & run
-----------
Build like other DART converters and run with an ``input.nml`` containing
``bmkg_to_obs_nml`` (see below). One of ``file_in`` or ``file_list`` must be set.

.. code-block:: bash

   ./bmkg_to_obs

The program writes a single output sequence to ``file_out``.

Namelist
--------
This namelist is added to the rest of DART program namelists in file ``input.nml``. 
Namelists start with an ampersand '&' and terminate with a slash '/'.

.. code-block:: fortran

   &bmkg_to_obs_nml
      file_in          = '',                 ! single file (or '')
      file_list        = 'bmkg_file_list',   ! text file of paths (or '')
      file_out         = 'obs_seq.bmkg',
      avg_obs_per_file = 500000,             ! pre-allocation hint
      convert_svp      = .true.,
      convert_arvorc   = .true.,
      convert_arvori   = .true.,
      obs_error_tmp    = 0.02,               ! float temperature   (C)
      obs_error_sal    = 0.02,               ! float salinity      (PSU)
      obs_error_sst    = 0.20,               ! drifter SST         (C)
      obs_error_vel    = 0.10,               ! drifter U,V         (m/s)
      debug            = .true.
   /

.. list-table::
   :header-rows: 1
   :widths: 22 12 26 40

   * - **Namelist item**
     - **Type**
     - **Default**
     - **Description**

   * - ``file_in``
     - character(256)
     - ``''``
     - Path to a single ASCII input file. Mutually exclusive with ``file_list``.

   * - ``file_list``
     - character(256)
     - ``''``
     - Path to a text file containing one input path per line. Use for batch runs.

   * - ``file_out``
     - character(256)
     - ``'obs_seq.bmkg'``
     - Name of the output DART observation sequence file (replaced if it exists).

   * - ``avg_obs_per_file``
     - integer
     - ``500000``
     - Pre-allocation hint for obs sequence size. Set near your typical per-file count.

   * - ``convert_svp``, ``convert_arvorc``, ``convert_arvori``
     - logical
     - ``.true.``
     - Enable/disable ingestion per product type.

   * - ``obs_error_tmp``
     - real(r8)
     - ``0.02``
     - Temperature SD (C) for ARVOR floats.

   * - ``obs_error_sal``
     - real(r8)
     - ``0.02``
     - Salinity SD (PSU) for ARVOR floats.

   * - ``obs_error_sst``
     - real(r8)
     - ``0.20``
     - SST SD (C) for SVP drifters.

   * - ``obs_error_vel``
     - real(r8)
     - ``0.10``
     - U,V SD (m/s) for SVP drifters.

   * - ``debug``
     - logical
     - ``.true.``
     - Verbose logging of parsing, unit conversions, and per-file summaries.

Key routines
------------
* **detect_product**: Opens the file, detects delimiter (semicolon/comma),
  tokenizes the header, locates ``process``, and classifies the data as
  SVP / ARVOR-C / ARVOR-I.

* **split_fields / normalize_delims**: Normalizes delimiters and preserves
  empty fields so the header column order matches data rows. Delegates
  tokenization to DART’s ``get_args_from_string`` parser.

* **find_field**: Case-insensitive header lookup for required fields
  (e.g., ``dat``, ``lat``, ``lon``, ``sea_water_pressure``,
  ``sea_water_temperature``, ``salinity``, ``sea_temperature``, ``u``, ``v``).

* **parse_svp**: Reads surface drifter rows; converts SST from Kelvin to C;
  sets surface depth (0 m) and emits SST, U, V observations.

* **parse_arvor_profile**: Reads profile rows; converts pressure (dbar) to
  depth (m, +down), temperature (K to C), and emits temperature and salinity
  observations for each level.

* **add_svp / add_arvor_profile**: Creates DART observations with appropriate
  kinds (see *Forward operators*), error SDs from the namelist, and append
  them to the observation sequence in time order.

* **depth_from_pressure**: UNESCO-style approximation with latitude-dependent
  gravity; converts pressure (dbar) to depth (m, positive downward).

Observation errors
^^^^^^^^^^^^^^^^^^
By default, the converter assigns constant per-type SDs from the namelist:

* **ARVOR**: ``obs_error_tmp`` (C) for ``FLOAT_TEMPERATURE`` and
  ``obs_error_sal`` (PSU) for ``FLOAT_SALINITY``.
* **SVP**: ``obs_error_sst`` (C) for ``DRIFTER_TEMPERATURE`` and
  ``obs_error_vel`` (m/s) for ``DRIFTER_U/V_CURRENT_COMPONENT``.

These are meant as reasonable starting points. You can tune them by region,
vendor, or season by adjusting the namelist values.

Output
------
The program writes a single ``obs_seq`` file (``file_out``) containing all
observations from the provided input(s). A sample output would look like this: 

.. code-block:: text

 obs_sequence
 obs_type_definitions
            5    
           15 FLOAT_SALINITY     
           16 FLOAT_TEMPERATURE     
           17 DRIFTER_U_CURRENT_COMPONENT    
           18 DRIFTER_V_CURRENT_COMPONENT    
           20 DRIFTER_TEMPERATURE             
   num_copies:       1       num_qc:            1    
   num_obs:          358  max_num_obs:          358  
 In-situ observation                                    
 In-situ QC              
   first:          1  last:           358 
  OBS            1    
    26.1200000000000     
   0.000000000000000E+000
          269           2          -1   
 obdef
 loc3d
     0.7560975759149676        0.8255407361933162E-02     0.000000000000000     -1   
 kind
           20   
      0     155141
   4.000000000000001E-002
  OBS            2    
  -6.184206827891370E-003
   0.000000000000000E+000
            1           3          -1   
 obdef
 loc3d
     0.7560975759149676        0.8255407361933162E-02     0.000000000000000     -1   
 kind
           17   
      0     155141
   1.000000000000000E-002


Troubleshooting
^^^^^^^^^^^^^^^
* **“Unknown product; skipping.”**  
  Ensure the input has a ``process`` column with one of:
  ``SVP``, ``ARVORC``, or ``ARVORI``.
* **Header index drift or wrong columns parsed**  
  Files with empty cells are supported, but if a row is badly malformed the
  converter may skip it. Use ``debug = .true.`` to print the offending line.
* **“Unable to read time variable.”**  
  The converter expects ISO timestamps like ``YYYY-MM-DDThh:mm:ssZ``.
* **Unexpected units**  
  SVP ``u/v`` must be m/s; temperatures must be Kelvin in the input (the
  converter subtracts 273.15). Salinity is assumed in PSU. If your feeds use
  different units, convert upstream.

Performance & tips
^^^^^^^^^^^^^^^^^^
* Use ``file_list`` for many files; it improves logging and pre-allocation.
* Set ``avg_obs_per_file`` near your typical count to minimize reallocations.
* Turn off verbose prints (``debug = .false.``) for production runs.
* Longitudes are normalized to [0, 360).

Further reading
---------------
For background on the observing platforms supported by this converter:

**Surface Velocity Program (SVP) drifters**

* `Global Drifter Program (GDP): <https://gdp.ucsd.edu/ldl/svpb/>`_ *"SVP drifter overview"* -- NOAA and Scripps-maintained description of the standard SVP and SVP-B drifters used globally.  
* `Lumpkin & Pazos (2007): <https://www.aoml.noaa.gov/phod/docs/LumpkinPazos.pdf>`_ *"Measuring surface currents with Surface Velocity Program drifters"* -- Classic technical paper detailing instrument design, calibration, and data characteristics.  
* `Ocean Observers: <https://www.oceanobservers.org/Learn-about-observing-the-ocean/Drifting-buoys-DBCP>`_ *"Drifting buoys (DBCP)"* -- General introduction to the international drifter network and its coordination under WMO/IOC.

**ARVOR profiling floats**

* `Argo Program: <https://argo.ucsd.edu/how-do-floats-work/float-types/>`_ *"Float types"* -- Official overview of the various Argo float families, including ARVOR, PROVOR, and DEEP-Arvor.  
* `Le Reste et al. (2016): <https://archimer.ifremer.fr/doc/00318/42969/42473.pdf>`_ *"Deep-Arvor: A new profiling float to extend the Argo observations to 4000 m"* -- Detailed design and performance study of the deep variant of the ARVOR float.  
* `NKE Instrumentation: <https://nke-instrumentation.com/standard-profiling-floats/>`_ *"Profiling floats"* -- Manufacturer specifications for the ARVOR and PROVOR families (pressure range, sensors, telemetry).
