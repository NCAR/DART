===
SVP
===
This utility converts **Surface Velocity Program (SVP)** drifter data
(SST + surface currents) into a DART observation sequence (``obs_seq``) file.

.. contents::
   :depth: 2
   :local:

Overview
--------
SVP drifters measure near–surface ocean temperature and currents. For the
InaCAWO project, BMKG/Baron provide ASCII (CSV-like) files with:

- Header line of field names  
- One row per surface observation time

The converter reads:
 - ``dat``  – timestamp  
 - ``lat``  – latitude (degrees)  
 - ``lon``  – longitude (degrees)  
 - ``sea_temperature`` – SST (Kelvin; converted to °C)  
 - ``u`` / ``v`` – eastward / northward surface velocity (m/s)

Each valid row becomes up to three DART observations at the surface
(``VERTISSURFACE``): ``DRIFTER_TEMPERATURE``, ``DRIFTER_U_CURRENT_COMPONENT``,
``DRIFTER_V_CURRENT_COMPONENT``

.. note::

   **Time handling:** The converter reads the ``dat`` column, parses the timestamp, and
   constructs a DART ``time_type`` in the Gregorian calendar.
   Expected format: ``YYYY-MM-DDThh:mm:ssZ``

Build & run
-----------
Build like other DART converters, then run with an ``input.nml`` that includes
``svp_to_obs_nml``:

.. code-block:: bash

   ./svp_to_obs

The program writes a single output sequence defined by ``file_out``.

Namelist
********
This namelist is added to ``input.nml``:

.. code-block:: fortran

   &svp_to_obs_nml
      file_in          = ''                         ! single ASCII file (or '')
      file_list        = ''                         ! text file of paths (or '')
      file_out         = 'obs_seq.svp',             ! output obs_seq file
      obs_error_sst    = 0.20_r8,                   ! SST error (°C)
      obs_error_vel    = 0.10_r8,                   ! U/V error (m/s)
      avg_obs_per_file = 500000,                    ! pre-allocation hint
      debug            = .true.
   /

.. list-table::
   :header-rows: 1
   :widths: 22 12 23 43

   * - **Namelist item**
     - **Type**
     - **Default**
     - **Description**

   * - ``file_in``
     - character(len=256)
     - ``''`` (empty)
     - Path to a single SVP ASCII input file. Mutually exclusive with ``file_list``.

   * - ``file_list``
     - character(len=256)
     - ``''`` (empty)
     - Text file listing SVP ASCII files, one per line. Use instead of ``file_in`` for batches.

   * - ``file_out``
     - character(len=256)
     - ``'obs_seq.svp'``
     - Name of the DART output observation sequence file. Overwritten if it exists.

   * - ``obs_error_sst``
     - real(r8)
     - ``0.20_r8``
     - Observation error standard deviation for SST (in °C).

   * - ``obs_error_vel``
     - real(r8)
     - ``0.10_r8``
     - Observation error standard deviation for U and V (in m/s).

   * - ``avg_obs_per_file``
     - integer
     - ``500000``
     - Estimated number of valid observations per input file. Used only for pre-allocation.

   * - ``debug``
     - logical
     - ``.true.``
     - If true, prints detailed information for each file and observation.

Output
------
A new observation sequence file named according to `file_out`.  
Each observation record includes location (lat, lon), time, and 
one of the three kinds (SST, U, or V), with error variance 
set to the square of `obs_error_sst` or `obs_error_vel`.

With ``debug = .true.``, the converter prints a summary of accepted/filtered
observations for each input file.

.. code-block:: text

   Output file: obs_seq.svp exists. Replacing it ...
    
   Input file: #1 T_HCSV00_C_BMKG_20251006000000_SVP_300534064005670.csv
        * lat: 0.4730, lon: 43.3212, SST: 26.1200, U: -0.0062, V: -0.0061, date: 2025 Oct 06 00:00:00
    
   Input file: #2 T_HCSV00_C_BMKG_20251006010000_SVP_300534064005670.csv
        * lat: 0.4732, lon: 43.3212, SST: 26.0100, U: 0.0000, V: 0.0061, date: 2025 Oct 06 01:00:00
    
   Input file: #3 T_HCSV00_C_BMKG_20251006020000_SVP_300534064005670.csv
        * lat: 0.4730, lon: 43.3214, SST: 26.0400, U: 0.0062, V: -0.0061, date: 2025 Oct 06 02:00:00
    
   Input file: #4 T_HCSV00_C_BMKG_20251006030000_SVP_300534064005670.csv
        * lat: 0.4732, lon: 43.3214, SST: 26.1400, U: 0.0000, V: 0.0061, date: 2025 Oct 06 03:00:00
    
   Input file: #5 T_HCSV00_C_BMKG_20251006040000_SVP_300534064005670.csv
        * lat: 0.4730, lon: 43.3214, SST: 26.3300, U: 0.0000, V: -0.0061, date: 2025 Oct 06 04:00:00
    
   ...
    
   > Ready to write 72 observations:
     write_obs_seq  opening formatted observation sequence file "obs_seq.svp"
     write_obs_seq  closed observation sequence file "obs_seq.svp"
     svp_to_obs Finished successfully.

Further reading
---------------
For background on the observing platform supported by this converter:

* `Global Drifter Program (GDP): <https://gdp.ucsd.edu/ldl/svpb/>`_ *"SVP drifter overview"* -- NOAA and Scripps-maintained description of the standard SVP and SVP-B drifters used globally.  
* `Lumpkin & Pazos (2007): <https://www.aoml.noaa.gov/phod/docs/LumpkinPazos.pdf>`_ *"Measuring surface currents with Surface Velocity Program drifters"* -- Classic technical paper detailing instrument design, calibration, and data characteristics.  
* `Ocean Observers: <https://www.oceanobservers.org/Learn-about-observing-the-ocean/Drifting-buoys-DBCP>`_ *"Drifting buoys (DBCP)"* -- General introduction to the international drifter network and its coordination under WMO/IOC.
