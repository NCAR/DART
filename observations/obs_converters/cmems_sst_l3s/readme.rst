.. _cmems_sst_to_obs:

CMEMS L3S SST observations
==========================

.. contents:: 
   :depth: 3
   :local:

Overview
--------
This converter reads Sea Surface Temperature (SST) observations from the
Copernicus Marine Environment Monitoring Service (CMEMS) Level-3S (L3S)
blended SST product and converts them into a DART ``obs_seq`` file.

The input data consist of bias-corrected, multi-sensor SST observations
from the ODYSSEA system. These products are distributed in netCDF format
and include per-pixel uncertainty estimates and quality control flags.

The converter extracts:

- SST observations (``adjusted_sea_surface_temperature``)
- Observation uncertainties (``sses_standard_deviation``)
- Quality flags (``quality_level``)
- Observation time and location

and writes them into a DART observation sequence file using the
``SATELLITE_BLENDED_SST`` observation type. The unit of 
temperature is C.  

For more information about the SST product, visit the 
Copernicus Marine Service: https://marine.copernicus.eu

Input data
----------
The converter expects netCDF files from the CMEMS SST L3S product, e.g., ``SST_GLO_PHY_L3S_*.nc``

Required variables:

- ``adjusted_sea_surface_temperature`` (Kelvin)
- ``sses_standard_deviation`` (Kelvin)
- ``quality_level`` (0–5)
- ``latitude``
- ``longitude``
- ``time``

Only observations with quality level greater than 3 are retained.

Obs Error handling
------------------
The input uncertainty field ``sses_standard_deviation`` is used as the
observation error.

To ensure numerical stability and prevent unrealistically small errors
from dominating the assimilation, a lower bound is applied:

::

   err = max(sses_standard_deviation, OBS_ERROR_SD_MIN)

with default:

::

   OBS_ERROR_SD_MIN = 0.05 °C

NaN values in the uncertainty field are replaced with this minimum value.
This lower bound represents a practical safeguard against over-weighting
individual observations and does not replace more sophisticated error
modeling or inflation during assimilation.

Namelist
--------

The converter uses the following namelist:

::

   &cmems_sst_to_obs_nml
      file_list         = 'sst_file_list.txt'  ! text file listing input netCDF files
      file_out          = 'obs_seq.sst'        ! output obs_seq filename
      avg_obs_per_file  = 500000               ! pre-allocation hint
      debug             = .true.               ! print diagnostic output
   /

Usage
-----
1. Create a file list:

::

   ls *.nc > sst_file_list.txt

2. Edit ``input.nml`` to configure options

3. Build the converter:

::

   cd work
   ./quickbuild.sh

4. Run:

::

   ./cmems_sst_to_obs

5. Output:

::

   obs_seq.sst
