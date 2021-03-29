Important capabilities of DART
==============================

In this section we discuss the capabilities of DART that may be of interest to
the user. This is a partial list of all of the functionality that is available
in DART, and additional capabilities and improvements are continually being
added.

As mentioned above, DART allows for both OSSE and OSE systems of models large
and small. This allows users to test both theoretical limits of DA, models, and
observations with idealized experiments as well as to improve actual real-world
forecasts of chaotic systems with real observations.

Models supported by DART
^^^^^^^^^^^^^^^^^^^^^^^^

A full list of models can be found `here <Models.md>`__, but in brief the models
supported by DART include:

============ ============== ================ ==============
Model        Latest version Model            Latest version
============ ============== ================ ==============
lorenz_63    Manhattan      lorenz_84        Manhattan
lorenz_96    Manhattan      lorenz_96_2scale Manhattan
lorenz_04    Manhattan      simple_advection Manhattan
bgrid_solo   Manhattan      WRF              Manhattan
MPAS         Manhattan      ATM              Manhattan
ROMS         Manhattan      CESM             Manhattan
CAM-FV       Manhattan      CAM-CHEM         Manhattan
WACCM        Manhattan      WACCM-X          Manhattan
CICE         Manhattan      CM1              Manhattan
FESOM        Manhattan      NOAH-MP          Manhattan
WRF-Hydro    Manhattan      GCCOM            Lanai
LMDZ         Lanai          MITgcm_ocean     Lanai
NAAPS        Lanai          AM2              Lanai
CAM-SE       Lanai          CLM              Lanai
COAMPS       Lanai          COSMO            Lanai
Dynamo       Lanai          GITM             Lanai
Ikeda        Lanai          JULES            Lanai
MPAS_ocean   Lanai          null_model       Lanai
openggcm     Lanai          PARFLOW          Lanai
sqg          Lanai          TIE-GCM          Lanai
WRF-CHEM     Lanai          ECHAM            Prior to Lanai
PBL_1d       Prior to Lanai MITgcm_annulus   Prior to Lanai
forced_barot Prior to Lanai pe2lyr           Prior to Lanai
ROSE         Prior to Lanai CABLE            Prior to Lanai
============ ============== ================ ==============

The models listed as “Prior to Lanai” will take some additional work to
integrate with a supported version of DART; please contact the dart @ ucar.edu
team for more information. The versions listed as “Lanai” will be ported to the
Manhattan version of DART depending on the needs of the user community as well
as the availablity of resources on the DART team.

` <#obsConverters>`__

Observation converters provided by DART
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given a way to compute the expected observation value from the model state, in
theory any and all observations can be assimilated by DART through the
``obs_seq.out`` file. In practice this means a user-defined observation
converter is required. DART provides many observation converters to make this
process easier for the user. Under the directory
``DARTHOME/observations/obs_converters`` there are multiple subdirectories, each
of which has at least one observation converter. The list of these directories
is as follows:

+----------------------------+-------------------+----------------------------+
| Observation                | Directory         | Format                     |
+============================+===================+============================+
| `Atmospheric Infrared      | AIRS              | HDF-EOS                    |
| Sounder <http              |                   |                            |
| s://airs.jpl.nasa.gov/>`__ |                   |                            |
| satellite retrievals       |                   |                            |
+----------------------------+-------------------+----------------------------+
| `Avi                       | Aviso             | netCDF                     |
| so <https://www.aviso.alti |                   |                            |
| metry.fr/en/home.html>`__: |                   |                            |
| satellite derived sea      |                   |                            |
| surface height             |                   |                            |
+----------------------------+-------------------+----------------------------+
| Level 4 Flux Tower data    | Ameriflux         | Comma-separated text       |
| from                       |                   |                            |
| `AmeriFlux <ht             |                   |                            |
| tp://ameriflux.lbl.gov>`__ |                   |                            |
+----------------------------+-------------------+----------------------------+
| Level 2 soil moisture from | COSMOS            | Fixed-width text           |
| `COSMOS <http://c          |                   |                            |
| osmos.hwr.arizona.edu/>`__ |                   |                            |
+----------------------------+-------------------+----------------------------+
| Doppler wind lidar         | DWL               | ASCII text                 |
+----------------------------+-------------------+----------------------------+
| GPS retrievals of          | GPSPW             | netCDF                     |
| precipitable water         |                   |                            |
+----------------------------+-------------------+----------------------------+
| GSI observation file       | GSI2DART          | Fortran binary             |
+----------------------------+-------------------+----------------------------+
| Global                     | GTSPP             | netCDF                     |
| Temperature-Salinity       |                   |                            |
| Profile Program            |                   |                            |
| (`G                        |                   |                            |
| TSPP <http://www.nodc.noaa |                   |                            |
| .gov/GTSPP/index.html>`__) |                   |                            |
+----------------------------+-------------------+----------------------------+
| Meteorological             | MADIS             | netCDF                     |
| Assimilation Data Ingest   |                   |                            |
| System                     |                   |                            |
| (`MADIS <h                 |                   |                            |
| ttp://madis.noaa.gov/>`__) |                   |                            |
+----------------------------+-------------------+----------------------------+
| `                          | MIDAS             | netCDF                     |
| MIDAS <https://www.science |                   |                            |
| direct.com/science/article |                   |                            |
| /pii/S0273117712001135>`__ |                   |                            |
| ionospheric obs            |                   |                            |
+----------------------------+-------------------+----------------------------+
| `MODIS <https:             | MODIS             | Comma-separated text       |
| //modis.gsfc.nasa.gov/>`__ |                   |                            |
| satellite retrievals       |                   |                            |
+----------------------------+-------------------+----------------------------+
| `NCEP                      | NCEP/prep_bufr    | PREPBUFR                   |
| PREPBUFR <h                |                   |                            |
| ttps://www.emc.ncep.noaa.g |                   |                            |
| ov/mmb/data_processing/pre |                   |                            |
| pbufr.doc/document.htm>`__ |                   |                            |
+----------------------------+-------------------+----------------------------+
| NCEP ASCII observations    | NCEP/ascii_to_obs | NCEP text files            |
+----------------------------+-------------------+----------------------------+
| `ROMS <h                   | ROMS              | netCDF                     |
| ttps://www.myroms.org/>`__ |                   |                            |
| verification observations  |                   |                            |
+----------------------------+-------------------+----------------------------+
| Satellite winds from       | SSEC              | ASCII text                 |
| `SSEC <https://w           |                   |                            |
| ww.ssec.wisc.edu/data/>`__ |                   |                            |
+----------------------------+-------------------+----------------------------+
| Sea surface temperature    | SST               | netCDF                     |
+----------------------------+-------------------+----------------------------+
| Special Sensor Ultraviolet | SSUSI             | netCDF                     |
| Spectrographic Imager      |                   |                            |
| (`SSUSI <http              |                   |                            |
| s://ssusi.jhuapl.edu/>`__) |                   |                            |
| retrievals                 |                   |                            |
+----------------------------+-------------------+----------------------------+
| World Ocean Database       | WOD               | World Ocean Database       |
| (`WOD <h                   |                   | packed ASCII               |
| ttp://www.nodc.noaa.gov/OC |                   |                            |
| 5/WOD09/pr_wod09.html>`__) |                   |                            |
+----------------------------+-------------------+----------------------------+
| `National Snow and Ice     | cice              | Binary sea ice             |
| Data                       |                   |                            |
| Cen                        |                   |                            |
| ter <http://nsidc.org/>`__ |                   |                            |
| sea ice obs                |                   |                            |
+----------------------------+-------------------+----------------------------+
| VTEC                       | gnd_gps_vtec      | ASCII text                 |
| `                          |                   |                            |
| Madrigal <http://millstone |                   |                            |
| hill.haystack.mit.edu/>`__ |                   |                            |
| upper atmospheric obs      |                   |                            |
+----------------------------+-------------------+----------------------------+
| GPS obs from               | gps               | netCDF                     |
| `COSMIC <http              |                   |                            |
| ://www.cosmic.ucar.edu>`__ |                   |                            |
+----------------------------+-------------------+----------------------------+
| Oklahoma                   | ok_mesonet        | Oklahoma Mesonet MDF files |
| `Mesonet <h                |                   |                            |
| ttp://www.mesonet.org/>`__ |                   |                            |
| MDF obs                    |                   |                            |
+----------------------------+-------------------+----------------------------+
| `QuikSCAT <http:           | quikscat          | HDF 4                      |
| //winds.jpl.nasa.gov/missi |                   |                            |
| ons/quikscat/index.cfm>`__ |                   |                            |
| scatterometer winds        |                   |                            |
+----------------------------+-------------------+----------------------------+
| Radar reflectivity/radial  | Radar             | WSR-88D (NEXRAD)           |
| velocity obs               |                   |                            |
+----------------------------+-------------------+----------------------------+
| `MODIS Snowcover           | snow              | General text               |
| Fraction <htt              |                   |                            |
| ps://modis.gsfc.nasa.gov/d |                   |                            |
| ata/dataprod/mod10.php>`__ |                   |                            |
| obs                        |                   |                            |
+----------------------------+-------------------+----------------------------+
| Text file                  | Text              | General text               |
| (e.g. spreadsheet) obs     |                   |                            |
+----------------------------+-------------------+----------------------------+
| Total precipitable water   | tpw               | HDF-EOS                    |
| from AQUA                  |                   |                            |
+----------------------------+-------------------+----------------------------+
| Automated Tropical Cyclone | Tropical Cyclones | Fixed width text           |
| Forecast                   |                   |                            |
| (`ATCF <https://www.nrlm   |                   |                            |
| ry.navy.mil/atcf_web/>`__) |                   |                            |
| obs                        |                   |                            |
+----------------------------+-------------------+----------------------------+
| `LITT                      | var               | little-r                   |
| LE_R <http://www2.mmm.ucar |                   |                            |
| .edu/mm5/On-Line-Tutorial/ |                   |                            |
| little_r/little_r.html>`__ |                   |                            |
| obs                        |                   |                            |
+----------------------------+-------------------+----------------------------+
| `MM5 <http://              | var               | MM5 3D-VAR 2.0 Radar data  |
| www2.mmm.ucar.edu/mm5/>`__ |                   | files                      |
| 3D-VAR radar obs           |                   |                            |
+----------------------------+-------------------+----------------------------+

Data assimilation algorithms available in DART
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DART allows users to test the impact of using multiple different types of
algorithms for filtering, inflation/deflation, and covariance localization.

DART offers numerous **filter algorithms**. These determine how the posterior
distribution is updated based on the observations and the prior ensemble. The
following table lists the filters supported in DART along with their type (set
by *filter_kind* in ``input.nml`` under the “assim_tools_nml” section):

+--------------------+----------------------------+----------------------------+
| Filter #           | Filter Name                | References                 |
+====================+============================+============================+
| 1                  | EAKF (Ensemble Adjustment  | **Anderson, J. L.**, 2001. |
|                    | Kalman Filter)             | An Ensemble Adjustment     |
|                    |                            | Kalman Filter for Data     |
|                    |                            | Assimilation. *Monthly     |
|                    |                            | Weather Review*, 129,      |
|                    |                            | 2884-2903.                 |
|                    |                            | https://doi.org/10.1175/   |
|                    |                            | 1520-0493%282001%29129%3C2 |
|                    |                            | 884%3AAEAKFF%3E2.0.CO%3B2. |
|                    |                            | \ **Anderson, J. L.**,     |
|                    |                            | 2003. A local least        |
|                    |                            | squares framework for      |
|                    |                            | ensemble filtering.        |
|                    |                            | *Monthly Weather Review*,  |
|                    |                            | 131, 634-642.              |
|                    |                            | https://do                 |
|                    |                            | i.org/10.1175/1520-0493%28 |
|                    |                            | 2003%29131%3C0634%3AALLSFF |
|                    |                            | %3E2.0.CO%3B2\ **Anderson, |
|                    |                            | J., Collins, N.**, 2007:   |
|                    |                            | Scalable Implementations   |
|                    |                            | of Ensemble Filter         |
|                    |                            | Algorithms for Data        |
|                    |                            | Assimilation. *Journal of  |
|                    |                            | Atmospheric and Oceanic    |
|                    |                            | Technology*, 24,           |
|                    |                            | 1452-1463.                 |
|                    |                            | https://d                  |
|                    |                            | oi.org/10.1175/JTECH2049.1 |
+--------------------+----------------------------+----------------------------+
| 2                  | ENKF (Ensemble Kalman      | **Evensen, G.**, 2003. The |
|                    | Filter)                    | Ensemble Kalman Filter:    |
|                    |                            | Theoretical Formulation    |
|                    |                            | and Practical              |
|                    |                            | Implementation. *Ocean     |
|                    |                            | Dynamics*. 53(4), 343–367. |
|                    |                            | https://doi.org/1          |
|                    |                            | 0.1007%2Fs10236-003-0036-9 |
+--------------------+----------------------------+----------------------------+
| 3                  | Kernel filter              |                            |
+--------------------+----------------------------+----------------------------+
| 4                  | Observation Space Particle |                            |
|                    | filter                     |                            |
+--------------------+----------------------------+----------------------------+
| 5                  | Random draw from posterior | None. :exclamation:        |
|                    |                            | *IMPORTANT*: (contact dart |
|                    |                            | @ ucar.edu before using)   |
+--------------------+----------------------------+----------------------------+
| 6                  | Deterministic draw from    | None. :exclamation:        |
|                    | posterior with fixed       | *IMPORTANT*: (contact dart |
|                    | kurtosis                   | @ ucar.edu before using)   |
+--------------------+----------------------------+----------------------------+
| 7                  | Boxcar kernel filter       |                            |
+--------------------+----------------------------+----------------------------+
| 8                  | Rank Histogram filter      | **Anderson, J. L.,** 2010. |
|                    |                            | A Non-Gaussian Ensemble    |
|                    |                            | Filter Update for Data     |
|                    |                            | Assimilation. *Monthly     |
|                    |                            | Weather Review*, 139,      |
|                    |                            | 4186-4198.                 |
|                    |                            | https://doi                |
|                    |                            | .org/10.1175/2010MWR3253.1 |
+--------------------+----------------------------+----------------------------+
| 9                  | Particle filter            | **Poterjoy, J.**, 2016. A  |
|                    |                            | localized particle filter  |
|                    |                            | for high-dimensional       |
|                    |                            | nonlinear systems.         |
|                    |                            | *Monthly Weather Review*,  |
|                    |                            | 144 59-76.                 |
|                    |                            | https://doi.o              |
|                    |                            | rg/10.1175/MWR-D-15-0163.1 |
+--------------------+----------------------------+----------------------------+

DART also has several **inflation algorithms** available for both prior (the
first value in the namelist) and posterior (the second value in the namelist).
The following table lists the inflation “flavors” supported in DART along with
their type number (set by *inf_flavor* in ``input.nml`` under the “filter_nml”
section):

+--------------------+----------------------------+----------------------------+
| Flavor #           | Inflation flavor name      | References                 |
+====================+============================+============================+
| 0                  | No inflation               | n/a                        |
+--------------------+----------------------------+----------------------------+
| 1                  | (Not Supported)            | n/a                        |
+--------------------+----------------------------+----------------------------+
| 2                  | Spatially-varying          | **Anderson, J. L.**, 2009. |
|                    | state-space (Gaussian)     | Spatially and temporally   |
|                    |                            | varying adaptive           |
|                    |                            | covariance inflation for   |
|                    |                            | ensemble filters. *Tellus  |
|                    |                            | A*, 61, 72-83,             |
|                    |                            | https://doi.org/10.111     |
|                    |                            | 1/j.1600-0870.2008.00361.x |
+--------------------+----------------------------+----------------------------+
| 3                  | Spatially-fixed            | **Anderson, J. L.**, 2007. |
|                    | state-space (Gaussian)     | An adaptive covariance     |
|                    |                            | inflation error correction |
|                    |                            | algorithm for ensemble     |
|                    |                            | filters. *Tellus A*, 59,   |
|                    |                            | 210-224,                   |
|                    |                            | https://doi.org/10.111     |
|                    |                            | 1/j.1600-0870.2006.00216.x |
+--------------------+----------------------------+----------------------------+
| 4                  | Relaxation to prior spread | **Whitaker, J.S. and T.M.  |
|                    | (posterior inflation only) | Hamill**, 2012. Evaluating |
|                    |                            | Methods to Account for     |
|                    |                            | System Errors in Ensemble  |
|                    |                            | Data Assimilation.         |
|                    |                            | *Monthly Weather Review*,  |
|                    |                            | 140, 3078–3089,            |
|                    |                            | https://doi.or             |
|                    |                            | g/10.1175/MWR-D-11-00276.1 |
+--------------------+----------------------------+----------------------------+
| 5                  | Enhanced spatially-varying | **El Gharamti M.**, 2018.  |
|                    | state-space (inverse       | Enhanced Adaptive          |
|                    | gamma)                     | Inflation Algorithm for    |
|                    |                            | Ensemble Filters. *Monthly |
|                    |                            | Weather Review*, 2,        |
|                    |                            | 623-640,                   |
|                    |                            | https://doi.o              |
|                    |                            | rg/10.1175/MWR-D-17-0187.1 |
+--------------------+----------------------------+----------------------------+

DART also offers the ability to correct for sampling errors. DART’s localization
and sampling error correction algorithm is described in > > **Anderson, J.L.**,
2012. Localization and Sampling Error Correction in Ensemble Kalman Filter Data
Assimilation. *Monthly Weather Review*, 140, 2359–2371.
https://doi.org/10.1175/MWR-D-11-00013.1

This behavior can be turned on or off via the *sampling_error_correction* in
``input.nml`` under the “assim_tools_nml” section. The following covariance
localization options are available (set by *select_localization* in
``input.nml`` under the “cov_cutoff_nml” section):

+--------------------+----------------------------+----------------------------+
| Loc #              | Localization type          | References                 |
+====================+============================+============================+
| 1                  | Gaspari-Cohn eq. 4.10      | **Gaspari, G. and Cohn, S. |
|                    |                            | E.**, 1999. Construction   |
|                    |                            | of correlation functions   |
|                    |                            | in two and three           |
|                    |                            | dimensions. *Quarterly     |
|                    |                            | Journal of the Royal       |
|                    |                            | Meteorological Society*,   |
|                    |                            | 125, 723-757.              |
|                    |                            | https://doi.               |
|                    |                            | org/10.1002/qj.49712555417 |
+--------------------+----------------------------+----------------------------+
| 2                  | Boxcar                     | None                       |
+--------------------+----------------------------+----------------------------+
| 3                  | Ramped boxcar              | None                       |
+--------------------+----------------------------+----------------------------+

The following image depicts all three of these options:

|cutoff_fig|

.. |cutoff_fig| image:: images/cutoff_fig.png
   :width: 100%