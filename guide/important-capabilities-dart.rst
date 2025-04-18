.. _important-capabilities-dart:

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

A full list of models can be found :doc:`here <../models/README>`, but in brief the models
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
CAM-SE       Manhattan      CLM              Manhattan
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


Observation converters provided by DART
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given a way to compute the expected observation value from the model state, in
theory any and all observations can be assimilated by DART through the
``obs_seq.out`` file. In practice this means a user-defined observation
converter is required. DART provides many observation converters to make this
process easier for the user. Under the directory
``DART/observations/obs_converters`` there are multiple subdirectories, each
of which has at least one observation converter. The list of these directories
is as follows:



+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Observation                                                                                          | Directory         | Format                            |
+======================================================================================================+===================+===================================+
| `Atmospheric Infrared Sounder <https://airs.jpl.nasa.gov/>`__ satellite retrievals                   | AIRS              | HDF-EOS                           |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| `Advanced Microwave Sounding Unit <https://aqua.nasa.gov/content/amsu>`__ brightness temperatures    | AIRS              | netCDF                            |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| `Aviso <https://www.aviso.altimetry.fr/en/home.html>`__: satellite derived sea surface height        | Aviso             | netCDF                            |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Level 4 Flux Tower data                                                                              | Ameriflux         | Comma-separated text              |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Ameriflux Fullset Flux Tower data from `AmeriFlux <https://ameriflux.lbl.gov/data/download-data>`__  | Ameriflux         | Comma-separated text              |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Level 2 soil moisture from `COSMOS <http://cosmos.hwr.arizona.edu/>`__                               | COSMOS            | Fixed-width text                  |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Doppler wind lidar                                                                                   | DWL               | ASCII text                        |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| GPS retrievals of precipitable water                                                                 | GPSPW             | netCDF                            |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| GSI observation file                                                                                 | GSI2DART          | Fortran binary                    |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Global Temperature-Salinity Profile Program (`GTSPP <http://www.nodc.noaa.gov/GTSPP/index.html>`__)  | GTSPP             | netCDF                            |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Meteorological Assimilation Data Ingest System (`MADIS <http://madis.noaa.gov/>`__)                  | MADIS             | netCDF                            |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| `MIDAS <https://www.sciencedirect.com/science/article/pii/S0273117712001135>`__ ionospheric obs      | MIDAS             | netCDF                            |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| `MODIS <https://modis.gsfc.nasa.gov/>`__ satellite retrievals                                        | MODIS             | Comma-separated text              |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| `NCEP PREPBUFR <https://www.emc.ncep.noaa.gov/mmb/data_processing/prepbufr.doc/document.htm>`__      | NCEP/prep_bufr    | PREPBUFR                          |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| NCEP ASCII observations                                                                              | NCEP/ascii_to_obs | NCEP text files                   |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| `ROMS <https://www.myroms.org/>`__ verification observations                                         | ROMS              | netCDF                            |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Satellite winds from `SSEC <https://www.ssec.wisc.edu/data/>`__                                      | SSEC              | ASCII text                        |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Sea surface temperature                                                                              | SST               | netCDF                            |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Solar-Induced Fluorescence                                                                           | SIF               | netCDF                            |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Special Sensor Ultraviolet Spectrographic Imager `(SSUSI) <https://ssusi.jhuapl.edu/>`__ retrievals  | SSUSI             | netCDF                            |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| World Ocean Database `(WOD) <http://www.nodc.noaa.gov/OC5/WOD09/pr_wod09.html>`__                    | WOD               | World Ocean Database packed ASCII |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| `National Snow and Ice Data Center <http://nsidc.org/>`__ sea ice obs                                | cice              | Binary sea ice                    |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| VTEC `Madrigal <http://millstone hill.haystack.mit.edu/>`__ upper atmospheric obs                    | gnd_gps_vtec      | ASCII text                        |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| GPS obs from `COSMIC <http://www.cosmic.ucar.edu/>`__                                                | gps               | netCDF                            |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Oklahoma `Mesonet <http://www.mesonet.org/>`__ MDF obs                                               | ok_mesonet        | Oklahoma Mesonet MDF files        |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| `QuikSCAT <http://winds.jpl.nasa.gov/missions/quikscat/index.cfm>`__ scatterometer winds             | quikscat          | HDF 4                             |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Radar reflectivity/radial velocity obs                                                               | Radar             | WSR-88D (NEXRAD)                  |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| `MODIS <https://modis.gsfc.nasa.gov/data/dataprod/mod10.php>`__ Snowcover Fraction obs               | snow              | General text                      |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Text file (e.g. spreadsheet) obs                                                                     | Text              | General text                      |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Total precipitable water from AQUA                                                                   | tpw               | HDF-EOS                           |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| Automated Tropical Cyclone Forecast (`ATCF <https://www.nrlmry.navy.mil/atcf_web/>`__) obs           | Tropical Cyclones | Fixed width text                  |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| `LITTLE_R <http://www2.mmm.ucar.edu/mm5/On-Line-Tutorial/little_r/little_r.html>`__ obs              | var               | little-r                          |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+
| `MM5 <http://www2.mmm.ucar.edu/mm5/>`__ 3D-VAR radar obs                                             | var               | MM5 3D-VAR 2.0 Radar data files   |
+------------------------------------------------------------------------------------------------------+-------------------+-----------------------------------+


Data assimilation algorithms available in DART
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DART allows users to test the impact of using multiple different types of
algorithms for filtering, inflation/deflation, and covariance localization.

DART offers numerous **filter algorithms** in the :ref:`QCEFF`. These determine how the prior
distribution is updated based on the observations. The
following table lists the filters supported in DART. For other types of filters,
such as the particle filter (PF), please contact the DART team.

+--------------------+----------------------------+--------------------------------------------+
| Filter #           | Filter Name                | References                                 |
+====================+============================+============================================+
| 1                  | EAKF (Ensemble Adjustment  | **Anderson, J. L.**, 2001. [1]_            |
|                    | Kalman Filter)             | **Anderson, J. L.**, 2003. [2]_            |
+--------------------+----------------------------+--------------------------------------------+
| 2                  | EnKF (Ensemble Kalman      | **Evensen, G.**, 2003. [3]_                |
|                    | Filter)                    |                                            |
+--------------------+----------------------------+--------------------------------------------+
| 3                  | Rank Histogram filter      | **Anderson, J. L.,** 2010. [4]_            |
|                    | (Unbounded)                |                                            |
+--------------------+----------------------------+--------------------------------------------+
| 4                  | Gamma filter               | **Anderson, J. L.** 2022. [11]_            |
|                    | (aka GIGG filter)          | **Bishop, C. H.** 2016. [15]_              |
+--------------------+----------------------------+--------------------------------------------+
| 5                  | Bounded Normal Rank        | **Anderson, J. L.** 2023. [12]_            |
|                    | Histogram filter (BNRHF)   | **Anderson, J. L.** 2024. [13]_            |
+--------------------+----------------------------+--------------------------------------------+
| 6                  | Kernel Density Estimation  | **Grooms, I. and C. Riedel** 2024. [14]_   |
|                    | (KDE) filter               |                                            |
+--------------------+----------------------------+--------------------------------------------+

DART also has several **inflation algorithms** (i.e. "flavors") available for both prior 
(``inf_flavor``, first value) and posterior (``inf_flavor``, second value) inflation. The
``inf_flavor`` setting is located under the ``&filter_nml`` heading within ``input.nml``
The following table lists the inflation “flavors” supported in DART along with
their integer value:

+-------------+-----------------------------+----------------------------------+
| inf_flavor #| Inflation flavor name       | References                       |
+=============+=============================+==================================+
| 0           | No inflation                | n/a                              |
+-------------+-----------------------------+----------------------------------+
| 1           | (Not Supported)             | n/a                              |
+-------------+-----------------------------+----------------------------------+
| 2           | Spatially-varying           | **Anderson, J. L.**, 2009. [6]_  |
|             | state-space (Gaussian)      |                                  |
+-------------+-----------------------------+----------------------------------+
| 3           | Spatially-fixed             | **Anderson, J. L.**, 2007. [5]_  |
|             | state-space (Gaussian)      |                                  |
+-------------+-----------------------------+----------------------------------+
| 4           | Relaxation to prior spread  | **Whitaker, J. S.**              |
|             | (posterior inflation only)  | **and T. M. Hamill**, 2012. [7]_ |
+-------------+-----------------------------+----------------------------------+
| 5           | Enhanced spatially-varying  | **El Gharamti M.**, 2018. [8]_   |
|             | state-space (inverse gamma) |                                  |
+-------------+-----------------------------+----------------------------------+

DART has the ability to correct for sampling errors in the regression
caused by finite ensemble sizes. DART’s sampling error correction algorithm
(and localization algorithm) is described in **Anderson, J. L.**, 2012 [9]_
Sampling error correction can be turned on or off via the *sampling_error_correction*
variable in the ``input.nml`` under the “assim_tools_nml” section.

The following covariance localization options are available
(set by *select_localization* in ``input.nml`` under the “cov_cutoff_nml” section):

+--------+----------------------------+----------------------------------+
| Loc #  | Localization type          | References                       |
+========+============================+==================================+
| 1      | Gaspari-Cohn eq. 4.10      | **Gaspari, G.**                  |
|        |                            | **and S. E. Cohn**, 1999. [10]_  |
+--------+----------------------------+----------------------------------+
| 2      | Boxcar                     | None                             |
+--------+----------------------------+----------------------------------+
| 3      | Ramped boxcar              | None                             |
+--------+----------------------------+----------------------------------+

The following image depicts all three of these options:

|cutoff_fig|

.. |cutoff_fig| image:: images/loc_types.png
   :width: 100%


References
----------

.. [1] Anderson, J. L., 2001:
       An Ensemble Adjustment Kalman Filter for Data Assimilation.
       *Monthly Weather Review*, **129**, 2884-2903.
       `doi:10.1175/1520-0493(2001)129<2884:AEAKFF>2.0.CO;2 <https://doi.org/10.1175/1520-0493(2001)129\<2884:AEAKFF\>2.0.CO;2>`__

.. [2] Anderson, J. L., 2003:
       A local least squares framework for ensemble filtering.
       *Monthly Weather Review*, **131**, 634-642.
       `doi:10.1175/1520-0493(2003)131<0634:ALLSFF>2.0.CO;2 <https://doi.org/10.1175/1520-0493(2003)131\<0634:ALLSFF\>2.0.CO;2>`__

.. [3] Evensen, G., 2003:
       The Ensemble Kalman Filter: Theoretical Formulation and Practical Implementation.
       *Ocean Dynamics*. **53(4)**, 343–367.
       `doi:10.1007%2Fs10236-003-0036-9 <https://doi.org/10.1007%2Fs10236-003-0036-9>`__

.. [4] Anderson, J. L., 2010:
       A Non-Gaussian Ensemble Filter Update for Data Assimilation.
       *Monthly Weather Review*, **139**, 4186-4198.
       `doi:10.1175/2010MWR3253.1 <https://doi.org/10.1175/2010MWR3253.1>`__

.. [5] Anderson, J. L., 2007:
       An adaptive covariance inflation error correction algorithm for ensemble filters.
       *Tellus A*, **59**, 210-224,
       `doi:10.1111/j.1600-0870.2006.00216.x <https://doi.org/10.1111/j.1600-0870.2006.00216.x>`__

.. [6] Anderson, J. L., 2009:
       Spatially and temporally varying adaptive covariance inflation for ensemble filters.
       *Tellus A*, **61**, 72-83,
       `doi:10.1111/j.1600-0870.2008.00361.x <https://onlinelibrary.wiley.com/doi/10.1111/j.1600-0870.2008.00361.x>`__

.. [7] Whitaker, J. S. and T. M.  Hamill, 2012:
       Evaluating Methods to Account for System Errors in Ensemble Data Assimilation.
       *Monthly Weather Review*, **140**, 3078–3089,
       `doi:10.1175/MWR-D-11-00276.1 <https://doi.org/10.1175/MWR-D-11-00276.1>`__

.. [8] El Gharamti M., 2018:
       Enhanced Adaptive Inflation Algorithm for Ensemble Filters.
       *Monthly Weather Review*, **2**, 623-640,
       `doi:10.1175/MWR-D-17-0187.1 <https://doi.org/10.1175/MWR-D-17-0187.1>`__

.. [9] Anderson, J. L., 2012:
       Localization and Sampling Error Correction in Ensemble Kalman Filter Data Assimilation.
       *Monthly Weather Review*, 140, 2359–2371.
       `doi:10.1175/MWR-D-11-00013.1 <https://doi.org/10.1175/MWR-D-11-00013.1>`__

.. [10] Gaspari, G. and S. E. Cohn, 1999:
        Construction of correlation functions in two and three dimensions.
        *Quarterly Journal of the Royal Meteorological Society*, **125**, 723-757.
        `doi:10.1002/qj.49712555417 <https://doi.org/10.1002/qj.49712555417>`__

.. [11] Anderson, J. L., 2022:
        A quantile-conserving ensemble filter framework. Part I: Updating an observed variable.
        *Monthly Weather Review*, **150**, 1061-1074.
        `doi:10.1175/MWR-D-21-0229.1 <https://doi.org/10.1175/MWR-D-21-0229.1>`__

.. [12] Anderson, J. L., 2023:
        A quantile-conserving ensemble filter framework. Part II: Regression of
        observation increments in a Probit and probability integral transformed space.
        *Monthly Weather Review*, **151**, 2759–2777.
        `doi:10.1175/MWR-D-23-0065.1 <https://doi.org/10.1175/MWR-D-23-0065.1>`__

.. [13] Anderson, J. L., C. Riedel, M. Wieringa, F. Ishraque, M. Smith and H. Kershaw, 2024:
        A quantile-conserving ensemble filter framework. Part III: Data assimilation
        for mixed distributions with application to a low-order tracer advection model.
        *Monthly Weather Review*, **152**, 2111–2127.
        `doi:10.1175/MWR-D-23-0255.1 <https://doi.org/10.1175/MWR-D-23-0255.1>`__

.. [14] Grooms, I. and C. Riedel, 2024:
        A quantile-conserving ensemble filter based on kernel-density estimation.
        *Remote Sensing*, **16**, 2377.
        `doi:10.3390/rs16132377 <https://doi.org/10.3390/rs16132377>`__

.. [15] Bishop, C. H., 2016:
        The GIGG-EnKF: ensemble Kalman filtering for highly skewed non-negative
        uncertainty distributions.
        *Quarterly Journal of the Royal Meteorological Society*, **142**, 1395-1412.
        `doi:10.1002/qj.2742 <https://doi.org/10.1002/qj.2742>`__
