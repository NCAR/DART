PROGRAM ``fluxnetfull_to_obs``
==============================

Overview
--------

FLUXNET FULLSET data to DART observation sequence file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| This routine is designed to convert Ameriflux FLUXNET data using the FULLSET format into a DART obs_seq file.  
  The data can be downloaded from the `AmeriFlux server. <https://ameriflux.lbl.gov/data/download-data>`__ 
  or `FLUXNET server. <https://fluxnet.org/data/fluxnet2015-data  set/fullset-data-product/>`__  The AmeriFlux network
  is part of `FLUXNET <http://fluxnet.org>`__ and the converter is suitable for either network given the similarity in
  the data format.  This code is not compatible with the SUBSET format of Ameriflux data, however, this code could be
  modified for that purpose.  The FULLSET format contains additional uncertainty information as described in the ONEflux
  processing technique as described by `Pastorello et al., (2020) <https://www.nature.com/articles/s41597-020-0534-3>`__ 
  which is required to formally propogate error into aggregated estimates of tower flux data.  
 
| The AmeriFlux/FLUXNET data products use local time, therefore the code converts the time into UTC such that it is compatible
  with DART and CLM. More information about the AmeriFlux data product is provided below and also within the `flux variable 
  description webpage <https://fluxnet.org/data/fluxnet2015-dataset/fullset-data-product/>`__.

The steps required to prepare Ameriflux data for an assimilation usually include:

#. Download the Ameriflux or FLUXNET FULLSET  data for the towers and years in question (see DATA SOURCES below)
#. Record the TIME ZONE, latitude, longitude, and elevation and tower height at each site. This tower metadata can be found 
   `here <https://fluxnet.org/sites/site-list-and-pages/>`__ or `here <https://ameriflux.lbl.gov/sites/site-search/>`__.
#. Manually provide tower metadata information via the ``fluxnetfull_to_obs_nml`` namelist as this information is
   not contained in the data file itself.
#. Build the DART executables with support for the tower observations. This is done by running ``preprocess`` with
   ``obs_def_tower_mod.f90`` in the list of ``input_files`` for ``preprocess_nml``.
#. Convert each Ameriflux data file individually using ``fluxnetfull_to_obs``
#. If necessary, combine all output files for the region and timeframe of interest into one file using
   :doc:`../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool`


.. Note::
   The Ameriflux data typically have all years combined into one file. However, for some models (CLM, for example), 
   it is required to reorganize the observation sequence files into a series of files that contains ONLY the observations
   for each assimilation time step. This can be achieved with the `makedaily.sh <makedaily.sh>`__ script.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &fluxnetfull_to_obs_nml
      text_input_file = 'textdata.input',
      obs_out_file    = 'obs_seq.out',
      timezoneoffset  = -1,
      latitude        = -1.0,
      longitude       = -1.0,
      elevation       = -1.0,
      flux_height     = -1.0,
      maxgoodqc       = 3,
      gap_filled      = .true.
      energy_balance  = .false.
      time_resolution = 'HR'
      verbose         = .false.
      /

.. container::

   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | Contents        | Type               | Description                                                                 |
   +=================+====================+=============================================================================+
   | text_input_file | character(len=128) | Name of the Ameriflux ASCII file of comma-separated values. This may be a   |
   |                 |                    | relative or absolute filename.                                              |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | obs_out_file    | character(len=128) | Name of the output observation sequence file.                               |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | timezoneoffset  | real               | The local time zone difference (in hours) from UTC of the flux tower site.  |
   |                 |                    | The code will convert from local time to UTC.                               |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | latitude        | real               | Latitude (in degrees N) of the tower.                                       |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | longitude       | real               | Longitude (in degrees E) of the tower. For internal consistency, DART uses  |
   |                 |                    | longitudes in the range [0,360]. An input value of -90 will be converted to |
   |                 |                    | 270, for example.                                                           |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | elevation       | real               | Surface elevation (in meters) of the tower.                                 |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | flux_height     | real               | Height (in meters) of the flux instrument on the tower.                     |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | maxgoodqc       | real               | Maximum value of the observation (obs) quality control (QC) value to pass   |
   |                 |                    | to the output observation sequence file. Obs that have a QC exceeding this  |
   |                 |                    | value will be rejected from the observation sequence file. The 'filter' step|
   |                 |                    | is also capable of discriminating obs based on the QC value, therefore it   |
   |                 |                    | is recommended to keep all obs during the conversion unless you             |
   |                 |                    | are certain they are bad/missing.                                           |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | gap_filled      | logical            | If set to .false. excludes all observations except for measured NEE, SH and |
   |                 |                    | and LE.  Must be set to .true. for time resolutions greater than hourly     |
   |                 |                    | (HH) or half-hourly (HR).                                                   |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | energy_balance  | logical            | If .true. applies energy balance closure to SH and LE fluxes, otherwise     |
   |                 |                    | applies standard approach.                                                  |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | time_resolution | logical            | Defines the flux observation time resolution. Possible values are native    |
   |                 |                    | resolution of hourly (HH) or half-hourly (HR), or coarser gap-filled        |
   |                 |                    | resolutions including daily (DD), weekly (WW) or monthly (MM).              |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | verbose         | logical            | If .true. outputs diagnostic data to help with troubleshooting.             |
   +-----------------+--------------------+-----------------------------------------------------------------------------+

Data sources
------------

| The data was acquired from https://ameriflux.lbl.gov/data/download-data/
  by choosing the data product ``Ameriflux FLUXNET`` and variable set ``FULLSET``.
  The filenames are organized by site and time (HH,HR,WW,DD,MM,YY) and look like:
|  ``AMF_US-Ha1_FLUXNET_FULLSET_HR_1991-2020_3-5.csv``, 
|  ``AMF_US-xBR_FLUXNET_FULLSET_MM_2017-2021_3-5.csv``

| The Ameriflux product in question are ASCII files of comma-separated values taken from all years the data is available.
  The first line is a comma-separated list of column descriptors, and all subsequent lines 
  are comma-separated numerical values. The converter presently searches for the columns pertaining to carbon, water
  and energy fluxes, corresponding quality control fields, uncertainty values, and the time of the observation. Data is available
  at varying time resolutions incuding: native resolution, hourly (HR) or half-hourly (HH), and aggregated resolution, daily (DD),
  weekly (WW), and monthly (MM). The source data does include yearly (YY) time resolution as well, but the coarse nature of yearly flux
  observations poorly constrain fast changing ecological process, thus are not supported by this converter. The required column
  headers depend upon the namelist definitions, including the ``time_resolution``, ``energy_balance`` and ``gap_filled`` settings.  
  These variables are defined as follows:



.. container::

   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | Ameriflux Units | Ameriflux Variable   | Description                   | DART type                | DART kind                   | DART units    |
   +=================+======================+===============================+==========================+=============================+===============+
   |  YYYYMMDDHHMM   | TIMESTAMP_START      | start of time window          | N/A                      | N/A                         | Gregorian     |
   |                 | TIMESTAMP_END        | end of time window            |                          |                             |               |
   |                 |                      | (HH,HR,WW only)               |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   |  YYYYMMDDHHMM   | TIMESTAMP            | time (DD and MM only)         | N/A                      | N/A                         | Gregorian     |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | W/m^2           | LE_F_MDS             | Latent Heat (LE) Flux         | TOWER_LATENT_HEAT_FLUX   | QTY_LATENT_HEAT_FLUX        | W/m^2         |
   |                 |                      | energy_balance = .false.      |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | W/m^2           | LE_RANDUNC           | Uncertainty for LE Flux       | N/A                      | N/A                         | W/m^2         |
   |                 |                      | energy_balance = .false.      |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | [0-3] integer   | LE_F_MDS_QC          | QC for LE Flux                | N/A                      | N/A                         | [0-3] integer |
   |                 |                      | energy_balance = .false.      |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | W/m^2           | LE_CORR              | Latent Heat (LE) Flux         | TOWER_LATENT_HEAT_FLUX   | QTY_LATENT_HEAT_FLUX        | W/m^2         |
   |                 |                      | energy_balance = .true.       |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | W/m^2           | LE_CORR_JOINTUNC     | Uncertainty for LE Flux       | N/A                      | N/A                         | W/m^2         |
   |                 |                      | energy_balance = .true.       |                          |                             |               |
   |                 |                      | Random and Ustar contributions|                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | W/m^2           | H_F_MDS              | Sensible Heat (SH) Flux       | TOWER_SENSIBLE_HEAT_FLUX | QTY_SENSIBLE_HEAT_FLUX      | W/m^2         |
   |                 |                      | energy_balance = .false.      |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | W/m^2           | H_RANDUNC            | Uncertainty for SH Flux       | N/A                      | N/A                         | W/m^2         |
   |                 |                      | energy_balance = .false.      |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | [0-3] integer   | H_F_MDS_QC           | QC for SH Flux                | N/A                      | N/A                         | [0-3] integer |
   |                 |                      | energy_balance = .false.      |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | W/m^2           | H_CORR               | Sensible Heat (SH) Flux       | TOWER_SENSIBLE_HEAT_FLUX | QTY_SENSIBLE_HEAT_FLUX      | W/m^2         |
   |                 |                      | energy_balance = .true.       |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | W/m^2           | H_CORR_JOINTUNC      | Uncertainty for SH Flux       | N/A                      | N/A                         | W/m^2         |
   |                 |                      | energy_balance = .true.       |                          |                             |               |
   |                 |                      | Random and Ustar contributions|                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | umolCO2/m^2/s   | NEE_VUT_REF          | Net Ecosystem Exchange (NEE)  | TOWER_NETC_ECO_EXCHANGE  | QTY_NET_CARBON_PRODUCTION   | gC/m^2/s      |
   |                 |                      | Variable Ustar, reference     |                          |                             |               |
   |                 |                      | flux approach                 |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | umolCO2/m^2/s   | NEE_VUT_REF_JOINTUNC | Uncertainty for NEE Flux      | N/A                      | N/A                         | gC/m^2/s      |
   |                 |                      | Variable Ustar, reference     |                          |                             |               |
   |                 |                      | flux approach. Random and     |                          |                             |               |
   |                 |                      | and Ustar contributions       |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | [0-3] integer   | NEE_VUT_REF_QC       | QC for NEE Flux               | N/A                      | N/A                         | [0-3] integer |
   |                 |                      | Variable Ustar, reference     |                          |                             |               |
   |                 |                      | flux approach.                |                          |                             |               |
   |                 |                      | (HH or HR only)               |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | [0-1] fraction  | NEE_VUT_REF_QC       | QC for NEE Flux               | N/A                      | N/A                         | [0-3] integer |
   |                 |                      | Variable Ustar, reference     |                          |                             |               |
   |                 |                      | flux approach.                |                          |                             |               |
   |                 |                      | (DD, WW and MM only)          |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | umolCO2/m^2/s   | GPP_DT_VUT_REF       | Gross Primary Production (GPP)| TOWER_GPP_FLUX           | QTY_GROSS_PRIMARY_PROD_FLUX | gC/m^2/s      |
   |                 |                      | Day partition, Variable       |                          |                             |               |
   |                 |                      | Ustar, reference approach     |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | umolCO2/m^2/s   | GPP_NT_VUT_REF       | Gross Primary Production (GPP)| TOWER_GPP_FLUX           | QTY_GROSS_PRIMARY_PROD_FLUX | gC/m^2/s      |
   |                 |                      | Night partition, Variable     |                          |                             |               |
   |                 |                      | Ustar, reference approach     |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | umolCO2/m^2/s   | GPP_DT_VUT_16        | 16th percentile uncertainty   |                          |                             |               |
   |                 |                      | estimate for GPP Flux         | N/A                      | N/A                         | gC/m^2/s      |
   |                 |                      | Day partition, Variable       |                          |                             |               |
   |                 |                      | Ustar, reference approach     |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | umolCO2/m^2/s   | GPP_DT_VUT_84        | 84th percentile uncertainty   |                          |                             |               |
   |                 |                      | estimate for GPP Flux         | N/A                      | N/A                         | gC/m^2/s      |
   |                 |                      | Day partition, Variable       |                          |                             |               |
   |                 |                      | Ustar, reference approach     |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | umolCO2/m^2/s   | GPP_NT_VUT_16        | 16th percentile uncertainty   |                          |                             |               |
   |                 |                      | estimate for GPP Flux         | N/A                      | N/A                         | gC/m^2/s      |
   |                 |                      | Night partition, Variable     |                          |                             |               |
   |                 |                      | Ustar, reference approach     |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | umolCO2/m^2/s   | GPP_NT_VUT_84        | 84th percentile uncertainty   |                          |                             |               |
   |                 |                      | estimate for GPP Flux         | N/A                      | N/A                         | gC/m^2/s      |
   |                 |                      | Night partition, Variable     |                          |                             |               |
   |                 |                      | Ustar, reference approach     |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | umolCO2/m^2/s   | RECO_DT_VUT_REF      | Ecosystem Respiration (ER)    | TOWER_ER_FLUX            | QTY_ER_FLUX                 | gC/m^2/s      |
   |                 |                      | Day partition, Variable       |                          |                             |               |
   |                 |                      | Ustar, reference approach     |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | umolCO2/m^2/s   | RECO_NT_VUT_REF      | Ecosystem Respiration (ER)    | TOWER_ER_FLUX            | QTY_ER_FLUX                 | gC/m^2/s      |
   |                 |                      | Night partition, Variable     |                          |                             |               |
   |                 |                      | Ustar, reference approach     |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | umolCO2/m^2/s   | RECO_DT_VUT_16       | 16th percentile uncertainty   |                          |                             |               |
   |                 |                      | estimate for ER Flux          | N/A                      | N/A                         | gC/m^2/s      |
   |                 |                      | Day partition, Variable       |                          |                             |               |
   |                 |                      | Ustar, reference approach     |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | umolCO2/m^2/s   | RECO_DT_VUT_84       | 84th percentile uncertainty   |                          |                             |               |
   |                 |                      | estimate for ER Flux          | N/A                      | N/A                         | gC/m^2/s      |
   |                 |                      | Day partition, Variable       |                          |                             |               |
   |                 |                      | Ustar, reference approach     |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | umolCO2/m^2/s   | RECO_NT_VUT_16       | 16th percentile uncertainty   |                          |                             |               |
   |                 |                      | estimate for ER Flux          | N/A                      | N/A                         | gC/m^2/s      |
   |                 |                      | Night partition, Variable     |                          |                             |               |
   |                 |                      | Ustar, reference approach     |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+
   | umolCO2/m^2/s   | RECO_NT_VUT_84       | 84th percentile uncertainty   |                          |                             |               |
   |                 |                      | estimate for ER Flux          | N/A                      | N/A                         | gC/m^2/s      |
   |                 |                      | Night partition, Variable     |                          |                             |               |
   |                 |                      | Ustar, reference approach     |                          |                             |               |
   +-----------------+----------------------+-------------------------------+--------------------------+-----------------------------+---------------+


Carbon Fluxes
-------------
The flux data files come with several different approaches to estimate the carbon fluxes (NEE, GPP and ER). In the conversion code we choose the
**variable Ustar with reference flux approach (_VUT_REF)**. This choice was based on guidance from Pastorello et al., (2020) which states:

`"The variable proposed in the SUBSET (or FULLSET) product is NEE_VUT_REF since it maintains the temporal variability (as opposed to the MEAN NEE), 
it is representative of the ensemble, and the VUT method is sensitive to possible changes of the canopy (density and height) and site setup,
which can have an impact on the turbulence and consequently on the USTAR threshold. The RECO and GPP products in SUBSET (or FULLSET) are calculated 
from the corresponding NEE variables filtered with the VUT method, generating RECO_NT_VUT_REF and RECO_DT_VUT_REF for ER, and 
GPP_NT_VUT_REF and GPP_DT_VUT_REF for GPP. It is important to use both daytime (DT) and nighttime (NT) variables, and consider their 
difference as uncertainty."`

The reference NEE (_REF) is selected on the basis of model efficiency (MEF) thus is the single NEE data set (out of all ensemble members that sample
the Ustar uncertainty) that is most representative. The mean of the night and day partitioning methods for GPP and ER are used as the observation within
the converter.  


Uncertainty
-----------
Multiple methods are used to estimate uncertainty within the conversion code. There are, in general, three separate sources of uncertainty
in flux data. First, **random uncertainty**, represents the random movement of eddies within the atmosphere where smaller eddies are sampled
more frequently and larger eddies less frequently.  Random uncertainty (_RANDC) estimates are based on `Hollinger, D. Y. & Richardson, A. D. 
Uncertainty in eddy covariance measurements and its application to physiological models. Tree Physiol. 25, 873–885 (2005)`.    
Second, **Ustar uncertainty**, represents the uncertainty contribution from low turbulence conditions as calculated from the Ustar (friction velocity)
threshold.  The ONEflux method uses a bootstrap sampling method to generate a 200 member ensemble from which the flux percentiles are estimated.
The third source of uncertainty, **partitioning uncertainty**, applies to GPP and ER only.  The night (Reichsten et al., 2005) and day (Lasslop et al.,) partitioning methods
estimate both the contributions of photosynthesis (GPP) and ecosystem respiration (ER) as measured from the net carbon exchange (NEE). 

Within the conversion code, the flux uncertainty values (_RANDUNC) account for random uncertainty (SH and LE), where uncertainty denoted with (_JOINTUNC) indicates
combined uncertainty of random and energy balance closure uncertainty (SH and LE). The NEE uncertainty (NEE_VUT_REF_JOINTUNC) accounts for both random and Ustar contributions, 
whereas the GPP and ER uncertainty combine both Ustar and partitioning method uncertainty as:

     | ``Day method Ustar uncertainty (sigma_fluxdt)   = (((fluxDTUNC84-fluxDTUNC16) / 2)^2)^0.5``  
     | ``Night method Ustar uncertainty (sigma_fluxnt) = (((fluxNTUNC84-fluxNTUNC16) / 2)^2)^0.5``

     | ``Ustar and partitioning uncertainty (sigma)  = (0.25 * (sigma_fluxdt)^2 + 0.25 * (sigma_fluxnt)^2)^0.5``  

where ``flux`` stands for either GPP or ER. The 84th and 16th percentile estimates are used to generate 1 sigma estimates
for day and night Ustar uncertainty respectively. Then the contributions of the day and night Ustar uncertainty
are propogated together through standard technqiues assuming gaussian uncertainty distributions (Taylor et al,
An Introduction to Error Analysis).

.. Note::
   
   In practice the relative uncertainty reduces as the flux time resolution increases.  This is likely a result of random uncertainty
   decreasing with coarser time resolutions as the sample size of measurements increases.  This reduced relative
   uncertainty will cause a stronger impact of the observations on the prior model state (i.e. larger increments) during the ``filter``
   step.  To prevent overconfident observations the ONEflux method attempts to account for as many sources of uncertainty as possible
   and that is reflected in this converter code.

   If an observation has a missing uncertainty value, the code estimates an uncertainty based on an empirically-based relative uncertainty
   value that reduces with increasing time resolution.  The default relative uncertainty values are 20%, 10% and 5% for HH/HR, DD/WW and MM
   time resolution respectively.  The user can adjust these values within the source code.



Quality Control
---------------
The general QC naming convention uses an integer system (0-3) defined as the following:

   ``0 = measured``, ``1 = good quality gapfill``, ``2 = medium quality gapfill``, ``3 = poor quality gapfill``. (Refer to Pastorello et al., (2020) or
   Reichstein et al. 2005 Global Change Biology for more information)


The QC values **do not** follow this convention for NEE, SH and LE fluxes for time resolutions coarser than the native resolution of 
HH and HR.  For DD, WW, and MM time resolution observations, the QC value is based on a fraction from 0-1 that indicates the fraction of the time period
that consists of measured or good quality gap-filled data.  Because it is more straightforward to reject observations in DART
based on an integer value scale, the conversion code converts these fractional QC values (0-1) to integer QC values (0-3) as follows:

#. ``QC(integer)=1 when QC(fraction) > 0.90``; 
#. ``QC(integer)=2 when 0.90 >= QC(fraction) >= 0.60``;
#. ``QC(integer)=3 when 0.60 > QC(fraction) >=0``.

.. Note::

   The fraction QC to integer QC conversion approach was based on a qualitative assessment of flux data from Harvard Forest.
   Depending on location and topography not all flux tower data will have a comparably high % of gap filled data.. The user 
   can change these thresholds within the source code.

There are times when a QC value is missing or does not exist for an observation.  In these cases the converter code does the following:

#. Missing QC values (-9999) where the associated flux observation looks physically reasonable are assigned a QC = 2.
#. GPP and ER observations do not have a QC. If the observations are physically reasonable a QC = 2 is assigned.
#. There are situations where +100 is added to an existing QC value such that the observations are purposely rejected
   during the conversion (assuming maxgoodQC = < 100).  These situations are:

      a) When gap_filled = .false., all observations that are not measured (QC = 0)
      b) When GPP or RE give non-physical negative values.
      c) If NEE QC is missing. This is rare.



Gap-Filling
-----------
Gap-filled data is only available for the native resolution format (HH, HR) for flux observations of 
NEE, SH and LE.  For all other situations choosing gap-filled data is mandatory, because measurements of fluxes at
the native resolution are frequently violated that cause the eddy covariance method to fail. These include
situations where the friction velocity (Ustar) falls below a certain value that prevents adequate land-atmosphere mixing, or when
their is instrumentation failure.  In these cases, gap-filling methods (essentially models) are required to calculate daily and coarser 
time resolutions.  Because gap-filled data is technically modeled data (e.g. Marginal Distribution Method (MDS) which relies on
met condtions physcially and temporally similar to missing data) a user may desire only real observations.

In general, we recommend to turn gap-filled data to ``.true.`` during the conversion process, and then use the QC value as a way to discriminate
against lower quality observations during the ``filter`` step.




Energy-Balance
--------------

We provide the option to use LE and SH observations that have been corrected for energy-balance closure.  
In these cases the correction is based on comparing the latent and sensible heat flux against other sources
of energy loss/gain including net incoming radiation and energy radiated through the ground.

Data Policy
-----------
It is important to recognize the flux data providers who have made their research publically available to advance
scientific research.  Please see the Ameriflux data policy `here <https://ameriflux.lbl.gov/data/data-policy/>`__  
and the FLUXNET 2015 data policy provided `here <https://fluxnet.org/data/data-policy/>`__. 


Programs
--------

The ``fluxnetfull_to_obs.f90`` file is the source for the main converter program. Look at the source code where it reads the
example data file. The example code reads each text line into a character buffer and then reads from that buffer to parse up the data items.

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
