Introduction to DART’s support for RTTOV
========================================

DART supports satellite radiance assimilation through an interface to 
the radiative transfer model RTTOV. 
RTTOV is a fast radiative transfer model that is widely used in research 
and operations. RTTOV contains observation operators for visible/infrared and microwave 
radiances/brightness temperatures. Observations from a wide range of satellites 
(e.g. GOES, FY, METOP, ...) and sensors (e.g. ABI, AMSU-A, SEVIRI, ...) are supported 
(`see the complete list here <https://nwp-saf.eumetsat.int/site/software/rttov/documentation/platforms-supported/>`__).
For more detail on RTTOV see the `RTTOV user guide <https://www.nwpsaf.eu/site/software/rttov/documentation/>`__.

This documentation describes 
1. Compilation and setup
2. High-level workflow
3. Input data to RTTOV
4. Tips for the assimilation of visible/infrared radiances  
5. Converting real observations to DART format
6. Current list of known issues


The run-time behavior of RTTOV is mostly controlled by the 
``&obs_def_rttov_nml`` section in the ``input.nml`` namelist file and
the model variables that are passed to RTTOV (``&model_nml`` section).
See :ref:`obs_def_rttov_mod`.

.. note::
   The RTTOV support in DART is experimental and may not work under all circumstances.
   Moreover, satellite radiance assimilation is an active area of research. 
   If you encounter issues, please submit them through `Github
   Issues <https://github.com/NCAR/DART/issues>`__.


Compilation and setup
---------------------

At present, DART supports RTTOV v12.3 and v13.  
In version 12.3, RTTOV-direct for visible/infrared/microwave without 
scattering as well as RTTOV-scatt for microwave computations with full 
scattering are supported. The interface to v13 is limited to RTTOV-direct.

If you haven't compiled DART before, it is recommended to compile DART
without RTTOV first, to confirm that everything is working.
To compile DART with RTTOV, you will need to follow these steps:


1. Install RTTOV

Download the RTTOV code and coefficients (for the sensors you need) from this page:
https://www.nwpsaf.eu/site/software/rttov
You will need to register for a free account before downloading the code.
Read the RTTOV user guide carefully as DART primarily acts as
a pass through. Refer to the setup instructions included with the RTTOV
documentation.
Place the coeffient files in the respective directories.
Be aware that there are more coefficient files available once you
download the RTTOV package. There is a
``rtcoef_rttov12/rttov_coef_download.sh`` script that assists in the
process and you can select specific coefficient files or large batches.
There is also a website
https://nwp-saf.eumetsat.int/site/software/rttov/download/coefficients/rttov-v12-coefficient-download/
Build RTTOV as per the instructions.

2. Include RTTOV to the DART build system

Once you have successfully installed RTTOV, you should customize the
``mkmf.template.rttov.gfortran`` file to your own build system, possibly
referring to the other mkmf.template examples for additional information
if you are not using gfortran.

3. Modify input.nml

Go into the model's work directory (``models/wrf/work``) for your model of choice
Add a selection of the observation types listed in
``obs_def_rttov_mod.f90`` to the ``input.nml`` namelist 
(either to ``assimilate_these_obs_types`` or ``evaluate_these_obs_types``).

Add the observation operators to the build by 
adding ``obs_def_rttov_mod.f90`` to the ``&preprocess`` section of the ``input.nml`` file:

.. code-block:: bash

   &preprocess_nml
      input_files              = '../../../observations/forward_operators/obs_def_reanalysis_bufr_mod.f90',
                                 '../../../observations/forward_operators/obs_def_radar_mod.f90',
                                 '../../../observations/forward_operators/obs_def_metar_mod.f90',
                                 '../../../observations/forward_operators/obs_def_dew_point_mod.f90',
                                 '../../../observations/forward_operators/obs_def_rel_humidity_mod.f90',
                                 '../../../observations/forward_operators/obs_def_gts_mod.f90',
                                 '../../../observations/forward_operators/obs_def_rttov_mod.f90',

Use ``obs_def_rttov13_mod.f90`` to compile DART with RTTOV v13.

4. In your model of choice, run ``./quickbuild.sh``.


High-level workflow
-------------------

To assimilate radiances in an OSSE (Observing System Simulation Experiment)
with synthetic observations, you will need to do the following:

   -  Create an observation sequence file using ``./create_obs_sequence``
      and ``./create_fixed_network_seq`` as detailed in the DART
      Getting_Started documentation
   -  Run ``./perfect_model_obs``
   -  Setup your ensemble as appropriate
   -  Run filter and analyze the results in the usual way

To assimilate radiances in an OSE (Observing System Experiment) with real
observations, you will need to do the following:

   -  Run the observation converter for your desired observations
   -  Setup your ensemble as appropriate
   -  Run filter and analyze the results in the usual way


Input data to RTTOV
-------------------

RTTOV simulates radiances by taking in a set of atmospheric and surface
variables to simulate the radiances that would be observed by a
satellite instrument. 

The DART interface basically passes through model variables to RTTOV.
Besides pressure, temperature, and humidity, this can include aerosols, 
trace gases, cloud hydrometeor mixing ratios, and surface variables.

A particular model may not have all of the variables necessary
for RTTOV depending on the model and model setup. 
Although a model may not have the necessary inputs by itself,
in some cases, the defaults in RTTOV based on climatology can be used, 
but at a minimum the following quantities must be defined as state variables:

+-----------------------------+----------------------------------------+
| Quantity                    | Description                            |
+=============================+========================================+
| **QTY_PRESSURE**            | atmospheric pressure in hPa at the     |
|                             | model levels                           |
+-----------------------------+----------------------------------------+
| **QTY_TEMPERATURE**         | atmospheric temperature in K at the    |
|                             | model levels                           |
+-----------------------------+----------------------------------------+
| **QTY_VAPOR_MIXING_RATIO**  | atmospheric humidity mixing ratio in   |
|                             | kg/kg at the model levels              |
+-----------------------------+----------------------------------------+
| **QTY_SURFACE_PRESSURE**    | the surface pressure in hPa            |
+-----------------------------+----------------------------------------+
| **QTY_SURFACE_ELEVATION**   | the surface elevation in km            |
+-----------------------------+----------------------------------------+
| **QTY_2M_TEMPERATURE**      | the atmospheric temperature in K at 2  |
|                             | m above the surface                    |
+-----------------------------+----------------------------------------+
| **QTY_SKIN_TEMPERATURE**    | the surface (skin) temperature in K    |
+-----------------------------+----------------------------------------+
| **QTY_SURFACE_TYPE**        | 0 = land, 1 = water, 2 = sea ice       |
+-----------------------------+----------------------------------------+


If a DART model_mod cannot provide these required quantities, the RTTOV
forward operator will fail and cannot be used. It may be possible to
look up surface elevation or surface type through an look-up table or
“atlas,” although DART does not yet provide such functionality. 2M
temperature in theory could be interpolated based on skin temperature
and the lowest-level model temperature.

Beyond these fields, there are many other optional fields (such as
clouds, trace gases, and aerosols) that can be specified. See
:ref:`obs_def_rttov_mod` for a complete list of values.


Tips for the assimilation of visible/infrared radiances 
-------------------------------------------------------

A good overview over the most important parameters for the radiative transfer
can be found in the RTTOV user guide section "Simulation of UV, visible and IR cloud-affected radiances".

In general, the representation of clouds differs among microphysics parameterizations, which can lead
to biases in comparison with observed radiances.
Moreover, the representation might not be entirely compatible with RTTOV.  
For example, the Thompson microphysics has five cloud hydrometeor categories (cloud water, ice, snow, graupel, and rain), 
while RTTOV only accepts liquid water and ice mixing ratio (plus snow for RTTOV-scatt).

Since cloud optical properties are often not provided by the model, 
RTTOV provides parameterizations for liquid and ice clouds (see the RTTOV user guide for details).
For liquid water clouds there are (abbreviated) "OPAC" and "Deff".

*  The Deff scheme (`clw_scheme=2`) computes optical properties from an effective particle diameter as input.
   By default, DART accesses the model state variable associated to ``QTY_CLOUDWATER_DE`` in the DART namelist.
   Alternatively, users can modify to code to provide a constant value.
*  The OPAC scheme computes optical properties from based on the cloud type 
   (marine/continental, stratus/cumulus, clean/dirty). 
   If the user selects the OPAC scheme (`clw_scheme=1`), DART classifies the cloud type based 
   on maximum vertical velocity in the column and land type.
   In case of cumulus over land, DART currently assigns "Cumulus Continental Clean" , 
   as we lack of aerosol information and cannot differentiate between clean and dirty cumulus.
   This may have some impact on the forward calculations - but in experience the difference 
   in cloud phase (ice versus water) makes a much larger difference. 

Trace gases and aerosols may be important for actual observation system experiments 
using visible/infrared; this may depend on the precise frequencies you wish to use.


Converting real observations to DART format
-------------------------------------------

Note that currently obervation converters are only provided for AIRS,
AMSU/A, GOES, and GMI. These converters can be found in the
observations/obs_converters directories. The L1 converters are the
appropriate converters for the radiance or brightness temperatures
(rather than retrievals). If you need real L1 data for another satellite
(as opposed to running an OSSE with perfect_model_obs where you can
generate your own data), you may be able to use one of these converters
to get you started. We welcome your contributions back to the DART
public repository. Please issue a pull request to
https://github.com/NCAR/DART.

Note that some of the observation converters may require the HDF-EOS
libraries. See the BUILDME script in each directory for help in building
these observation converters.


Current list of known issues
----------------------------

DART support for satellite radiances cannot be considered 100% complete.
The following details the known issues that are being considered with
DART’s support for satellite radiances.

-  DART does not yet provide satellite bias correction capabilities. 
   It may be appropriate to preprocess your radiance
   observations to remove systematic  bias before assimilation, 
   using techniques such as cumulative distribution function (CDF) matching.
-  Cross-channel error correlations are not yet supported. A principal
   component approach has been discussed. For now, the best bet is to
   use a subset of channels that are nearly independent of one another.
-  Vertical localization is an issue for satellite radiances. The main
   choices are to turn off vertical localization, use the maximum peak
   of the weighting function or the cloud-top may be appropriate, or
   explore other options. We consider this an open research problem.
