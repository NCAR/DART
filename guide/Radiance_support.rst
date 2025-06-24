Introduction to DART’s support for RTTOV
========================================

DART supports satellite radiance assimilation through an interface to 
the radiative transfer model RTTOV. 
RTTOV is a fast radiative transfer model that is widely used in research 
and operations. RTTOV contains observation operators for visible/infrared and microwave 
radiances/brightness temperatures. RTTOV uses the atmospheric state to generate the radiance 
observation. The minimum list of atmospheric state variables required is listed in the section below
:ref:`Input data to RTTOV<inpudata>`.  Observations from a wide range of satellites 
(e.g. GOES, FY, METOP, ...) and sensors (e.g. ABI, AMSU-A, SEVIRI, ...) are supported 
(`see the complete list here <https://nwp-saf.eumetsat.int/site/software/rttov/documentation/platforms-supported/>`__).
For more detail on RTTOV see the `RTTOV user guide <https://www.nwpsaf.eu/site/software/rttov/documentation/>`__.

This documentation describes:
 
1. `Compilation and setup<compilationandsetup>`
2. `High-level workflow<workflow>`
3. `Input data to RTTOV<inputdata>`
4. `Tips for the assimilation of visible/infrared radiances<tips>`  
5. `Converting real observations to DART format<realobs>`
6. `Current list of known issues<knownissues>`


.. seealso::
   The page :ref:`obs_def_rttov_mod` describes the 
   ``&obs_def_rttov_nml`` section of the DART namelist.

.. note::
   The RTTOV support in DART has been tested successfully for some applications 
   (e.g. convective storm cells), however,  may not work under all circumstances.
   Moreover, satellite radiance assimilation is an active area of research. 
   If you encounter issues, please submit them through `Github
   Issues <https://github.com/NCAR/DART/issues>`__.

.._compilationandsetup:

Compilation and setup
---------------------

At present, DART supports RTTOV v12.3 and v13.  
In version 12.3, RTTOV-direct for visible/infrared/microwave without 
scattering as well as RTTOV-scatt for microwave computations with full 
scattering are supported. The interface to v13 is limited to RTTOV-direct.

If you haven't compiled DART before, it is recommended to compile DART
without RTTOV first, to confirm that everything is working. Refer to the 
DART compile instructions :doc:`compiling-dart`.
To compile DART with RTTOV, you will need to follow these steps:


1. Install RTTOV

Download the RTTOV code and coefficients (for the sensors you need) from this page:
https://www.nwpsaf.eu/site/software/rttov
You will need to register for a free account before downloading the code.
Read the RTTOV user guide carefully as DART primarily passes the atmospheric model state variables
to RTTOV. Refer to the setup instructions included with the RTTOV documentation.
Place the coeffient files in the respective directories.
Be aware that there are more coefficient files available once you
download the RTTOV package. There is a
``rtcoef_rttov12/rttov_coef_download.sh`` script that assists in the
process and you can select specific coefficient files or large batches.
There is also a website
for `RTTOVv12 <https://nwp-saf.eumetsat.int/site/software/rttov/download/coefficients/rttov-v12-coefficient-download/>`__ and `RTTOVv13 <https://nwp-saf.eumetsat.int/site/software/rttov/download/coefficients/rttov-v13-coefficient-download/>`__ coefficent files.
Build RTTOV as per the instructions.

2. Include RTTOV to the DART build system

Once you have successfully installed RTTOV, you should customize the
``mkmf.template.rttov.gfortran`` file to your own build system.  
Refer to the other DART mkmf.template examples that best matches your compile
setup if you are not using gfortran.

3. Modify input.nml

Go into the work directory (e.g. ``{DART_path}/models/wrf/work``) for your atmospheric model of choice
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

.._workflow:

High-level workflow
-------------------

To assimilate radiances in an OSSE (Observing System Simulation Experiment)
with synthetic observations, you will need to do the following:

   -  Create an observation sequence file using ``./create_obs_sequence``
      and ``./create_fixed_network_seq`` as detailed in the DART
      :doc: `documentation <creating_obs_seq_synthetic>` to generate an ``obs_seq.in``
   -  Run ``./perfect_model_obs`` to generate synthetic obs within the ``obs_seq.out``
   -  Setup your ensemble as appropriate
   -  Run ``./filter`` and analyze the results in the usual way

To assimilate radiances in an OSE (Observing System Experiment) with real
observations, you will need to do the following:

   -  Run the :doc: `observation converter <creating-obs_seq-real>` for your desired radiance observation.
   -  Setup your ensemble as appropriate
   -  Run ``./filter`` and analyze the results in the usual way

.._inputdata:

Input data to RTTOV
-------------------

RTTOV simulates radiances by taking in a set of atmospheric and surface
variables to simulate the radiances that would be observed by a
satellite instrument. 

The DART interface basically passes through model variables to RTTOV.
Besides mandatory inputs such as pressure, temperature, and humidity, the
user can specify information on aerosols, trace gases, and cloud hydrometeor mixing ratios 
depending on the application of interest.

A particular atmospheric model may not have all of the variables necessary
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
specify surface elevation or surface type directly to RTTOV through a look-up table,
indpendent of DART. The 2M temperature in theory could be interpolated based on 
skin temperature and the lowest-level model temperature.

Beyond these fields, there are many other optional fields (such as
clouds, trace gases, and aerosols) that can be specified. See
:ref:`obs_def_rttov_mod` for a complete list of values.

.._tips:

Tips for the assimilation of visible/infrared radiances 
-------------------------------------------------------

A good overview over the most important parameters for the radiative transfer
can be found in the RTTOV user guide section "Simulation of UV, visible and IR cloud-affected radiances".

In general, the representation of clouds differs among microphysics parameterizations, which can lead
to biases in comparison with observed radiances.
Moreover, the representation might not be entirely compatible with RTTOV.  
For example, the Thompson microphysics has five cloud hydrometeor categories (cloud water, ice, snow, graupel, and rain), 
while RTTOV only accepts liquid water and ice mixing ratio (plus snow for RTTOV-scatt).


**Specifying liquid and ice cloud optical properties:**

#. Liquid water clouds

   *  The Deff scheme (`clw_scheme=2`) computes optical properties from an effective particle diameter as input.
      By default, DART accesses the model state variable associated with ``QTY_CLOUDWATER_DE`` in the DART namelist.
      Alternatively, users can modify the code to specify a constant value.
   *  The OPAC scheme computes optical properties from based on the cloud type 
      (marine/continental, stratus/cumulus, clean/dirty). 
      If the user selects the OPAC scheme (`clw_scheme=1`), DART classifies the cloud type based 
      on the maximum vertical velocity (``QTY_VERTICAL_VELOCITY``) in the column and land type. 
      In case of cumulus over land, DART currently assigns "Cumulus Continental Clean" , 
      as we lack of aerosol information and cannot differentiate between clean and dirty cumulus.
      This may have some impact on the forward calculations - but in practice the difference 
      in cloud phase (ice versus water) makes a much larger difference. 

#. Ice clouds

   *  See the RTTOV user guide.


**Specifying `addsolar` namelist option:**

The `addsolar` option allows the user to specify the azimuth and zenith angle of the sun such that the
expected radiance values account for scattering of solar radiation.  It should be noted that specifying the
azimth and zenith angle are not mandatory metadata to account for solar. Alternatively,  RTTOV can also 
calculate the impact of solar based on the latitude, longitude, date and time associated with the observation.

**Specifying `cfrac_data` namelist option:**

The default setting in DART is **not** to use `cfrac_data` (.false.) to account for the impact of clouds
on radiation.  This may seem counter-intuitive given that RTTOV uses a weighted linear combination of cloudy 
and clear sky fraction to calculate radiance, where the cloudy fraction is specified by the 
hydrometeor data (e.g. clw_data, rain_data, ciw_data, snow_data, graupel_data, hail_data). However, when 
`cfrac_data` is not specified DART will automatically prescribe a cloud fraction of 1 for all locations.  
Therefore, for high resolution simulations (e.g. several kms) the clouds are much larger than the grid resolution.  
In general, the recommendation is to not include the `cfrac_data` for high resolution and/or convection 
permitting simulations.  On the other hand, for coarse and/or parameterized convection simulations specifying 
`cfrac_data` is recommended.     


.._realobs:

Converting real observations to DART format
-------------------------------------------

Note that currently observation converters are only provided for AIRS,
AMSU/A, GOES, and GMI. These converters can be found in the
{DART_path}/observations/obs_converters directories. The L1 converters are the
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

.._knownissues:

Current list of known issues
----------------------------

DART support for satellite radiances may not include all the features required
for your application. For example, the end user should consider how to best
address the following challenges in satellite DA.

-  DART does not automatically provide satellite bias correction capabilities. 
   It may be appropriate to preprocess your radiance
   observations to remove systematic  bias before assimilation, 
   using techniques such as cumulative distribution function (CDF) matching.
-  Cross-channel error correlations are not accounted for in DART. 
   It is recommended to use a subset of channels that are nearly independent 
   of one another.
-  Vertical localization is an ongoing research challenge for satellite radiances, 
   given it is an integrated measure of atmospheric properties.  
   One option is to turn off vertical localization altogether.  
   Another option is to assign a vertical location based on the maximum peak of 
   the weighting function or the cloud-top as appropriate.
