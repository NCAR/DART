.. index:: radiance

Introduction to DART’s support for RTTOV
========================================

DART supports satellite radiance assimilation through an interface to 
the radiative transfer model RTTOV.  RTTOV is a fast radiative transfer model
that is widely used in research and operations. RTTOV contains 
observation operators for UV/visible/infrared/microwave radiances/brightness temperatures. 
RTTOV uses the atmospheric state provided by DART to generate the expected radiance 
observation. The minimum list of atmospheric state variables required is listed in the section below
:ref:`Input data to RTTOV<inputdata>`.  Observations from a wide range of satellites 
(e.g. GOES, FY, METOP, ...) and sensors (e.g. ABI, AMSU-A, SEVIRI, ...) are supported 
(`see the complete list here <https://nwp-saf.eumetsat.int/site/software/rttov/documentation/platforms-supported/>`__).
For more detail on RTTOV see the `RTTOV user guide <https://www.nwpsaf.eu/site/software/rttov/documentation/>`__.

This documentation describes:
 
1. :ref:`Compilation and setup <compilationandsetup>`
2. :ref:`High-level workflow for performing OSSE and OSE <workflow>`
3. :ref:`Input data to RTTOV <inputdata>`
4. :ref:`Tips for the assimilation of visible/infrared radiances <tips>`  
5. :ref:`Converting real observations to DART format <realobs>`
6. :ref:`Current list of known issues <knownissues>`


.. Important::
   The :ref:`obs_def_rttov_mod` describes the 
   ``&obs_def_rttov_nml`` section of the DART namelist.

.. note::
   The RTTOV support in DART has been tested successfully for some applications 
   (e.g. convective storm cells), however, it may not work under all circumstances.
   Moreover, satellite radiance assimilation is an active area of research. 
   If you encounter problems, please submit them through `Github
   Issues <https://github.com/NCAR/DART/issues>`__.

.. _compilationandsetup:

Compilation and setup
---------------------

New versions of RTTOV are released regularly.
At present, DART supports RTTOV v12.3 and v13 for
the simulation of UV/visible/infrared and microwave radiances.

If you haven't compiled DART before, it is recommended to compile DART
first without RTTOV, to confirm that everything is working. Refer to the 
DART compile instructions :doc:`compiling-dart`.
To compile DART with RTTOV, perform the following steps:


**1. Install RTTOV**

Download the RTTOV package after registering for a free account 
`at the EUMETSAT website <https://www.nwpsaf.eu/site/software/rttov>`__. 
Detailed installation instructions can be found in the 'quick start beginners guide' 
under `Software-RTTOV-Documentation <https://nwp-saf.eumetsat.int/site/software/rttov/documentation/>`__.

Apart from the RTTOV source code, you will also need to download the
RTTOV coefficient files, available for different satellite instruments and sensors.
You can select individual files or download all using the ``rttov_coef_download.sh`` script 
which comes with the RTTOV package in the ``rtcoef_rttov?/`` directory.

Lastly, read the RTTOV user guide carefully as DART primarily passes the 
atmospheric model state variables from the atmospheric model to RTTOV.

**2. Include RTTOV to the DART build system**

Once you have successfully installed RTTOV, you should customize the
DART ``mkmf.template.rttov.gfortran`` file to your own build system.  
Refer to the other DART ``mkmf.template`` examples that best matches your compile
setup if you are not using gfortran.

**3. Modify input.nml**

Go into the work directory (e.g. ``${DART_install}/models/wrf/work``) for your 
atmospheric model of choice. These instructions use the WRF model as an example.
Add the appropriate  observation types for your application listed in
``obs_def_rttov_mod.f90`` to the ``input.nml`` namelist 
(either to ``assimilate_these_obs_types`` or ``evaluate_these_obs_types``).

Include the RTTOV observation operators within the build process by editing
the ``input_files`` and ``quantity_files`` namelist variables in the
``&preprocess`` section of the ``input.nml`` file:


.. code-block:: bash

   &preprocess_nml
      quantity_files          =  '../../../assimilation_code/modules/observations/atmosphere_quantities_mod.f90',
                                 '../../../assimilation_code/modules/observations/ocean_quantities_mod.f90',
                                 '../../../assimilation_code/modules/observations/chemistry_quantities_mod.f90',
                                 '../../../assimilation_code/modules/observations/land_quantities_mod.f90'
      input_files              = '../../../observations/forward_operators/obs_def_reanalysis_bufr_mod.f90',
                                 '../../../observations/forward_operators/obs_def_radar_mod.f90',
                                 '../../../observations/forward_operators/obs_def_metar_mod.f90',
                                 '../../../observations/forward_operators/obs_def_dew_point_mod.f90',
                                 '../../../observations/forward_operators/obs_def_rel_humidity_mod.f90',
                                 '../../../observations/forward_operators/obs_def_gts_mod.f90',
                                 '../../../observations/forward_operators/obs_def_rttov_mod.f90',


**4. For your model of choice, run ./quickbuild.sh.**

.. _workflow:

High-level workflow
-------------------

Prior to performing a perfect model experiment (PMO) or real data assimilation using
the DART executables ``./perfect_model_obs`` or ``./filter``, you need to link 
the RTTOV coefficient files to the expected coefficient filename in 
the work directory of your model.
An example of coefficent files is provided below. These should be customized
for your application:

.. code-block:: bash

   export DART=${DART_install}
   export RTTOV=${RTTOV_install}
   cd $DART/models/wrf/work/
   ln -sf $DART/observations/forward_operators/rttov_sensor_db.csv .
   ln -sf $RTTOV/rtcoef_rttov13/cldaer_visir/sccldcoef_msg_4_seviri.dat  .
   ln -sf $RTTOV/rtcoef_rttov13/mfasis_lut/rttov_mfasis_cld_msg_4_seviri_deff.H5 rttov_mfasis_cld_msg_4_seviri.H5


Edit ``input.nml`` to define the:

   *  RTTOV namelist variables within ``&obs_def_rttov_nml``
   *  Model state variables necessary to support the forward operators within ``&model_nml``
   *  All DART specifc namelist variables including ``&filter_nml`` or ``&perfect_model_obs_nml`` 


**You can test the DART-RTTOV interface by performing an OSSE (e.g. perfect model experiment (PMO)).** 
This uses the forward operators, exactly as is done in a real assimilation using
``filter``, but at a much lower cost (i.e. only 1 model member required).  To perform a PMO 
complete the following steps:

   -  Generate a nature run state file from which to harvest the synthetic observation(s).
      If you are using WRF, you can copy any wrfout file to the work directory 
      and rename it to ``wrfinput_d01``. If you use WRF in ideal mode, make sure 
      that the file contains valid geographical coordinates.
   -  Create an observation sequence file using ``./create_obs_sequence``
      and ``./create_fixed_network_seq`` as detailed in the DART
      :doc:`creating-obs-seq-synthetic` documentation to generate an ``obs_seq.in``
      for the same valid time as the nature run state file.
   -  Run ``./perfect_model_obs`` to generate the synthetic observation (``obs_seq.out``)
   -  Check the ``obs_seq.out`` file to identify if there is a QC value not equal 0 (see below).
    
.. Important::

   If the DART QC code in the obs_seq.out is not 0, something went wrong. If the forward 
   operator fails in pmo, it gets assigned a DART QC value of 1000 in addition to any error code 
   from the model. For instance, an error code of 1003 means that the WRF 
   ``model_mod:model_interpolate()`` routine returned an error code of 3 ... 
   "3 = unsupported obs kind". You can determine what variables were needed by the 
   ``$DART/observations/forward_operators/obs_def_rttov_mod.f90`` 
   get_expected_radiance() routine and check to see that they are specified to be part of the DART 
   WRF state and that the WRF model_interpolate() routine supports them.
..

  -  Run ``./filter`` to assimilate the synthetic observation
  -  Check the ``obs_seq.final`` file to confirm a successful assimlation (QC value = 0)

**Perform a real observation assimilation (e.g. OSE)**

   - Run an observation converter following the :doc:`creating-obs-seq-real` documentation.
     to generate an ``obs_seq.out`` file.  At present, there are three observation converters
     for radiance: AIRS, GMI, and AMSU/A.
     Be advised that the units of the forward operator must match the units of the observations 
     in the observation sequence files. Presently, the DART/RTTOV implementation is such that 
     all observations of QTY_BRIGHTNESS_TEMPERATURE are in degrees Kelvin, all observations of 
     QTY_RADIANCE are as described in the RTTOV v12 user guide V1.3 (p54): "mW/cm-1/sr/sq.m" 
   -  Run ``./filter`` to assimilate the real observation
   -  Check the ``obs_seq.final`` file to confirm a successful assimlation (QC value = 0) 


.. _inputdata:

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
for RTTOV depending on the user application. 
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
independent of DART. The 2M temperature in theory could be interpolated based on 
skin temperature and the lowest-level model temperature.

Beyond these fields, there are many other optional fields (such as
clouds, trace gases, and aerosols) that can be specified. See
:ref:`obs_def_rttov_mod` for a complete list of values.

.. _tips:

Tips for the assimilation of visible/infrared radiances 
-------------------------------------------------------

We recommended to study the user guide for `RTTOV-v12 <https://nwp-saf.eumetsat.int/site/software/rttov/rttov-v12/>`__ 
or `RTTOV-v13 <https://nwp-saf.eumetsat.int/site/software/rttov/rttov-v13/>`__,
especially section 8.5. "Simulation of UV, visible and IR cloud-affected radiances".

In general, the representation of clouds differ between microphysics parameterizations, which can lead
to biases in comparison with observed radiances.
Moreover, the representation might not be entirely compatible with RTTOV.  
For example, the Thompson microphysics scheme has five cloud hydrometeor categories (cloud water, ice, snow, graupel, and rain), 
while RTTOV only accepts liquid water and ice mixing ratio (plus snow for RTTOV-scatt).


**Specifying liquid and ice cloud optical properties:**

#. Liquid water clouds

   *  The Deff scheme (`clw_scheme=2`) computes optical properties from an effective particle diameter as input.
      If the model state variable associated with ``QTY_CLOUDWATER_DE`` exists in the DART state (&model_nml) then
      DART will pass these modeled values onto RTTOV. Alternatively, if that model state does not exist, then
      the code defaults to a contant value (e.g. clwde = 20.0 microns) as defined in the ``obs_def_rttov_mod.f90``.
      Be aware that the code does not account for unit conversion.  WRF effective particle diameter, for example, 
      requires a meter to micron conversion.  If using WRF you should provide the ``&model_nml`` with the 
      effective liquid water radius variable (``RE_QC``).  This is not a default WRF output, but can be added through
      the WRF diagnostic namelist setting (``&diags solar_diagnostics = 1``). Be aware that the units conversion
      (radius (meters) to diameter (microns)) is applied within the ``model_mod.f90`` and not the RTTOV related code.
   *  The OPAC scheme computes optical properties based on the cloud type 
      (marine/continental, stratus/cumulus, clean/dirty). 
      If the user selects the OPAC scheme (`clw_scheme=1`), DART classifies the cloud type based 
      on the maximum vertical velocity (``QTY_VERTICAL_VELOCITY``) in the column and land type. 
      In case of cumulus over land, DART currently assigns "Cumulus Continental Clean" , 
      as we lack aerosol information and cannot differentiate between clean and dirty cumulus.
      This may have some impact on the forward calculations but in practice the difference 
      in cloud phase (ice versus water) makes a much larger difference. 

#. Ice clouds

   *  There is a large uncertainty in the representation of ice-phase clouds in forecast models and 
      radiative transfer models (see Li et al. 2022), due to different assumptions in particle size distributions
      and particle shape.
   *  For a realistic simulation of infrared radiances, include at least the hydrometeor categories
      `ice` and `snow` in the DART state vector. Only the variables in the DART state vector 
      will be passed to RTTOV to compute the expected radiance. If you want to prescribe the ice particle
      diameter directly from the model to RTTOV, assign the appropriate model variable to 
      ``QTY_CLOUD_ICE_DE`` within ``&model_nml`` and set ``ice_scheme = 1`` and ``use_icede = .true.``.  
      If ``QTY_CLOUD_ICE_DE`` is not defined then DART defaults to a fixed particle 
      size diameter (e.g. icede = 60.0 microns). If using WRF you should provide the ``&model_nml`` with the 
      effective ice radius variable (``RE_QI``).  This is not a default WRF output, but can be added through
      the WRF diagnostic namelist setting (``&diags solar_diagnostics = 1``). Be aware that the units conversion
      (radius (meters) to diameter (microns)) is applied within the ``model_mod.f90`` and not the RTTOV related code.
   *  Regarding visible reflectance, DART follows Kostka et al. (2014), section 3.a, 
      and counts only 10% of 'snow' towards the cloud ice concentration.
      This is because large particles tend to have a smaller scattering cross-section than 
      many small particles of the same total mass.
      The percentage value can be used as a tuning parameter in real applications.



**Specifying addsolar namelist option:** 
See :ref:`obs_def_rttov_mod` for namelist options.

The ``addsolar`` option allows the user to specify the azimuth and zenith angle of the sun such that the
expected radiance values account for scattering of solar radiation.  It should be noted that specifying the
azimth and zenith angle are not mandatory metadata to account for solar. Alternatively,  RTTOV can also 
calculate the impact of solar based on the latitude, longitude, date and time associated with the observation.

**Specifying cfrac_data namelist option:**  
See :ref:`obs_def_rttov_mod` for namelist options.

The default setting in DART is **not** to use cloud fraction data (``cfrac_data = false``) 
to account for the impact of clouds on radiation.  This may seem counter-intuitive given 
that RTTOV uses a weighted linear combination of cloudy 
and clear sky fraction to calculate radiance, where the cloudy fraction is specified by the 
hydrometeor data (e.g. clw_data, rain_data, ciw_data, snow_data, graupel_data, hail_data). However, when 
``cfrac_data`` is not specified DART will automatically prescribe a cloud fraction of 1 for all locations.  
Therefore, for high resolution simulations (e.g. several kms) the clouds are much larger than the grid resolution.  
In general, the recommendation is to not include the ``cfrac_data`` for high resolution and/or convection 
permitting simulations.  On the other hand, for coarse and/or parameterized convection simulations setting 
``cfrac_data = true`` is recommended.     


.. _realobs:

Converting real observations to DART format
-------------------------------------------

DART provides observation converters for AIRS,
AMSU/A, GOES, and GMI satellite sensors. These converters can be found in the
${DART_install}/observations/obs_converters directories. The L1 converters are the
appropriate converters for the radiance or brightness temperatures
(rather than L2 retrievals, i.e. derived physical properties). If you need real L1 data for another satellite
(as opposed to running an OSSE with perfect_model_obs where you can
generate your own data), you may be able to use one of these converters
as a template to get you started. We welcome your contributions back to the DART
public repository. Please issue a pull request to
https://github.com/NCAR/DART.

Note that some of the observation converters may require the HDF-EOS
libraries, available from https://hdfeos.org/.

.. _knownissues:

Current list of known issues
----------------------------

DART support for satellite radiances may not include all the features required
for your application. For example, the end user should consider how to best
address the following challenges in satellite DA.

-  DART does not automatically provide satellite radiance bias correction. 
   It may be appropriate to preprocess your radiance
   observations to remove systematic bias before assimilation, 
   using techniques such as cumulative distribution function (CDF) matching.
-  Cross-channel error correlations are not accounted for in DART. 
   To account for correlations, It is recommended to use a subset of channels that are 
   nearly independent of one another. Be sure to use channels most sensitive to the atmospheric
   property(s) of interest.
-  Applying vertical localization is an ongoing research challenge for satellite radiances, 
   given radiance is a spatially integrated measure of atmospheric properties without
   a single location.  One approach is to turn off vertical localization altogether.  
   Another approach is to assign a vertical location based on the maximum peak of 
   the weighting function (i.e. vertical location of highest sensitivity to property of interest) 
   or the cloud-top as appropriate.


References
----------
  *  Li et al. (2022) “Satellite All-Sky Infrared Radiance Assimilation: Recent Progress and Future Perspectives.” Advances in Atmospheric Sciences 39(1): 9–21. doi:10.1007/s00376-021-1088-9.
  *  Kostka et al. (2014) “Observation Operator for Visible and Near-Infrared Satellite Reﬂectances.” Journal of Atmospheric and Oceanic Technology 31(6): 1216–33. doi:10.1175/JTECH-D-13-00116.1.
