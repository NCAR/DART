Introduction to DART’s support for RTTOV
========================================

DART supports satellite radiance assimilation through an interface to 
the radiative transfer model RTTOV. 
RTTOV is a fast radiative transfer model that is widely used in research 
and operations to simulate microwave, infrared, and visible radiances.
It simulates radiances by taking in a set of atmospheric and surface
variables to simulate the radiances that would be observed by a
satellite instrument. 


The DART interface to RTTOV passes through model variables to RTTOV.
This includes aerosols, trace gases, clouds, and atmospheric variables.
A particular model may not have all of the variables necessary
for RTTOV depending on the model and model setup. In some
cases RTTOV default climatologies can be used, but at a minimum the
following quantities must be defined as state variables:

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

.. note::
   The RTTOV support DART support cannot be considered 100% complete and 
   may not work under all circumstances.
   Moreover, satellite radiance assimilation is an active area of research. 
   If you encounter issues, please submit them through `Github
   Issues <https://github.com/NCAR/DART/issues>`__.

Setting up DART+RTTOV
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
adding obs_def_rttov_mod.f90 to the ``&preprocess`` section of the ``input.nml`` file:

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


Assimilating radiances
----------------------

To assimilate radiances in an OSSE (Observing System Simulation Experiment)
with synthetic observations, you will need to do the following:

   -  Create an observation sequence file using ./create_obs_sequence
      and ./create_fixed_network_seq as detailed in the DART
      Getting_Started documentation
   -  Run ./perfect_model_obs
   -  Setup your ensemble as appropriate
   -  Run filter and analyze the results in the usual way

To assimilate radiances in an OSE (Observing System Experiment) with real
observations, you will need to do the following:

   -  Run the observation converter for your desired observations
   -  Setup your ensemble as appropriate
   -  Run filter and analyze the results in the usual way


Advice regarding RTTOV
----------------------

The run-time behavior of RTTOV is mostly controlled by the 
``&obs_def_rttov_nml`` section in the ``input.nml`` namelist file and
the model variables that are passed to RTTOV (``&model_nml`` section).
See :ref:`obs_def_rttov_mod`.

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
-  Cross-channel error correlations are not yet supported. A principal
   component approach has been discussed. For now, the best bet is to
   use a subset of channels that are nearly independent of one another.
-  Vertical localization is an issue for satellite radiances. The main
   choices are to turn off vertical localization, use the maximum peak
   of the weighting function or the cloud-top may be appropriate, or
   explore other options. We consider this an open research problem.
