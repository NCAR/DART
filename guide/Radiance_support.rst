Introduction to DART’s support for RTTOV
========================================

This document serves as an orientation for DART’s support for satellite
radiance assimilation. At the current time, only ECMWF’s RTTOV radiative
transfer model is supported.

DART now includes the ability to use the RTTOV forward operators for
satellite radiance assimilation. This is a new capability for DART,
please submit issues with the `DART
Issues <https://github.com/NCAR/DART/issues>`__ facility.

Note that DART support for RTTOV does not mean that all issues regarding
satellite data assimilation with an ensemble system have been solved.
Rather, the DART team hopes to provide the tools necessary for
researchers to investigate the relevant issues with multiple models and
data assimilation methodologies.

DART supports RTTOV version 12.3. Support for RTTOV v13 has been added to
DART by Lukas Kugler. Use obs_def_rttov13_mod.f90 to compile DART with RTTOV v13.

Both RTTOV-direct for visible/infrared/microwave without scattering as well as RTTOV-scatt for
microwave computations with full scattering are supported for v12.3. DART supports
all features of RTTOV 12.3 as a pass-through from the models to RTTOV.
This includes aerosols, trace gases, clouds, and atmospheric variables.
It also includes directly specifying scattering properties.

However, a particular model may not have all of the variables necessary
for these functions depending on the model and model setup. In some
cases RTTOV default climatologies can be used, but at a minimum the
following quantities must be supplied by the model_mod interpolate:

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


Setting up DART+RTTOV
---------------------

The RTTOV code and coefficients can be downloaded from this page:

https://www.nwpsaf.eu/site/software/rttov

Be aware that there are more coefficient files available once you
download the RTTOV package. There is a
``rtcoef_rttov12/rttov_coef_download.sh`` script that assists in the
process and you can select specific coefficient files or large batches.
There is also a website
https://nwp-saf.eumetsat.int/site/software/rttov/download/coefficients/rttov-v12-coefficient-download/

You will need to register for a free account before downloading the
code.

You should read the RTTOV user guide carefully as DART primarily acts as
a pass through. Refer to the setup instructions included with the RTTOV
documentation.

It may also be useful to refer to:

https://github.com/NCAR/DART/wiki/Getting-Started-with-DART-RTTOV

Once you have successfully installed RTTOV, you should customize the
mkmf.template.rttov.gfortran file to your own build system, possibly
referring to the other mkmf.template examples for additional information
if you are not using gfortran.

There are many namelist options available through input.nml that control
the run-time behavior of the RTTOV model. These are documented in
:ref:`obs_def_rttov_mod`.

To get RTTOV to work with your model, you will need to follow these
steps:

1. Install RTTOV as above
2. Customize your mkmf.template to include the RTTOV libraries and
   include directories
3. Go into the models//work directory for your model of choice
4. Add your observation types (which are listed in
   obs_def_rttov_mod.f90) to the input.nml namelist (assimilate\_ /
   evaluate_these_obs_types)
5. Include observations/forward_operators/obs_def_rttov_mod.f90 in the
   input_files section under &preprocess
6. In your model of choice, run ./quickbuild.sh and ensure the RTTOV
   libraries are built
7. For OSSE runs with perfect_model_obs:

   -  Create an observation sequence file using ./create_obs_sequence
      and ./create_fixed_network_seq as detailed in the DART
      Getting_Started documentation
   -  Run perfect_model_obs
   -  Setup your ensemble as appropriate
   -  Run filter and analyze the results in the usual way

8. For OSE runs:

   -  Run the observation converter for your desired observations
   -  Setup your ensemble as appropriate
   -  Run filter and analyze the results in the usual way

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
   This will be released in the near future.
-  Cross-channel error correlations are not yet supported. A principal
   component approach has been discussed. For now, the best bet is to
   use a subset of channels that are nearly independent of one another.
-  Vertical localization is an issue for satellite radiances. The main
   choices are to turn off vertical localization, use the maximum peak
   of the weighting function or the cloud-top may be appropriate, or
   explore other options. We consider this an open research problem.
