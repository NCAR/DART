.. _available_observation_converters:

Available observation converter programs
========================================

The ``DART/observations/obs_converters`` directory contains a variety of
converter programs to read various external formats and convert the observations
into the format required by DART.

Each directory has at least one converter:


-  ``AIRS``: :ref:`airs`
-  ``ARVOR``: :ref:`arvor`
-  ``AURA``: See ``DART/observations/obs_converters/AURA``
-  ``Aviso+/CMEMS``: :ref:`aviso`
-  ``Ameriflux``: 
    
    - :ref:`fluxnetfull_to_obs`
    - :ref:`level4_to_obs`
-  ``BATS``: :ref:`bats`
-  ``CHAMP``: :ref:`champ`
-  ``cice``: :ref:`cice_to_obs`
-  ``CNOFS``: See ``DART/observations/obs_converters/CNOFS``
-  ``CONAGUA``: :ref:`conagua`
-  ``COSMOS``: :ref:`cosmos`
-  ``CrocoLake``: :ref:`crocolake`
-  ``DWL``: :ref:`dwl`
-  ``GMI``: :ref:`gmi`
-  ``GOES``: :ref:`goes`
-  ``GPSPW``: :ref:`gpspw`
-  ``GRACE``: See ``DART/observations/obs_converters/GRACE``
-  ``GSI2DART``: :ref:`gsi2dart`
-  ``GTSPP``: :ref:`gtspp`
-  ``HFradar``: :ref:`hfradar`
-  ``IODA``: :ref:`ioda2obsq <ioda2obsq>`
-  ``MADIS``: :ref:`madis`
-  ``MIDAS``: :ref:`midas`
-  ``MODIS``: :ref:`modis15` 
-  ``MODIS``: :ref:`modis29`
-  ``MPD``: See ``DART/observations/obs_converters/MPD``
-  ``NASA_Earthdata``: :ref:`nasa_earthdata`
-  ``NCEP``:

    - (prepbufr -> ascii) :ref:`ncep_prepbufr`
    - (ascii -> obs_seq) :ref:`ncep_ascii`
-  ``NSIDC``: :ref:`nsidc_smap_l2`
-  ``ocean color``: :ref:`ocean_color`
-  ``ROMS``: :ref:`roms`
-  ``SIF``: :ref:`sif`
-  ``SSEC``: :ref:`ssec`
-  ``SST``: :ref:`sst`
-  ``SSUSI``: :ref:`ssusi`
-  ``WOD``: :ref:`wod`
-  ``gnd_gps_vtec``: :ref:`gnd_gps_vtec`
-  ``GPS``: :ref:`gps`
-  ``ok_mesonet``: :ref:`ok_mesonet`
-  ``QuikSCAT``: :ref:`quikscat`
-  ``Radar``: :ref:`radar`
-  ``snow``: :ref:`snow`
-  ``SVP``: :ref:`svp`
-  ``Text``: :ref:`text`
-  ``text_GITM``: See ``DART/observations/obs_converters/text_GITM``
-  ``tpw``: :ref:`tpw`
-  ``Tropical Cyclones``: :ref:`tropical_cyclone`
-  ``Var (3D/4D)``: :ref:`var`

   - ``Var (little-r)``: :ref:`littler_tf_dart`
   - ``Var (radar)``: :ref:`rad_3dvar_to_dart`


In addition the following external program produces DART observation sequence
files:

-  `Observation Processing And Wind Synthesis
   (OPAWS) <http://code.google.com/p/opaws/>`__: OPAWS can process NSF NCAR Dorade
   (sweep) and NSF NCAR EOL Foray (netCDF) radar data. It analyzes (grids) data in
   either two-dimensions (on the conical surface of each sweep) or
   three-dimensions (Cartesian). Analyses are output in netCDF, Vis5d, and/or
   DART (Data Assimilation Research Testbed) formats.

Contact the `DART development group <mailto:dart@ucar.edu>`__ if you
have observations in a different format that you want to convert. We can
give you advice and pointers on how to approach writing the code.

.. _synthetic_observations:

Synthetic observations
--------------------------

For generating synthetic observations, see the documentation for the 
:ref:`create_obs_sequence`.
You can also generate observation files based on text input. See the
documentation for the :ref:`text`.
Or for simulating a large complex observing system, you can use the DART
library routines in a Fortran program to compute the observation information
and have the DART routines write the output file.

There are also a couple utilities of note:

-  :ref:`even_sphere` - a utility for generating evenly-spaced
   observation locations that can then be used in a perfect model experiment.
-  :ref:`obs_error` - modules that specify observation errors
   based on what is used by ECMWF and NCEP

To learn how to run a model with a set of observations that have only
locations, types, and times, and have the forward operators compute the
observation values, see the documentation for the
:ref:`pmo`.
