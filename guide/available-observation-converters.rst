Adding your observations to DART
================================

First, you should understand that DART already supports a tremendous variety of
observations. To fully support an observation means that the observation can be
converted from its native format to the DART observation sequence format and
that the observation forward operator is already implemented. Keep in mind that
forward operators are not specific to any one model.

The observation converters are in the *observations/obs_converter* directory and
you should look there for the documentation describing which converters are
available.

The forward operators are functionally or logically grouped into Fortran modules
in the *observations/forward_operator* directory. DART employs a ‘contractual’
style of programming in that the forward operator requests information from the
model, and if the model cannot provide it, the forward operator may request
different information in an attempt to collect the information needed to apply
the operator. If the model cannot provide any of the required information, the
forward operator fails, the DART QC for that observation is set to the
appropriate value, and the program continues.

Available observation converter programs
========================================

The ``DART/observations/obs_converters`` directory contains a variety of
converter programs to read various external formats and convert the observations
into the format required by DART.

Each directory has at least one converter:

-  ``AIRS``: :doc:`/observations/obs_converters/AIRS/README`
-  ``AURA``: See ``DART/observations/obs_converters/AURA``
-  ``Aviso+/CMEMS``: :doc:`../observations/obs_converters/AVISO/AVISO`
-  ``Ameriflux``: :doc:`../observations/obs_converters/Ameriflux/level4_to_obs`
-  ``CHAMP``: :doc:`../observations/obs_converters/CHAMP/work/README`
-  ``cice``: :doc:`../observations/obs_converters/cice/cice_to_obs`
-  ``CNOFS``: See ``DART/observations/obs_converters/CNOFS``
-  ``CONAGUA``: :doc:`../observations/obs_converters/CONAGUA/README`
-  ``COSMOS``: :doc:`../observations/obs_converters/COSMOS/COSMOS_to_obs`
-  ``DWL``: :doc:`../observations/obs_converters/DWL/dwl_to_obs`
-  ``GMI``: :doc:`../observations/obs_converters/GMI/README`
-  ``GOES``: :doc:`../observations/obs_converters/GOES/README`
-  ``GPSPW``: :doc:`../observations/obs_converters/GPSPW/README`
-  ``GRACE``: See ``DART/observations/obs_converters/GRACE``
-  ``GSI2DART``: :doc:`../observations/obs_converters/GSI2DART/readme`
-  ``GTSPP``: :doc:`../observations/obs_converters/GTSPP/GTSPP`
-  ``MADIS``: :doc:`../observations/obs_converters/MADIS/MADIS`
-  ``MIDAS``: :doc:`../observations/obs_converters/MIDAS/MIDAS_to_obs`
-  ``MODIS``: :doc:`../observations/obs_converters/MODIS/MOD15A2_to_obs`
-  ``MPD``: See ``DART/observations/obs_converters/MPD``
-  ``NCEP``: (prepbufr -> ascii) :doc:`../observations/obs_converters/NCEP/prep_bufr/prep_bufr`
-  ``NCEP``: (ascii -> obs_seq) :doc:`../observations/obs_converters/NCEP/ascii_to_obs/create_real_obs`
-  ``ROMS``: :doc:`../observations/obs_converters/ROMS/ROMS`
-  ``SIF``: :doc:`../observations/obs_converters/SIF/SIF_to_obs_netcdf`
-  ``SSEC``: :doc:`../observations/obs_converters/SSEC/SSEC`
-  ``SST``: :doc:`../observations/obs_converters/SST/SST`
-  ``SSUSI``: :doc:`../observations/obs_converters/SSUSI/convert_f16_edr_dsk`
-  ``WOD``: :doc:`../observations/obs_converters/WOD/WOD`
-  ``gnd_gps_vtec``: :doc:`../observations/obs_converters/gnd_gps_vtec/README`
-  ``GPS``: :doc:`../observations/obs_converters/gps/gps`
-  ``ok_mesonet``: :doc:`../observations/obs_converters/ok_mesonet/ok_mesonet`
-  ``QuikSCAT``: :doc:`../observations/obs_converters/quikscat/QuikSCAT`
-  ``Radar``: :doc:`../observations/obs_converters/radar/README`
-  ``snow``: :doc:`../observations/obs_converters/snow/snow_to_obs`
-  ``Text``: :doc:`../observations/obs_converters/text/text_to_obs`
-  ``text_GITM``: See ``DART/observations/obs_converters/text_GITM``
-  ``tpw``: :doc:`../observations/obs_converters/tpw/tpw`
-  ``Tropical Cyclones``: :doc:`../observations/obs_converters/tropical_cyclone/tc_to_obs`
-  ``Var (little-r)``: :doc:`../observations/obs_converters/var/littler_tf_dart`
-  ``Var (radar)``: :doc:`../observations/obs_converters/var/rad_3dvar_to_dart`

There are also a couple utilities of note:

-  :doc:`../observations/obs_converters/even_sphere/README` - a utility for generating evenly-spaced
   observation locations that can then be used in a perfect model experiment.
-  :doc:`../observations/obs_converters/obs_error/README` - modules that specify observation errors
   based on what is used by ECMWF and NCEP

In addition the following external program produces DART observation sequence
files:

-  `Observation Processing And Wind Synthesis
   (OPAWS) <http://code.google.com/p/opaws/>`__: OPAWS can process NSF NCAR Dorade
   (sweep) and NSF NCAR EOL Foray (netCDF) radar data. It analyzes (grids) data in
   either two-dimensions (on the conical surface of each sweep) or
   three-dimensions (Cartesian). Analyses are output in netCDF, Vis5d, and/or
   DART (Data Assimilation Research Testbed) formats.

For generating synthetic observations, see the documentation for the 
:doc:`../assimilation_code/programs/create_obs_sequence/create_obs_sequence`.
You can also generate observation files based on text input. See the
documentation for the :doc:`../observations/obs_converters/text/text_to_obs`.
Or for simulating a large complex observing system, you can use the DART
library routines in a Fortran program to compute the observation information
and have the DART routines write the output file.

To learn how to run a model with a set of observations that have only
locations, types, and times, and have the forward operators compute the
observation values, see the documentation for the
:doc:`/assimilation_code/programs/perfect_model_obs/perfect_model_obs`.
