Available observation converter programs
========================================

The **DART/observations/obs_converters** directory contains a variety of
converter programs to read various external formats and convert the observations
into the format required by DART.

Each directory has at least one converter:

-  `AIRS <../observations/obs_converters/AIRS/AIRS.html>`__
-  `Aviso+/CMEMS <../observations/obs_converters/AVISO/AVISO.html>`__
-  `Ameriflux <../observations/obs_converters/Ameriflux/level4_to_obs.html>`__
-  `cice <../observations/obs_converters/cice/cice_to_obs.html>`__
-  `COSMOS <../observations/obs_converters/COSMOS/COSMOS_to_obs.html>`__
-  `DWL <../observations/obs_converters/DWL/dwl_to_obs.html>`__
-  `GOES <../observations/obs_converters/GOES/README.md>`__
-  `GPSPW <../observations/obs_converters/GPSPW/README>`__
-  `GSI2DART <../observations/obs_converters/GSI2DART/README>`__
-  `GTSPP <../observations/obs_converters/GTSPP/GTSPP.html>`__
-  `MADIS <../observations/obs_converters/MADIS/MADIS.html>`__
-  `MIDAS <../observations/obs_converters/MIDAS/MIDAS_to_obs.html>`__
-  `MODIS <../observations/obs_converters/MODIS/MOD15A2_to_obs.htm>`__
-  `NCEP (prepbufr -> ascii) <../observations/obs_converters/NCEP/prep_bufr/prep_bufr.html>`__
-  `NCEP (ascii -> obs_seq) <../observations/obs_converters/NCEP/ascii_to_obs/create_real_obs.html>`__
-  `ROMS <../observations/obs_converters/ROMS/ROMS.htm>`__
-  `SSEC <../observations/obs_converters/SSEC/SSEC.html>`__
-  `SST <../observations/obs_converters/SST/SST.html>`__
-  `SSUSI <../observations/obs_converters/SSUSI/convert_f16_edr_dsk.html>`__
-  `WOD <../observations/obs_converters/WOD/WOD.html>`__
-  `gnd_gps_vtec <../observations/obs_converters/gnd_gps_vtec/README>`__
-  `GPS <../observations/obs_converters/gps/gps.html>`__
-  `ok_mesonet <../observations/obs_converters/ok_mesonet/ok_mesonet.html>`__
-  `QuikSCAT <../observations/obs_converters/quikscat/QuikSCAT.html>`__
-  `Radar <../observations/obs_converters/radar/radar.html>`__
-  `snow <../observations/obs_converters/snow/snow_to_obs.html>`__
-  `Text <../observations/obs_converters/text/text_to_obs.html>`__
-  `tpw <../observations/obs_converters/tpw/tpw.html>`__
-  `Tropical Cyclones <../observations/obs_converters/tropical_cyclone/tc_to_obs.html>`__
-  `Var (little-r) <../observations/obs_converters/var/littler_tf_dart.html>`__
-  `Var (radar) <../observations/obs_converters/var/rad_3dvar_to_dart.html>`__

There are also a couple utilities of note:

-  `even_sphere <../observations/obs_converters/even_sphere/README.html>`__ - a utility for generating evenly-spaced
   observation locations that can then be used in a perfect model experiment.
-  `obs_error <../observations/obs_converters/obs_error/README.html>`__ - modules that specify observation errors
   based on what is used by ECMWF and NCEP

In addition the following external program produces DART observation sequence
files:

-  `Observation Processing And Wind Synthesis
   (OPAWS) <http://code.google.com/p/opaws/>`__: OPAWS can process NCAR Dorade
   (sweep) and NCAR EOL Foray (netCDF) radar data. It analyzes (grids) data in
   either two-dimensions (on the conical surface of each sweep) or
   three-dimensions (Cartesian). Analyses are output in netCDF, Vis5d, and/or
   DART (Data Assimilation Research Testbed) formats.

For generating synthetic observations, see the
`create_obs_sequence <../../assimilation_code/programs/create_obs_sequence/create_obs_sequence.html>`__
program documentation. You can also generate observation files based on text
input. See the `text_to_obs <text/text_to_obs.html>`__ program documentation. Or
for simulating a large complex observing system, you can use the DART library
routines in a Fortran program to compute the observation information and have
the DART routines write the output file.

See the
`perfect_model <../../assimilation_code/programs/perfect_model_obs/perfect_model_obs.html>`__
program documentation on how to run a model with a set of observations that have
only locations, types, and times, and have the forward operators compute the
observation values.
