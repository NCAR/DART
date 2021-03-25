Available observation converter programs
========================================

The **DART/observations/obs_converters** directory contains a variety of
converter programs to read various external formats and convert the observations
into the format required by DART.

Each directory has at least one converter:

-  `AIRS <AIRS/AIRS.html>`__
-  `Aviso+/CMEMS <AVISO/AVISO.html>`__
-  `Ameriflux <Ameriflux/level4_to_obs.html>`__
-  `COSMOS <COSMOS/COSMOS_to_obs.html>`__
-  `DWL <DWL/dwl_to_obs.html>`__
-  `GOES <GOES/README.md>`__
-  `GPSPW <GPSPW/README>`__
-  `GSI2DART <GSI2DART/README>`__
-  `GTSPP <GTSPP/GTSPP.html>`__
-  `MADIS <MADIS/MADIS.html>`__
-  `MIDAS <MIDAS/MIDAS_to_obs.html>`__
-  `MODIS <MODIS/MOD15A2_to_obs.htm>`__
-  `NCEP (prepbufr -> ascii) <NCEP/prep_bufr/prep_bufr.html>`__
-  `NCEP (ascii -> obs_seq) <NCEP/ascii_to_obs/create_real_obs.html>`__
-  `ROMS <ROMS/ROMS.htm>`__
-  `SSEC <SSEC/SSEC.html>`__
-  `SST <SST/SST.html>`__
-  `SSUSI <SSUSI/convert_f16_edr_dsk.html>`__
-  `WOD <WOD/WOD.html>`__
-  `cice <cice/cice_to_obs.html>`__
-  `gnd_gps_vtec <gnd_gps_vtec/README>`__
-  `GPS <gps/gps.html>`__
-  `ok_mesonet <ok_mesonet/ok_mesonet.html>`__
-  `QuikSCAT <quikscat/QuikSCAT.html>`__
-  `Radar <radar/radar.html>`__
-  `snow <snow/snow_to_obs.html>`__
-  `Text <text/text_to_obs.html>`__
-  `tpw <tpw/tpw.html>`__
-  `Tropical Cyclones <tropical_cyclone/tc_to_obs.html>`__
-  `Var (little-r) <var/littler_tf_dart.html>`__
-  `Var (radar) <var/rad_3dvar_to_dart.html>`__

There are also a couple utilities of note:

-  `even_sphere <even_sphere/README>`__ - a utility for generating evenly-spaced
   observation locations that can then be used in a perfect model experiment.
-  `obs_error <obs_error/README>`__ - modules that specify observation errors
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
