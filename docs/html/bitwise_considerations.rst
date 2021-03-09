Bitwise Considerations
======================

By bitwise we mean bit for bit identical results in the output obs_sequence file and the restarts (netcdf-to-netcdf or
DART format-to-dart format) when comparing one version of the code to another. For testing the code to be bitwise with
Lanai there are several things to change/set in Manhattan and your mkmf:

Important
~~~~~~~~~

The *CAM* and *bgrid_solo* model_mods have been altered so the state is in a different order inside filter. Thus DART
format restarts will **not** be bitwise with Lanai DART format restarts, but netcdf files will be (after running
dart_to_cam).
