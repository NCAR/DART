3DVAR/4DVAR Observation Converters
==================================

Overview
--------

The programs in the ``obs_converters/var`` directory help convert data which is 
formatted for input into the 3DVAR/4DVAR programs into DART
obs_seq observation files.

The directory contains conversion programs for various
obs formats related to 3D-Var, WRF-Var, and MM5:

- :doc:`./littler_tf_dart` to and back from little-r format, temperature and winds only.
- :doc:`./rad_3dvar_to_dart` the radar 3d-var obs only to dart format.
- ``gts_to_dart.f90`` from GTS to dart format.

You need to add some WRF-Var source files to the 3DVAR_OBSPROC
directory, and then you can go into the work directory and
run the 'quickbuild.sh' script. The required WRF-Var source files are
listed in ``3DVAR_OBSPROC/README``.

The little-r converter may need changes to the code to convert
from the original quality control flags into QC flags compatible
with DART.  (in DART, 0 is good data.)

The GTS converter does not support SATEM thickness data but
there are versions around which do; write dart@ucar.edu if you
are interested in more about this.

And a final disclaimer:
Whether these work with the latest 3D-Var format is untested
at this point.  Please contact the DART Development group if
you are interested in using these tools.

