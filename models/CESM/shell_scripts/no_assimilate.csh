#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# The DART sourcemods for CAM require that the cam 'initial' files are used
# instead of the restart files. The name of the initial file is read from
# the cam namelist file. To avoid having to change the namelist setting, DART
# uses a static filename in the namelist. After each successful model advance,
# the new initial file must be linked to the static file name.
# This is done in ${CASEROOT}/cam_no_assimilate.csh
#
# No action is required for the other model components.

echo "`date` -- BEGIN CESM-DART NO-ASSIMILATE"

${CASEROOT}/cam_no_assimilate.csh
if ( $status != 0 ) exit $status

echo "`date` -- END CESM-DART NO-ASSIMILATE"

exit 0


