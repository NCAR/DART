#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

# if you are skipping the POP or CLM assimilations you can simply
# comment out the calls to the xxx_assimilate.csh scripts.  however
# if you are skipping the CAM assimilation you must replace the
# call below with a call to ${CASEROOT}/cam_no_assimilate.csh

echo "`date` -- BEGIN CESM ASSIMILATE"

${CASEROOT}/cam_assimilate.csh
if ( $status != 0 ) exit $status
${CASEROOT}/pop_assimilate.csh
if ( $status != 0 ) exit $status
${CASEROOT}/clm_assimilate.csh
if ( $status != 0 ) exit $status

echo "`date` -- END CESM ASSIMILATE"

exit 0


