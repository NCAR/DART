#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

echo "`date` -- BEGIN CESM ASSIMILATE"

${CASEROOT}/cam_assimilate.csh
if ( $status != 0 ) exit $status
${CASEROOT}/pop_assimilate.csh
if ( $status != 0 ) exit $status
${CASEROOT}/clm_assimilate.csh
if ( $status != 0 ) exit $status

echo "`date` -- END CESM ASSIMILATE"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

