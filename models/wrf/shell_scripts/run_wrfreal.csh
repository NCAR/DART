#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# Purpose: Run WRFV2 real.exe.

#-----------------------------------------------------------------------
# [1] Required environment variables:
#-----------------------------------------------------------------------

setenv NC                 $1
setenv WRF_DIR            $2

#-----------------------------------------------------------------------
# [2] Run real.exe:
#-----------------------------------------------------------------------

cd ${WRF_DIR}/test/em_real

echo "   Running real.exe"
real.exe >>& ./run_real_${NC}.out

mv namelist.input namelist.input_${NC}
mv rsl.out.0000 rsl.out.0000_${NC}
mv rsl.error.0000 rsl.error.0000_${NC}

echo ""

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

