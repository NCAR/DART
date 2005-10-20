#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$
#
#-----------------------------------------------------------------------
# Script run_wrfreal.csh
#
# Purpose: Run WRFV2 real.exe.
#
#-----------------------------------------------------------------------

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

exit (0)
