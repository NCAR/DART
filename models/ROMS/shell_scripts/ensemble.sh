#!/bin/csh -f
#
# DART $Id$
#######################################################################
# Copyright (c) 2002-2016 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
# Script to set-up ROMS ensemble runs for data assimilation.          #
#                                                                     #
# Usage:                                                              #
#                                                                     #
#   ensemble.csh <Esize>                                              #
#                                                                     #
#######################################################################

# Set ROMS root directory.

set ROMS_ROOT=/Users/arango/ocean/repository/trunk

# Set project path to one directory up in the tree.

set ProjectDir=`dirname ${PWD}`

# Set string manipulations perl script.

set SUBSTITUTE=${ROMS_ROOT}/ROMS/Bin/substitute

# Set ensemble size, executable, input scripts, and observations.

echo
if ( ($#argv) > 0 ) then
  set Esize=$1
else
  echo "Error: you need to provide ensemble number argument..."
  exit 1
endif

set ROMS_EXE=oceanM
set ROMS_STDIN=ocean_wc13.in
set ROMS_STDIN_TMP=ocean_wc13.tmp
set ROMS_DAPAR=verification.in
set ROMS_DAPAR_TMP=s4dvar.in

set ROMS_INI=wc13_ini.nc
set ROMS_RST=wc13_rst.nc
set ROMS_OBS=wc13_obs.nc
set ROMS_MOD=wc13_mod.nc

# Set ROMS standard input parameters needed in template scripts.

set VARNAME=${ROMS_ROOT}/ROMS/External/varinfo.dat
set NtileI=1
set NtileJ=2
set NCPUS=2
set NTIMES=48
set NRST=48

# On the first pass, copy initial conditions to project directory.

cp -f -p ${ProjectDir}/Data/${ROMS_INI} ${PWD}

# Loop over all the ensemble members.

set i=1

while ( $i <= $Esize )

# Create ensemble directory, rXXX.

  set MyRun=`printf "%03d" ${i}`
  set MyDir=r${MyRun}

  echo "Processing Ensemble Run: ${MyRun}"

  if ( ! -d ${MyDir} ) then
    mkdir ${MyDir}
  endif

# Copy all required configuration files (templates) to ensemble directory.

  cp -f ${ROMS_STDIN_TMP} ${MyDir}/${ROMS_STDIN}
  cp -f ${ROMS_DAPAR_TMP} ${MyDir}/${ROMS_DAPAR}

# Copy observations to ensemble directory.

  cp -f -p ${ProjectDir}/Data/${ROMS_OBS} ${MyDir}

# Copy ROMS executable to ensemble directory.

  cp -f -p ${ROMS_EXE} ${MyDir}

# Move initial conditions from project to ensemble directory.

  mv -f ${ROMS_INI} ${MyDir}

# Go to ensemble working directory.

  cd ${MyDir}

# Set DSTART for the current ensemble.

  set TIME_SEC=`ncdump -v 'ocean_time' ${ROMS_INI}|grep "ocean_time ="|tail -1|grep -oE '[[:digit:]]+'`
  @ DSTART=$TIME_SEC / 86400

  echo "Initial Conditions time (days): ${DSTART}"

# Modify template configuration files with the approprite parameters.

 $SUBSTITUTE $ROMS_STDIN MyVARNAME $VARNAME
 $SUBSTITUTE $ROMS_STDIN MyNtileI  $NtileI
 $SUBSTITUTE $ROMS_STDIN MyNtileJ  $NtileJ
 $SUBSTITUTE $ROMS_STDIN MyNTIMES  $NTIMES
 $SUBSTITUTE $ROMS_STDIN MyNRST    $NRST
 $SUBSTITUTE $ROMS_STDIN MyDSTART  $DSTART
 $SUBSTITUTE $ROMS_STDIN MyININAME $ROMS_INI
 $SUBSTITUTE $ROMS_STDIN MyRSTNAME $ROMS_RST
 $SUBSTITUTE $ROMS_STDIN MyAPARNAM $ROMS_DAPAR

 $SUBSTITUTE $ROMS_DAPAR ocean_obs.nc $ROMS_OBS
 $SUBSTITUTE $ROMS_DAPAR ocean_mod.nc $ROMS_MOD

# Run ROMS. Add ampersand to run them all at (aproximately) the same time.

 echo "Running ROMS ensemble: ${MyRun}"

 /opt/bin/mpirunO -np ${NCPUS} ${ROMS_EXE} ${ROMS_STDIN} >& log.${MyRun}

# Generate new initial conditions from restart file in the project directory.

 cp -pv ${ROMS_RST} ${ProjectDir}/Forward/${ROMS_INI}

# Change back to project directory.

 cd ${ProjectDir}/Forward

 @ i += 1
end

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

