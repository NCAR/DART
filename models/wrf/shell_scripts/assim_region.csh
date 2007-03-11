#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$

# Standard script for use in assimilation applications
# where regions are assimilated by separate executables.

# This script copies the necessary files into the temporary directory
# and then executes the fortran program assim_region.

set  CENTRALDIR = $1
set     element = $2
set    temp_dir = $3

set REMOVE = 'rm -rf'
set   COPY = 'cp -p'
set   MOVE = 'mv -f'

${REMOVE} $temp_dir
mkdir -p  $temp_dir
cd        $temp_dir

# Copy the executable, initial condition file, and other inputs to the temp directory.
# wrfinput_d0? provides auxilliary info not available from DART state vector.

${COPY} ${CENTRALDIR}/assim_region .
${MOVE} ${CENTRALDIR}/filter_assim_region__in$element filter_assim_region_in
${COPY} ${CENTRALDIR}/input.nml .
${COPY} ${CENTRALDIR}/filter_assim_obs_seq .
${COPY} ${CENTRALDIR}/wrfinput_d0? .
${COPY} ${CENTRALDIR}/namelist.input .

./assim_region >& assim_region.out

# writes out filter_assim_region_out
# for filter...assim_tools_mod/filter_assim() to read in with
# ensemble member #s tacked onto the end

# Move the updated region state to the working directory
${MOVE} filter_assim_region_out ${CENTRALDIR}/filter_assim_region_out$element

# Change back to working directory and get rid of temporary directory
cd ${CENTRALDIR}
#${REMOVE} $temp_dir
