#!/bin/csh -f
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$
#
# Standard script for use in assimilation applications
# where regions are assimilated by separate executables.

# This script copies the necessary files into the temporary directory
# and then executes the fortran program assim_region.

set WORKDIR = $1
set element = $2
set temp_dir = $3

set REMOVE = 'rm -rf'
set   COPY = 'cp -p'
set   MOVE = 'mv -f'

${REMOVE} $temp_dir
mkdir -p  $temp_dir
cd        $temp_dir

# Copy the executable, initial condition file, and other inputs to the temp directory.
# wrfinput_d0? provides auxilliary info not available from DART state vector.

${COPY} ${WORKDIR}/assim_region .
${MOVE} ${WORKDIR}/filter_assim_region__in$element filter_assim_region_in
${COPY} ${WORKDIR}/input.nml .
${COPY} ${WORKDIR}/filter_assim_obs_seq .
${COPY} ${WORKDIR}/wrfinput_d0? .
${COPY} ${WORKDIR}/namelist.input .

./assim_region >& assim_region.out

# writes out filter_assim_region_out
# for filter...assim_tools_mod/filter_assim() to read in with
# ensemble member #s tacked onto the end

${MOVE} filter_assim_region_out $WORKDIR/filter_assim_region_out$element

cd $WORKDIR
#${REMOVE} $temp_dir
