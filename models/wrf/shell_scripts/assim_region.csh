#!/bin/csh -f
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
# Standard script for use in assimilation applications
# where regions are assimilated by separate executables.

# This script copies the necessary files into the temporary directory
# and then executes the fortran program assim_region.

set WORKDIR = $1
set element = $2
set temp_dir = $3

rm -rf   $temp_dir
mkdir -p $temp_dir
cd       $temp_dir

# Copy the executable, initial condition file, and other inputs to the temp directory

cp -pv ${WORKDIR}/assim_region .
mv -v  ${WORKDIR}/filter_assim_region__in$element filter_assim_region_in
cp -pv ${WORKDIR}/input.nml .
cp -pv ${WORKDIR}/filter_assim_obs_seq .
cp -pv ${WORKDIR}/wrfinput_d0? .
                   # Provides auxilliary info not avail. from DART state vector
cp -pv ${WORKDIR}/namelist.input .

./assim_region >& assim_region.out
# writes out filter_assim_region_out
# for filter...assim_tools_mod/filter_assim() to read in with
# ensemble member #s tacked onto the end

mv -v filter_assim_region_out $WORKDIR/filter_assim_region_out$element

cd $WORKDIR
#rm -rf $temp_dir
