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
# Standard script for use in assimilation applications
# where regions are assimilated by separate executables.

# This script copies the necessary files into the temporary directory
# and then executes the fortran program assim_region.

set WORKDIR = $1
set element = $2
set temp_dir = $3

\rm -rf  $temp_dir
mkdir -p $temp_dir
cd       $temp_dir

# Copy the executable, initial condition file, and other inputs to the temp directory

cp ${WORKDIR}/assim_region .
cp ${WORKDIR}/filter_assim_region__in$element filter_assim_region_in
cp ${WORKDIR}/input.nml .
cp ${WORKDIR}/filter_assim_obs_seq .
cp ${WORKDIR}/caminput.nc .
cp ${WORKDIR}/clminput.nc .


# assim_region writes out filter_assim_region_out
# assim_tools_mod:filter_assim() reads in 
# filter_assim_region_outXXX where the XXX is a region number
./assim_region           > assim_region.stout
echo "element $element" >> assim_region.stout
ls -lt                  >> assim_region.stout
cat assim_region.stout  >> $WORKDIR/cam_reg_temp$element

mv filter_assim_region_out $WORKDIR/filter_assim_region_out$element

cd $WORKDIR          ;# simply get out of this directory
\rm -rf $temp_dir    ;# clean out 'this' directory
