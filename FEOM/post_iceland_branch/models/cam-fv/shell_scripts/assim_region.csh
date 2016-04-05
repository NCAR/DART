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

set CENTRALDIR = $1
set element = $2
set temp_dir = $3

# Originally, we ensured temp_dir was empty, this proved to be a bit more 
# overhead and, in fact, screwed up on a GPFS system - bug reported (and fixed?)

if ( -d $temp_dir ) then
   cd   $temp_dir
else
   echo "FATAL ERROR assim_region.csh ... temp_dir( ${temp_dir} ) does not exist."
   echo "FATAL ERROR assim_region.csh ... temp_dir( ${temp_dir} ) does not exist." >> assim_region.stout

   exit 99
endif

echo "starting assim_region.csh for region $element at"`date` > assim_region.stout
echo "CENTRALDIR is $CENTRALDIR"                             >> assim_region.stout
echo "temp_dir is $temp_dir"                                 >> assim_region.stout

# Copy or link the initial conditions files and other inputs to the this
# UNIQUE temp directory -- one for each region.

if ( -s  ${CENTRALDIR}/input.nml ) then
   cp -p ${CENTRALDIR}/input.nml .
else
   # using defaults, i suppose.
   echo "WARNING assim_region.csh ... unable to copy ${CENTRALDIR}/input.nml" >> assim_region.stout
endif

if ( -s  ${CENTRALDIR}/filter_assim_obs_seq ) then
   cp -p ${CENTRALDIR}/filter_assim_obs_seq .
else
   echo "FATAL ERROR assim_region.csh ... unable to copy ${CENTRALDIR}/filter_assim_obs_seq" >> assim_region.stout
   exit 100
endif

if ( -s  ${CENTRALDIR}/caminput.nc ) then
   ln -s ${CENTRALDIR}/caminput.nc .
else
   echo "FATAL ERROR assim_region.csh ... unable to link ${CENTRALDIR}/caminput.nc" >> assim_region.stout
   exit 101
endif

if ( -s  ${CENTRALDIR}/clminput.nc ) then
   ln -s ${CENTRALDIR}/clminput.nc .
else
   echo "FATAL ERROR assim_region.csh ... unable to link ${CENTRALDIR}/clminput.nc" >> assim_region.stout
   exit 102
endif

set FILE = ${CENTRALDIR}/filter_assim_region__in$element 
if ( -s $FILE ) then
   ln -s $FILE filter_assim_region_in
else
   echo "FATAL ERROR assim_region.csh ... MISSING $FILE" >> assim_region.stout
   echo "FATAL ERROR assim_region.csh ... MISSING $FILE" >> assim_region.stout
   echo "FATAL ERROR assim_region.csh ... MISSING $FILE" >> assim_region.stout
   exit 101
endif

echo "starting assim_region at "`date` >> assim_region.stout

# assim_region writes out filter_assim_region_out
# assim_tools_mod:filter_assim() reads in 
# filter_assim_region_outXXX where the XXX is a region number
${CENTRALDIR}/assim_region            >> assim_region.stout
echo "element $element"               >> assim_region.stout
ls -lt                                >> assim_region.stout
echo "dart_log.out contents follows:" >> assim_region.stout
cat dart_log.out                      >> assim_region.stout
cat assim_region.stout                >> $CENTRALDIR/cam_reg_temp$element

# This is the piece we have been waiting for -- save it.
mv filter_assim_region_out $CENTRALDIR/filter_assim_region_out$element

cd $CENTRALDIR       ;# simply get out of this directory

# If the region output does not exist or is zero length, we save everything we can
# to a dead directory for a post-mortem. If the region is 'full', we carry on.

if (-z filter_assim_region_out$element || ! -e filter_assim_region_out$element) then
   echo "NO filter_assim_region_out$element; filter_server should stop " >> $CENTRALDIR/cam_reg_temp$element
   mkdir ${temp_dir}_dead
   mv ${temp_dir}/*  ${temp_dir}_dead
else
   \rm  $temp_dir/*    ;# clean out 'this' directory
   echo "finished assim_region.csh for region $element at "`date` >> $CENTRALDIR/cam_reg_temp$element
endif

