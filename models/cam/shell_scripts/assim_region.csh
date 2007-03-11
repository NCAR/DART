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

#----------------------------------------------------------------------
# Standard script for use in assimilation applications
# where regions are assimilated by separate executables.
#
# This script copies the necessary files into the temporary directory
# and then executes the fortran program assim_region.
#----------------------------------------------------------------------

set     myname = $0
set CENTRALDIR = $1
set    element = $2
set   temp_dir = $3

# People have the craziest aliases. These prevent the obsessive-compulsive
# from causing themselves no end of angst.
if ( ! $?REMOVE ) then
  set REMOVE = 'rm -rf'
endif
if ( ! $?COPY ) then
  set COPY = 'cp -p'
endif
if ( ! $?MOVE ) then
  set MOVE = 'mv -f'
endif

echo "starting ${myname} for region $element at"`date` >! assim_region.stout
echo "CENTRALDIR is ${CENTRALDIR}"                     >> assim_region.stout
echo "temp_dir is $temp_dir"                           >> assim_region.stout

# Originally, we ensured temp_dir was empty, this proved to be a bit more 
# overhead and, in fact, screwed up on a GPFS system - bug reported (and fixed?)

if ( -d $temp_dir ) then
   cd   $temp_dir
   ${REMOVE} ${temp_dir}/*
else
   echo "FATAL ERROR ${myname} ... temp_dir( ${temp_dir} ) does not exist."
   echo "FATAL ERROR ${myname} ... temp_dir( ${temp_dir} ) does not exist." >> assim_region.stout 
   exit 99
endif

# Copy or link the initial conditions files and other inputs to the this
# UNIQUE temp directory -- one for each region.

if (  -s   ${CENTRALDIR}/input.nml ) then
   ${COPY} ${CENTRALDIR}/input.nml .
else
   # using defaults, i suppose.
   echo "WARNING ${myname} ... unable to copy ${CENTRALDIR}/input.nml" >> assim_region.stout
endif

if (  -s   ${CENTRALDIR}/filter_assim_obs_seq ) then
   ${COPY} ${CENTRALDIR}/filter_assim_obs_seq .
else
   echo "FATAL ERROR ${myname} ... unable to copy ${CENTRALDIR}/filter_assim_obs_seq" >> assim_region.stout
   exit 100
endif

if (  -s ${CENTRALDIR}/caminput.nc ) then
   ln -s ${CENTRALDIR}/caminput.nc .
else
   echo "FATAL ERROR ${myname} ... unable to link ${CENTRALDIR}/caminput.nc" >> assim_region.stout
   exit 101
endif

if ( -s  ${CENTRALDIR}/clminput.nc ) then
   ln -s ${CENTRALDIR}/clminput.nc .
else
   echo "FATAL ERROR ${myname} ... unable to link ${CENTRALDIR}/clminput.nc" >> assim_region.stout
   exit 102
endif

set FILE = ${CENTRALDIR}/filter_assim_region__in$element 
if ( -s $FILE ) then
   ln -s $FILE filter_assim_region_in
else
   echo "FATAL ERROR ${myname} ... MISSING $FILE" >> assim_region.stout
   echo "FATAL ERROR ${myname} ... MISSING $FILE" >> assim_region.stout
   echo "FATAL ERROR ${myname} ... MISSING $FILE" >> assim_region.stout
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
cat assim_region.stout                >> ${CENTRALDIR}/cam_reg_temp$element

# This is the piece we have been waiting for -- save it.
${MOVE} filter_assim_region_out ${CENTRALDIR}/filter_assim_region_out$element

cd ${CENTRALDIR}       ;# simply get out of this directory

# If the region output does not exist or is zero length, we save everything we can
# to a dead directory for a post-mortem. If the region is 'full', we carry on.

if (-z filter_assim_region_out$element || ! -e filter_assim_region_out$element) then
   set DEADDIR = ${temp_dir}_dead
   echo "WARNING - NO filter_assim_region_out$element; check $DEADDIR"
   echo "WARNING - NO filter_assim_region_out$element; check $DEADDIR"
   echo "WARNING - NO filter_assim_region_out$element; check $DEADDIR" >> ${CENTRALDIR}/cam_reg_temp$element
   echo "WARNING - NO filter_assim_region_out$element; check $DEADDIR" >> ${CENTRALDIR}/cam_reg_temp$element
   mkdir $DEADDIR
   ${MOVE} ${temp_dir}/* $DEADDIR
   exit $element
else
   ${REMOVE} ${temp_dir}/*    ;# clean out 'this' directory
   echo "finished ${myname} for region $element at "`date` >> ${CENTRALDIR}/cam_reg_temp$element
endif
