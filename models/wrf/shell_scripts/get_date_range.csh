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

#-----------------------------------------------------------------------
# Script get_date_range.csh
#
# Purpose: Set environment variables associated with date range.
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# [1] Set arguments:
#-----------------------------------------------------------------------

if ($#argv < 2) then
   echo "Usage: get_date_range <YYYYMMDDHH> <FCST_RANGE>"
   exit(1)
endif

setenv START_DATE $1
setenv FCST_RANGE $2

#-----------------------------------------------------------------------
# [2] Set environment variables:
#-----------------------------------------------------------------------

set END_DATE = `advance_cymdh $START_DATE $FCST_RANGE`

setenv START_YEAR `echo $START_DATE | cut -c1-4`
setenv START_MONTH `echo $START_DATE | cut -c5-6`
setenv START_DAY `echo $START_DATE | cut -c7-8`
setenv START_HOUR `echo $START_DATE | cut -c9-10`
setenv END_YEAR `echo $END_DATE | cut -c1-4`
setenv END_MONTH `echo $END_DATE | cut -c5-6`
setenv END_DAY `echo $END_DATE | cut -c7-8`
setenv END_HOUR `echo $END_DATE | cut -c9-10`

exit (0)
