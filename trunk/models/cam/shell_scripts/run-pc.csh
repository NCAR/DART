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

# #sl marks streamlining of 8/31/05; whole thing is rewritten

#  code added to specify an arbitrary forecast duration
#  using [START,STOP]_[YMD,TOD] passed as a -namelist argument

#-----------------------------------------------------------------------
# PC-linux
#------------

if ($#argv == 0) then
   echo 'Usage:  run-pc.csh CASE MODEL CENTRALDIR, '
   echo '        where the namelist is in $CENTRALDIR/namelist'
   echo '        CASE and MODEL are set in long_run.csh or file casemodel'
   echo '        and CENTRALDIR is passed from advance_ens.csh'
   exit
endif

# Do our best to get sufficient stack memory
limit stacksize unlimited

# Default namelist settings:
# $case is the case identifier for this run. It will be placed in the namelist.
# $runtype is the run type: initial, restart, or branch. (and now hybrid?)
# initialize the datasets for iteration by starting with standard CAM data 
# in namelist ~raeder/CamPathRel/Caminput/${case}/namelist0, 
# (but use case = $case, no '0'; must move namelists around)

set case         = $1
set runtype      = initial

echo case = $case 

# ROOT OF CAM DISTRIBUTION 
# Directory which contains the CAM configuration to be used
# (resolution, optimization, etc); has files cam and config_cache.xml
# Should be full pathname, passed from advance_model.csh
# ie /gpfs/lightning/raeder/Cam3/cam3_0_7_brnchT_assim01/models/atm/cam/bld/T85-O1
set camroot      = $2
echo camroot = $camroot 

# work directory (on anchorage, not nodes) where filter is run
set CENTRALDIR = $3
echo CENTRALDIR = $CENTRALDIR 

# $wrkdir is a working directory where the model will be built and run.
# in DART; its the temp directory in which advance_model runs.
set wrkdir   = `pwd`

# $cfgdir is the directory containing the CAM configuration scripts.
#         and subdirectories with various configurations of this CAM version
set cfgdir   = $camroot:h

echo wrkdir  = $wrkdir
echo cfgdir  = $cfgdir

# obtain cam executable and matching config_cache.xml
if ( ! -x cam ) then
    if ( -x $camroot/cam ) then
       cp $camroot/config_cache.xml .
       echo 'use cam and config_cache.xml from ' $camroot
    else
       echo 'cam is not found; must be pre-built and stored in ' $camroot
       exit
    endif
else
   echo 'cam exists in ' $wrkdir
endif

set times = `cat $wrkdir/times`
echo run-pc times $times

# Create the namelist
echo "build-namelist ..."
ls -lt 

$cfgdir/build-namelist -v 2 -case ${camroot:t}-$case -runtype $runtype \
  -o $wrkdir/namelist -infile $CENTRALDIR/namelistin \
  -namelist "&camexp START_YMD=$times[3] START_TOD=$times[4] \
                     STOP_YMD=$times[1] STOP_TOD=$times[2] NHTFRQ=$times[5] /" \
  || echo "build-namelist failed" && exit 1

echo "finished build-namelist ..."

# Run CAM
echo "running CAM in $wrkdir"
ls -lt 

if ($4 == true) then
   echo "with command mpirun.lsf $camroot/cam < namelist"
   mpirun.lsf $camroot/cam < namelist
else
   echo "with command $camroot/cam < namelist"
   $camroot/cam < namelist 
endif 

# Iter; these are in the run directory;
ls -l *.i.*
# in DART caminput and clminput need to be saved for each "element" of the ensemble
mv *cam2.i.* $wrkdir/caminput.nc
mv *clm2.i.* $wrkdir/clminput.nc
rm *.h* *.r*
echo ' '
ls -l *input*
exit 0
