#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL: $
# $Id: $
# $Revision: $
# $Date: $

# #sl marks streamlining of 8/31/05; whole thing is rewritten

#  code added to specify an arbitrary forecast duration
#  using [START,STOP]_[YMD,TOD] passed as a -namelist argument

#-----------------------------------------------------------------------
## PC-linux
##------------

if ($#argv == 0) then
   echo 'Usage:  run-cam.csh CASE MODEL CENTRALDIR, '
   echo '        where '
   echo '        CASE and MODEL are set in job_mpi.csh or file casemodel'
   echo '        and CENTRALDIR is passed from advance_model.csh'
   exit
endif

## Do our best to get sufficient stack memory
limit stacksize unlimited

## Default namelist settings:
## $case is the case identifier for this run. It will be placed in the namelist.
# initialize the datasets for iteration by starting with standard CAM data 
# in namelist ~raeder/CamPathRel/Caminput/${case}/namelist0, 
# (but use case = $case, no '0'; must move namelists around)

set case         = $1

echo case = $case 

# sl
## ROOT OF CAM DISTRIBUTION 
# Directory which contains the CAM configuration to be used
# (resolution, optimization, etc); has files cam and config_cache.xml
# Should be full pathname, passed from advance_model.csh
# ie /gpfs/lightning/raeder/Cam3/cam3_0_7_brnchT_assim01/models/atm/cam/bld/T85-O1
set camroot      = $2
echo camroot = $camroot 

# work directory (on anchorage, not nodes) where filter is run
# sl; changed name from PBS_O_WORKDIR
set CENTRALDIR = $3
echo CENTRALDIR = $CENTRALDIR 

## ROOT OF CAM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CAM distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
# Orig setenv CSMDATA     /fs/cgd/csm/inputdata
# This works for me for some reason;  setenv CSMDATA     $PWD
# sl
# This does not need to be set because build-namelist gets an input namelist
#    argument, which has full pathnames of all the CAM input files.

# sl
## $wrkdir is a working directory where the model will be built and run.
#  DART; it's the temp directory in which advance_model runs.
# set wrkdir       = $PWD       
set wrkdir       = `pwd`
## $cfgdir is the directory containing the CAM configuration scripts.
#          and subdirectories with various configurations of this CAM version
set cfgdir       = $camroot:h

echo wrkdir  = $wrkdir
echo cfgdir  = $cfgdir

# obtain cam executable and matching config_cache.xml
if ( ! -x cam ) then
    if ( -x $camroot/cam ) then
#       cp $camroot/cam .
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
echo run-cam times $times
# Can't use \ in the middle of a string and then add more to the string later.
# The \ doesn't appear in the string, but the "" become unbalanced.
set namelist_string =         "&camexp START_YMD=$times[3] START_TOD=$times[4] "
set namelist_string = "$namelist_string STOP_YMD=$times[1]  STOP_TOD=$times[2] NHTFRQ=$times[5] " 

#--------------------------------------------------------------
## Create the namelist
# Extract the relevant line from input.nml
grep model_version input.nml >! ensstring.$$

# Replace the ,s with nothings and put the words into a list
set  STRING = "1,$ s#,##g"
set ensstring = `sed -e "$STRING" ensstring.$$`
# Output the cam version to a file
echo $ensstring[3] >! ensstring.$$

# Replace the 's with nothings and put the words into a list
set  STRING = "1,$ s#'##g"
set ensstring = `sed -e "$STRING" ensstring.$$`
# Output the cam version to a file
echo $ensstring >! ensstring.$$

echo "build-namelist for cam $ensstring ..."

# replace the .s with spaces and put the resulting numbers into a string
set  STRING = "1,$ s#\.# #g"
set ensstring = `sed -e "$STRING" ensstring.$$`

set cam_version = multi-namelist
if (($ensstring[1] == 3) && ($ensstring[2] < 5)) set cam_version = single-namelist
if  ($ensstring[1] <= 2)                         set cam_version = single-namelist

#--------------------------------------------------------------
ls -lt 

# Figure out the 2D CAM domain decomposition for this number of processors.
# lat/height and lon/lat decomp pairs npr_yz == (npr_y, npr_z, nprxy_x, nprxy_y) 
# must satisfy                        nprocs =  npr_y*npr_z = nprxy_x*nprxy_y,
# and helps to have                    npr_y =  nprxy_y     and  npr_z = nprxy_x
# job_mpi.csh has calculated num_procs to make nprocs be correct,
# and the helpful condition is satisfied in the namelist below.
set length_casemodel = `wc -l ${CENTRALDIR}/casemodel`
if ($length_casemodel[1] == 7) then
   set list = `head -7 ${CENTRALDIR}/casemodel | tail -1`
   set num_procs  = $list[1]
   set lev_blocks = $list[2]
   set lat_blocks = $list[3]
   set lon_blocks = $lev_blocks 
   set namelist_string = \
       "$namelist_string npr_yz=$lat_blocks, $lev_blocks, $lon_blocks, $lat_blocks /" 
else
   set namelist_string = "$namelist_string /" 
endif


if ($cam_version == 'single-namelist') then
   $cfgdir/build-namelist -v 2 \
     -case     ${camroot:t}-$case \
     -runtype  initial \
     -o        $wrkdir/namelist \
     -infile   $CENTRALDIR/namelistin \
     -namelist "$namelist_string" \
     || echo   "build-namelist failed" && exit 1
else if ($cam_version == 'multi-namelist') then
   # This builds all the *_in namelists CAM3.5 needs

   $cfgdir/build-namelist -v 2 \
     -case     ${camroot:t}-$case \
     -runtype  startup \
     -d        $wrkdir \
     -infile   $CENTRALDIR/namelistin \
     -namelist "$namelist_string" \
     || echo   "build-namelist failed" && exit 1
     # For advance_model to copy back to CENTRALDIR for archiving; won't be used here.
     cat *_in >! namelist
endif

echo "finished build-namelist ..."

# Run CAM
# run_command is how *filter* is run, and may not be how CAM is run.
set parallel_cam = `head -5 ${CENTRALDIR}/casemodel | tail -1`
if ($parallel_cam == true) then
   # async=4;  filter is parallel and CAM is too
   set run_command = `head -6 ${CENTRALDIR}/casemodel | tail -1`
else
   # async=2;  filter is run parallel, but CAM is not
   set run_command = ' '
endif
echo "running CAM in $wrkdir"

if ($cam_version == 'single-namelist') then
   echo "with command $run_command $camroot/cam < namelist"
   $run_command $camroot/cam < namelist
else if ($cam_version == 'multi-namelist') then
   # CAM 3.5 explicitly opens the files?  No redirect needed?
   echo "with command $run_command $camroot/cam  (no < namelist)"
   $run_command $camroot/cam 
endif
ls -lt 

# Iter; these are in the run directory;
ls -l *.i.*
# in DART caminput and clminput need to be saved for each "element" of the ensemble
mv *cam2.i.* $wrkdir/caminput.nc
mv *clm2.i.* $wrkdir/clminput.nc
rm *.h* *.r*
echo ' '
ls -l *input*
exit 0
