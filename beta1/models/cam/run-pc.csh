#!/bin/tcsh -f
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
#  #td means code added to specify an arbitrary forecast duration
#      using [START,STOP]_[YMD,TOD] passed as a -namelist argument

#-----------------------------------------------------------------------
## PC-linux
##------------

if ($#argv == 0) then
   echo 'Usage:  run-pc.csh CASE MODEL PBS_O_WORKDIR, '
   echo '        where the namelist is in ~raeder/CamPathRel/Caminput/CASE/namelist'
   echo '        CASE and MODEL are set in setup_advance_model.csh'
   echo '        and PBS_O_WORKDIR is passed from advance_ens.csh'
   exit
endif

## nthreads is the number of Open-MP threads.  On the PC we assume a 
## pure shared-memory configuration which assigns 1 thread per processor.
set nthreads = 2

## Do our best to get sufficient stack memory
limit stacksize unlimited

## ROOT OF CAM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CAM distribution.
## (the root directory contains the subdirectory "models")
# original; set camroot      = /fs/cgd/data0/$LOGNAME/cam2.0.1
# PWD should be tempdir#, which resides under cam_test, along with cam2.0.1
set camroot      = $2
# set camroot      = $model

## Default namelist settings:
## $case is the case identifier for this run. It will be placed in the namelist.
## $runtype is the run type: initial, restart, or branch.
# initialize the datasets for iteration by starting with standard CAM data 
# in namelist ~raeder/CamPathRel/Caminput/${case}/namelist0, 
# (but use case = $case, no '0'; must move namelists around)

set case         = $1
set runtype      = initial
set relocate     = ' '
# set relocate     = "-fflags '-Kpic -Bstatic'"

# work directory (on anchorage, not nodes) where filter is run
set PBS_O_WORKDIR = $3

## ROOT OF CAM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CAM distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
# Orig setenv CSMDATA     /fs/cgd/csm/inputdata
# in DART caminput and clminput need to be saved for each "element" of
# the ensemble, but only when running them in parallel.
# For now keep caminput.nc and clminput.nc in 'cam_test'.
# not 'cam_test/tempdir#'

setenv CSMDATA     $PWD

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CAM configuration scripts.
set wrkdir       = $PWD
set blddir       = $wrkdir/bld
set rundir       = $wrkdir
set cfgdir       = $camroot/models/atm/cam/bld

## Ensure that run and build directories exist
# mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## If an executable doesn't exist, build one.
if ( ! -x $blddir/cam ) then
    if ( -x $PBS_O_WORKDIR/$cfgdir/cam ) then
       cp $PBS_O_WORKDIR/$cfgdir/cam $blddir/cam
       cp $PBS_O_WORKDIR/$cfgdir/config_cache.xml $blddir/config_cache.xml
    else
# This will be somewhat inefficient when running on multiple nodes;
# recompiling for each element of the first batch
# But this will be fixed (most of the time) when iteration 0 run of cam
# is installed in advance_init.csh(?)
       cd $blddir                  || echo "cd $blddir failed" && exit 1
       echo 'configuring CAM for case '$case 
#      config_cache_defaults is hard-wired into configure; -defaults doesn't change it
#         so copy the defaults file I want into the name configure needs
       $PBS_O_WORKDIR/$cfgdir/configure -nosmp $relocate  \
            || echo "configure failed" && exit 1
       echo "building CAM in $blddir ..." 
       rm -f Depends
       gmake -j2 >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
       echo "finished gmake of bld ..."
# store compiled cam for other "elements" in cam_async.csh
       cp $blddir/cam $PBS_O_WORKDIR/$cfgdir/cam
       cp $blddir/config_cache.xml $PBS_O_WORKDIR/$cfgdir/config_cache.xml
    endif
endif

#td added -namelist argument
## read in start and end dates and times
set times = `cat $PBS_O_WORKDIR/times`
echo $times

## Create the namelist
cd $blddir                      || echo "cd $blddir failed" && exit 1
echo "build-namelist ..."
$PBS_O_WORKDIR/$cfgdir/build-namelist -v 2 -case $case -runtype $runtype \
  -o $rundir/namelist -infile $PBS_O_WORKDIR/namelistin \
  -namelist "&camexp START_YMD=$times[3] START_TOD=$times[4] \
                     STOP_YMD=$times[1] STOP_TOD=$times[2] NHTFRQ=$times[5] /" \
  || echo "build-namelist failed" && exit 1
cp $rundir/namelist /scratch/cluster/raeder/SPMD/namelisttemp
#td added -namelist argument
# before move to CamPathRel;  -infile ../../$cfgdir/$case.nl \
# my namelist input   
echo "finished build-namelist ..."

## Run CAM
cd $rundir                      || echo "cd $rundir failed" && exit 1
echo "running CAM in $rundir"
env OMP_NUM_TASKS=$nthreads MPSTKZ="128M" $blddir/cam < namelist  || echo "CAM run failed" && exit 1

# Iter; these are in the run directory;
ls -l *.i.*
mv *cam2.i.* $wrkdir/caminput.nc
mv *clm2.i.* $wrkdir/clminput.nc
rm *.h* *.r*
echo ' '
ls -l *input*

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

