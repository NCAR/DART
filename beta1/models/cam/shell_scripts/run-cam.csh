#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# code added to specify an arbitrary forecast duration
# using [START,STOP]_[YMD,TOD] passed as a -namelist argument

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

set case  = $1
echo case = $case 

## ROOT OF CAM DISTRIBUTION 
# Directory which contains the CAM configuration to be used
# (resolution, optimization, etc); has files cam and config_cache.xml
# Should be full pathname, passed from advance_model.csh
# ie /gpfs/lightning/raeder/Cam3/cam3_0_7_brnchT_assim01/models/atm/cam/bld/T85-O1
set camroot  = $2
echo camroot = $camroot 

# work directory where filter is run
set CENTRALDIR  = $3
echo CENTRALDIR = $CENTRALDIR 

# ROOT OF CAM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
# Contains the initial and boundary data for the CAM distribution.
# (the root directory contains the subdirectories "atm" and "lnd")
# Up through 3.5.30 this does not need to be set if build-namelist gets an input namelist
#    argument (namelistin), which has full pathnames of all the CAM input files.
# After 3.5.30 it seems that the variable needs to be set, regardless of namelistin.
setenv CSMDATA     /fis/cgd/cseg/csm/inputdata

## $wrkdir is a working directory where the model will be built and run.
#  DART; it's the temp directory in which advance_model runs.
set wrkdir       = `pwd`
## $cfgdir is the directory containing the CAM configuration scripts.
#          and subdirectories with various configurations of this CAM version
set cfgdir       = $camroot:h

echo wrkdir  = $wrkdir
echo cfgdir  = $cfgdir

# obtain cam executable and matching config_cache.xml
if ( ! -x cam ) then
    if ( -x $camroot/cam ) then
#      camroot is in the call, so don't need to copy every time
#       cp $camroot/cam .
       cp $camroot/config_cache.xml .
       echo 'use cam and config_cache.xml from ' $camroot
    else
       echo 'cam is not found; must be pre-built and stored in ' $camroot
       exit 3
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
set verbosity = ''
set dir_arg = '-d'
set cice_nl = ''
if ($ensstring[1] == 3) then
   if ($ensstring[2] < 5) then
      set cam_version = single-namelist
      set verbosity = 2
      echo 'version < 3.5'
   else if ($ensstring[2] == 5 && $#ensstring == 3 ) then
      if ($ensstring[3] < 30) then
         set verbosity = 2
         echo 'version < 3.5.30'
      endif
      echo 'version 3.5.x'
   endif
endif
if (($ensstring[1] == 3 && $ensstring[2] > 6)   ||   $ensstring[1] >= 4 ) then
      echo 'version > 3.6.0'
      set dir_arg = '-dir'
      # Get info from job_mpi.csh/casemodel about the cice_nl for build-namelist for CAM 3.6.26
      # /ptmp/raeder/Cam3.6/Dev20-1/obs_0029/ICE/iceinput_ 1949 2007 ......./sst_file.nc
      # NO! advance_model puts the ice IC file in a generic name for this script to find
      set list = `head -5 ${CENTRALDIR}/casemodel | tail -1`
      if (-e iceinput.nc) then
         set ice_ic = iceinput.nc
      else
         set ice_ic = default
      endif
#     Single line if-test isn't bomb proof in csh; the resultant is evaluated before the conditional.
#     If resultant has a non-existent entity, an error results, even though it it won't try to use it.

      set str_yr_first  = $list[2]
      set str_yr_last   = $list[3]
      set sst           = $list[4]
      set cice_nl = "&ice     stream_year_first=$str_yr_first stream_year_last=$str_yr_last" 
      set cice_nl = "$cice_nl  model_year_align=$str_yr_first ice_ic='$ice_ic'"
      set cice_nl = "$cice_nl stream_domfilename='$sst' stream_fldfilename='$sst' /"
#     Fix the default value, which derails the gregorian calendar in ice_calendar.F90
      set cice_nl = "$cice_nl &setup_nml year_init=1 /"
      echo $cice_nl
   endif
   if (($ensstring[1] == 3 && $ensstring[2] == 5 && $ensstring[3] >= 30) || \
       ($ensstring[1] == 3 && $ensstring[2] >= 6)                        || \
       ($ensstring[1] >= 4                                                   )) then
      if (! $?CSMDATA ) then
         echo "run-cam.csh; CSMDATA must be defined in here for > Cam3.5.30"
         exit 2
      endif
      if (! -d $CSMDATA ) then
         echo "run-cam.csh; CSMDATA = $CSMDATA must exist for > Cam3.5.30"
         exit 2
      endif
   endif
endif
if ($ensstring[1] <= 2) then
   set cam_version = single-namelist
   set verbosity = 2
endif
rm ensstring.$$


echo "multi namelist? $cam_version"

#--------------------------------------------------------------
ls -lt 

# Figure out the 2D CAM domain decomposition for this number of processors.
# lat/height and lon/lat decomp pairs npr_yz == (npr_y, npr_z, nprxy_x, nprxy_y) 
# must satisfy                        nprocs =  npr_y*npr_z = nprxy_x*nprxy_y,
# and helps to have                    npr_y =  nprxy_y     and  npr_z = nprxy_x
# job_mpi.csh has calculated num_procs to make nprocs be correct,
# and the helpful condition is satisfied in the namelist below.
set length_casemodel = `wc -l ${CENTRALDIR}/casemodel`
if ($length_casemodel[1] == 8) then
   set list = `head -8 ${CENTRALDIR}/casemodel | tail -1`
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

   if ("$cice_nl" == '') then
      $cfgdir/build-namelist -v $verbosity \
        -case     ${camroot:t}-$case \
        -runtype  startup \
        $dir_arg  $wrkdir \
        -infile   $CENTRALDIR/namelistin \
        -namelist "$namelist_string" \
        || echo   "build-namelist failed" && exit 1
   else
      $cfgdir/build-namelist -v $verbosity \
        -case     ${camroot:t}-$case \
        -runtype  startup \
        $dir_arg  $wrkdir \
        -infile   $CENTRALDIR/namelistin \
        -cice_nl "$cice_nl" \
        -namelist "$namelist_string" \
        || echo   "build-namelist failed" && exit 1
      # For advance_model to copy back to CENTRALDIR for archiving; won't be used here.
      cat *_in >! namelist
   endif
endif

echo "finished build-namelist ..."

# Run CAM
# run_command is how *filter* is run, and may not be how CAM is run.
set parallel_cam = `head -6 ${CENTRALDIR}/casemodel | tail -1`
if ($parallel_cam[1] == 'true') then
   # async=4;  filter is parallel and CAM is too
   set run_command = `head -7 ${CENTRALDIR}/casemodel | tail -1`
else
   # async=2;  filter is run parallel, but CAM is not
   set run_command = ' '
endif
echo "running CAM in $wrkdir"

if ($cam_version == 'single-namelist') then
   echo "with command; eval $run_command "'"$camroot/cam < namelist"'"
#   $run_command $camroot/cam < namelist
# Bluefire
   eval $run_command "$camroot/cam < namelist"
else if ($cam_version == 'multi-namelist') then
   # CAM 3.5 explicitly opens the files?  No redirect needed?
   echo "with command $run_command $camroot/cam  (no < namelist)"
   eval $run_command $camroot/cam 
endif
ls -lt 

# Iter; these are in the run directory;
echo ' '
ls -l *\.[hir]*

# in DART caminput, clminput and maybe iceinput and a history file
# need to be saved for each "element" of the ensemble
ls *cam2\.i\.*
if ($status == 0) then
   mv *cam2\.i\.* caminput.nc
else
   echo "No cam2.i. (caminput.nc) file available from CAM; ABORTING"
   exit 4
endif

if ($cam_version == 'single-namelist') then
   ls *clm2\.i\.*
   if ($status == 0) then
      mv *clm2\.i\.*\.nc clminput.nc
   else
      echo "No clm2.i. (caminput.nc) file available from CAM; ABORTING"
      exit 5
   endif
else
   ls *clm2\.r\.*\.nc
   if ($status == 0) then
      mv *clm2\.r\.*\.nc     clminput.nc
   else
      echo "No clm2.r. (clminput.nc) file available from CAM; ABORTING"
      exit 5
   endif

   ls *\.cice\.r\.*.nc
   if ($status == 0) then
      mv *\.cice\.r\.*.nc  iceinput.nc
   else
      echo "No cice.r. (iceinput.nc) file available from CAM; ABORTING"
      exit 6
   endif


# 
# move only the youngest h0 file
   set hist = `ls -t *\.h0\.*`
   echo hist is $hist
   if ($#hist > 0) then
      mv $hist[1]    hide_hist
   endif
endif

# mv *\.h0\.* ${CENTRALDIR}/cam_phis_good.nc
# ls *\.[hir]*\.[0-9]*
# Remove 'leftovers' ONLY if CAM completed correctly (old ice initials need to
# hang around for the retry)
grep 'END OF MODEL RUN' cam_out_temp > /dev/null
if ($status == 0) then
   rm *\.[hir]*\.[0-9]*
   echo "removed .[hir]."
endif
# restore file to meaningful name advance_model expects
if (-e hide_hist ) then
   mv hide_hist $hist[1]
endif

echo ' '
ls -l ${wrkdir}/*input* ${wrkdir}/*\.h*

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

