#!/bin/csh
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

#----------------------------------------------------------------------
# advance_model.csh
#
# Script to advance one ensemble member one filter "time step"
# when the model advance is executed as a separate process.
# Called by filter_server.csh.
# Calls run-pc.csh, the CAM execution script.
# Calls 3 translation routines to translate time and model state.
# Runs on one of the compute nodes allotted to filter_server.csh.
#
# Has either 3 or 4 input arguments.
# arg#1  is the name of the CENTRALDIR
# arg#2  is the ensemble member number
# arg#3  is the name of the directory that is used to advance this particular 
#        instance of the model. Each ensemble member (i.e. each execution of CAM)
#        runs in its own, unique directory that is emptied upon entry.
# arg#4  is optional, and contains the 'machine file' for multi-threaded 
#        instances of CAM. 
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

echo "starting ${myname} for ens member $element at "`date` >! cam_out_temp
echo "CENTRALDIR is ${CENTRALDIR}"                          >> cam_out_temp
echo "temp_dir is $temp_dir"                                >> cam_out_temp

# Originally, we ensured temp_dir was empty, this proved to be a bit more 
# overhead and, in fact, screwed up on a GPFS system - bug reported (and fixed?)

if ( -d $temp_dir ) then
   cd   $temp_dir
   ${REMOVE} ${temp_dir}/*
   echo "advance_model; cd to an empty $temp_dir" >> cam_out_temp
else
   mkdir -p ${temp_dir}
   cd       ${temp_dir}
#  echo "FATAL ERROR ${myname} ... temp_dir( ${temp_dir} ) does not exist."
#  echo "FATAL ERROR ${myname} ... temp_dir( ${temp_dir} ) does not exist." >> cam_out_temp
#  exit 99
endif

# Get information about this experiment from file "casemodel", 
# created by the main controlling script (job.csh)
set caseinfo = `cat ${CENTRALDIR}/casemodel`

# set $case = the case and 
# $model = the directory name (in the central CAM directory) 
# where CAM executable will be found. 
# set locations of the CAM and CLM input files
set     case = $caseinfo[1]
set    model = $caseinfo[2]
set cam_init = $caseinfo[3]
set clm_init = $caseinfo[4]

# output diagnostic information to the same file as the CAM list-directed output
echo "case $case model $model"    >> cam_out_temp
echo "cam init is $cam_init"      >> cam_out_temp
echo "clm init is $clm_init"      >> cam_out_temp

# get model state initial conditions for this ensemble member
ln -s ${CENTRALDIR}/assim_model_state_ic$element temp_ic

# get filter namelists for use by 
${COPY} ${CENTRALDIR}/input.nml input.nml

# this just creates a file that helps you figure out which member is
# being advanced in this directory. FYI only, you don't need it.
echo "junk" >! element$element

echo "ls $temp_dir for element $element" >> cam_out_temp
ls -lRt                                  >> cam_out_temp
 
# Need a base CAM initial file into which to copy state vector from filter.
# c[al]minput_$element also carry along CAM/CLM fields which are not updated
#      by the filter (not part of the filter model state).
# First look for c[al]minput.nc resulting from the previous advance of this ensemble
#      member from within the same day/obs_seq.out time span (in CENTRALDIR)
# Failing that, look for the results of the last advance of this ensemble member
#      of the previous obs_seq.out (i.e. in CENTRALDIR/exp_name/day/CAM)
# Failing that (when starting an experiment which has no spun up set of members)
#      get a copy of a single CAM initial file (usually from somewhere independent
#      of this experiment, i.e. /scratch/.../New_state/T42_GWD/CAM/caminput_0.nc)

if (-e     ${CENTRALDIR}/caminput_${element}.nc) then
   ${COPY} ${CENTRALDIR}/caminput_${element}.nc caminput.nc
   echo "CENTRALDIR caminput comes from ${CENTRALDIR}/caminput_${element}.nc" >> cam_out_temp
else if (-e ${cam_init}${element}.nc) then
   ${COPY}  ${cam_init}${element}.nc caminput.nc
   echo "cam_init caminput comes from ${cam_init}${element}.nc" >> cam_out_temp
else
   ${COPY}  ${cam_init}0.nc caminput.nc
   echo "DEFAULT caminput comes from ${cam_init}0.nc" >> cam_out_temp
endif

if ( -e     ${CENTRALDIR}/clminput_${element}.nc) then
    ${COPY} ${CENTRALDIR}/clminput_${element}.nc clminput.nc
else if (-e ${clm_init}${element}.nc) then
    ${COPY} ${clm_init}${element}.nc clminput.nc
else
    ${COPY} ${clm_init}0.nc clminput.nc
endif

# create 'times' file for CAM from DART times in assim_model_state_ic#
# This info is passed to CAM through the creation of its namelist
if (-e temp_ic && -e ${CENTRALDIR}/trans_time) then
   echo 'advance_model; executing trans_time '`date` >> cam_out_temp
   ${CENTRALDIR}/trans_time
   ls -lt >> cam_out_temp
   ${COPY} times ${CENTRALDIR}
else
   echo "ERROR: either ic file $element or trans_time not available for trans_time"
   exit 1
endif

echo "before trans_sv_pv"`ncdump -h caminput.nc | grep -i phis`

# Create an initial CAM.nc file from the DART state vector
# Times are handled separately in trans_time
${CENTRALDIR}/trans_sv_pv
ls -ltR >> cam_out_temp

echo "before run-pc.csh "`ncdump -h caminput.nc | grep -i phis`

# advance cam 
# 'machine_file' is NOT USED by run-pc.csh in single-threaded executions 
if ($#argv == 4) then
   set machine_file = $4
   ${model:h}/run-pc.csh ${case}-$element $model ${CENTRALDIR} $machine_file >>& cam_out_temp
else
   ${model:h}/run-pc.csh ${case}-$element $model ${CENTRALDIR} >>& cam_out_temp
endif

echo "after run-pc.csh "`ncdump -h caminput.nc | grep -i phis`

grep 'END OF MODEL RUN' cam_out_temp > /dev/null
if ($status == 0) then
   # Extract the new state vector information from the new caminput.nc and
   # put it in temp_ud (time followed by state)
   ${CENTRALDIR}/trans_pv_sv

   echo "after trans_pv_sv "`ncdump -h caminput.nc | grep -i phis`

   # Move updated state vector and new CAM/CLM initial files back to experiment
   # directory for use by filter and the next advance.
   ${MOVE} temp_ud     ${CENTRALDIR}/assim_model_state_ud$element
   ${MOVE} clminput.nc ${CENTRALDIR}/clminput_${element}.nc
   ${MOVE} caminput.nc ${CENTRALDIR}/caminput_${element}.nc
   ${MOVE} namelist    ${CENTRALDIR}

   echo "finished ${myname} for ens member $element at "`date` >> cam_out_temp
   ${MOVE} cam_out_temp ${CENTRALDIR}/cam_out_temp$element
   cd ${CENTRALDIR}

   # lightnings filesystem was misbehaving if you removed the directory.
   # it was more reliable if you just left the directory empty. It's bogus, but true.
   ${REMOVE} $temp_dir/*
else
   set DEADDIR = ${temp_dir}_dead
   echo "WARNING - CAM $element stopped abnormally; might be retried, see $DEADDIR"
   echo "WARNING - CAM $element stopped abnormally; might be retried, see $DEADDIR"
   echo "WARNING - CAM $element stopped abnormally; might be retried, see $DEADDIR" >> cam_out_temp
   echo "WARNING - CAM $element stopped abnormally; might be retried, see $DEADDIR" >> cam_out_temp
   ${COPY} cam_out_temp ${CENTRALDIR}/cam_out_temp$element
   mkdir $DEADDIR
   ${MOVE} * $DEADDIR
   exit $element
endif
