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

# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.
# Can be used with most low-order models and the bgrid model which
# can be advanced using the integrate_model executable.

# This script copies the necessary files into the temporary directory
# and then executes the fortran program integrate_model.

# Arguments are the 
# 1) process number of caller, 
# 2) the number of state copies belonging to that process, and 
# 3) the name of the filter_control_file for that process
set process = $1
set num_states = $2
set control_file = $3

# Get unique name for temporary working directory for this process's stuff
set temp_dir = 'advance_temp'${process}

# Create a clean temporary directory and go there
\rm -rf  $temp_dir
mkdir -p $temp_dir
cd       $temp_dir

# Get files needed to run the ocean model
cp -p ../eedata ../topog.bin ../theta.bin ../salt.bin ../SST.bin ../SSS.bin .  || exit 1

# Get namelist files controlling run-time behavior
cp -p ../data ../data.cal ../data.exf ../data.kpp \
                         ../data.obcs ../data.pkg . || exit 2

# Get files needed to run DART input.nml ... is it used by trans_sv_pv ...?
cp ../input.nml . || exit 3

# Loop through each state
set state_copy = 1
set ensemble_member_line = 1
set input_file_line = 2
set output_file_line = 3
while($state_copy <= $num_states)
   
   set ensemble_member = `head -$ensemble_member_line ../$control_file | tail -1`
   set input_file      = `head -$input_file_line      ../$control_file | tail -1`
   set output_file     = `head -$output_file_line     ../$control_file | tail -1`
   
   # Get the ics file for this state_copy and
   # convert them to the form needed to cold-start the ocean model
   # trans_sv_pv  creates the following files:
   #   S.YYYYMMDD.HHMMSS.[data,meta],
   #   T.YYYYMMDD.HHMMSS.[data,meta],
   #   U.YYYYMMDD.HHMMSS.[data,meta],
   #   V.YYYYMMDD.HHMMSS.[data,meta],
   # Eta.YYYYMMDD.HHMMSS.[data,meta], and 
   # data.cal.new  ... which contains the appropriate startdate_1, startdate_2
   # so data&PARM05 will specify the input data files.
   mv ../$input_file assim_model_state_ic

   ../trans_sv_pv

   # Update the MIT namelist output ... 
   # and rename the input files to those defined in the data&PARM05 namelist.

   mv data.cal.new data.cal

   set FNAME = `grep -i hydrogSaltFile data | sed -e "s#=##"`
   set FNAME = `echo  $FNAME | sed -e "s#hydrogSaltFile##"`
   set FNAME = `echo  $FNAME | sed -e "s#,##g"`
   set FNAME = `echo  $FNAME | sed -e "s#'##g"`
   mv -v S.*.*.data $FNAME

   set FNAME = `grep -i hydrogThetaFile data | sed -e "s#=##"`
   set FNAME = `echo  $FNAME | sed -e "s#hydrogThetaFile##"`
   set FNAME = `echo  $FNAME | sed -e "s#,##g"`
   set FNAME = `echo  $FNAME | sed -e "s#'##g"`
   mv -v T.*.*.data $FNAME

   set FNAME = `grep -i uVelInitFile data | sed -e "s#=##"`
   set FNAME = `echo  $FNAME | sed -e "s#uVelInitFile##"`
   set FNAME = `echo  $FNAME | sed -e "s#,##g"`
   set FNAME = `echo  $FNAME | sed -e "s#'##g"`
   mv -v U.*.*.data $FNAME

   set FNAME = `grep -i vVelInitFile data | sed -e "s#=##"`
   set FNAME = `echo  $FNAME | sed -e "s#vVelInitFile##"`
   set FNAME = `echo  $FNAME | sed -e "s#,##g"`
   set FNAME = `echo  $FNAME | sed -e "s#'##g"`
   mv -v V.*.*.data $FNAME

   set FNAME = `grep -i hydrogSaltFile data | sed -e "s#=##"`
   set FNAME = `echo  $FNAME | sed -e "s#hydrogSaltFile##"`
   set FNAME = `echo  $FNAME | sed -e "s#,##g"`
   set FNAME = `echo  $FNAME | sed -e "s#'##g"`
   mv -v Eta.*.*.data $FNAME

#  set  hydrogSaltFile=`sed -n -e  's/hydrogSaltFile=.\(.*\).,/\1/p' data`
#  set hydrogThetaFile=`sed -n -e 's/hydrogThetaFile=.\(.*\).,/\1/p' data`
#  set    uVelInitFile=`sed -n -e    's/uVelInitFile=.\(.*\).,/\1/p' data`
#  set    vVelInitFile=`sed -n -e    's/vVelInitFile=.\(.*\).,/\1/p' data`
#  set   thetaClimFile=`sed -n -e   's/thetaClimFile=.\(.*\).,/\1/p' data`
#  mv   S.*.*.data  $hydrogSaltFile
#  mv   T.*.*.data  $hydrogThetaFile
#  mv   U.*.*.data  $uVelInitFile
#  mv   V.*.*.data  $vVelInitFile
#  mv Eta.*.*.data  $thetaClimFile

   # Advance the model saving standard out
   ./mitgcmuv >! integrate_model_out_temp

   # Extract the timestep from the ocean model output files.
   set TIMESTEP = `ls -1 S.*.data`
   set TIMESTEP = $TIMESTEP:r
   set TIMESTEP = $TIMESTEP:e

   echo $TIMESTEP | ../trans_pv_sv

   # Append the output from the advance to the file in the working directory
   #cat integrate_model_out_temp >> ../integrate_model_out_temp$process

   # Move the updated state vector back to 'centraldir'
   mv assim_model_state_ud ../$output_file

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3
end

# Change back to original directory and get rid of temporary directory
cd ..
\rm -rf $temp_dir

# Remove the filter_control file to signal completion
# Is there a need for any sleeps to avoid trouble on completing moves here?
\rm -rf $control_file

