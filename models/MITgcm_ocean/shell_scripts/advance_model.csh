#!/bin/tcsh
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

# Script for use in assimilation applications
# where the model advance is executed as a separate process.

# This script copies the necessary files into the temporary directory
# and then executes the model

# Arguments are the 
# 1) process number of caller, 
# 2) the number of state copies belonging to that process, and 
# 3) the name of the filter_control_file for that process
set process = $1
set num_states = $2
set control_file = $3

echo "process      is $process"
echo "num_states   is $num_states"
echo "control_file is $control_file"

# Get unique name for temporary working directory for this process's stuff
set temp_dir = 'advance_temp'${process}
echo "temp_dir is $temp_dir"

# Create a clean temporary directory and go there
\rm -rf  $temp_dir
mkdir -p $temp_dir
cd       $temp_dir

# Copy the namelist files - these are small, so we really copy them
foreach FILE ( data data.cal data.exf data.kpp \
               data.obcs data.pkg eedata )
   cp -pv ../inputs/$FILE .
end

# copy the files used by data&PARM05 - input datasets
# These get overwritten ... so maybe we dont actually copy them
foreach FILE ( bathymetry.bin gom_H_199601.bin gom_S_199601.bin \
             gom_T_199601.bin gom_U_199601.bin gom_V_199601.bin )
#   cp -pv ../inputs/$FILE .
end

# link the files used by data.exf&EXF_NML_02 - external forcings
foreach FILE ( lev05_monthly_sss_relax.bin \
               lev05_monthly_sst_relax.bin \
               run-off.bin_1x1             \
               ncep_air_19960101.bin       \
               ncep_dlwrf_19960101.bin     \
               ncep_dswrf_19960101.bin     \
               ncep_nswrs_19960101.bin     \
               ncep_prate_19960101.bin     \
               ncep_shum_19960101.bin      \
               ncep_uwnd_19960101.bin      \
               ncep_vwnd_19960101.bin      )
   ln -sf ../inputs/$FILE .
end

# link the files used by data.obcs&OBCS_PARM01 - open boundaries
foreach FILE ( Rs_SobcsE_52_01_nPx1.bin    Rs_SobcsN_52_01_nPy1.bin \
               Rs_TobcsE_52_01_nPx1.bin    Rs_TobcsN_52_01_nPy1.bin \
               Rs_UobcsE_52_01_nPx1_c1.bin Rs_UobcsN_52_01_nPy1.bin \
               Rs_VobcsE_52_01_nPx1.bin    Rs_VobcsN_52_01_nPy1_c1.bin)
   ln -sf ../inputs/$FILE .
end


# Get files needed to run DART
cp ../input.nml .

echo 'listing now that the table has been set ...'
ls -l

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

   mv -v ../$input_file assim_model_state_ic

   ../trans_sv_pv

   # Update the MIT namelist output ... 
   # and rename the input files to those defined in the data&PARM05 namelist.
   # FIXME - somehow the last ensemble member needs to copy the 
   # data.cal.new back to the CENTRALDIR ... I think

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

   set FNAME = `grep -i pSurfInitFile data | sed -e "s#=##"`
   set FNAME = `echo  $FNAME | sed -e "s#pSurfInitFile##"`
   set FNAME = `echo  $FNAME | sed -e "s#,##g"`
   set FNAME = `echo  $FNAME | sed -e "s#'##g"`
   mv -v Eta.*.*.data $FNAME

# set pattern="s/.*'\(.*\)'.*/\1/"
# set zfilename=`echo $zoneline |  sed -e $pattern`
#
#  set  hydrogSaltFile=`sed -n -e  's/hydrogSaltFile=.\(.*\).,/\1/p' data`
#  set hydrogThetaFile=`sed -n -e 's/hydrogThetaFile=.\(.*\).,/\1/p' data`
#  set    uVelInitFile=`sed -n -e    's/uVelInitFile=.\(.*\).,/\1/p' data`
#  set    vVelInitFile=`sed -n -e    's/vVelInitFile=.\(.*\).,/\1/p' data`
#  set   pSurfInitFile=`sed -n -e   's/pSurfInitFile=.\(.*\).,/\1/p' data`
#  mv   S.*.*.data  $hydrogSaltFile
#  mv   T.*.*.data  $hydrogThetaFile
#  mv   U.*.*.data  $uVelInitFile
#  mv   V.*.*.data  $vVelInitFile
#  mv Eta.*.*.data  $pSurfInitFile

   # Advance the model saving standard out
   mpirun.lsf ../mitgcmuv_20p

   # Extract the timestep from the ocean model output files.
   set TIMESTEP = `ls -1 S.*.data`
   set TIMESTEP = $TIMESTEP:r
   set TIMESTEP = $TIMESTEP:e

   echo $TIMESTEP | ../trans_pv_sv

   # Move the updated state vector back to 'centraldir'
   mv assim_model_state_ud ../$output_file

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3
end

# Change back to original directory and get rid of temporary directory
ls -l > ../directory_contents_$ensemble_member
cd ..
\rm -rf $temp_dir

# Remove the filter_control file to signal completion
# Is there a need for any sleeps to avoid trouble on completing moves here?
\rm -rf $control_file

