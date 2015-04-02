#!/bin/tcsh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This script has 4 logical 'blocks':
# 1) creates a clean, temporary directory in which to run a model instance 
#    and copies the necessary files into the temporary directory
# 2) converts the DART output to input expected by the ocean model
# 3) runs the ocean model
# 4) converts the ocean model output to input expected by DART
#
# The error code from the script reflects which block it failed.
#
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

#-------------------------------------------------------------------------
# Block 1: populate a run-time directory with the bits needed to 
# run the ocean model.
#-------------------------------------------------------------------------

# Get unique name for temporary working directory for this process's stuff
set temp_dir = 'advance_temp'${process}
echo "temp_dir is $temp_dir"

# Create a clean temporary directory and go there
\rm -rf  $temp_dir
mkdir -p $temp_dir
cd       $temp_dir

# Copy the unchanging namelist files - these are small, so we really copy them
foreach FILE ( data.exf data.kpp data.obcs data.pkg eedata )
   cp -pv ../inputs/$FILE . || exit 1
end

# Get the 'changing' namelist files from CENTRALDIR
# Only the namelists in CENTRALDIR have the updated information about
# the state of the model partway through an assimilation experiment.
foreach FILE ( data data.cal input.nml )
   cp -pv ../$FILE . || exit 1
end

# copy the files used by data&PARM05 - input datasets
# These get overwritten ... so maybe we dont actually copy them
#foreach FILE ( bathymetry.bin gom_H_199601.bin gom_S_199601.bin \
#             gom_T_199601.bin gom_U_199601.bin gom_V_199601.bin )
foreach FILE ( bathymetry.bin )
    cp -pv ../inputs/$FILE . || exit 1
end

# link the files used by data.exf&EXF_NML_02 - external forcings
foreach FILE ( lev05_monthly_sss_relax.bin \
               lev05_monthly_sst_relax.bin \
               run-off.bin_1x1             \
               ncep_air_2008.bin      \
               ncep_dlwrf_2008.bin    \
               ncep_dswrf_2008.bin    \
               ncep_nswrs_2008.bin    \
               ncep_prate_2008.bin    \
               ncep_shum_2008.bin     \
               ncep_uwnd_2008.bin     \
               ncep_vwnd_2008.bin      )
   ln -sf ../inputs/$FILE . || exit 1
end

# link the files used by data.obcs&OBCS_PARM01 - open boundaries
foreach FILE ( Rs_SobcsE_2008_7DAY_GOM_Jan01_nPx1.bin Rs_SobcsN_2008_7DAY_GOM_Jan01_nPy1.bin \
               Rs_TobcsE_2008_7DAY_GOM_Jan01_nPx1.bin Rs_TobcsN_2008_7DAY_GOM_Jan01_nPy1.bin \
               Rs_UobcsE_2008_7DAY_GOM_Jan01_nPx1_c1.bin Rs_UobcsN_2008_7DAY_GOM_Jan01_nPy1.bin \
               Rs_VobcsE_2008_7DAY_GOM_Jan01_nPx1.bin Rs_VobcsN_2008_7DAY_GOM_Jan01_nPy1_c1.bin)
   ln -sf ../inputs/$FILE . || exit 1
end


echo 'listing now that the table has been set ...'
ls -l

# Loop through each state
set state_copy = 1
set ensemble_member_line = 1
set input_file_line = 2
set output_file_line = 3
while($state_copy <= $num_states)
   
   set ensemble_member = `head -n $ensemble_member_line ../$control_file | tail -n 1`
   set input_file      = `head -n $input_file_line      ../$control_file | tail -n 1`
   set output_file     = `head -n $output_file_line     ../$control_file | tail -n 1`

   #----------------------------------------------------------------------
   # Block 2: Convert the DART output file to form needed by ocean model
   #----------------------------------------------------------------------
   #
   # trans_sv_pv  creates the following files:
   #
   # data.cal.DART  ... containing the appropriate startdate_1, startdate_2
   # data.DART     ... containing the appropriate endTime, dumpFreq, and taveFreq
   #   S.YYYYMMDD.HHMMSS.[data,meta],
   #   T.YYYYMMDD.HHMMSS.[data,meta],
   #   U.YYYYMMDD.HHMMSS.[data,meta],
   #   V.YYYYMMDD.HHMMSS.[data,meta],
   # Eta.YYYYMMDD.HHMMSS.[data,meta], 

   ln -sfv ../$input_file assim_model_state_ic || exit 2

   ../trans_sv_pv || exit 2

   cp -pv data.cal.DART data.cal
   cp -pv data.DART     data

   # The DART output files must be renamed to match the data&PARM05 namelist.
   # This is pretty gory, but it works.

   set FNAME = `grep -i hydrogSaltFile data | sed -e "s#=##"`
   set FNAME = `echo  $FNAME | sed -e "s#hydrogSaltFile##"`
   set FNAME = `echo  $FNAME | sed -e "s#,##g"`
   set FNAME = `echo  $FNAME | sed -e "s#'##g"`
   mv -v S.*.*.data $FNAME  || exit 2

   set FNAME = `grep -i hydrogThetaFile data | sed -e "s#=##"`
   set FNAME = `echo  $FNAME | sed -e "s#hydrogThetaFile##"`
   set FNAME = `echo  $FNAME | sed -e "s#,##g"`
   set FNAME = `echo  $FNAME | sed -e "s#'##g"`
   mv -v T.*.*.data $FNAME  || exit 2

   set FNAME = `grep -i uVelInitFile data | sed -e "s#=##"`
   set FNAME = `echo  $FNAME | sed -e "s#uVelInitFile##"`
   set FNAME = `echo  $FNAME | sed -e "s#,##g"`
   set FNAME = `echo  $FNAME | sed -e "s#'##g"`
   mv -v U.*.*.data $FNAME  || exit 2

   set FNAME = `grep -i vVelInitFile data | sed -e "s#=##"`
   set FNAME = `echo  $FNAME | sed -e "s#vVelInitFile##"`
   set FNAME = `echo  $FNAME | sed -e "s#,##g"`
   set FNAME = `echo  $FNAME | sed -e "s#'##g"`
   mv -v V.*.*.data $FNAME  || exit 2

   set FNAME = `grep -i pSurfInitFile data | sed -e "s#=##"`
   set FNAME = `echo  $FNAME | sed -e "s#pSurfInitFile##"`
   set FNAME = `echo  $FNAME | sed -e "s#,##g"`
   set FNAME = `echo  $FNAME | sed -e "s#'##g"`
   mv -v Eta.*.*.data $FNAME  || exit 2

   # The DART output files must be renamed to match the data&PARM05 namelist.
   # This is succinct, but it does not work.

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

   #----------------------------------------------------------------------
   # Block 3: Run the ocean model
   #----------------------------------------------------------------------

   # Must determine if we are running in a queueing environment or not
   # so we know the form of the advance command.

   if ($?LS_SUBCWD) then

      mpirun.lsf ../mitgcmuv || exit 3
   
   else if ($?PBS_O_WORKDIR) then

      mpirun     ../mitgcmuv || exit 3

   else

      # This is a (temporary) violation of the DART standard of separating
      # the architecture-specific quantities and the model-specific things.
      # At some point in the future, the MPIRUN variable should not be hardwired
      # to an architecture-specific value.

      if ( -e ../nodelist ) then

         setenv NUM_PROCS `cat ../nodelist | wc -l`

# compas
#        set MPIRUN = /opt/mpich/myrinet/pgi/bin/mpirun
# atlas
         set MPIRUN = /share/apps/mpich1/pgi/bin/mpirun

         $MPIRUN -np $NUM_PROCS -nolocal -machinefile ../nodelist ../mitgcmuv || exit 3

      else
         echo "ERROR - there is no CENTRALDIR/nodelist for this execution."
         echo "ERROR - there is no CENTRALDIR/nodelist for this execution."
         echo "        The current working directory is: "`pwd`
         echo "        The contents of the directory are: "
         ls -l
         exit 3
      endif

   endif

   #----------------------------------------------------------------------
   # Block 4: Convert the ocean model output to form needed by DART
   #----------------------------------------------------------------------

   # trans_pv_sv reads the snapshot file after the model advance and must
   # convert them to an array that is prefaced by the time at which the 
   # model state is valid. This requires converting the timestep info
   # contained in the snapshot file into a real time. This requires
   # reading the 'data.cal' that _started_ the model advance, and the 
   # timestep increment from the 'data' namelist.

   # This code block makes the STRONG ASSUMPTION that there are precisely
   # two sets of snapshot files - the unavoidable one at time zero - and
   # the set for the END of the model advance. Now that trans_sv_pv 
   # sets data:endTime and dumpFreq, this should not be a problem. 

   # Remove the snapshot file at time zero.
   # We are interested in the snapshot file at the end of the advance.
   # Daily advances ... with a timestep of 900 seconds  86400/900 = 96
   rm *.0000000000.*

   # Extract the timestep from the ocean model output files.
   set TIMESTEP = `ls -1 S.*.data`
   set TIMESTEP = $TIMESTEP:r
   set TIMESTEP = $TIMESTEP:e

   echo $TIMESTEP | ../trans_pv_sv || exit 4

   # Move the updated state vector back to 'centraldir'
   mv assim_model_state_ud ../$output_file || exit 4

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3
end

echo "old data.cal is"
cat ../data.cal
echo ""
echo "new data.cal is"
cat data.cal.DART
echo ""
 
cp -pv data.cal.DART ../data.cal

# Change back to original directory and get rid of temporary directory
cd ..
# \rm -rf $temp_dir

# Remove the filter_control file to signal completion
# Is there a need for any sleeps to avoid trouble on completing moves here?
\rm -rf $control_file

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

