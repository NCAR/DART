#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$
#
# Shell script to do repeated segment integrations
# This runs on the master node (c#, where # is the largest #
# allotted to this job).  submit_cam.csh runs CAM on 
# all the other allotted nodes.

# set echo verbose

set PBS_O_WORKDIR = $1
set PBS_NODEFILE = $2

set output_dir = $PBS_O_WORKDIR/Exp9_
set obs_seq = $PBS_O_WORKDIR/obs_seq_jan
set input = $PBS_O_WORKDIR/input_
set chunk1 = 1
set chunkn = 2
echo output_dir,nodefile = $output_dir $PBS_NODEFILE 

# cd to working directory for this master node
ls /scratch/local 
if (-d /scratch/local/dartcam) rm -rf /scratch/local/dartcam
mkdir /scratch/local/dartcam
cd /scratch/local/dartcam

cp $PBS_O_WORKDIR/caminput.nc .
cp $PBS_O_WORKDIR/clminput.nc .
cp $PBS_O_WORKDIR/final* .
# cp $PBS_O_WORKDIR/perfect_ics .
# if starting from restart; cp $PBS_O_WORKDIR/filter_ics .  BELOW
# otherwise
# cp perfect_ics filter_ics

pwd  
ls -l  

# Have an overall outer loop
set i = $chunk1
while($i <= $chunkn)

   echo ' '
   echo ' '
   echo outer loop i is $i  
# >> $PBS_O_WORKDIR/filter_cam.err

# kdr debug temp_ic
# copy namelist into generic name for cam nodes to use
   cp ${input}$i.nml $PBS_O_WORKDIR/input.nml
# copy input files for this chunk to filter node working directory
   cp $PBS_O_WORKDIR/input.nml .
   cp ${obs_seq}$i.out .

   if ($i == 1) then
      cp $PBS_O_WORKDIR/filter_restart_31 filter_ics
   else if ($i == $chunk1) then
      @ j = $i - 1
      cp $PBS_O_WORKDIR/filter_restart_$j filter_ics
   else
      mv filter_restart filter_ics
   endif

echo before submit_cam.csh  > $1/submit_cam.err
# Run perfect model obs to get truth series.
#    csh $PBS_O_WORKDIR/submit_cam.csh $PBS_O_WORKDIR $PBS_NODEFILE | \
#        $PBS_O_WORKDIR/perfect_model_obs
# echo finished perfect_model_obs   
# >> $PBS_O_WORKDIR/filter_cam.err

# Run filter
   csh $PBS_O_WORKDIR/submit_cam.csh $PBS_O_WORKDIR $PBS_NODEFILE | \
       $PBS_O_WORKDIR/filter 
echo finished filter   
# >> $PBS_O_WORKDIR/filter_cam.err

# Move the netcdf files to an output directory
   mkdir ${output_dir}$i
   mv Prior_Diag.nc Posterior_Diag.nc ${output_dir}$i
#   True_State.nc 
   mv prior_obs_diagnostics posterior_obs_diagnostics ${output_dir}$i
   mv data_cam_prob ${output_dir}$i
   mv input.nml ${output_dir}$i


# Copy the perfect model and filter restarts to start files
#   cp perfect_restart perfect_ics
   if ($i == $chunkn) then
      mv filter_restart $PBS_O_WORKDIR/filter_restart_$i
   else 
      cp filter_restart $PBS_O_WORKDIR/filter_restart_$i
   endif
   if ($i > 1) then
      @ j = $i - 1
      rm $PBS_O_WORKDIR/filter_restart_$j 
   endif

# Move along to next iteration
   @ i++

end

# remove working directory on master node
cd /scratch/local
rm -rf dartcam
