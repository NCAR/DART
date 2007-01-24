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

# Shell script to do repeated segment integrations
# BUT, don't know how to do MATLAB from script yet

set exp = T85_noPSQ

# for January experiment
set output_dir = 01
set chunk1 = 1
set chunkn = 2

set save_freq = 2
set mod_save = 1

set input = input_

if (-d ${exp}) then
   echo ${exp} already exists
else
   mkdir ${exp}
   cp namelistin ${exp}/namelistin
endif

# Have an overall outer loop
set i = $chunk1
while($i <= $chunkn)
   echo ' '
   echo ' '
   echo starting iteration $i
   echo ' '
   echo ' '

# remove clminput files from previous experiment
   if (-e clminput_1.nc && $i == $chunk1) rm clminput_[1-9]*.nc
   if (-e caminput_1.nc && $i == $chunk1) rm caminput_[1-9]*.nc

   cp $input$i.nml input.nml

   if ($i == 1) then
#     Fill standard names with data for the beginning of an experiment
#     filter_ics may be filled in input.nml, instead.
      cp T85.clm3.0.i.2003-01-01-all-fields.nc clminput.nc
#      cp caminput_T85_1.nc caminput.nc
#      cp filter_T42_1-1-03_ics filter_ics
#      cp perfect_31_ics perfect_ics
   else if ($i == $chunk1) then
#     get restart data from experiment archive
      @ j = $i - 1
      cp ${exp}/${output_dir}_$j/filter_restart filter_ics
      cp ${exp}/${output_dir}_$j/CLM/clminput*.nc .
      cp ${exp}/${output_dir}_$j/CAM/caminput*.nc .
   else
      mv filter_restart filter_ics
      mv perfect_restart perfect_ics
   endif

# Run perfect_model_obs with async=2 in input.nml
#   ./perfect_model_obs >&! perfect_$i.out

# Run filter
   ./filter >&! dart_$i.log

# Move the netcdf files to an output directory
   mkdir ${exp}/${output_dir}_$i
# cleanupgrade; see above.  
   mkdir ${exp}/${output_dir}_$i/CLM
   mkdir ${exp}/${output_dir}_$i/CAM
   mv Prior_Diag.nc Posterior_Diag.nc    ${exp}/${output_dir}_$i
   mv obs_seq.final                      ${exp}/${output_dir}_$i
#   mv data_cam_prob                     ${exp}/${output_dir}_$i
   cp input.nml                          ${exp}/${output_dir}_$i
   mv cam_out_temp1                      ${exp}/${output_dir}_$i

# Copy the filter restart to start files
   cp filter_restart ${exp}/${output_dir}_$i
   cp clminput_[1-9]*.nc ${exp}/${output_dir}_$i/CLM
   cp caminput_[1-9]*.nc ${exp}/${output_dir}_$i/CAM
   if ($i % $save_freq != $mod_save) then
      @ j = $i - 1
      if (-e ${exp}/${output_dir}_$i/CLM/clminput_20.nc) then
         cd ${exp}/${output_dir}_$j
         gzip -r CLM
         tar cf clminput.gz.tar CLM
         rm -rf CLM
         cd ../..
      else
         echo 'NO clminput.20;  ABORTING' >> dump
         stop
      endif       
      if (-e ${exp}/${output_dir}_$i/CAM/caminput_20.nc) then
         cd ${exp}/${output_dir}_$j
         gzip -r CAM
         tar cf caminput.gz.tar CAM
         rm -rf CAM
         cd ../..
      else
         echo 'NO caminput.20;  ABORTING' >> dump
         stop
      endif       
   endif

# Move along to next iteration
   echo ' '
   echo ending iteration $i
   echo ' '
   @ i++

end

mv input.nml    ${exp}
mv dart_log.out $exp
mv namelist $exp
