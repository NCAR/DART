#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# Shell script to do repeated segment integrations
# BUT, don't know how to do MATLAB from script yet

set exp = Ex28
set output_dir = 01
set chunk1 = 1
set chunkn = 7

set save_freq = 7
set mod_save = 1
# save 1/7 ... 1/21    by 7  mod_save 1      1...27
#      1/27... 1/31    by 2  mod_save 0      28...31
#      2/2 ... 2/10    by 2  mod_save 1      1...11
#      2/13... 2/25    by 2  mod_save 0      12...25  (skipped a day!)
#   if ($i % $save_freq != $mod_save) then

# Exp  grps #/grp  conf.slope  cut-off  #days
#   1    1   20      0.8        0.2      9
#   2    4   20      0.0        0.2      7
#   3    4   20      0.8        0.2      7
#  (4    1   20      0.0        0.2      7   postponed , but use in SingleBatch)
#   5    1   80      0.8        0.2      7
#   8    8   80      0.8        0.2      4
#   9    4   80      0.8        0.2     31 + 26
#    00 ; creates a restart valid at 00Z on 1/29, rather than 12Z
#  10    4   80      0.8        0.2      7
#    first clm fix; use Sep(July?) clm file for Jan 1, but keep an updated clminput
#          for each ens member.
#  11    4   80      0.8        0.2      7
#    use CamOnly/Archive/T42.clm2.i.2003-01-01-00000.nc as first land file, 
#    and use updated clminputs
#  12    4   80      0.8        0.2      7
#    use data between lowest CAM level and 50 hPa (instead of 200 hPa in previous exp)
#  13    4   80      0.8        0.2      7
#    like 11, but using 4x/day data (>200 hPa).
#  14    4   80      0.8        0.2      7
#    like 13, but using modified assim_tools_mod;update_from_obs_inc from Jeff
#             AND I mistakenly used the filter that includes data up to 50 hPa (like 12)
#  16    4   80      0.8        0.2      7
#    like 13, but using modified assim_tools_mod;update_from_obs_inc from Jeff
#  17    4   80      0.8        0.2      6
#    like 13, but using data up to 50 hPa (standard assim_tools)
#    network problems may have eliminated restart ability at day 6

# Exp  grps #/grp  conf.slope  cut-off  cov_infl  #days  (SSTs still climatlgcl)
#  18    4   80      0.0        0.2       1.01      6
#  19    4   80      0.0        0.2       1.03      2
#  20    4   80      0.0        0.2       1.4       7
#  21    4   80      0.0        0.2       2.25      2
#  22    4   80      0.0        0.2       1.8       7
#  23    4   80      0.0        0.2       1.6       7
#
# Resume branch from 13
#  24    4   80      0.2        0.2       1.0       7     like 16, but conf_slp
#  25    4   80      0.6        0.2       1.0       2     corrupted restarts
#                                                         but no apparent blowup
#  26    4   80      0.6        0.3       1.0       7     like 25 but cut-off
#  27    4   80      0.9        0.2       1.0       7     
#  28    4   80      0.7        0.2       1.0       7     
# 
-------------------------------------------------------------------------

set input = input_
# echo input = $input > input.log

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

   cp $input$i.nml input.nml
#    echo iteration $i >> input.log
#    ls -l input* >> input.log

   if ($i == 1) then
      cp filter_restart_31 filter_ics
#      echo filter_restart already copied
   else if ($i == $chunk1) then
      @ j = $i - 1
      cp ${exp}/${output_dir}_$j/filter_restart filter_ics
      cp ${exp}/${output_dir}_$j/CLM/clminput*.nc .
   else
      mv filter_restart filter_ics
   endif


# Run filter
   csh ./sync_submit.csh | ./filter

# Move the netcdf files to an output directory
   mkdir ${exp}/${output_dir}_$i
   mkdir ${exp}/${output_dir}_$i/CLM
   mv Prior_Diag.nc Posterior_Diag.nc                 ${exp}/${output_dir}_$i
   mv posterior_obs_diagnostics prior_obs_diagnostics ${exp}/${output_dir}_$i
#   mv data_cam_prob                                   ${exp}/${output_dir}_$i
   mv input.nml                                       ${exp}/${output_dir}_$i
   cp obs_seq_jan$i.out                               ${exp}/${output_dir}_$i/obs_seq.out

# Copy the filter restart to start files
   cp filter_restart ${exp}/${output_dir}_$i
   cp clminput_*.nc ${exp}/${output_dir}_$i/CLM
#   if ($i > 1) then
   if ($i % $save_freq != $mod_save) then
      @ j = $i - 1
      if (-e ${exp}/${output_dir}_$i/filter_restart)    rm ${exp}/${output_dir}_$j/filter_restart 
      if (-e ${exp}/${output_dir}_$i/CLM/clminput_80.nc) rm ${exp}/${output_dir}_$j/CLM/clminput*.nc 
   endif

# Move along to next iteration
   echo ' '
   echo ending iteration $i
   echo ' '
   @ i++

end
