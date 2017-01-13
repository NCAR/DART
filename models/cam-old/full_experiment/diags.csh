#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# DART source directory on this machine

set DART = ~/Pre-J/DART
set exp_dir = `pwd`
set innov_freq = 1

if ($#argv == 0) then
   echo 
   echo '1) Sets up environment for and runs the obs_diag program for'
   echo '   observation space diagnostics in the $exp directory.'
   echo '2) Runs the suite of matlab diagnostics.'
   echo '3) Generates some innovation files and optionally (securely)'
   echo 'copies the results to another machine.'
   echo 'Takes several arguments, but must be edited to change'
   echo 'other characteristics.'
   echo 
   echo 'Usage: diags.csh last_obs_seq ps_direct [exp_mach:/directory]'
   exit
else if ($#argv == 3) then
   set destin = $3
endif

set last_dir = $1
set ps_dir = $2

ln -s ../topog_file.nc topog_file.nc
ln -s ../caminput.nc caminput.nc
cp $last_dir/input.nml .
vi input.nml
$DART/models/cam/work/obs_diag >& diag.out

# matlab batch job
# addpath $DART/diagnostics/matlab
# to do from home; 
# matlab -nodisplay >&! matlab.out << EOF
# orig 
matlab >&! matlab.out << EOF
fit_ens_mean_time
fit_ens_spread_time
obs_num_time
obs_num_vertical
fit_mean_spread_time
fit_ens_mean_vertical
fit_ens_bias_vertical
exit
EOF

mkdir $ps_dir
mv *.m *.dat *.ps input.nml $ps_dir
tar c -f {$exp_dir:t}_${last_dir}_Diags.tar $ps_dir/*.ps $ps_dir/input.nml

set n = 1
set more = true
while ($more == true)
   set obs_seq = 01_0$n
   if ($n > 9) set obs_seq = 01_$n
   if ($obs_seq == $last_dir) set more = false

   if (! -e $obs_seq/Innov.nc && $n % $innov_freq == 0) then
      cd $obs_seq
      innov
      cd ..
#      tar r -f {$exp_dir:t}_${last_dir}_Diags.tar  $obs_seq/*.nc 
#      tar r -f {$exp_dir:t}_${last_dir}_Diags.tar  $obs_seq/*.nc $obs_seq/input.nml
   endif
# orig
   tar r -f {$exp_dir:t}_${last_dir}_Diags.tar  $obs_seq/P*.nc 
   @ n++
end

if ($#argv == 3) then
   scp {$exp_dir:t}_${last_dir}_Diags.tar $destin 
   scp ../explist $destin
endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

