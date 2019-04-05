#!/bin/tcsh

#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -W 1:30
#BSUB -q casper
#BSUB -P P86850054
#BSUB -N
#BSUB -u raeder@ucar.edu
#BSUB -o obs_diags.%J
#BSUB -e obs_diags.%J
#BSUB -J obs_diags
#-----------------------------------------
#PBS  -N obs_diags
#PBS  -A P86850054
#PBS  -q casper
# Resources I want:
#    select=#nodes
#    ncpus=#CPUs/node
#    mpiprocs=#MPI_tasks/node
#PBS  -l select=1:ncpus=1:mpiprocs=1
#PBS  -l walltime=02:00:00
# Send email after a(bort) or e(nd)
#PBS  -m ae
#PBS  -M raeder@ucar.edu
# Send standard output and error to this file.
# It's helpful to use the $casename here.
#PBS  -o obs_diags.eo
#PBS  -j oe 
#--------------------------------------------

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# DART $Id$

# This script feeds a series of obs_seq.final files to obs_diag,
# which processes them according to its namelist.
# The resulting obs_diag_output.nc is fed to matlab to generate observation space plots
# according to the script which is defined in this script.
# It can be run as a batch job, or called by a batch job,
# because it can take longer than is allowed on many systems' login nodes. 
# It can be run interactively for smaller jobs.

if ($#argv == 0) then
   # User submitted, independent batch job (not run by another batch job).
   # CASE could be replaced by setup_*, as is done for DART_config.
   set arch_dir = /glade/scratch/${USER}/CASE/archive/esp/hist
   set date_string = '2017-03'
   set output_dir = Diags_NTrS_2017.3.1-2017.4.1H0_s0
#   Diags_NTrS_2017.3.1-2017.4.1H0_s0
#         ^    ^        ^          ^
#         |    |        |          Skip first N cycles for vert profiles
#         |    |        |          (or other details)
#         |    Start,   end dates
#         Domains

else if ($#argv == 1) then
   # Request for help; any argument will do.
   echo "Usage: call by user or script:"
   echo "   obs_diags.csh arch_dir 'date_string' output_dir [dart_dir]"
   echo "      arch_dir    = directory where $case.dart.e.cam_obs_seq_final.$date.nc are"
   echo "      date_string = cshell regular expression which selects the files, e.g. '2017-03-[01]',"
   echo "                    or set manually in the script."
   echo "      output_dir  = directory where output from obs_diag and matlab will be put."
   echo "                    I give it a name that will distinguish it from other obs space output"
   echo "                    in the same arch_dir."
   echo "      dart_dir    = optional local DART code root."

else
   # Script run by another (batch) script or interactively.
   set arch_dir    = $1
   set date_string = "$2"
   set output_dir  = $3
   # DART could be filled by setup_*, as is done for DART_config.
   if ($#argv == 4) 
      set DART = $4
   else
      set DART = ~/DART/reanalysis
   endif

#       echo 'Two scripts are required to generate observation space diagnostics '  
#       echo 'from obs_seq.[dates].final files; this one and ~/Scripts/matlab_cesm.csh.'  
#       echo '1) Go to the $EXEDIR directory, which contains the run directory.'  
#       echo '2) If needed, make a directory, e.g. "Obs_seq_final", to hold the obs_seq.final files.'  
#       echo '3) Link (or copy: >hsi  >cd hsi_dir  >get obs*)   the obs_seq.final files into Obs_seq_final.'  
#       echo '4) In EXEDIR run this script, giving it the directory where'  
#       echo '   the diagnostics will be created and Obs_seq_final (no ../):'  
#       echo '    > > >  ~/Scripts/diags_cesm.csh Diag_[details]_[dates] Obs_seq_final  < < < '
#       echo '5) Go to the diag_[dates] directory. '  
#       echo '6) Execute ~/Scripts/matlab_cesm.csh'  
#       exit
#    endif
endif

if (! -d $arch_dir:t) then
   echo "Missing arch_dir"
   exit 10
endif
cd $arch_dir

if (-d $output_dir) then
   echo "$output_dir exists; choose another name"
   exit 20
endif
mkdir $output_dir

# ls -1 does not work; unusable formatting.
# Separate the obs_seq.final files we want to process.
set files = `ls *cam_obs_seq_final.*${date_string}*`
if ($#files == 0) then
   echo "No files matching $date_string"
   exit 30
endif

if (! -d Obs_seqs) mkdir Obs_seqs
cd Obs_seqs
if ($files[1]:t == 'gz') then
   foreach f ($files)
      gunzip --stdout ../$f > $f:r
   end
   set compr = 'true'  
else
   foreach f ($files)
      ln -s ../$f $f:r
   endif
   set compr = 'false'  
endif
ls *obs_seq_final* >! ../${output_dir}/obs.list 

# Run obs_diag in the directory where we want the obs space output to be.
cd ../$output_dir

>>> Update env_archive.xml to copy input.nml to $archive/esp/hist

set yr = ()
set mo = ()
set dy = ()
set t = 1
while ($t <= 2) 
   set this = 1
   if ($t == 2) set this = $#files

   set date =  $files[$this]:r
   if ($comp == 'true') set date = $date:r

   set parts = (`echo $date:t | sed -e 's/-[0]/ /g'`)
   set yr = ($yr $parts[1])
   set mo = ($mo $parts[2])
   set dy = ($dy $parts[3])

   @ t++
end
echo "Begin and end dates = $yr[1]-$mo[1]-$dy[1], $yr[2]-$mo[2]-$dy[2]"
ex ../input.nml<< ex_end
/obs_diag_nml/
/obs_sequence_name/
s;= '.*';= "";
/obs_sequence_list/
s;= '.*';= "./obs.list";
/first_bin_center/
s;= .*;= $yr[1], $mo[1], $dy[1], 0, 0, 0 ,
/last_bin_center/
s;= .*;= $yr[2], $mo[2], $dy[2], 0, 0, 0 ,
wq
ex_end

# /time_to_skip/
# s;= '.*';=  $yr_skip, $mo_skip, $day_skip, $hr_skip, 0, 0 ,

echo "Running ${DART}/models/cam-fv/work/obs_diag"
${DART}/models/cam-fv/work/obs_diag >&! obs_diag.out 

if ($status != 0) exit 50

#-----------------------------------------
# Plot selected data from obs_seq_output.nc using matlab.

# Plot all of the obs_types listed in input.nml:{assimilate,evaluate}_these_obs_types.
# Line number of the assimilate_these_obs_types variable.
set line_as = (`grep -n assimilate_these_obs_types ../input.nml`)
set line_assim = `echo $line_as[1] | sed -e "s/://"`

# Line number of the evaluate_these_obs_types variable.
set line_ev = (`grep -n evaluate_these_obs_types ../input.nml`)
if ($status != 0) then
   set line_eval = 1000000
else
   set line_eval = `echo $line_ev[1] | sed -e "s/://"`
endif

# Line number of the end of the obs_kind_nml
# -m 1: stop looking after finding 1 '/'
set line = (`tail -n +$line_assim ../input.nml | grep -n -m 1 '/'`)
set nml_len = `echo $line[1] | sed -e "s/://"`
@ line_end = ($line_assim - 1) + $nml_len

set obsnames = ()
# It's 4 because of the extra word added by grep to list the line number.
if ($#line_as == 4) then
   set obsnames = (`echo $line_as[4] | sed -e "s/,//"`)
endif

# Figure where assimilate_these_obs_types ends.
# The first line is handled separately, so no '+1' on the end.
if ($line_eval < $line_end) then
   @ last = $line_eval - 1
   @ lines = $last - $line_assim
   # Extract the (rest of) the values in assimilate_these_obs_types
   set obsnames = ($obsnames `head -n $last ../input.nml | tail -n $lines | sed -e "s/,//"`)

   # Include the evaluate_these_obs_types values too
   if ($#line_ev == 4) then
      set obsnames = ($obsnames `echo $line_ev[4] | sed -e "s/,//"`)
   endif
   @ last = $line_end - 1
   @ lines = $last - $line_eval
   set obsnames = ($obsnames `head -n $last ../input.nml | tail -n $lines | sed -e "s/,//"`)

else
   # No evaluate_ variable; use the ending / as the end of the values to harvest.
   @ last = $line_end - 1
   @ lines = $last - $line_assim
   # Extract the the values in assimilate_these_obs_types
   set obsnames = ($obsnames `head -n $last ../input.nml | tail -n $lines | sed -e "s/,//"`)
endif


# Create the matlab script
echo "fname = 'obs_diag_output.nc';"   >! script.m
foreach obs ($obsnames)
foreach copy (totalspread bias)
foreach func (plot_rmse_xxx_evolution plot_rmse_xxx_profile)
   echo "$func(fname,$copy,'obsname',$obs)"         >> script.m
end
end
end

echo "exit" >> script.m


# Alternatively, hard-wire the specific output you'd like to see.
# cat << EndOfFile > script.m
# 
# fname   = 'obs_diag_output.nc';
# copy    = 'totalspread';
# obsname = 'RADIOSONDE_TEMPERATURE';
# plotdat = plot_rmse_xxx_evolution(fname, copy, 'obsname', obsname);
# obsname = 'RADIOSONDE_TEMPERATURE';
# plotdat = plot_rmse_xxx_evolution(fname, copy);
# 
# copy    = 'bias';
# plot_profile(fname, copy, 'obsname', 'RADIOSONDE_TEMPERATURE');
# plot_profile(fname, copy, 'obsname', 'RADIOSONDE_SPECIFIC_HUMIDITY');
# 
# plot_evolution(fname, copy, 'obsname', obsname, 'level', 4,  'verbose', 'no');
# 
# plot_rmse_xxx_profile(fname, copy, 'obsname', obsname);
# obsname = 'RADIOSONDE_SPECIFIC_HUMIDITY';
# plot_rmse_xxx_profile(fname, copy, 'obsname', obsname);
# 
# copy    = 'rmse';
# obsname = 'RADIOSONDE_TEMPERATURE';
# plot_bias_xxx_profile(fname, copy, 'obsname', obsname);
# 
# obsname = 'RADIOSONDE_SPECIFIC_HUMIDITY';
# plot_bias_xxx_profile(fname, copy, 'obsname', obsname);
# 
# exit
# EndOfFile
# 
#---------------------------------------------------------------------
# Note that Matlab requires the last line of the script to 
# be 'exit' or 'quit' and that the file extension be LEFT OFF the 
# command line. Strange. But true.
#---------------------------------------------------------------------

matlab -r script

exit
