#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# DART $Id$

#==========================================================================

# Script to process the obs_seq_final files from multiple cases.
# and compare them in obs space.  The stages are:
# 1) Run each set through obs_sequence_tool to make the sets 
#    as consistent as they need to be.  E.g. remove obs types from one set 
#    which aren't in the other.
# 2) Run all the sets through obs_common_subset
# 3) Process each common obs_seq_final set through obs_diag
# 4) Feed the resulting obs_diag_output.nc files to 
#    matlab:two_experiments_overview{_Breck16}.m

# >>> INSTRUCTIONS:

# >>> Run st_archive on all cases before running this script.          <<<
# >>> Create a scratch directory separate from the cases.              <<<
# >>> Copy input.nml into it and edit:                                 <<<
# >>>    obs_sequence_tool_nml; e.g. to remove obs types               <<<
# >>>    obs_common_subset_nml; to create obs_seq_finals containing    <<<
# >>>                           only obs used by all cases.            <<<
# >>>    obs_diag_nml;          generate the obs_diag_output.nc        <<<
# >>> Edit diags_batch.csh to match obs_common_subset_nml              <<<
#       This may be all done by the arguments this script passes.
# >>> Submit this job from casper from the new scratch directory.      <<<
# >>>    (or turn off the matlab:two_experiments_overview call)        <<<

#-----------------------------------------
# Submitting the job to casper (or other NCAR DAV clusters?) requires using slurm.

# Important things to know about slurm:
#
# sinfo     information about the whole slurm system
# squeue    information about running jobs
# sbatch    submitting a job
# scancel   killing a job
#
#==========================================================================
#
#SBATCH --job-name=obs_to_two_experiments
#SBATCH --ntasks=1 
#SBATCH --time=02:00:00
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=raeder@ucar.edu
#SBATCH --account=NCIS0006
#SBATCH --partition=dav
#SBATCH --ignore-pbs
# 
#-----------------------------------------
#PBS  -N obs_to_two_experiments
#PBS  -A NCIS0006
#PBS  -q share
# Resources I want:
#    select=#nodes
#    ncpus=#CPUs/node
#    mpiprocs=#MPI_tasks/node
#PBS  -l select=1:ncpus=1:mpiprocs=1
#  This can take a while; 
#     obs_sequence_tool on weeks of files
#     multiple obs_diags, 
#     matlab batch processing of obs_diag_output.nc files     
#PBS  -l walltime=02:00:00
# Send email after a(bort) or e(nd)
#PBS  -m ae
#PBS  -M raeder@ucar.edu
# Send standard output and error to this file.
# It's helpful to use the $CASE name here.
#PBS  -o obs_to_two_experiments.eo
#PBS  -j oe 
#--------------------------------------------

module list

set my_dir   = `pwd`
echo my_dir = $my_dir
set DART = ~/DART/reanalysis
set tool_dir = ${DART}/models/cam-fv/work_casper
set scr_dir  = ${DART}/models/cam-fv/shell_scripts/cesm2_1

#--------------------------------------------
# This first iteration of this script uses some existing infrastructure,
# Namely some file and directory naming, and an existing scripts.  
# This restricts the case names for now: ${case_base}.###  

# 001 = the 2017/Jan-June assim using NCEP+ACARS+GPS
# 004 = the 2017Jan 1-14  assim using NCEP+ACARS+GPS+AIRS (T only)
set case_base = f.e21.FHIST_BGC.f09_025.CAM6assim
set cases     = (001 004)

# We may want common obs from more than 2 cases,
# but not compare all of the cases.  Select the pairs here.
# What will be plotted are the rmse and bias of the first minus the second,
# so that positive means the 2nd (later, newer, ...) is better.
set compare   = (001-004 )

# A csh pattern that will find the days of desired obs_seq files.
# set obs_times_pattern = '2017-{01}-01'
# full set; 
set obs_times_pattern = '2017-{01}-{0*,1[0-4]}'

#--------------------------------------------
# Possibly install a section to set up input.nml (from a template)
# to be sure that it's consistent with what's expected in this script.

#--------------------------------------------
set matlab = true
if ($?PBS_O_WORKDIR) then
   cd $PBS_O_WORKDIR
   set matlab = false
else if ($?SLURM_SUBMIT_DIR) then
   cd $SLURM_SUBMIT_DIR
endif

#--------------------------------------------
foreach ext ($cases)
   set case = ${case_base}.$ext

   # Create the obs_seq file list.
   cd ../$case/archive/esp/hist
   ls -1 ${case}.dart.e.cam_obs_seq_final.${obs_times_pattern}* >! obs_${ext}.list
   cd $my_dir
   ln -sf ../$case/archive/esp/hist/obs_${ext}.list .

   # Modify the obs_seq_finals to be consistent with each other.
   # E.g. remove the AIRS obs from some.
   # Hopefully, if it finds no AIRS, it won't die.

   # WARNING: obs_common_subset was designed to use cases which used the same
   #          obs_seq_output file set.  There's no guarantee that 
   #          obs_sequence_tool can make 2 different sets have the same obs 
   #          in the same order.
   #          NCEP+ACARS+GPS and NCEP+ACARS+GPS+AIRS worked,
   #          probably because the latter is just the former with AIRS added in.

   $scr_dir/obs_seq_tool_series.csh $case_base $ext
   if ($status != 0) exit 10

end

# Use the first case's list of files to harvest dates 
# for the diagnostics directories used by diags_batch, below. 
set ofile = `head -n 1 obs_$cases[1].list `
set first = `echo $ofile:e | sed -e "s#-# #g"`
set ofile = `tail -n 1 obs_$cases[1].list `
set last  = `echo $ofile:e | sed -e "s#-# #g"`
set date_span = $first[1].$first[2].$first[3]-$last[1].$last[2].$last[3]


#--------------------------------------------
# Create the 'common obs' obs_seq files.
# They will appear in the same directories as the input obs_seq files,
# but with extensions provided by obs_common_subset_nml.

$tool_dir/obs_common_subset
if ($status != 0) then
   echo obs_common_subset failed with status = $status
   exit 30
endif

set suffix = `grep filename_out_suffix input.nml | cut -d"=" -f 2 | sed -e "s#'##g"`

#--------------------------------------------
# Create an obs_diag_output.nc from the common obs files of each case.

foreach ext ($cases)
   # For now this naming is hard-wired into obs_seq_tool_series.csh,
   # so it's used here.

   # Cd into directory and add ../ to the ls path 
   # so that the resulting filenames have ../ at the beginning,
   # which diags_batch.csh needs.
   cd Obs_${ext}_noAIRS 
   ls -1 ../Obs_${ext}_noAIRS/*$suffix >! ../Obs_${ext}_noAIRS.list
   cd ..

   set diagnostics_dir = Diags.${ext}_NTrS_${date_span}_s0
   set case = ${case_base}.$ext

   $scr_dir/diags_batch.csh ${case} Obs_${ext}_noAIRS $diagnostics_dir

end

#--------------------------------------------
# Compare the 2 cases

if ($matlab == true) then
foreach pair ($compare)
   set first = `echo $pair | cut -d'-' -f 1`
   set last  = `echo $pair | cut -d'-' -f 2`

   echo "addpath('$DART/diagnostics/matlab','-BEGIN')"                                      >! script.m
   echo "file1 = './Diags.${first}_NTrS_${date_span}_s0/obs_diag_output.nc'"                  >> script.m
   echo "file2 = './Diags.${last}_NTrS_${date_span}_s0/obs_diag_output.nc'"                   >> script.m
   echo "two_experiments_overview_Breck2016(file1,'$first',file2,'$last','FlagLevel',.05) " >> script.m

   matlab -r script

   tar -c -f ${case_base}_${pair}_comparison.tar \
       rmse_diff.pdf abs_bias_diff.pdf num_obs_used.pdf num_obs_used_diff.pdf
   rm  rmse_diff.pdf abs_bias_diff.pdf num_obs_used.pdf num_obs_used_diff.pdf

end
endif
#--------------------------------------------

#==========================================================================
# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
