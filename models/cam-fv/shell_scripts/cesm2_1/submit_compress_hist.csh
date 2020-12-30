UNTESTED

#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# ------------------------------------------------------------------------------

#SBATCH --job-name=compress_hist
#SBATCH -o %x_%j.eo 
#SBATCH -e %x_%j.eo 
# 80 members
# Each type is done as a separate cmdfile
#SBATCH --ntasks=80 
#SBATCH --time=01:00:00
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=raeder@ucar.edu
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=dav
#SBATCH --ignore-pbs
# 
#-----------------------------------------
#PBS  -N compress_hist.csh
#PBS  -A YOUR_ACCOUNT
#PBS  -q regular
# #PBS  -q share
#PBS  -l select=5:ncpus=36:mpiprocs=36
# #PBS  -l select=1:ncpus=1:mpiprocs=1
#PBS  -l walltime=00:10:00
#PBS  -o compress_hist.out
#PBS  -j oe 
#PBS  -k eod 

# Get CASE environment variables from the central variables file.
# This should make environment variables available in compress_hist.csh too.
source YOUR_CASEROOT/data_scripts.csh

# 'sets' performs better when ordered by decreasing size (clm2 cpl cam cice hist dart)
# but this can handle only 1 entry in $sets.
# -k means keep the input files after creation of the output file.
set comp_cmd      = 'gzip -k'
# No -k option on casper.
# set comp_cmd      = 'gzip'
set ymds          = 2010-07-17-64800
# set ymds          = 2012

# set data_dir      = ${pr}/${data_CASE}/cpl/hist
set sets          = (cpl)
set types         = ( ha2x1d hr2x ha2x3h ha2x1h ha2x1hi )

# set data_dir      = /glade/p/nsc/ncis0006/Reanalyses/${data_CASE}/rof/hist
# set sets          = (mosart)
# set types         = ( h0 )

# set data_dir      = /glade/scratch/${USER}/${data_CASE}/run
# set sets          = (hist dart)
# set sets          = (clm2 cpl cam cice)

set stages        = (none)

if ($?PBS_O_WORKDIR) then
   cd $PBS_O_WORKDIR
else if ($?SLURM_SUBMIT_DIR) then
   cd $SLURM_SUBMIT_DIR
endif

${data_CASEROOT}/compress_hist.csh "$comp_cmd" $ymds "$sets" "$types" "$stages"
