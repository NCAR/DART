#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
#==========================================================================
# this script submits the jobs to advance the ROMS ensemble members using
# a job array. Each ensemble member gets executed on its own node.
#
#SBATCH --array=1-ENSEMBLESIZE
#SBATCH --ntasks 16 
#SBATCH --time=10:00
#SXXXXX --dependency=afterok:FILTERJOBID
#
#==========================================================================
#
# Things to note: many strings are intended to be replaced when this
# template gets copied and ultimately submitted. Anything that starts
# with 'My' is a string that gets replaced in the "normal" ROMs fashion.
# The next few are not standard. 
#
# ENSEMBLESIZE           gets replaced in the 'stage_experiment.csh' script
# SXXXXX                 gets replaced if this is a dependent job (with SBATCH)
# FILTERJOBID            gets replaced by the job ID that must finish first
# EXPERIMENT_DIRECTORY   gets replaced in the 'stage_experiment.csh' script

# machine-specific dereferencing
if ($?LS_SUBCWD) then
   set LAUNCHCMD = "mpirun.lsf"
else if ($?SLURM_JOB_ID) then
   set LAUNCHCMD = "mpirun -np $SLURM_NTASKS -bind-to core"
else if ($?PBS_O_WORKDIR) then
   set LAUNCHCMD = "mpirun -np $SLURM_NTASKS -bind-to core"
else
   set LAUNCHCMD = "aprun -n 16"
endif

# set some useful environment variables
set   DSTART = MyDSTART
set instance = $SLURM_ARRAY_TASK_ID
echo "STARTING ENSEMBLE MEMBER $instance at "`date`

set INSTANCE_DIRECTORY = `printf "instance_%04d" $instance`

# get to work
cd EXPERIMENT_DIRECTORY/$INSTANCE_DIRECTORY

rm -f log_$instance.txt

echo "advancing instance $instance at ..."`date`

\cp ../s4dvar.in.template s4dvar.in
set OBS_PREF = ../Obs/obs
set NEW_OBS     = `printf %s_%d.nc ${OBS_PREF} $DSTART`
MySUBSTITUTE s4dvar.in MyOBSname   $NEW_OBS

${LAUNCHCMD} ../MyROMS_EXE MyROMS_STDIN >& log_$instance.txt

# Check for successful completion - log file should NOT have something like:
# Blowing-up: Saving latest model state into  RESTART file
grep -i blow log_$instance.txt > /dev/null
if ($status == 0) then
   echo "ROMS instance $instance FAILED."
   echo "ROMS instance $instance FAILED."
   echo "ROMS instance $instance FAILED."
   exit 1
endif

# sometimes we need the full name, sometimes we need it without the extension
set RST_FILE = MyRSTNAME
set DAI_FILE = MyDAINAME
set OBS_FILE = MyMODname
set RST_ROOT = $RST_FILE:r
set DAI_ROOT = $DAI_FILE:r
set OBS_ROOT = $OBS_FILE:r

# The ROMS restart file will be treated as the DART prior.
# Create a ROMS POSTERIOR file that will be updated by DART and
# tag the output with the model time.

set DSTART_STRING = `ncdump -v dstart ${DAI_FILE} | grep '^ dstart = '`
set DSTART = `echo $DSTART_STRING | sed -e "s#[=;a-z_ ]##g"`

set ROMS_PRIOR     = `printf %s_%04d_%d.nc ${RST_ROOT} $instance $DSTART`
set ROMS_POSTERIOR = `printf roms_posterior_%04d_%d.nc $instance $DSTART`
set ROMS_OBSFILE   = `printf %s_%04d_%d.nc ${OBS_ROOT} $instance $DSTART`
set SAFETY         = `printf roms_dai_original_%04d_%d.nc $instance $DSTART`

\cp -v ${DAI_FILE} ${SAFETY}          || exit 1
\mv -v ${RST_FILE} ${ROMS_PRIOR}      || exit 1
\mv -v ${DAI_FILE} ${ROMS_POSTERIOR}  || exit 1
\mv -v ${OBS_FILE} ${ROMS_OBSFILE}    || exit 1

echo
echo "#---------------------------------------------------------------------"
echo "# ROMS instance $instance completed at "`date`
echo "#---------------------------------------------------------------------"
echo

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

