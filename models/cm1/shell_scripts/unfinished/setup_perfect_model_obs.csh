#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#-----------------------------------------------------------------------------
#
# This script is meant to provide an example of what is necessary to set up
# a very simple experiment that advances a single instance of CM1 (i.e. the 'truth')
# and harvests a set of synthetic observations from the truth.
#
# ASSUMPTIONS/REQUIREMENTS
#
# *) CM1 is required to be compiled *without* MPI so we can run in async==2
#
# *) CM1 is required to use netCDF restart files.
#
# *) A single CM1 model state will be available as well as a 
#    namelist.input, LANDUSE.TBL, etc.
#
# *) The DART obs_seq.in file will totally define the duration
#    of the truth run.
#
# *) The DART input.nml file has some required values as defined below.
#
# *) Each time CM1 is advanced, it will start from the same filename,
#    and the restart number in that filename will be 000001 - ALWAYS.
#    That filename will be a link to the most current model state.
#
# *) Each time CM1 has finished advancing the most current restart file
#    will be renamed to have the valid date/time of the model state.
#
#-----------------------------------------------------------------------------
# These are the DART namelist variables that must have these values.
# &io_filenames_nml
#   rpointer         = .true.
#   rpointer_file    = 'input_filelist.txt'
# /
# &perfect_model_obs_nml
#   adv_ens_command          = "./advance_model.csh"
#   async                    = 2
#   direct_netcdf_read       = .true.
#   direct_netcdf_write      = .true.
#   init_time_days           = -1
#   init_time_seconds        = -1
#   obs_seq_in_file_name     = "obs_seq.in"
#   start_from_restart       = .true.
# /
#-----------------------------------------------------------------------------

setenv JOBNAME     cm1_perfect
setenv JOBID       $$

# some systems don't like the -v option to any of the following

setenv COPY   'cp -pv'
setenv LINK   'ln -sv'
setenv REMOVE 'rm -rf'

#-----------------------------------------------------------------------------
# Define the directory that will hold the experiment (i.e. CENTRALDIR).
# All the necessary resources get staged here.
# All the run-time input should be defined/modified here.
#-----------------------------------------------------------------------------

set  CENTRALDIR = ${home}/temp/${JOBNAME}/job_${JOBID}

echo "${JOBNAME} ($JOBID) CENTRALDIR == $CENTRALDIR"

mkdir -p ${CENTRALDIR}
cd ${CENTRALDIR}

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
# DARTDIR      The location of the DART cm1 model directory
# CM1DIR       The location of the cm1 executable
# TRUTHDIR     The location of the initial model true state
#-----------------------------------------------------------------------------

set     DARTDIR = ${home}/svn/DART/cm1/models/cm1
set      CM1DIR = ${home}/svn/DART/cm1/models/cm1/src/cm1r18.3/run
set    TRUTHDIR = ${home}/temp/CM1R18.3

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
# Get the cm1 executable, control files, and data files.
#-----------------------------------------------------------------------------

${COPY} ${DARTDIR}/work/perfect_model_obs            . || exit 1
${COPY} ${DARTDIR}/work/advance_time                 . || exit 1
${COPY} ${DARTDIR}/work/input.nml  		     . || exit 1
${COPY} ${DARTDIR}/shell_scripts/advance_model.csh   . || exit 1
${COPY} ${DARTDIR}/work/obs_seq.in                   . || exit 1

${COPY} ${CM1DIR}/cm1.exe             	             . || exit 1

# perfect_model will read the contents of input_filelist.txt for the names
# of the CM1 restart file. After each model advance this same
# filename will be queried for the new model state.
# So, the actual file associated with this name must be updated.
# advance_model.csh is used by both perfect_model_obs and filter.
# Each model advance is performed in a unique model directory that
# has an instance or ensemble member number in the directory name.
# So, even though there is only 1 instance for a perfect_model_obs
# experiment, there must still be a directory name compatible with
# the advance_model logic.

set filename = "cm1out_rst_000001.nc"
set  dirname = dir_model_001

mkdir $dirname

${COPY} ${TRUTHDIR}/cm1out_rst_000001.nc    $dirname/$filename  || exit 2
${COPY} ${TRUTHDIR}/input_sounding          $dirname            || exit 2
${COPY} ${TRUTHDIR}/LANDUSE.TBL             $dirname            || exit 2
${COPY} ${TRUTHDIR}/namelist.input          $dirname/namelist.input.template || exit 2

# Enforce certain assumptions.
# 1) we are always starting from a netCDF file ... irst == 1
# 2) restart always be named cm1_our_rst_000001.nc ... rstnum == 1
# 3) the run_time value will be a dummy string to facilitate
#    the use of 'sed' in the advance_model.csh script.

ex $dirname/namelist.input.template <<ex_end
g;irst     ;s;= .*;= 1;
g;rstnum   ;s;= .*;= 1;
g;run_time ;s;= .*;= CM1_FORECAST_LENGTH;
g;rstfrq   ;s;= .*;= CM1_FORECAST_LENGTH;
wq
ex_end

grep CM1_FORECAST_LENGTH $dirname/namelist.input.template
if ( $status != 0 ) then
echo "The CM1 namelist file 'namelist.input' did not get updated correctly."
   echo "Check the syntax and try again. Aborting."
   exit 2
endif

${REMOVE} input_filelist.txt

echo "$dirname/$filename" > input_filelist.txt   || exit 2

# When perfect_model_obs starts, it needs to read the grid information.
# The DART/CM1 interface specifies the file for the grid information by
# the input.nml:&model_nml:cm1_template_file  variable, which defaults
# to cm1out_rst_000001.nc
#
# Consequently, we can get simply link to the true model state.

${LINK} $dirname/$filename                  cm1out_rst_000001.nc || exit 3
${COPY} dir_model_001/namelist.input.template     namelist.input || exit 3

echo ""
echo "cd ${CENTRALDIR}"
echo "Make sure that input.nml contents are correct."
echo "Make sure the obs_seq.in and the model state are appropriate."
echo "Make sure that namelist.input contents are correct."
echo "Make sure that input_filelist.txt contains the right file and the file exists."
echo "Launch ./perfect_model_obs"
echo ""

exit 0


