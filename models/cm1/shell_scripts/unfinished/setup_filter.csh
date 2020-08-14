#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#-----------------------------------------------------------------------------
#
# This script is meant to provide an example of what is necessary to set up
# a very simple experiment that involves DART and a very small ensemble of
# CM1 model states. This is guaranteed to be scientifically bogus - it is just
# a test of the file motion and mechanics.
#
# ASSUMPTIONS/REQUIREMENTS
#
# *) CM1 is required to be compiled *without* MPI so we can run in async==2
#
# *) CM1 is required to use netCDF restart files.
#
# *) A collection of CM1 model states will be available. Each of those
#    model states is required to have its own namelist.input, LANDUSE.TBL,
#    etc. Each set of states and resources will be in a separate directory
#    named for each ensemble member ... i.e.  xxx/m1/* xxx/m2/* xxx/m3/* ...
#
# *) The DART obs_seq.out file will totally define the duration
#    of the assimilation.
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
# &filter_nml
#   adv_ens_command          = "./advance_model.csh"
#   async                    = 2
#   direct_netcdf_read       = .true.
#   direct_netcdf_write      = .true.
#   init_time_days           = -1
#   init_time_seconds        = -1
#   obs_seq_in_file_name     = "obs_seq.out"
#   start_from_restart       = .true.
# /
#-----------------------------------------------------------------------------

setenv JOBNAME     cm1_filter
setenv JOBID       $$

setenv ENS_SIZE   3

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
# ENSEMBLEDIR  The location of the initial ensemble of cm1 files
#-----------------------------------------------------------------------------

set     DARTDIR = ${home}/svn/DART/cm1/models/cm1
set      CM1DIR = ${home}/svn/DART/cm1/models/cm1/src/cm1r18.3/run
set ENSEMBLEDIR = ${home}/temp/CM1R18.3

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
# Get the cm1 executable, control files, and data files.
#-----------------------------------------------------------------------------

${COPY} ${DARTDIR}/work/filter                       . || exit 1
${COPY} ${DARTDIR}/work/advance_time                 . || exit 1
${COPY} ${DARTDIR}/work/input.nml  		     . || exit 1
${COPY} ${DARTDIR}/shell_scripts/advance_model.csh   . || exit 1
${COPY} ${DARTDIR}/work/obs_seq.out                  . || exit 1

${COPY} ${CM1DIR}/cm1.exe             	             . || exit 1

# filter will read the contents of input_filelist.txt for the names
# of the CM1 restart files. After each model advance these same
# filenames will be queried for the new model state.
# So, the actual file associated with this name must be updated.
# advance_model.csh is used by both perfect_model_obs and filter.
# Each model advance is performed in a unique model directory that
# has an instance or ensemble member number in the directory name.

${REMOVE} input_filelist.txt
touch     input_filelist.txt

@ instance = 1
while ( $instance <= $ENS_SIZE )

   # filter will read the contents of input_filelist.txt for
   # the file containing the most up-to-date model state.
   # Each ensemble member (instance) will get advanced in its own directory.

   set filename = "cm1out_rst_000001.nc"
   set  dirname = `printf dir_model_%03d $instance`

   mkdir $dirname

   ${COPY} ${ENSEMBLEDIR}/m${instance}/cm1out_rst_000001.nc $dirname/$filename || exit 2
   ${COPY} ${ENSEMBLEDIR}/m${instance}/input_sounding       $dirname || exit 2
   ${COPY} ${ENSEMBLEDIR}/m${instance}/LANDUSE.TBL          $dirname || exit 2
   ${COPY} ${ENSEMBLEDIR}/m${instance}/namelist.input       $dirname/namelist.input.template || exit 2

   # Enforce certain assumptions.
   # 1) we are always starting from a netCDF file ... irst == 1
   # 2) restart always be named cm1_our_rst_000001.nc ... rstnum == 1
   # 3) the run_time value will be a dummy string to facilitate
   #    the use of 'sed' in the advance_model.csh script.

   # Need to modify namelist.input to be a template namelist
   ex $dirname/namelist.input.template <<ex_end
g;irst     ;s;= .*;= 1;
g;rstnum   ;s;= .*;= 1;
g;run_time ;s;= .*;= CM1_FORECAST_LENGTH;
g;rstfrq   ;s;= .*;= CM1_FORECAST_LENGTH;
wq
ex_end

   grep CM1_FORECAST_LENGTH $dirname/namelist.input.template
   if ( $status != 0 ) then
      echo "<$dirname/namelist.input> did not get updated correctly."
      echo "Check the syntax and try again. Aborting."
      exit 2
   endif

   # Might as well put the known ensemble size in the input.nml
   ex input.nml <<ex_end
g;ens_size ;s;= .*;= ${ENS_SIZE};
g;num_output_obs_members ;s;= .*;= ${ENS_SIZE};
g;num_output_state_members ;s;= .*;= ${ENS_SIZE};
wq
ex_end

   # append the name of the restart file into the list
   echo "$dirname/$filename" >> input_filelist.txt

   @ instance++
end

# When filter starts, it needs to read the grid information.
# The DART/CM1 interface specifies the file for the grid information by
# the input.nml:&model_nml:cm1_template_file  variable, which defaults
# to cm1out_rst_000001.nc
#
# Consequently, we can get simply link to any model state.
# Until such time as the 'time' variable is correct, we also need a namelist.input

${LINK} dir_model_001/cm1out_rst_000001.nc  cm1out_rst_000001.nc || exit 3
${COPY} dir_model_001/namelist.input.template     namelist.input || exit 3

echo ""
echo "cd ${CENTRALDIR}"
echo "Make sure that input.nml contents are correct."
echo "Make sure the obs_seq.in and the model state are appropriate."
echo "Make sure that namelist.input contents are correct."
echo "Make sure that input_filelist.txt contains the right file and the file exists."
echo "Launch (mpirun?) ./filter"
echo ""

exit 0


