#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# The DART program 'filter' has the ability to take a single model state and 
# perturb it to create an ensemble of model states.
#
# The intent of the this script is to take a single coamps state and make N ensemble
# members using the DART program 'filter'. Each ensemble member will differ from the 
# parent in that noise will be added to each of the variables defined in the state.vars.
# The nature of the noise is defined in the filter namelist. It is not expected that 
# these files are valid for coamps, but they will enable testing of the assimilation.
#
# FIXME: Someone should take the logic from perturb_init.f90 and put it in 
# model_mod:pert_model_copies(). At present, pert_model_copies() uses a default routine.
#
# The steps are basically:
# 1) create an empty directory for the ensemble and copy all the bits needed to
#    create the ensemble into it.
# 2) run trans_coamps_to_dart to create a netCDF file - required by DART
# 3) create an observation sequence file (required by filter) by running pmo
# 4) create a bunch of output files that filter will update.
# 5) convert all those output files into coamps hdf5 restart files
# 6) while we are at it, create a file containing inflation values.
# 7) clean up the garbage?
#
#==========================================================================
# This block of directives constitutes the preamble for the SLURM queuing system.
# the normal way to submit with slurm:  sbatch create_ensemble.csh
# and to check on the job status     :  squeue --users=thoar
#
#SBATCH --ignore-pbs
#SBATCH --job-name=ensemble
#SBATCH --output=ensemble%A.log
#SBATCH --error=ensemble%A.err
#SBATCH --ntasks=16
#SBATCH --time=00:30:00
#SBATCH --error=ensemble%A.err
#SBATCH --output=ensemble%A.log
#==========================================================================
# PBS directives                        qsub test_batch.csh
#                                       qstat -u $USER
#                                       qdel <jobnumber>
#PBS -N filter
#PBS -e filter.err
#PBS -o filter.log
#PBS -l select=1:ncpus=16:mpiprocs=16
#PBS -l walltime=00:30:00
#PBS -A P8685nnnn
#PBS -q economy
#PBS -r n
#
#==========================================================================
# LSF directives                        bsub < advance_ensemble.csh
#                                       bjobs
#                                       bkill <jobnumber>
#BSUB -J filter
#BSUB -o filter.%J.log
#BSUB -P P8685nnnn
#BSUB -q regular
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -W 00:30
#BSUB -N -u ${USER}@ucar.edu
# 
#==========================================================================

if ($?SLURM_JOB_ID) then

   set ORIGINALDIR = $SLURM_SUBMIT_DIR
   set     JOBNAME = $SLURM_JOB_NAME
   set       JOBID = $SLURM_JOBID
   set     MYQUEUE = $SLURM_JOB_PARTITION
   set      MYHOST = $SLURM_SUBMIT_HOST
   set   LAUNCHCMD = "mpirun -np $SLURM_NTASKS -bind-to core"

else if ($?PBS_O_WORKDIR) then

   set ORIGINALDIR = $PBS_O_WORKDIR
   set     JOBNAME = $PBS_JOBNAME
   set       JOBID = $PBS_JOBID
   set      MYHOST = $PBS_O_HOST
   set     MYQUEUE = $PBS_QUEUE
   set   LAUNCHCMD = "mpiexec_mpt"

else if ($?LS_SUBCWD) then

   set ORIGINALDIR = $LS_SUBCWD
   set     JOBNAME = $LSB_JOBNAME
   set       JOBID = $LSB_JOBID
   set     MYQUEUE = $LSB_QUEUE
   set      MYHOST = $LSB_SUB_HOST
   set   LAUNCHCMD = "mpirun.lsf"

else
   echo ""
   echo "This script requires an mpi-aware environment."
   echo ""
   exit 1
endif

#==========================================================================
# Just an echo of job attributes
#==========================================================================

echo
echo "${JOBNAME} ($JOBID) submit directory ${ORIGINALDIR}"
echo "${JOBNAME} ($JOBID) submit      host ${MYHOST}"
echo "${JOBNAME} ($JOBID) running in queue ${MYQUEUE}"
echo "${JOBNAME} ($JOBID) started at "`date`
echo

#==========================================================================
# Provide the required locations - you will have to edit these.
#==========================================================================

set DARTDIR = /home/$USER/DART/coamps
set ORIGIN = /home/nopp/COAMPS_hdf5_files
set DESTINATION = /home/$USER/COAMPS_hdf5_files/Ensemble2
set ENSEMBLE_SIZE = 4
set DTG = 2013011000
set TEMPLATE = coamps_${DTG}.hdf5

#==========================================================================
# Must convert a single state to a netcdf file for filter, perturb it,
# take the output netcdf files and update the hdf5 files ...
#==========================================================================

mkdir -p ${DESTINATION}
cd ${DESTINATION}
\rm -f dart_log.out dart_log.nml

\ln -sf ${DARTDIR}/models/coamps_nest/work/trans_coamps_to_dart .   || exit 1
\ln -sf ${DARTDIR}/models/coamps_nest/work/trans_dart_to_coamps .   || exit 1
\ln -sf ${DARTDIR}/models/coamps_nest/work/perfect_model_obs    .   || exit 1
\ln -sf ${DARTDIR}/models/coamps_nest/work/filter               .   || exit 1
\cp     ${DARTDIR}/models/coamps_nest/work/input.nml            .   || exit 1
\cp     ${DARTDIR}/models/coamps_nest/work/convert.nml          .   || exit 1
\cp     ${DARTDIR}/models/coamps_nest/work/state.vars.small     state.vars  || exit 1
\cp     ${DARTDIR}/models/coamps_nest/work/obs_seq.2obs.in      obs_seq.in  || exit 1
\cp     ${DARTDIR}/assimilation_code/programs/gen_sampling_err_table/work/sampling_error_correction_table.nc .

# Enforce the assumptions.

ex input.nml <<ex_end
g;input_state_files ;s;= .*;= 'dart_vector.nc', 'dart_vector.nc';
g;input_state_file_list ;s;= .*;= '', '';
g;output_state_file_list ;s;= .*;= 'file_list_domain_1.txt', 'file_list_domain_2.txt';
g;ens_size ;s;= .*;= ${ENSEMBLE_SIZE};
g;num_output_obs_members ;s;= .*;= ${ENSEMBLE_SIZE};
g;num_output_state_members ;s;= .*;= ${ENSEMBLE_SIZE};
g;stages_to_write ;s;= .*;= 'input', 'preassim', 'output';
g;output_mean ;s;= .*;= .TRUE.;
g;output_sd ;s;= .*;= .TRUE.;
g;perturb_from_single_instance ;s;= .*;= .TRUE.;
g;sampling_error_correction ;s;= .*;= .FALSE.;
g;inf_flavor ;s;= .*;= 2, 0;
g;inf_initial_from_restart ;s;= .*;= .FALSE., .FALSE.;
g;inf_sd_initial_from_restart ;s;= .*;= .FALSE., .FALSE.;
g;cutoff ;s;= .*;= 0.10;
g;debug ;s;= .*;= 0;
wq
ex_end

# trans_coamps_to_dart creates a netCDF file (dart_vector.nc) from 
# an HDF5 file (coamps.hdf5).  DART reads the variable shapes, etc. from a 
# file called 'coamps.hdf5' but the actual input data files are specified 
# via namelist mechanisms.

\cp  ${ORIGIN}/${TEMPLATE} .
ln -s          ${TEMPLATE} coamps.hdf5

./trans_coamps_to_dart  || exit 2

# Create an observation sequence file (needed for filter)
# obs_seq.2obs.in has precisely 2 observations - one is identically
# part of the COAMPS state (at whatever index 2101 happens to be)
# and the other is a potential temperature on LEVEL 50.0 (wherever that is)
#
# The namelist for pmo reads the two nests, each from the same file: dart_vector.nc
# pmo is not set up to WRITE both nests to the same file when creating the output
# file from scratch, but will happily write to an existing file - so we copy
# the input file to an output file and let pmo overwrite it.

\cp dart_vector.nc perfect_output.nc

./perfect_model_obs || exit 3

#  Create the list of output files for DART - one list per domain. In this case
#  they are both the same since both domains exist in the same input file.
\rm -f file_list_domain_?.txt

set instance = 1
while ( $instance <= ${ENSEMBLE_SIZE} )
   set OUTFILE = `printf dart_%04d_output.nc $instance`
   \cp dart_vector.nc ${OUTFILE}
   echo ${OUTFILE} >> file_list_domain_1.txt
   echo ${OUTFILE} >> file_list_domain_2.txt
   @ instance ++
end

${LAUNCHCMD} ./filter || exit 4

# We now have an ensemble of netCDF files. We need an ensemble of HDF5 files.
# At present, the best way to do this is to copy the template HDF5 file to a 
# bunch of output files and let 'trans_dart_to_coamps' update them.
# 'trans_dart_to_coamps' uses COAMPS write routines that expect a DTG in
# the file name, the dart call provides the base filename of 'CoampsUpdate'

set instance = 1
while ( $instance <= ${ENSEMBLE_SIZE} )

   echo " "
   echo "Converting instance $instance of $ENSEMBLE_SIZE at "`date`

   set  INFILE = `printf dart_%04d_output.nc     $instance`
   set OUTFILE = CoampsUpdate_${DTG}.hdf5
   set RESTART = `printf coamps_%04d_${DTG}.hdf5 $instance`

   ln -sf ${INFILE}  dart_vector.nc
   \cp -v ${TEMPLATE} ${OUTFILE}

   ./trans_dart_to_coamps || exit 5

   \mv -v ${OUTFILE} ${RESTART}

   echo "Finished Converting instance $instance of $ENSEMBLE_SIZE at "`date`
   echo " "

   @ instance ++
end

mv obs_seq.out obs_seq_${DTG}.out

#==========================================================================
# The inflation files have a lot of extra cruft that should not be there.
# The cruft is written by nc_write_prognostic_atts. When it becomes a problem,
# that routine should be rewritten. Until then, the files are a bit bigger
# than necessary. They are intentionally renamed to be output files to mimic
# the behavior during a cycling experiment. Some scripting looks for the
# output of the previous cycle and renames it to be the expected input names
# for the current cycle.

foreach FILE ( preassim_priorinf_mean_d01.nc preassim_priorinf_sd_d01.nc \
               preassim_priorinf_mean_d02.nc preassim_priorinf_sd_d02.nc )

ncks -C -x -v nx_g01,nxm1_g01,ny_g01,nym1_g01,nx_g02,nxm1_g02,ny_g02,nym1_g02,sigw,sigm,latT_g01,lonT_g01,latU_g01,lonU_g01,latV_g01,lonV_g01,latT_g02,lonT_g02,latU_g02,lonU_g02,latV_g02,lonV_g02 \
       $FILE bob.nc
       \mv bob.nc $FILE
end

mv preassim_priorinf_mean_d01.nc output_priorinf_mean_d01.nc
mv preassim_priorinf_mean_d02.nc output_priorinf_mean_d02.nc

# The EXNER functions are static and should not be inflated, so the inflation
# standard deviations should be set to zero.

ncap2 -O -s 'EXBM_g01(:,:,:)=0.0;EXBW_g01(:,:,:)=0.0;THBM_g01(:,:,:)=0.0' \
   preassim_priorinf_sd_d01.nc output_priorinf_sd_d01.nc

ncap2 -O -s 'EXBM_g02(:,:,:)=0.0;EXBW_g02(:,:,:)=0.0;THBM_g02(:,:,:)=0.0' \
   preassim_priorinf_sd_d02.nc output_priorinf_sd_d02.nc

\rm preassim_priorinf_sd_d01.nc preassim_priorinf_sd_d02.nc

#==========================================================================

cat << README >! README.txt

The original file that was perturbed: ${ORIGIN}/${TEMPLATE}

The inflation files are intentionally called 'output_*' because of the logic in
the run_filter.csh script. It looks for the output of the previous cycle and 
renames it to be the input of the next cycle.

README

cat README.txt

#==========================================================================
echo
echo "${JOBNAME} ($JOBID) ended   at "`date`
echo
#==========================================================================

# cleanup block

# \rm dart*.nc dart_log.* file_list*txt filter perfect_model_obs 
# \rm trans_dart_to_coamps trans_coamps_to_dart dart.kstart
# \rm obs_seq.in obs_seq.final sampling_error_correction_table.nc
# \rm coamps.hdf5 perfect_output.nc
# \rm ${TEMPLATE} 
# \rm preassim*.nc

mkdir diags
mv preassim_mean_d??.nc preassim_sd_d??.nc diags
mv   output_mean_d??.nc   output_sd_d??.nc diags

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
