#!/bin/csh
#
# Example script to run the gsi_to_dart converter
# The converter requires more mpi tasks than ensemble members, so choose
# the number of nodes and tasks accordingly.
# 
# Each queuing system will have their own way of launching an mpi executable
# so you will surely need to make sure the examples work on your system.
#
# This is designed to be submitted from the directory that has the
# gsi_to_dart executable and namelist.
#==========================================================================
# SLURM directives              sbatch submit.csh
#
# sinfo     information about the whole slurm system
# squeue    information about running jobs
# sbatch    submitting a job
# scancel   killing a job
#
#SBAXXX --ignore-pbs
#SBAXXX --job-name=gsi2dart
#SBAXXX --ntasks=81
#SBAXXX --ntasks-per-node=36
#SBAXXX --time=00:03:00
#SBAXXX -A P########
#SBAXXX -p dav
#SBAXXX -e testslurm.%j.err
#SBAXXX -o testslurm.%j.out
#SBAXXX --mail-type=END
#SBAXXX --mail-type=FAIL
#SBAXXX --mail-user=${USER}@ucar.edu
#
#==========================================================================
# PBS directives                qsub submit.csh
#
#PBS -N gsi2dart         
#PBS -l walltime=00:03:00
#PBS -q economy
#PBS -l select=3:ncpus=36:mpiprocs=36
#PBS -A P########
#PBS -j oe
#PBS -k eod
#PBS -m abe
#PBS -M ${USER}@ucar.edu
#
#==========================================================================
# LSF directives                bsub < submit.csh
#
#BSUB -n 81
#BSUB -J gsi2dart
#BSUB -o output.driver
#BSUB -e output.driver
#BSUB -q regular
#BSUB -P P########
#BSUB -W 15
#BSUB -R "span[ptile=16]"
#BSUB -N -u ${USER}@ucar.edu

# machine-specific dereferencing
if ($?SLURM_JOB_ID) then

   setenv ORIGINALDIR $SLURM_SUBMIT_DIR
   setenv     JOBNAME $SLURM_JOB_NAME
   setenv       JOBID $SLURM_JOBID
   setenv     MYQUEUE $SLURM_JOB_PARTITION
   setenv      MYHOST $SLURM_SUBMIT_HOST
   setenv   NODENAMES $SLURM_NODELIST
   setenv   LAUNCHCMD "mpirun -np $SLURM_NTASKS -bind-to core"

else if ($?PBS_NODEFILE) then

   # PBS environment variables:
   # env | grep PBS | sort

   setenv ORIGINALDIR $PBS_O_WORKDIR
   setenv     JOBNAME $PBS_JOBNAME
   setenv       JOBID $PBS_JOBID
   setenv     MYQUEUE $PBS_O_QUEUE
   setenv      MYHOST $PBS_O_HOST

   setenv NPROCS    `cat $PBS_NODEFILE | wc -l`
   setenv LAUNCHCMD "mpirun -np $NPROCS"

else if ($?LSB_HOSTS) then

   setenv ORIGINALDIR $LS_SUBCWD
   setenv     JOBNAME $LSB_OUTPUTFILE:ar
   setenv       JOBID $LSB_JOBID
   setenv ARRAY_INDEX $LSB_JOBINDEX
   setenv  LAST_INDEX $LSB_JOBINDEX_END
   setenv     MYQUEUE $LSB_QUEUE
   setenv      MYHOST $LSB_SUB_HOST
   setenv   NODENAMES ${LSB_HOSTS}
   setenv   LAUNCHCMD "mpirun.lsf"

endif

#----------------------------------------------------------------------
# Just an echo of job attributes
#----------------------------------------------------------------------

echo
echo "${JOBNAME} ($JOBID) submitted   from $ORIGINALDIR"
echo "${JOBNAME} ($JOBID) submitted   from $MYHOST"
echo "${JOBNAME} ($JOBID) running in queue $MYQUEUE"
echo "${JOBNAME} ($JOBID) running       on $HOST"

cd ${ORIGINALDIR}

${LAUNCHCMD} ./gsi_to_dart

exit 0


