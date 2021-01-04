#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
#-------------------------------------------------------------------------------
# DESCRIPTION:
#
# Driver script to consolidate the wrf_hydro run-time output.
# Submit this script to your batch system and it will invoke 
# the 'wrf_hydro_cleanup_worker.csh' script once for each ensemble member.
#
#-------------------------------------------------------------------------------

#===============================================================================
# SLRUM directives             sbatch script.csh
#
# sinfo     information about the whole slurm system
# squeue    information about running jobs
# sbatch    submit a job
# scancel   kill a job
#
#SBATCH --ignore-pbs
#SBATCH --job-name=summ
#SBATCH -n 32
#SBATCH --ntasks-per-node=32
#SBATCH --time=01:00:00
#SBATCH -A project
#SBATCH -p dav
#SBATCH -C caldera
#SBATCH -e summ.%j.err
#SBATCH -o summ.%j.out
#
#===============================================================================
# PBS directives                qsub script.csh
#
# qstat    information on running jobs
# qsub     submit a job
# qdel     kill a job
# qpeek    see output from a running job
#
#PBS -N summ         
#PBS -l walltime=01:00:00
#PBS -q regular
#PBS -l select=1:ncpus=36:mpiprocs=36
#PBS -A NRAL0017
#
#===============================================================================
# LSF directives                bsub < script.csh
#
# bstat    information on running jobs
# bsub     submit a job
# bdel     kill a job
# bpeek    see output from a running job
#
#BSUB -J summ
#BSUB -o summ.%J.log
#BSUB -q regular
#BSUB -n 16
#BSUB -W 0:10:00
#BSUB -P project
#
#===============================================================================

# USER SETTINGS HERE

set DARTdir = /glade/work/${USER}/DART/wrf_hydro_dart_git/wrf_hydro_dart/models/wrf_hydro
set HydroDARTdir = /glade/scratch/${USER}/wrfhydro_dart/sixmile/runs/test1

# END USER SETTINGS


# set things that vary between batch systems here.

if ($?SLURM_JOB_ID) then
  echo running SLURM
  setenv SLURM true
  source /glade/u/apps/opt/slurm_init/init.csh
  setenv RUNCMD "srun --multi-prog"
  setenv NJOBS 32
else if ($?PBS_NODEFILE) then
  echo running PBS
  setenv PBS true
  setenv MPI_SHEPHERD true
  setenv RUNCMD "mpiexec_mpt launch_cf.sh"
  setenv NJOBS 27
else if ($?LSB_HOSTS) then
  echo running LSF
  setenv LSF true
  setenv MP_PGMMODEL mpmd
  setenv RUNCMD "mpirun.lsf -cmdfile"
  setenv NJOBS 16
else
  echo running without a batch system
  setenv NOBATCH true
  setenv RUNCMD "csh"
  setenv NJOBS 1
endif

setenv TMPDIR /glade/scratch/${USER}/temp
mkdir -p $TMPDIR

#===============================================================================

cd ${HydroDARTdir} || exit 1

ln -sf ${DARTdir}/work/advance_time                          . || exit 2
ln -sf ${DARTdir}/shell_scripts/wrf_hydro_cleanup_worker.csh . || exit 2

echo "summarizing started at "`date`

# determine the ensemble size
@ ens_size = `find . -name "member_*" | wc -l`

# loop over each member ... these are referenced [0,N-1]
@ member = 0
while ( $member < $ens_size )

  \rm -f mycmdfile

  @ j = 1
  while ( $j <= $NJOBS && $member < $ens_size )

    # create a command file where each line calls a script with
    # unique arguments, including a unique work directory name.
    echo "csh ./wrf_hydro_cleanup_worker.csh $member" >> mycmdfile
  
    # advance the ensemble member index
    @ member ++

    # advance the concurrent job counter
    @ j ++
  end

  echo "running jobs:"

  cat mycmdfile

  # the system seems to want the same number of commands in each
  # invocation of the command file as there are cpus on the node.
  # if we aren't running an even multiple of real tasks compared
  # to the cpu count, pad the rest of the script with a call to 'date'

  # avoid echoing the filename by making wc read stdin
  set j=`cat mycmdfile | wc -l`
  while ( $j < $NJOBS )
    echo "date " >> mycmdfile
    @ j ++
  end

  # actually launch the jobs here
  ${RUNCMD} ./mycmdfile 

end

\rm -f mycmdfile

echo "summarizing ended at "`date`

#===============================================================================
# Concatenate all ensemble members together
#===============================================================================

foreach STAGE ( CHANOBS_DOMAIN1 CHRTOUT_DOMAIN1 LAKEOUT_DOMAIN1 )

   ls -1 ${STAGE}.*.nc | sort >! wrf_hydro_cleanup_files.txt
   if ( $status == 0 ) then
      
      # Consolidates all the ensemble members into one - unfortunately the
      # 'time' variable has invalid 'valid_min','valid_max' attributes. 

      cat wrf_hydro_cleanup_files.txt | \
          ncecat -O -h -H -u member -o ${STAGE}.nc || exit 1

      # This deletes the nuisance 'time' attributes and then permutes the
      # dimensions so the variables are shaped (time, member, links)
      # and restores 'time' as the unlimited dimension.

      ncatted -O -a valid_min,time,d,,, ${STAGE}.nc bob.nc || exit 5
      ncatted -O -a valid_max,time,d,,, bob.nc bill.nc     || exit 6
      ncpdq   -O -a time,member bill.nc ${STAGE}.nc        || exit 7

      \rm -f bob.nc bill.nc ${STAGE}.*.nc

   endif

   \rm -f wrf_hydro_cleanup_files.txt

end

# HYDRO_RST is a bit different ... cannot use ncatted commands 
# because there are no 'time:valid_[min,max]' attributes
foreach STAGE ( HYDRO_RST )

   ls -1 ${STAGE}.*.nc | sort >! wrf_hydro_cleanup_files.txt
   if ( $status == 0 ) then
      
      # Consolidates all the ensemble members into one - this time there
      # are no 'time:valid_[min,max]' attributes to delete. 

      cat wrf_hydro_cleanup_files.txt | \
          ncecat -O -h -H -u member -o bob.nc || exit 1

      # This permutes the variables to be shaped (time, member, links)
      # and restores 'time' as the unlimited dimension.

      ncpdq -O -a time,member bob.nc ${STAGE}.nc || exit 7

      \rm -f bob.nc ${STAGE}.*.nc

   endif

   \rm -f wrf_hydro_cleanup_files.txt

end

echo job ended at `date`

exit 0


