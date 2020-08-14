#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
#-------------------------------------------------------------------------------
# DESCRIPTION:
#
# Driver script to consolidate the DART run-time output.
# Modify a couple directory names in this script and
# submit to the batch system and it will invoke the
# 'DART_cleanup_add_time.csh' script for each output timestep and
# 'DART_cleanup_pack_members.csh' to consolidate all ensemble members.
#
#-------------------------------------------------------------------------------
# SLRUM directives             sbatch script.csh
#
# sinfo     information about the whole slurm system
# squeue    information about running jobs
# sbatch    submit a job
# scancel   kill a job
#
#SBATCH --ignore-pbs
#SBATCH --job-name=DART_cleanup
#SBATCH -n 32
#SBATCH --ntasks-per-node=32
#SBATCH --time=01:00:00
#SBATCH -A project
#SBATCH -p dav
#SBATCH -C caldera
#SBATCH -e DART_cleanup.%j.err
#SBATCH -o DART_cleanup.%j.out
#
#-------------------------------------------------------------------------------
# PBS directives                qsub script.csh
#
# qstat    information on running jobs
# qsub     submit a job
# qdel     kill a job
# qpeek    see output from a running job
#
#PBS -N DART_cleanup         
#PBS -l walltime=00:10:00
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
#BSUB -J DART_cleanup
#BSUB -o DART_cleanup.%J.log
#BSUB -q regular
#BSUB -n 16
#BSUB -W 0:30:00
#BSUB -P project
#
#===============================================================================

# USER SETTINGS HERE

set DARTdir = ${HOME}/WRF_Hydro/wrf_hydro_dart/models/wrf_hydro
set HydroDARTdir = /glade/scratch/${USER}/wrfhydro_dart/flo_cut/runs/da_ln1_n80_op2_lp4_sp4_noloc_ga_ana_sample_var
set MARKER = all

# END USER SETTINGS

set nonomatch

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
  setenv NJOBS 36
else if ($?LSB_HOSTS) then
  echo running LSF
  setenv LSF true
  setenv MP_PGMMODEL mpmd
  setenv RUNCMD "mpirun.lsf -cmdfile"
  setenv NJOBS 16
else
  echo running without a batch system
  setenv NOBATCH true
  setenv RUNCMD "echo"
  setenv NJOBS 36
endif

setenv TMPDIR /glade/scratch/$USER/temp
mkdir -p $TMPDIR

cd ${HydroDARTdir} || exit 1

ln -sf ${DARTdir}/shell_scripts/DART_cleanup_add_time.csh     . || exit 2
ln -sf ${DARTdir}/shell_scripts/DART_cleanup_pack_members.csh . || exit 2

echo "DART time support started at "`date`

#===============================================================================
# There are some awful situations that come up if some partial results
# from previous attempts are lying about. Clean up!

\rm -f `find . -name "*.out.nc"`
\rm -f `find . -name "*.out.nc*.tmp"`
\rm -f ${MARKER}*

#===============================================================================
# how many dates are going to be processed 
# This is the number of directories in the output directory

find output -type d | grep 'output/' | sort >! dart_output_file_list.txt

@ timesteps = `cat dart_output_file_list.txt | wc -l`

# loop over all the timesteps
set d=1
while ( $d <= $timesteps )

  \rm -f concatcommands

  set j=1
  while ( $j <= $NJOBS && $d <= $timesteps)

    set DIRECTORY = `head -n $d dart_output_file_list.txt | tail -n 1`
    set TIMESTAMP = ${DIRECTORY:t}

    # create a command file where each line calls a script with
    # unique arguments, including a unique work directory name.
    echo "csh ./DART_cleanup_add_time.csh $TIMESTAMP" >> concatcommands
  
    # advance the loop counter
    @ d ++

    # advance the concurrent job counter
    @ j ++
  end

  # need the same number of commands as there are cpus on the node.
  # pad the script with a call to 'date'

  set j=`cat concatcommands | wc -l`
  while ( $j < $NJOBS )
    echo "date " >> concatcommands
    @ j ++
  end

  echo "running jobs:"

  cat concatcommands

  # actually launch the jobs here
  ${RUNCMD} ./concatcommands 

end

\rm -f concatcommands dart_output_file_list.txt

echo "DART time support ended at "`date`

#===============================================================================
# Concatenate individual files from each stage together
#===============================================================================

foreach STAGE ( input preassim analysis output )
   foreach TYPE (mean     sd     priorinf_mean     priorinf_sd     postinf_mean     postinf_sd     \
                 mean_d01 sd_d01 priorinf_mean_d01 priorinf_sd_d01 postinf_mean_d01 postinf_sd_d01 \
                 mean_d02 sd_d02 priorinf_mean_d02 priorinf_sd_d02 postinf_mean_d02 postinf_sd_d02 )

      set   OUTPUT = ${MARKER}_${STAGE}_${TYPE}.nc     
      ls -1 output/*/${STAGE}_${TYPE}.*.out.nc | sort >! concatlist.txt
      if ( $status == 0 ) then
         echo -n "Trying to create $OUTPUT ..."
         cat concatlist.txt | ncrcat -O -h -H ${OUTPUT} || exit 3
         echo " done."
      endif

   end
end

\rm -f concatlist.txt

echo "DART diagnostics summarized at "`date`

#===============================================================================
# Concatenate all ensemble members together
#===============================================================================
# This consolidates all the timesteps for each ensemble member into one.
# There will be ens_size files at the end of this

@ ens_size = `find . -name "member_*" | wc -l`
@ member = 1

while ( $member <= $ens_size )

  \rm -f concatcommands

  set j=1
  while ( $j <= $NJOBS && $member <= $ens_size)

    # create a command file where each line calls a script with
    # unique arguments, including a unique work directory name.
    echo "csh ./DART_cleanup_pack_members.csh $member" >> concatcommands
  
    @ member ++
    @ j ++
  end

  # need the same number of commands as there are cpus on the node.
  # pad the script with a call to 'date'

  set j=`cat concatcommands | wc -l`
  while ( $j < $NJOBS )
    echo "date " >> concatcommands
    @ j ++
  end

  echo "running jobs:"

  cat concatcommands

  # actually launch the jobs here
  ${RUNCMD} ./concatcommands 

end

\rm -f concatcommands

echo "DART members created at "`date`

#===============================================================================
# So by here we have a bunch of *member_????.nc files that have the entire
# timeseries for each member.
#===============================================================================

cat << EndOfFile >! member.cdl
netcdf member.template {
dimensions:
        member = 100 ;
variables:
        int member(member) ;
                member:long_name = "ensemble member number" ;
// global attributes:
                :version = "none" ;
data:
 member =      1,  2,  3,  4,  5,  6,  7,  8,  9,
          10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
          20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
          30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
          40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
          50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
          60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
          70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
          80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
          90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100 ;
}
EndOfFile
ncgen -o nuisance.nc member.cdl || exit 4

@ nuisance = $ens_size - 1
ncks -O -d member,0,$nuisance -v member nuisance.nc member.nc

\rm -f nuisance.nc member.cdl

#===============================================================================
# This consolidates all the DART ensemble members into one 
# and necessarily changes the unlimited dimension to be 'member'.
# The resulting variable shape is (member, time, links)

foreach STAGE ( input preassim analysis output )

   ls -1 ${STAGE}_member_*.nc | grep _d01 | sort >! DART_cleanup_sorted_list.txt
   if ($status == 0) then

      set FILE = ${MARKER}_${STAGE}_ensemble_d01.nc
      cat DART_cleanup_sorted_list.txt | ncecat -O -h -H -u member -o ${FILE} || exit 5

      # This permutes the variable shape to be (time, member, links)
      # and restores 'time' as the unlimited dimension.
      ncpdq -a time,member ${FILE} ${FILE}.temp || exit 6
      ncks -h -H -A member.nc ${FILE}.temp      || exit 7
      \mv ${FILE}.temp ${FILE}

      \rm -f `cat DART_cleanup_sorted_list.txt`

   endif


   ls -1 ${STAGE}_member_*.nc | grep _d02 | sort >! DART_cleanup_sorted_list.txt
   if ($status == 0) then

      set FILE = ${MARKER}_${STAGE}_ensemble_d02.nc
      cat DART_cleanup_sorted_list.txt | ncecat -O -h -H -u member -o ${FILE} || exit 5

      ncpdq -a time,member ${FILE} ${FILE}.temp || exit 6
      ncks -h -H -A member.nc ${FILE}.temp      || exit 7
      \mv ${FILE}.temp ${FILE}

      \rm -f `cat DART_cleanup_sorted_list.txt`

   endif


   ls -1 ${STAGE}_member_*.nc | grep -v _d0 | sort >! DART_cleanup_sorted_list.txt
   if ($status == 0) then

      set FILE = ${MARKER}_${STAGE}_ensemble.nc
      cat DART_cleanup_sorted_list.txt | ncecat -O -h -H -u member -o ${FILE} || exit 5

      ncpdq -a time,member ${FILE} ${FILE}.temp || exit 6
      ncks -h -H -A member.nc ${FILE}.temp      || exit 7
      \mv ${FILE}.temp ${FILE}

      \rm -f `cat DART_cleanup_sorted_list.txt`

   endif

end

\rm -f member.nc DART_cleanup_sorted_list.txt

echo "job ended at "`date`

exit 0


