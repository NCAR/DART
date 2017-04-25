#!/bin/tcsh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Based on scripts to control WRF-DART provided by 
# William James Schouler Miller of the University of Maryland.
# Thanks William!
#
#==========================================================================
#
#  this script submits the jobs to run the WRF ensemble members 
#
#  this script is dependent on successful completion of Filter   
#
#SBATCH --array=1-60%60
#SBATCH -n 20 
#SBATCH --mem-per-cpu=6000
#SBATCH -t 12:00:00
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wmiller@atmos.umd.edu

#load compiler configuration used to build DART
module unload openmpi
module unload gnu
sleep 1
module load intel/2015.0.3.032
module load openmpi/intel/1.8.6

set nprocs = 20

set paramfile = init_param.csh
source $paramfile
#EDIT#
  #  First determine the appropriate analysis date
  if ( $#argv > 0 ) then
    set initial_date   = ${1}
    echo "submitted jobs for running ensembles initialized at $initial_date"
  else
    echo "please enter a date: yyyymmddhh"
    exit
  endif
# initial_date is of the form yyyymmddhh

cd $RUN_DIR

#  Integrate ensemble members to next analysis time

set NUM_ENS = 60 

set n = 1
while ( $n <= $NUM_ENS )

 if ( $SLURM_ARRAY_TASK_ID == $n ) then

   echo "  STARTING ENSEMBLE MEMBER $n"

   set ensstring = `printf "%04d" $n`

     echo `module list`

     set gdate  = (`echo $initial_date 0h -g | ${DART_DIR}/models/wrf/work/advance_time`)
     set gdatef = (`echo $initial_date ${ASSIM_INT_HOURS}h -g | ${DART_DIR}/models/wrf/work/advance_time`)
     set wdate  =  `echo $initial_date 0h -w | ${DART_DIR}/models/wrf/work/advance_time`
     set yyyy   = `echo $initial_date | cut -b1-4`
     set mm     = `echo $initial_date | cut -b5-6`
     set dd     = `echo $initial_date | cut -b7-8`
     set hh     = `echo $initial_date | cut -b9-10`

   if ( -e ${RUN_DIR}/run_ensmem_${n}.csh ) then
     ${REMOVE} ${RUN_DIR}/run_ensmem_${n}.csh
   endif
   touch ${RUN_DIR}/run_ensmem_${n}.csh

   cat >> ${RUN_DIR}/run_ensmem_${n}.csh << EOF
#!/bin/csh

   mv ${RUN_DIR}/filter_restart.${ensstring} ${RUN_DIR}/filter_ic_new.${ensstring}

   cd $RUN_DIR
   ./advance_mem_restart.csh $initial_date $n $paramfile $nprocs

   if ( -e ${RUN_DIR}/assim_model_state_ud.${ensstring} ) then

      ${MOVE} assim_model_state_ud.${ensstring}             ${OUTPUT_DIR}/${initial_date}/DART/filter_ic.${ensstring}
      ${MOVE} WRFOUT/wrf.out_${gdatef[1]}_${gdatef[2]}_${n} ${OUTPUT_DIR}/${initial_date}/.
      ${REMOVE} filter_ic_new.${ensstring} start_member_${n} run_ensmem_${n}.csh
      rm -rf advance_temp${n}  
 
   endif
EOF

   chmod u+x run_ensmem_${n}.csh
   ./run_ensmem_${n}.csh

   endif

   @ n++

sleep 1

end


exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

