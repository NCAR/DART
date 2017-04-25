#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Based on scripts to control WRF-DART provided by 
# William James Schouler Miller of the University of Maryland.
# Thanks William!
#
#=========================================================================
#
#SBATCH -n 40
#SBATCH --mem-per-cpu=2000
#SBATCH -t 01:00:00
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wmiller@atmos.umd.edu

set nprocs = 40

#load compiler configuration used to build DART
module unload gnu
module unload openmpi
module load intel/2015.0.3.032
module load openmpi/intel/1.8.6

  set paramfile = init_param.csh
  source $paramfile

  #  First determine the appropriate analysis date
  if ( $#argv > 0 ) then
    set datea   = ${1}
    setenv restore 1   #  set the restore variable
    echo 'starting a restore'
  else
    echo "please enter a date: yyyymmddhh"
    setenv restore 0   #  set the restore variable
    exit
  endif

cd $RUN_DIR

     # If the current output directory does not exist, and this is a restore,
     #   then the analyses have been caught up, so exit.
     if ( ! -d ${OUTPUT_DIR}/${datea} && $restore == 1 ) then
        ${REMOVE} ${RUN_DIR}/ABORT_RETRO
        echo 'exiting because output directory does not exist and this is a restore'
        exit
     endif

     set datep  = `echo $datea -${ASSIM_INT_HOURS}   | ${DART_DIR}/models/wrf/work/advance_time`
     set gdate  = `echo $datea 0 -g                  | ${DART_DIR}/models/wrf/work/advance_time`
     set gdatef = `echo $datea ${ASSIM_INT_HOURS} -g | ${DART_DIR}/models/wrf/work/advance_time`
     set wdate  = `echo $datea 0 -w                  | ${DART_DIR}/models/wrf/work/advance_time`
     set hh     = `echo $datea | cut -b9-10`

     echo 'ready to check inputs'
     #  Check to make sure all input data exists
     foreach infile ( wrfinput_d01_${gdate[1]}_${gdate[2]}_mean \
                      wrfbdy_d01_${gdatef[1]}_${gdatef[2]}_mean obs_seq.out )

        if ( ! -e ${OUTPUT_DIR}/${datea}/${infile} ) then

           echo "${OUTPUT_DIR}/${datea}/${infile} is missing!  Stopping the system"
           touch ABORT_RETRO
           exit

        endif
     #  Copy the appropriate LSM, inflation, input.nml files
     set n = 1
     while ( $n <= $NUM_ENS )

        set ensstring = `echo $n + 10000 | bc | cut -b2-5`
        if ( -e ${OUTPUT_DIR}/${datep}/DART/filter_ic.${ensstring} ) then
           ${LINK} ${OUTPUT_DIR}/${datep}/DART/filter_ic.${ensstring} filter_ic.${ensstring}
        else
           echo "${OUTPUT_DIR}/${datep}/DART/filter_ic.${ensstring} is missing!  Stopping the system"
           touch ABORT_RETRO
           exit
        endif
        @ n++

     end

     #  Get wrfinput source information
     if ( -e ${RUN_DIR}/wrfinput_d01 ) rm -f ${RUN_DIR}/wrfinput_d01 
     ${COPY} ${OUTPUT_DIR}/${datea}/wrfinput_d01_${gdate[1]}_${gdate[2]}_mean wrfinput_d01

     #  Copy the inflation files from the previous time, update for domains
     if ( -e ${OUTPUT_DIR}/${datep}/DART/prior_inf_ic ) then
        ${COPY} ${OUTPUT_DIR}/${datep}/DART/prior_inf_ic prior_inf_ic_old
    endif

     ${LINK} ${OUTPUT_DIR}/${datea}/obs_seq.out .

     ${REMOVE} ${RUN_DIR}/WRF
     ${LINK} ${OUTPUT_DIR}/${datea} ${RUN_DIR}/WRF
     #  run filter to generate the analysis

#  run data assimilation system
mpirun -n $nprocs filter
if ( -e ${RUN_DIR}/obs_seq.final )  touch ${RUN_DIR}/filter_done

## verfify successful Filter run and set up run directory for ensemble advances
set initial_date = $datea 
set gdate  = (`echo $initial_date 0h -g | ${DART_DIR}/models/wrf/work/advance_time`)
set gdatef = (`echo $initial_date ${ASSIM_INT_HOURS}h -g | ${DART_DIR}/models/wrf/work/advance_time`)
set wdate  =  `echo $initial_date 0h -w | ${DART_DIR}/models/wrf/work/advance_time`
set yyyy   = `echo $initial_date | cut -b1-4`
set mm     = `echo $initial_date | cut -b5-6`
set dd     = `echo $initial_date | cut -b7-8`
set hh     = `echo $initial_date | cut -b9-10`

${REMOVE} ${RUN_DIR}/WRF
${LINK} ${OUTPUT_DIR}/${initial_date} WRF

if ( -e wrfinput_d01 ) rm -f wrfinput_d01
ln -sf ${OUTPUT_DIR}/${initial_date}/wrfinput_d01_${gdate[1]}_${gdate[2]}_mean wrfinput_d01

#-------------------------------------
#   first clean up ${RUN_DIR} after previous Filter run
     rm -f filter_ic* prior_inf_ic_old obs_seq.out

#    verify that Filter has run successfully for current analysis time

     echo "Listing contents of rundir before archiving"
     ls -l *.nc blown* dart_log* filter_* input.nml obs_seq* *inf_ic*
     mkdir -p ${OUTPUT_DIR}/${initial_date}/{DART}

     #  Move diagnostic and obs_seq.final data to storage directories
     foreach FILE ( Prior_Diag.nc Posterior_Diag.nc obs_seq.final )
        if ( -e $FILE && ! -z $FILE ) then
           ${MOVE} $FILE ${OUTPUT_DIR}/${initial_date}/$FILE
           if ( ! $status == 0 ) then
              echo "failed moving ${RUN_DIR}/${FILE}"
              touch BOMBED
           endif
        else
           echo "${OUTPUT_DIR}/${FILE} does not exist and should."
           ls -l
           touch BOMBED
        endif
     end

     #  Move Filter log files to storage directories
     mv dart_log* ${OUTPUT_DIR}/${initial_date} 

     #  Move inflation files to storage directories
     foreach FILE ( prior_inf_ic post_inf_ic )
        if ( -e ${FILE}_new && ! -z ${FILE}_new ) then
           ${MOVE} ${FILE}_new ${OUTPUT_DIR}/${inintial_date}/DART/${FILE}
           if ( ! $status == 0 ) then
              echo "failed moving ${RUN_DIR}/${FILE}"
              touch BOMBED
           endif
        endif
     end

     if ( -e BOMBED ) then

        echo "FATAL SYSTEM ERROR"
        touch ABORT_RETRO
        ${REMOVE} BOMBED
        exit

     endif

#  Create time info file for advance_model.csh
if ( -e ${RUN_DIR}/initial_adv_info ) rm -f ${RUN_DIR}/initial_adv_info
touch ${RUN_DIR}/initial_adv_info
echo ${gdatef[2]} ${gdatef[1]} >> ${RUN_DIR}/initial_adv_info
echo ${gdate[2]} ${gdate[1]} >> ${RUN_DIR}/initial_adv_info
echo ${yyyy} ${mm} ${dd} ${hh} "00" "00"  >> ${RUN_DIR}/initial_adv_info
echo ${domains}  >> ${RUN_DIR}/initial_adv_info
echo "./wrf.exe" >> ${RUN_DIR}/initial_adv_info

echo "------------------------------------------"
echo "finished ensemble update for $initial_date"   
#==========================================================================

