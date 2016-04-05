#!/bin/ksh -x
###############################################################################
#
#  Wrapper script to run WRF/CHEM-DART cycling experiment
#
############################################################################### 
#
# Define experiment parameters
export W_START_DATE=2008060106
export W_START_DATE=2008060600
export W_END_DATE=2008063018
#export W_END_DATE=2008072418
#export W_END_DATE=2008070300
export W_INITIAL_DATE=2008060106
export W_INITIAL_FILE_DATE=2008-06-01_06:00:00
export W_FIRST_FILTER_DATE=2008060112
export W_CYCLE_PERIOD=6
export W_FCST_PERIOD=6
export W_RUN_INITIAL=false
export W_EXPERIMENT=Ex2_Mg_p10p30_100km_DBL_Rev
export NL_APM_SCALE=1.
export NL_APM_SCALE_SW=.FALSE.
export W_PROJ_NUMBER=P19010000
export W_PROJ_NUMBER=NNES0003
export W_PROJ_NUMBER=NACD0002
#   
# Define directories
export W_PROJECT_DIR=/glade/p/work/mizzi
export W_SCRATCH_DIR=/glade/scratch/mizzi
export W_ACD_DIR=/glade/p/acd/mizzi
export W_HSI_DIR=/MIZZI
export W_TRUNK_DIR=${W_PROJECT_DIR}/TRUNK
export W_DART_DIR=${W_TRUNK_DIR}/DART_CHEM
export W_RUN_SCRIPT_DIR=${W_DART_DIR}/models/wrf_chem/run_scripts/RUN_WRF_CHEM/EXPERIMENT_SCRIPTS
export W_RUN_DIR=${W_SCRATCH_DIR}/RUN_${W_EXPERIMENT}
#
# Setup wrapper run directory
if [[ ! -d ${W_RUN_DIR} ]]; then
   mkdir -p ${W_RUN_DIR}
   cd ${W_RUN_DIR}
else
   cd ${W_RUN_DIR}
fi
#
# Copy files to wrapper run directory
cp ${W_DART_DIR}/models/wrf_chem/work/input.nml ./.
cp ${W_DART_DIR}/models/wrf_chem/work/advance_time ./.
cp ${W_RUN_SCRIPT_DIR}/run_Exp2_Mig_DA_100km_DBL_Revised.ksh run_script.ksh
chmod +x run_script.ksh
#
export W_DATE=${W_START_DATE}
while [[ ${W_DATE} -le ${W_END_DATE} ]]; do
    if ${W_RUN_INITIAL}; then
       export W_RUN_WARM=false
       export W_WARM_FILTER=false
       export W_WARM_WRFCHEM=false
       export W_WARM_ARCHIVE=false 
    elif [[ ${W_DATE} -eq ${W_START_DATE} ]]; then
       export W_RUN_WARM=true
       export W_WARM_FILTER=true
       export W_WARM_WRFCHEM=false
       export W_WARM_ARCHIVE=false
    else
       export W_RUN_WARM=true
       export W_WARM_FILTER=true
       export W_WARM_WRFCHEM=false
       export W_WARM_ARCHIVE=false
    fi
#
    ./run_script.ksh > ${W_EXPERIMENT}_${W_DATE} 2>&1 
#
   export W_RUN_INITIAL=false
   export W_NEXT_DATE=`echo ${W_DATE} +${W_CYCLE_PERIOD}h | ./advance_time`
   export W_DATE=${W_NEXT_DATE}
done
rm -rf input.nml
rm -rf advance_time
rm -rf dart_log.out
exit
