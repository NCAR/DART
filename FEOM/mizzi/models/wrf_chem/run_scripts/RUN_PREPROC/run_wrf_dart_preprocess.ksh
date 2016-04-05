#!/bin/ksh -aeux 
#
# set accounting info
export PROJ_NUMBER=P19010000
#
# set switches
  export DO_WRF_DART_PREPROC=true
  export NL_APM_SCALE=1.
  export NL_APM_SCALE_SW=.FALSE.
#
# set time data (greg day, sec: 148805, 64800)
  export START_DATE=2008060106
  export END_DATE=2008063018
#  export START_DATE=2008070100
#  export END_DATE=2008073118
  export TIME_INC=6
  export ASIM_WINDOW=3
#
# set use switches
  export USE_HSI=false
#
# set versions
  export DART_VER=DART_CHEM
  export WRFDA_VER=WRFDAv3.4_dmpar
  export WRF_VER=WRFv3.4_dmpar
#
# independent directories
  export RUN_DIR=/glade/scratch/mizzi/WRF_OBS_PREPROCESS
  export CODE_DIR=/glade/p/work/mizzi/TRUNK
  export DATA_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA
  export HSI_DATA_DIR=/MIZZI/AVE_TEST_DATA
#
# dependent directories
#  export INPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Std_SVD
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Std_SVD
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Std_SVD_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Std_SVD_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Trc_Ret
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Trc_Ret
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Trc_Ret_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Trc_Ret_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_Rev
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_Rev
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_Rev_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_Rev_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_no_rot1
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_no_rot1
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_no_rot1_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_no_rot1_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_no_rot1_bloc
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_no_rot1_bloc
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_no_rot1_bloc_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_no_rot1_bloc_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_bloc
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_bloc
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_bloc_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_bloc_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_scale
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_scale
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_scale_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_scale_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_DBL
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_DBL
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_DBL_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_DBL_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_DBL_bloc
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_bloc_DBL
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_DBL_bloc_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_DBL_bloc_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_DBL_Rev
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_DBL_Rev
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_DBL_Rev_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_DBL_Rev_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_Rev_bloc
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_Rev_bloc
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPCOMB_Mig_DA_Rev_bloc_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPCOMB_Mig_DA_Rev_bloc_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_MOPnIAS_CO_Mig_DA
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPnIAS_CO_Mig_DA
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPnIAS_CO_Mig_DA_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPnIAS_CO_Mig_DA_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_IASCOMB_O3_Mig_DA
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_IASCOMB_O3_Mig_DA
#  export OUTPUT_DIR=${DATA_DIR}/obs_IASCOMB_O3_Mig_DA_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_IASCOMB_O3_Mig_DA_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_MOPnIAS_CO_Mig_DA_bloc
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPnIAS_CO_Mig_DA_bloc
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPnIAS_CO_Mig_DA_bloc_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPnIAS_CO_Mig_DA_bloc_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_MOPnIAS_CO_Mig_DA_DBL
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPnIAS_CO_Mig_DA_DBL
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPnIAS_CO_Mig_DA_DBL_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPnIAS_CO_Mig_DA_DBL_filt
#
#  export INPUT_DIR=${DATA_DIR}/obs_MOPnIAS_CO_Mig_DA_DBL_bloc
#  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPnIAS_CO_Mig_DA_DBL_bloc
#  export OUTPUT_DIR=${DATA_DIR}/obs_MOPnIAS_CO_Mig_DA_DBL_bloc_filt
#  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPnIAS_CO_Mig_DA_DBL_bloc_filt
#
  export INPUT_DIR=${DATA_DIR}/obs_MOPnIAS_CO_O3_Mig_DA
  export HSI_INPUT_DIR=${HSI_DATA_DIR}/obs_MOPnIAS_CO_O3_Mig_DA
  export OUTPUT_DIR=${DATA_DIR}/obs_MOPnIAS_CO_O3_Mig_DA_filt
  export HSI_OUTPUT_DIR=${HSI_DATA_DIR}/obs_MOPnIAS_CO_O3_Mig_DA_filt
#
  export WRF_DIR=${CODE_DIR}/${WRF_VER}
  export WRFDA_DIR=${CODE_DIR}/${WRFDA_VER}
  export BUILD_DIR=${WRFDA_DIR}/var/build
  export DART_DIR=${CODE_DIR}/${DART_VER}
  export DART_WRF_CHEM=${DART_DIR}/models/wrf_chem/work
  mkdir -p ${OUTPUT_DIR}
#  hsi "mkdir -p ${HSI_OUTPUT_DIR}"
#
# make run directory and go to it
  if [[ ! -d ${RUN_DIR} ]]; then 
     mkdir -p ${RUN_DIR}
     cd ${RUN_DIR}
  else
     cd ${RUN_DIR}
  fi
  cp ${DART_WRF_CHEM}/advance_time ./.
  cp ${DART_WRF_CHEM}/input.nml ./.
#
# begin time loop
  export L_DATE=${START_DATE}
  export YYYY=$(echo $L_DATE | cut -c1-4)
  export MM=$(echo $L_DATE | cut -c5-6)
  export DD=$(echo $L_DATE | cut -c7-8)
  export HH=$(echo $L_DATE | cut -c9-10)
  export FILE_DATE=${YYYY}-${MM}-${DD}_${HH}:00:00
#
# copy template wrfinput data
  if ${USE_HSI}; then
     hsi get wrfinput_d01 : ${HSI_DATA_DIR}/wpb_rc_100km/${L_DATE}/wrfinput_d01_${FILE_DATE}.e001
  else
     cp ${DATA_DIR}/wpb_rc_100km/${L_DATE}/wrfinput_d01_${FILE_DATE}.e001 wrfinput_d01
  fi
#
# loop over analsis dates
  while [[ ${L_DATE} -le ${END_DATE} ]]; do
     export YYYY=$(echo $L_DATE | cut -c1-4)
     export MM=$(echo $L_DATE | cut -c5-6)
     export DD=$(echo $L_DATE | cut -c7-8)
     export HH=$(echo $L_DATE | cut -c9-10)
     export FILE_DATE=${YYYY}-${MM}-${DD}_${HH}:00:00
#
     set -A GREG_DATA `echo $L_DATE 0 -g | ./advance_time`
     export DAY_GREG=${GREG_DATA[0]}
     export SEC_GREG=${GREG_DATA[1]}
#
# copy DART utilities
     cp ${DART_DIR}/models/wrf_chem/work/wrf_dart_obs_preprocess ./.
     cp ${DART_DIR}/models/wrf_chem/WRF_DART_utilities/wrf_dart_obs_preprocess.nml ./.
     cp ${DART_DIR}/models/wrf_chem/work/input.nml ./.
#
# copy input data
     if ${USE_HSI}; then
        hsi get obs_seq.old : ${HSI_INPUT_DIR}/${L_DATE}/obs_seq_comb_${L_DATE}.out
     else
        cp ${INPUT_DIR}/${L_DATE}/obs_seq_comb_${L_DATE}.out obs_seq.old
     fi
#
# Make obs_def_apm_nml for apm_scale to adjust observation error variance
        rm -rf obs_def_apm.nml
        cat <<EOF > obs_def_apm.nml
&obs_def_apm_nml
apm_scale=${NL_APM_SCALE}
apm_scale_sw=${NL_APM_SCALE_SW}
/
EOF
#
# Create job script 
     if [[ -f job.ksh ]]; then rm -rf job.ksh; fi
     touch job.ksh
     RANDOM=$$
     export JOBRND=pre_$RANDOM
     cat << EOF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -n 1                                  # number of total (MPI) tas    758 ks
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W 00:10                              # wallclock time (minutes)
#BSUB -q geyser
#
# Run wrf_to_dart
rm -rf pre_*.err
rm -rf pre_*.out
./wrf_dart_obs_preprocess ${DAY_GREG} ${SEC_GREG} > index_preprocess 2>&1 
#
export RC=\$?     
if [[ -f PRE_SUCCESS ]]; then rm -rf PRE_SUCCESS; fi     
if [[ -f PRE_FAILED ]]; then rm -rf PRE_FAILED; fi          
if [[ \$RC = 0 ]]; then
   touch PRE_SUCCESS
else
   touch PRE_FAILED 
   exit
fi
EOF
#
# Submit convert file script for each and wait until job completes
     bsub -K < job.ksh 
#
# save output
     if ${USE_HSI}; then
#        hsi "mkdir -p ${HSI_OUTPUT_DIR}/${L_DATE}"
#        hsi put obs_seq.new : ${HSI_OUTPUT_DIR}/${L_DATE}/obs_seq_comb_filtered_${L_DATE}.out
        mkdir -p ${OUTPUT_DIR}/${L_DATE}
        cp obs_seq.new ${OUTPUT_DIR}/${L_DATE}/obs_seq_comb_filtered_${L_DATE}.out 
     else
        mkdir -p ${OUTPUT_DIR}/${L_DATE}
        cp obs_seq.new ${OUTPUT_DIR}/${L_DATE}/obs_seq_comb_filtered_${L_DATE}.out 
     fi
#
# clean up work directory
     rm dart_log.*
     rm obs_seq.*
#
# loop to next cycle time
     export P_DATE=${L_DATE}
     export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${TIME_INC} 2>/dev/null)  
  done 
  rm advance_time
  rm wrf_dart_obs_preprocess*
  rm input.nml
  rm wrfinput*
exit
