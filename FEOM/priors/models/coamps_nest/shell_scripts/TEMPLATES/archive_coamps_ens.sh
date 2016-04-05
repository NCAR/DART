#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
###############################################################################
#
# AUTHOR:   P. A. Reinecke
#           Naval Research Laboratory
######
#
usage="Usage: `basename $0` -d dtg -e dtg_end -i icycle -b arch_bkgrnd -n arch_init -f arch_fcst -r arch_restart [-W depend]]" 
while getopts ":d:e:i:b:f:r:n:W:" option
do
  case $option in
      d  ) dtg=$OPTARG;;
      e  ) dtg_end=$OPTARG;;
      i  ) icycle_loc=$OPTARG;;
      b  ) arch_bkgrnd=$OPTARG;;
      f  ) arch_fcst=$OPTARG;;
      r  ) arch_restart=$OPTARG;;
      n  ) arch_init=$OPTARG;;
      W  ) depend=$OPTARG;;
      \? ) echo $usage ; exit 65;;
      *  ) echo $usage ; exit 65;;
  esac
done

echo "#"
echo "#################### BEGIN archive_coamps_ens.sh ####################"
echo "#"

BATCH='qsub'
if [[ $depend ]]; then depend="-W depend=afterok:${depend}"; fi

PATH_CONFIG=CONFIG/FILE/GOES/HERE
. ${PATH_CONFIG}

: ${dtg:=$cdtg_beg}
: ${dtg_end:=$dtg}
: ${icycle_loc:=$icycle}
: ${icycle_ocn:=24}
: ${arch_fcst:='f'}
: ${arch_restart:='f'}
: ${arch_init:='t'}
: ${arch_bkgrnd:='t'}

scrDir=${COAMPS_RUN_BASE_DIR}
LOG_DIR=${scrDir}/log
LOG_DTG=${LOG_DIR}/${dtg}
scriptDir=${scrDir}/scripts
DART_BIN=${scrDir}/bin

stageOUT=${scriptDir}/stageOUT_ens.${dtg}.sh

cd ${scrDir}
nnest=`grep nnest ${scrDir}/namelist | sed 's/,//g' | awk -F= '{print $2}'`
dtgp1=`echo ${dtg} ${icycle_loc} | ${DART_BIN}/advance_time`
dtgm1=`echo ${dtg} -${icycle_loc} | ${DART_BIN}/advance_time`
dtgm1_ocn=`echo ${dtgp1} -${icycle_ocn} | ${DART_BIN}/advance_time`
cdtg_beg_p1=`echo ${cdtg_beg}  ${fcstlen_init} | ${DART_BIN}/advance_time`

# Get the machine layout
. ${scriptDir}/HPC_CONFIG.sh
MACH_LAYOUT


if [ -e ${stageOUT} ]; then rm ${stageOUT}; fi
cat > ${stageOUT} << EOF
#!/bin/bash
`PBS_DIRECTIVE STAGE_OUT_ENS ${dtg} ${TRANS_PROCS} ${MIN_PROCS} ${WALLTIME_TRANS} ${LOG_DIR} y n` 

DART_ARC=${DART_ARCH}
ARCMACH=${ARCMACH}
scrDir=${scrDir}
scriptDir=${scriptDir}

dtgm1_ocn=${dtgm1_ocn}
dtgm1=${dtgm1}
dtg=${dtg}
dtg_end=${dtg_end}
dtgp1=${dtgp1}
icycle=${icycle_loc}
icycle_ocn=${icycle_ocn}
digits=5

ENSEMBLE_SIZE=${ENSEMBLE_SIZE}

arch_restart=${arch_restart}
arch_fcst=${arch_fcst}

arch_init=${arch_init}
arch_bkgrnd=${arch_bkgrnd}

cycle_fmt=\`printf %04d0000 \${icycle} \`
cycle_zero=\`printf %04d0000 0 \`

archive mkdir -p \${DART_ARC}

cd \${scrDir}
EOF

cat >> ${stageOUT} << 'EOF'
#-----------------------------------------------------------------------------------
# RESTART FILES
#-----------------------------------------------------------------------------------
if [ $arch_restart = 't' ]; then
(
tar_file=${scrDir}/${dtg}/`printf ${dtg}_restart%03d.tar ${icycle}`
rm -f ${tar_file}
echo $tar_file
for element in `seq 1 ${ENSEMBLE_SIZE}`; do
  dataDir_mbr=`printf %0${digits}d/data ${element}`
  echo $dataDir_mbr
  target_files=`printf ./${dataDir_mbr}/restart*${dtg}%03d* ${icycle}`
  tar rf ${tar_file} ${target_files}
done

archive put -C ${DART_ARC} ${tar_file}
rm -f ${tar_file}
) < /dev/null &
fi

#-----------------------------------------------------------------------------------
# BACKGROUND FIELD FILES
#-----------------------------------------------------------------------------------
if [ $arch_bkgrnd = 't' ]; then
(
tar_file=${scrDir}/${dtg}/${dtgp1}_bkgrnd.tar
rm -f ${tar_file}
echo $tar_file
for element in `seq 1 ${ENSEMBLE_SIZE}`; do

  dataDir_mbr=`printf %0${digits}d/data ${element}`
echo $dataDir_mbr
  target_files="./${dataDir_mbr}/*_sfc_*_[1-9]a*_${dtg}_${cycle_fmt}_fcstfld  \
                ./${dataDir_mbr}/*_msl_*_[1-9]a*_${dtg}_${cycle_fmt}_fcstfld  \
			    ./${dataDir_mbr}/*_zht_*_${dtg}_${cycle_fmt}_fcstfld  \
			    ./${dataDir_mbr}/uuwind_sig_*_[1-9]a*_${dtg}_${cycle_fmt}_fcstfld  \
			    ./${dataDir_mbr}/vvwind_sig_*_[1-9]a*_${dtg}_${cycle_fmt}_fcstfld  \
			    ./${dataDir_mbr}/wwwind_sig_*_[1-9]a*_${dtg}_${cycle_fmt}_fcstfld  \
			    ./${dataDir_mbr}/pottmp_sig_*_[1-9]a*_${dtg}_${cycle_fmt}_fcstfld  \
			    ./${dataDir_mbr}/perprs_sig_*_[1-9]a*_${dtg}_${cycle_fmt}_fcstfld  \
			    ./${dataDir_mbr}/wvapor_sig_*_[1-9]a*_${dtg}_${cycle_fmt}_fcstfld  \
			    ./${dataDir_mbr}/turbke_sig_*_[1-9]a*_${dtg}_${cycle_fmt}_fcstfld  \
			    ./${dataDir_mbr}/???mix_sig_*_[1-9]a*_${dtg}_${cycle_fmt}_fcstfld  \
			    ./${dataDir_mbr}/datahd_*_1a*${dtg}_*  \
			    ./${dataDir_mbr}/terrht_*_[1-9]a*_${dtg}_*"
  tar rf ${tar_file} ${target_files}

  target_files="./${dataDir_mbr}/??wind_pre_*_[1-9]a*_${dtg}_${cycle_fmt}_fcstfld  \
                ./${dataDir_mbr}/geopht_pre_*_[1-9]a*_${dtg}_${cycle_fmt}_fcstfld  \
                ./${dataDir_mbr}/airtmp_pre_*_[1-9]a*_${dtg}_${cycle_fmt}_fcstfld  \
                ./${dataDir_mbr}/vpress_pre_*_[1-9]a*_${dtg}_${cycle_fmt}_fcstfld  \
                ./${dataDir_mbr}/dwptdp_pre_*_[1-9]a*_${dtg}_${cycle_fmt}_fcstfld"
  tar rf ${tar_file} ${target_files}

  target_files="./${dataDir_mbr}/*_sfc_*_[1-9]o*_${dtgm1_ocn}_${cycle_zero}_*"
  tar rf ${tar_file} ${target_files}

done
target_files="./${dtg}/prior_inflate_restart \
              ./input.nml \
              ./namelist \
  			  ./state.vars \
			  ./convert.vars "
tar rf ${tar_file} ${target_files}

archive put -C ${DART_ARC} ${tar_file}
rm -f ${tar_file}
) < /dev/null &
fi

#-----------------------------------------------------------------------------------
  # INIT FILES (FILES TO INTIALIZE THE MODEL AFTER DA)
#-----------------------------------------------------------------------------------
if [ $arch_init = 't' ]; then
(

tar_file=${scrDir}/${dtg}/${dtg}_init.tar
rm -f ${tar_file}
echo $tar_file
for element in `seq 1 ${ENSEMBLE_SIZE}`; do
  dataDir_mbr=`printf %0${digits}d/data ${element}`
  echo $dataDir_mbr

  target_files="./${dataDir_mbr}/*_sfc_*_[1-9]a*_${dtg}_00000000_fcstfld 	
                ./${dataDir_mbr}/analyses/*_msl_*_[1-9]a*_${dtg}_00000000_analfld
		     	./${dataDir_mbr}/analyses/*_zht_*_[1-9]a*_${dtg}_00000000_fcstfld
			    ./${dataDir_mbr}/analyses/*_sig_*_[1-9]a*_${dtg}_00000000_analfld
			    ./${dataDir_mbr}/datahd_*_${dtg}_*"
  tar rf ${tar_file} ${target_files}

  target_files="./${dataDir_mbr}/*_${dtg}_*_bndyfld"
  tar rf ${tar_file} ${target_files}

  target_files="./${dataDir_mbr}/*_sfc_*_[1-9]o*_${dtg}_${cycle_zero}_*"
  tar rf ${tar_file} ${target_files}
done
target_files="./input.nml \
              ./namelist \
  			  ./state.vars \
			  ./convert.vars "
tar rf ${tar_file} ${target_files}

archive put -C ${DART_ARC} ${tar_file}
rm -f ${tar_file}
) < /dev/null &

fi

#-----------------------------------------------------------------------------------
# SFC, ZHT, and MSL FILES
#-----------------------------------------------------------------------------------
if [ $arch_fcst = 't' ]; then
(

tar_file=${scrDir}/${dtg}/${dtg}_sfc.tar
rm -f ${tar_file}
echo $tar_file
for element in `seq 1 ${ENSEMBLE_SIZE}`; do
  dataDir_mbr=`printf %0${digits}d/data ${element}`
  echo $dataDir_mbr
  target_files="./${dataDir_mbr}/*_sfc_*_[1-9]a*_${dtg}_*_fcstfld \
                ./${dataDir_mbr}/*_zht_*_[1-9]a*_${dtg}_*_fcstfld \
                ./${dataDir_mbr}/*_msl_*_[1-9]a*_${dtg}_*_fcstfld \
				./${dataDir_mbr}/datahd_*_${dtg}_*"
  tar rvf ${tar_file} ${target_files}
done

archive put -C ${DART_ARC} ${tar_file}
rm -f ${tar_file}
)  < /dev/null &

#-----------------------------------------------------------------------------------
# SIGMA FILES
#-----------------------------------------------------------------------------------
(
tar_file=${scrDir}/${dtg}/${dtg}_sgl.tar
rm -f ${tar_file}
for element in `seq 1 ${ENSEMBLE_SIZE}`; do
  dataDir_mbr=`printf %0${digits}d/data ${element}`
  target_files="./${dataDir_mbr}/*_sig_*_${dtg}_*_fcstfld ./${dataDir_mbr}/datahd_*_${dtg}_* ./${dataDir_mbr}/terrht_*_${dtg}_*"
  tar rvf ${tar_file} ${target_files}
done

archive put -C ${DART_ARC} ${tar_file}
rm -f ${tar_file}
)  < /dev/null &

#-----------------------------------------------------------------------------------
# PRESSURE LEVEL FILES
#-----------------------------------------------------------------------------------
(
tar_file=${scrDir}/${dtg}/${dtg}_pre.tar
rm -f ${tar_file}
for element in `seq 1 ${ENSEMBLE_SIZE}`; do
  dataDir_mbr=`printf %0${digits}d/data ${element}`
  target_files="./${dataDir_mbr}/*_pre_*_${dtg}_*_fcstfld ./${dataDir_mbr}/datahd_*_${dtg}_* ./${dataDir_mbr}/terrht_*_${dtg}_*" 
  tar rvf ${tar_file} ${target_files}
done

archive put -C ${DART_ARC} ${tar_file}
rm -f ${tar_file}
)  < /dev/null &

fi

wait

if [ ${dtgp1} -le ${dtg_end} ]; then
${scriptDir}/archive_coamps_ens.sh -d ${dtgp1} -e ${dtg_end} -i ${icycle} \
		-n ${arch_init} -f ${arch_fcst} -r ${arch_restart} -b {arch_bkgrnd}

fi

exit 0

EOF

chmod u+x ${stageOUT}

# Submit archive job.
TRNS_JOB_ID=`${BATCH} ${depend} ${stageOUT}`
rm -f ${stageOUT}

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

