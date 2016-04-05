#!/bin/ksh -x
unalias ls
export DIR=${WPB_RC_CHEM_DIR}
export HSIDIR=${HSI_WPB_RC_CHEM_DIR}
if ${USE_HSI}; then
   hsi "mkdir ${HSIDIR}"
else
   mkdir ${DIR}
fi
#
export YEAR=${YYYY}
export MONTH=${MM}
export DAY=${DD}
export HOUR=${HH}
if ${USE_HSI}; then
   export METDIR=${HSI_WPB_RC_DIR}/${YEAR}${MONTH}${DAY}${HOUR}
else
   export METDIR=${WPB_RC_DIR}/${YEAR}${MONTH}${DAY}${HOUR}
fi
export CHEMDIR=${WPB_RC_CHEM_DIR}/${YEAR}${MONTH}${DAY}${HOUR}
mkdir -p ${CHEMDIR}
for ((i = 1; i <= ${NUM_MEMBERS}; i += 1))
do
   if [ "$i" -lt "10"  ]; then 
      export IENS="00"${i}
   fi
   if [ "$i" -lt "100"  ]; then
      if [ "$i" -ge "10"  ]; then
         export IENS="0"${i}
      fi
   fi
   if [ "$i" -lt "1000"  ]; then
      if [ "$i" -ge "100"  ]; then
         export IENS=${i}
      fi
   fi
   export WRFINP=wrfinput_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e${IENS}
   export WRFBDY=wrfbdy_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e${IENS}
   echo 'get inp'
   if ${USE_HSI}; then
      hsi get ${METDIR}/${WRFINP} 
      hsi get ${METDIR}/${WRFBDY}
   else
      cp ${METDIR}/${WRFINP} ./.
      cp ${METDIR}/${WRFBDY} ./.
   fi
   ls set${i}
   rm -f mozbc.ic.inp.set${i}
   cat mozbc.ic.inp set${i} > mozbc.ic.inp.set${i}
   rm -f mozbc.bc.inp.set${i}
   cat mozbc.bc.inp set${i} > mozbc.bc.inp.set${i}
   echo  'run mozbc'
   ./run_mozbc.csh type=ic mozbc_inp=mozbc.ic.inp.set${i} ens=${IENS}
   ./run_mozbc.csh type=bc mozbc_inp=mozbc.bc.inp.set${i} ens=${IENS}
   echo 'ok'
   cp ${WRFINP} ${CHEMDIR}/.
   cp ${WRFBDY} ${CHEMDIR}/.
   echo 'put files'
#   hsi "mkdir -p ${HSIDIR}/${YEAR}${MONTH}${DAY}${HOUR}; cd ${HSIDIR}/${YEAR}${MONTH}${DAY}${HOUR}; put ${WRFINP}; put ${WRFBDY}"
   echo 'OK'
done
#################parent files################
export DIR=${RC_CHEM_DIR}
export HSIDIR=${HSI_RC_CHEM_DIR}
if ${USE_HSI}; then
   export METDIR=${HSI_RC_DIR}/${YEAR}${MONTH}${DAY}${HOUR}
else
   export METDIR=${RC_DIR}/${YEAR}${MONTH}${DAY}${HOUR}
fi
export CHEMDIR=${RC_CHEM_DIR}/${YEAR}${MONTH}${DAY}${HOUR}
mkdir -p ${CHEMDIR}
export WRFINP=wrfinput_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00
export WRFBDY=wrfbdy_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00
if ${USE_HSI}; then
   hsi get ${METDIR}/${WRFINP}
   hsi get ${METDIR}/${WRFBDY}
else
   cp ${METDIR}/${WRFINP} ./.
   cp ${METDIR}/${WRFBDY} ./.
fi
mv ${WRFINP} wrfinput_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000
mv ${WRFBDY} wrfbdy_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000
rm -f mozbc.ic.inp.set00
cat mozbc.ic.inp set00 > mozbc.ic.inp.set00
rm -f mozbc.bc.inp.set00
cat mozbc.bc.inp set00 > mozbc.bc.inp.set00
./run_mozbc.csh type=ic mozbc_inp=mozbc.ic.inp.set00 ens=000
./run_mozbc.csh type=bc mozbc_inp=mozbc.bc.inp.set00 ens=000
mv wrfinput_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000 ${WRFINP}
mv wrfbdy_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000 ${WRFBDY}
cp ${WRFINP} ${CHEMDIR}/.
cp ${WRFBDY} ${CHEMDIR}/.
#hsi "cd ${HSIDIR}/${YEAR}${MONTH}${DAY}${HOUR}; put ${WRFINP}; put ${WRFBDY}"
################parent files#################
