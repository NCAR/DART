#!/bin/ksh -x
unalias ls
export DIR=${RC_CHEM_DIR}
export HSIDIR=${HSI_RC_CHEM_DIR}
if ${USE_HSI}; then
   hsi "mkdir -p ${HSIDIR}"
else
   mkdir ${DIR}
fi
#
export YEAR=${YYYY}
export MONTH=${MM}
export DAY=${DD}
export HOUR=${HH}
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
mv ${WRFBDY} wrfbdy_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000
mv ${WRFINP} wrfinput_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000
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
rm -f ${WRFBDY} ${WRFINP}
#echo "cd ${HSIDIR} ;  put -R ${CHEMDIR}"
#hsi "cd ${HSIDIR} ;  put -R ${CHEMDIR}"
#rm -rf ${CHEMDIR}
