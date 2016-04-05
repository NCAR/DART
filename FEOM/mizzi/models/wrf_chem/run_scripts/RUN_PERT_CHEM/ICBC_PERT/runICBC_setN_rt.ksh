#!/bin/ksh -x
unalias ls
export YEAR=${YYYY}
export MONTH=${MM}
export DAY=${DD}
export HOUR=${HH}
export IC_METDIR=${WRFCHEM_MET_IC_DIR}
export BC_METDIR=${WRFCHEM_MET_BC_DIR}
export CHEMDIR=${RUN_DIR}/${YEAR}${MONTH}${DAY}${HOUR}/wrfchem_chem_icbc
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
   cp ${IC_METDIR}/${WRFINP} ./.
   cp ${BC_METDIR}/${WRFBDY} ./.
   ls set${i}
   rm -f mozbc.ic.inp.set${i}
   cat mozbc.ic.inp set${i} > mozbc.ic.inp.set${i}
   rm -f mozbc.bc.inp.set${i}
   cat mozbc.bc.inp set${i} > mozbc.bc.inp.set${i}
   echo  'run mozbc'
   ./run_mozbc_rt.csh type=ic mozbc_inp=mozbc.ic.inp.set${i} ens=${IENS}
   ./run_mozbc_rt.csh type=bc mozbc_inp=mozbc.bc.inp.set${i} ens=${IENS}
   echo 'ok'
#   cp ${WRFINP} ${CHEMDIR}/.
#   cp ${WRFBDY} ${CHEMDIR}/.
   echo 'put files'
   echo 'OK'
done
#################parent files################
export METDIR=${REAL_DIR}
export CHEMDIR=${RUN_DIR}/${YEAR}${MONTH}${DAY}${HOUR}/wrfchem_chem_icbc
export WRFINP=wrfinput_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00
export WRFBDY=wrfbdy_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00
cp ${METDIR}/${WRFINP} ./.
cp ${METDIR}/${WRFBDY} ./.
mv ${WRFINP} wrfinput_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000
mv ${WRFBDY} wrfbdy_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000
rm -f mozbc.ic.inp.set00
cat mozbc.ic.inp set00 > mozbc.ic.inp.set00
rm -f mozbc.bc.inp.set00
cat mozbc.bc.inp set00 > mozbc.bc.inp.set00
./run_mozbc_rt.csh type=ic mozbc_inp=mozbc.ic.inp.set00 ens=000
./run_mozbc_rt.csh type=bc mozbc_inp=mozbc.bc.inp.set00 ens=000
mv wrfinput_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000 ${WRFINP}
mv wrfbdy_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000 ${WRFBDY}
#cp ${WRFINP} ${CHEMDIR}/.
#cp ${WRFBDY} ${CHEMDIR}/.
################parent files#################
