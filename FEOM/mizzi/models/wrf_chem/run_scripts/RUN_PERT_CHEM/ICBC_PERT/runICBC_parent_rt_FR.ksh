#!/bin/ksh -x
export YEAR=${YYYY}
export MONTH=${MM}
export DAY=${DD}
export HOUR=${HH}
export METDIR=${REAL_DIR}
export CHEMDIR=${RUN_DIR}/${DATE}/wrfchem_chem_icbc
export WRFINP_FR=wrfinput_d02_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00
#export WRFBDY_FR=wrfbdy_d02_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00
cp ${METDIR}/${WRFINP_FR} ./.
#cp ${METDIR}/${WRFBDY_FR} ./.
mv ${WRFINP_FR} wrfinput_d02_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000
#mv ${WRFBDY_FR} wrfbdy_d02_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000
rm -f mozbc.ic.inp.set00
cat mozbc.ic.inp set00 > mozbc.ic.inp.set00
#rm -f mozbc.bc.inp.set00
#cat mozbc.bc.inp set00 > mozbc.bc.inp.set00
./run_mozbc_rt_FR.csh type=ic mozbc_inp=mozbc.ic.inp.set00 ens=000
#./run_mozbc_rt_FR.csh type=bc mozbc_inp=mozbc.bc.inp.set00 ens=000
mv wrfinput_d02_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000 ${WRFINP_FR}
#mv wrfbdy_d02_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000 ${WRFBDY_FR}

