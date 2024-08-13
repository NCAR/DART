#!/bin/csh


 set paramfile = ${1}
 source $paramfile

 cd ${ICBC_DIR}
 mpiexec -n 128 -ppn 128 ${RUN_DIR}/WRF_RUN/real.exe

#if ( `grep "Successful completion of program real.exe" ./rsl.out.0000 | wc -l ` == 1  )  touch ${ICBC_DIR}/real_done

 touch ${ICBC_DIR}/real_done

exit 0

