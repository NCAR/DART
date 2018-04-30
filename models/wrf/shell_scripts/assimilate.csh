#!/bin/csh

set datea     = ${1}
set paramfile = ${3}

source $paramfile

set start_time = `date +%s`
echo "host is " `hostname`

cd $RUN_DIR
echo $start_time >& ${RUN_DIR}/filter_started

#  run data assimilation system
if ( $SUPER_PLATFORM == 'yellowstone' ) then
## Yellowstone
 setenv TARGET_CPU_LIST -1
 setenv FORT_BUFFERED true
 mpirun.lsf ./filter
else if ( $SUPER_PLATFORM == 'cheyenne' ) then
## Cheyenne
 setenv TMPDIR  /dev/shm
 setenv MPI_SHEPHERD true
 limit stacksize unlimited
 module load mpt
 mpiexec_mpt dplace -s 1 ./filter
# module load openmpi
# module load peak_memusage
# mpirun ./filter
endif

if ( -e ${RUN_DIR}/obs_seq.final )  touch ${RUN_DIR}/filter_done

set end_time = `date  +%s`
@ length_time = $end_time - $start_time
echo "duration = $length_time"

exit 0
