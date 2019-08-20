#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

set datea     = ${1}
set paramfile = ${3}

source $paramfile

set start_time = `date +%s`
echo "host is " `hostname`

cd ${RUN_DIR}
echo $start_time >& ${RUN_DIR}/filter_started

#  run data assimilation system
if ( $SUPER_PLATFORM == 'yellowstone' ) then

   setenv TARGET_CPU_LIST -1
   setenv FORT_BUFFERED true
   mpirun.lsf ./filter || exit 1

else if ( $SUPER_PLATFORM == 'cheyenne' ) then

# TJH MPI_SHEPHERD may be a very bad thing 
# TJH setenv MPI_SHEPHERD true
# TJH module load mpt
   setenv TMPDIR  /dev/shm
   limit stacksize unlimited
   mpiexec_mpt dplace -s 1 ./filter || exit 1

endif

if ( -e ${RUN_DIR}/obs_seq.final )  touch ${RUN_DIR}/filter_done

set end_time = `date  +%s`
@ length_time = $end_time - $start_time
echo "duration = $length_time"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
