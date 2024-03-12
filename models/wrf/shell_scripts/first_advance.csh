#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

set datea     = ${1}
set emember   = ${2}
set paramfile = ${3}

source $paramfile

set start_time = `date +%s`
echo "host is " `hostname`

set domains = $NUM_DOMAINS

cd ${RUN_DIR}

set gdate = (`echo $datea 0 -g | ${RUN_DIR}/advance_time`)
set gdatef = (`echo $datea $ASSIM_INT_HOURS -g | ${RUN_DIR}/advance_time`)
set yyyy  = `echo $datea | cut -b1-4`
set mm    = `echo $datea | cut -b5-6`
set dd    = `echo $datea | cut -b7-8`
set hh    = `echo $datea | cut -b9-10`
set nn    = "00"
set ss    = "00"

echo $start_time >! ${RUN_DIR}/start_member_${emember}

# go into member directory and generate the needed wrf.info file
cd $RUN_DIR/advance_temp${emember}

set icnum = `echo $emember + 10000 | bc | cut -b2-5`
if ( -e $RUN_DIR/advance_temp${emember}/wrf.info ) then
   ${REMOVE} $RUN_DIR/advance_temp${emember}/wrf.info
endif

touch wrf.info

if ( $SUPER_PLATFORM == 'LSF queuing system' ) then

   cat >! $RUN_DIR/advance_temp${emember}/wrf.info << EOF
 ${gdatef[2]}  ${gdatef[1]}
 ${gdate[2]}   ${gdate[1]}
 $yyyy $mm $dd $hh $nn $ss
           1
 mpirun.lsf ./wrf.exe
EOF

else if ( $SUPER_PLATFORM == 'derecho' ) then

   setenv MPI_SHEPHERD false

   cat >! $RUN_DIR/advance_temp${emember}/wrf.info << EOF
 ${gdatef[2]}  ${gdatef[1]}
 ${gdate[2]}   ${gdate[1]}
 $yyyy $mm $dd $hh $nn $ss
           $domains
 mpiexec -n 128 -ppn 128 ./wrf.exe
EOF

endif

cd $RUN_DIR

echo $emember                      >! ${RUN_DIR}/filter_control${icnum}
echo filter_restart_d01.${icnum}   >> ${RUN_DIR}/filter_control${icnum}
echo prior_d01.${icnum}            >> ${RUN_DIR}/filter_control${icnum}

#  integrate the model forward in time
${RUN_DIR}/new_advance_model.csh ${emember} 1 filter_control${icnum} $paramfile
${REMOVE} ${RUN_DIR}/filter_control${icnum}

# move the output to the appropriate directory
mkdir -p ${OUTPUT_DIR}/${datea}/PRIORS
mv $RUN_DIR/prior_d01.${icnum} ${OUTPUT_DIR}/${datea}/PRIORS/prior_d01.${icnum}

set end_time   = `date  +%s`
@ length_time  = $end_time - $start_time
echo "duration = $length_time"

exit 0

