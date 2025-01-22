#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# datea, emember, paramfile are command-line arguments - OR -
# are set by a string editor (sed) command.

set datea     = ${1}
set emember   = ${2}
set paramfile = ${3}

source $paramfile

set domains = $NUM_DOMAINS

set start_time = `date +%s`
echo "host is " `hostname`
echo "assim_advance.csh is running in `pwd`"

cd ${RUN_DIR}

set gdate = (`echo $datea 0 -g | ${RUN_DIR}/advance_time`)

if ( $ASSIM_INT_MINUTES <= 0 ) then
  set gdatef = (`echo $datea $ASSIM_INT_HOURS -g | ${RUN_DIR}/advance_time`)
else
  set gdatef = (`echo $datea ${ASSIM_INT_MINUTES}m -g | ${RUN_DIR}/advance_time`)
endif

set yyyy  = `echo $datea | cut -b1-4`
set mm    = `echo $datea | cut -b5-6`
set dd    = `echo $datea | cut -b7-8`
set hh    = `echo $datea | cut -b9-10`
set nn    = "00"
set ss    = "00"

#  copy files to appropriate location
echo $start_time >! ${RUN_DIR}/start_member_${emember}

# if ( -d ${RUN_DIR}/advance_temp${emember} )  ${REMOVE} ${RUN_DIR}/advance_temp${emember}
#mkdir -p ${RUN_DIR}/advance_temp${emember}

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
          $domains
mpirun.lsf ./wrf.exe
EOF

else if ( $SUPER_PLATFORM == 'derecho' ) then

# module load openmpi
cat >! $RUN_DIR/advance_temp${emember}/wrf.info << EOF
${gdatef[2]}  ${gdatef[1]}
${gdate[2]}   ${gdate[1]}
$yyyy $mm $dd $hh $nn $ss
           $domains
 mpiexec -n 128 -ppn 128  ./wrf.exe
EOF

endif

cd $RUN_DIR

echo $emember                      >! ${RUN_DIR}/filter_control${icnum}
echo filter_restart_d01.${icnum}   >> ${RUN_DIR}/filter_control${icnum}
echo prior_d01.${icnum}            >> ${RUN_DIR}/filter_control${icnum}

#  integrate the model forward in time
${RUN_DIR}/new_advance_model.csh ${emember} $domains filter_control${icnum} $paramfile
${REMOVE} ${RUN_DIR}/filter_control${icnum}

set end_time   = `date  +%s`
@ length_time  = $end_time - $start_time
echo "duration = $length_time"

exit 0

