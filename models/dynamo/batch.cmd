#!/bin/ksh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
#BSUB -J dyn_dart_test
#BSUB -n 20
#BSUB -R "span[ptile=32]"
#BSUB -o dyn_dart_out.%J
#BSUB -e dyn_dart_err.%J
#BSUB -W 1:05
#BSUB -P 22104000
#BSUB -q share
###BSUB -N

LSB_MCPU_HOSTS_STORE=${LSB_MCPU_HOSTS}

export MP_LABELIO=yes

start=1  #-- Starting iteration
stop=200 #-- Stop at

diff=`expr ${stop} - ${start} + 1`
store_dir=${start}_${stop}

#-- Generate input file create_fixed_network_seq and run that
#-- usually done by create_obs_sequence
cat << EOF > inp.obs
set_def.out
1
${diff}
${start} 0
1 0
obs_seq.in
EOF
./create_fixed_network_seq < inp.obs

sh clean #-- Remove the debri from previous run
#-- This is just one task parallel job
export LSB_MCPU_HOSTS=$(echo ${LSB_MCPU_HOSTS_STORE}|awk '{printf "%s 1", $1}')
rm -f env.lsf
env | egrep  '^(XL|L)S(F|B)' > env.lsf #-- Store all the LSF env for mpirun.lsf 
cp env.lsf env.lsf.tmp                 #-- to be spawned by advance_model.csh
printf "%s " export >> env.lsf
awk -F= '{printf "%s ", $1}' env.lsf.tmp >> env.lsf
printf "\n" >> env.lsf

env PERFECT=1 ./perfect_model_obs 
if [ $? != 0 ] ; then
   exit
fi

#-- 20 task parallel spawn by advance_model.csh
sh clean #-- Remove the debri from perfect_model run
if [ ! -z $prev_dir ] ; then
  cp $prev_dir/inr* $prev_dir/filter_restart .
  chmod 600 inr*
  cp filter_restart perfect_ics
fi

export LSB_MCPU_HOSTS=${LSB_MCPU_HOSTS_STORE}
rm -f env.lsf
env | egrep  '^(XL|L)S(F|B)' > env.lsf #-- Store all the LSF env for mpirun.lsf 
cp env.lsf env.lsf.tmp                 #-- to be spawned by advance_model.csh
printf "%s " export >> env.lsf
awk -F= '{printf "%s ", $1}' env.lsf.tmp >> env.lsf
printf "\n" >> env.lsf

env PERFECT=0 ./filter
if [ $? != 0 ] ; then
   exit
fi

#-- Store everything in a separate directory
mkdir $store_dir
cp inr* perfect_ics filter_restart *.nc obs_seq.* $store_dir
chmod 400 $store_dir/*

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
