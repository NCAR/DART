#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$
#
# script for copying 1 day/obs_seq of output diagnostics to mass store.

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/local/dcs/lib
echo $LD_LIBRARY_PATH

# set  echo verbose

set compress = false

#   echo "usage; from case/experiment/obs_seq directory; "
#   echo "       diag2ms compress(optnl)"
#   echo "       optionally compresses diagnostic files."
#   exit

if ($#argv == 1) then
   set compress = true
endif 

set direct = `pwd`
set obs_seq = $direct:t

cd ..
set direct = `pwd`
set case = $direct:t

cd ..
set direct = `pwd`
set exp_dir = $direct:t

cd $case/${obs_seq}
set ms_dir = /RAEDER/DAI/${exp_dir}/$case/${obs_seq}
echo files will be written to ${ms_dir}/diagnostics.tar

tar cv -f diagnostics.tar --exclude saved [^CD]* >& saved
if ($compress == true) then
   gzip diagnostics.tar
   msrcp -pe 365 diagnostics.tar.gz mss:${ms_dir}/diagnostics.tar.gz 
else
   msrcp -pe 365 diagnostics.tar mss:${ms_dir}/diagnostics.tar
endif

msls -l "${ms_dir}/diagnostics.tar*"
if ($status == 0) then
   rm diagnostics.tar.gz 
   echo ${ms_dir} was copied to mass store >> saved
else
   echo ${ms_dir} was NOT copied to mass store
endif

echo "msrcp done; check files on mass store, and then delete local copies"

exit
