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
# script for copying multiple day/obs_seq of diagnostics files to mass store.
# For lack of better bookkeeping, the diag2ms.csh script needs to be copied
# to $HOME ... at some point, this may be streamlined by integrating diag2ms.csh
# into this script.  TJH 9/30/2005

#### LSF options for BSUB
### -J      job name    (master script job.csh presumes filter.xxxx.log)
### -o      output listing filename
### -P      account number
### -q      queue
### -n      number of tasks (processors)
### -x      exclusive use of node
### -R "span[ptile=(num procs you want on each node)]"
#
#BSUB -J diags2ms
#BSUB -o diags2ms.%J.log
###BSUB -P 86850054
#BSUB -P 39510050
#BSUB -q standby
#BSUB -n 1
#xxxx -x
#xxxx -R "span[ptile=1]"

set compress = true
set obs_root = 01_
set obs_seq_first = 1
set obs_seq_last = 58


date
echo 'diags2ms starts in'
pwd
cd $LS_SUBCWD


setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/local/dcs/lib
echo $LD_LIBRARY_PATH

# set  echo verbose

# if ($#argv < 2) then
#    echo "usage; from case/experiment directory; "
#    echo "       diags2ms obs_seq_first obs_seq_last compress(optnl)"
#    echo "       Optionally Compresses  tar files."
#    echo "       ASSUMES obs_root of 01_; change for other months/roots"
#    exit
# else if ($#argv == 3) then
#    set compress = true
# endif 

set obs_seq = $obs_seq_first
while ($obs_seq <= $obs_seq_last)
   set obs = ${obs_root}
   if ($obs_seq < 10) set obs = ${obs_root}0
   set obs = $obs$obs_seq
   cd $obs
   if ($compress == true) then
      ~/diag2ms.csh $compress
   else
      ~/diag2ms.csh
   endif
   cd ..
   @ obs_seq++
end

exit
