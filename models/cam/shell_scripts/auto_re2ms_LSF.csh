#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source: /home/thoar/CVS.REPOS/DART/models/cam/shell_scripts/auto_re2ms_LSF.csh,v $
# $Name:  $

#
# Saves restart files to MSS. CAM/CLM, and filter. This is the batch driver
# for auto_re2ms.csh ...
#

#### LSF options for BSUB
### -J      job name    (master script job.csh presumes filter.xxxx.log)
### -o      output listing filename 
### -P      account number
### -q      queue
### -n      number of tasks (processors)
### -x      exclusive use of node
### -R "span[ptile=(num procs you want on each node)]"
#
#BSUB -J restart2ms
#BSUB -o restart2ms.%J.log
#BSUB -P 86850054
#BSUB -q standby
#BSUB -W 1:00
#BSUB -n 1
#xxxx -x
#xxxx -R "span[ptile=1]"

date
echo 'auto_re2ms starts in'
pwd
cd $LS_SUBCWD

set num_ens = 1
while (-s CAM/caminput_${num_ens}.nc)
   @ num_ens++
end
@ num_ens = $num_ens - 1
echo "num_ens = $num_ens"

# leaving off a factor of 100 here to handle limit on size of numbers (2^31 -1)

if (-e DART/filter_ic.0001) then
   set list = `du -bc DART/filter_ic.0001 CAM/caminput_1.nc CLM/clminput_1.nc`
   @ size_element = $list[7] / 100
else
   set list = `du -bc CAM/caminput_1.nc CLM/clminput_1.nc`
   @ size_element = $list[5] / 100
endif
echo "size_element = $size_element"
# This assumes compression factor of .8, from CAM and CLM only
# (2^31 - 1) / .8 = 2680000000
# orig; 
@ test_element =    26500000 / $size_element
# gz 
# @ test_element =    21000000 / $size_element
set div = 1
set size_file = ${num_ens}
echo div, size_file are $div $size_file     
while ( $size_file > $test_element)
    @ div = $div * 2
    @ size_file = $size_file / 2
    echo div, size_file are $div $size_file 
end
if ($div != 0) then
   @ num_per_batch = ${num_ens} / $div
else
   set num_per_batch = 666
endif
echo "For MS backup num_per_batch = $num_per_batch"        
echo "   for size_element(uncomp) = $size_element * 100"   
echo "                for num_ens = $num_ens"              


../../auto_re2ms.csh $num_ens $num_per_batch comp

echo "finished with auto_re2ms at " `date` 

exit
