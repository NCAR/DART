#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$
#
# script for copying diagnostics files to mass store.

#### LSF options for BSUB
### -J      job name    (master script job.csh presumes filter.xxxx.log)
### -o      output listing filename
### -P      account number
### -q # Queue name    regular   economy  standby     long   
#        proclim            32        16        8        2      
#        timelim         6 hrs        18       48   5 days  
#        # jobs/person       2         -        -        2
### -n      number of tasks (processors)
### -x      exclusive use of node
### -R "span[ptile=(num procs you want on each node)]"
#
#BSUB -J auto_diag2ms
#BSUB -o auto_diag2ms.%J.log
#BSUB -e auto_diag2ms.%J.err
#BSUB -P account_number
#BSUB -q premium
#BSUB -W 3:00
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#xxxx -x

# set  echo verbose

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/local/dcs/lib

set compress = true
set proj_num = 12345678

set diag_name = diagnostics.tar
set saved = saved_diagnostics
set write_pass = da$$

if ($?LS_SUBCWD) then
   cd $LS_SUBCWD
endif

echo 'auto_diag2ms_LSF starts in' >! $saved
pwd                               >> $saved
date                              >> $saved

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

# IBM tar requires 1 entry/line in list of things to exclude from tar
echo DART     >! tar_excl_list
echo CAM      >> tar_excl_list
echo CLM      >> tar_excl_list
# batch* are the files into which DART,CAM,CLM are archived by auto_re2ms_LSF.csh,
# which is run at the same time as this script.
echo 'batch*' >> tar_excl_list
echo $saved   >> tar_excl_list

tar -c -v -f $diag_name -X tar_excl_list * >>& $saved
rm tar_excl_list

if ($compress == true) then
   gzip $diag_name
   set diag_name = ${diag_name}.gz
endif

echo "files will be written to ${ms_dir}/${diag_name}" >> $saved
echo "with write password $write_pass" >> $saved

msrcp -pe 1000 -pr $proj_num -wpwd $write_pass -comment "write password $write_pass" \
      ${diag_name} mss:${ms_dir}/${diag_name} >>& $saved


# check to see if it's okay to remove diagnostics.tar and P*.nc (leave obs_seq.final)
# Better check would be to add before the gzip
#   tar -t -f diagnostics.tar > table
#   foreach name (Posterior_Diag.nc Prior_Diag.nc obs_seq.final)
#      grep $name table
#      if ($status != 0) then
#         set tarok = false
#      endif
#   end
#   
set list = `ls -l $diag_name`
set local_size = $list[5]
set list = `msls -l ${ms_dir}/${diag_name}`
set ms_size = $list[5]
echo " ${diag_name} local_size = $local_size, ms_size = $ms_size" >> $saved

if ($local_size == $ms_size) then
   echo "Archived files with write password $write_pass" >> $saved
   echo "msrcp of $ms_dir/$obs_seq succeeded; REMOVING $diag_name and P*.nc " >> $saved
   rm $diag_name P*.nc
else
   echo "msrcp of ${ms_dir}/$obs_seq  failed; " >> $saved
   echo "NOT removing $diag_name and P*.nc"      >> $saved
endif

chmod 444 $saved

exit
