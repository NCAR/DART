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
#BSUB -P xxxxxxxx
#BSUB -q share
#BSUB -W 2:00
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#xxxx -x

# set  echo verbose

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/local/dcs/lib

set compress = true
set proj_num = 12345678

set diag_name = diagnostics.tar
set saved = saved_diagnostics
set write_pass = $$

if ($?LS_SUBCWD) then
   cd $LS_SUBCWD
endif

touch $saved
echo '------------------------------------------------------------' >> $saved
echo 'auto_diag2ms_LSF starts in' >> $saved
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
echo DART                >! tar_excl_list
echo CAM                 >> tar_excl_list
echo CLM                 >> tar_excl_list
echo ICE                 >> tar_excl_list
# batch* are the files into which DART,CAM,CLM are archived by auto_re2ms_LSF.csh,
# which is run at the same time as this script.
echo 'batch*'            >> tar_excl_list
echo $saved              >> tar_excl_list

## Added to make mean easily accessible in the form of a CAM initial file
echo 'cam_analyses.tar'  >> tar_excl_list
echo 'H[0-9]*'           >> tar_excl_list

#-----------------------------
# Stuff the Posterior mean fields into CAM initial files.
# Tar together with a CLM initial file.
# Arguments to analyses2initial.csh are
#   set ms_file   = $1   script searches for local Posterior_Diag.nc first, so give a dummy MS name.
#   set local_dir = $2
#   set kind      = $3
#   set dim       = $4
#   set element1  = $5
#   set yrmoday   = $6   set to $obs_seq instead of yyyymmdd, since obs_seq is easily available

# Save out history files from H* directories before calculating CAM/CLM analyses and deleting H*
# compressing saves a factor of 12.
   ls H[012]*/*.h0.* >& /dev/null
   if ($status == 0) then
      gzip H[012]*/*.h0.*
      tar -c -f H_all.h0.gz.tar H[012]*/*.h0.*
      set ar_status = $status
      if ($ar_status == 0 && -e H_all.h0.gz.tar) then
         msrcp -pe 1000 -pr $proj_num -wpwd $write_pass -comment "write password $write_pass" \
               H_all.h0.gz.tar mss:${ms_dir}/H_all.h0.gz.tar                    >>& $saved  &
         set ar_status = $status
         if ($ar_status == 0) rm H[012]*/*.h0.* 
      endif
      if ($ar_status != 0) echo 'ARCHIVING of H[012]*.h0.gz.tar FAILED' >> & $saved
   else
       echo 'ARCHIVING of H[012]*.h0.gz.tar FAILED; no files available' >> & $saved
   endif 

if (-e ../../analyses2initial.csh) then
   # analyses2initial.csh needs CAM initial files to average and receive the analyses 
   # from Posterior_Diag.nc.
   # They should have been saved during the assimilation and be living in the 
   # exp/obs_####/H## directories.
   ls H*/caminput_1.nc > /dev/null
   set stat = $status
   if ($stat == 0) then
      ls H*/clminput_1.nc > /dev/null
      set stat = $status
   endif
   if ($stat != 0) then
      echo "H*/c[al]minput_* not available"                           >>& $saved
      echo "H*/c[al]minput_* not available" >&  ANALYSES_NOT_SAVED
   else
#   if ($status == 0) ls H*/clminput_1.nc > /dev/null
#   if ($status != 0) then
#      echo "H*/c[al]minput_* not available; exiting"                           >>& $saved
#      echo "H*/c[al]minput_* not available; exiting" >&  ANALYSES_NOT_SAVED
#      exit
#   endif
   
      set num_anal = `ls H[0-2]*/cam_init_*`
#      ls H[0-2]*/c*_init_* > /dev/null
#      if ($status != 0) then
      set tar_stat = 0
      if (! -e cam_analyses.tar) then
         if ($#num_anal < 4) then
            echo " "                                                                >>& $saved
            ../../analyses2initial.csh no_MS '.' Posterior copy 1 ${obs_seq}        >>& $saved
            set num_anal = `ls H[0-2]*/cam_init_*`
         endif
   
         tar -c -f cam_analyses.tar H[0-2]*/{c,ice}*_init_* 
         set tar_stat = $status
      endif
      if ($tar_stat == 0 )  \
         msrcp -pe 1000 -pr $proj_num -wpwd $write_pass -comment "write password $write_pass" \
               cam_analyses.tar mss:${ms_dir}/cam_analyses.tar                    >>& $saved  
      set list = `ls -l cam_analyses.tar`
      set local_size = $list[5]
      set list = `msls -l ${ms_dir}/cam_analyses.tar`
      set ms_size = $list[5]
      echo " cam_analyses.tar local_size = $local_size, ms_size = $ms_size"       >> $saved
   
      if ($local_size == $ms_size) then
         echo "Archived $ms_dir/cam_analyses.tar with write password $write_pass" >> $saved
         echo '    REMOVING H[0-9]* and cam_analyses.tar '                        >> $saved
         rm -rf H[0-9]*/[ci]* cam_analyses.tar
      else
         echo "msrcp of ${ms_dir}/cam_analyses.tar  failed; "                     >> $saved
         echo 'NOT removing H[0-9]* and cam_analyses.tar '                        >> $saved
      endif
   endif
else
   echo "NO analyses2initial.csh, so no CAM initial file format analyses created"
endif

#-----------------------------

if (! -e $diag_name && ! -e ${diag_name}.gz) tar -c -v -f $diag_name -X tar_excl_list * >>& $saved
rm tar_excl_list

if ($compress == true) then
   if (-e $diag_name) then
      gzip $diag_name                         >>& $saved
      set diag_name = ${diag_name}.gz
   else if (-e ${diag_name}.gz) then
      set diag_name = ${diag_name}.gz
   else
      echo "$diag_name does not exist at gzip" >> $saved
      exit
   endif
endif

echo "files will be written to ${ms_dir}/${diag_name}" >> $saved
echo "with write password $write_pass" >> $saved

msrcp -pe 1000 -pr $proj_num -wpwd $write_pass -comment "write password $write_pass" \
      ${diag_name} mss:${ms_dir}/${diag_name} >>& $saved

set list = `ls -l $diag_name`
set local_size = $list[5]
set list = `msls -l ${ms_dir}/${diag_name}`
set ms_size = $list[5]
echo " ${diag_name} local_size = $local_size, ms_size = $ms_size" >> $saved

if ($local_size == $ms_size) then
   echo "Archived files with write password $write_pass" >> $saved
   echo "msrcp of $ms_dir/$diag_name succeeded; REMOVING $diag_name and P*Diag.nc " >> $saved
   rm $diag_name P*Diag.nc
else
   echo "msrcp of ${ms_dir}/$obs_seq  failed; " >> $saved
   echo "NOT removing $diag_name and P*.nc"      >> $saved
endif

chmod 444 $saved

wait
if ($ar_status == 0) rm H_all.h0.gz.tar

exit
