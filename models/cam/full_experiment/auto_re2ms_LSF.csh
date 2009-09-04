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

# Saves restart files to MSS. CAM/CLM/ICE, and filter. This is the batch driver
# for auto_re2ms.csh ...

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
#BSUB -P [project number]
#BSUB -q share
#BSUB -W 2:00
#BSUB -n 1
#xxxx -x
#xxxx -R "span[ptile=1]"

date
echo 'auto_re2ms starts in'
pwd
cd $LS_SUBCWD

set comp = ''
set ret_period = 1000
set proj_num = #####
set write_pass = $$

set num_ens = 1
while (-e CAM/caminput_${num_ens}.nc && ! -z CAM/caminput_${num_ens}.nc)
   @ num_ens++
end
@ num_ens = $num_ens - 1
echo "num_ens = $num_ens"

# leaving off a factor of 100 here to handle limit on size of numbers (2^31 -1)

set mem_parts = 'CAM/caminput_1.nc CLM/clminput_1.nc'
set mem_entry = 5
if (-e DART/filter_ic.0001) then
   set mem_parts = "$mem_parts DART/filter_ic.0001 "
   @ mem_entry = $mem_entry + 2
endif
if (-e ICE/iceinput_1.tar) then
   set mem_parts = "$mem_parts ICE/iceinput_1.tar "
   @ mem_entry = $mem_entry + 2
endif
set list = `du -bc $mem_parts`
@ size_element = $list[$mem_entry] / 100

echo "size_element = $size_element"
# This assumes compression factor of .8, from CAM and CLM only
# (2^31 - 1) / .8 = 2680000000
@ test_element =    70000000 / $size_element
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


touch saved_restart
echo "------------------------------------------------" >> saved_restart
echo "with write password $write_pass" >> saved_restart

if ($num_per_batch > $num_ens) exit "restart2ms; num_per_batch > num_ens"

# Parse parts of the path name for construction of MSS name.
set direct = `pwd`
set obs_seq = $direct:t

cd ..
set direct = `pwd`
set exp_dir = $direct:t

cd ..
set direct = `pwd`
set case = $direct:t

cd ${exp_dir}/${obs_seq}

set ms_root = /RAEDER/DAI/$case/${exp_dir}/${obs_seq}/${num_ens}x${num_per_batch}
set ms_dir = mss:$ms_root
echo "files will be written to ${ms_root}/batch#" >> saved_restart
echo "files will be written to ${ms_root}/batch#"

# Figure out how many files to divide ensemble members among
@ nbatch = $num_ens / $num_per_batch
if ($num_ens % $num_per_batch != 0 ) @ nbatch++

set opts = "-pe $ret_period -pr $proj_num -wpwd $write_pass -comment "
set p1 = "write password $write_pass"
echo $opts $p1 >> saved_restart

set ok_to_remove = true

# Tally up and list the extra (non-ensemble) ic files
# Add more 'if' blocks as needed.
set DART_list = ' '
set num_files = 0
set do_filter = false
if (-e DART/prior_inf_ic)    then
   set DART_list = ($DART_list prior_inf_ic)
   @ num_files++
endif
if (-e DART/post_inf_ic)     then
   set DART_list = ($DART_list post_inf_ic)
   @ num_files++
endif
if (-e DART/filter_ic) then
   set DART_list = ($DART_list filter_ic)
   @ num_files++
else if (-e DART/filter_ic*0$num_ens || -e DART/filter_ic*.$num_ens) then
   set do_filter = true
else
   echo 'NOT ENOUGH ENSEMBLE MEMBERS IN .../DART' >> saved_restart
   exit
endif
echo "DART_list = ($DART_list) " >> saved_restart
echo "DART_list = ($DART_list) "

# Convert wordlist into a single character string (insert ,s)
set n = 1
if ($n <= $num_files) then
   set DART_files = $DART_list[1]
   while ($n < $num_files)
      @ n++
      set DART_files = "$DART_files,$DART_list[$n]"
   end
else
   set DART_files = 'none'
endif

if ($DART_files != 'none' && $do_filter == false) then
   tar -c -f ic_files.tar DART/{$DART_files}
   msrcp $opts "$p1" ic_files.tar ${ms_dir}/DART/ic_files.tar &
endif

# Compress here, instead of within tar; -z not available on some machines.
# It's not effective to compress filter_ic files (big and dense)
set fil = CAM/caminput_1.*
if ($comp != '' &&  $fil:e != 'gz') then
   gzip -r CAM
   gzip -r CLM
   echo "compressed CAM and CLM" >> saved_restart
endif

set batch = 1
while($batch <= $nbatch)

   @ base = (($batch - 1) * $num_per_batch) + 1
   @ member = $base
   if ($batch == 1 && (-e CAM/caminput_0.nc || -e CAM/caminput_0.nc.gz)) set member = 0
# BUT will there be a filter_ic.0000 ?
    echo base and member $base $member >> saved_restart
    echo base and member $base $member

# create the tar file using the first ensemble member of this batch
   set mem_parts = "CAM/caminput_$member.nc* CLM/clminput_$member.nc*"
   if ($do_filter == true)                         set mem_parts = "$mem_parts DART/filter_ic*[.0]$member"
   if (-e ICE/iceinput_${member}.tar)              set mem_parts = "$mem_parts ICE/iceinput_${member}.tar "
   if ($DART_files != 'none' && $member <= $base)  set mem_parts = "$mem_parts DART/{$DART_files}"
   tar -c -f batch${batch}${comp} $mem_parts
   echo "tar -c -f batch${batch}${comp} $mem_parts" 
   echo "tar -c -f batch${batch}${comp} $mem_parts" >> saved_restart

   # tack on additional ens members until this batch is complete
   set n = 2
   while ($n <= $num_per_batch)
     @ member++
      if ($member > $num_ens) then
         @ batch = $nbatch + 1
         @ n = $num_per_batch + 1
      else
         set mem_parts = "CAM/caminput_$member.nc* CLM/clminput_$member.nc*"
         if ($do_filter == true)      set mem_parts = "$mem_parts DART/filter_ic*[.0]$member"
         if (-e ICE/iceinput_${member}.tar) set mem_parts = "$mem_parts ICE/iceinput_${member}.tar "
         tar -r -f batch${batch}${comp} $mem_parts
         echo "tar -r -f batch${batch}${comp} $mem_parts" 

      endif
      @ n++
   end

#  Doing all at once, outside loop, would give more efficient MS access
#  but
#     msrcp $opts "$p1"  batch*${comp} ${ms_dir}
#     $opts and "p1" contents are being interpretted as files, 
#     maybe because of the multi-file copy batch*...
   msrcp $opts "$p1"  batch${batch}${comp} ${ms_dir}/batch${batch}${comp}
   if ($batch < $nbatch) rm batch${batch}${comp}

   @ batch++
end
#
# Correct batch back to value of last batch, for use below.
@ batch = $batch - 1

# ls -l batch*${comp}
# echo "msrcp $opts $p1  batch*${comp} ${ms_dir}" >> saved_restart
# msrcp $opts "$p1"  batch*${comp} ${ms_dir}

wait

# Check to see if it's okay to remove DART directory
if ($do_filter != true) then
   set list = `ls -l DART/filter_ic`
   set local_size = $list[5]
   set list = `msls -l ${ms_root}/DART/filter_ic`
   set ms_size = $list[5]
   echo "local_size ms_size = $local_size $ms_size" >> saved_restart
   if ($local_size != $ms_size) set ok_to_remove = false
endif

# Check to see if it's okay to remove CAM/CLM directories
set list = `ls -l batch${batch}${comp}`
set local_size = $list[5]
set list = `msls -l ${ms_root}/batch${batch}${comp}`
set ms_size = $list[5]
echo " CAM/CLM(/DART/ICE) local_size ms_size = $local_size $ms_size" >> saved_restart
if ($local_size != $ms_size) set ok_to_remove = false

if ($ok_to_remove == true) then
   echo "Archived files with write password $write_pass" >> saved_restart
   echo "msrcp of ${ms_root} succeeded; REMOVING $obs_seq DART,CAM,CLM,ICE" >> saved_restart
   rm -rf DART CAM CLM ICE
   # rm batch${batch}${comp}
   rm batch*${comp}
else
   echo "msrcp of ${ms_root} failed; NOT removing $obs_seq DART,CAM,CLM,ICE" >> saved_restart
endif

chmod 444 saved_restart
#================================================

echo "finished with auto_re2ms at " `date` 

exit
