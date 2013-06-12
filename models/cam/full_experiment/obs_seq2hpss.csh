#!/usr/local/bin/tcsh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Script for archiving months of obs_seq.final files 
# (for easier access than diagnostics.tar.gz)
#
#BSUB -J obs_seq2hpss
#BSUB -o obs_seq2hpss.%J.log
#BSUB -P xxxxxxxx
#BSUB -q share
#BSUB -W 3:00
#BSUB -n 1
#BSUB -R "span[ptile=1]"

if ($?LS_SUBCWD) then
   cd $LS_SUBCWD
endif

set obs_seq_first = 1
set mo_first = 8
set mo_last  = 8
set input_root = obs_
set digits = 4
set obs_seq_freq = 1

set saved = saved_obs_seq_${mo_first}-${mo_last}

set remove_finals = true

set year = 2006

set months =  (Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec)
set mo_days = (31  28  31  30  31  30  31  31  30  31  30  31)

#-----------------------------------------------------------------------------
# Parse parts of the path name for construction of HPSS name.
set direct = `pwd`
set exp_dir = $direct:t

cd ..
set direct = `pwd`
set case = $direct:t

cd ${exp_dir}

set hpss_root = /RAEDER/DAI/$case/${exp_dir}
set hpss_dir = $hpss_root
echo "files will be written to ${hpss_root}/MM_obs_seq.tar" >> $saved

#
#-----------------------------------------
# Figure out how many months to archive
if ($mo_last >= $mo_first) then
   @ mo_num = ($mo_last - $mo_first) + 1
else
   @ mo_num = $mo_last + (12 - $mo_first) + 1
endif

#-----------------------------------------

set ok_to_remove = true

set mo = 1
set mo_cal = $mo_first
@ obs_seq_n = $obs_seq_first - 1
echo "starting mo loop with mo mo_cal mo_num obs_seq_n = $mo $mo_cal $mo_num $obs_seq_n" >> $saved
# set echo verbose
while($mo <= $mo_num)
   if ($mo_cal == 13) then
      @ year++
      set mo_cal = 1
   endif
      
   @ leap = $year % 4
   if ($leap == 0) set mo_days[2] = 29

   @ obs_seq_1 = $obs_seq_n + 1
   @ obs_seq_n = ($obs_seq_1 + ($mo_days[${mo_cal}] * $obs_seq_freq)) - 1

   # create the tar file using the first obs_seq/day of this month
   set obs_seq = $obs_seq_1
   set input_dir = `printf "%s%0${digits}d" ${input_root} $obs_seq`

   set out_file = $months[${mo_cal}]_obs_seq.tar 
   echo "out_file and input_dir = "$out_file $input_dir >> $saved

   if (-e ${input_dir}/obs_seq.final) then
      tar -c -f $out_file ${input_dir}/{obs_seq.final,input.nml}
   else
      echo "${input_dir}/obs_seq.final does not exist; exiting" >> $saved
      exit
   endif

   # tack on additional ens members until this month is complete
   @ obs_seq++
   echo "starting obs_seq loop " >> $saved

   while ($obs_seq <= $obs_seq_n)
     set input_dir = `printf "%s%0${digits}d" ${input_root} $obs_seq`
     if ( -e ${input_dir}/obs_seq.final) then
        tar -r -f $out_file ${input_dir}/{obs_seq.final,input.nml}
     else
        # finish up
        echo "obs_seq $obs_seq does not exist; finishing" >> $saved
        @ mo = $mo_num + 1
        @ obs_seq = $obs_seq_n + 1
        if ($mo_cal != $mo_last) then
           echo "DID NOT FINISH ARCHIVE; $mo_cal is last month done" >> $saved
        endif
     endif
     @ obs_seq++
   end

   hsi mkdir ${hpss_dir}
   hsi put $out_file : ${hpss_dir}/$out_file

# Check to see if it's okay to remove obs_seq.final, etc
   set list = `ls -l $out_file`
   set local_size = $list[5]
   set list = `hsi -P ls -l ${hpss_root}/$out_file`
   set hpss_size = $list[5]
   echo "for $months[${mo_cal}] local_size, hpss_size = $local_size $hpss_size" >> $saved
   if ($local_size == $hpss_size) then
      rm $out_file &
      echo "hsi put of ${hpss_root}/$out_file SUCCEEDED " >> $saved
      echo "Archived obs_seq.final files                " >> $saved
      if ($remove_finals == true) then
         set obs_seq = $obs_seq_1
         while ($obs_seq <= $obs_seq_n)
            set input_dir = `printf "%s%0${digits}d" ${input_root} $obs_seq`
            rm ${input_dir}/[^is]* &
            @ obs_seq++
         end
         echo "REMOVED $out_file, ./obs_seq.final " >> $saved
      endif
   else
      echo "hsi put failed; NOT removing $months[${mo_cal}]: $obs_seq_1 - $obs_seq_n " >> $saved
      set mo = $mo_num + 1
   endif

   @ mo++
   @ mo_cal++
end

wait

chmod 444 $saved

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

