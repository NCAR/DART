#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# script for copying 1 day/obs_seq of restart files to the mass store.
# CAM,CLM, and possibly filter_ic for each ensemble member are lumped together
# so that we can retrieve a subset of the ensemble members for a new experiment.
# Then it lumps together ensemble members into batches to reduce the number of files.

# updated for HPSS which replaced MSS

# set  echo verbose


if ($#argv < 2) then
   echo "usage; from case/experiment/obs_seq directory; "
   echo "       auto_ICs2hpss.csh num_ens_members num_per_batch compress(optnl)"
   echo "       Compresses caminput and clminput files."
   echo "       tars ensemble members together into batches"
   exit
else if ($#argv == 2) then
   set comp = ' '
else if ($#argv == 3) then
   set comp = .cmp
endif 

set num_ens = $1
set num_per_batch = $2

if ($num_per_batch > $num_ens) then
  echo "auto_ICs2hpss; num_per_batch > num_ens"
  exit -1
endif

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

# New section for archiving initial ensembles from a spin-up
set months = (Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec)
echo $obs_seq >! obs_seq_name

# For ICs every 5 days
set STRING = "1,$ s#-# #g"
set parts = `sed -e "$STRING" obs_seq_name`
\rm obs_seq_name 
set month = $parts[1]
set day   = $parts[2]
set date = $months[$month]_$day
# END of spin-up section

set hpss_root = /RAEDER/DAI/CAM_init/$case/${date}/${1}x${2}
set hpss_dir = $hpss_root
echo "files will be written to ${hpss_root}/batch#" >! saved_restart
echo "files will be written to ${hpss_root}/batch#" 


# Figure out how many files to divide ensemble members among
@ nbatch = $num_ens / $num_per_batch
if ($num_ens % $num_per_batch != 0 ) @ nbatch++

set ok_to_remove = true

# Tally up and list the extra (non-filter) ic files
# Add more if blocks as needed.
set DART_list = ' '
set num_files = 0
set do_filter = false
if (-d DART) then
   if (-e DART/assim_tools_ics) then
      set DART_list = ($DART_list assim_tools_ics)
      @ num_files++
   endif
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
endif

echo "DART_list = ($DART_list) " >> saved_restart
echo "DART_list = ($DART_list) "

# make sure directory exists first
hsi mkdir ${hpss_dir}

# convert wordlist into a single character string (insert ,s)
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

if ($DART_files != 'none' && do_filter == false) then
   tar -c -f ic_files.tar DART/{$DART_files}
   hsi mkdir ${hpss_dir}/DART
   hsi put ic_files.tar : ${hpss_dir}/DART/ic_files.tar 
endif

# Do this here, manually, or within tar; -z not available on tempest
# also, not effective to compress filter_ic files (big and dense)
set fil = CAM/caminput_1.*
if ($comp != ' ' &&  $fil:e != 'gz') then
   gzip -r CAM
   gzip -r CLM
   echo "compressed CAM and CLM" >> saved_restart
endif

set batch = 1
while($batch <= $nbatch)

   @ base = (($batch - 1) * $num_per_batch) + 1
   @ member = $base
   if ($batch == 1 && (-e CAM/caminput_0.nc || -e CAM/caminput_0.nc.gz)) \
      set member = 0
# BUT will there be a filter_ic.0000 ?
   echo base and member $base $member >> saved_restart
   echo base and member $base $member 

   # create the tar file using the first ensemble member of this batch
   if ($do_filter == true) then
      # if (-e DART/assim_tools_ics && $member <= $base) then
      if ($DART_files != 'none' && $member <= $base) then
# won't work for $member < $base?            
        tar -c -f batch${batch}${comp} CAM/caminput_$member.nc* \
                 CLM/clminput_$member.nc* DART/filter_ic*[.0]$member DART/{$DART_files}
        echo "tar -c -f batch${batch}${comp} CAM/caminput_$member.nc* \
                 CLM/clminput_$member.nc* DART/filter_ic*[.0]$member DART/{$DART_files} "
      else
        tar -c -f batch${batch}${comp} CAM/caminput_$member.nc* \
                 CLM/clminput_$member.nc* DART/filter_ic*[.0]$member
        echo "tar -c -f batch${batch}${comp} CAM/caminput_$member.nc* \
                 CLM/clminput_$member.nc* DART/filter_ic*[.0]$member "
      endif
   else
        tar -c -f batch${batch}${comp} CAM/caminput_${member}.nc* \
                  CLM/clminput_${member}.nc* ICE/iceinput_${member}.tar
   endif
   # tack on additional ens members until this batch is complete
   set n = 2
   while ($n <= $num_per_batch)
     @ member++ 
# This looks wrong; $n will never be big enough, unless something crazy happens.
#     if ($n > $num_ens) then
     if ($member > $num_ens) then
        @ batch = $nbatch + 1
        @ n = $num_per_batch + 1
     else
        if ($do_filter == true) then
           tar -r -f batch${batch}${comp} CAM/caminput_$member.nc* \
               CLM/clminput_$member.nc* DART/filter_ic*[.0]$member
           echo "tar -r -f batch${batch}${comp} CAM/caminput_$member.nc* \
               CLM/clminput_$member.nc* DART/filter_ic*[.0]$member "
        else
           tar -r -f batch${batch}${comp} CAM/caminput_${member}.nc* \
                     CLM/clminput_${member}.nc* ICE/iceinput_${member}.tar
        endif
     endif
     @ n++
   end

   hsi put batch${batch}${comp} : ${hpss_dir}/batch${batch}${comp}
   if ($batch < $nbatch) rm batch${batch}${comp}

   @ batch++
end
# Correct batch back to value of last batch, for use below.
@ batch = $batch - 1

# put them one at a time?
foreach file ( batch*${comp} )
  hsi put ${file} : ${hpss_dir}/${file}
end


wait

# Check to see if it's okay to remove DART directory
if ($do_filter != true && $num_files > 0) then
   set list = `ls -l DART/filter_ic`
   set local_size = $list[5]
   set list = `hsi -P ls -l ${hpss_root}/DART/filter_ic`
   set hpss_size = $list[5]
   echo "local_size hpss_size = $local_size $hpss_size" >> saved_restart
   if ($local_size != $hpss_size) set ok_to_remove = false
endif

# Check to see if it's okay to remove CAM/CLM directories
set list = `ls -l batch${batch}${comp}`
set local_size = $list[5]
set list = `hsi -P ls -l ${hpss_root}/batch${batch}${comp}`
set hpss_size = $list[5]
echo " CAM/CLM/ICE(/DART) local_size hpss_size = $local_size $hpss_size" >> saved_restart
if ($local_size != $hpss_size) set ok_to_remove = false


if ($ok_to_remove == true) then
   echo "Archived files " >> saved_restart
   echo "archive of ${hpss_root} succeeded; REMOVING $obs_seq DART,CAM,CLM,ICE" >> saved_restart
   rm -rf DART CAM CLM ICE
   rm batch${batch}${comp}
else
   echo "archive of ${hpss_root} failed; NOT removing $obs_seq DART,CAM,CLM,ICE" >> saved_restart
endif

chmod 444 saved_restart

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

