#!/bin/csh

# script for copying 1 day/obs_seq of restart files to the mass store.
# CAM,CLM, and possibly filter_ic for each ensemble member are lumped together
# so that we can retrieve a subset of the ensemble members for a new experiment.
# Then it lumps together ensemble members into batches to reduce the number of files.

# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$

# change for lightning?
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/local/dcs/lib
echo $LD_LIBRARY_PATH

# set  echo verbose


if ($#argv < 2) then
   echo "usage; from case/experiment/obs_seq directory; "
   echo "       restart2ms num_ens_members num_per_batch compress(optnl)"
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

# set proj_num = 86850054
set proj_num = 39510050
set ret_period = 1000
set write_pass = da$$
echo "with write password $write_pass" > saved_restart

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
set ms_root = /RAEDER/DAI/$case/${exp_dir}/${obs_seq}/${1}x${2}
set ms_dir = mss:$ms_root
echo "files will be written to ${ms_root}/batch#" >> saved_restart

# Figure out how many files to divide ensemble members among
@ nbatch = $num_ens / $num_per_batch
if ($num_ens % $num_per_batch != 0 ) @ nbatch++

set ok_to_remove = true

if (-e DART/filter_ic) then
   set do_filter = false
   if (-e DART/assim_tools_ics) then
      msrcp -pe $ret_period -pr $proj_num -wpwd $write_pass \
            DART/filter_ic DART/assim_tools_ics ${ms_dir}/DART &
   else
      msrcp -pe $ret_period -pr $proj_num -wpwd $write_pass DART/filter_ic ${ms_dir}/DART/filter_ic &
   endif
else if (-e DART/filter_ic*0$num_ens || -e DART/filter_ic*.$num_ens) then
   set do_filter = true
else
   echo 'NOT ENOUGH ENSEMBLE MEMBERS IN .../DART' >> saved_restart
   exit
endif

# Do this here, manually, or within tar; -z not available on tempest
# also, not effective to compress filter_ic files (big and dense)
if ($comp != ' ') then
   gzip -r CAM
   gzip -r CLM
endif

set batch = 1
while($batch <= $nbatch)

   @ base = (($batch - 1) * $num_per_batch) + 1
   @ member = $base
   if ($batch == 1 && (-e CAM/caminput_0.nc || -e CAM/caminput_0.nc.gz)) \
      set member = 0
# BUT will there be a filter_ic.0000 ?
   echo base and member $base $member >> saved_restart

   # create the tar file using the first ensemble member of this batch
   if ($do_filter == true) then
      if (-e DART/assim_tools_ics && $member <= $base) then
# won't work for $member < $base?            
        tar c -f batch${batch}${comp} CAM/caminput_$member.nc* \
                 CLM/clminput_$member.nc* DART/filter_ic*[.0]$member DART/assim_tools_ics 
      else
        tar c -f batch${batch}${comp} CAM/caminput_$member.nc* \
               CLM/clminput_$member.nc* DART/filter_ic*[.0]$member
      endif
   else
        tar c -f batch${batch}${comp} CAM/caminput_$member.nc* CLM/clminput_$member.nc*
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
           tar r -f batch${batch}${comp} CAM/caminput_$member.nc* \
               CLM/clminput_$member.nc* DART/filter_ic*[.0]$member
        else
           tar r -f batch${batch}${comp} CAM/caminput_$member.nc* CLM/clminput_$member.nc*
        endif
     endif
     @ n++
   end

   msrcp -pe $ret_period -pr $proj_num -wpwd $write_pass batch${batch}${comp} ${ms_dir}/batch${batch}${comp}
   if ($batch < $nbatch) rm batch${batch}${comp}

   @ batch++
end
# correct batch back to value of last batch, for use below.
@ batch = $batch - 1

# msrcp -pe $ret_period -pr $proj_num -wpwd $write_pass batch*${comp} ${ms_dir}
mscomment -R -wpwd $write_pass -c "write password $write_pass" ${ms_root}

wait

# check to see if it's okay to remove DART directory
if ($do_filter != true) then
   set list = `ls -l DART/filter_ic`
   set local_size = $list[5]
   set list = `msls -l ${ms_root}/DART/filter_ic`
   set ms_size = $list[5]
   echo "local_size ms_size = $local_size $ms_size" >> saved_restart
   if ($local_size != $ms_size) set ok_to_remove = false
endif

# check to see if it's okay to remove CAM/CLM directories
set list = `ls -l batch${batch}${comp}`
set local_size = $list[5]
set list = `msls -l ${ms_root}/batch${batch}${comp}`
set ms_size = $list[5]
echo " CAM/CLM(/DART) local_size ms_size = $local_size $ms_size" >> saved_restart
if ($local_size != $ms_size) set ok_to_remove = false


if ($ok_to_remove == true) then
   echo "Archived files with write password $write_pass" >> saved_restart
   echo "msrcp of ${ms_root} succeeded; REMOVING $obs_seq DART,CAM,CLM" >> saved_restart
   rm -rf DART CAM CLM
   rm batch${batch}${comp}
else
   echo "msrcp of ${ms_root} failed; NOT removing $obs_seq DART,CAM,CLM" >> saved_restart
endif

chmod 444 saved_restart

exit
