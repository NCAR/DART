#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# script for copying 1 day/obs_seq of restart files from the mass store.
# CAM,CLM, and possibly filter_ic for each ensemble member are lumped together
# so that we can retrieve a subset of the ensemble members for a new experiment.
# Each batch has several ensemble member groups.

# set  echo verbose

if ($#argv < 3) then
   echo "usage; from case/experiment/obs_seq directory; "
   echo "       hpss2restart.csh HPSSpathname(no mss:) batch1 batchn uncompress(optnl)"
   echo "       untars ensemble members and puts them where they belong"
   echo "       Uncompresses caminput and clminput files, if told to."
   exit
endif 

set hpss_root = $1
set hpss_dir = $hpss_root

hsi -P ls -l $hpss_root
if ($status == 0) then
   echo "files will be read from $hpss_root/batch#"
else
   echo "$hpss_root does not exist.  Check name and try again"
   exit
endif

set uncomp = false
if ($#argv == 4) set uncomp = true

set batchn = $3
set batch = $2

while($batch <= $batchn)

   echo "getting ${hpss_dir}/batch${batch} to batch$batch"
   if ($uncomp == true) then
      hsi get batch$batch : ${hpss_dir}/batch${batch}.cmp
   else
      hsi get batch$batch : ${hpss_dir}/batch$batch
   endif

   # This makes CAM CLM DART directories, as necessary, to place its files.
   echo "untarring"
   tar -x -f batch$batch
   if ($status == 0) then
      rm batch$batch
   else
      echo "Could not untar batch$batch"
   endif
   @ batch++
end

# Get filter ics when all were written to one file.
# Also get assim_tools_ics, if it exists.
# DART may exist from untarring, above.
if (! -d DART) then
   mkdir DART
endif
if (! -f DART/filter_ic.0001) then
   echo DART has no filter_ic.0001
   ls DART
   hsi -P ls ${hpss_root}/DART/filter_ic 
   if ($status == 0) 
      echo "Getting single filter_ic"
      cd DART
      hsi get filter_ic : ${hpss_dir}/DART/filter_ic 
      cd ..
   else
      echo 'no filter_ic (nor filter_ic.0001) found'
   endif 

   hsi -P ls ${hpss_root}/DART/assim_tools_ics
   if ($status == 0) then
      cd DART
      hsi get assim_tools_ics : ${hpss_dir}/DART/assim_tools_ics  
   else
      echo 'no assim_tools_ics found'
   endif
endif

if ($uncomp == true) then
   gunzip -r CAM
   gunzip -r CLM
endif 

echo "Finished fetching files; check local files "

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


