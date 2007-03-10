#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source: /home/thoar/CVS.REPOS/DART/models/cam/shell_scripts/ms2restart.csh,v $
# $Name:  $

# script for copying 1 day/obs_seq of restart files from the mass store.
# CAM,CLM, and possibly filter_ic for each ensemble member are lumped together
# so that we can retrieve a subset of the ensemble members for a new experiment.
# Each batch has several ensemble member groups.

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/local/dcs/lib
echo $LD_LIBRARY_PATH

# set  echo verbose

if ($#argv < 3) then
   echo "usage; from case/experiment/obs_seq directory; "
   echo "       ms2restart MSpathname(no mss:) batch1 batchn uncompress(optnl)"
   echo "       untars ensemble members and puts them where they belong"
   echo "       Uncompresses caminput and clminput files, if told to."
   exit
endif 

set ms_root = $1
set ms_dir = mss:$ms_root

msls -l $ms_root
if ($status == 0) then
   echo "files will be read from $ms_root/batch#"
else
   echo "$ms_root does not exist.  Check name and try again"
   exit
endif

set uncomp = false
if ($#argv == 4) set uncomp = true

set batchn = $3
set batch = $2

while($batch <= $batchn)

   echo "msrcping ${ms_dir}/batch${batch}"
   if ($uncomp == true) then
      msrcp ${ms_dir}/batch${batch}.cmp batch$batch
   else
      msrcp ${ms_dir}/batch$batch .
   endif

#  This makes CAM CLM DART directories, as necessary, to place its files.
   echo "untarring"
   tar x -f batch$batch
   rm batch$batch
   @ batch++
end

# Get filter ics when all were written to one file.
# Also get assim_tools_ics, if it exists.
# DART may exist from untarring, above.
mkdir DART
if (! -e DART/filter_ic.0001) then
   echo "msrcping single filter_ic and/or inflation ic file(s)"
   msls ${ms_root}/DART/ic_files.tar
   if ($status == 0) then
      msrcp ${ms_dir}/DART/ic_files.tar DART
   endif
endif

if ($uncomp == true) then
   gunzip -r CAM
   gunzip -r CLM
endif 

echo "msrcp done; check local files "

exit
