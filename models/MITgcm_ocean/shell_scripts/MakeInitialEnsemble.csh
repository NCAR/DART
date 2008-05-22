#!/bin/tcsh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2008, Data Assimilation Research Section, 
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# $Id: runmodel_1x 2799 2007-04-04 23:17:51Z thoar $
#
#=============================================================================
# So the background is that we ran the model for 14 days with a 900 second
# timestep - starting 1996 01 01 ??Z 
#
# I wanted 20 ensemble members.
#
# I output snapshot files every 12 hours to generate 29 sets of snapshots.
# We're just going to make the first snapshot into ensemble member 1. 
# Second snapshot will be ensemble member 2 ... I'm just not going to use
# the last snapshots. This may not be the greatest idea - I don't know
# how much spread we will have to start. There must be SOME ...
# depends on the nonlinearity of the model.
#
# I modified trans_pv_sv.f90 to ignore the timestamps 
#
# Repeatedly call trans_pv_sv to turn each snapshot into a 
# DART initial condition. The first (line - for ascii files) or 
# (4+4+4+4 = 16 bytes - for binary files) contain the timestamp.
# We need the timestamp for the first one to be saved and used for all
# subsequent 'members'. Simply 'cat' them together at the end.
# There is no problem with having more ensemble members in a DART
# initial condition file ... but there is a big problem if you don't have
# enough ...
#
# TJH - 21 May 2008

# trans_pv_sv needs three namelist files:  data, data.cal and input.nml

set DARTROOT = /fs/image/home/thoar/SVN/DART/models/MITgcm_ocean
cp -p ${DARTROOT}/inputs/data      .
cp -p ${DARTROOT}/inputs/data.cal  .
cp -p ${DARTROOT}/work/input.nml   .

@ memcount = 0

foreach FILE ( T.*.data )
   @ memcount = $memcount + 1

   set FILEBASE = $FILE:r
   set TIMESTEP = $FILEBASE:e
   echo "Converting snapshot timestep $TIMESTEP ..."

   echo $TIMESTEP | ${DARTROOT}/work/trans_pv_sv

   if ( $memcount == 1) then
   #  head -1 assim_model_state_ud >! timestamp
   else
   #  set nlines = `wc -l assim_model_state_ud`
   endif

   if ( $memcount < 10 ) then
      set OFNAME = ens_mem_00$memcount
   else if ( $memcount < 100 ) then
      set OFNAME = ens_mem_0$memcount
   else
      set OFNAME = ens_mem_$memcount
   endif

   mv assim_model_state_ud $OFNAME
end
