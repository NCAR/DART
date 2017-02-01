#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# KDR This is a copy of /glade/p/image/Observations/ACARS/progs/make6h.sh, 
#     adapted for SABER.  This script needs 
#        program 'advance_time' (e.g. from $dart_cam/work)
#        program 'obs_sequence_tool' (e.g. from $dart/observations/SABER/work)
#        input file input.nml.template.h (e.g. from $obs/ACARS/progs) 
#        obs_seq files (daily) in a neighbor directory named with format ../YYYYMM.
#     Submit/run this script in a directory which is a neighbor 
#     of the obs_seq files directory, generally a 'progs' directory.
#     On yellowstone
#     > bsub < make6h.sh
#     The 6 hourly files will be put in a new directory with the same name
#     as where the daily files are, but with _6H appended; e.g. ../YYYYMM_6H.
#     Make sure the the daily files are in place, and uncompressed.

# split a set of daily files into 6 hr files which start at:
# 03:00:01Z, 09:00:01Z, 15:00:01Z, and 21:00:01Z each day.
# the daily input files run from 03:00:01Z to 03:00:00 the next day.
# for those files which cross 0Z, the date in the filename is
# the current day at the start of the file.

#BSUB -J splitobs
#BSUB -W 0:20
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -q caldera
#BSUB -P P86850054
#BSUB -n 1


# set the first and last days to be split.  can roll over
# month and year boundaries now!   note that for the last day
# the data from the following day must be available.

let start_year=2008
let start_month=8
let start_day=2

let end_year=2008
let end_month=8
let end_day=3

# end of things you should have to set in this script

# convert the start and stop times to gregorian days, so we can
# compute total number of days including rolling over month and
# year boundaries.  make sure all values have leading 0s if they
# are < 10.  do the end time first so we can use the same values
# to set the initial day while we are doing the total day calc.

# make sure there is an initial input.nml for advance_time
cp -f input.nml.template.h input.nml

# these outputs from advance time (with the -g flag) are
# 2 integers: gregorian_day_number seconds
# for now this script is doing whole days at a time so for
# the start/stop limits we are ignoring the hours here.
# we are setting 3Z to 3Z on the limits.
mon2=`printf %02d $end_month`
day2=`printf %02d $end_day`
end_d=(`echo ${end_year}${mon2}${day2}03 0 -g | ./advance_time`)
#echo $end_d

mon2=`printf %02d $start_month`
day2=`printf %02d $start_day`
start_d=(`echo ${start_year}${mon2}${day2}03 0 -g | ./advance_time`)
#echo $start_d

# these are a string in the format YYYYMMDDHH
# do them here to prime the loop below which first takes them apart.
currtime=(`echo ${start_year}${mon2}${day2}03  0  | ./advance_time`)
nexttime=(`echo ${start_year}${mon2}${day2}03 +6h | ./advance_time`)
middtime=(`echo ${start_year}${mon2}${day2}03 +3h | ./advance_time`)

# how many total days are going to be split (for the loop counter)
# and we are going to make 4 files per day, so times 4.  if we change
# the increment (6h currently) the multiple will have to change as well.
let totaldays=${end_d}-${start_d}+1
let totaltimes=${totaldays}*4
#echo totaltimes = $totaltimes

# loop over each new output file
let t=1
while (( t <= totaltimes)) ; do
  echo t=$t
#echo top of loop
#echo $currtime
  # parse out the parts from a string which is YYYYMMDDHH
  # both for the current time, the end time, and the middle of
  # the interval (the middle is for the filename; the others
  # for the obs time limits.)
  cyear=${currtime:0:4}
  cmonth=${currtime:4:2}
  cday=${currtime:6:2}
  chour=${currtime:8:2}
  nyear=${nexttime:0:4}
  nmonth=${nexttime:4:2}
  nday=${nexttime:6:2}
  nhour=${nexttime:8:2}
  myear=${middtime:0:4}
  mmonth=${middtime:4:2}
  mday=${middtime:6:2}
  mhour=${middtime:8:2}
#echo curr $cyear $cmonth $cday $chour
#echo next $nyear $nmonth $nday $nhour
#echo midd $myear $mmonth $mday $mhour

  # compute the equivalent gregorian days here for the current time 
  # plus 1 second (the interval start) and the end time.
  g=(`echo ${cyear}${cmonth}${cday}${chour} +1s -g | ./advance_time`)
  greg0=${g[0]}
  let secs0=${g[1]}
  g=(`echo ${cyear}${cmonth}${cday}${chour} +6h -g | ./advance_time`)
  greg1=${g[0]}
  let secs1=${g[1]}
#echo $greg0 $secs0
#echo $greg1 $secs1

  echo next output file ${myear}${mmonth}${mday}${mhour} 
  echo "starting time" ${cyear}${cmonth}${cday}${chour} plus 1 sec
  echo "ending   time" ${nyear}${nmonth}${nday}${nhour} 
  #echo gregorian range = $greg0 $secs0 thru $greg1 $secs1

  # make sure output dir exists
  OUTDIR=../${myear}${mmonth}_6H
  if [[ ! -d ${OUTDIR} ]] ; then
       mkdir ${OUTDIR}
  fi

  # BUFR:
  # in that case we've got input files that wrap around to 3Z the 
  # next day so there's only 1 input file ever.
  # echo "../${cyear}${cmonth}/obs_seq${cyear}${cmonth}${cday}" > olist
  # SABER: 
  # Input files range from 0Z to 0Z; need 2 input files for 0Z windows.
  # obs_seq.saber_v2.0_20080803
  if [[ -f olist ]] ; then
     rm olist;
  fi
  touch olist;

  if [[ -f ../${currtime:0:6}/obs_seq.saber_v2.0_${currtime:0:8} ]] ; then
     echo "../${currtime:0:6}/obs_seq.saber_v2.0_${currtime:0:8}" >> olist;
  else
     echo "WARNING; ${currtime:0:6}/obs_seq.saber_v2.0_${currtime:0:8}";
     echo "         is not available.  Final obs_seq file will not have hours 0-3";
  fi
  if [[ ${mhour} -eq 0 ]] ; then
     # nday2=`printf %02d ${nday}` ;
     # echo "../${cyear}${cmonth}/obs_seq.saber_v2.0_${cyear}${cmonth}${nday2}" >> olist;
     echo "../${nexttime:0:6}/obs_seq.saber_v2.0_${nexttime:0:8}" >> olist;
  fi

  echo "olist = "
  cat olist

  sed -e "s#YYYY#${myear}#g"    \
      -e "s#MM#${mmonth}#g"     \
      -e "s#DD#${mday}#g"       \
      -e "s#HH#${mhour}#g"      \
      -e "s#OUTDIR#${OUTDIR}#g" \
      -e "s#GREG0#${greg0}#g"   \
      -e "s#GREG1#${greg1}#g"   \
      -e "s#SECS0#${secs0}#g"   \
      -e "s#SECS1#${secs1}#g"   < ./input.nml.template.h > input.nml

  # cat input.nml

  # do the extract here
  #echo calling ./obs_sequence_tool here
  ./obs_sequence_tool
  

  # advance 6 hours; the output is YYYYMMDDHH
  currtime=$nexttime
  nexttime=(`echo ${nyear}${nmonth}${nday}${nhour} +6h | ./advance_time`)
  middtime=(`echo ${nyear}${nmonth}${nday}${nhour} +3h | ./advance_time`)
#echo $currtime $nexttime $middtime

  # advance the loop counter
  let t=t+1
 
done

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

