#!/bin/ksh

# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
#--------------------------------------------------------------
# DESCRIPTION:
# Merge NCEP BUFR obs (inc ACARS), GPS RO obs, and AIRS T/Q obs.
#
# This script merges a set of input obs files into 3-hour files.  
# Each file is named by the center time of the data and uses
# the CESM time convention for the name format.  For example, 
# a file with a name ending in YYYY-MM-DD-21600 will contain data
# starting at 04:30:01Z and ending at 07:30:00Z.  Note that a file
# with a name ending in YYYY-MM-DD-00000 contains data from
# 22:30:01Z the previous day to 01:30:00Z of the day YYYY-MM-DD.
#
# this version assumes the NCEP+ACARS input files were created
# as 6H files with filenames ending in 06, 12, 18, 00.
# it assumes daily GPS obs, and daily AIRS obs.
#
# be *very* careful to check that you have the right files for the
# first and last day -- many datasets that are distributed as collections 
# of weekly or monthly files, and you cannot generate obs for midnight
# on days where you don't have the previous day's data.
#
# this script runs one conversion at a time.  it has not been
# updated to run multiple merge jobs at once.
#

#--------------------------------------------------------------
# SLRUM directives             sbatch script.csh
#
# sinfo     information about the whole slurm system
# squeue    information about running jobs
# sbatch    submit a job
# scancel   kill a job
#
#SBATCH --ignore-pbs
#SBATCH --job-name=mergeobs
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH -A P86850054
#SBATCH -p dav
#SBATCH -C casper
#SBATCH -e mergeobs.%j.err
#SBATCH -o mergeobs.%j.out
#
#--------------------------------------------------------------
# PBS directives                qsub script.csh
#
# qstat    information on running jobs
# qsub     submit a job
# qdel     kill a job
# qpeek    see output from a running job
#
#PBS -N mergeobs
#PBS -l walltime=06:00:00
#PBS -q share
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -A P86850054
#
#--------------------------------------------------------------
# LSF directives                bsub < script.csh
#
# bstat    information on running jobs
# bsub     submit a job
# bdel     kill a job
# bpeek    see output from a running job
#
#BSUB -J mergeobs
#BSUB -o mergeobs.%J.log
#BSUB -q share
#BSUB -n 1
#BSUB -W 0:10:00
#BSUB -P P86850054
#
#--------------------------------------------------------------

# set the first and last days/hours.  can roll over month and year boundaries.
# if start_hour is 12Z and step is 3H, then the first file will be centered 
# at 12Z on this day, with data starting at 10.5Z and ending at 13.5Z.
#
# set these WITHOUT leading 0s.  they are numeric values, not strings.
#
# if you set the start hour to 0, you must have data files for
# the day immediately preceeding this date.

# on cheyenne, this is taking about 2 hours for a month.

let start_year=2011
let start_month=9
let start_day=1 
let start_hour=0


let end_year=2011
let end_month=9
let end_day=3
let end_hour=21

# string to control the number of hours in each file.
# set the forward and backwards half steps to be consistent.
# (the +1s will be done when setting the namelist seconds)
step='+3h'
halfstep='+1h30m'
halfback='-1h30m'

# end of things you should have to set in this script


# before running advance_time there must be an input.nml
# with at least a utilities namelist.  it will be overwritten
# in the first real loop by the input.nml.template, so make
# any permanent changes the template and not to input.nml.

rm -f input.nml
echo '&utilities_nml' > input.nml
echo '/' >> input.nml

# convert the start and stop times to gregorian days, so we can
# compute total number of days including rolling over month and
# year boundaries.  make sure all values have leading 0s if they
# are < 10.  do the end time first so we can use the same values
# to set the initial day while we are doing the total day calc.

# the output of advance time with the -g input is:
#   gregorian_day_number  seconds
# use ${var[0]} to return the day number, ${var[1]} for the seconds

mon2c=`printf %02d $end_month`
day2c=`printf %02d $end_day`
hor2c=`printf %02d $end_hour`
endstring=${end_year}${mon2c}${day2c}${hor2c}00
set -A end_t `echo $endstring 0 -g | ./advance_time`

mon2c=`printf %02d $start_month`
day2c=`printf %02d $start_day`
hor2c=`printf %02d $start_hour`
startstring=${start_year}${mon2c}${day2c}${hor2c}00
set -A start_t `echo $startstring 0 -g | ./advance_time`

# count up how many steps we are going to be executing

curtime=$startstring
let totalsteps=0

while true ; do

  # compute the current gregorian day here.  
  #  this_t[0] is day, this_t[1] is seconds
  set -A this_t `echo $curtime 0 -g | ./advance_time`

  # if the current day is beyond the end day, we are done.  break out of loop.
  # if the current day is equal to end day, check the seconds to see if we quit.
  # otherwise, do the loop again.
  if [[ ${this_t[0]} -gt ${end_t[0]} ]] ; then break ; fi
  if [[ ${this_t[0]} -eq ${end_t[0]}  &&  ${this_t[1]} > ${end_t[1]} ]] ; then break ; fi

  # advance to the next time and increment the step count
  curtime=`echo $curtime $step | ./advance_time`
  let totalsteps=$totalsteps+1

done

echo This script is going to create $totalsteps DART obs sequence files.
echo Starting at $start_year $start_month $start_day $start_hour 
echo Ending at $end_year $end_month $end_day $end_hour
echo ''

# the output of this call is a string YYYYMMDDHH
# see below for help in how to easily parse this up into words
  curstep=`echo $startstring         0 | ./advance_time`
 nextstep=`echo $startstring     $step | ./advance_time`
stepstart=`echo $startstring $halfback | ./advance_time`
  stepend=`echo $startstring $halfstep | ./advance_time`


# ok, let's actually do something.  up to now it has
# been computing fiddly time bits.

let s=1
while (( s <= totalsteps )) ; do

  # parse out the parts from a string which is YYYYMMDDHHMM
  # use cut with the byte option to pull out columns 1-4, 5-6, 7-8, 9-10, and 11-12
  # c = current middle-of-step, n = next middle, ss = step start, se = step end

   cyear=`echo $curstep | cut -b1-4`
  cmonth=`echo $curstep | cut -b5-6`
    cday=`echo $curstep | cut -b7-8`
   chour=`echo $curstep | cut -b9-10`
    cmin=`echo $curstep | cut -b11-12`

   nyear=`echo $nextstep | cut -b1-4`
  nmonth=`echo $nextstep | cut -b5-6`
    nday=`echo $nextstep | cut -b7-8`
   nhour=`echo $nextstep | cut -b9-10`
    nmin=`echo $nextstep | cut -b11-12`

  ssyear=`echo $stepstart | cut -b1-4`
 ssmonth=`echo $stepstart | cut -b5-6`
   ssday=`echo $stepstart | cut -b7-8`
  sshour=`echo $stepstart | cut -b9-10`
   ssmin=`echo $stepstart | cut -b11-12`

  seyear=`echo $stepend | cut -b1-4`
 semonth=`echo $stepend | cut -b5-6`
   seday=`echo $stepend | cut -b7-8`
  sehour=`echo $stepend | cut -b9-10`
   semin=`echo $stepend | cut -b11-12`

  # compute the equivalent gregorian days/secs here.
  set -A g `echo ${cyear}${cmonth}${cday}${chour}${cmin} 0 -g | ./advance_time`
  cgregday=${g[0]}; cgregsec=${g[1]}

  # special for dart: step start needs to be +1 second 
  set -A g `echo ${ssyear}${ssmonth}${ssday}${sshour}${ssmin} +1s -g | ./advance_time`
  ssgregday=${g[0]}; ssgregsec=${g[1]}

  set -A g `echo ${seyear}${semonth}${seday}${sehour}${semin} 0 -g | ./advance_time`
  segregday=${g[0]}; segregsec=${g[1]}

  # compute the CESM-style time string for the output filename
  cesmtime=`echo ${cyear}${cmonth}${cday}${chour} 0 -c | ./advance_time`

  # the next day, needed for observations centered at 21Z
  nextday=`echo ${cyear}${cmonth}${cday}00 +1d | ./advance_time | cut -b1-10`

  # compute the DOYs, and make sure they're 3 digits long
  set -A jan1 `echo ${cyear}010100 0 -g | ./advance_time`
  let cdy=${cgregday}-${jan1[0]}+1
  cdoy=`printf %03d $cdy`

  set -A jan1 `echo ${ssyear}010100 0 -g | ./advance_time`
  let ssdy=${ssgregday}-${jan1[0]}+1
  ssdoy=`printf %03d $ssdy`

  set -A jan1 `echo ${seyear}010100 0 -g | ./advance_time`
  let sedy=${segregday}-${jan1[0]}+1
  sedoy=`printf %03d $sedy`
  
  # status/debug - comment in or out as desired.
  echo "starting processing for "  ${cyear}  ${cmonth}  ${cday}  ${chour}  ${cmin}, gregorian day/sec:  $cgregday  $cgregsec #, DOY $cdoy 
  echo "          step start is " ${ssyear} ${ssmonth} ${ssday} ${sshour} ${ssmin}, gregorian day/sec: $ssgregday $ssgregsec #, DOY $ssdoy 
  echo "            step end is " ${seyear} ${semonth} ${seday} ${sehour} ${semin}, gregorian day/sec: $segregday $segregsec #, DOY $sedoy 
  echo ''


  # construct often-used strings below
     ym=${cyear}${cmonth}   ;   ymd=${ym}${cday}      ;   ymdh=${ymd}${chour}
   ssym=${ssyear}${ssmonth} ; ssymd=${ssym}${ssday}   ; ssymdh=${ssymd}${sshour}
   seym=${seyear}${semonth} ; seymd=${seym}${seday}   ; seymdh=${seymd}${sehour}

#echo
#echo     ym=$ym ; echo   ymd=${ymd}   ; echo   ymdh=${ymdh}
#echo   ssym=$ym ; echo ssymd=${ssymd} ; echo ssymdh=${ssymdh}
#echo   seym=$ym ; echo seymd=${seymd} ; echo seymdh=${seymdh}
#echo

  # ----------------------------------------
  # start of this merge-specific code

  # add whatever input observation files you need to the 'olist' file.
  # that sets the input for the merge.  do not list the same file twice or
  # those obs will be silently duplicated.  todo: could put in a 'sort|uniq'
  # to catch that.

  rm -f olist

  # the base NCEP + ACARS - 3Z to 3Z ; use the mid time to
  # get the right filename.  these are 6h files so we will
  # need to generate filenames to match the right intervals.
  case ${chour} in

   03 )
     echo /glade/p/cisl/dares/Observations/NCEP+ACARS/${ym}_6H/obs_seq${ymd}00 >> olist
     echo /glade/p/cisl/dares/Observations/NCEP+ACARS/${ym}_6H/obs_seq${ymd}06 >> olist ;;

   09 )
     echo /glade/p/cisl/dares/Observations/NCEP+ACARS/${ym}_6H/obs_seq${ymd}06 >> olist
     echo /glade/p/cisl/dares/Observations/NCEP+ACARS/${ym}_6H/obs_seq${ymd}12 >> olist ;;

   15 )
     echo /glade/p/cisl/dares/Observations/NCEP+ACARS/${ym}_6H/obs_seq${ymd}12 >> olist
     echo /glade/p/cisl/dares/Observations/NCEP+ACARS/${ym}_6H/obs_seq${ymd}18 >> olist ;;

   21 )
     echo /glade/p/cisl/dares/Observations/NCEP+ACARS/${ym}_6H/obs_seq${ymd}18   >> olist
     echo /glade/p/cisl/dares/Observations/NCEP+ACARS/${ym}_6H/obs_seq${nextday} >> olist ;;

    * )
     echo /glade/p/cisl/dares/Observations/NCEP+ACARS/${ym}_6H/obs_seq${ymdh} >> olist ;;

  esac

  # GPS (local operator)  
  #  old files were - 06,12,18,24 +/- 3H, and named obs_seq.gpsroYYYYMMMDDHH
  #  new files are 0Z to 0Z, daily, and named obs_seq.gpsro_YYYYMMDD
  #  for 0Z files, start time will be previous day.
  if [[ $ssymd != $ymd ]]; then
    echo /glade/p/cisl/dares/Observations/GPS/local-allsats-2019/${ssym}/obs_seq.gpsro_${ssymd} >> olist
  fi
  echo /glade/p/cisl/dares/Observations/GPS/local-allsats-2019/${ym}/obs_seq.gpsro_${ymd} >> olist


#  ## quikscat - named by day-of-year, 0Z to 0Z. 
#  # for 0Z files, start time will be previous day.
#  if [[ $ssymd != $ymd ]]; then
#    echo /glade/p/cisl/dares/Observations/obsolete/QuikSCAT_24_subx2_ascii/${ssym}/qscatL2B_${ssyear}_${ssdoy}_obs_seq.out >> olist
#  fi
#  echo /glade/p/cisl/dares/Observations/obsolete/QuikSCAT_24_subx2_ascii/${ym}/qscatL2B_${cyear}_${cdoy}_obs_seq.out >> olist

#  ## FIXME!  assumes daily files, but you can copy the 6H section if you want
#  ## to start there as well.
#  ## MOPITT - fill this in as you need with the right dir name and file format.
#  if [[ $ssymd != $ymd ]]; then
#    echo /glade/../MOPITT/${ssym}/obs_seq.${ssymd}.out >> olist
#  fi
#  echo /glade/../MOPITT/${ym}/obs_seq.${ymd}.out >> olist

  ## AIRS - named by YYYYMMDD, 0Z to 0Z (granule 240 has scans from next day)
  # for 0Z files, start time will be previous day.  Last file of a day 
  # (granule 240) contains data for next day, but this gets taken care 
  # of since we are already including the previous day's data for 0Z.
  if [[ $ssymd != $ymd ]]; then
    echo /glade/p/cisl/dares/Observations/AIRS/AIRS_24_sub9x10_nceperrs/${ssym}/obs_seq.AIRS.${ssymd}.out >> olist
  fi
  echo /glade/p/cisl/dares/Observations/AIRS/AIRS_24_sub9x10_nceperrs/${ym}/obs_seq.AIRS.${ymd}.out >> olist
 

  # make sure the output dir exists before running merge.
  # if you want to change the naming convention for the output location,
  # change 'outdir' here and it will get changed in the template file.

  # FIX THIS TO YOUR OWN LOCATION!
  outdir=/glade/../where_you_want/${ym}_3H_CESM
  if [ ! -d $outdir ]; then mkdir -p $outdir; fi

  sed -e "s/YYYY/${cyear}/g"       \
      -e "s/MM/${cmonth}/g"        \
      -e "s/DD/${cday}/g"          \
      -e "s/HH/${chour}/g"         \
      -e "s/GREG1/${ssgregday}/g"  \
      -e "s/GREG2/${ssgregsec}/g"  \
      -e "s/GREG3/${segregday}/g"  \
      -e "s/GREG4/${segregsec}/g"  \
      -e "s;OUTDIR;${outdir};g"    \
      -e "s/OUTTIME/${cesmtime}/g"  < ./input.nml.template > input.nml

  # FIXME:
  # actually run the tool, or just debug script 
  #  change "doit" to "" to skip; change "" to "doit" to make files
  if [ "" ] ; then
  
   ./obs_sequence_tool
  
  else
  
   echo current time is: $curstep
   echo input files are:
   cat olist
   echo input.nml is:
   fgrep obs_sequence_tool_nml -A 7 input.nml
   echo would be calling ./obs_sequence_tool here
  
  fi  

  # ----------------------------------------

  # Make the output files read-only so they don't get overwritten by accident.
  chmod 444 $outdir/*


  # advance the step; the output is YYYYMMDDHH
  curstep=$nextstep
  nextstep=`echo $curstep     $step | ./advance_time`
 stepstart=`echo $curstep $halfback | ./advance_time`
   stepend=`echo $curstep $halfstep | ./advance_time`


  # advance the loop counter
  let s=$s+1
 
done

echo job finished.

exit 0

