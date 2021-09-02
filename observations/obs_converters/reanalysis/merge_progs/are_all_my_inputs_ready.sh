#!/bin/ksh

# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
#--------------------------------------------------------------
# DESCRIPTION:
# validate the input files before you try to
# Merge NCEP BUFR obs (inc ACARS), GPS RO obs, and AIRS T/Q obs.
#
# this script does the work to compute which files are needed 
# for the merge period and lists any which aren't available.
#
# run mergeit6h.sh afterwards.
#
#--------------------------------------------------------------

# set the first and last days/hours.  can roll over month and year boundaries.
# if start_hour is 12Z and step is 6H, then the first file will be centered 
# at 12Z on this day, with data starting at 9Z and ending at 15Z.
#
# set these WITHOUT leading 0s.  they are numeric values, not strings.
#
# if you set the start hour to 0, you must have data files for
# the day immediately preceeding this date.


let start_year=2016
let start_month=1
let start_day=1
let start_hour=0

let end_year=2016
let end_month=12
let end_day=31
let end_hour=18

# string to control the number of hours in each file.
# set the forward and backwards half steps to be consistent.
# (the +1s will be done when setting the namelist seconds)
step='+6h'
halfstep='+3h'
halfback='-3h'

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
endstring=${end_year}${mon2c}${day2c}${hor2c}
set -A end_t `echo $endstring 0 -g | ./advance_time`

mon2c=`printf %02d $start_month`
day2c=`printf %02d $start_day`
hor2c=`printf %02d $start_hour`
startstring=${start_year}${mon2c}${day2c}${hor2c}
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

echo This script is going to validate inputs for $totalsteps DART obs sequence files.
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

  # parse out the parts from a string which is YYYYMMDDHH
  # use cut with the byte option to pull out columns 1-4, 5-6, 7-8, and 9-10
  # c = current middle-of-step, n = next middle, ss = step start, se = step end

   cyear=`echo $curstep | cut -b1-4`
  cmonth=`echo $curstep | cut -b5-6`
    cday=`echo $curstep | cut -b7-8`
   chour=`echo $curstep | cut -b9-10`

   nyear=`echo $nextstep | cut -b1-4`
  nmonth=`echo $nextstep | cut -b5-6`
    nday=`echo $nextstep | cut -b7-8`
   nhour=`echo $nextstep | cut -b9-10`

  ssyear=`echo $stepstart | cut -b1-4`
 ssmonth=`echo $stepstart | cut -b5-6`
   ssday=`echo $stepstart | cut -b7-8`
  sshour=`echo $stepstart | cut -b9-10`

  seyear=`echo $stepend | cut -b1-4`
 semonth=`echo $stepend | cut -b5-6`
   seday=`echo $stepend | cut -b7-8`
  sehour=`echo $stepend | cut -b9-10`

  # compute the equivalent gregorian days/secs here.
  set -A g `echo ${cyear}${cmonth}${cday}${chour} 0 -g | ./advance_time`
  cgregday=${g[0]}; cgregsec=${g[1]}

  # special for dart: step start needs to be +1 second 
  set -A g `echo ${ssyear}${ssmonth}${ssday}${sshour} +1s -g | ./advance_time`
  ssgregday=${g[0]}; ssgregsec=${g[1]}

  set -A g `echo ${seyear}${semonth}${seday}${sehour} 0 -g | ./advance_time`
  segregday=${g[0]}; segregsec=${g[1]}

  # compute the CESM-style time string for the output filename
  cesmtime=`echo ${cyear}${cmonth}${cday}${chour} 0 -c | ./advance_time`

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
  #echo "starting processing for "  ${cyear}  ${cmonth}  ${cday}  ${chour}, gregorian day/sec:  $cgregday  $cgregsec #, DOY $cdoy 
  #echo "          step start is " ${ssyear} ${ssmonth} ${ssday} ${sshour}, gregorian day/sec: $ssgregday $ssgregsec #, DOY $ssdoy 
  #echo "            step end is " ${seyear} ${semonth} ${seday} ${sehour}, gregorian day/sec: $segregday $segregsec #, DOY $sedoy 
  #echo ''


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
  # get the right filename.  since these are already 6h files
  # we will only every need one.
  ls /glade/p/cisl/dares/Observations/NCEP+ACARS/${ym}_6H/obs_seq${ymdh} >> olist

  # GPS (local operator)  
  #  old files were - 06,12,18,24 +/- 3H, and named obs_seq.gpsroYYYYMMMDDHH
  #  new files are 0Z to 0Z, daily, and named obs_seq.gpsro_YYYYMMDD
  #  for 0Z files, start time will be previous day.
  ls /glade/p/cisl/dares/Observations/GPS/local-allsats-2019/${ym}/obs_seq.gpsro_${ymd} >> olist
  if [[ $ssymd != $ymd ]]; then
    ls /glade/p/cisl/dares/Observations/GPS/local-allsats-2019/${ssym}/obs_seq.gpsro_${ssymd} >> olist
  fi


  ## quikscat - named by day-of-year, 0Z to 0Z. 
  #ls /glade/p/cisl/dares/Observations/obsolete/QuikSCAT_24_subx2_ascii/${ym}/qscatL2B_${cyear}_${cdoy}_obs_seq.out >> olist
  # for 0Z files, start time will be previous day.
  #if [[ $ssymd != $ymd ]]; then
  #  ls /glade/p/cisl/dares/Observations/obsolete/QuikSCAT_24_subx2_ascii/${ssym}/qscatL2B_${ssyear}_${ssdoy}_obs_seq.out >> olist
  #fi

  ## AIRS - named by YYYYMMDD, 0Z to 0Z (granule 240 has scans from next day)
  # for 0Z files, start time will be previous day.  Last file of a day 
  # (granule 240) contains data for next day, but this gets taken care 
  # of since we are already including the previous day's data for 0Z.
  ls /glade/p/cisl/dares/Observations/AIRS/AIRS_24_sub9x10_nceperrs/${ym}/obs_seq.AIRS.${ymd}.out >> olist
  if [[ $ssymd != $ymd ]]; then
    ls /glade/p/cisl/dares/Observations/AIRS/AIRS_24_sub9x10_nceperrs/${ssym}/obs_seq.AIRS.${ssymd}.out >> olist
  fi


  # validate the input files.  list any that are not present.

  echo current time is: $curstep
  ls -l `cat olist`
  echo ""

  # ----------------------------------------

  # advance the step; the output is YYYYMMDDHH
  curstep=$nextstep
  nextstep=`echo $curstep     $step | ./advance_time`
 stepstart=`echo $curstep $halfback | ./advance_time`
   stepend=`echo $curstep $halfstep | ./advance_time`


  # advance the loop counter
  let s=s+1
 
done

echo job finished.

exit 0

