#!/bin/bash
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# THIS VERSION IS USEFUL FOR FLUX OBSERVATIONS FOR CLM.
# split file(s) into "daily" files which start at 00:00Z
# the previous day and end the same (previous) day at 23:59Z.
# The date in the filename is the date/time at which CLM stops.
# The CLM history file has fluxes for the PREVIOUS 24 hours.

# -----------------------------------------------------------------------------
# set the first and last days to be split.
# depending on the window and the input file,
# the data from outside these bounds may be needed.

start_year=2011
start_month=01
start_day=05

end_year=2011
end_month=1
end_day=08

SOURCE_OBS_BASE=/Users/thoar/Desktop/EASE_Grid/2011_north/DART_OBS_SEQ
OUTPUT_OBS_BASE=/Users/thoar/Desktop/EASE_Grid/2011_north/DART_OBS_CLM

# end of things you should have to set in this script IFF you are
# content to have 'daily' files with observations 
# date in the filename.

# -----------------------------------------------------------------------------
# convert the start and stop times to gregorian days, so we can compute
# total number of days including rolling over month and year boundaries.
# do the end time first so we can use the same values to set the
# initial day while we are doing the total day calculation.

# bc strips off any preceeding zeros to prevent interpretation as octal.
year1=`echo  $start_year  | bc`
month1=`echo $start_month | bc`
day1=`echo   $start_day   | bc`

year2=`echo  $end_year    | bc`
month2=`echo $end_month   | bc`
day2=`echo   $end_day     | bc`

# make sure there is an initial input.nml for advance_time
# input.nml gets overwritten in the subsequent loop.
cp -f  ../work/input.nml.template input.nml || exit -1

# these outputs from advance time (with the -g flag) are 2 integers:
# gregorian_day_number seconds
# since we are concerned with daily files at 00Z,
# we can hardwire hours, minutes, and seconds
mon2=`printf %02d $month2`
day2=`printf %02d $day2`
end_d=(`echo     ${end_year}${mon2}${day2}00 0 -g | ../work/advance_time`)
echo "last  day,seconds is ${end_d[0]} ${end_d[1]}"

mon2=`printf %02d $month1`
day2=`printf %02d $day1`
start_d=(`echo ${start_year}${mon2}${day2}00 0 -g | ../work/advance_time`)
echo "first day,seconds is ${start_d[0]} ${start_d[1]}"

# how many total days are going to be split (for the loop counter)
let totaldays=${end_d[0]}-${start_d[0]}+1

# -----------------------------------------------------------------------------
# form some strings for logging.
# time_one    .... the first time in the file
# time_end    .... the last  time in the file
# filetime    .... the time in the file NAME ... usually the center.
# with no -g option advance_time returns strings in the format YYYYMMDDHH

time_one=(`echo ${start_year}${mon2}${day2}00 -12h | ../work/advance_time`)
time_end=(`echo ${start_year}${mon2}${day2}00 +12h | ../work/advance_time`)
filetime=(`echo ${start_year}${mon2}${day2}00    0 | ../work/advance_time`)

# loop over each output file time (day)

let d=1
while (( d <= totaldays)) ; do

   echo "subsetting file $d of $totaldays ..."
   #echo $filetime $time_end

   # string for first time in the file
    year1=${time_one:0:4}
   month1=${time_one:4:2}
     day1=${time_one:6:2}
    hour1=${time_one:8:2}

   # string for last time in the file
    yearN=${time_end:0:4}
   monthN=${time_end:4:2}
     dayN=${time_end:6:2}
    hourN=${time_end:8:2}

   # string for time for the file NAME
    fyear=${filetime:0:4}
   fmonth=${filetime:4:2}
     fday=${filetime:6:2}
    fhour=${filetime:8:2}

   # compute the equivalent DART timestamps here - seconds and days.
   # first time in file ... usually 1 second after the preceeding end time.
   g=(`echo ${fyear}${fmonth}${fday}${fhour} -12h -g | ../work/advance_time`)
   dart1d=${g[0]}
   let dart1s=${g[1]}+1

   # last time in file ... this time IS included in the file.
   g=(`echo ${fyear}${fmonth}${fday}${fhour} +12h -g | ../work/advance_time`)
   dartNd=${g[0]}
   dartNs=${g[1]}

   # time in file name ... must have all the zeros for seconds
   g=(`echo ${fyear}${fmonth}${fday}${fhour}    0 -g | ../work/advance_time`)
   dartFd=${g[0]}
   dartFs=`printf %05d ${g[1]}`

   echo "first $year1 $month1 $day1 $hour1 which is dart $dart1d $dart1s"
   echo "file  $fyear $fmonth $fday $fhour which is dart $dartFd $dartFs"
   echo "last  $yearN $monthN $dayN $hourN which is dart $dartNd $dartNs"

   # Create the file containing the input observation sequences.
   # The input.nml.template explicitly references this file.
   # The input observation sequences generally live in directories
   # with names YYYYMM.

   echo "${SOURCE_OBS_BASE}/${year1}${month1}/obs_seq.${year1}-${month1}-${day1}-*.out"
   echo "${SOURCE_OBS_BASE}/${yearN}${monthN}/obs_seq.${yearN}-${monthN}-${dayN}-*.out"

   ls -1 ${SOURCE_OBS_BASE}/${year1}${month1}/obs_seq.${year1}-${month1}-${day1}-*.out >  olist
   ls -1 ${SOURCE_OBS_BASE}/${yearN}${monthN}/obs_seq.${yearN}-${monthN}-${dayN}-*.out >> olist

   # make sure output dir exists
   OUTDIR=${OUTPUT_OBS_BASE}/${fyear}${fmonth}
   if [[ ! -d ${OUTDIR} ]] ; then
        mkdir -p ${OUTDIR}
   fi

   sed -e "s#OUTDIR#${OUTDIR}#g" \
       -e "s#YYYY#${fyear}#g"    \
       -e "s#MM#${fmonth}#g"     \
       -e "s#DD#${fday}#g"       \
       -e "s#SSSSS#${dartFs}#g"  \
       -e "s#DART1D#${dart1d}#g" \
       -e "s#DART1S#${dart1s}#g" \
       -e "s#DARTND#${dartNd}#g" \
       -e "s#DARTNS#${dartNs}#g" < ../work/input.nml.template > input.nml

   # do the extract here
   ../work/obs_sequence_tool

   # advance the day; the output is YYYYMMDD00
   time_one=(`echo ${year1}${month1}${day1}${hour1} +1d | ../work/advance_time`)
   time_end=(`echo ${yearN}${monthN}${dayN}${hourN} +1d | ../work/advance_time`)
   filetime=(`echo ${fyear}${fmonth}${fday}${fhour} +1d | ../work/advance_time`)
     echo "next set of times are: $time_one $time_end $filetime"

   # advance the loop counter
   let d=d+1

done

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

