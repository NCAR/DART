#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# THIS VERSION IS USEFUL FOR FLUX OBSERVATIONS FOR CLM.
# split file(s) into "daily" files which start at 00:00Z
# the previous day and end the same (previous) day at 23:59Z.
# The date in the filename is the date/time at which CLM stops.
# The CLM history file has fluxes for the PREVIOUS 24 hours.

# -----------------------------------------------------------------------------
# set the first and last days to be split.
# depending on the window and the input file,
# the data from outside these bounds may be needed.

start_year=2004
start_month=01
start_day=01

end_year=2004
end_month=12
end_day=31

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
end_d=(`echo ${end_year}${mon2}${day2}00 0 -g | ../work/advance_time`)
echo "last  day,seconds is day ${end_d[0]} ${end_d[1]}"

mon2=`printf %02d $month1`
day2=`printf %02d $day1`
start_d=(`echo ${start_year}${mon2}${day2}00 0 -g | ../work/advance_time`)
echo "first day,seconds is day ${start_d[0]} ${start_d[1]}"

# how many total days are going to be split (for the loop counter)
let totaldays=${end_d[0]}-${start_d[0]}+1

# -----------------------------------------------------------------------------
# form some strings for logging.
# time_one    .... the first time in the file
# time_end    .... the last  time in the file
# filetime    .... the time in the file NAME ... usually the center.
# with no -g option advance_time returns strings in the format YYYYMMDDHH

time_one=(`echo ${start_year}${mon2}${day2}00 -1d | ../work/advance_time`)
time_end=(`echo ${start_year}${mon2}${day2}00 -1s | ../work/advance_time`)
filetime=(`echo ${start_year}${mon2}${day2}00   0 | ../work/advance_time`)

# loop over each day
let d=1
while (( d <= totaldays)) ; do

   echo "subsetting $d of $totaldays ..."
   #echo $filetime $time_end

   # string for first time in the file
    pyear=${time_one:0:4}
   pmonth=${time_one:4:2}
     pday=${time_one:6:2}
    phour=${time_one:8:2}

   # string for last time in the file
    nyear=${time_end:0:4}
   nmonth=${time_end:4:2}
     nday=${time_end:6:2}
    nhour=${time_end:8:2}

   # string for time for the file NAME
    fyear=${filetime:0:4}
   fmonth=${filetime:4:2}
     fday=${filetime:6:2}
    fhour=${filetime:8:2}

   # compute the equivalent DART timestamps here - seconds and days.
   g=(`echo ${fyear}${fmonth}${fday}${fhour} -1d -g | ../work/advance_time`)
   dart1d=${g[0]}
   darts1=${g[1]}
#  dart1s=${darts1}+1     for non-flux tower observations
   dart1s=${darts1}

   # flux tower observations do not actually contain the time in the file name.
   g=(`echo ${fyear}${fmonth}${fday}${fhour} -1s -g | ../work/advance_time`)
   dartNd=${g[0]}
   dartNs=${g[1]}

   g=(`echo ${fyear}${fmonth}${fday}${fhour}   0 -g | ../work/advance_time`)
   dartFd=${g[0]}
   dartFs=`printf %05d ${g[1]}`

   echo "first $pyear $pmonth $pday $phour which is dart $dart1d $dart1s"
   echo "last  $nyear $nmonth $nday $nhour which is dart $dartNd $dartNs"
   echo "file  $fyear $fmonth $fday $fhour which is dart $dartFd $dartFs"

   # Create the file containing the input observation sequences.
   # The input.nml.template explicitly references this file.
   # The input observation sequences generally live in directories
   # with names YYYYMM. FIXME ... year boundaries ...

   echo "/glade/p/image/Observations/FluxTower/obs_seq.USBar.2004" >  olist
   echo "/glade/p/image/Observations/FluxTower/obs_seq.USHa1.2004" >> olist
   echo "/glade/p/image/Observations/FluxTower/obs_seq.USNR1.2004" >> olist
   echo "/glade/p/image/Observations/FluxTower/obs_seq.USSP3.2004" >> olist
   echo "/glade/p/image/Observations/FluxTower/obs_seq.USSRM.2004" >> olist
   echo "/glade/p/image/Observations/FluxTower/obs_seq.USWcr.2004" >> olist
   echo "/glade/p/image/Observations/FluxTower/obs_seq.USWrc.2004" >> olist

   # make sure output dir exists
   OUTDIR=../${fyear}${fmonth}
   if [[ ! -d ${OUTDIR} ]] ; then
        mkdir ${OUTDIR}
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
   time_one=(`echo ${pyear}${pmonth}${pday}${phour} +1d | ../work/advance_time`)
   filetime=(`echo ${fyear}${fmonth}${fday}${fhour} +1d | ../work/advance_time`)
   time_end=(`echo ${nyear}${nmonth}${nday}${nhour} +1d | ../work/advance_time`)
   echo "next set of times are: $time_one $filetime $time_end"

   # advance the loop counter
   let d=d+1

done

exit 0

