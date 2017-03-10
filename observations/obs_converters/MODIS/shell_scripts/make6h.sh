#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# Split a long file into "6hourly" files which start 6 hours before the
# time in the filename and end exactly one second before the day in the
# filename. This is only logical if you understand the CLM obs operator.
#
# The only flux tower forward observation operator is for CLM,
# which simply grabs values from a history file.
# The CLM h1 files contain everything STARTING with the time in their filename.
# all output intervals from that run are simply appended to that file.
# Consequently, we need to know the filename from the START of the model advance
# that resulted in the current model state.
#
#  | start of the model advance ... *.h1.* file starts getting written
#  |
#  X==X==X==X==X==X==X==X==X==X==X==[O] (CLM model advance)
#  |<------- hist_nhtfrq ------->|   |
#                                |   | END of the model advance
#                                |
#                                | END of the data in the *.h1.* file

# -----------------------------------------------------------------------------
# set the first and last days to be split.
# depending on the window and the input file,
# the data from outside these bounds may be needed.

start_year=2004
start_month=1
start_day=1

end_year=2004
end_month=12
end_day=31

# end of things you should have to set in this script IFF you are
# content to have 6hourly files with observations from the 6 hour
# period BEFORE the date in the filename.

# -----------------------------------------------------------------------------
# convert the start and stop times to DART days/seconds, so we can compute
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
cp -f ../work/input.nml.template input.nml || exit -1

# these outputs from advance time (with the -g flag) are 2 integers:
# gregorian_day_number seconds
# we can hardwire hours, minutes, and seconds
mon2=`printf %02d $month2`
day2=`printf %02d $day2`
end_d=(`echo     ${end_year}${mon2}${day2}00 0 -g | ../work/advance_time`)
echo "last  ${year2} ${month2} ${day2} is DART (days,seconds) ${end_d[0]} ${end_d[1]}"

mon2=`printf %02d $month1`
day2=`printf %02d $day1`
start_d=(`echo ${start_year}${mon2}${day2}00 0 -g | ../work/advance_time`)
echo "first ${year1} ${month1} ${day1} is DART (days,seconds) ${start_d[0]} ${start_d[1]}"

# how many total days are going to be split (for the loop counter)
# and we are going to make 4 files per day, so times 4.  if we change
# the increment (6h currently) the multiple will have to change as well.
let totaldays=${end_d[0]}-${start_d[0]}+1
let totaltimes=${totaldays}*4

# -----------------------------------------------------------------------------
# form some strings for logging.
# time_one    .... the first time in the file
# time_end    .... the last  time in the file
# filetime    .... the time in the file NAME ... usually the center.
# with no -g option advance_time returns strings in the format YYYYMMDDHH

time_one=(`echo ${start_year}${mon2}${day2}00 -6h | ../work/advance_time`)
time_end=(`echo ${start_year}${mon2}${day2}00 -1s | ../work/advance_time`)
filetime=(`echo ${start_year}${mon2}${day2}00   0 | ../work/advance_time`)

# loop over each output file time.

let itime=1
while ((itime <= totaltimes)) ; do

   echo ''
   echo "subsetting file $itime of $totaltimes"
   # The FluxTower observations are one-sided ...

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
   # first time in file ...
   g=(`echo ${fyear}${fmonth}${fday}${fhour} -6h -g | ../work/advance_time`)
   dart1d=${g[0]}
   dart1s=${g[1]}

   # last time in file ...
   g=(`echo ${fyear}${fmonth}${fday}${fhour} -1s -g | ../work/advance_time`)
   dartNd=${g[0]}
   dartNs=${g[1]}

   # time in file name ... must have all the zeros for seconds
   g=(`echo ${fyear}${fmonth}${fday}${fhour}   0 -g | ../work/advance_time`)
   dartFd=${g[0]}
   dartFs=`printf %05d ${g[1]}`

   echo "first $year1 $month1 $day1 $hour1 which is dart $dart1d $dart1s"
   echo "last  $yearN $monthN $dayN $hourN which is dart $dartNd $dartNs"
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
   OUTDIR=${fyear}${fmonth}_6H
   if [[ ! -d ${OUTDIR} ]] ; then
        mkdir ${OUTDIR}
   fi

   sed -e "s/OUTDIR/${OUTDIR}/g" \
       -e "s/YYYY/${fyear}/g"    \
       -e "s/MM/${fmonth}/g"     \
       -e "s/DD/${fday}/g"       \
       -e "s/SSSSS/${dartFs}/g"  \
       -e "s/DART1D/${dart1d}/g" \
       -e "s/DART1S/${dart1s}/g" \
       -e "s/DARTND/${dartNd}/g" \
       -e "s/DARTNS/${dartNs}/g" < ../work/input.nml.template > input.nml

   # do the extract here
   ../work/obs_sequence_tool

   # advance 6 hours; the output is YYYYMMDDHH
   time_one=(`echo ${year1}${month1}${day1}${hour1} +6h | ../work/advance_time`)
   time_end=(`echo ${yearN}${monthN}${dayN}${hourN} +6h | ../work/advance_time`)
   filetime=(`echo ${fyear}${fmonth}${fday}${fhour} +6h | ../work/advance_time`)
   # echo "next set of times are (first,last,file): $time_one $time_end $filetime"

   # advance the loop counter
   let itime=itime+1

done

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

