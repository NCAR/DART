#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# This version splits files into chunks that are compatible with the
# Flux Tower forward operator, which relies on using the CLM *.h1.* 
# history file with timestamps that do not actually include the time
# used in the filename.  
# 
# set the first and last days to be split.  can roll over
# month and year boundaries now! note that for the last day
# the data from the following day must be available.

start_year=2004
start_month=1
start_day=1

end_year=2004
end_month=1
end_day=2

# end of things you should have to set in this script

# convert the start and stop times to DART days/seconds, so we can
# compute total number of days including rolling over month and
# year boundaries.
# are < 10.  do the end time first so we can use the same values
# to set the initial day while we are doing the total day calc.

# make sure there is an initial input.nml for advance_time
cp -f ../work/input.nml.template input.nml || exit -1

# advance_time (with the -g flag) outputs 2 integers:
# days seconds
year1=`echo  $start_year  | bc`
month1=`echo $start_month | bc`
day1=`echo   $start_day   | bc`
dart_1=`printf %04d%02d%02d%02d $year1 $month1 $day1 0`

year2=`echo  $end_year  | bc`
month2=`echo $end_month | bc`
day2=`echo   $end_day   | bc`
dart_N=`printf %04d%02d%02d%02d $year2 $month2 $day2 0`

T1=(`echo ${dart_1} -6h -g | ../work/advance_time`)
TN=(`echo ${dart_N}  0  -g | ../work/advance_time`)

echo "start of period of interest is $dart_1 and ${T1[0]} ${T1[1]}"
echo "end   of period of interest is $dart_N and ${TN[0]} ${TN[1]}"

# how many total days are going to be split (for the loop counter)
# and we are going to make 4 files per day, so times 4.  if we change
# the increment (6h currently) the multiple will have to change as well.
let totaldays=${TN[0]}-${T1[0]}+1
let totaltimes=${totaldays}*4

# these are a string in the format YYYYMMDDHH
# do them here to prime the loop below which first takes them apart.
prevtime=(`echo ${dart_1} -6h | ../work/advance_time`)
middtime=(`echo ${dart_1}   0 | ../work/advance_time`)
nexttime=(`echo ${dart_1}   0 | ../work/advance_time`)

# loop over each new output file.

let itime=1
while ((itime <= totaltimes)) ; do

   echo ''
   # The FluxTower observations are one-sided ...
   # 6 hours before the time of the filename.

   # parse out the parts from a string which is YYYYMMDDHH
   # both for the current time, the end time, and the middle of
   # the interval (the middle is for the filename; the others
   # for the obs time limits.)

    pyear=${prevtime:0:4}
   pmonth=${prevtime:4:2}
     pday=${prevtime:6:2}
    phour=${prevtime:8:2}

    myear=${middtime:0:4}
   mmonth=${middtime:4:2}
     mday=${middtime:6:2}
    mhour=${middtime:8:2}

    nyear=${nexttime:0:4}
   nmonth=${nexttime:4:2}
     nday=${nexttime:6:2}
    nhour=${nexttime:8:2}

   # echo prev $pyear $pmonth $pday $phour
   # echo midd $myear $mmonth $mday $mhour
   # echo next $nyear $nmonth $nday $nhour

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

   # compute the equivalent DART days here for the current time
   # plus 1 second (the interval start) and the end time.

   g=(`echo ${myear}${mmonth}${mday}${mhour} -6h -g | ../work/advance_time`)
   dart1d=${g[0]}
   darts1=${g[1]}
   let dart1s=${darts1}+1

   g=(`echo ${myear}${mmonth}${mday}${mhour}   0 -g | ../work/advance_time`)
   dartMd=${g[0]}
   dartMs=`printf %05d ${g[1]}`

   g=(`echo ${myear}${mmonth}${mday}${mhour}   0 -g | ../work/advance_time`)
   dartNd=${g[0]}
   dartNs=`printf %05d ${g[1]}`

   echo prev $pyear $pmonth $pday $phour which is dart $dart1d $dart1s
   echo midd $myear $mmonth $mday $mhour which is dart $dartMd $dartMs
   echo next $nyear $nmonth $nday $nhour which is dart $dartNd $dartNs

   echo "starting    time" ${pyear}${pmonth}${pday}${phour} plus 1 sec
   echo "output      file" ${myear}${mmonth}${mday}${mhour}
   echo "ending      time" ${nyear}${nmonth}${nday}${nhour}

   # make sure output dir exists
   OUTDIR=../${myear}${mmonth}_6H
   if [[ ! -d ${OUTDIR} ]] ; then
        mkdir ${OUTDIR}
   fi

   sed -e "s#OUTDIR#${OUTDIR}#g" \
   sed -e "s#YYYY#${myear}#g"    \
       -e "s#MM#${mmonth}#g"     \
       -e "s#DD#${mday}#g"       \
       -e "s#SSSSS#${dartMs}#g"  \
       -e "s#DART1D#${dart1d}#g" \
       -e "s#DART1S#${dart1s}#g" \
       -e "s#DARTND#${dartNd}#g" \
       -e "s#DARTNS#${dartNs}#g" < ../work/input.nml.template > input.nml

   # do the extract here
   #echo calling ./obs_sequence_tool here
   ../work/obs_sequence_tool

   # advance 6 hours; the output is YYYYMMDDHH
   prevtime=(`echo ${pyear}${pmonth}${pday}${phour} +6h | ../work/advance_time`)
   middtime=(`echo ${myear}${mmonth}${mday}${mhour} +6h | ../work/advance_time`)
   nexttime=(`echo ${nyear}${nmonth}${nday}${nhour} +6h | ../work/advance_time`)
   echo $prevtime $middtime $nexttime

   # advance the loop counter
   let itime=itime+1
   echo "itime = $itime"

done

exit 0

