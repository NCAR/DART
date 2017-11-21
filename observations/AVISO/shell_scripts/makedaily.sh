#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Create a series of "daily" files which start at 12:01Z the previous day and 
# end at 12:00Z on the day that matches the day in the filename.
# The assumption is that the the list of input files for each daily file may
# contain observation sequence files with predictable file names. Simply
# specify the start and end dates, a geographic region of interest, and
# (for the AVISO obs converter
#
#BSUB -J splitobs
#BSUB -W 12:00
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -q geyser
#BSUB -P CMIA0001
#BSUB -n 1
#
#PBS -N splitobs      
#PBS -o splitobs.out
#PBS -e splitobs.err
#PBS -l walltime=00:00:10
#PBS -q economy 
#PBS -l select=1:ncpus=1
#PBS -A P86850054 

# account for the fact that some PBS/SLURM implementations need
# to be told where to run.
if [[ -n ${!PBS_O_WORKDIR} ]] ; then
  cd ${PBS_O_WORKDIR}
fi

# set the first and last days to be split.  can roll over
# month and year boundaries now!   note that for the first day
# you need the previous day data available.

let start_year=2015
let start_month=7
let start_day=1

let end_year=2015
let end_month=7
let end_day=31

EXEDIR='../work'
DATDIR='/glade2/scratch2/fredc/Observations/AVISO/data'
OUTDIR="/glade/scratch/${USER}/Observations/AVISO"

LONGITUDE_WEST=0.0
LONGITUDE_EAST=360.0
LATITUDE_SOUTH=-90.0
LATITUDE_NORTH=90.0

# end of things you should have to set in this script

# make sure there is an initial input.nml for advance_time                                            
# input.nml gets overwritten in the subsequent loop.                                                  
cp -f input.nml.template input.nml || exit -1

# convert the start and stop times to gregorian days, so we can
# compute total number of days including rolling over month and
# year boundaries.  make sure all values have leading 0s if they
# are < 10.  do the end time first so we can use the same values
# to set the initial day while we are doing the total day calc.

# these outputs from advance time (with the -g flag) are
# 2 integers: gregorian_day_number seconds
# and since we don't set hours, minutes, or seconds, the second
# number is always 0 and uninteresting for us.
mon2=`printf %02d $end_month`
day2=`printf %02d $end_day`
echo ${end_year}${mon2}${day2}00
end_d=(`echo ${end_year}${mon2}${day2}00 0 -g | ${EXEDIR}/advance_time`)
echo $end_d

mon2=`printf %02d $start_month`
day2=`printf %02d $start_day`
echo ${start_year}${mon2}${day2}00
start_d=(`echo ${start_year}${mon2}${day2}00 0 -g | ${EXEDIR}/advance_time`)
echo $start_d

# these are a string in the format YYYYMMDDHH
# do them here to prime the loop below which first takes them apart.
prevday=(`echo ${start_year}${mon2}${day2}00 -1d | ${EXEDIR}/advance_time`)
currday=(`echo ${start_year}${mon2}${day2}00   0 | ${EXEDIR}/advance_time`)

# how many total days are going to be merged (for the loop counter)
# (pull out the first of the 2 numbers which are output from advance_time)
let totaldays=${end_d}-${start_d}+1
echo $totaldays

# loop over each day
let d=1
while (( d <= totaldays)) ; do

  # parse out the parts from a string which is YYYYMMDDHH
  # both for the current day and tomorrow
  pyear=${prevday:0:4}
  pmonth=${prevday:4:2}
  pday=${prevday:6:2}

  cyear=${currday:0:4}
  cmonth=${currday:4:2}
  cday=${currday:6:2}

  echo ' '
  echo day $d ... prev $pyear $pmonth $pday
  echo day $d ... curr $cyear $cmonth $cday

  # compute the equivalent gregorian days here.
  g=(`echo ${cyear}${cmonth}${cday}00 -12h -g | ${EXEDIR}/advance_time`)
  greg0=${g[0]}
  g=(`echo ${cyear}${cmonth}${cday}00 +12h -g | ${EXEDIR}/advance_time`)
  greg1=${g[0]}

  echo "start = $greg0"
  echo "end   = $greg1"

  echo starting AVISO obs for ${cyear}${cmonth}${cday} gregorian= $greg1

  # make sure we start with an empty list - no other days/files to process.
  rm -f olist input.nml

  # These are the character strings defining the satellites in the filenames in use by AVISO
  for sat in al c2 e1 e2 en enn g2 h2 j1 j1g j1n j2 tp tpn
  do

     # To ensure that we have all the observations, we need the previous day and the current day.
     # Each 'daily' output file has observations from 12:00:01 Z the previous day to 12:00:00 Z of the current day.
     # The time in the output file name is the 00:00:00 Z of the current day.

     if [ -f ${DATDIR}/obs_seq.${sat}.${pyear}${pmonth}${pday} ] && [ -f ${DATDIR}/obs_seq.${sat}.${cyear}${cmonth}${cday} ]; then
        echo "${DATDIR}/obs_seq.${sat}.${pyear}${pmonth}${pday}" >> olist
        echo "${DATDIR}/obs_seq.${sat}.${cyear}${cmonth}${cday}" >> olist
     fi

  done

  # There are also observations in the world ocean database we might want.
  # Use them if they exist.

  WOD_PDAY="/glade/p/image/Observations/WOD09/${pyear}${pmonth}/obs_seq.0Z.${pyear}${pmonth}${pday}"
  WOD_CDAY="/glade/p/image/Observations/WOD09/${cyear}${cmonth}/obs_seq.0Z.${cyear}${cmonth}${cday}"

  if [ -f  ${WOD_PDAY} ]; then
     echo "${WOD_PDAY}" >> olist
  fi
  if [ -f  ${WOD_CDAY} ]; then
     echo "${WOD_CDAY}" >> olist
  fi

  # update run-time input

  FILENAMEOUT=obs_seq.0Z.${cyear}${cmonth}${cday}

  sed -e "s/<FILENAMEOUT>/${FILENAMEOUT}/g"       \
      -e "s/<LONGITUDE_WEST>/${LONGITUDE_WEST}/g" \
      -e "s/<LONGITUDE_EAST>/${LONGITUDE_EAST}/g" \
      -e "s/<LATITUDE_SOUTH>/${LATITUDE_SOUTH}/g" \
      -e "s/<LATITUDE_NORTH>/${LATITUDE_NORTH}/g" \
      -e "s/<FIRSTDAY>/${greg0}/g"  \
      -e "s/<LASTDAY>/${greg1}/g"  < input.nml.template > input.nml

  # make sure output dir exists
  OBSDIR=${OUTDIR}/${cyear}${cmonth}
  if [[ ! -d $OBSDIR ]] ; then
     mkdir -p $OBSDIR || exit -2
  fi

  # do the extract here
  ${EXEDIR}/obs_sequence_tool || exit -3

  \mv -v ${FILENAMEOUT} ${OBSDIR}

  # advance the day; the output is YYYYMMDD00
  prevday=$currday
  currday=(`echo ${cyear}${cmonth}${cday}00 +1d | ${EXEDIR}/advance_time`)

  # advance the loop counter
  let d=d+1
echo d=$d
 
done

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

