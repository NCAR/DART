#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# this version gets the tar file from the mass store first.
# unpack one day of tar files at a time, convert them into
# individual obs_seq files.  this program also does the merge
# of the 240 individual daily swaths into a single obs_seq file.
#
# this program should be started from the work directory.
# it assumes ../data, ../tars, the output dir, etc
# exist relative to starting from AIRS/work.

# set the first and last days to be converted.  can roll over
# month and year boundaries now!
let start_year=2006
let start_month=10
let start_day=1

let end_year=2007
let end_month=1
let end_day=31

# relative to work dir
output_dir=../output.thin

# whether to download the tar file from the mass store or not 
# set to one of: true or false
download=true

# end of things you should have to set in this script

# convert the start and stop times to gregorian days, so we can
# compute total number of days including rolling over month and
# year boundaries.  make sure all values have leading 0s if they
# are < 10.  do the end time first so we can use the same values
# to set the initial day while we are doing the total day calc.
mon2=`printf %02d $end_month`
day2=`printf %02d $end_day`
end_d=(`echo ${end_year}${mon2}${day2}00 0 -g | ./advance_time`)

mon2=`printf %02d $start_month`
day2=`printf %02d $start_day`
start_d=(`echo ${start_year}${mon2}${day2}00 0 -g | ./advance_time`)

curday=(`echo ${start_year}${mon2}${day2}00 0 | ./advance_time`)

# how many total days are going to be converted (for the loop counter)
let totaldays=${end_d[0]}-${start_d[0]}+1

# loop over each day
let d=1
while (( d <= totaldays)) ; do

  # parse out the parts from a string which is YYYYMMDDHH
  year=${curday:0:4}
  month=${curday:4:2}
  day=${curday:6:2}

  # compute the equivalent gregorian day here.
  g=(`echo ${year}${month}${day}00 0 -g | ./advance_time`)
  greg=${g[0]}

  echo starting AIRS to obs ${year}${month}${day}
  echo gregorian: $greg

  # download the tar file from the hpss first
  if [[ "$download" = "true" ]]; then
    echo getting ${year}${month}${day}.tar from mass store
    (cd ../tars; hsi get /MIJEONG/AIRS/V5/L2/${year}${month}/${year}${month}${day}.tar )
  fi

  # assume the original collection of data (hdf files, one per swath)
  # are in ../tars and that the filenames inside the tar files are named
  # YYYYMM/YYYYMMDD/*.hdf
  (cd ../data; tar -xvf ../tars/${year}${month}${day}.tar >> tarlog)
  
  # construct the input list of files for the converter.
  # cd there first in a subshell so the ls just contains simple file names
  (cd ../data/${year}${month}/${year}${month}${day}; ls AIR*hdf > flist)

  # get back to work dir and edit a template file to set the 
  # values that change in the namelists.
  sed -e "s/YYYY/${year}/g" \
      -e "s/MM/${month}/g"  \
      -e "s/DD/${day}/g"    \
      -e "s/GREG/${greg}/g" < ./input.nml.template > input.nml

  # actually make the obs_seq files, one per input.  these still need to
  # be merged if you want daily files.
  ./convert_airs_L2
  
  # do the merge now
  ls ${output_dir}/AIRS.${year}.${month}.${day}.*.out > olist
  ./obs_sequence_tool

  # start local mods
  # ok, this is a local mod - to try to keep from running out of disk space
  remote_dir=/gpfs/ptmp/dart/Obs_sets/AIRS_24_subx4_ascii/${year}${month}/
  cp -f ${output_dir}/AIRS.${year}${month}${day}.out $remote_dir
  # and clean up so we don't run out of disk space
  (cd ../data/${year}${month}/${year}${month}${day}; rm AIR*hdf)
  (cd ${output_dir}; rm AIRS.${year}.${month}.${day}.*.out)
  (cd ../tars; rm ${year}${month}${day}.tar)
  # end local mods

  # advance the day; the output is YYYYMMDD00
  curday=(`echo ${year}${month}${day}00 +1d | ./advance_time`)

  # advance the loop counter
  let d=d+1
 
done

exit 0

