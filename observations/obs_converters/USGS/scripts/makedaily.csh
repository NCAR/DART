#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# split file(s) into "daily" files which start at 23:30:01 Z
# the previous day and end at 23:30:00 of the supplied date.

# requires 2 args:
#    $1 - analysis date (yyyymmdd format)
#    $2 - directory containing the raw input TimeSlice.ncdf files
#    $3 - base directory for output files 

# set to no if you don't want as much output.
set chatty=yes

if ($# != 3) then
   echo usage: $0 date rawdir outputdir
   exit -1
endif

set     datea = $1     # target date, YYYYMMDD
set    rawdir = $2     # directory of raw input data
set outputdir = $3     # output directory

if ( $chatty == 'yes' ) then
   echo 'starting gpsro script run at ' `date`
endif

# record where we started running
set origindir = `pwd`

# do a bunch of bullet-proofing here before going on.
# e.g. verify the dirs all exist, the input.nml is in place.

setenv DART_DIR       /glade/p/work/thoar/DART/wrf_hydro_dart_git/wrf_hydro_dart
setenv DART_WORK_DIR  $DART_DIR/observations/obs_converters/USGS/work
setenv CONV_PROG      convert_streamflow

if ( ! -d $DART_WORK_DIR ) then
  echo 'DART work directory not found: ' $DART_WORK_DIR
  exit 1
endif

if ( ! -d $outputdir ) then
  echo 'creating data processing directory: ' $outputdir
  mkdir -p $outputdir
endif

cd $outputdir
echo 'current dir is now '`pwd`

# copy DART parts from to the data processing directory.
foreach FILE ( input.nml RouteLink.nc wanted_gages_list.txt $CONV_PROG )
   if ( ! -e $FILE ) then
     echo "data processing directory does not contain $FILE"
     echo 'copying from DART dir to data proc dir'
     echo  $DART_WORK_DIR/$FILE '->' $outputdir/$FILE
     cp -f $DART_WORK_DIR/$FILE .
   endif
end

set yyyy    = `echo $datea | cut -b1-4`
set mm      = `echo $datea | cut -b5-6`
set dd      = `echo $datea | cut -b7-8`
echo 'converting obs for date: ' $datea  $yyyy $mm $dd

ls -1 $rawdir/${yyyy}-${mm}-${dd}* >! list_of_obs_files.txt

./$CONV_PROG || exit 2

mv -v obs_seq.out obs_seq.$datea

exit 0


