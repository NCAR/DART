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

if ($# != 2) then
    echo "usage: $0 date rawdir outputdir"
    echo "date: YYYYMMDD"
    echo "rawdir: directory of timeslice input data"
    exit -1
endif

set     datea = $1     # target date, YYYYMMDD
set    rawdir = $2     # directory of raw input data
set outputdir = $3     # output directory

if ( $chatty == 'yes' ) then
   echo 'Starting USGS obs converter script: ' `date`
endif

# record where we started running
set origindir = `pwd`
setenv CONV_PROG    CONVPROG_TEMPLATE

# Every day creates files in its own dir.
set datedir = day_proc.${datea}
mkdir -p $datedir
cd $datedir
echo 'current dir: '`pwd`

ln -s $origindir/$CONV_PROG
ln -s $origindir/input.nml
ln -s $origindir/hydro.namelist
TOP_LEVEL_CONFIG_TEMPLATE
HYDRO_RST_TEMPLATE
ln -s $origindir/ROUTELINK_TEMPLATE
ln -s $origindir/GAGES_LIST_FILE_TEMPLATE

set yyyy    = `echo $datea | cut -b1-4`
set mm      = `echo $datea | cut -b5-6`
set dd      = `echo $datea | cut -b7-8`
echo "Converting obs for date:  $datea = $yyyy $mm $dd"

ls -1 $rawdir/${yyyy}-${mm}-${dd}* >! list_of_obs_files.txt
./$CONV_PROG || exit 2
mv -v obs_seq.out obs_seq.$datea

cd $origindir

ln -s $datedir/obs_seq.$datea .

exit 0


