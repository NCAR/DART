#!/bin/csh 
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id: gpsro_to_obsseq.csh 12574 2018-05-03 17:15:02Z nancy@ucar.edu $
#
########################################################################
#
#   convert_airs_L2.csh - script that downloads AIRS L2 retrievals
#               of temp/moisture, converts them into a series of
#               DART observations, and outputs a daily DART obs_seq file.
#
# can be used standalone but intended to be called by one of the two
# convert_xxx_airs_L2.csh scripts, or multi_parallel.batch (not .ksh), in this same directory.
#
# requires 5 args:
#    $1 - analysis date (yyyymmdd format)
#    $2 - base directory where files will be processed
#    $3 - yes to download AIRS L2 data files.  no if already downloaded.
#    $4 - yes to convert AIRS L2 data to an obs_seq file.  no to skip convert.
#    $5 - yes to delete AIRS L2 data when finished.  no to keep original files.
#
# update DART_DIR in this script to match your system, 
# and if downloading, set your web page user name and password.
# (No longer needed in 2021).
#
# edit the input.nml in the work directory (copied from ./input.nml(.template))
# [KDR; the template file has output_file = OUTDIR/YYYYMM/obs_seq.AIRS.YYYYMMDD.out,
#       which is expecting the scripts to replace the CAPs, which they don't do.]
# to select any options;
# it will be copied to the various places it is needed.
#
# the processing directory name is relative to the 'work' directory.
#
#     created June 2008, Ryan Torn NCAR/MMM
#     updated Nov  2008, nancy collins ncar/cisl
#     updated Aug  2009, nancy collins ncar/cisl
#     updated Oct  2010, Ryan and nancy
#     updated Jul  2011, nancy 
#     stolen from GPS converter, Aug 2018, nancy
#
########################################################################

########################################################################

# should only have to set the DART_DIR, and if downloading, set the
# web site user name and password for access.  expects to use the 'wget'
# utility to download files from the web page. 'curl' is another web 
# download program you can try although this script doesn't have any
# examples because it's not worked for me.

# top level directory (root where observations/AIRS dir is found)
#setenv DART_DIR    /home/user/DART
setenv DART_DIR    /glade/work/nancy/subversion/cleanrmatrunk.cheyenne

# NASA download site for AIRS data, individual files:
set AIRS_repository_path = 'https://airsl2.gesdisc.eosdis.nasa.gov/data/Aqua_AIRS_Level2/AIRX2RET.006'

# example single file name.  there are 240 files/day.
# the 3 digit number after the DD (day) and the ".L2.xxx" goes 
# from 001 to 240.
#
# https://airsl2.gesdisc.eosdis.nasa.gov/data/Aqua_AIRS_Level2/AIRX2RET.006/2010/349/AIRS.2010.12.15.062.L2.RetStd.v6.0.7.0.G13072192312.hdf

setenv DART_WORK_DIR  $DART_DIR/observations/obs_converters/AIRS/work
setenv CONV_PROG      convert_airs_L2
setenv DATE_PROG      advance_time
#setenv DOWNLOAD_DIR   rawdata
setenv DOWNLOAD_DIR   /glade/work/nancy/AIRS_data

# their example wget command:
# 2021; check the cookies requirement
#       Make the date_time part of the file name an variable, set with AIRS_repository_path?
#       Context; the get_files and downloads* filename in the README refer to 
#       $obs/AIRS/Daily_raw_files/get_files:       
#         wget
#         x  --background
#            --load-cookies /glade/u/home/nancy/.urs_cookies
#            --save-cookies /glade/u/home/nancy/.urs_cookies
#         x  --no-verbose
#            --auth-no-challenge=on
#            --keep-session-cookies
#            --content-disposition
#            -i downloads*.txt

touch ~/.urs_cookies
set wget_cmd = "wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies --content-disposition -i ./subset_AIRX2RET_V006_20180813_220243.txt"
# that last filename is a text file with a list of files
# to be downloaded.  FIXME here.

# set to no if you don't want as much output.
set chatty=yes

if ($# != 5) then
   echo usage: $0 date workdir downld convert cleanup 
   exit -1
endif

set datea   = $1     # target date, YYYYMMDD
set datadir = $2     # under what directory to process the files
set downld  = $3     # download?  'yes' or 'no'
set convert = $4     # convert?   'yes' or 'no'
set cleanup = $5     # delete downloaded files at end? 'yes' or 'no'


if ( $chatty == 'yes' ) then
   echo 'starting AIRS L2 script run at ' `date`
endif

# record where we started running
set origindir = `pwd`

# i've done this wrong enough times and wasted a lot of download 
# time, so do a bunch of bullet-proofing here before going on.
# e.g. verify the dirs all exist, the input.nml is in place.
# make sure the user has not edited the download/processing
# version of the input.nml - we use the one from the work dir.

if ( ! -d $DART_WORK_DIR ) then
  echo 'work directory not found: ' $DART_WORK_DIR
  exit 1
endif

echo 'current dir is ' `pwd`
if ( `pwd` != $DART_WORK_DIR ) then
   echo 'if not already there, changing directory to work dir.'
   cd $DART_WORK_DIR
   echo 'current dir now ' `pwd`
endif

if ( ! -d $datadir ) then
  echo 'creating data processing directory: ' $datadir
  mkdir -p $datadir
endif
if ( ! -e $datadir/input.nml ) then
  echo 'data processing directory does not contain an input.nml'
  echo 'copying from work dir to data proc dir'
  cp -fv ./input.nml $datadir
else
  diff -q ./input.nml $datadir/input.nml
  if ( $status == 1 ) then
     echo 'the input.nml file in the work directory is different'
     echo 'than the one in the data processing directory.'
     echo 'update them to be consistent, or remove the one in the'
     echo 'data processing directory and a new one will be copied'
     echo 'over from the work directory.'
     exit -1
  endif
endif

# copy over the date program and the converter from
# the work dir to the data processing directory.
# also the merge program and template for that.
# Make this conditional on $convert ?
cp -fv ./$DATE_PROG ./$CONV_PROG $datadir

########################################################################
# DOWNLOAD SECTION
########################################################################

if ( $downld == 'yes' ) then
   cd $datadir
   echo 'current dir now ' `pwd`

   if ( $chatty == 'yes' ) then
      echo 'starting raw file download at' `date`
   endif
   
   if ( ! -d $DOWNLOAD_DIR ) then
     echo 'creating new download directory: ' $DOWNLOAD_DIR
     mkdir -p $DOWNLOAD_DIR
   else
     echo 'existing directory found, cleaning up before new download'
     rm -fr $DOWNLOAD_DIR
     mkdir -p $DOWNLOAD_DIR
   endif
   
   cd $DOWNLOAD_DIR
   ln -sf $DART_WORK_DIR/input.nml .
   
   echo 'downloading tar file of obs for date: ' $datea 
   set yyyy    = `echo $datea | cut -b1-4`
   set mm      = `echo $datea | cut -b5-6`
   set dd      = `echo $datea | cut -b7-8`
   set dstring = ${yyyy}.${mm}.${dd}

   # get and untar, leaving all files in the current dir
   set day_file = daily_${datea}
# FIXME: put in the right wget command here
# each day has 240 files
   $wget_cmd -O $day_file $gps_repository_path/atmPrf/${dstring}
   if ( -z $day_file ) then
      echo no data available/downloaded for $datea
      rm $day_file
   else
      tar -xf $day_file --strip-components=3    # remove long path
      endif
   end
   rm input.nml dart_log.*
   
   set nfiles = `/bin/ls . | grep .hdf | wc -l`
   if ( $chatty == 'yes' ) then
      echo $nfiles total raw files downloaded at `date`
   endif
   if ( "$nfiles" == 0 ) then
      echo exiting: no profiles found for date $datea
      exit 1
   endif
   
   
   cd $DART_WORK_DIR
else
   echo 'not downloading data; assume it is already on local disk'
endif

######################################################################## 
# CONVERSION SECTION
########################################################################

if ( $convert == 'yes') then
   cd $datadir
   echo 'current dir now ' `pwd`
   
   if ( $chatty == 'yes' ) then
      echo 'starting AIRS conversion at ' `date`
   endif
   
   set yyyy    = `echo $datea | cut -b1-4`
   set mm      = `echo $datea | cut -b5-6`
   set dd      = `echo $datea | cut -b7-8`
   echo 'converting obs for date: ' $datea
   set dstring = ${yyyy}.${mm}.${dd}
   
   rm -f flist
   /bin/ls -1 ${DOWNLOAD_DIR}/${dstring}/AIRS.${dstring}.*.hdf >! flist
   
   set nfiles = `cat flist | wc -l`
   if ( $chatty == 'yes' ) then
      echo $nfiles ' to process for file ' $datea 
   endif
   
   ./$CONV_PROG >>! convert_output_log

   # if it worked, add a date string
   if ( -e obs_seq.out ) then
      mv obs_seq.out obs_seq.${datea}
      /bin/ls -l obs_seq.${datea}
   else
      if ( $chatty == 'yes' ) then
         echo "no obs found for date $datea, or conversion failed."
      endif
   endif

   cd $DART_WORK_DIR
else
   echo 'not converting data'
endif

########################################################################
# FILE REMOVAL SECTION
########################################################################

if ( $cleanup == 'yes' ) then
   cd $datadir
   echo 'current dir now ' `pwd`
   
   if ( $chatty == 'yes' ) then
      echo 'cleaning up files starting at ' `date`
   endif

   echo 'removing original AIRS L2 data files for date: ' $datea

   # just remove the whole subdir here.  trying to list individual
   # files can cause problems with long command line lengths.
   rm -fr $DOWNLOAD_DIR

   cd $DART_WORK_DIR
else
   echo 'not removing original AIRS L2 data files'
endif

########################################################################

if ( $chatty == 'yes' ) then
   echo 'finished AIRS L2 script run at ' `date`
   echo ' '
endif

exit 0

# <next few lines under version control, do not edit>
# $URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/rma_trunk/observations/obs_converters/gps/shell_scripts/gpsro_to_obsseq.csh $
# $Revision: 12574 $
# $Date: 2018-05-03 11:15:02 -0600 (Thu, 03 May 2018) $

