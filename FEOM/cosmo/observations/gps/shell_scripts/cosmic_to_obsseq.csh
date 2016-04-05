#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
########################################################################
#
#   cosmic_to_obsseq.csh - script that downloads COSMIC observations 
#               and converts them to a DART observation sequence file.
#
# requires 5 args:
#    $1 - analysis date (yyyymmdd format)
#    $2 - base directory where files will be processed
#    $3 - yes to download raw COSMIC data 
#    $4 - yes to convert COSMIC data to an obs_seq file
#    $5 - yes to delete raw COSMIC data when finished
#
# update DART_DIR in this script to match your system, and if
# downloading, set your COSMIC web page user name and password.
#
# edit the input.nml in the work directory to select any options;
# it will be copied to the various places it is needed.
#
# the processing directory name is relative to the 'work' directory.
#
#
#     created June 2008, Ryan Torn NCAR/MMM
#     updated Nov  2008, nancy collins ncar/cisl
#     updated Aug  2009, nancy collins ncar/cisl
#     updated Oct  2010, Ryan and nancy
# 
#
# ------- 
# from the cosmic web site about the use of 'wget' to download
# the many files needed to do this process:
#
#   Hints for using wget for fetching CDAAC files from CDAAC:
#   
#   Here is one recipe for fetching all cosmic real time atmPrf files for one day:
#   
#   wget -nd -np -r -l 10 -w 2 --http-user=xxxx --http-passwd=xxxx \
#         http://cosmic-io.cosmic.ucar.edu/cdaac/login/cosmicrt/level2/atmPrf/2009.007/
#   
#   The option -np (no parents) is important. Without it, all manner of 
#   files from throughout the site will be loaded, I think due to the 
#   links back to the main page which are everywhere.
#   
#   The option -r (recursive fetch) is necessary unless there is just 
#   one file you want to fetch.
#   
#   The option -l 10 (limit of recursive depth to 10 levels) is necessary 
#   in order to get around the default 5 level depth.
#   
#   The option -nd dumps all fetched files into your current directory. 
#   Otherwise a directory hierarchy will be created: 
#     cosmic-io.cosmic.ucar.edu/cdaac/login/cosmic/level2/atmPrf/2006.207
#   
#   The option -w 2 tells wget to wait two seconds between each file fetch 
#   so as to have pity on the poor web server.
# ------- 
# 
# note: there are between 1000 and 3000 files per day.  without the -w
# flag, i was getting about 5 files per second (each individual file is
# relatively small).  but with -w 1 obviously we get slightly less than
# a file a second, -w 2 is half that again.  this script uses -w 1 by
# default, but if you are trying to download a lot of days i'd recommend
# removing it.
#
########################################################################

########################################################################

# should only have to set the DART_DIR, and if downloading, set the
# web site user name and password for access.  expects to use the 'wget'
# utility to download files from the web page.  

# top level directory (where observations dir is found)
setenv DART_DIR    /home/user/dart
set cosmic_user    = xxx
set cosmic_pw      = yyy

set gps_repository_path = 'http://cosmic-io.cosmic.ucar.edu/cdaac/login'

setenv DART_WORK_DIR  ${DART_DIR}/observations/gps/work
setenv CONV_PROG      convert_cosmic_gps_cdf
setenv DATE_PROG      advance_time

# if you are in a hurry, this has no delay between requests:
#set wget_cmd           = 'wget -q -nd -np -r -l 10 -w 0'
set wget_cmd            = 'wget -q -nd -np -r -l 10 -w 1'

set chatty=yes

set datea   = ${1}     # target date, YYYYMMDD
set datadir = ${2}     # under what directory to process the files
set downld  = ${3}     # download?  'yes' or 'no'
set convert = ${4}     # convert?   'yes' or 'no'
set cleanup = ${5}     # delete COSMIC files at end? 'yes' or 'no'


if ( $chatty == 'yes' ) then
   echo 'starting gps script run at ' `date`
endif

# i've done this wrong enough times and wasted a lot of download 
# time, so do a bunch of bullet-proofing here before going on.
# e.g. verify the dirs all exist, the input.nml is in place.
# make sure the user has not edited the download/processing
# version of the input.nml - we use the one from the work dir.

if ( ! -d ${DART_WORK_DIR} ) then
  echo 'work directory not found: ' ${DART_WORK_DIR}
  exit 1
endif

echo 'current dir is ' `pwd`
if ( `pwd` != ${DART_WORK_DIR} ) then
   echo 'if not already there, changing directory to work dir.'
   cd ${DART_WORK_DIR}
   echo 'current dir now ' `pwd`
endif

if ( ! -d ${datadir} ) then
  echo 'creating data processing directory: ' ${datadir}
  mkdir ${datadir}
endif
if ( ! -e ${datadir}/input.nml ) then
  echo 'data processing directory does not contain an input.nml'
  echo 'copying from work dir to data proc dir'
  echo `pwd`/input.nml '->' ${datadir}/input.nml
  cp -f ./input.nml ${datadir}
else
  diff -q ./input.nml ${datadir}/input.nml
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
cp -f ./${DATE_PROG} ./${CONV_PROG} ${datadir}

if ( $downld == 'yes' ) then
   echo 'if not already there, changing dir to data proc directory'
   cd ${datadir}
   echo 'current dir now ' `pwd`

   if ( ! $?cosmic_user || ! $?cosmic_pw ) then
      echo "You must set cosmic_user to your username for the cosmic web site"
      echo "and set cosmic_pw to your password, then rerun this script. "
      exit -1
   endif

   if ( $chatty == 'yes' ) then
      echo 'starting raw file download at' `date`
   endif
   
   set get = "${wget_cmd} --http-user=${cosmic_user} --http-passwd=${cosmic_pw}" 
   set yyyy   = `echo $datea | cut -b1-4`
   
   if ( ! -d ${datea} ) then
     echo 'year/month/day directory not found: ' ${datea}
     echo 'creating now.'
     mkdir ${datea}
   else
     echo 'existing directory found, cleaning up before new download'
     rm -fr ${datea}
     mkdir ${datea}
   endif
   
   cd ${datea}
   ln -sf ${DART_WORK_DIR}/input.nml .
   
   set jyyyydd = `echo $datea 0 -j | ../${DATE_PROG}` 
   @ mday = $jyyyydd[2] + 1000  ;  set mday = `echo $mday | cut -b2-4`
   echo 'downloading obs for date: ' $datea ', which is julian day: ' $jyyyydd

   ${get} ${gps_repository_path}/cosmic/level2/atmPrf/${yyyy}.${mday}/
   rm -f *.html *.txt
   ${get} ${gps_repository_path}/champ/level2/atmPrf/${yyyy}.${mday}/
   rm -f *.html *.txt input.nml
   
   if ( $chatty == 'yes' ) then
      echo `/bin/ls . | grep _nc | wc -l` 'raw files downloaded at ' `date`
   endif
   
   cd ${DART_WORK_DIR}
else
   echo 'not downloading data; assume it is already on local disk'
endif

if ( $convert == 'yes') then
   echo 'if not already there, changing dir to data proc directory'
   cd ${datadir}
   echo 'current dir now ' `pwd`
   
   if ( $chatty == 'yes' ) then
      echo 'starting gps conversion at ' `date`
   endif
   
   rm -f flist
   set yyyy    = `echo $datea | cut -b1-4`
   set jyyyydd = `echo ${datea}00 0 -j | ./${DATE_PROG}`
   @ mday = $jyyyydd[2] + 1000  ;  set mday = `echo $mday | cut -b2-4`
   echo 'converting obs for date: ' $datea
   
   /bin/ls -1 ${datea}/*.${yyyy}.${mday}.*_nc >! flist
   
   set nfiles = `cat flist | wc -l`
   if ( $chatty == 'yes' ) then
      echo $nfiles ' to process for file ' $datea 
   endif
   
   ./${CONV_PROG} >>! convert_output_log

   rm -rf cosmic_gps_input.nc flist
   if ( -e obs_seq.gpsro ) then
       mv obs_seq.gpsro obs_seq.gpsro_${datea}
   
      if ( $chatty == 'yes' ) then
         echo "all observations for day in file obs_seq.gpsro_${datea} at " `date`
      endif
   else
      if ( $chatty == 'yes' ) then
         echo "no obs found for date ${datea}, or conversion failed.'
      endif
   endif

   cd ${DART_WORK_DIR}
else
   echo 'not converting data'
endif

if ( $cleanup == 'yes' ) then
   echo 'if not already there, changing dir to data proc directory'
   cd ${datadir}
   echo 'current dir now ' `pwd`
   
   if ( $chatty == 'yes' ) then
      echo 'cleaning up files at ' `date`
   endif

   echo 'removing original cosmic data files for date: ' $datea

   # just remove the whole subdir here.  trying to list individual
   # files can cause problems with long command line lengths.
   rm -fr ${datea}

   cd ${DART_WORK_DIR}
else
   echo 'not removing original cosmic data files'
endif

if ( $chatty == 'yes' ) then
   echo 'finished gps script run at ' `date`
   echo ' '
endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

