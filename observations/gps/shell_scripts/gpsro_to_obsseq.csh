#!/bin/csh 
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
########################################################################
#
#   gpsro_to_obsseq.csh - script that downloads GPS radio occultation
#               observation profiles and converts them into a series of
#               DART observations and outputs a daily DART obs_seq file.
#
# can be used standalone but intended to be called by the do_convert.csh
# script in this same directory.
#
# requires 6 args:
#    $1 - analysis date (yyyymmdd format)
#    $2 - base directory where files will be processed
#    $3 - yes to download RO data profiles.  no if already downloaded.
#    $4 - yes to convert RO data to an obs_seq file.  no to skip convert.
#    $5 - yes to delete RO data when finished.  no to keep original files.
#    $6 - file containing the names of the satellites to get data from
#
# update DART_DIR in this script to match your system, and if
# downloading, set your CDAAC web page user name and password.
#
# edit the input.nml in the work directory to select any options;
# it will be copied to the various places it is needed.
#
# if the processing directory name is a relative path, it should
# be relative to the gps converter 'work' directory.
#
# CHECK these against the cdaac web site for latest info.  as an
# example, a snapshot of available data as of May 2016 is:
#
# champ2014  2001.138 - 2008.279 
# cnofs      2010.060 - 2011.365
# cnofsrt    2012.001 - 2015.193
# cosmic2013 2006.112 - 2014.120
# cosmic     2014.121 - 2015.364
# cosmicrt   2014.181 - 2016.123
# gpsmet     1995.111 - 1997.047
# gpsmetas   1995.237 - 1997.016
# grace      2007.059 - 2015.364
# kompsat5rt 2015.305 - 2016.123
# metopa2016 2007.274 - 2015.365
# metopb     2013.032 - 2015.059
# sacc       2006.068 - 2011.215
# saccrt     2011.329 - 2013.226
# tsx        2008.041 - 2015.333
#
#     created June 2008, Ryan Torn NCAR/MMM
#     updated Nov  2008, nancy collins ncar/cisl
#     updated Aug  2009, nancy collins ncar/cisl
#     updated Oct  2010, Ryan and nancy
#     updated Jul  2011, nancy (added support for other satellites)
#     revised May  2016, nancy (full support for daily tar files)
# 
#
# This script uses the 'wget' command to download data from the CDAAC
# web site.  Data is available as a single tar file per satellite per day.  
# The profiles are stored one per file after untarring the archive.  
# The current filename scheme includes the satellite and YYYY.DOY, e.g.:
#
# wget http://cdaac-www.cosmic.ucar.edu/cdaac/rest/tarservice/data/cosmic/atmPrf/2012.304
#        -O cosmic_atmPrf_2012.304.tar 
#   
# If this changes (again) this script will need to be updated.
#
########################################################################

# start of user-settable section - should only need to set these once

# top level DART directory
setenv DART_DIR    /glade/p/home/$USER/DART

# your CDAAC web site user name and password
set cdaac_user    = nscollins
set cdaac_pw      = xxxxxxx

# end of user-settable section

########################################################################

setenv DART_WORK_DIR  ${DART_DIR}/observations/gps/work
setenv CONV_PROG      convert_cosmic_gps_cdf
setenv DATE_PROG      advance_time

# CDAAC web site path:
set gps_repository_path = 'http://cdaac-www.cosmic.ucar.edu/cdaac/rest/tarservice/data'

# wget seems to print a lot of output showing progress.
# this flag tries to print less often but i haven't found a way to turn 
# it off completely.  --progress:none was discussed in a wget forum 
# but i see no indication it was implemented
set wget_cmd            = 'wget --progress=dot:mega '

# this helps with debugging and isn't really that verbose
set chatty=yes

# start of executable stuff

if ($# != 6) then
   echo usage: $0 date workdir downld convert cleanup satlist
   exit -1
endif

set datea   = ${1}     # target date, YYYYMMDD
set datadir = ${2}     # under what directory to process the files
set downld  = ${3}     # download?  'yes' or 'no'
set convert = ${4}     # convert?   'yes' or 'no'
set cleanup = ${5}     # delete downloaded files at end? 'yes' or 'no'
set satlist = ${6}     # file with list of satellites to use


if ( $chatty == 'yes' ) then
   echo 'starting gpsro script run at ' `date`
endif

# record where we started running
set origindir = `pwd`

# i've done this wrong enough times and wasted a lot of download 
# time, so do a bunch of bullet-proofing here before going on.
# e.g. verify the dirs all exist, the input.nml is in place.
# make sure the user has not edited the download/processing
# version of the input.nml - we use the one from the work dir.

if ( ! -d ${DART_WORK_DIR} ) then
  echo 'work directory not found: ' ${DART_WORK_DIR}
  exit 1
endif

# copy satlist to workdir
if ( ! -e ${DART_WORK_DIR}/${satlist} ) then
  cp -f $satlist $DART_WORK_DIR
else
  diff -q $satlist ${DART_WORK_DIR}/${satlist}
  if ( $status == 1 ) then
     echo "the satellite list file ${satlist} in the work directory is different"
     echo "than the one in the $origindir directory.  "
     echo "update them to be consistent, or remove the one in the"
     echo "work directory and a new one will be copied"
     echo "over from the $origindir directory."
     exit -1
  endif
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

# copy satlist to datadir
if ( ! -e ${datadir}/${satlist} ) then
  echo 'data processing directory does not contain a satellite list file'
  echo 'copying from work dir to data proc dir'
  echo `pwd`/${satlist}'->' ${datadir}/${satlist}
  cp -f ${satlist} ${datadir}/
else
  diff -q ${satlist} ${datadir}/${satlist}
  if ( $status == 1 ) then
     echo "the satellite list file ${satlist} in the work directory is different"
     echo "than the one in the ${datadir} directory.  "
     echo "update them to be consistent, or remove the one in the"
     echo "${datadir} directory and a new one will be copied"
     echo "over from the ${DART_WORK_DIR} directory."
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

   if ( ! $?cdaac_user || ! $?cdaac_pw ) then
      echo "You must set cdaac_user to your username for the cdaac web site"
      echo "and set cdaac_pw to your password, then rerun this script. "
      exit -1
   endif

   if ( $chatty == 'yes' ) then
      echo 'starting raw file download at' `date`
   endif
   
   set get = "${wget_cmd} --http-user=${cdaac_user} --http-passwd=${cdaac_pw}" 
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

   ## here is where you select which satellites to download
   foreach sat ( `cat ../$satlist` )
      if ( $chatty == 'yes' ) then
         echo 'copying data files from satellite: ' ${sat}
      endif
      # new, tar file per day!
      echo ${get} ${gps_repository_path}/${sat}/atmPrf/${yyyy}.${mday} -O ${sat}_atmPrf_${yyyy}.${mday}.tar
      ${get} ${gps_repository_path}/${sat}/atmPrf/${yyyy}.${mday} -O ${sat}_atmPrf_${yyyy}.${mday}.tar
   end
   rm input.nml
   
   if ( $chatty == 'yes' ) then
      echo `/bin/ls *_atmPrf_${yyyy}.${mday}.tar | wc -l` 'tar files downloaded at ' `date`
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
      echo 'starting gpsro conversion at ' `date`
   endif
     
   set yyyy    = `echo $datea | cut -b1-4`
   set yyyymm  = `echo $datea | cut -b1-6`
   set jyyyydd = `echo ${datea}00 0 -j | ./${DATE_PROG}`
   @ mday = $jyyyydd[2] + 1000  ;  set mday = `echo $mday | cut -b2-4`
   echo 'converting obs for date: ' $datea

   cd ${datea}
   cp ../input.nml .
   echo 'current dir now ' `pwd`

   # make sure directory has no leftovers from before
   rm -fr atmPrf*_nc flist obs_seq.gpsro

   foreach sat ( `cat ../$satlist` )
     echo 'converting obs for satellite: ' $sat

     set next_tarfile = ${sat}_atmPrf_${yyyy}.${mday}.tar 

     if ( ! -e $next_tarfile || -z $next_tarfile ) then
       echo $next_tarfile NOT DOWNLOADED or ZERO LENGTH, SKIPPING
       continue
     else
       echo untarring $next_tarfile into daily profiles
     endif

     tar --strip-components=3 -xf $next_tarfile
     /bin/ls -1 atmPrf_*.${yyyy}.${mday}.*_nc >! flist
     
     set nfiles = `cat flist | wc -l`
     if ( $chatty == 'yes' ) then
      echo $nfiles $sat ' profiles to process for day ' $datea 
     endif
   
     ../${CONV_PROG} >>! ../convert_output_log

     # keep tar files but remove individual profile files so next
     # pass doesn't add them into the output file.
     rm -rf flist atmPrf_*.${yyyy}.${mday}.*_nc 
   end

   if ( -e obs_seq.gpsro ) then
     mv obs_seq.gpsro ../obs_seq.gpsro_${datea}
     
     if ( $chatty == 'yes' ) then
       echo "all observations for day in file obs_seq.gpsro_${datea} at " `date`
     endif
   else
     if ( $chatty == 'yes' ) then
       echo "no obs found for date ${datea}, or conversion failed."
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

   echo 'removing original gpsro data files for date: ' $datea

   # just remove the whole subdir here.  trying to list individual
   # files can cause problems with long command line lengths.
   rm -fr ${datea}

   cd ${DART_WORK_DIR}
else
   echo 'not removing original gpsro data files'
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

