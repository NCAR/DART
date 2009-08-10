#!/bin/csh 
########################################################################
#
#   cosmic_download.csh - script that downloads COSMIC observations.
#         then you can run cosmic_to_obsseq.csh with a third argument of
#         'no' so it does not re-download the same files.
#
# requires 2 args:
#    $1 - analysis date (yyyymmdd format)
#    $2 - base observation directory
#
# and update the 3 env settings at the top to match your system.
#
#     created May  2009, nancy collins ncar/cisl
#                        converted from cosmic_to_obsseq.csh
#     updated Aug  2009, nancy collins ncar/cisl
#
# from the cosmic web site about the use of 'wget' to download
# the many files needed to do this process:
# ------- 
# Hints for using wget for fetching CDAAC files from CDAAC:
# 
# Here is one recipe for fetching all cosmic real time atmPrf files for one day:
# 
# wget -nd -np -r -l 10 -w 2 --http-user=xxxx --http-passwd=xxxx http://cosmic-io.cosmic.ucar.edu/cdaac/login/cosmicrt/level2/atmPrf/2009.007/
# 
# The option -np (no parents) is important. Without it, all manner of 
# files from throughout the site will be loaded, I think due to the 
# links back to the main page which are everywhere.
# 
# The option -r (recursive fetch) is necessary unless there is just 
# one file you want to fetch.
# 
# The option -l 10 (limit of recursive depth to 10 levels) is necessary 
# in order to get around the default 5 level depth.
# 
# The option -nd dumps all fetched files into your current directory. 
# Otherwise a directory hierarchy will be created: 
#   cosmic-io.cosmic.ucar.edu/cdaac/login/cosmic/level2/atmPrf/2006.207
# 
# The option -w 2 tells wget to wait two seconds between each file fetch 
# so as to have pity on our poor web server.
# ------- 
# 
########################################################################

# should only have to set the DART_DIR and the rest should be in the
# right place relative to it.
setenv DART_DIR      /home/user/dart
setenv cosmic_user   xxx
setenv cosmic_pw     yyy

setenv DART_WORK_DIR  ${DART_DIR}/observations/gps/work
setenv CONV_PROG      convert_cosmic_gps_cdf
setenv DATE_PROG      advance_time

set chatty=yes
set downld=yes

set datea   = ${1}     # target date, YYYYMMDD
set datadir = ${2}     # where to process the files

set assim_freq          = 6  # hours, sets centers of windows.
set download_window     = 3  # window half-width (some users choose 2 hours)
set gps_repository_path = 'http://cosmic-io.cosmic.ucar.edu/cdaac/login'
set wget_cmd            = 'wget -q -nd -np -r -l 10 -w 1'

# i've done this wrong enough times and wasted a lot of download time,
# so do a bunch of bullet-proofing here before going on.

# verify the dirs all exist, the input.nml is in place.
if ( ! -d ${DART_WORK_DIR} ) then
  echo 'work directory not found: ' ${DART_WORK_DIR}
  exit
endif

echo 'current dir is ' `pwd`
if ( `pwd` != ${DART_WORK_DIR} ) then
   echo 'if not already there, changing directory to work dir.'
   cd ${DART_WORK_DIR}
   echo 'current dir now ' `pwd`
endif

if ( ! -d ${datadir} ) then
  echo 'data processing directory not found: ' ${datadir}
  echo 'creating now.'
  mkdir ${datadir}
  ls -ld ${datadir}
endif

if ( ! -e ${datadir}/${DATE_PROG} ) then
  echo 'data processing directory does not contain the time converter'
  echo 'copying from work dir to data proc dir'
  echo `pwd`/${DATE_PROG} '->' ${datadir}/${DATE_PROG}
  cp -f ./${DATE_PROG} ${datadir}
else
  echo 'using time conversion program already found in data proc dir'
endif

echo 'changing dir to data proc directory'
cd ${datadir}
echo 'current dir now ' `pwd`

 if ( ! $?cosmic_user || ! $?cosmic_pw ) then
    echo "You must setenv cosmic_user to your username for the cosmic web site"
    echo "and setenv cosmic_pw to your password. (or export cosmic_user=val for"
    echo "ksh/bash users) then rerun this script. "
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
endif

cd ${datea}

echo $datea
set jyyyydd = `echo $datea 0 -j | ../${DATE_PROG}` 
echo $jyyyydd
@ mday = $jyyyydd[2] + 1000  ;  set mday = `echo $mday | cut -b2-4`
${get} ${gps_repository_path}/cosmic/level2/atmPrf/${yyyy}.${mday}/
rm -f *.html *.txt
${get} ${gps_repository_path}/champ/level2/atmPrf/${yyyy}.${mday}/
rm -f *.html *.txt

set jyyyydd = `echo $datea 24 -j | ../${DATE_PROG}`
@ mday = $jyyyydd[2] + 1000  ;  set mday = `echo $mday | cut -b2-4` 
${get} ${gps_repository_path}/cosmic/level2/atmPrf/${yyyy}.${mday}/
rm -f *.html *.txt
${get} ${gps_repository_path}/champ/level2/atmPrf/${yyyy}.${mday}/
rm -f *.html *.txt

if ( $chatty == 'yes' ) then
   # the ls arg list line gets too long in some cases
   echo `/bin/ls | grep _nc | wc -l` 'raw files'
   echo 'all raw files download at ' `date`
endif

cd ..

exit 0


