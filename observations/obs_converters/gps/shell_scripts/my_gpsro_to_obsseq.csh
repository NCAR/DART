#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
########################################################################
#
#   gpsro_to_obsseq.csh - script that downloads GPS radio occultation
#               observation profiles and converts them into a series of
#               DART observations and outputs a daily DART obs_seq file.
#
# can be used standalone but intended to be called by one of the two
# convert_xxx_gpsro.csh scripts in this same directory.
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
# the processing directory name is relative to the 'work' directory.
#
# the options for satellite names, and available data times are:
#
#    cosmic:      all 6 COSMIC : 2006.194 - now*
#    cosmic-2013: reprocessed COSMIC : 2006.194 - now*
#    cosmicrt:    COSMIC : realtime
#    champ:       CHAMP : 2001.139 - 2008.274
#    grace:       Grace-A : 2007.059 - now*
#    tsx:         German TerraSAR-X : 2008.041 - now*
#    metopa:      Metop-A/GRAS : 2008.061 - now*
#    sacc:        Argentinan SAC-C : 2006.068 - now*
#    saccrt:      SAC-C : realtime
#    ncofs:       new Air Force C/NOFS : 2010.335 - now*
#    ncofsrt:     C/NOFS : realtime
#
#  - dates are YYYY.DDD where DDD is day number in year
#  - now* means current date minus 3-4 months.  reprocessed data
#    lags that much behind.  realtime means up to today's date,
#    with less quality control.
#  - only select one of reprocessed or realtime for a satellite 
#    or you will get duplicated observations.
#
#
#     created June 2008, Ryan Torn NCAR/MMM
#     updated Nov  2008, nancy collins ncar/cisl
#     updated Aug  2009, nancy collins ncar/cisl
#     updated Oct  2010, Ryan and nancy
#     updated Jul  2011, nancy (added support for other satellites)
# 
#
# ------- 
# From the CDAAC web site about the use of 'wget' to download
# the many files needed to do this process:
#
#   Hints for using wget for fetching CDAAC files from CDAAC:
#   
#   Here is one recipe for fetching one cosmic realtime atmPrf tar file for one day:
#   
#   wget -q --http-user=xxxx --http-passwd=xxxx \
#         http://cdaac-www.cosmic.ucar.edu/cdaac/rest/tarservice/data/cosmicrt/atmPrf/2018.017
#   
# at various times we've used some of these flags, although most of them
# are unneeded at this time.
#
#   The option -q (quiet) suppresses the progress output which can fill
#   your log files quickly if you're downloading lots of files.
#
#   The option -O file (direct output to given file) lets you rename the
#   file or put it in a different location easily if you are only downloading
#   a single file.
#
#   The option -np (no parents) tells wget not to follow back links to
#   higher level directories.
#   
#   The option -r (recursive fetch) is needed if there is more than 
#   one file you want to fetch.
#   
#   The option -l 10 (limit of recursive depth to 10 levels) allows you
#   to download a deep directory structure if you're also using -r.
#   
#   The option -nd dumps all fetched files into your current directory. 
#   Otherwise a directory hierarchy will be created: 
#     cosmic-io.cosmic.ucar.edu/cdaac/login/cosmic/level2/atmPrf/2006.207
#   
#   The option -w 2 tells wget to wait two seconds between each file fetch 
#   so as to have pity on the poor web server.
#
########################################################################

########################################################################

# should only have to set the DART_DIR, and if downloading, set the
# web site user name and password for access.  expects to use the 'wget'
# utility to download files from the web page. 'curl' is another web 
# download program you can try although this script doesn't have any
# examples because it's not worked for me.

# top level directory (root where observations/gps dir is found)
setenv DART_DIR   /glade/scratch/$USER
set cdaac_user    = username
set cdaac_pw      = password

# cdaac download site, daily tar files:
set gps_repository_path = 'http://cdaac-www.cosmic.ucar.edu/cdaac/rest/tarservice/data'

# example daily tar file name:
#http://cdaac-www.cosmic.ucar.edu/cdaac/rest/tarservice/data/tsx/atmPrf/2017.213

setenv DART_WORK_DIR  $DART_DIR/gps_conversion/work
setenv CONV_PROG      convert_cosmic_gps_cdf
setenv DATE_PROG      advance_time
setenv DOWNLOAD_DIR   rawdata

# options on wget command:
set wget_cmd = "wget -q --http-user=$cdaac_user --http-passwd=$cdaac_pw" 

# set to no if you don't want as much output.
set chatty=yes

if ($# != 6) then
   echo usage: $0 date workdir downld convert cleanup satlist
   exit -1
endif

set datea   = $1     # target date, YYYYMMDD
set datadir = $2     # under what directory to process the files
set downld  = $3     # download?  'yes' or 'no'
set convert = $4     # convert?   'yes' or 'no'
set cleanup = $5     # delete downloaded files at end? 'yes' or 'no'
set satlist = $6     # file with list of satellites to use


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

if ( ! -d $DART_WORK_DIR ) then
  echo 'work directory not found: ' $DART_WORK_DIR
  exit 1
endif

# copy satlist to workdir
# (override any existing satlist files at destination)
#if ( ! -e $DART_WORK_DIR/$satlist ) then
  cp -f $satlist $DART_WORK_DIR
#else
#  diff -q $satlist $DART_WORK_DIR/$satlist
#  if ( $status == 1 ) then
#     echo "the satellite list file $satlist in the work directory is different"
#     echo "than the one in the $origindir directory.  "
#     echo "update them to be consistent, or remove the one in the"
#     echo "work directory and a new one will be copied"
#     echo "over from the $origindir directory."
#     exit -1
#  endif
#endif

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
  echo `pwd`/input.nml '->' $datadir/input.nml
  cp -f ./input.nml $datadir
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

# copy satlist to datadir
# (override any existing satlist file)
#if ( ! -e $datadir/$satlist ) then
#  echo 'data processing directory does not contain a satellite list file'
  echo 'copying from work dir to data proc dir'
  echo `pwd`/$satlist'->' $datadir/$satlist
  cp -f $satlist $datadir/
#else
#  diff -q $satlist $datadir/$satlist
#  if ( $status == 1 ) then
#     echo "the satellite list file $satlist in the work directory is different"
#     echo "than the one in the $datadir directory.  "
#     echo "update them to be consistent, or remove the one in the"
#     echo "$datadir directory and a new one will be copied"
#     echo "over from the $DART_WORK_DIR directory."
#     exit -1
#  endif
#endif

# copy over the date program and the converter from
# the work dir to the data processing directory.
cp -f ./$DATE_PROG ./$CONV_PROG $datadir

########################################################################
# DOWNLOAD SECTION
########################################################################

if ( $downld == 'yes' ) then
   cd $datadir
   echo 'current dir now ' `pwd`

   if ( ! $?cdaac_user || ! $?cdaac_pw ) then
      echo "You must set cdaac_user to your username for the cdaac web site"
      echo "and set cdaac_pw to your password, then rerun this script. "
      exit -1
   endif

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
   
   set jyyyydd = `echo $datea 0 -j | ../$DATE_PROG` 
   set yyyy    = `echo $datea | cut -b1-4`
   set mday    = `printf %03d $jyyyydd[2]`
   echo 'downloading tar file of obs for date: ' $datea ', which is julian day: ' $jyyyydd

   # do this loop once for each satellite source you've selected
   foreach sat ( `cat ../$satlist` )
      if ( $chatty == 'yes' ) then
         echo 'copying data files from satellite: ' $sat
      endif
       
      # get and untar, leaving all files in the current dir
      set day_file = daily_${sat}_${yyyy}_${mday}
      $wget_cmd -O $day_file $gps_repository_path/$sat/atmPrf/${yyyy}.${mday}
      if ( -z $day_file ) then
         echo no data available/downloaded for satellite $sat on $yyyy $mday
         rm $day_file
      else
         tar -xf $day_file --strip-components=3    # remove long path
      endif
   end
   rm input.nml dart_log.*
   
   set nfiles = `/bin/ls . | grep _nc | wc -l`
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
      echo 'starting gpsro conversion at ' `date`
   endif
   
   set jyyyydd = `echo ${datea}00 0 -j | ./$DATE_PROG`
   set yyyy    = `echo $datea | cut -b1-4`
   set mday    = `printf %03d $jyyyydd[2]`
   echo 'converting obs for date: ' $datea
   
   rm -f flist
   /bin/ls -1 $DOWNLOAD_DIR/*.${yyyy}.${mday}.*_nc >! flist
   
   set nfiles = `cat flist | wc -l`
   if ( $chatty == 'yes' ) then
      echo $nfiles ' to process for file ' $datea 
   endif
   
   ./$CONV_PROG >>! convert_output_log

   rm -rf flist
   if ( -e obs_seq.gpsro ) then
       mv obs_seq.gpsro obs_seq.gpsro_${datea}
   
      if ( $chatty == 'yes' ) then
         echo "all observations for day in file obs_seq.gpsro_${datea} at " `date`
      endif
   else
      if ( $chatty == 'yes' ) then
         echo "no obs found for date $datea, or conversion failed.'
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

   echo 'removing original gpsro data files for date: ' $datea

   # just remove the whole subdir here.  trying to list individual
   # files can cause problems with long command line lengths.
   rm -fr $DOWNLOAD_DIR

   cd $DART_WORK_DIR
else
   echo 'not removing original gpsro data files'
endif

########################################################################

if ( $chatty == 'yes' ) then
   echo 'finished gpsro script run at ' `date`
   echo ' '
endif

exit 0


