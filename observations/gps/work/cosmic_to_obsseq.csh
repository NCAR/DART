#!/bin/csh 
########################################################################
#
#   cosmic_to_obsseq.csh - script that downloads COSMIC observations 
#               and converts them to a DART observation sequence file.
#
# requires 3 args:
#    $1 - analysis date (yyyymmdd format)
#    $2 - base observation directory
#    $3 - yes to download COSMIC data (calls cosmic_download.csh)
#    $4 - yes to delete raw COSMIC data when finished
#
# update the DART_DIR setting at the top to match your system.
#
#     created June 2008, Ryan Torn NCAR/MMM
#     updated Nov  2008, nancy collins ncar/cisl
#     updated Aug  2009, nancy collins ncar/cisl
#
# 
########################################################################

# should only have to set the DART_DIR and the rest should be in the
# right place relative to this directory.
setenv DART_DIR      /home/user/dart

setenv DART_WORK_DIR  ${DART_DIR}/observations/gps/work
setenv CONV_PROG      convert_cosmic_gps_cdf
setenv DATE_PROG      advance_time

set chatty=yes

set datea   = ${1}     # target date, YYYYMMDD
set datadir = ${2}     # where to process the files
set downld  = ${3}     # download?  'yes' or 'no'
set cleanup = ${4}     # delete COSMIC files at end? 'yes' or 'no'

set assim_freq          = 6  # hours, sets centers of windows.
set download_window     = 3  # window half-width (some users choose 2 hours)

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
if ( ! -e ${datadir}/input.nml ) then
  echo 'data processing directory does not contain an input.nml'
  echo 'copying from work dir to data proc dir'
  echo `pwd`/input.nml '->' ${datadir}/input.nml
  cp -f ./input.nml ${datadir}
else
  echo 'data processing directory already contains an input.nml'
  echo 'which will be used for this run.'
  diff -q ./input.nml ${datadir}/input.nml
endif

if ( ! -e ${datadir}/${DATE_PROG} ) then
  echo 'data processing directory does not contain the time converter'
  echo 'copying from work dir to data proc dir'
  echo `pwd`/${DATE_PROG} '->' ${datadir}/${DATE_PROG}
  cp -f ./${DATE_PROG} ${datadir}
else
  echo 'using time conversion program already found in data proc dir'
endif

if ( ! -e ${datadir}/${CONV_PROG} ) then
  echo 'data processing directory does not contain the data converter'
  echo 'copying from work dir to data proc dir'
  echo `pwd`/${CONV_PROG} '->' ${datadir}/${CONV_PROG}
  cp -f ./${CONV_PROG} ${datadir}
else
  echo 'using data conversion program already found in data proc dir'
endif

if ( $downld == 'yes' ) then

   ./cosmic_download.csh ${datea} ${datadir}

endif

echo 'changing dir to data proc directory'
cd ${datadir}
echo 'current dir now ' `pwd`

# save original
set dates = $datea

@ hh = $assim_freq + 100  ;  set hh = `echo $hh | cut -b2-3`
set datea = `echo $datea | cut -b1-8`${hh}
set datee = `echo ${datea}00 24 | ./${DATE_PROG}`

if ( $chatty == 'yes' ) then
   echo 'starting gps conversion at ' `date`
endif

while ( $datea < $datee )

  set yyyy = `echo $datea | cut -b1-4`
  set   mm = `echo $datea | cut -b5-6`
  set   dd = `echo $datea | cut -b7-8`
  set   hh = `echo $datea | cut -b9-10` 
  echo ${yyyy}-${mm}-${dd}_${hh}:00:00  > convert_cosmic_input

  rm -f flist
  @ nhours = 2 * $download_window  ;  set n = 1
  set datef = `echo $datea -$download_window | ./${DATE_PROG}`
  while ( $n <= $nhours ) 

    set    yyyy = `echo $datef | cut -b1-4`
    set      hh = `echo $datef | cut -b9-10`
    set jyyyydd = `./${DATE_PROG} $datef 0 -j`
    @ mday = $jyyyydd[2] + 1000  ;  set mday = `echo $mday | cut -b2-4`

    /bin/ls ${dates}/*.${yyyy}.${mday}.${hh}.*_nc >>! flist

    if ( $chatty == 'yes' ) then
      echo $datea ' hour ' $hh
    endif

    set datef = `echo $datef 1 | ./${DATE_PROG}`
    @ n += 1
  end
 
  set nfiles = `cat flist | wc -l`
  if ( $chatty == 'yes' ) then
      echo $nfiles ' to process for file ' $datea 
  endif

  ./convert_cosmic_gps_cdf < convert_cosmic_input >>! convert_output_log
  rm -rf cosmic_gps_input.nc convert_cosmic_input  flist

  set datef = $datea
  if ( `echo $datea | cut -b9-10` == '00' ) then
    set datef = `echo $datea -24 | ./${DATE_PROG}`
    set datef = `echo $datef | cut -b1-8`24
  endif
  if ( -e obs_seq.gpsro )  mv obs_seq.gpsro ../obs_seq.gpsro_${datef}

   if ( $chatty == 'yes' ) then
      echo "moved accumulated obs to ../obs_seq.gpsro_${datef} at " `date`
   endif

  set datea = `echo $datea $assim_freq | ./${DATE_PROG}`

end

if ( $cleanup == 'yes' ) then

   if ( $chatty == 'yes' ) then
      echo 'cleaning up files at ' `date`
   endif

  # do in chunks because some systems have problems with the command
  # line getting too long (large numbers of files here).
  rm -f ${dates}/atmPrf_C001*
  rm -f ${dates}/atmPrf_C002*
  rm -f ${dates}/atmPrf_C003*
  rm -f ${dates}/atmPrf_C004*
  rm -f ${dates}/atmPrf_C005*
  rm -f ${dates}/atmPrf_C006*
  rm -f ${dates}/atmPrf_CHAM*
  rm -f ${dates}/atmPrf_*

endif

if ( $chatty == 'yes' ) then
   echo 'finished gps conversion at ' `date`
endif


