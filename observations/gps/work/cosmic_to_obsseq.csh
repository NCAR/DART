#!/bin/csh 
########################################################################
#
#   cosmic_to_obsseq.csh - script that downloads COSMIC observations 
#               and converts them to a DART observation sequence file.
#
# requires 3 args:
#    $1 - analysis date (yyyymmdd format)
#    $2 - base observation directory
#    $3 - yes to download COSMIC data
#
# and update the 3 env settings at the top to match your system.
#
#     created June 2008, Ryan Torn NCAR/MMM
#     updated Nov  2008, nancy collins ncar/cisl
#
########################################################################

setenv DART_DIR      /home/user/dart
setenv DART_OBS_DIR ${DART_DIR}/observations/gps/work
setenv DATE_PROG    ./advance_time

set chatty=yes

set datea   = ${1}
set datadir = ${2}
set downld  = ${3}

set assim_freq      = 6
set download_window = 3   # was 2 hours
set gps_repository_path = 'http://cosmic-io.cosmic.ucar.edu/cdaac/login'

if ( ! -e ${datadir}/cosmic )  exit
cd ${datadir}/cosmic

if ( $downld == 'yes' ) then

   if ( ! $?cosmic_user || ! $?cosmic_pw ) then
      echo "You must setenv cosmic_user to your username for the cosmic web site"
      echo "and setenv cosmic_pw to your password. (or export cosmic_user=val for"
      echo "ksh/bash users) then rerun this script. "
      exit -1
   endif

  if ( $chatty == 'yes' ) then
     echo 'starting raw file download at' `date`
  endif

  set yyyy = `echo $datea | cut -b1-4`

  set jyyyydd = `${DATE_PROG} $datea 0 -j` 
  @ mday = $jyyyydd[2] + 1000  ;  set mday = `echo $mday | cut -b2-4`
  wget -q -nd -np -r -l 10 --http-user=${cosmic_user} --http-passwd=${cosmic_pw} ${gps_repository_path}/cosmic/level2/atmPrf/${yyyy}.${mday}/
  rm -f *.html *.txt
  wget -q -nd -np -r -l 10 --http-user=${cosmic_user} --http-passwd=${cosmic_pw} ${gps_repository_path}/champ/level2/atmPrf/${yyyy}.${mday}/
  rm -f *.html *.txt

  set jyyyydd = `${DATE_PROG} $datea 24 -j`
  @ mday = $jyyyydd[2] + 1000  ;  set mday = `echo $mday | cut -b2-4` 
  wget -q -nd -np -r -l 10 --http-user=${cosmic_user} --http-passwd=${cosmic_pw} ${gps_repository_path}/cosmic/level2/atmPrf/${yyyy}.${mday}/
  rm -f *.html *.txt
  wget -q -nd -np -r -l 10 --http-user=${cosmic_user} --http-passwd=${cosmic_pw} ${gps_repository_path}/champ/level2/atmPrf/${yyyy}.${mday}/
  rm -f *.html *.txt

  if ( $chatty == 'yes' ) then
     # the ls arg list line gets too long in some cases
     echo `/bin/ls | grep _nc | wc -l` 'raw files'
     echo 'all raw files download at ' `date`
  endif

endif

@ hh = $assim_freq + 100  ;  set hh = `echo $hh | cut -b2-3`
set datea = `echo $datea | cut -b1-8`${hh}
set datee = `${DATE_PROG} ${datea}00 24`

if ( $chatty == 'yes' ) then
   echo 'starting gps conversion at ' `date`
endif

while ( $datea < $datee )

  set yyyy = `echo $datea | cut -b1-4`
  set   mm = `echo $datea | cut -b5-6`
  set   dd = `echo $datea | cut -b7-8`
  set   hh = `echo $datea | cut -b9-10` 
  echo ${yyyy}-${mm}-${dd}_${hh}:00:00  > convert_cosmic_input

  @ nhours = 2 * $download_window  ;  set n = 1
  set datef = `${DATE_PROG} $datea -$download_window`
  while ( $n <= $nhours ) 

    set    yyyy = `echo $datef | cut -b1-4`
    set      hh = `echo $datef | cut -b9-10`
    set jyyyydd = `${DATE_PROG} $datef 0 -j`
    @ mday = $jyyyydd[2] + 1000  ;  set mday = `echo $mday | cut -b2-4`

    set f = 1  ;  set nfiles=`/bin/ls *.${yyyy}.${mday}.${hh}.*_nc | wc -l`
    if ( $chatty == 'yes' ) then
        echo $nfiles ' to process for file ' $datea ' hour ' $hh
    endif
    while ( $f <= $nfiles )

      ln -sf `/bin/ls *.${yyyy}.${mday}.${hh}.*_nc | head -n $f | tail -1` cosmic_gps_input.nc
      ${DART_OBS_DIR}/convert_cosmic_gps_cdf < convert_cosmic_input >>! convert_output_log
      rm -rf cosmic_gps_input.nc
      @ f += 1

    end
    set datef = `${DATE_PROG} $datef 1`
    @ n += 1
 
  end
  rm -rf convert_cosmic_input

  set datef = $datea
  if ( `echo $datea | cut -b9-10` == '00' ) then
    set datef = `${DATE_PROG} $datea -24`
    set datef = `echo $datef | cut -b1-8`24
  endif
  if ( -e obs_seq.gpsro )  mv obs_seq.gpsro ../obs_seq.gpsro_${datef}

   if ( $chatty == 'yes' ) then
      echo "moved accumulated obs to ../obs_seq.gpsro_${datef} at " `date`
   endif

  set datea = `${DATE_PROG} $datea $assim_freq`

end

if ( $downld == 'yes' ) then

   if ( $chatty == 'yes' ) then
      echo 'cleaning up files at ' `date`
   endif

  rm -f atmPrf_C001*
  rm -f atmPrf_C002*
  rm -f atmPrf_C003*
  rm -f atmPrf_C004*
  rm -f atmPrf_C005*
  rm -f atmPrf_C006*
  rm -f atmPrf_CHAM*
  rm -f atmPrf_*

endif

if ( $chatty == 'yes' ) then
   echo 'finished gps conversion at ' `date`
endif


