#!/usr/bin/env bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

fdate=${1}
otype=${2}
outfile=${3}

if [ $fdate == '-h' ]; then

  echo '  '
  echo '   get_madis_archive - script that ftps observation data from FSL   '
  echo '               MADIS archive database.                              '
  echo '  '
  echo '    arg1 - yyyymmddhh date of observation '
  echo '    arg2 - type of observation to get:'
  echo '             HDW - GOES cloud track winds'
  echo '             acars - ACARS observations'
  echo '             maritime - buoys and ships'
  echo '             metar - metar data'
  echo '             profiler - NOAA wind profiler'
  echo '             radiometer -' 
  echo '             raob - radiosonde observations'
  echo '             sao - surface observations foreign'
  echo '             satrad - microwave satellite radiances'
  echo '             snow  - snow reports'
  echo '    arg3 - name for observation file'
  echo '  '

  exit
fi
		     
yyyy=`echo $fdate | cut -b1-4`
mm=`echo $fdate | cut -b5-6`
dd=`echo $fdate | cut -b7-8`
hh=`echo $fdate | cut -b9-10`

cat > ftp_madis << END_FTP 
open rftp.madis-data.noaa.gov
################# below change to match your madis acct info ##########################
user your_madis_username your_madis_password
#######################################################################################
binary
passive
get /research/archive/${yyyy}/${mm}/${dd}/point/${otype}/netcdf/${yyyy}${mm}${dd}_${hh}00.gz ${outfile}.gz
get /research/archive/${yyyy}/${mm}/${dd}/LDAD/${otype}/netCDF/${yyyy}${mm}${dd}_${hh}00.gz  ${outfile}.gz
get /research2/point/${otype}/netcdf/${yyyy}${mm}${dd}_${hh}00.gz ${outfile}.gz
quit
END_FTP

ftp -n < ftp_madis > /dev/null
\rm ftp_madis

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

