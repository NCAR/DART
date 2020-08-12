#!/bin/csh 
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
#--------------------------------------------------------------
# DESCRIPTION:
#
#  This script is used to generate 1 day's worth of observations
#  in an intermediate text file format.  it can output either a
#  single daily (3:01Z to 3:00Z of next day) file, or 4 6-hr files.
#  the text file contains decoded NCEP PREPBUFR text/ascii data.
#
#  It requires 3 arguments: year month day 
#
#--------------------------------------------------------------
# USER SET PARAMETERS

# if daily is 'yes', 4 6-hour files will be processed and a single, 1-day
# output file will be created.  if daily is set to anything else, then
# each 6-hour input file will be converted to a single 6-hour output file.
# this script still processes a day's worth of files at a time even if
# daily is 'no' - it just makes 4 individual files per day.

set    daily = no

# if daily is 'no' and zeroZ is 'yes', input files at 0Z will be translated 
# into output files also marked 0Z.  otherwise, they will be named with the 
# previous day number and 24Z (chose here based on what script will be 
# processing these files next.  the 'create_real_obs' script assumes 
# filenames with the pattern 6Z,12Z,18Z,24Z, so 'no' is correct for it.)
# this variable is ignored completely if daily is 'yes'.

set zeroZ = no

# if convert is 'yes', then the big-endian BUFR files will be converted
# to little-endian files before processing. this is needed if you are running
# on a machine that uses Intel chips (e.g. linux clusters, altix, pcs, etc).
# it is not needed for ibm power systems.  any value other than 'yes' will
# skip the convert step.

set  convert = yes

# if block is 'yes', then the cword program will be run to convert an
# unblocked file into a blocked one.  this is not required for recent
# prepbufr files, but older ones may require it.

set block = no

# starting year, month, day, and ending day.  this script does not allow
# you to do more than a single month at a time, but does handle the last
# day of the month, leap day in feb, and the last day of the year correctly.
# this version of the conversion tool takes up to 3 hours of observations
# from the day *following* the end day, so you must have at least the 6Z
# file from one day beyond the last day for this to finish ok (to get obs
# at exactly 3Z from the next file).

if ( $#argv == 3 ) then
   set year   = $argv[1]
   set month  = $argv[2]
   set day    = $argv[3]
else
   echo usage: $0 year month day
   exit -1
endif

# all output files are put into a single directory with names that follow
# the filename pattern 'prepqmYYMMDDHH' for 6-hr files, or 'prepqmYYMMDD' 
# for daily files. ("daily" here means from 3Z+1sec to 3Z on day+1)
#
# the format of the input prepqm filenames have changed over the years; if 
# this script no longer matches what is untarred from the files, you will 
# have to edit this script below.
#
# there are several shell variables in the loop you can use to construct
# alternate names:
# 'year' is 4 digits; 'yr' is 2.
# 'mn', 'dy', and 'hr' are 0 padded so they are always 2 digits.
# 'oyear', 'omn', 'ody', 'ohr' are the original date, if the day, month, 
# and/or year have rolled over.

#set BUFR_dir = ../data
set BUFR_dir  = /glade/p/image/Observations/bufr

set BUFR_idir = ${BUFR_dir}/prepqm
set BUFR_odir = ${BUFR_dir}/prepout


# directory where DART prepbufr programs are located, relative to here.
set DART_exec_dir = .

# END USER SET PARAMETERS
#--------------------------------------------------------------

# make sure these variables are 2 chars wide
# which means adding leading 0s if less than 10.

set yr = `echo $year | cut -c3-4`
set mn = `printf %02d $month` 
set dy = `printf %02d $day` 
set hr = "06"

# date-time group, and 'short' date-time group where
# year is only last 2 digits (matches filename conventions)

set dtg  = ${year}${mn}${dy}${hr}
set sdtg = ${yr}${mn}${dy}${hr}

# these files come every 6 hours.

set inc = +6h
set files_per_day = 4

# daily files include obs at exactly 3Z from day+1
# so they have to read 5 files instead of 4.
if ($daily == 'yes') then
   set files_per_day = 5 
endif

set hr_file = 1
while ( $hr_file <= $files_per_day )

   # clear any old intermediate (text) files
   rm -f temp_obs prepqm.in prepqm.out 

   # the prepqm input files.  match general naming pattern with the
   # data time encoded in the filename.  if the pattern of the filename
   # is different (4 digit year vs 2, extra fixed text in the name, etc)
   # fix the BUFR_in line below to match what you have.  if the file is
   # gzipped, you can leave it and this program will unzip it before
   # processing it.
   set BUFR_in = ${BUFR_idir}/prepqm${sdtg}.no_ship_id

   if ( -e ${BUFR_in} ) then
      echo "copying ${BUFR_in} into prepqm.in"
      rm -f prepqm.in
      cp -f ${BUFR_in} prepqm.in
   else if ( -e ${BUFR_in}.gz ) then
      echo "unzipping ${BUFR_in}.gz into prepqm.in"
      rm -f prepqm.in
      gunzip -c -f ${BUFR_in}.gz >! prepqm.in
   else
      echo "MISSING INPUT FILE: cannot find either"
      echo ${BUFR_in}
      echo   or 
      echo ${BUFR_in}.gz
      echo "Script will abort now."
      exit -1
   endif

   # blocking
   if ($block == 'yes') then
      echo "blocking prepqm.in"
      mv -f prepqm.in prepqm.unb
      echo 'block' >! in
      echo 'prepqm.unb' >> in
      echo 'prepqm.blk' >> in
      ${DART_exec_dir}/cword.x < in
      mv -f prepqm.blk prepqm.in
      rm -f prepqm.unb in
   endif

   # byte swapping
   if ($convert == 'yes') then
      echo "byteswapping bigendian to littleendian prepqm.in"
      mv -f prepqm.in prepqm.big
      ${DART_exec_dir}/grabbufr.x prepqm.big prepqm.lit
      mv -f prepqm.lit prepqm.in
      rm -f prepqm.big
   endif

   if ($hr_file == 5) then
      # get any obs exactly at 3Z from the 6Z file of the next day using a
      # modified prepbufr program, since 6 hour assimilation windows 
      # centered on 6Z, 12Z, etc would be 03:01Z-09:00Z, 09:01Z-15:00Z, etc.
      # obs exactly at 03Z need to be part of the same file which spans 
      # 21:01Z-03:00Z.   both of these prepbufr programs also add 24h to 
      # any time after midnight, so the hours in the output files run
      # from 3.016 hours to 27.000 hours.  (the ascii intermediate files
      # do not contain day numbers, but probably should and then the hours
      # can run from a normal 0Z to 23:59Z.)
      ${DART_exec_dir}/prepbufr_03Z.x
   else
      ${DART_exec_dir}/prepbufr.x
   endif

   if ($daily == 'yes') then
      cat prepqm.out >>! temp_obs
      rm -f prepqm.out
   else
      if ($hr == "00") then   # last file of the day; possibly handle naming special
         if ($zeroZ == 'yes') then
            # if 0Z, output named with current day and 0Z
            set BUFR_out = ${BUFR_odir}/temp_obs.${dtg}
         else
            # if not 0Z, output named with previous day and 24Z
            set BUFR_out = ${BUFR_odir}/temp_obs.${oyear}${omn}${ody}24
         endif
      else
         set BUFR_out = ${BUFR_odir}/temp_obs.${dtg}
      endif
      mkdir -p ${BUFR_odir}/
      echo "moving output to $BUFR_out"
      mv -fv prepqm.out $BUFR_out
   endif

   set oyear = $year
   set oyr   = $yr
   set omn   = $mn
   set ody   = $dy
   set ohr   = $hr

   @ hr_file++

   set dtg  = `echo $dtg $inc | ./advance_time`

   set year = `echo $dtg | cut -c 1-4`
   set yr   = `echo $dtg | cut -c 3-4`
   set mn   = `echo $dtg | cut -c 5-6`
   set dy   = `echo $dtg | cut -c 7-8`
   set hr   = `echo $dtg | cut -c 9-10`

   set sdtg = ${yr}${mn}${dy}${hr}

end

if ($daily == 'yes') then
   # use the original dates without rollover
   set BUFR_out = ${BUFR_odir}/temp_obs.${oyear}${omn}${ody}
   mkdir -p ${BUFR_odir}/
   echo "moving output to ${BUFR_out}"
   mv -fv temp_obs ${BUFR_out}
endif

rm -f prepqm.in

exit 0


