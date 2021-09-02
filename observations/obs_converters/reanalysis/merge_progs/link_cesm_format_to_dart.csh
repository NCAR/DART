#!/bin/tcsh

# a month at a time, create symbolic links in a parallel directory.
# this version of the script should be used when the files already
# exist in an old-style dart directory and you want to make links 
# in a CESM named directory.  see link_dart_to_cesm.csh if you
# have files in a CESM dir and you want to make old-style links.

# usage: $0 yyyy mm

# or comment out the argument processing and hardcode a
# year and month.

# set a 4 digit year and a 2 digit month.
# make sure the month has a leading 0 if it is less than 10.
if ( $# != 2 ) then
   echo usage: $0 YYYY MM
   echo month must be 2 digits, with a leading 0 if less than 10.
   exit -1
endif

set yr = $1
set mo = $2

echo creating links for year, month: $1 $2

# or if you don't want to deal with arguments, comment out the
# code above and set them here.

#set yr = 2012
#set mo = 08

set THISDIR = ../${yr}${mo}_6H
set CESMDIR = ../${yr}${mo}_6H_CESM

if ( ! -d $THISDIR ) then
   echo $THISDIR must exist and contain original dart named obs_seq files 
   echo before running this script.  exiting with an error.
   exit -1
endif

if ( ! -d $CESMDIR ) then
   mkdir $CESMDIR
endif

cd $CESMDIR

set d = 1
while ($d <= 31) 
   set dd = `printf %02d $d`
   if (-f $THISDIR/obs_seq${yr}${mo}${dd}00) then
      ln -s $THISDIR/obs_seq${yr}${mo}${dd}00 obs_seq.${yr}-${mo}-${dd}-00000
   endif
   if (-f $THISDIR/obs_seq${yr}${mo}${dd}06) then
      ln -s $THISDIR/obs_seq${yr}${mo}${dd}06 obs_seq.${yr}-${mo}-${dd}-21600
   endif
   if (-f $THISDIR/obs_seq${yr}${mo}${dd}12) then
      ln -s $THISDIR/obs_seq${yr}${mo}${dd}12 obs_seq.${yr}-${mo}-${dd}-43200
   endif
   if (-f $THISDIR/obs_seq${yr}${mo}${dd}18) then
      ln -s $THISDIR/obs_seq${yr}${mo}${dd}18 obs_seq.${yr}-${mo}-${dd}-64800
   endif

   @ d++
end

exit 0
