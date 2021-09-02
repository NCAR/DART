#!/bin/tcsh

# a month at a time, create symbolic links in a parallel directory.
# this version of the script should be used when the files already
# exist in the CESM directory and you want to make links in an
# old-style name directory.  see link_cesm_to_dart.csh if you
# have files in an old-style directory and you want to make
# CESM links.

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

set THISDIR = ../${yr}${mo}_6H_CESM
set DARTDIR = ../${yr}${mo}_6H

if ( ! -d $THISDIR ) then
   echo $THISDIR must exist and contain CESM format obs_seq files before
   echo running this script.  exiting with an error.
   exit -1
endif

if ( ! -d $DARTDIR ) then
   mkdir $DARTDIR
endif

cd $DARTDIR

set d = 1
while ($d <= 31) 
   set dd = `printf %02d $d`
   if (-f $THISDIR/obs_seq.${yr}-${mo}-${dd}-00000) then
      ln -s $THISDIR/obs_seq.${yr}-${mo}-${dd}-00000 obs_seq${yr}${mo}${dd}00
   endif
   if (-f $THISDIR/obs_seq.${yr}-${mo}-${dd}-21600) then
      ln -s $THISDIR/obs_seq.${yr}-${mo}-${dd}-21600 obs_seq${yr}${mo}${dd}06
   endif
   if (-f $THISDIR/obs_seq.${yr}-${mo}-${dd}-43200) then
      ln -s $THISDIR/obs_seq.${yr}-${mo}-${dd}-43200 obs_seq${yr}${mo}${dd}12 
   endif
   if (-f $THISDIR/obs_seq.${yr}-${mo}-${dd}-64800) then
      ln -s $THISDIR/obs_seq.${yr}-${mo}-${dd}-64800 obs_seq${yr}${mo}${dd}18
   endif

   @ d++
end

exit 0
