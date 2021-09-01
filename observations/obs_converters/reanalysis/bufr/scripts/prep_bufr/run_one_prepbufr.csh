#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id: prepbufr.csh 10982 2017-02-01 23:43:10Z thoar@ucar.edu $
#


#
# run a single instance of prepbufr.  
# assumes you are starting this from one dir down
# from the work dir under prep_bufr (so ../../exe 
# contains the prepbufr # converter programs, etc)
#
# usage: $0 unique_subdir_name year month day
#

if ( $# != 4 ) then
  echo usage: $0 unique_subdir_name year month day
  exit -1
endif

# make a unique directory to work in

mkdir -p $1
cd $1

# copy the executables and namelist files we need

foreach i ( ../*.x )
  ln -sf $i .
end

ln -sf ../advance_time .

cp -f ../my_prepbufr.csh .
cp -f ../input.nml .

# and now call the script that actually does the work

./my_prepbufr.csh $2 $3 $4 

# clean up things which don't have any log info in them

rm *.x
rm advance_time
rm input.nml 
rm my_prepbufr.csh
if ( -z dart_log.out ) rm dart_log.out
cd ..

exit 0

# <next few lines under version control, do not edit>
# $URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/rma_trunk/observations/obs_converters/NCEP/prep_bufr/work/prepbufr.csh $
# $Revision: 10982 $
# $Date: 2017-02-01 16:43:10 -0700 (Wed, 01 Feb 2017) $

