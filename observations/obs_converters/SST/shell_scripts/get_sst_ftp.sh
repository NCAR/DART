#!/bin/bash
#
# This code was donated to DART by Romain Escudier, who was at Rutgers at the time.
# It is not protected by the DART copyright agreement. Thanks Romain!
#
# DART $Id$

if [ $# -gt 0 ]; then
   YEAR=$1
else
   echo "must supply a year as an argument"
   exit 1
fi

. ./parameters_SST

mkdir -p ${DIR_IN}/${YEAR}
cd ${DIR_IN}/${YEAR}

for i_day in $( seq 1 366 ) ; do
   printf -v DDAY "%03d" $i_day

   wget ${FTP_ADDRESS}/${YEAR}/${DDAY}/*.nc*
   cmd_status=$?

   if (( cmd_status==0 )); then
      file_tmp=$(ls *.nc.md5 | awk -F ".md5" '{print $1}')
      md5file=$(md5sum ${file_tmp} | awk '{print $1}')
      md5veri=$(cat ${file_tmp}.md5 | awk '{print $1}')
      if [ "$md5file" != "$md5veri" ]; then
         echo "Error! Files md5 not matching!"
         break
      fi
      rm ${file_tmp}.md5
      bunzip2 ${file_tmp}
   fi

done

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

